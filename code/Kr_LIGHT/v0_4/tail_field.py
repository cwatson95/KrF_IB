"""v0_4: the two-temperature closure's tail ratio for field-driven pure
krypton, from this model's own kinetic EEDF.

v0_1 to v0_3 read r = T_t / Te from ML's table of beam-driven KrF states
(data/kr_tail_ratio.txt, its F2 = 0 column). Along the reference run, a
field-driven discharge, the kinetic EEDF's own r is 0.12 to 0.16 over the
pulse body where that column gives 0.18 to 0.20, so the closure's tail was
too hot and its 147-nm share 12 % over the kinetic one (the 2026-10-15
groundwork, on v0_3). This module builds the table v0_4's closure reads,
r(n_e/N, Te), from field-driven states solved as the run solves them.

Each state: ground-state Kr at TG_K and the density of 1 atm at TG_K, as
many cations as electrons (their Coulomb collisions), a DC field at the
mean energy 1.5 Te (the mean match), electron-electron collisions at n_e/N
and the growth term where it holds (Solver.eedf with the run's processes,
cells and tolerance), no excited population. r is the T_t / Te at which the
two-temperature EEDF at Te removes the same energy through the tail rows
(the 51 levels and the ionization, each at its threshold) as the kinetic
EEDF does, the definition of ML's beam table (ML/ml/tail.py). Where the
kinetic tail is colder than the two-temperature form can follow (low n_e/N
and low Te: no r down to R_RANGE[0] removes as little as the kinetic EEDF,
or its loss is below 1e-30 of the Maxwellian's), r is R_RANGE[0], the
closest the form comes; xsec holds T_t at its grid's floor, 0.01 eV, in any
case. The closure reads the table bilinear in (log10 n_e/N, ln Te) on
ln r.

    python3 -m v0_4.tail_field          # from LIGHT_Kr: writes the table

prints the file's md5, which xsec.TAIL_FIELD_FILE pins (regenerate, never
edit the file)."""

import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import boltzmann, xsec

# the table's axes: log10(n_e/N) and Te [eV]
X_GRID = tuple(round(x, 3) for x in np.arange(-9.0, -2.49, 0.5))
TE_GRID = tuple(round(2.0 ** e, 6) for e in np.arange(-1.0, 3.01, 0.25))
TG_K = 300.0
P_ATM = 1.0
R_RANGE = (1e-3, 3.0)
MATCH_FLOOR = 1e-30      # kinetic loss / Maxwellian loss below it: held


def _models():
    from .model import Model, Params
    kin = Model(Params(eedf="kinetic"))
    two = Model(Params(eedf="two_t"))
    return kin, two


def _tail_rows(kin, two):
    """(threshold [eV], kinetic Process, two-temperature table) of each
    tail row, in the run's order, with its factor (1 at the defaults)."""
    proc = {j: p for j, p in kin._kin_xsec}
    out = []
    for c in two.compiled:
        if c.kind != "xsec_tail":
            continue
        if c.j not in proc:              # a row at factor 0 (HM_ at hm_kr 0)
            if c.factor != 0.0:
                raise RuntimeError(f"{c.row.rid}: not in the kinetic set")
            continue
        p = proc[c.j]
        out.append((p.threshold, p, c.rate, c.factor))
    return out


def kinetic_state(ne_over_n, n_col_m3, tg=TG_K):
    ne = ne_over_n * n_col_m3
    return boltzmann.State(x={"Kr": 1.0}, tg=tg, n_col=n_col_m3, ne=ne,
                           n_ion=ne)


def solve_r(solver, rows, ne_over_n, te, n_col_m3, tg=TG_K, warm=None):
    """r at one state, the kinetic solution and the matched quantities:
    (r or nan, sol, kinetic loss, Maxwellian loss) [eV cm^3/s]."""
    s = kinetic_state(ne_over_n, n_col_m3, tg)
    sol = solver.eedf(s, target_mean=1.5 * te, warm=warm)
    g, F = sol.grid, sol.F
    te_k = sol.te
    target = sum(thr * g.rate(p, F) * 1e6 for thr, p, _, _ in rows)

    def loss(r):
        return sum(thr * f * tab(te_k, r * te_k) for thr, _, tab, f in rows)
    maxw = loss(1.0)
    if not target > MATCH_FLOOR * maxw:
        return float("nan"), sol, target, maxw
    fn = lambda lr: math.log(loss(math.exp(lr)) / target)
    lo, hi = (math.log(v) for v in R_RANGE)
    if fn(lo) > 0.0 or fn(hi) < 0.0:
        return float("nan"), sol, target, maxw
    return math.exp(brentq(fn, lo, hi, xtol=1e-10)), sol, target, maxw


def build(x_grid=X_GRID, te_grid=TE_GRID, tg=TG_K, p_atm=P_ATM,
          verbose=False):
    """The table r[x, te], with what each point matched."""
    kin, two = _models()
    rows = _tail_rows(kin, two)
    n_col = p_atm * 101325.0 / (1.380649e-23 * tg)
    r = np.full((len(x_grid), len(te_grid)), np.nan)
    info = []
    for i, x in enumerate(x_grid):
        for j, te in enumerate(te_grid):
            # each point from a cold start on a fresh solver at the run's
            # first grid, so a point solved alone is the table's (a warm
            # sweep along Te landed on other grids, 1e-3 apart in r)
            solver = boltzmann.Solver(
                kin.kin.procs, n_cells=kin.p.kin_cells,
                emax=kin.p.kin_emax_eV, tol=kin.p.kin_tol, elastic=("Kr",))
            t0 = time.time()
            val, sol, target, maxw = solve_r(solver, rows, 10.0 ** x, te,
                                             n_col, tg)
            r[i, j] = val
            info.append(dict(x=x, te=te, r=val, te_kin=sol.te,
                             en_td=sol.en_td, growth_off=sol.growth_off,
                             loss_kin=target, loss_maxw=maxw,
                             emax=sol.grid.emax, s=time.time() - t0))
            if verbose:
                print(f"x {x:6.2f} Te {te:6.3f} r {val:8.5f} E/N "
                      f"{sol.en_td:9.3f} Td kin/maxw "
                      f"{target / maxw if maxw > 0 else float('nan'):9.3e} "
                      f"({time.time() - t0:.2f} s)", flush=True)
    held = np.isnan(r)
    r[held] = R_RANGE[0]          # colder than the form follows: its floor
    return r, held, info


def write(path, r, held, x_grid=X_GRID, te_grid=TE_GRID, tg=TG_K,
          p_atm=P_ATM):
    lines = [
        "# LIGHT_Kr v0_4's two-temperature tail: r = T_tail / Te for the Kr "
        "excitation and ionization rows",
        "# (xsec.two_t_table). Rows: log10(n_e/N). Columns: Te [eV]. "
        "Generated by v0_4/tail_field.py",
        f"# (python3 -m v0_4.tail_field) on {time.strftime('%Y-%m-%d')}. "
        "Do not edit: regenerate.",
        f"# Each state: field-driven pure Kr at {tg:g} K and {p_atm:g} atm, "
        "as many cations as electrons, no excited",
        "# population; its kinetic EEDF (boltzmann.py, the run's processes, "
        "e-e collisions and growth) at the mean",
        "# energy 1.5 Te sets r: the two-temperature EEDF at Te removes the "
        "same energy through the 51 levels",
        "# and the ionization. Read bilinear in (log10 n_e/N, ln Te) on ln r "
        f"and held at its edges. {int(held.sum())} of {held.size} points",
        f"# (low n_e/N and Te) are at r's floor, {R_RANGE[0]:g}: the "
        "kinetic tail is colder than the form follows.",
        "TE " + " ".join(f"{t:g}" for t in te_grid),
    ]
    for x, row in zip(x_grid, r):
        lines.append(f"{x:.3f} " + " ".join(f"{v:.5f}" for v in row))
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv=None):
    path = xsec.DATA_DIR / xsec.TAIL_FIELD_FILE[0]
    t0 = time.time()
    r, held, _ = build(verbose=True)
    write(path, r, held)
    print(f"wrote {path} ({time.time() - t0:.0f} s), md5 "
          f"{xsec.md5_of(path)}: pin it in xsec.TAIL_FIELD_FILE")


if __name__ == "__main__":
    sys.exit(main())
