"""v0_4: the two-temperature closure's tail ratio for field-driven pure
krypton, from this model's own kinetic EEDF; v0_5 (2026-10-16): with an
excited-population axis.

v0_1 to v0_3 read r = T_t / Te from ML's table of beam-driven KrF states
(data/kr_tail_ratio.txt, its F2 = 0 column). v0_4 built r(n_e/N, Te) from
field-driven states with no excited atom, and its r came out 16 % below
the kinetic run's own along the reference run: the run's excited atoms
shape its tail (the 2026-10-15 entry, Known-bad: the fix is an
excited-fraction axis, never a factor on r). v0_5's table is
r(n_e/N, Te, x_s), x_s = [Kr(5s)] / [Kr], the excited populations at each
point set by the model's own rates:

* the four 5s levels at x_s by statistical weight, g/12;
* the 5p rung in its steady state against them at the point's Te and n_e:
  the electrons' pumping out of each 5s level (the K76 pairs, Maxwellian at
  Te) over its losses (relaxation by Kr, radiation, the electrons' return
  and ionization), its levels at the gas temperature (excited.py);
* Kr2* in its steady state: the 5s levels' formation at 300 K over its
  radiation and its ionization by electrons at Te.

Each state: ground-state Kr at TG_K and the density of P_ATM at TG_K, as
many cations as electrons (their Coulomb collisions), those populations, a
DC field at the mean energy 1.5 Te (the mean match), electron-electron
collisions and the growth term where it holds (Solver.eedf with the run's
processes, its population processes included, cells and tolerance). r is
the T_t / Te at which the two-temperature EEDF at Te removes the same
energy through the tail rows (the 51 levels and the ionization, each at its
threshold) as the kinetic EEDF does, the definition of ML's beam table
(ML/ml/tail.py). Where the kinetic tail is colder than the two-temperature
form can follow (no r down to R_RANGE[0] removes as little, or its loss is
below 1e-30 of the Maxwellian's), r is R_RANGE[0]. Where it is hotter than
the form can follow, r is R_RANGE[1] (v0_5): the 5s atoms' superelastic
electrons feed the tail at 10 eV directly, so at low Te with excited atoms
the kinetic tail removes up to some 400 times the Maxwellian's energy,
which a tail tied to the Maxwellian bulk at the break cannot reach (the
first build held those 175 points at the floor, the wrong end). The
closure reads the
table trilinear in (log10 n_e/N, ln Te, log10(x_s + XS0)) on ln r, held at
its edges.

    python3 -m v0_5.tail_field [workers]   # from LIGHT_Kr: writes the table

prints the file's md5, which xsec.TAIL_FIELD_FILE pins (regenerate, never
edit the file). Every point is solved from a cold start on a fresh
solver, so the table does not depend on the order or the workers."""

import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import boltzmann, excited, xsec

# the table's axes: log10(n_e/N), Te [eV] and x_s = [Kr(5s)]/[Kr]
X_GRID = tuple(round(x, 3) for x in np.arange(-9.0, -2.49, 0.5))
TE_GRID = tuple(round(2.0 ** e, 6) for e in np.arange(-1.0, 3.01, 0.25))
XS_GRID = (0.0, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
XS0 = 1e-7                  # the axis coordinate log10(x_s + XS0)
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
    """(threshold [eV], kinetic Process, two-temperature table, factor) of
    each tail row, in the run's order."""
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


def composition(two, ne_over_n, te, x_s, n_col_cm3, tg=TG_K) -> dict:
    """The excited populations of one table state as fractions of the
    colliding gas: the 5s levels at x_s by weight, the rung and Kr2* in
    their steady states at Te and n_e (the closure's own Maxwellian rates,
    two._k), the rung's levels at tg."""
    from .network import G_LEVEL_5S, G_5S, LEVELS_5S, kr5s
    pops = {}
    if x_s <= 0.0:
        return pops
    ne = ne_over_n * n_col_cm3
    n = np.zeros(two.ns)
    n[two.idx["Kr"]] = n_col_cm3
    for lv in LEVELS_5S:
        pops[kr5s(lv)] = x_s * G_LEVEL_5S[lv] / G_5S
        n[two.idx[kr5s(lv)]] = pops[kr5s(lv)] * n_col_cm3

    def k(c):
        v = two._k(c, te, tg, n, ne)
        return v * (excited.rung_fraction(c.rung, tg) if c.rung else 1.0)
    i_pp, i_2 = two.idx["Kr**"], two.idx["Kr2*"]
    make_pp = loss_pp = make_2 = loss_2 = 0.0
    for c in two.compiled:
        r = c.row
        heavy = [two.species[i] for i in c.heavy]
        if "Kr**" in r.products and "Kr**" not in r.reactants and \
                r.rid.startswith("K76_"):
            make_pp += k(c) * ne * n[c.heavy[0]]
        elif heavy and heavy[0] == "Kr**" and "Kr**" not in r.products:
            f = k(c) * (ne if c.ne_pow else 1.0)
            for i in c.heavy[1:]:
                f *= n[i]
            loss_pp += f
        if "Kr2*" in r.products and r.rid.startswith("AK8_"):
            f = k(c)
            for i in c.heavy:
                f *= n[i]
            make_2 += f
        elif heavy and heavy[0] == "Kr2*" and "Kr2*" not in r.products:
            f = k(c) * (ne if c.ne_pow else 1.0)
            for i in c.heavy[1:]:
                f *= n[i]
            loss_2 += f
    n_pp = make_pp / loss_pp if loss_pp > 0 else 0.0
    n_2 = make_2 / loss_2 if loss_2 > 0 else 0.0
    pops["Kr**"] = n_pp / n_col_cm3
    for x in excited.RUNG:
        pops[f"Kr**:{x}"] = pops["Kr**"] * excited.rung_fraction(x, tg)
    pops["Kr2*"] = n_2 / n_col_cm3
    pops["Kr2*(0u+)"] = 0.0
    return pops


def kinetic_state(ne_over_n, n_col_m3, tg=TG_K, pops=None):
    ne = ne_over_n * n_col_m3
    return boltzmann.State(x={"Kr": 1.0}, tg=tg, n_col=n_col_m3, ne=ne,
                           n_ion=ne, pops=dict(pops or {}))


def solve_r(solver, rows, ne_over_n, te, n_col_m3, tg=TG_K, warm=None,
            pops=None):
    """r at one state, the kinetic solution and the matched quantities:
    (r, sol, kinetic loss, Maxwellian loss) [eV cm^3/s]; r is nan where
    the kinetic tail is colder than the form follows and inf where it is
    hotter (v0_5)."""
    s = kinetic_state(ne_over_n, n_col_m3, tg, pops)
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
    if fn(lo) > 0.0:
        return float("nan"), sol, target, maxw
    if fn(hi) < 0.0:
        return float("inf"), sol, target, maxw
    return math.exp(brentq(fn, lo, hi, xtol=1e-10)), sol, target, maxw


_CACHE = {}


def _point(args):
    """One table point from a cold start on a fresh solver (a worker's)."""
    i, j, k, x, te, x_s, tg, p_atm = args
    if "models" not in _CACHE:
        kin, two = _models()
        _CACHE["models"] = (kin, two, _tail_rows(kin, two))
    kin, two, rows = _CACHE["models"]
    n_col = p_atm * 101325.0 / (1.380649e-23 * tg)
    pops = composition(two, 10.0 ** x, te, x_s, n_col * 1e-6, tg)
    solver = boltzmann.Solver(
        kin.kin.procs, n_cells=kin.p.kin_cells, emax=kin.p.kin_emax_eV,
        tol=kin.p.kin_tol, elastic=("Kr",), pops=kin.kin.pops)
    t0 = time.time()
    val, sol, target, maxw = solve_r(solver, rows, 10.0 ** x, te, n_col,
                                     tg, pops=pops)
    return (i, j, k, val, dict(x=x, te=te, x_s=x_s, r=val, te_kin=sol.te,
                               en_td=sol.en_td, growth_off=sol.growth_off,
                               loss_kin=target, loss_maxw=maxw,
                               emax=sol.grid.emax, s=time.time() - t0,
                               x_pp=pops.get("Kr**", 0.0),
                               x_2=pops.get("Kr2*", 0.0)))


def build(x_grid=X_GRID, te_grid=TE_GRID, xs_grid=XS_GRID, tg=TG_K,
          p_atm=P_ATM, workers=1, verbose=False):
    """The table r[x, te, x_s], with what each point matched."""
    jobs = [(i, j, k, x, te, xs, tg, p_atm)
            for k, xs in enumerate(xs_grid)
            for i, x in enumerate(x_grid) for j, te in enumerate(te_grid)]
    r = np.full((len(x_grid), len(te_grid), len(xs_grid)), np.nan)
    info = []
    if workers > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers,
                                 mp_context=mp.get_context("spawn")) as ex:
            results = ex.map(_point, jobs, chunksize=4)
            for n_done, (i, j, k, val, d) in enumerate(results):
                r[i, j, k] = val
                info.append(d)
                if verbose and n_done % 100 == 0:
                    print(f"{n_done + 1}/{len(jobs)}", flush=True)
    else:
        for job in jobs:
            i, j, k, val, d = _point(job)
            r[i, j, k] = val
            info.append(d)
            if verbose:
                print(f"x {d['x']:6.2f} Te {d['te']:6.3f} x_s {d['x_s']:.0e} "
                      f"r {val:8.5f} ({d['s']:.2f} s)", flush=True)
    cold, hot = np.isnan(r), np.isinf(r)
    r[cold] = R_RANGE[0]          # colder than the form follows: its floor
    r[hot] = R_RANGE[1]           # v0_5: hotter than it follows: its top
    return r, cold | hot, info


def write(path, r, held, x_grid=X_GRID, te_grid=TE_GRID, xs_grid=XS_GRID,
          tg=TG_K, p_atm=P_ATM):
    lines = [
        "# LIGHT_Kr v0_5's two-temperature tail: r = T_tail / Te for the Kr "
        "excitation and ionization rows",
        "# (xsec.two_t_table). Blocks: x_s = [Kr(5s)]/[Kr] (line 'XS'). "
        "Rows: log10(n_e/N). Columns: Te [eV].",
        "# Generated by v0_5/tail_field.py (python3 -m v0_5.tail_field) on "
        f"{time.strftime('%Y-%m-%d')}. Do not edit: regenerate.",
        f"# Each state: field-driven pure Kr at {tg:g} K and {p_atm:g} atm, "
        "as many cations as electrons, the 5s",
        "# levels at x_s by weight, the 5p rung and Kr2* in their steady "
        "states at Te and n_e (the model's rates);",
        "# its kinetic EEDF (boltzmann.py: the run's processes, population "
        "processes, e-e and growth) at the mean",
        "# energy 1.5 Te sets r: the two-temperature EEDF at Te removes the "
        "same energy through the 51 levels",
        "# and the ionization. Read trilinear in (log10 n_e/N, ln Te, "
        f"log10(x_s + {XS0:g})) on ln r, held at its edges.",
        f"# {int(np.sum(held & (r == R_RANGE[0])))} of {held.size} points "
        f"are at r's floor, {R_RANGE[0]:g} (the kinetic tail colder than the "
        f"form follows), {int(np.sum(held & (r == R_RANGE[1])))} at its top, "
        f"{R_RANGE[1]:g} (hotter).",
    ]
    for k, xs in enumerate(xs_grid):
        lines.append(f"XS {xs:g}")
        lines.append("TE " + " ".join(f"{t:g}" for t in te_grid))
        for i, x in enumerate(x_grid):
            lines.append(f"{x:.3f} " + " ".join(f"{v:.5f}"
                                                for v in r[i, :, k]))
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    workers = int(argv[0]) if argv else 1
    path = xsec.DATA_DIR / xsec.TAIL_FIELD_FILE[0]
    t0 = time.time()
    r, held, _ = build(workers=workers, verbose=True)
    write(path, r, held)
    print(f"wrote {path} ({time.time() - t0:.0f} s), md5 "
          f"{xsec.md5_of(path)}: pin it in xsec.TAIL_FIELD_FILE")


if __name__ == "__main__":
    sys.exit(main())
