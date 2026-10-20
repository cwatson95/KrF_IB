"""The two-temperature tail: the kinetic EEDF's tail in the compact form
LIGHT uses (step 3 of the owner's plan, 2026-09-28).

LIGHT's basic network integrates its Kr excitation and ionization rows
(G1, C1a, C1b, C1c, C2) over a two-temperature EEDF: the Maxwellian at Te
below the lowest Kr threshold E_b = 9.915 eV, and above it a tail that
falls with temperature T_t = r Te, continuous at E_b (LIGHT/v0/xsec.py).
This module makes r, a table in (log10 n_e/N, Te) for LIGHT/v0/data/.

Each point of the table rests on one kinetic solution:
  1. the state. LIGHT under the prototype beam drive at one pump rate and
     F2 fraction gives the gas, the pair rate and the recombination
     frequency of the 200-400 ns window. It runs on its Maxwellian closure
     (forward.MAXWELL), so the table never depends on itself, and stops
     at the end of the pulse (the afterglow is not needed, and at low pump
     with F2 attachment drives n_e to LIGHT's floor there).
  2. the EEDF. ml/eedf.py's beam-sustained Boltzmann solution with e-e
     collisions at that state gives the EEDF, its n_e and its mean energy.
  3. r. The two-temperature EEDF at Te = 2/3 <eps> (LIGHT's convention,
     U_e = 3/2 n_e Te) removes the same energy per electron through the
     five Kr rows as the kinetic EEDF does. G1 (Kr 5s) carries nearly all
     of it and is matched to 3 %. The minor rows are not matched one by
     one (C1a x0.2-1.2, C2 x30-800), and each is below 3 % of G1 there.
The table's axes are the states' own: n_e/N (the e-e collisions that
refill the tail) and the F2 fraction (how fast F2 slows the injected
electrons below the break). Te does not carry it. Without F2 the injected
electrons are slowed only by elastic and e-e collisions, and e-e lifts
them over the break: pure Kr at 300 W/cm3 has r 0.27 where 0.1 % F2 at
the same n_e/N has 0.035. A polynomial surface in (n_e/N, Te) fitted to
the states put r = 1 in corners no state reaches (2026-09-28). Each
column is one series of states, and nothing is fitted.

The scored outputs never enter. The table follows from the cross sections
and the collision physics alone, as ml/eedf.py does.
"""
from __future__ import annotations

import datetime
import hashlib
import math
import signal

import numpy as np
from scipy.optimize import brentq

from . import ROOT, eedf, forward, score
import ee  # noqa: E402  (ml.eedf puts the vendored solver on sys.path)

TAIL = ("G1", "C1a", "C1b", "C1c", "C2")
E_B = eedf.E_TOP                       # 9.915 eV
R_LO, R_HI = 0.002, 2.0                # search bounds for r
PATH = ROOT / "LIGHT" / "v0" / "data" / "kr_tail_ratio.txt"
# the states: F2 fraction x pump [W/cm3] at 1 atm, and a pressure check
X_F2 = (0.0, 0.0003, 0.001, 0.002, 0.004, 0.008, 0.016, 0.03)
PUMPS = (3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6)
PRESSURE_CHECK = ((380.0, 0.004, 3e4), (380.0, 0.004, 3e5),
                  (1520.0, 0.004, 3e4), (1520.0, 0.004, 3e5))
# the table grid
GRID_X = np.round(np.arange(-8.5, -2.99, 0.125), 3)      # log10(n_e/N)


def two_t(ec, du, te, tt, e_b=E_B):
    """The two-temperature F0 on the solver's grid (int sqrt(e) F0 = 1)."""
    f = (np.exp(-np.minimum(ec, e_b) / te)
         * np.exp(-np.maximum(ec - e_b, 0.0) / tt))
    return f / float(np.sum(np.sqrt(ec) * f) * du)


def tail_rates(F, ec, du, kr, f2):
    s = dict(F=F, ec=ec, du=du)
    return {r: ee.rate(eedf.collision(kr, f2, *eedf.ROWS[r]), s) for r in TAIL}


def tail_loss(k):
    """Energy per electron per gas particle into the Kr rows [eV m3/s]."""
    return sum(k[r] * eedf.ROWS[r][1] for r in TAIL)


def _kinetic(b, tg):
    sol = b["sol"]
    kr, f2 = eedf.gases(tg)
    return (sol["ec"], sol["du"], kr, f2, 2.0 / 3.0 * b["mean"],
            tail_rates(sol["F"], sol["ec"], sol["du"], kr, f2))


def rows_at(b, tg, r):
    """Each tail row's two-temperature rate at ratio r over its kinetic
    rate, for one kinetic solution (eedf.beam_eedf output)."""
    ec, du, kr, f2, te, kin = _kinetic(b, tg)
    k2 = tail_rates(two_t(ec, du, te, te * r), ec, du, kr, f2)
    return {row: k2[row] / kin[row] for row in TAIL}


def ratio(b, tg):
    """r for one kinetic solution: the two-temperature EEDF at Te = 2/3
    <eps> removes the kinetic energy through the tail rows. Returns r and
    rows_at(b, tg, r)."""
    ec, du, kr, f2, te, kin = _kinetic(b, tg)
    target = tail_loss(kin)

    def gap(lr):
        k = tail_rates(two_t(ec, du, te, te * math.exp(lr)), ec, du, kr, f2)
        return math.log(tail_loss(k) / target)
    lo, hi = math.log(R_LO), math.log(R_HI)
    if gap(lo) > 0.0:
        r = R_LO
    elif gap(hi) < 0.0:
        r = R_HI
    else:
        r = math.exp(brentq(gap, lo, hi, xtol=1e-7))
    return r, rows_at(b, tg, r)


class _Late(Exception):
    pass


def _alarm(signum, frame):
    raise _Late()


def point(c, x_f2, pump, torr=760.0, alarm_s=240):
    """One state: LIGHT's window state, its kinetic EEDF and r. None if the
    run fails or takes longer than alarm_s."""
    cond = {"pressure_torr": torr, "f2_torr": x_f2 * torr,
            "power_W_cm3": pump, "pulse_s": 450e-9}
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(alarm_s)
    try:
        dev = forward.run_device(c, cond, closure=forward.MAXWELL,
                                 afterglow_s=0.0)
        b = eedf.beam_eedf(dev, tol=1e-7)
    except (_Late, RuntimeError):
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    r, rows = ratio(b, dev["tg"])
    n_gas = dev["n"]["Kr"] + dev["n"]["F2"] + dev["n"]["F"]
    return {"torr": torr, "x_f2": x_f2, "pump": pump,
            "x_win": dev["n"]["F2"] / n_gas,
            "ne_N": b["ne"] / n_gas, "te": 2.0 / 3.0 * b["mean"], "r": r,
            "rows": rows, "tg": dev["tg"], "f_G1": b["f"]["G1"],
            "beam": b, "dev": dev}


def scan(c=None, states=None, echo=print):
    """Every state of the grid (X_F2 x PUMPS at 1 atm) plus the pressure
    check. Returns the points that ran."""
    c = c or score.candidate_sets()["LIGHT now"]
    states = states or ([(760.0, x, p) for x in X_F2 for p in PUMPS]
                        + list(PRESSURE_CHECK))
    pts = []
    for torr, x, pump in states:
        o = point(c, x, pump, torr)
        if o is None:
            echo(f"  {torr:6.0f} Torr  F2 {x:6.4f}  {pump:8.1e} W/cm3  "
                 "no solution")
            continue
        pts.append(o)
        echo(f"  {torr:6.0f} Torr  F2 {x:6.4f}  {pump:8.1e} W/cm3  n_e/N "
             f"{o['ne_N']:8.2e}  Te {o['te']:5.3f}  r {o['r']:6.4f}  "
             f"G1 {o['rows']['G1']:5.3f}  C1a {o['rows']['C1a']:6.3f}")
    return pts


def series(pts, x_f2, torr=760.0):
    """One solved series at a fill F2 fraction: rows (log10 n_e/N, ln r),
    in order of n_e/N."""
    return np.array(sorted((math.log10(p["ne_N"]), math.log(p["r"]))
                           for p in pts if p["torr"] == torr
                           and p["x_f2"] == x_f2 and R_LO < p["r"] < R_HI))


def table(pts, grid=GRID_X, columns=X_F2):
    """r on the grid: each column one series, linear in (log10 n_e/N,
    ln r) between its states and held at its ends."""
    cols = []
    for x in columns:
        s = series(pts, x)
        cols.append(np.exp(np.interp(grid, s[:, 0], s[:, 1])))
    return np.stack(cols, axis=1)


def lookup(r, ne_N, x_f2, grid=GRID_X, columns=X_F2):
    """LIGHT's bilinear read of the table (the same rule as
    LIGHT/v0/xsec.TailRatio)."""
    u = min(max(math.log10(ne_N), grid[0]), grid[-1])
    f = min(max(x_f2, columns[0]), columns[-1])
    i = min(max(int(np.searchsorted(grid, u)) - 1, 0), len(grid) - 2)
    j = min(max(int(np.searchsorted(columns, f)) - 1, 0), len(columns) - 2)
    a = (u - grid[i]) / (grid[i + 1] - grid[i])
    b = (f - columns[j]) / (columns[j + 1] - columns[j])
    return float((1 - a) * (1 - b) * r[i, j] + a * (1 - b) * r[i + 1, j]
                 + (1 - a) * b * r[i, j + 1] + a * b * r[i + 1, j + 1])


def readback(pts, r):
    """Each state's table r, read as LIGHT reads it (its n_e/N and its
    window's F2, which burns down at high pump) over its kinetic r."""
    return [lookup(r, p["ne_N"], p["x_win"]) / p["r"] for p in pts]


def write_table(pts, path=PATH):
    """Write the table. Returns (path, md5, the pressure check as
    [(torr, table r / kinetic r)]). The table comes from v0's states, and
    v0_1 carries a pinned copy of it, so only a v0 run may write it."""
    from . import LIGHT_VERSION
    if LIGHT_VERSION != "v0" and path == PATH:
        raise RuntimeError("the tail table is made from v0's states; "
                           f"refusing to write it from {LIGHT_VERSION}")
    r = table(pts)
    used = [p for p in pts if p["torr"] == 760.0]
    dense = [p for p in used if p["ne_N"] >= 1e-6]
    back = np.log(readback(dense, r))
    w = dense[int(np.argmax(np.abs(back)))]
    worst = (f"{w['pump']:g} W/cm3 at {100 * w['x_f2']:g} % F2, "
             f"{100 * (1 - w['x_win'] / w['x_f2']) if w['x_f2'] else 0:.0f} "
             "% of it burned")
    check = [(p["torr"], lookup(r, p["ne_N"], p["x_win"]) / p["r"])
             for p in pts if p["torr"] != 760.0]
    lines = [
        "# LIGHT's two-temperature tail: r = T_tail / Te for the Kr "
        "excitation and ionization rows",
        "# (LIGHT/v0/xsec.py). Rows: log10(n_e/N). Columns: the F2 "
        "fraction [F2]/N.",
        f"# Generated by ML/ml/tail.py (python3 ML/run.py tail) on "
        f"{datetime.date.today().isoformat()}. Do not edit: regenerate.",
        "# Each state: LIGHT basic network + ions + f2e + kr2 on its "
        "Maxwellian closure under the prototype",
        "# e-beam drive at 1 atm, 200-400 ns window. Its kinetic EEDF "
        "(ml/eedf.py, e-e collisions on) sets r:",
        "# the two-temperature EEDF at Te = 2/3 <eps> removes the same "
        "energy through G1, C1a, C1b, C1c, C2.",
        f"# Each column is one series of states (pump {PUMPS[0]:g}-"
        f"{PUMPS[-1]:g} W/cm3, {len(used)} solved in all), linear",
        "# in (log10 n_e/N, ln r) between them and held at its ends. "
        "Columns are the fill F2. Read at each",
        "# state's own n_e/N and window F2 (burned down at high pump), the "
        "table gives its kinetic r to",
        f"# {100 * float(np.sqrt(np.mean(back ** 2))):.0f} % rms where "
        f"n_e/N >= 1e-6 (largest {100 * float(np.max(np.abs(back))):.0f} "
        f"%: {worst}).",
        "# Pressure check (0.5 and 2 atm, not in the table): table r / "
        "kinetic r = "
        + ", ".join(f"{q:.3f}" for _, q in check) + ".",
        "F2 " + " ".join(f"{x:g}" for x in X_F2),
    ]
    for u, row in zip(GRID_X, r):
        lines.append(f"{u:.3f} " + " ".join(f"{v:.5f}" for v in row))
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    return path, hashlib.md5(text.encode()).hexdigest(), check


def write_states(pts, path):
    """The kinetic states behind the table, one per line (a results file:
    run.py tail rewrites it)."""
    cols = ("torr", "x_f2", "x_win", "pump", "ne_N", "te", "r", "tg", "f_G1")
    lines = [",".join(cols + tuple(f"{r}_2T_over_kinetic" for r in TAIL))]
    for p in pts:
        lines.append(",".join([f"{p[c]:.6g}" for c in cols]
                              + [f"{p['rows'][r]:.6g}" for r in TAIL]))
    path.write_text("\n".join(lines) + "\n")
    return path
