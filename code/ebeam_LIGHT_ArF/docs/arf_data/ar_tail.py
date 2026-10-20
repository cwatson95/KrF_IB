"""The Ar two-temperature tail table for LIGHT_ArF v0_2 (ML's method on Ar).

    python3 LIGHT_ArF_ebeam/docs/arf_data/ar_tail.py            # the whole table
    python3 LIGHT_ArF_ebeam/docs/arf_data/ar_tail.py --quick    # three states

LIGHT_ArF's basic network integrates its Ar excitation and ionization rows
(every IST-Lisbon Ar block) over a two-temperature EEDF: the Maxwellian at
Te below the lowest Ar threshold E_b = 11.548 eV, and above it a tail that
falls with temperature T_t = r Te, continuous at E_b (v0_2/xsec.py). This
script makes r, the table v0_2/data/ar_tail_ratio.txt, exactly the way
ML/ml/tail.py makes KrF's (LIGHT/v0/data/kr_tail_ratio.txt):

  1. the state. LIGHT_ArF v0_2 on its Maxwellian closure (eedf maxwellian,
     (3/2)Te per captured electron) under a prototype e-beam drive at one
     pump and F2 fraction, 1 atm, the 200-400 ns window of a 450-ns pulse.
     The drive (BeamModel) is the survey's KrF drive
     (LIGHT/docs/survey_2026-09-27/scripts/ebeam_screen.py) with REBUILD's
     Spencer-Fano per-pair yields for Ar/0.5 % F2 in place of Kr's. It is
     part of this generator only, never of LIGHT_ArF (LIGHT has no drive).
  2. the EEDF. ML's beam-sustained Boltzmann solution with e-e collisions
     (ML/ml/boltzmann/ee.py solve_beam, imported read-only) at that state:
     electrons injected at the pair rate with the sub-excitation spectrum
     on [0, E_b] (mean EPS_SEC), colliding with Ar and F2 (LoKI-MC's
     IST-Lisbon and SIGLO files, the ones LIGHT_ArF vendors) and with one
     another, and lost by attachment and by recombination (LIGHT_ArF's R11p
     law as the sink).
  3. r. The two-temperature EEDF at Te = 2/3 <eps> removes the same energy
     per electron through the Ar tail rows as the kinetic EEDF does.
The columns are the fill F2 fractions, the rows log10(n_e/N); each column is
one series of states, linear in (log10 n_e/N, ln r) between them and held at
its ends. Nothing is fitted, and no scored output enters.

KM_ROOT (default: the directory holding LIGHT_ArF) locates ML/ and
LoKI-MC/. Nothing there is written.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import math
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve()
LIGHT_ARF = next(p for p in HERE.parents if (p / "run_v0_2.py").is_file())
KM = Path(os.environ.get("KM_ROOT", LIGHT_ARF.parent)).expanduser().resolve()
BOLTZ = KM / "ML" / "ml" / "boltzmann"
INP = KM / "LoKI-MC" / "Code" / "Input"
if not (BOLTZ / "ee.py").is_file() or not INP.is_dir():
    raise SystemExit(f"KM_ROOT={KM}: needs ML/ml/boltzmann/ee.py and "
                     "LoKI-MC/Code/Input (set KM_ROOT to Kinetic_Modeling)")
sys.path.insert(0, str(LIGHT_ARF))
sys.path.insert(0, str(BOLTZ))          # the vendored solver imports flat

import ee  # noqa: E402
import gases_twin as GT  # noqa: E402

from v0_2 import network as NET  # noqa: E402
from v0_2 import xsec  # noqa: E402
from v0_2.model import QE, Params  # noqa: E402
from v0_2.species import RG  # noqa: E402

GT.load_tables(str(INP))

# --- the prototype beam drive (beam_drive.py, shared with score_arf.py) ----
sys.path.insert(0, str(HERE.parent))
from beam_drive import (E_B, EPS_SEC, W_EFF, BeamModel)  # noqa: E402


WINDOW = (200e-9, 400e-9)
PULSE = 450e-9


def run_state(x_f2, pump, torr=760.0):
    """The drive's window state: mean densities, event rates, n_e, Te, Tg."""
    p = Params(pressure_atm=torr / 760.0, x_f2=x_f2, power_W_cm3=pump,
               pulse_s=PULSE, t_end_s=PULSE, ne_cm3=1e10, n_save=1201,
               eedf="maxwellian", capture="three_halves")
    m = BeamModel(p)
    r = m.run()
    t = r.t
    win = (t >= WINDOW[0]) & (t <= WINDOW[1])
    return {"n": {s: float(np.mean(r.n[s][win])) for s in r.n},
            "events": {c.row.rid: float(np.polyfit(
                t[win], r.events[c.row.rid][win], 1)[0])
                for c in m.compiled},
            "ne": float(np.mean(r.ne[win])), "te": float(np.mean(r.te[win])),
            "tg": float(np.mean(r.tg[win])),
            "pairs": pump / QE / W_EFF, "model": m}


# --- the kinetic EEDF (ML/ml/eedf.py beam_eedf, on Ar) -----------------------
_GASES: dict = {}


def gases(tg):
    key = round(tg, 1)
    if key not in _GASES:
        _GASES[key] = (GT.build("Ar", tg, str(INP)),
                       GT.build("F2", tg, str(INP)))
    return _GASES[key]


def tail_rows():
    """LIGHT_ArF's tail rows: (row id, threshold) of every Ar excitation and
    the ionization, read from the network (xsec.in_tail)."""
    out = []
    for r in NET.build():
        if r.rate[0] == "xsec" and xsec.in_tail(*r.rate[1:4]):
            out.append((r.rid, xsec.find(*r.rate[1:4]).threshold))
    return out


def collision(gas, thr):
    hits = [c for c in gas["colls"] if abs(c["thr"] - thr) < 1e-6]
    if len(hits) != 1:
        raise KeyError(f"threshold {thr}: {len(hits)} collisions")
    return hits[0]


def r11p_pow() -> float:
    row = {r.rid: r for r in NET.build()}["R11p"]
    if row.rate[0] != "te_pow":
        raise ValueError(f"R11p: rate kind {row.rate[0]!r}")
    return float(row.rate[3])


def beam_eedf(dev, N=600, emax=30.0, tol=1e-7):
    n, tg = dev["n"], dev["tg"]
    ar, f2 = gases(tg)
    n_gas = (n[RG] + n["F2"] + n["F"]) * 1e6
    mix = [(ar, n[RG] * 1e6 / n_gas), (f2, n["F2"] * 1e6 / n_gas)]
    nu_rec = dev["events"]["R11p"] / dev["ne"]
    sol = ee.solve_beam(mix, tg, n_gas, N, emax, dev["pairs"] * 1e6,
                        lambda e: ee.injection(e, E_B, EPS_SEC), nu_rec,
                        dev["ne"] * 1e6, ee_on=True, rec_te=dev["te"],
                        tol=tol, rec_pow=r11p_pow())
    return {"ne": sol["ne"] / 1e6, "mean": sol["mean"], "sol": sol,
            "n_gas": n_gas / 1e6}


def two_t(ec, du, te, tt, e_b=E_B):
    f = (np.exp(-np.minimum(ec, e_b) / te)
         * np.exp(-np.maximum(ec - e_b, 0.0) / tt))
    return f / float(np.sum(np.sqrt(ec) * f) * du)


def tail_loss(F, ec, du, ar, rows):
    """Energy per electron per gas particle into the tail rows."""
    s = dict(F=F, ec=ec, du=du)
    return sum(ee.rate(collision(ar, thr), s) * thr for _, thr in rows)


R_LO, R_HI = 0.002, 2.0


def ratio(b, tg, rows):
    sol = b["sol"]
    ar, _ = gases(tg)
    ec, du = sol["ec"], sol["du"]
    te = 2.0 / 3.0 * b["mean"]
    target = tail_loss(sol["F"], ec, du, ar, rows)

    def gap(lr):
        return math.log(tail_loss(two_t(ec, du, te, te * math.exp(lr)), ec,
                                  du, ar, rows) / target)
    lo, hi = math.log(R_LO), math.log(R_HI)
    if gap(lo) > 0.0:
        return R_LO
    if gap(hi) < 0.0:
        return R_HI
    return math.exp(brentq(gap, lo, hi, xtol=1e-7))


class _Late(Exception):
    pass


def _alarm(signum, frame):
    raise _Late()


def point(x_f2, pump, torr=760.0, rows=None, alarm_s=300):
    rows = rows or tail_rows()
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(alarm_s)
    try:
        dev = run_state(x_f2, pump, torr)
        b = beam_eedf(dev)
    except (_Late, RuntimeError, ValueError) as exc:
        return {"fail": f"{type(exc).__name__}: {exc}"}
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    r = ratio(b, dev["tg"], rows)
    n_gas = dev["n"][RG] + dev["n"]["F2"] + dev["n"]["F"]
    return {"torr": torr, "x_f2": x_f2, "pump": pump,
            "x_win": dev["n"]["F2"] / n_gas, "ne_N": b["ne"] / n_gas,
            "te": 2.0 / 3.0 * b["mean"], "r": r, "tg": dev["tg"],
            "te_light": dev["te"], "ne_light": dev["ne"]}


# --- the table -----------------------------------------------------------------
X_F2 = (0.0, 0.0003, 0.001, 0.002, 0.004, 0.008, 0.016, 0.03)
PUMPS = (3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6)
PRESSURE_CHECK = ((380.0, 0.004, 3e4), (380.0, 0.004, 3e5),
                  (1520.0, 0.004, 3e4), (1520.0, 0.004, 3e5))
GRID_X = np.round(np.arange(-8.5, -2.99, 0.125), 3)
PATH = LIGHT_ARF / "v0_2" / "data" / "ar_tail_ratio.txt"
STATES = HERE.parent / "ar_tail_states.csv"


def series(pts, x_f2, torr=760.0):
    return np.array(sorted((math.log10(p["ne_N"]), math.log(p["r"]))
                           for p in pts if p["torr"] == torr
                           and p["x_f2"] == x_f2 and R_LO < p["r"] < R_HI))


def table(pts):
    cols = []
    for x in X_F2:
        s = series(pts, x)
        if len(s) == 0:
            raise RuntimeError(f"F2 {x}: no solved state in the column")
        cols.append(np.exp(np.interp(GRID_X, s[:, 0], s[:, 1])))
    return np.stack(cols, axis=1)


def lookup(r, ne_N, x_f2):
    """v0_2's bilinear read (xsec.TailRatio), on this grid."""
    tr = xsec.TailRatio(GRID_X, X_F2, r, "")
    return tr(ne_N, x_f2)


def write_table(pts, path=PATH):
    r = table(pts)
    used = [p for p in pts if p["torr"] == 760.0]
    dense = [p for p in used if p["ne_N"] >= 1e-6]
    back = np.log([lookup(r, p["ne_N"], p["x_win"]) / p["r"] for p in dense])
    w = dense[int(np.argmax(np.abs(back)))]
    worst = (f"{w['pump']:g} W/cm3 at {100 * w['x_f2']:g} % F2, "
             f"{100 * (1 - w['x_win'] / w['x_f2']) if w['x_f2'] else 0:.0f} "
             "% of it burned")
    check = [(p["torr"], lookup(r, p["ne_N"], p["x_win"]) / p["r"])
             for p in pts if p["torr"] != 760.0]
    lines = [
        "# LIGHT_ArF's two-temperature tail: r = T_tail / Te for the Ar "
        "excitation and ionization rows",
        "# (LIGHT_ArF_ebeam/v0_2/xsec.py). Rows: log10(n_e/N). Columns: the F2 "
        "fraction [F2]/N.",
        f"# Generated by LIGHT_ArF_ebeam/docs/arf_data/ar_tail.py on "
        f"{datetime.date.today().isoformat()} (ML/ml/tail.py's method on Ar). "
        "Do not edit: regenerate.",
        "# Each state: LIGHT_ArF v0_2 basic network on its Maxwellian "
        "closure under the prototype ArF",
        "# e-beam drive at 1 atm, 200-400 ns window. Its kinetic EEDF "
        "(ML ee.solve_beam, e-e on) sets r:",
        f"# the two-temperature EEDF at Te = 2/3 <eps> removes the same "
        f"energy through the {len(tail_rows())} Ar tail rows.",
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
        "kinetic r = " + ", ".join(f"{q:.3f}" for _, q in check) + ".",
        "F2 " + " ".join(f"{x:g}" for x in X_F2),
    ]
    for u, row in zip(GRID_X, r):
        lines.append(f"{u:.3f} " + " ".join(f"{v:.5f}" for v in row))
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    return path, hashlib.md5(text.encode()).hexdigest(), check


def write_states(pts, path=STATES):
    cols = ("torr", "x_f2", "x_win", "pump", "ne_N", "te", "r", "tg",
            "te_light", "ne_light")
    lines = [",".join(cols)]
    for p in pts:
        lines.append(",".join(f"{p[c]:.6g}" for c in cols))
    path.write_text("\n".join(lines) + "\n")
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="three states only; writes nothing")
    a = ap.parse_args(argv)
    rows = tail_rows()
    print(f"W_eff {W_EFF:.3f} eV/pair, E_b {E_B:.3f} eV, {len(rows)} tail "
          f"rows, R11p Te^{r11p_pow():g}")
    states = ([(760.0, 0.004, 3e4), (760.0, 0.004, 3e5), (760.0, 0.0, 3e3)]
              if a.quick else
              [(760.0, x, p) for x in X_F2 for p in PUMPS]
              + list(PRESSURE_CHECK))
    pts, t0 = [], time.time()
    for torr, x, pump in states:
        o = point(x, pump, torr, rows)
        if "fail" in o:
            print(f"  {torr:6.0f} Torr  F2 {x:6.4f}  {pump:8.1e} W/cm3  "
                  f"no solution ({o['fail']})", flush=True)
            continue
        pts.append(o)
        print(f"  {torr:6.0f} Torr  F2 {x:6.4f}  {pump:8.1e} W/cm3  n_e/N "
              f"{o['ne_N']:8.2e}  Te {o['te']:5.3f}  r {o['r']:6.4f}  "
              f"(LIGHT's Maxwellian state: Te {o['te_light']:5.3f} eV, "
              f"n_e {o['ne_light']:8.2e} cm-3)  "
              f"{time.time() - t0:6.0f} s", flush=True)
    if a.quick:
        return 0
    path, md5, check = write_table(pts)
    write_states(pts)
    print(f"wrote {path} md5 {md5}; pressure check {check}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
