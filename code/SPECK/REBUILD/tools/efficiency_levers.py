#!/usr/bin/env python3
"""Which operating levers raise the KrF intrinsic efficiency? — one-at-a-time
scan around the HYPER reference point, then the best-found combination.

Owner ask (2026-08-21): "figure out a way to increase the efficiency."
Beam-only runs (no microwave) on krf_hyper_mw; every lever varied alone
around the reference (2 atm, 100 A cm^-2, 60 ns flat, 10 cm cube, R_OC
0.95, L_op 0.02, Ar/Kr/F2 0.949/0.05/0.001), reporting intrinsic
efficiency (laser / deposited), laser energy, deposited energy, peak g0,
threshold margin g0/g_th, and F2 burn. Then the top levers are combined.

Levers and their honesty limits:
  pressure          0.5 - 3 atm (DE-AR range)
  pump density      j 50 - 400 A cm^-2 (600 kV)
  pulse length      flat 30 - 150 ns (longer extraction window)
  gain length       Lz 0.1 - 0.3 m (optical axis; transverse 10 x 10 cm
                    kept; beam along x so the beam area grows with Lz —
                    same current density, a longer cell)
  optics            R_OC 0.90 - 0.99, L_op 0.005 - 0.05
  F2 fraction       0.05 - 0.3 % (Ar renormalized) — the electron-impact
                    tables are the 0.1 % fill's: attachment scales with
                    n_F2 linearly, the EEDF itself is Ar-dominated, so a
                    x3 bracket is defensible; stated
  Kr fraction       2.5 - 10 % (Ar renormalized) — same caveat, weaker:
                    Kr shapes the EEDF more; read as a bracket, not a curve

Predictions scored against nothing; the engine's anchored formation
deficit (x0.8-0.9) makes RELATIVE slopes the robust output.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/efficiency_levers.py
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.MODULES.engine import Device  # noqa: E402
from REBUILD.MODULES.params import engine_parameters  # noqa: E402
from REBUILD.run import PRESETS  # noqa: E402

REF = {"mw_power_W": 0.0, "t_end_ns": 250.0}
LEVERS = [
    ("pressure_atm", [0.5, 1.0, 1.5, 2.0, 3.0]),
    ("beam_j_A_cm2", [50.0, 100.0, 200.0, 400.0]),
    ("beam_flat_ns", [30.0, 60.0, 100.0, 150.0]),
    ("Lz_m", [0.1, 0.2, 0.3]),
    ("r_oc", [0.90, 0.95, 0.98, 0.99]),
    ("l_op", [0.005, 0.01, 0.02, 0.05]),
    ("F2_frac", [0.0005, 0.001, 0.002, 0.003]),
    ("Kr_frac", [0.025, 0.05, 0.10]),
]


def run_point(overrides, mixture=None):
    spec = PRESETS["krf_hyper_mw"]
    P = engine_parameters()
    P.set_many(spec["params"], source="problem")
    ov = dict(REF)
    ov.update(overrides)
    if "beam_flat_ns" in ov and ov["beam_flat_ns"] > 150.0:
        ov["t_end_ns"] = ov["beam_flat_ns"] + 150.0
    if ov.get("beam_flat_ns", 60.0) > 80.0:
        ov["t_end_ns"] = ov["beam_flat_ns"] + 150.0
    P.set_many(ov, source="cli")
    dev = Device(P, mixture=mixture)
    res = dev.run(n_save=121)
    if not res["success"]:
        return None
    g0 = float(np.max(res["g0_m1"]))
    cav = dev.cav
    g_th = (-np.log(cav.R_OC) - np.log(1.0 - cav.single_pass_loss(0.0))) \
        / (2.0 * cav.l_g_m)
    fin = res["final_channels_J"]
    return {"eta": res["eta_intrinsic"], "E_laser_J": fin["laser_out"],
            "E_dep_J": fin["drive_in"], "g0_peak_m1": g0, "g_th_m1": g_th,
            "margin": g0 / g_th if g_th > 0 else float("inf"),
            "f2_burn": res["validity"]["f2_burn_frac"],
            "regime": res["regime"]}


def mixture_for(f2=None, kr=None):
    base = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
    if f2 is not None:
        base["F2"] = f2
    if kr is not None:
        base["Kr"] = kr
    base["Ar"] = 1.0 - base["Kr"] - base["F2"]
    return base


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--levers", default=",".join(k for k, _ in LEVERS))
    ap.add_argument("--no-combo", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    want = {s.strip() for s in args.levers.split(",") if s.strip()}
    rows = []
    hdr = (f"{'lever':>13} {'value':>8} {'eta%':>7} {'E_las J':>9} "
           f"{'E_dep J':>8} {'g0/m':>6} {'g0/gth':>7} {'F2burn%':>8}")
    print(hdr)
    print("-" * len(hdr))
    ref = run_point({})
    print(f"{'reference':>13} {'':>8} {100 * ref['eta']:7.3f} "
          f"{ref['E_laser_J']:9.4f} {ref['E_dep_J']:8.2f} "
          f"{ref['g0_peak_m1']:6.2f} {ref['margin']:7.2f} "
          f"{100 * ref['f2_burn']:8.1f}", flush=True)
    rows.append({"lever": "reference", "value": None, **ref})
    best = {}
    for lever, values in LEVERS:
        if lever not in want:
            continue
        for val in values:
            if lever == "F2_frac":
                r = run_point({}, mixture=mixture_for(f2=val))
            elif lever == "Kr_frac":
                r = run_point({}, mixture=mixture_for(kr=val))
            else:
                r = run_point({lever: val})
            if r is None:
                print(f"{lever:>13} {val:8g}   FAILED", flush=True)
                continue
            rows.append({"lever": lever, "value": val, **r})
            print(f"{lever:>13} {val:8g} {100 * r['eta']:7.3f} "
                  f"{r['E_laser_J']:9.4f} {r['E_dep_J']:8.2f} "
                  f"{r['g0_peak_m1']:6.2f} {r['margin']:7.2f} "
                  f"{100 * r['f2_burn']:8.1f}", flush=True)
            if r["eta"] > best.get(lever, (None, -1))[1]:
                best[lever] = (val, r["eta"])
    combo = None
    if not args.no_combo:
        # combine each lever's best value (only those that beat the
        # reference) into one run — a greedy combination, stated as such
        ov, mix_kw = {}, {}
        for lever, (val, eta) in best.items():
            if eta <= ref["eta"] * 1.001:
                continue
            if lever == "F2_frac":
                mix_kw["f2"] = val
            elif lever == "Kr_frac":
                mix_kw["kr"] = val
            else:
                ov[lever] = val
        combo_mix = mixture_for(**mix_kw) if mix_kw else None
        r = run_point(ov, mixture=combo_mix)
        combo = {"overrides": ov, "mixture": combo_mix, **(r or {})}
        if r:
            print(f"\nGREEDY COMBINATION of the winning levers {ov} "
                  f"{'+ mixture ' + str(combo_mix) if combo_mix else ''}:"
                  f"\n  eta {100 * r['eta']:.3f}% (reference "
                  f"{100 * ref['eta']:.3f}%), laser {r['E_laser_J']:.4f} J, "
                  f"deposited {r['E_dep_J']:.2f} J, g0/g_th {r['margin']:.2f}")
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "efficiency")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "efficiency_levers.json"
    path.write_text(json.dumps({"reference": ref, "rows": rows,
                                "best_per_lever": best, "combo": combo},
                               indent=2, default=float))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
