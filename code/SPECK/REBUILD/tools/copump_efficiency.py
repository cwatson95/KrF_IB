#!/usr/bin/env python3
"""Efficiency range WITH co-pumping — post-processing of the co-pump map.

Owner ask (2026-08-21): "what is the range of efficiency with this
co-pumping?" Reads rebuild_out/copump_scan/copump_scan_report.json (every
point carries the baseline laser energy, the deposited beam energy, the
absorbed microwave energy and dE_laser) and reports, per point:

  eta_dep0   = E_laser0 / E_beam_dep                 intrinsic, beam only
  eta_dep    = (E_laser0 + dE) / (E_beam_dep + E_MW) intrinsic, co-pumped
               (the DE-AR efficiency-metric note: microwave input
               SUBTRACTED — it sits in the denominator)
  eta_inj0/eta_inj   the same per INJECTED beam energy (600 kV x j x
               aperture x pulse) — the 10 cm cell passes most of a
               600 kV beam, so this is the scale a beam operator sees;
               microwave counted at the coupled (absorbed) value
  P_+1pp     the coupled microwave power that would lift eta_dep by ONE
               percentage point at this point's marginal m and window:
               from (E_L0 + m E)/(E_d + E) = eta0 + d  =>
               E = d E_d / (m - eta0 - d), P = E / t_window
               (m <= eta0 + d means no power can do it — reported inf)

What is NOT in these numbers (stated): beam-diode / pulser / hibachi
wall-plug factors and the gyrotron's own efficiency (DE-AR: 2 MW source
-> ~250 kW coupled, ~12.5 %); the intrinsic numbers here are per joule
DEPOSITED in the gas, which is where the two pumps can be compared
like-for-like. Predictions, scored against nothing.

Usage (from the Kinetic_Modeling root, after tools/copump_scan.py):
    python3 SPECK/REBUILD/tools/copump_efficiency.py
"""

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

BEAM_KV = 600.0
APERTURE_M2 = 0.1 * 0.1           # krf_hyper_mw: 10 cm cube, beam along x
PULSE_S = (60.0 + 5.0) * 1e-9      # flat + rise (trapezoid area), preset
REF_J = 100.0


def e_injected_J(j_A_cm2):
    return BEAM_KV * 1e3 * j_A_cm2 * 1e4 * APERTURE_M2 * PULSE_S


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    path = (Path(args.report) if args.report
            else ROOT / "rebuild_out" / "copump_scan" / "copump_scan_report.json")
    d = json.loads(path.read_text())
    rows = []
    hdr = (f"{'axis':>17} {'point':>26} {'eta_dep0%':>9} {'eta_dep%':>8} "
           f"{'d_pp':>7} {'eta_inj0%':>9} {'eta_inj%':>8} {'marg%':>6} "
           f"{'P_+1pp MW':>10}  note")
    print(hdr)
    print("-" * len(hdr))
    for r in d["points"]:
        e_d = r["E_beam_dep_J"]
        e_l0 = r["E_laser_baseline_J"]
        e_mw = r["E_MW_J"]
        d_l = r["dE_laser_J"]
        m = r["marginal_laser"]
        j = r["beam_overrides"].get("beam_j_A_cm2", REF_J)
        e_inj = e_injected_J(j)
        eta0 = e_l0 / e_d
        eta = (e_l0 + d_l) / (e_d + e_mw)
        eta_i0 = e_l0 / e_inj
        eta_i = (e_l0 + d_l) / (e_inj + e_mw)
        win = r["window"]
        t_win = (win["mw_flat_ns"] + win["mw_rise_ns"]) * 1e-9
        dpp = 0.01
        denom = m - eta0 - dpp
        p_1pp = (dpp * e_d / denom / t_win) if denom > 0 else float("inf")
        note = "" if r["baseline_lases"] else "BASELINE DOES NOT LASE"
        rows.append({**{k: r[k] for k in ("axis", "point", "P_mw_W")},
                     "eta_dep_baseline": eta0, "eta_dep_copump": eta,
                     "delta_pp": 100 * (eta - eta0),
                     "eta_inj_baseline": eta_i0, "eta_inj_copump": eta_i,
                     "E_injected_beam_J": e_inj, "E_beam_dep_J": e_d,
                     "P_plus_1pp_W": p_1pp, "marginal_laser": m})
        print(f"{r['axis']:>17} {r['point']:>26} {100 * eta0:9.3f} "
              f"{100 * eta:8.3f} {100 * (eta - eta0):+7.4f} "
              f"{100 * eta_i0:9.3f} {100 * eta_i:8.3f} {100 * m:6.2f} "
              f"{p_1pp / 1e6:10.1f}  {note}")
    lasing = [x for x in rows if x["eta_dep_baseline"] * 0 == 0
              and x["E_beam_dep_J"] > 0]
    good = [x for x in lasing if x["eta_dep_baseline"] > 1e-3]
    if good:
        lo0 = min(x["eta_dep_baseline"] for x in good)
        hi0 = max(x["eta_dep_baseline"] for x in good)
        lo = min(x["eta_dep_copump"] for x in good)
        hi = max(x["eta_dep_copump"] for x in good)
        dmax = max(x["delta_pp"] for x in good)
        print(f"\nlasing points: intrinsic eta (per deposited J) baseline "
              f"{100 * lo0:.2f}–{100 * hi0:.2f}% -> co-pumped "
              f"{100 * lo:.2f}–{100 * hi:.2f}%; largest lift "
              f"{dmax:+.3f} percentage points")
    out_dir = Path(args.out) if args.out else path.parent
    out = out_dir / "copump_efficiency.json"
    out.write_text(json.dumps({"source": str(path), "rows": rows}, indent=2))
    print(f"report -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
