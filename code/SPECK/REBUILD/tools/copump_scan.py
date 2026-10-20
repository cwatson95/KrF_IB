#!/usr/bin/env python3
"""Co-pump operating map — where does microwave DURING the beam pay?

Owner ask (2026-08-20): "what range of operating parameters allows for
co-pumping?" The 0.1.21 timing study found the one laser-positive scheme:
the microwave window ON during the e-beam pulse (+4.1% marginal at the
DE-AR point: 2 atm, 250 kW, 100 A cm^-2, window over the whole beam). This
driver maps that marginal over the DE-AR operating ranges
(docs/DE-AR_targets.txt: 0.5-3 atm; booster up to ~250 kW coupled from a
2 MW source; beam 600 kV) and over the two knobs the beam operator owns:
pump rate (current density) and the window's placement inside the pulse.

Every point is a twin pair (identical beam/gas, microwave the only
difference) on the krf_hyper_mw preset; reported per point:
  marginal_laser  = dE_laser / E_MW_absorbed   (the payoff per microwave J)
  dE_laser, E_MW, the baseline laser energy (a point that does not lase in
  the baseline is flagged — a marginal on ~0 is noise), E/N settled in the
  window (max / mean from the history), peak Te in the window, and the
  fluorescence marginal (the DE-AR observable).

Axes:
  power x pressure     at fixed beam hardware (j = 100 A cm^-2): the
                       deposition falls with density, so low-pressure
                       points may sit below threshold — that IS an
                       operating-range fact, flagged not hidden
  pressure @ const dep j scaled 1/p so the deposited power density is the
                       reference value: isolates the pressure physics
  beam pump rate       j 50 / 100 / 200 at 2 atm, 250 kW
  window placement     full beam / first half / second half / beam +
                       early afterglow / late-beam + afterglow
  optics margin        R_OC 0.9 vs 0.95 (threshold-margin amplification)

Caveats carried: the LOKI electron-impact tables are the v2 fill's
(Ar/Kr/F2 0.949/0.05/0.001) — a Kr- or F2-fraction axis would need new
LoKI-MC tables (LoKI-MC is in-repo) and is NOT scanned here; uniform MW
field over the cube (at the beam-on densities the 95 GHz wave is
collisionally absorbed in a cm-class skin — the uniform idealization
overstates penetration during the pulse, best in the late window);
fluid-Te keying; internal optics stated. Predictions, scored against
nothing.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/copump_scan.py
    python3 SPECK/REBUILD/tools/copump_scan.py --axes power_pressure,timing
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

from REBUILD.tools.hybrid_timing_study import (build_device, fluor_J,  # noqa: E402
                                                 run_pair)
from REBUILD.tools.mw_reheat_experiment import mw_power_history  # noqa: E402

REF_P_ATM = 2.0
REF_J = 100.0
REF_POWER = 2.5e5
FULL_WINDOW = {"mw_t0_ns": 0.0, "mw_flat_ns": 70.0, "mw_rise_ns": 2.0}
LASING_FLOOR_J = 1e-3       # baseline below this: "does not lase" flag


def axes_points():
    pts = []
    for p in (0.5, 1.0, 2.0, 3.0):
        for pw in (5e4, 2.5e5, 1e6, 3e6):
            pts.append(("power_pressure", f"p{p:g}atm_P{pw:g}W", pw,
                        dict(FULL_WINDOW),
                        {"pressure_atm": p, "beam_j_A_cm2": REF_J,
                         "t_end_ns": 250.0}))
    for p in (0.5, 1.0, 3.0):
        pts.append(("pressure_constdep", f"p{p:g}atm_j{REF_J * REF_P_ATM / p:g}",
                    REF_POWER, dict(FULL_WINDOW),
                    {"pressure_atm": p, "beam_j_A_cm2": REF_J * REF_P_ATM / p,
                     "t_end_ns": 250.0}))
    for j in (50.0, 200.0):
        pts.append(("pump_rate", f"j{j:g}", REF_POWER, dict(FULL_WINDOW),
                    {"pressure_atm": REF_P_ATM, "beam_j_A_cm2": j,
                     "t_end_ns": 250.0}))
    for tag, win in (("first_half", {"mw_t0_ns": 0.0, "mw_flat_ns": 33.0,
                                     "mw_rise_ns": 2.0}),
                     ("second_half", {"mw_t0_ns": 35.0, "mw_flat_ns": 35.0,
                                      "mw_rise_ns": 2.0}),
                     ("beam_plus_early_afterglow",
                      {"mw_t0_ns": 0.0, "mw_flat_ns": 110.0,
                       "mw_rise_ns": 2.0}),
                     ("late_beam_plus_afterglow",
                      {"mw_t0_ns": 40.0, "mw_flat_ns": 70.0,
                       "mw_rise_ns": 2.0})):
        pts.append(("timing", tag, REF_POWER, win,
                    {"pressure_atm": REF_P_ATM, "beam_j_A_cm2": REF_J,
                     "t_end_ns": 250.0}))
    pts.append(("optics", "r_oc0.90", REF_POWER, dict(FULL_WINDOW),
                {"pressure_atm": REF_P_ATM, "beam_j_A_cm2": REF_J,
                 "r_oc": 0.90, "t_end_ns": 250.0}))
    return pts


def window_stats(dev, res):
    p_hist, en_hist = mw_power_history(dev, res)
    win = en_hist > 0.0
    if not win.any():
        return 0.0, 0.0, 0.0
    te = res["Te_eV"][win]
    return float(en_hist[win].max()), float(en_hist[win].mean()), float(te.max())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_hyper_mw")
    ap.add_argument("--axes", default="power_pressure,pressure_constdep,"
                    "pump_rate,timing,optics")
    ap.add_argument("--n-save", type=int, default=161)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    want = {a.strip() for a in args.axes.split(",") if a.strip()}
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "copump_scan")
    rows = []
    hdr = (f"{'axis':>17} {'point':>22} {'P_mw W':>8} {'E_MW mJ':>8} "
           f"{'E_las0 J':>9} {'dE mJ':>8} {'marg%':>7} {'flmarg%':>8} "
           f"{'EN max/mean Td':>15} {'Te_pk':>6}  note")
    print(hdr)
    print("-" * len(hdr))
    for axis, tag, pw, win, beam_ov in axes_points():
        if axis not in want:
            continue
        (dev0, res0), (dev1, res1) = run_pair(args.preset, f"copump_{tag}",
                                              win, beam_ov, pw, args.n_save)
        e_mw = (res1["final_channels_J"]["drive_in"]
                - res0["final_channels_J"]["drive_in"])
        e_las0 = res0["final_channels_J"]["laser_out"]
        d_las = res1["final_channels_J"]["laser_out"] - e_las0
        d_fl = fluor_J(dev1, res1) - fluor_J(dev0, res0)
        marg = d_las / e_mw if e_mw > 0 else 0.0
        fl_marg = d_fl / e_mw if e_mw > 0 else 0.0
        en_max, en_mean, te_pk = window_stats(dev1, res1)
        note = "" if e_las0 >= LASING_FLOOR_J else "BASELINE DOES NOT LASE"
        g0_pk = float(np.max(res0["g0_m1"]))
        rows.append({
            "axis": axis, "point": tag, "P_mw_W": pw, "window": win,
            "beam_overrides": beam_ov, "E_MW_J": e_mw,
            "E_laser_baseline_J": e_las0, "dE_laser_J": d_las,
            "marginal_laser": marg, "marginal_fluor": fl_marg,
            "g0_peak_baseline_m1": g0_pk,
            "en_window_max_td": en_max, "en_window_mean_td": en_mean,
            "Te_peak_window_eV": te_pk,
            "baseline_lases": bool(e_las0 >= LASING_FLOOR_J),
            "E_beam_dep_J": res0["final_channels_J"]["drive_in"],
        })
        print(f"{axis:>17} {tag:>22} {pw:8.3g} {1e3 * e_mw:8.2f} "
              f"{e_las0:9.4f} {1e3 * d_las:8.3f} {100 * marg:7.2f} "
              f"{100 * fl_marg:8.2f} {en_max:7.3f}/{en_mean:6.3f} "
              f"{te_pk:6.2f}  {note}", flush=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "copump_scan_report.json"
    path.write_text(json.dumps({"preset": args.preset, "points": rows},
                               indent=2))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
