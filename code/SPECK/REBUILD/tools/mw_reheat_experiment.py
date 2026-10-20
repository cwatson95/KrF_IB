#!/usr/bin/env python3
"""The ARPA-chamber microwave-reheat design study: baseline vs reheat.

THE QUESTION (owner ask, 2026-08-14): pump the chamber's KrF fill with the
e-beam, let it lase; then, in the afterglow as the electron gas cools,
reheat it through the three WR90 feeds — does the microwave energy buy
laser output at a useful marginal efficiency?

WHAT THIS SCRIPT DOES (a driver, not an engine — feeds nothing back):
runs the SAME preset at mw_power_W = 0 (baseline) and at each scanned
injected power, and reports per run

  E_MW      = drive_in - drive_in(baseline)   [J]  the ABSORBED
              microwave energy (exact: the beam waveform and gas density
              are identical by construction; cross-checked against a
              trapezoid quadrature of the settled-field swarm balance
              P_abs(t) = e <n_e> v_d(E/N) E fac V from saved history)
  dE_laser  = laser_out - laser_out(baseline) [J]
  marginal  = dE_laser / E_MW                 [the answer]
  absorbed fraction, and the settled E/N at window start/end (the field
  solver sags the field into dense plasma and rides the mw_en_max_td
  ceiling in dilute plasma — absorbed never exceeds injected). The
  STATED feed budget for context: three WR90 feeds, pulsed X-band
  magnetron class ~0.3 MW each -> ~1 MW total injected.

Construction (matched by construction, never scored): chamber geometry from
freecad/arpa_chamber_ports.json; beam 500 kV / 140 A cm^-2 / 100 ns; fill
2 atm Ar/Kr/F2 0.949/0.05/0.001. Outputs (scored against nothing — this is
a PREDICTION): every energy, density and temperature history.

Stated model limits (the engine blocks carry them too): uniform MW field
over the pumped column (overmoded cavity; at pulse-peak n_e the collisional
absorption length is cm-class, so early-window heating is really skin-
localized — the uniform idealization is best late in the afterglow);
fluid-Te rate keying; no wall cooling of the gas (ms-scale, off-window).

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/mw_reheat_experiment.py
    python3 SPECK/REBUILD/tools/mw_reheat_experiment.py \
        --power-scan 3e5,1e6,1e7 --t0 130 --preset krf_arpa_mw
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

FEED_BUDGET_W = 1.0e6      # stated context: 3x WR90 pulsed magnetrons


def build_device(preset: str, overrides: dict):
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    from REBUILD.run import PRESETS
    spec = PRESETS[preset]
    P = engine_parameters()
    P.set_many(spec["params"], source="problem")
    P.set_many(overrides, source="cli")
    return Device(P, mixture=spec.get("mixture"))


def mw_power_history(dev, res):
    """(P_abs(t) [W], E/N(t) [Td]) reconstructed from the saved history —
    exact for a 0-D run (one cell), a volume-average approximation on a
    grid. Mirrors the engine's rhs MW block term for term, including the
    settled-field solver."""
    from REBUILD.MODULES.engine import _mw_envelope
    from REBUILD.MODULES.units import E_CHARGE, M_E
    p = dev.p
    t = res["t_ns"]
    if float(p["mw_power_W"]) <= 0.0:
        return np.zeros_like(t), np.zeros_like(t)
    v_tot = float(dev.cell_vols.sum())
    p_out = np.zeros_like(t)
    en_out = np.zeros_like(t)
    for j, t_ns in enumerate(t):
        env = _mw_envelope(t_ns * 1e-9, p["mw_t0_ns"], p["mw_rise_ns"],
                           p["mw_flat_ns"])
        if env <= 0.0:
            continue
        ne = float(res["ne_m3"][j])
        en = dev._mw_field_td(ne * v_tot, float(p["mw_power_W"]) * env)
        en_out[j] = en
        v_d = dev._mw_swarm.drift_velocity(en)
        if v_d <= 0.0:
            continue
        e_f = en * 1e-21 * dev.n_tot
        nu = E_CHARGE * e_f / (M_E * v_d)
        fac = nu * nu / (nu * nu + dev._mw_omega ** 2)
        p_out[j] = E_CHARGE * ne * v_d * e_f * fac * v_tot
    return p_out, en_out


def window_slice(res, t0_ns, t1_ns):
    t = res["t_ns"]
    return (t >= t0_ns) & (t <= t1_ns)


def run_study(preset, p_values, t0_ns, t_end_ns, n_save, out_dir,
              en_max_td=None):
    rows = []
    devs = {}
    results = {}
    for pw in [0.0] + [x for x in p_values if x > 0.0]:
        ov = {"mw_power_W": pw}
        if en_max_td is not None:
            ov["mw_en_max_td"] = en_max_td
        if t0_ns is not None:
            ov["mw_t0_ns"] = t0_ns
        if t_end_ns is not None:
            ov["t_end_ns"] = t_end_ns
        dev = build_device(preset, ov)
        res = dev.run(n_save=n_save)
        devs[pw], results[pw] = dev, res
        print(res["regime"])
        if not res["success"]:
            sys.exit(f"integration FAILED at mw_power_W={pw}")
    base = results[0.0]
    e_beam = base["final_channels_J"]["drive_in"]
    e_laser0 = base["final_channels_J"]["laser_out"]
    dev0 = devs[0.0]
    t0 = float(dev0.p["mw_t0_ns"])
    t1 = t0 + 2.0 * float(dev0.p["mw_rise_ns"]) + float(dev0.p["mw_flat_ns"])
    i_t0 = int(np.searchsorted(base["t_ns"], t0))
    i_t1 = min(int(np.searchsorted(base["t_ns"], t1)), len(base["t_ns"]) - 1)
    print()
    print(f"BASELINE (beam only): deposited {e_beam:.1f} J, laser "
          f"{e_laser0:.2f} J (eta {100 * base['eta_intrinsic']:.2f}%), "
          f"g0 peak {np.max(base['g0_m1']):.2f} /m, "
          f"n_e at window start {base['ne_m3'][i_t0]:.2e} m^-3, "
          f"at window end {base['ne_m3'][i_t1]:.2e} m^-3")
    print()
    hdr = (f"{'P_inj W':>9} {'E_MW J':>9} {'quad J':>9} {'abs %':>6} "
           f"{'EN0/EN1 Td':>12} {'E_las J':>9} {'dE_las J':>10} "
           f"{'marg %':>8} {'Te_pk eV':>8} {'ne_end':>9}  vs 1MW budget")
    print(hdr)
    print("-" * len(hdr))
    for pw in p_values:
        if pw <= 0.0:
            continue
        dev, res = devs[pw], results[pw]
        e_mw = res["final_channels_J"]["drive_in"] - e_beam
        p_hist, en_hist = mw_power_history(dev, res)
        trap = getattr(np, "trapezoid", None) or np.trapz
        e_quad = float(trap(p_hist, res["t_ns"] * 1e-9))
        e_inj = pw * (float(dev.p["mw_flat_ns"])
                      + float(dev.p["mw_rise_ns"])) * 1e-9
        d_laser = res["final_channels_J"]["laser_out"] - e_laser0
        marg = d_laser / e_mw if e_mw > 0 else 0.0
        w = window_slice(res, t0, t1)
        te_pk = float(np.max(res["Te_eV"][w])) if w.any() else 0.0
        ne_end = float(res["ne_m3"][i_t1])
        win = np.nonzero(en_hist > 0.0)[0]
        en0 = float(en_hist[win[0]]) if win.size else 0.0
        en1 = float(en_hist[win[-1]]) if win.size else 0.0
        feas = ("within" if pw <= FEED_BUDGET_W
                else f"x{pw / FEED_BUDGET_W:.0f} OVER")
        rows.append({
            "mw_power_W": pw, "E_MW_absorbed_J": e_mw,
            "E_MW_quad_J": e_quad, "E_injected_J": e_inj,
            "absorbed_frac": e_mw / e_inj if e_inj > 0 else 0.0,
            "en_window_start_td": en0, "en_window_end_td": en1,
            "E_laser_J": res["final_channels_J"]["laser_out"],
            "dE_laser_J": d_laser, "marginal_eff": marg,
            "eta_intrinsic": res["eta_intrinsic"],
            "Te_peak_window_eV": te_pk, "ne_window_end_m3": ne_end,
            "within_1MW_budget": bool(pw <= FEED_BUDGET_W),
            "validity": {k: v for k, v in res["validity"].items()
                         if str(k).startswith("mw_")},
        })
        print(f"{pw:9.3g} {e_mw:9.3g} {e_quad:9.3g} "
              f"{100 * e_mw / max(e_inj, 1e-300):6.1f} "
              f"{en0:5.2f}/{en1:5.2f} "
              f"{res['final_channels_J']['laser_out']:9.3f} "
              f"{d_laser:10.4f} {100 * marg:8.3f} "
              f"{te_pk:8.2f} {ne_end:9.2e}  {feas}")
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "preset": preset,
        "baseline": {
            "E_beam_dep_J": e_beam, "E_laser_J": e_laser0,
            "eta_intrinsic": base["eta_intrinsic"],
            "g0_peak_m1": float(np.max(base["g0_m1"])),
            "ne_at_window_start_m3": float(base["ne_m3"][i_t0]),
            "ne_at_window_end_m3": float(base["ne_m3"][i_t1]),
            "Tg_end_K": float(base["Tg_K"][-1]),
            "regime": base["regime"],
        },
        "window_ns": [t0, t1],
        "feed_budget_W": FEED_BUDGET_W,
        "scan": rows,
    }
    path = out_dir / "mw_reheat_report.json"
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nreport -> {path}")
    return payload


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_arpa_mw")
    ap.add_argument("--power-scan", default="3e5,1e6,3e6,1e7",
                    help="comma list of mw_power_W values [W injected]")
    ap.add_argument("--t0", type=float, default=None, metavar="NS",
                    help="window start (defaults to the preset's)")
    ap.add_argument("--t-end", type=float, default=None, metavar="NS")
    ap.add_argument("--n-save", type=int, default=241)
    ap.add_argument("--out", default=None)
    ap.add_argument("--en-max", type=float, default=None, metavar="TD",
                    help="cavity-field ceiling for the scan (the "
                         "empty-cavity buildup scales as sqrt(P) x "
                         "sqrt(Q); raise it for what-would-it-take "
                         "brackets)")
    args = ap.parse_args(argv)
    p_values = [float(x) for x in args.power_scan.split(",") if x.strip()]
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "mw_reheat")
    run_study(args.preset, p_values, args.t0, args.t_end, args.n_save,
              out_dir, en_max_td=args.en_max)
    return 0


if __name__ == "__main__":
    sys.exit(main())
