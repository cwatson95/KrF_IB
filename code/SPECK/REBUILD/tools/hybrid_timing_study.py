#!/usr/bin/env python3
"""HYPER-LASER hybrid-pumping TIMING study — when should the microwave fire?

THE QUESTIONS (owner ask 2026-08-14, following docs/DE-AR_targets.txt):
reheat between lasing pulses? prime the gas before the beam to lower the
e-beam requirements? DE-AR context: e-beam primary + 95 GHz gyrotron
booster (<= 250 kW coupled), and the project's recorded prior mapping
"afterglow reheat helps efficiency; during the discharge it hurts" — whose
observable is 248.5 nm FLUORESCENCE (prop. to small-signal gain), while
REBUILD's 0.1.20 study scored extracted LASER energy. This driver scores
BOTH metrics for every timing scheme, each against its no-microwave twin
with identical beam timing (the microwave is the only difference), so the
marginal numbers subtract the microwave input exactly as the DE-AR
efficiency-metric note demands.

SCHEMES (all on the krf_hyper_mw preset unless overridden):
  afterglow   window after beam-off (the 0.1.20 scheme, DE-AR frequency)
  copump      window ON during the beam (DE-AR expects this to hurt)
  preprime    window BEFORE a delayed beam (beam_t0_ns) — can priming
              lower the beam requirement?
  keeper      two identical beam pulses; the window spans the inter-pulse
              gap (rep-rate keeper). Run at the field-buildup ceiling
              (30 Td ~ the self-sustainment crossing) AND at the
              traveling-wave ceiling (~6 Td) — steady hold needs the
              field to reach ~28 Td (ionization = attachment).

Fluorescence metric: E_fluor = integral n_B/tau_B x h nu x V dt from the
saved KrF*(B) history — the 248.5 nm B->X emission the experiment images;
reported total and post-beam. A gap: the sim's inter-pulse gap is
microseconds (compressed; a real rep-rate gap is ms with gas FLOW — no
flow in the model, so F2 burn carries into pulse 2 in BOTH twins and
differences isolate the microwave).

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/hybrid_timing_study.py
    python3 SPECK/REBUILD/tools/hybrid_timing_study.py --schemes afterglow,copump
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

GYROTRON_W = 2.5e5     # DE-AR: 95 GHz booster, up to ~250 kW coupled

# scheme -> (mw overrides, beam overrides shared by BOTH twins)
SCHEMES = {
    "afterglow": ({"mw_t0_ns": 80.0, "mw_flat_ns": 200.0,
                   "mw_rise_ns": 10.0}, {"t_end_ns": 400.0}),
    "copump": ({"mw_t0_ns": 0.0, "mw_flat_ns": 70.0, "mw_rise_ns": 2.0},
               {"t_end_ns": 400.0}),
    "preprime": ({"mw_t0_ns": 0.0, "mw_flat_ns": 185.0, "mw_rise_ns": 5.0},
                 {"beam_t0_ns": 200.0, "t_end_ns": 600.0}),
    # gap COMPRESSED to 500 ns (stated): the UNDRIVEN twin decays nine
    # decades by ~400 ns (n_e ~1e8 m^-3 — vacuum for every practical
    # purpose) and past ~550 ns the BDF numerical Jacobian overflows on
    # the frozen rows even with restarts; a real rep-rate gap is ms with
    # gas flow, so the keeper numbers are per-gap RATES to scale, not
    # totals. Held-vs-unheld contrast at pulse 2 is already 8 decades.
    "keeper": ({"mw_t0_ns": 80.0, "mw_flat_ns": 390.0, "mw_rise_ns": 10.0},
               {"beam2_t0_ns": 500.0, "t_end_ns": 800.0}),
    "keeper_tw": ({"mw_t0_ns": 80.0, "mw_flat_ns": 390.0,
                   "mw_rise_ns": 10.0, "mw_en_max_td": 6.0},
                  {"beam2_t0_ns": 500.0, "t_end_ns": 800.0}),
}


def build_device(preset, overrides):
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    from REBUILD.run import PRESETS
    spec = PRESETS[preset]
    P = engine_parameters()
    P.set_many(spec["params"], source="problem")
    P.set_many(overrides, source="cli")
    return Device(P, mixture=spec.get("mixture"))


def fluor_J(dev, res, t_from_ns=None):
    """248.5 nm B->X fluorescence energy [J]: n_B/tau_B x h nu x V."""
    upper = dev.optics["upper"]
    tau = float(dev.optics["tau_upper_s"])
    v = float(dev.cell_vols.sum())
    hv = dev.cav.photon_J
    t_s = res["t_ns"] * 1e-9
    p_fl = res["dens_avg"][upper] * v / tau * hv          # W
    if t_from_ns is not None:
        m = res["t_ns"] >= t_from_ns
        t_s, p_fl = t_s[m], p_fl[m]
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(p_fl, t_s))


def laser_after_J(res, t_from_ns):
    m = res["t_ns"] >= t_from_ns
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(res["P_out_W"][m], res["t_ns"][m] * 1e-9))


def run_one_segmented(preset, ov, bounds_ns, n_save, tag, out_dir):
    """One run integrated in checkpoint-restart segments (the checkpoint
    module's stated tolerance-exact semantics). Needed because a SINGLE
    BDF solve across the nine-decade vacuum decay of the inter-pulse gap
    overflows its numerical-Jacobian factor state (measured: chunked
    restarts integrate the same trajectory cleanly). Returns (device,
    merged history dict) shaped like Device.run's result for the metric
    helpers."""
    from REBUILD.MODULES import checkpoint as ckpt
    ov0 = dict(ov)
    ov0["t_end_ns"] = bounds_ns[0]
    dev = build_device(preset, ov0)
    n_of = lambda dt: max(int(n_save * dt / bounds_ns[-1]), 30)
    res = dev.run(n_save=n_of(bounds_ns[0]))
    if not res["success"]:
        sys.exit(f"segmented run FAILED in first segment: {ov0}")
    keys = ("t_ns", "ne_m3", "Te_eV", "P_out_W")
    merged = {k: res[k] for k in keys}
    dens = {sp: res["dens_avg"][sp] for sp in res["dens_avg"]}
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"ckpt_{tag}.npz"
    t_last = bounds_ns[0]
    for b in bounds_ns[1:]:
        ckpt.save(path, res["device"], res["t_ns"][-1] * 1e-9,
                  res["y_final"])
        res = ckpt.resume(path, t_end_ns=b, n_save=n_of(b - t_last))
        if not res["success"]:
            sys.exit(f"segmented run FAILED resuming to {b} ns: {ov}")
        for k in keys:
            merged[k] = np.concatenate([merged[k], res[k][1:]])
        for sp in dens:
            dens[sp] = np.concatenate([dens[sp],
                                       res["dens_avg"][sp][1:]])
        t_last = b
    merged["dens_avg"] = dens
    merged["final_channels_J"] = res["final_channels_J"]
    merged["regime"] = res["regime"]
    return res["device"], merged


def run_pair(preset, name, mw_ov, beam_ov, power_W, n_save, rtol=None,
             out_dir=None):
    """(no-MW twin, MW run) with identical beam timing. Keeper schemes
    integrate in checkpoint-restart segments (see run_one_segmented)."""
    off = dict(beam_ov)
    if rtol is not None:
        off["rtol"] = rtol
    off["mw_power_W"] = 0.0
    on = dict(beam_ov)
    if rtol is not None:
        on["rtol"] = rtol
    on.update(mw_ov)
    on["mw_power_W"] = power_W
    out = []
    for tag, ov in (("off", off), ("on", on)):
        if name.startswith("keeper"):
            t_end = float(ov["t_end_ns"])
            t2 = float(ov["beam2_t0_ns"])
            # 100-ns restarts through the gap (the probe-proven chunk),
            # a boundary just before pulse 2, then the pulse-2 segment
            bounds = [float(b) for b in np.arange(100.0, t2 - 10.0, 100.0)]
            bounds += [t2 - 10.0, t_end]
            dev, res = run_one_segmented(preset, ov, bounds, n_save,
                                         f"{name}_{tag}", out_dir)
        else:
            dev = build_device(preset, ov)
            res = dev.run(n_save=n_save)
            if not res["success"]:
                sys.exit(f"integration FAILED in scheme {name}: {ov}")
        print(res["regime"])
        out.append((dev, res))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_hyper_mw")
    ap.add_argument("--power", type=float, default=GYROTRON_W,
                    help="injected microwave power [W]")
    ap.add_argument("--schemes", default="afterglow,copump,preprime,"
                    "keeper,keeper_tw")
    ap.add_argument("--n-save", type=int, default=241)
    ap.add_argument("--rtol", type=float, default=None,
                    help="BDF rtol override (robustness checks)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    names = [s.strip() for s in args.schemes.split(",") if s.strip()]
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "hybrid_timing")
    rows = []
    print(f"[hybrid timing study] preset={args.preset} "
          f"P_mw={args.power:g} W (DE-AR gyrotron class)\n")
    for name in names:
        mw_ov, beam_ov = SCHEMES[name]
        (dev0, res0), (dev1, res1) = run_pair(
            args.preset, name, mw_ov, beam_ov, args.power, args.n_save,
            rtol=args.rtol, out_dir=out_dir)
        e_mw = (res1["final_channels_J"]["drive_in"]
                - res0["final_channels_J"]["drive_in"])
        d_las = (res1["final_channels_J"]["laser_out"]
                 - res0["final_channels_J"]["laser_out"])
        beam_end = (float(dev0.p["beam_t0_ns"]) + float(dev0.p["beam_flat_ns"])
                    + 2.0 * float(dev0.p["beam_rise_ns"]))
        d_fl = fluor_J(dev1, res1) - fluor_J(dev0, res0)
        d_fl_post = (fluor_J(dev1, res1, beam_end)
                     - fluor_J(dev0, res0, beam_end))
        row = {
            "scheme": name, "P_mw_W": args.power, "E_MW_J": e_mw,
            "E_laser_off_J": res0["final_channels_J"]["laser_out"],
            "E_laser_on_J": res1["final_channels_J"]["laser_out"],
            "dE_laser_J": d_las,
            "marginal_laser": d_las / e_mw if e_mw > 0 else 0.0,
            "E_fluor_off_J": fluor_J(dev0, res0),
            "dE_fluor_J": d_fl, "dE_fluor_postbeam_J": d_fl_post,
            "marginal_fluor": d_fl / e_mw if e_mw > 0 else 0.0,
            "drive_in_off_J": res0["final_channels_J"]["drive_in"],
        }
        if name.startswith("keeper"):
            t2 = float(dev0.p["beam2_t0_ns"])
            row["ne_at_pulse2_off_m3"] = float(
                res0["ne_m3"][int(np.searchsorted(res0["t_ns"], t2)) - 1])
            row["ne_at_pulse2_on_m3"] = float(
                res1["ne_m3"][int(np.searchsorted(res1["t_ns"], t2)) - 1])
            row["pulse2_laser_off_J"] = laser_after_J(res0, t2)
            row["pulse2_laser_on_J"] = laser_after_J(res1, t2)
        if name == "preprime":
            t_b = float(dev0.p["beam_t0_ns"])
            j = max(int(np.searchsorted(res1["t_ns"], t_b)) - 1, 0)
            row["ne_at_beam_on_m3"] = float(res1["ne_m3"][j])
            row["Te_at_beam_on_eV"] = float(res1["Te_eV"][j])
        rows.append(row)
        print(f"  -> {name}: E_MW {e_mw:.3g} J | dE_laser {d_las:+.4f} J "
              f"(marg {100 * row['marginal_laser']:+.2f}%) | dE_fluor "
              f"{d_fl:+.4f} J (post-beam {d_fl_post:+.4f} J, marg "
              f"{100 * row['marginal_fluor']:+.2f}%)\n")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "hybrid_timing_report.json"
    path.write_text(json.dumps(
        {"preset": args.preset, "P_mw_W": args.power, "schemes": rows},
        indent=2))
    print(f"report -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
