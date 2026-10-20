#!/usr/bin/env python3
"""Microwave reheat across the benchmarks — does it raise efficiency anywhere?

Owner ask (2026-08-21): "test the reheating on the different benchmarks to
see if it increases efficiency for all of the cases."

Every BEAM-drive benchmark configuration the package carries (the discharge
anchors cannot take the window — it is a beam-drive overlay) is run as twin
pairs for the two microwave windows, each against its no-microwave twin:

  afterglow   window from beam-off + 10 ns, 150 ns flat (the reheat)
  copump      window on the lasing plateau (second half of the beam pulse,
              ending AT beam-off — the 0.1.22 rule)

at ONE stated coupling rule so the cases compare: the DE-AR booster's
250 kW per litre of pumped cell (2.5e8 W/m^3 x cell volume), 95 GHz, field
ceiling 30 Td — except the ARPA chamber, which keeps its own X-band
hardware (9.5 GHz WR90, 1 MW). Efficiency with the microwave counts the
absorbed microwave energy in the denominator (the DE-AR metric note).

Cases:
  krf_beam_0d     DE-AR-class 0-D cell: 3 atm Ar/Kr/F2, 600 kV, 0.75 m,
                  Electra-class optics (R_OC 0.1, L_op 0.088)
  krf_ne_0d       the Ne-buffered fill, same cell
  kannari_R20     Kannari 1985 device anchor (2.5 atm Ar/Kr/F2 89.8/10/0.2,
                  660 kV, 4x4x44 cm, lossy CaF2 L_op 0.64, R_OC 0.2;
                  j iterated so 66 J deposit — matched by construction)
  electra_arf_1p09  NRL Electra ArF anchor at the 1.09 atm peak point
                  (arf_beam, 193 nm, 2 m, R_OC 0.08, L_op 0.37; j iterated
                  to 0.8 MW/cc flat-top)
  krf_hyper_mw    the DE-AR HYPER cell (2 atm, 10 cm cube)
  krf_arpa_mw     the ARPA chamber column (X-band hardware)
  ar_optimized    the 0.1.24 Ar-optimized point (1 atm, 200 A cm^-2,
                  100 ns, 0.2 m, L_op 0.005, F2 0.2 %)
  ne_optimized    the 0.1.24 Ne-optimized point (2 atm Ne, same knobs)

Reported per case and scheme: E_MW absorbed, dE_laser, marginal, eta
without / with microwave (and the change in percentage points), dFluor
post-beam (the DE-AR observable). Summary: in how many cases each window
RAISES the intrinsic efficiency. 0.1.23 caveat carried: the co-pump sign
is model-fragile (KSTEP bracket, rate tier) — `--tier kinetic` re-runs the
sweep on the PIC_MCC_3D-measured rates for the sign-sensitivity column.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/mw_benchmarks.py
    python3 SPECK/REBUILD/tools/mw_benchmarks.py --cases krf_beam_0d,kannari_R20 --tier kinetic
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
from REBUILD.MODULES.laser import Cavity  # noqa: E402
from REBUILD.MODULES.network import mode_optics  # noqa: E402
from REBUILD.MODULES.params import engine_parameters  # noqa: E402
from REBUILD.run import PRESETS  # noqa: E402

POWER_PER_M3 = 2.5e8          # DE-AR: 250 kW coupled per litre of cell
AFTERGLOW = {"delay_ns": 10.0, "flat_ns": 150.0, "rise_ns": 10.0}


def _preset_params(name, overrides, tier):
    P = engine_parameters()
    P.set_many(PRESETS[name]["params"], source="problem")
    P.set_many(overrides, source="cli")
    P.set("electron_tier", tier, source="cli")
    return P


def _case_configs():
    """name -> (builder(mw_ov, tier) -> Device, notes). Builders take the
    microwave overrides and return a Device; the beam timing is read back
    from the built device's params."""
    cases = {}

    def simple(preset, base_ov=None, mixture=None, cavity_fn=None):
        def build(mw_ov, tier):
            ov = dict(base_ov or {})
            ov.update(mw_ov)
            P = _preset_params(preset, ov, tier)
            cav = cavity_fn(P) if cavity_fn else None
            return Device(P, mixture=mixture or PRESETS[preset].get("mixture"),
                          cavity=cav)
        return build

    cases["krf_beam_0d"] = (simple("krf_beam_0d"), "DE-AR class, 3 atm, 0.75 m, R_OC 0.1/L_op 0.088")
    cases["krf_ne_0d"] = (simple("krf_ne_0d"), "Ne-buffered fill, same cell")
    cases["krf_hyper_mw"] = (simple("krf_hyper_mw", {"mw_freq_GHz": 95.0}),
                             "DE-AR HYPER cell, 2 atm")
    cases["krf_arpa_mw"] = (simple("krf_arpa_mw"), "ARPA chamber column (X-band 9.5 GHz kept)")
    ar_opt = {"pressure_atm": 1.0, "beam_j_A_cm2": 200.0, "beam_flat_ns": 100.0,
              "Lz_m": 0.2, "l_op": 0.005, "t_end_ns": 400.0, "mw_freq_GHz": 95.0}
    cases["ar_optimized"] = (simple("krf_hyper_mw", ar_opt,
                                    {"Ar": 0.948, "Kr": 0.05, "F2": 0.002}),
                             "0.1.24 Ar-optimized (1 atm, F2 0.2%, L_op 0.005)")
    ne_opt = dict(ar_opt)
    ne_opt["pressure_atm"] = 2.0
    cases["ne_optimized"] = (simple("krf_hyper_mw", ne_opt,
                                    {"Ne": 0.948, "Kr": 0.05, "F2": 0.002}),
                             "0.1.24 Ne-optimized (2 atm Ne, F2 0.2%)")

    # --- anchors: rebuilt the way scorecard.py builds them --------------
    def kannari(mw_ov, tier, j_holder=[301.0]):
        opt = mode_optics("krf")
        mix = {"Ar": 0.898, "Kr": 0.10, "F2": 0.002}
        P = engine_parameters()
        P.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 2.5,
                    "beam_kV": 660.0, "beam_j_A_cm2": j_holder[0],
                    "beam_flat_ns": 50.0, "beam_rise_ns": 5.0,
                    "t_end_ns": 120.0, "nx": 1, "ny": 1, "nz": 1,
                    "Lx_m": 0.04, "Ly_m": 0.04, "Lz_m": 0.44,
                    "beam_axis": 0}, source="problem")
        P.set_many(mw_ov, source="cli")
        P.set("electron_tier", tier, source="cli")
        cav = Cavity(l_c_m=0.44, l_g_m=0.44, K_axial=1e-5, R_OC=0.2,
                     k_c_loss=0.0, L_op_fixed=0.64, area_m2=0.04 * 0.04,
                     wavelength_nm=248.0, sigma_se_m2=opt["sigma_se_m2"],
                     absorbers=dict(opt["absorbers"]))
        return Device(P, mixture=mix, cavity=cav)
    cases["kannari_R20"] = (kannari, "Kannari 1985 anchor cell, R_OC 0.2 (66 J matched)")

    def electra(mw_ov, tier, j_holder=[80.0]):
        opt = mode_optics("arf_beam")
        mix = {"Ar": 0.997, "F2": 0.003}
        lx, ly, lz = 0.105, 0.11, 2.0
        P = engine_parameters()
        P.set_many({"mode": "arf_beam", "drive": "beam",
                    "pressure_atm": 1.09, "beam_kV": 500.0,
                    "beam_j_A_cm2": j_holder[0], "beam_flat_ns": 100.0,
                    "beam_rise_ns": 20.0, "t_end_ns": 200.0,
                    "nx": 1, "ny": 1, "nz": 1, "Lx_m": lx, "Ly_m": ly,
                    "Lz_m": lz, "beam_axis": 0}, source="problem")
        P.set_many(mw_ov, source="cli")
        P.set("electron_tier", tier, source="cli")
        cav = Cavity(l_c_m=lz, l_g_m=lz, K_axial=1e-5, R_OC=0.08,
                     k_c_loss=0.0, L_op_fixed=0.37, area_m2=lx * ly,
                     wavelength_nm=193.0, sigma_se_m2=opt["sigma_se_m2"],
                     absorbers=dict(opt["absorbers"]))
        return Device(P, mixture=mix, cavity=cav)
    cases["electra_arf_1p09"] = (electra, "NRL Electra ArF anchor, 1.09 atm peak point (0.8 MW/cc matched)")
    return cases


ANCHOR_TARGETS = {"kannari_R20": ("deposit_J", 66.0),
                  "electra_arf_1p09": ("deposit_J", 0.8e12 * 0.105 * 0.11 * 2.0 * 120e-9)}


def match_anchor_j(name, build, tier):
    """Reproduce the scorecard's input matching: scale j until the
    no-microwave run deposits the anchor's recorded energy."""
    if name not in ANCHOR_TARGETS:
        return None
    _kind, target = ANCHOR_TARGETS[name]
    j = None
    for _attempt in range(4):
        dev = build({"mw_power_W": 0.0} if j is None else
                    {"mw_power_W": 0.0, "beam_j_A_cm2": j}, tier)
        res = dev.run(n_save=25)
        dep = res["final_channels_J"]["drive_in"]
        j = float(dev.p["beam_j_A_cm2"])
        if abs(dep - target) < 0.005 * target:
            break
        j *= target / dep
    return j


def run_segmented(dev, bounds_ns, n_save):
    """Integrate one device in restart segments (the 0.1.21 cure for the
    BDF numerical-Jacobian overflow across a deep vacuum afterglow — the
    Ne fill empties fastest): Device.run to the first bound, then
    solve_ivp continuations on the SAME device (no rebuild — the anchors
    carry custom cavities the checkpoint restore would not reproduce),
    histories merged into one run()-shaped result."""
    import time as _time
    from scipy.integrate import solve_ivp
    n_of = lambda dt: max(int(n_save * dt / bounds_ns[-1]), 20)
    dev.p.set("t_end_ns", float(bounds_ns[0]), source="api")
    res = dev.run(n_save=n_of(bounds_ns[0]))
    if not res["success"]:
        return res
    keys = ("t_ns", "ne_m3", "Te_eV", "P_out_W", "g0_m1", "Tg_K")
    merged = {k: res[k] for k in keys}
    dens = {sp: res["dens_avg"][sp] for sp in res["dens_avg"]}
    t_last = float(bounds_ns[0])
    y = res["y_final"]
    atol = np.full(dev.lay.size, 1e3)
    ca = dev.lay.cells_view(atol)
    ca[dev.lay.i_ue] = 1e-10
    ca[dev.lay.i_tg] = 1e-6
    atol[dev.lay.i_ph] = 1e-3
    atol[dev.lay.i_ledger0:] = 1e-12
    for b in bounds_ns[1:]:
        t_eval = np.linspace(t_last * 1e-9, b * 1e-9, n_of(b - t_last))
        w0 = _time.time()
        sol = solve_ivp(dev.rhs, (t_last * 1e-9, b * 1e-9), y, method="BDF",
                        rtol=dev.p["rtol"], atol=atol, t_eval=t_eval,
                        max_step=2e-9)
        seg = dev._collect(sol, y, _time.time() - w0, False)
        if not seg["success"]:
            return seg
        for k in keys:
            merged[k] = np.concatenate([merged[k], seg[k][1:]])
        for sp in dens:
            dens[sp] = np.concatenate([dens[sp], seg["dens_avg"][sp][1:]])
        y = seg["y_final"]
        t_last = float(b)
    merged["dens_avg"] = dens
    merged["final_channels_J"] = seg["final_channels_J"] if len(bounds_ns) > 1 else res["final_channels_J"]
    merged["regime"] = res["regime"]
    merged["success"] = True
    return merged


def segment_bounds(beam_end_ns, t_end_ns):
    """First segment to beam-off + 50 ns (the pulse + early afterglow),
    then 100-ns restarts to the end."""
    b = [beam_end_ns + 50.0]
    while b[-1] + 100.0 < t_end_ns - 1e-9:
        b.append(b[-1] + 100.0)
    if b[-1] < t_end_ns - 1e-9:
        b.append(t_end_ns)
    return b


def fluor_post(dev, res, t_from_ns):
    upper = dev.optics["upper"]
    tau = float(dev.optics["tau_upper_s"])
    v = float(dev.cell_vols.sum())
    m = res["t_ns"] >= t_from_ns
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(res["dens_avg"][upper][m] * v / tau * dev.cav.photon_J,
                      res["t_ns"][m] * 1e-9))


def run_case(name, build, tier, n_save):
    j = match_anchor_j(name, build, tier)
    base_ov = {"mw_power_W": 0.0}
    if j is not None:
        base_ov["beam_j_A_cm2"] = j
    probe = build(dict(base_ov), tier)
    p = probe.p
    beam_t0 = float(p["beam_t0_ns"])
    pulse = float(p["beam_flat_ns"]) + 2.0 * float(p["beam_rise_ns"])
    beam_end = beam_t0 + pulse
    vol = float(probe.cell_vols.sum())
    power = (float(PRESETS["krf_arpa_mw"]["params"]["mw_power_W"])
             if name == "krf_arpa_mw" else POWER_PER_M3 * vol)
    t_end = max(float(p["t_end_ns"]), beam_end + AFTERGLOW["delay_ns"]
                + AFTERGLOW["flat_ns"] + 2 * AFTERGLOW["rise_ns"] + 80.0)
    schemes = {
        "afterglow": {"mw_t0_ns": beam_end + AFTERGLOW["delay_ns"],
                      "mw_flat_ns": AFTERGLOW["flat_ns"],
                      "mw_rise_ns": AFTERGLOW["rise_ns"]},
        "copump": {"mw_t0_ns": beam_t0 + 0.5 * pulse,
                   "mw_flat_ns": max(0.5 * pulse - 4.0, 1.0),
                   "mw_rise_ns": 2.0},
    }
    common = dict(base_ov)
    common["t_end_ns"] = t_end
    if name != "krf_arpa_mw":
        common["mw_en_max_td"] = 30.0
    bounds = segment_bounds(beam_end, t_end)
    d0 = build(dict(common), tier)
    r0 = run_segmented(d0, bounds, n_save)
    if not r0["success"]:
        sys.exit(f"{name}: baseline failed")
    rows = {}
    for sch, win in schemes.items():
        ov = dict(common)
        ov.update(win)
        ov["mw_power_W"] = power
        d1 = build(ov, tier)
        r1 = run_segmented(d1, bounds, n_save)
        if not r1["success"]:
            sys.exit(f"{name}/{sch}: run failed")
        e_mw = r1["final_channels_J"]["drive_in"] - r0["final_channels_J"]["drive_in"]
        d_las = r1["final_channels_J"]["laser_out"] - r0["final_channels_J"]["laser_out"]
        e_dep = r0["final_channels_J"]["drive_in"]
        eta0 = r0["final_channels_J"]["laser_out"] / e_dep
        eta1 = r1["final_channels_J"]["laser_out"] / (e_dep + e_mw)
        d_fl = fluor_post(d1, r1, beam_end) - fluor_post(d0, r0, beam_end)
        rows[sch] = {"E_MW_J": e_mw, "dE_laser_J": d_las,
                     "marginal": d_las / e_mw if e_mw > 0 else 0.0,
                     "eta0": eta0, "eta_mw": eta1,
                     "d_eta_pp": 100 * (eta1 - eta0),
                     "dFluor_postbeam_J": d_fl, "window": win}
    return {"case": name, "tier": tier, "P_mw_W": power, "cell_m3": vol,
            "beam_end_ns": beam_end, "E_dep_J": r0["final_channels_J"]["drive_in"],
            "E_laser0_J": r0["final_channels_J"]["laser_out"],
            "eta0": rows["afterglow"]["eta0"], "j_A_cm2": float(p["beam_j_A_cm2"]),
            "schemes": rows, "regime": r0["regime"]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    cfg = _case_configs()
    ap.add_argument("--cases", default=",".join(cfg))
    ap.add_argument("--tier", default="loki", choices=("loki", "kinetic"))
    ap.add_argument("--n-save", type=int, default=161)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    names = [c.strip() for c in args.cases.split(",") if c.strip()]
    out = []
    hdr = (f"{'case':>17} {'P_mw W':>8} {'E_dep J':>8} {'eta0%':>6} | "
           f"{'AFTERGLOW: E_MW mJ':>18} {'dE mJ':>8} {'marg%':>6} {'d_eta pp':>8} "
           f"{'dFl mJ':>7} | {'COPUMP: E_MW mJ':>15} {'dE mJ':>8} {'marg%':>6} "
           f"{'d_eta pp':>8}")
    print(f"[mw_benchmarks] tier={args.tier}, rule {POWER_PER_M3:g} W/m^3 "
          f"(250 kW/L), 95 GHz (ARPA keeps X-band)\n")
    print(hdr)
    print("-" * len(hdr))
    for name in names:
        build, note = cfg[name]
        r = run_case(name, build, args.tier, args.n_save)
        r["note"] = note
        out.append(r)
        a, c = r["schemes"]["afterglow"], r["schemes"]["copump"]
        print(f"{name:>17} {r['P_mw_W']:8.3g} {r['E_dep_J']:8.2f} "
              f"{100 * r['eta0']:6.2f} | {1e3 * a['E_MW_J']:18.2f} "
              f"{1e3 * a['dE_laser_J']:8.3f} {100 * a['marginal']:6.2f} "
              f"{a['d_eta_pp']:+8.4f} {1e3 * a['dFluor_postbeam_J']:7.3f} | "
              f"{1e3 * c['E_MW_J']:15.2f} {1e3 * c['dE_laser_J']:8.3f} "
              f"{100 * c['marginal']:6.2f} {c['d_eta_pp']:+8.4f}", flush=True)
    n = len(out)
    up_a = sum(r["schemes"]["afterglow"]["d_eta_pp"] > 0 for r in out)
    up_c = sum(r["schemes"]["copump"]["d_eta_pp"] > 0 for r in out)
    print(f"\nSUMMARY ({args.tier} tier): afterglow reheat RAISES intrinsic "
          f"efficiency in {up_a}/{n} cases; plateau co-pump in {up_c}/{n} "
          f"cases; largest |d_eta| "
          f"{max(abs(r['schemes'][s]['d_eta_pp']) for r in out for s in r['schemes']):.4f} pp")
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "mw_benchmarks")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"mw_benchmarks_{args.tier}.json"
    path.write_text(json.dumps({"tier": args.tier, "rule_W_m3": POWER_PER_M3,
                                "cases": out}, indent=2, default=float))
    print(f"report -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
