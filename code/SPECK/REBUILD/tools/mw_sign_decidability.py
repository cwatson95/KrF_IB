"""Can the models TELL whether any microwave reheat is real? ([0.1.27])

Two discriminators, both reusing recorded configurations verbatim:

A) THE CO-PUMP SIGN'S DECIDABILITY. 0.1.23 measured that the plateau
   co-pump marginal at the HYPER point is hostage to exactly ONE rate:
   KSTEP (e + Kr* -> Kr+ + 2e), negative at BOTH x1/3 and x3 around a
   positive nominal — so two zero crossings exist. This part maps the
   marginal vs the KSTEP multiplier and locates the crossings. The
   positive window's width against the EST tier's own x3 uncertainty
   convention IS the answer to "can the model tell": if the window sits
   inside the uncertainty, the sign is undecidable in-model and the map
   states exactly how well KSTEP must be known to decide it.

B) IS THE ONE AFTERGLOW-POSITIVE CASE SINGLE-RATE-HOSTAGE TOO? The
   Ne-optimized cell is the only benchmark whose afterglow-reheat
   marginal is laser-positive (LOKI +10.0% post-errata). This part runs
   the full 0.1.23-style EST x1/3..x3 + LIT-conflict bracket pass on
   THAT twin (segmented integration, the 0.1.21 cure). If the sign
   survives every variant, the "high-margin still-lasing afterglow"
   class is the models' one bracket-robust positive; if a row flips it,
   that row is named.

Run from the Kinetic_Modeling root:
    python3 SPECK/REBUILD/tools/mw_sign_decidability.py
Writes SPECK/rebuild_out/mw_benchmarks/sign_decidability.json.
"""
import json
import sys
import time
from dataclasses import replace as _rep
from pathlib import Path

import numpy as np

for base in Path(__file__).resolve().parents:
    if (base / "REBUILD" / "MODULES").is_dir():
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.MODULES.chemistry import CompiledNetwork          # noqa: E402
from REBUILD.MODULES.scorecard import (CONFLICT_ALTERNATES,    # noqa: E402
                                       _bracket_rate,
                                       est_active_and_bound)
from REBUILD.tools.copump_verify import (BEAM, PLATEAU, POWER,  # noqa: E402
                                         solve)
from REBUILD.tools.hybrid_timing_study import build_device      # noqa: E402
from REBUILD.tools.mw_benchmarks import (AFTERGLOW,             # noqa: E402
                                         POWER_PER_M3,
                                         _case_configs,
                                         run_segmented,
                                         segment_bounds)

OUT = Path(__file__).resolve()
for base in OUT.parents:
    if (base / "REBUILD" / "MODULES").is_dir():
        OUT = base / "rebuild_out" / "mw_benchmarks" / "sign_decidability.json"
        break

KSTEP_FACTORS = [1.0 / 3.0, 0.45, 0.6, 0.8, 1.0, 1.25, 1.6, 2.2, 3.0]


def _swap(devs, rows):
    for d in devs:
        d.compiled = CompiledNetwork(rows, d.lay, d.tiers, d.optics)
        d.compiled.set_photon_energy(d.hv_eV)


def _marginal_pair(r0, r1):
    e_mw = r1["final_channels_J"]["drive_in"] - r0["final_channels_J"]["drive_in"]
    d_las = r1["final_channels_J"]["laser_out"] - r0["final_channels_J"]["laser_out"]
    return e_mw, d_las, (d_las / e_mw if e_mw > 0 else 0.0)


def part_kstep_map():
    print("[A] co-pump sign vs the KSTEP multiplier (HYPER plateau point, "
          "the 0.1.23 configuration verbatim; LOKI tier)")
    off = dict(BEAM)
    off["mw_power_W"] = 0.0
    on = dict(BEAM)
    on.update(PLATEAU)
    on["mw_power_W"] = POWER
    d0 = build_device("krf_hyper_mw", off)
    d1 = build_device("krf_hyper_mw", on)
    kstep = [r for r in d0.net if r.rid == "KSTEP"]
    if len(kstep) != 1:
        sys.exit(f"expected exactly one KSTEP row, found {len(kstep)}")
    print(f"    KSTEP row: {kstep[0].rid}  (the one sign-flipping bracket "
          f"of 0.1.23)")
    rows = []
    print(f"{'x':>8} {'E_off J':>10} {'dE mJ':>9} {'marg%':>8}")
    for f in KSTEP_FACTORS:
        net2 = [_rep(x, rate=_bracket_rate(x.rate, f))
                if x.rid == "KSTEP" else x for x in d0.net]
        _swap((d0, d1), net2)
        r0 = solve(d0)
        r1 = solve(d1)
        e_mw, d_las, marg = _marginal_pair(r0, r1)
        rows.append({"factor": f, "E_laser_off_J":
                     r0["final_channels_J"]["laser_out"],
                     "E_MW_J": e_mw, "dE_laser_J": d_las, "marginal": marg})
        print(f"{f:8.3f} {r0['final_channels_J']['laser_out']:10.6f} "
              f"{1e3 * d_las:9.4f} {100 * marg:8.3f}", flush=True)
    _swap((d0, d1), d0.net)

    fs = np.array([r["factor"] for r in rows])
    ms = np.array([r["marginal"] for r in rows])
    crossings = []
    for i in range(len(fs) - 1):
        if ms[i] == 0.0 or ms[i] * ms[i + 1] >= 0:
            continue
        lf = np.log(fs[i]) + (np.log(fs[i + 1]) - np.log(fs[i])) * (
            -ms[i] / (ms[i + 1] - ms[i]))
        crossings.append(float(np.exp(lf)))
    verdict = {}
    if len(crossings) == 2:
        f_lo, f_hi = crossings
        width = f_hi / f_lo
        need = float(np.sqrt(width))
        verdict = {"f_lo": f_lo, "f_hi": f_hi, "window_ratio": width,
                   "decisive_knowledge_factor": need}
        print(f"  positive window: KSTEP x{f_lo:.2f} .. x{f_hi:.2f} "
              f"(ratio {width:.1f}); nominal marginal "
              f"{100 * ms[fs == 1.0][0]:+.2f}%")
        print(f"  EST uncertainty convention x3 straddles BOTH crossings -> "
              f"the sign is UNDECIDABLE in-model at current rate knowledge; "
              f"pinning KSTEP to better than ~x{need:.2f} decides it.")
    else:
        print(f"  crossings found: {crossings} (expected 2) — map the "
              f"curve before interpreting")
    return {"config": {"preset": "krf_hyper_mw", "window": PLATEAU,
                       "beam": BEAM, "power_W": POWER, "tier": "loki"},
            "rows": rows, "crossings": crossings, **verdict}


def part_neopt_afterglow_brackets():
    print("\n[B] EST x1/3..x3 + LIT-conflict brackets on the ONE afterglow-"
          "positive case (ne_optimized, segmented; LOKI tier)")
    build, note = _case_configs()["ne_optimized"]
    d_probe = build({"mw_power_W": 0.0}, "loki")
    p = d_probe.p
    beam_end = float(p["beam_rise_ns"] + p["beam_flat_ns"])
    t_end = float(p["t_end_ns"])
    cell_m3 = float(p["Lx_m"] * p["Ly_m"] * p["Lz_m"])
    p_mw = POWER_PER_M3 * cell_m3
    mw_ov = {"mw_power_W": p_mw,
             "mw_t0_ns": beam_end + AFTERGLOW["delay_ns"],
             "mw_flat_ns": AFTERGLOW["flat_ns"],
             "mw_rise_ns": AFTERGLOW["rise_ns"]}
    bounds = segment_bounds(beam_end, t_end)
    print(f"    {note}; window {mw_ov['mw_t0_ns']:.0f}"
          f"+{mw_ov['mw_flat_ns']:.0f} ns at {p_mw / 1e3:.0f} kW, "
          f"segments {bounds}")
    d0 = build({"mw_power_W": 0.0}, "loki")
    d1 = build(mw_ov, "loki")
    active, bound = est_active_and_bound(d0.net, d0.mixture)
    by_rid = {r.rid: r for r in d0.net}
    targets = [(r, (1.0 / 3.0, 3.0)) for r in active]
    for rid, alt in CONFLICT_ALTERNATES.items():
        if rid in by_rid:
            targets.append((by_rid[rid], (alt,)))
    rows = []
    print(f"{'variant':>22} {'E_off J':>10} {'dE mJ':>9} {'marg%':>8}")
    t0 = time.time()
    for r, factors in targets:
        for f in factors:
            net2 = [_rep(x, rate=_bracket_rate(x.rate, f))
                    if x.rid == r.rid else x for x in d0.net]
            _swap((d0, d1), net2)
            r0 = run_segmented(d0, bounds, 41)
            r1 = run_segmented(d1, bounds, 41)
            if not (r0["success"] and r1["success"]):
                sys.exit(f"segmented solve failed at {r.rid} x{f:g}")
            e_mw, d_las, marg = _marginal_pair(r0, r1)
            rows.append({"variant": f"{r.rid} x{f:g}", "rid": r.rid,
                         "factor": f, "E_laser_off_J":
                         r0["final_channels_J"]["laser_out"],
                         "dE_laser_J": d_las, "marginal": marg})
            print(f"{rows[-1]['variant']:>22} "
                  f"{r0['final_channels_J']['laser_out']:10.6f} "
                  f"{1e3 * d_las:9.4f} {100 * marg:8.3f}", flush=True)
    _swap((d0, d1), d0.net)
    ms = np.array([r["marginal"] for r in rows])
    surv = int((ms > 0).sum())
    print(f"  {len(rows)} variants in {time.time() - t0:.0f} s: marginal "
          f"min {100 * ms.min():+.2f}% / max {100 * ms.max():+.2f}%; "
          f"POSITIVE in {surv}/{len(rows)}"
          + ("" if surv == len(rows) else " — flipped by: " + ", ".join(
              sorted({r['variant'] for r in rows if r['marginal'] <= 0}))))
    print(f"  element-bound EST rows (structurally zero impact, not run): "
          f"{len(bound)}")
    return {"case": "ne_optimized", "note": note, "window": mw_ov,
            "segments": bounds, "rows": rows,
            "marginal_min": float(ms.min()),
            "marginal_max": float(ms.max()),
            "sign_survival": f"{surv}/{len(rows)}",
            "bound_est_rows": sorted(r.rid for r in bound)}


def main():
    out = {"A_kstep_map": part_kstep_map(),
           "B_neopt_afterglow_brackets": part_neopt_afterglow_brackets()}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print(f"\nreport -> {OUT}")


if __name__ == "__main__":
    main()
