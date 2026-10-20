#!/usr/bin/env python3
"""Is the co-pump gain numerically real, and how sure? — the verification.

Owner ask (2026-08-21): "Is there a way to know that the positive increase of
efficiency is numerically correct? Is there a way to verify statistically to
ensure no false negatives or false positives?"

The engine is deterministic — there is no shot noise — so "statistics" here
means three things, each measured on the reference co-pump point (preset
krf_hyper_mw: 2 atm, 100 A cm^-2, 250 kW on the lasing plateau 35–74 ns):

  1. NUMERICAL CONVERGENCE of the twin difference dE_laser: the same twin
     pair re-solved at rtol 1e-6 / 1e-8 / 1e-10, max_step 2 / 1 ns,
     tighter atol, and with an INDEPENDENT integrator (Radau vs the
     engine's BDF). A real effect is invariant; an artifact moves.
  2. THE NULL FLOOR: the spread of the BASELINE laser energy across those
     solver settings is the absolute numerical noise; the spread of dE
     across settings is the DIFFERENTIAL noise (the twin structure cancels
     common-mode error). z = dE / sigma_diff is the false-positive
     statistic; the minimum detectable dE = 3 sigma_diff sets the
     false-negative bound — the smallest microwave power whose effect the
     engine could resolve at this point (P_min = 3 sigma_diff / (m t_win)).
  3. PHYSICS ROBUSTNESS: the marginal re-measured under the network's own
     rate uncertainty — every active EST row x1/3 and x3 and every recorded
     LIT-conflict alternate (the S5 bracket machinery of scorecard.py, a
     variant CompiledNetwork per run, BOTH twins on the same variant) —
     and on the INDEPENDENT electron-impact rate tier (KINETIC =
     PIC_MCC_3D-measured rates vs LOKI = LoKI-MC). Sign survival across
     all variants is the physics false-positive check; a variant flipping
     the sign would name the rate the conclusion hangs on.
  4. EXPERIMENTAL POWER ANALYSIS (arithmetic): interleaved on/off shot
     pairs needed for alpha 0.05 (two-sided) / power 0.8 at a given shot
     jitter — N = (z_a + z_b)^2 x 2 sigma^2 / delta^2.

Outputs a table + rebuild_out/copump_scan/copump_verify.json. Predictions
verified against themselves, not against data — a numerical-correctness
and robustness statement, stated as such.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/copump_verify.py           (~15-25 min)
    python3 SPECK/REBUILD/tools/copump_verify.py --parts numeric,tier
"""

import argparse
import json
import sys
import time
from dataclasses import replace as _rep
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

from REBUILD.tools.hybrid_timing_study import build_device  # noqa: E402

PLATEAU = {"mw_t0_ns": 35.0, "mw_flat_ns": 35.0, "mw_rise_ns": 2.0}
BEAM = {"pressure_atm": 2.0, "beam_j_A_cm2": 100.0, "t_end_ns": 250.0}
POWER = 2.5e5
T_WIN_S = (35.0 + 2.0) * 1e-9          # trapezoid area of the window


def solve(dev, method="BDF", rtol=None, max_step=2e-9, atol_scale=1.0,
          n_save=161):
    """Device.run's solve with the method/tolerance exposed (mirrors run()
    exactly otherwise — same atol construction, same _collect)."""
    from scipy.integrate import solve_ivp
    p = dev.p
    t_end = p["t_end_ns"] * 1e-9
    y0 = dev.initial_state()
    t_eval = np.linspace(0.0, t_end, int(n_save))
    atol = np.full(dev.lay.size, 1e3)
    ca = dev.lay.cells_view(atol)
    ca[dev.lay.i_ue] = 1e-10
    ca[dev.lay.i_tg] = 1e-6
    atol[dev.lay.i_ph] = 1e-3
    atol[dev.lay.i_ledger0:] = 1e-12
    atol *= atol_scale
    t0 = time.time()
    sol = solve_ivp(dev.rhs, (0.0, t_end), y0, method=method,
                    rtol=(p["rtol"] if rtol is None else rtol), atol=atol,
                    t_eval=t_eval, max_step=max_step)
    res = dev._collect(sol, y0, time.time() - t0, False)
    if not res["success"]:
        sys.exit(f"solve failed: {method} rtol={rtol}")
    return res


def twin(preset, mw_ov, beam_ov, power, tier="loki", **solve_kw):
    off = dict(beam_ov)
    off["mw_power_W"] = 0.0
    off["electron_tier"] = tier
    on = dict(beam_ov)
    on.update(mw_ov)
    on["mw_power_W"] = power
    on["electron_tier"] = tier
    d0 = build_device(preset, off)
    d1 = build_device(preset, on)
    r0 = solve(d0, **solve_kw)
    r1 = solve(d1, **solve_kw)
    e_mw = r1["final_channels_J"]["drive_in"] - r0["final_channels_J"]["drive_in"]
    d_las = r1["final_channels_J"]["laser_out"] - r0["final_channels_J"]["laser_out"]
    return {"E_laser_off_J": r0["final_channels_J"]["laser_out"],
            "E_laser_on_J": r1["final_channels_J"]["laser_out"],
            "E_MW_J": e_mw, "dE_laser_J": d_las,
            "marginal": d_las / e_mw if e_mw > 0 else 0.0,
            "devs": (d0, d1)}


def part_numeric(preset):
    settings = [
        ("BDF rtol 1e-6 (engine default)", dict(method="BDF", rtol=1e-6)),
        ("BDF rtol 1e-8", dict(method="BDF", rtol=1e-8)),
        ("BDF rtol 1e-10", dict(method="BDF", rtol=1e-10)),
        ("BDF rtol 1e-6, max_step 1 ns", dict(method="BDF", rtol=1e-6,
                                              max_step=1e-9)),
        ("BDF rtol 1e-6, atol x0.1", dict(method="BDF", rtol=1e-6,
                                          atol_scale=0.1)),
        ("Radau rtol 1e-6 (independent integrator)",
         dict(method="Radau", rtol=1e-6)),
        ("Radau rtol 1e-8", dict(method="Radau", rtol=1e-8)),
    ]
    rows = []
    print("\n[1] numerical convergence of the twin difference (plateau point)")
    print(f"{'setting':>42} {'E_off J':>10} {'dE mJ':>9} {'marg%':>7}")
    for name, kw in settings:
        t = twin(preset, PLATEAU, BEAM, POWER, **kw)
        rows.append({"setting": name, **{k: v for k, v in t.items()
                                         if k != "devs"}})
        print(f"{name:>42} {t['E_laser_off_J']:10.6f} "
              f"{1e3 * t['dE_laser_J']:9.4f} {100 * t['marginal']:7.3f}",
              flush=True)
    de = np.array([r["dE_laser_J"] for r in rows])
    e0 = np.array([r["E_laser_off_J"] for r in rows])
    m = np.array([r["marginal"] for r in rows])
    sig_diff = float(de.std(ddof=1))
    sig_abs = float(e0.std(ddof=1))
    z = float(de.mean() / sig_diff) if sig_diff > 0 else float("inf")
    de_min = 3.0 * sig_diff
    p_min = de_min / (m.mean() * T_WIN_S) if m.mean() > 0 else float("inf")
    print(f"  dE mean {1e3 * de.mean():.4f} mJ, spread (sigma_diff) "
          f"{1e3 * sig_diff:.5f} mJ -> z = {z:.0f}; baseline spread "
          f"(sigma_abs) {1e3 * sig_abs:.5f} mJ; all signs "
          f"{'POSITIVE' if (de > 0).all() else 'MIXED'}; "
          f"min detectable dE {1e3 * de_min:.5f} mJ -> P_min "
          f"{p_min:.3g} W at this marginal/window")
    return {"rows": rows, "sigma_diff_J": sig_diff, "sigma_abs_J": sig_abs,
            "z": z, "dE_min_J": de_min, "P_min_W": p_min,
            "all_positive": bool((de > 0).all())}


def part_brackets(preset):
    from REBUILD.MODULES.chemistry import CompiledNetwork
    from REBUILD.MODULES.scorecard import (CONFLICT_ALTERNATES, _bracket_rate,
                                           est_active_and_bound)
    print("\n[3] physics robustness: marginal under EST x1/3..x3 and "
          "LIT-conflict alternates (both twins on the same variant)")
    off = dict(BEAM)
    off["mw_power_W"] = 0.0
    on = dict(BEAM)
    on.update(PLATEAU)
    on["mw_power_W"] = POWER
    d0 = build_device(preset, off)
    d1 = build_device(preset, on)
    active, bound = est_active_and_bound(d0.net, d0.mixture)
    by_rid = {r.rid: r for r in d0.net}
    targets = [(r, (1.0 / 3.0, 3.0)) for r in active]
    for rid, alt in CONFLICT_ALTERNATES.items():
        if rid in by_rid:
            targets.append((by_rid[rid], (alt,)))
    rows = []
    print(f"{'variant':>22} {'E_off J':>10} {'dE mJ':>9} {'marg%':>7}")
    for r, factors in targets:
        for f in factors:
            rows2 = [_rep(x, rate=_bracket_rate(x.rate, f))
                     if x.rid == r.rid else x for x in d0.net]
            for d in (d0, d1):
                d.compiled = CompiledNetwork(rows2, d.lay, d.tiers, d.optics)
                d.compiled.set_photon_energy(d.hv_eV)
            r0 = solve(d0)
            r1 = solve(d1)
            e_mw = (r1["final_channels_J"]["drive_in"]
                    - r0["final_channels_J"]["drive_in"])
            d_las = (r1["final_channels_J"]["laser_out"]
                     - r0["final_channels_J"]["laser_out"])
            marg = d_las / e_mw if e_mw > 0 else 0.0
            name = f"{r.rid} x{f:g}"
            rows.append({"variant": name, "rid": r.rid, "factor": f,
                         "E_laser_off_J": r0["final_channels_J"]["laser_out"],
                         "dE_laser_J": d_las, "marginal": marg})
            print(f"{name:>22} {r0['final_channels_J']['laser_out']:10.6f} "
                  f"{1e3 * d_las:9.4f} {100 * marg:7.3f}", flush=True)
    for d in (d0, d1):
        d.compiled = CompiledNetwork(d.net, d.lay, d.tiers, d.optics)
        d.compiled.set_photon_energy(d.hv_eV)
    ms = np.array([r["marginal"] for r in rows])
    print(f"  {len(rows)} variants: marginal min {100 * ms.min():.2f}% / "
          f"max {100 * ms.max():.2f}%, sign survival "
          f"{int((ms > 0).sum())}/{len(rows)}; element-bound EST rows "
          f"(structurally zero impact, not run): {len(bound)}")
    return {"rows": rows, "marginal_min": float(ms.min()),
            "marginal_max": float(ms.max()),
            "sign_survival": f"{int((ms > 0).sum())}/{len(rows)}",
            "bound_est_rows": sorted(r.rid for r in bound)}


def part_tier(preset):
    print("\n[3b] independent rate tier: KINETIC (PIC_MCC_3D-measured) vs "
          "LOKI (LoKI-MC)")
    out = {}
    for tier in ("loki", "kinetic"):
        t = twin(preset, PLATEAU, BEAM, POWER, tier=tier)
        out[tier] = {k: v for k, v in t.items() if k != "devs"}
        print(f"  {tier:>8}: E_off {t['E_laser_off_J']:.6f} J  dE "
              f"{1e3 * t['dE_laser_J']:+.4f} mJ  marginal "
              f"{100 * t['marginal']:+.3f}%", flush=True)
    return out


def part_power_analysis():
    print("\n[4] experimental power analysis (interleaved on/off shot "
          "pairs, alpha 0.05 two-sided, power 0.8)")
    za, zb = 1.96, 0.84
    rows = []
    for sigma in (0.01, 0.02, 0.05):
        for label, delta in (("250 kW plateau, +0.15%", 0.0015),
                             ("1 MW, +0.6%", 0.006),
                             ("3 MW @1 atm, +3%", 0.03)):
            n = (za + zb) ** 2 * 2.0 * sigma ** 2 / delta ** 2
            rows.append({"shot_jitter": sigma, "case": label,
                         "effect": delta, "pairs": n})
            print(f"  jitter {100 * sigma:3.0f}%  {label:>24}: "
                  f"{n:8.0f} shot pairs", flush=True)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_hyper_mw")
    ap.add_argument("--parts", default="numeric,brackets,tier,power")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    parts = {s.strip() for s in args.parts.split(",") if s.strip()}
    out = {"preset": args.preset, "point": {"window": PLATEAU,
                                            "beam": BEAM, "P_W": POWER}}
    if "numeric" in parts:
        out["numeric"] = part_numeric(args.preset)
    if "brackets" in parts:
        out["brackets"] = part_brackets(args.preset)
    if "tier" in parts:
        out["tier"] = part_tier(args.preset)
    if "power" in parts:
        out["power_analysis"] = part_power_analysis()
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "copump_scan")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "copump_verify.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
