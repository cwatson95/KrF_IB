"""Is the afterglow-reheat gain a numerical artifact? A placebo-calibrated
profile-likelihood / CLs certificate ([0.1.31]).

SCOPE — stated before any number: this test rejects (or fails to reject)
the hypothesis that the twin-pair laser-energy difference is drawn from
the CODE'S OWN NUMERICS-NOISE distribution, calibrated in situ. It
cannot certify the absence of a conceptual error in the physics model —
that is covered by the independent-implementation agreement (REBUILD vs
1DModelKrF), the closed energy ledger, the gate suite, and the
rate-bracket program, not by any statistic of one model's outputs.

Design:
  * SIGNAL ensemble: the [0.1.28]-B ne_optimized afterglow twin
    (115+150 ns window, 500 kW) re-measured under 7 solver settings
    (BDF rtol 1e-6/1e-8/1e-10, atol x0.1, max_step 1 ns, and the
    independent Radau integrator at rtol 1e-6/1e-8), segmented
    integration throughout.
  * NULL calibration: per setting, a PLACEBO window — identical but
    1 W absorbed (linear physical response ~2e-6 of signal, and the
    first-order physical part is additionally subtracted) — whose twin
    difference measures pure numerical differencing noise.
  * Model: x_i ~ N(nu, sigma^2) (signal), y_j ~ N(0, sigma^2) (null),
    common sigma profiled as a nuisance. Exact results used:
      q0   = (n+m) ln(sigma0_hat^2 / sigma_hat^2), Z_asym = sqrt(q0)
      T    = xbar / (s_p / sqrt(n)) ~ t_{n+m-1} under H0 (EXACT)
      s_p^2 = [SS_x + SS_y] / (n + m - 1)
    plus toy Monte Carlo confirmation, a CLs-allowed band for the
    signal strength, the worst case over the 41 [0.1.28] bracket
    variants, and the noise-inflation factor that would degrade the
    discovery to exactly 3 sigma.
  * Self-tests run FIRST (q0 algebra vs direct optimization, t/CLs
    closed forms, synthetic discovery and synthetic null) and abort on
    any failure. A NEGATIVE CONTROL on real data (placebo split
    treated as signal) must come out insignificant.
  * 1-D (SPECKLE): the power ladder 250 kW -> 250 W bounds that code's
    deterministic differencing floor; reported as a floor-ratio bound.

Run from the Kinetic_Modeling root:
    python3 SPECK/REBUILD/tools/artifact_significance.py
Writes SPECK/rebuild_out/mw_benchmarks/artifact_significance.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import optimize, stats

for base in Path(__file__).resolve().parents:
    if (base / "REBUILD" / "MODULES").is_dir():
        sys.path.insert(0, str(base))
        SPECK = base
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.tools.mw_benchmarks import (AFTERGLOW,            # noqa: E402
                                         POWER_PER_M3,
                                         _case_configs,
                                         segment_bounds)

OUT = SPECK / "rebuild_out" / "mw_benchmarks" / "artifact_significance.json"
SPECKLE_RUNS = SPECK.parent / "SPECKLE" / "runs"

P_PLACEBO_W = 1.0
THREE_SIGMA_P = 1.0 - stats.norm.cdf(3.0)          # one-sided, 1.3499e-3


# ===================================================== statistics core
def fit_stats(x, y):
    """Exact Gaussian profile-likelihood + Student-t results for
    x_i ~ N(nu, s^2), y_j ~ N(0, s^2)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n, m = len(x), len(y)
    nu_hat = x.mean()
    ss_x = float(((x - nu_hat) ** 2).sum())
    ss_y = float((y ** 2).sum())
    sig2_hat = (ss_x + ss_y) / (n + m)
    sig2_0 = (float((x ** 2).sum()) + ss_y) / (n + m)
    q0 = (n + m) * np.log(sig2_0 / sig2_hat) if nu_hat > 0 else 0.0
    z_asym = float(np.sqrt(max(q0, 0.0)))
    dof = n + m - 1
    s_p = float(np.sqrt((ss_x + ss_y) / dof))
    T = nu_hat / (s_p / np.sqrt(n))
    p_t = float(stats.t.sf(T, dof))                 # exact under H0
    z_t = float(stats.norm.isf(max(p_t, 1e-300)))
    z_single = nu_hat / s_p                          # no sqrt(n) benefit
    return {"n": n, "m": m, "nu_hat": float(nu_hat), "s_p": s_p,
            "se": s_p / np.sqrt(n), "dof": dof, "q0": float(q0),
            "Z_asym": z_asym, "T": float(T), "p_t_exact": p_t,
            "Z_from_t": z_t, "z_single_measurement": float(z_single)}


def toy_p0(x, y, fit, n_toys=200000, seed=20260827):
    """Toy-MC confirmation of the discovery p-value: generate null
    datasets at sigma = sigma0_hat, recompute T, count exceedances."""
    rng = np.random.default_rng(seed)
    n, m, dof = fit["n"], fit["m"], fit["dof"]
    sig0 = np.sqrt((np.sum(np.asarray(x) ** 2) + np.sum(np.asarray(y) ** 2))
                   / (n + m))
    xt = rng.normal(0.0, sig0, size=(n_toys, n))
    yt = rng.normal(0.0, sig0, size=(n_toys, m))
    nu = xt.mean(axis=1)
    ssx = ((xt - nu[:, None]) ** 2).sum(axis=1)
    ssy = (yt ** 2).sum(axis=1)
    sp = np.sqrt((ssx + ssy) / dof)
    Tt = nu / (sp / np.sqrt(n))
    exceed = int((Tt >= fit["T"]).sum())
    p = exceed / n_toys if exceed else 1.0 / n_toys
    return {"n_toys": n_toys, "exceedances": exceed,
            "p0_toy_bound": p, "Z_toy_bound": float(stats.norm.isf(p)),
            "toy_T_matches_t": float(stats.kstest(
                Tt[:5000], "t", args=(dof,)).pvalue)}


def cls_band(fit, alpha=0.05):
    """95% CLs-allowed band for the signal strength nu (exact central-t
    pivots: under hypothesis nu, (nu_hat-nu)/se ~ t_dof)."""
    tcrit = stats.t.ppf(1.0 - alpha, fit["dof"])
    lo = fit["nu_hat"] - tcrit * fit["se"]
    hi = fit["nu_hat"] + tcrit * fit["se"]
    # CLs denominators (protection terms) at the observed point:
    p_b_up = float(stats.t.cdf(fit["T"], fit["dof"]))    # ~1: no penalty
    return {"nu_lo_95CLs": float(lo), "nu_hi_95CLs": float(hi),
            "zero_below_band_by_se": float((lo - 0.0) / fit["se"]),
            "cls_denominator_up": p_b_up}


def self_tests():
    rng = np.random.default_rng(11)
    # 1) q0 algebra vs direct likelihood optimization
    x = rng.normal(3.0, 1.0, 6)
    y = rng.normal(0.0, 1.0, 5)

    def nll(pars, nu_fixed=None):
        if nu_fixed is None:
            nu, ls = pars
        else:
            nu, ls = nu_fixed, pars[0]
        s2 = np.exp(2 * ls)
        return (0.5 * (len(x) + len(y)) * np.log(2 * np.pi * s2)
                + (np.sum((x - nu) ** 2) + np.sum(y ** 2)) / (2 * s2))

    f = fit_stats(x, y)
    g = optimize.minimize(nll, [1.0, 0.0], method="Nelder-Mead",
                          options={"xatol": 1e-10, "fatol": 1e-12})
    g0 = optimize.minimize(lambda p: nll(p, nu_fixed=0.0), [0.0],
                           method="Nelder-Mead",
                           options={"xatol": 1e-10, "fatol": 1e-12})
    q0_direct = 2 * (g0.fun - g.fun)
    assert abs(q0_direct - f["q0"]) < 1e-5, (q0_direct, f["q0"])
    # 2) t machinery: sf/ppf consistency
    assert abs(stats.t.sf(stats.t.ppf(0.9, 10), 10) - 0.1) < 1e-12
    # 3) synthetic discovery at a known scale: nu = 5 sigma
    xs = rng.normal(5.0, 1.0, 7)
    ys = rng.normal(0.0, 1.0, 7)
    fs = fit_stats(xs, ys)
    assert 5.0 < fs["T"] < 30.0, fs["T"]
    # 4) synthetic null: no manufactured significance
    xn = rng.normal(0.0, 1.0, 7)
    yn = rng.normal(0.0, 1.0, 7)
    fn = fit_stats(xn, yn)
    assert abs(fn["T"]) < 4.0, fn["T"]
    # 5) CLs band width equals the exact central-t formula
    cb = cls_band(fs)
    t95 = stats.t.ppf(0.95, fs["dof"])
    width = cb["nu_hi_95CLs"] - cb["nu_lo_95CLs"]
    assert abs(width - 2 * t95 * fs["se"]) < 1e-12, width
    print("[self-tests] q0-vs-optimizer, t closed forms, synthetic "
          "discovery, synthetic null, CLs band: ALL PASS", flush=True)


# ============================================== 0-D ensemble machinery
SETTINGS = [
    ("BDF rtol 1e-6 (default)", dict(method="BDF", rtol=1e-6)),
    ("BDF rtol 1e-8", dict(method="BDF", rtol=1e-8)),
    ("BDF rtol 1e-10", dict(method="BDF", rtol=1e-10)),
    ("BDF rtol 1e-6, atol x0.1", dict(method="BDF", rtol=1e-6,
                                      atol_scale=0.1)),
    ("BDF rtol 1e-6, max_step 1 ns", dict(method="BDF", rtol=1e-6,
                                          max_step=1e-9)),
    ("Radau rtol 1e-6", dict(method="Radau", rtol=1e-6)),
    ("Radau rtol 1e-8", dict(method="Radau", rtol=1e-8)),
]


def solve_segmented(dev, bounds_ns, method="BDF", rtol=1e-6,
                    atol_scale=1.0, max_step=2e-9):
    """Segmented integration (the 0.1.21 cure) with the solver exposed;
    ledger channels ride as state, so the final segment carries the
    totals. Returns final_channels_J."""
    from scipy.integrate import solve_ivp
    atol = np.full(dev.lay.size, 1e3)
    ca = dev.lay.cells_view(atol)
    ca[dev.lay.i_ue] = 1e-10
    ca[dev.lay.i_tg] = 1e-6
    atol[dev.lay.i_ph] = 1e-3
    atol[dev.lay.i_ledger0:] = 1e-12
    atol *= atol_scale
    y = dev.initial_state()
    y0 = y
    t_last = 0.0
    sol = None
    for b in bounds_ns:
        t_eval = np.linspace(t_last * 1e-9, b * 1e-9, 5)
        sol = solve_ivp(dev.rhs, (t_last * 1e-9, b * 1e-9), y,
                        method=method, rtol=rtol, atol=atol,
                        t_eval=t_eval, max_step=max_step)
        if not sol.success:
            raise RuntimeError(f"{method} rtol={rtol} failed in segment "
                               f"ending {b} ns: {sol.message}")
        y = sol.y[:, -1]
        t_last = b
    res = dev._collect(sol, y0, 0.0, False)
    return res["final_channels_J"]


def part_0d():
    print("\n[0-D] ne_optimized afterglow twin under the solver ensemble "
          "(signal 500 kW + placebo 1 W per setting)", flush=True)
    build, note = _case_configs()["ne_optimized"]
    probe = build({"mw_power_W": 0.0}, "loki")
    p = probe.p
    beam_end = float(p["beam_rise_ns"] + p["beam_flat_ns"])
    t_end = float(p["t_end_ns"])
    cell = float(p["Lx_m"] * p["Ly_m"] * p["Lz_m"])
    p_sig = POWER_PER_M3 * cell
    win = {"mw_t0_ns": beam_end + AFTERGLOW["delay_ns"],
           "mw_flat_ns": AFTERGLOW["flat_ns"],
           "mw_rise_ns": AFTERGLOW["rise_ns"]}
    bounds = segment_bounds(beam_end, t_end)
    print(f"    {note}; window {win['mw_t0_ns']:.0f}+{win['mw_flat_ns']:.0f}"
          f" ns, signal {p_sig / 1e3:.0f} kW, placebo {P_PLACEBO_W:g} W, "
          f"segments {bounds}", flush=True)
    rows = []
    for name, kw in SETTINGS:
        t0 = time.time()
        try:
            offs = solve_segmented(build({"mw_power_W": 0.0}, "loki"),
                                   bounds, **kw)
            ons = solve_segmented(build({**win, "mw_power_W": p_sig},
                                        "loki"), bounds, **kw)
            pls = solve_segmented(build({**win,
                                         "mw_power_W": P_PLACEBO_W},
                                        "loki"), bounds, **kw)
        except RuntimeError as exc:
            print(f"    {name}: EXCLUDED ({exc})", flush=True)
            rows.append({"setting": name, "failed": str(exc)})
            continue
        d_sig = ons["laser_out"] - offs["laser_out"]
        d_pl = pls["laser_out"] - offs["laser_out"]
        e_mw = ons["drive_in"] - offs["drive_in"]
        d_pl_corr = d_pl - (P_PLACEBO_W / p_sig) * d_sig
        rows.append({"setting": name, "dE_signal_J": float(d_sig),
                     "E_MW_J": float(e_mw),
                     "dE_placebo_J": float(d_pl),
                     "dE_placebo_linear_corrected_J": float(d_pl_corr),
                     "wall_s": time.time() - t0})
        print(f"    {name}: dE_sig {1e3 * d_sig:+.5f} mJ, dE_placebo "
              f"{1e3 * d_pl:+.3e} mJ (corr {1e3 * d_pl_corr:+.3e}) "
              f"[{time.time() - t0:.0f} s]", flush=True)
    return {"window": win, "p_signal_W": p_sig, "segments": bounds,
            "rows": rows}


def analyze_0d(part):
    ok = [r for r in part["rows"] if "failed" not in r]
    if len(ok) < 5:
        sys.exit("fewer than 5 solver settings succeeded — cannot certify")
    x = [1e3 * r["dE_signal_J"] for r in ok]                     # mJ
    y = [1e3 * r["dE_placebo_linear_corrected_J"] for r in ok]   # mJ
    fit = fit_stats(x, y)
    toys = toy_p0(x, y, fit)
    band = cls_band(fit)
    t3 = stats.t.ppf(1.0 - THREE_SIGMA_P, fit["dof"])
    inflation_to_3sigma = fit["T"] / t3
    # worst case over the 41 bracket variants of [0.1.28]:
    sd = json.loads((SPECK / "rebuild_out" / "mw_benchmarks" /
                     "sign_decidability.json").read_text())
    d_min = 1e3 * min(r["dE_laser_J"]
                      for r in sd["B_neopt_afterglow_brackets"]["rows"])
    z_worst = d_min / fit["s_p"]
    # negative control on REAL placebo data (split): must be insignificant
    neg = fit_stats(y[: (len(y) + 1) // 2], y[(len(y) + 1) // 2:])
    return {"signal_mJ": x, "null_mJ": y, "fit": fit, "toys": toys,
            "cls_band_mJ": band,
            "t_crit_3sigma": float(t3),
            "noise_inflation_factor_to_3sigma": float(inflation_to_3sigma),
            "worst_bracket_variant_dE_mJ": float(d_min),
            "z_worst_bracket_single_measurement": float(z_worst),
            "negative_control_T": float(neg["T"]),
            "negative_control_p": float(neg["p_t_exact"])}


# ===================================================== 1-D floor bound
def part_1d():
    print("\n[1-D] SPECKLE power ladder: deterministic differencing floor",
          flush=True)
    pts = []
    for name in ("ag250k", "ag27k", "pb_p2k5", "pb_p250"):
        f = SPECKLE_RUNS / f"{name}_afterglow_report.json"
        if not f.exists():
            print(f"    {name}: report missing, skipped", flush=True)
            continue
        d = json.loads(f.read_text())
        pts.append({"name": name, "P_W": d["config"]["mw_power"],
                    "dE_mJ": 1e3 * d["dE_laser_J"],
                    "marginal": d["marginal"]})
    pts.sort(key=lambda r: -r["P_W"])
    for r in pts:
        print(f"    P {r['P_W']:9.3e} W: dE {r['dE_mJ']:+.6f} mJ "
              f"(marginal {100 * r['marginal']:+.4f}%)", flush=True)
    out = {"points": pts}
    if len(pts) >= 3 and pts[-1]["P_W"] < 1e3:
        floor = abs(pts[-1]["dE_mJ"])
        sig = pts[0]["dE_mJ"]
        out["floor_bound_mJ"] = floor
        out["signal_mJ"] = sig
        out["signal_over_floor_bound"] = sig / floor if floor > 0 else float("inf")
        scaled = all(abs(pts[i + 1]["dE_mJ"]) < abs(pts[i]["dE_mJ"])
                     for i in range(len(pts) - 1))
        out["monotone_decrease_with_P"] = bool(scaled)
    return out


def main():
    self_tests()
    part = part_0d()
    ana = analyze_0d(part)
    oned = part_1d()

    f = ana["fit"]
    print("\n================ CERTIFICATE (scope in the docstring) "
          "================")
    print(f"0-D signal: nu_hat = {f['nu_hat']:+.4f} mJ over n={f['n']} "
          f"solver settings; placebo-calibrated noise s_p = "
          f"{f['s_p']:.2e} mJ (m={f['m']} placebo twins)")
    print(f"  exact Student-t discovery: T = {f['T']:.1f} "
          f"(dof {f['dof']}), p0 = {f['p_t_exact']:.3g}  ->  "
          f"Z = {f['Z_from_t']:.1f} sigma")
    print(f"  profile-likelihood q0 = {f['q0']:.1f} -> Z_asym = "
          f"{f['Z_asym']:.1f};  conservative single-measurement z = "
          f"{f['z_single_measurement']:.1f}")
    print(f"  toys: {ana['toys']['exceedances']} exceedances in "
          f"{ana['toys']['n_toys']} null toys -> p0 < "
          f"{ana['toys']['p0_toy_bound']:.1e} (Z > "
          f"{ana['toys']['Z_toy_bound']:.2f} at toy resolution)")
    print(f"  95% CLs-allowed signal band: [{ana['cls_band_mJ']['nu_lo_95CLs']:.4f}, "
          f"{ana['cls_band_mJ']['nu_hi_95CLs']:.4f}] mJ; nu=0 sits "
          f"{ana['cls_band_mJ']['zero_below_band_by_se']:.1f} s.e. below it")
    print(f"  worst of the 41 bracket variants: dE = "
          f"{ana['worst_bracket_variant_dE_mJ']:+.3f} mJ -> z = "
          f"{ana['z_worst_bracket_single_measurement']:.1f}")
    print(f"  the noise scale would have to be underestimated x"
          f"{ana['noise_inflation_factor_to_3sigma']:.0f} to reduce the "
          f"discovery to 3 sigma")
    neg_ok = abs(ana["negative_control_T"]) < 3.5
    print(f"  negative control (placebo split as signal): T = "
          f"{ana['negative_control_T']:+.2f}, p = "
          f"{ana['negative_control_p']:.2f} -> "
          f"{'PASS (insignificant, as required)' if neg_ok else 'FAIL — noise model questioned; rely on the single-measurement z'}")
    if "signal_over_floor_bound" in oned:
        print(f"1-D: differencing floor bounded by |dE(250 W)| = "
              f"{oned['floor_bound_mJ']:.2e} mJ; the 250-kW signal "
              f"exceeds it x{oned['signal_over_floor_bound']:.0f}; "
              f"monotone with P: {oned['monotone_decrease_with_P']}")
    print("NOT covered by this statistic: shared conceptual model error "
          "— see the independent-implementation agreement, ledger "
          "closure, suite, and bracket program.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"part_0d": part, "analysis_0d": ana,
                               "part_1d": oned}, indent=1))
    print(f"\nreport -> {OUT}")


if __name__ == "__main__":
    main()
