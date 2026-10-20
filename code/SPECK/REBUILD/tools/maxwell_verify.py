"""The MAXWELL verification — how much of each electron-impact rate, and of
the beam-mode device answers, is the EEDF-shape assumption (0.1.37).

Rate level (every lane): for a fill's LOKI table set, at each tabulated
<eps> point of each bucket, k_LOKI, k_MAXWELL(Te = <eps>/1.5), their ratio,
and the two-Maxwellian bracket [min, max] over a STATED (tail fraction,
tail temperature) grid at fixed <eps>. The discharge lane (Chen ArF) is
covered the only way a field-free Maxwellian can be: the same ratios
indexed by E/N through the swarm table's own <eps>(E/N) (k_by_en against
the Maxwellian at that mean energy). The tail grid is a sweep coordinate
for the bracket, stated on the command line, never an engine parameter.

Device level (beam drive only — the engine refuses the circuit lane on
this tier): presets and anchors under LOKI (default, relax bridge), LOKI
with the taper stand-in (no bridge, so the bridge's own share is visible),
and MAXWELL. Populations, E_L, eta, peaks: misses stated, nothing tuned.

    python3 SPECK/REBUILD/tools/maxwell_verify.py                # rate level, 3 fills
    python3 SPECK/REBUILD/tools/maxwell_verify.py --device       # + presets/ne_optimized
    python3 SPECK/REBUILD/tools/maxwell_verify.py --population   # + NRL population rows
    python3 SPECK/REBUILD/tools/maxwell_verify.py --kannari      # + Kannari anchor (slow)
    python3 SPECK/REBUILD/tools/maxwell_verify.py --tail 0.05:5 --tail 0.2:3

Writes rebuild_out/maxwell_verify/report.json. Anchors score; they never
steer — this tool changes nothing in the engine.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
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

from REBUILD.MODULES.loki_tier import (RUN_NAMES_AR_FILL, RUN_NAMES_NE_FILL,  # noqa: E402
                                       RUN_NAMES_ARF_CHEN, LokiRates)
from REBUILD.MODULES.maxwell_tier import MaxwellRates  # noqa: E402

FILLS = {
    "ar": ({"Ar": 0.949, "Kr": 0.05, "F2": 0.001}, RUN_NAMES_AR_FILL, "krf"),
    "ne": ({"Ne": 0.949, "Kr": 0.05, "F2": 0.001}, RUN_NAMES_NE_FILL, "krf"),
    "chen": ({"Ne": 0.964, "Ar": 0.035, "F2": 0.001}, RUN_NAMES_ARF_CHEN, "arf_chen"),
}
DEFAULT_TAILS = [(a, t) for a in (0.01, 0.03, 0.1) for t in (2.0, 4.0, 8.0)]


# ------------------------------------------------------------------ rate level
def rate_rows(fill: str, tails) -> list:
    mixture, run_names, _mode = FILLS[fill]
    lok = LokiRates(run_names)
    mx = MaxwellRates(mixture)
    rows = []
    for b in lok.buckets:
        xs, ys = lok.curves[b]
        for x, y in zip(xs, ys):
            eps = float(np.exp(x))
            te = eps / 1.5
            k_l = float(np.exp(y))
            k_m = float(mx.k(b, te))
            ratios = []
            for a, tt in tails:
                if a * tt < te:
                    ratios.append(float(mx.k_two_t(b, te, a, tt)) / k_l)
            rows.append({"fill": fill, "index": "eps", "bucket": b,
                         "eps_eV": eps, "k_loki": k_l, "k_maxwell": k_m,
                         "ratio": k_m / k_l,
                         "two_t_min": (min(ratios) if ratios else None),
                         "two_t_max": (max(ratios) if ratios else None),
                         "n_two_t": len(ratios)})
        # the discharge view: the same bucket indexed by E/N
        cur = lok.en_curves.get(b)
        if cur is None:
            continue
        for x, y in zip(*cur):
            en = float(np.exp(x))
            eps = lok.mean_energy_at_en(en)
            te = eps / 1.5
            k_l = float(np.exp(y))
            k_m = float(mx.k(b, te))
            rows.append({"fill": fill, "index": "en", "bucket": b,
                         "en_td": en, "eps_eV": eps, "k_loki": k_l,
                         "k_maxwell": k_m, "ratio": k_m / k_l})
    return rows


def _fmt_ratio(r):
    return "   --   " if r is None else f"{r:8.3g}"


def format_rate_rows(rows: list, per_bucket: int = 5) -> str:
    out = []
    fills = sorted({r["fill"] for r in rows}, key=list(FILLS).index)
    for fill in fills:
        mixture = FILLS[fill][0]
        out.append(f"== {fill} fill {mixture}: k_MAXWELL / k_LOKI at the table's "
                   "own <eps> points; [two-Maxwellian bracket at fixed <eps>]")
        out.append(f"  {'bucket':8s} {'<eps> eV':>9s} {'k_LOKI':>10s} {'M/L':>8s} "
                   f"{'2T min':>8s} {'2T max':>8s}")
        by_b = {}
        for r in rows:
            if r["fill"] == fill and r["index"] == "eps":
                by_b.setdefault(r["bucket"], []).append(r)
        for b, rs in by_b.items():
            idx = np.unique(np.linspace(0, len(rs) - 1, min(per_bucket, len(rs))).astype(int))
            for i in idx:
                r = rs[i]
                out.append(f"  {b:8s} {r['eps_eV']:9.3g} {r['k_loki']:10.3e} "
                           f"{r['ratio']:8.3g} {_fmt_ratio(r['two_t_min'])} "
                           f"{_fmt_ratio(r['two_t_max'])}")
        en_rows = [r for r in rows if r["fill"] == fill and r["index"] == "en"]
        if en_rows:
            out.append(f"  -- E/N view (the discharge lane): M/L at the table's E/N "
                       "points through <eps>(E/N)")
            by_b = {}
            for r in en_rows:
                by_b.setdefault(r["bucket"], []).append(r)
            for b, rs in by_b.items():
                idx = np.unique(np.linspace(0, len(rs) - 1, min(4, len(rs))).astype(int))
                out.append("  " + b.ljust(8) + "  " + "  ".join(
                    f"{rs[i]['en_td']:.3g} Td/{rs[i]['eps_eV']:.2g} eV: {rs[i]['ratio']:.3g}"
                    for i in idx))
    return "\n".join(out)


def rate_summary(rows: list) -> dict:
    """Per fill/bucket: the shape ratio at the LOW and HIGH ends of the
    table and the two-Maxwellian envelope — the numbers a changelog quotes."""
    s = {}
    for r in rows:
        if r["index"] != "eps":
            continue
        key = f"{r['fill']}:{r['bucket']}"
        d = s.setdefault(key, {"eps_lo": None, "ratio_lo": None, "eps_hi": None,
                               "ratio_hi": None, "two_t_min": None, "two_t_max": None})
        if d["eps_lo"] is None or r["eps_eV"] < d["eps_lo"]:
            d["eps_lo"], d["ratio_lo"] = r["eps_eV"], r["ratio"]
        if d["eps_hi"] is None or r["eps_eV"] > d["eps_hi"]:
            d["eps_hi"], d["ratio_hi"] = r["eps_eV"], r["ratio"]
        for k in ("two_t_min", "two_t_max"):
            if r[k] is not None:
                d[k] = r[k] if d[k] is None else (min if k == "two_t_min" else max)(d[k], r[k])
    return s


# ---------------------------------------------------------------- device level
TIERS = (("loki", None), ("loki-taper", "taper"), ("maxwell", None))


def _observables(dev, res) -> dict:
    fc = res["final_channels_J"]
    dens = res["dens_avg"]
    ob = {"E_L_J": float(fc["laser_out"]), "E_dep_J": float(fc["drive_in"]),
          "eta": float(res["eta_intrinsic"]),
          "ne_peak_m3": float(np.max(res["ne_m3"])),
          "Te_peak_eV": float(np.max(res["Te_eV"])),
          "Te_end_eV": float(res["Te_eV"][-1]),
          "g0_peak_m1": float(np.max(res["g0_m1"])),
          # the engine's own gauge: fraction of the pulse with Te inside the
          # LOKI tables' <eps> range — below it BOTH tiers extrapolate by the
          # same flat+taper rule, so a low fraction says the device answer is
          # the extrapolation's, not the in-table shape's
          "te_in_table_frac": float(res.get("validity", {}).get("te_in_table_frac", float("nan"))),
          "wall_s": float(res["wall_s"])}
    for sp in ("F-", "Ar*", "Ar**", "Kr*", "Kr**", "KrF*(B)", "ArF*", "Ne*"):
        if sp in dens:
            ob[f"{sp}_peak_m3"] = float(np.max(dens[sp]))
    if "F2" in dens and dens["F2"][0] > 0:
        ob["F2_burn"] = float(1.0 - dens["F2"][-1] / dens["F2"][0])
    tier = dev.tiers
    ob["tier"] = type(tier).__name__ + (
        f" [{len(getattr(tier, 'maxwell_buckets', []))} Maxwell buckets]"
        if hasattr(tier, "maxwell_buckets") else "")
    return ob


def device_cases(names) -> dict:
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    from REBUILD.run import PRESETS
    from REBUILD.tools.mw_benchmarks import _case_configs
    cfg = _case_configs()
    out = {}
    for name in names:
        out[name] = {}
        for tier, sub in TIERS:
            t0 = time.time()
            base_tier = "loki" if tier.startswith("loki") else tier
            if name in PRESETS:
                P = engine_parameters()
                P.set_many(PRESETS[name]["params"], source="problem")
                P.set("electron_tier", base_tier, source="cli")
                if sub:
                    P.set("beam_subtable", sub, source="cli")
                dev = Device(P, mixture=PRESETS[name].get("mixture"))
            else:
                build, _notes = cfg[name]
                ov = {"mw_power_W": 0.0}
                if sub:
                    ov["beam_subtable"] = sub
                dev = build(ov, base_tier)
            try:
                res = dev.run(n_save=121)
            except Exception as e:                        # report, do not hide
                out[name][tier] = {"error": repr(e)[:300]}
                print(f"  {name:14s} {tier:11s} ERROR {e!r:.200}", flush=True)
                continue
            ob = _observables(dev, res)
            out[name][tier] = ob
            print(f"  {name:14s} {tier:11s} E_L {ob['E_L_J']:.4g} J  eta {ob['eta']:.4g}  "
                  f"ne_pk {ob['ne_peak_m3']:.3g}  Te_pk {ob['Te_peak_eV']:.3g} eV  "
                  f"({time.time() - t0:.0f} s)", flush=True)
    return out


def format_device(out: dict) -> str:
    keys = ["E_L_J", "eta", "ne_peak_m3", "Te_peak_eV", "te_in_table_frac", "g0_peak_m1",
            "F-_peak_m3", "Ar**_peak_m3", "Ar*_peak_m3", "Kr**_peak_m3", "F2_burn"]
    lines = []
    for name, per in out.items():
        lines.append(f"== {name}: ratio to LOKI (default) per tier "
                     "(te_in_table_frac printed as a value)")
        ref = per.get("loki", {})
        head = f"  {'observable':14s}" + "".join(f"{t:>14s}" for t, _ in TIERS)
        lines.append(head)
        for k in keys:
            if k not in ref:
                continue
            row = f"  {k:14s}"
            for t, _ in TIERS:
                v = per.get(t, {}).get(k)
                if v is None or (ref[k] == 0 and k != "te_in_table_frac"):
                    row += f"{'--':>14s}"
                elif k == "te_in_table_frac" or t == "loki":
                    row += f"{v:14.4g}"
                else:
                    row += f"{v / ref[k]:14.4g}"
            lines.append(row)
    return "\n".join(lines)


def anchor_rows(fn_name: str) -> dict:
    from REBUILD.MODULES import scorecard
    fn = getattr(scorecard, fn_name)
    out = {}
    for tier, sub in TIERS:
        t0 = time.time()
        base_tier = "loki" if tier.startswith("loki") else tier
        try:
            rows = fn(tier=base_tier, subtable=sub)
        except Exception as e:
            out[tier] = {"error": repr(e)[:300]}
            print(f"  {fn_name} {tier}: ERROR {e!r:.200}", flush=True)
            continue
        out[tier] = [list(r) if isinstance(r, (list, tuple)) else r for r in rows]
        print(f"  {fn_name} under {tier} ({time.time() - t0:.0f} s):", flush=True)
        for r in rows:
            print("    " + (" | ".join(str(x) for x in r) if isinstance(r, (list, tuple))
                            else str(r)), flush=True)
    return out


# ------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fill", choices=sorted(FILLS), action="append",
                    help="fill(s) for the rate-level comparison (default: all)")
    ap.add_argument("--tail", action="append", metavar="FRACTION:T_eV",
                    help="two-Maxwellian sweep point at fixed <eps> (repeatable; "
                         "default grid 0.01/0.03/0.1 x 2/4/8 eV)")
    ap.add_argument("--device", action="store_true",
                    help="beam-drive presets + ne_optimized under loki / "
                         "loki-taper / maxwell")
    ap.add_argument("--population", action="store_true",
                    help="the NRL population benchmark rows under the three tiers")
    ap.add_argument("--kannari", action="store_true",
                    help="the Kannari anchor rows under the three tiers (slow)")
    ap.add_argument("--out", default=str(ROOT / "rebuild_out" / "maxwell_verify"))
    args = ap.parse_args(argv)

    tails = DEFAULT_TAILS
    if args.tail:
        tails = []
        for s in args.tail:
            a, t = s.split(":")
            tails.append((float(a), float(t)))
    fills = args.fill or list(FILLS)
    report = {"tails": tails, "fills": fills}

    rows = []
    for f in fills:
        rows += rate_rows(f, tails)
    print(format_rate_rows(rows))
    report["rate_rows"] = rows
    report["rate_summary"] = rate_summary(rows)
    print()
    print("summary (M/L at table low end -> high end; two-Maxwellian envelope):")
    for key, d in report["rate_summary"].items():
        print(f"  {key:14s} {d['eps_lo']:6.3g} eV: {d['ratio_lo']:8.3g}   "
              f"{d['eps_hi']:6.3g} eV: {d['ratio_hi']:8.3g}   "
              f"2T [{_fmt_ratio(d['two_t_min'])}, {_fmt_ratio(d['two_t_max'])}]")

    if args.device:
        print("\n== device level (beam drive): presets + the ne_optimized lasing cell")
        report["device"] = device_cases(["krf_beam_0d", "krf_ne_0d", "ne_optimized"])
        print(format_device(report["device"]))
    if args.population:
        print("\n== NRL Orestes population benchmark (electra_population_rows)")
        report["population"] = anchor_rows("electra_population_rows")
    if args.kannari:
        print("\n== Kannari 1985 anchor (kannari_device_rows)")
        report["kannari"] = anchor_rows("kannari_device_rows")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=1, default=str))
    print(f"\nreport -> {out / 'report.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
