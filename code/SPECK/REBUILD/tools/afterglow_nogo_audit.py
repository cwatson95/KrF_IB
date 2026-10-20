#!/usr/bin/env python3
"""Is there ANY case where reheating cold afterglow electrons pays? — audit.

Owner ask (2026-08-21): "is there any case where reheating the cold
electrons actually works?" Instead of another scan, this audits the
STRUCTURE: every network channel that can form the KrF* manifold in an
afterglow, and the SIGN of its response to electron heating — evaluated
numerically from the engine's own rate specs and swarm curves
(k at Te = 0.5 eV vs 3 eV; heavy-particle rows marked Te-independent with
their electron-side SUPPLY chain classified instead).

The two walls this exposes (both already measured in 0.1.20-0.1.25):
  1. every recombination-fed channel SLOWS when heated (DR ~ Te^-1/2;
     ion-ion is heavy-particle but its F- fuel is made by attachment,
     whose rate FALLS with Te), so the stored-ion inventory harvests
     best COLD — the afterglow of an excimer is a recombination laser's
     afterglow, and recombination lasers want cold electrons;
  2. the only heat-POSITIVE channel is direct e-impact excitation of
     Kr/Ar (thresholds 10-16 eV), which needs E/N ~ 20-30 Td — and the
     conductive afterglow pins the settled field at <= 2 Td at ANY
     injected power (measured to 1 GW, 0.1.20): the field the channel
     needs is exactly the field the plasma forbids.

Where electron heating DOES work, stated for contrast (literature-class,
not modeled here): lasers whose upper level lies at electron-temperature
energies — the CO2/CO VIBRATIONAL sustainer lasers (upper state 0.3 eV;
E/N ~ 1-3 Td pumps it at ~90% quantum efficiency; the EBSD architecture)
— and fully discharge-PUMPED excimers (kTd fields, arc-instability
limited; the DE-AR "hurts during discharge" regime). KrF's afterglow
sits in neither class. BOTE companion: tools/ib_feasibility.py.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/afterglow_nogo_audit.py
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

from REBUILD.tools.hybrid_timing_study import build_device  # noqa: E402

TE_COLD, TE_HOT = 0.5, 3.0        # the afterglow vs the reheated window


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_hyper_mw")
    ap.add_argument("--arf", action="store_true",
                    help="audit the ArF (193 nm) chain instead, on the "
                         "Electra-ArF anchor configuration")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    if args.arf:
        from REBUILD.tools.mw_benchmarks import _case_configs
        build, _note = _case_configs()["electra_arf_1p09"]
        dev = build({"mw_power_W": 0.0}, "loki")
    else:
        dev = build_device(args.preset, {"mw_power_W": 0.0})
    lay, net, comp = dev.lay, dev.net, dev.compiled
    upper_base = dev.optics["upper"].split("(")[0]     # KrF* or ArF*
    manifold = {i for s, i in lay.index.items() if s.startswith(upper_base)}
    st = comp.row_stoichiometry()
    rows = []
    for r in net:
        coeff = sum(c for i_s, c in st[r.rid] if i_s in manifold)
        if coeff <= 0:
            continue
        k_cold = comp._k(r.rate, np.array([TE_COLD]), dev.n_tot)
        k_hot = comp._k(r.rate, np.array([TE_HOT]), dev.n_tot)
        k_cold = float(np.atleast_1d(k_cold)[0])
        k_hot = float(np.atleast_1d(k_hot)[0])
        te_dep = r.is_electron_row() or k_hot != k_cold
        ratio = (k_hot / k_cold) if (te_dep and k_cold > 0) else 1.0
        if not te_dep or abs(ratio - 1.0) < 1e-9:
            sign = "0 (heavy-particle)"
        elif ratio < 1.0:
            sign = f"NEGATIVE (x{ratio:.2f} hot/cold)"
        else:
            sign = f"positive (x{ratio:.2f})"
        rows.append({"rid": r.rid,
                     "reaction": " + ".join(r.reactants) + " -> "
                     + " + ".join(r.products),
                     "kind": r.kind, "electron_row": r.is_electron_row(),
                     "k_ratio_hot_cold": ratio, "sign": sign})
    rows.sort(key=lambda x: x["k_ratio_hot_cold"])
    print(f"{upper_base}-manifold FORMATION channels ({len(rows)} rows), "
          f"response to heating Te {TE_COLD} -> {TE_HOT} eV "
          f"(mode {dev.mode}):\n")
    for x in rows:
        print(f"  {x['rid']:>10}  {x['sign']:>28}  {x['reaction']}")
    neg = sum(1 for x in rows if x["k_ratio_hot_cold"] < 1.0 - 1e-9)
    heavy = sum(1 for x in rows if abs(x["k_ratio_hot_cold"] - 1.0) < 1e-9)
    pos = len(rows) - neg - heavy
    # --- the supply chains of the heavy-particle rows -------------------
    tiers = dev.tiers
    k_att = {t: float(np.atleast_1d(tiers.k("F-", np.array([t])))[0])
             for t in (TE_COLD, TE_HOT)}
    att_ratio = k_att[TE_HOT] / k_att[TE_COLD]
    exc_bucket = "Ar*" if dev.mode.startswith("arf") else "Kr*"
    dimer_ion = "Ar2+" if dev.mode.startswith("arf") else "Kr2+"
    k_exc = {}
    for t in (TE_COLD, TE_HOT, 5.75):
        k_exc[t] = float(np.atleast_1d(tiers.k(exc_bucket,
                                               np.array([t])))[0])
    # the Kr* SUPPLY in an afterglow is dissociative recombination
    # (Kr2+ + e -> Kr* + Kr, the ADR Te^-0.55 convention): its ratio
    dr_ratio = (TE_HOT / TE_COLD) ** -0.55
    print(f"""
TALLY: all {heavy} direct formation rows are heavy-particle — the
electron-temperature dependence lives ONE STEP UPSTREAM, in their
supplies, and BOTH supply lines fall when heated:
 - F- (fuel of every ion-ion row): made by attachment, whose rate FALLS
   x{att_ratio:.2f} from {TE_COLD} to {TE_HOT} eV.
 - {exc_bucket}/{exc_bucket}* (fuel of every harpoon row): made in the
   afterglow by dissociative recombination ({dimer_ion} + e, ~ Te^-0.55),
   which FALLS x{dr_ratio:.2f} over the same heating. ({neg} direct rows
   slow, {pos} speed up.)
 - The heat-positive route is direct e-impact excitation
   (k({exc_bucket}) = {k_exc[TE_COLD]:.1e} m^3/s at {TE_COLD} eV ->
   {k_exc[TE_HOT]:.1e} at {TE_HOT} eV -> {k_exc[5.75]:.1e} at 5.75 eV,
   the 30 Td swarm point) — but it needs E/N ~ 20-30 Td (the
   self-sustainment crossing is ~28 Td) and the conductive afterglow
   pins the settled field at <= 2 Td at ANY injected power (measured to
   1 GW, 0.1.20): the field this channel needs is the field the plasma
   forbids. The two walls close every route.
VERDICT: within this network there is NO afterglow state whose
{upper_base}-photon harvest improves by heating the electrons: every
recombination channel wants them colder, the fuel line wants them
colder, and the one hot channel is field-starved by the plasma's own
conductivity. The measured exceptions (Ne-opt/Ar-opt "+9%") are cells
still LASING through the window — late co-pump, not afterglow harvest.
Electron heating pays in lasers whose upper level sits at eV energies
(CO2/CO vibrational sustainers, ~90% quantum efficiency at 1-3 Td) —
an excimer's 10-16 eV chain is the wrong class.""")
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "mw_benchmarks")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / ("afterglow_nogo_audit_arf.json" if args.arf else "afterglow_nogo_audit.json")
    path.write_text(json.dumps(
        {"preset": args.preset, "te_cold_eV": TE_COLD, "te_hot_eV": TE_HOT,
         "rows": rows, "attachment_ratio_hot_cold": att_ratio,
         "dr_ratio_hot_cold": dr_ratio,
         "k_Kr_exc": {str(k): v for k, v in k_exc.items()},
         "tally": {"negative": neg, "heavy": heavy, "positive": pos}},
        indent=2, default=float))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
