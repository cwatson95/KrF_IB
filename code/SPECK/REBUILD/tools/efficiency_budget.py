#!/usr/bin/env python3
"""Where does the efficiency go? — the KrF energy/channel budget of one run.

Owner ask (2026-08-21): "figure out a way to increase the efficiency or
see if it is just an ionization issue." This tool answers the second half
by decomposition: it runs a preset (0-D) and reconstructs, from the saved
history and the engine's own per-row fluxes (CompiledNetwork.row_fluxes,
gated to reproduce the RHS), the full chain

  deposited energy  -> beam partition (ion pairs / excitation / electron
                       heat / sub-threshold) from the Spencer-Fano yields
                    -> KrF*(B) FORMATION by channel class, integrated over
                       the run: ion-ion recombination (Kr2+/Kr+/... + F-),
                       neutral harpoon (Kr*/Kr** + F2), transfer from
                       ArF*/Ar2F*/Kr2F*/ArKrF*, other
                    -> KrF*(B) LOSS by channel class: stimulated emission
                       into the cavity (the laser), ASE, spontaneous B->X,
                       quench by partner (e / Ar / Kr / F2 / three-body),
                       other
                    -> the ledger's end sinks (laser / cavity loss /
                       radiated escape / stores)

and the derived figures: eV of deposited energy per KrF*(B) formed, the
formation efficiency (KrF* formed x h nu / E_dep), the extraction
efficiency (laser photons / KrF* formed), and the ION-PAIR RECOVERY — the
fraction of beam-made ion pairs that end up as a KrF*(B) through the
ion-ion channels. "Ionization issue" has a number: the energy the beam
sinks into ion pairs (W x pairs) vs what the ion channels return as KrF*.

Everything here is a RE-DESCRIPTION of the engine's solution (0-D: the
volume-averaged histories are the state), not new physics; predictions
scored against nothing.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/efficiency_budget.py
    python3 SPECK/REBUILD/tools/efficiency_budget.py --preset krf_hyper_mw \
        --set pressure_atm=1.0 --set beam_j_A_cm2=200
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

ANIONS = ("F-",)
ION_PARTNERS = ("Kr+", "Kr2+", "Ar+", "Ar2+", "ArKr+", "Ne+", "Ne2+", "F2+",
                "NeKr+")
HARPOON_DONORS = ("Kr*", "Kr**", "Kr2*")
TRANSFER_DONORS = ("ArF*", "Ar2F*", "Kr2F*", "ArKrF*", "NeKrF*", "NeF*")
QUENCH_PARTNER_ORDER = ("e", "F2", "Kr", "Ar", "Ne")


def classify_formation(r, upper):
    rs = set(r.reactants)
    if rs & set(ANIONS) and rs & set(ION_PARTNERS):
        return "ion-ion recombination"
    if rs & set(HARPOON_DONORS) and "F2" in rs:
        return "neutral harpoon (Kr*/Kr**/Kr2* + F2)"
    if rs & set(TRANSFER_DONORS):
        return "transfer from ArF*/Ar2F*/Kr2F*/ArKrF*"
    if "e" in rs and any(x in rs for x in ION_PARTNERS):
        return "electron recombination"
    return "other formation"


def classify_loss(r, upper):
    rs = list(r.reactants)
    partners = [x for x in rs if not x.startswith("KrF*")]
    if r.kind == "rad" or (not partners and len(rs) == 1):
        return "spontaneous B->X (+ other rad)"
    if "e" in partners:
        return "quench by electrons"
    heavy = [x for x in partners if x != "e"]
    if len(heavy) >= 2:
        return "three-body (Ar/Kr 2-body partners)"
    for pn in ("F2", "Kr", "Ar", "Ne"):
        if pn in heavy:
            return f"quench by {pn}"
    if any(x in heavy for x in ANIONS + ION_PARTNERS):
        return "loss to ions"
    return "other loss"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="krf_hyper_mw")
    ap.add_argument("--set", action="append", default=[],
                    help="param=value overrides (repeatable)")
    ap.add_argument("--n-save", type=int, default=401)
    ap.add_argument("--mixture", default="ar", choices=("ar", "ne"),
                    help="buffer: ar = Ar/Kr/F2 0.949/0.05/0.001 (preset), "
                         "ne = Ne/Kr/F2 0.949/0.05/0.001 (the Ne-fill "
                         "LOKI tables, krf_ne_0d's fill)")
    ap.add_argument("--f2", default=None,
                    help="F2 mole fraction override (Kr 5 %, buffer "
                         "renormalized) — e.g. 0.002")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    ov = {"mw_power_W": 0.0}
    for kv in args.set:
        k, v = kv.split("=", 1)
        ov[k] = float(v) if k not in ("mode", "drive", "electron_tier",
                                     "beam_subtable",
                                     "discharge_mode") else v
    if args.mixture == "ne" or args.f2 is not None:
        from REBUILD.MODULES.engine import Device
        from REBUILD.MODULES.params import engine_parameters
        from REBUILD.run import PRESETS
        buf = "Ne" if args.mixture == "ne" else "Ar"
        f2 = 0.001 if args.f2 is None else float(args.f2)
        mix = {buf: 1.0 - 0.05 - f2, "Kr": 0.05, "F2": f2}
        P = engine_parameters()
        P.set_many(PRESETS[args.preset]["params"], source="problem")
        P.set_many(ov, source="cli")
        dev = Device(P, mixture=mix)
        ov["mixture"] = f"{buf}{mix[buf]:.4f}_Kr0.05_F2{f2:g}"
    else:
        dev = build_device(args.preset, ov)
    res = dev.run(n_save=args.n_save)
    if not res["success"]:
        sys.exit("run failed")
    lay, net = dev.lay, dev.net
    upper = dev.optics["upper"]
    i_up = lay.index[upper]
    v = float(dev.cell_vols.sum())
    t_s = res["t_ns"] * 1e-9
    hv_J = dev.cav.photon_J
    hv_eV = dev.hv_eV
    st = dev.compiled.row_stoichiometry()
    by_rid = {r.rid: r for r in net}
    # the KrF* MANIFOLD (B, C, vibrational levels) is one reservoir for
    # this accounting: rows that only move population inside it (B<->C
    # collisional mixing, vibrational relaxation) are neither formation
    # nor loss — they dominate the raw per-state flux and would swamp the
    # real channels (found on the first run: CB_Ar/BC_Ar at 1e20 each way)
    manifold = {i for s, i in lay.index.items() if s.startswith("KrF*")}
    net_manifold = {rid: sum(c for i_s, c in pairs if i_s in manifold)
                    for rid, pairs in st.items()}

    # --- integrate per-row fluxes over the saved history ------------------
    n_t = len(t_s)
    form = {}
    loss = {}
    ste_rate = np.zeros(n_t)
    se_rate = np.zeros(n_t)
    ase_rate_t = np.zeros(n_t)
    for j in range(n_t):
        dens = np.array([res["dens_avg"][s][j] for s in lay.species])
        dens = np.maximum(dens, 0.0).reshape(lay.n_species, 1)
        ne = np.atleast_1d(lay.electron_density(dens[:, 0]))
        te = np.atleast_1d(res["Te_eV"][j])
        tg = np.atleast_1d(res["Tg_K"][j])
        fl = dev.compiled.row_fluxes(dens, te, tg, ne, dev.n_tot)
        for rid, flux in fl.items():
            coeff = net_manifold.get(rid, 0.0)
            if coeff == 0.0:
                continue
            f = float(flux[0])
            (form if coeff > 0 else loss).setdefault(rid, np.zeros(n_t))[j] = \
                abs(coeff) * f
        n_ph = float(res["n_ph"][j])
        from REBUILD.MODULES.units import C_LIGHT
        ste_rate[j] = (n_ph * C_LIGHT * dev.optics["sigma_se_m2"]
                       * float(dens[i_up, 0]))
        ase_rate_t[j] = ((dev._ase_rate(dens) if dev._ase_on else 0.0)
                         * float(dens[i_up, 0]))
    trap = getattr(np, "trapezoid", None) or np.trapz

    def integ(arr):
        return float(trap(arr, t_s)) * v            # molecules (or events)

    form_by_class = {}
    form_rows = {}
    for rid, arr in form.items():
        c = classify_formation(by_rid[rid], upper)
        n = integ(arr)
        form_by_class[c] = form_by_class.get(c, 0.0) + n
        form_rows[rid] = n
    loss_by_class = {}
    loss_rows = {}
    for rid, arr in loss.items():
        c = classify_loss(by_rid[rid], upper)
        n = integ(arr)
        loss_by_class[c] = loss_by_class.get(c, 0.0) + n
        loss_rows[rid] = n
    n_ste = integ(ste_rate)
    n_ase = integ(ase_rate_t)
    # --- where the stimulated photons go: extracted / mirror+window loss /
    # photo-absorbed by species (the cavity absorber inventory) -----------
    from REBUILD.MODULES.units import C_LIGHT as _C
    absorbed = {}
    for sp, sig in dev.cav.absorbers.items():
        if sp not in res["dens_avg"]:
            continue
        r_a = np.maximum(res["n_ph"], 0.0) * _C * sig * res["dens_avg"][sp]
        absorbed[sp] = integ(r_a)
    n_abs = sum(absorbed.values())
    loss_by_class["stimulated emission (cavity)"] = n_ste
    loss_by_class["ASE"] = n_ase
    n_formed = sum(form_by_class.values())
    n_lost = sum(loss_by_class.values())

    # --- beam partition from the yields --------------------------------
    y = dev.beam_yields
    e_dep = res["final_channels_J"]["drive_in"]
    pairs = e_dep / (y.W_eV * 1.602176634e-19)
    from REBUILD.MODULES.species import charge_of, e_form_of
    ion_eV = sum(n * e_form_of(s) for s, n in y.per_pair.items()
                 if charge_of(s) > 0)
    exc_eV = sum(n * e_form_of(s) for s, n in y.per_pair.items()
                 if charge_of(s) == 0)
    part = {"W_eV": y.W_eV, "ion_pairs": pairs,
            "ionization_eV_per_pair": ion_eV,
            "excitation_eV_per_pair": exc_eV,
            "eps_sec_eV_per_pair": y.eps_sec_eV,
            "subfloor_eV_per_pair": y.subfloor_eV}
    fin = res["final_channels_J"]
    laser_photons = fin["laser_out"] / hv_J
    n_ion_pairs = pairs
    recov = form_by_class.get("ion-ion recombination", 0.0) / max(n_ion_pairs, 1)
    report = {
        "preset": args.preset, "overrides": ov, "regime": res["regime"],
        "E_dep_J": e_dep, "E_laser_J": fin["laser_out"],
        "eta_intrinsic": res["eta_intrinsic"],
        "beam_partition": part,
        "KrF_B_formed": n_formed, "KrF_B_lost_accounted": n_lost,
        "formation_by_class": form_by_class,
        "loss_by_class": loss_by_class,
        "top_formation_rows": dict(sorted(form_rows.items(),
                                          key=lambda kv: -kv[1])[:8]),
        "top_loss_rows": dict(sorted(loss_rows.items(),
                                     key=lambda kv: -kv[1])[:8]),
        "eV_deposited_per_KrF_formed": (e_dep / 1.602176634e-19) / max(n_formed, 1),
        "formation_efficiency": n_formed * hv_J / e_dep,
        "extraction_efficiency": laser_photons / max(n_formed, 1),
        "ion_pair_recovery_as_KrF": recov,
        "ledger_sinks_J": {k: fin[k] for k in ("laser_out", "cavity_loss",
                                               "radiated_escape")},
        "ledger_stores_J": res["closure"][-1]["stores_delta_J"],
    }
    # --- print ----------------------------------------------------------
    print(res["regime"])
    print(f"\nDEPOSITED {e_dep:.2f} J -> LASER {fin['laser_out']:.4f} J "
          f"(intrinsic {100 * res['eta_intrinsic']:.2f}%)")
    print(f"beam partition per ion pair (W = {y.W_eV:.1f} eV): ionization "
          f"{ion_eV:.1f} eV ({100 * ion_eV / y.W_eV:.0f}%), excitation "
          f"{exc_eV:.1f} eV ({100 * exc_eV / y.W_eV:.0f}%), electron heat "
          f"{y.eps_sec_eV:.1f} eV ({100 * y.eps_sec_eV / y.W_eV:.0f}%), "
          f"sub-threshold {y.subfloor_eV:.2f} eV; ion pairs {pairs:.3e}")
    print(f"\nKrF*(B) FORMED {n_formed:.3e} ({(e_dep / 1.602176634e-19) / max(n_formed, 1):.0f} eV deposited per KrF*; "
          f"formation efficiency {100 * n_formed * hv_J / e_dep:.1f}% of E_dep as KrF* quanta)")
    for c, n in sorted(form_by_class.items(), key=lambda kv: -kv[1]):
        print(f"   {100 * n / max(n_formed, 1):5.1f}%  {c}")
    print(f"   ion-pair recovery as KrF* via ion-ion channels: "
          f"{100 * recov:.1f}% of {pairs:.2e} pairs")
    print(f"\nKrF*(B) LOST {n_lost:.3e} (accounted; formed-lost mismatch "
          f"{100 * (n_formed - n_lost) / max(n_formed, 1):+.2f}% = stored/unresolved)")
    for c, n in sorted(loss_by_class.items(), key=lambda kv: -kv[1]):
        print(f"   {100 * n / max(n_lost, 1):5.1f}%  {c}")
    print(f"   extraction: {laser_photons:.3e} laser photons = "
          f"{100 * laser_photons / max(n_formed, 1):.1f}% of KrF* formed "
          f"({100 * n_ste / max(n_formed, 1):.1f}% stimulated incl. "
          f"intracavity loss)")
    cav_loss_photons = fin["cavity_loss"] / hv_J
    print(f"   stimulated photons {n_ste:.3e}: extracted "
          f"{100 * laser_photons / max(n_ste, 1):.1f}%, HR/window loss "
          f"{100 * cav_loss_photons / max(n_ste, 1):.1f}%, photo-absorbed "
          f"{100 * n_abs / max(n_ste, 1):.1f}% — by species: "
          + ", ".join(f"{sp} {100 * n / max(n_abs, 1):.0f}%"
                      for sp, n in sorted(absorbed.items(),
                                          key=lambda kv: -kv[1])[:6]))
    report["stimulated_photons"] = n_ste
    report["photoabsorbed_by_species"] = absorbed
    print("\nledger end sinks: " + ", ".join(
        f"{k} {fin[k]:.3f} J" for k in ("laser_out", "cavity_loss",
                                        "radiated_escape")))
    print("ledger stores (delta): " + ", ".join(
        f"{k} {v_:+.3f} J" for k, v_ in res["closure"][-1]["stores_delta_J"].items()))
    print("\ntop formation rows:", {k: f"{v_:.2e}" for k, v_ in report["top_formation_rows"].items()})
    print("top loss rows:", {k: f"{v_:.2e}" for k, v_ in report["top_loss_rows"].items()})
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "efficiency")
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = "_".join(f"{k}{v_}" for k, v_ in ov.items() if k != "mw_power_W") or "ref"
    path = out_dir / f"budget_{args.preset}_{tag}.json"
    path.write_text(json.dumps(report, indent=2, default=float))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
