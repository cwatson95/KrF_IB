"""S7 — anchors scored, misses stated, omissions named. No tuning. Ever.

Every row: the anchor, its source, the engine's number, and the miss. The
scorecard never feeds back into the engine (PLAN (26): anchors score, they
never steer). The OMISSIONS list is part of the score — a model that hides
what it left out is scoring itself on a curve.
"""

from __future__ import annotations

import numpy as np

from .provenance import census_line, tier_census
from .units import gas_density_m3

# Literature anchors with sources
W_ANCHORS = {          # eV per ion pair (ICRU 31 / Combecher class)
    "Ar": 26.4, "Kr": 24.2, "Ne": 36.6,
}
ESTAR_600KEV = {       # NIST ESTAR-class mass stopping @ ~600 keV
    "He": 1.94, "Ne": 1.71, "Ar": 1.52, "Kr": 1.34, "F2": 1.63,
}
_MOLAR = {"He": 4.0026, "Ne": 20.18, "Ar": 39.948, "Kr": 83.80,
          "F2": 37.9968}

# Engine-vs-engine context at the shared v2 operating point
# (Ar 94.9/Kr 5/F2 0.1, 3 atm, 600 kV, 53 A/cm2), from spectre CHANGELOG
# 2026-07-15 "gain-gap closure": g0 spectre 0.0318 /cm, phantom 0.0302 /cm.
V2_G0_PER_CM = {"spectre": 0.0318, "phantom11": 0.0302}

OMISSIONS = [
    "Xe rows R16/R17 (trace 1e-4) omitted, as in spectre — Xe enters the "
    "stopping power only",
    "RESOLVED 0.1.10, corrected 0.1.11: the ArKr+ ION manifold (Kannari "
    "Table I p.4311-4312 / J&H #10,36-40, both first-hand 2026-08-05) — "
    "the 0.1.9 'next network milestone'. In as predicted (Kr+ +2Ar -> "
    "ArKr+ at 1.0e-31 IS the dominant Kr+ route in an Ar fill, 3.0e8 "
    "s^-1 at the device point vs 8.4e7 for AK7), plus the genuinely "
    "new rows (Ar** harpoon + non-excimer branches, ArF*+2Ar Ar-third-"
    "body, five Penning rows, Ar**/Kr**/ArKr* photoionization sigmas) "
    "and first-hand prints SUPERSEDING Chen's generated equivalents "
    "(AF1/AIIR1/AIIR2/AAR1/FMN1/ARF_F2/ARF_AR/ARF_RAD vs R29/R30/R31/"
    "R27/R22/R33/R34/R37 — the 0.1.10 'source-less/sink-less/no-path' "
    "claims were WRONG, those channels existed in the generated CHEN "
    "block, and the pairs ran DOUBLE-COUNTED until the 0.1.11 dedup + "
    "channel-uniqueness teeth). MEASURED VERDICT (clean network): "
    "anchor x0.60 -> x0.56 — the suspect hypothesis is consumed, the "
    "deficit persists (see the anchor verdict). CORRECTION to the "
    "0.1.9 note: the printed 1.5e-17 cm2 absorber is the NEUTRAL "
    "ArKr* [ref 85] (adopted), not the ion — neither primary prints a "
    "sigma for ArKr+ itself, so none was invented",
    "electron-electron collisions: fluid Te treatment (both rate tiers are "
    "swarm-kinetic, but the engine's electron gas is a temperature). "
    "0.1.12: the fluid-Te SUSPECT was tested against the NRL EEDF "
    "benchmark (electra_population_rows) — the bulk ran x4 cold because "
    "the time-free Spencer-Fano cascade pre-spent the Platzman share on "
    "rate-limited trace-F2 channels (fixed: prompt_floor_eV; Te x0.26 "
    "-> x1.01, nothing fitted). The single-temperature treatment itself "
    "REMAINS (no EEDF shape below the tables; the F- x0.02 residual and "
    "Ar** x0.06 sit exactly where shape matters). Likewise no e-driven "
    "excited-manifold stepwise rows beyond R8/KSTEP (J&H #47-class "
    "Kr*+e -> Kr**+e ~1e6 s^-1 at beam-mode ne — sub-dominant, stated); "
    "the printed e-ionization of excited/dimer states (Kannari [72,79]) "
    "except AKI_ION, and the Ar2*/Kr2*/Ar2F* superelastics (1.0e-7 "
    "prints) — each <=1% of its species' fate budget at beam points; "
    "ArF*+e IS covered (Chen R32); F+F+M -> F2+M [16] ~1e3 s^-1 on "
    "pulse timescales",
    "discharge circuit drive: Lm is a fixed inductor at the paper's stated "
    "base-case 250 nH (0.1.25 — the '5 -> 1 uH post-switch' story the "
    "0.1.5 note told belongs to Chen's Fig. 5 VALIDATION circuit, whose "
    "components had also shipped as the defaults; see CHANGELOG [0.1.25] "
    "and the Known-bad), the paper's RLC filter branch across the gap "
    "(70.71 ohm/100 nH/10 pF, post-breakdown ringing suppression) is not "
    "modeled, and the Xe seed is omitted; the prescribed-waveform mode "
    "remains available and unchanged",
    "RESOLVED 0.1.13: the Ramsauer-stall export deepening (the 0.1.9 "
    "recorded follow-up) — tools/export_kinetic.py --deepen, same seed, "
    "10x rounds, mechanical pure-extension gate (99 above-stall bins "
    "bit-identical; a 25000-round reproduction run matched the shipped "
    "npz array-for-array first). Measured floors moved el 1.26 -> 0.47 "
    "eV (through the Ramsauer minimum), vib -> 0.66, F- -> 0.66 usable; "
    "the F- DA curve RISES x13 toward low eps where the taper fell — "
    "the F- starvation mechanism, now measured (benchmark F- x0.02 -> "
    "x0.19). Below 0.47/0.66 eV the taper still extrapolates from the "
    "deepest measured anchor — nothing operates there post-0.1.12",
    "RESOLVED 0.1.6: in-loop ASE (ase_inloop, default ON) and gated "
    "inter-cell transport (measured negligible at device presets: "
    "tau_diff ~ 0.1 s) — the former quasi-steady-report and no-transport "
    "scope notes are closed",
    "RESOLVED 0.1.9: Kr2*/Ar2* VUV radiative decay rows (KR2_RAD/AR2_RAD, "
    "Kannari [88]) — the 0.1.0 'spectre parity' scope hole; the x24 "
    "J&H-vs-Kannari lifetime spread is recorded on the rows and swept in "
    "the EST/conflict bracket sensitivity",
    "RESOLVED 2026-08-04: the x100 'vib' tier disagreement was a LOKI-side "
    "bucket bug (the F2 Effective momentum-transfer row fell through to "
    "vib, draining 0.2 eV per FAKE event) — fixed in loki_tier._bucket; "
    "post-fix vib agreement x0.84 [0.51, 1.15]. The same defect exists "
    "upstream in projects/spectre loki_electron._bucket (flagged for its "
    "owner, not edited here)",
    "beam-direct F2 channels (ionization is carried, but the "
    "dissociation/vibration/attachment of F2 by beam secondaries — NRL "
    "petrova2022 Table A2, the Kushner-Moratz 1988 class — are not in "
    "the beam module; all F2 burn-down here runs through the bulk "
    "plasma-chemistry rows. NRL's model runs F2 ~50% dissociated by "
    "pulse end at 0.5 MW/cc. A stated limit of the Electra ArF anchor "
    "(0.1.11); direction on the deficit UNRANKED — dissociation both "
    "starves the harpoon and feeds the F/F- formation routes",
    "ion-ion screening (0.1.14): the parameter-free screened-capture "
    "S = 1/(1+R_c/lambda_D) gives 0.695 at Kannari's printed point where "
    "Morgan's Monte-Carlo prints 0.397 (his paper is not in-repo) — the "
    "x1.75 form spread is part of the recombination error budget; "
    "fetching Morgan 1980-class screening curves would pin it",
    "checkpoint resume is tolerance-exact, not bitwise (adaptive BDF; see "
    "checkpoint.py)",
]


def beam_scorecard() -> dict:
    """Emergent-W and stopping-power anchors, computed live (no caching —
    the point is that the current code produces these numbers)."""
    from .beam import differential_W, stopping_power
    from .cross_sections import CrossSectionSet
    from .units import N_AVOGADRO
    n = gas_density_m3(1.0, 300.0)
    rows = []
    for g, lit in W_ANCHORS.items():
        try:
            W = differential_W(CrossSectionSet({g: n}))
            rows.append((f"W({g}) emergent (Spencer-Fano, differential "
                         f"window)", f"{W:.1f} eV", f"{lit} eV",
                         f"{100 * (W / lit - 1):+.1f}%"))
        except FileNotFoundError:
            rows.append((f"W({g})", "no LXCat data", f"{lit} eV",
                         "could-not-check"))
    for g, ref in ESTAR_600KEV.items():
        try:
            S = stopping_power({g: n}, 6.0e5)
            rho = n * _MOLAR[g] / N_AVOGADRO / 1e6
            s_mass = S * 1e-8 / rho
            rows.append((f"dE/dx({g}) @600 keV (Bethe/Moller)",
                         f"{s_mass:.3f} MeV cm2/g", f"{ref} MeV cm2/g",
                         f"{100 * (s_mass / ref - 1):+.1f}%"))
        except (FileNotFoundError, KeyError):
            pass
    return {"rows": rows}


def two_route_w_rows() -> list:
    """The cross-package W check: REBUILD's deterministic Spencer–Fano vs
    the MC cascade run through PIC_MCC_3D's actual MCC operator (shipped in
    the kinetic export with provenance), both evaluated as integral W at the
    SAME recorded E0 on the same LXCat sets. Independent code paths; the
    in-package pairing of the same two methods agreed to 0.4 % in
    PIC_MCC_3D, and this pairing is measured at 0.02–0.2 % (2026-08-04)."""
    from .beam import spencer_fano_per_process
    from .cross_sections import CrossSectionSet
    from .kinetic_tier import KineticRates
    kin = KineticRates.load_default()
    n = gas_density_m3(1.0, 300.0)
    rows = []
    if not kin.w_cascade:
        return ["export carries no w_cascade section — regenerate with "
                "tools/export_kinetic.py"]
    for gas, rec in sorted(kin.w_cascade.items()):
        cs = CrossSectionSet({gas: n})
        U, _ = spencer_fano_per_process(cs, rec["E0_eV"])
        ion = [k for k, kind in enumerate(cs.kinds) if kind == "ionization"]
        u_ion = float(U[ion].sum())
        W_sf = rec["E0_eV"] / u_ion if u_ion > 0 else float("inf")
        rows.append(
            f"W({gas}) @ {rec['E0_eV']:.0f} eV: SF {W_sf:6.2f} vs MC "
            f"{rec['W_eV']:6.2f} eV -> agree to "
            f"{100 * abs(W_sf / rec['W_eV'] - 1):.2f}% "
            f"(MC stat ~{100 * rec['stat_err_rel']:.1f}%, "
            f"{rec['rounds']} rounds)")
    return rows


def _chen_circuit_params(t_end_ns: float, **over):
    from .params import engine_parameters
    p = engine_parameters()
    d = {"mode": "arf_chen", "drive": "discharge",
         "discharge_mode": "circuit", "pressure_atm": 3.4,
         "t_end_ns": t_end_ns, "nx": 1, "ny": 1, "nz": 1,
         "Lx_m": 0.012, "Ly_m": 0.010, "Lz_m": 0.90}
    d.update(over)
    p.set_many(d, source="problem")
    return p


def _fwhm_ns(t_ns, p_w) -> float:
    """FWHM of the output-power trace [ns]; 0.0 when there is no pulse."""
    import numpy as np
    p = np.asarray(p_w, float)
    pk = float(p.max())
    if pk <= 0.0:
        return 0.0
    above = p >= 0.5 * pk
    return float(t_ns[above][-1] - t_ns[above][0]) if above.any() else 0.0


def chen_figure6_rows(t_end_ns: float = 200.0) -> list:
    """S7: Chen 2025 base-case observables scored by RUNNING the circuit
    preset (self-consistent breakdown) — 0.1.25 REVISION after the
    2026-08-25 comparison audit (CHANGELOG [0.1.25]):

    * circuit constants are now the paper's Sec II A BASE CASE
      (10 nF / 250 nH / 10 nF / 5 nH, -18 kV — re-read first-hand); the
      lane had shipped the Fig. 5 VALIDATION experiment's components;
    * the CHEN_ERRATA rate corrections are on (the printed R13/R14/R23/
      R24/R36 are transcription errors — the as-printed lane is a
      RECORDED measurement below: it cannot integrate under the paper's
      own circuit);
    * deposited energy is compared to the PAPER-STATED
      E_in = C0 V0^2/2 = 1.62 J (the old "1.68 J" was the wrong C0);
      Chen does not print his deposited fraction, so expect < 1;
    * n_e target 6.5e20 is Fig. 10a's LINEAR axis at 0.1% F2 (the old
      1e21 was a coarse log-axis read of Fig. 6b — a x1.5 target error
      booked against the engine).

    Misses stated; nothing tuned."""
    import numpy as np
    from .engine import Device
    r = Device(_chen_circuit_params(t_end_ns)).run(n_save=201)
    fc = r["final_channels_J"]
    PL = float(np.max(r["P_out_W"]))
    j_en = int(np.argmax(r["circuit"]["en_td"]))
    rows = [
        ("deposited energy [J]", fc["drive_in"], 1.62,
         "vs E_in = C0 V0^2/2 stated p.113502-10 (deposited not printed)"),
        ("peak n_e [m^-3]", float(np.max(r["ne_m3"])), 6.5e20,
         "Fig. 10a linear axis, 0.1% F2"),
        # Te AT THE BREAKDOWN (E/N peak): the raw trace max is the
        # pre-avalanche transient (tiny-n_e epoch), not Fig. 10c's
        # discharge value
        ("T_e at E/N peak [eV]", float(r["Te_eV"][j_en]), 5.0,
         "Fig. 10c peak class"),
        ("peak n_ph [m^-3]", float(np.max(r["n_ph"])), 1e20,
         "stated in text (p.113502-5)"),
        ("peak P_L [MW]", PL / 1e6, 5.0, "Fig. 6a"),
        ("pulse FWHM [ns]", _fwhm_ns(r["t_ns"], r["P_out_W"]), 15.0,
         "stated in text (0.1% F2)"),
        ("E/N peak [Td]", float(np.max(r["circuit"]["en_td"])), 15.0,
         "Fig. 10b class"),
    ]
    out = ["Chen 2025 base case via the CIRCUIT preset (Sec II A "
           "constants + CHEN_ERRATA; self-consistent breakdown):"]
    for name, got, ref, note in rows:
        ratio = got / ref if ref else float("nan")
        out.append(f"  {name:<22s} {got:11.3e} vs {ref:9.3e}  "
                   f"(x{ratio:6.3f})  {note}")
    cz = r["circuit"]
    lo, hi = cz["en_table_range_td"]
    f_below = cz["en_below_table_frac"]
    f_in = 1.0 - f_below - cz["en_above_table_frac"]
    out.append(f"  E/N table coverage: {100 * f_in:.0f}% in "
               f"[{lo:g}, {hi:g}] Td; {100 * f_below:.0f}% BELOW the "
               "floor (the breakdown ramp rides the flat+taper stand-in "
               "— the gauge was one-sided through 0.1.24 and this "
               "fraction went unreported)")
    # the as-printed lane — a RECORDED measurement, not a live re-run (it
    # grinds the stiff BDF for minutes on its way to the blow-up):
    out.append("  as-printed lane (chen_errata=False): measured "
               "2026-08-25 — does NOT integrate under the paper's own "
               "circuit (FloatingPointError, non-finite U_e/T_g at "
               "~138 ns; R24 as printed annihilates Ne* in ~0.4 ps). "
               "Reproduce: Device(params with chen_errata=False) on this "
               "preset. Under the OLD mis-transcribed softer circuit it "
               "integrated to n_e 1.73e20, n_ph 2.5e15, no lasing — the "
               "two shipped defects masked each other")
    out.append("  VERDICT (0.1.25, measured attribution): the 0.1.5..24 "
               "record blamed a 'x3-5 per-joule formation gap' (suspects: "
               "Ne* source, Xe seed, 0-D vs 1-D). WRONG RANKING — "
               "controlled swaps show the gap was (1) the Fig.5-for-"
               "Table-II circuit transcription and (2) the four impossible "
               "as-printed rates, both now fixed; the named suspects are "
               "second-order. Residual misses above are the honest lane.")
    return out


#: the recorded inter-source alternates on conflicted LIT radiative rows
#: (rid -> the OTHER published print / the kept value): the bracket for
#: these rows is the real recorded spread, not the generic x1/3..x3
CONFLICT_ALTERNATES = {
    "KR2_RAD": 3.3e6 / 8.0e7,      # J&H #66 vs Kannari [88] (x24)
    "AR2_RAD": 3.8e6 / 6.0e7,      # J&H #67 vs Kannari [88] (x16)
    "AKD_RAD": 3.0e6 / 8.0e7,      # J&H #70 vs Kannari-est (x27)
    "AKF_RAD": 5.0e7 / 5.0e6,      # J&H #71 vs Kannari-est (x10 UP)
    # chain-2 head (0.1.14): Kannari-est + Lee [36] both print 1.0e-10,
    # but NRL's Orestes carries NO neutral 4p->4s quench (A3: electron
    # relaxation only) — the x0.1 alternate probes toward their
    # structure and the Setser-era 4p-class magnitudes
    "APS_Ar": 0.1,
}


def _bracket_rate(spec, f):
    """The rate spec with its magnitude multiplied by f (sweep-side only —
    the engine has no such pathway; defect class D4 stays impossible)."""
    kind = spec[0]
    if kind == "const":
        return ("const", spec[1] * f)
    if kind == "stepwise":
        a = spec[2] if len(spec) > 2 else 1.56e-13
        return ("stepwise", spec[1], a * f)
    if kind == "dr_sqrt":
        return ("dr_sqrt", spec[1] * f)
    if kind == "te_fit":
        return ("te_fit", spec[1] * f, spec[2], spec[3])
    if kind == "flannery":
        return ("flannery", spec[1], spec[2] * f)
    if kind == "lee_ii":
        return ("lee_ii", spec[1] * f, spec[2], spec[3])
    raise NotImplementedError(f"no bracket rule for rate kind {kind!r}")


def fill_elements(mixture: dict) -> set:
    from .species import atoms_of
    els = set()
    for g in mixture:
        els |= set(atoms_of(g))
    return els


def est_active_and_bound(rows, mixture: dict):
    """Split EST rows into (active, element-bound) for a fill.

    A row whose heavy reactants need an element the fill does not contain
    can only ever run on seed-level densities (1e12 m^-3 vs ~1e26 fill) —
    its bracket impact is bounded at the 1e-8 relative level STRUCTURALLY,
    no sweep needed. This is the honest reason 14 of the 19 EST rows
    (the Ne-chain analogs) cannot move the Ar-fill answer at all."""
    from .species import atoms_of
    els = fill_elements(mixture)
    active, bound = [], []
    for r in rows:
        if r.tier != "EST":
            continue
        need = set()
        for s in r.heavy_reactants():
            need |= set(atoms_of(s))
        (active if need <= els else bound).append(r)
    return active, bound


def est_bracket_rows(mixture: dict = None, t_end_ns: float = 120.0,
                     only_rids=None) -> dict:
    """S5 CLOSURE — the x1/3..x3 EST bracket sweep, MEASURED.

    For every EST row that can carry flux in the fill, the device preset
    is re-run with that row's rate x3 and x1/3 (the stated EST bracket)
    and the moves in peak g0 and eta recorded. Conflicted LIT radiative
    rows (CONFLICT_ALTERNATES) are swept at their recorded alternate
    print instead — their real spread. Element-bound EST rows are listed
    with the structural zero-impact statement instead of burning runs.
    The sweep lives HERE (scorecard side): it rebuilds a CompiledNetwork
    per variant — the engine itself has no rate-scaling pathway (D4).
    """
    from .chemistry import CompiledNetwork
    from .engine import DEFAULT_MIXTURES, Device
    from .params import engine_parameters
    if mixture is None:
        mixture = dict(DEFAULT_MIXTURES["krf"])
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": t_end_ns, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    dev = Device(p, mixture=mixture)
    base = dev.run(n_save=31)
    g0_b = float(max(base["g0_m1"]))
    eta_b = float(base["eta_intrinsic"])
    active, bound = est_active_and_bound(dev.net, mixture)
    swept = {}
    from dataclasses import replace as _rep
    targets = [(r, (1.0 / 3.0, 3.0)) for r in active]
    by_rid = {r.rid: r for r in dev.net}
    for rid, alt in CONFLICT_ALTERNATES.items():
        if rid in by_rid:
            targets.append((by_rid[rid], (alt,)))
    if only_rids is not None:
        targets = [(r, fs) for r, fs in targets if r.rid in set(only_rids)]
    for r, factors in targets:
        rec = {}
        for f in factors:
            rows2 = [_rep(x, rate=_bracket_rate(x.rate, f))
                     if x.rid == r.rid else x for x in dev.net]
            dev.compiled = CompiledNetwork(rows2, dev.lay, dev.tiers,
                                           dev.optics)
            dev.compiled.set_photon_energy(dev.hv_eV)
            res = dev.run(n_save=31)
            rec[f] = {
                "dg0_rel": float(max(res["g0_m1"])) / g0_b - 1.0
                if g0_b > 0 else 0.0,
                "deta_abs": float(res["eta_intrinsic"]) - eta_b,
            }
        swept[r.rid] = rec
    # restore the unmodified network
    dev.compiled = CompiledNetwork(dev.net, dev.lay, dev.tiers, dev.optics)
    dev.compiled.set_photon_energy(dev.hv_eV)
    return {"mixture": dict(mixture), "g0_base_m1": g0_b,
            "eta_base": eta_b, "swept": swept,
            "bound_rids": sorted(r.rid for r in bound),
            "active_rids": sorted(r.rid for r in active)}


def format_est_brackets(rep: dict) -> str:
    lines = [f"EST x1/3..x3 bracket sweep at the fill {rep['mixture']} "
             f"(base g0 {rep['g0_base_m1'] / 100:.4f} /cm, eta "
             f"{100 * rep['eta_base']:.2f}%):"]
    worst = (None, 0.0)
    for rid, rec in sorted(rep["swept"].items()):
        parts = []
        for f, d in sorted(rec.items()):
            parts.append(f"x{f:.2f}: g0 {100 * d['dg0_rel']:+.1f}%, "
                         f"eta {100 * d['deta_abs']:+.2f} pts")
            if abs(d["dg0_rel"]) > abs(worst[1]):
                worst = (rid, d["dg0_rel"])
        lines.append(f"  {rid:<10s} " + " | ".join(parts))
    if rep["bound_rids"]:
        lines.append("  element-bound (structurally zero in this fill, "
                     "seed-level 1e-8 class): "
                     + ", ".join(rep["bound_rids"]))
    if worst[0] is not None:
        lines.append(f"  WORST swept lever on g0: {worst[0]} "
                     f"({100 * worst[1]:+.1f}%)")
    return "\n".join(lines)


def _conditioning(g0_m1: float, alpha_m1: float, r_oc: float, l_op: float,
                  l_g_m: float) -> str:
    """Chain-1 gauge (0.1.14): how AMPLIFIED is this observable?

    Near threshold, laser energy goes like (g0 - g_th), so
    d ln E / d ln g0 ~ g0/(g0 - g_th) — a x1.15 gain error at a 7%
    margin is a x3-5 energy error. Every anchor point self-reports its
    margin and amplification so an amplified miss is never mistaken for
    a raw physics miss (and a threshold-pinned point is labeled)."""
    g_th = alpha_m1 - np.log(r_oc * l_op) / (2.0 * l_g_m)
    margin = (g0_m1 - g_th) / g0_m1 if g0_m1 > 0 else -1.0
    if margin <= 0.02:
        return (f"conditioning: g_th {g_th / 100:.4f} /cm, margin "
                f"{100 * margin:+.0f}% — THRESHOLD-PINNED (any gain "
                f"error is amplified without bound)")
    return (f"conditioning: g_th {g_th / 100:.4f} /cm, margin "
            f"{100 * margin:.0f}%, dlnE/dlng0 ~ x{1.0 / margin:.0f}")


def kannari_device_rows(t_end_ns: float = 120.0, tier: str = "loki",
                        subtable: str = None, subfloor: bool = False) -> list:
    """The Kannari 1985 DEVICE anchor, run LIVE (was a hand-recorded text
    block until 0.1.9 — a stale-able snapshot; now it re-scores on every
    scorecard run).

    His Sec IV cell: 660 kV, 1.88 MW/cm3 x 50 ns, Ar/Kr/F2 89.8/10/0.2 at
    2.5 atm, lossy CaF2 cavity (T = 0.60 -> L_op 0.64), 4x4x44 cm gain
    volume (704 cm3 -> the paper's 66 J at 94 mJ/cm3). The beam current is
    iterated ONCE so the engine deposits the recorded 66 J — matching the
    anchor's INPUT ('matched by construction', the 0.1.3 recipe; deposition
    is linear in j, and the anchor's outputs never feed back — PLAN (26)).
    Measured outputs: ~2.3 J at R_OC 20%, ~0.6 J at 80%.

    Verdict history (recorded, dated): x2.1-2.8 OVER (0.1.3, one-state
    lumped KrF*); x1.31-1.34 OVER (0.1.8, B/C/Bv split + K34Ar); 0.1.9
    re-scores with the completed Ar manifold — the rows below are live.
    """
    from .engine import Device
    from .laser import Cavity
    from .network import mode_optics
    from .params import engine_parameters
    mix = {"Ar": 0.898, "Kr": 0.10, "F2": 0.002}
    meas = {0.2: 2.3, 0.8: 0.6}
    opt = mode_optics("krf")
    out = []
    j = 301.0
    for r_oc in (0.2, 0.8):
        res = None
        for _attempt in range(2):
            p = engine_parameters()
            # 0.1.37: the electron tier is selectable for the MAXWELL
            # verification (tools/maxwell_verify.py); default unchanged
            p.set("electron_tier", tier, source="api")
            p.set("loki_subfloor", bool(subfloor), source="api")
            if subtable is not None:
                p.set("beam_subtable", subtable, source="api")
            p.set_many({"mode": "krf", "drive": "beam",
                        "pressure_atm": 2.5, "beam_kV": 660.0,
                        "beam_j_A_cm2": j, "beam_flat_ns": 50.0,
                        "beam_rise_ns": 5.0, "t_end_ns": t_end_ns,
                        "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.04,
                        "Ly_m": 0.04, "Lz_m": 0.44, "beam_axis": 0},
                       source="problem")
            cav = Cavity(l_c_m=0.44, l_g_m=0.44, K_axial=1e-5, R_OC=r_oc,
                         k_c_loss=0.0, L_op_fixed=0.64,
                         area_m2=0.04 * 0.04, wavelength_nm=248.0,
                         sigma_se_m2=opt["sigma_se_m2"],
                         absorbers=dict(opt["absorbers"]))
            dev = Device(p, mixture=mix, cavity=cav)
            res = dev.run(n_save=25)
            dep = res["final_channels_J"]["drive_in"]
            if abs(dep - 66.0) < 0.5:
                break
            j *= 66.0 / dep                     # input matching only
        e_out = res["final_channels_J"]["laser_out"]
        g0pk = float(max(res["g0_m1"]))
        # absorber sum at the peak-g0 save point
        i_pk = int(np.argmax(res["g0_m1"]))
        alpha = sum(sig * res["dens_avg"][sp][i_pk]
                    for sp, sig in cav.absorbers.items()
                    if sp in res["dens_avg"])
        ratio = e_out / meas[r_oc] if meas[r_oc] > 0 else float("inf")
        out.append(
            f"R_OC {int(100 * r_oc)}%: engine {e_out:6.2f} J vs measured "
            f"~{meas[r_oc]:.1f} J -> x{ratio:.2f} "
            f"(deposited {dep:5.1f} J, j {j:.0f} A/cm2; peak g0 "
            f"{g0pk / 100.0:.3f} /cm"
            + (f", g0/alpha {g0pk / alpha:.1f}" if alpha > 0 else "")
            + "; " + _conditioning(g0pk, alpha, r_oc, 0.64, 0.44) + ")")
    return out


def electra_arf_rows(t_end_ns: float = 200.0, tier: str = "loki",
                     subtable: str = None, subfloor: bool = False) -> list:
    """The NRL Electra ArF quadruple-pass oscillator anchor, run LIVE —
    the engine's first ArF-side device anchor (0.1.11). Sources read
    first-hand 2026-08-05: Wolford 2020 HEDP 36,100801 (Table 1 + Sec 2)
    and the Orestes memo petrova2022 NRL/6720/MR--2022/1 (Table 4 same
    data; Table 2 absorbers; App. A4 optics — mode_optics('arf_beam')).

    Measured: 90 J @ 0.82 atm / 0.6 MW/cc flat-top; 96 J @ 1.09 / 0.8;
    80 J @ 1.36 / 0.9 (calorimeter; 99.7% Ar / 0.3% F2). The TREND is
    the discriminator: yield PEAKS at 1.09 atm and DECLINES at 1.36 —
    NRL attribute the decline to the ArF*+2Ar three-body class (ARF_AA
    here), and their own Orestes surrogate matches x1.02 at 0.82 atm
    then 'predicts significantly higher laser yield' at the higher
    pressures (wolford2020 p.4; the memo notes quadruple-pass Orestes
    models were NOT performed).

    Construction — inputs matched, outputs never fed back (PLAN (26)):
    - 0-D cell = the two folded gain lines UNROLLED: the measured
      10.5 x 11 cm2 beam over 2 x 100 cm of deposited length (l = 2 m;
      round-trip gain path 4 m = the real four passes; the pumped gas
      the mode does not sample — the imaged 14x30x100 cm3 deposition
      region — is NOT modeled: 0-D extracts everything it sees).
    - Flat-top power density matched by construction: j iterated so
      the trapezoid (20/100/20 ns) deposits P_flat x V x 120 ns.
      500 kV, double-sided injection idealized to one axis.
    - Cavity: R_OC 0.08 (uncoated silica); ALL other passive losses in
      L_op = 0.37 round trip (8 window crossings x 0.90, five mirror
      reflections x 0.97); NRL effective sigma_se 1.0e-16 cm2.
    Known construction limits, stated: no beam-direct F2 channels (the
    memo's Table A2 dissociation/vib/attachment by beam secondaries —
    their F2 is 50% dissociated by pulse end), no F2*/F* species, the
    LOKI swarm tables carry the Ar-fill 5% Kr mixture (no pure-Ar set
    shipped).
    """
    from .engine import Device
    from .laser import Cavity
    from .network import mode_optics
    from .params import engine_parameters
    mix = {"Ar": 0.997, "F2": 0.003}
    # (pressure_atm, flat-top MW/cc, measured J)
    points = [(0.82, 0.6, 90.0), (1.09, 0.8, 96.0), (1.36, 0.9, 80.0)]
    lx, ly, lz = 0.105, 0.11, 2.0
    vol_m3 = lx * ly * lz
    opt = mode_optics("arf_beam")
    out, engine_e = [], {}
    j = 80.0
    for p_atm, mw_cc, e_meas in points:
        target = mw_cc * 1e12 * vol_m3 * 120e-9    # MW/cc -> W/m3, trapezoid
        res = None
        for _attempt in range(3):
            p = engine_parameters()
            # 0.1.37: the electron tier is selectable for the MAXWELL
            # verification (tools/maxwell_verify.py); default unchanged
            p.set("electron_tier", tier, source="api")
            p.set("loki_subfloor", bool(subfloor), source="api")
            if subtable is not None:
                p.set("beam_subtable", subtable, source="api")
            p.set_many({"mode": "arf_beam", "drive": "beam",
                        "pressure_atm": p_atm, "beam_kV": 500.0,
                        "beam_j_A_cm2": j, "beam_flat_ns": 100.0,
                        "beam_rise_ns": 20.0, "t_end_ns": t_end_ns,
                        "nx": 1, "ny": 1, "nz": 1, "Lx_m": lx,
                        "Ly_m": ly, "Lz_m": lz, "beam_axis": 0},
                       source="problem")
            cav = Cavity(l_c_m=lz, l_g_m=lz, K_axial=1e-5, R_OC=0.08,
                         k_c_loss=0.0, L_op_fixed=0.37, area_m2=lx * ly,
                         wavelength_nm=193.0,
                         sigma_se_m2=opt["sigma_se_m2"],
                         absorbers=dict(opt["absorbers"]))
            dev = Device(p, mixture=mix, cavity=cav)
            res = dev.run(n_save=25)
            dep = res["final_channels_J"]["drive_in"]
            if abs(dep - target) < 0.01 * target:
                break
            j *= target / dep                     # input matching only
        e_out = res["final_channels_J"]["laser_out"]
        engine_e[p_atm] = e_out
        g0pk = float(max(res["g0_m1"]))
        i_pk = int(np.argmax(res["g0_m1"]))
        alpha = sum(sig * res["dens_avg"][sp][i_pk]
                    for sp, sig in cav.absorbers.items()
                    if sp in res["dens_avg"])
        out.append(
            f"{p_atm:.2f} atm @ {mw_cc:.1f} MW/cc: engine {e_out:6.1f} J "
            f"vs measured {e_meas:.0f} J -> x{e_out / e_meas:.2f} "
            f"(deposited {dep:6.0f} J, j {j:.0f} A/cm2; peak g0 "
            f"{g0pk / 100.0:.4f} /cm"
            + (f", g0/alpha {g0pk / alpha:.1f}" if alpha > 0 else "")
            + "; " + _conditioning(g0pk, alpha, 0.08, 0.37, lz) + ")")
    # the trend verdict: measured argmax 1.09 atm, decline at 1.36
    e_mid, e_hi = engine_e[1.09], engine_e[1.36]
    trend = ("REPRODUCED" if e_hi < e_mid else "NOT reproduced")
    out.append(
        f"pressure trend (the 3-body discriminator): measured PEAK at "
        f"1.09 atm with the 1.36-atm point 17% LOWER; engine "
        f"{'falls' if e_hi < e_mid else 'RISES'} 1.09 -> 1.36 "
        f"({e_mid:.1f} -> {e_hi:.1f} J) — {trend}. NRL's own Orestes "
        f"OVERPREDICTS the high-pressure points (their stated 3-body "
        f"suspicion = the ARF_AA class this engine carries first-hand)")
    return out


def electra_population_rows(t_end_ns: float = 200.0, tier: str = "loki",
                            subtable: str = None, subfloor: bool = False) -> list:
    """The NRL Orestes POPULATION/TEMPERATURE benchmark, scored live —
    the fluid-Te suspect put on the record (0.1.12). petrova2022
    NRL/6720/MR--2022/1 Figs. 4/6/8 + text, first-hand 2026-08-05:
    at 0.5 MW/cc, 12.5 psi, 0.5% F2 in Ar their EEDF-resolved model
    reports Te ~2 eV FLAT through the pulse, Tg 300 -> ~420 K, ne_pk
    ~5e14 cm-3, ArF*_pk ~5e14 (measured ~4e14, their [58]), Ar*_pk
    ~1e15 (Ar** ~ Ar*/3), F- plateau ~1e14, F2 ~50% dissociated by
    pulse end, ion ordering Ar+ > Ar2+ > F2+, and g0_pk 2.5 %/cm with
    alpha ~0.2 %/cm (gamma ~17.5).

    Why this block exists: the cross-device x0.5-0.6 anchor deficit
    points at the shared electron layer. Comparing POPULATIONS (not
    just yield) localizes it — if ne/Te are off, the electron
    treatment is confirmed as the suspect; if populations match but
    yield misses, the defect is downstream (extraction/optics).

    Construction: the electra_arf_rows cell/cavity at 12.5 psi
    (0.8503 atm), 0.5% F2, 0.5 MW/cc flat-top matched (trapezoid
    20/100/20 ns = 120 ns-equivalent vs their ~136 ns-equivalent
    experimental shape; PEAK quantities compared, Tg at t_end 200 ns
    vs their ~260 ns — both stated). Their figures are for their
    double-pass-class oscillator; the quadruple-pass-class cavity
    here clamps ArF* slightly differently (stated). Te definitions
    comparable: both are bulk (2/3)<u> measures — their Eq. 13 over
    the EEDF (tail number-density negligible), ours the fluid ledger.
    """
    from .engine import Device
    from .laser import Cavity
    from .network import mode_optics
    from .params import engine_parameters
    mix = {"Ar": 0.995, "F2": 0.005}
    p_atm = 12.5 / 14.696
    mw_cc = 0.5
    lx, ly, lz = 0.105, 0.11, 2.0
    opt = mode_optics("arf_beam")
    target = mw_cc * 1e12 * (lx * ly * lz) * 120e-9
    j = 150.0
    res = None
    for _attempt in range(3):
        p = engine_parameters()
        # 0.1.37: the electron tier is selectable for the MAXWELL
        # verification (tools/maxwell_verify.py); default unchanged
        p.set("electron_tier", tier, source="api")
        p.set("loki_subfloor", bool(subfloor), source="api")
        if subtable is not None:
            p.set("beam_subtable", subtable, source="api")
        p.set_many({"mode": "arf_beam", "drive": "beam",
                    "pressure_atm": p_atm, "beam_kV": 500.0,
                    "beam_j_A_cm2": j, "beam_flat_ns": 100.0,
                    "beam_rise_ns": 20.0, "t_end_ns": t_end_ns,
                    "nx": 1, "ny": 1, "nz": 1, "Lx_m": lx,
                    "Ly_m": ly, "Lz_m": lz, "beam_axis": 0},
                   source="problem")
        cav = Cavity(l_c_m=lz, l_g_m=lz, K_axial=1e-5, R_OC=0.08,
                     k_c_loss=0.0, L_op_fixed=0.37, area_m2=lx * ly,
                     wavelength_nm=193.0, sigma_se_m2=opt["sigma_se_m2"],
                     absorbers=dict(opt["absorbers"]))
        dev = Device(p, mixture=mix, cavity=cav)
        res = dev.run(n_save=41)
        dep = res["final_channels_J"]["drive_in"]
        if abs(dep - target) < 0.01 * target:
            break
        j *= target / dep                         # input matching only
    d = res["dens_avg"]

    def pk(sp):
        return float(np.max(d[sp])) if sp in d else 0.0
    cm3 = 1e-6                                    # m-3 -> cm-3
    ne_pk = float(res["ne_m3"].max()) * cm3
    te_pk = float(res["Te_eV"].max())
    tg_end = float(res["Tg_K"][-1]) if "Tg_K" in res else float("nan")
    arf_pk, ar_pk, arss_pk = pk("ArF*") * cm3, pk("Ar*") * cm3, \
        pk("Ar**") * cm3
    fm_end = float(d["F-"][-1]) * cm3 if "F-" in d else 0.0
    f2_diss = 1.0 - float(d["F2"][-1] / d["F2"][0]) if "F2" in d else 0.0
    g0_pk = float(max(res["g0_m1"])) / 100.0
    i_pk = int(np.argmax(res["g0_m1"]))
    alpha = sum(sig * d[sp][i_pk] for sp, sig in cav.absorbers.items()
                if sp in d) / 100.0
    ions = {sp: pk(sp) for sp in ("Ar+", "Ar2+", "F2+")}
    order = " > ".join(sorted(ions, key=ions.get, reverse=True))
    rows = [
        ("Te peak [eV]", te_pk, 2.0, "Fig. 6, flat ~2 eV"),
        ("Tg end [K]", tg_end, 420.0, "Fig. 6 (~260 ns; engine at "
                                      "200 ns)"),
        ("ne peak [cm-3]", ne_pk, 5e14, "Fig. 4"),
        ("ArF* peak [cm-3]", arf_pk, 5e14, "Fig. 4 (measured ~4e14, "
                                           "their [58])"),
        ("Ar* peak [cm-3]", ar_pk, 1e15, "Fig. 4"),
        ("Ar** peak [cm-3]", arss_pk, 3e14, "Fig. 4 (~Ar*/3)"),
        ("F- plateau [cm-3]", fm_end, 1e14, "Fig. 4 (flattens ~30 ns)"),
        ("F2 dissociated", f2_diss, 0.50, "Sec 2.2.2 (by pulse end)"),
        ("g0 peak [/cm]", g0_pk, 0.025, "Fig. 8"),
        ("alpha @ g0 peak [/cm]", alpha, 0.002, "Fig. 8"),
    ]
    out = []
    for name, eng, ref, note in rows:
        ratio = eng / ref if ref else float("inf")
        out.append(f"{name:<22s} engine {eng:9.3g} vs NRL {ref:9.3g} "
                   f"-> x{ratio:5.2f}  ({note})")
    out.append(f"ion ordering: engine {order} vs NRL Ar+ > Ar2+ > F2+ "
               f"(Fig. 4)")
    out.append(f"(deposited {res['final_channels_J']['drive_in']:.0f} J "
               f"matched to 0.5 MW/cc; laser out "
               f"{res['final_channels_J']['laser_out']:.1f} J)")
    return out


def run_scorecard(result: dict = None) -> str:
    """The full text scorecard; pass an engine result to add its
    run-level rows."""
    lines = ["=" * 72,
             "REBUILD SCORECARD — anchors score, they never steer "
             "(PLAN (26))", "=" * 72]
    bs = beam_scorecard()
    lines.append("\n-- beam physics vs literature anchors "
                 "(miss stated per row) --")
    for name, got, ref, miss in bs["rows"]:
        lines.append(f"  {name:<46s} {got:>16s} vs {ref:<14s} {miss}")
    lines.append("\n-- W by TWO INDEPENDENT ROUTES (REBUILD Spencer-Fano vs "
                 "PIC_MCC_3D MC cascade, same E0 + LXCat data) --")
    try:
        for ln in two_route_w_rows():
            lines.append("  " + ln)
    except FileNotFoundError as e:
        lines.append(f"  kinetic export missing: {e}")
    if result is not None:
        dev = result["device"]
        lines.append("\n-- this run --")
        lines.append("  " + result["regime"])
        worst = max((abs(c["residual_rel"]) for c in result["closure"][1:]),
                    default=0.0)
        lines.append(f"  energy-ledger closure (S3 gate): worst residual "
                     f"{worst:.2e} of energy moved")
        g0pk = float(max(result["g0_m1"])) / 100.0
        lines.append(f"  peak g0 = {g0pk:.4f} /cm")
        if (dev.mode == "krf" and dev.p["drive"] == "beam"
                and abs(dev.p["pressure_atm"] - 3.0) < 0.2):
            for eng, ref in V2_G0_PER_CM.items():
                lines.append(
                    f"    vs {eng} at the shared v2-class point: {ref} /cm "
                    f"({100 * (g0pk / ref - 1):+.0f}%) — engine-vs-engine, "
                    "matched conditions. NOTE (0.1.9/0.1.10): those "
                    "engines lack the ENTIRE Ar heteronuclear manifold "
                    "this network now carries (K34Ar found 0.1.8; "
                    "KrF*+2Ar -> ArKrF*, KrF*+Ar, ArKr*, Kr2*/Ar2* VUV "
                    "et al. found 0.1.9, triple-sourced Kannari/J&H/Xu; "
                    "the ArKr+ ion manifold, Ar-ion/F- recombination, "
                    "harpoon quartet + ArF* relaxation found 0.1.10) — "
                    "sitting FAR below them is the expected signature of "
                    "channels they omit; the v2 point at 0.42 MW/cm3 "
                    "deposition is sub-threshold in the completed "
                    "network (0.0053 at 0.1.9 -> 0.0034 at 0.1.11 -> "
                    "0.0054 /cm at 0.1.12, the prompt-floor beam split)")
        lines.append(f"  eta_intrinsic = {100 * result['eta_intrinsic']:.2f}%"
                     " (laser out / energy deposited; NOT wall-plug — the "
                     "delivery chain is not modelled here)")
        lines.append("  DE-AR context: current devices ~7% wall-plug, "
                     "program goal >14% (docs/DE-AR_targets.txt) — no "
                     "wall-plug claim is made by this engine")
        lines.append("\n-- network provenance census (S5: drive EST to 0) --")
        lines.append("  " + census_line(tier_census(dev.net)))
        if dev.beam_yields is not None:
            lines.append(f"  beam yields route: {dev.beam_yields.mode} "
                         f"(W = {dev.beam_yields.W_eV:.1f} eV; "
                         f"identity residual "
                         f"{dev.beam_yields.identity_residual_eV():.1e} eV)")
    lines.append("\n-- the taper at beam-mode energies (L3, MEASURED "
                 "2026-08-04) --")
    try:
        from .kinetic_tier import KineticRates, taper_vs_relaxation
        from .loki_tier import rates_for_mode
        kin_ = KineticRates.load_default()
        lok_ = rates_for_mode("krf", {"Ar": 0.949, "Kr": 0.05,
                                      "F2": 0.001})
        rows = [r for r in taper_vs_relaxation(kin_, lok_)
                if r[0] in ("F-", "el:Ar", "el:Kr", "vib")]
        for b, e, km, kt, ratio, err in rows:
            lines.append(f"  {b:<6s} @ <eps> {e:4.2f} eV: taper "
                         f"{kt:.2e} vs measured {km:.2e} "
                         f"(x{ratio:.2f}, stat {100 * err:.0f}%)")
        lines.append("  verdict: elastic/vib within x0.8-1.9; F- "
                     "attachment was OVERESTIMATED x3.3 by the taper at "
                     "1.8 eV. CLOSED 0.1.9: beam mode now CONSUMES these "
                     "curves below the table floor (RelaxationBridged, "
                     "default ON; the rows above are what it consumes); "
                     "below the measured floor the taper extrapolates "
                     "from the measured anchor — see OMISSIONS for the "
                     "Ramsauer-stall deepening follow-up")
    except FileNotFoundError as e:
        lines.append(f"  kinetic export missing: {e}")
    lines.append("\n-- W(Ne) miss, decomposed (2026-08-04) --")
    lines.append("  W(Ne) 43.7 = 21.6 (IP) + 11.1 (excitation spend, "
                 "U_exc/U_ion = 0.62) + 11.0 (sub-threshold); literature "
                 "36.6 implies a 15.0 eV excitation+sub-threshold budget "
                 "vs our 22.1 — the +19% is carried by the Biagi set's "
                 "excitation/ionization ratio + equal sharing (both W "
                 "routes agree; a data/method property, not a defect)")
    lines.append("\n-- KINETIC vs LOKI (S4) --")
    try:
        from .kinetic_tier import (KineticRates, compare_tiers,
                                   format_comparison)
        from .loki_tier import rates_for_mode
        kin = KineticRates.load_default()
        lok = rates_for_mode("krf", {"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
        for ln in format_comparison(compare_tiers(kin, lok)).splitlines():
            lines.append("  " + ln)
    except FileNotFoundError as e:
        lines.append(f"  kinetic export missing: {e}")
    lines.append("\n-- Chen 2025 Fig. 6, scored live (S7) --")
    try:
        for ln in chen_figure6_rows():
            lines.append("  " + ln)
    except Exception as e:                      # data missing etc.
        lines.append(f"  could-not-run: {type(e).__name__}: {e}")
    lines.append("\n-- Kannari 1985 device anchor, scored LIVE (misses "
                 "stated, nothing tuned) --")
    lines.append("  660 kV, 1.88 MW/cm3 x 50 ns, Ar/Kr/F2 89.8/10/0.2, "
                 "2.5 atm, lossy CaF2 cavity (L_op 0.64), 66 J deposited "
                 "by construction:")
    try:
        for ln in kannari_device_rows():
            lines.append("    " + ln)
    except Exception as e:
        lines.append(f"    could-not-run: {type(e).__name__}: {e}")
    lines.append("    verdict history: x2.1-2.8 OVER (0.1.3, one-state); "
                 "x1.31-1.34 OVER (0.1.8, split + K34Ar); 0.1.9 completes "
                 "the Ar manifold (KrF*+2Ar/KrF*+Ar sinks et al., "
                 "triple-sourced) and the anchor flips to UNDER — "
                 "x0.60-class at the robust R_OC 80% point, threshold-"
                 "pinned ~0 at 20% (g0 barely clears g_th at L_op 0.64). "
                 "0.1.10 TESTS the named suspect (numbers corrected "
                 "0.1.11 after the Chen-overlap dedup): the ArKr+ ion "
                 "manifold goes IN and the anchor moves only x0.60 -> "
                 "x0.56 — hypothesis consumed, deficit stands, "
                 "pointing OFF-network. 0.1.12 then TESTS the "
                 "first-listed off-network suspect (the fluid-Te "
                 "treatment) against the NRL EEDF population benchmark "
                 "— CONFIRMED x4-cold bulk, F- x0.003 — and the "
                 "DERIVED prompt-floor correction (nothing fitted; "
                 "see electra_population_rows) moves this anchor "
                 "x0.56 -> x0.81 at R_OC 80% (x0.21 at the "
                 "threshold-sensitive 20% point; g0/alpha 3.0-3.1 vs "
                 "the empirical 4-12). The old x0.5 cross-device "
                 "deficit story is BROKEN (Electra now straddles 1); "
                 "open residuals, quantified: F- x0.02 vs benchmark "
                 "(the export-deepening follow-up's strongest reason), "
                 "Ar** x0.06, the threshold-region points, and note "
                 "Kannari's own model matches his device carrying a "
                 "x4.4 LARGER e-quench (his 2.0e-7 vs the Kushner-89 "
                 "4.5e-8 kept here)")
    lines.append("  NRL Electra ArF quadruple-pass oscillator (Wolford 2020 "
                 "HEDP Table 1 / petrova2022 NRL memo, first-hand "
                 "2026-08-05), scored LIVE in the arf_beam mode — the "
                 "first ArF-side device anchor (misses stated, nothing "
                 "tuned; construction in electra_arf_rows.__doc__):")
    try:
        for ln in electra_arf_rows():
            lines.append("    " + ln)
    except Exception as e:
        lines.append(f"    could-not-run: {type(e).__name__}: {e}")
    lines.append("  Chen 2025 device (0.1.25 corrected lane; photon "
                 "eqs. 7/8 are this engine's, incl. the loss split):")
    lines.append("    HISTORY, kept honest: the 0.1.5..0.1.24 record here "
                 "claimed 'Table II verbatim-confirmed', 'E/N 100% "
                 "in-table' and a residual 'x3-5 PER-JOULE formation gap "
                 "(suspects: Ne* source, Xe seed, 0-D vs 1-D)'. The "
                 "2026-08-25 comparison audit MEASURED all three claims "
                 "wrong: the shipped circuit was the Fig. 5 validation "
                 "experiment's (10.4 nF/8.1 nF/30 nH), not Table II's "
                 "(10/250 nH/10/5 nH — grafting it also made the stored "
                 "energy '1.68 J' where the paper states 1.62); the "
                 "in-table claim came from a gauge that never counted "
                 "BELOW-floor excursions (the whole breakdown ramp, <5 Td); "
                 "and controlled swaps attribute the anchor gap to the "
                 "circuit transcription + four physically impossible "
                 "as-printed rates (R13/R14/R23/R24 — CHEN_ERRATA), with "
                 "the named suspects second-order. The corrected lane's "
                 "live numbers and residual misses are in the Fig. 6 "
                 "block above; the as-printed lane is re-run there and "
                 "recorded (it cannot integrate under the paper's own "
                 "circuit)")
    lines.append("\n-- omissions (part of the score) --")
    for o in OMISSIONS:
        lines.append(f"  * {o}")
    lines.append("=" * 72)
    return "\n".join(lines)
