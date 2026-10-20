"""network.py gates — D2's mechanical balance, WITH TEETH.

The S2 gate demands a deliberately broken row FAIL (ORESTES_CLASS_REBUILD
§6: "prove the gate has teeth"). Every teeth test here constructs the exact
defect class the history recorded and asserts the mechanical check kills it.
"""

import pytest

from REBUILD.MODULES.network import (NetworkError, Reaction, build_network,
                                     check_network, check_reaction,
                                     kannari_ps_relaxation, mode_optics)
from REBUILD.MODULES.provenance import tier_census

NE_BUCKETS = {"Ne*", "Ne**", "Ne+", "F-", "F2(a3)", "F2(A1)", "F2+",
              "Kr*", "Kr**", "Kr+"}
AR_BUCKETS = NE_BUCKETS | {"Ar*", "Ar**", "Ar+"}


def test_networks_build_and_balance():
    for mode, buckets in (("krf", NE_BUCKETS), ("krf", AR_BUCKETS),
                          ("arf_chen", AR_BUCKETS - {"Kr*", "Kr**", "Kr+"})):
        rows = build_network(mode, buckets)
        assert len(rows) >= (80 if mode == "krf" else 38)
        for r in rows:
            assert r.cite, r.rid


def test_chen_errata_lane_default_and_as_printed():
    """0.1.25: the errata lane (default) carries the adjudicated values on
    the CHEN rows with provenance in the cite; errata=False reproduces the
    printed table exactly (the record of the PAPER). The corrections apply
    in every mode that includes chen_rows — in krf, R24 as printed
    (4e-40 m6/s = 4e-28 cm6/s) annihilated Ne* in ~0.4 ps in any Ne fill."""
    arf_b = AR_BUCKETS - {"Kr*", "Kr**", "Kr+"}
    cor = {r.rid: r for r in build_network("arf_chen", arf_b)}
    raw = {r.rid: r for r in build_network("arf_chen", arf_b,
                                           chen_errata=False)}
    from REBUILD.MODULES.chen2025 import CHEN_ERRATA
    for rid, (fixed, printed, _basis) in CHEN_ERRATA.items():
        assert cor[rid].rate[1] == fixed, rid
        assert raw[rid].rate[1] == printed, rid
        assert "ERRATA" in cor[rid].cite and "ERRATA" not in raw[rid].cite
    # exponents/products untouched by the errata (magnitude-only)
    assert cor["R13"].rate[2:] == raw["R13"].rate[2:]
    # krf: R13/R14 are dropped by the repulsive reroute either way, and
    # the Ne-partnered three-body rows carry the corrected magnitudes
    kr = {r.rid: r for r in build_network("krf", AR_BUCKETS)}
    assert kr["R24"].rate[1] == CHEN_ERRATA["R24"][0]
    assert kr["R23"].rate[1] == CHEN_ERRATA["R23"][0]


def test_teeth_mass_imbalance_fails():
    """The r115 class: a term with wrong species content cannot be written.
    (The original sign error survived phantom11's parity harness because
    every parity case had Ar_g = 0 — root CHANGELOG Known-bad.)"""
    bad = Reaction("BAD1", ("Kr*", "F2"), ("KrF*",), ("const", 1e-16),
                   "EST", "deliberately drops an F atom")
    with pytest.raises(NetworkError, match="element balance"):
        check_reaction(bad)


def test_teeth_charge_imbalance_fails():
    bad = Reaction("BAD2", ("Kr+", "F2"), ("Kr", "F2"), ("const", 1e-16),
                   "EST", "deliberately loses a charge")
    with pytest.raises(NetworkError, match="charge balance"):
        check_reaction(bad)


def test_teeth_over_radiating_row_fails():
    bad = Reaction("BAD3", ("Kr2F*",), ("Kr", "Kr", "F"), ("const", 1e6),
                   "EST", "radiates more than it stores", kind="rad",
                   hv_eV=99.0)
    with pytest.raises(NetworkError, match="exceeds"):
        check_reaction(bad)


def test_teeth_missing_citation_fails():
    bad = Reaction("BAD4", ("Kr*", "F2"), ("KrF*", "F"), ("const", 1e-16),
                   "EST", "")
    with pytest.raises(NetworkError, match="citation"):
        check_reaction(bad)


def test_teeth_duplicate_rid_fails():
    r = Reaction("DUP", ("Kr*", "F2"), ("KrF*", "F"), ("const", 1e-16),
                 "EST", "x")
    with pytest.raises(NetworkError, match="duplicate"):
        check_network([r, r])


def test_kannari_ps_rows_are_unconditional():
    """PLAN (72): the p->s relaxation is IN, with no flag to turn it off.
    In phantom11 it defaults OFF as a parity artifact worth 3.4x in baseline
    light output and a breached W-value bound."""
    rows = build_network("krf", AR_BUCKETS)
    rids = {r.rid for r in rows}
    for rid in ("KPS_Ar", "KPS_Ne", "KPS_Kr", "APS_Ar", "APS_Ne"):
        assert rid in rids, rid
    for r in kannari_ps_relaxation():
        assert r.rate == ("const", 1.0e-16)
        # cite is Kannari's Table I or (for APS_Ne since the 0.1.7 pass)
        # Lee 1989 [37] — the direct measurement class outranks the analog
        assert "Kannari" in r.cite or "Lee 1989" in r.cite


def _code_only(path):
    """Source with strings and comments stripped — docstrings that DESCRIBE
    a defect must not trip a scan for the defect."""
    import io
    import tokenize
    out = []
    with open(path, "rb") as f:
        for tok in tokenize.tokenize(f.readline):
            if tok.type in (tokenize.STRING, tokenize.COMMENT):
                continue
            out.append(tok.string)
    return " ".join(out)


def test_no_parity_flag_exists_in_source():
    """D1: no physics default may be gated by a code-equivalence flag.
    (The scan strips docstrings/comments: network.py's Kannari docstring
    NAMES LIT_QUENCH as the defect it closes, which is documentation, not
    a flag.)"""
    from pathlib import Path
    mod_dir = Path(__file__).resolve().parents[1] / "MODULES"
    for f in mod_dir.glob("*.py"):
        code = _code_only(f)
        assert "LIT_QUENCH" not in code, f.name
        assert "og_pyfiles" not in code, f.name


def test_tier_census_recorded():
    rows = build_network("krf", AR_BUCKETS)
    c = tier_census(rows)
    assert c["CHEN"] >= 13 and c["LOKI"] >= 10 and c["REPO"] >= 15
    # CHEN 23 -> 15 in krf at 0.1.11: 8 generated Chen rows are
    # superseded by the first-hand set (CHEN_SUPERSEDED — the 0.1.10
    # double-count fix); 15 -> 13 at 0.1.14: R13/R14 dropped by the
    # repulsive-F2 reroute (excitation IS dissociation; the two-step
    # bookkeeping stays in arf_chen only). arf_chen keeps all of
    # Chen's verbatim.
    # REPO fell 20+ -> 19 in 0.1.9: K34/K32/AK8 upgraded REPO -> LIT
    # (Xu 1993 / Kannari first-hand) — the intended direction. 19 -> 22
    # in 0.1.10: AK7/R11p left for LIT (two-primary consensus) while the
    # four ion-ion rows whose Flannery params mirror K26/K108 arrived
    # honestly REPO (AKI_F x2 split, AKI_FA, AIIR1, AIIR2).
    # S5 target: EST -> 0. Record the current count so a silent GROWTH
    # fails; driving it down is progress and only moves this pin.
    # 22 -> 18 (Kannari) -> 23 (B/C/Bv split additions) -> 19 on
    # 2026-08-04: the Lee/Boichenko first-hand pass retired NK2, HK1,
    # HK4, APS_Ne (CHANGELOG [0.1.7]). Still 19 after 0.1.9: the entire
    # Ar-manifold completion (30 rows) arrived fully cited — LIT 40 -> 78
    # — and the x1/3..x3 bracket sweep now PRICES the remainder
    # (scorecard.est_bracket_rows; measured verdict in CHANGELOG 0.1.9).
    # Still 19 after 0.1.10: all 28 ion-manifold/completion rows carry
    # primary prints — LIT 78 -> 103 (census 149 -> 177 rows; 169 after
    # the 0.1.11 dedup removed the 8 double-counted Chen channels).
    # 19 -> 24 at 0.1.38: five Kr-partner analogs (BC_Kr, CB_Kr, VREL_Kr,
    # APS_Kr, AK8_Kr) — the buffer-partner policy of kannari_ps_relaxation
    # applied to the rows the 0.1.4 split left Ar/Ne-only. Found by the
    # Mandl-1985 pure-Kr/F2 check (a fill with no Ar/Ne could not lase);
    # each row names its primary (Gadomski 1992, Kvaran 1988, Kannari
    # 1985) for the EST -> LIT move. CHANGELOG [0.1.38].
    assert c["EST"] <= 24, "EST tier grew — every new row needs a citation"
    assert c["LIT"] >= 100


def test_ar_manifold_complete_and_sourced():
    """The 0.1.9 completion (found via the Xu 1993 re-verification pass):
    the KrF*+2Ar/KrF*+Ar sinks, the ArKr*/ArKrF* chains, the Ar2F* fate
    block and the dimer VUV rows exist per state, fully cited — the
    formerly-missing dominant channels cannot silently vanish again."""
    rows = build_network("krf", AR_BUCKETS)
    by_rid = {r.rid: r for r in rows}
    # the headline sinks, per B/C/Bv state
    for base in ("KAA", "KQAr"):
        for st in ("_B", "_C", "_Bv"):
            assert base + st in by_rid, base + st
    # KAA magnitude: 7.0e-32 cm6/s (Kannari [62,63]; J&H #14 8.0e-32;
    # Xu (2-8)e-32) — a transcription slip of one decade would be the
    # single largest error available in this network
    assert abs(by_rid["KAA_B"].rate[1] / 7.0e-44 - 1.0) < 1e-12
    # every chain link present
    for rid in ("AKD1", "AKD2", "AKD3_fBv", "AKD3_fC", "AKD4", "AKD_RAD",
                "AKF_F2", "AKF_KR", "AKF_AR", "AKF_RAD", "A2F_F2",
                "A2F_KR1", "A2F_KR2_fBv", "A2F_KR2_fC", "A2F_AR",
                "A2F_RAD", "KR2_RAD", "AR2_RAD", "K2F_BR", "A2_BR",
                "K2F_F_fBv", "K2F_F_fC", "A2_F", "AK9", "ADR"):
        assert rid in by_rid, rid
        assert by_rid[rid].tier == "LIT" and by_rid[rid].cite
    # recorded inter-source conflicts stay recorded (the x24/x27/x10
    # radiative spreads are load-bearing honesty)
    for rid in ("AKD_RAD", "KR2_RAD", "AR2_RAD", "AKF_RAD"):
        assert "CONFLICT" in by_rid[rid].cite, rid
    # Ar2F* is no longer a sink-less absorber: it has consumers
    consumers = [r for r in rows if "Ar2F*" in r.reactants]
    assert len(consumers) >= 4


def test_ar_ion_manifold_complete_and_sourced():
    """The 0.1.10 completion — the 0.1.9 'next network milestone': the
    ArKr+ ion manifold (Kannari 1985 Table I p.4311-4312 CT/ion-ion/DR/
    Penning/superelastic blocks; J&H 1980 #10,36-40), the Ar-ion + F-
    recombination pair, the source-less-Ar2+/sink-less-F2+ fixes, the
    Ar*/Ar** harpoon quartet and the ArF* relaxation block. Every rid
    pinned here so none can silently vanish."""
    from REBUILD.MODULES.species import REGISTRY, KRF_MODE_SPECIES
    rows = build_network("krf", AR_BUCKETS)
    by_rid = {r.rid: r for r in rows}
    # the species itself: registered, +1, Ar+Kr, in the krf state vector
    assert "ArKr+" in REGISTRY and "ArKr+" in KRF_MODE_SPECIES
    sp = REGISTRY["ArKr+"]
    assert sp.charge == +1 and dict(sp.atoms) == {"Ar": 1, "Kr": 1}
    # AKI1 magnitude: 1.0e-31 cm6/s (Kannari [16]; J&H #37 same) — the
    # dominant Kr+ route in an Ar fill; a decade slip here would be the
    # largest formation-side error available
    assert abs(by_rid["AKI1"].rate[1] / 1.0e-43 - 1.0) < 1e-12
    # every manifold link present, tiered, cited
    lit = ("AKI1", "AKI2", "AKI3", "AKI_DR", "AKI_ION", "AKI_SE",
           "AKF_SE", "AAR1", "FMN1", "PEN_KR", "PEN_AR", "PEN_ARSS",
           "PEN_KR2", "PEN_AR2", "AF1", "AF2", "AF3", "AF4", "AKD5",
           "ARF_AA", "ARF_F2", "ARF_AR", "ARF_RAD")
    for rid in lit:
        assert rid in by_rid, rid
        assert by_rid[rid].tier == "LIT" and by_rid[rid].cite, rid
    # the ion-ion rows: existence printed, params mirrored (REPO), and
    # the 50/50 ArKr+ branch split really sums to the K108 form
    # (compared PRE-split — K108 itself forms KrF* and splits to _fBv/_fC)
    for rid in ("AKI_F_fBv", "AKI_F_fC", "AKI_FA", "AIIR1", "AIIR2"):
        assert rid in by_rid and by_rid[rid].tier == "REPO", rid
    from REBUILD.MODULES.network import ar_kr_ion_manifold, kr_core
    raw = {r.rid: r for r in ar_kr_ion_manifold()}
    k108 = next(r for r in kr_core() if r.rid == "K108")
    total = raw["AKI_F"].rate[2] + raw["AKI_FA"].rate[2]
    assert abs(total / k108.rate[2] - 1.0) < 1e-12
    assert raw["AKI_F"].rate[1] == raw["AKI_FA"].rate[1] == k108.rate[1]
    # ArKr+ has both sources and sinks (no source-less/sink-less repeat)
    srcs = [r for r in rows if "ArKr+" in r.products]
    snks = [r for r in rows if "ArKr+" in r.reactants]
    assert len(srcs) >= 3 and len(snks) >= 4
    # the ArKr* photoionization absorber: adopted, routed, energetically
    # open (E_form + 5.0 eV >= E_form(ArKr+)); and NO invented sigma for
    # the ion itself — neither primary prints one
    from REBUILD.MODULES.network import mode_optics, PHOTO_PRODUCTS
    opt = mode_optics("krf")
    assert abs(opt["absorbers"]["ArKr*"] / 1.5e-21 - 1.0) < 1e-12
    assert PHOTO_PRODUCTS["ArKr*"] == ("ArKr+", "e")
    assert REGISTRY["ArKr*"].e_form_eV + 5.0 >= REGISTRY["ArKr+"].e_form_eV
    assert "ArKr+" not in opt["absorbers"]
    # the two-primary value re-pins hold (K15-precedent consensus)
    assert abs(by_rid["AK7"].rate[1] / 2.5e-43 - 1.0) < 1e-12
    assert abs(by_rid["AK2"].rate[1] / 7.5e-16 - 1.0) < 1e-12
    # R11p: dominant printed branch (Kr*), ADR-convention total
    assert by_rid["R11p"].products == ("Kr*", "Kr")
    assert abs(by_rid["R11p"].rate[1] / 2.24e-13 - 1.0) < 1e-12


def test_chen_supersede_and_channel_uniqueness():
    """The 0.1.11 double-count fix, mechanically enforced both ways:
    in krf the 8 superseded Chen rows are OUT (their first-hand
    successors in), in arf_chen Chen's set is verbatim (that anchor
    scores THEIR model), and check_network refuses two rows on one
    physical channel — the defect class that ran 8 channels at 2x
    (found 2026-08-05: row existence had been grepped from SOURCE while
    the CHEN block is GENERATED from chen2025.TABLE_III)."""
    from REBUILD.MODULES.network import (CHEN_SUPERSEDED, NetworkError,
                                         Reaction, check_network)
    rows = build_network("krf", AR_BUCKETS)
    rids = {r.rid for r in rows}
    for gone, kept in CHEN_SUPERSEDED.items():
        assert gone not in rids, f"{gone} still in krf"
        base = kept if kept in rids else kept + "_fBv"
        assert base in rids, f"successor {kept} missing in krf"
    arf = build_network("arf_chen", AR_BUCKETS - {"Kr*", "Kr**", "Kr+"})
    arf_rids = {r.rid for r in arf}
    for gone in CHEN_SUPERSEDED:
        assert gone in arf_rids, f"{gone} missing from arf_chen"
    # the teeth: an identical-channel pair must raise
    dup = [Reaction("T1", ("Ar*", "F2"), ("ArF*", "F"), ("const", 1e-16),
                    "EST", "teeth probe"),
           Reaction("T2", ("F2", "Ar*"), ("F", "ArF*"), ("const", 2e-16),
                    "EST", "teeth probe, permuted order")]
    try:
        check_network(dup)
        assert False, "duplicate channel not caught"
    except NetworkError as e:
        assert "duplicate channel" in str(e)


def test_repulsive_f2_reroute_and_screening():
    """The 0.1.14 chain-3 gates. (a) F2's repulsive a3/A1 states: in the
    Kannari-chemistry modes the measured excitation rows deliver F+F
    directly (electron pays the full threshold via e_cost_override) and
    the two-step R13/R14 are gone; arf_chen keeps Chen VERBATIM. The
    R11/R12 pair is the sole ALLOWED_SAME_CHANNEL exemption — two
    repulsive states, two measured rates, same fragments. (b) The
    ion-ion Debye screened-capture factor: parameter-free, 0.695 at
    Kannari's printed point (Morgan's MC prints 0.397 — the x1.75
    mildness is the RECORDED spread), -> 1 as the plasma vanishes."""
    import numpy as np
    from REBUILD.MODULES.chemistry import ion_ion_screening
    rows = build_network("krf", AR_BUCKETS)
    by = {r.rid: r for r in rows}
    assert "R13" not in by and "R14" not in by
    assert by["R11"].products == ("F", "F", "e")
    assert abs(by["R11"].e_cost_override_eV - 1.56) < 1e-12
    assert abs(by["R12"].e_cost_override_eV - 2.74) < 1e-12
    arf = {r.rid: r for r in
           build_network("arf_chen", AR_BUCKETS - {"Kr*", "Kr**", "Kr+"})}
    assert arf["R11"].products == ("F2(a3)", "e") and "R13" in arf
    # chain-2 rows present with superelastic routing
    assert by["AR_EDX1"].energy_to == "electrons"
    assert by["AR_EDX2"].tier == "LIT"
    from REBUILD.MODULES.scorecard import CONFLICT_ALTERNATES
    assert CONFLICT_ALTERNATES.get("APS_Ar") == 0.1
    # the screening factor: printed-point check + limits
    s_kannari = float(ion_ion_screening(np.array([1.0e20]),
                                        np.array([0.0]),
                                        np.array([300.0]))[0])
    assert abs(s_kannari - 0.695) < 0.01     # vs Morgan 0.397, recorded
    s_lo = float(ion_ion_screening(np.array([1.0e10]), np.array([0.0]),
                                   np.array([300.0]))[0])
    assert s_lo > 0.999                       # no plasma, no screening
    assert s_kannari < 1.0


def test_arf_beam_mode_and_electra_optics():
    """The Electra ArF anchor mode (0.1.11): the arf_beam assembly is the
    Kannari pool mechanically species-filtered — no Kr leakage, the
    krf dedup rule applied — and its optics are NRL Orestes' PRINTED
    193-nm set (petrova2022 Table 2 / App. A4, first-hand 2026-08-05)."""
    from REBUILD.MODULES.network import CHEN_SUPERSEDED, PHOTO_PRODUCTS
    from REBUILD.MODULES.species import ARF_MODE_SPECIES, REGISTRY
    rows = build_network("arf_beam", AR_BUCKETS - {"Kr*", "Kr**", "Kr+"})
    allowed = set(ARF_MODE_SPECIES)
    for r in rows:
        for s in r.reactants + r.products:
            assert (s in allowed or s == "e" or s.startswith("hv")), \
                f"{r.rid}: species {s} leaked past the arf filter"
    rids = {r.rid for r in rows}
    for gone, kept in CHEN_SUPERSEDED.items():
        assert gone not in rids
        assert kept in rids, f"successor {kept} missing in arf_beam"
    for rid in ("AF1", "AK5", "AK6", "ADR", "PEN_AR", "A2_BR", "AR2_RAD",
                "ARF_AA", "AAR1", "AIIR1", "AIIR2", "FMN1"):
        assert rid in rids, rid
    # NRL Table 2 optics, pinned (cm2 prints -> m2)
    opt = mode_optics("arf_beam")
    assert abs(opt["sigma_se_m2"] / 1.0e-20 - 1.0) < 1e-12
    assert abs(opt["absorbers"]["Ar**"] / 3.0e-22 - 1.0) < 1e-12
    assert abs(opt["absorbers"]["F-"] / 6.5e-22 - 1.0) < 1e-12
    assert set(opt["absorbers"]) == {"F2", "Ar*", "Ar**", "Ar2*",
                                     "Ar2F*", "Ar2+", "F-"}
    # every absorber routed and energetically open at the 193 nm quantum
    hv = 1239.841984 / opt["wavelength_nm"]
    for sp in opt["absorbers"]:
        prods = PHOTO_PRODUCTS[sp]
        de = (sum(REGISTRY[q].e_form_eV for q in prods if q != "e")
              - REGISTRY[sp].e_form_eV)
        assert de <= hv + 1e-9, f"{sp}: {de} eV > {hv}"


# EST rows whose reactants are all present in the Ar/Kr/F2 fill: the S5 bracket
# levers. 0.1.9 had four; 0.1.38 added the five Kr-partner rows (priced in
# CHANGELOG 0.1.38 — they carry krypton, so no fill binds them).
EST_ACTIVE_AR_FILL = {"R8", "KSTEP", "KPS_Kr", "URAD_Bv",
                      "BC_Kr", "CB_Kr", "VREL_Kr", "APS_Kr", "AK8_Kr"}


def test_est_bracket_machinery():
    """S5 closure mechanics: every EST row's rate spec is bracketable,
    and the element-bound classification is structurally right — in the
    Ar fill exactly the Ne-chain analogs are bound (their reactants need
    an element the fill lacks), leaving the active levers: the four of
    0.1.9 plus the five Kr-partner rows of 0.1.38 (krypton is in every
    fill, so they are active everywhere)."""
    from REBUILD.MODULES.scorecard import (CONFLICT_ALTERNATES,
                                           _bracket_rate,
                                           est_active_and_bound)
    rows = build_network("krf", AR_BUCKETS)
    for r in rows:
        if r.tier == "EST":
            up = _bracket_rate(r.rate, 3.0)
            dn = _bracket_rate(r.rate, 1.0 / 3.0)
            assert up[0] == r.rate[0] == dn[0], r.rid
    active, bound = est_active_and_bound(
        rows, {"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
    assert {r.rid for r in active} == EST_ACTIVE_AR_FILL
    assert all("Ne" in str(r.reactants) or "NeF*" in r.reactants
               for r in bound)
    by_rid = {r.rid: r for r in rows}
    for rid in CONFLICT_ALTERNATES:
        assert rid in by_rid, rid


def test_est_bracket_spot_measurement():
    """One live bracket point (the URAD_Bv lever, x3): the sweep really
    moves the observable, in the physically-required direction (a faster
    hot-manifold radiative drain lowers the relaxed-B population and
    with it peak g0). Full sweep: run.py --est-sweep (measured verdict
    recorded in CHANGELOG 0.1.9)."""
    from REBUILD.MODULES.scorecard import est_bracket_rows
    rep = est_bracket_rows(t_end_ns=60.0, only_rids={"URAD_Bv"})
    assert set(rep["active_rids"]) == EST_ACTIVE_AR_FILL
    assert len(rep["bound_rids"]) == 15
    d = rep["swept"]["URAD_Bv"][3.0]
    assert d["dg0_rel"] < 0.0
    assert abs(d["dg0_rel"]) < 0.6
    from REBUILD.MODULES.chemistry import CompiledNetwork
    from REBUILD.MODULES.species import KRF_MODE_SPECIES, StateLayout
    from REBUILD.MODULES.loki_tier import rates_for_mode
    import numpy as np
    lok = rates_for_mode("krf", {"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
    rows = build_network("krf", lok.available_buckets())
    lay = StateLayout(KRF_MODE_SPECIES, 1)
    cn = CompiledNetwork(rows, lay, lok, mode_optics("krf"))
    te = np.array([1.0])
    for r in rows:
        k = cn._k(r.rate, te, 7.4e25)
        val = float(k if np.isscalar(k) else k[0])
        assert val >= 0.0 and val == val, r.rid       # finite, non-negative


def test_mode_optics_complete():
    for mode in ("krf", "arf_chen"):
        o = mode_optics(mode)
        assert o["sigma_se_m2"] > 0 and o["tau_upper_s"] > 0
        assert o["upper"] and o["absorbers"] and o["fragments"]


def test_krf_buffer_partner_families_complete():
    """0.1.38 teeth: every catalytic-partner family that feeds or relaxes the
    lasing state must carry ALL three buffer partners (Ar, Ne, Kr). The 0.1.4
    split shipped B<->C mixing and Bv->B relaxation with Ar and Ne partners
    only, so a fill with little or no Ar/Ne (Electra's 40 % Kr, Mandl 1985's
    99.6 % Kr) had no channel into KrF*(B) — a pure Kr/F2 fill could not lase
    at all. Construct the defect (drop the Kr rows) and assert the check
    catches it; then assert the shipped network passes."""
    rows = build_network("krf", AR_BUCKETS)
    families = {
        "B->C mixing": (("KrF*(B)",), ("KrF*(C)",)),
        "C->B mixing": (("KrF*(C)",), ("KrF*(B)",)),
        "Bv->B relaxation": (("KrF*(Bv)",), ("KrF*(B)",)),
        "Ar** p->s relaxation": (("Ar**",), ("Ar*",)),
        "Kr** p->s relaxation": (("Kr**",), ("Kr*",)),
    }

    def partners(rows_, rea, pro):
        found = set()
        for r in rows_:
            for m in ("Ar", "Ne", "Kr"):
                if (tuple(sorted(r.reactants)) == tuple(sorted(rea + (m,)))
                        and tuple(sorted(r.products)) == tuple(sorted(pro + (m,)))):
                    found.add(m)
        return found

    for name, (rea, pro) in families.items():
        assert partners(rows, rea, pro) == {"Ar", "Ne", "Kr"}, name
    # detailed balance holds for every partner: k(C->B)/k(B->C) = 0.682
    for m in ("Ar", "Ne", "Kr"):
        bc = next(r for r in rows if r.rid == f"BC_{m}")
        cb = next(r for r in rows if r.rid == f"CB_{m}")
        assert abs(cb.rate[1] / bc.rate[1] - 0.682) < 1e-3, m
    # teeth: the pre-0.1.38 network (no Kr partners) must FAIL this gate
    broken = [r for r in rows if r.rid not in ("BC_Kr", "CB_Kr", "VREL_Kr")]
    assert partners(broken, ("KrF*(Bv)",), ("KrF*(B)",)) == {"Ar", "Ne"}
    # and Kr2* formation exists with Kr as the third body (AK8_Kr)
    assert any(tuple(sorted(r.reactants)) == ("Kr", "Kr", "Kr*")
               and "Kr2*" in r.products for r in rows), "Kr*+Kr+Kr -> Kr2*+Kr missing"
