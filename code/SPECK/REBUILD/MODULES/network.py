"""The reaction network AS DATA — tiers, citations, mechanical balance.

Defect class D2's prevention (ORESTES_CLASS_REBUILD.txt §2): "the network is
DATA; the RHS is assembled from it. Then every reaction is individually
testable, mass/charge balance is checkable mechanically per row, and no term
can hide." Every row here carries (tier, cite); rate laws are DECLARATIVE
tuples evaluated by ``chemistry.py`` (no lambdas — a spec can be printed,
diffed, and audited; a closure cannot).

Rate spec vocabulary (SI: m^3/s two-body, m^6/s three-body, s^-1 radiative):

  ("const", k)                  constant
  ("te_fit", A, B, C)           A * Te^B * exp(C/Te)   [Chen Table III form]
  ("stepwise", E_bind_eV)       Chen R9 template te_fit(1.56e-13, 0.71, -E)
  ("dr_sqrt", k0)               k0 / sqrt(max(Te, 0.02))   [molecular-ion DR]
  ("flannery", a, scale_cm3)    Kannari/Petrova pressure-turnover ion-ion:
                                [a*N/(N^2+1.7)] * scale_cm3 * 1e-6, N amagat
  ("lee_ii", A_cm3, b, c)       Lee 1989 [49] ion-ion form:
                                A*N/(N^2+b*N+c) * 1e-6, N amagat
  ("loki", bucket)              electron-impact k from the LOKI/KINETIC tier

ENERGY ROUTING IS MECHANICAL (see species.py): each row's release_eV =
sum(E_form reactants) - sum(E_form products). Endothermic rows draw from the
electron gas when an electron is a reactant, else from the gas. Exothermic
release goes to radiation (kind='rad', quantum ``hv_eV``), to the electron
gas when the row ejects a new free electron from heavy reactants (Penning
class), else to the gas. ``e_cost_override_eV`` exists for the dissociation
rows whose activation exceeds the potential-energy change (Chen R13/R14) —
the excess is released to the fragments, so closure is unaffected.

THE PLAN (72) CLOSURE: the Kannari p→s manifold relaxation rows (KPS*/APS*)
are IN this network unconditionally. No flag gates them; the class of defect
"physics default chosen by a parity requirement" cannot arise because there
is no parity obligation and no flag.

Row provenance: CHEN rows are chen2025.TABLE_III verbatim; REPO/LIT/EST rows
are the spectre network (projects/spectre modules/spectre/rates.py, audited
2026-07-14..15) transcribed to declarative specs with their citation strings
intact — including the AK7/AK8 Ar-third-body rows whose ABSENCE once left
Kr2+ formation dead in Ar fills (spectre CHANGELOG 2026-07-15, gain-gap
closure), and the K68Ne/NK5 first-hand Xu/Gadomski/Setser fetches.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from .chen2025 import CHEN_ERRATA, TABLE_III
from .provenance import check_tier
from .species import (ELECTRON_NAME, PHOTON_PREFIX, atoms_of, charge_of,
                      e_form_of)

_RATE_KINDS = ("const", "te_fit", "stepwise", "dr_sqrt", "flannery",
               "lee_ii", "loki")


@dataclass(frozen=True)
class Reaction:
    rid: str
    reactants: tuple
    products: tuple
    rate: tuple                     # declarative spec, see module docstring
    tier: str
    cite: str
    kind: str = "k"                 # 'k' collisional | 'rad' s^-1
    tg_pow: float = 0.0             # k *= (Tg/300)^tg_pow (Boichenko -0.75 set)
    # (0.1.40) T-independent share of the 300 K rate: with tg_pow set,
    # k *= tg_floor + (1 - tg_floor)(Tg/300)^tg_pow — the a + b Tg^n form
    # NRL prints for ArF*+2Ar (petrova2022 Table A5), normalised at 300 K.
    # 0.0 (the default) is the pure power law, bit-identical to pre-0.1.40.
    tg_floor: float = 0.0
    e_cost_override_eV: float = None
    energy_to: str = None           # override: 'gas' | 'radiation' | 'electrons'
    hv_eV: float = None             # rad rows: the emitted quantum

    # ------------------------------------------------------------ derived
    def release_eV(self) -> float:
        """sum(E_form reactants) - sum(E_form products) [eV]; >0 exothermic."""
        return (sum(e_form_of(s) for s in self.reactants)
                - sum(e_form_of(s) for s in self.products))

    def ne_delta(self) -> int:
        """Free electrons created (products) minus consumed (reactants)."""
        return (sum(1 for s in self.products if s == ELECTRON_NAME)
                - sum(1 for s in self.reactants if s == ELECTRON_NAME))

    def heavy_reactants(self) -> tuple:
        return tuple(s for s in self.reactants if s != ELECTRON_NAME
                     and not s.startswith(PHOTON_PREFIX))

    def is_electron_row(self) -> bool:
        return ELECTRON_NAME in self.reactants


class NetworkError(ValueError):
    """A row failed a mechanical check."""


# The ONLY channel-collision exemptions (see check_network): physically
# distinct channels whose element signatures collapse. Each entry must
# say why.
ALLOWED_SAME_CHANNEL = frozenset({
    # R11/R12: excitation of F2's REPULSIVE a3Piu / A1Piu states — two
    # measured LOKI rates at two thresholds (3.16 / 4.34 eV). After the
    # 0.1.14 repulsive reroute both deliver e+F2 -> F+F+e; they remain
    # two physical channels (per electronic state), not a double count.
    frozenset({"R11", "R12"}),
})


def check_reaction(r: Reaction) -> None:
    """Mechanical per-row audit: balance, provenance, spec validity."""
    check_tier(r.tier)
    if not r.cite:
        raise NetworkError(f"{r.rid}: no citation")
    if not r.rate or r.rate[0] not in _RATE_KINDS:
        raise NetworkError(f"{r.rid}: unknown rate spec {r.rate!r}")
    if r.kind not in ("k", "rad"):
        raise NetworkError(f"{r.rid}: unknown kind {r.kind!r}")
    if r.kind == "rad":
        if len(r.heavy_reactants()) != 1 or r.rate[0] != "const":
            raise NetworkError(f"{r.rid}: a radiative row is one species with "
                               f"a constant s^-1 rate")
    # --- elemental mass balance (the r115-class killer) ---------------------
    lhs, rhs = {}, {}
    for side, acc in ((r.reactants, lhs), (r.products, rhs)):
        for s in side:
            for el, n in atoms_of(s).items():
                acc[el] = acc.get(el, 0) + n
    if lhs != rhs:
        raise NetworkError(f"{r.rid}: element balance fails: {lhs} -> {rhs}")
    # --- charge balance -----------------------------------------------------
    qL = sum(charge_of(s) for s in r.reactants)
    qR = sum(charge_of(s) for s in r.products)
    if qL != qR:
        raise NetworkError(f"{r.rid}: charge balance fails: {qL} -> {qR}")
    # --- energy routing sanity ---------------------------------------------
    if r.hv_eV is not None:
        if r.kind != "rad":
            raise NetworkError(f"{r.rid}: hv_eV on a non-radiative row")
        if r.hv_eV > r.release_eV() + 1e-9:
            raise NetworkError(
                f"{r.rid}: emitted quantum {r.hv_eV} eV exceeds the row's "
                f"release {r.release_eV():.3f} eV — radiating energy that "
                f"was never stored")
    if r.energy_to not in (None, "gas", "radiation", "electrons"):
        raise NetworkError(f"{r.rid}: unknown energy_to {r.energy_to!r}")
    # --- gas-temperature form (0.1.40) -------------------------------------
    if r.tg_floor != 0.0:
        if r.tg_pow == 0.0:
            raise NetworkError(f"{r.rid}: tg_floor without tg_pow — the "
                               "floor is a share of a Tg-dependent rate")
        if not 0.0 < r.tg_floor < 1.0:
            raise NetworkError(f"{r.rid}: tg_floor {r.tg_floor} outside "
                               "(0, 1) — it is the T-independent SHARE of "
                               "the 300 K rate")


def check_network(rows) -> None:
    seen = set()
    channels = {}
    for r in rows:
        if r.rid in seen:
            raise NetworkError(f"duplicate rid {r.rid}")
        seen.add(r.rid)
        check_reaction(r)
        # --- channel uniqueness (the 0.1.10 double-count killer): two rows
        #     with the same reactant AND product multisets are one physical
        #     channel counted twice. Distinct product branches of one
        #     reactant pair remain legal (AKD3/AKD5-class branching).
        #     Found 2026-08-05: 8 channels double-counted because row
        #     existence was checked by grepping SOURCE while the CHEN block
        #     is GENERATED from data — this makes that class impossible.
        #     ALLOWED_SAME_CHANNEL lists the sole exemptions: rid pairs
        #     that are PHYSICALLY DISTINCT channels whose element-level
        #     signatures collapse (R11/R12: two repulsive F2 states, two
        #     measured rates and thresholds, same fragments after the
        #     0.1.14 reroute).
        key = (tuple(sorted(r.reactants)), tuple(sorted(r.products)))
        if key in channels:
            pair = frozenset((r.rid, channels[key]))
            if pair not in ALLOWED_SAME_CHANNEL:
                raise NetworkError(
                    f"{r.rid}: duplicate channel with {channels[key]}: "
                    f"{key[0]} -> {key[1]} — one physical channel, two rows")
        channels[key] = r.rid


R = Reaction

# --- excimer three-body temperature form (0.1.40) ---------------------------
# Every excimer three-body association row (KrF* + 2Ar / Kr+Ar / 2Kr /
# Kr+Ne, ArF* + 2Ar) is a 300 K value whose SOURCE model carried a steep
# negative gas-temperature dependence that this engine dropped on adoption:
#   - Kannari, Obara & Fujioka, JAP 57, 4309 (1985) Table I footnote b:
#     "Expressed for 300 K gas temperature. Temperature dependence is
#     available from the calculation method in Refs. 49-51" (= V. H. Shui,
#     APL 31, 50 (1977); APL 34, 203 (1979); Shui & Duzy, APL 36, 135
#     (1980)), and Sec. II.D: the three-body rates were "treated here as a
#     temperature dependent rate, which was based on the theoretical
#     calculation by Shui et al.";
#   - Lee, Matsui, Kannari & Obara, IEEE TED 36, 2053 (1989) Table I
#     footnote b ("At a gas temperature of 300 K") and text: the ArF*+2Ar
#     quench "is reduced with increasing excitation rate due to the gas
#     temperature rise";
#   - NRL Orestes 2022 (petrova2022 Table A5, red), the one PRINTED ArF
#     form: ArF* + 2Ar -> Ar2F* + Ar = 1.62e-32 + 2.245e-21 Tg^-3.9 cm6/s
#     [k28, Shui-1979, Chen-1978] (all three read first-hand 2026-09-18).
# KrF: Shui's calculation itself, in the Avco contract report that carries
# it (Mangano et al., "One Meter KrF Laser System", DTIC ADA069053, Sec.
# V.B, p. 73, page image read 2026-09-18): "a change in temperature from
# 300°K to 600°K lowers kr(Ar) by a factor of 8 and kr(Kr) by a factor of
# 10" (R = the atom that binds: KrF*+2Ar -> ArKrF*, KrF*+Kr+M -> Kr2F*),
# i.e. (Tg/300)^-3.00 and ^-3.32; Fig. 22 on the page agrees (5.1e-31 ->
# 5.3e-32). MEASURED: Klimek & Hsia, JAP 52, 5361 (1981), 315-420 K at
# 1-3 amagat, R, M = Ar, Kr: "agrees to within 15% of that calculated by
# using classical phase-space theory". Applied as pure power laws on each
# row's own 300 K value; the Kr-bound exponent covers K34Ar (Hsia/Jacob/
# Duzy take k_ArKr ~ k_2Kr) and K68Ne (Ne third body, by the same R rule).
KRF_3B_TG_POW_AR = math.log(1.0 / 8.0) / math.log(2.0)      # -3.000
KRF_3B_TG_POW_KR = math.log(1.0 / 10.0) / math.log(2.0)     # -3.322
# ArF: NRL's form normalised at 300 K on the row's own 300 K value,
# k(Tg) = k300 [f + (1 - f)(Tg/300)^-3.9], f = 1.62e-32 / k_NRL(300) = 0.0320
ARF_3B_TG_POW = -3.9
ARF_3B_TG_FLOOR = 1.62e-32 / (1.62e-32 + 2.245e-21 * 300.0 ** -3.9)

# ---------------------------------------------------------------------------
# CHEN tier — chen2025.TABLE_III mapped to rows (bolsig -> LOKI buckets).
# ---------------------------------------------------------------------------
_LOKI_BUCKET = {"R1": "Ne*", "R2": "Ne**", "R3": "Ne+", "R10": "F-",
                "R11": "F2(a3)", "R12": "F2(A1)", "R15": "F2+",
                "R5": "Ar*", "R6": "Ar**", "R7": "Ar+"}

# radiative quanta for the Chen rad rows [eV]
_CHEN_HV = {"R37": 6.42, "R39": 13.05}   # 193 nm; Ar** ~95 nm resonance
# activation overrides (fit exponent C is the activation the electron pays;
# the excess over the potential-energy change is released to the fragments)
_E_COST = {"R13": 4.69, "R14": 4.94}


def chen_rows(available_buckets, errata: bool = True) -> list:
    """Chen Table III as Reactions. 'bolsig' rows map to tier buckets when
    the loaded electron-rate tables provide them (spectre's rule: Ar ground
    rows live only when an Ar-mixture table is present); R4 uses the cited
    Lee 1989 fit; R8 the R9-template EST fallback. R16/R17 (trace Xe,
    1e-4 mole fraction) are omitted exactly as spectre omits them — recorded
    in the scorecard as an omission, not silently. 'sigma' rows (R38/R40/R41)
    are optics data, not reactions — see mode_optics().

    ``errata`` (0.1.25, default True) applies chen2025.CHEN_ERRATA — the
    adjudicated corrections to the five NUMERIC transcription errors of the
    printed table (R13/R14 cm3-misprint, R23 R27-paste, R24 x1e6 slip, R36
    exponent typo; forensics on the chen2025 module docstring). The errata
    lane reproduces Chen's MODEL (the printed values cannot describe what
    he ran — measured: they sever the Ne bath's Penning chain and NaN the
    integration under his own base-case circuit); ``errata=False`` runs the
    table AS PRINTED, kept as the record of the PAPER. The R19 exponent and
    R26 product adjudications predate the errata block and are
    unconditional, as before."""
    avail = set(available_buckets)
    out = []
    for rid, rea, pro, kind, val, ref, note in TABLE_III:
        cite = f"Chen 2025 Table III {rid} [ref {ref}]" + (
            f"; {note}" if note else "")
        fix = CHEN_ERRATA.get(rid) if errata else None
        if fix is not None:
            cite += (f" ; ERRATA 0.1.25: printed {fix[1]:g} -> {fix[0]:g} "
                     f"({fix[2]})")
        if kind == "bolsig":
            b = _LOKI_BUCKET.get(rid)
            if b is not None and b in avail:
                out.append(R(rid, rea, pro, ("loki", b), "LOKI",
                             cite + " ; k(<eps>) from the curated LoKI-MC "
                             "mixture tables (PLAN (45))"))
            elif rid == "R4":
                out.append(R(rid, rea, pro, ("te_fit", 3.04e-14, 0.54, -7.50),
                             "REPO", cite + " ; Lee 1989 Table I [28]: "
                             "3.04e-8 Te^0.54 exp(-7.50/Te) cm3/s"))
            elif rid == "R8":
                out.append(R(rid, rea, pro, ("stepwise", 4.21), "EST",
                             cite + " ; R9-template, E_bind(Ar*) = 4.21 eV"))
            continue
        if rid == "R19":
            # exponent triple-confirmed NEGATIVE (Ma 2024 #46, Lee 1989 [45],
            # Boichenko 2000 #13); the printed Te^+0.5 is a typography flag
            # preserved in chen2025.py.
            out.append(R(rid, rea, pro, ("te_fit", 3.7e-14, -0.43, 0.0),
                         "REPO", cite + " ; exponent -0.43 (Ma #46, 3-source "
                         "adjudication, spectre 2026-07-14)"))
            continue
        if kind == "te_fit":
            A, B, C = val
            if fix is not None:
                A = fix[0]                 # prefactor only; exponents stand
            out.append(R(rid, rea, pro, ("te_fit", A, B, C), "CHEN", cite,
                         e_cost_override_eV=_E_COST.get(rid)))
        elif kind in ("k2", "k3"):
            tgp = -0.75 if rid in ("R23", "R24", "R26", "R27", "R28") else 0.0
            out.append(R(rid, rea, pro,
                         ("const", fix[0] if fix is not None else val),
                         "CHEN", cite, tg_pow=tgp))
        elif kind == "rad":
            # strip the photon markers: the quantum is carried by hv_eV
            pro_clean = tuple(p for p in pro if not p.startswith(PHOTON_PREFIX))
            out.append(R(rid, rea, pro_clean, ("const", val), "CHEN", cite,
                         kind="rad", hv_eV=_CHEN_HV.get(rid)))
        # 'sigma' rows: mode_optics()
    return out


# ---------------------------------------------------------------------------
# REPO/LIT/EST tiers — the spectre network, transcribed to declarative specs
# with citation strings intact (rates.py functions kr_core/ne_kr_couplings/
# ne_refinements/nef_channels/ar_kr_couplings, audited 2026-07-14..15).
# ---------------------------------------------------------------------------
def kr_core() -> list:
    return [
        R("K25", ("Kr*", "F2"), ("KrF*", "F"), ("const", 8.1e-16), "REPO",
          "J&H 1980 Table I row 1; VERIFIED FIRST-HAND in Kannari 1985 "
          "Table I p.4311: 8.1e-10 [refs 52,59,93] (2026-08-04); "
          "config.py:261"),
        R("K27", ("Kr**", "F2"), ("KrF*", "F"), ("const", 8.1e-16), "REPO",
          "J&H 1980 Table I row 6; config.py:262 K27"),
        R("K26", ("Kr+", "F-"), ("KrF*",), ("flannery", 2.4, 3.417e-6),
          "REPO", "Flannery turnover, Kannari 1985 Sec II.B / Petrova 2022 "
          "A5; config.py:496 k26_flannery"),
        R("K108", ("Kr2+", "F-"), ("KrF*", "Kr"), ("flannery", 2.65, 3.096e-6),
          "REPO", "same Flannery form; config.py:499 k108_flannery (flat "
          "1e-6 placeholder rejected)"),
        R("K15", ("Kr+", "Kr", "Kr"), ("Kr2+", "Kr"), ("const", 2.5e-43),
          "LIT", "Kr+ +2Kr -> Kr2+ +Kr 2.5e-31 cm6/s by THREE first-hand "
          "sources (2026-08-04): J&H 1980 Table I #41 [30], Kannari 1985 "
          "CT block, Boichenko 2000 Table 1 #18 (1.62e-32 Tg[eV]^-0.75 = "
          "2.5e-31 at 300 K). Chow's 2.7e-31 (phantom config.py:263) sits "
          "in-spread and is superseded by the consensus",
          tg_pow=-0.75),
        R("R11p", ("e", "Kr2+"), ("Kr*", "Kr"), ("dr_sqrt", 2.24e-13),
          "LIT", "Kannari 1985 Table I p.4312 DR block (first-hand "
          "2026-08-05): Kr2+ +e -> Kr*+Kr 1.2e-6 (Te/300K)^-0.55 [71] "
          "+ -> Kr**+Kr 1.9e-7 [68-70], lumped to the dominant printed "
          "branch at the 300 K total 1.39e-6 cm3/s on the Te^-1/2 "
          "template (the ADR convention; 1.39e-12/6.219 = 2.24e-13). "
          "Supersedes the REPO reactions.py:65 2.4e-7/sqrt(Te) -> "
          "Kr**+Kr (x1.07 above this mapping with the product on the "
          "SUB-dominant branch; the 0.1.0 citation debt is now paid). "
          "RECORDED: J&H #60 prints a single Kr**-product row, typical "
          "1.1e-7 — the branch question; Kr** relaxes to Kr* via "
          "KPS_Ar-class rows either way"),
        R("K67", ("KrF*", "e"), ("Kr", "F", "e"), ("const", 4.5e-14), "REPO",
          "Kushner-Hanson-Schneider 1989 APL 55,2482; config.py:298"),
        R("K32", ("KrF*", "F2"), ("Kr", "F", "F", "F"), ("const", 3.0e-16),
          "LIT", "Xu/Gadomski/Setser JCP 99,2591 (1993) Table I + p.2594: "
          "KrF(B,C)+F2 = (3.0+-0.3)e-10 cm3/s, 'the average of all three "
          "determinations ... We strongly support this smaller value' "
          "(VERIFIED FIRST-HAND 2026-08-04; Hakuta's 2.5+-1.0e-10 "
          "agrees). Superseded era: Kannari prints 7.8e-10; J&H #12-era "
          "5-6e-10 class — both recorded"),
        R("K33", ("KrF*", "Kr"), ("Kr", "Kr", "F"), ("const", 4.0e-18),
          "LIT", "Xu/Gadomski/Setser JCP 99,2591 (1993) Table I + p.2596: "
          "KrF(B,C)+Kr recommended (4.0+-0.5)e-12 cm3/s — the dedicated "
          "later MEASUREMENT (read first-hand 2026-08-04), adopted per "
          "the AK1 precedent. Source history: Kannari 1985 prints 2.0e-12 "
          "marked 'estimated'; phantom's 1.0e-12 was attributed to 'J&H "
          "Table 1' but J&H #12 actually prints 2.0e-11 (a x20 mis-cite, "
          "found 2026-08-04). Spread [2.0e-12 .. 2.0e-11] recorded"),
        R("K34", ("KrF*", "Kr", "Kr"), ("Kr2F*", "Kr"), ("const", 6.0e-43),
          "LIT", "Xu/Gadomski/Setser JCP 99,2591 (1993) Table II + "
          "abstract: KrF(B,C)+2Kr = (6.0+-2.0)e-31 cm6/s (VERIFIED "
          "FIRST-HAND 2026-08-04; was a repo fetch). Kannari 1985 "
          "Table I prints 6.7e-31 [62,63] — in-band. 300 K value; Tg "
          "form (0.1.40): Shui's R = Kr calculation, x1/10 from 300 to "
          "600 K, measured to 15 % by Klimek & Hsia 1981 "
          "(KRF_3B_TG_POW_KR)",
          tg_pow=KRF_3B_TG_POW_KR),
        R("K68Ne", ("KrF*", "Kr", "Ne"), ("Kr2F*", "Ne"), ("const", 4.9e-43),
          "LIT", "Xu/Gadomski/Setser JCP 99,2591 (1993) Table II: Ne+Kr "
          "mixed 3-body = (4.9+-1.5)e-31 cm6/s (VERIFIED FIRST-HAND "
          "2026-08-04; was a repo fetch). 300 K value (Xu measured at "
          "300 K); Tg form (0.1.40) by Shui's R rule: the bound atom is "
          "Kr, so the R = Kr exponent (KRF_3B_TG_POW_KR); no Ne-third-"
          "body temperature study exists in the tree — stated EST",
          tg_pow=KRF_3B_TG_POW_KR),
        # FOUND MISSING 2026-08-04: the Ar-third-body KrF* three-body
        # quench — dominant in Ar fills (Ar is ~95% of the DE-AR mixture;
        # ~1.7e8 s^-1 at the v2 point, comparable to the B radiative
        # rate) — was absent while the Kr- and Ne-partner rows existed.
        R("K34Ar", ("KrF*", "Kr", "Ar"), ("Kr2F*", "Ar"), ("const", 6.5e-43),
          "LIT", "TRIPLE-SOURCED, all first-hand 2026-08-04: J&H 1980 "
          "Table I #18 6.5e-31 cm6/s [their 13]; Kannari 1985 Table I "
          "6.5e-31 [62,63]; Xu 1993 Table II mixed Ar+Kr = (5.5+-2.0)e-31 "
          "— 6.5 sits inside Xu's band, so the recorded follow-up "
          "(Xu re-verification) is CLOSED as confirmed. Kannari's row "
          "carries footnote b ('Expressed for 300 K gas temperature. "
          "Temperature dependence is available from ... Refs. 49-51' = "
          "Shui): the Tg form, dropped until 0.1.40, is restored with "
          "Shui's R = Kr exponent (KRF_3B_TG_POW_KR; Hsia/Jacob/Duzy "
          "take k_ArKr ~ k_2Kr)",
          tg_pow=KRF_3B_TG_POW_KR),
        R("K169", ("Kr2F*",), ("Kr", "Kr", "F"), ("const", 4.5e6), "LIT",
          "Xu/Gadomski/Setser JCP 99,2591 (1993) p.2597: tau(Kr2F*) "
          "'established as 220+-20 ns to within 10%' by their dedicated "
          "re-measurement -> 4.5e6 s^-1 (read first-hand 2026-08-04; "
          "adopted per the AK1 later-dedicated-measurement precedent). "
          "Source history: Quigley-era 181+-12 ns = the 5.6e6 Kannari "
          "prints [90,93] (verified first-hand; phantom's 5.4e6 was the "
          "Ar2F* transposition); Ueda 210 ns; earlier Setser-lab "
          "240+-35 ns; J&H 1980 #68 prints 6.7e6. Spread recorded",
          kind="rad", hv_eV=3.0),
        R("K125", ("Kr2F*", "F2"), ("Kr", "Kr", "F", "F2"), ("const", 1.5e-16),
          "LIT", "Kannari 1985 Table I p.4311: Kr2F*+F2 -> 2Kr+F+F2 "
          "1.5e-10 cm3/s [refs 89,94] (verified first-hand 2026-08-04 — "
          "phantom's 2.0e-10, config.py:695, is the Ar2F*+F2 rate: "
          "transposed). Xu 1993 Fig.9/p.2597 re-measures (1.2+-0.2)e-10 "
          "(first-hand 2026-08-04) — in-spread, kept at the Kannari/"
          "Morgan-Szoke value; J&H #22's 1.0e-9 is the superseded era"),
        R("K_2Fe", ("Kr2F*", "e"), ("Kr", "Kr", "F", "e"), ("const", 1.0e-13),
          "LIT", "Kannari 1985 Table I p.4312 superelastic: Kr2F*+e "
          "1.0e-7 cm3/s, marked 'estimated' there (verified first-hand "
          "2026-08-04; config.py:695)"),
        R("KSTEP", ("e", "Kr*"), ("Kr+", "e", "e"), ("stepwise", 4.09), "EST",
          "R9-template, E_bind(Kr*) = 4.09 eV. Corroborated 2026-08-04: "
          "J&H 1980 Table I #48 gives the discharge-typical value 1.5e-8 "
          "cm3/s (2 kV/cm/atm); the template evaluates 1.4e-8 at "
          "Te = 1.5 eV — same class"),
    ]


def kannari_ps_relaxation() -> list:
    """THE PLAN (72) ROWS — p→s manifold relaxation, unconditionally on.

    phantom11 gates these behind LIT_QUENCH, which defaults False because 0.0
    reproduces og-pyfiles bit-for-bit — "a parity artifact, not a physics
    choice" (root CHANGELOG Known-bad). At 1 atm the 1e-10 cm3/s class rate
    gives tau ~ 0.45 ns, an order of magnitude faster than the KrF* radiative
    lifetime, so holding the manifolds isolated is not defensible at these
    pressures. Here the rows exist with NO flag; their uncertainty is
    expressed the honest way, as tier + bracket, not as an off switch.
    Buffer-partner analogs (Ne/Kr third bodies) are EST at the same class
    value pending a direct measurement.
    """
    cite = ("Kannari 1985 Table I p.4311: Kr**+Ar -> Kr*+Ar and "
            "Ar**+Ar -> Ar*+Ar, both 1.0e-10 cm3/s and both marked "
            "'estimated' THERE (VERIFIED FIRST-HAND 2026-08-04); the "
            "phantom11 K_PS_RELAX value, PLAN (72) root-cause entry")
    return [
        R("KPS_Ar", ("Kr**", "Ar"), ("Kr*", "Ar"), ("const", 1.0e-16),
          "LIT", cite),
        R("KPS_Ne", ("Kr**", "Ne"), ("Kr*", "Ne"), ("const", 1.0e-16),
          "EST", cite + "; Ne-partner analog at the class value"),
        R("KPS_Kr", ("Kr**", "Kr"), ("Kr*", "Kr"), ("const", 1.0e-16),
          "EST", cite + "; Kr-partner analog at the class value"),
        R("APS_Ar", ("Ar**", "Ar"), ("Ar*", "Ar"), ("const", 1.0e-16),
          "LIT", cite + "; independently Lee 1989 Table I [36] 1.0e-10 "
          "(second source, first-hand 2026-08-04). CHAIN-2 HEAD "
          "(0.1.14): NRL's Orestes carries NO neutral 4p->4s quench at "
          "all (petrova2022 A3 relaxes Ar** by electrons only) — a "
          "STRUCTURAL disagreement worth x15 on Ar** and with it the "
          "193-nm absorber budget; the two-print value is KEPT and the "
          "x0.1 alternate is SWEPT every est-sweep "
          "(CONFLICT_ALTERNATES)"),
        R("APS_Ne", ("Ar**", "Ne"), ("Ar*", "Ne"), ("const", 1.0e-16),
          "LIT", "Lee 1989 Table I p.2055 [37]: Ar** +Ne -> Ar* +Ne "
          "1.0e-10 cm3/s — the ACTUAL process, first-hand 2026-08-04"),
        R("APS_Kr", ("Ar**", "Kr"), ("Ar*", "Kr"), ("const", 1.0e-16),
          "EST", cite + "; Kr-partner analog at the class value (0.1.38: "
          "the Kr-rich fills had Ar** relaxation on Ar/Ne only)"),
    ]


def ne_kr_couplings() -> list:
    return [
        R("NK1", ("Ne*", "Kr"), ("Kr+", "Ne", "e"), ("const", 6.5e-17), "EST",
          "Ne*+Ar Penning analog 6.5e-11 cm3/s: Kannari 1985 [Chen R21] = "
          "Lee 1989 [38] = Ma 2024 #33 (Lee [38] digits verified "
          "first-hand 2026-08-04; the Kr process itself is unmeasured)"),
        R("NK1b", ("Ne*", "Kr"), ("NeKr+", "e"), ("const", 2.3e-17), "EST",
          "associative-ionization branch, Ne*+Ar->NeAr+ +e analog 2.3e-11: "
          "Lee 1989 [38] (digits verified first-hand 2026-08-04)"),
        R("NK2", ("Ne2+", "Kr"), ("Kr+", "Ne", "Ne"), ("const", 5e-19), "LIT",
          "Boichenko/Yakovlenko QE 30,681 p.682: Ne2+ +Kr <= 5e-13 cm3/s "
          "[their 12]; value AT the published bound (RE-VERIFIED first-hand "
          "2026-08-04). Recorded conflict: Lee 1989 Table I [43] lists the "
          "Ar twin Ne2+ +Ar at 3.0e-11 — far above Boichenko's Ar bound; "
          "the Kr bound is the direct source and is kept"),
        R("NK3", ("Ne+", "Kr"), ("Kr+", "Ne"), ("const", 1e-17), "EST",
          "slow atomic CT: Ne++Ar and Ne++Xe both 1.0e-11 cm3/s "
          "(Ma #75/#76; Chen R25)"),
        R("NK4", ("Kr+", "Kr", "Ne"), ("Kr2+", "Ne"), ("const", 2.5e-43),
          "EST", "Ne-third-body dimerization, adjudicated first-hand "
          "2026-07-15 (Boichenko #39 vs measured-analog cluster); keep "
          "2.5e-31 cm6/s, bracket [1.0, 3.0]e-31", tg_pow=-0.75),
        R("NK5", ("KrF*", "Ne"), ("Kr", "Ne", "F"), ("const", 2.1e-18), "LIT",
          "Xu/Gadomski/Setser JCP 99,2591 (1993) Table I: KrF(B,C)+Ne = "
          "(2.1+-0.5)e-12 cm3/s (VERIFIED FIRST-HAND 2026-08-04, with "
          "their p.2596 caveat now recorded: rare-gas two-body values "
          "are impurity-limited UPPER BOUNDS — value kept AT the bound, "
          "the NK2 precedent)"),
        R("NK6", ("Kr*", "Kr", "Ne"), ("Kr2*", "Ne"), ("const", 2.5e-43),
          "EST", "analog of Chen R23 (Ar*+Ar+Ne->Ar2*+Ne)", tg_pow=-0.75),
        R("AK8_Kr", ("Kr*", "Kr", "Kr"), ("Kr2*", "Kr"), ("const", 2.7e-44),
          "EST", "Kr-third-body analog of AK8 (Kannari 1985 Kr*+Kr+Ar -> "
          "Kr2*+Ar 2.7e-32 cm6/s) at the class value (0.1.38): the "
          "krypton-rich and buffer-free fills had NO Kr2* formation "
          "channel with Kr as the third body (K15 covers Kr2+, nothing "
          "covered Kr2*)"),
        R("NK7", ("Kr2*", "F2"), ("KrF*", "Kr", "F"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Kr2*+F2 -> KrF*+Kr+F "
          "3.0e-10 cm3/s [ref 61] (first-hand 2026-08-04; replaces the "
          "K25-class EST analog, which sat x2.7 high)"),
        R("NK8", ("Ne2*", "Kr"), ("Kr+", "Ne", "Ne", "e"), ("const", 6.5e-17),
          "EST", "Ne2*+Ar Penning analog 6.5e-11 (Ma #36, Lee [38])"),
        R("HK1", ("Ne+", "Kr", "Ne"), ("NeKr+", "Ne"), ("const", 1.0e-43),
          "LIT", "Boichenko 2000 Table 1 #23: Ne+ +Kr+Ne -> NeKr+ +Ne "
          "1e-31 cm6/s — the ACTUAL process, first-hand 2026-08-04 (their "
          "text marks the universal 1e-31 as assumed; Lee [40] NeAr+ "
          "analog 1.0e-31 verified first-hand corroborates)",
          tg_pow=-0.75),
        R("HK2", ("NeKr+", "Kr"), ("Kr+", "Ne", "Kr"), ("const", 1.1e-15),
          "EST", "NeAr++Ar->Ar++Ne+Ar analog 1.1e-9 cm3/s (Lee [44], "
          "digits verified first-hand 2026-08-04)"),
        R("HK3", ("NeKr+", "Kr"), ("Kr2+", "Ne"), ("const", 1.2e-17), "EST",
          "NeAr++Ar->Ar2++Ne analog 1.2e-11 cm3/s (Lee [44], digits "
          "verified first-hand 2026-08-04)"),
        R("HK4", ("e", "NeKr+"), ("Kr**", "Ne"), ("dr_sqrt", 1.0e-13), "LIT",
          "Boichenko 2000 Table 1 #8: NeKr+ +e -> Kr** +Ne, 1e-7 Te^-1/2 "
          "cm3/s — the ACTUAL process, RE-VERIFIED first-hand 2026-08-04. "
          "Recorded conflict: Lee [46]'s NeAr+ DR form (1.5e-9 (Tg/Te)^"
          "0.67) sits far lower; Boichenko's typical-magnitude direct row "
          "is kept"),
    ]


def ne_refinements() -> list:
    return [
        R("NR1", ("Ne2*",), ("Ne", "Ne"), ("const", 7.5e7), "REPO",
          "Ne2* radiative decay 7.5e7 s-1 (Rauf&Kushner 1999 Table I; "
          "Lee [47] 8.9e7)", kind="rad", hv_eV=16.1),
        R("NR2", ("e", "Ne**"), ("Ne+", "e", "e"),
          ("te_fit", 1.0e-13, 0.61, -4.65), "REPO",
          "Lee 1989 [28]: 1.0e-7 Te^0.61 exp(-4.65/Te) cm3/s"),
        R("NR3", ("Ne**", "Ne"), ("Ne*", "Ne"), ("const", 7.1e-17), "LIT",
          "Lee 1989 Table I p.2055 [35]: 7.1e-11 cm3/s (verified "
          "first-hand 2026-08-04)"),
        R("NR4", ("Ne2*", "Ne2*"), ("Ne2+", "Ne", "Ne", "e"),
          ("const", 5.0e-16), "REPO", "Ma 2024 Table 5 #29: 5.0e-16 m3/s "
          "dimer pooling"),
    ]


def nef_channels() -> list:
    return [
        R("NF1", ("Ne+", "F-"), ("NeF*",), ("lee_ii", 2.80e-5, 1.44, 10.9),
          "LIT", "Lee 1989 Table I p.2055 [49]: Ne+ +F- +Ne form "
          "(2.80e-5, 1.44, 10.9) — digits verified first-hand 2026-08-04"),
        R("NF2", ("Ne2+", "F-"), ("NeF*", "Ne"), ("lee_ii", 3.66e-5, 1.74, 17.3),
          "LIT", "Lee 1989 Table I p.2055 [49]: Ne2+ +F- +Ne form "
          "(3.66e-5, 1.74, 17.3) — digits verified first-hand 2026-08-04"),
        R("NF3", ("NeF*",), ("Ne", "F"), ("const", 1e9), "EST",
          "predissociation-dominated fate (Lee [62]: ~2/3 to Ne+F*); lumped "
          "fast decay — fragments, NOT a photon, hence energy_to='gas'",
          kind="rad", energy_to="gas"),
        R("HK5", ("NeKr+", "F-"), ("KrF*", "Ne"), ("lee_ii", 5.27e-5, 6.81, 20.3),
          "EST", "NeAr+ +F- +Ne -> ArF* +2Ne analog (Lee 1989 Table I "
          "[49] (5.27e-5, 6.81, 20.3) — analog digits verified first-hand "
          "2026-08-04; the NeKr+ process itself is unmeasured)"),
        # FOUND 2026-08-04, first-hand: Lee 1989 Table I p.2056 "Neutral
        # Formation" DOES list the Ne*/Ne** harpoon channels to NeF* —
        # CONTRADICTING the recorded spectre 2026-07-14 adjudication
        # ("the row does not exist on the page"), which shaped
        # nef_channels() to ion-ion only. At 0.1% F2 in a Ne fill this is
        # a ~1e7-1e8 s^-1 Ne*-drain competing with the Kr Penning chain.
        R("LNF1", ("Ne*", "F2"), ("NeF*", "F"), ("const", 5.3e-16), "LIT",
          "Lee 1989 Table I p.2056 [32]: Ne* +F2 -> NeF* +F 5.3e-10 cm3/s "
          "(first-hand 2026-08-04; corrects the spectre 2026-07-14 "
          "absence claim)"),
        R("LNF2", ("Ne**", "F2"), ("NeF*", "F"), ("const", 7.5e-16), "LIT",
          "Lee 1989 Table I p.2056 [54]: Ne** +F2 -> NeF* +F 7.5e-10 "
          "cm3/s (first-hand 2026-08-04)"),
    ]


def ar_kr_couplings() -> list:
    return [
        R("AK1", ("ArF*", "Kr"), ("KrF*", "Ar"), ("const", 1.6e-15), "LIT",
          "Kannari 1985 Table I p.4311: ArF*+Kr -> KrF*+Ar 1.6e-9 cm3/s "
          "[ref 62] (first-hand 2026-08-04). RECORDED CONFLICT: J&H 1980 "
          "Table I #2 prints 3.0e-10 [32,33] — a x5.3 inter-source "
          "spread; the later dedicated measurement [62] is kept"),
        R("AK2", ("Ar2+", "Kr"), ("Kr+", "Ar", "Ar"), ("const", 7.5e-16),
          "LIT", "Ar2+ +Kr -> Kr+ +2Ar 7.5e-10 cm3/s by TWO primary "
          "prints (both first-hand 2026-08-05): Kannari 1985 Table I "
          "p.4312 CT block [ref 64] and J&H 1980 #35 [31]. Boichenko & "
          "Yakovlenko QE 30,681 p.682's 'well-known' 7.8e-10 [their 12] "
          "sits x1.04 in-spread and is superseded by the primary "
          "consensus (the K15 precedent; was this row's sole source "
          "through 0.1.9)"),
        R("AK3", ("Ar+", "Kr"), ("Kr+", "Ar"), ("const", 3.0e-17), "LIT",
          "Kannari 1985 Table I p.4312: Ar+ +Kr -> Kr+ +Ar 3.0e-11 cm3/s "
          "[ref 64] (first-hand 2026-08-04; the EST 5e-10 guess was x17 "
          "high — non-resonant atomic CT is slow)"),
        R("AK4", ("Ar*", "Kr"), ("Kr*", "Ar"), ("const", 6.2e-18), "REPO",
          "excitation transfer 6.2e-12 cm3/s: J&H 1980 r132 (phantom "
          "config.py:540 K132); was a 2e-12 EST (x3.1 up, 2026-07-15)"),
        R("AK5", ("Ar*", "Ar", "Ar"), ("Ar2*", "Ar"), ("const", 1.14e-44),
          "LIT", "Kannari 1985 Table I p.4311: Ar*+2Ar -> Ar2*+Ar "
          "1.14e-32 cm6/s [ref 52] (verified first-hand 2026-08-04)"),
        R("AK6", ("Ar2*", "F2"), ("ArF*", "Ar", "F"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar2*+F2 -> ArF*+Ar+F "
          "3.0e-10 cm3/s [ref 61] (verified first-hand 2026-08-04; "
          "phantom config.py:436 K123B)"),
        R("AK7", ("Kr+", "Kr", "Ar"), ("Kr2+", "Ar"), ("const", 2.5e-43),
          "LIT", "Kr+ +Kr+Ar -> Kr2+ +Ar 2.5e-31 cm6/s by TWO primary "
          "prints (both first-hand 2026-08-05): Kannari 1985 Table I "
          "p.4312 CT block [refs 43,67] and J&H 1980 #43 [32,40] — the "
          "same consensus value as the 2Kr channel K15. Chow 1978's "
          "r60-Ar 3.25e-31 (refs [5,6]; phantom FIX-E1, this row's "
          "source through 0.1.9) sits x1.3 above and is superseded "
          "(K15 precedent). Tg^-0.75 kept per the K15/Boichenko "
          "ion-association family template. ABSENCE of this row once "
          "left Kr2+ formation dead in Ar fills (spectre CHANGELOG "
          "2026-07-15)", tg_pow=-0.75),
        R("AK8", ("Kr*", "Kr", "Ar"), ("Kr2*", "Ar"), ("const", 2.7e-44),
          "LIT", "Kannari 1985 Table I p.4311: Kr*+Kr+Ar -> Kr2*+Ar "
          "2.7e-32 cm6/s [ref 52] (first-hand 2026-08-04). RECORDED "
          "CONFLICT: the phantom REPO chain carried 1.0e-32 attributed "
          "to J&H r130 — x2.7 below the Kannari print; the first-hand "
          "later table is kept (AK1 precedent)", tg_pow=-0.75),
    ]


def ar_kr_manifold() -> list:
    """The Ar-fill heteronuclear manifold — FOUND MISSING 2026-08-04.

    The Xu 1993 re-verification pass exposed it: Xu p.2595 prices the
    KrF*+2Ar three-body at (2-8)e-32 cm6/s "because of ArKrF* formation",
    and Kannari 1985 Table I p.4311 (read first-hand) carries the ENTIRE
    chain — KrF*+2Ar -> ArKrF* (7.0e-32 [62,63]), the KrF*+Ar two-body
    quench (1.8e-12), ArKr*, and the Ar2F* fate block — every row of which
    was absent here. J&H 1980 Table I p.2408 (#11,13,14,16,17,19,20,23,
    24,26,27,28,30,33) double-sources the set; their Fig. 2 literally
    draws the ArKr*/ArKrF* boxes. Root cause of the omission: the REPO
    source network (spectre) is Ne-buffer-oriented, and the earlier
    per-row passes filtered tables by species already in the network, so
    rows of UNTRACKED species were invisible. At the v2 point the missing
    KrF* sinks total ~4.6e8 s^-1 per state — x3.2 the B radiative rate;
    this was the largest single defect in the Ar-fill answer.

    Values are Kannari's prints (300 K), SI. Where J&H print a different
    value the conflict is recorded on the row and the later table is
    kept (the AK1 precedent). The ArKr+ ION manifold — the 0.1.9
    OMISSIONS follow-up — is closed in ar_kr_ion_manifold() (0.1.10),
    and the 0.1.10 re-walk of the same two tables picked up the neutral
    rows this pass missed (ar_fill_completion()).
    """
    return [
        # --- ArKr*: formation, harpoon, displacement, VUV decay ----------
        R("AKD1", ("Kr*", "Ar", "Ar"), ("ArKr*", "Ar"), ("const", 1.0e-44),
          "LIT", "Kannari 1985 Table I p.4311: Kr*+2Ar -> ArKr*+Ar "
          "1.0e-32 cm6/s [ref 52]; J&H 1980 #28 same value [33] "
          "(both first-hand 2026-08-04). At the v2 point this diverts "
          "~45% of the Kr* harpoon flux through the dimer chain"),
        R("AKD2", ("Ar*", "Kr", "Ar"), ("ArKr*", "Ar"), ("const", 1.0e-44),
          "LIT", "Kannari 1985 Table I p.4311: Ar*+Kr+Ar -> ArKr*+Ar "
          "1.0e-32 cm6/s [ref 52]; J&H 1980 #30 same value [33] "
          "(both first-hand 2026-08-04)"),
        R("AKD3", ("ArKr*", "F2"), ("KrF*", "Ar", "F"), ("const", 6.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: ArKr*+F2 -> KrF*+Ar+F "
          "6.0e-10 cm3/s [ref 52] (first-hand 2026-08-04). RECORDED "
          "CONFLICT: J&H 1980 #11 prints 1.0e-9 [33] — x1.7; the later "
          "table is kept (AK1 precedent)"),
        R("AKD4", ("ArKr*", "Kr"), ("Kr2*", "Ar"), ("const", 1.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: ArKr*+Kr -> Kr2*+Ar "
          "1.0e-10 cm3/s [ref 52]; J&H 1980 #27 same value [33] "
          "(both first-hand 2026-08-04) — the dominant ArKr* fate at "
          "5% Kr (3.7e8 s^-1 at the v2 point)"),
        R("AKD_RAD", ("ArKr*",), ("Ar", "Kr"), ("const", 8.0e7), "LIT",
          "Kannari 1985 Table I p.4312: ArKr* -> Ar+Kr+hv 8.0e7 s^-1, "
          "marked 'estimated' there (first-hand 2026-08-04). RECORDED "
          "CONFLICT: J&H 1980 #70 prints 3.0e6 [33] — x27 apart, the "
          "triplet-vs-mixed-manifold question (same spread class as "
          "KR2_RAD/AR2_RAD); Kannari's model value kept for internal "
          "consistency, swept in the bracket sensitivity. Radiates the "
          "full formation energy (no band printed in either source)",
          kind="rad"),
        # --- the KrF* sinks Xu's pass exposed ----------------------------
        R("KAA", ("KrF*", "Ar", "Ar"), ("ArKrF*", "Ar"), ("const", 7.0e-44),
          "LIT", "Kannari 1985 Table I p.4311: KrF*+2Ar -> ArKrF*+Ar "
          "7.0e-32 cm6/s [refs 62,63]; J&H 1980 #14 prints 8.0e-32 "
          "[32,38]; Xu 1993 p.2595 quotes the class at (2-8)e-32 "
          "'because of ArKrF* formation' (all first-hand 2026-08-04). "
          "~3.4e8 s^-1 per KrF* state at the v2 point — x2.4 the B "
          "radiative rate; THE dominant formerly-missing sink. 0.1.40, "
          "read first-hand 2026-09-18. MAGNITUDE: Kannari's 7.0e-32 is "
          "the UNREDUCED k*tau / 6.5 ns of Rokni, Jacob & Mangano (PRA "
          "16, 2216 (1977); manuscript Table 3, DTIC ADA047221 p.27: "
          "k*tau 4.6e-40 cm6 -> 7e-32 at tau_R = 6.5 ns — equal to "
          "Kannari's print, as are his F2/2Kr/Kr+Ar rows), so his "
          "footnote-a 0.433 is applied INSIDE his model, not in the "
          "table. Where a DIRECT decay measurement exists the unreduced "
          "scale is the one it confirms (2Kr 6.7 vs Gadomski/Xu 6.0+-2.0; "
          "Kr+Ar 6.5 vs Xu 5.5+-2.0 e-31) while k*tau/15 ns (Kannari's "
          "model) and k*tau/25 ns (Xu's footnote-d adjustment, and the "
          "2026-09-04 audit's x0.297) miss it x2.1 and x3.3 — so the "
          "300 K value is KEPT. Tg FORM (Kannari footnote b, dropped on "
          "adoption, restored): Shui's R = Ar calculation, x1/8 from 300 "
          "to 600 K (KRF_3B_TG_POW_AR), measured to 15 % by Klimek & "
          "Hsia 1981",
          tg_pow=KRF_3B_TG_POW_AR),
        R("KQAr", ("KrF*", "Ar"), ("Kr", "Ar", "F"), ("const", 1.8e-18),
          "LIT", "Kannari 1985 Table I p.4311: KrF*+Ar -> Kr+Ar+F "
          "1.8e-12 cm3/s, marked 'estimated' there. THREE-SOURCE RECORD "
          "(all first-hand 2026-08-04): J&H 1980 #13 prints 0.0 "
          "(Estimated); Xu 1993 Table I/p.2596 measures (2.9+-0.5)e-12 "
          "as an impurity-limited UPPER BOUND. Kannari's estimate sits "
          "inside [0, 2.9e-12] and is kept; bracket = that interval"),
        # --- ArKrF*: fates -----------------------------------------------
        R("AKF_F2", ("ArKrF*", "F2"), ("Ar", "Kr", "F", "F2"),
          ("const", 1.0e-15), "LIT",
          "Kannari 1985 Table I p.4311: ArKrF*+F2 -> Ar+Kr+F+F2 1.0e-9 "
          "cm3/s [ref 52]; J&H 1980 #19 same value [33] (both "
          "first-hand 2026-08-04)"),
        R("AKF_KR", ("ArKrF*", "Kr"), ("Kr2F*", "Ar"), ("const", 2.0e-17),
          "LIT", "Kannari 1985 Table I p.4311: ArKrF*+Kr -> Kr2F*+Ar "
          "2.0e-11 cm3/s [ref 43]; J&H 1980 #16 same value [32] (both "
          "first-hand 2026-08-04)"),
        R("AKF_AR", ("ArKrF*", "Ar"), ("Ar2F*", "Kr"), ("const", 2.0e-17),
          "LIT", "Kannari 1985 Table I p.4311: ArKrF*+Ar -> Ar2F*+Kr "
          "2.0e-11 cm3/s [ref 43]; J&H 1980 #17 same value [32] (both "
          "first-hand 2026-08-04) — the dominant ArKrF* fate in an Ar "
          "fill (1.4e9 s^-1 at the v2 point)"),
        R("AKF_RAD", ("ArKrF*",), ("Ar", "Kr", "F"), ("const", 5.0e6),
          "LIT", "Kannari 1985 Table I p.4312: ArKrF* -> Ar+Kr+F+hv "
          "5.0e6 s^-1, marked 'estimated' there (first-hand 2026-08-04; "
          "Kr2F*-class, consistent with Xu's tau(Kr2F*) = 220 ns). "
          "RECORDED CONFLICT: J&H 1980 #71 prints 5.0e7 [33] — x10, "
          "typography-suspect (R36 class); the triatomic-class value is "
          "kept. Band unprinted in both: hv taken at the Kr2F*/Ar2F* "
          "EST interpolation 3.7 eV, remainder to gas",
          kind="rad", hv_eV=3.7),
        # --- Ar2F*: the fate block (absent since 0.1.0 — Ar2F* was a
        #     tracked absorber with NO consumers) ------------------------
        R("A2F_F2", ("Ar2F*", "F2"), ("Ar", "Ar", "F", "F2"),
          ("const", 2.0e-16), "LIT",
          "Kannari 1985 Table I p.4311: Ar2F*+F2 -> 2Ar+F+F2 2.0e-10 "
          "cm3/s [ref 61] (first-hand 2026-08-04; the value phantom once "
          "transposed onto K125). RECORDED CONFLICT: J&H 1980 #25 "
          "prints 1.0e-9 [33] — x5, the later table kept"),
        R("A2F_KR1", ("Ar2F*", "Kr"), ("ArKrF*", "Ar"), ("const", 1.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar2F*+Kr -> ArKrF*+Ar "
          "1.0e-10 cm3/s [refs 52,59]; J&H 1980 #26 same value [33] "
          "(both first-hand 2026-08-04) — the displacement pair's "
          "return direction"),
        R("A2F_KR2", ("Ar2F*", "Kr"), ("KrF*", "Ar", "Ar"),
          ("const", 1.0e-16), "LIT",
          "Kannari 1985 Table I p.4311: Ar2F*+Kr -> KrF*+2Ar 1.0e-10 "
          "cm3/s, marked 'estimated' there (first-hand 2026-08-04) — "
          "the chain's KrF*-return branch"),
        R("A2F_AR", ("Ar2F*", "Ar"), ("Ar", "Ar", "Ar", "F"),
          ("const", 2.2e-20), "LIT",
          "Lee, Matsui, Kannari & Obara, IEEE TED 36, 2053 (1989) Table I "
          "(Ar2F* relaxation): Ar2F*+Ar -> 2Ar+F 2.2e-14 cm3/s [their "
          "58] (read first-hand on the page image 2026-09-18) — Kannari's "
          "OWN group replacing the 1.0e-10 his 1985 Table I marks "
          "'estimated' (first-hand 2026-08-04; RECORDED, kept reachable "
          "in the 0.1.40 changelog). NRL Orestes carries the same "
          "2.2e-14 (petrova2022 Table A5: 2021 'Lee-1989'; the 2022 red "
          "row re-routes products to ArF*+2Ar, 'estimated products', "
          "Boewering-1982 — not adopted, Lee's products kept). "
          "CONSEQUENCE, measured 0.1.40: Ar2F* now lives ~1e-6 s in an Ar "
          "fill instead of flashing through, so in KrF mixtures the "
          "estimated A2F_KR2 return (and the AKF_AR feed into Ar2F*) "
          "become live — see the changelog's dependency statement"),
        R("A2F_RAD", ("Ar2F*",), ("Ar", "Ar", "F"), ("const", 5.4e6),
          "LIT", "Kannari 1985 Table I p.4312: Ar2F* -> 2Ar+F+hv "
          "5.4e6 s^-1 [refs 61,93] (first-hand 2026-08-04; the rate "
          "phantom once transposed onto Kr2F*). J&H 1980 #69 prints "
          "2.0e6 [33], recorded. hv at the registry's 285 nm band",
          kind="rad", hv_eV=4.35),
        # --- dimer VUV decay: the 0.1.0 'spectre parity' scope hole,
        #     closed with the printed rates ------------------------------
        R("KR2_RAD", ("Kr2*",), ("Kr", "Kr"), ("const", 8.0e7), "LIT",
          "Kannari 1985 Table I p.4312: Kr2* -> 2Kr+hv 8.0e7 s^-1 "
          "[ref 88] (first-hand 2026-08-04). RECORDED CONFLICT: J&H "
          "1980 #66 prints 3.3e6 [33] — x24: the triplet lifetime vs "
          "the mixed-manifold effective rate; Kannari's model value "
          "kept, swept in the bracket sensitivity. Without ANY such "
          "row (the 0.1.0 scope), 100% of Kr2* reached KrF* via the "
          "harpoon; with it, the VUV branch takes ~2/3 at the v2 "
          "point. hv at the registry's 147 nm band",
          kind="rad", hv_eV=8.43),
        R("AR2_RAD", ("Ar2*",), ("Ar", "Ar"), ("const", 6.0e7), "LIT",
          "Kannari 1985 Table I p.4312: Ar2* -> 2Ar+hv 6.0e7 s^-1 "
          "[ref 88] (first-hand 2026-08-04). RECORDED CONFLICT: J&H "
          "1980 #67 prints 3.8e6 [33] — same x16 spread class as "
          "KR2_RAD, same adjudication. hv at the registry's 126 nm "
          "band", kind="rad", hv_eV=9.84),
        # --- the second harpoon branches (both tables list them) --------
        R("K2F_BR", ("Kr2*", "F2"), ("Kr2F*", "F"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Kr2*+F2 -> Kr2F*+F "
          "3.0e-10 cm3/s [ref 52]; J&H 1980 #20 same value [33] (both "
          "first-hand 2026-08-04). Runs beside NK7's KrF*-forming "
          "branch [61] at the same k — the branch NK7's 0.1.3 "
          "adjudication did not cover"),
        R("A2_BR", ("Ar2*", "F2"), ("Ar2F*", "F"), ("const", 2.5e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar2*+F2 -> Ar2F*+F "
          "2.5e-10 cm3/s [ref 52]; J&H 1980 #23 same value [33] (both "
          "first-hand 2026-08-04). Runs beside AK6's ArF*-forming "
          "branch [61]"),
        R("K2F_F", ("Kr2*", "F"), ("KrF*", "Kr"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Kr2*+F -> KrF*+Kr "
          "3.0e-10 cm3/s; J&H 1980 #21 same value [33] (both "
          "first-hand 2026-08-04; grows with F2 burn-down)"),
        R("A2_F", ("Ar2*", "F"), ("ArF*", "Ar"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar2*+F -> ArF*+Ar "
          "3.0e-10 cm3/s [ref 52]; J&H 1980 #24 same value [33] (both "
          "first-hand 2026-08-04)"),
        # --- Ar2* -> Kr* collisional transfer (dominant Ar2* fate) ------
        R("AK9", ("Ar2*", "Kr"), ("Kr*", "Ar", "Ar"), ("const", 4.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar2*+Kr -> Kr*+2Ar "
          "4.0e-10 cm3/s [ref 43]; J&H 1980 #33 same value [32] (both "
          "first-hand 2026-08-04) — x20 every previously-present Ar2* "
          "channel at 5% Kr (1.5e9 s^-1 at the v2 point)"),
        # --- Ar2+ dissociative recombination (absent while Kr2+ had
        #     R11p) ------------------------------------------------------
        R("ADR", ("e", "Ar2+"), ("Ar*", "Ar"), ("dr_sqrt", 1.14e-13),
          "LIT", "Kannari 1985 Table I p.4312 DR block (first-hand "
          "2026-08-04): Ar2+ +e -> Ar*+Ar 0.6e-6 (Te/300K)^-0.66 [71] "
          "+ -> Ar**+Ar 1.1e-7 [68-70], lumped to the dominant printed "
          "branch at the 300 K total 7.1e-7 cm3/s on the R11p Te^-1/2 "
          "template (J&H #61 7.7e-7 corroborates; Ar** would relax to "
          "Ar* at 7e9 s^-1 anyway — downstream-equivalent)"),
    ]


def ar_fill_completion() -> list:
    """Neutral-block rows the 0.1.9 walk missed — 0.1.10, CORRECTED 0.1.11.

    Closing the ion sections (below) forced a fresh row-by-row re-walk of
    Kannari 1985 Table I pp.4311-4313 and J&H 1980 Table I p.2408. The
    0.1.10 framing ('krf had NO Ar-atom harpoon; AK1 was ArF*'s ONLY
    consumer') was WRONG: the CHEN block — generated from
    chen2025.TABLE_III, invisible to the source grep that made those
    claims — already carried R29 (harpoon), R30/R31 (ion-ion) and
    R32-R37 (ArF* fates), and the 0.1.10 additions ran double-counted
    beside them until the 0.1.11 dedup (CHEN_SUPERSEDED + the
    channel-uniqueness check). What this block truly contributes:

    - First-hand PRINTED values superseding Chen's model constants where
      they overlap (AF1 7.5e-10 vs R29's 4.0e-10; ARF_F2/ARF_AR/ARF_RAD
      vs R33/R34/R37 — conflicts recorded on the rows).
    - Genuinely absent rows: the Ar*/Ar** non-excimer branches (AF2/AF4),
      the Ar** harpoon (AF3), the Ar-third-body ArF*+2Ar -> Ar2F* channel
      (ARF_AA — R36 is the Ne-third-body variant), and AKD5.

    Kannari's ArF* quench prints are footnote-a REDUCED by the B/C factor
    0.418 already — consumed as printed, never re-reduced (the 0.433 trap,
    Kr side, applies here with the ArF number).
    """
    return [
        # --- the Ar*/Ar** + F2 harpoons (formation side) -----------------
        R("AF1", ("Ar*", "F2"), ("ArF*", "F"), ("const", 7.5e-16),
          "LIT", "Kannari 1985 Table I p.4311 row 1: Ar*+F2 -> ArF*+F "
          "7.5e-10 cm3/s [refs 52,59,92] (first-hand 2026-08-05). "
          "RECORDED CONFLICTS: J&H 1980 #4 prints 8.5e-10 [33] (x1.13); "
          "Chen 2025 R29 carries 4.0e-10 (x1.9 below) — this row "
          "SUPERSEDES R29 in this mode (CHEN_SUPERSEDED; the 0.1.10 "
          "claim that no Ar-atom harpoon existed was WRONG — R29 is "
          "generated from Chen's table and the source grep missed it; "
          "both ran double-counted until 0.1.11). 5.6e7 s^-1 per Ar* "
          "at the v2 point"),
        R("AF2", ("Ar*", "F2"), ("Ar", "F", "F"), ("const", 3.1e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar*+F2 -> Ar+2F 3.1e-10 "
          "cm3/s [ref 92] (first-hand 2026-08-05) — the non-excimer "
          "branch beside AF1, 29% of the pair"),
        R("AF3", ("Ar**", "F2"), ("ArF*", "F"), ("const", 4.7e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar**+F2 -> ArF*+F "
          "4.7e-10 cm3/s [refs 52,60] (first-hand 2026-08-05)"),
        R("AF4", ("Ar**", "F2"), ("Ar", "F", "F"), ("const", 3.1e-16),
          "LIT", "Kannari 1985 Table I p.4311: Ar**+F2 -> Ar+2F "
          "3.1e-10 cm3/s, marked 'analogous to Ar*' there (first-hand "
          "2026-08-05)"),
        # --- the third ArKr* row (0.1.9 adopted only two of three) -------
        R("AKD5", ("ArKr*", "F2"), ("ArKrF*", "F"), ("const", 3.0e-16),
          "LIT", "Kannari 1985 Table I p.4311: ArKr*+F2 -> ArKrF*+F "
          "3.0e-10 cm3/s, marked 'analogous to Ar2*' there (first-hand "
          "2026-08-05) — the adduct branch beside AKD3's 6.0e-10 "
          "KrF*-forming branch; MISSED by the 0.1.9 walk"),
        # --- the ArF* relaxation block (AK1 was the sole consumer) -------
        R("ARF_AA", ("ArF*", "Ar", "Ar"), ("Ar2F*", "Ar"),
          ("const", 4.0e-43), "LIT",
          "Kannari 1985 Table I p.4311 ArF* relaxation: ArF*+2Ar -> "
          "Ar2F*+Ar 4.0e-31 cm6/s [refs 62,63,90, footnote a/b] "
          "(first-hand 2026-08-05). 1.9e9 s^-1 at 3 atm — 24% "
          "of the ArF* budget, diverting it from AK1's KrF* transfer. "
          "CORRECTED READING (0.1.40): the printed value is NOT footnote-a "
          "reduced — it is Rokni's k*tau 1.6e-39 cm6 at tau = 4 ns "
          "(Avco manuscript Table 2, DTIC ADA047221 p.15, as are his "
          "ArF*+F2/Kr/Ar rows); the 0.418 lives in Kannari's model. For "
          "this engine's LUMPED ArF* radiating at tau_B (ARF_RAD) the "
          "tau_B-based value is the consistent one (Lee 1989's footnote c "
          "states the same convention), so it is kept. Footnote b "
          "('Expressed for 300 K gas temperature', Tg dependence per "
          "Shui) restored with NRL's printed form normalised at 300 K "
          "(ARF_3B_TG_POW/FLOOR). RECORDED CONFLICT on the 300 K "
          "magnitude: Lee 1989 Table I prints 5.0e-31 (300 K, tau_B "
          "based) and NRL 2022's form gives 5.06e-31; Kannari's 4.0e-31 "
          "kept (anchor effect of the swap measured x0.98 on Electra)",
          tg_pow=ARF_3B_TG_POW, tg_floor=ARF_3B_TG_FLOOR),
        R("ARF_F2", ("ArF*", "F2"), ("Ar", "F", "F2"), ("const", 1.9e-15),
          "LIT", "Kannari 1985 Table I p.4311 ArF* relaxation: ArF*+F2 "
          "-> Ar+F+F2 1.9e-9 cm3/s [refs 62,63, footnote a] "
          "(first-hand 2026-08-05). SUPERSEDES Chen 2025 R33 — the "
          "same 1.9e-9 with the F2 fragment DISSOCIATED (Ar+3F); "
          "Kannari's product routing kept, the variant recorded "
          "(CHEN_SUPERSEDED; double-counted until 0.1.11)"),
        R("ARF_AR", ("ArF*", "Ar"), ("Ar", "Ar", "F"), ("const", 9.0e-18),
          "LIT", "Kannari 1985 Table I p.4311 ArF* relaxation: ArF*+Ar "
          "-> 2Ar+F 9.0e-12 cm3/s [refs 62,63, footnote a] (first-hand "
          "2026-08-05). SUPERSEDES Chen 2025 R34 — the IDENTICAL "
          "9.0e-12 (CHEN_SUPERSEDED; double-counted until 0.1.11)"),
        # --- Ar** electron de-excitation (chain 2 of the amplification
        #     audit, 0.1.14): NRL's Orestes relaxes Ar** by ELECTRONS
        #     (petrova2022 App. A3) and carries NO neutral 4p->4s
        #     quench, while Kannari/Lee both print the neutral quench
        #     at 1.0e-10 (APS_Ar) — the structural disagreement heading
        #     the Ar** x0.06 / alpha x0.23 residuals. Both printed
        #     electron rows adopted; the neutral quench stays at its
        #     two-print value and is SWEPT (CONFLICT_ALTERNATES) ------
        R("AR_EDX1", ("Ar**", "e"), ("Ar*", "e"),
          ("te_fit", 2.69e-13, -0.17, 0.0), "LIT",
          "petrova2022 NRL/6720/MR--2022/1 App. A3 (first-hand "
          "2026-08-05): Ar**+e -> Ar*+e 2.69e-7 Te^-0.17 cm3/s "
          "[their k5, Vriens-1980]; superelastic — the 1.6 eV returns "
          "to the electron gas", energy_to="electrons"),
        R("AR_EDX2", ("Ar**", "e"), ("Ar", "e"),
          ("te_fit", 9.52e-16, -0.17, 0.0), "LIT",
          "petrova2022 NRL/6720/MR--2022/1 App. A3 (first-hand "
          "2026-08-05): Ar**+e -> Ar+e 9.52e-10 Te^-0.17 cm3/s "
          "[McGeoch-2018]; superelastic", energy_to="electrons"),
        R("ARF_RAD", ("ArF*",), ("Ar", "F"), ("const", 2.50e8),
          "LIT", "Kannari 1985 Table I p.4312: ArF*(B) -> Ar+F+hv "
          "2.50e8 s^-1 [refs 47,48] (first-hand 2026-08-05). "
          "SUPERSEDES Chen 2025 R37 (2.381e8, x0.95 — recorded) in "
          "this mode (CHEN_SUPERSEDED; the 0.1.10 claim that krf "
          "carried no ArF* spontaneous decay was WRONG — R37 was in "
          "the krf assembly via chen_rows, and both radiated until "
          "0.1.11). hv at the 193 nm quantum; the lumped-state "
          "registry reference (7.52 eV) covers it",
          kind="rad", hv_eV=6.42),
    ]


def ar_kr_ion_manifold() -> list:
    """The ArKr+ ION manifold + the ion-block completeness pass (0.1.10).

    Closes the 0.1.9 OMISSIONS follow-up ('the next network milestone'):
    the board's one open warn was the formation-side deficit exposed when
    the neutral Ar manifold landed, with the omitted ArKr+ manifold named
    the leading suspect (Kr+ +2Ar -> ArKr+ at 1.0e-31 is the DOMINANT
    Kr+ route in an Ar fill — 3.1e8 s^-1 at the Kannari device point vs
    1.1e8 for the tracked AK7 route). Both primaries print the full
    manifold and were read first-hand 2026-08-05: Kannari 1985 Table I
    p.4311-4312 (CT block, ion-ion block, DR block, Penning block,
    superelastic block) and J&H 1980 Table I p.2408 (#10, 36-40).

    Also here, same sections, same walk (framing CORRECTED 0.1.11): the
    AAR1/FMN1/AIIR1/AIIR2 rows SUPERSEDE Chen's generated R27/R22/R30/
    R31 — the 0.1.10 claims 'source-less Ar2+', 'sink-less F2+' and 'the
    primary beam ion had no F- path' were WRONG (those channels existed
    in the CHEN block, which is generated from chen2025.TABLE_III and
    invisible to the source grep; the pairs ran double-counted until the
    0.1.11 dedup, CHEN_SUPERSEDED). Genuinely new: the five printed
    Penning rows (the Ne network already carried its analog, NR4).

    NOT adopted, recorded: the e-impact ionization ladder of excited
    states [refs 72,79] except ArKr* (see AKI_ION), the Ar2*/Kr2*
    superelastics, F+F+M -> F2+M [16] (~1e3 s^-1 at pulse densities),
    and any ArKr+ 248-nm photoabsorption — NEITHER primary prints a
    sigma for the heteronuclear ion (the 0.1.9 OMISSIONS note that
    said 'ArKr+ photoabsorption 1.5e-17' had a species typo: Kannari's
    printed 1.5e-17 cm2 row [85] is the NEUTRAL ArKr*, adopted in
    mode_optics).
    """
    return [
        # --- ArKr+ formation (the dominant Kr+ route in an Ar fill) ------
        R("AKI1", ("Kr+", "Ar", "Ar"), ("ArKr+", "Ar"), ("const", 1.0e-43),
          "LIT", "Kannari 1985 Table I p.4312 CT block: Kr+ +2Ar -> "
          "ArKr+ +Ar 1.0e-31 cm6/s [ref 16]; J&H 1980 #37 same value "
          "(Estimated) (both first-hand 2026-08-05). THE dominant Kr+ "
          "route in an Ar fill — the 0.1.9 OMISSIONS quantification "
          "this row closes. Tg^-0.75 per the K15/Boichenko "
          "ion-association family template", tg_pow=-0.75),
        R("AKI2", ("Ar+", "Kr", "Ar"), ("ArKr+", "Ar"), ("const", 1.0e-43),
          "LIT", "Kannari 1985 Table I p.4312 CT block: Ar+ +Kr+Ar -> "
          "ArKr+ +Ar 1.0e-31 cm6/s [ref 66]; J&H 1980 #36 same value "
          "[39] (both first-hand 2026-08-05) — runs beside AK3's "
          "two-body CT", tg_pow=-0.75),
        R("AKI3", ("ArKr+", "Kr"), ("Kr2+", "Ar"), ("const", 3.2e-16),
          "LIT", "Kannari 1985 Table I p.4312 CT block: ArKr+ +Kr -> "
          "Kr2+ +Ar 3.2e-10 cm3/s [ref 16]; J&H 1980 #39 same value "
          "(Estimated) (both first-hand 2026-08-05) — displacement "
          "toward the lower-IP homonuclear dimer, the HK3 (NeKr+) "
          "analog; 5.9e7 s^-1 at 5% Kr, 3 atm"),
        # --- ArKr+ + F- ion-ion recombination (both printed branches) ----
        R("AKI_F", ("ArKr+", "F-"), ("KrF*", "Ar"),
          ("flannery", 2.65, 1.548e-6), "REPO",
          "existence LIT-printed: Kannari 1985 Table I p.4312 ion-ion "
          "block ArKr+ +F- (+Ar,Kr) -> KrF* (+M) [refs 28-35, "
          "Natanson-class 'see text']; J&H 1980 #10 prints 1.0e-6 "
          "cm3/s effective two-body at 2 atm (Estimated) (both "
          "first-hand 2026-08-05). RATE: the K108 (Kr2+ +F-) Flannery "
          "turnover form reused — Natanson/Flannery recombination is "
          "buffer-collision-dominated at these densities, weakly "
          "ion-specific — hence REPO tier; split 50/50 with AKI_FA per "
          "J&H's equal branch prints (total = the K108 form, 2.1e-6 "
          "cm3/s at 2 atm — in-band of #10+#38)"),
        R("AKI_FA", ("ArKr+", "F-"), ("ArKrF*",),
          ("flannery", 2.65, 1.548e-6), "REPO",
          "existence LIT-printed: Kannari 1985 Table I p.4312 ion-ion "
          "block ArKr+.F- (+M) -> ArKrF* (+M) — the adduct-stabilized "
          "branch; J&H 1980 #38 -> ArKrF* 1.0e-6 (Estimated) (both "
          "first-hand 2026-08-05). 50/50 split with AKI_F per the "
          "equal prints; params mirrored from K108 (REPO)"),
        # --- the missing Ar-ion + F- recombination pair ------------------
        R("AIIR1", ("Ar+", "F-"), ("ArF*",), ("flannery", 2.4, 3.417e-6),
          "REPO", "existence LIT-printed: Kannari 1985 Table I p.4312 "
          "ion-ion block Ar+ +F- (+Ar,Kr) -> ArF* (+M) [refs 28-35]; "
          "J&H 1980 #5 prints 1.0e-6 [44] (both first-hand "
          "2026-08-05). SUPERSEDES Chen 2025 R30 (9.1e-7 const) in "
          "this mode: the Flannery turnover adds the pressure "
          "dependence Chen's operating-point constant lacks "
          "(CHEN_SUPERSEDED; the 0.1.10 claim that the primary beam "
          "ion had no F- path was WRONG — R30 existed, generated from "
          "Chen's table; both ran double-counted until 0.1.11). "
          "Params mirrored from K26 (atomic-ion class, REPO)"),
        R("AIIR2", ("Ar2+", "F-"), ("ArF*", "Ar"),
          ("flannery", 2.65, 3.096e-6), "REPO",
          "existence LIT-printed: Kannari 1985 Table I p.4312 ion-ion "
          "block Ar2+ +F- (+Ar,Kr) -> ArF*+Ar (+M) [refs 28-35]; J&H "
          "1980 #9 prints 1.0e-6 [44] (both first-hand 2026-08-05). "
          "SUPERSEDES Chen 2025 R31 (9.26e-7 const) in this mode "
          "(CHEN_SUPERSEDED, as AIIR1). Params mirrored from K108 "
          "(molecular-ion class, REPO)"),
        # --- ArKr+ dissociative recombination ----------------------------
        R("AKI_DR", ("e", "ArKr+"), ("Kr**", "Ar"), ("dr_sqrt", 1.93e-13),
          "LIT", "Kannari 1985 Table I p.4312 DR block: ArKr+ +e -> "
          "Kr**+Ar 1.2e-6 (Te/300K)^-0.55 cm3/s [ref 71], mapped at "
          "the 300 K print onto the dr_sqrt Te^-1/2 template (the "
          "ADR/R11p convention; 1.2e-12/6.219 = 1.93e-13). J&H 1980 "
          "#40 prints typical 1.0e-7 — in-band at ~1 eV (both "
          "first-hand 2026-08-05). Single printed branch, product "
          "Kr** as printed"),
        # --- ArKr* <-> ArKr+ closure rows --------------------------------
        R("AKI_ION", ("ArKr*", "e"), ("ArKr+", "e", "e"),
          ("te_fit", 9.00e-14, 0.70, -3.60), "LIT",
          "Kannari 1985 Table I p.4311: ArKr*+e -> ArKr+ +2e "
          "9.00e-8 Te^0.70 exp(-3.60/Te) cm3/s [ref 72] (first-hand "
          "2026-08-05) — completes the printed fate set of the 0.1.9 "
          "ArKr* species; its [72] siblings for Ar2*/Kr2* are recorded "
          "in OMISSIONS (sub-1% of those species' fate budgets)"),
        R("AKI_SE", ("ArKr*", "e"), ("Ar", "Kr", "e"), ("const", 1.0e-13),
          "LIT", "Kannari 1985 Table I p.4312 superelastic: ArKr*+e -> "
          "Ar+Kr+e 1.0e-7 cm3/s, marked 'estimated' there (first-hand "
          "2026-08-05; the K_2Fe class)"),
        R("AKF_SE", ("ArKrF*", "e"), ("Ar", "Kr", "F", "e"),
          ("const", 1.0e-13), "LIT",
          "Kannari 1985 Table I p.4312 superelastic: ArKrF*+e -> "
          "Ar+Kr+F+e 1.0e-7 cm3/s, marked 'estimated' there "
          "(first-hand 2026-08-05)"),
        # --- ion-block source-less/sink-less fixes -----------------------
        R("AAR1", ("Ar+", "Ar", "Ar"), ("Ar2+", "Ar"), ("const", 2.5e-43),
          "LIT", "Kannari 1985 Table I p.4312 CT block: Ar+ +2Ar -> "
          "Ar2+ +Ar 2.5e-31 cm6/s [ref 65] (first-hand 2026-08-05). "
          "SUPERSEDES Chen 2025 R27 — the IDENTICAL 2.5e-31, same "
          "lineage (CHEN_SUPERSEDED; the 0.1.10 claim that Ar2+ was "
          "source-less was WRONG — R26/R27 existed, generated from "
          "Chen's table; the channel ran double-counted until "
          "0.1.11). R26 (Ne third body) remains beside this. "
          "Tg^-0.75 per the K15 family template", tg_pow=-0.75),
        R("FMN1", ("F2+", "F-"), ("F", "F", "F"), ("const", 4.0e-14),
          "LIT", "Kannari 1985 Table I p.4312 ion-ion block: F2+ +F- "
          "-> 3F 4.0e-8 cm3/s [ref 91] (first-hand 2026-08-05). "
          "SUPERSEDES Chen 2025 R22 — the IDENTICAL 4.0e-8 "
          "(CHEN_SUPERSEDED; the 0.1.10 claim that F2+ was sink-less "
          "was WRONG — R22 existed, generated; double-counted until "
          "0.1.11)"),
        # --- Penning ionization (the five printed rows; NR4 precedent) ---
        R("PEN_KR", ("Kr*", "Kr*"), ("Kr+", "Kr", "e"), ("const", 5.0e-16),
          "LIT", "Kannari 1985 Table I p.4312 Penning block: Kr*+Kr* "
          "-> Kr+ +Kr+e 5.0e-10 cm3/s [ref 86] (first-hand "
          "2026-08-05; quadratic in [Kr*] — the Ne network already "
          "carried its NR4 analog)"),
        R("PEN_AR", ("Ar*", "Ar*"), ("Ar+", "Ar", "e"), ("const", 5.0e-16),
          "LIT", "Kannari 1985 Table I p.4312 Penning block: Ar*+Ar* "
          "-> Ar+ +Ar+e 5.0e-10 cm3/s [ref 86] (first-hand 2026-08-05)"),
        R("PEN_ARSS", ("Ar**", "Ar**"), ("Ar+", "Ar", "e"),
          ("const", 5.0e-16), "LIT",
          "Kannari 1985 Table I p.4312 Penning block: Ar**+Ar** -> "
          "Ar+ +Ar+e 5.0e-10 cm3/s, marked 'estimated' there "
          "(first-hand 2026-08-05)"),
        R("PEN_KR2", ("Kr2*", "Kr2*"), ("Kr2+", "Kr", "Kr", "e"),
          ("const", 5.0e-16), "LIT",
          "Kannari 1985 Table I p.4312 Penning block: Kr2*+Kr2* -> "
          "Kr2+ +2Kr+e 5.0e-10 cm3/s [ref 87] (first-hand 2026-08-05)"),
        R("PEN_AR2", ("Ar2*", "Ar2*"), ("Ar2+", "Ar", "Ar", "e"),
          ("const", 5.0e-16), "LIT",
          "Kannari 1985 Table I p.4312 Penning block: Ar2*+Ar2* -> "
          "Ar2+ +2Ar+e 5.0e-10 cm3/s [ref 87] (first-hand 2026-08-05)"),
    ]


def kr_loki_rows(available_buckets) -> list:
    """e + Kr ground-manifold rows from the mixture tables.

    A COMPLETENESS DELTA OVER SPECTRE, recorded: spectre's KrF-mode network
    carries no e+Kr ground excitation/ionization rows at all (its Kr
    manifolds are populated by the beam only). In a discharge — and in any
    avalanche — those channels exist; the curated V2point Kr/F2/Ne(Ar)
    tables price them, so the rows belong in the network."""
    avail = set(available_buckets)
    rows = []
    for rid, pro, b in (("LK1", ("Kr*", "e"), "Kr*"),
                        ("LK2", ("Kr**", "e"), "Kr**"),
                        ("LK3", ("Kr+", "e", "e"), "Kr+")):
        if b in avail:
            rows.append(R(rid, ("e", "Kr"), pro, ("loki", b), "LOKI",
                          "k(<eps>) from the curated LoKI-MC V2point mixture "
                          "tables (PLAN (45)); bucket " + b))
    return rows


def upper_level_radiation(mode: str) -> list:
    """Spontaneous decay of the upper laser level.

    KrF mode (the Kannari 1985 B/C/Bv split, Table I verified first-hand
    2026-08-04): B radiates at 1.43e8 s^-1 (tau_B = 6.99 ns) and lases;
    C radiates C->A at 1.33e7 s^-1 (tau_C = 75.19 ns) and does NOT lase;
    the hot-B manifold Bv radiates at the B rate (per-v lifetime taken
    v-independent, Kannari Sec II C). This replaces the lumped 25 ns /
    sigma_eff pair (Gadomski/Xu), which encoded a DIFFERENT B fraction
    (0.206) than Kannari's Boltzmann split (0.405) — the two-state model
    carries the split explicitly instead of folding it into constants.
    ArF mode: Chen R37 is already in the table."""
    if mode == "krf":
        return [
            R("URAD_B", ("KrF*(B)",), ("Kr", "F"), ("const", 1.43e8),
              "LIT", "Kannari 1985 Table I p.4312: KrF*(B)->Kr+F+hv "
              "1.43e8 s^-1 [refs 47,48] (first-hand 2026-08-04)",
              kind="rad", hv_eV=5.00),
            R("URAD_C", ("KrF*(C)",), ("Kr", "F"), ("const", 1.33e7),
              "LIT", "Kannari 1985 Table I p.4312: KrF*(C)->Kr+F+hv "
              "1.33e7 s^-1 [refs 47,48] (first-hand 2026-08-04; broadband "
              "C->A, booked at ~the 248 nm quantum)",
              kind="rad", hv_eV=4.99),
            R("URAD_Bv", ("KrF*(Bv)",), ("Kr", "F"), ("const", 1.43e8),
              "EST", "hot-B manifold at the B radiative rate (Kannari "
              "Sec II C: per-v lifetimes treated v-independent)",
              kind="rad", hv_eV=5.00),
        ]
    return []


def krf_bc_mixing() -> list:
    """The Kannari 1985 B<->C mixing and vibrational relaxation.

    Mixing 5e-10 cm3/s (two-body with Ar; abstract + Sec III), detailed
    balance from E_B - E_C = 80 cm^-1: k(C->B) = k(B->C) exp(-DE/kT) =
    0.682 at 300 K (Tg dependence unmodeled — Tg moves < 10% in runs, a
    stated approximation). Vibrational relaxation Bv -> B at 4e-11 cm3/s
    (abstract: "estimated ... 4e-11 ... for a two-body collision rate
    with argon gas"). Ne-partner analogs EST at the class value.

    0.1.38: Kr-partner analogs added. The 0.1.4 split carried Ar (LIT) and
    Ne (EST) partners only, so a fill with little or no argon/neon had no
    channel from C or from the hot manifold Bv into the lasing state B: a
    pure Kr/F2 fill could not lase at all (Mandl 1985 measured 0.82 of the
    10 %-Kr optimum there) and a 40 %-Kr Electra fill was starved of
    relaxation (327 -> 402 J at fixed beam once the rows exist). Found by
    the PIC/Neon replication (NRL_REPLICATION.md §6d). Same policy as
    kannari_ps_relaxation: buffer-partner analogs EST at the class value,
    bracket direction recorded, primary named."""
    return [
        R("BC_Ar", ("KrF*(B)", "Ar"), ("KrF*(C)", "Ar"), ("const", 5.0e-16),
          "LIT", "Kannari 1985 abstract/Sec III: KrF*(B,C) mixing 5e-10 "
          "cm3/s with Ar (first-hand 2026-08-04)"),
        R("CB_Ar", ("KrF*(C)", "Ar"), ("KrF*(B)", "Ar"), ("const", 3.41e-16),
          "LIT", "detailed balance of BC_Ar: x exp(-80 cm^-1 / kT_300) = "
          "0.682 (Kannari p.4314 energy gap)"),
        R("BC_Ne", ("KrF*(B)", "Ne"), ("KrF*(C)", "Ne"), ("const", 5.0e-16),
          "EST", "Ne-partner analog of BC_Ar at the class value. Bracket "
          "direction recorded 2026-08-04: Xu 1993 p.2596 — 'the k_BC "
          "values for He, Ne, and Ar may be much smaller than for Kr' "
          "(sensitization studies [17,24]); the analog likely "
          "OVERestimates"),
        R("CB_Ne", ("KrF*(C)", "Ne"), ("KrF*(B)", "Ne"), ("const", 3.41e-16),
          "EST", "Ne-partner analog of CB_Ar at the class value (detailed "
          "balance partner of BC_Ne; same Xu 1993 p.2596 bracket "
          "direction)"),
        R("VREL_Ar", ("KrF*(Bv)", "Ar"), ("KrF*(B)", "Ar"), ("const", 4.0e-17),
          "LIT", "Kannari 1985 abstract: vibrational relaxation 4e-11 "
          "cm3/s, two-body with Ar (first-hand 2026-08-04; the finite "
          "rate his paper exists to establish, 2e-11..5e-11 band)"),
        R("VREL_Ne", ("KrF*(Bv)", "Ne"), ("KrF*(B)", "Ne"), ("const", 4.0e-17),
          "EST", "Ne-partner analog of VREL_Ar at the class value"),
        # --- 0.1.38: the krypton partners the split left out --------------
        R("BC_Kr", ("KrF*(B)", "Kr"), ("KrF*(C)", "Kr"), ("const", 5.0e-16),
          "EST", "Kr-partner analog of BC_Ar at the class value (0.1.38). "
          "Bracket direction: Xu 1993 p.2596 — the k_BC values for He, Ne "
          "and Ar 'may be much smaller than for Kr', so the class value "
          "likely UNDERestimates the Kr partner; primary for the Kr value "
          "= Gadomski/Xu/Setser/Nelson CPL 189,153 (1992), adopt when read "
          "first-hand"),
        R("CB_Kr", ("KrF*(C)", "Kr"), ("KrF*(B)", "Kr"), ("const", 3.41e-16),
          "EST", "detailed-balance partner of BC_Kr: x exp(-80 cm^-1 / "
          "kT_300) = 0.682, same gap as CB_Ar (0.1.38)"),
        R("VREL_Kr", ("KrF*(Bv)", "Kr"), ("KrF*(B)", "Kr"), ("const", 4.0e-17),
          "EST", "Kr-partner analog of VREL_Ar at the class value (0.1.38); "
          "primary = Kvaran/Shaw/Simons Appl. Phys. B 46,95 (1988), "
          "level-resolved relaxation of KrF* by rare gases, adopt when read"),
    ]


def _halve_rate(spec):
    kind = spec[0]
    if kind == "const":
        return ("const", spec[1] * 0.5)
    if kind == "flannery":
        return ("flannery", spec[1], spec[2] * 0.5)
    if kind == "lee_ii":
        return ("lee_ii", spec[1] * 0.5, spec[2], spec[3])
    raise NotImplementedError(
        f"cannot halve a {kind!r} rate spec for the B/C formation split")


def split_krf_bc(rows) -> list:
    """Rewrite every lumped-KrF* row into the Kannari B/C/Bv states.

    MECHANICAL, not hand-listed: any row CONSUMING KrF* becomes three rows
    (same k per state — the quench constants in this network are the
    Xu/Gadomski/Kushner measurements, already effective-lifetime-consistent,
    so no Kannari 0.433 rescale is applied: doing both would double-count);
    any row PRODUCING KrF* splits 0.5 -> Bv (nascent KrF* is born
    vibrationally hot) and 0.5 -> C (Kannari: "the formation rate for the
    KrF*(B) and KrF*(C) states is the same value" [41])."""
    B, C, BV = "KrF*(B)", "KrF*(C)", "KrF*(Bv)"
    out = []
    for r in rows:
        in_rea = "KrF*" in r.reactants
        in_pro = "KrF*" in r.products
        if not (in_rea or in_pro):
            out.append(r)
            continue
        if in_rea and in_pro:
            raise NetworkError(f"{r.rid}: KrF* on both sides — the split "
                               "transform does not define this case")
        if in_rea:
            for st, tag in ((B, "_B"), (C, "_C"), (BV, "_Bv")):
                rea = tuple(st if s == "KrF*" else s for s in r.reactants)
                out.append(replace(r, rid=r.rid + tag, reactants=rea,
                                   cite=r.cite + " ; applied per B/C/Bv "
                                   "state (Kannari split, same k)"))
        else:
            for st, tag in ((BV, "_fBv"), (C, "_fC")):
                pro = tuple(st if s == "KrF*" else s for s in r.products)
                out.append(replace(r, rid=r.rid + tag, products=pro,
                                   rate=_halve_rate(r.rate),
                                   cite=r.cite + " ; formation split 0.5 "
                                   "-> " + st + " (Kannari [41]: equal "
                                   "B/C formation; B-side born hot)"))
    return out


# ---------------------------------------------------------------------------
# mode assembly
# ---------------------------------------------------------------------------
# CHEN rows whose physical channel the adopted first-hand set supersedes in
# the Kannari-chemistry modes. Found 2026-08-05: the 0.1.10 completion pass
# added AF1/AIIR1/AIIR2/AAR1/FMN1/ARF_F2/ARF_AR/ARF_RAD blind to these —
# row existence had been checked by grepping network.py SOURCE, but the
# CHEN block is GENERATED from chen2025.TABLE_III, so 8 channels ran
# double-counted. In krf (and any future Kannari-chemistry mode) the
# first-hand LIT rows win and these are dropped; arf_chen keeps Chen's set
# VERBATIM — that anchor scores THEIR model, not a hybrid. Conflicts are
# recorded on the superseding rows.
CHEN_SUPERSEDED = {
    "R22": "FMN1",     # F2+ +F- -> 3F (identical value, Kannari [91] lineage)
    "R27": "AAR1",     # Ar+ +2Ar -> Ar2+ (identical 2.5e-31)
    "R29": "AF1",      # Ar*+F2 harpoon (Chen 4.0e-10 vs the 7.5e-10 print)
    "R30": "AIIR1",    # Ar+ +F- (Chen const vs Flannery turnover)
    "R31": "AIIR2",    # Ar2+ +F- (same)
    "R33": "ARF_F2",   # ArF*+F2 quench (same 1.9e-9; product routing differs)
    "R34": "ARF_AR",   # ArF*+Ar (identical 9.0e-12)
    "R37": "ARF_RAD",  # ArF* spontaneous (Chen 2.381e8 vs Kannari 2.50e8)
}

# F2's a3Piu and A1Piu are REPULSIVE states: excitation IS dissociation
# within a vibrational period. Chen's Table III books them as bound
# species that wait for a SECOND electron collision (R13/R14) — at
# beam-mode ne that parks excited F2 for ~us and froze the burn-down
# x0.41 low vs the NRL benchmark (found 2026-08-05, chain 3 of the
# amplification audit). In the Kannari-chemistry modes the measured
# excitation rows deliver the fragments directly: same LOKI rate, the
# electron still pays the full threshold (D0(F2) = 1.6 eV stored in the
# fragments via release, the excess as fragment heat via
# e_cost_override), and R13/R14 are dropped (source-less). arf_chen
# keeps Chen's two-step bookkeeping VERBATIM.
_REPULSIVE_F2 = {"R11": 3.16, "R12": 4.34}     # rid -> threshold [eV]


def reroute_repulsive_f2(rows) -> list:
    out = []
    for r in rows:
        if r.rid in ("R13", "R14"):
            continue
        if r.rid in _REPULSIVE_F2:
            thr = _REPULSIVE_F2[r.rid]
            r = R(r.rid, r.reactants, ("F", "F", "e"), r.rate, r.tier,
                  r.cite + " ; REROUTED 0.1.14: the target state is "
                  "REPULSIVE — excitation IS prompt dissociation "
                  "(fragments direct; threshold paid in full, excess "
                  f"over D0 = {thr - 1.6:.2f} eV to fragment heat)",
                  e_cost_override_eV=thr - 1.6)
        out.append(r)
    return out


def build_network(mode: str, available_buckets,
                  chen_errata: bool = True) -> list:
    """The full network for a mode. ``available_buckets`` comes from the
    loaded electron-rate tier (LOKI or KINETIC), so a missing table degrades
    to the documented fallback rather than a silent zero.

    ``chen_errata`` (0.1.25) is threaded to every chen_rows() call — see
    its docstring. In krf/arf_beam the affected rows are either dropped
    (R13/R14, repulsive reroute) or Ne-partnered (R23/R24/R36: flux ~0 in
    an Ar fill, live in a Ne fill — where the as-printed R24 annihilated
    Ne* in ~0.4 ps), so the Ar-fill anchors are unchanged by the default."""
    if mode == "arf_chen":
        rows = (chen_rows(available_buckets, errata=chen_errata)
                + ne_refinements())
    elif mode == "arf_beam":
        # The Electra ArF anchor mode (0.1.11): the SAME Kannari-chemistry
        # pool as krf, mechanically restricted to the ArF species list —
        # no hand-picked row set to drift. Chen's e-impact block included;
        # the CHEN_SUPERSEDED channels excluded (the krf dedup rule —
        # arf_chen alone keeps Chen verbatim); every row whose species all
        # exist in ARF_MODE_SPECIES survives, the rest (Kr/NeF chains)
        # fall away by the filter, not by judgment.
        from .species import ARF_MODE_SPECIES
        allowed = set(ARF_MODE_SPECIES)
        pool = (chen_rows(available_buckets, errata=chen_errata) + kr_core()
                + kannari_ps_relaxation() + ne_kr_couplings()
                + ne_refinements() + ar_kr_couplings() + ar_kr_manifold()
                + ar_fill_completion() + ar_kr_ion_manifold()
                + nef_channels() + kr_loki_rows(available_buckets))
        pool = [r for r in pool if r.rid not in CHEN_SUPERSEDED]
        pool = reroute_repulsive_f2(pool)

        def _fits(r):
            return all(s in allowed or s == ELECTRON_NAME
                       or s.startswith(PHOTON_PREFIX)
                       for s in r.reactants + r.products)
        rows = [r for r in pool if _fits(r)]
        rids = {r.rid for r in rows}
        missing = [v for v in CHEN_SUPERSEDED.values() if v not in rids]
        if missing:
            raise NetworkError(f"arf_beam lost superseding rows {missing} "
                               f"to the species filter — dedup broken")
    elif mode == "krf":
        rows = (chen_rows(available_buckets, errata=chen_errata) + kr_core()
                + kannari_ps_relaxation() + ne_kr_couplings()
                + ne_refinements() + ar_kr_couplings() + ar_kr_manifold()
                + ar_fill_completion() + ar_kr_ion_manifold()
                + nef_channels() + kr_loki_rows(available_buckets))
        rids = {r.rid for r in rows}
        missing = [v for v in CHEN_SUPERSEDED.values() if v not in rids]
        if missing:
            raise NetworkError(f"CHEN_SUPERSEDED replacements absent: "
                               f"{missing} — refusing to drop Chen rows "
                               f"without their first-hand successors")
        rows = [r for r in rows if r.rid not in CHEN_SUPERSEDED]
        rows = reroute_repulsive_f2(rows)
        rows = split_krf_bc(rows) + krf_bc_mixing()
    else:
        raise ValueError(f"unknown mode {mode!r}")
    rows = rows + upper_level_radiation(mode)
    check_network(rows)
    return rows


# ---------------------------------------------------------------------------
# optics data (the 'sigma' rows of Chen Table III + the REPO 248 nm set)
# ---------------------------------------------------------------------------
# Photo-product routing: an absorbed photon TRANSFORMS its absorber (the
# species-side twin of the photon equation's loss term, spectre 2026-07-15).
PHOTO_PRODUCTS = {
    "Kr2+": ("Kr+", "Kr"), "Ar2+": ("Ar+", "Ar"), "F-": ("F",),
    "F2": ("F", "F"), "Kr2F*": ("Kr", "Kr", "F"),
    "Ar2F*": ("Ar", "Ar", "F"),      # Kannari 1985 Table I p.4313
    "Ar**": ("Ar+", "e"), "Kr**": ("Kr+", "e"),
    "Ar*": ("Ar+", "e"), "Kr*": ("Kr+", "e"),   # Kannari Table I photoion.
    "ArKr*": ("ArKr+", "e"),                    # 0.1.10, see mode_optics
    # Ar2* photoionization at 193 nm (Sauerbrey-1987 via petrova2022
    # Table 2; 10.4 + 6.42 >= E(Ar2+) 14.5 — open by 2.3 eV)
    "Ar2*": ("Ar2+", "e"),
}


def mode_optics(mode: str) -> dict:
    """Stimulated-emission cross section, upper level, absorber inventory.

    KrF (REPO): sigma_eff = sigma_B * f_B = 5.354e-21 m2 (Tellinghuisen 1976
    2.6e-16 cm2 x f_B = 0.206) with tau_eff 25 ns — the CONSISTENT lumped
    pair (spectre rates.py). ArF (CHEN): Chen R38 sigma = 4e-20 m2, R37
    tau = 1/2.381e8; absorbers R40 (Ar** 1e-22) and R41 (F- 1e-21) verbatim.
    """
    if mode == "krf":
        return {
            # Kannari 1985 split: only B lases; sigma_B = 2.6e-16 cm2
            # (Table I [ref 46], first-hand 2026-08-04), tau_B = 6.99 ns.
            # The lumped (sigma_eff 5.354e-21, 25 ns) pair is retired with
            # the lumped species.
            "upper": "KrF*(B)", "wavelength_nm": 248.0,
            "sigma_se_m2": 2.6e-20, "tau_upper_s": 6.993e-9,
            "absorbers": {
                # F-: 5.6e-18 cm2 CONFIRMED first-hand, Kannari 1985
                # Table I p.4312 [ref 82] (2026-08-04)
                "F-": 5.6e-22,
                "Kr2+": 1.5e-21 * 0.182,  # J&H 1980 R74 x thermal scale
                "Ar2+": 1.0e-21 * 0.126,  # J&H 1980 R75 x thermal scale
                # Kannari Table I p.4313 [ref 16]: 6.0e-18 cm2 — the
                # phantom-derived EST 2.2e-19 cm2 sat x27 LOW (2026-08-04)
                "Kr2F*": 6.0e-22,
                # was MISSING: Kannari Table I p.4313 [ref 16] 1.0e-18 cm2
                "Ar2F*": 1.0e-22,
                # excited-atom photoionization, Kannari Table I p.4312/4313
                # [ref 16], first-hand 2026-08-04 (energetically valid:
                # 11.6+5.0 > IP(Ar), 9.92+5.0 > IP(Kr))
                "Ar*": 4.5e-22,
                "Kr*": 1.0e-23,
                # Ar**/Kr** photoionization: Kannari Table I p.4312
                # [ref 84] 6.0e-18 cm2 each — MISSING through 0.1.9
                # (the routing existed for the arf_chen R40 absorber,
                # the krf sigmas were never adopted); first-hand
                # 2026-08-05. J&H #76/#77 print 2.3e-18/4.5e-18
                "Ar**": 6.0e-22,
                "Kr**": 6.0e-22,
                # ArKr* photoionization: Kannari Table I p.4313 [ref 85]
                # 1.5e-17 cm2, first-hand 2026-08-05 — adoptable only
                # once ArKr+ exists (0.1.10). The scan's product print
                # is ambiguous; read as ArKr+ +e, the only channel open
                # at 248 nm (9.5+5.0 = 14.5 >= 13.5; Ar+ +Kr would need
                # 15.76). NOTE: the 0.1.9 OMISSIONS text called this
                # 'the ArKr+ photoabsorption' — species TYPO, the
                # printed absorber is the NEUTRAL ArKr*; neither
                # primary prints a sigma for the ion itself
                "ArKr*": 1.5e-21,
                # Kannari Table I p.4312 [ref 81]: F2+hv -> 2F
                # 1.2e-20 cm2 (first-hand 2026-08-05; supersedes the
                # 1.5e-24 m2 EST, which sat x1.25 above the print)
                "F2": 1.2e-24,
            },
            "fragments": ("Kr", "F"),
        }
    if mode == "arf_chen":
        return {
            "upper": "ArF*", "wavelength_nm": 193.0,
            "sigma_se_m2": 4e-20, "tau_upper_s": 1.0 / 2.381e8,
            "absorbers": {"Ar**": 1e-22, "F-": 1e-21},   # Chen R40, R41
            "fragments": ("Ar", "F"),
        }
    if mode == "arf_beam":
        # The Electra ArF anchor mode (0.1.11): NRL Orestes' PRINTED
        # 193-nm optical set (petrova2022 NRL/6720/MR--2022/1 Table 2 +
        # Appendix A4, first-hand 2026-08-05). sigma_se: their
        # vibrationally-averaged EFFECTIVE 1.0e-16 cm2 (memo p.16 Fig.10
        # + text); the B(v=0) line prints 1.9e-16 [McGeoch-2017/
        # Dunning-1978] and Lee-1989 adopted 2.75e-16 — RECORDED, the
        # effective value is what their laser model runs. tau matches
        # ARF_RAD (Kannari 2.50e8; NRL print 2.67e8+2.1e7 [Dunning] is
        # x1.15 — recorded on the row).
        return {
            "upper": "ArF*", "wavelength_nm": 193.0,
            "sigma_se_m2": 1.0e-20, "tau_upper_s": 4.0e-9,
            "absorbers": {                    # Table 2 (cm2 -> m2)
                "F2": 1.0e-25,                # Arguello-1995
                "Ar*": 3.5e-24,               # Duzy-1980
                "Ar**": 3.0e-22,              # Duzy-1980
                "Ar2*": 6.3e-23,              # Sauerbrey-1987
                "Ar2F*": 5.0e-23,             # McGeoch-2017
                "Ar2+": 1.0e-23,              # NRL estimate
                "F-": 6.5e-22,                # Mandl-1971
            },
            "fragments": ("Ar", "F"),
        }
    raise ValueError(f"unknown mode {mode!r}")
