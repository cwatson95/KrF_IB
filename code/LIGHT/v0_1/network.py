"""The V0 reaction set AS DATA, plus the optional Kr/F2 groups.

Version 0 is the reaction list in KrF_model_notes_with_Gary.rtf (group
'gary', G1..G6, numbered as in the notes) plus the three rows the notes'
own output list needs to be computable (group 'closure'):

  G1  e + Kr -> e + Kr*            (5s manifold, Biagi 'Kr(s)' 9.915 eV)
  G2  Kr* + Kr* -> Kr+ + Kr + e    (notes: 'associative ionization')
  G3  Kr* + F2 -> KrF(B) + F       (the harpoon reaction)
  G4  KrF(B) -> Kr + F + hv 248 nm (spontaneous emission, 1/tau_sp)
  G5  KrF(B) + Kr + Kr -> Kr2F* + Kr (three-body trimer 4^2Gamma; M = Kr)
  G6  e + Kr elastic               (energy loss to the gas)
  C1  e + Kr -> Kr(5p, 4d, higher) -> Kr* + heat   (prompt cascade to 5s)
  C2  e + Kr -> 2e + Kr+           (direct ionization: the 'Kr+' output)
  C3  Kr2F* -> 2Kr + F + hv ~415 nm (so trimer energy has a fate)

Optional groups (OFF in V0; switch on with Params.groups / --add) are the
rest of the pure Kr/F2 chemistry that SPECK/REBUILD carries, copied with
REBUILD 0.1.40's row ids and values (the drift test compares them). They
exist so the cost of V0's omissions can be MEASURED (run.py --omissions)
rather than argued.

Units: cm^3/s (two-body), cm^6/s (three-body), 1/s (radiative).

Rate specs (declarative, evaluated by model.py):
  ("const", k)                      constant
  ("inv_tau", param[, frac])        frac / Params.<param> (frac 1 if absent)
  ("xsec", gas, kind, product)      Maxwellian <sigma v>(Te) of one LXCat block
  ("elastic", gas)                  G6: energy loss only, no density change
  ("stepwise", E_bind)              1.56e-7 Te^0.71 exp(-E/Te)  (REBUILD)
  ("dr_sqrt", k0)                   k0 / sqrt(max(Te, 0.02))    (REBUILD)
  ("te_pow", k0, Te0, n)            k0 (Te/Te0)^n, Te held above 0.02 eV;
                                    a captured electron takes (3/2 + n) Te
  ("flannery", a, scale[, n3])      [a N/(N^2 + 1.7)] scale, N in amagat,
                                    times the Debye screening S (REBUILD).
                                    It is the harmonic mean of a three-body
                                    part a N/1.7 and a mobility part a/N;
                                    n3 puts (Tg/300)^n3 on the three-body
                                    part (v0_1: n3 = -3, Thomson theory)
  ("superelastic", gas, product, g_ratio)
                                    de-excitation from the Klein-Rosseland
                                    cross section of an LXCat excitation
                                    (= k_fwd g_ratio exp(E_thr/Te))
  ("balance", k_fwd, K_ref, upper, lower)
                                    reverse of a two-state exchange by
                                    detailed balance: k_fwd / K(Tg), with
                                    K(Tg) = K_ref exp(dE/k (1/Tg - 1/300)),
                                    dE = E_form(upper) - E_form(lower), the
                                    gap the energy ledger books
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from .species import ELECTRON, atoms, charge, e_form


@dataclass(frozen=True)
class Row:
    rid: str
    reactants: tuple
    products: tuple
    rate: tuple
    group: str
    cite: str
    tg_pow: float = 0.0            # k *= (Tg/300)^tg_pow
    hv_eV: float = None            # radiative rows: the emitted quantum
    band: str = None               # emission band label (radiative rows)
    e_cost_eV: float = None        # xsec rows: the electron's energy loss
    energy_to: str = None          # 'electrons': exothermic release -> U_e
    rebuild: str = None            # REBUILD 0.1.40 row id this is copied from

    @property
    def release_eV(self) -> float:
        """sum E_form(reactants) - sum E_form(products); > 0 exothermic."""
        return (sum(e_form(s) for s in self.reactants)
                - sum(e_form(s) for s in self.products))

    @property
    def ne_delta(self) -> int:
        return self.products.count(ELECTRON) - self.reactants.count(ELECTRON)

    @property
    def heavy(self) -> tuple:
        return tuple(s for s in self.reactants if s != ELECTRON)

    @property
    def n_electrons_in(self) -> int:
        return self.reactants.count(ELECTRON)

    @property
    def is_rad(self) -> bool:
        return self.hv_eV is not None


class NetworkError(ValueError):
    pass


def energy_routing(r: Row) -> dict:
    """Per-event energy flows [eV] (the REBUILD routing rule):

    * radiative: the quantum escapes, the rest (lower-state repulsion) heats;
    * electron-impact with an LXCat threshold: the electron pays the
      threshold, the products store dE_form, the difference heats (C1's
      cascade, F2 vibration);
    * exothermic: to the electron gas when the row ejects a new electron from
      heavy reactants (Penning) or says energy_to='electrons'
      (superelastic), otherwise to the gas;
    * endothermic without a threshold: from the electron gas if an electron
      collides, otherwise from the gas.
    The captured electron's (3/2)Te (DR, attachment) moves electron gas ->
    gas at run time ('capture')."""
    dE = r.release_eV
    out = {"from_e": 0.0, "to_e": 0.0, "heat": 0.0, "rad": 0.0,
           "capture": max(-r.ne_delta, 0)}
    if r.rate[0] == "elastic":
        return out
    if r.is_rad:
        out["rad"] = r.hv_eV
        out["heat"] = dE - r.hv_eV
    elif r.e_cost_eV is not None:
        out["from_e"] = r.e_cost_eV
        out["heat"] = r.e_cost_eV + dE
    elif dE >= 0.0:
        penning = r.n_electrons_in == 0 and r.ne_delta > 0
        if r.energy_to == "electrons" or penning:
            out["to_e"] = dE
        else:
            out["heat"] = dE
    elif r.n_electrons_in:
        out["from_e"] = -dE
    else:
        out["heat"] = dE                  # endothermic: drawn from the gas
    return out


def check_row(r: Row) -> None:
    if not r.cite:
        raise NetworkError(f"{r.rid}: no citation")
    lhs, rhs = {}, {}
    for side, acc in ((r.reactants, lhs), (r.products, rhs)):
        for s in side:
            for el, n in atoms(s).items():
                acc[el] = acc.get(el, 0) + n
    if lhs != rhs:
        raise NetworkError(f"{r.rid}: elements {lhs} -> {rhs}")
    q_in = sum(charge(s) for s in r.reactants)
    if q_in != sum(charge(s) for s in r.products):
        raise NetworkError(f"{r.rid}: charge not balanced")
    route = energy_routing(r)
    if route["heat"] < -1e-9 and (r.is_rad or r.e_cost_eV is not None):
        raise NetworkError(f"{r.rid}: quantum/threshold exceeds the energy "
                           f"available (heat {route['heat']:.4f} eV)")


# ---------------------------------------------------------------------------
# Version 0: the notes' six reactions + three closures
# ---------------------------------------------------------------------------
def v0_rows() -> list:
    return [
        Row("G1", ("e", "Kr"), ("e", "Kr*"),
            ("xsec", "Kr", "excitation", "Kr(s)(9.915eV)"), "gary",
            "notes (1) e + Kr -> Kr* 5s[3/2]. Biagi-v7.1 'Kr(s)' = "
            "the whole 5s manifold (5s[3/2]1,2 and 5s'[1/2]0,1) lumped at "
            "9.915 eV; k = <sigma v> over a Maxwellian (notes: Q = k n_e[Kr])",
            e_cost_eV=9.915),
        Row("G2", ("Kr*", "Kr*"), ("Kr+", "Kr", "e"), ("const", 5.0e-10),
            "gary", "notes (2) Kr* + Kr* 'associative ionization'. Rate and "
            "products: Kannari 1985 Table I Penning block, 5.0e-10 cm3/s "
            "[ref 86] (REBUILD PEN_KR, first-hand 2026-08-05), which prints "
            "Kr+ + Kr + e; an associative Kr2+ + e branch would move 1.15 "
            "eV per event from 'ion' to 'electron energy', nothing else",
            rebuild="PEN_KR"),
        Row("G3", ("Kr*", "F2"), ("KrF(B)", "F"), ("const", 8.1e-10),
            "gary", "notes (3) Kr*(5s) + F2 -> KrF(B) (the harpoon). Kannari "
            "1985 Table I p.4311, 8.1e-10 cm3/s [52,59,93] (REBUILD K25, "
            "first-hand 2026-08-04); all product booked to B (V0 has no C). "
            "Untraceable: Kannari's ref. 52, SRI MP-76-99, gives 7.2e-10 "
            "(docs/krf_data/FORMATION_LINKS_SEARCH.md). The basic network "
            "uses the first-hand mean, 8.47e-10", rebuild="K25"),
        Row("G4", ("KrF(B)",), ("Kr", "F"), ("inv_tau", "tau_sp_s"),
            "gary", "notes (4) KrF* -> Kr + F + hv, spontaneous emission. "
            "Default tau_sp = 6.99 ns: Kannari 1985 Table I p.4312, "
            "1.43e8 1/s [47,48] (REBUILD URAD_B). The notes write "
            "'tau_sp^-1 ~ 4 ns'; set Params.tau_sp_s = 4e-9 to use it",
            hv_eV=5.00, band="248 nm (KrF B-X)", rebuild="URAD_B"),
        Row("G5", ("KrF(B)", "Kr", "Kr"), ("Kr2F*", "Kr"), ("const", 6.0e-31),
            "gary", "notes (5) KrF(B) + Kr + M -> Kr2F* (4^2Gamma) (the notes "
            "call it harpooning; it is three-body trimer formation). M = Kr "
            "in a Kr/F2 fill (F2 < 1% of M, no measured F2 third body). "
            "Xu/Gadomski/Setser JCP 99, 2591 (1993) Table II, "
            "(6.0 +- 2.0)e-31 cm6/s at 300 K (REBUILD K34); Tg form "
            "(Tg/300)^-3.32 = Shui's R = Kr calculation (x1/10 from 300 to "
            "600 K, Klimek & Hsia 1981 measured it to 15 %)",
            tg_pow=-3.321928094887362, rebuild="K34"),
        Row("G6", ("e", "Kr"), ("e", "Kr"), ("elastic", "Kr"),
            "gary", "notes (6) e + Kr elastic. Biagi-v7.1 momentum-transfer "
            "set; power to the gas (2m/M)(1 - Tg/Te) n_e N <sigma_m v eps>"),
        Row("C1a", ("e", "Kr"), ("e", "Kr*"),
            ("xsec", "Kr", "excitation", "Kr(p)(11.30eV)"), "closure",
            "Biagi-v7.1 5p group (11.30 eV). Not in the notes: "
            "without it the electrons lose ~half their inelastic power "
            "nowhere. Cascade to 5s taken prompt: p->s collisional "
            "relaxation, Kannari 'estimated' 1e-10 cm3/s (REBUILD KPS_Kr), "
            "is <= 4 ns at >= 0.1 atm; the 1.385 eV excess heats the gas",
            e_cost_eV=11.30),
        Row("C1b", ("e", "Kr"), ("e", "Kr*"),
            ("xsec", "Kr", "excitation", "Kr(d+p)(11.998eV)"), "closure",
            "Biagi-v7.1 d+p group (11.998 eV); prompt cascade to 5s as C1a",
            e_cost_eV=11.998),
        Row("C1c", ("e", "Kr"), ("e", "Kr*"),
            ("xsec", "Kr", "excitation", "Kr*(12.75eV)"), "closure",
            "Biagi-v7.1 higher group (12.75 eV); prompt cascade "
            "to 5s as C1a", e_cost_eV=12.75),
        Row("C2", ("e", "Kr"), ("e", "e", "Kr+"),
            ("xsec", "Kr", "ionization", "Kr^+"), "closure",
            "Biagi-v7.1 ionization (13.996 eV): the notes' 'Kr+' "
            "output needs electron-impact ionization", e_cost_eV=13.996),
        Row("C3", ("Kr2F*",), ("Kr", "Kr", "F"), ("const", 4.5e6), "closure",
            "Kr2F* -> 2Kr + F + hv: Xu/Gadomski/Setser 1993 p.2597, "
            "tau = 220 +- 20 ns (REBUILD K169); gives G5's trimer a fate",
            hv_eV=3.0, band="~415 nm (Kr2F*)", rebuild="K169"),
    ]


# ---------------------------------------------------------------------------
# Optional groups: the rest of SPECK's pure Kr/F2 chemistry (OFF in V0)
# ---------------------------------------------------------------------------
def optional_rows() -> list:
    return [
        # --- quench: collisional losses of the emitters -------------------
        Row("K33", ("KrF(B)", "Kr"), ("Kr", "Kr", "F"), ("const", 4.0e-12),
            "quench", "Xu/Gadomski/Setser 1993 Table I + p.2596, "
            "recommended (4.0 +- 0.5)e-12 cm3/s (REBUILD K33)", rebuild="K33"),
        Row("K32", ("KrF(B)", "F2"), ("Kr", "F", "F", "F"), ("const", 3.0e-10),
            "quench", "Xu/Gadomski/Setser 1993: (3.0 +- 0.3)e-10 cm3/s "
            "(REBUILD K32)", rebuild="K32"),
        Row("K67", ("KrF(B)", "e"), ("Kr", "F", "e"), ("const", 4.5e-8),
            "quench", "Kushner-Hanson-Schneider APL 55, 2482 (1989) "
            "(REBUILD K67)", rebuild="K67"),
        Row("K125", ("Kr2F*", "F2"), ("Kr", "Kr", "F", "F2"),
            ("const", 1.5e-10), "quench", "Kannari 1985 Table I p.4311, "
            "1.5e-10 cm3/s [89,94] (REBUILD K125)", rebuild="K125"),
        Row("K_2Fe", ("Kr2F*", "e"), ("Kr", "Kr", "F", "e"), ("const", 1.0e-7),
            "quench", "Kannari 1985 Table I p.4312, 'estimated' there "
            "(REBUILD K_2Fe)", rebuild="K_2Fe"),
        # --- kr2: the Kr2* excimer channel --------------------------------
        Row("AK8_Kr", ("Kr*", "Kr", "Kr"), ("Kr2*", "Kr"), ("const", 2.7e-32),
            "kr2", "Kr-third-body analog of Kannari's Kr* + Kr + Ar 2.7e-32 "
            "cm6/s, EST (REBUILD AK8_Kr)", rebuild="AK8_Kr"),
        Row("KR2_RAD", ("Kr2*",), ("Kr", "Kr"), ("const", 8.0e7), "kr2",
            "Kannari 1985 Table I p.4312, 8.0e7 1/s [88] (REBUILD KR2_RAD; "
            "J&H 1980 prints 3.3e6, recorded conflict)", hv_eV=8.43,
            band="147 nm (Kr2* VUV)", rebuild="KR2_RAD"),
        Row("NK7", ("Kr2*", "F2"), ("KrF(B)", "Kr", "F"), ("const", 3.0e-10),
            "kr2", "Kannari 1985 Table I p.4311 [61] (REBUILD NK7)",
            rebuild="NK7"),
        Row("K2F_BR", ("Kr2*", "F2"), ("Kr2F*", "F"), ("const", 3.0e-10),
            "kr2", "Kannari 1985 Table I p.4311 [52]; J&H #20 (REBUILD "
            "K2F_BR)", rebuild="K2F_BR"),
        Row("K2F_F", ("Kr2*", "F"), ("KrF(B)", "Kr"), ("const", 3.0e-10),
            "kr2", "Kannari 1985 Table I p.4311; J&H #21 (REBUILD K2F_F). "
            "An estimate: Kannari's ref. 52, SRI MP-76-99, tabulates no such "
            "rate and only interprets a slow KrF* component at very low F2 "
            "as this reaction (p.42). No measurement exists "
            "(docs/krf_data/FORMATION_LINKS_SEARCH.md). Re-read 2026-09-29: "
            "the report's model leaves Ar2* and Kr2* + F out (p.45), and its "
            "one number is conditional and for argon, 'if the reaction Ar2* "
            "+ F -> ArF* + Ar is sufficiently rapid (3-10 x 10^-10)' it could "
            "explain an Ar2F* yield trend (p.65). 3.0e-10 is the low end of "
            "that range, carried to Kr, and is no upper bound. Chen & Payne, "
            "IEEE JQE 15, 149 (1979) Fig. 31 assign a slow KrF* component to "
            "ion-ion recombination instead. The owner's rate tables and Ong "
            "1996 carry the same 3.0e-10, from the same report or from "
            "Kannari", rebuild="K2F_F"),
        Row("PEN_KR2", ("Kr2*", "Kr2*"), ("Kr2+", "Kr", "Kr", "e"),
            ("const", 5.0e-10), "kr2", "Kannari 1985 Table I p.4312 [87] "
            "(REBUILD PEN_KR2). Ref. 87 is Lorents & Olson, SRI semiannual "
            "technical report N00014-72-C-0557 (1973, unpublished), and "
            "Kannari gives Ar2* + Ar2* the same 5.0e-10 from it (read "
            "first-hand 2026-09-29). No measurement of this rate is among "
            "the sources held", rebuild="PEN_KR2"),
        # --- ions: conversion and recombination ---------------------------
        Row("K15", ("Kr+", "Kr", "Kr"), ("Kr2+", "Kr"), ("const", 2.5e-31),
            "ions", "J&H 1980 #41, Kannari 1985, Boichenko 2000 (REBUILD K15)",
            tg_pow=-0.75, rebuild="K15"),
        Row("R11p", ("e", "Kr2+"), ("Kr*", "Kr"), ("dr_sqrt", 2.24e-7),
            "ions", "Kannari 1985 Table I p.4312 DR block, 1.39e-6 cm3/s at "
            "300 K on the Te^-1/2 template (REBUILD R11p)", rebuild="R11p"),
        # --- electron: electron collisions with Kr* -----------------------
        Row("KSTEP", ("e", "Kr*"), ("Kr+", "e", "e"), ("stepwise", 4.09),
            "electron", "Chen R9 template, E_bind 4.09 eV, EST; J&H #48 "
            "discharge value 1.5e-8 cm3/s agrees (REBUILD KSTEP)",
            rebuild="KSTEP"),
        Row("SUPER", ("e", "Kr*"), ("e", "Kr"),
            ("superelastic", "Kr", "Kr(s)(9.915eV)", 1 / 12),
            "electron", "detailed balance of G1, g(ground)/g(5s manifold) = "
            "1/12 (statistical over the four 5s levels), EST — LIGHT's own "
            "row, REBUILD carries superelastics inside its swarm tables",
            energy_to="electrons"),
        # --- bc: the KrF B<->C manifold -----------------------------------
        Row("BC_Kr", ("KrF(B)", "Kr"), ("KrF(C)", "Kr"), ("const", 5.0e-10),
            "bc", "Kr-partner analog of Kannari's B<->C mixing 5e-10, EST "
            "(REBUILD BC_Kr)", rebuild="BC_Kr"),
        Row("CB_Kr", ("KrF(C)", "Kr"), ("KrF(B)", "Kr"), ("const", 3.41e-10),
            "bc", "detailed balance of BC_Kr, exp(-80 cm^-1/kT) at 300 K "
            "(REBUILD CB_Kr)", rebuild="CB_Kr"),
        Row("URAD_C", ("KrF(C)",), ("Kr", "F"), ("const", 1.33e7), "bc",
            "Kannari 1985 Table I p.4312, C->A 1.33e7 1/s, booked at 4.99 "
            "eV as REBUILD does (REBUILD URAD_C)", hv_eV=4.99,
            band="KrF C-A", rebuild="URAD_C"),
        # --- f2e: electron collisions with F2 and the ions they make ------
        Row("ATT", ("e", "F2"), ("F-", "F"),
            ("xsec", "F2", "attachment", None), "f2e",
            "dissociative attachment, SIGLO set (the F2 file REBUILD reads)"),
        Row("DIS_a3", ("e", "F2"), ("e", "F", "F"),
            ("xsec", "F2", "excitation", "F2(a3)(3.16eV)"), "f2e",
            "SIGLO a3Pi_u (repulsive: excitation = dissociation, REBUILD "
            "0.1.14 reroute)", e_cost_eV=3.16),
        Row("DIS_A1", ("e", "F2"), ("e", "F", "F"),
            ("xsec", "F2", "excitation", "F2(A1)(4.34eV)"), "f2e",
            "SIGLO A1Pi_u (repulsive, as DIS_a3)", e_cost_eV=4.34),
        Row("VIB1", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V1)(0.111eV)"), "f2e",
            "SIGLO v=1, the quantum relaxes to heat", e_cost_eV=0.111),
        Row("VIB2", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V2)(0.219eV)"), "f2e",
            "SIGLO v=2", e_cost_eV=0.219),
        Row("VIB3", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V3)(0.324eV)"), "f2e",
            "SIGLO v=3", e_cost_eV=0.324),
        Row("VIB4", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V4)(0.421eV)"), "f2e",
            "SIGLO v=4", e_cost_eV=0.421),
        Row("K26", ("Kr+", "F-"), ("KrF(B)",), ("flannery", 2.4, 3.417e-6),
            "f2e", "ion-ion, Flannery turnover (Kannari 1985 Sec II.B) with "
            "Debye screening (REBUILD K26)", rebuild="K26"),
        Row("K108", ("Kr2+", "F-"), ("KrF(B)", "Kr"),
            ("flannery", 2.65, 3.096e-6), "f2e", "ion-ion, same form "
            "(REBUILD K108)", rebuild="K108"),
    ]


GROUPS = ("quench", "kr2", "ions", "electron", "bc", "f2e")


def c_state_copies(rows) -> list:
    """With the 'bc' group on, every non-radiative row that consumes KrF(B)
    also consumes KrF(C) at the same k (REBUILD's split_krf_bc rule: the
    quench constants are manifold-effective). Without this the C state
    would be a quench-free reservoir."""
    out = []
    for r in rows:
        if "KrF(B)" in r.reactants and not r.is_rad and r.group != "bc":
            rea = tuple("KrF(C)" if s == "KrF(B)" else s for s in r.reactants)
            out.append(replace(r, rid=r.rid + "_C", reactants=rea, group="bc",
                               cite=r.cite + " ; applied to KrF(C) at the "
                               "same k (REBUILD split_krf_bc)"))
    return out


# ---------------------------------------------------------------------------
# The basic network: the six notes' processes plus the KrF data that pin them
# ---------------------------------------------------------------------------
# KrF(B,C) as two coupled levels (Gadomski, Xu, Setser & Nelson, Chem.
# Phys. Lett. 189, 153 (1992) eq. 5), with the constants of ML's joint fit
# (ML/README.md, 2026-09-27) to the three measurement sets that agree with
# one another:
#   Eden, Waynant, Searles & Burnham 1978: KrF(B) alone, decay against Kr
#     and F2 (APL 32, 733 Figs. 3 and 5);
#   Gadomski et al. 1992 Figs. 1-2: the Kansas State B/C decays against F2
#     and Kr;
#   Yu & Setser 1990: the equilibrium C/B emission ratio at 300 K (J. Phys.
#     Chem. 94, 2934, eq. 18).
# Hakuta et al. 1990 (tau_C 60 ns, k_BC 1.2e-10) disagrees with all three
# (chi2 12 per point) and is held out. B and C carry their own quench
# constants. The MCMC gives k_BC 9.8e-12 [4.6e-12, 1.7e-11], K_eq 1.07
# [0.84, 1.38] and tau_C 50 [41, 63] ns (68 %). The B + Kr two-body quench
# sits at the fit's lower bound: Eden's two-body term is B -> C transfer,
# and no other B + Kr loss is seen. tau_B stays 6.99 ns (Params; the fit
# gives 7.03 under a 6.9 +- 0.3 ns prior). The set predicts Tahoma's
# held-out Isat, 1.49 against 1.4 +- 0.1 MW/cm2. Adopted 2026-09-28.
# K_eq follows the gas temperature by detailed balance ("balance" rows).
KRF_K_BC = 1.05e-11
KRF_K_EQ = 1.10
KRF_TAU_C = 50.0e-9
KRF_K3_KR = 1.90e-30        # KrF(B) + 2Kr -> Kr2F* + Kr, 300 K
KRF_K3_KR_C = 3.79e-31      # KrF(C) + 2Kr
KRF_K_KR = 1.0e-14          # KrF(B) + Kr, the fit's lower bound
KRF_K_KR_C = 8.37e-13       # KrF(C) + Kr
KRF_K_F2 = 4.75e-10         # KrF(B) + F2
KRF_K_F2_C = 3.04e-10       # KrF(C) + F2
KRF_K_E = 6.6e-8            # KrF(B,C) + e, Trainor & Jacob re-read (6.58)
ML_FIT = ("ML joint fit 2026-09-27 to Eden 1978, Gadomski 1992 and Yu & "
          "Setser 1990 (ML/README.md, Results for the KrF(B,C) block)")
# Kr* + F2 forms KrF in B and in C (2026-09-28, step 2): k_B/k_C is about
# 1.6 for Kr(3P1) + F2 (Yu & Setser 1990, Sec. III.A, "1.6 is a good
# estimate"), and Setser's 1979 report finds the same B/C for 3P1 and 3P2.
# With the fit's slow B-C mixing the C share matters for the gain.
# The Kr* + F2 rate is the mean of the three first-hand measurements, 7.2
# (Velazco, Kolts & Setser 1976), 8.8 +- 0.5 and 9.4 +- 1.0 (Ueda 1984
# Table V) e-10 cm3/s (2026-09-28, step 4). v0 keeps Kannari's 8.1e-10,
# which he and Ueda credit to SRI MP-76-99, a report that gives 7.2e-10.
KRF_K_HARPOON = 8.47e-10    # Kr* + F2, cm3/s
KRF_F_C_HARPOON = 0.38      # share of the harpoon KrF formed in C
# The KrF* yield of Kr* + F2 (2026-09-28, from the papers the owner
# supplied, read first-hand). Velazco, Kolts & Setser, J. Chem. Phys. 65,
# 3468 (1976) measured the KrF* formation rate constant absolutely, 'the
# measured value of the KrF* formation rate constant (65 x 10^-11 cm3
# molecule-1 sec-1)' (p.3476), 0.88 of the KrCl* one (Table III, +-10 %),
# against the total Kr(3P2) + F2 quenching rate constant 72 x 10^-11
# (Table II). The yield is 65/72 = 0.90. The paper rounds it to unity
# (Table IV, 'The uncertainties are estimated as +-0.20'), so the measured
# range is 0.70 to 1.0. Velazco's dissertation (1977, p.148) quotes the
# Kr* + Cl2 reference at 0.90. The rest of Kr* + F2 dissociates (G3X).
KRF_Y_HARPOON = 0.90        # KrF* per Kr* + F2 reaction
# v0_1 (2026-09-28): the B -> A branch of KrF(B) and the printed band
# positions (the provenance is in the rows that use them)
KRF_BA_OVER_BX = 0.035      # A(B-A)/A(B-X), Dunning & Hay's calculation
HC_EV_NM = 1239.84193       # eV nm, photon energy = HC_EV_NM / lambda
KR_CONV_300 = 2.4e-31       # Kr+ + 2Kr -> Kr2+ + Kr, cm6/s at 300 K
FF_KR_300 = 2.8e-34         # F + F + Kr -> F2 + Kr, cm6/s
BASIC_GROUPS = ("kr2", "electron")     # groups that add to it (v0_1: the
                                        # ions and f2e rows are basic rows)
TRIMER_TG = ("Tg form (Tg/300)^-3.32, a two-point fit to Shui's classical phase-space calculation for R = Kr (V. H. Shui, APL 34, 203 (1979), computed 300-800 K; the power law stays within +16/-11 % of his curve there, Avco ADA069053 Fig. 22). Measured temperature dependence: 315-420 K, Klimek & Hsia, JAP 52, 5361 (1981), agreeing with Shui to 15 % (Ar-rich mixtures, Avco report 10.2172/5268869 Sec. IV.F). No data above 425 K. Salesky & Kimura (APL 47, 774 (1985), 48, 897 (1986), abstracts) inferred a weaker T^-1/2 for pure Kr/F2 by fitting a laser model to efficiencies, so it bounds the law and never sets it")


def basic_rows() -> list:
    """The KrF(B,C) network: v0's rows with KrF(B) given its partner state
    KrF(C), and the two-level constants of ML's joint fit to the measured
    KrF(B), KrF(B,C) and C/B data (the header above; provenance in
    ML/README.md and docs/krf_data/BC_AND_ELECTRON_QUENCH.md). tau_B is
    Kannari's 6.99 ns. The electron quench is Trainor & Jacob's measured
    product re-read with this network's own Stern-Volmer lifetime."""
    v0 = {r.rid: r for r in v0_rows()}
    v0["G5"] = Row(
        "G5", ("KrF(B)", "Kr", "Kr"), ("Kr2F*", "Kr"), ("const", KRF_K3_KR),
        "basic", "KrF(B) + 2Kr -> Kr2F* + Kr, 1.90e-30 cm6/s at 300 K, "
        + ML_FIT + ". " + TRIMER_TG + ". Xu 1993 Table II printed 6.0 "
        "+- 2.0e-31, read with K_eq = 5.0 and fast B-C mixing, which the "
        "joint fit rejects",
        tg_pow=-3.321928094887362)
    v0["G3"] = Row(
        "G3", ("Kr*", "F2"), ("KrF(B)", "F"),
        ("const", KRF_K_HARPOON * KRF_Y_HARPOON * (1.0 - KRF_F_C_HARPOON)),
        "basic",
        "Kr* + F2 -> KrF(B) + F: 62 % of the 90 % of 8.47e-10 cm3/s that "
        "forms KrF*. The rate is the mean of the first-hand 7.2 (Velazco, "
        "Kolts & Setser 1976 via COO-2807-2 Table 3), 8.8 +- 0.5 and 9.4 "
        "+- 1.0e-10 (Ueda 1984 Table V). The KrF* yield 0.90 is the "
        "measured KrF* formation rate constant, 65e-11, over the total "
        "quenching, 72e-11 (Velazco, Kolts & Setser, JCP 65, 3468 (1976) "
        "p.3476 and Table II), which the paper rounds to unity +- 0.20 "
        "(Table IV). The 62/38 B/C split is Yu & Setser's k_B/k_C 1.6 "
        "(docs/krf_data/FORMATION_LINKS_SEARCH.md)")
    v0["G4"] = Row(
        "G4", ("KrF(B)",), ("Kr", "F"),
        ("inv_tau", "tau_sp_s", 1.0 / (1.0 + KRF_BA_OVER_BX)), "basic",
        "KrF(B) -> Kr + F + hv, 248 nm (B-X): the 96.6 % share of B's "
        "radiative decay, tau_B = 6.99 ns in all (Kannari 1985 Table I "
        "p.4312, 1.43e8 1/s; Eden 1978 measured 6.8 +- 0.2 ns). The rest "
        "is B -> A (G4A)", hv_eV=5.00, band="248 nm (KrF B-X)")
    v0["C3"] = Row(
        "C3", ("Kr2F*",), ("Kr", "Kr", "F"), ("const", 4.5e6), "basic",
        "Kr2F* -> 2Kr + F + hv: Xu/Gadomski/Setser 1993 p.2597, tau = 220 "
        "+- 20 ns (REBUILD K169). Band maximum 420 +- 35 nm, printed "
        "(Huestis, Marowsky & Tittel 1984 Table 6.1 p.185)",
        hv_eV=round(HC_EV_NM / 420.0, 2), band="420 nm (Kr2F*)")
    rows = [v0[k] for k in ("G1", "G2", "G3", "G4", "G5", "G6",
                            "C1a", "C1b", "C1c", "C2", "C3")]
    gad = ("Gadomski, Xu, Setser & Nelson, CPL 189, 153 (1992) eq. 5 and "
           "p.156")
    rows += [
        Row("G3C", ("Kr*", "F2"), ("KrF(C)", "F"),
            ("const", KRF_K_HARPOON * KRF_Y_HARPOON * KRF_F_C_HARPOON),
            "basic",
            "Kr* + F2 -> KrF(C) + F: 38 % of the KrF* from 8.47e-10 cm3/s "
            "(yield 0.90, see G3). Yu & Setser, J. "
            "Phys. Chem. 94, 2934 (1990) Sec. III.A: 'the true kB/kC ratio "
            "will be somewhat larger than the apparent values, 1.3 +- 0.1, "
            "and 1.6 is a good estimate' (Kr(3P1) + F2, nascent, a few "
            "Torr; read first-hand). C -> B transfer is fast in high v, so "
            "the relaxed split at 1 atm may lie between this and all-B"),
        Row("G3X", ("Kr*", "F2"), ("Kr", "F", "F"),
            ("const", KRF_K_HARPOON * (1.0 - KRF_Y_HARPOON)), "basic",
            "Kr* + F2 without KrF*: the 10 % the measured yield (0.90, see "
            "G3) leaves. The products are not measured, and 'no other "
            "emissions were observed from the F2 reactions' (Velazco, Kolts "
            "& Setser 1976 p.3476), so dissociation to Kr + 2F, its energy "
            "to heat, is the assignment"),
        Row("G4C", ("KrF(C)",), ("Kr", "F"), ("const", 1.0 / KRF_TAU_C),
            "basic", "KrF(C) -> Kr + F + hv (C-A). tau_C = 50 ns [41, 63], "
            + ML_FIT + ". The one measurement, Hakuta, Fujino, Nakayama & "
            "Takuma, JJAP 29, 1194 (1990) Table I, 60 +- 10 ns, lies inside "
            "that interval; " + gad + " assumed Hay & Dunning's calculated "
            "75 ns. C -> X is negligible (Yu & Setser 1990 p.2941). Band "
            "maximum 275 nm, printed (Huestis, Marowsky & Tittel, Excimer "
            "Lasers 2nd ed. 1984, Table 6.1 p.185; Hay & Dunning, JCP 66, "
            "1306 (1977) abstract), FWHM about 25 nm at 300 K (Yu & Setser "
            "1990 p.2937)",
            hv_eV=round(HC_EV_NM / 275.0, 2), band="275 nm (KrF C-A)"),
        Row("G4A", ("KrF(B)",), ("Kr", "F"),
            ("inv_tau", "tau_sp_s", KRF_BA_OVER_BX / (1.0 + KRF_BA_OVER_BX)),
            "basic", "KrF(B) -> Kr + F + hv (B-A), the "
            "3.4 % branch of B's radiative decay: A(B-A)/A(B-X) = 0.035, "
            "calculated (Dunning & Hay 1978, quoted by Yu & Setser, JPC 94, "
            "2934 (1990) p.2941). 272 nm, a theoretical position under the "
            "C-A band (Huestis, Marowsky & Tittel 1984 Table 6.1 p.185, "
            "'Theoretical values in parenthesis'). It shares tau_B with "
            "G4, so B's total radiative rate is unchanged",
            hv_eV=round(HC_EV_NM / 272.0, 2), band="272 nm (KrF B-A)"),
        Row("BC", ("KrF(B)", "Kr"), ("KrF(C)", "Kr"), ("const", KRF_K_BC),
            "basic", "B -> C mixing by Kr, k_BC = 1.05e-11 cm3/s [4.6e-12, "
            "1.7e-11], " + ML_FIT + ". Slow mixing, as Eden's B-only decays "
            "demand (at most about 5e-11) and as Yu & Setser describe for "
            "low v. Hakuta et al. 1990 measured 1.2 +- 0.3e-10, which "
            "conflicts with the three fitted sets (held out)"),
        Row("CB", ("KrF(C)", "Kr"), ("KrF(B)", "Kr"),
            ("balance", KRF_K_BC, KRF_K_EQ, "KrF(B)", "KrF(C)"), "basic",
            "C -> B by detailed balance, k_BC / K_eq(Tg). K_eq = 1.10 at 300 "
            "K [0.84, 1.38], " + ML_FIT + " (Yu & Setser's measured I_r "
            "0.16 +- 0.02 pins it). K_eq(Tg) uses the B-C gap the ledger "
            "books, Yu & Setser's measured 120 +- 50 cm^-1 between 300 and "
            "230 K (species.py). " + gad + " fitted 5.0 +- 1.0 with fast "
            "mixing"),
        Row("G5C", ("KrF(C)", "Kr", "Kr"), ("Kr2F*", "Kr"),
            ("const", KRF_K3_KR_C), "basic",
            "KrF(C) + 2Kr -> Kr2F* + Kr, 3.79e-31 cm6/s at 300 K, " + ML_FIT
            + ". Tg form as G5", tg_pow=-3.321928094887362),
        Row("QF2", ("KrF(B)", "F2"), ("Kr", "F", "F", "F"),
            ("const", KRF_K_F2), "basic",
            "KrF(B) + F2, 4.75e-10 cm3/s, " + ML_FIT + " (Eden's B-only "
            "decays against F2 pin it). Xu, Gadomski & Setser JCP 99, 2591 "
            "(1993) Table I printed 3.0 +- 0.3e-10 for the B/C manifold"),
        Row("QF2C", ("KrF(C)", "F2"), ("Kr", "F", "F", "F"),
            ("const", KRF_K_F2_C), "basic",
            "KrF(C) + F2, 3.04e-10 cm3/s, " + ML_FIT),
        Row("QKR", ("KrF(B)", "Kr"), ("Kr", "Kr", "F"), ("const", KRF_K_KR),
            "basic", "KrF(B) + Kr, 1.0e-14 cm3/s, the lower bound of " + ML_FIT
            + ": Eden's two-body term is B -> C transfer, and no other B + "
            "Kr loss is seen. Xu 1993 printed 4.0 +- 1.0e-12 for the "
            "manifold with fast mixing"),
        Row("QKRC", ("KrF(C)", "Kr"), ("Kr", "Kr", "F"),
            ("const", KRF_K_KR_C), "basic",
            "KrF(C) + Kr, 8.37e-13 cm3/s, " + ML_FIT),
        Row("QE", ("KrF(B)", "e"), ("Kr", "F", "e"), ("const", KRF_K_E),
            "basic", "KrF* + e, 6.6e-8 cm3/s: Trainor & Jacob's measured "
            "product k_e tau / k_att = 2.9e-7 s (Kushner, Hanson & "
            "Schneider APL 55, 2482 (1989), recombination factor 0.62) "
            "re-read with this network's two-level Stern-Volmer lifetime "
            "at 190 Torr (docs/krf_data/electron_quench.py, ML forward."
            "with_rederived_k_e, and a test)"),
        Row("QEC", ("KrF(C)", "e"), ("Kr", "F", "e"), ("const", KRF_K_E),
            "basic", "as QE for the C state"),
        Row("Q2F", ("Kr2F*", "F2"), ("Kr", "Kr", "F", "F2"),
            ("const", 1.2e-10),
            "basic", "Kr2F* + F2, (1.2 +- 0.2)e-10 cm3/s measured, Xu 1993 "
            "Table I p.2594 (REBUILD keeps Kannari's 1.5e-10)"),
    ]
    return rows + basic_charge_rows()


def basic_charge_rows() -> list:
    """v0_1: the charge chemistry and fluorine recombination belong to the
    basic network (2026-09-28). No electrons are injected, so n_e follows
    ionization, attachment and recombination (Params.ne_mode 'auto')."""
    opt = {r.rid: r for r in optional_rows()}
    rows = [
        Row("K15", ("Kr+", "Kr", "Kr"), ("Kr2+", "Kr"), ("const", KR_CONV_300),
            "basic", "Kr+ + 2Kr -> Kr2+ + Kr, 2.4e-31 cm6/s at 300 K, the "
            "mean of three measurements: 2.3e-31 (Bhattacharya, PRA 1, 1196 "
            "(1970)), (2.2 +- 0.2)e-31 (Grossl et al., JCP 74, 1728 (1981)) "
            "and 342 +- 15 Torr^-2 s^-1 (Smith, Dean & Plumb, J. Phys. B 5, "
            "2134 (1972)), which is 2.74e-31 with Torr taken at 273 K as "
            "Brau's tables read it (3.30e-31 at 300 K); abstracts. Tg law "
            "(Tg/300)^-0.75: the classical theory of Perez-Rios & Greene, "
            "PRA 98, 062707 (2018), 'k3 prop. to T^-3/4' for rare-gas ion "
            "three-body recombination (read first-hand); Smith et al. saw "
            "the rate fall by 510 K. No measurement above 510 K, and "
            "Kr2+ + Kr redissociation (D0 1.15 eV) is not in the model",
            tg_pow=-0.75),
        Row("R11p", ("e", "Kr2+"), ("Kr*", "Kr"),
            ("te_pow", 1.6e-6, 300.0 * 8.617333262e-5, -0.55), "basic",
            "Kr2+ + e -> Kr* + Kr, 1.6e-6 (300 K/Te)^0.55 cm3/s +- 10 %, "
            "Shiu & Biondi, PRA 16, 1817 (1977) (abstract, Te 300-8400 K, "
            "Tgas 300 K); Mikus, J. Phys. D 11, L39 (1978) measured "
            "1.48e-6 (Te/300 K)^-0.53 over 300-21300 K, within 8 %. The "
            "products are mostly 5p states (Shiu & Biondi), which cascade "
            "to 5s as C1 assumes; the 2.93 eV difference heats the gas. No "
            "gas-temperature dependence has been measured"),
    ]
    rows += [replace(opt[rid], group="basic") for rid in
             ("ATT", "DIS_a3", "DIS_A1", "VIB1", "VIB2", "VIB3", "VIB4")]
    ion_ion = ("Flannery form of Kannari 1985 Sec. II.B (the modified "
               "Natanson calculation of Hoffman & Moreno, SAND80-1486, which "
               "Kannari used gas-temperature dependent), with Debye "
               "screening. Measured with F- only below 300 Torr: Jacob, "
               "Rokni, Mangano & Brochu, APL 32, 109 (1978) fitted KrF* "
               "pulse shapes in e-beam Kr/F2 to Kr+ + F- (+M) = 2e-7 at 60 "
               "Torr and 3e-7 cm3/s at 120 Torr (this form, screened at "
               "their 1e14 cm-3 ions, gives x1.2 and x1.6), and Rokni, Jacob "
               "& Mangano, PRA 16, 2216 (1977) found by changing the beam "
               "current that Kr2+ + F- forms KrF* 'with a branching ratio "
               "near unity', as efficiently as Kr+ + F- (relative "
               "intensities only). The absolute KrF* yield, the B/C split "
               "and the rate near 1 atm have no measurement (Brau 1984 p.98; "
               "docs/krf_data/FORMATION_LINKS_SEARCH.md). v0_1 restores the "
               "temperature dependence from Thomson theory: the three-body "
               "part scales as (Tg/300)^-3 at fixed density (b^3 prop. to "
               "T^-3, polarization mean free paths), the mobility part does "
               "not (reduced mobilities are constant); Mezyk et al., JPC 95, "
               "3152 (1991) use the same (300/T)^3 in their Bates scaling")
    rows += [
        Row("K26", ("Kr+", "F-"), ("KrF(B)",),
            ("flannery", 2.4, 3.417e-6, -3.0), "basic",
            "Kr+ + F- -> KrF(B), " + ion_ion),
        Row("K108", ("Kr2+", "F-"), ("KrF(B)", "Kr"),
            ("flannery", 2.65, 3.096e-6, -3.0), "basic",
            "Kr2+ + F- -> KrF(B) + Kr, " + ion_ion),
        Row("FF", ("F", "F", "Kr"), ("F2", "Kr"), ("const", FF_KR_300),
            "basic", "F + F + Kr -> F2 + Kr, 2.8e-34 cm6/s, the evaluation "
            "of Cohen & Bott (Aerospace TR-0076(6603)-2 (1976), DTIC "
            "ADA024332, pp.12-16: '10^13.7 exp(-35,100/RT) is recommended', "
            "10^14.0 cm6 mol-2 s-1 at 300 K; their 1982 update gives every "
            "third body but F, F2 and He the Ar value). No Kr measurement "
            "exists; measured at 295-298 K: (8.0 +- 0.6)e-35 for Ar "
            "(Ganguli & Kaufman 1974) and (6.0 +- 1.1)e-34 for He (Ultee "
            "1977), both via Baulch et al., JPCRD 10 Suppl. 1 (1981) "
            "p.1-22. Held at its 300 K value: the evaluation falls with Tg, "
            "so this overstates a channel that re-forms F2 on a millisecond "
            "scale (d[F2]/dt = k[F]^2[M], Baulch's convention)"),
    ]
    return rows


def basic_kr2_rows() -> dict:
    """The basic network's kr2 group: the Kr2* route with its measured
    constants (2026-09-28, step 2 of the owner's plan). They replace the
    REBUILD copies of the same rows, which v0 keeps (the drift gate guards
    them there)."""
    rows = [
        Row("AK8_Kr", ("Kr*", "Kr", "Kr"), ("Kr2*", "Kr"), ("const", 1.7e-32),
            "kr2", "Kr*(3P2) + 2Kr -> Kr2* + Kr, 1.7 +- 0.2e-32 cm6/s at 300 "
            "K, Ueda, Kanada, Sasaki & Takuma, Rev. Laser Eng. 12, 357 "
            "(1984) Table V (e-beam, Kr* absorption at 0.5-1 atm; read "
            "first-hand). Low-pressure afterglows gave 2.6-4.4e-32 (Tracy & "
            "Oskam 1976, Smith & Turner 1963, Bouciqué & Mortier 1970, "
            "abstracts). Tg form: Ea = -0.020 eV (Wieme & Lenaerts, JCP 74, "
            "483 (1981), abstract), (Tg/300)^-0.64 over 300-451 K. REBUILD's "
            "AK8_Kr is an estimate, 2.7e-32", tg_pow=-0.64),
        Row("KR2_RAD", ("Kr2*",), ("Kr", "Kr"), ("const", 3.7e6), "kr2",
            "Kr2* -> 2Kr + hv, 3.7e6 1/s = 270 ns, the slow (1u, 0u-) "
            "component that carries the kinetics: Morikawa et al., JCP 91, "
            "1469 (1989) 270 +- 5 ns; Salamero 1979 250 +- 20, Quigley & "
            "Hughes 1978 280 +- 30, Audouard 1988 265, Oka 1974 353 ns, "
            "Leichner & Ericson 1974 about 0.30 us at 1000 Torr (abstracts; "
            "J&H 1980's 3.3e6 agrees). REBUILD's KR2_RAD keeps Kannari's "
            "8.0e7, which none of the measurements supports", hv_eV=8.43,
            band="147 nm (Kr2* VUV)"),
        Row("NK7", ("Kr2*", "F2"), ("KrF(B)", "Kr", "F"), ("const", 1.05e-10),
            "kr2", "half of Kr2* + F2 = 2.1 +- 0.2e-10 cm3/s (Quigley & "
            "Hughes, APL 32, 649 (1978), abstract; products not measured). "
            "The even KrF*/Kr2F* split is Kannari's assignment (Table I "
            "p.4311 [52, 61], 3e-10 each)"),
        Row("K2F_BR", ("Kr2*", "F2"), ("Kr2F*", "F"), ("const", 1.05e-10),
            "kr2", "the other half of Quigley & Hughes' 2.1e-10 (as NK7); "
            "J&H #20. Chen & Payne, IEEE JQE 15, 149 (1979) find this "
            "direct channel minor (Fig. 25, f2'/f3 << 1: Kr2F* forms "
            "'through KrF* as an intermediate species', p.158) but do not "
            "give the KrF* share. The split moves formation by 1 % at most, "
            "so it stays even. Their Kr2* + F2 total is 5e-10 and SRI's "
            "3e-10: the measurements disagree beyond their errors"),
    ]
    return {r.rid: r for r in rows}


def build(groups=(), net="basic") -> list:
    """The network for ``net`` ('basic', the default, or 'v0', the notes'
    list as written) plus optional groups ('all' = every group that fits the
    network). The basic network already carries measured quench and B/C
    rows, so the 'quench' and 'bc' groups are refused there."""
    if net not in ("basic", "v0"):
        raise NetworkError(f"unknown network {net!r}; have basic, v0")
    allowed = BASIC_GROUPS if net == "basic" else GROUPS
    groups = allowed if "all" in groups else tuple(groups)
    bad = set(groups) - set(allowed)
    if bad:
        why = (" (the basic network already has measured quench, B/C, "
               "ion and F2-electron rows)"
               if net == "basic" and bad <= {"quench", "bc", "ions", "f2e"}
               else "")
        raise NetworkError(f"groups {sorted(bad)} not available for "
                           f"{net!r}{why}; have {allowed}")
    base = basic_rows() if net == "basic" else v0_rows()
    opt = optional_rows()
    if net == "basic":
        own = basic_kr2_rows()
        opt = [own.get(r.rid, r) for r in opt]
    rows = base + [r for r in opt if r.group in groups]
    if "bc" in groups:
        rows += c_state_copies(rows)
    seen = set()
    for r in rows:
        if r.rid in seen:
            raise NetworkError(f"duplicate row id {r.rid}")
        seen.add(r.rid)
        check_row(r)
    return rows
