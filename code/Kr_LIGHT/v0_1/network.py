"""The pure-krypton reaction set AS DATA (LIGHT_Kr v0_1, 2026-10-13).

v0_1 is LIGHT v0_5's basic network with every row that touches fluorine
taken out: the 18 rows among its 59 whose reactants and products are
krypton and electrons alone, each with v0_5's rate, citation and energy
routing unchanged, in v0_5's order. The rows taken out are the fluorine
chemistry (F2's electron collisions, the harpoon reactions of Kr*, Kr**
and Kr2*, the KrF(B,C) block, the trimer Kr2F*, ion-ion recombination with
F-, F + F + Kr) and their bands. With no F2 in the fill those rows carry
nothing, so v0_1 runs v0_5's own pure-Kr chemistry (tests/test_v0_1.py).

  G1      e + Kr -> e + Kr*               (5s manifold, Biagi 'Kr(s)' 9.915 eV)
  G2      Kr* + Kr* -> Kr+ + Kr + e       (Penning, band kr_pen)
  G6      e + Kr elastic                  (energy loss to the gas)
  C1a-c   e + Kr -> e + Kr**              (the 5p, d+p and higher groups,
                                           to the 5p rung)
  C2      e + Kr -> 2e + Kr+              (direct ionization)
  K15     Kr+ + 2Kr -> Kr2+ + Kr
  R11p    Kr2+ + e -> Kr** + Kr           (dissociative recombination)
  AK8_Kr  Kr* + 2Kr -> Kr2* + Kr          (band kr2_form)
  KR2_RAD Kr2* -> 2Kr + hv, 147 nm        (the second continuum)
  PEN_KR2 Kr2* + Kr2* -> Kr2+ + 2Kr + e   (band kr2_pen, 0: no measurement)
  KPS     Kr** + Kr -> Kr* + Kr           (band kps_kr)
  PRAD    Kr** -> Kr* + hv, near IR       (the 5p -> 5s lines)
  K76     e + Kr* -> e + Kr**             (Jung et al.'s measured cross
  K76r    e + Kr** -> e + Kr*              sections, band k_5s5p; reverse
                                           by detailed balance)
  KSTEP_S e + Kr* -> Kr+ + 2e             (bands kstep_s, kstep_p, 0: no
  KSTEP_P e + Kr** -> Kr+ + 2e             measurement)

The optional group 'electron' (Params.groups / --add) adds KSTEP, REBUILD's
stepwise-ionization template, and SUPER, the superelastic return of G1.

Units: cm^3/s (two-body), cm^6/s (three-body), 1/s (radiative).

Rate specs (declarative, evaluated by model.py):
  ("const", k)                      constant
  ("xsec", gas, kind, product)      Maxwellian <sigma v>(Te) of one LXCat block
  ("elastic", gas)                  G6: energy loss only, no density change
  ("stepwise", E_bind)              1.56e-7 Te^0.71 exp(-E/Te)  (REBUILD)
  ("te_pow", k0, Te0, n)            k0 (Te/Te0)^n, Te held above 0.02 eV;
                                    a captured electron takes (3/2 + n) Te
  ("superelastic", gas, product, g_ratio)
                                    de-excitation from the Klein-Rosseland
                                    cross section of an LXCat excitation
                                    (= k_fwd g_ratio exp(E_thr/Te))
  ("param", name[, frac])           frac x the value of band `name`
                                    (BANDS): Params.rates if set there,
                                    else the band's default
  ("te_fit", A, B, C)               A Te^B exp(-C/Te), Kannari's printed
                                    form for the electron rows of the 5p
                                    rung. A may be a band name, read as
                                    ("param", name) is. Te is held above
                                    0.026 eV
  ("te_fit", A, B, C, "reverse", g_ratio, upper, lower)
                                    the reverse of a te_fit row by detailed
                                    balance over the ledger's gap dE =
                                    E_form(upper) - E_form(lower): A g_ratio
                                    Te^B exp(-(C - dE)/Te)
  ("pop_xsec", "jung2006"[, "reverse", g_ratio, upper, lower])
                                    the 5s -> 5p excitation on Jung et al.'s
                                    measured cross sections (jung2006.py):
                                    its Maxwellian table, and the reverse its
                                    Klein-Rosseland table over the ledger's
                                    gap; under the kinetic EEDF a population
                                    process of the Boltzmann equation

A row's `frac` multiplies its rate by band values: (("name", "x"), ...)
by the band's value, (("name", "1-x"), ...) by one minus it.
"""

from __future__ import annotations

from dataclasses import dataclass

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
    frac: tuple = ()               # ((band, "x" | "1-x"), ...) rate factors

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
      cascade to the 5p rung);
    * exothermic: to the electron gas when the row ejects a new electron from
      heavy reactants (Penning) or says energy_to='electrons'
      (superelastic), otherwise to the gas;
    * endothermic without a threshold: from the electron gas if an electron
      collides, otherwise from the gas.
    The captured electron's energy (dissociative recombination) moves
    electron gas -> gas at run time ('capture')."""
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
# constants (LIGHT v0_1 to v0_5, unchanged)
# ---------------------------------------------------------------------------
KR_CONV_300 = 2.4e-31       # Kr+ + 2Kr -> Kr2+ + Kr, cm6/s at 300 K
BASIC_GROUPS = ("electron",)           # the groups that add to the network
GROUPS = BASIC_GROUPS

# Low-field reduced mobilities K0 [cm2/V/s] at the Loschmidt density, the
# drift-tube compilation of Ellis, Pai, McDaniel, Mason & Viehland, At.
# Data Nucl. Data Tables 17, 177 (1976) (Kr+, Kr2+ in Kr). The same table
# Fluid1D's 1-D transport reads
# (nathaniel/1DModelKrF/data/ion_transport/ion_mobility.toml, +-3 to 5 %).
# mu = K0 N0 / N. They set the ions' drift to the walls (model.py,
# Params.gap_cm).
ION_K0 = {"Kr+": 0.84, "Kr2+": 1.22}
AMAGAT_CM3 = 2.6867811e19         # N0, Loschmidt

# The 5p rung. Statistical weights of the lumped manifolds (sum of 2J+1):
# 5s (4 levels) 12, 5p (10 levels) 36.
G_5S, G_5P = 12, 36
# Kannari, Obara & Fujioka, JAP 57, 4309 (1985) Table I, the electron rows
# of the 5p rung, A Te^B exp(-C/Te) cm3/s (read from the rendered page,
# IB/quant_IB/qib3_rows.py checks the digits against the PDF text):
#   [76] Kr* + e -> Kr** + e   9.00e-7 Te^0.50 exp(-1.51/Te)  (Jacob &
#        Mangano, APL 28, 724 (1976): a fit to their e-beam KrF data)
#   [79] Kr* + e -> Kr+ + 2e   7.00e-8 Te^0.68 exp(-4.04/Te)  (Vriens,
#        Physica 31, 385 (1965), a binary-encounter calculation)
#   [79] Kr** + e -> Kr+ + 2e  2.00e-7 Te^0.62 exp(-2.53/Te)  (Vriens)
KANNARI_K76 = (9.00e-7, 0.50, 1.51)
KANNARI_KSTEP_S = (7.00e-8, 0.68, 4.04)
KANNARI_KSTEP_P = (2.00e-7, 0.62, 2.53)
# the preset Params.ladder = "kannari": every ladder band at Kannari's
# printed value at once (run_v0_1.py --ladder kannari), a sensitivity
# configuration and never a default (the unmeasured rows default to 0).
# LIGHT's preset also set Kr** + F2 (k5p_f2), which has no row here
LADDER_KANNARI = {"kstep_s": KANNARI_KSTEP_S[0],
                  "kstep_p": KANNARI_KSTEP_P[0], "kps_kr": 1.0e-10}
# the ladder bands with no measurement, 0 by default (BANDS): the preset
# turns these on; kps_kr carries a measurement and the preset moves it to
# Kannari's mark
LADDER_UNMEASURED = ("kstep_s", "kstep_p")
# the 5p -> 5s radiation: one lifetime for the lumped manifold, the
# statistical-weight mean of the ten measured 5p radiative rates, 29 ns.
# Measured lifetimes [ns] (Dzierzega, Volz, Nave & Griesmann, PRA 62,
# 022505 (2000), Table 1, their collected experiments; the beam-gas-laser
# values of Schmitt et al. where measured, else Chang, Horiguchi & Setser,
# JCP 73, 778 (1980), Whitehead et al., JCP 102, 1965 (1995), Cannon et
# al., PRA 47, 147 (1993), Ernst & Schulz-Gulde, Physica 93C, 136 (1978)):
# 2p1 24.6, 2p2 26.9-32.0, 2p3 28.1, 2p4 29.4, 2p5 23.5, 2p6 27.3, 2p7
# 29.6, 2p8 32.0, 2p9 27.7, 2p10 40.9. Weighted by 2J+1 (1,5,3,3,1,5,3,5,
# 7,3) the mean decay rate is 1/29.0 ns. The lines (760-830 nm) end on the
# 5s levels and are not trapped
KR_5P_TAU_NS = {"2p1": 24.6, "2p2": 28.6, "2p3": 28.1, "2p4": 29.4,
                "2p5": 23.1, "2p6": 27.3, "2p7": 29.6, "2p8": 32.0,
                "2p9": 27.7, "2p10": 38.1}
KR_5P_G = {"2p1": 1, "2p2": 5, "2p3": 3, "2p4": 3, "2p5": 1, "2p6": 5,
           "2p7": 3, "2p8": 5, "2p9": 7, "2p10": 3}
KR_5P_TAU_S = 1e-9 * sum(KR_5P_G.values()) / sum(
    g / KR_5P_TAU_NS[k] for k, g in KR_5P_G.items())


# ---------------------------------------------------------------------------
# The basic network: krypton and electrons
# ---------------------------------------------------------------------------
def excitation_rows() -> list:
    """The notes' krypton rows and their closures (KrF_model_notes_with_
    Gary.rtf, as LIGHT's v0 wrote them, and LIGHT's basic network since):
    G1 makes the 5s manifold, G2 is Kr* + Kr* Penning ionization on its
    band, G6 the elastic loss, the higher excitations C1a-c make the 5p
    rung Kr** (LIGHT v0_4) and C2 ionizes. The notes' G3, G4, G5 and C3
    are fluorine's and are not here."""
    return [
        Row("G1", ("e", "Kr"), ("e", "Kr*"),
            ("xsec", "Kr", "excitation", "Kr(s)(9.915eV)"), "gary",
            "notes (1) e + Kr -> Kr* 5s[3/2]. Biagi-v7.1 'Kr(s)' = "
            "the whole 5s manifold (5s[3/2]1,2 and 5s'[1/2]0,1) lumped at "
            "9.915 eV; k = <sigma v> over a Maxwellian (notes: Q = k n_e[Kr])",
            e_cost_eV=9.915),
        Row("G2", ("Kr*", "Kr*"), ("Kr+", "Kr", "e"), ("param", "kr_pen"),
            "basic", "notes (2) Kr* + Kr* 'associative ionization', 5.0e-10 "
            "cm3/s, Kannari 1985 Table I Penning block from his ref. 86, "
            "Leichner & Ericson, PRA 9, 251 (1974), with Kr+ + Kr + e as the "
            "products (v0's value). Band kr_pen, since Jancaitis's LLNL model "
            "(UCRL-53465, 1983, Table II p.86) carries 1.5e-10 from Werner's "
            "MIT thesis (1975). Neither original is among the sources held"),
        Row("G6", ("e", "Kr"), ("e", "Kr"), ("elastic", "Kr"),
            "gary", "notes (6) e + Kr elastic. Biagi-v7.1 momentum-transfer "
            "set; power to the gas (2m/M)(1 - Tg/Te) n_e N <sigma_m v eps>"),
        Row("C1a", ("e", "Kr"), ("e", "Kr**"),
            ("xsec", "Kr", "excitation", "Kr(p)(11.30eV)"), "closure",
            "Biagi-v7.1 5p group (11.30 eV). Not in the notes: without it "
            "the electrons lose ~half their inelastic power nowhere. v0_4: "
            "the product is the 5p rung Kr** (11.30 eV, the group's own "
            "threshold), which relaxes to Kr* by KPS, radiates (PRAD) or is "
            "moved on by the electrons (ladder_rows) (v0_4: product Kr**)",
            e_cost_eV=11.30),
        Row("C1b", ("e", "Kr"), ("e", "Kr**"),
            ("xsec", "Kr", "excitation", "Kr(d+p)(11.998eV)"), "closure",
            "Biagi-v7.1 d+p group (11.998 eV); cascade to the 5p rung (Kr**) "
            "taken prompt, the excess heats the gas (v0_4: product Kr**, "
            "0.698 eV to heat)", e_cost_eV=11.998),
        Row("C1c", ("e", "Kr"), ("e", "Kr**"),
            ("xsec", "Kr", "excitation", "Kr*(12.75eV)"), "closure",
            "Biagi-v7.1 higher group (12.75 eV); cascade to the 5p rung "
            "(Kr**) taken prompt, the excess heats the gas (v0_4: product "
            "Kr**, 1.450 eV to heat)", e_cost_eV=12.75),
        Row("C2", ("e", "Kr"), ("e", "e", "Kr+"),
            ("xsec", "Kr", "ionization", "Kr^+"), "closure",
            "Biagi-v7.1 ionization (13.996 eV): the notes' 'Kr+' "
            "output needs electron-impact ionization", e_cost_eV=13.996),
    ]


def charge_rows() -> list:
    """LIGHT v0_1's charge chemistry, its krypton part: no electrons are
    injected, so n_e follows ionization and recombination (quasineutral,
    n_e = Kr+ + Kr2+). Kr+ converts to Kr2+ in three-body collisions and
    Kr2+ recombines dissociatively onto the 5p rung (v0_4)."""
    return [
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
        Row("R11p", ("e", "Kr2+"), ("Kr**", "Kr"),
            ("te_pow", 1.6e-6, 300.0 * 8.617333262e-5, -0.55), "basic",
            "Kr2+ + e -> Kr** + Kr, 1.6e-6 (300 K/Te)^0.55 cm3/s +- 10 %, "
            "Shiu & Biondi, PRA 16, 1817 (1977) (abstract, Te 300-8400 K, "
            "Tgas 300 K); Mikus, J. Phys. D 11, L39 (1978) measured "
            "1.48e-6 (Te/300 K)^-0.53 over 300-21300 K, within 8 %. The "
            "products are mostly 5p states (Shiu & Biondi observed the 5p "
            "-> 5s lines in the afterglow), so v0_4 makes the 5p rung Kr** "
            "(11.30 eV) and books 1.546 eV of heat per event, where v0_3 "
            "made Kr* and booked 2.93 eV. No gas-temperature dependence "
            "has been measured"),
    ]


def kr2_rows() -> list:
    """LIGHT v0_2's Kr2* route, its krypton part: Kr* + 2Kr makes the
    excimer, which radiates the 147-nm second continuum. The measured rows
    carry their measurements and read their bands (BANDS) where the
    published values spread. Kr2* + Kr2* has no measurement and is 0
    unless a band run sets it."""
    return [
        Row("AK8_Kr", ("Kr*", "Kr", "Kr"), ("Kr2*", "Kr"),
            ("param", "kr2_form"), "basic",
            "Kr*(3P2) + 2Kr -> Kr2* + Kr, 1.7 +- 0.2e-32 cm6/s at 300 K, "
            "Ueda, Kanada, Sasaki & Takuma, Rev. Laser Eng. 12, 357 (1984) "
            "Table V (e-beam, Kr* absorption at 0.5-1 atm, read "
            "first-hand). Band kr2_form carries the low-pressure "
            "afterglows, 1.79 to 4.4e-32 (Salamero 1979, Tracy & Oskam 1976, "
            "Smith & Turner 1963, Bouciqué & Mortier 1970, abstracts), and "
            "Hughes 1977, 5.5e-32 (via Quigley & Hughes 1978 Table I). "
            "Jancaitis's model uses Tracy & Oskam's 2.6e-32. Tg form Ea = "
            "-0.020 eV (Wieme & Lenaerts, JCP 74, 483 (1981), abstract), "
            "(Tg/300)^-0.64 over 300-451 K. REBUILD's AK8_Kr is an "
            "estimate, 2.7e-32", tg_pow=-0.64),
        Row("KR2_RAD", ("Kr2*",), ("Kr", "Kr"), ("const", 3.7e6), "basic",
            "Kr2* -> 2Kr + hv, 3.7e6 1/s = 270 ns, the slow (1u, 0u-) "
            "component that carries the kinetics: Morikawa et al., JCP 91, "
            "1469 (1989) 270 +- 5 ns; Salamero 1979 250 +- 20, Quigley & "
            "Hughes 1978 280 +- 30, Audouard 1988 265, Oka 1974 353 ns, "
            "Leichner & Ericson 1974 about 0.30 us at 1000 Torr (abstracts; "
            "J&H 1980's 3.3e6 agrees). REBUILD's KR2_RAD keeps Kannari's "
            "8.0e7, which none of the measurements supports", hv_eV=8.43,
            band="147 nm (Kr2* VUV)"),
        Row("PEN_KR2", ("Kr2*", "Kr2*"), ("Kr2+", "Kr", "Kr", "e"),
            ("param", "kr2_pen"), "basic",
            "Kr2* + Kr2* -> Kr2+ + 2Kr + e, no measurement, so 0 unless a "
            "band run sets it (kr2_pen, 0 to 5.0e-10 cm3/s). Jancaitis's "
            "model carries 1.0e-10 from Werner et al., IEEE J. Quant. Elect. "
            "13, 769 (UCRL-53465 Table II p.86 and ref. 1 p.109). Kannari "
            "1985 Table I p.4312 [87] (REBUILD PEN_KR2). Ref. 87 is Lorents "
            "& Olson, SRI semiannual technical report N00014-72-C-0557 "
            "(1973, unpublished), and Kannari gives Ar2* + Ar2* the same "
            "5.0e-10 from it (read first-hand 2026-09-29). No measurement of "
            "this rate is among the sources held"),
    ]


def kannari_k76_rows() -> list:
    """LIGHT v0_4's K76 pair, Kannari's Table I [76] fit A Te^0.5
    exp(-C/Te) and its detailed-balance reverse, for Params.k76 =
    'kannari' (a comparison configuration; the band k_5s5p scales it as it
    scales Jung's)."""
    return [
        Row("K76", ("e", "Kr*"), ("e", "Kr**"),
            ("te_fit",) + KANNARI_K76, "basic",
            "e + Kr(5s) -> e + Kr(5p), Kannari 1985 Table I [76], 9.0e-7 "
            "Te^0.5 exp(-1.51/Te) cm3/s, Jacob & Mangano's fit to e-beam "
            "KrF data (a fit, not a direct measurement; v0_4's default). "
            "The 1.385 eV comes from the electron gas",
            frac=(("k_5s5p", "x"),)),
        Row("K76r", ("e", "Kr**"), ("e", "Kr*"),
            ("te_fit",) + KANNARI_K76
            + ("reverse", G_5S / G_5P, "Kr**", "Kr*"), "basic",
            "e + Kr(5p) -> e + Kr(5s), the reverse of Kannari's K76 by "
            "detailed balance, g(5s)/g(5p) = 12/36 and exp(+1.385/Te) over "
            "the ledger's gap; its 1.385 eV returns to the electron gas",
            energy_to="electrons", frac=(("k_5s5p", "x"),)),
    ]


def ladder_rows() -> list:
    """LIGHT v0_4: the 5p rung Kr** and its rows. Kr** is made by the
    higher excitations (C1a-c) and by Kr2+ + e (R11p). It relaxes to Kr*
    in collisions with Kr (KPS), radiates to the 5s levels (PRAD, near-IR
    lines that are not trapped), and the electrons pump it from Kr* (K76),
    return it (K76r, detailed balance) and ionize it (KSTEP_P); KSTEP_S
    ionizes Kr*. Measured rows carry their measurements (PRAD, KPS, and
    K76 on Jung et al.'s cross sections since v0_5); the two stepwise
    ionizations have no measurement and are bands at 0 (BANDS).
    Params.ladder = 'kannari' sets every ladder band to Kannari's printed
    value at once."""
    return [
        Row("KPS", ("Kr**", "Kr"), ("Kr*", "Kr"), ("param", "kps_kr"),
            "basic", "Kr(5p) + Kr -> Kr(5s) + Kr, the p -> s relaxation. "
            "Band kps_kr, default (1.7 +- 0.2)e-10 cm3/s: the one Kr "
            "measurement in hand, the total two-body loss of Kr 5p'[3/2]2 "
            "(2p2) in Kr at 300 K by ps-TALIF (Gazeli et al., Phys. "
            "Plasmas 28, 043301 (2021), Stern-Volmer over 0.1-10 mbar, "
            "'in the range of values reported in the literature'; the "
            "product channel is not resolved). Whitehead et al., JCP 102, "
            "1965 (1995) measured the 5p + Kr constants level by level "
            "(their table is not among the sources held) and Chang, "
            "Horiguchi & Setser, JCP 73, 778 (1980) the Ar partner. The "
            "band runs down to Jancaitis's 3.0e-11, the Ar analogue "
            "(UCRL-53465 Table II, from Chang & Setser 1978), past "
            "Kannari's 1.0e-10 'estimated' for Ar (REBUILD KPS_Kr). The "
            "1.385 eV gap heats the gas"),
        Row("PRAD", ("Kr**",), ("Kr*",), ("const", 1.0 / KR_5P_TAU_S),
            "basic", "Kr(5p) -> Kr(5s) + hv, the near-infrared 5p -> 5s "
            "lines (760-830 nm), which end on the 5s levels and are not "
            "trapped. One lifetime for the manifold, 29.0 ns, the "
            "statistical-weight mean of the ten measured 5p radiative "
            "rates (23 to 41 ns; Dzierzega et al. 2000 Table 1, collecting "
            "Schmitt et al.'s beam-gas-laser values, Chang, Horiguchi & "
            "Setser 1980, Whitehead et al. 1995, Cannon et al. 1993 and "
            "Ernst & Schulz-Gulde 1978; KR_5P_TAU_NS)",
            hv_eV=11.30 - 9.915, band="5p-5s (near IR)"),
        Row("K76", ("e", "Kr*"), ("e", "Kr**"),
            ("pop_xsec", "jung2006"), "basic",
            "e + Kr(5s) -> e + Kr(5p) on the measured cross sections out "
            "of the two metastables into the ten 5p levels (Jung, Stone, "
            "Boffard, Anderson & Lin, PRA 73, 022722 (2006) Tables II-IV; "
            "PRL 94, 163202 (2005)), the statistical-weight mean over the "
            "5s manifold with the J = 1 levels given their same-core "
            "metastable's set (jung2006.py). v0_5: replaces Kannari's "
            "Table I [76] fit (Jacob & Mangano), kept as Params.k76 = "
            "'kannari'. Band k_5s5p is a scale on the cross section, 1 by "
            "default, 0.5 (the metastables alone) to 1.4 (the measurement's "
            "+40 %). The 1.385 eV comes from the electron gas",
            frac=(("k_5s5p", "x"),)),
        Row("K76r", ("e", "Kr**"), ("e", "Kr*"),
            ("pop_xsec", "jung2006", "reverse", G_5S / G_5P, "Kr**", "Kr*"),
            "basic",
            "e + Kr(5p) -> e + Kr(5s), the reverse of K76 by the "
            "Klein-Rosseland relation on the same cross section, g(5s)/"
            "g(5p) = 12/36 over the ledger's gap; its 1.385 eV returns to "
            "the electron gas", energy_to="electrons",
            frac=(("k_5s5p", "x"),)),
        Row("KSTEP_S", ("e", "Kr*"), ("Kr+", "e", "e"),
            ("te_fit", "kstep_s") + KANNARI_KSTEP_S[1:], "basic",
            "e + Kr(5s) -> Kr+ + 2e, stepwise ionization. Band kstep_s: no "
            "measurement; Kannari 1985 Table I [79] prints Vriens's "
            "binary-encounter calculation 7.0e-8 Te^0.68 exp(-4.04/Te), and "
            "REBUILD's template (LIGHT's optional KSTEP) runs 2.2x above it "
            "at 1.3 eV, the band's top. The 4.08 eV comes from the electron "
            "gas (the ledger's IP - E(5s))"),
        Row("KSTEP_P", ("e", "Kr**"), ("Kr+", "e", "e"),
            ("te_fit", "kstep_p") + KANNARI_KSTEP_P[1:], "basic",
            "e + Kr(5p) -> Kr+ + 2e. Band kstep_p: no measurement; Kannari "
            "1985 Table I [79], Vriens's calculation 2.0e-7 Te^0.62 "
            "exp(-2.53/Te). The 2.70 eV comes from the electron gas"),
    ]


def basic_rows() -> list:
    """The pure-Kr network, in LIGHT v0_5's order: the excitations, the
    charge chemistry, the Kr2* route and the 5p rung."""
    return excitation_rows() + charge_rows() + kr2_rows() + ladder_rows()


# ---------------------------------------------------------------------------
# The optional group 'electron': electron collisions with Kr* (OFF by default)
# ---------------------------------------------------------------------------
def optional_rows() -> list:
    return [
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
    ]


# ---------------------------------------------------------------------------
# Bands (LIGHT v0_2, 2026-09-29): the rates known only as a published
# spread, or not at all
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Band:
    """A rate the basic network carries as a band. Rows read it through a
    ("param", name[, frac]) rate, a ("te_fit", name, ...) rate or a frac
    factor, and run_v0_1.py --bands runs each end. A measured rate keeps
    its measurement as the default and spans the other published values.
    A rate with no measurement is 0 by default, and its band runs from 0 to
    the largest published value, which is no bound. Device data and
    analogs never set a default or narrow a band."""
    name: str
    what: str
    default: float
    lo: float
    hi: float
    unit: str
    cite: str
    marks: tuple = ()          # ((value, who), ...) published values inside


BANDS = {b.name: b for b in (
    Band("kr2_pen", "Kr2* + Kr2* -> Kr2+ + 2Kr + e", 0.0, 0.0, 5.0e-10,
         "cm3/s",
         "No measurement. Kannari's 5.0e-10 is Lorents & Olson's "
         "unpublished 1973 SRI report, which gives Ar2* the same value, and "
         "Jancaitis's model carries 1.0e-10 from Werner et al. 1977",
         marks=((1.0e-10, "Werner et al. via Jancaitis"),
                (5.0e-10, "Kannari 1985"))),
    Band("kr2_form", "Kr* + 2Kr -> Kr2* + Kr at 300 K", 1.7e-32, 1.7e-32,
         5.5e-32, "cm6/s",
         "Ueda 1984 measured 1.7 +- 0.2e-32 at 0.5-1 atm, the default. "
         "Low-pressure afterglows give 1.79 to 4.4e-32 and Hughes 1977 "
         "5.5e-32",
         marks=((2.59e-32, "Tracy & Oskam 1976"),
                (4.4e-32, "Bouciqué & Mortier 1970"))),
    Band("kr_pen", "Kr* + Kr* -> Kr+ + Kr + e", 5.0e-10, 1.5e-10, 5.0e-10,
         "cm3/s",
         "Kannari 1985 gives 5.0e-10 from Leichner & Ericson 1974, v0's and "
         "the default. Jancaitis's model carries 1.5e-10 from Werner's "
         "thesis (1975)"),
    # --- LIGHT v0_4 (2026-10-08) --------------------------------------------
    Band("kps_kr", "Kr(5p) + Kr -> Kr(5s) + Kr", 1.7e-10, 3.0e-11, 1.7e-10,
         "cm3/s",
         "Gazeli et al. 2021 measured (1.7 +- 0.2)e-10 for the total loss "
         "of Kr 5p'[3/2]2 in Kr at 300 K, the default (one level, product "
         "unresolved). Jancaitis's 3.0e-11 is the Ar analogue (Chang & "
         "Setser 1978), the low end; Kannari's 1.0e-10 'estimated' for Ar "
         "(REBUILD KPS_Kr) lies inside",
         marks=((3.0e-11, "Jancaitis 1984, Ar analogue"),
                (1.0e-10, "Kannari 1985 Ar estimate, REBUILD KPS_Kr"))),
    Band("k_5s5p", "e + Kr(5s) -> e + Kr(5p), a scale on the measured "
         "lumped cross section (jung2006.py)", 1.0, 0.5, 1.4, "factor",
         "Jung et al., PRA 73, 022722 (2006): the cross sections out of "
         "the two metastables, +-35 to 40 %; nothing measured out of the "
         "two J = 1 resonance levels, given their same-core metastable's "
         "set. 1.4 is the measurement's upper error; 0.5 is the "
         "metastables alone (the J = 1 levels, half the manifold's weight, "
         "contributing nothing). Kannari's fit is Params.k76 = 'kannari'",
         marks=((1.0, "Jung 2006, the manifold on the proxies"),
                (0.5, "the measured metastables alone"))),
    Band("kstep_s", "e + Kr(5s) -> Kr+ + 2e, the pre-exponential of "
         "Kannari's Te^0.68 exp(-4.04/Te)", 0.0, 0.0, 1.5e-7, "cm3/s",
         "No measurement. Kannari 1985 Table I [79] prints Vriens's "
         "calculation 7.0e-8; REBUILD's template 1.56e-7 Te^0.71 "
         "exp(-4.09/Te) (LIGHT's optional KSTEP) equals 1.5e-7 on this "
         "form at 1.3 eV, the top",
         marks=((7.0e-8, "Kannari 1985 [79], Vriens"),)),
    Band("kstep_p", "e + Kr(5p) -> Kr+ + 2e, the pre-exponential of "
         "Kannari's Te^0.62 exp(-2.53/Te)", 0.0, 0.0, 2.0e-7, "cm3/s",
         "No measurement. Kannari 1985 Table I [79], Vriens's calculation "
         "2.0e-7, the top",
         marks=((2.0e-7, "Kannari 1985 [79], Vriens"),)),
    Band("ib_scale", "the electron-ion (Turnbull) share of the EM wave's "
         "absorbed power, a factor", 1.0, 0.9, 1.1, "factor",
         "Turnbull et al. 2023-2024 state their inverse-bremsstrahlung "
         "model reproduces the measured absorption to about 10 %, the "
         "band (ib.py TURNBULL_MODEL_ERROR). It scales the electron-ion "
         "part of the power LIGHT books, under either closure. "
         "Params.rates sets it; no row reads it"),
)}
# bands no row reads: Params.rates may still set them (model.py)
FREE_BANDS = ("ib_scale",)


def build(groups=(), k76="jung") -> list:
    """The pure-Kr network plus optional groups ('all' = every group; the
    one group is 'electron'). k76: 'jung' (the measured 5s -> 5p cross
    sections, the default) or 'kannari' (LIGHT v0_4's fit)."""
    if k76 not in ("jung", "kannari"):
        raise NetworkError(f"k76 {k76!r}: use 'jung' or 'kannari'")
    groups = GROUPS if "all" in groups else tuple(groups)
    bad = set(groups) - set(GROUPS)
    if bad:
        raise NetworkError(f"groups {sorted(bad)} not available; have "
                           f"{GROUPS} (the fluorine groups of LIGHT are "
                           "not in the pure-Kr model)")
    base = basic_rows()
    if k76 == "kannari":
        swap = {r.rid: r for r in kannari_k76_rows()}
        base = [swap.get(r.rid, r) for r in base]
    rows = base + [r for r in optional_rows() if r.group in groups]
    seen = set()
    for r in rows:
        if r.rid in seen:
            raise NetworkError(f"duplicate row id {r.rid}")
        seen.add(r.rid)
        check_row(r)
        if r.rate[0] == "param" and r.rate[1] not in BANDS:
            raise NetworkError(f"{r.rid}: no band {r.rate[1]!r}")
        if r.rate[0] == "te_fit" and isinstance(r.rate[1], str) \
                and r.rate[1] not in BANDS:
            raise NetworkError(f"{r.rid}: no band {r.rate[1]!r}")
        if r.rate[0] == "pop_xsec" and r.rate[1] != "jung2006":
            raise NetworkError(f"{r.rid}: pop_xsec source {r.rate[1]!r}")
        for name, how in r.frac:
            if name not in BANDS or how not in ("x", "1-x"):
                raise NetworkError(f"{r.rid}: frac ({name!r}, {how!r})")
    return rows


def bands_read(rows) -> set:
    """The bands the rows read, through 'param' and 'te_fit' rates and
    through `frac` factors."""
    out = set()
    for r in rows:
        if r.rate[0] == "param" or (r.rate[0] == "te_fit"
                                    and isinstance(r.rate[1], str)):
            out.add(r.rate[1])
        out.update(name for name, _ in r.frac)
    return out
