"""The pure-krypton reaction set AS DATA (LIGHT_Kr v0_4, 2026-10-15).

v0_1 was LIGHT v0_5's basic network with every row that touches fluorine
taken out: v0_5's 18 krypton rows, in its order. v0_2 books the electron's
collisions with ground-state krypton level by level, on Biagi's Magboltz
8.97 cross sections (xsec.py): each of the set's 51 excitation levels is
its own row at its own threshold, where v0_1's five rows (G1, C1a-c, C2)
read Biagi-v7.1's four groups at the groups' lowest thresholds. v0_2 adds
the Hornbeck-Molnar branch of the levels from which an excited atom can
ionize associatively on a ground-state one (HM_*, band hm_kr, 0 by
default: no measurement for krypton). v0_4 resolves the 5s manifold: each
of the four levels is its own species, Kr*(1s5), Kr*(1s4), Kr*(1s3) and
Kr*(1s2), with its own rows, and the associations have their reverses by
detailed balance (thermo.py).

  G1_<1s>   e + Kr -> e + Kr*(<1s>)       (the four 5s levels, 9.92-10.64 eV)
  G2_<a>_<b> Kr*(a) + Kr*(b) -> Kr+ + Kr + e  (Penning, band kr_pen; a pair
                                           of two levels counts both orders)
  G6        e + Kr elastic                (energy loss to the gas)
  C1a_<2p>  e + Kr -> e + Kr**            (the ten 5p levels, 11.30-12.26 eV)
  C1b_<lvl> e + Kr -> e + Kr**            (the 4d and 6s levels, 12.00-12.39
                                           eV, cascading to the 5p rung)
  C1c_<lvl> e + Kr -> e + Kr**            (every level from 12.76 eV up,
                                           cascading to the 5p rung)
  HM_<lvl>  e + Kr (+ Kr) -> 2e + Kr2+    (the levels at or above E(Kr2+) =
                                           12.85 eV ionizing associatively,
                                           Hornbeck-Molnar; band hm_kr)
  C2        e + Kr -> 2e + Kr+            (direct ionization)
  K15, K15r Kr+ + 2Kr <-> Kr2+ + Kr       (and its redissociation)
  R11p      Kr2+ + e -> Kr** + Kr         (dissociative recombination; band
                                           dr_tg, its gas-temperature factor)
  AK8_<1s>  Kr*(<1s>) + 2Kr -> Kr2* + Kr  (bands kr2_form, kr2_form_1s4,
                                           kr2_form_1s3; 1s2 forms through
                                           1s4)
  AK8f_1s4  Kr*(1s4) + 2Kr -> Kr2*(0u+) + Kr  (band kr2_fast, 0)
  AK8r_<1s>, AK8fr_1s4                    (the redissociations, detailed
                                           balance; band kr2_d0)
  KR2_RAD   Kr2* -> 2Kr + hv, 147 nm      (the slow 1u, 0u- states, 270 ns)
  KR2F_RAD  Kr2*(0u+) -> 2Kr + hv, 147 nm (the fast state, 3.4 ns)
  QKR2      Kr2* + Kr -> 3Kr              (band kr2_quench, 0)
  MIX_*     the 5s levels' transfers by Kr (Turner's 1s4 -> 1s5, Audouard's
            1s2 -> 1s4) and their reverses by detailed balance
  RES_1s4, RES_1s2  the escaping resonance lines (band res_escape, 0)
  PEN_KR2   Kr2* + Kr2* -> Kr2+ + 2Kr + e (band kr2_pen, 0: no measurement)
  KPS_<1s>, KPS_low  Kr** + Kr -> Kr*(<1s>) + Kr  (band kps_kr; band kps_5s
                                           splits the products)
  PRAD_<1s> Kr** -> Kr*(<1s>) + hv, near IR  (the 5p -> 5s lines, NIST's
                                           branching)
  K76_<1s>  e + Kr*(<1s>) -> e + Kr**     (Jung et al.'s measured cross
  K76r_<1s> e + Kr** -> e + Kr*(<1s>)      sections, bands k_5s5p and
                                           k_5s5p_j1; reverse by detailed
                                           balance)
  KSTEP_S_<1s> e + Kr*(<1s>) -> Kr+ + 2e  (bands kstep_s, kstep_p, 0: no
  KSTEP_P   e + Kr** -> Kr+ + 2e           measurement)

A level's row is named by its group and its Paschen label as Biagi prints
it, lower case, with Biagi's '!' marks (levels sharing a label) written
_1, _2, ... ('3D1!!' is C1b_3d1_2). group_rows() gathers a group's rows.

v0_1 to v0_3's optional group 'electron' (KSTEP and SUPER, on the lumped
Kr*) is gone: v0_4 carries no optional group.

Units: cm^3/s (two-body), cm^6/s (three-body), 1/s (radiative).

Rate specs (declarative, evaluated by model.py):
  ("const", k)                      constant
  ("xsec", gas, kind, product)      Maxwellian <sigma v>(Te) of one LXCat block
  ("assoc", gas, kind, product)     v0_2: the Hornbeck-Molnar branch of an
                                    excitation, the level's own <sigma v>
                                    per electron and per atom excited. The
                                    row's second Kr is the ground-state
                                    atom the excited one ionizes on, whose
                                    density the band's share stands for, so
                                    the rate is first order in Kr
                                    (model.py). Under the kinetic EEDF the
                                    electron pays the level's threshold
                                    and the new electron is born at rest
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
  ("pop_xsec", "jung2006", level[, "reverse", g_ratio, upper, lower])
                                    the 5s -> 5p excitation out of one 5s
                                    level on Jung et al.'s measured cross
                                    sections (jung2006.py): its Maxwellian
                                    table, and the reverse its
                                    Klein-Rosseland table over the ledger's
                                    gap; under the kinetic EEDF a population
                                    process of the Boltzmann equation
  ("eq_reverse", forward, fwd_pow, molecule, g_atoms[, d0 band])
                                    v0_4: the reverse of an association by
                                    detailed balance, k_f (Tg/300)^fwd_pow /
                                    K_eq(Tg) (thermo.py) on the ledger's
                                    binding energy (shifted by the band)
  ("level_reverse", forward, fwd_pow, g_ratio)
                                    v0_4: the reverse of a transfer between
                                    two levels by detailed balance at the
                                    gas temperature, k_f (Tg/300)^fwd_pow
                                    g_ratio exp(-dE/kTg)

A row's tg_fac (v0_4), (band, ((Tg, F), ...)), multiplies its rate by
F(Tg)^x, x the band's value, F log-log interpolated between the nodes.

A row's `frac` multiplies its rate by band values: (("name", "x"), ...)
by the band's value, (("name", "1-x"), ...) by one minus it.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import xsec
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
    # v0_4: (band, ((Tg [K], F), ...)): the rate times F(Tg)^x, F read
    # log-log between the points and held at their ends, x the band's value
    tg_fac: tuple = ()

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
BASIC_GROUPS = ()            # v0_4: no optional group (optional_rows)
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
# printed value at once (run_v0_4.py --ladder kannari), a sensitivity
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
# v0_4: the 5p -> 5s radiative branching. Einstein A [1/s] of every 2p ->
# 1s line, NIST Atomic Spectra Database (Kramida, Ralchenko, Reader & NIST
# ASD Team, doi:10.18434/T4W30F; Kr I lines 740-900 nm, read 2026-10-15),
# accuracy grades AAA to B+. Weighted by 2J + 1 over the ten 5p levels
# (the statistical distribution v0_1 to v0_3's lifetime assumes), their
# mean decay is 1/29.5 ns, against the measured lifetimes' 29.0 ns, and
# the lumped rung's radiation splits over the 5s levels as KR_5P_BRANCH
KR_5P_5S_A = {
    ("2p1", "1s2"): 4.064e7, ("2p2", "1s2"): 3.416e7,
    ("2p3", "1s2"): 1.418e7, ("2p3", "1s3"): 2.041e7,
    ("2p4", "1s2"): 1.811e7, ("2p4", "1s3"): 1.583e7,
    ("2p5", "1s4"): 4.310e7,
    ("2p6", "1s4"): 8.940e6, ("2p6", "1s5"): 2.732e7,
    ("2p7", "1s4"): 2.931e7, ("2p7", "1s5"): 4.270e6,
    ("2p8", "1s4"): 2.217e7, ("2p8", "1s5"): 8.960e6,
    ("2p9", "1s5"): 3.610e7, ("2p10", "1s5"): 2.289e7,
}
_A_W = {lv: sum(KR_5P_G[u] * a for (u, low), a in KR_5P_5S_A.items()
                if low == lv) for lv in ("1s5", "1s4", "1s3", "1s2")}
KR_5P_BRANCH = {lv: w / sum(_A_W.values()) for lv, w in _A_W.items()}


# ---------------------------------------------------------------------------
# The basic network: krypton and electrons
# ---------------------------------------------------------------------------
# v0_2 (2026-10-14): the ground-state excitations level by level. A level
# below the 5p rung is a 5s level (G1, product Kr*); the 2p levels are the
# 5p manifold (C1a, Kr**); the 4d and 6s levels below C1B_TOP_EV are C1b;
# every level above is C1c. C1b and C1c cascade to the 5p rung, taken
# prompt, the excess over Kr**'s energy heating the gas, as v0_1's d+p and
# 5d groups did. C1B_TOP_EV lies in the 0.37-eV gap between the 6s levels
# (2s4, 12.3853 eV) and the 6p (3p10, 12.7564 eV)
C1B_TOP_EV = 12.5
# the excited atoms that can ionize on a ground-state Kr, Kr(nl) + Kr ->
# Kr2+ + e: those at or above Kr2+'s energy, IP(Kr) - D0(Kr2+) (species.py,
# D0 measured by photoionization of Kr2, Ng et al. 1977), with thermal
# collision energy neglected. In Biagi 8.97 that is the twenty levels from
# 3s1_3 (12.8573 eV) up; 3s1_4 at 12.8253 eV lies 0.024 eV below
E_ASSOC_EV = e_form("Kr2+")


def level_tag(product: str) -> str:
    """'Kr(3D1!!)' -> '3d1_2', 'Kr(4F SUM' -> '4f', 'Kr(HIGH)' -> 'high'."""
    lab = product.strip()
    if lab.startswith("Kr("):
        lab = lab[3:]
    lab = lab.rstrip(")").replace(" SUM", "").strip()
    bangs = lab.count("!")
    lab = lab.replace("!", "").lower()
    return f"{lab}_{bangs}" if bangs else lab


def level_group(p) -> str:
    """The group a Biagi 8.97 excitation level is booked in."""
    if p.threshold < e_form("Kr**"):
        return "G1"
    if level_tag(p.product).startswith("2p"):
        return "C1a"
    return "C1b" if p.threshold < C1B_TOP_EV else "C1c"


def levels() -> list:
    """The excitation blocks of the vendored file, in its order (rising
    threshold)."""
    return [p for p in xsec.processes("Kr") if p.kind == "excitation"]


_PRODUCT = {"C1a": "Kr**", "C1b": "Kr**", "C1c": "Kr**"}
# v0_4: the four 5s levels, each its own species (species.py), by
# statistical weight 2J + 1 and the 5p -> 5s radiative branching (PRAD)
LEVELS_5S = ("1s5", "1s4", "1s3", "1s2")
G_LEVEL_5S = {"1s5": 5, "1s4": 3, "1s3": 1, "1s2": 3}


def kr5s(level: str) -> str:
    """The species of a 5s level, 'Kr*(1s5)'."""
    return f"Kr*({level})"


KR_5S = tuple(kr5s(lv) for lv in LEVELS_5S)
_WHAT = {
    "G1": "a 5s level, its own species since v0_4 (v0_1 to v0_3 lumped "
          "the four at the 1s5 energy)",
    "C1a": "a 5p level; Kr** is the 5p rung at the 2p10 energy, the excess "
           "over it heats the gas. It relaxes to Kr* by KPS, radiates "
           "(PRAD) or is moved on by the electrons (ladder_rows)",
    "C1b": "a 4d or 6s level; the cascade to the 5p rung (Kr**) is taken "
           "prompt and the excess heats the gas",
    "C1c": "a level from 12.76 eV up (6p, 5d, 7s, ...); the cascade to the "
           "5p rung (Kr**) is taken prompt and the excess heats the gas",
}


def excitation_rows() -> list:
    """The notes' krypton rows and their closures (KrF_model_notes_with_
    Gary.rtf, as LIGHT's v0 wrote them, and LIGHT's basic network since),
    level by level since v0_2: G1_* make the 5s manifold, G2 is Kr* + Kr*
    Penning ionization on its band, G6 the elastic loss, C1a_* to C1c_*
    make the 5p rung Kr** (LIGHT v0_4), HM_* ionize associatively on the
    band hm_kr and C2 ionizes. The notes' G3, G4, G5 and C3 are
    fluorine's and are not here."""
    by_group = {"G1": [], "C1a": [], "C1b": [], "C1c": []}
    hm = []
    for p in levels():
        g = level_group(p)
        tag = level_tag(p.product)
        label = p.product.strip()
        assoc = g == "C1c" and p.threshold >= E_ASSOC_EV
        cite = (f"Biagi (Magboltz 8.97, LXCat 'Biagi') {label} at "
                f"{p.threshold:.4f} eV: {_WHAT[g]}")
        if g == "G1" and tag == "1s5":
            cite = ("notes (1) e + Kr -> Kr* 5s[3/2]. " + cite
                    + "; k = <sigma v> over the EEDF (notes: Q = k n_e[Kr])")
        if assoc:
            cite += (". At or above E(Kr2+): the share hm_kr ionizes "
                     f"associatively instead (HM_{tag})")
        product = kr5s(tag) if g == "G1" else _PRODUCT[g]
        by_group[g].append(Row(
            f"{g}_{tag}", ("e", "Kr"), ("e", product),
            ("xsec", "Kr", "excitation", p.product),
            "gary" if g == "G1" else "closure", cite,
            e_cost_eV=p.threshold,
            frac=(("hm_kr", "1-x"),) if assoc else ()))
        if assoc:
            hm.append(Row(
                f"HM_{tag}", ("e", "Kr", "Kr"), ("e", "e", "Kr2+"),
                ("assoc", "Kr", "excitation", p.product), "basic",
                f"{label} + Kr -> Kr2+ + e, the Hornbeck-Molnar ionization "
                f"of the level Biagi (Magboltz 8.97) excites at "
                f"{p.threshold:.4f} eV, {p.threshold - E_ASSOC_EV:.3f} eV "
                "above Kr2+ (species.py): the share hm_kr of the level's "
                "excitations (band, 0 by default: no measurement for "
                "krypton). The electron pays the threshold, the new "
                "electron is born at rest and the excess over Kr2+'s "
                "energy heats the gas",
                e_cost_eV=p.threshold, frac=(("hm_kr", "x"),)))
    ion = xsec.find("Kr", "ionization")
    return (by_group["G1"] + penning_rows() + [
        Row("G6", ("e", "Kr"), ("e", "Kr"), ("elastic", "Kr"),
            "gary", "notes (6) e + Kr elastic. Biagi (Magboltz 8.97) "
            "momentum-transfer cross section; power to the gas (2m/M)(1 - "
            "Tg/Te) n_e N <sigma_m v eps>")]
        + by_group["C1a"] + by_group["C1b"] + by_group["C1c"] + hm + [
        Row("C2", ("e", "Kr"), ("e", "e", "Kr+"),
            ("xsec", "Kr", "ionization", ion.product), "closure",
            f"Biagi (Magboltz 8.97) ionization ({ion.threshold:.4f} eV): "
            "the notes' 'Kr+' output needs electron-impact ionization",
            e_cost_eV=ion.threshold),
    ])


def penning_rows() -> list:
    """G2, Kr* + Kr* -> Kr+ + Kr + e on the band kr_pen, between every pair
    of 5s levels (v0_4). The lumped row's rate k [Kr*]^2 is kept in total:
    a pair of one level reads k, a pair of two levels 2k (the event counts
    once for each order)."""
    out = []
    for i, a in enumerate(LEVELS_5S):
        for b in LEVELS_5S[i:]:
            out.append(Row(
                f"G2_{a}_{b}", (kr5s(a), kr5s(b)), ("Kr+", "Kr", "e"),
                ("param", "kr_pen", 1.0 if a == b else 2.0), "basic",
                f"Kr*({a}) + Kr*({b}) -> Kr+ + Kr + e, the notes' (2) Kr* "
                "+ Kr* 'associative ionization' on the band kr_pen (5.0e-10 "
                "cm3/s, Kannari 1985 Table I from Leichner & Ericson, PRA 9, "
                "251 (1974); 1.5e-10 in Jancaitis's model from Werner's "
                "1975 thesis), the same coefficient for every pair of 5s "
                "levels: the one thermal measurement found, Kr 3P2-3P2 "
                "(Blagoev et al., J. Phys. Conf. Ser. 44, 80 (2006)), "
                "prints its value only in the full text"
                + ("" if a == b else "; 2k for a pair of two levels")))
    return out


def group_rows(rows, group: str) -> list:
    """The rows of one level group ('G1', 'C1a', 'C1b', 'C1c', 'HM'),
    by their ids."""
    return [r for r in rows if r.rid.startswith(group + "_")]


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
            "the rate fall by 510 K. No measurement above 510 K; v0_4 adds "
            "the reverse, K15r", tg_pow=-0.75),
        Row("K15r", ("Kr2+", "Kr"), ("Kr+", "Kr", "Kr"),
            ("eq_reverse", ("const", KR_CONV_300), -0.75, "Kr2+", 4),
            "basic", "Kr2+ + Kr -> Kr+ + 2Kr, the redissociation of the "
            "molecular ion (v0_4): K15's rate over the equilibrium "
            "constant of Kr+ + Kr <-> Kr2+ (thermo.py), on the measured "
            "D0 = 1.15 eV (Ng et al. 1977; 1.149 eV, Signorell, "
            "Hollenstein & Merkt, JCP 114, 9840 (2001), abstract). Its "
            "1.15 eV comes from the gas. It matters above about 1000 K"),
        Row("R11p", ("e", "Kr2+"), ("Kr**", "Kr"),
            ("te_pow", 1.6e-6, 300.0 * 8.617333262e-5, -0.55), "basic",
            "Kr2+ + e -> Kr** + Kr, 1.6e-6 (300 K/Te)^0.55 cm3/s +- 10 %, "
            "Shiu & Biondi, PRA 16, 1817 (1977) (abstract, Te 300-8400 K, "
            "Tgas 300 K); Mikus, J. Phys. D 11, L39 (1978) measured "
            "1.48e-6 (Te/300 K)^-0.53 over 300-21300 K, within 8 %. The "
            "products are mostly 5p states (Shiu & Biondi observed the 5p "
            "-> 5s lines in the afterglow), so v0_4 makes the 5p rung Kr** "
            "(11.30 eV) and books 1.546 eV of heat per event, where v0_3 "
            "made Kr* and booked 2.93 eV. v0_4: the band dr_tg on a "
            "gas-temperature factor, 0 by default (Shiu & Biondi's Tgas = "
            "300 K law): at its top the rate is times F(Tg), F = 4.80 at "
            "900 K and 2.01 at 2700 K, Cunningham & Hobson's shock-tube "
            "coefficients at Te = Tg (J. Phys. B 5, 1773 (1972), abstract: "
            "4.2e-6 cm3/s at 900 K, falling as T^-3/2 to 9.6e-7 at 2700 K) "
            "over this law at Te = Tg, and 1 at 300 K",
            tg_fac=("dr_tg", ((300.0, 1.0), (900.0, 4.80), (2700.0, 2.01)))),
    ]


# Turner 1967's 1s4 -> 1s5 transfer by Kr, 1.3e-13 cm3/s at 295 K on his
# two-point law (T/295)^-1.29, as the rows' 300-K form takes it
MIX_1S4_300 = 1.3e-13 * (295.0 / 300.0) ** 1.29


def kr2_rows() -> list:
    """The Kr2* route (LIGHT v0_2), level by level since v0_4: each 5s
    level makes the excimer at its own rate (AK8_<level>), Kr2* radiates
    the 147-nm second continuum, and every association has its reverse by
    detailed balance (AK8r_<level>, thermo.py), the redissociation that hot
    gas drives. The fast state Kr2*(0u+), made from Kr*(1s4) on the band
    kr2_fast, radiates the same continuum in 3.4 ns. The 5s levels mix on
    krypton at their measured rates (MIX_*), the two resonance levels
    radiate what escapes the trapping (RES_*, band res_escape), and Kr2*
    may be quenched by Kr (QKR2, band kr2_quench) and ionizes on itself
    (PEN_KR2, band kr2_pen). Unmeasured rows are bands at 0."""
    form_1s5 = ("param", "kr2_form")
    form_1s4 = ("param", "kr2_form_1s4")
    form_1s3 = ("param", "kr2_form_1s3")
    a_note = ("The reverse by detailed balance (thermo.py): the rate over "
              "the equilibrium constant of the excimer's A state on its "
              "measured well (band kr2_d0, the ledger's 0.709 eV to "
              "Killeen & Eden's 0.62 eV), the well drawn from the gas")
    rows = [
        Row("AK8_1s5", (kr5s("1s5"), "Kr", "Kr"), ("Kr2*", "Kr"), form_1s5,
            "basic",
            "Kr*(1s5, 3P2) + 2Kr -> Kr2* + Kr, band kr2_form, default 3.81 "
            "+- 0.08e-32 cm6/s: Marchal, Sewraj, Gardou, Merbahi & Yousfi, "
            "'Study of formation and decay of rare-gas excimers by "
            "laser-induced fluorescence' (IntechOpen 2018, "
            "doi:10.5772/intechopen.71942) Table 2, k32 for pure Kr by "
            "selective excitation (temperature not stated). The band spans "
            "the other direct 1s5 measurements: Ueda et al. 1984 1.7 +- "
            "0.2e-32 (e-beam, 0.5-1 atm; v0_1-v0_3's default), Tracy & "
            "Oskam 1976 2.59e-32, Turner 1967 4e-32 (Phys. Rev. 158, 121; "
            "6.9e-32 at 196 K), Smith & Turner 1963 4.1-4.2e-32. Salamero "
            "1979, Hughes 1977 and Bouciqué & Mortier 1970, in v0_3's band, "
            "do not measure it (2026-10-15 survey). Tg form (Tg/300)^-0.64: "
            "Wieme & Lenaerts' Ea = -0.020 eV, measured over 200-400 K",
            tg_pow=-0.64),
        Row("AK8_1s4", (kr5s("1s4"), "Kr", "Kr"), ("Kr2*", "Kr"), form_1s4,
            "basic",
            "Kr*(1s4, 3P1) + 2Kr -> Kr2* + Kr, band kr2_form_1s4, default "
            "9.24 +- 0.80e-32 cm6/s: Marchal et al. 2018 Table 2, k31 "
            "(pure Kr, selective excitation). Audouard, Laporte, Subtil & "
            "Damany, JCP 89, 6176 (1988) (abstract) found the dimer "
            "formation times from 3P1 and 3P2 in the ratio 3/1, the band's "
            "low end at a third of the 1s5 default, 1.27e-32. Janssens et "
            "al., JCP 86, 4935 (1987) (abstract) give 82.0p^2 /s (p in "
            "mbar), 1.36e-31 at 295 K, from a quasiclassical calculation "
            "compared with their experiments, above the band (a "
            "calculation; its 3P2 term, 3.85e-32, lies inside kr2_form's "
            "band). The share kr2_fast goes to the fast state (AK8f_1s4). "
            "Tg form as 1s5's",
            tg_pow=-0.64, frac=(("kr2_fast", "1-x"),)),
        Row("AK8f_1s4", (kr5s("1s4"), "Kr", "Kr"), ("Kr2*(0u+)", "Kr"),
            form_1s4, "basic",
            "Kr*(1s4) + 2Kr -> Kr2*(0u+) + Kr, the share kr2_fast of the "
            "1s4 association that radiates from the fast B 1Sigma_u+ (0u+) "
            "state (band, 0 by default: all to the slow state, as v0_1 to "
            "v0_3). At its top, 0.17, Audouard et al. 1988's intensity "
            "ratio of the fast to the slow component after 3P1 excitation "
            "above 6e19 cm^-3, F/S = 0.21",
            tg_pow=-0.64, frac=(("kr2_fast", "x"),)),
        Row("AK8_1s3", (kr5s("1s3"), "Kr", "Kr"), ("Kr2*", "Kr"), form_1s3,
            "basic",
            "Kr*(1s3, 3P0) + 2Kr -> Kr2* + Kr, band kr2_form_1s3. Kholin, "
            "L'dov, Ustinovskii & Zayarnyi measured the quenching of 3P0 by "
            "Kr (J. Russ. Laser Res. 36, 74 (2015)), the value not in a "
            "source in hand; the default is the 1s5 rate as a proxy, the "
            "band 0 to the largest measured 5s rate (9.24e-32)",
            tg_pow=-0.64),
        Row("AK8r_1s5", ("Kr2*", "Kr"), (kr5s("1s5"), "Kr", "Kr"),
            ("eq_reverse", form_1s5, -0.64, "Kr2*", 5, "kr2_d0"), "basic",
            "Kr2* + Kr -> Kr*(1s5) + 2Kr. " + a_note),
        Row("AK8r_1s4", ("Kr2*", "Kr"), (kr5s("1s4"), "Kr", "Kr"),
            ("eq_reverse", form_1s4, -0.64, "Kr2*", 3, "kr2_d0"), "basic",
            "Kr2* + Kr -> Kr*(1s4) + 2Kr, AK8_1s4's reverse. " + a_note,
            frac=(("kr2_fast", "1-x"),)),
        Row("AK8fr_1s4", ("Kr2*(0u+)", "Kr"), (kr5s("1s4"), "Kr", "Kr"),
            ("eq_reverse", form_1s4, -0.64, "Kr2*(0u+)", 3), "basic",
            "Kr2*(0u+) + Kr -> Kr*(1s4) + 2Kr, AK8f_1s4's reverse by "
            "detailed balance on the B state's measured well (0.684 eV, "
            "LaRocque et al. 1986)", frac=(("kr2_fast", "x"),)),
        Row("AK8r_1s3", ("Kr2*", "Kr"), (kr5s("1s3"), "Kr", "Kr"),
            ("eq_reverse", form_1s3, -0.64, "Kr2*", 1, "kr2_d0"), "basic",
            "Kr2* + Kr -> Kr*(1s3) + 2Kr, AK8_1s3's reverse. " + a_note),
        Row("KR2_RAD", ("Kr2*",), ("Kr", "Kr"), ("const", 3.7e6), "basic",
            "Kr2* -> 2Kr + hv, 3.7e6 1/s = 270 ns, the slow (1u, 0u-) "
            "component that carries the kinetics: Morikawa et al., JCP 91, "
            "1469 (1989) 270 +- 5 ns; Salamero 1979 250 +- 20, Quigley & "
            "Hughes 1978 280 +- 30, Audouard 1988 265, Marchal et al. 2018 "
            "261.6 +- 5.7, Oka 1974 353 ns, Leichner & Ericson 1974 about "
            "0.30 us at 1000 Torr (abstracts; J&H 1980's 3.3e6 agrees). "
            "REBUILD's KR2_RAD keeps Kannari's 8.0e7, which none of the "
            "measurements supports", hv_eV=8.43, band="147 nm (Kr2* VUV)"),
        Row("KR2F_RAD", ("Kr2*(0u+)",), ("Kr", "Kr"), ("const", 1.0 / 3.4e-9),
            "basic",
            "Kr2*(0u+) -> 2Kr + hv, the fast (0u+) component of the second "
            "continuum, 3.4 +- 0.1 ns (Morikawa et al. 1989, the low-density "
            "limit; abstract)", hv_eV=8.43, band="147 nm (Kr2* VUV)"),
        Row("QKR2", ("Kr2*", "Kr"), ("Kr", "Kr", "Kr"),
            ("param", "kr2_quench"), "basic",
            "Kr2* + Kr -> 3Kr, collisional quenching of the slow excimer: "
            "band kr2_quench, 0 by default (Oka, Rama Rao, Redpath & "
            "Firestone, JCP 61, 4740 (1974), 200-1400 Torr: 'no evidence "
            "for collisional steps'). Firestone & Selander, JCP 88, 717 "
            "(1988) (abstract) measured the 1u, 0u- decay at 984 nm as "
            "2.53e6 + 194(+-84) P_Kr s-1 (Torr, 298 K), 6.0e-15 cm3/s, "
            "which the band's top, 8.6e-15, spans with its error. The well "
            "and the photon's energy heat the gas"),
        Row("MIX_1s4", (kr5s("1s4"), "Kr"), (kr5s("1s5"), "Kr"),
            ("const", MIX_1S4_300), "basic",
            "Kr*(1s4) + Kr -> Kr*(1s5) + Kr, 1.3e-13 cm3/s at 295 K: "
            "Turner, Phys. Rev. 158, 121 (1967) (abstract: 13e-14 at 295 "
            "K, 22e-14 at 196 K, from the decay of the 123.6- and 125-nm "
            "emission). Tg law (Tg/295)^-1.29 through the two points "
            "(carried as its 300-K value, 1.272e-13, on (Tg/300)^-1.29), "
            "held to them in form above 295 K (no measurement there). Its "
            "reverse by detailed balance (MIX_1s5r) agrees with Turner's "
            "measured 1s5 -> 1s4, 7.6e-16 at 295 K, to 2 %", tg_pow=-1.29),
        Row("MIX_1s5r", (kr5s("1s5"), "Kr"), (kr5s("1s4"), "Kr"),
            ("level_reverse", ("const", MIX_1S4_300), -1.29, 3, 5), "basic",
            "Kr*(1s5) + Kr -> Kr*(1s4) + Kr, MIX_1s4's reverse by detailed "
            "balance at the gas temperature, g(1s4)/g(1s5) exp(-0.117 "
            "eV/kTg); the gap comes from the gas"),
        Row("MIX_1s2", (kr5s("1s2"), "Kr", "Kr"), (kr5s("1s4"), "Kr", "Kr"),
            ("const", 0.2e-32), "basic",
            "Kr*(1s2, 1P1) + 2Kr -> Kr*(1s4, 3P1) + 2Kr, about 0.2e-32 "
            "cm6/s: Audouard, Laporte, Subtil & Reininger, JCP 95, 3283 "
            "(1991) (abstract: above 4e20 cm^-3 the 1P1 decay occurs fully "
            "through 3P1, with a reservoir effect of rate constant ~0.2e-32 "
            "for the transfer). The 0.611 eV heats the gas"),
        Row("MIX_1s4r2", (kr5s("1s4"), "Kr", "Kr"), (kr5s("1s2"), "Kr", "Kr"),
            ("level_reverse", ("const", 0.2e-32), 0.0, 3, 3), "basic",
            "Kr*(1s4) + 2Kr -> Kr*(1s2) + 2Kr, MIX_1s2's reverse by "
            "detailed balance, g(1s2)/g(1s4) exp(-0.611 eV/kTg)"),
        Row("RES_1s4", (kr5s("1s4"),), ("Kr",), ("param", "res_escape"),
            "basic",
            "Kr*(1s4) -> Kr + hv, 123.6 nm, the resonance radiation that "
            "escapes the trapping (A = 2.98e8 /s untrapped, NIST ASD): band "
            "res_escape, 0 by default (a 0-D box has no size). Its top, "
            "1.7e5 /s, is the pressure-independent term of the 3P1 loss "
            "rate Janssens, Vanmarcke, Desoppere, Bouciqué & Wieme give "
            "with diffusion and resonance radiation, 1.7e5 + 2050p + "
            "82.0p^2 /s (p in mbar; JCP 86, 4935 (1987), abstract: a "
            "quasiclassical theory compared with their alpha-excitation "
            "experiments)",
            hv_eV=e_form(kr5s("1s4")), band="123.6 nm (resonance)"),
        Row("RES_1s2", (kr5s("1s2"),), ("Kr",), ("param", "res_escape"),
            "basic",
            "Kr*(1s2) -> Kr + hv, 116.5 nm, the escaping resonance "
            "radiation (A = 3.09e8 /s untrapped, NIST ASD), band "
            "res_escape as RES_1s4's", hv_eV=e_form(kr5s("1s2")),
            band="116.5 nm (resonance)"),
        Row("PEN_KR2", ("Kr2*", "Kr2*"), ("Kr2+", "Kr", "Kr", "e"),
            ("param", "kr2_pen"), "basic",
            "Kr2* + Kr2* -> Kr2+ + 2Kr + e, no measurement, so 0 unless a "
            "band run sets it (kr2_pen, 0 to 5.0e-10 cm3/s). Jancaitis's "
            "model carries 1.0e-10 from Werner et al., IEEE J. Quant. Elect. "
            "13, 769 (UCRL-53465 Table II p.86 and ref. 1 p.109). Kannari "
            "1985 Table I p.4312 [87] (REBUILD PEN_KR2). Ref. 87 is Lorents "
            "& Olson, SRI semiannual technical report N00014-72-C-0557 "
            "(1973, unpublished), and Kannari gives Ar2* + Ar2* the same "
            "5.0e-10 from it (read first-hand 2026-09-29). The 2026-10-15 "
            "survey found no krypton measurement (Xe2* + Xe2*: 8e-11, Zamir "
            "et al., APL 27, 56 (1975))"),
    ]
    return rows


def kannari_k76_rows() -> list:
    """LIGHT v0_4's K76 pair, Kannari's Table I [76] fit A Te^0.5
    exp(-C/Te) and its detailed-balance reverse, for Params.k76 =
    'kannari' (a comparison configuration; the band k_5s5p scales it as it
    scales Jung's). LIGHT_Kr v0_4: the same fit out of every 5s level, its
    reverse over each level's gap."""
    out = []
    for lv in LEVELS_5S:
        gap = e_form("Kr**") - e_form(kr5s(lv))
        out += [
            Row(f"K76_{lv}", ("e", kr5s(lv)), ("e", "Kr**"),
                ("te_fit",) + KANNARI_K76, "basic",
                f"e + Kr*({lv}) -> e + Kr(5p), Kannari 1985 Table I [76], "
                "9.0e-7 Te^0.5 exp(-1.51/Te) cm3/s, Jacob & Mangano's fit "
                "to e-beam KrF data (a fit, not a direct measurement), the "
                f"same out of every 5s level. The {gap:.3f} eV comes from "
                "the electron gas", frac=(("k_5s5p", "x"),)),
            Row(f"K76r_{lv}", ("e", "Kr**"), ("e", kr5s(lv)),
                ("te_fit",) + KANNARI_K76
                + ("reverse", G_LEVEL_5S[lv] / G_5P, "Kr**", kr5s(lv)),
                "basic",
                f"e + Kr(5p) -> e + Kr*({lv}), the reverse of Kannari's K76 "
                f"by detailed balance, g({lv})/g(5p) = "
                f"{G_LEVEL_5S[lv]}/{G_5P} over the ledger's gap; its "
                f"{gap:.3f} eV returns to the electron gas",
                energy_to="electrons", frac=(("k_5s5p", "x"),)),
        ]
    return out


def ladder_rows() -> list:
    """LIGHT v0_4's 5p rung Kr** and its rows, on the resolved 5s levels
    since LIGHT_Kr v0_4. Kr** is made by the higher excitations (C1a-c) and
    by Kr2+ + e (R11p). It relaxes to the 5s levels in collisions with Kr
    (KPS_*: the measured total loss, its products over the levels on the
    band kps_5s), radiates to them on the measured branching (PRAD_*,
    near-IR lines that are not trapped), and the electrons pump it from
    each 5s level (K76_*, Jung et al.'s cross sections), return it (K76r_*,
    detailed balance) and ionize it (KSTEP_P); KSTEP_S_* ionize the 5s
    levels. The two stepwise ionizations have no measurement and are bands
    at 0 (BANDS). Params.ladder = 'kannari' sets every ladder band to
    Kannari's printed value at once."""
    rows = []
    kps_cite = ("Kr(5p) + Kr -> Kr(5s) + Kr, the p -> s relaxation. Band "
                "kps_kr, default (1.7 +- 0.2)e-10 cm3/s: the one Kr "
                "measurement in hand, the total two-body loss of Kr "
                "5p'[3/2]2 (2p2) in Kr at 300 K by ps-TALIF (Gazeli et al., "
                "Phys. Plasmas 28, 043301 (2021), Stern-Volmer over 0.1-10 "
                "mbar; the product channel is not resolved). The band runs "
                "down to Jancaitis's 3.0e-11, the Ar analogue (UCRL-53465 "
                "Table II, from Chang & Setser 1978), past Kannari's "
                "1.0e-10 'estimated' for Ar (REBUILD KPS_Kr). v0_4: the "
                "products, which nothing measures, go over the 5s levels "
                "by statistical weight (the share kps_5s, 1 by default) or "
                "all to 1s5, the lowest (the rest). The gap heats the gas")
    for lv in LEVELS_5S:
        rows.append(Row(
            f"KPS_{lv}", ("Kr**", "Kr"), (kr5s(lv), "Kr"),
            ("param", "kps_kr", G_LEVEL_5S[lv] / G_5S), "basic",
            f"{kps_cite}. This row: the statistical share to Kr*({lv}), "
            f"{G_LEVEL_5S[lv]}/{G_5S}", frac=(("kps_5s", "x"),)))
    rows.append(Row(
        "KPS_low", ("Kr**", "Kr"), (kr5s("1s5"), "Kr"),
        ("param", "kps_kr"), "basic",
        f"{kps_cite}. This row: the share 1 - kps_5s, all to Kr*(1s5)",
        frac=(("kps_5s", "1-x"),)))
    for lv in LEVELS_5S:
        rows.append(Row(
            f"PRAD_{lv}", ("Kr**",), (kr5s(lv),),
            ("const", KR_5P_BRANCH[lv] / KR_5P_TAU_S), "basic",
            f"Kr(5p) -> Kr*({lv}) + hv, the near-infrared 5p -> 5s lines "
            "(760-893 nm), which end on the 5s levels and are not trapped. "
            "One lifetime for the rung, 29.0 ns, the statistical-weight "
            "mean of the ten measured 5p radiative rates (23 to 41 ns; "
            "Dzierzega et al. 2000 Table 1; KR_5P_TAU_NS), times the "
            f"rung's branching to {lv}, {KR_5P_BRANCH[lv]:.3f}, from the "
            "NIST ASD A-values weighted by 2J + 1 (KR_5P_5S_A). The "
            "photon carries the ledger's gap from the rung's 2p10 energy",
            hv_eV=e_form("Kr**") - e_form(kr5s(lv)),
            band="5p-5s (near IR)"))
    for lv in LEVELS_5S:
        gap = e_form("Kr**") - e_form(kr5s(lv))
        proxy = lv in ("1s4", "1s2")
        how = ("the J = 1 level's own cross section is not measured: its "
               "same-core metastable's set at its own energy (1s4 the "
               "1s5 set, 1s2 the 1s3 set; jung2006.py), scaled by the band "
               "k_5s5p_j1 as well" if proxy else
               f"the measured cross sections out of {lv} into the ten 5p "
               "levels")
        frac = (("k_5s5p", "x"),) + ((("k_5s5p_j1", "x"),) if proxy else ())
        rows += [
            Row(f"K76_{lv}", ("e", kr5s(lv)), ("e", "Kr**"),
                ("pop_xsec", "jung2006", lv), "basic",
                f"e + Kr*({lv}) -> e + Kr(5p) on {how} (Jung, Stone, "
                "Boffard, Anderson & Lin, PRA 73, 022722 (2006) Tables II-"
                "IV; PRL 94, 163202 (2005)). Band k_5s5p is a scale on the "
                "measurement, 1 by default, 0.6 to 1.4 (its +-35 to 40 %). "
                f"The electron pays the ledger's gap, {gap:.3f} eV, the "
                "rung booked at its 2p10 energy", frac=frac),
            Row(f"K76r_{lv}", ("e", "Kr**"), ("e", kr5s(lv)),
                ("pop_xsec", "jung2006", lv, "reverse",
                 G_LEVEL_5S[lv] / G_5P, "Kr**", kr5s(lv)), "basic",
                f"e + Kr(5p) -> e + Kr*({lv}), the reverse of K76_{lv} by "
                "the Klein-Rosseland relation on the same cross section, "
                f"g({lv})/g(5p) = {G_LEVEL_5S[lv]}/{G_5P} over the ledger's "
                f"gap; its {gap:.3f} eV returns to the electron gas",
                energy_to="electrons", frac=frac),
        ]
    for lv in LEVELS_5S:
        c_ev = KANNARI_KSTEP_S[2] - (e_form(kr5s(lv)) - e_form(kr5s("1s5")))
        rows.append(Row(
            f"KSTEP_S_{lv}", ("e", kr5s(lv)), ("Kr+", "e", "e"),
            ("te_fit", "kstep_s", KANNARI_KSTEP_S[1], c_ev), "basic",
            f"e + Kr*({lv}) -> Kr+ + 2e, stepwise ionization. Band kstep_s: "
            "no measurement (the 2026-10-15 survey found none; one 1977 "
            "attempt failed); Kannari 1985 Table I [79] prints Vriens's "
            "binary-encounter calculation 7.0e-8 Te^0.68 exp(-4.04/Te) for "
            f"the lumped 5s, here with its exponent at {lv}'s own "
            f"ionization energy, {c_ev:.3f} eV. REBUILD's template runs "
            "2.2x above it at 1.3 eV, the band's top"))
    rows.append(Row(
        "KSTEP_P", ("e", "Kr**"), ("Kr+", "e", "e"),
        ("te_fit", "kstep_p") + KANNARI_KSTEP_P[1:], "basic",
        "e + Kr(5p) -> Kr+ + 2e. Band kstep_p: no measurement; Kannari 1985 "
        "Table I [79], Vriens's calculation 2.0e-7 Te^0.62 exp(-2.53/Te). "
        "The 2.70 eV comes from the electron gas"))
    return rows


def basic_rows() -> list:
    """The pure-Kr network, in LIGHT v0_5's order: the excitations, the
    charge chemistry, the Kr2* route and the 5p rung."""
    return excitation_rows() + charge_rows() + kr2_rows() + ladder_rows()


# ---------------------------------------------------------------------------
# The optional group 'electron' of v0_1 to v0_3 (none in v0_4)
# ---------------------------------------------------------------------------
def optional_rows() -> list:
    """v0_1 to v0_3's optional group 'electron' (REBUILD's stepwise
    template KSTEP and the lumped superelastic SUPER) read the lumped Kr*,
    which v0_4 resolves: it carries no optional group. The stepwise
    ionization of each 5s level is KSTEP_S_* on its band; the superelastic
    return of the 5s levels, measured in an afterglow (Blagoev, Mishonov &
    Popov, J. Phys. B 17, 435 (1984): 1.5e-10 to 9.8e-10 cm3/s for thermal
    electrons), is not carried."""
    return []


# ---------------------------------------------------------------------------
# Bands (LIGHT v0_2, 2026-09-29): the rates known only as a published
# spread, or not at all
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Band:
    """A rate the basic network carries as a band. Rows read it through a
    ("param", name[, frac]) rate, a ("te_fit", name, ...) rate or a frac
    factor, and run_v0_4.py --bands runs each end. A measured rate keeps
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
    # --- LIGHT_Kr v0_4 (2026-10-15): the 5s levels resolved -----------------
    Band("kr2_form", "Kr*(1s5) + 2Kr -> Kr2* + Kr at 300 K", 3.81e-32,
         1.7e-32, 4.2e-32, "cm6/s",
         "Marchal et al. 2018 (doi:10.5772/intechopen.71942) Table 2, "
         "3.81 +- 0.08e-32 by selective excitation in pure Kr, the default. "
         "The band spans the other direct 1s5 measurements: Ueda et al. "
         "1984 1.7e-32 (e-beam at 0.5-1 atm, the low end), Tracy & Oskam "
         "1976 2.59e-32, Turner 1967 4e-32, Smith & Turner 1963 4.1-4.2e-32 "
         "(the top)",
         marks=((1.7e-32, "Ueda 1984"), (2.59e-32, "Tracy & Oskam 1976"),
                (4.0e-32, "Turner 1967"))),
    Band("kr2_form_1s4", "Kr*(1s4) + 2Kr -> Kr2* + Kr at 300 K", 9.24e-32,
         1.27e-32, 9.24e-32, "cm6/s",
         "Marchal et al. 2018 Table 2, 9.24 +- 0.80e-32 (3P1, selective "
         "excitation), the default and top. Audouard et al. 1988 found the "
         "formation times from 3P1 and 3P2 in the ratio 3/1: a third of "
         "the 1s5 default, 1.27e-32, the low end",
         marks=((1.27e-32, "Audouard 1988, 3P2 rate / 3"),)),
    Band("kr2_form_1s3", "Kr*(1s3) + 2Kr -> Kr2* + Kr at 300 K", 3.81e-32,
         0.0, 9.24e-32, "cm6/s",
         "No value in a source in hand (Kholin et al. 2015 measured the "
         "quenching of 3P0 by Kr; its number is not readable here). The "
         "default is the 1s5 rate as a proxy; 0 to the largest measured "
         "5s rate"),
    Band("kr2_fast", "the share of the Kr*(1s4) association that forms the "
         "fast 0u+ state", 0.0, 0.0, 0.17, "fraction",
         "0 by default: every excimer in the slow (1u, 0u-) state, as v0_1 "
         "to v0_3. Audouard et al. 1988: after 3P1 excitation the fast "
         "component appears from 1.5e19 cm^-3 and stays at F/S = 0.21 "
         "above 6e19, F/(F + S) = 0.17, the top"),
    Band("kr2_quench", "Kr2* + Kr -> 3Kr", 0.0, 0.0, 8.6e-15, "cm3/s",
         "0 by default: Oka et al. 1974 saw no collisional step at 200-1400 "
         "Torr. Firestone & Selander 1988 measured 194 +- 84 Torr^-1 s^-1 "
         "at 298 K, 6.0e-15 cm3/s, the band to its upper error, 8.6e-15",
         marks=((6.0e-15, "Firestone & Selander 1988"),)),
    Band("kr2_d0", "the excimer's well for its redissociation (thermo.py)",
         0.709, 0.62, 0.709, "eV",
         "LaRocque et al. 1986: D_e = 5807.0 cm^-1 of the A state, less half "
         "its low-v quantum, 0.709 eV, the default and the ledger's well "
         "(species.py). Killeen & Eden 1986 estimate D0 = 0.62 +- 0.02 eV, "
         "the low end. The band moves the redissociation, not the ledger"),
    Band("dr_tg", "the gas-temperature factor on Kr2+ + e (R11p)", 0.0, 0.0,
         1.0, "fraction",
         "0 by default: Shiu & Biondi's coefficient, measured at Tgas = 300 "
         "K. 1: Cunningham & Hobson's shock-tube coefficients at Te = Tg "
         "(4.2e-6 cm3/s at 900 K, 9.6e-7 at 2700 K) over that law, a "
         "factor 4.80 at 900 K and 2.01 at 2700 K (network.R11p)"),
    Band("res_escape", "the escaping share of the resonance radiation, a "
         "rate per Kr*(1s4) and Kr*(1s2)", 0.0, 0.0, 1.7e5, "1/s",
         "0 by default: a 0-D box has no size to escape from. 1.7e5 /s is "
         "the pressure-independent term of the 3P1 loss rate Janssens et "
         "al. (JCP 86, 4935 (1987), abstract) give with diffusion and "
         "resonance radiation for their cell, a quasiclassical theory "
         "compared with their alpha-excitation experiments"),
    Band("kps_5s", "the share of Kr(5p) + Kr whose products go over the 5s "
         "levels by statistical weight, the rest to Kr*(1s5)", 1.0, 0.0,
         1.0, "fraction",
         "No measurement of the products (Gazeli et al. 2021 measured the "
         "total loss). 1, statistical, by default; 0, all to the lowest "
         "level, the slowest excimer formation"),
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
         "cross sections out of every 5s level (jung2006.py)", 1.0, 0.6,
         1.4, "factor",
         "Jung et al., PRA 73, 022722 (2006): the cross sections out of "
         "the two metastables, +-35 to 40 %, the band. Kannari's fit is "
         "Params.k76 = 'kannari'"),
    Band("k_5s5p_j1", "a scale on the J = 1 levels' 5s -> 5p cross "
         "sections, their same-core metastable's sets (jung2006.py)", 1.0,
         0.0, 1.0, "factor",
         "Nothing measured out of 1s4 and 1s2: 1 by default, the proxy; 0, "
         "the J = 1 levels pumped by no electron (v0_3's k_5s5p = 0.5 end)"),
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
    # --- LIGHT_Kr v0_3 (2026-10-14) -----------------------------------------
    Band("hm_kr", "the share of the excitations of the levels at or above "
         "E(Kr2+) (Biagi 8.97's twenty, 12.857 to 13.6 eV) that ionize "
         "associatively on a ground-state atom, Kr(nl) + Kr -> Kr2+ + e "
         "(Hornbeck-Molnar), where the rest cascade to the 5p rung", 0.0,
         0.0, 1.0, "fraction",
         "No measurement for krypton of the share against the levels' "
         "radiation and quenching. Hornbeck & Molnar, Phys. Rev. 84, 621 "
         "(1951) found Kr2+ made by electrons below the atomic ionization "
         "energy, which shows the process and not its probability at a "
         "pressure. 0 by default, as every unmeasured row; 1, every such "
         "excitation ionizing, is the top, a probability's own bound"),
)}
# bands no row reads: Params.rates may still set them (model.py)
FREE_BANDS = ("ib_scale",)
# v0_4: bands at 0 by default although a measurement bounds them, each for
# the reason its cite gives (a measurement at another density, measurements
# that disagree, a law measured at 300 K, no geometry in a 0-D box); the
# report lists them apart from the unmeasured ones
OFF_BY_DEFAULT = ("dr_tg", "kr2_fast", "kr2_quench", "res_escape")


def build(groups=(), k76="jung") -> list:
    """The pure-Kr network (v0_4 carries no optional group; 'all' is
    accepted and adds none). k76: 'jung' (the measured 5s -> 5p cross
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
        missing = _row_bands(r) - set(BANDS)
        if missing:
            raise NetworkError(f"{r.rid}: no band {sorted(missing)}")
        if r.rate[0] == "pop_xsec" and r.rate[1] != "jung2006":
            raise NetworkError(f"{r.rid}: pop_xsec source {r.rate[1]!r}")
        for name, how in r.frac:
            if name not in BANDS or how not in ("x", "1-x"):
                raise NetworkError(f"{r.rid}: frac ({name!r}, {how!r})")
    return rows


def _row_bands(r) -> set:
    """The bands one row reads: its 'param' or 'te_fit' rate, the forward
    rate and the well of a detailed-balance reverse (v0_4), its `frac`
    factors and its gas-temperature factor (v0_4)."""
    out = set()
    spec = r.rate
    if spec[0] in ("eq_reverse", "level_reverse"):
        if spec[1][0] == "param":
            out.add(spec[1][1])
        if spec[0] == "eq_reverse" and len(spec) > 5:
            out.add(spec[5])
    elif spec[0] == "param" or (spec[0] == "te_fit"
                                and isinstance(spec[1], str)):
        out.add(spec[1])
    out.update(name for name, _ in r.frac)
    if r.tg_fac:
        out.add(r.tg_fac[0])
    return out


def bands_read(rows) -> set:
    """The bands the rows read, through 'param' and 'te_fit' rates, the
    reverses' forward rates and wells, `frac` factors and gas-temperature
    factors."""
    out = set()
    for r in rows:
        out |= _row_bands(r)
    return out
