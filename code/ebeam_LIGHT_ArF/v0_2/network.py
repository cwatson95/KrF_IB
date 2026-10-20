"""The ArF basic network AS DATA (LIGHT_ArF v0_2, 2026-09-30).

The Ar/F2 version of LIGHT/v0_2/network.py: one rare gas, argon, no
buffer. It keeps KrF v0_2's design. Every rate is a measurement, cited in
its row. A rate known only as a published spread is a band (BANDS) that
keeps its measurement as the default. A rate with no measurement is a band
with default 0, so it is off unless a band run sets it. The first-hand
checks of the numbers are in docs/arf_data/SOURCES_CHECKED.md, and the
research reading behind them is docs/arf_data/RATES_RESEARCH.md.

What differs from KrF, and why:
  * ArF* is one state, ArF(B). No C position, B-C mixing rate or
    equilibrium is measured; the calculated B/C splitting is about 0.25 eV
    near R_e (Krauss & Mies, Excimer Lasers 2nd ed. 1984, p.15), so C holds
    ~1e-4 of the manifold at 300 K. Every measurement-based ArF model
    (Rokni, Chen & Payne, SRI, Jancaitis) lumps them the same way.
  * The ArF(B) radiative lifetime is not measured (4.2 ns is Dunning &
    Hay's calculation). The quench measurements are Stern-Volmer products
    k tau, so the rows carry k tau and divide by Params.tau_sp_s: changing
    the calculated lifetime keeps every measured product.
  * Ar* + F2 forms ArF* with yield about 0.5, measured (Setser 1977), not
    the unity every published ArF model books.
  * The electron-impact set is IST-Lisbon's: 37 excitation blocks, each its
    own row with its own threshold (the 4s rows G1_*, the higher rows C1_*
    cascading promptly to the 4s manifold).

Units: cm^3/s (two-body), cm^6/s (three-body), 1/s (radiative).

Rate specs (declarative, evaluated by model.py):
  ("const", k)                      constant
  ("inv_tau", param[, frac])        frac / Params.<param>: 1/tau, or with
                                    param 'tau_sp_s' and frac a measured k
                                    tau product, k = (k tau) / tau_B
  ("xsec", gas, kind, product)      <sigma v>(Te) of one LXCat block
  ("elastic", gas)                  G6: energy loss only
  ("stepwise", E_bind)              1.56e-7 Te^0.71 exp(-E/Te)  (REBUILD)
  ("te_pow", k0, Te0, n)            k0 (Te/Te0)^n, Te held above 0.02 eV;
                                    a captured electron takes (3/2 + n) Te
  ("hm", A, B, C[, n3])             Hoffman & Moreno's ion-ion form
                                    A N / (N^2 + B N + C / f3), N in amagat,
                                    times the Debye screening S; f3 =
                                    (Tg/300)^n3 acts on the three-body limit
                                    A N f3 / C (Thomson theory, n3 = -3)
  ("superelastic", gas, product, g_ratio)
                                    de-excitation from the Klein-Rosseland
                                    cross section of an LXCat excitation
  ("param", name[, frac[, mods]])   frac x the value of band `name`
                                    (Params.rates if set there, else the
                                    band's default), and each modifier in
                                    mods: ("share", b) multiplies by band b,
                                    ("rest", b) by 1 - band b, ("per_tau",)
                                    divides by Params.tau_sp_s
"""

from __future__ import annotations

from dataclasses import dataclass

from .species import (ELECTRON, HC_EV_NM, atoms, charge, e_form,
                      E_AR_STAR)


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
    """Per-event energy flows [eV] (the REBUILD routing rule, as KrF):

    * radiative: the quantum escapes, the rest (lower-state repulsion) heats;
    * electron-impact with an LXCat threshold: the electron pays the
      threshold, the products store dE_form, the difference heats (the
      cascade of the upper Ar levels to 4s, F2 vibration);
    * exothermic: to the electron gas when the row ejects a new electron from
      heavy reactants (Penning) or says energy_to='electrons'
      (superelastic), otherwise to the gas;
    * endothermic without a threshold: from the electron gas if an electron
      collides, otherwise from the gas.
    The captured electron's energy moves electron gas -> gas at run time
    ('capture')."""
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
# Electron impact on Ar: the IST-Lisbon set, one row per block
# ---------------------------------------------------------------------------
IST = ("IST-Lisbon (Yanguas-Gil, Cotrino & Alves, J. Phys. D 38, 1588 "
       "(2005); LoKI-MC/Code/Input/Argon/Ar_LXCat.txt, md5-pinned in "
       "xsec.py)")
# (row id, LXCat product label, threshold eV), in threshold order; the test
# holds this table to the vendored file
AR_LEVELS = (
    ("G1_3P2", "Ar(3P2)", 11.548),
    ("G1_3P1", "Ar(3P1)", 11.623),
    ("G1_3P0", "Ar(3P0)", 11.723),
    ("G1_1P1", "Ar(1P1)", 11.828),
    ("C1_4p01", "Ar(4p[1/2]1)", 12.907),
    ("C1_4p02", "Ar(4p[5/2]3)", 13.076),
    ("C1_4p03", "Ar(4p[5/2]2)", 13.095),
    ("C1_4p04", "Ar(4p[3/2]1)", 13.153),
    ("C1_4p05", "Ar(4p[3/2]2)", 13.172),
    ("C1_4p06", "Ar(4p[1/2]0)", 13.273),
    ("C1_4p07", "Ar(4p'[3/2]1)", 13.283),
    ("C1_4p08", "Ar(4p'[3/2]2)", 13.302),
    ("C1_4p09", "Ar(4p'[1/2]1)", 13.328),
    ("C1_4p10", "Ar(4p'[1/2]0)", 13.48),
    ("C1_3d01", "Ar(3d[1/2]0 + 3d[1/2]1)", 13.84),
    ("C1_3d02", "Ar(3d[3/2]2)", 13.9),
    ("C1_3d03", "Ar(3d[7/2]4)", 13.98),
    ("C1_3d04", "Ar(3d[7/2]3)", 14.01),
    ("C1_3d05", "Ar(3d[5/2]2 + 5s[3/2]2)", 14.06),
    ("C1_3d06", "Ar(3d[5/2]3 + 5s[3/2]1)", 14.09),
    ("C1_3d07", "Ar(3d[3/2]1)", 14.15),
    ("C1_3d08", "Ar(3d'[5/2]2)", 14.21),
    ("C1_3d09", "Ar(3d'[5/2or3/2]3or2+5s'[1/2]0or1)", 14.23),
    ("C1_3d10", "Ar(3d'[3/2]1)", 14.3),
    ("C1_5p01", "Ar(5p[1/2]1)", 14.464),
    ("C1_5p02", "Ar(5p[5/2]3)", 14.499),
    ("C1_5p03", "Ar(5p[5/2]2)", 14.506),
    ("C1_5p04", "Ar(5p[3/2]1)", 14.525),
    ("C1_5p05", "Ar(5p[3/2]2)", 14.529),
    ("C1_5p06", "Ar(5p[1/2]0)", 14.564),
    ("C1_5p07", "Ar(5p'[3/2]1)", 14.681),
    ("C1_5p08", "Ar(5p'[1/2]1)", 14.687),
    ("C1_5p09", "Ar(5p'[3/2]2)", 14.688),
    ("C1_5p10", "Ar(5p'[1/2]0)", 14.738),
    ("C1_4d", "Ar(4d)", 14.78),
    ("C1_6s", "Ar(6s)", 14.842),
    ("C1_4dp", "Ar(4d')", 14.967),
)
CASCADE = ("cascade to the 4s manifold taken prompt, as KrF's C1 rows: "
           "Ar(4p) relaxes by Ar at 6e-11 to 1e-10 cm3/s (Chang & Setser "
           "1978 via NRL's Orestes memo; Kannari's estimate), about 1e9 1/s "
           "at 1 atm, and Setser's ERDA report COO-2807-11 (1977) p.-4- "
           "finds more than half of the Ar(4p) quenching goes directly to "
           "Ar(4s); the excess over 11.548 eV heats the gas")


def electron_rows() -> list:
    rows = []
    for rid, label, thr in AR_LEVELS:
        if rid.startswith("G1"):
            what = (f"e + Ar -> Ar({label[3:-1]}), a 4s level lumped into "
                    f"Ar* at 11.548 eV")
            if thr > E_AR_STAR:
                what += (f"; the {thr - E_AR_STAR:.3f} eV above 3P2 heats "
                         "the gas (collisional mixing within 4s)")
            if label in ("Ar(3P1)", "Ar(1P1)"):
                what += ("; a resonance level, its 107/105-nm radiation "
                         "trapped at these pressures")
        else:
            what = f"e + Ar -> {label}, {CASCADE}"
        rows.append(Row(rid, ("e", "Ar"), ("e", "Ar*"),
                        ("xsec", "Ar", "excitation", label), "basic",
                        f"{what}. {IST} block '{label}', threshold {thr} "
                        "eV (the 4s resonance cross sections are Khakoo's "
                        "x0.5 there)", e_cost_eV=thr))
    rows += [
        Row("C2", ("e", "Ar"), ("e", "e", "Ar+"),
            ("xsec", "Ar", "ionization", "Ar+"), "basic",
            f"e + Ar -> 2e + Ar+ (15.76 eV), {IST}, Rapp & Englander-Golden's "
            "measurement", e_cost_eV=15.76),
        Row("G6", ("e", "Ar"), ("e", "Ar"), ("elastic", "Ar"), "basic",
            f"e + Ar elastic, {IST} momentum-transfer set (Phelps' "
            "effective, via Yamabe et al.); power to the gas (2m/M)(1 - "
            "Tg/Te) n_e N <sigma_m v eps>"),
    ]
    return rows


SIGLO = "SIGLO set (the F2 file LIGHT and REBUILD read)"


def f2_electron_rows() -> list:
    """Electron collisions with F2, as KrF's basic network carries them."""
    return [
        Row("ATT", ("e", "F2"), ("F-", "F"),
            ("xsec", "F2", "attachment", None), "basic",
            f"dissociative attachment, {SIGLO}"),
        Row("DIS_a3", ("e", "F2"), ("e", "F", "F"),
            ("xsec", "F2", "excitation", "F2(a3)(3.16eV)"), "basic",
            "SIGLO a3Pi_u (repulsive: excitation = dissociation, REBUILD "
            "0.1.14 reroute)", e_cost_eV=3.16),
        Row("DIS_A1", ("e", "F2"), ("e", "F", "F"),
            ("xsec", "F2", "excitation", "F2(A1)(4.34eV)"), "basic",
            "SIGLO A1Pi_u (repulsive, as DIS_a3)", e_cost_eV=4.34),
        Row("VIB1", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V1)(0.111eV)"), "basic",
            "SIGLO v=1, the quantum relaxes to heat", e_cost_eV=0.111),
        Row("VIB2", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V2)(0.219eV)"), "basic",
            "SIGLO v=2, the quantum relaxes to heat", e_cost_eV=0.219),
        Row("VIB3", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V3)(0.324eV)"), "basic",
            "SIGLO v=3, the quantum relaxes to heat", e_cost_eV=0.324),
        Row("VIB4", ("e", "F2"), ("e", "F2"),
            ("xsec", "F2", "excitation", "F2(V4)(0.421eV)"), "basic",
            "SIGLO v=4, the quantum relaxes to heat", e_cost_eV=0.421),
    ]


# ---------------------------------------------------------------------------
# The measured constants (first-hand checks: docs/arf_data/SOURCES_CHECKED.md)
# ---------------------------------------------------------------------------
# ArF(B) radiative lifetime: not measured. 4.2 ns is Dunning & Hay's
# calculation (JCP 69, 134 (1978), via Jancaitis, UCRL-53465 p.97, and
# Brau's Table 4.2, Excimer Lasers 2nd ed. p.90). Rokni and Chen & Payne
# reduced their Stern-Volmer products with 4 ns.
ARF_TAU_B = 4.2e-9
# Stern-Volmer products k tau (the measured quantities)
KT_ARF_F2 = 7.55e-18        # cm3: Rokni 7.6 +- 0.7, Chen & Payne 7.5 (mean)
KT_ARF_AR = 3.6e-20         # cm3: Rokni 3.6 +- 1
# the Ar2F* constants, each the mean of the three first-hand measurements
AR2F_TAU = (185e-9 + 205e-9 + 220e-9) / 3.0          # 203 ns
K_AR2F_F2 = (2.1e-10 + 1.79e-10 + 1.9e-10) / 3.0      # 1.93e-10 cm3/s
AR_CONV_300 = 2.2e-31       # Ar+ + 2Ar -> Ar2+ + Ar, cm6/s at 300 K
DR_AR2 = 7.5e-7             # Ar2+ + e, cm3/s at Te = 300 K
FF_AR_300 = 8.0e-35         # F + F + Ar -> F2 + Ar, cm6/s, measured
TE_300K = 300.0 * 8.617333262e-5
# the photons: band peaks measured by Chen & Payne, IEEE JQE 15, 149 (1979)
# p.150 (1933 A ArF*, 2840 A Ar2F*, 1270 A Ar2*(1u))
HV_BX = HC_EV_NM / 193.3
HV_AR2F = HC_EV_NM / 284.0
HV_AR2 = HC_EV_NM / 127.0
BX_BAND = "193 nm (ArF B-X)"
TRIMER_BAND = "284 nm (Ar2F*)"
VUV_BAND = "127 nm (Ar2* VUV)"
TRIMER_RIDS = ("G5",)
HARPOON_RIDS = ("G3",)

RK77 = ("Rokni, Jacob & Mangano, PRA 16, 2216 (1977) Table II p.2218, "
        "e-beam Stern-Volmer, read first-hand")
CP79 = ("Chen & Payne, IEEE JQE 15, 149 (1979) Table I p.155, "
        "proton-excited, read first-hand")
JANC = "Jancaitis, UCRL-53465 (1983) ch. V, read first-hand"
UEDA84 = ("Ueda et al., Rev. Laser Eng. 12, 357 (1984) Table V p.365, "
          "pulse radiolysis, read first-hand")
BRAU = "Brau, in Excimer Lasers 2nd ed. (Rhodes ed., Springer 1984)"


def basic_rows() -> list:
    """The ArF basic network: electron impact on Ar and F2, the harpoon, the
    ArF(B) block, the trimer, the ions, and the Ar2* route."""
    rows = electron_rows() + f2_electron_rows()
    rows += [
        # --- Ar* ----------------------------------------------------------
        Row("G2", ("Ar*", "Ar*"), ("Ar+", "Ar", "e"), ("param", "ar_pen"),
            "basic", "Ar* + Ar* -> Ar+ + Ar + e. No measurement among the "
            "sources held, so 0 unless a band run sets it (ar_pen, 0 to "
            "6.2e-10 cm3/s). Kannari 1985's 5.0e-10 cites Leichner & "
            "Ericson 1974, a Kr paper, and Lee & Lieberman's global model "
            "uses 6.2e-10 (via NRL's Orestes memo)"),
        Row("G3", ("Ar*", "F2"), ("ArF(B)", "F"),
            ("param", "ar_f2", 1.0, (("share", "arf_yield"),)), "basic",
            "Ar* + F2 -> ArF(B) + F: the ArF* share (band arf_yield, 0.5) "
            "of the total Ar(3P2) + F2 rate (band ar_f2, 8.1e-10 cm3/s). "
            "The yield is Setser's flowing-afterglow measurement, 'the "
            "branching ratio is ~0.5 rather than unity' (ERDA COO-2807-11 "
            "(1977) pp.-5-/-6-, read first-hand); Kolts et al. 1978 give "
            "4.0e-10 for the ArF* channel (via Ueda 1984 Table III), 0.53 "
            "of Velazco's total, and Chen & Payne's relative f1/f2 = 0.55 "
            "+- 0.25 agrees. Every published ArF model books all of it to "
            "ArF* (the band's top)"),
        Row("G3X", ("Ar*", "F2"), ("Ar", "F", "F"),
            ("param", "ar_f2", 1.0, (("rest", "arf_yield"),)), "basic",
            "Ar* + F2 without ArF*: the rest of the total (see G3). The "
            "products are not measured; dissociation to Ar + 2F, its energy "
            "to heat, is the assignment (as KrF's G3X)"),
        Row("AK8_Ar", ("Ar*", "Ar", "Ar"), ("Ar2*", "Ar"),
            ("param", "ar2_form"), "basic",
            "Ar*(3P2) + 2Ar -> Ar2* + Ar, (6.0 +- 0.5)e-33 cm6/s by Ar* "
            "absorption at 1-4 atm, " + UEDA84 + " (also Ueda 1983 Table "
            "III). Band ar2_form carries the low-pressure measurements "
            "Brau tabulates, 2.8e-33 (Keto 1974) to 3.3e-32 (Kolts & "
            "Setser 1977) (" + BRAU + " Table 4.7 p.102), and Chen & "
            "Payne's 1.1e-32. No gas-temperature law is measured"),
        # --- ArF(B) ---------------------------------------------------------
        Row("G4", ("ArF(B)",), ("Ar", "F"), ("inv_tau", "tau_sp_s"),
            "basic",
            "ArF(B) -> Ar + F + hv (B-X), 1/tau_B with tau_B = 4.2 ns, "
            "Dunning & Hay's calculation (JCP 69, 134 (1978), via "
            "Jancaitis p.97 and " + BRAU + " Table 4.2 p.90); no "
            "measurement is among the sources held (McCown & Greene, "
            "LA-UR-90-291, still planned one in 1990). Band peak 1933 A, "
            "measured (Chen & Payne p.150)",
            hv_eV=HV_BX, band=BX_BAND),
        Row("QF2", ("ArF(B)", "F2"), ("Ar", "F", "F", "F"),
            ("inv_tau", "tau_sp_s", KT_ARF_F2), "basic",
            "ArF(B) + F2, k tau_B = 7.55e-18 cm3, the mean of the measured "
            "Stern-Volmer products (7.6 +- 0.7)e-18 (" + RK77 + ") and "
            "7.5e-18 (" + CP79 + ", beta3 = 1.9e-9 +- 25 % at 4 ns); k = "
            "1.80e-9 cm3/s at tau_B = 4.2 ns"),
        Row("QAR", ("ArF(B)", "Ar"), ("Ar", "Ar", "F"),
            ("inv_tau", "tau_sp_s", KT_ARF_AR), "basic",
            "ArF(B) + Ar, k tau_B = (3.6 +- 1)e-20 cm3 (" + RK77 + "); k "
            "= 8.6e-12 cm3/s at tau_B = 4.2 ns"),
        Row("G5", ("ArF(B)", "Ar", "Ar"), ("Ar2F*", "Ar"),
            ("param", "arf_trimer", 1.0, (("per_tau",),)), "basic",
            "ArF(B) + 2Ar -> Ar2F* + Ar, k tau_B = 2.27e-39 cm6 (band "
            "arf_trimer), the mean of the three first-hand Stern-Volmer "
            "products: (1.6 +- 0.3)e-39 (" + RK77 + "), 2.1e-39 (" + CP79
            + ") and 3.1e-39 +- 25 % (" + JANC + " p.150); k = 5.4e-31 "
            "cm6/s at tau_B = 4.2 ns. Marowsky et al. 1982 give (1.2 +- "
            "0.2)e-30 (via Huestis, Marowsky & Tittel, Excimer Lasers 2nd "
            "ed. Table 6.4 p.191), 5.0e-39 at 4.2 ns, the band's top. No "
            "gas-temperature law is measured, so none is applied (band "
            "trimer_tg_pow)"),
        Row("QE", ("ArF(B)", "e"), ("Ar", "F", "e"),
            ("param", "arf_e", 1.0, (("per_tau",),)), "basic",
            "ArF(B) + e, k_e tau_B = 5.75e-16 cm3 (band arf_e): Duzy et al., "
            "Avco Everett final report 1980 (OSTI 5268869) Table III p.-23- "
            "measured 1.6e-7 and 3.0e-7 cm3/s at 44.4 and 10.1 A/cm2, 0.25 "
            "atm, read first-hand. Their Stern-Volmer slope is k_e tau_rad "
            "/ k_att, reduced with k_att = 4.5e-9 cm3/s (stated for their "
            "KrF data) and, for ArF*, tau = 2.5 ns (p.-21-), so the measured "
            "products k_e tau are 4.0e-16 and 7.5e-16 cm3 (the band) and "
            "k_e = 1.37e-7 cm3/s at tau_B = 4.2 ns. Trainor & Jacob 1980's "
            "1.6e-7 (via Lee 1989) is from the same laboratory"),
        # --- Ar2F* ----------------------------------------------------------
        Row("C3", ("Ar2F*",), ("Ar", "Ar", "F"), ("const", 1.0 / AR2F_TAU),
            "basic",
            "Ar2F* -> 2Ar + F + hv, tau = 203 ns, the mean of the three "
            "first-hand measurements: 185 +- 25 % (" + CP79 + "), 205 +- 20 "
            "(" + JANC + " Table I p.158) and 220 ns (" + UEDA84 + "). "
            "Marowsky 1982's 219 +- 15 and Boewering 1982's 230 ns (via "
            "Huestis Table 6.4) lie inside; Wadt & Hay calculated 132 ns. "
            "Band peak 2840 A, measured (Chen & Payne p.150)",
            hv_eV=HV_AR2F, band=TRIMER_BAND),
        Row("Q2F", ("Ar2F*", "F2"), ("Ar", "Ar", "F", "F2"),
            ("const", K_AR2F_F2), "basic",
            "Ar2F* + F2, 1.93e-10 cm3/s, the mean of the three first-hand "
            "measurements: 2.1e-10 +- 25 % (" + CP79 + "), (1.79 +- "
            "0.06)e-10 (" + JANC + " Table I p.158) and (1.9 +- 0.1)e-10 ("
            + UEDA84 + "). Marowsky 1982's (2.05 +- 0.06)e-10 (via "
            "Jancaitis) disagrees with Jancaitis's beyond both errors"),
        Row("Q2FF", ("Ar2F*", "F"), ("Ar", "Ar", "F", "F"),
            ("param", "ar2f_f"), "basic",
            "Ar2F* + F, 1.45e-10 cm3/s (band ar2f_f): Jancaitis found the "
            "trimers quenched after the pulse at the rate the initial F2 "
            "would give with the F2 nearly all dissociated, so 'the "
            "quenching rate coefficient of the RGH trimers by fluorine "
            "atoms is between .5 and 1 times the corresponding value for "
            "F2' (" + JANC + " p.153). The band is that range times Q2F, the "
            "default its middle. The products are not measured: "
            "dissociation to 2Ar + 2F, its energy to heat"),
        Row("Q2AR", ("Ar2F*", "Ar"), ("Ar", "Ar", "Ar", "F"),
            ("param", "ar2f_ar"), "basic",
            "Ar2F* + Ar: below 1e-14 cm3/s, measured (" + JANC + " Table I "
            "p.158, and Marowsky 1982 via it), so 0 unless a band run sets "
            "it; Boewering et al. 1982 give (2.2 +- 0.8)e-14 (via Huestis "
            "Table 6.4), the band's top"),
        # --- ions ------------------------------------------------------------
        Row("K15", ("Ar+", "Ar", "Ar"), ("Ar2+", "Ar"),
            ("const", AR_CONV_300), "basic",
            "Ar+ + 2Ar -> Ar2+ + Ar, 2.2e-31 cm6/s at 300 K, the median of "
            "the eight measurements " + BRAU + " Table 4.8 p.104 lists (5e-32 "
            "Kebarle to 4.4e-31 Knewstubb; Liu & Conway's drift-tube values "
            "2.07 and 2.3e-31, Smith & Cromey 2.5e-31). Tg law (Tg/300)^"
            "-0.75, the classical theory of Perez-Rios & Greene, PRA 98, "
            "062707 (2018) for rare-gas ion three-body recombination (as "
            "KrF's K15)", tg_pow=-0.75),
        Row("R11p", ("e", "Ar2+"), ("Ar*", "Ar"),
            ("te_pow", DR_AR2, TE_300K, -0.67), "basic",
            "Ar2+ + e -> Ar* + Ar, 7.5e-7 (Te/300 K)^-0.67 cm3/s (Oskam & "
            "Mittelstadt 1963, Mehr & Biondi 1968, via " + BRAU + " Table "
            "4.9 p.105). The products are 4p states (McCusker, Excimer "
            "Lasers ch. 3 pp.49-50), which cascade to 4s as C1 assumes; "
            "the 2.98 eV difference heats the gas. No gas-temperature law "
            "is measured"),
        Row("K26", ("Ar+", "F-"), ("ArF(B)",),
            ("hm", 9.32e-6, 6.21e-2, 1.69, -3.0), "basic",
            "Ar+ + F- (+Ar) -> ArF(B), 9.32e-6 N/(N^2 + 6.21e-2 N + 1.69) "
            "cm3/s at 300 K, N in amagat: Hoffman & Moreno's modified "
            "Natanson calculation (SAND80-1486) as Lee et al., IEEE Trans. "
            "Electron Devices 36, 2053 (1989) Table I p.2055 print it (read "
            "first-hand), with the Debye screening and the Thomson (Tg/300)"
            "^-3 on the three-body limit that KrF's K26 carries. No "
            "absolute measurement: 'has not been directly measured' ("
            + BRAU + " p.98). Rokni, Jacob & Mangano 1977 found by changing "
            "the beam current that the ions form ArF* 'with unit branching' "
            "(p.2219, Fig. 5)"),
        Row("K108", ("Ar2+", "F-"), ("ArF(B)", "Ar"),
            ("hm", 1.51e-5, 1.81, 2.25, -3.0), "basic",
            "Ar2+ + F- (+Ar) -> ArF(B) + Ar, 1.51e-5 N/(N^2 + 1.81 N + "
            "2.25) cm3/s at 300 K, the same calculation and printing as "
            "K26. Rokni et al. 1977: 'almost all the Ar2+ forms ArF*' "
            "(p.2219)"),
        Row("FF", ("F", "F", "Ar"), ("F2", "Ar"), ("const", FF_AR_300),
            "basic",
            "F + F + Ar -> F2 + Ar, (8.0 +- 0.6)e-35 cm6/s at 295 K, "
            "Ganguli & Kaufman 1974 (via Baulch et al., JPCRD 10 Suppl. 1 "
            "(1981) p.1-22, as KrF's FF cites it); d[F2]/dt = k[F]^2[M]"),
        # --- the Ar2* route ----------------------------------------------------
        Row("AR2_RAD", ("Ar2*",), ("Ar", "Ar"), ("param", "ar2_rad"),
            "basic",
            "Ar2* -> 2Ar + hv, 1/3.2 us (band ar2_rad): the slow (1u) "
            "component that carries the kinetics, Keto, Gleason & Walters, "
            "PRL 33, 1365 (1974), 3.2 +- 0.3 us (via " + BRAU + " Table "
            "4.10 p.105; Thonnard & Hurst 2.8, Boucique & Mortier 3.7 us). "
            "The singlet (0u+) radiates in 4.2 ns; the band's top, 6.0e7 "
            "1/s, is Kannari's full singlet-triplet mixing average of those "
            "two measurements. Band peak 1270 A, measured (Chen & Payne "
            "p.150)", hv_eV=HV_AR2, band=VUV_BAND),
        Row("NK7", ("Ar2*", "F2"), ("ArF(B)", "Ar", "F"),
            ("param", "ar2_f2", 1.0, (("share", "ar2_f2_arf"),)), "basic",
            "Ar2*(1u) + F2 -> ArF(B) + Ar + F: the ArF* share (band "
            "ar2_f2_arf, 0.91) of 5.2e-10 cm3/s (band ar2_f2; Chen et al., "
            "via " + CP79 + " and Table II). The share is Setser's 0.5 over "
            "Chen & Payne's measured f1/f2 = 0.55 +- 0.25 (0.62 to 1)"),
        Row("K2F_BR", ("Ar2*", "F2"), ("Ar2F*", "F"),
            ("param", "ar2_f2", 1.0, (("rest", "ar2_f2_arf"),)), "basic",
            "Ar2*(1u) + F2 -> Ar2F* + F, the rest of ar2_f2 (see NK7). Chen "
            "& Payne find this direct channel minor, f2'/f3 < 0.3 (Table "
            "I); SRI MP-76-99 Table 2 assigned all of its 2.5e-10 to it"),
        Row("K2F_F", ("Ar2*", "F"), ("ArF(B)", "Ar"), ("param", "ar2_f"),
            "basic",
            "Ar2* + F -> ArF(B) + Ar, no measurement, so 0 unless a band "
            "run sets it (ar2_f, 0 to 1e-9 cm3/s): SRI MP-76-99 p.65's "
            "conditional 'if the reaction ... is sufficiently rapid (3-10 x "
            "10^-10)'"),
        Row("PEN_AR2", ("Ar2*", "Ar2*"), ("Ar2+", "Ar", "Ar", "e"),
            ("param", "ar2_pen"), "basic",
            "Ar2* + Ar2* -> Ar2+ + 2Ar + e, no measurement, so 0 unless a "
            "band run sets it (ar2_pen, 0 to 5e-10 cm3/s: Lorents & Olson's "
            "unpublished SRI reports, via Kannari 1985 and McCusker)"),
    ]
    return rows


def electron_group_rows() -> list:
    """Optional group 'electron' (off by default): electron collisions with
    Ar*, estimates only, so their cost can be measured."""
    rows = [Row("KSTEP", ("e", "Ar*"), ("Ar+", "e", "e"),
                ("stepwise", round(15.76 - E_AR_STAR, 3)), "electron",
                "Ar* + e -> Ar+ + 2e, REBUILD's Chen R9 template at E_bind "
                "4.212 eV, EST. The IST-Lisbon file has no stepwise "
                "ionization, and the literature fits differ x100 at Te = 1 "
                "eV (Kannari 6.80e-9 Te^0.67 exp(-4.20/Te), Lee 5.18e-6 "
                "Te^0.67 exp(-6.27/Te))")]
    for rid, label, thr in AR_LEVELS[:4]:
        rows.append(Row(
            "SUPER_" + rid[3:], ("e", "Ar*"), ("e", "Ar"),
            ("superelastic", "Ar", label, 1.0 / 12.0), "electron",
            f"Ar* + e -> Ar + e through {label}: detailed balance of "
            f"{rid} with the 4s levels populated by their statistical "
            "weights (g 5, 3, 1, 3 of 12), so each block's weight over its "
            "own g is 1/12, EST", energy_to="electrons"))
    return rows


BASIC_GROUPS = ("electron",)
GROUPS = BASIC_GROUPS


# ---------------------------------------------------------------------------
# Bands: the rates known only as a published spread, or not at all
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Band:
    """A rate the basic network carries as a band. Rows read it through a
    ("param", name, ...) rate, and run_v0_2.py --bands runs each end. A
    measured rate keeps its measurement as the default and spans the other
    published values. A rate with no measurement is 0 by default, and its
    band runs from 0 to the largest published value, which is no bound.
    Device data and analogs never set a default or narrow a band."""
    name: str
    what: str
    default: float
    lo: float
    hi: float
    unit: str
    cite: str
    marks: tuple = ()          # ((value, who), ...) published values inside


BANDS = {b.name: b for b in (
    Band("ar_f2", "Ar* + F2, total", 8.1e-10, 6.4e-10, 1.1e-9, "cm3/s",
         "The mean of the three groups' Ar(3P2) removal measurements: "
         "Velazco, Kolts & Setser 1976 7.5e-10 (flowing afterglow, a lower "
         "limit by their caveat), Chen et al. 9.4e-10 (via Chen & Payne "
         "Table II) and Ueda 1984 (7.4 +- 1.0)e-10 (absorption, 1-4 atm). "
         "Ueda 1984's ArF* fluorescence decay gives (6.4 +- 0.3)e-10, the "
         "bottom, and Ueda 1983's absorption at 2.5 atm (1.1 +- 0.1)e-9, "
         "the top",
         marks=((7.5e-10, "Velazco 1976"), (9.4e-10, "Chen et al."))),
    Band("arf_yield", "ArF* yield of Ar* + F2", 0.5, 0.25, 1.0, "fraction",
         "Setser 1977 measured about 0.5 (Kolts et al. 1978: 4.0e-10, 0.53 "
         "of Velazco's total). Chen & Payne measured the resonance levels "
         "3P1 and 1P1 below 0.25, the bottom; every published ArF model "
         "books 1.0, the top",
         marks=((0.53, "Kolts et al. 1978"),)),
    Band("ar2_form", "Ar* + 2Ar -> Ar2* + Ar at 300 K", 6.0e-33, 2.8e-33,
         3.3e-32, "cm6/s",
         "Ueda 1984 measured (6.0 +- 0.5)e-33 at 1-4 atm, the default. "
         "Brau's Table 4.7 spans 2.8e-33 (Keto 1974) to 3.3e-32 (Kolts & "
         "Setser 1977); Chen & Payne use 1.1e-32",
         marks=((1.1e-32, "Chen & Payne"), (1.2e-32, "Thonnard & Hurst"))),
    Band("ar2_rad", "Ar2* -> 2Ar + hv", 1.0 / 3.2e-6, 1.0 / 3.2e-6, 6.0e7,
         "1/s",
         "The 1u (triplet) 3.2 us of Keto 1974, the default. Full "
         "singlet-triplet mixing (Kannari's 6.0e7, 0.25/4.2 ns + 0.75/3.2 "
         "us) is the top; SRI's effective value with F2 is 3.8e6",
         marks=((3.8e6, "SRI 1976, effective"),)),
    Band("ar2_f2", "Ar2* + F2, both channels", 5.2e-10, 2.5e-10, 5.2e-10,
         "cm3/s",
         "Chen et al. measured 5.2e-10 (via Chen & Payne), the default. SRI "
         "MP-76-99 gives 2.5e-10 from the Ar2* and Ar2F* decays"),
    Band("ar2_f2_arf", "ArF* share of Ar2* + F2", 0.91, 0.62, 1.0,
         "fraction",
         "f2 = f1 / (f1/f2) with Setser's f1 = 0.5 and Chen & Payne's "
         "measured f1/f2 = 0.55 +- 0.25; the rest forms Ar2F* + F"),
    Band("ar2_f", "Ar2* + F -> ArF(B) + Ar", 0.0, 0.0, 1.0e-9, "cm3/s",
         "No measurement. SRI MP-76-99's conditional 3-10e-10 (p.65) sets "
         "the top",
         marks=((3.0e-10, "Kannari 1985"),)),
    Band("ar2_pen", "Ar2* + Ar2* -> Ar2+ + 2Ar + e", 0.0, 0.0, 5.0e-10,
         "cm3/s",
         "No measurement. Kannari's 5.0e-10 is Lorents & Olson's "
         "unpublished SRI report"),
    Band("ar_pen", "Ar* + Ar* -> Ar+ + Ar + e", 0.0, 0.0, 6.2e-10, "cm3/s",
         "No measurement among the sources held. Kannari's 5.0e-10 cites a "
         "Kr paper; Lee & Lieberman's global model uses 6.2e-10",
         marks=((5.0e-10, "Kannari 1985"),)),
    Band("arf_trimer", "ArF* + 2Ar -> Ar2F* + Ar, k tau_B",
         (1.6e-39 + 2.1e-39 + 3.1e-39) / 3.0, 1.6e-39, 5.0e-39, "cm6",
         "The mean of Rokni 1977 (1.6 +- 0.3), Chen & Payne 1979 (2.1) and "
         "Jancaitis 1983 (3.1 +- 25 %) e-39 cm6. Marowsky 1982's (1.2 +- "
         "0.2)e-30 cm6/s, 5.0e-39 at 4.2 ns, is the top",
         marks=((2.1e-39, "Chen & Payne"), (3.1e-39, "Jancaitis"))),
    Band("trimer_tg_pow", "ArF* + 2Ar -> Ar2F* + Ar, its (Tg/300)^n", 0.0,
         -3.9, 0.0, "exponent",
         "No measurement: the rate is held at its 300 K value. NRL's "
         "Orestes form (Petrova et al., NRL/MR/6720 2022 p.25, k28, "
         "'Shui-1979, Chen-1978', read first-hand), 1.62e-32 + 2.245e-21 "
         "Tg^-3.9 cm6/s, sets the bottom as a pure power. "
         "Params.trimer_tg_pow sets it"),
    Band("arf_e", "ArF* + e, k_e tau_B", 5.75e-16, 4.0e-16, 7.5e-16, "cm3",
         "Duzy et al. 1980's two measurements, 1.6e-7 and 3.0e-7 cm3/s "
         "with their tau = 2.5 ns, as the products they measured; the "
         "default is their mean"),
    Band("ar2f_f", "Ar2F* + F -> 2Ar + 2F", 0.75 * K_AR2F_F2,
         0.5 * K_AR2F_F2, K_AR2F_F2, "cm3/s",
         "Jancaitis 1983: between 0.5 and 1 times the F2 rate (Q2F), the "
         "default its middle"),
    Band("ar2f_ar", "Ar2F* + Ar", 0.0, 0.0, 2.2e-14, "cm3/s",
         "Jancaitis 1983 and Marowsky 1982 measured below 1e-14; Boewering "
         "1982 gives (2.2 +- 0.8)e-14, the top",
         marks=((1.0e-14, "Jancaitis, Marowsky: upper bound"),)),
)}

# Tg ranges the rows' laws were measured over, for the report's warnings:
# (row id, K, text)
TG_VALIDITY = (
    ("G5", 400.0, "ArF* + 2Ar is measured near 300 K and its gas-"
     "temperature law is not measured: the model holds it at its 300 K "
     "value. NRL's form, 1.62e-32 + 2.245e-21 Tg^-3.9 (after Shui's "
     "calculation), keeps 0.07 of it at 700 K, and the band's pure "
     "(Tg/300)^-3.9 (--trimer-tg-pow -3.9) 0.036"),
    ("AK8_Ar", 400.0, "Ar* + 2Ar and the other Ar* and Ar2* rates are "
     "measured near 300 K, with no gas-temperature law"),
)


def bands_of(r: Row) -> tuple:
    """The band names a 'param' row reads (its own and any modifier's)."""
    if r.rate[0] != "param":
        return ()
    names = [r.rate[1]]
    mods = r.rate[3] if len(r.rate) > 3 else ()
    names += [m[1] for m in mods if m[0] in ("share", "rest")]
    return tuple(names)


def param_rate(spec: tuple, value, tau_sp_s: float) -> float:
    """The constant of a ("param", name[, frac[, mods]]) rate, with
    value(name) the band's value for this run."""
    frac = float(spec[2]) if len(spec) > 2 else 1.0
    k = frac * float(value(spec[1]))
    for m in (spec[3] if len(spec) > 3 else ()):
        if m[0] == "share":
            k *= float(value(m[1]))
        elif m[0] == "rest":
            k *= 1.0 - float(value(m[1]))
        elif m[0] == "per_tau":
            k /= float(tau_sp_s)
        else:
            raise NetworkError(f"unknown param modifier {m!r}")
    return k


def build(groups=(), net="basic") -> list:
    """The ArF basic network plus optional groups ('all' = every group)."""
    if net != "basic":
        raise NetworkError(f"unknown network {net!r}: LIGHT_ArF has only "
                           "'basic' (the notes' list 'v0' is KrF's)")
    groups = GROUPS if "all" in groups else tuple(groups)
    bad = set(groups) - set(GROUPS)
    if bad:
        raise NetworkError(f"groups {sorted(bad)} not available; have "
                           f"{GROUPS}")
    rows = basic_rows()
    if "electron" in groups:
        rows += electron_group_rows()
    seen = set()
    for r in rows:
        if r.rid in seen:
            raise NetworkError(f"duplicate row id {r.rid}")
        seen.add(r.rid)
        check_row(r)
        for b in bands_of(r):
            if b not in BANDS:
                raise NetworkError(f"{r.rid}: no band {b!r}")
    return rows
