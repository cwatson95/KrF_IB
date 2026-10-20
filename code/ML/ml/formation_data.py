"""Measurements of the KrF formation chain, link by link (ml/formation.py).

Every entry names its source and where in it the number sits, and says how
it was read:
  checked     re-read first-hand in the source for this file
  secondhand  quoted by another paper (named), the original not read
  abstract    from an abstract only
Laser performance (Tahoma) only scores, as everywhere in ML. Paths named
sources/... are in LIGHT/docs/krf_data/sources/ (checktext CLEAN).
"""
from __future__ import annotations

from dataclasses import dataclass

SETSER_1976 = ("Velazco, Kolts & Setser, preprint in D. W. Setser, ERDA "
               "report COO-2807-2 (1976), OSTI 7361313, LIGHT/docs/krf_data/"
               "sources/Setser_ERDA_COO-2807-2_1976_OSTI7361313.pdf")
UEDA_1984 = ("Ueda, Kanada, Sasaki & Takuma, Rev. Laser Eng. 12, 357 (1984), "
             "LIGHT/docs/krf_data/sources/Ueda_RevLaserEng12_357_1984.pdf")
SRI_1976 = ("Nakano, Hill, Lorents, Huestis & McCusker, SRI MP-76-99 (Dec "
            "1976), DTIC AD-A033997, LIGHT/docs/krf_data/sources/SRI_MP-76-"
            "99_Nakano_1976_DTIC_ADA033997_text.pdf")
# supplied by the owner 2026-09-28 (LIGHT/docs/krf_data/sources/)
CHEN_PAYNE_1979 = ("Chen & Payne, IEEE J. Quantum Electron. 15, 149 (1979), "
                   "sources/Kinetics_of_proton-excited_Ar-F2and_Kr-F2mixtures"
                   ".pdf")
MEZYK_1991 = ("Mezyk, Cooper & Sherwell, J. Phys. Chem. 95, 3152 (1991), "
              "sources/mezyk1991.pdf")
VKS_1976 = ("Velazco, Kolts & Setser, J. Chem. Phys. 65, 3468 (1976), "
            "sources/velazco1976.pdf")
# supplied by the owner 2026-09-29 (Kinetic_Modeling/Recent_KrF_Papers/)
QH_1978 = ("Quigley & Hughes, Appl. Phys. Lett. 32, 649 (1978), "
           "Recent_KrF_Papers/")
JACOB_1978 = ("Jacob, Rokni, Mangano & Brochu, Appl. Phys. Lett. 32, 109 "
              "(1978), Recent_KrF_Papers/")
ROKNI_1977 = ("Rokni, Jacob & Mangano, Phys. Rev. A 16, 2216 (1977), "
              "Recent_KrF_Papers/")
JANCAITIS_1983 = ("Jancaitis, 'Kinetic Modelling of Krypton Fluoride "
                  "Laser Systems', Ph.D. thesis, LLNL UCRL-53465 (1983), "
                  "Recent_KrF_Papers/janacaitis.pdf (supplied 2026-09-29)")
VELAZCO_1977 = ("Velazco, 'A study of the total and reactive electronic "
                "quenching of the lowest metastable state of xenon', Ph.D. "
                "dissertation, Kansas State University (1977), sources/"
                "A_STUDY_OF_THE_TOTAL_AND_REACT.pdf (a scan: PDF page = "
                "printed page + 11)")


@dataclass(frozen=True)
class Value:
    link: str            # a Chain field
    value: float
    err: float           # 1 sigma, same units; 0 = not stated
    source: str
    where: str
    status: str          # checked | secondhand | abstract
    note: str = ""


VALUES = [
    # Kr* + F2, total rate [cm3/s]
    Value("k_H", 7.2e-10, 0.0, SETSER_1976, "Table 3 (Kr 3P2 + F2: 72, "
          "units 1e-11 cm3/s)", "checked", "published in Velazco, Kolts & "
          "Setser, JCP 65, 3468 (1976) Table II (72). Their two F2 "
          "mixtures gave 61 and 72 (p.3473), the mixture being 'the most "
          "likely source of systematic error'"),
    Value("k_H", 8.8e-10, 0.5e-10, UEDA_1984, "Table V, KrF* fluorescence",
          "checked", "the KrF* rise, i.e. the Kr* removal rate"),
    Value("k_H", 9.4e-10, 1.0e-10, UEDA_1984, "Table V, Kr*(3P2) absorption",
          "checked"),
    Value("k_H", 6.8e-10, 0.0, CHEN_PAYNE_1979, "p.157 (Figs. 22-23): 'The "
          "quenching rate constants of Kr(3P1), Kr(1P1), and Kr2*(1u) by F2 "
          "are obtained from Figs. 22 and 23 as 6.8 x 10^-10, 7.6 x 10^-10, "
          "and 5 x 10^-10 cm3 s-1'", "checked", "Kr(3P1), a resonance "
          "state. Collisions at laser pressures keep the 5s manifold mostly "
          "in 3P2 (3P1 lies 0.12 eV higher), so the 3P2 values set k_H"),
    Value("k_H", 7.6e-10, 0.0, CHEN_PAYNE_1979, "p.157, as above", "checked",
          "Kr(1P1), a resonance state (see the 3P1 entry). Their Kr(3P1) "
          "and Kr(1P1) make KrF* 'with efficiency comparable to Kr(3P2)' "
          "(p.159), a relative statement"),
    Value("k_H", 8.1e-10, 0.0, "H. H. Nakano et al., via Ueda 1984 Table II "
          "and Kannari 1985 Table I", "Ueda 1984 Table II", "secondhand",
          "untraceable, so not used: the report both cite (" + SRI_1976 +
          ") gives 7.2e-10 from Velazco, Kolts & Setser in Table 2, "
          "reaction (14), p.46 and in the text on p.67, and '8.1' appears "
          "nowhere in it or in the midterm report MP-76-46 (checked in two "
          "OCR texts, 2026-09-28)"),
    # KrF* share of Kr* + F2
    Value("gamma_H", 0.90, 0.2, VKS_1976, "p.3476: 'the measured value of "
          "the KrF* formation rate constant (65 x 10^-11 cm3 molecule-1 "
          "sec-1)', over the total Kr(3P2) + F2 quenching rate constant 72 "
          "x 10^-11 (Table II)", "checked", "the KrF* rate constant is 0.88 "
          "of the KrCl* one (Table III, reproducible to +-10 %) times k_KrCl "
          "= 74e-11 (Sec. III C, whose absolute scale carries +-20 % from "
          "the ends of the spectral response and +-30 % from [Kr*]/[Ar*]). "
          "The +-0.20 is Table IV's. The ratio 0.7 in the 1976 ERDA report "
          "(COO-2807-2 p.17) came before this recalibration. LIGHT's value "
          "since 2026-09-28"),
    Value("gamma_H", 1.0, 0.2, VKS_1976, "Table IV: F2 with Kr(3P2), 1.0 "
          "('The uncertainties are estimated as +-0.20')", "checked",
          "the published value, an assignment: 'branching ratios for "
          "formation of XeF* and KrF* from F2 were assigned as unity' "
          "(p.3476), because 65e-11 equals the total quenching within "
          "experimental uncertainty, the XeF*/XeCl* rate ratio matches the "
          "total quenching ratio, and 'No other emissions were observed "
          "from the F2 reactions'. The Kr* + Cl2 reference there is KrCl* "
          "0.90 and Cl2* 0.10, which Velazco 1977 p.148 quotes. The group's "
          "ArF* check came out 0.5 'rather than unity' for Ar(3P2) + F2 "
          "(COO-2807-11 (1977), OSTI 6681231, p.6, sources/Setser_ERDA_"
          "COO-2807-11_1977_OSTI6681231.pdf)"),
    # share of the new KrF formed in C (nascent, Kr(3P1) + F2)
    Value("f_C_H", 0.38, 0.03, "Yu & Setser, J. Phys. Chem. 94, 2934 (1990)",
          "Sec. III.A after Fig. 2: 'the true kB/kC ratio will be somewhat "
          "larger than the apparent values, 1.3 +- 0.1, and 1.6 is a good "
          "estimate'", "checked", "the nascent split at a few Torr; C->B "
          "transfer is fast in high v and slow in low v (same page), so "
          "the relaxed split at 1 atm may be smaller. Setser's final ERDA "
          "report (COO-2807-19 (1979), OSTI 5572256, pp.4-9, sources/"
          "Setser_ERDA_COO-2807-19_1979_OSTI5572256.pdf) shows B/C 0.2 "
          "near 1 Torr and the same B/C from 3P1 and 3P2"),
    # Kr* + 2Kr -> Kr2* + Kr [cm6/s]
    Value("k3", 1.7e-32, 0.2e-32, UEDA_1984, "Table V, Kr*(3P2) absorption",
          "checked"),
    Value("k3", 1.79e-32, 0.0, "Salamero et al., J. Phys. B 12, 419 "
          "(1979)", "abstract: 1/tau3 = 18.6 p^2 + 3.9e3 p s^-1 (p in Torr) "
          "-> 18.6/(3.22e16)^2", "abstract", "formation of the radiative "
          "state from an intermediate state, a few hundred Torr"),
    Value("k3", 2.59e-32, 0.0, "Tracy & Oskam, J. Chem. Phys. 65, 1666 "
          "(1976)", "abstract: 'k2 = 2.59e-32 cm6/s' (5s[3/2]2, afterglow)",
          "abstract"),
    Value("k3", 4.2e-32, 0.0, "Smith & Turner, Can. J. Phys. 41, 1949 "
          "(1963)", "abstract: 1/tau = 30/p + 75p + 44p^2 (p in mm Hg) -> "
          "44/(3.22e16)^2", "abstract", "0.2-1.8 Torr: the p^2 term is a "
          "small part of the decay"),
    Value("k3", 4.4e-32, 0.0, "Boucique & Mortier, J. Phys. D 3, 1905 "
          "(1970)", "abstract: '44e-33 cm6/s for Kr_m'", "abstract",
          "indirect: the rise of delayed VUV in Townsend discharges"),
    Value("k3", 5.5e-32, 0.0, "W. M. Hughes, in Electronic Transition "
          "Lasers II (MIT Press, 1977) p.81, via " + QH_1978, "Table I: "
          "'Kr* + 2Kr -> Kr2* + Kr, k2 = 5.5 x 10^-32 cm6 s-1' (their ref. "
          "5)", "secondhand", "the original is not read"),
    Value("k3_tpow", -0.64, 0.0, "Wieme & Lenaerts, J. Chem. Phys. 74, 483 "
          "(1981)", "abstract: Ea(Kr) = -0.020 eV, 200-400 K", "abstract",
          "as (Tg/300)^n between 300 and 451 K: n = ln(0.772)/ln(1.503)"),
    # the Kr2* excimer: radiative rate [1/s] and F2 quenching [cm3/s]
    Value("kr2_rad", 3.7e6, 0.3e6, "Morikawa et al., J. Chem. Phys. 91, "
          "1469 (1989)", "abstract: 270 +- 5 ns (slow component, low "
          "density)", "abstract", "Salamero 1979 250 +- 20 ns, Quigley & "
          "Hughes 1978 280 +- 30 ns, Audouard 1988 265 ns, Oka 1974 353 ns, "
          "Leichner & Ericson 1974 about 0.30 us at 1000 Torr; the fast "
          "0u+ component (3.4 ns) carries F/S = 0.21 of the light"),
    Value("k_kr2_f2", 2.1e-10, 0.2e-10, QH_1978, "p.650 and Fig. 1: the "
          "slope of the Kr2F* decay against F2 (photolytic, 'no charged "
          "particles'), 'k3 = 2.1 +- 0.2 x 10^-10 cm3 sec-1'", "checked",
          "the total of both product routes (see ABSENT). Their harpoon "
          "estimate, the Kr2* + F2 crossing at 6.5 A scaled from Kr* + F2, "
          "gives 2.9e-10. REBUILD's 6e-10 is Kannari's 3e-10 to KrF* plus "
          "3e-10 to Kr2F*"),
    Value("k_kr2_f2", 5.0e-10, 1.25e-10, CHEN_PAYNE_1979, "p.157 (Fig. 23) "
          "and Fig. 25: Kr2*(1u) + F2 '5 x 10^-10 cm3 s-1', beta2 = 1.6e7 "
          "Torr-1 s-1 +- 25 %", "checked", "from the decay of the 1425 A "
          "continuum at 600 Torr Kr. The three measurements, 2.1 +- 0.2, 3 "
          "+- 0.9 and 5 +- 1.25 (e-10), disagree beyond their errors"),
    Value("kr2_rad", 1.0e6, 0.0, CHEN_PAYNE_1979, "Fig. 23 caption: 'The "
          "lifetime of Kr2*(1u) is determined as ~1 us'", "checked",
          "the zero-F2 intercept at 600 Torr Kr. The six other "
          "measurements give 250-353 ns"),
    Value("k_kr2_f2", 3.0e-10, 0.9e-10, SRI_1976, "Table 2, reaction (21), "
          "p.46-47: Kr2* + F2 -> Kr2F* + F, 3e-10, uncertainty 0.3, 'From "
          "decay of Kr2* and Kr2F* vs F2'", "checked", "the total removal "
          "rate with the product assigned, agreeing with Quigley & Hughes "
          "within its 30 %"),
]


def of(link):
    return [v for v in VALUES if v.link == link]


# The ionic channel. Its rate is measured only below 300 Torr (Jacob et al.
# 1978) and the Kr2+ channel's KrF* yield only relative to Kr+ (Rokni et
# al. 1977), both read first-hand 2026-09-29. The absolute KrF* yield, the
# B/C split and the rate near 1 atm have no measurement. The rest is model
# input or theory:
ION_NOTES = [
    Value("ion_scale", 0.63, 0.0, JACOB_1978, "pp.110-111: KrF* pulse "
          "shapes in e-beam Kr/F2 give alpha_i = 2 x 10^-7 at 60 Torr (Fig. "
          "1) and '3 x 10^-7 cm3/sec' at 120 Torr, 10^-7 extrapolated to "
          "zero pressure, and 'the effective two-body ion-ion recombination "
          "rate at 1 atm is estimated to be 10^-6 cm3/sec'", "checked",
          "Kr+ + F- (+M) below 300 Torr. Over LIGHT's Flannery form "
          "screened at their 1e14 cm-3 ions: 0.83 at 60 Torr, 0.63 at 120 "
          "Torr. The 1 atm value is an estimate"),
    Value("y_ion", 1.0, 0.0, ROKNI_1977, "p.2221: 'Experiments performed "
          "with different E-beam currents ... show that Kr2+ recombines "
          "with F- to form KrF* with a branching ratio near unity'",
          "checked", "relative to Kr+ + F-: the beam current changed x25 "
          "shifts Kr2+ against Kr+, and only relative intensities were "
          "used (the Ar case, Fig. 5, against branchings 1, 0.8 and 0.6). "
          "It rules out Kr2F* from Kr2+ + F-, not a loss common to both "
          "ions"),
    Value("ion_scale", 0.38, 0.0, "Avco Everett, 'One meter KrF laser "
          "system', May 1978, DTIC ADA069053 (32 MB, not kept: https://"
          "archive.org/download/DTIC_ADA069053/DTIC_ADA069053.pdf)",
          "Table 1 p.42: Kr+ + F- + (M) "
          "-> KrF* + (M): (1e-7 + 1e-6 P) cm3/s for P < 1 atm, 1.1e-6 above",
          "checked", "a model input built on a sub-300 Torr measurement "
          "(Jacob, Rokni, Mangano & Brochu, APL 32, 109 (1978)) and theory; "
          "1.1e-6 / LIGHT's Flannery 2.9e-6 = 0.38"),
    Value("y_ion", 0.1, 0.0, CHEN_PAYNE_1979, "p.159: 'the ratio of the "
          "intensity of ArF* or KrF* from ion-ion recombination to that from "
          "the slow energy transfer process is about 10 percent for F2 "
          "pressures higher than a few torr, although the population of "
          "argon ions is always more than the population of excited argon "
          "atoms'", "checked", "not a KrF* yield. At their low ion density "
          "F- is lost to F2 first (the slow part decays at 2e-10 cm3/s per "
          "F2, and 'the quenching of F- by F2 plays an important role in "
          "the ion kinetics', mechanism unknown). LIGHT has no such row. At "
          "Tahoma recombination outruns it: 0.7 % of F- loss, about 0.3 % "
          "of formation (our estimate from LIGHT's window state)"),
    Value("y_ion", 1.0, 0.0, "Avco Everett, 'One meter KrF laser system', "
          "Feb 1977, DTIC ADA047684, sources/Avco_OneMeterKrF_Feb1977_DTIC_"
          "ADA047684.pdf", "p.31: 'with sufficient Kr the "
          "channeling of energy from e-beam ionization to KrF* proceeds "
          "with essentially unit branching for pressures up to 4 atm'",
          "checked", "Ar-rich mixtures, inferred with a kinetic model from "
          "fluorescence; not a measured recombination yield"),
]

# Analogs (2026-09-28, step 4 of the owner's plan, LIGHT/docs/krf_data/
# FORMATION_LINKS_SEARCH.md). Measurements of neighbouring processes, not of
# the links. They set the UNMEASURED brackets and never a link.
ANALOGS = [
    Value("gamma_H", 0.5, 0.0, "Setser, ERDA COO-2807-11 (1977), OSTI "
          "6681231", "p.6: Ar(3P2) + F2 gives ArF* 0.5 'rather than unity'",
          "checked", "ArF, the same group and method"),
    Value("gamma_H", 0.6, 0.0, "Nakano et al. (SRI), via Brau in Rhodes (ed.), "
          "Excimer Lasers (1984), Table 4.4 p.98", "Ar* + F2 -> ArF*, <= 0.6",
          "secondhand", "an upper bound for ArF; SRI's ArF fluorescence "
          "yields were 'surprisingly small (30% being the greatest "
          "observed)' (" + SRI_1976 + ", p.40)"),
    Value("gamma_H", 1.0, 0.2, VELAZCO_1977, "Table 4.4 p.164: Xe(3P2) + "
          "F2 -> XeF* 1.0, relative to Xe* + Cl2 -> XeCl* assigned 1.0 "
          "('uncertainties are +-0.20')", "checked", "XeF, the same method. "
          "XeCl* was assigned 1.0 because KrCl* came out 0.90 and the XeCl* "
          "rate is the larger (p.148)"),
    Value("gamma_H", 0.31, 0.05, "Cooper, Mezyk & Young, J. Phys. Chem. 100, "
          "16195 (1996)", "abstract: Kr* + SF6 -> KrF*, G = 0.70 +- 0.10 per "
          "100 eV", "abstract", "SF6, not F2. Per Kr*: 0.70 over the "
          "0.55 x 100/24.2 = 2.27 Kr* per 100 eV of Saito's N_ex/N_i and W "
          "(our conversion)"),
    Value("y_ion", 0.64, 0.0, "Tsuji et al., J. Chem. Phys. 97, 1079 (1992)",
          "abstract: KrCl(B,C) : Kr(5s 3P1) : Cl(4s) = 0.64 : 0.23 : 0.13 "
          "from Kr+(2P3/2) + Cl- (+ He)", "abstract", "Cl-, not F-, in a "
          "flowing afterglow"),
    Value("y_ion", 0.41, 0.03, "Cooper, Mezyk & Young 1996", "abstract: Kr2+ "
          "+ SF6- + Kr -> KrF*, G = 1.70 +- 0.10 per 100 eV", "abstract",
          "SF6-, not F-. Per ion pair: 1.70 over 100/24.2 (our conversion)"),
    Value("f_C_ion", 0.38, 0.0, "Tsuji, Furusawa & Nishimura, J. Chem. Phys. "
          "92, 6502 (1990)", "abstract: Kr+(2P3/2) + SF6- -> KrF B:C = "
          "0.62:0.38", "abstract", "SF6-, not F-, flowing afterglow at 300 K"),
    Value("f_C_ion", 0.40, 0.04, "Tsuji et al., J. Chem. Phys. 94, 4291 "
          "(1991)", "abstract: Kr+(2P3/2) + Cl- + He -> KrCl B:C = 0.60 +- "
          "0.06 : 0.40 +- 0.04", "abstract", "Cl-, not F-"),
    Value("ion_scale", 1.0, 0.0, MEZYK_1991, "p.3155 and Fig. 4a: Kr2+ + "
          "SF6- (+ Kr) rises to 'a maximum value (2.9 x 10^15 M-1 s-1, ~4.8 "
          "x 10^-6 cm3 s-1) at approximately 300 Torr', about 3e-6 cm3/s "
          "near 665 Torr (read from the figure)", "checked", "SF6-, not F-, "
          "pulse radiolysis at 30-1600 Torr Kr. LIGHT's Flannery form gives "
          "2.9e-6 at 665 Torr, so the ratio is about 1. The method cancels "
          "the KrF* yield, so the paper gives no yield and no B/C split"),
]
# Links with no measurement at all, and why the numbers in use are not ones
ABSENT = [
    Value("k_kr2_f", 3.0e-10, 0.0, SRI_1976, "p.42: 'The interpretation is "
          "that this component results from the reaction Kr2* + F -> KrF* + "
          "Kr' (a slow KrF* component at very low F2), and no rate in Table 2",
          "checked", "Kannari's 3e-10 (Table I p.4311) cites this report as "
          "ref. 52, which tabulates no such rate. It is an estimate"),
    Value("kr2_krf", 0.5, 0.0, SRI_1976, "Table 2, reaction (21): the "
          "measured Kr2* + F2 rate is assigned to Kr2F* + F", "checked",
          "no measured KrF*/Kr2F* split exists"),
    Value("kr2_krf", 0.5, 0.0, QH_1978, "p.650: 'The relative probability "
          "for the occurrence of one route to Kr2F* or the other is not "
          "known and cannot be discerned in this experiment'", "checked",
          "their k3 is the total over KrF* + F + Kr and Kr2F* + F"),
    Value("k_kr2_f", 3.0e-10, 0.0, SRI_1976, "p.45: 'An important class of "
          "reactions has not yet been included; namely, reactions of Ar2* "
          "and Kr2* with F atoms'; p.65: 'if the reaction Ar2* + F -> ArF* "
          "+ Ar is sufficiently rapid (3-10 x 10^-10)'", "checked", "the "
          "model leaves Rg2* + F out and its one number is conditional and "
          "for argon, so Kannari's 3e-10 is its low end carried to Kr and "
          "no upper bound. Quigley & Hughes 1978, Jacob et al. 1978, Rokni "
          "et al. 1977 and Rokni et al., IEEE JQE 14, 464 (1978) give no "
          "Kr2* + F either"),
    Value("kr2_krf", 0.5, 0.0, CHEN_PAYNE_1979, "Fig. 25: f2'/f3 << 1, and "
          "p.158: 'we conclude that the dominant process to produce Kr2F* "
          "is through KrF* as an intermediate species'", "checked", "the "
          "direct Kr2F* channel of Kr2* + F2 is minor, which contradicts "
          "the half LIGHT sends there. The KrF* share itself is not given, "
          "and the whole split moves formation by 1 % at most"),
    Value("k_kr2_f", 3.0e-10, 0.0, JANCAITIS_1983, "Table II pp.94-95, KrF* "
          "Formation: the ion pairs, Kr** + F2, Kr* + F2, Kr2* + F2 and Kr + "
          "F2*, and no Kr2* + F", "checked", "the LLNL model leaves the "
          "reaction out, as SRI's does. LIGHT v0_2 carries it as a band "
          "from 0 (its default) to 1e-9"),
    Value("kr2_krf", 1.0, 0.0, JANCAITIS_1983, "Table II p.95: 'Kr2* + F2 "
          "-> KrF* + Kr + F', 3.0e-10 (his refs. 93, 95, Brau and Quigley & "
          "Hughes)", "checked", "all of Kr2* + F2 to KrF*, an assignment, "
          "at SRI's total rather than Quigley & Hughes' measured 2.1e-10"),
]


# The chain with every link that has a measurement set to it (LIGHT's own
# values since 2026-09-28, step 2). Choices:
#   k_H      the mean of the three first-hand values (7.2, 8.8, 9.4), as
#            LIGHT's basic network has used since 2026-09-28 (step 4). The
#            8.1e-10 credited to Nakano et al. is not in their report
#   gamma_H  0.90: the measured KrF* formation rate constant over the total
#            quenching (Velazco, Kolts & Setser 1976). LIGHT's since
#            2026-09-28
#   f_C_H    Yu & Setser's nascent split
#   k3       Ueda's, the one direct measurement at high pressure, with
#            Wieme & Lenaerts' temperature dependence
#   kr2_rad, k_kr2_f2  the measured slow Kr2* lifetime and F2 quenching
# Unmeasured links keep LIGHT's or REBUILD's values: kr2_krf, k_kr2_f and
# everything ionic. UNMEASURED brackets each one.
MEASURED = dict(k_H=8.47e-10, gamma_H=0.90, f_C_H=0.38, k3=1.7e-32,
                k3_tpow=-0.64, kr2_rad=3.7e6, k_kr2_f2=2.1e-10)
UNMEASURED = [
    ("harpoon KrF* yield 0.70 (0.90 less the stated 0.20)",
     dict(gamma_H=0.70)),
    ("harpoon KrF* yield 1.0 (the published assignment)", dict(gamma_H=1.0)),
    ("harpoon KrF* yield 0.5 (the ArF* analog)", dict(gamma_H=0.5)),
    ("nascent C relaxes into B at 1 atm (C share 0)", dict(f_C_H=0.0)),
    ("ionic KrF* formed 38 % in C (the SF6- analog)", dict(f_C_ion=0.38)),
    ("ionic KrF* yield 0.64 (the Cl- analog)", dict(y_ion=0.64)),
    ("ion-ion rate as Avco's 1.1e-6 (x0.38)", dict(ion_scale=0.38)),
    ("no Kr2* + F -> KrF", dict(k_kr2_f=0.0)),
    ("Kr2* + F2 all to Kr2F*", dict(kr2_krf=0.0)),
    ("Kr2* + F2 all to KrF*", dict(kr2_krf=1.0)),
    ("k3 4.4e-32 (Boucique)", dict(k3=4.4e-32)),
]


# Laser performance at Tahoma (score only): ml/data.py _tahoma, and the
# report's own theoretical formation efficiency for the Kr mixture
TAHOMA_THEORY = Value("eta", 0.28, 0.0, "Spectra Technology 1989, OSTI "
                      "10136715, Eq. 3-5 ('from earlier analysis (3-3) under "
                      "ideal conditions')", "p.3-14", "checked",
                      "a model value, not a measurement")
