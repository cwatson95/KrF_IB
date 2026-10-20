"""Chen 2025 Table III — verbatim transcription (the CHEN tier's source data).

COPY PROVENANCE (REBUILD): byte-verbatim copy of
``projects/spectre/modules/spectre/chen2025_table3.py`` taken 2026-08-04;
only this provenance note, the module name, the flag list below (extended
2026-08-04 R36; 2026-08-25 R13/R14/R23/R24) and the appended CHEN_ERRATA
block changed. The transcription itself, its R19/R26 typography
adjudications and the reference list below are the 2026-07-14 spectre pass
and are preserved exactly — re-transcribing would discard an audited
artifact. TABLE_III stays AS PRINTED; corrections live ONLY in CHEN_ERRATA.

Source: Jiayi Chen, Liyang Zhang, Zhihang Zhao, Yangyang Fu, "Fluid modeling
of discharge-pumped ArF excimer lasers: Effects of F2 concentration on pulsed
discharge dynamics," Phys. Plasmas 32, 113502 (2025), doi:10.1063/5.0287524.
Table III, p. 113502-11 (docs/Recent_KrF_Papers/chen2025.pdf), transcribed
2026-07-14 from the page image.  Their system: Ne/Ar/F2/Xe (96.4/3.5/0.1/0.01)
at 3.4 atm, 26 species, 41 reactions; heavy-particle rates primarily from
Ma et al. [13]; f(sigma) rows are computed from LXCat cross sections by a
Boltzmann solver (their BOLSIG+; spectre uses the in-repo LoKI-MC tables).

UNITS ARE SI AS PRINTED: m3/s (two-body), m6/s (three-body), s-1 (radiative),
m2 (cross sections).  Te in eV.

Chen's reference numbers used below:
  [13] X. Ma, L. Bai, Y. Zhu, X. Jiang, Y. Wu, Plasma Sources Sci. Technol.
       33, 075012 (2024)   (their primary heavy-particle source)
  [22] www.lxcat.net, SIGLO database (electron cross sections)
  [30] H. Akashi, Y. Sakai, T. Tagashira, J. Phys. D: Appl. Phys. 28, 445
       (1995)   (self-sustained discharge ArF model)
  [38] J. Tennyson et al., Plasma Sources Sci. Technol. 31, 095020 (2022)
       (Quantemol QDB)
  [39] S. Rauf, M. J. Kushner, J. Appl. Phys. 85, 3460 (1999)  (Ne chemistry)
  [40] F. Kannari, M. Obara, T. Fujioka, J. Appl. Phys. 57, 4309 (1985)
       (in-repo: docs/Recent_KrF_Papers/kannari1985.pdf)
  [41] A. Demyanov, L. Feenstra, P. J. M. Peters, A. M. Napartovich,
       W. J. Witteman, Appl. Phys. B 72, 823 (2001)  (discharge ArF model)

TRANSCRIPTION FLAGS (typography issues in the printed table; the first two
are recorded in PLAN (46), the 2026-08-25 set in CHANGELOG [0.1.25];
resolutions are marked on the rows or carried in CHEN_ERRATA, originals
preserved here):
  * R19 rate prints as "3.7e-14 * Te^0.5".  Dissociative recombination
    physically falls with Te and the cited source [39] uses Te^-0.5-class
    scaling; spectre implements Te^-0.5 with a config switch to run
    as-printed (SENSITIVITY: both bracketed in tests).
  * R36 (ArF* + Ar + Ne -> Ar2F* + Ne) prints as 3.5e-45 m6/s =
    3.5e-33 cm6/s. REBUILD flag 2026-08-04 (first-hand reads): Lee 1989
    Table I [58]/[59] gives the ArF* three-body class at 5.0e-31 (2Ar) /
    1.2e-31 (Ar+Ne) cm6/s and Kannari 1985 Table I gives 4.0e-31 —
    BOTH ~x100 above the printed value. Kept verbatim here (this file is
    the as-printed record; the arf_chen mode is Chen-verbatim by
    definition) — a candidate term in the recorded x3-5 per-joule
    formation gap, flagged for the comparison work.
  * R26 product prints as "Ne+ + 2Ar", which is endothermic by ~5.8 eV
    (IP Ne 21.56 > IP Ar 15.76) and involves no species consistent with the
    reactants at 1.6e-44 m6/s; by exothermicity, magnitude (Ne third body
    ~6x less efficient than Ar, cf. R27), and their Fig. 2 arrows, the
    intended row is Ar+ + Ar + Ne -> Ar2+ + Ne.  spectre implements that.
  * R13/R14 print 2.08e-8 / 5.61e-9 "m3/s" — PHYSICALLY IMPOSSIBLE as
    printed: ~1e4 above the electron gas-kinetic ceiling (~1e-12 m3/s);
    these are the cm3/s numbers of the cited QDB fits [38] mislabeled SI.
    Flagged 2026-08-25 (first-hand ceiling check + the 1DModelKrF
    primary-source pass); correction /1e6 in CHEN_ERRATA.
  * R23 prints 2.5e-43 m6/s — byte-identical to their own R27 (Ar+
    dimerization, a genuinely 2.5e-31 cm6/s process) pasted onto the
    EXCIMER row; the Ar* + Ar + M formation class measures 1.1e-32 cm6/s
    (SRI/Kannari lineage, and Ma [13] — Chen's own cited source — prints
    it). Flagged 2026-08-25; correction x1/22.7 in CHEN_ERRATA.
  * R24 prints 4e-40 m6/s = 4e-28 cm6/s — ~4 orders above ANY measured
    three-body rate (it would annihilate Ne* in ~0.4 ps and sever the
    Penning chain of the 96% Ne bath). The primary value is
    triple-confirmed at 4.0e-34 cm6/s (Rauf & Kushner [39] — Chen's own
    cited source — 4.0; Phelps 1959 5; Lee 1989 4.1). Flagged 2026-08-25;
    correction /1e6 in CHEN_ERRATA.
  NOTE (0.1.25): none of the numeric errata can describe the model Chen
  RAN — with the printed R13/R14/R23/R24 the network cannot reproduce his
  own figures and (measured) cannot even integrate under his Sec. II A
  base-case circuit — so the errata lane, not the as-printed lane, is the
  faithful reproduction of his MODEL; the as-printed lane remains the
  faithful record of his PAPER (network.chen_rows(errata=False)).

Each row: (rid, reactants, products, kind, rate, ref, note)
  kind: 'bolsig' (f(sigma): rate from EEDF tables), 'te_fit'
        (A*Te**B*exp(C/Te), params (A, B, C)), 'k2' (constant, m3/s),
        'k3' (constant, m6/s), 'rad' (s-1), 'sigma' (m2, photon process).
"""

# fmt: off
TABLE_III = [
    # --- A. Electron collisions -------------------------------------------
    ("R1",  ("e", "Ne"),        ("Ne*", "e"),         "bolsig", None,                        "22", ""),
    ("R2",  ("e", "Ne"),        ("Ne**", "e"),        "bolsig", None,                        "22", ""),
    ("R3",  ("e", "Ne"),        ("Ne+", "e", "e"),    "bolsig", None,                        "22", ""),
    ("R4",  ("e", "Ne*"),       ("Ne+", "e", "e"),    "bolsig", None,                        "22", "stepwise ionization"),
    ("R5",  ("e", "Ar"),        ("Ar*", "e"),         "bolsig", None,                        "22", ""),
    ("R6",  ("e", "Ar"),        ("Ar**", "e"),        "bolsig", None,                        "22", ""),
    ("R7",  ("e", "Ar"),        ("Ar+", "e", "e"),    "bolsig", None,                        "22", ""),
    ("R8",  ("e", "Ar*"),       ("Ar+", "e", "e"),    "bolsig", None,                        "22", "stepwise ionization"),
    ("R9",  ("e", "Ar**"),      ("Ar+", "e", "e"),    "te_fit", (1.56e-13, 0.71, -2.63),     "22", ""),
    ("R10", ("e", "F2"),        ("F-", "F"),          "bolsig", None,                        "22", "dissociative attachment"),
    ("R11", ("e", "F2"),        ("F2(a3)", "e"),      "bolsig", None,                        "22", ""),
    ("R12", ("e", "F2"),        ("F2(A1)", "e"),      "bolsig", None,                        "22", ""),
    ("R13", ("e", "F2(a3)"),    ("F", "F", "e"),      "te_fit", (2.08e-8, -0.51, -4.69),     "38", ""),
    ("R14", ("e", "F2(A1)"),    ("F", "F", "e"),      "te_fit", (5.61e-9, -0.23, -4.94),     "38", ""),
    ("R15", ("e", "F2"),        ("F2+", "e", "e"),    "bolsig", None,                        "22", ""),
    ("R16", ("e", "Xe"),        ("Xe*", "e"),         "bolsig", None,                        "22", "trace Xe (off in spectre KrF mode)"),
    ("R17", ("e", "Xe"),        ("Xe+", "e", "e"),    "bolsig", None,                        "22", "trace Xe (off in spectre KrF mode)"),
    ("R18", ("e", "Ne2*"),      ("Ne2+", "e", "e"),   "te_fit", (9.75e-15, 0.71, -3.4),      "39", ""),
    ("R19", ("e", "Ne2+"),      ("Ne", "Ne*"),        "te_fit", (3.7e-14, -0.5, 0.0),        "39", "PRINTED Te^+0.5; implemented Te^-0.5 per [39] physics -- see module docstring"),
    # --- B. Penning ionization --------------------------------------------
    ("R20", ("Ne*", "Ne*"),     ("Ne+", "Ne", "e"),   "k2",     5e-16,                       "13", "metastable pooling"),
    ("R21", ("Ne*", "Ar"),      ("Ar+", "Ne", "e"),   "k2",     6.5e-17,                     "40", "Penning; Kannari 1985 (in-repo PDF)"),
    # --- B'. Recombination -------------------------------------------------
    ("R22", ("F2+", "F-"),      ("F", "F", "F"),      "k2",     4e-14,                       "13", "ion-ion neutralization"),
    # --- C. Three-body heavy-particle collisions ---------------------------
    ("R23", ("Ar*", "Ar", "Ne"), ("Ar2*", "Ne"),      "k3",     2.5e-43,                     "13", ""),
    ("R24", ("Ne*", "Ne", "Ne"), ("Ne2*", "Ne"),      "k3",     4e-40,                       "39", ""),
    # --- D. Charge transfer -------------------------------------------------
    ("R25", ("Ne+", "Ar"),      ("Ar+", "Ne"),        "k2",     1e-17,                       "13", "slow non-resonant atomic CT"),
    ("R26", ("Ar+", "Ar", "Ne"), ("Ar2+", "Ne"),      "k3",     1.6e-44,                     "13", "PRINTED product 'Ne+ + 2Ar' is endothermic; corrected to Ar2+ + Ne -- see module docstring"),
    ("R27", ("Ar+", "Ar", "Ar"), ("Ar2+", "Ar"),      "k3",     2.5e-43,                     "41", ""),
    ("R28", ("Ne+", "Ne", "Ne"), ("Ne2+", "Ne"),      "k3",     4.4e-44,                     "39", ""),
    # --- D'. ArF* formation -------------------------------------------------
    ("R29", ("Ar*", "F2"),      ("ArF*", "F"),        "k2",     4e-16,                       "30", "harpoon"),
    ("R30", ("Ar+", "F-"),      ("ArF*",),            "k2",     9.1e-13,                     "30", "ion-ion recombination"),
    ("R31", ("Ar2+", "F-"),     ("ArF*", "Ar"),       "k2",     9.26e-13,                    "30", ""),
    # --- E. ArF* quenching ---------------------------------------------------
    ("R32", ("ArF*", "e"),      ("Ar", "F", "e"),     "k2",     2e-13,                       "41", ""),
    ("R33", ("ArF*", "F2"),     ("Ar", "F", "F", "F"), "k2",    1.9e-15,                     "41", ""),
    ("R34", ("ArF*", "Ar"),     ("Ar", "Ar", "F"),    "k2",     9e-18,                      "41", ""),
    ("R35", ("ArF*", "Ne"),     ("Ar", "Ne", "F"),    "k2",     1.6e-18,                    "13", ""),
    ("R36", ("ArF*", "Ar", "Ne"), ("Ar2F*", "Ne"),    "k3",     3.5e-45,                    "13", ""),
    # --- F. UV radiation ------------------------------------------------------
    ("R37", ("ArF*",),          ("Ar", "F", "hv193"), "rad",    2.381e8,                    "30", "spontaneous emission"),
    ("R38", ("ArF*", "hv193"),  ("Ar", "F", "hv193", "hv193"), "sigma", 4e-20,              "30", "stimulated emission cross section"),
    ("R39", ("Ar**",),          ("Ar", "hv95"),       "rad",    1e9,                        "30", ""),
    # --- G. UV absorption (photoionization/photodetachment) -------------------
    ("R40", ("Ar**", "hv193"),  ("Ar+", "e"),         "sigma",  1e-22,                      "30", "photoionization of Ar**"),
    ("R41", ("F-", "hv193"),    ("F", "e"),           "sigma",  1e-21,                      "30", "photodetachment"),
]
# fmt: on

# ---------------------------------------------------------------------------
# CHEN_ERRATA (0.1.25) — the numeric transcription errors of the PRINTED
# table, adjudicated against primary sources (each is physically impossible
# or provably mis-pasted AS PRINTED; the flag list in the module docstring
# carries the forensics, CHANGELOG [0.1.25] the measured consequences).
# rid -> (corrected magnitude, printed magnitude, one-line basis).
# For te_fit rows the magnitude is the prefactor A (exponents unchanged);
# for k3 rows it is the constant [m6/s]. TABLE_III itself stays AS PRINTED;
# network.chen_rows applies these when errata=True (the default — the
# errata lane is the reproduction of Chen's MODEL, the as-printed lane the
# record of his PAPER). The R19 exponent and R26 product adjudications
# predate this block and remain in-row/unconditional (PLAN (46)).
# ---------------------------------------------------------------------------
CHEN_ERRATA = {
    "R13": (2.08e-14, 2.08e-8,
            "cm3/s misprinted as m3/s (x1e4 above the gas-kinetic "
            "ceiling as printed); QDB [38] class"),
    "R14": (5.61e-15, 5.61e-9,
            "cm3/s misprinted as m3/s (same class as R13)"),
    "R23": (1.1e-44, 2.5e-43,
            "their R27 value pasted onto the excimer row; Ar*+Ar+M "
            "formation = 1.1e-32 cm6/s (Ma [13]/SRI/Kannari lineage)"),
    "R24": (4.0e-46, 4e-40,
            "x1e6 slip; primary triple-confirmed 4.0e-34 cm6/s "
            "(Rauf & Kushner [39] 4.0, Phelps 1959 5, Lee 1989 4.1)"),
    "R36": (3.5e-43, 3.5e-45,
            "exponent typo vs his cited source; ArF* three-body class "
            "1.2-5.0e-31 cm6/s (Lee 1989 Table I, Kannari 1985, Ma [13]) "
            "— the 2026-08-04 REBUILD flag, now applied"),
}

# Discharge/laser parameters from the paper (Table II + Sec. II C), for the
# Chen-verbatim ArF scoring mode:
CHEN_PARAMS = {
    "mixture": {"Ne": 0.964, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
    "pressure_atm": 3.4,
    "gas_temperature_K": 300.0,
    "applied_voltage_V": -18e3,
    "gap_m": 12e-3,
    "discharge_area_m2": 90e-4,
    # laser gain module (their Eqs. 7-18)
    "cavity_length_m": 1.05,      # l_c = 105 cm
    "gain_length_m": 0.90,        # l_g = 90 cm
    "K_axial": 1e-5,              # axial spontaneous-emission fraction
    "R_OC": 0.08,                 # output coupler reflectivity
    "L_op_slope_kc": 80.0,        # L_op = k_c * C_F2 (0.08 at 0.1% F2)
    "photon_wavelength_nm": 193.0,
}


def te_fit(A: float, B: float, C: float, te_eV: float) -> float:
    """A * Te^B * exp(C/Te) [m3/s], the Table III fit form (C <= 0)."""
    if te_eV <= 0.0:
        return 0.0
    import math
    return A * te_eV ** B * math.exp(C / te_eV)
