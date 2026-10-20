"""Pure rare gases: measurements of the half of the KrF chain that needs no
F2 (ml/pure.py). Each entry names its source, the place in it, and how it
was read (checked = re-read first-hand, abstract, secondhand), as in
ml/formation_data.py. Tahoma's pure-Kr absorption only scores.
"""
from __future__ import annotations

from .formation_data import Value

IAEA_799 = ("IAEA-TECDOC-799 (1995), ch. 8, OSTI ETDE 98221 "
            "(https://www.osti.gov/etdeweb/servlets/purl/98221)")
SAITO_2003 = ("Saito et al., IEEE Trans. Nucl. Sci. (2003), doi "
              "10.1109/TNS.2003.820615: 5.49 MeV alphas, gas at 1-10 atm")
MOROZOV_2008 = ("Morozov, Heindl, Krucken, Ulrich & Wieser, J. Appl. Phys. "
                "103, 103301 (2008): continuous 8-15 keV electron beams, "
                "75-1400 hPa, LIGHT/docs/krf_data/sources/morozov2008.pdf")

# What Morozov's efficiency counts (read in full 2026-09-28, the owner's
# copy in LIGHT/docs/krf_data/sources/morozov2008.pdf): all the VUV power
# emitted into 4 pi over the power deposited in the gas. The integration
# is not restricted to one continuum, so the 22 % shortfall below the
# one-photon ideal is not a definition effect.
MOROZOV_COUNTS = (
    "p.103301-3: 'Energy conversion efficiency is defined here as the ratio "
    "of the power emitted by the target gas into 4 pi solid angle in the VUV "
    "(or XUV) wavelength range and the power deposited in the gas by the "
    "electron beam'. The VUV photocurrent is the difference with and "
    "without a glass filter, converted with the excimer spectrum. Kr above "
    "500 hPa: apparent 0.35 +- 0.02 (12 keV), divided by the deposited "
    "fraction 0.67 from GEANT4 (+-10 %, backscatter included) and corrected "
    "for reflection on the foil (0.26), gives 0.42 +- 0.05 (systematic: "
    "emission 3 %, current 3 %, reflection 5 %, GEANT4 10 %)")

# energy per ion pair W [eV]
W = [
    Value("W_Kr", 24.2, 0.0, "Mano et al., arXiv:2102.00032 (2021), citing "
          "Krajcar Bronic, Hoshasen 24, 101 (1998)", "Sec. 3: 'The w-value "
          "in Kr is 24.2 eV'", "secondhand"),
    Value("W_Kr", 24.1, 0.0, IAEA_799, "p.610: 'values for gases (26.3 eV "
          "for Ar, 24.1 eV for Kr, and 21.9 eV for Xe)'", "checked",
          "the particle is not named there"),
    Value("W_Ar", 26.4, 0.0, IAEA_799, "Fig. 8.4: 'The ICRU [1] reference "
          "value of 26.4 eV' (electrons in argon)", "checked"),
    Value("W_Ne", 37.0, 0.0, "PDG 2010 review 'Particle detectors at "
          "accelerators', Table 28.5", "minimum-ionizing particles, called "
          "approximate", "secondhand"),
]
# excitations per ion pair, measured in the gas
NEX = [
    Value("nex_Ar", 0.52, 0.0, SAITO_2003, "abstract: 'found to be 0.52, "
          "0.55, and 0.60 in argon, krypton and xenon'", "abstract",
          "photons and electrons counted together; 'one scintillation "
          "photon is emitted from one recombination process'"),
    Value("nex_Kr", 0.55, 0.0, SAITO_2003, "abstract", "abstract"),
]
# energy per VUV photon if every pair recombines, W_s [eV]
WS = [
    Value("ws_Ar", 17.5, 0.0, SAITO_2003, "abstract: W_s 17.5, 15.4, 13.0 "
          "eV (Ar, Kr, Xe)", "abstract"),
    Value("ws_Kr", 15.4, 0.0, SAITO_2003, "abstract", "abstract"),
]
# fraction of deposited e-beam energy emitted as excimer VUV (> 500 hPa)
VUV = [
    Value("vuv_Ne", 0.31, 0.06, MOROZOV_2008, "p.103301-7 and Fig. 10",
          "checked"),
    Value("vuv_Ar", 0.33, 0.04, MOROZOV_2008, "p.103301-7 and Fig. 10",
          "checked"),
    Value("vuv_Kr", 0.42, 0.05, MOROZOV_2008, "p.103301-7 and Fig. 10",
          "checked",
          "low pump: a continuous beam of a few keV, so no Penning or "
          "excimer-excimer losses. " + MOROZOV_COUNTS),
]
# Kr ion chemistry
KR_DR = Value("dr_Kr2p", 1.6e-6, 0.16e-6, "Shiu & Biondi, Phys. Rev. A 16, "
              "1817 (1977)", "abstract: alpha(Kr2+) = 1.6e-6 (300/Te)^0.55 "
              "cm3/s +-10 %, 300-8400 K; the 5p states are the main "
              "products", "abstract")
KR_DR_POW = 0.55
KR_CONV = Value("conv_Kr", 342.0, 15.0, "Smith, Dean & Plumb, J. Phys. B 5, "
                "2134 (1972)", "abstract: beta = (342 +- 15) Torr^-2 s^-1 "
                "at 300 K, the same at 180 K, a 'marked decrease' by 510 K",
                "abstract", "Torr at 300 K gives 3.3e-31 cm6/s, Torr reduced "
                "to 273 K gives 2.7e-31; the abstract does not say which")
# VUV photon energies of the second continua [eV] (approximate peaks)
HV_VUV = {"Kr": 8.43, "Ar": 9.84}
# Tahoma's pure-Kr run (device: scores only)
TAHOMA_PURE_KR = Value("absorb_257", 0.67, 0.06, "Spectra Technology 1989, "
                       "OSTI 10136715", "Table 3-2: pure Kr, 257.3 nm, "
                       "-0.67 %/cm (beam stop in), -0.72 (stop out)",
                       "checked", "needs absorption cross sections of Kr2+ "
                       "and Kr2*, which are theory or assumption (UCRL-53465 "
                       "Table III)")

ALL = W + NEX + WS + VUV + [KR_DR, KR_CONV, TAHOMA_PURE_KR]


def value(name, pick=0):
    return [v for v in ALL if v.link == name][pick].value
