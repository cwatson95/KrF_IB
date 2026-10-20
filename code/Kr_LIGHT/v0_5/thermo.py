"""Detailed balance for the dimers' redissociation (LIGHT_Kr v0_4,
2026-10-15): the equilibrium constant of A + B <-> M from the measured
binding energies, with the dimer's vibrational and rotational partition
functions, so that a measured association rate gives its reverse,
k_rev = k_fwd / K_eq.

    K_eq = [M] / ([A][B]) = (g_M / (g_A g_B)) (h^2 / 2 pi mu k T)^(3/2)
           q_rot q_vib exp(D0 / kT)                                [cm^3]

D0 is the well from the dimer's v = 0 level to the asymptote A + B at rest
(the species' formation energies, species.py, so that the ledger and the
equilibrium agree), q_vib sums the harmonic levels below D0 from v = 0,
q_rot = kT / (sigma B) with sigma = 2 (identical nuclei) and B from the
internuclear distance. g_M / (g_A g_B) counts the bound molecular states
over the atomic pair's (sigma = 2 and the g/u pair of an excited and a
ground atom cancel, so this is the plain ratio of degeneracies).

The dimers' constants (MOLECULES):
* Kr2* (A 3Sigma_u+: 1u and 0u-, g = 3). D0 0.709 eV (LaRocque et al.
  1986, D_e = 5807.0 cm^-1, less half Delta G_1/2), with the band kr2_d0
  down to Killeen & Eden's 0.62 eV; the vibrational quantum Delta G_1/2 =
  172 cm^-1 measured at low v by Killeen & Eden (JCP 84, 6048 (1986),
  abstract), where LaRocque's omega_e = 244.6 cm^-1 is extrapolated from
  v' = 32 to 38.
* Kr2*(0u+) (B 1Sigma_u+, g = 1): D0 0.684 eV below Kr*(1s4) + Kr and
  omega_e = 219.50 cm^-1 (LaRocque et al., abstract).
* Kr2+ (I(1/2u), g = 2; Kr+ 2P3/2, g = 4): D0 9267.8 +- 2.8 cm^-1 =
  1.1491 eV (Signorell, Hollenstein & Merkt, JCP 114, 9840 (2001), the
  ledger's since v0_5; Ng et al. 1977's 1.15 eV before); its vibrational
  quantum is not in the sources in hand,
  and the excimer's 172 cm^-1 stands in: the 5s Rydberg state is bound by
  the same A 2Sigma_u+ ion core (Killeen & Eden).
* re 2.6 A for all three: Huber & Herzberg's Kr2+ X re, '(2.6)', an
  estimate from Mittman & Weise's elastic scattering of Kr+ by Kr (NIST
  Chemistry WebBook, SRD 69, Kr2+ diatomic constants), the core the
  Rydberg excimers share. K_eq goes as re^2.
"""

from __future__ import annotations

import math
from functools import lru_cache

KB = 1.380649e-23            # J/K
H = 6.62607015e-34           # J s
C_CM = 2.99792458e10         # cm/s
AMU = 1.66053906660e-27      # kg
EV = 1.602176634e-19         # J
M_KR = 83.798                # amu, natural krypton
CM_PER_EV = 8065.543937      # cm^-1 per eV

# (g_M, omega [cm^-1], re [A]) per dimer; D0 comes from the species
MOLECULES = {
    "Kr2*": (3, 172.0, 2.6),
    "Kr2*(0u+)": (1, 219.50, 2.6),
    "Kr2+": (2, 172.0, 2.6),
}


@lru_cache(maxsize=4096)
def _keq(mol: str, d0_eV: float, g_atoms: int, tg: float) -> float:
    g_m, we, re_a = MOLECULES[mol]
    mu = 0.5 * M_KR * AMU
    kt_cm = KB * tg / (H * C_CM)
    b_cm = H / (8.0 * math.pi ** 2 * C_CM * mu * (re_a * 1e-10) ** 2)
    q_rot = kt_cm / (2.0 * b_cm)
    d0_cm = d0_eV * CM_PER_EV
    q_vib, v = 0.0, 0
    while v * we < d0_cm:
        q_vib += math.exp(-v * we / kt_cm)
        v += 1
    lam3 = (H ** 2 / (2.0 * math.pi * mu * KB * tg)) ** 1.5 * 1e6    # cm^3
    return (g_m / g_atoms) * lam3 * q_rot * q_vib * math.exp(
        d0_eV * EV / (KB * tg))


def k_eq(mol: str, d0_eV: float, g_atoms: int, tg: float) -> float:
    """[M] / ([A][B]) [cm^3] for the dimer mol bound by d0_eV below the
    asymptote whose atoms' degeneracies multiply to g_atoms, at tg [K]."""
    return _keq(mol, round(float(d0_eV), 9), int(g_atoms),
                round(max(float(tg), 50.0), 6))
