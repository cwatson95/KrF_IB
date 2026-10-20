"""Species table: formation energy, composition, charge, heat capacity (ArF).

LIGHT_ArF (2026-09-30): the Ar/F2 version of LIGHT/v0_2/species.py. One
rare gas, argon, and no buffer. The KrF original is LIGHT/v0_2/species.py.

Formation energy E_form [eV] is the energy above ground-state Ar + F2 at
rest. It is what makes the energy ledger close BY CONSTRUCTION: every
reaction's energy release is sum(E_form reactants) - sum(E_form products),
so no row can create or destroy energy silently (the REBUILD rule,
SPECK/REBUILD/MODULES/species.py, adopted unchanged).

The electron-impact thresholds (Ar* 11.548 eV, Ar+ 15.76 eV) are taken
from the vendored IST-Lisbon cross-section file itself, so the energy an
electron pays in a collision and the energy the product stores are the same
number. The molecular energies are the measured band positions plus the
calculated well depths below; which ones are measured and which are
calculated is stated per species. The first-hand checks are in
docs/arf_data/SOURCES_CHECKED.md.

cv is the constant-volume heat capacity per particle in units of k_B
(monatomic 3/2, diatomic 5/2 with vibration frozen, nonlinear triatomic 3).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    name: str
    e_form_eV: float
    atoms: tuple            # (("Ar", 1), ("F", 1))
    charge: int             # units of e
    cv: float               # k_B per particle
    cite: str


# the rare gas and the roles the model and the report look up by name
RG = "Ar"                   # the rare gas (the only buffer)
RG_STAR = "Ar*"             # the lumped 4s manifold
RG_ION = "Ar+"
EXCIMER = "ArF(B)"          # the lasing state
TRIMER = "Ar2F*"
MANIFOLD = ("ArF(B)",)      # the ArF excimer states the model carries

HC_EV_NM = 1239.84193       # eV nm, photon energy = HC_EV_NM / lambda
E_F = 0.80                  # eV, D0(F2)/2 (REBUILD species.py)
EA_F = 3.40                 # eV, as F- books it (3.4011895 eV, Blondel)
E_AR_STAR = 11.548          # IST-Lisbon Ar(3P2) threshold
IP_AR = 15.76               # IST-Lisbon ionization threshold
# ArF(B). Its well depth is not measured (Krauss & Mies, in Excimer Lasers
# 2nd ed., Rhodes ed., Springer 1984, Table 2.2 p.14, leave D_e(B) blank for
# ArF and give only the calculated R_e(B) = 2.4 A). Chen & Payne, IEEE JQE
# 15, 149 (1979) p.150 measured the B-X band peak at 1933 A. B lies at least
# that photon above Ar + F, where the X state is flat, so E_form(B) is set
# at that spectroscopic lower bound: 0.80 + 6.414 = 7.214 eV, D_e(B) = IP(Ar)
# - EA(F) - 6.414 = 5.95 eV at most. A 193-nm photon then leaves no X-state
# repulsion to heat; Jancaitis's model value 7.24 eV (UCRL-53465 Table I,
# p.85) would leave 0.03 eV.
ARF_BX_NM = 193.3
E_ARF_B = E_F + HC_EV_NM / ARF_BX_NM                    # 7.2141 eV
# Ar2F* (4^2 Gamma), bound by 0.65 +- 0.08 eV below ArF* + Ar, calculated
# (Huestis, Marowsky & Tittel, Excimer Lasers 2nd ed. 1984, ch. 6 Table 6.3
# p.185; Wadt & Hay give 0.6 eV, via Rokni, Jacob & Mangano 1977 p.2223)
D_AR2F = 0.65
# Ar2+, bound by D_e = 1.23 eV (photoionization, Ng, Trevor, Mahan & Lee,
# JCP 66, 446 (1977), via Krauss & Mies Table 2.1 p.9; Wadt's calculation
# 1.19 eV). Ar2* (the lowest triplet Rydberg dimer), D_e = 0.68 eV below
# Ar* + Ar, calculated (Saxon & Liu, JCP 64, 3291 (1976), via Krauss & Mies
# Table 2.1 p.9)
D_AR2_ION = 1.23
D_AR2_STAR = 0.68

SPECIES = {s.name: s for s in [
    Species("Ar", 0.0, (("Ar", 1),), 0, 1.5, "ground state"),
    Species("F2", 0.0, (("F", 2),), 0, 2.5, "ground state"),
    Species("F", E_F, (("F", 1),), 0, 1.5,
            "D0(F2)/2, D0 = 1.60 eV (REBUILD species.py, spectre "
            "photolysis adjudication 2026-07-15)"),
    Species("F-", E_F - EA_F, (("F", 1),), -1, 1.5,
            "E_form(F) - EA(F), EA = 3.40 eV (REBUILD species.py)"),
    Species("Ar*", E_AR_STAR, (("Ar", 1),), 0, 1.5,
            "the 4s manifold (3P2, 3P1, 3P0, 1P1) lumped at its lowest "
            "level, the IST-Lisbon Ar(3P2) threshold 11.548 eV (NIST "
            "4s[3/2]2 11.5484 eV); excitation of the upper three pays their "
            "own thresholds and the excess heats the gas"),
    Species("Ar+", IP_AR, (("Ar", 1),), +1, 1.5,
            "IST-Lisbon ionization threshold 15.76 eV (NIST 15.7596 eV)"),
    Species("Ar2+", IP_AR - D_AR2_ION, (("Ar", 2),), +1, 2.5,
            "IP(Ar) - D_e(Ar2+), D_e = 1.23 eV measured by photoionization "
            "(Ng, Trevor, Mahan & Lee, JCP 66, 446 (1977), via Krauss & "
            "Mies, Excimer Lasers 2nd ed. 1984, Table 2.1 p.9)"),
    Species("Ar2*", E_AR_STAR - D_AR2_STAR, (("Ar", 2),), 0, 2.5,
            "Ar2* (3Sigma_u+, the 126-nm excimer), 0.68 eV below Ar* + Ar, "
            "calculated (Saxon & Liu, JCP 64, 3291 (1976), via Krauss & "
            "Mies Table 2.1 p.9); no measured well depth is among the "
            "sources held"),
    Species("ArF(B)", E_ARF_B, (("Ar", 1), ("F", 1)), 0, 2.5,
            "0.80 (F) + the 1933-A B-X band peak (Chen & Payne, IEEE JQE "
            "15, 149 (1979) p.150), 7.214 eV: the spectroscopic lower bound "
            "with a flat X state (D_e(B) is not measured, Krauss & Mies "
            "Table 2.2 p.14). The C state is not carried: its calculated "
            "splitting from B is about 0.25 eV (Krauss & Mies p.15), and no "
            "C position, mixing rate or equilibrium is measured"),
    Species("Ar2F*", E_ARF_B - D_AR2F, (("Ar", 2), ("F", 1)), 0, 3.0,
            "Ar2F* (4^2 Gamma) bound by 0.65 +- 0.08 eV below ArF* + Ar, "
            "calculated (Huestis, Marowsky & Tittel, Excimer Lasers 2nd ed. "
            "1984, Table 6.3 p.185)"),
]}

ELECTRON = "e"


def e_form(name: str) -> float:
    """Formation energy [eV]; the electron's kinetic energy lives in U_e."""
    return 0.0 if name == ELECTRON else SPECIES[name].e_form_eV


def atoms(name: str) -> dict:
    return {} if name == ELECTRON else dict(SPECIES[name].atoms)


def charge(name: str) -> int:
    return -1 if name == ELECTRON else SPECIES[name].charge
