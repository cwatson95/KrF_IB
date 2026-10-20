"""Species table: formation energy, composition, charge, heat capacity.

Formation energy E_form [eV] is the energy above ground-state Kr at rest.
It is what makes the energy ledger close BY CONSTRUCTION: every reaction's
energy release is sum(E_form reactants) - sum(E_form products), so no row
can create or destroy energy silently (the REBUILD rule,
SPECK/REBUILD/MODULES/species.py, adopted unchanged).

LIGHT_Kr v0_1 (2026-10-13) carries LIGHT v0_5's krypton species with
their values and citations unchanged, and none of its fluorine ones (F2,
F, F-, KrF(B), KrF(C), Kr2F*). v0_2 (2026-10-14) takes the electron-impact
thresholds from its own cross-section file, Biagi's Magboltz 8.97 set
(xsec.py): Kr* 9.9152 eV (1s5), Kr** 11.3035 eV (2p10) and Kr+ 13.9996 eV,
which are the NIST levels to the digits printed, so the energy an electron
pays in a collision and the energy the product stores are still the same
number. v0_1 had Biagi-v7.1's 9.915, 11.30 and 13.996 (REBUILD rounds them
to 9.92 / 14.00). Kr2+ is IP(Kr) - D0(Kr2+), the energy at and above which
an excited atom can ionize associatively with a ground-state one (the
Hornbeck-Molnar rows, network.py).

cv is the constant-volume heat capacity per particle in units of k_B
(monatomic 3/2, diatomic 5/2 with vibration frozen).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    name: str
    e_form_eV: float
    atoms: tuple            # (("Kr", 2),)
    charge: int             # units of e
    cv: float               # k_B per particle
    cite: str


_WELL = ("EST well depth (REBUILD species.py); it moves energy between "
         "'stored' and 'heat', never the ledger closure")

IP_KR = 13.9996             # eV, Biagi 8.97's ionization threshold (Kr+ below)

SPECIES = {s.name: s for s in [
    Species("Kr", 0.0, (("Kr", 1),), 0, 1.5, "ground state"),
    Species("Kr*", 9.9152, (("Kr", 1),), 0, 1.5,
            "Biagi (Magboltz 8.97) 1s5 threshold 9.9152 eV = the 5s[3/2]2 "
            "level, the lowest of the 4p^5 5s manifold (NIST ASD 9.9152 "
            "eV; v0_1's Biagi-v7.1 group 9.915, REBUILD 9.92)"),
    Species("Kr+", IP_KR, (("Kr", 1),), +1, 1.5,
            "Biagi (Magboltz 8.97) ionization threshold 13.9996 eV (NIST "
            "ASD 13.99961; v0_1 13.996, REBUILD 14.00)"),
    # LIGHT v0_4 (2026-10-08): the 5p rung of the ladder. One lumped state
    # for the 4p^5 5p manifold (ten levels, 11.3035 to 12.2565 eV, g = 36)
    # and the levels above it that cascade through it. Its energy is the
    # lowest 5p level, 5p[1/2]_1 (2p10), so the electron's cost of that
    # level and the stored energy are one number, as for Kr* and Kr+; the
    # manifold's statistical mean lies about 0.4 eV higher, which the
    # ledger books as heat at the excitation of every higher level (each
    # level is its own row since v0_2) and never loses. REBUILD and
    # Kannari carry 11.6 eV.
    Species("Kr**", 11.3035, (("Kr", 1),), 0, 1.5,
            "Biagi (Magboltz 8.97) 2p10 threshold 11.3035 eV = the "
            "5p[1/2]_1 level, the lowest of the 4p^5 5p manifold (NIST ASD "
            "11.3035 eV; v0_1's Biagi-v7.1 'Kr(p)' group 11.30); REBUILD "
            "and Kannari 1985 use 11.6 eV"),
    # --- the excimer and the molecular ion ---------------------------------
    Species("Kr2*", 9.2, (("Kr", 2),), 0, 2.5,
            "REBUILD Kr2* (147 nm VUV excimer), 0.72 eV below Kr* + Kr. "
            "No measured or calculated Kr2* binding energy is among the "
            "sources held (Krauss & Mies, Excimer Lasers 2nd ed. 1984, Table "
            "2.1 p.9, leave the Kr Rydberg dimer blank), so v0_2 keeps it. "
            + _WELL),
    Species("Kr2+", IP_KR - 1.15, (("Kr", 2),), +1, 2.5,
            "IP(Kr) - D0(Kr2+), D0 = 1.15 eV, measured by photoionization "
            "(Ng, Trevor, Mahan & Lee, JCP 66, 446 (1977), via Krauss & Mies, "
            "Excimer Lasers 2nd ed. 1984, Table 2.1 p.9). REBUILD carries the "
            "same value as an estimated well"),
]}

ELECTRON = "e"


def e_form(name: str) -> float:
    """Formation energy [eV]; the electron's kinetic energy lives in U_e."""
    return 0.0 if name == ELECTRON else SPECIES[name].e_form_eV


def atoms(name: str) -> dict:
    return {} if name == ELECTRON else dict(SPECIES[name].atoms)


def charge(name: str) -> int:
    return -1 if name == ELECTRON else SPECIES[name].charge
