"""Species table: formation energy, composition, charge, heat capacity.

Formation energy E_form [eV] is the energy above ground-state Kr at rest.
It is what makes the energy ledger close BY CONSTRUCTION: every reaction's
energy release is sum(E_form reactants) - sum(E_form products), so no row
can create or destroy energy silently (the REBUILD rule,
SPECK/REBUILD/MODULES/species.py, adopted unchanged).

LIGHT_Kr v0_1 (2026-10-13) carries LIGHT v0_5's krypton species with
their values and citations unchanged, and none of its fluorine ones (F2,
F, F-, KrF(B), KrF(C), Kr2F*). The electron-impact thresholds (Kr* 9.915
eV, Kr** 11.30 eV, Kr+ 13.996 eV) are taken from the vendored Biagi-v7.1
cross-section file itself, so the energy an electron pays in a collision
and the energy the product stores are the same number (REBUILD rounds them
to 9.92 / 14.00).

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

IP_KR = 13.996              # eV, Biagi's ionization threshold (Kr+ below)

SPECIES = {s.name: s for s in [
    Species("Kr", 0.0, (("Kr", 1),), 0, 1.5, "ground state"),
    Species("Kr*", 9.915, (("Kr", 1),), 0, 1.5,
            "Biagi-v7.1 'Kr(s)' group threshold 9.915 eV = the 5s[3/2]2 "
            "level (REBUILD 9.92)"),
    Species("Kr+", 13.996, (("Kr", 1),), +1, 1.5,
            "Biagi-v7.1 ionization threshold 13.996 eV (REBUILD 14.00)"),
    # LIGHT v0_4 (2026-10-08): the 5p rung of the ladder. One lumped state
    # for the 4p^5 5p manifold (ten levels, 11.304 to 12.257 eV, g = 36)
    # and the levels above it that cascade through it. Its energy is the
    # Biagi-v7.1 'Kr(p)' threshold, the lowest 5p level 5p[1/2]_1, so the
    # electron's cost in C1a and the stored energy are one number, as for
    # Kr* and Kr+; the manifold's statistical mean lies about 0.4 eV
    # higher, which the ledger books as heat at the cascade (C1b, C1c)
    # and never loses. REBUILD and Kannari carry 11.6 eV.
    Species("Kr**", 11.30, (("Kr", 1),), 0, 1.5,
            "Biagi-v7.1 'Kr(p)' group threshold 11.30 eV = the 5p[1/2]_1 "
            "level, the lowest of the 4p^5 5p manifold (NIST ASD: 11.304 "
            "eV); REBUILD and Kannari 1985 use 11.6 eV"),
    # --- the excimer and the molecular ion ---------------------------------
    Species("Kr2*", 9.2, (("Kr", 2),), 0, 2.5,
            "REBUILD Kr2* (147 nm VUV excimer), 0.72 eV below Kr* + Kr. "
            "No measured or calculated Kr2* binding energy is among the "
            "sources held (Krauss & Mies, Excimer Lasers 2nd ed. 1984, Table "
            "2.1 p.9, leave the Kr Rydberg dimer blank), so v0_2 keeps it. "
            + _WELL),
    Species("Kr2+", 13.996 - 1.15, (("Kr", 2),), +1, 2.5,
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
