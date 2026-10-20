"""Species table: formation energy, composition, charge, heat capacity.

Formation energy E_form [eV] is the energy above ground-state Kr + F2 at
rest. It is what makes the energy ledger close BY CONSTRUCTION: every
reaction's energy release is sum(E_form reactants) - sum(E_form products),
so no row can create or destroy energy silently (the REBUILD rule,
SPECK/REBUILD/MODULES/species.py, adopted unchanged).

Values are REBUILD 0.1.40's, with one deliberate difference: the
electron-impact thresholds (Kr* 9.915 eV, Kr+ 13.996 eV) are taken from
the vendored Biagi-v7.1 cross-section file itself, so the energy an
electron pays in a collision and the energy the product stores are the
same number (REBUILD rounds them to 9.92 / 14.00).

cv is the constant-volume heat capacity per particle in units of k_B
(monatomic 3/2, diatomic 5/2 with vibration frozen, nonlinear triatomic 3).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    name: str
    e_form_eV: float
    atoms: tuple            # (("Kr", 1), ("F", 1))
    charge: int             # units of e
    cv: float               # k_B per particle
    cite: str


_WELL = ("EST well depth (REBUILD species.py); it moves energy between "
         "'stored' and 'heat', never the ledger closure")

SPECIES = {s.name: s for s in [
    # --- the V0 set (Gary's six reactions + closures) ----------------------
    Species("Kr", 0.0, (("Kr", 1),), 0, 1.5, "ground state"),
    Species("F2", 0.0, (("F", 2),), 0, 2.5, "ground state"),
    Species("F", 0.80, (("F", 1),), 0, 1.5,
            "D0(F2)/2, D0 = 1.60 eV (REBUILD species.py, spectre "
            "photolysis adjudication 2026-07-15)"),
    Species("Kr*", 9.915, (("Kr", 1),), 0, 1.5,
            "Biagi-v7.1 'Kr(s)' group threshold 9.915 eV = the 5s[3/2]2 "
            "level (REBUILD 9.92)"),
    Species("Kr+", 13.996, (("Kr", 1),), +1, 1.5,
            "Biagi-v7.1 ionization threshold 13.996 eV (REBUILD 14.00)"),
    Species("KrF(B)", 6.15, (("Kr", 1), ("F", 1)), 0, 2.5,
            "REBUILD KrF*(B): 0.80 (F) + 5.00 (248 nm quantum) + 0.35 "
            "X-state repulsion; " + _WELL),
    Species("Kr2F*", 4.10, (("Kr", 2), ("F", 1)), 0, 3.0,
            "REBUILD Kr2F* (4^2 Gamma): 0.80 + ~3.0 (~415 nm band) + ~0.3; "
            + _WELL),
    # --- species that only the optional groups create -----------------------
    Species("Kr2*", 9.2, (("Kr", 2),), 0, 2.5,
            "REBUILD Kr2* (147 nm VUV excimer); " + _WELL),
    Species("Kr2+", 13.996 - 1.15, (("Kr", 2),), +1, 2.5,
            "IP(Kr) - D0(Kr2+), D0 = 1.15 eV (REBUILD EST well)"),
    Species("F-", 0.80 - 3.40, (("F", 1),), -1, 1.5,
            "E_form(F) - EA(F), EA = 3.40 eV (REBUILD species.py)"),
    Species("KrF(C)", 6.15 - 0.0099, (("Kr", 1), ("F", 1)), 0, 2.5,
            "80 cm^-1 below B (Kannari 1985 p.4314, REBUILD KrF*(C))"),
]}

ELECTRON = "e"


def e_form(name: str) -> float:
    """Formation energy [eV]; the electron's kinetic energy lives in U_e."""
    return 0.0 if name == ELECTRON else SPECIES[name].e_form_eV


def atoms(name: str) -> dict:
    return {} if name == ELECTRON else dict(SPECIES[name].atoms)


def charge(name: str) -> int:
    return -1 if name == ELECTRON else SPECIES[name].charge
