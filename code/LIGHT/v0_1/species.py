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

# v0_1 (2026-09-28): KrF(B), KrF(C) and Kr2F* from measurement. KrF(B)
# lies D_e below the ion pair Kr+ + F- (Tellinghuisen, Hays, Hoffman &
# Tisone, JCP 65, 4473 (1976): D_e(B) = 5.54 +- 0.03 eV, quoted by Krauss &
# Mies in Excimer Lasers, 2nd ed. 1984, Table 2.2 p.14). With the ledger's
# own IP(Kr) and EA(F) that puts it 5.056 eV above Kr + F, and the 248-nm
# photon leaves 0.06 eV of X-state repulsion, where REBUILD's estimate
# books 0.35 eV. The value is the potential minimum; v = 0 lies about
# 0.02 eV higher, inside the +-0.03 eV.
IP_KR = 13.996              # eV, Biagi's ionization threshold (Kr+ below)
EA_F = 3.40                 # eV, as F- books it
DE_KRF_B = 5.54             # eV, D_e(B) relative to Kr+ + F-
E_KRF_B = 0.80 + IP_KR - EA_F - DE_KRF_B          # 5.856 eV above Kr + F2
GAP_BC_CM1 = 120.0          # C below B, Yu & Setser 1990 p.2941 consensus
CM1_EV = 1.23984193e-4      # eV per cm^-1
D_KR2F = 0.58               # Kr2F* below KrF* + Kr, eV (Huestis et al.)

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
    Species("KrF(B)", E_KRF_B, (("Kr", 1), ("F", 1)), 0, 2.5,
            "0.80 (F) + IP(Kr) - EA(F) - D_e(B), D_e(B) = 5.54 +- 0.03 eV "
            "(Tellinghuisen et al. 1976, via Krauss & Mies Table 2.2); v0 "
            "keeps REBUILD's estimated 6.15 eV"),
    Species("Kr2F*", E_KRF_B - D_KR2F, (("Kr", 2), ("F", 1)), 0, 3.0,
            "Kr2F* (4^2 Gamma) bound by 0.58 +- 0.08 eV below KrF* + Kr, "
            "calculated (Huestis, Marowsky & Tittel, in Excimer Lasers, 2nd "
            "ed., Rhodes ed., Springer 1984, ch. 6 Table 6.3 p.185; Wadt & "
            "Hay 1978 give 0.49 eV with spin-orbit, and the DIM value of "
            "SRI MP-78-07 is 0.66 eV). v0 keeps REBUILD's 4.10 eV, an "
            "estimated well that put 2.05 eV of heat into trimer formation; "
            "with the calculated binding the 420-nm emission leaves about "
            "1.5 eV of lower-state repulsion"),
    # --- species that only the optional groups create -----------------------
    Species("Kr2*", 9.2, (("Kr", 2),), 0, 2.5,
            "REBUILD Kr2* (147 nm VUV excimer); " + _WELL),
    Species("Kr2+", 13.996 - 1.15, (("Kr", 2),), +1, 2.5,
            "IP(Kr) - D0(Kr2+), D0 = 1.15 eV (REBUILD EST well)"),
    Species("F-", 0.80 - 3.40, (("F", 1),), -1, 1.5,
            "E_form(F) - EA(F), EA = 3.40 eV (REBUILD species.py)"),
    Species("KrF(C)", E_KRF_B - GAP_BC_CM1 * CM1_EV, (("Kr", 1), ("F", 1)),
            0, 2.5, "120 +- 50 cm^-1 below B, 'the consensus value from our "
            "work' (Yu & Setser, JPC 94, 2934 (1990) p.2941), the one "
            "equilibrium measurement; their 230/300 K ratio route alone "
            "gives 180 +- 30. v0 keeps REBUILD's 80 cm^-1 (Julienne & "
            "Krauss's low-pressure reanalysis, via Kannari 1985 p.4313)"),
]}

ELECTRON = "e"


def e_form(name: str) -> float:
    """Formation energy [eV]; the electron's kinetic energy lives in U_e."""
    return 0.0 if name == ELECTRON else SPECIES[name].e_form_eV


def atoms(name: str) -> dict:
    return {} if name == ELECTRON else dict(SPECIES[name].atoms)


def charge(name: str) -> int:
    return -1 if name == ELECTRON else SPECIES[name].charge
