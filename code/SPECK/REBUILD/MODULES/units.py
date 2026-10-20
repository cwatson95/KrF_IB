"""SI constants and the validity scales the gates are written against.

SI THROUGHOUT IS A DESIGN RULE, NOT A PREFERENCE. phantom11's internal
cm / 1e16 cm^-3 / ns unit system is itself a documented bug source
(ORESTES_CLASS_REBUILD.txt §4: "phantom11's cm / 1e16 / ns internal units are
themselves a bug source"); spectre went SI and this package inherits that.
Densities m^-3, rates m^3/s and m^6/s, energies J internally with eV at the
interfaces (always suffixed ``_eV``), lengths m, times s.

Constant values are CODATA 2018 exact/recommended (the 2019 SI redefinition
made e, h, k_B exact). The Loschmidt constant is the derived ideal-gas value
2.686780111...e25 m^-3 — NOT ARPAE's rounded 2.6867811e25, whose 7th-digit
slip desynchronised a bit-exact oracle once (PIC_MCC_3D/CHANGELOG.md [0.3.0]).
"""

from __future__ import annotations

import math

# --- exact SI (2019 redefinition) ------------------------------------------
E_CHARGE = 1.602176634e-19      # C — exact
EV = E_CHARGE                   # J per eV — exact
K_B = 1.380649e-23              # J/K — exact
H_PLANCK = 6.62607015e-34       # J s — exact
C_LIGHT = 2.99792458e8          # m/s — exact
N_AVOGADRO = 6.02214076e23      # 1/mol — exact

# --- CODATA 2018 recommended ------------------------------------------------
M_E = 9.1093837015e-31          # kg
AMU = 1.66053906660e-27         # kg
EPS0 = 8.8541878128e-12         # F/m
R_E = 2.8179403262e-15          # classical electron radius, m
MEC2_EV = 510998.95000          # electron rest energy, eV

K_B_EV_PER_K = K_B / EV         # 8.617333...e-5 eV/K
ATM_PA = 101325.0               # Pa — exact by definition

# Loschmidt (273.15 K, 1 atm), derived: p/(k_B T). Kept as the derived
# expression so it cannot drift from the constants above.
LOSCHMIDT_M3 = ATM_PA / (K_B * 273.15)
# One amagat: the Loschmidt density (used by Flannery/Lee pressure-turnover
# ion-ion recombination forms, which are parameterised in amagat).
AMAGAT_M3 = LOSCHMIDT_M3


def gas_density_m3(pressure_atm: float, T_K: float) -> float:
    """Ideal-gas number density [m^-3]."""
    if pressure_atm < 0.0 or T_K <= 0.0:
        raise ValueError(f"unphysical gas state: {pressure_atm} atm, {T_K} K")
    return pressure_atm * ATM_PA / (K_B * T_K)


def photon_energy_J(wavelength_nm: float) -> float:
    """E = h c / lambda [J]."""
    return H_PLANCK * C_LIGHT / (wavelength_nm * 1e-9)


# --- plasma scales (validity/regime reporting, not stepping) ----------------
def plasma_frequency(n_e_m3: float) -> float:
    """omega_pe [rad/s]."""
    return math.sqrt(max(n_e_m3, 0.0) * E_CHARGE ** 2 / (EPS0 * M_E))


def debye_length(n_e_m3: float, Te_eV: float) -> float:
    """lambda_D [m]."""
    if n_e_m3 <= 0.0:
        return math.inf
    return math.sqrt(EPS0 * max(Te_eV, 0.0) * EV / (n_e_m3 * E_CHARGE ** 2))


def critical_density_m3(f_hz: float) -> float:
    """Microwave critical density n_c = eps0 m_e (2 pi f)^2 / e^2 [m^-3]."""
    w = 2.0 * math.pi * f_hz
    return EPS0 * M_E * w * w / (E_CHARGE * E_CHARGE)
