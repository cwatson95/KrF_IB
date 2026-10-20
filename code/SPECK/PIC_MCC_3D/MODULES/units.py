"""SI constants and the plasma scales this package's gates are written against.

Self-contained on purpose: this package must be runnable if it is moved out of
the repo (the ``ARPAE/pic3d`` precedent).  Values are CODATA 2018, matching
``ARPAE/pic3d/units.py`` so numbers are comparable across the two packages.
"""

from __future__ import annotations

import numpy as np

# --- CODATA 2018 -----------------------------------------------------------
C_LIGHT = 2.99792458e8          # m/s (exact)
E_CHARGE = 1.602176634e-19      # C (exact)
EPS0 = 8.8541878128e-12         # F/m
MU0 = 1.25663706212e-6          # H/m
M_E = 9.1093837015e-31          # kg
M_P = 1.67262192369e-27         # kg
K_B = 1.380649e-23              # J/K (exact)
AMU = 1.66053906660e-27         # kg
EV = E_CHARGE                   # J per eV
ATM_PA = 101325.0               # Pa per atm (exact)
# Loschmidt constant: ATM_PA / (K_B * 273.15) = 2.6867801118e25 m^-3, which is
# also the CODATA 2018 value.  NOTE it differs from ARPAE/pic3d/units.py's
# 2.6867811e25 by 3.7e-7 relative -- a rounding slip in ARPAE's 7th digit
# (…7801… written as …7811…).  Physically irrelevant, but it is why the
# ARPAE oracle gate in tests/test_mcc3d.py builds both operators on ONE
# explicitly shared density instead of letting each package compute its own:
# the gate is there to test the collision OPERATOR, not an ideal-gas constant.
N_LOSCHMIDT_M3 = 2.686780111e25  # m^-3, ideal gas at 1 atm and 273.15 K


def gas_density_m3(pressure_atm: float, T_gas_K: float = 300.0) -> float:
    """Ideal-gas number density [m^-3] at ``pressure_atm`` and ``T_gas_K``."""
    return N_LOSCHMIDT_M3 * float(pressure_atm) * (273.15
                                                   / max(float(T_gas_K), 1.0))


def thermal_speed(Te_eV: float, mass: float = M_E) -> float:
    """sqrt(kT/m) [m/s] -- the PER-AXIS thermal speed (not the mean speed).

    This is the convention every loader and gate in this package uses: a
    Maxwellian is loaded as ``rng.normal(0, thermal_speed(...))`` per axis.
    """
    return float(np.sqrt(max(Te_eV, 0.0) * EV / mass))


def plasma_frequency(n_e_m3: float) -> float:
    """Electron plasma ANGULAR frequency omega_pe [rad/s]."""
    return float(np.sqrt(max(n_e_m3, 0.0) * E_CHARGE ** 2 / (EPS0 * M_E)))


def debye_length(n_e_m3: float, Te_eV: float) -> float:
    """Electron Debye length [m] = sqrt(eps0 kTe / (n e^2)) = v_th / omega_pe."""
    if n_e_m3 <= 0.0:
        return np.inf
    return float(np.sqrt(EPS0 * max(Te_eV, 0.0) * EV / (n_e_m3 * E_CHARGE ** 2)))


def courant_dt_em(h: float) -> float:
    """The 3-D Yee/FDTD stability limit dt <= h / (c sqrt(3)) [s].

    Not used by any solver here -- this package is electrostatic.  It exists
    so the ES/EM timestep ratio quoted in ``../PIC_MCC_3D.txt`` section 2.1
    can be recomputed from code rather than trusted from a document.
    """
    return float(h / (C_LIGHT * np.sqrt(3.0)))


def es_over_em_dt_ratio(Te_eV: float, safety: float = 0.2) -> float:
    """dt_ES / dt_EM at dx = lambda_D -- independent of density.

    dt_ES = safety / omega_pe, dt_EM = lambda_D / (c sqrt 3), and
    lambda_D = v_th / omega_pe, so the ratio is safety * c * sqrt(3) / v_th.
    Documents why a discharge does not belong on the FDTD path.
    """
    return float(safety * C_LIGHT * np.sqrt(3.0) / thermal_speed(Te_eV))
