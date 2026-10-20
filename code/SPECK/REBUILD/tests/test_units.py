"""units.py gates: exact constants, derived scales."""

import math

from REBUILD.MODULES import units as U


def test_exact_si_constants():
    assert U.E_CHARGE == 1.602176634e-19        # exact by SI definition
    assert U.K_B == 1.380649e-23
    assert U.H_PLANCK == 6.62607015e-34
    assert U.C_LIGHT == 2.99792458e8


def test_loschmidt_is_derived_not_rounded():
    """The ARPAE units.py Loschmidt slip (7th digit) is the recorded reason
    this constant must be DERIVED (PIC_MCC_3D CHANGELOG [0.3.0])."""
    assert abs(U.LOSCHMIDT_M3 / 2.6867801118e25 - 1.0) < 1e-9
    assert U.gas_density_m3(1.0, 273.15) == U.LOSCHMIDT_M3


def test_photon_energy_248nm():
    ev = U.photon_energy_J(248.0) / U.EV
    assert abs(ev - 4.9997) < 1e-3


def test_plasma_scales():
    # n_e = 1e21 m^-3 (= 1e15 cm^-3), Te = 1 eV: the PIC_MCC_3D.txt §3
    # table row — lambda_D = 0.235 um, f_pe = 283.9 GHz
    assert abs(U.plasma_frequency(1e21) / (2 * 3.14159265 * 283.9e9) - 1) \
        < 0.01
    assert abs(U.debye_length(1e21, 1.0) / 2.35e-7 - 1) < 0.01
    # 95 GHz critical density = 1.12e14 cm^-3 (ARPAE: peak 7.3e16 = 651x
    # critical -> 7.3e16/651 = 1.12e14)
    assert abs(U.critical_density_m3(95e9) / 1.12e20 - 1) < 0.02


def test_unphysical_gas_state_raises():
    import pytest
    with pytest.raises(ValueError):
        U.gas_density_m3(1.0, -3.0)
    assert U.gas_density_m3(0.0, 300.0) == 0.0
    assert math.isinf(U.debye_length(0.0, 1.0))
