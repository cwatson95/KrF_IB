"""Gates for the neutral background.

Small module, but every collision rate in the package is proportional to what
it returns, so the arithmetic is worth pinning.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.gas import GAS_MASSES_AMU, GasMixture, pure
from PIC_MCC_3D.MODULES.units import AMU, K_B, N_LOSCHMIDT_M3


def test_density_is_the_ideal_gas_law():
    g = pure("Ar", 1.0, 273.15)
    assert abs(g.total_density_m3 / N_LOSCHMIDT_M3 - 1.0) < 1e-12
    # and it scales as p/T
    assert abs(pure("Ar", 2.0, 273.15).total_density_m3
               / g.total_density_m3 - 2.0) < 1e-12
    assert abs(pure("Ar", 1.0, 546.30).total_density_m3
               / g.total_density_m3 - 0.5) < 1e-9


def test_pressure_round_trips():
    """n k T must give back the pressure it was built from."""
    g = pure("Ar", 1.0, 300.0)
    assert abs(g.pressure_Pa / 101325.0 - 1.0) < 2e-3   # Loschmidt vs exact atm


def test_component_densities_follow_the_fractions():
    mix = GasMixture({"Ar": 0.9495, "Kr": 0.05, "F2": 0.0005}, 1.0)
    n = mix.total_density_m3
    assert abs(mix.density_of("Kr") / n - 0.05) < 1e-12
    assert abs(mix.density_of("F2") / n - 0.0005) < 1e-12
    assert mix.density_of("He") == 0.0
    assert abs(sum(mix.densities.values()) - n) < 1e-6 * n


def test_fractions_are_not_renormalised():
    """A 2 % admixture on top of a full carrier stays what the user wrote."""
    mix = GasMixture({"Ar": 1.0, "F2": 0.02}, 1.0)
    assert mix.fractions["Ar"] == 1.0
    assert mix.fractions["F2"] == 0.02


def test_masses_and_mean_mass():
    mix = GasMixture({"Ar": 0.5, "Kr": 0.5}, 1.0)
    assert abs(mix.mass_of("Ar") / AMU - GAS_MASSES_AMU["Ar"]) < 1e-12
    want = 0.5 * (GAS_MASSES_AMU["Ar"] + GAS_MASSES_AMU["Kr"]) * AMU
    assert abs(mix.mean_mass - want) / want < 1e-12


def test_thermal_speed_matches_kT_over_M():
    g = pure("Ar", 1.0, 300.0)
    want = np.sqrt(K_B * 300.0 / (GAS_MASSES_AMU["Ar"] * AMU))
    assert abs(g.thermal_speed("Ar") - want) / want < 1e-12


def test_mean_free_path():
    g = pure("Ar", 1.0, 300.0)
    sigma = 1e-20
    assert abs(g.mean_free_path(sigma) - 1.0 / (g.total_density_m3 * sigma)) < 1e-30
    assert np.isinf(g.mean_free_path(0.0))


def test_depletion_gauges_the_static_gas_assumption():
    g = pure("Ar", 1.0, 300.0)
    assert g.depletion(0.0) == 0.0
    assert abs(g.depletion(0.01 * g.total_density_m3) - 0.01) < 1e-12


def test_scaled_to_keeps_composition():
    mix = GasMixture({"Ar": 0.95, "Kr": 0.05}, 1.0, 300.0)
    hot = mix.scaled_to(T_gas_K=600.0)
    assert hot.fractions == mix.fractions
    assert abs(hot.total_density_m3 / mix.total_density_m3 - 0.5) < 1e-9


@pytest.mark.parametrize("bad,msg", [
    ({}, "every fraction is zero"),
    ({"Ar": 0.0}, "every fraction is zero"),
    ({"Ar": -0.1}, "negative"),
    ({"Ar": 95.0}, "not a mixture"),
    ({"Xe": 1.0}, "unknown gas"),
])
def test_constructor_rejects_nonsense(bad, msg):
    with pytest.raises((ValueError, KeyError), match=msg):
        GasMixture(bad, 1.0)


def test_constructor_rejects_bad_state():
    with pytest.raises(ValueError, match="pressure_atm"):
        GasMixture({"Ar": 1.0}, 0.0)
    with pytest.raises(ValueError, match="T_gas_K"):
        GasMixture({"Ar": 1.0}, 1.0, -5.0)
