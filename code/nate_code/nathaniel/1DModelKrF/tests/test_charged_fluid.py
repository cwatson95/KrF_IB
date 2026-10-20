"""Charged-species fluid container: densities and electron energy."""

from __future__ import annotations

import unittest

import numpy as np

from Fluid1D.charged_fluid import ChargedFluid
from Fluid1D.components import Preionization
from particles import excimer_system_registry
from particles.species import KR


def _charged_species():
    reg = excimer_system_registry()
    return tuple(s for s in reg if s.is_charged)


class TestChargedFluid(unittest.TestCase):
    def setUp(self) -> None:
        self.cf = ChargedFluid(_charged_species(), 8, 0.1, area_m2=1.0)

    def test_rejects_neutral_species(self) -> None:
        with self.assertRaises(ValueError):
            ChargedFluid((KR,), 8, 0.1)

    def test_seed_and_total(self) -> None:
        self.cf.seed_uniform("e-", 1.0e18)
        np.testing.assert_allclose(self.cf.density_of("e-"), 1.0e18)
        self.assertAlmostEqual(self.cf.total_number("e-"), 1.0e18 * self.cf.cell_volume * 8)

    def test_seed_profile_sets_column(self) -> None:
        prof = np.linspace(1.0e13, 1.0e14, self.cf.n_cells)
        self.cf.seed_profile("e-", prof)
        np.testing.assert_allclose(self.cf.density_of("e-"), prof)

    def test_seed_profile_shape_checked(self) -> None:
        with self.assertRaises(ValueError):
            self.cf.seed_profile("e-", np.zeros(self.cf.n_cells + 1))


    def test_electron_energy_and_temperature(self) -> None:
        self.cf.seed_electron_energy(6.0)  # <eps> = 6 eV
        np.testing.assert_allclose(self.cf.electron_mean_energy(), 6.0)
        np.testing.assert_allclose(self.cf.electron_temperature_ev(), 4.0)  # Te = 2/3 <eps>

    def test_net_charge_density_quasineutral(self) -> None:
        self.cf.seed_uniform("e-", 1.0e18)
        self.cf.seed_uniform("Kr+", 1.0e18)
        np.testing.assert_allclose(self.cf.net_charge_density(), 0.0, atol=1e-6)

    def test_net_charge_density_electron_excess(self) -> None:
        from particles.constants import E_CHARGE

        self.cf.seed_uniform("e-", 1.0e18)
        rho = self.cf.net_charge_density()
        np.testing.assert_allclose(rho, -E_CHARGE * 1.0e18)


class TestPreionization(unittest.TestCase):
    """t=0 electron pre-ionization profiles n_{e,0}(x)."""

    def test_uniform(self) -> None:
        x = np.linspace(0.0, 12e-3, 16)
        np.testing.assert_allclose(Preionization(mode="uniform", density=1e18).profile(x), 1e18)

    def test_tanh_matches_paper_formula(self) -> None:
        # Chen et al. Eq. 6: n_e0(x) = 1e14*tanh(100 x) + 1e13 m^-3.
        x = np.linspace(0.0, 12e-3, 64)
        p = Preionization(mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0)
        np.testing.assert_allclose(p.profile(x), 1e14 * np.tanh(100.0 * x) + 1e13, rtol=1e-12)

    def test_tanh_endpoints_and_monotone(self) -> None:
        p = Preionization(mode="tanh")  # paper defaults
        self.assertAlmostEqual(float(p.profile(np.array([0.0]))[0]), 1e13)  # x=0: floor
        # x = 12 mm: 1e14*tanh(1.2) + 1e13
        self.assertAlmostEqual(float(p.profile(np.array([12e-3]))[0]),
                               1e14 * np.tanh(1.2) + 1e13, delta=1e9)
        v = p.profile(np.linspace(0.0, 12e-3, 32))
        self.assertTrue(np.all(np.diff(v) > 0.0))  # rises toward the cathode (x=length)

    def test_unknown_mode_raises(self) -> None:
        with self.assertRaises(ValueError):
            Preionization(mode="bogus").profile(np.array([0.0]))


if __name__ == "__main__":
    unittest.main()
