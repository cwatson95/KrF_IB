"""0D excimer gain/output model (optics.gain.GainCavity, Chen Eqs. 10/14/18)."""

from __future__ import annotations

import math
import unittest

import numpy as np

from Fluid1D.charged_fluid import ChargedFluid
from kinetics.fluid import NeutralFluid
from optics.gain import GainCavity
from particles.constants import C_LIGHT, EV_TO_J
from particles.species import (
    AR, ARF_EXCIMER, AR_PLUS, AR_2STAR, ELECTRON, F, F_MINUS, photon_energy_ev,
)


def _cavity(**kw):
    base = dict(
        photon_energy_ev=photon_energy_ev(193.0),
        radiatives=[
            ("ArF*", 4.2e-9, ((1, "Ar"), (1, "F"))),
            ("Ar**", 1.0e-9, ((1, "Ar"),)),
        ],
        lasing=("ArF*", 4.0e-20, ((1, "Ar"), (1, "F")), "hv193"),
        absorbers=[("Ar**", 1.0e-22), ("F-", 1.0e-21)],
        gain_length_m=0.90, cavity_length_m=1.05, beam_area_m2=2.0e-3,
        r_oc=0.08, l_op=0.08, axial_se_fraction=1.0e-5,
    )
    base.update(kw)
    return GainCavity(**base)


class TestOutputRelations(unittest.TestCase):
    def test_output_fraction_zero_gain_limit(self) -> None:
        gc = _cavity()
        self.assertAlmostEqual(gc.output_fraction(0.0), (1.0 - 0.08) / (1.0 + 0.08), places=9)

    def test_output_fraction_matches_eq18(self) -> None:
        gc = _cavity()
        g = 3.0  # 1/m
        x = g * gc.l_g
        expected = x * math.exp(x) * (1 - 0.08) / ((1 + math.exp(x) * 0.08) * (math.exp(x) - 1))
        self.assertAlmostEqual(gc.output_fraction(g), expected, places=9)

    def test_output_power_units_eq10(self) -> None:
        gc = _cavity()
        gc._pd.n_ph = 1.0e20  # set a cavity density directly
        g = 2.0
        p = gc.output_power(g)
        e_photon_j = photon_energy_ev(193.0) * EV_TO_J
        expected = 1.0e20 * C_LIGHT * e_photon_j * gc.beam_area_m2 * gc.output_fraction(g)
        self.assertAlmostEqual(p / expected, 1.0, places=9)
        self.assertGreater(p, 0.0)  # W

    def test_total_photons_is_density_times_mode_volume(self) -> None:
        gc = _cavity()
        gc._pd.n_ph = 5.0e18
        # PhotonDensity is averaged over the full optical cavity (Chen Eq. 7), not only
        # the pumped gain length.
        expected = 5.0e18 * gc.l_c * gc.beam_area_m2
        self.assertAlmostEqual(gc.total_photons() / expected, 1.0, places=12)


class TestSpeciesBookkeeping(unittest.TestCase):
    def _fluid(self, n_arf=1.0e20, n_ar2s=1.0e19):
        fluid = NeutralFluid((AR, ARF_EXCIMER, AR_2STAR, F), 8, 1.0e-4, area_m2=2.0e-3)
        fluid.seed_uniform("ArF*", n_arf)
        fluid.seed_uniform("Ar**", n_ar2s)
        return fluid

    def test_spontaneous_decay_conserves_particles(self) -> None:
        # Pure spontaneous (lasing disabled -> no stimulated depletion): ArF* decays by
        # 1-e^{-dt/tau}, creating Ar + F 1:1.
        gc = _cavity(lasing=None, absorbers=[])
        fluid = self._fluid()
        arf0 = fluid.total_number("ArF*")
        ar0, f0 = fluid.total_number("Ar"), fluid.total_number("F")
        dt = 1.0e-9
        gc.advance(dt, fluid=fluid)
        lost = arf0 - fluid.total_number("ArF*")
        self.assertAlmostEqual(lost / arf0, -math.expm1(-dt / 4.2e-9), places=6)
        # Ar gets ArF* decay AND Ar** decay; F gets only ArF* decay.
        self.assertAlmostEqual((fluid.total_number("F") - f0) / lost, 1.0, places=6)
        self.assertGreater(fluid.total_number("Ar") - ar0, lost)  # + Ar** channel

    def test_advance_returns_power_and_gain(self) -> None:
        gc = _cavity()
        fluid = self._fluid(n_arf=5.0e20)
        p_out, g = gc.advance(1.0e-9, fluid=fluid)
        self.assertTrue(math.isfinite(p_out) and p_out >= 0.0)
        # g = n_ArF*·σ_STE − n_Ar**·σ (F- absent from this neutral-only fluid → 0).
        self.assertAlmostEqual(g, 5.0e20 * 4.0e-20 - 1.0e19 * 1.0e-22, delta=abs(g) * 1e-6 + 1e-9)

    def test_huge_gain_saturates_instead_of_overflowing(self) -> None:
        """A 100×-threshold inventory burns down in one finite co-integrated advance."""
        gc = _cavity()
        fluid = self._fluid(n_arf=2.0e22, n_ar2s=0.0)  # g·l_g ≈ 720 at start
        p_out, g = gc.advance(1.0e-9, fluid=fluid)
        self.assertAlmostEqual(g, 2.0e22 * 4.0e-20, delta=1e-3)  # pinned start-of-step gain
        self.assertTrue(math.isfinite(p_out) and p_out >= 0.0)
        self.assertTrue(math.isfinite(gc.n_ph) and gc.n_ph > 0.0)
        n_arf_end = fluid.total_number("ArF*") / (fluid.cell_volume * 8)
        self.assertLess(n_arf_end, 2.0e21)  # inventory burned toward saturation
        # and the output-fraction guard itself:
        self.assertEqual(gc.output_fraction(1.0e4), gc.output_fraction(1.0e4))  # finite, no raise
        self.assertGreaterEqual(gc.output_fraction(1.0e4), 0.0)

    def test_stimulated_growth_cannot_exceed_excimer_inventory(self) -> None:
        gc = _cavity(r_oc=0.999999, l_op=0.0, axial_se_fraction=0.0)
        fluid = self._fluid(n_arf=2.0e21, n_ar2s=0.0)
        gc._pd.n_ph = 1.0e16
        n_up0 = float(fluid.density[:, fluid.column["ArF*"]].mean())
        photons0 = gc.total_photons()
        gc.advance(2.0e-9, fluid=fluid)
        lost_upper = (
            n_up0 - float(fluid.density[:, fluid.column["ArF*"]].mean())
        ) * gc.l_g * gc.beam_area_m2
        created_photons = gc.total_photons() - photons0
        # Cavity loss and spontaneous decay can only reduce the photon side of this bound.
        self.assertLessEqual(created_photons, lost_upper * (1.0 + 2.0e-3))

    def test_no_lasing_channel_gives_zero_output(self) -> None:
        gc = _cavity(lasing=None, absorbers=[])
        fluid = self._fluid()
        p_out, g = gc.advance(1.0e-9, fluid=fluid)
        self.assertEqual(p_out, 0.0)
        self.assertEqual(g, 0.0)
        self.assertEqual(gc.total_photons(), 0.0)
        # the non-lasing Ar** radiative channel still decays it
        self.assertLess(fluid.total_number("Ar**"), 1.0e19 * 8 * fluid.cell_volume + 1.0)

    def test_radiative_decay_preserves_product_multiplicity(self) -> None:
        gc = _cavity(
            radiatives=[("ArF*", 1.0e-9, ((2, "Ar"),))],
            lasing=None,
            absorbers=[],
        )
        fluid = self._fluid(n_arf=5.0e20, n_ar2s=0.0)
        upper_before = fluid.density[:, fluid.column["ArF*"]].copy()
        ar_before = fluid.density[:, fluid.column["Ar"]].copy()
        gc.advance(1.0e-9, fluid=fluid)
        lost = upper_before - fluid.density[:, fluid.column["ArF*"]]
        np.testing.assert_allclose(fluid.density[:, fluid.column["Ar"]] - ar_before, 2.0 * lost)

    def test_photo_hazards_conserve_atoms_charge_photons_and_electron_energy(self) -> None:
        # Zero declared birth energy: photoelectrons dilute the pool and the electron
        # ENERGY DENSITY is conserved exactly across the photo update.
        channels = [
            ("R40", "Ar**", 1.0e-22, ((1, "Ar**"), (1, "hv193")),
             ((1, "Ar+"), (1, "e-")), 0.0, "hv193"),
            ("R41", "F-", 1.0e-21, ((1, "F-"), (1, "hv193")),
             ((1, "F"), (1, "e-")), 0.0, "hv193"),
        ]
        gc = _cavity(radiatives=[], absorbers=[], photo_channels=channels)
        fluid = NeutralFluid((AR, AR_2STAR, F), 4, 1.0e-4, area_m2=2.0e-3)
        charged = ChargedFluid(
            (ELECTRON, AR_PLUS, F_MINUS), 4, 1.0e-4, area_m2=2.0e-3,
        )
        fluid.seed_uniform("Ar**", 2.0e20)
        charged.seed_uniform("F-", 3.0e19)
        charged.seed_uniform("e-", 5.0e19)
        charged.seed_electron_energy(4.0)
        atom_ar0 = fluid.density[:, fluid.column["Ar**"]] + fluid.density[:, fluid.column["Ar"]]
        atom_ar0 = atom_ar0 + charged.density[:, charged.column["Ar+"]]
        atom_f0 = fluid.density[:, fluid.column["F"]] + charged.density[:, charged.column["F-"]]
        charge0 = charged.net_charge_density().copy()
        electron_energy0 = (
            charged.electron_density() * charged.electron_mean_energy_ev
        ).copy()

        fluence = 1.0e8
        absorber_before = (2.0e20, 3.0e19)
        sigmas = (1.0e-22, 1.0e-21)
        exact_extent = sum(
            n * (-math.expm1(-C_LIGHT * sigma * fluence))
            for n, sigma in zip(absorber_before, sigmas)
        )
        frozen_extent = sum(
            C_LIGHT * sigma * n * fluence
            for n, sigma in zip(absorber_before, sigmas)
        )
        photons0 = 1.0e21
        gc._pd.n_ph = photons0 - (gc.l_g / gc.l_c) * frozen_extent
        gc._advance_photo_channels(fluence, fluid, charged)

        atom_ar1 = fluid.density[:, fluid.column["Ar**"]] + fluid.density[:, fluid.column["Ar"]]
        atom_ar1 = atom_ar1 + charged.density[:, charged.column["Ar+"]]
        atom_f1 = fluid.density[:, fluid.column["F"]] + charged.density[:, charged.column["F-"]]
        np.testing.assert_allclose(atom_ar1, atom_ar0, rtol=1e-14)
        np.testing.assert_allclose(atom_f1, atom_f0, rtol=1e-14)
        np.testing.assert_allclose(charged.net_charge_density(), charge0, atol=1.0e-12)
        np.testing.assert_allclose(
            charged.electron_density() * charged.electron_mean_energy_ev,
            electron_energy0,
            rtol=1e-14,
        )
        self.assertAlmostEqual(
            gc._pd.n_ph,
            photons0 - (gc.l_g / gc.l_c) * exact_extent,
            delta=1.0e-13 * photons0,
        )
        self.assertEqual(gc.last_photoelectron_birth_energy_j, 0.0)

    def test_photoelectron_birth_energy_enters_the_pool_and_the_ledger(self) -> None:
        birth_r40, birth_r41 = 3.76, 3.02
        channels = [
            ("R40", "Ar**", 1.0e-22, ((1, "Ar**"), (1, "hv193")),
             ((1, "Ar+"), (1, "e-")), birth_r40, "hv193"),
            ("R41", "F-", 1.0e-21, ((1, "F-"), (1, "hv193")),
             ((1, "F"), (1, "e-")), birth_r41, "hv193"),
        ]
        gc = _cavity(radiatives=[], absorbers=[], photo_channels=channels)
        fluid = NeutralFluid((AR, AR_2STAR, F), 4, 1.0e-4, area_m2=2.0e-3)
        charged = ChargedFluid(
            (ELECTRON, AR_PLUS, F_MINUS), 4, 1.0e-4, area_m2=2.0e-3,
        )
        fluid.seed_uniform("Ar**", 2.0e20)
        charged.seed_uniform("F-", 3.0e19)
        charged.seed_uniform("e-", 5.0e19)
        charged.seed_electron_energy(4.0)
        energy0 = (charged.electron_density() * charged.electron_mean_energy_ev).copy()

        fluence = 1.0e8
        gc._pd.n_ph = 1.0e21
        gc._advance_photo_channels(fluence, fluid, charged)

        lost_r40 = 2.0e20 * (-math.expm1(-C_LIGHT * 1.0e-22 * fluence))
        lost_r41 = 3.0e19 * (-math.expm1(-C_LIGHT * 1.0e-21 * fluence))
        expected_gain_ev = lost_r40 * birth_r40 + lost_r41 * birth_r41
        energy1 = charged.electron_density() * charged.electron_mean_energy_ev
        # rtol allows the ~1e-12 cancellation noise of differencing two ~1e20-scale
        # energy densities to recover the ~1e15-scale birth-energy gain.
        np.testing.assert_allclose(energy1 - energy0, expected_gain_ev, rtol=1e-11)
        self.assertAlmostEqual(
            gc.last_photoelectron_birth_energy_j / (expected_gain_ev * EV_TO_J),
            1.0, places=12,
        )


class TestRadiativeSinkLedger(unittest.TestCase):
    """SE/STE upper-state sink ledger on advance() (feeds the ExcimerRates diagnostic)."""

    def _fluid(self, n_arf=1.0e20, n_ar2s=1.0e19):
        fluid = NeutralFluid((AR, ARF_EXCIMER, AR_2STAR, F), 8, 1.0e-4, area_m2=2.0e-3)
        fluid.seed_uniform("ArF*", n_arf)
        fluid.seed_uniform("Ar**", n_ar2s)
        return fluid

    def test_pure_spontaneous_ledger_is_exact(self) -> None:
        gc = _cavity(lasing=None, absorbers=[])
        fluid = self._fluid()
        dt = 1.0e-9
        gc.advance(dt, fluid=fluid)
        self.assertEqual(gc.last_dt_s, dt)
        self.assertAlmostEqual(
            gc.last_se_sink_m3["ArF*"] / (1.0e20 * -math.expm1(-dt / 4.2e-9)), 1.0,
            places=12,
        )
        self.assertAlmostEqual(
            gc.last_se_sink_m3["Ar**"] / (1.0e19 * -math.expm1(-dt / 1.0e-9)), 1.0,
            places=12,
        )
        self.assertEqual(gc.last_ste_sink_m3, {})

    def test_se_ste_split_sums_to_booked_upper_loss(self) -> None:
        # Bright cavity: the joint per-substep exponential loss is attributed by hazard
        # share, so SE + STE must reproduce the booked ArF* depletion exactly.
        gc = _cavity()
        fluid = self._fluid(n_arf=2.0e21, n_ar2s=0.0)
        gc._pd.n_ph = 1.0e16
        n0 = float(fluid.density[:, fluid.column["ArF*"]].mean())
        gc.advance(2.0e-9, fluid=fluid)
        n1 = float(fluid.density[:, fluid.column["ArF*"]].mean())
        se = gc.last_se_sink_m3["ArF*"]
        ste = gc.last_ste_sink_m3["ArF*"]
        self.assertGreater(se, 0.0)
        self.assertGreater(ste, 0.0)
        self.assertAlmostEqual((se + ste) / (n0 - n1), 1.0, places=10)
        # a fresh advance resets the ledger
        gc.advance(1.0e-12, fluid=fluid)
        self.assertLess(gc.last_se_sink_m3["ArF*"], se)


if __name__ == "__main__":
    unittest.main()
