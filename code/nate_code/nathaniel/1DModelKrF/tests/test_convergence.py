"""Refinement-gate acceptance semantics."""

from __future__ import annotations

import unittest

from Fluid1D.convergence import ConvergenceMetrics, compare_refinement


def _metric(**changes) -> ConvergenceMetrics:
    values = dict(
        cells=256,
        dt_s=50e-12,
        duration_s=200e-9,
        peak_electron_density_m3=1.0e20,
        deposited_energy_j=1.4,
        output_energy_j=0.08,
        peak_output_power_w=5.0e6,
        pulse_fwhm_s=12e-9,
    )
    values.update(changes)
    return ConvergenceMetrics(**values)


class TestConvergenceGate(unittest.TestCase):
    def test_accepts_refinement_inside_all_limits(self) -> None:
        result = compare_refinement(_metric(), _metric(cells=512, output_energy_j=0.085))
        self.assertTrue(result.passed)
        self.assertEqual(result.failures, ())

    def test_rejects_material_output_and_fwhm_changes(self) -> None:
        result = compare_refinement(
            _metric(), _metric(cells=512, output_energy_j=0.05, pulse_fwhm_s=10e-9),
        )
        self.assertFalse(result.passed)
        self.assertTrue(any("output_energy_rel" in item for item in result.failures))
        self.assertTrue(any("pulse_fwhm_abs_s" in item for item in result.failures))

    def test_nonlasing_pair_still_gates_plasma_and_deposition(self) -> None:
        result = compare_refinement(
            _metric(output_energy_j=0.0, peak_output_power_w=0.0, pulse_fwhm_s=0.0),
            _metric(cells=512, output_energy_j=0.0, peak_output_power_w=0.0,
                    pulse_fwhm_s=0.0, peak_electron_density_m3=5.0e19),
        )
        self.assertFalse(result.passed)
        self.assertTrue(any("peak_electron_density_rel" in item for item in result.failures))


if __name__ == "__main__":
    unittest.main()
