"""Shared ArF/KrF sweep diagnostic and figure presets."""

from __future__ import annotations

import unittest

from Fluid1D.diagnostics import CircuitWaveforms, SpeciesDensities
from Fluid1D.sweep_presets import (
    STANDARD_FWHM_CE_OPTIONS,
    STANDARD_OVERLAY_FIGURES,
    standard_point_diagnostics,
    standard_point_plots,
)


class TestSweepPresets(unittest.TestCase):
    def test_conversion_efficiency_uses_bank_energy(self) -> None:
        self.assertEqual(STANDARD_FWHM_CE_OPTIONS["ce_denominator"], "bank")

    def test_standard_overlay_datasets_are_unique(self) -> None:
        filenames = [str(figure.filename) for figure in STANDARD_OVERLAY_FIGURES]
        self.assertEqual(len(filenames), len(set(filenames)))

    def test_diagnostics_include_panel_species_and_optional_circuit(self) -> None:
        panel_species = ("Ar+", "Ar*", "ArF*", "F-")
        diagnostics = standard_point_diagnostics(
            0.001,
            excimer="ArF*",
            circuit_backend="mpc",
            panel_species=panel_species,
        )
        self.assertTrue(any(isinstance(item, CircuitWaveforms) for item in diagnostics))
        densities = next(item for item in diagnostics if isinstance(item, SpeciesDensities))
        self.assertEqual(densities.track, panel_species)

    def test_point_plot_names_are_system_specific(self) -> None:
        plots = standard_point_plots(
            0.001,
            system="krf",
            formation_ids=("R28",),
            quenching_ids=("R33", "SE", "STE"),
            circuit_backend="stub",
        )
        self.assertEqual(
            [str(plot.filename) for plot in plots],
            ["krf_formation_rates.png", "krf_quenching_rates.png", "laser_output_power.png"],
        )


if __name__ == "__main__":
    unittest.main()
