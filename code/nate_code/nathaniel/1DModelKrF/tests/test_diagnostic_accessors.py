"""Public in-memory diagnostic accessors and pulse measurements."""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from Fluid1D.convergence import pulse_fwhm
from Fluid1D.diagnostics import (
    CircuitWaveforms,
    EnergyClosure,
    ExcimerRates,
    PlasmaScalarTraces,
    SpeciesDensities,
)


class TestDiagnosticAccessors(unittest.TestCase):
    def test_scalar_trace_arrays_are_copies_and_unknown_names_are_explicit(self) -> None:
        trace = PlasmaScalarTraces()
        trace._t.extend((1.0, 2.0))
        trace._ne.extend((3.0, 4.0))
        trace._en.extend((5.0, 6.0))
        trace._te.extend((7.0, 8.0))
        trace._pdis.extend((9.0, 10.0))
        trace._pl.extend((11.0, 12.0))

        times = trace.times_ns
        density = trace.series("n_e_m3")
        times[0] = -1.0
        density[0] = -1.0

        np.testing.assert_array_equal(trace.times_ns, (1.0, 2.0))
        np.testing.assert_array_equal(trace.series("n_e_m3"), (3.0, 4.0))
        with self.assertRaisesRegex(KeyError, "unknown diagnostic series.*available"):
            trace.series("missing")

    def test_species_and_circuit_accessors_preserve_units(self) -> None:
        species = SpeciesDensities()
        species._t.extend((20.0, 40.0))
        species._series["ArF*"] = [1.0e15, 2.0e15]
        np.testing.assert_array_equal(species.times_ps, (20.0, 40.0))
        np.testing.assert_array_equal(species.series("ArF*"), (1.0e15, 2.0e15))

        circuit = CircuitWaveforms()
        circuit._t.extend((1.0, 2.0))
        circuit._id.extend((-1.5, 2.5))
        np.testing.assert_array_equal(circuit.current_a, (-1500.0, 2500.0))
        current = circuit.current_a
        current[:] = 0.0
        np.testing.assert_array_equal(circuit.current_a, (-1500.0, 2500.0))

    def test_excimer_split_accessors_return_independent_arrays(self) -> None:
        rates = ExcimerRates(excimer="ArF*")
        rates._t.extend((1.0, 2.0))
        rates._ids = ("R_form", "R_quench")
        rates._nu_x = np.array((2.0, -1.0))
        rates._rows.extend((np.array((3.0, 4.0)), np.array((5.0, 6.0))))
        rates._se.extend((7.0, 8.0))
        rates._ste.extend((9.0, 10.0))

        formation = rates.formation_rates()
        quenching = rates.quenching_rates()
        np.testing.assert_array_equal(formation["R_form"], (6.0, 10.0))
        np.testing.assert_array_equal(quenching["R_quench"], (4.0, 6.0))
        np.testing.assert_array_equal(quenching["SE"], (7.0, 8.0))
        formation["R_form"][:] = 0.0
        np.testing.assert_array_equal(rates.formation_rates()["R_form"], (6.0, 10.0))

    def test_closure_totals_are_copied_and_expose_chemistry_energy(self) -> None:
        closure = EnergyClosure()
        sim = SimpleNamespace(_energy_closure_totals={
            "accepted_reaction_cost_j": 1.25,
            "coefficient_bound_adjustment_j": -0.5,
        })
        closure.export(sim, archive=None)

        totals = closure.totals
        totals["accepted_reaction_cost_j"] = 0.0
        self.assertEqual(closure.chemistry_energy_cost_j, 1.25)
        self.assertEqual(closure.chemistry_bound_adjustment_j, -0.5)
        self.assertEqual(closure.totals["accepted_reaction_cost_j"], 1.25)


class TestPulseFwhm(unittest.TestCase):
    def test_uses_the_time_axis_units_and_ignores_nonfinite_samples(self) -> None:
        times_ns = np.array((0.0, 1.0, np.nan, 2.0, 3.0, 4.0))
        power_w = np.array((0.0, 0.5, 100.0, 1.0, 0.5, 0.0))
        self.assertEqual(pulse_fwhm(times_ns, power_w), 2.0)

    def test_returns_zero_for_nonpositive_or_single_sample_peaks(self) -> None:
        self.assertEqual(pulse_fwhm((0.0, 1.0), (0.0, 0.0)), 0.0)
        self.assertEqual(pulse_fwhm((0.0, 1.0, 2.0), (0.0, 1.0, 0.0)), 0.0)

    def test_rejects_incompatible_shapes(self) -> None:
        with self.assertRaisesRegex(ValueError, "one-dimensional arrays of equal length"):
            pulse_fwhm((0.0, 1.0), (1.0,))


if __name__ == "__main__":
    unittest.main()
