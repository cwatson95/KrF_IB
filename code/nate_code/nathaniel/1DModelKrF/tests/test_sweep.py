"""F2-concentration sweep driver — the top-level use case."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    Preionization,
    Simulation,
    StubPlasmaModel,
    SweepFailure,
    sweep_f2_parallel,
)
from Fluid1D.sweep import (
    SweepResult,
    _base_point_payload,
    _simulation_for_point,
    f2_range,
    sweep_f2,
)

KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"


class TestF2Sweep(unittest.TestCase):
    def test_sweep_runs_and_collects_scalars(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            base = Simulation(
                name="sweep", grid=Grid1D(cells=8, length=2.0e-2, area_m2=1.0e-4),
                gas=GasMixture({"Kr": 0.99, "F2": 0.01}), plasma=Plasma(density=1e18, ion="Kr+"),
                cavity=Cavity(), boltzmann=BoltzmannSpec(backend="stub"), circuit=CircuitSpec(),
                reactions=KRF_TOML, diagnostics=[], output_dir=Path(d),
            )
            f2_values = [0.005, 0.01, 0.02]
            result = sweep_f2(f2_values, base,
                              numerics=Numerics(steps=20, dt=1e-12, fluid_every=10),
                              cache_dir=Path(d) / "lut_cache")
            self.assertIsInstance(result, SweepResult)
            self.assertEqual(len(result.points), 3)
            np.testing.assert_allclose(result.f2(), f2_values)
            self.assertTrue(np.all(np.isfinite(result.efficiency())))
            self.assertTrue(np.all(np.isfinite(result.output_energy_j())))
            best = result.best()
            self.assertIn(best.f2_fraction, f2_values)
            # one cached LUT per distinct mixture
            self.assertEqual(len(list((Path(d) / "lut_cache").glob("*.npz"))), 3)

    def test_with_fraction_rebalances_into_kr(self) -> None:
        gas = GasMixture({"Kr": 0.99, "F2": 0.01})
        g2 = gas.with_fraction("F2", 0.05, rebalance="Kr")
        self.assertAlmostEqual(g2.fractions["F2"], 0.05)
        self.assertAlmostEqual(g2.fractions["Kr"], 0.95)  # delta absorbed into Kr

    def test_production_sweep_is_exported_at_package_top_level(self) -> None:
        self.assertTrue(callable(sweep_f2_parallel))
        self.assertTrue(issubclass(SweepFailure, object))


class TestSweepConfiguration(unittest.TestCase):
    def _base(self, tmp: Path, **kwargs) -> Simulation:
        return Simulation(
            name="clone",
            grid=Grid1D(cells=8, length=2.0, area_m2=3.0),
            gas=GasMixture({"Kr": 0.99, "F2": 0.01}, pressure_atm=2.5,
                           temperature_k=345.0),
            plasma=Plasma(density=4.0e12, ion="Kr+", electron_mean_energy_ev=5.0),
            preionization=Preionization(mode="uniform", amplitude=4.0e12),
            cavity=Cavity(wavelength_nm=248.0, gain_length_m=4.0, cavity_length_m=5.0),
            boltzmann=BoltzmannSpec(backend="stub"),
            circuit=CircuitSpec(backend="stub", phi_right=-123.0),
            reactions=KRF_TOML,
            volume_units="cm",
            drive=object(),
            output_dir=tmp / "output",
            cross_sections_dir=tmp / "xs",
            ion_mobility_path=tmp / "ions.toml",
            diagnostics=[],
            **kwargs,
        )

    def test_shared_point_builder_preserves_all_simulation_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = self._base(Path(directory))
            diagnostics = object()
            plot = object()
            plasma_model = object()
            coeff_table = object()
            payload = _base_point_payload(
                base,
                rebalance="Kr",
                diagnostics_factory=lambda _f2: [diagnostics],
                plots_factory=lambda _f2: [plot],
                plasma_model_factory=lambda _f2: plasma_model,
                coeff_table_factory=lambda _f2: coeff_table,
            )
            with patch("Fluid1D.sweep.Simulation") as simulation:
                _simulation_for_point(payload, 0.02, name="point")

            kwargs = simulation.call_args.kwargs
            self.assertEqual(kwargs["name"], "point")
            self.assertIs(kwargs["grid"], base.grid_spec)
            self.assertEqual(kwargs["gas"].fractions, {"Kr": 0.98, "F2": 0.02})
            self.assertIs(kwargs["plasma"], base.plasma_spec)
            self.assertIs(kwargs["preionization"], base.preionization_spec)
            self.assertIs(kwargs["cavity"], base.cavity)
            self.assertEqual(kwargs["reactions"], base.reactions_path)
            self.assertIs(kwargs["boltzmann"], base.boltzmann_spec)
            self.assertIs(kwargs["circuit"], base.circuit_spec)
            self.assertEqual(kwargs["volume_units"], "cm")
            self.assertIs(kwargs["drive"], base.drive)
            self.assertEqual(kwargs["diagnostics"], [diagnostics])
            self.assertEqual(kwargs["plots"], [plot])
            self.assertIs(kwargs["plasma_model"], plasma_model)
            self.assertIs(kwargs["coeff_table"], coeff_table)
            self.assertEqual(kwargs["output_dir"], base.output_dir)
            self.assertEqual(kwargs["cross_sections_dir"], base.cross_sections_dir)
            self.assertEqual(kwargs["ion_mobility_path"], base.ion_mobility_path)

    def test_injected_runtime_objects_require_per_point_factories(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            model_base = self._base(tmp, plasma_model=StubPlasmaModel())
            with self.assertRaisesRegex(ValueError, "plasma_model_factory"):
                sweep_f2([0.01], model_base)
            table_base = self._base(tmp, coeff_table=object())
            with self.assertRaisesRegex(ValueError, "coeff_table_factory"):
                sweep_f2([0.01], table_base)

    def test_parallel_runtime_factory_must_be_picklable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = self._base(Path(directory))
            captured = object()
            with self.assertRaisesRegex(TypeError, "must be picklable"):
                sweep_f2_parallel(
                    base,
                    values=(0.01,),
                    plasma_model_factory=lambda _f2: captured,
                    progress=False,
                )


class TestSweepValidation(unittest.TestCase):
    def test_range_requires_finite_ascending_bounds_and_positive_step(self) -> None:
        self.assertEqual(f2_range(0.0, 0.2, 0.1), (0.0, 0.1, 0.2))
        for args in ((0.2, 0.1, 0.01), (-0.1, 0.1, 0.01), (0.9, 1.1, 0.01),
                     (0.0, 0.1, 0.0),
                     (0.0, 0.1, float("nan"))):
            with self.subTest(args=args), self.assertRaises(ValueError):
                f2_range(*args)

    def test_parallel_selection_and_worker_validation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = TestSweepConfiguration()._base(Path(directory))
            invalid_calls = (
                {"values": ()},
                {"values": (0.01, 0.01)},
                {"values": (-0.01,)},
                {"values": (1.01,)},
                {"values": (float("inf"),)},
                {"values": (0.01,), "start": 0.01, "stop": 0.02, "step": 0.01},
                {"start": 0.02, "stop": 0.01, "step": 0.01},
                {"values": (0.01,), "max_workers": 0},
            )
            for kwargs in invalid_calls:
                with self.subTest(kwargs=kwargs), self.assertRaises((TypeError, ValueError)):
                    sweep_f2_parallel(base, progress=False, **kwargs)


if __name__ == "__main__":
    unittest.main()
