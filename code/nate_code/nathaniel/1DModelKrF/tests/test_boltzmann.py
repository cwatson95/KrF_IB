"""Boltzmann coefficient table + analytic stub solver + offline LUT cache."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D.boltzmann import (
    AnalyticStubSolver,
    CoefficientTable,
    build_lut,
    en_grid_from_spec,
    reduced_field_td,
)
from reactions import load_reaction_network
from particles import excimer_system_registry

KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"


class TestReducedField(unittest.TestCase):
    def test_townsend_conversion(self) -> None:
        en = reduced_field_td(1.0e6, 2.4e25)
        self.assertAlmostEqual(float(en), (1.0e6 / 2.4e25) / 1e-21, places=6)

    def test_zero_density_safe(self) -> None:
        self.assertEqual(float(reduced_field_td(1.0e6, 0.0)), 0.0)


class TestAnalyticStub(unittest.TestCase):
    def setUp(self) -> None:
        self.compiled = load_reaction_network(KRF_TOML, excimer_system_registry()).compile()
        self.grid = en_grid_from_spec((1.0, 1000.0, 32))
        self.table = AnalyticStubSolver().solve({"Kr": 0.99, "F2": 0.01}, self.grid,
                                                reactions=self.compiled)

    def test_table_shapes_and_monotone_te(self) -> None:
        self.assertIsInstance(self.table, CoefficientTable)
        self.assertEqual(self.table.electron_temp_ev.shape, self.grid.shape)
        self.assertTrue(np.all(np.diff(self.table.electron_temp_ev) >= 0))  # Te rises with E/N

    def test_interpolation_and_missing_rate(self) -> None:
        te = self.table.te(50.0)
        self.assertGreater(float(te), 0.0)
        self.assertEqual(float(self.table.rate("no_such_reaction", 50.0)), 0.0)

    def test_electron_impact_rates_present(self) -> None:
        # the KrF network has electron-impact reactions -> the stub gives them a rate curve
        self.assertTrue(len(self.table.reaction_ids()) > 0)
        for rid in self.table.reaction_ids():
            self.assertTrue(np.all(self.table.rate_coeff[rid] > 0.0))


class TestLutCache(unittest.TestCase):
    def test_cache_round_trip(self) -> None:
        grid = en_grid_from_spec((1.0, 500.0, 16))
        solver = AnalyticStubSolver()
        with tempfile.TemporaryDirectory() as d:
            t1 = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=grid, cache_dir=d)
            # second call must hit the cache and return identical arrays
            t2 = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=grid, cache_dir=d)
            np.testing.assert_allclose(t1.electron_temp_ev, t2.electron_temp_ev)
            np.testing.assert_allclose(t1.mobility_N, t2.mobility_N)
            self.assertEqual(len(list(Path(d).glob("*.npz"))), 1)

    def test_distinct_grid_interiors_do_not_collide(self) -> None:
        solver = AnalyticStubSolver()
        with tempfile.TemporaryDirectory() as d:
            a = np.array([1.0, 2.0, 10.0, 100.0])
            b = np.array([1.0, 50.0, 75.0, 100.0])
            ta = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=a, cache_dir=d)
            tb = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=b, cache_dir=d)
            np.testing.assert_array_equal(ta.key_grid, a)
            np.testing.assert_array_equal(tb.key_grid, b)
            self.assertEqual(len(list(Path(d).glob("*.npz"))), 2)


class TestMeanEnergyKeying(unittest.TestCase):
    """LMEA closure: re-key an E/N table onto mean energy + the energy-transport pair."""

    def setUp(self) -> None:
        net = load_reaction_network(KRF_TOML, excimer_system_registry()).compile()
        self.en = AnalyticStubSolver().solve(
            {"Kr": 0.99, "F2": 0.01}, en_grid_from_spec((1.0, 1000.0, 32)), reactions=net,
        )

    def test_rekey_to_mean_energy(self) -> None:
        me = self.en.to_mean_energy_key()
        self.assertEqual(me.key_kind, "mean_energy_eV")
        # abscissa is now the (sorted, ascending) mean-energy grid
        np.testing.assert_allclose(me.key_grid, np.sort(self.en.mean_energy_ev))
        self.assertTrue(np.all(np.diff(me.key_grid) >= 0.0))
        # idempotent
        self.assertIs(me.to_mean_energy_key(), me)

    def test_energy_transport_pair_default_closure(self) -> None:
        # μ_ε = (5/3)μ_e and D_ε = (5/3)D_e when the solver omits the pair.
        me = self.en.to_mean_energy_key()
        eps = 5.0
        self.assertAlmostEqual(me.energy_mobility(eps) / me.mobility(eps), 5.0 / 3.0, places=9)
        self.assertAlmostEqual(me.energy_diffusion(eps) / me.diffusion(eps), 5.0 / 3.0, places=9)

    def test_build_lut_mean_energy_key(self) -> None:
        net = load_reaction_network(KRF_TOML, excimer_system_registry()).compile()
        t = build_lut({"Kr": 0.99, "F2": 0.01}, solver=AnalyticStubSolver(),
                      en_grid_td=en_grid_from_spec((1.0, 1000.0, 16)), reactions=net,
                      key_kind="mean_energy_eV")
        self.assertEqual(t.key_kind, "mean_energy_eV")

    def test_energy_pair_survives_npz(self) -> None:
        me = self.en.to_mean_energy_key()
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "t.npz"
            me.save_npz(p)
            r = CoefficientTable.load_npz(p)
            self.assertEqual(r.key_kind, "mean_energy_eV")
            np.testing.assert_allclose(r.energy_mobility(5.0), me.energy_mobility(5.0))

    def test_saturated_mean_energy_ties_are_consolidated(self) -> None:
        table = AnalyticStubSolver(te_ceiling_ev=1.0).solve(
            {"Kr": 1.0}, np.array([1.0, 10.0, 100.0, 1000.0]),
        ).to_mean_energy_key()
        self.assertTrue(np.all(np.diff(table.key_grid) > 0.0))

    def test_mean_energy_table_can_be_requested_on_en_key(self) -> None:
        class NativeMeanSolver(AnalyticStubSolver):
            def solve(self, *args, **kwargs):
                return super().solve(*args, **kwargs).to_mean_energy_key()

        en = np.geomspace(1.0, 100.0, 8)
        table = build_lut({"Kr": 1.0}, solver=NativeMeanSolver(), en_grid_td=en)
        self.assertEqual(table.key_kind, "EN_Td")
        np.testing.assert_allclose(table.key_grid, en)


if __name__ == "__main__":
    unittest.main()
