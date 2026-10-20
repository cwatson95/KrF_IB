"""Custom Maxwellian-EEDF coefficient solver: analytic limits, transport, breakdown physics.

This is the validated replacement for the placeholder ``AnalyticStubSolver`` — it computes real
rate and transport coefficients as EEDF moments of the project's measured SIGLO cross-sections.
"""

from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D.boltzmann import CoefficientTable, build_lut
from Fluid1D.eedf_solver import (
    GAMMA_V_PER_SQRT_EV,
    MaxwellianEedfSolver,
    maxwellian_f0,
)
from particles import excimer_system_registry
from reactions import ReactionNetwork, load_reaction_network
from reactions.cross_sections import CrossSectionRecord, TabulatedSigma
from reactions.types import RateForm, ReactionKind

XS_DIR = _ROOT / "data" / "cross_sections"
ARF_TOML = _ROOT / "data" / "reactions" / "arf.toml"
KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"
ARF_MIX = {"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001}


def _const_sigma(sigma_m2: float, kind: ReactionKind, *, threshold_ev=None) -> TabulatedSigma:
    """A synthetic energy-independent cross-section, for closed-form checks."""
    energy = np.geomspace(1.0e-4, 2.0e3, 400)
    record = CrossSectionRecord(
        process_id="synthetic", data_kind="total", reaction="e + X -> e + X",
        reaction_kind=kind, projectile="e", target="X", threshold_ev=threshold_ev,
        metadata={}, energy_ev=energy, sigma_m2=np.full_like(energy, sigma_m2),
    )
    return TabulatedSigma.from_record(record)


class TestMaxwellianEedf(unittest.TestCase):
    """The distribution itself: normalisation and its first moment."""

    def test_normalisation_and_mean_energy(self) -> None:
        eps = np.geomspace(1e-6, 400.0, 20000)
        for te in (0.5, 2.0, 7.0):
            f0 = maxwellian_f0(eps, te)
            self.assertAlmostEqual(float(np.trapezoid(f0 * np.sqrt(eps), eps)), 1.0, places=4)
            mean = float(np.trapezoid(f0 * eps**1.5, eps))
            self.assertAlmostEqual(mean / (1.5 * te), 1.0, places=4)  # ⟨ε⟩ = 3Tₑ/2


class TestAnalyticLimits(unittest.TestCase):
    """Closed forms for an energy-independent cross-section validate the quadrature."""

    def setUp(self) -> None:
        self.sigma0 = 3.0e-20
        reg = excimer_system_registry()
        net = ReactionNetwork(reg)
        net.add("Kr + e", "Kr+ + 2e", sigma="synthetic",
                rate_form=RateForm.ELECTRON_IMPACT, id="SYN")
        self.compiled = net.compile()
        self.solver = MaxwellianEedfSolver(XS_DIR)
        # Inject synthetic σ for both the reaction and the momentum transfer.
        self.solver._sigma_cache["synthetic"] = _const_sigma(self.sigma0, ReactionKind.IONIZATION)
        self.solver._elastic_sigma_for = lambda gas: _const_sigma(self.sigma0, ReactionKind.ELASTIC)

    def test_rate_matches_closed_form(self) -> None:
        # k = γ·σ0·(2/√π)·√Tₑ  for constant σ (∫ε e^{−ε/Tₑ}dε = Tₑ²)
        table = self.solver.solve({"Kr": 1.0}, np.zeros(32), reactions=self.compiled)
        for i in (5, 15, 25):
            te = table.electron_temp_ev[i]
            expected = GAMMA_V_PER_SQRT_EV * self.sigma0 * (2.0 / np.sqrt(np.pi)) * np.sqrt(te)
            self.assertAlmostEqual(table.rate_coeff["SYN"][i] / expected, 1.0, places=3)

    def test_mobility_matches_closed_form(self) -> None:
        # μ·N = 2γ/(3σ0√π)·Tₑ^{−1/2}  for constant σ_m
        table = self.solver.solve({"Kr": 1.0}, np.zeros(32), reactions=self.compiled)
        for i in (5, 15, 25):
            te = table.electron_temp_ev[i]
            expected = 2.0 * GAMMA_V_PER_SQRT_EV / (3.0 * self.sigma0 * np.sqrt(np.pi)) / np.sqrt(te)
            self.assertAlmostEqual(table.mobility_N[i] / expected, 1.0, places=3)


class TestRealCrossSections(unittest.TestCase):
    """Transport + rates from the shipped SIGLO datasets."""

    @classmethod
    def setUpClass(cls) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.net = load_reaction_network(ARF_TOML, excimer_system_registry()).compile()
            cls.table = MaxwellianEedfSolver(XS_DIR).solve(
                ARF_MIX, np.zeros(48), reactions=cls.net)

    def test_einstein_relation_is_exact(self) -> None:
        # D/μ = Tₑ by construction — the defining consistency check of the transport moments.
        ratio = self.table.diffusion_N / self.table.mobility_N
        np.testing.assert_allclose(ratio, self.table.electron_temp_ev, rtol=1e-12)

    def test_transport_is_physical(self) -> None:
        # Reduced mobility for rare gases sits in ~1e23-1e27 1/(V·m·s).
        self.assertTrue(np.all(self.table.mobility_N > 1.0e23))
        self.assertTrue(np.all(self.table.mobility_N < 1.0e27))
        self.assertTrue(np.all(np.diff(self.table.electron_temp_ev) > 0.0))
        self.assertTrue(np.all(self.table.diffusion_N > 0.0))

    def test_mobility_falls_as_momentum_transfer_grows(self) -> None:
        # Xe has a far larger momentum-transfer σ than Ne, so it is much less mobile.
        solver = MaxwellianEedfSolver(XS_DIR)
        ne = solver.solve({"Ne": 1.0}, np.zeros(24), reactions=None)
        xe = solver.solve({"Xe": 1.0}, np.zeros(24), reactions=None)
        self.assertLess(xe.mobility(3.0), ne.mobility(3.0))

    def test_rates_respect_thresholds(self) -> None:
        # Ne ionization (21.56 eV threshold) must be steeply rising and negligible when cold.
        cold = float(self.table.rate("R3", 0.3))
        warm = float(self.table.rate("R3", 12.0))
        self.assertLess(cold, 1.0e-30)
        self.assertGreater(warm, 1.0e-18)
        series = self.table.rate_coeff["R3"]
        self.assertTrue(np.all(np.diff(series) >= 0.0))  # monotone in ⟨ε⟩
        self.assertTrue(all(np.all(v >= 0.0) for v in self.table.rate_coeff.values()))

    def test_missing_cross_sections_are_reported_not_silent(self) -> None:
        # The warning machinery: a reaction whose dataset key resolves to nothing must be
        # reported, not silently zeroed. The shipped ArF network has no such gaps, so exercise
        # the warning synthetically.
        reg = excimer_system_registry()
        net = ReactionNetwork(reg)
        net.add("Ar + e", "Ar+ + 2e", sigma="no_such_dataset",
                rate_form=RateForm.ELECTRON_IMPACT, id="RX")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            MaxwellianEedfSolver(XS_DIR).solve(ARF_MIX, np.zeros(16), reactions=net.compile())
        msgs = " ".join(str(w.message) for w in caught)
        self.assertIn("without a usable cross-section", msgs)

    def test_arf_network_has_no_cross_section_gaps(self) -> None:
        # Every electron-impact row in arf.toml now resolves to a dataset: the last gap (R6,
        # Ar**) closed with e_Ar_excitation_2 — the NGFSRDW ground->4p sum (2p1..2p10).
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            MaxwellianEedfSolver(XS_DIR).solve(ARF_MIX, np.zeros(16), reactions=self.net)
        msgs = " ".join(str(w.message) for w in caught)
        self.assertNotIn("without a usable cross-section", msgs)

    def test_gas_without_momentum_transfer_warns(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            MaxwellianEedfSolver(XS_DIR).solve({"Ne": 0.5, "He": 0.5}, np.zeros(16))
        self.assertTrue(any("momentum-transfer" in str(w.message) for w in caught))

    def test_no_momentum_transfer_at_all_raises(self) -> None:
        with self.assertRaises(ValueError):
            MaxwellianEedfSolver(XS_DIR).solve({"He": 1.0}, np.zeros(16))


class TestBreakdownPhysics(unittest.TestCase):
    """The competition that decides whether a discharge can break down at all."""

    @staticmethod
    def _budget(table, mixture, eps, n_gas=8.3e25):
        ionize = sum(float(table.rate(rid, eps)) * n_gas * mixture.get(gas, 0.0)
                     for rid, gas in (("R3", "Ne"), ("R7", "Ar"), ("R17", "Xe")))
        attach = float(table.rate("R10", eps)) * n_gas * mixture.get("F2", 0.0)
        return ionize, attach

    @classmethod
    def setUpClass(cls) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.net = load_reaction_network(ARF_TOML, excimer_system_registry()).compile()
            cls.solver = MaxwellianEedfSolver(XS_DIR)
            cls.table = cls.solver.solve(ARF_MIX, np.zeros(48), reactions=cls.net)

    def test_attachment_dominates_cold_ionization_dominates_hot(self) -> None:
        ion_cold, att_cold = self._budget(self.table, ARF_MIX, 3.0)   # Tₑ = 2 eV
        self.assertLess(ion_cold, att_cold)                            # no avalanche possible
        ion_hot, att_hot = self._budget(self.table, ARF_MIX, 9.0)      # Tₑ = 6 eV
        self.assertGreater(ion_hot, 10.0 * att_hot)                    # runaway avalanche

    def test_more_f2_raises_the_breakdown_threshold(self) -> None:
        """Chen's central effect: F₂ suppresses breakdown, so the crossover moves up in Tₑ."""
        def crossover(mix):
            table = self.solver.solve(mix, np.zeros(48), reactions=self.net)
            for eps in np.linspace(1.0, 12.0, 200):
                ion, att = self._budget(table, mix, eps)
                if ion > att:
                    return eps
            return np.inf

        lean = crossover({**ARF_MIX, "F2": 0.001})
        rich = crossover({**ARF_MIX, "F2": 0.01})
        self.assertLess(lean, rich)
        self.assertTrue(np.isfinite(lean) and np.isfinite(rich))


class TestLocalFieldClosure(unittest.TestCase):
    """⟨ε⟩_eq(E/N): the local-field closure that repairs LMEA where nₑ → 0."""

    @classmethod
    def setUpClass(cls) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            net = load_reaction_network(ARF_TOML, excimer_system_registry()).compile()
            cls.table = MaxwellianEedfSolver(XS_DIR).solve(ARF_MIX, np.zeros(48), reactions=net)

    def test_reduced_field_supplied_and_monotone(self) -> None:
        self.assertIsNotNone(self.table.reduced_field_td)
        self.assertTrue(np.all(np.diff(self.table.reduced_field_td) > 0.0))

    def test_round_trips_against_the_tables_own_pairing(self) -> None:
        for i in (5, 20, 40):
            en = self.table.reduced_field_td[i]
            self.assertAlmostEqual(
                float(self.table.equilibrium_mean_energy(en)) / self.table.key_grid[i], 1.0, places=6)

    def test_hotter_field_gives_hotter_equilibrium(self) -> None:
        eps = [float(self.table.equilibrium_mean_energy(en)) for en in (0.5, 5.0, 15.4, 50.0)]
        self.assertTrue(all(b > a for a, b in zip(eps, eps[1:])))
        # At Chen's breakdown field the closure must sit above the ionization crossover.
        self.assertGreater(float(self.table.equilibrium_mean_energy(15.4)) * 2.0 / 3.0, 2.6)

    def test_mean_energy_bounds_match_the_grid(self) -> None:
        lo, hi = self.table.mean_energy_bounds_ev
        self.assertAlmostEqual(lo, float(self.table.key_grid.min()))
        self.assertAlmostEqual(hi, float(self.table.key_grid.max()))

    def test_absent_when_a_table_omits_it(self) -> None:
        """A table with no ⟨ε⟩↔E/N pairing reports None, so callers can fall back cleanly."""
        bare = CoefficientTable(
            key_kind="mean_energy_eV", key_grid=np.array([1.0, 2.0]),
            mobility_N=np.array([1e24, 1e24]), diffusion_N=np.array([1e24, 2e24]),
            drift_velocity=np.array([1.0, 2.0]), mean_energy_ev=np.array([1.0, 2.0]),
            electron_temp_ev=np.array([2 / 3, 4 / 3]), rate_coeff={}, mixture={"Kr": 1.0},
            gas_temp_k=300.0,
        )
        self.assertIsNone(bare.reduced_field_td)
        self.assertIsNone(bare.equilibrium_mean_energy(10.0))


class TestCacheAndIdentity(unittest.TestCase):
    def test_solver_id_is_stable_and_configuration_dependent(self) -> None:
        a = MaxwellianEedfSolver(XS_DIR)
        b = MaxwellianEedfSolver(XS_DIR)
        c = MaxwellianEedfSolver(XS_DIR, mean_energy_range_ev=(0.1, 40.0))
        self.assertEqual(a.solver_id, b.solver_id)
        self.assertNotEqual(a.solver_id, c.solver_id)

    def test_lut_cache_round_trip_preserves_everything(self) -> None:
        solver = MaxwellianEedfSolver(XS_DIR)
        with tempfile.TemporaryDirectory() as d, warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = np.zeros(24)
            t1 = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=grid, cache_dir=d,
                           key_kind="mean_energy_eV")
            t2 = build_lut({"Kr": 1.0}, solver=solver, en_grid_td=grid, cache_dir=d,
                           key_kind="mean_energy_eV")
            self.assertEqual(len(list(Path(d).glob("*.npz"))), 1)  # second call hit the cache
            np.testing.assert_allclose(t1.mobility_N, t2.mobility_N)
            np.testing.assert_allclose(t1.reduced_field_td, t2.reduced_field_td)
            self.assertIsInstance(t2, CoefficientTable)

    def test_energy_pair_left_to_the_five_thirds_closure(self) -> None:
        # Consumers rely on CoefficientTable applying (5/3) when the pair is absent.
        table = MaxwellianEedfSolver(XS_DIR).solve({"Kr": 1.0}, np.zeros(16))
        self.assertIsNone(table.energy_mobility_N)
        self.assertIsNone(table.energy_diffusion_N)
        self.assertAlmostEqual(table.energy_mobility(4.0) / table.mobility(4.0), 5.0 / 3.0, places=9)


if __name__ == "__main__":
    unittest.main()
