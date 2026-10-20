"""Chemistry/electron-energy co-integration, bounds ledger, and call cadence."""

from __future__ import annotations

import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]

from Fluid1D.boltzmann import AnalyticStubSolver, build_lut, en_grid_from_spec
from Fluid1D.simulation import Simulation
from kinetics import FluidChemistry, FluidIntegrator
from particles import excimer_system_registry
from particles.constants import E_CHARGE
from reactions import load_reaction_network

_KRF = _ROOT / "data" / "reactions" / "krf.toml"
_ARF = _ROOT / "data" / "reactions" / "arf.toml"


def _plasma_chem_with_rates(key: float = 200.0):
    """Combined charged+neutral KrF chemistry with stub-LUT rates filled at ``key``."""
    reg = excimer_system_registry()
    net = load_reaction_network(_KRF, reg).compile()
    chem = FluidChemistry.for_plasma(net, reg, reg.names)
    coeffs = build_lut({"Kr": 0.999, "F2": 0.001}, solver=AnalyticStubSolver(),
                       en_grid_td=en_grid_from_spec((1.0, 1000.0, 32)), reactions=net)
    chem.set_dynamic_rates(coeffs, key)
    return reg, chem


def _plasma_chem_and_mean_energy_table():
    reg = excimer_system_registry()
    net = load_reaction_network(_KRF, reg).compile()
    chem = FluidChemistry.for_plasma(net, reg, reg.names)
    coeffs = build_lut(
        {"Kr": 0.999, "F2": 0.001},
        solver=AnalyticStubSolver(),
        en_grid_td=en_grid_from_spec((1.0, 1000.0, 32)),
        reactions=net,
    ).to_mean_energy_key()
    return chem, coeffs


def _arf_chem_and_mean_energy_table():
    reg = excimer_system_registry()
    net = load_reaction_network(_ARF, reg).compile()
    chem = FluidChemistry.for_plasma(net, reg, reg.names)
    coeffs = build_lut(
        {"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
        solver=AnalyticStubSolver(),
        en_grid_td=en_grid_from_spec((1.0, 1000.0, 32)),
        reactions=net,
    ).to_mean_energy_key()
    return chem, coeffs


def _seed(chem, n_e=1.0e18):
    n0 = np.zeros(chem.n_species)
    n0[chem.index["Kr"]] = 2.4e25 * 0.999
    n0[chem.index["F2"]] = 2.4e25 * 0.001
    n0[chem.index["e-"]] = n_e
    n0[chem.index["Kr+"]] = n_e
    return n0


class TestKernelEnergyCost(unittest.TestCase):
    def test_per_cell_requested_power_matches_reaction_progress(self) -> None:
        _, chem = _plasma_chem_with_rates()
        state = np.vstack((_seed(chem, 1.0e18), _seed(chem, 2.0e18)))
        k = np.tile(chem.k, (2, 1))
        requested = chem.energy_consumption_rate_percell(state, k)
        expected = np.array(
            [np.sum(chem.threshold_ev * chem.rates_of_progress(row)) for row in state]
        )
        np.testing.assert_allclose(requested, expected, rtol=2e-15)

    def test_cost_matches_radau_augmented_reference(self) -> None:
        """Σ ε_th·∫R dt from the kernel tracks an exact augmented Radau integration."""
        from scipy.integrate import solve_ivp

        _, chem = _plasma_chem_with_rates()
        integ = FluidIntegrator(chem)
        n0 = _seed(chem)
        eth = chem.threshold_ev
        interval = 1.0e-9

        def rhs_aug(_t, y):
            n = np.maximum(y[:-1], 0.0)
            r = chem.rates_of_progress(n)
            return np.append(chem.nu @ r, float(np.sum(eth * r)))

        sol = solve_ivp(rhs_aug, (0.0, interval), np.append(n0, 0.0),
                        method="Radau", rtol=1e-9, atol=1e-6)
        cost_ref = float(sol.y[-1, -1])
        self.assertGreater(cost_ref, 0.0)  # the scenario must actually spend energy

        for substeps, rtol in ((20, 0.15), (400, 0.02)):
            nn = np.array([n0.copy()])
            integ.advance(nn, interval, substeps=substeps)
            cost = float(integ.last_energy_cost_ev_m3[0])
            self.assertAlmostEqual(cost / cost_ref, 1.0, delta=rtol,
                                   msg=f"substeps={substeps}: {cost:.4g} vs {cost_ref:.4g}")

    def test_zero_thresholds_give_zero_cost(self) -> None:
        """A neutral-only chemistry (no threshold_ev) reports exactly zero spend."""
        reg = excimer_system_registry()
        net = load_reaction_network(_KRF, reg).compile()
        neutrals = tuple(s.name for s in reg if not s.is_charged and not s.is_photon)
        chem = FluidChemistry.for_fluid(net, reg, neutrals)
        integ = FluidIntegrator(chem)
        n0 = np.zeros((2, chem.n_species))
        n0[:, chem.index["Kr"]] = 2.4e25
        integ.advance(n0, 1e-10, substeps=4)
        np.testing.assert_array_equal(integ.last_energy_cost_ev_m3, 0.0)

    def test_hard_budget_cannot_be_overspent_inside_one_substep(self) -> None:
        _, chem = _plasma_chem_with_rates(key=1000.0)
        integ = FluidIntegrator(chem)
        budget = 2.5e19
        for substeps in (1, 8):
            state = np.array([_seed(chem, n_e=1.0e20)])
            integ.advance(
                state, 2.0e-9, substeps=substeps,
                energy_budget_ev_m3=np.array([budget]),
            )
            self.assertLessEqual(integ.last_energy_cost_ev_m3[0], budget * (1.0 + 1.0e-10))
            self.assertTrue(np.all(state >= 0.0))


class TestEnergyCoupledChemistry(unittest.TestCase):
    def test_source_partitions_are_disjoint_and_complete(self) -> None:
        reg = excimer_system_registry()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            compiled = load_reaction_network(_ARF, reg).compile()
        all_chem = FluidChemistry.for_plasma(compiled, reg, reg.names)
        electron = FluidChemistry.for_plasma(
            compiled, reg, reg.names, partition="electron_source",
        )
        heavy = FluidChemistry.for_plasma(
            compiled, reg, reg.names, partition="heavy_source",
        )
        electron_ids = set(electron.reaction_ids)
        heavy_ids = set(heavy.reaction_ids)
        self.assertFalse(electron_ids & heavy_ids)
        self.assertEqual(electron_ids | heavy_ids, set(all_chem.reaction_ids))
        self.assertLessEqual({"R1", "R9", "R10", "R18", "R20"}, electron_ids)

    def test_augmented_stage_uses_one_extent_for_species_and_energy(self) -> None:
        from kinetics.integrator import _augmented_linearly_implicit_trial

        # A + e -> B + e, with a 2 eV transfer per accepted event.  A large enough step
        # exercises the single positivity/energy damping factor.
        state = np.array([10.0, 0.0, 10.0])
        nu = np.array([[-1.0], [1.0], [0.0]])
        k = np.array([1.0])
        rs = np.array([[0, 2]], dtype=np.int64)
        ro = np.array([[1, 1]], dtype=np.int64)
        rc = np.array([2], dtype=np.int64)
        candidate, cost, alpha = _augmented_linearly_implicit_trial(
            state, 100.0, nu, k, rs, ro, rc, 1.0, np.array([2.0]), 2, 0.0,
        )
        extent = state[0] - candidate[0]
        self.assertGreater(extent, 0.0)
        self.assertLessEqual(alpha, 1.0)
        self.assertAlmostEqual(candidate[1], extent, places=12)
        self.assertAlmostEqual(cost, 2.0 * extent, places=12)

    def test_reaction_progress_and_energy_debit_close_per_cell(self) -> None:
        chem, coeffs = _plasma_chem_and_mean_energy_table()
        integ = FluidIntegrator(chem)
        state = np.vstack((_seed(chem, 1.0e18), _seed(chem, 2.0e18)))
        initial_j = np.array((8.0e18, 16.0e18)) * E_CHARGE
        initial_ev = initial_j / E_CHARGE

        integ.advance_energy_coupled(
            state,
            2.0e-10,
            electron_energy_j_m3=initial_j,
            electron_col=chem.index["e-"],
            coeffs=coeffs,
            energy_bounds_ev=coeffs.mean_energy_bounds_ev,
            substeps=16,
        )

        final_ev = integ.last_electron_energy_j_m3 / E_CHARGE
        cost = integ.last_energy_cost_ev_m3
        adjustment = integ.last_energy_bound_adjustment_ev_m3
        np.testing.assert_allclose(final_ev + cost - adjustment, initial_ev, rtol=2e-13)
        self.assertGreater(float(np.max(cost)), 0.0)
        mean_energy = final_ev / np.maximum(state[:, chem.index["e-"]], 1.0)
        self.assertTrue(np.all(mean_energy >= coeffs.mean_energy_bounds_ev[0] - 1e-12))
        self.assertTrue(np.all(mean_energy <= coeffs.mean_energy_bounds_ev[1] + 1e-12))

    def test_live_energy_refresh_cools_the_dynamic_rates(self) -> None:
        chem, coeffs = _plasma_chem_and_mean_energy_table()
        integ = FluidIntegrator(chem)
        initial = _seed(chem, 1.0e19)
        coupled = np.array([initial.copy()])
        initial_mean_ev = 12.0
        initial_j = np.array([initial_mean_ev * initial[chem.index["e-"]] * E_CHARGE])
        integ.advance_energy_coupled(
            coupled,
            5.0e-10,
            electron_energy_j_m3=initial_j,
            electron_col=chem.index["e-"],
            coeffs=coeffs,
            energy_bounds_ev=coeffs.mean_energy_bounds_ev,
            substeps=40,
        )
        final_mean_ev = (
            integ.last_electron_energy_j_m3[0]
            / (max(coupled[0, chem.index["e-"]], 1.0) * E_CHARGE)
        )
        self.assertLess(final_mean_ev, initial_mean_ev)
        # At least one threshold reaction must have a lower coefficient at the cooled state.
        costly_dynamic = np.flatnonzero(chem.dynamic_mask & chem.energy_costly_mask)
        ratios = []
        for r in costly_dynamic:
            rid = chem.reaction_ids[r]
            k0 = float(coeffs.rate(rid, initial_mean_ev))
            k1 = float(coeffs.rate(rid, final_mean_ev))
            if k0 > 0.0:
                ratios.append(k1 / k0)
        self.assertTrue(ratios)
        self.assertLess(min(ratios), 1.0)

    def test_energy_floor_preserves_threshold_free_penning_progress(self) -> None:
        """The energy limiter must not erase R20 while all costly channels are stopped."""
        chem, coeffs = _arf_chem_and_mean_energy_table()
        integ = FluidIntegrator(chem)
        state = np.zeros((1, chem.n_species))
        initial_ne = 1.0e15
        state[0, chem.index["e-"]] = initial_ne
        state[0, chem.index["Ne+"]] = initial_ne
        state[0, chem.index["Ne*"]] = 1.0e20
        floor_ev = coeffs.mean_energy_bounds_ev[0]
        initial_j = np.array([floor_ev * initial_ne * E_CHARGE])

        integ.advance_energy_coupled(
            state,
            1.0e-9,
            electron_energy_j_m3=initial_j,
            electron_col=chem.index["e-"],
            coeffs=coeffs,
            energy_bounds_ev=coeffs.mean_energy_bounds_ev,
            substeps=8,
        )

        final_ne = state[0, chem.index["e-"]]
        self.assertGreater(final_ne, initial_ne)
        self.assertGreater(integ.last_energy_bound_adjustment_ev_m3[0], 0.0)
        self.assertAlmostEqual(integ.last_energy_cost_ev_m3[0], 0.0, delta=1.0e-10)
        np.testing.assert_allclose(
            integ.last_electron_energy_j_m3[0] / E_CHARGE,
            floor_ev * final_ne,
            rtol=2e-13,
        )

    def test_power_source_extends_the_per_cell_identity(self) -> None:
        """With an in-stage source: final + cost - adjustment - delivered == initial."""
        chem, coeffs = _plasma_chem_and_mean_energy_table()
        integ = FluidIntegrator(chem)
        state = np.vstack((_seed(chem, 1.0e18), _seed(chem, 2.0e18)))
        initial_ev = np.array((8.0e18, 16.0e18))
        interval = 2.0e-10
        # One strongly heated cell, one strongly drained cell (negative source).
        source_w = np.array((5.0e8, -2.0e8)) * E_CHARGE  # eV·m⁻³·s⁻¹ -> W/m³

        integ.advance_energy_coupled(
            state,
            interval,
            electron_energy_j_m3=initial_ev * E_CHARGE,
            electron_col=chem.index["e-"],
            coeffs=coeffs,
            energy_bounds_ev=coeffs.mean_energy_bounds_ev,
            substeps=16,
            power_source_j_m3_s=source_w,
        )

        final_ev = integ.last_electron_energy_j_m3 / E_CHARGE
        cost = integ.last_energy_cost_ev_m3
        adjustment = integ.last_energy_bound_adjustment_ev_m3
        delivered_ev = source_w / E_CHARGE * interval
        np.testing.assert_allclose(
            final_ev + cost - adjustment - delivered_ev, initial_ev, rtol=1e-12,
        )
        # The drained cell's overdraw is absorbed and ledgered at the floor, never hidden.
        mean_energy = final_ev / np.maximum(state[:, chem.index["e-"]], 1.0)
        self.assertTrue(np.all(mean_energy >= coeffs.mean_energy_bounds_ev[0] - 1e-12))

    def test_power_source_holds_the_stage_at_its_balance_point(self) -> None:
        """Heating inside the stage: the state settles where demand matches the source
        instead of draining the stored buffer to the coefficient floor."""
        chem, coeffs = _plasma_chem_and_mean_energy_table()
        interval = 2.0e-9
        e_col = chem.index["e-"]

        def run(source_w: float) -> np.ndarray:
            integ = FluidIntegrator(chem)
            state = np.array([_seed(chem, 1.0e18)])
            initial_j = np.array([6.0 * state[0, e_col] * E_CHARGE])
            integ.advance_energy_coupled(
                state, interval,
                electron_energy_j_m3=initial_j,
                electron_col=e_col,
                coeffs=coeffs,
                energy_bounds_ev=coeffs.mean_energy_bounds_ev,
                substeps=400,
                power_source_j_m3_s=(
                    None if source_w == 0.0 else np.array([source_w])
                ),
            )
            eps = integ.last_electron_energy_j_m3[0] / (
                max(state[0, e_col], 1.0) * E_CHARGE
            )
            return np.array([eps, *state[0]])

        unsourced = run(0.0)
        # Source sized to the initial demand scale so a balance exists inside the table.
        k0 = chem.dynamic_rates_percell(coeffs, np.array([6.0]))
        demand0 = float(chem.energy_consumption_rate_percell(
            np.array([_seed(chem, 1.0e18)]), k0,
        )[0]) * E_CHARGE
        self.assertGreater(demand0, 0.0)
        sourced = run(0.5 * demand0)

        eps_unsourced, eps_sourced = unsourced[0], sourced[0]
        # Without heating the stage drains toward the floor; with it, it holds hot.
        self.assertGreater(eps_sourced, eps_unsourced + 0.5)
        # And the final state sits near the demand == source balance point.
        state_final = sourced[1:].reshape(1, -1)
        k_final = chem.dynamic_rates_percell(coeffs, np.array([eps_sourced]))
        demand_final = float(chem.energy_consumption_rate_percell(
            state_final, k_final,
        )[0]) * E_CHARGE
        self.assertLess(abs(demand_final - 0.5 * demand0), 0.35 * 0.5 * demand0)


class TestElectronEnergyHandoff(unittest.TestCase):
    def _resolved_model(self):
        from tests.test_plasma_model import _dd_setup
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        cf.seed_uniform("e-", 1e16)
        cf.seed_uniform("Kr+", 1e16)
        cf.seed_electron_energy(3.0)
        e = np.zeros((n, 3))
        e[:, 0] = 1.0e6
        model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        return cf, model, n

    def test_adopts_coupled_energy_and_redistributes(self) -> None:
        cf, model, n = self._resolved_model()
        # Emulate a coupled chemistry call: electrons doubled and the returned energy state
        # lies inside the coefficient table's resolved range.
        cf.density[:, model._e_col] *= 2.0
        n_e_post = cf.density[:, model._e_col]
        target_ev = 0.5 * (model._eps_bounds[0] + min(3.0, model._eps_bounds[1]))
        final_energy = target_ev * n_e_post * E_CHARGE
        model.adopt_electron_energy_density(final_energy, charged_fluid=cf)
        expected = final_energy / (n_e_post * E_CHARGE)
        np.testing.assert_allclose(cf.electron_mean_energy(), expected, rtol=1e-12)
        self.assertLess(float(cf.electron_mean_energy().mean()), 3.0)

    def test_rejects_energy_outside_the_resolved_bounds(self) -> None:
        cf, model, _ = self._resolved_model()
        excessive = cf.density[:, model._e_col] * (model._eps_bounds[1] + 1.0) * E_CHARGE
        with self.assertRaisesRegex(ValueError, "outside coefficient-table bounds"):
            model.adopt_electron_energy_density(excessive, charged_fluid=cf)


class TestAdaptiveChemistryCadence(unittest.TestCase):
    def _due(self, i, *, nu, fluid_every=20, dt=50e-12, last=0):
        ns = SimpleNamespace(
            numerics=SimpleNamespace(fluid_every=fluid_every, dt=dt),
            _last_nu_max=nu, _chem_efolds_per_call=0.05, _last_chem_step=last,
        )
        return Simulation._chemistry_due(ns, i)

    def test_quiescent_plasma_keeps_the_nominal_cadence(self) -> None:
        # ν·dt = 1e6·5e-11 = 5e-5 e-folds/step → target 1000 steps: only the nominal cadence.
        self.assertFalse(self._due(7, nu=1.0e6))
        self.assertTrue(self._due(20, nu=1.0e6))

    def test_conducting_phase_grades_to_intermediate_cadence(self) -> None:
        # ν = 1e8 → 0.05/(1e8·5e-11) = 10 steps: attachment churn must not wait 20 steps
        # (the F2-sweep sawtooth), nor force every-step cost.
        self.assertFalse(self._due(7, nu=1.0e8))
        self.assertTrue(self._due(12, nu=1.0e8))

    def test_breakdown_tightens_to_every_step(self) -> None:
        # ν = 5e9 → target 0.2 steps: every step.
        self.assertTrue(self._due(7, nu=5.0e9))
        self.assertTrue(self._due(8, nu=5.0e9))

    def _due_heavy(self, i, *, heavy_nu, fluid_every=20, dt=50e-12, last=0):
        ns = SimpleNamespace(
            numerics=SimpleNamespace(
                fluid_every=fluid_every, dt=dt,
                electron_energy_closure="source_coupled",
            ),
            _real_transport=True,
            _last_heavy_nu=heavy_nu, _chem_efolds_per_call=0.05, _last_chem_step=last,
        )
        return Simulation._chemistry_due(ns, i)

    def test_source_coupled_heavy_cadence_grades_on_ion_stiffness(self) -> None:
        # Quiescent heavy rows: only the nominal cadence.
        self.assertFalse(self._due_heavy(7, heavy_nu=0.0))
        self.assertTrue(self._due_heavy(20, heavy_nu=0.0))
        # Conduction-phase ion-ion recombination ν = 1e9 → 0.05/(1e9·5e-11) = 1 step.
        self.assertTrue(self._due_heavy(1, heavy_nu=1.0e9))
        # ν = 1e8 → 10-step grade.
        self.assertFalse(self._due_heavy(7, heavy_nu=1.0e8))
        self.assertTrue(self._due_heavy(12, heavy_nu=1.0e8))


class TestChargedChangingMask(unittest.TestCase):
    def test_heavy_partition_flags_ion_channels_not_neutral_quench(self) -> None:
        reg = excimer_system_registry()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            compiled = load_reaction_network(_ARF, reg).compile()
        heavy = FluidChemistry.for_plasma(
            compiled, reg, reg.names, partition="heavy_source",
        )
        mask = heavy.charged_changing_mask
        self.assertIsNotNone(mask)
        by_id = dict(zip(heavy.reaction_ids, mask))
        # Ion-ion recombination consumes Ar+/F-: charged-changing.
        self.assertTrue(by_id["R30"])
        # Neutral harpoon/quench rows change no charged species.
        self.assertFalse(by_id["R29"])
        self.assertFalse(by_id["R33"])


if __name__ == "__main__":
    unittest.main()
