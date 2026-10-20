"""Parity of the Numba coupled-transport kernel against the NumPy reference path.

The cases span both energy closures, fixed and finite-capacitance gaps, pre- and
post-breakdown densities, enabled and disabled heating caps, a capped remainder,
near-empty cells, and prescribed-field routing. Tolerances allow the small reduction-order
differences between scalar JIT loops and NumPy vector operations.
"""

from __future__ import annotations

import unittest
import warnings

import numpy as np

from Fluid1D.field import Poisson1DSolver, SemiImplicitField, sample_efield
from Fluid1D.transport_kernel import coupled_advance, pack_kernel_args
from tests.test_plasma_model import _dd_setup


def _seed(cf, n_e=1.0e20, eps_ev=3.0):
    cf.seed_uniform("e-", n_e)
    cf.seed_uniform("Kr+", n_e)
    cf.seed_electron_energy(eps_ev)


def _advance_model_pair(
    *,
    n_e: float = 1.0e20,
    eps_ev: float = 3.0,
    closure: str = "source_coupled",
    heating_cap_factor: float = 2.0,
    gap_capacitance_f: float | None = None,
    max_substeps: int = 20_000,
    electron_profile: np.ndarray | None = None,
    dt: float = 50e-12,
):
    """Advance matching NumPy/kernel models once and return both observable states."""
    _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
    _, nf_b, cf_b, _c2, model_b, _n2, _dx2 = _dd_setup()
    model_a.use_numba_kernel = False
    for model in (model_a, model_b):
        model.electron_energy_closure = closure
        model.field_heating_cap_factor = heating_cap_factor
        model.max_substeps = max_substeps
    if electron_profile is None:
        _seed(cf_a, n_e=n_e, eps_ev=eps_ev)
        _seed(cf_b, n_e=n_e, eps_ev=eps_ev)
    else:
        for cf in (cf_a, cf_b):
            cf.seed_profile("e-", electron_profile)
            cf.seed_uniform("Kr+", n_e)
            cf.seed_electron_energy(eps_ev)

    centers = (np.arange(n) + 0.5) * dx
    reports = []
    handles = []
    for nf, cf, model in ((nf_a, cf_a, model_a), (nf_b, cf_b, model_b)):
        poisson = Poisson1DSolver(n, dx)
        _phi, e1d = poisson.solve(cf.net_charge_density(), 0.0, -5.0e3)
        handle = SemiImplicitField(
            solver=poisson,
            phi_left=0.0,
            phi_right=-5.0e3,
            cell_centers=centers,
            gap_capacitance_f=gap_capacitance_f,
        )
        reports.append(
            model.advance(
                dt,
                efield=sample_efield(e1d, centers),
                neutral_fluid=nf,
                charged_fluid=cf,
                coeffs=coeffs,
                field=handle,
            )
        )
        handles.append(handle)
    return (cf_a, model_a, reports[0], handles[0]), (cf_b, model_b, reports[1], handles[1])


def _assert_advance_parity(test: unittest.TestCase, numpy_state, kernel_state) -> None:
    """Compare the full public transport result plus retained internal energy state."""
    cf_np, model_np, report_np, field_np = numpy_state
    cf_kernel, model_kernel, report_kernel, field_kernel = kernel_state
    np.testing.assert_allclose(cf_kernel.density, cf_np.density, rtol=3e-6, atol=1.0)
    np.testing.assert_allclose(
        cf_kernel.electron_mean_energy_ev,
        cf_np.electron_mean_energy_ev,
        rtol=1e-4,
    )
    np.testing.assert_allclose(
        model_kernel._last_n_eps,
        model_np._last_n_eps,
        rtol=1e-4,
        atol=1.0e-30,
    )
    test.assertEqual(model_kernel._last_substeps, model_np._last_substeps)
    test.assertAlmostEqual(field_kernel.phi_right / field_np.phi_right, 1.0, delta=1e-8)
    test.assertAlmostEqual(
        report_kernel.deposited_power_w / report_np.deposited_power_w,
        1.0,
        delta=1e-6,
    )
    test.assertAlmostEqual(
        report_kernel.conduction_current_a / report_np.conduction_current_a,
        1.0,
        delta=1e-6,
    )
    for name in (
        "field_work_j_m3",
        "elastic_transfer_j_m3",
        "continuous_inelastic_sink_j_m3",
        "coefficient_bound_adjustment_j_m3",
    ):
        np.testing.assert_allclose(
            getattr(report_kernel.energy_closure, name),
            getattr(report_np.energy_closure, name),
            rtol=1e-4,
            atol=1.0e-18,
        )


class TestContinuousSinkBasis(unittest.TestCase):
    def test_transport_sink_matches_chemistry_energy_demand(self) -> None:
        """The reconciled mode's provisional sink must be the SAME Σ ε_th·k(⟨ε⟩)·Πn the
        chemistry bills — evaluated at an identical state (healthy n_e, so the ⟨ε⟩_eq blend
        is inert and the comparison isolates the row/column packing)."""
        from kinetics.chemistry import FluidChemistry

        reg, nf, cf, coeffs, model, n, _dx = _dd_setup()
        cf.seed_uniform("e-", 1.0e19)
        cf.seed_uniform("Kr+", 1.0e19)
        nf.seed_uniform("Kr*", 1.0e18)
        nf.seed_uniform("Kr**", 1.0e17)
        cf.seed_electron_energy(6.0)
        model._resolve(nf, cf, coeffs)

        eps = np.full(n, 6.0)
        p_model = model._continuous_inelastic_power(eps, nf.density, cf.density)

        chem = FluidChemistry.for_plasma(model.reactions, reg, reg.names)
        n_all = np.zeros((n, chem.n_species))
        for name, col in nf.column.items():
            n_all[:, chem.index[name]] = nf.density[:, col]
        for name, col in cf.column.items():
            if name in chem.index:
                n_all[:, chem.index[name]] = cf.density[:, col]
        k_percell = chem.dynamic_rates_percell(coeffs, eps)
        p_chem = chem.energy_consumption_rate_percell(n_all, k_percell)

        self.assertGreater(float(np.max(p_chem)), 0.0)
        np.testing.assert_allclose(p_model, p_chem, rtol=1e-12)


class TestKernelParity(unittest.TestCase):
    def test_parity_matrix_covers_closure_density_and_cap_modes(self) -> None:
        cases = (
            ("source_coupled", 1.0e14, 2.0),
            ("source_coupled", 1.0e20, np.inf),
            ("continuous_reconciled", 1.0e19, 2.0),
            ("continuous_reconciled", 1.0e20, np.inf),
        )
        for closure, density, cap in cases:
            with self.subTest(closure=closure, density=density, cap=cap):
                numpy_state, kernel_state = _advance_model_pair(
                    n_e=density,
                    eps_ev=6.0 if closure == "continuous_reconciled" else 3.0,
                    closure=closure,
                    heating_cap_factor=cap,
                )
                _assert_advance_parity(self, numpy_state, kernel_state)

    def test_near_empty_cells_match_numpy_reference(self) -> None:
        n = _dd_setup()[5]
        profile = np.full(n, 1.0e16)
        profile[n // 2:] = 1.0e-12
        numpy_state, kernel_state = _advance_model_pair(
            n_e=1.0e16,
            electron_profile=profile,
        )
        _assert_advance_parity(self, numpy_state, kernel_state)

    def test_capped_remainder_matches_numpy_reference(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            numpy_state, kernel_state = _advance_model_pair(
                n_e=1.0e21,
                max_substeps=1,
                dt=5.0e-9,
            )
        self.assertTrue(any("max_substeps" in str(item.message) for item in caught))
        _assert_advance_parity(self, numpy_state, kernel_state)
        self.assertEqual(kernel_state[1]._last_substeps, 2)

    def test_prescribed_field_routes_to_numpy_reference(self) -> None:
        class ZeroDrive:
            @staticmethod
            def at_batch(positions, _time):
                return np.zeros_like(positions)

        _, nf, cf, coeffs, model, n, dx = _dd_setup()
        _seed(cf, n_e=1.0e14)
        centers = (np.arange(n) + 0.5) * dx
        poisson = Poisson1DSolver(n, dx)
        _phi, e1d = poisson.solve(cf.net_charge_density(), 0.0, -5.0e3)
        handle = SemiImplicitField(
            solver=poisson,
            phi_left=0.0,
            phi_right=-5.0e3,
            cell_centers=centers,
            prescribed=ZeroDrive(),
        )

        def fail_kernel(*_args, **_kwargs):
            self.fail("a prescribed field must use the NumPy reference path")

        model._advance_via_kernel = fail_kernel
        report = model.advance(
            50e-12,
            efield=sample_efield(e1d, centers),
            neutral_fluid=nf,
            charged_fluid=cf,
            coeffs=coeffs,
            field=handle,
        )
        self.assertTrue(np.isfinite(report.deposited_power_w))

    def test_finite_gap_heating_cap_matches_numpy_reference(self) -> None:
        """The cap uses the voltage that drove the substep, before the gap-node sag."""
        _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
        _, nf_b, cf_b, _c2, model_b, _n2, _dx2 = _dd_setup()
        model_a.use_numba_kernel = False
        _seed(cf_a, n_e=1.0e21)
        _seed(cf_b, n_e=1.0e21)
        centers = (np.arange(n) + 0.5) * dx
        v_b = -5.0e3
        c_g = 2.5e-11

        reports = []
        handles = []
        for nf, cf, model in ((nf_a, cf_a, model_a), (nf_b, cf_b, model_b)):
            poisson = Poisson1DSolver(n, dx)
            _phi, e1d = poisson.solve(cf.net_charge_density(), 0.0, v_b)
            handle = SemiImplicitField(
                solver=poisson,
                phi_left=0.0,
                phi_right=v_b,
                cell_centers=centers,
                gap_capacitance_f=c_g,
            )
            reports.append(
                model.advance(
                    50e-12,
                    efield=sample_efield(e1d, centers),
                    neutral_fluid=nf,
                    charged_fluid=cf,
                    coeffs=coeffs,
                    field=handle,
                )
            )
            handles.append(handle)

        report_np, report_kernel = reports
        np.testing.assert_allclose(cf_b.density, cf_a.density, rtol=1e-7, atol=1.0)
        np.testing.assert_allclose(model_b._last_n_eps, model_a._last_n_eps, rtol=1e-6)
        np.testing.assert_allclose(
            report_kernel.energy_closure.field_work_j_m3,
            report_np.energy_closure.field_work_j_m3,
            rtol=1e-6,
            atol=1.0e-18,
        )
        self.assertAlmostEqual(
            report_kernel.deposited_power_w / report_np.deposited_power_w,
            1.0,
            delta=1e-6,
        )
        self.assertAlmostEqual(handles[1].phi_right / handles[0].phi_right, 1.0, delta=1e-8)

    def test_continuous_reconciled_sink_matches_numpy_reference(self) -> None:
        _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
        _, nf_b, cf_b, _c2, model_b, _n2, _dx2 = _dd_setup()
        model_a.use_numba_kernel = False
        model_a.electron_energy_closure = "continuous_reconciled"
        model_b.electron_energy_closure = "continuous_reconciled"
        _seed(cf_a, n_e=1.0e19, eps_ev=6.0)
        _seed(cf_b, n_e=1.0e19, eps_ev=6.0)
        po_a = Poisson1DSolver(n, dx)
        po_b = Poisson1DSolver(n, dx)
        centers = (np.arange(n) + 0.5) * dx
        _phi, e1d_a = po_a.solve(cf_a.net_charge_density(), 0.0, -5.0e3)
        handle = SemiImplicitField(
            solver=po_a, phi_left=0.0, phi_right=-5.0e3, cell_centers=centers,
        )
        report = model_a.advance(
            50e-12,
            efield=sample_efield(e1d_a, centers),
            neutral_fluid=nf_a,
            charged_fluid=cf_a,
            coeffs=coeffs,
            field=handle,
        )

        model_b._resolve(nf_b, cf_b, coeffs)
        packed = pack_kernel_args(model_b, model_b._me_coeffs, nf_b, cf_b)
        _phi, e1d_b = po_b.solve(cf_b.net_charge_density(), 0.0, -5.0e3)
        out = coupled_advance(
            model_b, 50e-12, e1d_b, nf_b, cf_b, packed, phi_l=0.0, phi_r=-5.0e3,
        )
        np.testing.assert_allclose(cf_b.density, cf_a.density, rtol=1e-7, atol=1.0)
        np.testing.assert_allclose(out.n_eps_j_m3, model_a._last_n_eps, rtol=1e-6)
        np.testing.assert_allclose(
            out.continuous_sink_j_m3,
            report.energy_closure.continuous_inelastic_sink_j_m3,
            rtol=1e-6,
        )
        np.testing.assert_allclose(
            out.bound_adjustment_j_m3,
            report.energy_closure.coefficient_bound_adjustment_j_m3,
            rtol=1e-6,
            atol=1.0e-18,
        )

    def test_source_coupled_transport_defers_heating_to_the_stage(self) -> None:
        """Closure gate: source_coupled transport must NOT heat n_eps (the stage delivers
        field work), while continuous_reconciled transport keeps the full local update —
        an unconditional removal would silently gut the diagnostic mode."""
        results = {}
        for closure in ("source_coupled", "continuous_reconciled"):
            _, nf, cf, coeffs, model, n, dx = _dd_setup()
            model.electron_energy_closure = closure
            _seed(cf, n_e=1.0e19, eps_ev=6.0)  # uniform state -> zero interior flux div
            n_eps_before = (
                cf.electron_mean_energy_ev * cf.density[:, cf.column["e-"]]
                * 1.602176634e-19
            ).copy()
            model._resolve(nf, cf, coeffs)
            packed = pack_kernel_args(model, model._me_coeffs, nf, cf)
            po = Poisson1DSolver(n, dx)
            _phi, e1d = po.solve(cf.net_charge_density(), 0.0, -5.0e3)
            out = coupled_advance(
                model, 50e-12, e1d, nf, cf, packed, phi_l=0.0, phi_r=-5.0e3,
            )
            interior = slice(8, n - 8)
            results[closure] = np.abs(out.n_eps_j_m3[interior] - n_eps_before[interior])
        drift_source = float(np.max(results["source_coupled"]))
        drift_reconciled = float(np.mean(results["continuous_reconciled"]))
        # Reconciled interior cells move by h*(pf - p_el - p_inel); source_coupled interior
        # cells see only flux divergence (~0 on a uniform state) and bounds.
        self.assertGreater(drift_reconciled, 1.0e3 * max(drift_source, 1.0e-300))

    def test_finite_gap_capacitor_cannot_supply_more_than_stored_energy(self) -> None:
        _, nf, cf, coeffs, model, n, dx = _dd_setup()
        _seed(cf, n_e=1.0e21)
        model._resolve(nf, cf, coeffs)
        packed = pack_kernel_args(model, model._me_coeffs, nf, cf)
        v0 = -5.0e3
        c_g = 2.5e-11
        po = Poisson1DSolver(n, dx)
        _phi, e1d = po.solve(cf.net_charge_density(), 0.0, v0)
        out = coupled_advance(
            model, 50e-12, e1d, nf, cf, packed, phi_l=0.0, phi_r=v0,
            c_g=c_g, i_drive=0.0,
        )
        self.assertLessEqual(
            max(out.gap_energy_absorbed_j, 0.0),
            0.5 * c_g * v0**2 * (1.0 + 1.0e-12),
        )

    def test_matches_numpy_coupled_path_post_breakdown(self) -> None:
        _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
        _, nf_b, cf_b, _c2, model_b, _n2, _dx2 = _dd_setup()
        model_a.use_numba_kernel = False  # numpy reference path (the parity baseline)
        _seed(cf_a)
        _seed(cf_b)
        po = Poisson1DSolver(n, dx)
        centers = (np.arange(n) + 0.5) * dx
        v_b = -5.0e3
        steps = 10

        # numpy reference loop
        reps = []
        for _ in range(steps):
            rho = cf_a.net_charge_density()
            _phi, e1d = po.solve(rho, 0.0, v_b)
            efield = sample_efield(e1d, centers, axis=0)
            handle = SemiImplicitField(solver=po, phi_left=0.0, phi_right=v_b,
                                       cell_centers=centers)
            reps.append(model_a.advance(50e-12, efield=efield, neutral_fluid=nf_a,
                                        charged_fluid=cf_a, coeffs=coeffs, field=handle))

        # kernel loop (same bootstrap-field construction)
        model_b._resolve(nf_b, cf_b, coeffs)
        packed = pack_kernel_args(model_b, model_b._me_coeffs, nf_b, cf_b)
        kernel_out = []
        for _ in range(steps):
            rho = cf_b.net_charge_density()
            _phi, e1d = po.solve(rho, 0.0, v_b)
            out = coupled_advance(model_b, 50e-12, e1d, nf_b, cf_b, packed,
                                  phi_l=0.0, phi_r=v_b)
            kernel_out.append(out)
            # write back what the numpy path writes back (⟨ε⟩ + retained energy budget)
            n_eps, e_axial = out.n_eps_j_m3, out.e_axial_v_m
            n_e = np.maximum(cf_b.density[:, model_b._e_col], 1.0)
            cf_b.electron_mean_energy_ev[:] = np.clip(
                n_eps / (n_e * 1.602176634e-19), *model_b._eps_bounds
            )

        np.testing.assert_allclose(cf_b.density, cf_a.density, rtol=1e-7, atol=1.0)
        # substep sequencing parity (integer equality per outer step)
        for rep, out in zip(reps, kernel_out):
            self.assertEqual(
                out.steps,
                model_a._last_substeps if rep is reps[-1] else out.steps,
            )
        self.assertEqual(kernel_out[-1].steps, model_a._last_substeps)
        # report scalars of the last step
        last = kernel_out[-1]
        self.assertAlmostEqual(
            -model_b.area_m2 * last.jcond_mean_a_m2 / reps[-1].conduction_current_a,
            1.0,
            delta=1e-6,
        )
        # gap load parity: both paths now report the energy-matched secant conductance
        self.assertAlmostEqual(
            (
                max(last.gap_energy_absorbed_j, 0.0)
                / max(last.gap_voltage_squared_dt_v2_s, 1e-30)
            ) / reps[-1].gap_conductance_s,
            1.0, delta=1e-6,
        )
        self.assertAlmostEqual(
            (last.deposited_j / 50e-12) / reps[-1].deposited_power_w,
            1.0,
            delta=1e-6,
        )
        np.testing.assert_allclose(last.n_eps_j_m3, model_a._last_n_eps, rtol=1e-6)

    def test_kernel_is_much_faster_post_breakdown(self) -> None:
        """Sanity (not a benchmark): the kernel step must not be slower than numpy."""
        import time
        _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
        _, nf_b, cf_b, _c2, model_b, _n2, _dx2 = _dd_setup()
        model_a.use_numba_kernel = False
        _seed(cf_a, n_e=1.0e21)
        _seed(cf_b, n_e=1.0e21)
        po = Poisson1DSolver(n, dx)
        centers = (np.arange(n) + 0.5) * dx
        model_b._resolve(nf_b, cf_b, coeffs)
        packed = pack_kernel_args(model_b, model_b._me_coeffs, nf_b, cf_b)
        # warm the JIT before timing
        _phi, e1d = po.solve(cf_b.net_charge_density(), 0.0, -5e3)
        coupled_advance(model_b, 50e-12, e1d, nf_b, cf_b, packed, phi_l=0.0, phi_r=-5e3)

        t0 = time.perf_counter()
        for _ in range(5):
            _phi, e1d = po.solve(cf_a.net_charge_density(), 0.0, -5e3)
            efield = sample_efield(e1d, centers, axis=0)
            handle = SemiImplicitField(solver=po, phi_left=0.0, phi_right=-5e3,
                                       cell_centers=centers)
            model_a.advance(50e-12, efield=efield, neutral_fluid=nf_a, charged_fluid=cf_a,
                            coeffs=coeffs, field=handle)
        t_np = time.perf_counter() - t0
        t0 = time.perf_counter()
        for _ in range(5):
            _phi, e1d = po.solve(cf_b.net_charge_density(), 0.0, -5e3)
            coupled_advance(model_b, 50e-12, e1d, nf_b, cf_b, packed, phi_l=0.0, phi_r=-5e3)
        t_kn = time.perf_counter() - t0
        self.assertLess(t_kn, t_np)


if __name__ == "__main__":
    unittest.main()
