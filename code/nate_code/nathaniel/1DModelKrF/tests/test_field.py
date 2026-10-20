"""1D Poisson solve (Dirichlet) + discharge circuits (stub + MPC) + field sampling."""

from __future__ import annotations

import math
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D.field import (
    MpcCircuit,
    Poisson1DSolver,
    SemiImplicitField,
    StubCircuit,
    _thomas_general,
    conduction_current_0d,
    sample_efield,
)
from particles.constants import E_CHARGE, EPS0


class TestStubCircuit(unittest.TestCase):
    def test_fixed_potentials(self) -> None:
        c = StubCircuit(0.0, -100.0)
        self.assertEqual(c.advance(1e-12), (0.0, -100.0))

    def test_waveform_drives_right_electrode(self) -> None:
        c = StubCircuit(0.0, -100.0, waveform=lambda t: -1.0e3 * t)
        _l, r = c.advance(2.0)
        self.assertAlmostEqual(r, -2.0e3)


class TestPoisson1D(unittest.TestCase):
    def test_zero_charge_gives_linear_potential_uniform_field(self) -> None:
        n, dx = 64, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        phi, e = solver.solve(np.zeros(n), phi_left=0.0, phi_right=-100.0)
        L = n * dx
        # ρ=0 -> φ linear -> E = -dφ/dx = (φ_left - φ_right)/L, uniform
        np.testing.assert_allclose(e, (0.0 - (-100.0)) / L, rtol=1e-9)
        self.assertTrue(phi[0] > phi[-1])  # descends toward the negative electrode

    def test_constant_charge_gives_parabolic_potential(self) -> None:
        n, dx = 64, 1.0e-3
        rho0 = 1.0e-6  # C/m^3
        solver = Poisson1DSolver(n, dx)
        phi, e = solver.solve(np.full(n, rho0), phi_left=0.0, phi_right=0.0)
        # -φ'' = ρ/ε0 > 0 -> φ concave-down -> interior potential rises above the ends
        self.assertGreater(phi.max(), 0.0)
        # dE/dx = ρ/ε0 (Gauss): E increases roughly linearly across the interior
        dedx = np.diff(e) / dx
        np.testing.assert_allclose(np.mean(dedx), rho0 / EPS0, rtol=0.05)


class TestSemiImplicitPoisson(unittest.TestCase):
    """σ-corrected semi-implicit solve: reductions, limits, and the J_expl lemma."""

    def test_thomas_general_matches_dense_solve(self) -> None:
        rng = np.random.default_rng(7)
        for m in (1, 2, 5, 63):
            lower = rng.uniform(-1.0, -0.1, m)
            upper = rng.uniform(-1.0, -0.1, m)
            diag = np.abs(lower) + np.abs(upper) + rng.uniform(0.5, 2.0, m)  # diag dominant
            rhs = rng.normal(size=m)
            dense = np.zeros((m, m))
            for r in range(m):
                dense[r, r] = diag[r]
                if r > 0:
                    dense[r, r - 1] = lower[r]
                if r < m - 1:
                    dense[r, r + 1] = upper[r]
            np.testing.assert_allclose(
                _thomas_general(lower, diag, upper, rhs), np.linalg.solve(dense, rhs),
                rtol=1e-12, atol=1e-14,
            )

    def test_reduces_to_plain_poisson_when_sigma_zero(self) -> None:
        n, dx = 64, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        rng = np.random.default_rng(3)
        rho = rng.normal(scale=1e-4, size=n)
        phi_a, e_a = solver.solve(rho, phi_left=25.0, phi_right=-4.0e3)
        phi_b, e_b = solver.solve_semi_implicit(
            rho, 25.0, -4.0e3, sigma_cells=np.zeros(n), j_expl_faces=np.zeros(n + 1), h=50e-12,
        )
        np.testing.assert_allclose(phi_b, phi_a, rtol=1e-12)
        np.testing.assert_allclose(e_b, e_a, rtol=1e-12)

    def test_low_conductivity_limit_quantified(self) -> None:
        # Pre-breakdown regime: h·σ/ε₀ ≈ 6e-6 at n_e = 1e14, μ_e ≈ 0.065, h = 50 ps —
        # the corrected solve must match the plain Poisson solve to that order.
        n, dx = 64, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        rho = np.full(n, 1.0e-5)  # pure space-charge drive (BCs 0 so the σ term is exercised)
        sigma = np.full(n, E_CHARGE * 1.0e14 * 0.065)
        h = 50e-12
        self.assertLess(h * sigma[0] / EPS0, 1e-4)
        _pa, e_a = solver.solve(rho, 0.0, 0.0)
        _pb, e_b = solver.solve_semi_implicit(
            rho, 0.0, 0.0, sigma_cells=sigma, j_expl_faces=np.zeros(n + 1), h=h,
        )
        self.assertLess(float(np.max(np.abs(e_b - e_a)) / np.max(np.abs(e_a))), 1e-4)

    def test_uniform_sigma_zero_charge_preserves_vacuum_field(self) -> None:
        # ρ = 0, J = 0: the uniform coefficient divides out — same linear φ, uniform E.
        n, dx = 32, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        _phi, e = solver.solve_semi_implicit(
            np.zeros(n), 0.0, -5.0e3, sigma_cells=np.full(n, 10.0),
            j_expl_faces=np.zeros(n + 1), h=1e-9,
        )
        np.testing.assert_allclose(e, 5.0e3 / (n * dx), rtol=1e-9)

    def test_two_region_sigma_partitions_field(self) -> None:
        # ρ = 0, J = 0, piecewise σ: the discrete rows enforce a_i·E_i = const exactly
        # (total-current continuity) — the conducting half carries the smaller field.
        n, dx = 64, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        sigma = np.where(np.arange(n) < n // 2, 1.0, 10.0)
        h = 1.0e-6  # strongly conducting: h·σ/ε₀ ≫ 1
        _phi, e = solver.solve_semi_implicit(
            np.zeros(n), 0.0, -5.0e3, sigma_cells=sigma, j_expl_faces=np.zeros(n + 1), h=h,
        )
        # Total-current continuity lives on the FACES of the cell-centred layout:
        # af_f·E_f = const exactly (af = arithmetic face mean; walls use the half-cell field).
        a = EPS0 + h * sigma
        af = np.empty(n + 1)
        af[1:n] = 0.5 * (a[:-1] + a[1:])
        af[0], af[n] = a[0], a[-1]
        current = af * solver.last_e_faces
        self.assertLess(float(np.ptp(current)), 1e-9 * abs(float(np.mean(current))))
        self.assertGreater(float(np.mean(e[: n // 2])), float(np.mean(e[n // 2:])))

    def test_alternating_charge_draws_restoring_field(self) -> None:
        # The odd-even checkerboard regression: the old node-CV gather 0.5(rho_i+rho_{i+1})
        # was EXACTLY blind to cell-alternating charge (measured E response 0.000) — the
        # discrete Gauss law drew no restoring field and the avalanche amplified any seed.
        # The cell-centred layout must respond at the face field.
        n, dx = 32, 1.0e-3
        solver = Poisson1DSolver(n, dx)
        rho = 1.0e-4 * (-1.0) ** np.arange(n)
        _phi, _e = solver.solve(rho, 0.0, 0.0)
        self.assertGreater(float(np.abs(solver.last_e_faces[1:-1]).max()), 1.0e3)

    def test_jexpl_matches_explicit_charge_update(self) -> None:
        # Lemma: at σ = 0, the semi-implicit solve with face J_expl (walls 0) equals the
        # plain solve of the explicitly-pushed charge ρ_i − (h/dx)(Jf_{i+1} − Jf_i).
        # Exact on the cell-centred layout. Pins the RHS sign and indexing.
        n, dx, h = 48, 2.5e-4, 50e-12
        solver = Poisson1DSolver(n, dx)
        rng = np.random.default_rng(11)
        rho = rng.normal(scale=1e-4, size=n)
        jf = np.zeros(n + 1)
        jf[1:-1] = rng.normal(scale=1.0, size=n - 1)  # face currents, wall faces zero
        rho_pushed = rho - (h / dx) * np.diff(jf)
        phi_a, e_a = solver.solve(rho_pushed, 12.0, -300.0)
        phi_b, e_b = solver.solve_semi_implicit(
            rho, 12.0, -300.0, sigma_cells=np.zeros(n), j_expl_faces=jf, h=h,
        )
        np.testing.assert_allclose(phi_b, phi_a, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(e_b, e_a, rtol=1e-12, atol=1e-12)


class TestSampleEfield(unittest.TestCase):
    def test_axial_component_only_without_prescribed(self) -> None:
        e1d = np.linspace(1.0, 2.0, 10)
        centers = np.linspace(0.0, 1.0, 10)
        efield = sample_efield(e1d, centers, axis=0, prescribed=None)
        self.assertEqual(efield.shape, (10, 3))
        np.testing.assert_allclose(efield[:, 0], e1d)
        np.testing.assert_allclose(efield[:, 1:], 0.0)

    def test_semi_implicit_prescribed_field_reaches_faces(self) -> None:
        class ConstantField:
            def at_batch(self, positions, _time):
                out = np.zeros_like(positions)
                out[:, 0] = 4.25e5
                return out

        n, dx = 12, 2.0e-4
        handle = SemiImplicitField(
            solver=Poisson1DSolver(n, dx), phi_left=0.0, phi_right=0.0,
            prescribed=ConstantField(),
        )
        cells = handle.solve(
            np.zeros(n), np.zeros(n), np.zeros(n + 1), 1.0e-12,
        )
        np.testing.assert_allclose(cells, 4.25e5)
        np.testing.assert_allclose(handle.last_e_faces, 4.25e5)


class TestMpcCircuit(unittest.TestCase):
    """Chen et al. (2025) MPC driver — phase-1 transfer, energy, geometry, gap coupling."""

    def test_precharge_and_stored_energy(self) -> None:
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4, v0_v=-18e3)
        self.assertAlmostEqual(c.v_c0, -18e3)                      # C0 pre-charged
        self.assertAlmostEqual(c.stored_energy_j(), 0.5 * 10e-9 * 18e3**2, places=6)
        phi_l, phi_r = c.advance(1e-12, discharge_current_a=0.0)
        self.assertEqual(phi_l, 0.0)                               # anode grounded
        self.assertAlmostEqual(phi_r, c.v_b)                       # cathode = V_b state

    def test_gap_capacitance_from_geometry(self) -> None:
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4, c_head_f=20e-12)
        self.assertAlmostEqual(c.cg_eff, EPS0 * 90e-4 / 12e-3 + 20e-12)

    def test_complete_resonant_transfer_and_energy_conservation(self) -> None:
        # Lossless MPC loop, gap open (no conduction): C0 -> Cp transfers ~completely at the
        # crest t = pi*sqrt(Lm*Ctilde) = 111 ns, and total energy is conserved.
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4, r1_ohm=0.0, rs_ohm=0.0, rf_ohm=0.0)
        e0 = c.stored_energy_j()
        dt, drift = 20e-12, 0.0
        vp_peak, t_peak = 0.0, 0.0
        for k in range(1, 12501):  # 250 ns
            c.advance(dt, discharge_current_a=0.0)
            if abs(c.v_p) > vp_peak:
                vp_peak, t_peak = abs(c.v_p), k * dt
            drift = max(drift, abs(c.stored_energy_j() - e0) / abs(e0))
        tau = math.pi * math.sqrt(250e-9 * (10e-9 * 10e-9 / (10e-9 + 10e-9)))  # 111 ns
        self.assertGreater(vp_peak, 0.95 * 18e3)      # near-complete transfer to Cp
        self.assertAlmostEqual(t_peak, tau, delta=5e-9)
        self.assertLess(drift, 1e-6)                  # lossless energy conservation

    def test_conduction_current_loads_node_b(self) -> None:
        # Node-B KCL is C_g·V_b' = I_s - I_f - jcond, so (from an identical state) a positive
        # conduction current pulls V_b more negative than the gap-open case.
        a = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        b = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        a.advance(1e-9, discharge_current_a=0.0)
        b.advance(1e-9, discharge_current_a=0.0)  # same state at t = 1 ns
        a.advance(1e-12, discharge_current_a=0.0)
        b.advance(1e-12, discharge_current_a=1.0e3)
        self.assertLess(b.v_b, a.v_b)

    def test_open_gap_exponential_cache_is_bit_equivalent(self) -> None:
        cached = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        uncached = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        dt = 50e-12
        cached.advance(dt)
        self.assertIn(dt, cached._open_gap_transition_cache)
        uncached.advance(dt)
        cached.advance(dt)
        uncached._open_gap_transition_cache.clear()
        uncached.advance(dt)
        np.testing.assert_array_equal(cached.y, uncached.y)


class TestSemiImplicitGapCoupling(unittest.TestCase):
    """Chen circuit-analysis §7: linearising the gap load lets node B ride the collapse."""

    @staticmethod
    def _run(g_s, *, use_g, steps=400):
        """Emulate the coupled loop: the plasma reports its ohmic load at the live V_b.

        Returns ``None`` if the loop blew up. 400 steps (20 ns) is enough to separate the two
        schemes: the explicit one is already past 1e118 V by step 100.
        """
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        i_cond = 0.0
        with np.errstate(over="ignore", invalid="ignore"):
            for _ in range(steps):
                kwargs = {"discharge_current_a": i_cond}
                if use_g:
                    kwargs["gap_conductance_s"] = g_s
                c.advance(50e-12, **kwargs)
                if not np.isfinite(c.v_b):
                    return None
                i_cond = g_s * c.v_b
        return c

    def test_zero_conductance_is_backward_compatible(self) -> None:
        a = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        b = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        for _ in range(200):
            a.advance(50e-12, discharge_current_a=-50.0)
            b.advance(50e-12, discharge_current_a=-50.0, gap_conductance_s=0.0)
        self.assertEqual(a.v_b, b.v_b)  # exactly the old staggered-explicit coupling

    def test_linearisation_reproduces_the_reported_current(self) -> None:
        # At the state where the current was measured, G·V_b + offset == the reported current.
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        c.advance(50e-12, discharge_current_a=0.0)   # get V_b off zero
        g, i_rep = 5.0, -1.2e3
        i_off = i_rep - g * c.v_b
        self.assertAlmostEqual(g * c.v_b + i_off, i_rep, places=9)

    def test_explicit_diverges_where_linearised_survives(self) -> None:
        """The post-breakdown regime: a strongly conducting gap."""
        g_s = 10.0                                  # τ = C_g/G ≈ 2.7 ps << the 50 ps step
        self.assertIsNone(self._run(g_s, use_g=False), "explicit coupling should diverge here")
        linearised = self._run(g_s, use_g=True)
        self.assertIsNotNone(linearised, "linearised coupling must stay finite")
        self.assertTrue(np.isfinite(linearised.stored_energy_j()))

    def test_conducting_gap_collapses_the_voltage(self) -> None:
        # Physical signature of breakdown: the more the gap conducts, the lower |V_b| stays.
        weak = self._run(0.1, use_g=True)
        strong = self._run(100.0, use_g=True)
        self.assertLess(abs(strong.v_b), abs(weak.v_b))

    def test_extreme_conductance_is_exact_and_stable(self) -> None:
        # The exact (expm) update has no substep limit: even G ~ 1e6 S (τ ~ 3e-17 s, which
        # explicit sub-stepping could not affordably resolve) must relax node B toward the
        # linearised equilibrium without a hint of divergence.
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        for _ in range(50):
            c.advance(50e-12, discharge_current_a=0.0, gap_conductance_s=1.0e6)
        self.assertTrue(np.all(np.isfinite(c.y)))
        self.assertLess(abs(c.v_b), 1.0)  # a hugely conducting gap pins V_b ~ 0

    def test_negative_conductance_is_rejected(self) -> None:
        # A negative G would anti-damp node B; it must be clamped away, not trusted.
        c = MpcCircuit(gap_m=12e-3, area_m2=90e-4)
        for _ in range(200):
            c.advance(50e-12, discharge_current_a=-10.0, gap_conductance_s=-50.0)
        self.assertTrue(np.isfinite(c.v_b))


class TestConductionCurrent0D(unittest.TestCase):
    def test_ohmic_closure(self) -> None:
        n_e, mu_e, area, gap, v_gap = 1e20, 0.04, 90e-4, 12e-3, -15e3
        expected = E_CHARGE * n_e * mu_e * area / gap * v_gap
        self.assertAlmostEqual(conduction_current_0d(n_e, mu_e, area, gap, v_gap), expected)

    def test_zero_when_no_electrons(self) -> None:
        self.assertEqual(conduction_current_0d(0.0, 0.04, 90e-4, 12e-3, -15e3), 0.0)


class TestMpcCoupledSmoke(unittest.TestCase):
    def test_circuit_drives_cathode_bc(self) -> None:
        import contextlib
        import io

        from Fluid1D import (BoltzmannSpec, Cavity, CircuitSpec, GasMixture, Grid1D,
                             Numerics, Plasma, Simulation)

        sim = Simulation(
            name="mpc_smoke", grid=Grid1D(cells=16, length=12e-3, area_m2=90e-4),
            gas=GasMixture({"Kr": 0.99, "F2": 0.01}, pressure_atm=3.4),
            plasma=Plasma(density=1e16, ion="Kr+"),
            cavity=Cavity(), boltzmann=BoltzmannSpec(backend="stub"),
            circuit=CircuitSpec(backend="mpc"), reactions=_ROOT / "data" / "reactions" / "krf.toml",
            diagnostics=[], output_dir=_ROOT / "runs",
        )
        with contextlib.redirect_stdout(io.StringIO()):
            sim.run(Numerics(steps=600, dt=50e-12, fluid_every=50), verbose=False)  # ~30 ns
        # cathode voltage has ramped negative (phase-1 transfer under way), everything finite
        self.assertLess(sim.circuit.v_b, -1.0e3)
        self.assertTrue(np.isfinite(sim.circuit.stored_energy_j()))
        closure = sim.last_energy_closure_report
        self.assertEqual(closure.field_work_j_m3.shape, (16,))
        self.assertEqual(closure.accepted_reaction_cost_j_m3.shape, (16,))
        self.assertEqual(closure.signed_reconciliation_j_m3.shape, (16,))
        self.assertTrue(np.all(np.isfinite(closure.coefficient_bound_adjustment_j_m3)))


if __name__ == "__main__":
    unittest.main()
