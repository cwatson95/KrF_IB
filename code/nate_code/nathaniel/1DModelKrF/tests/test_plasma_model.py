"""Charged-species plasma model interface + stub (fixed-Te, Joule-heating power)."""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D.boltzmann import (
    AnalyticStubSolver,
    build_lut,
    en_grid_from_spec,
    reduced_field_td,
)
from Fluid1D.charged_fluid import ChargedFluid
from Fluid1D.plasma_model import (
    DriftDiffusionPlasmaModel,
    StepReport,
    StubPlasmaModel,
)
from kinetics.fluid import NeutralFluid
from particles import excimer_system_registry
from reactions import load_ion_mobility, load_reaction_network

_KRF = _ROOT / "data" / "reactions" / "krf.toml"
_ION = _ROOT / "data" / "ion_transport" / "ion_mobility.toml"


def _fluids(n=8, dx=1.0e-3):
    reg = excimer_system_registry()
    neutrals = tuple(s for s in reg if not s.is_charged and not s.is_photon)
    charged = tuple(s for s in reg if s.is_charged)
    nf = NeutralFluid(neutrals, n, dx, area_m2=1.0e-4)
    nf.seed_uniform("Kr", 2.4e25)
    cf = ChargedFluid(charged, n, dx, area_m2=1.0e-4)
    cf.seed_uniform("e-", 1.0e18)
    cf.seed_uniform("Kr+", 1.0e18)
    return nf, cf


class TestStubPlasmaModel(unittest.TestCase):
    def setUp(self) -> None:
        self.nf, self.cf = _fluids()
        self.coeffs = AnalyticStubSolver().solve({"Kr": 1.0}, en_grid_from_spec((1.0, 1000.0, 32)))
        self.efield = np.zeros((self.cf.n_cells, 3))
        self.efield[:, 0] = 1.0e6  # V/m

    def test_pins_te_from_lut_and_reports_power(self) -> None:
        model = StubPlasmaModel()
        report = model.advance(1e-12, efield=self.efield, neutral_fluid=self.nf,
                               charged_fluid=self.cf, coeffs=self.coeffs)
        self.assertIsInstance(report, StepReport)
        self.assertGreater(report.mean_electron_energy_ev, 0.0)
        self.assertTrue(np.isfinite(report.deposited_power_w))
        self.assertGreaterEqual(report.deposited_power_w, 0.0)
        # <eps> pinned to 1.5 * Te (Te from the LUT at this |E|/N)
        np.testing.assert_allclose(self.cf.electron_mean_energy(),
                                   self.cf.electron_mean_energy()[0])

    def test_fixed_te_override(self) -> None:
        model = StubPlasmaModel(fixed_te_ev=5.0)
        model.advance(1e-12, efield=self.efield, neutral_fluid=self.nf,
                      charged_fluid=self.cf, coeffs=self.coeffs)
        np.testing.assert_allclose(self.cf.electron_mean_energy(), 7.5)  # 1.5 * 5.0

    def test_does_not_grow_charged_density(self) -> None:
        before = self.cf.density.copy()
        StubPlasmaModel().advance(1e-12, efield=self.efield, neutral_fluid=self.nf,
                                  charged_fluid=self.cf, coeffs=self.coeffs)
        np.testing.assert_array_equal(self.cf.density, before)  # stub holds densities


def _dd_setup(n=48, dx=None):
    reg = excimer_system_registry()
    L = 0.012
    dx = dx or L / n
    area = 9.0e-3
    neutrals = tuple(s for s in reg if not s.is_charged and not s.is_photon)
    charged = tuple(s for s in reg if s.is_charged)
    nf = NeutralFluid(neutrals, n, dx, area_m2=area)
    nf.seed_uniform("Kr", 8.3e25 * 0.999)
    nf.seed_uniform("F2", 8.3e25 * 0.001)
    cf = ChargedFluid(charged, n, dx, area_m2=area)
    net = load_reaction_network(_KRF, reg).compile()
    coeffs = build_lut({"Kr": 0.999, "F2": 0.001}, solver=AnalyticStubSolver(),
                       en_grid_td=en_grid_from_spec((1.0, 1000.0, 32)), reactions=net,
                       key_kind="mean_energy_eV")
    model = DriftDiffusionPlasmaModel(
        reactions=net, ion_mobility=load_ion_mobility(_ION), mixture={"Kr": 0.999, "F2": 0.001},
        gas_temp_k=300.0, area_m2=area, gap_m=L,
    )
    return reg, nf, cf, coeffs, model, n, dx


class TestDriftDiffusionPlasmaModel(unittest.TestCase):
    def test_runs_stable_finite(self) -> None:
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        x = (np.arange(n) + 0.5) * cf.dx
        cf.seed_profile("e-", 1e14 * np.tanh(100 * x) + 1e13)
        cf.seed_profile("Kr+", 1e14 * np.tanh(100 * x) + 1e13)
        cf.seed_electron_energy(2.0)
        e = np.zeros((n, 3)); e[:, 0] = 1.0e6
        rep = model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        self.assertIsInstance(rep, StepReport)
        self.assertTrue(np.all(np.isfinite(cf.density)) and np.all(cf.density >= 0.0))
        self.assertGreaterEqual(rep.deposited_power_w, 0.0)   # Joule heating >= 0
        self.assertNotEqual(rep.conduction_current_a, 0.0)     # real flux integral fed back

    def test_electrons_drift_toward_anode(self) -> None:
        # E points +x (anode at x=0); electrons (z=-1) drift toward x=0.
        _, nf, cf, coeffs, model, n, dx = _dd_setup()
        x = (np.arange(n) + 0.5) * dx
        blob = 1e16 * np.exp(-((x - 0.006) / 0.0015) ** 2)
        cf.seed_profile("e-", blob.copy()); cf.seed_profile("Kr+", blob.copy())
        cf.seed_electron_energy(2.0)
        e = np.zeros((n, 3)); e[:, 0] = 5.0e5
        com0 = float(np.sum(x * cf.density_of("e-")) / np.sum(cf.density_of("e-")))
        for i in range(200):
            model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        com1 = float(np.sum(x * cf.density_of("e-")) / np.sum(cf.density_of("e-")))
        self.assertLess(com1, com0)  # electron cloud moved toward the anode

    def test_conduction_current_loads_the_gap_dissipatively(self) -> None:
        """The 1D flux integral must agree in sign with the 0D closure and dissipate power."""
        from Fluid1D.field import conduction_current_0d
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        cf.seed_uniform("e-", 1e16); cf.seed_uniform("Kr+", 1e16)
        cf.seed_electron_energy(3.0)
        v_gap = -15.0e3                     # cathode negative -> E_x > 0 on our mesh
        e = np.zeros((n, 3)); e[:, 0] = 1.0e6
        rep = model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        i_1d = rep.conduction_current_a
        i_0d = conduction_current_0d(1e16, 0.012, model.area_m2, model.gap_m, v_gap)
        self.assertEqual(np.sign(i_1d), np.sign(i_0d))   # same convention as the 0D closure
        self.assertGreaterEqual(i_1d * v_gap, 0.0)       # dissipative: I and V_gap share a sign

    def test_local_field_closure_rescues_emptied_cells(self) -> None:
        """Where attachment empties a cell, ⟨ε⟩ follows the local field, not 0/0."""
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        # Half the domain has electrons, half is emptied to a pathological trace density.
        n_e = np.full(n, 1.0e16); n_e[n // 2:] = 1.0e-12
        cf.seed_profile("e-", n_e); cf.seed_uniform("Kr+", 1.0e16)
        cf.seed_electron_energy(3.0)
        e = np.zeros((n, 3)); e[:, 0] = 1.0e6
        model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        eps = cf.electron_mean_energy()
        self.assertTrue(np.all(np.isfinite(eps)))
        # The emptied interior must sit at the local-field equilibrium, not pinned to a clamp.
        # The last cell is excluded: cathode secondary emission reseeds it with electrons, so
        # LMEA legitimately applies there.
        lo, hi = model._eps_bounds
        interior = slice(n // 2 + 2, n - 1)
        emptied = eps[interior]
        self.assertTrue(np.all(emptied > lo * 1.01), f"pinned to the floor: {emptied[:4]}")
        self.assertTrue(np.all(emptied < hi * 0.99), f"pinned to the ceiling: {emptied[:4]}")
        n_gas = nf.density.sum(axis=1)
        eps_eq = coeffs.equilibrium_mean_energy(reduced_field_td(np.abs(e[:, 0]), n_gas))
        np.testing.assert_allclose(emptied, eps_eq[interior], rtol=1e-6)

    def test_mean_energy_never_leaves_the_table_range(self) -> None:
        """An out-of-range ⟨ε⟩ would silently clamp inside every accessor — never allow it."""
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        cf.seed_uniform("e-", 1e16); cf.seed_uniform("Kr+", 1e16)
        cf.seed_electron_energy(3.0)
        e = np.zeros((n, 3)); e[:, 0] = 5.0e6  # deliberately violent field
        for _ in range(20):
            model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        lo, hi = coeffs.mean_energy_bounds_ev
        eps = cf.electron_mean_energy()
        self.assertTrue(np.all(eps >= lo - 1e-12) and np.all(eps <= hi + 1e-12))

    def test_energy_pde_bounds_and_evolves(self) -> None:
        _, nf, cf, coeffs, model, n, _ = _dd_setup()
        cf.seed_uniform("e-", 1e16); cf.seed_uniform("Kr+", 1e16)
        cf.seed_electron_energy(2.0)
        e = np.zeros((n, 3)); e[:, 0] = 1.0e6
        rep = model.advance(5e-11, efield=e, neutral_fluid=nf, charged_fluid=cf, coeffs=coeffs)
        eps = cf.electron_mean_energy()
        self.assertTrue(np.all(eps >= model.eps_floor - 1e-12))
        self.assertTrue(np.all(eps <= model.eps_ceiling + 1e-9))
        self.assertGreater(rep.mean_electron_energy_ev, 0.0)


class TestSemiImplicitFieldCoupling(unittest.TestCase):
    """The coupled (σ-corrected, per-substep) field solve: stability + reference parity.

    The frozen-field staggered coupling is linearly unstable once σ > ~2ε₀/dt
    (nₑ ≳ 3e19 m⁻³ at dt = 50 ps): with a −5 kV Dirichlet gap and nₑ = 1e20 it drives
    |E| → ~1e10 V/m within 1 ns while manufacturing energy. The semi-implicit handle must keep
    the same scenario bounded.
    """

    @staticmethod
    def _coupled_loop(nf, cf, coeffs, model, n, dx, *, steps, v_b=-5.0e3, dt=50e-12):
        from Fluid1D.field import Poisson1DSolver, SemiImplicitField, sample_efield
        po = Poisson1DSolver(n, dx)
        centers = (np.arange(n) + 0.5) * dx
        rep, e_max_hist = None, []
        for _ in range(steps):
            rho = cf.net_charge_density()
            _phi, e1d = po.solve(rho, 0.0, v_b)                # bootstrap/sizing field
            efield = sample_efield(e1d, centers, axis=0)
            handle = SemiImplicitField(solver=po, phi_left=0.0, phi_right=v_b,
                                       cell_centers=centers)
            rep = model.advance(dt, efield=efield, neutral_fluid=nf, charged_fluid=cf,
                                coeffs=coeffs, field=handle)
            e_max_hist.append(float(np.max(np.abs(handle.last_efield[:, 0]))))
        return rep, handle, np.asarray(e_max_hist)

    def test_semi_implicit_survives_post_breakdown_density(self) -> None:
        """nₑ = 1e20 for 5 ns at dt = 50 ps: bounded fields, finite state, sane energy.

        The frozen-field coupling passes 1e8 V/m within ~10 steps here and reaches ~1e10
        within 1 ns; the coupled path instead builds a cathode fall that plateaus around
        7e7 V/m and then relaxes — assert boundedness *and* the non-growing trend.
        """
        _, nf, cf, coeffs, model, n, dx = _dd_setup()
        cf.seed_uniform("e-", 1.0e20); cf.seed_uniform("Kr+", 1.0e20)
        cf.seed_electron_energy(3.0)
        rep, handle, e_hist = self._coupled_loop(nf, cf, coeffs, model, n, dx, steps=250)
        self.assertTrue(np.all(np.isfinite(cf.density)) and np.all(cf.density >= 0.0))
        self.assertTrue(np.all(np.isfinite(cf.electron_mean_energy())))
        # Cell-centred field (2026-07-31): the wall half-cell face reads the cathode-fall
        # gradient ~2x steeper than the old node scheme's cell average, and the fall forms
        # over ~150 steps instead of ~40 (measured plateau 1.77e8 at steps 150-200, then
        # relaxing to 1.39e8 by step 400) — bounds re-pinned accordingly.
        self.assertLess(float(e_hist.max()), 5.0e8)  # sheath-scale, well below blowup (1e10+)
        # Bounded orbit: once the sheath has formed (~step 150) the peak field stops growing.
        self.assertLess(float(e_hist[-40:].max()), 1.2 * float(e_hist[120:160].max()))
        self.assertLess(float(np.mean(cf.electron_mean_energy())), 10.0)
        self.assertTrue(np.isfinite(rep.deposited_power_w))
        self.assertTrue(np.isfinite(rep.conduction_current_a))

    def test_coupled_path_matches_frozen_field_prebreakdown(self) -> None:
        """At nₑ = 1e14 (h·σ/ε₀ ~ 1e-5) the coupled and staggered loops must agree."""
        from Fluid1D.field import Poisson1DSolver, sample_efield
        _, nf_a, cf_a, coeffs, model_a, n, dx = _dd_setup()
        _, nf_b, cf_b, _cb, model_b, _n2, _dx2 = _dd_setup()
        for cf in (cf_a, cf_b):
            cf.seed_uniform("e-", 1.0e14); cf.seed_uniform("Kr+", 1.0e14)
            cf.seed_electron_energy(3.0)
        # coupled loop
        self._coupled_loop(nf_a, cf_a, coeffs, model_a, n, dx, steps=10)
        # staggered loop (frozen field per outer step)
        po = Poisson1DSolver(n, dx)
        centers = (np.arange(n) + 0.5) * dx
        for _ in range(10):
            rho = cf_b.net_charge_density()
            _phi, e1d = po.solve(rho, 0.0, -5.0e3)
            efield = sample_efield(e1d, centers, axis=0)
            model_b.advance(50e-12, efield=efield, neutral_fluid=nf_b, charged_fluid=cf_b,
                            coeffs=coeffs)
        np.testing.assert_allclose(cf_a.density, cf_b.density, rtol=1e-3, atol=1.0)
        np.testing.assert_allclose(cf_a.electron_mean_energy(),
                                   cf_b.electron_mean_energy(), rtol=1e-3)


if __name__ == "__main__":
    unittest.main()
