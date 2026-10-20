"""Two-term Boltzmann solver validation: exact limits, tail depletion, table contract.

Anchors chosen to be *provably* correct rather than remembered numbers: the thermal limit
(detailed balance with the gas), the Druyvesteyn limit (closed-form EEDF for constant σ,
elastic-only, cold gas), Maxwellian-moment closed forms (Einstein exact; energy-pair ratio
4/3 for constant σ), and the tail-depletion trend against the Maxwellian solver.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D.boltzmann import en_grid_from_spec
from Fluid1D.eedf_solver import MaxwellianEedfSolver, maxwellian_f0
from Fluid1D.two_term import TwoTermBoltzmannSolver
from particles import excimer_system_registry
from particles.constants import E_CHARGE, K_B, molecular_mass
from reactions import build_cross_section_index, load_reaction_network

_XS = _ROOT / "data" / "cross_sections"
_ARF = _ROOT / "data" / "reactions" / "arf.toml"


class _ConstSigma:
    """Duck-typed stand-in for TabulatedSigma: constant σ, arbitrary threshold."""

    def __init__(self, sigma_m2: float, threshold_ev: float = 0.0):
        self._s = float(sigma_m2)
        self.threshold_ev = float(threshold_ev)

    def __call__(self, eps):
        eps = np.asarray(eps, dtype=np.float64)
        return np.where(eps > self.threshold_ev, self._s, 0.0)


class TestMoments(unittest.TestCase):
    """Closed-form Maxwellian moments — pins the integral conventions exactly."""

    def setUp(self) -> None:
        self.te = 3.0
        self.sigma = 5.0e-20
        n, eps_max = 6000, 90.0
        self.d_eps = eps_max / n
        self.eps = (np.arange(n) + 0.5) * self.d_eps
        self.f0 = maxwellian_f0(self.eps, self.te)
        self.sig_c = np.full(n, self.sigma)

    def test_einstein_relation_exact(self) -> None:
        m = TwoTermBoltzmannSolver._moments(self.eps, self.d_eps, self.f0, self.sig_c, [])
        self.assertAlmostEqual(m["dif"] / m["mob"], self.te, delta=self.te * 2e-3)

    def test_diffusion_closed_form_constant_sigma(self) -> None:
        # DN = (γ/3σ)·(2/√π)·√Tₑ for a Maxwellian with constant σ.
        from Fluid1D.eedf_solver import GAMMA_V_PER_SQRT_EV as g
        m = TwoTermBoltzmannSolver._moments(self.eps, self.d_eps, self.f0, self.sig_c, [])
        ref = (g / (3.0 * self.sigma)) * (2.0 / np.sqrt(np.pi)) * np.sqrt(self.te)
        self.assertAlmostEqual(m["dif"] / ref, 1.0, delta=2e-3)

    def test_energy_pair_ratio_constant_sigma_is_4_3(self) -> None:
        # For constant σ the Maxwellian energy-transport ratio is exactly 4/3 (5/3 holds
        # for constant collision frequency — the table's default closure).
        m = TwoTermBoltzmannSolver._moments(self.eps, self.d_eps, self.f0, self.sig_c, [])
        self.assertAlmostEqual(m["dif_e"] / m["dif"], 4.0 / 3.0, delta=5e-3)
        self.assertAlmostEqual(m["mob_e"] / m["mob"], 4.0 / 3.0, delta=2e-2)

    def test_rate_closed_form_constant_sigma_with_threshold(self) -> None:
        # k = γσ·(2/√π)·Tₑ^{-3/2}·∫_u^∞ ε e^{-ε/Tₑ} dε = γσ(2/√π)√Tₑ·(1+u/Tₑ)e^{-u/Tₑ}.
        from Fluid1D.eedf_solver import GAMMA_V_PER_SQRT_EV as g
        u = 10.0
        procs = [("RX", 1.0, u, _ConstSigma(self.sigma, u), False, False)]
        m = TwoTermBoltzmannSolver._moments(self.eps, self.d_eps, self.f0, self.sig_c, procs)
        ref = g * self.sigma * (2.0 / np.sqrt(np.pi)) * np.sqrt(self.te) \
            * (1.0 + u / self.te) * np.exp(-u / self.te)
        self.assertAlmostEqual(m["rates"]["RX"] / ref, 1.0, delta=5e-3)


class TestExactLimits(unittest.TestCase):
    def test_thermal_limit_detailed_balance(self) -> None:
        """E/N → 0 in pure Kr: F₀ relaxes to the gas temperature (⟨ε⟩ → (3/2)kT_g)."""
        solver = TwoTermBoltzmannSolver(_XS, energy_max_floor_ev=1.0)
        tab = solver.solve({"Kr": 1.0}, np.array([1.0e-3]), gas_temp_k=300.0, reactions=None)
        kt = K_B * 300.0 / E_CHARGE
        self.assertAlmostEqual(tab.key_grid[0] / (1.5 * kt), 1.0, delta=0.15)
        # Einstein at thermal equilibrium: D/μ = kT_g/e.
        eps0 = tab.key_grid[0]
        self.assertAlmostEqual(
            float(tab.diffusion(eps0)) / float(tab.mobility(eps0)) / kt, 1.0, delta=0.1
        )

    def test_druyvesteyn_limit(self) -> None:
        """Constant σ, elastic-only, cold gas: ln F₀ is exactly linear in ε²."""
        sigma = 1.0e-19
        solver = TwoTermBoltzmannSolver(_XS)
        solver._elastic_sigma_for = lambda gas: _ConstSigma(sigma)  # type: ignore[method-assign]
        en_td = 10.0
        res = solver._solve_one({"Kr": 1.0}, en_td, gas_temp_k=1.0, processes=[])
        f0, eps = res["f0"], res["eps"]
        m_over_M = 9.109e-31 / molecular_mass("Kr")
        slope_ref = -3.0 * m_over_M * sigma**2 / (en_td * 1.0e-21) ** 2
        sel = f0 > f0.max() * 1e-6
        slope = np.polyfit(eps[sel] ** 2, np.log(f0[sel]), 1)[0]
        self.assertAlmostEqual(slope / slope_ref, 1.0, delta=0.03)

    def test_normalisation_and_positivity(self) -> None:
        solver = TwoTermBoltzmannSolver(_XS)
        res = solver._solve_one({"Kr": 1.0}, 50.0, gas_temp_k=300.0, processes=[])
        f0, eps = res["f0"], res["eps"]
        d_eps = eps[1] - eps[0]
        self.assertTrue(np.all(f0 >= 0.0))
        self.assertAlmostEqual(float(np.sum(f0 * np.sqrt(eps)) * d_eps), 1.0, delta=1e-9)


class TestChenMixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        reg = excimer_system_registry()
        cls.net = load_reaction_network(_ARF, reg).compile()
        f2 = 0.001
        cls.mix = {"Ne": 1.0 - 0.035 - 0.0001 - f2, "Ar": 0.035, "F2": f2, "Xe": 0.0001}
        loss = tuple(p for p in build_cross_section_index(_XS) if "_vib_" in p)
        grid = en_grid_from_spec((0.5, 300.0, 24))
        cls.tab2 = TwoTermBoltzmannSolver(_XS, extra_loss_process_ids=loss).solve(
            cls.mix, grid, 300.0, reactions=cls.net
        )
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.tabM = MaxwellianEedfSolver(_XS, extra_loss_process_ids=loss).solve(
                cls.mix, grid, 300.0, reactions=cls.net
            )

    def test_table_contract(self) -> None:
        t = self.tab2
        self.assertEqual(t.key_kind, "mean_energy_eV")
        self.assertTrue(np.all(np.diff(t.key_grid) > 0.0))
        self.assertIsNotNone(t.reduced_field_td)
        for rid in ("R1", "R3", "R7", "R10"):
            self.assertIn(rid, t.rate_coeff)
        # The low-E/N extension must reach well below the Ne-buffer elastic plateau (~2.4 eV
        # at 0.1 Td), so cooled post-collapse cells aren't clamped onto hot coefficients.
        self.assertLess(float(t.key_grid[0]), 1.0)

    def test_real_energy_pair_not_the_5_3_closure(self) -> None:
        # The solver supplies the true ε-weighted pair; ratio sits between the constant-σ
        # (4/3) and constant-ν (5/3) Maxwellian limits, and differs from the 5/3 closure.
        eps = 5.0
        ratio = float(self.tab2.energy_mobility(eps)) / float(self.tab2.mobility(eps))
        self.assertGreater(ratio, 1.15)
        self.assertLess(ratio, 5.0 / 3.0)

    def test_tail_depletion_cuts_ionization_at_fixed_mean_energy(self) -> None:
        """The physical point of the solver: k_ion(⟨ε⟩) well below the Maxwellian value."""
        for eps in (5.0, 6.0, 7.0):
            k2 = float(self.tab2.rate("R3", eps))
            kM = float(self.tabM.rate("R3", eps))
            self.assertLess(k2, 0.3 * kM, f"<eps>={eps}: {k2:.3e} !< 0.3*{kM:.3e}")

    def test_attachment_less_depleted_than_ionization(self) -> None:
        # F2 attachment peaks at low ε — the depleted tail hits it far less than ionization.
        eps = 5.0
        att_ratio = float(self.tab2.rate("R10", eps)) / float(self.tabM.rate("R10", eps))
        ion_ratio = float(self.tab2.rate("R3", eps)) / float(self.tabM.rate("R3", eps))
        self.assertGreater(att_ratio, 5.0 * ion_ratio)

    def test_growth_map_vs_maxwellian_compensation(self) -> None:
        """ν_ion(E/N) closure sensitivity — re-pinned after the R6 Ar(4p) channel landed.

        Empirical finding pinned as a regression: from ~30 Td up, the two-term ⟨ε⟩(E/N) sits
        well ABOVE the Maxwellian power-balance inversion (Druyvesteyn-like bulk), which
        largely compensates the tail depletion when both are read at their own self-consistent
        (E/N → ⟨ε⟩) point — the breakdown *field* barely moves between closures; guard both
        directions within a factor of 5 there. BELOW ~30 Td the tail-depletion direction
        holds (two-term ν_ion < Maxwellian) but its depth depends on the R5/R6 excitation
        set: the raw NGFSRDW sums (×4–13 SIGLO at 13–20 eV) gave 0.017–0.049; the
        SIGLO-total × RDW-branching records (2026-07-31) give ~0.4 — guard direction only.
        """
        n_gas = 8.317e25

        def nu_at(tab, en_want):
            en_grid = np.asarray(tab.reduced_field_td, dtype=np.float64)
            j = int(np.argmin(np.abs(en_grid - en_want)))
            eps = tab.key_grid[j]
            return n_gas * (
                float(tab.rate("R3", eps)) * self.mix["Ne"]
                + float(tab.rate("R7", eps)) * self.mix["Ar"]
            )

        for en in (30.0, 50.0, 100.0):
            r = nu_at(self.tab2, en) / nu_at(self.tabM, en)
            self.assertGreater(r, 0.2, f"E/N={en}: ratio {r:.3g}")
            self.assertLess(r, 5.0, f"E/N={en}: ratio {r:.3g}")
        for en in (10.0, 20.0):
            r = nu_at(self.tab2, en) / nu_at(self.tabM, en)
            self.assertGreater(r, 1.0e-3, f"E/N={en}: ratio {r:.3g}")
            self.assertLess(r, 1.0, f"E/N={en}: ratio {r:.3g}")

        # …while at FIXED ⟨ε⟩ the two-term ⟨ε⟩(E/N) map is genuinely hotter:
        j = int(np.argmin(np.abs(np.asarray(self.tab2.reduced_field_td) - 10.0)))
        eps_2t = float(self.tab2.key_grid[j])
        eps_mx = float(self.tabM.equilibrium_mean_energy(10.0))
        self.assertGreater(eps_2t, 1.2 * eps_mx)


if __name__ == "__main__":
    unittest.main()
