"""Tests for the fluid chemistry: analytic Jacobian, numba step vs SciPy Radau, conservation."""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from kinetics import FluidChemistry, FluidIntegrator
from reactions import load_reaction_network

KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"


def _seed(chem: FluidChemistry) -> np.ndarray:
    n = np.zeros(chem.n_species)
    n[chem.index["Kr"]] = 1.0e25
    n[chem.index["F2"]] = 2.0e23
    n[chem.index["Kr*"]] = 5.0e21
    return n


class TestFluidChemistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        net = load_reaction_network(KRF_TOML)
        cls.chem = FluidChemistry.from_compiled(net.compile(), net.registry)

    def test_analytic_jacobian_matches_finite_difference(self) -> None:
        # Validate at evolved states (5 ns, 50 ns): the raw seed has one row dominated
        # by a ~1e28 three-body term that swamps ~1e7 cross-derivatives, so the
        # finite-difference reference (not the analytic Jacobian) breaks there.
        chem = self.chem
        n = _seed(chem)
        for t_evolve in (5e-9, 50e-9):
            state = chem.integrate_reference(n, (0.0, t_evolve)).y[:, -1]
            Ja = chem.jacobian(state)
            f0 = chem.rhs(state)
            Jn = np.zeros_like(Ja)
            for j in range(chem.n_species):
                # Number-density RHS values are O(1e28--1e36). A one-particle/m3 absolute
                # perturbation for an unpopulated intermediate is swallowed by subtraction
                # cancellation once the corrected two-body Penning channels are present.
                # 1e12 m^-3 remains negligible relative to the 1e25 m^-3 gas while yielding
                # a resolvable finite-difference column.
                h = max(1e-6 * abs(state[j]), 1.0e12)
                pert = state.copy()
                pert[j] += h
                Jn[:, j] = (chem.rhs(pert) - f0) / h
            mask = np.abs(Ja) > 1e-3 * np.abs(Ja).max()
            rel = np.abs(Ja - Jn)[mask] / np.abs(Ja)[mask]
            self.assertLess(rel.max(), 1e-3)

    def test_atom_conservation_reference(self) -> None:
        chem = self.chem
        n0 = _seed(chem)
        a0 = chem.atom_counts(n0)
        n_end = chem.integrate_reference(n0, (0.0, 200e-9)).y[:, -1]
        a_end = chem.atom_counts(n_end)
        seeded = a0 > 0  # He/Ar/Ne/Xe are absent in this fill -> skip 0/0
        rel = np.abs(a_end - a0)[seeded] / a0[seeded]
        self.assertLess(rel.max(), 1e-8)


class TestNumbaVsReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        net = load_reaction_network(KRF_TOML)
        cls.chem = FluidChemistry.from_compiled(net.compile(), net.registry)
        cls.integ = FluidIntegrator(cls.chem)
        cls.integ.warmup()

    def test_numba_step_matches_radau(self) -> None:
        chem, integ = self.chem, self.integ
        n0 = _seed(chem)
        t_end = 100e-9
        ref = chem.integrate_reference(n0, (0.0, t_end)).y[:, -1]
        cells = n0[np.newaxis, :].copy()
        integ.advance(cells, t_end, substeps=4000)  # 1st-order Euler -> many substeps
        got = cells[0]
        # compare the physically meaningful populations to the stiff reference
        for name in ("KrF*", "Kr2F*", "Kr*", "F"):
            i = chem.index[name]
            scale = max(abs(ref[i]), 1e-3 * abs(ref).max())
            self.assertLess(abs(got[i] - ref[i]) / scale, 2e-2, msg=f"{name} mismatch")

    def test_conservation_and_parallel_cells(self) -> None:
        chem, integ = self.chem, self.integ
        n0 = _seed(chem)
        a0 = chem.atom_counts(n0)
        cells = np.repeat(n0[np.newaxis, :], 32, axis=0).copy()  # 32 identical cells
        integ.advance(cells, 100e-9, substeps=1000)
        # every cell identical (deterministic, independent) and atom-conserving
        self.assertTrue(np.allclose(cells, cells[0][np.newaxis, :]))
        seeded = a0 > 0
        rel = np.abs(chem.atom_counts(cells[0]) - a0)[seeded] / a0[seeded]
        self.assertLess(rel.max(), 1e-6)


if __name__ == "__main__":
    unittest.main()
