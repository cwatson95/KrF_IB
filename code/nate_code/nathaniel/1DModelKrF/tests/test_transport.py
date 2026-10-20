"""Scharfetter–Gummel drift-diffusion flux kernel (Fluid1D.transport)."""

from __future__ import annotations

import unittest

import numpy as np

from Fluid1D.transport import (
    bernoulli,
    face_average,
    flux_divergence,
    sg_interior_flux,
)


class TestBernoulli(unittest.TestCase):
    def test_limit_at_zero(self) -> None:
        self.assertAlmostEqual(bernoulli(np.array([0.0]))[0], 1.0, places=12)
        self.assertAlmostEqual(bernoulli(np.array([1e-12]))[0], 1.0, places=9)

    def test_asymptotics(self) -> None:
        # B(x) -> 0 for large +x; B(x) -> -x (>0) for large -x.
        self.assertLess(bernoulli(np.array([50.0]))[0], 1e-10)
        self.assertAlmostEqual(bernoulli(np.array([-50.0]))[0], 50.0, places=6)

    def test_identity(self) -> None:
        # B(-x) - B(x) = x  (a defining Bernoulli identity).
        x = np.array([-3.0, -0.5, 0.7, 4.0])
        np.testing.assert_allclose(bernoulli(-x) - bernoulli(x), x, atol=1e-9)


class TestSGFlux(unittest.TestCase):
    def setUp(self) -> None:
        self.dx = 0.1
        self.n = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])  # dn/dx = +10

    def test_diffusion_only_limit(self) -> None:
        # v=0 -> Γ = -D dn/dx (constant for a linear profile).
        d = 2.0
        j = sg_interior_flux(self.n, np.zeros(5), np.full(5, d), self.dx)
        np.testing.assert_allclose(j, -d * 10.0, atol=1e-9)

    def test_drift_only_upwind(self) -> None:
        # D->0, v>0 -> upwind from the left cell; v<0 -> from the right.
        d = np.full(5, 1e-25)
        jp = sg_interior_flux(self.n, np.full(5, 3.0), d, self.dx)
        np.testing.assert_allclose(jp, 3.0 * self.n[:-1], rtol=1e-6)
        jm = sg_interior_flux(self.n, np.full(5, -3.0), d, self.dx)
        np.testing.assert_allclose(jm, -3.0 * self.n[1:], rtol=1e-6)

    def test_conservation_no_flux_walls(self) -> None:
        # With zero boundary fluxes, transport conserves total density.
        j = sg_interior_flux(self.n, np.full(5, 1.5), np.full(5, 0.7), self.dx)
        dndt = flux_divergence(j, 0.0, 0.0, self.dx)
        self.assertLess(abs(float(np.sum(dndt)) * self.dx), 1e-9)

    def test_face_average(self) -> None:
        np.testing.assert_allclose(face_average(np.array([1.0, 3.0, 5.0])), [2.0, 4.0])


if __name__ == "__main__":
    unittest.main()
