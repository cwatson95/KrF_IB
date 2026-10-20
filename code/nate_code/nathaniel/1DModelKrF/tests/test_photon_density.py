"""Cavity photon-density solver (optics.photon_density.PhotonDensity, Chen Eq. 7)."""

from __future__ import annotations

import math
import unittest

from optics.photon_density import PhotonDensity
from particles.constants import C_LIGHT


def _pd(**kw):
    base = dict(l_g=0.90, l_c=1.05, tau_s=4.2e-9, sigma_ste=4.0e-20,
                absorber_sigmas=(1.0e-22, 1.0e-21), r_oc=0.08, l_op=0.08, k_axial=1.0e-5)
    base.update(kw)
    return PhotonDensity(**base)


class TestPhotonDensity(unittest.TestCase):
    def test_gain_eq14(self) -> None:
        pd = _pd()
        g = pd.gain(5.0e20, (1.0e19, 1.0e18))
        expected = 5.0e20 * 4.0e-20 - (1.0e19 * 1.0e-22 + 1.0e18 * 1.0e-21)
        self.assertAlmostEqual(g, expected, places=6)

    def test_matches_constant_coefficient_analytic(self) -> None:
        # Fixed densities -> constant (a, b); n_ph(T) = (b/a)(e^{aT}-1) from n_ph(0)=0.
        pd = _pd()
        n_up = 2.0e20
        a = (pd.l_g / pd.l_c) * C_LIGHT * pd.gain(n_up, ()) - pd._loss_coeff
        b = pd.k_axial * n_up / pd.tau_s
        dt, n = 1.0e-12, 300
        for _ in range(n):
            pd.step(dt, n_up=n_up)
        analytic = (b / a) * (math.exp(a * n * dt) - 1.0)
        self.assertAlmostEqual(pd.n_ph / analytic, 1.0, places=6)

    def test_net_gain_grows_net_loss_saturates(self) -> None:
        # Strong upper-state density -> net gain (a>0) -> n_ph grows from 0.
        pd = _pd()
        for _ in range(50):
            pd.step(1.0e-12, n_up=5.0e20, absorber_densities=(1.0e19, 1.0e18))
        self.assertGreater(pd.n_ph, 0.0)
        a_gain = (pd.l_g / pd.l_c) * C_LIGHT * pd.last_gain - pd._loss_coeff
        self.assertGreater(a_gain, 0.0)
        # Weak upper-state -> net loss (a<0) -> n_ph saturates at the spontaneous floor b/|a|.
        pd2 = _pd(absorber_sigmas=())
        n_up = 1.0e17
        for _ in range(200000):
            pd2.step(1.0e-12, n_up=n_up)
        a2 = (pd2.l_g / pd2.l_c) * C_LIGHT * pd2.gain(n_up, ()) - pd2._loss_coeff
        self.assertLess(a2, 0.0)
        ss = (pd2.k_axial * n_up / pd2.tau_s) / (-a2)
        self.assertAlmostEqual(pd2.n_ph / ss, 1.0, places=4)

    def test_agrees_with_fine_step_ode(self) -> None:
        # The exponential-integrator step == a fine forward-Euler integration of Eq. 7.
        pd = _pd()
        n_up = 3.0e20
        pd.step(1.0e-11, n_up=n_up)  # one coarse exact step
        # reference: many tiny explicit Euler steps of dn/dt = a n + b
        a = (pd.l_g / pd.l_c) * C_LIGHT * pd.gain(n_up, (0.0, 0.0)) - pd._loss_coeff
        b = pd.k_axial * n_up / pd.tau_s
        n_ph = 0.0
        h = 1.0e-11 / 200000
        for _ in range(200000):
            n_ph += h * (a * n_ph + b)
        self.assertAlmostEqual(pd.n_ph / n_ph, 1.0, places=4)

    def test_lossy_cavity_required(self) -> None:
        with self.assertRaises(ValueError):
            _pd(r_oc=1.0, l_op=0.0)  # R_OC(1-L_op)=1 -> not lossy


if __name__ == "__main__":
    unittest.main()
