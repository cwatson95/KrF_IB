"""laser.py gates: the ledgered loss split and the Rigrod overlay (S6)."""

import math

import numpy as np

from REBUILD.MODULES.laser import (Cavity, eta_out_chen, gain_m1, loss_rates,
                                   output_power_W, photon_rhs,
                                   rigrod_intensity)
from REBUILD.MODULES.units import C_LIGHT


def _cav(**kw):
    d = dict(l_c_m=0.75, l_g_m=0.75, K_axial=1e-5, R_OC=0.1, k_c_loss=0.0,
             L_op_fixed=0.088, area_m2=2e-4, wavelength_nm=248.0,
             sigma_se_m2=5.354e-21, absorbers={})
    d.update(kw)
    return Cavity(**d)


def test_loss_split_reproduces_chen_eq8():
    """out + internal must equal the single Chen Eq. (8) magnitude
    (c n_ph / 2 l_c) |ln(R_OC (1-L_op))|."""
    cav = _cav()
    n_ph = 3e17
    out, internal = loss_rates(cav, n_ph, 0.0)
    total_ref = (C_LIGHT * n_ph / (2 * cav.l_c_m)) * (
        -math.log(cav.R_OC * (1.0 - 0.088)))
    assert abs((out + internal) / total_ref - 1.0) < 1e-12
    assert out > 0 and internal > 0


def test_photon_rhs_signs():
    cav = _cav()
    # no gain, no spont: photons decay
    assert photon_rhs(cav, 0.0, 0.0, 0.0, 1e17, 0.0) < 0.0
    # strong gain: photons grow
    assert photon_rhs(cav, 1e30, 0.0, 0.0, 1e17, 0.0) > 0.0
    # spontaneous seeding from zero photons
    assert photon_rhs(cav, 0.0, 0.0, 1e25, 0.0, 0.0) > 0.0


def test_gain_and_output_power():
    cav = _cav(absorbers={"F-": 5.6e-22})
    dens = {"KrF*": 1e20, "F-": 1e20}
    g = gain_m1(cav, dens, "KrF*")
    assert abs(g - (1e20 * 5.354e-21 - 1e20 * 5.6e-22)) < 1e-12
    p = output_power_W(cav, 1e17)
    out, _ = loss_rates(cav, 1e17, 0.0)
    assert abs(p - out * cav.photon_J * cav.cavity_vol_m3) < 1e-6 * p


def test_rigrod_threshold_and_linearity():
    Isat = 1e8
    lg, R1, R2 = 0.75, 0.99, 0.7
    g_th = -0.5 * math.log(R1 * R2) / lg
    assert rigrod_intensity(g_th * 0.999, lg, R1, R2, Isat) == 0.0
    I1 = rigrod_intensity(g_th + 0.5, lg, R1, R2, Isat)
    I2 = rigrod_intensity(g_th + 1.0, lg, R1, R2, Isat)
    assert I1 > 0.0
    # linear in (g0 l - loss): doubling the excess doubles the output
    assert abs(I2 / I1 - 2.0) < 1e-9


def test_eta_out_chen_form():
    """Chen Eq. (18) is a photon-to-output COUPLING factor, not an
    efficiency bounded by 1 (measured: 1.31 at gl = 1.5, R_OC = 0.08). The
    engine's LEDGERED output deliberately uses loss_rates instead; Eq. 18
    stays a reporting overlay. Gates: the small-gl limit is (1-R)/(1+R),
    values are positive and finite."""
    cav = _cav(R_OC=0.08)
    small = eta_out_chen(1e-6, cav)
    assert abs(small - (1 - 0.08) / (1 + 0.08)) < 1e-3
    for g in (0.5, 2.0, 5.0):
        e = eta_out_chen(g, cav)
        assert 0.0 < e < 10.0
    assert eta_out_chen(0.0, cav) == 0.0
