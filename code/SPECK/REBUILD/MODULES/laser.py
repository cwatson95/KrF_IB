"""Laser cavity — Chen 2025 Eqs. (7)–(18), multi-zone, energy-ledgered.

ADAPTED from spectre/laser.py (the Chen photon-equation port) with two
upgrades, both stated:

* MULTI-ZONE: the gain medium is the device grid's z-column; every rate that
  Chen writes with a single density uses the gain-volume average here, which
  reduces exactly to Chen's form for one cell. The l_g/l_c factor is the
  cavity/gain volume ratio and is what makes the photon-energy bookkeeping
  close (see ledger.py).

* THE LOSS TERM IS SPLIT for the ledger: Chen Eq. (8)'s
  L = (c n_ph / 2 l_c)[ln R_OC + ln(1 - L_op)] is one number in the photon
  equation but TWO different fates for the energy — |ln R_OC| leaves as the
  laser output, |ln(1-L_op)| dies at the HR/windows. They are integrated
  into separate named channels.

The Rigrod saturated-extraction analytic (the S6 gate's oracle) is
``rigrod_intensity``: for a homogeneously broadened two-mirror oscillator,
the round-trip-steady-state circulating intensity satisfies the standard
Rigrod relation; the gate compares the engine's quasi-steady extraction
against it in the flat-gain limit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .units import C_LIGHT, photon_energy_J


@dataclass(frozen=True)
class Cavity:
    """Resonator + gain-column geometry (Chen Table II vocabulary)."""

    l_c_m: float = 1.05
    l_g_m: float = 0.90
    K_axial: float = 1e-5
    R_OC: float = 0.08
    k_c_loss: float = 80.0        # L_op = k_c * C_F2 (Chen Eq. 9)
    L_op_fixed: float = 0.0       # if > 0: fixed single-pass loss (e-beam cells)
    area_m2: float = 90e-4
    wavelength_nm: float = 248.0
    sigma_se_m2: float = 0.0
    absorbers: dict = field(default_factory=dict)

    @property
    def photon_J(self) -> float:
        return photon_energy_J(self.wavelength_nm)

    @property
    def cavity_vol_m3(self) -> float:
        return self.area_m2 * self.l_c_m

    @property
    def gain_vol_m3(self) -> float:
        return self.area_m2 * self.l_g_m

    def single_pass_loss(self, c_f2: float) -> float:
        if self.L_op_fixed > 0.0:
            return self.L_op_fixed
        return min(self.k_c_loss * max(c_f2, 0.0), 0.9999)


def gain_m1(cav: Cavity, dens_avg: dict, upper: str) -> float:
    """g [1/m] = n_upper sigma_se - sum absorbers (Chen Eq. 14), on the
    gain-volume averaged densities."""
    g = dens_avg.get(upper, 0.0) * cav.sigma_se_m2
    for sp, sig in cav.absorbers.items():
        g -= dens_avg.get(sp, 0.0) * sig
    return g


def loss_rates(cav: Cavity, n_ph: float, c_f2: float):
    """(output-coupling, internal-loss) photon drain rates [m^-3 s^-1],
    both >= 0 — the ledger split of Chen Eq. (8)."""
    pre = C_LIGHT * max(n_ph, 0.0) / (2.0 * cav.l_c_m)
    out = pre * (-math.log(max(cav.R_OC, 1e-9)))
    internal = pre * (-math.log(max(1.0 - cav.single_pass_loss(c_f2), 1e-9)))
    return out, internal


def photon_rhs(cav: Cavity, r_ste_avg: float, r_x_avg: float, r_se_avg: float,
               n_ph: float, c_f2: float) -> float:
    """dn_ph/dt (Chen Eq. 7) from gain-volume-averaged rates."""
    out, internal = loss_rates(cav, n_ph, c_f2)
    return ((cav.l_g_m / cav.l_c_m) * (r_ste_avg - r_x_avg)
            + cav.K_axial * r_se_avg - out - internal)


def output_power_W(cav: Cavity, n_ph: float) -> float:
    """Power through the output coupler [W] — the ledgered quantity:
    P = (OC drain rate) x h nu x cavity volume."""
    out, _ = loss_rates(cav, n_ph, 0.0)   # OC term has no c_f2 dependence
    return out * cav.photon_J * cav.cavity_vol_m3


def eta_out_chen(g_m1: float, cav: Cavity) -> float:
    """Chen Eq. (18) pulsed output-coupling efficiency (reporting overlay)."""
    gl = g_m1 * cav.l_g_m
    if gl <= 1e-9:
        return 0.0
    egl = math.exp(min(gl, 50.0))
    return gl * egl * (1.0 - cav.R_OC) / ((1.0 + egl * cav.R_OC)
                                          * (egl - 1.0))


# ------------------------------------------------------------------ Rigrod
def rigrod_intensity(g0_m1: float, l_g_m: float, R1: float, R2: float,
                     I_sat_W_m2: float, alpha_m1: float = 0.0) -> float:
    """Rigrod (1965) steady-state output intensity [W/m^2] of a homogeneous
    two-mirror oscillator, mirror-2 output; the space-independent
    approximation with distributed loss:

        I_out = (1-R2)/( (1+sqrt(R2/R1)) (1-sqrt(R1 R2)) )
                * I_sat * (g0 l - |ln sqrt(R1 R2)| - alpha l)

    Valid above threshold; returns 0 below. This is the ANALYTIC OVERLAY the
    extraction gate scores against — it never feeds back into the engine.
    """
    loss = -0.5 * math.log(max(R1 * R2, 1e-30)) + alpha_m1 * l_g_m
    net = g0_m1 * l_g_m - loss
    if net <= 0.0:
        return 0.0
    s1, s2 = math.sqrt(max(R1, 1e-30)), math.sqrt(max(R2, 1e-30))
    pre = (1.0 - R2) / ((1.0 + s2 / s1) * (1.0 - s1 * s2))
    return pre * I_sat_W_m2 * net
