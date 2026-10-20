"""Analytic overlays that SCORE this code -- they never steer a number.

Policy (PLAN (26)): every comparison target here is derived from theory, not
fitted to the code's output.  Where the DISCRETE scheme has a known dispersion
error, the honest target is the discrete prediction (which converges to the
continuum answer as h -> 0), not the continuum value with a loose tolerance.

THE COLD-PLASMA MODE, DISCRETELY
--------------------------------
A cold, uniform, unmagnetised plasma oscillates at exactly omega_pe in the
continuum.  A PIC code does not reproduce that exactly, and the deviation is
NOT an error to be tuned away -- it is a computable property of the scheme.
For this package (trilinear/CIC deposit and gather, centred-difference
gradient, spectrally exact Poisson) the linearised cold-fluid analysis gives

    omega^2 / omega_pe^2  =  k * k_g(k) * S(k)^2 / K(k)^2

with, for cell size h and a mode k:

    K(k)^2  = (2 - 2 cos(kh)) / h^2      the discrete Laplacian eigenvalue
    k_g(k)  = sin(kh) / h                the centred-difference gradient
    S(k)    = [sin(kh/2) / (kh/2)]^2     the CIC shape factor (deposit AND
                                         gather, hence S^2)

Continuity uses the TRUE k because particle motion is exact between steps.
Each factor is < 1, so a PIC cold mode always runs slightly SLOW, and the
deficit scales as (kh)^2 -- i.e. the measured frequency converges to omega_pe
at second order in h, which is the gate.

Verified against this package's own runs (2026-08-03, periodic 3-D box, cold
quiet start, mode 1): predicted 0.99141 vs measured 0.99107 at nx = 24;
predicted 0.99523 vs measured 0.99391 at nx = 32.
"""

from __future__ import annotations

import numpy as np


def cic_shape_factor(kh: float) -> float:
    """CIC (trilinear) assignment-function transform S(k), normalised S(0)=1."""
    u = 0.5 * np.asarray(kh, dtype=float)
    s = np.where(np.abs(u) < 1e-12, 1.0, np.sin(u) / np.where(u == 0, 1.0, u))
    return float(s ** 2) if np.ndim(kh) == 0 else s ** 2


def discrete_laplacian_eigenvalue(kh: float, h: float) -> float:
    """K^2(k) of the 7-point stencil for a single-axis mode [1/m^2]."""
    return float((2.0 - 2.0 * np.cos(kh)) / (h * h))


def centred_gradient_wavenumber(kh: float, h: float) -> float:
    """k_g(k) = sin(kh)/h -- the centred-difference gradient's effective k."""
    return float(np.sin(kh) / h)


def cold_mode_frequency_ratio(kh: float) -> float:
    """omega / omega_pe for a cold plasma mode of this scheme (h-independent).

    Depends only on the dimensionless kh, since every factor scales as 1/h^2.
    Returns 1.0 in the limit kh -> 0.
    """
    kh = float(kh)
    if kh == 0.0:
        return 1.0
    h = 1.0                                   # cancels; kept for readability
    k = kh / h
    k_g = centred_gradient_wavenumber(kh, h)
    K2 = discrete_laplacian_eigenvalue(kh, h)
    S = cic_shape_factor(kh)
    ratio2 = k * k_g * S * S / K2
    return float(np.sqrt(max(ratio2, 0.0)))


def cold_mode_frequency(omega_pe: float, k: float, h: float) -> float:
    """The scheme's predicted cold-mode angular frequency [rad/s]."""
    return omega_pe * cold_mode_frequency_ratio(k * h)


def zero_crossing_frequency(t: np.ndarray, signal: np.ndarray) -> float:
    """Frequency [Hz] of an oscillating signal, from interpolated crossings.

    Far more precise than an FFT peak on a short record: an FFT over N periods
    resolves only 1/N of the frequency, while a straight-line fit through the
    interpolated zero-crossing times uses every crossing.  Returns 0.0 if
    fewer than three crossings are present.
    """
    t = np.asarray(t, dtype=float)
    s = np.asarray(signal, dtype=float)
    s = s - s.mean()
    idx = np.nonzero(np.sign(s[:-1]) != np.sign(s[1:]))[0]
    if idx.size < 3:
        return 0.0
    denom = s[idx + 1] - s[idx]
    tc = t[idx] + (t[idx + 1] - t[idx]) * (-s[idx] / denom)
    half_period = float(np.polyfit(np.arange(tc.size), tc, 1)[0])
    if half_period <= 0.0:
        return 0.0
    return 1.0 / (2.0 * half_period)


def child_langmuir_current(V: float, d: float, mass: float,
                           charge: float) -> float:
    """Space-charge-limited current density [A/m^2] -- the M4 sheath target.

    j = (4/9) eps0 sqrt(2 q / M) V^{3/2} / d^2.  Present now so the M4 gate
    has its target already written down and un-tunable; nothing in this
    package produces a sheath yet (missing module 4).
    """
    from .units import EPS0
    return float((4.0 / 9.0) * EPS0 * np.sqrt(2.0 * abs(charge) / mass)
                 * V ** 1.5 / (d * d))


def bohm_speed(Te_eV: float, ion_mass: float) -> float:
    """sqrt(q Te / M_i) [m/s] -- the M4 Bohm-criterion target."""
    from .units import EV
    return float(np.sqrt(max(Te_eV, 0.0) * EV / ion_mass))
