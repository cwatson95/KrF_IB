"""Cavity photon-density solver — Chen et al. (2025) Eq. 7, closed-quadrature form.

The excimer laser module of Chen et al., *Phys. Plasmas* **32**, 113502 (2025) §II.C evolves a
single scalar cavity photon density ``n_ph`` (m⁻³) along the optical (gain) axis — which is
**transverse to the 1D discharge axis** — driven by the plasma densities. Eq. 7 is

    dn_ph/dt = (l_g/l_c)·(R_STE − R_z) + K·R_SE − L,

    R_STE = n_ph·n_up·σ_STE       (stimulated emission),
    R_z   = n_ph·Σ_i n_i·σ_i      (UV photo-absorption, e.g. Ar** and F⁻),
    R_SE  = n_up/τ                (spontaneous emission of the upper state),
    K     ≈ 1e-5                  (fraction of isotropic SE launched along the axis),
    L     = (c·n_ph/2l_c)·ln[1/(R_OC·(1−L_op))]   (cavity loss, Eq. 8, as a positive loss rate).

This is **linear in n_ph** with plasma-driven coefficients, ``dn_ph/dt = a(t)·n_ph + b(t)``
(R_STE/R_z carry the photon speed c, so the growth term is a temporal rate in 1/s):

    a(t) = (l_g/l_c)·c·G(t) − (c/2l_c)·ln[1/(R_OC·(1−L_op))],   G = n_up·σ_STE − Σ_i n_i·σ_i,
    b(t) = K·n_up(t)/τ.

The integrating-factor (closed-quadrature) solution with χ = (c/2)·ln[1/(R_OC(1−L_op))] − l_g·G
(so a = −χ/l_c) is
    n_ph(t) = e^{−(1/l_c)∫χ ds}·(K/τ)·∫ e^{+(1/l_c)∫χ ds′}·n_up dτ ,   n_ph(0)=0.
We evaluate it **per fluid step** as the exact exponential-integrator update (a, b frozen over a
step where the plasma densities are ~constant):

    n_ph⁺ = n_ph·e^{a·Δt} + b·(e^{a·Δt} − 1)/a          (→ n_ph + b·Δt  as a → 0).

Stable against the ODE stiffness and the e^{gl_g} overflow of the raw quadrature.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from particles.constants import C_LIGHT


@dataclass
class PhotonDensity:
    """Scalar cavity photon-density solver (Chen Eq. 7, exponential-integrator step).

    ``sigma_ste`` is the stimulated-emission cross-section σ (m², Chen R38); ``absorber_sigmas``
    is the ordered tuple of UV photo-absorption cross-sections σ_i (m², e.g. R40/R41). Lengths
    ``l_g``/``l_c`` in metres, ``tau_s`` the upper-state radiative lifetime (s), ``r_oc`` the
    output-coupler reflectivity R_OC, ``l_op`` the single-pass optical loss L_op, ``k_axial`` the
    axial-SE fraction K.
    """

    l_g: float
    l_c: float
    tau_s: float
    sigma_ste: float
    absorber_sigmas: tuple[float, ...] = ()
    r_oc: float = 0.08
    r_back: float = 1.0
    l_op: float = 0.08
    k_axial: float = 1.0e-5
    c: float = C_LIGHT
    n_ph: float = 0.0
    _last_gain: float = field(default=0.0, repr=False)

    def __post_init__(self) -> None:
        # Cavity-loss rate coefficient (per n_ph): (c/2l_c)·ln[1/(R_OC(1−L_op))] ≥ 0.
        arg = self.r_oc * self.r_back * (1.0 - self.l_op)
        if not (0.0 < arg < 1.0):
            raise ValueError(
                f"R_OC·R_back·(1−L_op) must be in (0,1) for a lossy cavity; got {arg:g} "
                f"(R_OC={self.r_oc}, R_back={self.r_back}, L_op={self.l_op})"
            )
        self._loss_coeff = (self.c / (2.0 * self.l_c)) * math.log(1.0 / arg)

    def gain(self, n_up: float, absorber_densities: tuple[float, ...]) -> float:
        """Net gain coefficient g = n_up·σ_STE − Σ_i n_i·σ_i (Chen Eq. 14) [1/m]."""
        g = n_up * self.sigma_ste
        for n_i, sig_i in zip(absorber_densities, self.absorber_sigmas):
            g -= n_i * sig_i
        return g

    def growth_rate(self, g: float) -> float:
        """Temporal photon growth rate a = (l_g/l_c)·c·g − loss (1/s) for a gain ``g``."""
        return (self.l_g / self.l_c) * self.c * g - self._loss_coeff

    def step(self, dt: float, *, n_up: float, absorber_densities: tuple[float, ...] = ()) -> float:
        """Advance n_ph by one step of ``dt`` at fixed plasma densities; return the new n_ph.

        ``n_up`` = upper-state (ArF*/KrF*) density; ``absorber_densities`` aligned with
        ``absorber_sigmas``. Uses the exact exponential-integrator update of Eq. 7.
        """
        value, _integral = self.step_with_integral(
            dt, n_up=n_up, absorber_densities=absorber_densities,
        )
        return value

    def step_with_integral(
        self, dt: float, *, n_up: float, absorber_densities: tuple[float, ...] = (),
    ) -> tuple[float, float]:
        """Advance and return ``(n_ph_new, integral(n_ph dt))`` for frozen coefficients.

        The exact integral drives stimulated depletion and interval-averaged output. Using
        only the end-of-step density over-depletes the upper state during exponential growth.
        """
        if dt < 0.0:
            raise ValueError("photon step dt must be non-negative")
        g = self.gain(n_up, absorber_densities)
        self._last_gain = g
        # Temporal growth rate: photons traverse the gain region at c, so the STE/absorption
        # rate coefficient is (l_g/l_c)·c·g [1/s] (Chen R_STE = c·σ·n_ph·n_up); loss already has c.
        a = self.growth_rate(g)
        b = self.k_axial * n_up / self.tau_s if self.tau_s > 0.0 else 0.0
        ax = a * dt
        if abs(ax) < 1.0e-6:
            growth = math.exp(ax)
            src = dt * (1.0 + ax / 2.0 + ax**2 / 6.0 + ax**3 / 24.0)
            psi = dt * dt * (0.5 + ax / 6.0 + ax**2 / 24.0 + ax**3 / 120.0)
        else:
            # Exponent clamp: a caller that keeps a·dt ≲ few (GainCavity co-integrates) never
            # hits this; it only stops a pathological frozen-gain step from overflowing exp.
            growth = math.exp(min(ax, 700.0))
            src = (growth - 1.0) / a
            psi = (src - dt) / a
        n0 = self.n_ph
        integral = n0 * src + b * psi
        self.n_ph = n0 * growth + b * src
        if self.n_ph < 0.0:
            self.n_ph = 0.0
        return self.n_ph, max(float(integral), 0.0)

    @property
    def last_gain(self) -> float:
        """Most recent net gain coefficient g (1/m) from :meth:`step`."""
        return self._last_gain
