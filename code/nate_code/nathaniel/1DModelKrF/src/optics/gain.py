"""0D excimer gain/output optical model — Chen et al. (2025) §II.C, Eqs. 10/14/18.

The same lumped-gain cavity formulation serves KrF* (248 nm) and ArF* (193 nm). The beam runs
along the transverse optical axis, so the model has one cavity photon density ``n_ph`` (evolved by
:class:`optics.photon_density.PhotonDensity`, Chen Eq. 7) plus the algebraic output-coupling
relations.

- **Gain (Eq. 14):** ``g = n_up·σ_STE − Σ_i n_i·σ_i`` [1/m] over domain-mean densities.
- **Output-coupling efficiency (Eq. 18):**
  ``η_out = g·l_g·e^{g l_g}(1−R_OC) / [(1 + e^{g l_g}·R_OC)(e^{g l_g} − 1)]``  (→ (1−R_OC)/(1+R_OC) as g→0).
- **Output power (Eq. 10):** ``P_out = n_ph·c·hν·S·η_out``.

Species bookkeeping (the radiative reactions are excluded from the chemistry, so the optical model
owns them): each step the excited upper states decay spontaneously (``1−e^{−dt/τ}`` per radiative
channel) and the lasing upper additionally by stimulated emission (``1−e^{−c σ n_ph dt}``), with the
lower products created. UV PHOTO rows use the same frozen-fluence substeps to consume their
absorbers and create the declared ion/electron products with an exact exponential hazard.
"""

from __future__ import annotations

import math

import numpy as np

from particles.constants import C_LIGHT, EV_TO_J

from .photon_density import PhotonDensity


class GainCavity:
    """Lumped-gain excimer output over a :class:`PhotonDensity` cavity solver.

    ``radiatives`` rows are ``(upper, lifetime_s, lower_terms)`` and ``lasing`` is
    ``(upper, sigma_ste_m2, lower_terms, photon_mode)``. Each lower term is an integer
    ``(multiplicity, species_name)`` pair. ``photo_channels`` rows follow
    :func:`reactions.network.photo_channels`. Density mutations are applied in place;
    ``last_*`` ledgers are reset and populated on every :meth:`advance` call.
    """

    def __init__(
        self,
        *,
        photon_energy_ev: float,
        radiatives: list[tuple[str, float, tuple[tuple[int, str], ...]]],
        lasing: tuple[str, float, tuple[tuple[int, str], ...], str] | None,
        absorbers: list[tuple[str, float]],
        gain_length_m: float,
        cavity_length_m: float,
        beam_area_m2: float,
        r_oc: float,
        l_op: float,
        axial_se_fraction: float,
        r_back: float = 1.0,
        photo_channels: list[tuple] | tuple = (),
    ) -> None:
        self.photon_energy_ev = photon_energy_ev
        self.radiatives = list(radiatives)
        self.lasing = lasing
        self.absorbers = list(absorbers)
        self.l_g = gain_length_m
        self.l_c = cavity_length_m
        self.beam_area_m2 = beam_area_m2
        self.r_oc = r_oc
        self.l_op = l_op
        # n_ph is a cavity-averaged density in Chen Eq. 7, hence its inventory occupies l_c*S.
        self.v_mode = cavity_length_m * beam_area_m2
        # Diagnostic ledger of the last advance(): domain-mean upper-state density lost to
        # spontaneous (per upper state) and stimulated (lasing upper) emission, in m^-3.
        # Joint losses inside a substep split by competing-risks hazard shares, so
        # SE + STE equals the booked upper-state depletion exactly.
        self.last_dt_s: float = 0.0
        self.last_se_sink_m3: dict[str, float] = {}
        self.last_ste_sink_m3: dict[str, float] = {}
        self.photo_channels = list(photo_channels)
        self.last_photo_extent_m3: dict[str, float] = {}
        self.last_photoelectron_birth_m3: float = 0.0
        # Domain-mean electron energy density added by photoelectron birth (J/m^3):
        # each channel's declared electron_energy_eV (hv minus binding/detachment).
        self.last_photoelectron_birth_energy_j: float = 0.0
        self._pd: PhotonDensity | None = None
        if lasing is not None:
            _upper, sigma_ste, _lower, _mode = lasing
            self._pd = PhotonDensity(
                l_g=gain_length_m, l_c=cavity_length_m, tau_s=self._lasing_tau(),
                sigma_ste=sigma_ste, absorber_sigmas=tuple(s for _n, s in absorbers),
                r_oc=r_oc, r_back=r_back, l_op=l_op, k_axial=axial_se_fraction,
            )

    def _lasing_tau(self) -> float:
        """Radiative lifetime of the lasing upper state (from the radiative channels)."""
        if self.lasing is None:
            return math.inf
        upper = self.lasing[0]
        for u, tau, _lower in self.radiatives:
            if u == upper:
                return tau
        return math.inf

    @property
    def n_ph(self) -> float:
        return self._pd.n_ph if self._pd is not None else 0.0

    def total_photons(self) -> float:
        """Cavity photon count = n_ph · cavity mode volume (l_c·S)."""
        return self.n_ph * self.v_mode

    # -- density helpers (species live in the neutral or charged fluid) -------------------
    @staticmethod
    def _column(name: str, fluid, charged):
        if name in fluid.column:
            return fluid.density, fluid.column[name]
        if charged is not None and name in charged.column:
            return charged.density, charged.column[name]
        return None, -1

    def _mean(self, name: str, fluid, charged) -> float:
        arr, col = self._column(name, fluid, charged)
        return float(arr[:, col].mean()) if arr is not None else 0.0

    def _decay(
        self,
        upper: str,
        frac,
        products: tuple[tuple[int, str], ...],
        fluid,
        charged,
    ) -> None:
        """Apply one radiative extent per lost upper with declared product stoichiometry."""
        arr, col = self._column(upper, fluid, charged)
        if arr is None:
            return
        lost = arr[:, col] * frac
        arr[:, col] -= lost
        for count, prod in products:
            parr, pcol = self._column(prod, fluid, charged)
            if parr is not None:
                parr[:, pcol] += count * lost

    def _advance_photo_channels(self, fluence: float, fluid, charged) -> None:
        """Apply PHOTO stoichiometry with the exact exponential photon hazard."""
        if fluence <= 0.0:
            return
        photon_density_correction = 0.0
        for (reaction_id, absorber, sigma, reactants, products,
             birth_ev, _mode) in self.photo_channels:
            arr, col = self._column(absorber, fluid, charged)
            if arr is None:
                continue
            absorber_order = sum(count for count, name in reactants if name == absorber)
            if absorber_order != 1:
                raise ValueError(
                    f"PHOTO channel {reaction_id} requires first-order absorber stoichiometry"
                )
            absorber_before = arr[:, col].copy()
            hazard = C_LIGHT * float(sigma) * fluence
            lost = absorber_before * (-math.expm1(-hazard))
            arr[:, col] -= lost

            electron_before = None
            electron_energy_density_ev = None
            if charged is not None and "e-" in charged.column:
                electron_before = charged.density[:, charged.column["e-"]].copy()
                electron_energy_density_ev = electron_before * charged.electron_mean_energy_ev
            for count, product in products:
                parr, pcol = self._column(product, fluid, charged)
                if parr is not None:
                    parr[:, pcol] += count * lost
                    if product == "e-":
                        born = count * lost
                        self.last_photoelectron_birth_m3 += float(np.mean(born))
                        if electron_energy_density_ev is not None and birth_ev > 0.0:
                            # Photoelectrons carry hν minus the binding/detachment energy
                            # (declared per data row) into the pool.
                            electron_energy_density_ev = (
                                electron_energy_density_ev + born * birth_ev
                            )
                            self.last_photoelectron_birth_energy_j += float(
                                np.mean(born) * birth_ev * EV_TO_J
                            )
            if electron_before is not None:
                electron_after = charged.density[:, charged.column["e-"]]
                live = electron_after > 0.0
                charged.electron_mean_energy_ev[live] = (
                    electron_energy_density_ev[live] / electron_after[live]
                )

            extent_mean = float(np.mean(lost))
            self.last_photo_extent_m3[reaction_id] = (
                self.last_photo_extent_m3.get(reaction_id, 0.0) + extent_mean
            )
            # PhotonDensity used a frozen absorber during this substep.  Replace its linear
            # absorption extent with the exponential species extent so one reacted absorber
            # removes exactly one cavity photon after gain/cavity-volume conversion.
            frozen_extent = C_LIGHT * float(sigma) * float(np.mean(absorber_before)) * fluence
            photon_density_correction += (self.l_g / self.l_c) * (frozen_extent - extent_mean)
        if self._pd is not None and photon_density_correction != 0.0:
            self._pd.n_ph = max(self._pd.n_ph + photon_density_correction, 0.0)

    # -- Chen output relations ------------------------------------------------------------
    def output_fraction(self, g: float) -> float:
        """η_out (Eq. 18); the removable g→0 singularity limits to (1−R_OC)/(1+R_OC)."""
        x = g * self.l_g
        if abs(x) < 1.0e-9:
            return (1.0 - self.r_oc) / (1.0 + self.r_oc)
        if x > 350.0:
            # Large-gain asymptote of Eq. 18 (den → R_OC·e^{2x} would overflow):
            # η → x·(1−R_OC)·e^{−x}/R_OC → 0 (the intra-cavity intensity dwarfs the output).
            return x * (1.0 - self.r_oc) * math.exp(-x) / self.r_oc
        ex = math.exp(x)
        num = x * ex * (1.0 - self.r_oc)
        den = (1.0 + ex * self.r_oc) * (ex - 1.0)
        return num / den if den != 0.0 else 0.0

    def output_power(self, g: float) -> float:
        """Laser output power P_out = n_ph·c·hν·S·η_out (Eq. 10) [W]."""
        e_photon_j = self.photon_energy_ev * EV_TO_J
        return self.n_ph * C_LIGHT * e_photon_j * self.beam_area_m2 * self.output_fraction(g)

    # -- per-step advance -----------------------------------------------------------------
    def advance(self, dt: float, *, fluid, charged=None) -> tuple[float, float]:
        """Advance n_ph + deplete the excited states over ``dt``; return ``(P_out_W, g)``.

        The (n_ph, lasing-upper) pair is **co-integrated**: internal substeps keep the photon
        e-folding |a|·h ≲ a few, re-reading the upper-state density each substep, so above
        threshold the photon field saturably burns the excimer inventory *within* the call
        instead of ratcheting against a frozen gain (whose e^{g·l_g} would overflow). The
        non-lasing radiative uppers decay once over the whole ``dt`` (uncoupled to n_ph).
        Returns the interval-averaged output power and the **start-of-interval** gain.
        The per-call SE/STE upper-state sinks are ledgered in ``last_se_sink_m3`` /
        ``last_ste_sink_m3`` (domain-mean m^-3 over ``last_dt_s``) for diagnostics.
        """
        if dt < 0.0:
            raise ValueError("optical advance dt must be non-negative")
        self.last_dt_s = dt
        self.last_se_sink_m3 = {}
        self.last_ste_sink_m3 = {}
        self.last_photo_extent_m3 = {}
        self.last_photoelectron_birth_m3 = 0.0
        self.last_photoelectron_birth_energy_j = 0.0
        # Start-of-interval gain from the untouched densities (the returned, pinned value).
        g_first = 0.0
        if self._pd is not None:
            g_first = self._pd.gain(
                self._mean(self.lasing[0], fluid, charged),
                tuple(self._mean(sp, fluid, charged) for sp, _ in self.absorbers),
            )
        lasing_upper = self.lasing[0] if self.lasing is not None else None
        for upper, tau, lower in self.radiatives:
            if upper == lasing_upper:
                continue  # handled inside the co-integration loop
            if tau > 0.0 and math.isfinite(tau):
                frac = -math.expm1(-dt / tau)
                self.last_se_sink_m3[upper] = (
                    self.last_se_sink_m3.get(upper, 0.0)
                    + frac * self._mean(upper, fluid, charged)
                )
                self._decay(upper, frac, lower, fluid, charged)
        if self._pd is None or dt == 0.0:
            return 0.0, 0.0

        upper, sigma, lower, _mode = self.lasing
        tau_las = self._lasing_tau()
        energy_j = 0.0
        remaining = dt
        for _guard in range(65536):
            n_up = self._mean(upper, fluid, charged)
            abs_means = tuple(self._mean(sp, fluid, charged) for sp, _ in self.absorbers)
            g_now = self._pd.gain(n_up, abs_means)
            a_now = self._pd.growth_rate(g_now)
            h = remaining if a_now == 0.0 else min(remaining, 0.1 / abs(a_now))
            # Limit upper-state depletion per frozen-coefficient substep. Trial steps are
            # cheap exact scalar updates and avoid an arbitrary global timestep floor.
            while True:
                n_before = self._pd.n_ph
                _n_after, fluence = self._pd.step_with_integral(
                    h, n_up=n_up, absorber_densities=abs_means,
                )
                hazard_ste = C_LIGHT * sigma * fluence
                hazard_se = (h / tau_las
                             if tau_las > 0.0 and math.isfinite(tau_las) else 0.0)
                hazard = hazard_ste + hazard_se
                photo_hazard = max(
                    (C_LIGHT * channel[2] * fluence for channel in self.photo_channels),
                    default=0.0,
                )
                if max(hazard, photo_hazard) <= 0.05 or h <= 1.0e-15 * dt:
                    break
                self._pd.n_ph = n_before
                h *= 0.5
            frac = -math.expm1(-hazard)
            self._advance_photo_channels(fluence, fluid, charged)
            if hazard > 0.0:
                self.last_se_sink_m3[upper] = (
                    self.last_se_sink_m3.get(upper, 0.0)
                    + frac * n_up * (hazard_se / hazard)
                )
                self.last_ste_sink_m3[upper] = (
                    self.last_ste_sink_m3.get(upper, 0.0)
                    + frac * n_up * (hazard_ste / hazard)
                )
            self._decay(upper, frac, lower, fluid, charged)
            energy_j += (
                C_LIGHT * self.photon_energy_ev * EV_TO_J * self.beam_area_m2
                * self.output_fraction(g_now) * fluence
            )
            remaining -= h
            if remaining <= 1.0e-12 * dt:
                break
        if remaining > 1.0e-12 * dt:
            raise RuntimeError("optical co-integration exceeded its substep limit")
        return energy_j / dt, g_first
