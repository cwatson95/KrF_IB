"""Maxwellian-EEDF electron coefficients from the shipped measured cross-sections.

This diagnostic coefficient backend takes moments of the SIGLO-derived datasets in
``data/cross_sections/`` (Chen et al. 2025 ref. 22). It uses a Maxwellian electron energy
distribution; the production ``two_term`` backend solves the in-tree two-term equation over the
same datasets.

With ``F₀(ε) = (2/√π)·Tₑ^{-3/2}·e^{−ε/Tₑ}`` (normalised ``∫F₀√ε dε = 1``) and
``γ = √(2e/mₑ) ≈ 5.931e5 m·s⁻¹·eV^{-1/2}`` (Hagelaar & Pitchford 2005, Eqs. 61–63):

- **rate coefficients**   ``k_r(Tₑ) = γ ∫ ε·σ_r(ε)·F₀(ε) dε``                        [m³/s]
- **reduced mobility**    ``μ_e·N = −(γ/3) ∫ (ε/σ_m(ε))·(∂F₀/∂ε) dε
                                  = (γ/3Tₑ) ∫ (ε/σ_m(ε))·F₀(ε) dε``                [1/(V·m·s)]
- **reduced diffusion**   ``D_e·N = (γ/3) ∫ (ε/σ_m(ε))·F₀(ε) dε = Tₑ·μ_e·N``        [1/(m·s)]

so the Einstein relation ``D/μ = Tₑ`` holds *exactly*, by construction. ``σ_m`` is the
mixture-weighted momentum-transfer cross-section ``Σ_g x_g σ_{m,g}(ε)``.

The table is keyed directly on the **mean energy** ``⟨ε⟩ = 3Tₑ/2`` (``key_kind="mean_energy_eV"``),
which is what the LMEA transport model looks up — no E/N→Tₑ guesswork anywhere. A reduced field is
*derived* for reporting from the steady-state power balance (see :meth:`_reduced_field`).

**Approximation.** A Maxwellian EEDF is exact only when electron–electron collisions dominate; real
discharge EEDFs are depleted in the inelastic tail, so ionization rates here are upper estimates at
low Tₑ. It is nonetheless a standard fluid-model closure, and — unlike the stub — it is built from
measured cross-sections and reproduces the attachment/ionization competition that sets breakdown.
"""

from __future__ import annotations

import hashlib
import warnings
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from particles.constants import E_CHARGE, M_E, molecular_mass
from reactions.cross_sections import (
    TabulatedSigma,
    build_cross_section_index,
    load_cross_section,
    resolve_cross_section_path,
)
from reactions.types import RateForm, ReactionKind

from .boltzmann import BoltzmannSolver, CoefficientTable, _dynamic_reactions, _threshold_ev

# γ = √(2e/mₑ): electron speed (m/s) per √eV.
GAMMA_V_PER_SQRT_EV = float(np.sqrt(2.0 * E_CHARGE / M_E))
# Floor on σ_m so the ε/σ_m integrand stays finite where a table has a deep Ramsauer minimum.
_SIGMA_M_FLOOR_M2 = 1.0e-24


def _cross_section_fingerprint(path: Path) -> str:
    """Content fingerprint of the JSON/CSV inputs that determine a coefficient solve."""
    digest = hashlib.sha256()
    for item in sorted((*path.rglob("*.json"), *path.rglob("*.csv"))):
        digest.update(str(item.relative_to(path)).encode("utf-8"))
        digest.update(item.read_bytes())
    return digest.hexdigest()[:16]


def maxwellian_f0(energy_ev: NDArray[np.float64], te_ev: float) -> NDArray[np.float64]:
    """Normalised Maxwellian EEDF F₀(ε) (eV^{-3/2}); ``∫F₀·√ε dε = 1``."""
    return (2.0 / np.sqrt(np.pi)) * te_ev**-1.5 * np.exp(-energy_ev / te_ev)


class MaxwellianEedfSolver(BoltzmannSolver):
    """Electron rate + transport coefficients by Maxwellian-EEDF moments of measured σ(ε).

    ``cross_sections_dir`` is the JSON+CSV cross-section directory. Per-reaction rates use each
    reaction's own ``cross_section`` dataset; transport uses the per-gas **elastic /
    momentum-transfer** datasets (``ReactionKind.ELASTIC``), which are not reactions and are found
    by target. ``extra_loss_process_ids`` names further inelastic channels (e.g. F₂ vibrational)
    that carry energy but change no species — they enter only the power balance behind the derived
    reduced field.

    ``solve`` ignores the *values* of ``en_grid_td`` and uses its length to size the mean-energy
    abscissa (this solver is natively mean-energy-keyed); the LUT cache key still tracks it.
    """

    def __init__(
        self,
        cross_sections_dir: str | Path,
        *,
        mean_energy_range_ev: tuple[float, float] = (0.05, 30.0),
        energy_max_ev: float = 400.0,
        quad_points: int = 4000,
        extra_loss_process_ids: tuple[str, ...] = (),
        sigma_lut_points: int = 1024,
    ) -> None:
        self.cross_sections_dir = Path(cross_sections_dir)
        self.mean_energy_range_ev = mean_energy_range_ev
        self.energy_max_ev = float(energy_max_ev)
        self.quad_points = int(quad_points)
        self.extra_loss_process_ids = tuple(extra_loss_process_ids)
        self.sigma_lut_points = int(sigma_lut_points)
        self._index: dict[str, Path] | None = None
        self._sigma_cache: dict[str, TabulatedSigma] = {}

    @property
    def solver_id(self) -> str:
        lo, hi = self.mean_energy_range_ev
        return (f"maxwell-eedf:v1:eps{lo:g}-{hi:g}:emax{self.energy_max_ev:g}"
                f":nq{self.quad_points}:nlut{self.sigma_lut_points}"
                f":loss{','.join(self.extra_loss_process_ids)}"
                f":xs{_cross_section_fingerprint(self.cross_sections_dir)}")

    # ------------------------------------------------------------------ cross-section access
    def _process_index(self) -> dict[str, Path]:
        if self._index is None:
            self._index = build_cross_section_index(self.cross_sections_dir)
        return self._index

    def _sigma(self, key: str) -> TabulatedSigma | None:
        """Load σ(ε) for a dataset id/path as a clamped, non-negative lookup (cached)."""
        if key in self._sigma_cache:
            return self._sigma_cache[key]
        try:
            record = load_cross_section(resolve_cross_section_path(key, self.cross_sections_dir))
            sigma = TabulatedSigma.from_record(record, n_points=self.sigma_lut_points)
        except (FileNotFoundError, ValueError, KeyError):
            return None
        self._sigma_cache[key] = sigma
        return sigma

    def _elastic_sigma_for(self, gas: str) -> TabulatedSigma | None:
        """The momentum-transfer (elastic) dataset for a background gas, found by target."""
        for pid, path in sorted(self._process_index().items()):
            try:
                record = load_cross_section(path)
            except (FileNotFoundError, ValueError):
                continue
            if record.target == gas and record.reaction_kind is ReactionKind.ELASTIC:
                sigma = TabulatedSigma.from_record(record, n_points=self.sigma_lut_points)
                self._sigma_cache[pid] = sigma
                return sigma
        return None

    # ------------------------------------------------------------------ solve
    def solve(self, mixture, en_grid_td, gas_temp_k=300.0, *, reactions=None) -> CoefficientTable:
        n_key = max(8, int(np.asarray(en_grid_td).size))
        lo, hi = self.mean_energy_range_ev
        mean_energy = np.geomspace(lo, hi, n_key)
        te = (2.0 / 3.0) * mean_energy

        # Quadrature grid in ε (log-spaced: the integrands are exponentially weighted).
        eps = np.geomspace(1.0e-3, self.energy_max_ev, self.quad_points)

        sigma_m = self._mixture_momentum_transfer(mixture, eps)
        rate_coeff = self._reaction_rates(reactions, eps, te)

        # Transport moments (vectorised over the Tₑ grid).
        mobility_N = np.empty_like(te)
        diffusion_N = np.empty_like(te)
        integrand_pref = eps / sigma_m
        for i, t in enumerate(te):
            f0 = maxwellian_f0(eps, t)
            moment = float(np.trapezoid(integrand_pref * f0, eps))  # ∫(ε/σ_m)F₀ dε
            diffusion_N[i] = (GAMMA_V_PER_SQRT_EV / 3.0) * moment
            mobility_N[i] = diffusion_N[i] / t  # = (γ/3Tₑ)·moment → Einstein exact

        en_td, drift = self._reduced_field(mixture, eps, te, mobility_N, rate_coeff, reactions)

        return CoefficientTable(
            key_kind="mean_energy_eV",
            key_grid=mean_energy,
            mobility_N=mobility_N,
            diffusion_N=diffusion_N,
            drift_velocity=drift,
            mean_energy_ev=mean_energy,
            electron_temp_ev=te,
            rate_coeff=rate_coeff,
            mixture=dict(mixture),
            gas_temp_k=float(gas_temp_k),
            # energy pair left None: consumers apply the (5/3) closure.
            reduced_field_td=en_td,  # enables the local-field ⟨ε⟩_eq(E/N) closure
        )

    # ------------------------------------------------------------------ pieces
    def _mixture_momentum_transfer(self, mixture, eps) -> NDArray[np.float64]:
        """σ_m(ε) = Σ_g x_g·σ_{m,g}(ε) over the background gases (m²), floored."""
        sigma_m = np.zeros_like(eps)
        missing: list[str] = []
        for gas, frac in mixture.items():
            if frac <= 0.0:
                continue
            sig = self._elastic_sigma_for(gas)
            if sig is None:
                missing.append(gas)
                continue
            sigma_m += float(frac) * np.asarray(sig(eps), dtype=np.float64)
        if missing:
            warnings.warn(
                f"no elastic/momentum-transfer cross-section for {sorted(missing)} in "
                f"{self.cross_sections_dir}; electron mobility/diffusion will be OVERESTIMATED "
                f"(those gases contribute no momentum transfer)", stacklevel=3,
            )
        if not np.any(sigma_m > 0.0):
            raise ValueError(
                f"no momentum-transfer cross-section found for any gas in mixture "
                f"{sorted(mixture)}; cannot compute electron transport"
            )
        return np.maximum(sigma_m, _SIGMA_M_FLOOR_M2)

    def _reaction_rates(self, reactions, eps, te) -> dict[str, NDArray[np.float64]]:
        """k_r(Tₑ) = γ∫ε·σ_r·F₀ dε for every electron-impact reaction carrying a dataset."""
        rate_coeff: dict[str, NDArray[np.float64]] = {}
        if reactions is None:
            return rate_coeff
        no_sigma: list[str] = []
        for rxn in _dynamic_reactions(reactions):
            if rxn.rate_form is not RateForm.ELECTRON_IMPACT:
                # TE_DEPENDENT reactions are evaluated analytically from their own (K0, p, E0)
                # by the chemistry layer, so they need no LUT entry.
                continue
            key = getattr(rxn, "cross_section", None)
            sig = self._sigma(str(key)) if key else None
            if sig is None:
                no_sigma.append(rxn.reaction_id or rxn.describe())
                continue
            sigma_eps = np.asarray(sig(eps), dtype=np.float64)
            weight = eps * sigma_eps
            rate_coeff[rxn.reaction_id] = np.array(
                [GAMMA_V_PER_SQRT_EV * float(np.trapezoid(weight * maxwellian_f0(eps, t), eps))
                 for t in te]
            )
        if no_sigma:
            warnings.warn(
                f"electron-impact reactions without a usable cross-section get k=0: "
                f"{sorted(no_sigma)}", stacklevel=3,
            )
        return rate_coeff

    def _reduced_field(self, mixture, eps, te, mobility_N, rate_coeff, reactions):
        """Derived E/N (Td) and drift velocity from the steady-state power balance.

        Equating per-electron field heating to collisional losses,
        ``e·(μ_e·N)·(E/N)² = Σ_r ε_r·k_r·x_r + elastic``, which is independent of N. Diagnostic
        only — the fluid model keys on ⟨ε⟩ — but it lets an E/N-keyed consumer use this table.
        """
        # Inelastic loss per electron per unit target density: Σ_r ε_th,r·k_r·x_r.
        loss = np.zeros_like(te)
        if reactions is not None:
            for rxn in _dynamic_reactions(reactions):
                rid = rxn.reaction_id
                if rid not in rate_coeff:
                    continue
                eth = _threshold_ev(rxn)
                if eth <= 0.0:
                    continue
                target_frac = sum(
                    float(mixture.get(name, 0.0)) for _c, name in rxn.reactants
                )
                loss += eth * E_CHARGE * rate_coeff[rid] * target_frac
        # Extra loss-only channels (e.g. F₂ vibrational): energy sinks, no species change.
        for pid in self.extra_loss_process_ids:
            sig = self._sigma(pid)
            if sig is None:
                continue
            eth = float(sig.threshold_ev)
            frac = float(mixture.get(sig.record.target, 0.0))
            if eth <= 0.0 or frac <= 0.0:
                continue
            weight = eps * np.asarray(sig(eps), dtype=np.float64)
            k = np.array([GAMMA_V_PER_SQRT_EV * float(np.trapezoid(weight * maxwellian_f0(eps, t), eps))
                          for t in te])
            loss += eth * E_CHARGE * k * frac
        # Elastic (recoil) loss: fraction 2mₑ/M per momentum-transfer collision.
        loss += self._elastic_power_loss(mixture, eps, te)

        with np.errstate(divide="ignore", invalid="ignore"):
            en_sq = np.where(mobility_N > 0.0, loss / (E_CHARGE * mobility_N), 0.0)
        en_v_m2 = np.sqrt(np.maximum(en_sq, 0.0))  # E/N in V·m²
        drift = mobility_N * en_v_m2               # μ·N·(E/N) = drift velocity (m/s)
        return en_v_m2 / 1.0e-21, drift            # E/N in townsend

    def _elastic_power_loss(self, mixture, eps, te) -> NDArray[np.float64]:
        """Elastic recoil power loss per electron per unit gas density (W·m³)."""
        out = np.zeros_like(te)
        for gas, frac in mixture.items():
            if frac <= 0.0:
                continue
            sig = self._elastic_sigma_for(gas)
            if sig is None:
                continue
            try:
                mass = molecular_mass(*_atoms_of(gas))
            except (KeyError, ValueError):
                continue
            sigma_eps = np.asarray(sig(eps), dtype=np.float64)
            # ⟨ν_m·ε⟩ moment: γ∫ ε·σ_m·F₀·ε dε ; times 2mₑ/M for the fractional recoil loss.
            weight = eps * sigma_eps * eps
            mom = np.array([GAMMA_V_PER_SQRT_EV * float(np.trapezoid(weight * maxwellian_f0(eps, t), eps))
                            for t in te])
            out += float(frac) * (2.0 * M_E / mass) * mom * E_CHARGE
        return out


def _atoms_of(gas: str) -> tuple[str, ...]:
    """Constituent atoms of a simple background-gas formula (``Ne`` → (Ne,), ``F2`` → (F,F))."""
    if len(gas) > 1 and gas[-1].isdigit():
        return tuple([gas[:-1]] * int(gas[-1]))
    return (gas,)
