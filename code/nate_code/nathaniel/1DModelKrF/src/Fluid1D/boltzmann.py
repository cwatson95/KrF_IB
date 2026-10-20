"""Electron Boltzmann coupling: transport + rate-coefficient lookup tables.

The 1D charged-species fluid does not solve the electron Boltzmann equation in the time
loop. Instead, for each gas mixture (i.e. each F2-concentration sweep point) a two-term
Boltzmann solve is run **once** over a reduced-field (E/N) grid to produce a
:class:`CoefficientTable`: electron transport coefficients (reduced mobility/diffusion,
drift velocity, mean energy, Te) and per-reaction electron-impact rate coefficients, as
functions of the reduced field (or mean energy). The fluid interpolates this table at
runtime — see :func:`build_lut` for the offline precompute + on-disk cache.

Backends implement :class:`BoltzmannSolver`:

- :class:`AnalyticStubSolver` — closed-form placeholder coefficients so the program and
  tests run with no external binary. **Not physical** — a baseline stand-in.
- :class:`Fluid1D.eedf_solver.MaxwellianEedfSolver` and
  :class:`Fluid1D.two_term.TwoTermBoltzmannSolver` — the real coefficient sources, built
  from the measured cross-sections in ``data/cross_sections/``.
"""

from __future__ import annotations

import hashlib
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from particles.constants import TD
from reactions.network import CompiledNetwork
from reactions.types import RateForm

# Rate forms whose rate constant comes from the Boltzmann LUT (not a fixed TOML value).
_DYNAMIC_RATE_FORMS = (RateForm.ELECTRON_IMPACT, RateForm.TE_DEPENDENT)


def reduced_field_td(
    e_field_mag: NDArray[np.floating] | float, gas_number_density: NDArray[np.floating] | float
) -> NDArray[np.float64]:
    """Reduced field |E|/N in townsend (Td) — the ``EN_Td`` LUT lookup key."""
    e = np.asarray(e_field_mag, dtype=np.float64)
    n = np.asarray(gas_number_density, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        en = np.where(n > 0.0, e / n, 0.0) / TD
    return en


def _dynamic_reactions(compiled: CompiledNetwork):
    """Reactions whose rate constant is Boltzmann-derived (electron-impact / Te-dependent)."""
    return [
        r for r in (*compiled.fluid, *compiled.kinetic)
        if r.rate_form in _DYNAMIC_RATE_FORMS
    ]


def _threshold_ev(reaction) -> float:
    params = getattr(reaction, "params", None) or {}
    for key in ("threshold_eV", "threshold_ev", "threshold", "E0", "Eth"):
        if key in params:
            try:
                return float(params[key])
            except (TypeError, ValueError):
                pass
    return 0.0


@dataclass(frozen=True)
class CoefficientTable:
    """Interpolable electron transport + per-reaction rate coefficients vs a lookup key.

    ``key_kind`` is ``"EN_Td"`` (reduced field, townsend) or ``"mean_energy_eV"``; ``key_grid``
    is the monotone-increasing lookup axis. Transport arrays and each ``rate_coeff`` entry
    are sampled on ``key_grid`` and linearly interpolated (clamped outside the range).
    ``rate_coeff`` is keyed by ``Reaction.reaction_id`` so a network reaction's runtime rate
    can be looked up directly.
    """

    key_kind: str
    key_grid: NDArray[np.float64]
    mobility_N: NDArray[np.float64]  # reduced mobility μ·N [1/(V·m·s)]
    diffusion_N: NDArray[np.float64]  # reduced diffusion D·N [1/(m·s)]
    drift_velocity: NDArray[np.float64]  # electron drift speed [m/s]
    mean_energy_ev: NDArray[np.float64]  # <eps>(key) [eV]
    electron_temp_ev: NDArray[np.float64]  # Te = 2/3 <eps> [eV]
    rate_coeff: dict[str, NDArray[np.float64]]  # reaction_id -> k(key) [SI]
    mixture: dict[str, float]
    gas_temp_k: float
    # Energy-transport pair for the electron energy PDE (LMEA). Reduced energy mobility μ_ε·N
    # and energy diffusion D_ε·N (Hagelaar & Pitchford 2005 §3.2, Eqs. 61-62). When a solver
    # does not emit them, the textbook/COMSOL closure μ_ε=(5/3)μ_e, D_ε=(5/3)D_e is used (see
    # the accessors) — sufficient for the fluid energy equation.
    energy_mobility_N: NDArray[np.float64] | None = None
    energy_diffusion_N: NDArray[np.float64] | None = None
    # Reduced field E/N (Td) corresponding to each key — the steady-state (local-field) pairing
    # of ⟨ε⟩ with E/N. Lets a consumer invert the relation via :meth:`equilibrium_mean_energy`,
    # which is how the LMEA closure is repaired where nₑ → 0 (see DriftDiffusionPlasmaModel).
    reduced_field_td: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        if self.key_kind not in {"EN_Td", "mean_energy_eV"}:
            raise ValueError(f"unsupported coefficient-table key kind {self.key_kind!r}")
        key = np.asarray(self.key_grid, dtype=np.float64)
        if key.ndim != 1 or key.size == 0 or not np.all(np.isfinite(key)):
            raise ValueError("coefficient-table key_grid must be a non-empty finite vector")
        if np.any(np.diff(key) <= 0.0):
            raise ValueError("coefficient-table key_grid must be strictly increasing")
        for name in (
            "mobility_N", "diffusion_N", "drift_velocity", "mean_energy_ev",
            "electron_temp_ev", "energy_mobility_N", "energy_diffusion_N",
            "reduced_field_td",
        ):
            value = getattr(self, name)
            if value is not None and np.asarray(value).shape != key.shape:
                raise ValueError(f"coefficient-table {name} must match key_grid")
        for rid, value in self.rate_coeff.items():
            if np.asarray(value).shape != key.shape:
                raise ValueError(f"rate coefficient {rid!r} must match key_grid")

    def _interp(self, arr: NDArray[np.float64], key) -> NDArray[np.float64] | float:
        return np.interp(key, self.key_grid, arr)

    def mobility(self, key):
        return self._interp(self.mobility_N, key)

    def diffusion(self, key):
        return self._interp(self.diffusion_N, key)

    def energy_mobility(self, key):
        """Reduced energy mobility μ_ε·N (1/(V·m·s)); (5/3)·μ_e·N if the solver omits it."""
        arr = self.energy_mobility_N if self.energy_mobility_N is not None else (5.0 / 3.0) * self.mobility_N
        return self._interp(arr, key)

    def energy_diffusion(self, key):
        """Reduced energy diffusion D_ε·N (1/(m·s)); (5/3)·D_e·N if the solver omits it."""
        arr = self.energy_diffusion_N if self.energy_diffusion_N is not None else (5.0 / 3.0) * self.diffusion_N
        return self._interp(arr, key)

    @property
    def mean_energy_bounds_ev(self) -> tuple[float, float]:
        """(min, max) mean energy the table actually resolves — lookups clamp outside it."""
        arr = np.asarray(self.mean_energy_ev, dtype=np.float64)
        return float(arr.min()), float(arr.max())

    def equilibrium_mean_energy(self, en_td):
        """Local-field-equilibrium mean energy ⟨ε⟩_eq (eV) at reduced field ``en_td`` (Td).

        The inverse of the table's own ⟨ε⟩ → E/N steady-state pairing. This is the closure the
        LMEA transport model falls back on where the electron density is too small for
        ``⟨ε⟩ = n_ε/nₑ`` to be meaningful (an isolated electron equilibrates to the local field).
        Requires a solver that supplied :attr:`reduced_field_td`; ``None`` otherwise.
        """
        if self.reduced_field_td is None:
            return None
        en = np.asarray(self.reduced_field_td, dtype=np.float64)
        eps = np.asarray(self.mean_energy_ev, dtype=np.float64)
        order = np.argsort(en, kind="stable")  # np.interp needs an increasing abscissa
        return np.interp(en_td, en[order], eps[order])

    def to_mean_energy_key(self) -> "CoefficientTable":
        """Return an equivalent table **keyed on mean electron energy ⟨ε⟩** (LMEA closure).

        The solve is naturally on an E/N grid; a fluid model *with* an electron energy equation
        looks up coefficients on the local mean energy instead (Hagelaar & Pitchford 2005 §1).
        Re-keys every sampled array against ``mean_energy_ev`` (sorted ascending), leaving the
        physics untouched. Idempotent if already mean-energy-keyed.
        """
        if self.key_kind == "mean_energy_eV":
            return self
        order = np.argsort(self.mean_energy_ev, kind="stable")
        sorted_energy = np.asarray(self.mean_energy_ev, dtype=np.float64)[order]
        # Saturating closures (notably the analytic stub's Te ceiling) produce repeated
        # energies.  A lookup abscissa must be strictly increasing; retain the first sample
        # at each tie, which is deterministic and preserves the lower-field branch.
        keep = np.concatenate(([True], np.diff(sorted_energy) > 0.0))
        order = order[keep]
        g = sorted_energy[keep]

        def _reorder(arr):
            return None if arr is None else np.asarray(arr, dtype=np.float64)[order]

        return CoefficientTable(
            key_kind="mean_energy_eV",
            key_grid=g,
            mobility_N=_reorder(self.mobility_N),
            diffusion_N=_reorder(self.diffusion_N),
            drift_velocity=_reorder(self.drift_velocity),
            mean_energy_ev=g,  # identity under mean-energy keying
            electron_temp_ev=_reorder(self.electron_temp_ev),
            rate_coeff={rid: np.asarray(arr, dtype=np.float64)[order] for rid, arr in self.rate_coeff.items()},
            mixture=dict(self.mixture),
            gas_temp_k=self.gas_temp_k,
            energy_mobility_N=_reorder(self.energy_mobility_N),
            energy_diffusion_N=_reorder(self.energy_diffusion_N),
            reduced_field_td=_reorder(self.reduced_field_td),
        )

    def to_reduced_field_key(self) -> "CoefficientTable":
        """Return an equivalent table keyed on reduced field E/N in townsend.

        Native mean-energy solvers retain their steady-state E/N pairing in
        :attr:`reduced_field_td`; without it, the requested conversion is undefined.
        """
        if self.key_kind == "EN_Td":
            return self
        if self.reduced_field_td is None:
            raise ValueError("cannot re-key a coefficient table to E/N without reduced_field_td")
        order = np.argsort(self.reduced_field_td, kind="stable")
        sorted_en = np.asarray(self.reduced_field_td, dtype=np.float64)[order]
        keep = np.concatenate(([True], np.diff(sorted_en) > 0.0))
        order = order[keep]
        key = sorted_en[keep]

        def _reorder(arr):
            return None if arr is None else np.asarray(arr, dtype=np.float64)[order]

        return CoefficientTable(
            key_kind="EN_Td",
            key_grid=key,
            mobility_N=_reorder(self.mobility_N),
            diffusion_N=_reorder(self.diffusion_N),
            drift_velocity=_reorder(self.drift_velocity),
            mean_energy_ev=_reorder(self.mean_energy_ev),
            electron_temp_ev=_reorder(self.electron_temp_ev),
            rate_coeff={rid: np.asarray(arr, dtype=np.float64)[order]
                        for rid, arr in self.rate_coeff.items()},
            mixture=dict(self.mixture),
            gas_temp_k=self.gas_temp_k,
            energy_mobility_N=_reorder(self.energy_mobility_N),
            energy_diffusion_N=_reorder(self.energy_diffusion_N),
            reduced_field_td=key,
        )

    def drift(self, key):
        return self._interp(self.drift_velocity, key)

    def mean_energy(self, key):
        return self._interp(self.mean_energy_ev, key)

    def te(self, key):
        return self._interp(self.electron_temp_ev, key)

    def rate(self, reaction_id: str, key):
        arr = self.rate_coeff.get(reaction_id)
        if arr is None:
            return np.zeros_like(np.asarray(key, dtype=np.float64))
        return self._interp(arr, key)

    def reaction_ids(self) -> tuple[str, ...]:
        return tuple(self.rate_coeff)

    # -- on-disk cache (npz) --------------------------------------------------
    def save_npz(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = {
            "key_grid": self.key_grid,
            "mobility_N": self.mobility_N,
            "diffusion_N": self.diffusion_N,
            "drift_velocity": self.drift_velocity,
            "mean_energy_ev": self.mean_energy_ev,
            "electron_temp_ev": self.electron_temp_ev,
        }
        if self.energy_mobility_N is not None:
            arrays["energy_mobility_N"] = self.energy_mobility_N
        if self.energy_diffusion_N is not None:
            arrays["energy_diffusion_N"] = self.energy_diffusion_N
        if self.reduced_field_td is not None:
            arrays["reduced_field_td"] = self.reduced_field_td
        rate_ids = list(self.rate_coeff)
        for i, rid in enumerate(rate_ids):
            arrays[f"rate__{i}"] = self.rate_coeff[rid]
        meta = {
            "key_kind": self.key_kind,
            "mixture": self.mixture,
            "gas_temp_k": self.gas_temp_k,
            "rate_ids": rate_ids,
        }
        arrays["_meta"] = np.frombuffer(json.dumps(meta).encode("utf-8"), dtype=np.uint8)
        np.savez(path, **arrays)
        return path

    @classmethod
    def load_npz(cls, path: Path) -> "CoefficientTable":
        with np.load(path) as data:
            meta = json.loads(bytes(data["_meta"]).decode("utf-8"))
            rate_coeff = {
                rid: data[f"rate__{i}"] for i, rid in enumerate(meta["rate_ids"])
            }
            return cls(
                key_kind=meta["key_kind"],
                key_grid=data["key_grid"],
                mobility_N=data["mobility_N"],
                diffusion_N=data["diffusion_N"],
                drift_velocity=data["drift_velocity"],
                mean_energy_ev=data["mean_energy_ev"],
                electron_temp_ev=data["electron_temp_ev"],
                rate_coeff=rate_coeff,
                mixture=meta["mixture"],
                gas_temp_k=meta["gas_temp_k"],
                energy_mobility_N=data["energy_mobility_N"] if "energy_mobility_N" in data else None,
                energy_diffusion_N=data["energy_diffusion_N"] if "energy_diffusion_N" in data else None,
                reduced_field_td=data["reduced_field_td"] if "reduced_field_td" in data else None,
            )


class BoltzmannSolver(ABC):
    """Produce a :class:`CoefficientTable` for a gas mixture over an E/N grid (offline)."""

    @abstractmethod
    def solve(
        self,
        mixture: dict[str, float],
        en_grid_td: NDArray[np.float64],
        gas_temp_k: float = 300.0,
        *,
        reactions: CompiledNetwork | None = None,
    ) -> CoefficientTable:
        """Compute transport + per-reaction rate coefficients for ``mixture`` on ``en_grid_td``."""

    @property
    @abstractmethod
    def solver_id(self) -> str:
        """Stable identity string for this backend (folds into the LUT cache key)."""


class AnalyticStubSolver(BoltzmannSolver):
    """Closed-form placeholder coefficients so tests need no external binary.

    NOT physical: transport uses simple E/N scalings; each dynamic reaction gets an
    Arrhenius-shaped ``k = k0·exp(−Eth/Te)`` from its threshold. A deterministic stand-in
    that keeps the chemistry + optics path exercised across the F2 sweep.
    """

    def __init__(
        self,
        *,
        mobility_N: float = 1.0e24,
        te_coeff_ev_per_td: float = 0.03,
        te_floor_ev: float = 0.3,
        te_ceiling_ev: float = 40.0,
        default_k: float = 1.0e-15,
    ) -> None:
        self._mobility_N = mobility_N
        self._te_coeff = te_coeff_ev_per_td
        self._te_floor = te_floor_ev
        self._te_ceiling = te_ceiling_ev
        self._default_k = default_k

    @property
    def solver_id(self) -> str:
        return (f"analytic-stub:mu{self._mobility_N:g}:tec{self._te_coeff:g}"
                f":floor{self._te_floor:g}:ceiling{self._te_ceiling:g}"
                f":k{self._default_k:g}")

    def solve(self, mixture, en_grid_td, gas_temp_k=300.0, *, reactions=None) -> CoefficientTable:
        key = np.asarray(en_grid_td, dtype=np.float64)
        te = np.clip(self._te_coeff * key + self._te_floor, self._te_floor, self._te_ceiling)
        mean_energy = 1.5 * te
        mobility_N = np.full_like(key, self._mobility_N)
        diffusion_N = mobility_N * te  # Einstein-like D·N = μ·N · Te
        drift_velocity = 1.0e2 * key  # monotone placeholder [m/s]

        rate_coeff: dict[str, NDArray[np.float64]] = {}
        if reactions is not None:
            for rxn in _dynamic_reactions(reactions):
                eth = _threshold_ev(rxn)
                rate_coeff[rxn.reaction_id] = self._default_k * np.exp(-eth / te)
        return CoefficientTable(
            key_kind="EN_Td",
            key_grid=key,
            mobility_N=mobility_N,
            diffusion_N=diffusion_N,
            drift_velocity=drift_velocity,
            mean_energy_ev=mean_energy,
            electron_temp_ev=te,
            rate_coeff=rate_coeff,
            mixture=dict(mixture),
            gas_temp_k=gas_temp_k,
            # The stub's own Tₑ(E/N) closure doubles as the ⟨ε⟩ ↔ E/N pairing, so the
            # local-field energy fallback works on this backend too.
            reduced_field_td=key,
        )




def _reaction_fingerprint(reactions: CompiledNetwork | None) -> list[dict]:
    """Canonical dynamic-reaction signature for LUT cache invalidation."""
    if reactions is None:
        return []
    return [
        {
            "id": rxn.reaction_id,
            "form": str(rxn.rate_form),
            "reactants": rxn.reactants,
            "products": rxn.products,
            "cross_section": rxn.cross_section,
            "rate_si": rxn.rate_si,
            "params": rxn.params,
        }
        for rxn in _dynamic_reactions(reactions)
    ]


def _mixture_cache_key(
    mixture, en_grid_td, gas_temp_k, solver_id: str, key_kind: str = "EN_Td", *, reactions=None,
) -> str:
    """Stable hash of normalized mixture + E/N grid + gas temperature + solver id + key kind."""
    en = np.asarray(en_grid_td, dtype=np.float64)
    payload = {
        "mixture": {k: round(float(v), 12) for k, v in sorted(mixture.items())},
        # Hash every sample.  Min/max/count collided for grids with different interiors and
        # returned a numerically unrelated table from disk.
        "en": en.tolist(),
        "gas_temp_k": round(float(gas_temp_k), 6),
        "solver": solver_id,
        "key_kind": key_kind,
        "reactions": _reaction_fingerprint(reactions),
        "cache_schema": 2,
    }
    digest = hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()
    return f"lut_{digest[:16]}"


def build_lut(
    mixture: dict[str, float],
    *,
    solver: BoltzmannSolver,
    en_grid_td: NDArray[np.float64],
    gas_temp_k: float = 300.0,
    reactions: CompiledNetwork | None = None,
    cache_dir: Path | str | None = None,
    key_kind: str = "EN_Td",
) -> CoefficientTable:
    """Offline precompute + on-disk cache of a :class:`CoefficientTable` for ``mixture``.

    The F2 sweep calls this once per mixture; the ``.npz`` cache (keyed by mixture + grid +
    temperature + ``solver.solver_id`` + ``key_kind``) is reused across sweep points sharing a
    mixture. ``key_kind="mean_energy_eV"`` re-keys the solve onto the mean electron energy
    (LMEA) — required by the drift-diffusion model's electron energy equation.
    """
    if key_kind not in {"EN_Td", "mean_energy_eV"}:
        raise ValueError(f"unsupported coefficient-table key kind {key_kind!r}")

    def _solve() -> CoefficientTable:
        table = solver.solve(mixture, en_grid_td, gas_temp_k, reactions=reactions)
        if key_kind == "mean_energy_eV":
            return table.to_mean_energy_key()
        return table.to_reduced_field_key()

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_path = cache_dir / (
            f"{_mixture_cache_key(mixture, en_grid_td, gas_temp_k, solver.solver_id, key_kind, reactions=reactions)}.npz"
        )
        if cache_path.exists():
            return CoefficientTable.load_npz(cache_path)
        table = _solve()
        table.save_npz(cache_path)
        return table
    return _solve()


def en_grid_from_spec(spec: tuple[float, float, int]) -> NDArray[np.float64]:
    """Log-spaced reduced-field grid from a ``(min_td, max_td, points)`` spec."""
    lo, hi, n = spec
    return np.geomspace(lo, hi, int(n))
