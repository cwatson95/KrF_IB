"""Quantitative refinement gates for production discharge/laser runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike


def pulse_fwhm(times: ArrayLike, signal: ArrayLike) -> float:
    """Return a pulse's sampled full width at half maximum in the time-axis units.

    The width spans the first and last finite samples at or above half of the finite
    positive peak. A non-positive or insufficient trace has zero width. Inputs must be
    one-dimensional and have matching lengths.
    """
    t = np.asarray(times, dtype=np.float64)
    y = np.asarray(signal, dtype=np.float64)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
        raise ValueError("times and signal must be one-dimensional arrays of equal length")
    finite = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(finite) < 2:
        return 0.0
    t_finite = t[finite]
    y_finite = y[finite]
    peak = float(np.max(y_finite))
    if peak <= 0.0:
        return 0.0
    above = np.flatnonzero(y_finite >= 0.5 * peak)
    if above.size < 2:
        return 0.0
    return float(t_finite[above[-1]] - t_finite[above[0]])


@dataclass(frozen=True)
class ConvergenceMetrics:
    """Scalar observables compared between two mesh or timestep refinements."""

    cells: int
    dt_s: float
    duration_s: float
    peak_electron_density_m3: float
    deposited_energy_j: float
    output_energy_j: float
    peak_output_power_w: float
    pulse_fwhm_s: float
    chemistry_energy_cost_j: float = 0.0
    chemistry_bound_adjustment_j: float = 0.0
    f2_mole_fraction: float = 0.0
    breakdown_time_s: float = 0.0
    peak_current_a: float = 0.0
    peak_electron_temperature_ev: float = 0.0
    peak_excimer_density_m3: float = 0.0
    excimer_formation_integral_mol_m3: float = 0.0
    excimer_quenching_integral_mol_m3: float = 0.0
    max_energy_closure_residual_j_m3: float = 0.0
    electron_energy_closure: str = "source_coupled"

    def as_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


@dataclass(frozen=True)
class ConvergenceTolerance:
    """Acceptance limits for successive refinements."""

    peak_electron_density_rel: float = 0.10
    deposited_energy_rel: float = 0.05
    output_energy_rel: float = 0.10
    peak_output_power_rel: float = 0.10
    pulse_fwhm_abs_s: float = 0.5e-9
    lasing_energy_floor_j: float = 1.0e-6


@dataclass(frozen=True)
class ConvergenceResult:
    """One coarse/fine comparison and its failed acceptance clauses."""

    passed: bool
    changes: dict[str, float]
    failures: tuple[str, ...]


def _relative_change(coarse: float, fine: float) -> float:
    scale = max(abs(fine), 1.0e-300)
    return abs(coarse - fine) / scale


def compare_refinement(
    coarse: ConvergenceMetrics,
    fine: ConvergenceMetrics,
    tolerance: ConvergenceTolerance = ConvergenceTolerance(),
) -> ConvergenceResult:
    """Apply mesh/timestep refinement gates to two otherwise identical runs."""
    changes = {
        "peak_electron_density_rel": _relative_change(
            coarse.peak_electron_density_m3, fine.peak_electron_density_m3,
        ),
        "deposited_energy_rel": _relative_change(
            coarse.deposited_energy_j, fine.deposited_energy_j,
        ),
    }
    failures: list[str] = []
    for name in ("peak_electron_density_rel", "deposited_energy_rel"):
        limit = getattr(tolerance, name)
        if changes[name] > limit:
            failures.append(f"{name}={changes[name]:.3g} > {limit:.3g}")

    lasing = max(coarse.output_energy_j, fine.output_energy_j) >= tolerance.lasing_energy_floor_j
    if lasing:
        changes["output_energy_rel"] = _relative_change(
            coarse.output_energy_j, fine.output_energy_j,
        )
        changes["peak_output_power_rel"] = _relative_change(
            coarse.peak_output_power_w, fine.peak_output_power_w,
        )
        changes["pulse_fwhm_abs_s"] = abs(coarse.pulse_fwhm_s - fine.pulse_fwhm_s)
        for name in ("output_energy_rel", "peak_output_power_rel", "pulse_fwhm_abs_s"):
            limit = getattr(tolerance, name)
            if changes[name] > limit:
                failures.append(f"{name}={changes[name]:.3g} > {limit:.3g}")
    else:
        changes["output_energy_rel"] = 0.0
        changes["peak_output_power_rel"] = 0.0
        changes["pulse_fwhm_abs_s"] = 0.0

    return ConvergenceResult(not failures, changes, tuple(failures))


__all__ = (
    "ConvergenceMetrics",
    "ConvergenceResult",
    "ConvergenceTolerance",
    "compare_refinement",
    "pulse_fwhm",
)
