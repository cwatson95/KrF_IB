"""Small, auditable electron-energy reconciliation helpers.

The production source-coupled path does not use an interval reconciliation: its accepted
reaction extent and energy debit come from one augmented source stage.  The diagnostic
continuous path does, so its signed correction and any subsequent coefficient-table bound
adjustment live here instead of being open-coded in the orchestrator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ReconciledEnergy:
    """Per-cell result of the diagnostic continuous-cooling reconciliation."""

    final_j_m3: NDArray[np.float64]
    signed_reconciliation_j_m3: NDArray[np.float64]
    bound_adjustment_j_m3: NDArray[np.float64]


def reconcile_continuous_energy(
    transport_j_m3: NDArray[np.float64],
    prepaid_sink_j_m3: NDArray[np.float64],
    accepted_cost_j_m3: NDArray[np.float64],
    electron_density_m3: NDArray[np.float64],
    energy_bounds_ev: tuple[float, float],
    *,
    ev_to_j: float,
) -> ReconciledEnergy:
    """Apply ``+S_PDE-C_accepted`` and expose, rather than hide, energy bounds.

    ``signed_reconciliation_j_m3`` is positive for a refund and negative for an additional
    debit.  Bounds are a separate numerical ledger term and never change that sign.
    """
    transport = np.asarray(transport_j_m3, dtype=np.float64)
    prepaid = np.asarray(prepaid_sink_j_m3, dtype=np.float64)
    accepted = np.asarray(accepted_cost_j_m3, dtype=np.float64)
    electron_density = np.maximum(np.asarray(electron_density_m3, dtype=np.float64), 0.0)
    signed = prepaid - accepted
    raw = transport + signed
    lower = energy_bounds_ev[0] * electron_density * ev_to_j
    upper = energy_bounds_ev[1] * electron_density * ev_to_j
    final = np.minimum(np.maximum(raw, lower), upper)
    return ReconciledEnergy(
        final_j_m3=final,
        signed_reconciliation_j_m3=signed,
        bound_adjustment_j_m3=final - raw,
    )


__all__ = ("ReconciledEnergy", "reconcile_continuous_energy")
