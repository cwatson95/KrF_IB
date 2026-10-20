"""Canonical electron-energy transfer metadata for compiled reaction networks.

Every electron-energy closure consumes this module's table.  Rate expressions remain owned by
the coefficient/chemistry solvers; this table only answers how much electron energy is transferred
by one accepted reaction event.  An explicit ``energy_loss_eV`` wins over a cross-section threshold
because an analytic activation energy is not necessarily the physical energy transfer.
"""

from __future__ import annotations

from .network import Reaction


def reaction_energy_loss_ev(reaction: Reaction) -> float:
    """Return the declared energy transfer per event in eV, or zero when absent."""
    params = reaction.params or {}
    for key in ("energy_loss_eV", "energy_loss_ev", "threshold_eV", "threshold_ev"):
        if key in params:
            try:
                return max(float(params[key]), 0.0)
            except (TypeError, ValueError):
                return 0.0
    return 0.0


__all__ = ("reaction_energy_loss_ev",)
