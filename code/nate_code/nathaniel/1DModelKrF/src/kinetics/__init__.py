"""Per-cell fluid neutral chemistry: dn/dt = ν @ R(n), numba implicit + SciPy reference.

Assemble a :class:`FluidChemistry` from a compiled reaction network, then advance
per-cell neutral densities with :class:`FluidIntegrator` (numba, embarrassingly
parallel over cells). ``FluidChemistry.integrate_reference`` is the SciPy-Radau oracle.
"""

from __future__ import annotations

from .chemistry import FluidChemistry
from .efficiency import LaserEnergetics
from .fluid import NeutralFluid
from .integrator import HAS_NUMBA, FluidIntegrator

__all__ = (
    "FluidChemistry",
    "FluidIntegrator",
    "HAS_NUMBA",
    "LaserEnergetics",
    "NeutralFluid",
)
