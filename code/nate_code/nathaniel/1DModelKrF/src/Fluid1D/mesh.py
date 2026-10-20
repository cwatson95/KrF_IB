"""1D discharge/cavity mesh geometry + boundary helpers.

A :class:`Mesh1D` is the single source of truth for the axial discretisation shared by the neutral
fluid, charged-species fluid, and 1D Poisson solver: ``n_cells`` cells over ``length`` metres along
the single axis, with a transverse ``area_m2`` so absolute densities (m^-3) map to counts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Mesh1D:
    """Uniform 1D mesh: ``n_cells`` cells over ``length`` m, transverse area ``area_m2`` m²."""

    n_cells: int
    length: float
    area_m2: float = 1.0

    @property
    def dx(self) -> float:
        return self.length / self.n_cells

    @property
    def cell_volume(self) -> float:
        return self.dx * self.area_m2

    @property
    def cell_centers(self) -> NDArray[np.float64]:
        """Axial coordinate (m) of each cell centre, shape ``(n_cells,)``."""
        return (np.arange(self.n_cells, dtype=np.float64) + 0.5) * self.dx

    @property
    def edges(self) -> NDArray[np.float64]:
        """Axial coordinate (m) of each cell edge, shape ``(n_cells + 1,)``."""
        return np.arange(self.n_cells + 1, dtype=np.float64) * self.dx

    def cell_of(self, x: NDArray[np.floating] | float) -> NDArray[np.int64]:
        """Flat (periodic) cell index for axial coordinate(s) ``x`` (m)."""
        xa = np.asarray(x, dtype=np.float64)
        return np.floor(xa / self.dx).astype(np.int64) % self.n_cells


