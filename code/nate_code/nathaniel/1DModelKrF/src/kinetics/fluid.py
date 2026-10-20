"""Per-cell neutral fluid: number densities of neutral species on the 1D mesh.

The neutral species (Kr, Kr*, F, F2, KrF*, …) are represented as a per-cell density field
``density[cell, species]`` (SI, m⁻³) over the ``N`` cells
of the 1D discharge/cavity axis. The charged-species fluid reads the *local* cell density
to weight electron-impact rates and writes reaction products back into it; the fluid
chemistry integrator (:mod:`kinetics.integrator`) evolves each cell's densities
independently.
"""

from __future__ import annotations

import numpy as np

from particles.species import Species


class NeutralFluid:
    """Cell-indexed neutral-species number densities on a uniform 1D mesh (``N`` cells)."""

    def __init__(
        self,
        species: tuple[Species, ...],
        n_cells: int,
        dx: float,
        area_m2: float = 1.0,
    ) -> None:
        if any(s.is_charged or s.is_photon for s in species):
            raise ValueError("NeutralFluid holds only neutral species (no ions/photons)")
        self.species = species
        self.names = tuple(s.name for s in species)
        self.column = {name: i for i, name in enumerate(self.names)}
        self.n_cells = int(n_cells)
        self.dx = dx
        self.cell_volume = dx * area_m2
        # (n_cells, n_species) so cells are the parallel axis for the fluid integrator.
        self.density = np.zeros((self.n_cells, len(species)), dtype=np.float64)

    @property
    def n_species(self) -> int:
        return len(self.species)

    def seed_uniform(self, name: str, density: float) -> None:
        """Set a uniform density (m⁻³) of ``name`` in every cell."""
        self.density[:, self.column[name]] = density

    def total_number(self, name: str) -> float:
        """Total real-particle count of ``name`` across the domain (density·volume summed)."""
        return float(self.density[:, self.column[name]].sum()) * self.cell_volume
