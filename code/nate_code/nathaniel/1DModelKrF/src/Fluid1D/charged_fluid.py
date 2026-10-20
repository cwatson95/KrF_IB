"""Per-cell charged-species fluid: densities + electron energy on the 1D mesh.

The charged species (e⁻, Kr⁺, Kr₂⁺, F⁻, …) are represented as per-cell number-density vectors
``density[cell, species]`` (SI, m⁻³) over the ``N`` cells of the 1D axis, plus a per-cell
**electron mean-energy** field ``electron_mean_energy_ev[cell]`` (eV) from which
``Tₑ = ⅔⟨ε⟩``. Electron-distribution diagnostics reconstruct a Maxwellian from this mean
energy; the fluid does not store a resolved cell-by-energy EEDF.

``ChargedFluid`` mirrors :class:`kinetics.fluid.NeutralFluid` so the two share cell geometry
and the same ``FluidIntegrator`` can advance a combined charged+neutral chemistry.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from particles.species import Species, SpeciesCategory


class ChargedFluid:
    """Cell-indexed charged-species number densities + electron energy on the 1D mesh."""

    def __init__(
        self,
        species: tuple[Species, ...],
        n_cells: int,
        dx: float,
        area_m2: float = 1.0,
    ) -> None:
        if any(not s.is_charged for s in species):
            raise ValueError("ChargedFluid holds only charged species (electrons + ions)")
        self.species = species
        self.names = tuple(s.name for s in species)
        self.column = {name: i for i, name in enumerate(self.names)}
        self.n_cells = int(n_cells)
        self.dx = dx
        self.area_m2 = area_m2
        self.cell_volume = dx * area_m2
        self.density = np.zeros((self.n_cells, len(species)), dtype=np.float64)
        # Per-cell electron mean energy <eps> (eV); Te = 2/3 <eps>.
        self.electron_mean_energy_ev = np.zeros(self.n_cells, dtype=np.float64)
        # Electron column (by category, fallback to the conventional "e-" name).
        self._electron_col = next(
            (i for i, s in enumerate(species) if s.category is SpeciesCategory.ELECTRON),
            self.column.get("e-"),
        )
        # Per-species charge (C), used by net_charge_density.
        self._charges = np.array([s.charge for s in species], dtype=np.float64)

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def electron_col(self) -> int:
        if self._electron_col is None:
            raise KeyError("ChargedFluid has no electron species")
        return self._electron_col

    def density_of(self, name: str) -> NDArray[np.float64]:
        """1D vector n_s(x) (m⁻³) for charged species ``name`` (a view into the column)."""
        return self.density[:, self.column[name]]

    def electron_density(self) -> NDArray[np.float64]:
        """Per-cell electron density nₑ(x) (m⁻³)."""
        return self.density[:, self.electron_col]

    def electron_mean_energy(self) -> NDArray[np.float64]:
        """Per-cell electron mean energy ⟨ε⟩(x) (eV)."""
        return self.electron_mean_energy_ev

    def electron_temperature_ev(self) -> NDArray[np.float64]:
        """Per-cell electron temperature Tₑ(x) = ⅔⟨ε⟩(x) (eV)."""
        return (2.0 / 3.0) * self.electron_mean_energy_ev

    def seed_uniform(self, name: str, density: float) -> None:
        """Set a uniform density (m⁻³) of ``name`` in every cell."""
        self.density[:, self.column[name]] = density

    def seed_profile(self, name: str, values: NDArray[np.floating]) -> None:
        """Set a per-cell density profile (m⁻³, length ``n_cells``) for ``name``."""
        col = self.column[name]
        arr = np.asarray(values, dtype=np.float64)
        if arr.shape != (self.n_cells,):
            raise ValueError(
                f"seed_profile for {name!r} expects shape ({self.n_cells},), got {arr.shape}"
            )
        self.density[:, col] = arr

    def seed_electron_energy(self, mean_energy_ev: float) -> None:
        """Set a uniform electron mean energy (eV) in every cell."""
        self.electron_mean_energy_ev[:] = mean_energy_ev

    def total_number(self, name: str) -> float:
        """Total real-particle count of ``name`` across the domain (density·volume summed)."""
        return float(self.density[:, self.column[name]].sum()) * self.cell_volume

    def net_charge_density(self) -> NDArray[np.float64]:
        """Per-cell net charge density ρ(x) = Σ_s q_s n_s(x) (C/m³) — the Poisson source."""
        return self.density @ self._charges
