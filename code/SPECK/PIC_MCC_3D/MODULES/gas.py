"""The neutral background -- the collision partner every MCC needs.

A DEPENDENCY THE ROADMAP HAD FOLDED INTO "cross sections" AND IS NOT.  Until
now this package had no gas at all: no composition, no number density, no
temperature.  Every collisional module (M3 electron-neutral and ion-neutral,
M5 emission yields, M9 beam stopping) prices its rates against a neutral
target density, so that object needs to exist once, explicitly, rather than
being re-derived as a loose float in each caller.

Scope, stated up front: the gas here is a STATIC, UNIFORM, IDEAL background.
It is not transported, not depleted by ionization, and not heated by the
discharge.  Those are real effects with real modules behind them elsewhere in
the repo (``pic1d/gas_heating.py`` is the 0-D lumped heating + rarefaction
model, README_LIMITS 2.7), and they are deliberately NOT implemented here --
a 3-D run at the ns timescales this package can afford does not move the
neutrals.  When that assumption is broken (long pulses, strong depletion), it
must be revisited rather than quietly relied on; :meth:`GasMixture.depletion`
exists so a caller can measure how far it has strayed instead of assuming.

Masses are here too, because ion-neutral kinematics (M3) needs the TARGET
mass and the electron-neutral operator needs m_e/M -- and the LXCat files
carry the latter only as a per-process ratio.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .units import AMU, K_B, gas_density_m3

# Atomic/molecular masses [amu] for the gas set this project uses.  These are
# kinematics, not chemistry: the 26-species reaction network stays in
# phantom11.  Values are standard atomic weights (CIAAW), which is ample --
# the mass enters as m_e/M ~ 1e-5 and as reduced masses.
GAS_MASSES_AMU = {
    "Ar": 39.948,
    "Kr": 83.798,
    "He": 4.002602,
    "F2": 37.99681,
    "Ne": 20.1797,
    "N2": 28.0134,
    "O2": 31.9988,
}


@dataclass
class GasMixture:
    """A static, uniform neutral background.

    ``fractions`` maps a gas name to its mole fraction.  They are used AS
    GIVEN, not renormalised -- so a 0.2 % F2 admixture on top of Ar/Kr is
    expressed the way an experimentalist would state it, and a set summing to
    1.002 is not silently rescaled into something nobody asked for.  A sum far
    from 1 is a mistake, though, so the constructor refuses it.
    """

    fractions: dict
    pressure_atm: float = 1.0
    T_gas_K: float = 300.0
    _densities: dict = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        clean = {}
        for name, frac in self.fractions.items():
            f = float(frac)
            if f < 0.0:
                raise ValueError(f"gas fraction for {name!r} is negative ({f})")
            if f > 0.0:
                if name not in GAS_MASSES_AMU:
                    raise KeyError(
                        f"unknown gas {name!r}; known: {sorted(GAS_MASSES_AMU)}")
                clean[name] = f
        if not clean:
            raise ValueError("GasMixture: every fraction is zero")
        total = sum(clean.values())
        if not 0.5 < total < 1.5:
            raise ValueError(
                f"mole fractions sum to {total:.3f}, which is not a mixture -- "
                "give fractions, not percentages or densities")
        if self.pressure_atm <= 0.0:
            raise ValueError(f"pressure_atm = {self.pressure_atm}: must be > 0")
        if self.T_gas_K <= 0.0:
            raise ValueError(f"T_gas_K = {self.T_gas_K}: must be > 0")
        self.fractions = clean
        n_tot = gas_density_m3(self.pressure_atm, self.T_gas_K)
        self._densities = {g: f * n_tot for g, f in clean.items()}

    # ------------------------------------------------------------- densities
    @property
    def total_density_m3(self) -> float:
        """Total neutral number density [m^-3] (ideal gas at p, T)."""
        return float(sum(self._densities.values()))

    def density_of(self, name: str) -> float:
        """Number density of one component [m^-3]; 0 if absent."""
        return float(self._densities.get(name, 0.0))

    @property
    def densities(self) -> dict:
        return dict(self._densities)

    def mass_of(self, name: str) -> float:
        """Target mass [kg] -- what ion-neutral kinematics needs."""
        return GAS_MASSES_AMU[name] * AMU

    @property
    def mean_mass(self) -> float:
        """Fraction-weighted mean neutral mass [kg]."""
        tot = sum(self.fractions.values())
        return float(sum(f * self.mass_of(g)
                         for g, f in self.fractions.items()) / tot)

    # ------------------------------------------------------------- diagnostics
    @property
    def pressure_Pa(self) -> float:
        return float(self.total_density_m3 * K_B * self.T_gas_K)

    def mean_free_path(self, sigma_m2: float) -> float:
        """lambda = 1 / (n sigma) for a given total cross section [m]."""
        n = self.total_density_m3
        if n <= 0.0 or sigma_m2 <= 0.0:
            return np.inf
        return float(1.0 / (n * sigma_m2))

    def thermal_speed(self, name: str) -> float:
        """sqrt(kT/M) of a neutral component [m/s] (per-axis convention)."""
        return float(np.sqrt(K_B * self.T_gas_K / self.mass_of(name)))

    def depletion(self, n_ionized_m3: float) -> float:
        """Fraction of the neutrals converted -- the static-gas assumption's gauge.

        This module holds the neutral density FIXED.  That is defensible while
        the ionized fraction is small; this returns the number so a caller can
        check rather than assume.  Past ~1 % the static background is a stated
        approximation, not a detail.
        """
        n = self.total_density_m3
        return float(max(n_ionized_m3, 0.0) / n) if n > 0 else np.inf

    def scaled_to(self, pressure_atm: float = None, T_gas_K: float = None
                  ) -> "GasMixture":
        """A copy at a different p/T (the composition is unchanged)."""
        return GasMixture(
            dict(self.fractions),
            self.pressure_atm if pressure_atm is None else pressure_atm,
            self.T_gas_K if T_gas_K is None else T_gas_K)

    def summary(self) -> str:
        parts = ", ".join(f"{g} {f:.4g}" for g, f in sorted(self.fractions.items()))
        return (f"{parts} @ {self.pressure_atm:g} atm, {self.T_gas_K:g} K "
                f"-> N = {self.total_density_m3:.3e} m^-3")


# convenience mixtures used by the gates and the presets.  Not "the" operating
# points -- those are a device-design choice recorded in PIC_MCC_3D.txt.
def pure(name: str, pressure_atm: float = 1.0, T_gas_K: float = 300.0
         ) -> GasMixture:
    return GasMixture({name: 1.0}, pressure_atm, T_gas_K)
