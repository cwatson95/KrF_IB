"""Convert tagged reaction rates and cross-sections to internal SI.

Reaction rows retain their published cgs or SI unit tag. The loader normalizes each
value once so the engine runs in m^-3, m^2, and s.
"""

from __future__ import annotations

# Multiplicative factor to convert a tagged unit to SI.
#   density in cm^-3 -> m^-3 is x1e6, so an n-body rate constant (cm^{3(n-1)}/s)
#   converts by x1e-6^(n-1): 2-body cm^3/s -> x1e-6, 3-body cm^6/s -> x1e-12.
#   cross-section cm^2 -> m^2 is x1e-4. "eV^k" factors are dimensionless here (T_e in eV).
_UNIT_TO_SI: dict[str, float] = {
    "cm3/s": 1e-6,
    "cm^3/s": 1e-6,
    "cm6/s": 1e-12,
    "cm^6/s": 1e-12,
    "cm2": 1e-4,
    "cm^2": 1e-4,
    "m3/s": 1.0,
    "m6/s": 1.0,
    "m2": 1.0,
    "1/s": 1.0,
    "s^-1": 1.0,
    "s": 1.0,  # a lifetime tau (converted trivially)
    "cm3 ev^0.5/s": 1e-6,  # e.g. K0*Te^-1/2 dissociative recombination
    "cm3 ev^-0.5/s": 1e-6,
    "m3 ev^0.5/s": 1.0,
}


def to_si(value: float, units: str) -> float:
    """Convert a tagged cgs or SI ``value`` to SI."""
    key = " ".join(units.lower().split())
    if key not in _UNIT_TO_SI:
        raise KeyError(
            f"unknown rate/cross-section unit {units!r}; known: {sorted(_UNIT_TO_SI)}"
        )
    return value * _UNIT_TO_SI[key]
