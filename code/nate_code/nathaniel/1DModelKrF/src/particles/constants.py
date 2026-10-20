"""SI physical constants for the 1D fluid plasma model."""

E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
M_U = 1.66053906660e-27
M_H = 1.00784 * M_U
M_KR = 83.798 * M_U
M_F = 18.998403163 * M_U
M_HE = 4.002602 * M_U
M_AR = 39.948 * M_U
M_NE = 20.1797 * M_U
M_XE = 131.293 * M_U
K_B = 1.380649e-23
EV_TO_J = E_CHARGE
C_LIGHT = 299_792_458.0  # m/s, for photon (hν) kinematics and photo-process rates
EPS0 = 8.8541878128e-12  # F/m vacuum permittivity (seeds the Poisson solve)
AVOGADRO = 6.02214076e23  # mol^-1
ATM_PA = 101_325.0  # standard atmosphere in pascals
TD = 1.0e-21  # V*m^2 per townsend

# Atomic masses (kg) keyed by element symbol, used to compose molecular species
# (e.g. F2, Kr2, KrF) from their constituents. "e" lets electrons share the lookup.
ATOMIC_MASS: dict[str, float] = {
    "e": M_E,
    "H": M_H,
    "F": M_F,
    "He": M_HE,
    "Ne": M_NE,
    "Ar": M_AR,
    "Kr": M_KR,
    "Xe": M_XE,
}


def molecular_mass(*symbols: str) -> float:
    """Sum atomic masses (kg) for a molecule given its element symbols.

    ``molecular_mass("Kr", "F")`` -> mass of KrF; ``molecular_mass("F", "F")`` -> F2.
    """
    try:
        return sum(ATOMIC_MASS[symbol] for symbol in symbols)
    except KeyError as exc:  # pragma: no cover - guards catalog typos
        raise KeyError(f"unknown element symbol {exc.args[0]!r}; add it to ATOMIC_MASS") from exc
