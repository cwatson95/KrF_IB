"""Species metadata for the KrF and ArF low-temperature plasma systems.

:class:`Species` represents electrons, neutral ground and excited states, atomic and
molecular ions, excimers, and photon modes with one immutable metadata record. Runtime
fluids store cell-density arrays keyed by these records; reaction, transport, chemistry,
and optics layers share the same registry.

Orthogonal physical attributes (rather than a deep type hierarchy) keep the model
flexible:

- ``charge_number`` — signed integer charge; ``0`` for neutrals, ``-1`` for e- and F-.
- ``constituents`` — atomic symbols making up the species; ``len >= 2`` is a molecule.
- ``excitation_energy_ev`` — electronic energy above the *ground state of the same
  composition and charge*; ``0`` for every ground species, ``> 0`` for Kr*, KrF*, ...
- ``category`` — a coarse label for grouping/diagnostics, derivable info made explicit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

from .constants import E_CHARGE, molecular_mass

# --- Representative electronic-state energies (eV above the ground state) -----------
# Kr I 5s manifold (4p^5 5s): metastable 3P2 ~9.915, 3P1 ~10.033, 1P1 ~10.644 eV.
# A single lumped "Kr*" effective state takes the metastable energy.
KR_5S_ENERGY_EV = 9.92
# Kr I 5p manifold (4p^5 5p): ~12.0-12.9 eV; lumped "Kr**" effective state.
KR_5P_ENERGY_EV = 12.3
# KrF* (B -> X) bound-free emission near 248 nm carries ~5.0 eV; used as the excimer's
# stored electronic energy above the (dissociative) ground KrF configuration.
KRF_EMISSION_EV = 5.0

# Planck relation constant hc in eV·nm, so photon_energy_ev(λ_nm) = _HC_EV_NM / λ_nm.
_HC_EV_NM = 1239.841984


def photon_energy_ev(wavelength_nm: float) -> float:
    """Photon energy (eV) for a vacuum wavelength in nm (E = hc/λ)."""
    return _HC_EV_NM / wavelength_nm


class SpeciesCategory(StrEnum):
    """Coarse classification used for grouping, diagnostics, and field coupling.

    Charge is tracked separately on :attr:`Species.charge_number`, so a charged
    molecule (Kr2+) is a :attr:`DIMER` with non-zero charge rather than a distinct
    category. Use the :class:`Species` helper properties for derived questions.
    """

    ELECTRON = "electron"  # e-
    ION = "ion"  # charged atomic species (Kr+, F-)
    NEUTRAL = "neutral"  # ground-state neutral atom (Kr, F)
    EXCITED = "excited"  # electronically excited neutral atom (Kr*, Kr**)
    DIMER = "dimer"  # ground-state diatomic molecule (F2, Kr2, Kr2+)
    EXCIMER = "excimer"  # electronically excited diatomic/polyatomic (KrF*, Kr2*, Kr2F*)
    PHOTON = "photon"  # radiation quantum hν (a laser/optical mode); no atoms, no charge


@dataclass(frozen=True)
class Species:
    """One physical species tracked (or referenced) by the plasma model.

    A :class:`Species` is immutable identity and physics metadata. Runtime state is held
    as cell-density arrays in the neutral and charged fluids.
    """

    name: str
    mass: float
    charge_number: int = 0
    category: SpeciesCategory = SpeciesCategory.NEUTRAL
    excitation_energy_ev: float = 0.0
    constituents: tuple[str, ...] = ()
    metastable: bool = False
    relaxes_to: str | None = None  # name of the de-excitation product (same composition)
    wavelength_nm: float | None = None  # set only for PHOTON modes (hν); None otherwise

    @property
    def charge(self) -> float:
        """Physical charge in coulombs (``charge_number * e``)."""
        return self.charge_number * E_CHARGE

    @property
    def is_charged(self) -> bool:
        return self.charge_number != 0

    @property
    def is_photon(self) -> bool:
        """A radiation quantum (optical mode) — carries energy but no atoms/charge."""
        return self.category == SpeciesCategory.PHOTON

    @property
    def is_field_coupled(self) -> bool:
        """Whether the species responds to the electrostatic field."""
        return self.charge_number != 0

    @property
    def is_molecule(self) -> bool:
        return len(self.constituents) >= 2

    @property
    def is_excited(self) -> bool:
        return self.excitation_energy_ev > 0.0

@dataclass(frozen=True)
class IonSpecies(Species):
    """A :class:`Species` whose defaults describe a singly charged ion."""

    charge_number: int = 1
    category: SpeciesCategory = SpeciesCategory.ION


@dataclass
class SpeciesRegistry:
    """Lookup table of :class:`Species` by name, for assembling a reaction system."""

    _by_name: dict[str, Species] = field(default_factory=dict)

    @classmethod
    def from_species(cls, species: tuple[Species, ...]) -> SpeciesRegistry:
        """Create a registry in tuple order, rejecting duplicate names."""
        registry = cls()
        for spec in species:
            registry.add(spec)
        return registry

    def add(self, species: Species) -> None:
        """Register one species, preserving insertion order."""
        if species.name in self._by_name:
            raise ValueError(f"duplicate species name in registry: {species.name!r}")
        self._by_name[species.name] = species

    def get(self, name: str) -> Species:
        """Return the named species or raise a registry-aware :class:`KeyError`."""
        try:
            return self._by_name[name]
        except KeyError as exc:
            raise KeyError(f"unknown species {name!r}; registered: {sorted(self._by_name)}") from exc

    def __contains__(self, name: object) -> bool:
        return name in self._by_name

    def __iter__(self):
        return iter(self._by_name.values())

    def __len__(self) -> int:
        return len(self._by_name)

    @property
    def names(self) -> tuple[str, ...]:
        """Registered species names in insertion order."""
        return tuple(self._by_name)

    def by_category(self, category: SpeciesCategory) -> tuple[Species, ...]:
        """Return species carrying exactly ``category`` in insertion order."""
        return tuple(s for s in self._by_name.values() if s.category == category)

    def charged(self) -> tuple[Species, ...]:
        """Return every species with nonzero charge number."""
        return tuple(s for s in self._by_name.values() if s.is_charged)

    def excited(self) -> tuple[Species, ...]:
        """Return species with positive stored energy, including photon modes."""
        return tuple(s for s in self._by_name.values() if s.is_excited)


# --- Kr/F excimer system catalog ----------------------------------------------------
# The species relevant to a KrF* excimer plasma. Excited states and excimers start
# with zero population and are produced during a run (e.g. Kr* + F2 -> KrF* + F); they
# are defined here so the model has a stable identity to create them under.

ELECTRON = Species("e-", molecular_mass("e"), -1, SpeciesCategory.ELECTRON, constituents=("e",))

KR = Species("Kr", molecular_mass("Kr"), 0, SpeciesCategory.NEUTRAL, constituents=("Kr",))
KR_STAR = Species(
    "Kr*",
    molecular_mass("Kr"),
    0,
    SpeciesCategory.EXCITED,
    excitation_energy_ev=KR_5S_ENERGY_EV,
    constituents=("Kr",),
    metastable=True,
    relaxes_to="Kr",
)
KR_2STAR = Species(
    "Kr**",
    molecular_mass("Kr"),
    0,
    SpeciesCategory.EXCITED,
    excitation_energy_ev=KR_5P_ENERGY_EV,
    constituents=("Kr",),
    relaxes_to="Kr",
)
KR_PLUS = IonSpecies("Kr+", molecular_mass("Kr"), constituents=("Kr",))
KR2_PLUS = Species(
    "Kr2+",
    molecular_mass("Kr", "Kr"),
    1,
    SpeciesCategory.DIMER,
    constituents=("Kr", "Kr"),
)

F = Species("F", molecular_mass("F"), 0, SpeciesCategory.NEUTRAL, constituents=("F",))
F2 = Species("F2", molecular_mass("F", "F"), 0, SpeciesCategory.DIMER, constituents=("F", "F"))
F_MINUS = Species("F-", molecular_mass("F"), -1, SpeciesCategory.ION, constituents=("F",))

KRF_EXCIMER = Species(
    "KrF*",
    molecular_mass("Kr", "F"),
    0,
    SpeciesCategory.EXCIMER,
    excitation_energy_ev=KRF_EMISSION_EV,
    constituents=("Kr", "F"),
)
KR2_EXCIMER = Species(
    "Kr2*",
    molecular_mass("Kr", "Kr"),
    0,
    SpeciesCategory.EXCIMER,
    excitation_energy_ev=KR_5S_ENERGY_EV,
    constituents=("Kr", "Kr"),
)

KR2F_EXCIMER = Species(
    "Kr2F*",
    molecular_mass("Kr", "Kr", "F"),
    0,
    SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(415.0),  # broadband triatomic emission ~415 nm
    constituents=("Kr", "Kr", "F"),
)
F2_STAR = Species(
    "F2*",
    molecular_mass("F", "F"),
    0,
    SpeciesCategory.EXCIMER,
    excitation_energy_ev=4.0,  # lumped F2 excited (A^1Π / a^3Π) manifold
    constituents=("F", "F"),
)
# The two resolved F2 excited states of Chen et al. Table III (R11-R14), distinct from the
# lumped F2* above, plus the F2 molecular cation (R15/R22).
F2_A3 = Species(
    "F2(a3)", molecular_mass("F", "F"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=3.16, constituents=("F", "F"),  # a^3Π_u
)
F2_A1 = Species(
    "F2(A1)", molecular_mass("F", "F"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=4.34, constituents=("F", "F"),  # A^1Π_u
)
F2_PLUS = Species(
    "F2+", molecular_mass("F", "F"), 1, SpeciesCategory.DIMER, constituents=("F", "F"),
)

# --- Buffer / rare-gas partners (He, Ne, Ar, Xe) and their excited/ionic states -----
# Excitation energies use representative metastable levels (eV); thresholds for the
# actual reactions come from the reaction data, so these are diagnostic metadata.
HE = Species("He", molecular_mass("He"), 0, SpeciesCategory.NEUTRAL, constituents=("He",))
HE_STAR = Species(
    "He*", molecular_mass("He"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=19.82, constituents=("He",), metastable=True, relaxes_to="He",
)
HE_PLUS = IonSpecies("He+", molecular_mass("He"), constituents=("He",))

NE = Species("Ne", molecular_mass("Ne"), 0, SpeciesCategory.NEUTRAL, constituents=("Ne",))
NE_STAR = Species(
    "Ne*", molecular_mass("Ne"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=16.62, constituents=("Ne",), metastable=True, relaxes_to="Ne",
)
NE_2STAR = Species(
    "Ne**", molecular_mass("Ne"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=18.38, constituents=("Ne",), relaxes_to="Ne",  # higher Ne manifold (Chen ArF R2/R6)
)
NE_PLUS = IonSpecies("Ne+", molecular_mass("Ne"), constituents=("Ne",))
NE2_PLUS = Species("Ne2+", molecular_mass("Ne", "Ne"), 1, SpeciesCategory.DIMER, constituents=("Ne", "Ne"))
NE2_STAR = Species(
    "Ne2*", molecular_mass("Ne", "Ne"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(80.0), constituents=("Ne", "Ne"),
)

AR = Species("Ar", molecular_mass("Ar"), 0, SpeciesCategory.NEUTRAL, constituents=("Ar",))
AR_STAR = Species(
    "Ar*", molecular_mass("Ar"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=11.55, constituents=("Ar",), metastable=True, relaxes_to="Ar",
)
AR_2STAR = Species(
    "Ar**", molecular_mass("Ar"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=13.2, constituents=("Ar",), relaxes_to="Ar",
)
AR_PLUS = IonSpecies("Ar+", molecular_mass("Ar"), constituents=("Ar",))
AR2_PLUS = Species("Ar2+", molecular_mass("Ar", "Ar"), 1, SpeciesCategory.DIMER, constituents=("Ar", "Ar"))
AR2_STAR = Species(
    "Ar2*", molecular_mass("Ar", "Ar"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(126.0), constituents=("Ar", "Ar"),
)

XE = Species("Xe", molecular_mass("Xe"), 0, SpeciesCategory.NEUTRAL, constituents=("Xe",))
XE_STAR = Species(
    "Xe*", molecular_mass("Xe"), 0, SpeciesCategory.EXCITED,
    excitation_energy_ev=8.32, constituents=("Xe",), metastable=True, relaxes_to="Xe",
)
XE_PLUS = IonSpecies("Xe+", molecular_mass("Xe"), constituents=("Xe",))

# --- Rare-gas-halide excimers and mixed dimer ions (Ar/Kr systems) ------------------
ARF_EXCIMER = Species(
    "ArF*", molecular_mass("Ar", "F"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(193.0), constituents=("Ar", "F"),
)
AR2F_EXCIMER = Species(
    "Ar2F*", molecular_mass("Ar", "Ar", "F"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(290.0), constituents=("Ar", "Ar", "F"),
)
ARKR_PLUS = Species(
    "ArKr+", molecular_mass("Ar", "Kr"), 1, SpeciesCategory.DIMER, constituents=("Ar", "Kr"),
)
ARKRF_EXCIMER = Species(
    "ArKrF*", molecular_mass("Ar", "Kr", "F"), 0, SpeciesCategory.EXCIMER,
    excitation_energy_ev=photon_energy_ev(248.0), constituents=("Ar", "Kr", "F"),
)

# Lumped, composition-agnostic positive ion used to close charge for reactions whose
# ionized partner is otherwise untracked (e.g. Penning He* + M -> He + M+ + e). It
# carries +1 charge but no ``constituents``, so it participates in charge balance while
# being (intentionally) transparent to atom-conservation bookkeeping.
M_PLUS = Species("M+", molecular_mass("Ar"), 1, SpeciesCategory.ION, constituents=())

# --- Photon (optical) modes: hν quanta tracked as their own species -----------------
# Laser modes plus the non-lasing continua used as radiative sinks in the reduced chemistry.
# They are massless and carry energy hc/λ; only hv248/hv193 are cavity modes.
PHOTON_80 = Species(
    "hv80", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(80.0), wavelength_nm=80.0,
)
PHOTON_147 = Species(
    "hv147", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(147.0), wavelength_nm=147.0,
)
PHOTON_248 = Species(
    "hv248", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(248.0), wavelength_nm=248.0,
)
PHOTON_193 = Species(
    "hv193", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(193.0), wavelength_nm=193.0,
)
PHOTON_95 = Species(
    "hv95", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(95.0), wavelength_nm=95.0,
)
PHOTON_415 = Species(
    "hv415", 0.0, 0, SpeciesCategory.PHOTON,
    excitation_energy_ev=photon_energy_ev(415.0), wavelength_nm=415.0,
)

KRF_SYSTEM_SPECIES: tuple[Species, ...] = (
    ELECTRON,
    KR,
    KR_STAR,
    KR_2STAR,
    KR_PLUS,
    KR2_PLUS,
    F,
    F2,
    F_MINUS,
    F2_A3,
    F2_PLUS,
    KRF_EXCIMER,
    KR2_EXCIMER,
    KR2F_EXCIMER,
    NE,
    NE_STAR,
    NE_2STAR,
    NE_PLUS,
    NE2_PLUS,
    NE2_STAR,
    PHOTON_80,
    PHOTON_95,
    PHOTON_147,
    PHOTON_248,
    PHOTON_415,
)

# Full KrF/ArF excimer-laser species set (Kr + Ar + He + Ne + Xe systems, excimers,
# mixed-dimer ions, the lumped M+, and the photon modes). Superset of KRF_SYSTEM_SPECIES.
EXCIMER_SYSTEM_SPECIES: tuple[Species, ...] = (
    *KRF_SYSTEM_SPECIES,
    F2_STAR, F2_A1,
    HE, HE_STAR, HE_PLUS,
    AR, AR_STAR, AR_2STAR, AR_PLUS, AR2_PLUS, AR2_STAR,
    XE, XE_STAR, XE_PLUS,
    ARF_EXCIMER, AR2F_EXCIMER, ARKR_PLUS, ARKRF_EXCIMER,
    M_PLUS,
    PHOTON_193,
)


def krf_system_registry() -> SpeciesRegistry:
    """A fresh registry for the reduced Ne/Kr/F2 KrF* reaction system."""
    return SpeciesRegistry.from_species(KRF_SYSTEM_SPECIES)


def excimer_system_registry() -> SpeciesRegistry:
    """A fresh :class:`SpeciesRegistry` with the full KrF/ArF excimer-laser species set."""
    return SpeciesRegistry.from_species(EXCIMER_SYSTEM_SPECIES)
