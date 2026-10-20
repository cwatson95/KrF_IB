"""Enums for the reaction network: rate forms and execution layers."""

from __future__ import annotations

from enum import StrEnum


class RateForm(StrEnum):
    """How a reaction's rate constant is specified (and evaluated)."""

    CONSTANT = "constant"  # tagged two-/three-body coefficient; SI-converted on load
    TE_DEPENDENT = "te_dependent"  # k = f(T_e): K0*exp(-E0/kTe), K0*Te^p, ...
    RADIATIVE = "radiative"  # unimolecular decay, k = 1/tau (s^-1)
    ELECTRON_IMPACT = "electron_impact"  # cross-section-driven coefficient table
    PHOTO = "photo"  # photon absorption/ionization from a bare sigma_m2 cross-section
    STIMULATED = "stimulated"  # stimulated emission: sigma_eff (cm^2), gain channel


class Layer(StrEnum):
    """Which physics layer executes a reaction (assigned at compile time by charge)."""

    KINETIC = "kinetic"  # at least one charged reactant -> combined plasma chemistry
    PHOTON = "photon"  # >=1 photon reactant -> scalar cavity/optics layer
    FLUID = "fluid"  # all-neutral reactants -> per-cell mass-action rate equation


class ReactionKind(StrEnum):
    """Cross-section process kind, mapped from Cross-Section Browser ``type``/``category``.

    The cross-section loader tags each dataset with its process kind, which the Boltzmann rate
    builder consumes.
    """

    ELASTIC = "elastic"
    EXCITATION = "excitation"
    VIBRATIONAL = "vibrational"
    IONIZATION = "ionization"
    ATTACHMENT = "attachment"
    RECOMBINATION = "recombination"
    OTHER = "other"


def reaction_kind_from_metadata(type_name: str, category: str | None = None) -> ReactionKind:
    """Map JSON ``type`` / ``category`` strings to :class:`ReactionKind`."""
    key = (type_name or category or "").strip().lower()
    if "ioniz" in key:
        return ReactionKind.IONIZATION
    if "vibrat" in key:
        return ReactionKind.VIBRATIONAL
    if "excit" in key:
        return ReactionKind.EXCITATION
    if "elastic" in key or "effective" in key:
        # LXCat "EFFECTIVE" is the total momentum-transfer cross-section (elastic + inelastic
        # folded in) and is used exactly where an elastic momentum transfer would be; without
        # this it would fall through to OTHER and be silently dropped by every consumer.
        return ReactionKind.ELASTIC
    if "attach" in key:
        return ReactionKind.ATTACHMENT
    if "recomb" in key:
        return ReactionKind.RECOMBINATION
    return ReactionKind.OTHER
