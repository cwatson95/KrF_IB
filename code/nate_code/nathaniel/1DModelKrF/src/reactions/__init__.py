"""Reaction network: declarative reactions -> layer-routed, SI, mass-action matrices.

The single source of truth for the KrF/ArF chemistry, consumed by the charged-species
fluid (electron-impact / charged reactants), the scalar cavity model (photon
reactants), and the per-cell neutral fluid integrator (all-neutral reactants). Declare
with :meth:`ReactionNetwork.add` or load a TOML file with :func:`load_reaction_network`,
then :meth:`ReactionNetwork.compile`. Cross-section datasets (LXCat-derived JSON+CSV)
are loaded here too — see :mod:`reactions.cross_sections`.
"""

from __future__ import annotations

from .cross_sections import (
    CrossSectionCatalog,
    CrossSectionImport,
    CrossSectionImportReport,
    CrossSectionRecord,
    TotalCrossSectionInterpolator,
    build_cross_section_index,
    import_cross_sections,
    interpolate_sigma_linear,
    load_cross_section,
    load_cross_section_directory,
    match_record_to_reaction,
    resolve_cross_section_path,
)
from .data import load_reaction_network
from .energy_transfer import reaction_energy_loss_ev
from .ion_transport import IonMobilityRecord, IonMobilityTable, load_ion_mobility
from .network import (
    CompiledNetwork,
    Reaction,
    ReactionNetwork,
    photo_absorbers,
    photo_channels,
    radiative_channels,
    stimulated_channel,
    stimulated_sigma,
)
from .parsing import ALIASES, normalize_name, parse_side
from .types import Layer, RateForm, ReactionKind, reaction_kind_from_metadata
from .units import to_si

__all__ = (
    "ALIASES",
    "CompiledNetwork",
    "CrossSectionCatalog",
    "CrossSectionImport",
    "CrossSectionImportReport",
    "CrossSectionRecord",
    "IonMobilityRecord",
    "IonMobilityTable",
    "Layer",
    "Reaction",
    "ReactionKind",
    "ReactionNetwork",
    "RateForm",
    "TotalCrossSectionInterpolator",
    "build_cross_section_index",
    "import_cross_sections",
    "interpolate_sigma_linear",
    "load_ion_mobility",
    "load_cross_section",
    "load_cross_section_directory",
    "load_reaction_network",
    "reaction_energy_loss_ev",
    "match_record_to_reaction",
    "normalize_name",
    "parse_side",
    "photo_absorbers",
    "photo_channels",
    "stimulated_channel",
    "radiative_channels",
    "reaction_kind_from_metadata",
    "resolve_cross_section_path",
    "stimulated_sigma",
    "to_si",
)
