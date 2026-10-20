"""The reaction network: the single source of truth over all physics layers.

``ReactionNetwork.add(reactants, products, r=... | sigma=... | tau_s=...)`` is the
user-facing (Keras-like) declaration; ``.compile()`` resolves species against a
:class:`~particles.species.SpeciesRegistry`, **routes each reaction to a layer by
reactant charge** (photon reactant -> PHOTON, any charged reactant -> KINETIC, else
FLUID), and builds the fluid mass-action matrices:

- ``orders`` O (S x R_fluid): reactant multiplicity of species s in fluid reaction r.
- ``stoichiometry`` Sigma (S x R_fluid): net (products - reactants) of species s in r.

The fluid RHS is then ``dn/dt = Sigma @ (k * prod)`` with ``prod[r] = prod_s n_s**O[s,r]``
(assembled in the kinetics integrator). Electron-driven reactions are handed to coefficient/
chemistry consumers and photon reactions to the cavity model. Rates are converted to SI at add time.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np

from particles.species import SpeciesRegistry, excimer_system_registry

from .parsing import parse_side
from .types import Layer, RateForm
from .units import to_si

Term = tuple[int, str]  # (multiplicity, species name)


@dataclass(frozen=True)
class Reaction:
    """One declared reaction with SI-normalized rate metadata (layer assigned at compile)."""

    reactants: tuple[Term, ...]
    products: tuple[Term, ...]
    rate_form: RateForm
    rate_si: float | None = None  # CONSTANT/PHOTO/STIMULATED k (SI); RADIATIVE 1/tau; TE_DEPENDENT K0
    cross_section: str | None = None  # dataset process-id/path for ELECTRON_IMPACT
    params: dict = field(default_factory=dict)  # E0_eV, te_power, threshold_eV, sigma_eff_cm2, ...
    photon_mode: str | None = None
    reaction_id: str = ""
    reference: str = ""
    enabled: bool = True

    @property
    def reactant_names(self) -> tuple[str, ...]:
        return tuple(n for _, n in self.reactants)

    @property
    def product_names(self) -> tuple[str, ...]:
        return tuple(n for _, n in self.products)

    def describe(self) -> str:
        def side(terms: tuple[Term, ...]) -> str:
            return " + ".join((f"{c} {n}" if c != 1 else n) for c, n in terms)
        return f"{side(self.reactants)} -> {side(self.products)}"


@dataclass
class CompiledNetwork:
    """Compiled, layer-routed reaction set + fluid mass-action matrices."""

    species_names: tuple[str, ...]
    index: dict[str, int]
    fluid: tuple[Reaction, ...]
    kinetic: tuple[Reaction, ...]
    photon: tuple[Reaction, ...]
    stoichiometry: np.ndarray  # Sigma, shape (S, R_fluid)  net products-reactants
    orders: np.ndarray  # O, shape (S, R_fluid)  reactant multiplicities
    rate_si: np.ndarray  # (R_fluid,) base SI rate constant (K0 for te_dependent)
    is_dynamic: np.ndarray  # (R_fluid,) bool: rate depends on T_e/T_gas (runtime eval)

    @property
    def n_species(self) -> int:
        return len(self.species_names)

    @property
    def n_fluid(self) -> int:
        return self.stoichiometry.shape[1]


class ReactionNetwork:
    """Declare reactions with :meth:`add`, then :meth:`compile` for the layered engine."""

    def __init__(self, registry: SpeciesRegistry | None = None) -> None:
        self.registry: SpeciesRegistry = registry if registry is not None else excimer_system_registry()
        self.reactions: list[Reaction] = []

    def add(
        self,
        reactants: str | Iterable[str],
        products: str | Iterable[str],
        *,
        r: float | None = None,
        r_units: str = "cm3/s",
        sigma: str | None = None,
        tau_s: float | None = None,
        rate_form: RateForm | str | None = None,
        params: dict | None = None,
        photon_mode: str | None = None,
        id: str = "",
        reference: str = "",
        enabled: bool = True,
    ) -> Reaction:
        """Declare and append one reaction.

        Without ``rate_form``, ``tau_s``, ``sigma``, or ``r`` selects radiative,
        electron-impact, or constant form respectively. Explicit forms use their relevant
        argument and ``params`` metadata: for example TE-dependent rows use ``r`` as K0,
        PHOTO rows use an SI cross-section in ``r``, and stimulated rows carry their
        effective cross-section in ``params``.
        """
        rx = parse_side(reactants)
        px = parse_side(products)
        if not rx or not px:
            raise ValueError("reaction needs at least one reactant and one product")
        params = dict(params or {})

        if rate_form is None:
            if tau_s is not None:
                rate_form = RateForm.RADIATIVE
            elif sigma is not None:
                rate_form = RateForm.ELECTRON_IMPACT
            elif r is not None:
                rate_form = RateForm.CONSTANT
            else:
                raise ValueError("provide one of r=, sigma=, or tau_s=")
        rate_form = RateForm(rate_form)

        if rate_form is RateForm.RADIATIVE:
            if tau_s is None or tau_s <= 0:
                raise ValueError("radiative reaction requires a positive tau_s")
            rate_si: float | None = 1.0 / tau_s
        elif rate_form is RateForm.ELECTRON_IMPACT:
            rate_si = None  # cross-section driven coefficient lookup
        else:  # CONSTANT / TE_DEPENDENT / PHOTO / STIMULATED
            rate_si = to_si(r, r_units) if r is not None else None

        reaction = Reaction(
            reactants=rx,
            products=px,
            rate_form=rate_form,
            rate_si=rate_si,
            cross_section=sigma,
            params=params,
            photon_mode=photon_mode,
            reaction_id=id,
            reference=reference,
            enabled=enabled,
        )
        self.reactions.append(reaction)
        return reaction

    def _layer_of(self, rxn: Reaction) -> Layer:
        reactants = [self.registry.get(n) for n in rxn.reactant_names]
        if any(s.is_photon for s in reactants):
            return Layer.PHOTON
        if any(s.is_charged for s in reactants):
            return Layer.KINETIC
        return Layer.FLUID

    def compile(self) -> CompiledNetwork:
        """Validate and route active reactions, then build neutral-fluid matrices.

        Photon-reactant rows enter ``photon``; remaining charged-reactant rows enter
        ``kinetic``; all-neutral rows enter ``fluid``. The returned ``stoichiometry`` and
        ``orders`` arrays have shape ``(n_species, n_fluid_reactions)`` in registry order,
        while ``rate_si`` and ``is_dynamic`` have one entry per fluid reaction. Compilation
        does not mutate declared reactions.
        """
        names = self.registry.names
        index = {n: i for i, n in enumerate(names)}
        n_species = len(names)

        active = [r for r in self.reactions if r.enabled]
        # validate every species is known up front (clear error before matrix build)
        for rxn in active:
            for name in (*rxn.reactant_names, *rxn.product_names):
                if name not in index:
                    raise KeyError(
                        f"reaction {rxn.describe()!r} references unknown species {name!r}"
                    )

        fluid: list[Reaction] = []
        kinetic: list[Reaction] = []
        photon: list[Reaction] = []
        for rxn in active:
            layer = self._layer_of(rxn)
            (fluid if layer is Layer.FLUID else kinetic if layer is Layer.KINETIC else photon).append(rxn)

        n_fluid = len(fluid)
        orders = np.zeros((n_species, n_fluid), dtype=np.float64)
        stoich = np.zeros((n_species, n_fluid), dtype=np.float64)
        rate_si = np.zeros(n_fluid, dtype=np.float64)
        is_dynamic = np.zeros(n_fluid, dtype=bool)
        for r, rxn in enumerate(fluid):
            for count, name in rxn.reactants:
                orders[index[name], r] += count
                stoich[index[name], r] -= count
            for count, name in rxn.products:
                stoich[index[name], r] += count
            rate_si[r] = rxn.rate_si if rxn.rate_si is not None else 0.0
            is_dynamic[r] = rxn.rate_form is RateForm.TE_DEPENDENT

        self._warn_orphan_charge(active, index)

        return CompiledNetwork(
            species_names=names,
            index=index,
            fluid=tuple(fluid),
            kinetic=tuple(kinetic),
            photon=tuple(photon),
            stoichiometry=stoich,
            orders=orders,
            rate_si=rate_si,
            is_dynamic=is_dynamic,
        )

    def _warn_orphan_charge(self, active: list[Reaction], index: dict[str, int]) -> None:
        """Warn if a charged species is produced but never consumed (charge would leak).

        The lumped ``M+`` sink is exempt (it exists precisely to absorb the paper's
        untracked-ion channels).
        """
        produced: set[str] = set()
        consumed: set[str] = set()
        for rxn in active:
            produced.update(rxn.product_names)
            consumed.update(rxn.reactant_names)
        for name in produced - consumed:
            if name == "M+":
                continue
            if self.registry.get(name).is_charged and name not in _warned_orphans:
                _warned_orphans.add(name)  # once per species per process
                note = _KNOWN_ORPHANS.get(name, "add a recombination/loss channel or route to M+.")
                warnings.warn(
                    f"charged species {name!r} is produced but never consumed "
                    f"(charge may not close); {note}",
                    OrphanChargeWarning,
                    stacklevel=2,
                )


class OrphanChargeWarning(UserWarning):
    """A charged species has no loss channel in the loaded network (see _KNOWN_ORPHANS)."""


#: Documented orphans inherited from the literature source — not project omissions.
_KNOWN_ORPHANS: dict[str, str] = {
    "Xe+": (
        "known orphan inherited from Chen Table III / Ma et al. PSST 33 075012 Table 2 "
        "(neither paper gives Xe+ a loss channel); in-pulse effect negligible (~3% of "
        "late-time ions). An opt-in Xe+ + F- -> Xe + F ion-ion channel would close it "
        "if afterglow charge balance ever matters."
    ),
}
_warned_orphans: set[str] = set()


def radiative_channels(
    compiled: CompiledNetwork, registry: SpeciesRegistry
) -> list[tuple[str, float, tuple[Term, ...]]]:
    """(upper state, lifetime τ_s, non-photon lower products) per radiative channel.

    The engine's spontaneous-emission driver consumes this; ``rate_si`` for a
    RADIATIVE reaction is 1/τ, so a missing/zero rate maps to an infinite lifetime.
    """
    out: list[tuple[str, float, tuple[Term, ...]]] = []
    for rxn in compiled.fluid:
        if rxn.rate_form is RateForm.RADIATIVE and rxn.photon_mode:
            upper = rxn.reactant_names[0]
            tau = 1.0 / rxn.rate_si if rxn.rate_si else np.inf
            lower = tuple(term for term in rxn.products if not registry.get(term[1]).is_photon)
            out.append((upper, tau, lower))
    return out


def stimulated_sigma(compiled: CompiledNetwork) -> float:
    """Effective stimulated-emission cross-section σ_eff (m²), 0.0 if undeclared."""
    for rxn in compiled.photon:
        if rxn.rate_form is RateForm.STIMULATED and "sigma_eff_m2" in rxn.params:
            return float(rxn.params["sigma_eff_m2"])
    return 0.0


def stimulated_channel(
    compiled: CompiledNetwork, registry: SpeciesRegistry
) -> tuple[str, float, tuple[Term, ...], str] | None:
    """The stimulated-emission channel as ``(upper_state, σ_eff_m2, lower_products, photon_mode)``.

    ``upper_state`` is the excited-species reactant (e.g. ``KrF*``/``ArF*``); ``lower_products``
    are the non-photon products; ``photon_mode`` the emitted photon species. ``None`` if the
    network declares no stimulated-emission reaction.
    """
    for rxn in compiled.photon:
        if rxn.rate_form is RateForm.STIMULATED and "sigma_eff_m2" in rxn.params:
            upper = next((n for n in rxn.reactant_names if not registry.get(n).is_photon), None)
            lower = tuple(term for term in rxn.products if not registry.get(term[1]).is_photon)
            if upper is not None:
                return upper, float(rxn.params["sigma_eff_m2"]), lower, rxn.photon_mode or ""
    return None


def photo_absorbers(
    compiled: CompiledNetwork, registry: SpeciesRegistry
) -> list[tuple[str, float, str]]:
    """UV photo-absorption channels as ``(absorber_species, σ_m2, photon_mode)``.

    Enumerates the PHOTO reactions (e.g. Chen R40 ``Ar**+hν→Ar⁺+e``, R41 ``F⁻+hν→F+e``); the
    absorber is the non-photon reactant and σ is the bare cross-section (m²) from the reaction's
    ``sigma_m2`` param.
    """
    out: list[tuple[str, float, str]] = []
    for rxn in compiled.photon:
        if rxn.rate_form is RateForm.PHOTO and "sigma_m2" in rxn.params:
            absorber = next((n for n in rxn.reactant_names if not registry.get(n).is_photon), None)
            if absorber is not None:
                out.append((absorber, float(rxn.params["sigma_m2"]), rxn.photon_mode or ""))
    return out


def photo_channels(
    compiled: CompiledNetwork, registry: SpeciesRegistry,
) -> list[tuple[str, str, float, tuple[Term, ...], tuple[Term, ...], float, str]]:
    """Declared PHOTO stoichiometry for optics-owned species feedback.

    Rows are ``(reaction_id, absorber, sigma_m2, reactants, products,
    electron_energy_ev, photon_mode)`` — the photoelectron birth energy is 0.0 when the
    data row declares none.  ``photon_mode`` stays LAST (downstream filters select on
    ``channel[-1]``).  Keeping the declared multiplicities lets the optical hazard update
    preserve atoms and charge without maintaining a second hard-coded photo-chemistry table.
    """
    out = []
    for reaction in compiled.photon:
        if reaction.rate_form is not RateForm.PHOTO or "sigma_m2" not in reaction.params:
            continue
        absorber = next(
            (name for _count, name in reaction.reactants if not registry.get(name).is_photon),
            None,
        )
        if absorber is not None:
            out.append((
                reaction.reaction_id,
                absorber,
                float(reaction.params["sigma_m2"]),
                reaction.reactants,
                reaction.products,
                float(reaction.params.get("electron_energy_eV", 0.0)),
                reaction.photon_mode or "",
            ))
    return out
