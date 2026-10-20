"""Mass-action plasma chemistry: dn/dt = ν @ R(n), with an analytic Jacobian.

Built from a :class:`~reactions.network.CompiledNetwork`, this assembles the
per-cell charged-and-neutral right-hand side and its Jacobian from the stoichiometry
matrix ν and the reactant-order matrix O:

    R[r](n) = k[r] · Πₛ nₛ**O[s,r]          (rate of progress of reaction r)
    dn/dt   = ν @ R(n)
    J[i,j]  = Σ_r ν[i,r] · k[r] · ∂R[r]/∂n[j]

``∂R/∂n`` is built by the zero-safe *product-of-the-others* form (never divides by a
density), so species that start at zero (excited states, excimers, photons) are
handled correctly. This class is the reference/oracle path (single state vector, NumPy
+ SciPy Radau); the numba per-cell parallel step lives in :mod:`kinetics.integrator`
and reuses the padded reactant tables built here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from particles.species import SpeciesRegistry
from reactions.energy_transfer import reaction_energy_loss_ev
from reactions.network import CompiledNetwork, Reaction
from reactions.types import RateForm


@dataclass
class FluidChemistry:
    """Assemble and evaluate the fluid-layer mass-action ODE from a compiled network."""

    nu: NDArray[np.float64]  # ν stoichiometry (S, R)
    orders: NDArray[np.float64]  # O reactant orders (S, R)
    k: NDArray[np.float64]  # SI rate constants (R,)
    species_names: tuple[str, ...]
    index: dict[str, int]
    # Padded reactant tables (for the numba per-cell path): shape (R, max_reactants)
    react_species: NDArray[np.int64]
    react_order: NDArray[np.int64]
    react_count: NDArray[np.int64]  # number of reactant species per reaction (R,)
    atom_matrix: NDArray[np.float64]  # (n_elements, S) heavy-atom counts for conservation
    elements: tuple[str, ...]
    # Populated by for_plasma(): per-reaction id + a mask of Boltzmann-derived (dynamic)
    # rate constants whose k is filled at runtime from the LUT (see set_dynamic_rates).
    reaction_ids: tuple[str, ...] = ()
    dynamic_mask: NDArray[np.bool_] | None = None
    # Reactions whose rate is an explicit analytic function of Tₑ (RateForm.TE_DEPENDENT):
    # reaction index -> (K0 [SI], te_power, E0_eV) evaluated as K0·Tₑ^p·exp(−E0/Tₑ). These are
    # NOT Boltzmann-LUT quantities — their coefficients are published with the reaction (e.g. Chen
    # et al. 2025 Table III R9/R13/R14/R18/R19), so they must be evaluated from their own
    # parameters rather than looked up.
    te_analytic: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    # Per-reaction inelastic energy cost ε_th (eV), 0 where a reaction consumes no electron
    # energy. ``energy_costly_mask`` selects the electron-driven, energy-consuming reactions —
    # the ones whose rate must be capped by the electron power actually available.
    threshold_ev: NDArray[np.float64] | None = None
    energy_costly_mask: NDArray[np.bool_] | None = None
    electron_multiplying_mask: NDArray[np.bool_] | None = None
    electron_changing_mask: NDArray[np.bool_] | None = None
    # Reactions that change ANY charged species (ions included). The heavy partition's
    # electron masks are all-False by construction, but its ion-ion recombination and
    # dimerization rows still reach ν ~ 1e9 s⁻¹ during conduction — this mask lets the
    # orchestrator grade the heavy-chemistry cadence on that stiffness.
    charged_changing_mask: NDArray[np.bool_] | None = None

    @classmethod
    def from_compiled(
        cls, compiled: CompiledNetwork, registry: SpeciesRegistry
    ) -> "FluidChemistry":
        """Build contiguous mass-action arrays from ``compiled`` fluid-layer rows.

        Species preserve ``compiled.species_names`` order. ``nu`` and ``orders`` have
        shape ``(n_species, n_reactions)``; the padded reactant tables have shape
        ``(n_reactions, max_distinct_reactants)``. Arrays are newly allocated, so later
        coefficient updates do not mutate the compiled network.
        """
        nu = np.ascontiguousarray(compiled.stoichiometry, dtype=np.float64)
        orders = np.ascontiguousarray(compiled.orders, dtype=np.float64)
        k = np.ascontiguousarray(compiled.rate_si, dtype=np.float64)
        n_species, n_react = orders.shape

        # Sparse reactant structure -> padded (R, max_reactants) integer tables.
        per_reaction = [
            [(s, int(orders[s, r])) for s in range(n_species) if orders[s, r] > 0]
            for r in range(n_react)
        ]
        max_reactants = max((len(t) for t in per_reaction), default=1)
        react_species = np.full((n_react, max_reactants), -1, dtype=np.int64)
        react_order = np.zeros((n_react, max_reactants), dtype=np.int64)
        react_count = np.zeros(n_react, dtype=np.int64)
        for r, terms in enumerate(per_reaction):
            react_count[r] = len(terms)
            for a, (s, o) in enumerate(terms):
                react_species[r, a] = s
                react_order[r, a] = o

        # Heavy-atom conservation matrix (exclude the electron pseudo-element "e").
        elements = sorted(
            {
                el
                for name in compiled.species_names
                for el in registry.get(name).constituents
                if el != "e"
            }
        )
        atom_matrix = np.zeros((len(elements), n_species), dtype=np.float64)
        for j, name in enumerate(compiled.species_names):
            constituents = registry.get(name).constituents
            for e_i, el in enumerate(elements):
                atom_matrix[e_i, j] = constituents.count(el)

        return cls(
            nu=nu,
            orders=orders,
            k=k,
            species_names=compiled.species_names,
            index=compiled.index,
            react_species=react_species,
            react_order=react_order,
            react_count=react_count,
            atom_matrix=atom_matrix,
            elements=tuple(elements),
        )

    @classmethod
    def for_fluid(
        cls,
        compiled: CompiledNetwork,
        registry: SpeciesRegistry,
        fluid_species: tuple[str, ...],
        *,
        exclude_radiative: bool = True,
    ) -> "FluidChemistry":
        """Build over a neutral-species subset (the :class:`~kinetics.fluid.NeutralFluid`).

        ν/O span exactly ``fluid_species`` (same order as the fluid's columns). Reactions
        with a ``RADIATIVE`` rate form are excluded when ``exclude_radiative`` — that decay
        (KrF*→Kr+F+hν) is driven by the photon layer's spontaneous emission, so keeping it
        here too would double-count the loss. Cross-layer products (charged/photon) simply
        fall outside the subset and are dropped from the fluid update (handled elsewhere).
        """
        reactions: list[Reaction] = [
            r for r in compiled.fluid
            if not (exclude_radiative and r.rate_form is RateForm.RADIATIVE)
        ]
        index = {name: i for i, name in enumerate(fluid_species)}
        n_species, n_react = len(fluid_species), len(reactions)
        orders = np.zeros((n_species, n_react), dtype=np.float64)
        stoich = np.zeros((n_species, n_react), dtype=np.float64)
        rate_si = np.zeros(n_react, dtype=np.float64)
        is_dynamic = np.zeros(n_react, dtype=bool)
        for r, rxn in enumerate(reactions):
            for count, name in rxn.reactants:
                if name in index:  # fluid reactions have all-neutral reactants -> always in subset
                    orders[index[name], r] += count
                    stoich[index[name], r] -= count
            for count, name in rxn.products:
                if name in index:  # non-fluid products (ions/photons) are cross-layer -> dropped
                    stoich[index[name], r] += count
            rate_si[r] = rxn.rate_si if rxn.rate_si is not None else 0.0
            is_dynamic[r] = rxn.rate_form is RateForm.TE_DEPENDENT
        compiled_subset = CompiledNetwork(
            species_names=tuple(fluid_species), index=index, fluid=tuple(reactions),
            kinetic=(), photon=(), stoichiometry=stoich, orders=orders,
            rate_si=rate_si, is_dynamic=is_dynamic,
        )
        return cls.from_compiled(compiled_subset, registry)

    @classmethod
    def for_plasma(
        cls,
        compiled: CompiledNetwork,
        registry: SpeciesRegistry,
        all_species: tuple[str, ...],
        *,
        exclude_radiative: bool = True,
        partition: str = "all",
    ) -> "FluidChemistry":
        """Build the combined charged+neutral chemistry (fluid + electron-impact reactions).

        ν/O span **all** species. Both ``compiled.fluid`` (all-neutral) and
        ``compiled.kinetic`` (≥1 charged reactant — electron-impact ionization/excitation/
        attachment and charged recombination) are included; ``RADIATIVE`` decay is excluded
        (driven by the photon layer's spontaneous emission — keeping it here double-counts).
        Photon-layer reactions are excluded (the scalar cavity optics handles them).

        Electron-impact / Te-dependent reactions get ``k = 0`` here — their runtime rate
        constant is filled from the Boltzmann LUT by :meth:`set_dynamic_rates` (keyed on the
        local reduced field / mean energy). ``reaction_ids``/``dynamic_mask`` record which.
        """
        if partition not in {"all", "electron_source", "heavy_source"}:
            raise ValueError(
                "plasma chemistry partition must be 'all', 'electron_source', or 'heavy_source'"
            )
        reactions_all: list[Reaction] = [
            r for r in (*compiled.fluid, *compiled.kinetic)
            if not (exclude_radiative and r.rate_form is RateForm.RADIATIVE)
        ]
        if partition == "all":
            reactions = reactions_all
        else:
            def is_electron_source(reaction: Reaction) -> bool:
                electron_net = (
                    sum(count for count, name in reaction.products if name == "e-")
                    - sum(count for count, name in reaction.reactants if name == "e-")
                )
                return (
                    reaction.rate_form in (RateForm.ELECTRON_IMPACT, RateForm.TE_DEPENDENT)
                    or electron_net != 0
                    or reaction_energy_loss_ev(reaction) > 0.0
                )

            want_electron = partition == "electron_source"
            reactions = [r for r in reactions_all if is_electron_source(r) == want_electron]
        index = {name: i for i, name in enumerate(all_species)}
        n_species, n_react = len(all_species), len(reactions)
        orders = np.zeros((n_species, n_react), dtype=np.float64)
        stoich = np.zeros((n_species, n_react), dtype=np.float64)
        rate_si = np.zeros(n_react, dtype=np.float64)
        dynamic_mask = np.zeros(n_react, dtype=bool)
        reaction_ids: list[str] = []
        te_analytic: dict[int, tuple[float, float, float]] = {}
        thresholds = np.zeros(n_react, dtype=np.float64)
        for r, rxn in enumerate(reactions):
            # Shared with transport's diagnostic closure through reactions.energy_transfer.
            # E0_eV is an analytic-rate activation, not an implicit energy-transfer value.
            thresholds[r] = reaction_energy_loss_ev(rxn)
            for count, name in rxn.reactants:
                orders[index[name], r] += count
                stoich[index[name], r] -= count
            for count, name in rxn.products:
                stoich[index[name], r] += count
            dyn = rxn.rate_form in (RateForm.ELECTRON_IMPACT, RateForm.TE_DEPENDENT)
            dynamic_mask[r] = dyn
            if rxn.rate_form is RateForm.TE_DEPENDENT:
                # Capture (K0, p, E0) before rate_si is zeroed below: this reaction's rate is
                # published analytically, not derived from a Boltzmann solve.
                te_analytic[r] = (
                    float(rxn.rate_si or 0.0),
                    float(rxn.params.get("te_power", 0.0)),
                    float(rxn.params.get("E0_eV", 0.0)),
                )
            rate_si[r] = 0.0 if dyn else (rxn.rate_si if rxn.rate_si is not None else 0.0)
            reaction_ids.append(rxn.reaction_id)
        compiled_subset = CompiledNetwork(
            species_names=tuple(all_species), index=index, fluid=tuple(reactions),
            kinetic=(), photon=(), stoichiometry=stoich, orders=orders,
            rate_si=rate_si, is_dynamic=dynamic_mask,
        )
        chem = cls.from_compiled(compiled_subset, registry)
        chem.reaction_ids = tuple(reaction_ids)
        chem.dynamic_mask = dynamic_mask
        chem.te_analytic = te_analytic
        chem.threshold_ev = thresholds
        chem.energy_costly_mask = thresholds > 0.0
        # Electron-multiplying reactions: net electron production > 0. These are the avalanche
        # drivers and set the stiffness — note this INCLUDES the analytic TE_DEPENDENT stepwise
        # ionizations (Chen R9/R18: e + Ar**/Ne2* -> ion + 2e).
        e_row = index.get("e-")
        chem.electron_multiplying_mask = (
            np.asarray(stoich[e_row, :] > 0.0) if e_row is not None
            else np.zeros(n_react, dtype=bool)
        )
        # Electron-changing reactions include both production and consumption so the cadence
        # resolves ionization and attachment rather than grading only avalanche growth.
        chem.electron_changing_mask = (
            np.asarray(stoich[e_row, :] != 0.0) if e_row is not None
            else np.zeros(n_react, dtype=bool)
        )
        charged_rows = [
            index[name] for name in all_species
            if name in registry and registry.get(name).is_charged
        ]
        chem.charged_changing_mask = (
            np.any(stoich[charged_rows, :] != 0.0, axis=0)
            if charged_rows else np.zeros(n_react, dtype=bool)
        )
        return chem

    def progress_frequency(self, densities: NDArray[np.float64], k_percell, columns) -> float:
        """Fastest per-cell rate of progress / reactant density over ``columns`` (s⁻¹).

        Used to size the chemistry sub-cycling: for ``e + X → …`` the electron growth frequency is
        ``k·n_X``, so we form ``k·Π n^O`` with the *actual* reactant densities and divide out the
        scarcest reactant, which bounds the fastest fractional reactant change for that row.
        """
        n = np.asarray(densities, dtype=np.float64)
        k = np.asarray(k_percell, dtype=np.float64)
        worst = 0.0
        for r in np.flatnonzero(columns):
            cnt = int(self.react_count[r])
            dens = [n[:, self.react_species[r, a]] for a in range(cnt)]
            if not dens:
                continue
            prod = np.ones(n.shape[0], dtype=np.float64)
            orders = [self.react_order[r, a] for a in range(cnt)]
            for d, o in zip(dens, orders):
                prod = prod * (d ** o if o != 1 else d)
            # Divide out the SCARCEST reactant: for e + X -> X+ + 2e the electrons are the
            # multiplying species, so their growth frequency is R/nₑ = k·n_X — dividing by the
            # abundant target instead would understate the stiffness by many orders.
            scarcest = np.maximum(np.min(np.stack(dens), axis=0), 1.0)
            worst = max(worst, float(np.max(k[:, r] * prod / scarcest)))
        return worst

    def energy_consumption_rate_percell(
        self, densities: NDArray[np.float64], k_percell: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Requested threshold-energy consumption rate in each cell (eV m⁻³ s⁻¹)."""
        n = np.asarray(densities, dtype=np.float64)
        k = np.asarray(k_percell, dtype=np.float64)
        if n.ndim != 2 or n.shape[1] != self.n_species:
            raise ValueError("densities must have shape (n_cells, n_species)")
        if k.shape != (n.shape[0], self.n_reactions):
            raise ValueError("k_percell must have shape (n_cells, n_reactions)")
        eth = self.threshold_ev
        if eth is None:
            return np.zeros(n.shape[0], dtype=np.float64)
        requested = np.zeros(n.shape[0], dtype=np.float64)
        for r in np.flatnonzero(eth > 0.0):
            progress = k[:, r].copy()
            for a in range(int(self.react_count[r])):
                col = self.react_species[r, a]
                order = self.react_order[r, a]
                progress *= n[:, col] if order == 1 else n[:, col] ** order
            requested += eth[r] * progress
        return requested

    @staticmethod
    def _arrhenius_te(params: tuple[float, float, float], te):
        """``k = K0·Tₑ^p·exp(−E0/Tₑ)`` (SI), safe at Tₑ → 0. Scalar or array ``te``."""
        k0, power, e0 = params
        t = np.asarray(te, dtype=np.float64)
        safe = np.where(t > 0.0, t, np.inf)  # Tₑ<=0 -> exp(-E0/inf)=1 * 0**p -> 0 below
        out = k0 * safe**power * np.exp(-e0 / safe)
        return np.where(t > 0.0, out, 0.0)

    def set_dynamic_rates(self, coeffs, key) -> None:
        """Fill the Boltzmann-derived rate constants ``k`` from the LUT at ``key``.

        ``coeffs`` is a ``boltzmann.CoefficientTable``; ``key`` is a scalar lookup value
        (uniform |E|/N or mean energy) — the stub uses a domain-representative scalar, so a
        single ``k`` per reaction is exact. TE_DEPENDENT reactions are evaluated from their own
        published ``K0·Tₑ^p·exp(−E0/Tₑ)`` (Tₑ read from the table, so this works for either key
        kind); only ELECTRON_IMPACT reactions come from the LUT.
        """
        if self.dynamic_mask is None:
            return
        te = float(np.asarray(coeffs.te(key))) if self.te_analytic else 0.0
        for r, is_dyn in enumerate(self.dynamic_mask):
            if not is_dyn:
                continue
            if r in self.te_analytic:
                self.k[r] = float(self._arrhenius_te(self.te_analytic[r], te))
            else:
                self.k[r] = float(coeffs.rate(self.reaction_ids[r], key))

    def dynamic_rates_percell(self, coeffs, keys_percell) -> NDArray[np.float64]:
        """Per-cell rate-coefficient matrix ``k[cell, reaction]`` for the LMEA path.

        Static (non-dynamic) reactions broadcast their constant ``rate_si``; ELECTRON_IMPACT
        reactions are interpolated from the LUT at the per-cell lookup key ``keys_percell``
        (shape ``(n_cells,)`` — mean energy ⟨ε⟩ or E/N per cell); TE_DEPENDENT reactions are
        evaluated per cell from their own published ``K0·Tₑ^p·exp(−E0/Tₑ)``. Does not mutate
        ``self.k``; returns a fresh ``(n_cells, n_reactions)`` array for
        :meth:`FluidIntegrator.advance(..., k_percell=...)`.
        """
        keys = np.asarray(keys_percell, dtype=np.float64)
        k = np.tile(self.k, (keys.shape[0], 1))  # (n_cells, R); static rates broadcast
        if self.dynamic_mask is None:
            return k
        te = np.asarray(coeffs.te(keys), dtype=np.float64) if self.te_analytic else None
        for r, is_dyn in enumerate(self.dynamic_mask):
            if not is_dyn:
                continue
            if r in self.te_analytic:
                k[:, r] = self._arrhenius_te(self.te_analytic[r], te)
            else:
                k[:, r] = np.asarray(coeffs.rate(self.reaction_ids[r], keys), dtype=np.float64)
        return k

    @property
    def n_species(self) -> int:
        return self.nu.shape[0]

    @property
    def n_reactions(self) -> int:
        return self.nu.shape[1]

    def rates_of_progress(self, n: NDArray[np.float64]) -> NDArray[np.float64]:
        """R(n): rate of progress of each reaction, ``k[r]·Πₛ nₛ**O[s,r]``."""
        prod = np.ones(self.n_reactions, dtype=np.float64)
        for r in range(self.n_reactions):
            for a in range(int(self.react_count[r])):
                s = self.react_species[r, a]
                prod[r] *= n[s] ** self.react_order[r, a]
        return self.k * prod

    def rhs(self, n: NDArray[np.float64]) -> NDArray[np.float64]:
        """dn/dt = ν @ R(n)."""
        return self.nu @ self.rates_of_progress(n)

    def jacobian(self, n: NDArray[np.float64]) -> NDArray[np.float64]:
        """Analytic Jacobian ∂(dn/dt)/∂n via the zero-safe product-of-others form."""
        n_react = self.n_reactions
        deriv = np.zeros((self.n_species, n_react), dtype=np.float64)  # ∂R̃[r]/∂n[j]
        for r in range(n_react):
            cnt = int(self.react_count[r])
            for a in range(cnt):
                j = self.react_species[r, a]
                oj = self.react_order[r, a]
                others = 1.0
                for b in range(cnt):
                    if b != a:
                        s = self.react_species[r, b]
                        others *= n[s] ** self.react_order[r, b]
                deriv[j, r] = oj * (n[j] ** (oj - 1)) * others
        weighted = self.nu * self.k[np.newaxis, :]  # (S, R)
        return weighted @ deriv.T  # (S, S)

    def atom_counts(self, n: NDArray[np.float64]) -> NDArray[np.float64]:
        """Total heavy-atom count per element (conserved by balanced neutral chemistry)."""
        return self.atom_matrix @ n

    def integrate_reference(
        self,
        n0: NDArray[np.float64],
        t_span: tuple[float, float],
        *,
        method: str = "Radau",
        rtol: float = 1e-8,
        atol: float = 1.0,
        t_eval: NDArray[np.float64] | None = None,
    ):
        """Single-cell reference integration via SciPy (validation oracle for the numba step)."""
        from scipy.integrate import solve_ivp

        return solve_ivp(
            lambda _t, y: self.rhs(y),
            t_span,
            np.asarray(n0, dtype=np.float64),
            method=method,
            jac=lambda _t, y: self.jacobian(y),
            rtol=rtol,
            atol=atol,
            t_eval=t_eval,
        )
