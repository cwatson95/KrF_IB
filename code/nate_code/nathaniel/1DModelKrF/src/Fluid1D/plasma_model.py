"""Charged-species plasma model interface, production drift-diffusion model, and test stub.

:class:`PlasmaModel` is implemented by the 1D drift-diffusion charged-fluid model: given the
per-cell electric field and the neutral/charged fluids, it
advances the charged densities + electron energy one timestep and reports what it did.

:class:`StubPlasmaModel` is a minimal, bounded stand-in for fixed-voltage tests:
it pins the electron mean energy to Tₑ from the Boltzmann LUT (or a fixed value) and reports
a Joule-heating deposited power for the efficiency ledger. It does **not** transport or grow
the charged densities; the orchestrator holds those quasi-steady while chemistry exercises the
remaining coupling path.
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from kinetics.fluid import NeutralFluid
from particles.constants import E_CHARGE, EPS0, K_B, M_E
from reactions.network import CompiledNetwork
from reactions.energy_transfer import reaction_energy_loss_ev
from reactions.types import RateForm

from .boltzmann import CoefficientTable, reduced_field_td
from .charged_fluid import ChargedFluid
from .transport import face_average, flux_divergence, sg_interior_flux

# Numerical floor on the electron density used to recover the mean energy ⟨ε⟩ = n_ε/(nₑ·e).
# 1 electron per m³ is utterly negligible physically (seeds are ≥1e13 m⁻³) but keeps nₑ·e far
# above the denormal range, so an attachment-emptied cell cannot produce 0/0.
_N_E_FLOOR_M3 = 1.0

# Default reduced ion mobility μ·N (1/(V·m·s)) for ions absent from the mobility table
# (≈1 cm²/V·s at the Loschmidt density). Only touches inactive, zero-density ions in practice.
_DEFAULT_ION_MU_N = 1.0e-4 * 2.6867805e25


@dataclass
class EnergyClosureReport:
    """Per-cell, interval-integrated electron-energy terms from one transport advance.

    Arrays use J/m^3.  Positive field work adds electron energy; elastic and continuous
    inelastic terms are positive sinks; ``coefficient_bound_adjustment`` is signed.
    Chemistry fills accepted cost and reconciliation at the orchestration layer.
    """

    field_work_j_m3: np.ndarray
    elastic_transfer_j_m3: np.ndarray
    continuous_inelastic_sink_j_m3: np.ndarray
    accepted_reaction_cost_j_m3: np.ndarray
    signed_reconciliation_j_m3: np.ndarray
    coefficient_bound_adjustment_j_m3: np.ndarray
    # source_coupled: energy actually delivered into the chemistry stage (signed; the
    # orchestration layer fills it — transport advances leave it None).
    stage_delivered_j_m3: np.ndarray | None = None
    maximum_energy_jump_ev: float = 0.0
    limiter_activity: int = 0
    max_ampere_residual: float = 0.0
    max_energy_residual_j_m3: float = 0.0


@dataclass
class StepReport:
    """What a charged-fluid advance did this step: energy deposition + bookkeeping."""

    deposited_power_w: float  # power coupled into electrons this step (feeds LaserEnergetics)
    mean_electron_energy_ev: float = 0.0
    dt: float = 0.0
    conduction_current_a: float = 0.0  # A·⟨J_cond⟩_x drawn at the cathode node (feeds the circuit)
    # Differential gap conductance G = ∂I_cond/∂V_b (S). Lets the circuit linearise the plasma
    # load about the current state (semi-implicit gap coupling) instead of freezing the current.
    gap_conductance_s: float = 0.0
    gap_energy_absorbed_j: float | None = None
    energy_closure: EnergyClosureReport | None = None


@dataclass(frozen=True)
class _TransportCoefficients:
    """Per-cell transport data evaluated at one substep's recovered mean energy."""

    mu_e: np.ndarray
    d_e: np.ndarray
    mu_eps: np.ndarray
    d_eps: np.ndarray
    mu_all: tuple[np.ndarray, ...]
    d_face_all: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class _SpeciesTransportResult:
    """Flux information retained after one conservative species update."""

    gamma_e_cell: np.ndarray
    electron_flux_left: float
    electron_flux_right: float
    j_cond_cell: np.ndarray
    max_ampere_residual: float


@dataclass(frozen=True)
class _EnergyTransportResult:
    """Electron-energy state and ledger increments from one transport substep."""

    n_eps: np.ndarray
    field_work_j_m3: np.ndarray
    elastic_transfer_j_m3: np.ndarray
    continuous_sink_j_m3: np.ndarray
    bound_adjustment_j_m3: np.ndarray
    deposited_j: float
    maximum_energy_jump_ev: float
    limiter_activity: int
    max_energy_residual_j_m3: float


class PlasmaModel(ABC):
    """Interface for the 1D charged-species fluid (transport + energy + reactions)."""

    @abstractmethod
    def advance(
        self,
        dt: float,
        *,
        efield: np.ndarray,
        neutral_fluid: NeutralFluid,
        charged_fluid: ChargedFluid,
        coeffs: CoefficientTable,
        field=None,
    ) -> StepReport:
        """Advance the charged fluid one ``dt`` in place; return a :class:`StepReport`.

        ``efield`` is the per-cell vec-3 field ``(N, 3)``; ``coeffs`` supplies transport and
        rate coefficients keyed on |E|/N (or mean energy). ``field`` is an optional per-step
        field handle (:class:`Fluid1D.field.SemiImplicitField`-like, duck-typed:
        ``solve(rho, sigma_cells, j_expl_cells, h) -> axial E``): when given, the model
        re-solves the semi-implicit field every internal substep and ``efield`` serves only
        as the bootstrap/sizing field; when ``None`` the given field stays frozen over the
        whole ``dt`` (the frozen-field staggered-explicit coupling, unstable past breakdown).
        """


class StubPlasmaModel(PlasmaModel):
    """Fixed-Tₑ, quasi-steady stand-in: pins electron energy, reports Joule-heating power.

    ``fixed_te_ev`` pins Tₑ directly; when ``None`` (default) Tₑ is read from the Boltzmann
    LUT at the per-cell |E|/N, so the circuit voltage → E → E/N → Tₑ chain (and the F2 sweep)
    actually drives the electron energy. Does not mutate charged densities.
    """

    def __init__(self, *, fixed_te_ev: float | None = None) -> None:
        self.fixed_te_ev = fixed_te_ev

    def advance(self, dt, *, efield, neutral_fluid, charged_fluid, coeffs, field=None) -> StepReport:
        n_gas = neutral_fluid.density.sum(axis=1)  # (N,) total neutral density per cell
        emag = np.linalg.norm(np.asarray(efield, dtype=np.float64), axis=1)  # (N,)
        en_td = reduced_field_td(emag, n_gas)  # (N,) E/N in townsend
        if coeffs.key_kind == "mean_energy_eV":
            equilibrium = coeffs.equilibrium_mean_energy(en_td)
            if equilibrium is None:
                key = charged_fluid.electron_mean_energy_ev
            else:
                key = np.asarray(equilibrium, dtype=np.float64)
        else:
            key = en_td

        if self.fixed_te_ev is not None:
            te = np.full(charged_fluid.n_cells, float(self.fixed_te_ev))
            if coeffs.key_kind == "mean_energy_eV":
                key = 1.5 * te
        else:
            te = np.asarray(coeffs.te(key), dtype=np.float64)
        # Pin the electron mean energy field (<eps> = 3/2 Te).
        charged_fluid.electron_mean_energy_ev[:] = 1.5 * te

        # Joule heating P = e·nₑ·μ·|E|² per cell (μ = μ·N / N), summed over the domain.
        mobility_n = np.asarray(coeffs.mobility(key), dtype=np.float64)
        mu = np.where(n_gas > 0.0, mobility_n / n_gas, 0.0)
        n_e = charged_fluid.electron_density()
        p_density = E_CHARGE * n_e * mu * emag**2  # W/m³
        deposited = float(np.sum(p_density)) * charged_fluid.cell_volume

        return StepReport(
            deposited_power_w=deposited,
            mean_electron_energy_ev=float(charged_fluid.electron_mean_energy_ev.mean()),
            dt=dt,
        )


class DriftDiffusionPlasmaModel(PlasmaModel):
    """Real 1D drift-diffusion-reaction charged-fluid transport (Chen et al. 2025, §II.B).

    Advances the charged-species densities by the continuity equation ``∂n/∂t + ∂ₓΓ = S`` with
    the drift-diffusion flux ``Γ = −D∂ₓn + zμnE`` (Eqs. 1–2), and the electron mean energy by the
    energy-conservation PDE ``∂n_ε/∂t + ∂ₓΓ_ε + eE·Γ_e = S_ε`` (Eqs. 3–4). Fluxes use the
    Scharfetter–Gummel scheme (:mod:`Fluid1D.transport`). Reactions ``S`` are handled by the
    orchestrator's combined chemistry integrator; this model does **transport + energy** and
    returns the plasma **conduction current** ``A·⟨J_cond⟩`` for the discharge circuit.

    Coefficient closure is **LMEA** (Hagelaar & Pitchford 2005 §1): electron μ_e, D_e and the
    energy-transport pair μ_ε, D_ε are looked up on the local mean energy ⟨ε⟩ from ``coeffs``
    (which must be mean-energy-keyed — auto-rekeyed if not). Ion μ_i, D_i come from
    ``ion_mobility`` (constant reduced mobility + Einstein diffusion). Electrode boundary fluxes
    are drift-limited ion outflux + cathode secondary emission (γ), building the cathode fall.

    Energy-closure contract: under ``electron_energy_closure="source_coupled"`` (production)
    ``advance`` does NOT apply field work or elastic exchange to n_ε — it accumulates them
    into the per-interval ``field_work``/``elastic_transfer`` ledgers for the orchestrator to
    deliver inside the co-integrated chemistry stage, and n_ε evolves here by energy-flux
    divergence and bounds only. The ``continuous_reconciled`` diagnostic mode retains the full
    in-transport local update (it has no stage source).

    Internally sub-stepped. Two field-coupling modes (selected per ``advance`` call):

    - ``field=None`` — the caller's ``efield`` is **frozen** for the whole ``dt``
      (staggered-explicit coupling; substeps sized by the drift/diffusion CFL *and* a
      fraction of the dielectric-relaxation time ε₀/σ). Linearly unstable once
      σ > ~2ε₀/dt (nₑ ≳ 3×10¹⁹ m⁻³ at dt = 50 ps) — kept for stub-path tests and as the
      frozen-field baseline.
    - ``field=SemiImplicitField`` — the model **re-solves the semi-implicit σ-corrected
      Poisson every substep** (ohmic response implicit, diffusion explicit), which removes the
      ε₀/σ limit entirely; substeps are sized by the drift/diffusion CFL alone, so the
      fully-collapsed phase (nₑ~10²¹ m⁻³) integrates in a handful of substeps instead of
      hitting the 20 000-substep dielectric cap. The mean energy's local-field closure
      ⟨ε⟩_eq(E/N) is refreshed from each fresh field (the recovery uses it with a one-substep
      lag). This is the production coupling with the MPC circuit.
    """

    def __init__(
        self,
        *,
        reactions: CompiledNetwork | None,
        ion_mobility,
        mixture: dict[str, float],
        gas_temp_k: float = 300.0,
        area_m2: float = 1.0,
        gap_m: float = 1.0,
        secondary_emission_gamma: float = 0.05,
        energy_floor_ev: float = 0.05,
        energy_ceiling_ev: float = 60.0,
        lmea_ref_density_m3: float = 1.0e10,
        cfl_safety: float = 0.3,
        max_substeps: int = 20000,
        default_ion_mu_N: float = _DEFAULT_ION_MU_N,
        use_numba_kernel: bool = True,
        field_heating_cap_factor: float | None = 2.0,
        electron_energy_closure: str = "source_coupled",
    ) -> None:
        self.reactions = reactions
        self.ion_mobility = ion_mobility
        self.mixture = dict(mixture)
        self.gas_temp_k = float(gas_temp_k)
        self.area_m2 = float(area_m2)
        self.gap_m = float(gap_m)
        self.gamma = float(secondary_emission_gamma)
        self.eps_floor = float(energy_floor_ev)
        self.eps_ceiling = float(energy_ceiling_ev)
        # Electron density below which ⟨ε⟩ is set by the local field rather than by n_ε/nₑ
        # (far below any physical seed, far above the values where the ratio degenerates).
        self.lmea_ref_density_m3 = float(lmea_ref_density_m3)
        self.cfl_safety = float(cfl_safety)
        self.max_substeps = int(max_substeps)
        self.default_ion_mu_N = float(default_ion_mu_N)
        # Coupled-path substep loop runs in the numba kernel (transport_kernel) by default;
        # False keeps the numpy reference path (the kernel's parity baseline).
        self.use_numba_kernel = bool(use_numba_kernel)
        if electron_energy_closure not in {"source_coupled", "continuous_reconciled"}:
            raise ValueError(
                "electron_energy_closure must be 'source_coupled' or 'continuous_reconciled'"
            )
        self.electron_energy_closure = electron_energy_closure
        if field_heating_cap_factor is not None and (
            not np.isfinite(field_heating_cap_factor) or field_heating_cap_factor <= 0.0
        ):
            raise ValueError("field_heating_cap_factor must be finite and positive or None")
        self.field_heating_cap_factor = (
            np.inf if field_heating_cap_factor is None else float(field_heating_cap_factor)
        )
        self._kernel_packed = None
        self._resolved = False
        self._me_coeffs: CoefficientTable | None = None
        self._last_substeps = 0
        self._inelastic_packed = None

    # ------------------------------------------------------------------ setup (lazy)
    def _resolve(self, neutral_fluid: NeutralFluid, charged_fluid: ChargedFluid, coeffs) -> None:
        cf = charged_fluid
        self._dx = cf.dx
        self._e_col = cf.electron_col
        self._z = np.array([s.charge_number for s in cf.species], dtype=np.float64)
        self._q = self._z * E_CHARGE
        self._ion_cols = [i for i in range(cf.n_species) if i != self._e_col]
        self._pos_ion_cols = [i for i in self._ion_cols if self._z[i] > 0.0]
        self._ion_masses = np.array([s.mass for s in cf.species], dtype=np.float64)
        # Per-ion reduced mobility invariant μ·N (mixture is fixed for the run).
        self._ion_mu_N = np.zeros(cf.n_species, dtype=np.float64)
        missing = []
        for i in self._ion_cols:
            name = cf.names[i]
            if self.ion_mobility is not None and self.ion_mobility.has_ion(name):
                self._ion_mu_N[i] = self.ion_mobility.reduced_mobility(name, self.mixture)
            else:
                self._ion_mu_N[i] = self.default_ion_mu_N
                missing.append(name)
        if missing:
            warnings.warn(
                f"no ion-mobility data for {sorted(set(missing))}; using default reduced "
                f"mobility (these ions are typically inactive/zero-density)", stacklevel=2,
            )
        # Neutral masses (for the elastic energy-loss term's mixture-mean mass).
        self._neutral_masses = np.array([s.mass for s in neutral_fluid.species], dtype=np.float64)
        # Coefficients must be mean-energy-keyed for the LMEA lookup.
        self._me_coeffs = coeffs if coeffs.key_kind == "mean_energy_eV" else coeffs.to_mean_energy_key()
        # Never let ⟨ε⟩ leave the range the table resolves: outside it every accessor silently
        # clamps, so an out-of-range ⟨ε⟩ is a wrong coefficient dressed up as a valid one.
        lo_tbl, hi_tbl = self._me_coeffs.mean_energy_bounds_ev
        self._eps_bounds = (max(self.eps_floor, lo_tbl), min(self.eps_ceiling, hi_tbl))
        self._last_n_eps = cf.density[:, self._e_col] * cf.electron_mean_energy_ev * E_CHARGE
        self._inelastic_packed = self._pack_inelastic_closure(
            neutral_fluid, charged_fluid, self._me_coeffs,
        )
        self._resolved = True

    def _pack_inelastic_closure(self, neutral_fluid, charged_fluid, coeffs) -> dict[str, np.ndarray]:
        """Pack the canonical loss table and compiled reactants for the diagnostic PDE sink."""
        charged_index = charged_fluid.column
        neutral_index = neutral_fluid.column
        selected = []
        if self.reactions is not None:
            for reaction in (*self.reactions.fluid, *self.reactions.kinetic):
                loss = reaction_energy_loss_ev(reaction)
                if loss > 0.0:
                    selected.append((reaction, loss))
        count = len(selected)
        max_reactants = max((len(reaction.reactants) for reaction, _ in selected), default=1)
        react_kind = np.full((count, max_reactants), -1, dtype=np.int64)
        react_col = np.full((count, max_reactants), -1, dtype=np.int64)
        react_order = np.zeros((count, max_reactants), dtype=np.int64)
        react_count = np.zeros(count, dtype=np.int64)
        loss_ev = np.zeros(count)
        base_k = np.zeros(count)
        rate_kind = np.zeros(count, dtype=np.int64)
        analytic = np.zeros((count, 3))
        rate_lut = np.zeros((count, coeffs.key_grid.size))
        for r, (reaction, loss) in enumerate(selected):
            loss_ev[r] = loss
            react_count[r] = len(reaction.reactants)
            for a, (order, name) in enumerate(reaction.reactants):
                if name in charged_index:
                    react_kind[r, a] = 0
                    react_col[r, a] = charged_index[name]
                elif name in neutral_index:
                    react_kind[r, a] = 1
                    react_col[r, a] = neutral_index[name]
                else:
                    raise KeyError(f"energy-transfer reaction {reaction.reaction_id} has untracked reactant {name}")
                react_order[r, a] = order
            if reaction.rate_form is RateForm.ELECTRON_IMPACT:
                rate_kind[r] = 1
                rate_lut[r] = np.asarray(
                    coeffs.rate_coeff.get(reaction.reaction_id, np.zeros_like(coeffs.key_grid)),
                    dtype=np.float64,
                )
            elif reaction.rate_form is RateForm.TE_DEPENDENT:
                rate_kind[r] = 2
                analytic[r] = (
                    float(reaction.rate_si or 0.0),
                    float(reaction.params.get("te_power", 0.0)),
                    float(reaction.params.get("E0_eV", 0.0)),
                )
            else:
                base_k[r] = float(reaction.rate_si or 0.0)
        return {
            "loss_ev": np.ascontiguousarray(loss_ev),
            "base_k": np.ascontiguousarray(base_k),
            "rate_kind": np.ascontiguousarray(rate_kind),
            "analytic": np.ascontiguousarray(analytic),
            "rate_lut": np.ascontiguousarray(rate_lut),
            "react_kind": np.ascontiguousarray(react_kind),
            "react_col": np.ascontiguousarray(react_col),
            "react_order": np.ascontiguousarray(react_order),
            "react_count": np.ascontiguousarray(react_count),
        }

    def _continuous_inelastic_power(self, eps, neutral_density, charged_density) -> np.ndarray:
        """Requested continuous inelastic sink in eV m^-3 s^-1 from canonical reactions."""
        packed = self._inelastic_packed
        out = np.zeros(charged_density.shape[0])
        if packed is None or packed["loss_ev"].size == 0:
            return out
        key = np.asarray(eps, dtype=np.float64)
        te = np.asarray(self._me_coeffs.te(key), dtype=np.float64)
        for r, loss in enumerate(packed["loss_ev"]):
            kind = packed["rate_kind"][r]
            if kind == 1:
                k = np.interp(key, self._me_coeffs.key_grid, packed["rate_lut"][r])
            elif kind == 2:
                k0, power, e0 = packed["analytic"][r]
                with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                    k = np.where(te > 0.0, k0 * te**power * np.exp(-e0 / te), 0.0)
            else:
                k = np.full(key.shape, packed["base_k"][r])
            progress = k
            for a in range(int(packed["react_count"][r])):
                col = int(packed["react_col"][r, a])
                density = (charged_density[:, col]
                           if packed["react_kind"][r, a] == 0 else neutral_density[:, col])
                progress = progress * density ** int(packed["react_order"][r, a])
            out += loss * progress
        return np.nan_to_num(out, nan=0.0, posinf=np.finfo(np.float64).max)

    # ------------------------------------------------------------------ advance
    def advance(self, dt, *, efield, neutral_fluid, charged_fluid, coeffs, field=None) -> StepReport:
        """Resolve static data once, then route to the coupled kernel or NumPy driver."""
        if not self._resolved:
            self._resolve(neutral_fluid, charged_fluid, coeffs)
        if (field is not None and self.use_numba_kernel
                and getattr(field, "prescribed", None) is None):
            return self._advance_via_kernel(dt, efield, neutral_fluid, charged_fluid, field)
        return self._advance_numpy(dt, efield, neutral_fluid, charged_fluid, field)

    def _advance_numpy(
        self,
        dt: float,
        efield: np.ndarray,
        neutral_fluid: NeutralFluid,
        charged_fluid: ChargedFluid,
        field,
    ) -> StepReport:
        """Reference transport driver for frozen and semi-implicit field coupling."""
        me = self._me_coeffs
        cf = charged_fluid
        dx = cf.dx
        e_axial = np.asarray(efield, dtype=np.float64)[:, 0]  # (N,) axial field V/m
        n_gas = neutral_fluid.density.sum(axis=1)
        n_gas_safe = np.where(n_gas > 0.0, n_gas, 1.0)
        # Mixture-mean neutral mass per cell (for elastic energy loss).
        m_gas = np.where(
            n_gas > 0.0, (neutral_fluid.density @ self._neutral_masses) / n_gas_safe, self._neutral_masses.mean()
        )
        eps_gas_ev = 1.5 * K_B * self.gas_temp_k / E_CHARGE  # gas thermal energy (eV)

        # Local-field-equilibrium mean energy ⟨ε⟩_eq(E/N) per cell — the closure the LMEA
        # relation ⟨ε⟩ = n_ε/nₑ degenerates to where the electron density vanishes. In the
        # coupled mode it is refreshed from each substep's fresh field.
        eps_eq = self._equilibrium_eps(me, e_axial, n_gas_safe)

        density = cf.density  # (N, S) — mutated in place
        e_col = self._e_col
        # Electron energy density n_ε [J/m³] evolved internally; recovered to ⟨ε⟩ at the end.
        n_e0 = density[:, e_col]
        eps0 = self._recover_eps(cf.electron_mean_energy_ev * n_e0 * E_CHARGE, n_e0, eps_eq)
        n_eps = n_e0 * eps0 * E_CHARGE  # J/m³

        coupled = field is not None
        n_sub = None
        if coupled:
            # Adaptive driver: h from the drift/diffusion CFL alone — the semi-implicit field
            # solve removes the dielectric-relaxation (ε₀/σ) limit entirely.
            remaining = float(dt)
            h = self.cfl_safety * self._cfl_dt_target(e_axial, n_gas_safe, me, eps0)
            if not np.isfinite(h) or h <= 0.0:
                h = remaining
            h = min(h, remaining)
        else:
            n_sub = self._plan_substeps(dt, density[:, e_col], e_axial, n_gas_safe, me, eps0)
            h = dt / n_sub

        deposited_j = 0.0
        jcond_mean_accum = 0.0
        sigma_mean_accum = 0.0
        gap_e_abs = 0.0
        gap_v2 = 0.0
        steps_done = 0
        capped_warned = False
        field_work_j_m3 = np.zeros(cf.n_cells)
        elastic_transfer_j_m3 = np.zeros(cf.n_cells)
        continuous_sink_j_m3 = np.zeros(cf.n_cells)
        bound_adjustment_j_m3 = np.zeros(cf.n_cells)
        maximum_energy_jump_ev = 0.0
        limiter_activity = 0
        max_ampere_residual = 0.0
        max_energy_residual_j_m3 = 0.0
        while True:
            n_e = density[:, e_col]
            n_e_safe = np.maximum(n_e, _N_E_FLOOR_M3)
            eps = self._recover_eps(n_eps, n_e, eps_eq)
            # Relax the energy density onto the recovered ⟨ε⟩ so a near-empty cell's n_ε cannot
            # drift away from the local-field closure and re-poison the next substep.
            n_eps = n_e_safe * eps * E_CHARGE

            coeff = self._evaluate_transport_coefficients(eps, n_gas_safe, cf)

            if coupled:
                e_axial, e_face, eps_eq, sigma_cell = self._solve_coupled_field(
                    charged_fluid=cf,
                    coeff=coeff,
                    field=field,
                    h=h,
                    n_gas_safe=n_gas_safe,
                )
            else:
                e_face = None

            species = self._transport_species(
                charged_fluid=cf,
                coeff=coeff,
                e_axial=e_axial,
                e_face=e_face,
                eps=eps,
                h=h,
            )
            max_ampere_residual = max(
                max_ampere_residual, species.max_ampere_residual,
            )

            j_cond_cell = species.j_cond_cell
            jcond_mean_accum += float(np.mean(j_cond_cell)) * (h / dt)
            if coupled:
                v_before, gap_energy = self._advance_gap_node(
                    field=field,
                    h=h,
                    j_cond_cell=j_cond_cell,
                    sigma_cell=sigma_cell,
                )
                gap_v2 += v_before * v_before * h
                gap_e_abs += gap_energy
            # Differential gap conductance for the circuit's semi-implicit load linearisation:
            # only the drift part responds to V_b, and σ = e·Σ_k μ_k n_k (every carrier drifts
            # so as to carry current *with* the field, whatever its charge sign).
            if not coupled:
                sigma_cell = self._conductivity(density, coeff)
            sigma_mean_accum += float(np.mean(sigma_cell)) * (h / dt)

            energy = self._advance_electron_energy(
                n_eps=n_eps,
                charged_fluid=cf,
                neutral_fluid=neutral_fluid,
                coeff=coeff,
                species=species,
                e_axial=e_axial,
                e_face=e_face,
                eps=eps,
                m_gas=m_gas,
                eps_gas_ev=eps_gas_ev,
                h=h,
                gap_voltage_v=v_before if coupled else None,
            )
            n_eps = energy.n_eps
            field_work_j_m3 += energy.field_work_j_m3
            elastic_transfer_j_m3 += energy.elastic_transfer_j_m3
            continuous_sink_j_m3 += energy.continuous_sink_j_m3
            bound_adjustment_j_m3 += energy.bound_adjustment_j_m3
            deposited_j += energy.deposited_j
            maximum_energy_jump_ev = max(
                maximum_energy_jump_ev, energy.maximum_energy_jump_ev,
            )
            limiter_activity += energy.limiter_activity
            max_energy_residual_j_m3 = max(
                max_energy_residual_j_m3, energy.max_energy_residual_j_m3,
            )

            steps_done += 1
            if coupled:
                remaining -= h
                if remaining <= 1.0e-12 * dt:
                    break
                h = self.cfl_safety * self._cfl_dt_target(e_axial, n_gas_safe, me, eps)
                if not np.isfinite(h) or h <= 0.0:
                    h = remaining
                h = min(h, remaining)
                if steps_done >= self.max_substeps:
                    if not capped_warned:
                        warnings.warn(
                            "semi-implicit field coupling hit max_substeps; finishing the "
                            "interval in one capped substep", stacklevel=2,
                        )
                        capped_warned = True
                    h = remaining
            else:
                if steps_done >= n_sub:
                    break
        self._last_substeps = steps_done

        # Write ⟨ε⟩ back to the fluid.
        eps_final = self._recover_eps(n_eps, density[:, e_col], eps_eq)
        cf.electron_mean_energy_ev[:] = eps_final
        self._last_n_eps = np.maximum(n_eps, 0.0)

        return StepReport(
            deposited_power_w=deposited_j / dt if dt > 0 else 0.0,
            mean_electron_energy_ev=float(eps_final.mean()),
            dt=dt,
            # Sign convention: ⟨J_cond⟩ above is measured along **our** mesh axis (anode at x=0,
            # cathode at x=d), but the circuit's node-B KCL (`field.py`,
            # `C_g·V̇_b = i_s − i_f − jcond`) defines jcond as the current flowing **out of node B
            # into the gap** — the discharge-circuit convention, whose axis runs cathode→anode.
            # The two are opposite, so negate: with V_b < 0 the gap must *load* node B (relax V_b
            # toward 0, dissipative I·V_gap ≥ 0), matching `conduction_current_0d`'s G·V_gap.
            conduction_current_a=-self.area_m2 * jcond_mean_accum,
            # Coupled path: energy-matched secant conductance (dissipation == real absorption);
            # frozen path keeps G = A·⟨σ⟩/d — the gap's ohmic (drift) conductance, ≥0.
            gap_conductance_s=(max(gap_e_abs, 0.0) / max(gap_v2, 1e-30)) if coupled
            else self.area_m2 * sigma_mean_accum / self.gap_m,
            gap_energy_absorbed_j=max(gap_e_abs, 0.0) if coupled else None,
            energy_closure=EnergyClosureReport(
                field_work_j_m3=field_work_j_m3,
                elastic_transfer_j_m3=elastic_transfer_j_m3,
                continuous_inelastic_sink_j_m3=continuous_sink_j_m3,
                accepted_reaction_cost_j_m3=np.zeros(cf.n_cells),
                signed_reconciliation_j_m3=np.zeros(cf.n_cells),
                coefficient_bound_adjustment_j_m3=bound_adjustment_j_m3,
                maximum_energy_jump_ev=maximum_energy_jump_ev,
                limiter_activity=limiter_activity,
                max_ampere_residual=max_ampere_residual,
                max_energy_residual_j_m3=max_energy_residual_j_m3,
            ),
        )

    def _advance_via_kernel(self, dt, efield, neutral_fluid, charged_fluid, field) -> StepReport:
        """Coupled-path advance through the numba substep kernel (parity-pinned vs numpy)."""
        from .transport_kernel import coupled_advance, pack_kernel_args

        cf = charged_fluid
        if self._kernel_packed is None or self._kernel_packed[0] is not self._me_coeffs:
            self._kernel_packed = (
                self._me_coeffs,
                pack_kernel_args(self, self._me_coeffs, neutral_fluid, cf),
            )
        result = coupled_advance(
            self, dt, np.asarray(efield, dtype=np.float64)[:, 0], neutral_fluid, cf,
            self._kernel_packed[1], phi_l=field.phi_left, phi_r=field.phi_right,
            c_g=getattr(field, "gap_capacitance_f", None),
            i_drive=getattr(field, "gap_drive_a", 0.0),
        )
        field.v_b_final = result.v_b_final_v
        field.phi_right = result.v_b_final_v
        if result.capped:
            warnings.warn(
                "transport kernel hit max_substeps; its final remainder step exceeded the "
                "requested CFL target", stacklevel=2,
            )
        eps_final = self._recover_eps(
            result.n_eps_j_m3,
            cf.density[:, self._e_col],
            result.equilibrium_energy_ev,
        )
        cf.electron_mean_energy_ev[:] = eps_final
        self._last_n_eps = np.maximum(result.n_eps_j_m3, 0.0)
        self._last_substeps = result.steps
        lef = np.zeros((result.e_axial_v_m.shape[0], 3))
        lef[:, 0] = result.e_axial_v_m
        field.last_efield = lef
        return StepReport(
            deposited_power_w=result.deposited_j / dt if dt > 0 else 0.0,
            mean_electron_energy_ev=float(eps_final.mean()),
            dt=dt,
            conduction_current_a=-self.area_m2 * result.jcond_mean_a_m2,
            # Energy-matched secant conductance: the passive G whose circuit-side dissipation
            # over the interval equals the plasma's actual absorbed energy (the bulk-sigma
            # conductance over-dissipated ~10x into a fictional resistor at the collapse,
            # annihilating the bank energy the plasma never received).
            gap_conductance_s=(
                max(result.gap_energy_absorbed_j, 0.0)
                / max(result.gap_voltage_squared_dt_v2_s, 1e-30)
            ),
            gap_energy_absorbed_j=max(result.gap_energy_absorbed_j, 0.0),
            energy_closure=EnergyClosureReport(
                field_work_j_m3=result.field_work_j_m3,
                elastic_transfer_j_m3=result.elastic_transfer_j_m3,
                continuous_inelastic_sink_j_m3=result.continuous_sink_j_m3,
                accepted_reaction_cost_j_m3=np.zeros(cf.n_cells),
                signed_reconciliation_j_m3=np.zeros(cf.n_cells),
                coefficient_bound_adjustment_j_m3=result.bound_adjustment_j_m3,
                maximum_energy_jump_ev=result.maximum_energy_jump_ev,
                limiter_activity=result.limiter_activity,
                max_ampere_residual=result.max_ampere_residual,
                max_energy_residual_j_m3=result.max_energy_residual_j_m3,
            ),
        )

    # ------------------------------------------------------------------ energy handoff
    def adopt_electron_energy_density(self, energy_j_m3, *, charged_fluid) -> None:
        """Adopt the chemistry integrator's final co-integrated electron-energy state."""
        if not self._resolved:
            return
        cf = charged_fluid
        energy = np.asarray(energy_j_m3, dtype=np.float64)
        if energy.shape != (cf.n_cells,) or np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("energy_j_m3 must be one finite non-negative value per cell")
        n_e = np.maximum(cf.density[:, self._e_col], 0.0)
        eps = np.full(cf.n_cells, self._eps_bounds[0], dtype=np.float64)
        live = n_e > _N_E_FLOOR_M3
        eps[live] = energy[live] / (n_e[live] * E_CHARGE)
        tol = 1.0e-10 * max(1.0, self._eps_bounds[1])
        if np.any(eps[live] < self._eps_bounds[0] - tol) or np.any(
            eps[live] > self._eps_bounds[1] + tol
        ):
            raise ValueError("co-integrated electron energy lies outside coefficient-table bounds")
        eps = np.clip(eps, *self._eps_bounds)
        cf.electron_mean_energy_ev[:] = eps
        self._last_n_eps = n_e * eps * E_CHARGE

    # ------------------------------------------------------------------ energy closure
    def _equilibrium_eps(self, me, e_axial, n_gas_safe) -> np.ndarray | None:
        """Per-cell local-field-equilibrium mean energy ⟨ε⟩_eq(E/N), or ``None`` if unavailable."""
        eq = me.equilibrium_mean_energy(reduced_field_td(np.abs(e_axial), n_gas_safe))
        if eq is None:
            return None
        return np.clip(np.asarray(eq, dtype=np.float64), *self._eps_bounds)

    def _recover_eps(self, n_eps, n_e, eps_eq) -> np.ndarray:
        """Mean energy ⟨ε⟩ from the energy density, repaired where the electron density vanishes.

        LMEA defines ``⟨ε⟩ = n_ε/(nₑ·e)``, which is ill-posed as nₑ → 0: an attachment-emptied
        cell divides ~0 by ~0 and the recovered energy slams between its floor and ceiling,
        corrupting every coefficient looked up on it. Physically, a lone electron simply
        equilibrates to the **local field**, so we blend toward ⟨ε⟩_eq(E/N) with weight
        ``w = nₑ/(nₑ + n_ref)``: the LMEA value where electrons are plentiful, the local-field
        value where they are not. The result is clamped to the coefficient table's own range so
        no lookup ever extrapolates.
        """
        n_e = np.asarray(n_e, dtype=np.float64)
        n_e_safe = np.maximum(n_e, _N_E_FLOOR_M3)
        eps_lmea = np.asarray(n_eps, dtype=np.float64) / (n_e_safe * E_CHARGE)
        if eps_eq is None:
            eps = eps_lmea
        else:
            w = n_e / (n_e + self.lmea_ref_density_m3)
            eps = w * eps_lmea + (1.0 - w) * eps_eq
        eps = np.where(np.isfinite(eps), eps, self._eps_bounds[0])
        return np.clip(eps, *self._eps_bounds)

    # ------------------------------------------------------------------ helpers
    def _evaluate_transport_coefficients(
        self,
        eps: np.ndarray,
        n_gas_safe: np.ndarray,
        charged_fluid: ChargedFluid,
    ) -> _TransportCoefficients:
        """Evaluate electron and ion transport data for the current substep state."""
        me = self._me_coeffs
        mu_e = np.asarray(me.mobility(eps), dtype=np.float64) / n_gas_safe
        d_e = np.asarray(me.diffusion(eps), dtype=np.float64) / n_gas_safe
        mu_eps = np.asarray(me.energy_mobility(eps), dtype=np.float64) / n_gas_safe
        d_eps = np.asarray(me.energy_diffusion(eps), dtype=np.float64) / n_gas_safe

        mu_all: list[np.ndarray] = []
        d_face_all: list[np.ndarray] = []
        for s in range(charged_fluid.n_species):
            if s == self._e_col:
                mu_s = mu_e
                d_s = d_e
            else:
                mu_s = self._ion_mu_N[s] / n_gas_safe
                d_s = mu_s * (K_B * self.gas_temp_k / E_CHARGE)
            mu_all.append(mu_s)
            d_face_all.append(face_average(d_s))
        return _TransportCoefficients(
            mu_e=mu_e,
            d_e=d_e,
            mu_eps=mu_eps,
            d_eps=d_eps,
            mu_all=tuple(mu_all),
            d_face_all=tuple(d_face_all),
        )

    def _conductivity(
        self,
        density: np.ndarray,
        coeff: _TransportCoefficients,
    ) -> np.ndarray:
        """Cell conductivity from the current charged state and frozen substep mobility."""
        sigma_cell = E_CHARGE * density[:, self._e_col] * coeff.mu_e
        for i in self._ion_cols:
            sigma_cell = sigma_cell + E_CHARGE * density[:, i] * coeff.mu_all[i]
        return sigma_cell

    def _solve_coupled_field(
        self,
        *,
        charged_fluid: ChargedFluid,
        coeff: _TransportCoefficients,
        field,
        h: float,
        n_gas_safe: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]:
        """Solve the semi-implicit field using the pre-transport charged state."""
        sigma_cell = self._conductivity(charged_fluid.density, coeff)
        j_expl = self._diffusive_current(
            charged_fluid.density, coeff.d_face_all, charged_fluid.dx,
        )
        e_axial = field.solve(
            charged_fluid.net_charge_density(), sigma_cell, j_expl, h,
        )
        e_face = field.last_e_faces
        eps_eq = self._equilibrium_eps(self._me_coeffs, e_axial, n_gas_safe)
        return e_axial, e_face, eps_eq, sigma_cell

    def _transport_species(
        self,
        *,
        charged_fluid: ChargedFluid,
        coeff: _TransportCoefficients,
        e_axial: np.ndarray,
        e_face: np.ndarray | None,
        eps: np.ndarray,
        h: float,
    ) -> _SpeciesTransportResult:
        """Build SG fluxes, apply electrode fluxes, and advance every charged species."""
        density = charged_fluid.density
        n_cells = charged_fluid.n_cells
        coupled = e_face is not None
        j_left = np.zeros(charged_fluid.n_species)
        j_right = np.zeros(charged_fluid.n_species)
        interior_all: list[np.ndarray] = []
        for s in range(charged_fluid.n_species):
            n_s = density[:, s]
            v_cell = self._z[s] * coeff.mu_all[s] * e_axial
            if coupled:
                v_face = (
                    self._z[s]
                    * 0.5
                    * (coeff.mu_all[s][:-1] + coeff.mu_all[s][1:])
                    * e_face[1:-1]
                )
                v_cell = v_cell.copy()
                v_cell[0] = self._z[s] * coeff.mu_all[s][0] * e_face[0]
                v_cell[-1] = self._z[s] * coeff.mu_all[s][-1] * e_face[-1]
            else:
                v_face = face_average(v_cell)
            interior_all.append(
                sg_interior_flux(n_s, v_face, coeff.d_face_all[s], charged_fluid.dx)
            )
            j_left[s], j_right[s] = self._boundary_fluxes(s, n_s, v_cell, eps)

        ion_cathode_flux = sum(max(0.0, j_right[i]) for i in self._pos_ion_cols)
        j_right[self._e_col] -= self.gamma * ion_cathode_flux

        rho_before = density @ self._q
        j_face_total = np.zeros(n_cells + 1)
        for s in range(charged_fluid.n_species):
            j_face_total[0] += self._q[s] * j_left[s]
            j_face_total[-1] += self._q[s] * j_right[s]
            j_face_total[1:-1] += self._q[s] * interior_all[s]

        gamma_cell: list[np.ndarray] = []
        for s in range(charged_fluid.n_species):
            faces = np.empty(n_cells + 1)
            faces[0] = j_left[s]
            faces[-1] = j_right[s]
            faces[1:-1] = interior_all[s]
            gamma_cell.append(0.5 * (faces[:-1] + faces[1:]))
            dndt = -(faces[1:] - faces[:-1]) / charged_fluid.dx
            density[:, s] = np.maximum(density[:, s] + h * dndt, 0.0)

        charge_rate = (density @ self._q - rho_before) / h
        current_div = (j_face_total[1:] - j_face_total[:-1]) / charged_fluid.dx
        ampere_abs = float(np.max(np.abs(charge_rate + current_div)))
        ampere_scale = max(
            float(np.max(np.abs(charge_rate))),
            float(np.max(np.abs(current_div))),
            1.0e-300,
        )
        j_cond_cell = sum(
            self._q[s] * gamma_cell[s] for s in range(charged_fluid.n_species)
        )
        return _SpeciesTransportResult(
            gamma_e_cell=gamma_cell[self._e_col],
            electron_flux_left=float(j_left[self._e_col]),
            electron_flux_right=float(j_right[self._e_col]),
            j_cond_cell=j_cond_cell,
            max_ampere_residual=ampere_abs / ampere_scale,
        )

    def _advance_gap_node(
        self,
        *,
        field,
        h: float,
        j_cond_cell: np.ndarray,
        sigma_cell: np.ndarray,
    ) -> tuple[float, float]:
        """Advance the finite-capacitance gap and return driving voltage and absorbed energy."""
        v_before = float(field.phi_right)
        i_gap = -self.area_m2 * float(np.mean(j_cond_cell))
        if field.gap_capacitance_f is None:
            return v_before, i_gap * v_before * h

        c_g = float(field.gap_capacitance_f)
        g_sub = self.area_m2 * float(np.mean(sigma_cell)) / self.gap_m
        off = i_gap - g_sub * v_before
        if g_sub > 0.0:
            v_eq = (field.gap_drive_a - off) / g_sub
            v_after = v_eq + (v_before - v_eq) * np.exp(-g_sub * h / c_g)
        else:
            v_after = v_before + (field.gap_drive_a - off) * h / c_g
        field.phi_right = float(v_after)
        field.v_b_final = float(v_after)
        energy_absorbed = 0.5 * c_g * (v_before**2 - v_after**2)
        return v_before, energy_absorbed

    def _advance_electron_energy(
        self,
        *,
        n_eps: np.ndarray,
        charged_fluid: ChargedFluid,
        neutral_fluid: NeutralFluid,
        coeff: _TransportCoefficients,
        species: _SpeciesTransportResult,
        e_axial: np.ndarray,
        e_face: np.ndarray | None,
        eps: np.ndarray,
        m_gas: np.ndarray,
        eps_gas_ev: float,
        h: float,
        gap_voltage_v: float | None,
    ) -> _EnergyTransportResult:
        """Advance energy flux and local closure terms after species and gap updates."""
        dx = charged_fluid.dx
        v_eps_cell = -coeff.mu_eps * e_axial
        if e_face is None:
            v_eps_face = face_average(v_eps_cell)
        else:
            v_eps_face = -0.5 * (coeff.mu_eps[:-1] + coeff.mu_eps[1:]) * e_face[1:-1]
        j_eps_int = sg_interior_flux(
            n_eps, v_eps_face, face_average(coeff.d_eps), dx,
        )
        j_eps_left = eps[0] * E_CHARGE * species.electron_flux_left
        j_eps_right = eps[-1] * E_CHARGE * species.electron_flux_right
        deps_transport = flux_divergence(j_eps_int, j_eps_left, j_eps_right, dx)

        p_field = -E_CHARGE * species.gamma_e_cell * e_axial
        limiter_activity = 0
        if gap_voltage_v is not None:
            pos_tot = float(np.sum(np.maximum(p_field, 0.0))) * dx * self.area_m2
            neg_tot = -float(np.sum(np.minimum(p_field, 0.0))) * dx * self.area_m2
            cap_w = (
                self.field_heating_cap_factor
                * abs(
                    self.area_m2
                    * float(np.mean(species.j_cond_cell))
                    * gap_voltage_v
                )
                + 1.0e3
            )
            magnitude = max(pos_tot, neg_tot)
            if magnitude > cap_w:
                p_field = p_field * (cap_w / magnitude)
                limiter_activity += 1

        n_e = charged_fluid.density[:, self._e_col]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            nu_m = np.where(coeff.mu_e > 0.0, E_CHARGE / (M_E * coeff.mu_e), 0.0)
        p_elastic = (
            (2.0 * M_E / m_gas)
            * nu_m
            * n_e
            * (eps - eps_gas_ev)
            * E_CHARGE
        )
        if self.electron_energy_closure == "continuous_reconciled":
            p_inelastic = self._continuous_inelastic_power(
                eps, neutral_fluid.density, charged_fluid.density,
            ) * E_CHARGE
            local_power = p_field - p_elastic - p_inelastic
        else:
            p_inelastic = np.zeros_like(p_field)
            local_power = np.zeros_like(p_field)

        raw_n_eps = n_eps + h * (deps_transport + local_power)
        lo_bound = self.eps_floor * n_e * E_CHARGE
        hi_bound = self.eps_ceiling * np.maximum(n_e, 0.0) * E_CHARGE + 1.0e-300
        next_n_eps = np.clip(raw_n_eps, lo_bound, hi_bound)
        adjustment = next_n_eps - raw_n_eps
        limiter_activity += int(np.count_nonzero(adjustment != 0.0))
        density_safe = np.maximum(n_e, _N_E_FLOOR_M3)
        maximum_jump = float(
            np.max(np.abs(next_n_eps - n_eps) / (density_safe * E_CHARGE))
        )
        residual = next_n_eps - n_eps - h * (deps_transport + local_power) - adjustment
        return _EnergyTransportResult(
            n_eps=next_n_eps,
            field_work_j_m3=p_field * h,
            elastic_transfer_j_m3=p_elastic * h,
            continuous_sink_j_m3=p_inelastic * h,
            bound_adjustment_j_m3=adjustment,
            deposited_j=float(np.sum(p_field)) * charged_fluid.cell_volume * h,
            maximum_energy_jump_ev=maximum_jump,
            limiter_activity=limiter_activity,
            max_energy_residual_j_m3=float(np.max(np.abs(residual))),
        )

    def _boundary_fluxes(self, s, n_s, v_cell, eps) -> tuple[float, float]:
        """(+x)-directed flux at x=0 and x=L for species ``s`` (drift-limited + thermal wall)."""
        if s == self._e_col:
            v_th0 = _thermal_speed_ev(eps[0], M_E)
            v_thL = _thermal_speed_ev(eps[-1], M_E)
        else:
            kt_ev = K_B * self.gas_temp_k / E_CHARGE
            m = self._ion_masses[s]
            v_th0 = _thermal_speed_ev(kt_ev * 1.5, m)  # (3/2)kT ~ mean energy
            v_thL = v_th0
        n0, nL = float(n_s[0]), float(n_s[-1])
        # Left wall (x=0, outward −x): thermal escape + drift only if drifting toward the wall.
        j_left = -0.25 * n0 * v_th0 + min(0.0, float(v_cell[0])) * n0
        # Right wall (x=L, outward +x).
        j_right = 0.25 * nL * v_thL + max(0.0, float(v_cell[-1])) * nL
        return j_left, j_right

    def _cfl_dt_target(self, e_axial, n_gas_safe, me, eps) -> float:
        """Drift/diffusion CFL time (unscaled min of dx/v_max and dx²/2D_max; may be inf).

        This is the density-independent part of the substep limit — the coupled
        (semi-implicit-field) driver uses it alone, since the implicit field solve removes
        the dielectric-relaxation clause.
        """
        mu_e = np.asarray(me.mobility(eps), dtype=np.float64) / n_gas_safe
        d_e = np.asarray(me.diffusion(eps), dtype=np.float64) / n_gas_safe
        d_eps = np.asarray(me.energy_diffusion(eps), dtype=np.float64) / n_gas_safe
        emag = np.abs(e_axial)
        # Fastest carrier is the energy field (μ_ε ~ 5/3 μ_e); include ion drift for completeness.
        v_max = float(np.max((5.0 / 3.0) * mu_e * emag)) if mu_e.size else 0.0
        for i in self._ion_cols:
            mu_i = self._ion_mu_N[i] / n_gas_safe
            v_max = max(v_max, float(np.max(mu_i * emag)))
        d_max = float(max(np.max(d_e) if d_e.size else 0.0, np.max(d_eps) if d_eps.size else 0.0))
        # Boundary thermal outflux 0.25·n·v_th must respect the CFL too (see the kernel's
        # matching guard): an over-emitting wall cell mints particles through the clip.
        v_max = max(v_max, 0.25 * _thermal_speed_ev(float(eps[0]), M_E),
                    0.25 * _thermal_speed_ev(float(eps[-1]), M_E))
        dx = self._dx
        dt_drift = dx / v_max if v_max > 0.0 else np.inf
        dt_diff = dx * dx / (2.0 * d_max) if d_max > 0.0 else np.inf
        return min(dt_drift, dt_diff)

    def _plan_substeps(self, dt, n_e, e_axial, n_gas_safe, me, eps) -> int:
        """Frozen-field sub-step count: drift/diffusion CFL + dielectric-relaxation time."""
        dt_cfl = self._cfl_dt_target(e_axial, n_gas_safe, me, eps)
        mu_e = np.asarray(me.mobility(eps), dtype=np.float64) / n_gas_safe
        sigma = E_CHARGE * float(np.max(n_e * mu_e)) if n_e.size else 0.0
        tau_diel = EPS0 / sigma if sigma > 0.0 else np.inf
        dt_target = self.cfl_safety * min(dt_cfl, tau_diel)
        n = 1 if (not np.isfinite(dt_target) or dt_target <= 0.0) else int(np.ceil(dt / dt_target))
        n = min(max(1, n), self.max_substeps)
        self._last_substeps = n
        return n

    def _diffusive_current(self, density, d_face_all, dx) -> np.ndarray:
        """Face explicit (diffusion-only) current density Σ_s q_s·Γ_diff,s [A/m²] (N+1,).

        The explicit part of the semi-implicit field predictor: per species the pure-diffusion
        face flux ``−D_face·Δn/dx`` (the SG flux at v = 0, since B(0) = 1), wall faces zero —
        the physical wall fluxes stay in the explicit transport step; the Dirichlet rows
        dominate the boundary in the field solve anyway. Returned at the N+1 faces (the
        cell-centred field solve consumes face currents natively).
        """
        n_cells = density.shape[0]
        jf = np.zeros(n_cells + 1, dtype=np.float64)
        for s in range(density.shape[1]):
            n_s = density[:, s]
            jf[1:-1] += self._q[s] * (-d_face_all[s] * (n_s[1:] - n_s[:-1]) / dx)
        return jf


def _thermal_speed_ev(energy_ev: float, mass_kg: float) -> float:
    """Mean thermal speed √(8·(2/3)ε·e/(π m)) from a mean energy ε (eV) and particle mass."""
    kt_j = (2.0 / 3.0) * max(energy_ev, 0.0) * E_CHARGE
    return float(np.sqrt(8.0 * kt_j / (np.pi * mass_kg)))
