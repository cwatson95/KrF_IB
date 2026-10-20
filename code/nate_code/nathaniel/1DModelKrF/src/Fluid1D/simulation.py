"""1D charged-species fluid, chemistry, discharge circuit, and 0D excimer optics.

``Simulation`` is composed declaratively from the component specs in
:mod:`Fluid1D.components` (1D grid, gas, plasma, cavity, Boltzmann source, discharge
circuit), numeric :mod:`Fluid1D.diagnostics`, and archive-backed plot declarations, then run with a
:class:`~Fluid1D.components.Numerics`. Per step it: advances the discharge circuit for the
electrode potentials, solves the 1D Poisson field, samples the per-cell vec-3 E field,
advances the charged-species :class:`~Fluid1D.plasma_model.PlasmaModel` (electron energy +
deposited power), and — on the fluid cadence — integrates the combined charged+neutral
chemistry (electron-impact rates from the coefficient LUT), and the transverse scalar cavity model
(spontaneous emission, saturated stimulated gain, absorption, and output). MPC runs use the real
drift-diffusion transport; fixed-voltage and analytic stubs remain fast test configurations.
"""

from __future__ import annotations

import sys
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np

import backend
from diagnostics import Plot
from kinetics import FluidChemistry, FluidIntegrator, LaserEnergetics, NeutralFluid
from optics import GainCavity
from particles import EV_TO_J, excimer_system_registry, photon_energy_ev
from reactions import (
    build_cross_section_index,
    import_cross_sections,
    load_ion_mobility,
    load_reaction_network,
    photo_absorbers,
    photo_channels,
    radiative_channels,
    stimulated_channel,
)

from .boltzmann import AnalyticStubSolver, build_lut, en_grid_from_spec, reduced_field_td
from .charged_fluid import ChargedFluid
from .components import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    Preionization,
    ProgressReport,
    Units,
)
from .diagnostics import Diagnostic
from .eedf_solver import MaxwellianEedfSolver
from .field import (
    MpcCircuit,
    Poisson1DSolver,
    SemiImplicitField,
    StubCircuit,
    conduction_current_0d,
    sample_efield,
)
from .mesh import Mesh1D
from .plasma_model import (
    DriftDiffusionPlasmaModel,
    EnergyClosureReport,
    PlasmaModel,
    StepReport,
    StubPlasmaModel,
)
from .progress import ProgressReporter


class Simulation:
    """Declarative 1D fluid excimer simulation with opt-in collectors and plots."""

    def __init__(
        self,
        *,
        name: str,
        grid: Grid1D,
        gas: GasMixture,
        plasma: Plasma,
        cavity: Cavity,
        reactions: str | Path,
        boltzmann: BoltzmannSpec | None = None,
        circuit: CircuitSpec | None = None,
        preionization: Preionization | None = None,
        volume_units: str = "m",
        drive=None,
        diagnostics: list[Diagnostic] | None = None,
        plots: list[Plot] | None = None,
        plasma_model: PlasmaModel | None = None,
        coeff_table=None,
        output_dir: str | Path = "runs",
        cross_sections_dir: str | Path | None = None,
        ion_mobility_path: str | Path | None = None,
    ) -> None:
        self.name = name
        self.grid_spec = grid
        self.gas = gas
        self.plasma_spec = plasma
        self.preionization_spec = preionization
        # User-facing volume/length unit system ("cm" or "m"); geometry + seed densities are
        # interpreted in it (→ SI) at build, and densities/lengths are reported back in it.
        self.units = volume_units if isinstance(volume_units, Units) else Units(volume_units)
        self.cavity = cavity
        self.boltzmann_spec = boltzmann or BoltzmannSpec()
        self.circuit_spec = circuit or CircuitSpec()
        self.drive = drive  # optional provider of prescribed cell and face electric fields
        self.reactions_path = Path(reactions)
        self.cross_sections_dir = (
            Path(cross_sections_dir)
            if cross_sections_dir is not None
            else self.reactions_path.parent.parent / "cross_sections"
        )
        self.ion_mobility_path = (
            Path(ion_mobility_path)
            if ion_mobility_path is not None
            else self.reactions_path.parent.parent / "ion_transport" / "ion_mobility.toml"
        )
        self.diagnostics = list(diagnostics or [])
        self.plots = list(plots or [])
        self._plasma_model = plasma_model
        self._coeff_table = coeff_table
        self.output_dir = Path(output_dir)
        self.numerics: Numerics | None = None
        self.domain_length = grid.length

    # ------------------------------------------------------------------ build
    def _build(self, numerics: Numerics) -> None:
        """Resolve declarative specifications into the mutable runtime model."""
        backend.configure()
        self._build_fluids()
        self._build_chemistry()
        self._build_gain_cavity()
        self._initialize_runtime_state()
        self._build_plasma_transport(numerics)
        self._initialize_energetics()
        self._validate_cavity_absorbers()

    def _build_fluids(self) -> None:
        """Build the SI mesh and seed neutral and charged fluid state."""
        cfg = self.grid_spec
        length_si = self.units.length_to_si(cfg.length)
        self.area_si = self.units.area_to_si(cfg.area_m2)
        self.mesh = Mesh1D(cfg.cells, length_si, self.area_si)
        self.dx = self.mesh.dx
        self.cell_volume = self.mesh.cell_volume
        self.domain_length = length_si

        self.registry = excimer_system_registry()
        neutral_species = tuple(s for s in self.registry if not s.is_charged and not s.is_photon)
        charged_species = tuple(s for s in self.registry if s.is_charged)
        self.fluid = NeutralFluid(neutral_species, cfg.cells, self.dx, self.area_si)
        for name in self.gas.fractions:
            if name not in self.fluid.column:
                raise ValueError(f"GasMixture species {name!r} is not a neutral in the registry")
            self.fluid.seed_uniform(name, self.gas.partial_density(name))

        self.charged = ChargedFluid(charged_species, cfg.cells, self.dx, self.area_si)
        ps = self.plasma_spec
        preion = self.preionization_spec or Preionization(mode="uniform", density=ps.density)
        preion_si = replace(
            preion,
            density=self.units.density_to_si(preion.density),
            amplitude=self.units.density_to_si(preion.amplitude),
            floor=self.units.density_to_si(preion.floor),
            steepness=self.units.inv_length_to_si(preion.steepness),
        )
        n_e0 = preion_si.profile(self.mesh.cell_centers)
        if "e-" in self.charged.column:
            self.charged.seed_profile("e-", n_e0)
        if ps.ion not in self.charged.column:
            raise ValueError(f"Plasma.ion {ps.ion!r} is not a registered charged species")
        ion_charge = self.registry.get(ps.ion).charge_number
        if ion_charge <= 0:
            raise ValueError(f"Plasma.ion {ps.ion!r} must be positively charged")
        self.charged.seed_profile(ps.ion, n_e0 / ion_charge)
        self.charged.seed_electron_energy(ps.electron_mean_energy_ev)
        self._charged_seed = self.charged.density.copy()

    def _build_chemistry(self) -> None:
        """Compile reactions, coefficients, and the chemistry partitions."""
        net = load_reaction_network(self.reactions_path, self.registry)
        self.compiled = net.compile()
        self.cross_section_report = import_cross_sections(self.compiled, self.cross_sections_dir)
        for xs in self.cross_section_report.imports:
            if not xs.ok:
                warnings.warn(f"cross-section: {xs.detail}", stacklevel=2)

        self.coeffs = self._coeff_table or self._build_coeff_table()

        all_species = self.registry.names
        self.plasma_chem = FluidChemistry.for_plasma(self.compiled, self.registry, all_species)
        self.plasma_integ = FluidIntegrator(self.plasma_chem)
        self.plasma_integ.warmup()
        self.electron_chem = FluidChemistry.for_plasma(
            self.compiled, self.registry, all_species, partition="electron_source",
        )
        self.heavy_chem = FluidChemistry.for_plasma(
            self.compiled, self.registry, all_species, partition="heavy_source",
        )
        self.electron_integ = FluidIntegrator(self.electron_chem)
        self.heavy_integ = FluidIntegrator(self.heavy_chem)
        self.electron_integ.warmup()
        self.heavy_integ.warmup()
        all_index = self.plasma_chem.index
        self._fluid_all_cols = np.array([all_index[n] for n in self.fluid.names], dtype=np.int64)
        self._charged_all_cols = np.array([all_index[n] for n in self.charged.names], dtype=np.int64)

    def _build_gain_cavity(self) -> None:
        """Build the scalar cavity from reaction channels and cavity settings."""
        lasing = stimulated_channel(self.compiled, self.registry)
        mode = lasing[3] if lasing else None
        absorbers = [(sp, sig) for sp, sig, m in photo_absorbers(self.compiled, self.registry) if m == mode]
        photo_feedback = [channel for channel in photo_channels(self.compiled, self.registry)
                          if channel[-1] == mode]
        absorbers += list(self.cavity.absorption.items())
        x_f2 = self.gas.normalized_fractions.get("F2", 0.0)
        l_op = min(self.cavity.loss_coeff_per_f2 * x_f2, 0.99)
        beam_area = (self.cavity.beam_area_m2 if self.cavity.beam_area_m2 is not None
                     else self.domain_length * self.area_si / self.cavity.gain_length_m)
        self.gain = GainCavity(
            photon_energy_ev=photon_energy_ev(self.cavity.wavelength_nm),
            radiatives=radiative_channels(self.compiled, self.registry),
            lasing=lasing, absorbers=absorbers,
            gain_length_m=self.cavity.gain_length_m, cavity_length_m=self.cavity.cavity_length_m,
            beam_area_m2=beam_area, r_oc=self.cavity.r_out, l_op=l_op,
            axial_se_fraction=self.cavity.axial_se_fraction, r_back=self.cavity.r_back,
            photo_channels=photo_feedback,
        )

    def _initialize_runtime_state(self) -> None:
        """Initialize field, coupling, chemistry-cadence, and closure ledgers."""
        cfg = self.grid_spec
        self.circuit = self._build_circuit()
        self.poisson = Poisson1DSolver(cfg.cells, self.dx)
        self.cell_centers = self.mesh.cell_centers
        self._conduction_current = 0.0
        self._gap_conductance = 0.0
        self._chem_energy_cost_j = 0.0
        self._chem_energy_bound_adjustment_j = 0.0
        self._continuous_prepaid_j_m3 = np.zeros(self.grid_spec.cells)
        self._last_interval_field_work_j_m3 = np.zeros(self.grid_spec.cells)
        self._last_interval_elastic_j_m3 = np.zeros(self.grid_spec.cells)
        self._last_interval_h_s = 0.0
        self._lead_delivered_j_m3 = np.zeros(self.grid_spec.cells)
        self._last_deliverable_j_m3 = np.zeros(self.grid_spec.cells)
        self._last_heavy_nu = 0.0
        closure_energy_terms = (
            "field_work_j", "elastic_transfer_j", "continuous_inelastic_sink_j",
            "accepted_reaction_cost_j", "signed_reconciliation_j",
            "coefficient_bound_adjustment_j", "photoelectron_birth_energy_j",
            "delivered_source_j", "undelivered_field_work_j",
        )
        self._energy_closure_totals = {name: 0.0 for name in closure_energy_terms}
        self._energy_closure_step = {name: 0.0 for name in closure_energy_terms}
        self._energy_closure_step.update({
            "equilibrium_te_residual_ev": 0.0,
            "bulk_equilibrium_te_residual_ev": 0.0,
            "maximum_source_step_energy_jump_ev": 0.0,
            "limiter_activity": 0.0,
            "max_energy_closure_residual_j_m3": 0.0,
            "max_ampere_residual": 0.0,
        })
        self.last_energy_closure_report = self._empty_energy_closure_report()
        self._chem_efolds_per_substep = 0.25
        self._chem_efolds_faststep = 0.25
        self._chem_max_substeps = 8000
        self._last_chem_step = 0
        self._last_nu_max = 0.0
        self._chem_efolds_per_call = 0.05
        self._exchange_efolds_max = 2.0
        self._exchange_max_subdiv = 256
        self._rates_cache = None

    def _build_plasma_transport(self, numerics: Numerics) -> None:
        """Resolve the injected, drift-diffusion, or deterministic plasma model."""
        self.ion_mobility = self._load_ion_mobility()
        if self._plasma_model is not None:
            self.plasma = self._plasma_model
        elif self.circuit_spec.backend == "mpc":
            if self.coeffs.key_kind != "mean_energy_eV":
                self.coeffs = self.coeffs.to_mean_energy_key()
            self.plasma = DriftDiffusionPlasmaModel(
                reactions=self.compiled, ion_mobility=self.ion_mobility,
                mixture=self.gas.normalized_fractions, gas_temp_k=self.gas.temperature_k,
                area_m2=self.area_si, gap_m=self.domain_length,
                field_heating_cap_factor=numerics.field_heating_cap_factor,
                electron_energy_closure=numerics.electron_energy_closure,
                secondary_emission_gamma=self.plasma_spec.secondary_emission_gamma,
            )
        else:
            self.plasma = StubPlasmaModel()
        self._real_transport = isinstance(self.plasma, DriftDiffusionPlasmaModel)
        if self._real_transport and self.coeffs.key_kind != "mean_energy_eV":
            self.coeffs = self.coeffs.to_mean_energy_key()
        if self._real_transport:
            self.plasma.electron_energy_closure = numerics.electron_energy_closure
        if self._real_transport and not self.plasma._resolved:
            self.plasma._resolve(self.fluid, self.charged, self.coeffs)

    def _initialize_energetics(self) -> None:
        """Initialize run-level power and charged-energy accounting."""
        n_e_total = self.charged.total_number("e-") if "e-" in self.charged.column else 0.0
        e_init = self.plasma_spec.electron_mean_energy_ev * EV_TO_J * n_e_total
        self.energetics = LaserEnergetics(initial_energy_j=e_init)
        self._p_in_accum = 0.0
        self._gain_coeff = 0.0
        self._last_efield = np.zeros((self.grid_spec.cells, 3), dtype=np.float64)

    def _validate_cavity_absorbers(self) -> None:
        """Reject cavity absorbers that are absent from both fluid stores."""
        for name in self.cavity.absorption:
            if name not in self.fluid.column and name not in self.charged.column:
                raise ValueError(
                    f"Cavity.absorption species {name!r} is neither a neutral fluid species "
                    f"nor a registered charged species"
                )

    def _build_coeff_table(self):
        spec = self.boltzmann_spec
        en_grid = en_grid_from_spec(spec.en_grid_td)
        if spec.backend in ("eedf", "two_term"):
            # Coefficients come from measured cross-sections. "eedf" takes Maxwellian-EEDF
            # moments; "two_term" solves the two-term Boltzmann equation per E/N (H&P 2005).
            # Extra energy-loss channels (e.g. F2 vibrational) shape the EEDF/power balance.
            loss_ids = tuple(
                pid for pid in build_cross_section_index(self.cross_sections_dir)
                if "_vib_" in pid
            )
            if spec.backend == "two_term":
                from .two_term import TwoTermBoltzmannSolver

                solver = TwoTermBoltzmannSolver(
                    self.cross_sections_dir, extra_loss_process_ids=loss_ids
                )
            else:
                solver = MaxwellianEedfSolver(
                    self.cross_sections_dir, extra_loss_process_ids=loss_ids
                )
        else:
            solver = AnalyticStubSolver()
        return build_lut(
            self.gas.normalized_fractions, solver=solver, en_grid_td=en_grid,
            gas_temp_k=self.gas.temperature_k, reactions=self.compiled,
            cache_dir=spec.cache_dir, key_kind=spec.key_kind,
        )

    def _build_circuit(self):
        spec = self.circuit_spec
        if spec.backend == "mpc":
            return MpcCircuit(
                gap_m=self.domain_length, area_m2=self.area_si,
                c0_f=spec.c0_f, lm_h=spec.lm_h, cp_f=spec.cp_f, ls_h=spec.ls_h,
                rf_ohm=spec.rf_ohm, lf_h=spec.lf_h, cf_f=spec.cf_f, v0_v=spec.v0_v,
                r1_ohm=spec.r1_ohm, rs_ohm=spec.rs_ohm, c_head_f=spec.c_head_f,
            )
        return StubCircuit(spec.phi_left, spec.phi_right)

    def _load_ion_mobility(self):
        """Load the cited ion reduced-mobility table (None if absent — model uses defaults)."""
        if not self.ion_mobility_path.exists():
            warnings.warn(
                f"ion-mobility table {self.ion_mobility_path} not found; ion transport falls "
                "back to a default reduced mobility", stacklevel=2,
            )
            return None
        return load_ion_mobility(self.ion_mobility_path)

    def _gap_conduction_current(self, v_gap: float) -> float:
        """Return the stub model's 0D closure current ``G·V_gap``."""
        n_gas_mean = float(self.fluid.density.sum(axis=1).mean())
        if n_gas_mean <= 0.0:
            return 0.0
        mu_e = float(self.coeffs.mobility(self._representative_key())) / n_gas_mean
        n_e_mean = float(self.charged.electron_density().mean())
        return conduction_current_0d(
            n_e_mean, mu_e, self.area_si, self.domain_length, v_gap
        )

    # ------------------------------------------------------------------ read interface (diagnostics)
    def electron_temperature(self) -> float:
        """Domain-mean electron temperature Tₑ (eV)."""
        return float(self.charged.electron_temperature_ev().mean())

    def electron_mean_energy(self) -> np.ndarray:
        """Per-cell electron mean energy ⟨ε⟩(x) (eV)."""
        return self.charged.electron_mean_energy()

    def reaction_rates_mean(self) -> tuple[tuple[str, ...], np.ndarray]:
        """Per-reaction volume-mean rates of progress R_r [m⁻³ s⁻¹] at the live state.

        Uses the same combined charged+neutral network and the same per-cell LMEA rate
        constants as the runtime chemistry (``_advance_fluid_and_photons``), so a rate
        sampled here is exactly what the integrator would apply at this instant.
        Returns ``(reaction_ids, rates)`` aligned with ``plasma_chem.reaction_ids``.
        """
        chem = self.plasma_chem
        n_all = np.zeros((self.fluid.density.shape[0], chem.n_species), dtype=np.float64)
        n_all[:, self._fluid_all_cols] = self.fluid.density
        n_all[:, self._charged_all_cols] = self.charged.density
        if self._real_transport and self.coeffs.key_kind == "mean_energy_eV":
            k = chem.dynamic_rates_percell(self.coeffs, self.charged.electron_mean_energy_ev)
        else:
            chem.set_dynamic_rates(self.coeffs, self._representative_key())
            k = np.broadcast_to(chem.k, (n_all.shape[0], chem.n_reactions))
        prod = np.ones((n_all.shape[0], chem.n_reactions), dtype=np.float64)
        for a in range(chem.react_species.shape[1]):
            s = chem.react_species[:, a]
            valid = s >= 0
            prod[:, valid] *= n_all[:, s[valid]] ** chem.react_order[valid, a]
        return chem.reaction_ids, np.asarray((k * prod).mean(axis=0), dtype=np.float64)

    def species_total(self, name: str) -> float:
        if name in self.charged.column:
            return self.charged.total_number(name)
        if name in self.fluid.column:
            return self.fluid.total_number(name)
        if name in self.registry:
            species = self.registry.get(name)
            if species.is_photon and abs((species.wavelength_nm or 0.0) - self.cavity.wavelength_nm) < 1.0:
                return self.gain.total_photons()
        return 0.0

    def gain_coefficient(self) -> float:
        """Net optical gain g = n_up·σ_STE − Σ n_i·σ_i (1/m) from the last optical advance."""
        return self._gain_coeff

    def energy_closure_snapshot(self) -> dict[str, float | str]:
        """Latest outer-step closure terms for the opt-in HDF5 collector."""
        snapshot: dict[str, float | str] = dict(self._energy_closure_step)
        snapshot["closure_mode"] = self.numerics.electron_energy_closure
        if self._real_transport:
            n_gas = self.fluid.density.sum(axis=1)
            eq = self.coeffs.equilibrium_mean_energy(
                reduced_field_td(np.abs(self._last_efield[:, 0]), n_gas)
            )
            if eq is not None:
                residual = (2.0 / 3.0) * (
                    self.charged.electron_mean_energy_ev - np.asarray(eq, dtype=np.float64)
                )
                snapshot["equilibrium_te_residual_ev"] = float(np.max(np.abs(residual)))
                # Density-weighted companion: the max-over-cells metric is dominated by the
                # near-empty sheath cells (huge E/N -> table-ceiling eps_eq); the bulk
                # residual measures how far the CONDUCTING population sits from the
                # local-field balance — the closure-health gate.
                n_e = self.charged.electron_density()
                weight = float(np.sum(n_e))
                snapshot["bulk_equilibrium_te_residual_ev"] = (
                    float(np.sum(np.abs(residual) * n_e) / weight) if weight > 0.0 else 0.0
                )
        return snapshot

    def _reset_step_energy_closure(self) -> None:
        for name in self._energy_closure_step:
            self._energy_closure_step[name] = 0.0
        report = self.last_energy_closure_report
        for values in (
            report.field_work_j_m3,
            report.elastic_transfer_j_m3,
            report.continuous_inelastic_sink_j_m3,
            report.accepted_reaction_cost_j_m3,
            report.signed_reconciliation_j_m3,
            report.coefficient_bound_adjustment_j_m3,
            report.stage_delivered_j_m3,
        ):
            values.fill(0.0)
        report.maximum_energy_jump_ev = 0.0
        report.limiter_activity = 0
        report.max_ampere_residual = 0.0
        report.max_energy_residual_j_m3 = 0.0

    def _empty_energy_closure_report(self) -> EnergyClosureReport:
        def zeros() -> np.ndarray:
            return np.zeros(self.grid_spec.cells, dtype=np.float64)

        return EnergyClosureReport(
            field_work_j_m3=zeros(),
            elastic_transfer_j_m3=zeros(),
            continuous_inelastic_sink_j_m3=zeros(),
            accepted_reaction_cost_j_m3=zeros(),
            signed_reconciliation_j_m3=zeros(),
            coefficient_bound_adjustment_j_m3=zeros(),
            stage_delivered_j_m3=zeros(),
        )

    def _book_closure_energy(self, name: str, value_j: float) -> None:
        value = float(value_j)
        self._energy_closure_step[name] += value
        self._energy_closure_totals[name] += value

    def _book_transport_closure(self, report: StepReport) -> None:
        closure = report.energy_closure
        if closure is None:
            return
        volume = self.cell_volume
        self._book_closure_energy("field_work_j", np.sum(closure.field_work_j_m3) * volume)
        self._book_closure_energy(
            "elastic_transfer_j", np.sum(closure.elastic_transfer_j_m3) * volume,
        )
        self._book_closure_energy(
            "continuous_inelastic_sink_j",
            np.sum(closure.continuous_inelastic_sink_j_m3) * volume,
        )
        self._book_closure_energy(
            "coefficient_bound_adjustment_j",
            np.sum(closure.coefficient_bound_adjustment_j_m3) * volume,
        )
        self._continuous_prepaid_j_m3 += closure.continuous_inelastic_sink_j_m3
        # Capture the JUST-COMPLETED interval's per-cell heating ledgers for the stage
        # delivery (copies: the live report is a per-outer-step accumulator, wrong under
        # climax subdivision, and the advance report itself is discarded by the caller).
        self._last_interval_field_work_j_m3 = closure.field_work_j_m3.copy()
        self._last_interval_elastic_j_m3 = closure.elastic_transfer_j_m3.copy()
        self._last_interval_h_s = float(report.dt)
        live = self.last_energy_closure_report
        live.field_work_j_m3 += closure.field_work_j_m3
        live.elastic_transfer_j_m3 += closure.elastic_transfer_j_m3
        live.continuous_inelastic_sink_j_m3 += closure.continuous_inelastic_sink_j_m3
        live.coefficient_bound_adjustment_j_m3 += closure.coefficient_bound_adjustment_j_m3
        live.maximum_energy_jump_ev = max(
            live.maximum_energy_jump_ev, closure.maximum_energy_jump_ev,
        )
        live.limiter_activity += closure.limiter_activity
        live.max_energy_residual_j_m3 = max(
            live.max_energy_residual_j_m3, closure.max_energy_residual_j_m3,
        )
        live.max_ampere_residual = max(live.max_ampere_residual, closure.max_ampere_residual)
        self._energy_closure_step["maximum_source_step_energy_jump_ev"] = max(
            self._energy_closure_step["maximum_source_step_energy_jump_ev"],
            closure.maximum_energy_jump_ev,
        )
        self._energy_closure_step["limiter_activity"] += closure.limiter_activity
        self._energy_closure_step["max_energy_closure_residual_j_m3"] = max(
            self._energy_closure_step["max_energy_closure_residual_j_m3"],
            closure.max_energy_residual_j_m3,
        )
        self._energy_closure_step["max_ampere_residual"] = max(
            self._energy_closure_step["max_ampere_residual"], closure.max_ampere_residual,
        )

    def _deliverable_heating_j_m3(self) -> np.ndarray:
        """Rectified per-cell stage heating of the last transport interval (J/m³, ≥ 0).

        The signed (field_work − elastic) ledger books large phantom NEGATIVE work in the
        one-cell electrode falls (sub-Debye dx: the cell-centred Γ_e·E sign is unresolved
        there — measured to cancel the entire bulk heating post-breakdown). Physically
        deliverable heating is its positive part; the discarded negative part stays visible
        in the ``undelivered_field_work_j`` ledger and in the signed ``field_work`` books.
        """
        signed = (self._last_interval_field_work_j_m3
                  - self._last_interval_elastic_j_m3)
        deliverable = np.maximum(signed, 0.0)
        self._book_closure_energy(
            "undelivered_field_work_j",
            float(np.sum(signed - deliverable)) * self.cell_volume,
        )
        return deliverable

    def _source_exchange(self, t: float, interval: float) -> None:
        """Strang-split local electron source around one field/SG transport interval.

        The stage's power source is transport's ledgered field work minus elastic exchange,
        rectified per cell (see ``_deliverable_heating_j_m3``): the LEAD half-step uses the
        previous interval's mean power (lagged predictor, zeroed on the first interval and
        across climax-subdivision h jumps where an unguarded lead would deliver one
        spike-power interval), and the TRAIL half-step delivers the exact residual — so the
        energy delivered over each interval equals that interval's rectified ledger
        identically, and heating coexists with the stiff inelastic sink inside the stage
        instead of being drained from a stored buffer.
        """
        half = 0.5 * interval
        h_prev = self._last_interval_h_s
        if h_prev > 0.0 and interval <= 2.0 * h_prev:
            p_lead = self._last_deliverable_j_m3 / h_prev
        else:
            p_lead = np.zeros(self.grid_spec.cells)
        self._lead_delivered_j_m3 = p_lead * half
        self._advance_electron_source(half, power_source_j_m3_s=p_lead)
        self._exchange(t, interval)
        self._last_deliverable_j_m3 = self._deliverable_heating_j_m3()
        p_trail = (self._last_deliverable_j_m3 - self._lead_delivered_j_m3) / half
        self._advance_electron_source(half, power_source_j_m3_s=p_trail)

    # ------------------------------------------------------------------ step
    def _step(self, i: int) -> None:
        num = self.numerics
        dt = num.dt
        self._reset_step_energy_closure()

        # Climax subdivision: once the avalanche e-folds faster than the exchange cadence,
        # the circuit↔plasma handshake (I_cond, G → node-B linearisation) at dt granularity
        # feeds currents that are stale by whole e-foldings — node B then rings unphysically
        # and the books explode. Subdivide the WHOLE exchange (circuit + field + transport +
        # chemistry) so its granularity follows the physics through the spike.
        k = 1
        if self._last_nu_max > 0.0:
            need = self._last_nu_max * dt / self._exchange_efolds_max
            if need > 1.0:
                k = int(min(2.0 ** np.ceil(np.log2(need)), self._exchange_max_subdiv))
        if k == 1:
            if num.electron_energy_closure == "source_coupled" and self._real_transport:
                self._source_exchange(i * dt, dt)
            else:
                self._exchange(i * dt, dt)
            if i > self._last_chem_step and self._chemistry_due(i):
                self._advance_fluid_and_photons(i * dt, interval=(i - self._last_chem_step) * dt)
                self._last_chem_step = i
        else:
            h = dt / k
            t0 = (i - 1) * dt
            for j in range(1, k + 1):
                if num.electron_energy_closure == "source_coupled" and self._real_transport:
                    self._source_exchange(t0 + j * h, h)
                else:
                    self._exchange(t0 + j * h, h)
                    self._advance_fluid_and_photons(t0 + j * h, interval=h,
                                                    reuse_rates=(j > 1))
            if num.electron_energy_closure == "source_coupled" and self._real_transport:
                if i > self._last_chem_step and self._chemistry_due(i):
                    self._advance_fluid_and_photons(
                        i * dt, interval=(i - self._last_chem_step) * dt,
                    )
                    self._last_chem_step = i
            else:
                self._last_chem_step = i

    def _exchange(self, t: float, h: float) -> None:
        """One circuit → field → charged-transport handshake over ``h``."""
        # Circuit advances on the plasma load; returns the electrode Dirichlet potentials
        # (anode=0 at x=0, cathode=V_b at x=d). The load is handed over as the PASSIVE
        # conductance form I = G·V_b (offset-free). This keeps the linearized load passive
        # and dissipative, matching Chen §VI's closure; the non-ohmic diffusive contribution
        # is negligible in the circuit frame.
        # For real MPC transport, a Lie split assigns the network subflow to MpcCircuit and
        # the gap-draw subflow to transport; node B is then pinned to the latter's final value.
        split_gap = self._real_transport and isinstance(self.circuit, MpcCircuit)
        v_b_now = getattr(self.circuit, "v_b", None)
        if split_gap:
            # Lie split: first advance the MPC network and gap capacitor with the plasma open;
            # transport then advances the same node under gap draw only and pins it back.
            i_load, g_load = 0.0, 0.0
        elif self._real_transport and v_b_now is not None:
            i_load, g_load = self._gap_conductance * v_b_now, self._gap_conductance
        else:
            i_load, g_load = self._conduction_current, self._gap_conductance
        phi_l, phi_r = self.circuit.advance(
            h, discharge_current_a=i_load, gap_conductance_s=g_load,
        )
        kernel_path = (self._real_transport
                       and getattr(self.plasma, "use_numba_kernel", False)
                       and self.drive is None)
        if kernel_path:
            # The kernel re-solves the semi-implicit field internally every substep; the
            # previous final field supplies the bootstrap and CFL-sizing input.
            efield = self._last_efield
        else:
            rho = self.charged.net_charge_density()
            _phi, e1d = self.poisson.solve(rho, phi_l, phi_r)
            efield = sample_efield(e1d, self.cell_centers, axis=0, prescribed=self.drive, t=t)
            self._last_efield = efield

        # Real transport re-solves the semi-implicit σ-corrected field inside its own substep
        # loop (stable past breakdown); the plain solve above bootstraps it and feeds the stub.
        handle = None
        if self._real_transport:
            handle = SemiImplicitField(
                solver=self.poisson, phi_left=phi_l, phi_right=phi_r,
                cell_centers=self.cell_centers, prescribed=self.drive, t=t,
                gap_capacitance_f=getattr(self.circuit, "cg_eff", None),
                gap_drive_a=0.0 if split_gap else (
                    getattr(self.circuit, "i_s", 0.0) - getattr(self.circuit, "i_f", 0.0)
                ),
            )
        report = self.plasma.advance(
            h, efield=efield, neutral_fluid=self.fluid, charged_fluid=self.charged,
            coeffs=self.coeffs, field=handle,
        )
        self._book_transport_closure(report)
        if handle is not None and handle.last_efield is not None:
            self._last_efield = handle.last_efield  # diagnostics see the final-substep field
        if split_gap and handle is not None and handle.v_b_final is not None:
            self.circuit.pin_gap_voltage(handle.v_b_final)
        # Deposited-power ledger: the conduction power |I_gap*V_gap| is the honest electrical
        # deposition that closes against the circuit. The signed field-heating integral is
        # retained separately in the energy-closure ledger.
        if self._real_transport:
            if report.gap_energy_absorbed_j is not None:
                self._p_in_accum += max(report.gap_energy_absorbed_j, 0.0)
            else:
                self._p_in_accum += abs(report.conduction_current_a * phi_r) * h
        else:
            self._p_in_accum += report.deposited_power_w * h
        # Conduction current the circuit sees next step: the real 1D flux integral
        # A·⟨J_cond⟩ from the transport model, or the 0D closure G·V_gap for the stub.
        if self._real_transport:
            self._conduction_current = report.conduction_current_a
            self._gap_conductance = report.gap_conductance_s
        else:
            self._conduction_current = self._gap_conduction_current(phi_r)

    def _representative_key(self) -> float:
        """Domain-representative lookup key in the coefficient table's declared units."""
        if self.coeffs.key_kind == "mean_energy_eV":
            return float(np.mean(self.charged.electron_mean_energy_ev))
        n_gas = self.fluid.density.sum(axis=1)
        emag = np.linalg.norm(self._last_efield, axis=1)
        key = reduced_field_td(emag, n_gas)
        return float(np.mean(key))

    def _chemistry_due(self, i: int) -> bool:
        """Chemistry cadence: every ``fluid_every`` steps, tightened to every step at breakdown.

        Between calls the electron-impact rates k(⟨ε⟩) are frozen, so once the last-seen
        electron-multiplying frequency ν makes the nominal interval span more than a fraction
        of an e-folding, the cadence drops to the transport step — the stale-k window then
        matches the cadence the energy PDE (and the field) already run at.
        """
        num = self.numerics
        if i % num.fluid_every == 0:
            return True
        if (getattr(self, "_real_transport", False)
                and getattr(num, "electron_energy_closure", "continuous_reconciled")
                == "source_coupled"):
            # Electron rows already run every transport interval; grade the HEAVY cadence
            # on its own charged-changing stiffness (ion-ion recombination, dimerization
            # reach ν ~ 1e9 s⁻¹ during conduction — a fixed 1 ns cadence lands them as one
            # impulsive kick per call).
            nu = self._last_heavy_nu
            if nu <= 0.0:
                return False
            return (i - self._last_chem_step) >= max(
                1.0, self._chem_efolds_per_call / (nu * num.dt),
            )
        # Graded cadence: target <= _chem_efolds_per_call electron e-folds per call. The
        # binary 20-or-1 trigger left the conducting phase (nu ~ 0.1 e-folds/ns, below the
        # 0.25 in-call threshold) on the 1 ns cadence, where each call's accumulated
        # attachment landed as one impulsive kick (the F2-sweep sawtooth).
        nu = self._last_nu_max
        if nu <= 0.0:
            return False
        steps_target = self._chem_efolds_per_call / (nu * num.dt)
        return (i - self._last_chem_step) >= max(1.0, steps_target)

    def _chemistry_substeps(self, interval, k_percell, n_all, floor_substeps) -> int:
        """Substeps needed to resolve the fastest reaction timescale over ``interval``.

        The stiff integrator is L-stable, so a coarse step never oscillates — but it is only
        first-order, so a step spanning many e-foldings of an avalanche saturates to "everything
        reacted" instead of tracking the growth. At breakdown the ionization frequency
        ``ν ≈ k·n_target`` reaches ~1e12 s⁻¹ while the fluid interval is ~1 ns, i.e. ~1e3
        e-foldings per call: a fixed substep count cannot follow that. Size the count from
        ``ν_max`` instead, keeping each substep to a fraction of an e-folding.
        """
        # Cadence stiffness watches every electron-CHANGING channel (ionization AND
        # attachment): watching only multipliers let the conducting phase relax to the 1 ns
        # cadence while attachment kept running -- each call then kicked n_e/sigma/V_b
        # impulsively, which creates a phase-locked F2-sweep sawtooth.
        mask = getattr(self.plasma_chem, "electron_changing_mask", None)
        if mask is None:
            mask = self.plasma_chem.electron_multiplying_mask
        if mask is None or not mask.any():
            self._last_nu_max = 0.0
            return floor_substeps
        nu_max = self.plasma_chem.progress_frequency(n_all, k_percell, mask)
        if not np.isfinite(nu_max) or nu_max <= 0.0:
            self._last_nu_max = 0.0
            return floor_substeps
        self._last_nu_max = float(nu_max)  # drives the adaptive call cadence (_chemistry_due)
        # At the per-step cadence the stale-k window already equals dt, and the L-stable
        # integrator is accurate to ~6 % at 1.5 e-foldings — resolving 0.25 e-foldings on top
        # of that is double-paying (it was the chemistry cost driver at breakdown). Keep the
        # fine target only for the long (fluid_every) intervals.
        efolds = (self._chem_efolds_per_substep
                  if interval > 2.0 * self.numerics.dt else self._chem_efolds_faststep)
        needed = int(np.ceil(interval * nu_max / efolds))
        n = int(np.clip(needed, floor_substeps, self._chem_max_substeps))
        return n

    def _combined_densities(self) -> np.ndarray:
        """Gather neutral and charged fluids into the all-species chemistry ordering."""
        n_all = np.zeros(
            (self.fluid.density.shape[0], self.plasma_chem.n_species), dtype=np.float64,
        )
        n_all[:, self._fluid_all_cols] = self.fluid.density
        n_all[:, self._charged_all_cols] = self.charged.density
        return n_all

    def _adopt_combined_densities(self, n_all: np.ndarray) -> None:
        self.fluid.density[:] = n_all[:, self._fluid_all_cols]
        self.charged.density[:] = n_all[:, self._charged_all_cols]

    def _book_accepted_source(self, integrator: FluidIntegrator) -> None:
        cost_ev = integrator.last_energy_cost_ev_m3
        if cost_ev is not None:
            cost_j = float(np.sum(cost_ev)) * self.cell_volume * EV_TO_J
            self._chem_energy_cost_j += cost_j
            self._book_closure_energy("accepted_reaction_cost_j", cost_j)
            self.last_energy_closure_report.accepted_reaction_cost_j_m3 += cost_ev * EV_TO_J
        adjustment_ev = integrator.last_energy_bound_adjustment_ev_m3
        if adjustment_ev is not None:
            adjustment_j = float(np.sum(adjustment_ev)) * self.cell_volume * EV_TO_J
            self._chem_energy_bound_adjustment_j += adjustment_j
            self._book_closure_energy("coefficient_bound_adjustment_j", adjustment_j)
            self.last_energy_closure_report.coefficient_bound_adjustment_j_m3 += (
                adjustment_ev * EV_TO_J
            )
        self._energy_closure_step["maximum_source_step_energy_jump_ev"] = max(
            self._energy_closure_step["maximum_source_step_energy_jump_ev"],
            integrator.last_maximum_energy_jump_ev,
        )
        self._energy_closure_step["limiter_activity"] += integrator.last_limiter_activity

    def _advance_electron_source(
        self, interval: float, power_source_j_m3_s: np.ndarray | None = None,
    ) -> None:
        """Advance production electron reactions and their debit on the transport cadence.

        ``power_source_j_m3_s`` (per cell, W/m³, signed) is the transport-ledgered local
        heating delivered inside the stage's substep loop, so field work and the stiff
        inelastic sink share one balance point.
        """
        if interval <= 0.0 or not self._real_transport:
            return
        n_all = self._combined_densities()
        keys = self.charged.electron_mean_energy_ev
        k_percell = self.electron_chem.dynamic_rates_percell(self.coeffs, keys)
        mask = self.electron_chem.electron_changing_mask
        nu_max = (
            self.electron_chem.progress_frequency(n_all, k_percell, mask)
            if mask is not None and mask.any() else 0.0
        )
        requested = self.electron_chem.energy_consumption_rate_percell(n_all, k_percell)
        stored_ev = np.maximum(self.plasma._last_n_eps / EV_TO_J, 1.0)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            energy_frequency = float(np.max(np.abs(requested) / stored_ev))
        source_frequency = max(float(nu_max), energy_frequency)
        if power_source_j_m3_s is not None:
            # Resolve the delivered heating too — against the PHYSICAL energy scale
            # eps_ceiling·n_e rather than the stored buffer: field work does not scale
            # with n_e, so trace-density sheath cells would otherwise pin the global
            # substep count at the cap every half-stage (the clamp window bounds their
            # attainable change anyway).
            n_e = np.maximum(n_all[:, self.electron_chem.index["e-"]], 1.0e10)
            denom_ev = np.maximum(stored_ev, self.plasma._eps_bounds[1] * n_e)
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                source_frequency = max(
                    source_frequency,
                    float(np.max(np.abs(power_source_j_m3_s) / EV_TO_J / denom_ev)),
                )
        # Resolve predicted source-step changes to at most 5% in either electron population
        # or stored electron energy.  The accepted augmented stage enforces the floor exactly.
        needed = max(1, int(np.ceil(interval * source_frequency / 0.05)))
        substeps = min(needed, self._chem_max_substeps)
        self._last_nu_max = float(nu_max)
        self.electron_integ.advance_energy_coupled(
            n_all,
            interval,
            substeps=substeps,
            electron_energy_j_m3=self.plasma._last_n_eps,
            electron_col=self.electron_chem.index["e-"],
            coeffs=self.coeffs,
            energy_bounds_ev=self.plasma._eps_bounds,
            power_source_j_m3_s=power_source_j_m3_s,
        )
        self._adopt_combined_densities(n_all)
        energy_j = self.electron_integ.last_electron_energy_j_m3
        if energy_j is None:
            raise RuntimeError("electron-source integrator did not return electron energy")
        self.plasma.adopt_electron_energy_density(energy_j, charged_fluid=self.charged)
        self._book_accepted_source(self.electron_integ)
        if power_source_j_m3_s is not None:
            delivered_j_m3 = power_source_j_m3_s * interval
            self._book_closure_energy(
                "delivered_source_j", float(np.sum(delivered_j_m3)) * self.cell_volume,
            )
            self.last_energy_closure_report.stage_delivered_j_m3 += delivered_j_m3

    def _advance_source_coupled_interval(
        self,
        n_all: np.ndarray,
        interval: float,
    ) -> None:
        """Advance the heavy partition after transport-cadence electron sources."""
        num = self.numerics
        self.heavy_integ.advance(n_all, interval, substeps=num.fluid_substeps)
        mask = getattr(self.heavy_chem, "charged_changing_mask", None)
        if mask is None or not mask.any():
            self._last_heavy_nu = 0.0
            return
        k_heavy = np.broadcast_to(
            self.heavy_chem.k,
            (n_all.shape[0], self.heavy_chem.n_reactions),
        )
        self._last_heavy_nu = float(
            self.heavy_chem.progress_frequency(n_all, k_heavy, mask)
        )

    def _interval_rate_table(self, reuse_rates: bool) -> np.ndarray:
        """Return per-cell chemistry rates for the current electron energy."""
        if reuse_rates and self._rates_cache is not None:
            return self._rates_cache.copy()
        keys = self.charged.electron_mean_energy_ev
        k_percell = self.plasma_chem.dynamic_rates_percell(self.coeffs, keys)
        self._rates_cache = k_percell.copy()
        return k_percell

    def _rates_for_substep_sizing(
        self,
        n_all: np.ndarray,
        k_percell: np.ndarray,
        interval: float,
        prepaid_j_m3: np.ndarray,
    ) -> np.ndarray:
        """Limit only the rate estimate used to choose chemistry substeps."""
        costly = self.plasma_chem.energy_costly_mask
        if costly is None or not costly.any():
            return k_percell
        requested = self.plasma_chem.energy_consumption_rate_percell(n_all, k_percell)
        electron_col = self.plasma_chem.index["e-"]
        reserve_ev = self.plasma._eps_bounds[0] * n_all[:, electron_col]
        available_ev = np.maximum(
            (self.plasma._last_n_eps + prepaid_j_m3) / EV_TO_J - reserve_ev,
            0.0,
        )
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            factor = np.where(
                requested * interval > available_ev,
                available_ev / np.maximum(requested * interval, 1.0e-300),
                1.0,
            )
        factor = np.clip(np.nan_to_num(factor, nan=0.0, posinf=1.0), 0.0, 1.0)
        k_for_sizing = k_percell.copy()
        k_for_sizing[:, costly] *= factor[:, np.newaxis]
        return k_for_sizing

    def _advance_reconciled_interval(
        self,
        n_all: np.ndarray,
        interval: float,
        reuse_rates: bool,
    ) -> np.ndarray:
        """Advance continuous-reconciled chemistry and return its prepaid sink."""
        num = self.numerics
        k_percell = self._interval_rate_table(reuse_rates)
        prepaid_j_m3 = self._continuous_prepaid_j_m3.copy()
        k_for_sizing = self._rates_for_substep_sizing(
            n_all,
            k_percell,
            interval,
            prepaid_j_m3,
        )
        substeps = self._chemistry_substeps(
            interval,
            k_for_sizing,
            n_all,
            num.fluid_substeps,
        )
        self.plasma_integ.advance_energy_coupled(
            n_all,
            interval,
            substeps=substeps,
            electron_energy_j_m3=self.plasma._last_n_eps.copy(),
            electron_col=self.plasma_chem.index["e-"],
            coeffs=self.coeffs,
            energy_bounds_ev=self.plasma._eps_bounds,
            power_source_j_m3_s=(
                prepaid_j_m3 / interval if interval > 0.0 else None
            ),
        )
        return prepaid_j_m3

    def _advance_stub_interval(self, n_all: np.ndarray, interval: float) -> None:
        """Advance chemistry using one representative coefficient-table key."""
        self.plasma_chem.set_dynamic_rates(self.coeffs, self._representative_key())
        self.plasma_integ.advance(
            n_all,
            interval,
            substeps=self.numerics.fluid_substeps,
        )

    def _adopt_reconciled_chemistry(self, prepaid_j_m3: np.ndarray) -> None:
        """Adopt the coupled electron energy and book its closure terms."""
        energy_j = self.plasma_integ.last_electron_energy_j_m3
        if energy_j is None:
            raise RuntimeError("energy-coupled chemistry did not return electron energy")
        self.plasma.adopt_electron_energy_density(energy_j, charged_fluid=self.charged)

        cost_ev = self.plasma_integ.last_energy_cost_ev_m3
        if cost_ev is not None:
            cost_j = float(np.sum(cost_ev)) * self.cell_volume * EV_TO_J
            self._chem_energy_cost_j += cost_j
            self._book_closure_energy("accepted_reaction_cost_j", cost_j)
            reconciliation_j_m3 = prepaid_j_m3 - cost_ev * EV_TO_J
            self._book_closure_energy(
                "signed_reconciliation_j",
                float(np.sum(reconciliation_j_m3)) * self.cell_volume,
            )

        adjustment_ev = self.plasma_integ.last_energy_bound_adjustment_ev_m3
        if adjustment_ev is not None:
            adjustment_j = float(np.sum(adjustment_ev)) * self.cell_volume * EV_TO_J
            self._chem_energy_bound_adjustment_j += adjustment_j
            self._book_closure_energy("coefficient_bound_adjustment_j", adjustment_j)

        self._energy_closure_step["maximum_source_step_energy_jump_ev"] = max(
            self._energy_closure_step["maximum_source_step_energy_jump_ev"],
            self.plasma_integ.last_maximum_energy_jump_ev,
        )
        self._energy_closure_step["limiter_activity"] += self.plasma_integ.last_limiter_activity

        live = self.last_energy_closure_report
        if cost_ev is not None:
            live.accepted_reaction_cost_j_m3 += cost_ev * EV_TO_J
            live.signed_reconciliation_j_m3 += reconciliation_j_m3
        if adjustment_ev is not None:
            live.coefficient_bound_adjustment_j_m3 += adjustment_ev * EV_TO_J
        live.maximum_energy_jump_ev = max(
            live.maximum_energy_jump_ev,
            self.plasma_integ.last_maximum_energy_jump_ev,
        )
        live.limiter_activity += self.plasma_integ.last_limiter_activity
        self._continuous_prepaid_j_m3.fill(0.0)

    def _adopt_interval_state(
        self,
        n_all: np.ndarray,
        *,
        source_coupled: bool,
        prepaid_j_m3: np.ndarray | None,
    ) -> None:
        """Scatter chemistry state back into neutral and charged fluid stores."""
        self.fluid.density[:] = n_all[:, self._fluid_all_cols]
        if not self._real_transport:
            self.charged.density[:] = self._charged_seed
            return
        self.charged.density[:] = n_all[:, self._charged_all_cols]
        if source_coupled:
            return
        if prepaid_j_m3 is None:
            raise RuntimeError("continuous-reconciled chemistry has no prepaid energy")
        self._adopt_reconciled_chemistry(prepaid_j_m3)

    def _resync_transport_energy_after_optics(self) -> None:
        """Synchronize transport electron energy after photo reactions."""
        n_e = np.maximum(self.charged.density[:, self.charged.column["e-"]], 0.0)
        eps = np.clip(self.charged.electron_mean_energy_ev, *self.plasma._eps_bounds)
        clipped_j = float(np.sum(
            (eps - self.charged.electron_mean_energy_ev) * n_e,
        )) * EV_TO_J * self.cell_volume
        if clipped_j != 0.0:
            self._book_closure_energy("coefficient_bound_adjustment_j", clipped_j)
        self.charged.electron_mean_energy_ev[:] = eps
        self.plasma._last_n_eps = n_e * eps * EV_TO_J

    def _advance_optics(self, t: float, interval: float) -> None:
        """Advance the scalar cavity and record interval-average input/output power."""
        p_out, self._gain_coeff = self.gain.advance(
            interval,
            fluid=self.fluid,
            charged=self.charged,
        )
        self._book_closure_energy(
            "photoelectron_birth_energy_j",
            self.gain.last_photoelectron_birth_energy_j
            * self.cell_volume
            * self.grid_spec.cells,
        )
        if self._real_transport:
            self._resync_transport_energy_after_optics()
        p_in = self._p_in_accum / interval if interval > 0 else 0.0
        self.energetics.record(t, 0.0, p_out, p_in, interval_s=interval)
        self._p_in_accum = 0.0

    def _advance_fluid_and_photons(
        self,
        t: float,
        interval: float | None = None,
        reuse_rates: bool = False,
    ) -> None:
        """Advance cadence-level chemistry, optics, and energy accounting."""
        num = self.numerics
        if interval is None:
            interval = num.fluid_every * num.dt

        n_all = self._combined_densities()
        source_coupled = (
            self._real_transport and num.electron_energy_closure == "source_coupled"
        )
        prepaid_j_m3: np.ndarray | None = None
        if source_coupled:
            self._advance_source_coupled_interval(n_all, interval)
        elif self._real_transport and self.coeffs.key_kind == "mean_energy_eV":
            prepaid_j_m3 = self._advance_reconciled_interval(
                n_all,
                interval,
                reuse_rates,
            )
        else:
            self._advance_stub_interval(n_all, interval)
        self._adopt_interval_state(
            n_all,
            source_coupled=source_coupled,
            prepaid_j_m3=prepaid_j_m3,
        )
        self._advance_optics(t, interval)

    # ------------------------------------------------------------------ run
    def run(
        self,
        numerics: Numerics | None = None,
        *,
        progress: ProgressReport | None = None,
        verbose: bool = True,
    ) -> Path:
        """Build and advance the run; return its output directory (``output_dir/name``).

        ``verbose`` gates the startup banner, the default-cadence progress lines, and the
        closing energy summary. ``progress`` customizes the progress line (cadence, style,
        fields — see :class:`~Fluid1D.components.ProgressReport`) and reports even when
        ``verbose`` is False, so quiet embedders can keep a heartbeat.
        """
        self.numerics = numerics or Numerics()
        self._build(self.numerics)
        # An explicit spec always reports, even on quiet runs; ``verbose`` alone reports with
        # the default fields. Fully silent runs skip reporter construction/validation.
        reporter = (
            ProgressReporter(progress or ProgressReport(), self)
            if verbose or progress is not None
            else None
        )
        out = self.output_dir / self.name
        out.mkdir(parents=True, exist_ok=True)
        if reporter is not None:
            reporter.start(out)
        if verbose:
            n_e = self.charged.total_number("e-") if "e-" in self.charged.column else 0.0
            print(f"{self.name}: {self.grid_spec.cells}-cell 1D mesh, "
                  f"{n_e:.3e} electrons, {self.numerics.steps} steps @ dt={self.numerics.dt:.1e}s")
            closure_label = self.numerics.electron_energy_closure
            if closure_label == "continuous_reconciled":
                closure_label += " (diagnostic only)"
            print(f"  electron-energy closure: {closure_label}")
            for xs in self.cross_section_report.imports:
                print(f"  {xs.describe()}")
            # Land the banner before the (possibly long) first JIT'd step, even when
            # stdout is a block-buffered pipe (nohup/tee).
            sys.stdout.flush()
        for i in range(1, self.numerics.steps + 1):
            self._step(i)
            # Flush a partial final cadence. Short runs and a steps % fluid_every remainder
            # must still advance chemistry, photons, and the accumulated energy ledger.
            if i == self.numerics.steps and i > self._last_chem_step:
                interval = (i - self._last_chem_step) * self.numerics.dt
                self._advance_fluid_and_photons(i * self.numerics.dt, interval=interval)
                self._last_chem_step = i
            for diag in self.diagnostics:
                diag.record(self, i, i * self.numerics.dt)
            if reporter is not None and reporter.due(i):
                reporter.emit(i)

        written: list[Path] = []
        if self.diagnostics:
            from diagnostics import DiagnosticsArchive

            archive_path = out / "diagnostics.h5"
            with DiagnosticsArchive(archive_path) as archive:
                for diag in self.diagnostics:
                    diag.export(self, archive)
            written.append(archive_path)
            for plot in self.plots:
                written.extend(plot.render(archive_path, out))
        elif self.plots:
            raise ValueError("plots require at least one diagnostic collector")
        if verbose:
            e = self.energetics
            e_out = e.output_energy_j()
            stim_frac = e.output_energy_stim_j() / e_out if e_out > 0.0 else 0.0
            print(f"  energy ledger: E_init={e.initial_energy_j:.3e} J, "
                  f"drive ∫P_in dt={e.input_energy_j():.3e} J, "
                  f"E_out={e_out:.3e} J (stimulated fraction {stim_frac:.3f})")
            print(f"  intrinsic efficiency = {e.intrinsic_efficiency():.3e}")
            print(f"  peak cavity photons = {self.gain.total_photons():.3e}, net gain g = "
                  f"{self._gain_coeff:.3e} 1/m")
            print(f"  wrote {len(written)} outputs to {out}/")
        if reporter is not None:
            # After the archive/plot writes, so a "done" heartbeat means outputs exist.
            reporter.finalize()
        return out
