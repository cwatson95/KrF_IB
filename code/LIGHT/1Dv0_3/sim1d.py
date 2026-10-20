"""The 1-D model: LIGHT v0_2's chemistry across a gap, on Fluid1D, on a mesh
graded toward the walls.

1Dv0_3 is 1Dv0_2 with OSIRIS's grid scaling (mesh.py). OSIRIS divides its
grid among processors by equal shares of a load; 1Dv0_3 places its cells
by equal shares of a load that is dense at the electrodes, where the
sheaths are, so the cells are small there and grow by a fixed ratio into
the uniform interior. ``mesh="uniform"`` is 1Dv0_2's mesh, and there
1Dv0_3 runs 1Dv0_2's model to the bit, cell by cell.

The input is v0_2's: a given power into the electron gas for ``pulse_s``,
deposited directly, uniform across the gap unless a profile is given (no
drive, beam or cavity model; LIGHT's rule). The gas fills a gap between
two grounded electrodes. Each cell runs v0_2's rate equations (its Kr2*
route included, and its bands at their defaults unless ``rates`` sets
them), electron energy and gas temperature (fluid1d.CellChemistry), and
transport.py moves the charged species and the electron energy between
the cells and into the walls under the space-charge field it solves, with
Fluid1D's scheme on each cell's own width.

One step of dt (Fluid1D's source-coupled split, Simulation._source_exchange):

    chemistry + electron energy, dt/2   (this module's stage: every row of
                                         the network, per cell, v0_2's
                                         rates refreshed every substep)
    transport, dt                       (transport.py: field, drift-
                                         diffusion, electron energy
                                         fluxes, walls, on the mesh)
    the gas, dt                         (neutrals.py, neutrals='move':
                                         flow, heat conduction into the
                                         electrodes, F2 and F diffusion,
                                         every species riding with it)
    chemistry + electron energy, dt/2

The stage's energy source per cell is the given power minus v0_2's
elastic loss n_e N_Kr kel(Te)(1 - Tg/Te), which heats the gas. The
electric field's work (field_work='conserve', the default) goes to the
electrons in the transport step, either sign, and to the ions, whose share
heats the gas, and the two are kept to the field's energy, which the books
carry. With field_work='fluid1d' the stage also takes the field's work on
the electrons where it is positive, as input (1Dv0_2's). Heat from the
reactions, the captured electrons' energy and the elastic loss raise each
cell's Tg at constant volume within the stage, as in v0_2 (U_th = sum n_s
cv_s k Tg), and the gas step then lets the cells expand and conduct.
Every total is a sum over the cells' own volumes.

What Fluid1D does not take from LIGHT: the ion mobilities (its cited
table: Kr+, Kr2+, F- in Kr), the electron mobility and diffusion (a
Maxwellian at the mean energy, from 1Dv0_3's own cross sections), and the
gas temperature its transport uses (held at tg0_K; it enters the ions'
diffusion and wall speed only). neutrals='move' moves the neutrals
(neutrals.py), and 'still' holds them in place, as 1Dv0_2 does (1 us at
1 atm spreads them by microns, which is safe on 39-um cells and not on
cells of a few microns).

    sim1d = importlib.import_module("1Dv0_3.sim1d")   # '1Dv0_3' is not an
    res = sim1d.run(sim1d.Params1D(mesh="graded", cells=140))   # identifier
    print(res.summary())
"""

from __future__ import annotations

import json
import time
import warnings
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import __version__
from . import fluid1d as bridge
from . import mesh as meshes
from . import neutrals as gasflow
from .model import KB_EV, QE, Model, Params

RUNS_DIR = Path(__file__).resolve().parent / "runs"
CHEMISTRY_FILES = ("model.py", "network.py", "species.py", "xsec.py",
                   "two_level.py", "data/Kr_Biagi71_LXCat.txt",
                   "data/F2_SIGLO_LXCat.txt", "data/kr_tail_ratio.txt")


@dataclass(frozen=True)
class Params1D:
    # --- the cell: v0_2's inputs (model.Params), same names and defaults
    pressure_atm: float = 1.0
    x_f2: float = 0.005
    tg0_K: float = 300.0
    power_W_cm3: float = 3.0e5
    pulse_s: float = 1.0e-6
    t_end_s: float = 2.0e-6
    ne_cm3: float = 1.0e15            # seed, as v0_2: Kr+ = n_e, uniform
    tau_sp_s: float = 1.0 / 1.43e8
    trimer_tg_pow: float = None
    # v0_2's bands (network.BANDS) as ((name, value), ...); one left out
    # keeps its default, and the unmeasured ones default to 0. For
    # sensitivity runs, never to match a device (v0_2's rule)
    rates: tuple = ()
    groups: tuple = ()
    eedf: str = "auto"                # 'two_t' (basic's) | 'maxwellian'
    capture: str = "auto"             # 'exact' (basic's) | 'three_halves'
    # --- the gap
    gap_cm: float = 1.0
    area_cm2: float = 1.0             # scales totals only
    cells: int = 128                  # with the wall cell below, the 64-cell
                                      # mesh from 2 um refined once (the
                                      # moving gas needs it, README)
    mesh: str = "graded"              # 'graded' (mesh.py: OSIRIS's equal-
                                      # load cuts) | 'uniform' (1Dv0_2's)
    wall_cell_um: float = 0.9508138901170884   # graded: the cell at each
                                      # electrode
    max_cell_um: float = None         # graded: a cap on any cell (None: none)
    neutrals: str = "move"            # 'move': the gas flows, conducts heat to
                                      # the electrodes (held at tg0_K) and its
                                      # F2 and F diffuse (neutrals.py);
                                      # 'still': every neutral held in place,
                                      # as in 1Dv0_2
    field_work: str = "conserve"      # 'conserve': the field's work on the
                                      # electrons and on the ions, both signs,
                                      # kept to the field's energy, the ions'
                                      # to the gas (transport.py); 'fluid1d':
                                      # Fluid1D's, its positive part given as
                                      # input, as in 1Dv0_2
    dt_s: float = 50e-12              # transport step (Fluid1D's dt)
    power_profile: tuple = None       # relative power per cell, volume-
                                      # weighted mean 1; None = uniform
    # --- output
    n_save: int = 801                 # time records
    n_profiles: int = 21              # profile snapshots
    name: str = "1Dv0_3"
    output_dir: str = None            # default 1Dv0_3/runs/
    threads: int = None               # numba threads (None: Fluid1D's pick)
    verbose: bool = True

    def replace(self, **kw) -> "Params1D":
        return replace(self, **kw)

    def cell_params(self) -> Params:
        """The 0-D Params of one cell: the chemistry 1Dv0_3 compiles, and the
        run v0_2 makes of the same inputs (the 0-D limit)."""
        return Params(pressure_atm=self.pressure_atm, x_f2=self.x_f2,
                      tg0_K=self.tg0_K, power_W_cm3=self.power_W_cm3,
                      pulse_s=self.pulse_s, t_end_s=self.t_end_s,
                      ne_cm3=self.ne_cm3, ne_mode="free",
                      tau_sp_s=self.tau_sp_s,
                      trimer_tg_pow=self.trimer_tg_pow,
                      rates=tuple((str(k), float(v)) for k, v in self.rates),
                      network="basic", groups=tuple(self.groups),
                      eedf=self.eedf, capture=self.capture)

    def make_mesh(self, length_m: float, area_m2: float) -> meshes.Mesh:
        """The mesh over a gap of length_m [m] (Fluid1D's own SI length)."""
        mx = None if self.max_cell_um is None else self.max_cell_um * 1e-6
        return meshes.build(self.mesh, self.cells, length_m, area_m2,
                            wall_cell_m=self.wall_cell_um * 1e-6,
                            max_cell_m=mx)

    def profile(self, volumes=None) -> np.ndarray:
        """The relative power per cell. Its mean is 1 over the gap's volume,
        so a profile moves power without changing the total."""
        if self.power_profile is None:
            return np.ones(self.cells)
        f = np.asarray(self.power_profile, dtype=float)
        if f.shape != (self.cells,) or np.any(f < 0) or f.mean() <= 0:
            raise ValueError("power_profile needs one non-negative value "
                             "per cell")
        if volumes is None or self.mesh == "uniform":
            return f / f.mean()
        v = np.asarray(volumes, dtype=float)
        return f / (np.sum(f * v) / np.sum(v))


def chemistry_md5() -> dict:
    here = Path(__file__).resolve().parent
    return {f: bridge.md5_of(here / f) for f in CHEMISTRY_FILES}


# ---------------------------------------------------------------------------
# The simulation (defined once Fluid1D is importable)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def simulation_class():
    """LightSimulation, a subclass of Fluid1D's Simulation."""
    bridge.load()
    import Fluid1D.simulation as fsim
    from Fluid1D import (BoltzmannSpec, Cavity, CircuitSpec, GasMixture,
                         Grid1D, Plasma, Preionization, Simulation,
                         build_lut)
    from Fluid1D.eedf_solver import MaxwellianEedfSolver
    from particles.constants import E_CHARGE as EV_TO_J
    from reactions import load_ion_mobility
    from reactions.network import OrphanChargeWarning

    from . import transport

    class LightSimulation(Simulation):
        """Fluid1D's Simulation with LIGHT's chemistry in the stage and
        1Dv0_3's mesh under the transport.

        Overrides (tests hold the parent methods they rely on to a
        recorded hash): the species registry, the mesh (after Fluid1D has
        built its uniform one, whose length it keeps), the electron
        transport table (Fluid1D's Maxwellian solver on 1Dv0_3's cross
        sections), the plasma model (transport.MeshPlasmaModel: grounded
        walls, real transport on the mesh), the field solver, the stage
        (_advance_electron_source), its energy source
        (_deliverable_heating_j_m3), the heavy-species chemistry call,
        which does nothing here (every row runs in the stage), and the
        three books Fluid1D keeps with one cell volume
        (_book_transport_closure, _resync_transport_energy_after_optics,
        _advance_optics), which here sum over the cells' own."""

        def __init__(self, p: Params1D):
            if p.neutrals not in ("move", "still"):
                raise ValueError("neutrals is 'move' or 'still'")
            if p.field_work not in transport.MeshPlasmaModel.FIELD_WORK:
                raise ValueError("field_work is 'conserve' or 'fluid1d'")
            self.p1 = p
            self.cell_model = Model(p.cell_params())
            self.chem = bridge.CellChemistry(self.cell_model)
            root = Path(p.output_dir) if p.output_dir else RUNS_DIR
            self.inputs_dir = root / p.name / "fluid1d_inputs"
            xs_dir, toml = bridge.write_inputs(self.inputs_dir,
                                               self.cell_model)
            fractions = {"Kr": 1.0 - p.x_f2}
            if p.x_f2 > 0.0:
                fractions["F2"] = p.x_f2
            gas = GasMixture(fractions, pressure_atm=p.pressure_atm,
                             temperature_k=p.tg0_K)
            eps0 = 1.5 * KB_EV * p.tg0_K
            plasma_model = transport.MeshPlasmaModel(
                reactions=None,
                ion_mobility=load_ion_mobility(bridge.ion_mobility_path()),
                mixture=gas.normalized_fractions, gas_temp_k=p.tg0_K,
                area_m2=p.area_cm2 * 1e-4, gap_m=p.gap_cm * 1e-2,
                secondary_emission_gamma=0.0,
                energy_floor_ev=bridge.EPS_LO,
                energy_ceiling_ev=bridge.EPS_HI,
                field_heating_cap_factor=None,
                electron_energy_closure="source_coupled",
                field_work=p.field_work)
            super().__init__(
                name=p.name,
                grid=Grid1D(cells=p.cells, length=p.gap_cm,
                            area_m2=p.area_cm2),
                gas=gas,
                plasma=Plasma(density=p.ne_cm3, ion="Kr+",
                              electron_mean_energy_ev=eps0,
                              secondary_emission_gamma=0.0),
                preionization=Preionization(mode="uniform",
                                            density=p.ne_cm3),
                cavity=Cavity(),
                reactions=toml,
                boltzmann=BoltzmannSpec(backend="eedf",
                                        key_kind="mean_energy_eV"),
                circuit=CircuitSpec(backend="stub", phi_left=0.0,
                                    phi_right=0.0),
                volume_units="cm",
                plasma_model=plasma_model,
                output_dir=root,
                cross_sections_dir=xs_dir,
                ion_mobility_path=bridge.ion_mobility_path(),
            )

        # ------------------------------------------------------- build
        def _build(self, numerics):
            import backend
            if self.p1.threads:
                backend.configure(threads=self.p1.threads)
            super()._build(numerics)
            bridge.stage_kernel()              # define the stage (numba)
            self._start_ledger()

        def _build_fluids(self):
            original = fsim.excimer_system_registry
            fsim.excimer_system_registry = (
                lambda: bridge.registry(self.chem.heavy))
            try:
                super()._build_fluids()
            finally:
                fsim.excimer_system_registry = original
            # v0_2 seeds Kr+ out of the krypton
            kr_plus = self.charged.density[:, self.charged.column["Kr+"]]
            self.fluid.density[:, self.fluid.column["Kr"]] -= kr_plus
            # 1Dv0_3's mesh over Fluid1D's own gap (its SI length to the
            # bit); the uniform one is Fluid1D's Mesh1D exactly
            fm = self.mesh
            m = self.p1.make_mesh(fm.length, fm.area_m2)
            if m.uniform and not (np.array_equal(m.widths,
                                                 np.full(m.n_cells, fm.dx))
                                  and np.array_equal(m.centres,
                                                     fm.cell_centers)):
                raise RuntimeError("the uniform mesh is not Fluid1D's")
            self.mesh = m
            self.vols = m.volumes
            self._plasma_model.set_mesh(m)
            if not m.uniform:
                # nothing may read one cell size from here on: Fluid1D's
                # paths that would are overridden, and any other fails
                vols = self.vols
                self.dx = self.cell_volume = None
                for fl in (self.fluid, self.charged):
                    fl.dx = fl.cell_volume = None
                    fl.total_number = (
                        lambda name, fl=fl: float(np.sum(
                            fl.density[:, fl.column[name]] * vols)))

        def _initialize_runtime_state(self):
            super()._initialize_runtime_state()
            self.poisson = transport.MeshPoisson(self.mesh)
            self.cell_centers = self.mesh.cell_centers

        def _build_chemistry(self):
            with warnings.catch_warnings():
                # the file holds the electron-impact rows alone, so F- and
                # Kr+ look unconsumed to Fluid1D; the stage consumes them
                warnings.simplefilter("ignore", OrphanChargeWarning)
                super()._build_chemistry()
            if tuple(self.plasma_chem.species_names) != self.chem.names:
                raise RuntimeError("Fluid1D's species order is not "
                                   "CellChemistry's")

        def _build_coeff_table(self):
            solver = MaxwellianEedfSolver(
                self.cross_sections_dir,
                mean_energy_range_ev=(bridge.EPS_LO, bridge.EPS_HI),
                sigma_lut_points=4096)
            return build_lut(self.gas.normalized_fractions, solver=solver,
                             en_grid_td=np.geomspace(0.1, 1000.0, 400),
                             gas_temp_k=self.gas.temperature_k,
                             reactions=self.compiled, cache_dir=None,
                             key_kind="mean_energy_eV")

        # ------------------------------------------- the books, per cell
        def _total(self, per_m3) -> float:
            """A per-cell density of something [x m^-3] summed over the
            cells' volumes [x]."""
            return float(np.sum(np.asarray(per_m3) * self.vols))

        def _book_transport_closure(self, report):
            """Fluid1D's, with the totals over the cells' own volumes."""
            closure = report.energy_closure
            if closure is None:
                return
            book = self._book_closure_energy
            book("field_work_j", self._total(closure.field_work_j_m3))
            book("elastic_transfer_j",
                 self._total(closure.elastic_transfer_j_m3))
            book("continuous_inelastic_sink_j",
                 self._total(closure.continuous_inelastic_sink_j_m3))
            book("coefficient_bound_adjustment_j",
                 self._total(closure.coefficient_bound_adjustment_j_m3))
            self._continuous_prepaid_j_m3 += \
                closure.continuous_inelastic_sink_j_m3
            self._last_interval_field_work_j_m3 = \
                closure.field_work_j_m3.copy()
            self._last_interval_elastic_j_m3 = \
                closure.elastic_transfer_j_m3.copy()
            self._last_interval_bound_j_m3 = \
                closure.coefficient_bound_adjustment_j_m3.copy()
            self._last_interval_h_s = float(report.dt)
            live = self.last_energy_closure_report
            live.field_work_j_m3 += closure.field_work_j_m3
            live.elastic_transfer_j_m3 += closure.elastic_transfer_j_m3
            live.continuous_inelastic_sink_j_m3 += \
                closure.continuous_inelastic_sink_j_m3
            live.coefficient_bound_adjustment_j_m3 += \
                closure.coefficient_bound_adjustment_j_m3
            live.maximum_energy_jump_ev = max(live.maximum_energy_jump_ev,
                                              closure.maximum_energy_jump_ev)
            live.limiter_activity += closure.limiter_activity
            live.max_energy_residual_j_m3 = max(
                live.max_energy_residual_j_m3,
                closure.max_energy_residual_j_m3)
            live.max_ampere_residual = max(live.max_ampere_residual,
                                           closure.max_ampere_residual)
            step = self._energy_closure_step
            step["maximum_source_step_energy_jump_ev"] = max(
                step["maximum_source_step_energy_jump_ev"],
                closure.maximum_energy_jump_ev)
            step["limiter_activity"] += closure.limiter_activity
            step["max_energy_closure_residual_j_m3"] = max(
                step["max_energy_closure_residual_j_m3"],
                closure.max_energy_residual_j_m3)
            step["max_ampere_residual"] = max(step["max_ampere_residual"],
                                              closure.max_ampere_residual)

        def _resync_transport_energy_after_optics(self):
            """Fluid1D's, with the clipped energy over the cells' volumes."""
            n_e = np.maximum(
                self.charged.density[:, self.charged.column["e-"]], 0.0)
            eps = np.clip(self.charged.electron_mean_energy_ev,
                          *self.plasma._eps_bounds)
            clipped_j = self._total(
                (eps - self.charged.electron_mean_energy_ev) * n_e) * EV_TO_J
            if clipped_j != 0.0:
                self._book_closure_energy("coefficient_bound_adjustment_j",
                                          clipped_j)
            self.charged.electron_mean_energy_ev[:] = eps
            self.plasma._last_n_eps = n_e * eps * EV_TO_J

        # ------------------------------------------------------ ledger
        def _start_ledger(self):
            p, c = self.p1, self.grid_spec.cells
            self.tg = np.full(c, float(p.tg0_K))
            self._power = p.power_W_cm3 * 1e6 * p.profile(self.vols)  # W/m^3
            z = np.zeros
            self.led = {          # per cell, eV m^-3, cumulative
                "input": z(c), "field_work": z(c), "elastic": z(c),
                "heat_reactions": z(c), "adjust": z(c), "to_e": z(c),
                "delivered": z(c),
                **{f"rad:{b}": z(c) for b in self.chem.bands},
                "krf_formed": z(c),                   # m^-3, not energy
                "wall_heat": z(c),        # conducted into the electrodes,
                                          # booked in the two wall cells
                "field_e": z(c),          # the field's work on the electrons
                "field_i": z(c)}          # and on the ions, to the gas
                                          # (field_work='conserve')
            self._conserve = p.field_work == "conserve"
            self.field_defect = 0.0       # eV: the theorem's residual
            self.field_moved = 0.0        # eV: the field's energy the gas
                                          # changed, carrying the charge
            self._positivity_limits = 0
            self.u = z(c)                 # the gas's velocity [m/s]
            self.flow = None
            self._gas_substeps = []
            if p.neutrals == "move":
                chem = self.chem
                masses = np.array([bridge.xsec.ME] + [
                    gasflow.mass_kg(bridge.SPECIES[s].atoms)
                    for s in chem.heavy])
                heavy = np.array([n != "e-" for n in chem.names])
                self.flow = gasflow.NeutralGas(
                    self.mesh, chem.names, masses, chem.cv, heavy,
                    e_col=chem.e_col,
                    transport=gasflow.BathTransport(chem.names, masses),
                    t_wall=float(p.tg0_K))
                self._m_heavy = np.where(heavy, masses, 0.0)
            self.e_other = 0.0            # eV: transport + optics, whole gap
            self.field_dropped = 0.0      # eV: negative field work not given
            self._t_interval = (0.0, 0.0)
            self._half_start = 0.0
            self._substeps = []
            self._capped = []     # (t, species, energy, source frequency)
            self._transport_substeps = []
            self.e_total0 = self.energy_total()
            steps = self.numerics.steps
            self._rec_steps = set(np.unique(np.round(np.linspace(
                0, steps, max(p.n_save, 2))).astype(int)).tolist())
            self._prof_steps = set(np.unique(np.round(np.linspace(
                0, steps, max(p.n_profiles, 2))).astype(int)).tolist())
            self.records, self.profiles = [], []
            self._wall0 = time.time()
            self._next_print = 0.0
            self._record(0)

        def energy_total(self) -> float:
            """Chemical + electron + gas thermal energy in the gap [eV],
            the gas's kinetic energy where it moves, and the field's energy
            where its work is kept to it (field_work='conserve')."""
            dens = self._combined_densities()
            u = (self.chem.chemical_energy(dens)
                 + self.plasma._last_n_eps / QE
                 + self.chem.heat_capacity(dens) * self.tg)
            if self.flow is not None:
                rho = (dens * self._m_heavy).sum(axis=1)
                u = u + 0.5 * rho * self.u * self.u / QE
            total = self._total(u)
            if self._conserve:
                total += self._field_energy_ev()
            return total

        def _field_energy_ev(self) -> float:
            """The field's energy in the gap [eV] of the charge as it
            stands (the exact field between the grounded electrodes)."""
            return (self.plasma.field_energy_j_m2(self.charged)
                    * self.mesh.area_m2 / QE)

        def _heat_gas(self, q_ev_m3, dens_before, dens_after):
            """Raise Tg by q [eV m^-3] at constant volume while the
            densities change from before to after (v0_2's rule: the
            particles a reaction makes take their cv k Tg from the heat)."""
            u = self.chem.heat_capacity(dens_before) * self.tg + q_ev_m3
            self.tg = u / self.chem.heat_capacity(dens_after)

        # --------------------------------------------------- the split
        def _source_exchange(self, t, interval):
            self._t_interval = (t - interval, t)
            self._half_start = t - interval
            super()._source_exchange(t, interval)

        def _exchange(self, t, h):
            e0 = self.energy_total()
            super()._exchange(t, h)
            self._transport_substeps.append(self.plasma._last_substeps)
            booked = self._take_field_work() if self._conserve else 0.0
            self.e_other += self.energy_total() - e0 - booked
            if self.flow is not None:
                self._move_gas(h)

        def _take_field_work(self) -> float:
            """field_work='conserve': the advance applied the electrons'
            share of the field's work to their energy; the ions' share
            heats the gas here. Both are booked, the kernel's mean-energy
            bounds go with the bounds (returned [eV], so the exchange's
            other energy leaves them out), and the theorem's residual,
            round-off, is kept."""
            pl = self.plasma
            wi = pl._last_ion_work_j_m3 / QE                     # eV m^-3
            dens = self._combined_densities()
            self._heat_gas(wi, dens, dens)
            self.led["field_i"] += wi
            self.led["field_e"] += self._last_interval_field_work_j_m3 / QE
            bound = self._last_interval_bound_j_m3 / QE
            self.led["adjust"] += bound
            w0, w1 = pl._last_field_energy_j_m2
            self.field_defect += ((w1 - w0) + pl._last_field_work_j_m2) \
                * self.mesh.area_m2 / QE
            self._positivity_limits += pl._last_positivity_limits
            return self._total(bound)

        def _move_gas(self, h):
            """The neutral gas over h (neutrals.py): every species and the
            electrons' mean energy ride with it, F2 and F diffuse, heat
            flows into the electrodes. Its energy is conserved but for that
            heat, which the wall cells book, and the electrons' mean energy
            is held in the transport table's range (what that moves is
            booked with the bounds)."""
            chem = self.chem
            w_field = self._field_energy_ev() if self._conserve else 0.0
            n_all = self._combined_densities()
            # the electrons' energy as the books hold it, Fluid1D's energy
            # density, per electron (its mean energy is blended toward the
            # local field's where the electrons are few, and the transport
            # step books that difference itself)
            ne0 = n_all[:, chem.e_col]
            eps = np.array(self.charged.electron_mean_energy_ev, float)
            live = ne0 > 0.0
            eps[live] = self.plasma._last_n_eps[live] / (ne0[live] * QE)
            r = self.flow.step(n_all, self.tg, self.u, eps, h)
            self._gas_substeps.append(r.substeps)
            self._adopt_combined_densities(r.dens)
            self.tg, self.u = r.T, r.u
            ne = np.maximum(r.dens[:, chem.e_col], 0.0)
            carried = ne * r.eps                          # eV/m^3, conserved
            lo, hi = self.plasma._eps_bounds
            self.plasma.adopt_electron_energy_density(
                np.clip(carried, lo * ne, hi * ne) * QE,
                charged_fluid=self.charged)
            # what Fluid1D's adoption changed (the table's bounds, an empty
            # cell) is booked with the bounds
            self.led["adjust"] += self.plasma._last_n_eps / QE - carried
            w = self.mesh.widths
            self.led["wall_heat"][0] += r.wall_heat[0] / (w[0] * QE)
            self.led["wall_heat"][-1] += r.wall_heat[1] / (w[-1] * QE)
            if self._conserve:
                # the charge the gas carried changed the field's energy,
                # booked with the other energy the transport moves, and
                # kept, so the field's own books close with the gas moving
                moved = self._field_energy_ev() - w_field
                self.e_other += moved
                self.field_moved += moved

        def _advance_optics(self, t, interval):
            """Fluid1D's, with the photoelectrons' energy over the gap's
            volume, and the energy it moves booked as e_other."""
            e0 = self.energy_total()
            p_out, self._gain_coeff = self.gain.advance(
                interval, fluid=self.fluid, charged=self.charged)
            self._book_closure_energy(
                "photoelectron_birth_energy_j",
                self.gain.last_photoelectron_birth_energy_j
                * float(np.sum(self.vols)))
            if self._real_transport:
                self._resync_transport_energy_after_optics()
            p_in = self._p_in_accum / interval if interval > 0 else 0.0
            self.energetics.record(t, 0.0, p_out, p_in, interval_s=interval)
            self._p_in_accum = 0.0
            self.e_other += self.energy_total() - e0

        def _advance_source_coupled_interval(self, n_all, interval):
            self._last_heavy_nu = 0.0          # every row ran in the stage

        def _deliverable_heating_j_m3(self):
            """The interval's field and elastic energy for the stage
            [J m^-3]: the field's positive work on the electrons minus
            v0_2's elastic loss (whose energy heats the gas now). Fluid1D
            splits it between the two half-stages (a predictor from the
            last interval, then the remainder). The given power is not in
            it: each half-stage takes the power over its own time window
            (_advance_electron_source), so the pulse's end needs no
            correction."""
            if self._conserve:
                # the field's work went to the electrons in the advance and
                # its ions' share to the gas: v0_2's elastic loss alone
                h = self._last_interval_h_s
                dens = self._combined_densities()
                te = np.clip(self.charged.electron_mean_energy_ev,
                             *self.plasma._eps_bounds) / 1.5
                el = self.chem.elastic_power(dens, te, self.tg) * h
                self._heat_gas(el, dens, dens)
                self.led["elastic"] += el
                return -el * QE
            field = self._last_interval_field_work_j_m3
            pos = np.maximum(field, 0.0)
            dropped = self._total(field - pos)
            self._book_closure_energy("undelivered_field_work_j", dropped)
            self.field_dropped += dropped / QE
            h = self._last_interval_h_s
            dens = self._combined_densities()
            te = np.clip(self.charged.electron_mean_energy_ev,
                         *self.plasma._eps_bounds) / 1.5
            el = self.chem.elastic_power(dens, te, self.tg) * h  # eV/m^3
            self._heat_gas(el, dens, dens)
            self.led["field_work"] += pos / QE
            self.led["elastic"] += el
            return pos - el * QE

        def _power_ev(self, interval):
            """The given power over the next half-stage [eV m^-3 s^-1],
            averaged over its window (which may hold the pulse's end)."""
            t0 = self._half_start
            t1 = t0 + interval
            self._half_start = t1
            on = max(0.0, min(t1, self.p1.pulse_s) - max(t0, 0.0))
            p_in = self._power * on / QE                        # eV/m^3
            self.led["input"] += p_in
            return p_in / interval

        def _advance_electron_source(self, interval,
                                     power_source_j_m3_s=None):
            """Every row of the network, in every cell, over ``interval``,
            with the electron energy (the stage)."""
            if interval <= 0.0:
                return
            chem = self.chem
            lo, hi = self.plasma._eps_bounds
            n_all = self._combined_densities()
            eps = np.clip(self.charged.electron_mean_energy_ev, lo, hi)
            source = self._power_ev(interval)
            if power_source_j_m3_s is not None:
                source = source + np.asarray(power_source_j_m3_s, float) / QE
            st = chem.stage(n_all, self.plasma._last_n_eps / QE, eps,
                            self.tg, source, interval, lo, hi)
            self._substeps.append(st.n_sub)
            if st.n_sub == bridge.MAX_SUBSTEPS:
                self._capped.append((self._t_interval[1],) + st.freqs)
            self._heat_gas(st.heat, n_all, st.dens)
            self._adopt_combined_densities(st.dens)
            self.plasma.adopt_electron_energy_density(
                st.energy * QE, charged_fluid=self.charged)
            led = self.led
            led["heat_reactions"] += st.heat
            led["adjust"] += st.adjust
            led["to_e"] += st.counters["to_e"]
            led["delivered"] += source * interval
            led["krf_formed"] += st.counters["krf_formed"]
            for b in chem.bands:
                led[f"rad:{b}"] += st.counters[f"rad:{b}"]
            self._last_nu_max = st.nu_e

        # ------------------------------------------------------ records
        def _step(self, i):
            super()._step(i)
            if i in self._rec_steps or i in self._prof_steps:
                self._record(i)
            if self.p1.verbose:
                frac = i / self.numerics.steps
                if frac >= self._next_print:
                    self._print(i)
                    self._next_print += 0.1

        def _record(self, i):
            t = i * self.numerics.dt
            chem = self.chem
            dens = self._combined_densities()
            eps = self.charged.electron_mean_energy_ev
            ne = dens[:, chem.e_col]
            mid = dens.shape[0] // 2
            vols = self.vols
            vol = float(np.sum(vols))
            if i in self._rec_steps:
                rad = sum(self._total(self.led[f"rad:{b}"])
                          for b in chem.bands)
                e_in = self._total(self.led["input"]) \
                    + self._total(self.led["field_work"])
                adjust = self._total(self.led["adjust"])
                e_tot = self.energy_total()
                w = vols / vol
                self.records.append({
                    "t": t,
                    "mean": ((w @ dens) * 1e-6).tolist(),
                    "centre": (dens[mid] * 1e-6).tolist(),
                    "te_mean": float(np.sum(ne * eps * vols)
                                     / max(np.sum(ne * vols), 1e-300) / 1.5),
                    "te_centre": float(eps[mid] / 1.5),
                    "tg_mean": float(np.sum(self.tg * w)),
                    "tg_centre": float(self.tg[mid]),
                    "ledger": {k: self._total(v) / vol * 1e-6
                               for k, v in self.led.items()},
                    "ledger_centre": {k: float(v[mid]) * 1e-6
                                      for k, v in self.led.items()},
                    "e_other": self.e_other / vol * 1e-6,
                    "e_total": (e_tot - self.e_total0) / vol * 1e-6,
                    "closure": ((e_tot - self.e_total0)
                                - (e_in - rad + adjust + self.e_other
                                   - self._total(self.led["wall_heat"])))
                    / vol * 1e-6,
                })
            if i in self._prof_steps:
                prof = {"t": t, "dens": (dens * 1e-6).T.tolist(),
                        "te": (eps / 1.5).tolist(), "tg": self.tg.tolist()}
                if self.flow is not None:
                    heavy = self._m_heavy > 0.0
                    prof["u"] = self.u.tolist()
                    prof["p"] = (dens[:, heavy].sum(axis=1) * gasflow.KB
                                 * self.tg).tolist()
                self.profiles.append(prof)

        def _print(self, i):
            r = self.records[-1] if self.records else None
            if r is None:
                return
            e_in = r["ledger"]["input"]
            key = "rad:248 nm (KrF B-X)"
            share = (100 * r["ledger"].get(key, 0.0) / e_in
                     if e_in > 0 else 0.0)
            ne = r["mean"][self.chem.e_col]
            print(f"  t={r['t'] * 1e9:7.1f} ns  <n_e>={ne:.2e}  "
                  f"Te={r['te_mean']:.2f} eV  Tg={r['tg_mean']:.0f} K  "
                  f"248 nm {share:.2f} %  "
                  f"wall {time.time() - self._wall0:.0f} s", flush=True)

    return LightSimulation


# ---------------------------------------------------------------------------
# Running and the result
# ---------------------------------------------------------------------------
def run(params: Params1D = None, **kw) -> "Result1D":
    """run(cells=64) or run(Params1D(...)): build, run, save
    <output_dir>/<name>/result.json, return the Result1D."""
    p = params if params is not None else Params1D()
    if kw:
        p = p.replace(**kw)
    sim_cls = simulation_class()
    from Fluid1D import Numerics
    steps = int(round(p.t_end_s / p.dt_s))
    numerics = Numerics(steps=steps, dt=p.dt_s,
                        fluid_every=max(1, int(round(1e-9 / p.dt_s))),
                        fluid_substeps=1, field_heating_cap_factor=None,
                        electron_energy_closure="source_coupled")
    t0 = time.time()
    sim = sim_cls(p)
    if p.verbose:
        print(f"LIGHT 1Dv0_3: {p.cells} {p.mesh} cells over {p.gap_cm:g} cm, "
              f"{steps} steps of {p.dt_s * 1e12:g} ps", flush=True)
    out = sim.run(numerics, verbose=False)
    res = Result1D.from_simulation(sim, time.time() - t0)
    res.save(Path(out) / "result.json")
    return res


class Result1D:
    """A 1Dv0_3 run: time records, profiles and the energy ledger.

    Densities in cm^-3, energies in eV per cm^3 of gap (averaged over the
    gap's volume, or the centre cell's own), temperatures in eV (Te) and K
    (Tg). ``x_cm`` holds the cell centres and ``widths_cm`` their widths."""

    def __init__(self, data: dict):
        self.data = data
        self.p = data["params"]
        self.species = data["species"]
        self.t = np.array([r["t"] for r in data["records"]])
        self.x_cm = np.array(data["x_cm"])
        self.widths_cm = np.array(data["widths_cm"])

    @classmethod
    def from_simulation(cls, sim, wall_s):
        chem = sim.chem
        vol = float(np.sum(sim.vols))
        data = {
            "model": f"LIGHT 1Dv0_3 {__version__} (LIGHT v0_2 chemistry on "
                     "Fluid1D's 1-D transport, on 1Dv0_3's mesh)",
            "params": {k: (list(v) if isinstance(v, tuple) else v)
                       for k, v in asdict(sim.p1).items()},
            "species": list(chem.names),
            "bands": list(chem.bands),
            "mesh": sim.mesh.describe(),
            "x_cm": (sim.mesh.centres * 100).tolist(),
            "widths_cm": (sim.mesh.widths * 100).tolist(),
            "records": sim.records,
            "profiles": sim.profiles,
            "ledger_cells": {k: (v * 1e-6).tolist()
                             for k, v in sim.led.items()},
            "field_dropped": sim.field_dropped / vol * 1e-6,
            # field_work='conserve': the field's energy at the end, the
            # theorem's residual summed over the advances (round-off), what
            # the gas's carrying the charge changed in the field's energy,
            # and the cells whose outflow the positive transport scaled
            "field_energy": (sim._field_energy_ev() / vol * 1e-6
                             if sim._conserve else 0.0),
            "field_defect": sim.field_defect / vol * 1e-6,
            "field_moved": sim.field_moved / vol * 1e-6,
            "positivity_limits": int(sim._positivity_limits),
            "substeps": {"mean": float(np.mean(sim._substeps)),
                         "max": int(np.max(sim._substeps))},
            "transport_substeps": {
                "mean": float(np.mean(sim._transport_substeps)),
                "max": int(np.max(sim._transport_substeps))},
            "gas_substeps": ({"mean": float(np.mean(sim._gas_substeps)),
                              "max": int(np.max(sim._gas_substeps))}
                             if sim._gas_substeps else None),
            "eps_bounds_eV": list(sim.plasma._eps_bounds),
            # the bands the network's rows read, and Kr* + 2Kr's rate (the
            # summary's bands line and warnings, as v0_2's report makes them)
            "bands_read": sorted({c.row.rate[1]
                                  for c in sim.cell_model.compiled
                                  if c.row.rate[0] == "param"}),
            "ak8_rate": float(sum(c.rate for c in sim.cell_model.compiled
                                  if c.row.rid == "AK8_Kr")),
            "wall_s": wall_s,
            "provenance": {"fluid1d_src": str(bridge.src_dir()),
                           "chemistry_md5": chemistry_md5()},
        }
        return cls(data)

    # --------------------------------------------------------- access
    def series(self, key, where="mean"):
        """Time series of a species density [cm^-3] ('mean' over the gap
        or 'centre')."""
        i = self.species.index(key)
        return np.array([r[where][i] for r in self.data["records"]])

    def ledger(self, key, where="ledger"):
        return np.array([r[where][key] for r in self.data["records"]])

    def fate(self, where="ledger") -> dict:
        """Percent of the energy put in (power + field work) by t_end:
        radiated by band, gas heat, still stored, lost to the walls. For
        the whole gap ('ledger') the entries sum to 100 %: the closure
        residual is the bookkeeping test (it recomputes the gap's energy
        from its densities, Te and Tg). A cell ('ledger_centre') exchanges
        energy with its neighbours, so only its own flows are listed."""
        r = self.data["records"][-1]
        led = r[where]
        e_in = led["input"] + led["field_work"]
        if e_in <= 0:
            return {}
        out = {f"radiated {b}": 100 * led[f"rad:{b}"] / e_in
               for b in self.data["bands"]}
        out["heat (elastic)"] = 100 * led["elastic"] / e_in
        out["heat (reactions + captured electrons)"] = (
            100 * led["heat_reactions"] / e_in)
        conserve = self.p.get("field_work") == "conserve"
        f_i = led.get("field_i", 0.0) if conserve else 0.0
        if conserve:
            # the ions' share of the field's work heats the gas
            out["heat (ions in the field)"] = 100 * f_i / e_in
        if where == "ledger":
            # e_total is the change of chemical + electron + gas thermal
            # (and kinetic, and the field's) energy; the gas heat above is
            # what reached the gas, of which the electrodes took wall_heat
            # by conduction
            stored = (r["e_total"] - led["elastic"] - led["heat_reactions"]
                      + led.get("wall_heat", 0.0))
            if conserve:
                stored -= f_i
            out["stored (species + electrons)"] = 100 * stored / e_in
            out["lost to the walls"] = -100 * r["e_other"] / e_in
            # positive: energy the mean-energy bounds took out (the
            # drained electrons beside the walls pinned at the ceiling)
            out["cut by the mean-energy bounds"] = -100 * led["adjust"] / e_in
            out["closure residual"] = -100 * r["closure"] / e_in
        return out

    def share_248(self, where="ledger") -> float:
        led = self.data["records"][-1][where]
        if led["input"] <= 0:
            return float("nan")
        return 100 * led["rad:248 nm (KrF B-X)"] / led["input"]

    def summary(self) -> str:
        p, d = self.p, self.data
        last = d["records"][-1]
        lines = [d["model"],
                 f"  {p['pressure_atm']:g} atm, F2 {100 * p['x_f2']:g} %, "
                 f"{p['power_W_cm3'] / 1e3:g} kW/cm^3 for "
                 f"{p['pulse_s'] * 1e9:g} ns, end {p['t_end_s'] * 1e9:g} ns, "
                 f"seed n_e {p['ne_cm3']:.1e}",
                 f"  gap {p['gap_cm']:g} cm, {d['mesh']}",
                 f"  dt {p['dt_s'] * 1e12:g} ps, transport substeps mean "
                 f"{d['transport_substeps']['mean']:.1f} max "
                 f"{d['transport_substeps']['max']}, chemistry substeps "
                 f"mean {d['substeps']['mean']:.1f} max "
                 f"{d['substeps']['max']}, wall {d['wall_s']:.0f} s",
                 f"  deposited {last['ledger']['input'] * QE:.3f} J/cm^3 "
                 f"(gap mean)" + (
                     "" if p.get("field_work") == "conserve" else
                     f", field work "
                     f"{last['ledger']['field_work'] * QE:.2e} J/cm^3")]
        read = d.get("bands_read", [])
        if read:
            # v0_2's report header: the bands set, and the unmeasured ones
            # left at 0
            from .network import BANDS
            rates = [tuple(r) for r in p["rates"]]
            set_ = [f"{k} {v:.3g}" for k, v in rates]
            if p["trimer_tg_pow"] is not None:
                set_.append(f"trimer_tg_pow {p['trimer_tg_pow']:.3g}")
            given = {k for k, _ in rates}
            off = [k for k in read if BANDS[k].default == 0.0
                   and k not in given]
            lines.append("  bands: " + (", ".join(set_) if set_ else
                                        "every band at its default")
                         + (";  no measurement, 0: " + ", ".join(off)
                            if off else ""))
        lines.append("  percent of the energy in    gap    centre cell")
        gap, mid = self.fate("ledger"), self.fate("ledger_centre")
        for k, v in gap.items():
            if k.startswith("_"):
                continue
            m = mid.get(k)
            lines.append(f"    {k:40s} {v:7.3f}"
                         + (f"  {m:7.3f}" if m is not None else ""))
        e_in = last["ledger"]["input"]
        if e_in > 0 and p.get("field_work") == "conserve":
            led = last["ledger"]
            moved = (f"the gas carrying the charge changed it by "
                     f"{100 * d.get('field_moved', 0.0) / e_in:.2g} %, "
                     if p.get("neutrals") == "move" else "")
            lines.append(
                f"  the field's work, kept to its energy: on the electrons "
                f"{100 * led['field_e'] / e_in:.3g} %, on the ions "
                f"{100 * led['field_i'] / e_in:.3g} % of the input, its "
                f"energy at the end {100 * d['field_energy'] / e_in:.2g} %, "
                f"{moved}the books' residual "
                f"{100 * d['field_defect'] / e_in:.2g} %")
        elif e_in > 0:
            # Fluid1D's split gives the electrons the field's positive work
            # and never its negative part. On a uniform mesh the sheath is
            # one unresolved cell and the negative part a phantom that grows
            # as the mesh coarsens (1Dv0_2's README, "The wall cells")
            lines.append(
                f"  field work given {100 * last['ledger']['field_work'] / e_in:.3g} "
                f"% of the input, negative work kept out "
                f"{-100 * d['field_dropped'] / e_in:.3g} %")
            # the energy in above counts the field's work, which in the
            # grounded gap is the coupled step's error (README, "The
            # field's work"); the given power alone is the input
            lines.append(
                f"  248 nm of the given power alone: gap "
                f"{self.share_248():.3f} %, centre cell "
                f"{self.share_248('ledger_centre'):.3f} %")
        prof = d["profiles"][-1] if d["profiles"] else {}
        if "u" in prof and e_in > 0:
            u = np.array(prof["u"])
            pr = np.array(prof["p"])
            sub = d.get("gas_substeps") or {"mean": 0.0}
            lines.append(
                f"  gas: |u| up to {np.max(np.abs(u)):.3g} m/s at the end, "
                f"pressure {pr.min() / 1e3:.4g} to {pr.max() / 1e3:.4g} kPa, "
                f"heat conducted into the electrodes "
                f"{100 * last['ledger'].get('wall_heat', 0.0) / e_in:.3g} % "
                f"of the input ({sub['mean']:.2f} gas substeps per step)")
        elif "u" not in prof:
            lines.append("  gas: held still (neutrals='still')")
        lines += [f"  Tg at the end: centre {last['tg_centre']:.0f} K, "
                  f"gap mean {last['tg_mean']:.0f} K",
                  f"  F2 left: {100 * self.series('F2')[-1] / self.series('F2')[0]:.2f} % "
                  f"(gap mean)"]
        lines += [f"  WARNING {w}" for w in self.warnings()]
        return "\n".join(lines)

    def warnings(self) -> list:
        """v0_2's report warnings that concern the gas temperature, read in
        the centre cell (the interior), and the seed that did not break
        down."""
        d, out = self.data, []
        recs = d["records"][1:]
        if not recs:
            return out
        tg = max(r["tg_centre"] for r in recs)
        if tg > 800.0:
            out.append(f"Tg reached {tg:.0f} K in the centre cell. The "
                       "trimer's (Tg/300)^-3.32 follows Shui's calculation "
                       "for 300-800 K (measured 294-425 K) and is "
                       "extrapolated above it")
        if d.get("ak8_rate", 0.0) > 0.0 and tg > 451.0:
            out.append(f"Tg reached {tg:.0f} K in the centre cell. Kr* + "
                       "2Kr's (Tg/300)^-0.64 is measured over 300-451 K "
                       "(Wieme & Lenaerts 1981) and extrapolated above it")
        if tg > 1000.0:
            out.append("Tg above 1000 K: Kr2+ and Kr2F* redissociate on these "
                       "time scales, which the model does not include")
        te_top = 0.99 * d["eps_bounds_eV"][1] / 1.5
        if max(r["te_centre"] for r in recs) >= te_top and \
                self.share_248("ledger_centre") < 1e-3:
            out.append("the centre cell sat at the mean-energy ceiling and "
                       "made no KrF: the seed n_e was too low for the power "
                       "to be absorbed below the ceiling (README, the seed)")
        return out

    # ------------------------------------------------------------ files
    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.data) + "\n")
        return path

    @classmethod
    def load(cls, path) -> "Result1D":
        return cls(json.loads(Path(path).read_text()))
