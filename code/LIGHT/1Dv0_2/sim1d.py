"""The 1-D model: LIGHT v0_2's chemistry across a gap, on Fluid1D.

The input is v0_2's: a given power into the electron gas for ``pulse_s``,
deposited directly, uniform across the gap unless a profile is given (no
drive, beam or cavity model; LIGHT's rule). The gas fills a gap between
two grounded electrodes. Each cell runs v0_2's rate equations (its Kr2*
route included, and its bands at their defaults unless ``rates`` sets
them), electron energy and gas temperature (fluid1d.CellChemistry), and
Fluid1D moves the charged species and the electron energy between the
cells and into the walls under the space-charge field it solves.

One step of dt (Fluid1D's source-coupled split, Simulation._source_exchange):

    chemistry + electron energy, dt/2   (this module's stage: every row of
                                         the network, per cell, v0_2's
                                         rates refreshed every substep)
    transport, dt                       (Fluid1D: field, drift-diffusion,
                                         electron energy fluxes, walls)
    chemistry + electron energy, dt/2

The stage's energy source per cell is the given power, plus the electric
field's work on the electrons where it is positive, minus v0_2's elastic
loss n_e N_Kr kel(Te)(1 - Tg/Te), which heats the gas. Heat from the
reactions, the captured electrons' energy and the elastic loss raise each
cell's Tg at constant volume, as in v0_2 (U_th = sum n_s cv_s k Tg).

What Fluid1D does not take from LIGHT: the ion mobilities (its cited
table: Kr+, Kr2+, F- in Kr), the electron mobility and diffusion (a
Maxwellian at the mean energy, from 1Dv0_2's own cross sections), and the gas
temperature its transport uses (held at tg0_K; it enters the ions'
diffusion and wall speed only). The neutrals do not move (1 us at 1 atm
spreads them by microns).

    sim1d = importlib.import_module("1Dv0_2.sim1d")   # '1Dv0_2' is not an
    res = sim1d.run(sim1d.Params1D(cells=64))           # identifier
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
    cells: int = 128
    dt_s: float = 50e-12              # transport step (Fluid1D's dt)
    power_profile: tuple = None       # relative power per cell, mean 1;
                                      # None = uniform (v0_2's rule)
    # --- output
    n_save: int = 801                 # time records
    n_profiles: int = 21              # profile snapshots
    name: str = "1Dv0_2"
    output_dir: str = None            # default 1Dv0_2/runs/
    threads: int = None               # numba threads (None: Fluid1D's pick)
    verbose: bool = True

    def replace(self, **kw) -> "Params1D":
        return replace(self, **kw)

    def cell_params(self) -> Params:
        """The 0-D Params of one cell: the chemistry 1Dv0_2 compiles, and the
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

    def profile(self) -> np.ndarray:
        if self.power_profile is None:
            return np.ones(self.cells)
        f = np.asarray(self.power_profile, dtype=float)
        if f.shape != (self.cells,) or np.any(f < 0) or f.mean() <= 0:
            raise ValueError("power_profile needs one non-negative value "
                             "per cell")
        return f / f.mean()


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
    from Fluid1D import (BoltzmannSpec, Cavity, CircuitSpec,
                         DriftDiffusionPlasmaModel, GasMixture, Grid1D,
                         Plasma, Preionization, Simulation, build_lut)
    from Fluid1D.eedf_solver import MaxwellianEedfSolver
    from reactions import load_ion_mobility
    from reactions.network import OrphanChargeWarning

    class LightSimulation(Simulation):
        """Fluid1D's Simulation with LIGHT's chemistry in the stage.

        Overrides (tests hold the parent methods they rely on to a
        recorded hash): the species registry, the electron transport table
        (Fluid1D's Maxwellian solver on 1Dv0_2's cross sections), the plasma
        model (grounded walls, real transport), the stage
        (_advance_electron_source), its energy source
        (_deliverable_heating_j_m3), and the heavy-species chemistry call,
        which does nothing here: every row runs in the stage."""

        def __init__(self, p: Params1D):
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
            plasma_model = DriftDiffusionPlasmaModel(
                reactions=None,
                ion_mobility=load_ion_mobility(bridge.ion_mobility_path()),
                mixture=gas.normalized_fractions, gas_temp_k=p.tg0_K,
                area_m2=p.area_cm2 * 1e-4, gap_m=p.gap_cm * 1e-2,
                secondary_emission_gamma=0.0,
                energy_floor_ev=bridge.EPS_LO,
                energy_ceiling_ev=bridge.EPS_HI,
                field_heating_cap_factor=None,
                electron_energy_closure="source_coupled")
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

        # ------------------------------------------------------ ledger
        def _start_ledger(self):
            p, c = self.p1, self.grid_spec.cells
            self.tg = np.full(c, float(p.tg0_K))
            self._power = p.power_W_cm3 * 1e6 * p.profile()      # W/m^3
            z = np.zeros
            self.led = {          # per cell, eV m^-3, cumulative
                "input": z(c), "field_work": z(c), "elastic": z(c),
                "heat_reactions": z(c), "adjust": z(c), "to_e": z(c),
                "delivered": z(c),
                **{f"rad:{b}": z(c) for b in self.chem.bands},
                "krf_formed": z(c)}                   # m^-3, not energy
            self.e_other = 0.0            # eV: transport + optics, whole gap
            self.field_dropped = 0.0      # eV: negative field work not given
            self._t_interval = (0.0, 0.0)
            self._half_start = 0.0
            self._substeps = []
            self._capped = []     # (t, species, energy, source frequency)
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
            """Chemical + electron + gas thermal energy in the gap [eV]."""
            dens = self._combined_densities()
            u = (self.chem.chemical_energy(dens)
                 + self.plasma._last_n_eps / QE
                 + self.chem.heat_capacity(dens) * self.tg)
            return float(np.sum(u)) * self.cell_volume

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
            self.e_other += self.energy_total() - e0

        def _advance_optics(self, t, interval):
            e0 = self.energy_total()
            super()._advance_optics(t, interval)
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
            field = self._last_interval_field_work_j_m3
            pos = np.maximum(field, 0.0)
            dropped = float(np.sum(field - pos)) * self.cell_volume
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
            vol = self.cell_volume * dens.shape[0]
            if i in self._rec_steps:
                rad = sum(float(np.sum(self.led[f"rad:{b}"]))
                          for b in chem.bands) * self.cell_volume
                e_in = (float(np.sum(self.led["input"]))
                        + float(np.sum(self.led["field_work"]))) \
                    * self.cell_volume
                adjust = float(np.sum(self.led["adjust"])) * self.cell_volume
                e_tot = self.energy_total()
                self.records.append({
                    "t": t,
                    "mean": (dens.mean(axis=0) * 1e-6).tolist(),
                    "centre": (dens[mid] * 1e-6).tolist(),
                    "te_mean": float(np.sum(ne * eps)
                                     / max(np.sum(ne), 1e-300) / 1.5),
                    "te_centre": float(eps[mid] / 1.5),
                    "tg_mean": float(np.mean(self.tg)),
                    "tg_centre": float(self.tg[mid]),
                    "ledger": {k: float(np.sum(v)) * self.cell_volume / vol
                               * 1e-6 for k, v in self.led.items()},
                    "ledger_centre": {k: float(v[mid]) * 1e-6
                                      for k, v in self.led.items()},
                    "e_other": self.e_other / vol * 1e-6,
                    "e_total": (e_tot - self.e_total0) / vol * 1e-6,
                    "closure": ((e_tot - self.e_total0)
                                - (e_in - rad + adjust + self.e_other))
                    / vol * 1e-6,
                })
            if i in self._prof_steps:
                self.profiles.append({
                    "t": t, "dens": (dens * 1e-6).T.tolist(),
                    "te": (eps / 1.5).tolist(), "tg": self.tg.tolist()})

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
        print(f"LIGHT 1Dv0_2: {p.cells} cells over {p.gap_cm:g} cm, "
              f"{steps} steps of {p.dt_s * 1e12:g} ps", flush=True)
    out = sim.run(numerics, verbose=False)
    res = Result1D.from_simulation(sim, time.time() - t0)
    res.save(Path(out) / "result.json")
    return res


class Result1D:
    """A 1Dv0_2 run: time records, profiles and the energy ledger.

    Densities in cm^-3, energies in eV per cm^3 of gap (averaged over the
    gap, or the centre cell's own), temperatures in eV (Te) and K (Tg)."""

    def __init__(self, data: dict):
        self.data = data
        self.p = data["params"]
        self.species = data["species"]
        self.t = np.array([r["t"] for r in data["records"]])
        self.x_cm = np.array(data["x_cm"])

    @classmethod
    def from_simulation(cls, sim, wall_s):
        chem = sim.chem
        data = {
            "model": f"LIGHT 1Dv0_2 {__version__} (LIGHT v0_2 chemistry on "
                     "Fluid1D's 1-D transport)",
            "params": {k: (list(v) if isinstance(v, tuple) else v)
                       for k, v in asdict(sim.p1).items()},
            "species": list(chem.names),
            "bands": list(chem.bands),
            "x_cm": (sim.mesh.cell_centers * 100).tolist(),
            "records": sim.records,
            "profiles": sim.profiles,
            "ledger_cells": {k: (v * 1e-6).tolist()
                             for k, v in sim.led.items()},
            "field_dropped": sim.field_dropped / (sim.cell_volume
                                                  * sim.grid_spec.cells)
            * 1e-6,
            "substeps": {"mean": float(np.mean(sim._substeps)),
                         "max": int(np.max(sim._substeps))},
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
        if where == "ledger":
            # e_total is the change of chemical + electron + gas thermal
            # energy; the gas heat above is its thermal part
            stored = r["e_total"] - led["elastic"] - led["heat_reactions"]
            out["stored (species + electrons)"] = 100 * stored / e_in
            out["lost to the walls"] = -100 * r["e_other"] / e_in
            # positive: energy the mean-energy bounds took out (the wall
            # cells' drained electrons pinned at the ceiling)
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
                 f"  gap {p['gap_cm']:g} cm, {p['cells']} cells, dt "
                 f"{p['dt_s'] * 1e12:g} ps, chemistry substeps mean "
                 f"{d['substeps']['mean']:.1f} max {d['substeps']['max']}, "
                 f"wall {d['wall_s']:.0f} s",
                 f"  deposited {last['ledger']['input'] * QE:.3f} J/cm^3 "
                 f"(gap mean), field work "
                 f"{last['ledger']['field_work'] * QE:.2e} J/cm^3"]
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
        if e_in > 0:
            # the wall cells hold the unresolved sheath: Fluid1D's split
            # gives their electrons the field's positive work only and
            # never applies the negative part (its "phantom" work), which
            # grows as the mesh coarsens (README, "The wall cells")
            lines.append(
                f"  wall cells (unresolved sheaths): field work given "
                f"{100 * last['ledger']['field_work'] / e_in:.3g} %, "
                f"phantom negative work kept out "
                f"{-100 * d['field_dropped'] / e_in:.3g} % of the input")
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
