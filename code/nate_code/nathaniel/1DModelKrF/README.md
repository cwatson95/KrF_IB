# F2-sweep-KrF

A one-dimensional charged-species fluid model for discharge-pumped **KrF\*** (248 nm) and
**ArF\*** (193 nm) excimer lasers, with F2 concentration as the principal sweep parameter.

Electrons and ions are density fluids on a 1D electrode-gap mesh. Electron coefficients are
computed from the shipped cross sections by an in-tree two-term Boltzmann solver (or a faster
Maxwellian approximation), cached as lookup tables, and consumed by the transport and chemistry
models. Optics is transverse to the discharge mesh and is represented by one cavity-averaged
photon density.

See the [full codebase guide](docs/index.html) for the current implementation map and
[CLAUDE.md](CLAUDE.md) for the architecture charter and numerical invariants.

## Model

- **Charged transport:** per-cell electron and ion densities solve
  `dn_s/dt + dGamma_s/dx = S_s` with Scharfetter-Gummel drift-diffusion fluxes. A separate
  electron-energy density solves transport, field heating, and elastic cooling.
- **Electrostatics:** a cell-centred 1D Poisson solve with Dirichlet electrode potentials. The
  conducting phase uses a semi-implicit conductivity correction and face-native electric fields.
- **Discharge circuit:** either fixed electrode potentials (`StubCircuit`) or the seven-state MPC
  network (`MpcCircuit`). The MPC network and finite gap capacitance are advanced with a split
  circuit/plasma coupling.
- **Chemistry:** combined charged and neutral mass-action kinetics, advanced per cell with a
  linearly implicit stiff integrator. Electron-impact rates are refreshed from the live local mean
  energy inside each chemistry substep, and reaction progress is accepted with its exact
  threshold-energy debit from the same trial.
- **Electron coefficients:** `TwoTermBoltzmannSolver` is the production backend;
  `MaxwellianEedfSolver` is a faster measured-cross-section approximation;
  `AnalyticStubSolver` exists only for tests and fast smoke runs. No external LoKI-B or BOLSIG+
  executable is required.
- **Optics:** `PhotonDensity` and `GainCavity` evolve a 0D cavity density, spontaneous seeding,
  saturated stimulated depletion, absorption, mirror loss, and interval-averaged output power.
  The optical gain axis is transverse to the simulated discharge axis.

The engine computes in SI. `Simulation(volume_units="cm" | "m")` changes only geometry and
density input/output conventions.

## Setup

Python 3.11+ is required. Runtime dependencies are NumPy, SciPy, Numba, Matplotlib, and h5py;
Linux installations may also use TBB for Numba threading.

```bash
./MacOS-setup.sh       # macOS
./Unix-setup.sh        # Linux
# or .\Windows-setup.ps1 on Windows
source .venv/bin/activate
```

The setup scripts install the dependencies and then install this repository as an editable
`f2-sweep-krf` package from `pyproject.toml`. A manual environment can use:

```bash
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
```

## Run

Run these commands from the repository root after activating the installed environment.

```bash
# Fast three-point ArF development sweep
python examples/arf_f2_sweep.py --tiny

# Production ArF F2 sweep: two-term coefficients + MPC transport
python examples/arf_f2_sweep.py

# One ArF point or one KrF point
python examples/example-arf.py
python examples/01_krf_excimer.py

# KrF sweep, with the same tiny/production modes
python examples/krf_f2_sweep.py --tiny
python examples/krf_f2_sweep.py

# Test suite
python -m unittest discover -s tests -p "test_*.py"

# Production ArF mesh/timestep gate (128/256/512 cells; 100/50/25 ps)
python examples/arf_convergence.py

# Resume a partially completed production gate from per-member checkpoints
python examples/arf_convergence.py --resume
```

Production coefficient builds are cached per mixture, complete E/N grid, solver configuration,
cross-section contents, and reaction-network signature. F2 sweep points can run concurrently with
`sweep_f2_parallel`; each worker is pinned to one Numba thread.

### Live reporting

Single runs print an opt-in progress line (`Simulation.run(progress=ProgressReport(...))`,
default fields `excimer, photons, temp:e-, efficiency` plus trailing `elapsed=`/`eta=` columns;
the `excimer` field auto-resolves to the network's lasing upper state, KrF\*/ArF\*). Parallel
sweeps report live by default (`sweep_f2_parallel(..., progress=True)`): a start banner, one
completion line per finished point, a periodic `[sweep] done k/N running … eta~…` rollup, and an
end summary. A point that raises is recorded on `SweepResult.failures` (message + worker
traceback) while the remaining points complete and the overlay figures render from the
survivors; `progress=False` restores the fully silent behavior.

Each sweep worker also rewrites `runs/<name>/f2_<value>/progress.json` atomically on its report
cadence (status `running`/`done`/`failed`, step, %, sim time, wall time, field values). Watch a
live or crashed sweep from another terminal:

```bash
python -m Fluid1D.sweep_status runs/krf_f2_sweep --watch 5
```

`running` heartbeats older than `--stale` seconds (default 120) are flagged `stalled` — the
post-mortem trail for dead or interrupted workers.

## Declarative API

```python
from Fluid1D import (
    BoltzmannSpec, Cavity, CircuitSpec, DataRef, GasMixture, Grid1D,
    LinePlot, Numerics, OpticalGain, Plasma, Preionization, Simulation,
    sweep_f2_parallel,
)

sim = Simulation(
    name="arf_baseline",
    grid=Grid1D(cells=256, length=12e-3, area_m2=90e-4),
    gas=GasMixture(
        {"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
        pressure_atm=3.4,
    ),
    plasma=Plasma(density=1e16, ion="Ne+", electron_mean_energy_ev=2.0),
    preionization=Preionization(
        mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0,
    ),
    cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=0.99),
    boltzmann=BoltzmannSpec(
        backend="two_term", key_kind="mean_energy_eV", cache_dir="runs/lut_cache",
    ),
    circuit=CircuitSpec(backend="mpc", v0_v=-18e3),
    reactions="data/reactions/arf.toml",
    diagnostics=[OpticalGain()],
    plots=[LinePlot(
        filename="gain_history.png",
        x=DataRef("optical_gain/times_ps", "time (ns)", scale=1e-3),
        y=(DataRef("optical_gain/gain_per_m", "net gain"),),
        ylabel="optical gain (1/m)",
    )],
)
sim.run(Numerics(steps=4000, dt=50e-12, fluid_every=20))
```

`sweep_f2` and `sweep_f2_parallel` clone every physical configuration field from a base
`Simulation`, replace only the F2 fraction (and the selected balance gas), run independent states,
and return output energy plus intrinsic efficiency for each point. Both are exported from the
top-level `Fluid1D` package. Select parallel points with either an inclusive
`start`/`stop`/`step` range or a unique explicit `values` sequence; empty, duplicate, out-of-range,
or ambiguous selections are rejected before workers start. Diagnostics and plots are constructed
per point by optional factories. If the base simulation contains an injected mutable
`plasma_model` or `coeff_table`, supply the corresponding `plasma_model_factory(f2)` or
`coeff_table_factory(f2)`; parallel factories must be picklable, and the sweep rejects an injected
object that cannot be recreated independently for every point.

## Outputs

Diagnostics are opt-in numeric collectors that write `diagnostics.h5` under `runs/<name>/`.
Plots are separately opt-in and declared in the run file with mandatory filenames and explicit
x/y dataset references. Run authors can select linear, log, or symlog axes and fixed ranges; an
unpinned axis fits all finite visible data with a small margin instead of rounding to a decade.
This keeps plot selection, naming, and scale policy out of the simulation runtime.

Available datasets include domain-mean species histories, charged-density x-t maps, field
profiles, circuit traces, E/N and electron-temperature traces, reaction-rate histories, gain,
laser power, and energy efficiency. Pulse width is derived from archived laser power with
`pulse_fwhm`. The `EEDF` collector exports a qualitative
Maxwellian reconstructed from the mean-energy field; the runtime state does not retain a resolved
energy distribution in every cell. `EnergyClosure` adds the per-interval field, elastic,
continuous-inelastic, accepted-reaction, signed-reconciliation, and bound-adjustment ledgers plus
residual/limiter histories and the selected closure mode. See `examples/example-arf.py` and
`examples/01_krf_excimer.py` for complete collector and plot lists.

The supported outputs are static figures and the HDF5 diagnostic archive.

## Current Limits

- PHOTO feedback uses the same exact frozen-fluence hazard for photon attenuation and species
  mutation. R40/R41 consume Ar**/F- and create their declared products; their photoelectrons enter
  the electron-energy density at 3.76/3.02 eV. The analogous KrF R43/R44 values are 3.10/1.60 eV.
  These are sourced point-energy assignments (`hv` minus binding/detachment), not resolved
  photoelectron energy distributions.
- Ion mobility data cover cited ion/background-gas pairs. Missing pairs use the documented
  fallback reduced mobility; mixture values are only as complete as the available pair data.
- `Numerics.electron_energy_closure` defaults to production `source_coupled`;
  `continuous_reconciled` is causal diagnostics only and is marked as such in HDF5/terminal output.
  Historical failures that motivated the current in-stage heating source are not current
  validation results. Production claims still require a fresh mesh/timestep gate at the operating
  point.
- `AnalyticStubSolver` and fixed-voltage smoke runs test plumbing only; they are not predictive.
- ArF reaction data are a Chen-derived, cited network with documented corrections and additions,
  not a byte-for-byte transcription of Table III. `arf_chen_2025.toml` is a reproduction-only
  overlay of the printed rows and is never the production default.

## Layout

```text
src/Fluid1D/    simulation orchestration, specs, fields, transport, coefficient solvers, sweeps
src/kinetics/   combined chemistry, implicit integrator, energy ledger
src/optics/     0D photon-density and saturated-gain model
src/reactions/  reaction parsing, cross-section loading, ion transport
src/particles/  species registry and physical constants
src/diagnostics HDF5 archive and plotting helpers
src/backend/    CPU topology and Numba threading selection
data/           ArF/KrF networks, cross sections, cited ion mobility data, test fixtures
examples/       single-point and parallel F2-sweep programs
tests/          unit, conservation, parity, and end-to-end regression tests
```
