# F2-sweep-KrF: project charter

> Canonical repository instructions. `AGENTS.md` is a symlink to this file. Keep this document,
> `README.md`, `NUMERICAL_AUDIT.md`, `changelog.md`, and `bug-fixes.md` consistent with code changes.

## Purpose

This repository models discharge-pumped KrF* (248 nm) and ArF* (193 nm) excimer lasers while
sweeping the initial F2 mole fraction. The spatial model is one-dimensional along the electrode
gap. The optical gain axis is transverse and is represented by a zero-dimensional cavity model.

Cross sections feed coefficient solvers and reaction-rate tables. Keep the simulation focused on
the 1D charged-fluid, chemistry, circuit, and scalar-cavity formulation described below.

## Live physics

### Charged fluid and field

- `src/Fluid1D/plasma_model.py`: electron/ion continuity equations and electron-energy density.
- `src/Fluid1D/transport.py`: Scharfetter-Gummel drift-diffusion fluxes.
- `src/Fluid1D/transport_kernel.py`: Numba production loop; maintain parity with the NumPy path.
- `src/Fluid1D/field.py`: cell-centred Poisson solve, face-native E field, semi-implicit
  conductivity correction, prescribed-field superposition, fixed-voltage and MPC circuits.
- The MPC circuit is a seven-state linear ODE advanced by an augmented matrix exponential. The
  circuit network and finite-capacitance plasma draw use a first-order Lie split; transport returns
  the final node-B voltage, which is pinned back into the circuit.
- The production transport path is selected with `CircuitSpec(backend="mpc")`.

### Electron coefficients

`src/Fluid1D/boltzmann.py` defines `CoefficientTable` and the `BoltzmannSolver` interface.

- `backend="two_term"`: in-tree two-term energy-space Boltzmann solve using shipped cross sections.
- `backend="eedf"`: Maxwellian-EEDF moment approximation using the same measured data.
- `backend="stub"`: deterministic analytic coefficients for tests and smoke runs only.

No external LoKI-B or BOLSIG+ executable is part of the live architecture. Production LMEA tables
are keyed on local mean electron energy. Tables retain the paired E/N values for low-density local
field closure and can be re-keyed explicitly. Disk-cache identities include the full E/N grid,
mixture, gas temperature, key kind, complete solver configuration, cross-section contents, and
dynamic reaction signature.

Ion mobility is separate electron-independent data in `data/ion_transport/ion_mobility.toml`.
Mixtures use the available cited pairs through Blanc's law; absent pairs use a visible fallback.

### Chemistry

- `data/reactions/krf.toml`: 45-channel reduced KrF network.
- `data/reactions/arf.toml`: Chen-derived ArF network with cited corrections/additions.
- `src/reactions/`: data loading, layer routing, SI conversion, cross sections, and ion transport.
- `src/kinetics/`: combined charged/neutral mass action and per-cell linearly implicit integration.

Radiative and stimulated channels are owned by optics and excluded from the chemistry integrator.
Electron-impact and analytic Te-dependent reactions are evaluated at runtime. Energy-consuming
channels carry explicit loss/threshold data from one canonical reaction-energy table. The default
`Numerics(electron_energy_closure="source_coupled")` partitions electron/energy-transfer rows from
heavy chemistry, refreshes live-energy coefficients on the transport cadence, and accepts reaction
progress with its debit in one augmented implicit stage; the transport-ledgered local heating
(field work minus elastic, rectified per cell against the phantom-negative sub-Debye sheath
booking) is delivered INSIDE that stage as a per-substep power source so heating and the stiff
inelastic sink share one balance point, and the heavy partition's cadence grades on its
charged-species-changing frequency. The alternative `"continuous_reconciled"` is diagnostic only:
transport applies provisional continuous cooling and the chemistry interval refunds the prepaid
sink as an in-stage rate while recording a signed reconciliation. Numerical energy-bound
adjustments remain separate in both modes. Preserve stoichiometric update directions when
enforcing positivity; do not independently clip species unless no conservative alternative exists.

Optics applies PHOTO species feedback with the same exact frozen-fluence hazard that attenuates the
photon field. R40/R41 consume Ar**/F- and create their declared ion/electron products.
Photoelectron birth energies come from the data rows' `electron_energy_eV` (hv minus
binding/detachment, NIST/Blondel anchored) and are ledgered; the transport-owned electron energy
density is re-synced after every optics call.

### Optics and energetics

`src/optics/PhotonDensity` evolves Chen Eq. 7 as a scalar cavity-averaged density.
`GainCavity` co-integrates the lasing upper state with exact frozen-coefficient photon fluence,
saturated depletion, mirror/optical losses, and interval-integrated output power. Photon inventory
uses the full cavity volume `l_c*S`; upper-state inventory uses the gain volume `l_g*S`. Both output
and input powers recorded by `LaserEnergetics` are interval averages and must be integrated as
`sum(P_avg*dt_interval)`, including the first and final partial cadence.

`GainCavity.advance` ledgers its per-call upper-state sinks (`last_se_sink_m3` /
`last_ste_sink_m3`, hazard-share attribution so SE + STE equals the booked depletion exactly);
`ExcimerRates` exports them as the `SE`/`STE` quenching series beside the chemistry channels, so
the rate figures are faithful Chen Fig. 7 analogs even though optics owns the radiative rows.

The `EEDF` diagnostic reconstructs a Maxwellian from mean energy for visualization. It is not a
stored cell-resolved EEDF and must not be described as one.

The opt-in `EnergyClosure` diagnostic archives field work, elastic transfer, continuous sink,
accepted reaction cost, signed reconciliation, coefficient-bound adjustment, residuals, source
jumps, limiter activity, and closure-mode metadata. Never present a `continuous_reconciled` archive
as a production-conservative result.

## Declarative API

Experiments compose frozen specs and call `Simulation.run`:

```python
Simulation(
    name="arf",
    grid=Grid1D(cells=256, length=12e-3, area_m2=90e-4),
    gas=GasMixture({"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
                   pressure_atm=3.4),
    plasma=Plasma(density=1e16, ion="Ne+", electron_mean_energy_ev=2.0),
    preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13,
                                steepness=100.0),
    cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=0.99),
    boltzmann=BoltzmannSpec(backend="two_term", key_kind="mean_energy_eV"),
    circuit=CircuitSpec(backend="mpc", v0_v=-18e3),
    reactions="data/reactions/arf.toml",
    diagnostics=[OpticalGain()],
    plots=[LinePlot(
        filename="gain_history.png",
        x=DataRef("optical_gain/times_ps", "time (ns)", scale=1e-3),
        y=(DataRef("optical_gain/gain_per_m", "net gain"),),
        ylabel="optical gain (1/m)",
    )],
).run(Numerics(steps=4000, dt=50e-12, fluid_every=20))
```

`sweep_f2_parallel` is the production sweep entry point. Serial and parallel sweeps clone all
physical configuration, paths, units, drives, and preionization from the base simulation. Mutable
injected runtime objects are never shared between points: an injected plasma model or coefficient
table requires a per-F2 construction factory, and process-parallel factories must be picklable.

## Numerical invariants

- Internal arrays and coupling constants are SI.
- `CoefficientTable.key_grid` is finite and strictly increasing; consumers use its declared units.
- Seed electrons are balanced only by a known positive ion, with `n_i=n_e/z_i`.
- Prescribed axial E fields must be applied at cells and faces.
- Finite gap-capacitor energy transfer cannot exceed its stored-energy change during a gap-only
  split.
- Accepted chemistry progress and its electron-energy debit must come from the same substep trial;
  any coefficient-table-bound energy adjustment must remain separately visible. The per-cell
  identity is `E_final + E_cost − E_bound_adjustment − E_delivered = E_initial`, with delivered
  stage heating ledgered (`delivered_source_j`) and the rectification discard ledgered
  (`undelivered_field_work_j`).
- `Numerics(field_heating_cap_factor=None)` is a sensitivity control for the audited `J.E` cap,
  not a production shortcut around field/sheath energy closure.
- Stimulated cavity photon production cannot exceed the corresponding gain-volume upper-state
  loss, apart from explicit spontaneous sources.
- Every fluid interval, including a final remainder, contributes to chemistry, optics, and energy.
- Kernel changes require the NumPy parity test and `NUMBA_BOUNDSCHECK=1` coverage.

## Diagnostics

Diagnostics are opt-in numeric collectors that write `diagnostics.h5` under `runs/<name>/`.
Collectors never choose figures or filenames. Run files own all archive-backed plot declarations,
including mandatory names, x/y datasets, scale flags, and optional fixed ranges. Unpinned axes use
finite-data-aware limits with small padding. Domain-mean histories and spatial x-t maps are
distinct diagnostics; describe them accurately.

## Repository map

```text
src/Fluid1D/    orchestrator, specs, field/circuit, transport, coefficient solvers, sweeps
src/kinetics/   combined chemistry, implicit integration, energy ledger
src/optics/     scalar photon density and saturated gain/output
src/reactions/  network parser, cross-section loader, ion transport
src/particles/  species and constants
src/diagnostics archive and plotting helpers
src/backend/    Numba threading and topology
data/           networks, cross sections, mobility data, fixtures
examples/       ArF/KrF single runs and parallel F2 sweeps
tests/          unit, parity, conservation, and smoke tests
```

## Development

```bash
python -m pip install --no-deps -e .
python -m unittest discover -s tests -p "test_*.py"
NUMBA_BOUNDSCHECK=1 python -m unittest tests.test_transport_kernel
python examples/arf_f2_sweep.py --tiny
python examples/arf_convergence.py --resume
```

Use writable `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` in restricted environments. Add focused tests for
numerical fixes and preserve unrelated working-tree changes.
