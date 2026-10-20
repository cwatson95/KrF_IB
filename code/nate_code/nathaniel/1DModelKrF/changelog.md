# Changelog

Major commits, newest first. One entry per significant change set. Keep entries terse; link to
`bug-fixes.md` for defect fixes.

Format:
```
## YYYY-MM-DD — <short title>
- <what changed and why, at a high level>
```

---

## 2026-08-01 — Dead-code removal, runtime decomposition, and public-contract cleanup

- Removed unreferenced transport/circuit/fluid/spec/reaction/cross-section compatibility surface
  and dead telemetry, plus unused vector/elastic-cross-section fixtures, the unwired Ar 4p
  excitation dataset, and tracked macOS metadata. Each removal was committed independently before
  the aggregate removal test gate.
- Decomposed simulation construction/interval coupling and the NumPy/Numba transport stages while
  preserving parity contracts. Centralized shared physical constants, normalized reaction-unit
  whitespace, preserved radiative/stimulated product multiplicities, and exposed copy-returning
  diagnostic accessors used by the analysis programs.
- Added editable `pyproject.toml` packaging and a locked environment description; platform setup
  scripts now install the repository itself, so examples and tests import the installed package
  instead of mutating `sys.path`.
- Serial and parallel F2 sweeps now share one point-construction path, preserve every representable
  base configuration field, validate range/value selections and worker counts up front, and require
  per-point factories for injected plasma models or coefficient tables (picklable for process
  workers). `sweep_f2_parallel` is part of the top-level package API. The ArF/KrF tiny entrypoints
  pass their explicit three values without also forwarding the production range defaults.
- Refreshed current-code comments and canonical docs, including PHOTO product/birth-energy behavior,
  while retaining older numerical results only as labeled historical checkpoints. No final test
  count is recorded here; use the acceptance-run output for the completed change set.

## 2026-08-01 — Live reporting: sweep parent lines, heartbeat files, status CLI

- `sweep_f2_parallel(..., progress=True)` (new default) now reports live: start banner, one
  completion line per finished point (E_out/eff/wall/ETA), a periodic in-flight rollup built
  from worker heartbeats, and an end summary. Failures are isolated: a raising point is
  recorded on the new `SweepResult.failures` (message + worker traceback, sentinel-returned so
  it can never poison the pool; `BrokenProcessPool` folds in as "victim or collateral") and the
  rest of the sweep completes — previously one bad point aborted `pool.map` and lost everything.
  Ctrl-C cancels pending futures cleanly. `progress=False` restores full silence; a
  `ProgressReport` customizes worker cadence/fields (forced heartbeat-only — the parent owns
  the terminal).
- `ProgressReport` gains `timing` (trailing `elapsed=`/`eta=` columns, default on), `terminal`,
  and `heartbeat` sinks; workers/single runs with `heartbeat=True` atomically rewrite
  `<run_dir>/progress.json` (schema 1: status/step/pct/t_sim/wall/fields) each cadence plus a
  final `done`/`failed` record. `DEFAULT_PROGRESS_FIELDS` is now system-aware: the new
  `excimer` field resolves to the network's lasing upper state (KrF*/ArF*) instead of
  hard-coding KrF*.
- New `python -m Fluid1D.sweep_status runs/<name> [--watch N] [--stale S]` renders a per-point
  status table from the heartbeats (`stalled` flags dead workers).
- Tests: 3 new suites/extensions (progress heartbeats + timing/legacy formats, sweep failure
  isolation incl. a BrokenProcessPool probe, status-CLI scan/render/return-codes).

## 2026-08-01 — Progress-line cleanup: explicit spec reports on quiet runs

- `Simulation.run` progress lines are decoupled from `verbose`: an explicit `ProgressReport`
  now reports even under `verbose=False` (banner/summary stay quiet), and fully silent runs no
  longer construct-and-validate a default reporter. `ProgressReporter.due` always reports the
  final step, so uneven cadences still close on a 100% / N-of-N line. Removed the write-only
  `Simulation._verbose` attribute; hoisted `sweep_f2_parallel`'s function-local imports.
- Tests: uneven-cadence final line, quiet-run explicit progress, fully-quiet run prints nothing.

## 2026-08-01 — In-stage heating source: the closure ignites at Chen scale

- Delivered the transport-ledgered local heating (field work − elastic, RECTIFIED per cell
  against the phantom-negative sub-Debye sheath booking) inside the co-integrated electron
  stage as a per-substep power source, Strang-split with exact per-interval conservation
  (`delivered_source_j` / `undelivered_field_work_j` ledgers; identity extends to
  `E_final + E_cost − E_bound − E_delivered = E_initial`). The stage now relaxes to the
  J·E-balanced ε_eq(E/N) instead of draining to the floor: 256c/50 ps/0.15% two-term went from
  n_e 8.2e19 / 0.1–2.2 kA / µJ (with +0.277 J of floor-minted clamp energy) to
  **n_e 1.6e20 / 4.0 kA / 1.45 J deposited (Chen ≈1.36) / mJ-scale lasing**, with bound
  adjustments collapsed ×60 to +4.5 mJ (0.8% of the 0.56 J accepted cost) and the avalanche
  riding T_e 4.35–4.60 eV.
- Repaired `continuous_reconciled`: the prepaid sink is refunded as an in-stage rate (the old
  lump-sum add was destroyed above the ceiling by the pre-clamp — the ×1278 divergence); the
  sink's row/column packing is pinned correct against `energy_consumption_rate_percell`.
  Cross-mode agreement test added (pre-fix collapse was ×13).
- Graded the heavy-partition cadence on a new `charged_changing_mask` (ion-ion recombination
  ran as 1 ns impulses); photoelectron birth energies (`electron_energy_eV` on PHOTO rows,
  NIST/Blondel anchored) now enter the pool and the ledger with `_last_n_eps` re-synced after
  optics; γ exposed on the `Plasma` spec; examples aligned to Chen's 100% HR (`r_back=1.0`).
  See `bug-fixes.md` (five entries, 2026-08-01).

## 2026-07-31 — Explicit ArF electron-energy closure modes and failed ignition gate

- Added production-default `source_coupled` and diagnostic-only `continuous_reconciled` modes,
  one canonical compiled reaction-energy table, electron/heavy source partitions, augmented-stage
  extent/debit acceptance, signed reconciliation, and the opt-in `EnergyClosure` HDF5 ledger.
- Preserved R9/R18 energy metadata, documented R10's zero-threshold limitation, cached open-gap
  MPC exponentials, corrected ArF sweep CE to the 1.62 J bank denominator, and added a
  reproduction-only printed Chen network plus a non-fabricated published-reference fixture.
- Added exact R40/R41 PHOTO species feedback with photon/atom/charge conservation. This checkpoint's
  explicit zero-photoelectron-birth-energy assumption was superseded on 2026-08-01 by sourced and
  ledgered ArF 3.76/3.02 eV and KrF 3.10/1.60 eV point birth energies.
- The staged 128-cell/50 ps/two-term production gate reaches `8.20e19 m^-3`, 2.24 kA, and 0.888 J
  deposition but only 1.21 microjoules output and a +0.277 J bound adjustment. The diagnostic
  reconciles +1.99 J against 1.56 mJ accepted cost. Both fail; no 256/512 or calibration run was
  launched. This is a historical, pre-heating-source checkpoint rather than a current result.

## 2026-07-31 — Co-integrated chemistry/electron energy + executable convergence gates

- Removed electron-impact reaction losses from the transport PDE and added a Numba-parallel
  chemistry path that refreshes live-energy rates, accepts reaction progress with its exact
  threshold debit, and exposes coefficient-bound energy adjustments separately.
- Added per-cell closure/dynamic-cooling regressions and `examples/arf_convergence.py`, whose
  production matrix gates 128/256/512 cells and 100/50/25 ps, checkpoints/restarts expensive
  members, and writes a JSON audit record. Added an explicit field-heating-cap sensitivity control.
- Fixed the zero-cost limiter endpoint so threshold-free Penning electron production is retained
  and any energy-floor adjustment is visible. The full 0.15% matrix fails all mesh/timestep gates:
  peak `n_e` ranges `3.34e18`--`2.22e19 m^-3` and output 5.58--48.7 nJ. Even the 512-cell mesh has
  `dx/lambda_D ~= 5.3`; chemistry-cadence and heating-cap sensitivities are near-null. This matrix
  is retained here as historical defect evidence and does not validate the current closure.
- A corrected-closure 128-cell sensitivity at 0.20/0.233/0.267/0.30% F2 remains entirely
  nonlasing (`n_e = 2.70e18`--`3.43e18 m^-3`, output 3.16--4.79 nJ). Its apparent 11--14 ns
  widths are spontaneous leakage and are explicitly excluded from Chen FWHM validation.

## 2026-07-31 — SE/STE quench series + Chen Fig. 15 density panels

- `GainCavity.advance` now ledgers its per-call spontaneous/stimulated upper-state sinks
  (`last_se_sink_m3`/`last_ste_sink_m3`, competing-risks hazard split — SE + STE equals the
  booked depletion exactly, unit-tested). `ExcimerRates` samples the ledger and exports `SE`/
  `STE` beside the chemistry quench channels, making the rate figures faithful Chen Fig. 7(b)
  analogs (the two radiative sinks are the #1/#2 quench channels and were previously invisible,
  which made the collisional curves read as if they were SE/STE).
- Added `OverlayPanel`/`PanelOverlayFigure`/`plot_sweep_panels` to `sweep_plots` and a
  `sweep_species_densities.png` declaration (Chen Fig. 15: 2×2 Ar⁺/Ar\*/ArF\*/F⁻ — Kr analogs
  for KrF — one curve per F₂ value, shared legend) to both F₂ sweep run files, fed by the
  existing `SpeciesDensities` collector.

## 2026-07-31 — Data-aware, run-owned diagnostic plots

- Split diagnostics into numeric collectors and archive-backed `LinePlot`/`HeatmapPlot`
  declarations. Experiment files now own every figure filename, x/y dataset selection,
  linear/log/symlog scale, and optional fixed range; the simulation runtime supplies no plots.
- Replaced decade-rounded and symmetric axes with finite-data limits plus small padding, including
  y-fitting inside explicit x windows. Updated both single-run examples and F2 sweep plot lists.
- Corrected circuit/rate archive time keys from `times_ps` to `times_ns` and field position from a
  time-like dataset to the experiment's declared length unit. Added line, heatmap, naming, scaling,
  and sweep regressions. See `bug-fixes.md`.

## 2026-07-31 — Numerical audit and fluid-only documentation reset

- Audited transport, field/circuit coupling, coefficient caching/key units, cross-section
  interpolation, chemistry energy bounds, saturated optics, power integration, final cadences,
  sweeps, and spec validation. Fixed the Numba Thomas-workspace overrun, within-step chemistry
  overspend, optical inventory minting, double node-B advance, prescribed-face omission, LUT cache
  collisions, key-unit confusion, interval-energy loss, and partial-cadence loss.
- Added conservation and parity regressions for bounds-checked transport, gap-capacitor energy,
  chemistry budgets, high-gain photon inventory, cache identities, and coefficient re-keying.
- Rewrote the README, canonical charter, HTML guide, and active source guides around the current
  1D fluid + 0D optics architecture. Removed source-catalog material outside the active model and
  kept the provenance map focused on current physics. See `NUMERICAL_AUDIT.md` and `bug-fixes.md`.

## 2026-07-31 — Ne-buffered KrF* network and full F2 sweep

- Replaced the KrF starter subset with a cited, atom/charge-balanced 45-reaction Ne/Kr/F2
  network: complete Ne/Kr/F2 electron channels, Penning and dimer kinetics, three KrF* source
  paths, five measured/estimated quench paths, molecular-intermediate radiative sinks, and
  species-resolved 248 nm absorption. Corrected the F2+ ionization products, Kr2+/F- branch,
  F2 quench products, two-body Penning units, and several mis-cited rates (`bug-fixes.md`).
- Expanded the KrF species registry for Ne, F2(a3), F2+, Ne2*/Ne2+, Kr2F*, and non-lasing
  radiative photon modes. The 248 nm upper-state model now consistently pairs the measured
  KrF(B) 9 ns lifetime with its 2.6e-16 cm2 stimulated-emission cross section.
- Converted `examples/01_krf_excimer.py` to the three-gas Ne-buffered case and added
  `examples/krf_f2_sweep.py`, the Chen-style 0.05--0.50% F2 production sweep with per-point
  archives, rate/circuit plots, sweep overlays, FWHM/efficiency, and alpha/eta figures.
- Completed the 10-point real two-term/MPC production sweep. The lasing window peaks at 0.10%
  F2 (20.32 mJ, 1.382% intrinsic efficiency, 14.9 ns FWHM), remains open at 0.15%
  (17.17 mJ), and closes by 0.20% as attachment suppresses avalanche (peak n_e drops from
  1.13e20 to 2.64e16 m^-3). Results and seven aggregate figures are in `runs/krf_f2_sweep/`.

## 2026-07-31 — Ma et al. verification: R26 confirmed, three Chen exponent typos fixed (R23/R24/R36)

Chen's Table III cross-checked row-by-row against Ma 2024 Table 2 (PSST **33** 075012 — Chen's
declared primary heavy-particle source, now in `sources/kinetics/`). **R26 verified** as
Ar₂⁺ + Ne dimerization (Ma prints the same "Ne⁺+2Ar" typo, but his own network forms Ar₂⁺
nowhere else while consuming it twice, and Chen's Fig. 2 draws the Ar⁺→Ar₂⁺ arrow). Three
Chen prints corrected to their cited sources' values: R23 2.5×10⁻⁴³→1.1×10⁻⁴⁴ (Ar\* donor
drained ×23 too fast), R24 4×10⁻⁴⁰→8×10⁻⁴⁶ m⁶/s (cm⁶ slip — Ne\* lifetime was 0.4 ps,
killing the Ne-bath Penning/stepwise chain), R36 3.5×10⁻⁴⁵→3.5×10⁻⁴³ (real ArF\* sink).
Kept where Chen cites non-Ma sources: R29/R30/R31. 200 ns: deposited 1.364 J (Chen 1.36),
P_L 0.635 MW, **FWHM 15.0 ns (Chen 15)**, E_out 11.3 mJ (was 2.6), η 0.70% of stored
(Chen 5.8%), n_e 4.2×10²⁰.

## 2026-07-31 — Cross-section sensitivity campaign: BSR branching, Ne manifolds, R39 cascade

Autonomous hypothesis testing against the ranked lasing-gap list (all data derived from the
crossbrowser LXCat exports, no fabrication): (1) **R5/R6 branching** rebuilt with the BSR
(B-spline R-matrix) 4p/(4s+4p) ratio — f = 0.25–0.36 in the 13.5–30 eV window vs RDW's 0.68
(beam-experiment range), ArF\* +29 %, laser ×26; a fixed f = 0.25 bound run confirms the
lever. (2) **R1/R2 Ne records** = SIGLO manifold sums (previously only the first level of
each manifold — half the 96 %-bath sink discarded); near-null on lasing, kept for
faithfulness. (3) **R39 read as the 4p→4s cascade** (printed products E1-forbidden,
R26-class typo; `bug-fixes.md`) — the Ar\*\* flux returns to Ar\*: gain 3.62 m⁻¹, E_out
2.6 mJ (×500 over the campaign start). R4's SIGLO record cross-validated against the
independent TRINITI state-mean (≤2 % above 8 eV); R8 has no alternative dataset (TRINITI
unique). Remaining Chen gap concentrated in ArF\* ×3.4.

## 2026-07-31 — R6 Ar\*\* cross-section: NGFSRDW ground→4p sum closes the last ArF data gap

- `data/cross_sections/e_Ar_excitation_2.{json,csv}`: e + Ar → Ar\*\*(4p), the sum of the ten
  NGFSRDW (LXCat, retrieved 2026-07-30) state-resolved ground→4p blocks (Paschen 2p1–2p10,
  thresholds 12.907–13.480 eV; RDW — Kaur et al., J. Phys. B **31**, 4833 (1998)) on the union
  energy grid (29 pts, 12.9–500 eV, σ_peak 3.44×10⁻²⁰ m² at 13.4 eV). SIGLO (Chen ref 22) has no
  second Ar excitation level. R6 `threshold_eV` 13.1 → 12.907 (data). Every electron-impact row in
  `arf.toml` now imports a dataset (0 failures). Caveats documented in the record + `arf.toml`:
  R5's SIGLO set is "all excitation grouped into one level" so R5+R6 overlap in the 4p share, and
  RDW is threshold-peaked (σ_max at 13.4 eV vs beam-experiment 4p peaks near ~20 eV).
- Physics shift on the 200 ns two-term ArF run (arf-fixcheck2 → arf-r6): the 4p sink eats the
  12.9–20 eV tail and reroutes deposited energy to Ar\*\*→hv95 — breakdown-deposition peak moves
  88 → 115 ns (Chen 100–110 ns), same 1.65 J deposited, but ArF\* peak 9.2×10¹⁹ → 3.6×10¹⁹ m⁻³
  (Chen 2.5×10²⁰), peak gain 3.91 → 1.40 %/cm: below cavity threshold, photon build-up dies —
  E_out 19.4 mJ → 0.47 µJ, η 1.19×10⁻² → 2.9×10⁻⁷ (Chen 75 mJ / 5.8×10⁻²).
- Tests re-pinned (234 pass): two-term vs Maxwellian ν_ion compensation now holds ≥30 Td only
  (ratio 0.017/0.049 at 10/20 Td — the new channel's tail depletion); the missing-cross-section
  warning machinery is exercised synthetically and ArF full coverage asserted.

## 2026-07-31 — Energy-matched secant gap conductance; books close end-to-end; ArF* over threshold

The circuit's gap load was `G = A·σ_mean/d` (bulk) while the real gap current is
sheath-limited after collapse — the circuit dissipated 1.64 J into a fictional resistor while
the plasma absorbed 0.154 J (the ×9 deposited-energy deficit vs Chen that starved ArF\* below
the lasing threshold). Both coupled paths now hand the circuit the energy-matched secant
conductance `G = max(∫i_gap·v_b dt, 0)/∫v_b² dt` per advance — passive, and circuit-side
dissipation equals true plasma absorption by construction (`bug-fixes.md`). 200 ns two-term
gate: delivered +1.599 J / absorbed +1.607 J / bank left 0.001 J (Chen deposits 1.36 of
1.620 J); collapse from −18.7 kV @ ~89 ns; I_s peak 7.7 kA; mean n_e peak 5.5×10²⁰;
T_e(mean) peak 5.28 eV (Chen ~5 eV); ArF\* peak 7.3×10¹⁹ m⁻³, above the 3.6×10¹⁹ threshold
for 30 ns. Full 200 ns run: ~3 min wall.

## 2026-07-30 — Rectified-heating cap closes the collapse-window mint; discharge circuit-clamped

The last energy leak is closed: discrete J·E work in sub-Debye sheath cells oscillated
±10¹¹–10¹³ W at the collapse and was *rectified* into licence (per-cell `max(p,0)` + n_ε
accumulation), funding runaway ionization to 10²⁶ m⁻³. Both rectified totals Σ±p_field are
now capped per substep at the gap's true electrical power 2·|I_gap·V_b| (exact identity
pre-breakdown; kills the mint at the collapse) in the numba kernel and the numpy coupled
path (`bug-fixes.md` entry). 95 ns two-term verification: nₑ clamps at 7.9×10²¹ (Chen's
circuit-limited scale), first collapse −17.8 kV @ 84.4 ns (Chen band −16…−19 kV), I_s
through −2.9 kA, ionization spend 0.020 J on the 1.62 J bank, 106 s wall for 1900 steps.
Known residual: ~−1.5 J of booked negative sheath work across the ringing (capped at supply
scale; cannot license ionization; λ_D-resolution refinement documented).

## 2026-07-30 — Energy-ledger closure + numba kernels + two-term solver + dead-code purge

The end-to-end energy books now CLOSE: circuit→gap +8.5 mJ vs plasma field-heating +7.1 mJ
over a 65 ns breakdown (previously the ledger blew up by 3–9 orders through a chain of
coupling leaks, each found by an instrumented probe and fixed in turn — full chain in
`bug-fixes.md`). A 65 ns run costs ~60 s wall (was multiple hours).

- **Semi-implicit field coupling** (`Poisson1DSolver.solve_semi_implicit`, `SemiImplicitField`):
  σ-corrected tridiagonal re-solved per transport substep — removes the ε₀/σ explicit-coupling
  instability (nₑ ≳ 3e19 at 50 ps) and the dielectric-relaxation substep cost.
- **Chemistry pays its bill**: per-cell inelastic spend accumulated in-kernel (trapezoidal),
  debited from n_ε (`debit_electron_energy`); adaptive per-step cadence at breakdown; and a
  **hard per-cell energy budget inside the integrator** — costly channels stop the substep
  their licence (field power + stored energy) runs out. The throttle shapes rates; the budget
  makes the licence exact.
- **Exact circuit**: the MPC 7-state advance is an 8×8 augmented matrix exponential —
  unconditionally stable at any gap conductance (RK4 sub-stepping was hours-slow resolved or
  V_b→∞ capped). The gap enters as a PASSIVE conductance load (an offset current source
  pumped unbounded energy through stale exchanges), and the kernel co-advances a local node-B
  copy per substep so the plasma's Dirichlet boundary sags self-consistently.
- **Climax exchange subdivision**: the whole circuit↔field↔transport↔chemistry handshake
  subdivides (≤256×) when the avalanche e-folds faster than 50 ps.
- **Numba transport kernel** (`Fluid1D/transport_kernel.py`): the entire coupled substep loop
  in one @njit kernel (~40× per substep), parity-pinned against the numpy path
  (`tests/test_transport_kernel.py`); numpy path retained as the reference.
- **Two-term Boltzmann solver** (`Fluid1D/two_term.py`, `backend="two_term"`, now the ArF
  example default): H&P Eq. 39 energy-space SG solve, conservative inelastic deposits,
  temporal-growth **Perron eigenpair** by shifted inverse-power iteration (mandatory in an
  attaching mixture), sparse Δε≈0.05 eV grid (the 16.6–21.6 eV wall is exponentially
  load-bearing), column-scaled tail. 12 exact-limit validation tests. Empirical finding
  (pinned): its hotter ⟨ε⟩(E/N) largely compensates the tail depletion — ν_ion(E/N) is nearly
  closure-invariant; the measured breakdown delay vs the Maxwellian backend is ~15 ns.
- **Optics co-integration**: GainCavity sub-steps the (n_ph, upper-state) pair (growth-limited),
  Eq. 18 large-gain asymptote + exp clamps — the 100×-past-threshold overflow scenario is now a
  regression test.
- **Parallel F₂ sweep**: `sweep_f2_parallel(start, stop, step)` — process-pool batches, workers
  pinned to one numba thread. Threaded two-term LUT build.
- **arf.toml R26 → Ar₂⁺ + Ne** (paper's printed products endothermic by +5.8 eV; documented
  deviation, Ma et al. PSST 33 075012 (2024) verification pending); Xe⁺ orphan noted.
- **Dead-code purge (~3,400 LOC + 8 MB)**: unused external-solver, field-I/O, optics,
  diagnostics, adapter, and setup-script components; package inits slimmed. Tests 246 → 223
  (dead subjects retired; new kernel/energy/two-term suites added).

## 2026-07-30 — Semi-implicit field coupling + chemistry energy debit (the discharge survives breakdown)

The staggered-explicit Poisson↔transport coupling was linearly unstable past nₑ ≈ 3×10¹⁹ m⁻³
at dt = 50 ps — it manufactured field energy (measured 1.45×10⁹ J against the 1.62 J bank),
pinned ⟨ε⟩ at the LUT ceiling, fully ionized the domain, and NaN'd. Fixed in two parts
(diagnosis + full measurement chain in `bug-fixes.md`):

- **Semi-implicit σ-corrected Poisson per transport substep** — `Fluid1D/field.py` gains
  `_thomas_general`, `Poisson1DSolver.solve_semi_implicit` (node-CV form
  `−d/dx[(ε₀+h·σ)dφ/dx] = ρ − h·dJ_expl/dx`; reduces exactly to the plain solve at σ=0) and
  the `SemiImplicitField` per-step handle; `DriftDiffusionPlasmaModel.advance(..., field=…)`
  re-solves the field each substep (drift/diffusion CFL only — no ε₀/σ clause, ~10³× cheaper
  post-breakdown; frozen-field path preserved bit-for-bit for `field=None`). The orchestrator
  hands the handle to the real-transport path and diagnostics see the final-substep field.
- **Chemistry energy debit + adaptive cadence** — the numba stiff integrator now accumulates
  each cell's actual inelastic spend Σ_r ε_th,r·∫R_r dt (trapezoidal;
  `FluidIntegrator.last_energy_cost_ev_m3`), which
  `DriftDiffusionPlasmaModel.debit_electron_energy` deducts from the electron energy density
  over the post-chemistry electron count; the chemistry call cadence tightens from
  `fluid_every` steps to every transport step once ν·interval > 0.25 e-foldings
  (`Simulation._chemistry_due`), with the optics advanced on the true elapsed interval.
- **arf.toml R26 corrected to the physical products** `Ar+ + Ar + Ne → Ar2+ + Ne` (the printed
  `Ne+ + 2Ar` is endothermic by +5.8 eV — a table typo; documented deviation at the block,
  verification against Ma et al. PSST 33, 075012 (2024) pending). Xe⁺ orphan noted as the
  paper's own omission.
- Tests 246 → 260: `TestSemiImplicitPoisson` (Thomas-vs-dense, σ=0 reduction, low-σ limit,
  uniform-σ invariance, two-region current continuity, the J_expl lemma),
  `TestSemiImplicitFieldCoupling` (the nₑ=10²⁰ scenario that used to explode now plateaus at a
  ~7×10⁷ V/m cathode fall; pre-breakdown parity with the frozen-field path),
  `tests/test_energy_debit.py` (kernel spend vs augmented Radau, debit semantics, cadence).

## 2026-07-30 — Breakdown unlocked: real EEDF coefficients + LMEA closure + energy-limited chemistry (A1–A8)

The audit-driven series that lets the discharge actually break down (the stub's zero-threshold
attachment made avalanche impossible by construction) and keeps the avalanche paid for:

- **A1 `Fluid1D/eedf_solver.py` `MaxwellianEedfSolver`** — electron rate + transport coefficients
  as Maxwellian-EEDF moments of the measured SIGLO cross-sections (`data/cross_sections/`),
  mean-energy-keyed `CoefficientTable` (native LMEA), Einstein-exact D/μ = Tₑ, derived E/N from
  power balance. 20 dedicated tests (analytic closed forms, benchmarks, the ionization/attachment
  crossover near ⟨ε⟩≈3.7 eV and its F₂ dependence — Chen's central effect, emergent).
  `BoltzmannSpec(backend="eedf")`; the ArF example now runs on real coefficients.
- **A2** SIGLO elastic/momentum-transfer exports for Ne/Ar/Xe/F₂ + F₂ vibrational v1–4
  (loss-only channels), + Ne*/Ar* stepwise-ionization sets (R4/R8 now resolve);
  `reaction_kind_from_metadata` maps `effective → ELASTIC`.
- **A3** TE_DEPENDENT reactions (R9/R13/R14/R18/R19) evaluated from their **published**
  `K0·Tₑ^p·exp(−E0/Tₑ)` per cell (`FluidChemistry.te_analytic`, `dynamic_rates_percell`) instead
  of being wrongly sourced from the LUT. **A4** R13/R14 prefactors read as cm³/s (÷1e6),
  documented in `arf.toml` (as printed they exceed the gas-kinetic ceiling ×10⁴).
- **A5** conduction-current sign: the 1D flux integral now returns the **gap-loading** current in
  the circuit's node-B convention (−A·⟨J_x⟩); 1D and 0D closures agree, gap dissipative.
- **A6** LMEA local-field closure: ⟨ε⟩ recovered as `w·(n_ε/nₑe) + (1−w)·⟨ε⟩_eq(E/N)` with
  `w = nₑ/(nₑ+n_ref)` — repairs the nₑ→0 degeneracy; `CoefficientTable.equilibrium_mean_energy`
  inverts the table's E/N curve; ⟨ε⟩ clamped to the table's resolved range.
- **A7** semi-implicit **gap** coupling in `MpcCircuit.advance` (Chen circuit-analysis §7):
  `I_cond(V_b) ≈ G·V_b + i_off` with G = A⟨σ⟩/d from the plasma; node B rides the collapse where
  the frozen-current form diverges; substeps resolve τ = C_g/G.
- **A8** energy-limited chemistry: `energy_throttle` caps inelastic electron-impact rates at the
  available electron power; adaptive chemistry sub-cycling sized by the true electron-multiplying
  progress frequency (`electron_multiplying_mask` — includes the analytic R9/R18);
  `_N_E_FLOOR_M3` kills a sub-denormal 0/0.

246 tests green. Known remaining defect at this commit (diagnosed, fix planned): the staggered
explicit field↔charge coupling is unstable for nₑ ≳ 3×10¹⁹ m⁻³ at dt=50 ps — the discharge still
runs away and NaNs shortly after breakdown; see `bug-fixes.md`.

## 2026-07-30 — Unified 0D excimer gain/output optical model (Chen §II.C)

Replaced the KrF-only two-stream (Rigrod) photon march with Chen et al.'s lumped-gain laser
module, valid for both KrF* (248 nm) and ArF* (193 nm) — the ArF example now produces a real
193 nm laser pulse. The optical (gain) axis is transverse to the discharge axis, so the model is
0D optically.

- `optics/photon_density.py` `PhotonDensity`: cavity photon density n_ph via Chen Eq. 7, written
  as the linear ODE `dn_ph/dt = a·n_ph + b` and advanced by the exact **exponential-integrator**
  step `n_ph⁺ = n_ph·e^{aΔt} + b·(e^{aΔt}−1)/a` (the closed-quadrature/integrating-factor solution,
  stable against stiffness/overflow). Gain term carries the photon speed c; cavity loss from Eq. 8.
- `optics/gain.py` `GainCavity`: wraps `PhotonDensity` + the output layer — gain `g` (Eq. 14),
  output-coupling efficiency η_out (Eq. 18, with the correct g→0 limit (1−R_OC)/(1+R_OC)), power
  `P_out = n_ph·c·hν·S·η_out` (Eq. 10); depletes the excited states (spontaneous + stimulated) and
  creates products. Absorbers come from the network's PHOTO reactions ∪ `Cavity.absorption`.
- `reactions.network`: new `stimulated_channel` + `photo_absorbers` accessors; `data.py` PHOTO rows
  accept a bare `sigma_m2`. `arf.toml` R38 retagged `stimulated`, R40/R41 `photo`.
- `Cavity` gains `gain_length_m` (l_g), `cavity_length_m` (l_c), `beam_area_m2` (S),
  `axial_se_fraction` (K), `loss_coeff_per_f2` (k_c ⇒ L_op = k_c·x_F2, Chen Eq. 9); R_OC = `r_out`.
- Orchestrator rewired to `GainCavity`; `gain_coefficient()`/`species_total(photon)` now generic
  (no KrF hardcode). The retired two-stream implementation was subsequently removed.
- `examples/example-arf.py` gains OpticalGain/LaserOutput/Efficiency + R_OC=8 %; full MPC run shows
  ArF* avalanche → 193 nm photon burst (n_ph rises ~8 decades) → non-zero output.
- Tests (217): `PhotonDensity` (analytic + ODE agreement, gain/loss), `GainCavity` (Eqs. 14/18/10,
  depletion conservation), network accessors; two-stream test retired.

## 2026-07-29 — Real 1D drift-diffusion-reaction plasma model + energy PDE + LoKI-B integration

Replaced the charged-species stub with the real 1D transport of Chen et al. (2025) §II.B, and
built the LoKI-B coefficient pipeline it consumes.

- `DriftDiffusionPlasmaModel` (`Fluid1D/plasma_model.py`) + `Fluid1D/transport.py`: continuity
  `∂n/∂t+∂ₓΓ=S` with the drift-diffusion flux `Γ=−D∂ₓn+zμnE` (Eqs. 1–2) via **Scharfetter–Gummel**
  exponential fitting; the **electron energy PDE** (Eqs. 3–4: field heating − elastic − inelastic
  losses); electrode BCs (drift-limited ion outflux + cathode secondary emission γ). Internally
  sub-stepped (drift/diffusion CFL + dielectric relaxation). Returns the **exact conduction current**
  `A·⟨J_cond⟩=A·(1/d)∫Σ q_kΓ_k dx` to the circuit — the real node-B load (replaces the 0D closure).
- **LMEA closure** (Hagelaar & Pitchford 2005): `CoefficientTable` gains mean-energy keying +
  `to_mean_energy_key()` + the energy-transport pair μ_ε/D_ε (default (5/3)·electron pair). The
  combined chemistry integrates with **per-cell** electron-impact rates at the local ⟨ε⟩
  (`FluidChemistry.dynamic_rates_percell`, `FluidIntegrator.advance(k_percell=…)`).
- **Ion transport** (neither LoKI-B nor BOLSIG+ produces it): `data/ion_transport/ion_mobility.toml`
  + `reactions.ion_transport` (constant reduced mobility → μ_i·N/N, Einstein D_i, Blanc's-law
  mixtures), every K₀ cited in `sources/kinetics/ion-mobility-data.md` (Ellis/Viehland compilations).
- **LoKI-B**: vendored at `extern/loki-b/`; `LxCatAdapter` (`reactions.lxcat_adapter`) converts the
  project's Kr/F₂ JSON+CSV cross-sections to LoKI-B LXCat `.txt`; `LokiBSolver` defaults to the
  vendored binary + adapter. (The vendored binary currently errors on its own example — a C++-port
  setup-schema issue — so the stub backend stays the validated default; the model is solver-agnostic.)
- Orchestrator: selects `DriftDiffusionPlasmaModel` with the MPC circuit, stops pinning charged
  densities, feeds back the real conduction current. New `ChargedDensityXT` diagnostic exports the
  spatio-temporal charge profiles `n_s(x,t)` (`(n_times, n_cells)` per species) to `diagnostics.h5`.
- Tests (220 pass): SG flux limits/conservation, ion-mobility loader (Einstein/Blanc), LXCat block
  format, DriftDiffusion (drift direction, conduction feedback, energy bounds), mean-energy re-keying,
  real-transport end-to-end + archive.
- **ArF gas system** (toward the Chen replication): `data/reactions/arf.toml` — the complete,
  faithfully-transcribed + per-value-cited Chen Table III (R1–R41); Ne/Xe/F₂ electron-impact
  cross-sections added to `data/cross_sections/`; resolved species `Ne**`, `F2(a3)`, `F2(A1)`, `F2+`
  registered in `particles/species.py`. Loads + compiles (41 reactions). Argon cross-sections
  (R5–R8) are the one data gap (absent from the local LXCat mirror), flagged in the file.
- Full example with the MPC driver now reproduces a genuine avalanche/breakdown pulse (Tₑ→2.4 eV,
  KrF*→1.4e11, 8-decade photon burst) — the phase-2 behavior the 0D closure could not produce.

## 2026-07-29 — Unit-consistency audit + `--volume-units` (cm|m)

Audited every species-density + coupled-unit touchpoint across the program (6-way parallel
audit). Finding: the engine is uniformly **SI (m⁻³) internally and in all coupling constants**
(rate constants m³/s & m⁶/s, cross-sections m², Poisson ρ, reduced field E/N Td); the only
discrepancy was the **output layer** — `SpeciesDensities`/progress reported cm⁻³ while `FieldLine`
reported m⁻³ (same quantity, 10⁶ apart).

Added a global volume/length unit system so users can work entirely in cm or m:

- `Units` spec (`Fluid1D/components.py`): ``volume ∈ {"cm","m"}`` with a single length scale +
  conversion helpers (length/area/density/inverse-length, user↔SI) and labels.
- `Simulation(volume_units="m")` interprets the user's mesh geometry + seed densities in the
  chosen units (→ SI at build) and reports densities/lengths/gain/absorption back in them.
  Internal computation stays SI; **rate constants & cross-sections are untouched** (they carry
  their own units in the data files). Neutral fill (P/k_BT) is intrinsically SI, unaffected.
- Diagnostics (`SpeciesDensities`, `FieldLine`, `OpticalGain`, `Absorption`) + progress now report
  in `sim.units` with matching axis/archive labels — resolving the FieldLine discrepancy.
- Example `01_krf_excimer.py` gains `--volume-units {cm,m}` (default cm); geometry/densities are
  expressed via the paper's SI values so physics is unit-invariant.
- Tests (194 pass): `Units` conversions/labels/round-trips; cm and m inputs → identical SI state;
  density reported in the chosen unit.

## 2026-07-29 — Electron pre-ionization profiles (1D charged-species seed)

Added spatially-resolved t=0 electron pre-ionization for the 1D charged fluid (replacing the
flat seed).

- `Preionization` spec (`Fluid1D/components.py`): `mode="uniform"` or `mode="tanh"`, the latter
  `n_{e,0}(x) = amplitude·tanh(steepness·x) + floor` reproducing Chen et al. Eq. 6
  (`1e14·tanh(100·x) + 1e13` m⁻³ by default), rising from the anode (x=0) toward the cathode.
- `ChargedFluid.seed_profile(name, values)` seeds a per-cell density vector; `Simulation` takes a
  `preionization=` spec and seeds electrons + the balancing ion quasi-neutrally with n_{e,0}(x).
- Example `01_krf_excimer.py` gains `--preionization {tanh,uniform}` (tanh default) and applies it.
- Tests: profile matches the paper formula + endpoints/monotonicity; `seed_profile` shape check;
  end-to-end 1D quasi-neutral seed. 187 tests pass.

## 2026-07-29 — MPC discharge-circuit driver

Added the magnetic-pulse-compression (MPC) discharge circuit that pumps the plasma, replicating
the circuit module of Chen et al., *Phys. Plasmas* **32**, 113502 (2025) (reference + ODE
derivation in `sources/ArF-model/`). Gas-invariant, solved concurrently with the plasma.

- `MpcCircuit` (`src/Fluid1D/field.py`): 7-state RK4 ODE [V_C0, I_m, V_p, I_s, V_b, I_f, V_Cf] —
  C₀→L_m→C_p→L_s→gap with an R_f–L_f–C_f snubber; gap capacitance C_g=ε₀A/d at the cathode node
  (V_b a proper state, no DAE). Internally sub-stepped for stiffness; all component values
  configurable (defaults = paper's). Cross-checked ODE signs via an independent nodal+mesh
  derivation workflow.
- `CircuitSpec` extended with backend `"mpc"` + component values (C0, Lm, Cp, Ls, Rf, Lf, Cf, V0,
  bus R1/Rs, C_head, max_substep).
- Wired into `Simulation`: circuit → cathode Dirichlet BC for the 1D Poisson solve; plasma →
  conduction current back to the circuit (0D closure `G·V_b`, `conduction_current_0d`, pending the
  real 1D flux integral). Example `01_krf_excimer.py` gains `--circuit {mpc,stub}` (mpc default).
- Tests: complete phase-1 C₀→C_p resonant transfer (peak V_p→V₀ at ~111 ns), lossless energy
  conservation, gap-capacitance geometry, node-B conduction loading, coupled smoke. 180 tests pass.
- Validated: MPC phase-1 voltage compression reproduced (peak |V_b|≈17.9 kV at 110 ns). Phase-2
  breakdown awaits the real 1D charged-species transport (avalanche vs. F₂ attachment).
