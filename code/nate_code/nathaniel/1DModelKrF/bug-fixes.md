# Known bugs & fixes

Track each non-trivial bug and the implementation that fixed it. Newest first. Use this to avoid
re-introducing regressions and to document subtle physics/numerics traps.

Format:
```
## <short bug title>  (status: open | fixed YYYY-MM-DD)
- **Symptom:** what was observed.
- **Cause:** the root cause.
- **Fix:** the change that resolved it (file/function), or the current workaround if still open.
```

---

## Tiny sweep presets supplied two point selections  (status: fixed 2026-08-01)

- **Symptom:** `examples/arf_f2_sweep.py --tiny` and `examples/krf_f2_sweep.py --tiny`
  stopped before launching workers with `values=... is mutually exclusive with start/stop/step`.
- **Cause:** each tiny preset supplied an explicit three-value sequence while the entrypoint also
  forwarded the parser's default range to the sweep API. The stricter ambiguity validation exposed
  the duplicate selection.
- **Fix:** each entrypoint now constructs one selection mapping: explicit values for `--tiny`, or
  `start`/`stop`/`step` for a normal run. Three-point ArF and KrF process-pool smoke runs cover the
  corrected path.

## Source-coupled stage starved of heating: floor-minted energy, no lasing branch  (status: fixed 2026-08-01)

- **Symptom:** the 128-cell/50 ps/0.15% two-term run recovered a 2.24 kA, 0.888 J conducting
  branch but peaked at `8.20e19 m^-3`, emitted only 1.21 microjoules, and booked +0.277 J of
  coefficient-bound energy against 0.288 J of accepted reaction cost; run-mean T_e 1.63 eV.
- **Cause:** the electron stage had NO heating source — it drained a stored buffer while J.E
  lived only in transport. Inelastic cooling is stiff (tau_eps ~ 10-90 ps at 7-12 eV, at/below
  the 25 ps half-steps), so the stage relaxed toward the floor instead of the J.E-balanced
  eps_eq(E/N); the floor clamp then minted the difference. Ionization is exponential across the
  eps_eq knife-edge (nu_net -1.9e7 at 6.0 eV -> +7.7e8 at 7.4 eV), so the bias compounded.
- **Fix:** deliver the transport-ledgered local heating INSIDE the stage's substep loop
  (`_energy_coupled_substeps` gains `source_ev_s`; `advance_energy_coupled` gains
  `power_source_j_m3_s`). Delivery is Strang-consistent (lead half-step: previous interval's
  power, guarded across climax-k h jumps; trail half-step: exact residual — per-interval
  conservation exact) and RECTIFIED per cell (see next entry). Substep sizing uses a physical
  denominator `max(stored, eps_ceiling*max(n_e,1e10))` so trace-density sheath cells cannot pin
  the global substep count. Per-cell identity extends to
  `E_final + E_cost - E_bound_adjustment - E_delivered = E_initial` (tests:
  tests/test_energy_debit.py, tests/test_energy_closure.py). Measured after the fix (256c/50ps
  0.15%): bound adjustments 0.277 J -> +4.5 mJ (0.8% of 0.563 J accepted cost), avalanche rides
  T_e 4.35-4.60 eV, peak n_e `1.6e20 m^-3`, 4.0 kA, 1.45 J deposited, output jumped
  discontinuously from 1 uJ to ~2 mJ.

## Signed sheath field-work booking cancels bulk heating when delivered as a source  (status: fixed 2026-08-01)

- **Symptom:** first in-stage delivery attempt made the books WORSE (+0.482 J bound ledger) with
  `total_field_work_j = -0.153 J` — net-negative electron field work while the gap absorbed
  +0.883 J.
- **Cause:** the cell-centred Γ_e·E booking in the one-cell electrode falls (sub-Debye dx) books
  large phantom NEGATIVE work — measured to cancel the entire positive bulk heating post-
  breakdown (24c probe: cathode cell -266 units vs bulk +274 in 80-100 ns). The old in-transport
  scheme silently rectified this at the per-substep floor; a signed source re-imported it.
- **Fix:** `Simulation._deliverable_heating_j_m3` delivers `max(field_work - elastic, 0)` per
  cell; the discarded phantom-negative part stays visible in `undelivered_field_work_j` and in
  the signed `field_work` books. The supply cap still bounds the positive lobe at accumulation.

## Continuous reconciliation dwarfed accepted chemistry (x1278)  (status: fixed 2026-08-01)

- **Symptom:** the 128-cell diagnostic accumulated a 1.994 J provisional sink but only 1.561 mJ
  accepted reaction cost (+1.993 J signed refund, 34 eV equilibrium residual).
- **Cause:** two-part. The sink's row/column packing is CORRECT (pinned by
  tests/test_transport_kernel.py::TestContinuousSinkBasis against
  `energy_consumption_rate_percell`, rtol 1e-12); the divergence came from the lump-sum refund
  `electron_energy_j_m3 = transport_energy + prepaid`, which the stage pre-clamp destroyed above
  the ceiling in exactly the cells that prepaid the most — the refund became a bounds-ledger
  entry instead of usable energy, so the mode starved and its books diverged.
- **Fix:** the prepaid sink is refunded as an in-stage RATE (`power_source_j_m3_s =
  prepaid/interval`), reusing the source-coupled delivery mechanism; sizing includes the prepaid
  energy. Cross-mode agreement is pinned by
  tests/test_energy_closure.py::TestCrossModeAgreement (n_e within x2, T_e within 2 eV on a
  short kick; the pre-fix ratio was ~x13 collapse). Mode remains diagnostic-only.

## Heavy-partition chemistry landed as 1 ns impulses in source_coupled  (status: fixed 2026-08-01)

- **Symptom:** ion-ion recombination/dimerization (nu ~ 1e9 1/s at conduction ion densities) ran
  strictly every `fluid_every` (1 ns) because `_chemistry_due` returned False for graded
  triggers in source_coupled — the same impulsive-kick class as the old F2-sweep sawtooth.
- **Cause:** heavy_chem's electron masks are all-False by construction (heavy = complement of
  the electron-source predicate), so the old grading had nothing to watch.
- **Fix:** new `charged_changing_mask` on `FluidChemistry.for_plasma` (any charged species with
  nonzero stoichiometry, registry `is_charged`); `_chemistry_due` grades the heavy cadence on
  `progress_frequency` over that mask (`_last_heavy_nu`, kept strictly separate from
  `_last_nu_max` so climax subdivision is not polluted).

## Photoelectrons born at zero energy  (status: fixed 2026-08-01)

- **Symptom:** R40/R41 (and KrF R43/R44) photo feedback diluted the electron pool — a cooling
  bias applied exactly during the pulse, when the eps_eq knife-edge matters most.
- **Cause:** no birth-energy data path from the reaction rows to the optics feedback.
- **Fix:** optional `electron_energy_eV` on PHOTO rows (hv minus binding/detachment; NIST ASD +
  Blondel EA(F) anchored: ArF 3.76/3.02 eV, KrF 3.10/1.60 eV), carried through
  `photo_channels()` (inserted BEFORE `photon_mode` — downstream filters select on
  `channel[-1]`), added to the electron energy density in `_advance_photo_channels`, ledgered
  in `photoelectron_birth_energy_j` (domain-mean J/m^3 at the cavity, converted to joules at
  booking), and the transport-owned `_last_n_eps` is re-synced after every optics call.

## PHOTO absorption omitted declared species feedback  (status: fixed 2026-07-31)

- **Symptom:** R40/R41 removed cavity photons but left Ar**, F-, Ar+, F, and electrons unchanged.
- **Cause:** optics owned photon attenuation while the chemistry integrator deliberately excluded
  PHOTO rows, leaving no layer responsible for their stoichiometry.
- **Fix:** `GainCavity` applies the declared products with the exact exponential frozen-fluence
  hazard. Tests conserve atoms, charge, photons, and existing electron energy. At this checkpoint,
  new photoelectrons were explicitly born at zero energy pending sourced yields.
- **Superseded 2026-08-01:** the later "Photoelectrons born at zero energy" fix above supplies
  and ledgers the active ArF 3.76/3.02 eV and KrF 3.10/1.60 eV birth energies. The zero-energy
  sentence is retained only as the historical behavior of the 2026-07-31 fix.

## Reaction energy was debited independently of accepted chemistry  (status: fixed 2026-07-31)

- **Symptom:** the F2 sweep could charge electrons for a full frozen-rate inelastic sink while a
  later rate throttle or chemistry budget accepted only part of the corresponding ionization or
  excitation. The reverse split also let chemistry advance against a stale mean energy.
- **Cause:** `plasma_model.py`/`transport_kernel.py` owned a continuous reaction sink, while
  `FluidIntegrator` independently owned reaction progress and a second energy constraint.
- **Fix:** transport now owns energy flux, field work, and elastic exchange only.
  `FluidIntegrator.advance_energy_coupled` refreshes electron-impact/analytic coefficients from
  the live mean energy and accepts each implicit reaction trial together with its exact threshold
  debit. A post-reaction energy-floor constraint reduces only costly channels; any unavoidable
  table-bound adjustment is a separate ledger. Its zero-cost endpoint retains threshold-free
  Penning electron production rather than reverting the entire trial. `examples/arf_convergence.py`
  supplies restartable 128/256/512-cell and 100/50/25-ps production gates.
- **Verification:** focused energy/transport/field/simulation suites pass; the per-cell identity
  `E_final + E_cost - E_bound_adjustment = E_initial` closes to floating-point tolerance, and an
  ArF R20 regression pins the Penning endpoint. The full 0.15% matrix is stable but fails every
  spatial/timestep gate; peak `n_e` spans `3.34e18`--`2.22e19 m^-3` and output remains nJ-scale.
  Bound adjustments are microjoules, and disabling the `J.E` cap is near-null, confirming that
  matching Chen is still a mesh/coupling/model-validation issue rather than a reason to relax the
  corrected energy closure.

## Diagnostic plots hid data and imposed source-owned filenames  (status: fixed 2026-07-31)

- **Symptom:** `runs/arf-lee` plots placed peaks at only 2–40% of the vertical extent, symmetric
  optical-gain limits reserved an unused negative half-axis, and the linear species plot hid
  active species under the Ar background. Every collector also selected its own `.png` name and
  x/y composition, preventing an experiment author from splitting or rescaling those series.
- **Cause:** `linear_y_limits` rounded small positive peaks up to a full decade, signed data were
  always made symmetric, and collection/rendering were coupled in each diagnostic's `finalize`
  method. Circuit and reaction nanosecond arrays were also mislabeled `times_ps`; field position
  was exported through the archive's default time key.
- **Fix:** collectors now only record/export arrays. Run files declare mandatory plot filenames,
  x/y `DataRef`s, optional linear/log/symlog scales, and fixed ranges through archive-backed
  `LinePlot`/`HeatmapPlot` objects. Unpinned limits retain every finite plottable value with a
  small margin and fit y only over the visible x window. The mislabeled archive axes were fixed,
  and sweep plot helpers likewise require run-owned names and dataset specifications.

## 2026 numerical audit: bounds, energy, keys, cadences, and coupling  (status: fixed 2026-07-31)
- **Symptom:** bounds-checked transport failed; distinct coefficient grids reused one cache file;
  equivalent E/N and mean-energy tables diverged; prescribed fields heated but did not drift
  species; single fluid intervals integrated to zero energy; short runs skipped chemistry/optics;
  high-gain optics and implicit chemistry exceeded their inventories; node B was advanced twice.
- **Cause:** an `N-1` Thomas workspace, incomplete cache identities, implicit lookup-unit
  assumptions, centre-only field superposition, trapezoids applied to interval averages, no final
  cadence flush, end-state optical depletion, pre-step-only chemistry budgets, and overlapping
  circuit/plasma node subflows.
- **Fix:** see `NUMERICAL_AUDIT.md` for the full defect table and regression anchors. Core changes
  are in `transport_kernel.py`, `boltzmann.py`, `field.py`, `plasma_model.py`, `integrator.py`,
  `photon_density.py`, `gain.py`, `simulation.py`, and `sweep.py`.

## F2+ transport uses an uncited default mobility  (status: open)
- **Symptom:** the Ne/Kr/F2 production sweep warns that `F2+` has no entry in
  `data/ion_transport/ion_mobility.toml`; its drift/diffusion therefore uses the engine's generic
  reduced-mobility fallback. (`ArKr+`, `He+`, and `M+` appear in the same registry-wide warning
  but remain unpopulated in the KrF network.)
- **Cause:** the Boltzmann solvers are electron-only and the current cited swarm table has no
  F2+ measurements in Ne, Kr, or F2. The new R12/R21 pair makes F2+ an active, short-lived ion.
- **Fix:** pending a primary F2+ swarm-mobility source (or a documented collision-integral model)
  for the three background gases. Do not silently promote the generic fallback to a cited value.

## KrF starter network mixed wrong units, products, state models, and citations  (status: fixed 2026-07-31)
- **Symptom:** the starter `krf.toml` had only a small Kr/F2 subset, labeled a two-body
  `2 Kr*` Penning reaction as cm6/s, treated a rate `sigma*c` as a photo cross section, used
  an uncited second Kr excitation, omitted all Ne-buffer paths, and left several excited and
  charged products without physical sinks. The proposed `Kr2+ + F- -> KrF* + F` route also
  violated both Kr and F atom balance, while the copied `KrF* + F2 -> Kr + 3F` products belonged
  to an ArF dissociative branch rather than the measured KrF electronic quench.
- **Cause:** the file was explicitly a schema-demonstration starter assembled from secondary
  Cooper-Watson/Chow rows; cgs reaction order, B/C state reduction, and primary DOI provenance
  had not been audited as a coupled network.
- **Fix:** replaced it with the 45-row, cited Ne/Kr/F2 reduction in `data/reactions/krf.toml`.
  Rates are tagged as measured/model/transfer in comments and references; every row conserves
  atoms and charge; molecular intermediates have radiative or collisional exits; F2+, Ne2+,
  Kr2+, Ne**, and Kr** have explicit losses. The KrF(B) lifetime/gain-cross-section pair is
  internally consistent, and network photo absorbers replace (rather than duplicate) the
  empirical cavity F2 loss in the KrF examples.
- **Verification:** 45 rows compile as 18 fluid + 23 kinetic + 4 photon reactions with no orphan
  charge warning; chemistry atom-conservation, Boltzmann, optics/simulation, and photon-channel
  tests pass. The 10-point real two-term/MPC sweep completed with all point archives and seven
  aggregate figures under `runs/krf_f2_sweep/`.

## R39 printed products are E1-forbidden — Ar** radiative exit read as the 4p→4s cascade  (status: fixed 2026-07-31)
- **Symptom:** with the R5/R6 branching corrected, ~all Ar**(4p) inventory still exited as
  95 nm photons (printed R39: `Ar** → Ar + hv95`, 10⁹ s⁻¹) — a one-way drain that fed
  neither the harpoon pool nor stepwise ionization.
- **Cause:** the printed transition is parity-forbidden (3p⁵4p is even, like the 3p⁶ ground
  state — E1 requires a parity change); the 95 nm resonance line belongs to Ar*(4s ¹P₁)→ground,
  and the fast 4p decay is the ~810 nm 4p→4s cascade INTO Ar*. Same typo class as R26
  (endothermic printed products). Akashi 1995 (Chen ref 30) models the cascade explicitly.
- **Fix:** R39 products read as `Ar* + hv` (documented deviation in `arf.toml`; Chen's printed
  A = 10⁹ s⁻¹ kept). The radiative-channel driver credits non-photon products, so the 4p flux
  now lands on Ar*.
- **Verification (200 ns two-term):** gain 2.91 → 3.62 m⁻¹, ArF* 6.0×10¹⁹ → 7.3×10¹⁹,
  n_e → 3.4×10²⁰, E_out 106 µJ → 2.6 mJ — every observable moves toward Chen. Suite green.

## F2-sweep sawtooth — impulsive chemistry-call kicks at the 1 ns cadence  (status: fixed 2026-07-31)
- **Symptom (user-reported):** every F2-sweep overlay trace (n_e, T_e, E/N, P_dis) — and V_b
  itself — jagged/saw-toothed vs Chen's smooth Fig. 7/12 curves; survived the dilution-only
  debit fix.
- **Cause (measured):** FFT + phase-profiling of the raw h5 series put the dominant period at
  exactly 1.00 ns (the `fluid_every` chemistry cadence; the physical L-C ring is 1.4-2.5 ns)
  with the whole sample-to-sample jump concentrated at the call phase (x154 the median step).
  Two roots: the cadence trigger watched only electron-MULTIPLYING channels, so the conducting
  phase (ionization slow; thresholdless F2 attachment ~1e8/s still churning) stayed on the
  1 ns cadence while each call applied the accumulated electron loss as one impulsive kick to
  n_e -> sigma -> V_b/E; and the binary 20-or-1 trigger's 0.25 e-fold threshold sat above the
  conduction-phase nu ~ 0.1 e-folds/ns.
- **Fix:** `electron_changing_mask` (net electron production OR consumption) now drives the
  cadence stiffness estimate, and the cadence is graded — chemistry runs every
  `ceil(0.05/(nu*dt))` steps (`_chem_efolds_per_call`), from every step at breakdown to the
  nominal `fluid_every` when quiet.
- **Verification (200 ns two-term):** 1 ns phase-lock ratio n_e 153.8 -> 2.0, T_e 11.0 -> 1.5,
  E/N 6.4 -> 1.4 (~1 = uniform); residual V_b modes 1.4-2.5 ns = the physical ring. +33% wall.
- **Honest consequence:** the sawtooth's impulsive reheating was inflating excimer yield
  (pulse energy drops with the smooth chemistry) — energy numbers from sawtooth-era runs carry
  that bias; re-baseline via a sweep rerun.

## Odd–even checkerboard — discrete Gauss law blind to cell-alternating charge  (status: fixed 2026-07-31)
- **Symptom (user-reported):** `field_line.png` shows cell-alternating n_e/Ar⁺ (orders of
  magnitude) and a ρ(x) zigzag (±600 C/m³) while |E|(x) stays smooth. Audit (x–t maps):
  onset at avalanche start (80–82 ns), order-unity through the WHOLE conducting/laser
  window — not just afterglow cosmetics; the F⁻/ArF\* split carried a mode-induced bias.
- **Cause (measured):** the node-CV charge gather `ρ_node = 0.5(ρ_i + ρ_{i+1})` in the
  Poisson solves (and the kernel's in-loop Thomas block) annihilates the cell-alternating
  charge mode EXACTLY (E response 0.000 vs 8.9×10⁷ V/m for a smooth ρ of equal magnitude):
  charge separation between adjacent cells drew no restoring field, so the avalanche's
  differential gain amplified any seed unopposed. Two earlier hypotheses disproven: the
  debit double-count only seeded it, and SG face-velocity averaging damps (not pumps) it.
- **Fix:** cell-centred φ formulation throughout — unknowns at cell centres, Dirichlet walls
  via dx/2 half-cells, arithmetic face coefficients, E native at the N+1 faces
  (`last_e_faces`), and the SG drift terms + boundary fluxes + energy-PDE drift consume the
  FACE field directly (any cell-averaging of the face field re-blinds the operator — the
  restoring flux needs face-native E). `field.py` solve/solve_semi_implicit, the kernel's
  in-loop solve, and the numpy coupled path, in parity; J_expl now passes as face currents.
  Standalone validation: a seeded 2 % charge checkerboard neutralizes ×140 in 7.5 ns (node
  scheme: untouched). Regression: `test_alternating_charge_draws_restoring_field`.
- **Re-baselines:** wall half-cell resolves a ~2× deeper cathode fall forming over ~150
  steps (was ~40) — sheath-test bounds re-pinned; expect small shifts in breakdown numbers.

## T_e sawtooth at the chemistry cadence — inelastic energy debited twice  (status: fixed 2026-07-31)
- **Symptom (user-reported):** ⟨ε⟩ stable at ~6 eV (T_e 3–4 eV) collapses to the 0.2 eV floor
  **domain-wide** at exactly the chemistry cadence (every `fluid_every·dt` = 1 ns), re-heating
  within ~0.4 ns — a sawtooth the T_e diagnostics sample as "drops to ~0.18 eV then spikes back".
- **Cause:** double subtraction. The transport energy PDE already sinks inelastic losses
  continuously (the ELECTRON_IMPACT eth-channel table, per substep, in both kernel and numpy
  paths); `debit_electron_energy` then subtracted the *same interval's* chemistry spend again
  at each chemistry call. Pre-breakdown the interval's inelastic loss ≈ the whole stored n_ε
  (heating pays for it continuously), so the second subtraction floored every cell. Its
  non-uniform per-cell kick also seeded cell-to-cell density noise (see field_line entry).
- **Fix:** both transport paths accumulate the PDE-sunk inelastic energy per advance
  (`inel_sunk`, J/m³); the simulation sums it between chemistry calls and debits only
  `max(chem_cost − sunk, 0)` — Phase 2's intent (chemistry pays what the PDE under-counted;
  newborn electrons don't inherit free energy) without the double count.
- **Verification:** 70 ns two-term probe — 0 of 720 samples below 0.5 eV in the 30–66 ns
  window (was ~1 in 6); ⟨ε⟩ steady at 5.9 eV.

## Intrinsic efficiency blows up (∞ / 0) — ledger drive was the ± field-heating integral  (status: fixed 2026-07-31)
- **Symptom (user-reported):** `intrinsic efficiency` prints ∞ (or 0) and the cumulative
  efficiency curve is garbage; the checkpoint run's ledger showed `drive ∫P_in dt = −6.7 J`.
- **Cause:** the energetics ledger recorded `report.deposited_power_w` — the discrete J·E
  field-heating integral, which swings ± with the under-resolved sheath work (booked −34 J
  over 200 ns) — so the denominator E_init + ∫P_in could cross zero. E_init is the µJ-scale
  seed-plasma energy, not a stabilizer.
- **Fix:** with real transport the ledger records the conduction power |I_cond·V_gap| — the
  quantity whose books close against the circuit (gap_j = 1.607 J vs delivered 1.599 J).
- **Verification:** 200 ns run — `∫P_in dt = +1.631 J`, efficiency finite.

## Laser output ×75 energy mint — mode volume vs discharge volume (plan B2)  (status: fixed 2026-07-31)
- **Symptom:** once the sawtooth fix let the cavity lase properly, E_out = 1.452 J from a
  1.63 J deposition (η = 89 %) — needing ~10¹⁸ ArF\* de-excitations against a ~10¹⁶ inventory.
- **Cause:** `beam_area` defaulted to the **discharge area** A = 90 cm², making the mode
  volume l_g·S = 75× the real discharge volume. Photon pumping treated the domain-mean n_up
  as filling that fictional volume while stimulated depletion drew on the real one — a ×75
  photon/energy mint, invisible below threshold (mJ scale), joules once lasing.
- **Fix:** default beam face S = d·A/l_g (the optical axis runs along l_g, so S = gap ×
  discharge width and l_g·S ≡ the discharge volume — the plan's deferred B2 geometry item,
  now required). `Cavity(beam_area_m2=...)` still overrides explicitly.

## Circuit gap load annihilates the bank — bulk-σ conductance vs sheath-limited gap  (status: fixed 2026-07-31)
- **Symptom:** with the heating mint closed and the discharge circuit-clamped, the 200 ns
  two-term gate deposited only 0.154 J into the plasma against Chen's ≈1.36 J, while the
  1.62 J bank still drained (1.4 J unaccounted). Downstream: ArF\* peaked at 2.0×10¹⁹ m⁻³ —
  below the 3.6×10¹⁹ lasing threshold — F⁻ ×3 low, no laser.
- **Cause (measured from the gate series):** the circuit's passive gap load used
  `G = A·σ_mean/d` — the **bulk** ohmic conductance. During/after collapse the bulk is hugely
  conducting but the gap current is sheath/electrode-limited (I_cond ≪ G·V_b), so the circuit
  dissipated ∫G·V_b² dt = 1.64 J into a fictional resistor the plasma never saw
  (∫I_cond·V_b dt = 0.154 J). Passivity (the fix for the earlier i_off active-source pump)
  was preserved, but the *magnitude* of the load was an order too large exactly when it
  mattered.
- **Fix:** energy-matched secant conductance. Both coupled paths accumulate the actual gap
  absorption ∫i_gap·v_b dt (circuit convention, dissipative ≥ 0) and ∫v_b² dt per advance;
  the report hands the circuit `G = max(∫i·v, 0)/∫v²` — still passive, and its circuit-side
  dissipation over the next interval equals the plasma's true absorption by construction
  (one-interval lag, bounded by the climax subdivision). Kernel + numpy coupled paths in
  parity; the frozen path keeps the ohmic form.
- **Verification (200 ns two-term gate):** books close end-to-end — circuit→node +1.599 J,
  plasma absorbed +1.607 J, bank left 0.001 J (Chen: 1.36 of 1.620 J deposited). ArF\* peak
  7.3×10¹⁹ m⁻³, **above threshold for 30 ns**; F⁻ 3.8×10²⁰ (Chen 2×10²⁰); mean-density n_e
  peak 5.5×10²⁰ (Chen ~10²¹); T_e(mean) peak 5.28 eV (Chen ~5 eV); collapse from −18.7 kV
  (Chen −16…−19 kV) at ~89 ns (Chen 100–110 ns — remaining gap, see checkpoint notes).

## Collapse-window heating mint — discrete J·E work in sub-Debye sheaths  (status: fixed 2026-07-30)
- **Symptom:** with every licence fix in place (hard budget, network cap, boundary-thermal CFL
  guard), the collapse still ran away: at ~86.5 ns (V_b collapsed to −40…−600 V, nₑ crossing
  8×10²², sheath |E| → 2.3×10⁹ V/m) the recorded field-heating power exploded to 10¹¹–5×10¹³ W
  against an honest supply |i_s·V_b| of 10⁵–10⁶ W, and the cumulative deposition swung
  ±1.9×10⁷ J — funding ionization to nₑ ~ 10²⁶.
- **Cause (measured on the 84–91 ns probe window):** discrete J·E work is not conservative in
  cells whose Debye length is far below Δx = 47 µm (λ_D ~ nm at nₑ > 10²²): the space-charge
  field exchange oscillates with enormous ± lobes that nearly cancel in the net. Per-cell
  *rectification* converts that numerical noise into real licence: `max(p_field, 0)` feeds the
  spendable-power licence, and the energy PDE accumulates the positive lobes into n_ε, which
  the global cap then counts as stored "spendable" energy.
- **Fix:** rectified-heating cap in the kernel (and the numpy coupled path, for parity): per
  substep, compute the rectified totals Σmax(p,0) and Σmax(−p,0) over the domain; if either
  exceeds `2·|I_gap·V_b| + 1 kW` — the gap's true instantaneous electrical power, which bounds
  ∫J·E dV in 1D electrostatics — scale the whole p_field profile down to that bound. Exact
  pre-breakdown (∫J·E dV ≈ I·V there, ratio ~1, cap never engages — verified 83.5–86.1 ns);
  kills the mint by construction at the collapse.
- **Verification (95 ns two-term probe):** stable throughout, nₑ clamped at 7.9×10²¹ (Chen's
  circuit-clamped scale; was 10²⁶), first collapse −17.8 kV @ 84.4 ns (Chen −16…−19 kV band),
  I_s ramping through −2.9 kA, ionization licence spend 0.020 J on the 1.62 J bank, positive
  heating 26 mJ vs 136 mJ delivered. Wall time 106 s for 1900 steps (was 9–20 min).
- **Residual (known artifact, conservative direction):** the deposition *accumulator* still
  books ~−1.5 J of net-negative sheath work across the post-breakdown ringing (starts
  ~88 ns) — the negative lobes are capped at the supply scale, not zero. It cannot license
  ionization (only positive lobes do, and they are capped + n_ε is floor-clipped), it
  under-heats wall cells slightly, and it would vanish with a λ_D-resolved sheath — the same
  documented next refinement as the semi-implicit coupling note.

## Vendored LoKI-B binary errors on its own shipped example  (status: superseded 2026-07-31)
- **Symptom:** `extern/loki-b/build/app/loki <setup.json>` prints `[json.exception.type_error.304]
  cannot use at() with number` and produces no output — including on the *shipped* `default_lokib_setup.json`
  (N₂), so it is not specific to our generated Kr/F₂ setup.
- **Cause:** a setup-schema mismatch in the C++ LoKI-B port's JSON reader (a field it accesses with
  `.at()` is a scalar in the input) — a LoKI-B tool/build issue, not our integration.
- **Fix:** the external adapter/vendor path was deleted. Production coefficients now come from the
  in-tree `TwoTermBoltzmannSolver`; `MaxwellianEedfSolver` is the measured-data approximation and
  `AnalyticStubSolver` is test-only. No external LoKI-B executable is required.

## Energy-ledger leaks past breakdown — the full chain  (status: fixed 2026-07-30)
- **Symptom:** with the field-coupling and chemistry fixes in place, runs still reached
  near-full ionization (nₑ → 7e25) by ~62 ns with the ledger violated by 3–9 orders
  (kJ–GJ against a 1.62 J bank), books diverging at the avalanche climax.
- **Causes (each measured on an instrumented probe, fixed in sequence):**
  1. circuit RK4 sub-stepping ∝ G — hours of wall time resolved honestly, or V_b → ∞ when
     capped → replaced by the exact augmented-matrix-exponential advance;
  2. the reported-current linearisation offset `i_off` is an ACTIVE source under staleness —
     a 1.62 J bank "delivered" 16 kJ → gap handed to the circuit as a PASSIVE G·V_b load;
  3. the plasma's pinned Dirichlet V_b is an infinite-impedance source between exchanges
     (0.085 J drawn on a 91 µJ licence in one call) → kernel co-advances a local node-B copy
     (network drive + own measured draw) per substep; exchange subdivides ≤256× at climax;
  4. the A8 throttle is a RATE cap — in-call growth outran the licence ×930 → hard per-cell
     energy BUDGET inside the integrator kernel (costly channels stop when spent);
  5. (regression, reverted same day) relaxing chemistry substeps to 2 e-folds overshot the
     avalanche catastrophically — 0.25 e-foldings restored.
- **Verification:** 65 ns eedf probe: circuit→gap +8.5 mJ vs plasma heating +7.1 mJ
  (ledger agreement ~20 %), ionization cost 1.2 mJ ≤ delivered, nₑ a paid-for avalanche
  (1.5e17 → 1.9e19 by 64 ns, V_b −10 kV still standing), 60 s wall for 1300 steps.

## Explicit field↔charge coupling unstable past breakdown — manufactures energy, then NaN  (status: fixed 2026-07-30)
- **Symptom:** the ArF run breaks down, then runs away — nₑ climbs decades past the
  circuit-clamped ~10²¹ m⁻³, ⟨ε⟩ pins at the coefficient-table ceiling (30 eV) domain-wide,
  V_b rings to ±MV, and `charged.density` goes NaN one step later (measured at t = 68.1 ns with
  the example defaults; 53.25 ns in an earlier configuration). "Peak nₑ" numbers from such runs
  (e.g. 3.2×10²³) are snapshots of the runaway in flight, not physics.
- **Cause (measured, not inferred):** Poisson is solved once per outer Δt = 50 ps with the field
  then **frozen** across all transport substeps; the internal ε₀/σ substep guard stabilizes
  transport against the frozen field but not the field↔charge feedback loop, which is linearly
  unstable once σ > ~2ε₀/Δt, i.e. nₑ ≳ 3.4×10¹⁹ m⁻³ (measured onset 2.8×10¹⁹). The oscillating
  field *creates energy*: 1.45×10⁹ J of Joule heating against the 1.62 J capacitor bank, the MPC
  circuit pumped from 1.61 → 379 J by garbage ±10⁵ A conduction currents, one 1 ns chemistry call
  then ionizing 99.6 % of the domain's neon. Isolation repro (chemistry off, V_b frozen −5 kV,
  nₑ=10²⁰): Δt = 50 ps explodes ×30 000 in 1 ns; Δt = 5 ps stays physical. Amplifier: chemistry
  k(⟨ε⟩) frozen over the 1 ns fluid interval with a start-of-interval throttle under-counts the
  in-interval growth (×2.3 measured at a healthy ⟨ε⟩=5 eV state; ×15–30 at breakdown states).
  The chemistry integrator itself is accurate (+5.8 % vs Radau at the substep floor) and
  conservative (atom error ≤ 4×10⁻¹⁶ through the crash).
- **Fix (2026-07-30):** two-part.
  (1) **Semi-implicit σ-corrected Poisson** re-solved inside the transport substep loop:
  ``−d/dx[(ε₀+h·σ)dφ/dx] = ρ − h·dJ_expl/dx`` (`Poisson1DSolver.solve_semi_implicit` +
  `SemiImplicitField` handle in `Fluid1D/field.py`; `DriftDiffusionPlasmaModel.advance(...,
  field=handle)` re-solves per substep, sized by the drift/diffusion CFL alone — the ε₀/σ
  substep clause is gone, so the collapsed phase integrates in a handful of substeps instead
  of the 20 000-cap). The frozen-field path (`field=None`) is preserved bit-for-bit for the
  stub/test surface. Regression: `TestSemiImplicitFieldCoupling` — the nₑ=10²⁰/−5 kV scenario
  that blew to 10¹⁰ V/m now builds a ~7×10⁷ V/m cathode fall that plateaus and relaxes.
  (2) **Chemistry pays its energy bill**: the numba integrator accumulates each cell's actual
  inelastic spend Σ_r ε_th,r·∫R_r dt (trapezoidal; `FluidIntegrator.last_energy_cost_ev_m3`),
  and `DriftDiffusionPlasmaModel.debit_electron_energy` deducts it from the electron energy
  density over the **post-chemistry** electron count (no free ⟨ε⟩ for newly created
  electrons); the chemistry call cadence tightens from `fluid_every` to every transport step
  once ν·interval > 0.25 e-foldings (`Simulation._chemistry_due`), so k(⟨ε⟩) staleness during
  the collapse matches the energy PDE's own cadence. Tests: `tests/test_energy_debit.py`
  (kernel vs augmented-Radau reference, debit semantics, cadence trigger).
