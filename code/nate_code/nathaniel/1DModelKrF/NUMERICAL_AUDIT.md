# Numerical audit

Audit date: 2026-07-31; current-contract update: 2026-08-01. Scope: coefficient construction and caching, cross-section interpolation,
charged transport, Poisson/circuit coupling, stiff chemistry, optical gain, sweep orchestration,
diagnostic integration, and configuration validation.

The audit combined source inspection with targeted reproductions, NumPy/Numba parity checks,
closed-form comparisons, conservation bounds, and the repository test suite. It is not a claim of
experimental validation or mesh/timestep convergence for every production sweep.

## Corrected defects

| Severity | Defect | Correction and regression |
|---|---|---|
| Critical | The Numba Thomas solve allocated its forward RHS workspace with `N-1` entries but wrote index `N-1`. Normal Numba execution could silently corrupt memory. | Allocate `N`; run the transport-kernel suite with `NUMBA_BOUNDSCHECK=1`. |
| Critical | Transport removed the full frozen-rate inelastic sink before chemistry, while chemistry could independently reduce the corresponding reaction progress. Electrons could therefore pay for reactions that never occurred, or reaction progress could use a stale energy. | Remove reaction losses from transport. Co-integrate live-energy coefficient refresh, reaction progress, and its threshold-energy debit inside every implicit chemistry substep. Costly trials are reduced against the post-reaction electron-energy floor; bound adjustments are ledgered separately. The zero-cost limiter endpoint explicitly preserves Penning progress. Per-cell closure, dynamic-cooling, and floor-Penning tests pin the contract. |
| Critical | The co-integrated stage had no heating source: it drained a stored buffer while J.E lived in transport. At stiff inelastic timescales (tau_eps 10-90 ps at 7-12 eV) the stage relaxed to the floor instead of eps_eq(E/N); the floor minted +0.277 J and the discharge never ignited past `8e19 m^-3`/microjoules. The signed sheath field-work booking (phantom-negative in the one-cell electrode falls) cancelled bulk heating when delivered naively. | Deliver the rectified (positive-part) field-work-minus-elastic ledger inside the stage as a per-substep power source, Strang-split with exact per-interval conservation; `delivered_source_j`/`undelivered_field_work_j` ledgers; identity extends to `E_final + E_cost − E_bound − E_delivered = E_initial`. 50→25 ps sensitivity fell 70%/80% → 14%/3%; bound adjustments fell to +4.5 mJ vs 0.56 J accepted; the 256-cell 0.15% case reaches `1.6e20 m^-3`, 4 kA, 1.45 J, mJ lasing. See `bug-fixes.md` 2026-08-01. |
| Critical | One closure could not serve production conservation and Chen's published continuous source interpretation without reintroducing double billing. | Add explicit `source_coupled` (default) and `continuous_reconciled` (diagnostic-only) modes. Both consume one canonical compiled reaction-energy table. Production partitions electron sources from heavy chemistry and derives debit from the accepted augmented stage; the diagnostic path retains provisional continuous cooling and exposes a signed interval reconciliation. |
| High | UV absorption attenuated photons without consuming Ar**/F- or creating the declared charged products. | Apply PHOTO stoichiometry with the same exact exponential photon hazard used for attenuation. R40/R41 now consume their absorbers, create the declared ion/electron products, and add sourced 3.76/3.02 eV photoelectron birth energies; KrF R43/R44 use 3.10/1.60 eV. Atom, charge, photon, and electron-energy-inventory regressions cover the update. |
| Critical | Frozen-upper optical growth and end-density depletion created photons far beyond the excimer inventory. | Integrate exact frozen-coefficient photon fluence, limit internal growth/depletion, combine spontaneous and stimulated hazards, and integrate output over fluence. A gain-volume/cavity-volume conservation bound covers the high-gain case. |
| High | The MPC node B was advanced once with a plasma load in the circuit and again in transport; the transport result was discarded. | Use a Lie split: open-gap MPC network advance, gap-only finite-capacitance draw, then pin the returned voltage into the circuit. |
| High | Finite-capacitance gap energy used `I(V_old)V_old dt`, which is unbounded for `G dt/C >> 1`. | Use the exact capacitor-energy change `C(V_old^2-V_new^2)/2`; test that the gap cannot supply more than its stored energy. |
| High | Prescribed E fields were added at cell centres but omitted at faces, so they heated electrons without driving SG transport or wall flux. | Sample and superpose the prescribed axial field on all faces as well as cells. |
| High | Mean-energy coefficient tables were sometimes queried with E/N values. Native mean-energy solver output also ignored a requested E/N key. | Enforce the table key contract, implement both re-key directions, make the stub key-aware, and normalize injected drift-diffusion simulations to LMEA. |
| High | LUT disk keys used only E/N minimum, maximum, and count and omitted network/data identity. Distinct grids could return the same cached table. | Hash every grid point, the dynamic reaction signature, complete solver configuration, and JSON/CSV contents. |
| High | Saturated analytic temperature curves produced duplicate mean-energy abscissas, invalid for interpolation. | Consolidate exact ties deterministically and require strictly increasing table keys. |
| High | Interval-average output/input power was sampled at interval endpoints and trapezoid-integrated, dropping the first interval. One recorded interval produced zero energy. | Store each interval duration and integrate `sum(P_avg*dt)`. Cumulative efficiency uses the same convention. |
| High | Runs ending between chemistry cadences silently discarded the final chemistry, optics, and accumulated input energy. | Flush the final partial interval before final diagnostics. |
| High | Serial and parallel F2 sweeps omitted physical configuration and serial sweeps silently substituted the analytic solver. | Preserve backend, key kind, units, preionization, drive, mobility and cross-section paths, output root, and cache configuration. Reject ambiguous/invalid point selections. Injected mutable runtime objects require fresh per-point factories; process-parallel factories must be picklable. Only an explicitly supplied serial solver overrides the base backend. |
| Medium | Cross-section linear interpolation extrapolated outside the measured range; a differential out-of-range lookup returned `10^0 = 1 m^2`. | Clamp total tables to endpoints, return zero outside differential axes, reject negative/non-finite data, and consolidate duplicate total-energy rows. |
| Medium | The two-term process list could add an elastic/effective reaction dataset on top of the mixture momentum-transfer cross section. | Exclude elastic datasets from the inelastic energy-shift process list. |
| Medium | Rear-mirror reflectivity was accepted by the public cavity spec but omitted from photon loss. | Include `R_back` in round-trip loss and pass it from `Simulation`. |
| Medium | Cavity photon count used the pumped gain length although Eq. 7 stores a cavity-averaged density. | Use `n_ph*l_c*S`; upper-state depletion continues to use `l_g*S`. |
| Medium | Unknown or negative seed ions produced a non-neutral initial state without an error. | Require a registered positive ion and seed `n_i=n_e/z_i`. |
| Medium | Invalid grids, timesteps, pressures, optical geometry, and circuit values often failed later as NaNs or divisions by zero. | Add finite/range validation to the declarative specs and the chemistry integrator. |
| Medium | Decade-rounded and symmetric diagnostic limits hid pulse structure in whitespace; collectors also imposed figure names/composition. | Fit finite visible data with small padding and move mandatory names, x/y references, scales, and optional ranges into run-file plot declarations. |
| Cleanup | Cross-section APIs exposed channel classes and descriptions unrelated to coefficient/rate construction. | Restrict cross-section APIs to coefficient and reaction-rate use. |

The ArF Te-dependent stepwise-ionization channels R9 and R18 now retain their explicit energy
losses through TOML loading/compilation; their Arrhenius `E0` values are no longer implicitly
confused with energy bookkeeping. R10 attachment remains at its current zero-threshold treatment
until a sourced attachment-energy moment is available.

## Verification anchors

- Thomas solver: dense linear solve comparison, Numba/NumPy transport parity, bounds-checked kernel.
- Poisson: vacuum linear field, constant-charge Gauss law, alternating-charge response, and
  semi-implicit current-continuity identity.
- Circuit: lossless resonant-transfer energy conservation, passive conducting-gap collapse, finite
  gap-capacitor energy bound.
- Coefficients: Maxwellian closed forms, two-term thermal/Druyvesteyn limits, Einstein relation,
  key round trips, cache collision regression.
- Chemistry: SciPy Radau augmented-energy reference, co-integrated per-cell energy identity,
  live-energy rate refresh, atom/positivity checks.
- Optics: constant-coefficient analytic photon solution, fine-step ODE comparison, high-gain
  saturation, upper-state/photon inventory bound.
- Orchestration: final cadence flush, units equivalence, sweep smoke tests, diagnostic archive,
  run-owned plot naming, and data-aware line/heatmap scaling.

Run the principal checks with:

```bash
python -m unittest discover -s tests -p "test_*.py"
NUMBA_BOUNDSCHECK=1 python -m unittest tests.test_transport_kernel
python examples/arf_f2_sweep.py --tiny
```

The principal suite covers augmented-stage extent/debit identity, signed reconciliation in both
directions, PHOTO conservation and birth-energy booking, cache equivalence, reaction metadata,
Chen bank energy, sweep configuration cloning, and NumPy/Numba parity for both closure modes. Run
the kernel suite with `NUMBA_BOUNDSCHECK=1` after every kernel change. Test success is not a
convergence claim.

## Open numerical and modeling risks

These are visible limitations, not silently accepted implementation claims.

1. **Photoelectron energy:** the active ArF and KrF ionizing/detaching PHOTO rows now assign
   sourced point birth energies (`hv` minus binding/detachment): ArF 3.76/3.02 eV and KrF
   3.10/1.60 eV. These are not resolved birth-energy distributions; distribution-sensitive work
   needs sourced state-resolved data.
2. **Ion mixture coverage:** the mobility table lacks some ion/background-gas pairs. The fallback is
   finite and visible, but a quantitative mixed-gas result needs cited data for every active pair.
3. **Operator-split order:** source-coupled electron chemistry is Strang-split around SG
   transport with the transport-ledgered heating delivered inside the stage (rectified,
   lead/trail with exact per-interval conservation). Measured 50→25 ps sensitivity at 256
   cells/0.15%: n_e 14%, deposition 3% (was 70%/80% before the in-stage source) — but output
   remains threshold-sensitive (×2 across that step), so pulse-energy claims still need the
   dt gate run at the operating point, and the 12.5 ps member has not been run.
4. **Conducting-collapse closure:** the semi-implicit field removes the dielectric-relaxation
   instability, but field heating uses a supply-based cap in under-resolved sheath states. A
   cap-disabled 128-cell/50-ps sensitivity is effectively unchanged, so this cap does not cause
   the present pre-breakdown low-density branch. Treat post-collapse local `J.E` as unresolved;
   `Numerics(field_heating_cap_factor=None)` exists for explicit sensitivity runs only.
5. **Transport positivity:** charged transport still floors individual negative densities after a
   finite-volume update. Normal CFL-selected runs should not rely on this; cap warnings and species
   conservation should be monitored in severe collapse cases.
6. **Reaction energy completeness:** threshold losses are now co-integrated exactly with accepted
   progress and active photoionization/photodetachment rows declare point birth energies, but the
   network does not yet declare kinetic-energy yields for electrons born in exothermic Penning
   channels. Those yields require sourced state-resolved data.
7. **Effective momentum data:** when a source labels a base momentum-transfer set “effective,” adding
   explicit inelastic momentum contributions can still depend on the source convention. The audit
   removed direct elastic double counting, but dataset-level convention review remains necessary.
8. **Diagnostic reconstruction:** `EEDF` plots a Maxwellian reconstructed from mean energy. It must
   not be used as evidence that the runtime two-term distribution has that shape. Reaction-rate
   diagnostics are reconstructed from sampled states and are not an exact per-substep event ledger.
9. **Reaction-network provenance:** ArF data include cited corrections/additions beyond Chen Table
   III. `arf_chen_2025.toml` preserves the printed network for reproduction only. Sensitivity to
   uncertain branches and orphan charged species remains a physics-data question.
10. **Current convergence status:** the quoted 128-cell failure in the 2026-07-31 checkpoint was
    superseded by the in-stage heating-source correction; it remains useful defect evidence, not a
    current production result. The corrected 256-cell sensitivity reduced the major energy-book
    discrepancy, but output remains timestep-sensitive and no completed current mesh/timestep
    matrix authorizes predictive sweep claims.

## Required convergence work

Before treating sweep trends as predictive, record at least:

- mesh refinement in the electrode gap;
- outer timestep and chemistry-cadence refinement;
- E/N and electron-energy grid refinement for the two-term table;
- optics internal-tolerance refinement in the saturated pulse;
- sensitivity to missing ion mobility pairs and uncertain reaction branches;
- energy closure among circuit storage/loss, gap deposition, electron/chemical storage, cavity
  photons, and extracted optical output.

`examples/arf_convergence.py` implements the first two as 128/256/512-cell and 100/50/25-ps
successive-refinement gates, checkpoints every member, supports `--resume`, and writes
`runs/analysis/arf_convergence/convergence.json`. The 2026-07-31 pre-heating-source matrix failed
every comparison and is historical defect evidence rather than validation of the current model.
At that checkpoint, even 512 uniform cells had `dx/lambda_D ~= 5.3`; the next current-model mesh
study must resolve the electrode layer rather than interpreting 512 cells as a terminal grid.
