# COMMERCIAL_AUDIT.md — what stops SPECK, SPECKLE and SPECKLE3D from operating at commercial-grade confidence

Date: 2026-09-03. Read-only audit; no engine code was changed. The
question answered here is technical: what in the code, the numerics, the
physics representation and the verification prevents its results from
carrying the confidence a commercial simulation code's results carry.
Licensing, packaging-as-a-product and market questions are out of scope
(a five-line note at the very end records the two data-rights facts found
in passing, because they bear on whether the data the code runs on can be
shipped, and nothing more).

Scope: `SPECK/REBUILD` 0.1.36 (2026-09-01), `SPECK/PIC_MCC_3D` 0.8.4
(2026-08-14), `SPECK/tools/couple0d.py`, `SPECKLE/` (2026-08-27c),
`SPECKLE3D/` (2026-09-02).

Method: every README, CHANGELOG and plan document read (`PIC_MCC_3D.txt`,
`ORESTES_CLASS_REBUILD.txt`, `missing_modules.txt`, `MESH.txt`,
`EXPERIMENTS.md`, `README_LIMITS.md`); four code-level passes with
file:line evidence (engine numerics and API, PIC internals, the study
layer, packaging/QA); both suites run today; and a numerical-confidence
probe run today on the engine's own presets (rtol ladders, repeatability
within and across processes and interpreters, save-grid, max-step and
restart-schedule sensitivity, rate-tier sensitivity — Appendix A, scripts
under the job's tmp directory, re-runnable). Findings are VERIFIED (read
or executed today) unless marked INFERRED. Status words follow the repo's
vocabulary: **ok** (at the bar), **warn** (partial), **gap** (absent).

The modeling policy is not in question: predictive-only (anchors score,
never steer) stays, and this audit does not propose calibration.

--------------------------------------------------------------------------------
CONTENTS
  0.  Verdict
  1.  What "commercial-grade confidence" means, operationally (six pillars)
  2.  Measured today: how converged and repeatable the engine actually is
  3.  Findings by pillar
        P1  numerical convergence and error control
        P2  solver robustness and predictability
        P3  model form: the physics not represented that moves the answer
        P4  validated accuracy and the missing accuracy envelope
        P5  reproducibility and result provenance
        P6  guarded operation: does the code know when it is outside its envelope
        P7  what the test suites gate, and what they do not
        P8  defects found in passing
  4.  Scorecard matrix
  5.  What already meets the bar — keep it
  6.  What it would take, ordered by confidence gained per unit of work
  Appendix A. The probe (tables, scripts)
  Appendix B. The engine's own stated omissions
  Appendix C. In-house seeds elsewhere in the tree
  Appendix D. Limits of this audit
  Note.       Two data-rights facts, recorded and not discussed
--------------------------------------------------------------------------------


================================================================================
0.  VERDICT
================================================================================

The engine's arithmetic is trustworthy; its physics envelope and its
guard rails are not yet at the level where a user can take a number
without knowing the code's history.

At commercial standard already (measured today, Appendix A): a 0-D or
single-cell REBUILD run is converged in output energy to about one part
per million at the default tolerance, is bit-repeatable within a process,
across processes and across restart schedules, agrees between two scipy
generations to a part per million, and closes its energy ledger on every
saved step. The structural tests (import layering, no module state, no
calibration knob, per-row mass/charge balance, channel uniqueness) are
stronger than most commercial codes carry.

What stops the results from carrying commercial-grade confidence, in
order of weight:

1. **Model form (P3).** The physics the engine does not represent is
   first-order for the questions it is asked: the electron gas is a single
   temperature reading fixed-composition swarm tables (the F⁻ ×0.19 and
   Ar** ×0.06 population misses sit exactly where distribution shape
   matters, by the code's own diagnosis); 3-D cells are independent
   reactors, so filamentation and arcing — the physics behind the
   "over-ionization" question — cannot occur in the model; the gas density
   is frozen while the gas temperature evolves; the foil is a transmission
   of 1.0; VUV radiation escapes untrapped; all 248-nm absorbers are one
   aggregate; the KINETIC tier covers one mixture at seven points. The
   PIC's scattering basis is known to overheat the tail, and it has no
   Coulomb or ion-neutral collisions.
2. **Validation without an envelope (P4).** The anchors miss by ×0.23 to
   ×2.0, the Electra pressure trend is not reproduced, the g0/α ratio sits
   at 3 against an empirical 4–12, the Turner CCP benchmark is not claimed,
   and no document states for which operating region and to what tolerance
   the code's outputs can be trusted. The scorecard is a manual CLI flag.
3. **Robustness on anything beyond one cell (P2).** The finite-difference
   Jacobian's factor state overflows (scipy warnings appear even in the 0-D
   probe); on multi-cell Ne fills it kills the run, and the cure — restart
   segments — lives in four copies in the study scripts, not in the
   engine. Wall clock is unpredictable: the same 2×2×1 trajectory took
   2,535 s in one run and 344 s in another; a 2×2×1 reheat twin took 22.8 h
   with 99 % of it in two 2-ns segments; the KINETIC-tier Ne-optimized case
   has never completed (7-day runs terminated; today's 10-minute probe
   budget expired on it too).
4. **Silent operation outside the envelope (P6).** Te, Tg, gain and
   density are clipped without a counter; swarm tables are extrapolated by
   a flat clamp above and a taper below; validity gauges are printed, never
   enforced; the PIC's stability limits are advisory; a relocated checkout
   silently skips the collision physics tests and reports green.
5. **Tests gate structure, not device-level convergence (P7).** No suite
   halves a tolerance, a timestep or a grid on a lasing case; the anchors
   are not a regression suite; determinism is checked on one 0-D 20-ns
   run; SPECKLE and SPECKLE3D, whose scripts produced the 3-D verdicts,
   have zero tests; a circuit resume defect sits untested.
6. **Provenance (P5).** No result file records the code version, commit,
   interpreter, dependency versions, seed or wall clock; the data the
   engine runs on is resolved from sibling directories; SPECKLE cannot be
   re-run in any environment on this machine; a 3-D verdict is checked
   against a hand-typed literal.

Section 6 orders the fixes by confidence gained per unit of work. The
single largest item is not physics: an analytically generated sparse
Jacobian from the network data removes the overflow class, the restart
hacks and most of the wall-clock unpredictability at once.


================================================================================
1.  WHAT "COMMERCIAL-GRADE CONFIDENCE" MEANS, OPERATIONALLY
================================================================================

A user of a commercial simulation code trusts a number because six things
hold, and each is demonstrable rather than asserted:

  P1  NUMERICAL CONVERGENCE. The answer does not move when the tolerance,
      the timestep, the grid or the save grid is refined, and the code can
      show that on the cases it is sold for.
  P2  SOLVER ROBUSTNESS. Runs finish, in a predictable time, without
      hand-placed restarts; failures are diagnosed by the code, not by the
      user.
  P3  MODEL FORM. The physics that moves the requested observable at first
      order is in the model, or its absence is bounded.
  P4  VALIDATED ACCURACY. Benchmarks are reproduced within a stated
      tolerance, and a document states the operating envelope and the
      expected error inside it.
  P5  REPRODUCIBILITY AND PROVENANCE. The same input gives the same output
      on any supported platform, and every output records what produced
      it.
  P6  GUARDED OPERATION. The code knows when a run has left the validated
      envelope — a clamp activated, a table extrapolated, a stability
      limit crossed — and says so loudly or stops.

Two supporting pillars: P7 the test suites gate all of the above, and P8
no known defects are left in the path a user's run takes.


================================================================================
2.  MEASURED TODAY
================================================================================

All runs: `/usr/bin/python3` 3.9.6, numpy 2.0.2, scipy 1.13.1, two BLAS
threads, `nice 19`, from `Kinetic_Modeling/`. Cases: the three CLI presets
`krf_beam_0d`, `krf_ne_0d` (both sub-threshold at 3 atm: output energy
~7e-7 J, so they measure the fluorescence-level output), the lasing
`arf_chen_circuit` discharge anchor (2.04 mJ), and the lasing
`ne_optimized` cell behind the SPECKLE3D verdicts (7.773 J, η 6.19 %),
built exactly as `tools/mw_benchmarks.py` builds it.

Relative change of laser output energy E_L when rtol is tightened one
decade (default rtol = 1e-6):

| Case | 1e-4→1e-5 | 1e-5→1e-6 | 1e-6→1e-7 | 1e-7→1e-8 |
|---|---|---|---|---|
| krf_beam_0d | 4.5e-6 | 5.7e-7 | −2.8e-8 | 7.5e-9 |
| krf_ne_0d | 7.8e-6 | −7.3e-7 | −5.5e-8 | 6.0e-9 |
| arf_chen_circuit | 3.3e-4 | 3.1e-5 | −1.7e-6 | — |
| ne_optimized (segmented) | — | −4.2e-6 | −9.1e-7 | −1.3e-7 |

Other sensitivities at the default tolerance (relative change in E_L
unless stated):

| Probe | Result |
|---|---|
| Same process twice; fresh process | identical to all 16 digits (both presets) |
| Save grid n_save 50 / 200 / 800 (krf_beam_0d) | E_L identical; the reported *peak* n_e moves by up to 0.11 % because peaks are sampled on the save grid |
| Save grid 121 → 401 (ne_optimized) | E_L identical; sampled peak n_e +0.31 % |
| max_step 2e-9 → 1e-9 → 5e-10 | −5.0e-8, +1.8e-9 (krf_beam_0d); −4.8e-8, −2.9e-9 (ne_optimized) |
| Restart schedule, ne_optimized: one segment vs standard 5 vs fine 9 | single = plain run to all digits; fine vs standard +3.7e-7 |
| Rate tier LOKI → KINETIC (krf_beam_0d) | E_L −3.8 %, peak n_e +11.0 % |
| Rate tier KINETIC, ne_optimized | did not complete within the 10-minute probe budget; CHANGELOG 0.1.35 records the same case never completing in 7-day runs |
| Interpreter 3.9.6 / numpy 2.0.2 / scipy 1.13.1 vs 3.14.6 / numpy 2.5.0 / scipy 1.18.0 | 9.6e-9 (krf_beam_0d, 2269 vs 2284 RHS evaluations); 1.8e-6 (arf_chen_circuit) |
| scipy `num_jac` factor overflow warnings during the 0-D ne_optimized probe | 4 RuntimeWarnings (silent in a normal run) |

Reading: at the default tolerance a single-cell run is converged to ~1e-6
in E_L on lasing cases and ~1e-7 on sub-threshold ones; the machine
arithmetic is deterministic; restart placement and max_step are
sub-ppm effects in 0-D. The two numerical quantities that are NOT
converged by construction are the sampled peaks (a diagnostics defect,
not a solver one) and anything that depends on the electron-rate tier,
which is a modelling choice worth several percent in output and tens of
percent in density even where the tier runs. These 0-D numbers do not
transfer to the multi-cell grid or to the 1-D study, where the recorded
timestep and segment sensitivities are large (P1).


================================================================================
3.  FINDINGS BY PILLAR
================================================================================

Evidence paths are relative to `Kinetic_Modeling/`.

--------------------------------------------------------------------------------
 P1  NUMERICAL CONVERGENCE AND ERROR CONTROL           REBUILD 0-D ok · grid/1-D warn
--------------------------------------------------------------------------------
- 0-D: converged and repeatable (section 2). The tolerance floor is
  policed (`SPECK/REBUILD/tests/test_defect_history.py:87-91` refuses
  rtol looser than 1e-3) and analytic trajectories are checked at rtol
  1e-11 (`tests/test_analytic_trajectories.py`).
- Multi-cell and 1-D: the recorded sensitivities are first-order, and no
  gate exists for them. SPECKLE's co-pump marginal went +16.6 % → −24.8 %
  → −5.8 % at dt 50 → 25 → 12.5 ps while the baseline itself moved
  12.8 → 16.9 → 19.4 mJ, "~+15 % per halving" — the discharge phase is
  unconverged at practical timesteps (SPECKLE CHANGELOG 2026-08-26/27).
  SPECKLE3D's 1-D verification: +5.74 % (dt 50) → +10.06 % (dt 25),
  "magnitude unconverged UPWARD"; a timing-ladder point went −4.37 →
  +0.76 → +2.79 % across three dt values, the coarse dip being "another
  discretization artifact" (CHANGELOG 2026-08-27b). These are honest
  records, but they mean the headline 1-D and 3-D magnitudes are quoted
  from unconverged ladders.
- Save-grid dependence of reported peaks: `Device._collect` takes
  `max()` over `t_eval` samples (`SPECK/REBUILD/MODULES/engine.py:757-790`),
  so peak n_e, peak P_out and any FWHM depend on n_save (0.1–0.3 % measured
  today). The Chen anchor scores peak P_L and FWHM this way
  (`scorecard.py:202-278`).
- `max_step = 2e-9` and the per-block atol are literals in `run()`
  (`engine.py:736-752`), not parameters, so a user cannot run a
  convergence ladder without editing source or monkey-patching (as the
  probe did). The atol table is copied verbatim into five other files
  (`tools/mw_benchmarks.py:199-204`, `tools/artifact_significance.py:195-200`,
  `tools/copump_verify.py:81-84`, `SPECKLE3D/mw_3d_study.py:69-74`,
  `SPECKLE3D/cold_plasma_vis.py:98-103`) — any future change diverges.
- The differencing floor of twin-pair marginals was measured once, by
  hand, for the 1-D model (7.8e-4 mJ, SPECKLE 2026-08-27c) and certified
  once for 0-D (`tools/artifact_significance.py`); it is not computed per
  study.
- PIC: order-of-accuracy is gated for the Poisson solver (MMS slope
  1.97/2.03, `PIC_MCC_3D/tests/test_poisson3d.py:93-107`) and the Langmuir
  frequency converges with h (`tests/test_es_sim.py:55`); nothing gates
  particles-per-cell or timestep convergence of a collisional quantity.

The bar: convergence ladders (rtol, max_step, dt, grid, particles per
cell) on the lasing cases as CI gates; peaks taken from the solver's
dense output; tolerances and step limits as declared parameters; one
atol builder.

--------------------------------------------------------------------------------
 P2  SOLVER ROBUSTNESS AND PREDICTABILITY                                   warn
--------------------------------------------------------------------------------
- One `solve_ivp(method="BDF")` with a finite-difference Jacobian
  (`engine.py:750-752`), `jac_sparsity` only when `n_cells > 1`. The
  diagnosed multi-cell failure is scipy's `num_jac` factor-state overflow
  on BOTH the sparse and dense paths (SPECKLE3D CHANGELOG 2026-08-28,
  correcting 0.1.34's attribution); the 0-D probe today produced the same
  overflow warnings without failing. The engine's own guard is a
  non-finite check that raises with the offending state named
  (`engine.py:597-621`) and a 1 m⁻³ density floor added because "the bare
  0/0 NaN'd the numerical Jacobian" (`engine.py:470-475`).
- The cure is outside the engine: restart segments implemented four times
  (`tools/hybrid_timing_study.py:111-138`, `tools/artifact_significance.py:189-215`,
  `SPECKLE3D/mw_3d_study.py:60-93`, `SPECKLE3D/cold_plasma_vis.py:81-171`),
  with the fine-prefix logic duplicated and one variant lacking the
  circuit-state tolerances (`mw_3d_study.py:69-74`; CHANGELOG 2026-08-31).
  `Device.run` itself has no segmentation, no events (beam end, threshold),
  no warm restart.
- Wall clock is not predictable: 2×2×1 baseline 2,535 s vs 344 s for the
  same trajectory ("BDF wall swings with machine state", SPECKLE3D
  CHANGELOG 2026-08-28c); the 2×2×1 port twin 82,001 s with 51,261 s and
  30,241 s spent in the 392→394 and 394→396 ns segments
  (`SPECKLE3D/runs/coldplasma_port_2x2x1.log:247-251`), six attempts over
  five days; 3×3 and 4×4 never completed. A cold resume into the deep
  afterglow "rebuilds BDF/numjac state from scratch and can grind minutes
  where a warm continuous run passes in seconds" (CHANGELOG 2026-08-31).
- KINETIC tier: the Ne-optimized case never completes (CHANGELOG 0.1.35:
  7-day and 2-day runs terminated; "kinetic 5/8 is the standing record");
  reproduced today within the 10-minute budget.
- Checkpoint resume is tolerance-exact, not bitwise (`checkpoint.py:8-15`),
  and its atol table covers five auxiliary states where the circuit has
  seven (`checkpoint.py:80-82` vs `circuit.py:57`): a circuit resume gives
  two filter states atol 1e3 — untested (`tests/test_checkpoint.py` has no
  circuit case).
- Per-RHS cost: a Python loop over 169 rows (`chemistry.py:393-435`,
  duplicated at 321-359) and 48 ASE chords per call; single-core.
- PIC: explicit leapfrog only; `ESSim.stability()` reports dx/λ_D and
  ω_pe·dt but nothing enforces them (`PIC_MCC_3D/MODULES/es_sim.py:302-315`);
  no merge/split, so an avalanche grows unbounded (the recorded 13-CPU-hour
  e^37-particle runaway, `PIC_MCC_3D/CHANGELOG.md:104-111`); no Coulomb
  collisions; the plan's own conclusion is that an implicit mover is "the
  only lever" and it is not built (`PIC_MCC_3D.txt:443-455`).

The bar: an analytic or generated sparse Jacobian (the network is data;
this is code generation); segmentation, events and warm restarts inside
`run()`; a wall-clock model and an ETA; stability limits that stop or
adapt; merge/split and an implicit option in the PIC.

--------------------------------------------------------------------------------
 P3  MODEL FORM — PHYSICS NOT REPRESENTED THAT MOVES THE ANSWER              warn
--------------------------------------------------------------------------------
The engine states most of these itself (`scorecard.py:31-115` OMISSIONS,
Appendix B; `missing_modules.txt` P1–P8). Listed by the observable they
move.

Electron kinetics (moves every rate):
- One electron temperature per cell; rates read from LoKI-MC swarm tables
  at three FIXED compositions (Ar/Kr/F2 94.9/5/0.1, Ne/Kr/F2 94.9/5/0.1,
  Ne/Ar/F2 96.4/3.5/0.1), chosen by mode and the sign of Ar ≥ Ne only
  (`loki_tier.py:293-300`) — the Electra 99.7 % Ar anchor runs on the 5 %
  Kr table (`scorecard.py:545-547`); lookup by ⟨ε⟩ = 1.5 Te, a Maxwellian
  mapping of a non-Maxwellian table (`loki_tier.py:242`); as F2 burns
  (tracked as `f2_burn_frac`) the distribution cannot respond. No in-loop
  Boltzmann solver in either engine. The population misses F⁻ ×0.19 and
  Ar** ×0.06 "sit exactly where shape matters" (OMISSIONS).
- KINETIC tier: seven mean-energy points 3.4–7.7 eV for Ar/Kr/F2 only
  (`kinetic_tier.py:34-90`), measured by the PIC whose own swarm chain
  agrees with LoKI-MC to ×1.6 (`EXPERIMENTS.md` A10, "consistency"
  provenance) and bridged below its floor by stated rules
  (`MAX_REL_ERR = 0.5`, `BLEND = 0.1`, `kinetic_tier.py:158-165`). The
  two tiers flip the co-pump sign (CHANGELOG 0.1.23/0.1.25) and differ by
  3.8 % in E_L and 11 % in peak n_e on the sub-threshold cell today.
- No electron–electron relaxation, no electron-impact stepwise rows
  beyond R8/KSTEP, and KSTEP alone (an EST row) decides the co-pump sign:
  the positive window is ×0.49–1.99 of its value (CHANGELOG 0.1.28).
- PIC: isotropic scattering on momentum-transfer cross sections by
  default, stated as "inherited, not solved" (`mcc3d.py:36-38`;
  `cross_sections.py:28-34` "do not present isotropic scattering here as
  validated"); the 1-D sibling measured the basis ~10× too hot in the tail
  (`README_LIMITS.md:120-124`). Anisotropic options exist with only
  CDF-shape tests.

- ADDENDUM (same day, REBUILD 0.1.37): the MAXWELL verification tier
  measured the shape sensitivity from the same cross sections. Threshold
  rates under a Maxwellian are ×1e3–1e6 the swarm values at the tables'
  low end and ×0.9–2.8 at the high end; sub-threshold rates agree within
  ×2. On the NRL benchmark the swarm closure holds Te ×1.15 where a
  Maxwellian collapses to ×0.20, so a Maxwellian family is not the
  accuracy path; and the LOKI runs themselves sit inside the tables'
  range for only 1–33 % of a beam pulse (`te_in_table_frac`), which
  makes the sub-floor rule, not the in-table shape, the closure that
  decides beam-mode answers (REBUILD CHANGELOG [0.1.37]).
- ADDENDUM 2 (REBUILD 0.1.38, the cross-section atlas): the sub-floor
  region now has a named owner. Ar and Kr have Ramsauer-Townsend windows
  (×165 at 0.25 eV, ×185 at 0.51 eV) in which the buffer's elastic
  energy-relaxation time reaches 998 ns, ten times a beam pulse, so the
  bulk cannot thermalise and sits hot and non-thermal. What cools it there
  is trace F2 vibrational excitation, exceeding the elastic term ×187 at
  0.3 eV. Both that channel and F2 attachment physically RISE toward low
  energy while the engine's flat+taper rule falls, leaving the engine low
  by ×142 to ×3670 at the operating point; exactly 2 of 14 buckets carry
  that signature and both are F2. The published-database spread is only
  ×1.6 to ×2.1, so a different evaluation is not the fix. The fix is
  measured swarm data below 0.5 Td, a data extension with the existing
  provenance rather than a model change, and it is now the highest-value
  item in this section (REBUILD CHANGELOG [0.1.38]).

Transport and gas state (moves inhomogeneity, arcing, E/N):
- 3-D cells are independent reactors sharing one photon density and the
  deposition profile (`engine.py:264-275, 568-578`); the gated diffusion
  (`params.py:371-375`, default off) is a scalar and cannot be combined
  with a nonuniform axis (`engine.py:160-162`). Filamentation and arcing
  cannot occur in the model; SPECKLE3D's Known-bad says so, and it is the
  physics behind the "during the discharge it hurts" question.
- Total density is frozen at t = 0 (`engine.py:91`) while Tg evolves
  (clipped 200–3000 K): no expansion, no E/N drift with density, no
  pressure state, no acoustics. Heavy-particle rates depend on Tg only as
  (Tg/300)^n on 13 krf / 5 arf rows (`chemistry.py:347-348`); the Ar C↔B
  detailed balance is fixed at 300 K (`network.py:1057-1068`).

Deposition and drive (moves the energy in):
- Straight-ray CSDA, single energy, uniform transverse profile, no
  scattering, straggling or backscatter (`beam.py:174-198`,
  `engine.py:264-277`); `foil_transmission` defaults to 1.0 and is not
  plumbed through the parameters (`beam.py:103,113`, `engine.py:224-229`).
  The beam-EEDF F2 channels (dissociation, vibration, attachment by
  secondaries) exist since 0.1.36 but default OFF; with them ON the engine
  "lands where Orestes lands" and the pressure trend is still not
  reproduced — the remaining owner is shared model physics (CHANGELOG
  0.1.36).
- One circuit netlist, 0-D, LOKI tier only, fixed Lm (`circuit.py:10-98`,
  `engine.py:122-141`); preionization is an initial condition; the PIC has
  no circuit, constant energy-independent γ, no field emission
  (`emission3d.py:44-47`).

Photons (moves g0, α, E_L):
- All absorbers folded into one r_x (`chemistry.py`; `missing_modules.txt`
  M7/P3), so the open g0/α miss (3.0–3.1 vs the empirical 4–12,
  `scorecard.py:848-849`) cannot be attributed to a species.
- VUV rows (Kr2*, Ar2* second continua) escape untrapped (P2 of
  `missing_modules.txt`): radiative loss overcounted.
- 0-D photon density with one axial mode shared by every cell
  (`laser.py:84-89`, `species.py:225-226`); the Electra four-pass is
  emulated by unrolling to 2 m with a lumped loss (`scorecard.py:532-542`).
- 48 closed-form ASE chords on volume-averaged gain in the loop
  (`engine.py:334-358`); wall/window reflectivities are literals
  (`ase.py:185-213`).

Rates themselves:
- 19 EST rows in krf mode (`R8, KSTEP, KPS_Ne, KPS_Kr, NK1, NK1b, NK3,
  NK4, NK6, NK8, HK2, HK3, NF3, HK5_fBv, HK5_fC, BC_Ne, CB_Ne, VREL_Ne,
  URAD_Bv`), bracketed ×1/3..×3; the ion-ion screening form spread is
  ×1.75 (OMISSIONS); the Kr2*/Ar2* lifetime spread ×24 between sources.
- Absent classes: wall/surface chemistry, impurities, F2 vibrational
  manifold beyond one lumped 0.2 eV sink (`chemistry.py:62, 481-494`).

PIC geometry and collisions:
- Two boundary conditions (`poisson3d.py:52-54, 98-99`), full-plane
  electrodes, no internal electrodes or mesh import; neutrals static,
  uniform, unheated (`gas.py:10-12`); no ion-neutral or charge-exchange
  collisions ("ions arrive at the cathode by free fall",
  `emission3d.py:44-47`); the resolution budget (`PIC_MCC_3D.txt:401-461`)
  puts a device-density 3-D run beyond any hardware; the affordable box is
  ~64³–96³ cells.

The bar: an in-loop Boltzmann solver with composition and Tg dependence;
real transport (drift–diffusion with sheaths) on the device grid; gas
density and pressure as state; foil and hibachi deposition; per-absorber
photon ledger; radiation trapping; the EST count driven down; in the PIC,
Coulomb and ion-neutral collisions and validated anisotropic scattering.

--------------------------------------------------------------------------------
 P4  VALIDATED ACCURACY AND THE ACCURACY ENVELOPE                             warn
--------------------------------------------------------------------------------
- Verification against the analytic canon is real and gated: MMS order
  1.97/2.03, FFT-vs-LU 1e-15, cold Langmuir 0.99158 vs the discrete
  0.99145, two-stream ±20 %, Child–Langmuir J/J_CL 1.023, W-values
  emergent (Ar +5 %, Kr +8 %, He +33 %), ledger closure 1e-8 to 1e-16
  class, bitwise PIC restart (`EXPERIMENTS.md` Tier A/C; suites).
- Validation against devices, as the scorecard itself reports (CHANGELOG
  0.1.36; `scorecard.py:846-856`): Kannari ×0.86 (R_OC 80 %) and ×0.23
  (20 %, threshold-pinned); Electra ArF ×1.06 / ×1.62 / ×2.03 with the
  pressure trend NOT reproduced; NRL population benchmark F⁻ ×0.19,
  Ar** ×0.06, g0 ×0.56, n_e ×1.53, Te ×1.22, F2 burn ×0.53; Chen n_e ×0.79
  class, Te ×0.69; Paschen minimum ×1.33 in voltage (α deficit ×0.80,
  pinned). The conditioning gauge (`scorecard.py:430-446`, dlnE/dlng0)
  correctly labels the ×0.23 point as ×9–11 amplified — good practice,
  and it also says that near threshold the code's error is unbounded.
- Not claimed: Turner 2013 CCP ("BEGUN, not claimed", `tools/rf_ccp.py:1-22`;
  case table transcribed from a reproduction repo, M6 Coulomb open);
  Landau damping, Bohm–Gross, Bohm sheath, ambipolar and ion-mobility
  gates do not exist; the PIC swarm chain vs LoKI-MC is ×1.6.
- No accuracy envelope exists anywhere: no statement of the form "for
  e-beam KrF amplifiers at 1–2.5 atm and ≥40 % gain margin, output energy
  within ×0.7–×1.3". The scorecard is `run.py --scorecard`, its results
  are changelog prose, and `status.html` (the hand-kept board) is nine
  releases stale (0.1.27 vs 0.1.36). The `tests/test_missing_modules.py`
  Chen gap pin is the only anchor that runs in the suite, and it pins the
  miss rather than a tolerance.

The bar: a V&V document generated by a script from the suite, with the
anchor table, provenance class, conditioning, and an explicit envelope;
anchors as a CI regression suite with tolerance bands; the model-form
omissions bounded with a number each.

--------------------------------------------------------------------------------
 P5  REPRODUCIBILITY AND RESULT PROVENANCE                                     gap
--------------------------------------------------------------------------------
- Reproducible where measured: bitwise within and across processes;
  1e-6-class across scipy generations; the 2×2×1 baseline reproduced "to
  the digit" through a different capture path (SPECKLE3D 2026-08-28c);
  `test_repeat_run_is_deterministic` exists but covers one 0-D 20-ns run
  (`tests/test_architecture.py:140-154`) — no multi-cell, circuit or
  KINETIC-tier determinism check; BLAS thread count is set by hand.
- No output carries provenance: `run.py` prints and returns a dict
  (`run.py:221-229`); tool JSON/npz files record no package version,
  commit, interpreter, dependency versions, thread count, seed, or wall
  clock. SPECKLE's report is `vars(args)` plus numbers; the upstream
  model's commit hash is never recorded although the study's validity
  rests on it.
- The data the engine runs on is not part of the engine: both
  `cross_sections.py` modules resolve LXCat files by walking parent
  directories for `ARPAE/pic3d/lxcat` or `LoKI-MC/Code/Input`
  (`SPECK/REBUILD/MODULES/cross_sections.py:35-62`,
  `SPECK/PIC_MCC_3D/MODULES/cross_sections.py:47-98`); the Ne file exists
  only under the vendored LoKI-MC tree. The Spencer–Fano route falls back
  silently to the W-table if the files are missing (`engine.py:253-258`).
- SPECKLE cannot be re-run: it imports a colleague's model from a
  home-relative path through two private methods (`SPECKLE/mw_1d_study.py:57-61,
  124-141`); the model's pinned environment (h5py, numba, Python ≥3.11)
  does not exist on this machine; the recorded runs used unpinned
  Homebrew versions.
- SPECKLE3D's comparison tool asserts the verdict against a hand-typed
  literal `9.565` (`SPECKLE3D/compare_runs.py:83-84`), the grid against
  `(2,2,1)` (`:412`), and identifies the heated layer by a threshold
  (`:463-465`); the run artefacts (103 MB) and the 511 files of the package
  are untracked in git; REBUILD 0.1.28–0.1.36 are uncommitted.
- Copied numerics diverge silently: deposition weights are max-normalised
  in SPECKLE3D and mean-normalised in SPECKLE (`mw_3d_study.py:127` vs
  `mw_1d_study.py:143-162`); the microwave envelope is duplicated
  (`mw_1d_study.py:85-96` = `engine.py:51-62`).

The bar: a provenance stamp on every result file; data shipped with the
engine and hashed in the stamp; determinism gates on the grid and circuit
paths; study inputs and verdicts as data, not literals.

--------------------------------------------------------------------------------
 P6  GUARDED OPERATION                                                          gap
--------------------------------------------------------------------------------
- Silent clamps in the RHS (`engine.py:469-490`): Te clipped to
  0.026–50 eV, Tg to 200–3000 K, n_e floored at 1 m⁻³, gain exponent
  clipped to ±60 (`engine.py:638`). None increments a counter or is
  reported; a run that spends its afterglow on the Te floor looks like any
  other run.
- Swarm tables are extrapolated: flat then exponential taper below the
  table, flat clamp above (`loki_tier.py:244-249, 262-268`); drift velocity
  clamped above the table (`:271-283`). The Known-bad rule "never compare
  rate tiers outside their tabulated overlap" exists because the clamps
  change the physics; the engine reports `te_in_table_frac` and
  `mw_cap_in_table` (`engine.py:803-871`) but proceeds regardless.
- Validity gauges (ionization fraction, F2 burn, in-table fraction,
  bridged buckets) are printed on the regime line (`engine.py:913-917`);
  no threshold stops or flags the run, and the regime line prints "no
  inter-cell transport" even when transport is on (`engine.py:919`).
- PIC: stability limits advisory (P2); σ = 0 above the last tabulated
  energy (`cross_sections.py:133`); electrons above the ceiling are a
  Known-bad, not a guard.
- Tests skip silently when the LXCat data is not found (nine sites, e.g.
  `PIC_MCC_3D/tests/test_mcc3d.py:29`, `test_integration.py:35`): a
  relocated copy reports green without its collision physics.
- Dead or misleading controls: `SPECKLE3D/mw_3d_study.py --n-save` parsed
  and never used (`:193`); `--help` description truncated; `MODULES/params.py`
  in the PIC declares a `wpe_dt_max` nothing reads; `PIC_MCC_3D/README.md:188`
  advertises a numba path that does not exist; `__version__ = "0.2.0"`
  and "There is no MCC here yet" in `PIC_MCC_3D/__init__.py:17-35`.

The bar: every clamp and extrapolation counted and reported; validity
thresholds that fail the run or mark the result; stability limits
enforced; tests that fail, not skip, without data; controls that do what
their names say.

--------------------------------------------------------------------------------
 P7  WHAT THE SUITES GATE, AND WHAT THEY DO NOT                                  warn
--------------------------------------------------------------------------------
Gated today (417 tests green, Appendix A): import layering and no sibling
imports; no module-state mutation; banned calibration names; per-row
mass/charge balance with a teeth test; channel uniqueness; ledger closure
on 0-D and one 6×1×4 grid; analytic trajectories at rtol 1e-11; the
Poisson MMS order; the cold Langmuir frequency; two-stream growth; MCC
per-process ledger and CDF shapes; Child–Langmuir; Paschen band; W-value
routes; bitwise PIC restart; tolerance-exact REBUILD resume (0-D);
determinism (0-D, 20 ns); the Chen gap pinned as a miss.

Not gated:
- device-level convergence in rtol, max_step, dt, grid or particles per
  cell on any lasing case (section 2 is the first such ladder);
- the anchors as a regression suite with bands (the scorecard is manual);
- determinism on multi-cell grids, with the circuit, or on the KINETIC
  tier;
- cross-interpreter numerical parity (the two-interpreter rule is an
  import smoke test);
- checkpoint resume with the circuit (the atol defect);
- anything in SPECKLE or SPECKLE3D (zero tests; the "26 falsification
  checks" in SPECKLE3D's changelog have no file in the tree), although
  the 3-D verdicts rest on those scripts;
- the collision physics when data is absent (skips);
- runtime regressions (the suite's slowest test is 283 s; no timing gate).

--------------------------------------------------------------------------------
 P8  DEFECTS FOUND IN PASSING
--------------------------------------------------------------------------------
- `checkpoint.py:80-82` tolerances 5 auxiliary states; the circuit has 7
  (`circuit.py:57`) → two filter states resume at atol 1e3; untested.
- Sampled-peak diagnostics depend on `n_save` (0.1–0.3 % today); the Chen
  anchor scores peak P_L and FWHM from them.
- `engine.py:919` regime line asserts "no inter-cell transport"
  unconditionally.
- `scorecard.py:31-115` OMISSIONS still says the RLC filter branch "is
  not modeled" although `circuit.py` carries it (7 states).
- `SPECKLE3D/mw_3d_study.py:193` dead `--n-save`; `:69-74` atol block lacks
  circuit states (scope-bounded, but a trap).
- `SPECKLE3D/compare_runs.py:74,519-522` hard-codes `/opt/homebrew/bin/ffmpeg`
  and fails hard without it, while `visualize_3d.py` uses `shutil.which`.
- `PIC_MCC_3D/__init__.py:17-35` version and MCC statements stale;
  `tools/rf_ccp.py:6` claims data "ships in MODULES/lxcat/" (does not
  exist); `README.md:178` cites a guide that does not exist;
  `SPECK/REBUILD/status.html` stale at 0.1.27.
- `beam.py:103,113` foil transmission unreachable from parameters
  (effectively 1.0 everywhere).
- `tools/mw_benchmarks.py:199-204` and four other copies of the engine's
  atol table; four copies of the segmented loop.


================================================================================
4.  SCORECARD MATRIX
================================================================================

ok = at the bar · warn = partial · gap = absent · — = not applicable

| Pillar | REBUILD | PIC_MCC_3D | SPECKLE (driver) | SPECKLE3D |
|---|---|---|---|---|
| P1 numerical convergence | ok (0-D) / warn (grid) | warn | gap (dt-unconverged, no gate) | warn (unconverged upward, no gate) |
| P2 solver robustness | warn | warn | — | gap (restart hacks, 22.8 h) |
| P3 model form | warn | warn | (upstream 1-D fluid, 45 rows) | warn (independent reactors) |
| P4 validated accuracy | warn | warn | gap (scored against nothing) | gap |
| P5 reproducibility / provenance | warn | warn | gap (cannot re-run) | gap (literal-checked verdict) |
| P6 guarded operation | gap | gap | — | gap |
| P7 tests gate the above | warn | warn | gap (0 tests) | gap (0 tests) |


================================================================================
5.  WHAT ALREADY MEETS THE BAR — KEEP IT
================================================================================

- The energy ledger integrated as state and closed on every saved step;
  residuals normalised by energy movement.
- Provenance tiers on every rate with the census printed; EST rows
  bracketed by tests; the completeness-pass rule.
- The network as data with mechanical balance, uniqueness as a hard error,
  and RHS assembled from rows (the r115-class sign error is structurally
  impossible).
- The architecture suite (layering, no state, determinism, no-knob scan),
  the teeth tests, `verify.py` isolation ladders, the two-interpreter rule
  as a written gate, bitwise PIC restart.
- The conditioning gauge on anchors (amplification stated, never hidden),
  the "predictions, scored against nothing" labels, the resolution-wall
  computation, and the Known-bad lists with root causes.
- Measured today: 0-D convergence and determinism at a level a commercial
  code would advertise.


================================================================================
6.  WHAT IT WOULD TAKE — ORDERED BY CONFIDENCE GAINED PER UNIT OF WORK
================================================================================

Effort classes are judgment: S days, M weeks, L months, XL a quarter or
more.

 1. Generate the sparse analytic Jacobian from the network rows (code
    generation over the `Reaction` data; cells block-diagonal, photon and
    ledger rows global). Removes the `num_jac` overflow class, the restart
    hacks and most of the wall-clock unpredictability; speeds BDF.     [M]
 2. Move segmentation, events (beam end, threshold, window edges) and warm
    restarts into `Device.run`; one atol builder; declared `max_step`;
    fix the circuit-resume atol; peaks from dense output.               [S]
 3. Convergence gates in CI on the lasing cases: rtol, max_step, grid,
    n_save; dt and segment ladders for the study drivers; determinism on
    the grid, the circuit and the KINETIC tier; cross-interpreter parity
    with a numeric band.                                                 [M]
 4. Guard rails: count and report every clamp and table extrapolation;
    validity thresholds that fail or mark a run; PIC stability enforced;
    tests fail without data.                                             [S]
 5. Anchors as a regression suite with bands, and a generated V&V
    document with an explicit accuracy envelope and a bounded number for
    each model-form omission.                                            [M]
 5b. EXTEND THE SWARM TABLES BELOW 0.5 Td (new, 0.1.38 — placed ahead of
     item 6 by measurement). LoKI-MC runs at 0.01–0.5 Td carry the tables
     from <eps> 2 eV down toward 1.5 kT_gas = 0.039 eV and cover the whole
     beam-mode operating band, removing both the extrapolation the engine
     runs on for 67–99 % of a pulse and the wrong-way taper on the two F2
     channels. Same provenance, no model change.                      [M]
 6. Electron kinetics: an in-loop two-term Boltzmann solver with
    composition and Tg dependence (Appendix C names an uncalibrated
    in-house seed), tables as a fast path only; extend or retire the
    seven-point KINETIC tier; validate anisotropic scattering against
    swarm data in the PIC.                                               [L]
 7. Model form on the device grid: drift–diffusion transport with
    sheaths (the filamentation question needs at least this), density and
    pressure as state, Tg-dependent rates with detailed balance at live
    Tg, per-absorber photon ledger (settles g0/α), radiation trapping,
    foil and hibachi deposition, F2 vibrational manifold, a decision on
    the beam-EEDF F2 default.                                            [XL]
 8. PIC fidelity: Coulomb, ion-neutral and charge exchange, merge/split,
    an energy-conserving implicit mover, circuit boundary; then claim
    Turner 2013.                                                         [XL]
 9. Provenance: a stamp on every result (version, commit, interpreter,
    dependency versions, seed, threads, wall clock, data hashes); data
    shipped with the engine; study inputs and verdicts as data; tests for
    the study drivers; SPECKLE either pinned to a recorded upstream commit
    in a locked environment or retired as evidence.                     [M]
10. Performance: compiled RHS/Jacobian kernels and cell-parallel
    chemistry; PIC compiled kernels and cell sorting; a wall-clock model
    with an ETA before each run.                                         [L]


================================================================================
APPENDIX A.  THE PROBE
================================================================================

Suites (today): `SPECK/PIC_MCC_3D/tests` + `SPECK/tools/tests` 214 passed
in 131.7 s; `SPECK/REBUILD/tests` 203 passed in 540.9 s (slowest
`test_engine.py::test_3d_gradient_case_runs_and_closes` 283.4 s, next
`test_3d_uniform_limit_matches_0d` 44.8 s). Homebrew 3.14.6 imports both
packages but has no pytest.

Scripts (re-runnable, read-only w.r.t. the repo):
`~/.claude/jobs/b2c334e1/tmp/conv/probe.py` (presets; `python3 -B probe.py
<tag> [--quick]`) and `probe2.py` (ne_optimized via
`tools/mw_benchmarks._case_configs`); results in `results_*.jsonl` with
E_L, deposited energy, η, sampled peaks, closure residual, RHS-evaluation
counts, wall time and interpreter versions per run. `max_step` was varied
by wrapping `scipy.integrate.solve_ivp` from the probe; nothing in the
repo was edited.

Raw E_L values at the default tolerance (J): krf_beam_0d 6.778310597629288e-07;
krf_ne_0d 7.974392905730435e-07; arf_chen_circuit 2.0399965575306923e-03;
ne_optimized 7.773223032343338 (plain) / 7.773223030440069 (standard
segments) / 7.773225877196952 (fine segments). Homebrew: krf_beam_0d
6.778310532459566e-07; arf_chen_circuit 2.0399928336188363e-03.

Recorded wall clocks (SPECKLE3D CHANGELOG, `runs/*/report.json`, logs):
1×1×1 ne_optimized 6 s; 2×2×1 off 2,535 s and 344 s (same trajectory);
2×2×1 native/shell/port twins 387/407/254 s; 2×2×1 port via
`cold_plasma_vis.py` 82,001 s in the final of six attempts; 3×3 and 4×4
never completed; SPECKLE 1-D twins 30–115 s (median 50 s, 51 runs).


================================================================================
APPENDIX B.  THE ENGINE'S OWN STATED OMISSIONS (scorecard.py:31-115)
================================================================================

Kept verbatim in spirit; these are the code's confidence caveats and the
starting list for bounding model-form error:
- Xe rows omitted (trace); Xe enters stopping power only.
- Electron–electron collisions: "fluid Te treatment ... no EEDF shape
  below the tables; the F- x0.02 residual and Ar** x0.06 sit exactly where
  shape matters"; no e-driven excited-manifold stepwise rows beyond
  R8/KSTEP; superelastics omitted.
- Discharge circuit: fixed Lm; Xe seed omitted (the RLC-filter sentence
  is stale — the branch is now modeled).
- Below 0.47/0.66 eV the KINETIC tier "still extrapolates from the deepest
  measured anchor".
- Kr2*/Ar2* lifetime spread ×24 between sources, swept as a bracket.
- Beam-direct F2 channels: now built, default OFF (0.1.36).
- Ion-ion screening: ×1.75 form spread, "part of the recombination error
  budget".
- Checkpoint resume tolerance-exact, not bitwise.


================================================================================
APPENDIX C.  IN-HOUSE SEEDS ELSEWHERE IN THE TREE
================================================================================

Verified to exist; each carries the legacy line's parity history and would
need the same gates before use.
- `KrF_and_LoKI/KrF Code/factored-pyfiles/light_unified/eedf_multiterm.py`
  (507 lines): two-term solver with a parameter-free anisotropy limiter,
  "without any empirical calibration" — the candidate for item 6.
  `boltzmann.py` (185) adds the AC/microwave Joule term but carries a
  calibrated `FIELD_K`, disallowed as-is.
- `light_unified/gas_dynamics.py`, `preionization.py`, `pfl.py`,
  `circuit.py`, `light_chem.py` (Holstein trapping): seeds for item 7.
- `pic1d/implicit.py` (energy-conserving implicit mover), `radiation.py`
  (axially resolved Rigrod resonator, never wired in), `ion_transport.py`:
  seeds for items 7–8.
- `HIBACHI/` (pulsed power → diode → foil/hibachi ledger): item 7.
- `ARPAE/pic3d` (Yee/CPML electromagnetic PIC, 41 gates): the EM path.


================================================================================
APPENDIX D.  LIMITS OF THIS AUDIT
================================================================================

- The probe covers 0-D and single-cell runs; multi-cell and 1-D
  convergence statements come from the changelogs' recorded ladders, not
  from new runs (a 2×2×1 ladder is a day-class computation on this
  machine).
- The KINETIC-tier ne_optimized run was bounded at 10 minutes.
- Physics verdicts in the changelogs were not re-derived.
- Effort classes are judgment.
- One tool-generated claim was rejected: that pytest is unavailable here
  (the suites ran).


--------------------------------------------------------------------------------
NOTE — two data-rights facts, recorded and not discussed
--------------------------------------------------------------------------------
Every cross-section file the engines resolve is an LXCat download (Ar/He
IST-Lisbon, Kr Biagi-v7.1, F2 SIGLO, Ne Biagi) whose site terms restrict
redistribution and commercial use pending written permission per
database; and `SPECKLE` links at runtime to a colleague's model whose
LICENSE file is AGPL-3. Both bear only on whether the data and that one
driver can be shipped, not on the confidence of any result.
