# PIC_MCC_3D changelog

All notable changes to the self-consistent 3-D PIC-MCC package.
Dated entries, newest first (Added / Changed / Fixed / Removed), per the
repo-wide changelog policy (`Agent/CLAUDE.md`).

**Read this file before editing anything in this directory**, especially the
Known-bad list — it is permanent, and it is where the traps that already cost
time are recorded.

The roadmap this package implements is `../PIC_MCC_3D.txt` (the nine missing
modules of its section 2.1, milestones M1–M10). Every module here closes a
numbered gap in that file; when a module lands, its board row there moves from
`[gap]` to `[ok]`.

## Known-bad — do not reintroduce

- **Do not solve the variable-coefficient (implicit/D1) Poisson equation with
  the separable FFT solver.** `Poisson3D` is fast *because* the operator is
  separable: a 2-D FFT in (x, y) turns the 7-point stencil into an
  independent tridiagonal solve in z. The direct-implicit screened equation
  `div[(1 + chi) grad phi] = -rho/eps0` has a spatially varying `chi`, which
  breaks separability outright. Use conjugate gradient *preconditioned by*
  `Poisson3D` (planned M7, `implicit3d.py`); do not "extend" the FFT path.
- **A fully periodic Poisson solve requires a net-neutral charge density.**
  The k=0 mode of the discrete Laplacian is exactly singular, so a net charge
  has no periodic solution — the potential would diverge. `solve()` raises on
  a net charge above tolerance rather than silently zeroing the mode and
  returning a wrong answer. Do not "fix" that by dropping the check.
- **Do not compare a Dirichlet solve against an analytic free-space Coulomb
  field.** Boundary planes at finite distance are image charges; the discrete
  answer is *correct* and the free-space formula is the wrong reference. Gate
  the discrete operator against the sparse-LU reference (exact) and against
  configurations with genuine analytic solutions (parallel-plate, vacuum
  ramp, MMS) instead.
- **Gather and deposit must use the SAME shape function.** Mixing (e.g. CIC
  deposit with a nearest-grid-point gather) introduces a self-force that
  accelerates isolated particles out of nothing. `es_field.py` uses trilinear
  (CIC) for both; measured, the sum of internal forces is 5.7e-15 of the
  largest single force and the self-force is 1.4e-16 of the one-cell Coulomb
  force, and gates hold both there.
- **Equilibrate the LU oracle's boundary rows.** `reference_solve_lu` is the
  gate oracle, and its Dirichlet rows (`phi = value`, entry 1.0) sit beside
  interior rows carrying 1/h². At h = 1e-3 that is a 1e6 row-scaling
  imbalance and the sparse LU loses ~6 digits — which presents as *the fast
  solver* disagreeing at 2.4e-10 when the fast solver is in fact the accurate
  one (it honours the planes to 2e-16). Scaling those rows by 1/h² is the
  same equation and restores 1e-15 agreement. Cost 2026-08-03: this was
  diagnosed only by checking which solution had the smaller residual — do
  that before believing an oracle.
- **Do NOT wire this package into `phantom11`'s pipeline.** Export and
  compare; never inject. Two defects are OPEN in the repo's own records:
  PLAN (70) — "the engine's electron-heating → light conversion is
  internally inconsistent by ~18x", with its own conclusion that no
  hybrid-pumping efficiency claim can be made until it is resolved; and
  PLAN (72) — the Kannari p→s relaxation defaults OFF as a *parity artifact*
  (τ ≈ 0.45 ns at 1 atm Ar, faster than the KrF* radiative lifetime), under
  which the bulk-deposited conversion is 54.2% and **breaches the ≤25%
  W-value bound**, while the baseline light output moves **3.4×** on that
  flag alone. Feeding a first-principles EEDF into that accounting does not
  improve the answer — it lends this package's provenance to a wrong one.
  Also note the parity harness is structurally blind to Ar_g-weighted terms
  (a sign error worth 40% of dArF survived there), so parity green is not
  evidence of correctness (PLAN (26) policy 4). Owner warning, 2026-08-03.
- **Never launch an electron above `CrossSectionSet.energy_ceiling_eV`.**
  `Process.sigma_at` returns ZERO above its last tabulated point — it does not
  extrapolate and does not raise. An electron there collides with nothing,
  never slows, and a degradation cascade returns "zero ionizations", which
  reads as a physical answer rather than a broken one. Measured ceilings on
  this repo's files: **Ar and He 1000 eV, Kr 965 eV, F2 937 eV** — well under
  the 2 keV the lookup table spans, and eight orders below a 600 keV beam
  primary. Fast electrons must be transported by stopping power
  (`beam3d`), never as MCC particles. Cost 2026-08-03: the first W-value
  measurement returned infinity.
- **Quote the DIFFERENTIAL W (dE/dU_ion), not E0/U_ion(E0).** The integral
  ratio includes the bottom of the cascade, where sub-excitation losses make
  no ion pairs, so it always sits above the asymptotic value the literature
  reports (argon: 35 eV at E0 = 100 eV). And it must be evaluated away from
  the table edge, where the yield saturates and bends W back up (argon 28.4 eV
  at 600 eV → 31.7 eV at 950 eV — an artefact of the data ending).
- **Do not locate the repo by counting parent directories.**
  `cross_sections.py` first resolved its data root as `_HERE.parents[2]`,
  which is correct in the main checkout and silently wrong inside a git
  worktree under `.claude/worktrees/<name>/` — two extra levels, so it
  pointed at a directory that does not exist and every LXCat lookup failed.
  Search **upward** for a known marker instead (`ARPAE/pic3d/lxcat`,
  `LoKI-MC/Code/Input`). Same rule the repo-wide `CLAUDE.md` states for the
  repo root: never hardcode a depth.
- **Do not fold the gas-reservoir term into the electron energy closure.**
  A newborn ion carries thermal energy drawn from the neutral gas. The
  electron ledger is an exact invariant *because* that term is booked
  separately (`gas_reservoir_eV`); merging them would produce an invariant
  that always "closes" while hiding a source, and the closure gate would stop
  detecting anything.
- **Do not deposit with `np.add.at`.** It is an unbuffered ufunc call and ran
  ~20–50× slower than the equivalent `np.bincount` over flat indices; at 1e5
  particles × 8 corners × 1e3 steps that is the difference between a gate
  that runs in seconds and one that times out. `es_field.deposit_rho`
  concatenates the eight corners into a single `bincount`.
- **`float()` on a ndim>0 numpy array is a hard TypeError on numpy ≥ 2.5**
  (the owner's Homebrew interpreter) and only a warning on numpy 2.0
  (`/usr/bin/python3`). This class of bug is invisible to the pytest run.
  Use `.item()`, and honour the two-interpreter rule below.
- **Never evaluate an avalanche counter at a fixed high bracket voltage.**
  M0 = ions per seed grows like exp(alpha d); a "safe" upper bracket of
  1500 V at pd = 8 Torr cm is alpha d ~ 37, i.e. e^37 particles — the first
  bisection evaluation ran 13 CPU-hours and produced nothing (2026-08-14).
  `tools/paschen.py` now climbs the bracket geometrically from a low
  voltage (overshoot bounded to x1.3 in V) and `measure_M0` caps the live
  macro count (a capped point is a LOWER bound, read as above-threshold).
  Any future threshold search over an exponential observable must do both.
- **The two-interpreter rule (PLAN (68)).** Run pytest under
  `/usr/bin/python3` AND import-smoke under the Homebrew `python3` before
  claiming anything green.

## [0.8.4] — 2026-08-14 — M5 emission + the C-tier experiments: G4.1 met, Ar Paschen vs Raizer, RF drive, EBSD i-V

Owner ask: add the Paschen curve and begin all five "one-build-away"
experiments of `SPECK/EXPERIMENTS.md` Tier C. One shared build unlocked two
of them outright (cathode emission), and the rest are begun with their
blockers named, not hidden.

### Added
- **`MODULES/emission3d.py` — MISSING MODULE 5 (partial), the gamma process
  + a flux emitter.** Exactly as the plan words M5: per-step ion areal flux
  to each electrode × γ → WHOLE super-particles through a fractional
  accumulator (nothing rounded away), emission energy distribution STATED
  (flux-weighted half-Maxwellian: v_n = v_th sqrt(-ln U), thermal
  tangentials). γ is a stated surface parameter, never a prediction; an
  ion-energy-resolved γ stays blocked on ion-neutral collisions (§2.1
  item 3). `FluxEmitter` is the constant-j cathode source G4.1 needs.
  Gates (`tests/test_emission3d.py`, 8): sampled moments match the stated
  distribution, accumulator emits whole particles and loses nothing,
  emitter delivers j·A·t/e to 1e-9, conservation emitted = collected +
  in-flight to 1e-12, γ bookkeeping exact (γ·21 = 8.4 → 8 macros + 0.4
  carried), wall attribution, determinism.
- **`es_sim.ESSim`: three small hooks** — `injectors` (run at the top of
  every step; `BeamInjector.inject(sim)` fits unchanged), `wall_hooks`
  (called with the indices about to be absorbed, per electrode, before
  removal) + `wall_flux_w[(species, 'lo'|'hi')]` per-electrode attribution,
  and `_bias()` so `phi_lo/phi_hi` may be a phi(t) callable (the C3/RF
  path). Old global tallies unchanged; architecture layering gated
  (`emission3d` at layer 2, imports species3d/units only).
- **`tools/diode.py` — Child–Langmuir, plan gate G4.1 MET.** Emit 3× J_CL;
  transmitted current measured from the anode's absorbed-weight series:
  J/J_CL = **1.023** at 250 V / 5 mm (virtual cathode −0.24 V at 39 µm),
  1.072 → 1.040 → 1.023 → 1.010 over 100 → 400 V — the residual is the
  stated O(√(T/V)) Langmuir thermal bias shrinking exactly as it should —
  and 1.017 / 1.012 / 1.023 / 1.043 across nz 33 → 257 (the virtual-cathode
  layer is ~1e-5 m; the band is flat at the few-% level). Artifacts
  `picmcc_out/diode/`. Gate `tests/test_diode.py` (3): within 10 %,
  virtual cathode exists, limiting real (3× offered, ~1× crossed).
- **`tools/paschen.py` — Townsend breakdown in argon vs Raizer, read
  FIRST-HAND** (`books/gas_discharge_physics_Raizer.pdf`, scanned — pages
  read as images): criterion Eq. (7.1) γ[e^{αd} − 1] = 1 (p. 131);
  similarity law (7.3)/(7.4) with Table 4.1 argon A = 12 /(cm Torr),
  B = 180 V/(cm Torr), E/p 100–600 (p. 56), inert-gas form (4.6) C = 29.2,
  D = 26.6; the printed experimental anchor "argon with an iron cathode,
  (pd)_min = 1.5 Torr cm, V_min = 265 V, (E/p)_min = 176" (p. 135) — the
  law + anchor imply γ_iron ≈ 1.5e-3 (our transcription reproduces 265 V,
  gated). M0 is MEASURED by direct avalanche counting (seeds at the
  cathode, every ionization counted, ions immobile — their return trip
  does not change how many there are); γ stays a stated cathode
  parameter; V_b(pd) is the bisection root of γ·M0 = 1. FINDING: α/p by
  the source's own two-length method (§4.1.4, slope of ln M0 vs d —
  cancels seed/back-diffusion/equilibration; single-length ln(1+M0)/pd
  under-reads by ×0.68) is **×0.80 of the Table 4.1 fit and ×0.87 of the
  inert-gas form** at 176 V/(cm Torr) (the two fits differ ×1.10 between
  themselves); ionization energy-sharing model moves it only 6 %. Gates
  `tests/test_paschen.py` (5): transcription reproduces the printed
  anchor, full electron ledger (seeds + births = absorbed at both walls),
  determinism, monotone in V, two-length α in a stated band, bracketing.
  Dynamic demo (`--demo`): real Ar⁺ macro-ions free-fall back and knock
  out real secondaries — the γ loop closed in time (ions free-fall: ion
  MCC open, stated; the criterion itself is transit-time-independent).
- **Paschen curve RESULTS (production run 2026-08-14, `picmcc_out/paschen/`,
  npz of record).** γ_iron = 1.5e-3: V_b = 436 / 376 / 373 / **352** / 365 /
  423 / 546 V at pd = 0.7 / 0.9 / 1.2 / **1.5** / 2.5 / 4 / 8 Torr cm —
  the minimum sits at pd ≈ 1.5 exactly where Raizer prints it, at
  **×1.33** the printed 265 V (the α deficit, ×0.80, propagated through
  e^{αd}: predicted ×1.2–1.3). γ = 1e-2: 278 / 260 / 259 / **264** / 300 /
  337 / 456 V — this curve passes THROUGH the measured Ar/iron point
  (264 vs 265 V at pd = 1.5). Two readings of the same physics, both
  stated: with the PIC's own α the iron anchor is reproduced by an
  effective γ ≈ 1e-2 rather than the 1.5e-3 the Table-4.1 α fit implies.
  Both PIC branches are flatter on the left than the similarity law
  (which diverges at (pd)_lim = 0.54); points with E/p outside the fit's
  100–600 window are flagged in the output (pd 0.7 and 8 at γ_iron; 4 and
  8 at γ = 1e-2). Dynamic demo at γ = 0.1, pd = 1.0: V_b = 155 V, the
  population GREW at 1.3 V_b and DIED at 0.7 V_b over three ion-transit
  cycles — the γ loop closed in time. α/p two-length table: ×0.82 / 0.76 /
  0.80 / 0.81 of Raizer (4.5) at E/p = 120 / 150 / 176 / 220; the first
  pass read ×0.47 at 300 and the 220/300 points were CAPPED long-gap
  avalanches (lower bounds, not physics) → `alpha_two_length` is now
  cap-aware (shrinks both gaps and re-runs; lengths used reported). FINAL
  cap-aware table (`alpha_ar.png`): α/p = ×1.02 / 0.90 / 0.85 / 0.79 /
  0.82 / 0.91 of Raizer (4.5) at E/p = 100 / 120 / 150 / 176 / 220 / 300
  V/(cm Torr) (gaps auto-shrunk to 0.36/0.72 cm at 300) — monotone, and
  closer still to the inert-gas form (4.6). So: the PIC argon chain
  under-reads α by ~20 % around the Paschen minimum and agrees at the
  window edges; that ×0.8 is what moves V_b ×1.33 at γ_iron.
- **`tools/rf_ccp.py` — the Turner 2013 CCP benchmark BEGUN, not claimed.**
  RF phi(t) drive gated against the closed-form single-electron trajectory
  (`tests/test_rf_bias.py`, 3: cos drive, < 1 % of amplitude over 3
  cycles; a sin drive from rest carries the physical DC drift u·t and
  walks out of any small gap — documented); case table (He, 6.7 cm,
  13.56 MHz, 450/200/150/120 V, n 9.64e20…3.21e22) transcribed from the
  lase-unb/ccp-benchmark reproduction repo — PROVENANCE STATED, verify
  against Phys. Plasmas 20, 013507 (DOI 10.1063/1.4775084) before
  claiming G5.1; scaled case-1 smoke closes the loop (both species reach
  both electrodes, 4k collisions in 2 RF cycles, stability honestly
  reported unresolved at the smoke grid). Blockers named: first-hand case
  spec, M6 Coulomb.
- **`tools/ebsd.py` — e-beam-sustained conduction, the non-self-sustaining
  i–V** (Raizer Fig. 7.1 regimes, read first-hand): `BeamInjector` in
  kinetic pair mode (W emerges) drives a 38 Torr Ar/Kr/F₂ gap along a
  periodic axis; with NO γ registered the gap cannot self-sustain by
  construction — exactly Raizer's statement, and the whole point of EBSD.
  Measured: **5.6 → 49.8 A/m² over 20 → 240 V**, ×5.8 amplification from
  60 → 240 V (consistent with e^{αd} at that E/p), ion space charge ≤ 10 %
  of the applied field and reported per point. Gates `tests/test_ebsd.py`
  (2). Artifacts `picmcc_out/ebsd/`.
- Suite 188 → **211** (`/usr/bin/python3`), plus the Homebrew smoke under
  `-W error::DeprecationWarning` over every new module (two-interpreter
  rule). The coupler that these experiments feed is `SPECK/tools/couple0d.py`
  (its own suite, 3; reproduce-limit parity 19 species / worst 4 %) —
  recorded in the root changelog, not here.

### Changed
- `es_sim` module docstring: what is still NOT in (circuit, dielectric
  charging, Coulomb, ion-neutral, implicit) vs the M4/M5 slices now in.
- `status.html`: M5 γ row → ok, flux-emitter row added, Child–Langmuir
  target → ok with the measured band, electrode-potential row notes phi(t),
  Turner row → warn ("begun, not claimed") with the blockers, suite 211.

## [0.8.3] — 2026-08-12 — two-stream demo: self-generated bunching, gated against the closed form

Owner feedback on the first particle views: "perfectly uniform with
absolutely no deviation or local filamentation/grouping." The Langmuir
demo is a single coherent mode by design; grouping that GROWS from the
dynamics is the two-stream instability, and it is exactly what this
package's honest regime is for.

### Added
- `--demo twostream` in `tools/visualize_particles.py`: two equal
  electron beams counter-streaming at ±v0 along z in the periodic
  jelly, v0 set so the box's mode-1 wavenumber sits at the
  fastest-growing κ = k·v0/ω_b = √3/2 — where the cold symmetric
  dispersion (closed form in `twostream_growth_rate`, nothing tunable)
  gives γ = ω_b/2. Anti-phase displacement seeds start the growing
  eigenmode; the macro-particles bunch, trap, and roll into the
  textbook phase-space cat's-eye. Frames show the 3-D beams (colored
  by beam), the (z, v_z) phase space where the vortex lives, and the
  beam-0 bunching trace with the closed-form γ as a dashed slope
  guide — falsification printed on every frame. Viewer.html colors by
  beam and shows the bunching readout.
- Recorded signal: the per-beam BUNCHING parameter |⟨e^{−ikz}⟩| —
  the anti-phase seeds cancel in the total-charge mode at t = 0 and
  the counter-drift makes it beat at k·v0, so the total mode is the
  wrong growth observable; the per-beam bunching is monotone through
  the linear phase (recorded per snapshot, in the npz).
- A small SEEDED Gaussian v_z spread (2 % of v0, deterministic)
  breaks a real degeneracy found while shipping: with the regular
  transverse lattice and an x,y-independent mode-1 field, every
  particle in a z-plane shares ONE (z, v_z) point forever — the first
  render's "vortex" was 12 points per beam. The spread is stated in
  the meta and far below v0, so the cold rate still applies — and is
  still gated.
- Gates (tests, 7 → 10): the fitted bunching growth rate matches the
  closed-form γ within ±20 % (measured on the 8³/2×216 config); the
  quiet start's z-histogram modulation is <0.05 and grows above 0.5
  (bins = lattice planes — more bins alias the lattice into fake
  modulation, learned here); schema round-trip incl. beam/bunch/mode_im
  and the beam-colored viewer payload.

### Measured
- Default demo (16³, 2×1728 macro-particles × 2.89e4 real e⁻, 10
  plasma periods, 127 snapshots): saturation at t ≈ 2 T_pe with the
  bunching rise parallel to the γ = ω_b/2 guide; run ~8 s; headless
  viewer check renders all 3456 dots + 12 box edges.
- Suite 188 passed under `/usr/bin/python3` (30.0 s; 185 + 3 new).

## [0.8.2] — 2026-08-12 — see the macro-particles: tools/visualize_particles.py

Owner feedback on REBUILD's 3-D viewer ("I do not see any particles or
macro-particles of the plasma, I instead see the cube just change
hues"): REBUILD's device grid is a rate-equation model — it has no
particles. The particles live HERE, so this package gets the viewer
that shows its own state.

### Added
- `tools/` package + `tools/visualize_particles.py`: snapshots REAL
  macro-particle state (positions + v_z) from an ES PIC run by reading
  the public particle arrays between steps — `MODULES/` untouched, and
  nothing in `MODULES/` imports the tool (gated). Stock demo is the
  cold-Langmuir frequency-gate configuration made visible (`cold_wave`
  quiet start, periodic neutralising jelly; amplitude 0.05 L vs the
  gate's 0.004, stated, so the motion is visible): dots colored by v_z
  show the mode-1 standing slosh, and the strip below is
  `density_mode`'s amplitude — the same signal the frequency gate
  measures. Outputs: `particles_data.npz` (the data of record),
  frame PNGs → `movie.gif` / `movie.mp4` (`--mp4`), and a
  self-contained `viewer.html` (uint8-quantized positions; time/yaw
  sliders, play). Scale honesty on every figure: reduced-density cm
  box — the plan §3 regime, not the laser chamber. Outputs default to
  `<root>/picmcc_out/` (gitignored).
- `tests/test_visualize_particles.py` (7 gates): snapshot shapes and
  bounds match the loop; the cold start carries only the leapfrog's
  retarded half-step kick (<1e-2 of the slosh peak); the mode amplitude
  swings through zero repeatedly; runs are deterministic; npz
  round-trips exactly; frame + GIF render headless (Agg); `viewer.html`
  is self-contained and its quantized positions invert within one
  8-bit bin. Plus the MODULES-never-import-tools gate.

### Measured
- Default demo (16³ grid, 1728 macro-particles × 5.79e4 real electrons
  each, 3 plasma periods, 96 snapshots): run ~2 s, 96 frames + GIF +
  MP4 + viewer; the board's "Run outputs" row moves gap → warn (the
  particle-state path is closed; flare's per-step .dat/.log discipline
  is still open).
- Suite 185 passed under `/usr/bin/python3` (27.4 s; 178 + 7 new);
  Homebrew 3.14.6 / numpy 2.5 smoke with `-W error::DeprecationWarning`
  covering import + a tiny run + every output kind.

## [0.8.1] — 2026-08-03 — retract the phantom11 coupling recommendation

Owner warning: "phantom was seriously flawed and would need careful
examination before implementation." Correct, and documented in the repo — the
0.8.0 recommendation to hand this package's EEDF/rates to phantom11 was
unsafe as written and is retracted.

### Changed
- `../PIC_MCC_3D.txt` §10.1 rewritten and §10.1a added: the Orestes item is
  now a **comparison, not a pipeline**. Export the EEDF, per-process rates
  and deposition profile with provenance; SCORE phantom11's assumed swarm
  tables against them; inject nothing.
- §10.6 priority 2 reworded to match, with the blocking conditions stated.
- Known-bad gains the rule, because this is exactly the trap a future
  session would walk into while "finishing the integration".

### Why (both OPEN in the repo's records, not speculation)
- **PLAN (70)**: the engine's electron-heating → light conversion is
  internally inconsistent by **~18×**, the CO₂ inverse-bremsstrahlung term
  takes the same path so the PLAN (67) numbers inherit it, and the entry's
  own conclusion is that no hybrid-pumping efficiency claim can be made until
  it is resolved.
- **PLAN (72)**: the Kannari p→s relaxation defaults OFF because 0.0
  reproduces og-pyfiles bit-for-bit — "a parity artifact, not a physics
  choice". At 1 atm Ar, τ ≈ 0.45 ns, faster than the KrF* radiative lifetime.
  Under the default, bulk-deposited conversion is 54.2%, **breaching the
  ≤25% W-value bound**; with it on, 1.7%. Baseline light output moves 3.4×,
  device efficiency 3.05% → 0.89%. Deliberately not flipped, because picking
  the flag by which value lands nearer the measured ~7% would be steering by
  an anchor (PLAN (26)). Undecided owner call.

The general principle, worth stating once: **a better input to broken
accounting is a more convincing wrong answer, not a better one** — and this
package's first-principles provenance would be what makes it convincing.

## [0.8.0] — 2026-08-03 — external survey; scattering-angle gates

Owner ask: look up what other codes have that this one may be missing (NRL
Orestes, Hydra/xRAGE, other PIC codes). Full findings with sources are in
`../PIC_MCC_3D.txt` §10; the headlines:

### Fixed — a verification hole the literature named and we had
- **Six scattering-angle distribution gates.** There was **no test on the
  angular distributions at all**. That is precisely the hole the PIC-MCC
  verification literature warns about (arXiv:2602.15130 — such tests catch
  "coding errors that do not modify the expected outcome from the collision
  algorithm"): a wrong angular model still conserves energy exactly, still
  tallies correctly, still reproduces n·sigma·v, and still matches the ARPAE
  oracle. Every gate this package had would have stayed green.
  - empirical CDF of cos(chi) vs the analytic one for all three models
    (isotropic / screened / Surendra), by KS statistic;
  - isotropic directions uniform on the sphere (⟨cos⟩ = 0, ⟨cos²⟩ = 1/3, and
    uniform azimuth);
  - **the rotation deflects by EXACTLY the sampled angle** (1e-12) — an
    exact check of the basis construction in `_scatter_dirs`, independent of
    any distribution, including its near-x special case.
  - The samplers were **correct**; two of my first CDFs were wrong (both
    models map p → cos *decreasing*, so the CDF is 1 − p(x), not p(x)).
    Recorded because the failure looked exactly like a broken sampler.

### Found — where this package sits relative to Orestes
NRL's **Orestes** is e-beam deposition + 23-species non-equilibrium chemistry
(116 reactions, Kr/Ar/F2) + 1-D lasing + 3-D ASE transport. **This repo
already has that architecture — in `phantom11` + `unified/ase3d` +
`unified/ebeam`, not here.** PIC_MCC_3D sits one layer *below* Orestes: it
prices kinetically what an Orestes-class code assumes. So the gap is not
chemistry (out of scope by design) but **the coupling**: handing this
package's EEDF, per-process rates and deposition profile to phantom11 in
place of its assumed swarm tables. That is the highest-value new item the
survey produced and it was not previously on the M-list.

Hydra/xRAGE: rad-hydro, not PIC, and not publicly documented at a useful
level — as anticipated. Their transferable ideas are architectural, and those
already arrived from `projects/MHD` (flare), which is readable.

### Still open, now with sources (§10.5)
Self-consistent neutrals via DSMC coupling (PICLas, Pantera) — the principled
retirement of our static-gas caveat; energy-preserving explicit PIC+MCC
coupling (arXiv:2410.14859), directly relevant to our split
electron-ledger/gas-reservoir bookkeeping; **particle management changes
answers, not just cost** (CPC 2016), which promotes M8 from a performance
chore to a correctness item; complex geometry (a hibachi/foil is not a box);
AMR and GPU/SIMD with particle sorting. Also: an **open reference
implementation of the Turner 2013 benchmark exists**
(github.com/lase-unb/ccp-benchmark), so the M5 comparison data need not be
re-derived from the paper.

### Verified
- **178 passed in 24.3 s** (`/usr/bin/python3`), suite 172 → 178.

## [0.7.0] — 2026-08-03 — verification: each module alone, then together

Owner ask: is there a way to verify each module works independently, and then
together?

The suite ran everything in ONE process, which is exactly the condition under
which cross-contamination hides — a module that leaks state, or a file that
only passes because an earlier one ran first, still shows green. Split into
the two questions it could not answer.

### Added
- **`verify.py`** — the runner, `python3 PIC_MCC_3D/verify.py` (`--quick`
  skips the slow sweep). Exit code 0 only if everything passes, so it is
  usable as a gate. Six sections: import graph, module isolation, per-file
  gates, together, composition ladder, summary.
- **`tests/test_architecture.py` (30 gates)** — the package's *shape*:
  - the import graph is **acyclic** (real topological sort) and **layered**,
    with the intended layering written down so an import that inverts the
    architecture fails here rather than surfacing later as a circular-import
    crash;
  - third-party imports are declared — scipy is allowed **only** in the
    Poisson LU oracle, which is gate-only code, so the package still runs
    without it;
  - every module imports **alone in a fresh interpreter**, and importing the
    package prints nothing, touches no data file, and writes no stderr;
  - **module-level state does not move.** This is the PLAN (74) rule
    mechanised: that entry traced **10 of 17** suite failures to `krf0d`
    mutating `phantom11` globals without restoring them, where one leaked
    value (`CHOW_EN_TD_DEFAULT`, 51 → 90 Td) was enough on its own to break a
    parity gate. The test diffs every module-level scalar, array and container
    across a representative run;
  - the same seed gives the same digest twice in one process — a second net
    for anything caching across runs.
- **`tests/test_integration.py` — the composition LADDER (5 rungs)**, built
  that way deliberately: each rung adds exactly one module and re-checks the
  invariants that already held, so a failure names the module that introduced
  it. A single everything-at-once test would leave eight suspects.
  1. fields only · 2. + kinetic particles · 3. + gas/cross-sections/collisions
  · 4. + beam injection & preionization · 5. + persistence.
  - Rung 3 checks the **charge budget across three modules at once**: seeded
    charge + created pairs − collected at walls must still balance.
  - Rung 5 is the capstone: a **bitwise-exact restart of the full live
    stack** — two species, immobile-charge grid, collision tallies, fields,
    and the single RNG stream the MCC and the beam both draw from. It can
    only pass if every module put all of its state where the checkpoint could
    see it.
  - The beam injector is handed the **sim's** rng, not its own. One stream per
    run is what makes rung 5 repeatable; two would have made it unrepeatable
    for no benefit.

### Verified
- `verify.py` → **ALL CHECKS PASS** under `/usr/bin/python3`:

| check | result |
|---|---|
| import graph | acyclic, **5 layers** |
| module isolation | 12/12 import alone, no side effects |
| per-file gates (each file its own process) | **172 passed**, 11 files |
| combined suite, one process | **172 passed** |
| isolated total vs combined | **equal — nothing is order-dependent** |
| reversed file order | **172 passed** |
| composition ladder | 5/5 rungs |

- The suite is **137 → 172** gates. All 12 modules also import alone under
  Homebrew 3.14.6 / numpy 2.5.0 (which has no pytest — the reason the repo's
  two-interpreter rule pairs pytest on one with an import smoke on the other).

The equality of the isolated and combined totals is the number that answers
the question: every test passes both by itself and in company, and the suite
is order-independent in both directions.

## [0.6.0] — 2026-08-03 — infrastructure adapted from projects/MHD (flare)

Owner ask: look at `projects/MHD` for architecture this package is missing.

`flare` is a clean-room recreation of FLASH's *architecture*, and its
infrastructure is the part this package had skipped entirely: it had physics
and gates but no configuration system, no persistence, and no run identity.
The two highest-value pieces are adapted here; the rest is named below rather
than quietly left out.

### Added
- **`params.py` — declared runtime parameters** (from `flare/params.py`,
  FLASH guide ch. 10). Every knob declared with a doc, type, range and
  choices; overrides applied in a fixed precedence
  (`default < problem < file < cli < api`); every value remembers where it
  came from. Three concrete failures this closes:
  - **a misspelled override was silently ignored** and the run used the
    default — `set()` now rejects unknown names outright;
  - out-of-range values failed late and obscurely (`ppw = 4` surfaced as a
    strange dispersion result, not an error) — ranges are checked at
    assignment;
  - "why is `dt` this value?" was unanswerable — `source("dt")` now says.
  A lower-precedence source cannot overwrite a higher one, so a problem
  definition cannot clobber what the caller passed.

- **`checkpoint.py` — save a run, resume it bitwise-exactly** (from
  `flare/io.py`, FLASH ch. 9). Full state to `.npz`: particles, fields,
  `rho_heavy`, tallies, history, and the parameter set with its provenance.
  - **The RNG bit-generator state is saved, not the seed.** Reseeding would
    rewind the stream and silently change every collisional answer.
  - **The MCC is deliberately NOT serialised** — large, immutable, and
    rebuildable from the saved gas description, which `restore(mcc=…)` takes
    as an argument. Pickling the object graph would have been less code and
    much worse: it would freeze the cross-section tables into every
    checkpoint and break the moment a module was edited.
  - The FIELDS are saved rather than recomputed. They are a cache, but the
    cache from *before* the last push/move — re-solving on restore uses
    post-move positions and yields a different field. (Found by the exactness
    gate, which failed on `phi` alone while all 22 other items matched.)

### Verified
- **137 passed in 23.9 s** (`/usr/bin/python3`); params + checkpoint smoke
  clean under Homebrew 3.14.6 / numpy 2.5.0.
- `test_restart_is_bitwise_exact` — checkpoint, restore, and continue is
  identical to never having stopped. This gate is strong because **any**
  unsaved mutable state breaks it.
- `test_restart_diverges_if_the_rng_is_not_restored` — proves the gate has
  teeth: a restart that reseeds instead of restoring diverges, while the
  correct one matches. Without this, a checkpoint that forgot the RNG would
  pass every collisionless test.
- `state_digest` is built from **raw particle arrays**, not summary scalars,
  and a gate moves one particle by 1e-12 m to prove the digest notices.

### Still missing from flare's architecture (named, not silently skipped)
- **Run outputs**: a per-step integrals table (`<base>.dat`), a logfile
  echoing the resolved parameters with `[CHANGED]` markers, quicklook plots,
  and a CLI (`python3 -m flare run sod nx=512 --plot`). This package has no
  entry point at all and keeps its history in memory.
- **A `problems/` registry** — flare makes each benchmark a class with
  `DEFAULTS`, so a case is reproducible by name. Here every test still builds
  its configuration inline.
- **A driver with dt polling** — flare takes the min over each unit's
  `compute_dt` with `dtinit/dtmin/dtmax/tstep_change_factor` and lands exactly
  on `tmax`. That structure fits a PIC-MCC unusually well, where omega_pe·dt,
  the CFL condition and the MCC collision probability each want a say; today
  `dt` is a fixed input.
- **The Unit registry with no-op stubs** — flare's way of making absent
  physics safe to call. Here the same job is done by `if self.mcc is not None`
  checks, which will not scale to missing modules 4–8.
- **`flare diff a.npz b.npz`** (FLASH's `sfocu`) — checkpoint comparison for
  regression testing. `state_digest` is the seed of this.

## [0.5.1] — 2026-08-03 — flow diagram: readable routing, visible feedback loop

Owner report: "many lines that overlap and are hard to discern", and the
feedback loop was not apparent.

### Fixed
- **The 8 → 1 feedback loop was drawn but invisible.** It existed as an edge
  (`mcc → dep`) routed up a channel at x = 352 — which is exactly where five
  input Béziers crossed it, so it was lost in the tangle. It now **wraps the
  entire diagram**: out of the bottom of stage 8, under everything, up the far
  right, back across the top, into the top of stage 1. It crosses *nothing*,
  and it is drawn dashed in the accent colour and labelled "FEEDBACK ...
  next step" so it cannot be mistaken for another data tap.
- **Orthogonal routing replaces the Bézier curves.** Curves crossing at
  shallow angles are unreadable; right angles are not. Each dogleg now gets
  its **own** vertical channel, and the boxes were repositioned so most edges
  are dead straight — inputs sit at the height of the stage they feed
  (electrodes level with the Poisson solve, cross sections level with the
  MCC), and outputs sit at the height of their source.

### Verified — measured, not eyeballed
Ran the **real** drawing code under a stubbed DOM (JXA), recovered every
emitted path, and analysed the geometry:

| check | before | after |
|---|---|---|
| edge crossings | many | **0** |
| collinear overlap (edges drawn on top of each other) | present | **0 px** |
| edges passing through node boxes | — | **0** |
| content inside the viewBox | — | x 108–958, y 44–692 of 1000×726 |

That harness executes the shipped renderer rather than a reimplementation of
it, so the numbers describe what the page actually draws. What is still NOT
machine-checked is typography — text metrics are unknown without a browser,
so label placement relies on the 24 px inter-box gaps being wider than the
9.5 px labels.

## [0.5.0] — 2026-08-03 — the data-flow diagram

### Added
- **`status.html` gains a second view: DATA FLOW.** The board says what
  exists; this says where a number came from. A generated SVG of the whole
  chain — inputs (gas, cross sections, beam, preionization seed, electrodes)
  → the eight-stage step loop (deposit → solve → E → gather → push → move →
  boundaries → collide, with the collision products looping back) → outputs
  (fields, particle diagnostics, ledgers, validity, scoring).
  - **Built for debugging, which is what makes it worth having.** Click any
    box or arrow and it reports what crosses that edge, the **invariant** that
    must hold there, the **gate** that enforces it, and the **symptom you
    would see if it broke** — e.g. "isolated particles accelerating from
    nothing = deposit/gather shape mismatch", "Te climbing with no source =
    dx > lambda_D, the wall not a bug", "closure failing while ions look fine
    = reservoir energy folded into the electron ledger". 13 of the 19 edges
    carry the full annotation; the six pure output taps carry a label only.
  - Data-driven like the tree: `FLOW_NODES` / `FLOW_EDGES` are the only edit
    point, and the SVG lays itself out from them.
  - Still self-contained — no CDN, no build step. (The one `w3.org` string is
    the SVG XML namespace for `createElementNS`, not a network request.)
  - Verified by executing the graph through JXA: 18 nodes, 19 edges, every
    edge resolving to a real node, no duplicate ids, no overlapping boxes.
    As before the *rendering* is not machine-checked — there is no headless
    browser here.

## [0.4.0] — 2026-08-03 — M9: the e-beam, and the surface-discharge decision

Owner ask: do the e-beam next, and "I do not have any idea how to implement
the surface discharge without many bugs."

### The surface discharge: not simulated, and that is the answer

`seed_preionization` enters the preionized sheet as an **initial condition** —
matched, co-located electron/ion populations at a stated density — and does
not model the flashover that made it. Simulating a surface discharge needs
dielectric surface charging, triple-point field enhancement, field and photo
emission, and streamer propagation across the insulator: **missing modules 4
and 5 plus physics beyond either**, i.e. several unvalidated models at once,
whose errors would compound straight into the beam-pumped answer the project
actually wants. The repo has taken this same decision before for the same
reason — PLAN (67)'s Argon_UI enters "the no-arc surface-discharge sheet ...
as an INITIAL CONDITION seeded on Ar+". What the seed density *is* remains an
input the caller must justify from hardware; the function makes that
assumption explicit and checkable rather than buried.

The two populations are seeded **co-located** on purpose: sampling them
independently would leave a per-cell charge imbalance from noise alone and
launch a spurious space-charge field on step one.

### Added — `beam3d.py`

- **The fast primary is NOT a PIC particle**, and that is the module's central
  design decision. At 600 kV it is relativistic (gamma 2.17) *and* far above
  the cross-section tables, where every sigma reads zero — a macro-particle
  there would collide with nothing and fly through the gas unstopped, with no
  error raised. It is transported as a ray losing energy by stopping power;
  only its low-energy progeny become particles, which is the regime the LXCat
  data and the MCC are valid in. Same split the repo's 0-D `unified/ebeam.py`
  uses, so the two are comparable by construction.
- **Stopping power is first-principles**, not a table lookup: the relativistic
  Bethe/Møller collisional formula with ICRU mean excitation energies.
  Unlike the single tabulated constant in `unified/ebeam.py` it carries the
  right energy dependence, which is what shapes the deposition profile.
- `BeamSpec` (voltage, current density, pulse envelope, foil transmission),
  `deposition_profile` (CSDA ray march), `sample_opal_peterson` (exact
  inverse-CDF sampling of the Opal–Peterson–Beaty secondary spectrum),
  `BeamInjector` (deposition → particles, per step).
- **The W-value is an OUTPUT, by two independent routes** —
  `spencer_fano_yields` (deterministic Fowler recursion) and `cascade_yields`
  (Monte-Carlo through the actual MCC operator). `BeamInjector(pair_mode="w")`
  exists for cheap runs and **refuses to default the W-value**: it must be
  passed explicitly, so no run can quietly turn W into an input and then be
  quoted for it.

### Verified

- **121 passed in 24.3 s** (`/usr/bin/python3`); beam physics + preionization
  smoke clean under Homebrew 3.14.6 / numpy 2.5.0 with `DeprecationWarning` as
  an error.

| gate | result |
|---|---|
| Bethe/Møller stopping power vs NIST ESTAR @600 keV | He 2.9%, Ar 5.6%, Kr 5.7%, F2 1.9% |
| CSDA range, 600 keV in 1 atm Ar | 1.94 m; scales exactly as 1/p |
| **the two W routes agree** | **0.4% (Ar), 1.0% (Kr), 1.7% (He)** |
| W vs literature — **scored anchor, miss reported** | Ar 27.7 vs 26.4 (+5%), Kr 26.2 vs 24.2 (+8%), **He 55.1 vs 41.3 (+33%)** |
| injector energy conservation (kinetic mode) | secondaries carry the deposited energy to 1e-9 |
| preionization seed | exactly net-neutral, populations exactly co-located |

**Nothing is tuned to the literature W.** The Ar/Kr offsets are the expected
direction for a Spencer–Fano treatment with equal ionization sharing on a
momentum-transfer basis. **Helium's +33% is a real miss**, and it is pinned by
a gate so it cannot drift silently in either direction; candidate causes are
the equal-sharing assumption and the He set's channel coverage.

### Found while building

- **The cross-section tables end at ~1 keV, and above that every sigma is
  silently ZERO.** Measured ceilings: Ar and He 1000 eV, Kr 965 eV, F2 937 eV
  — far below the 2 keV the lookup table spans. This is exactly the failure
  the beam module was designed around, and it bit anyway: the first W-value
  run returned **infinity**, from a cascade launched at 1000 eV that never
  collided. `CrossSectionSet.energy_ceiling_eV` now measures it and the
  cascades refuse to start above it.
- **W is energy-dependent, and the table edge corrupts it.** Argon's
  E0/U_ion falls 35.3 → 28.4 eV as E0 goes 100 → 600 eV, then bends back **up**
  to 31.7 eV at 950 eV as the tables run out. So the reported quantity is now
  the *differential* W = dE/dU_ion over a clean 25–55% window of the ceiling —
  which is also what the literature value means.
- **The cascade's stopping criterion was wrong and slow.** Grinding every
  electron below the lowest *excitation* threshold chases a straggler cooled
  elastically at 2m/M ~ 3e-5 per collision: >6000 rounds for an answer settled
  after a few hundred. Stopping when no further *ionization* is possible gives
  the identical W in 260 rounds. `U_exc` is then a lower bound, reported as
  `exc_complete: False`.
- **`BeamInjector.inject` returned one number called "energy"** where two are
  needed: a 600 keV beam ranges metres, so a centimetre cell stops well under
  1% of what enters. It now returns `energy_in_J` and `energy_deposited_J`
  separately — the single name invited exactly the wrong comparison, and did.

## [0.3.0] — 2026-08-03 — M3 (electron half): gas, cross sections, MCC

Surveyed `GUI_v2.py`'s device chain (its tab-⑥ map: LXCat → LoKI-MC →
embedded tables → phantom11 → e-beam / discharge drive / in-loop EEDF →
cavity → ASE) and `projects/spectre/` (the Ne-buffered sibling; its unique
code is the 0-D `modules/spectre/` kinetics — chemistry that stays out of a
3-D PIC by scope) for what this package still needs. The survey's finding is
the entry below.

### Added

- **`gas.py` — the neutral background. A DEPENDENCY THAT WAS MISSING FROM THE
  ROADMAP.** `PIC_MCC_3D.txt` §2.1 folded it into "cross sections", and it is
  not the same thing: until now this package had no gas at all — no
  composition, no number density, no temperature — while every collision rate
  is proportional to one. `GasMixture` owns it once: mole fractions (used as
  given, never silently renormalised), ideal-gas density from (p, T),
  per-component densities and masses, mean free path, thermal speeds.
  - Scope is stated in the module, not assumed: the gas is **static, uniform
    and ideal** — not transported, not depleted by ionization, not heated.
    That holds at the ns windows this package can afford; `depletion()`
    exists so a caller can *measure* how far a run has strayed from it
    instead of trusting it. The 0-D heating/rarefaction model, when it is
    needed, is `pic1d/gas_heating.py` (README_LIMITS 2.7).

- **`cross_sections.py` — the LXCat loader.** Adapted from
  `ARPAE/pic3d/cross_sections.py` (itself adapted from the tested
  `pic1d/cross_sections.py` @ 4e2a9fc). Parser, `Process`/`CrossSectionSet`
  semantics, the effective→elastic synthesis and the fixed log-energy lookup
  tables are deliberately identical — the point of adapting a gated module is
  that its gates still mean something. Densities now come from a
  `GasMixture`; data is located by a search path (local `lxcat/` → ARPAE's
  copy → the repo's LoKI-MC tree, plus a `PIC_MCC_3D_LXCAT` override), so no
  third copy of ~6000 lines of tabulated data.

- **`mcc3d.py` — MISSING MODULE 3, the electron half.** Null-collision
  (Vahedi–Surendra) MCC, adapted from `ARPAE/pic3d/mcc3d.py`: same selection,
  kinematics, angular models, ionization sharing and energy-grouped adaptive
  substepping. **The one real upgrade: ionization can create a REAL ION
  MACRO-PARTICLE.** ARPAE's ions are immobile, so its MCC deposits frozen
  charge; this package has kinetic ions, so an event can spawn into an actual
  `ParticleSpecies` — which is what a sheath, an ion flux to a wall and
  therefore secondary emission (M5) all require. The immobile-charge path is
  kept for callers who do not want to pay for ion particles, and both are
  gated.
  - **The energy audit is split, deliberately.** The electron ledger is
    exact. A newborn ion, however, inherits its neutral's thermal velocity,
    so it arrives carrying energy from the **gas reservoir** — energy that
    was in the system but not in this model's particle set. That is booked
    as `gas_reservoir_eV`, *outside* the electron closure. Folding it in
    would leave an "exact" invariant that quietly absorbs a source term,
    which is how an energy audit stops meaning anything.

- **The MCC wired into the loop.** `ESSim(mcc=…)` collides once per field
  step, after push and move; per-process tallies and the energy ledger
  accumulate on the run; unregistered products deposit as immobile charge so
  the space-charge accounting closes either way. The step is now
  deposit → solve → gather → push → move → boundaries → **collide**.

- **`tests/test_gas.py`, `test_cross_sections.py`, `test_mcc3d.py`** and three
  loop-integration gates — the suite goes 50 → **96**.

### Verified

- `python3 -m pytest PIC_MCC_3D/tests -q` → **96 passed in 23.1 s**
  (`/usr/bin/python3`, numpy 2.0.2). Full driven-gap run *including the
  collision path* under Homebrew `python3` 3.14.6 / numpy 2.5.0 with
  `DeprecationWarning` promoted to an error — clean.

| gate | result |
|---|---|
| **ARPAE provenance oracle** — adapted MCC vs `ARPAE/pic3d/mcc3d.py` on a shared RNG stream | **exact equality** of tallies, energy ledger and post-collision velocities (20k electrons, ~128k collisions) |
| electron energy ledger closure | **2.0e-16** (gated at 1e-10) |
| ionization bookkeeping | 1 event → **exactly** 1 new electron + 1 ion |
| gas-reservoir term vs the born ions' kinetic energy | **1e-12** |
| measured collision rate vs n·sigma·v (monoenergetic beam) | within Poisson error on ~1e4 events |
| no inelastic events below threshold | **0** excitation/ionization at 0.05 eV in Ar |
| nu_m(1 eV, 1 atm Ar) | **2.0e11 /s** (the band ARPAE's own map probes) |
| Ar ionization threshold from the parse | **15.76 eV** |
| ideal-gas density at 1 atm, 273.15 K | Loschmidt to **1e-12** |

The oracle gate is the one that matters here: it is what carries the
pic1d/LoKI-MC validation across to this copy, instead of this package
re-asserting a lineage it has not earned.

### Found while building (reported, not silently patched)

- **`ARPAE/pic3d/units.py` has a rounding slip in its Loschmidt constant**:
  `2.6867811e25` where the exact ideal-gas value (and CODATA 2018) is
  `2.6867801118e25` — the 7th digit, 3.7e-7 relative. Physically irrelevant,
  but it desynchronises a bit-exact oracle comparison, and it is why the
  oracle gate hands **one shared density** to both operators instead of
  letting each compute its own. **ARPAE was not modified**: it is a separate
  package with 41 gates of its own and its measured numbers are recorded
  against that constant; changing another package's physics constant as a
  side effect of building this one is not this session's call. Flagged here
  for its owner.

## [0.2.0] — 2026-08-03 — `MODULES/` layout + the status board

### Changed
- **Every dependency moved into `MODULES/`** (owner instruction). The package
  root now holds only the things *about* the code — `CHANGELOG.md`,
  `status.html` — plus the harness that checks it (`tests/`). Moved with
  `git mv` so history follows the files. Relative imports between the modules
  were unaffected (they are all in the same subpackage); the test imports
  moved to the canonical path `PIC_MCC_3D.MODULES.<module>`.
- **No root-level re-export, deliberately.** Aliasing the modules on the
  package root would make `from PIC_MCC_3D import poisson3d` work while
  `from PIC_MCC_3D.poisson3d import Poisson3D` still failed — one working and
  one broken path to the same module is worse than a single explicit one.
  The reasoning is recorded in `MODULES/__init__.py` so it is not "tidied"
  away later.

### Added
- **`status.html` — the system map.** A grouped, collapsible board:
  top level `PIC_MCC_3D`; 1st level **inputs / black box / outputs**; 2nd
  level the sub-systems (field solve, particle↔grid, movers, collisions,
  boundaries, ledgers, scoring, gates); 3rd level the individual items, each
  with a status, its source file, and the measured number backing it where
  one exists. Filters by status, expand/collapse, light+dark. Self-contained
  — no CDN, no build step, open it directly in a browser.
  - The tree is **data-driven**: one `const MAP` literal at the top of the
    script is the only thing to edit when the code changes, and sub-groups
    nest arbitrarily deep by adding another `children` array.
  - Status vocabulary is `GUI_v2.py` tab ⑥'s (PLAN (53)/(60)), so this board
    reads the same way as the rest of the repo's maps.
  - **It is a hand-maintained snapshot, not a live probe** — unlike the Tk
    system maps in `unified/`, `ARPAE/pic3d` and `Argon_UI`, which re-run
    their import/file probes on demand. The page says so at the top, and the
    gate figures on it name the dated run they come from. If the code and
    this file disagree, the file is wrong.
  - Contents as shipped: 73 nodes, **51 statused items — 32 working, 3 open,
    16 missing**. Verified by executing the map and tally through
    `osascript -l JavaScript` (JXA), which also confirms the whole script
    parses; there is no headless browser on this machine, so the page's
    visual rendering has NOT been machine-checked.

## [0.1.0] — 2026-08-03 — M1 + M2: the electrostatic core and mobile species

First code in this package. Closes missing modules 1 and 2 of
`../PIC_MCC_3D.txt` §2.1 — the two that everything else waits on.

### Added

- **`units.py`** — SI constants and the plasma scales the gates are written
  against: `debye_length`, `plasma_frequency`, `thermal_speed`,
  `courant_dt_em` (used only to document the ES/EM timestep ratio the plan
  computes), and the `EV`/`E_CHARGE`/`EPS0`/`M_E`/`M_P` set.

- **`poisson3d.py` — MISSING MODULE 1, the blocker.** The 3-D electrostatic
  Poisson solve, `div grad phi = -rho/eps0`, in two boundary modes:
  - `"periodic_xy_dirichlet_z"` — the device geometry (electrodes at
    z = 0 and z = Lz, each a scalar or an (nx, ny) profile);
  - `"periodic"` — the verification geometry, fully spectral.

  Method: a 2-D FFT over the periodic (x, y) plane diagonalises those two
  directions, reducing the 7-point stencil to one **independent tridiagonal
  system in z per (kx, ky) mode**, solved by a Thomas algorithm vectorised
  over all modes at once. Cost O(N log N), exact to round-off, no iteration,
  no preconditioner, no new dependency. A `reference_solve_lu()` built on
  `scipy.sparse` assembles the same operator explicitly and is used only by
  the gates.

- **`es_field.py`** — the electrostatic field/particle coupling:
  `deposit_rho` (trilinear/CIC to nodes, periodic wrap in x/y, volume
  weighted), `efield_from_phi` (E = −grad phi; centred in the periodic
  directions, one-sided at the z electrodes), and `gather_E` (trilinear, the
  **same** shape function as the deposit — the momentum-conserving pairing).

- **`species3d.py` — MISSING MODULE 2.** A `Species` record (name, charge,
  mass) plus a structure-of-arrays `ParticleSpecies` store carrying
  (x, y, z, vx, vy, vz, w) with grow/compact/spawn/kill, Maxwellian and
  cold-uniform loaders, and per-species diagnostics. Movers: `push_es`
  (leapfrog kick, electrostatic) and `push_boris` (E and B), both taking
  `q/m` from the species — so **ions are ordinary kinetic particles here**,
  not the immobile grid charge of `ARPAE/pic3d`. `SpeciesSet` holds the
  mixture and carries `subcycle` per species (heavy species may move on
  k×dt); `ion_subcycle_ratio` suggests k from the mass ratio.

- **`es_sim.py`** — a minimal electrostatic PIC loop composing the above
  (deposit → solve → gather → push → move → boundaries), enough to make the
  physics gates real. Periodic and Dirichlet-electrode geometries; per-step
  energy/momentum history; absorbing electrodes with wall-charge tallies.

- **`theory.py`** — the analytic overlays that SCORE the code: the scheme's
  own **discrete cold-plasma dispersion relation** (derived in the module
  docstring), a precise zero-crossing frequency estimator, and the
  Child–Langmuir and Bohm targets written down now so the M4 gates cannot be
  authored to fit whatever the code later produces.

- **`tests/`** — 50 falsification gates, ~21 s under `/usr/bin/python3`.

### Verified

Measured on this checkout, 2026-08-03 — every number below was produced by a
run, not estimated:

- `python3 -m pytest PIC_MCC_3D/tests -q` → **50 passed in 20.59 s**
  (`/usr/bin/python3`, numpy 2.0.2 / scipy 1.13.1).
- Full-loop smoke under the Homebrew `python3` (**3.14.6**, numpy 2.5.0 /
  scipy 1.18.0) with `DeprecationWarning` promoted to an error — clean, no
  scalar-conversion failures.

| gate | result |
|---|---|
| separable solve vs sparse-LU oracle, Dirichlet-z (h over 3 decades) | **1.2e-15** relative |
| separable solve vs sparse-LU oracle, fully periodic | **1.3e-14** relative |
| MMS order of accuracy, Dirichlet-z | slope **1.973** (L2 over 4 refinements) |
| MMS order of accuracy, fully periodic | slope **2.025** |
| electrode planes carry their prescribed potential (incl. a profile) | **< 1e-12** |
| vacuum capacitor → linear ramp | **1.5e-16** relative |
| parallel plate, uniform rho → analytic parabola (stencil is exact on quadratics) | **6.3e-16** relative |
| solution satisfies the 7-point operator (discrete Gauss law), re-applied independently | **3.6e-15** |
| deposit conserves total charge | **3.2e-18** relative |
| gather of a linear field | **1.3e-15** absolute (field scale ~10 V/m) |
| sum of internal forces (momentum conservation) | **5.7e-15** of the largest single force |
| self-force on an isolated particle | **1.4e-16** of the one-cell Coulomb force |
| Boris: speed conserved in pure B over 20 steps | **7.8e-16** |
| Boris: rotation angle vs −2 atan(ω_c dt/2) | **1.6e-16** relative |
| **cold plasma oscillation, nx=24** | measured **f/f_pe = 0.99158** vs the scheme's predicted **0.99145** → **1.2e-4** |
| **convergence toward the continuum** | f/f_pe = **0.97438** (nx=16) → **0.99774** (nx=32) |

The oscillation gate is the one that matters: it is the first end-to-end
confirmation that deposit, solve, gather and push are mutually consistent in
3-D. Note what it is scored against — **the discrete dispersion relation of
this scheme**, not omega_pe with a generous tolerance. A PIC cold mode always
runs slightly slow by a computable amount, so scoring against the continuum
value would have hidden real error behind a loose band; scoring against the
derived prediction agrees to 1.2e-4, and a separate gate shows the measured
frequency converging to omega_pe as h → 0.

### Scope — what this package does NOT have yet

Stated so no one mistakes an early package for a finished code. Missing
modules 3–9 of `../PIC_MCC_3D.txt` §2.1 are all still open:
ion–neutral collisions (M3), electrodes/external circuit/surface charge
(M4, only fixed-potential planes exist), secondary emission (M5), Coulomb
e–e/e–i (M6), the implicit/energy-conserving mover (M7, so **dx < lambda_D
and omega_pe dt < 2 are hard requirements today**), particle merge/split
(M8), and beam injection (M9). There is no MCC in this package yet — the
electron–neutral operator to adapt is `ARPAE/pic3d/mcc3d.py`, which is
already oracle-matched against `pic1d.mcc_ref`.

### Note on location

`../PIC_MCC_3D.txt` §7 originally recommended extending `ARPAE/pic3d` rather
than starting a second 3-D package. The owner created this directory
(2026-08-03) and directed the work here, so §7 has been updated to match.
The two packages are complementary and neither imports the other yet:
`ARPAE/pic3d` is the **electromagnetic** instrument (Yee FDTD, horn, wave
coupling, prescribed plasma); this package is the **electrostatic,
self-consistent** one. The MCC, cross-section and deposition modules are the
natural things to share, and that consolidation is tracked as an open item
rather than done silently.
