# REBUILD — the Orestes-class excimer-laser engine

The rebuild `ORESTES_CLASS_REBUILD.txt` plans, housed where the owner directed
it (2026-08-04): `Kinetic_Modeling/REBUILD/`. A clean, self-contained,
fully-3-D device engine for the e-beam-pumped KrF/ArF cell — beam deposition,
non-equilibrium plasma chemistry, laser extraction, 3-D ASE — that runs on a
laptop and is structured to scale.

**POLICY (inherited, non-negotiable — PLAN (26)): predictive-only.** Every
rate comes from theory or a published measurement of that process. Anchors
SCORE; they never steer. Misses are reported as misses. No calibration knob
exists anywhere in this package, and a test asserts that.

## What this is (and what it is not)

Orestes (NRL) is e-beam deposition + 23-species/116-reaction non-equilibrium
chemistry + 1-D lasing + 3-D ASE transport. This package is that architecture,
built to be *more checkable* (C7 of the plan): every number traceable to a
citation, every assumption gated, every miss reported.

It is NOT a 3-D PIC code, and it deliberately does not contain one.
`PIC_MCC_3D/` is the kinetic layer — it prices, from first principles, the
electron behaviour this engine's rate tiers assume — and it stays where it is.
The coupling is DATA with provenance (the KINETIC tier below), never an
import. `PIC_MCC_3D.txt` §3 proves a λ_D-resolved 3-D PIC of the full cell is
unaffordable on any hardware this project will have; the affordable fully-3-D
device model is exactly this: a coarse device grid (cm-scale cells) carrying
chemistry + deposition + radiation, with kinetics supplied as data.

**On the root Known-bad "do not propose a per-cell 26-species network in
3-D":** that rule bans the network on a *PIC* grid (λ_D-resolved, ≥1e6 cells,
inside the PIC step loop) and it still stands. The grid here is the *device*
grid — O(1e3–1e4) cells of ~mm–cm size, integrated by a stiff ODE solver with
a block-diagonal Jacobian. Cost scales as cells × reactions per RHS call and
is measured in the changelog, not assumed.

## The five defect classes (ORESTES_CLASS_REBUILD.txt §2) — status here

- **D1 parity-driven physics defaults** — structurally impossible: this
  package owes bit-equivalence to nothing. There is no parity harness, no
  frozen reference, and no flag whose default is chosen by code equivalence.
  The PLAN (72) lesson is closed the only honest way: the Kannari p→s
  manifold relaxation rows are IN the network, unconditionally, as cited
  rows (no flag exists to turn them off).
- **D2 unauditable hand-written RHS** — the network is DATA (`network.py`
  rows with tier + citation); the RHS is assembled from it. Per-row elemental
  mass balance and charge balance are checked mechanically, and the teeth
  test proves a deliberately unbalanced row fails.
- **D3 no enforced energy budget** — the whole-engine energy ledger is
  integrated as state alongside the physics by the same solver, so closure
  holds analytically and any residual is integrator error. The suite fails
  if the budget does not close. (S3 — the headline gate.)
- **D4 calibration paths** — no fittable scale factor exists in live code;
  an architecture test scans the declared parameters and the source for
  banned patterns.
- **D5 mutable module-level state** — immutable config objects passed
  explicitly; the architecture suite diffs module-level state across a run
  (the PLAN (74) rule, ported from PIC_MCC_3D).

## The three layers, and the tiers

```
PIC_MCC_3D  (kinetic, 3-D, first-principles electrons)
     |  exports DATA with provenance (EEDF, per-process k, W yields)
     v
REBUILD chemistry+beam  (this package: tiered network, SI, no parity debt)
     |
     v
REBUILD radiation  (1-D lasing + voxel-field 3-D ASE / parasitics)
```

Electron-impact rate tiers, ranked: **KINETIC** (PIC_MCC_3D export, read as
data, no import) > **LOKI** (curated LoKI-MC mixture tables, copied under
`MODULES/loki_tables/` with provenance) > **CHEN** (Chen 2025 Table III
verbatim) > **REPO** (value + citation lifted from phantom11/spectre, never
re-tuned) > **EST** (stated analogy, flagged, sensitivity-bracketed ×1/3..×3
by tests; S5 drives the count to zero). The KINETIC and LOKI tiers run side
by side and their difference is REPORTED, never averaged away (S4).

## Milestones (from ORESTES_CLASS_REBUILD.txt §6) → what implements them

- **S1 read the history, write the tests** → `tests/test_defect_history.py`
  (docstrings cite PLAN entries; the Known-bad lists of the root, ARPAE,
  PIC_MCC_3D and spectre were walked 2026-08-04).
- **S2 structural gates** → `tests/test_architecture.py` (acyclic layered
  imports, import-alone, no module-state mutation, determinism, no-knob scan)
  + `tests/test_network.py` (balance checks with teeth).
- **S3 energy budget as invariant** → `MODULES/ledger.py` + gates in
  `tests/test_ledger.py`: every joule in is accounted to a named sink each
  step; beam-vs-discharge marginal conversions delivered side by side with
  the W-value bound gate (the audit PLAN (70) could not do).
- **S4 kinetic electrons** → `MODULES/kinetic_tier.py` (+
  `tools/export_kinetic.py`, which imports PIC_MCC_3D *outside* the engine
  import graph) with the LOKI-vs-KINETIC comparison report.
- **S5 network completeness** → tier census in the scorecard; EST rows carry
  brackets; sensitivity sweep hooks in `tests/test_network.py`.
- **S6 lasing + ASE** → `MODULES/laser.py` (Rigrod-gated) +
  `MODULES/ase.py` (voxel gain field; uniform-field oracle vs the closed
  form; parasitic margins).
- **S7 score against anchors** → `MODULES/scorecard.py` — every anchor, its
  source, and the miss stated. No tuning. Ever.

## Layout

```
REBUILD/
  README.md          this file
  CHANGELOG.md       dated entries + Known-bad — READ FIRST
  status.html        THE SYSTEM MAP: inputs → engine → outputs → gates with
                     per-item status + measured numbers, plus the data-flow
                     view (one RHS evaluation, gates on every edge). A
                     hand-kept snapshot — if it and the code disagree, the
                     page is wrong. Open directly in a browser.
  verify.py          each module alone, then together; exit 0 = all green
  run.py             CLI entry: presets, run report, scorecard
  tools/             scripts allowed to import OTHER packages (never imported
                     by MODULES/ — gated): export_kinetic.py;
                     visualize_3d.py (per-cell fields over time → npz;
                     movie.gif = the GLOW view, what the cell visibly
                     emits, in the vessel; movie_particles/_fields.gif =
                     the analysis views; self-contained HTML viewer with
                     glow/particles/voxels modes; REAL macro-particles
                     are PIC_MCC_3D/tools/)
  MODULES/           the engine (import path REBUILD.MODULES.<name>)
    units.py         SI constants + plasma/validity scales
    provenance.py    tiers, citations, the provenance report
    params.py        declared parameters (flare/FLASH pattern, via PIC_MCC_3D)
    species.py       registry as data: mass, charge, atoms, E_form + tier
    chen2025.py      Chen 2025 Table III verbatim (R19/R26 flags preserved)
    network.py       the tiered reaction network as data + balance checks
    loki_tier.py     LOKI electron rates from MODULES/loki_tables/ (copied,
                     provenance in loki_tables/README.txt)
    kinetic_tier.py  PIC_MCC_3D exports as data (KINETIC tier) + comparison
    maxwell_tier.py  MAXWELL verification tier (0.1.37): Maxwellian / two-
                     Maxwellian EEDF rates integrated from the LXCat
                     sets; scores the LOKI shape assumption side by
                     side (beam drive only; never a default)
    cross_sections.py LXCat loader (adapted from PIC_MCC_3D, + Ne)
    beam.py          Bethe/Møller stopping, CSDA profile, Spencer–Fano
                     yields (W emergent), conserving per-pair ledger split
    ledger.py        the named-sink energy budget, integrated as state
    chemistry.py     RHS assembled from the network; vectorized over cells
    laser.py         cavity photon equation (Chen Eqs. 7–18 form) + Rigrod
    ase.py           voxel-field 3-D ASE chords + parasitic margins
    engine.py        Device: modes (arf_chen | krf), drives (beam |
                     discharge), geometry (0-D | 3-D grid), the 3-D
                     deposition extrusion, BDF integration with
                     block-diagonal Jacobian sparsity
    checkpoint.py    save/resume (tolerance-exact — adaptive BDF; the
                     honest semantics are in its docstring)
    scorecard.py     anchors scored, misses stated
  tests/             one test_*.py per module + architecture + defect history
```

Import path is single and explicit — `from REBUILD.MODULES.engine import ...`
— run everything from the `Kinetic_Modeling/` root:

```bash
python3 -m pytest SPECK/REBUILD/tests -q        # the gates
python3 SPECK/REBUILD/verify.py [--quick]       # isolation + composition ladder
python3 SPECK/REBUILD/run.py --preset arf_chen_0d
python3 SPECK/REBUILD/run.py --preset krf_beam_3d
python3 SPECK/REBUILD/tools/visualize_3d.py --preset krf_beam_3d   # SEE the 3-D
#   run: per-cell fields over time → SPECK/rebuild_out/viz_krf_beam_3d/
python3 SPECK/REBUILD/tools/mw_reheat_experiment.py   # the ARPA-chamber
#   MW-afterglow design study: baseline vs reheat, marginal efficiency per
#   injected power (see CHANGELOG [0.1.20] for the measured verdict)
python3 SPECK/REBUILD/tools/hybrid_timing_study.py    # HYPER-LASER (DE-AR)
#   timing study: afterglow / copump / preprime / inter-pulse keeper, each
#   vs its no-MW twin, laser AND 248.5 nm fluorescence metrics ([0.1.21])
python3 SPECK/REBUILD/tools/copump_scan.py            # the co-pump operating
#   map: marginal payoff of microwave DURING the beam over pressure, power,
#   pump rate, window placement, optics margin ([0.1.22])
python3 SPECK/REBUILD/tools/copump_efficiency.py      # the efficiency range
#   with co-pumping, from the map (intrinsic + per injected beam energy, and
#   the MW power per +1 percentage point)
python3 SPECK/REBUILD/tools/copump_verify.py          # is the co-pump gain
#   real? solver/integrator convergence + null floor, EST x1/3..x3 and
#   LIT-conflict brackets, LOKI vs KINETIC tier, shot-pair power analysis
#   ([0.1.23]: numerically real, sign model-fragile — KSTEP + tier flip it)
python3 SPECK/REBUILD/tools/efficiency_budget.py      # where the efficiency
#   goes: beam partition -> KrF* formation/loss channels (per-row fluxes)
#   -> fate of the stimulated photons by absorber ([0.1.24]; --mixture ne,
#   --f2, --set k=v)
python3 SPECK/REBUILD/tools/efficiency_levers.py      # the levers that raise
#   it: pressure / pump / pulse / cell / optics / F2 / Kr one-at-a-time +
#   greedy combination (reference 1.54% -> 4.95%; Ne optimized 5.76%)
python3 SPECK/REBUILD/tools/mw_benchmarks.py          # microwave reheat and
#   plateau co-pump on EVERY beam-drive benchmark (DE-AR class, Ne fill,
#   Kannari + Electra-ArF anchors, HYPER, ARPA, the optimized points) at the
#   250 kW/L rule; --tier kinetic for the sign-sensitivity pass ([0.1.26])
python3 SPECK/REBUILD/tools/ib_feasibility.py         # BOTE: inverse-
#   bremsstrahlung coupling vs frequency (1 GHz - CO2) with the engine's own
#   nu_m -- below ~0.3 THz the plasma can't tell frequencies apart; the
#   wall-plug inequality eta_src x marginal > 7% decides ([0.1.26])
python3 SPECK/REBUILD/tools/afterglow_nogo_audit.py   # every KrF* formation
#   row's response to heating Te 0.5 -> 3 eV (all heavy-particle; both supply
#   lines FALL when heated); --arf for the 193 nm chain ([0.1.26])
python3 SPECK/REBUILD/tools/mw_sign_decidability.py   # can the model TELL?
#   (A) co-pump marginal vs the KSTEP multiplier -> the positive window vs
#   the EST x3 uncertainty (the sign's decidability, quantified);
#   (B) full EST/LIT bracket pass on the one afterglow-positive case
#   (ne_optimized) ([0.1.28])
python3 SPECK/REBUILD/tools/artifact_significance.py  # is it a numerical
#   artifact? placebo-calibrated profile likelihood + exact Student-t +
#   toys + CLs band: Z = 16.5 sigma against the numerics-noise null,
#   self-testing, with the scope stated ([0.1.31])
python3 SPECK/REBUILD/tools/derated_operating_point.py # what survives
#   realistic coupling: eta_c > 62% for intrinsic gain, wall-plug
#   excluded; Fresnel favors 95 GHz; the blended-gain table ([0.1.33])
python3 SPECK/REBUILD/tools/maxwell_verify.py --device --population
#   the MAXWELL verification ([0.1.37]): k_MAXWELL / k_LOKI per
#   electron-impact bucket at the tables' own <eps> points, the
#   two-Maxwellian bracket at fixed <eps> (stated tail grid), the E/N
#   view for the discharge lane, and the beam-mode presets / NRL
#   population rows under loki / loki-taper / maxwell (--kannari adds
#   the Kannari anchor). Scores the shape assumption; steers nothing.
python3 SPECK/REBUILD/tools/sigma_atlas.py
#   the cross-section atlas ([0.1.38]): elastic sigma vs eV per gas with
#   the Ramsauer-Townsend windows and the swarm-table floors marked,
#   energy-relaxation times, which channel actually cools the bulk, the
#   published-database spread, and the SUB-FLOOR AUDIT (which buckets the
#   flat+taper rule extrapolates the WRONG WAY: F- and vib, in every
#   fill). Seven figures + report.json in rebuild_out/sigma_atlas/
python3 SPECK/REBUILD/tools/excimer_anchors.py [--quick] [--set k=v]
#   the MEASURED excimer anchors ([0.1.40]): Electra ArF small-signal
#   gain at printed deposition (matched exactly) and the Klimek KrF
#   gain/absorption grid, plus the pump-independent discriminators
#   (g0/alpha per pressure, the pressure exponent). Refuses to score on
#   the W-table fallback. Scores; never steers. ~2 min (--quick 30 s);
#   rebuild_out/excimer_anchors/<tag>.json
```

Two-interpreter rule (PLAN (68)): pytest under `/usr/bin/python3`, import +
run smoke under the Homebrew `python3` (numpy ≥ 2.5) before claiming green.

## Laptop honesty / scaling

Presets print their regime (grid, cells, reactions, wall-clock class). The
3-D grid is embarrassingly parallel over cells everywhere except the shared
cavity photon field and the ASE chords; the scaling path (cluster) is domain
decomposition over cells plus per-rank chord bundles — recorded here so the
single-core numpy v1 cannot quietly be mistaken for the cluster code.
MESH.txt's adaptive-mesh ask is tracked in the CHANGELOG as a roadmap item;
v1 supports a nonuniform *stated* grid, not automatic refinement.
