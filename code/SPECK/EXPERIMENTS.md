# SPECK experiment matrix — known experiments as verification targets

Purpose: one table of real, published experiments (plus the analytic canon) that
SPECK either recreates today or can take on next, each mapped to the **geometry**
it needs and the **modules** it verifies. The improvement path this encodes:
build one geometry per experiment class, and every module gets exercised against
something measured.

Covenant (PLAN (26), carried from REBUILD's scorecard): anchors **score, they
never steer**. No experiment on this list may feed a number back into an engine.

**Provenance classes** — stated per row, because they mean different things:

- **blind** — the engine's inputs do not descend from the measured result; the
  score is a genuine prediction (the REBUILD device-anchor lane).
- **consistency** — the compared data shares provenance with the engine's inputs
  (e.g. swarm data the LXCat σ-sets were built against; DR rates measured in the
  very afterglow experiments being recreated). Still a real module verification
  — "does the chain reproduce the data its inputs came from" — but never to be
  presented as a blind prediction.
- **analytic** — closed-form theory target; for discretized schemes the honest
  target is the *discrete* prediction (see `PIC_MCC_3D/MODULES/theory.py`).

Status legend: ✅ live now (ran/verified 2026-08-14) · 🟢 ready — existing code,
add a preset/scorecard row/demo · 🟡 needs a named small build (named in the row)
· 🔴 needs the PIC↔chemistry coupler (the named next build) · ⛔ out of lane.

## The geometry ladder

| ID | Geometry | Exists? | Where |
|----|----------|---------|-------|
| G0 | 0-D reactor cell | ✅ | REBUILD `krf_beam_0d`, `arf_chen_0d` (nx=ny=nz=1) |
| G1 | Periodic cube | ✅ | PIC_MCC_3D periodic mode (Langmuir/two-stream demos) |
| G2 | Plane-parallel gap: periodic (x,y), biased electrodes at z=0/Lz | ✅ solver+walls | `poisson3d` `periodic_xy_dirichlet_z` ("THE DEVICE GEOMETRY"), settable `phi_lo`/`phi_hi`, absorbing electrodes tallying `wall_charge`/`wall_energy_J` |
| G3 | Coarse 3-D device grid (independent-reactor cells) | ✅ small/mid | REBUILD `krf_beam_3d` (Electra-class, 240 cells), `krf_beam_3d_bragg` (100 kV Bragg slab); **to add:** large-aperture Nike/LAM-class preset |
| G4 | Sub-chamber coupled box (PIC + sub-cycled chemistry) | ⬜ future | the `tools/` coupler — the named next build |

## Tier A — live today (the existing anchor set)

### A-chem: device anchors (REBUILD + HIBACHI)

| # | Experiment | Compared | Modules verified | Geometry | Provenance | Status / score |
|---|-----------|----------|------------------|----------|------------|----------------|
| A1 | Kannari et al. 1985 e-beam KrF amplifier (Table I; 4×4×44 cm, 66 J deposited matched by construction) | laser output vs output coupling (2 pts) | chemistry network + cavity + deposition | G3 (`kannari_device_rows`) | blind | ✅ ×0.83 robust pt (47% margin) / ×0.20 fragile pt (self-labels ×9-amplified, 11% margin) |
| A2 | NRL Electra ArF quadruple-pass oscillator — Wolford 2020 HEDP 36, 100801 (Table 1 + §2) | output energy / populations | ArF chain + cavity | G3 (`electra_arf_rows`) | blind | ✅ ×0.84–×1.76 band |
| A3 | NRL Electra KrF population benchmark | KrF\* formation populations | formation chain diagnosis | G3 (`electra_population_rows`) | blind | ✅ recorded finding: bulk ran ×4 cold — the cross-device ×0.5–0.6 deficit diagnosis |
| A4 | Chen 2025 ArF **discharge** device (Fig. 6 waveforms; Table II network verbatim-confirmed) | n_e, n_ph, peak P_L, E/N class, deposited energy | discharge **circuit drive** (v0.1.5) + ArF chain + self-consistent breakdown | G0 (`arf_chen_circuit`, `chen_figure6_rows`) | blind | ✅ waveform/E-N class + deposition reproduced; n_e/pulse low by the recorded ×3–5 per-joule formation gap |
| A5 | W-values, ICRU 31 / Combecher class (Ar 26.4, Kr 24.2, Ne 36.6 eV/ip) | emergent W per ion pair | beam slow-down / deposition partition | G0 | blind (emergent, computed live) | ✅ |
| A6 | NIST ESTAR-class stopping powers @ ~600 keV (He/Ne/Ar/Kr/F2) | mass stopping | stopping module | G0 | blind | ✅ |
| A7 | Electra hibachi drive chain (NAS study, s14 numbers) | deposited-in-gas / beam fraction | HIBACHI energy ledger (emitter→foil→gas) | HIBACHI chain | blind | ✅ 0.649 vs 0.75 measured → ×0.86 |

### A-kin: analytic canon (PIC_MCC_3D, gated in the suite)

| # | Benchmark | Compared | Modules verified | Geometry | Provenance | Status |
|---|-----------|----------|------------------|----------|------------|--------|
| A8 | Cold Langmuir oscillation | ω vs the **discrete** dispersion (CIC S², spectral K², centred k_g — `theory.py`) | deposit/gather + Poisson + leapfrog | G1 | analytic | ✅ gated |
| A9 | Two-stream instability | growth of per-beam bunching vs closed-form γ=ω_b/2 (±20%) | collective instability path | G1 | analytic | ✅ gated |
| A10 | LoKI-MC oracle match | swarm/EEDF chain ×1.6 agreement | σ-chain end-to-end | G1 | consistency (shared isotropic-on-MT basis — stated in `cross_sections.py`) | ✅ recorded |

## Tier B — ready now: new geometry/preset + a scorecard row or gated demo

| # | Experiment | Compared | Modules verified | Geometry | Provenance | What to build |
|---|-----------|----------|------------------|----------|------------|---------------|
| B1 | Ar / Kr drift velocity + characteristic energy vs E/N (Phelps-era swarm compilations, the LXCat sets' own reference data; Ar/Kr/F2/He σ-files ship in `MODULES/lxcat/`) | v_d(E/N), D/μ | null-collision MCC + σ tables + pusher under uniform E | G2 drift tube: bias electrodes, seed via `beam3d` preionization, drift, absorb | consistency | preset + measurement rows |
| B2 | Townsend first ionization coefficient α/N vs E/N — Kruithof 1940 (Physica 7: Ne/Ar/Kr/Xe) | ionization events per drift length | ionization spawning (real e⁻+ion macro-particles, `mcc3d`) + energy balance | G2 drift tube, count births vs z | consistency | same preset, birth bookkeeping |
| B3 | Bohm sheath / floating-wall potential (Bohm 1949 probe theory) | φ_float ≈ (T_e/2e)·ln(m_i/2πm_e); sheath scale ~ λ_D | bounded Poisson + absorbing walls + ion (Ar⁺) dynamics | G2 grounded slab | analytic | demo + gate |
| B4 | Debye shielding of a biased electrode | shielded φ(z) vs analytic profile | Poisson + thermal equilibrium | G2, small bias | analytic | demo + gate |
| B5 | Landau damping of a warm Langmuir wave | γ_L(kλ_D) | velocity-space fidelity (quiet start) | G1 | analytic | demo + gate |
| B6 | Single-particle E×B / gyro drifts (Boris takes B already — `species3d.push_boris`) | drift speeds, μ conservation | Boris rotation correctness | G1 + static B | analytic | demo + gate |
| B7 | Microwave-afterglow dissociative recombination, Shiu–Biondi class (Ar₂⁺/Kr₂⁺, PRA 1977–78) | 1/n_e(t) slope → α_DR(T_e) | REBUILD recombination chain in decay | G0 afterglow (pulse→decay; the movie machinery) | consistency (DR rates in the network descend from these experiments) | scorecard row |
| B8 | F₂ attachment-dominated decay (Christophorou-school attachment data) | n_e e-fold time vs [F₂] | attachment rows | G0 afterglow, F₂ scan | consistency | scorecard row |
| B9 | Electron depth–dose (Bragg) profile vs Tabata–Ito semi-empirical fits | dose(z) shape + practical range | deposition shape module | G3 Bragg slab (`krf_beam_3d_bragg` exists) | data-fit target | comparison rows |
| B10 | NRL Nike large-aperture KrF amplifier (Obenschain-era facility papers) | output energy; ASE fraction vs aperture | **ASE geometry scaling** — the module only large apertures stress | G3 large-aperture preset (the S/M/L ladder: Kannari → Electra → Nike) | blind | sourcing pass in `books/` for deposition+geometry numbers, then preset + rows |
| B11 | LANL Aurora / LAM (Rosocha-era) | same as B10 at ~1 m class | same | G3 | blind | same sourcing caveat |

## Tier C — the one-build-away set (build status 2026-08-14: all five begun; C1/C2 done)

| # | Experiment | Status 2026-08-14 | Geometry | Result |
|---|-----------|-------------------|----------|--------|
| C1 | Paschen curve V_b(pd) in Ar vs Raizer (printed Ar/iron anchor: 265 V @ 1.5 Torr·cm, p. 135) | ✅ built + gated (`PIC_MCC_3D/tools/paschen.py`; M5 `emission3d` γ process landed) | G2 gap | Criterion = Raizer Eq. (7.1), γ[e^{αd}−1]=1, with M0 measured by direct avalanche counting (γ = stated cathode parameter). Two-length α/p (the source's own method, §4.1.4) = **×0.79–1.02 of the Table 4.1 fit over 100–300 V/(cm·Torr)** (×0.8 near the Paschen minimum, ×1.0 at the window edges; closer still to the inert-gas form (4.6)) — recorded; energy-sharing sensitivity only 6%. **Curve**: γ_iron=1.5e-3 minimum at pd≈1.5 (352 V, ×1.33 the printed 265 V — the α deficit through e^{αd}); γ=1e-2 passes THROUGH the measured point (264 V @ 1.5). Dynamic γ-loop demo: grew at 1.3 V_b, died at 0.7 V_b. Artifacts `picmcc_out/paschen/` |
| C2 | Child–Langmuir diode, J ∝ V^{3/2}/d² (plan gate G4.1) | ✅ **G4.1 met** (`PIC_MCC_3D/tools/diode.py`) | G2 gap | Emit 3× J_CL → transmitted **J/J_CL = 1.023** at 250 V (virtual cathode −0.24 V @ 39 µm); 1.072→1.010 over 100→400 V — the stated O(√(T/V)) thermal bias shrinking exactly as it should; 1.01–1.04 across nz 33→257 |
| C3 | Turner et al. 2013 CCP benchmark (Phys. Plasmas 20, 013507, DOI 10.1063/1.4775084; He cases) | 🟡 begun, **not claimed** (`PIC_MCC_3D/tools/rf_ccp.py`) | G2 gap + RF | φ(t) electrode waveform landed, leapfrog gated on the closed-form RF trajectory (<1% of amplitude); case table transcribed from a reproduction repo (provenance stated in the tool — verify against the paper before claiming G5.1); case-1 smoke closes the loop. Blockers: first-hand case spec, M6 Coulomb |
| C4 | PIC↔chemistry coupler (`SPECK/tools/couple0d.py`) — the streamer-chemistry prerequisite | 🟢 begun + **reproduce-limit gated** | G0 core (G4 later) | Clamped-BDF sub-cycling of REBUILD's compiled network under PIC-supplied (n_e, T_e, T_g, E/N, n_ph). Parity vs the engine's own afterglow: **19 species, worst 4.0%** at 1 ns clamps (16% at 2 ns — clamp resolution is the error knob, measured). Demo: 38 Torr Ar/Kr/F₂ gap → chemistry handoff. Return leg (gas update into the MCC) = named next step |
| C5 | E-beam-sustained conduction — the non-self-sustaining i–V (Raizer Fig. 7.1 regimes; EBSD-class) | 🟢 kinetic skeleton built (`PIC_MCC_3D/tools/ebsd.py`) | G2 gap + beam | BeamInjector (kinetic pair mode, W emerges) + bias scan: **5.6→49.8 A/m² over 20→240 V**, ×5.8 amplification 60→240 V (consistent with e^{αd} at that E/p), ion space charge ≤10% and reported, structurally non-self-sustaining (no γ registered). Coupled EBSD rides on C4 |

## ⛔ Out of lane (structural, stated)

- Magnetized-confinement / electromagnetic devices (tokamak, helicon, ICP wave
  heating): the field solve is electrostatic (`poisson3d` names why the Yee/FDTD
  path cannot substitute); Boris accepts only a prescribed static B.
- Full-chamber kinetic simulation at device density: the λ_D resolution wall
  (`PIC_MCC_3D.txt` §3) — sub-chamber boxes only, by design.
- Resolving the relativistic primary beam as PIC particles: excluded on purpose
  (`beam3d` header — γ=2.17 at 600 keV, above the 2 keV σ-table ceiling); the
  primary is a source model, and that stays.
- Discharge-pumped **KrF** device anchors: the circuit drive exists (A4), but
  the discharge lane is currently the ArF/Chen network; a KrF discharge anchor
  would need its own vetted source first.

## Adding an anchor — the recipe (unchanged from REBUILD's scorecard)

1. A construction docstring: which published numbers are INPUT ("matched by
   construction") and which are OUTPUT (scored, never fed back).
2. Run LIVE in the scorecard/demo — no hand-recorded result text.
3. Misses stated per row (×ratio), with conditioning gauges where the operating
   point amplifies (g_th margin, dlnE/dlng0 — the ×0.20 Kannari point self-labels).
4. Provenance class stated (blind / consistency / analytic).
5. Omissions named. "Could-not-check" is an acceptable verdict; a rubber-stamp is not.

---
*Status audit 2026-08-14: ✅ rows re-verified against the live code this date
(Kannari rows re-run; solver/walls/σ-roots/Chen block read from source). 🟢/🟡
assessments name their blocking module from the code's own stated exclusions.
Citations are given at the confidence level the repo records them (scorecard
source strings verbatim); B10/B11 explicitly need a sourcing pass before use.*
