# SPECKLE — Spatial Plasma and Excimer Chemistry Kinetic Light-Engine

A **fully three-dimensional, fully kinetic** particle-in-cell / Monte-Carlo-collision
(PIC-MCC) simulation of an excimer-laser discharge: the plasma, the excimer
chemistry, and the laser light are all produced **from particles**, not from a
fluid closure. Where the sibling engines model the electron gas as a temperature
reading a swarm table, SPECKLE carries the electrons as macro-particles and reads
the reaction rates off the actual velocity distribution.

> Status: the electrostatic PIC-MCC **plasma core is built and verified**
> (v0.8.4). The **excimer chemistry** and **light-engine** subsystems that turn a
> discharge code into a *laser* code are the current frontier. See
> [The frontier](#the-frontier--what-makes-it-a-laser).

This code currently lives at `SPECK/PIC_MCC_3D/`; the plan of record is
`SPECK/PIC_MCC_3D.txt` (milestones M0–M10, each with a falsification gate). This
README is the human-facing entry point to that plan. **On the name, see
[Naming](#naming--how-speckle-relates-to-the-other-codes) — "SPECKLE" is already
used for a 1-D fluid study, and this is a different, kinetic thing.**

---

## The name, decoded

Each word of the backronym is a subsystem, with its honest build state:

| Word | Subsystem | Built from | Status |
|---|---|---|---|
| **S**patial **P**lasma | 3-D electrostatic PIC: charge-conserving deposition, FFT+tridiagonal Poisson, Boris push, mobile multi-species ions | `poisson3d`, `es_field`, `es_sim`, `species3d` | **working, gated** |
| **E**xcimer **C**hemistry | the KrF\*/ArF\* formation–quenching network run per PIC cell, driven by the kinetic electron/ion populations + null-collision MCC | `mcc3d`, `gas`, `cross_sections` (e-neutral done); reaction-network coupling **frontier** | core done, coupling open |
| **K**inetic | the defining choice: real macro-particles with a self-consistent EEDF — no fluid T<sub>e</sub>, no swarm-table lookup | the whole particle store + MCC | **by design** |
| **L**ight-**E**ngine | 3-D laser photon field: cavity, ASE/parasitics, and photon transport from the excimer upper states | **frontier** (not yet in the PIC core) | open |

The through-line: **kinetic**. SPECKLE's reason to exist is that the electron
energy distribution is *computed from the particles*, so nothing downstream — not
the ionization rate, not the excimer formation rate, not the gain — depends on a
fluid-temperature assumption.

---

## Why a fully-kinetic engine (the gap it fills)

Three excimer codes already exist in this program, and each is missing something
the other two have:

- **REBUILD** (SPECK / SPECKLE3D) — full excimer chemistry, e-beam physics, a
  laser light-engine, an energy ledger, and device anchors — but its electron gas
  is a **single fluid temperature** on precomputed swarm tables. Its documented
  residuals (F⁻ ×0.19, Ar\*\* ×0.06) sit exactly where the EEDF *shape* matters.
- **Fluid1D** (1DModelKrF) — a genuine in-loop two-term Boltzmann EEDF and real
  1-D transport/sheaths — but **discharge-only, 1-D, no beam, no ASE, no gas
  heating**, and no experimental anchors.
- **PIC_MCC_3D** — the kinetic plasma done right in 3-D — but today it produces
  the *plasma*, and the excimer chemistry and laser optics are read from it as
  data (the KINETIC tier REBUILD consumes), not run inside it.

SPECKLE is the union: the kinetic 3-D plasma of PIC_MCC_3D, **plus** the excimer
chemistry and the light-engine coupled in, so a single code makes the electrons,
the excimers, and the photons self-consistently. It is the first-principles
reference the fluid engines are meant to be scored against — the KINETIC tier,
grown up into a laser.

---

## What already runs (the honest inventory)

The plan is a **promotion, not a rewrite**: every missing piece existed first in
the `pic1d`/`pic2d` lineage and is lifted into the 3-D core. Done and gated
(under both interpreters, per the two-interpreter rule):

| Milestone | Subsystem | Verification (recorded) |
|---|---|---|
| **M1** | 3-D electrostatic Poisson (`poisson3d`, `es_field`) | MMS slope 1.97/2.03; FFT-vs-LU 1e-15; cold-plasma f/f<sub>pe</sub> to <1% |
| **M2** | mobile, subcycled multi-species particles (`species3d`, `es_sim`) | cold-plasma oscillation 0.99158 vs predicted 0.99145; ionization births real ions |
| **M3a** | electron-neutral MCC (`mcc3d`, `gas`, `cross_sections`) | Vahedi–Surendra null-collision; per-process ledger closure 2e-16; matches the EM instrument's operator on a shared RNG stream |
| **M5** | cathode secondary emission — the γ process (`emission3d`) | self-sustaining discharge; G4.1 met; Ar Paschen vs Raizer |
| **M9** | e-beam injection + preionization (`beam3d`) | fast primary via Bethe/Møller; W-value emergent by two routes to 0.4% |

The electromagnetic sibling it grew out of — `ARPAE/pic3d` (Yee FDTD + CPML, a
wave-plasma instrument) — stays as the EM path (`field_mode="em"`), byte-identical
and guarded by its 41 gates. Current version: **v0.8.4** (2026-08-14).

---

## The construction — the PIC-MCC cycle

One timestep of the electrostatic core, and the module that owns each stage:

```
   ┌─────────────────────────────────────────────────────────────┐
   │  1. PUSH      Boris advance of every macro-particle          │  species3d
   │               (electrons every dt; ions subcycled k_i·dt)    │
   │      ↓                                                        │
   │  2. WEIGHT    charge-conserving CIC deposition  ρ(x)         │  es_field
   │      ↓                                                        │
   │  3. SOLVE     ∇²φ = −ρ/ε₀   via 2-D FFT in (x,y)             │  poisson3d
   │               + independent tridiagonal solve in z           │
   │               → O(N log N), exact, no iteration              │
   │      ↓                                                        │
   │  4. GATHER    E = −∇φ, gathered with the SAME shape          │  es_field
   │               function as the deposit (energy-conserving)    │
   │      ↓                                                        │
   │  5. COLLIDE   null-collision MCC: elastic / excitation /     │  mcc3d
   │               ionization / attachment; ionization makes ions │
   │      ↓                                                        │
   │  6. BOUNDARY  electrodes + external circuit + surface charge;│  (M4, open)
   │               secondary emission at the cathode              │  emission3d
   └─────────────────────────────────────────────────────────────┘
```

The **design decision that makes 3-D affordable** is stage 3: the periodic-in-(x,y),
Dirichlet-in-z geometry lets a 2-D FFT reduce the 7-point stencil to one
tridiagonal solve per (k<sub>x</sub>,k<sub>y</sub>) mode — O(N log N), round-off
exact, no preconditioner, no new dependency. Ion subcycling (stage 1) keeps the
ion timescale affordable without an implicit mover.

---

## The frontier — what makes it a laser

The plasma core is a **discharge** code. Three things turn it into an **excimer
laser** code, and they are what's being built now:

1. **Excimer chemistry coupling.** Run the KrF\*/ArF\* formation–quenching network
   (the 169-row network that today lives in REBUILD) inside each PIC cell, sourced
   by the MCC-produced electrons and ions. The kinetic EEDF drives the excimer
   kinetics directly — the payoff over the fluid engines.
2. **The light-engine.** A 3-D laser photon field: cavity build-up, ASE /
   parasitic transport, and stimulated emission from the excimer upper states —
   the "L·E" of the name. (Note: `emission3d` is *secondary-electron* emission,
   not light; the optical engine is new.)
3. **The remaining plasma modules**, in the forced dependency order from the plan:
   **M3b** ion-neutral collisions (elastic + charge exchange), **M4**
   electrodes / external circuit / surface charge, **M6** Coulomb collisions
   (e–e / e–i, Takizuka–Abe), **M7** an implicit / energy-conserving mover, **M8**
   particle merge/split. Each already exists in `pic1d`/`pic2d`.

---

## Policy and verification (non-negotiable)

- **Predictive-only.** Every cross section, rate, and branching ratio comes from
  theory or a published measurement *of that process*. Analytic overlays and
  published benchmarks **score** the code; they never steer a number. A gate that
  fails is reported as a failure.
- **The verification ladder** (named, published targets, no tunable in sight):
  MMS order-of-accuracy (Riva 2017), cold-plasma oscillation at ω<sub>pe</sub>,
  Bohm–Gross / Landau dispersion, two-stream growth rate, Paschen curve
  (vs Raizer), Child–Langmuir, and the Turner et al. RF-discharge benchmark once
  M5 is in.
- **Two-interpreter rule.** No milestone is "done" until its gate passes under
  both interpreters and the changelog carries the entry.
- **Energy ledgers everywhere.** The MCC closes to 2e-16, Gauss's law holds through
  ionization/attachment/wall events — the diagnostics are the product, not an
  afterthought.

---

## Naming — how SPECKLE relates to the other codes

**Heads-up on the name, because it is reused in this program:**

- **`SPECKLE/`** (existing) is a **1-D microwave-reheat study** — one driver on
  Nathaniel's *fluid* 1DModelKrF. Fluid electrons, 1-D, reheat-only.
- **`SPECKLE3D/`** (existing) is a **3-D device-grid study on REBUILD** — *fluid*
  electrons, independent-reactor cells, plus the cold-plasma visualization.
- **This SPECKLE** is the **kinetic PIC-MCC engine** — real particles, a computed
  EEDF, in 3-D. It is a *different kind of thing* from the two above (kinetic vs
  fluid), and it currently ships as `SPECK/PIC_MCC_3D/`.

**Decision for the owner:** the three "SPECKLE" names now span two fluid studies
and one kinetic engine. Options: (a) keep the kinetic engine as `PIC_MCC_3D` and
let this README carry the SPECKLE identity/vision; (b) rename the fluid studies
(e.g. `SPECKLE-1D`, `SPECKLE3D-fluid`) and promote this to the `SPECKLE` name; or
(c) pick a distinct name for one side. This README does **not** rename anything —
it names the vision and flags the collision for you to resolve.

Relationship to the fluid engines (see the *"Two Engines, Four Wrappers"* guide):
SPECKLE is the first-principles reference REBUILD and Fluid1D are scored against.
At the full vision, the "KINETIC tier" REBUILD reads from PIC as data becomes
SPECKLE running the chemistry and light itself — the same physics, no longer
handed across a coupler.

---

## Running

Python 3 with `numpy` (a numba fast path is optional; the reference stays numpy).
From the `SPECK/` root:

```bash
/usr/bin/python3 -m pytest PIC_MCC_3D/tests -q     # the gated suite
python3 PIC_MCC_3D/verify.py                        # the live verification board
```

The plan of record, milestone gates, the nine-gap dependency order, and the
resolution budget for a full 3-D device run are all in `SPECK/PIC_MCC_3D.txt`.
History is in `SPECK/PIC_MCC_3D/CHANGELOG.md`.
