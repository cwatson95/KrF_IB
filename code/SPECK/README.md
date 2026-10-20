# SPECK — the plasma-chemistry simulation stack

SPECK is the two-layer PIC-MCC plasma-chemistry effort. The two packages
here are the two halves of one machine, split by TIMESCALE and coupled by
DATA — never by imports (gated in both test suites):

```
SPECK/
  PIC_MCC_3D/          the KINETIC half — electrons as macro-particles
                       (3-D electrostatic PIC + null-collision MCC on
                       LXCat cross sections). Resolves the
                       non-Maxwellian EEDF at ps scales in small boxes
                       (the chamber is unreachable for PIC — the
                       resolution budget in PIC_MCC_3D.txt §3).
  REBUILD/             the CHEMISTRY half — the Orestes-class engine:
                       169-reaction / 35-species KrF-ArF network with
                       provenance tiers, energy ledger, cavity photons,
                       ASE, anchors. Rate equations at ns scales, 0-D
                       or a coarse 3-D device grid, up to chamber
                       scale. Predictive-only: anchors score, never
                       steer.

  the coupling (today, offline):
    PIC_MCC_3D  ──(tools/export_kinetic.py: swarm-measured rate
                   coefficients → kinetic_tables/*.npz)──►  REBUILD's
                   KINETIC rate tier
  the coupling (the next build, online):
    a tools/ coupler that sub-cycles REBUILD's compiled network inside
    the PIC loop of a small box — MCC event tallies in, updated gas
    composition back. tools/ may import both packages; the engines
    still never import each other.

  PIC_MCC_3D.txt              the PIC roadmap (the nine gaps, the
                              resolution wall, milestones M1–M10)
  ORESTES_CLASS_REBUILD.txt   the engine roadmap
  missing_modules.txt         observables the models cannot yet predict
  MESH.txt                    the adaptive-mesh ask (static rung closed
                              in REBUILD 0.1.14-era; dynamic open)
  COMMERCIAL_AUDIT.md         what stops the results from carrying
                              commercial-grade confidence (2026-09-03):
                              convergence measured, solver robustness,
                              model form, accuracy envelope, guards
  rebuild_out/                generated artifacts (gitignored):
                              movie.gif = the GLOW camera view,
                              movie_particles/_fields.gif, viewer.html,
                              viz_data.npz per preset
  picmcc_out/                 PIC demos: Langmuir slosh, two-stream
                              instability (real macro-particles)
```

Run everything from the `Kinetic_Modeling/` root:

```bash
python3 -m pytest SPECK/REBUILD/tests SPECK/PIC_MCC_3D/tests -q
python3 SPECK/REBUILD/run.py --preset krf_beam_3d
python3 SPECK/REBUILD/tools/visualize_3d.py --preset krf_beam_3d_bragg
python3 SPECK/PIC_MCC_3D/tools/visualize_particles.py --demo twostream
```

Two-interpreter rule (PLAN (68)) applies to both halves: pytest under
`/usr/bin/python3`, import + run smoke under Homebrew `python3` with
`-W error::DeprecationWarning`, before claiming green.

Read each package's `CHANGELOG.md` (Known-bad first) before editing
anything, and each `status.html` for the live board. `HIBACHI/`
(the e-beam driver chain, one level up) imports REBUILD as a consumer
and knows to look here.

`EXPERIMENTS.md` (this directory) is the experiment matrix: known
experiments and the analytic canon as verification targets, tiered by
readiness (live anchors / ready / named small builds / out-of-lane),
each mapped to the geometry it needs and the modules it verifies,
provenance class stated per row.

`tools/` is the ONLY place allowed to import both halves (the engines
never import each other): `couple0d.py` is the PIC↔chemistry coupler —
the data contract (n_e, T_e, T_g, E/N, n_ph clamps → heavy-species
inventory) plus a clamped-BDF sub-cycling core, gated in the
reproduce-limit against a REBUILD afterglow (worst species 4% at 1 ns
clamps). Its suite: `pytest SPECK/tools/tests` (runs with the two
package suites from the Kinetic_Modeling root).
