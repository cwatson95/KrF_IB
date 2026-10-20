# SPECK/PAPER — the in-silico microwave-reheat article draft

Quick draft of a computational article:

> **Microwave reheating of electron-beam-pumped rare-gas–halide laser
> plasmas: a robust observable, a real but marginal channel**

Thesis, in the owner's framing: reheating is **feasible as an
observable** (fluorescence rises in 8/8 benchmarks on both rate tiers
and in the independent 1-D model), it **can be done and is not an
artifact** (the high-margin afterglow gain survives all 41
rate/literature brackets, a tier swap in the sibling cell, and timestep
refinement), and it is **very difficult** (five quantitative walls:
the ≤0.006-pp blended ceiling, ≤2-Td conductivity pinning, the
supply-chain inversion, the wall-plug inequality, and the one-rate
co-pump sign).

## Contents

- `mw_reheat_draft.tex` — REVTeX 4.2 (`aps,pra,preprint`); compiles
  clean, 12 pp. `[TODO: …]` markers (red in the PDF) are the
  placeholders: citations, affiliation line, funding text.
- `make_figs.py` — regenerates both figures FROM THE RECORDED RUN
  ARTIFACTS (sign_decidability / copump_verify / ib_feasibility JSONs +
  the SPECKLE co-pump reports). No hand-typed numbers.
- `mw_reheat_slides.tex` — the group-talk Beamer deck (Metropolis,
  16:9; 11 frames + backup); `make_slide_figs.py` regenerates its
  single-panel figures (`s_*.pdf`) from the same recorded JSONs.
- `fig_during_pulse.pdf`, `fig_robust_and_wallplug.pdf`,
  `fig_1d_verification.pdf`, `s_*.pdf`, `mw_reheat_slides.pdf`,
  `mw_reheat_draft.pdf` — generated; regenerate with the script +
  `pdflatex` (×2).

## Before circulating ANYWHERE

- **Export-control / data-rights review** (see the header comment in
  the .tex): some operating parameters trace to a controlled program
  document. The draft uses them only as modeling choices and never
  quotes or cites that document — review anyway.
- Funding/acknowledgment text: insert the program-approved statement
  VERBATIM (repo rule: never paraphrase funding text).
- The 1-D model is another author's code — obtain permission and a
  citation before naming it.

## Build

```bash
python3 SPECK/PAPER/make_figs.py       # from the Kinetic_Modeling root
cd SPECK/PAPER && pdflatex mw_reheat_draft.tex && pdflatex mw_reheat_draft.tex
```

Provenance: REBUILD CHANGELOG [0.1.20]–[0.1.29]; data of record listed
in the .tex header.
