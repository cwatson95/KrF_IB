# docs/arf_data: the ArF data behind LIGHT_ArF

- `SOURCES_CHECKED.md`: every number LIGHT_ArF's species table and network
  use, and the page it was read from first-hand (2026-09-30).
- `RATES_RESEARCH.md`: the research reading of the local library that
  found them (a subagent's report, kept verbatim). It is a reading, not a
  verification: the table above is the check.
- `SCORING_DATA.md`: the ArF measurements found to score against (a
  subagent's report, kept verbatim), with the traps in the secondary
  sources (its section C).
- `fig_readings/`: the figures read into numbers, one CSV per figure, each
  with its source, page and reading method in its header.
- `ar_tail.py`: makes `v0_2/data/ar_tail_ratio.txt`, the two-temperature
  tail of the Ar rows, the way `ML/ml/tail.py` makes KrF's table. It reads
  ML's vendored Boltzmann solver and LoKI-MC's cross sections (set
  `KM_ROOT` to `Kinetic_Modeling/` when running from a copy) and writes
  only the table and `ar_tail_states.csv`.
- `ar_tail_states.csv`: the kinetic states behind the table (a results file:
  the generator rewrites it).
- `beam_drive.py`: the prototype ArF e-beam drive (Spencer-Fano yields for
  Ar/0.5 % F2, from REBUILD) that `ar_tail.py` and the scoring's 'beam'
  column use. It is an evidence tool, never part of LIGHT_ArF, which keeps
  LIGHT's own input.
- `score_arf.py`, `make_score_report.py`: LIGHT_ArF v0_2 and 1Dv0_2 against
  the ArF measurements. The results are `score_0d.json` and `score_1d.json`
  (each with the md5 of the chemistry it ran), the report `SCORE.md`; the
  1-D run outputs go to `runs/` (not kept).

From `Kinetic_Modeling/`:

    KM_ROOT=$PWD python3 LIGHT_ArF_ebeam/docs/arf_data/ar_tail.py            # the tail table
    python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py                        # 0-D score, about 13 min
    python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --only rokni_vs_ar     # one test, merged in
    ~/.venvs/f2sweep/bin/python LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --one-d   # 1-D, 64 cells, about 2.5 min
    python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --report               # SCORE.md
