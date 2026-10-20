# LIGHT_Kr: the lightweight pure-krypton model

LIGHT_Kr is the owner's copy of `LIGHT/` (2026-10-13), made to model pure
krypton instead of Kr/F2. The KrF files the copy came with were removed
the same day. v0_1 was made from LIGHT's v0_5, v0_2 from v0_1 on
2026-10-14, v0_3 from v0_2 and v0_4 from v0_3 on 2026-10-15, and v0_5
from v0_4 on 2026-10-16 (`CHANGELOG.md`). It holds:

- **`v0_5/`**: the current pure-krypton model. v0_4 with the electrons'
  collisions with the excited atoms and with Kr2* taken in: the stepwise
  ionization of the 5s levels and of the 5p rung and the ionization of
  Kr2* on Vriens & Smeets's calculation (bands 0 to twice it), Kr2*'s
  quenching and the mixing of its two states by electrons as bands at 0,
  the 5s levels' superelastic return and their 1s5 <-> 1s4 mixing by
  electrons; each 5s -> 5p pumping charged the energy of the 5p level it
  reaches, the rung's ten levels at the gas temperature with NIST's lines
  one by one and its relaxation's reverse; the two-temperature closure's
  tail table given an excited-population axis; Signorell's Kr2+. A kinetic
  table point with no solution is left out instead of stopping the run.
  `v0_5/README.md` has what changed and what it gives, `run_v0_5.py`
  launches it.
- **`v0_4/`**: v0_3 with its four 5s
  levels resolved: each is its own species and forms the excimer at its
  own measured rate (Marchal et al. 2018: 1s5 3.81e-32, 1s4 9.24e-32
  cm6/s, where v0_3's one Kr* used Ueda's 1.7e-32), is pumped to the 5p
  rung on its own cross sections, and is fed by the 5p lines on NIST's
  branching; the levels mix on krypton (Turner, Audouard). The excimer
  and the molecular ion redissociate by detailed balance on their
  measured wells (`thermo.py`); the fast excimer state, its quenching,
  the recombination's gas-temperature factor and the resonance lines are
  bands. The two-temperature closure reads a field-driven tail table from
  the model's own kinetic EEDF (`tail_field.py`), and the kinetic solver
  falls back to the growth term where a table point has no steady
  solution. `v0_4/README.md` has what changed and what it gives,
  `run_v0_4.py` launches it.
- **`v0_3/`**: v0_2 with its electron
  cross sections trained on krypton's measured swarm data by ML (`ML/`,
  "Pure krypton" in its README). Biagi's file is multiplied by seven
  factors, fitted to the eight drift-velocity and ionization-coefficient
  datasets that agree with one another, and the model is scored as each
  experiment measures it, through a Monte Carlo's ratios. Everything else
  is v0_2's. `v0_3/README.md` has the training and its scores, and
  `run_v0_3.py` launches it.
- **`v0_2/`**: v0_1 with the electron's collisions on ground-state
  krypton booked level by level. Every one of the 51 levels of Biagi's
  Magboltz 8.97 cross sections is a row at its own threshold (v0_1 read
  Biagi-v7.1's four groups), and the Hornbeck-Molnar ionization of the
  levels above Kr2+'s energy is a band at 0. `v0_2/README.md` has what it
  gives and why it was made; `run_v0_2.py` launches it.
- **`v0_1/`**: the first 0-D pure-krypton model. LIGHT's v0_5 with every
  fluorine species, row, band and parameter taken out, v0_5's krypton row
  for row, and two entries of its Jung et al. tables corrected to the
  paper. `v0_1/README.md` has what it gives at the reference point, what
  decides the answer and what it leaves out; `v0_1/SOURCES.md` the
  literature. `run_v0_1.py` launches it.
- **`docs/`**: each version's accuracy overview
  (`accuracy_overview_v0_5.pdf`, `_v0_4.pdf`, `_v0_3.pdf`, `_v0_2.pdf`,
  `_v0_1.pdf`) and system map (`system_map_v0_5.pdf`, ...), each generated
  from its model (`make_overview_v0_*.py`, `make_system_map_v0_*.py`). The
  runs go through `kr_docs_v0_5.py`, `kr_docs_v0_4.py`, `kr_docs_v0_3.py`,
  `kr_docs_v0_2.py` and
  `kr_docs.py` (v0_1), which run and cache the documents' runs
  (`docs/cache/`, keyed by the md5 of that version's code) and carry the
  e-beam drive the overview scores the model under.
- **`improvement_for_EEDF.txt`**: a note on the kinetic EEDF's steady-shape
  assumption in the first nanosecond and the afterglow, written for LIGHT
  v0_3. Since LIGHT's 2026-10-01 the pulse's EEDF is time-stepped; v0_4
  bounds the afterglow's steady shape (`Params.kin_afterglow_memory`:
  2e-5 of the reference point's 147-nm share).

The copy of LIGHT's v0_5 that v0_1 was made from (`v0_5/`, with
`run_v0_5.py`) is no longer in this directory (found absent on
2026-10-15). Its files were LIGHT's v0_5 byte for byte, less the KrF-only
ones, so `LIGHT/v0_5` holds the same model. v0_1's conversion tests skip
without the copy.

```bash
python3 LIGHT_Kr/run_v0_5.py --eedf two_t           # pure Kr, seconds
python3 LIGHT_Kr/run_v0_5.py                        # pure Kr, kinetic EEDF, ~20 min
python3 ML/run.py kr-train                          # retrain v0_3's factors (v0_4 and v0_5 keep them)
/usr/bin/python3 -m pytest LIGHT_Kr/v0_5/tests -q   # from Kinetic_Modeling
python3 LIGHT_Kr/docs/kr_docs_v0_5.py --workers 12       # the documents' 193 runs, ~4 h
/usr/bin/python3 LIGHT_Kr/docs/make_system_map_v0_5.py   # then latexmk -pdf in docs/
/usr/bin/python3 LIGHT_Kr/docs/make_overview_v0_5.py     # ~15 min on a current cache
```

At LIGHT's reference drive (1 atm, 300 kW/cm3 for 1 us) pure krypton
radiates 6.08 % of the input at 147 nm on v0_5's kinetic EEDF
(5.58 % on the two-temperature closure; v0_4 31.22 and 32.92 %,
v0_3 14.05 and 15.78 %, v0_2 14.20 and 16.00 %, v0_1 14.06 and 15.05 %).
The electrons now ionize 56.9 % of the excimers formed. The three
calculated rows are most of the change: with them at v0_4's 0 v0_5 gives
28.02 %, and without the largest of them, the ionization of Kr2*,
15.43 %. The propagated uncertainty is +-7.0
points (+-115 %), and its three largest terms are the
electrons' collisions with Kr2*: its ionization (a calculation), its
quenching and the mixing of its two states (nothing measures or
calculates them; at Kannari's estimate the share is 1.27
and 12.44 %) (`v0_5/README.md`).

## Where the KrF material is

Everything removed is in `LIGHT/` at the same path, byte for byte: the
older 0-D models `v0/` to `v0_4/`, the 1-D models `1Dv0/` to `1Dv0_3/` on
Fluid1D, `docs/` (the accuracy overviews, system maps, KrF data, scorers
and source papers), `IMPROVEMENT_PLAN.md`, and LIGHT's README and
CHANGELOG. LIGHT's Known-bad list in `LIGHT/CHANGELOG.md` applies to the
code here, with LIGHT_Kr's own additions in `CHANGELOG.md`.

The accuracy overviews score each version against krypton's swarm data
and three e-beam measurements (Morozov 2008, Eckstrom 1988, Saito 2003);
none of them tests the field-driven 5s <-> 5p cycle. v0_3 trains its
cross sections on the 8 swarm datasets that agree with one another (the
drift velocity 4.8 -> 2.8 % rms, the ionization coefficient at 50-300 Td
14.0 -> 9.9 % rms, scored as the experiments measure them). The tests
that neither the training nor the choice of its factors saw move with
it. The ionization coefficient below 50 Td goes from -44 % (Biagi's file,
as measured) to -36 %; every published cross-section set falls below
those measurements, which Bordage et al. 2013 ascribe to impurities. The
beam's partition moves toward Saito et al.'s measurement (W +16 -> +10 %,
excitations per ion pair +44 -> +24 %). v0_3's propagated uncertainty on
the reference point's 147-nm share is +-84 % (quadrature), three
quarters of it the spread of the Kr* + 2Kr measurements; the trained
cross sections' posterior adds +-0.04 points. Morozov's paper and
LIGHT's reading of it stay with the KrF data in `LIGHT/docs/krf_data/`
(`sources/morozov2008.pdf`, `FORMATION_LINKS_SEARCH.md` section 4), as
does Tahoma's pure-Kr absorption at 257.3 nm (`README.md` there), which
needs optics v0_1 does not carry.
