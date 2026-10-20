# LIGHT_ArF_ebeam: LIGHT_ArF with an electron-beam drive

LIGHT_ArF_ebeam is the owner's copy of `LIGHT_ArF/` (2026-09-30), made to
add an electron-beam drive: "Look through Kinetic_Modeling/SPECK/REBUILD/
and see if adding the electron beam drive in LIGHT_ArF_ebeam/ improves our
models and if it agrees with Morozov as well". The answer, in
`docs/arf_data/SCORE.md` and under Accuracy below: it improves the model
overall, and it does not reach Morozov. It holds LIGHT_ArF's two models:

- **v0_2** (`v0_2/`): the 0-D model, one well-mixed box of argon and
  fluorine, with two ways to deposit the given power: LIGHT's own input
  into the electron gas (`--drive power`, the default) and an electron beam
  (`--drive beam`, `v0_2/beam.py`). `v0_2/README.md`.
- **1Dv0_2** (`1Dv0_2/`): the 1-D model, v0_2's chemistry in every cell of
  a gap on Nathaniel's Fluid1D, on the same two drives: the beam's products
  enter every cell of its stage. `1Dv0_2/README.md`.

LIGHT_ArF's history is in `LIGHT_ArF/CHANGELOG.md`. This copy's starts with
the drive (`CHANGELOG.md`, entry (a), and entry (b) for the 1-D model),
under the same Known-bad list.

```bash
python3 LIGHT_ArF_ebeam/run_v0_2.py                    # the 0-D reference point
python3 LIGHT_ArF_ebeam/run_v0_2.py --drive beam       # the same under the e-beam
python3 LIGHT_ArF_ebeam/run_v0_2.py --bands            # each band's two ends
~/.venvs/f2sweep/bin/python LIGHT_ArF_ebeam/run_1Dv0_2.py --compare   # the 1-D model (numba)
~/.venvs/f2sweep/bin/python LIGHT_ArF_ebeam/run_1Dv0_2.py --drive beam --compare   # 1-D, beam
python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py     # the 0-D score (about 14 min)
python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --report   # SCORE.md
python3 LIGHT_ArF_ebeam/docs/make_overview.py          # the overview's numbers and figures
FLUID1D_SRC=.../nathaniel/1DModelKrF/src ~/.venvs/f2sweep/bin/python \
    LIGHT_ArF_ebeam/docs/make_overview.py --runs       # + its runs, ~20 min
python3 LIGHT_ArF_ebeam/docs/make_energy_ladder.py     # the energy ladder
python3 LIGHT_ArF_ebeam/docs/make_system_map.py        # the system map's numbers, ~7 min
(cd LIGHT_ArF_ebeam/docs && /opt/homebrew/bin/latexmk -pdf system_map.tex)          # the map page
(cd LIGHT_ArF_ebeam/docs && /opt/homebrew/bin/latexmk -pdf accuracy_overview.tex)   # the PDF
python3 -m pytest LIGHT_ArF_ebeam/v0_2/tests -q        # the 0-D tests
```

## The electron-beam drive

LIGHT's own input puts the given power into the electron gas, which then
excites and ionizes the gas as a discharge would: in pure argon it makes
9.5 excited atoms per ion at 153 eV per ion. A beam's fast electrons make
ions and excited atoms directly and leave each secondary electron a few
eV. `v0_2/beam.py` is that drive, adapted from
`SPECK/REBUILD/MODULES/beam.py`:

- REBUILD's Bethe/Moller stopping power (a printed beam voltage and current
  density give the deposited power) and its Spencer-Fano cascade, run on
  the model's own cross sections (the IST-Lisbon Ar and SIGLO F2 files,
  md5-identical to REBUILD's);
- the secondary's energy from the measured Opal-Peterson-Beaty spectrum
  (Ebar(Ar) = 10.0 eV, LoKI-MC's database), and the yields per pair as the
  cascade's asymptotic slopes;
- per ion pair at 0.5 % F2: W 25.2 eV (4.5 % under the measured 26.4), 0.39
  excited atoms, 4.4 eV to the electron gas, 0.55 eV of prompt heat.

REBUILD's own rule, equal sharing read at one source energy, makes the
yields oscillate with that energy (0.31 to 1.28 excited atoms per ion
between 100 and 900 eV). LIGHT_ArF's prototype drive (`docs/arf_data/
beam_drive.py`, which only the tail-table generator still uses) sat on one
point of it: W 30.85 eV and 0.85 excited atoms per pair.

## The ArF network

Every rate cites a measurement whose page was read first-hand
(`docs/arf_data/SOURCES_CHECKED.md`: many are review tables quoting the
measurement, and F + F + Ar is KrF's reading of Baulch 1981) or is a band:
a measured rate known only as a spread keeps its measurement as the
default; a rate with no measurement defaults to 0. What differs from KrF:

- **One excimer state.** No ArF C position, B-C mixing rate or
  equilibrium is measured (the calculated splitting is about 0.25 eV), so
  ArF* is ArF(B) alone.
- **Stern-Volmer products, not rate constants.** The ArF(B) lifetime has
  never been measured (4.2 ns is Dunning & Hay's calculation). Its
  quenching by F2, Ar, 2Ar and electrons is measured as k tau, and the
  rows divide by tau_B.
- **Ar* + F2 makes ArF* with yield 0.5** (Setser's measurement), not the
  1.0 every published ArF model books (the band's top).
- **Each IST-Lisbon Ar level is its own row** (37 excitations), the 4s ones
  making Ar*, the higher ones cascading to it.

`v0_2/README.md` has the whole network and the reference point; the
research reading behind it is `docs/arf_data/RATES_RESEARCH.md`, and the
measurements it is scored on `docs/arf_data/SCORING_DATA.md`.

## Accuracy

`docs/accuracy_overview.pdf` is the accuracy overview: every score, the
reference point and its bands, the 1-D model, the energy ladder, the system
map, the reaction rates and the data scored. Its numbers are written by
`docs/make_overview.py`, the ladder by `docs/make_energy_ladder.py` and the
map's numbers by `docs/make_system_map.py`.
`docs/arf_data/SCORE.md` holds every point, beside LIGHT_ArF's record
(`score_0d_LIGHT_ArF.json`). No rate was set from these data. Each 0-D test
runs on LIGHT's own input, on the same input with a Maxwellian EEDF (a
sensitivity) and under the electron-beam drive. The seeded SRI starts take
the drive's excited atoms per ion. The 1-D model (64 cells) runs the device
points and Rokni's pressure curve on both drives, at five of the pressures
read off his figure (0.243 to 4.003 atm). The numbers are model / measured:
geometric means, with the rms of the logs in brackets where there are three
points or more.

| test | n | own input | own, Maxwellian | e-beam drive | LIGHT_ArF prototype beam | 1-D centre | 1-D centre, beam |
|---|---|---|---|---|---|---|---|
| Electra formation efficiency g0 Isat / P (1.09, 0.90 MW/cm3) | 2 | 0.64, 0.55 | 0.71, 0.61 | 0.84, 0.72 | 0.81, 0.70 | 0.64, 0.55 | 0.84, 0.72 |
| the same at 1.09 MW/cm3 with the ArF* yield at 1.0 (sensitivity) | 1 | 1.19 | | 1.08 | 1.12 | | |
| Electra gain at 1.09 and 1.36 atm, ratio to 0.82 atm | 2 | 1.03, 1.31 | | 1.03, 1.32 | 1.03, 1.32 | 1.03, 1.31 | 1.03, 1.32 |
| Keio formation efficiency at its printed 0.51 MW/cm3 | 1 | 0.86 | 0.98 | 1.14 | 1.10 | 0.86 | 1.14 |
| Sandia formation efficiency (abstract only) | 1 | 1.38 | 1.55 | 1.82 | 1.76 | 1.38 | 1.82 |
| Rokni 193 nm against Ar (shape) | 10 | 0.97 (0.15) | 1.05 (0.07) | 0.99 (0.10) | 1.00 (0.08) | 1.00 (0.24), 4 points | 1.03 (0.14), 4 points |
| Rokni 193 nm for a x25 current change | 4 | 1.67 (0.58) | 0.96 (0.09) | 0.97 (0.05) | 0.97 (0.06) | | |
| Jancaitis peak ArF* and Ar2F* against current (shape) | 14 | 1.06 (0.28) | 1.18 (0.32) | 1.17 (0.35) | 1.18 (0.35) | | |
| SRI 193 nm per initial excitation against Ar | 6 | 0.66 (0.47) | | | | | |
| SRI 290 nm (Ar2F*) per initial excitation against Ar | 9 | 0.90 (0.54) | | | | | |
| SRI 193 nm per initial excitation against F2 | 20 | 0.84 (0.24) | | | | | |
| Morozov pure-Ar VUV efficiency | 16 | 2.15 (0.77) | | 1.75 (0.56) | 1.89 (0.64) | | |
| the constants' own sources: decays, Stern-Volmer plots | 69, 15 | 0.97 (0.12), 1.06 (0.17) | | | | | |

- **The drive improves the model overall.** Over the 11 held-out tests the
  mean |ln gm| goes from 0.314 (LIGHT_ArF's own input) to 0.256 under the
  drive, and the rms of the logs over their 87 points from 0.476 to 0.392.
  It removes the own input's failure at low pumps: under the beam Te sits
  at 1.5 to 2.2 eV at every device and shape point (21 to 94 eV at
  Rokni's pumps on the own input), and Rokni's current ratio goes from
  x1.67 to x0.97. Electra rises from x0.64 to x0.84.
- **It overshoots at high pump.** Sandia goes to x1.82, Keio at its printed
  pump to x1.14, and Jancaitis's ArF* peak rises x5.5 from 74 to 562
  A/cm2 against x9.1 measured. Under the beam the model makes ArF(B) at
  0.25 to 0.30 of the input at every device pump, where the measured share
  falls from 0.36 (Electra, 0.90 MW/cc) to 0.16 (Sandia, 5 MW/cm3).
- **It does not reach Morozov.** Pure argon improves from x2.15 to x1.75
  (x1.70 at 500-1400 hPa, where his 0.33 +- 0.04 is printed). That is the
  network's one-photon ceiling, every ion and excitation one 127-nm photon,
  0.538 with the drive's yields, near the 0.53 theoretical maximum Morozov
  quotes (Lorents 1976). His argon sits at 62 % of that maximum while his
  Kr and Xe sit near theirs, and run as a continuous beam the model is
  flat at 0.54 from 151 to 1000 hPa. No input closes this: the network
  lacks a measured argon VUV loss, or his argon value is low.
- **1Dv0_2 is v0_2** on either drive. On the own input its centre cell
  gives v0_2's numbers to 4e-4 (2e-5 at the reference point), except at
  Rokni's 2 and 4 atm, where v0_2's own-input Te passes Fluid1D's 40-eV
  ceiling (x0.90 and x0.94). Under the beam it gives v0_2's beam numbers to
  1.1e-4 at every device point, to 5e-5 at Rokni's five pressures, where
  v0_2's Te stays at 1.5 to 2.2 eV, and to 5e-6 at the reference point.
  The gap average lies 2.1 to 2.4 % lower (3.1 % on the own input), the
  wall cells, which shrink with the mesh.
- **What the data say is missing** is LIGHT_ArF's: an Ar2F* loss that grows
  with excitation (Jancaitis's Ar2F* x1.34 under the beam at the top
  currents, SRI's 290-nm yield x0.35 at 1 atm and x1.4 at 9 atm), an argon
  VUV loss, and Electra's gain falling to x0.70 at 1.36 atm where the
  model's falls to x0.92.

## Files

```
LIGHT_ArF_ebeam/
  README.md, CHANGELOG.md      this file; the changelog (Known-bad first)
  run_v0_2.py                  the 0-D ArF model
  run_1Dv0_2.py                the 1-D ArF model (Python >= 3.11 with numba)
  v0_2/                        the 0-D ArF model (README.md, tests/),
                               beam.py the electron-beam drive
  1Dv0_2/                      the 1-D ArF model (README.md, tests/)
  docs/accuracy_overview.pdf   the accuracy overview (accuracy_overview.tex,
                               written by make_overview.py with its
                               overview_*.tex files and figures/)
  docs/energy_ladder.tex       the energy ladder, by make_energy_ladder.py
  docs/system_map.pdf          the system map page (system_map.tex, with
                               system_map_numbers.tex by make_system_map.py)
  docs/arf_data/               the ArF evidence: SOURCES_CHECKED.md,
                               RATES_RESEARCH.md, SCORING_DATA.md and
                               fig_readings/, ar_tail.py (the tail table),
                               beam_drive.py (LIGHT_ArF's prototype, for
                               ar_tail.py), score_arf.py, SCORE.md and
                               score_0d_LIGHT_ArF.json
```
