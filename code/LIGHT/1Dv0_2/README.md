# LIGHT 1Dv0_2: the 1-D model on v0_2's chemistry

1Dv0_2 runs LIGHT v0_2's chemistry in every cell of a gap, on Nathaniel's
Fluid1D (`Kinetic_Modeling/nathaniel/1DModelKrF`). It is 1Dv0's 1-D layer
on v0_2, the version the review calls the best of the three
(`docs/REVIEW/v0_vs_v0_1_vs_v0_2.pdf`). The owner asked for it on
2026-09-29 ("Then make a 1-D model with v0_2 chemistry"). 1Dv0, the same
layer on v0_1's chemistry, and v0_2 itself stay as they are.

The chemistry files (`model.py`, `network.py`, `species.py`, `xsec.py`,
`two_level.py`, `report.py`, `data/`) are v0_2's, byte for byte, as they
stood when 1Dv0_2 was made, and `model.py` still runs them in 0-D. What
v0_2 changed and why is in `LIGHT/v0_2/README.md`, its literature in
`SOURCES.md`. `docs/accuracy_overview_1Dv0_2.pdf` says how accurate
1Dv0_2 is: its chemistry scores as v0_2's, and the 1-D layer is scored
against v0_2's 0-D run (2026-09-30).

```bash
PY=~/.venvs/f2sweep/bin/python   # Python >= 3.11 with numba, as Fluid1D needs
$PY LIGHT/run_1Dv0_2.py                      # v0_2's reference point, 1-cm gap, 128 cells, 80 s
$PY LIGHT/run_1Dv0_2.py --compare            # and v0_2's 0-D run of the same inputs
$PY LIGHT/run_1Dv0_2.py --rate kr2_f=3e-10   # one band set (repeatable)
$PY LIGHT/run_1Dv0_2.py --trimer-tg-pow -0.5 # the trimer's weak law
$PY LIGHT/run_1Dv0_2.py --profile power.txt  # relative power per cell
$PY -m pytest LIGHT/1Dv0_2/tests -q          # the 1-D tests and the 0-D copies
LIGHT_1D_SLOW=1 $PY -m pytest LIGHT/1Dv0_2/tests/test_1d.py -q   # + the reference point
$PY LIGHT/docs/make_overview_1Dv0_2.py       # the accuracy overview's numbers, ~10 min
open LIGHT/docs/accuracy_overview_1Dv0_2.pdf # how accurate 1Dv0_2 is
```

## What it is

The same model as 1Dv0 (`1Dv0/README.md`: the input, the walls, how a step
works, what comes from Fluid1D), with v0_2's chemistry in the stage:

- The basic network's 44 rows: v0_1's 37 with the Kr2* route (Kr* + 2Kr,
  Kr2* radiation at 147 nm, Kr2* + F2 in two channels, Kr2* + F,
  Kr2* + Kr2*) and Kr2F* + F.
- v0_2's bands. The eight rows v0_2 carries as a band (rate kind
  `("param", name[, frac])`) run at the value v0_2's model resolves:
  `Params1D.rates` (`--rate NAME=VALUE`) or the band's default, which is a
  measurement or 0. They are constants of the run, the same in every cell,
  times their Tg law. The trimer's law is `trimer_tg_pow`, as in v0_2.
- The 147-nm band has its own emission counter, and the summary prints
  v0_2's bands line (the bands set, and the unmeasured ones left at 0) and
  its gas-temperature warnings (Kr* + 2Kr's law above 451 K, the trimer's
  above 800 K, redissociation above 1000 K), read in the centre cell.

v0_2's `--bands` (each band's two ends, one run per end) stays with
`run_v0_2.py`. With a uniform input the gap's interior is the 0-D run, so
a band sweep in 1-D would repeat v0_2's own.

## What it gives

v0_2's reference point: 1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us,
seed n_e 1e15, a 1-cm gap, dt = 50 ps.

| | v0_2, 0-D | centre cell | gap, 64 cells | gap, 128 cells |
|---|---|---|---|---|
| 248 nm | 4.66672 % | 4.66665 % | 4.521 % | 4.594 % |
| 147 nm | 20.0184 % | 20.0185 % | 19.391 % | 19.703 % |
| Tg at the end | 649.35 K | 649.35 K | 638 K | 644 K |
| cut at the walls | | | 3.10 % | 1.54 % |
| to the walls | | | 0.042 % | 0.041 % |
| wall time | | | 59 s | 79 s |

(percent of the input; the centre cell is the same on both meshes to 1e-6;
four numba threads)

The interior is v0_2's 0-D run, to 1.5e-5 in the 248-nm share and 5e-6 in
the 147-nm share. With Kr2* + F at the top of its band (`--rate
kr2_f=1e-9`) it stays so: 9.49826 % at 248 nm against v0_2's 9.49835 %,
and 726.35 K at the end in both. The gap averages carry the wall cells,
which are unresolved sheaths and shrink with the mesh, as in 1Dv0
(`1Dv0/README.md`, "The wall cells are not physics").

At Tahoma on v0_2's own power input (665 Torr, 2.8 Torr F2, 330 kW/cm^3
for 450 ns, ML's constants imposed, g0 from KrF(B) over 200 to 400 ns, 64
cells), the centre cell's g0 is v0_2's to 3e-4: 5.528 %/cm against 5.530
from a seed of 1e13 cm^-3, and 5.554 against 5.555 from 1e15, x0.865 and
x0.869 of the measured 6.39. The review's v0_2 column is 1Dv0_2's.

## Limits

1Dv0's, the same model with the same layer: the wall cells, the gas
temperature transport reads (tg0), Maxwellian electron transport, first
order in dt, the mean energy held to 0.03-60 eV, and no discharge. See
`1Dv0/README.md`. Two apply to v0_2 in numbers:

- The seed. Give 1Dv0_2 a seed of 1e13 cm^-3 or more at these powers.
  From ML's Tahoma seed, 1e10, v0_2's 0-D Te runs to 73 eV before its tail
  reaches the Kr thresholds (66 eV from 1e12, 38 eV from 1e13). 1Dv0_2
  holds Te below 40 eV, so from 1e10 or 1e12 it never breaks down: the
  ceiling cuts 99.9 % of the input or more, no KrF forms, and the summary
  says so. v0_2's own g0 barely moves with the seed (5.530 %/cm from 1e10,
  5.555 from 1e15).
- A band is for sensitivity, never to match a device (v0_2's rule). The
  1-D model carries it through unchanged.

Two more, found for the accuracy overview (2026-09-30):

- Above 30 eV, the end of v0_2's Te tables, the stage holds Kr2+ + e
  (R11p, a power law) at its 30-eV value, where v0_2's `Model._k` follows
  the law. At the 40-eV ceiling the stage's rate is 17 % high. The wall
  cells sit there, and the interior passes 30 eV only while a low seed
  breaks down, and at Tahoma from 1e13 the centre cell still gives v0_2's
  g0 to 3e-4.
- The cells beside the wall cell do not converge. On 32 to 256 cells the
  first of them sends 4.668 to 4.751 % out at 248 nm against the
  interior's 4.667 %, and on 256 cells they alternate about the interior
  out to about 0.02 cm, with n_e and Te. The cells are 39 um wide where
  the Debye length is about 1 um, so no mesh here resolves the sheath. The
  gap average still converges first order.

## Files

```
1Dv0_2/
  fluid1d.py     the bridge (1Dv0's, plus v0_2's 'param' rows as constants
                 of the run)
  sim1d.py       Params1D (with v0_2's rates), LightSimulation, run(),
                 Result1D (with v0_2's bands line and Tg warnings)
  run.py         the command line (LIGHT/run_1Dv0_2.py), with --rate
  model.py ...   v0_2's chemistry and its 0-D model, unchanged
  tests/test_1d.py   the 1-D tests (skip without Python >= 3.11 + numba)
  tests/test_*.py    v0_2's 0-D tests, pointed at 1Dv0_2's own modules
  runs/          run output (git-ignored)
```

## Tests

`test_1d.py` runs 22 tests in about 12 s, and a slow one. The bridge tests
(species, cross sections, every rate against `Model._k`, the tail, the
elastic term, the stage's energy, atom and charge conservation) run twice:
with the bands at their defaults and with the three unmeasured ones set,
so every band row runs at a nonzero rate. One test holds each of the
eight band rows to frac x (rate or default) x its Tg law. The runs: the
centre cell against v0_2's 0-D run (its starting state, emission by band,
KrF formed, species, Tg, Te), the top of the Kr2* + F band reaching the
1-D chemistry (the centre cell follows v0_2's run with the same band,
which differs from the default by far more), the books closing, the wall
cells, Te converging first order in dt, the result file, the command line
with `--rate`, and the Fluid1D hashes. The slow test runs the reference
point on 32 cells against v0_2. A mutation check planted 19 bugs (1Dv0's
16 and three in the bands: the band rows' kind, the rates reaching the
cell, the bands line), and the fast tests caught all 19.

The 0-D tests are v0_2's, pointed at 1Dv0_2's copies, and run under
LIGHT's usual `python3`. v0_2's `test_docs.py` and its two tests of
`run_v0_2.py`'s command line (`--rate`, `--bands`) were not kept: 1Dv0_2's
`run.py` is the 1-D command line. The summary and bands-line checks of the
first stay, on the 0-D report. `tests/conftest.py` puts this LIGHT/ on
sys.path, as in 1Dv0.

`test_docs.py` (2026-09-30) holds the accuracy overview to the model. Under
any Python it checks that the chemistry files still have v0_2's md5 sums
and that the overview prints only numbers its generated files define,
v0_2's energy ladder among them.
With Fluid1D it compares the overview's 1-D numbers with a fresh run, the
slow ones (about 10 minutes) with `LIGHT_DOCS_SLOW=1`. The rate test above
draws its random states below 30 eV, so it cannot see the Kr2+ + e
difference under Limits.
