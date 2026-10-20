# LIGHT_ArF_ebeam 1Dv0_2: the 1-D ArF model, on both drives

1Dv0_2 runs LIGHT_ArF_ebeam v0_2's ArF chemistry in every cell of a gap, on
Nathaniel's Fluid1D (`Kinetic_Modeling/nathaniel/1DModelKrF`, imported and
never edited). It is LIGHT's 1Dv0_2 (the 1-D layer on KrF v0_2) converted
to Ar/F2 on 2026-09-30, with the electron-beam drive added the same day
(LIGHT_ArF_ebeam). The KrF original stays in `LIGHT/1Dv0_2/`, the ArF one
without the beam in `LIGHT_ArF/1Dv0_2/`.

The chemistry files (`model.py`, `network.py`, `species.py`, `xsec.py`,
`report.py`, `beam.py`, `data/`) are LIGHT_ArF_ebeam v0_2's, byte for
byte, and `model.py` still runs them in 0-D. The network, its sources, its
closure and the drive are in `../v0_2/README.md`.

```bash
PY=~/.venvs/f2sweep/bin/python   # Python >= 3.11 with numba, as Fluid1D needs
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py                       # the reference point, 1-cm gap, 128 cells
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py --compare             # and v0_2's 0-D run of the same inputs
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py --drive beam --compare # the same under the electron beam
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py --rate arf_yield=1.0  # one band set (repeatable)
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py --trimer-tg-pow -3.9  # NRL's trimer law
$PY LIGHT_ArF_ebeam/run_1Dv0_2.py --profile power.txt   # relative power per cell
$PY -m pytest LIGHT_ArF_ebeam/1Dv0_2/tests -q           # the 1-D tests and the 0-D copies
LIGHT_1D_SLOW=1 $PY -m pytest LIGHT_ArF_ebeam/1Dv0_2/tests/test_1d.py -q   # + the reference point
LIGHT_DOCS_SLOW=1 $PY -m pytest LIGHT_ArF_ebeam/1Dv0_2/tests/test_docs.py -q # the overview's runs, ~20 min
```

## What it is

LIGHT's 1-D layer (`LIGHT/1Dv0/README.md`: the input, the walls, how a
step works, what comes from Fluid1D) with the ArF chemistry in the stage.
The input is a given power, uniform across the gap unless `--profile` sets
a shape, deposited one of v0_2's two ways (`--drive`). LIGHT's own input
(`power`) puts it into the electron gas. The electron beam (`beam`,
`beam.py`) makes ion pairs at W eV each: the stage adds their electrons,
Ar+ ions, excited atoms and prompt F2 dissociations to every cell in every
substep, the secondaries' share to the electron energy beside the field's
work, and the cascade's prompt heat to the cell's gas (ledger
`heat_beam`). With no beam the stage's arithmetic is the own input's, bit
for bit (the same run on LIGHT_ArF's 1Dv0_2 and on this one agrees in every
recorded number). The gas fills a gap between two grounded
electrodes. Each cell runs every row of the ArF network (69 rows: the 38 Ar
electron-impact rows on the two-temperature tail, the F2 rows, the harpoon,
ArF(B), the trimer, the ions, the Ar2* route), its electron energy and its
gas temperature, while Fluid1D moves the charged species and the electron
energy between the cells and into the walls.

What changed from the KrF layer:

- the rare gas is argon everywhere the layer names it (the gas mixture,
  the seed ion, the seed taken out of the argon, the elastic loss);
- Fluid1D's cited ion mobilities for Ar+, Ar2+ and F- in Ar (Ellis et al.:
  1.52, 1.83, 3.33 cm2/V s);
- the ion-ion rows are Hoffman & Moreno's form ("hm"), evaluated per cell
  with the Debye screening and the Thomson Tg law;
- the band rows carry the band modifiers (a share, a rest, a k tau product
  over tau_B), resolved by the 0-D model to constants of the run;
- the ledger counts ArF formed and the 193-, 284- and 127-nm bands, and
  the summary prints the ArF warnings.

## What it gives

v0_2's reference point: 1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to
2 us, seed n_e 1e15, a 1-cm gap, dt = 50 ps (2026-09-30), on LIGHT's own
input:

| | v0_2, 0-D | centre cell | gap, 64 cells | gap, 128 cells |
|---|---|---|---|---|
| 193 nm | 5.0522 % | 5.0523 % | 4.896 % | 4.974 % |
| 127 nm | 4.6040 % | 4.6039 % | 4.459 % | 4.530 % |
| Tg at the end | 675.19 K | 675.19 K | 664 K | 669 K |
| cut at the walls | | | 3.11 % | 1.55 % |
| phantom negative work kept out | | | 12.7 % | 2.4 % |
| to the walls | | | 0.040 % | 0.036 % |
| wall time | | | 136 s | 164 s |

(percent of the input; the centre cell is the same on both meshes)

The same under the electron beam (`--drive beam`):

| | v0_2, 0-D | centre cell | gap, 64 cells | gap, 128 cells |
|---|---|---|---|---|
| 193 nm | 6.3398 % | 6.3398 % | 6.196 % | 6.264 % |
| 127 nm | 2.0907 % | 2.0907 % | 2.027 % | 2.060 % |
| Tg at the end | 731.25 K | 731.25 K | 722 K | 726 K |
| cut at the walls | | | 0 | 0 |
| phantom negative work kept out | | | 1522 % | 225 % |
| to the walls | | | 2.28 % | 1.18 % |
| wall time | | | 125 s | 169 s |

The interior is v0_2's 0-D run, to 2e-5 in the 193- and 127-nm shares and
in Tg on the own input and to 5e-6 under the beam.

The wall cells change their channel with the drive, not their size. On
the own input the walls drain their electrons and the mean-energy ceiling
cuts what the power leaves them. Under the beam they keep receiving pairs,
nothing is cut, and the walls take the pairs made in them instead. Both
halve with each mesh doubling, as does the gap average's shortfall from
the centre. Fluid1D's phantom negative work is far larger under the beam,
15 times the input on 64 cells, and falls x0.15 to x0.17 per doubling (89,
15 and 2.2 times on 32, 64 and 128 cells). It never enters the books,
which close to 3e-11 % of the input on either drive.

With a uniform input the interior is v0_2's 0-D run, and the gap average
carries the wall cells, unresolved sheaths that shrink with the mesh
(`LIGHT/1Dv0/README.md`, "The wall cells are not physics").

In the scoring (`../docs/arf_data/SCORE.md`, 64 cells) the centre cell
gives v0_2's Electra, Keio and Sandia numbers to 4e-4 on the own input and
to 1.1e-4 under the beam. On the own input Rokni's curve follows to 2e-4
up to 1 atm, and at 2 and 4 atm, where v0_2's Te passes Fluid1D's 40-eV
ceiling, the centre falls to x0.90 and x0.94 of v0_2. Under the beam
v0_2's Te stays at 1.5 to 2.2 eV and the centre follows v0_2 at all five
pressures to 5e-5. The gap average is 3.1 % lower on the own input and 2.1
to 2.4 % lower under the beam. On ArF data the 1-D model is as accurate as
the 0-D one and no more: nothing in these experiments varies across the
gap that a uniform input could follow.

## Limits

LIGHT's 1-D layer's: the wall cells, the gas temperature Fluid1D's
transport reads (held at tg0), Maxwellian electron transport, first order
in dt, the mean energy held to 0.03-60 eV, and no discharge. The seed rule
of KrF (1e13 cm^-3 or more) does not bind here: v0_2's 0-D Te stays below
28 eV from a 1e10 seed at 300 kW/cm3, inside the 40-eV ceiling. A band is
for sensitivity, never to match a device. The beam deposits uniformly
unless `--profile` shapes it, and every cell takes the fill's yields per
pair.

## Files

```
1Dv0_2/
  fluid1d.py     the bridge (the ArF names, the "hm" kind per cell, the
                 beam's products per cell in the stage)
  sim1d.py       Params1D (with v0_2's rates and drive), LightSimulation,
                 run(), Result1D (bands line, ArF warnings, share_193)
  run.py         the command line (LIGHT_ArF_ebeam/run_1Dv0_2.py), with
                 --rate and --drive
  model.py ...   LIGHT_ArF_ebeam v0_2's chemistry, its drive (beam.py) and
                 its 0-D model, byte for byte
  tests/test_1d.py   the 1-D tests (skip without Python >= 3.11 + numba)
  tests/test_docs.py the accuracy overview's runs file against fresh runs
  tests/test_*.py    v0_2's 0-D tests, pointed at 1Dv0_2's own modules
  runs/          run output (git-ignored)
```

## Tests

`test_1d.py`: the bridge tests (species, the 50 cross-section processes
round-tripped through Fluid1D, every rate against `Model._k`, the 13 band
rows, the trimer law, the tail, the elastic term, the stage's energy, atom
and charge conservation) with the bands at their defaults and with five
set, then the runs: the centre cell against v0_2's 0-D run, two bands
reaching the 1-D chemistry, the books closing, the wall cells, Te first
order in dt, the result file, the command line, and the Fluid1D hashes.
A mutation check planted 12 bugs in the new code (the band modifiers, the
"hm" kind in the model and the bridge, the k tau products, the ArF(B)
energy, the tail gas, the elastic column, the seed, the trimer law in the
cells) and the tests caught all 12. The beam's tests: the stage's energy,
atom and charge conservation with the beam's sources (and the same stage
with a zero source equal to the own input's to the last bit), the centre
cell against v0_2's 0-D beam run (its emission by band, ArF formed, the
beam's prompt heat, species, Tg and Te), the books closing under the beam,
the beam's electrons entering a cell at the ceiling inside the energy
bounds, and the command line's `--drive beam`. A second mutation check
planted 8 bugs in the beam's path (no electron per pair, the kernel
ignoring the species source or adding it with the wrong sign, the
secondaries' energy dropped, the prompt heat left out of Tg or of the
books, pairs at the measured W in place of the drive's, the bounds taken
on the electrons before the source) and the tests caught all 8. The 0-D
tests are v0_2's, pointed at
1Dv0_2's copies (the two that drive v0_2's 0-D command line are dropped),
`test_beam.py` among them.
`test_docs.py` (with `LIGHT_DOCS_SLOW=1`) regenerates the runs of
`docs/accuracy_overview.pdf`, v0_2's and 1Dv0_2's, and compares them with
`docs/overview_runs.tex` byte for byte.
