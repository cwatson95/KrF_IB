# LIGHT 1Dv0: the 1-D model

1Dv0 runs LIGHT v0_1's chemistry in every cell of a gap, on Nathaniel's
Fluid1D (`Kinetic_Modeling/nathaniel/1DModelKrF`). v0_1 is the 0-D model
and stays as it is. The owner made 1Dv0 as a copy of v0_1 on 2026-09-29
and asked for it to become the 1-D model on Fluid1D. Its chemistry files
(`model.py`, `network.py`, `species.py`, `xsec.py`, `two_level.py`,
`report.py`, `data/`) are v0_1's, unchanged, and `model.py` still runs
them in 0-D. What v0_1 changed and why is in `LIGHT/v0_1/README.md`, its
literature in `SOURCES.md`.

```bash
PY=~/.venvs/f2sweep/bin/python   # Python >= 3.11 with numba, as Fluid1D needs
$PY LIGHT/run_1Dv0.py                       # v0_1's reference point, 1-cm gap, 128 cells, 80 s
$PY LIGHT/run_1Dv0.py --compare             # and v0_1's 0-D run of the same inputs
$PY LIGHT/run_1Dv0.py --cells 256 --dt-ps 25
$PY LIGHT/run_1Dv0.py --pressure 2 --power 500 --f2-pct 0.3
$PY LIGHT/run_1Dv0.py --profile power.txt   # relative power per cell
$PY -m pytest LIGHT/1Dv0/tests -q           # the 1-D tests and the 0-D copies
LIGHT_1D_SLOW=1 $PY -m pytest LIGHT/1Dv0/tests/test_1d.py -q   # + the reference point
```

Fluid1D is read from `nathaniel/1DModelKrF/src` beside LIGHT, or from
`FLUID1D_SRC`. The f2sweep venv also holds an editable install of another
checkout (`~/coding/nathaniel/1DModelKrF`), and `fluid1d.load()` refuses
a Fluid1D imported from anywhere but the chosen path. A run writes
`1Dv0/runs/<name>/result.json` and the files Fluid1D read
(`fluid1d_inputs/`: the cross sections and the reaction file generated
from this package).

## What it is

The input is v0_1's: a given power into the electron gas for the pulse,
deposited directly and uniformly across the gap. `--profile` gives it a
shape across the gap (one relative number per cell) and keeps its mean.
There is no drive, beam or cavity model, which is LIGHT's rule.

The gas fills a gap between two grounded electrodes. Fluid1D moves the
charged species (e, Kr+, Kr2+, F-) by drift and diffusion in the field
their own space charge makes, carries the electron energy with them, and
lets the walls absorb what reaches them. The neutrals stay where they
are. In 1 us at 1 atm they diffuse a few microns.

Every cell runs v0_1's whole basic network with v0_1's rates: the
two-temperature tail read at the cell's own n_e/N and F2/N, the captured
electrons' exact energies, the Tg laws, the screened ion-ion
recombination, the elastic loss n_e N_Kr kel(Te)(1 - Tg/Te), and the gas
temperature at constant volume.

## How a step works

One step of dt is Fluid1D's source-coupled split:

1. chemistry and electron energy over dt/2 (the stage),
2. transport over dt (Fluid1D: field, drift-diffusion, electron energy
   fluxes, walls),
3. chemistry and electron energy over dt/2.

The stage (`fluid1d.CellChemistry.stage`) runs every row of the network in
every cell with Fluid1D's linearly implicit trial step, which moves the
species and debits each event's electron energy together. It sits in an
outer substep loop that follows Fluid1D's own, with two changes. The
rates are v0_1's, refreshed every substep (the tail rows also read the
cell's tail ratio, which Fluid1D's loop cannot). And they are read at the
substep's starting energy, which makes a steady Te exact. Read after the
power is added, as Fluid1D's loop does, the centre's Te sat 3 to 4 % low
at dt = 50 ps.

The rates come from v0_1's own tables. The tables' Te rows sit on v0_1's
Te nodes, so an LXCat row is v0_1's table itself. The tail rows sit on
v0_1's (Te, T_tail) nodes with -E_break/Te added back exactly. The Tg,
density and ion rows are v0_1's formulas per cell. Against `Model._k` the
rates agree to rounding, and the tail to 8e-6.

Counter columns in the stage record each cell's reaction heat, emission
by band, energy returned to the electrons (Penning), captured electrons
and KrF formed, as extents times v0_1's per-event energies
(`network.energy_routing`). Tg rises by the heat at constant volume, as
in v0_1. The stage's energy source is the given power over its own half
step, plus the field's positive work on the electrons and minus the
elastic loss, which Fluid1D's split delivers.

The books close. The gap's chemical, electron and thermal energy,
recomputed from its state, equals what the ledger booked (input minus
radiation, plus the bounds and the transport) to 1e-12 of the input.

## What it gives at v0_1's reference point

1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us, seed n_e 1e15, a
1-cm gap. v0_1's 0-D run gives 4.73236 % at 248 nm and ends at 762.49 K.

| cells | dt | centre 248 nm | gap 248 nm | cut at the walls | field work given / not | to the walls | wall time |
|---|---|---|---|---|---|---|---|
| 32 | 50 ps | 4.73229 % | 4.438 % | 6.21 % | 0 / 87.7 % | 0.059 % | 47 s |
| 64 | 50 ps | 4.73227 % | 4.585 % | 3.09 % | 0 / 17.1 % | 0.053 % | 56 s |
| 128 | 50 ps | 4.73227 % | 4.659 % | 1.81 % | 0.31 / 5.1 % | 0.077 % | 80 s |
| 128 | 25 ps | 4.73234 % | 4.659 % | 1.50 % | 0 / 4.7 % | 0.073 % | 144 s |
| 256 | 50 ps | 4.73226 % | 4.696 % | 0.95 % | 0.24 / 1.6 % | 0.057 % | 127 s |

(percent of the input; four numba threads; the centre cell ends at
762.48 to 762.49 K at every mesh)

The interior is v0_1. Every cell away from the walls is the same cell,
and it is v0_1's 0-D run to 2e-5 in the 248-nm share and 0.01 K in Tg.
With uniform power the walls change nothing inside the gap. The ions
spread about ten microns in a microsecond at 1 atm (ambipolar), and under
0.1 % of the energy reaches the walls. This is the comparison's finding
for power (e-beam-like) pumping, now measured in 1-D. The 1-D model earns
its keep with a non-uniform input (`--profile`).

At Tahoma, the one device the LIGHT versions are scored on (665 Torr,
2.8 Torr F2, 330 kW/cm^3 for 450 ns, ML's constants imposed, g0 from
KrF(B) over 200 to 400 ns), the centre cell on 64 cells gives v0_1's g0 to
2e-4 with the default network, 5.93 to 5.96 %/cm or x0.93 of the measured
6.39, and to 7e-4 with the Kr2* group, 6.63 %/cm or x1.04, from seeds of
1e13 to 1e15 cm^-3. So 1Dv0 scores as v0_1 does, and the review
(`docs/REVIEW/v0_vs_v0_1_vs_v0_2.pdf`) is its comparison with v0_2. The
gap average, x0.90, is the wall cells' (below).

The wall cells are not physics. The sheath at a wall is about a micron
thick, far below any mesh here, so it lives inside the cell next to each
wall. There the walls drain the electrons, and v0_1's power per volume
pins what is left at the mean-energy ceiling (60 eV), which cuts the rest.
Fluid1D's split gives the electrons only the positive part of the field's
work there, and keeps out the large negative part that its own notes call
phantom (a sub-Debye cell cannot resolve the sign of the electron flux
times the field). The cut halves with each mesh doubling, the negative
field work falls faster (x0.2 to x0.3), and the gap average's shortfall
from the interior halves (0.29, 0.15, 0.07, 0.04 points). Read a gap
average with its mesh, and read the interior for v0_1's physics.

## What comes from Fluid1D, and what does not

From Fluid1D, imported and never edited: the mesh, the Scharfetter-Gummel
drift-diffusion of the charged species and of the electron energy, the
semi-implicit field, the wall fluxes, the source-coupled split, the
linearly implicit trial step, its ion mobilities (Kr+ 0.84, Kr2+ 1.22,
F- 2.29 cm^2/V s in Kr, Ellis et al., `data/ion_transport/ion_mobility.toml`),
and the electron mobility and diffusion, from its Maxwellian solver
(Hagelaar and Pitchford 2005 moments) run on 1Dv0's own cross sections.

Not from Fluid1D: its reaction network (`krf.toml`), its species list,
its rate tables, its two-term Boltzmann solver, its circuit and its
cavity. `fluid1d.write_inputs` converts 1Dv0's two md5-pinned LXCat files
(Biagi-v7.1 Kr, SIGLO F2) to Fluid1D's JSON and CSV (exact to rounding)
and writes a reaction file of the 12 electron-impact rows. Fluid1D builds its
electron transport table from those, and runs no chemistry from them.

`sim1d.LightSimulation` subclasses Fluid1D's `Simulation`. It swaps in the
species registry, the transport table, a plasma model with grounded walls
and the stage, and leaves the rest of the step to Fluid1D. A test holds
the 22 Fluid1D functions it overrides or relies on to recorded hashes, so
a change upstream fails loudly instead of silently.

## Limits

- The wall cells, above. They converge away with the mesh, first order.
- Transport reads the gas at tg0 (300 K): the ions' diffusion and their
  speed at the wall. The chemistry reads each cell's own Tg.
- The electron mobility and diffusion are Maxwellian at the mean energy.
  Below the tail's break (9.9 eV) v0_1's two-temperature EEDF is the same
  Maxwellian, and transport is set there.
- First order in dt. In a 20-ns pulse at 50 ps the centre's Te stays
  within 1.5 % of v0_1's (the largest in the first 0.1 ns of the ramp),
  and 25 ps halves it. At the reference point the 248-nm share moves by
  1.5e-5 between 50 and 25 ps.
- The mean energy is held to 0.03-60 eV (Te 0.02-40 eV), the tables'
  range. At the reference point only the wall cells reach the ceiling,
  and no cell the floor.
- The seed. Give 1Dv0 a seed of 1e13 cm^-3 or more at these powers. v0_1
  breaks down from any seed because its 0-D Te has no ceiling: from ML's
  1e10 cm^-3 at Tahoma's 330 kW/cm^3 it runs to 518 eV in 0.1 ns before
  the two-temperature tail (T_t = r Te, r about 0.004 at that ionization
  degree) reaches the Kr thresholds. 1Dv0 holds Te below 40 eV, so from
  that seed it never breaks down: the ceiling cuts the whole input
  ("cut by the mean-energy bounds" reads 100 %) and no KrF forms. From
  1e13 up it breaks down as v0_1 does. v0_1's scores barely depend on the
  seed (Tahoma's g0 is 5.94 %/cm from 1e10 to 1e14 and 5.96 at 1e15).
- No discharge. Fluid1D's own drive is a circuit, which LIGHT's rule
  leaves out. Running one would also need a field-driven EEDF in place of
  v0_1's power-driven tail, and the comparison found Fluid1D's discharge
  not mesh-converged (`comparisons/README.md`).

## Files

```
1Dv0/
  fluid1d.py     the bridge: Fluid1D import, species, cross sections and
                 reaction file, the network as per-cell arrays, the stage
  sim1d.py       Params1D, LightSimulation (Fluid1D's Simulation with the
                 stage), run(), Result1D (records, profiles, ledger, fate)
  run.py         the command line (LIGHT/run_1Dv0.py)
  model.py ...   v0_1's chemistry and its 0-D model, unchanged
  tests/test_1d.py   the 1-D tests (skip without Python >= 3.11 + numba)
  tests/test_*.py    v0_1's 0-D tests, pointed at 1Dv0's own modules
  runs/          run output (git-ignored)
```

## Tests

`test_1d.py` has 13 tests, about 10 s, and a slow one, about 50 s. The
bridge: the species against `species.py`, the cross sections through
Fluid1D's own loader, every rate against `Model._k`, the vectorized tail
against `xsec.py`, the elastic term against `Model.rhs`. The stage: from
random states it conserves energy per cell, Kr and F atoms and charge.
The runs: the centre cell against v0_1's 0-D run (its starting state to
1e-12; emission by band, KrF formed and the 248-nm share to 1e-3; species
to 1e-2; Tg to 0.5 K; Te to 1 % after 2 ns), the books closing, a uniform
interior with the cut in the wall cells only, Te converging first order
in dt, the result file and the command line. And the Fluid1D hashes. The
slow test runs the reference point on 32 cells against v0_1 (248 nm to
1e-4). A mutation check planted 16 bugs in `fluid1d.py` and `sim1d.py`
(unit factors, the rate order, the energy flows, the Tg update, the
counters, the seed), and the fast tests caught all 16.

`tests/conftest.py` puts this LIGHT/ first on sys.path. pytest stops its
package search at a directory name that is not a Python identifier, so
for '1Dv0' it never does so itself. Without it the tests import no
package, or the one in the current directory.

The 0-D tests are v0_1's, pointed at 1Dv0's copies, and run under
LIGHT's usual `python3` as v0_1's do. v0_1's `test_docs.py`, which checks
v0_1's documents, was not kept.
