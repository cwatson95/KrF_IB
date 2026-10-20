# LIGHT 1Dv0_3: the 1-D model on a mesh graded to the walls

1Dv0_3 is 1Dv0_2 with OSIRIS's grid scaling. The owner asked on 2026-10-01:
"In Kinetic_Modeling/osiris/ the software uses grid scaling for where
important interactions occur. Can we apply that to the new copy of 1Dv0_2
titled 1Dv0_3". 1Dv0_3 is the owner's copy of 1Dv0_2, made that day. Its
chemistry is v0_2's, byte for byte (`model.py`, `network.py`, `species.py`,
`xsec.py`, `two_level.py`, `report.py`, `data/`), and 1Dv0_2 stays as it
is. What changes is where the cells are, since 2026-10-01 (c) that the
neutral gas moves, and since 2026-10-06 that the electric field's work is
kept to the field's energy. The resolved walls showed that the gas must
move there, and `neutrals.py` lets it flow, conduct heat into the
electrodes and its F2 and F diffuse. The field's work now moves energy
from the electrons to the ions through the field, and none of it is input
("The field's work", below). On its uniform mesh with the gas still and
Fluid1D's field work (`mesh="uniform"`, `neutrals="still"`,
`field_work="fluid1d"`) 1Dv0_3 runs 1Dv0_2's model to the bit.

```bash
PY=~/.venvs/f2sweep/bin/python   # Python >= 3.11 with numba, as Fluid1D needs
$PY LIGHT/run_1Dv0_3.py                      # reference point, 128 cells from 0.95 um, the gas moving, 13 min
$PY LIGHT/run_1Dv0_3.py --compare            # and v0_2's 0-D run of the same inputs
$PY LIGHT/run_1Dv0_3.py --cells 64 --wall-cell 2                   # the quick mesh, 6 min
$PY LIGHT/run_1Dv0_3.py --neutrals still                           # the gas held still, as in 1Dv0_2
$PY LIGHT/run_1Dv0_3.py --field-work fluid1d                       # Fluid1D's field work, as before 2026-10-06
$PY LIGHT/run_1Dv0_3.py --mesh uniform --cells 128 --neutrals still --field-work fluid1d  # 1Dv0_2 to the bit
$PY -m pytest LIGHT/1Dv0_3/tests -q          # the 1-D, gas, transport and mesh tests and the 0-D copies
LIGHT_1D_SLOW=1 $PY -m pytest LIGHT/1Dv0_3/tests/test_1d.py -q   # + the reference point, twice
```

## What OSIRIS does

OSIRIS (`Kinetic_Modeling/osiris/`, the particle-in-cell code) has one cell
size, dx = (xmax - xmin)/nx_p on each axis (`os-simulation.f03`, line 330).
It has no stretched or refined grid. Its grid scaling is load balancing. It
divides the grid among its processors so that each holds an equal share of
a load, the particles in its cells or a load density given as an expression
(`grid{ lb_type, spatial_loaddensity }`). The cuts come from recursive
bisection of the running load total (`rec_tree_dist`,
`os-grid-parallel.f90`, line 236). A range of n parts is cut where the load
reaches its left end's value plus nl/n of its own, nl = n/2, and each side
is cut again. On a continuous load that puts the k-th cut at the k/N
quantile of the total. `SOURCES.md` section 6 has the lines.

## What 1Dv0_3 takes from it

The rule, applied to the cells instead of the processors. `mesh.py` puts
the N + 1 edges at the k/N quantiles of a load, so every cell holds an
equal share of it and the cells are small where the load is dense. The load
is the number of cells wanted per unit length,

```text
w = 1/h(d),   h(d) = h0 + s d,   d = the distance to the nearer wall,
```

so the cells grow geometrically from each electrode, the k-th h_wall r^k
with r = e^s. Given the cell count and the wall cell, s is solved so that
the load totals N cells. OSIRIS sets no limit on how fast its ranks' widths
change, but a mesh needs one, and a mesh whose cells would grow by more
than 1.3 a cell is refused. `refine` doubles the cells on the same load.
Every edge stays and each cell splits at its load midpoint, so successive
meshes are nested.

The load sits at the walls because the important interactions of the gap
are there. The interior is v0_2's 0-D run (1Dv0_2's finding), while the
sheath, a few Debye lengths thick (about 1 um in the pulse), held all of
1Dv0_2's trouble. On 32 to 256 uniform cells of 312 to 39 um it was one
drained cell at the 40-eV ceiling, and the cells beside it did not
converge.

| mesh | cells | at the walls | in the middle | grows by |
|---|---|---|---|---|
| quick | 64 | 2 um | 895 um | 1.218 |
| default | 128 | 0.95 um | 470 um | 1.103 |
| refined | 256 | 0.46 um | 241 um | 1.050 |

The default is the quick mesh refined once (`--wall-cell
0.9508138901170884`, which rebuilds its edges exactly), so the three stay
nested.

## How it works

`mesh.py` places the cells. `transport.py` moves the charged species and
the electron energy on them with Fluid1D's scheme, cell by cell. Fluid1D
reads one dx in every Scharfetter-Gummel flux, divergence and substep limit
and in its semi-implicit Poisson solve (`Fluid1D/transport_kernel.py`), and
it is never edited. `transport.py` is that kernel's substep loop with the
mesh in place of dx. A face flux takes the centre-to-centre spacing, a
cell's update its own width, a face average weighs the two cells by their
half-widths (two half-cells' conductances in series, to first order), the
field solve's wall faces sit half a wall cell from the first centres, and
the substep limit is each cell's own, on its narrowest length. It calls
Fluid1D's own geometry-free kernel functions for the rest. As a derivative
of Fluid1D's kernel it carries the AGPL. On a uniform mesh every one of
these is Fluid1D's arithmetic exactly (the weights are 0.5 and 0.5, and 0.5
a + 0.5 b is 0.5 (a + b) to the bit), and the step is Fluid1D's to the bit.

`neutrals.py` moves the gas after each transport step. It solves the 1-D
Navier-Stokes equations of the mixture on the mesh. Every species and the
electrons' energy ride with the gas, F2 and F diffuse in the krypton (whose
flux balances theirs, so diffusion carries no mass), heat conducts into the
electrodes, held at tg0 (300 K), and viscosity slows the flow. Its
coefficients are Chapman-Enskog's on Lennard-Jones potentials, with
Svehla's force constants for Kr, F2 and F (`SOURCES.md` section 7). The
gas is 99.5 % Kr, so its viscosity and conductivity are krypton's. The
scheme is finite volumes, MUSCL-Hancock with van Leer's limiter and the
HLLC Riemann solver, all explicit and conservative. Each species' count,
the gas's energy and the electrons' energy are kept to round-off but for
the heat the electrodes take, which the books carry. Sound crosses the
narrowest cell in about a nanosecond, so one gas substep per 50-ps step
is enough in every run here.

`sim1d.py` builds Fluid1D's simulation as 1Dv0_2 does, then puts the mesh
under its transport, its field solver and every total, which sums over the
cells' own volumes (the input, the radiation, the ledger, the gap means, a
power profile's mean of 1). On a graded mesh Fluid1D's one cell size and
cell volume are unset, so any Fluid1D code that would read them fails
loudly, and the three books Fluid1D keeps with one volume are redone per
cell. The gas step comes after the transport in each step and hands the
gas's densities, temperature and velocity back to the chemistry.

## What it gives

v0_2's reference point: 1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us,
seed n_e 1e15, a 1-cm gap, dt = 50 ps, four numba threads (several runs at
a time on a shared machine, so the times are rough), in percent of the
given power, which is all the input. The field's work is kept to its
energy (the default).

| | uniform 64 | uniform 128 | still, 64 | still, 128 | still, 256 | moving, 64 | moving, 128 | moving, 256 |
|---|---|---|---|---|---|---|---|---|
| 248 nm, the gap | 4.528 | 4.598 | 4.621 | 4.608 | 4.593 | 4.642 | 4.632 | 4.618 |
| the same at dt = 25 ps | | | | | | 4.643 | 4.632 | |
| 248 nm, the centre cell | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 |
| 147 nm, the gap | 19.249 | 19.642 | 19.834 | 19.747 | 19.682 | 19.704 | 19.589 | 19.525 |
| the field's work on the ions | 5.71 | 2.07 | 0.37 | 0.28 | 0.21 | 0.50 | 0.49 | 0.29 |
| the same at dt = 25 ps | | | | | | 0.52 | 0.50 | |
| cut by the mean-energy bounds | -2.27 | -0.39 | 0.17 | 0.18 | 0.19 | 0.10 | 0.14 | 0.14 |
| conducted into the electrodes | | | | | | 0.052 | 0.049 | 0.050 |
| the walls' hottest gas, K | 1317 | 1030 | 1777 | 1957 | 1901 | 1234 | 1225 | 1180 |
| wall time | 141 s | 1409 s | 312 s | 938 s | 5192 s | 368 s | 781 s | 3486 s |

The centre cell is v0_2's 0-D run (4.66672 %) on every mesh, with either
gas and either field work. On a graded mesh the walls are resolved. The
first 7 to 12 um at each wall are the drained sheath, at the mean-energy
ceiling, and the bounds cut only what those cells take of the uniform
input. The field's work moves 0.2 to 0.5 % of the given power from the
electrons to the ions near the walls, whose share heats the gas there ("The
field's work", below). The walls' hottest gas is the hottest cell within
200 um of either wall in the profiles saved every 100 ns. With the gas
moving the graded meshes converge in dt to 0.02 %, and the gap's 248-nm
share falls as the walls are resolved, 4.642, 4.632 and 4.618 % with the
gas moving and 4.621, 4.608 and 4.593 % with it still. The moving gas
raises it by about 0.5 %, and every graded mesh agrees with every other to
1.1 %. The uniform meshes' wall cells hold the whole sheath, whose work on
their electrons is more than they hold, so the floor puts the difference
back there (the negative cut), and that costs the stage its substeps
(Limits, below).

With Fluid1D's field work (`--field-work fluid1d`, 1Dv0_3 before
2026-10-06) the same runs gave these. Its uniform columns are 1Dv0_2's own
numbers.

| | uniform 64 | uniform 128 | still, 64 | still, 128 | still, 256 | moving, 64 | moving, 128 | moving, 256 |
|---|---|---|---|---|---|---|---|---|
| 248 nm, the gap | 4.521 | 4.594 | 4.637 | 4.608 | 4.591 | 4.662 | 4.637 | 4.616 |
| the same at dt = 25 ps | | | 4.645 | 4.607 | | 4.661 | 4.645 | |
| 248 nm, the centre cell | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 | 4.667 |
| 147 nm, the gap | 19.391 | 19.703 | 20.039 | 19.802 | 19.724 | 19.970 | 19.705 | 19.576 |
| field work given | 2e-8 | 2e-4 | 1.09 | 0.033 | 5e-5 | 3.08 | 0.36 | 1e-4 |
| the same at dt = 25 ps | | | 2.58 | 0.095 | | 5.32 | 0.81 | |
| cut by the mean-energy bounds | 3.10 | 1.54 | 0.12 | 0.20 | 0.21 | 0.08 | 0.15 | 0.16 |
| conducted into the electrodes | | | | | | 0.095 | 0.046 | 0.045 |
| the walls' hottest gas, K | 302 | 650 | 2562 | 1949 | 1950 | 2810 | 1375 | 1182 |
| wall time | 81 s | 100 s | 309 s | 960 s | 4959 s | 429 s | 803 s | 3365 s |

Its field work given was the coupled step's error ("The field's work",
below). Its books also keep the transport step's own mean-energy bounds
in the walls' line, as 1Dv0_2's do, where the conserving books put them
in the cut. The run's summary then lists each energy's share of everything
put in, the given power and that work, and where the work is percents
(the 64-cell meshes) those shares differ from the ones here. Read per the
given power these graded meshes converge in dt to 0.2 % (the 2026-10-01
(b) entry's "4.587 to 4.528 %" on 64 cells was the field work's growth in
the denominator) and agree with each other to 1.5 %.

## What the resolved walls show

Inside the drained sheath, from about 10 um to 200 um, the graded meshes
show what the uniform cells held in one or two. The F2 burns out from the
wall inward during the pulse, a band at a time. With the gas still, on 128
cells the cell at 29 um has no F2 left by 300 ns (3.6e13 cm^-3 against the
centre's 6.5e16), the cell at 75 um by 500 ns, and those at 116 and 159 um
by 600 ns. A burned band stops attaching electrons, its n_e rises to 100
times the centre's, and its gas ends near 2000 K while its neighbours keep
their F2. This is the halogen-depletion (attachment) instability of
discharge-pumped excimer lasers, seeded at the wall. On 256 cells the bands
sit at 28, 73, 119 and 163 um and peak at 1850 and 1950 K nearest the walls
and 1710 and 1430 K further in, so the structure converges in place. These
are the runs with Fluid1D's field work. With the work kept to the field's
energy the still gas's bands on 128 cells burn out at the same times and
peak within 70 K of them (1957 K at 29 um, 1686 K at 75 um), and on 256
cells the four sit in the same places and peak 10 to 50 K cooler (1840,
1900, 1690 and 1410 K).

What heats a band is the energy the charges bring it. With the gas still
and Fluid1D's field work, the hottest band cells on 128 and 256 cells
receive none of it, yet their gas takes 2.4 to 2.8 times their given power
as heat. Only the charges move in the still gas, so the extra is the
ionization energy of the ions and electrons that drift and diffuse into a
burned band and recombine there. With the work kept to the field's energy
the hottest cell on 128 cells still takes 2.80 times its given power, of
which the ions bring 0.22 from the field.

With the gas moving the bands are cooler, since a hot band now expands and
conducts. On 256 cells they sit where the still gas put them, at 25, 69,
119 and 163 um, and peak near 1.1 us at 1050 to 1180 K, against the still
gas's 1410 to 1950 K, with n_e 10 to 15 times the centre's, while the
200 um at each wall hold up to 11 % more pressure than the centre. On 128
cells the two nearest each wall form, at 26 and 75 um, peaking at 1190 and
1230 K, and two weaker ones at 116 and 159 um (980 and 820 K). On 64 cells
the two nearest form at 24 and 80 um, at 1230 and 1080 K. So the bands'
places converge, the nearest band's peak agrees to 60 K from 64 to 256
cells and the second's to 150 K, and on every mesh the walls hold up to 10
to 11 % more pressure than the centre and the gas moves at 15 m/s at
most. These are the runs with the field's work
kept to its energy. Fluid1D's field work gives the same bands on 256 cells
(to 15 K) and on 128 cells the second one 150 K hotter, but on 64 cells it
made them hotter instead. There its near-wall cells took the field work
given, 3 to 4 times their given power, and the expanding gas fed it (N
falls, so E/N and the electrons' mobility rise). One layer 15 to 65 um
from each wall reached 2000 to 2800 K and 290 kPa against the centre's
211 kPa, and drove a compression wave into the gap at up to 73 m/s. That
layer was the step's error given as input ("The field's work", below), and
with the work kept to the field's energy it is gone.

Three causes were tried with the gas still and ruled out (2026-10-01 (b),
CHANGELOG):

- the given power in the drained sheath. With none within 10 um of the
  walls the bounds cut nothing, and the burned cells move outward.
- Fluid1D's field work, -e Gamma E of a cell's mean flux and mean field,
  replaced by the average of its two faces' products and still given as
  input where positive. The near-wall cells ran away sooner (the gap's mean
  n_e 2.8e16 cm^-3 by 1.2 us). The conserving work below also takes the
  faces, but with its sign and kept to the field's energy.
- v0_2's power per unit volume. Shared in proportion to n_e instead, the
  gap contracts at once (the gap's mean n_e 3.7e16 cm^-3 by 200 ns on 64
  cells, against 9e13).

## The field's work

The electrodes are grounded and nothing outside drives the gap, so the
field can only give the charges what its own energy loses. From a start
with no field, its net work on all of them over the run is minus its energy
at the end, about 1e-6 % of the input here.

Fluid1D's step (`field_work="fluid1d"`) takes the work on the electrons as
-e Gamma E of each cell's mean flux and mean field, gives its positive part
to the stage as input, drops its negative part and leaves out the ions. A
copy of the kernel that also booked the work on every charge measured the
books to 600 ns on 64 cells (2026-10-01 (c)).

| percent of the given power, to 600 ns | gas still | gas moving |
|---|---|---|
| on the electrons, as booked | 0.260 | 1.763 |
| on the ions | 0.037 | 0.114 |
| on all charges, the cells' means | 0.297 | 1.877 |
| on all charges, the faces' products | 0.69 | |

All of it is positive, where the gap allows none, so the field work given
was the coupled step's error, never a physical input. It sat within 200 um
of the walls, fell as they were resolved (1.09, 0.033 and 5e-5 % of the
given power with the gas still, 3.08, 0.36 and 1.2e-4 % with it moving),
grew as dt fell, and on 64 cells the moving gas fed it, as above.

Since 2026-10-06 (`field_work="conserve"`, the default) `transport.py`
keeps the work to the field's energy. Gauss's law in each cell and each
species' continuity across its faces give, between grounded electrodes,
the discrete energy theorem of each substep h,

```text
W' - W = -h sum_f J_f (E_f + E'_f)/2 l_f,    W = sum_f eps0 E_f^2 l_f / 2,
```

where J_f is the current through face f that moved the charge, E_f and
E'_f are the exact fields of the charge before and after, and l_f is the
face's length (the centres' spacing, half a wall cell at the walls). Each
substep solves those exact fields, one tridiagonal solve, takes each
face's work on the electrons and on the ions apart, and splits it between
the face's two cells by the half-cell each holds of its length (all of a
wall face's to the wall cell). The electrons' share goes into their energy
in the transport step, with its sign, the ions' share heats the gas, and
the books carry the field's energy. None of it is input. Fluid1D's species
update sets a density that goes negative to zero, which makes particles
and charge, and on 48 cells from 10 um that missed the theorem by 44 % of
the work. The conserving step instead scales what a cell's faces would
carry out over a substep beyond its content, each face by its upwind
cell's factor, so the counts stay non-negative and their continuity exact.
With both, the theorem holds to round-off.

On the graded meshes the electrons give the field 0.21 to 0.52 % of the
given power and the ions take it (the table above). The field's energy at
the end is about 1e-6 % and the theorem's residual, summed over the run,
3e-11 % at most. No reference run needed an outflow scaled, so there the
species update is Fluid1D's (the tests' 48 cells from 10 um scale 303 in 40
ns). Within 10 um of the walls, the sheath, the ions take 0.007 to 0.010 %
on every graded mesh. The rest is at the bands' edges, where n_e jumps up
to 160, 45 and 20 times from one cell to the next on 64, 128 and 256 cells,
and it falls as the edges are resolved (0.50, 0.49 and 0.29 % with the gas
moving, 0.37, 0.28 and 0.21 % with it still). Over a run the electrons lose
and the ions gain in every cell, and within one substep some cells'
electrons gain far less than the others lose. The moving gas carries charge
too, which changes the field's energy, by -2.5e-6 of the input over the
default run (the field's work on all charges is -2.5e-6 there, and the
field ends with 1e-8). That is booked with the energy the transport moves
and kept as `field_moved`, and with it the field's own books close to 3e-13
of the electrons' work (5e-4 without).

## Limits

1Dv0_2's, except the wall cells. They are the gas temperature transport
reads (tg0), Maxwellian electron transport, first order in dt, the mean
energy held to 0.03-60 eV, the seed (1e13 cm^-3 or more at these powers),
no discharge, and R11p held at its 30-eV value above 30 eV
(`1Dv0_2/README.md`). And:

- The 200 um at each wall, above. From 128 to 256 cells the gap's 248-nm
  share moves 0.3 % (4.632 to 4.618 %) and the interior is the 0-D run,
  but 128 cells resolve two of the four bands at each wall, and the
  field's work at the bands' edges falls as they are resolved (0.49 to
  0.29 % of the given power).
- The field's work, above. It is kept to the field's energy only between
  grounded electrodes with nothing driving the gap, and
  `field_work="conserve"` refuses a bias, a drive current or a gap
  capacitance. Where a wall cell holds the whole sheath, as on the uniform
  mesh, the sheath's work on its electrons is more than they hold, and the
  mean-energy floor puts the difference back (2.3 % of the given power on
  64 uniform cells, 0.4 % on 128, none on the graded meshes, whose wall
  cells are thin and hold few electrons). The stage then takes up to
  Fluid1D's 8000 substeps a step to raise those electrons' energy again,
  so a run on 128 uniform cells takes 23 min, against 100 s with Fluid1D's
  field work.
- The gas. Krypton's transport stands for the mixture's, the excited
  species and the ions only ride with the gas (the ions' own drift and
  diffusion are Fluid1D's), the electrodes are held at tg0, and there is
  no thermal diffusion. Svehla notes that the (12-6) potential gives high
  viscosities at high temperature (`SOURCES.md` section 7).
- A run costs more. The explicit transport's substep follows the smallest
  cell's diffusion limit, h^2/2D, and the burned cells are stiff. On 128
  cells a 50-ps step takes about 320 transport substeps and a run 13 min,
  against 5 substeps and 100 s on 128 uniform cells with Fluid1D's field
  work.
- The transport runs the coupled semi-implicit step only. Fluid1D's
  frozen-field NumPy path, and a prescribed drive with it, read one dx and
  are refused.

## Files

```
1Dv0_3/
  mesh.py        the cells: OSIRIS's equal-load rule on a wall load,
                 uniform(), graded(), refine()
  transport.py   Fluid1D's coupled transport step on the mesh (AGPL), the
                 field's work kept to its energy, MeshPlasmaModel and
                 MeshPoisson
  neutrals.py    the moving gas: Navier-Stokes of the mixture, its
                 Chapman-Enskog coefficients, NeutralGas
  fluid1d.py     the bridge (1Dv0_2's)
  sim1d.py       Params1D (with the mesh and the gas), LightSimulation,
                 run(), Result1D (with the cells' widths)
  run.py         the command line (LIGHT/run_1Dv0_3.py), with --mesh,
                 --wall-cell, --max-cell, --neutrals and --field-work
  model.py ...   v0_2's chemistry and its 0-D model, unchanged
  tests/test_mesh.py       the mesh (NumPy only)
  tests/test_neutrals.py   the gas (NumPy and SciPy)
  tests/test_transport.py  the transport against Fluid1D's
  tests/test_1d.py         the 1-D model (Fluid1D's Python)
  tests/test_*.py          v0_2's 0-D tests, pointed at 1Dv0_3's modules
  runs/          run output (git-ignored)
```

## Tests

`test_mesh.py` (9 tests, NumPy only, any Python) holds the mesh to its
rule. The uniform mesh is Fluid1D's to the bit, the graded cells grow from
the wall cell by one ratio, every cell holds an equal share of the load,
OSIRIS's recursive bisection gives the same edges, the rule for any load
density gives the analytic ones, the face weights are the half-cells in
series, refinement is nested, a cap holds, and impossible meshes are
refused. `test_neutrals.py` (18 tests, NumPy and SciPy, any Python) holds
the gas to its sources and to exact solutions. The collision integrals are
a direct quadrature's, the viscosities Svehla's tables', the diffusion
Brokaw's worked examples and their measurements. Sod's shock tube follows the
exact Riemann solution, on uniform and graded cells, a species front
rides the exact contact, a contact at rest stays at rest, sound rings at c
and decays at Stokes and Kirchhoff's rate, heat spreads at kappa
k^2/(rho c_p) between adiabatic walls and into isothermal electrodes, a
trace species spreads as 2Dt, diffusion carries no mass and mixing
releases no heat, a closed box keeps every species and its energy, the
electrodes take the heat the books say, the substep is the tighter of the
sound and diffusion limits, and a uniform gas stays put to the bit.
`test_transport.py` (12 tests) holds the transport to Fluid1D's. From a
live state on a uniform mesh one advance gives Fluid1D's densities, energy,
field and ledgers to the bit. The field solve is Fluid1D's on a uniform
mesh and second order on graded ones (a uniform charge between grounded
walls), the kernel's solve is the NumPy one, the Scharfetter-Gummel fluxes
hold the Boltzmann equilibrium on a graded mesh, the substep is each cell's
own limit, and a graded step keeps every species' count and the electrons'
energy. Five more hold the field's work to its energy. The exact field
holds Gauss's law in every cell and its energy is the faces' sum. From a
live state with its sheaths formed one conserving advance keeps the faces'
work to the field's energy to round-off, the cells' shares sum to the
faces' total, and no density goes negative. Twin advances from one state,
conserving and Fluid1D's, move the charges alike to the bit and the
electrons' energy apart by their share of the work (to 2e-13). An ion
drifted out of cells faster than they hold it keeps its count with the
positive transport, the last cell giving the wall all it holds, where
Fluid1D's floor makes particles. `test_1d.py` adds to 1Dv0_2's tests. The
same short run on the uniform mesh with the gas still is 1Dv0_2's to the
bit, the input is the power over the gap's volume, the gap means and a
power profile weigh the cells by volume, the graded mesh's centre is the
0-D model, on the same 48 cells the graded mesh's walls take a small part
of what the uniform mesh's do, and the moving gas keeps the books and the
interior. With the field's work kept to its energy, on 48 graded cells
with the gas still and moving, none of it is input, the electrons' work,
the ions' and the field's energy close to round-off (less what the moving
gas changed), the ions' heat has its line in the fate, the walls take no
more than with Fluid1D's work, and the summary and `--field-work` say
which work a run kept. It also pins the Fluid1D functions transport.py
follows by hash. `LIGHT_1D_SLOW=1` adds the reference point on 32 uniform
cells and on the default mesh.

Three mutation checks planted bugs one at a time, each in a copy of its
own. In `mesh.py`, `transport.py` and `sim1d.py` the tests catch 19 of 20.
The 20th, a uniform mesh's dx read as edges[-1]/N where edges[1] is dx to
the bit, changes no result for any gap or cell count tried. In
`neutrals.py` and its coupling they catch 28 of 29. The 29th leaves
unbooked what Fluid1D's adoption changes in the energy the gas carried,
which is 5e-16 of the input in every run tried (the gas carries the
electrons' mean energy between its neighbours' values). In the field's
work kept to its energy (`transport.py`, `sim1d.py` and `run.py`) they
catch all 26. The first pass let one through, the ions' share booked but
left out of the gas, until the walls' loss was held below Fluid1D's, and
the twin advances and the wall cell's content were written for two bugs no
check then covered (a limited cell giving nothing, the electrons' step
leaving out their share).
