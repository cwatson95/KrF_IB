# LIGHT v0_3

v0_3 is v0_2 with the kinetic EEDF in the loop. The owner asked on
2026-09-30 for v0_3 as v0_2 with ML added, and chose the route that takes
ML's Boltzmann solver into the run. The chemistry, the bands and the input
are v0_2's (`LIGHT/v0_2/README.md`). What changes is where the electron
rates come from. v0_2 integrates every electron rate over a Maxwellian at
its Te with a tabulated tail. v0_3 solves the two-term Boltzmann equation
with electron-electron collisions during the run, at LIGHT's own mean
energy, and takes every electron rate, the elastic loss and the energy a
captured electron carries from that solution.

`docs/accuracy_overview_v0_3.pdf` scores v0_3 against the measured data,
beside v0_2.

```bash
python3 LIGHT/run_v0_3.py                        # the reference point, several minutes
python3 LIGHT/run_v0_3.py --eedf two_t           # v0_2's closure, which is v0_2
python3 -m pytest LIGHT/v0_3/tests -q            # the tests, about 8 minutes
LIGHT_KINETIC_SLOW=1 python3 -m pytest LIGHT/v0_3/tests/test_kinetic.py -q   # and the convergence runs
LIGHT_DOCS_SLOW=1 python3 -m pytest LIGHT/v0_3/tests/test_docs.py -q         # and the overview's slow numbers
```

## The kinetic equation

`boltzmann.py` is ML's scheme (`ML/ml/boltzmann`, the website LXCat build's
two-term solver and ML's e-e term) written again on sparse matrices. It
reads LoKI-MC's cross sections the way LIGHT's `xsec.py` does, and it
carries exactly the electron processes the basic network books, so its
energy balance and LIGHT's ledger close the same budget.

- **Collisions.** Kr's excitations and ionization, and F2's attachment,
  dissociation and vibration, which are the network's cross-section rows
  (G1, C1a-c, C2, ATT, DIS_a3, DIS_A1, VIB1-4), with Kr's elastic loss
  (G6). F2's other processes and its elastic collisions, which no row
  books, stay in its momentum transfer and nowhere else.
- **Sources and sinks.** Penning ionization (G2, PEN_KR2) injects its
  electrons at the energy the network gives them, 5.834 and 5.554 eV.
  Recombination (R11p) removes them with sigma v ~ eps^n, the template of
  its Te^n law. The growth term carries the net, so every electron is
  accounted for.
- **Electrons born with a spectrum.** A drive that injects electrons gives
  them through `Model.kinetic_injection`, as a rate with a spectrum linear
  on [0, top] of a given mean. The prototype beam drive
  (`docs/survey_2026-09-27/scripts/ebeam_screen.py`) gives its secondary
  electrons this way, one per ion pair, on [0, 9.915] eV with a mean of
  4.506 eV, as ML's `solve_beam` injects them. A given power injects none
  and heats as a field.
- **Electron-electron collisions.** ML's Landau term, with its integrals
  taken on ln F linear between cell centres. A Maxwellian of any
  temperature is then left alone on the grid, as the exact operator leaves
  it. ML's quadrature (F constant in each cell) misses that at the lowest
  boundaries by more than the elastic coupling to the gas, and an undriven
  EEDF condenses into the first cell.
- **Its own powers.** The fluxes through the cell boundaries are
  Scharfetter-Gummel fluxes. The elastic and field powers LIGHT reads are
  the parts of those fluxes, the central W F - D dF/deps with the diffusive
  part times (z/2) coth(z/2), so they are the energy the solver itself
  exchanges.

## How LIGHT reads it

Every solution has LIGHT's mean energy, (3/2) Te, held by a field E/N >= 0.
In the pulse the field stands for the input power. After it, the field
stands for the energy LIGHT's electrons still carry, which is the local
mean energy approximation of fluid models with LIGHT's own energy equation
choosing the mean. Where even the no-field solution is hotter than LIGHT's
electrons (the Penning electrons or a beam's secondaries alone heat it),
the field is 0. The field is found inside every step of the iteration,
since at a fixed field the problem can have two solutions.

LIGHT reads a table of the kinetic rates, the captured energies and the
elastic power at five mean energies, Te e^r for r = -0.16, -0.08, 0, 0.08
and 0.16, interpolated in ln Te. Past the table's ends it carries them with
its own Maxwellian dependence on Te, a rate by its logarithmic slope, a
captured energy by its ratio and the elastic power by its difference, which
passes through 0 at the gas temperature. The table is renewed when Te has
moved 0.12 in ln Te, when n_e, the F2 fraction, Tg, the Penning or the
recombination frequency per electron or the Kr* fraction has moved 0.10 in
its ln (each above a floor below which it cannot matter), and at least
every 10 ns in the pulse. After the pulse that interval grows to a tenth of
the time since the pulse ended. Below n_e = 1e3 cm^-3 the EEDF is not
solved and the Maxwellian closure stands in.

In the pulse the EEDF is advanced in time (`Params.kin_memory_s`, the whole
pulse by default since 2026-10-01). Each update is one backward-Euler step
from the last update's EEDF over the time since it, to LIGHT's mean energy,
`boltzmann.Solver.step`. The growth term there is the population's own
change, nu from F, with no rule to drop it, since in a step it moves the
shape by nu dt at most. Where the step with no field is hotter than LIGHT's
electrons, the no-field step stands, unreached. In Morozov's low-pump
pure-Kr run under the beam that is 52 of the pulse's 139 steps, 0.11 to
0.23 % above LIGHT's Te (solved steady, 56 updates of its first 2 ns
missed, by up to a factor 115), and the VUV yield moves by 3.6e-4 of
itself. The run starts from the gas's Maxwellian, with LIGHT's mean energy
on the grid, under LIGHT's own closure until the first step. A history
moved onto another grid keeps its mean energy (`Grid.remap`), since a step
reads that mean as its start.

After the pulse the EEDF is solved steady. The growth term keeps its shape
steady while n_e changes. It is dropped where it cannot give LIGHT's mean
energy, where the energy it moves outweighs the heating (field and sources)
ten times over while the gas relaxes the shape less than ten times faster
than the population turns over, and where it has no positive solution.
Without it the first cell takes up the net gain or loss of electrons. Where
even its no-field solution is hotter than LIGHT's electrons, that solution
stands, unreached, when the population grows or the gas relaxes the shape
ten times faster than the turnover (a recombining afterglow). A solution
without the term reaches a colder mean there only by putting the lost
electrons back cold. A cold start tries the growth term first, from the net
of the sources, on a budget.

## At the reference point

1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us.

| Model | 248 nm | 272 nm | 275 nm | 420 nm | 147 nm | all bands | Tg end |
|---|---|---|---|---|---|---|---|
| v0_2 | 4.67 % | 0.149 % | 0.712 % | 4.87 % | 20.0 % | 30.4 % | 649 K |
| **v0_3** | **4.79 %** | **0.153 %** | **0.735 %** | **4.96 %** | **19.7 %** | **30.3 %** | **650 K** |

The kinetic EEDF raises the laser line by 2.7 % of itself and lowers the
147-nm band by 1.8 %. Its electrons spend less on elastic heating, 5.3 % of
the deposited energy against v0_2's 6.2 %, and on ionization, 26.0 %
against 26.8 %, and more on Kr excitation, 99.4 % against 98.3 % (gross
flows, since the Penning electrons bring part of it back to the electrons).
KrF(B) forms from 24.6 % of the deposited energy, against 23.8 %.

- **The two energy balances agree.** In the pulse's second half the field
  brings 1.006 to 1.031 times the input power. The excess is the table's
  renewal. Between renewals LIGHT runs on a table made at the last
  renewal's electron density and gas, and n_e grows by up to 5.7 % between
  renewals there. So at each renewal LIGHT's electrons sit off the fresh
  table's balance by the excess (their n_e d<eps>/dt is -0.7 to -3.2 % of
  the input power there), and they regain it with an e-folding time of
  about 0.2 ns (at 721 ns, -2.8 % falls to -0.12 % in 1.5 ns). Renewed
  every 5 ns, the excess is at most 1.015. With LIGHT's own n_e d<eps>/dt
  at the renewal, less the step's own (-0.03 to -0.21 %), the two balances
  agree to 2.9e-5 of the input power. (Until 2026-10-01 this README read
  the excess as the electrons cooling through the pulse. They do cool, but
  by 6 to 44 times less.)
- **The two particle balances agree.** The kinetic growth frequency equals
  LIGHT's d n_e/dt / n_e wherever the growth term is kept.
- **Converged.** Renewing the table at half the drift moves no band by more
  than 3.9e-5 of itself, renewing it every 5 ns by 6.8e-5, and doubling the
  cells to 1200 by 1.6e-4.
- **Its cost.** 355 kinetic solutions. 61 of them come in the first 46 ps,
  where Kr2+, the ion R11p recombines, builds up from nothing and its
  recombination frequency climbs through the trigger's floor (1e4 per
  second) a tenth at a time. Solved steady they took 114 s of a 391-s run
  (2026-09-30). As time steps they are cheap, and the first 2 ns take 28 s
  against 111 s solved steady (2026-10-01). A floor tied to the solution's
  own frequencies would save most of them.

## At Tahoma

Under the prototype beam drive, in ML's Tahoma configuration (v0_2's
default network with every band at its default and n_e free from a seed of
1e10 cm^-3, 665 Torr with 2.8 Torr of F2, 330 kW/cm^3 for 450 ns), v0_3
gives g0 = 5.64 %/cm, x0.88 the measured 6.39 (v0_2 gives 5.69, x0.89),
Isat = 1.56 MW/cm^2 (x1.11) and g0 Isat/Pin = 0.266 (x0.99). Over the 200
to 400 ns window its Te, n_e, gain and formation are those of ML's
self-consistent LIGHT-Boltzmann loop to 0.08, 0.47, 0.08 and 0.07 %, where
v0_2's closure is off by 0.7, 2.5, 1.0 and 0.5 %. LIGHT's mean energy is
the kinetic one at every update in the window, held by a field of at most
0.13 Td, and the two energy balances agree to 6.0e-4 of the beam's electron
power.

## Found on the way

1. **At a fixed field the problem has two solutions.** 0.12 ns into the
   pulse, at 0.4 Td, a hot start gives a mean energy of 0.282 eV and a cold
   start 0.042 eV. Kr's Ramsauer minimum under strong e-e collisions makes
   the mean energy against E/N an S curve. At a fixed mean energy the
   solution is one, and the field it needs falls from 0.444 Td at 0.06 eV
   to 0.388 Td at 0.2 eV before it rises again. So the field is found
   inside the iteration, and a search on E/N (the first design) could not
   converge there.
2. **ML's e-e quadrature condenses a cold EEDF.** With no drive, pure Kr at
   300 K and n_e/N = 1e-5 went to a mean energy of 0.029 eV on 600 cells
   and 0.007 eV on 2400, falling with the cell width. With the integrals on
   ln F linear between cells it is the Maxwellian at 300 K to 1e-14.
3. **Factors over a Maxwellian fail where the two mean energies part.** The
   first coupling carried the kinetic rates after the pulse as factors over
   the Maxwellian at the kinetic mean energy, applied at LIGHT's Te (ML's
   method). 325 ns into a short run the kinetic EEDF's temperature was 0.29
   eV and LIGHT's 1.0 eV, the elastic factor was 22 and F2's dissociation
   factor 6e4, and LIGHT's Te fell from 1.1 to 0.41 eV in 0.12 ns. The
   table at LIGHT's own mean energy replaced it.
4. **The grid.** A top scaled by the decades F falls across it oscillated
   between 41 and 13 eV, since the fall is slow in the bulk and steep past
   the Kr thresholds, and the 13-eV grid, with Kr's ionization threshold
   off its end, was returned. The top is now read off F.
5. **A beam's secondaries belong in the kinetic equation.** A field carries
   a given power, but a beam's secondaries are born with the sub-excitation
   spectrum, as ML's beam-sustained solver injects them. With them in, at
   Tahoma's start the sources outrun every sink in a seed of 1e10 cm^-3.
   The solution without the growth term then does not exist (the first cell
   would have to take up a net gain), and the growth term started from nu =
   0 lost positivity at once. A cold start now takes nu from the net of the
   sources, a positivity reset raises a gain, and the no-field growth
   solution stands where it is hotter than LIGHT's electrons.
6. **A table held flat past its ends stalls LIGHT's Te.** At Tahoma LIGHT's
   Te sat at 1.3 eV where the kinetic mean energy gives 1.6. Below the
   table's bottom the rates and powers were held at the bottom's values,
   which left LIGHT's energy equation no loss that falls with its Te, so
   its Te crept toward the kinetic one at the electrons' turnover rate
   alone. Past its ends the table now follows LIGHT's own Maxwellian.
7. **A pure-Kr afterglow renewed the table 2000 times.** The table was
   renewed at least every 10 ns, so a 20-us afterglow took 2000 renewals,
   and the triggers watched Te, n_e and the F2 fraction only, so the
   interval was the only guard on the other inputs. The triggers now watch
   every input of the kinetic state, and after the pulse the interval grows
   with the time since the pulse ended. Morozov's low-pump run takes 235
   renewals.
8. **A trigger crossed at t = 0.** Under the beam the secondaries heat a
   seed of 1e10 cm^-3 so fast that Te passes its trigger 1.6e-16 s into the
   run, inside the 8.9e-16 s to which scipy (1.13) locates an event (brentq
   with an absolute tolerance of 4 machine epsilons). At Tahoma the event
   came back at 3.0e-16 s and the run moved on. In pure Kr at 300 kW/cm^3
   it came back at t = 0 itself, and the run renewed its table there 12,405
   times without moving. A segment whose event lands where it began now
   takes a short step without the triggers, doubled while it repeats.
9. **The elastic power was not the solver's.** The table's elastic power
   was taken by central differences, F + (kTg/e) dF/deps at the cell
   boundaries, while the solver's fluxes are Scharfetter-Gummel. On 600
   cells to 12 eV a cell is 0.77 kTg at 300 K, and central differences gave
   the solver's own equilibrium a loss to the gas of 4.7 % of its cooling
   term. Late in Morozov's low-pump afterglow LIGHT read that loss and its
   Te fell to 0.930 of the gas temperature, below the kinetic solution's
   0.98. The powers are now the solver's own flux split, and LIGHT's Te
   there is the kinetic one, 0.977 of the gas temperature at 20 us (the
   grid's own thermal bias, below).
10. **The late afterglow went below the gas temperature twice more.** The
    elastic power past the table's bottom was carried as a ratio of the
    Maxwellian's powers, which divides by a number passing through 0 near
    the gas temperature. It is now carried by their difference. And
    dropping the growth term in a recombining afterglow put the recombined
    electrons back cold, in the first cell. The term now stays where the
    gas relaxes the shape ten times faster than the population turns over.
11. **Rows dropped after construction.** ML drops rows from a model it has
    built (`run_beam`'s drop, its formation fits) and rescales others, and
    the kinetic equation read its processes once, when the model was built.
    It now reads them at the start of every run. The table's Maxwellian
    values at its ends looked rows up by their place in the list, which a
    dropped row shifts, and now look them up by row (found reading the
    code, before any run used it).
12. **No steady EEDF holds the first nanosecond.** At the reference point
    the electrons heat from the gas temperature through the first
    nanosecond, and most of the input power raises their mean energy.
    Solved steady, the EEDF dropped the growth term at 86 of its 88 updates
    there, its field carried 0.02 to 0.19 of the input power, and its
    energy balance missed LIGHT's by up to 0.22 of it. A backward-Euler
    step from the last update's EEDF (option 2 of
    `LIGHT/improvement_for_EEDF.txt`, the memory term) holds it. Its field
    carries 0.98 to 1.07 of the input power, 0 updates drop the growth
    term, and after 0.1 ns the two energy balances agree to 1.5e-3. Marched
    instead with the physical drive, the Joule power set to the input power
    and the mean left free (option 3, run offline), the EEDF's mean energy
    stays within 1.7 % of LIGHT's (1.1 % after 0.1 ns), so holding LIGHT's
    mean loses nothing there. The bands barely move. Steady against
    stepped, no band moves by more than 1.9e-4 of itself, and LIGHT's own
    Maxwellian closure in the first nanosecond moves none by more than
    9.6e-4, since only 1.0e-3 of the pulse's energy enters then. At Tahoma
    a steady start moves the window's gain by 3.6e-4 and n_e by 9.2e-4.
13. **The flux scheme's energy split is coarse at thermal energies.** Where
    a cell is as wide as kT (thermal electrons on 600 cells to 12 eV), the
    e-e part of the Scharfetter-Gummel flux takes 12 % of a 1e9 eV/s
    field's power in a picosecond's step from the gas's Maxwellian, and the
    step's energy balance holds only to 1.0e-1 (7.6e-3 on 2400 cells,
    2.8e-5 with e-e left out). Rescaling the e-e drift to force its own
    energy exchange to zero did nothing (lambda 0.999994). The error lies
    in splitting a flux whose c(z) couples its parts, not in the e-e term
    itself, and only resolution cures it.
14. **A time step in the afterglow is not yet the answer.** Run through the
    afterglow of Morozov's low-pump pure-Kr run as well, the steps left 87
    of its 95 afterglow updates short of LIGHT's mean energy, the no-field
    step up to 8 % hotter than LIGHT's Te (LIGHT's Te falls faster than a
    step with no field can follow), and moved the VUV yield by 4.6e-4
    against the steady run. The afterglow stays steady until that is
    understood.

## Checks

`tests/test_kinetic.py` and `tests/test_docs.py`, with v0_2's own tests (on
v0_2's closure, as `tests/test_model.py` says) beside them.

- **ML's solver.** With ML's e-e quadrature the solver is ML's dense one
  (`ee.solve_field`) to 2.2e-9 in F and 4.2e-8 in every rate, at three
  states with e-e and the growth term (the test holds them to 1e-7 and
  1e-6). v0_3's own quadrature moves F by 1.1e-5 at most at those driven
  states.
- **Balances.** The growth frequency is the solution's own net frequency,
  and the energy the field and the sources bring equals what the
  collisions, the captures and the growth term take, to 2.0e-6 in the pulse
  and 5.6e-5 with no field. Without the growth term the first cell's
  injection closes it to 2.9e-4, a residual that falls with the cell width.
  The elastic and field powers are the parts of the solver's own fluxes,
  and the undriven equilibrium exchanges no energy with the gas to 1e-12 of
  its cooling term.
- **The limits.** No drive gives the Maxwellian at the gas temperature. At
  the transient's state the fixed-field problem has its two solutions and
  the fixed-mean problem one, from three different starts. The growth term
  is kept in the pulse and dropped where it should be, with the heating of
  Penning electrons and of a beam's secondaries alike. A recombining
  afterglow keeps it, from a cold start and from one without it, though the
  solution without it could reach LIGHT's colder mean. The grid comes out
  the same from a top far too long and one far too short.
- **Beam electrons.** The injected spectrum reaches the grid with its mean
  to 1.6e-5 on grids of 12 to 40 eV. A seed under the beam grows on its
  secondaries. The beam drive's run carries them, closes its ledger, holds
  LIGHT's mean energy after 10 ns and agrees with LIGHT's energy balance to
  6.5e-4 of the drive's electron power. A run steps past a trigger at its
  start.
- **The table and the rows.** Past its ends the table follows LIGHT's
  Maxwellian, every rate, the elastic power and the captured energies, for
  a steady table and a time step's. Rows dropped after construction leave
  the kinetic equation.
- **The time step.** A step too long to matter is the steady growth
  solution, and one too short to matter, with no field, leaves F where it
  was. With e-e left out its energy balance closes to 1e-4 at 1 and 100 ps,
  and its field takes a set Joule power to 1e-9. Where the target is colder
  than the step can reach, the step with no field stands, converged. A
  history keeps its mean energy on another grid. The run's first nanosecond
  reaches LIGHT's mean energy at every update with the growth term on, and
  after 0.1 ns its energy balance agrees with LIGHT's to 2e-3.
- **The run.** A 0.3-us run closes its ledger, every EEDF has LIGHT's mean
  energy, and the two energy balances agree to 1e-3. Under v0_2's closure
  v0_3 gives v0_2's numbers exactly. `LIGHT_KINETIC_SLOW=1` adds the
  reference point's convergence.
- **The overview.** `test_docs.py` holds every number the overview prints
  to its generated files, the files to a fresh run (`LIGHT_DOCS_SLOW=1` for
  the slow ones), its words to its numbers, and ML's Tahoma configuration
  to the one it names.
- **Mutations.** Twelve bugs planted in `boltzmann.py` one at a time
  (2026-09-30 (d)), and the fast tests fail on eleven. The survivor removes
  the guard that stops the grid search from revisiting a grid, which the
  tests never reach. Thirty-four more in the code that scored v0_3
  (2026-09-30 (e)), and the tests fail on thirty. Two survivors change only
  the cost, and two guard fallbacks no test reaches (`LIGHT/CHANGELOG.md`).
  Twelve more in the time step (2026-10-01), and the fast tests fail on all
  twelve. One, the no-field step where a target is out of reach, failed
  only once its check moved from a 1-ps step, where the field's floor exit
  leaves F 3e-12 off the converged step, to a 0.1-ns one (5e-7).

## Limits

- After the pulse the EEDF is a steady solution, while the electrons' shape
  is set by their history, and v0_3 takes the steady shape at their mean
  energy. In the pulse it follows the history (found 12, and 14 for why the
  afterglow stays steady).
- After the pulse, where the growth term is dropped, the net loss of
  electrons comes back in the first cell. At low field with strong
  attachment that cools the solution, and since v0_3 holds LIGHT's mean
  energy regardless, it moves only the field and the shape below a few
  tenths of an eV.
- Thermal electrons. On the default 600 cells the grid's Maxwellian at 300
  K has a mean energy of 0.0377 eV, below the 0.0388 eV of the gas, so late
  in an afterglow LIGHT's electrons sit a few percent below the gas
  temperature (0.977 of it at the end of Morozov's low-pump run). 2400
  cells give 0.0386 eV. There a time step's energy exchange holds only to
  1.0e-1 on the default cells (found 13).
- QE and QEC, the electron quench of KrF(B, C), keep their constant rate,
  and their energy goes to heat.
- A run takes minutes against seconds for v0_2, so the bands' scans and the
  scan over the notes' ranges (`--bands`, `--scan`, `--omissions`) take
  hours under the kinetic closure. `--eedf two_t` is fast.
- ML runs v0_3 when `LIGHT_VERSION=v0_3` names it. The overview runs ML's
  forward model, its formation observables, its kinetic loop and its
  pure-Kr check on it. ML's fits have not been rerun on it.
