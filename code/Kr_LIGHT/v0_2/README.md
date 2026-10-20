# LIGHT_Kr v0_2: the pure-krypton model, its electron collisions level by level

v0_2 (2026-10-14) is v0_1 with the electron's collisions on ground-state
krypton booked level by level. The owner asked for it on 2026-10-14:
"LIGHT_Kr appears to have an error of -54% in some cases, let us proceed
with making v0_2 to improve this". The -54 % is v0_1's ionization
coefficient below 50 Td in its accuracy overview, against Kruithof 1940
(-73 %) and Heylen 1971 (-33 %). Offered the choice, the owner took the
level-resolved cross sections with the Hornbeck-Molnar ionization as a
band, and both documents.

What changed from v0_1, and nothing else did:

- **Cross sections.** v0_1 read Biagi-v7.1, which groups krypton's
  excitations into four effective levels (5s, 5p, the 4d and 6s group,
  5d and above) at the groups' lowest thresholds. v0_2 reads Biagi's
  level-resolved set (Magboltz 8.97, LXCat's 'Biagi' database, the
  owner's copy of 2026-06-18): 51 levels, each its own row at its own
  threshold (`network.py`). The rows keep their groups' names: `G1_<1s>`
  the four 5s levels into Kr*, `C1a_<2p>` the ten 5p levels into Kr**,
  `C1b_*` the 4d and 6s levels and `C1c_*` every level from 12.76 eV up,
  cascading to Kr**; each pays its own threshold and heats the gas with
  its excess over its product's energy. The species sit at the file's
  thresholds, which are NIST's levels (Kr* 9.9152, Kr** 11.3035, Kr+
  13.9996 eV; v0_1 had 9.915, 11.30 and 13.996). The optional SUPER row
  reads the four 5s levels' summed cross section at the 1s5 energy
  (`xsec.lumped_5s`), as v0_1 read Biagi-v7.1's 5s group.
- **The Hornbeck-Molnar rows.** An atom excited at or above Kr2+'s energy
  (IP - D0 = 12.85 eV; 20 of the levels) can ionize on a ground-state
  atom, Kr(nl) + Kr -> Kr2+ + e. Each such level has a twin row `HM_<lvl>`
  that takes the share `hm_kr` of its excitations (the cascade keeps 1 -
  hm_kr). No measurement of that share exists for krypton, so the band is
  0 by default, 1 its top. The rows' rate is first order in Kr (the
  partner is in the share), and under the kinetic EEDF their electron is
  born at rest (`boltzmann.py`, collision kind 'assoc').
- **`swarm.py`**: the kinetic EEDF's solver run as a swarm (no e-e, the
  growth term on), with the band's share, which the accuracy overview
  scores against krypton's swarm data.
- **Speed.** The stiff solver gets the Jacobian's structure
  (`Model.ivp_options`: the ledger and the 85 event counters feed
  nothing) and the two-temperature closure computes its shared tail parts
  once per state. A two-temperature run takes 3.5 s (v0_1 2.3 s; 36 s
  without the two), the same numbers to 6e-9.

```bash
python3 LIGHT_Kr/run_v0_2.py --eedf two_t           # the reference point, seconds
python3 LIGHT_Kr/run_v0_2.py                        # on the kinetic EEDF, minutes
python3 LIGHT_Kr/run_v0_2.py --rate hm_kr=1         # every Hornbeck-Molnar level ionizing
python3 LIGHT_Kr/run_v0_2.py --eedf two_t --scan    # pressure x power
python3 LIGHT_Kr/run_v0_2.py --eedf two_t --bands   # each band's two ends
/usr/bin/python3 -m pytest LIGHT_Kr/v0_2/tests -q   # from Kinetic_Modeling
```

## The -54 %: what it is and what v0_2 does about it

- **Not the solver.** LoKI-MC, the vendored Monte Carlo solver (no
  two-term approximation), on the same cross sections gives an
  ionization coefficient the two-term solution matches to 3-10 % on v0_1's
  file and 2.5-4 % on v0_2's at 20-50 Td, inside the Monte Carlo's
  counting error; the mean energies agree to 0.2 %
  (`tests/test_v0_2.py`). ML's dense solver matches v0_2's sparse one to
  round-off on v0_2's file (`tests/test_kinetic.py`).
- **Every cross-section set does it.** Bordage et al., J. Phys. D 46,
  334003 (2013) compared every published Kr set: calculated alpha/N falls
  below the low-field measurements, "a general trend", which they ascribe
  to Penning ionization of impurities in the experiments. Kruithof's and
  Heylen's values differ by up to 2.2x near 30 Td, and Jacques et al.
  (1986) also measured below Kruithof at the low end.
- **The level-resolved set** moves the mean error below 50 Td from -54 to
  -47 % (Heylen -33 -> -23 %) and the drift velocity's rms from 9.6 to
  10.8 % (the newer set's momentum transfer).
- **The Hornbeck-Molnar band.** Its top raises alpha/N 17-fold at 14 Td,
  where the electrons' tail falls by orders of magnitude between the
  levels (12.86-13.6 eV) and the ionization threshold (14.0 eV), and
  1.7-fold at 100 Td. At its top 18 of the 21 measurements below 50 Td
  lie between the band's ends (Kruithof's three at and below 17 Td lie
  above even its top), but the same top puts the model's mean error at
  100-300 Td at +39 %, where direct ionization alone is -7 % (7.6 % rms):
  a share that large cannot hold at every field. The share depends on the
  pressure (association against radiation and quenching), and a Townsend
  experiment's low-field points come from higher pressures; no
  measurement in hand gives either for krypton, so it stays 0 by default
  under the predictive policy.
- **The answer does not move.** At the reference point direct ionization
  makes 0.01 % of the ions and Penning ionization the rest; hm_kr at its
  top moves the 147-nm share by -0.02 % (kinetic) and -0.01 %
  (two-temperature).

## At the reference point

LIGHT's drive in pure krypton: 1 atm at 300 K, 300 kW/cm3 for 1 us into
the electrons, a 1e15 cm-3 seed, 2 us in all.

| Fate at 2 us, % of the input | kinetic EEDF (default) | two-temperature | v0_1 kinetic | v0_1 two-temperature |
|---|---|---|---|---|
| 147 nm (Kr2*) | 14.20 | 16.00 | 14.06 | 15.05 |
| near IR (5p -> 5s) | 0.64 | 0.63 | 0.65 | 0.63 |
| heat, reactions and relaxation | 83.11 | 81.40 | 83.13 | 81.81 |
| heat, elastic collisions | 2.58 | 2.49 | 2.70 | 3.03 |
| gas temperature at 2 us | 807 K | 797 K | 808 K | 802 K |

Under the kinetic EEDF the field carries 1.0001-1.0002 of the input power
and the kinetic and LIGHT's energy balances agree to 9e-6 in the pulse's
second half (381 solutions, 1.3-15.5 Td); the ledger closes to 3e-15.

The electrons pump Kr* to the 5p rung (K76) faster than three-body
collisions turn it into Kr2*, and krypton relaxes the rung back (KPS)
with the 1.388-eV gap as heat: 34 pumpings for each Kr2* formed on the
kinetic EEDF, 29 on the two-temperature closure, where that heat is
75.2 % of the input, 40 eV for each 147-nm photon. So the band on the
measured 5s -> 5p cross section decides the answer:

| k_5s5p | 0 (no cycle) | 0.5 | 1 (Jung et al.) | 1.4 |
|---|---|---|---|---|
| 147 nm, % of the input (two-temperature) | 52.2 | 21.1 | 16.00 | 13.9 |
| 147 nm, % of the input (kinetic) | 50.9 | 18.8 | 14.20 | 12.4 |

Kannari's fit in its place (`--k76 kannari`) gives 17.5 % (two-temperature;
15.4 % kinetic).
Pressure and power move it as much (`--scan`, two-temperature, % of the
input at 147 nm):

| kW/cm3 | 0.1 atm | 0.25 atm | 0.5 atm | 1 atm | 2 atm |
|---|---|---|---|---|---|
| 100 | 0.2 | 1.6 | 8.7 | 29.2 | 61.9 |
| 300 | 0.1 | 0.6 | 3.4 | 16.0 | 46.1 |
| 500 | 0.0 | 0.3 | 2.1 | 11.0 | 37.7 |

At 0.1 atm a microsecond of this drive heats the gas to 2000-8800 K,
beyond the measured range of several rates (the report warns).

## Checks

`tests/test_v0_2.py` holds the change from v0_1: every row but the
levels is v0_1's field for field in v0_1's order, v0_1's bands are
unchanged with hm_kr the one more, jung2006's data and the tail table are
v0_1's, and the species sit at Biagi 8.97's thresholds; at hm_kr = 0 the
right-hand side is the one without the HM rows, number for number, and
the kinetic equation leaves them out; at hm_kr = 1 every such excitation
ionizes and the ledger closes on both closures; the kinetic equation's
associative ionization counts in the growth term and keeps the energy
balance; the swarm against LoKI-MC and the band's effect on it. The rest
are v0_1's tests carried to the levels: `test_network.py` (every level
one row at its own threshold, the groups, the twins), `test_xsec.py` (the
51 blocks, the lumped 5s manifold for SUPER), `test_model.py`,
`test_kinetic.py` (with ML's dense solver on v0_2's own file),
`test_rung_walls_wave.py` (row_scale on a level and on an HM row),
`test_jung2006.py`, `test_ib.py` and `test_docs.py` (the two documents).
From Kinetic_Modeling on 2026-10-14: 389 passed and 3 skipped (the slow
reference-point test and test_docs' two regeneration tests,
LIGHT_KINETIC_SLOW=1) in 14 minutes; with LIGHT_KINETIC_SLOW=1
`test_docs.py` passes all 7 (both generators reproduce their files byte
for byte).

## Found on the way

- The kinetic energy-balance gap of a short pulse is the grid's: 1.1e-5
  on 600 cells with 51 thresholds (v0_1 7.8e-6), 3.2e-6 on 1200 cells.
- v0_1's strict-xfail grid search (from a grid five times too long the
  tail reached the double-precision floor and the search kept the long
  grid) does not reproduce at that state on the 51 levels: from 60, 40 and
  30 eV the search lands on one 21.3-eV grid. The solver is unchanged, so
  the rule stays (runs start at 30 eV).
- One of the solver's mean matches lands 2.2e-9 from its target, inside
  its design (the inner match is held to 1e-3 of the outer change); v0_1's
  test asked 1e-9.

## Not in v0_2

- v0_1's list: Kr3+ and its recombination; electron collisions with
  Kr2*; Penning ionization between Kr* and Kr2*, and of Kr**; Kr2+ and
  Kr2* redissociation in hot gas; the resonance lines of the two J = 1 5s
  levels (v0_2 resolves their excitation but books them into the lumped
  Kr*).
- The Hornbeck-Molnar share as a measured rate: 0 until a krypton
  measurement of the levels' association against their radiation and
  quenching sets it.
- The two-temperature closure reads LIGHT's kinetic tail table at its F2
  = 0 column, made (ML/ml/tail.py) with LIGHT's 2026-09-28 Kr/F2 network
  on Biagi-v7.1's groups. The default closure does not use it.
- A pure-Kr measurement of v0_2's own regime, a field-heated gas. The
  accuracy overview (`LIGHT_Kr/docs/accuracy_overview_v0_2.pdf`) scores
  its parts; the system map (`system_map_v0_2.pdf`) draws the reference
  point's energy flow.

## Files

```
LIGHT_Kr/
  run_v0_2.py         launcher (puts LIGHT_Kr/ on sys.path)
  v0_2/
    species.py        the six krypton species, at Biagi 8.97's thresholds
    network.py        the 85 rows (51 levels, 20 HM twins), the optional group, the bands
    model.py          the 0-D ODE (densities, U_e, Tg, ledger), Params, run()
    report.py         flows, fate, reservoirs (Kr*, Kr**, Kr2*), warnings
    run.py            command line: single point, --bands, --scan, --omissions
    swarm.py          the solver as a swarm, with the Hornbeck-Molnar share
    xsec.py           LXCat parser, Maxwellian tables, the two-temperature tail
    boltzmann.py      the kinetic EEDF (two-term Boltzmann with e-e)
    jung2006.py       the measured 5s -> 5p cross sections (v0_1's)
    ib.py             the EM wave's inverse bremsstrahlung (Turnbull)
    data/             Kr_Biagi897_LXCat.txt and kr_tail_ratio.txt, md5-pinned
    tests/            the tests above
    SOURCES.md        the literature, with how each value was read (section 7 new)
```
