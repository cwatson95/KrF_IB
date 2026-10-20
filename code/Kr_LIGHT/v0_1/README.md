# LIGHT_Kr v0_1: the pure-krypton model

v0_1 (2026-10-13) is LIGHT's v0_5 for a fill of pure krypton. It is one
well-mixed box (0-D) of krypton into whose electrons a given power is put,
and it says where that energy goes: the 5s and 5p manifolds, the ions, the
excimer Kr2* and its 147-nm second continuum, the near-infrared 5p -> 5s
lines, heat, and the electron gas. The owner asked for it on 2026-10-13:
"Make v0_1 of pure Kr from the remaining v0_5 directory".

v0_1 is v0_5 with every fluorine species, row, band and parameter taken
out, and the 248-nm optics with them. What is left is v0_5's krypton,
unchanged row for row, and v0_5's machinery: the energy ledger, the
kinetic EEDF in the loop (the default closure), LIGHT v0_2's
two-temperature closure (seconds), the walls, the EM wave with Turnbull's
inverse bremsstrahlung, and the bands. One thing is not v0_5's: two
entries of jung2006's tables, corrected to the paper (`SOURCES.md`
section 6).

```bash
python3 LIGHT_Kr/run_v0_1.py --eedf two_t           # the reference point, seconds
python3 LIGHT_Kr/run_v0_1.py                        # on the kinetic EEDF, minutes
python3 LIGHT_Kr/run_v0_1.py --eedf two_t --scan    # pressure x power
python3 LIGHT_Kr/run_v0_1.py --eedf two_t --bands   # each band's two ends
/usr/bin/python3 -m pytest LIGHT_Kr/v0_1/tests -q   # from Kinetic_Modeling
```

## The network

Eighteen rows, each v0_5's (`network.py`, sources in `SOURCES.md`):

| Rows | What |
|---|---|
| G1, C1a-c, C2 | electron-impact excitation of the 5s manifold (Kr*) and of the 5p rung (Kr**, the higher groups cascading through it), and ionization, on Biagi-v7.1 |
| G6 | elastic loss to the gas |
| G2 | Kr* + Kr* Penning ionization (band kr_pen) |
| K15, R11p | Kr+ + 2Kr -> Kr2+ + Kr, and Kr2+ + e -> Kr** + Kr |
| AK8_Kr, KR2_RAD | Kr* + 2Kr -> Kr2* + Kr (band kr2_form), and Kr2* -> 2Kr + hv at 147 nm (270 ns) |
| PEN_KR2 | Kr2* + Kr2* Penning (band kr2_pen, 0: no measurement) |
| KPS, PRAD | Kr** + Kr -> Kr* + Kr (band kps_kr), and Kr** -> Kr* + hv in the near IR (29 ns) |
| K76, K76r | e + Kr* <-> e + Kr** on Jung et al.'s measured cross sections (band k_5s5p) |
| KSTEP_S, KSTEP_P | stepwise ionization of Kr* and Kr** (bands at 0: no measurement) |

The optional group `electron` (`--add electron`) adds REBUILD's stepwise
template (KSTEP) and the superelastic return of G1 (SUPER). n_e is always
solved: no electrons are injected, and n_e = Kr+ + Kr2+.

## At the reference point

LIGHT's drive in pure krypton: 1 atm at 300 K, 300 kW/cm3 for 1 us into
the electrons, a 1e15 cm-3 seed, 2 us in all.

| Fate at 2 us, % of the input | kinetic EEDF (default) | two-temperature |
|---|---|---|
| 147 nm (Kr2*) | 14.06 | 15.05 |
| near IR (5p -> 5s) | 0.65 | 0.63 |
| heat, reactions and relaxation | 83.13 | 81.81 |
| heat, elastic collisions | 2.70 | 3.03 |
| gas temperature at 2 us | 808 K | 802 K |

Under the kinetic EEDF the field carries the input power to 1e-3 and the
kinetic and LIGHT's energy balances agree to 7e-6 in the pulse's second
half; the ledger closes to 4e-15.

The electrons pump Kr* to the 5p rung (K76) faster than three-body
collisions turn it into Kr2*, and krypton relaxes the 5p rung back (KPS)
with the 1.385-eV gap as heat. Under the two-temperature closure the
electrons pump Kr* to the 5p rung 32 times for each Kr2* formed, and that
relaxation heat is 76.5 % of the input, 43 eV for each 147-nm photon. So
the band on the measured 5s -> 5p cross section decides the answer:

| k_5s5p | 0 (no cycle) | 0.5 | 1 (Jung et al.) | 1.4 |
|---|---|---|---|---|
| 147 nm, % of the input (two-temperature) | 51.6 | 20.1 | 15.05 | 13.0 |

Kannari's fit in its place (`--k76 kannari`) gives 16.3 %. Pressure and
power move it as much (`--scan`, two-temperature, % of the input at 147
nm):

| kW/cm3 | 0.1 atm | 0.25 atm | 0.5 atm | 1 atm | 2 atm |
|---|---|---|---|---|---|
| 100 | 0.2 | 1.6 | 8.3 | 27.2 | 59.3 |
| 300 | 0.1 | 0.5 | 3.3 | 15.1 | 43.4 |
| 500 | 0.0 | 0.3 | 2.0 | 10.5 | 35.5 |

At 0.1 atm a microsecond of this drive heats the gas to 2000-8800 K,
beyond the measured range of several rates (the report warns).

## Checks

`tests/test_v0_1.py` holds the conversion: the species, rows and bands are
v0_5's krypton ones field for field, no fluorine is left, jung2006 differs
from v0_5's by its two corrections alone, the kinetic EEDF table and the
right-hand side from one state are v0_5's bit for bit and to round-off,
and whole runs agree with v0_5 (no F2 in its fill, the same tables) to the
solver's tolerance in nine configurations, the difference shrinking a
hundredfold from rtol 1e-8 to 1e-10. The rest are LIGHT's tests carried to
krypton, each assertion measured again on pure-Kr states: `test_model.py`,
`test_kinetic.py`, `test_rung_walls_wave.py` (LIGHT v0_4's 5p rung, walls,
wave and row_scale), `test_jung2006.py` (v0_5's cross sections, now
against the paper's tables entry by entry), `test_network.py`,
`test_xsec.py` and `test_ib.py`; `test_docs.py` holds the accuracy
overview and the system map (`LIGHT_Kr/docs`) to the model. From
Kinetic_Modeling on 2026-10-13: 346 passed, 4 skipped (the slow runs,
LIGHT_KINETIC_SLOW=1) and 1 expected failure (the long-grid search below)
in 10 minutes. With LIGHT_KINETIC_SLOW=1 the slow tests pass too: the
documents' generators reproduce their files exactly, and (33 minutes) the
kinetic runs of v0_1 and v0_5 agree within the closure's own spread, and
at the reference point half the renewal step, a renewal every 5 ns or
twice the energy cells moves neither emission band by 5e-4 of itself.

## Found on the way

- `Params.row_scale` reached no row the kinetic table carries (G1, C1a-c,
  C2, R11p, SUPER) under the kinetic EEDF: the table's rate was read
  without the row's factor. Found by the accuracy overview (R11p at +-10 %
  moved the kinetic run by 1e-10) and fixed 2026-10-13: the factor goes
  into the row's cross section or recombination template, as K76's band
  always did (`test_row_scale_reaches_the_kinetic_table`). A row_scale on
  G6 is refused, and K76 and K76r take one factor (2026-10-14). LIGHT v0_4
  and v0_5, LIGHT_ArF v0_3 and v0_4 and LIGHT_Kr's v0_5 copy carried the
  defect too and were fixed the same way on 2026-10-14.
- jung2006 carried two misreadings of Jung et al. 2006 (`SOURCES.md`
  section 6): 21 for Table II's 23 (1s5 -> 2p6, 4 eV), and Table IV's
  1s5 -> 2p5 booked as a core-changing 1s3 -> 2p5. v0_1 corrects both;
  LIGHT/v0_5 still carries them (147 nm at the reference point 15.066 %
  with them, 15.051 % without). LIGHT's test of the tables compared
  through pytest's default absolute tolerance, which passes any two cross
  sections of 1e-16 cm2, and its expectation for the 1s5 -> 2p8
  extrapolation took the 8-eV point where the column ends at 10 eV.
- The kinetic closure's own spread: the same run at rtol 1e-8 and 1e-9
  differs by 3e-4 on the 147-nm share and 1.2e-3 on the elastic share (the
  EEDF table is renewed where the solver's steps put the drift triggers).
  v0_1 and v0_5 differ by that much: by the same two numbers on the
  same tables, and within three times the spread with v0_1's corrected
  ones (the slow test).
- Pure krypton takes the kinetic solver where LIGHT's Kr/F2 states did
  not. The S curve of the Ramsauer minimum sits near 0.04 Td (in LIGHT's
  state at 0.4 Td, where F2 held a cold branch): a fixed-field iteration
  there finds the upper solution or none, while the mean match gives one
  solution from any start. From a grid far too long (40 or 60 eV at the
  plasma state of the tests) the EEDF's tail falls to the double-precision
  floor and the grid search keeps the long grid, its ionization rate 31 %
  off (a strict xfail test). Runs start at 30 eV and warm-start; at the
  reference point no solution met that, and the steady ones outside the
  search's window (35 of 105) are on its 12-eV minimum grid, where cold
  electrons fall 25 to 34 decades.

## Not in v0_1

- Processes v0_5 never carried, each with more room in a fill without F2
  (no harpoon removes Kr*, no attachment removes the electrons), and each
  needing a measured rate before it can enter: Kr3+ (Kr2+ + 2Kr) and its
  recombination; electron collisions with Kr2* (ionization, quenching,
  mixing of its two states, which v0_1 carries as one with the slow
  component's 270 ns); Penning ionization between Kr* and Kr2*, and of
  Kr**; Kr2+ and Kr2* redissociation in hot gas; the resonance lines of the
  two J = 1 5s levels.
- The two-temperature closure reads LIGHT's kinetic tail table at its F2 =
  0 column; the table was generated (ML/ml/tail.py) with LIGHT's
  2026-09-28 network under the prototype e-beam drive. The default
  closure, the kinetic EEDF, does not use it.
- A pure-Kr measurement of v0_1's own regime, a field-heated gas. The
  accuracy overview (`LIGHT_Kr/docs/accuracy_overview_v0_1.pdf`) scores
  its parts: the solver against krypton's swarm data (drift velocity 9.6 %
  rms; alpha/N 6.5 % at 100-300 Td, -54 % below 50 Td) and the chain to
  147 nm under an e-beam (+31 % against Morozov et al. 2008 at low pump;
  at 3e5 W/cm3 and 1 atm 37 % below its ceiling to Penning ionization,
  where Eckstrom et al. 1988 report little loss). Its propagated
  uncertainty on the reference point's 147-nm share is +-85 %, led by Kr*
  + 2Kr. The system map (`system_map_v0_1.pdf`) draws the reference
  point's energy flow.

## Files

```
LIGHT_Kr/
  run_v0_1.py         launcher (puts LIGHT_Kr/ on sys.path)
  v0_1/
    species.py        the six krypton species (v0_5's values)
    network.py        the 18 rows, the optional group, the bands
    model.py          the 0-D ODE (densities, U_e, Tg, ledger), Params, run()
    report.py         flows, fate, reservoirs (Kr*, Kr**, Kr2*), warnings
    run.py            command line: single point, --bands, --scan, --omissions
    xsec.py           LXCat parser, Maxwellian tables, the two-temperature tail
    boltzmann.py      the kinetic EEDF (two-term Boltzmann with e-e)
    jung2006.py       the measured 5s -> 5p cross sections (two entries corrected)
    ib.py             the EM wave's inverse bremsstrahlung (Turnbull)
    data/             Kr_Biagi71_LXCat.txt and kr_tail_ratio.txt, md5-pinned
    tests/            the tests above
    SOURCES.md        the literature, with how each value was read
```
