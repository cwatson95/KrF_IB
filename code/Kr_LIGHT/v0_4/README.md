# LIGHT_Kr v0_4: the pure-krypton model with its 5s levels resolved

v0_4 (2026-10-15) is v0_3 with section 7 of v0_3's accuracy overview
taken up. The owner asked: "Begin approaching the limitations of v0_3
such as from section 7 of the accuracy overview to begin making the
improved v0_4", and, shown the groundwork (`../CHANGELOG.md`, 2026-10-15
(b)), chose the scope: the 5s levels resolved (each forming the excimer
at its own measured rate), redissociation by detailed balance, bands for
the excimer's two states, its quenching and the recombination's
gas-temperature factor, and numerics: a field-driven tail table for the
two-temperature closure, Mano et al.'s electroluminescence yield as a new
score, the afterglow's bound, and the 451-K and first-nanosecond
statements corrected.

What changed from v0_3:

- **The four 5s levels** (`species.py`, `network.py`). Kr*(1s5),
  Kr*(1s4), Kr*(1s3) and Kr*(1s2) at Biagi's thresholds, each made by its
  own G1 row. v0_3's one Kr* sat at the 1s5 energy.
  - The excimer forms out of 1s5 at 3.81e-32 and out of 1s4 at 9.24e-32
    cm6/s (Marchal et al. 2018, selective excitation in pure Kr; bands
    kr2_form 1.7-4.2e-32, kr2_form_1s4 1.27-9.24e-32); out of 1s3 at the
    1s5 rate as a proxy (band 0 to 9.24e-32). v0_3 used Ueda's 1.7e-32.
  - 1s4 -> 1s5 on Kr at Turner's 1.3e-13 cm3/s at 295 K (through his
    196-K point), 1s2 -> 1s4 three-body at Audouard's 0.2e-32, both with
    detailed-balance reverses (MIX_*).
  - The 5s -> 5p pair per level (K76_<level>, K76r_<level>) on Jung's
    sets, the J = 1 levels on their same-core metastable's (band
    k_5s5p_j1); per level they average to v0_3's lump by weight.
  - PRAD branches by NIST's A-values weighted by 2J + 1, KPS by the
    levels' weights (band kps_5s); Penning per pair of levels (G2_a_b);
    stepwise ionization per level (KSTEP_S_<level>, band 0).
  - The resonance lines (RES_1s4, RES_1s2; band res_escape, 0).
- **Redissociation by detailed balance** (`thermo.py`). AK8r_* and
  K15r: the forward rate over the equilibrium constant on the measured
  wells (Kr2* 0.709 eV, band kr2_d0 to Killeen & Eden's 0.62; Kr2+ 1.15
  eV). At the fill's density the excimer's bound share is
  27 times its atoms' at 700 K and 1 at 988 K.
- **Excimer and ion bands.** The fast 0u+ state (AK8f_1s4, KR2F_RAD,
  3.4 ns; band kr2_fast 0-0.17), Kr2* + Kr quenching (QKR2; band
  kr2_quench 0-8.6e-15), the recombination's gas-temperature factor
  (band dr_tg 0-1, Cunningham & Hobson over Shiu & Biondi). The report
  lists these "off by default" apart from the unmeasured bands.
- **The two-temperature closure's tail** (`tail_field.py`,
  `data/kr_tail_field.txt`, md5-pinned). r(n_e/N, Te) from this model's
  own kinetic EEDF for a field-driven discharge, every point from a cold
  start; `Params.tail = 'beam'` reads v0_3's table. Along the kinetic
  reference run it gives 0.82-0.85 of the kinetic
  EEDF's own tail ratio (v0_3's table 0.95-1.30), and
  0.73-1.24 over the scan's corners: its states
  carry no excited atoms, whose 5s <-> 5p collisions shape the run's tail
  (hotter than the table's at the reference point and at 2 atm).
  Never tune it toward the kinetic run (`../CHANGELOG.md`, Known-bad).
- **`Params.kin_afterglow_memory`** time-steps the kinetic EEDF after the
  pulse too: it moves the reference point's 147-nm share by 2e-5.
- **`swarm.el_yield`**: the swarm's electroluminescence yield against
  Mano et al. 2022 (1.1 bar), a score nothing trained on.
- **`boltzmann.py`**: where a kinetic table point told to leave the
  growth term out has no steady solution (the afterglow of very hot gas,
  where Penning ionization outruns recombination), `Solver.eedf` takes
  the growth term on the same grid, then a grid of its own
  (`Solver.fallbacks`). Before, the run stopped (the scan's 0.25-atm,
  500-kW corner).
- **`run_v0_4.py`**: `--tail beam`, `--afterglow-memory`. No optional
  group ('electron' booked the lumped Kr*).

```bash
python3 LIGHT_Kr/run_v0_4.py --eedf two_t           # the reference point, seconds
python3 LIGHT_Kr/run_v0_4.py                        # on the kinetic EEDF, minutes
python3 -m v0_4.tail_field                          # from LIGHT_Kr: regenerate the tail table (pin its md5)
/usr/bin/python3 -m pytest LIGHT_Kr/v0_4/tests -q   # from Kinetic_Modeling
```

## At the reference point

LIGHT's drive in pure krypton: 1 atm at 300 K, 300 kW/cm3 for 1 us into
the electrons, a 1e15 cm-3 seed, 2 us in all.

| % of the input at 2 us | kinetic EEDF (default) | two-temperature | v0_3 kinetic | v0_3 two-temperature |
|---|---|---|---|---|
| 147 nm (Kr2*) | 31.22 | 32.92 | 14.05 | 15.78 |
| near IR (5p -> 5s) | 0.47 | 0.45 | 0.64 | 0.63 |
| heat, reactions and relaxation | 65.45 | 63.16 | 83.09 | 81.37 |
| heat, elastic collisions | 3.21 | 3.81 | 2.75 | 2.74 |
| gas temperature at 2 us | 707 K | 697 K | 808 K | 798 K |

What each change does, undone alone (kinetic EEDF; v0_4 31.22 %):

| undone | 147 nm |
|---|---|
| every level forms the excimer at Ueda's 1.7e-32 (v0_3's) | 11.92 |
| 1s4 forms it at the 1s5 rate | 23.05 |
| Kr** + Kr all to 1s5 | 27.52 |
| no redissociation of Kr2* | 31.68 |
| no redissociation of Kr2+ | 31.22 |
| no transfer between the 5s levels | 31.17 |
| the J = 1 levels pumped by no electron | 54.81 |

Over the run 5.8 % of the excimers formed redissociate;
52.6 % form out of 1s4. Propagated uncertainty of the
kinetic 147-nm share: +-20.1 points in quadrature
(+-64 %); the largest term is the J = 1 levels' pumping.

## The scores

- Swarm (v0_3's cross sections, unchanged): train W 2.8 %
  rms, train alpha/N (50-300 Td) 9.9 % rms, alpha/N below
  50 Td mean -36 %.
- Electroluminescence (new, Mano et al. 2022): 0.99 to 1.04 of
  the fit from 6 to 13.5 Td.
- Beam: Morozov +29 % (v0_3 +29 %), Eckstrom
  at 1 atm +1 % (v0_3 -26 %).
- The two-temperature closure against the kinetic EEDF: +5 %
  at the reference point (v0_3's table +8 %), within
  5 % at the scan's corners.

## Checks

- `/usr/bin/python3 -m pytest -p no:cacheprovider LIGHT_Kr/v0_4/tests`
  from Kinetic_Modeling: 403 passed, 1 skipped without test_docs.py (29 min), and test_docs.py with LIGHT_KINETIC_SLOW=1 7 passed (10 min: both generators reproduce their numbers files). The gate, `tests/test_v0_4.py`:
  everything v0_4 does not change is v0_3's; the 5s levels are Biagi's
  and each G1 row makes its own; the per-level 5s -> 5p rates average to
  v0_3's lump; NIST's branching; Turner's transfer and its reverse; the
  equilibrium constant against the classical integral; the reverses
  balance their forwards; the closure reads the field table; the
  afterglow's time step; Mano's yield; the solver's restart and fallback.
- The overview's generator raises when a claim its text makes fails.

## Files

```
LIGHT_Kr/
  run_v0_4.py         launcher
  v0_4/
    species.py        the four 5s levels, the two excimer states
    network.py        the rows (126), the bands (OFF_BY_DEFAULT apart)
    thermo.py         equilibrium constants for the redissociations
    tail_field.py     the field-driven tail table's generator
    model.py, report.py, run.py, boltzmann.py, jung2006.py, xsec.py, swarm.py, swarm_mc.py, ib.py
    data/kr_tail_field.txt          the field tail table (md5-pinned; tail_field.py writes it)
    data/Kr_Biagi897_LXCat.txt, kr_swarm_trained.json, kr_tail_ratio.txt   (v0_3's)
    tests/            test_v0_4.py (the gate) and v0_3's tests carried over
  docs/
    accuracy_overview_v0_4.tex/.pdf, system_map_v0_4.tex/.pdf
    make_overview_v0_4.py, make_system_map_v0_4.py, kr_docs_v0_4.py
```
