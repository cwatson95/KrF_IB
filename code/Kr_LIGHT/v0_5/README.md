# LIGHT_Kr v0_5: the pure-krypton model with the electrons on its excited atoms

v0_5 (2026-10-16) is v0_4 with the limitations of v0_4's accuracy
overview taken up. The owner asked: "Begin approaching the limitations of
v0_4 to make v0_5", and, shown the groundwork (`../CHANGELOG.md`,
2026-10-16 (a)), chose the scope (2026-10-16 (b)): published calculations
where nothing is measured, the pumps' energy pair by pair, the electrons'
collisions with the excited atoms and with Kr2*, two data fixes, and the
two-temperature closure's tail table given an excited-population axis;
no paywalled or bot-checked full texts.

What changed from v0_4:

- **Calculations where nothing is measured** (`excited.py`,
  `network.py`). The stepwise ionization of the four 5s levels
  (`KSTEP_S_<level>`) and of the 5p rung's ten levels (`KSTEP_P_<2p>`),
  held at 0 by v0_1 to v0_4, and the ionization of Kr2* (`EI_KR2`,
  `EI_KR2F`), not carried before, run on Vriens & Smeets's cross section
  (PRA 22, 940 (1980)) at each level's own ionization energy; the bands
  kstep_s, kstep_p and ei_kr2 are scales, 1 by default, 0 to 2. For Kr2*
  the formula's peak is Flannery et al.'s ~1e-15 cm2. e + Kr2* -> e + 2Kr
  (`EQ_KR2`) and the mixing of its two states (`EM_KR2`, `EM_KR2r`) were
  found neither measured nor calculated: bands at 0, to Kannari's
  estimate of the quenching, 1e-7 cm3/s.
- **The pumps' energy.** Each of Jung et al.'s 26 5s -> 5p
  pairs is its own row (`K76_<1s>_<2p>`), charged its own threshold on
  Biagi's level energies, the level's excess over 2p10 to heat; on average
  1.537 eV out of 1s5 where v0_4 charged 1.388.
  The reverses (`K76r_<1s>_<2p>`) leave from the rung's levels, which are
  populated at the gas temperature (`excited.rung_fraction`); the lines
  are NIST's one by one (`PRAD_<2p>_<1s>`, 28), and the
  rung's relaxation has its detailed-balance reverse (`KPSr_*`).
- **The 5s levels and the electrons.** Their superelastic return to the
  ground state (`SUPER_<level>`, the Klein-Rosseland reverse of each G1
  cross section) and 1s5 <-> 1s4 by electrons (`EMIX_1s5`, `EMIXr_1s4`,
  Blagoev & Pilosov's 4.3e-8 cm3/s; band emix_5s to 2.3e-7).
- **The kinetic EEDF** (`boltzmann.py`): every electron row on an excited
  population is a process of its own (74 processes), the
  stepwise ionizations sharing the energy between the two electrons and
  counting in the growth term.
- **The tail table's excited axis** (`tail_field.py`,
  `data/kr_tail_field.txt`, md5-pinned): r(n_e/N, Te, x_s), x_s the 5s
  fraction, 1666 points, the rung and Kr2* in steady state at
  each; its x_s = 0 block is v0_4's table. `Params.tail = 'field0'`
  reads that block alone.
- **Data**: Kr2+ bound by Signorell et al.'s 9267.8 cm-1; Kholin et al.'s
  1s3 rate cited with its argon third body.
- **A kinetic table point with no solution no longer stops the run**
  (`model._kinetic_update`): after the restart and v0_4's fallbacks it is
  left out of the table (`dropped`), and an update whose centre has none
  leaves LIGHT's Maxwellian until the next renewal (`failed`). The scan's
  0.25-atm, 500-kW corner needed it (n_e 1e17 cm-3 at 2243 K in its
  afterglow): 11 points dropped, no failed update.

```bash
python3 LIGHT_Kr/run_v0_5.py --eedf two_t           # the reference point, seconds
python3 LIGHT_Kr/run_v0_5.py                        # on the kinetic EEDF, about 21 minutes
python3 -m v0_5.tail_field 10                       # from LIGHT_Kr: regenerate the tail table on 10 processes (pin its md5)
/usr/bin/python3 -m pytest LIGHT_Kr/v0_5/tests -q   # from Kinetic_Modeling
```

## At the reference point

LIGHT's drive in pure krypton: 1 atm at 300 K, 300 kW/cm3 for 1 us into
the electrons, a 1e15 cm-3 seed, 2 us in all.

| % of the input at 2 us | kinetic EEDF (default) | two-temperature | v0_4 kinetic | v0_4 two-temperature |
|---|---|---|---|---|
| 147 nm (Kr2*) | 6.08 | 5.58 | 31.22 | 32.92 |
| near IR (5p -> 5s) | 0.48 | 0.48 | 0.47 | 0.45 |
| heat, reactions and relaxation | 90.69 | 90.67 | 65.45 | 63.16 |
| heat, elastic collisions | 3.23 | 3.75 | 3.21 | 3.81 |
| gas temperature at 2 us | 856 K | 859 K | 707 K | 697 K |

At the reference point the electrons ionize 56.9 % of the
excimers formed (28.5 % radiate). Counted by number, the
excimer's ionization makes 45.2 % of the ions, the 5s levels'
33.8 % and the 5p rung's 16.0 %, and every
ion recombines onto the 5p rung. The electrons' collisions with Kr2* decide
the answer: without its ionization the share is 15.43 %, at
twice the formula 4.47 %; its quenching and the mixing of
its two states, at Kannari's estimate of 1e-7 cm3/s, give
1.27 and 12.44 %.

What each change does, undone alone (kinetic EEDF; v0_5 6.08 %):

| undone | 147 nm |
|---|---|
| the three calculated rows at 0 (v0_4's rule) | 28.02 |
| no stepwise ionization (5s and 5p) | 7.50 |
| no ionization of Kr2* by electrons | 15.43 |
| no superelastic return of the 5s levels | 6.11 |
| no 1s5 -> 1s4 by electrons | 6.07 |
| no reverse of the 5p rung's relaxation | 6.07 |
| every level forms the excimer at Ueda's 1.7e-32 (v0_3's) | 3.16 |
| no redissociation of Kr2* | 6.42 |
| the J = 1 levels pumped by no electron | 5.55 |

Propagated uncertainty of the kinetic 147-nm share: +-7.0
points in quadrature (+-115 %); the largest term is
the ionization of Kr2* by electrons.

## The scores

- Swarm (v0_3's cross sections, unchanged): train W 2.8 %
  rms, train alpha/N (50-300 Td) 9.9 % rms, alpha/N below
  50 Td mean -36 %.
- Electroluminescence (Mano et al. 2022): 0.99 to 1.04 of the
  fit from 6 to 13.5 Td.
- Superelastic return (new, Blagoev et al. 1984, thermal electrons):
  5.3-6.6 times the measurement for 1s5,
  5.5-9.5 for 1s4, 0.5-1.2
  for 1s3.
- Beam: Morozov +29 % (v0_4 +29 %), Eckstrom
  at 1 atm -3 % (v0_4 +1 %).
- The two-temperature closure against the kinetic EEDF:
  -8 % at the reference point (on v0_4's table, without
  the excited axis, +4 %), within 6 % at
  the scan's corners.

## Checks

- `/usr/bin/python3 -m pytest -p no:cacheprovider LIGHT_Kr/v0_5/tests`
  from Kinetic_Modeling: 417 passed, 1 skipped without test_docs.py (1 h 52 min), and test_docs.py with
  LIGHT_KINETIC_SLOW=1 7 passed (10 min: both generators reproduce their numbers files). The gate, `tests/test_v0_5.py`:
  everything v0_5 does not change is v0_4's; Vriens & Smeets's cross
  section integrates to their rate formula; the pseudo cross sections
  carry their rate laws exactly; every population operator exchanges
  rate x energy with the electrons; the pairs sum to Jung's level rates
  and balance in detail; the rung's partition; the KPS reverse in hot gas;
  Signorell's Kr2+; the Blagoev ratios pinned; the table survives a point
  or a centre with no EEDF.
- The overview's generator raises when a claim its text makes fails.

## Files

```
LIGHT_Kr/
  run_v0_5.py         launcher
  v0_5/
    excited.py        Vriens & Smeets, the pseudo cross sections, the 5p rung at Tg
    network.py        the rows (219), the bands (CALCULATED and OFF_BY_DEFAULT apart)
    boltzmann.py      the population processes (one direction, ionization)
    jung2006.py       Jung et al.'s pairs, their thresholds and reverses
    tail_field.py     the tail table's generator, with the excited axis
    model.py, report.py, run.py, species.py, thermo.py, xsec.py, swarm.py, swarm_mc.py, ib.py
    data/kr_tail_field.txt          the field tail table (md5-pinned; tail_field.py writes it)
    data/Kr_Biagi897_LXCat.txt, kr_swarm_trained.json, kr_tail_ratio.txt   (v0_3's)
    tests/            test_v0_5.py (the gate), test_carried_v0_4.py and v0_4's tests carried over
  docs/
    accuracy_overview_v0_5.tex/.pdf, system_map_v0_5.tex/.pdf
    make_overview_v0_5.py, make_system_map_v0_5.py, kr_docs_v0_5.py
```
