# LIGHT v0_5: the 5s -> 5p excitation measured

v0_5 is v0_4 with the one rate the 2026-10-09 accuracy overview named as
deciding the partition replaced by a measurement, and with the new
articles the owner added (`Kinetic_Modeling/new_articles/`) read against
the overview's Section 7 ("Where v0_4 is not yet reliable"). Everything
else is v0_4's (`v0_4/README.md`): the 5p rung, the 248-nm optics, KrF* +
F, the ion-ion limbs, the walls, the wave, the kinetic EEDF.

```bash
python3 LIGHT/run_v0_5.py --eedf two_t                  # the reference point (seconds)
python3 LIGHT/run_v0_5.py --k76 kannari --eedf two_t    # v0_4's fit in its place
python3 LIGHT/run_v0_5.py --bands                       # each band's two ends
LIGHT_VERSION=v0_5 python3 LIGHT/docs/krf_data/tahoma_v0_4.py [--kinetic]   # Tahoma g0, alpha0
LIGHT_VERSION=v0_5 python3 LIGHT/docs/krf_data/kleikamp_v0_5.py            # the Kr/F2 laser's optima
/usr/bin/python3 -m pytest LIGHT/v0_5/tests -q          # test_kinetic.py on its own
```

## What changed

1. **K76, e + Kr(5s) -> e + Kr(5p), is a measured cross section.** Jung,
   Stone, Boffard, Anderson and Lin measured the electron-impact
   excitation out of the two metastable levels of Kr into the ten 4p5 5p
   levels (Phys. Rev. A 73, 022722 (2006), Tables II-IV; the Letter Phys.
   Rev. Lett. 94, 163202 (2005)). `jung2006.py` carries the tables and
   builds the one cross section LIGHT's lumped Kr* -> Kr** needs: the
   statistical-weight mean over the 5s manifold of each level's sum over
   the ten 5p levels, with the two J = 1 resonance levels (nothing
   measured; radiation-trapped and populated with the metastables at 1
   atm) given the set of the metastable with the same ion core. Linear
   between the tabulated energies, from 0 at each pair's own threshold,
   and the 1s5 -> 2p9 shape (tabulated to 300 eV) above the tables. The
   electron loses the ledger's gap, 1.385 eV, as in v0_4. Under v0_2's
   closure the row is the Maxwellian table of that cross section and its
   reverse the Klein-Rosseland table (detailed balance over the gap to
   2e-4); under the kinetic EEDF the population process carries the
   measured cross section itself (the solver's rates match the tables to
   1.5e-3). Kannari's Table I [76] fit (Jacob and Mangano) stays as
   `Params.k76 = "kannari"`, which reproduces v0_4 to four digits.
2. **The band k_5s5p is now a scale on the measurement**, 1 by default,
   0.5 to 1.4: 1.4 the measurement's upper error (+35 to 40 %), 0.5 the
   metastables alone, the J = 1 levels contributing nothing (the one
   assumption, bracketed). `--ladder kannari` no longer touches it.
3. **The finding.** The measured set confirms the 1976 fit: its Maxwellian
   rate is 1.37 times Kannari's at 1.16 eV (Tahoma's electrons), 1.11 at 3
   eV and 0.91 at 5.8 eV (the reference point's). The lumped cross section
   is 1.14e-14 cm2 at 6 eV against the fit's effective 1.34e-14. So the
   partition v0_4 found rests on a measurement now, and moves little.
4. **A TALIF lifetime as a cross-check.** The TALIF review Gazeli et al.,
   Plasma 4, 145 (2021), Table 3, gives the natural lifetime of the Kr
   5p'[3/2]_2 calibration level as 34.1 ns (quoting Niemi et al.); PRAD's
   2J+1 mean of the ten measured level lifetimes is 29.0 ns (Dzierzega et
   al. 2000 give that level 28.6 ns). No change; the review tabulates no
   Kr self-quenching, so KPS stays on Gazeli et al., Phys. Plasmas 28,
   043301 (2021).
5. **A second laser, held out.** Kleikamp, Peters and Witteman, IEEE J.
   Quantum Electron. 22, 219 (1986), ran a coaxial e-beam KrF laser with
   krypton as the only buffer: with F2 at 10 mbar the output peaks at
   about 2 atm Kr, with Kr at 2 atm at about 10 mbar F2, and their
   steady-state cavity balance deduces a three-body KrF* + Kr quenching
   of 4.6e-31 cm6/s. `docs/krf_data/kleikamp_v0_5.py` scores the optima's
   locations (below). It is a score, never a target.

## At the reference point

1 atm Kr, 0.5 % F2, 300 kW/cm3 for 1 us, run to 2 us, percent of E_dep
unless marked; every band at its default.

| | v0_5 (two_t) | v0_4 (two_t) | v0_5 kinetic | v0_4 kinetic |
|---|---|---|---|---|
| 248 nm, KrF B-X | 5.68 | 5.79 | 5.73 | 5.68 |
| 147 nm, Kr2* | 4.56 | 4.53 | 4.82 | 5.3 |
| T_g at the end, K | 742 | 742 | 738 | 735 |
| T_e, pulse mean, eV | 5.68 | 5.81 | 2.96 | 2.97 |
| Kr* pumped to the rung, net, % of Kr* formed | | 81.0 | 79.4 | 78.4 |
| g0 / alpha, 500-1000 ns, %/cm | 4.12 / 0.63 | 4.12 / 0.63 | 3.79 / 0.67 | 3.73 / 0.68 |
| k_5s5p 0.5 / 1.4, 248 nm | 5.42 / 5.75 | (0: 4.67) | | |

`Params.k76 = "kannari"` gives v0_4's 5.787 / 4.53 / 741.7 exactly. Under
v0_2's closure the electrons sit at 5.7 eV, where the measured cycle is
9 % weaker than the fit's, so the 248-nm share falls by 0.1 point and the
Kr2* route gains it; under the kinetic EEDF they sit at 3 eV, where it is
10 % stronger, and the share rises by 0.05 point. The kinetic run: 425
solutions, the field carrying 1.001 to 1.003 of the input through the
pulse's second half, the two energy balances agreeing to 7e-6, 17 min.

## At Tahoma

Under the prototype beam drive, ML's configuration, the 200-400 ns window
(`tahoma_v0_4.py`, which scores any version; measured 6.39 / 1.11 / 5.28
%/cm, Spectra Technology 1989 Table 3-3):

| | g0 | alpha0 | g0 - alpha0 | T_e, eV | n_e, cm-3 |
|---|---|---|---|---|---|
| v0_5, two_t | 5.93 (x0.93) | 0.82 (x0.74) | 5.11 (x0.97) | 1.10 | 4.2e14 |
| v0_4, two_t | 5.86 (x0.92) | 0.82 (x0.74) | 5.04 (x0.95) | 1.18 | 4.5e14 |
| v0_5, kinetic | 5.99 (x0.94) | 0.83 (x0.75) | 5.16 (x0.98) | 1.08 | 4.1e14 |
| v0_4, kinetic | 5.91 (x0.93) | 0.82 (x0.74) | 5.09 (x0.96) | 1.16 | 4.4e14 |

At Tahoma's 1.1-eV electrons the measured cycle is 1.4 times the fit's,
cools the electrons a little more and lifts the gain by 1 % under either
closure (the kinetic EEDF and v0_2's closure agree to 1 % here, as in
v0_4). The absorption does not move: its 73 % is Kr2F*. The errors
against Tahoma, read as v0_4's overview reads them: gain -6 %, absorption
-25 %, net gain -2 %.

## The Kr/F2 laser of Kleikamp et al. (held out)

LIGHT under the beam drive at their conditions (250 A/cm2, 30 ns,
electrons of at most 170 keV, 1-cm tube; the deposited power 4.7 MW/cm3
at 2 atm proportional to the density, from Berger and Seltzer's stopping
power 2.2 MeV cm2/g and from their own B_Kr with W = 24.1 eV, 3.8 and 5.6
MW/cm3, run at the mean and at 0.7 and 1.4 times it), two_t closure, the
pulse's second half (15-30 ns). Two proxies for the output: the net
small-signal gain, and the paper's own saturated balance (their Eqs.
(1)-(3), the inversion clamped at threshold (gamma + alpha)/sigma_s with
gamma = 5.75e-2 1/cm their outcoupling), output ~ gamma/(gamma + alpha)
(P - Q [KrF]_thr) on LIGHT's formation rate P, loss rate Q and absorption.

| | LIGHT v0_5 | measured |
|---|---|---|
| Kr pressure of the output maximum, F2 10 mbar | 1.26 atm (saturated balance; 1.02 at 0.7x power, 1.48 at 1.4x); 0.67 atm (small-signal net gain) | ~2 atm |
| F2 pressure of the output maximum, Kr 2 atm | above 30 mbar (both proxies still rising) | ~10 mbar |
| KrF(B) loss rate Q at 2 atm, 10 mbar, in the window | 3.3e9 1/s (the trimer 2.8e9 of it) | |
| three-body KrF* + Kr quenching at the window's T_g (362 K) | 1.02e-30 cm6/s (G5 1.9e-30 at 300 K on the (T_g/300)^-3.32 law) | their deduced 4.6e-31 |
| the same with G5, G5C scaled to their value (x0.45, a sensitivity run, not adopted) | the saturated optimum at 1.83 atm (small-signal 0.99), the effective quenching 4.57e-31 | ~2 atm |

What it says: at the laser's 2 atm LIGHT quenches the threshold inversion
about twice as fast as the laser can bear, and puts the pressure optimum
at 1.0 to 1.5 atm where the laser found 2. The term is the trimer's rate
at the laser's gas temperature, 350 to 390 K here, which the propagated
uncertainty of v0_4 already named as the largest; the laser's deduced
constant is a model-dependent number (their Q is three-body only, their
output insensitive to the outcoupling), so it bounds rather than fixes
it. The F2 optimum is not reproduced: LIGHT's absorption falls with more
F2 (Kr2F* is quenched away faster than F2 and F- absorb), so neither
proxy turns over by 30 mbar where the laser peaked at 10. That points at
the F2-side absorbers (F2's 1.2e-20 cm2 is a compilation value; F-'s
5.6e-18 a 1971 measurement) and at KrF* + F2 quenching, none of them
moved here.

## Found on the way

- **The band's old meaning did not survive.** k_5s5p was Kannari's
  pre-exponential in cm3/s; as a scale on a measured cross section it is
  a factor, and `--ladder kannari` must not set it. Tests that set it to
  9e-7 meant 1.0.
- **A proxy for laser output must saturate.** The small-signal net gain
  peaks at 0.67 atm because it quenches the whole unsaturated inversion;
  a laser clamps the inversion at threshold and quenches only that, so
  its optimum sits at twice the pressure. The paper's own balance is the
  honest proxy, and even it falls short of 2 atm.
- **Rates are events per volume.** The first Kleikamp run summed rate
  constants; the formation and loss rates need the reactants' densities
  (`Compiled.heavy`, `ne_pow`).

## Limits

Section 7 of v0_4's overview, revisited:

- The 5s <-> 5p rate: measured now (item 1 above); the J = 1 proxies are
  the remaining assumption, bracketed by the band.
- KPS one level's total loss, K27 a search snippet: unchanged (no new
  source).
- The stepwise ionizations and the 5p harpoon as LIGHT-side rates: no new
  cross sections in the articles; unchanged, 0 by default.
- The absorbers as a diagnostic, the wall closure, no propagation: unchanged.
- New: the Kr/F2 laser's optima are not reproduced (above); the trimer's
  rate at 350-400 K and the F2-side absorbers are the terms to measure.

## Checks

`tests/test_v0_5.py` (the tables as printed, the lumped cross section's
sums, the Maxwellian table against a direct quadrature, detailed balance,
the kinetic population process against the tables, the measured rate
against the fit, the fit switch reproducing v0_4, the reference point
pinned, the band's ends, the ledger under the kinetic EEDF) and v0_4's
tests re-pointed with their K76 expectations rewritten. On 2026-10-10 the
fast files give 326 passed, 11 skipped (no v0_5 overview yet) and 1
strict expected failure (Hakuta's decays); test_kinetic.py 24 passed and
1 skipped. Details in the CHANGELOG entry of 2026-10-09 (c).

## Documents

`docs/accuracy_overview_v0_5.pdf` (from `accuracy_overview_v0_5.tex`;
numbers by `make_overview_v0_5.py` in the `\wg` namespace, v0_4's
generated numbers beside them as `\wf`, the map's as `\wm`), with the
errors table, the propagated uncertainty, the measured cross section's
checks, the Kleikamp laser and a reference list; `docs/system_map_v0_5.pdf`
(`make_system_map_v0_5.py`), the energy flow with the 5s <-> 5p cycle on
the measured cross sections; `docs/energy_ladder_v0_5.tex`.
`tests/test_docs.py` holds them to the model (keys defined, citations,
words against numbers, the ladder current; `LIGHT_DOCS_FAST=1` and
`LIGHT_DOCS_SLOW=1` regenerate and compare).

## Files

```
v0_5/
  README.md, SOURCES.md    this file; the sources (section 8 is v0_5's)
  jung2006.py    Jung et al. 2006 Tables I-IV, the lumped Kr(5s) -> Kr(5p)
                 cross section, its Maxwellian, superelastic and captured-
                 energy tables, Kannari's fit for comparison
  network.py     K76/K76r on ("pop_xsec", "jung2006"), the band k_5s5p as a
                 scale (Row.frac), kannari_k76_rows, build(k76=...)
  model.py       + Params.k76; the pop_xsec rate kind under both closures;
                 the kinetic population process from the measured cross
                 section (or the fit's pseudo cross section)
  run.py         + --k76
  tests/         test_v0_5.py and v0_4's, re-pointed
docs/krf_data/kleikamp_v0_5.py   the Kr/F2 laser's optima, two proxies
docs/accuracy_overview_v0_5.tex/.pdf, make_overview_v0_5.py,
     system_map_v0_5.tex/.pdf, make_system_map_v0_5.py,
     energy_ladder_v0_5.tex   the documents above
```
