# LIGHT v0_1

v0_1 is v0's working copy. It was copied from v0 on 2026-09-28 to work on
the KrF trimer, ion recombination, fluorine recombination, and the KrF C-A
band and C-state decay. v0 and the accuracy overview stay as
they were. The literature behind every change is in `SOURCES.md`.

```bash
python3 LIGHT/run_v0_1.py                        # the reference point
python3 LIGHT/run_v0_1.py --scan                 # the notes' pressure x power grid
python3 LIGHT/run_v0_1.py --trimer-tg-pow -0.5   # the trimer's temperature band
python3 LIGHT/run_v0_1.py --network v0           # the notes' list as written
python3 -m pytest LIGHT/v0_1/tests -q            # 77 tests, about 2 minutes
open LIGHT/docs/accuracy_overview_v0_1.pdf       # how accurate v0_1 is
open LIGHT/docs/REVIEW/v0_vs_v0_1.pdf            # v0 or v0_1, the same tests
```

## The input

The only energy input is a given power, deposited uniformly and directly
into the electron gas for the pulse. No drive is modeled, and v0's EM-wave
option is gone. No electrons are injected either. n_e starts from a seed
(`--ne`, 1e15 cm^-3 by default) and then follows ionization, attachment and
recombination, with the plasma quasineutral. The seed only sets the first
nanoseconds. A 1e12 seed gives the same 248-nm share as 1e15 to within
0.01 of a percent. Holding n_e fixed is still possible (`--ne-mode
clamped`), but with attachment and recombination in the network it makes
charge from nothing, and the report says so.

## What changed

**Ion recombination.** The charge chemistry is now part of the basic network
rather than optional groups:
- Kr2+ + e -> Kr* + Kr at Shiu and Biondi's measured 1.6e-6 (300 K/Te)^0.55
  cm3/s. Mikus measured the same to within 8 % up to 1.8 eV. The captured
  electron takes 0.95 Te from the electron gas.
- Kr+ + 2Kr -> Kr2+ + Kr at 2.4e-31 cm6/s, the mean of three measurements,
  with the T^-3/4 law of Perez-Rios and Greene's classical theory.
- F2 attachment, dissociation and vibration from the SIGLO cross sections.
- Kr+ and Kr2+ + F- -> KrF(B) on Kannari's density profile, which is a
  harmonic mean of a three-body part and a mobility part. Kannari used it
  gas-temperature dependent, and REBUILD had dropped that. v0_1 puts Thomson
  theory's (Tg/300)^-3 back on the three-body part. No measurement with F-
  exists.

The other F- channels are closed or negligible. Associative detachment
F- + F is endothermic by 1.8 eV, and electron detachment needs more than
13 eV.

**Fluorine recombination.** F + F + Kr -> F2 + Kr at 2.8e-34 cm6/s, Cohen
and Bott's evaluation, which gives every rare gas the Ar value. It is held
at its 300 K value, which overstates it at higher Tg. No Kr measurement
exists, and the Ar and He measurements are 8.0e-35 and 6.0e-34. At 1 atm
the F atoms recombine in milliseconds. At the reference point this re-forms
0.5 % of the F2 fill during the run, and the F2 still burns out.

**The trimer.** The rate and its (Tg/300)^-3.32 law stay, but the law's
stated range is now Shui's calculation, 300-800 K. It is measured only to
425 K, and nothing exists above that. Kr2F* now sits 0.58 eV below KrF(B)
+ Kr, the calculated binding (Huestis et al. 1984). REBUILD's estimate put
it 2.05 eV below, which moved 1.5 eV per trimer from its emission into its
formation. The trimer band is at the printed 420 nm.

**The C-A band and the C-state decay.**
- The C-A band sits at the printed 275 nm.
- B -> A emission at 272 nm, calculated at 3.5 % of B -> X, lies under the
  C-A band, and v0_1 carries it (row G4A). B's total radiative lifetime is
  unchanged.
- C keeps the joint fit's 50 ns lifetime, and the one measurement,
  Hakuta's 60 +- 10 ns, lies inside the fit's 41-63 ns.
- The B-C gap is Yu and Setser's measured 120 cm^-1, where v0 had 80.
- KrF(B) now lies at its measured well depth, 5.54 eV below the ion pair
  (Tellinghuisen 1976), which puts it 5.856 eV above ground Kr + F2. v0
  booked REBUILD's estimated 6.15 eV. With the measured energy a 248-nm
  photon leaves 0.06 eV of heat and a C-A photon 0.53 eV, as the
  spectroscopy implies.

## At the reference point

1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us. Each row changes one
thing from v0_1 as built (248-nm share of the deposited energy):

| Model | 248 nm | all bands | Tg end |
|---|---|---|---|
| v0 basic network, n_e held at 1e15 (the overview's model) | 5.24 % | 10.1 % | 604 K |
| v0's values with the ion and F2-electron groups, n_e free | 4.85 % | 10.2 % | 762 K |
| **v0_1** | **4.73 %** | **10.2 %** | **762 K** |
| v0_1 without the B -> A branch | 4.90 % | 10.2 % | 762 K |
| v0_1 without F + F + Kr | 4.69 % | 10.1 % | 762 K |
| v0_1 with v0's Kr+ + 2Kr and Kr2+ + e values | 4.75 % | 10.2 % | 762 K |
| v0_1 without the ion-ion (Tg/300)^-3 | 4.73 % | 10.2 % | 763 K |
| v0_1 with REBUILD's KrF and Kr2F* energies | 4.73 % | 10.2 % | 763 K |

The change from v0 comes mostly from letting n_e find its own value. With no
injection the plasma runs near 1e14 cm^-3 while F2 lasts, and 3.7e15 at its
peak after the F2 burns out at 613 ns. Electrons then quench 7 % of the KrF
instead of 19 %. More Kr* ends in Penning ionization, and recombination
turns that ionization back into heat, so the gas ends hotter. The B -> A
branch takes 3.4 % of the B-X light. The recombination values, the ion-ion
temperature law and the measured energies each move the 248-nm share by
0.02 or less here. The energies change where heat is booked, not how much
light leaves.

## Across the notes' ranges, and the trimer band

248-nm share of the deposited energy. The last columns replace the trimer's
temperature law by the pure-Kr T^-1/2, the weakest the data allow.

| p (atm) | P (kW/cm^3) | v0 basic, n_e held | v0_1 | Tg end (K) | v0_1, trimer T^-1/2 |
|---|---|---|---|---|---|
| 0.1 | 100 | 3.73 | 3.46 | 1176 | 3.39 |
| 0.1 | 300 | 1.25 | 0.80 | 2437 | 0.79 |
| 0.1 | 500 | 0.75 | 0.39 | 3461 | 0.39 |
| 0.25 | 100 | 8.43 | 9.22 | 852 | 8.06 |
| 0.25 | 300 | 2.83 | 2.39 | 2083 | 2.15 |
| 0.25 | 500 | 1.71 | 1.21 | 3141 | 1.10 |
| 0.5 | 100 | 11.31 | 14.21 | 546 | 9.54 |
| 0.5 | 300 | 4.38 | 4.11 | 1313 | 2.93 |
| 0.5 | 500 | 2.63 | 2.24 | 2081 | 1.64 |
| 1 | 100 | 4.79 | 5.73 | 441 | 3.43 |
| 1 | 300 | 5.24 | 4.73 | 762 | 2.28 |
| 1 | 500 | 3.09 | 2.72 | 1154 | 1.35 |
| 2 | 100 | 1.00 | 1.25 | 379 | 0.88 |
| 2 | 300 | 2.30 | 2.45 | 524 | 0.97 |
| 2 | 500 | 2.81 | 2.31 | 688 | 0.78 |

- The trimer's temperature law is now the largest uncertainty in the model.
  At the reference point the two laws give 4.73 % and 2.28 %, and at 2 atm
  they differ by up to three times. Shui's calculation is the default. The
  T^-1/2 came from fitting a laser model to laser efficiencies, so under the
  predictive-only policy it only bounds the law.
- With recombination the ionization energy comes back as heat. Below 0.5
  atm and at 300 kW/cm^3 or more the gas ends above 1000 K, and at 0.1 atm
  it passes 3000 K. There Kr2+ and Kr2F* would redissociate, which the model
  does not include, and every rate law is extrapolated. The report warns.
- The best 248-nm share in the grid is 14.2 %, at 0.5 atm and 100
  kW/cm^3.

## What stays open

1. The trimer's temperature dependence above 425 K, and Kr2F*
   redissociation above about 600 K. Neither has a measurement.
2. Ion-ion recombination with F-. Its rate was measured only below 300
   Torr. Jacob et al. (APL 32, 109, 1978) fitted KrF* pulse shapes to Kr+
   + F- (+M) = 2e-7 at 60 Torr and 3e-7 cm3/s at 120 Torr, where v0_1's
   form, screened at their 1e14 cm-3 ions, gives x1.2 and x1.6. Rokni et
   al. (PRA 16, 2216, 1977) found that Kr2+ + F- forms KrF* as efficiently
   as Kr+ + F-, from relative intensities. The absolute KrF* yield, its B/C
   split and the rate near 1 atm have no measurement.
3. Hakuta's C-A decays still conflict with the three datasets the B/C
   constants were fitted to.
4. Kr2+ redissociation above about 1000 K has no published rate.
5. Kr2* + F -> KrF(B) + Kr, Kannari's estimate of 3.0e-10 cm3/s. The
   optional Kr2* group raises the reference point's 248-nm share from 4.73
   to 9.17 %, and without this one row it gives 4.65 %. 80 % of the gain
   comes after F2 burns out, when the F atoms left behind turn Kr2* into
   KrF 17 times faster than Kr2* radiates. At Tahoma the row moves g0 from
   5.69 to 6.10 %/cm. The 3.0e-10 traces to SRI MP-76-99, which has no
   such rate. Its model leaves Rg2* + F out, and its one number is a
   conditional 3-10e-10 for Ar2* + F. Chen & Payne (1979) assign the slow
   KrF* component SRI saw to ion-ion recombination instead. 3.0e-10 is no
   upper bound, but the share saturates, 9.49 % at 1e-9 and already 5.85 %
   at 1e-11. Until the rate is measured the group is a band from 4.65 to
   9.49 %, like the trimer law. Kr2* + Kr2* Penning ionization rests on an
   unpublished 1973 report, and the even split of Kr2* + F2 is Kannari's.

## Accuracy

`docs/accuracy_overview_v0_1.pdf` scores v0_1 on every test v0's overview
used, with v0_1's map and every rate. `docs/REVIEW/v0_vs_v0_1.pdf` compares
the two versions. They score the same wherever they run the same
chemistry. Under Tahoma's beam v0_1's default network gives the gain at
x0.93 of the measured value, where v0's default gives x0.32, since v0's
holds n_e fixed and has no charge chemistry. The trimer's weak T^-1/2 law
gives Tahoma's gain at x0.47 (the window sits at 448 K), so Tahoma counts
against it, as a score.

## Tests

`python3 -m pytest LIGHT/v0_1/tests -q` runs 77 tests (74 pass, Hakuta's
C-A decays are a strict expected failure, and 2 docs checks run only with
`LIGHT_DOCS_SLOW=1`). `test_docs.py` holds the map, the overview and the
comparison to fresh runs, the ML-backed ones in a process of their own, the energy ladder to the species and network, and every printed number to a generated file. They cover v0's checks,
adapted, plus the new rates against the numbers their sources print, the
B-A branching, detailed balance over the measured gap, the Tg law of
ion-ion recombination, F + F + Kr's size, the trimer switch, the
charge-conservation warning and a reference-point regression. A mutation
check planted 12 bugs in the new code, one at a time, and every one failed
the suite.

Every comparison in these tests uses `approx(..., abs=0)`.
`pytest.approx` adds an absolute tolerance of 1e-12 by default, even when
a relative one is given, so any two rate constants below 1e-12 compare
equal. v0's tests use it as it comes. There, for example, the check that
G5 is 6.0e-31 passes against the network's 1.9e-30, and the REBUILD drift
gate cannot see a change in a three-body rate. A strict re-run found no
drift in v0's copied values. v0's step-cross-section check claims 0.2 %,
where the quadrature is good to half a grid cell of the Boltzmann factor
(0.7 % at 0.8 eV), and v0_1's check uses that bound.
