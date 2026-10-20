# LIGHT v0_2

v0_2 is v0_1's working copy. The owner copied it from v0_1 on 2026-09-29
to take the Kr2* route into the basic network and to carry every rate that
is known only as a published spread, or not at all, as a band. v0, v0_1
and their documents stay as they were. What v0_1 changed from v0 is in
`LIGHT/v0_1/README.md`, and the literature behind both versions is in
`SOURCES.md` (section 5 for v0_2).

```bash
python3 LIGHT/run_v0_2.py                        # the reference point
python3 LIGHT/run_v0_2.py --bands                # each band's two ends, about 45 s
python3 LIGHT/run_v0_2.py --rate kr2_f=3e-10     # one band at one value
python3 LIGHT/run_v0_2.py --scan                 # the notes' pressure x power grid
python3 LIGHT/run_v0_2.py --network v0           # the notes' list as written
python3 -m pytest LIGHT/v0_2/tests -q            # 82 tests, about 5 minutes
open LIGHT/docs/REVIEW/v0_vs_v0_1_vs_v0_2.pdf    # v0, v0_1 or v0_2, the same tests
open LIGHT/docs/accuracy_overview_v0_2.pdf       # how accurate v0_2 is, with its map
```

The input is v0_1's. A given power goes uniformly and directly into the
electron gas, no electrons are injected, and n_e follows ionization,
attachment and recombination from its seed.

## What changed

**The Kr2* route is basic.** Kr* + 2Kr -> Kr2* + Kr, Kr2* radiation at 147
nm, Kr2* + F2 in its two channels, Kr2* + F and Kr2* + Kr2* are rows of the
basic network, where v0_1 had them as an optional group. The measured rows
keep their measurements, Ueda's 1.7e-32 cm6/s on Wieme and Lenaerts'
(Tg/300)^-0.64, the 270-ns lifetime and Quigley and Hughes' 2.1e-10
cm3/s split evenly. Kr2* + F and Kr2* + Kr2* have no measurement, so they
are 0 unless a band run sets them. With Kannari's estimates for those two,
3.0e-10 and 5.0e-10, v0_2 gives v0_1's optional group to one part in
1e8. `--add kr2` is refused on the basic network and still works on v0's.

**The bands.** A rate carried as a band keeps its measurement as the
default and spans the other published values. A rate with no measurement
is 0 by default and runs to the largest published value, which is no
bound. Analogs and device data never set a default or narrow a band. The
bands live in `network.BANDS`, and rows read them through the new rate kind
`("param", name[, frac])`.

| Band | Process | Default | Low | High |
|---|---|---|---|---|
| kr2_f | Kr2* + F -> KrF(B) + Kr | 0 | 0 | 1e-9 (SRI's conditional 3-10e-10 for Ar2*) |
| kr2_pen | Kr2* + Kr2* -> Kr2+ + 2Kr + e | 0 | 0 | 5e-10 (Kannari, from an unpublished report) |
| kr2f_f | Kr2F* + F -> 2Kr + 2F | 0 | 0 | 1.8e-10 (the Ar2F* analog) |
| kr2_form | Kr* + 2Kr -> Kr2* + Kr, 300 K | 1.7e-32 (Ueda) | 1.7e-32 | 5.5e-32 (Hughes) |
| kr2_f2 | Kr2* + F2, both channels | 2.1e-10 (Quigley & Hughes) | 2.1e-10 | 5.0e-10 (Chen & Payne) |
| kr_pen | Kr* + Kr* -> Kr+ + Kr + e | 5.0e-10 (v0's) | 1.5e-10 (Werner) | 5.0e-10 |
| kr2f_f2 | Kr2F* + F2 | 1.2e-10 (Xu) | 1.2e-10 | 4.3e-10 (Quigley & Hughes) |
| trimer_tg_pow | the trimer's (Tg/300)^n | -3.32 (Shui) | -3.32 | -0.5 (pure-Kr efficiencies) |

**Kr2F* + F is new.** Jancaitis found in 1983 that after the pulse, with
the F2 nearly all dissociated, Ar2F* was quenched as fast as the initial F2
would quench it, so F atoms quench the trimer at 0.5 to 1 times the F2
rate. That was argon (his Fig. 11). Carried to Kr2F* it is an analog, so it
sets the top of the band and the default stays 0.

**The warnings and the report.** The report warns when the gas passes 451
K, the top of Kr* + 2Kr's measured temperature law. Its header names the
bands set and the unmeasured ones at 0. `row_dict` gains the 147-nm share.

**Citations.** Kr2+'s 1.15 eV binding is a photoionization measurement (Ng
et al. 1977), where the species table called it an estimate. Kr2* keeps
REBUILD's estimated well, since no measured or calculated binding is among
the sources held. It only moves energy between heat and storage. The
ion-ion citations and `SOURCES.md` carry the owner's 1977-78 papers. Jacob
et al. measured Kr+ + F- below 300 Torr, so the rate near 1 atm is still
theory.

## At the reference point

1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us.

| Model | 248 nm | 147 nm | all bands | Tg end |
|---|---|---|---|---|
| v0_1 | 4.73 % | none | 10.2 % | 762 K |
| **v0_2** | **4.67 %** | **20.0 %** | **30.4 %** | **649 K** |
| v0_2 with Kannari's Kr2* + Kr2*, 5.0e-10 | 4.65 % | 10.8 % | 21.1 % | 704 K |
| v0_2 with both Kannari estimates (v0_1's kr2 group) | 9.17 % | 0.5 % | 18.3 % | 724 K |

The measured part of the route barely moves the laser line. It takes 9 %
of Kr* before F2 burns out at 640 ns, and after that the Kr* that v0_1
turned into Penning ionization and heat radiates at 147 nm instead. The
gas ends 113 K cooler.

`--bands` at the reference point, each band alone at its two ends:

| Band | End | 248 nm | change | 147 nm | Tg end |
|---|---|---|---|---|---|
| kr2_f | 1e-9 | 9.50 % | +4.83 | 0.5 % | 726 K |
| trimer_tg_pow | -0.5 | 2.22 % | -2.44 | 19.5 % | 658 K |
| kr_pen | 1.5e-10 | 4.87 % | +0.21 | 25.0 % | 616 K |
| kr2f_f2 | 4.3e-10 | 4.88 % | +0.21 | 19.9 % | 663 K |
| kr2f_f | 1.8e-10 | 4.87 % | +0.20 | 19.9 % | 670 K |
| kr2_form | 5.5e-32 | 4.52 % | -0.15 | 28.0 % | 601 K |
| kr2_f2 | 5.0e-10 | 4.60 % | -0.06 | 19.1 % | 654 K |
| kr2_pen | 5e-10 | 4.65 % | -0.01 | 10.8 % | 704 K |

- Two bands matter for the laser line, Kr2* + F (4.67 to 9.50 %) and the
  trimer law (4.67 to 2.22 %). Each of the others moves it by 0.21 points
  or less.
- The Kr2F* quench bands raise the 248-nm share through the gas
  temperature. Their energy goes to heat instead of 420-nm light, and the
  hotter gas slows the trimer.

## Across the notes' ranges

248-nm share of the deposited energy, and the same at the top of the Kr2*
+ F band.

| p (atm) | P (kW/cm^3) | v0_1 | v0_2 | Tg end (K) | F2 left (%) | 147 nm | v0_2, Kr2* + F at 1e-9 |
|---|---|---|---|---|---|---|---|
| 0.1 | 100 | 3.46 | 3.47 | 1173 | 0.0 | 0.8 | 3.55 |
| 0.1 | 300 | 0.80 | 0.80 | 2434 | 0.0 | 0.3 | 0.81 |
| 0.1 | 500 | 0.39 | 0.39 | 3458 | 0.0 | 0.2 | 0.40 |
| 0.25 | 100 | 9.22 | 9.23 | 826 | 0.0 | 4.9 | 10.33 |
| 0.25 | 300 | 2.39 | 2.40 | 2044 | 0.0 | 2.7 | 2.61 |
| 0.25 | 500 | 1.21 | 1.21 | 3102 | 0.0 | 1.8 | 1.30 |
| 0.5 | 100 | 14.21 | 13.84 | 532 | 2.2 | 5.7 | 16.13 |
| 0.5 | 300 | 4.11 | 4.13 | 1191 | 0.0 | 11.1 | 5.86 |
| 0.5 | 500 | 2.24 | 2.26 | 1928 | 0.0 | 8.3 | 3.03 |
| 1 | 100 | 5.73 | 5.42 | 440 | 44.4 | 1.5 | 5.84 |
| 1 | 300 | 4.73 | 4.67 | 649 | 0.2 | 20.0 | 9.50 |
| 1 | 500 | 2.72 | 2.71 | 915 | 0.2 | 24.7 | 7.18 |
| 2 | 100 | 1.25 | 1.17 | 378 | 72.7 | 1.0 | 1.28 |
| 2 | 300 | 2.45 | 2.17 | 519 | 22.2 | 2.5 | 2.64 |
| 2 | 500 | 2.31 | 2.12 | 604 | 0.9 | 17.3 | 4.48 |

- Kr2* + F matters where Kr2* forms fast and F2 burns out, at 1 atm and
  above with 300 kW/cm^3 or more. There its band spans up to a factor 2.7.
  At 0.1 atm three-body formation is too slow for it to matter.
- The measured route lowers the 248-nm share by up to 0.37 points from
  0.5 atm up, most at 0.5 atm and 100 kW/cm^3, and leaves 0.1 and 0.25 atm
  unchanged.

## Pure krypton (the owner's step 4)

ML's pure-Kr check runs LIGHT with no F2 under the prototype beam drive
(`ML/ml/pure.py kr_vuv`, `LIGHT_VERSION=v0_2`). At Morozov's pump, 330
W/cm^3 at 665 Torr, v0_2 gives a VUV yield of 0.571 against the measured
0.42 +- 0.05. Every Kr2* band moves it by 0.001 or less (Kr2* + Kr2* at
5e-10 gives 0.570, Kr* + 2Kr at 5.5e-32 gives 0.571), and without the
thermal-electron excitation it is 0.552, at the ideal of one photon per
excitation and ion. So the route cannot explain the miss, and Morozov's
point cannot score it. At 300 kW/cm^3 and the same 665 Torr the bands
move the pure-Kr yield from 0.13 (Kr2* + Kr2* at 5e-10) through 0.27 (the defaults) to 0.49
(Kr* + 2Kr at 5.5e-32). A pure-Kr VUV yield measured against pump rate,
or time-resolved Kr2* decays at high density, would score them. None is
among the sources held.

## What stays open

1. Kr2* + F. No measurement, and at 1-2 atm and 300-500 kW/cm^3 its band
   spans up to a factor 2.7 in the 248-nm share. A measurement, or a bound,
   would matter between 1e-12 and 1e-10 cm3/s.
2. The trimer's temperature law above 425 K, and Kr2F* redissociation.
3. Kr2* + Kr2* and Kr* + Kr* Penning ionization. Neither has a measurement
   among the sources held, and at high pump they set the pure-Kr VUV yield.
4. Kr* + 2Kr above 451 K. The reference point ends at 649 K and the
   report warns.
5. Ion-ion recombination with F- near 1 atm, its rate, KrF* yield and B/C
   split.
6. Hakuta's C-A decays, which still conflict with the fitted B/C
   constants, and Kr2+ redissociation above about 1000 K.
7. Kr2*'s binding energy.

## Accuracy

`docs/accuracy_overview_v0_2.pdf` gives v0_2's accuracy on every test,
with its system map (`docs/system_map_v0_2.pdf`), each band at both ends at
the reference point and at Tahoma, and its rates and bands in Appendix A.
`docs/REVIEW/v0_vs_v0_1_vs_v0_2.pdf` scores v0_2 beside v0 and v0_1 on the
same tests. v0_1 and v0_2 score alike on every kinetic dataset, on the
pure-gas checks and on the electron closure, since they share the KrF(B,C)
block and every formation link. Under Tahoma's beam v0_2's default network
gives the gain at x0.89 and the formation efficiency at x0.99. v0_1 as ML
scores it, with Kannari's two Kr2* estimates, gives x0.95 and x1.06, and
v0_2 with the same two estimates gives v0_1's numbers. Tahoma is device
data, so it scores those two rates and never sets them. The review's
bottom line is that v0_2 is the best of the three, since v0_1 comes closer
on Tahoma's gain only by leaving out the measured Kr* + 2Kr or by adding
the two estimates.

## Tests

`python3 -m pytest LIGHT/v0_2/tests -q` runs 82 tests (79 pass, Hakuta's
C-A decays are a strict expected failure, and the review's and the accuracy
overview's slow numbers are checked against a fresh run only with
`LIGHT_DOCS_SLOW=1`). Since 2026-09-30 `test_docs.py` also holds the
accuracy overview and the system map to the model: the map's numbers, the
rate tables and the fast numbers against fresh runs, every short source to
its row's citation, every printed number to a generated file, and the
energy ladder to v0_2's species and network. v0_1's
docs tests stay with v0_1. New since v0_1 are the Kr2* route in the basic network, the
bands held to their published values, the param rows, the band runner,
the 451-K warning and a reference-point regression that also pins v0_1's
kr2 group through Kannari's two estimates. The harpoon-seed test now
switches Kr* + 2Kr off to isolate the harpoon, and checks the route's
analytic share (9 % at 1 atm) with it on. A mutation check planted 13
bugs in the new code, one at a time, and every one failed the suite. Every
comparison uses `approx(..., abs=0)`, as in v0_1.
