# ML: fitting LIGHT's channels to the measurements

Owner direction (2026-09-27): a place to train the best models to the
available data. ML holds four things:
- the datasets;
- the forward models (LIGHT itself, its `v0` package or the version the
  environment names, `LIGHT_VERSION=v0_1` since 2026-09-29);
- the fitters: least squares, an ensemble MCMC and a neural-network
  estimator;
- a scorecard.

LIGHT stays the physics model, and ML fits its channels. Three models so far,
all in LIGHT since the owner's four steps of 2026-09-28:
- The KrF(B,C) block of LIGHT's basic network (10 channels), adopted as
  step 1.
- KrF formation: the chain from the beam's Kr+ and Kr* to KrF(B,C), each
  link against its measurements. The measured links went into LIGHT as
  step 2, and step 4 searched for the rest.
- The electron side: the LoKI-MC cross sections against swarm
  measurements, and the EEDF of the laser plasma. Its exact captured
  energies and kinetic tail are LIGHT's closure since step 3
  (`ml/tail.py`).

Since 2026-10-15 ML also trains LIGHT_Kr, the pure-krypton model: its
electron cross sections against krypton's swarm data (`ml/kr_*.py`,
"Pure krypton" below), pinned by LIGHT_Kr v0_3.

## The rule for data

The repo's predictive-only policy (`KrF_and_LoKI/agents/CLAUDE.md`,
`PLAN_physical_validity.txt`) allows a rate to come from "a published
measurement of that process". Anchors "score the model after the fact; they
never steer a number". `ml/data.py` gives every dataset one of three roles:

| Role | What it is | What it does |
|---|---|---|
| train | kinetic measurements of the channels being fitted | fits the channels |
| held_out | kinetic measurements left out | tests the fit |
| device | laser performance (gain, efficiency) | scores only, never trains |

- `tests/test_data.py` fails if a device dataset is given the train role,
  and `fit.fit` refuses one. Changing that is an owner decision, and the
  policy text would need changing first.
- The train set is the largest group of kinetic datasets that agree with one
  another within their errors (`python3 ML/run.py consistency`).
- Electron swarm measurements (drift velocity, ionization and attachment
  coefficients, attachment rate constants in a buffer gas) measure the
  cross sections themselves, so they may train them. `ml/swarm.py` gives
  them the roles train, held_out and check.

## Quick start

```
python3 ML/run.py data          # the registry: role, points, error, source
python3 ML/run.py consistency   # which datasets agree (leave-one-out), 2 s
python3 ML/run.py fit           # least squares + MCMC on the train set, 1.5 min
python3 ML/run.py nn            # the neural-network estimator, tested, 1 min
python3 ML/run.py score         # every candidate set against every dataset, 5 s
python3 ML/run.py electron      # cross sections vs swarm data, the EEDF, Tahoma, 5 min
python3 ML/run.py formation     # the KrF formation chain against Tahoma, 2 min
python3 ML/run.py pure          # the F2-free half of the chain on pure-gas data, 5 s
python3 ML/run.py all           # all of it, 10 min
python3 ML/run.py tail          # regenerate LIGHT's tail table, 4 min (not in all)
python3 ML/run.py kr-data       # krypton's swarm datasets and their roles
python3 ML/run.py kr-consistency  # which krypton datasets agree (joint fit)
python3 ML/run.py kr-fit        # LIGHT_Kr's factors: least squares, MCMC, network
python3 ML/run.py kr-train      # the same, and the factors file LIGHT_Kr pins
                                # (the first kr command simulates the solver
                                # for the emulator: 16 min on 6 workers)
python3 -m pytest ML/tests -q   # 45 tests, 80 s (/usr/bin/python3 has pytest)
```

ML runs v0 unless `LIGHT_VERSION` names another LIGHT version. The version
is fixed when `ml` is first imported, so a run for v0_1 needs its own
process (`LIGHT_VERSION=v0_1 python3 ...`). v0_1 carries the ion and
F2-electron rows in its basic network, shares B's radiative rate between
B-X and B-A, and puts a Tg law on ion-ion recombination, and ML follows
each (`tests/test_light_version.py`). The tail table is made from v0's
states, and `ml/tail.py` refuses to write it from another version.

Outputs also go to `ML/results/<command>.txt`, which is gitignored.

## What is here

```
ML/
  run.py          the command line
  ml/data.py      the dataset registry (values from LIGHT/docs/krf_data/)
  ml/params.py    the channels, their bounds and priors
  ml/forward.py   the forward models: LIGHT's closed-form B/C block
                  (LIGHT/v0/two_level.py) for kinetic data, and LIGHT under
                  the prototype beam drive for device data
  ml/fit.py       least squares, the consistency check, the MCMC
  ml/nn.py        neural posterior estimation (numpy MLP, Gaussian head)
  ml/score.py     the scorecard
  ml/swarm.py     electron swarm datasets and the cross-section fit
  ml/eedf.py      the laser plasma's EEDF and LIGHT with it
  ml/electron.py  changes to LIGHT's electron rows (scales, tables, capture)
  ml/formation.py the formation chain on LIGHT, and the probe-flux
                  observables g0, Isat and g0 Isat / Pin
  ml/formation_data.py  the chain's measurements, link by link
  ml/pure.py      verification on pure rare gases (the F2-free half)
  ml/pure_data.py pure-gas measurements: W, N_ex/N_i, VUV yields, Kr ions
  ml/tail.py      LIGHT's two-temperature tail table from kinetic EEDFs
  ml/boltzmann/   two-term Boltzmann solver (vendored, md5 in __init__.py)
                  and ee.py (electron-electron collisions, beam drive)
  ml/kr_light.py  LIGHT_Kr's version under its own name (LIGHT_KR_VERSION)
  ml/kr_data.py   krypton's swarm datasets, roles and DOIs
  ml/kr_swarm.py  LIGHT_Kr's factors (FREE: the seven that train): forward
                  model, error model, fits, consistency, MCMC, the record
                  LIGHT_Kr pins
  ml/kr_nn.py     the swarm emulator and neural posterior estimation
  ml/kr_train.py  the krypton training end to end, and the comparison of
                  which factors train (free_variants)
  data/swarm/     the swarm measurements, provenance in each file's header
                  (Kr_*_LAPLACE_LXCat.txt: LXCat's output, md5-pinned)
  tests/          registry, forward models, fitter, network, electron side
```

## The datasets (2026-09-27)

| Name | Role | What |
|---|---|---|
| eden_kr, eden_f2 | train | Eden 1978: cold KrF(B) alone, decay against Kr and F2 |
| ksu_fig1, ksu_fig2 | train | Gadomski 1992: coupled B/C decays against F2 and Kr |
| ys_ir300 | train | Yu & Setser 1990: equilibrium C/B ratio at 300 K |
| xu_figs5_7 | held_out | Xu 1993: pure Kr/F2 decays at 2 Torr F2 (same lab as Gadomski) |
| hakuta_fig2, hakuta_fig3 | held_out | Hakuta 1990: C-A decays. Inconsistent with the train set |
| ys_ratio_230_300 | held_out | Yu & Setser: the 230/300 K ratio (tests the B-C gap) |
| bs_tau_b | held_out | Burnham & Searles 1977: tau_B = 9.0 ± 0.6 ns |
| tahoma_g0, tahoma_formation, tahoma_isat | device | the Tahoma laser, 99.6 % Kr/F2: g0, g0·Isat/Pin, Isat |

The sources, page references and digitizations are in `LIGHT/docs/krf_data/`.

## Results for the KrF(B,C) block

**Consistency.** chi2 per point, where about 1 means a dataset agrees with
its own stated error:

| Left out | Eden | Kansas State | Hakuta | Yu & Setser |
|---|---|---|---|---|
| none | 0.25 | 2.00 | 1.97 | 2.45 |
| Eden | 20.1 | 1.89 | 0.91 | 0.53 |
| Kansas State | 0.17 | 18.5 | 0.50 | 0.00 |
| **Hakuta** | **0.07** | **0.20** | 11.9 | **0.00** |
| Yu & Setser | 0.34 | 1.65 | 1.92 | 99.4 |

Only Eden, Kansas State and Yu & Setser agree, so they form the train set.

**The fit** (MCMC on the train set: median and 68 % interval):
- k_BC 9.8e-12 [4.6e-12, 1.7e-11] cm³/s. B-C mixing is slow.
- K_eq 1.07 [0.84, 1.38].
- tau_C 50 [41, 63] ns.
- F2 quench: B 4.6e-10, C 3.0e-10 cm³/s.
- Not pinned by these data: the Kr two-body quench (wide posteriors) and
  the B-state three-body rate.
- **Lasing B share at 1 atm: 0.69 [0.57, 0.80].** LIGHT gave 0.27 until
  2026-09-27 and gives 0.68 since adopting the fit (2026-09-28, step 1).

**The neural network** (30 000 simulated experiments):
- Network 1 (train set). On held-out simulations: R² 0.94 for k_BC, 0.97
  for K_eq, 0.86 for tau_C and 0.90-0.93 for F2, with ±1σ coverage
  0.64-0.65 (ideal 0.68). On the real data it gives the MCMC's answer:
  share 0.73 [0.61, 0.85].
- Network 2 (train set plus Hakuta, which conflicts). It returns a confident
  compromise, share 0.56 ± 0.04, that misfits Kansas State (chi2 15 per
  point). Nothing in its output flags this. The consistency check, not the
  network, is what separates conflicting data.

**Scorecard** (`python3 ML/run.py score`, 2026-09-28). rms or signed %
against each dataset. The Tahoma rows give the value and its ratio to the
measurement, on LIGHT's network and electron closure after the four steps
and the measured harpoon yield (0.90).

| Dataset (role) | LIGHT now (the fit) | LIGHT to 09-27 | Gadomski 1992 | Fit (train set) |
|---|---|---|---|---|
| eden_kr, eden_f2 (train) | 2.3, 3.2 % | 26.1, 10.5 % | 25.4, 10.5 % | 2.0, 3.2 % |
| ksu_fig1, ksu_fig2 (train) | 6.1, 4.2 % | 4.0, 7.5 % | 7.4, 3.5 % | 6.0, 4.3 % |
| ys_ir300 (train) | −0.5 % | +277 % | +202 % | −0.1 % |
| xu_figs5_7 (held out) | 11.7 % | 10.5 % | 12.7 % | 11.6 % |
| hakuta_fig3, hakuta_fig2 (held out) | 34, 19 % | 32, 19 % | 39, 26 % | 34, 19 % |
| bs_tau_b (held out) | −22 % | −22 % | −22 % | −22 % |
| ys_ratio_230_300 (held out) | −14 % | −14 % | −14 % | −14 % |
| **tahoma_g0 (device)** | **6.11 (×0.96)** | 3.17 (×0.50) | 3.28 (×0.51) | 6.13 (×0.96) |
| tahoma_formation (device) | **0.284 (×1.05)** | 0.291 (×1.08) | 0.289 (×1.07) | 0.284 (×1.05) |
| tahoma_isat (device) | **1.54 (×1.10)** | 3.03 (×2.16) | 2.91 (×2.08) | 1.53 (×1.09) |

The formation and Isat rows use the model's own g0·Isat/Pin and Isat (the
formation section below). Isat was in no fit.

Over 40 posterior draws the fit predicts Tahoma's gain as 5.83 %/cm
[5.42, 6.63] on today's network and closure (7.19 [6.39, 7.98] on
2026-09-27's). The measured value is 6.39, and no fit used it.

## Results for the electron side (2026-09-27)

Owner question: does training LoKI-MC's cross sections and EEDF improve
the model? `python3 ML/run.py electron` prints all of this.

**The solver.** A two-term Boltzmann solver (Hagelaar & Pitchford), vendored
from the website's LXCat sim. For He + 1 % F2 it matches LoKI-MC to 1 % at
6-15 Td. It fails at 3 Td, and for Kr/F2 at 0.3 Td, where attachment
outruns energy relaxation. Those points are not modelled, and LoKI-MC's
0.3 Td value stands.

**Cross sections against swarm data** (chi2 per point, about 1 = agrees):

| Dataset (role) | SIGLO or Biagi | TRINITI | Model/data |
|---|---|---|---|
| Kr W, Nakamura 1991 (check) | 0.12 | | 0.96-1.03 |
| Kr W, Pack 1962 (check) | 4.7 | | 0.80-1.19 |
| Kr alpha/N, Heylen 1971 (check) | 6.8 | | 0.54-1.05 |
| Kr alpha/N, Kruithof 1940 (check) | 96 | | 0.02-1.50 |
| F2 k_att in He, Nygaard 1978 (train) | 0.98 | 5.3 | 0.93-1.51 |
| F2 k_att in N2 300 K, Trainor & Jacob 1979 (train) | 0.44 | 0.07 | 0.69-0.88 |
| F2 k_att in N2 500 K, Trainor & Jacob (train) | 2.0 | 0.67 | 0.56-0.73 |
| F2 k_att in N2, Chen 1977 (train) | 3.0 | 4.5 | 0.94-2.4 |
| F2 k_att, Schneider & Brau 1978 (held out) | 47 | 25 | 0.48-0.53 |
| F2 thermal, Sides 1976 and McCorkle 1986 (held out) | 12 | 5.0 | 0.98-5.1 |

- Kr: Biagi's drift velocity matches Nakamura. Its ionization coefficient
  is low near threshold, where the two measurements differ by x2 from each
  other. Neither reaches the laser, whose thermal tail is empty (below) and
  whose beam yields come from Spencer-Fano.
- F2: near 1 eV the measured rate constants spread x3.5 between groups
  (Chen 2.1e-9, Trainor & Jacob 4.5e-9, Schneider & Brau 7.5e-9 cm3/s).
  SIGLO sits in the middle.

**Training the F2 attachment.** One scale on the main resonance (0-2.5 eV).
The 6 eV resonance carries under 10 % of any dataset's rate and 0.3 % of
the laser's, so no dataset pins it.
- Leave-one-out: Nygaard and Trainor & Jacob at 300 K agree. Chen conflicts
  with them. The 500 K points rise 25 % above the 300 K ones at the same
  mean energy, a gas-temperature effect no cross-section set here carries.
- Train set Nygaard + Trainor & Jacob 300 K: **s_low = 0.921 +- 0.045**.
  SIGLO as it is lies 1.8 sigma away.

**The EEDF of the laser plasma.** LoKI-MC and the two-term solver describe a
swarm, with electrons too few to collide with one another. At Tahoma the
ionization degree is 2e-5. `ml/boltzmann/ee.py` adds the e-e collision
term (checked against a Maxwellian's invariance, the Spitzer slowing-down
rate and energy conservation) and a beam drive: electrons appear with the
sub-excitation spectrum, 4.5 eV each, and leave by attachment and
recombination. At LIGHT's Tahoma state (fit set, 200-400 ns):
- The bulk is Maxwellian. Every rate below 9 eV is within 4 % of the
  Maxwellian with the same mean energy.
- The tail above the Kr thresholds is emptied by excitation faster than
  e-e collisions refill it. Kr(5s) x0.07, 5p x0.005, 4d+6s x6e-4, 5d
  x1e-5, ionization x2e-7.
- An attaching electron carries off 0.54 eV. LIGHT books (3/2)Te = 1.69 eV,
  which overstates the electrons' cooling by 21 % of their heating.
- Without e-e (a swarm code's physics) attachment falls to x0.29 of the
  Maxwellian value and n_e nearly triples. LoKI-MC's EEDF is the wrong
  EEDF for the laser.

`ml/eedf.py` scales LIGHT's electron rows by the ratio to the same-mean
Maxwellian, books the exact capture energies, and iterates LIGHT with the
Boltzmann solution to a fixed point. There LIGHT's Te and n_e agree with
the kinetic solution to 0.4 % and 5 %.

**Tahoma** on 2026-09-27's network, a score only (g0 measured 6.39 %/cm,
formation 0.27):

| B/C set | Electron side | g0 | Formation |
|---|---|---|---|
| fit (train set) | LIGHT today | 7.28 (x1.14) | 0.341 (x1.26) |
| fit (train set) | Beam EEDF + exact capture | **6.78 (x1.06)** | **0.333 (x1.23)** |
| fit (train set) | ... + trained attachment | 6.74 (x1.05) | 0.332 (x1.23) |
| LIGHT now | LIGHT today | 3.67 (x0.57) | 0.341 (x1.26) |
| LIGHT now | Beam EEDF + exact capture | 3.46 (x0.54) | 0.333 (x1.23) |

- The EEDF correction improves g0 by 7 % and formation by 2 %.
- The trained cross section moves g0 by 0.6 %. The whole spread of the
  conflicting attachment data (Chen's scale to Trainor & Jacob's 500 K
  scale) moves it by -3 % to +4 %, and formation by 1.5 %.
- The capture fix alone makes both scores worse (g0 7.40, formation
  0.354). It only helps together with the emptied tail.
- Formation stays 23 % high. The rest of that error is not on the
  electron side.

**Adopted in LIGHT (2026-09-28, step 3).** LIGHT's basic network now books
the exact captured energy of each electron it removes, and integrates its
Kr rows over a two-temperature EEDF whose tail temperature comes from
`ml/tail.py`'s table of kinetic solutions (66 LIGHT states, 0 to 3 % F2,
300 W/cm3 to 3 MW/cm3). The table's axes are n_e/N and the F2 fraction.
Te does not carry the tail: pure Kr at 300 W/cm3 has r 0.27 where 0.1 %
F2 at the same n_e/N has 0.035. On its own closure LIGHT lands on this
loop's fixed point at Tahoma:

| Tahoma, fit set | Te (eV) | n_e (cm^-3) | g0 (%/cm) | formation |
|---|---|---|---|---|
| LIGHT, Maxwellian closure | 1.126 | 4.66e14 | 6.83 | 0.341 |
| LIGHT, two-temperature closure | 1.595 | 6.38e14 | 6.39 | 0.335 |
| this loop | 1.606 | 6.53e14 | 6.33 | 0.334 |

The loop runs LIGHT on its Maxwellian closure (`forward.MAXWELL`) and stays
the reference that `tests/test_electron.py` holds LIGHT's closure to.

## Results for KrF formation (2026-09-28)

Owner direction: set the cross-section and EEDF training aside, and make KrF
formation itself the target. `python3 ML/run.py formation` prints all of
this.

**The observable.** Tahoma's formation efficiency is g0·Isat/Pin = 0.27. Its
Isat, 1.4 ± 0.1 MW/cm², comes from a Rigrod fit to the extraction data
(Spectra Technology 1989, p. 3-16). LIGHT's formation rate × 5 eV / Pin is a
different number, because KrF that forms in C, or mixes into it, and decays
there never meets the laser field. `ml/formation.py` puts a 248 nm probe
flux on the model and fits its Isat at fixed g0, as Tahoma did.
- On LIGHT's chain and closure of 2026-09-28 every B/C set gives nearly the
  same formation, 0.298-0.306 (×1.10-1.13), while their g0 differ ×2. The
  excess over 0.27 is in the formation chain, not in B and C.
- Isat is a second test, and no fit used it. The adopted set gives 1.54
  MW/cm² (×1.10). The constants until 2026-09-27 give 3.03 (×2.16).

**The chain against its measurements.** Checked means re-read first-hand.
LIGHT's column is its basic network since the four steps of 2026-09-28.

| Link | LIGHT | Measured | Source (status) |
|---|---|---|---|
| Kr* + F2 rate | 8.47e-10, the mean of the three measurements | 7.2, 8.8 ± 0.5, 9.4 ± 1.0 (e-10) | Velazco 1976 preprint, Ueda 1984 (checked) |
| same, the 8.1e-10 LIGHT used until 2026-09-28 | | not in the report Kannari and Ueda credit it to, which gives 7.2e-10 | SRI MP-76-99 Table 2 (checked) |
| KrF* share of Kr* + F2 | 0.90 | 0.90: the KrF* formation rate constant, 65e-11, over the total quenching, 72e-11. Published as 1.0 ± 0.20 | Velazco, Kolts & Setser, JCP 65, 3468 (1976) p.3476, Tables II and IV (checked) |
| same, the absolute scale | | from Ar(3P2) + Kr and absorption-measured [Ar*] and [Kr*]: ± 20 % spectral response, ± 30 % concentration ratio. KrCl* 0.90, Cl2* 0.10 | same, Sec. III C (checked), quoted by Velazco's dissertation (1977) p.148 |
| Share of that KrF formed in C | 0.38 | 0.38 (k_B/k_C ≈ 1.6), nascent | Yu & Setser 1990 (checked) |
| Kr* + 2Kr → Kr2* | 1.7e-32 | 1.7 ± 0.2 (e-32), and 1.8, 2.6, 4.2, 4.4 | Ueda 1984 (checked), four abstracts |
| its Tg dependence | (Tg/300)^-0.64 | Ea = −0.020 eV | Wieme & Lenaerts 1981 (abstract) |
| Kr2* radiative | 3.7e6 s⁻¹ | 3.7e6 s⁻¹ (270 ns), and five more at 250-353 ns. Chen & Payne ~1 µs | Morikawa 1989 and others (abstracts), Chen & Payne 1979 (checked) |
| Kr2* + F2 | 2.1e-10, split evenly | 2.1 ± 0.2, 3 ± 0.9 and 5 ± 1.25 (e-10). The direct Kr2F* channel is minor (Chen & Payne) | Quigley & Hughes 1978 (abstract), SRI MP-76-99 and Chen & Payne 1979 (checked) |
| Kr2* + F → KrF | 3e-10 (Kannari's estimate) | none exists. SRI only interprets a slow KrF* component as this reaction | SRI MP-76-99 p.42 (checked) |
| ion-ion rate | Flannery ≈ 2.9e-6 | none with F-. Analog: Kr2+ + SF6- about 3e-6 near 665 Torr. Avco's model input 1.1e-6 | Mezyk, Cooper & Sherwell 1991 (checked), Avco 1978 Table 1 (checked) |
| KrF* share of ion-ion | 1 | none with F-. Analogs: KrCl* 0.64 of Kr+ + Cl- products, 0.41 KrF* per pair from Kr2+ + SF6- | Tsuji 1992, Cooper 1996 (abstracts) |
| its B/C split | all B | none with F-. Analogs: 0.62:0.38 (SF6-), 0.60:0.40 (Cl-) | Tsuji 1990, 1991 (abstracts) |

**LIGHT's chain, every measured link at its measured value**
(`formation_data.MEASURED`):

| B/C set | g0 (6.39) | Isat (1.4) | g0·Isat/Pin (0.27) |
|---|---|---|---|
| LIGHT now (the fit) | **6.11 (×0.96)** | **1.54 (×1.10)** | **0.284 (×1.05)** |
| LIGHT to 09-27 | 3.17 | 3.03 | 0.291 (×1.08) |

- The measured change that mattered most is the nascent split. 38 % of the
  KrF from Kr* + F2 forms in C, and with the fitted set's slow B-C mixing
  little of it reaches the laser.
- The Kr2* route is not a loss once its measured 270 ns lifetime is used.
  Kr2* reacts before it radiates. REBUILD's 12.5 ns made it look like a 4 %
  loss.
- LIGHT's electron closure (step 3) moved g0 from ×1.07 to ×1.00 and the
  formation from ×1.13 to ×1.10.
- The measured harpoon yield (0.90, Velazco, Kolts & Setser 1976,
  supplied by the owner) moved g0 to ×0.96 and the formation to ×1.05.

**What remains.** Step 4 searched for the links with no measurement
(`LIGHT/docs/krf_data/FORMATION_LINKS_SEARCH.md`), and the papers the
owner supplied then gave the harpoon yield (Velazco, Kolts & Setser 1976).
One at a time from LIGHT's chain (fit set):

| Link | g0·Isat/Pin |
|---|---|
| KrF* share of Kr* + F2 at 0.70 (0.90 less the stated 0.20) | 0.256 (−0.028) |
| KrF* share of Kr* + F2 at 1.0 (the published assignment) | 0.298 (+0.014) |
| KrF* share of Kr* + F2 at 0.5 (the ArF* analog) | 0.228 (−0.057) |
| KrF* share of ion-ion recombination 0.64 (the Cl- analog) | 0.241 (−0.043) |
| ionic KrF formed 38 % in C (the SF6- analog) | 0.266 (−0.018) |
| no Kr2* + F → KrF | 0.268 (−0.017) |
| ion-ion rate at Avco's 1.1e-6 | 0.284 (no change) |
| nascent C relaxes into B at 1 atm | 0.307 (+0.023) |

- The measured harpoon yield, 0.90 with 0.70 to 1.0, brackets Tahoma's
  0.27. The
  still-unmeasured links (the ionic yield and B/C split with F-, Kr2* + F)
  have analogs only, and every analog points toward the measurement.
  Under the predictive-only policy none is chosen by its fit to Tahoma,
  and an analog never sets a link.
- Other formation measurements agree with Tahoma: 24 % (Kannari, Shaw &
  O'Neill 1987, 99.7 % Kr, abstract) and about 25 % (Ueda 1988, Ar-rich,
  g0·Is/Pin constant from 0.1 to 1.5 MW/cm³, checked).

## Verification on pure rare gases (2026-09-28)

Owner question: how can we verify this is done correctly, using pure gases
such as Ar, Ne or Kr? KrF formation runs through three segments:
1. beam energy into Kr+ and Kr* (W and the excitations per ion);
2. Kr+ and Kr* into Kr2* (ion conversion, recombination, three-body
   excimer formation, VUV emission);
3. the F2 steps.
Pure gases test the first two on their own measurements.
`python3 ML/run.py pure` prints all of this (5 s).

| Check | LIGHT | Pure-gas measurement | Verdict |
|---|---|---|---|
| Kr* + Kr+ made per eV, Kr | 0.0656 | 0.0640 (Saito 2003 N_ex/N_i 0.55, W 24.2) | within 3 % |
| same, Ar | 0.0592 | 0.0576 (0.52, W 26.4) | within 3 % |
| Kr2+ + e, Te 0.03-2 eV | Kannari, Te^-1/2 | Shiu & Biondi 1977 | x0.87-1.08 |
| Kr+ + 2Kr -> Kr2+ | 2.5e-31 | 2.7-3.3e-31 (Smith et al. 1972) | within 25 % |
| thermal-electron Kr excitation, low pump | 0.060 Kr* per pair (0.14 on the Maxwellian closure) | no measurement. The kinetic loop gives about half of LIGHT's | about x2 |
| pure-Kr e-beam VUV yield, low pump | 0.570 | 0.42 +- 0.05 (Morozov 2008) | x1.36 |

- REBUILD's Spencer-Fano gets W 11 % high for Kr and Ar (19 % Ne, 38 %
  He) and N_ex/N_i about 40 % high. The two errors cancel in what
  formation uses, the Kr* + Kr+ made per eV. Kannari's older ledger
  (N_ex/N_i 0.28) is 15-17 % low.
- Without the thermal excitation, LIGHT's pure-Kr yield at low pump is
  0.551, the one-photon-per-excitation-and-ion ideal (0.54 with the measured
  W and N_ex/N_i, and Saito's "one scintillation photon ... from one
  recombination process"). No rate in LIGHT removes the further 22 % that
  Morozov's e-beam measurement shows. The e-beam yields fall below the
  ideal in every gas: Ne about 52 %, Ar 59 %, Kr 78 %, Xe 80 %. A shortfall
  growing toward the light gases looks like a common definition or
  measurement effect more than Kr kinetics. The full text (the owner's
  copy, read 2026-09-28) defines the efficiency as all the VUV power
  emitted into 4 pi over the power deposited in the gas. For Kr the
  apparent 0.35 +- 0.02, divided by GEANT4's deposited fraction 0.67
  (backscatter included) and corrected for reflection on the foil, gives
  0.42 +- 0.05 (`pure_data.MOROZOV_COUNTS`). So the shortfall is not a
  definition effect. It is physics or a systematic error (GEANT4's 10 % is
  the largest term).
- Thermal excitation is where LIGHT's Maxwellian closure broke. Since
  2026-09-28 (step 3) LIGHT takes the tail from the kinetic EEDF, and the
  thermal excitation at 0.33 kW/cm3 fell from 0.14 to 0.060 per pair. It is
  not zero, as this section said on 2026-09-28 (b): with no F2 to slow
  them, e-e collisions lift some injected electrons over the break. The
  kinetic loop puts the Kr(5s) rate at about half of LIGHT's new value. In
  pure Kr at this density the bulk itself is not Maxwellian (LIGHT's Te
  0.88 eV against the kinetic 0.78-0.82), so LIGHT's one-temperature energy
  equation cannot follow it there. At Tahoma it agrees to 0.3 %.
- What the thermal excitation does to KrF formation (LIGHT's chain, fitted
  set): without the thermal 5p/4d/5d/ionization nothing changes (0.298), as
  the tail already empties them. Without any thermal excitation 0.261,
  which is not physical at Tahoma's density.
- Tahoma's pure-Kr run (absorption 0.67 %/cm at 257.3 nm) cannot be used
  yet. The Kr2+ and Kr2* cross sections there are theory or assumption.
- Ar and Ne have no chemistry in LIGHT. They enter through segment 1 (W and
  N_ex/N_i) and the pattern of Morozov's yields.

## Pure krypton: LIGHT_Kr's cross sections trained (2026-10-15)

Owner direction (2026-10-14): "Let us try to improve like we did for
LIGHT/, we can also use ML to train it to improve the accuracy", for
LIGHT_Kr (the pure-krypton model). LIGHT_Kr v0_2's largest scored error
was the electron swarm. Swarm data measure the cross sections themselves,
so they may train them (the rule above). `python3 ML/run.py kr-train`
trains LIGHT_Kr's cross sections and writes the factors file that
LIGHT_Kr v0_3 pins (`LIGHT_Kr/v0_3/data/kr_swarm_trained.json`).

**The data** (`ml/kr_data.py`, `python3 ML/run.py kr-data`): 28 krypton
datasets, 441 points. They come from Dutton's printed tables (Pack 1962
drift velocity, Kruithof 1940 and Heylen 1971 alpha/N) and from LXCat's
LAPLACE database, read with the owner's approval of LXCat's terms
(2026-10-14) and stored unedited and md5-pinned. LAPLACE holds:
- drift velocities: Hunter 1988, Pack 1992, Nakamura 1990, Dall'Armi 1992,
  Bowe 1960, English & Hanna 1953, Kucukarpaci 1981;
- characteristic energies: Koizumi 1986, Al-Amin 1987;
- ionization coefficients: Specht 1980 (direct and total), Jacques 1986,
  Bhattacharya 1979, Dall'Armi 1992.

No file prints an error, so each dataset carries a class value: 3 % for a
drift velocity measured after 1985, 5 % before, 7 % for a characteristic
energy, 10 % for alpha/N. A 2 % error in E/N enters through the model's
slope. Roles:
- **check**: points above 300 Td, Specht's total coefficient, and
  Al-Amin's D_T/mu above 1 Td.
- **held out by rule**: alpha/N below 50 Td (LIGHT_Kr's Known-bad list).
- **held out**: D_T/mu below 1 Td (Koizumi), and the sets the consistency
  check drops. Koizumi's and Pack's two 195-K sets chose which factors
  train (below), so they are validation sets, not independent tests.
- **train**: the consistency check's choice among the rest.

**The forward model** is LIGHT_Kr's own two-term swarm (`swarm.py`) on the
factors (`xsec.using`), times LoKI-MC's ratios (`LIGHT_Kr/v0_3/swarm_mc.py`).
The bulk drift velocity is what a drift tube measures, the SST coefficient
what a Townsend experiment measures, and the Monte Carlo's D_T/mu replaces
the two-term one, which is 46-79 % high from 10 to 300 Td. The two-term
flux drift velocity matches the Monte Carlo's to 0.3 % from 10 to 300 Td.
At 0.2 Td it is 1.8 % off, a two-term error that the factors absorb (they
are effective for LIGHT_Kr's two-term solver). On the trained set the
Monte Carlo gives the same ratios to 1.65 % from 30 to 300 Td (the most
is alpha/N at 50 Td, where its counting error is about 1 %).

**The factors** (natural logs z, prior sd 0.5): the elastic momentum
transfer at six knots (0.03, 0.2, 0.6, 2, 6, 20 eV, log-linear between
them) and the excitation of the 1s, 2p and remaining levels. Ionization
is not trained. `kr_swarm.FREE` trains seven: the knots at 0.03 and
0.2 eV stay at 1 (Biagi's file).

**Which factors train** (`kr_train.free_variants`): each variant's least
squares on the train set, chi2 per point on the true solver.

| factors trained | train sets, pooled | Hunter W | Koizumi D_T/mu | Pack 1962 W, 195 K | Pack 1992 W, 195 K |
|---|---|---|---|---|---|
| none (Biagi's file) | 1.75 | 0.79 | 2.52 | 3.07 | 6.12 |
| the three excitations | 1.16 | 0.79 | 2.52 | 3.07 | 6.12 |
| these + momentum transfer at 2-20 eV | 0.66 | 0.75 | 2.43 | 2.89 | 5.29 |
| these + 0.6 eV (**v0_3**) | 0.64 | 0.70 | 0.84 | 2.77 | 5.27 |
| all nine | 0.49 | 0.13 | 8.06 | 5.45 | 8.70 |

With all nine free the fit lowers the momentum transfer at the Ramsauer
minimum (0.6 eV) to 0.88 and raises it at 0.03 eV to 1.13. The train
sets allow that trade, and the data that weigh those energies differently
reject it. v0_3's variant fits the three held-out columns best of the
five and keeps every train set within the threshold. Those three sets
made the choice, so they validate it. The beam's partition and alpha/N
below 50 Td saw neither the training nor the choice.

**The emulator** (`ml/kr_nn.py`). A 3x160 tanh network maps (z, ln E/N,
Tg) to the change from Biagi's file. It is trained on the true solver at
81 nodes for 2000 Latin-hypercube designs in |z| < 0.6 (16 min on 6
workers). Held-out rms in ln: W 0.10 %, D_T/mu 0.13 %, alpha/N 0.74 %.
Every answer is re-checked on the true solver.

**Consistency** (`kr_swarm.consistency`). Fit every candidate together,
drop the worst while its chi2 per point at the joint fit exceeds 3, and
refit. Of 16 candidates it dropped, in order (chi2 when dropped):
- English 1953 (74), Kucukarpaci 1981 (45);
- Pack 1992 at 300 K (15), Dall'Armi 1992 W (16), Pack 1992 at 195 K
  (5.4), Pack 1962 at 300 K (5.0);
- Pack 1962 at 195 K (3.6), Bowe 1960 (3.5).

Train: Hunter 1988, Nakamura 1990 and Pack 1962 at 368 K (W), and
Bhattacharya, Dall'Armi, Jacques, Kruithof and Heylen at 50-300 Td
(alpha/N). Threshold 2 gives the same eight. At 4, Pack 1962 at 195 K and
Bowe 1960 join. Bowe's W lies 9 % below Hunter's on average where both
exist, and it pulls m_6 to 1.31 and x_1s to 1.16. A leave-one-out
criterion was tried first and dropped Hunter 1988 (see the Known-bad
list).

**The fit** (least squares; MCMC median [16, 84 %]; neural posterior
estimate):

| factor on | LSQ | MCMC | network |
|---|---|---|---|
| momentum transfer, 0.03 eV | 1 (not trained) | | |
| 0.2 eV | 1 (not trained) | | |
| 0.6 eV (Ramsauer minimum) | 1.157 | 1.157 [1.062, 1.260] | 1.125 |
| 2 eV | 0.926 | 0.926 [0.891, 0.962] | 0.909 |
| 6 eV | 1.039 | 1.040 [1.008, 1.071] | 1.038 |
| 20 eV and up | 1.227 | 1.221 [1.106, 1.354] | 1.240 |
| 5s (1s) excitation | 0.886 | 0.871 [0.776, 0.963] | 0.838 |
| 5p (2p) excitation | 0.925 | 0.954 [0.693, 1.293] | 1.047 |
| higher levels | 0.753 | 0.738 [0.631, 0.854] | 0.703 |

MCMC acceptance was 0.43. The network's held-out R^2 is 0.15 for the 5p
factor (the data barely pin it) and 0.83-0.99 for the rest. Its
+-1 sigma coverage is 0.60-0.65.

**Scores** on the true solver (chi2 per point, Biagi's file -> trained):

| | Biagi's file | trained |
|---|---|---|
| Hunter W | 0.79 | 0.70 |
| Nakamura W | 2.93 | 0.39 |
| Pack 1962 W at 368 K | 0.82 | 0.61 |
| Jacques alpha | 3.19 | 1.16 |
| Kruithof alpha | 3.28 | 0.76 |
| Heylen alpha | 0.91 | 0.17 |
| Dall'Armi alpha | 1.76 | 0.76 |
| Bhattacharya alpha | 0.43 | 1.25 |
| Koizumi D_T/mu (validation) | 2.52 | 0.84 |
| Pack 1962 W at 195 K (validation) | 3.07 | 2.77 |
| Pack 1992 W at 195 K (validation) | 6.12 | 5.27 |
| Pack 1992 W at 300 K (held out) | 14.2 | 17.7 |
| Dall'Armi W (held out) | 10.3 | 13.4 |
| Kucukarpaci W (held out) | 39.6 | 48.4 |
| Heylen alpha below 50 Td (held out) | 3.10 | 1.08 |
| Specht direct (held out) | 15.0 | 10.4 |
| Kruithof alpha below 50 Td (held out) | 69 | 60 |

Training moves five of the eight held-out drift velocities closer and
three (Pack 1992 at 300 K, Dall'Armi, Kucukarpaci) further: the train
sets, which agree with one another, disagree with those.

Under an electron beam the trained excitation also moves the
Spencer-Fano partition toward Saito et al.'s measurement, which neither
the training nor the choice of factors saw: W 27.6 -> 26.2 eV (23.9
implied) and N_ex/N_i 0.79 -> 0.68 (0.55).
LIGHT_Kr's accuracy overview v0_3 has the full scores and the reference
point.

## Next

The owner's four steps of 2026-09-28 are done, and the five papers the
search could not reach have been read. LIGHT carries the fitted B/C set,
the measured formation links with the measured harpoon yield (0.90,
published as 1.0 +- 0.20), and the exact captured energies with the
kinetic tail. Tahoma now scores g0 ×0.96, Isat ×1.10 and g0·Isat/Pin
×1.05, and the harpoon yield's measured range brackets the formation.

1. **Still unmeasured:** the KrF* yield and B/C split of ion-ion
   recombination with F- (every measurement found uses SF6- or Cl-), and
   Kr2* + F -> KrF*.
2. **New datasets** go into `ml/data.py` (kinetics) or `ml/swarm.py`
   (electron swarms) with source, error and role, and the consistency check
   shows whether they agree. After any change to LIGHT's network or cross
   sections, `python3 ML/run.py tail` regenerates its tail table.

Pure Kr at low pump stays outside LIGHT's one-temperature closure (the bulk
is not Maxwellian there). It affects the pure-gas check only, not KrF.
