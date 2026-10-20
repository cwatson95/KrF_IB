# Changelog: ML (fitting LIGHT's channels to the measurements)

Dated entries, newest first. Read the Known-bad list before editing.

## Known-bad: do not reintroduce

- **Device data never trains.** Laser gain and efficiency score the model
  and never fit it (the repo's predictive-only policy,
  `KrF_and_LoKI/agents/CLAUDE.md`). `ml/data.py` gives them the role
  `device`, `fit.fit` refuses them, and `tests/test_data.py` fails if one
  is given `train`. Change the policy text first, and only on the owner's
  decision.
- **A fitter cannot settle conflicting experiments.** Fed the train set plus
  Hakuta, the neural network returned a confident compromise that misfits
  Kansas State at chi2 15 per point, with no warning. Run
  `run.py consistency` before choosing a train set.
- **Sets built with k_e = 0 break the Stern-Volmer lifetime.** `sv_tau`
  divides by k_e, so it returned 0/0 = NaN, and the device run crashed on a
  NaN Jacobian. `forward.with_rederived_k_e` computes the lifetime with a
  unit placeholder, which is valid because the lifetime is independent of
  k_e's size.
- **numpy 2 on macOS Accelerate** raises spurious divide/overflow warnings
  in matmuls whose operands hold exact zeros, and the network's tanh layers
  hit this. `ml/nn.py` silences them inside the network only and asserts
  every weight is finite.
- **A swarm EEDF is not the laser's EEDF.** LoKI-MC and the two-term solver
  have no electron-electron collisions. At Tahoma's ionization degree
  (2e-5) those collisions make the bulk Maxwellian. The swarm EEDF of the
  same plasma puts attachment at x0.29 and n_e at x2.3. Use
  `ml/boltzmann/ee.py` (e-e on) for any laser-plasma EEDF.
- **The two-term solver fails where attachment outruns energy relaxation.**
  He + 1 % F2 at 3 Td gave a negative mean energy, and Kr/F2 99.5/0.5 at
  0.3 Td a mean energy x7 below LoKI-MC's. `ml/swarm.py` models He/F2 only
  from 6 Td, where it matches LoKI-MC to 1 %.
- **LIGHT's formation rate x 5 eV / Pin is not Tahoma's formation
  efficiency.** Tahoma's is g0 Isat / Pin with a measured Isat. KrF that
  forms in, or mixes into, C and decays there never meets the laser field.
  Compare with `formation.observables` (probe flux, Isat fitted at fixed
  g0), which gives the same 0.328 for B/C sets whose g0 differ x2.
- **REBUILD's Kr2* radiative rate, 8e7 1/s (Kannari), is 20x the measured
  3.7e6 (270 ns, and five more measurements at 250-353 ns).** With it the
  Kr2* route looks like a 4 % formation loss. Measured, Kr2* reacts first.
- **Setser's "unity" KrF* branching is relative.** It was measured against
  Kr* + Cl2 -> KrCl*, taken as 1. That reference was itself measured at 0.7
  "+-50 %" (COO-2807-3, p.21), and the group's ArF* check came out 0.5.
  Treat the absolute KrF* yield of Kr* + F2 as unmeasured.
- **REBUILD's Spencer-Fano W and N_ex/N_i are each wrong for pure gases**
  (W +11 % for Kr and Ar, +19 % Ne, +38 % He, and N_ex/N_i +40 % for Kr
  and Ar),
  but their combination, Kr* + Kr+ per eV, is within 3 % of the measured.
  Use only the combination. Do not feed W or N_ex/N_i separately into
  anything.
- **A two-temperature EEDF with a chosen tail temperature is a knob.** The
  tail temperature that reproduces the kinetic EEDF (ml/eedf.py) is 0.25 Te
  at Tahoma and 0.10 Te at a tenth of the pump. Take it from the kinetic
  solution at each state, never fix it (REBUILD 0.1.37's two-T tier is on
  REBUILD's Known-bad list for the same reason).
- **Count a VUV yield to the end of the afterglow.** At low pump Kr2+
  recombines over microseconds and Kr2* lives 270 ns, so a 200-400 ns
  window read 0.14-0.17 where the full count is 0.59.
- **The tail table's states run on LIGHT's Maxwellian closure.**
  `ml/tail.py` computes each state with `forward.MAXWELL`. On LIGHT's
  two-temperature closure the table would depend on itself. The kinetic
  loop (`eedf.self_consistent`) runs on it too, or its Boltzmann factors
  would correct LIGHT's tail a second time.
- **Te does not carry the tail ratio.** A polynomial surface in (n_e/N,
  Te) fitted to the kinetic states put r = 1 in corners no state reaches,
  and pure Kr (r 0.27 at 300 W/cm3) sits far from 0.1 % F2 (0.035) at the
  same n_e/N. The table's axes are n_e/N and the F2 fraction, each column
  one series of states, nothing fitted.
- **F2 attachment data conflict between groups.** Near 1 eV Chen (1977),
  Trainor & Jacob (1979) and Schneider & Brau (1978) spread x3.5. Train on
  the consistent pair (Nygaard, Trainor & Jacob 300 K) and hold the rest
  out. The N2 data's mean energies were computed by the authors, not
  measured.
- **ML binds its LIGHT version when first imported.** `LIGHT_VERSION` is
  read once, in `ml/__init__.py`. Setting it after `ml` is imported does
  nothing, and setting it at another module's import time rebinds ML for
  every later caller in that process. Run a second version's numbers in a
  process of their own.
- **ML reaches LIGHT_Kr only through `ml/kr_light.py` (2026-10-15).**
  `from v0_3 import ...` inside ML loads LIGHT's KrF v0_3, because `ml`
  puts LIGHT/ first on sys.path. A scratch script that mixed the two got
  LIGHT's xsec (no `Factors`).
- **Choose a train set by the joint fit's chi2, not by leave-one-out
  (2026-10-15).** With six momentum-transfer knots, a dataset that alone
  pins some knots predicts badly when left out. The leave-one-out pruning
  dropped Hunter 1988's drift velocity first and kept 1962's.
  `kr_swarm.consistency` drops by the chi2 at the joint fit and reports
  leave-one-out for information.
- **The emulator is valid only inside its box (|z| < kr_swarm.BOX).**
  The first fit (box 0.45) put the momentum transfer at 0.6 eV on the
  box's edge. The box was widened to 0.6 and resimulated, and the factor
  came back to 0.88 inside it. A least-squares answer on the bound means
  widen and resimulate, never accept.
- **A swarm's D_T/mu needs the Monte Carlo above 1 Td.** The two-term
  transverse diffusion is 46-79 % high from 10 to 300 Td, so Al-Amin
  1987 is a check of the correction, never a training set.
- **Freeing every factor lets the train set trade knots it cannot
  separate (2026-10-15).** With all nine free, least squares raised the
  momentum transfer at 0.03 eV to 1.13 and lowered it at the Ramsauer
  minimum (0.6 eV) to 0.88. The train sets fit best (pooled chi2 0.49),
  but held-out data got worse: Koizumi's D_T/mu went from chi2 2.5 to
  8.1, and Pack's 195-K drift velocities from 3.1 to 5.5 (1962) and from
  6.1 to 8.7 (1992). `kr_swarm.FREE` holds the knots at 0.03 and 0.2 eV
  at 1, and `kr_train.free_variants` keeps the comparison in the record.
- **A set that chose the factors is a validation set, not a test
  (2026-10-15).** Koizumi's D_T/mu and Pack's 195-K drift velocities
  decided `FREE`. Never report their improvement as independent
  evidence. The tests nothing chose on are the beam's partition and
  alpha/N below 50 Td.
- **A script that uses kr_swarm's process pool needs an
  `if __name__ == "__main__":` guard (2026-10-15).** macOS starts
  workers by spawn, which re-imports the main module. An unguarded
  scratch script re-ran itself in every worker, recursively, until the
  machine's load reached 33.

## 2026-10-15: LIGHT_Kr's cross sections trained on krypton's swarm data

The owner (2026-10-14), on LIGHT_Kr: "Let us try to improve like we did
for LIGHT/, we can also use ML to train it to improve the accuracy".
LIGHT_Kr v0_3 pins the result (LIGHT_Kr/CHANGELOG.md 2026-10-15). The
README has a new section, "Pure krypton".

### Added
- `ml/kr_light.py`: loads a LIGHT_Kr version under its own name,
  light_kr_<v> (LIGHT_KR_VERSION, default v0_3). LIGHT_Kr's v0_1 to v0_3
  share names with LIGHT's packages, which `ml` puts on sys.path.
- `ml/kr_data.py`: 28 krypton datasets with roles (train / held_out /
  check) and DOIs, from Dutton's tables and LXCat's LAPLACE database.
- `data/swarm/Kr_*_LAPLACE_LXCat.txt`: LXCat's text output, unedited
  and md5-pinned. Read with the owner's approval of LXCat's terms of use
  (cite, no redistribution, non-commercial) on 2026-10-14.
- `ml/kr_swarm.py`: the factors on LIGHT_Kr's cross sections, the true
  solver (parallel), the Monte Carlo ratios (from LIGHT_Kr's
  `swarm_mc.py`), the error model, least squares, the consistency check,
  the MCMC and the record LIGHT_Kr pins. `FREE` names the factors that
  train (seven of nine: the momentum transfer at 0.03 and 0.2 eV stays
  at 1). `fit`, `mcmc` and `kr_nn.npe` work in the free components and
  return the full vectors.
- `ml/kr_nn.py`: the emulator, simulations cached under results/ by a
  code key, and neural posterior estimation.
- `ml/kr_train.py`: the pipeline (setup, consistency, fit, record,
  threshold scan, and `free_variants`: five choices of free factors
  fitted and scored on the true solver).
- `run.py`: kr-data, kr-consistency, kr-fit, kr-train.
- `tests/test_kr.py` (7, fast) and `tests/test_kr_slow.py` (4, ML_KR_SLOW=1).

### Found
- The train set, the factors and the scores are in the README.
- Which factors train. The first training freed all nine. It fitted the
  train sets best but degraded three held-out sets (Koizumi's D_T/mu
  chi2 2.5 -> 8.1, Pack's 195-K drift velocities 3.1 -> 5.5 and
  6.1 -> 8.7). Of five variants, holding the knots at 0.03 and 0.2 eV
  at 1 fits those three best (0.84, 2.77, 5.27) and keeps every train
  set within the threshold, so `FREE` is that variant. The three sets
  are validation sets from here on (Known-bad list). A thermal-limit
  alternative for the D_T correction below 0.2 Td (ratio 1 at 0.005 Td)
  moves Koizumi's chi2 by +0.1 to +1.1 per point and keeps the ranking.
- LoKI-MC on Biagi's file: the two-term W matches the Monte Carlo's flux
  W to 0.3 % from 10 to 300 Td (1.8 % at 0.2 Td, absorbed by the factors),
  and its k_i/W lies 1.4-3.9 % below the flux coefficient.
  The bulk W and the SST coefficient part from them as ionization grows.
  The two-term D_T is 46-79 % high from 10 to 300 Td.
- The two-term swarm's D_T/mu obeys Einstein's relation at 1e-4 Td
  (kTg/e to 2.6e-5).

### Verified
- `ML/tests` on the seven-factor record: 53 passed, 4 skipped (the slow
  krypton tests), 1 xfailed, in 2 min. LIGHT's 45 are unchanged and 8
  are new, among them `test_the_pinned_record_trains_the_free_factors_only`
  (the record's free mask is `FREE`, its fixed factors exactly 1 with no
  spread).
- `ML_KR_SLOW=1 ML/tests/test_kr_slow.py`: 4 passed, in 66 s from a warm
  cache. The emulator's held-out errors are under 1 % (W, D_T/mu) and
  2 % (alpha/N); it matches the true solver within 3 % at a random trial.
  The train set is the consistency check's, and the pinned factors give
  every train set chi2 at most 3 on a fresh true-solver run, equal to
  the record's to 1e-3.
- LIGHT_Kr's overview generator scores the five variants live on the
  solver and agrees with the record's chi2 to 1 % for every W and alpha
  dataset (it raises otherwise).
- The toy tests: the consistency check drops a set planted x1.6 off
  (leave-one-out chi2 above 50) and keeps all when none is.

## 2026-09-29 (c): LIGHT v0_2, and Jancaitis 1983 registered

LIGHT/CHANGELOG.md 2026-09-29 (e) has v0_2 and the findings.
`LIGHT_VERSION=v0_2` runs it, and v0 and v0_1 run as before.

### Changed
- `forward.py`. `tahoma_params` asks for the kr2 group only where it is
  optional. In v0_2 the Kr2* route is basic, so kr2=False cannot remove it.
- `formation.py`. `_impose_chain` sets the Kr2* links whenever their rows
  are in the model, so a chain with k3 = 0 zeroes v0_2's route. The
  chain's Kr2* + F default is the LIGHT version's own (`_light_default`),
  0 in v0_2, where it is an unmeasured band, and Kannari's 3.0e-10 in v0
  and v0_1.
- `pure.py`. `kr_vuv(rates=...)` passes band values to a version that
  carries bands.
- `formation_data.py`. JANCAITIS_1983 and two ABSENT entries. His model's
  KrF* formation has no Kr2* + F (Table II pp.94-95), and it sends all of
  Kr2* + F2 to KrF* at SRI's 3.0e-10 (p.95), an assignment. MEASURED and
  UNMEASURED do not change.
- `tests/test_light_version.py`. The probe reads the chain's Kr2* rows,
  and a v0_2 probe checks no group, v0_1's electron side and 0 for Kr2* +
  F and Kr2* + Kr2*.

### Found
- The owner's step 4, pure krypton under v0_2. At Morozov's pump the VUV
  yield is 0.571 (0.570 in v0 and v0_1), and every Kr2* band moves it by
  0.001 or less, so the 0.42 +- 0.05 miss is not the Kr2* route. At 300
  kW/cm^3 and 665 Torr the bands move it from 0.13 (Kr2* + Kr2* at 5e-10)
  through 0.27 to 0.49 (Kr* + 2Kr at 5.5e-32).

### Verified
- LIGHT and ML together, 258 passed, 3 skipped and 4 xfailed. The v0 and
  v0_1 probes read the same groups and Kr2* rates as before, and ML's v0
  regressions pass unchanged.

## 2026-09-29 (b): the owner's four 1977-78 papers registered

LIGHT/CHANGELOG.md 2026-09-29 (d) has the findings.

### Changed (`ml/formation_data.py`)
- ION_NOTES gain two first-hand entries. Jacob, Rokni, Mangano & Brochu
  1978 measured Kr+ + F- (+M) at 60 and 120 Torr (0.83 and 0.63 of
  LIGHT's screened Flannery form). Rokni, Jacob & Mangano 1977 found the
  Kr2+ + F- -> KrF* branching near unity relative to Kr+ + F-. The comment
  above the list says what is and is not measured.
- Quigley & Hughes 1978 read first-hand. The k_kr2_f2 entry is checked,
  and ABSENT records their statement that the Kr2* + F2 product split
  cannot be discerned.
- ABSENT records SRI MP-76-99's p.45 and p.65 on Kr2* + F.
- VALUES gain k3 = 5.5e-32 (Hughes 1977, secondhand via Quigley &
  Hughes).
- MEASURED and UNMEASURED do not change, so no generated number moves.

### Verified
- ML, 45 tests as before, with the registry test over the new entries.

## 2026-09-29: LIGHT_VERSION, one ML for v0 and v0_1

The owner asked for v0_1's accuracy on v0's tests and a comparison of the
two (LIGHT/CHANGELOG.md 2026-09-29). ML imported `v0` by name, so it now
runs the LIGHT version the environment names.

### Changed
- `ml/__init__.py`: `LIGHT_VERSION` (default v0, refused unless
  `LIGHT/run_<v>.py` exists) and `light(module)`. Every `from v0 import`
  in `ml/` goes through it.
- `forward.py`. `impose` gives each row that runs on tau_sp_s its own
  share of 1/tau_B (v0's G4 alone, v0_1's G4 and G4A). `tahoma_params` asks
  for the ions and f2e groups only where they are optional (in v0_1 they
  are basic rows). `_drive` loads the survey's beam script with its two v0
  imports pointed at the version, so the drive is the same file.
- `formation.py`. The ion-ion share rows copy the version's own K26 and
  K108 rate form, and `flan` keeps the rest of the rate tuple, v0_1's Tg
  exponent included.
- `pure.py`. `kr_rates` reads R11p and K15 from the version's rows (v0's
  values are the ones it had typed), and `kr_vuv` asks for groups the
  version offers.
- `eedf.py` and `boltzmann/ee.py`. The loop reads R11p's Te exponent
  (`r11p_pow`), books (3/2 + n) Te per recombination, and gives the
  Boltzmann solve the matching sigma v ~ eps^n sink (`rec_pow`). v0's
  -1/2 runs the code it ran before.
- `tail.py`. `write_table` refuses to write the v0-made table from another
  version.

### Added
- `tests/test_light_version.py`, three tests in fresh processes. The
  default builds v0's Tahoma model exactly as before, `LIGHT_VERSION=v0_1`
  builds v0_1's, and an unknown version is refused.

### Verified
- v0's numbers did not move. LIGHT's v0 overview tests, fast and slow,
  regenerate their files byte for byte.
- Six bugs planted in the switch, one at a time. Every one failed
  `test_light_version.py`.

## 2026-09-28 (h): Velazco, Kolts & Setser 1976 registered

LIGHT/CHANGELOG.md 2026-09-28 (g) has the findings.

### Changed (`ml/formation_data.py`)
- gamma_H 0.90 now cites the paper first-hand: the KrF* formation rate
  constant, 65e-11 (p.3476), over the total quenching, 72e-11 (Table II),
  +-0.20 from Table IV. It is MEASURED and LIGHT's value, unchanged.
- gamma_H 1.0 is the paper's Table IV assignment (checked), replacing the
  abstract entry.
- The 0.63 entry is gone. The ratio 0.7 it rested on came before the
  recalibration (the published ratio is 0.88).
- The k_H 7.2e-10 entry records the paper's two F2 mixtures (61 and 72).
- The brackets are 0.70 (0.90 less the stated 0.20) and 1.0 (the published
  assignment).
- README: the chain table, the brackets and Next. The papers are all in,
  and the remaining unmeasured links are the ionic yield and B/C split
  with F-, and Kr2* + F.

### Verified
- ML: 41 passed, 1 xfailed. LIGHT: 67 passed, 1 xfailed. No regression
  value moved, since the yield is unchanged.

## 2026-09-28 (g): the owner's four papers registered

The owner supplied the four papers step 4 asked for
(`LIGHT/docs/krf_data/sources/`). LIGHT/CHANGELOG.md 2026-09-28 (f) has
what they change in LIGHT.

### Added (`ml/formation_data.py`)
- VALUES:
  - gamma_H 0.90 from Velazco 1977 p.148, the absolute Kr* + Cl2
    reference that the KrF* "unity" is relative to (checked). It is
    MEASURED and LIGHT's value now.
  - gamma_H 0.63, the raw pre-recalibration ratio 0.7 on the same scale
    (replaces the 0.7 entry).
  - k_H 6.8e-10 and 7.6e-10 for Kr(3P1) and Kr(1P1), a third k_kr2_f2
    (5e-10 +- 25 %) and kr2_rad (~1 us), all from Chen & Payne (checked).
- ANALOGS: Velazco's XeF* 1.0 (relative to XeCl*), and Mezyk's Kr2+ +
  SF6- rate as ion_scale about 1 (about 3e-6 cm3/s near 665 Torr).
- ION_NOTES: Chen & Payne's ionic KrF* light, about 10 % of the neutral at
  low ion density, set by F- lost to F2 (0.7 % of the F- loss at Tahoma).
- ABSENT: Chen & Payne on Kr2* + F2. The direct Kr2F* channel is minor,
  and the KrF* share is not given.
- `ml/pure_data.py`: MOROZOV_COUNTS from the full text (all the VUV into 4
  pi over the power deposited in the gas, 0.35 apparent -> 0.42 physical),
  and the Morozov values are now checked.

### Changed
- The Chain default gamma_H is 0.90 (LIGHT's). The brackets are 0.63 (the
  raw low end), 1.0 (LIGHT until today) and 0.5 (the ArF* analog).
- The tail table was regenerated with `run.py tail` after LIGHT's network
  changed.
- Regression values:
  - The device g0 with LIGHT's constants is 6.11 (was 6.37), and 3.17 with
    the constants until 2026-09-27.
  - The measured chain gives g0 Isat / Pin 0.284, g0 6.13 and Isat 1.53.
  - The kinetic loop gives g0 6.07 and formation 0.317. LIGHT's own
    closure lands within 1 % of it (2.3 % in n_e).
  - The Isat test now drops the yield to 0.63, the raw low end.

### Verified
- ML: 41 passed, 1 xfailed. LIGHT: 67 passed, 1 xfailed.

## 2026-09-28 (f): the unmeasured formation links searched (step 4)

Owner decision (2026-09-28, "Proceed with all 4 steps"). The search record
is `LIGHT/docs/krf_data/FORMATION_LINKS_SEARCH.md`.

### Added
- `ml/formation_data.py`:
  - ANALOGS: seven values from neighbouring processes (ArF*, Kr* + SF6,
    Kr+ + SF6- and Cl-, Kr2+ + SF6-). They set brackets, never a link.
  - ABSENT: the two links with no measurement, and why (Kr2* + F -> KrF*,
    and the KrF*/Kr2F* split of Kr2* + F2).
  - SRI's Kr2* + F2 value, 3e-10 +- 30 % (checked).
- `ml/pure_data.py`: MOROZOV_COUNTS, what Morozov's efficiency counts.

### Changed
- MEASURED k_H and the Chain default are 8.47e-10, the mean of the three
  first-hand values. The Nakano entry says why 8.1e-10 is not used.
- The UNMEASURED brackets use the analogs:
  - The harpoon KrF* yield at 0.5 (ArF*) is new.
  - The ionic C share is 0.38 (SF6-). It was 0.5 (statistical).
  - The ionic KrF* yield is 0.64 (Cl-). It was 0.8.
- README: every results table regenerated (score, formation, pure, fit).
  The electron-side and formation sections now describe LIGHT's adoption,
  and Next is rewritten.

### Found
- Every analog points below LIGHT's value, toward Tahoma's 0.27. One at a
  time from LIGHT's chain, g0 Isat / Pin becomes 0.228 with the harpoon
  yield at 0.5, 0.255 with the ionic yield at 0.64, 0.280 with the ionic C
  share at 0.38, and 0.281 with no Kr2* + F. None of these sets a link.
- Morozov 2008 counts the second continuum over the energy deposited in
  the gas, PTB-traceable to 14 % at k = 2. The integration limits and the
  backscatter handling are in no open source.
- The fit's posterior predictive for Tahoma's g0 on today's network and
  closure is 6.07 %/cm [5.66, 6.93]. It was 7.19 on 2026-09-27's.

### Verified
- ML: 41 passed, 1 xfailed. LIGHT: 67 passed, 1 xfailed.
- The research agent used 3 web searches and sent no email in any
  request. The 16 kept files are checktext CLEAN, rechecked before they
  were copied. It departed from the checktext-first rule twice: it parsed
  one bot-check page, and read 73 API responses, before scanning them (all
  CLEAN afterwards). Three HTML responses are quarantined.

## 2026-09-28 (e): LIGHT's two-temperature tail and exact captured energies (step 3)

Owner decision (2026-09-28, "Proceed with all 4 steps"). LIGHT/CHANGELOG.md
2026-09-28 (d) has the closure itself. ML makes its table and is its
reference.

### Added
- `ml/tail.py`: the tail table for LIGHT.
  - `point`: one state. A LIGHT beam-drive run on the Maxwellian closure
    (pulse only) is followed by its kinetic EEDF and r, the ratio at which
    the two-temperature EEDF removes the kinetic energy through G1, C1a,
    C1b, C1c and C2.
  - `scan`, `table`, `lookup`, `write_table`: 66 states (0, 0.03, 0.1, 0.2,
    0.4, 0.8, 1.6 and 3 % F2 at 300 W/cm3 to 3 MW/cm3, 1 atm) plus four at
    0.5 and 2 atm, written to `LIGHT/v0/data/kr_tail_ratio.txt`.
- `run.py tail` regenerates the table (about 4 min). It is not part of
  `all`, because it rewrites a pinned LIGHT file.
- `tests/test_tail.py` (4). At Tahoma's own state (665 Torr, 2.8 Torr F2,
  not a table state), at its pump and a tenth of it, the table gives the
  kinetic r within 10 % and Kr(5s) within 15 %. The generator reads the
  pinned table exactly as LIGHT does.
- `tests/test_electron.py`: LIGHT on its own closure against the kinetic
  loop at Tahoma (Te 0.993, n_e 0.977, g0 1.009, formation 1.005 of the
  loop's values).

### Changed
- `ml/forward.py`: `MAXWELL` (LIGHT's closure until 2026-09-28).
  `tahoma_params` and `run_device` take `closure` and `afterglow_s`.
  `b_share` is NaN when no KrF forms (pure Kr).
- `ml/eedf.py`: the kinetic loop runs LIGHT on `forward.MAXWELL`, so its
  numbers are unchanged. `beam_eedf` takes `tol`.
- Regression values on LIGHT's new closure:
  - Tahoma g0 with LIGHT's constants is 6.37 %/cm (was 6.81), and 3.32
    with the constants until 2026-09-27.
  - The measured chain gives g0 Isat / Pin 0.298 (x1.10 of 0.27, was
    0.306), g0 6.39 and Isat 1.54.
  - Pure Kr at 330 W/cm3 gives a VUV yield of 0.570 (was 0.595), with
    thermal Kr excitation 0.060 per pair (was 0.14).

### Found
- **Thermal excitation in pure Kr at low pump is not zero** (this
  corrects 2026-09-28 (b)). With no F2 to slow them, e-e collisions lift
  some injected electrons over the break. The kinetic loop gives a Kr(5s)
  window rate of 0.0165 per pair at 330 W/cm3, and LIGHT's new closure
  gives 0.0305. In this regime even the loop does not close LIGHT's
  budget: LIGHT's Te 0.90 eV against the Boltzmann solution's 0.82, and
  n_e 25 % apart. The bulk is not Maxwellian there, while at Tahoma the
  two agree to 0.3 %.

### Verified
- ML: 41 passed, 1 xfailed. LIGHT: 67 passed, 1 xfailed.

## 2026-09-28 (d): LIGHT takes the fitted B/C set and the measured chain (steps 1 and 2)

Owner decision (2026-09-28, "Proceed with all 4 steps"): LIGHT adopts the
fitted KrF(B,C) constants (step 1) and the measured formation links (step
2). LIGHT/CHANGELOG.md 2026-09-28 (b) and (c) have the network changes.
ML follows them.

### Changed
- `ml/params.py`: `gadomski_1992()` is Gadomski's published set (K_eq 5.0,
  fast mixing), and `light_to_0927()` is LIGHT's set until 2026-09-27.
  `light_now()` reads LIGHT, which now carries the fitted set.
- `ml/score.py`: the scorecard compares "LIGHT now" (the fit, rounded),
  "LIGHT to 09-27", "Gadomski 1992" and "fit (train set)".
- `ml/formation.py`: `Chain()` defaults are LIGHT's values, which are now
  the measured ones (`formation_data.MEASURED`). `build` adds a chain row
  only if LIGHT lacks it (LIGHT has its own G3C now), and the chain still
  sets every rate by row id.
- `ml/formation_data.py`: MEASURED takes LIGHT's k_H 8.1e-10 (Nakano's
  measurement). The mean of the three first-hand values, 8.47e-10, moves
  g0 Isat / Pin by +0.1 %.
- `ml/forward.py`: `tahoma_params` runs the kr2 group by default.

### Regression values (Tahoma, today's network)
- `forward.run_device`: LIGHT's constants give g0 6.81 %/cm with B share
  0.53. The constants until 2026-09-27 give 3.50 (3.67 on that day's
  network).
- `eedf.self_consistent` (fit set): g0 6.33 %/cm, formation 0.334. LIGHT's
  Te and n_e stay within 1 % and 4 % of the Boltzmann solution's.
- `formation.observables` (fit set, LIGHT's chain): g0 6.83 %/cm, Isat
  1.48 MW/cm2, g0 Isat / Pin 0.306 (x1.13 of 0.27).

### Verified
- ML: 36 passed, 1 xfailed. LIGHT: 55 passed, 1 xfailed.

## 2026-09-28 (c): would a two-temperature EEDF fix the tail? (no code change)

Owner question: add two tail Maxwellian temperatures to improve this?
Tested by fitting a two-T EEDF (bulk Te, break at 9.915 eV, tail T) to
the kinetic beam EEDF with e-e (ml/eedf.py) in the Tahoma mix:
- 330 kW/cm3 (n_e/N 2.9e-5): tail T 0.40 eV = 0.25 Te. Kr(5s) and Kr(5p)
  rates within 2 % of the kinetic ones, the bulk rates within 4 %.
- 33 kW/cm3 (n_e/N 2.9e-6): tail T 0.19 eV = 0.10 Te. Kr(5s) x1.5, and
  the bulk is no longer Maxwellian (attachment x1.29).
- The far tail (4d+6s, 5d, ionization) is wrong by x2 to x1e5 either way,
  but those rates are negligible.
Verdict: useful only as a compact form of the kinetic EEDF, with the tail
temperature computed at each state (a table in n_e/N and Te) and the exact
capture energies alongside (without them LIGHT reroutes the energy into
Kr(5s), as removing only the 5p-5d rows showed). Then the pure-Kr low-pump
excess goes, and KrF formation at Tahoma moves about 2 % (the kinetic EEDF
result). It does not close the x1.13.

## 2026-09-28 (b): verification on pure rare gases

Owner question: how can we verify this is done correctly, and would pure
gases such as Ar, Ne or Kr show how to improve the modeling?

### Added
- **`ml/pure.py`**:
  - `supply`: REBUILD's Spencer-Fano against the measured W and N_ex/N_i,
    and the Kr* + Kr+ per eV that formation uses;
  - `kr_rates`: LIGHT's Kr2+ recombination and Kr+ conversion against
    pure-Kr afterglows;
  - `kr_vuv`: LIGHT on pure Kr under the beam drive, with the VUV yield
    counted to the end of the afterglow.
- **`ml/pure_data.py`**: the pure-gas measurements with source, place and
  status:
  - W from IAEA-TECDOC-799, Mano 2021 and PDG;
  - N_ex/N_i and W_s from Saito 2003;
  - e-beam VUV yields from Morozov 2008;
  - Kr ion chemistry from Shiu & Biondi 1977 and Smith, Dean & Plumb 1972;
  - Tahoma's pure-Kr absorption (Table 3-2, score only).
- **`ml/formation.py`**: `Chain.thermal_upper` (drop only the thermal
  5p/4d/5d excitation and ionization).
- **`run.py pure`** (5 s) and 6 tests in `tests/test_pure.py`, one of them
  a strict xfail (LIGHT's pure-Kr VUV yield against Morozov's).

### Verified (first-hand)
- The abstracts of Saito 2003 (N_ex/N_i 0.52, 0.55, 0.60 for Ar, Kr, Xe,
  W_s 17.5, 15.4, 13.0 eV, and one photon per recombination), Morozov 2008
  (Ne 0.31, Ar 0.33, Kr 0.42, Xe 0.42 above 500 hPa), Shiu & Biondi 1977
  and Smith, Dean & Plumb 1972.
- IAEA-TECDOC-799 (W: Ar 26.4 eV for electrons, and 26.3, 24.1, 21.9 eV for
  Ar, Kr, Xe gases), Mano 2021 ("The w-value in Kr is 24.2 eV") and
  Tahoma's Table 3-2 (pure Kr, -0.72 and -0.67 %/cm at 257.3 nm).
- Found by a research agent (3 web searches, no email in any request,
  every file checktext'd, 3 unclean files quarantined unopened).
- REBUILD's `beam.py` and `cross_sections.py` in this worktree are
  byte-identical to main's. The Spencer-Fano numbers are the same.

### Results (the README has the table)
- Segment 1 passes: Kr* + Kr+ per eV within 3 % (Kr, Ar), although W and
  N_ex/N_i are each off. Kannari's ledger is 15-17 % low.
- Segment 2 links pass: Kr2+ recombination x0.87-1.08, Kr+ conversion within
  25 %.
- LIGHT's pure-Kr e-beam VUV yield is x1.42 of Morozov's 0.42 +- 0.05.
  - 0.045 of it is thermal-electron excitation, impossible at that electron
    density.
  - The remaining 22 % below the one-photon ideal is in no LIGHT rate.
    Across gases it grows toward the light gases, so it looks like a common
    definition or measurement effect.
- For KrF formation, thermal excitation is recycled sub-excitation energy.
  It is physical at Tahoma's density (kinetic EEDF, about 2 %), so the pure
  gases do not close the remaining x1.13.

## 2026-09-28: KrF formation as the target

Owner direction: ignore the cross-section and EEDF training, and focus on
KrF formation itself as the target.

### Added
- **`ml/formation.py`**:
  - `Chain`, the formation links with LIGHT's values as defaults. They are
    the Kr* + F2 rate with its KrF* and C shares, Kr* + 2Kr -> Kr2* with
    its Tg exponent and the Kr2* fates, the ion-ion rate with its KrF* and
    C shares, and thermal-electron excitation;
  - the chain's rows are added to LIGHT's own network at build time (split
    rows checked by LIGHT's element, charge and energy rules), and LIGHT is
    unchanged;
  - `observables`, a 248 nm probe flux (stimulated emission sigma_B I/h nu
    from KrF(B)) giving g0, Isat (fitted at fixed g0 over 0.5-4 MW/cm2, as
    Tahoma did) and g0 Isat / Pin, with the formation budget per ion pair.
- **`ml/formation_data.py`**: the chain's measurements with source,
  place and status (checked, secondhand or abstract), the measured chain
  (`MEASURED`) and brackets for the unmeasured links (`UNMEASURED`).
- **`ml/data.py`**: `tahoma_isat` (device, 1.4 +- 0.1 MW/cm2, p.3-16). The
  scorecard now scores formation as the model's g0 Isat / Pin.
- **`ml/forward.py`**: `tahoma_params`, `pair_rate`, `WINDOW`, shared by
  `run_device` (unchanged output) and the formation model.
- **`run.py formation`** (2 min), and 8 tests in `tests/test_formation.py`:
  - the default chain reproduces `run_device` exactly;
  - probe power equals the stimulated events times h nu;
  - lasing cools the gas;
  - the harpoon split keeps the total rate;
  - Isat is a lifetime, not a formation number;
  - the measured chain's regression values;
  - the Kr2* lifetime check;
  - the registry is complete.
- **LIGHT/docs/krf_data/sources/**: Setser's ERDA reports COO-2807-3
  (1976), -11 (1977) and -19 (1979) (OSTI), and Avco's Feb 1977 "One meter
  KrF laser system" report (DTIC ADA047684), all checktext CLEAN.

### Verified (first-hand)
- Setser COO-2807-3 p.21: Kr* + Cl2 branching "0.7 ... as much as +-50 %".
- Setser COO-2807-11 pp.5-6: ArF* from F2 "~0.5 rather than unity".
- Setser COO-2807-19 pp.4-9: B/C from 3P1 and 3P2 "very similar", and
  KrF B/C about 0.2 near 1 Torr.
- Yu & Setser 1990, Sec. III.A: k_B/k_C "1.6 is a good estimate".
- Ueda 1984 Tables I, II and V: k3 1.7 +- 0.2e-32, k_H 8.8 and 9.4e-10.
- Avco 1977 p.31 ("essentially unit branching", Ar-rich, model) and Avco
  1978 Table 1 (ion-ion 1.1e-6 as model input).
- Ueda 1988 p.211: formation efficiency about 25 %, constant over 0.1-1.5
  MW/cm3 (Ar-rich).
- The Crossref and OpenAlex abstracts of Velazco 1976, Tracy & Oskam 1976,
  Smith & Turner 1963, Bouciqué & Mortier 1970, Salamero 1979, Wieme &
  Lenaerts 1981, Audouard 1988, Morikawa 1989, Oka 1974, Leichner &
  Ericson 1974 and Quigley & Hughes 1978.
- Found by two research agents (8 + 8 web searches, no email in any
  request, every download checktext'd, the non-clean pages quarantined
  unopened, or opened by mistake and then re-read from clean copies).

### Results (the README has the tables)
- Like for like, LIGHT's formation is x1.22 of Tahoma's 0.27 for every
  B/C set.
- The fitted B/C set predicts Tahoma's held-out Isat, 1.49 against
  1.4 +- 0.1 MW/cm2. LIGHT's current constants give 2.95.
- With every measured link at its measured value, formation is x1.13
  (0.306), g0 x1.07 and Isat x1.06 (fitted set). The measured nascent C
  share (0.38) is the change that matters. The Kr2* route is not a loss
  with its measured lifetime.
- What remains lies in three unmeasured links: the absolute KrF* yield of
  Kr* + F2, and the KrF* yield and B/C split of ion-ion recombination.
  Each alone closes most of it. None is chosen by fit to Tahoma.
- A W-value of 24.1 eV (Avco's table) in place of the beam drive's 27.04
  moves formation by -2 % only (0.322): the beam partition is not the
  lever.

## 2026-09-27 (b): the electron side (cross sections and the EEDF)

Owner question: does training LoKI-MC's cross sections and EEDF improve
the model?

### Added
- **`ml/boltzmann/`**: the two-term Boltzmann solver of the website's
  LXCat sim (`lx_lib.py`, `gases_twin.py`, `boltz_twin.py`), copied
  byte-identical with md5s in `__init__.py`. It reads LoKI-MC's own input
  files, which are byte-identical to LIGHT's pinned copies.
- **`ml/boltzmann/ee.py`**:
  - the electron-electron collision term (the isotropic Landau operator in
    energy space, as in Hagelaar & Pitchford 2005, re-derived);
  - `solve_field` (swarm drive with e-e) and `solve_beam` (the beam drive:
    sub-excitation electrons in, attachment and recombination out).
- **`ml/swarm.py`** and **`data/swarm/`**: 11 swarm datasets with
  provenance.
  - Kr: W from Pack and from Nakamura, alpha/N from Kruithof and from
    Heylen.
  - F2 attachment: Nygaard in He, Trainor & Jacob (300 and 500 K) and Chen
    in N2, Schneider & Brau, thermal values.
  - The forward model is the two-term solver. Trace-F2 rates are linear
    functionals over cached EEDFs.
  - `fit_low` fits one scale on the main attachment resonance.
- **`ml/eedf.py`**: the laser plasma's EEDF at LIGHT's Tahoma state.
  - LIGHT's electron rows are scaled by the ratio to the same-mean
    Maxwellian, and the capture energies are booked exactly.
  - It iterates LIGHT and the Boltzmann solution to a fixed point.
- **`ml/electron.py`**: row scales, trained tables (`att_table`) and
  capture factors for a compiled LIGHT model.
- **`ml/forward.py`**: `run_device` also returns the window's densities,
  event rates and pair rate.
- **`run.py electron`**, and 12 tests in `tests/test_electron.py`:
  - the e-e term's Maxwellian invariance, Spitzer limit and energy
    conservation;
  - the solver against LoKI-MC (He + 1 % F2, 10 Td);
  - the bulk and tail factors, the capture energy, the fixed point's
    regression values and self-check;
  - linearity and the registry.

### Verified
- He + 1 % F2 LoKI-MC grid (3, 6, 10, 15 Td, 2e4 electrons, 6 min). At
  6-15 Td the two-term solver agrees to 1 % in <eps>, W and k_att. It
  fails at 3 Td, and this also explains the Kr/F2 0.3 Td disagreement
  (LoKI-MC's value stands).
- e-e term: the Maxwellian flux residual is under 1e-3, the Spitzer rate
  holds to 0.03 %, and the energy error at Tahoma's ionization degree is
  0.07 % of the input power.
- The fixed point closes LIGHT's energy budget against the kinetic one:
  Te 1.603 vs 1.597 eV, n_e 6.70e14 vs 6.37e14.

### Results (the README has the tables)
- Kr (Biagi): W matches Nakamura (0.96-1.03). alpha/N is low near
  threshold, where the data themselves differ x2. No effect on the laser.
- F2 attachment: only Nygaard and Trainor & Jacob 300 K agree. Trained
  scale 0.921 +- 0.045. The 500 K data rise 25 %, which no cross-section
  set carries.
- The laser's EEDF has a Maxwellian bulk (within 4 %) and an emptied tail
  (Kr(5s) x0.07 down to ionization x2e-7). An attaching electron carries off
  0.54 eV where LIGHT books 1.69 eV, 21 % of the electrons' heating.
- Tahoma, fit set: g0 7.28 -> 6.78 (x1.14 -> x1.06), formation 0.341 ->
  0.333 (x1.26 -> x1.23). The trained attachment adds -0.6 % on g0. So the
  EEDF correction improves the model, and cross-section training barely
  moves it.
- LIGHT itself is unchanged. Adopting the corrections is the owner's call
  (README, Next).

## 2026-09-27: ML created (v0: the KrF(B,C) block)

Owner direction: a subdirectory to train the best models to the available
data.

### Added
- **`ml/data.py`**, the registry of 13 datasets, each with source, error
  model and role:
  - train: Eden Figs. 3 and 5, Gadomski Figs. 1 and 2, Yu & Setser
    I_r(300 K);
  - held out: Xu Figs. 5-7, Hakuta Figs. 2 and 3, Yu & Setser's
    temperature ratio, Burnham & Searles' tau_B;
  - device: Tahoma g0 and formation efficiency.
- **`ml/params.py`**: 10 channels.
- **`ml/forward.py`**:
  - LIGHT's closed-form B/C block for kinetic data;
  - LIGHT under the prototype beam drive for device data;
  - k_e re-derived per set.
- **`ml/fit.py`**: least squares, leave-one-group-out consistency, and an
  ensemble MCMC (Goodman-Weare).
- **`ml/nn.py`**: neural posterior estimation (numpy MLP, Gaussian head).
- **`ml/score.py`** and **`run.py`**:
  - commands: data, consistency, fit, nn, score, all;
  - the whole pipeline runs in 2.5 min.
- **11 tests** (7 s), covering:
  - registry integrity and the device-never-trains guard;
  - Yu & Setser's eq. 18 returning their I_r (0.160), and Gadomski's eq. 5
    returning 4.45e7 at 50 Torr;
  - the re-derived k_e matching LIGHT's QE, and the device forward model's
    regression value (g0 3.67 %/cm);
  - recovery of known constants from synthetic data;
  - the consistency check singling out Hakuta;
  - a small network learning K_eq.
- **Moved here from `LIGHT/docs/krf_data/`:** `bc_infer.py`, `bc_nn.py` and
  `accuracy.py`, now `ml/fit.py`, `ml/nn.py` and `ml/score.py` on the
  registry. The closed form itself moved into LIGHT as `v0/two_level.py`.

### Results (the README has the tables)
- Only Eden, Kansas State and Yu & Setser agree.
- MCMC on them:
  - k_BC 9.8e-12 [4.6e-12, 1.7e-11], K_eq 1.07, tau_C 50 ns;
  - F2 quench: B 4.6e-10, C 3.0e-10;
  - lasing B share at 1 atm 0.69 [0.57, 0.80], against 0.27 in LIGHT now.
- The network on the train set gives the same answer.
- Tahoma's gain (a score only): the fit gives 7.28 %/cm (x1.14), and 7.19
  [6.39, 7.98] over the posterior, against a measured 6.39. LIGHT now
  gives x0.57.
