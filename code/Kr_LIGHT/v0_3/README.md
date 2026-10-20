# LIGHT_Kr v0_3: the pure-krypton model, its cross sections trained by ML

v0_3 (2026-10-15) is v0_2 with its electron cross sections trained on
krypton's measured swarm data by machine learning. The owner asked for
it on 2026-10-14: "Let us try to improve like we did for LIGHT/, we can
also use ML to train it to improve the accuracy". LIGHT's channels were
fitted to their own measurements in `ML/`. Swarm measurements (drift
velocity, characteristic energy, ionization coefficient) measure the
cross sections themselves, so under the predictive-only policy they may
train them. Device and beam data never do.

What changed from v0_2, and nothing else did:

- **Cross sections.** Biagi's level-resolved file (Magboltz 8.97) stays
  the source, md5-pinned, and v0_3 multiplies it by seven trained factors.
  - The elastic momentum transfer gets a factor read log-linearly between
    0.03, 0.2, 0.6, 2, 6 and 20 eV. The two lowest stay at 1, a choice
    made on held-out data (below), so it trains at 0.6 to 20 eV.
  - The four 5s (1s) levels, the ten 5p (2p) levels and the 37 higher
    levels get one factor each.
  - Ionization, which beam experiments measure, is not trained.
  - `xsec.processes()` applies the factors pinned in
    `data/kr_swarm_trained.json`, `xsec.raw_processes()` is the file as
    published, and `xsec.using(f)` runs with other factors.
  - `xsec.UNIT` is v0_2 number for number (`tests/test_v0_3.py`).
- **`boltzmann.py`**: where no EEDF is found from a warm start,
  `Solver.eedf` starts the grid's solve again cold, once, and counts it
  (`Solver.restarts`, logged per table update as "restarts"). The mean
  match has one solution from any start, so no answer changes; a run
  that never failed is v0_2's path exactly. A beam's afterglow on the
  trained set needed it (`CHANGELOG.md`, 2026-10-15).
- **`swarm.py`** also gives the transverse diffusion N D_T and D_T/mu.
  It obeys Einstein's relation in the thermal limit (kTg/e to 3e-5 at
  1e-4 Td).
- **`swarm_mc.py`** (new) holds LoKI-MC's ratios: what a swarm experiment
  measures, over the two-term swarm's own numbers.
  - The bulk drift velocity, which a drift tube times, is up to 23 %
    above the flux one at 300 Td.
  - The SST coefficient, which a Townsend experiment measures, is up to
    9 % below the two-term k_i/W at 300 Td.
  - The two-term transverse diffusion is 46-79 % high from 10 to 300 Td.
  - The two-term flux W matches the Monte Carlo's to 0.3 % from 10 to
    300 Td; at 0.2 Td it is 1.8 % off, which the factors absorb.
  - On the trained set the ratios hold to 1.65 % (30 to 300 Td; the most
    is alpha/N at 50 Td, where the Monte Carlo's counting error is about
    1 %).

```bash
python3 LIGHT_Kr/run_v0_3.py --eedf two_t           # the reference point, seconds
python3 LIGHT_Kr/run_v0_3.py                        # on the kinetic EEDF, minutes
python3 ML/run.py kr-train                          # retrain the factors (the file to pin)
/usr/bin/python3 -m pytest LIGHT_Kr/v0_3/tests -q   # from Kinetic_Modeling
```

## The training (ML/ml/kr_*.py)

- **The data.** 28 krypton datasets, 441 points. They come from Dutton's
  printed tables and from LXCat's LAPLACE database, the latter read with
  the owner's approval of LXCat's terms. Roles:
  - **train:** what the consistency check keeps.
  - **held out:** what it drops, and D_T/mu below 1 Td. Koizumi's D_T/mu
    and Pack's two 195-K drift velocities chose which factors train, so
    they are **validation** sets, not independent tests.
  - **held out by rule:** alpha/N below 50 Td (the Known-bad list).
  - **check:** above 300 Td, Specht's total coefficient, and D_T/mu above
    1 Td.
- **The error model.** No file prints an error, so each dataset carries a
  class value: 3 % W after 1985, 5 % before, 7 % D_T/mu, 10 % alpha/N.
  On top of that, 2 % in E/N enters through the model's own slope.
- **The emulator.** A neural network reproduces the true solver: 2000
  Latin-hypercube designs at 81 nodes, held-out rms 0.10 % (W), 0.13 %
  (D_T/mu), 0.74 % (alpha/N). Every answer is re-checked on the solver.
- **Consistency.** Fit all candidates and drop the worst while its chi2
  per point at the joint fit exceeds 3. That drops 8 of 16 candidates
  (English 1953, Kucukarpaci 1981, Pack 1992 at 300 K, Dall'Armi's W,
  Pack 1992 at 195 K, Pack 1962 at 300 K, Pack 1962 at 195 K, Bowe 1960).
  - Train: Hunter 1988, Nakamura 1990 and Pack 1962 at 368 K (W), and
    Bhattacharya, Dall'Armi, Jacques, Kruithof and Heylen at 50-300 Td
    (alpha/N).
  - Threshold 2 gives the same set. At 4, Pack 1962 at 195 K and Bowe
    1960 (9 % below Hunter where both exist) join, and the 6-eV factor
    moves to 1.31.
- **Which factors train.** Five variants were fitted and scored on the
  true solver (chi2 per point):

  | factors trained | train sets | Koizumi D_T/mu | Pack 1962, 195 K | Pack 1992, 195 K |
  |---|---|---|---|---|
  | none (Biagi's file) | 1.75 | 2.52 | 3.07 | 6.12 |
  | the three excitations | 1.16 | 2.52 | 3.07 | 6.12 |
  | these + momentum transfer at 2-20 eV | 0.66 | 2.43 | 2.89 | 5.29 |
  | these + 0.6 eV (**v0_3**) | 0.64 | 0.84 | 2.77 | 5.27 |
  | all nine | 0.49 | 8.06 | 5.45 | 8.70 |

  All nine fit the train sets best, but they trade the 0.03-eV knot
  (1.13) against the Ramsauer minimum (0.88), which the held-out sets
  reject. v0_3's variant fits those best of the five.
- **The factors** (least squares; the MCMC median and the neural
  posterior estimate agree, `ML/README.md`):
  - momentum transfer: 1 at 0.03 and 0.2 eV (not trained), 1.16 at
    0.6 eV, 0.93 at 2 eV, 1.04 at 6 eV, 1.23 at 20 eV and up;
  - excitation: 0.89 (5s), 0.92 (5p), 0.75 (higher levels).
  - The data barely pin the 5p factor (the network's R^2 0.15).
- **Scores**, the model as each experiment measures it, Biagi's file ->
  trained:
  - train W: 4.8 -> 2.8 % rms (2.1 % with all nine factors free).
  - train alpha/N (50-300 Td): 14.0 -> 9.9 % rms (mean -11 -> 0 %).
  - validation, Koizumi's D_T/mu: 15.8 -> 7.1 % rms (mean +11 -> -2 %).
  - held out, the eight drift velocities the check drops: five closer
    (the two 195-K sets, validation, among them) and three further:
    Pack 1992 at 300 K (chi2 14.2 -> 17.7 per point), Dall'Armi
    (10.3 -> 13.4) and Kucukarpaci (39.6 -> 48.4).
  - untouched tests, which neither the training nor the choice of factors
    saw:
    - alpha/N below 50 Td: mean -44 -> -36 %. Those datasets part by up to
      x4.0 from one another where they overlap.
    - the beam's partition against Saito et al.: W +16 -> +10 % and
      N_ex/N_i +44 -> +24 %.

## At the reference point

LIGHT's drive in pure krypton: 1 atm at 300 K, 300 kW/cm3 for 1 us into
the electrons, a 1e15 cm-3 seed, 2 us in all.

| Fate at 2 us, % of the input | kinetic EEDF (default) | two-temperature | v0_2 kinetic | v0_2 two-temperature |
|---|---|---|---|---|
| 147 nm (Kr2*) | 14.05 | 15.78 | 14.20 | 16.00 |
| near IR (5p -> 5s) | 0.64 | 0.63 | 0.64 | 0.63 |
| heat, reactions and relaxation | 83.09 | 81.37 | 83.11 | 81.40 |
| heat, elastic collisions | 2.75 | 2.74 | 2.58 | 2.49 |
| gas temperature at 2 us | 808 K | 798 K | 807 K | 797 K |

The trained excitation is smaller, so a little less of the input reaches
Kr* and the 147-nm band. The answer is still decided by the 5s <-> 5p
cycle on Jung et al.'s measured cross sections. With the cycle off half
the input is 147-nm light (50.6 % kinetic, 51.8 % two-temperature). At
the band's ends the share runs from 18.6 % (x0.5) to 12.2 % (x1.4) on the
kinetic EEDF, and Kannari's fit in K76's place gives 15.3 %.

Propagated uncertainty of the kinetic 147-nm share: +-11.8 points in
quadrature (+-84 %). The trained cross sections' posterior (12 draws on
the kinetic EEDF, 40 on the two-temperature closure) adds +-0.04 points
(+-0.16 on the two-temperature closure). Kr* + 2Kr is still the largest
term, 75 % of the variance.

## Checks

- `/usr/bin/python3 -m pytest -p no:cacheprovider LIGHT_Kr/v0_3/tests`
  from Kinetic_Modeling: 387 passed, 1 skipped without `test_docs.py`
  (19 min), and `test_docs.py` with LIGHT_KINETIC_SLOW=1 7 passed
  (8 min). The gate, `tests/test_v0_3.py`, holds the following.
  - Everything but the factors and the solver's restart is v0_2's.
  - With unit factors v0_3 is v0_2 number for number.
  - The trained file is the pinned one, inside the training's box.
  - The swarm's D_T/mu obeys Einstein's relation in the thermal limit.
  - A failed warm start restarts cold, once, with the cold answer.
- `ML/tests`: 53 passed, 4 skipped (the slow ones), 1 xfailed.
  `ML_KR_SLOW=1 ML/tests/test_kr_slow.py`: 4 passed (the emulator against
  the true solver, the train set, and the pinned factors' chi2 on a fresh
  true-solver run).
- The overview's generator raises when a claim its text makes fails, for
  example "training improves the train sets", "v0_3's variant fits the
  validation sets best" or "the beam's partition moves toward Saito's".
  It also scores the five variants live and checks the training record
  against them to 1 %.
- LoKI-MC on the trained set: the bulk-W and SST ratios lie within
  1.65 % of Biagi's file's, 30 to 300 Td.

## Files

```
LIGHT_Kr/
  run_v0_3.py         launcher (puts LIGHT_Kr/ on sys.path)
  v0_3/
    xsec.py           LXCat parser, the trained factors (Factors, TRAINED_FILE, using), Maxwellian tables, the two-temperature tail
    swarm.py          the solver as a swarm: W, alpha/N, N D_T, D_T/mu, with the Hornbeck-Molnar share
    swarm_mc.py       LoKI-MC's ratios: the bulk W, the SST coefficient and D_T/mu over the two-term's
    species.py, network.py, model.py, report.py, run.py, boltzmann.py, jung2006.py, ib.py   (v0_2's; boltzmann.py and model.py with the cold restart)
    data/Kr_Biagi897_LXCat.txt      Biagi's file (md5-pinned)
    data/kr_swarm_trained.json      the trained factors and the training's record (md5-pinned; ML/run.py kr-train writes it)
    data/kr_tail_ratio.txt          the two-temperature tail table (v0_1's)
    tests/            test_v0_3.py (the gate) and v0_2's tests carried over
  docs/
    accuracy_overview_v0_3.tex/.pdf, system_map_v0_3.tex/.pdf
    make_overview_v0_3.py, make_system_map_v0_3.py, kr_docs_v0_3.py
ML/
  ml/kr_light.py, kr_data.py, kr_swarm.py, kr_nn.py, kr_train.py
  data/swarm/Kr_*_LAPLACE_LXCat.txt
```
