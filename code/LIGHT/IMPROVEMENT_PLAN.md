# LIGHT: what to build next

> **Status, 2026-09-27 (later the same day).** Cooper chose the basic path.
> Keep the basic processes, pin each to a KrF measurement, and perfect the
> Kr/F2 gas first. That is now LIGHT's default `basic` network (see
> README.md, "The basic network"). It settles decision 1 below in favour of
> the measured B/C structure. The rest of this document is a reference menu,
> not a to-do list.

Survey of Kinetic_Modeling, 2026-09-27. Cooper asked which of the many
reactions in the repo, and which missing pieces such as recombination, would
make an improved LIGHT most accurate against experimental data. This plan
answers from two kinds of evidence. Three survey notes walk the rest of the
repo, and measurements made with LIGHT itself back every number below that
comes from LIGHT. Both live in `docs/survey_2026-09-27/`, and the
README there lists what was re-checked first-hand.

## The short answer

1. **The pool is large.** The 15 rate tables printed in
   `docs/Recent_KrF_Papers` hold 893 rows, which come to about 429 distinct
   reactions. With REBUILD, phantom11 and Nathaniel's model added there are
   about 460, and about 650 across every gas family. About 310 fall inside
   the Kr/Ar/Ne/F family, and 130 of the literature ones there are not in
   REBUILD. LIGHT carries 38 today.

2. **Pick reactions by importance and by source, never by fit.** Choosing the
   subset that best matches the data would be calibration, which the repo's
   predictive-only policy forbids. With dozens of uncertain rows it would also
   match almost anything. Instead, keep every reaction that measurably moves
   an output in the operating range, take each rate from its best source, and
   then score against experiments blind.

3. **The largest lever is the KrF B/C structure, not a missing reaction.** Xu
   1993 measured the coupled KrF(B,C) lifetime as 25 ns. That forces a B/C
   equilibrium constant of 4 to 5, where LIGHT and REBUILD use Kannari's 1.47.
   The Kr quench constants both codes use were extracted assuming 4 to 5, so
   the present combination is inconsistent. It misses the only experimental
   data LIGHT can score today, the Xu and Gadomski photolysis decays in pure
   Kr/F2, by ×1.5 to ×2.9. The measured structure reproduces them to ×1.00.
   The choice halves the 248-nm output, 5.9 % against 2.9 % at 1 atm.

4. **The drive decides whether recombination matters.** An electron beam
   puts 52 % of its energy into ionization, 32 % into excitation and 17 %
   into slow electrons. Power into a Maxwellian electron gas, as in v0, puts
   92 % into excitation and 10 % into ionization. Under a beam the ionic
   channel forms about half of KrF(B). Switching off attachment, ion
   conversion and recombination cuts the 248-nm output ×9 to ×14. Under v0's
   drive the same switch moves it by less than one point. Every KrF gain
   dataset in the tree is e-beam pumped.

5. **One omitted electron process matters under v0's drive.** Electrons lift
   Kr*(5s) to Kr**(5p) at Kannari's 9.0e-7 Te^0.5 exp(−1.51/Te). At
   n_e = 1e15 that runs about three times faster than the harpoon. REBUILD
   left it out because it is negligible at beam densities. In LIGHT it cuts
   the 248-nm output by 37 %.

6. **n_e has to follow charge balance.** Once electrons attach, holding n_e
   fixed makes F⁻ out of nothing. The ion store went to −25 % of the
   deposited energy at 2 atm in the prototype.

7. **Gary's Version 0 list overstates the 248-nm output ×3 to ×6** against
   the complete Kr/F2 network (16.8 % against 5.9 % or 2.9 %, depending on
   the B/C choice). Its partition is a statement about those six reactions,
   not yet about a KrF gas.

## How to pick from the pool

Three steps, applied to every candidate reaction.

1. **Importance.** Start from the most complete network available and remove
   one reaction at a time across the operating range (pressure, power, n_e
   and drive). A reaction that moves no output anywhere can be dropped. Any
   reaction that moves the 248-nm share, KrF(B) formation or the heat share
   by more than about one point stays. Parallel paths can hide each other,
   so recombination is also tested as a group. `loo_screen.py` and
   `recomb_isolate.py` do this now.
2. **Source.** For each reaction that stays, take the best printed value: a
   direct measurement before a model fit, and a fit before an estimate. Where
   sources disagree, carry the spread as a bracket and report the spread in
   the outputs, as REBUILD does.
3. **Scoring.** Compare against experiments without changing a rate in
   response. Start with data that isolate kinetics (decays and quench
   rates), then gain, then device energies. Device energies amplify every
   error by about one over the lasing margin.

## What matters, measured in LIGHT

### One row at a time

Largest change in the 248-nm share of E_dep, in percentage points, when one
row is removed from v0 plus every optional group. The power-drive column is
the largest change over five points (0.1 to 2 atm, n_e 1e14 to 1e16). The
e-beam column is the prototype beam drive at 1 atm.

| Row | Reaction | Source status | Power drive | E-beam |
|---|---|---|---|---|
| CB_Kr, BC_Kr | KrF B↔C mixing by Kr | estimated class value | 6.8, 6.1 | 6.3, 5.8 |
| K67 | KrF(B) + e quench | repo value, band 3–6e-8 | 5.0 | 0.6 |
| K33 | KrF(B) + Kr quench | measured (Xu) | 4.4 | 3.7 |
| G3 | Kr* + F2 harpoon | measured | 4.1 | 1.4 |
| G5 | KrF(B) + 2Kr trimer | measured, Mandl hints ×0.5 | 3.0 | 2.4 |
| KR2_RAD | Kr2* radiative decay | sources differ ×24 | 2.1 | 0.9 |
| K2F_F | Kr2* + F → KrF | literature | 1.4 | 0.7 |
| ATT | e + F2 attachment | SIGLO cross sections | 1.1 | 1.2 |
| AK8_Kr | Kr* + 2Kr → Kr2* | estimated | 0.8 | 0.5 |
| KSTEP | e + Kr* → Kr⁺ | estimated | 0.7 | 0.1 |
| K15, R11p | Kr⁺ → Kr2⁺, Kr2⁺ + e | literature | 0.2, 0.2 | 0.4, 0.4 |

The recombination rows barely move the light one at a time, because the
parallel paths take over. They move the heat share by 15 to 19 points.

### Recombination as a group

The `ions` and `f2e` groups (ion conversion, dissociative recombination,
attachment, ion-ion recombination) switched off and on, every other group
on. Numbers are percent of E_dep.

| Drive | p (atm) | KrF(B) formed, off → on | 248 nm, off → on | Ions left at the end, off |
|---|---|---|---|---|
| e-beam | 0.5 | 7.6 → 17.4 | 0.48 → 5.02 | 63.6 |
| e-beam | 1.0 | 10.7 → 31.2 | 0.47 → 6.65 | 54.8 |
| e-beam | 2.0 | 11.5 → 33.0 | 0.28 → 2.48 | 52.1 |
| power | 0.5 | 16.4 → 17.8 | 4.98 → 5.78 | 52.8 |
| power | 1.0 | 28.3 → 27.8 | 5.08 → 5.94 | 24.1 |
| power | 2.0 | 32.5 → 23.1 | 2.14 → 2.16 | 11.2 |

Under the beam, the harpoon forms 39 to 42 % of KrF(B), ion-ion
recombination 48 to 53 %, and the Kr2* route the rest. Without electron
losses, n_e grows until electron quenching of KrF(B) takes over, which is why
the light collapses.

### The B/C structure against the only data LIGHT can score

KrF decay constants (s⁻¹) in pure Kr with 1 Torr F2 at 300 K, no power, no
electrons. The measured values at 50 and 100 Torr are printed on Xu 1993
p.2594. The 500-Torr value was read from Gadomski's fitted curve.

| Kr (Torr) | Measured | v0 | + quench | + quench + bc | + measured B/C structure |
|---|---|---|---|---|---|
| 50 | 5.0e7 | ×2.89 | ×3.21 | ×1.61 | ×1.01 |
| 100 | 6.2e7 | ×2.41 | ×2.77 | ×1.50 | ×1.01 |
| 500 | 2.63e8 | ×1.13 | ×1.42 | ×1.12 | ×1.00 |

This is a consistency test rather than a blind prediction, since Xu and
Gadomski extracted their constants from these same decays. It still settles
that LIGHT's present mix of sources is inconsistent. With every group on,
the measured structure takes the 248-nm share from 5.94 to 2.93 % at 1 atm
(power drive) and from 6.65 to 3.29 % (e-beam).

### The Kr** (5p) coupling

An explicit Kr** with Kannari's electron up-pumping ([76]), its reverse,
stepwise ionization from 5p ([79]), the 5p harpoon, radiative decay and p→s
relaxation. All other groups on, 1 atm.

| Drive, n_e | 248 nm, without → with | KrF(B) formed | Heat |
|---|---|---|---|
| power, 1e15 | 5.94 → 3.75 | 27.8 → 17.8 | 71.0 → 84.6 |
| power, 1e16 | 2.46 → 1.59 | 27.4 → 17.8 | 74.6 → 88.5 |
| e-beam | 6.65 → 6.54 | 31.2 → 30.2 | 74.1 → 76.1 |

The Kr* ↔ Kr** cycle turns electron power into heat and cools the electrons,
which starves direct ionization as well.

## What LIGHT can be scored against

From `notes_B_experimental_data.md`, ranked by how well a 0-D model can use
them.

| Rank | Data | What it measures | What LIGHT needs |
|---|---|---|---|
| 1 | Xu 1993 and Gadomski 1992, pure Kr/F2 | KrF(B,C) and Kr2F* decays against Kr and F2 | Seed mode (exists), the B/C decision |
| 2 | Xu 1993, He/Ne/Ar series | Buffer quench rates | Ar and Ne (v1) |
| 3 | Trainor & Jacob 1980 | KrF emission in e-beam Kr/F2 99.5/0.5 at 190 Torr, about 70 kW/cm³ | The paper (not in the tree), the beam drive |
| 4 | Klimek 1981 via Kannari Tables II/III | 18 gain and 18 absorption values, g0/α and pressure exponents | v1, beam drive, gain and absorption outputs |
| 5 | Aurora and Tacoma | Gain at printed deposition | v1, beam drive |
| 6 | Johnson & Hunter 1980 | g0/α for four devices | v1, beam drive |
| 7 | Mandl 1985 | Efficiency at 4, 10 and 99.6 % Kr | v1, beam drive, extraction |
| 8 | Device energies (Kannari, Electra, Rhodes) | Laser output | A cavity. Errors grow as one over the margin |

The gain-to-absorption ratio is the safest device-level target, because it
does not depend on the unknown pump calibration.

## The plan

### Step 1. A consistent, complete Kr/F2 box (v0.2)

1. **Settle the B/C structure** (decision 1 below). Measured K_eq with
   Hakuta's k_BC = 1.2e-10 is the only structure consistent with the quench
   constants LIGHT already uses, and the only one that fits the data LIGHT
   can score.
2. **Promote the six optional groups** into the default network. All are
   drift-gated REBUILD copies, and quench and B/C mixing each move the
   248-nm share by about 6 points.
3. **Add Kr** (5p) with its electron coupling** (Kannari [76] and [79]), the
   5p harpoon (REBUILD K27), radiative decay and p→s relaxation. Kr**
   replaces the prompt cascade in C1.
4. **Replace estimated rows where a printed value exists**, in REBUILD first
   so the drift gate keeps the two codes in step:
   - KSTEP: Kannari prints 7.00(−08) Te^0.68 exp(−4.04/Te). The estimated
     template runs ×2.2 above it at 1.3 eV.
   - BC_Kr: 1.2e-10 (Xu p.2591, from Hakuta).
   - AK8_Kr: the estimate matches Tracy & Oskam's 2.59e-32 (via Rhodes
     Table 4.7), so this is a change of tier with a source, not of value.
   - KPS_Kr: 3.0e-11 appears second-hand in the Janacaitis list. Keep the
     estimate as a bracket until the primary is found.
5. **Add KrF* + F** at Kannari's 7.8e-10 ("analogous to F2", not measured),
   bracketed down to 1e-12 (decision 4). It outruns KrF(B) radiation once
   the F2 has burnt out, which the notes' drive reaches in 830 ns.
6. **Let n_e follow charge balance by default**, with injection as a source
   term. Stabilise the electron energy at low n_e, since the prototype's
   Jacobian went non-finite when attachment drove n_e to its floor.
7. **Add gain and absorption outputs**, g0(t) = σ_B n_B with σ_B = 2.6e-16 cm²
   and α(t) species by species from REBUILD's 248-nm absorber list. They are
   what the experiments measure.
8. **Add the Xu/Gadomski decays as a test.** `decay_run.py` already runs it.

**Done when:** the decay test passes, the drift gate and ledger closure hold,
and the leave-one-out screen shows no row of the complete Kr/F2 set that
moves an output by more than a point while resting on an estimate that has a
printed alternative.

### Step 2. The e-beam drive (v0.3)

1. **Port REBUILD's Spencer–Fano yields** for the fill, giving Kr⁺, Kr* and
   Kr** per ion pair plus the slow-electron energy handed to the electron
   gas. Refuse the fixed-W fallback, which moved REBUILD's anchors ×2 to ×7.
2. **Complete the recombination set** now that it matters:
   - route part of Kr2⁺ recombination to Kr** (Kannari 1.9(−07))
   - add three-body electron–ion recombination for atomic Kr⁺.

   At 0.1 atm, Kr⁺ takes about 700 ns to convert to Kr2⁺. The tree has no Kr
   value for three-body recombination, only Ar analogs (Kossyi 7.2e-27
   Te^−4.5 cm⁶/s), so the Kr primary is needed first.
3. **Score against Trainor & Jacob 1980** once the paper is in hand. It is
   the one pure-Kr e-beam dataset with LIGHT-native conditions.

### Step 3. Ar and Ne (v1, the notes' Version 1)

1. **Bring in the buffer chemistry.** Add species and about 130 rows from
   REBUILD's Ar and Ne sets, plus md5-pinned copies of the Biagi Ar and Ne
   cross sections already in `LoKI-MC/Code/Input`.
2. **Fix these rows before relying on them:**
   - NK6 copies Chen's misprinted R23 and runs ×9.3 high in Ne fills.
   - NK3 is estimated at 1e-11. Nathaniel's rate table cites a measured
     bound of ≤1e-14 (Johnsen 1978).
   - NK1 is estimated at 6.5e-11. The same table cites a measurement of
     9.0e-11 (Neynaber & Magnuson).
   - A2F_KR2 and AKF_AR are estimates that carry about half of REBUILD's
     Klimek gain.
3. **Score** against Klimek's g0/α and pressure exponents, Xu's buffer quench
   series, Aurora and Tacoma, Johnson & Hunter, and the Mandl 4/10/99.6 % Kr
   ratios.

### Step 4. Electrons, afterglow and extraction (later)

1. **Replace the Maxwellian for field and EM drives** with a swarm table for
   LIGHT's own mixture from LoKI-MC, or the two-term Boltzmann solver in
   `~/coding/PIC/Neon/scripts/boltzmann2t.py`. The repo measured Maxwellian
   threshold rates ×1e3 to ×1e6 off at low E/N.
2. **Complete the afterglow chemistry:**
   - F + F + M, printed at 1.0e-33 and negligible within a pulse
   - F⁻ detachment, for which no primary source exists in the tree
   - the afterglow tests the repo planned and never built (EXPERIMENTS.md
     B7 and B8).
3. **Add a cavity only if device energies are wanted** (port REBUILD's
   `laser.py`).

## Decisions for Cooper and Gary

1. **The B/C equilibrium constant**, Kannari's 1.47 or the measured 4 to 5. It
   halves the 248-nm output and decides whether LIGHT matches the decay data.
   In REBUILD it would lower the Klimek gain score from ×0.64 to about ×0.47,
   and that anchor's pump is borrowed from a different device.
2. **The drive.** The notes say the deposition route does not matter, but it
   moves ionization from 10 % to 52 % of the deposited energy. I recommend
   carrying both, with the beam used for validation.
3. **The trimer rate K34.** Xu measured (6.0 ± 2.0)e-31. The Mandl ratio
   prefers about half. Carry both until pure-Kr data are in hand.
4. **KrF* + F.** Printed values span 1e-12 to 7.8e-10 and none is a
   measurement.
5. **τ_sp**, Kannari's 6.99 ns (Xu uses 7.0) or the notes' 4 ns.
6. **Order of changes.** Recommend changing REBUILD first and letting LIGHT
   follow through the drift gate, rather than letting LIGHT carry values
   REBUILD does not.
7. **Primary sources not in the tree:** Trainor & Jacob 1980, Mandl 1985,
   Klimek 1981 and Suda 1984. They hold the only pure-Kr device data.

## Found outside LIGHT

These were found during the survey and not changed. They belong to REBUILD
or to the frozen legacy code.

- **NK6** (Kr* + Kr + Ne → Kr2* + Ne) is set to 2.5e-43 m⁶/s as an "analog of
  Chen R23". REBUILD's own errata correct R23 to 1.1e-44, so NK6 runs ×9.3
  high and touches every Ne-buffered REBUILD result.
- **BC_Kr and CB_Kr** name Gadomski 1992 as their primary, but Gadomski
  assumed k_BC. The measured value (Xu p.2591, from Hakuta) is 1.2e-10,
  ×0.24 of the row. The row's note says the class value likely
  underestimates.
- **VREL_Kr** names Kvaran 1988 as its primary. Kvaran avoided a Kr buffer
  on purpose and has no Kr data.
- **KSTEP**'s estimated template runs ×2.2 above Kannari's printed rate at
  1.3 eV.
- **Kr* + e → Kr\*\*** is left out of REBUILD as negligible at beam
  densities, which holds there. Any REBUILD discharge or microwave study at
  n_e ≥ 1e15 needs it.
- **The phantom11 Stevefelt row** (`factored-pyfiles/phantom11/config.py`
  657–669) gives 1.3e-42 cm⁶/s at 0.3 eV where its comment says 2.4e-22,
  and its third body should be an electron. Do not port it.
- **Reported by the survey but not re-checked:**
  - `SPECK/EXPERIMENTS.md` line 46 labels an ArF model-output benchmark as a
    blind KrF anchor.
  - `excimer_anchors.py` cites p.4314 for a number printed on p.4315.
  - phantom's r120 codes Kr2⁺ + F2 where J&H print Kr2* + F2.

## Traps to carry forward

- Never choose the reaction subset by fit to the data.
- Kannari prints unreduced quench constants. Do not reduce them.
- B/C brackets move as a pair, or detailed balance breaks.
- A silent fixed-W fallback moved REBUILD's anchors ×2 to ×7.
- The parity harness is blind to Ar-weighted terms, which matters for v1.
- Check every phantom row against its primary before copying it.
