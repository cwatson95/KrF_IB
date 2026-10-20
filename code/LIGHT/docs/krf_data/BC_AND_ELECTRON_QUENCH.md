# KrF(B), KrF(C) and the electron quench (2026-09-27)

Owner direction, 2026-09-27: get the KrF interactions right before any
drive or device module. Focus on the B and C states and on the electron
quench rate.

**checked** means re-read first-hand in the source. **agent** means a
search agent quoted or digitized it and nobody has re-read it.

## Latest: Yu & Setser in full, and the channels inferred from the data

**Yu & Setser, J. Phys. Chem. 94, 2934 (1990), read in full (checked, pp. 2940-2941,
`sources/Yu_Setser_JPC94_2934_1990.pdf`).**
- **Method.**
  - They measured the equilibrium C/B emission ratio at high pressure (up to 5 atm of Ar, N2 or CF4), at 300 and 230 K.
  - They removed the B→A emission under the C-A band (branching 0.035).
  - They take the partition-function ratio as "nearly unity".
- **Result.**
  - "I_r(300 K) = 0.16 ± 0.02, K_eq(KrF) = 1.7 ± 0.2", using the calculated τ_B 7.0 and τ_C 75.2 ns.
  - ΔE_BC = −110 ± 20 cm⁻¹ from the absolute ratio, and −180 ± 30 cm⁻¹ from the 230/300 K change. Their "consensus value" is −120 ± 50 cm⁻¹.
  - Xu 1993's "2-3" was a loose paraphrase.
- **On mixing.**
  - "more than 1000 Torr of Ar is needed to obtain equilibrium B and C concentrations"
  - B↔C mixing by Ar ≈ 1e-11
  - a coupled lifetime of "15 ± 5 ns", against Gadomski's later 25 ns
  - This is slow mixing, from the Kansas State lab itself.

**The channels inferred from the four kinetic datasets** (now `python3 ML/run.py consistency` and `fit`).
- **The data:** Eden, Kansas State, Hakuta, and Yu & Setser's I_r(300 K). The ratio is modelled as their eq. 18, I_r = 1.035 (τ_B/τ_C) K_eq, which returns 0.160 with their own constants.
- **The channels (10):** k_BC, K_eq, τ_C, τ_B, and Kr, 2Kr and F2 quench for B and C separately.
- **Part 1 (consistency).** Least-squares fits with each dataset left out in turn. Only Eden + Kansas State + Yu & Setser fit together within their errors. Hakuta is the odd one, whichever way the leave-outs are cut.

| Fitted to | χ² per point: Eden | KSU | Hakuta | Yu & Setser | k_BC | lasing share, 1 atm |
|---|---|---|---|---|---|---|
| all four | 0.26 | 2.22 | 2.54 | 2.45 | 1.5e-11 (τ_C 30 ns) | 0.58 |
| **Eden + KSU + Yu & Setser** | **0.07** | **0.18** | 14.6 | **0.00** | **1.1e-11** | **0.68** |
| KSU + Hakuta + Yu & Setser | 22.4 | 2.11 | 1.15 | 0.53 | 1.6e-10 | 0.41 |
| Eden + Hakuta + Yu & Setser | 0.19 | 20.0 | 0.56 | 0.00 | 4.2e-11 | 0.56 |
| Eden + KSU + Hakuta | 0.37 | 1.85 | 2.49 | 99.4 | 3.3e-11 | 0.43 |

- **Part 2 (MCMC on Eden + KSU + Yu & Setser).** Median and 68 % interval:
  - k_BC = 9.8e-12 [4.6e-12, 1.7e-11]
  - K_eq = 1.07 [0.84, 1.38]
  - τ_C = 50 [41, 63] ns (Hakuta measured 60 ± 10)
  - τ_B = 7.0 [6.8, 7.2] ns
  - F2 quench: B 4.6e-10 (Eden's B-only 4.8e-10), C 3.0e-10 (Xu's coupled 3.0e-10). The old "F2 conflict" was the two states.
  - 2Kr: B 1.9e-30, C 3.6e-31. The Kr two-body quench is near zero and poorly determined.
  - **Lasing B share at 1 atm, 300 K: 0.69 [0.58, 0.79], against 0.27 in LIGHT now.**

**The neural-network test** (now `python3 ML/run.py nn`, numpy only, amortized neural posterior estimation):
- **Setup.** An MLP with a Gaussian output was trained on 30 000 simulated experiments. Each simulation drew channel values from wide priors, pushed them through LIGHT's two-level model, and added each dataset's noise.
- **Network 1 (Eden + KSU + Yu & Setser).**
  - On held-out simulations: R² 0.94 (k_BC), 0.97 (K_eq), 0.86 (τ_C), 0.90-0.93 (F2 quench). Coverage of its ±1σ is 0.64-0.66 (ideal 0.68).
  - Kr two-body (R² 0.08-0.24) and B-state 2Kr (0.35) are not identifiable, and it says so with wide ranges.
  - On the real data: k_BC 7.5e-12, K_eq 1.15, τ_C 49 ns, share 0.73 [0.61, 0.85]. This is the MCMC's answer.
- **Network 2 (all four datasets).**
  - It returns a confident compromise: k_BC 2.9e-11, τ_C 35 ns, share 0.59 ± 0.04.
  - That compromise misfits Kansas State at χ² 8.6 per point, and nothing in the network's output flags it.
- **What this means.** A network reads constants off consistent data as well as a direct fit. It cannot decide between experiments that disagree, and that decision came from the consistency check in Part 1.

**Accuracy against data the fits never saw** (now `python3 ML/run.py score`). Each set's electron quench is re-derived from Trainor & Jacob with its own lifetime.

| Data | LIGHT now | Gadomski 1992 | Inferred (Eden + KSU + Yu & Setser) |
|---|---|---|---|
| Eden B-only decays (fitted for Inferred) | 21 % rms | 21 % | 2.5 % |
| Kansas State Figs. 1-2 (fitted) | 6.4 % | 5.4 % | 5.1 % |
| Yu & Setser I_r(300 K) (fitted for Inferred) | +277 % | +202 % | 0.1 % |
| **Xu 1993 Figs. 5-7, 2 Torr F2 (held out)** | 10.5 % | 12.7 % | 11.6 % |
| **Hakuta Figs. 2-3 (held out)** | 29 % | 36 % | 30 % |
| **Burnham & Searles τ_B 9.0 ns (held out)** | −22 % | −22 % | −22 % |
| **Yu & Setser 230/300 K ratio (held out)** | −14 % | −14 % | −14 % |
| **Tahoma small-signal gain 6.39 %/cm (device, held out)** | 3.67 (×0.57) | 3.82 (×0.60) | 7.28 (×1.14). Over 40 posterior draws: 7.15 [6.44, 7.89] |
| Tahoma formation efficiency 0.27 | 0.341 (×1.26) | 0.341 | 0.341 |

How to read it:
- **The held-out kinetic data** are matched to 10-12 % by every set. That is inside the data's own reproducibility: the Kansas State lab's two 1992-1993 series differ by up to 25 % at one condition.
- **Misses common to all three sets**, whose cause is in inputs the sets share:
  - Hakuta, by about 1.8× in absolute rate
  - Burnham & Searles' τ_B (Eden measured 6.8 ns)
  - Yu & Setser's temperature dependence (LIGHT's ledger gap is 80 cm⁻¹, their two-temperature value 180 ± 30).
- **Tahoma's gain.** This is the one device measurement, and no fit used it. The inferred set predicts it to +14 %, with the measurement at the edge of the 68 % band. LIGHT now misses it by −43 %.
- **Formation efficiency** is 26 % high for every set. The excess is upstream of B and C, in the prototype beam drive and formation chemistry, and is the next thing to understand.
- **Policy.** The device score is after the fact only, and no constant was moved toward it.

**What adopting Eden + KSU + Yu & Setser would change in LIGHT.**
- **Mixing and lifetimes:** k_BC 1.2e-10 → about 1e-11, K_eq 5.0 → about 1.1, τ_C 60 → about 50 ns.
- **Quench:** the F2, Kr and 2Kr quench constants become state-specific. The B and C rows already exist separately.
- **Result:** the lasing B share at 1 atm goes from 0.27 to about 0.7.
- **What it reverses:** it would undo Hakuta's k_BC, adopted one step earlier.
- **What it leaves unexplained:** Hakuta's own B-X tail at 9 Torr.

This is the owner's call, and LIGHT's constants are unchanged until then.

## In short (updated after Hakuta 1990 and Yu & Setser 1990 were read)

1. **LIGHT now uses the measured B/C values where they exist.**
   - τ_C = 60 ns and k_BC = 1.2e-10 come from Hakuta et al., JJAP 29, 1194 (1990), read in full (checked). They replace Gadomski's assumed 75 ns and 1.0e-10.
   - Gadomski's reduction of the Kansas State decays, redone with those two values, gives K_eq 5.4 ± 10 % and k_Kr 3.6e-12. Both are inside the published 5.0 ± 1.0 and 4.0 ± 1.0, so K_eq and the quench constants stay.
   - K_eq now follows the gas temperature by detailed balance over the gap the energy ledger books (80 cm⁻¹). It is 5.0 at 300 K and 4.4 at 460 K. Before, the kinetics used a fixed 5.0 while the ledger used 80 cm⁻¹.
   - The effect is small. The lasing B share at 1 atm goes from 0.29 to 0.27 at 300 K, and from 0.23 to 0.24 at 460 K. The 248-nm share at the reference point goes from 3.47 % to 3.52 %.
2. **The electron quench stays at 4.5e-8.** Re-read with the new constants, Trainor & Jacob's product gives 4.6e-8.
3. **Two conflicts remain open, and each would raise the lasing share.**
   - **Eden versus Hakuta on B→C mixing.** Eden's B-only decays want k_BC ≤ about 5e-11, and would put the lasing share at 0.4 to 0.8. Hakuta's B-X decay at 9 Torr (his Fig. 1) has a slow tail that "is coincident with the decay time of the C-A fluorescence". With slow mixing that tail would be about 10⁻⁵ of the peak, so Hakuta's data say mixing is fast.
   - **Yu & Setser versus Kansas State on K_eq.** Superseded by the section above: Yu & Setser print K_eq 1.7 ± 0.2 (with τ_C 75.2), and the joint inference gives about 1.1 (with τ_C about 50).
4. **Hakuta and Kansas State disagree on absolute decay rates**, by about 1.8× at 20–50 Torr Kr with 1 Torr F2. No single two-level constant set fits both (Hakuta 29 % rms with LIGHT's constants). A B→A contribution to Hakuta's 280-nm signal (Julienne & Krauss 1979) makes it worse, not better.

## What LIGHT uses and where it comes from

| Constant | LIGHT | Where it comes from | What was measured |
|---|---|---|---|
| τ_B | 6.99 ns | Hay & Dunning's calculation, via Kannari 1985 | 6.8 ± 0.2 ns, Eden 1978, cold B (checked). 9.0 ± 0.6, 9, 9.5 ± 0.6 ns, from the late decay or hot B (agent) |
| τ_C | **60 ns** | **Hakuta 1990 Table I (checked)** | the one measurement. Gadomski assumed Hay & Dunning's 75 ns |
| k_BC (B→C by Kr) | **1.2e-10** | **Hakuta 1990 Table I (checked)** | Gadomski assumed 1.0e-10. Eden's B-only data allow ≤ about 5e-11 |
| K_eq = [C]/[B] | 5.0 at 300 K, **K_eq(Tg) by detailed balance** | Gadomski's fit. Redone with Hakuta's τ_C and k_BC it gives 5.4 ± 10 % | Yu & Setser "2-3" (via Xu, checked) and ΔE −120 ± 50 cm⁻¹ (abstract, checked). Xu prints "4-5" and "3-4". Kannari 1985 and Hakuta used 1.47 (a gap of 80 cm⁻¹) |
| B-C gap in the ledger | 80 cm⁻¹ | REBUILD and Kannari 1985 (Julienne & Krauss) | Yu & Setser 120 ± 50 cm⁻¹, from ratios at 300 and 230 K |
| Kr, 2Kr, F2 quench | 4.0e-12, 6.0e-31, 3.0e-10 | fitted by Gadomski and Xu with B and C set equal | nobody measured B and C separately (Gadomski p.155, checked). Hakuta: Kr 7 ± 3e-12, F2 2.5 ± 1e-10 |
| k_e | 4.5e-8 | Kushner, Hanson & Schneider 1989, a one-state reading of Trainor & Jacob (checked) | two-level re-reading 4.6e-8 (below) |

## The electron quench

Only two direct measurements exist (agent's search, first-hand where marked).

1. **Trainor & Jacob.**
   - Duzy et al., AVCO final report 1980, pp. 19–23 (checked). Journal version APL 37, 675 (1980).
   - Pure Kr/F2, e-beam, 190 and 95 Torr Kr, F2 0.5–4 Torr. Steady B-X sidelight measured against [F2], with n_e = S/(k_att[F2]).
   - The measured product is k_e·τ/k_att = 2.9e-7 s (p. 21). With τ = 6.5 ns and k_att = 4.5e-9 this gives 2.0 ± 0.2e-7.
2. **Nishioka et al.**
   - Rev. Laser Eng. 16, 571 (1988): 4.5e-7 (checked).
   - Measured from the saturation intensity against pump and pressure in Ar/Kr/F2 = 94/6/0.284, using a one-state formula with τ = 6.5 ns and k_att = 7.5e-9.
   - The gas is mostly Ar, outside LIGHT's.

Kushner, Hanson & Schneider, APL 55, 2482 (1989) re-read Trainor & Jacob for KrF(B) alone (checked). Their changes:
- k_a of 1–2e-9, from their computed EEDF at T_e = 1.58 eV
- Kr₂⁺ recombination added as an electron loss
- τ = 5.5 ns (6.5 ns less the Kr quench and the trimer).

That gives 3.7–6.5e-8, recommended as "3-6×10⁻⁸". They never mention KrF(C). Hazi's theory gives 3e-8 for B.

In a two-level model the τ in the product is the Stern–Volmer lifetime of B-X emission. It depends on the B/C constants and on where KrF is formed. `electron_quench.py` re-reads the product with LIGHT's lifetime (8.3 ns at 190 Torr), LIGHT's own F2 attachment rate, and KHS's recombination factor (0.62):

| Reading | τ (ns) | k_att | k_e (cm³/s) |
|---|---|---|---|
| Trainor & Jacob as printed | 6.5 | 4.5e-9 | 2.0e-7 |
| KHS, k_a 2e-9 | 5.5 | 2.0e-9 | 6.5e-8 |
| LIGHT two-level, KrF formed in B | 8.3 | 2.1e-9 (T_e 1.58 eV) | **4.6e-8** |
| LIGHT two-level, formed half in B, half in C | 10.5 | 2.1e-9 | 3.6e-8 |
| same as the bold row, k_BC 1e-11 | 4.7 | 2.1e-9 | 8.0e-8 |
| same as the bold row, T_e 0.67 eV (a 1-eV mean) | 8.3 | 5.2e-9 | 1.1e-7 |

The remaining uncertainties in the measurement itself (agent):
- **The F2 subtraction.** Trainor & Jacob's "AB[F2]" is never printed. A 10 % change in it moves the result by −22 % or +25 %, per the agent's reading of their Fig. 9.
- **The ionization rate S.** Their axis implies about 60 kW/cm³, where a straight-path estimate gives 12–20.
- **The EEDF.**

The test `test_electron_quench_is_trainor_jacob_under_this_b_c_model` holds LIGHT's QE row to this two-level reading within 10 %.

## B→C mixing

**Eden, Waynant, Searles & Burnham, APL 32, 733 (1978) (checked).**
- Photolysis of KrF2 at 193 nm makes KrF(B) alone. Their He test says it is in low vibrational levels.
- The pump is about 3 ns FWHM with a 2.5-ns tail, and it was deconvolved. A 249-nm filter isolates B-X.
- Their signals "were observed to decay exponentially over several e-foldings".
- Their printed fit is Σ = 1/(6.8 ns) + 8.6e-12[Kr] + 9.7e-31[Kr]². They themselves attribute the Kr term to "mixing collisions between the B and C states".

**Kansas State (agent readings, Fig. 2 spot-checked by eye).**
- Gadomski et al. 1992 and Xu et al. 1993. A two-photon Kr + F2 reaction makes hot B and C.
- The decay constant is taken from the latest part of the waveform.
- Their own figures scatter by up to 30 % at one condition: 4.1e7 s⁻¹ in Fig. 2a against about 5.4e7 in Fig. 1, at 50 Torr Kr and 1 Torr F2. Xu prints 5.0e7.

`bc_check.py` compares both with LIGHT's closed-form two-level model:
- Eden's Σ against the model's B(t) after a B-only start, fitted with one exponential down to 5 %
- Kansas State's points against the slow eigenvalue.

Part 2 frees every two-level constant (B and C quench separately, τ_B, τ_C, K_eq) and profiles over k_BC. Deviations are rms.

| k_BC | K_eq (fitted) | τ_C (fitted) | Eden Kr | Eden F2 | KSU Fig. 2 | KSU Fig. 1 | B share, 1 atm |
|---|---|---|---|---|---|---|---|
| **LIGHT: 1.0e-10** | 5 | 75 ns | **25 %** | 10 % | 3.5 % | 7.4 % | **0.29** |
| 3e-12 | 0.38 | 46 ns | 1.5 % | 3.2 % | 4.2 % | 6.3 % | 0.87 |
| 1e-11 | 1.00 | 50 ns | 1.7 % | 3.2 % | 4.2 % | 6.1 % | 0.68 |
| 3e-11 | 2.21 | 60 ns | 3.4 % | 4.0 % | 4.3 % | 5.9 % | 0.47 |
| 5e-11 | 3.18 | 60 ns | 7.5 % | 4.2 % | 4.5 % | 5.7 % | 0.38 |
| 1e-10 | 4.56 | 49 ns | 18 % | 4.0 % | 7.0 % | 5.9 % | 0.29 |

How to read the table:
- Eden quotes ≤10 % error per point, so his data allow k_BC up to about 5e-11 and exclude 1e-10 with any other constants.
- The same holds when Eden's measurement is modelled differently (raw pump-convolved trace, or a late fit window). There 1e-10 misses by 18–44 %.
- The Kansas State data fit equally well anywhere on the valley. With k_BC 1e-10 the free fit returns K_eq 4.6, Gadomski's own answer.
- The B share is the steady [B]/([B]+[C]) at 1 atm, 0.5 % F2 and 300 K, for KrF formed in B. The gain is proportional to it.

**Hakuta, Fujino, Nakayama & Takuma, JJAP 29, 1194 (1990), full text (checked).**
- **Method.** 248-nm photoassociation of Kr + F makes KrF(B) in low v. The 10-ns laser pulse sets the time resolution (p.1194). B-X is monitored at 245 nm and C-A at 280 nm.
- **Fig. 1.** At 9 Torr Kr and 1 Torr F2 the B-X decay is double-exponential: fast for 30 ns, then a tail with the C-A decay time (25 ns). This is direct evidence of fast B↔C coupling.
- **Figs. 2-3.** They fit only the slow C-A decay, against F2 at 9 Torr Kr and against Kr at 1-70 Torr with 1 Torr F2. It rises steeply and saturates near 20-30 Torr.
- **Held fixed in the fit:** τ_B 6.8 ns, 2Kr 2.9e-31 and a gap of 80 cm⁻¹ with C below B (K_eq 1.47), with equal B and C quench.
- **Table I:** k_BC 1.2 ± 0.3e-10 (Kr), τ_C 60 ± 10 ns, F2 2.5 ± 1e-10, Kr 7 ± 3e-12.
- **Conflicts.**
  - With Eden, on the mixing rate.
  - With Kansas State, on absolute rates at the same conditions. Hakuta's C-A decay is about 6e7 s⁻¹ at 20 Torr against Kansas State's 3.4e7.

The strict expected-failure test `test_basic_network_against_eden_b_only_decays` records the Eden conflict. With k_BC = 1.2e-10, LIGHT's fitted rate is ×1.35, ×0.95 and ×0.82 of Eden's at 50, 100 and 200 Torr, and it falls with pressure where Eden's rises. The test will start passing, and so must be updated, if the mixing rate is ever revised down.

The table above uses LIGHT's constants of 2026-09-27 before Hakuta was read (k_BC 1.0e-10, τ_C 75 ns). `bc_check.py` now prints the current ones as well. Eden Kr 26 %, KSU 7.5 % and 4.0 %, Hakuta 29 %, B share 0.27.

## What would settle the rest

1. **Yu & Setser, J. Phys. Chem. 94, 2934-2942 (1990), the full article.** It would show how they turned C/B intensities into K_eq: which τ_C they used, and whether the B→A emission under the C-A band was removed. That decides K_eq between 2.5 and 5.
2. **A B-only decay with sub-ns time resolution in 10-200 Torr Kr.** Only this would settle Eden against Hakuta, and nothing in the literature at hand does.

With 1 in hand:
- add its observable to `bc_check.py`
- if K_eq changes, re-derive k_e with `electron_quench.py`
- change `KRF_K_EQ` once, with the tests as the gate.

## Files

| File | What it does |
|---|---|
| `../../v0/two_level.py` | the two-level model in closed form, held to LIGHT's run within 1e-5 by `v0/tests/test_two_level.py` |
| `bc_check.py` | LIGHT against Eden and Kansas State, and the k_BC profile (2 s) |
| `electron_quench.py` | the Trainor & Jacob product re-read (under 1 s) |
| `digitized/` | figure readings, with uncertainties in each header. Hakuta Figs. 2 and 3 were read 2026-09-27 at 600 dpi by fitting equal disks to the data circles. Hakuta's printed F2 slope comes back as 2.38e-10 against 2.5 ± 1e-10 |
