<!-- Survey note A of three, compiled 2026-09-27 by a read-only survey agent
for LIGHT/IMPROVEMENT_PLAN.md. The claims the plan relies on were re-checked
first-hand (see README.md in this directory). Everything else here is a lead
with its file and page, to verify before use. Scratch paths under
~/.claude/jobs/ were temporary and no longer exist. -->

# Slice A: reaction sources outside REBUILD, and what REBUILD lacks

I only read files and wrote nothing under `/Users/cooper/coding`. My scratch files, including the transcriptions and the dedup script, are in `/Users/cooper/.claude/jobs/89373d08/tmp/survey/sliceA/`. Three helper agents inventoried the phantom family, Nathaniel's model with the PIC trees, and the SPECK side (HIBACHI, ARPAE, LXCat). I spot-checked their key claims against the files.

Conventions:
- Paths are relative to `Kinetic_Modeling/` unless shown in full. `Recent/` means `docs/Recent_KrF_Papers/`.
- "p." is the printed journal page. "PDF n" is the page index in the file.
- "Ref pt" is LIGHT's reference point (1 atm, Te 1.31 eV, n_e 1e15, 0.5 % F2). The rates quoted at the ref pt are my own arithmetic from the printed forms.

## Headlines

1. **Your absence claims check out.** In the krf rows:
   - There is no F⁻ detachment of any kind, no three-body electron–ion recombination, no F + F + M, and exactly one attachment row (R10).
   - Every cation has an F⁻ partner, but F2⁺ + F⁻ is a flat 4.0e-8 (FMN1), not a Flannery or Lee form.
   - Dissociative recombination (DR) covers Kr2⁺, Ar2⁺, Ne2⁺, NeKr⁺ and ArKr⁺. There is no e + F2⁺ row, and no repo source prints one.
2. **REBUILD left several of these out on purpose, but judged them at e-beam points.** `SPECK/REBUILD/MODULES/scorecard.py:61-67` and `network.py:976-984` record four omissions: Kr* + e → Kr** ("~1e6 s^-1 at beam-mode ne"), Kannari's excited- and dimer-state ionization ladder, the Kr2*/Ar2*/Ar2F* superelastics ("<=1 %"), and F + F + M ("~1e3 s^-1"). **That judgment does not carry over to LIGHT.** At the ref pt:
   - Kr* + e → Kr** is 3.3e8 s⁻¹, 3.3 times the harpoon rate.
   - Kr** + e → Kr⁺ is 3.4e7 s⁻¹ per Kr**.
   - Kr2* + e → 2Kr + e (Kannari's estimate) is 1e8 s⁻¹, against 8e7 s⁻¹ for Kr2* radiation.
3. **e + Ne excitation and ionization (Chen R1–R3) only enter krf when Ne LoKI buckets are loaded** (`network.py:231, 270-275`). They are absent from the dumped fill. In krf, Ne*/Ne⁺ otherwise come only from the beam (`MODULES/beam.py:68-69`). LIGHT v1 drives the electron gas, so it needs its own versions.
4. **For v0's F2 burn-out, the row that matters is KrF* + F, not F + F + M.**
   - F + F + M runs at about 6e3 s⁻¹ per F atom, i.e. on a 0.2 ms timescale.
   - KrF* + F at Kannari's 7.8(−10) gives 1.9e8 s⁻¹ once the 0.5 % fill has burnt, which is above KrF(B) radiation at 1.43e8.
   - The printed values for KrF* + F span 1e-12 to 7.8e-10.
5. **No printed value exists anywhere in the repo for:**
   - F⁻ detachment by electrons, by F (associative), by collision, or by Kr*/Ar*;
   - e + F2⁺;
   - any Kr⁺ + e recombination;
   - 248-nm absorption by Kr2*, the Ne-family species, or KrF*.
   
   For F⁻ detachment the only numbers are "Estimated" rows in `docs/KrF_Reaction_Rate_Tables.pdf` p.7, and phantom copies of them. For the rest there are only analogs (Ar⁺, or other wavelengths).

## 1. Inventory

Literature (read first-hand; values checked on the page images where the text layer was garbled):

| Source (file) | Table, pages | Rows | Gases | Provenance |
|---|---|---|---|---|
| Johnson & Hunter 1980 (`Recent/Johnson.pdf`) | T.I p.2408 (PDF 4) | 82 | Kr/Ar/F2 | Primary compilation. 11 rows "Estimated". Electron rows are single "typical" values at 2 kV/cm/atm. |
| Kannari 1985 (`Recent/kannari1985.pdf`) | T.I pp.4311-4313 (PDF 4-6) | 119, incl. 5 e-beam W rows | Kr/Ar/F2 | Primary, refs per row. 17 rows "estimated" or "analogous". Footnote a: the quench rates are lifetime-rescaled (×0.433 for KrF*). |
| Lee 1989 (`Recent/lee1989.pdf`) | T.I pp.2055-2056 (PDF 3-4) | 97 | Ne/Ar/F2 | Primary ArF model. 6 rows "Guessed". Some exponent misprints (e.g. "Ne2* → 2Ne + hν 8.9 × 10^-7 s^-1"). |
| Petrova 2022, NRL Orestes (`Recent/petrova2022_orestes.pdf`) | App. T.A1-A5 pp.20-26 (PDF 27-33) | 65 | Ar/F2 | Refs per row, plus 2022 updates in red. Some rows "this work" or "estimate". |
| Boichenko 2000 (`Recent/Boichenko_2000.pdf`) | T.1 pp.683-684 (PDF 4-5) | 83 | Ne/Ar/Kr, no F | Ion chemistry. Most three-body conversions are set to an assumed 10^-31. Absorption is at 210-224 nm. |
| Chow 1978 (`docs/KrF_Chow.pdf`) | T.1 p.23 | 29, + 3 cavity terms | He/Kr/F2 | Report. Linear-threshold Maxwellian fits, no per-row refs. |
| Ong 1996 (`Recent/ong1996.pdf`) | T.1 pp.2588-2589 (PDF 4-5) | 93 | He/Kr/F2 | Refs per row. 21 rows come from a Boltzmann solver. |
| Ma 2024 (`Recent/Ma_2024.pdf`) | T.5 pp.17-20 (PDF 18-21), SI units | 152 | Ne/Ar/Xe/F2 | Refs per row, 2 "guess". |
| Chen 2025 (`Recent/chen2025.pdf`) | T.III (PDF 11-12) | 41 | Ne/Ar/Xe/F2 | Already carried in REBUILD `MODULES/chen2025.py:89-140`, with errata at :156-172. |
| Xu 1993 (`Recent/xu1993.pdf`) | T.I p.2594, T.II p.2596 | 40 + 4 | KrF(B,C) and Kr2F* against 20 quenchers | Measurements |
| Gadomski 1992 (`Recent/gadomski1992.pdf`) | pp.153-158 | 4 | Kr/Xe/F2 | Measured; superseded by Xu 1993 |
| Kvaran 1988 (`Recent/Kvaran1988.pdf`) | T.3 p.101 | KrF(B,v') relaxation by He, Ne and Ar, per v' | | Measured |
| Johnson 1989, XeCl (`Recent/johnson1989.pdf`) | T.I pp.5714-5718 (PDF 9-13) | ~200 (rough count; noisy text layer), 43 of them Ar/Ne-only | Xe/Ar/Ne/HCl | Refs and notes per row |
| Rauf 1999 (`Recent/rauf1999.pdf`) | T.I p.3464 (PDF 6) | ~65 (rough), 10 Ne-only | He/Ne/Xe | Refs per row |
| Rhodes/Brau, *Excimer Lasers* (`docs/Excimer_Lasers_Topics_Rhodes.pdf`) | T.4.3-4.12 (PDF 105-118), T.6.4-6.5 (PDF 199-201) | not counted | all rare gases | Review of measurements |
| Kushner & Moratz 1988 (`Recent/kushner1988.pdf`) | p.1857 | 2 fitted values: k_d ≈ 4.6e-9, β ≈ 1.4e-9 cm³/s (XeF, Te ≈ 2 eV) | F2 | Fitted |
| `docs/KrF_Reaction_Rate_Tables.pdf` (Watson, Jan 2026) | 17 pp. | ~120 | He/Kr/F2 + Ar | Secondary compilation: Chow T.1 + a "Janacaitis" list (rows 53-71, p.4) + J&H renumbered +100 + "model stability" rows "Estimated from [34]" (p.7). Several "[?]" citations. |

No rate tables in: `nighan1978`, `petrov2017` (Ar/F2 e-beam yields only), `ArF_report`, `AD1329084`, `obenschain2020`, `wolford2020`, `turnbull*`, `Langdon1980`, `Matte1988` and `dudey1967`. `kushner_krfs.pdf` is image-only and I could not check it.

Code and data:

| Source | Rows | Gases | Provenance |
|---|---|---|---|
| **Phantom monolith** `KrF_and_LoKI/KrF Code/og-pyfiles/KrF_phantom11.py` (constants :324-471, RHS :1893-2181) | 118 | He/Kr/F2 + Ar | Chow, Jancaitis, J&H and Kannari. The literature completions sit behind LIT_* flags that are off by default for parity. **Uncited:** k26B, K163, k163_inv0, k10_inv, k_F_detach, r59_Ar, K151, K158. **Calibrated:** K26_3BODY, K162_ALPHA_F2DD. |
| **Phantom factored** `factored-pyfiles/phantom11/config.py:222-776`, `reactions.py:33-304`, `_rhs_u_body.py:171-683` | 128 | same | Monolith + 10 gated rows |
| Phantom copies and subsets | legacy_11 28, ode_13 46, spatial 118, pic1d `fluid.py` 118, light_chem 40, unified kinetics 30 (flux replay), krf0d 0 | | |
| spectre `factored-pyfiles/spectre/rates.py:69-279` | 78 | Kr/F2/Ne/Ar | REBUILD's ancestor (same row ids) |
| `presentation/architecture/chemistry.txt:49-213` | 105 + 8 | | Documentation of phantom11 |
| **Nathaniel** `/Users/cooper/coding/nathaniel/1DModelKrF/data/reactions/krf.toml:39-486` | 45 | Ne/Kr/F2 | Independent; a DOI on every row |
| Nathaniel `arf.toml` | 44 | Ne/Ar/F2/Xe | same |
| Nathaniel audits `rate-references/kannari1985_KrF_rates_final.csv` | 119 | | Status per row. The Kannari electron rows are marked "OCR-transcribed (exponents uncertain)". |
| Nathaniel audits `rate-references/ne_kr_f_rates.csv` | 69 | | Status per row (MEASURED/MODEL/ESTIMATE/THEORY) |
| `/Users/cooper/coding/PIC/Chen2025/scripts/chen_table3.py:36-78` | 41 + 7 errata + 10 Ma rows (:173-179) | | Same transcription as REBUILD's |
| `/Users/cooper/coding/PIC/Neon/patches/rebuild-0.1.38-kr-partner-rows.patch:9-65` | 5 | | All 5 are now in 0.1.40 |
| `LIGHT/v0/network.py:149-331` | 11 + 27 optional | Kr/F2 | REBUILD copies, plus its own SUPER row |
| LXCat pool `LoKI-MC/Code/Input/*` | see below | | Only the Ne BSR set is used by the krf0d decks and REBUILD's V2point_NeKrF2 table |
| PIC_MCC_3D and ARPAE (`SPECK/PIC_MCC_3D/MODULES/cross_sections.py:93-98`) | e-neutral only | Ar/Kr/F2/He, no Ne | No ion-neutral processes |
| HIBACHI, SPECKLE3D, SPECKLE, Argon_UI, KrF_and_LoKI/docs | none of their own | | Reuse REBUILD, or read Nathaniel's krf.toml. Argon_UI has only Ar 4p→4s A-values (`emission.py:38-45`). |

LXCat pool contents:
- **Kr:** Biagi-v7.1, 4 lumped groups (5s, 5p, 4d+6s, 5d; `Krypton/Kr_LXCat.txt:275, 489, 703, 917`) + ionization. There is no level-resolved Kr set anywhere in the repo.
- **Ar:** 37 excitations.
- **Ne:** BSR 30 excitations, or Biagi 45.
- **F2:** SIGLO, 11 blocks.
- **F atom:** BSR, 353 blocks (unused).

`old_files/` exists and was not used.

## 2. Candidates not in REBUILD's krf rows

Values are as printed, in cm³/s, cm⁶/s, s⁻¹ or cm² unless stated. Te is in eV.

| # | Reaction | Printed value | Where | When it matters |
|---|---|---|---|---|
| A1 | e + e + Kr⁺ → Kr** + e | **No Kr row.** Ar analogs: 7.2e-27 Te^-4.5 cm⁶ s⁻¹ (Kossyi-1992); 8.75e-39 (Te)^-4.5 m⁶/s → Ar*; 5.0e-39 → Ar** | Petrova p.23 (PDF 30); Ma #41, #42 p.18 (PDF 19) | Afterglow with free n_e. Per ion: 2.6e3 s⁻¹ at Te 1.31 eV, but 2.8e8 s⁻¹ at Te 0.1 eV and n_e 1e15, beating Kr⁺ → Kr2⁺ (1.5e8 at 1 atm) |
| A2 | e + e + Rg2⁺ | Only He2⁺: 2.8(−20)(Tg/Te)^4 and 8.0(−21)(Tg/Te)^4 cm⁶/s | Ong #38, #39 p.2588 | Adds to DR in a cold, dense afterglow |
| A3 | e + Kr⁺ + Kr → Kr* + Kr | **No Kr row.** "Ar⁺ + e + Ar → Ar* + Ar 1.0 × 10^-38" m⁶/s, with no Te dependence | Ma #43 p.18 | Taken by analogy it is 2.4e8 s⁻¹ per ion (1 atm, n_e 1e15), competing with Kr2⁺ formation. Check before use. |
| A4 | Kr2⁺ + e → Kr** + Kr (a DR branch) | 1.9(−07) [68-70], beside the Kr* branch 1.2(−6)[Te(K)/300]^−0.55 [71] | Kannari p.4312 (PDF 5) | REBUILD's R11p sends all DR to Kr*. Routing to Kr** opens photoionization (6.0(−18) cm²) and C2. |
| | | J&H #60: 1.1e-7, all to Kr** | J&H p.2408 | |
| | | Chow #11: K0 = 2.4e-7 Te^-1/2 | Chow p.23 | |
| | | 1.6e-7 Te^-1/2 | Boichenko #15 p.683 | |
| | | 1.6e-6 (300/Te)^0.55 total, all to Kr** (Shiu & Biondi 1977) | `nathaniel/.../ne_kr_f_rates.csv:38` | |
| A5 | e + F2⁺ → F + F | none found | — | could not find |
| A6 | Ne⁺/Ne2⁺ + F⁻ + M → neutral products (no NeF*) | 4.19e-11 × N and 4.82e-11 × N. N is not defined; I could not check the units. | Ma #66, #67 p.18 | Ne buffer. REBUILD sends every pair to NeF* (NF1, NF2). |
| B1 | F + F + M → F2 + M | J&H #44: 1.0e-33 cm⁶/sec "Estimated" | J&H p.2408 | About 6e3 s⁻¹ per F at 1 atm. Matters for interpulse and rep-rate F2 recovery only. |
| | | 1.0(−33) [16] | Kannari p.4311 | |
| | | F + F + He: 1.0(−33) [22] | Ong #81 p.2589 | |
| | | Phantom K144 = K3(1.0e-32) "Baulch 1981"; K144_LIT = K3(1.0e-33); low edge 8.0e-35 | `factored-pyfiles/phantom11/config.py:384-401` | |
| B2 | F + e → F⁻ | J&H #63: 1.0e-12, flagged "effective two-body at 2 atm", Estimated | J&H p.2408 | About 2.4e5 s⁻¹ on a fully burnt fill, vs 5.5e8 attachment on intact F2. Matters only after >99.9 % burn-out. |
| | | 1.0e-12 [31] | Lee p.2055 | |
| | | 1.0e-12 ("Lee-1989") | Petrova p.23 | |
| | | 1.0(−12) [6] | Ong #26 | |
| | | Phantom K163 | `config.py:384` | |
| B3 | F⁻ + e → F + 2e | "3.0 × 10^-7 cm3 s-1 × e^-0.2 eV/kT", "Estimated from [34]". [34] is Flannery & Yang, a paper on ion-ion recombination. | rate tables p.7, p.13; phantom `config.py:736` | No primary source. 0.2 eV contradicts the code's own EA(F) = 3.45 eV (`config.py:223`). |
| B4 | F⁻ + F → F2 + e | 1.0e-9, "Estimated from [34]" | p.7; `config.py:737` | No primary source; roughly 1.8 eV endothermic |
| B5 | F⁻ + M → F + M + e | Rate table: "F⁻ + He → F + He + e 5.0 × 10^-12"; phantom `k_F_detach = 1e-6` for Ar | p.7; `config.py:738` | Three inconsistent estimates. Under phantom's K2 convention 1e-6 is 1e-13 cm³/s; `chemistry.txt:119` prints it as "1e-6 cm3/s". |
| B6 | F⁻ + Kr*/Ar* (detachment) | none | — | could not find |
| B7 | **KrF* + F → Kr + 2F** | Kannari: 7.8(−10), "analogous to F2^a" | Kannari p.4311 | Burn-out phase: 1.9e8 s⁻¹ at 7.8e-10, above KrF(B) radiation. The 780× spread among the printed values is itself the finding. |
| | | Ong #67: 7.8(−10) [22]; #73 for KrF** | Ong p.2589 | |
| | | Phantom K_KRF_F_LIT = K2(3.0e-10), re-based on the measured KrF* + F2 | `config.py:731` | |
| | | Analog: "ArF* + F → Ar + 2F 1.0 × 10^-12, Rokni-1977, k=0, not known" | Petrova p.25 | |
| B8 | F2* chemistry | Lee: F* + F2 → F2* + F 5.1e-10 [32]; F2* + F2 → products 3.5e-10; F2* → F2 + hν 2.7e8 | Lee p.2056 | Only with beam or hot-electron F2 excitation (F2* at 11.57 eV) |
| | | Janacaitis list row 64: "Kr + F2* → KrF* + F 7.6 × 10^-10 [7,9]" | rate tables p.4 | |
| B9 | NeF* → Ne + F* | "2/3 of formed NeF*"; F* + Ar → Ar* + F 1.0e-10 "Guessed" | Lee pp.2055-2056 | Ne buffer. REBUILD's NF3 dumps NeF* to heat, losing the 12.7-eV F* (energy from Petrova p.3, PDF 10). |
| C1 | **Kr* + e ↔ Kr** + e** | Kannari: 9.00(−07)Te^0.50 exp(−1.51/Te) [76], two-way | Kannari p.4311 | 3.3e8 s⁻¹ at the ref pt, 3.3× the harpoon |
| | | J&H #47: 6.4e-7 and #50: 8.0e-7 | J&H p.2408 | |
| | | Chow #3: E0 = 0.65 eV, α = 1.5 × 10^-15 cm²/eV; #8 is its reverse | Chow p.23 | |
| | | Ong #11: 6.4(−7); #27: 8.0(−7) | Ong p.2588 | |
| C2 | **Kr** + e → Kr⁺ + 2e** | Kannari: 2.00(−07)Te^0.62 exp(−2.53/Te) [79] | Kannari p.4311 | 3.4e7 s⁻¹ per Kr** at the ref pt, 9× Kannari's Kr* stepwise rate. Scales up by 10× at 1e16. |
| | | J&H #52: 1.8e-7; Ong #14: 1.8(−7); Chow #7: E0 = 3.4 eV | same pages | |
| C3 | Kr** → Kr* + hν | Chow #18: τ_Kr** = 2.76 × 10^-8 sec; phantom r18 | Chow p.23 | 13 % of Kr** loss at 0.1 atm; 1.5 % at 1 atm |
| C4 | Kr** + Kr → Kr* + Kr; Kr** + 2Kr → Kr* + 2Kr; Kr** + 2Kr → Kr2* + Kr | 3.0e-11 (ref 22); 5.4e-32 [4]; 2.6e-32 [4] | Janacaitis list, rate tables p.4 rows 55, 57, 58 | Sourced (second-hand) alternative to the EST KPS_Kr (1e-10): 7.6e8 vs 2.4e9 s⁻¹ at 1 atm |
| C5 | Kr* + e → Kr + e | J&H #51: 8.0e-11; Ong #28: 1.0(−8) [20] | J&H p.2408; Ong p.2588 | LIGHT's SUPER already covers it. LIGHT claims REBUILD carries it inside its swarm tables (`LIGHT/v0/network.py:270-271`); I could not check that. |
| C6 | Kr2* + e → Kr2⁺ + 2e | Kannari: 9.00(−08)Te^0.70 exp(−3.53/Te) [72]; Ong #17: 1.0(−10) | Kannari p.4311; Ong p.2588 | 7.4e6 s⁻¹ at the ref pt, 9 % of the Kr2* radiative rate; comparable to it at 1e16 |
| C7 | Kr2* + e → 2Kr + e | Kannari: 1.0(−07) "estimated" | Kannari p.4312 | Kannari's value gives 1e8 s⁻¹ at the ref pt, above radiation. The sources spread 100×. |
| | | Ong #30: 1.0(−9) [22] | Ong p.2588 | |
| | | Janacaitis list row 62: 5.7 × 10^-9 (kTe/1 eV)^1/2 | rate tables p.4 | |
| C8 | e + Kr → Kr***(6s, 4d), Kr**** | Boltzmann-solver rates | Ong #6, #7 | Beyond 5s/5p. The only in-repo data is Biagi-v7.1's lumped groups. |
| D1 | KrF(Bv) → KrF(B) vibrational relaxation. REBUILD's VREL_Ne and VREL_Kr are EST 4e-11. | Ong #60: KrF** + Kr → KrF* + Kr 1.0(−10) [13] | Ong p.2589 | Sets the B(v) feed rate. Only matters if Bv is tracked. |
| | | Kvaran: Ne (1.2±0.5)–(3.5±0.8), Ar (42±12)–(58±12), in units of 10^-11 | Kvaran T.3 p.101 | |
| | | Kvaran text: Ar "2.5-5 × 10^-11", "similar … Kr-buffered" | Kvaran p.95 | |
| D2 | Kr2F* + 2Kr → 4Kr + F | Ong #75: 2.9(−31) cm⁶ s⁻¹ [40] | Ong p.2589 | **Do not adopt.** It would be 1.7e8 s⁻¹ (40× Kr2F* radiation) and contradicts Xu p.2599 ("Kr does not quench Kr2F* even at 600 Torr") and Xu T.I (< 0.6 × 10^-13). |
| D3 | Kr2⁺ + F⁻ (+M) → Kr2F* + M | ~1e-6, THEORY (Flannery & Yang) | `ne_kr_f_rates.csv:32` | A trimer branch that REBUILD's K108 lacks |
| D4 | KrF(C) + F2; B↔C mixing by Kr | 2.5e-10; 1.2e-7 "MEASURED" (Hakuta 1990) | `ne_kr_f_rates.csv:43-44` | 1.2e-7 is about 10³× gas-kinetic as written, so the units are suspect. Could not check (paper not in repo). |
| D5 | Impurity quenching of KrF(B,C) and Kr2F* | HF (3.1±0.4)e-10 / (3.7±0.4)e-11; Xe (1.0±0.2)e-10; N2, O2, CF4 and others | Xu T.I p.2594 | Only if impurities are modelled |
| E1 | Kr2* + hν → Kr2⁺ + e | 10^-18, but only at 210-224 nm | Boichenko #79 p.684 | No 248-nm value in the repo |
| E2 | Ne-family absorbers | Lee (193 nm): Ne2⁺ 2.3e-19, Ne2* 9.2e-19, Ne** 2.5e-19, Ne* 1.3e-19 | Lee p.2056 | No 248-nm values |
| E3 | ArKr⁺ + hν → Ar⁺ + Kr | 1.5e-17: Nathaniel's reading of Kannari's [85] as the ion. Kannari prints "ArKr* + hν → Ar⁺ + Kr + hν". | `kannari1985_KrF_rates_final.csv:120`; Kannari p.4313 | Conflicts with REBUILD's ArKr* reading (`network.py:1449-1458`) |
| E4 | Preionizer VUV | Kr: σ 3.5e-17 (Samson 1966); F2: 5.0e-20 / 3.0e-19; Ar → Ar*: 1.0e-19 | phantom `config.py:45-68` | Discharge start-up only |
| F1 | e + Ne → Ne*, Ne**, Ne⁺ | 2.55 × 10^-9 Te^1.69 exp(−24.9/Te); 4.82 × 10^-10 Te^0.48 exp(−27.8/Te); 1.23 × 10^-9 Te^0.72 exp(−32.4/Te) | Lee p.2055; SIGLO in Nathaniel `krf.toml:39-64` | Needed for LIGHT v1 |
| F2 | Ne* + e ↔ Ne**; Ne2* + e → 2Ne + e | 3.83 × 10^-7 Te^0.32 exp(−2.85/Te) | Lee p.2055 | Ne buffer at n_e ≳ 1e15 |
| | | Ne2* + e: 1.0e-7 "Guessed" | Lee | |
| | | Ne2* + e: 3.0e-13 m³/s | Ma #57 | |
| F3 | Ne/Kr ion conversion | #22, #24, #34, #38, #47, #49, all "10^-31" | Boichenko p.683-684 | Ne buffer. The value is an assumed class. |
| | | Ne⁺ + 2Kr ~1e-31; Kr⁺ + 2Ne ≤1e-32 | `ne_kr_f_rates.csv:28-29` | |
| F4 | Ar buffer | Ar* + e → Ar** 8.90(−07)Te^0.51 exp(−1.590/Te); Ar2* + e → Ar2⁺ + 2e 9.00(−08)Te^0.70 exp(−3.66/Te); Ar2* + e → 2Ar + e 1.0(−07); Ar2F* + e 1.0(−07) [17]; Ar2⁺ + e → Ar** + Ar 1.1(−07) | Kannari pp.4311-4312 | Mirrors of C1–C7 |
| | | ArF* + 2Ar → 3Ar + F 5.0e-32 | Lee p.2056, Petrova p.25 | |
| | | Ar* + e → Ar + e 5.4e-12 | J&H #58 | |

## 3. Value problems on rows REBUILD already has

- **NK6** (Kr* + Kr + Ne → Kr2* + Ne) is 2.5e-43 m⁶/s, cited as an "analog of Chen R23" (`network.py:499-500`).
  - That is Chen's misprinted R23. REBUILD's own errata correct R23 to 1.1e-44 (`chen2025.py:64, 162`).
  - NK6 therefore runs 9.3× above AK8/AK8_Kr. Nathaniel uses 2.7e-32 cm⁶/s (`krf.toml:245-251`).
- **Printed values exist for several EST rows:**
  - AK8_Kr: Rhodes T.4.7 p.102 (PDF 112) gives 2.59 × 10^-32 (Tracy & Oskam) and 4.4 × 10^-32. Chow #16 gives 4.5e-32; Ong #55 gives 2.6(−32).
  - KSTEP: Kannari prints 7.00(−08)Te^0.68 exp(−4.04/Te).
  - NK1: 9.0e-11, measured (Neynaber & Magnuson; `krf.toml:196-204`).
  - NK3: ≤1e-14, a measured bound (Johnsen 1978; `krf.toml:269-275`), against REBUILD's EST 1e-11.
  - KPS_Kr: see C4.
  - VREL_Ne and VREL_Kr: see D1.
- **Phantom r120** is Kr2⁺ + F2 → Kr2F* + F (`config.py:510`, `_rhs_u_body.py:331`). J&H #20 prints Kr2* + F2. Do not copy rows from phantom without checking them.
- **The Stevefelt row** (Ar⁺ + e + Ar; rate tables p.14; `config.py:657-669`) does not agree with its own check values. The formula 1.09e-26 Te[K]^-4.5 gives 4.0e-48 at 5 eV, where the printed check says 6.7e-31. Its Te^-4.5 law also belongs to an electron third body, not Ar.

## 4. Rough count

**Dedup rule.** Two rows are the same reaction when they have the same reactants and products after:
- lumping KrF(B)/(C)/(Bv)/KrF** into KrF*, and ArF(B)/(C) into ArF*;
- dropping the photon from spontaneous emission;
- collapsing any quench of an excited species to ground-state products into one reaction per collision partner, however the products are written (so Kr + F + F2 and Kr + 3F count once);
- counting ion–ion recombination once per (cation, excimer product), whatever third body is printed. Other three-body rows with a different third body stay distinct.

Rate values and sources are ignored. Cavity loss terms are excluded. The e-beam W rows count as reactions.

**Results** (script `sliceA/dedup.py` and `count_all.py`):

- **Literature tables transcribed:** 15 sources, 893 rows. For Johnson 1989, Rauf and Rhodes I transcribed only their Ar/Ne or relevant rows.
  - That gives **429 distinct reactions**: 283 within Kr/Ar/Ne/F, and 70 within Kr/F only.
- **Adding REBUILD** (143 distinct after lumping) and about 20 code-only rows (phantom, Nathaniel):
  - **~459 distinct**, of which **~310 are within Kr/Ar/Ne/F** and **~79 are Kr/F-only**.
- **Of the 283 literature reactions within Kr/Ar/Ne/F, 130 are not in REBUILD:**

  | Group | Count |
  |---|---|
  | Kr/F-only | 26 |
  | involving Ar but not Ne | 36 |
  | involving Ne but not Ar | 29 |
  | involving both Ar and Ne | 39 |

  46 of the 130 come only from Boichenko: mostly assumed-10^-31 ion conversions and NeAr⁺/NeKr⁺ ionic-excimer rows.
- **All gas families** (adding the ~150 untranscribed Xe/Cl/HCl rows of Johnson 1989, Rauf's He/Xe rows and the Rhodes entries): roughly 650.
- **The LXCat pool adds level-resolved electron-impact processes on top:** about 130 from the ground state in the Kr/Ar/Ne/F family (Kr 6, Ar 39, Ne 47, F2 11, F 28), or about 450 counting the F-atom file's excited initial states. REBUILD collapses each gas to 2 buckets.
