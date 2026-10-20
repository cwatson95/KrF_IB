<!-- Survey note B of three, compiled 2026-09-27 by a read-only survey agent
for LIGHT/IMPROVEMENT_PLAN.md. The claims the plan relies on were re-checked
first-hand (see README.md in this directory). Everything else here is a lead
with its file and page, to verify before use. Scratch paths under
~/.claude/jobs/ were temporary and no longer exist. -->

# Slice B: measured data a KrF 0-D model can be scored against

Path keys: KM = /Users/cooper/coding/Agent/projects/Kinetic_Modeling, PAP = KM/docs/Recent_KrF_Papers, RB = KM/SPECK/REBUILD, LEG = "KM/KrF_and_LoKI/KrF Code/factored-pyfiles", VAL = /Users/cooper/coding/PIC/Validation, NEON = /Users/cooper/coding/PIC/Neon. "p." = printed journal page. Nothing under /Users/cooper/coding was run or written.

### Headline
1. **The only in-tree data LIGHT v0 can reproduce with no new physics** are the photolytic decay measurements of Gadomski 1992 and Xu 1993: pure Kr/F2, 300 K, no electrons, no cavity. LIGHT v0 misses them. Its KrF(B) decays ×2.4 to ×2.9 too fast at low pressure. With the quench and bc groups on it is still ×1.5 to ×1.6 too fast.
2. **Every KrF gain or device dataset in the tree is Ar-diluted (4–10 % Kr).** Using them needs v1 (Ar), an e-beam W/yield source term instead of Maxwellian heating, the stimulated-emission cross-section and the 248-nm absorbers, and a cavity for device energies.
3. **Pure-Kr device data exist only as abstracts or derived constants:** Mandl 1985 (99.6 % Kr), Suda 1984 (buffer-free) and Trainor & Jacob 1980 (Kr/F2 99.5/0.5). Their conditions and raw numbers are not on disk.
4. **REBUILD has only two firm gain anchors:** Klimek (KrF) and Electra (ArF). It never scores the direct decay data and scores no Kr-rich data. KrF gain points with printed deposition (Aurora, Tacoma) exist only in legacy docstrings.

### 1. Data that isolate kinetics (no electrons, no cavity)

| Dataset | Measured | Conditions | Where | REBUILD use | LIGHT needs |
|---|---|---|---|---|---|
| Gadomski, Xu, Setser & Nelson 1992 | KrF(B,C) decay constant vs F2 (Fig. 1) and vs Kr (Fig. 2). Fit: kF2 (2.8±0.2)e-10, kKr (4.0±1.0)e-12, k'Kr (6.0±2.0)e-31, K_eq 5.0±1.0, τ_eff ≈25 ns (τB 7.0 and τC 75 ns assumed). kBC was *assumed* 1.0e-10 | 7-ns photolysis at 240–242 nm, 300 K; Kr 50–500 Torr with F2 1 or 2 Torr; F2 0–12 Torr at Kr 50 and 100 Torr | PAP/gadomski1992.pdf p.153, p.155 (figures, eq. 5), p.156 | Rate source only (RB/MODULES/network.py:359-375). Measured K_eq recorded but not adopted (RB/CHANGELOG.md:345-352) | Seed mode (`init_extra_cm3`, LIGHT/v0/model.py:71) with power 0 and n_e ≈ 0; quench + bc groups; late-time exponential fit. Pure Kr/F2 |
| Xu, Gadomski & Setser 1993 | k' = 5.0e7 and 6.2e7 s⁻¹ at 50/100 Torr Kr + 1 Torr F2; intercepts 6.0–6.5e7 measured. Quench table: F2 (3.0±0.3)e-10. He/Ne/Ar series: two-body 2.5/2.1/2.9e-12, mixed three-body 4.6/4.9/5.5e-31. 2Ar: "no evidence … ≤400 Torr". Kr2F*: 220±20 ns, kF2 (1.2±0.2)e-10, Kr < 0.6e-13. Spectrum of B-X, C-A and Kr2F | 243 nm, 7 ns, 300 K; Kr 50–600 Torr; F2 1–2 Torr | PAP/xu1993.pdf p.2594 (Figs 3–4, Table I), p.2595–2596 (Figs 5–8, Table II), p.2597 (Figs 9–10), p.2592 (Fig. 1) | Rate source (network.py:352-416). K125 kept at 1.5e-10, not Xu's 1.2e-10 | Pure-Kr parts as above. The He/Ne/Ar series needs v1 (He is not in the v1 plan) |
| Kvaran, Shaw & Simons 1988 | Effective relaxation constant vs vibrational level, He/Ne/Ar: Ar 42–58e-11, He 0.8–4.0e-11, Ne 1.2–3.5e-11 | Kr resonance lamp + NF3, 0.25 + 0.25 Torr plus buffer. **Kr buffer deliberately avoided** (p.100) | PAP/Kvaran1988.pdf p.101, Table 3 | Named "primary" for VREL_Kr (network.py:1217-1220), but it has no Kr data | Needs a vibrational ladder. Not scorable by LIGHT |
| Rhodes (ed.) 1984, Tables 4.4 and 6.5 | Kr*+F2 → KrF*: 7.2e-10, branching 1.0. Kr2F* lifetime 145–185 ns. KrF*+2Kr 2.9–9.7e-31 | Compilations | KM/docs/Excimer_Lasers_Topics_Rhodes.pdf PDF p108, p201 | Unmined | A direct check of G3 (LIGHT uses 8.1e-10) |

**LIGHT against the pure-Kr decay data**, from LIGHT's constants at 300 K with n_e = 0 (the survey computed these from the printed constants, and LIGHT itself reproduces them; see README.md):

| Kr/F2 (Torr) | measured | v0 | v0+quench | v0+quench+bc |
|---|---|---|---|---|
| 50/1 | 5.0e7 (Xu p.2594) | 1.45e8 (×2.9) | 1.61e8 (×3.2) | 8.1e7 (×1.6) |
| 100/1 | 6.2e7 | 1.49e8 (×2.4) | 1.72e8 (×2.8) | 9.3e7 (×1.5) |
| 500/1 | 2.63e8 (Gadomski fit curve) | ×1.13 | ×1.42 | ×1.12 |

- The miss is the zero-pressure intercept. LIGHT's 1/τ is 1.43e8 s⁻¹ in v0 (B only) or 6.6e7 with bc (K_eq 1.47, τ_eff 15.2 ns). The measured coupled τ_eff is 25±2 ns (Xu p.2591, p.2594).
- LIGHT's K125 gives ×1.22 the measured Kr2F* decay constant at 300 Torr Kr and 10 Torr F2.
- LIGHT copies two-body and three-body constants that Gadomski and Xu extracted assuming K_eq ≈ 5, then runs them under a different B/C structure. At 1 atm Kr, 0.5 % F2 and 300 K the B-X share of KrF* decay is 0.285 (v0), 0.104 (bc) and 0.045 (Gadomski's structure).
- Provenance class: "consistency" (KM/SPECK/EXPERIMENTS.md:16-20). It tests structure, not blind prediction.
- The F2 came as 10 % F2 in He (gadomski1992 p.154). That adds at most 1.4e6 s⁻¹ via Xu's He quench, which is negligible.

### 2. Pure-Kr and Kr-rich data (the buffer-free set)
- **Trainor & Jacob 1980:** e-beam, Kr/F2 99.5/0.5, 190 Torr, ~70 kW/cm³, KrF(B) emission vs F2. These are LIGHT-native conditions and need no cavity. The tree has only derived electron-quench values: 2.0e-7 (theirs) and 3–6e-8 (Kushner's reassessment) (PAP/kushner_krfs.pdf p.2482).
- **Suda, Obara & Fujioka 1984:** Kr/F2 at 1 atm, 4.2 J/l at 1.7 MW/cm³, ~5 % intrinsic. Abstract-level only (NEON/LITERATURE.md:6, :68).
- **Mandl, Klimek & Salesky 1985:** output efficiency 9.3/11.2/9.2 % at 4/10/99.6 % Kr in Ar, 350-ns pulse (NEON/LITERATURE.md:69). Scored only in the PIC bundle, at assumed Electra conditions (1.36 atm, 0.3 % F2, 0.7 MW/cm³, output coupler 0.4; NEON/scripts/mandl_check.py:4-9, 23-38). The 99.6 %/10 % ratio came out 0.27 against 0.82 measured, and 4 %/10 % 0.74 against 0.83 (NEON/NRL_REPLICATION.md:378-381; RB/CHANGELOG.md:627-633). That was REBUILD 0.1.38, not re-run after 0.1.40. REBUILD itself only gates that a pure fill lases (RB/tests/test_engine.py:170-193).
- **Referenced, not in the tree:** Kannari, Shaw & O'Neill 1987 (efficiency flat up to 1.25 MW/cm³ at 10–99.7 % Kr; kushner_krfs.pdf p.2484). Suda, Kumagai & Obara 1987 (APL 51, 218) and Lee, Kannari & Obara, both Kr-rich amplifiers (PAP/lee1989.pdf reference list). Salesky & Kimura 1985: net gain 3.8 %/cm at 0.12 MW/cm³, 4 % Kr (KM/PLAN_physical_validity.txt:3118-3119).

### 3. Ar-buffered KrF gain and device data (need LIGHT v1, an e-beam source, cross-sections and absorbers)

| Dataset | Observable | Where | REBUILD score |
|---|---|---|---|
| Klimek 1981 grid (via Kannari Tables II/III) | 18 gain + 18 absorption values (F2 excluded) at 110/55/20 A/cm² × 1–3 amagat × 0.2/0.6 % F2. Kr %, beam energy and pump not printed | PAP/kannari1985.pdf p.4317; F2 convention p.4317–4318 | 0.1.40: gain ×0.64 (0.48–0.79), absorption ×0.68, g0/α ×0.96 (1 amagat) and ×0.86 (3 amagat). Pressure exponent +0.57 vs +0.45 at 110 A/cm², +0.01 vs +0.13 at 20 A/cm² (RB/CHANGELOG.md:288-293; RB/tools/excimer_anchors.py) |
| Aurora preamplifier; Tacoma | g0 1.96 %/cm at 0.080 MW/cm³ (10 % Kr, 0.5 % F2, 750 Torr); g0 5.45 %/cm at 0.285 MW/cm³ with α 0.76–0.79 (1.32 atm, 14.4 % Kr) | LEG/light_unified/benchmark.py:66-92; LEG/unified/validation_ebeam.py:220-280 | Not in REBUILD |
| Johnson & Hunter 1980, Tables II/III | AVCO 102 J, 9.8 %, α 0.27 %/cm, g0/α ~10. Maxwell 120 J, 7.2 %, g0/α ~6. Northrop 54 J, 15 %. E-beam-sustained discharge 36 J, 3 %, ~4 | PAP/Johnson.pdf p.2415 | Not scored. The authors' code was tuned on one case (p.2416) |
| Kannari 1985 device | Fig. 5 output vs output-coupler reflectivity (0.6/1.3/2.0/2.3 J at 80–20 %); Fig. 2 fluorescence; Fig. 6 depletion vs flux; Fig. 7 Tisone data; Table IV six measured intrinsic efficiencies | kannari1985 p.4315–4319 | Fig. 5 only, last recorded ×0.86/×0.23 (KM/SPECK/COMMERCIAL_AUDIT.md:433); reflectivity trend inverted (VAL/validation_report.tex:177-179) |
| Rhodes / Brau §4.3 | 350 J, ~10 % (1 atm Ar, 7 % Kr); 15 % and 13 % intrinsic | Rhodes PDF p123 | Not scored |
| Electra KrF | 730 J best, 400–700 J routine | KM/HIBACHI/README.md:70-94; NEON/NRL_REPLICATION.md:19-35 | ×0.68 (REBUILD 0.1.9), ×0.47 (0.1.36). Obenschain's "about 12 %" KrF intrinsic is a simulation (PAP/obenschain2020.pdf p.1) |

### 4. Other data in the tree
- **Electron density:** Kimura 1987 reproduced in PAP/kushner1988.pdf p.1857 Fig. 3 (XeF mix: Ne / 0.5 % Xe / F2, 3.8 atm, 180 kW/cm³, F2 0.05–0.2 %). Kushner's fit gives a direct F2 dissociation kd ≈ 4.6e-9 cm³/s, about 5× the SIGLO-based 9.3e-10 at 2 eV (KM/PLAN_physical_validity.txt:2682-2685). This tests LIGHT's SIGLO-based f2e group.
- **Discharge KrF (He buffer):** Chow 1978 MEPHISTO, 300 mJ and ~5.5 %, with efficiency vs F2 and energy vs reflectivity (KM/docs/KrF_Chow.pdf PDF p24–28). The legacy code was calibrated to Chow, which is Known-bad. Ong 1996: output vs pressure and charging voltage (PAP/ong1996.pdf p.2592–2593).
- **ArF only:** Electra gain and yields, the Orestes "population benchmark" (model outputs), Chen 2025.
- **Out of scope:** the inverse-bremsstrahlung papers (Turnbull, Langdon, Matte), dudey1967, rauf1999, Boichenko 2000.
- **Planned with no data yet:** the DE-AR microwave experiment (KM/docs/DE-AR_targets.txt).

### 5. Ranking by how well a 0-D model can be scored

| Rank | Dataset | Why it fits 0-D | What blocks LIGHT |
|---|---|---|---|
| 1 | Gadomski/Xu pure Kr/F2 decays | Uniform cell, printed conditions, no electrons or cavity. Tests structure plus G4, G5, K32, K33, C3 and K125 together | Seed mode, n_e ≈ 0, a K_eq decision |
| 2 | Xu He/Ne/Ar series | Same quality | v1 buffer rows |
| 3 | Trainor & Jacob | LIGHT-native e-beam fill, no cavity | Raw data not in tree; e-beam source |
| 4 | Klimek g0/α and pressure exponents | Steady state; the ratios are immune to the unknown pump | v1, e-beam yields, cross-sections |
| 5 | Aurora/Tacoma | Deposition printed | v1; primaries not in tree |
| 6 | Johnson & Hunter g0/α | Conditions printed | v1, e-beam |
| 7 | Kannari Figs 2 and 6 | Shapes | v1, prescribed flux |
| 8 | Kimura n_e | Electron balance | Ne/Xe, e-beam |
| 9 | Mandl/Suda | The only pure-Kr device numbers | Cavity; conditions unknown |
| 10 | Device energies (Kannari Fig. 5 and Table IV, Tisone, Electra, Rhodes) | Cavity knife-edge: error amplified ~1/margin (RB/MODULES/scorecard.py:430-446) | Cavity, geometry |
| 11 | Discharge (Chow, Ong) | none | He, circuit |
| 12 | Kvaran | none | Vibrational levels |

### 6. What REBUILD's scorecard misses or gets most wrong
- **The cleanest data are never scored.** The measured τ_eff of 25 ns (the engine gives 15.2 ns) is recorded but unresolved (RB/CHANGELOG.md:345-352). The validation audit recommends settling it first (VAL/validation_report.tex:569-584, 880-884).
- **No Kr-rich anchor.** The Mandl ratio lives in the PIC bundle, at assumed conditions, and is stale after 0.1.40. The "K34 ~×2 high" hypothesis (RB/CHANGELOG.md:668-674) conflicts with the 300 K measurement of 6.0±2.0e-31 unless temperature explains it.
- **The KrF level is unpinned.** Klimek's pump is borrowed from Kannari's device (excimer_anchors.py:29-38). Gain and absorption are low by the same factor, which points at the pump rather than the kinetics (RB/CHANGELOG.md:295-300).
- **Largest misses:** KrF pressure exponents are wrong at both ends, and the measured 55 A/cm² exponent (+0.235, the survey's arithmetic from Tables II/III) is unscored. The Kannari device scores ×0.23 at 20 % reflectivity with an inverted trend, using a 44-cm length and a "66 J" the paper never prints (scorecard.py:455-457, 491; Kannari p.4315 gives a 40×4 cm² excitation area). Only 2 of the 4 Fig. 5 points are scored. ArF: pressure trend ×1.13 vs ×0.70 measured, level ×0.63 at 0.82 atm (RB/CHANGELOG.md:286-287, 301-303).
- **Unscored measured sets in the tree:** Johnson & Hunter Tables II/III, Kannari Table IV and Figs 2/6/7, Rhodes §4.3, Wolford's efficiencies (VAL/validation_report.tex:350-351).
- **Mislabelled or stale entries (reported, not re-checked for the plan):** KM/SPECK/EXPERIMENTS.md:46 calls an ArF model-output benchmark a "KrF population benchmark" and labels it blind (the code runs Ar/F2 0.995/0.005, scorecard.py:657, and its targets are Orestes outputs, VAL/validation_report.tex:319-322). EXPERIMENTS.md:44-45 and scorecard.py:852-877 still carry pre-0.1.40 verdicts. The "empirical g0/α 4–12" band (scorecard.py:870) traces to a Thompson 1989 compilation cited in KM/PLAN_physical_validity.txt:3122-3125, and that citation "did not survive source verification" (LEG/light_unified/benchmark.py:77-78). In-tree measured g0/α spans 3.7–16.6 (Klimek) and ~4–10 (Johnson & Hunter).
- **Provenance of the rows these data test:** BC_Kr (5.0e-10; network.py:1207-1213, copied by LIGHT) names Gadomski as its primary, but Gadomski *assumed* kBC. The only measured Kr value in the tree is 1.2e-10 (Hakuta, quoted in Xu p.2591 and p.2598), ×0.24 of the row, contradicting the row's note that the class value "likely UNDERestimates". VREL_Kr names Kvaran as its primary, and Kvaran has no Kr data. Minor: excimer_anchors.py:33 and :93 cite the 0.75 MW/cm³·atm calibration as p.4314, where the survey reads it on p.4315.

### 7. Could not check
- Primary sources not in the tree: Klimek 1981 (pump and Kr %), Klimek & Hsia 1981, Trainor & Jacob, Suda 1984/1987, Mandl and Salesky & Kimura (conditions), Kannari/Shaw/O'Neill 1987, Aurora, Tacoma, Electra KrF (Wolford 2013, NAS slides), Kimura's KrF data.
- Figure reads were not digitised: Gadomski Figs 1–2, Kannari Figs 2/5/6/7, Chow, Ong.
