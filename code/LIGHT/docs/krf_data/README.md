# Measured Kr/F2 data for LIGHT (2026-09-27)

This folder holds the explicit measurements on the Kr/F2 gas that LIGHT is
checked against, and the source documents they come from. It follows the
owner direction of 2026-09-27. Keep a few basic processes, pin each one to
KrF data, and perfect the Kr/F2 gas before adding anything else.

The data score the model and never steer it (the repo's predictive-only
policy). A mismatch is a reason to look for a missing or wrong process with
its own measurement. It is never a reason to move a rate toward the data.

**checked** means the number was re-read first-hand in the source for this
file. **agent** means a search agent quoted it and nobody has re-read it.

## Where LIGHT stands

| Test | Measured | LIGHT basic network | Where it runs |
|---|---|---|---|
| KrF(B,C) decay, 50 Torr Kr, 1 Torr F2 | 5.0e7 s⁻¹ (Xu 1993 p.2594) | ×0.97 (×0.89 with Gadomski's own constants) | `v0/tests/test_model.py` |
| KrF(B,C) decay, 100 Torr Kr, 1 Torr F2 | 6.2e7 s⁻¹ (Xu 1993 p.2594) | ×1.01 | `v0/tests/test_model.py` |
| KrF(B,C) decay, 500 Torr Kr | 2.63e8 s⁻¹ (read off Gadomski 1992's fitted curve, not re-checked) | ×1.01 | by hand with the same test helper |
| KrF(B) alone, decay against Kr (Eden 1978) | Σ = 1/6.8 ns + 8.6e-12[Kr] + 9.7e-31[Kr]² | 26 % rms, a known conflict | strict xfail test, `bc_check.py` |
| KrF(C) decay against Kr and F2 (Hakuta 1990) | Figs. 2-3, `digitized/` | 29 % rms, rates about 1.8× Kansas State's | `bc_check.py` |
| Tahoma small-signal gain g0 (parked) | 6.39 %/cm | 3.67 %/cm (×0.57) | `tahoma_check.py` |
| Tahoma g0, with Gadomski's 300 cm⁻¹ B-C gap (parked) | 6.39 %/cm | 4.40 %/cm (×0.69) | `tahoma_check.py` |
| Tahoma formation efficiency | 0.27 | 0.341 (×1.26), an approximate comparison | `tahoma_check.py` |

The Tahoma rows use a prototype e-beam drive (below) and are parked by
owner direction (kinetics first, no drive modules in LIGHT). The decay rows
are exact runs of the basic network.

`BC_AND_ELECTRON_QUENCH.md` (2026-09-27, later) covers the KrF(B), KrF(C)
and electron-quench measurements in detail:
- Eden 1978's B-only decays, which LIGHT misses by 25 % rms
- the B-to-C mixing rate, which is assumed and open by a factor of 2 to 10
- the Trainor & Jacob product, which LIGHT's 4.5e-8 matches once it is
  read with the two-level lifetime.

## 1. Tahoma, the best independent test (checked)

Spectra Technology, *Experimental Support for KrF Laser Modeling*, final
report to Los Alamos, 19 July 1989, OSTI 10136715
(`sources/Tahoma_SpectraTech_1989_OSTI10136715.pdf`).

Conditions:
- 99.6 % Kr/F2 with 2.8 Torr F2 and 665 Torr total (Table 3-3 footnote)
- e-beam, 450 ns pulse (p. 2-2)
- 330 kW/cm³ deposited at the centre of the active volume, where gain was
  probed (p. 2-3). The volume average is 285 kW/cm³.
- The beam stop raises the pumping by about 20 % (p. 3-7). The 330 kW/cm³
  figure belongs to the beam-stop-in runs.
- Gas temperature is not stated.

Measured at 248.4 nm, in %/cm (Table 3-3):

| | g0 − α0 | α0 | g0 |
|---|---|---|---|
| beam stop out (lower pumping) | 5.05 | −1.25 | 6.30 |
| beam stop in (330 kW/cm³) | 5.28 | −1.11 | 6.39 |

- Saturation intensity Isat = 1.4 ± 0.1 MW/cm² (p. 3-16).
- Formation efficiency g0·Isat/Pin = 0.27 (p. 3-16). Their arithmetic is
  g0·Isat = 0.064 cm⁻¹ × 1.4 MW/cm² = 90 kW/cm³, over Pin = 330 kW/cm³.
- Intrinsic laser efficiency 9.9 % (p. 4-1). This needs a cavity model.
- Pure Kr with no F2, 257.3 nm: −0.72 %/cm (stop out) and −0.67 %/cm (stop
  in) (Table 3-2). This isolates the absorbers made from Kr alone.
- Gain uncertainty ±0.1 %/cm, loss ±0.06 %/cm (Table 3-2).

What LIGHT must compute to use it:
- g0 = σ_B·[KrF(B)] needs no cavity and compares directly. LIGHT uses
  σ_B = 2.6e-16 cm² (Tellinghuisen 1976 via Kannari 1985 Table I, as in
  REBUILD `mode_optics`). Only the B state lases.
- α0 needs the absorber densities times their cross sections (F⁻, Kr2⁺,
  Kr2F*, F2, excited Kr).
- Tahoma's formation efficiency is g0·Isat/Pin. The model's is the KrF
  formation rate times 5.0 eV over Pin. The two agree only for a simple
  two-level medium. With B and C mixing, a like-for-like number needs the
  model's own Isat.

Trap: pair 6.39 %/cm with 330 kW/cm³. The beam-stop-out gain (6.30 %/cm)
was taken at about 20 % lower pumping. The legacy benchmarks
(`light_unified/benchmark.py`, `unified/validation_ebeam.py`) pair 5.45 %/cm
with 285 kW/cm³. In Table 3-3, 5.45 is g0 − α0 for the 14.5 % Kr Ar-diluent
mix with the beam stop out.

## 2. KrF(B,C) decays after photolysis (checked)

- Xu, Gadomski & Setser 1993, p.2594: 5.0e7 and 6.2e7 s⁻¹ at 50 and 100 Torr
  Kr with 1 Torr F2. Their p.2591 sets the coupled lifetime at 25 ± 2 ns, so
  K_eq "must be 4-5".
- Gadomski, Xu, Setser & Nelson 1992, eq. 5: the two-level B/C model whose
  constants the basic network uses (τ_B 7.0 ns, τ_C 75 ns, K_eq 5.0, k_F2
  3.0e-10, k_Kr 4.0e-12, k'_Kr 6.0e-31). Their p.156 puts the B-C gap near
  300 cm⁻¹ (">200").
- Both are in `docs/Recent_KrF_Papers/` of this repo, not in `sources/`.

## 3. Eden et al. 1978 decays (checked in part)

J.G. Eden, R.W. Waynant, S.K. Searles, R. Burnham, APL 32, 733 (1978),
reprinted in the NRL annual report DTIC ADA072359
(`sources/NRL_1978_ADA072359_Eden_KrF_section.pdf`, extract p.9 Table I,
p.13 the reprint).

- KrF2 photolysis in 20-200 Torr Kr. τ(B-X) = 6.8 ± 0.2 ns (agent).
- Kr 8.6 ± 0.7e-12 cm³/s and 2Kr 9.7 ± 0.8e-31 cm⁶/s (checked). F2 4.8 ±
  0.3e-10 (checked).
- Their decay constant at 100 Torr (about 1.85e8 s⁻¹, agent's figure reading)
  is about 3× Xu's 6.2e7. Gadomski 1992 explains it: "the photolysis of KrF2
  initially gives only KrF(B)" (checked). A B-only start sees the fast B-to-C
  transfer. So this set tests the B/C mixing rate, and LIGHT can run it in
  seed mode with KrF(B) alone at t = 0.

## 4. Electron quenching of KrF (checked)

Trainor & Jacob in Duzy et al., *Rare Gas Halide Lasers for Fusion*, AVCO
final report 1980, OSTI 5268869 (`sources/Duzy_AVCO_RGHlasers_1980_OSTI5268869.pdf`).

- Kr/F2 binary, e-beam, 190 Torr Kr, F2 0.5-4 Torr.
- Electron quenching 2.0 ± 0.2e-7 cm³/s, reduced with τ_rad = 6.5 ns and
  k_att = 4.5e-9 cm³/s (checked).
- LIGHT uses 4.5e-8 (Kushner 1989's reassessment of the same data, REBUILD
  K67). The printed value is 4.4× higher. At Tahoma's n_e near 5e14 cm⁻³ the
  two give 2e7 and 1e8 s⁻¹. The model's total KrF loss there is 2.3e8 s⁻¹,
  so the choice matters for g0.
- Nishioka et al., Rev. Laser Eng. 16, 571 (1988) is a later e-beam
  measurement of the same rate (`sources/Nishioka_RevLaserEng16_571_1988.pdf`,
  not yet read).

## 5. Rate constants measured on Kr/F2

| Process | Measured | Source | LIGHT basic |
|---|---|---|---|
| Kr(³P₂) + F2 | 9.4 ± 1.0e-10 (absorption), 8.8 ± 0.5e-10 (fluorescence) (checked) | Ueda et al., Rev. Laser Eng. 12, 357 (1984) Table V | 8.1e-10 |
| Kr(³P₂) + F2 | 7.2e-10, KrF branching near 1 (agent) | Velazco, Kolts & Setser preprint in `sources/Setser_ERDA_COO-2807-2_1976_OSTI7361313.pdf` pp.34-35 | |
| Kr* + 2Kr → Kr2* | 1.7 ± 0.2e-32 cm⁶/s (agent) | Ueda 1984 | |
| KrF(B) radiative lifetime | 6.8-7 ns (Eden 1978, NRL) or 9-10 ns (Burnham & Searles 1977, Quigley & Hughes 1978, Tiee 1989) (agent) | `sources/NRL_MR3753_1978_DTIC_ADA055912.pdf`, `sources/Tiee_LA-UR-89-2337_OSTI5844798.pdf` | 7.0 ns |
| B↔C mixing by Kr | 1.2 ± 0.3e-10, τ_C 60 ± 10 ns (abstract only) | Hakuta et al., JJAP 29, 1194 (1990) | 1.0e-10, 75 ns |
| Kr2F* lifetime, F2 quench | 210 ± 20 ns, 1.8 ± 0.1e-10 (agent) | Ueda et al., Rev. Laser Eng. 11, 576 (1983) | Q2F 1.2e-10 (Xu Table I) |

None of the measured KrF(B) lifetimes supports the 4 ns in Gary's notes.

## 6. Cross-checks that need a cavity or are abstract-only (agent)

- Ueda, Rev. Laser Eng. 16, 209 (1988), Fig. 9: intrinsic efficiency against
  pressure for 99.5/0.5 Kr/F2, 9.5 % at 1 atm. Needs a cavity model and the
  deposition at each pressure.
- Kannari, Shaw & O'Neill, JAP 61, 476 (1987): formation efficiency 24 % at
  99.7 % Kr and 1.4 MW/cm³ (abstract). A direct 0-D check once found in full.
- Salesky & Kimura, APL 46, 927 (1985): g0 − α = 3.2 %/cm at 99.6 % Kr near
  120 kW/cm³ (abstract).
- Suda et al. (1984, 1987), Mandl, Klimek & Salesky (CLEO 1985): abstracts or
  no copy found.

## 7. KrF formation (2026-09-28, ML)

The formation chain's measurements are registered link by link in
`ML/ml/formation_data.py`, with the analysis in ML/README.md ("Results for
KrF formation"). What changed here:
- The Tahoma formation row above compares an approximate quantity. Like for
  like (the model's own g0·Isat/Pin with a probe flux), LIGHT gives 0.328
  (×1.22) for every B/C set. The model's Isat is a second test: 1.49 MW/cm²
  with ML's fitted B/C set against the measured 1.4 ± 0.1.
- Kr* + F2 forms KrF* partly in C: k_B/k_C ≈ 1.6 (Yu & Setser 1990,
  checked). The "unity" KrF* branching of Velazco, Kolts & Setser (1976) is
  relative to Kr* + Cl2 taken as 1, and that reference was measured at 0.7
  ± 50 % (COO-2807-3 p.21, checked).
- Kr* + 2Kr → Kr2*: 1.7 ± 0.2e-32 (Ueda 1984, checked) against 2.6-4.4e-32
  from low-pressure and discharge work (abstracts). Kr2* lives 250-353 ns
  in six measurements, not the 12.5 ns of Kannari's table.
- No measurement was found of the Kr⁺/Kr₂⁺ + F⁻ recombination rate near 1
  atm, or of its KrF* yield or B/C split.

## Findings that touch the rest of the repo

- **The Klimek grid has no measured pump (checked).** Duzy 1980: "it was not
  possible to make direct absolute energy deposition measurements ... the
  e-beam stopping power used in the computer code was adjusted to match the
  gain data at 20 and 55 A/cm²." REBUILD scores against that grid (4 % Kr in
  Ar).
- The legacy Tahoma benchmark pairing (section 1 trap).

## The prototype e-beam drive

Tahoma, like nearly all KrF data, was pumped by an electron beam, not by
heating the electrons. `tahoma_check.py` borrows the prototype source from
`../survey_2026-09-27/scripts/ebeam_screen.py`. It uses REBUILD's
Spencer-Fano yields for Kr with 0.5 % F2. Per ion pair (27.07 eV) that is
Kr⁺ 0.999, Kr** 0.512 (cascading to Kr*), Kr* 0.263, and 4.5 eV into the
slow electrons.

Two results from it:
- Recombination decides the answer under e-beam pumping. Without the `ions`
  and `f2e` groups, n_e runs to 2e16 cm⁻³ and g0 falls to 0.55 %/cm.
- K_eq(Tg) by detailed balance on the printed gap raises g0 from 3.60 to
  4.20 (200 cm⁻¹) and 4.53 %/cm (300 cm⁻¹) at the computed Tg near 460 K.
  The formation efficiency does not move.

Where the KrF goes at the Tahoma state (basic + ions + f2e, 200-400 ns,
Tg 462 K, n_e 4.8e14 cm⁻³). The total loss is 2.3e8 s⁻¹, a 4.35 ns lifetime
of the B/C manifold.

| Loss | Share |
|---|---|
| two-body quench by Kr (B and C) | 37 % |
| trimer Kr2F* (B and C) | 29 % |
| radiation (248 nm and C-A) | 19 % |
| electrons | 9 % |
| F2 | 6 % |

Only 23 % of the manifold sits in the lasing B state (C/B = 3.4). By the
report's own relation, g0·Isat/Pin is the formation efficiency, so the
measured g0 and Isat imply B holds 2.2 ns of population per KrF formed
(hν/(σ_B·Isat)). The model holds 0.99 ns. The model already makes more KrF
than measured, so with σ_B = 2.6e-16 cm² the gain shortfall is in how long
KrF stays in B, not in how much KrF is made.

## Sources

All files passed `~/text-sandbox/checktext` (CLEAN) on 2026-09-27 before
they were opened. The NRL annual report is cut to its cover and the
KrF/HgX kinetics section (original PDF pages 1 and 6-30), from 52 MB.

| File | What it is | Used for |
|---|---|---|
| `Tahoma_SpectraTech_1989_OSTI10136715.pdf` | Spectra Technology final report to LANL, 1989 | Section 1 |
| `NRL_1978_ADA072359_Eden_KrF_section.pdf` | NRL annual report, Eden et al. reprints | Section 3 |
| `NRL_MR3753_1978_DTIC_ADA055912.pdf` | NRL MR 3753, incl. the Burnham & Searles preprint (PDF p.61) | Lifetimes |
| `Duzy_AVCO_RGHlasers_1980_OSTI5268869.pdf` | AVCO final report 1980 | Section 4, the Klimek grid |
| `Tiee_LA-UR-89-2337_OSTI5844798.pdf` | Tiee et al. 1989, KrF(B) lifetime 9.9 ns | Lifetimes |
| `Setser_ERDA_COO-2807-2_1976_OSTI7361313.pdf` | Setser ERDA report, Velazco preprint | Kr* + F2 |
| `Ueda_RevLaserEng12_357_1984.pdf` | Pulse radiolysis of Kr/F2 (II) | Kr* + F2, Kr* + 2Kr |
| `Ueda_RevLaserEng11_576_1983.pdf` | Pulse radiolysis, de-excitation rates | Kr2F* |
| `Ueda_RevLaserEng16_209_1988.pdf` | E-beam KrF for fusion drivers | Section 6 |
| `Ueda_JPlasmaFusionRes58_464_1987.pdf` | E-beam KrF review | Section 6 |
| `Nishioka_RevLaserEng16_571_1988.pdf` | E-beam measurement of KrF electron quenching | Section 4, not yet read |
| `Takuma_Ueda_JPlasmaFusionRes53_276_1985.pdf` | High-efficiency e-beam KrF | Not yet read |
| `McCown_Greene_LA-UR-90-291.pdf` | Kr2F absorption 5.1e-18 cm² at 308 nm (the PDF lacks its Table 1) | Absorbers |
| `Hakuta_1990_Jpn._J._Appl._Phys._29_1194.pdf` | Hakuta et al., JJAP 29, 1194 (1990), added by Cooper 2026-09-27 | k_BC, τ_C in LIGHT |
| `jenkin1990.pdf` | Jenkin et al., J. Phys. Chem. 94, 2927 (1990), the article before Yu & Setser (added by mistake). Its last page carries Yu & Setser's title and abstract | none |
| `Setser_ERDA_COO-2807-3_1976_OSTI7334418.pdf` | Setser ERDA report, July 1976 | Kr* + Cl2 reference 0.7 ± 50 % (p.21) |
| `Setser_ERDA_COO-2807-11_1977_OSTI6681231.pdf` | Setser ERDA annual report 1977 (scan) | ArF* from F2 ~0.5 (pp.5-6) |
| `Setser_ERDA_COO-2807-19_1979_OSTI5572256.pdf` | Setser ERDA final report 1979 (scan) | KrF B/C against pressure (pp.4-9) |
| `Avco_OneMeterKrF_Feb1977_DTIC_ADA047684.pdf` | Avco Everett, One meter KrF laser system, Feb 1977 | "essentially unit branching" (p.31, Ar-rich) |
| `Yu_Setser_JPC94_2934_1990.pdf` | Yu & Setser, J. Phys. Chem. 94, 2934-2942 (1990), added by Cooper 2026-09-27 as `Kinetic_Modeling/yu1990.pdf` (checktext CLEAN) | I_r(300 K) = 0.16 ± 0.02, K_eq 1.7 ± 0.2, ΔE_BC −120 ± 50 cm⁻¹, slow B↔C mixing |
| `Svehla_1962_NASA_TR_R-132.pdf` | Svehla, Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures, NASA TR R-132 (1962), NTRS 19630012982, added 2026-10-01 (checktext CLEAN) | 1Dv0_3's moving gas, the Lennard-Jones force constants of Kr, F2 and F (`1Dv0_3/SOURCES.md` section 7) |
| `Brokaw_1961_NASA_TR_R-81.pdf` | Brokaw, Alignment Charts for Transport Properties, NASA TR R-81 (1961), NTRS 19980232907, added 2026-10-01 (checktext CLEAN) | 1Dv0_3's moving gas, the diffusion coefficient and its combining rules |
