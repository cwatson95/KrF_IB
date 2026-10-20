# ArF (Ar/F2) scoring data: the survey report (2026-09-30)

The report of the research subagent that surveyed the ArF measurements
LIGHT_ArF can be scored against, kept verbatim as returned. Its figure
readings are in `fig_readings/` (one CSV per figure, each with its source
and page in the header). "checked" below is the subagent's own first-hand
reading; the values `score_arf.py` scores against were re-read first-hand
before use (`SCORE.md` says which). The subagent's web and API calls are
listed at the end (all checktext CLEAN); nothing it fetched is relied on
unless read first-hand here.

How to read: "checked" = read first-hand in the PDF (printed page / PDF
page). "fig" = the subagent's own reading of a figure (rendered at 200-300
dpi, axes calibrated on frame and tick pixels, markers located by template
matching; the ± is the reading error, not the authors'). "abstract" =
OpenAlex/Crossref abstract only. "2nd-hand" = quoted by another paper.
Pressures convert to densities at 3.22e16 cm^-3 per Torr (300 K).

## Headlines

1. No direct measurement of the ArF(B) radiative lifetime was found, in the
   library or in web/API searches. Every value in use is a calculation or an
   assumption: 4.2 ns (Hay & Dunning, calc., Rhodes Table 4.2 p.90/PDF 100);
   A = 2.5e8 s^-1 (Dunning & Hay, calc., Rhodes Table 4.14 row 65, p.132/PDF
   142); "4 ns" (assumed by Rokni 1977 and Chen & Payne 1979); 2.5 ns
   (assumed by Trainor & Jacob, Duzy 1980 p.21/PDF 36); the Petrova memo
   contradicts itself (p.5 "A = 2.67 × 10^-8 s^-1 (lifetime τ = 2.7 ns)",
   p.15 2.67e8 s^-1 = "radiative lifetime 3.7 ns"); SRI's "30 nsec"
   (MP-76-99 Table 2 reaction 7) was inferred from yields and is superseded.
   Score against the measured kτ products and decay constants, not against k
   derived with an assumed τ.
2. The cleanest ArF analogue of the KrF photolysis set is Chen & Payne 1979
   (K1): proton-excited Ar/F2 with n_e < 1e12 cm^-3, decays plus three
   Stern-Volmer plots. Next, Ar2F* decay constants from three independent
   groups (K2).
3. The pure-Ar analogue of the Morozov Kr check is Morozov 2008 itself: 0.33
   ± 0.04 above 500 hPa, falling to about 0.23 at 150 hPa (fig). Already in
   ML/ml/pure_data.py as vuv_Ar.
4. Only three groups have gain data on Ar/F2 with no buffer gas: NRL
   Electra, Keio (Lee 1988) and Sandia (Adams & Dillon 1991). Each has a
   pump-rate ambiguity or a composition known only secondhand. There is a
   σ-free target: at Electra (0.82 atm, 1.09 MW/cm³), g0·Isat/P = 0.0321
   cm^-1 × 10 MW/cm² / 1.09 MW/cm³ = 0.29 (arithmetic from printed values),
   the ArF counterpart of Tahoma's 0.27.
5. REBUILD's ArF population check rests on NRL model outputs, not
   measurements. The one "measured" number, ArF* ≈ 4e14 cm^-3 attributed to
   Petrova's ref [58], comes from a He-buffered discharge (section C).
6. No ArF measurement on an Ar/F2-only discharge was found; every discharge
   dataset uses a He or Ne buffer.

## A. Ranked catalogue (best 0-D fit first)

**K1. Chen & Payne 1979: proton-excited Ar/F2, decays and Stern-Volmer
plots.** C.H. Chen & M.G. Payne, IEEE JQE 15, 149 (1979);
`LIGHT/docs/krf_data/sources/Kinetics_of_proton-excited_Ar-F2and_Kr-F2mixtures.pdf`
(printed page = PDF page + 148).
- Conditions (checked, pp.150–151): 2-MeV proton pulses of 10 ns, about
  1e6–1e7 protons per pulse; cell 2×2 in by 3.5 in; flowing Ar (99.9995 %)
  and F2 (99.8 %), F2 from a 1 % F2/Ar premix; room temperature; "electron
  density is always kept low (<10^12 cm^-3)" (p.151).
- (a) Table I p.155 (checked): Ar(3P2)+F2 = 9.4e-10; Ar2*(1u)+F2 =
  5.2e-10 cm³/s; ArF*+F2 = 1.9e-9 ± 25 %; ArF*+2Ar = 5.3e-31 ± 25 % cm⁶/s;
  Ar2F*+F2 = 2.1e-10 ± 25 %; τ(Ar2F*) = 185 ± 25 % ns; τ(ArF*) = 4 ns
  (calculated, their ref. 40); f1/f2 = 0.55 ± 0.25; f2'/f3 < 0.3.
- (b) The measured products (p.155, checked): τ3β3 = 7.5e-18 cm³, τ3β3' =
  2.1e-39 cm⁶; from Fig. 15, τ4β4 = 4.5e-17 cm³.
- (c) p.151 (checked): Ar(3P1), Ar(1P1) and Ar2*(1u) are quenched by F2 at
  8.9e-10, 13e-10 and 5.2e-10; "Ar(3P1) and Ar(1P1) do not produce
  significant amounts of ArF*".
- (d) Fig. 13 p.153 (fig, ±0.1e7 s^-1): ArF* decay frequency (1e7 s^-1)
  vs F2 (Torr). 600 Torr Ar: 0.60→1.81, 1.21→3.86, 2.20→5.87, 3.76→11.49,
  6.49→19.19. 300 Torr Ar: 0.61→1.32, 0.92→2.50, 1.70→4.62, 2.95→8.39.
  Follows the Ar(3P2)+F2 line; slope about 9.2e-10 (subagent's fit).
- (e) Fig. 14 p.154 (fig): Ar2F* decay (1e7 s^-1) vs F2, 300 and 600 Torr
  pooled: 0.76–0.88 Torr → 1.00–1.47; 1.17 → 1.41, 1.51; 1.71 → 1.88;
  2.15 → 1.79; 2.89 → 2.30, 2.40; 3.81 → 3.24; 6.55 → 4.57–4.74. Fit above
  1 Torr: 1.87e-10 and 140 ns (printed 2.1e-10 and 185 ns).
- (f) Figs. 16–17 p.155 (fig): 1/I(193.3 nm) vs F2 at 600 Torr Ar: 1.20→
  15.7, 2.20→18.3, 3.73→22.3, 4.88→23.6, 5.85→25.9 (drawn line intercept
  12.2, slope 2.53 per Torr). Vs P_Ar² (1e5 Torr²) at 0.84 Torr F2: 0.40→
  2.98, 1.11→2.94, 1.59→3.25, 2.57→3.84, 3.72→4.33 (drawn line 2.50 to 5.62
  at 6).
- How to score: seed mode (initial Ar(3P2), and Ar2* for the Ar2F* runs)
  with no electrons at the printed pressures; compare late-time decay
  constants of ArF(B) and Ar2F* emission; for the Stern-Volmer plots,
  compare ratios of time-integrated ArF(B→X) photon counts I(x)/I(ref),
  which need neither the collection constant nor τ.
- Caveats: the printed τ3β3 does not follow cleanly from Fig. 16
  (uncorrected intercept/slope 6.2 Torr → kτ 5.0e-18; drawn line 4.8 Torr →
  6.4e-18; correcting the intercept for the paper's own three-body Ar term
  at 600 Torr, ×1.78, gives 8.9e-18 and 1.15e-17), so score the raw 1/I,
  not 7.5e-18. Below 1 Torr F2 the decays follow the precursors (p.153). The
  ion-ion slow component is about 10 % (p.159).

**K2. Ar2F* decay after a pulse: three groups plus SRI.**
- (i) Jancaitis & Powell, appendix to Ch. V of K.S. Jancaitis, UCRL-53465
  (1983); `Recent_KrF_Papers/janacaitis.pdf` (printed page = PDF page − 8).
  Conditions (checked, pp.134–137): 600 keV, 50 ns, 2–500 A/cm² (Faraday
  cup); excited volume 2.5×1×1.5 cm; 200–2000 Torr; deposition "25 mJ/cm³
  per atmosphere of Ar at 100 A/cm²", ±25 % (SANDYL), about 0.5 MW/cm³ per
  atm per 100 A/cm². Table I p.158 (checked): τ_rad = 205 ± 20 ns; k_Ar <
  1e-14; k_F2 = (1.79 ± .06)e-10; k_NF3 = (1.07 ± .07)e-10; k_SF6 < 1e-13.
  p.150 (checked): kτ for Ar2F* formation from ArF* = 3.1e-39 cm⁶, ±25 %.
  Fig. 7 p.169 (fig, ±0.07e7): decay rate (1e7 s^-1) vs total pressure
  (Torr), at 11 and 33 A/cm²; (a) [Ar]/[F2] = 300: 298→1.09, 444→1.38,
  595→1.63/1.79, 750→1.87, 907→2.43, 1204→2.90, 1510→3.38/3.48,
  1909→3.69/3.81; (b) [Ar]/[F2] = 600: 298→0.80, 444→0.90, 595→1.00,
  750→1.21, 966→1.43, 1204→1.72, 1510→2.00, 1909→2.10/2.34. Fit about 165
  ns and 1.65e-10.
- (ii) Ueda, Kanada, Kitagawa & Takuma, Rev. Laser Eng. 11, 576 (1983);
  Febetron 706, 500 kV, 4.4 kA, 3-ns triangular pulse; Ar/F2 at 1–3 atm.
  Fig. 8 p.581 (fig ±0.1e7): [F2] (1e17 cm^-3) → decay (1e7 s^-1): 1 atm
  0.70→1.91, 1.21→2.85; 1.5 atm 1.04→2.50, 1.86→4.05; 2 atm 1.48→3.39,
  2.45→4.97; 2.5 atm 1.86→3.94, 3.06→6.17; 3 atm 2.17→4.92, 3.66→7.70.
  Printed k = 1.9e-10; Table II p.582: 220 ns. Ueda 1984 Table V (checked):
  220 ns; k_F2 = 1.9 ± 0.1e-10; Ar2F*+Kr = 1.35 ± 0.05e-10.
- (iii) SRI MP-76-99 Fig. 15 (p.58/PDF 62; Table 2 rows 4–5 p.46/PDF 50,
  checked): 290-nm decay vs F2 at 1–9 atm, slope 2.5e-10 (±20 %), intercept
  3.8e6 s^-1; SRI attributed this to Ar2* + F2; the work above shows it is
  the Ar2F* decay itself.
- (iv) Rhodes Table 6.4 (p.191/PDF 199, checked; originals not read):
  Marowsky et al. 219 ± 15 ns, (2.05 ± 0.06)e-10, k_Ar < 1e-14, ArF*+2Ar
  (1.2 ± 0.2)e-30 cm⁶; Böwering et al. (230 ± 25) ns, k_Ar (2.2 ±
  0.8)e-14; Gu Zhi-yu et al. 236 ns, 1.83e-10; Wadt & Hay 132 ns (calc.).
- How to score: fit the post-pulse Ar2F* decay against [F2] and [Ar] and
  compare with the points; for Jancaitis use the ≤33 A/cm² data or the
  post-pulse part.

**K3. Ar* removal by F2 and by Ar under e-beam (absorption).** Ueda 1983
Fig. 12 (p.584/PDF 9): Ar(3P2) 696.5-nm absorption decay in Ar-F2 at 2.5
atm (fig ±0.3e7): F2 % → decay (1e7 s^-1): 0→2.11, 0.049→5.80,
0.098→9.14, 0.148→11.85, 0.200→15.1; printed k = (1.1 ± 0.1)e-9 (text
p.583); ArF* fluorescence gave an apparent 0.9e-9. Ueda 1984 Table V
(checked): 7.4 ± 1.0e-10 by absorption; 6.4 ± 0.3e-10 from ArF*
fluorescence in 0.5 atm with 30 % Kr (not Ar/F2-only); Ar*+2Ar = 6.0 ±
0.5e-33. The two Ueda Ar*+F2 values conflict.

**P1. Morozov 2008: pure Ar VUV efficiency.** JAP 103, 103301 (2008),
`sources/morozov2008.pdf`. Continuous 8–15 keV beams at 5–50 mW;
efficiency independent of power (Fig. 5). Physical efficiency 0.33 ± 0.04
at 500–1400 hPa (abstract and p.103301-7), VUV power into 4π over power
deposited in the gas; apparent 0.31 ± 0.03 at 12 keV (p.-4); deposited
fraction 0.73 (p.-5); foil reflection 0.27 (p.-6). Fig. 10 Ar panel (fig),
hPa → efficiency: 151→0.228; 199→0.26/0.275; 248→0.283; 299→0.306;
400→0.30–0.31; 500→0.313; 600→0.30/0.32; 700→0.31/0.325; 800→0.328;
900–1000→0.32–0.33; 1100→0.326; 1200–1400→0.32–0.34. Low power only.

**P2. Ueda 1983 Fig. 10: Ar(3P2) decay in pure Ar** (Febetron 3-ns pulse,
p.583/PDF 8; fig): atm → decay (1e7 s^-1): 1.46→0.94/1.12, 1.99→1.50,
2.48→2.09/2.24, 3.0→2.94/3.09, 3.4–3.5→3.80, 4.0→4.88. Printed
k(Ar*+2Ar) = (6.0 ± 0.5)e-33 by "parameter fitting" against a simulation.
Recombination repopulates Ar*: score the model's post-pulse Ar(3P2) decay,
not k.

**S1. Rokni, Jacob & Mangano 1977: quasi-steady e-beam fluorescence.**
PRA 16, 2216 (1977), `Recent_KrF_Papers/"Dominant formation and quenching
processes...pdf"`. Conditions (checked, pp.2216–2217): 150 keV; 5 A/cm²
after the foil; 300-ns pulse "enabling the fluorescence amplitude to reach
a steady state"; Teflon cell 22×2×0.2 cm³; F2 98 % pure. Power density not
printed; about 20 kW/cm³ per atm of Ar (subagent's estimate from 5 A/cm²
and a 150-keV Ar stopping power of about 2.4 MeV cm²/g). Table II p.2218
(checked): kτ (k at τ = 4 ns): F2 7.6 ± 0.7e-18 cm³ (1.9e-9); Kr 6.1 ±
0.5e-18 (1.6e-9); Ar 3.6 ± 1e-20 (9e-12); 2Ar → Ar2F* 1.6 ± 0.3e-39 cm⁶
(4e-31); F2 half-pressure 4.0 Torr (p.2218). Fig. 3 p.2219 (fig): 1/I vs
F2 at 150 Torr Ar, 2–16 Torr (the y label "0.8" means 0.08);
intercept/slope 4.78 Torr, which divided by the Ar terms at 150 Torr from
Table II (1 + 0.174 + 0.037 = 1.21) gives 3.95 Torr, the printed 4.0.
Fig. 4 (fig): I vs Ar (atm) at 2 Torr F2, 5 A/cm², relative units:
0.24→103/109; 0.49→194/199; 0.75→236; 1.0→256/264; 1.26→260; 1.5→244/255;
2.0→227/234; 2.5→215; 3.0→205; 3.5→184; 4.0→179. Fig. 5 (fig): ratio of
fluorescence for a ×25 current change: 24.1–26.6 over 0.26–4 atm (1.01→
26.6, 2.0→24.5, 4.0→23.7), read by Rokni as Ar2+ + F- forming ArF* with
branching near unity. Their eq. (1): S = αN_Ar/[1 + (k_F2 N_F2 + k_Ar N_Ar
+ k_2Ar N_Ar²)τ]. Caveats: ionic formation dominates (Fig. 2 is a model);
absolute pump not printed; the Petrova memo's "300 ns e-beam, 1.25 MW/cm3"
for this work is inconsistent with 5 A/cm² at 150 keV.

**S2. SRI MP-76-99 (Nakano, Hill, Lorents, Huestis, McCusker, 1976):
absolute ArF and Ar2F yields.** `sources/SRI_MP-76-99_Nakano_1976_DTIC_ADA033997_text.pdf`.
Conditions (checked): Febetron 706, 500 kV, 2000 A/cm², 3 ns (p.17/PDF 22);
charge density about "10^15 [Ar(atm)] cm^-3"; recombination beats
attachment (p.37/PDF 42). Yield (checked): photons in the band divided by
the "total number of excited states produced initially" (p.25/PDF 30),
calibrated against N2(2+) emission from Ar + 5 % N2 with an assumed yield
of "17.6%/P" (p.21/PDF 26). Fig. 4 (p.26/PDF 31, fig ±0.01), F2 = 1 Torr,
no Kr: 193 nm 1 atm 0.30, 2 atm 0.147, 3 atm 0.124, 5 atm 0.036/0.027, 7 atm
0.023; 290 nm 1 atm 0.42, 2 atm 0.38/0.51, 3 atm 0.41, 5 atm
0.35/0.31/0.18, 7 atm 0.32, 9 atm 0.30. Fig. 11 (p.52/PDF 56, fig ±8 %),
ArF yield vs F2: Ar 1 atm 0.24 (0.5 Torr) → 0.30–0.32 (1–4 Torr) → 0.22 (6
Torr); Ar 2 atm 0.105 → 0.145 → 0.12 → 0.117; Ar 5 atm 0.019–0.033. Text
p.40 (checked): "30% being the greatest observed". Their model agrees
"within a factor of 2" (Fig. 12); their Table 2 ArF* constants (30-ns
radiative, k_F2 1e-9 "Guessed") are superseded. How to score: 3-ns pulse at
the stated ionization density; count 193- and 290-nm photons per initial
excitation (ions plus excited atoms). Caveat: the normalization rests on
their N2 actinometer.

**S3. Jancaitis & Powell: pump-rate scan and excitation-dependent Ar2F*
loss.** Fig. 3 (p.165/PDF 173, fig ±8 %): peak fluorescence (arb.) vs
current density, 5 Torr F2 + 1500 Torr Ar. ArF*: 4.3→0.049, 10.2→0.130,
14.0→0.179, 31→0.37, 74→1.32, 186→3.81, 449→9.0, 562→12.0. Ar2F*:
4.3→0.046, 10.2→0.119, 14.0→0.157, 31→0.34, 73→0.94, 186→1.9, 430→3.1,
540→3.4. ArF* log-log slope about 1.14. Text p.141: electron quenching
"does not cause a significant reduction in the dimer intensity" up to 500
A/cm². Fig. 12 (p.174/PDF 182, fig): extra Ar2F* loss vs current density
at 2.5, 5.0 and 7.5 Torr F2 in 1500 Torr Ar, about 1.3e8 s^-1 at about 550
A/cm² for 2.5 and 5 Torr and about 0.95e8 for 7.5 Torr; linear in J.
Estimated Ar2F*+e about 4e-8 cm³/s (p.152, model-dependent). How to score:
e-beam source of 0.04–5 MW/cm³ over 50 ns (from the printed deposition).

**S4. Trainor & Jacob 1980: ArF* + e.** APL 37, 675 (1980) (abstract, doi
10.1063/1.92047); data in Duzy et al., AVCO 1980, Table III p.23/PDF 38
(checked): Ar/F2 at 0.25 atm, k_el = 1.6e-7 cm³/s at 44.4 A/cm², 3.0e-7 at
10.1 A/cm²; the reduction assumed τ_rad = 2.5 ns (p.21); k_el ∝ 1/τ_rad,
so 4.2 ns gives 0.95e-7 and 1.8e-7. Raw Stern-Volmer data not printed.
Lee 1989's 2.4e-7 for Ar/F2 is a fit to their amplifier data (p.2056).

**G1. NRL Electra: gain (99.7 % Ar / 0.3 % F2).** Wolford 2020 (HEDP 36,
100801, `Recent_KrF_Papers/wolford2020.pdf` pp.2–4); NRL/6730/MR--2023/1
(`ArF_report.pdf` p.8/PDF 15 and Fig. 5 p.9/PDF 16); Obenschain 2020 (Phil.
Trans. R. Soc. A, p.7); Petrova memo Fig. 3 p.9/PDF 16 and Fig. 7 p.14/PDF
21 (all checked). Machine: 500 keV, 110 kA total, 100-ns flat-top in a
140-ns pulse; two counter-propagating beams; 25.4-µm stainless foil; 1.4 kG
guide field. Deposition measured by pressure rise; deposition volume imaged
as 14×30×100 cm³. Gain measured "at temporal peak power deposition"
(Obenschain p.7). Values (checked): single-pass, 30-cm cathode: 1.84 %/cm
at 12 psi (0.82 atm) and 0.395 MW/cc; 1.77 %/cm at 16 psi and 0.497 MW/cc
(the report prints "kW/cc"); 1.29 %/cm at 20 psi and 0.583 MW/cc.
Double-pass, 30-cm cathode: 1.37 and 1.10 %/cm at 12 and 20 psi (called low
because of the analysis). Double-pass, 10-cm cathode, 12 psi: g0 = 3.21
%/cm, α = 0.16 %/cm, Isat = 10 MW/cm² at 1.09 MW/cc (Isat "could be
estimated"). Intrinsic efficiency 17.7 % (Wolford p.2) or 16 %
(Obenschain). Fig. 5 adds 1.36 atm, 10-cm cathode, 2.90 %/cm at 1.38 MW/cc.
Measured pump waveform (Petrova Fig. 3, fig): peak 504 kW/cm³ at about 135
ns; integral 67.9 mJ/cm³. How to score: drive with the waveform scaled to
each power; compute g0 at peak deposition; prefer the σ-free g0·Isat/P =
0.29 (0.36 if the pump was 0.90 MW/cc) and g0/α ≈ 20; the drop of g0 with
pressure at rising pump (1.84 → 1.29) is the discriminator. Caveats: g0
needs σ_se (NRL uses an effective 1e-16 cm², memo p.16; Lee 2.75e-16); the
beam is non-uniform at low pressure (Wolford p.4); the 3.21 %/cm point is
plotted at 0.90 MW/cc in Petrova Fig. 7 but at 1.10 in the ArF report; the
seed timing is not given.

**G2. Keio: Lee et al., APL 52, 1294 (1988) and IEEE JQE 25, 2053 (1989).**
The APL paper as an abstract (doi 10.1063/1.99142) and
`docs/Recent_KrF_Papers/lee1989.pdf` (running header IEEE Trans. Electron
Devices 36; Crossref gives IEEE JQE 25, doi 10.1109/3.35232). Single pass,
L = 50 cm; e-beam 1.3 MV, 90 kA, 65 ns FWHM; Ar/F2 = 99.6/0.4 at 1 atm
(1989 Fig. 3 caption p.2057). Power efficiency 15.7 % and extraction
efficiency 59.4 % for 99.6 % Ar (APL abstract), inputs 6.7–9.5 MW/cm².
1989 Fig. 5(a) p.2058 (fig): g0 ≈ 4.18 %/cm, α ≈ 0.232 %/cm, Is ≈ 3.2
MW/cm². Text p.2053: Is "1.9–3.2 MW/cm2", "g0/α (>18)", excitation rate
"0.51 MW/cm3". Caveats: the pump is ambiguous (text 0.51 MW/cm³; NRL plots
the same point at 1.21, ArF report Fig. 5, and 0.69, Petrova Fig. 7); the
1989 Fig. 4 model pulse peaks near 1.05 MW/cm³ for the 70 % Ar mix; the
probe arrived 30 ns after the sidelight peak.

**G3. Sandia: Adams & Dillon, JAP 70, 4073 (1991).** Abstract only (doi
10.1063/1.349157; OSTI 5167262): longitudinal pumping at "∼5-MW/cm3";
963-Torr mix; 15-ns probe. At 193.31 nm: g0 = 12.1 %/cm, α ≈ 0.70 %/cm, Is
= 6.7 MW/cm²; at 193.68 nm: 10.8, 0.70, 7.9. Composition (953 Torr Ar + 8
Torr F2) 2nd-hand, Petrova p.13. g0·Is/P ≈ 0.16–0.17.

**B1. Avco: F2 burn-up.** Duzy 1980 pp.14–16/PDF 29–31 (checked):
return-current traces for Ar/F2 0.2 % and 0.4 % at 2 atm (Fig. 5), current
rising to about 6.5 A/cm² over about 600 ns. Burn-up numbers printed only
for 0.4 % F2 / 4 % Kr / 95.6 % Ar (Fig. 6, fig): about 4.9e17 F2 per J/cm³
(e.g. 9.9e16 cm^-3 at 0.185 J/cm³), about 12.6 eV per F2. The text says the
Ar/F2 burn-up is "the same" within scatter; the Ar/F2 numbers are in AERL
quarterlies, not obtained.

**O1. Oscillator and device energies (need a cavity).** Electra
quadruple-pass, Wolford Table 1 p.4 (checked): 80 J at 1.36 atm / 0.9 MW/cc
(12 %); 96 J at 1.09 / 0.8 (14 %); 90 J at 0.82 / 0.6 (16 %). Also 137 J
dual oscillator at 99.71/0.29, 0.95 atm (p.5); 200 J at 0.76 MW/cm³, 13 psi,
0.5 % F2 (memo p.20); oscillator/amplifier pairs at 0.5 % F2, 0.884 atm, for
8/16/23 % output couplers: 57/95, 54/126, 51/100 J (ArF report p.5). Keio
(abstracts): 95 J and 3.4 % from 650 Torr Ar/F2 at 2.3 MW/cm³ with 65 ns
(Suda, Obara & Noguchi, JAP 60, 3791 (1986)); 96 J (5 J/ℓ) and 3.4 % (Suda
et al., JAP 64, 1720 (1988)). Sandia 92 J (Hoffman, Hays & Tisone, APL 28,
538 (1976), abstract; mixture not in the abstract).

## B. Rate-constant measurements (inputs or single-process checks)

- Ar*+F2 total (cm³/s): 7.5e-10 (Velazco et al. 1976); 9.4e-10 (Chen &
  Payne); 1.1e-9 (Ueda 1983); 7.4e-10 (Ueda 1984).
- ArF* branching from Ar*+F2: "∼0.5 rather than unity" (Setser,
  COO-2807-11, p.6/PDF 8, checked); Kolts & Setser, J. Phys. Chem. 82, 1766
  (1978): 4.0e-10, 2nd-hand via Ueda 1984 Table III (0.53 of 7.5e-10); "≤0.6"
  (SRI, via Rhodes Table 4.4 p.98); "about 0.6" (Rokni 1978 review, citing
  unpublished Chen/Judish/Payne); resonance states give little ArF* (Chen &
  Payne). The ArF C/B ratio against pressure is in Kolts & Setser 1978, not
  read (paywalled).

## C. Problems in the repo's current ArF anchors

1. REBUILD `electra_population_rows` (scorecard.py): Te 2 eV, Tg 420 K, ne
   5e14, Ar* 1e15, F- 1e14, "F2 50 % dissociated", g0 2.5 %/cm and α 0.2
   %/cm are Orestes model outputs (memo Figs. 4, 6, 8), not measurements. The
   ArF* "measured ~4e14" is the memo's citation of ref [58], Nagai et al.
   (the memo prints IEEE JQE 34, 40 (1998); Crossref gives 34, 942–948, doi
   10.1109/3.678589), whose abstract describes a discharge at F2-Ar-He =
   0.2-10.0-89.8, 2.5 atm, measuring electron density (peak 6.7e15); the memo
   instead says "Ar/F2/He = 89.5%/0.5%/10%". An ArF* density is not in the
   abstract.
2. Memo p.13: "measured small-signal gain on Electra SPA was 12.1 %/cm, α
   ∼0.70 %/cm, Isat 6.6 MW/cm²" are Adams & Dillon's Sandia numbers (12.1 /
   0.70 / 6.7), not Electra.
3. Memo Fig. 7 shows an "Electra" point at 3.99 %/cm and 1.22 MW/cc; it
   coincides with the ArF report's placement of Lee's point, and the 1.36-atm
   10-cm-cathode point (2.90 %/cm) is missing.
4. REBUILD `electra_arf_rows`: the (atm, MW/cc, J) triples match Wolford
   Table 1 as printed.

## D. Not found or not read

- No ArF(B) lifetime measurement.
- No measured ArF*/Ar* densities in e-beam Ar/F2 (NRL says it measured
  species but reports no numbers, ArF report §2.4).
- No Ar/F2-only discharge data; the discharge sets use buffers (Nagai 1998
  He; Feenstra thesis, Twente 1999: 6 mbar F2 / 100 mbar Ar / 200 mbar Ne /
  3.7 bar He; chen2025 Ne/Ar/F2/Xe 96.4/3.5/0.1/0.01; Mochizuki 1989,
  composition not read).
- Buffered e-beam ArF lasers, excluded: Edwards 1980 (Ne/He/Ar/F2); Mandl
  1986 (Ne); Peters 1987 (Ne); Glessner 1991 (lased only with Ne); Murray &
  Powell 1976 (UCRL-78484; qualitative only, "best performance using neon").
- Abstract only: Chang, JCP 78, 2210 (1983), doi 10.1063/1.445064, measured
  ArF(B) fluorescence yield, formation yield and quenching by Ne, Ar, F2 and
  NF3 under a long-pulse e-beam. Worth obtaining.
- WebSearch listings seen but not opened and not relied on:
  https://osti.gov/biblio/5210436-ar-sub-radiative-lifetime-measurement ;
  https://pubs.acs.org/doi/10.1021/j100504a021 ;
  https://www.academia.edu/104506262/Rate_constants_for_argon_fluoride_ArF_formation_from_reactions_of_argon_3P2_0_with_fluorine_containing_molecules_and_the_pressure_dependence_of_the_C_to_B_state_ratios_for_argon_fluoride_ArF_krypton_fluoride_KrF_and_xenon_fluoride_XeF_ ;
  https://opg.optica.org/ol/abstract.cfm?uri=ol-9-12-539 ;
  https://www.researchgate.net/publication/231019493_Fluorescence_of_high-pressure_argon_excited_by_an_energetic_flash_X-ray_source

## E. Files and URLs

- The figure readings: `fig_readings/` (16 CSVs, copied from the
  subagent's scratch directory).
- URLs the subagent fetched with curl, all checktext CLEAN:
  api.openalex.org/works (DOI filters and title searches, 5 calls, one
  rate-limited); api.crossref.org/works?query.bibliographic=… (4 calls);
  https://www.osti.gov/api/v1/records?q=%22ArF%22%20laser%20electron%20beam%20gain… ;
  https://www.osti.gov/servlets/purl/7334373 (Murray & Powell);
  https://ris.utwente.nl/ws/files/284073806/thesis_L_Feenstra.pdf. 2
  WebSearch calls.
