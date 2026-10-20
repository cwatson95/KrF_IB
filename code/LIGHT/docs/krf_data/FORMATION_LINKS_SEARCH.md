# KrF formation: the search for the unmeasured links (2026-09-28)

Step 4 of the owner's plan (2026-09-28, "Proceed with all 4 steps"). The
formation chain (ML/ml/formation.py) had three links with no measurement:

1. the absolute KrF* yield of Kr* + F2,
2. the KrF* yield and B/C split of Kr+ and Kr2+ + F- recombination,
3. Kr2* + F -> KrF* + Kr.

A fourth question was what Morozov et al.'s e-beam VUV efficiency counts
(ML's pure-gas check).

A research agent searched on 2026-09-28. It used 3 web searches and free
APIs, with no email address in any request. The numbers that change LIGHT
were then re-read first-hand in the files kept in `sources/`. The registry
entries are in `ML/ml/formation_data.py` (VALUES, ANALOGS, ABSENT) and
`ML/ml/pure_data.py` (MOROZOV_COUNTS).

## Summary

| Question | Result | What changed |
|---|---|---|
| Absolute KrF* yield of Kr* + F2 | No measurement in the open sources. The published "unity" is relative to Kr* + Cl2, and Setser's promised recalibration never appears in his later reports (COO-2807-4 to -18). **Found later: Velazco, Kolts & Setser 1976 measured the KrF* formation rate constant absolutely, 65e-11 of the 72e-11 total, so the yield is 0.90, published as unity +-0.20** (below) | Analogs registered: ArF* 0.5 to 0.6, Kr* + SF6 about 0.31 per Kr*. LIGHT uses 0.90 since the papers came |
| Kr* + F2 rate | LIGHT's 8.1e-10 cannot be traced. The report Kannari and Ueda credit it to gives 7.2e-10 | The basic network uses the first-hand mean, 8.47e-10. v0 keeps REBUILD's 8.1e-10 |
| Ion-ion KrF* yield and B/C split | No measurement with F-. Analogs: B:C 0.62:0.38 with SF6-, 0.60:0.40 with Cl-. KrCl* is 0.64 of the Kr+ + Cl- products | Brackets now use the analogs (C share 0.38, yield 0.64) |
| Kr2* + F -> KrF* + Kr | No measurement exists. Kannari's 3e-10 cites a report with no such rate | Citation corrected. The row stays an estimate |
| Kr2* + F2 | A second measurement: 3e-10 +- 30 %, product assigned to Kr2F* (SRI) | Registered. It agrees with Quigley & Hughes' 2.1e-10 |
| Morozov 2008 | Second-continuum light over the energy deposited in the gas, PTB-traceable to +-14 % (k = 2) | Registered. The integration limits and backscatter handling are still open |

Status labels: checked = read first-hand in the full text (OCR text layer or
page scan), abstract = abstract only, secondhand = quoted by another document
(named). Only measurements of the process itself count. Analogs inform the
brackets and never set a link.

## 1. The absolute KrF* yield of Kr* + F2

**Setser's branching ratio.** Measured as 0.7, published as "unity" in
anticipation of a recalibration (COO-2807-2, p.17, checked, already in the
registry):

> "Although the Kr(3P2) + F2 branching ratio was measured as 0.7,
> recalibration measurements from 250-190 nm (which are in progr[ess] with
> the D2 lamp) strongly suggest that the final ratios between the KrF and
> KrCl intensities will lead to unity branching ratio for KrF formation
> from F2. This anticipated result was used in Table 5."

The Kr* + Cl2 reference itself "may have an uncertainty of as much as + 50%"
(COO-2807-3, p.21). None of the later reports gives the outcome. COO-2807-4
(July 1976) and -8 (July 1977) were read by text layer, and -12 to -18
(1978-79) by OCR. COO-2807-16 (March 1979): "Work was concentrated on the
Kr + Clo [sic] reaction, Though final results have not yet been obtained".
COO-2807-18, p.6, on 3P1 against 3P2: "the branching fractions for the
metastable and resonance states must be quite similar for a given halogen
donor. Although this has not been proven, we are fairly certain of this
result."

**The Kr* + F2 rate (checked, re-read 2026-09-28).** Ueda 1984 Table II and
Kannari 1985 Table I both credit 8.1e-10 to Nakano et al., SRI MP-76-99
(DTIC AD-A033997). That report lists reaction (14), Kr* + F2 -> KrF* + F, at
7.2e-10 with uncertainty 0.1, taken from Velazco, Kolts & Setser (Table 2,
p.46), and repeats "k14 = 7.2 x 10^-10 cm3/sec" in the text (p.67). "8.1"
appears in neither of two independent OCR texts of it, nor in the midterm
report MP-76-46. The first-hand values are 7.2 (Setser), 8.8 +- 0.5 and
9.4 +- 1.0 (Ueda 1984 Table V) e-10 cm3/s, mean 8.47e-10.

**Analogs (not the link).**
- Ar* + F2 -> ArF*: 0.5 "rather than unity" (Setser, COO-2807-11, p.6,
  checked) and <= 0.6 (SRI, secondhand via Brau in Rhodes (ed.), *Excimer
  Lasers* (1984), Table 4.4 p.98). SRI's own text (p.40): "The ArF
  fluorescence yields were surprisingly small (30% being the greatest
  observed)."
- Kr* + SF6 -> KrF*: G = 0.70 +- 0.10 per 100 eV (Cooper, Mezyk & Young,
  J. Phys. Chem. 100, 16195 (1996), abstract). With Saito's N_ex/N_i 0.55
  and W 24.2 eV that is about 0.31 KrF* per Kr* (our conversion).
- SRI's KrF* fluorescence yield "approaching 100% at low argon pressures"
  (MP-76-99, p.42) measures the whole chain in Ar-rich gas and does not
  isolate Kr* + F2.

**Could settle it (paywalled or unreachable):**
- J. E. Velazco, Ph.D. dissertation, Kansas State University (1977), the
  likeliest home of the recalibration. K-REx answered 403.
- Kolts, Velazco & Setser, J. Chem. Phys. 71, 1247 (1979),
  doi:10.1063/1.438480.
- Setser et al., Faraday Discuss. Chem. Soc. 67, 255 (1979),
  doi:10.1039/dc9796700255, and the General Discussion,
  doi:10.1039/dc9796700367.
- Kolts & Setser, in *Reactive Intermediates in the Gas Phase* (1979),
  pp. 151-231, doi:10.1016/b978-0-12-637450-6.50008-5.
- Chen & Payne, "Kinetics of proton-excited Ar-F2 and Kr-F2 mixtures", IEEE
  J. Quantum Electron. 15, 149 (1979), doi:10.1109/jqe.1979.1069984, the
  only high-pressure non-device Kr/F2 formation study found.
- Chen, Judish & Payne, J. Phys. B 11, 2189 (1978),
  doi:10.1088/0022-3700/11/12/017 (listed as open access, but IOP served
  the script an HTML page, quarantined unopened. It should open in a
  browser).

## 2. The KrF* yield and B/C split of Kr+ and Kr2+ + F-

No measurement with F- was found. Tsuji's group (228 works listed,
1985-2005) never used F- itself. Jacob, Rokni, Mangano & Brochu, Appl.
Phys. Lett. 32, 109 (1978) measured a rate from the KrF* fluorescence shape
and assumed the product. SRI in 1978 (MP-78-07, p.E-11, checked): "These
rates have not been measured for rare gas halides nor have the branching
ratios of reactions such as (8a) and (8b) been studied." Brau in 1984
(Rhodes p.99) calls the branching "believed to be very high, based on the
high efficiencies observed from rare gas halide lasers", a device inference
the predictive-only policy does not admit.

**Analogs (all abstracts):**

| Reaction | Value | Source |
|---|---|---|
| Kr+(2P3/2) + SF6- -> KrF B:C | 0.62:0.38 | Tsuji, Furusawa & Nishimura, J. Chem. Phys. 92, 6502 (1990) |
| Kr+(2P1/2) + SF6- -> KrF B:D | 0.05:0.95 | same |
| Kr+(2P3/2) + Cl- + He -> KrCl B:C | 0.60 +- 0.06 : 0.40 +- 0.04 | Tsuji et al., J. Chem. Phys. 94, 4291 (1991) |
| Kr+(2P3/2) + Cl- (+ He) -> KrCl(B,C) : Kr(5s 3P1) : Cl(4s) | 0.64 : 0.23 : 0.13 | Tsuji et al., J. Chem. Phys. 97, 1079 (1992) |
| Kr2+ + SF6- + Kr -> KrF* | G = 1.70 +- 0.10 per 100 eV, about 0.41 per ion pair (our conversion) | Cooper, Mezyk & Young 1996 |

**Could settle it:** Mezyk, Cooper & Sherwell, J. Phys. Chem. 95, 3152
(1991), doi:10.1021/j100161a037, the direct paper (ACS answered 403).
Companions: J. Phys. Chem. 93, 8187 (1989), 96, 8858 (1992), 97, 9413
(1993), and Radiat. Phys. Chem. 24, 545 (1984). Mezyk's thesis is not in
Melbourne's repository.

## 3. Kr2* + F -> KrF* + Kr

No measurement exists. Kannari's 3.0e-10 (Table I p.4311) cites ref. 52,
which is SRI MP-76-99. That report tabulates no such rate (checked, re-read
2026-09-28). It only interprets a slow component (p.42):

> "We also note that at all F2 pressures, there was a very slow component to
> the KrF* decay. At very low F2, most of the signal intensity is in the
> slow component. The interpretation is that this component results from
> the reaction Kr2* + F -> KrF* + Kr and that the decay of Kr2* is quite
> slow, as we see below."

Its measured 3e-10 is Kr2* + F2 -> Kr2F* + F (Table 2, reaction (21),
uncertainty 0.3, "From decay of Kr2* and Kr2F* vs F2"), now in the registry
next to Quigley & Hughes' 2.1 +- 0.2e-10. No measured KrF*/Kr2F* split of
Kr2* + F2 was found, and Chen & Payne 1979 is the likeliest place. SRI's
1978 report (MP-78-07, p.E-11): "the branching ratio between RgX and Rg2X
formation is not known."

**Re-read 2026-09-29, with the owner's papers.** The report's own model
leaves the reaction out (p.45): "An important class of reactions has not
yet been included; namely, reactions of Ar2* and Kr2* with F atoms produced
by previous dissociation of the F2." Its one number is for argon and
conditional (p.65): "if the reaction Ar2* + F -> ArF* + Ar is sufficiently
rapid (3-10 x 10^-10)" it could explain the Ar2F* yield falling at high
argon pressure. So 3.0e-10 is the low end of that range, carried to Kr. It
is no measurement and no upper bound. Chen & Payne, IEEE JQE 15, 149
(1979), Fig. 31, assign the slow KrF*(B) decay component in proton-excited
Kr/F2 to ion-ion recombination instead. Of the 76 PDFs searched (docs/,
docs/Recent_KrF_Papers/, docs/PIC_papers/, the sources here and the
owner's KrF_Reaction_Rate_Tables.pdf), every value is this 3.0e-10. The
owner's tables carry it (rows 121 and 124, ref. [16], the same SRI
report), and so do Ong 1996 (via Yukimura 1989) and Ma 2024 (via Kannari,
in m3/s). At LIGHT v0_1's
reference point the 248-nm share saturates in this rate: 4.65 % at zero,
5.85 % at 1e-11, 9.17 % at 3.0e-10 and 9.49 % at 1e-9.

## 4. What Morozov et al. 2008 counts

Morozov, Heindl, Krucken, Ulrich & Wieser, J. Appl. Phys. 103, 103301
(2008), doi:10.1063/1.2931000 (paywalled, and the mediaTUM copy is
restricted access).
- **Denominator:** "the energy deposited by the electron beam in the gas"
  (abstract). "Apparent" efficiencies for devices are reported separately.
- **Numerator:** the second continuum. Ulrich, Laser Part. Beams 30, 199
  (2012), arXiv:1511.08393 (checked): "Efficiencies for the conversion of
  electron beam power into VUV light of 31%, 32%, 42% and 42% have been
  measured for the so called second continuum radiation at 83, 126, 150,
  and 172 nm in Ne, Ar, Kr, and Xe, respectively (Morozov et al., 2008a)."
- **Scale:** a deuterium lamp traceable to PTB, with a coverage-factor-2
  uncertainty of 14 % at 122.8-170.0 nm (Heindl et al., JINST 6, P02011
  (2011), arXiv:1511.07720, checked).
- **Pressure:** 75-1400 hPa, with Kr at 0.42 +- 0.05 above 500 hPa.
- **Open:** whether the first-continuum shoulder, the atomic lines or the
  third continuum were integrated, and how backscatter was handled. Dandl
  et al., EPL 94, 53001 (2011), arXiv:1511.08395: "The first continuum is
  visible as a short wavelength 'shoulder' of the second continuum."

LIGHT's pure-Kr yield with no thermal excitation is 0.551, the ideal of one
photon per excitation and per ion (0.54 with the measured W and N_ex/N_i).
Morozov's 0.42 +- 0.05 is 2.4 sigma below the 0.54. The definition explains part
of that only if the second continuum misses a sizeable share of the light,
which the open sources do not settle.

## The four papers the owner supplied (2026-09-28, later)

The owner added the four papers this search could not reach:
- `sources/mezyk1991.pdf`
- `sources/Kinetics_of_proton-excited_Ar-F2and_Kr-F2mixtures.pdf` (Chen &
  Payne 1979)
- `sources/morozov2008.pdf`
- `sources/A_STUDY_OF_THE_TOTAL_AND_REACT.pdf` (Velazco's dissertation)

Each was checktexted CLEAN before it was opened. The text layers were read
first-hand, and the page images were used where the OCR garbled numbers.
The 318-page scanned dissertation was OCR'd by an agent. Its key passage,
pp.147-148, was re-read on the page images.

| Paper | What it measures | What it settles |
|---|---|---|
| Velazco 1977, Ph.D. dissertation (Kansas State) | Xe(3P2) quenching and excimer branching, relative to XeCl* | The Kr* + Cl2 reference is 0.90 absolute. The 1976 paper, supplied next, then gave the KrF* yield directly (below) |
| Chen & Payne 1979 | Proton-excited Ar/F2 and Kr/F2: rates and a kinetic scheme | A third Kr2* + F2 rate (5e-10). The direct Kr2* + F2 -> Kr2F* channel is minor |
| Mezyk, Cooper & Sherwell 1991 | Kr2+ + SF6- (+Kr) recombination rate, 30-1600 Torr | The ion-ion rate scale (about 3e-6 cm3/s near 665 Torr, LIGHT 2.9e-6). SF6-, not F-, and no yield |
| Morozov et al. 2008 | e-beam VUV efficiency of Ne, Ar, Kr, Xe | What the efficiency counts: all the VUV over the deposited power, so the 22 % shortfall is not a definition effect |

**Velazco 1977 (checked on the page images, pp.147-148).** Xe* + Cl2 ->
XeCl* "was assigned the value of 1.0 because: (i) comparison experiments
with Kr(3P2) for which the absolute metastable concentration could be
measured reliably, showed that the branching ratio for excimer formation,
k_KrCl/k^Kr_Cl2, was 0.90." The reference cited is Velazco, Kolts & Setser
1976, the paper whose abstract calls the Kr(3P2) + F2 branching "unity"
relative to Kr* + Cl2. p.147: "the only observed products from Kr(3P2) +
Cl2 are KrCl* and Cl2* 256 nm (trace)".
- With the 1976 abstract's relative "unity" this pointed to an absolute
  KrF* yield of 0.90. The 1976 paper measures it directly (next section).
- An earlier version of this file took the ratio 0.7 of the 1976 ERDA
  report (COO-2807-2 p.17) as a low end, 0.63. That ratio came before the
  group's spectral recalibration, and the published recalibrated ratio is
  0.88 (next section), so that low end no longer applies.
- No uncertainty is printed for the 0.90 here. The 1976 paper gives
  +-0.20.
- The dissertation has no Kr* + F2 rate or yield of its own. Its XeF*
  yield, 1.0 +- 0.20, is relative to XeCl* (Table 4.4). Its spectral
  simulation of KrF from Kr + F2 needs an "inverted ramp" vibrational
  distribution peaking near the highest allowed v' (p.193), so the nascent
  KrF is vibrationally hot, as Yu & Setser found.

**Chen & Payne 1979 (checked, pp.157-159 and Fig. 25).**
- Kr2*(1u) + F2 is 5 x 10^-10 cm3/s (beta2 = 1.6e7 Torr-1 s-1 +- 25 %),
  from the 1425 A continuum's decay at 600 Torr Kr. With Quigley & Hughes'
  2.1 +- 0.2 and SRI's 3 +- 0.9, the three measurements disagree beyond
  their errors.
- Fig. 25: f2'/f3 << 1, and p.158 says "the dominant process to produce
  Kr2F* is through KrF* as an intermediate species". The direct Kr2F*
  channel of Kr2* + F2 is minor, which contradicts the half LIGHT sends
  there. The KrF* share is not given, and the split moves formation by 1 %
  at most, so LIGHT keeps the even split with this note.
- Kr(3P1) + F2 is 6.8e-10 and Kr(1P1) + F2 is 7.6e-10 (p.157). They make
  KrF* "with efficiency comparable to Kr(3P2)" (p.159). At laser pressures
  collisions keep the 5s manifold mostly in 3P2.
- At their low ion density the ion-ion part of the KrF* light is "about 10
  percent" of the neutral part (p.159), which they put down to F- lost to
  F2 (2e-10 cm3/s per F2). At Tahoma recombination outruns that loss:
  0.7 % of the F- loss and about 0.3 % of formation.
- Kr2*(1u) lifetime "~1 us" (Fig. 23), against 250-353 ns from six other
  measurements.
- **Not used:** their KrF* + F2 (8.8e-10) and KrF* + 2Kr (7.3e-31). They
  are single-state Stern-Volmer products of time-integrated B emission,
  divided by a 7.5 ns rise-time lifetime. Their own trimer rate alone gives
  a decay rate of 2.7e8 1/s at 600 Torr (3.7 ns), so the pair is
  inconsistent, and the plots also fold in the F2 dependence of formation.

**Mezyk, Cooper & Sherwell 1991 (checked).** Kr/SF6 pulse radiolysis: "a
maximum value (2.9 x 10^15 M-1 s-1, ~4.8 x 10^-6 cm3 s-1) at approximately
300 Torr", then falling as the recombination becomes diffusion controlled.
The anion is SF6-, not F-. Their method cancels the KrF* yield per
recombination, so the paper gives no yield and no B/C split. Its Kr2+ + F-
curve (Fig. 4b) is a Bates calculation, not data.

**Morozov et al. 2008 (checked).** "Energy conversion efficiency is defined
here as the ratio of the power emitted by the target gas into 4 pi solid
angle in the VUV (or XUV) wavelength range and the power deposited in the
gas by the electron beam."
- For Kr above 500 hPa the apparent efficiency is 0.35 +- 0.02 (12 keV).
- Dividing by the deposited fraction 0.67 (GEANT4, +-10 %, backscatter
  included) and correcting for reflection on the foil (0.26) gives 0.42 +-
  0.05.
- The integration is not restricted to one continuum. The 22 % shortfall
  below the one-photon ideal is therefore physics or a systematic error,
  not a definition effect.

**What it does to LIGHT.** The basic network now uses the measured yield,
0.90 (row G3X takes the other 10 %). At Tahoma, a score only:

| KrF* yield of Kr* + F2 | g0 (6.39 %/cm) | Isat (1.4 MW/cm2) | g0 Isat / Pin (0.27) |
|---|---|---|---|
| 1.0 (the published assignment, LIGHT until 2026-09-28) | 6.39 (x1.00) | 1.54 (x1.10) | 0.298 (x1.10) |
| **0.90 (measured, adopted)** | **6.11 (x0.96)** | **1.53 (x1.09)** | **0.284 (x1.05)** |
| 0.80 | 5.83 (x0.91) | 1.53 (x1.09) | 0.270 (x1.00) |
| 0.70 (0.90 less the stated 0.20) | 5.54 (x0.87) | 1.52 (x1.09) | 0.256 (x0.95) |

The measured range brackets Tahoma's formation. Still unmeasured: the KrF*
yield and B/C split of ion-ion recombination with F-, and Kr2* + F.

## Velazco, Kolts & Setser 1976, supplied last (checked)

The owner then added the paper the dissertation cites,
`sources/velazco1976.pdf` (J. Chem. Phys. 65, 3468 (1976), checktext
CLEAN, text layer read first-hand). It settles the harpoon yield.
- **The relative rate.** Table III puts KrF* formation from Kr(3P2) + F2
  at 88 relative to KrCl* from Cl2 (100). "For the same gas mixture, the
  reproducibility of the relative rate constants was +- 10%", although
  "the comparison between KrCl* and KrF* could be subject to some
  systematic error" because the detection response "declines rapidly at
  short wavelength". This is the recalibrated ratio. The 0.7 of the ERDA
  report came before it.
- **The absolute scale.** Sec. III C ties the Cl2* 256 nm continuum to the
  Kr* 760.2 nm line from Ar(3P2) + Kr and measures [Ar*] and [Kr*] by
  absorption, giving k_KrCl = 74 x 10^-11 cm3/s. The comparison "uses the
  ends of our calibrated spectral response curve, which could be in error
  by +- 20%", and the [Kr(3P2)]/[Ar(3P2)] ratio "still could be +- 30%".
  "Within these experimental limitations, the branching ratios for
  formation of KrCl* and Cl2* are assigned as 0.90 and 0.10".
- **The KrF* yield.** "the measured value of the KrF* formation rate
  constant (65 x 10^-11 cm3 molecule-1 sec-1) equals the total quenching
  rate constant" within uncertainty (p.3476). The total is 72 x 10^-11
  (Table II). Their two F2 mixtures gave 61 and 72 (p.3473). The measured
  yield is therefore 65/72 = 0.90.
- **The published value.** Table IV lists Kr(3P2) + F2 as 1.0, and "The
  uncertainties are estimated as +-0.20". The branching ratios for XeF*
  and KrF* from F2 "were assigned as unity" for three reasons: the rate constants agree within uncertainty,
  the XeF*/XeCl* rate ratio matches the total quenching ratio, and "No
  other emissions were observed from the F2 reactions".
- LIGHT keeps the measured 0.90, inside the published 1.0 +- 0.20. The
  range 0.70 to 1.0 brackets Tahoma's formation.
- Also in the paper: "The secondary emission, which is thought to
  represent termination on the 2Pi, rather than the 2Sigma+ lower state,
  is very weak for XeF* and KrF*" at 1-2 Torr. That sits uneasily with Yu &
  Setser's 1990 nascent C share of 0.38. The 1990 study is the dedicated,
  later measurement, and LIGHT keeps it.

## The owner's 1977-78 papers (2026-09-29, checked)

Asked for after the Kr2* + F search came up empty, and supplied in
`Kinetic_Modeling/Recent_KrF_Papers/`, all checktext CLEAN. Jancaitis's
1983 LLNL report is still missing.

| Paper | What it is |
|---|---|
| Rokni, Jacob & Mangano, Phys. Rev. A 16, 2216 (1977) | e-beam ArF* and KrF* formation and quenching |
| Jacob, Rokni, Mangano & Brochu, Appl. Phys. Lett. 32, 109 (1978) | e-beam Kr/F2, ion-ion rate from pulse shapes |
| Rokni, Mangano, Jacob & Hsia, IEEE J. Quantum Electron. 14, 464 (1978) | Avco's review of rare-gas fluoride lasers |
| Quigley & Hughes, Appl. Phys. Lett. 32, 649 (1978) | photolytic Kr2* and Kr2F* lifetimes and quenching |

- **Kr2* + F.** None of the four gives it. Quigley & Hughes ran
  photolytically with 2e-4 to 7e-4 F2, so F atoms are far too few to
  bound it, and the Avco papers never mention the reaction.
- **Kr2* + F2.** Quigley & Hughes (p.650): "The relative probability for
  the occurrence of one route to Kr2F* or the other is not known and
  cannot be discerned in this experiment." Their 2.1 +- 0.2e-10 is the
  total, and a harpoon estimate from the 6.5 A crossing gives 2.9e-10.
- **Kr* + 2Kr.** Their Table I carries 5.5e-32 from Hughes 1977 (MIT
  Press), not read. The registry's spread is now 1.7 to 5.5e-32.
- **Ion-ion rate (section 2).** Jacob et al. (pp.110-111) fitted KrF*
  pulse shapes in e-beam Kr/F2 to Kr+ + F- (+M) = 2e-7 at 60 Torr and
  3e-7 cm3/s at 120 Torr, 1e-7 extrapolated to zero pressure, and
  estimated 1e-6 at 1 atm. v0_1's Flannery form, screened at their 1e14
  cm-3 ions, gives x1.2 and x1.6 of the two measurements.
- **Ion-ion KrF* yield (section 2).** Rokni et al. (p.2221) changed the
  beam current x25, which shifts Kr2+ against Kr+, and found that Kr2+ +
  F- forms KrF* "with a branching ratio near unity". The test used
  relative intensities only (the Ar case, Fig. 5, against branchings 1,
  0.8 and 0.6), so it rules out Kr2F* from Kr2+ + F- but not a loss
  common to both ions. The absolute yield and the B/C split stay
  unmeasured.
- **Kr2F*.** Quigley & Hughes measured 181 +- 12 ns and 4.3 +- 0.4e-10
  for Kr2F* + F2. Xu et al. 1993 weighed these against three other
  groups (1.5, 1.8 and 1.5e-10) and remeasured 1.2e-10 and 220 ns, the
  values LIGHT carries.

## Jancaitis 1983, supplied last (2026-09-29, checked)

K. S. Jancaitis, *Kinetic Modelling of Krypton Fluoride Laser Systems*,
Ph.D. thesis, LLNL UCRL-53465 (1983), in `Kinetic_Modeling/Recent_KrF_Papers/`
(checktext CLEAN). Printed page numbers, PDF page = printed + 8.

- **Kr2* + F.** His model's KrF* formation (Table II pp.94-95) is the ion
  pairs, Kr** + F2, Kr* + F2, Kr2* + F2 and Kr + F2*, with no Kr2* + F.
  So neither SRI's model nor the LLNL one carries it, and no measurement
  or bound is among the 81 sources now held. LIGHT v0_2 carries it as a
  band from 0 (its default) to 1e-9.
- **Kr2* + F2.** 3.0e-10, all to KrF* + Kr + F (p.95, his refs. 93 and 95,
  Brau and Quigley & Hughes), an assignment at SRI's total.
- **Kr* + 2Kr.** 2.6e-32 from Tracy & Oskam 1976 (p.86), inside the
  registry's 1.7 to 5.5e-32.
- **Penning.** Kr* + Kr* 1.5e-10 from Werner's MIT thesis (1975) and Kr2*
  + Kr2* 1.0e-10 from Werner et al., IEEE J. Quant. Elect. 13, 769 (p.86,
  refs. 21 and 1 on p.109). Kannari has 5.0e-10 for both. Neither original
  is held, so both stay bands in v0_2.
- **Kr2F*.** Measured 225 +- 21 ns and Kr2F* + F2 = 1.8 +- 0.2e-10 (Table
  II p.159). His model table (p.96) uses 1.6e-10. Xu 1993 weighed 1.5 to
  4.3e-10 and remeasured 1.2e-10, the value LIGHT carries.
- **Kr2F* + F.** With the F2 nearly all dissociated at high current, the
  trimer quench after the pulse matched the initial F2, "This suggests
  that the quenching rate coefficient of the RGH trimers by fluorine atoms
  is between .5 and 1 times the corresponding value for F2" (p.153). The
  analysis behind it is Fig. 11, Ar2F* in 2.5 Torr F2 and 1500 Torr Ar.
  For Kr2F* it is an analog, 0.9 to 1.8e-10 with his Kr2F* + F2, so it sets
  the top of v0_2's band for the row and never its default.

## Also noted (not asked)

COO-2807-18 (1979), p.15 (OCR): "For KrCl, KrF and XeI the high pressure
C/B ratios are 1.7, > 4.0 and » 1.0, respectively." For KrF this is a lower
limit, with all broadband emission assigned to C-A. It is older than the
Kansas State data in `BC_AND_ELECTRON_QUENCH.md`.

## Sources kept (`sources/`, checktext CLEAN, rechecked 2026-09-28)

- `SRI_MP-76-99_Nakano_1976_DTIC_ADA033997_text.pdf` and `_djvu.txt` (an
  independent OCR of the same report)
- `SRI_MP-76-46_McCusker_1976_DTIC_ADA027128.pdf` (the midterm report)
- `SRI_MP-78-07_Huestis_1978_DTIC_ADA056810_djvu.txt` (OCR text only, the
  PDF is 99.5 MB)
- `Setser_ERDA_COO-2807-{4,8,12,13,14,15,16,17,18}_*.pdf` (-12 to -18 are
  scans)
- `Ulrich_LaserPartBeams30_199_2012_arXiv1511.08393.pdf`,
  `Heindl_JINST6_P02011_2011_arXiv1511.07720.pdf`,
  `Dandl_EPL94_53001_2011_arXiv1511.08395.pdf`

## How the search ran

- 3 WebSearch calls. The rest went through Crossref, OpenAlex, Semantic
  Scholar, OpenAIRE and arXiv APIs, OSTI's API and full-text servlet, the
  archive.org DTIC mirrors, and Melbourne's repository API. No email
  address was sent in any request, and no paywall was worked around (ACS
  403, K-REx 403, mediaTUM restricted access behind a bot check).
- Every downloaded file was run through `~/text-sandbox/checktext`. Three
  HTML responses came back SUSPICIOUS (scripts) and were quarantined: an
  OSTI response, the IOP page, and the mediaTUM bot-check page.
- Two departures from the rule that checktext comes before any reading.
  The agent parsed the mediaTUM bot-check page as text before scanning it.
  It also read 73 API responses (JSON and XML) before checktexting them,
  all CLEAN afterwards. Nothing from those pages is used here.
- Quotes were byte-checked by the agent against the raw text. The SRI
  report's reactions (14) and (21) and its p.42 passage were re-read
  first-hand for this file.

## URLs fetched

Final HTTP status first. 000 means no connection.

Crossref: `200` api.crossref.org/works?query.bibliographic= for Kolts
Velazco Setser JCP 71 1247, Setser Faraday Discuss 67 255, Mezyk Cooper
Sherwell JPC 95 3152, Jacob Rokni Mangano Brochu APL 32 109, Morozov et al.
conversion efficiencies, Chen Payne IEEE JQE 15 149, Quigley Hughes APL 32
649, Johnson Hunter JAP 51 2406, Cooper Mezyk Armstrong Radiat. Phys. Chem.
1984, and "Cooper rare gas fluoride exciplex pulse radiolysis" (1980-95).
`200` api.crossref.org/works/ for 10.1021/jp961251q, 10.1021/jp960129o,
10.1021/j100139a027, 10.1063/1.2931000, 10.1063/1.89955, 10.1063/1.458284,
10.1063/1.460615, 10.1063/1.463287, 10.1109/jqe.1979.1069984,
10.1063/1.89882, 10.1063/1.328010, 10.1063/1.456349,
10.1016/0146-5724(84)90191-2.

OpenAlex (`200` unless marked): works?filter=doi: (four batches), authors?
search=Andrei Morozov (TU Munich), works?filter=cites:W2030971770, searches
for Velazco metastable xenon, ion-ion recombination KrF flowing afterglow,
Kr+ F- recombination KrF formation, branching fraction rare gas halide
formation, proton excited Kr F2 mixtures, Mezyk Cooper Sherwell, pulse
radiolysis rare gas fluoride excimer ion recombination, dissertations on
rare gas halide ion recombination pulse radiolysis, author works for
A5048445296, A5078517768, A5109969001, A5025212885, A5061869813,
A5035412340 and A5053559574, works/doi:10.15017/6583 and
works/doi:10.1016/0146-5724(84)90191-2. `400` for one search with a
literal "Kr2*".

Semantic Scholar: `200` graph/v1/paper/DOI: for 10.1021/j100161a037,
10.1021/j100362a010, 10.1021/j100201a033, 10.1063/1.2931000,
10.1063/1.438480, 10.1039/dc9796700255, 10.1021/jp961251q.

OpenAIRE and arXiv: `200` api.openaire.eu (doi 10.1063/1.2931000),
export.arxiv.org queries (au:Ulrich_A AND au:Wieser, ti:"conversion
efficiencies" AND ti:vacuum), arxiv.org/pdf/1511.08393v1, 1511.08395v1,
1511.07720v1.

OSTI: `200` api/v1/records for q=COO-2807, author=Setser (1975-85),
Lorents and Huestis (1974-82), Nakano with title=krypton, Payne (1975-82),
Chen with title=proton (1975-82), title=KrF author=Chen (1976-81). `000`
for records 6474864 and 5066237, `404` for record 5383187. `200` for servlets/purl/
7326141, 7289512, 5070777, 6935020 (first try returned HTML, quarantined),
6598867, 5675592, 6012435, 5738041, 5888225. `000` for servlets/purl/4089579
(COO-2807-1, never retrieved).

archive.org: `200` metadata for DTIC_ADA033997, DTIC_ADA049829,
DTIC_ADA027128, DTIC_ADA056810, and downloads DTIC_ADA033997_text.pdf,
DTIC_ADA033997_djvu.txt, DTIC_ADA049829.pdf, DTIC_ADA027128.pdf,
DTIC_ADA056810_djvu.txt.

Repositories and publishers: `200` mediatum.ub.tum.de/node?id=796917 (the
bot-check page) and its OAI-PMH record. `403` for the krex.k-state.edu
search API. `200` for the minerva-access.unimelb.edu.au search API (Mezyk,
exciplex rare gas, pulse radiolysis gas, rare gas halide, ion recombination
exciplex), `500` for its rest/items/find-by-metadata-field. `200` for
catalog.lib.kyushu-u.ac.jp/opac_download_md/6583/KJ00004507583.pdf (Tsuji
& Nishimura 1991, not relevant) and for
iopscience.iop.org/article/10.1088/0022-3700/11/12/017/pdf (served HTML,
quarantined).

WebFetch: mediatum.ub.tum.de/node?id=796917 (bot check), krex.k-state.edu
search (403), pubs.acs.org/doi/abs/10.1021/j100161a037 (403),
osti.gov/etdeweb/biblio/5383187 (200), osti.gov/biblio/5066237 (200).
