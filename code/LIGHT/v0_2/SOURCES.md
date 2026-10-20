# v0_2 sources, from v0_1's and for the Kr2* route and the bands

Sections 1 to 4 are v0_1's, the literature behind its changes of
2026-09-28, which v0_2 keeps. Section 2's ion-ion rows were corrected on
2026-09-29 from the owner's 1977-78 papers. Section 5 is v0_2's, the
Kr2* route and the rates carried as bands (2026-09-29). Four research agents
searched the local library (`docs/Recent_KrF_Papers/`, the Rhodes book
`docs/Excimer_Lasers_Topics_Rhodes.pdf`, `LIGHT/docs/krf_data/sources/`)
and free APIs (Crossref, OpenAlex, OSTI, arXiv, NIST, DTIC via
archive.org). Every download was scanned with `~/text-sandbox/checktext`
before it was read. The numbers that changed the model were re-read here
first-hand in the local sources, as the status column says.

The status column says how each value was read. **checked** means read in
the full text (page given), **abstract** means the abstract only,
**secondhand** means quoted by another source that was read (named),
**derived** means arithmetic on the cited numbers, and **estimate** means a
researcher's own estimate, never used to set a value.

## 1. The trimer, KrF* + 2Kr -> Kr2F* + Kr

| Quantity | Value | Source | Status |
|---|---|---|---|
| rate at 300 K, B/C mixture | (6.0 +- 2.0)e-31 cm6/s | Xu, Gadomski & Setser, JCP 99, 2591 (1993) Table II p.2596 | checked |
| B and C separately | not measured. Gadomski 1992 p.155 assumed them equal | Gadomski et al., CPL 189, 153 (1992) | checked |
| temperature law, theory | Shui's classical phase-space calculation, computed 300-800 K. 300 -> 600 K lowers k by 10 | V. H. Shui, APL 34, 203 (1979), and the Avco report DTIC ADA069053 Sec. V.B pp.68-73, Fig. 22 | abstract, and checked on page images by the agent |
| temperature law, measured | 315-420 K, agreeing with Shui to 15 % (Ar-rich mixtures) | Klimek & Hsia, JAP 52, 5361 (1981), and Duzy et al., Avco report OSTI 5268869 Sec. IV.F pp.52-60 (in `docs/krf_data/sources/`) | abstract, and checked |
| temperature law, pure Kr/F2 | "an inverse square root temperature dependence for the Kr2F* formation, which is weaker than earlier predictions", 294-425 K | Kimura & Salesky, APL 48, 897 (1986), and Salesky & Kimura, APL 47, 774 (1985) | abstract |
| data above 425 K | none found | | |
| Kr2F* binding below KrF* + Kr | 0.58 +- 0.08 eV (calculated) | Huestis, Marowsky & Tittel, in *Excimer Lasers*, 2nd ed., C. K. Rhodes ed. (Springer 1984), ch. 6, Table 6.3 p.185 | checked |
| same, other calculations | 0.49 eV with spin-orbit (Wadt & Hay, JCP 68, 3850 (1978)), and 0.66 eV DIM (SRI MP-78-07 App. G) | | abstract, and checked |
| Kr2F* emission band maximum | 420 +- 35 nm | Huestis et al. 1984, Table 6.1 p.185 | checked |
| Kr2F* radiative lifetime | 220 +- 20 ns | Xu 1993 p.2597 | checked |
| Kr2F* thermal redissociation | no published rate. An estimate puts it near the radiative rate from about 600 K | the trimer agent's detailed-balance sketch | estimate |

**Used in v0_1.** The rate and the (Tg/300)^-3.32 law stay (the joint fit's
300 K values, Shui's law), and the law's stated range becomes 300-800 K.
Kr2F* sits 0.58 eV below KrF(B) + Kr (it was REBUILD's estimated 2.05 eV).
C3 emits at 420 nm. The pure-Kr T^-1/2 was inferred by fitting a laser
model to laser efficiencies, so it bounds the law and never sets it
(`Params.trimer_tg_pow` runs the band).

## 2. Charged-particle recombination

| Quantity | Value | Source | Status |
|---|---|---|---|
| Kr2+ + e -> Kr* + Kr | 1.6e-6 (300/Te)^0.55 cm3/s +- 10 %, Te 300-8400 K, Tgas 300 K | Shiu & Biondi, PRA 16, 1817 (1977) | abstract (quote checked) |
| same | 1.48e-6 (Te/300 K)^-0.53, Te 300-21300 K, 3-12 Torr | Mikus, J. Phys. D 11, L39 (1978) | abstract (quote checked) |
| same | (5.5 +- 3.5)e-7 at Te 1200 K | Howorka & Kuen, JCP 70, 758 (1979) | abstract |
| DR products | mostly 5p states, and "excited states of Kr* having energies up to that of the Kr2+ ion" | Shiu & Biondi 1977 | abstract |
| DR gas-temperature dependence | not measured | | |
| Kr+ + 2Kr -> Kr2+ + Kr | 2.3e-31 cm6/s at 300 K | Bhattacharya, PRA 1, 1196 (1970) | abstract (quote checked) |
| same | (2.2 +- 0.2)e-31 | Grossl et al., JCP 74, 1728 (1981) | abstract (quote checked) |
| same | 342 +- 15 Torr^-2 s^-1 at 300 K (2.74e-31 with Torr at 273 K, 3.30e-31 at 300 K). "a marked decrease in beta as the gas temperature is increased to 510 K" | Smith, Dean & Plumb, J. Phys. B 5, 2134 (1972) | abstract (quote checked) |
| its temperature law, theory | "k3 prop. to T^-3/4" for rare-gas ion-atom-atom recombination | Perez-Rios & Greene, PRA 98, 062707 (2018), arXiv 1810.02108 (in `docs/krf_data/sources/`) | checked |
| Kr+ + e + e | 2.0e-27 n_e Te(eV)^-4.5 cm3/s, below 1e-4 of DR at LIGHT's conditions | Mansbach & Keck, PR 181, 275 (1969) | abstract |
| ion-ion with F- | Brau found no measurement for any rare-gas cation | Brau in *Excimer Lasers* 1984 p.98, and SRI MP-78-07 p.E-11 | checked |
| Kr+ + F- (+M), below 300 Torr | 2e-7 cm3/s at 60 Torr and 3e-7 at 120 Torr, fitted to KrF* pulse shapes in e-beam Kr/F2. LIGHT's form, screened at their 1e14 cm-3 ions, gives x1.2 and x1.6 of these | Jacob, Rokni, Mangano & Brochu, APL 32, 109 (1978) | checked |
| Kr2+ + F- against Kr+ + F- | forms KrF* "with a branching ratio near unity", from relative intensities as the beam current changed x25 | Rokni, Jacob & Mangano, PRA 16, 2216 (1977) | checked |
| ion-ion, analog | Kr2+ + SF6- rises to 4.8e-6 cm3/s near 300 Torr, then falls | Mezyk, Cooper & Sherwell, JPC 95, 3152 (1991) p.3155 (in `docs/krf_data/sources/`) | checked |
| ion-ion, theory | peak (2-3)e-6 cm3/s above 1 atm (Wadehra & Bardsley, APL 32, 76 (1978)), and about 3e-6 (Flannery & Yang, APL 33, 574 (1978)) | | abstract |
| Kannari's ion-ion form | "the gas-temperature-dependent recombination rate, the modified Natanson calculation [Hoffman & Moreno, SAND80-1486] is directly introduced and their profile is used" | Kannari et al., JAP 57, 4309 (1985) Sec. II.B p.4313 | checked |
| Thomson three-body part | proportional to N T^-3 at fixed density (critical radius b = 2e^2/3kT, polarization mean free paths), and the mobility part proportional to 1/N | Brau 1984 eqs. 4.11-4.12, and the same (300/T)^3 in Mezyk 1991's Bates scaling p.3156 | checked for Mezyk, secondhand for Brau (read by the agent) |
| F- + F -> F2 + e | endothermic by 1.80 eV, since EA(F) 3.4011895 eV (Blondel et al. 2001) exceeds D0(F2) | NIST WebBook and JANAF | checked (by the agent) |
| e + F- -> F + 2e | measured only 13-1000 eV (peak 2.74e-16 cm2 at 68 eV), negligible for Te <= 3 eV | Peart, Forrest & Dolder, J. Phys. B 12, L115 (1979) | abstract |
| F- + F2 | 2e-10 cm3/s effective two-body, products unknown | Chen & Payne, IEEE JQE 15, 149 (1979) pp.156, 159 (in `docs/krf_data/sources/`) | checked (by the agent) |

**Used in v0_1.** Charge-conserving n_e with no injection. Kr2+ + e at
Shiu & Biondi (R11p, rate kind te_pow, where a captured electron takes
0.95 Te).
Kr+ + 2Kr at 2.4e-31 cm6/s, the mean of the three measurements, on
Perez-Rios & Greene's T^-3/4 (K15). F2 attachment, dissociation and
vibration from the SIGLO set (ATT, DIS, VIB, as v0's f2e group). Ion-ion
recombination keeps Kannari's density profile at 300 K and gets back its
gas-temperature dependence through (Tg/300)^-3 on the three-body part
(K26, K108). Near 1 atm its rate, its absolute KrF* yield and its B/C
split have no measurement. The one measurement is below 300 Torr. Not included: F- + F2 (products unknown), Kr2+ and Kr2F*
redissociation above about 1000 K (no rates).

## 3. Fluorine-atom recombination, F + F + M -> F2 + M

Rates follow the convention d[F2]/dt = k[F]^2[M] (Baulch et al., JPCRD 10
Suppl. 1 (1981) p.1-11).

| Quantity | Value | Source | Status |
|---|---|---|---|
| M = Ar, 295 K | (8.0 +- 0.6)e-35 cm6/s | Ganguli & Kaufman, CPL 25, 221 (1974), via Baulch 1981 p.1-22 | secondhand |
| M = He, 298 K | (6.0 +- 1.1)e-34 | Ultee, CPL 46, 366 (1977), via Baulch 1981 p.1-22 | secondhand |
| M = Kr | not measured at any temperature | | |
| evaluation, M = Ar and every M but F, F2 and He | "10^13.7 exp(-35,100/RT) is recommended", which implies 10^14.0 cm6 mol^-2 s^-1 = 2.8e-34 cm6/s at 300 K, "uncertain by at least an order of magnitude" | Cohen & Bott, Aerospace TR-0076(6603)-2 = SAMSO-TR-76-82 (1976), DTIC ADA024332, pp.12-16 (the F2 pages are in `docs/krf_data/sources/`), and the 1982 update, DTIC ADA122175 App. A | checked |
| KrF models | 1.0e-33, "Estimated" | Johnson & Hunter, JAP 51, 2406 (1980) Table I p.2408, and Kannari 1985 Table I p.4311 | checked |
| D0(F2) | 1.602 eV | NIST-JANAF table F-001 | checked (by the agent) |
| KrF* + F quenching | no measurement. Kannari's 7.8e-10 is "analogous to F2" | Kannari 1985 p.4311 | checked |
| KrF(X) well | 13 meV (0.31 kcal/mol +- 10 %), not a sink for F | Becker, LBL-8589 (1978) Table I p.36 | checked (by the agent) |

**Used in v0_1.** FF: F + F + Kr -> F2 + Kr at 2.8e-34 cm6/s, held at its
300 K value (the evaluation falls with Tg, so this overstates the channel).
At 1 atm with F near 1e17 cm^-3 the recombination time is milliseconds. At
the reference point it re-forms 0.5 % of the F2 fill within the run and the
F2 still burns out.

## 4. The KrF(C) state, the C-A band and its decay

| Quantity | Value | Source | Status |
|---|---|---|---|
| C3/2 -> A3/2 band maximum | 275 nm | Huestis et al. 1984, Table 6.1 p.185 | checked |
| same | "The observed emissions at 220 nm and 275 nm are assigned to the IV1/2-I1/2 and II3/2-I3/2 transitions" | Hay & Dunning, JCP 66, 1306 (1977) | abstract |
| C-A band width | "ca. 25 nm for 300 K conditions", "~22 nm at 230 K" | Yu & Setser, JPC 94, 2934 (1990) p.2937 | checked (by the agent) |
| B1/2 -> A1/2 | (272) nm, "Theoretical values in parenthesis" | Huestis et al. 1984, Table 6.1 p.185 | checked |
| A(B-A)/A(B-X) | "calculated to be 0.035" | Dunning & Hay 1978, quoted by Yu & Setser 1990 p.2941 | checked |
| C -> X | "The B-X transition dominates I1 with negligible C-X transition" | Yu & Setser 1990 p.2941 | checked (by the agent) |
| tau_C, measured | 60 +- 10 ns, the only measurement | Hakuta, Fujino, Nakayama & Takuma, JJAP 29, 1194 (1990) Table I | checked |
| tau_C, calculated | 75.2 ns | Dunning & Hay 1978, via Yu & Setser p.2935 and Kannari p.4314 | secondhand |
| B-C gap | "the consensus value from our work seems to be -120 +- 50 cm-1" (C below B). -110 +- 20 from the absolute ratio, -180 +- 30 from the 230 and 300 K ratios | Yu & Setser 1990 p.2941 | checked |
| B-C gap, v0's 80 cm-1 | Julienne & Krauss's reanalysis of Kolts & Setser's low-pressure intensities | Kannari 1985 p.4313 | checked (by the agent) |
| KrF(B) well | D_e(B) = 5.54 +- 0.03 eV below Kr+ + F-, R_e 2.27 A | Tellinghuisen, Hays, Hoffman & Tisone, JCP 65, 4473 (1976), quoted by Krauss & Mies in *Excimer Lasers* 1984, Table 2.2 p.14 | checked |
| Hakuta's conditions | 248-nm pump at 10 MW/cm2, 10 ns time resolution not deconvolved, C-A detected at 280 nm, 1-70 Torr Kr with 1 Torr F2. The fit held tau_B 6.8 ns and equal B and C quench | Hakuta 1990 pp.1194-1196 | checked (by the agent) |

**Used in v0_1.** G4C at the printed 275 nm (4.51 eV). A B -> A branch
(G4A, 272 nm, 3.4 % of B's radiative decay) under the C-A band, sharing
tau_B with B -> X. C stays at the joint fit's 50 ns, and Hakuta's measured
60 +- 10 ns lies inside the fit's 41-63 ns. The B-C gap is Yu & Setser's
120 cm^-1. KrF(B) sits D_e(B) below the ion pair, 5.856 eV above Kr + F2
(REBUILD's estimate was 6.15 eV), so a 248-nm photon leaves 0.06 eV of
X-state repulsion and a C-A photon 0.53 eV, as the spectroscopy implies.

**Hakuta's conflict stays open.** Candidate causes, none established: B ->
A light in the 280-nm channel, the undeconvolved 10-ns response, and cold
versus hot B/C preparation. No later paper reconciles the absolute rates.

## 5. The Kr2* route and the bands (v0_2)

The owner's step 2 asked for Kr2* + F to be treated like the trimer law,
as a band that only a measurement can set, and step 3 for the measured
spreads to be carried as bands. The sources below are in
`Kinetic_Modeling/Recent_KrF_Papers/` and `LIGHT/docs/krf_data/sources/`.
Jancaitis is his LLNL thesis, *Kinetic Modelling of Krypton Fluoride Laser
Systems*, UCRL-53465 (1983), with printed page numbers.

| Quantity | Value | Source | Status |
|---|---|---|---|
| Kr* + 2Kr -> Kr2* + Kr, 0.5-1 atm | 1.7 +- 0.2e-32 cm6/s at 300 K, e-beam, Kr* absorption | Ueda, Kanada, Sasaki & Takuma, Rev. Laser Eng. 12, 357 (1984) Table V | checked |
| same, low pressure | 1.79e-32 (Salamero 1979), 2.59e-32 (Tracy & Oskam, JCP 65, 1666 (1976)), 4.2e-32 (Smith & Turner 1963), 4.4e-32 (Bouciqué & Mortier 1970) | abstracts | abstract |
| same | 5.5e-32 | Hughes 1977, quoted by Quigley & Hughes, APL 32, 649 (1978) Table I | secondhand |
| same, Jancaitis's model | 2.6e-32, citing Tracy & Oskam | Jancaitis Table II p.86 | checked |
| its temperature law | Ea = -0.020 eV, which is (Tg/300)^-0.64 over 300-451 K | Wieme & Lenaerts, JCP 74, 483 (1981) | abstract |
| Kr2* + F2, both channels | 2.1 +- 0.2e-10 cm3/s. Their products "cannot be discerned in this experiment" | Quigley & Hughes 1978 p.650 | checked |
| same | 3.0e-10 | SRI MP-76-99 Table 2, reaction (21) | checked |
| same | 5.0e-10 | Chen & Payne, IEEE JQE 15, 149 (1979) p.157 | checked |
| same, Jancaitis's model | 3.0e-10, all to KrF* + Kr + F | Jancaitis Table II p.95 | checked |
| Kr2* + F -> KrF* + Kr | no measurement. SRI's model leaves Rg2* + F out (p.45), and its one number is a conditional "3-10 x 10^-10" for Ar2* + F (p.65). Kannari carries 3.0e-10 | SRI MP-76-99, Kannari 1985 Table I p.4311 | checked |
| same, Jancaitis's model | no such row among the KrF* formation reactions | Jancaitis Table II pp.94-95 | checked |
| Kr2* + Kr2* -> Kr2+ + 2Kr + e | no measurement. Kannari's 5.0e-10 is Lorents & Olson's unpublished 1973 SRI report, which gives Ar2* the same value | Kannari 1985 Table I p.4312, ref. 87 | checked |
| same, Jancaitis's model | 1.0e-10, from Werner et al., IEEE J. Quant. Elect. 13, 769 | Jancaitis Table II p.86 and ref. 1 p.109 | secondhand |
| Kr* + Kr* -> Kr+ + Kr + e | 5.0e-10, from Leichner & Ericson, PRA 9, 251 (1974) | Kannari 1985 Table I p.4312, ref. 86 | secondhand |
| same, Jancaitis's model | 1.5e-10, from Werner's MIT thesis (1975) | Jancaitis Table II p.86 and ref. 21 p.109 | secondhand |
| Kr2F* + F2 | 1.2 +- 0.2e-10 cm3/s, and 220 +- 20 ns. Earlier, 4.3 +- 0.4e-10 (Quigley & Hughes), 1.5e-10 (Morgan & Szoke), 1.8e-10 (Ueda) and 1.5 +- 0.2e-10 (Setser's group). "Our new result supports the smaller range of values" | Xu, Gadomski & Setser, JCP 99, 2591 (1993) Table I p.2594 and p.2595 | checked |
| same | 1.8 +- 0.2e-10 cm3/s, and 225 +- 21 ns | Jancaitis Table II p.159 | checked |
| Kr2F* + F | no Kr measurement. In Ar/F2 at high current, with the F2 nearly all dissociated, Ar2F* was quenched after the pulse as the initial F2 would quench it, so F atoms quench "between .5 and 1 times the corresponding value for F2" | Jancaitis p.153 and Fig. 11 (Ar2F*, 2.5 Torr F2 in 1500 Torr Ar) | checked |
| Kr2+ binding | D0 = 1.15 eV, by photoionization | Ng, Trevor, Mahan & Lee, JCP 66, 446 (1977), via Krauss & Mies in *Excimer Lasers* 1984, Table 2.1 p.9 | secondhand |
| Kr2* binding | none among the sources held. Krauss & Mies' Table 2.1 gives the Rydberg dimer's D_e for He, Ne and Ar only | *Excimer Lasers* 1984, Table 2.1 p.9 | checked |

**Used in v0_2.** The Kr2* route is part of the basic network. Kr* + 2Kr
keeps Ueda's 1.7e-32 on Wieme & Lenaerts' law, Kr2* radiates at 270 ns,
and Kr2* + F2 keeps Quigley & Hughes' 2.1e-10 split evenly. Kr2* + F,
Kr2* + Kr2* and the new Kr2F* + F row have no measurement, so they are 0
unless a band run sets them. Every rate in this section with a published
spread reads a band (`network.BANDS`), Kr* + Kr* and Kr2F* + F2 among
them, and `run_v0_2.py --bands` runs each end. The Ar2F* inference sets
the top of Kr2F* + F's band and never its default, since an analog never
sets a link. Kr2+ keeps 1.15 eV, now cited to the measurement, and Kr2*
keeps REBUILD's estimated well, which only moves energy between heat and
storage.
