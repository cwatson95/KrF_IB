# LIGHT_ArF: the numbers the network uses, checked first-hand (2026-09-30)

Every number LIGHT_ArF v0_2's network and species table take as a default
or a band end, and the page it was read from. "First-hand" means the page
image (or, for a text-layer PDF, the page text) was read in this session,
not a secondary quotation. Paths are relative to `Kinetic_Modeling/`.
Status: **M** measured, **M-via** a measurement read in a source that
quotes it, **C** calculated, **D** derived here by arithmetic.

Abbreviations for the PDFs:
- *Rhodes* = `docs/Excimer_Lasers_Topics_Rhodes.pdf` (Excimer Lasers, 2nd
  ed., C. K. Rhodes ed., Springer 1984; printed page + 10 = PDF page)
- *Rokni 1977* = `Recent_KrF_Papers/Dominant formation and quenching
  processes in E-beam pumped ArF and KrF lasers.pdf` (Rokni, Jacob &
  Mangano, PRA 16, 2216 (1977))
- *Chen & Payne* = `LIGHT/docs/krf_data/sources/Kinetics_of_proton-excited_Ar-F2and_Kr-F2mixtures.pdf`
  (IEEE JQE 15, 149 (1979))
- *Jancaitis* = `Recent_KrF_Papers/janacaitis.pdf` (UCRL-53465, 1983)
- *Ueda 1983 / 1984* = `LIGHT/docs/krf_data/sources/Ueda_RevLaserEng11_576_1983.pdf`,
  `Ueda_RevLaserEng12_357_1984.pdf`
- *Setser 1977* = `LIGHT/docs/krf_data/sources/Setser_ERDA_COO-2807-11_1977_OSTI6681231.pdf`
  (scanned; read as page images)
- *Duzy 1980* = `LIGHT/docs/krf_data/sources/Duzy_AVCO_RGHlasers_1980_OSTI5268869.pdf`
- *Velazco 1976* = `LIGHT/docs/krf_data/sources/velazco1976.pdf` (JCP 65,
  3468 (1976))
- *Lee 1989* = `docs/Recent_KrF_Papers/lee1989.pdf` (IEEE Trans. Electron
  Devices 36, 2053 (1989))
- *Petrova 2022* = `docs/Recent_KrF_Papers/petrova2022_orestes.pdf` (NRL
  memo, the Orestes ArF model)
- *Kannari 1985* = `Recent_KrF_Papers/kannari1985.pdf` (JAP 57, 4309 (1985))
- *SRI MP-76-99* = `LIGHT/docs/krf_data/sources/SRI_MP-76-99_Nakano_1976_DTIC_ADA033997_djvu.txt`

## Energetics (species.py)

| Quantity | Value | Source, page | Status |
|---|---|---|---|
| Ar* (4s manifold, at 3P2) | 11.548 eV | IST-Lisbon file threshold (vendored, md5 pinned) | M |
| Ar+ | 15.76 eV | IST-Lisbon ionization threshold | M |
| D_e(Ar2+) | 1.23 eV | Rhodes Table 2.1 p.9 (Krauss & Mies), photoionization of Ng et al. 1977 [2.7] | M-via |
| D_e(Ar2*) | 0.68 eV | Rhodes Table 2.1 p.9, Saxon & Liu 1976 [2.15] | C |
| ArF B-X band peak | 1933 Å = 6.414 eV | Chen & Payne p.150 ("The peaks at 1933 Å and 2840 Å correspond to the emission from ArF* and Ar2F*") | M |
| E_form(ArF(B)) | 0.80 + 6.414 = 7.214 eV | the spectroscopic lower bound (flat X state); D_e(B) blank in Rhodes Table 2.2 p.14 | D |
| ArF B/C splitting | about 0.25 eV (so C is not carried) | Rhodes p.15 (Krauss & Mies, calculation of Michels [2.39]) | C |
| D(Ar2F* below ArF* + Ar) | 0.65 ± 0.08 eV | Rhodes Table 6.3 p.185 (Huestis, Marowsky & Tittel) | C |
| Ar2F* band peak | 2840 Å = 4.366 eV | Chen & Payne p.150 | M |
| Ar2*(1u) band peak | 1270 Å = 9.763 eV | Chen & Payne p.150 | M |

## Rates (network.py)

| Row / band | Value used | Source, page, what it says | Status |
|---|---|---|---|
| ar_f2 (Ar* + F2 total) | 8.1e-10 (mean) | Velazco 1976 Table II: Ar(3P2) + F2 kQ = 75e-11 cm3/s, σ 132 Å² | M |
| | | Chen & Payne Table II p.160: 9.4e-10 (their ref. 47), Velazco's 7.5 beside it | M-via |
| | | Ueda 1984 Table V p.365: Ar*(3P2) absorption, Ar* + F2 (7.4 ± 1.0)e-10; ArF* fluorescence (6.4 ± 0.3)e-10 (band bottom) | M |
| | band top 1.1e-9 | Ueda 1983 p.584 Fig. 12 and Table IV: (1.1 ± 0.1)e-9 at 2.5 atm | M |
| arf_yield | 0.5 | Setser 1977 pp.-5-/-6-: "Total quenching rate constants and the branching fractions for formation of ArF* ... were measured", "Molecular F2 is still the best donor, but the branching ratio is ∿ 0.5 rather than unity" | M |
| | | Ueda 1984 Table III p.365: 4.0e-10 "J. H. Kolts et al." | M-via |
| | f1/f2 = 0.55 ± 0.25 | Chen & Payne Table I p.155 | M (relative) |
| G4, tau_B | 4.2 ns | Rhodes Table 4.2 p.90 (Brau): τ(Σ→Σ) 4.2 ns [Hay & Dunning, calculation] | C |
| QF2 (k τ) | 7.55e-18 cm3 (mean) | Rokni 1977 Table II p.2218: 7.6 ± 0.7e-18 cm³ (k = 1.9e-9 at their assumed 4 ns) | M |
| | | Chen & Payne p.155: "Values of τ3β3 and τ3β3′ are obtained as 7.5 × 10⁻¹⁸ cm³ and 2.1 × 10⁻³⁹ cm⁶" | M |
| QAR (k τ) | 3.6e-20 cm3 | Rokni 1977 Table II p.2218: 3.6 ± 1e-20 cm³ | M |
| arf_trimer (k τ) | 2.27e-39 cm6 (mean) | Rokni 1977 Table II: 1.6 ± 0.3e-39 cm⁶; Chen & Payne p.155: 2.1e-39; Jancaitis p.150: "kτd products of 3.1×10⁻³⁹ cm⁶ for Ar2F* formation from ArF* in argon", ±25 % | M |
| | band top 5.0e-39 | Rhodes Table 6.4 p.191: (1.2 ± 0.2)e-30 cm⁶/s [Marowsky 1982] × 4.2 ns | M-via, D |
| trimer_tg_pow | 0 (bottom −3.9) | Petrova 2022 p.25 (PDF 32): "ArF* + 2Ar → Ar2F* + Ar, 1.62 × 10⁻³² + 2.245 × 10⁻²¹ Tg^−3.9, k28, Shui-1979, Chen-1978" | model form after C |
| arf_e (k_e τ) | 5.75e-16 cm3 (mean) | Duzy 1980 Table III p.-23-: ArF* + e 1.6e-7 (44.4 A/cm², 0.25 atm) and 3.0e-7 (10.1 A/cm²); p.-21-: the slope is k_el τ_rad/k_ATT, "values of τ_rad = 6.5 ns and k_ATT = 4.5 × 10⁻⁹ cm³/s were used" (KrF), "The ArF* results ... short radiative lifetime, τ = 2.5 ns" | M |
| C3 (Ar2F* radiative) | 203 ns (mean) | Chen & Payne Table I: 185 ± 25 % ns; Jancaitis Table I p.158: 205 ± 20 ns; Ueda 1984 Table V: 220 ns | M |
| Q2F (Ar2F* + F2) | 1.93e-10 (mean) | Chen & Payne Table I: 2.1e-10 ± 25 %; Jancaitis Table I p.158: (1.79 ± 0.06)e-10; Ueda 1984 Table V: (1.9 ± 0.1)e-10 | M |
| ar2f_f (Ar2F* + F) | 0.75 × Q2F, band 0.5–1 × | Jancaitis p.153: "the quenching rate coefficient of the RGH trimers by fluorine atoms is between .5 and 1 times the corresponding value for F2" | M (inferred) |
| ar2f_ar (Ar2F* + Ar) | 0, top 2.2e-14 | Jancaitis Table I p.158: kAr < 1e-14 (this work and Ref. 13); Rhodes Table 6.4 p.191: (2.2 ± 0.8)e-14 [Boewering 1982] | M bound, M-via |
| ar2_form (Ar* + 2Ar) | 6.0e-33 | Ueda 1984 Table V: (6.0 ± 0.5)e-33 cm⁶/s (absorption) | M |
| | band 2.8e-33 to 3.3e-32 | Rhodes Table 4.7 p.102 (Brau): 2.8e-33 (Keto 1974) ... 2.6–3.3e-32 (Kolts & Setser 1977); 1.0e-32 (Chen et al.); Chen & Payne Table I 1.1e-32 | M-via |
| ar2_rad (Ar2*) | 1/3.2 µs | Rhodes Table 4.10 p.105: Ar2*(1u) 3.2e-6 s, (0u+) 4.2e-9 s [Keto et al. 1974] | M-via |
| | band top 6.0e7 1/s | 0.25/4.2 ns + 0.75/3.2 µs, Kannari's full-mixing value | D |
| ar2_f2 (Ar2* + F2) | 5.2e-10 | Chen & Payne Table I p.155: β2 = 5.2e-10 (ref. [47]) | M-via |
| | band bottom 2.5e-10 | SRI MP-76-99 Table 2 reaction (5): Ar2* + F2 → Ar2F* + F 2.5e-10, "From decay of Ar2* and Ar2F* vs F2" | M |
| ar2_f2_arf | 0.91 | 0.5 / 0.55 (Setser's f1 over Chen & Payne's f1/f2) | D |
| K15 (Ar+ + 2Ar) | 2.2e-31 | Rhodes Table 4.8 p.104: eight measurements 5e-32 to 4.4e-31; median 2.2e-31 | M-via |
| R11p (Ar2+ + e) | 7.5e-7 (Te/300 K)^−0.67 | Rhodes Table 4.9 p.105 [Oskam & Mittelstadt 1963; Mehr & Biondi 1968] | M-via |
| K26, K108 (ion-ion) | Lee's Hoffman & Moreno forms | Lee 1989 Table I p.2055: Ar+ + F- + Ar 9.32e-6 N/(N² + 6.21e-2 N + 1.69); Ar2+ + F- + Ar 1.51e-5 N/(N² + 1.81 N + 2.25) cm³/s, N in amagat, 300 K, ref. [49] = Hoffman & Moreno SAND80-1486 | C |
| | branching to ArF* | Rokni 1977 p.2218–2219: the ion channel forms ArF* ("almost all") | M (relative) |
| FF (F + F + Ar) | 8.0e-35 | Ganguli & Kaufman 1974 via Baulch et al. 1981 p.1-22, as LIGHT/v0_2's FF cites it (KrF's first-hand check) | M-via |
| C1 cascade | prompt | Kannari 1985 Table I: Ar** + Ar → Ar* + Ar 1.0e-10 "estimated"; Petrova 2022 p.25: 6.0e-11 "k36, Chang-1978"; Setser 1977 p.-4- abstract: "More than 50% of the quenching leads directly to intermultiplet transfer to Ar(4s) states" | E, M-via, M |

## Not used as numbers, only as context

- Kannari 1985 Table I lists one ion-ion entry for all rare-gas ions with F-
  ("see text"); its text (p.4313) uses Hoffman & Moreno's profiles with a
  constant screening reduction for every recombination.
- Lee 1989 p.2056 (PDF 4) prints electron-quench values fitted to amplifier
  data (1.6, 1.9, 2.4 e-7); these are a device fit and are not used.
