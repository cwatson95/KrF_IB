# ArF (Ar/F2) rates and energetics: the research agent's report (2026-09-30)

This is the report of the research subagent that read the local library for
LIGHT_ArF, kept verbatim as it was returned (its tables and page references).
It is a reading, not a verification. The numbers LIGHT_ArF's network uses
were re-read first-hand from the PDFs before use; what was checked, and
where, is in `SOURCES_CHECKED.md` beside this file and in each row's
citation in `v0_2/network.py`. The subagent made no web requests, so
nothing here came from outside the local library.

Library shorthand used below: **LIGHT/src** = `LIGHT/docs/krf_data/sources/`;
**RKP** = `Kinetic_Modeling/Recent_KrF_Papers/`; **docs/RKP** =
`Kinetic_Modeling/docs/Recent_KrF_Papers/`; **Rhodes** =
`docs/Excimer_Lasers_Topics_Rhodes.pdf`. Page numbers: printed page / PDF
page.

**Key finds.**
1. Chen & Payne, IEEE JQE 15, 149 (1979) is in LIGHT/src as
   `Kinetics_of_proton-excited_Ar-F2and_Kr-F2mixtures.pdf`. It is the main
   Ar/F2 measurement paper.
2. Setser's 1977 ERDA report gives the ArF* branching of Ar(3P2) + F2 as
   about 0.5.
3. Ueda's two pulse-radiolysis papers carry measured Ar* and Ar2F*
   constants.
4. Jancaitis Chapter V has a full Ar2F* measurement.
5. Duzy/Avco 1980 Table III has measured ArF* + e values.

## 0. The coordinator's items 1–8, one line each

1. **ArF(B) energy and the C state.** Neither is measured in the held
   sources. Krauss & Mies (Rhodes Table 2.2, p.14/PDF 24) give only
   R_e(B) = 2.4 Å for ArF (calculated, their ref 2.39); D_e(B) is blank.
   Their p.15 (PDF 25) says the ab initio ArF* B and C curves are split by
   "about 0.25eV near R_e" (calculated), and that accurate ab initio results
   put the ²Σ₁/₂ state "always lower in energy than the ²Π₃/₂" (calculated).
   From the 193.3-nm band (Chen & Payne's peak is 1933 Å), B lies at least
   6.41 eV above Ar + F (7.21 eV above Ar + F2), so D_e(B) ≤ 5.95 eV, with
   IP(Ar) − EA(F) = 12.359 eV (derived). Jancaitis's model value
   (UCRL-53465, Table I p.85/PDF 93) is ArF* = 7.24 eV (model). C lifetime,
   calculated only: Kannari (RKP/kannari1985 Table I p.4312/PDF 5) gives
   C → Ar + F + hν = 2.27×10⁷ s⁻¹ [Hay & Dunning]; Lee 1989 (docs/RKP Table
   I p.2056/PDF 4) prints "2.08 × 10⁻⁷ s⁻¹", an exponent misprint for 10⁷.
   Nothing held gives a C position, a B↔C mixing rate or an equilibrium
   constant.
2. **ArF(B) radiative lifetime.** Not measured in any held source. McCown &
   Greene (LA-UR-90-291, final page) still describe it in 1990 as a planned
   measurement, "less than 10 ns". The values in use are all calculations:
   4 ns (Dunning & Hay, used by Rokni and by Chen & Payne), 4.2 ns (Dunning
   & Hay, JCP 69, 134 (1978), per Jancaitis p.97), 2.50e8 s⁻¹ (Kannari) and
   2.38e8 s⁻¹ (Lee).
3. **ArF* yield of Ar* + F2: about 0.5, measured.** Setser ERDA
   COO-2807-11 (1977), printed p.-6-/PDF 8, flowing afterglow, Ar(3P2):
   "Molecular F₂ is still the best donor, but the branching ratio is ∿ 0.5
   rather than unity." Its published form appears to be Kolts, Velazco &
   Setser, J. Phys. Chem. 82, 1766 (1978) (not held); Ueda 1984 Table III
   (p.365/PDF 9) lists "4.0 × 10⁻¹⁰, J. H. Kolts et al." for Ar* + F2 →
   ArF* + F, which is 0.53 of Velazco's 7.5e-10. Chen & Payne, JQE 1979
   Table I (p.155/PDF 7): f₁/f₂ = 0.55 ± 0.25, a relative measurement (ArF*
   fraction from Ar(3P2) + F2 over that from Ar2*(1u) + F2). Chen & Payne
   p.159/PDF 11: Ar(3P1) and Ar(1P1) + F2 give ArF* at "< 25 percent".
4. **Ar* + 2Ar → Ar2* + Ar, measured at 1–4 atm:** (6.0 ± 0.5)×10⁻³³ cm⁶/s
   (Ueda et al., Rev. Laser Eng. 12, 357 (1984) Table V p.365/PDF 9; same
   value in Rev. Laser Eng. 11, 576 (1983) Table III p.584/PDF 9). Other
   measurements span 2×10⁻³³ to about 2.4×10⁻³² (below). Chen & Payne use
   1.1×10⁻³² (their Table I) and quote β = 11 ± 2 torr⁻² s⁻¹ (p.152,
   Thonnard & Hurst), which converts to 1.06e-32 at 300 K (derived).
   **Penning:** no measurement in held sources. The 5.0×10⁻¹⁰ cm³/s in
   Kannari Table I p.4312 cites Leichner & Ericson, PRA 9, 251 (1974),
   apparently a Kr paper (McCusker lists it with Kr2* lifetime data in
   Rhodes Table 3.1, p.51). Treat it as an estimate.
5. **Ar2* radiative lifetime, measured:** Keto, Gleason & Walters, PRL 33,
   1365 (1974), via McCusker (Rhodes Table 3.1 p.51/PDF 61) and Brau (Table
   4.10 p.105/PDF 115): ¹Σu⁺ (0u⁺) 4.2 ± 0.13 ns; ³Σu⁺ (1u) 3.2 ± 0.3 µs.
   Other ³Σu values: 2.8 µs (Thonnard & Hurst 1972), 3.7 µs (Boucique &
   Mortier 1970), 3.4 µs (Colli 1954), 3.22 µs (Oka 1974), 4.0 ± 2.0 µs.
6. **Ar⁺ + 2Ar → Ar2⁺ + Ar, measured:** 2.07×10⁻³¹ and 2.3×10⁻³¹ (Liu &
   Conway 1974, 1975) and 2.5×10⁻³¹ (Smith & Cromey 1968), cm⁶/s (Brau
   Table 4.8, p.104/PDF 114). **Dissociative recombination, measured:**
   7.5×10⁻⁷ (Te/300 K)^−0.67 cm³/s (Oskam & Mittelstadt 1963; Mehr & Biondi
   1968; Brau Table 4.9, p.105/PDF 115). Shiu & Biondi 1978 (argon) is not
   held.
7. **ArF* + e, measured:** 1.6 and 3.0 ×10⁻⁷ cm³/s at 44.4 and 10.1 A/cm²,
   0.25 atm (Duzy et al., Avco 1980, Table III p.-23-/PDF 38). Lee p.2054
   (PDF 2) says Trainor & Jacob (APL 37, 675 (1980)) "reported it to be 1.6
   × 10⁻⁷ cm³/s", derived "by using a theoretical rate constant such as the
   electron dissociative attachment to F₂". Ar⁺ + F⁻ and Ar2⁺ + F⁻ have no
   absolute measurement ("has not been directly measured", Brau p.98/PDF
   108). Measured, relative only: Ar2⁺ + F⁻ gives ArF* with branching near
   unity; "almost all the Ar₂⁺ forms ArF*" (Rokni, Jacob & Mangano, PRA 16,
   2216 (1977), p.2219/PDF 4, Fig. 5).
8. **Ar⁺ + F2 charge transfer.** Nothing in any held source or model
   table. IP(F2) = 15.69 eV (Petrova Table A1, model value) is below IP(Ar)
   = 15.76 eV, so it is near-resonant and exothermic by about 0.07 eV
   (derived). Carry it as a band with default 0.

## 1. Process table

Units cm³/s, cm⁶/s, s⁻¹. **Status:** M = measured (primary source read
first-hand); M-via = a measurement quoted by a source read here; C =
calculation; E = estimate or guess; MT = model-table value; D = derived by
arithmetic.

**Lifetime convention.** Rokni, and Chen & Payne, measured k·τ
(Stern-Volmer). Their absolute k assume τ(ArF) = 4 ns. Lee's rows equal the
same k·τ divided by 4.2 ns: 7.6e-18/4.2 ns = 1.81e-9; 3.6e-20/4.2 ns =
8.6e-12; 2.1e-39/4.2 ns = 5.0e-31 (D). **Carry k·τ, not k.**

### A. Metastable and excited-atom chemistry

**Ar(3P2) + F2, total quench**
- 7.5×10⁻¹⁰, σ = 132 Å² (Velazco, Kolts & Setser, JCP 65, 3468 (1976)
  Table II p.3472/PDF 6, LIGHT/src/velazco1976.pdf). M, flowing afterglow
  at about 0.8 torr, "< ±10%" reproducibility; p.3469 notes wall loss of F2
  "would tend to make the measurements lower limits".
- 9.4×10⁻¹⁰ (Chen & Payne Table I p.155 and Table II p.160, from their ref.
  47). M, proton beam, N2(C) tracer.
- (7.4 ± 1.0)×10⁻¹⁰ by Ar* absorption and (6.4 ± 0.3)×10⁻¹⁰ by ArF*
  fluorescence decay, 0.5 atm with 30% Kr (Ueda 1984 Table V p.365). M. The
  fluorescence number tracks Ar* removal; it is not an ArF* yield.
- (1.1 ± 0.1)×10⁻⁹ by absorption at 2.5 atm, and 0.9×10⁻⁹ by fluorescence
  (Ueda 1983 p.583–584, Table IV). M.
- Tg law: none measured.
- Other: Brau Table 4.4 (9×10⁻¹⁰, Chen private communication). Ueda's
  comparison lists give 0.85e-9 ("Nakano"; SRI Table 2 prints 7.5e-10) and
  0.40e-9 (Morgan & Szoke, model). Models: 7.5e-10 (Rokni Table I, SRI,
  Kannari, Lee, Jancaitis) and 1.0e-9 (Petrova).

**Ar(3P2) + F2 → ArF* + F (yield)**
- ≈ 0.5 (Setser COO-2807-11, p.-6-/PDF 8). M.
- 4.0×10⁻¹⁰, "J. H. Kolts et al." = JPC 82, 1766 (1978) (Ueda 1984 Table
  III p.365). M-via.
- f₁/f₂ = 0.55 ± 0.25 (Chen & Payne Table I p.155). M, relative.
- Tg law: none.
- Notes: Brau Table 4.4 prints "≤ 0.6" from Chen private communication.
  Kannari prints 7.5e-10 to ArF* plus 3.1e-10 to Ar + 2F (Table I p.4311),
  which sums above the total he cites. Chen 2025 R29 uses 4e-10 cm³/s (from
  Akashi). With f₁ = 0.5, f₁/f₂ gives an Ar2* ArF* fraction f₂ ≈ 0.91,
  range 0.62–1 (D).

**Ar(3P1) + F2 and Ar(1P1) + F2:** 8.9 and 12.9 ×10⁻¹⁰ (text says 13); σ =
156 and 225 Å²; ArF* yield "< 25 percent" (Chen & Payne Table II p.160 and
p.151, p.159). M. These are the resonance 4s states.

**Ar(3P0) + F2:** no separate measurement (Velazco p.3469: 3P0 is about 1/6
of the metastables and was not measured).

**Ar(4p) = Ar** + F2:** no measurement. Model values: Kannari 4.7e-10 to
ArF* plus 3.1e-10 dissociative, "analogous to Ar*"; Lee 7.5e-10; Jancaitis
7.5e-10; Petrova 1.0e-9 ("estimated").

**Ar* + 2Ar → Ar2* + Ar**
- (6.0 ± 0.5)×10⁻³³ (Ueda 1984 Table V; Ueda 1983 Table III p.584). M,
  696.5-nm absorption, 1–4 atm.
- 1.1×10⁻³² (Chen & Payne Table I, citing Chen, Payne, Hurst & Judish JCP
  65, 4028 (1976)). M-via.
- Tg law: none held.
- Other: Brau Table 4.7 (p.102/PDF 112): 6.6e-33 (Boucique & Mortier
  1970), 9e-33 (Colli 1954), 1.14–1.70e-32 (Ellis & Twiddy 1969), 1.3e-32
  (Futch & Grant 1956), 1.2e-32 (Thonnard & Hurst 1972), 1–1.2e-32
  (Leichner), 1.2e-32 (Hurst 1974), 6e-33 (LeCalve & Bourene 1973), 8.7e-33
  (Phelps & Molnar 1953), 2.6–3.3e-32 (Kolts & Setser 1977), 2.8e-33 (Keto
  1974). Ueda 1983 Table III adds 17e-33 (Wieme), 4.1e-33 (Firestone),
  2e-33 (Suzuki) and Sadeghi, printed "23.  8×10⁻³³" (reads as 2.38e-32;
  ambiguous spacing). Models: SRI Table 2 1e-32 (±30%); Kannari 1.14e-32;
  Lee 1.1e-32.

**Ar* + Ar* → Ar⁺ + Ar + e (Penning/associative):** no measurement held.
E/MT: 5.0e-10 (Kannari Table I p.4312, ref 86 = Leichner & Ericson 1974);
5×10⁻¹⁰ (McCusker Table 3.7 p.76/PDF 86, same ref); 5.0e-10 (Kolokolov
1994) and 6.2e-10 (Lee & Lieberman 1995, a global model), both from Petrova
p.26 (PDF 33). Brau p.105: "generally occur at every collision". The
branching between Ar⁺ + Ar + e and Ar2⁺ + e is not measured.

**Ar** + Ar → Ar* + Ar (4p relaxation):** 1.0e-10 "estimated" (Kannari);
1.0e-10 (Lee, "Estimated from [29]"); 6.0e-11 (Petrova p.25,
"Chang-1978"). Setser COO-2807-11 p.-4- abstract: "More than 50% of the
quenching leads directly to intermultiplet transfer to Ar(4s)" (M,
qualitative). The measurement is Chang & Setser, JCP 69, 3885 (1978), not
held.

### B. The ArF* upper level

**ArF(B) → Ar + F + hν, 193 nm:** τ = 4 ns (Rokni Table II footnote; Chen
& Payne Table I, "ab initio"); 4.2 ns (Jancaitis Table II p.97, citing
Dunning & Hay JCP 69, 134 (1978)). C. Others: 2.50e8 (Kannari, Hay &
Dunning); 2.38e8 (Lee; Chen 2025 R37 2.381e8); 2.67e8 (Petrova p.15/PDF 22,
"Dunning-1978", = 3.7 ns); Petrova p.5 also prints "A = 2.67 × 10⁻⁸ s⁻¹
(lifetime τ = 2.7 ns)", internally inconsistent. SRI Table 2 prints 3.3e7,
a fitted yield ratio against a guessed ArF* + F2; ignore it.

**ArF* + F2 → products:** kτ = (7.6 ± 0.7)×10⁻¹⁸ cm³, k = 1.9e-9 at 4 ns
(Rokni, Jacob & Mangano Table II p.2218/PDF 3, RKP/"Dominant
formation…ArF and KrF lasers.pdf"). M. The same paper reports an F2
half-pressure of 4.0 torr at 150 torr Ar, which gives kτ = 7.77e-18 (D).
τ₃β₃ = 7.5×10⁻¹⁸ cm³, β₃ = 1.9×10⁻⁹ ± 25% (Chen & Payne p.155 and Table
I). M. Tg law: none measured (Jancaitis's model applies (Tg/300)^0.5).
Others: 8.0e-18 (Brau Table 4.5, Chen private, τ = 4.2 ns); 1.8e-9 (Lee);
2.0e-9·A (Jancaitis p.95); 1e-9 "Guessed" (SRI). Kannari prints 1.9e-9
unreduced; his footnote a applies ×0.418 inside his model.

**ArF* + Ar → 2Ar + F:** kτ = (3.6 ± 1)×10⁻²⁰ cm³, k = 9×10⁻¹² at 4 ns
(Rokni Table II). M. Rokni p.2218: Ar quenching is negligible at ≤ 200
torr. Others: 8.6e-12 (Lee); 1.3e-11·A (Jancaitis p.95); 9e-12 (Kannari;
Chen 2025 R34).

**ArF* + 2Ar → Ar2F* + Ar:** kτ = (1.6 ± 0.3)×10⁻³⁹ cm⁶, k = 4e-31 (Rokni
Table II), M. τ₃β₃′ = 2.1×10⁻³⁹ cm⁶, β₃′ = 5.3×10⁻³¹ ± 25% (Chen & Payne
p.155), M. kτ = 3.1×10⁻³⁹ cm⁶, ±25% (Jancaitis ch. V p.150/PDF 158), M. Tg
law: no measurement held. Others: 2.2e-39 (Brau Table 4.6, Chen private,
4.2 ns); Lee 5.0e-31 "at 300 K"; Kannari 4.0e-31 (footnote b: Tg per Shui);
Jancaitis model R×3.2e-31·A (M = Ar R = 1, Kr R = 0.8, He R = 0.4); SRI
2×10⁻³³ (factor-2 uncertainty; SRI called the reaction minor). The Petrova
form is 1.62e-32 + 2.245e-21·Tg^−3.9 ("k28, Shui-1979, Chen-1978",
p.25/PDF 32); it equals 5.07e-31 at 300 K and falls ×0.097 by 600 K (D).

**ArF* + e → Ar + F + e:** 1.6 and 3.0 ×10⁻⁷ (Duzy et al., Avco 1980 Table
III p.-23-/PDF 38, LIGHT/src). M; the KrF average there is 2.0 ± 0.2e-7.
Trainor & Jacob APL 37, 675 (1980), 1.6e-7 via Lee p.2054: M-via. Others:
Kannari 1.6e-7 (ref 45); Jancaitis 1e-7 (Mangano, private communication);
Petrova 2.4e-7. Lee's 1.6, 1.9 and 2.4 e-7 are fits to amplifier data, a
device fit. **Do not use them.**

**ArF* + F → Ar + 2F:** no measurement. Petrova lists 1.0e-12 as
"Rokni-1977" but also as "k=0, not known"; Rokni has no such row.

**ArF(C) → Ar + F + hν:** 2.27e7 (Kannari, Hay & Dunning; C); 2.08e7 (Lee,
printed 10⁻⁷); 2.1e7 (Petrova). Nothing measured.

**Stimulated emission cross section, ArF B–X:** 2.75×10⁻¹⁶ cm² (Lee Table
I p.2056 and text; C, from Dunning & Hay's lifetime and a measured 15 Å
bandwidth, Kudryavtsev & Kuzmina); 1.9×10⁻¹⁶ cm² (Petrova p.24/PDF 31,
"McGeoch-2017; Dunning-1978"; C/MT; Petrova attributes the same number to
Lee, which is wrong); 4e-16 cm² (Chen 2025 Table III R38, Akashi; MT). None
is measured.

### C. Ar2* and Ar2F*

**Ar2* → 2Ar + hν, 126 nm:** fast 4.2 ± 0.13 ns, slow 3.2 ± 0.3 µs (Keto
1974, via McCusker Table 3.1 and Brau Table 4.10). M-via. Model forms:
Kannari 6.0e7 (ref 88 = Keto), the 1:3 singlet:triplet equilibrium
average, 0.25/4.2 ns + 0.75/3.2 µs = 5.98e7 (D), which assumes fast
electron mixing (McCusker eqs. 3.6–3.7, p.50–51, k_e ≈ 2e-6 cm³/s); SRI
3.8e6 ("Effective value in the presence of F₂, 1–9 atm Ar", ±50%); Chen &
Payne τ₂ = 3 µs; Lee prints "6.0 × 10⁻⁷ s⁻¹", an exponent misprint;
Jancaitis 3,100 ns and 4.2 ns.

**Ar2*(1u) + F2, total:** 5.2×10⁻¹⁰, σ = 105 Å² (Chen & Payne Fig. 5 at 605
torr, Tables I and II). M. Brau Table 4.4 lists 5.2e-10 (Chen private) and
3×10⁻¹⁰ "(SRI)"; SRI's Table 2 prints 2.5e-10, for the Ar2F* channel.

**Ar2* + F2 split:** Chen & Payne p.156: "f₂′/f₃ … should be between 0 to
0.3" — direct Ar2* + F2 → Ar2F* + F is minor next to ArF* + 2Ar; ArF* is
the main product (f₂ ≈ 0.9, D). M, relative. SRI Table 2 (5): Ar2* + F2 →
Ar2F* + F = 2.5e-10, "From decay of Ar₂* and Ar₂F* vs F₂" (disputed).
Models: Kannari 3.0e-10 (ArF*) + 2.5e-10 (Ar2F*); Lee 5.2e-10 + 2.5e-10
(Lee's ref. 52 calls the 2.5e-10 "assumed to be analogous to Ar₂* + NF₃";
his two channels sum above the measured total).

**Ar2* + F → ArF* + Ar:** no measurement. SRI p.65 conditional: "3-10 ×
10⁻¹⁰" if the reaction is fast. Models: Kannari 3.0e-10; Petrova 1.0e-9
("k=0, not known").

**Ar2* + Ar2* → Ar2⁺ + 2Ar + e:** no Ar measurement held. 5×10⁻¹⁰ (Kannari
ref 87 and McCusker Table 3.7 ref 3.3: Lorents & Olson SRI reports
1972–73, unpublished). McCusker p.52 (PDF 62) quotes 3.5×10⁻¹⁰ (Gerardo &
Johnson 1973) and 8×10⁻¹¹ (Zamir et al. 1975) for excimer-excimer
annihilation, written as Ar2* in his eq. 3.10; the gas measured is not
stated there, likely Xe.

**Ar2* + e (quench and superelastic):** 1.0e-7 "estimated" (Kannari);
1×10⁻⁹ (McCusker Table 3.7); 5.0e-8 (Petrova, "estimated based on
Neeser-1997").

**Ar2F* → 2Ar + F + hν, about 285–290 nm, all M:** 185 ± 25% ns (Chen &
Payne Table I); 205 ± 20 ns (Jancaitis ch. V Table I p.158/PDF 166); 219 ±
15 ns (Marowsky et al., IEEE JQE 18, 898 (1982), via Jancaitis and Ueda);
220 ns (Ueda 1983 Table II p.582/PDF 7; Ueda 1984 Table V); 230 ns
(Boewering et al., JCP 76, 3524 (1982), via Ueda). Theory: 132 ns (Wadt &
Hay, JCP 68, 3850 (1978), via Jancaitis), C. Models: Kannari 5.4e6; Lee
4.4e6; Petrova 4.3e6. SRI 2×10⁸ is obsolete.

**Ar2F* + F2, all M:** 2.1×10⁻¹⁰ ± 25%, τ₄β₄ = 4.5e-17 cm³ (Chen & Payne
p.154–155); (1.79 ± 0.06)×10⁻¹⁰ (Jancaitis ch. V Table I p.158); (2.05 ±
0.06)×10⁻¹⁰ (Marowsky, via Jancaitis); (1.9 ± 0.1)×10⁻¹⁰ (Ueda 1984 Table
V). Others: 3.9e-17 cm³ kτ (Brau Table 4.5, Chen private, τ = 185 ns);
Morgan & Szoke 4e-10 (model); Kannari and Lee 2.0e-10; SRI 1e-9
"Guessed".

**Ar2F* + Ar:** < 10⁻¹⁴ (Jancaitis ch. V Table I; Marowsky, same bound). M,
upper bound. Others: Lee 2.2e-14 (Boewering, ref 58); Kannari 1.0e-10
"estimated", 10⁴ above the measured bound.

**Ar2F* + e:** 4×10⁻⁸ (Jancaitis p.152/PDF 160) against 2.×10⁻⁸ in his own
Table II p.96. E, derived from an assumed attachment rate of 5e-9. Kannari
1.0e-7 (ref 17); Lee 1.0e-7 "Guessed".

**Ar2F* + F:** no direct measurement. Jancaitis p.153 infers quenching by F
at 0.5–1 × the F2 rate (via LIGHT/v0_2/SOURCES.md §5, not re-read here).

### D. Ions and fluorine

**Ar⁺ + 2Ar → Ar2⁺ + Ar:** 2.07e-31 (Liu & Conway 1974), 2.3e-31 (Liu &
Conway 1975), 2.5e-31 (Smith & Cromey 1968) (Brau Table 4.8 p.104).
M-via. Tg law: none measured (T^−3/4 theory, Pérez-Ríos & Greene 2018, as
the KrF job). Others: 4.4e-31 (Knewstubb 1968), 1.46e-31 (Kritschmer 1963),
5e-32 (Kebarle 1966), 1.9e-31 (Bhattacharya 1970), 3.0e-31 (Cronin &
Sexton 1968). Models: Rokni 2.5e-31 (McDaniel); Lee 2.5e-31; Kannari
2.5e-31; Petrova 2.1e-31 and 2.07e-31; Chen 2025 R27 2.5e-43 m⁶/s.

**Ar2⁺ + e → Ar(*) + Ar:** 7.5×10⁻⁷ (Te/300 K)^−0.67 (Brau Table 4.9).
M-via. Products: McCusker p.49–50, the dissociating curves "lead to argon
atoms in 4p states". Models: Kannari 0.6×10⁻⁶ (Te/300)^−0.66 to Ar* plus
1.1e-7 to Ar**; Lee 1.5×10⁻⁹ (Tg/Te)^0.67 to Ar** (= 6.85e-8
Te[eV]^−0.67 against Brau's 6.48e-8, D); Petrova 7.3e-8 Te^−0.67 to Ar**
("Mehr-1968"); McCusker 1×10⁻⁶. No gas-temperature law measured.

**Ar⁺ + F⁻ (+M) → ArF*; Ar2⁺ + F⁻ (+M) → ArF* + Ar:** no absolute
measurement. Measured, relative: Ar2⁺ + F⁻ gives ArF* with branching near
1 (Rokni p.2219, Fig. 5; p.2223: "formed via the ion channel with unit
branching"). Models: Rokni Table I 10⁻⁷ + 10⁻⁷ p (p < 1 atm), 1.1×10⁻⁶
above 1 atm, "assumed" equal for Ar⁺ and Ar2⁺; Lee (Hoffman & Moreno)
9.32e-6 N/(N² + 6.21e-2 N + 1.69) for Ar⁺ and 1.51e-5 N/(N² + 1.81N +
2.25) for Ar2⁺, N in amagat (3.39e-6 and 2.98e-6 at 1 amagat, D); Petrova
(Flannery) 2.4N/(N² + 1.7)×3.417e-6 and 2.65N/(N² + 1.7)×3.096e-6,
identical to the KrF K26/K108 constants; Chen 2025 9.1e-13 and 9.26e-13
m³/s, constant.

**F⁻ + F2 (+Ar):** 2×10⁻¹⁰, effective two-body at 600 torr Ar, products
unknown (Chen & Payne p.156–157, Fig. 20). M.

**F + F + Ar → F2 + Ar:** (8.0 ± 0.6)×10⁻³⁵ at 295 K (Ganguli & Kaufman
1974, via Baulch 1981 p.1-22; M-via, via LIGHT/v0_2/SOURCES §3); 2.8e-34
at 300 K (Cohen & Bott evaluation, "uncertain by at least an order of
magnitude"); 1.0e-33 (Kannari, citing Johnson & Hunter 1980, "Estimated").
Convention d[F2]/dt = k[F]²[M].

### E. 193-nm absorbers (none measured first-hand in the library)

F2 1.5×10⁻²¹ cm² (Lee; estimated from the absorption spectrum over 2100 Å,
Steunenberg & Vogel) vs 1.0×10⁻²¹ (Petrova, Argüello 1995). F⁻ 1.0×10⁻¹⁸
(Lee, Robinson & Geltman 1967) vs 6.5×10⁻¹⁸ (Petrova, Mandl 1971). Ar*
(4s) 4.0×10⁻²⁰ (Lee) and 3.5×10⁻²⁰ (Petrova), both "Duzy-1980", a
calculated curve (Duzy Fig. 14 p.-39-/PDF 54), C. Ar** (4p) 3.0×10⁻¹⁸
(Lee, Petrova; Duzy Fig. 15 p.-40-), C. Ar2* 6.3×10⁻¹⁹ (Lee and Petrova,
citing Sauerbrey et al., J. Phys. B 19, 2279 (1986), not held). Ar2⁺
1.0×10⁻¹⁹ "Guessed" (Lee). Ar2F* 5.0×10⁻¹⁹ "Guessed" (Lee); McCown &
Greene (LA-UR-90-291, Fig. 4) measured only a relative Ar2F spectrum
(bands at 295 and 220 nm). ArF* self-absorption: nothing held.

### F. Electron-impact sets used by other models

Rokni/Lee 26.4 eV per ion pair; Lee 93 eV per Ar* (Lorents) plus Te fits
from the MSNW 1980 report. Kannari 26.2 eV per ion pair, Ar* = 0.28 ×
ionization. Chen 2025: SIGLO through BOLSIG+ (SIGLO lumps Ar excitation
into one level). NRL Orestes: its own Boltzmann solver (Petrov 2017).
**Ours:** the IST-Lisbon file `LoKI-MC/Code/Input/Argon/Ar_LXCat.txt`, each
4s, 4p, 3d/5s, 5p, 4d and 6s level and ionization from the ground state
only. It has no stepwise ionization from Ar* (literature fits disagree ×100
at Te = 1 eV: Lee 5.18e-6 Te^0.67 exp(−6.27/Te) against Kannari 6.80e-9
Te^0.67 exp(−4.20/Te)) and no Ar2* data.

## 2. Energetics (eV above ground-state Ar + F2)

- Ar 4s: 3P2 11.548, 3P1 11.623, 3P0 11.723, 1P1 11.828 (IST-Lisbon
  thresholds). 4p 12.907–13.480; 3d/5s 13.84–14.30; 5p 14.46–14.74. IP(Ar)
  15.760.
- Ar2⁺: 15.760 − 1.23 = 14.53 (D). D_e(Ar2⁺) = 1.23 eV (photoionization,
  Ng et al., via Krauss & Mies Table 2.1 p.9/PDF 19), M-via; Wadt's
  calculation 1.19. Models: Petrova 14.53, Jancaitis 14.50.
- Ar2* (³Σu⁺): about 10.87, D_e = 0.68 below Ar(3P2) + Ar (D); D_e 0.68 eV
  and R_e 2.43 Å calculated (Saxon & Liu, via Krauss & Mies Table 2.1 and
  p.11/PDF 21), C. Models: Jancaitis 11.06 (³Σ) and 11.135 (¹Σ); Petrova
  10.60; REBUILD 10.4.
- Ar2* emission 126 nm (SRI) = 9.84 eV; Chen & Payne's continua peak at
  1270 Å (1u, 9.76 eV) and 1100 Å (0u) (p.150).
- F 0.80 (D0(F2)/2); F⁻ 0.80 − 3.401 = −2.60.
- ArF(B): at least 7.21 (D, flat X state: 0.80 + 6.414 eV at 193.3 nm).
  The ion-pair asymptote lies 12.359 above Ar + F, so D_e(B) ≤ 5.95, not
  measured. R_e(B) 2.4 Å (Krauss & Mies Table 2.2), C. Models: Jancaitis
  7.24 (about 0.03 eV of X repulsion); REBUILD 7.52; Petrova 8.30 (not
  credible).
- ArF(C): not measured; calculated B/C splitting about 0.25 eV near R_e
  (Krauss & Mies p.15). Petrova's 120 cm⁻¹ is the KrF value copied over.
- Ar2F* (4²Γ): about 0.65 below ArF* + Ar, so about 6.56 if B is 7.21 (D).
  0.65 ± 0.08 eV (Huestis Table 6.3); ΔE = 0.6 eV for Ar2F* and Kr2F* (Wadt
  & Hay, via Rokni p.2223). C. Models: Jancaitis 6.10; Petrova 6.00;
  REBUILD 5.45.
- Photon bands: ArF B–X 1933 Å peak (Chen & Payne) = 6.414 eV, 193.0 nm =
  6.424 eV; Ar2F* 2840 Å (Chen & Payne) = 4.366 eV, 290 nm (SRI) = 4.275
  eV, Rokni "290-300 nm"; Ar2* 126 nm = 9.84 eV. ArF B–A and C–A: no
  position in the held sources.

## 3. No measurement: candidates for bands with default 0

- Ar* + Ar* → Ar⁺ + Ar + e: 6.2e-10 (Lee & Lieberman 1995, model); 5.0e-10
  (Leichner & Ericson, a Kr paper, via Kannari).
- Ar* + Ar** and Ar** + Ar**: 5.0e-10, "estimated".
- Ar2* + Ar2*: 5e-10 (Lorents & Olson, unpublished).
- Ar2* + F → ArF* + Ar: 1.0e-9 (Petrova, "not known"); SRI conditional
  3–10e-10.
- Ar** (4p) + F2: 1.0e-9 (Petrova, estimated); analog values 7.5e-10.
  Ar** relaxation by 1 atm Ar is about 1.5e9 s⁻¹ (6e-11 × 2.45e19), so at
  0.3 % F2 F2 would take only about 4 % of Ar**.
- Ar(3P0) + F2: no separate value.
- ArF(C) radiative rate, B↔C mixing, gap and quench: calculated only.
  Most measurement-based ArF models (Rokni, Chen & Payne, SRI, Jancaitis)
  lump B and C.
- ArF* + F: 1e-12 (Petrova, "k=0, not known").
- Ar2F* + F: 0.5–1 × the Ar2F* + F2 rate (Jancaitis's inference).
- Ar2F* + e: 1e-7 (guessed/estimated); Jancaitis's own estimates 2–4e-8.
- Ar2F* + Ar: measured upper bound 1e-14; band 0 to 2.2e-14 (Lee/Boewering).
- Ar2* + e: 1e-7 (estimated).
- Ar⁺ + F⁻ and Ar2⁺ + F⁻ absolute rates near 1 atm: theory only, about
  3e-6 at 1 amagat (Flannery; Hoffman & Moreno).
- Ar⁺ + F2 → F2⁺ + Ar: no value held; capture limit about 6e-10 (Langevin
  with α(F2) ≈ 1.4 Å³, the subagent's estimate, α not verified).
- Gas-temperature laws: none measured for ArF* + 2Ar, Ar* + 2Ar, Ar2⁺ + e,
  ArF* + F2 or Ar(3P2) + F2. Theory/model forms: Shui/NRL Tg^−3.9 for the
  trimer; Jancaitis's blanket (Tg/300)^0.5; T^−3/4 for Ar⁺ + 2Ar.
- Stepwise e + Ar* and e + Ar** ionization: absent from IST-Lisbon; fits
  differ ×100.

## 4. Open conflicts between sources

1. Ar(3P2) + F2 total: 6.4–11 ×10⁻¹⁰ across measurements (Velazco 7.5, a
   lower limit by his own caveat; Chen & Payne 9.4; Ueda 1984 7.4 ± 1.0 and
   6.4 ± 0.3; Ueda 1983 11 ± 1 by the same absorption method, not
   reconciled by the 1984 paper).
2. ArF* yield of Ar* + F2: measured about 0.5 (Setser; Kolts 4.0e-10), but
   every measurement-era model (Rokni, SRI, Lee, Jancaitis, Kannari, and
   nathaniel's arf.toml R29) books the full 7.5e-10 to ArF*.
3. Brau Table 4.4's "3×10⁻¹⁰ (SRI)" Ar2* + F2 row does not match SRI Table
   2, which prints 2.5e-10 for Ar2* + F2 → Ar2F* + F (its 3e-10 is the Kr2*
   row). Brau's "< 0.3" corresponds to Chen & Payne's f₂′/f₃ < 0.3, the
   direct Ar2F* share relative to trimer formation from ArF*, not an ArF*
   branching; Chen & Payne's data imply Ar2* + F2 gives mostly ArF* (f₂ ≈
   0.9).
4. Where Ar2F* comes from: SRI attributes it to Ar2* + F2 and calls ArF* +
   2Ar minor at 2e-33; Rokni, Chen & Payne and Jancaitis measure ArF* + 2Ar
   at kτ = 1.6–3.1e-39 cm⁶ and find the direct Ar2* route minor. Jancaitis
   shows trimer formation "independent of fluorine donor, mixture ratio, and
   e-beam current density" (p.149).
5. ArF* + 2Ar magnitude: kτ spread ×1.9 (1.6, 2.1, 2.2, 3.1 e-39 cm⁶); the
   absolute k depends on the unmeasured τ_B; the Tg law is unmeasured, and
   Jancaitis's model assumes the opposite sign.
6. Ar* + 2Ar: ×10 spread (2e-33 to 2.4e-32); Ueda's high-pressure 6.0e-33
   sits below the low-pressure afterglow cluster at about 1.1e-32.
7. Ar2F* + F2: Jancaitis (1.79 ± 0.06)e-10 and Marowsky (2.05 ± 0.06)e-10
   disagree beyond their stated errors; Ueda 1.9 and Chen 2.1 ± 25% lie
   between.
8. Ar2F* lifetime: 185–230 ns measured, against 132 ns calculated.
9. Ar2F* + Ar: measured bound < 1e-14, Kannari "estimated" 1.0e-10.
10. Ar2F* + e: Jancaitis 4e-8 (text) vs 2e-8 (table).
11. Effective Ar2* radiative rate: 6.0e7 (Kannari, full singlet/triplet
    mixing) vs 3.8e6 (SRI, effective) vs the 3.1–3.2 µs triplet (Chen &
    Payne; Jancaitis).
12. ArF(B) radiative rate: all calculated (2.38e8, 2.50e8, 2.67e8 s⁻¹).
13. ArF(B) and Ar2F* formation energies scatter by more than 1 eV across
    models; only Jancaitis's pair is near the spectroscopic lower bound.
14. Ar2⁺ + e products: 4p (McCusker); Kannari splits about 85/15 4s/4p;
    Lee and Petrova send all to Ar**. Magnitudes agree within 6 %.
15. F⁻ photodetachment at 193 nm: 1.0e-18 (Lee) vs 6.5e-18 (Petrova).
16. σ_SE(ArF): 1.9, 2.75 and 4.0 e-16 cm², all calculated or model values.
17. ArF* + e: Duzy's two points differ ×1.9 with current density; Trainor
    & Jacob's value depends on a calculated attachment rate; Lee's 2.4e-7
    for Ar/F2 is a device fit (excluded).
18. Print and attribution errors not to copy: Lee Table I exponent signs,
    "Ar₂R*" for Ar2F*, the σ_SE row's [50], the F2 absorption's [57], the
    Ar2⁺ photodissociation row's "+ e", the Ar2F* radiative row's "+ e";
    Petrova's σ_SE attribution to Lee; Ueda 1984 Table III's "H. L. Chen
    7.5e-10" and "Nakano 8.5e-10"; Ueda 1983 Table IV's units "cm² s⁻¹".

## 5. Not held and worth obtaining (priority order)

1. Kolts, Velazco & Setser, J. Phys. Chem. 82, 1766 (1978): the published
   ArF* branching (about 0.5, 4.0e-10) and the ArF B/C ratio against
   pressure (per Setser COO-2807-17).
2. C. H. Chen, M. G. Payne & J. P. Judish, JCP 69, 1625 (1978).
3. Chen, Judish & Payne, J. Phys. B (Ar-Xe and Ar-F2 energy transfer; the
   source of 9.4e-10 and 5.2e-10).
4. Boewering, Sauerbrey & Langhoff, JCP 76, 3524 (1982).
5. Marowsky et al., IEEE JQE 18, 898 (1982).
6. R. S. F. Chang, JCP 78, 2210 (1983), ArF*(B) in Ne.
7. Trainor & Jacob, APL 37, 675 (1980).
8. Shiu & Biondi 1978, Ar2⁺ recombination.
9. Chang & Setser, JCP 69, 3885 (1978), Ar(4p) relaxation.
10. Sauerbrey et al., J. Phys. B 19, 2279 (1986), Ar2* absorption.
11. Any measured ArF(B) radiative lifetime, and any measured Ar
    metastable-metastable ionization rate.
