# LIGHT_Kr v0_2 sources

v0_2 carries v0_1's rows, and their sources are v0_1's (sections 1 to 6,
copied from LIGHT/v0_5/SOURCES.md with their status column, the fluorine
rows left out). Section 7 is new: the level-resolved cross sections v0_2
reads, the Hornbeck-Molnar band, and the swarm literature behind the
accuracy overview's low-field ionization coefficient.

The status column says how each value was read. **checked** means read in
the full text (page given), **abstract** means the abstract only,
**secondhand** means quoted by another source that was read (named),
**derived** means arithmetic on the cited numbers, **calculation** a
calculation, **snippet** a search-engine rendering of an abstract (verify
before quoting), and **estimate** a researcher's own estimate, never used
to set a value.

## 1. Charged-particle recombination (LIGHT v0_1)

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
| its temperature law, theory | "k3 prop. to T^-3/4" for rare-gas ion-atom-atom recombination | Perez-Rios & Greene, PRA 98, 062707 (2018), arXiv 1810.02108 (in `LIGHT/docs/krf_data/sources/`) | checked |
| Kr+ + e + e | 2.0e-27 n_e Te(eV)^-4.5 cm3/s, below 1e-4 of DR at LIGHT's conditions | Mansbach & Keck, PR 181, 275 (1969) | abstract |

**Used.** Charge-conserving n_e with no injection. Kr2+ + e at Shiu &
Biondi (R11p, rate kind te_pow, where a captured electron takes 0.95 Te),
onto the 5p rung since LIGHT v0_4 (section 4). Kr+ + 2Kr at 2.4e-31
cm6/s, the mean of the three measurements, on Perez-Rios & Greene's
T^-3/4 (K15). Not included: Kr+ + e + e (below 1e-4 of DR), Kr2+
redissociation above about 1000 K (no rates).

## 2. The Kr2* route and the bands (LIGHT v0_2)

Jancaitis is his LLNL thesis, *Kinetic Modelling of Krypton Fluoride Laser
Systems*, UCRL-53465 (1983), with printed page numbers.

| Quantity | Value | Source | Status |
|---|---|---|---|
| Kr* + 2Kr -> Kr2* + Kr, 0.5-1 atm | 1.7 +- 0.2e-32 cm6/s at 300 K, e-beam, Kr* absorption | Ueda, Kanada, Sasaki & Takuma, Rev. Laser Eng. 12, 357 (1984) Table V | checked |
| same, low pressure | 1.79e-32 (Salamero 1979), 2.59e-32 (Tracy & Oskam, JCP 65, 1666 (1976)), 4.2e-32 (Smith & Turner 1963), 4.4e-32 (Bouciqué & Mortier 1970) | abstracts | abstract |
| same | 5.5e-32 | Hughes 1977, quoted by Quigley & Hughes, APL 32, 649 (1978) Table I | secondhand |
| same, Jancaitis's model | 2.6e-32, citing Tracy & Oskam | Jancaitis Table II p.86 | checked |
| its temperature law | Ea = -0.020 eV, which is (Tg/300)^-0.64 over 300-451 K | Wieme & Lenaerts, JCP 74, 483 (1981) | abstract |
| Kr2* -> 2Kr + hv, 147 nm | 3.7e6 1/s = 270 ns, the slow (1u, 0u-) component: Morikawa et al., JCP 91, 1469 (1989) 270 +- 5 ns; Salamero 1979 250 +- 20, Quigley & Hughes 1978 280 +- 30, Audouard 1988 265, Oka 1974 353 ns, Leichner & Ericson 1974 about 0.30 us at 1000 Torr | as the row KR2_RAD cites them | abstract |
| Kr2* + Kr2* -> Kr2+ + 2Kr + e | no measurement. Kannari's 5.0e-10 is Lorents & Olson's unpublished 1973 SRI report, which gives Ar2* the same value | Kannari 1985 Table I p.4312, ref. 87 | checked |
| same, Jancaitis's model | 1.0e-10, from Werner et al., IEEE J. Quant. Elect. 13, 769 | Jancaitis Table II p.86 and ref. 1 p.109 | secondhand |
| Kr* + Kr* -> Kr+ + Kr + e | 5.0e-10, from Leichner & Ericson, PRA 9, 251 (1974) | Kannari 1985 Table I p.4312, ref. 86 | secondhand |
| same, Jancaitis's model | 1.5e-10, from Werner's MIT thesis (1975) | Jancaitis Table II p.86 and ref. 21 p.109 | secondhand |
| Kr2+ binding | D0 = 1.15 eV, by photoionization | Ng, Trevor, Mahan & Lee, JCP 66, 446 (1977), via Krauss & Mies in *Excimer Lasers* 1984, Table 2.1 p.9 | secondhand |
| Kr2* binding | none among the sources held. Krauss & Mies' Table 2.1 gives the Rydberg dimer's D_e for He, Ne and Ar only | *Excimer Lasers* 1984, Table 2.1 p.9 | checked |

**Used.** Kr* + 2Kr keeps Ueda's 1.7e-32 on Wieme & Lenaerts' law, Kr2*
radiates at 270 ns, and Kr2* + Kr2* has no measurement, so it is 0 unless
a band run sets it. Kr* + 2Kr and Kr* + Kr* read bands (`network.BANDS`).
Kr2+ keeps 1.15 eV, and Kr2* keeps REBUILD's estimated well, which only
moves energy between heat and storage.

## 3. The kinetic EEDF (LIGHT v0_3)

The scheme and its conventions are the solver ML vendored from the website
LXCat build, which ML checked against LoKI-MC's Monte Carlo (`ML/ml/
boltzmann/__init__.py`: under 0.5 % in mean energy and 1 % in drift
velocity for Kr, Ar, He, Xe and a Kr/F/Ar mixture). No literature value
enters the model through it.

| Quantity | Value | Source | Status |
|---|---|---|---|
| the two-term equation, its temporal growth term and its Scharfetter-Gummel fluxes | the scheme LoKI-B solves | Hagelaar & Pitchford, PSST 14, 722 (2005), as cited by ML's vendored solver (`ML/ml/boltzmann/boltz_twin.py`) | secondhand |
| the electron-electron term | W_ee = -3 alpha A1, D_ee = 2 alpha (A2 + eps^3/2 A3), with ln Lambda at kTe = 2/3 <eps> | ML's `ee.py`, citing Hagelaar & Pitchford sec. 2.5 and re-derived there from the Rosenbluth coefficients | secondhand |
| a Maxwellian leaves the e-e term alone | 3 T A1 = 2 (A2 + eps^3/2 A3) for F = e^(-u/T) | integration by parts | derived, and held on the grid to 1e-10 (`tests/test_kinetic.py`) |
| cross sections | LoKI-MC's krypton file (Biagi-v7.1), as `xsec.py` reads it | `LoKI-MC/Code/Input/Krypton/Kr_LXCat.txt` | md5-pinned copy |
| the elastic and field powers | the parts of the Scharfetter-Gummel flux, W F_b - D c(z) dF/deps with c(z) = (z/2) coth(z/2), z = W du/D | the flux itself, since B(-z) - B(z) = z and B(-z) + B(z) = z coth(z/2) | derived, and held on the grid to 1e-10 (`tests/test_kinetic.py`) |
| the time step in the pulse | backward Euler on sqrt(eps) dF/dt, the term the steady scheme drops, with nu the population's own change | the two-term equation as Hagelaar & Pitchford write it (2005), and option 2 of `LIGHT_Kr/improvement_for_EEDF.txt` | derived. A long step is the steady solution to 1e-10 and a short one its start to 1e-8 (`tests/test_kinetic.py`) |

**Used.** Every electron rate the network books by a cross section, Kr's
elastic loss, the recombination template's rate and the energies the
captured electrons carry come from the kinetic EEDF at LIGHT's mean
energy. The Penning energies, 5.834 and 5.554 eV, are the network's own.

## 4. The 5p rung, the walls and the wave (LIGHT v0_4)

| Row / constant | Value | Source | Status |
|---|---|---|---|
| Kr** energy | 11.30 eV, the Biagi 'Kr(p)' threshold = 5p[1/2]_1 (NIST ASD 11.304) | Biagi-v7.1 file; NIST ASD | checked |
| PRAD, 5p -> 5s lifetime | 29.0 ns, the 2J+1 mean of the ten measured level lifetimes 23.5-40.9 ns | Dzierzega, Volz, Nave & Griesmann, PRA 62, 022505 (2000) Table 1, collecting Schmitt et al. (beam-gas-laser), Chang, Horiguchi & Setser, JCP 73, 778 (1980), Whitehead et al., JCP 102, 1965 (1995), Cannon et al., PRA 47, 147 (1993), Ernst & Schulz-Gulde, Physica 93C, 136 (1978) | checked (arXiv physics/0005057) |
| KPS, Kr(5p) + Kr | (1.7 +- 0.2)e-10 cm3/s, the total loss of 5p'[3/2]_2 in Kr at 300 K | Gazeli et al., Phys. Plasmas 28, 043301 (2021), Stern-Volmer 0.1-10 mbar | checked (arXiv 2012.14685); product channel unresolved. Band to 3.0e-11 (Jancaitis UCRL-53465 Table II, Ar analogue from Chang & Setser 1978) via Kannari's 1.0e-10 'estimated' |
| K76, e + Kr(5s) -> Kr(5p), the fit | 9.0e-7 Te^0.5 exp(-1.51/Te) cm3/s | Kannari 1985 Table I [76] = Jacob & Mangano, APL 28, 724 (1976), a fit to e-beam KrF data | checked (Kannari's page); a fit. Params.k76 = 'kannari' since LIGHT v0_5 (section 5) |
| K76r | detailed balance, g(5s)/g(5p) = 12/36 | statistical weights of the manifolds | derived |
| KSTEP_S, e + Kr(5s) -> Kr+ | band 0 to 1.5e-7 on Te^0.68 exp(-4.04/Te) | Kannari 1985 Table I [79] = Vriens, Physica 31 (1965), a calculation (7.0e-8); REBUILD's template 1.56e-7 Te^0.71 exp(-4.09/Te) (1.5e-7 on this form at 1.3 eV) | calculation; no measurement of metastable Kr ionization found (Ton-That & Flannery 1977, Hyman 1979, Margreiter 1990 are calculations) |
| KSTEP_P, e + Kr(5p) -> Kr+ | band 0 to 2.0e-7 on Te^0.62 exp(-2.53/Te) | Kannari 1985 Table I [79], Vriens | calculation |
| R11p products | Kr2+ + e -> Kr(5p) + Kr | Shiu & Biondi, PRA 16, 1817 (1977): the 5p -> 5s lines in the afterglow | abstract |
| ION_K0 | Kr+ 0.84, Kr2+ 1.22 cm2/V/s at the Loschmidt density | Ellis, Pai, McDaniel, Mason & Viehland, ADNDT 17, 177 (1976) | as Fluid1D's ion_mobility.toml reads them (+-3 to 5 %) |
| Turnbull IB model | Sommerfeld lnL (Pradler & Semmelrock's Gaunt fit), ion screening at a_i, the Langdon factor; model error 10 % | Turnbull et al., PRL 130, 145103 (2023) p.145103-4; Phys. Plasmas 31, 063304 (2024) Sec. I and IV | checked (local PDFs): 'all predicted to better than 10 %' |
| Electron-neutral IB kernel | the Lorentz gas on the Biagi momentum transfer | LIGHT v0's kernel; IB/electron_atom finds it within 1.6 % of the quantum rate at 0.3-4 THz | checked |
| HF field term | D = (gamma/3)(E/N)^2 eps sigma~/(sigma~^2 + (omega/N)^2/(gamma^2 eps)) | Hagelaar & Pitchford, PSST 14, 722 (2005) | the form IB/hf_boltzmann.py carries |
| Wall loss | D_a/Lambda^2, D_a = mu_i k(Te + Tg)/e, Lambda = gap/pi; 2 Te per lost electron to the wall, Te to the gas | the electropositive global-model closure (Lieberman & Lichtenberg ch. 5), stated, not fitted | closure |

## 5. The 5s -> 5p excitation measured (LIGHT v0_5)

| Row / constant | Value | Source | Status |
|---|---|---|---|
| K76, e + Kr(5s) -> e + Kr(5p) | the measured cross sections out of the two metastables into the ten 4p5 5p levels, Tables II (1s5, J = 2, +-35 %; 2p3 and 2p4 +-50 %) and III (1s3, J = 0, +-40 %); the lumped Kr* -> Kr** cross section is the statistical-weight mean over the 5s manifold with the J = 1 levels (1s4, 1s2, nothing measured) given their same-core metastable's set (`jung2006.py`). Sum at 6 eV: 1s5 106.9e-16, 1s3 127.1e-16 cm2 (LIGHT v0_5: 127.2, section 6) | Jung, Stone, Boffard, Anderson & Lin, PRA 73, 022722 (2006); the Letter PRL 94, 163202 (2005) Table I | checked (both PDFs, in `Kinetic_Modeling/new_articles/`). Apparent cross sections (cascades under 10 %, Sec. II F) used as printed. Replaces Kannari Table I [76], which stays as Params.k76 = 'kannari' |
| band k_5s5p | 1.0 [0.5, 1.4], a scale on the measured cross section | 1.4 the measurement's upper error; 0.5 the metastables alone (the J = 1 levels contributing nothing) | the one assumption (the J = 1 proxies) bracketed |
| K76r | the Klein-Rosseland reverse of the same cross section, g 12/36 over the ledger's gap (1.385 eV) | detailed balance | derived; the table equals g exp(gap/Te) times the forward to 2e-3 (`tests/test_jung2006.py`) |
| Kr 5p'[3/2]_2 natural lifetime | 34.1 ns (the TALIF calibration level; 2p2 in Paschen) | Gazeli, Lombardi, Aubert, Duluard, Prasanna & Hassouni, Plasma 4, 145 (2021) Table 3, quoting Niemi et al. | checked; a cross-check of PRAD's 2J+1 mean 29.0 ns (Dzierzega 2000 gives 2p2 28.6 ns) |

## 6. Two corrections of jung2006's tables (LIGHT_Kr v0_1, 2026-10-13)

Both were read on the paper's page images and in two text extractions
(pdftotext and pypdf). Both are in LIGHT/v0_5/jung2006.py as well.

| Entry | LIGHT v0_5 | The paper | Source | Effect in v0_1 |
|---|---|---|---|---|
| 1s5 -> 2p6 at 4 eV | 21 x 1e-16 cm2 | 23 x 1e-16 cm2 | Jung et al. 2006, Table II, p.022722-6 | lumped cross section +1.15 % at 4 eV; K76 rate +0.3 % at 1-3 eV |
| 1s3 -> 2p5 | 0.4 at 3.5 eV and 0.11 at 6 eV, booked as a core-changing excitation from Table IV | not measured. Table IV's last row is 1s5 -> 2p5 (J = 2 -> J = 0, dipole-forbidden), whose 0.11 at 6 eV is Table II's 2p5 entry; the Letter's Table I prints the same row. Of the core-changing excitations out of 1s3 the paper tried 1s3 -> 2p7 alone: it "could be as large as 5 x 10^-16 cm2 but could also be much smaller" | Table IV, p.022722-10; text, p.022722-7; PRL 94, 163202 (2005) Table I | lumped cross section -0.1 % at 2-4 eV; K76 rate -0.1 % |

With both, the reference point's 147-nm share under LIGHT v0_2's closure
moves from 15.066 % to 15.051 %, and the band k_5s5p (0.5 to 1.4) spans
20.1 % to 13.0 %.

## 7. The level-resolved cross sections, the Hornbeck-Molnar band and the swarm (LIGHT_Kr v0_2, 2026-10-14)

| Quantity | Value | Source | Status |
|---|---|---|---|
| e + Kr cross sections | 51 excitation levels, each its own block at its own threshold (1s5 9.9152 eV to 'HIGH' 13.6 eV), ionization 13.9996 eV, elastic momentum transfer; "Transcribed from S.F. Biagi's FORTRAN code Magboltz, version 8.97, Sept 2011. Data are based in part on results in BSR database" | LXCat 'Biagi' database, retrieved 2026-06-18 by the owner (`old_files/Kr_LXCat_new.txt`, Kinetic_Modeling commit 36f70d3), vendored as `data/Kr_Biagi897_LXCat.txt` (md5 faeab9cb...) | the file's own header |
| v0_1's set | Biagi-v7.1: "The excitation cross sections are grouped into a few effective excitation levels. For more detailed cross sections with multiple levels, see Biagi 8.9" | LXCat 'Biagi-v7.1' header (`v0_1/data/Kr_Biagi71_LXCat.txt`) | checked |
| Magboltz | S. F. Biagi, Nucl. Instrum. Methods A 421, 234 (1999), doi:10.1016/S0168-9002(98)01233-9 | Crossref | metadata checked |
| E(Kr2+), the associative-ionization threshold | IP(Kr) - D0(Kr2+) = 13.9996 - 1.15 = 12.85 eV; the 20 levels of the file from 12.8573 eV up | Ng, Trevor, Mahan & Lee, JCP 66, 446 (1977) (D0; species.py) | as in v0_1 |
| Kr2+ from excited atoms | "The appearance potentials of the molecular ions were less than those of the atomic ions by ... 0.7-0.3+0.7 volt in Kr"; formation by "an excitation by electron impact ... followed by the collision process, He*+He -> He2++e" | Hornbeck & Molnar, Phys. Rev. 84, 621 (1951), doi:10.1103/PhysRev.84.621 | abstract (APS page) |
| the share of a level's excitations that ionizes associatively, against radiation and quenching | not measured for krypton: the band hm_kr, 0 by default, 1 its top (a probability's bound) | searched 2026-10-14: no rate for Kr(nl) + Kr -> Kr2+ + e found | none |
| calculated alpha/N below the low-field measurements | "Figures 13 and 14 show a general trend of calculated ionization coefficients at low E/N being below the experimental values"; "we ascribe the general behaviour to the experimental measurements having a contribution from Penning ionization of impurities which explains the large effect seen at low values of E/N"; "The two-term approximation yields accurate results for alpha_i/N in Kr and Xe, when compared to Monte Carlo simulations"; for Kr below 100 Td Biagi-v7.1 lies slightly below Biagi-v8.9 | Bordage, Biagi, Alves, Bartschat, Chowdhury, Pitchford, Hagelaar, Morgan, Puech & Zatsarinny, J. Phys. D 46, 334003 (2013), doi:10.1088/0022-3727/46/33/334003 | full text through a web extraction (IOP HTML), quotes as extracted |
| the same, earlier | "Penning effects influence the ionization coefficient in the low E/N range and more experiments are needed here, too" | Bordage, Chowdhury, Hagelaar, Pitchford, Biagi, Puech & Morgan, GEC 2011 poster (LXCat notes, lxcat_team4) | read (poster image); not DOI-published, never cited in the documents |
| alpha/N in Kr, 40-250 Td | their results "differ from Kruithof's values, mainly at the lower end of the krypton range" (paraphrase of the abstract) | Jacques, Bruynooghe, Boucique & Wieme, J. Phys. D 19, 1731 (1986), doi:10.1088/0022-3727/19/9/017 | abstract (IOP page) |
| LoKI-MC | the Monte Carlo electron-kinetics solver the overview checks the two-term solution with | Dias, Tejero-del-Caz, Alves & Guerra, Comput. Phys. Commun. 282, 108554 (2023), doi:10.1016/j.cpc.2022.108554; the vendored binary `LoKI-MC/Code/lokimc` | metadata checked |

**Used.** Every level of the file is a row (network.py), named by its
group and its Paschen label: G1 the four 5s levels into Kr*, C1a the ten
5p levels into Kr**, C1b the 4d and 6s levels below 12.5 eV and C1c every
level above, cascading to Kr**; each pays its own threshold and the
excess over its product's energy heats the gas. The species sit at the
file's thresholds, which are NIST's levels (species.py). The twenty
levels at or above E(Kr2+) have a Hornbeck-Molnar twin (HM_*) on the band
hm_kr, 0 by default, which no measurement sets.

**Found.** The accuracy overview's -54 % below 50 Td (v0_1) is not the
solver's: LoKI-MC on the same cross sections agrees with the two-term
solution to 3-10 % (v7.1) and 2.5-4 % (8.97) at 20-50 Td, where its
counting error is 1-19 %. Resolving the levels moves it to -47 %. The
rest is what Bordage et al. 2013 call a general trend of every cross-
section set; the band hm_kr at its top brackets 18 of the 21 points
below 50 Td (not Kruithof's three at and below 17 Td) but puts the model
39 % above the measurements at 100-300 Td, where direct ionization is
within 7.6 % rms (`LIGHT_Kr/docs/accuracy_overview_v0_2.pdf`).

