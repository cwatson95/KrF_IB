# LIGHT_Kr v0_5 sources

v0_5 carries v0_4's rows that it did not rewrite, and their sources are
v0_4's (sections 1 to 9: v0_1's, copied from LIGHT/v0_5/SOURCES.md with
their status column, the fluorine rows left out; v0_2's level-resolved
cross sections, Hornbeck-Molnar band and swarm literature; v0_3's swarm
data, Monte Carlo and training; v0_4's 5s levels, excimer formation and
redissociations). Section 10 is new: the calculations v0_5 takes where
nothing is measured (Vriens & Smeets; Flannery et al.), the electrons'
collisions with the excited atoms and with Kr2*, the 5p rung's levels and
lines, the measured 5s mixing, the superelastic score, and the corrected
citations. Where a v0_4 row is replaced, section 10 says so; the older
sections' "Used" paragraphs describe v0_4. DOIs in section 10 were
checked against Crossref on 2026-10-16.

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
| its temperature law | Ea = -0.020 eV, which is (Tg/300)^-0.64, measured "between 200 and 400 K" (title and abstract; v0_1 to v0_3 said 300-451 K, the range of the power-law conversion) | Wieme & Lenaerts, JCP 74, 483 (1981), doi:10.1063/1.440855 | abstract |
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


## 8. The swarm data, the Monte Carlo and the training (LIGHT_Kr v0_3, 2026-10-15)

The data live in ML's registry (`ML/ml/kr_data.py`, `ML/data/swarm/`).
LXCat's LAPLACE database ("measurements after 1975", digitized by S.
Chowdhury, LAPLACE, from the papers' tables or figures) was read with the
owner's approval of LXCat's terms of use on 2026-10-14. Its text output
is stored unedited and md5-pinned. Dutton's tables were transcribed from
the printed page (ML, 2026-09-27). None of the data files prints an
error; the class values are ML's (3 % W after 1985, 5 % before, 7 %
D_T/mu, 10 % alpha/N, and 2 % in E/N).

| Dataset | What | Source | Status |
|---|---|---|---|
| Hunter, Carter & Christophorou 1988 | W, 0.002-3 Td, 301 K (train) | Phys. Rev. A 38, 5539, doi:10.1103/PhysRevA.38.5539 | LAPLACE, from the paper's table; Crossref metadata checked |
| Nakamura 1990 | W, 0.3-53 Td (train) | in Nonequilibrium Effects in Ion and Electron Transport (Springer 1990) p. 363, doi:10.1007/978-1-4613-0661-0_33 | LAPLACE, digitized from a figure; metadata checked |
| Pack, Voshall & Phelps 1962 | W, 0.0008-3 Td at 195, 300, 368 K (368 K trains; 195 K validation) | Phys. Rev. 127, 2084, doi:10.1103/PhysRev.127.2084 | Dutton Table 1.9, printed |
| Pack, Voshall, Phelps & Kline 1992 | W at 195 and 300 K (held out; 195 K validation) | J. Appl. Phys. 71, 5363, doi:10.1063/1.350555 | LAPLACE, from the paper's table; metadata checked |
| Dall'Armi, Brown, Purdie & Fletcher 1992 | W and alpha/N, 10-80 Td (alpha/N from 50 Td trains) | Aust. J. Phys. 45, 185, doi:10.1071/PH920185 | LAPLACE, digitized from figures; abstract (search snippet); metadata checked |
| Bowe 1960 | W, 0.025-7.4 Td (held out) | Phys. Rev. 117, 1411, doi:10.1103/PhysRev.117.1411 | LAPLACE, digitized; metadata checked |
| English & Hanna 1953 | W, 0.04-2 Td (held out) | Can. J. Phys. 31, 768, doi:10.1139/p53-070 | LAPLACE, digitized; metadata checked |
| Kucukarpaci & Lucas 1981 | W, 14-820 Td (held out; above 300 Td check) | J. Phys. D 14, 2001, doi:10.1088/0022-3727/14/11/008 | LAPLACE, digitized; metadata checked |
| Koizumi, Shirakawa & Ogawa 1986 | D_T/mu, 0.005-0.9 Td (held out; validation) | J. Phys. B 19, 2331, doi:10.1088/0022-3700/19/15/014 | LAPLACE, digitized; metadata checked |
| Al-Amin & Lucas 1987 | D_T/mu, 1.4-5650 Td (check) | J. Phys. D 20, 1590, doi:10.1088/0022-3727/20/12/007 | LAPLACE, digitized; metadata checked |
| Jacques, Bruynooghe, Boucique & Wieme 1986 | alpha/N, 40-260 Td (60-260 trains) | J. Phys. D 19, 1731, doi:10.1088/0022-3727/19/9/017 | LAPLACE, from the paper's table |
| Bhattacharya 1979 | alpha/N, 62-625 Td (to 300 trains) | J. Appl. Phys. 50, 6207, doi:10.1063/1.325754 (the pure-Kr end of a Ne-Kr Penning-mixture study) | LAPLACE, digitized; metadata checked |
| Specht, Lawton & DeTemple 1980 | alpha/N, direct (18-26 Td, held out) and total (4.5-24 Td, check) | J. Appl. Phys. 51, 166, doi:10.1063/1.327395. Abstract: "At low field values, the ionization coefficient was found to be anomalously large, a fact attributed to surface photoelectron emission from radiating metastables. This contribution also explains the discrepancy between earlier measurements and recent calculations based on the transport equation." | LAPLACE, digitized; abstract (Crossref) |
| Kruithof 1940, Heylen 1971 | alpha/N (50-300 Td trains; below 50 held out by rule) | as v0_2 (section 7) | Dutton Table 7.4, printed |
| LXCat | the platform | Pancheshnyi et al., Chem. Phys. 398, 148 (2012); Carbone et al., Atoms 9, 16 (2021) | as v0_2 |

| Quantity | Value | Source | Status |
|---|---|---|---|
| what a swarm experiment measures | the bulk W (a drift tube) and the SST coefficient (a Townsend experiment) against the two-term flux W and k_i/W; the Monte Carlo's D_T | LoKI-MC (Dias et al. 2023, section 7) run on Biagi's file (xsec.UNIT) and on the trained set, 300 K, 1 Torr, 1e5 electrons, equal sharing, 2026-10-14 (`swarm_mc.py`) | calculation |
| MCMC | the affine-invariant ensemble sampler (stretch move) | Goodman & Weare, Commun. Appl. Math. Comput. Sci. 5, 65 (2010), doi:10.2140/camcos.2010.5.65 | metadata checked |
| the emulator's designs | Latin hypercube sampling | McKay, Beckman & Conover, Technometrics 21, 239 (1979), doi:10.1080/00401706.1979.10489755 | metadata checked |
| the networks' optimizer | Adam | Kingma & Ba, arXiv:1412.6980 (2014), doi:10.48550/arXiv.1412.6980 | metadata checked (DataCite) |
| neural-network swarm analysis | cross sections from transport coefficients with deep networks (prior art for the estimator) | Stokes, Cocks, Brunger & White, Plasma Sources Sci. Technol. 29, 055009 (2020), doi:10.1088/1361-6595/ab85b6 | metadata checked |

**Used.** v0_3's cross sections are Biagi's file times the factors ML
trained on the consistent swarm data (`xsec.TRAINED_FILE`, its record
names the train set). Every other dataset only scores them. Seven of the
nine factors train. The momentum transfer at 0.03 and 0.2 eV stays
Biagi's, a choice made on Koizumi's D_T/mu and Pack's two 195-K drift
velocities (the record's `free_variants`). Those three are therefore
validation sets, not independent tests.

## 9. The 5s levels resolved, the excimer and the redissociations (LIGHT_Kr v0_4, 2026-10-15)

| Quantity | Value | Source | Status |
|---|---|---|---|
| the four 5s levels | 1s5 9.9152, 1s4 10.0324, 1s3 10.5624, 1s2 10.6436 eV (Biagi's thresholds), g = 5, 3, 1, 3 | Biagi, Magboltz 8.97 (section 7) | checked |
| Kr*(1s5) + 2Kr -> Kr2* + Kr | k32 = 3.81 +- 0.08e-32 cm6/s, pure Kr, selective excitation, TALIF (Table 2; the column header misprints k31; no temperature or pressure stated) | Marchal, Sewraj, Gardou, Merbahi & Yousfi, in *Photon Counting - Fundamentals and Applications* (IntechOpen 2018), doi:10.5772/intechopen.71942 | checked (chapter read online) |
| Kr*(1s4) + 2Kr -> Kr2* + Kr | k31 = 9.24 +- 0.80e-32 cm6/s | Marchal et al. 2018, Table 2 | checked |
| same, other values | 4e-32 for 3P2 (6.9e-32 at 196 K) | Turner, Phys. Rev. 158, 121 (1967), doi:10.1103/PhysRev.158.121 | abstract (OpenAlex copy; APS bot check not bypassed) |
| same | dimer formation times from 3P1 and 3P2 in the ratio 3/1; no pressure fall-off 0.1-15 atm; the fast 0u+ component appears above 1.5e19 cm^-3, F/S = 0.21 above 6e19 | Audouard, Laporte, Subtil & Damany, JCP 89, 6176 (1988), doi:10.1063/1.455434 | abstract |
| same, calculated | 3P1 loss nu2 = 1.7e5 + 2050p + 82.0p^2 s^-1, 3P2 nu1 = 35.3/(Lambda^2 p) + 15p + 23.2p^2 (p in mbar), with diffusion and resonance radiation: 82.0 p^2 is 1.36e-31 cm6/s at 295 K, 23.2 p^2 is 3.85e-32 | Janssens, Vanmarcke, Desoppere, Bouciqué & Wieme, JCP 86, 4935 (1987), doi:10.1063/1.452663: a quasiclassical theory compared with their alpha-excitation experiments | abstract |
| Kr*(1s3) + Kr | quenching of 3P0 by Kr measured; the value not in a source in hand (v0_5: the abstract gives Kr + Ar, section 10) | Kholin, L'dov, Ustinovskii & Zayarnyi, J. Russ. Laser Res. 36, 74 (2015), doi:10.1007/s10946-015-9479-y | title only |
| Kr*(1s4) + Kr -> Kr*(1s5) + Kr | 13e-14 cm3/s at 295 K, 22e-14 at 196 K; the reverse 7.6e-16 at 295 K | Turner 1967 | abstract |
| Kr*(1s2) + 2Kr -> Kr*(1s4) + 2Kr | above 4e20 cm^-3 the 1P1 decay runs fully through 3P1, a reservoir effect of about 0.2e-32 cm6/s | Audouard, Laporte, Subtil & Reininger, JCP 95, 3283 (1991), doi:10.1063/1.460834 | abstract |
| Kr2*(0u+) lifetime | 3.4 +- 0.1 ns (1u 270 +- 5) | Morikawa et al., JCP 91, 1469 (1989), doi:10.1063/1.457108 | abstract |
| Kr2* + Kr quenching | 194 +- 84 Torr^-1 s^-1 at 298 K (6.0e-15 cm3/s), from the 984-nm absorption's decay, 250-1100 Torr | Firestone & Selander, JCP 88, 717 (1988), doi:10.1063/1.454150 | abstract |
| same | "no evidence for collisional steps", 200-1400 Torr | Oka, Rama Rao, Redpath & Firestone, JCP 61, 4740 (1974), doi:10.1063/1.1681797 | abstract |
| Kr2* A 3Sigma_u+ well | D_e = 5807.0 +- 8.0 cm^-1 (0.720 eV), from v' = 30-38; B 1Sigma_u+ D_e 5629.1 cm^-1, we 219.50 | LaRocque, Lipson, Herman & Stoicheff, JCP 84, 6627 (1986), doi:10.1063/1.450716 | abstract |
| same | D0 = 0.62 +- 0.02 eV, Delta G_1/2 = 172 cm^-1 | Killeen & Eden, JCP 84, 6048 (1986), doi:10.1063/1.450794 | abstract |
| Kr2+ binding | D0 = 9267.8 +- 2.8 cm^-1 = 1.149 eV | Signorell, Hollenstein & Merkt, JCP 114, 9840 (2001), doi:10.1063/1.1370939 | abstract |
| Kr2+ bond length | r_e about 2.6 A, an estimate (Mittman & Weise 1974) | Huber & Herzberg via the NIST Chemistry WebBook | secondhand |
| Kr2+ + e at Te = Tg | 4.2e-6 cm3/s at 900 K, T^-3/2 to 9.6e-7 at 2700 K, shock tube | Cunningham & Hobson, J. Phys. B 5, 1773 (1972), doi:10.1088/0022-3700/5/9/020 | abstract (IOP page) |
| 5p -> 5s Einstein coefficients | the fifteen 2p -> 1s lines, 740-900 nm, grades AAA to B+; the 1s4 and 1s2 resonance lines' A = 2.98e8 and 3.09e8 /s | Kramida, Ralchenko, Reader & NIST ASD Team, NIST Atomic Spectra Database, version 5.12, doi:10.18434/T4W30F (read 2026-10-15) | checked (the database) |
| superelastic return of the 5s levels by thermal electrons | 3P2 1.5 +- 0.2e-10, 3P1 5.2e-10, 3P0 9.8e-10 cm3/s, afterglow | Blagoev, Mishonov & Popov, J. Phys. B 17, 435 (1984), doi:10.1088/0022-3700/17/3/013 | abstract |
| electroluminescence yield | Y/N = 0.113 E/N - 0.30 (1e-17 photons per electron per atom per cm^-2, E/N in Td), threshold 2.7 Td, 1.1 bar, room temperature, below 14 Td | Mano, Henriques, Amaro & Monteiro, Phys. Lett. B 824, 136782 (2022), doi:10.1016/j.physletb.2021.136782 (CC-BY) | checked |

**Used.** Each 5s level is a species (`species.py`). AK8_1s5 at Marchal's
k32 (band kr2_form, 1.7e-32 to 4.2e-32: Ueda to Smith & Turner), AK8_1s4
at k31 (band kr2_form_1s4 down to a third of k32), AK8_1s3 at k32 as a
proxy (band 0 to k31), all on Wieme & Lenaerts' law. Janssens et al.'s
3P1 term is a calculation and is cited, not a band end. The share
kr2_fast of the 1s4 formation goes to the fast state (0 to 0.17 =
F/(F+S)). MIX_1s4 on Turner's two points, MIX_1s2 at Audouard et al.'s
0.2e-32, both with detailed-balance reverses. The redissociations
(AK8r_*, AK8fr_1s4, K15r) are detailed balance on the measured wells
(`thermo.py`: K_eq = (g_M/g_atoms) lambda^3 (kT/2B) q_vib exp(D0/kT),
sigma = 2), Kr2* at 0.709 eV (LaRocque's D_e less half Killeen & Eden's
Delta G_1/2; band kr2_d0 to 0.62), the ion at 1.15 eV; r_e = 2.6 A for
both, and the ion's vibrational quantum taken as the excimer's (the
Rydberg core). QKR2 is 0 by default (band to 8.6e-15). dr_tg moves R11p
from Shiu & Biondi's law (0) to Cunningham & Hobson's coefficients over
it (1). The resonance lines' escape is 0 in a 0-D box (band to Janssens'
pressure-independent 1.7e5 /s). PRAD branches by NIST's A-values weighted
by 2J + 1 (their mean decay 1/29.5 ns against the measured 29.0 ns); KPS
over the levels by 2J + 1 (band kps_5s). Not carried: the superelastic
return of the 5s levels (Blagoev et al.'s rates are for thermal electrons
in an afterglow; the kinetic EEDF needs cross sections), electron-driven
transfers between the 5s levels, any mixing of the two excimer states.
Mano et al.'s yield scores the swarm (`swarm.el_yield`); nothing trains
on it.

## 10. Calculations where nothing is measured, the electrons on excited atoms and Kr2*, the 5p rung (LIGHT_Kr v0_5, 2026-10-16)

| Quantity | Value | Source | Status |
|---|---|---|---|
| ionization of an excited level by electrons | sigma(E) = 4 pi a0^2 R^2/(E + 3.25 E_n) (5/(3 E_n) - 1/E - 2 E_n/(3 E^2)); its Maxwellian rate 9.56e-6 Te^-1.5 exp(-e)/(e^2.33 + 4.38 e^1.72 + 1.32 e), e = E_n/Te | Vriens & Smeets, PRA 22, 940 (1980), doi:10.1103/PhysRevA.22.940 | calculation; the two formulas recalled and held to each other to 0.2 % at Te 2-5 eV, 0.8 % at 1 eV and 3 % at 0.5 eV (the paper behind APS's bot check) |
| the same rows as Kannari prints them | Kr* + e: 7.0e-8 Te^0.68 exp(-4.04/Te) [79, Vriens]; Kr** + e: 2.0e-7 Te^0.62 exp(-2.53/Te) [79]; Kr2* + e -> Kr2+ + 2e: 9.0e-8 Te^0.70 exp(-3.53/Te) [72, Vriens & Smeets] | Kannari, Obara & Fujioka, JAP 57, 4309 (1985), doi:10.1063/1.334590, Table I p.4311 | secondhand (LIGHT's 2026-09-27 reaction survey read p.4311; the repo's CSV transcription flags [72]'s exponents as uncertain, and [72] lies within 5 % of Vriens & Smeets at 2.6 eV) |
| ionization of Kr2* by electrons | binary-encounter cross sections about 1e-15 cm2, twice He2*'s, 40-50 % above the metastable atom's | Flannery, McCann & Winter, J. Phys. B 14, 3789 (1981), doi:10.1088/0022-3700/14/19/025 | abstract |
| ionization of Kr* (and Ne*, Ar*, Xe*) by electrons | calculated | Ton-That & Flannery, PRA 15, 517 (1977), doi:10.1103/PhysRevA.15.517 | title only (abstract unread) |
| e + Kr2* -> e + 2Kr | 1.0e-7 cm3/s "estimated" | Kannari et al. 1985, Table I, p.4312 | checked |
| e + Kr(3P2) -> e + Kr(3P1) | 4.3e-8 cm3/s in a krypton afterglow, the electron temperature not stated | Blagoev & Pilosov, Phys. Lett. A 83, 161 (1981), doi:10.1016/0375-9601(81)90875-6 | abstract (OSTI ETDEWeb) |
| e + Kr(1s5) -> 1s4, 1s3 -> 1s2 | measured over Te 0.07-1 eV; best agreement with the R-matrix cross sections of Allan et al. 2011 | Qiu, Cheng, Zhu & Pu, J. Phys. D 51, 135202 (2018), doi:10.1088/1361-6463/aab0f7 | abstract (numbers behind IOP's bot check) |
| excitation out of every 5s level into the 5p levels | fully relativistic R-matrix; within about 10 % of Jung et al.'s metastable values for four transitions after an energy shift | Zeng, Wu, Jin, Zhao & Yuan, PRA 72, 042707 (2005), doi:10.1103/PhysRevA.72.042707 | abstract (numbers behind APS's bot check) |
| superelastic return of the 5s levels by thermal electrons | as section 9: 3P2 1.5 +- 0.2e-10, 3P1 5.2e-10, 3P0 9.8e-10 cm3/s; to the ground state | Blagoev, Mishonov & Popov 1984 | abstract |
| Kr(3P0) + Kr + Ar -> Kr2* + Ar | (5.3 +- 0.8)e-33 cm6/s; Kr* + 2Ar -> ArKr* 4.6e-36; Kr* + Ar 4.1e-15 cm3/s; e-beam-excited Ar-Kr mixtures | Kholin, L'dov, Ustinovskii & Zayarnyi 2015 | abstract (Springer) |
| Kr2+ binding (the ledger's since v0_5) | D0 = 9267.8 +- 2.8 cm^-1 = 1.1491 eV | Signorell, Hollenstein & Merkt 2001 | abstract |
| Kr(5p) deactivation by Kr, per level | two-photon LIF, lifetimes, bimolecular quenching in Kr and Xe, the 2p8 <-> 2p9 mixing | Whitehead, Pournasr, Bruce, Cai, Kohel, Layne & Keto, JCP 102, 1965 (1995), doi:10.1063/1.468763 | title and abstract (numbers paywalled) |
| the rest of the 5p -> 5s lines | 2p10 -> 1s4 975.176 nm A = 3.13e6 (B+); 2p10 -> 1s3, 1s2 1.26e5, 7.4e4; 2p8 -> 1s2 8.1e4; 2p7 -> 1s3 7.6e4; 2p6 -> 1s2 3.10e5; 2p5 -> 1s2 1.50e5; 2p4 -> 1s5, 1s4 1.5e4, 5.0e4; 2p3 -> 1s5, 1s4 9.80e5, 5.5e4; 2p2 -> 1s5, 1s4 1.10e5, 7.1e5 /s (grades B+ to C) | NIST Atomic Spectra Database, doi:10.18434/T4W30F (Kr I 550-3000 nm, read 2026-10-16) | checked (the database) |

**Used.** KSTEP_S_<level> and KSTEP_P_<2p> are Vriens & Smeets's cross
section at each level's ionization energy (13.9996 eV less the level's)
times the bands kstep_s and kstep_p (1 by default, 0 to 2; v0_1 to v0_4
held them at 0, the owner's choice of 2026-10-16), and EI_KR2/EI_KR2F the
same at the excimer's (Kr2+ less Kr2*: 3.64 and 3.50 eV), band ei_kr2;
Flannery et al.'s ~1e-15 cm2 is the cross section's own peak. Kannari's
fits are cited as 10-40 % from the formula over Te 1-3 eV. EQ_KR2 (band
eq_kr2, 0 to Kannari's 1e-7) and EM_KR2 (band em_kr2, 0 to 1e-7) have no
published value for any rare-gas excimer. SUPER_<level> is the
Klein-Rosseland reverse of each level's G1 cross section; Blagoev et al.'s
rates score it (the overview). EMIX_1s5 is A exp(-0.117 eV/Te) with A
Blagoev & Pilosov's 4.3e-8 (band emix_5s to 2.3e-7, their value at a Te
of 0.07 eV); the other pairs' rates are not in a source in hand and are
not carried. The 5s -> 5p pumps are Jung et al.'s pairs (section 5) each at
its own threshold on Biagi's level energies, the reverses per atom in the
rung's level; the rung's levels are populated at Tg (Whitehead et al.
measured the 2p8 <-> 2p9 mixing in Kr), the lines are NIST's per level, and
KPS keeps Gazeli et al.'s total with its reverse by detailed balance at Tg.
Zeng et al.'s J = 1 cross sections would replace the proxies (band
k_5s5p_j1) if read; Kholin et al.'s value has argon for its third body, so
the 1s3 formation stays a proxy (band kr2_form_1s3).
