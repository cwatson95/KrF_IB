# Inverse-bremsstrahlung heating of KrF laser gases

Notes by Cooper K. Watson, October 2026. Each folder holds one note: the PDF, and the LaTeX
source with the tables and figures it reads, so the PDF can be rebuilt with
`latexmk -pdf <name>.tex` inside that folder. `compendium/` holds all six in one document.

## All in one document

`compendium/compendium.pdf` (163 pages) collects the six notes below as chapters 4 to 9. Its
Part I (chapters 1 to 3) is new: it asks whether inverse bremsstrahlung can reheat the electrons
of the afterglow, after the electron beam has stopped, and make the laser more efficient, and
answers link by link (coupling, absorption, heating, partition, laser and wall plug), each
verified, invalidated or open on the evidence of the notes, with what would settle the open ones.
Rebuild with `latexmk -pdf compendium.tex` inside `compendium/`.

`compendium/compendium_summary.pdf` (5 pages; the same text in `compendium_summary.md`) is the
short version: the quantum derivation of the absorption in seven steps, the bottom line, and what a
time-dependent kinetic model of the Kr/F2 laser gas gives for the efficiency when the afterglow is
reheated. Rebuild with `latexmk -pdf compendium_summary.tex` inside `compendium/`.

## The question

At what frequency does an electromagnetic wave best heat the free electrons of a KrF laser gas by
inverse bremsstrahlung, at electron densities of 10^14 to 10^16 cm^-3, 0.1 to 2 atm and 100 to
500 kW/cm^3 deposited for a microsecond, and what can that heating do for the laser?

## The notes, in reading order

| Folder | Note | Pages | What it answers |
|---|---|---|---|
| `envelope_kr/` | `envelope_kr_v6.pdf` | 17 | The best frequency for pure Kr and Kr + 2 % F2, derived step by step: absorption, electron temperature, GHz against THz, uncertainty, and a worked example |
| `quant_ib_v4/` | `quant_ib_v4.pdf` | 31 | The same question with the electron kinetics and the excited atoms followed: Ne- and Ar-buffered fills and Kr + 0.5 % F2, the measured Kr 5s to 5p cross section, and a 1 µs pulse |
| `nonlinear_ib/` | `nonlinear_ib_v2.pdf` | 20 | Whether the nonlinearity of inverse bremsstrahlung, or pumping by an auxiliary field, can make gas lasers more efficient |
| `ib_ladder_v3/` | `ib_ladder_v3.pdf` | 18 | The absorption model built up one physical effect at a time and tested against the measured absorption of Turnbull et al.; electron-atom absorption in argon |
| `derivations/` | `detailed_derivation.pdf` | 40 | The quantum theory of inverse bremsstrahlung on ions, every step written out |
| `derivations/electron_atom/` | `detailed_derivation_v2.pdf` | 23 | Free-free absorption by electron-atom collisions at THz frequencies, and how accurate the classical absorption kernel is there |

## Main results

- **Afterglow reheating** (`compendium`, Part I): the heating half stands and the efficiency half
  is open. At the afterglow's 10^20 m^-3 a wave of 0.64 to 3.8 THz is absorbed through 10 cm of gas,
  where microwaves heat a layer millimetres deep, and 1.4×10^3 to 1.5×10^4 W/cm^2 hold its
  electrons. A THz wave would form 1.4 to 1.8 times the e-beam's KrF* per absorbed joule, but that
  count leaves out the electron collisions of Kr* and KrF*, and no note computes the extraction or
  the wall plug. The nonlinearity of inverse bremsstrahlung and the laser's own field are ruled out
  as levers.
- **The efficiency it buys** (`compendium_summary`): reheating is possible, since quantum mechanics
  gives the afterglow's electrons the classical absorption to 1.6 % at 0.3 to 4 THz, but it does not
  make the laser meaningfully more efficient. In a time-dependent kinetic model of the Kr/F2 gas,
  the compendium's 3-MW/cm^2 wave in the afterglow lowers the 248-nm efficiency by 0.037 points with
  every rate at its default, and raises it by at most 0.14 points, about 2 % of the efficiency, only
  if the unmeasured Kr2* + F → KrF(B) + Kr runs faster than about 2×10^-12 cm^3/s and the wave's
  source is at least 0.49 times as efficient as the pump's.
- **Best frequency** (`envelope_kr`): pure Kr at 10^15 cm^-3, 0.5 atm and 300 kW/cm^3, in a gas 10 cm
  deep: f_best = 6.2 ± 0.6 THz. At 95 GHz the same gas absorbs within 0.17 mm and reflects 31 %.
- **With the kinetics** (`quant_ib_v4`): Kr + 0.5 % F2 wants 1.05 to 8.8 THz at 0.1 atm and 3.8 to
  26.9 THz at 2 atm.
- **Efficiency** (`nonlinear_ib`): the nonlinearity of inverse bremsstrahlung does not by itself
  improve a gas laser. In the Ne fill at beam-on density, where the field's power exceeds the beam's,
  field pumping forms KrF* at 1.5 to 1.8 times the e-beam's formation efficiency of 0.198.
- **Absorption model** (`ib_ladder_v3`): with no fitted parameter, 0.96 ± 0.04 of the measured
  absorption for ten nitrogen-hydrogen shots and 1.19 ± 0.06 for four hydrogen shots.
- **Atoms at THz** (`derivations/electron_atom`): from 0.3 to 4 THz the classical kernel is within
  1.6 % of the quantum result in either fill, for electron temperatures from 300 K to 1.08 eV.
- **Theory** (`derivations`): the quantum theory fixes both cutoffs of the Coulomb logarithm,
  b_max = v_t/ω and b_min = 0.445 ħ/(m v_t).

## Numbers and code

The numbers in these notes were computed and checked numerically by scripts; most notes list their
checks in an appendix. `derivations/` includes the two scripts of that note:
`detailed_derivation_check.py` (writes Appendix C) and `detailed_derivation_fig.py` (writes the
figures), which need Python 3 with numpy >= 2.0, scipy and matplotlib. The other scripts are
available on request.
