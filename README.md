# KrF_IB

Modelling a KrF (Kr/F2) laser gas pumped by electrons, and an inverse-
bremsstrahlung (IB) scheme that would reheat its afterglow with an EM
wave. Two 0-D kinetic models carry the work: **LIGHT**, which models Kr/F2,
and **LIGHT_Kr**, which models pure krypton.

**[`overview.pdf`](overview.pdf)** covers LIGHT v0_5 and LIGHT_Kr v0_5:
- how accurate they are against measurements;
- their propagated uncertainty and reference points;
- their energy ladders;
- every reaction rate, band and species energy (Appendices A and B).

## Layout

| path | what |
|---|---|
| `code/LIGHT/` | LIGHT, the Kr/F2 0-D model: `v0` to `v0_6` (the overview covers **`v0_5`**), 1-D variants `1Dv0*`, `run_v0_*.py` launchers, `CHANGELOG.md` (read its Known-bad list before editing), `docs/` (each version's generated accuracy overview and system map) |
| `code/Kr_LIGHT/` | LIGHT_Kr, the pure-krypton copy of LIGHT: `v0_1` to **`v0_5`**, `run_v0_*.py`, `CHANGELOG.md`, `docs/` |
| `code/ML/` | the fits LIGHT's KrF(B,C) constants and LIGHT_Kr's swarm-trained cross sections come from |
| `code/SPECK/` | SPECK/REBUILD, the source of LIGHT's copied rate values (a drift test holds them equal) |
| `code/nate_code/` | `nathaniel/`'s code (not covered by the overview) |
| `IB/` | the envelope notes (`envelope_kr_v6.tex/.pdf`) and the detailed derivation of the IB absorption |
| `doc/` | papers, books, notes and user guides |
| `overview/` | the source of `overview.pdf`: `make_overview.py` (generates the rate, band and species tables from the code), `overview.tex`, `generated/` |

## The models in the overview

| | LIGHT v0_5 | LIGHT_Kr v0_5 |
|---|---|---|
| gas | Kr/F2 | pure Kr |
| network | 59 rows, 12 species, 17 bands | 219 rows, 10 species, 22 bands |
| electrons | kinetic EEDF (two-term Boltzmann, e-e collisions) by default; two-temperature closure as the fast option | the same |
| reference point (1 atm, 300 kW/cm³ for 1 µs) | 248 nm: **5.73 %** of the input (kinetic), 5.7 % (two-temp.) | 147 nm: **6.08 %** (kinetic), 5.58 % (two-temp.) |
| **propagated uncertainty** | **±40 %** (Tahoma net gain; ±38 % on the 248-nm share) | **±115 %** on the 147-nm share |
| what the version adds | K76 (e + Kr* → e + Kr**) on Jung et al.'s measured cross sections; Kleikamp et al.'s laser held out | the electrons' collisions with excited Kr and Kr2* (Vriens–Smeets calculation); each 5s→5p pair at its own threshold |

### LIGHT v0_5: accuracy

From `code/LIGHT/docs/accuracy_overview_v0_5.pdf`.

| test | measured | LIGHT | error |
|---|---|---|---|
| KrF(B,C) decays its rate fit used (Eden 1978, Gadomski 1992) | four sets | | within 6.1 % rms |
| KrF(B,C) decays no fit used (Xu 1993) | three figures | | within 11.7 % rms |
| C/B emission ratio, 300 K (Yu & Setser 1990) | 0.16 ± 0.02 | | −0.5 % |
| Tahoma small-signal gain g0 | 6.39 %/cm | 5.99 | −6 % |
| Tahoma absorption α0 | 1.11 %/cm | 0.83 | −25 % |
| Tahoma net gain g0 − α0 | 5.28 %/cm | 5.16 | −2 % |
| Tahoma saturation intensity | 1.4 ± 0.1 MW/cm² | 1.45 | +4 % |
| Tahoma formation efficiency | 0.27 | 0.261 | −3 % |
| Kleikamp 1986 Kr pressure optimum | 2 atm | 1.26 atm | −37 % |
| Kleikamp 1986 F2 pressure optimum | 10 mbar | > 30 mbar | not reproduced |
| C–A decays (Hakuta 1990, in conflict with the fitted sets) | two figures | | 19.3 and 34.1 % rms |
| KrF(B) lifetime (Burnham & Searles 1977, in conflict) | 9.0 ± 0.6 ns | | −22.3 % |

- **Uncertainty: ±40 %.** Propagated from 22 inputs, it is **±40 %** on
  Tahoma's net gain, ±38 % on the 248-nm share (±2.17 points), and ±33
  and ±37 % on Tahoma's gain and absorption. Its largest term is the
  trimer's gas-temperature law, a calculation measured only to 425 K.

### LIGHT_Kr v0_5: accuracy

From `code/Kr_LIGHT/docs/accuracy_overview_v0_5.pdf`. Its run cache
matches the current code (md5 `eac04189da72`).

| test | measured | error |
|---|---|---|
| drift velocity, 3 sets that agree (train) | | 2.8 % rms |
| ionization coefficient 50–300 Td, 5 sets (train) | | 9.9 % rms |
| characteristic energy D_T/µ, Koizumi 1986 (validation) | | 7.1 % rms |
| ionization coefficient below 50 Td (held out by rule) | | mean −36 % |
| electroluminescence yield 6–13.5 Td, Mano 2022 | 0.113 E/N − 0.3 | 0.99 to 1.04 of it |
| superelastic 1s5, 1s4, 1s3 returns, Blagoev 1984 | 1.5, 5.2, 9.8e-10 cm³/s | ×5.3–6.6, ×5.5–9.5, ×0.5–1.2 |
| energy per ion pair / excitations per pair, Saito 2003 | 23.9 eV / 0.55 | +10 % / +24 % |
| VUV efficiency at low pump, Morozov 2008 | 0.42 ± 0.05 | +29 % |
| VUV efficiency at 3e5 W/cm³, 1 atm, Eckstrom 1988 | near 0.46 | −3 % |
| energy ledger closure | | 4e-15 of the input |

- **Uncertainty: ±115 %.** The 147-nm share's propagated uncertainty is
  **±115 %** (±7.0 points). Three rates dominate it: the electrons' ionization,
  quenching and state mixing of Kr2*. Nothing measures any of them.
- **Unscored regime.** No pure-krypton measurement tests the
  field-driven regime the model is built for.

## Energy ladders and reaction rates

`overview.pdf` holds both. Sections 1.5 and 2.4 give each model's energy
ladder figure and its species' formation energies. Appendices A and B list
every row of each default network in the code's order, with:
- its rate law;
- its rate constant at 300 K and 1 atm, with electron rows on a Maxwellian
  at Te = 1 and 3 eV;
- the bands it reads;
- the start of its citation.

They are followed by the band tables (default and range of every rate
carried with an error band). The rate constants are evaluated by each
model's own `Model._k`, not typed by hand. The models' default closure is
the kinetic EEDF, so the Maxwellian electron rates are a reference.

## Running

```bash
cd code
python3 LIGHT/run_v0_5.py --eedf two_t        # LIGHT reference point, seconds
python3 LIGHT/run_v0_5.py                     # kinetic EEDF, ~20 min
python3 Kr_LIGHT/run_v0_5.py --eedf two_t     # LIGHT_Kr reference point, seconds
/usr/bin/python3 -m pytest LIGHT/v0_5/tests -q --ignore=LIGHT/v0_5/tests/test_kinetic.py
```

Each version's own README (`code/LIGHT/v0_5/README.md`,
`code/Kr_LIGHT/v0_5/README.md`) has the full options.

## Rebuilding overview.pdf

```bash
/usr/bin/python3 overview/make_overview.py    # tables from the code, seconds
cd overview && latexmk -pdf -outdir=build overview.tex && cp build/overview.pdf ../overview.pdf
```

Rerun the generator after any change to either network. The accuracy
numbers in `overview.tex` are copied from the models' generated overviews.
If a model's overview is regenerated (`code/*/docs/make_overview_v0_*.py`),
update them by hand. The two ladder figures are embedded from those
overview PDFs (LIGHT page 9, LIGHT_Kr page 17), so a new overview PDF can
move them.
