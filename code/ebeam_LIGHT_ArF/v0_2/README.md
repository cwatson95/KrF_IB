# LIGHT_ArF_ebeam v0_2: the 0-D ArF model, with an electron-beam drive

v0_2 here is LIGHT's v0_2 converted from Kr/F2 to Ar/F2 (LIGHT_ArF,
2026-09-30), with an electron-beam drive added (LIGHT_ArF_ebeam,
2026-09-30). The equations, the energy ledger, the two-temperature electron
closure and the band machinery are KrF v0_2's. The gas, the cross sections,
the tail table and every rate are ArF's. One rare gas, argon, and no
buffer. The KrF original stays in `LIGHT/v0_2/`, the ArF one without the
drive in `LIGHT_ArF/v0_2/`.

```bash
python3 LIGHT_ArF_ebeam/run_v0_2.py                        # the reference point
python3 LIGHT_ArF_ebeam/run_v0_2.py --drive beam           # the same under the e-beam
python3 LIGHT_ArF_ebeam/run_v0_2.py --bands                # each band's two ends, about 6 min
python3 LIGHT_ArF_ebeam/run_v0_2.py --rate arf_yield=1.0   # one band at one value
python3 LIGHT_ArF_ebeam/run_v0_2.py --pressure 0.82 --f2-pct 0.3 --power 1090 --pulse 0.12
python3 LIGHT_ArF_ebeam/run_v0_2.py --eedf maxwellian      # no kinetic tail
python3 -m pytest LIGHT_ArF_ebeam/v0_2/tests -q            # 88 tests, about 7 minutes
```

The input is a given power, deposited uniformly, in one of two ways
(`Params.drive`). LIGHT's own input (`'power'`, the default) puts it into
the electron gas directly. No electrons are injected, and n_e follows
ionization, attachment and recombination from its seed. The electron beam
(`'beam'`, `beam.py`) makes ion pairs at W eV each, with their excited
atoms and prompt F2 dissociations, puts the secondaries' share into the
electron gas and the cascade's excess into the gas.

## The electron-beam drive

`beam.py` is REBUILD's beam physics (`SPECK/REBUILD/MODULES/beam.py`) run
on this model's own cross sections: the Bethe/Moller stopping power, and
the Spencer-Fano cascade with its prompt floor at the lowest Ar threshold,
the secondary's energy drawn from the measured Opal-Peterson-Beaty spectrum
(Ebar(Ar) = 10.0 eV) and the yields per pair read as the cascade's
asymptotic slopes. Per ion pair at 0.5 % F2 it gives W 25.2 eV (4.5 % under
the measured 26.4, scored and never set), 0.39 excited atoms, 4.4 eV to
the electron gas and 0.55 eV of prompt heat. Every ionization is one Ar+
(no F2+ here), every Ar excitation one Ar*, a prompt F2 C1 or H1 excitation
2 F, and the identity W = sum(y E_form) + heat + eps_sec closes exactly.
`beam_route='integral', beam_sharing='equal'` is REBUILD's own rule, whose
yields oscillate with the source energy (the Known-bad list).
`beam.deposited_power_W_cm3` gives a printed beam's deposited power in a
thin cell.

## The network

Every rate is a measurement, read first-hand from the page it cites
(`docs/arf_data/SOURCES_CHECKED.md`: many are review tables quoting the
measurement, and F + F + Ar is KrF's reading of Baulch 1981), or a band. A
measured rate known only as a published spread keeps its measurement as the
band's default. A rate with no measurement is a band with default 0.

- **Electron impact on Ar.** All 37 excitation blocks of the IST-Lisbon set
  (Yanguas-Gil, Cotrino & Alves 2005), each its own row with its own
  threshold. The four 4s blocks (G1_*) make Ar*, lumped at 3P2 (11.548 eV).
  The 33 higher blocks (C1_*) cascade promptly to Ar*, the excess heating
  the gas. Ionization (C2), elastic loss (G6). F2's SIGLO attachment,
  dissociation and vibration rows are KrF's.
- **Ar* + F2** (G3, G3X): the total 8.1e-10 cm3/s (band `ar_f2`, the mean
  of Velazco 1976, Chen et al. and Ueda 1984) makes ArF(B) with yield 0.5
  (band `arf_yield`), Setser's measurement. The rest dissociates.
- **ArF(B)** is one state. Its radiative lifetime has never been measured.
  4.2 ns is Dunning & Hay's calculation (`--tau-sp`). Its quenchers are the
  measured Stern-Volmer products k tau: F2 7.55e-18 cm3 (Rokni 1977, Chen
  & Payne 1979), Ar 3.6e-20 cm3 (Rokni), 2Ar -> Ar2F* 2.27e-39 cm6 (band
  `arf_trimer`, Rokni, Chen & Payne and Jancaitis), electrons 5.75e-16 cm3
  (band `arf_e`, Duzy 1980). The rows divide each product by tau_B, so a
  change of the calculated lifetime keeps every measured product.
- **Ar2F*** radiates in 203 ns and is quenched by F2 at 1.93e-10 cm3/s
  (means of Chen & Payne, Jancaitis and Ueda), by F atoms at 0.5-1 times
  that (band `ar2f_f`, Jancaitis) and by Ar below 1e-14 (band `ar2f_ar`,
  default 0).
- **Ions.** Ar+ + 2Ar -> Ar2+ (2.2e-31 cm6/s, the median of Brau's eight
  measurements), Ar2+ + e -> Ar* + Ar (7.5e-7 (Te/300 K)^-0.67), and ion-ion
  recombination to ArF(B) on Hoffman & Moreno's calculation as Lee 1989
  prints it for Ar+ and Ar2+, screened, with Thomson's Tg^-3 on the
  three-body limit (no absolute measurement exists).
- **The Ar2* route.** Ar* + 2Ar -> Ar2* (6.0e-33 cm6/s, Ueda 1984, band
  `ar2_form`), Ar2* radiation at 127 nm (the 3.2-us triplet, band
  `ar2_rad`), Ar2* + F2 (5.2e-10, band `ar2_f2`) making ArF(B) with share
  0.91 (band `ar2_f2_arf`) and Ar2F* otherwise. Ar2* + F and Ar2* + Ar2*
  have no measurement (bands at 0).
- **Energetics.** ArF(B) sits at the spectroscopic lower bound, 0.80 eV +
  the measured 1933-A band peak (7.214 eV above Ar + F2); Ar2F* 0.65 eV
  below it (calculated); Ar2* 0.68 eV below Ar* (calculated); Ar2+ 1.23 eV
  below Ar+ (measured). The ledger closes by construction (4e-14 at the
  reference point).

Not carried, with the reason: the ArF(C) state (no measured position or
mixing; calculated 0.25 eV above B); Ar+ + F2 charge transfer (no
measurement); ArF* + F and Ar2F* + e (estimates only); stepwise ionization
and superelastics from Ar* (the optional `electron` group, estimates, off
by default).

## The electron closure

The basic network integrates its 38 Ar rows over a two-temperature EEDF:
the Maxwellian at Te below 11.548 eV and a tail temperature r Te above it,
r read from `data/ar_tail_ratio.txt` at the run's n_e/N and F2 fraction.
The table is made by `docs/arf_data/ar_tail.py` from ML's kinetic
beam-sustained Boltzmann solution with e-e collisions, exactly as KrF's
table is made (70 states at 1 atm; the pressure check reads back 0.97 to
1.02 at 0.5 and 2 atm). `--eedf maxwellian` runs without it.

## At the reference point

1 atm, 0.5 % F2, 300 kW/cm3 for 1 us, to 2 us (the KrF reference point's
inputs):

| run | 193 nm | 127 nm | formed ArF(B) x 6.41 eV | Tg end |
|---|---|---|---|---|
| defaults | 5.05 % | 4.60 % | 18.7 % | 675 K |
| Maxwellian EEDF | 5.89 % | 3.55 % | 22.3 % | 696 K |
| ArF* yield 1.0 (the published models') | 8.89 % | 5.17 % | 32.9 % | 636 K |
| Ar2* + F at 1e-9 (no measurement) | 8.56 % | 0.03 % | 30.4 % | 773 K |
| trimer Tg^-3.9 (NRL's form) | 7.01 % | 4.82 % | 18.4 % | 661 K |

Te averages 6.8 eV over the pulse with the cold kinetic tail (1.8 eV on the
Maxwellian), and ionization is 1.6 % of the input. ArF(B) comes 83 % from
Ar* + F2, 13 % from Ar2* + F2 and 4 % from ion-ion recombination. Of the
ArF(B) formed, 27 % radiates, 36 % forms Ar2F*, 24 % is quenched by Ar and
11 % by F2. F2 burns out at 728 ns. The seed only sets
the first nanoseconds: 1e10 to 1e15 cm^-3 give 5.05-5.06 % (no seed trap,
unlike KrF).

`--bands` at the reference point (193 nm, percent of the input): the three
levers are Ar2* + F (0 -> 1e-9: 5.05 -> 8.56), the ArF* yield (0.25 -> 1:
3.01 -> 8.89) and the trimer (k tau 1.6e-39 -> 5e-39: 5.63 -> 3.56; Tg^-3.9:
7.01). Ar* + 2Ar moves it 4.90 -> 5.89 (2.8e-33 -> 3.3e-32), and every
other band by 0.3 or less.

## Accuracy

`docs/arf_data/SCORE.md` scores v0_2 and 1Dv0_2 on ArF measurements (Electra,
Keio and Sandia gain, Rokni's and Jancaitis's fluorescence scans, SRI's
absolute yields, Morozov's pure-Ar VUV efficiency, and the decay data the
constants came from). The accuracy overview is
`docs/accuracy_overview.pdf`, and the summary is in `LIGHT_ArF_ebeam/README.md`.

## Tests

`python3 -m pytest LIGHT_ArF_ebeam/v0_2/tests -q`: the network as data (balance,
citations, the IST-Lisbon table against the vendored file, the cited
means), the cross-section integrator against closed forms and the pinned
files, the two-temperature tail, and the model (ledger closures across
modes, the measured Stern-Volmer branching at any tau_B, the harpoon
branching, capture energies, the rate laws against their printed forms, the
bands, the reference-point regression), `test_beam.py`: the drive (its
identity and atom balance on every route, REBUILD parity where
SPECK/REBUILD is in the checkout, W against the measurement, the
equal-sharing oscillation and the Opal-Beaty convergence, the stopping
power, closures under the beam in five modes, and pure argon's one photon
per ion and excitation), and `test_docs.py`: the accuracy
overview's generated rate and score files against fresh ones, every number
it prints defined, its short sources held to the citations, the scores
made on this chemistry, the energy ladder against its generator, and every
node of the system map balanced on a fresh run. One test, the map's numbers
file against a fresh one, takes about 7 minutes and runs only with
LIGHT_DOCS_SLOW=1.
