# Changelog: LIGHT_ArF_ebeam (LIGHT_ArF with an electron-beam drive)

Dated entries, newest first. Read the Known-bad list before editing.

LIGHT_ArF_ebeam is the owner's copy of `LIGHT_ArF/` (2026-09-30), made to
add an electron-beam drive to its 0-D model and to score it (entry (a)),
then to its 1-D model (entry (b)). The network, the tail table and the
data are LIGHT_ArF's. LIGHT_ArF's own history is in
`LIGHT_ArF/CHANGELOG.md`, and its Known-bad list is carried below because
the code is the same.

## Known-bad: do not reintroduce

E-beam drive:

- **Never read the beam's yields per pair off one source energy with equal
  sharing.** REBUILD's `spencer_fano_yields` (and LIGHT_ArF's prototype
  drive, which copied its numbers) splits what an ionization leaves equally
  between the two electrons and books the ratio at 55 % of the tables'
  ceiling. On the IST-Lisbon set that makes the yields oscillate with the
  source energy: 0.31 to 1.28 excitations per ion and W 28 to 35 eV between
  100 and 900 eV, some points moving x2 with the grid. `v0_2/beam.py` gives
  the secondary the measured Opal-Peterson-Beaty spectrum and reads the
  asymptotic slopes, which agree over every window to 0.2 % (2026-09-30).
- **LIGHT's own input is not a beam.** It puts the power into the electron
  gas, which in pure Ar makes 9.5 excitations per ion at 153 eV per ion,
  a discharge. A beam makes 0.39 per pair at 25.2 eV. Read a beam
  experiment's score in the beam column (2026-09-30).
- **The pulsed Morozov test is not his experiment below 500 hPa.** It
  counts photons for 20 us after a 450-ns pulse, and at 151 hPa Ar* + 2Ar
  takes 12.5 us, so the model's efficiency falls with pressure there as an
  artifact of the window. His beam was continuous, and run to steady state
  the model is flat at 0.54 from 151 to 1000 hPa (2026-09-30).
- **Bound the 1-D electron energy with the electrons after the beam's
  source.** The stage adds the beam's pairs at each substep and then holds
  the mean energy to its bounds. Bounded with the electrons before the
  source, a cell at the ceiling that the beam fills loses the new
  electrons' energy to the cut. `test_1d.py`'s
  `test_beam_electrons_enter_inside_the_energy_bounds` catches it
  (2026-09-30, entry (b)).
- **Compare the 1-D Rokni record with 0-D runs at its own pressures.**
  `score_1d.json` runs Rokni at 0.243 to 4.003 atm, `score_0d.json` at the
  figure's pressures rounded to 0.24 to 4.0. Side by side the two differ
  by up to 1 % from the pressures alone, where at the same pressure the
  beam's 1-D centre follows v0_2 to 5e-5. A first draft of entry (b)
  quoted that 1 % as the 1-D's agreement (2026-09-30).

ArF-specific:

- **Carry the Stern-Volmer products, not k.** The ArF(B) radiative
  lifetime has never been measured (4.2 ns is Dunning & Hay's calculation;
  Rokni and Chen & Payne reduced with 4 ns, Trainor & Jacob with 2.5 ns).
  ArF* + F2, + Ar, + 2Ar and + e are measured as k tau products, and the
  network divides them by `Params.tau_sp_s` (`("inv_tau", "tau_sp_s", k tau)`
  or the `("per_tau",)` band modifier). Never freeze a k derived with an
  assumed lifetime: a change of tau_B would then change every measured
  product (2026-09-30).
- **The ArF* yield of Ar* + F2 is measured at about 0.5, not 1.** Setser's
  ERDA report COO-2807-11 (1977) pp.5-6; Kolts et al. 4.0e-10 of 7.5e-10.
  Every published ArF model books all of Ar* + F2 to ArF*. That value is
  the top of the `arf_yield` band, never the default (2026-09-30).
- **No ArF(C) state from KrF.** KrF's B/C constants are ML's joint fit to
  KrF data. No ArF C position, mixing rate or equilibrium is measured, and
  the calculated B/C splitting is about 0.25 eV (Krauss & Mies p.15). NRL's
  memo carries the KrF 120 cm^-1 over to ArF. Do not copy it (2026-09-30).
- **Read the secondary sources' ArF rows against the originals.** Brau's
  Table 4.4 lists Ar2* + F2 "3 x 10^-10 (SRI)", which is SRI's Kr2* row;
  SRI prints 2.5e-10 for Ar2* + F2 -> Ar2F* + F. Brau's "< 0.3" there is
  Chen & Payne's f2'/f3, the direct Ar2F* share, not an ArF* branching.
  Lee 1989's electron-quench values (1.6, 1.9, 2.4e-7) are fits to his
  amplifier data. NRL's memo quotes Adams & Dillon's Sandia gain as
  Electra's, and REBUILD's `electra_population_rows` are Orestes outputs
  (docs/arf_data/SCORING_DATA.md section C) (2026-09-30).
- **The tail table is generated, never edited.** `v0_2/data/
  ar_tail_ratio.txt` comes from `docs/arf_data/ar_tail.py` (ML's method on
  Ar) and is md5-pinned in `xsec.TAIL_FILE`. It is read far from the
  states that made it: those were beam-driven with Te 1.3-2 eV, and LIGHT's
  own input runs at Te 4-10 eV. The same holds for KrF's table, less
  strongly (2026-09-30).
- **Score device data with the device's own composition.** Every ArF
  discharge dataset found uses a He or Ne buffer, and several e-beam ArF
  lasers do too (SCORING_DATA.md section D). LIGHT_ArF is Ar/F2 alone; do
  not score it on buffered data (2026-09-30).
- **Emission is the emitter's density times its radiative rate, never
  d(events)/dt.** The gradient of the saved cumulative counters picks up
  output-grid spikes: Jancaitis's Ar2F* peaks at 430 and 540 A/cm2 came out
  identical that way (2026-09-30).
- **Seeded runs on LSODA, not BDF.** Where n_e decays to its 1 cm^-3 floor
  (SRI's proton-excited starts, 5 atm and 0.5 Torr F2), BDF's
  finite-difference Jacobian overflows and the run dies. `score_arf.py`
  runs its seeded starts on LSODA, equal to BDF to four digits where both
  finish (2026-09-30).
- **Never compare 1Dv0_2 with v0_2 where v0_2's Te passes 40 eV.** Fluid1D
  caps the mean energy at 60 eV (Te 40 eV), and there the 1-D centre
  leaves the 0-D model (Rokni's 2 and 4 atm: x0.90, x0.94). Above the
  30-eV table edge the two-temperature tail rows still follow T_tail = r
  Te, so v0_2's ionization keeps rising (at Rokni's 2-atm state, where
  T_tail/Te = 0.010, C2 is 3.2e-17 cm3/s at 30 eV, 1.8e-15 at 40 eV and
  1.3e-14 at 47 eV) while the ceiling holds the 1-D's (2026-09-30,
  numbers corrected in (c)).
- **Above 30 eV not every rate is held.** Rows on Maxwellian tables stop at
  the 30-eV edge, the tail rows follow r Te up to it, and power laws (R11p)
  follow their laws in v0_2 (the 1-D stage holds R11p at its 30-eV value,
  less than 0.1 % of Rokni's emission). `report.warnings` said "rates were
  held at the edge" and now says which (2026-09-30).
- **Never read an own-input score at a low pump without its Te and its
  Maxwellian twin.** Below the ionization fractions the tail table was
  checked at (n_e/N under 1e-6) its tail is cold, T_tail/Te = 0.009 to
  0.037 at Rokni's pumps, and the bulk runs to 21-94 eV where the Maxwellian
  EEDF and the beam drive sit at 1.5-2.2 eV. The failures there are the
  tail table's (entry (c)). `SCORE.md` and the overview print the Te and
  the `*_maxwellian` rows beside those scores (2026-09-30, corrected in
  (c)).
- **`Model._k` takes densities in cm^-3**, as `Model.rhs` passes them. A
  check that passes the normalised state with n_e in cm^-3 reads the tail
  table at the wrong n_e/N, pinned at its edge. That gave entry (a) its
  wrong C2 numbers (2026-09-30).
- **Maxwellian-EEDF runs on LSODA.** On the Maxwellian EEDF, BDF's
  finite-difference Jacobian stalled at a Jancaitis point for 17 minutes
  and more. `score_arf.py`'s `*_maxwellian` entries run on LSODA, which
  gives the same peaks as BDF to five digits where both finish
  (2026-09-30).
- **Measure the 1-D against the 0-D up to where their Te part, not to the
  ceiling.** Ar+, whose production is steep in Te, moves by 7 % while the
  two Te differ by 1 %. Entry (a)'s "1 % in every species" did not hold
  for Ar+ at the ceiling (2026-09-30).
- **Never place the map page with `\afterpage`.** In the accuracy
  overview it fires inside Appendix A's longtables: `\afterpage{\clearpage}`
  gives "Misplaced \noalign" and `\afterpage{\includepdf...}` hung pdfTeX
  for over 10 minutes. The map follows the ladder's float page, and Section
  13 follows the map (2026-09-30, entry (e)).

Inherited from LIGHT (the same code):

- **No tuning.** Never move a rate to bring an output nearer a measured
  efficiency, gain or yield. The scores in `docs/arf_data/SCORE.md` score;
  they never steer.
- **A band's default is a measurement or 0.** A rate with no measurement is
  0 by default, and its band runs from 0 to the largest published value.
  Never set a default from an estimate, an analog or a device score.
- **The vendored cross sections are pinned** (`xsec.FILES`). Re-vendor and
  re-pin deliberately.
- **dTg/dt = heat / C_v alone is wrong** when reactions change the particle
  count (F2 -> 2F); the model subtracts kTg sum(cv_s dn_s/dt).
- **No BLAS matmul on the Maxwellian tables** (`xsec._mean` multiplies and
  sums); **superelastic rates from the Klein-Rosseland cross section**;
  **no zero-length solve_ivp segments**.
- **Do not read Te after n_e hits the floor** (1 cm^-3): Te = U_e/(1.5 n_e)
  is solver noise there.
- **"Flows" are not a partition**; only "fate at t_end" sums to 100 %.
- **(3/2)Te per captured electron is the 'three_halves' closure only.** The
  default books the captured electron's exact energy.
- **`pytest.approx` passes any two numbers below 1e-12**; the tests compare
  through an `approx` helper with `abs=0`.
- **The 1-D model's tests need `1Dv0_2/tests/conftest.py`**: pytest never
  puts LIGHT_ArF_ebeam/ on sys.path for the non-identifier '1Dv0_2'.
- **Import Fluid1D only through `fluid1d.load()`.** The f2sweep venv holds an
  editable install of another checkout.
- **The 1-D stage reads rates at each substep's starting energy**, and the
  given power reaches each half-step over its own window.
- **Never quote a 1-D gap average without its mesh**: the wall cells are
  unresolved sheaths.

## 2026-09-30 (b): the electron-beam drive in the 1-D model

The owner asked: "Apply the e-beam drive to LIGHT_ArF_ebeam/1Dv0_2".

### Added
- `1Dv0_2/fluid1d.py`: the stage takes a species source per cell.
  `CellChemistry` builds the beam's products per pair in Fluid1D's columns
  from `beam.per_pair()`: one electron, one Ar+, the excited atoms and the
  F atoms, less the Ar and F2 they come from. `stage(...,
  species_source=)` passes it to the numba kernel, which adds it at each
  substep after reading that substep's mean energy, then holds the
  electron energy to its bounds with the electrons the cell now has. The
  beam's pairs per electron join the substep frequency. With a zero source
  the arithmetic is the own input's, bit for bit.
- `1Dv0_2/sim1d.py`: `Params1D.drive`, `beam_route` and `beam_sharing`,
  passed to every cell's v0_2. Under the beam the given power in a cell
  makes pairs = P / W. The stage takes their species, the electron energy
  their secondaries' share (4.4 eV per pair at the reference fill) and the
  gas their prompt heat (0.55 eV per pair, ledger `heat_beam`). The model
  line names the drive, the fate lists the beam's heat, and the summary
  says where the power went. `beam.py` joins the chemistry files, byte for
  byte v0_2's. Version 0.2.0-arf-ebeam+1d0.
- `1Dv0_2/run.py` (`run_1Dv0_2.py`): `--drive`, `--beam-route`,
  `--beam-sharing`.
- `1Dv0_2/tests/test_1d.py`: the stage's energy, atom and charge
  conservation with the beam's sources, the same stage with a zero source
  equal to the own input's, the beam's electrons entering a cell at the
  ceiling inside the bounds, the centre cell against v0_2's 0-D beam run
  (emission by band, ArF formed, prompt heat, species, Tg, Te), the books
  closing under the beam, and `--drive beam`. `tests/test_beam.py` is
  v0_2's, pointed at 1Dv0_2's copy.
- The score: `score_arf.py` runs the 1-D tests on both drives
  (`centre_beam` and `gap_beam` in `score_1d.json`). `SCORE.md`, the
  README's table (a "1-D centre, beam" column) and the overview's 1-D
  section show them. The overview adds the beam's Rokni points
  (`rkp.bdep2`, `rkp.bdep4`).
- `1Dv0_2/README.md`: the reference point under the beam on 64 and 128
  cells, and the scoring on both drives.

### Found
- **Under the beam the 1-D interior is v0_2's beam run.** The centre cell
  gives v0_2's beam numbers at every device point to 1.1e-4, at Rokni's
  five pressures to 5e-5 and at the reference point to 5e-6, on 64 and 128
  cells. On the own input the device points agree to 3.9e-4 and the
  reference point to 2e-5.
- **The own input's Rokni exception is gone.** Under the beam v0_2's Te
  stays at 1.5 to 2.2 eV, so Fluid1D's 40-eV ceiling, which takes the own
  input's 1-D centre to x0.90 and x0.94 of v0_2 at 2 and 4 atm, never
  binds. The 1-D Rokni rms goes from 0.24 to 0.14 over its four points.
- **The wall cells change their channel, not their size.** On the own
  input the walls drain the wall cells' electrons and the mean-energy
  ceiling cuts what the power leaves them, 3.1 % of the input at the
  reference point on 64 cells. Under the beam those cells keep receiving
  pairs, nothing is cut, and the walls take 2.3 % instead. Both halve with
  each mesh doubling (1.55 % and 1.18 % on 128 cells, 4.2 % under the beam
  on 32), as does the gap average's shortfall from the centre, 2.1 to
  2.4 % at the device points under the beam and 3.1 % on the own input on
  64 cells. Fluid1D's phantom negative work in those cells is 15 times the
  input under the beam on 64 cells (0.13 times on the own input) and falls
  x0.15 to x0.17 per doubling. It never enters the books, which close to
  3e-11 % of the input.
- The drive carries to 1-D as it does in 0-D: the 1-D centre under the
  beam has the 0-D beam column's scores (Electra x0.84 and x0.72, Keio
  x1.14, Sandia x1.82).
- The own input's 1-D runs are LIGHT_ArF's, bit for bit.

### Fixed
- `SCORE.md`'s introduction said the own input differs from LIGHT_ArF's
  at Rokni (18.8 in place of 20 kW/cm3 per atm). Entry (a) kept 20 so
  that the two records differ by the drive alone, and the own records
  agree to 5e-8 outside the seeded starts. Root cause: the sentence was
  written while the score ran Rokni at 18.8 and was not revised when 20
  came back. The introduction and the Rokni test's docstring lost their
  semicolons, and the two stored Rokni "doc" fields of `score_0d.json`
  were refreshed from the docstrings without a re-score (no other byte
  changed).
- `1Dv0_2/README.md` (from LIGHT_ArF) said the centre cell gives Rokni's
  numbers to 4 digits. On the own input it falls to x0.90 and x0.94 at 2
  and 4 atm. LIGHT_ArF's copy is corrected in its entry (g).
- `LIGHT_ArF/` paths where this copy's own belong: `1Dv0_2/run.py`'s usage
  text (which `--help` prints, its new `--drive beam` line included, though
  LIGHT_ArF's runner has no `--drive`), two test docstrings, the overview's
  build comments and the Known-bad line on `conftest.py`. Root cause:
  entry (a)'s rename pass did not cover `accuracy_overview.tex` or this
  changelog, and it did rename 1Dv0_2's usage text, but 1Dv0_2's files
  were then restored from LIGHT_ArF because the pass had touched
  `fluid1d.py`'s generated strings, and the restore took the renamed text
  back with them. Those generated strings still name
  `LIGHT_ArF/1Dv0_2/fluid1d.py` as the lineage, as does the tail table's
  header, which is never edited.

### Verified
- Tests in the worktree: v0_2 85 passed, 3 skipped, 1Dv0_2
  under the f2sweep venv 100 passed, 4 skipped. The skips are
  REBUILD's parity and the LoKI origin check, which need the main
  checkout, and the slow runs (`LIGHT_DOCS_SLOW`, and `LIGHT_1D_SLOW` in
  1Dv0_2).
- From the main checkout after landing: v0_2 87 passed and 1 skipped, and
  its docs tests with `LIGHT_DOCS_SLOW=1` 10 passed, the bands included.
  1Dv0_2 with `LIGHT_1D_SLOW=1` and `LIGHT_DOCS_SLOW=1` 104 passed, none
  skipped: REBUILD's parity, the LoKI origin check, the reference point in
  1-D, and the overview's runs file reproduced byte for byte, its
  own-input 1-D runs included. Both PDFs rebuilt from main's sources carry
  the text of main's PDFs.
- A mutation check planted 8 bugs in the beam's 1-D path (no electron per
  pair, the kernel ignoring the species source or adding it with the wrong
  sign, the secondaries' energy dropped, the prompt heat left out of Tg or
  of the books, pairs at the measured W, the bounds taken on the electrons
  before the source) and the tests caught all 8.
- Bit for bit: a short 8-cell own-input run on LIGHT_ArF's 1Dv0_2 and on
  this one agrees in every recorded number, and the 52 own-input values of
  `score_1d.json` equal LIGHT_ArF's.
- The README's accuracy table, all eight columns, recomputed from the
  records cell by cell. The 1-D Rokni record against v0_2 run at the 1-D's
  own pressures. The reference point on 32, 64 and 128 cells against
  v0_2's 0-D runs on both drives.
- The overview: rebuilt, 16 pages, no LaTeX error, undefined reference or
  box warning.

## 2026-09-30 (a): the electron-beam drive, and what it changes

The owner asked: "Look through Kinetic_Modeling/SPECK/REBUILD/ and see if
adding the electron beam drive in LIGHT_ArF_ebeam/ improves our models and
if it agrees with Morozov as well".

### Added
- `v0_2/beam.py`, the electron-beam drive, adapted from
  `SPECK/REBUILD/MODULES/beam.py`:
  - the Bethe/Moller stopping power with the ICRU 37 mean excitation
    energies, line for line, and from it the power a beam deposits in a
    thin cell (`deposited_power_W_cm3`) and its CSDA range.
  - REBUILD's Spencer-Fano recursion with its prompt floor at the lowest Ar
    threshold, run on the model's own vendored cross sections (the
    IST-Lisbon Ar and SIGLO F2 files, md5-identical to the copies REBUILD
    reads).
  - the secondary's energy from the Opal-Peterson-Beaty spectrum, Ebar(Ar)
    = 10.0 eV as LoKI-MC's `Databases/OPBParameter.txt` carries it (Opal,
    Peterson and Beaty, J. Chem. Phys. 55, 4100 (1971)), and the yields per
    pair as the least-squares slopes of one sweep over 20-90 % of the
    tables' ceiling. `yields(x, 'integral', 'equal')` is REBUILD's own rule.
  - the split onto LIGHT's species: every ionization one Ar+ and one
    electron (LIGHT has no F2+, and F2 makes 0.2 % of the pairs), every Ar
    excitation one Ar* with the rest of its threshold as prompt heat (the
    cascade the G1 and C1 rows book), a prompt F2 C1 or H1 excitation 2 F
    with the rest over D0 as heat. The secondaries' share closes
    W = sum(y E_form) + heat + eps_sec exactly and goes to the electron gas.
  At the reference fill: W 25.2 eV per pair, 0.39 excited atoms per pair,
  4.4 eV per pair to the electrons, 0.55 eV prompt heat. W is scored
  against the measured 26.4 eV (ICRU 31) and is 4.5 % under it.
- `v0_2/model.py`: `Params.drive = 'power'` (LIGHT's own input, the
  default, unchanged) or `'beam'`, with `beam_route` and `beam_sharing`,
  and the ledger's `heat_beam` term for the prompt heat. `report.py` books the
  beam's products (the Ar* reservoir counts its excited atoms as formed,
  pairs = E_dep / W exactly) and prints its split. `run.py`: `--drive`,
  `--beam-route`, `--beam-sharing`. Version 0.2.0-arf-ebeam.
- `v0_2/tests/test_beam.py`: the identity and atom balance for every route,
  REBUILD parity (runs where SPECK/REBUILD is in the checkout), W against
  the measurement, the equal-sharing oscillation and the Opal-Beaty
  convergence, the fit window, the stopping power, closures under the beam
  in five modes, and pure argon's one photon per ion and excitation.
- The score (`docs/arf_data/score_arf.py`): the beam column is the model's
  drive, the seeded starts (SRI, Chen & Payne) take its excited atoms per
  ion and secondaries' energy at their fill, and `beam.py` joins the
  chemistry fingerprint. `score_0d_LIGHT_ArF.json` keeps LIGHT_ArF's record
  as it was, and `SCORE.md` shows both.
- The overview (`docs/accuracy_overview.pdf`) names this copy, adds the
  drive to its bottom line, method, pure-argon section and open questions,
  and its numbers, figures, ladder and map are regenerated.

### Found
- **The drive improves the model overall.** Over the 11 held-out tests the
  mean |ln gm| is 0.314 for LIGHT_ArF's own input and 0.256 under the drive,
  and the rms of the logs over their 87 points 0.476 and 0.392 (LIGHT_ArF's
  prototype 0.275 and 0.429). Rokni's x25 current ratio goes from x1.67 to
  x0.97, his pressure curve's rms from 0.15 to 0.10, Electra from x0.64 to
  x0.84. Sandia goes from x1.39 to x1.82, Keio at its printed pump from
  x0.86 to x1.14, and Jancaitis's ArF* scan from rms 0.28 to 0.35. The drive
  makes ArF(B) at 0.25 to 0.30 of the input at every device pump, where the
  measured share falls from 0.36 to 0.16.
- **It does not reach Morozov.** Pure argon goes from x2.15 to x1.75 (x1.70
  at 500-1400 hPa). The model's one-photon ceiling with the drive is 0.538,
  near the 0.53 theoretical maximum Morozov quotes (Lorents 1976), and his
  argon sits at 62 % of it while his Kr and Xe sit near theirs. Run as a
  continuous beam the model is flat at 0.54 from 151 to 1000 hPa (x1.65 to
  x2.36). No input closes this: the network lacks a measured argon VUV
  loss, or his argon value is low.
- **REBUILD's yields oscillate** (Known-bad). The device scores move by at
  most 4 % between REBUILD's rule and the Opal-Beaty one (Electra x0.81 and
  x0.84), Morozov at 1000 hPa from x1.78 to x1.67. REBUILD is not changed
  here.
- Rokni's deposition, 5 A/cm2 at 150 keV, is 18.8 kW/cm3 per atm of Ar by
  the drive's stopping power on a straight ray. The survey's estimate, 20,
  is kept so that the two records differ by the drive alone.
- The tail table (`v0_2/data/ar_tail_ratio.txt`) was made under
  LIGHT_ArF's prototype drive. It is indexed by n_e/N and the F2 fraction,
  not by the drive, and the secondaries' share barely moved (4.40 to
  4.37 eV), so it is kept. `ar_tail.py` still imports the prototype.

### Verified
- The drive: its integral route with equal sharing reproduces REBUILD's
  `spencer_fano_yields` to round-off (W 29.345 and 30.851 eV, excited
  atoms 0.7364 and 0.8441), and its stopping power REBUILD's. The
  Opal-Beaty yields converge with the grid (W to 0.01 %, excitations to 1 %
  on 2000, 4000 and 8000 points).
- The model: every beam run closes to below 1e-12 of E_dep, and the Ar*
  reservoir to 1e-15. The own input's results are LIGHT_ArF's to 4.5e-8
  (the extra ledger state moves BDF's steps), and the Maxwellian runs are
  identical.
- Tests: v0_2 (see its README for the count) pass, with the parity test
  skipped in the worktree (no SPECK there) and run in main.
- The overview: rebuilt with no LaTeX error, undefined number or box
  warning.

### Not done
- 1Dv0_2 has no beam drive. Its stage adds power to the electron energy
  inside Fluid1D's numba kernel, and a drive there needs species sources
  in that kernel. Its interior equals the 0-D run, so the 0-D score stands
  for it. Entry (b) adds the drive.
