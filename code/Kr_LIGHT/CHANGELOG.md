# Changelog: LIGHT_Kr (the lightweight pure-Kr model, in preparation)

Dated entries, newest first. Read the Known-bad list before editing.

LIGHT_Kr is the owner's copy of `LIGHT/` (2026-10-13), made for a pure-Kr
model. The KrF models and their history are in `LIGHT/` and
`LIGHT/CHANGELOG.md`; this file starts with the removal of the KrF copies.

## Known-bad: do not reintroduce

v0_1's code is LIGHT v0_5's with the fluorine taken out, and v0_2's is
v0_1's with the ground-state levels resolved, so the Known-bad list in
`LIGHT/CHANGELOG.md` applies to both, except its entries on the fluorine
chemistry, the 1-D models and the documents. LIGHT_Kr's own:

- **Jung et al. 2006, as printed (2026-10-13).** Table II's 1s5 -> 2p6 at
  4 eV is 23 (x 1e-16 cm2), not 21, and Table IV's last row is 1s5 -> 2p5
  (dipole-forbidden, Table II's 2p5 column), not a core-changing 1s3 ->
  2p5: the paper extracts no core-changing cross section out of 1s3.
  LIGHT/v0_5/jung2006.py still carries both. `test_jung2006.py` holds the
  tables to the paper entry by entry.
- **Never compare a cross section through pytest.approx's defaults.** Its
  absolute tolerance of 1e-12 passes any two numbers of 1e-16 cm2; LIGHT's
  test of jung2006's tables did, which hid both misreadings above and a
  wrong expectation of its own (the 1s5 -> 2p8 extrapolation from the 8-eV
  point, where the column ends at 10 eV). Pass abs=0.
- **Compare v0_1 with v0_5 on the same tables.** v0_5 carries the two
  misreadings; `test_v0_1.py` gives it v0_1's tables in memory
  (`v0_5_on_v0_1s_tables`) and restores them.
- **Under the kinetic EEDF, run-level parity is the closure's own spread,
  not round-off.** The table is renewed where the solver's steps put the
  drift triggers, so a state vector of another size (v0_5 carries its idle
  fluorine) or another rtol moves the 147-nm share by 3e-4 and the elastic
  share by 1.2e-3. Test the machinery from one state, where v0_1's table is
  v0_5's bit for bit.
- **On pure krypton's S curve (near 0.04 Td) a fixed-field iteration finds
  the upper solution or none.** Use the mean match (Solver.eedf), which
  gives one solution from any start. LIGHT's bistability example at 0.4 Td
  needed F2's low-energy losses for its cold branch.
- **Never start the kinetic grid search from a grid far too long in pure
  Kr.** From 40 or 60 eV the EEDF's tail falls to the double-precision
  floor, the clipped top cell reads 300 decades, and the search keeps the
  long grid (ionization rate 31 % off; a strict xfail test). Runs start at
  30 eV and warm-start, and at the reference point the steady solutions
  outside the window are on the 12-eV minimum grid with cold electrons
  (25-34 decades), not at the noise floor.
- **row_scale under the kinetic EEDF (fixed 2026-10-13 (c)).** The rows
  whose rates the kinetic table carries (G1, C1a-c, C2, R11p, and SUPER
  when added) read the table's rate without their factor, so row_scale
  moved none of them on the default closure. Their factor now goes into
  their process in the Boltzmann equation, as K76's band always did.
  `test_row_scale_reaches_the_kinetic_table` holds it. G6, the elastic
  exchange, cannot be scaled (refused), and K76 and K76r take one factor
  (2026-10-14). LIGHT v0_4 and v0_5, LIGHT_ArF v0_3 and v0_4 and this
  directory's v0_5 copy carried the defect too and were fixed the same
  way on 2026-10-14.
- **The low-field ionization coefficient is not a rate to fix
  (2026-10-14).** Below 50 Td every published Kr cross-section set falls
  below Kruithof's and Heylen's alpha/N (Bordage et al., J. Phys. D 46,
  334003 (2013): impurity Penning ionization in the experiments), and the
  solver agrees with Monte Carlo there (LoKI-MC, within 4 % on v0_2's
  file). Never tune a cross section or the band hm_kr toward those
  points: hm_kr stays 0 until a krypton measurement of the association's
  share sets it, and at its top it overshoots the 100-300 Td data by 39 %.
- **The Hornbeck-Molnar rows are first order in Kr (v0_2).** For the
  'assoc' kind `Compiled.heavy` holds one Kr: the partner's density is in
  the share. Giving the rate both of the row's Kr reactants would make it
  N times too large.
- **v0_2's stiff solver takes the Jacobian's structure
  (`Model.ivp_options`).** Only the densities, U_e and Tg feed the
  right-hand side. A new state component that does must join `live`, or
  the Newton iteration's Jacobian misses its column.
- **The two-term swarm is not what a swarm experiment measures
  (2026-10-15).** A drift tube times the bulk drift velocity, which lies
  above the two-term flux one as ionization grows (+23 % at 300 Td). A
  Townsend experiment measures the SST coefficient, 9 % below the
  two-term k_i/W at 300 Td (11 % below the Monte Carlo's flux one). The
  two-term transverse diffusion is 46-79 % high from 10 to 300 Td
  (LoKI-MC on Biagi's file, `v0_3/swarm_mc.py`).
  Score and train against `swarm_mc.as_measured`, never the bare
  two-term numbers. v0_2's overview compared k_i/W with Townsend data
  and so put the model closer to them at high field than it was.
- **Never choose a swarm train set by leave-one-out (2026-10-15).** A
  dataset that alone pins a factor predicts badly when left out, however
  well it agrees, so leave-one-out reads "informative" as "inconsistent".
  It dropped Hunter 1988, the most precise drift velocity, first. ML's
  consistency check drops by the chi2 at the joint fit
  (`ML/ml/kr_swarm.consistency`) and shows leave-one-out only for
  information.
- **v0_3's cross sections are the trained ones, by default
  (2026-10-15).**
  `xsec.processes` applies the pinned factors (`TRAINED_FILE`) and
  `raw_processes` is Biagi's file. A comparison with Biagi's own file (a
  Monte Carlo, ML's dense solver, LoKI-MC's input) runs under
  `xsec.using(xsec.UNIT)`. With unit factors v0_3 is v0_2 number for
  number (`test_v0_3.py`).
- **Retrain, never hand-edit, the factors (2026-10-15).**
  `python3 ML/run.py kr-train` writes `data/kr_swarm_trained.json` and
  prints its md5 to pin. The file's record (train set, consistency
  rounds, emulator accuracy, posterior draws, variants) is what the
  documents read. A change to the datasets, the error model or
  `swarm_mc`'s W and alpha ratios means retraining. Its D_T ratios
  train nothing; they reach the overview's variant table, which scores
  the variants live.
- **v0_3 trains seven factors, not nine (2026-10-15).** The momentum
  transfer at 0.03 and 0.2 eV stays Biagi's. Freed, it traded against the
  Ramsauer minimum (1.13 at 0.03 eV, 0.88 at 0.6 eV), and Koizumi's
  held-out D_T/mu went from chi2 2.5 to 8.1 and Pack's 195-K drift
  velocities from 3.1 to 5.5 and from 6.1 to 8.7. Those sets chose this,
  so the documents call them validation sets, and the beam's partition
  and the low-field alpha/N are the tests nothing chose on. The variant
  comparison is in the factors file's record (`free_variants`) and the
  overview's "Which factors train" table. Change it through ML's
  `kr_swarm.FREE` and retrain.
- **A kinetic EEDF that fails from a warm start has not shown that the
  state has no solution (2026-10-15).** In a beam's afterglow (0.5 atm,
  3e5 W/cm3, n_e 1e16 cm-3) a table point solved without the growth term,
  warm-started from the last table's point that had it, found none. A
  cold start converges there, to the one solution the mean match has.
  `Solver.eedf` restarts cold once and counts it (`Solver.restarts`, the
  overview reports the total). Never remove the restart, and never take
  a failed warm start as the state's answer. LIGHT v0_5 and LIGHT_ArF
  v0_4 carry the same table code without it.

- **A table point told to leave the growth term out can have no steady
  solution (2026-10-15).** In the afterglow of very hot gas (the scan's
  0.25-atm, 500-kW corner, 3970 K) Penning ionization made electrons 4.6
  times faster than recombination took them; a side point of the kinetic
  table below the no-field mean, without the growth term because its
  centre had none, found no solution on its centre's grid from any start
  and stopped the run. `Solver.eedf` now falls back to the growth term on
  the same grid, then to a grid of its own (`Solver.fallbacks`, logged
  per update). Never remove the fallback, and never read a stop there as
  the state having no EEDF: with the growth term it has one.
- **A kinetic table point with no solution at all must not stop the run
  (2026-10-16, v0_5).** After the restart and v0_4's fallbacks (growth
  term, a grid of its own) a side point can still have none: the scan's
  0.25-atm, 500-kW corner, n_e 1.05e17 cm-3 at 2243 K in the afterglow,
  at 0.97 eV. `model._kinetic_update` leaves such a point out of the table
  (`dropped`) and, where the centre itself has none, lets LIGHT's
  Maxwellian stand in until the next renewal (`failed`); the documents
  count both. Never turn either back into a raise: the batch lost 4.8 h to
  one. LIGHT v0_5 and LIGHT_ArF carry the update without it.
- **The two-temperature closure's tail table is the drive's (v0_4).**
  v0_1 to v0_3 read ML's beam-driven KrF table (its F2 = 0 column) for a
  field-driven discharge. v0_4's default is `data/kr_tail_field.txt`
  (`tail_field.py`, field-driven pure Kr from this model's own kinetic
  EEDF); the beam table stays as Params.tail = 'beam' for comparison, never
  as the default. Regenerate the field table (`python3 -m v0_4.tail_field`,
  every point from a cold start) and re-pin its md5; never edit it.
- **Never tune the field tail table toward the kinetic run (2026-10-15).**
  Its states carry no excited atoms, so its r is 16 % below the kinetic
  EEDF's own along the reference run (the 5s <-> 5p collisions heat the
  tail), and the closure's 147-nm share comes closer to the kinetic one
  partly by compensation (with the run's own r it would be farther: the
  Maxwellian bulk errs the other way). The fix is an excited-fraction
  axis, not a factor on r.
- **Janssens et al. 1987 is a calculation (2026-10-15).** Their 3P1 loss
  rate (1.7e5 + 2050p + 82.0p^2 /s, p in mbar) comes from a quasiclassical
  theory compared with their experiments: its 82.0p^2 (1.36e-31 cm6/s) is
  cited beside kr2_form_1s4's band, not a band end, and its 1.7e5 /s is
  the escape band's top only as the pressure-independent term of that
  rate (diffusion and resonance radiation together).
- **Turner's 1s4 -> 1s5 law is (T/295)^-1.29 through his 295- and 196-K
  points (2026-10-15).** The rows carry its 300-K value, 1.272e-13, on
  (Tg/300)^-1.29 (`network.MIX_1S4_300`), not 1.3e-13 at 300 K.
  `test_the_1s4_1s5_transfer_reverse_is_turners` holds both points.
- **The Kr* + 2Kr law is measured over 200-400 K (2026-10-15).** v0_1 to
  v0_3's documents and SOURCES.md said 300-451 K, the range of the
  power-law conversion; `report.KR2_TG_VALID_K` is 400.

## 2026-10-16 (c): v0_5, the electrons on the excited atoms and Kr2*, and its documents

`v0_5/` is v0_4 with the scope of (b); `run_v0_5.py` launches it.

### Added (code)
- **`excited.py`**: Vriens & Smeets's ionization cross section and rate;
  pseudo cross sections that carry A Te^B exp(-C/Te) or a constant
  superelastic rate into the kinetic EEDF exactly; the 5p rung's ten
  levels (Biagi's energies, 2J + 1) populated at the gas temperature
  (`rung_fraction`, `rung_partition`).
- **`network.py`** (219 rows, v0_4 126): the 5s -> 5p pumps pair by pair
  (`K76_<1s>_<2p>`, 26, each charged its own threshold, the excess over
  2p10 to heat) and their reverses out of the rung's levels
  (`K76r_<1s>_<2p>`, `Row.rung_level`, `Row.e_gain_eV`, the levels'
  excess lent by the gas); the 5p lines one by one (`PRAD_<2p>_<1s>`, 28:
  v0_4's 15 plus NIST's 975-nm 2p10 -> 1s4 and the weak core-changing
  lines); the rung's relaxation reversed (`KPSr_*`); the stepwise
  ionizations on Vriens & Smeets (`KSTEP_S_<1s>`, `KSTEP_P_<2p>`, bands
  kstep_s/kstep_p now scales 1 by default, 0 to 2); e + Kr2* (`EI_KR2`,
  `EI_KR2F` on Vriens & Smeets, band ei_kr2; `EQ_KR2` band eq_kr2 0 to
  Kannari's 1e-7; `EM_KR2`/`EM_KR2r` band em_kr2 0 to 1e-7); the 5s
  levels' superelastic return (`SUPER_<1s>`, Klein-Rosseland of G1); e +
  1s5 <-> 1s4 (`EMIX_1s5`, `EMIXr_1s4`, Blagoev & Pilosov's 4.3e-8, band
  emix_5s to 2.3e-7). `network.CALCULATED`; `LADDER_KANNARI` = kps_kr
  alone.
- **`boltzmann.py`**: `PopProcess.lower`/`kind`: one direction alone, an
  ionization out of an excited population (both electrons share the
  rest; the growth term counts it); the operators summed per population
  (`Grid.pop_ops`: 74 processes, 16 populations).
- **`model.py`**: every electron row on an excited population is a
  population process of the kinetic EEDF; the rung's level shares in the
  RHS; the Kr2* fraction among the kinetic EEDF's drift triggers;
  `Params.tail = 'field0'` (the field table's x_s = 0 block, v0_4's).
- **`tail_field.py`, `data/kr_tail_field.txt`** (md5 c04a5d06...):
  r(n_e/N, Te, x_s), 14 x 17 x 7 = 1666 points, the excited populations
  of each point from the model's own rates (the 5s levels by weight, the
  rung and Kr2* in steady state); its x_s = 0 block is v0_4's table to
  the digit.
- **Data**: Kr2+ bound by Signorell et al.'s 9267.8 cm-1 (1.1491 eV);
  Kholin et al.'s 3P0 rate cited with its argon third body.

### Found while building
- **The tail generator held a hot point at the cold floor.** At low Te
  with excited atoms the superelastic electrons feed the tail at 10 eV and
  the kinetic tail removes up to 440 times the Maxwellian's energy, which
  a tail tied to the bulk at the break cannot reach; the first build gave
  those 175 points r = 1e-3. Fixed: hotter than the form follows is r's
  top, 3 (`tail_field.solve_r` returns inf; the file header counts both).
- **The kinetic power balance is the grid's**: 5.7e-5 on 600 cells (grids
  to 31 eV, where v0_4's ended near 19), 1.6e-5 on 1200 (a 40-ns pulse).
  The start-up transient (6.7e-2 at t = 5 fs, Te = Tg) is v0_4's too.
- **Each population operator gives or takes exactly rate x u** of the
  electrons' energy on a solved EEDF (checked for all 74; a test).
- **The superelastic returns against Blagoev et al. 1984**: 5.3-6.6 times
  their thermal rate for 3P2, 5.5-9.5 for 3P1, 0.55-1.24 for 3P0 (the
  cross sections' first 0.2 eV above threshold); the returns move the
  reference point by -0.5 % (first written -0.4 %).
- **Runs**: the reference point 6.08 % to 147 nm on the kinetic EEDF
  (v0_4 31.22 %), 5.58 % on the two-temperature closure (6.30 % on v0_4's
  table, without the excited axis); 28 % of the excimers radiate, the
  electrons ionize most of the rest; n_e 8.2e14 cm-3 in the pulse (v0_4
  2.8e14). Kinetic runs take about twice v0_4's time per update.
- **The hot corner stopped the documents' batch (fixed).** The scan's
  0.25-atm, 500-kW corner on the kinetic EEDF reached n_e 1.05e17 cm-3
  (1.7 % ionized) at 2243 K in its afterglow, where a side point of the
  table (0.97 eV) had no solution with or without the growth term, even
  cold on a grid of its own, and the RuntimeError ended the run after
  hours (batch 1, 191 of 193 runs done, 4.8 h; the in-flight 1200-cell run
  was lost with it). Root cause: v0_4's fallbacks end in a raise. Fix
  (`model._kinetic_update`): a table point with no solution is left out of
  the table (logged `dropped`; the table reads that side from its next
  point and LIGHT's Maxwellian beyond it), and an update whose centre has
  none leaves LIGHT's Maxwellian (the table None) until the next renewal
  (logged `failed`); `kr_docs_v0_5.summarize` counts both, the overview
  prints them. The change sits in exception handlers only, so every run
  that completed before runs identically after; all 193 are rerun anyway
  (the code hash), longest first (`batch 2`). Gate:
  `test_the_table_survives_a_point_or_a_centre_with_no_eedf`.
- **Vriens & Smeets's two formulas agree to 0.2 % at 2-5 eV, 0.8 % at 1
  eV** (an independent integration, 2026-10-16 (c)); (b) had said 0.1 % at
  1-5 eV. Corrected in (b), SOURCES.md and the gate's docstring.
- **What dominates v0_5's uncertainty: the electrons on Kr2*.** On the
  pre-fix runs (unchanged by the fix): the excimer's ionization 15.43 %
  at 0 to 4.47 % at twice the formula (term 5.5 points, 61 % of the
  quadrature's square); its quenching by electrons at Kannari's 1e-7, 1.27
  %; the mixing of its two states at 1e-7, 12.44 % (the fast state
  radiates in 3.4 ns, before the ~74 ns the electrons take to ionize it).
  Quadrature +-7.0 points (+-115 %). The J = 1 proxies now lower the
  share when removed (5.55 %; v0_4 54.81 %).
- **Batch 2 (the fixed code, md5 b392d33a7c61 over every module):** 193
  runs, 0 failed, 4.0 h on 12 workers longest first. Every one of the 191
  runs batch 1 had finished came out identical in every output (timing and
  the new counters aside), as the fix's confinement to exception handlers
  says it must. The hot corner ran to its end in 4.0 h: 0.0053 % to 147 nm
  (two-temperature 0.0050 %), 3281 K, closure 1.8e-15, 11 table points
  dropped, no failed update, 40 restarts and 42 fallbacks. The 1200-cell
  run moves the reference point by 1.3e-3 of itself (the documents'
  largest numerical term).

### Added (documents)
- `docs/kr_docs_v0_5.py` (193 jobs, 72 kinetic; `summarize` counts the
  dropped points and failed updates), `make_system_map_v0_5.py` (the
  electrons' share of every ionization drawn as C2's box, the excited
  species' own energy as their bands; the rung's heat band carries KPS net
  of KPSr, the pumps' excess and the lines' loan), `make_overview_v0_5.py`
  (`claims()` holds 26 checks of the new text's statements against the
  runs, the parts 27 more, and the generator raises when one fails; the
  v0_4 runs are read from v0_4's own cache).
- `docs/system_map_v0_5.tex/.pdf` (one page) and
  `docs/accuracy_overview_v0_5.tex/.pdf` (23 pages): v0_5's changes and
  what each does undone, the reference point, the scan, the closure's
  excited axis, a section on the electrons and the excited atoms (Vriens
  & Smeets against their rate and Kannari's fits; the superelastic score
  against Blagoev et al.; Blagoev & Pilosov's mixing; the rows at 0; the
  rung at Tg), the beam, the ladder with v0_5's rows, the limitations, and
  the rates in two tables. New sources, each with a Crossref-checked DOI:
  Vriens & Smeets 1980, Flannery, McCann & Winter 1981, Blagoev & Pilosov
  1981, Qiu et al. 2018, Zeng et al. 2005, Whitehead et al. 1995; Ng et
  al. 1977 dropped (Kr2+'s well is Signorell's).
- `v0_5/README.md` (its numbers filled from the overview's numbers file)
  and `LIGHT_Kr/README.md` (v0_5 the current model, its commands and
  reference point).

### Found (the documents' numbers)
- **Reference point:** 6.08 % to 147 nm on the kinetic EEDF, 5.58 % on the
  closure (v0_4 31.22 and 32.92 %), a factor 5.1. The electrons ionize
  56.9 % of the excimers formed (in ~74 ns against 270 ns of radiation)
  and 28.5 % radiate; by number the excimer's ionization makes 45.2 % of
  the ions, the 5s levels' 33.8 %, the rung's 16.0 %, Penning 4.9 %. With
  the three calculated rows at 0 the share is 28.02 %; without the
  excimer's ionization 15.43 %.
- **Uncertainty:** +-7.0 points (+-115 %), the three largest terms the
  electrons on Kr2*: its ionization 15.43 % at 0 to 4.47 % at twice the
  formula (61 % of the quadrature's square); its quenching at Kannari's
  1e-7, 1.27 %; the mixing of its two states at 1e-7, 12.44 %.
- **Scores:** swarm and electroluminescence as v0_4 (0.99-1.04 of Mano et
  al.'s fit); superelastic return 5.3-6.6 x Blagoev et al. for 1s5;
  Morozov +29 % (v0_4's 0.543 to the digit, 0.005 excimers ionized per
  pair); Eckstrom 1 atm -3 % (v0_4 +1 %), and at 0.5 atm the new rows cost
  29 % of v0_4's light.
- **Scan:** the share rises with pressure and falls with power everywhere;
  best 19.85 % at 2 atm and 100 kW (v0_4 73.6 %).
- **Closure:** its excited axis puts the tail at 0.96-0.99 of the kinetic
  one along the reference run (the zero block, v0_4's table, 1.16-1.30);
  the share then lies -8 % from the kinetic EEDF (zero block +4 %): what is
  left is the Maxwellian below the break, which the zero block's hotter
  tail had more than offset. The corners within 6 %.

### Verified
- `/usr/bin/python3 -m pytest -p no:cacheprovider LIGHT_Kr/v0_5/tests`
  (without test_docs) from Kinetic_Modeling on the fixed code: 417 passed,
  1 skipped (1 h 52 min); on the code before the fix 416 passed, 1 skipped.
- `LIGHT_KINETIC_SLOW=1 /usr/bin/python3 -m pytest LIGHT_Kr/v0_5/tests/test_docs.py`:
  7 passed (9 min 42 s): both generators reproduce their numbers files
  byte for byte, every number printed is generated from the current code,
  every citation is cited and has a DOI, no local file is named.
- Both PDFs built clean (no overfull box but the map's 1.6-pt one, which
  v0_4's map has too) and read page by page; the text's claims raise in
  the generator if they fail (all hold).
- The Vriens & Smeets cross section against their rate formula by an
  independent integration (above); every new source's metadata on
  Crossref.

### Not done (in scope, deferred) and limitations found
- The J = 1 levels' 5s -> 5p cross sections stay proxies: Zeng et al.'s
  R-matrix values, Qiu et al.'s 5s transfer rates and Whitehead et al.'s
  5p deactivations are behind paywalls or bot checks (the owner's choice:
  no full texts).
- e + Kr2* -> e + 2Kr and the electrons' mixing of Kr2*'s two states have
  no published value: bands at 0, and they are now two of the three
  largest uncertainty terms. A measurement (or calculation) of e + Kr2*
  collisions is what v0_5's answer most needs.
- The field-driven regime has no pure-krypton score.
- LIGHT v0_5 and LIGHT_ArF v0_3/v0_4 carry the kinetic table update
  without the dropped/failed handling (a task chip was offered).

## 2026-10-16 (b): v0_5 scope, the owner's choice

Shown the groundwork below, the owner chose for v0_5 (2026-10-16):
- **defaults from published calculations** where nothing is measured
  (recommended over v0_4's rule of 0): Vriens & Smeets 1980 for the
  stepwise ionization out of the 5s levels and the 5p rung, and the same
  form, which Flannery et al. 1981's ~1e-15 cm2 supports, for the
  ionization of Kr2*; each band from 0 to twice the calculation.
  e + Kr2* -> e + 2Kr, with no published value, enters as a band at 0;
- **the pumping energy**: each measured 5s -> 5p pair its own row, charged
  its own threshold, and the 5p rung's relaxation with its reverse;
- **electrons and excited species**: the e + Kr2* rows (ionization,
  quenching, the mixing of its two states), the stepwise ionizations as
  cross sections in the kinetic EEDF, the 5s levels' superelastic return
  (scored on Blagoev et al. 1984) and their mixing by electrons (Blagoev &
  Pilosov 1981);
- **data fixes**: Kholin et al.'s 3P0 rate has argon for its third body
  (cite corrected), Signorell et al.'s Kr2+ constants;
- **two-temperature numerics**: the tail table's excited-population axis;
- **no full texts behind paywalls or bot checks** (recommended): v0_5 works
  from the abstracts, the J = 1 sets stay proxies with their band.

Checked before building: the Vriens-Smeets ionization cross section used
for the kinetic EEDF, 4 pi a0^2 R^2/(E + 3.25 E_n) (5/(3 E_n) - 1/E -
2 E_n/(3 E^2)), integrates over a Maxwellian to their rate formula 9.56e-6
Te^-1.5 exp(-e)/(e^2.33 + 4.38 e^1.72 + 1.32 e), e = E_n/Te, to 0.2 % for
Te of 2-5 eV, 0.8 % at 1 eV and 3 % at 0.5 eV (first written "0.1 % for
1-5 eV"; corrected 2026-10-16 (c) by an independent integration), for
E_n = 4.085, 3.64 and 2.70 eV. Its peak
for Kr2* (E_n 3.64 eV), 1.0e-15 cm2 at 9.9 eV, is Flannery et al.'s order.
Kannari's printed fits of the same rows ([72], [79]) lie within 10-40 %
of the formula over Te 1-3 eV.

## 2026-10-16 (a): v0_5 groundwork, v0_4's limitations surveyed (no code change)

The owner asked: "Begin approaching the limitations of v0_4 to make v0_5".
Each limitation of v0_4's overview (its section "Where v0_4 is not yet
reliable") and of the (d) entry's "Not done" was measured on scratch
copies of v0_4 or searched in the literature (two search agents; every
DOI below checked against Crossref; abstracts unless stated). Nothing in
v0_4 changed. What v0_5 builds is the owner's choice.

### Found, measured on v0_4 (scratch copies; reference point, kinetic 31.22 %, two-temperature 32.92 %)
- **The 5s -> 5p pumps are undercharged.** The ledger books every pump at
  the rung's 2p10 energy, but the measured cross sections reach the ten
  5p levels at their own energies. The rate-weighted gap is 1.540 eV out of
  1s5 (ledger 1.388), 1.420 out of 1s4 (1.271), 1.561 out of 1s3 (0.741)
  and 1.481 out of 1s2 (0.660): about 18 % of the input missing from the
  cycle's heat over the reference run. One row per measured pair, each
  paying its own threshold, the excess over 2p10 heating the gas: 27.67 %
  on the kinetic EEDF (-11.4 %, closure 3e-16), 29.17 % on the
  two-temperature closure. No new data needed.
- **The rows v0_4 leaves at 0 or does not carry decide the field-driven
  answer** (two-temperature closure):
  - the stepwise ionizations at Vriens's calculation (Kannari's printed
    fits, which are Vriens & Smeets's formula to 10-40 %): out of 5s
    21.17 %, out of 5p 25.52 %, both 19.22 % (Te 2.65 -> 1.82 eV);
  - e + Kr2* -> Kr2+ + 2e, not carried, at Kannari's 9.0e-8 Te^0.70
    exp(-3.53/Te) (Vriens & Smeets's form): 9.48 %, 30 % of the excimers
    radiating;
  - e + Kr2* -> e + 2Kr, not carried, at the KrF models' estimates
    (1e-9 to 1e-7 cm3/s, 100x apart; in-repo LIGHT survey of 2026-09-27):
    31.34 % at 1e-9, 20.56 % at 1e-8, 4.13 % at 1e-7;
  - every estimate at once (per-pair pumps, the superelastic return, both
    stepwise ionizations, the excimer's ionization): 7.13 %.
  Under the beam drive the stepwise ionizations at Vriens's values change
  nothing (Morozov 54.17 %, Eckstrom at 0.5/1/2 atm 13.04/46.40/53.67 %):
  the beam's electrons are cold, so no score in hand tests these rows.
- **Small:** the superelastic return of the 5s levels to the ground state
  (Klein-Rosseland on each level's trained cross section), kinetic 31.22
  -> 31.09 %; electron transfers between the 5s levels, every pair at
  A exp(-dE/Te) with detailed-balance reverses, A = 2e-8 to 1e-6 cm3/s,
  +0.3 to +1.2 %; electron mixing of the excimer's slow and fast states,
  1e-8 to 1e-6 cm3/s, +0.7 to +1.5 %; the 5p rung's relaxation with its
  detailed-balance reverse, no 147-nm share on the scan moves, the hottest
  corners' heat does (0.1 atm, 500 kW: 82.9 -> 75.8 % of the input, Tg
  8549 -> 7866 K).

### Found, literature
- **Electrons between the 5s levels:** e + 3P2 -> 3P1 measured at
  4.3e-8 cm3/s (Blagoev & Pilosov, Phys. Lett. A 83, 161 (1981); Te not in
  the abstract); 1s5 -> 1s4 and 1s3 -> 1s2 at Te 0.07-1 eV (Qiu, Cheng,
  Zhu & Pu, J. Phys. D 51, 135202 (2018); numbers behind IOP's bot check,
  not bypassed). Calculations: Dasgupta et al., PRA 65, 042724 (2002);
  Allan, Zatsarinny & Bartschat, J. Phys. B 44, 065201 (2011).
- **Blagoev, Mishonov & Popov 1984's thermal rates** (3P2 1.5e-10, 3P1
  5.2e-10, 3P0 9.8e-10 cm3/s) are the superelastic return to the ground
  state, not 5s mixing: a score for the superelastic rows.
- **Out of the J = 1 levels:** one calculation covers every 5s -> 5p
  transition, the fully relativistic R-matrix of Zeng, Wu, Jin, Zhao &
  Yuan, PRA 72, 042707 (2005), within about 10 % of Jung's metastable
  values for four transitions after an energy shift; its numbers are behind
  APS's bot check. "Dipti & Srivastava, JQSRT 2016" is about Ar+, not Kr.
- **Ionization of Kr(5s):** no measurement for Kr or Ar (only He 2^3S).
  Calculations: Ton-That & Flannery, PRA 15, 517 (1977); Vriens & Smeets,
  PRA 22, 940 (1980).
- **Ionization of Kr2*:** Flannery, McCann & Winter, J. Phys. B 14, 3789
  (1981), binary-encounter theory: about 1e-15 cm2, 40-50 % above the
  metastable atom's. e + Kr2* -> e + 2Kr and the electron mixing of its
  two states: nothing, for Kr2*, Ar2* or Xe2*.
- **Kr(5p) + Kr per level:** Chang, Horiguchi & Setser, JCP 73, 778
  (1980) measured Kr(5p) in argon, not krypton (2p10 relaxing fast to 5s
  there); in krypton Whitehead et al., JCP 102, 1965 (1995) give per-level
  and state-to-state constants (paywalled).
- **Kr(1s3) + 2Kr:** Kholin et al. 2015 measured Kr(3P0) + Kr + Ar ->
  Kr2* + Ar, (5.3 +- 0.8)e-33 cm6/s, argon the third body: no pure-Kr
  value, and v0_4's cite ("the quenching of 3P0 by Kr") is inexact.
- **Kr* + Kr*:** Blagoev et al. 2006's number is in its PDF only (IOP's
  bot check); the Kolokolov & Blagoev review, Phys. Usp. 36, 152 (1993),
  is an open PDF, not fetched.
- **Nothing found:** Kr2* + Kr2*, Kr* + Kr2*, Kr3+ kinetics, Kr* + 2Kr
  above 400 K, the Kr2* A state's bond length.

## 2026-10-15 (d): v0_4, the 5s levels resolved, and its documents

The owner's request of 2026-10-15 ("Begin approaching the limitations of
v0_3 such as from section 7 of the accuracy overview to begin making the
improved v0_4"), built to the scope chosen in (c). `v0_4/` is v0_3 plus
what follows; nothing else changed (`tests/test_v0_4.py`'s gate).

### Added
- **The four 5s levels** as species (Biagi's thresholds), with per-level
  G1, G2 pairs, K76/K76r, KPS, PRAD, KSTEP_S, the excimer formation
  AK8_<level> at Marchal et al.'s k32/k31 (1s3 on k32), MIX_* (Turner,
  Audouard 1991) with reverses, RES_1s4/1s2 (band 0). 126 rows, 10
  species (v0_3: 85 rows).
- **Redissociation by detailed balance** (`thermo.py`; AK8r_*, AK8fr_1s4,
  K15r). Checked against the classical configuration integral (exact in
  the narrow-well limit; 5-6 % at 2000-3000 K with the real quantum).
- **Bands**: kr2_form_1s4, kr2_form_1s3, kr2_fast, kr2_quench, kr2_d0,
  dr_tg, res_escape, kps_5s, k_5s5p_j1; `network.OFF_BY_DEFAULT` lists the
  measured-but-off ones apart in the report.
- **`tail_field.py` and `data/kr_tail_field.txt`** (md5 57e88a86...):
  the field-driven tail table, every point from a cold start;
  Params.tail ('field' default, 'beam' v0_3's).
- **`Params.kin_afterglow_memory`**, **`swarm.el_yield`** (Mano et al.
  2022), `run_v0_4.py --tail/--afterglow-memory`.
- **docs**: `accuracy_overview_v0_4.tex/.pdf` (19 pp) and
  `system_map_v0_4.tex/.pdf`, `make_overview_v0_4.py`,
  `make_system_map_v0_4.py`, `kr_docs_v0_4.py` (174 runs: the undo table,
  the scan's corners on the kinetic EEDF, the EL seeds, v0_4 on Biagi's
  file, the tail ratio along every kinetic run).

### Found (the documents' numbers)
- **The reference point's 147-nm share: 31.22 % on the kinetic
  EEDF** (v0_3 14.05 %), 32.92 % on the two-temperature
  closure (v0_3 15.78 %). With every level at Ueda's 1.7e-32:
  11.92 %; the 1s4's own rate is worth the difference to
  23.05 %; the J = 1 levels pumped by no electron:
  54.81 %; no redissociation: 31.68 %.
  5.8 % of the excimers formed redissociate; 52.6 %
  form out of 1s4.
- **Uncertainty**: +-20.1 points in quadrature
  (+-64 %); the largest terms the J = 1 levels' pumping
  (11.79), Kr2* + Kr2* (9.66),
  e + Kr* -> Kr+, the stepwise ionization of the 5s levels (7.30).
- **Mano's EL yield**: 0.99 to 1.04 of the fit at 6-13.5 Td
  (bands' defaults; 0.93-0.98 with Kr2* + Kr at its top);
  0.89 at 5 Td, 0.65 at 4 Td.
- **Beam**: Morozov +29 % (v0_3 +29 %),
  Eckstrom at 1 atm +1 % (v0_3 -26 %).
- **The closure**: +5 % from the kinetic EEDF at the
  reference point on the field table, +8 % on v0_3's;
  within 5 % (field) and 9 % (beam) over the
  scan's corners.
- **Numerics**: ledger 7e-15; convergence 2e-5;
  the afterglow's time step moves the share by 2e-5; restarts
  11 times, in one run; fallbacks 6 times, in one run.

### Fixed
- The kinetic EEDF's stop in hot gas's afterglow (the fallback; Known-bad).
- The 451-K statement (200-400 K), and the overview's "first nanosecond"
  (the pulse has been time-stepped since 2026-10-01; v0_4 bounds the
  afterglow's steady shape).

### Verified
- `/usr/bin/python3 -m pytest -p no:cacheprovider LIGHT_Kr/v0_4/tests`
  from Kinetic_Modeling: 403 passed, 1 skipped without test_docs.py (29 min, before the fallback); after the last code changes test_network + test_xsec + test_v0_4 107 passed and test_kinetic + test_jung2006 36 passed, 1 skipped.
- Both documents built, read page by page; `test_docs.py`
  7 passed with LIGHT_KINETIC_SLOW=1 (10 min; both generators' --check reproduce their numbers files).
- The overview's generator raises when a claim its text makes fails; on
  this run every such claim held (the J = 1 pumping leads the uncertainty,
  three of the four largest terms are unmeasured, the formation rates are
  most of the change, v0_4 meets Eckstrom at 1 atm, the field table is
  closer at the reference point and no worse at the corners, its tail
  lies below the kinetic one at 2 atm and above at 0.25 atm, the swarm
  meets Mano's fit to 6 %).
- Both PDFs built (overview 19 pp, map 1 p; no overfull box but the map's
  1.6-pt picture width, inherited from v0_3's) and read page by page;
  the reading found and fixed a false "few percent" (the fit's yield is
  0.22 and 0.12 of its 13.5-Td value at 5 and 4 Td), a misdated MIX rate
  (Turner's 1.3e-13 is at 295 K), an unsupported "the beam table keeps the
  tail hot" (it runs 0.67-1.68 of the kinetic tail over the corners), and
  layout (the map's page overflow and overlaps, the rates table's page).

### Not done (in scope, deferred) and limitations found
- The 5p rung's relaxation (KPS) has no detailed-balance reverse; by the
  overview it would match 1s4's excimer formation at the fill's density
  near 1213 K. Below that, and at the reference point, it is negligible;
  the scan's 0.1- and 0.25-atm columns pass it. Documented, not added.
- The field tail table has no excited-population axis (the residual
  above). The beam table is kept for comparison only.
- The J = 1 levels' 5s -> 5p cross sections remain proxies; their band
  k_5s5p_j1 (0 to 1) is now the largest uncertainty term.

## 2026-10-15 (c): v0_4 scope, the owner's choice

Shown the groundwork below, the owner chose for v0_4 (2026-10-15):
- **the 5s levels resolved**, each forming Kr2* at its own measured rate
  (Marchal et al.: 3P2 3.81e-32, 3P1 9.24e-32 cm6/s), Biagi's level cross
  sections setting the split, measured mixing, bands where nothing is
  measured (recommended over keeping one Kr* with 3.81e-32 or Ueda's
  1.7e-32);
- **redissociation** of Kr2* and Kr2+ by detailed balance on the measured
  binding energies;
- **excimer and ion bands**: the two Kr2* states, Kr2* quenching by Kr,
  the gas-temperature dependence of Kr2+ recombination;
- **numerics, a new score and fixes**: a field-driven pure-Kr tail table
  for the two-temperature closure, Mano et al.'s electroluminescence yield
  as an independent score, the afterglow bound, and the 451-K and
  first-nanosecond statements corrected.

### Progress (2026-10-15, evening; the entry is finished when the documents are)
- Built in `v0_4/`: the four 5s levels (species, G1/G2/K76/KPS/PRAD/KSTEP
  per level, MIX, RES), detailed-balance reverses (`thermo.py`, AK8r_*,
  K15r), the excimer and ion bands, `tail_field.py` and its pinned table
  `data/kr_tail_field.txt` (Params.tail = 'field' default, 'beam' v0_3's),
  `Params.kin_afterglow_memory`, `swarm.el_yield` (Mano et al. 2022), and
  a solver fallback (below). Full suite 403 passed, 1 skipped (29 min,
  /usr/bin/python3) before the fallback; the fallback's own test passes.
- Found: the field table's r is 16 % below the kinetic run's own along
  the reference run (0.82-0.86 of it; v0_3's beam table 0.95-1.30). The
  run's 5s and 5p populations, through the 5s <-> 5p collisions, heat the
  tail (r 0.164 with them, 0.130 without, at 500 ns; Penning electrons and
  recombination move it by under 1 %). With the run's own r the closure's
  147-nm share would be 33.36 % (scratch run): the field table's 32.92 %
  against the kinetic 31.22 % is partly compensation by the Maxwellian
  bulk. Documented, not tuned.
- Fixed: the kinetic EEDF stopped the run at the scan's hottest corner
  (0.25 atm, 500 kW/cm3, the gas at 3970 K): a table point 52 ns into the
  afterglow, below the no-field mean, told to leave the growth term out
  (its centre had none), found no steady solution on its centre's grid,
  cold restart included. Root cause: Penning ionization made electrons
  4.6 times faster than recombination took them, so the no-growth problem
  has no steady solution there (a change in the fifth decimal of the
  Penning energies moves it across its convergence edge). `Solver.eedf`
  now falls back to the growth term on the same grid, then to a grid of
  its own (from the solver's first grid), counted in `Solver.fallbacks`
  and the kinetic log. `test_a_table_point_without_a_steady_solution_falls_back`.
- Corrected from sources (DOIs verified against Crossref): Janssens et
  al. 1987 is a quasiclassical theory compared with alpha-excitation
  experiments; its 1.7e5 /s is the pressure-independent term of their
  3P1 loss rate with diffusion and resonance radiation (the band
  res_escape's cite reworded), and its 3P1 three-body term, 82.0 p^2
  (1.36e-31 cm6/s at 295 K), is a calculation above kr2_form_1s4's band
  (cited, the band not widened); its 3P2 term, 3.85e-32, agrees with
  Marchal's 3.81e-32. Turner's 1s4 -> 1s5 law now passes through both his
  points: the rows carry its 300-K value, 1.272e-13, on (Tg/300)^-1.29.

## 2026-10-15 (b): v0_4 groundwork, section 7's limitations surveyed (no code change)

The owner asked: "Begin approaching the limitations of v0_3 such as from
section 7 of the accuracy overview to begin making the improved v0_4".
Each limitation was measured on scratch copies of v0_3 or searched in the
literature (four searches; every value below that would change a default
was checked against its source's abstract or open full text, the rest are
marked). Nothing in v0_3 changed. What v0_4 builds is the owner's choice.

### Found, measured on v0_3 (scratch copies)
- **The afterglow EEDF's steady shape does not matter.** v0_3 already
  time-steps the EEDF through the pulse (`kin_memory_s` = inf); with the
  step kept on after the pulse the reference point's 147-nm share moves by
  +1.5e-5 of itself (14.049435 -> 14.049639 %), and Morozov's low-pump runs
  by 3e-7 and 2e-6. The overview's "the first nanosecond" is outdated.
- **The two-temperature tail is too hot.** The tail ratio r = T_t/Te that
  reproduces the kinetic EEDF's 5s energy loss along the reference run is
  0.12-0.16 over the pulse body; the table (ML's beam-driven KrF states,
  F2 = 0 column) gives 0.18-0.20 there. The closure's 147-nm share is +12 %
  over the kinetic one (15.78 vs 14.05 %), which the pressure-power scan
  (two-temperature) inherits.
- **The formation default decides the answer.** Reference point, kinetic
  EEDF: 14.05 % at Ueda's 1.7e-32 (v0_3), 19.64 % at Tracy & Oskam's
  2.59e-32, 26.45 % at Marchal's 3.81e-32 (two-temperature 15.78, 21.71,
  28.64 %).
- **Redissociation of Kr2* is not negligible.** By detailed balance on the
  measured well depth and the formation rate, along the reference run's own
  history (Tg to 800 K, no feedback), 10-26 % as many Kr2* would break up
  as radiate (Ueda's rate) and 22-58 % (Marchal's); in the pulse a freed
  Kr* re-enters the 5s <-> 5p cycle.
- **The swarm-regime field test passes.** v0_3's swarm (one photon per
  excitation) against Mano et al.'s krypton electroluminescence yield (fit
  0.113 E/N - 0.30, 1e-17 photons cm2 per electron and atom, 1.1 bar):
  0.99-1.04 of it at 6-13.5 Td (Biagi's file 1.02-1.08), below near the
  2.7-Td threshold (0.30-0.89 at 3-5 Td). The training never saw it.

### Found, literature
- **Kr* + 2Kr.** Direct 1s5 values: 1.7 +- 0.2e-32 (Ueda 1984, e-beam,
  0.5-1 atm), 2.59 (Tracy & Oskam 1976), 3.81 +- 0.08 (Marchal et al., TALIF,
  pure Kr: IntechOpen chapter doi:10.5772/intechopen.71942, Table 2, read;
  no temperature stated), 4 (Turner, Phys. Rev. 158, 121 (1967); 6.9 at
  196 K), 4.1-4.2 (Smith & Turner 1963). The resonant 3P1: 9.24 +- 0.80e-32
  (Marchal). Three of the band's values do not measure atomic Kr(5s) + 2Kr:
  Salamero 1979 assigns its p^2 term to a molecular intermediate, Hughes
  1977 is secondhand only, Bouciqué & Mortier 1970 rest on a 1.7-us dimer
  lifetime later measurements contradict. No pressure fall-off from 0.1 to
  15 atm (Audouard et al., JCP 89, 6176 (1988)). At 1 atm 3P1 forms excimers
  17-26x faster than it relaxes to 3P2, so a lumped rate depends on how the
  excitation splits over the 5s levels.
- **Its temperature law** is measured over 200-400 K (Wieme & Lenaerts,
  JCP 74, 483 (1981), title and abstract), not "to 451 K" as v0_3's
  overview says: the sources table's "300-451 K" was the range of the
  power-law conversion (inherited from LIGHT v0_5's SOURCES.md). Nothing is
  measured above 400 K; the low-temperature data are steeper than -0.64.
- **Kr2* well depth.** D_e(A 3Sigma_u+) = 5807.0(8.0) cm-1 = 0.720 eV
  (LaRocque et al., JCP 84, 6627 (1986), abstract; D_0 = 0.705 eV from
  their constants, extrapolated from v' = 30-38); D_0 = 0.62 +- 0.02 eV
  (Killeen & Eden, JCP 84, 6048 (1986), abstract, Delta G_1/2 = 172 cm-1).
  v0_3's 0.715 eV had no source.
- **Kr2+.** D_0 = 9267.8 +- 2.8 cm-1 = 1.149 eV (Signorell et al., JCP
  114, 9840 (2001), abstract): v0_3's 1.15 confirmed. Recombination at
  Te = Tg = 900 K is 4.2e-6 cm3/s, falling as T^-3/2 to 9.6e-7 at 2700 K
  (Cunningham & Hobson, J. Phys. B 5, 1773 (1972), shock tube, abstract);
  Shiu & Biondi's Te law gives 0.87e-6 at 900 K. Kr3+: thermochemistry
  only (Hiraoka & Mori, JCP 92, 4408 (1990)), no formation or
  recombination rate.
- **Kr2* states.** 0u+ 3.4 +- 0.1 ns and 1u 270 +- 5 ns (Morikawa et al.
  1989); quenching of 1u by Kr: 194(+-84) Torr-1 s-1 (6.0e-15 cm3/s,
  Firestone & Selander, JCP 88, 717 (1988), abstract) against "no evidence
  for collisional steps" (Oka et al., JCP 61, 4740 (1974)). No mixing rate
  between the two states, by Kr or by electrons, is measured.
- **Resonance lines.** Six measurements of the 123.6- and 116.5-nm
  oscillator strengths (0.155-0.205, 0.139-0.193); trapped escape ~1.7e5
  /s against 3P1 excimer formation of 8e7-3e8 /s at 1-2 atm: negligible
  in cold gas.
- **The rows with no measurement stay so.** Stepwise ionization out of
  Kr(5s) or Kr(5p): none (one 1977 attempt failed, AFAPL-TR-77-84).
  Kr2* + Kr2*: no krypton measurement (Xe2*: 8e-11, Zamir et al., APL 27,
  56 (1975)); every Kr value traces to a model. Electron collisions with
  Kr2*: one calculation (Flannery et al. 1981). Kr* + Kr*: the Sofia group
  measured Kr 3P2-3P2 ionization in an afterglow (Blagoev et al., J. Phys.
  Conf. Ser. 44, 80 (2006), abstract), the number behind IOP's bot check
  (not bypassed); the same group's superelastic de-excitation by thermal
  electrons, 3P2 (1.5 +- 0.2)e-10, 3P1 5.2e-10, 3P0 9.8e-10 cm3/s
  (J. Phys. B 17, 435 (1984), abstract).
- **Field-driven scores.** Besides Mano et al. (Phys. Lett. B 824, 136782
  (2022), CC-BY), no pure-Kr measurement of VUV efficiency or 5s/5p
  densities at known E/N and high n_e exists: the 5s <-> 5p cycle stays
  unscored. On the Morozov gap: Saito's N_ex/N_i = 0.55 is the only
  measured one; other indications (Eckstrom's ~46 % limit, a preliminary
  W_sc) point to ~0.3; v0_3's beam gives 0.68.
- A search suggested capping the Kr2* + Kr2* band with Eckstrom's
  efficiency. That is anchor-steering and is not done.

## 2026-10-15: v0_3, the cross sections trained on krypton's swarm data by ML, and its documents

The owner asked (2026-10-14): "Let us try to improve like we did for
LIGHT/, we can also use ML to train it to improve the accuracy". LIGHT's
channels were fitted to their own measurements in `ML/`. In the same way,
v0_3's electron cross sections are trained on krypton's swarm
measurements, which measure the cross sections themselves. The owner
approved accepting LXCat's terms of use to read its krypton swarm data.
Device and beam data never train. The low-field ionization coefficient
stays out of the training by this file's Known-bad rule.

### Found
- **LXCat's LAPLACE database** holds 18 krypton swarm datasets LIGHT_Kr
  never used. With Dutton's tables that makes 28 datasets, 441 points
  (`ML/ml/kr_data.py`).
  - Against the most precise of them (Hunter 1988), Biagi's file already
    gives the drift velocity to 3.2 % rms. v0_2's 10.8 % came from
    scoring against Pack 1962 alone, where W climbs out of the Ramsauer
    minimum.
  - Jacques 1986 and Bhattacharya 1979 agree with v0_2's alpha/N to a
    few percent at 40-100 Td.
- **The two-term swarm is not what the experiments measure** (LoKI-MC on
  Biagi's file, `v0_3/swarm_mc.py`).
  - The bulk drift velocity is up to 23 % above the two-term flux one at
    300 Td.
  - The SST coefficient is 9 % below the two-term k_i/W at 300 Td.
  - The two-term transverse diffusion is 46-79 % high from 10 to 300 Td.
  - The two-term flux W matches the Monte Carlo's to 0.3 % from 10 to
    300 Td (1.8 % at 0.2 Td, absorbed by the trained factors), and k_i/W
    lies 1.4-3.9 % below the Monte Carlo's flux coefficient.
  - v0_2's overview compared k_i/W with Townsend data, which put it closer
    to them at high field than it was.
- **The train set** (the joint-fit consistency check, chi2 per point at
  most 3): Hunter 1988, Nakamura 1990 and Pack 1962 at 368 K (W), and
  Bhattacharya, Dall'Armi, Jacques, Kruithof and Heylen at 50-300 Td
  (alpha/N).
  - Dropped, in order: English 1953, Kucukarpaci 1981, Pack 1992 at
    300 K, Dall'Armi's W, Pack 1992 at 195 K, Pack 1962 at 300 K,
    Pack 1962 at 195 K (3.6), Bowe 1960 (3.5).
  - Threshold 2 keeps the same set. At 4, Pack 1962 at 195 K and Bowe
    1960 join, and the factors move to m_6 1.31 and x_1s 1.16.
  - A leave-one-out criterion, tried first, dropped Hunter 1988
    (Known-bad).
- **Which factors train: seven of nine, chosen on held-out data.** With
  all nine free (the first training) the fit matched the train sets best
  (pooled chi2 0.49 per point). But it traded the momentum transfer at
  0.03 eV (1.13) against the Ramsauer minimum at 0.6 eV (0.88), and the
  held-out data that weigh those energies differently got worse:
  Koizumi's D_T/mu went from chi2 2.5 on Biagi's file to 8.1, and Pack's
  195-K drift velocities from 3.1 to 5.5 (1962) and from 6.1 to 8.7
  (1992). Five variants were compared on the true solver
  (`ML/ml/kr_train.free_variants`). Holding the knots at 0.03 and 0.2 eV
  at 1 fits those three best of the five (0.84, 2.77, 5.27), with the
  train sets at 0.64 pooled and Hunter at 0.70 (0.13 with all nine).
  v0_3 trains that variant (`ML/ml/kr_swarm.FREE`). Those held-out sets
  made the choice, so the documents call them validation sets. The
  beam's partition and the low-field alpha/N saw neither the training nor
  the choice.
- **The factors** (least squares; MCMC and the neural estimator agree,
  `ML/README.md`).
  - Momentum transfer: 1 and 1 (not trained), then 1.157, 0.926, 1.039
    and 1.227 at 0.03, 0.2, 0.6, 2, 6 and 20+ eV.
  - Excitation: 0.886 (5s), 0.925 (5p), 0.753 (higher levels).
  - The first emulator box (|ln| < 0.45) put the 0.6-eV factor (then
    free with the two below it) on its edge. Widened to 0.6 and
    resimulated, it came back inside.
- **Scores** (as measured, Biagi's file -> trained):
  - train W 4.8 -> 2.8 % rms (2.1 % with all nine free); train alpha/N
    14.0 -> 9.9 % rms (mean -11 -> 0 %).
  - validation: Koizumi's D_T/mu 15.8 -> 7.1 % rms; Pack's 195-K drift
    velocities chi2 3.07 -> 2.77 (1962) and 6.12 -> 5.27 (1992).
  - held out: five of the eight dropped drift velocities move closer and
    three further (Pack 1992 at 300 K 14.2 -> 17.7, Dall'Armi 10.3 ->
    13.4, Kucukarpaci 39.6 -> 48.4 per point).
  - the untouched tests: alpha/N below 50 Td, mean -44 -> -36 %; the
    beam's partition against Saito et al., W +16 -> +10 % and N_ex/N_i
    +44 -> +24 %.
  - LoKI-MC on the trained set gives the bulk-W and SST ratios within
    1.65 % of Biagi's file's from 30 to 300 Td (the most: alpha at 50 Td,
    counting error about 1 %). A new 0.5-Td point (1e4 electrons) gives
    the D_T ratio 0.9052, which the table's log-linear reading between
    0.2 and 1 Td had to 5e-5.
- **The reference point** (kinetic / two-temperature): 147 nm 14.05 /
  15.78 % (v0_2 14.20 / 16.00 %). Propagated: +-11.8 points in quadrature
  (+-84 %), and the trained cross sections' posterior +-0.04 points
  (+-0.16 on the two-temperature closure). Under the beam: Morozov's low
  pump +29 % (v0_2 +31 %), Eckstrom's high pump -26 % at 1 atm (v0_2
  -25 %).
- **v0_5/ and run_v0_5.py are no longer in LIGHT_Kr** (found absent on
  2026-10-15, not removed by this work). The README no longer lists
  them, and v0_1's tests that compare with the copy skip.
- The kinetic grid search: on the trained set, as on Biagi's file, the
  60-, 40- and 30-eV starts land on one 21.3-eV grid with the same
  ionization rate (`test_grid_search_from_a_grid_five_times_too_long`).
  The nine-factor trial set had put the 30-eV start on a 19.3-eV grid,
  1.2 % off in that rate.

### Added
- `v0_3/` from `v0_2/`, and `run_v0_3.py`.
- `xsec.py`: `raw_processes` (Biagi's file), `Factors`, `apply_factors`,
  `exc_group`, `read_factors`, `trained`, `using` and `UNIT`, and
  `TRAINED_FILE = data/kr_swarm_trained.json`, md5-pinned. `processes`
  is the file times the factors in force.
- `swarm.py`: N D_T and D_T/mu.
- `swarm_mc.py`: LoKI-MC's ratios (RATIOS, DT_LOW, FLUX, TRAINED_CHECK)
  and `as_measured`.
- `data/kr_swarm_trained.json`: the factors and the training's record
  (`python3 ML/run.py kr-train`), including the variant comparison.
- `tests/test_v0_3.py` (the gate). v0_2's tests are carried over and
  point at v0_3.
- `docs/`: `accuracy_overview_v0_3.pdf` and `system_map_v0_3.pdf`, built
  by `make_overview_v0_3.py`, `make_system_map_v0_3.py` and
  `kr_docs_v0_3.py`.
  - New job kind "xs": the model on other factors, here the posterior
    draws.
  - The runs' cache keys on the files a run reads. The numbers files are
    stamped with `full_md5`.
  - The overview scores the five variants live on the solver (Table
    "Which factors train") and checks the record against them.
  - New references, each checked on Crossref: the 11 swarm datasets'
    papers, Goodman & Weare 2010, McKay et al. 1979, Kingma & Ba 2014,
    Stokes et al. 2020.

### Fixed
- **A kinetic beam run stopped: "kinetic EEDF: no solution with or
  without the growth term at mean 0.385 eV".**
  - Symptom: the documents' e-beam run at 0.5 atm and 3e5 W/cm3 raised
    in its afterglow (n_e 1.06e16 cm-3, Tg 1196 K) on the seven-factor
    set. The nine-factor trial set had run through.
  - Root cause: the EEDF table's point 16 % below the centre is solved
    without the growth term (the centre had none). It was warm-started
    from the last table's point, which had the term on (nu -5e6 /s).
    From that F the solve without the term converged nowhere, though the
    state has a solution: the no-field one, mean 0.449 eV, since the
    target lay below it. From a cold start, from this table's centre or
    its other points, it converges, to the same EEDF to 2e-7 in ln F
    (the failing state replayed from a capture, 2026-10-15).
  - Fix: `boltzmann.Solver.eedf` starts a grid's solve again cold, once,
    where the warm start finds no solution, and counts it
    (`Solver.restarts`). The model logs it per table update, the
    documents' runs sum it, and the overview reports it. Only the
    failure path changed: a run that never failed is the same number for
    number. `test_a_failed_warm_start_restarts_cold` holds it.
  - LIGHT v0_5 and LIGHT_ArF v0_4 carry the same table code without the
    restart; they were not changed here.
- **The documents' runner stalled 20 minutes with nothing in the log**
  after that failure. The exception stopped the main loop while the
  pool's jobs in flight ran on, and the traceback appeared only after
  them. `kr_docs_v0_3.run_many` now names the failed job at once and
  cancels the jobs not started.

### Changed (v0_3 against v0_2)
- The electron cross sections, and the solver's cold restart (above).
  Rows, bands, species, jung2006 and the tail table are v0_2's byte for
  byte, and with unit factors the swarm and a model run are v0_2's
  number for number.
- Tests pinned to v0_2's numbers now hold v0_3's, with v0_2's in their
  comments (`test_model.py`, `test_jung2006.py`). The tests of the solver
  against ML's dense solver and LoKI-MC run on Biagi's file
  (`xsec.using(xsec.UNIT)`).

### Verified
- `LIGHT_Kr/v0_3/tests` on the final code: 387 passed, 1 skipped (19 min,
  test_docs apart); `test_docs.py` with LIGHT_KINETIC_SLOW=1: 7 passed (8 min; both
  generators reproduce their files exactly).
- `ML/tests`: 53 passed, 4 skipped, 1 xfailed. `ML_KR_SLOW=1
  ML/tests/test_kr_slow.py`: 4 passed.
- The documents: 120 runs on the final code (36 min on 8 workers). Both
  PDFs build without a warning (the overview 15 pages, the map 1 page)
  and were read page by page. Every claim the generator checks holds:
  the training improves the train sets and the low-field alpha/N, v0_3's
  variant fits the validation sets best, the beam's partition moves
  toward Saito's, and the five variants' live scores agree with the
  training record to 1 % for every W and alpha set. The solver restarted
  once, in one run (the 0.5-atm beam).
- The failing beam state, captured and replayed: the patched solver
  returns the cold start's EEDF (mean 0.449366 eV) after one restart.
- Not run to completion: LoKI-MC at 0.05 Td (2e4 electrons) and the
  1e5-electron repeats at 0.5, 1 and 3 Td (stopped after about 5 hours).
  `DT_LOW` holds the 2e4 and 1e4 points.

## 2026-10-14 (b): v0_2, krypton's electron collisions level by level, and its documents

The owner asked: "LIGHT_Kr appears to have an error of -54% in some cases,
let us proceed with making v0_2 to improve this". The -54 % is v0_1's
accuracy overview's ionization coefficient below 50 Td (against Kruithof
1940 -73 %, Heylen 1971 -33 %). Shown the diagnosis below and asked which
v0_2 to build, the owner chose the level-resolved cross sections with the
Hornbeck-Molnar ionization as a band, and both documents. No rate was set
or moved by any score.

### Found (the -54 %)
- Not the solver. LoKI-MC (the vendored Monte Carlo solver, no two-term
  approximation) on v0_1's Biagi-v7.1 file: the two-term alpha/N is 10,
  6, 3 and 6 % below it at 20, 30, 40, 50 Td; on v0_2's 8.97 file +3.7,
  -4.0, -3.5, -2.5 %, inside the Monte Carlo's counting error (19 % at 20
  Td, 1 % at 50 Td); mean energies within 0.15 %. ML's dense solver
  matches v0_2's sparse one to round-off on v0_2's file.
- Every published Kr set falls below these low-field measurements:
  Bordage et al., J. Phys. D 46, 334003 (2013), "a general trend",
  ascribed to Penning ionization of impurities in the experiments;
  Jacques et al. 1986 measured below Kruithof at the low end too.
  Kruithof and Heylen differ by up to 2.2x near 30 Td.
- Resolving Biagi's levels moves the mean below 50 Td from -54 to -47 %
  (Kruithof -73 -> -68, Heylen -33 -> -23 %), 50-100 Td from -19 to -14 %,
  above 300 Td from +12 to +9 %; 100-300 Td 6.5 -> 7.6 % rms, the drift
  velocity 9.6 -> 10.8 % rms.
- Hornbeck-Molnar ionization of the 20 levels at or above E(Kr2+) =
  12.85 eV (Hornbeck & Molnar 1951 saw Kr2+ appear 0.7 (+0.7 -0.3) V
  below Kr+): at the band's top alpha/N rises 16.8x at 14 Td, 3.8x at 30,
  1.7x at 100 Td; 18 of the 21 points below 50 Td lie between the band's
  ends (Kruithof's three at and below 16.95 Td above even its top), but
  the top puts the model +39 % above the 100-300 Td data. No krypton
  measurement gives the share, so it is 0 by default.
- The answer does not move: at the reference point direct ionization
  makes 0.012 % of the ions (Penning the rest), and hm_kr at its top moves
  the 147-nm share by -0.02 % (kinetic), -0.01 % (two-temperature).

### Added
- `v0_2/` from a copy of `v0_1/`, and `run_v0_2.py`; `v0_2/README.md` and
  `v0_2/SOURCES.md` (section 7 new, each source with how it was read).
- `v0_2/data/Kr_Biagi897_LXCat.txt`: LXCat's 'Biagi' Kr file (Magboltz
  8.97, 51 levels), the owner's download of 2026-06-18
  (`old_files/Kr_LXCat_new.txt`, Kinetic_Modeling commit 36f70d3), byte
  for byte and md5-pinned (test_xsec checks it against the repo copy).
  v0_1's Biagi-v7.1 file is not carried.
- `network.py`: every level a row at its own threshold, `G1_<1s>` into
  Kr*, `C1a_<2p>`, `C1b_*` (4d, 6s) and `C1c_*` (12.76 eV up) into Kr**,
  the excess over the product's energy to heat (`level_tag`,
  `level_group`, `levels`, `group_rows`); twenty `HM_<lvl>` twins (rate
  kind 'assoc', first order in Kr: the partner is in the share) on the
  band `hm_kr` (0 default, 1 top); SUPER on the four 5s levels' summed
  cross section at the 1s5 energy (`xsec.lumped_5s`, outside processes()
  so the momentum transfer counts no level twice). 85 rows (87 with the
  optional group).
- `boltzmann.py`: collision kind 'assoc' (scattered at e - u, a second
  electron born at rest, counted in the growth term).
- `swarm.py`: the solver as a swarm with the Hornbeck-Molnar share.
- `model.py`: the 'assoc' kind; rows at factor 0 left out of the kinetic
  equation; `ivp_options` gives BDF/Radau the Jacobian's structure (only
  densities, U_e and Tg feed the RHS: 9 right-hand sides a Jacobian, not
  101) and the two-temperature closure computes its shared tail parts once
  per state (`_tail_memo`, `xsec._two_t_common`): a two-temperature run
  3.5 s (36 s without; v0_1 2.3 s), the same numbers to 6e-9.
- Tests: `test_v0_2.py` (the gate: v0_1's rows, bands, tables and jung2006
  data carried over; HM rows inert at 0 number for number; hm_kr = 1
  ionizes and closes on both closures; the 'assoc' collision's growth and
  energy balance; the swarm against LoKI-MC and the band). v0_1's tests
  carried to the levels; ML's dense-solver comparison now on v0_2's own
  file (rewritten for its parser in a tmp dir).
- `docs/`: `accuracy_overview_v0_2.pdf` (12 pages) and
  `system_map_v0_2.pdf` from `make_overview_v0_2.py`,
  `make_system_map_v0_2.py` and `kr_docs_v0_2.py` (68 runs, both closures,
  cache `docs/cache/runs_v0_2.json` keyed by v0_2's code md5; 40 min on 7
  workers at load 35-55). The overview's swarm section is new: v0_1's and
  v0_2's scores side by side, the Monte Carlo check, the literature, the
  band. Five new references, each checked on Crossref: Biagi 1999
  (Magboltz), Bordage et al. 2013, Dias et al. 2023 (LoKI-MC), Hornbeck &
  Molnar 1951, Jacques et al. 1986. `v0_2/tests/test_docs.py`.

### Changed (v0_2 against v0_1)
- `species.py`: Kr* 9.9152, Kr** 11.3035, Kr+ 13.9996 eV (8.97's
  thresholds, NIST's levels); Kr2+ = IP - 1.15; the 5s -> 5p gap 1.388 eV.
- `jung2006.py`: its docstring gives v0_2's energies; its data are v0_1's
  (test_v0_2 compares them).
- The reference point (1 atm, 300 kW/cm3 x 1 us), kinetic / two-
  temperature: 147 nm 14.20 / 16.00 % (v0_1 14.06 / 15.05); near IR 0.64 /
  0.63; reaction heat 83.11 / 81.40; elastic 2.58 / 2.49; Tg 807 / 797 K;
  Te pulse mean 2.24 / 2.00 eV (2.29 / 2.26). The 5s -> 5p cycle: 34 / 29
  pumpings per Kr2*, 77.6 / 75.2 % of the input to heat. Kinetic closure:
  381 solutions, 1.32-15.5 Td, field/input 1.0001-1.0002, balances to
  9e-6, ledger 3e-15.
- The overview's other scores: Morozov 2008 0.550 (+31 %, 2.6 sigma, at
  v0_2's ceiling 0.548); Eckstrom 1988 at 3e5 W/cm3: 0.113 / 0.345 /
  0.529 at 0.5 / 1 / 2 atm (79 / 37 / 4 % below the ceiling); the beam's
  W_s 15.38 eV (Saito 15.4, -0.1 %); propagated +-11.9 points (+-84 %) on
  the kinetic 147-nm share, Kr* + 2Kr 75 % of the variance.

### Verified
- `v0_2/tests` from Kinetic_Modeling on the final files: 389 passed, 3
  skipped (the slow reference-point test and test_docs' two regeneration
  tests, LIGHT_KINETIC_SLOW=1) in 14 min. `test_docs.py` with
  LIGHT_KINETIC_SLOW=1: 7 passed (both generators reproduce their files
  byte for byte). v0_1's `test_docs.py`: 5 passed, 2 skipped (its
  documents still match v0_1's code; nothing under v0_1/ changed).
- Both PDFs compiled with latexmk (one 1.6-pt overfull line in the map's
  footer) and read page by page; every claim the generators make about
  their numbers is asserted before they write.
- The Jacobian's structure changes no result: the kinetic reference run
  with the dense finite-difference Jacobian gives 14.200479608069 % at 147
  nm against 14.200479608068 % (two-temperature: 6e-9). The kinetic
  closure's convergence: half the renewal step 1e-5, a renewal every 5 ns
  3e-5, twice the cells 3e-5, rtol 1e-9 2e-7 (the same 381 renewals);
  v0_1's rtol change moved it 3e-4.

## 2026-10-14 (a): K76 and K76r take one row_scale factor

The owner asked: "Fix the bug under all LIGHT directories LIGHT/
LIGHT_ArF/" (entry (c)'s row_scale defect). Porting the fix found one
more case of it here: under the kinetic EEDF the 5s <-> 5p reverse
(K76r) comes from K76's scaled cross section, so a row_scale on K76r
alone reached neither direction's kinetic rate, and one on K76 alone
scaled both directions there but the forward alone on the two-temperature
closure.

### Fixed
- `v0_1/model.py`: `Model.__init__` refuses different row_scale factors
  on K76 and K76r, so a factor on the pair means the same on both
  closures (both directions scaled).
- `v0_5/model.py` (the KrF reference v0_1's parity tests run against):
  LIGHT/v0_5's fix applied, so the copy stays LIGHT's v0_5 byte for byte
  (`cmp`), with LIGHT's test added to `v0_5/tests/test_v0_4.py` (it
  failed on the copy before the fix).
- `v0_1/tests/test_rung_walls_wave.py`:
  `test_row_scale_reaches_the_kinetic_table` also holds K76's and K76r's
  kinetic rates at 1.1 times the default's with both scaled by 1.1, and
  refuses K76 alone, K76r alone and unequal factors.

### Changed
- `docs/`: both numbers files regenerated on the new code (its md5 is in
  them), the cache's 66 runs rerun (30 min on 6 workers), both PDFs
  rebuilt. The 66 runs' 7434 numbers are identical to the old cache's;
  the numbers files differ only in the code stamp and in the kinetic
  run's wall-clock minutes (3 -> 4, a loaded machine).

### Verified
- `v0_1/tests` from Kinetic_Modeling: 345 passed, 4 skipped, 1 xfailed,
  with `test_numbers_come_from_the_current_code` failing on the old stamp
  until the documents were regenerated; then `test_docs.py` with
  LIGHT_KINETIC_SLOW=1: 7 passed, and the full suite 346 passed, 4
  skipped, 1 xfailed (15 min).
- `v0_5/tests` on the fixed copy: 343 passed, 2 skipped (49 min; 342 and
  2 before, the new test the one more); v0_1's `test_v0_1.py` against it:
  19 passed, 1 skipped (the slow run).

## 2026-10-13 (c): the accuracy overview and system map of v0_1

The owner asked: "Then make the accuracy overview for v0_1, as well as the
system map of the pure Kr model". Both are in `docs/`, generated from the
model, and stand alone: no local file is named, and every source cited
has a DOI (each checked on Crossref, authors included). No rate was set
or moved by any score.

### Added
- `docs/system_map_v0_1.pdf` (one landscape page) from
  `system_map_v0_1.tex` and `make_system_map_v0_1.py`: pure krypton at
  LIGHT's reference point on the kinetic EEDF, every flow a percent of
  E_dep, LIGHT's palette and rules. 18 nodes are checked to conserve
  their flows (to 1e-8; they hold to 1e-13) before a number is written,
  every row is drawn or printed (PEN_KR2, KSTEP_S, KSTEP_P at 0), and the
  5s <-> 5p cycle is drawn as its two net bands with KPS's heat as a band
  to the ledger (78 % of the input). A bands box gives the 147-nm share at
  each end of every band the rows read; a footer gives the cycle, the
  kinetic EEDF's diagnostics and the sources of the printed rates.
- `docs/accuracy_overview_v0_1.pdf` (9 pages) from
  `accuracy_overview_v0_1.tex` and `make_overview_v0_1.py`: the scores,
  the propagated uncertainty, the numerical checks, what v0_1 is, the
  reference point on both closures with the 5s -> 5p scale, Kannari's fit
  and the pressure x power scan, the swarm figure, the e-beam scores, a
  simple energy ladder, the limits, the rates appendix and 35 DOI
  references.
- `docs/kr_docs.py`: the documents' runs (66 jobs, both closures) in
  parallel with a JSON cache (`docs/cache/`) keyed by the md5 of v0_1's
  code and data and by each job's Params, so a model change reruns
  everything; and the e-beam drive the overview scores under, a fixture of
  the documents and not of the model (BeamModel: ion pairs, Kr**, Kr* and
  injected secondaries from a Fowler/Spencer-Fano degradation on v0_1's
  own Biagi cross sections from 536 eV: W 26.9 eV, N_ex/N_i 0.76, 15.31
  eV per excitation or ion).
- `v0_1/tests/test_docs.py`: every printed number generated, both files
  stamped with the current code's md5, every citation a bibitem with a
  DOI and every bibitem cited, no local file or path in either text; with
  LIGHT_KINETIC_SLOW=1 both generators reproduce their files byte for byte.

### Found (the overview's scores; none moved a rate)
- Swarm (the solver as a swarm on Biagi v7.1, Dutton's printed values):
  Pack 1962's drift velocity, 41 points 0.0008-3 Td, rms 9.6 %; Kruithof
  1940 and Heylen 1971's alpha/N rms 6.5 % at 100-300 Td, mean +12 % above,
  -19 % at 50-100 Td and -54 % below 50 Td (Kruithof -73 %, Heylen -33 %).
  The reference point's ions are Penning's, so that does not reach it.
  Half and twice the 1000 cells move every rms by 0.17 points.
- E-beam: low pump (Morozov 2008, 0.42 +- 0.05) v0_1 0.551, +31 %, 2.6
  sigma, at its ceiling of one photon per excitation and per ion, which
  Saito 2003's W_s (15.4 eV, 0.547) supports. At 3e5 W/cm3 (Eckstrom
  1988, near their limit of ~46 %, pressure not in the abstract) v0_1
  loses 79, 37 and 4 % below its ceiling at 0.5, 1 and 2 atm, to Penning
  ionization of Kr* (0.655 photons a pair at 1 atm); kr_pen's low end or
  kr2_form's top cut the 1-atm loss to 18 or 7 %.
- Propagated (7 bands and 3 measured constants, both ends, both
  closures): the reference point's 147-nm share 14.06 % +- 11.9 points in
  quadrature (+-85 %; two-temperature +-78 %), excursions -16.9 / +31.4.
  Kr* + 2Kr leads with 75 % of the variance (its measurements spread
  x3.2: Ueda's 1.7e-32 at 0.5-1 atm is the default, Hughes's 5.5e-32 the
  top gives 34.7 %), then Kr2* + Kr2* and the 5s -> 5p scale (3.2 points
  each). The three unmeasured rows at their tops cut the share by 20-45 %.
- Kinetic closure at the reference point: 381 solutions, E/N 1.34-15.8 Td,
  field/input 1.0001-1.0002, balances to 7e-6; kin_step 0.05, a 5-ns
  renewal and 1200 cells move the 147-nm share by 3e-4, 3e-5 and 2e-4;
  rtol 1e-9 by 3e-4 (the closure's own spread).

### Fixed
- `v0_1/model.py`: `Params.row_scale` reached no row the kinetic table
  carries. Symptom: in the overview's propagation, R11p at +-10 % moved the
  kinetic run's R11p events by 1e-10 and its 147-nm share not at all,
  where the two-temperature closure moved the share by 1 %. Root cause:
  `Model._rate` read the table's rate without the row's `c.factor`, and
  past the table's ends the Maxwellian ratio cancelled it; the kinetic
  processes were built from the unscaled cross sections and recombination
  constant. Only K76 and K76r carried their factor (the band k_5s5p scales
  their cross section). Fix: `_kinetic_setup` builds a cross-section row's
  process on its cross section times the factor, multiplies the
  superelastic row's g ratio and the recombination template's constant by
  it, so the kinetic rate and the collision operator carry it as the
  closure's own rate does. A row_scale on G6 (the elastic exchange, which
  neither closure scales) is now an error. Default runs are unchanged (the
  factors are 1).
- `v0_1/tests/test_rung_walls_wave.py`:
  `test_row_scale_reaches_the_kinetic_table`, written first and failing
  before the fix: the scaled cross section and recombination constant in
  the kinetic processes, R11p's rate at one state 1.1 times the default's
  (to 1e-9), G1's at 0.8 its cross section's rate over the run's own EEDF
  (to 1e-6, the table's centre), the other rows untouched, G6 refused.
  After the fix R11p at +-10 % moves the kinetic 147-nm share by +-1.4 %
  (two-temperature +-1 %).

### Verified
- `v0_1/tests` from Kinetic_Modeling after the fix: 341 passed, 2 skipped,
  1 xfailed (13 min). With `test_docs.py`: 346 passed, 4 skipped (the slow
  runs, LIGHT_KINETIC_SLOW=1) and 1 xfailed (10 min).
- Both generators rerun from the cache reproduce their files exactly
  (`test_docs.py`, LIGHT_KINETIC_SLOW=1, 50 s). Every claim the overview's
  text makes about its numbers (the formation term leads, the unmeasured
  rows cut by a fifth to a half, v0_1 above Morozov, the order of the
  high-pump runs, the small 2-atm loss, Penning ionization making the
  reference point's ions, the sign of the swarm errors) is asserted by the
  generator, which refuses to write when one fails.
- Both PDFs compiled with latexmk, no overfull box in the overview, and
  read page by page.

## 2026-10-13 (b): v0_1, the pure-krypton model

The owner asked: "Make v0_1 of pure Kr from the remaining v0_5 directory".
v0_1 is v0_5 copied and converted: the fluorine is taken out and nothing
else changes but two table entries corrected to their paper. No rate was
set from a device score.

### Added
- `v0_1/`, from a copy of `v0_5/`, and `run_v0_1.py`.
- `v0_1/tests/test_v0_1.py`: the conversion's gate (species, rows and
  bands are v0_5's krypton ones field for field; no fluorine is left;
  jung2006 differs by the two corrections alone; the kinetic table and the
  right-hand side from one state are v0_5's; runs agree with v0_5 at x_f2 =
  0 on the same tables to the solver's tolerance in nine configurations).
- `v0_1/README.md` and `v0_1/SOURCES.md` (v0_5's krypton sources, and a
  section on the corrections).

### Changed (v0_1 against v0_5)
- `species.py`: the six krypton species, values and citations v0_5's.
- `network.py`: v0_5's 18 krypton rows of its 59, in its order; the
  optional group `electron` (KSTEP, SUPER) alone; the bands those rows read
  (kr2_pen, kr2_form, kr_pen, kps_kr, k_5s5p, kstep_s, kstep_p) and
  ib_scale; ION_K0 and LADDER_KANNARI without their fluorine entries;
  build() has no `net` (one network).
- `model.py`: Params loses x_f2, ne_mode, tau_sp_s, trimer_tg_pow and
  network; n_e is always solved (the cations' sum); eedf 'auto' is the
  kinetic EEDF with the exact captured energy; the two-temperature tail is
  read at the table's F2 = 0 column; the rate kinds only fluorine rows used
  (inv_tau, flannery, ionion, balance, dr_sqrt), the Debye screening and
  the kinetic attachment bookkeeping are gone; the kinetic drift triggers
  lose the F2 fraction (an event that never fired).
- `xsec.py` vendors krypton alone; `boltzmann.py` assembles the vendored
  gases (`GASES`, krypton) where it named Kr and F2; `ib.py` reads
  krypton's momentum transfer alone.
- `report.py`: the fluorine lines (248 nm, KrF*, Kr2F*, F atoms, F2
  burn-out, the trimer's law) and the 248-nm optics are gone; Kr2* has a
  flow ("to Kr2*") and a fate line; row_dict carries 147 nm, the near IR,
  the Kr* -> Kr2* and Kr2* -> 147 nm shares and n_e at the end. The bands
  line lists every 0-default band the rows read: v0_5's listed the 'param'
  bands alone and left out kstep_s and kstep_p.
- `run.py`: no --f2-pct, --ne-mode, --tau-sp, --trimer-tg-pow or
  --network; the scan is pressure x power; the tables report 147 nm and
  the near IR.
- The tests: LIGHT's carried to krypton, each assertion measured again on
  pure-Kr states. `test_v0_4.py` and `test_v0_5.py` are
  `test_rung_walls_wave.py` and `test_jung2006.py`. Dropped with the
  fluorine: the KrF(B,C) decays, the trimer, attachment, ion-ion and
  248-nm optics tests. Pure-Kr analogs added: a dilute Kr* seed ends as
  147-nm photons on the two-step chain's analytic course; a Kr** seed
  radiates its analytic share A/(A + k_ps N); recombination capture (0.95
  Te) in place of attachment capture. `test_ib.py` finds LIGHT's v0 kernel
  in `LIGHT/` beside LIGHT_Kr (appended to the path, never shadowing),
  where v0_5's copy skips it.

### Fixed
- `jung2006.py`, two misreadings of Jung et al., PRA 73, 022722 (2006),
  read on the page images and in two text extractions. Symptom: the
  ported table test failed once it compared without pytest's absolute
  tolerance. Root causes: Table II's 1s5 -> 2p6 at 4 eV transcribed as 21
  (printed 23, p.022722-6); Table IV's last row, 1s5 -> 2p5 (p.022722-10;
  the Letter's Table I), booked as a core-changing 1s3 -> 2p5 (0.4 at 3.5
  eV, 0.11 at 6 eV), which the paper never extracts (p.022722-7). Effect:
  the lumped cross section +1.15 % at 4 eV and -0.1 % at 2-4 eV; the
  reference point's 147-nm share 15.066 % -> 15.051 % (two-temperature).
  Both stay in LIGHT/v0_5 (LIGHT's to fix).

### Removed
- From the copy: `optics.py`, `data/F2_SIGLO_LXCat.txt`, and v0_5's
  README, SOURCES and the two renamed test files. Their originals are in
  `v0_5/`.

### Found
- The 5s -> 5p cycle decides pure krypton's answer at LIGHT's drive: the
  electrons pump Kr* to the 5p rung 32 times for each Kr2* formed, and
  krypton's relaxation of the rung returns 76.5 % of the input as heat
  (two-temperature). The band k_5s5p moves the 147-nm share from 20.1 %
  (0.5) to 13.0 % (1.4), 51.6 % with the cycle off. Pressure and power
  move it as much: 59 % at 2 atm and 100 kW/cm3, under 0.2 % at 0.1 atm,
  where a microsecond heats the gas to 2000-8800 K.
- The kinetic closure's run-to-run spread (Known-bad), the S curve near
  0.04 Td and the long-grid search (Known-bad), all of LIGHT's solver,
  unchanged.

### Verified
- `test_v0_1.py`: rows, species and bands field for field; the kinetic
  table from one state bit for bit, the right-hand side to 1e-14; runs at
  rtol 1e-8 within 3e-7 (trajectories) and 6e-7 (row_dict) of v0_5 at
  x_f2 = 0 in the two-temperature, Maxwellian, walls, wave, ladder, k76,
  electron-group and band-top configurations, falling to 3e-9 at rtol
  1e-10, with v0_5's fluorine densities and events exactly 0. Under the
  kinetic EEDF (0.3-us pulse) v0_1 and v0_5 differ by what v0_1 at rtol
  1e-8 and 1e-9 differ by: 3e-4 on the 147-nm share (slow test).
- The kinetic reference point: 147 nm 14.06 %, near IR 0.65 %, reaction
  heat 83.13 %, elastic 2.70 %, Tg 808 K; the field carries the input to
  1e-3, the two energy balances agree to 7e-6, the ledger closes to 4e-15.
  Its steady EEDF solutions all lie in the grid search's window, apart from
  those on the 12-eV minimum grid with cold electrons (25-34 decades).
- From Kinetic_Modeling: `v0_1/tests` 340 passed, 2 skipped (the slow runs)
  and 1 xfailed (the long-grid search) in 15 min (load average near 20).
  With LIGHT_KINETIC_SLOW=1 the two slow tests pass (33 min): v0_1 and v0_5
  under the kinetic EEDF within three times the closure's own rtol spread
  (the 147-nm share within 2e-3), and at the pure-Kr reference point
  kin_step 0.05, kin_max_dt_s 5 ns and kin_cells 1200 move the 147-nm and
  near-IR bands by less than 5e-4 of themselves (LIGHT's bound).
- Under Homebrew's python3 (3.14, numpy 2.5, scipy 1.18) the two-temperature
  reference point gives the same numbers, and a 20-ns kinetic pulse closes
  its ledger to 4e-16.

## 2026-10-13 (a): the KrF copies removed

The owner asked: "I copied LIGHT/ to make LIGHT_Kr/ can you please remove
the old files pertaining KrF, so that we can make a model that is for pure
Kr". The rule applied: remove every KrF-only file that v0_5 does not need
to run. v0_5 is LIGHT's newest model and the pure-Kr model's starting
point. The KrF parts it still runs on stay for the conversion
(`README.md`, "Still KrF inside v0_5").

### Removed
- `v0/` to `v0_4/`: the older 0-D KrF models. v0_5 carries their machinery.
- `1Dv0/`, `1Dv0_2/`, `1Dv0_3/`: the 1-D KrF models on Fluid1D.
- All of `docs/`: the accuracy overviews, system maps and energy ladders
  with their generators, `REVIEW/`, `survey_2026-09-27/`, and `krf_data/`
  (the KrF measurements, the Tahoma and Kleikamp scorers, the digitized
  figures and the source PDFs).
- `IMPROVEMENT_PLAN.md`, `improvements_next.txt` (the 1-D layer's next
  step), LIGHT's `README.md` and `CHANGELOG.md` (replaced by new ones), and
  the launchers `run_v0.py` to `run_v0_4.py` and `run_1Dv0.py` to
  `run_1Dv0_3.py`.
- From `v0_5/`: `two_level.py`, the KrF(B,C) pair in closed form, which
  nothing at run time imports; `tests/test_two_level.py`; and
  `tests/test_docs.py`, the tests of the KrF documents.
- From `v0_5/tests/test_model.py`:
  `test_basic_network_against_hakuta_c_decays`, a strict xfail scoring the
  KrF(C) decays with `two_level.py` against `docs/krf_data/digitized/`,
  and the `pathlib` import only it used.
- From `v0_5/tests/test_kinetic.py`: `test_two_t_closure_reproduces_v0_2`,
  v0_5 against KrF v0_2, which is gone.

The removed files went to the macOS Trash in one folder,
`LIGHT_Kr_removed_KrF_2026-10-13`, with their paths kept.

### Kept
- `v0_5/`, unchanged apart from the removals above, and `run_v0_5.py`.
- `improvement_for_EEDF.txt`, on the kinetic EEDF's steady-shape
  assumption, which v0_5's solver shares. Nothing in it is KrF's.

### Verified
- Before the move, each of the 374 removed files (91.5 MB) was
  byte-identical to `LIGHT/` at the same path, and the script would not
  move anything otherwise. Besides them went 27 .DS_Store files and
  bytecode caches. Nothing existed only here.
- After the removal, from Kinetic_Modeling: `v0_5/tests` 342 passed and 2
  skipped in 62 min (load average near 35). One skip is the removal's:
  `test_ib.py` compares against LIGHT/v0's kernel, which is gone here. The
  other is the slow reference-point test (LIGHT_KINETIC_SLOW=1).
- `python3 LIGHT_Kr/run_v0_5.py --f2-pct 0 --eedf two_t` runs: exit 0,
  ledger closure -1.4e-16, every fluorine and KrF line 0, and one
  RuntimeWarning from the report's "F2 left" line (0/0).
