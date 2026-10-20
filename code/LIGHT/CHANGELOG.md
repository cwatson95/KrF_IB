# Changelog: LIGHT (the lightweight 0-D KrF model)

Dated entries, newest first. Read the Known-bad list before editing.

## Known-bad: do not reintroduce

- **With the wave alone, solve the kinetic EEDF steady; time-step it
  under the given power only (v0_6, 2026-10-19).** Time-stepping it
  whenever the wave was on failed in the attaching afterglow of 2 % F2:
  at the envelope point two of six runs found no time step (the 1 V/cm
  control one of them), and a third found no steady solution as the
  window closed. A steady solution that cannot be found now falls back to
  a time step from the last one, logged `stepped`.
- **Call `rhs` and `_kinetic_update` with v0_5's three arguments wherever
  no window parts the wave from the given power (v0_6).** The beam drive
  in `docs/survey_2026-09-27/scripts/ebeam_screen.py` overrides
  `rhs(t, y, on)` and `test_kinetic.py`'s counters wrap
  `_kinetic_update(t, y, on)`; a fourth argument everywhere broke four
  tests (2026-10-19). `Model._legacy` keeps the old call;
  `test_a_subclass_written_for_v0_5_still_runs` holds it.
- **A twin pair's marginal means nothing where the wave absorbs nothing
  (v0_6).** The 1 V/cm control's dE_248 / E_wave came out at 108 % and
  -5600 %, ratios of 1e-12 J/cm3 of noise; `reheat.py`'s tables print a
  dash there.
- **A row's factor must reach the kinetic EEDF's rates (v0_4, v0_5;
  fixed 2026-10-14).** Under the kinetic EEDF the rows the kinetic table
  carries (G1, C1a-c, C2, the F2 electron rows ATT, DIS and VIB, R11p,
  and SUPER when added) read the table's rate, which was made from their
  unscaled cross sections and recombination template, so `row_scale`
  moved none of them; in v0_4 K76's pseudo cross section took the band's
  value without the row's factor too. Each row's factor now goes into its
  process in the Boltzmann equation. K76 and K76r take one factor (the
  reverse is built from the forward's cross section), and G6, the
  elastic exchange, cannot be scaled.
  `test_row_scale_reaches_the_kinetic_table` (in each version's
  `test_v0_4.py`) holds it.
- **A generated-numbers check must not compare wall-clock keys (v0_5
  docs, 2026-10-11).** The overview generators print how long their
  runs took (*.min); a `--check` that compares the whole file against a
  fresh run fails on them every time. Compare without them.
- **Never compare the solver's electron-ion wave power with the kernel's
  printed "ei" share, and never let the ib_scale band mean two things
  (v0_4, 2026-10-09).** The kernel's "ei" is nu_ei's share of the
  numerator of nu/(nu^2 + omega^2); the solver's, with the ions taken out
  of its momentum transfer, is the marginal part. They agree only for
  nu_en << omega; in 1 atm Kr at 95 GHz nu_en = 4.4 omega and the share
  is 1.5x the marginal part. A 50 %-tolerance test hid the 33 % for a
  day. Compare with the marginal part; the band acts on it under both
  closures.
- **Never apply the Langdon factor under the kinetic closure (v0_4).** It
  corrects a Maxwellian-averaged electron-ion absorption for the
  distortion a solved EEDF already holds; both would double count. ib.py
  applies it only where model.py's closure is two_t or maxwellian.
- **Every electron row of the ladder that moves energy must be in the
  Boltzmann equation (v0_4).** With K76/K76r as LIGHT-side rates alone the
  kinetic energy balance missed LIGHT's by 0.7 of the input (2026-10-08).
  Kannari's A Te^0.5 exp(-C/Te) is the Maxwellian rate of sigma0 (1 -
  C/eps): use that pseudo cross section (boltzmann.PopProcess), never a
  rate the solver cannot see. The stepwise rows are LIGHT-side still and
  are bands at 0 for that reason too.
- **An error band on an absorption acts on the booked power, never inside
  the collision frequency (v0_4).** nu_ei scaled inside nu/(nu^2 +
  omega^2) is nonlinear in the electron-ion share below omega ~ nu
  (10 % of itself at 1 THz, 1 atm Kr).
- **Kr2F*'s 248-nm cross section is Hakuta et al. 1987's measured 1.0e-18
  cm2, never the compilation's 6.0e-18** (Kannari Table I [16], Johnson &
  Hunter), which gives Tahoma's alpha0 3.4 times over (2026-10-08).
- **A 5p rung needs its collisional relaxation.** With Kr(5p) + Kr at 0
  the rung could only radiate and 17 % of the input left as near-IR
  light; the measured 1.7e-10 (Gazeli et al. 2021) is the default
  (2026-10-08).
- **Keep the event-bracket guard in `Model._run_kinetic` (v0_4).** scipy's
  solve_ivp raises "f(a) and f(b) must have different signs" where a
  trigger's sign change between a step's ends is not reproduced by its
  dense output (an input jumping, n_e at its floor); the guard takes a
  short step without the triggers and renews the table (2026-10-09).
- **No Te = U_e/(1.5 n_e) in a transport or loss term (v0_4).** Where n_e
  is a small difference of two large ion densities (the attaching
  afterglow) such a term's Jacobian explodes while its value stays small,
  and the stiff solver stalls (a 0.3-us run took over 10 minutes). Drive
  it by U_e over the cations, smooth and capped (`Model.wall_te`).
- **A population process scales each operator by its own population
  (v0_4).** The 5s -> 5p forward acts on Kr*, the reverse on Kr**; one
  fraction for both returned 7 times too much energy to the electrons.
- **Particles that leave the box take their thermal energy with them, once
  (v0_4).** Not through the reaction rule in the Tg equation (it raised
  Tg), and not as a second fate line (E_wall_th is inside the heat booked).
- **No tuning.** Never move a rate to bring an output nearer a measured
  efficiency (the repo's predictive-only policy: anchors score, never
  steer). LIGHT copies SPECK/REBUILD's sourced values, and
  `test_values_copied_from_rebuild_have_not_drifted` fails if they diverge.
  Change a value in REBUILD first, with its source, then here. Since
  2026-09-28 the basic network also carries its own measured values, in
  rows without a `rebuild=` tag. Each such row cites a measurement of that
  process, never a fit to device data, and v0 keeps REBUILD's copies.
- **The vendored cross sections are pinned.** `v0/data/*.txt` are md5-pinned
  copies of the LoKI-MC inputs. Never edit them in place. Re-vendor and
  re-pin deliberately (`xsec.FILES`).
- **dTg/dt = heat / C_v alone is wrong** when reactions change the particle
  count (F2 -> 2F): the new particles take cv k Tg from the heat. The first
  draft had it; the fix subtracts kTg sum(cv_s dn_s/dt), and
  `u_thermal` vs booked heat (to 1e-9 of E_dep over the whole trajectory)
  is the gate.
- **No BLAS matmul on the Maxwellian tables.** numpy 2 on macOS Accelerate
  raises spurious divide/overflow RuntimeWarnings in a matmul whose
  operands hold underflowed zeros (the numbers were right). `xsec._mean`
  multiplies and sums instead.
- **Superelastic rates come from the Klein-Rosseland cross section**, never
  from k_exc (g ratio) exp(E/Te): at room-temperature Te that multiplies an
  underflowed table by exp(+381).
- **A zero-length integration segment crashes solve_ivp** ("first_step
  exceeds bounds"). `Model.run` drops segments with no length (pulse_s = 0
  runs are the analytic-limit tests).
- **Do not read Te after n_e hits the floor.** With the `f2e` group and
  `ne_mode="free"`, attachment empties the electrons in nanoseconds
  (k_att = 6.6e-9 cm^3/s at 0.5 eV on 1e17 F2). Below about 1e3 cm^-3,
  Te = U_e/(1.5 n_e) is solver noise.
- **The drift gate reads the SPECK beside LIGHT.** In a checkout of commit
  4a70eb7 (REBUILD before 0.1.40, which is how a fresh worktree of this
  submodule starts on 2026-09-25) it fails on G5's Tg exponent, and rightly:
  LIGHT follows the live 0.1.40 (Shui's Tg form restored), which is
  uncommitted in the main checkout. Test from such a tree with
  `LIGHT_SPECK_DIR=<main>/SPECK`; never "fix" the gate by editing LIGHT's
  values back.
- **Two TeX traps in the docs .tex files.** The T1 typewriter font turns
  `--` into an en-dash, so `--omissions` must be written `-{}-omissions`.
  sansmath has no upright-sans uppercase Greek, so an uppercase Greek
  letter renders as a black box unless it comes from the normal math font
  via `\unsansmath`. The maps' 4^2 Gamma label and the overview's F2 state
  labels (a 3Pi_u, A 1Pi_u) do this.
- **The measured KrF constants belong to the two-level B/C model.** Xu and
  Gadomski fitted k_F2, k_Kr and k'_Kr jointly with K_eq = 5.0 in a
  two-level KrF(B,C) scheme (Gadomski 1992 eq. 5). Running them in a B-only
  model (v0) or with Kannari's K_eq = 1.47 misses their decays by x1.5 to
  x2.9. Keep them in the basic network's structure. Change them only with
  a new measurement, never to match another dataset.
- **Never run the `f2e` group with `ne_mode="clamped"` for results.**
  Attachment removes electrons that the clamp keeps replacing, so F⁻ is made
  from nothing and charge is not conserved. The ion store read −25 % of the
  deposited energy at 2 atm (survey 2026-09-27). Use free n_e with `f2e`.
  In v0_1 the basic network itself carries attachment and recombination,
  so its n_e is free by default and the report flags a held one.
- **Pair Tahoma's gain with its own deposition.** The 330 kW/cm³ at the
  probed centre belongs to the beam-stop-in runs, whose g0 is 6.39 %/cm
  (report Table 3-3, p. 3-16). The beam-stop-out 6.30 %/cm was taken at
  about 20 % less pumping. The legacy benchmarks pair 5.45 %/cm, which is
  g0 − α0 of the 14.5 % Kr Ar-diluent mix, with 285 kW/cm³. Never mix them.
- **The basic network's B/C constants are a joint fit, not single
  measurements.** Since 2026-09-28 tau_C, k_BC, K_eq and the B and C quench
  constants come from ML's fit to Eden 1978, Gadomski 1992 and Yu & Setser
  1990 (ML/README.md). tau_B is Kannari's 6.99 ns. Hakuta 1990's measured
  tau_C and k_BC disagree with those three and are held out
  (test_basic_network_against_hakuta_c_decays, strict xfail). Do not
  describe the constants as measured one by one, and do not refit them to
  device data.
- **K_eq is not a constant.** The C -> B row is a "balance" row: k_BC /
  K_eq(Tg), with the gap taken from the species energies. A fixed K_eq with
  a nonzero ledger gap breaks detailed balance, and the first basic network
  had exactly that (fixed 5.0 against a booked 80 cm^-1).
- **`sources/jenkin1990.pdf` is not Yu & Setser.** It is Jenkin et al.,
  J. Phys. Chem. 94, 2927 (1990). The real article is
  `sources/Yu_Setser_JPC94_2934_1990.pdf`.
- **Yu & Setser print K_eq 1.7 +- 0.2, not "2-3".** Xu 1993 p.2598
  paraphrased it as 2-3. Their value rests on tau_C 75.2 ns, so with a
  measured tau_C near 50-60 ns the same ratio gives about 1.1-1.4.
- **A fitter cannot settle conflicting experiments.** A neural network fed
  all four B/C datasets returned a confident compromise that misfits Kansas
  State at chi2 8.6 per point, with no warning. Check dataset consistency
  first (leave-one-out chi2 in `bc_infer.py`), then fit.
- **Never read a KrF Stern-Volmer product with a one-state lifetime here.**
  In the two-level network the lifetime in k_e tau / k_att is the
  Stern-Volmer lifetime of B-X emission (`v0/two_level.py` `sv_tau`). It depends on k_BC, K_eq and where KrF is formed. KHS's 3-6e-8
  used 5.5 ns for KrF(B) alone. Change k_BC or K_eq and k_e must be
  re-derived (`test_electron_quench_is_trainor_jacob_under_this_b_c_model`
  enforces it).
- **The Kansas State decays cannot fix k_BC or K_eq.** With the other
  two-level constants free they trade along a valley (K_eq 0.4 to 4.6 for
  k_BC 3e-12 to 1e-10). The x0.89 / x0.95 match is not support for K_eq = 5.
- **"Flows" are not a partition.** They overlap and can exceed 100 %
  (Penning returns energy to the electrons, which excite Kr again). Only
  the "fate at t_end" block sums to 100 %.
- **The tail table is generated, never edited.** `v0/data/kr_tail_ratio.txt`
  comes from `python3 ML/run.py tail` and is md5-pinned (`xsec.TAIL_FILE`).
  Never set the tail ratio by hand or fit it to device data. A chosen tail
  temperature is a knob (REBUILD 0.1.37's two-temperature tier is on
  REBUILD's Known-bad list for the same reason).
- **(3/2)Te per captured electron is v0's closure only.** F2 attaches slow
  electrons: an attached electron takes 0.54 eV at Te = 1.13 eV, not 1.69.
  Booking (3/2)Te put LIGHT's Tahoma Te at 1.13 eV where the kinetic EEDF
  gives 1.60 (2026-09-28).

- **`pytest.approx` passes any two numbers below 1e-12.** Its default
  absolute tolerance, 1e-12, applies even when a relative one is given, so
  every comparison of rate constants is vacuous unless the test passes
  `abs=0`. v0's tests use it as it comes. There the check that G5 is
  6.0e-31 passes against the network's 1.9e-30, and the REBUILD drift gate
  cannot
  see a change in a three-body rate (a strict re-run on 2026-09-28 found
  no drift). v0_1's tests compare through an `approx` helper with `abs=0`.
  Use it for any new test.
- **A bare % in generated LaTeX comments out the rest of its line.** In a
  table row that includes the row's \\, so the next row merges into it and
  LaTeX reports extra alignment tabs one table later. Two short sources
  ("96.6 %", "+- 10 %") did this in the v0_1 overview (2026-09-29). Escape
  with one substitution, `re.sub(r" ?%", r"\\,\\%", s)`. A second replace
  over its output escapes the % twice and makes `\\%`, a line break and a
  comment again.
- **ML binds its LIGHT version when it is first imported.** Never set
  `LIGHT_VERSION` at a module's import time, since a test that imports that
  module would rebind ML for every test after it. Generators set it in
  `main()`, ML-backed functions refuse another version, and tests run them
  in a process of their own (`--check`).
- **Say which configuration a device score used.** v0's overview scored
  Tahoma at x0.96 with v0's optional ion, F2-electron and Kr2* groups on and
  n_e free (ML's `forward.tahoma_params`) and did not say so. v0's default
  network, the one its map and reference point describe, gives x0.32 there.
  Every device score names its groups and n_e mode.
- **The Kr2* group's gain is one estimated rate.** Kr2* + F -> KrF(B) + Kr
  (K2F_F) is Kannari's estimate, 3.0e-10 cm3/s, and at the reference point
  it carries the group's whole gain, 4.73 -> 9.17 % at 248 nm and 4.65 %
  without it. Most of it comes after F2 burns out, where the F atoms left
  behind turn Kr2* into KrF and get their F back with each photon. Never
  call the Kr2* route measured, never make the group a default on the
  strength of its gain, and report it as a band until it is measured
  (2026-09-29). 3.0e-10 is no upper bound. It is the low end of a
  conditional range SRI gave for Ar2* + F. The share saturates, 4.65 % at
  zero, 9.49 % at 1e-9.
- **A band's default is a measurement or 0.** In v0_2's `network.BANDS` a
  rate with no measurement is 0 by default, and its band runs from 0 to
  the largest published value, which is no bound. Never set a default
  from an estimate, an analog or a device score, and never narrow a band
  to one. Jancaitis's F-atom quench of Ar2F* only sets the top of Kr2F* +
  F's band (2026-09-29).
- **The 1-D models' tests need their `tests/conftest.py` (1Dv0,
  1Dv0_2, 1Dv0_3).** pytest stops its package search at a directory
  name that is not an identifier, so for '1Dv0' it never puts LIGHT/ on sys.path, and `import_module("1Dv0.x")`
  finds no package, or the copy in the current directory. A mutation
  check ran 16 mutants against the unmutated package that way and passed
  them all (2026-09-29).
- **Import Fluid1D only through `fluid1d.load()` (1Dv0, 1Dv0_2,
  1Dv0_3).** The f2sweep venv
  holds an editable install of another checkout
  (`~/coding/nathaniel/1DModelKrF`), which a bare `import Fluid1D` finds
  first. `load()` puts `nathaniel/1DModelKrF/src` beside LIGHT (or
  `FLUID1D_SRC`) first and refuses a Fluid1D from anywhere else.
- **The 1-D stage (1Dv0, 1Dv0_2, 1Dv0_3) reads the rates at each substep's
  starting energy.**
  Fluid1D's own loop adds the source first, which with v0_1's power set a
  steady Te 3 to 4 % low at dt = 50 ps. And the given power reaches each
  half-step over its own time window, not through Fluid1D's lead/trail
  predictor, which after the pulse's end took the energy back in the next
  half-step (stages capped at 8000 substeps).
- **Never quote a 1-D gap average (1Dv0, 1Dv0_2, 1Dv0_3) without its
  mesh.** On a uniform mesh the cells next to the walls are unresolved
  sheaths: the mean-energy ceiling cuts their input and Fluid1D's split
  keeps out a phantom negative field work, both shrinking with dx. The
  interior is v0_1's 0-D run. On 1Dv0_3's graded mesh the sheath is
  resolved and the 200 um beside it are not converged (below).
- **Give the 1-D models (1Dv0, 1Dv0_2, 1Dv0_3) a seed of 1e13 cm^-3 or more at
  hundreds of kW/cm^3.**
  From ML's Tahoma seed, 1e10, it never breaks down: v0_1 does only because
  its 0-D Te runs to 518 eV, and 1Dv0's tables stop Te at 40 eV, below
  which v0_1's tail holds no ionizing electrons at that ionization degree.
  The fate then shows the whole input "cut by the mean-energy bounds".
- **An & in generated LaTeX is a column break.** v0_2's band citations
  name authors with "&" (Salesky & Kimura, Quigley & Hughes), and printed
  as they are they split the rows of the overview's bands table into extra
  columns (2026-09-30). `make_overview_v0_2._tex_text` writes "and", as
  the hand-written short sources do.
- **pdfpages does not flush floats.** A float still pending at an
  `\includepdf` lands after the included page, and with `fitpaper` on a
  page of the included PDF's size (2026-09-30: the energy ladder came out
  after the system map, on a landscape page that cut off its caption). A
  longtable holds pending floats until it ends. The overviews set
  `\clearpage` before the map, and 1Dv0_2's before its appendix.
- **No kinetic rates carried as factors from one mean energy to another.**
  v0_3's first coupling took ML's method, factors over the Maxwellian at
  the kinetic EEDF's own mean energy applied at LIGHT's Te. After the pulse
  the two means parted (0.29 against 1.0 eV), and the factors (22 on the
  elastic loss, 6e4 on F2's dissociation) drove LIGHT's Te from 1.1 to 0.41
  eV in 0.12 ns (2026-09-30). v0_3's table is at LIGHT's own mean energy.
  Past the table's ends LIGHT carries the end's values with its own
  Maxwellian's change between the end's Te and its own (2026-09-30 (e)).
  That moves no kinetic value by a kinetic-over-Maxwellian factor, meets
  the table at its end, and is needed, since held flat there the table
  stalled LIGHT's Te at Tahoma (1.3 against 1.6 eV).
- **No search on E/N for a mean energy.** At a fixed field the two-term
  problem with strong e-e collisions can have two solutions (Kr's Ramsauer
  minimum), and Brent's method on E/N finds no sign change. v0_3 matches
  the mean inside every step of the iteration (`boltzmann._mean_match`).
- **No e-e integrals with F constant in each cell where kT is near the
  cell width.** It condenses an undriven EEDF into the first cell, its mean
  falling with the cell width. `Grid.ee_integrals` takes ln F linear
  between cells, which leaves a Maxwellian alone.
- **No grid top scaled by the decades F falls across it.** It oscillates
  between grids too long and too short, and one returned with Kr's
  ionization threshold off its end. `Solver._emax_for` reads the top off F.
- **No elastic or field power by central differences on the kinetic grid
  (v0_3).** The solver's fluxes are Scharfetter-Gummel. Central differences
  gave its own undriven equilibrium a loss to the gas of 4.7 % of its
  cooling term on 600 cells to 12 eV at 300 K, and LIGHT's late-afterglow
  electrons fell to 0.93 Tg (2026-09-30 (e)). The powers are
  `Solver.elastic_power` and `field_power`, the solver's own flux split.
- **Keep the stall guard in `Model._run_kinetic` (v0_3).** solve_ivp
  locates an event only to 8.9e-16 s (scipy 1.13), so a trigger crossed in
  the first femtosecond can come back at t = 0, and the run then renews its
  table there without end (12,405 times in a pure-Kr beam run, 2026-09-30
  (e)).
- **No growth term dropped where the gas relaxes the shape faster than the
  population turns over (v0_3).** Without it the recombined electrons come
  back cold in the first cell, and a pure-Kr afterglow went below the gas
  temperature (2026-09-30 (e)). `growth_valid`'s relaxation clause keeps it.
- **Rank every score before naming the largest misses.** The bottom lines
  named Hakuta's decays and the pure-Kr yield, a sentence copied from v0's
  box, and left out Burnham and Searles' tau_B (-22.3 %) in all five
  overviews and v0's default network at Tahoma (32 % of the gain) in v0's
  (2026-09-30 (f)). v0_3's `test_docs.py` checks its box. The others were
  checked by hand.
- **Never move a time step's history onto another grid without keeping its
  mean (v0_3, `Grid.remap`).** Where a cell is as wide as kT, a plain
  interpolation of the gas's Maxwellian between the solver's grids moves
  its mean by 0.5 to 27 %, which over a picosecond's step reads as 2e8 to
  1e10 eV/s of heating per electron, against the test's Joule power of 1e9
  (2026-10-01).
- **No rescale of the e-e drift to force energy conservation (v0_3).**
  Where a cell is as wide as kT the energy error lies in splitting the
  Scharfetter-Gummel flux, whose c(z) couples its parts. A rescale of the
  e-e drift that zeroed its own exchange left the step's error at 10 %
  (lambda 0.999994). Resolution cures it (2026-10-01).
- **No time step after the pulse yet (v0_3).** Through the afterglow of
  Morozov's low-pump pure-Kr run it left 87 of the 95 afterglow updates
  short of LIGHT's mean energy (2026-10-01). `kin_memory_s` acts while
  the drive is on.
- **A step's field over the input power is not its energy check (v0_3).**
  It carries the step's average d<eps>/dt against LIGHT's instantaneous
  one, and the table's renewal. The energy-balance diagnostic is the
  check.
- **Never explain an energy ratio by a term without checking the term's
  size.** The field's 0.6 to 3.1 % excess over the input power was read
  for a day as the electrons cooling through the pulse. Their cooling is 6
  to 44 times smaller, and the excess is the table's renewal (2026-10-01).
- **On a graded mesh nothing may read Fluid1D's one dx or cell volume
  (1Dv0_3).** Fluid1D's transport, field solve and books assume one cell
  size. 1Dv0_3 unsets both on a graded mesh, so a read fails instead of
  quietly using a wrong size, and its transport is `transport.py`. A
  Fluid1D path new to the run is redone per cell, and its hash pinned in
  `test_1d.py` (2026-10-01 (b)).
- **Read nothing within 200 um of a wall on a graded mesh as converged
  (1Dv0_3).** On cells of a few microns the F2 burns out there a band at a
  time (the halogen-depletion instability). With the gas still the bands
  reach 2000 K beside cells at 600 K. With it moving (the default since
  2026-10-01 (c)) they are cooler, four at each wall near 1000 to 1200 K
  on 256 cells where the still gas put them and two of them on the default
  128. With Fluid1D's field work the step's error drove one layer on 64
  cells to 2000 to 2800 K, and with the work kept to the field's energy
  (the default since 2026-10-06) 64 cells put the two nearest bands where
  256 do. The field's work at the bands' edges moves with them (0.49 to
  0.29 % of the given power from 128 to 256 cells). The interior is the
  result. The gap's totals as shares of the given power move 0.3 % from
  128 to 256 cells and 0.02 % in dt with the work kept to its energy, 0.5
  and 0.2 % with Fluid1D's (2026-10-01 (b), corrected in (c), 2026-10-06).
- **No given power shared per electron, and no field work from the faces,
  in the 1-D models while the neutrals stand still.** Both were tried on
  1Dv0_3's graded mesh against the near-wall cells. Shared per electron,
  the gap contracts at once (its mean n_e 3.7e16 cm^-3 by 200 ns), and the
  faces' work made the near-wall cells run away sooner (2.8e16 by 1.2 us)
  (2026-10-01 (b)). 1Dv0_3's gas moves since 2026-10-01 (c), and neither
  has been tried with it. The faces' products do not conserve the field's
  work either (+0.69 % of the given power to 600 ns on 64 cells, (c)).
  The conserving coupling (2026-10-06) also takes the work at the faces,
  but with its sign, paired with the time-centred exact fields, the ions'
  share to the gas and the field's energy in the books. This item rules
  out the faces' work given as input, not that.
- **The field's work in the grounded gap is never an input (1Dv0_3).**
  From a start with no field its net work on all charges is minus its
  energy (9e-7 of the input at 600 ns). Fluid1D's step
  (`field_work="fluid1d"`) gives the electrons its positive part as input,
  and to 600 ns on 64 cells that is +0.30 % of the given power with the
  gas still and +1.9 % with it moving, all of it the coupled step's error.
  With `fluid1d` read the shares of the given power (`share_248()`,
  `--compare`, the summary's "of the given power alone"), never `fate()`'s
  shares of everything put in where the field work is percents
  (2026-10-01 (c)). The default since 2026-10-06 keeps the work to the
  field's energy and none of it is input. On the uniform mesh its wall
  cells hold the whole sheath, and the floor puts back what the sheath's
  work takes beyond the electrons' energy (2.3 % of the given power on 64
  cells), so read the field's work there as the unresolved sheath's.
- **No density floor in the conserving transport (1Dv0_3).** Fluid1D's
  species update sets a density that goes negative to zero, which makes
  particles and charge, and the field's work then misses the field's
  energy (by 44 % of the work on 48 cells from 10 um). With
  `field_work="conserve"` the outflow a cell cannot supply is scaled to its
  content instead, each face by its upwind cell's factor (2026-10-06).
- **The field's work is each face's current times the time-centred exact
  field (1Dv0_3).** The theorem needs the exact fields of the charge
  before and after each substep and the currents that moved it. One end's
  field alone, or a cell's mean flux and mean field (Fluid1D's), misses
  the field's energy. `test_conserving_work_keeps_the_field_s_energy` and
  the twin advances hold it (2026-10-06).
- **The gas carries Fluid1D's electron energy density, not its mean
  energy (1Dv0_3).** Fluid1D blends the mean energy toward the local
  field's where the electrons are few, so there it differs from the energy
  the books hold. Carried with the gas it left 2e-7 % of the input
  unbooked (2026-10-01 (c)).
- **Fluid1D's Simulation keeps its gas mixture in `self.gas` (1Dv0_3).** A
  LightSimulation attribute of that name hides it, so the gas solver is
  `self.flow` (2026-10-01 (c)).
- **1Dv0_3 on its uniform mesh stays 1Dv0_2 to the bit.**
  `test_uniform_mesh_is_1dv0_2_to_the_bit` holds it, and `transport.py`
  keeps Fluid1D's arithmetic exactly there (face weights 0.5 and 0.5,
  dx/dx = 1, the cells' sums in Fluid1D's order). A change that breaks it
  changes 1Dv0_2's model, not only the mesh (2026-10-01 (b)).

## 2026-10-19: v0_6, the wave's window, and the afterglow reheat tested

The owner asked to verify with FLASH 4.8 that "we can reheat the afterglow
of the KrF gas to improve efficiency" (with `IB/envelope_kr/envelope_kr_v6.*`,
this directory and `docs/accuracy_overview_v0_5.*`). FLASH cannot: its
laser heats an LTE plasma by electron-ion inverse bremsstrahlung, with no
electron-atom absorption, no Kr excitation, no attachment and no KrF
kinetics. LIGHT carries all of them, so the test was proposed here, and
the owner asked: "make a v0_6 that with the three things you would test
for".

### Added
- `v0_6/` from v0_5. `Params.em_t_on_s` and `em_t_off_s` (`--em-on`,
  `--em-off`) give the EM wave's window; the default is v0_5's, on with
  the pulse. `Model.segments` cuts the run at pulse_s and at both window
  edges. `rhs`, `_kinetic_state` and `_kinetic_update` take the given
  power's flag and the wave's apart; the wave's defaults to the given
  power's.
- `v0_6/reheat.py` (`python3 LIGHT/run_v0_6.py reheat`), the twin pairs:
  - at the reference point and at the envelope notes' worked example
    (0.5 atm, 2 % F2), with the notes' wave (E0 = 47.6 kV/cm, 33.7 kV/cm
    rms, at 6.2 and 6.6 THz);
  - three fields, two windows, both closures, and a 1 V/cm control;
  - for each pair: the marginal (dE_248 / E_wave), the efficiency with
    the wave counted as input, the break-even ratio of source
    efficiencies, where the wave's energy went, and the window's optics;
  - every band and UNC_ROWS constant at both ends: under two_t at both
    points, under the kinetic EEDF at the reference point
    (`--kinetic-bands`), and FURTHER's pairs (`--further`).

  The runs are resumable (`results/reheat_runs.json`), and the tables go
  to `results/reheat.md`.
- `run.py --twin`: one pair from the command line.
- `report.py`: the window in the summary, the wave-alone kinetic line,
  and `row_dict`'s given_J_cm3, E248_J_cm3, em_on_us and em_off_us.
- `run_v0_6.py`; `v0_6/tests/test_v0_6.py` (15 tests); `v0_6/SOURCES.md`
  section 9. `test_docs.py` was left out, since v0_6 has no documents of
  its own.

### Results (`v0_6/results/reheat.md`)
- **At LIGHT's defaults the reheat lowers the 248-nm efficiency** with
  the wave counted as input, at both points, under both closures, at
  all three fields and in both windows:
  - reference point, kinetic: -0.037 pp (marginal 4.5 % against the
    pump's 5.73 %);
  - reference point, two_t: -2.25 pp. The wave sustains a discharge,
    absorbs 75 % of the pump's energy, and adds 0.44 % of it as 248 nm;
  - envelope point, two_t: -0.033 pp (9.1 % against 14.75 %);
  - envelope point, kinetic: -0.0014 pp (13.9 % against 15.83 %).

  A gentle wave (1/100 of the power) has a negative marginal everywhere
  (-0.3 to -3.1 %).
- **One rate with no measurement decides the reference point under the
  kinetic EEDF.** Kr2* + F -> KrF(B) (kr2_f, 0 by default) matters
  there because the 0.5 % F2 burns out: the wave's Kr* becomes Kr2*,
  which at the default leaves as VUV (30 % of the wave). Scanned across
  its band, the change in efficiency:

  | kr2_f, cm3/s | change |
  |---|---|
  | 0 | -0.037 pp |
  | 1e-12 | -0.018 pp |
  | 3e-12 | +0.014 pp |
  | 1e-11 | +0.077 pp |
  | 3e-11 | +0.127 pp |
  | 1e-10 | +0.138 pp |
  | 3e-10 (Kannari) | +0.131 pp |
  | 1e-9 | +0.123 pp |

  The sign turns at about 2e-12 cm3/s. Above it the reheat pays, by at
  most 0.14 pp of about 7 %, with the notes' 3 MW/cm2 wave (the gentle
  one does not pay) and a source at least half as efficient as the
  pump.
- **The kinetic bands at the reference point** (all 20 inputs): -0.037
  +- 0.84 pp in quadrature, envelope -2.57 to +0.35 pp. Four ends raise
  it: kr2_f 1e-9 (+0.12), kstep_s 1.5e-7 (+0.071), kstep_p 2e-7 (+0.015)
  and k_5s5p 1.4 (+0.002).
- **Where it never pays:**
  - under two_t, where every band end lowers it, kr2_f's top included
    (-1.37 pp at the reference point);
  - at the envelope point under either closure. Its 2 % F2 attaches the
    electrons within about 10 ns and does not burn out, so kr2_f's top
    gives -0.0012 pp. Under the kinetic EEDF every end run there lowers
    it: -0.0014 +- 0.0004 pp over kr2_f, krf_f and the trimer law, and
    -0.0014 at k_5s5p's upper end.
- **The closures disagree 24-fold on the absorbed energy at the
  reference point.** two_t's Maxwellian with the e-beam-calibrated tail
  absorbs 5.2e9 eV/s per electron at 2.2 eV and holds the plasma; the
  solved EEDF absorbs 2.0e9 eV/s at 1.06 eV.

### How it ran
- The memory reaper stopped one batch: 28 workers ran beside another
  session's decompilers. A reboot ended the session that followed. Each
  finished run had been saved, and `python3 LIGHT/run_v0_6.py reheat
  --kinetic-bands --further --workers 8` finished the rest from the
  cache.
- One run was stopped unfinished: the envelope point's kinetic pair at
  k_5s5p = 0.5, after 3.4 hours of compute, twice its siblings'
  longest. Whether it stalls is not known. The same command runs it
  again.
- Fixed on the way: a failed run's error record reached the tables
  (`KeyError: 'given'`); `reheat._rec` now treats it as missing, and the
  table prints it as failed.

### Found on the way
- **A time-stepped EEDF under the wave alone failed** in the attaching
  afterglow of 2 % F2. At the envelope point two of six runs found no
  time step in the window (the 1 V/cm control one of them), and one
  found no steady solution as the window closed.
  - v0_6 time-steps the EEDF under the given power alone, v0_5's rule,
    and solves it steady with the wave alone, as after the pulse.
  - A steady solution that cannot be found is replaced by a time step
    from the last one, logged `stepped`: 2 of the study's 37,097
    kinetic solutions.
- **Passing the wave's flag everywhere broke four tests** in
  `test_kinetic.py` (the beam drive and the counters). The beam drive in
  `docs/survey_2026-09-27/scripts/ebeam_screen.py` overrides
  `rhs(t, y, on)`, and the tests wrap `_kinetic_update(t, y, on)`. The
  flag is now passed only where a window parts it from the given power
  (`Model._legacy`).
- **Shell lesson:** a `kill` of the study's workers from zsh left 14
  of them orphaned under launchd. zsh does not split an unquoted list of
  PIDs, so they were stopped one by one, and the next stop was checked
  with `ps`.

### Verified
- Default runs are v0_5's to the bit: densities, times and ledger
  against v0_5 itself under two_t, in three cases, with and without the
  wave. The kinetic reference point gives v0_5's 5.73 % at 248 nm,
  4.82 % VUV and 738 K.
- Every saved run of the study (250) closes its ledger to 1e-12 of its
  input.
  The 1 V/cm control moves 248 nm by 1e-12 J/cm3 under either closure.
- Tests on 2026-10-19, final code, `/usr/bin/python3`:
  - `test_v0_6.py`: 15 passed;
  - the other fast files: 327 passed and 1 xfailed (Hakuta's decays,
    strict, as in v0_5);
  - `test_kinetic.py`: 24 passed and 1 skipped.

## 2026-10-14: row_scale reaches the kinetic EEDF (v0_4, v0_5)

The owner asked: "Fix the bug under all LIGHT directories LIGHT/
LIGHT_ArF/". The bug was found in LIGHT_Kr v0_1 on 2026-10-13
(`LIGHT_Kr/CHANGELOG.md`, entry (c)), whose code is v0_5's.

### Fixed
- `v0_4/model.py`, `v0_5/model.py`: `Params.row_scale` reached no row
  the kinetic table carries. Symptom (LIGHT_Kr v0_1, the same code): R11p
  at +-10 % moved the kinetic run's R11p events by 1e-10 and its 147-nm
  share not at all, where the two-temperature closure moved the share by
  1 %. Root cause: `Model._rate` reads the kinetic table's rate without
  the row's `c.factor`, and past the table's ends the Maxwellian ratio
  cancels it; the processes the Boltzmann equation solves were built from
  the unscaled cross sections (`xsec` rows), g ratio (`superelastic`) and
  recombination constant (`te_pow`). K76 and K76r carried their factor
  in v0_5 (the band k_5s5p scales Jung et al.'s cross section) but not in
  v0_4, whose Kannari pseudo cross section took the band's pre-exponential
  alone. Fix: `_kinetic_setup` builds a cross-section row's process on
  its cross section times the factor, multiplies the superelastic g ratio
  and the recombination constant by it, and (v0_4) multiplies K76's
  pseudo cross section by K76's factor. `Model.__init__` refuses a
  row_scale on G6 (neither closure scales the elastic exchange) and
  different factors on K76 and K76r (under the kinetic EEDF the reverse
  comes from the forward's cross section, so a factor on one alone meant
  one thing on each closure).
- `v0_4/tests/test_v0_4.py`, `v0_5/tests/test_v0_4.py`:
  `test_row_scale_reaches_the_kinetic_table`, written first and failing on
  the old code in both: the scaled cross section and recombination
  constant in the kinetic processes; at one state R11p's rate 1.1 times
  the default's (1e-9), K76's and K76r's 1.1 times with both scaled, G1's
  0.8 of its cross section's rate over the run's own EEDF (1e-6); the
  other rows untouched; G6, K76 alone, K76r alone and unequal factors
  refused.

### What the defect touched
- No document number. Every band leaves the table's rows at a factor of
  1 at both its ends (checked band by band), so only a row_scale on G1,
  C1, C2, the F2 electron rows, R11p, SUPER or K76r was silenced, and
  only under the kinetic EEDF. The overviews' propagated uncertainty
  scaled R11p, K15 and C3 on v0_2's closure (`make_overview_v0_4.py`,
  `make_overview_v0_5.py`, `--unc`), and the Kleikamp sensitivity scales
  G5 and G5C, heavy rows the table does not carry.

### Verified
- Default runs are unchanged to the bit: in v0_4 and v0_5 the
  two-temperature reference and a 0.1-us kinetic pulse, run before and
  after the fix, give the same row_dict, final densities, U_e, Tg and RHS
  count (44 numbers each, compared as reprs). For every default run the
  new code takes the old path (no table row has a factor other than 1, and
  x 1.0 is exact).
- From Kinetic_Modeling (load near 25): `v0_4/tests` 348 passed, 5
  skipped, 1 xfailed (44 min); `v0_5/tests` 358 passed, 5 skipped, 1
  xfailed (51 min), both including `test_docs.py`, whose generators hash
  no code.

## 2026-10-11: v0_5's accuracy overview and system map

The owner (2026-10-11): "Make the accuracy overview and system map for
v0_5". Both on v0_4's generators, copied and rebound.

### Added

- `docs/make_system_map_v0_5.py` -> `system_map_v0_5_numbers.tex`
  (`--check`), `docs/system_map_v0_5.tex` -> `.pdf`: v0_4's map with the
  5s <-> 5p cycle on the measured cross sections (sig.76 at 6 eV, the
  rate's ratio to Kannari's fit at the pulse's Te, in the footer) and
  om.v04, the share with the fit in K76's place.
- `docs/make_overview_v0_5.py` (`--fast`, `--slow`, `--unc`, `--check`)
  -> `overview_v0_5_fast/slow/unc.tex` in a namespace of its own (`\wg`),
  so v0_4's generated numbers stay beside v0_5's as `\wf` in the
  overview: Tahoma's variants now include Kannari's fit in K76's place
  and the band's two ends; v0_3 and v0_4 scored through the version-
  agnostic Tahoma script; the solver's population process checked against
  Jung's tables (forward, reverse, detailed balance); Jung's numbers and
  the TALIF lifetime as constants; the v0_4-equivalent reference row
  (k76 = 'kannari'); the Kleikamp laser through `kleikamp_v0_5.py` in
  subprocesses (the optima of both proxies at three pump scales, the F2
  scan, the loss rate and the trimer's effective rate at 2 atm, the
  trimer sensitivity run) with a table.
- `docs/accuracy_overview_v0_5.tex` -> `.pdf`: the bottom line, the score
  table beside v0_4 with the Kleikamp rows, the errors table, the
  propagated uncertainty (k_5s5p now a scale), "What v0_5 changes" (the
  measured cross section, the finding, the TALIF lifetime, what v0_4 added
  and this keeps), Tahoma, the reference point (four columns: v0_5 and
  v0_4 under both closures), the Kleikamp section, the ladder, the new
  terms' checks, where v0_5 is not yet reliable, Appendix A of the two
  changed rows, the references (Jung 2005/2006, the TALIF review,
  Kleikamp 1986, Berger & Seltzer added).
- `v0_5/tests/test_docs.py` rewritten for the two namespaces and v0_5's
  claims (the measured-over-fit ratios, the share's two directions, the
  switch, the Kleikamp optima's ordering, the sensitivity run).

### Verified

- `make_system_map_v0_5.py`: every node balanced, every row drawn (its
  own checks); `system_map_v0_5.pdf` one page, 0 errors, rendered and
  read: the Kr** node, the two net bands (33.7 % pumped, 80.9 % of the
  Kr* formed), the bands box with "Kannari's fit in K76's place 5.79 %,
  v0_4's", the footer's measured cross section (1.14e-14 cm2 at 6 eV,
  0.92x Kannari's fit at the pulse's 5.68 eV), the optics (g0 3.96,
  alpha 0.64 %/cm); `make_system_map_v0_5.py --check` reports the file
  current. The first build ran to two pages (the footer's new clause),
  shortened; k.76.te first used the burn-out mean (6.24 eV), aligned to
  the map's pulse mean.
- The overview generator failed twice late in its fast numbers before
  writing anything, each time on a copy-and-rebind slip: the Tahoma
  scoring script's name (it is version-agnostic and keeps v0_4's name)
  and ML's two-argument wrapper of network.build (the k76 switch is now
  passed only when set), then the Kleikamp script's import (loaded by
  path). Each cost the twenty-minute fast run.
- The overview TeX compiled against numeric stubs of its three v0_5
  files: 0 errors, 10 pages (one bare-math value in the score table
  fixed).
- `make_overview_v0_5.py` (all three files, 61 min: fast 20, the two
  kinetic runs 20 on 2 workers, the propagation 17): the kinetic
  datasets equal v0_2's; Tahoma kinetic 5.99 / 0.83 / 5.16 %/cm; the
  reference point kinetic 5.73 % at 248 nm; the propagation +-38 / 33 /
  37 / 40 % in quadrature on the 248-nm share, g0, alpha0, g0 - alpha0
  (v0_4: 39 / 33 / 37 / 40), every Tahoma measurement within one unit;
  the Kleikamp optima 1.26 atm (1.02-1.48) and >30 mbar, the trimer
  x2.2, the sensitivity run 1.83 atm. The Kleikamp table's rows first
  ended in one backslash (an f-string escape in the generator); fixed,
  the fast file regenerated (20 min), the PDF rebuilt on it (10 pages,
  0 errors, 0 overfull) and test_docs.py rerun (7 passed, 4 skipped).
- `accuracy_overview_v0_5.pdf`: 10 pages, 0 errors, 0 missing glyphs, 0
  overfull lines after the Tahoma table was scaled to the line and the
  score table's unit moved to its row label; every page rendered and
  read against the numbers (the bottom line, both tables, the
  propagation, the four-column reference point, the Kleikamp section,
  the ladder, the references).
- `make_overview_v0_5.py --check unc` first reported the propagation
  file stale: the comparison read the wall-clock keys (*.min, how long
  the runs took), which differ between two runs of the same code. Both
  overview generators (v0_4's and v0_5's) now compare without them; the
  rerun (17 min) reports the file current. `make_system_map_v0_5.py
  --check` reports the map's numbers current.
- `v0_5/tests/test_docs.py`: 7 passed, 4 skipped (the gated
  regenerations), on the generated files, including the prose
  assertions (the measured-over-fit ratios, the share's two directions,
  the switch, the Kleikamp ordering, the sensitivity run, the largest
  errors).

## 2026-10-09 (c): v0_5, the 5s -> 5p excitation measured; a Kr/F2 laser held out

The owner (2026-10-09): new articles under `Kinetic_Modeling/new_articles/`
(Jung et al. PRA 73, 022722 (2006) and PRL 94, 163202 (2005); the TALIF
review Gazeli et al., Plasma 4, 145 (2021); Kleikamp, Peters & Witteman,
IEEE JQE 22, 219 (1986)), "see if it can help the analysis ... make the
new model v0_5 to improve upon v0_4 and the issues such as section 7 of
the accuracy_overview_v0_4.tex". v0_5 = a copy of v0_4 with K76 on the
measurement (`v0_5/README.md`, SOURCES section 8).

### Found

- **Jung et al.'s measured 5s -> 5p cross sections confirm the 1976 fit.**
  The lumped Kr(5s) -> Kr(5p) cross section (the statistical-weight mean
  over the 5s manifold, the J = 1 levels on their same-core metastable's
  set) is 1.14e-14 cm2 at 6 eV; its Maxwellian rate is 1.37x Kannari's
  Table I [76] fit at 1.16 eV, 1.11x at 3 eV and 0.91x at 5.8 eV. The
  partition v0_4 found therefore rests on a measurement now and moves
  little: the reference point's 248-nm share 5.79 -> 5.68 % (two_t,
  electrons at 5.7 eV where the measured rate is 0.91x) and 5.68 -> 5.73 %
  (kinetic, electrons at 3 eV where it is 1.1x; VUV 5.3 -> 4.8, g0 3.73 ->
  3.79 %/cm), Tahoma's g0 x0.92 -> x0.93 (two_t) and x0.93 -> x0.94
  (kinetic: 5.99 / 0.83 / 5.16 %/cm against v0_4's 5.91 / 0.82 / 5.09),
  alpha0 and the net gain within 1 %. The
  first item of the overview's Section 7 is closed except for the J = 1
  proxies, which the band k_5s5p (0.5 to 1.4, now a scale) brackets.
- **A Kr-buffered e-beam KrF laser is not reproduced, and the term is
  the trimer.** Kleikamp et al.'s two-component Kr-F2 laser peaks at ~2
  atm Kr (F2 10 mbar); LIGHT under their drive (250 A/cm2, 30 ns, 4.7
  MW/cm3 at 2 atm from Berger & Seltzer's stopping power and their own
  B_Kr, run at 0.7-1.4x) puts the optimum of the paper's own saturated
  balance (their Eqs. (1)-(3) on LIGHT's formation rate, loss rate and
  absorption) at 1.26 atm (1.02-1.48 over the power spread) and of the
  small-signal net gain at 0.67 atm. At 2 atm LIGHT's KrF(B) loss rate is
  3.3e9 1/s, 2.8e9 of it the trimer at the window's 362 K (1.02e-30
  cm6/s on the -3.32 law), where the laser's own balance deduces 4.6e-31.
  With G5/G5C scaled by 0.45 (a sensitivity run, NOT adopted) the
  saturated optimum moves to 1.83 atm and the effective quenching is
  their 4.57e-31. The F2 optimum (10 mbar measured) is not reproduced
  either: LIGHT's absorption falls with F2 (Kr2F* quenched away faster
  than F2 and F- absorb), both proxies still rising at 30 mbar; the
  F2-side absorbers (F2's compilation 1.2e-20, F-'s 1971 5.6e-18) and
  KrF* + F2 are the terms. The trimer's rate at 350-400 K was the
  largest term of (b)'s propagated uncertainty; this is the first
  device datum that bears on it. Nothing retuned (predictive-only).
- **A proxy for laser output must saturate.** The small-signal net gain
  quenches the whole inversion and peaks at 0.67 atm; a laser clamps the
  inversion at threshold and quenches only that. The paper's balance is
  the honest proxy. The first run summed rate constants instead of
  events per volume (fixed: `Compiled.heavy`, `ne_pow`).
- **The TALIF review adds a lifetime, not a quench.** Its Table 3 gives
  the Kr 5p'[3/2]_2 natural lifetime 34.1 ns (PRAD's 2J+1 mean 29.0 ns;
  Dzierzega gives that level 28.6); no Kr self-quenching is tabulated,
  so KPS stays on Gazeli et al., Phys. Plasmas 28, 043301 (2021).
- **The band's old meaning did not survive the measurement.** k_5s5p was
  a pre-exponential in cm3/s; it is a factor now and `--ladder kannari`
  must not set it (LADDER_KANNARI drops it; Params.k76 = 'kannari' is
  the fit). A K76 row captures no electron: the exciting electron pays
  the gap through the routing (from_e) and stays; the first draft gave
  the pop_xsec row a captured-energy table and a test expected it.

### Added

- `v0_5/` (a copy of v0_4, rebound): `jung2006.py` (Tables I-IV, the
  lumped cross section, Maxwellian/superelastic/captured-energy tables,
  `kannari_rate`, `summary`), `network.py` (K76/K76r on `("pop_xsec",
  "jung2006")` with `frac=(("k_5s5p", "x"),)`, `kannari_k76_rows`,
  `build(k76=...)`, the band as a scale, LADDER_KANNARI without it),
  `model.py` (`Params.k76`, the pop_xsec kind under both closures, the
  kinetic population process from the measured cross section or the
  fit's pseudo one, scaled by the band), `run.py --k76`,
  `tests/test_v0_5.py`, README, SOURCES section 8; `run_v0_5.py`.
- `docs/krf_data/kleikamp_v0_5.py` (--kinetic, --json, --trimer-scale,
  --only-nominal); `tahoma_v0_4.py` now takes the version's own optics;
  `make_energy_ladder.py` VERSIONS + v0_5 (`energy_ladder_v0_5.tex`).
- `LIGHT/README.md` v0_5 paragraph and quick-start lines.

### Verified

- `jung2006.py` against the papers: Table II/III/IV spot values, the
  1s5 and 1s3 sums at 6 eV (106.9e-16, 127.2e-16), the manifold mean
  1.14e-14 cm2 (the metastables alone 0.49 of it), the Maxwellian table
  against a direct 120 001-point trapezoid (2e-3), the reverse against
  g exp(gap/Te) (2e-4 at 1.3 eV), the kinetic population process against
  the tables for a Maxwellian (<= 1.5e-3 at 1 to 5.8 eV), the band
  scaling both.
- `Params.k76 = 'kannari'` reproduces v0_4's reference point to four
  digits (5.787 / 4.53 / 741.7); `run_v0_5.py --k76 kannari` prints 5.79.
- `v0_5/tests`: test_v0_5.py 10 passed (one test rewritten after the
  captured-energy finding above); the fast files (test_v0_5, test_network,
  test_model, test_v0_4, test_xsec, test_two_level, test_ib, test_docs)
  325 passed, 11 skipped (test_docs: no v0_5 overview yet), 1 strict
  expected failure (Hakuta's decays), and the one failure, the reference
  regression pinned to v0_4's default, repinned to v0_5's with v0_4's
  under the switch (passes on its own, 1 passed); 49 min under a machine
  load of 50 to 90 from other sessions' jobs. test_kinetic.py 24 passed,
  1 skipped (the slow convergence run), 21 min under the same load.
- The kinetic reference point (17 min): 425 solutions, 277 in the pulse,
  the field 1.001 to 1.003 of the input in the second half, the balances
  to 7e-6, the ledger to 1e-14. Tahoma under the kinetic EEDF (19 min
  under load): g0 5.99, alpha0 0.83, g0 - alpha0 5.16 %/cm, Te 1.08 eV,
  n_e 4.1e14, within 1 % of the two_t numbers.
- `kleikamp_v0_5.py`: 3 power scales x 9 Kr pressures + 7 F2 pressures
  under two_t, 167 s; the trimer sensitivity run (--trimer-scale 0.45
  --only-nominal). `make_energy_ladder.py v0_5` renders; `--check v0 ..
  v0_4` unchanged.

## 2026-10-09 (b): the overview's errors, propagated uncertainty and references

The owner (2026-10-09): the citations were not at the end of the page,
and the overview should state the overall accuracy with the largest
errors as errors ("within 6.1 % rms", "27 % error") and the uncertainty
propagated from the sources.

### Found

- **The overview had no reference list.** Its ten sources sat in small
  type inside the bottom-line box and the body named two dozen more
  (Gazeli, Dzierzega, Nelson & Setser, Shiu & Biondi, Ellis, Turnbull,
  Hagelaar & Pitchford, Jung, Vriens, Rokni, Mandl, Hyman, Johnson &
  Hunter, ...) with no citation anywhere. Now a numbered reference list
  of 33 entries at the end, every mention cited by key, and
  `test_docs.py` holds every \cite to a \bibitem and back. Pradler &
  Semmelrock's Gaunt fit is cited through Turnbull (its own reference is
  not among the sources held); Vriens 1965 through Kannari's ref. 79.
- **"73 % Kr2F*" read as an accuracy.** The bottom line said "the
  absorption is 73 % Kr2F*", which reads as 73 % of a measurement; the
  errors are now written as errors (model less measurement over the
  measurement, signed; rms over a set) in a table of their own, the
  largest named, and the Kr2F* sentence says "Kr2F* makes 73 % of that
  absorption".
- **The propagated uncertainty, first numbers.** See the overview's
  Section 1 and `overview_v0_4_unc.tex` (regenerated numbers; this entry
  does not copy them). Inputs: 17 bands at their two ends, 3 measured
  constants at their quoted errors (Kr2+ + e +-10 %, Shiu & Biondi; Kr+ +
  2Kr +-9 %, Grossl; the Kr2F* lifetime +-9 %, Xu), 2 cross sections
  analytically (sigma_se +-12 %, Rokni/Tellinghuisen; sigma(Kr2F*)
  +-20 %, Hakuta). The trimer's temperature law is Params.trimer_tg_pow,
  not a `rates` band: the first run died on it ("rates
  ['trimer_tg_pow']: no row reads them"); `_unc_kw` routes it.
- **At Tahoma the model's own uncertainty is larger than its errors.**
  In quadrature the propagation gives about a third to two fifths of
  each output (the 248-nm share, g0, alpha0, g0 - alpha0), and every
  Tahoma measurement lies within one quadrature unit of the model; the
  trimer's (Tg/300)^n law alone carries most of the gain's term (it has
  no measurement), the Kr2F* quenches and cross section the absorption's.
  The first draft's "envelope" (the linear sum of every input's
  excursion) went below 0 on the 248-nm share, so it is printed as the
  linear ceiling on the spread and not as a bound; the first draft's
  sentence "the absorption's error is not covered" was wrong (1.0
  quadrature units) and is gone.

### Added

- `Params.row_scale` (v0_4/model.py): ((row id, factor), ...) multiplying
  one row's compiled rate, every kind, validated (unknown id, negative
  factor). The sensitivity hook of the propagation; never a fit.
  `test_v0_4.py::test_row_scale_multiplies_one_rows_rate_and_nothing_else`.
- `make_overview_v0_4.py --unc` / `--check unc` -> `overview_v0_4_unc.tex`
  (UNC_ROWS, UNC_SIGMA, UNC_LABELS, `unc_numbers`, the `\wfUncTable`
  table; the Tahoma ends on 2 processes). `--fast`, `--slow`, `--unc`
  select; no flag writes all three.
- The overview: Section 1 "Overall accuracy and its uncertainty" (the
  errors table, the propagation's method, numbers and table), the bottom
  line rewritten around the errors and the propagated uncertainty, the
  `\errof` macro (a model-over-measurement factor as a signed percent),
  `\cite` keys throughout, the reference list.

### Verified

- `make_overview_v0_4.py --unc`: 40 reference-point runs (6 s each) and
  41 Tahoma runs under v0_2's closure (59 s on 2 processes), 6 min in
  all; `overview_v0_4_unc.tex` written, no stub left.
- `accuracy_overview_v0_4.pdf` rebuilt: 10 pages, 0 errors, 0 missing
  glyphs, one 0.5-pt overfull line; pages 1, 3, 4 and 10 rendered and
  read (the bottom line, the errors table, the propagation and its table,
  the reference list of 33 entries).
- `v0_4/tests/test_docs.py`: 7 passed, 4 skipped (the gated
  regenerations), including the new citation check and the propagation's
  assertions (quadrature <= linear, the default inside the excursions,
  every Tahoma deviation within one quadrature unit, the largest term a
  band, a third to two fifths of each output).
- `test_v0_4.py::test_row_scale_multiplies_one_rows_rate_and_nothing_else`
  passes; the four fast test files rerun after `Params.row_scale`:
  293 passed, 1 xfailed in 473.07s (0:07:53).

## 2026-10-09: v0_4's accuracy overview, system map and ladder

The owner (2026-10-09) asked for an accuracy overview and a system map
for v0_4. `docs/accuracy_overview_v0_4.pdf`, `docs/system_map_v0_4.pdf`
and `docs/energy_ladder_v0_4.tex`, in the repo's way: no model number
typed by hand, each written by a generator, and `v0_4/tests/test_docs.py`
holding the files to the model and the words to the numbers.

### Found

- **The wave's electron-ion part has two definitions, and the README
  and the test had the wrong one.** The overview's own check compared
  the solver's field term with the kernel for a Maxwellian: atoms to
  0.7 %, ions 33 % low, and the gap did not move from 600 to 4800 cells
  (not a grid effect). The kernel's printed "ei" is nu_ei's share of the
  numerator; the solver's, with the ions removed from the momentum
  transfer, is the marginal part (the power with the ions less the power
  without them). In 1 atm Kr at 95 GHz nu_en = 4.4 omega (the wave is
  collisional) and the share is 1.49x the marginal part; against the
  marginal part the solver agrees to 0.1 %. The README's "to within its
  grid" was wrong and the test's 50 % tolerance hid it. Fixed: the test
  compares with the marginal part at 1 %, and the Maxwellian closures'
  ib_scale band now acts on the marginal share, as the kinetic one's
  always did (nothing changes at ib_scale = 1; at the reference point
  the band's two ends now span 0.05 % of E_em, 0.09 % before). Known-bad
  line above.
- **The 1-D wall number the first draft cited was not a wall loss.** The
  generator carried a hand-typed 0.008 % as "1Dv0_3's resolved walls";
  1Dv0_3's README gives that as the field work the ions take within 10 um
  of the walls, and 1Dv0_3 prints no walls line at all. The overview now
  reads 1Dv0_2's walls line (0.036 % on 256 uniform cells) from its
  generated file and says that line also holds the transport step's
  mean-energy cut, so no 1-D number checks the 0-D closure yet. The 0-D
  ambipolar loss in a 1-cm gap is 9.4e-6 % of the input (63 1/s at
  mid-pulse), 0.069 % in a 0.1-cm gap at 0.1 atm.
- **What moves Tahoma's gain under v0_2's closure.** Each band alone from
  its default: the trimer's temperature law x0.44, the ion-ion KrF*
  yield at 0.5 x0.69, KrF(B,C) + F at Kannari's analogue x0.77, the
  ion-ion C share at 0.5 x0.84, Kr2* + F at the top x1.06, the cycle off
  x0.97; the rung's own bands (Kannari's values, Kr** + Kr at its low
  end) under 1 %. The kinetic EEDF moves g0, alpha0 and g0 - alpha0 by
  about 1 % from the two_t closure (5.91/0.82/5.09 against 5.86/0.82/
  5.04 %/cm; measured 6.39/1.11/5.28). Tahoma's absorption, scored for the
  first time, is x0.74 under both; at the compilation Kr2F* cross
  section it would be x3.43.
- **The reference point under the kinetic EEDF.** 248 nm 5.68 % of E_dep
  (two_t 5.8, v0_3 kinetic 4.79), Kr* formed 65.4 %, 78.4 % of it pumped
  to the rung net, 96.1 % of the rung relaxed on Kr, 10 min a run, 414
  solutions, balances to 8e-6.
- **The ladder generator's band parser assumed a wavelength label.**
  PRAD's "5p-5s (near IR)" is not "nnn nm"; `_band` now tolerates it.
- **The map needed its own conventions for a cycle.** K76/K76r circulate
  0.8 times the Kr* formed, so the map draws them as two net bands (the
  33.8 % the electrons pump net, the 4.2 % the rung returns to Kr*) and
  says so in the footer; a Kr** node of more than three lines overlapped
  C1 and the F2 rows; the footer needed \tiny to keep the map to one
  page.

### Added

- `docs/accuracy_overview_v0_4.tex` (inputs v0_2's, v0_3's, the map's and
  its own generated numbers plus the ladder; the score table with the
  absorption rows, the 5p rung, the Tahoma variants table, the reference
  point under both closures, the bands table with g0 and alpha, the
  wave's two definitions, the ladder, the reliability list, Appendix A of
  the new rows), `docs/make_overview_v0_4.py` (--fast, --slow, --check
  fast|slow; 2 workers for the kinetic runs; the kinetic datasets held
  equal to v0_2's), `docs/system_map_v0_4.tex` + `make_system_map_v0_4.py`
  (--check added; the Kr** node, the two net bands, the ion-ion limbs
  box, the 17-band table with g0/alpha, the optics footer),
  `docs/energy_ladder_v0_4.tex` (make_energy_ladder.py: the v0_4 entry,
  Kr**, the 5p rows), `v0_4/tests/test_docs.py`.
- `LIGHT/README.md` v0_4 docs paragraph and quick-start lines;
  `v0_4/README.md` item 5 corrected, Checks and Files extended.

### Verified

- `make_overview_v0_4.py --fast` and `--slow` regenerated after the
  fixes; `make_system_map_v0_4.py` (every node balanced, every row drawn);
  `make_energy_ladder.py --check v0 v0_1 v0_2 v0_3` unchanged; both PDFs
  build with latexmk, 0 errors, the map on one page, the overview's
  pages rendered and read.
- After the band change: `test_v0_4.py`, `test_model.py`,
  `test_network.py` and `test_ib.py` 292 passed, 1 strict expected
  failure (Hakuta's decays), 0 failed; the four wave tests against the
  tightened 1 % check among them. `test_kinetic.py` not rerun: the change
  is in `_em_maxwell` alone, which the kinetic path uses only as the
  table's between-points correction (`pem + pemM(te) - pemM(end)`),
  unchanged at ib_scale = 1.
- `v0_4/tests/test_docs.py`: 6 passed, 3 skipped (the gated
  regenerations: `--check fast` is the 15-min run just made, `--check
  slow` the 12-min one, and `make_system_map_v0_4.py --check` was run by
  hand and reports the numbers file current).
- `accuracy_overview_v0_4.pdf`: 8 pages, 0 errors, 0 missing glyphs, one
  0.5-pt overfull line; every page rendered and read (the score table on
  its own page 2, the bands table scaled to the line on page 5, the
  ladder on page 7, Appendix A on page 8).

## 2026-10-08: v0_4, the 5p rung, the optics, KrF* + F, ion-ion limbs, the walls, the wave

The owner (2026-10-08): "Implement resolutions to the bugs and missing
components to make v0_4. Add the absorbers, quenching, ladder, use
Turnbull's IB equations (note that it has a 10% error), add ion drift, and
attempt ion-ion recombination." v0_4 is the owner's copy of v0_3 with each
of these in (`v0_4/README.md`), on the audit of the same day that ranked
them (the audit is in the session's memory notes, not in the tree).

### Found

- **Kr2F*'s 248-nm cross section was the largest error of the absorption.**
  Kannari's Table I carries Johnson & Hunter's compilation value 6.0e-18
  cm2, and so does REBUILD. With it v0_4's alpha0 at Tahoma was 3.8 %/cm,
  3.4 times the measured 1.11 (Spectra Technology 1989 Table 3-3). Hakuta,
  Komori, Mukai & Takuma, JAP 61, 2113 (1987) measured (1.0 +- 0.2)e-18 by
  fluorescence suppression under KrF light and wrote that it is 'about 1/5
  of the value assumed in the current kinetic model'. With the measurement
  alpha0 is 0.82 %/cm (x0.74), Kr2F* 73 % of it, and g0 - alpha0 x0.95
  (two_t; kinetic: g0 x0.93, alpha0 x0.74, g0 - alpha0 x0.96, against
  v0_3's kinetic g0 x0.88).
- **A 5p rung that cannot relax radiates a sixth of the input.** The first
  build had Kr(5p) + Kr as a band at 0 (no Kr measurement in hand): 17 %
  of the input left as near-IR 5p -> 5s light. Gazeli et al., Phys.
  Plasmas 28, 043301 (2021) measured (1.7 +- 0.2)e-10 cm3/s for one 5p'
  level in Kr; with it the relaxation outruns the 29-ns radiation 120 to 1
  and the IR band is 0.3 %.
- **The 5s <-> 5p electron cycle has to be in the Boltzmann equation.**
  With K76/K76r as LIGHT-side rates the kinetic energy balance missed
  LIGHT's by 0.7 of the input (the DC field carried half the power, the
  cycle took the rest unseen). Kannari's A Te^0.5 exp(-C/Te) is exactly
  the Maxwellian rate of sigma(eps) = sigma0 (1 - C/eps) above C, so the
  pair is now a population process on Kr* with that pseudo cross section
  (`boltzmann.PopProcess`), its reverse the Klein-Rosseland superelastic;
  under a Maxwellian it gives the printed law to 0.1 % at 1 eV on the
  solver's grid. The first build scaled the pair's reverse by the Kr*
  fraction as well; the reverse acts on Kr**, seven times smaller at the
  reference point, so the solver returned seven times too much energy to
  the electrons and the gap only fell to 0.4. With each operator on its
  own population the field carries 1.001 to 1.004 of the input through
  the pulse and the two energy balances agree to 7.9e-6 (the reference
  run; 2.5e-7 over a 60-ns pulse). The kinetic reference point gives
  5.68 % at 248 nm (two_t 5.79, v0_3 4.79), VUV 5.29 %, heat 71.9 %.
- **The wall term must not read Te = U_e/(1.5 n_e).** In the attaching
  afterglow n_e is a small difference of two large ion densities; the
  ambipolar loss written with that Te has a Jacobian that explodes as
  n_e falls (its value stayed small), and BDF stalled re-evaluating it: a
  0.3-us run with a gap took longer than 10 minutes where 5 s is normal.
  The flow is now driven by U_e/(1.5 sum n+), the electron energy over
  the cations, smooth and capped at 30 eV (`Model.wall_te`).
- **scipy's event root-finder can fail inside a kinetic run.** In the
  0.3-us run of the kinetic tests, in the attaching afterglow where n_e
  falls to its floor, solve_ivp saw a trigger's sign change between a
  step's ends that its dense output could not bracket and raised
  "f(a) and f(b) must have different signs". Raising the Kr* and Kr**
  fraction floors to 1e-7 (the trigger had fired through the first
  nanoseconds while Kr** grew from nothing) was not enough.
  `Model._run_kinetic` now catches that error, takes a short step without
  the triggers (1e-3 of the renewal cap) and renews the table, as the
  stall guard does for a trigger crossed at a segment's start
  (2026-10-09).
- **Particles lost to the walls leave at Tg and change it not at all.**
  The reaction rule (new particles take cv k Tg from the heat) applied to
  the wall flux raised Tg spuriously; the lost particles' thermal energy
  is booked as E_wall_th and excluded from the Tg equation. Booking it a
  second time as a walls fate line double counted it (a residual of
  E_wall_th / E_dep, 2.3e-6 in a thin gap at 0.1 atm): the fate's walls
  line is E_wall alone.
- **The 5s <-> 5p cycle decides the partition at the reference point.** On
  Kannari's fit 81 % of the Kr* formed is pumped to 5p and relaxed to
  heat; the Kr2* VUV falls from 20 to 4.5 % and the 248-nm share rises
  from 4.67 to 5.79 % (two_t). The fit (Jacob & Mangano 1976) has no
  direct measurement behind it; Jung et al. 2005/2006 measured the cross
  sections out of the metastable but the absolute values are paywalled, so
  the band runs to 0 (-1.1 points).
- **The IB model's error band must act on the absorbed power, not inside
  the collision frequency**: scaling nu_ei in nu/(nu^2 + omega^2) is
  linear in the e-i share only for omega >> nu, and misses by 10 % of
  itself at 1 THz in 1 atm Kr (ib.py's tests). ib_scale scales the
  electron-ion share of the power under both closures.
- **Turnbull's 10 %, in their words.** PRL 130, 145103 (2023) p.145103-4:
  with screening 'Sommerfeld just 2% above' the data on average over 22
  shots; Phys. Plasmas 31, 063304 (2024) Sec. IV: the mixture shots with
  absorption > 1.5 % 'all predicted to better than 10%', and Sec. I: the
  models 'overestimated absorption by about 10% unless [screening] was
  included'. The 10 % is the per-shot scatter of the best data.
- **The Langevin limb from the measured mobilities is 30 % below
  REBUILD's** (8.2e-6 cm3/s amagat for both pairs against 5.66e-6 and
  6.35e-6); the ion-ion turnover at 1 amagat falls 10 to 14 %, the 248-nm
  share at the reference point by 0.03 %.
- **No measurement exists of KrF* + F**, nor of any RgX* + X analogue
  (Jancaitis 1984, Jacob 1978, Rokni 1978, Ong 1996 carry none or
  Kannari's 'analogous to F2'); nor of the electron-impact ionization of
  metastable Kr (Vriens 1965, Ton-That & Flannery 1977, Hyman 1979 and
  Margreiter 1990 are calculations). Both stay bands at 0.
- **At 1 THz a kilovolt per centimetre heats nothing**: 3 kV/cm rms
  absorbs 1.7e-5 J/cm3 in a microsecond at 1 atm (P ~ nu/omega^2).

### Added

- `v0_4/`: v0_3's copy with Kr** (5p, 11.30 eV, the Biagi group's own
  threshold) and the rows KPS, PRAD, K76, K76r, KSTEP_S, KSTEP_P, K27,
  K27C, K27X, QF, QFC, K26C, K26X, K108C, K108X; the rate kinds te_fit
  (Kannari's form, with its detailed-balance reverse) and ionion
  (Flannery's turnover with REBUILD's three-body limb and the Langevin
  limb of the Ellis mobilities, ION_K0); Row.frac for band fractions; the
  bands krf_f, kps_kr, k_5s5p, kstep_s, kstep_p, k5p_f2, ii_yield, ii_f_c
  and the free band ib_scale; Params.ladder ('kannari' = every ladder band
  at Kannari's value), gap_cm, em_E_V_cm, em_f_Hz; the ledger slots E_em,
  E_wall, E_ion_heat; `optics.py` (g0, alpha by species, Kannari's block
  with Hakuta's Kr2F*); `ib.py` (Turnbull's model and the Lorentz kernel,
  SI inside, 219 checks against the IB project's implementations);
  boltzmann's PopProcess, the ions' Debye-screened Coulomb collisions in
  the DC term, the high-frequency field term with Turnbull's per-energy
  cross section, em_power; the report's optics block, walls line and EM
  input; run.py's --ladder, --gap, --em-E, --em-f; `run_v0_4.py`.
- `docs/krf_data/tahoma_v0_4.py`: Tahoma g0, alpha0 by species and the
  pure-Kr loss for any LIGHT version under the beam drive (ML's
  configuration), two_t or kinetic.
- `v0_4/SOURCES.md` section 7: every new value with how it was read.
- Tests: `v0_4/tests/test_v0_4.py`, `test_ib.py`, the ladder and limb
  tests in `test_network.py`, the rate-law tests in `test_model.py`; the
  three inherited tests whose expectations the new rows move
  (reference-point regression, two_t reproduces v0_2, the flux split's
  five parts) rewritten to the neutralised configuration or the new
  signatures.

### Changed (v0_4 only; v0_3 untouched)

- C1a-c make Kr** (the cascade stops at 5p); R11p makes Kr** + Kr (Shiu &
  Biondi's 5p lines), 1.546 eV of heat per event where v0_3 booked 2.93.
- K26 and K108 split into B, C and dissociation rows through the ii bands.
- The kinetic state carries n_ion, the wave (omega, en_em_td) and the Kr*
  population; the table carries the wave's power (pem) and reads it as
  it reads the elastic power, carried past its ends by the Maxwellian
  kernel's change.

### Verified

- `/usr/bin/python3 -m pytest LIGHT/v0_4/tests` (2026-10-09): the fast
  files (test_v0_4, test_network, test_ib, test_xsec, test_two_level,
  test_model) 315 passed, 1 strict expected failure (Hakuta's decays, as
  in v0_3); test_kinetic.py 24 passed, 1 skipped (the slow convergence
  run, LIGHT_KINETIC_SLOW). Run the kinetic file alone: three kinetic
  jobs at once made each run 4 to 10 times slower.
- The ledger closes to 1e-13 in every new mode; the neutralised rung
  reproduces v0_3's reference point to 0.03 % (248 nm) and the two_t
  comparison with v0_2 holds the emission columns to 0.3 % and the ion
  flows to 1 %.
- The kinetic reference point: field/input power 1.001 to 1.004, the two
  energy balances to 7.9e-6; under the beam (Tahoma, 30 ns) the kinetic
  mean is LIGHT's Te to 1.9e-4 after 10 ns and 1e-9 after 17 ns, the
  balance 1.2e-3 falling to 5.4e-4.
- ib.py against the IB project's own Turnbull, qib_kinetics and
  nlib_kinetics implementations: 219 checks to 1e-10 (the Drude limit to
  8e-10, v0's kernel to 4e-7).
- Tahoma scored, never fitted: two_t g0 x0.92, alpha0 x0.74, g0 - alpha0
  x0.95; kinetic x0.93, x0.74, x0.96 (`docs/krf_data/tahoma_v0_4.py`).
- No accuracy overview or system map yet for v0_4 (docs/ generators are
  v0_3's); the README carries the numbers.

### Known-bad additions (also at the top)

See the list: the Langdon factor under the kinetic closure, ladder
electron rows outside the Boltzmann equation, an error band inside a
collision frequency, Kr2F*'s compilation cross section, a 5p rung with
no relaxation.

## 2026-10-06: 1Dv0_3 keeps the field's work to the field's energy

The owner asked to proceed with the conservation of energy in the field's
coupling, the next physics the 2026-10-01 (c) entry named. In 1Dv0_3 the
electric field's work now moves energy from the electrons through the
field to the ions, with its sign, and none of it is input.

### Found
- **Kept to its energy, the field's work settles the 64-cell walls.** At
  the reference point, the gas moving, percent of the given power:

  | | 64 cells | 128 cells | 256 cells |
  |---|---|---|---|
  | 248 nm, the gap, dt 50 ps | 4.642 | 4.632 | 4.618 |
  | the same at dt 25 ps | 4.643 | 4.632 | |
  | with Fluid1D's field work, dt 50 ps | 4.662 | 4.637 | 4.616 |
  | with Fluid1D's field work, dt 25 ps | 4.661 | 4.645 | |
  | the field's work on the ions | 0.498 | 0.486 | 0.285 |
  | the walls' hottest gas, K | 1234 | 1225 | 1180 |
  | with Fluid1D's field work, K | 2810 | 1375 | 1182 |

  The centre cell is v0_2's 0-D run in every run. 64 cells now put the two
  nearest bands at 24 and 80 um (1230 and 1080 K), where 256 put them at
  25 and 69 um (1180 and 1170 K). On every mesh the walls hold 10 to 11 %
  more pressure than the centre (51 % on 64 cells before) and the gas
  moves at 15 m/s at most (73 before), and the gap's share converges in dt
  to 0.02 % (0.2 before). With the gas still the gap's share is 4.621,
  4.608 and 4.593 % on 64, 128 and 256 cells. Root cause of the old
  64-cell layer at 2000 to 2800 K: the step's error, given as input, which
  the expanding gas fed (2026-10-01 (c)).
- **The field's work sits at the bands' edges, not in the sheath.** Within
  10 um of the walls the ions take 0.007 to 0.010 % of the given power on
  every graded mesh. The rest is where n_e jumps between neighbouring
  cells, up to 160, 45 and 20 times on 64, 128 and 256 cells, and it falls
  as the edges are resolved (0.50, 0.49 and 0.29 % with the gas moving,
  0.37, 0.28 and 0.21 % with it still). Over a run the electrons lose and
  the ions gain in every cell.
- **Fluid1D's density floor makes charge.** Its species update sets a
  density that goes negative to zero. With the work paired to the exact
  fields, the theorem missed by 44 % of the work on 48 cells from 10 um,
  whose wall cells' electrons went to -2.4e16 m^-3 before the floor, and
  held to 1e-10 with the floor removed. Root cause: the floor adds
  particles, so the charge the field sees is not the charge the currents
  moved. No reference run needed an outflow scaled, so there the species
  update is Fluid1D's (the tests' 48 cells from 10 um scale 303 in 40 ns).
- **The uniform mesh's wall cells hold the whole sheath.** There the
  sheath's work on the electrons is more than they hold, and the floor puts
  back 2.3 % of the given power on 64 uniform cells and 0.4 % on 128 (none
  on the graded meshes). The stage then takes up to Fluid1D's 8000
  substeps a step, since its substep follows the given power over the
  electrons' stored energy, so 128 uniform cells take 23 min against 100 s
  with Fluid1D's field work. On 64 uniform cells Fluid1D's step had
  dropped 16 % of the given power as negative work.
- **The moving gas changes the field's energy.** Carrying charge, it
  changed it by -2.5e-6 of the input over the default run, where the
  field's work on all charges is -2.5e-6 and the field ends with 1e-8. The
  field's books close only with that term, to 3e-13 of the electrons' work
  with it and 5e-4 without (2.9e-11 and 2.8e-4 on the tests' short run).

### Added
- `field_work` in `Params1D` and `MeshPlasmaModel` ('conserve', the
  default, or 'fluid1d'), and `--field-work`. In `transport.py`
  `_exact_field`, `_field_energy` and `_face_work`, the positive transport
  (`positive` in `_transport_species`), and closure mode 2, in which the
  electrons' energy takes their share of the faces' work in the step. The
  books carry `field_e` and `field_i` (the ions' heat, its own line in the
  fate) and the field's energy, and the results `field_energy`,
  `field_defect`, `field_moved` and `positivity_limits`. The summary
  reports the field's books. `field_work="conserve"` refuses a bias, a
  drive current and a gap capacitance, and books the transport step's own
  mean-energy bounds with the cut, where Fluid1D's mode keeps them in the
  walls' line as 1Dv0_2 does.
- Tests. In `test_transport.py` (12), Gauss's law and the field's energy,
  one conserving advance from a live state against the theorem, twin
  advances against Fluid1D's step, and the positive transport against the
  floor. In `test_1d.py`, the conserving books on 48 graded cells with the
  gas still and moving, and `--field-work`.

### Changed
- The default is `field_work="conserve"`. The parity with 1Dv0_2 needs
  `field_work="fluid1d"` beside the uniform mesh and the still gas, and
  the tests' short runs pin it.
- `test_reference_point_on_the_graded_mesh` (slow) holds the conserving
  books.
- 1Dv0_3's version is 0.1.2+1d3.
- `1Dv0_3/README.md` gives the conserving runs, the field's work, the
  walls with it and the uniform mesh's cost, and `LIGHT/README.md` says so
  in brief.

### Fixed
- `1Dv0_3/README.md` said the still gas's four bands on 256 cells peak at
  1850 to 1950 K. The two nearest the walls do (1848 and 1950 K), and the
  two further in peak at 1708 and 1426 K in the same run's saved profiles.
  The range was the two nearest bands' peaks.

### Verified
- The theorem, re-derived from Gauss's law and continuity (the comment
  above `_exact_field`). On a live 48-cell graded state one advance's work
  matches the field's energy change to 1e-8 of the work's scale, the cells'
  shares sum to the faces' total to 1e-12, and in the reference runs the
  residual summed over 2 us is 3e-11 % of the given power at most.
- The twin advances. The electrons' energy differs by their share to
  2.3e-13 and the densities agree to the bit.
- A mutation check, 26 bugs planted one at a time in `transport.py`,
  `sim1d.py` and `run.py`, all caught. The first pass let one through (the
  ions' share booked but left out of the gas) until the walls' loss was
  held below Fluid1D's, and the twin advances and the wall cell's content
  were written for two that no check then covered (a limited cell giving
  nothing, the electrons' step leaving out their share).
- The reference runs. 64, 128 and 256 graded cells with the gas moving (25
  ps too on 64 and 128) and still, 64 and 128 uniform cells, beside the
  2026-10-01 (b) and (c) runs with Fluid1D's field work.
- From the main checkout, on 2026-10-06. 1Dv0_3's tests under
  `~/.venvs/f2sweep/bin/python` 139 passed, 2 skipped (the slow two) and 1
  xfailed (Hakuta, as in v0_2), in 212 s, and with `LIGHT_1D_SLOW=1` the
  slow two pass (the default's reference point in 693 s, 32 uniform cells
  in 131 s). Under `/usr/bin/python3` 97 passed, the 1-D modules skipping.
  1Dv0_2's `test_1d.py` still passes (22 passed, 1 skipped). The default
  `LIGHT/run_1Dv0_3.py` from main gives the 128-cell reference run's
  numbers (4.632 %, the ions 0.486 %).

## 2026-10-01 (c): 1Dv0_3's neutral gas moves

1Dv0_3's resolved walls asked for moving neutrals (2026-10-01 (b)), and
the owner said to proceed with the next step. `neutrals.py` now lets the
gas flow, conduct heat into the electrodes and its F2 and F diffuse.

### Found
- **The moving gas does not settle the walls.** At the reference point,
  percent of the given power:

  | 248 nm, the gap | 64 cells | 128 cells | 256 cells |
  |---|---|---|---|
  | the gas still, dt 50 ps | 4.637 | 4.608 | 4.591 |
  | the gas still, dt 25 ps | 4.645 | 4.607 | |
  | the gas moving, dt 50 ps | 4.662 | 4.637 | 4.616 |
  | the gas moving, dt 25 ps | 4.661 | 4.645 | |

  The centre cell is v0_2's 0-D run (4.667 %) in every run. On 256 cells
  the moving gas keeps the still gas's four bands at each wall, at 25, 69
  to 73, 119 to 125 and 163 to 172 um, but cooler, peaking near 1.2 us at
  1030 to 1170 K against the still gas's 1850 to 1950 K, with n_e 11 to 19
  times the centre's. 128 cells resolve only the two nearest each wall (26
  and 75 um, 1170 and 1360 K). On 64 cells one layer 15 to 65 um from each
  wall takes 3 to 4 times its given power from the field, reaches 2000 to
  2800 K and 290 kPa against the centre's 211, and drives a compression
  wave into the gap at up to 73 m/s. Root cause: the field work's error
  (next item), which the moving gas feeds, since as the hot gas expands N
  falls and E/N and the electrons' mobility rise.
- **The field's work in the grounded gap is the coupled step's error.**
  Nothing outside drives the gap, so from a start with no field the
  field's net work on all charges is minus its energy at the end (9e-7 of
  the input at 600 ns). A copy of the kernel that also booked the work on
  every charge measured, to 600 ns on 64 cells, +0.297 % of the given power
  with the gas still (0.260 on the electrons, 0.037 on the ions) and
  +1.877 % with it moving, and the faces' products give +0.69 %. Root
  cause, most likely: the semi-implicit step solves each substep's field
  from a predicted charge, and the explicit species update and its floors
  then move the charge another way, so the field and the current never
  pair exactly. The electrons' share falls as the walls are resolved
  (1.09, 0.033 and 5e-5 % with the gas still, 3.08, 0.36 and 1.2e-4 %
  moving) and grows as dt falls (2.58 % on 64 cells and 0.095 % on 128 at
  25 ps, with the gas still).
- **The (b) entry's dt reading was the field work in the denominator.**
  Its "4.587 to 4.528 %" on 64 cells is the 248-nm share of the given power
  plus the field work, the books' total that `fate()` and the summary's
  table divide by. Per the given power, which `share_248()` and
  `--compare` use, it is 4.637 to 4.645 %, and every graded mesh converges
  in dt to 0.2 %. The field work went from 1.09 to 2.58 % there.
- **What heats the still gas's bands is the charges' energy.** The hottest
  band cells on 128 and 256 cells receive no field work, yet their gas
  takes 2.4 to 2.8 times their given power as heat. Only the charges move
  in the still gas, so this is the ionization energy of the ions and
  electrons that drift and diffuse into a burned band and recombine there.
- **The gas carries Fluid1D's electron energy density, not its mean
  energy.** Fluid1D blends the electrons' mean energy toward the local
  field's where they are few (`_recover_eps`), so it differs there from
  the energy the books hold. Carried with the gas, the mean energy left
  2e-7 % of the input unbooked on the short run. The gas now carries the
  energy density per electron (2e-11 %).
- **Svehla's eq. (B1) misprints its coefficient.** It prints 198.01 where
  Brokaw's eq. (6) and (15/4) R x 26.693 give 198.91, and his tables follow
  the relation (Kr at 300 K, 22.5e-6 cal/(cm s K)).

### Added
- `1Dv0_3/neutrals.py`, the gas after each transport step. The 1-D
  Navier-Stokes equations of the mixture on the mesh, every species and the
  electrons' energy riding with it, F2 and F diffusing in the krypton
  (whose flux balances theirs, so diffusion carries no mass), heat
  conducting into the electrodes at tg0, viscosity. Chapman-Enskog's
  coefficients on Lennard-Jones potentials with Svehla's force constants
  (Kr 3.655 A and 178.9 K, F2 3.357 A and 112.6 K, F an estimate),
  Brokaw's diffusion and combining rules, the collision integrals of
  Neufeld, Janzen and Aziz. Finite volumes, MUSCL-Hancock, HLLC, explicit
  and conservative, substeps under the sound and diffusion limits.
- `Params1D.neutrals` ('move' or 'still') and `--neutrals`. The books carry
  the heat conducted into the electrodes (`wall_heat`, in the two wall
  cells) and the gap's energy the gas's kinetic energy. The profiles carry
  the gas's velocity and pressure, the results its substeps, and the
  summary a line for the gas and the 248-nm share of the given power alone.
- Svehla 1962 (NASA TR R-132) and Brokaw 1961 (NASA TR R-81) in
  `docs/krf_data/sources/`, from NTRS, checktext CLEAN, and `SOURCES.md`
  section 7.
- Tests. `test_neutrals.py` (18, NumPy and SciPy, any Python) holds the
  coefficients to a quadrature, Svehla's tables and Brokaw's examples, and
  the flow to exact solutions (Sod's shock tube on uniform and graded
  cells, a species front on the contact, a contact at rest, sound and its
  classical decay, heat between adiabatic walls and into the electrodes,
  2Dt spreading, no mass in diffusion, no heat in mixing, a closed box's
  conservation, the electrodes' heat in the books, the substep's limits, a
  uniform gas to the bit). `test_1d.py` holds the coupled gas to the books
  and the interior, and the command line's `--neutrals`.

### Changed
- The default is the gas moving on 128 cells from 0.9508138901170884 um,
  the 64-cell mesh refined once (13 min at the reference point). On 64
  cells the moving gas's walls are the field work's (3.08 % of the given
  power, 5.32 % at 25 ps), and 64 cells stay the quick mesh
  (`--cells 64 --wall-cell 2`). `neutrals="still"` is 1Dv0_2's gas, and on
  the uniform mesh with it 1Dv0_3 is still 1Dv0_2 to the bit.
- 1Dv0_3's version is 0.1.2+1d2. The gas solver is `self.flow`, since
  Fluid1D's Simulation keeps its gas mixture in `self.gas`.
- `1Dv0_3/README.md` gives the shares of the given power, the moving gas,
  and the field's work.

### Verified
- From the main checkout, on 2026-10-05 (the Mac slept with its lid closed
  through the runs of 2026-10-01). 1Dv0_3's tests under
  `~/.venvs/f2sweep/bin/python` 133 passed, 2 skipped (the slow two) and 1
  xfailed (Hakuta, as in v0_2), in 117 s, and with `LIGHT_1D_SLOW=1` the
  slow two pass (the default's reference point in 587 s, 32 uniform cells
  in 80 s). Under `/usr/bin/python3` 97 passed (the 0-D copies, the mesh
  and the gas tests), the 1-D modules skipping. 1Dv0_2's `test_1d.py` still
  passes (22 passed, 1 skipped).
- The sources, read first-hand in the NTRS scans (`SOURCES.md` section 7).
  The collision-integral fits against a direct quadrature of the (12-6)
  definitions (0.2 % over T* = 0.5-50), the viscosities against Svehla's
  Table III (0.9 %), the diffusion against Brokaw's examples (argon
  0.1795 cm^2/s, measured 0.180, H2-N2 1.042, measured 1.052).
- The coupled books. The closure residual stays below 1e-10 % of the input
  in every reference run, and the short run's centre cell is the still
  gas's to 1e-9.
- The reference runs at 50 ps (64, 128 and 256 cells) and 25 ps (64 and
  128), the gas moving, beside 2026-10-01 (b)'s with it still. The field's
  work on every charge from an instrumented copy of the kernel (scratch,
  not landed).
- A mutation check planted 29 bugs in `neutrals.py` and its coupling in
  `sim1d.py`, one at a time, each in a copy of its own. Five tests were
  written first for slips the flow tests would have missed (a downwind
  composition, diffusion carrying mass, an enthalpy flux without its pV
  term, a full cell to the electrode, the substep's limits), each sized
  against its slip by a probe. The check then caught 28 of 29. The 29th
  leaves unbooked what Fluid1D's adoption changes in the energy the gas
  carried, 5e-16 of the input in every run tried.

## 2026-10-01 (b): 1Dv0_3, 1Dv0_2 on a mesh graded by OSIRIS's load rule

The owner made `LIGHT/1Dv0_3/`, a copy of 1Dv0_2, and asked: "In
Kinetic_Modeling/osiris/ the software uses grid scaling for where important
interactions occur. Can we apply that to the new copy of 1Dv0_2 titled
1Dv0_3". 1Dv0_2 and its chemistry stay as they are.

### Found
- **OSIRIS has one cell size.** dx = (xmax - xmin)/nx_p on each axis
  (`osiris/source/os-simulation.f03`, line 330), and it has no stretched,
  mapped or refined grid. Its grid scaling is load balancing. It divides
  the grid among its processors so each holds an equal share of a load (the
  particles in its cells, or a density given as `spatial_loaddensity`), by
  recursive bisection of the running load total (`rec_tree_dist`,
  `os-grid-parallel.f90`, line 236), which on a continuous load cuts at its
  quantiles. Read first-hand, `SOURCES.md` section 6.
- **The walls are where 1Dv0_2's gap needed cells.** At the reference point
  (1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, 2 us, 1-cm gap):

  | percent of the input | uniform 64 | uniform 128 | graded 64 | graded 128 | graded 256 |
  |---|---|---|---|---|---|
  | 248 nm, the gap | 4.521 | 4.594 | 4.587 | 4.606 | 4.591 |
  | cut by the mean-energy bounds | 3.100 | 1.536 | 0.120 | 0.204 | 0.210 |
  | lost to the walls | 0.042 | 0.041 | 0.015 | 0.013 | 0.014 |
  | negative field work kept out | 16 | 3.21 | 0.028 | 0.049 | 0.125 |

  The graded meshes start at 2, 0.95 and 0.46 um and are nested. The centre
  cell is v0_2's 0-D run (4.667 %) on all of them. On a uniform mesh the
  cut is the wall cell's whole input (2/64 of the gap holds 3.1 %), and the
  negative field work Fluid1D's split keeps out is a phantom of the
  one-cell sheath (1572 % of the input in a 20-ns pulse on 48 cells). On
  the graded meshes the cut is the drained sheath's input, the negative
  work is the electrons' real work against the sheath's field (it grows as
  the sheath is resolved, and the split still leaves it out), and the
  positive work the split gives falls away (1.09, 0.03 and 5e-5 %). 64
  graded cells give the gap about what 128 uniform ones give, and the three
  graded meshes agree to 0.4 %.
- **Resolved, the walls show the F2 burning out cell by cell.** Inside the
  drained sheath, from about 10 to 200 um, the cells burn their F2 one at a
  time from the wall inward during the pulse (on 128 cells the cell at
  29 um by 300 ns, at 75 um by 500 ns, at 116 and 159 um by 600 ns), and a
  burned cell stops attaching electrons, its n_e rises to 100 times the
  centre's, and it ends near 2000 K. It is the halogen-depletion
  instability of discharge-pumped excimer lasers, seeded at the wall. Root
  cause of its sharp bands, 2000 K beside 600 K: the 1-D models move no
  neutrals, no F2 diffusion and no gas expansion, which 1Dv0's README
  justified because 1 us at 1 atm spreads them by microns. That holds on
  39-um cells and not on cells of a few microns. On 64 and 128 cells each
  burned region is one cell wide, and on 256 cells they are bands of
  several cells at 28, 73, 119 and 163 um (29, 75, 116 and 159 on 128
  cells), peaking at 1850 to 1950 K, so the structure converges in place.
  On the default mesh, whose cells there are 3 to 25 um, it moves with dt.
  Halving dt to 25 ps takes the gap's 248-nm share from 4.587 to 4.528 %
  and the field work given from 1.09 to 2.58 % of the input (the centre
  cell stays at 4.667 %). On 128 cells 25 ps takes it from 4.606 to
  4.602 %, and the graded meshes and steps agree to 0.3 % without the
  default mesh's 25-ps run. The default mesh's gap totals are good to about
  1 %, 128 cells' to about 0.3 %, and the interior is the 0-D run. 1Dv0_2's
  "cells beside the wall cell do not converge" (2026-09-30) is this, held
  in one or two cells.

### Added
- `1Dv0_3/mesh.py`. The N + 1 edges sit at the k/N quantiles of a load w =
  1/h(d), h = h0 + s d from the nearer wall, so the cells grow by e^s from
  the wall cell. `graded(cells, length, wall_cell, max_cell)` solves s,
  refuses a mesh growing faster than 1.3 a cell, and `refine` doubles the
  cells on the same load (nested meshes). `uniform` is Fluid1D's Mesh1D to
  the bit.
- `1Dv0_3/transport.py`. Fluid1D's coupled transport kernel with the mesh
  in place of its one dx, with spacings at the faces, widths in the
  updates, half-width weights in the face averages (two half-cells in
  series), the field solve's wall faces half a cell out, each cell's own
  substep limit, length-weighted gap means. It calls Fluid1D's
  geometry-free kernel functions and, as a derivative of the kernel,
  carries the AGPL. Fluid1D is not edited. `MeshPlasmaModel` refuses the
  frozen-field NumPy path, which reads one dx, and `MeshPoisson` is the
  same rows in NumPy.
- `Params1D.mesh` ('graded' or 'uniform'), `wall_cell_um`, `max_cell_um`,
  and `--mesh`, `--wall-cell`, `--max-cell`. Every total sums over the
  cells' own volumes, a power profile has a volume-weighted mean of 1, and
  Fluid1D's three books kept with one volume (the transport's closure, the
  optics' resync and the photoelectrons) are redone per cell. On a graded
  mesh Fluid1D's dx and cell volume are unset so any read of them fails.
  Results carry the cells' widths and the transport's substeps.
- `run_1Dv0_3.py`, `1Dv0_3/README.md`, `SOURCES.md` section 6.
- Tests. `test_mesh.py` (9, NumPy only), `test_transport.py` (7) and, in
  `test_1d.py`, the uniform mesh against 1Dv0_2 to the bit, the input over
  the gap's volume, volume-weighted means and profiles, the graded centre
  against the 0-D model, the graded walls against the uniform ones on 48
  cells, the graded command line, the reference point on the default mesh
  (slow), and the Fluid1D functions transport.py follows, pinned by hash.

### Changed
- 1Dv0_3's copied files name 1Dv0_3, its version is 0.1.2+1d1, and its
  default mesh is graded, 64 cells from 2 um at the walls (5 min at the
  reference point against 1Dv0_2's 80 s on 128 uniform cells). 1Dv0_2's
  `tests/test_docs.py` was dropped from the copy, since it tests 1Dv0_2's
  overview.

### Tried and left
- **No given power in the drained sheath.** With none within 10 um of the
  walls the bounds cut nothing and the walls take 0.004 %, but the burned
  cells move outward (2160 K at 31 um).
- **The field's work from the faces.** Fluid1D gives a cell -e Gamma E of
  its mean flux and mean field. The average of its two faces' products,
  whose sum is the field's work between the walls, made the near-wall cells
  run away sooner (the gap's mean n_e 2.8e16 cm^-3 by 1.2 us).
- **The given power shared per electron.** In proportion to each cell's
  n_e, with the same total, the gap contracts at once (its mean n_e
  3.7e16 cm^-3 by 200 ns on 64 cells, 6.7e15 on 128, against 9e13).

### Verified
- From the main checkout, 1Dv0_3's tests under
  `~/.venvs/f2sweep/bin/python` 114 passed, 2 skipped (the slow two) and 1
  xfailed (Hakuta, as in v0_2), in 96 s, and with `LIGHT_1D_SLOW=1` the
  slow two pass (the default mesh's reference point in 227 s, 32 uniform
  cells in 61 s). Under `/usr/bin/python3` 79 passed (the 0-D copies and
  the mesh tests), the 1-D modules skipping. 1Dv0_2's `test_1d.py` still
  passes there (22 passed, 1 skipped).
- Parity. One transport advance from a live state on a uniform mesh is
  Fluid1D's to the bit on 8 and 64 cells (densities, energy, field, ledgers
  and substeps), and the same short run on 8 and 16 uniform cells is
  1Dv0_2's to the bit in every cell (the totals to 1e-12, their summation
  order). The reference point on 64 and 128 uniform cells gives 1Dv0_2's
  published numbers (4.521 and 4.594 % at 248 nm, the bounds' cut 3.100 and
  1.536 %).
- The convergence runs at the reference point. Graded 64, 128 and 256 cells
  (nested) at 50 ps, 64 and 128 at 25 ps, and uniform 64 and 128.
- A mutation check planted 20 bugs in `mesh.py`, `transport.py` and
  `sim1d.py`, one at a time, each in a copy of its own with its own pycache
  and numba cache. The first pass caught 13. Tests were added for six of
  the survivors (the face weights, both substep limits, the Peclet spacing,
  the drift's face weights, the energy divergence), and the second pass
  caught them all. The 20th, a uniform mesh's dx as edges[-1]/N, changes no
  result for any gap or cell count tried.
- Fluid1D is unedited. Its files were only read, and the hashes
  `test_1d.py` pins hold.

## 2026-10-01: the first nanosecond, as a time step

The owner added `LIGHT/improvement_for_EEDF.txt`, three ways to treat the
steady-shape limitation from cheapest to most complete, and asked: "Test if
we can fix the issues regarding the first nanosecond". The note's option 1
(measure first) came first, then its option 2 (a memory term in the
existing solve), with its option 3 (the physical drive) run offline as a
check.

### Found (measured before any change)
- **No steady EEDF holds the heating transient.** At the reference point
  the electrons heat from the gas temperature through the first nanosecond.
  Solved steady, the EEDF dropped the growth term at 86 of the 88 updates
  of the first nanosecond (and at none after it), its field carried 2 to
  19 % of the input power there, and its energy balance missed LIGHT's by
  up to 22 % of it. Root cause: most of the input power raises the mean
  energy, n_e d<eps>/dt, which a steady solution has no term for. With the
  field holding only the losses, the growth term's energy outweighed the
  heating, the term was dropped, and the first cell took up the attached
  electrons cold.
- **The bands do not depend on it.** With LIGHT's own Maxwellian closure in
  the first nanosecond no band moved by more than 6.7e-4 of itself against
  the steady run (248 nm, 4.79714 % against 4.79392 %), and the run took
  160 s against 263 s.
- **The first nanosecond was dear.** Solved steady, the first 2 ns took
  111 s, 110 s of it in 98 kinetic updates. Those were the cold start
  (13 s), 61 recombination-triggered renewals in the first 46 ps, 27
  Te-triggered and 9 Penning-triggered ones. As time steps the first 2 ns
  take 28 s, 27 s of it in 98 updates.

### Added
- `boltzmann.Solver.step`, one backward-Euler step of the time-dependent
  EEDF, the steady problem plus sqrt(eps) (F - F_prev)/(N dt) (a diagonal
  like the growth term's and a source). Its growth term is the population's
  own change, nu from F, with no rule to drop it. With it
  `solve_joint(memory=, target_power=)`, `_power_match` (the field from a
  Joule power, option 3's drive), `Grid.remap` (a history onto another grid
  with its mean kept), `Grid.maxwellian_at_mean` and `_sg_factor`.
- `Params.kin_memory_s`: while the drive is on and for this long from the
  start, every update is a time step from the last to LIGHT's mean energy.
  The run starts from the gas's Maxwellian with LIGHT's mean energy on the
  grid, under LIGHT's own closure until the first step. A step's
  energy-balance diagnostic takes its own d<eps>/dt.
- The overview's section "The first nanosecond", generated from three more
  runs (the reference point steady and with LIGHT's Maxwellian in the first
  nanosecond, Tahoma steady), and the step's checks.
- Tests for the step's limits and energy and its no-field step where a
  target is out of reach, the remap and the grid's Maxwellian at a mean,
  the first nanosecond as time steps, the afterglow left steady, and the
  table-end test on a step's table as well. `test_docs.py` holds the new
  section's words and the renewal sentence to their numbers.

### Changed
- **The default.** `kin_memory_s` is inf, the whole pulse. In the pulse the
  EEDF follows its history, and in its first nanosecond the field carries
  0.98 to 1.07 of the input power, 0 updates drop the growth term, and
  after 0.1 ns the two energy balances agree to 1.5e-3. Steady against
  stepped, no band at the reference point moves by more than 1.9e-4 of
  itself, and Tahoma's window gain by 3.6e-4 (n_e by 9.2e-4). After the
  pulse the EEDF stays steady. In Morozov's low-pump pure-Kr run under the
  beam the steps leave 52 of the pulse's 139 updates unreached, the
  no-field step 0.11 to 0.23 % hotter than LIGHT's Te between 10 and
  449 ns, where solved steady 56 updates of the first 2 ns missed, by up to
  a factor 115 at the start. Its VUV yield moves by 3.6e-4 of itself
  (0.55437 against 0.55457), and it still takes 235 renewals.

### Corrected
- **The field's excess over the input power was misread.** Since 2026-09-30
  (d) v0_3's README, and since (e) its accuracy overview, said that in the
  pulse's second half the field brings 1.007 to 1.032 times the input power
  because LIGHT's electrons cool through the pulse and the steady balance
  leaves out n_e d<eps>/dt. The excess is the table's renewal. Between
  renewals LIGHT runs on a table made at the last renewal's electron
  density and gas, and n_e grows by up to 5.7 % between renewals there. So
  at each renewal LIGHT's electrons sit off the fresh table's balance,
  their n_e d<eps>/dt -0.7 to -3.2 % of the input power, and they regain it
  with an e-folding time of about 0.2 ns (at 721 ns, -2.8 % falls to
  -0.12 % in 1.5 ns). At 422 ns LIGHT's dUe/dt is -1.35e22 eV cm^-3 s^-1
  right after a renewal and averages +3.6e21 over the next 10 ns. The
  electrons' real cooling, a step's own n_e d<eps>/dt, is -0.03 to -0.21 %
  of the input power, 6 to 44 times less than the excess. Renewed every
  5 ns the excess is at most 1.015, against 1.031, and the bands move by
  6.8e-5, so the excess is harmless. Root cause of the misreading: the
  ratio was explained by the term the steady balance lacks, without
  checking that term's size against the run. Over the 64 renewals of
  0.5-1 us the excess equals LIGHT's cooling at the renewal less the step's
  own to 2.9e-5 of the input power.

### Tried and left
- **Option 3 as the model.** Marched offline with the Joule power set to
  the input power and the mean left free, the EEDF's mean energy stays
  within 1.7 % of LIGHT's in the first nanosecond (1.1 % after 0.1 ns).
  Holding LIGHT's mean (option 2) loses nothing there and keeps LIGHT's
  ledger, so option 3's rebuild of the energy bookkeeping was not needed.
- **A time step in the afterglow.** Through the afterglow of Morozov's
  low-pump pure-Kr run it left 87 of the 95 afterglow updates short of
  LIGHT's mean energy, the no-field step up to 8 % hotter than LIGHT's Te,
  and moved the VUV yield by 4.6e-4 against the steady run. The afterglow
  stays steady.
- **An energy-conserving e-e drift.** At thermal energies on 600 cells the
  e-e part of the flux takes 12 % of a 1e9 eV/s field's power in a
  picosecond's step, and the step's energy balance holds to 1.0e-1 (7.6e-3
  on 2400 cells, 2.8e-5 with e-e left out). Rescaling the e-e drift so its
  own exchange is zero changed nothing (lambda 0.999994). The error is the
  split of a Scharfetter-Gummel flux whose c(z) couples its parts where a
  cell is as wide as kT. Documented, not fixed.

### Verified
- The v0_3 suite from the main checkout, 100 passed, 2 skipped and 1
  xfailed in 7.5 min. `test_docs.py` with `LIGHT_DOCS_SLOW=1`, 6 passed,
  its slow numbers against the overview's runs of this code (the code's md5
  is the same in the worktree and in main). The overview regenerated, 37
  runs in 1159 s on 14 processes, and rebuilt, 7 pages with no undefined
  reference or overfull box.
- The reference point's first 2 ns. Solved steady, 98 updates, 86 with the
  growth term dropped, the field 0.019 to 0.186 of the input power and the
  balance off by up to 0.222, in 111 s. As time steps, 98 updates, none
  dropped and none unreached, the field 0.982 to 1.070 and the balance off
  by up to 0.070 (at the thermal start), in 28 s.
- The renewal reading, from LIGHT's derivatives at every renewal of
  0.5-1 us and through one interval, and from the overview's 5-ns
  convergence run.
- Twelve bugs planted in the time step one at a time, and the fast tests
  fail on all twelve. The no-field step where a target is out of reach
  (returning the field's floor exit as it stands) failed only once its
  check moved from a 1-ps to a 0.1-ns step.
- The drift gate (`test_network.py`) fails in a checkout whose SPECK
  predates REBUILD 0.1.40, which gave K34 the Tg form LIGHT's G5 copies,
  and passes against main's SPECK.

## 2026-09-30 (f): the bottom lines named the wrong largest misses

The boxed bottom line atop each accuracy overview (2026-09-30 (c)) said the
model's largest misses were Hakuta's C-A decays (19.3 and 34.1 % rms) and
the pure-Kr VUV yield at low pump. Checking the claim against every score
it ranks (while scoring v0_3, (e) below) found it wrong in all five
overviews. Burnham and Searles' tau_B, missed by -22.3 % and in the same
scorecards, is larger than Hakuta's Fig. 2 and was left out. v0's box also
gives its default network 32 % of Tahoma's gain two sentences earlier, a
far larger miss.

### Fixed
- `accuracy_overview.tex`, `_v0_1`, `_v0_2`, `_1Dv0_2` and `_v0_3`: the
  sentence now names the misses on data that conflict with the fitted sets,
  Hakuta's and Burnham and Searles' (cited as J. Chem. Phys. 67, 5967
  (1977), "Radiative lifetime of KrF", which NIST's Chemistry WebBook lists
  and whose manuscript is in
  `docs/krf_data/sources/NRL_MR3753_1978_DTIC_ADA055912.pdf`), and the
  pure-Kr yield. v0's sentence sets its default network's Tahoma gain apart
  first. Burnham and Searles is [7], and Morozov moves to [8].
- Root cause: the sentence was written once, for v0's overview, and copied
  into the others with each one's numbers but without ranking its scores.

### Verified
- The four PDFs rebuilt (latexmk, no undefined reference), and the new
  sentence read back from each with pdftotext. Their doc tests pass (v0,
  v0_1 and v0_2: 23 passed, 5 slow skipped. 1Dv0_2: 2 passed, 2 slow
  skipped). For v0_3, `test_docs.py` now checks every other score of the
  default model against the named misses.

## 2026-09-30 (e): v0_3 scored against the measured data, and what that took

The owner said "Proceed" to the proposed next step, scoring v0_3 against
the measured data in a v0_3 accuracy overview. Running ML's scores on v0_3
needed the beam drive's electrons in the kinetic equation and turned up
faults in v0_3's coupling, fixed below with their root causes.

### Added
- `docs/accuracy_overview_v0_3.{tex,pdf}`, `docs/make_overview_v0_3.py` and
  its generated `overview_v0_3_fast.tex` and `overview_v0_3_slow.tex`, and
  `docs/energy_ladder_v0_3.tex` (`make_energy_ladder.py` now knows v0_3).
  No number in the overview is typed by hand. ML's ten kinetic datasets
  have no electrons, and the generator refuses to write a v0_3 score that
  differs from v0_2's. Tahoma, its bands, pure Kr, ML's kinetic loop and
  the reference point run in a pool of processes (about 25 minutes on 14).
  With `LIGHT_OVERVIEW_CACHE` set, the generator keeps each run under the
  md5 of v0_3, ML, the beam drive and itself.
- The scores (ML's Tahoma configuration, v0_2's default network with n_e
  free from 1e10 cm^-3): g0 5.64 %/cm, x0.88 the measured 6.39 (v0_2
  x0.89), Isat x1.11, g0 Isat/Pin x0.99, the pure-Kr VUV yield at low pump
  x1.32 (v0_2 x1.36). Over Tahoma's window v0_3 is ML's self-consistent
  LIGHT-Boltzmann loop to 0.06 % in Te, 0.56 % in n_e, 0.04 % in g0 and
  0.04 % in formation, where v0_2's closure is off by 0.7, 2.5, 1.0 and
  0.5 %.
- Electrons born with a spectrum: `boltzmann.State.inj`, `Grid.injection`
  (linear on [0, top] with a given mean, as ML's ee.injection, to 1.6e-5 in
  the mean), `Model.kinetic_injection` (none for a given power) and
  `BeamModel.kinetic_injection` in `ebeam_screen.py` (the secondaries, one
  per ion pair, on [0, 9.915] eV, mean 4.506 eV). Earlier versions never
  call it.
- `v0_3/tests/test_docs.py` (the overview's numbers defined, current and
  matching its words, ML's Tahoma configuration) and eight tests in
  `test_kinetic.py`: the injection on the grid, a beam seed growing on its
  secondaries, the beam drive's run, a trigger at the run's start, rows
  dropped after construction, the table past its ends, a recombining
  afterglow keeping its growth term, and the powers as the solver's own
  flux split. The threshold test now also runs on a beam's secondaries, and
  the undriven test checks that the equilibrium exchanges no energy.

### Changed
- The drift triggers watch every input of the kinetic state (Te, n_e, the
  F2 fraction, Tg, the Penning and recombination frequencies per electron
  and the Kr* fraction), each above a floor, and after the pulse the
  renewal interval grows to `kin_afterglow_frac` (0.1) of the time since
  the pulse ended.
- Past the table's ends LIGHT carries it with its own Maxwellian dependence
  on Te, the rates by their logarithmic slope, the captured energies by
  their ratio and the elastic power by its difference.
- The growth rules. The term also stays where the gas relaxes the shape ten
  times faster than the population turns over (`relax_freq`). A cold start
  tries it first, on a budget, with nu from the net of the sources, and a
  positivity reset raises a gain. Where even its no-field solution is
  hotter than the target, that solution stands, unreached, when the
  population grows or the relaxation dominates. `_without_growth` returns
  None where no solution exists, and the EEDF fails only when neither
  exists.
- The energy-balance diagnostic takes the drive's electron power from
  LIGHT's own equations, and a field carries only what is not injected. The
  report's electrons line gives no field ratio under a beam.

### Fixed
- **The table's elastic power was not the solver's.** Symptom: late in
  Morozov's low-pump pure-Kr afterglow LIGHT's Te fell to 0.930 of the gas
  temperature, below the kinetic solution's 0.978. Root cause:
  `Grid.elastic_power` took F + (kTg/e) dF/deps by central differences,
  while the solver's fluxes are Scharfetter-Gummel. On 600 cells to 12 eV a
  cell is 0.77 kTg at 300 K, and central differences gave the solver's own
  equilibrium a loss of 4.7 % of its cooling term. The additive extension
  below the table carried that loss down, putting LIGHT's elastic zero at
  0.93 Tg. Fix: `Solver.elastic_power` and `field_power` are the parts of
  the solver's own flux (the central form with the diffusive part times
  (z/2) coth(z/2), z from all of W and D), and the central form is gone.
  LIGHT's Te there is now the kinetic one, 0.977 of the gas temperature at
  20 us, the grid's own thermal bias.
- **A beam's secondaries had no place in the kinetic equation**, and once
  they had, Tahoma failed at t = 0. Root cause: in a seed of 1e10 cm^-3 the
  secondaries outrun every sink, so the solution without the growth term
  does not exist (the first cell would take up a net gain), and the growth
  term started from nu = 0 lost positivity. Fixed by the cold-start rules
  above.
- **LIGHT's Te stalled below the kinetic mean at Tahoma**, 1.3 against 1.6
  eV. Root cause: past the table's bottom the rates and powers were held
  flat, so LIGHT's energy equation had no loss that falls with its Te.
  Fixed by the Maxwellian extension.
- **A pure-Kr afterglow renewed the table 2000 times.** Root cause: a fixed
  10-ns interval over a 20-us afterglow, the only guard on the inputs the
  triggers did not watch. Fixed by the triggers and the growing interval.
  Morozov's low-pump run now takes 235 renewals.
- **A pure-Kr run under the beam looped at t = 0**, renewing its table
  12,405 times without moving. Root cause: Te passes its trigger 1.6e-16 s
  into the run, inside the 8.9e-16 s to which scipy 1.13 locates an event
  (brentq, absolute tolerance 4 machine epsilons), and in pure Kr the event
  came back at t = 0 itself (at Tahoma at 3.0e-16 s). Fix: a segment whose
  event lands where it began takes a short step without the triggers,
  doubled while it repeats.
- **The late afterglow went below the gas temperature** in two more ways.
  The elastic power past the table's bottom was carried as a ratio of the
  Maxwellian's powers, which divides by a number passing through 0 near Tg.
  It is carried by their difference. Dropping the growth term in a
  recombining afterglow put the recombined electrons back cold. The
  relaxation clause keeps it there.
- **Rows ML drops after construction stayed in the kinetic equation.** Root
  cause: `_kinetic_setup` ran once, when the model was built, and ML drops
  and rescales rows afterwards. It now runs at the start of every run. The
  table's Maxwellian end values looked rows up by list position, which a
  dropped row shifts (found reading, before any run used it). They are
  looked up by row.

### Verified
- The overview's numbers from fresh runs of the final code (34 runs in
  about 25 minutes on 14 processes), and `test_docs.py` green (6 tests, the
  slow check rendering the runs' cache). The PDF built and read page by
  page.
- The scores at printed precision did not move with the elastic fix. It
  moved the reference point's energy-balance gap from 1.9e-5 to 3.5e-5 and
  the worst ledger closure to 1e-13.
- A planted-bug check of the new code, each bug in a scratch copy of its
  own with its own bytecode cache: 34 bugs in `boltzmann.py`, `model.py`,
  `report.py` and the beam drive, and the fast tests fail on 30. Two
  survivors change only the cost (a cold start at nu = 0, which the first
  positivity reset corrects, and a cold start that tries the solution
  without the growth term first). Two guard fallbacks no test reaches
  (keeping the growth term's no-field solution where the solution without
  the term misses the target too, and a solution without the term that does
  not exist). The first round left 13 alive, and the tests for the table's
  ends, the recombining afterglow and the beam's heating were written for
  them.
- From the main checkout v0_3's tests pass, 96 with 2 slow ones skipped
  and 1 xfailed (Hakuta's C-A decays, as in v0_2), in 13 minutes. So do
  the other overviews' doc tests, 23 passed and 5 slow skipped for v0 to
  v0_2 and 3 passed and 1 skipped for 1Dv0_2.

## 2026-09-30 (d): v0_3, v0_2 with the kinetic EEDF in the loop

The owner asked: "Then let us proceed with making v0_3 where we add ML to
v0_2", and of the routes offered chose the kinetic EEDF in the loop. v0_3
is v0_2's copy in which every electron rate, the elastic loss and the
energy a captured electron carries come from the two-term Boltzmann
equation with electron-electron collisions, solved during the run at
LIGHT's own mean energy. The chemistry, the bands and the input stay
v0_2's.

### Added
- `v0_3/` and `run_v0_3.py`. `v0_3/boltzmann.py` is ML's two-term scheme
  (the website LXCat build's solver and ML's e-e term) on sparse matrices,
  with the network's Penning sources, its recombination template,
  attachment and the growth term, and only the processes the network
  books in its energy terms (Kr's elastic loss alone, as G6 books it).
  `model.py` reads it through a table of the kinetic rates, captured
  energies and elastic power at five mean energies around LIGHT's own,
  renewed as Te, n_e and the F2 fraction drift and at least every 10 ns.
  `eedf='two_t'` keeps v0_2's closure, and under it v0_3 is v0_2. The
  report's header says how the electrons were closed and, for the kinetic
  EEDF, how far its energy balance is from LIGHT's.
- `v0_3/README.md` (how it works, the reference point, what was found, the
  checks and the limits) and section 6 of `v0_3/SOURCES.md`.
- `v0_3/tests/test_kinetic.py`, 13 tests: the solver against ML's dense
  one, its particle and energy balances, the undriven Maxwellian, the
  fixed-field bistability and the one fixed-mean solution, the growth
  decisions and their threshold, the grid, a 0.3-us run's closure and its
  energy agreement with LIGHT, v0_2 under v0_2's closure, and (with
  `LIGHT_KINETIC_SLOW=1`) the reference point's convergence. v0_2's tests
  run on v0_2's closure, as `tests/test_model.py` says.

### Found
- **At a fixed field the problem has two solutions** where Kr's Ramsauer
  minimum meets strong e-e collisions. 0.12 ns into the reference pulse, at
  0.4 Td, a hot start gives a mean energy of 0.282 eV and a cold one 0.042
  eV, and the field a mean energy needs falls from 0.444 Td at 0.06 eV to
  0.388 Td at 0.2 eV before it rises. A search on E/N for a mean energy
  (Brent's method, the first design) then found no sign change. The field
  is now found inside every step of the iteration, at a fixed mean energy,
  where the solution is one.
- **ML's e-e quadrature condenses an undriven EEDF.** With F constant in
  each cell the identity that leaves a Maxwellian alone, 3T A1 = 2(A2 +
  eps^3/2 A3), fails at the lowest boundaries by more than the elastic
  coupling to the gas, and with no drive pure Kr at n_e/N = 1e-5 fell to a
  mean energy of 0.029 eV on 600 cells and 0.007 eV on 2400. The integrals
  are now taken on ln F linear between cells, exact for a Maxwellian, which
  then comes out to 1e-14. ML's own `ee.py` has the cell quadrature and was
  left as it is. Whether it moves ML's beam-driven states was not checked.
- **Factors over a Maxwellian fail where the two mean energies part.** The
  first coupling carried the kinetic rates after the pulse as ML does, as
  factors over the Maxwellian at the kinetic mean energy applied at LIGHT's
  Te. In a 0.3-us run the kinetic EEDF's temperature was 0.29 eV while
  LIGHT's was 1.0 eV, the elastic factor reached 22 and F2's dissociation
  factor 6e4, and LIGHT's Te fell from 1.1 to 0.41 eV in 0.12 ns. The table
  at LIGHT's own mean energy now serves the afterglow as well.
- **A grid top scaled by the decades F falls oscillated** between 41 and
  13 eV at the reference point, since the fall is slow in the bulk and
  steep past the Kr thresholds, and returned the 13-eV grid, three decades
  long with Kr's ionization threshold off its end. The top is now read off
  F, on a lattice of grids, and no grid is tried twice.
- **The growth term.** In the first nanosecond it moves more energy,
  replacing the electrons F2 attaches, than the field brings, and it is
  dropped there. In the Penning-fed afterglow its no-field solution has no
  positive solution while the one at LIGHT's mean energy exists, so that
  absence no longer counts as unreachable. A growth attempt that follows a
  solution without the term runs on a budget, 60 steps and 2 resets.
- **Without the growth term the first cell takes up the net loss of
  electrons**, cold. At low field with strong attachment that cools the
  no-field solution below the gas temperature, to a mean energy of 0.015
  eV at the reference run's start. v0_3 holds LIGHT's mean energy, so this
  moves only the field and the shape below a few tenths of an eV.
- **Near no field the secant on E/N never arrived.** Its slope was 1e-6
  there, where the Penning electrons and not the field set the mean, and
  one table took 58 s. The field is now bracketed, with the floor tried
  directly, and found by Brent's method.

### Verified
- The reference point: 4.79 % at 248 nm against v0_2's 4.67 %, 19.7 % at
  147 nm against 20.0 %, all bands 30.3 % against 30.4 %, 650 K against
  649 K at the end, the ledger closed to 1.5e-15, in 217 s. In the pulse's
  second half the field brings 1.007 to 1.032 times the input power, and
  with LIGHT's n_e d<eps>/dt put back the two energy balances agree to
  1.9e-5. The kinetic growth frequency equals LIGHT's d n_e/dt / n_e
  at every pulse table that keeps the growth term.
- Converged: renewing at half the drift or every 5 ns moves no band by
  more than 7e-5 of itself, and 1200 cells by 1.8e-4.
- The solver is ML's dense one with ML's quadrature to 2e-9 in F and 4e-8
  in every rate at three states, and its own quadrature moves F by 1.1e-5
  at most there. Its energy balance closes to 1.6e-6 in the pulse and
  5.7e-5 with no field, and without the growth term to 2.9e-4 with the
  first cell's injection, a residual first order in the cell width (5.9e-4,
  2.9e-4 and 6.6e-5 on cells of 0.048, 0.024 and 0.006 eV).
- The Maxwellian identity by quadrature at four (T, eps) pairs, to 5e-16.
- A planted-bug check: twelve bugs in `boltzmann.py`, one at a time, and
  the fast tests failed on eleven. The survivor removes the guard against
  trying a grid twice, which no test reaches. Each run used its own
  bytecode cache, and the file was restored after every one.
- From the main checkout: v0_3's tests 83 passed, 1 skipped (the slow
  convergence check) and 1 xfailed (Hakuta's C-A decays, as in v0_2), in
  5 minutes, and `run_v0_3.py` prints the reference point's summary with
  its electrons line.

## 2026-09-30 (c): a bottom line at the top of each accuracy overview

The owner asked: "Please have the bottom line model accuracy at the top of
the accuracy overviews. E.g. This model v0_2 replicates 99% of the
chemistry kinetics with a 90% agreement to the experimental gain of
[CITATION], etc." Each of the four overviews opens with a boxed bottom
line in that form, with the model's own numbers. The kinetics are given as
the rms deviation from each measured set, since no measurement gives a
share of the kinetics that a model replicates.

### Added
- The box under each title, in the form of the review's bottom line. Its
  first sentence gives the rms deviation from the KrF decays the model was
  fitted to (Eden 1978, Gadomski 1992) and from decays no fit used (Xu
  1993), and the share of Tahoma's measured small-signal gain the model
  gives (Spectra Technology 1989, OSTI 10136715). The rest gives
  Tahoma's saturation intensity and formation efficiency, the C/B
  emission ratio (Yu and Setser 1990) and the largest misses (Hakuta
  1990, Morozov 2008), against seven numbered references. v0_2 gives 2.3
  to 6.1 % and 11.7 % rms and 89 % of the gain. v0's box names both of its
  Tahoma configurations, as the Known-bad list requires: 32 % from its
  default network, and 96 % with its ion, F2-electron and Kr2* groups on
  and n_e free, ML's configuration, which its Table 1 uses. v0_1 gives
  93 % by default and 95 % with the Kr2* group, and 1Dv0_2's centre cell
  86.5 % on its own power input, its gap average on 64 cells 84 %.
- `\pctof` in each overview's preamble prints a factor as a percentage of
  the measurement, 100 times the generated factor to its own precision,
  so no score in a box is typed. v0's overview now inputs
  `REVIEW/review_numbers.tex` (`\rv`) for its default network's score.
- A test in v0's and v0_1's `test_docs.py` that every number the overview
  prints is defined in the generated file it comes from, as v0_2's and
  1Dv0_2's already had.

### Verified
- From the main checkout: v0's tests 74 passed, 1 skipped and 1 xfailed
  (Hakuta), v0_1's 74 passed, 2 skipped and 1 xfailed, v0_2's docs tests
  but the minutes-long fast and slow checks 7 passed and 1 skipped, and
  1Dv0_2's docs tests 2 passed and 2 skipped under `/usr/bin/python3`.
- A planted-bug check: six bugs, one at a time (a missing review key cited
  in v0's box, the review dropping the key v0 cites, a missing overview key
  in v0's box, and a missing key in v0_1's, v0_2's and 1Dv0_2's boxes),
  and each failed the test aimed at it. Every file was restored byte for
  byte, and no bytecode was written.
- The four PDFs build in the main checkout with no errors, overfull boxes
  or undefined references, and their text matches the branch's builds.
  Every box was read in a render and its numbers checked against the
  generated factors they print (0.89 as 89 %, 0.865 as 86.5 %). Each Table 1
  now sits on page 2, and each ladder page still falls just before its map
  and after it in 1Dv0_2's. The boxes' prose passes the style check for
  semicolons, dashes, colon-led clauses and parenthetical asides.

## 2026-09-30 (b): an energy ladder in each accuracy overview

The owner asked: "Add an energy ladder similar to the one on
~/coding/website/simulations/excimer.html for the accuracy_overview
files". All four overviews carry one, v0's and v0_1's included, since the
request names the overview files. No model changed, and no overview's
text changed beyond the ladder's section, the subtitle that names it
and the page breaks that place it.

### Added
- `docs/make_energy_ladder.py`, which draws the ladders of v0, v0_1 and
  v0_2 from each version's `species.py` and `network.py` and writes
  `energy_ladder_v0.tex`, `energy_ladder_v0_1.tex` and
  `energy_ladder_v0_2.tex` (`\laddef`, read as `\lad`, and the picture
  as `\ladpicture`). It follows excimer.html's LADDER view. Every species
  stands at its formation energy, the ion pairs dashed at the sum of
  theirs, and the electron-impact thresholds are ticks. Every row the
  version can run is an arrow with its id, solid in the default network,
  dashed in an optional group and dotted for a band whose default is 0,
  and a species that only an optional group makes is pale. A photon drops
  by its quantum, and the heat its lower state gives the gas as it falls
  apart is a red squiggle with its value. Quenches are red stubs. KrF(B)
  and KrF(C), 80 or 120 cm-1 apart, are a pair of bars that every B and C
  row meets. G6 and VIB1 to VIB4 move nothing between levels and the
  caption names them. The script refuses a row it has no drawing for, and
  every number on a ladder comes from the version's code.
- A section "The energy ladder" in each overview, before the system map
  in v0's, v0_1's and v0_2's and after it in 1Dv0_2's, which sets v0_2's
  ladder since its chemistry is v0_2's by md5. Its paragraph cites only
  `\lad` numbers: KrF(B) against the ion pair and Kr + F, the heat the
  248-nm and 420-nm photons leave, the heat of the harpoon, of K26 and of
  the trimer's formation, and in v0_2 the 147-nm band. v0's ladder shows
  REBUILD's estimated KrF(B), 5.25 eV below the ion pair, where v0_1 and
  v0_2 use the measured D_e(B), 5.54 eV, and its 248-nm photon leaves
  0.35 eV of X-state repulsion against their 0.06 eV. The subtitles name
  the ladder, and the v0 overview gained the `\Fm` macro.
- Two tests in each of v0's, v0_1's and v0_2's `test_docs.py` (the ladder
  file against a fresh render, every `\lad` the overview prints defined),
  and the `\lad` keys in 1Dv0_2's number check.

### Found
- pdfpages does not flush floats (Known-bad). The first build set the
  ladder after the map, on a page that took the map's landscape paper and
  cut off the caption, and in 1Dv0_2 the appendix's longtables held it
  three pages further.
- v0's harpoon heat is 9.915 - 6.95 = 2.965 eV, which binary floating
  point stores just below the tie, so `f"{x:.2f}"` printed 2.96. The
  ladder rounds half up on the decimal value.

### Verified
- From the main checkout: v0's tests 73 passed, 1 skipped (the overview's
  slow numbers) and 1 xfailed (Hakuta) in 69 s, v0_1's 73 passed, 2 skipped
  and 1 xfailed in 2.8 minutes, and v0_2's 79 passed, 2 skipped and 1
  xfailed in 4.7 minutes. 1Dv0_2's `test_docs.py` 2 passed and 2 skipped
  under `/usr/bin/python3`, and 3 passed and 1 skipped (the slow numbers)
  under `~/.venvs/f2sweep/bin/python`.
- A planted-bug check: eight bugs, one at a time (a ladder number edited, a
  ladder key dropped, the picture edited, a missing key cited in v0's and in
  1Dv0_2's overview, a species energy moved, a row left without a drawing,
  the rounding done on the binary value), and each failed the test aimed at
  it. Every file was restored byte for byte, and no bytecode was written.
  The row without a drawing also trips the script's own refusal ("rows with
  no drawing: ['ATT']").
- The four PDFs build in the main checkout with no errors, overfull boxes
  or undefined references, and their text matches the branch's builds. Each
  ladder page is portrait A4 with its caption, just before the map and just
  after it in 1Dv0_2's. Every ladder page was read in a render, and v0's and
  v0_2's KrF columns at 300 dpi. The new prose passes the style check for
  semicolons, dashes, colon-led clauses and parenthetical asides.

## 2026-09-30: how accurate v0_2 and 1Dv0_2 are, written up

The owner asked: "The subdirectory docs/ has an accuracy_overview please
replicate the file for the 1Dv0_2 as well as for v0_2". v0, v0_1, 1Dv0,
v0_2, 1Dv0_2, their documents and the review stay as they were. Only
documents, their generators and docs tests were added.

### Added
- `docs/accuracy_overview_v0_2.tex/pdf` (13 pages) with
  `make_overview_v0_2.py`, which writes `overview_v0_2_rates.tex`,
  `overview_v0_2_fast.tex` and `overview_v0_2_slow.tex` (`\widef`, read as
  `\wi`). It is v0_1's overview for v0_2, plus a section on the Kr2* route
  (Kr2* + F across its band at the reference point, the part of its gain
  that comes after the F2 burnout, pure krypton against the bands) and one
  on the bands (each at both ends at the reference point and at Tahoma).
  Tahoma is scored with v0_2's default network, which is ML's
  configuration now that the Kr2* route is basic, and with Kannari's two
  Kr2* estimates. The formation-links table tries Kr2* + F at Kannari's
  3.0e-10 and at its band's top, since ML's "no Kr2* + F" is v0_2's own
  default. Appendix A prints each band row with its default and band, and
  lists the bands (A.4).
- `docs/system_map_v0_2.tex/pdf` with `make_system_map_v0_2.py` (`\wmdef`,
  `\wm`). It is v0_1's map with a Kr2* node (Kr* + 2Kr drawn inside it, as
  K15 is inside the ions' node), one box for Kr2*'s fates (147-nm radiation
  and Kr2* + F2), the 147-nm band to a new ledger entry, Kr2* + F2's KrF(B)
  joining the ion-ion band and its Kr2F* rising into the trimer, and a box
  with each band's two ends where v0_1's showed the optional groups. All 44
  rows are drawn and 43 nodes balance to 1e-8 of E_dep. The three
  unmeasured rows run at 0 and are printed, and the script refuses to draw
  a run that sets them. Kr2* is #A63A83. The dataviz validator passed it
  against the map's other hues, all pairs on white, at a normal-vision
  Delta E of 16.2 or more and a CVD Delta E of 8.6 or more. The palette's
  own magenta, #E87BA4, failed against heat's red at 13.2.
- `docs/accuracy_overview_1Dv0_2.tex/pdf` (11 pages) with
  `make_overview_1Dv0_2.py`, which needs Fluid1D (Python >= 3.11 with
  numba) and writes `overview_1Dv0_2_fast.tex` and `overview_1Dv0_2_slow.tex`
  (`\didef`, `\di`). It checks that 1Dv0_2's nine chemistry files have v0_2's
  md5 sums, writes nothing if one differs, and takes the chemistry's scores
  from the v0_2 overview. It scores the 1-D layer itself: the stage's rates
  against v0_2's `Model._k` at 800 random cell states, v0_2's reference
  point on 32 to 256 cells, Kr2* + F at the top of its band, the time step
  at 50, 25 and 12.5 ps, and Tahoma on 1Dv0_2's own input from seeds of
  1e13, 1e15 and ML's 1e10. A pgfplots figure shows the cells beside the
  wall on each mesh and n_e and Te there, and v0_2's map stands for the
  interior. The four meshes' curves use one blue ramp, which passed the
  validator's ordinal checks.
- Docs tests. `v0_2/tests/test_docs.py` gains six beside its two review
  tests: the map, the rate tables and the fast numbers against fresh runs,
  the slow numbers with `LIGHT_DOCS_SLOW=1`, every short source held to its
  row's citation, and every number the map and the overview print defined
  in a generated file. `1Dv0_2/tests/test_docs.py` has four: the md5
  identity and the printed numbers under any Python, and with Fluid1D the
  1-D fast numbers and (slow) the runs against fresh ones.
- Pointers in `README.md`, `v0_2/README.md` and `1Dv0_2/README.md`, and the
  two findings below under 1Dv0_2's Limits.

### Found
- v0_2's overview reproduces the review's v0_2 column to the digit: at
  Tahoma g0 5.69 %/cm (x0.89), Isat x1.11 and g0 Isat/Pin x0.99, Kannari's
  estimates x0.95, the top of Kr2* + F x0.97, the trimer's weak law x0.44,
  the posterior 4.99 to 6.18 %/cm, and every kinetic score. With both of
  Kannari's estimates v0_2 gives v0_1's kr2 group, 9.17 %.
- At Tahoma the bands move the gain far less than at the reference point.
  Kr2* + F at its top takes g0 from 5.69 to 6.23 %/cm, where the reference
  point's 248-nm share goes from 4.67 to 9.50 %. Tahoma's window closes
  with 29 % of the F2 left, and 74 % of the band's gain at the reference
  point comes after the burnout. Among the six narrow bands, Kr* + 2Kr at
  5.5e-32 moves Tahoma's gain most, to 5.31 %/cm (x0.83).
- The v0_2 README's pure-Kr yields at 300 kW/cm^3 (0.13, 0.27 and 0.49) are
  for a 1-us pulse, which it does not say. With Tahoma's 450 ns they are
  0.18, 0.33 and 0.52. ML's `pure.kr_vuv` sets Kr* + 2Kr from its `k3`
  argument, so the kr2_form band passed through `rates=` is silently
  ignored. The overview passes `k3`.
- 1Dv0_2's stage holds Kr2+ + e (R11p, a power law) at its 30-eV table
  value above 30 eV, where v0_2's `Model._k` follows the law, so at the
  40-eV ceiling the stage's rate is 17 % high. Below 30 eV the stage
  matches `Model._k`, the density and Tg rows exactly, the Te tables to
  3e-4 (R11p again, just above its 0.02-eV floor, where the table
  interpolates across the law's kink) and the tail to 7.3e-6. Only the drained
  wall cells, a cell that fails to break down and a low seed's breakdown
  reach 30 eV, and at Tahoma from 1e13, where v0_2's Te peaks at 38 eV,
  the centre cell still gives v0_2's g0 to 2.8e-4. `test_1d.py`'s rate test
  draws 40 states whose hottest is at 20 eV, so it cannot see this. Not
  changed, since it is the owner's model: evaluating the power-law rows in
  the stage, as v0_2 does, would remove it.
- 1Dv0_2's cells beside the wall cell do not converge. The first sends out
  4.668, 4.672, 4.684 and 4.751 % at 248 nm on 32 to 256 cells against the
  interior's 4.667 %, and on 256 cells they alternate about the interior
  out to 195 um, with n_e and Te. The cells are 39 um wide and the
  interior's Debye length is 1.0 um, so no mesh here resolves the sheath.
  The gap average still converges first order (its deficit falls x2.00 to
  x2.02 per doubling), and beyond 0.05 cm from the walls every cell is the
  centre cell to 8e-7.
- The interior is v0_2's 0-D run to 1.7e-5 at 248 nm, 8.7e-6 at 147 nm and
  2.3e-6 in the final Tg on 32 to 256 cells, and to 5.8e-7 with Kr2* + F
  at the top of its band. The time step is first order: the Te transient
  is 1.3, 0.6 and 0.4 % at 50, 25 and 12.5 ps.
- `docs/REVIEW/make_review.py`'s comments still say v0_2 has no overview.
  The review is unchanged and still makes v0_2's column itself, which the
  new overview reproduces.

### Verified
- From the main checkout: v0_2's tests 77 passed, 2 skipped (the review's
  and the overview's slow checks) and 1 xfailed (Hakuta), in 5 minutes,
  and with `LIGHT_DOCS_SLOW=1` the overview's slow numbers regenerate byte
  for byte (5.4 minutes). 1Dv0_2's tests 72 passed, 3 skipped and 1
  xfailed under `/usr/bin/python3`, and `test_1d.py` with `test_docs.py`
  25 passed and 2 skipped under `~/.venvs/f2sweep/bin/python`, and with
  `LIGHT_DOCS_SLOW=1` the 1-D overview's slow numbers regenerate byte for
  byte (8.9 minutes on four numba threads).
- The map's numbers and both overviews' fast numbers regenerate byte for
  byte in fresh processes.
- A planted-bug check: nine bugs, one at a time (a short source quoting a
  number its row's citation lacks, undefined numbers in the map and in
  both overviews, a dropped generated table, a dropped profile, a stale
  rate table, a stale map number and a changed chemistry file), and each
  failed the test aimed at it. Every file was restored byte for byte, and
  each mutant changed its file's size, so no cached bytecode stood in.
- The three PDFs build with no warnings but the map's three overfull
  boxes, which v0_1's map has at the same sizes, and every page was read in
  a render. The overviews' prose passes the style check for semicolons,
  dashes, colon-led clauses and (i) lists.

## 2026-09-29 (j): 1Dv0_2, the 1-D model on v0_2's chemistry

After (i) the owner asked: "Then make a 1-D model with v0_2 chemistry".
1Dv0_2 is 1Dv0's 1-D layer on a copy of v0_2, the review's pick. 1Dv0 and
v0_2 stay as they are.

### Added
- `1Dv0_2/`, a copy of v0_2 (its chemistry byte for byte: model.py
  3ac682b6, network.py 7f818d4e, report.py 0d94799b, species.py c912fa5d,
  xsec.py and two_level.py and data/ as v0_1's) with 1Dv0's `fluid1d.py`,
  `sim1d.py`, `run.py`, `tests/test_1d.py` and `tests/conftest.py`, and
  `run_1Dv0_2.py`.
- v0_2's bands in the 1-D layer. The eight band rows (rate kind
  `("param", name[, frac])`) are constants of the run, as v0_2's model
  resolves them. `Params1D.rates` and `--rate NAME=VALUE` set them. The
  summary prints v0_2's bands line and its gas-temperature warnings, read
  in the centre cell, and warns when a seed too low left the centre cell
  at the mean-energy ceiling with no KrF.

### Changed
- 1Dv0_2's copied 0-D tests import its own modules. v0_2's `test_docs.py`
  and its two tests of `run_v0_2.py`'s command line are gone. The summary
  and bands-line checks of the first stay, on the 0-D report.
  `README.md` is rewritten, `SOURCES.md` says it is v0_2's, and
  `__version__` is 0.1.2+1d0.
- The Known-bad bullets written for 1Dv0 now name both 1-D models.

### Found
- At v0_2's reference point the interior is v0_2's 0-D run: 4.66665 % at
  248 nm against 4.66672 % (1.5e-5), 20.0185 % at 147 nm against 20.0184
  %, and 649.35 K at the end, on 64 and 128 cells. With Kr2* + F at the top
  of its band (1e-9) it gives 9.49826 % against 9.49835 %. The gap averages
  (4.521 % on 64 cells, 4.594 % on 128) carry the wall cells, as in 1Dv0.
- At Tahoma on v0_2's own input with ML's constants the centre cell's g0
  is v0_2's to 3e-4 (5.528 %/cm against 5.530 from a seed of 1e13, x0.865
  of the measured 6.39). So the review's v0_2 column is 1Dv0_2's.
- The seed rule holds for 1Dv0_2 too. From ML's 1e10 seed v0_2's 0-D Te
  runs to 73 eV (66 eV from 1e12), above the 40-eV ceiling, and 1Dv0_2
  never breaks down; from 1e13 it does.

### Verified
- `test_1d.py`: 22 passed (the bridge tests with the bands at their
  defaults and with the three unmeasured ones set), and the slow one. A
  mutation check planted 19 bugs, 1Dv0's 16 and three in the bands, and
  the fast tests caught all 19.
- From the main checkout: `test_1d.py` 23 passed with `LIGHT_1D_SLOW=1`
  (66 s) under `~/.venvs/f2sweep/bin/python`, and under `/usr/bin/python3`
  1Dv0_2's tests 70 passed, 1 skipped (`test_1d.py`) and 1 xfailed:
  v0_2's 0-D tests (55) less the `--bands` one, and less its document
  tests. `run_1Dv0_2.py --compare --rate kr2_f=3e-10` runs from main.

## 2026-09-29 (i): 1Dv0 at Tahoma, and the seed it needs

The owner asked how accurate 1Dv0 is against v0_2.

### Found
- 1Dv0 scores as v0_1 does. At Tahoma on its own power input, with ML's
  constants imposed, the centre cell's g0 over 200 to 400 ns is v0_1's to
  2e-4 with the default network (5.93 to 5.96 %/cm, x0.93 of the measured
  6.39) and to 7e-4 with the Kr2* group (6.63 %/cm, x1.04), from seeds of
  1e13 to 1e15 cm^-3 on 64 cells. The gap average is x0.90, lowered by the
  wall cells. So the review's v0_1 column is 1Dv0's, and against v0_2 it
  ties on every kinetic dataset, pure-gas check and the electron closure,
  and comes nearer on Tahoma's gain (x0.93 against v0_2's x0.89 under the
  beam, x0.93 against x0.87 on the power input), while v0_2 comes nearer on
  the formation efficiency (x0.99 against v0_1's x1.06). 1Dv0 cannot run
  the beam drive (LIGHT's rule), so its own-input numbers are the ones it
  has.
- From ML's seed, 1e10 cm^-3, 1Dv0 never breaks down at Tahoma: the
  ceiling cuts 99.999 % of the input and no KrF forms. v0_1 breaks down
  from that seed only because its 0-D Te runs to 518 eV in 0.1 ns (277 eV
  from 1e12, 38 eV from 1e13, 7.9 eV from 1e14). At the ionization degree
  of a 1e10 seed v0_1's tail ratio is about 0.004, so T_t = r Te reaches
  the Kr thresholds only at such a Te, and 1Dv0's tables stop Te at 40 eV.
  v0_1's g0 does not depend on the seed (5.94 %/cm from 1e10 to 1e14,
  5.96 at 1e15).

### Changed
- `1Dv0/README.md`: the Tahoma result, and the seed under Limits.

### Verified
- v0_1's 0-D runs in this check reproduce the review's numbers: 6.56 %/cm
  on its own input with the Kr2* group at ML's seed, and 5.94 with the
  default network.

## 2026-09-29 (h): 1Dv0, the 1-D model: v0_1's chemistry across a gap on Fluid1D

The owner copied v0_1 to `1Dv0/` and asked: "Make v0_1 0-D and under 1Dv0
use Fluid1D to make a 1-D model". v0_1 stays as it is and is the 0-D
model. 1Dv0 is the 1-D model, on Nathaniel's Fluid1D
(`../nathaniel/1DModelKrF/src`, imported and never edited). The input
stays v0_1's given power, with no drive, beam or cavity module.

### Added
- `1Dv0/fluid1d.py`, the bridge. It loads Fluid1D, builds its species
  records from `species.py`, converts the two pinned LXCat files to
  Fluid1D's JSON and CSV, and writes a reaction file of the 12
  electron-impact rows, from which Fluid1D builds its electron transport
  table. It holds the network as per-cell arrays, with counter columns for
  heat, emission by band, energy returned to the electrons, captured
  electrons and KrF formed, and v0_1's rates on v0_1's own Te nodes (the
  tail on its (Te, T_tail) nodes). The stage runs Fluid1D's linearly
  implicit trial step in an outer substep loop that reads v0_1's rates at
  each substep's starting energy.
- `1Dv0/sim1d.py`: `Params1D` (v0_1's inputs plus the gap, cells, dt and a
  power profile), `LightSimulation` (Fluid1D's `Simulation` with grounded
  walls, real transport, the stage, v0_1's elastic loss, Tg per cell and an
  energy ledger), `run()` and `Result1D`.
- `1Dv0/run.py`, the 1-D command line (it replaces the copied 0-D one),
  and `run_1Dv0.py`. `1Dv0/tests/test_1d.py` (13 tests and a slow one),
  `tests/conftest.py`, and `1Dv0/.gitignore` for `runs/`.

### Changed
- 1Dv0's copied 0-D tests import 1Dv0's own modules (they imported
  v0_1's). The copied `test_docs.py`, which checks v0_1's documents, and
  the 0-D command-line test are gone. `1Dv0/README.md` is rewritten for the
  1-D model, `SOURCES.md` says it is v0_1's, and `__version__` is
  0.1.1+1d0.
- The chemistry is v0_1's byte for byte (`model.py`, `network.py`,
  `species.py`, `xsec.py`, `two_level.py`, `report.py`, `data/`).
- `README.md` points to 1Dv0 (intro, quick start, files).

### Found
- At v0_1's reference point in a 1-cm gap the interior is v0_1's 0-D run.
  On 32 to 256 cells the centre cell's 248-nm share is 4.73226 to 4.73234 %
  against v0_1's 4.73236 %, and it ends at 762.48 to 762.49 K against
  762.49 K. Under 0.1 % of the energy reaches the walls, so a uniform power
  leaves the 1-D model nothing to add inside the gap.
- The cells next to the walls are unresolved sheaths. The walls drain
  their electrons, and the mean-energy ceiling cuts 6.21, 3.09, 1.81 and
  0.95 % of the input on 32, 64, 128 and 256 cells. Fluid1D's split gives
  their electrons the field's positive work only and keeps out a phantom
  negative part of 88, 17, 5.1 and 1.6 %. The gap average (4.438 to 4.696
  %) approaches the interior's value first order in dx.
- 1Dv0 has to carry the two-temperature tail. Without it v0_1's reference
  point gives 5.82 % at 248 nm and a mean Te of 1.50 eV, against 4.73 %
  and 4.05 eV. Fluid1D's rate tables are functions of the mean energy
  alone, so the stage reads its own, with the tail ratio per cell.
- Fluid1D's substep order, the source added before the rates are read,
  set a steady Te low by dt times the source over n_e (3 to 4 % at 50 ps).
  Read at the substep's start the steady state is exact, and the transient
  is first order in dt (1.5 % at 50 ps, 0.8 % at 25 ps and 0.4 % at 12.5
  ps in a 20-ns pulse).
- Fed through Fluid1D's lead/trail predictor, the given power left a large
  negative correction after the pulse, and stages hit the 8000-substep
  cap. It now reaches each half-step over that half-step's own window.
- The 2-us reference point takes 47 s on 32 cells and 127 s on 256 at
  dt = 50 ps, on four numba threads.

### Verified
- Every rate against `Model._k` at random states: the Tg, density and ion
  rows exactly, the Te tables to rounding, the tail to 7.8e-6 (2e5
  points). The cross sections through Fluid1D's own loader, exactly. The
  stage conserves energy per cell, Kr and F atoms and charge. The gap's
  energy, recomputed from its state, closes against the ledger to below
  1e-12 of the input.
- From the main checkout under `~/.venvs/f2sweep/bin/python` (Python
  3.14, numba 0.65.1, and pytest 9.1.1, installed there for this):
  `test_1d.py` 14 passed with `LIGHT_1D_SLOW=1`, in 59 s. Under
  `/usr/bin/python3`, 1Dv0's tests 67 passed, 1 skipped (`test_1d.py`) and
  1 xfailed, which is v0_1's 74 less its 6 document tests.
  `run_1Dv0.py --compare` runs from main.
- A mutation check planted 16 bugs in `fluid1d.py` and `sim1d.py`, and the
  fast tests caught all 16. Two of them needed new checks (the starting
  state and the KrF-formed counter).

## 2026-09-29 (g): the review's bottom line, v0_2 is the best model

The owner asked for a bottom line in the comparison review saying which
model is best and why.

### Changed
- `docs/REVIEW/v0_vs_v0_1_vs_v0_2.tex` opens with a boxed bottom line.
  v0_2 is the best of the three. It is more accurate than v0, it scores
  level with v0_1, and it is the only version whose default network rests
  on measurements and theory alone. v0_1 comes closer on Tahoma's gain,
  x0.93 by default and x0.95 as ML scores it against v0_2's x0.89, but only
  by leaving out the measured Kr* + 2Kr or by adding Kannari's two
  unmeasured rates, and with those two v0_2 gives v0_1's 6.10 %/cm. v0_2's
  formation efficiency, x0.99, is the nearest of the three, and its bands
  say how far the answer can move. The box replaces the opening paragraph,
  and the provenance note now rides with the last table, so it ends the
  document wherever the floats land.
- The LIGHT and v0_2 READMEs carry the bottom line.

### Verified
- Every number in the box is a generated macro, and v0_2's fast docs test
  passes. review_numbers.tex did not change. The PDF builds without
  warnings at 4 pages, and the changed pages were read in a render.

## 2026-09-29 (f): the review compares v0, v0_1 and v0_2

The owner asked for the REVIEW to compare all three versions.

### Changed
- `docs/REVIEW/v0_vs_v0_1.tex` is now `v0_vs_v0_1_vs_v0_2.tex` (4 pages),
  with v0_2 in every table. The verdict is that v0_1 and v0_2 are equally
  accurate and both more accurate than v0. Nothing measured separates v0_1
  from v0_2, and at Tahoma, the one test where their numbers differ, the
  gain lies nearer v0_1's and the formation efficiency nearer v0_2's.
- `make_review.py`. v0_2 at Tahoma with its default network, with
  Kannari's two Kr2* estimates, at the top of the Kr2* + F band, with the
  trimer's weak law and on its own input. v0_2's column of the accuracy
  tests, made with the v0_1 overview's ML calls (`--scores`), since v0_2
  has no overview. v0_2 at the reference point, its eight bands each alone
  at their two ends, and v0_2 in the notes' grid. The ML processes run side
  by side, about 6 minutes in all. v0's and v0_1's numbers regenerate
  unchanged.
- `v0_2/tests/test_docs.py` checks that the review prints only numbers its
  generated files define, and runs `make_review.py --check` with
  `LIGHT_DOCS_SLOW=1`. The READMEs point to the new file. v0_1's README,
  which the owner froze, still names `v0_vs_v0_1.pdf`.

### Fixed
- Entry (e) and v0_2's README gave Kr2* + F's largest factor on the grid
  as 2.6, from rounded table values. It is 2.652, so 2.7.

### Found
- At Tahoma v0_2's default gives g0 5.69 %/cm (x0.89), Isat x1.11 and a
  formation efficiency of x0.99. v0_1 as ML scores it gives x0.95, x1.10
  and x1.06, and v0_2 with Kannari's two estimates gives v0_1's 6.10. At
  the top of the Kr2* + F band v0_2 gives x0.97, and with the trimer's weak
  law x0.44.
- Over the posterior of the B/C constants v0_2's gain spans 4.99 to 6.18
  %/cm, below the measured 6.39, and v0_1's 5.38 to 6.63.
- The kinetic datasets, the pure-gas checks and the electron closure score
  identically in v0_1 and v0_2.
- Across the notes' grid v0 and v0_1 differ by up to 91 % and v0_1 and
  v0_2 by up to 11 %. v0_2's two wide bands reach farther than that spread
  at 10 of the 15 points.

### Verified
- The review's scores path, run on v0_1, reproduces all 30 of v0_1's
  overview numbers.
- `make_review.py --check` regenerates review_numbers.tex byte for byte,
  as v0_2's slow docs test and as v0_1's (10.7 minutes for both). The fast
  docs test fails on a planted wrong key.
- The PDF builds without warnings at 4 pages, and every page was read in
  a render.

## 2026-09-29 (e): v0_2, the Kr2* route in the basic network and the rates carried as bands

The owner supplied Jancaitis's 1983 LLNL thesis (Recent_KrF_Papers/,
checktext CLEAN) and asked for v0_1 to stay as it is and for the changes to
go to v0_2, their copy of v0_1. This carries out steps 2 to 6 of the side
discussion, the papers' findings and the citations v0_1 corrected after
the copy.

### Added
- `v0_2/` as its own version, 0.1.2, launched by `run_v0_2.py`. The
  owner's copy is commit a75dd95, and every change below is a diff on it.
- The Kr2* route in the basic network (step 6). Kr* + 2Kr, Kr2* radiation,
  Kr2* + F2 in two channels, Kr2* + F and Kr2* + Kr2* are basic rows, and
  `--add kr2` is refused there. v0 keeps REBUILD's copies in its group.
- `network.BANDS` and the rate kind `("param", name[, frac])` (steps 2 and
  3). A measured rate keeps its measurement as the default and spans the
  other published values. A rate with no measurement is 0 by default and
  runs to the largest published value. The eight bands are Kr2* + F,
  Kr2* + Kr2*, Kr2F* + F, Kr* + 2Kr, Kr2* + F2, Kr* + Kr*, Kr2F* + F2 and
  the trimer law. `Params.rates` sets them, and `run_v0_2.py --bands` runs
  each end one band at a time.
- Kr2F* + F -> 2Kr + 2F (Q2FF), 0 by default, band 0 to 1.8e-10.
- The report warns above 451 K, the top of Kr* + 2Kr's measured
  temperature law (step 5). Its header names the bands set and the
  unmeasured ones at 0, and `row_dict` has the 147-nm share.

### Changed
- `v0_2/species.py`. Kr2+'s 1.15 eV is cited to the photoionization
  measurement (Ng et al. 1977, via Krauss & Mies Table 2.1). Kr2* keeps
  REBUILD's estimated well, since no binding energy is among the sources
  held (step 5).
- The K2F_F and ion-ion citations v0_1 corrected in (c) and (d).
- `v0_2/README.md` and `SOURCES.md` (section 5 new, and section 2's
  ion-ion rows corrected).
- ML runs v0_2 (ML/CHANGELOG.md).

### Found
- At the reference point v0_2 gives 4.67 % at 248 nm and 20.0 % at 147
  nm, and the gas ends at 649 K (v0_1, 4.73 % and 762 K). After F2 burns
  out the Kr* that v0_1 turned into Penning ionization and heat radiates
  at 147 nm.
- Two bands move the 248-nm share, Kr2* + F (4.67 to 9.50 %) and the
  trimer law (4.67 to 2.22 %). Each of the other six moves it by 0.21
  points or less. Across the notes' ranges Kr2* + F spans up to a factor
  2.7, at 1-2 atm and 300-500 kW/cm^3, and hardly matters at 0.1 atm.
- Jancaitis's model has no Kr2* + F row (Table II pp.94-95). His F-atom
  quench of the trimers, 0.5 to 1 times F2 (p.153), was inferred for
  Ar2F* (Fig. 11), so for Kr2F* it is an analog and only sets the top of
  its band. He measured Kr2F* + F2 = 1.8 +- 0.2e-10 and 225 +- 21 ns
  (Table II p.159), and his model's Penning rates, 1.5e-10 for Kr* + Kr*
  and 1.0e-10 for Kr2* + Kr2*, come from Werner (thesis 1975, and IEEE
  JQE 13, 769).
- Step 4, pure krypton. At Morozov's pump every Kr2* band moves ML's VUV
  yield by 0.001 or less (0.571 against the measured 0.42 +- 0.05), so
  the route cannot explain the miss and that point cannot score it. At
  300 kW/cm^3 the bands move the yield from 0.13 to 0.49. A yield against
  pump rate would score them, and none is held.
- v0_1's harpoon-seed test assumed every Kr* is harpooned. With Kr* + 2Kr
  in the network 9 % forms Kr2* at 1 atm, and the test now checks that
  share analytically.

### Verified
- v0_2 with Kr2* + F at 3.0e-10 and Kr2* + Kr2* at 5.0e-10 reproduces v0_1
  with its kr2 group to 6e-9 (9.17282 %), and with Kr2* + Kr2* alone that
  group without K2F_F to 2e-9 (4.65249 %).
- v0_2's tests, 71 passed and 1 xfailed. A mutation check planted 13 bugs
  in the new code, and all 13 failed the suite.
- LIGHT and ML together, 258 passed, 3 skipped and 4 xfailed.

## 2026-09-29 (d): the owner's four 1977-78 papers read, still no Kr2* + F, and the ion-ion link corrected

The owner supplied four of the five papers asked for in (c), in
Kinetic_Modeling/Recent_KrF_Papers/, all checktext CLEAN. Jancaitis's 1983
LLNL report is still missing.
- Rokni, Jacob & Mangano, PRA 16, 2216 (1977)
- Jacob, Rokni, Mangano & Brochu, APL 32, 109 (1978)
- Rokni, Mangano, Jacob & Hsia, IEEE JQE 14, 464 (1978)
- Quigley & Hughes, APL 32, 649 (1978)

### Found
- None of the four gives Kr2* + F -> KrF* + Kr or a bound on it. Quigley
  & Hughes ran photolytically at 2e-4 to 7e-4 F2, far too few F atoms, and
  the Avco papers never mention the reaction. The band stays 4.65 to
  9.49 %.
- The ion-ion link is not unmeasured, as the v0_1 documents said. Jacob et
  al. fitted KrF* pulse shapes in e-beam Kr/F2 to Kr+ + F- (+M) = 2e-7 at
  60 Torr and 3e-7 cm3/s at 120 Torr. v0_1's Flannery form, screened at
  their 1e14 cm-3 ions, gives x1.2 and x1.6 of these. Rokni et al. changed
  the beam current x25 and found that Kr2+ + F- forms KrF* "with a
  branching ratio near unity", as efficiently as Kr+ + F-. That test used
  relative intensities, so the absolute KrF* yield, like the B/C split and
  the rate near 1 atm, is still unmeasured. v0_1's values do not change.
- Quigley & Hughes state that the KrF*/Kr2F* split of Kr2* + F2 "is not
  known and cannot be discerned in this experiment". Their Table I carries
  Kr* + 2Kr = 5.5e-32 from Hughes 1977, the high end of the spread.

### Changed
- `v0_1/network.py`. The ion-ion citation (K26, K108) records Jacob et al.
  and Rokni et al. No number moved.
- `docs/accuracy_overview_v0_1`. The reliability item on ion-ion
  recombination and K26's short source say what is and is not measured.
  A.5 keeps its heading with its table through needspace, 10 pages.
- `v0_1/README.md` item 2, and an addendum in
  `docs/krf_data/FORMATION_LINKS_SEARCH.md`.
- ML's registry (ML/CHANGELOG.md).

### Verified
- v0_1 and ML, 115 passed, 2 slow skipped, 2 xfailed. The v0_1 fast
  numbers regenerate unchanged, and only the rate table's K26 source
  moved. The PDF rebuilt clean and the changed passages were read in it.

## 2026-09-29 (c): Kr2* + F searched in the owner's papers, no measurement, and the band restated

The owner pointed to Recent_KrF_Papers/, PIC_papers/ and
docs/KrF_Reaction_Rate_Tables.pdf for the step-2 search on Kr2* + F. The
22 files in the two directories are byte-identical copies of files
already in docs/.

### Found
- None of the 76 local PDFs measures or bounds Kr2* + F -> KrF* + Kr.
  Kannari 1985, the owner's rate tables (rows 121 and 124, ref. [16]), Ong
  1996 and Ma 2024 all carry 3.0e-10, and all of them trace to SRI
  MP-76-99.
- That report tabulates no such rate. Its model leaves Rg2* + F out
  (p.45), and its one number is conditional and for argon, "if the
  reaction Ar2* + F -> ArF* + Ar is sufficiently rapid (3-10 x 10^-10)"
  (p.65). Chen & Payne 1979 (Fig. 31) assign the slow KrF* component SRI
  took for Kr2* + F to ion-ion recombination.
- 3.0e-10 is therefore no upper bound. The reference point's 248-nm share
  saturates in the rate, 4.65 % at zero, 5.85 % at 1e-11, 8.50 % at
  1e-10, 9.17 % at 3.0e-10 and 9.49 % at 1e-9. Half the gain is there at
  about 2.5e-11, so a measurement would matter between 1e-12 and 1e-10.

### Changed
- The v0_1 overview and the review state the band as 4.65 to 9.49 % and
  say that Kannari's value is no bound (`make_overview_v0_1.py` computes
  the shares at 1e-11 and 1e-9). The overview's A.5 starts on a page of its
  own, 10 pages.
- `v0_1/network.py`. K2F_F's citation records p.45, p.65 and Chen & Payne.
  No number moved.
- `v0_1/README.md` item 5, `docs/krf_data/FORMATION_LINKS_SEARCH.md`
  section 3 (a dated addendum).

### Verified
- v0_1's tests, 74 passed, 2 slow skipped, 1 xfailed. Both PDFs rebuilt
  clean, and the changed passages were read in the PDFs.

## 2026-09-29 (b): the Kr2* group's doubling taken apart, and its "measured constants" corrected

The owner brought in a side discussion. The review and the v0_1 overview
called the Kr2* route's constants measured, which is too strong, and the
doubling at the reference point might rest on the route's one unmeasured
formation row after F2 burns out. The first step asked for was to take the doubling
apart, at the reference point and at Tahoma.

### Found
- The doubling rests on Kr2* + F -> KrF(B) + Kr (K2F_F) alone. At the
  reference point the Kr2* group raises v0_1's 248-nm share from 4.73 to
  9.17 %. Without K2F_F it gives 4.65 %, and dropping any other row of the
  group leaves 9.20 to 9.44 %. 80 % of the gain, 3.56 of 4.44 points,
  comes after F2 burns out (665 ns with the group on). The 2.1e17 cm^-3 F
  atoms left then make Kr2* + F 17 times faster than Kr2* radiation, and
  each KrF photon returns its F atom, so the F atoms keep turning Kr2* into
  248-nm light for the rest of the pulse.
- At Tahoma, where 30 % of the F2 is left at 400 ns, the group moves v0_1's
  g0 from 5.94 to 6.10 %/cm. Without K2F_F it gives 5.69, x0.89, and v0
  gives 5.71. Both overviews' Tahoma scores include the row.
- Kr2* + Kr2* Penning ionization (PEN_KR2, 5.0e-10) cites Kannari's ref.
  87, Lorents & Olson, an unpublished 1973 SRI semiannual report that gives
  Ar2* the same value. No measurement of it is among the sources held.
- Four of the six Kr2* rows rest on measurements, Kr* + 2Kr, the 270-ns
  lifetime and the two halves of Kr2* + F2, and the even split of Kr2* + F2
  between its products is Kannari's assignment.

### Changed
- `docs/accuracy_overview_v0_1`. The introduction names the trimer law and
  Kr2* + F as the open questions and gives the Tahoma score without Kr2* +
  F. The scorecard has that row, the reliability item reports the group as
  a band, and Appendix A.2 no longer calls the route measured. The closing
  footer repeated page 1 and Section 1 and is gone, so it stays 9 pages.
- `docs/REVIEW/v0_vs_v0_1`. The Kr2* sentence reports the band, and the
  table has Tahoma without Kr2* + F for both versions.
- `docs/system_map_v0_1`. The groups box gives the Kr2* group without
  Kr2* + F.
- The generators. `make_overview_v0_1.py` has `_kr2_split` and the numbers
  kr2.* and tah.kr2noF.*, `make_system_map_v0_1.py` has om.kr2noF, and
  `make_review.py` runs each version at Tahoma without K2F_F.
- `v0_1/network.py`. PEN_KR2's citation names ref. 87. No number moved.
- The READMEs. v0_1's open items gain Kr2* + F, and LIGHT's comparison
  paragraph reports the group as a band.

### Verified
- v0_1's docs tests, 4 passed with 2 slow skipped, and the two slow checks
  with LIGHT_DOCS_SLOW=1 (the numbers above regenerate). The three PDFs
  rebuilt clean and were read where they changed.

## 2026-09-29: v0 and v0_1 compared, and v0_1's accuracy written up

The owner asked which version is more accurate, for a comparison of v0 and
v0_1 under docs/REVIEW/, and for an accuracy overview of v0_1. v0 and its
overview stay as they were.

### Added
- `docs/REVIEW/v0_vs_v0_1.tex/pdf` (3 pages) with `make_review.py`, which
  writes `review_numbers.tex`. It sets the two versions side by side on
  every test, takes the reference point apart one change at a time, and
  runs both defaults over the notes' grid. Its other numbers come from the
  two overviews' generated files.
- `docs/accuracy_overview_v0_1.tex/pdf` (9 pages) with
  `make_overview_v0_1.py`, which writes `overview_v0_1_rates.tex`,
  `overview_v0_1_fast.tex` and `overview_v0_1_slow.tex`. It is v0's
  overview for v0_1, plus Tahoma with v0_1's default network, with the
  trimer's weak law and with v0_1's own power input, and sections on the
  charge chemistry and on the input.
- `docs/system_map_v0_1.tex/pdf` with `make_system_map_v0_1.py`. It is the
  basic map with v0_1's ion cycle (Penning, conversion, recombination back
  to Kr*, ion-ion into KrF(B)), the F2-electron rows, F + F + Kr and the
  B-A band. All 37 rows are drawn and 35 nodes balance to 1e-8 of E_dep.
  F- lies below Kr + F2 on the ledger's scale, so its flows are printed and
  not drawn.
- `v0_1/tests/test_docs.py`, six tests. The map, the rate tables and the
  fast numbers are compared with fresh runs, the slow numbers and the
  review's with `LIGHT_DOCS_SLOW=1`, and every short source is held to its
  row's citation.

### Changed
- ML scores v0_1 when `LIGHT_VERSION=v0_1` (ML/CHANGELOG.md). v0's fast
  and slow overview numbers regenerate byte for byte under the switch.

### Fixed
- v0_1's CB citation still said the ledger books 80 cm^-1. It books Yu &
  Setser's 120, as species.py said. No number moved.

### Found
- v0_1 is the more accurate model, and the data show it only in the
  default network. Under the same beam at Tahoma, v0's default network (n_e
  held at 1e15, no charge chemistry) gives g0 2.05 %/cm, x0.32, because the
  beam's ions never recombine or reach KrF. Adding the Kr2* group changes
  nothing (x0.32). v0_1's default gives x0.93, and x0.95 with Kr2*.
- With the same chemistry the versions score alike. The 300 K kinetic
  datasets agree to the digit, and at Tahoma g0 is 6.11 and 6.10, Isat 1.54
  in both, and g0 Isat / Pin 0.284 and 0.285. Yu & Setser's 230/300 K
  ratio moves from -13.6 to -8.4 %, not independently. v0_1's gap is their
  consensus, and the ratio alone implies 180 +- 30 cm^-1.
- At the reference point 0.39 of the 0.51 points between v0 (5.24 %) and
  v0_1 (4.73 %) come from the default network and 0.12 from the changed
  values, mostly the B-A branch (0.17, less 0.04 from F + F + Kr).
- Tahoma's window sits at 448 K, just above the trimer law's measured
  425 K. The pure-Kr T^-1/2 law gives g0 there at x0.47, so Tahoma's gain
  counts against the weak law, as a score that sets nothing.
- Over the notes' grid the trimer band is wider than the gap between the
  versions at 10 of 15 points. The gap is wider only where v0_1's gas ends
  at 1176 K or more.
- v0_1's own power input at Tahoma gives g0 x1.03 with Te 3.8 eV, where the
  beam-sustained plasma sits at 1.6 eV. The agreement does not test the
  input.

### Verified
- From main, LIGHT (v0 and v0_1) and ML: 186 passed, 3 skipped (the slow
  docs checks), 3 xfailed. With the switch alone, before the new tests, the
  suites gave 179 passed, 1 skipped, 3 xfailed, as before it. The three
  slow checks pass with LIGHT_DOCS_SLOW=1, v0's overview (5 min) and v0_1's
  overview with the review (7.5 min).
- `ML/tests/test_light_version.py`. Six bugs planted in the switch, one at
  a time, and every one failed it.
- The three PDFs rendered and read page by page. The map has no label
  collisions, and each build is clean apart from the two overfull boxes
  the basic map already had.

## 2026-09-28 (i): v0_1, v0's working copy, with power the only input, recombination, fluorine, the trimer and the C-A band

The owner asked to ignore the drive and take the energy as uniformly and
directly deposited from a given input power, and to work on the trimer,
ion recombination, fluorine recombination and the C-A band and decay, in a
v0_1 copied from v0. v0 and the accuracy overview stay as they were.

### Added
- `v0_1/`, a copy of v0 with its imports renamed, and `run_v0_1.py`.
  `v0_1/README.md` says what changed and what it gives. `v0_1/SOURCES.md`
  holds the literature from four research agents, each download scanned
  with checktext. Every number that changed the model was re-read
  first-hand in the local sources.
- Two sources in `docs/krf_data/sources/`, the F2 pages of Cohen & Bott
  1976 (DTIC ADA024332) and Perez-Rios & Greene, PRA 98, 062707 (2018),
  both checktext CLEAN.

### Changed (v0_1 only)
- The input. The given power goes uniformly into the electron gas, and v0's
  EM-wave option is removed. No electrons are injected. n_e is free and
  quasineutral from a seed, which is the default for the basic network
  (`ne_mode` "auto"), while the notes' list keeps v0's held n_e. The scan
  runs pressure x power when n_e is free.
- Ion recombination joins the basic network. Kr2+ + e follows Shiu &
  Biondi's measured 1.6e-6 (300 K/Te)^0.55 (a new rate kind, te_pow, whose
  captured electron takes (3/2 + n) Te). Kr+ + 2Kr is 2.4e-31 cm6/s, the
  mean of three measurements, on Perez-Rios & Greene's T^-3/4. F2
  attachment, dissociation and vibration come from SIGLO as in v0's f2e
  group. Ion-ion recombination keeps Kannari's density profile and gets back
  the gas-temperature dependence he used, through (Tg/300)^-3 on its
  three-body part (the flannery spec takes an optional exponent).
- Fluorine recombination joins it, F + F + Kr -> F2 + Kr at 2.8e-34 cm6/s,
  Cohen & Bott's evaluation, held at its 300 K value.
- The trimer. Its law's stated range is Shui's 300-800 K. Kr2F* sits 0.58
  eV below KrF* + Kr (Huestis et al. 1984), where REBUILD's estimate had
  2.05 eV. C3 emits at the printed 420 nm. `Params.trimer_tg_pow` and
  `--trimer-tg-pow` run the unmeasured high-temperature band.
- The C-A band and decay. G4C emits at the printed 275 nm. A new row, G4A,
  carries B -> A at 272 nm with A(B-A)/A(B-X) = 0.035 and shares tau_B with
  B -> X. The B-C gap is Yu & Setser's 120 cm^-1. KrF(B) sits at its
  measured well depth, 5.54 eV below the ion pair (Tellinghuisen 1976),
  5.856 eV above Kr + F2 instead of 6.15.
- The report warns above 800 K (the trimer law) and above 1000 K
  (redissociation not modeled), flags a held n_e with attachment, and says
  F + F + Kr is far too slow to refill burnt-out F2.

### Found
- At the reference point v0_1 sends 4.73 % of the deposited energy out at
  248 nm, where v0's basic network gives 5.24 %. Letting n_e find its own
  value accounts for most of it. The plasma runs near 1e14 cm^-3 while F2
  lasts, so electrons quench 7 % of the KrF instead of 19 %, and more Kr*
  ends in Penning ionization, whose energy recombination returns as heat.
  The B -> A branch takes 0.17 of a point. F + F + Kr adds 0.04 and
  re-forms 0.5 % of the F2 fill. The measured recombination values, the
  ion-ion temperature law and the measured energies move it by 0.02 or less.
- The trimer's temperature law is now the largest uncertainty. At the
  reference point Shui's law gives 4.73 % and the pure-Kr T^-1/2 gives
  2.28 %. No measurement exists above 425 K.
- With recombination the gas gets hot at low pressure. It passes 1000 K at
  0.5 atm and below with 300 kW/cm^3 or more, and reaches 3461 K at 0.1 atm
  and 500 kW/cm^3, beyond every rate law's data.
- v0's tests: `pytest.approx` let rate comparisons pass unchecked (Known-bad
  above). A strict re-run found no drift in the copied REBUILD values. The
  step-cross-section check claims 0.2 % where the quadrature is good to half
  a grid cell of the Boltzmann factor (0.7 % at 0.8 eV).

### Verified
- v0_1: 67 passed, 1 xfailed (Hakuta, strict), every comparison strict. A
  mutation check planted 12 bugs in the new code, one at a time, and every
  one failed the suite.
- v0, ML and the docs are unchanged, and their tests still pass.

## 2026-09-28 (h): how accurate LIGHT is, written up with the basic network's map and every rate

The owner asked for an overview of the model's current accuracy, with the
system map and a list of the reaction rates in an appendix.

### Added
- `docs/accuracy_overview.tex` and `.pdf`, 8 pages. A scorecard of every
  test opens it. Sections follow on how the numbers were made, the KrF(B,C)
  block, KrF formation with each unmeasured link moved alone at Tahoma, the
  electron side, the system map and where LIGHT is not yet reliable.
  Appendix A lists every rate LIGHT can run with a short source. It covers
  the 24 basic rows, the 19 rows of the optional groups, the electron-impact
  rate coefficients from 0.5 to 3 eV with the tail factor at the reference
  point, v0's differing values and the formation energies.
- `docs/system_map_basic.tex` and `.pdf`, whose numbers
  `make_system_map_basic.py` writes to `system_map_basic_numbers.tex`. It
  draws the basic network's energy flow at the notes' reference point in
  the v0 map's style, adding KrF(C) with the net B-C mixing, the three
  branches of Kr* + F2, the quench boxes of B and C, the C-A band and the
  Kr2F* quench. The script checks 28 node balances, that every row of the
  network is drawn, and that the net mixing runs C to B as drawn. KrF(C) is
  green #008300, which passed the dataviz validator against the map's other
  hues with a normal-vision Delta E of 15.6 against KrF(B)'s aqua.
- `docs/make_overview.py`. It writes `overview_rates.tex` from the rows
  `network.build` returns, `overview_fast.tex` with the checks that take
  seconds, and `overview_slow.tex` with the MCMC posterior, the formation
  links, the electron closure against the kinetic loop and the swarm
  checks. A full run takes about 6 minutes, and `--fast` skips the slow
  part.
- `v0/tests/test_docs.py`, five tests. The basic map's numbers, the rate
  tables and the fast numbers are regenerated and compared byte for byte.
  The slow numbers are compared with `LIGHT_DOCS_SLOW=1`. Every row LIGHT
  can run must have a short source, and every number a short source quotes
  must appear in the row's full citation. A year may come from another
  row's citation.

### Found
The overview's numbers, all generated.
- At Tahoma, a score only, g0 is 6.11 %/cm (x0.96), Isat 1.54 MW/cm2
  (x1.10) and g0 Isat/Pin 0.284 (x1.05). Over 40 posterior draws g0 is 5.83
  [5.42, 6.63] against the measured 6.39.
- The fitted decay sets are matched to 2.3 to 6.1 % rms and the held-out
  Xu 1993 decays to 11.7 %. Hakuta's C-A decays (19 and 34 % rms), Burnham
  and Searles' tau_B (-22 %) and Yu and Setser's 230/300 K ratio (-14 %)
  stay in conflict.
- The measured range of the harpoon yield, 0.70 to 1.0, gives g0 Isat/Pin
  0.256 to 0.298 and brackets Tahoma's 0.27. The analogs for the two
  unmeasured ionic links give 0.241 and 0.266, and removing Kr2* + F gives
  0.268.
- LIGHT's electron closure lands within 0.7 % of the kinetic loop's Te and
  2.3 % of its n_e at Tahoma.
- The basic network reproduces Gadomski's eq. 5 to 1.0e-4 at six
  pressures, and the ledger closes below 1e-12 in nine runs across the
  modes.
- At the reference point the basic map sends 5.2 % of E_dep out at 248 nm,
  where the notes' list sends 16.8 %. The net B-C mixing runs C to B and
  carries 2.9 % of E_dep, with 11.5 % going to C and 14.4 % coming back.
- tau_B is the one KrF rate that is neither measured nor fitted. LIGHT
  keeps Kannari's 6.99 ns, inside the fit's 6.79 to 7.22 ns, and Eden
  measured 6.8 +- 0.2 ns. The overview says so.

### Verified
- LIGHT and ML: 112 passed, 1 skipped (the slow doc test), 2 xfailed.
  With LIGHT_DOCS_SLOW=1 the slow doc test passed in 5 minutes, so the
  seeded MCMC, the formation links and the kinetic loop reproduce byte for
  byte.
- The source guard caught a changed number in three short sources (G4C,
  C3, AK8_Kr) and a missing one (Q2F).
- Every page of both PDFs was checked by eye. The overview's prose passes
  the style grep for semicolons, dashes, colon-led clauses and (i) lists.

## 2026-09-28 (g): Velazco, Kolts & Setser 1976 read, the harpoon yield confirmed

The owner added the paper the dissertation cites, J. Chem. Phys. 65, 3468
(1976), as `docs/krf_data/sources/velazco1976.pdf` (checktext CLEAN, read
first-hand).

### Found
- The paper measures the KrF* yield directly. "the measured value of the
  KrF* formation rate constant (65 x 10^-11 cm3 molecule-1 sec-1)" (p.3476)
  is 0.88 of the KrCl* one (Table III, reproducible to +-10 %). The total
  Kr(3P2) + F2 quenching is 72 x 10^-11 (Table II), so the yield is 0.90.
- The absolute scale comes from Ar(3P2) + Kr and absorption-measured
  metastable densities, good to +-20 % (spectral response) and +-30 %
  (concentration ratio).
- The published value is an assignment: Table IV lists 1.0, and "The
  uncertainties are estimated as +-0.20".

### Changed
- `v0/network.py`: the provenance of `KRF_Y_HARPOON` (0.90, unchanged) and
  of G3 and G3X now cites the paper first-hand. No rate changed, so the
  tail table stands.

### Fixed
- Entry (f) gave the measured range as 0.63 to 0.90, taking the ratio 0.7
  of the 1976 ERDA report as the low end. That ratio came before the
  group's spectral recalibration, and the published ratio is 0.88. The
  range is 0.90 +- 0.20, that is 0.70 to 1.0. At Tahoma it gives g0 Isat /
  Pin 0.256 to 0.298 (x0.95 to x1.10), still bracketing 0.27.

### Verified
- LIGHT: 67 passed, 1 xfailed. ML: 41 passed, 1 xfailed.

## 2026-09-28 (f): the owner's four papers, and the measured KrF* yield of Kr* + F2

The owner supplied the four papers step 4 asked for, now in
`docs/krf_data/sources/`: Mezyk 1991, Chen & Payne 1979, Morozov 2008 and
Velazco's 1977 dissertation. Each was checktexted CLEAN before it was
opened. `docs/krf_data/FORMATION_LINKS_SEARCH.md` has what each measures
and settles.

### Found
- **The absolute scale of the harpoon yield.** Velazco's dissertation
  (p.148, read on the page image) puts the Kr* + Cl2 -> KrCl* reference at
  0.90 absolute, from comparison experiments with an absolutely measured
  Kr(3P2) concentration. The published KrF* "unity" is relative to it, so
  Kr* + F2 -> KrF* is 0.90. The raw ratio before the group's spectral
  recalibration, 0.7, gives 0.63, the low end.
- Chen & Payne 1979 add a third Kr2* + F2 rate (5e-10), and find the
  direct Kr2F* channel of Kr2* + F2 minor. Their KrF* quench constants are
  single-state Stern-Volmer readings that contradict their own trimer rate,
  so they are not used.
- Mezyk 1991 puts Kr2+ + SF6- recombination at about 3e-6 cm3/s near 665
  Torr, the scale of LIGHT's ion-ion rate. The anion is SF6-, and the
  paper gives no yield.
- Morozov 2008 counts all the VUV over the deposited power, so its 22 %
  shortfall is not a definition effect.

### Changed
- `v0/network.py`, basic network: the KrF* yield of Kr* + F2 is 0.90
  (`KRF_Y_HARPOON`). G3 and G3C carry 0.90 of the 8.47e-10, and the new row
  G3X sends the other 10 % to Kr + 2F, its energy to heat. K2F_BR's
  citation records Chen & Payne's finding. v0 is unchanged.
- `v0/data/kr_tail_ratio.txt` regenerated, because the network changed,
  and re-pinned.
- Tests: the harpoon seed test checks 90 % into KrF* (62/38 into B and C)
  and 10 % into G3X, and the basic network has 43 rows with every group.
- README: the formation row, the reference point, the Tahoma record and
  Next.

### What it does
At Tahoma (ML's like-for-like observable, a score only):
- g0 Isat / Pin goes from 0.298 to 0.284 (x1.10 -> x1.05 of 0.27).
- g0 goes from 6.39 to 6.11 %/cm (x1.00 -> x0.96).
- Isat goes from 1.54 to 1.53 MW/cm2 (measured 1.4 +- 0.1).
The measured yield range, 0.63 to 0.90, gives 0.246 to 0.284 and so
brackets Tahoma's 0.27.

### Verified
- LIGHT: 67 passed, 1 xfailed. ML: 41 passed, 1 xfailed.
- The dissertation's pp.147-148, Chen & Payne's p.157 and Fig. 25, and
  Mezyk's Fig. 4 were read on the page images.

## 2026-09-28 (e): the unmeasured formation links searched (step 4 of 4), and the docs brought up to date

Owner decision (2026-09-28, "Proceed with all 4 steps"): look for
measurements of the three formation links that have none.
`docs/krf_data/FORMATION_LINKS_SEARCH.md` has the search, the quotes and
the papers that would settle each link.

### Found
- No measurement exists for any of the three.
  - The absolute KrF* yield of Kr* + F2 is still relative to Kr* + Cl2.
    Setser's promised recalibration never appears in COO-2807-4 to -18.
  - Ion-ion recombination with F- has no measured KrF* yield or B/C split,
    only analogs with SF6- and Cl-.
  - Kr2* + F -> KrF* + Kr has no rate anywhere.
- **The Kr* + F2 rate LIGHT used was untraceable.** Kannari 1985 and Ueda
  1984 credit 8.1e-10 to SRI MP-76-99. That report gives Setser's 7.2e-10
  (Table 2, reaction (14), p.46, and p.67), and "8.1" is in neither of two
  OCR texts of it.
- **Kannari's Kr2* + F rate cites the same report, which has no such
  rate.** It only interprets a slow KrF* component at very low F2 as this
  reaction (p.42).
- The same report measured Kr2* + F2 at 3e-10 +- 30 % (reaction (21)),
  which agrees with Quigley & Hughes' 2.1e-10.

### Changed
- `v0/network.py`, basic network: the Kr* + F2 rate is the mean of the
  three first-hand measurements, 8.47e-10 cm3/s (7.2, 8.8 +- 0.5 and 9.4 +-
  1.0). It moves g0 Isat / Pin by +0.1 %. v0 keeps REBUILD's 8.1e-10, and
  its G3 and K2F_F citations now say what the report holds.
- `docs/krf_data/sources/`: 16 files added, all checktext CLEAN.
  - SRI MP-76-99 (with an independent OCR text), MP-76-46 and MP-78-07.
  - Setser's COO-2807-4, -8 and -12 to -18.
  - Ulrich 2012, Heindl 2011 and Dandl 2011, on what Morozov's efficiency
    counts.
- The docs are up to date with steps 1 to 4. `README.md` covers the basic
  network's table and constants, the electron closure, the decay checks,
  the reference point and scan, the Tahoma record, assumptions, tests,
  files and Next.
- Fixed: `docs/krf_data/bc_check.py` had crashed since step 1 with an
  OverflowError. Root cause: Part 3 built "Hakuta's inputs" from LIGHT's
  current constants, which since step 1 are the slow-mixing fit, so the
  K_eq and k_Kr refit diverged. It now builds LIGHT's until-09-27 set and
  Gadomski's set explicitly. `electron_quench.py` printed LIGHT's k_BC as a
  hard-coded 1.0e-10 and now reads it.

### Verified
- LIGHT: 67 passed, 1 xfailed. ML: 41 passed, 1 xfailed.
- `bc_check.py` and `electron_quench.py` run.
- The SRI report's reactions (14) and (21) and its p.42 passage were
  re-read first-hand.

## 2026-09-28 (d): exact captured energies and a two-temperature tail (step 3 of 4)

Owner decision (2026-09-28, "Proceed with all 4 steps"): fix the electron
energy bookkeeping and the tail above the Kr thresholds, with the tail
temperature taken from the kinetic EEDF at each state.

### Changed
- `v0/model.py`: two Params switches, `capture` and `eedf`. Left at "auto"
  they give the basic network the new closure and keep v0 on the notes'
  Maxwellian with (3/2)Te per captured electron, so the system map (v0) is
  unchanged.
  - `capture="exact"`: a captured electron takes <eps sigma v>/<sigma v>
    from the electron gas, which for a Maxwellian is (3/2 + d ln k/d ln Te)
    Te. F2 attachment reads its own table, R11p's Te^-1/2 template books
    Te, and a Te-independent rate books (3/2)Te.
  - `eedf="two_t"`: the five Kr rows (G1, C1a, C1b, C1c, C2) integrate over
    the Maxwellian at Te below the 9.915 eV break and a tail temperature
    T_t = r Te above it. r comes from a kinetic table in n_e/N and the F2
    fraction.
- `v0/xsec.py`: `capture_energy_table`, `two_t_table`, `TailRatio` and
  `tail_ratio_table`. The two-temperature rate separates into 1-D sums on
  the quadrature grid, so T_t = Te gives the Maxwellian table exactly.
- `v0/data/kr_tail_ratio.txt` (new, md5-pinned): made by `ML/ml/tail.py`
  (`python3 ML/run.py tail`).
  - Each point is one kinetic solution (ML/ml/eedf.py, e-e collisions on)
    at a LIGHT beam-drive state computed on the old closure. There are 66
    states, 300 W/cm3 to 3 MW/cm3 and 0 to 3 % F2.
  - r is set so the two-temperature EEDF removes the same energy through
    the five Kr rows as the kinetic one. Kr(5s) carries almost all of it
    and is matched to 3 %.
  - Each column is one series of states, and nothing is fitted. Read at
    each state's own n_e/N and F2, the table gives its r to 5 % rms where
    n_e/N >= 1e-6. States at 0.5 and 2 atm that it never saw come out
    within 3 %.
- Tests:
  - The capture identity holds on the table grid.
  - The two-temperature rate equals the Maxwellian at T_t = Te and a closed
    form for a linear-ramp cross section.
  - The table loads only at its pinned md5, and its reader is bilinear.
  - Attachment books the exact dU_e/dt. Te now rises as F2 removes the
    slow electrons, and v0's (3/2)Te case keeps its own test.
  - Each network runs its own closure.
- The quasi-neutrality test now clips densities at zero, as the model
  does. With the new closure n_e sits at 20-35 cm^-3 late in the
  afterglow instead of on the 1 cm^-3 floor, and Kr+ noise of -9 cm^-3
  (inside the solver's 25 cm^-3 absolute tolerance) showed there.

### Why the table's second axis is the F2 fraction
The plan named n_e/N and Te. Te does not carry the tail. Without F2 the
injected electrons are slowed only by elastic and e-e collisions, and e-e
lifts them over the break. Pure Kr at 300 W/cm3 has r 0.27 where 0.1 % F2
at the same n_e/N has 0.035. A polynomial surface in (n_e/N, Te) fitted to
the states put r = 1 in corners that no state reaches. n_e/N and the F2
fraction are the states' own axes.

### What it does
At Tahoma (the fit set under the prototype drive):

| closure | Te (eV) | n_e (cm^-3) | g0 (%/cm) | formation |
|---|---|---|---|---|
| Maxwellian, (3/2)Te capture | 1.126 | 4.66e14 | 6.83 | 0.341 |
| two-temperature, exact capture | 1.595 | 6.38e14 | 6.39 | 0.335 |
| ML's kinetic loop | 1.606 | 6.53e14 | 6.33 | 0.334 |

The kinetic loop runs LIGHT on its Maxwellian closure with every electron
row scaled by the Boltzmann factors. LIGHT on its own new closure lands
within 1 % of it in Te, g0 and formation, and within 3 % in n_e. The
measured g0 is 6.39 %/cm.
- ML's like-for-like formation g0 Isat / Pin goes from 0.306 to 0.298
  (x1.13 -> x1.10 of the measured 0.27).
- Isat goes from 1.48 to 1.54 MW/cm2 (measured 1.4 +- 0.1).
- Pure Kr at 330 W/cm3: thermal Kr excitation 0.14 -> 0.060 per ion pair,
  VUV yield 0.595 -> 0.570 (Morozov 0.42 +- 0.05). The kinetic loop puts
  the Kr(5s) excitation at about half of LIGHT's new value. At this density
  with no F2 the bulk itself is not Maxwellian (LIGHT's Te 0.88 eV against
  the kinetic 0.78-0.82), which no tail closure represents.

### Verified
- LIGHT: 67 passed, 1 xfailed (Hakuta).
- ML: 41 passed, 1 xfailed (Morozov).

## 2026-09-28 (c): the measured formation steps (step 2 of 4)

Owner decision (2026-09-28, "Proceed with all 4 steps"): put the measured
links of KrF formation into LIGHT (ML/README.md, Results for KrF
formation).

### Changed
- `v0/network.py`, basic network only (v0 keeps REBUILD's copies):
  - Kr* + F2 now forms KrF in B and in C. G3 takes 62 % of the 8.1e-10
    cm3/s and the new row G3C takes 38 % (Yu & Setser 1990 Sec. III.A,
    k_B/k_C about 1.6, read first-hand). The KrF* yield stays 1, relative
    to Kr* + Cl2 (Velazco, Kolts & Setser 1976).
  - The kr2 group carries measured values (`basic_kr2_rows`). Kr* + 2Kr
    -> Kr2* is 1.7e-32 cm6/s on (Tg/300)^-0.64 (Ueda 1984 Table V, Wieme &
    Lenaerts 1981). The Kr2* radiative rate is 3.7e6 1/s, the measured 270
    ns (Morikawa 1989 and five more at 250-353 ns). Kr2* + F2 is Quigley &
    Hughes' 2.1e-10, split evenly between KrF(B) and Kr2F* as Kannari
    assigned it.
  - Kr2* + F and Kr2* + Kr2* keep REBUILD's values, since no measurement
    replaces them.
- Tests:
  - The harpoon seed test checks the 62/38 split in the basic network and
    all-B in v0.
  - A new test pins the basic network's kr2 rows and v0's REBUILD copies.
  - The basic network has 42 rows with every group on.

### Why
Each link is its own measurement, and no Tahoma number chose any value.
At Tahoma with LIGHT's constants (ML's like-for-like observable):
- g0 Isat / Pin goes from 0.328 to 0.306 (measured 0.27, so x1.22 ->
  x1.13).
- g0 goes from 7.26 to 6.81 %/cm (measured 6.39).
- Isat goes from 1.49 to 1.48 MW/cm2 (measured 1.4 +- 0.1).

The nascent C share does it. With the slow B-C mixing, little of the KrF
formed in C reaches the laser. The Kr2* route is not a loss once it has the
measured 270 ns lifetime, because Kr2* reacts before it radiates (KrF made
per ion pair 1.845 -> 1.841).

### Verified
- LIGHT: 55 passed, 1 xfailed (Hakuta). The system map (v0) is unchanged,
  and the drift gate still checks all 24 REBUILD copies.
- ML: 36 passed, 1 xfailed. ML's measured chain is now LIGHT's own chain,
  and ML's Tahoma drive runs the kr2 group by default.

## 2026-09-28 (b): the joint-fit B/C constants adopted (step 1 of 4)

Owner decision (2026-09-28, "Proceed with all 4 steps"): adopt ML's fitted
KrF(B,C) constants in the basic network.

### Changed
- `v0/network.py`, basic network only (v0 unchanged):
  - k_BC 1.2e-10 -> 1.05e-11, K_eq 5.0 -> 1.10, tau_C 60 -> 50 ns;
  - B and C now carry their own quench constants. B: Kr 1.0e-14 (the fit's
    lower bound, so no B + Kr loss beyond B -> C transfer), 2Kr 1.90e-30,
    F2 4.75e-10. C: Kr 8.37e-13, 2Kr 3.79e-31, F2 3.04e-10;
  - k_e 4.5e-8 -> 6.6e-8, Trainor & Jacob's product re-read with the new
    Stern-Volmer lifetime (5.8 ns at 190 Torr);
  - the basic network has its own G5 row. v0 keeps REBUILD's 6.0e-31;
  - rows whose values are no longer REBUILD copies lost their `rebuild=`
    tag, as Q2F's and Hakuta's rows already had, so the drift gate still
    guards the true copies.
- `v0/two_level.py`: `Consts.from_light` reads the B and C rows separately.
- Tests:
  - Eden's B-only decays pass (x0.96, 0.95, 0.98), so the strict xfail came
    off;
  - Hakuta's C decays are the new strict xfail (19 % and 34 % rms low);
  - Xu's decay score reads x0.91 and x0.97;
  - the eigenvalue check fits the decay below 1e-3 of the start, because the
    fast eigenvalue at 50 Torr / 6 Torr F2 is only 2.7x the slow one.

### Why
The fit uses only measurements of these processes (Eden, Kansas State,
Yu & Setser). It predicts two Tahoma numbers no fit used: g0 x0.57 ->
x1.14 and Isat x2.11 -> x1.06 (1.49 against 1.4 +- 0.1 MW/cm2).

### Verified
- LIGHT: 54 passed, 1 xfailed (Hakuta).
- ML: 36 passed, 1 xfailed. ML's scorecard now names the old constants
  "LIGHT to 09-27" and keeps Gadomski's published set explicitly.

## 2026-09-28: KrF formation checked against its measurements (ML, no code change)

ML made KrF formation the target (`ML/ml/formation.py`,
`python3 ML/run.py formation`). Nothing in LIGHT changed. Adopting the
measured chain is the owner's call (ML/README.md, Next).
- **Tahoma's formation efficiency is g0 Isat / Pin, not LIGHT's formation
  rate x 5 eV / Pin.** With a probe flux on LIGHT and Isat fitted as
  Tahoma fitted it, LIGHT's formation is x1.22 of the measured 0.27 for
  every B/C set.
- **Two links of LIGHT's formation chain disagree with measurements.**
  38 % of the KrF from Kr* + F2 forms in C (Yu & Setser 1990), where LIGHT
  books all of it to B. The kr2 group's Kr2* radiative rate, 8e7 1/s
  (Kannari via REBUILD), is 20x the measured 3.7e6 (270 ns).
- **With every measured link at its measured value**, formation is x1.13
  (fitted B/C set, g0 x1.07, Isat x1.06). The rest lies in unmeasured
  links: the absolute KrF* yield of Kr* + F2, and the KrF* yield and B/C
  split of ion-ion recombination.
- New sources in `docs/krf_data/sources/`: Setser COO-2807-3, -11 and -19,
  and Avco DTIC ADA047684.

## 2026-09-27: the electron side checked against the kinetic EEDF (ML, no code change)

ML solved the steady Boltzmann equation of the beam-sustained Kr/F2 plasma
at LIGHT's Tahoma state, with electron-electron collisions
(`ML/ml/boltzmann/ee.py`, `ML/ml/eedf.py`, `python3 ML/run.py electron`).
It found two open limitations in LIGHT's electron side. Nothing in LIGHT
changed. Adopting the fixes is the owner's call (ML/README.md, Next).
- **The Maxwellian tail above the Kr thresholds is not physical.** The bulk
  is Maxwellian to 4 %, which the e-e collisions at n_e/N = 2e-5 enforce.
  Above 9.9 eV, excitation empties the tail faster than they refill it.
  At the same mean energy the rates fall: Kr(5s) x0.07, 5p x0.005, 4d+6s
  x6e-4, 5d x1e-5, ionization x2e-7.
- **The captured electron's (3/2)Te overstates the cooling.** F2 attaches
  below 1 eV, so an attaching electron carries off 0.54 eV where LIGHT
  books 1.69 eV. At Tahoma that is 21 % of the secondary-electron heating.
  The Maxwellian-exact energy per event is (3/2 + d ln k/d ln Te) Te,
  which is Te for R11p's Te^-1/2 template.
- With both fixed LIGHT agrees with the kinetic solution (Te 1.60 vs 1.60
  eV, n_e within 5 %). Tahoma with the fitted B/C set: g0 7.28 -> 6.78
  %/cm (measured 6.39), formation 0.341 -> 0.333 (measured 0.27).

## 2026-09-27: the fitting moves to ML/, and the closed form into v0

The owner made `Kinetic_Modeling/ML/` to train models to the available data.

### Changed
- `docs/krf_data/bc_lab.py` is now `v0/two_level.py`, a LIGHT module with
  package imports. `bc_check.py` and `electron_quench.py` import it from
  there.
- `docs/krf_data/bc_infer.py`, `bc_nn.py` and `accuracy.py` moved to
  `../ML/` as `ml/fit.py`, `ml/nn.py` and `ml/score.py`, running on a
  dataset registry. They are removed here, and the README and write-up
  point to `python3 ML/run.py`.

### Added
- `v0/tests/test_two_level.py` (6 tests): the closed form against the ODE,
  within 1e-5 of the seed at 300 K (four conditions) and from a 500 K start,
  plus `from_light` reading the network.

### Verified
- 54 tests: 53 pass, 1 is a strict expected failure.
- `bc_check.py` and `electron_quench.py` print the same numbers as before.

## 2026-09-27: accuracy against held-out and device data

Cooper asked how the model compares with other experimental data, and what
its accuracy is. No LIGHT constant changed.

### Added
- `docs/krf_data/accuracy.py`. It scores three B/C sets against fitted,
  held-out and device data:
  - the sets are LIGHT now, Gadomski 1992, and the set inferred from Eden +
    Kansas State + Yu & Setser;
  - each set's k_e is re-derived with its own lifetime.
  - `--posterior` propagates the MCMC posterior into the Tahoma gain.

### Found
- **Held-out Xu 1993 Figs. 5-7** (pure Kr/F2, 2 Torr F2, 50-200 Torr Kr):
  10.5 / 12.7 / 11.6 % rms, inside the Kansas State lab's own 25 % scatter.
- **Misses shared by all sets:**
  - Hakuta, 29-36 %;
  - Burnham & Searles' tau_B, -22 %;
  - Yu & Setser's 230/300 K ratio, -14 % (the ledger gap is 80 cm^-1,
    their two-temperature value 180 +- 30).
- **Tahoma small-signal gain** (6.39 %/cm, prototype beam drive, a score
  only):
  - LIGHT now 3.67 (x0.57), Gadomski 3.82 (x0.60), inferred 7.28 (x1.14);
  - over 40 posterior draws, 7.15 [6.44, 7.89];
  - no fit used any device data.
- **Formation efficiency** is 0.341 against 0.27 (x1.26) for every set, a
  miss upstream of B and C.

### Fixed (evidence scripts only)
- `k_e_for` in `accuracy.py`:
  - symptom: sets built by `bc_infer.make` carry k_e = 0, so `sv_tau`
    returned 0/0 = NaN and the Tahoma run crashed on a NaN Jacobian;
  - fix: the lifetime is independent of k_e's size, so a unit placeholder
    is used.

## 2026-09-27: Yu & Setser read, and the B/C channels inferred from the data (MCMC and a neural network)

Cooper asked whether a neural network could find the proper values for a
few channels, and added Yu & Setser (as `Kinetic_Modeling/yu1990.pdf`). No
LIGHT constant changed. The adoption decision is the owner's.

### Found
- **Yu & Setser (J. Phys. Chem. 94, 2934, pp. 2940-2941, first-hand).**
  - Equilibrium C/B emission ratio at high pressure, corrected for B -> A:
    I_r(300 K) = 0.16 +- 0.02, giving K_eq 1.7 +- 0.2 with tau_C 75.2 ns.
  - dE_BC: -110 +- 20 cm^-1 from the absolute ratio, -180 +- 30 cm^-1 from
    the 230/300 K change, consensus -120 +- 50.
  - They need more than 1000 Torr of Ar for B/C equilibrium, take Ar mixing
    as about 1e-11, and give a coupled lifetime of 15 +- 5 ns.
- **Consistency.** Least-squares fits with each dataset left out show that
  only Eden + Kansas State + Yu & Setser fit together within their errors
  (chi2 per point 0.07, 0.18, 0.00). Hakuta is inconsistent with them (14.6).
- **MCMC on those three.**
  - k_BC 9.8e-12 [4.6e-12, 1.7e-11], K_eq 1.07 [0.84, 1.38], tau_C 50
    [41, 63] ns.
  - F2 quench is 4.6e-10 for B and 3.0e-10 for C, which resolves the old F2
    spread as two states.
  - Lasing B share at 1 atm: 0.69 [0.58, 0.79], against LIGHT's 0.27.
- **Neural network (amortized posterior estimation, 30 000 simulations).**
  - On those three datasets it reproduces the MCMC (k_BC 7.5e-12, K_eq 1.15,
    tau_C 49 ns, share 0.73), with calibrated uncertainties (coverage
    0.64-0.66 of the ideal 0.68).
  - It cannot pin the Kr two-body quench (R^2 0.08-0.24) or the B-state
    2Kr (0.35).
  - Fed all four datasets, it returns a confident compromise that misfits
    Kansas State (chi2 8.6 per point).

### Added
- `docs/krf_data/bc_infer.py`: leave-one-out consistency and an ensemble
  MCMC (Goodman-Weare, numpy only), about 2 min.
- `docs/krf_data/bc_nn.py`: the neural-network test (numpy MLP with Adam
  and a Gaussian head), about 1 min.
  - The spurious Accelerate matmul warnings are silenced there.
  - Finiteness is asserted.
- `sources/Yu_Setser_JPC94_2934_1990.pdf`, a copy of Cooper's yu1990.pdf
  (checktext CLEAN).
- Docs: the write-up's new top section, the README "Next" (owner decision),
  and the catalogue source list.

### Verified
- The Yu & Setser formula returns their I_r = 0.160 with their constants.
- The network's held-out R^2 and coverage are as listed above. Its answer
  on the three consistent sets matches the MCMC within the intervals.
- LIGHT's 48 tests are unchanged (no code under `v0/` changed).

## 2026-09-27: Hakuta's measured tau_C and k_BC, and K_eq(Tg) by detailed balance

Cooper added Hakuta et al., JJAP 29, 1194 (1990), and a file meant to be
Yu & Setser, J. Phys. Chem. 94, 2934 (1990). The second is the article
before theirs, and its last page carries their abstract.

### Changed (basic network)
- **Hakuta's measured values replace two assumptions.**
  - tau_C: 75 -> 60 ns.
  - k_BC: 1.0e-10 -> 1.2e-10.
  - Source: Hakuta 1990 Table I, read in full. Gadomski had assumed both.
- **C -> B by detailed balance.** The row is now k_BC / K_eq(Tg), with
  K_eq(Tg) = 5.0 exp(dE/k (1/Tg - 1/300)). dE is the gap the ledger books,
  80 cm^-1 from species.py. Before, the kinetics used a fixed 5.0.
- **K_eq and the quench constants are unchanged.**
  - K_eq stays 5.0 at 300 K, and the Kr, 2Kr and F2 quench constants stay.
  - Gadomski's reduction redone with Hakuta's inputs gives K_eq 5.4 +- 10 %
    and k_Kr 3.6e-12, both inside his published errors.
- **k_e is unchanged at 4.5e-8.** The two-level re-reading gives 4.6e-8.
- **Effects.**
  - Decays: x0.97, x1.01 and x1.01 at 50, 100 and 500 Torr (were x0.89,
    x0.95, x1.00).
  - Reference-point 248 nm: 3.47 -> 3.52 %.
  - Lasing B share at 1 atm: 0.29 -> 0.27 at 300 K, 0.23 -> 0.24 at 460 K.

### Added
- `model.py` rate kind "balance" (k_fwd / K(Tg), with the gap from the
  species energies).
- `test_c_to_b_follows_detailed_balance`.
- `digitized/hakuta1990_fig2_F2.csv` and `hakuta1990_fig3_Kr.csv`.
  - The page was rendered at 600 dpi and the data circles fitted as equal
    disks.
  - The printed F2 slope comes back as 2.38e-10 against 2.5 +- 1e-10.
- `bc_check.py`:
  - Hakuta scoring;
  - Part 3, Gadomski's reduction redone (the basis of the constants).
- `tahoma_check.py` (parked): the K_eq column now shows the ledger-gap value.
- `sources/`: the two PDFs Cooper added, both checktext CLEAN.

### Found (open)
- **Eden against Hakuta on B -> C mixing.**
  - Hakuta's Fig. 1 B-X tail at 9 Torr, which tracks the C state, shows
    fast coupling.
  - Eden's B-only decays allow at most about 5e-11.
- **Yu & Setser against Kansas State on K_eq.**
  - Yu & Setser: 2-3 (via Xu), and dE -120 +- 50 cm^-1 (abstract).
  - Kansas State: 5.
  - At 2.5 the lasing share at 1 atm would be 0.35.
  - Needs the full Yu & Setser article.
- **Hakuta against Kansas State on absolute decay rates.** They differ by
  about 1.8x at 20-50 Torr. Adding a B -> A contribution to Hakuta's
  280-nm signal (Julienne & Krauss) makes the fit worse, not better.

### Verified
- 48 tests: 47 pass, 1 is a strict expected failure (Eden). Six mutants of
  the new rate kind and constants are all caught.
- LIGHT against the closed form: 3e-6 of the seed at 300 K, and 7e-8 for a
  seed run started at 500 K.
- The reference point and scan were re-run and the README numbers updated.
- Re-read first-hand:
  - Hakuta 1990, all pages;
  - Yu & Setser's abstract;
  - Xu 1993 p.2598 on Yu & Setser.

## 2026-09-27: KrF(B), KrF(C) and the electron quench against their measurements

Owner direction: kinetics first, and no drive or device modules in LIGHT
(the e-beam drivers live in SPECK and SPECKLE). Focus on the B and C
states and the electron quench rate.

### Found
- **The electron quench is tied to the B/C constants.**
  - Trainor & Jacob measured the product k_e tau / k_att = 2.9e-7 s (Duzy
    et al., AVCO 1980, p.21).
  - LIGHT's 4.5e-8 is Kushner, Hanson & Schneider's one-state (KrF(B) only)
    reading of that product.
  - Read with this network's two-level Stern-Volmer lifetime (8.3 ns at
    190 Torr) it gives 4.6e-8.
  - It moves to 8.0e-8 if k_BC is 1e-11, and spans 3.5e-8 to 1.1e-7 with
    the electron temperature assumed for F2 attachment.
- **k_BC = 1.0e-10 is an assumption, not a measurement.**
  - Eden et al. 1978's B-only decays, the one full-text direct test, allow
    at most about 5e-11. LIGHT misses them by 25 % rms.
  - Hakuta et al. 1990 (abstract only) gives 1.2e-10.
- **The Kansas State decays cannot separate k_BC from K_eq.**
  - With all two-level constants free, any k_BC from 3e-12 to 3e-11 fits
    Eden and Kansas State together within 1.5-6 %.
  - K_eq moves with k_BC (0.4 to 2.2 over that range).
  - The B share of KrF at 1 atm, which sets the gain, is 0.29 in LIGHT and
    0.38-0.87 across those fits.
- **tau_B and tau_C are calculated.** The one tau_C measurement is
  60 +- 10 ns (Hakuta, abstract).

### Added
- `docs/krf_data/`:
  - `BC_AND_ELECTRON_QUENCH.md`
  - `bc_lab.py`: the two-level model in closed form, equal to LIGHT's run
    to 3e-6 of the seed
  - `bc_check.py`
  - `electron_quench.py`
  - `digitized/`: a search agent's figure readings. Eden Fig. 3 and
    Gadomski Fig. 2 were checked by eye.
- `test_electron_quench_is_trainor_jacob_under_this_b_c_model`: the QE row
  must be within 10 % of the two-level reading.
- `test_basic_network_against_eden_b_only_decays`: a strict expected
  failure until the B/C constants are revised.

### Changed
- **Provenance text only, no numbers.** The `basic_rows` docstring, the
  G4C, BC, CB and QE row comments, the `Params` and CLI labels, and the
  README basic-network section.
- **Decay-check caveat.** The README now notes that Kansas State's own
  figures scatter by up to 30 % at one condition.
- **Correction.** The earlier entry today says the network carries "the
  constants they and Xu 1993 measured". That is wrong for tau_B and tau_C
  (calculated) and k_BC (assumed). The Known-bad list now says so.
- **README.**
  - The Tahoma section is parked.
  - "Next" is rewritten: settle k_BC and K_eq from Hakuta 1990 and Yu &
    Setser 1990, then change the network once, with k_e re-derived.

### Verified
- 47 tests: 46 pass, 1 is a strict expected failure. Three mutants are each
  caught by the new tests: QE x2, k_BC 1e-11 alone, and K_eq 1.5 alone.
- Re-read first-hand:
  - KHS 1989, all three pages
  - Nishioka's abstract and conclusion
  - Trainor & Jacob's equations
  - Eden 1978, APL pp. 733-735
  - Gadomski's k_BC assumption and Fig. 2
  - Xu's K_eq ranges
  - Kannari's 80 cm⁻¹
- Not re-checked:
  - the agents' other figure readings
  - Hakuta's and Yu & Setser's abstracts
  - the Trainor & Jacob F2-subtraction and ionization-rate estimates

## 2026-09-27: measured Kr/F2 data and the first Tahoma test

### Added
- `docs/krf_data/README.md`: the explicit Kr/F2 measurements a literature
  search found, each with its page and a checked or agent label. Tahoma
  (99.6 % Kr/F2, e-beam, g0 and formation efficiency measured) is the best
  independent test. The file also covers the Xu/Gadomski and Eden decays,
  Trainor & Jacob's electron quench and the measured Kr/F2 rate constants.
- `docs/krf_data/sources/`: 13 source PDFs (36 MB), all CLEAN under
  checktext. The 52 MB NRL annual report is cut to its KrF section.
- `docs/krf_data/tahoma_check.py`: the basic network against Tahoma under
  the prototype e-beam drive of `docs/survey_2026-09-27`, with an optional
  K_eq(Tg) from the B-C gap Gadomski 1992 prints.
- README: a Tahoma section, and a "Next" section that states the one path.

### Result (LIGHT itself unchanged)
- g0 3.60 %/cm against 6.39 measured (x0.56), and 4.53 (x0.71) with a
  300 cm⁻¹ gap. Formation efficiency 0.341 against 0.27 (x1.26).
- At that state the B/C manifold loses 2.3e8 s⁻¹. Kr quench takes 37 %,
  the trimer 29 %, radiation 19 %, electrons 9 % and F2 6 %. Only 23 % of
  the manifold is in B.

### Verified
- The Tahoma numbers were re-read in the report text (Table 3-3, pp. 2-2,
  2-3, 3-7, 3-16, 4-1). The beam-stop-in pairing is the report's own
  (g0·Isat = 0.064 × 1.4, over 330 kW/cm³).
- Four agent quotes were re-read in their sources: Trainor & Jacob's
  2.0e-7, the Klimek stopping-power sentence, Eden's 8.6e-12 and 9.7e-31,
  and Gadomski on KrF2 photolysis.
- `tahoma_check.py` reproduces the earlier prototype numbers. No code
  under `v0/` changed, and the 45 tests pass.

## 2026-09-27: the basic network (new default)

Owner direction: cut back to the basic processes (excitation, ionization,
scattering, emission, heat loss), use explicit KrF data, and perfect the
Kr/F2 gas before adding anything.

### Added
- `network="basic"` (`--network basic`), now the default. It has 15
  processes in 22 rows:
  - the notes' six reactions and the three v0 closures;
  - KrF(C) with its 75-ns emission near 275 nm;
  - B/C mixing by Kr (k_BC 1.0e-10, K_eq 5.0);
  - measured quenching of KrF(B,C) by F2 (3.0e-10), Kr (4.0e-12) and
    electrons (4.5e-8), plus the trimer from both states;
  - measured quenching of Kr2F* by F2 (1.2e-10).

  The KrF part is the two-level model of Gadomski, Xu, Setser & Nelson 1992
  (eq. 5), with the constants they and Xu 1993 measured. The electron quench
  is Kushner 1989's reassessment of Trainor & Jacob's measurement.
- The v0 network (the notes' list as written) stays available as
  `--network v0`.
- The `quench` and `bc` groups are refused with the basic network, which
  already carries measured versions of them.
- New tests:
  - the basic decay equals Gadomski's analytic eq. 5 at six pressures
    (rel 1e-3);
  - a score against Xu's measured decays;
  - the basic row list;
  - a v0 closure case.

  45 tests in total.

### Changed
- The reference point changes a lot:
  - 248-nm emission is 3.47 % of E_dep (v0: 16.8 %);
  - KrF formation is unchanged (31.0 % vs 32.5 %);
  - only 11 % of the KrF formed radiates at 248 nm. The trimer takes 34 %,
    Kr quenching 32 %, electron quenching 15 % and F2 quenching 5 %.

  Across the notes' grid the 248-nm share is 0.8 to 7.5 %, where v0 gave
  1.9 to 30 %.
- `docs/make_system_map.py` pins `network="v0"`. The system map still draws
  the notes' list as written.
- The report prints the network.
- The README leads with the basic network. Its v0 sections are labelled.
- `IMPROVEMENT_PLAN.md` is marked as a reference menu for the chosen path.

### Verified
- Measured KrF decays (pure Kr + 1 Torr F2, 300 K):

  | Kr (Torr) | basic | v0 |
  |---|---|---|
  | 50 | x0.89 | x2.89 |
  | 100 | x0.95 | x2.41 |
  | 500 | x1.00 | x1.13 |

  The 50-Torr residual equals that of Gadomski's own published fit.
- The KrF decay matches Gadomski's analytic eq. 5 at eight pressure and F2
  combinations to four decimal places.
- 45 tests pass, and the drift gate holds for the basic rows copied from
  REBUILD (K32, K33, K34, K67 values).
- The mutation check (21 deliberate bugs) catches 21/21. That includes
  reverting K_eq to 1.47 and a wrong tau_C. Homebrew 3.14 runs the basic,
  v0 and omissions modes clean with warnings as errors.

## 2026-09-27: improvement survey and plan

### Added
- `IMPROVEMENT_PLAN.md`: what an improved LIGHT should do, from a survey of
  the rest of Kinetic_Modeling (owner ask: the reactions and missing pieces,
  recombination among them, that would bring LIGHT closest to experiment).
- `docs/survey_2026-09-27/`: three survey notes (reaction sources,
  experimental data, missing components) and seven evidence scripts. The
  README there lists the claims re-checked first-hand.

### Findings (measured with LIGHT, scripts in the survey folder)
- The rate tables in `docs/Recent_KrF_Papers` hold about 429 distinct
  reactions, about 460 with the repo's codes, 130 of them in the Kr/Ar/Ne/F
  family and missing from REBUILD.
- The KrF B/C structure is the largest lever. Kannari's K_eq 1.47 with
  constants that Xu and Gadomski extracted under K_eq 4 to 5 misses their
  pure Kr/F2 decays ×1.5 to ×2.9. The measured structure reproduces them to
  ×1.00 and halves the 248-nm share (5.94 → 2.93 % at 1 atm).
- Recombination decides the answer only under an e-beam drive. There, the
  ionic channel forms about half of KrF(B), and switching it off cuts
  248-nm output ×9 to ×14. Under v0's drive it moves the output by less
  than one point.
- Kannari's Kr* + e → Kr** coupling, omitted by REBUILD as negligible at beam
  densities, cuts the 248-nm share 37 % at n_e 1e15 under v0's drive.
- The prototype beam drive exposed a solver failure: attachment drives free
  n_e to its floor and the Jacobian turns non-finite (one case at 2 atm).
- Outside LIGHT, not changed:
  - REBUILD's NK6 carries Chen's misprinted R23 and runs ×9.3 high in Ne
    fills.
  - BC_Kr names Gadomski as its primary, but Gadomski assumed the value, and
    the measured 1.2e-10 is ×0.24 of the row.
  - VREL_Kr's primary has no Kr data.
  - phantom11's Stevefelt row is broken.

### Verified
- The load-bearing literature values were read first-hand:
  - Xu 1993 p.2591 and p.2594;
  - Kannari 1985 Table I p.4311;
  - chen2025.py against network.py for NK6;
  - NRL_REPLICATION.md for the Mandl ratios.
- LIGHT's own decay runs reproduce the survey's hand calculation digit for
  digit.
- The evidence scripts run from their new location and give the numbers the
  plan quotes. No LIGHT code changed, so the suite (36 tests) is unaffected.

## 2026-09-27: the system map

### Added
- `docs/system_map.pdf` and its source `docs/system_map.tex`: a two-page
  A4-landscape map of the v0 kinetics at the reference point (owner ask: a
  system map of the kinetics, with power in feeding the uniform initial
  gas, then the reactions, then the outputs). Page 1 runs from the inputs
  through the uniform-gas container (electron gas, G1 to G6, C1 to C3, the
  Kr*, KrF(B), Kr2F* and Kr+ stores) to the fate ledger, with band widths
  proportional to energy flow and heat and F-atom losses as tags. A dashed
  box gives the 248-nm share with each optional group on. Page 2 is the
  key: how to read it, the reaction table with sources, the formation
  energies, the reference state and the not-in-v0 table.
- `docs/make_system_map.py` writes `docs/system_map_numbers.tex` from the
  model (reference point plus the seven omission runs, about 8 s). Before
  writing, it checks energy conservation at every box on the map to 1e-8
  and the ledger closure to 1e-12. No number on the map is typed by hand.
- `v0/tests/test_docs.py`: the numbers file must regenerate byte for byte,
  so a model change cannot leave the PDF stale.

### Fixed
- The first draft drew bands at 0.03 cm per percent, making the 100 %
  power band 3 cm wide. The scale is 0.03 mm per percent.

### Verified
- latexmk builds with 0 errors and 0 overfull boxes. Both pages were
  rendered at 110 and 200 dpi and inspected for overlaps and clipping.
- The five band hues pass the dataviz palette validator (all pairs, white
  surface). The two WARNs it reports, red against aqua in the 6 to 8 CVD
  band and amber and aqua below 3:1 contrast, are covered by a text label
  on every band and tag.
- The closest printed value to a rounding edge is 32.549, a 3e-5 relative
  margin against about 2e-6 cross-interpreter drift, so the byte-exact
  test is stable.
- `python3 -m pytest LIGHT/v0/tests -q`: 36 passed.

## 2026-09-25: v0 created (0.1.0)

### Added
- `LIGHT/`, the lightweight sibling of `SPECK/`, and `v0/`, Version 0 of
  `../KrF_model_notes_with_Gary.rtf`: a 0-D Kr/F2 model driven by uniform
  power into the electrons (100 to 500 kW/cm^3, microsecond pulse), an
  injected n_e (1e14 to 1e16 cm^-3, held or free) or an EM wave, at 0.1 to
  2 atm. Output: the notes' "percent to Kr+, KrF(B), Kr*, heat, elastic
  losses and emission losses", as flows and as a fate-at-t_end ledger that
  closes to about 1e-16.
- Network as data: the notes' six reactions (G1 to G6) plus three closures
  the output list needs (C1 5p/4d/higher cascade to 5s, C2 ionization, C3
  Kr2F* emission). Six optional groups (quench, kr2, ions, electron, bc,
  f2e) hold the rest of SPECK's pure Kr/F2 chemistry, off by default: 27
  rows, plus 4 KrF(C) copies when bc is on. `--omissions` measures what
  each costs.
- Maxwellian k(Te) = <sigma v> from vendored Biagi-v7.1 Kr and SIGLO F2
  files, the exact Maxwellian elastic loss (2m/M)(1 - Tg/Te)<sigma v eps>,
  and a Lorentz-gas EM absorption kernel with e-n and e-i collisions.
- `run_v0.py` CLI (single point, `--scan` over the notes' ranges,
  `--omissions`, `--add`, `--em-E/--em-f`, `--ne-mode`, `--tau-sp`, CSV and
  JSON output) and the README.

### Findings (v0 as written, reference point 1 atm, 0.5 % F2, 300 kW/cm^3,
1 us, n_e 1e15)
- 92.5 % of the power excites Kr (40 % of that through 5p and higher) and
  10.4 % ionizes it directly. KrF(B) formation is 32.5 % of E_dep in 5-eV
  quanta, and 248-nm spontaneous emission is 16.8 %.
- F2 burns out at 830 ns. The notes' drive (0.1 to 0.5 J/cm^3) exceeds the
  F2 inventory of a 0.5 % fill below 1 atm, and Kr* then ends in Penning
  ions (Kr+ up to 87 % at 0.1 atm).
- Switching on all optional groups cuts 248 nm from 16.8 % to 5.9 %.
  Quench and B/C mixing cost the most. v0's "more F2 is better" comes from
  the missing KrF(B) + F2 quench: with it, the share peaks near 1 % F2.
- Flags for the notes (README, "Where the notes and the literature
  differ"): tau_sp is 6.99 ns in Kannari, not 4 ns (4 ns raises 248 nm to
  20.7 %). Reaction 5 is trimer formation, not harpooning. Kannari prints
  Kr* + Kr* as Penning, not associative. The higher excitations are
  needed. And Q is in cm^-3 s^-1.

### Verified
- `python3 -m pytest LIGHT/v0/tests -q`: 35 passed in about 11 s under
  /usr/bin/python3 3.9.6 (scipy 1.13.1). That covers closed-form
  Maxwellian integrals, detailed balance, the Drude limit, balance of all
  42 rows, the drift gate against REBUILD 0.1.40 (24 rates plus 9 species
  energies), closure in every mode, gas thermal energy against heat,
  analytic branching limits, the cited rate-law values and rtol
  convergence.
- Mutation check: 19 deliberate bugs, one per fresh copy, all 19 caught.
- Homebrew python3 3.14.6 (scipy 1.18.0): reference, scan, omissions, EM
  and free-n_e runs clean with `-W error::DeprecationWarning -W
  error::RuntimeWarning`.
- Scan: 45 points, worst closure 4.4e-15.
