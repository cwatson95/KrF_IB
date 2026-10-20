# REBUILD changelog

All notable changes to the Orestes-class engine rebuild.
Dated entries, newest first (Added / Changed / Fixed / Removed), per the
repo-wide changelog policy (`Agent/CLAUDE.md`).

**Read this file before editing anything in this directory** — the
Known-bad list is permanent.

## Known-bad — do not reintroduce

- **Never hand-write an RHS term.** The RHS is assembled from `network.py`
  rows; a new process is a new ROW (balance-checked mechanically), never an
  edit to `chemistry.py`'s loop. The r115 class of defect (sign error worth
  40 % of dArF, invisible to parity) is impossible only while this holds.
- **Never route a radiative row's energy off its hand-stated quantum when a
  cavity is attached.** The upper level's spontaneous row must radiate the
  EXACT cavity photon energy (5.00 stated vs 4.9997 exact at 248 nm leaked
  0.3 meV/event and failed the per-event identity gate, found 2026-08-04).
  `CompiledNetwork.set_photon_energy` recomputes it; keep it that way.
- **Never filter a literature-table completeness pass by the species
  already in the network.** Rows of untracked species are exactly where
  whole missing MANIFOLDS hide: the ArKr*/ArKrF* chain (a ~4.6e8 s^-1
  KrF* sink class, x3.2 the B radiative rate) survived FIVE first-hand
  table passes because every pass looked up rows by known species
  (found 2026-08-04 via Xu 1993's citation trail). A completeness pass
  walks the SOURCE's table row by row and disposes of every row —
  adopt, or record why not.
- **Never establish a row's existence or absence from network.py source
  text.** The CHEN block is GENERATED from chen2025.TABLE_III — grep
  finds nothing for it. That is how 0.1.10 'discovered' four absences
  that existed (Ar* harpoon, Ar-ion+F- recombination, Ar2+ source,
  F2+ sink) and ran 8 channels DOUBLE-COUNTED (found 2026-08-05).
  Existence claims walk the BUILT network (build_network output);
  channel uniqueness — one (reactants, products) multiset, one row —
  is now a check_network hard error, and CHEN_SUPERSEDED records which
  generated rows the first-hand set replaces per mode.
- **Never compare rate tiers outside their tabulated overlap.** A clamped
  extrapolation against a table characterises the clamp, not the physics —
  the naive full-range KINETIC/LOKI ratios hit x1e7 on excitation buckets
  whose true-overlap agreement is x0.8–1.1 (found 2026-08-04;
  `kinetic_tier.compare_tiers` clips to the overlap).
- **Chen Eq. (18) `eta_out` is a coupling factor, NOT an efficiency** — it
  exceeds 1 near gl ~ 1.5 at R_OC = 0.08 (measured 1.31). The ledgered
  output power comes from `laser.loss_rates`; Eq. (18) is a reporting
  overlay only. Do not "fix" the budget by switching it to Eq. (18).
- **Normalize ledger residuals by total energy MOVEMENT, not drive alone.**
  Before the pulse the drive is ~0 while seed species rearrange picojoules;
  dividing by drive reports closure noise as failure (the first closure
  "2.6 residual" was exactly this, 2026-08-04).
- **npz bucket keys must round-trip.** `el:Ar` is stored as `el_Ar`; the
  reader maps it back. Breaking the round-trip silently kills the elastic
  drag term when running on the KINETIC tier.
- **Source scans in tests must strip strings/comments** (tokenize) and use
  the AST for imports — a docstring that NAMES a banned flag is
  documentation, not the defect.
- **Never drive the plasma with a prescribed FIELD window.** A fixed E/N
  held across a conducting plasma is an infinite generator: the first
  microwave-reheat cut (8 Td through the dense afterglow) avalanched the
  cell to ionization fraction 1.0 and deposited 7.2 MJ — the ledger
  closed at 1e-16 because the bookkeeping was right and the model
  assumption was wrong (found 2026-08-14). Field-class drives prescribe
  POWER and let the field settle on the swarm balance
  (engine._mw_field_td), or evolve the field as state (circuit.py).
- **Never transplant a paper's validation-experiment components onto its
  base-case conditions.** The 0.1.5..0.1.24 circuit shipped Chen's Fig. 5
  VALIDATION circuit (10.4 nF / 8.1 nF / 30 nH, the "5 -> 1 uH switch" —
  a 2 atm, 20 mm, -13 kV, mirror-less rig) under the Table II base case
  (-18 kV, 12 mm, 3.4 atm), and the graft propagated: the "1.68 J
  stored" energy target contradicted the paper's own stated
  E_in = 1.62 J, and the anchor verdict blamed the network for what the
  circuit was doing (found 2026-08-25, comparison audit; first-order in
  the Chen anchor miss by controlled swap). Constants come from the
  SECTION THAT DEFINES THE CASE BEING SCORED, and E_in recomputed from
  them must match the paper's stated value before the case is called
  transcribed.
- **A verbatim lane does not exempt a physically impossible printed value
  from flagging.** chen2025.py flagged R19/R26/R36 but carried R13/R14
  (cm3-misprint, x1e4 above the electron gas-kinetic ceiling) and
  R23/R24 (R27-paste; x1e6 slip, 0.4 ps Ne* lifetime) as "SI AS
  PRINTED" — and the engine ran them as real rates while the anchor
  verdict hunted elsewhere (found 2026-08-25; the printed table cannot
  even integrate under the paper's own circuit). Every CHEN-tier row
  gets the same ceiling/impossibility screen the LIT rows get; printed
  impossibilities are flagged in the module docstring and corrected via
  CHEN_ERRATA, with the as-printed lane kept as the record.
- **Never list a LOKI table directory twice.** The cross-directory overlap
  guard compares dir NAMES, so the same set passed twice sailed through and
  SUMMED its rates into every bucket a second time — measured x2.0000
  exactly (0.1.39, latent since the guard landed in 0.1.25, reachable once
  rates_for_mode began concatenating run-name lists). A duplicate-name
  check now refuses it.
- **Never "improve" a cross section, and never pick a database by which
  one lands nearer an anchor.** Choosing an evaluation is a PROVENANCE
  decision, recorded with its citation and its measured spread (0.1.38:
  F2 SIGLO vs TRINITI is x1.6-2.0 across channels, Ne Biagi vs BSR
  x1.1-2.1) — it is never a fit. The published spread is an error bar the
  engine carries, not a lever it turns.
- **The flat+taper sub-floor rule is monotone DOWN, and two buckets
  physically go UP there.** Measured 0.1.38 in all three fills: F2
  dissociative attachment (F-) and F2 vibrational excitation (vib) peak
  BELOW the swarm tables' floor, so the extrapolation runs the wrong way
  and the engine is low by x142-x3670 at the operating <eps> of 0.3 eV.
  Every other bucket is threshold-limited and genuinely falls. Do not
  "fix" this by flattening the taper (it would break the threshold
  buckets); the cure is measured data below 0.5 Td.
- **Never promote the MAXWELL tier to a default, blend it with LOKI by a
  weight, or run it on the circuit lane.** It is a verification tier
  (0.1.37): it scores the EEDF-shape assumption side by side. Its own
  verdict is that the swarm shape is the better beam-mode closure —
  a Maxwellian family cannot hold Te (NRL benchmark Te x0.20 vs LOKI
  x1.15) because its tail drains energy through the thresholds
  1e3-1e5x faster than the inelastic-wall swarm EEDF does.
- **Kannari 1985 Table I prints UNREDUCED quench constants.** They are
  Rokni/Jacob/Mangano's k*tau / tau_B (6.5 ns KrF*, 4 ns ArF*) — identical
  to the Avco tables (DTIC ADA047221 pp.15, 27) — and his footnote-a
  0.433/0.418 lives in his MODEL. The 0.1.x comments called them "already
  reduced"; the 2026-09-04 audit proposed a further x0.297. Neither: the
  direct decay measurements (Xu 1993, Gadomski 1992) confirm the
  unreduced scale where both exist (0.1.40). Do not reduce them, and do
  not describe them as reduced.
- **Never adopt a source's 300 K three-body excimer value without the
  temperature form the source carried.** Kannari's footnote b and Lee's
  footnote b both say "300 K"; the engine dropped the Tg dependence for
  0.1.4–0.1.39 and lost the measured pressure scaling of the gain
  (p-exponent +0.03 vs +0.45 until 0.1.40).
- **A solver change is not verified by a green suite.** PIC-copy [0.1.39]:
  a fixed-step Jacobian passed every correctness test at 7.4x the
  integrator work; gate solver changes on the integrator's own counters
  (steps, njev, nfev) against stock, back to back in one process.
- Inherited and binding: the root CHANGELOG's Known-bad list (no
  calibration to any anchor; predictive-only; the two-interpreter rule;
  loose-tolerance comparisons; `ne` is not `y[0]`), PIC_MCC_3D's list
  (cross-section ceilings, no parent-counting path resolution, differential
  W), and the coupling rule: **PIC_MCC_3D/phantom11/spectre are reachable
  as DATA only — `MODULES/` never imports them** (gated).

## [0.1.40] — 2026-09-18 — the fork merged; excimer survival restored at the source (Shui's temperature dependence, Lee/Böwering Ar2F*+Ar), and the measured gain anchors in-tree

Owner ask (2026-09-18): "I believe that my 0D SPECK should be sufficient for
Excimer simulation approximations, lets work on improving it so that the
Excimer kinetics no longer runs low." Rule applied throughout: anchors
SCORE, they never steer — every row changed below is changed because a
first-hand source shows the engine's value is not the source's, and the
anchor numbers are reported after the fact, including the ones that did
not move and the dependencies of the ones that did.

### Merged — the ~/coding/PIC fork (the two copies diverged after 0.1.36)

The PIC bundle's REBUILD (cut from this tree at 0.1.36; its committed HEAD
48d3c36 is byte-identical to this tree's 0.1.36 in every engine file) went
on to its own [0.1.38] (Kr-partner rows) and [0.1.39] (the num_jac cure)
while this tree took [0.1.37]–[0.1.39] (MAXWELL tier, sigma atlas,
sub-floor tables) — two histories with COLLIDING numbers. Both PIC changes
are now here: network.py / tests/test_network.py are byte-identical to
0.1.36 + `PIC/Neon/patches/rebuild-0.1.38-kr-partner-rows.patch`;
engine.py is a three-way merge (base PIC HEAD, ours this tree's 0.1.39,
theirs PIC's working copy) with ZERO conflicts (this tree touched params
plumbing only, PIC the Jacobian and the solve call); tests/test_engine.py
carries PIC's Kr-row gates + the three Jacobian gates. Two stale comments
in the imported `Device.jac` were corrected (the docstring said "FIXED
relative step" — that was the rejected first draft; the cost table marked
the 1e0 cap as shipped, the code ships 1e12; note also that PIC's entry
says the 0-D tolerance floor is ON for 0-D, while the code it shipped has
`_jac_floor = False` by default — opt in per device). Their two changelog
entries are carried verbatim at the end of this entry. The PIC copy does
NOT carry this tree's 0.1.37–0.1.40 (Cooper's call).

Measured on the merge alone: suite 463 passed (458 + the 5 imported gates,
14:07); krf_beam_0d E_L 6.778e-07 -> 6.844e-07 J (+1.0 %, the Kr rows);
arf_chen_circuit 2.040 mJ unchanged; ne_optimized 7.773 -> 7.818 J
(+0.6 %, Kr rows at 5 % Kr — inside the patch's "< 1 %" statement).

### Fixed — the excimer three-body rows had lost their temperature dependence

Every excimer three-body association row is a 300 K value whose source
carried a steep negative Tg dependence that was dropped when the row was
adopted. Read first-hand 2026-09-18 (page images):
- Kannari, Obara & Fujioka 1985 Table I: KrF*+2Ar 7.0e-32, KrF*+Kr+Ar
  6.5e-31, KrF*+2Kr 6.7e-31 (all "62,63^a,b"), ArF*+2Ar 4.0e-31;
  footnote b: "Expressed for 300 K gas temperature. Temperature dependence
  is available from the calculation method in Refs. 49–51" (= V. H. Shui,
  APL 31, 50 (1977); APL 34, 203 (1979); Shui & Duzy, APL 36, 135
  (1980)); Sec. II.D: the rates were "treated here as a temperature
  dependent rate, which was based on the theoretical calculation by Shui".
- Shui's KrF calculation itself, in the Avco contract report that carries
  it (Mangano et al., "One Meter KrF Laser System", DTIC ADA069053, Sec.
  V.B p.73; Fig. 22 read off the page, 5.1e-31 -> 5.3e-32 over 300–600 K):
  "a change in temperature from 300°K to 600°K lowers kr(Ar) by a factor
  of 8 and kr(Kr) by a factor of 10" -> (Tg/300)^-3.00 (R = Ar) and
  ^-3.32 (R = Kr), R = the atom that binds.
- MEASURED: Klimek & Hsia, JAP 52, 5361 (1981) (Crossref abstract): F2/Kr/
  Ar, 315–420 K, 1–3 amagat, KrF*+R+M with R, M = Ar, Kr: "The decrease in
  the three-body rates with increasing temperature deduced from the
  experimental data, agrees to within 15% of that calculated by using
  classical phase-space theory."
- ArF: Lee, Matsui, Kannari & Obara, IEEE TED 36, 2053 (1989) Table I
  (ArF*+2Ar 5.0e-31, footnote b "At a gas temperature of 300 K"; text:
  quench "reduced with increasing excitation rate due to the gas
  temperature rise") and NRL Orestes 2022 (petrova2022 Table A5, red), the
  one printed ArF form: 1.62e-32 + 2.245e-21 Tg^-3.9 cm6/s [k28, Shui-1979,
  Chen-1978]. (Aside: NRL's 2021 black row prints Lee's value as 5.0e-32;
  Lee's page prints 5.0e-31.)

Restored on each row's OWN 300 K value: KAA (R = Ar) x(Tg/300)^-3.00;
K34, K34Ar, K68Ne (R = Kr) x(Tg/300)^-3.32 (K34Ar: Hsia/Jacob/Duzy take
k_ArKr ~ k_2Kr; K68Ne: Ne third body, by Shui's R rule — no Ne study in
the tree, stated EST); ARF_AA x[f + (1 - f)(Tg/300)^-3.9], f = 0.0320
(NRL's form normalised at 300 K; `Reaction.tg_floor`, new). R36
(ArF*+Ar+Ne, generated CHEN lane, flux ~0 in these fills) left verbatim.

### Established — the 300 K magnitudes STAY; Kannari's table prints UNREDUCED values

Kannari's footnote a says his quench constants came from k*tau products
and "were modified using the effective radiative lifetime" (0.433 KrF*,
0.418 ArF*). The engine's comments (and the 2026-09-04 audit) read the
printed numbers as already modified. They are not: the Avco manuscripts
of Rokni, Jacob & Mangano (PRA 16, 2216 (1977); DTIC ADA047221, Table 3
p.27 KrF*, Table 2 p.15 ArF*, page images) print k*tau and k at tau_R =
6.5 ns (KrF*) / 4 ns (ArF*): F2 5e-18 -> 7.8e-10, 2Kr 4.4e-39 -> 6.7e-31,
Kr+Ar 4.2e-39 -> 6.5e-31, 2Ar 4.6e-40 -> 7e-32; ArF* F2 1.9e-9, Kr 1.6e-9,
Ar 9e-12, 2Ar 1.6e-39 -> 4e-31 — every one identical to Kannari's print.
The 0.433/0.418 lives in Kannari's MODEL. Which scale is physical? Where a
DIRECT decay measurement exists the unreduced one is what it confirms:
2Kr 6.7 vs Gadomski/Xu 6.0 +- 2.0, Kr+Ar 6.5 vs Xu 5.5 +- 2.0 (e-31);
k*tau/15 ns (Kannari's model) misses x2.1 and k*tau/25 ns (Xu's own
footnote-d adjustment, ~1.7e-31, and the audit's proposed KAA x0.297)
misses x3.3. So no magnitude changes: KAA 7.0e-32 (no direct measurement;
the cross-check carries it), K34Ar 6.5e-31, and the lumped ArF* rows
(tau_B-based, the consistent convention for a one-state ArF* radiating at
tau_B — Lee 1989's footnote c states the same). The audit's KAA x0.297 is
NOT applied. The ARF_AA comment's "printed value is footnote-a reduced by
0.418" was wrong and is corrected.

### Fixed — A2F_AR measured, was Kannari's "estimated"

Ar2F*+Ar -> 3Ar+F: 1.0e-10 cm3/s (Kannari 1985, "estimated") -> 2.2e-14:
Lee, Matsui, Kannari & Obara 1989 Table I [58] (page image) — Kannari's
own group; Böwering, Sauerbrey & Langhoff, JCP 76, 3524 (1982) measured
(2.2 +- 0.8)e-14 (second-hand, Gu et al., Acta Phys. Sin. 35, 228 (1986)
Table 1, who measured (3.07 +- 0.74)e-14 themselves); NRL Orestes carries
2.2e-14. Products kept as Lee prints them (NRL 2022's ArF*+2Ar re-route is
labelled "estimated products" — not adopted).

DEPENDENCY, measured and stated plainly: in ArF fills the change scores
flat (as the audit found). In KrF fills it matters because an Ar2F* that no
longer flashes through lives long enough for Kannari's ESTIMATED
Ar2F*+Kr -> KrF*+2Ar (A2F_KR2) to return it, fed by ArKrF*+Ar -> Ar2F*+Kr
(AKF_AR, 2.0e-11, Kannari/J&H via Mangano 1977 [43]). That feed moves a
Kr+-cored triatomic to an Ar2+-cored one, which the ion-core energies
(IP(Ar) - IP(Kr) = 1.76 eV; Ar2+ bound ~1.2 eV vs ArKr+ ~0.1 eV) suggest
is ~0.6 eV uphill — an argument, not a source, so the row is kept as
printed. With either row zeroed the KrF gain from this correction
vanishes: full grid, AKF_AR = 0 -> Klimek x0.53 (p-exp @110
+0.39, g0/alpha x0.88 / x0.68); A2F_KR2 = 0 -> x0.52 (+0.38, x0.84 /
x0.65). The Tg restoration ALONE therefore carries Klimek x0.42 -> ~x0.52
and the high-current pressure exponent +0.03 -> ~+0.38; the further
x0.52 -> x0.64 rides on the two weak rows.

### Added

- `tools/excimer_anchors.py`: the 2026-09-04 audit's gain anchors as an
  engine tool — Electra ArF small-signal gain at PRINTED deposition
  (matched exactly: the drive is linear in j and independent of the
  plasma, so one analytic rescale, verified to 1e-9) and the Klimek KrF
  grid (Kannari Tables II/III, F2 excluded from absorption per his
  p.4316–17 convention; pump on Kannari's printed 0.75 MW/cm3 atm at
  187.5 A/cm2, p.4314 — a DIFFERENT device, stated), plus the
  pump-independent discriminators (g0/alpha per pressure, the pressure
  exponent). Refuses to score on the W-table fallback. Reproduces the
  audit's pre-0.1.40 numbers (Electra x0.68, Klimek x0.42). ~2 min,
  --quick 30 s.
- `Reaction.tg_floor` + check_reaction rule (needs tg_pow, 0 < f < 1);
  f = 0 takes the untouched `tg_fac ** tg_pow` branch, so every existing
  tg_pow row is bit-identical.
- `tests/test_excimer_anchors.py` (6): analytic deposition == the ledger's
  drive_in (rel 1e-6); linearity + exact matching; the W-table refusal;
  printed-target pins; the tg_floor form at 300/400 K and f = 0; misuse
  refused. No anchor ratio is a test target.

### Measured — the anchors, full grid (22 points), before -> after

| | pre-0.1.40 | 0.1.40 | measured |
|---|---|---|---|
| Electra ArF gain, peak / flat top | x0.68 / x0.60 | x0.75 / x0.67 | 1 |
| Electra per point, peak (0.82 / 1.09 / 1.36 atm; 10-cm) | 0.59 / 0.64 / 0.87; 0.63 | 0.63 / 0.72 / 1.01; 0.67 | 1 |
| Electra 0.82 -> 1.36 atm | x1.04 | x1.13 | x0.70 |
| Klimek KrF gain (18) | x0.42 (0.25–0.67) | x0.64 (0.48–0.79) | 1 |
| Klimek absorption (18) | x0.66 | x0.68 | 1 |
| g0/alpha, 1 amagat | x0.78 | x0.96 | 1 |
| g0/alpha, 3 amagat | x0.51 | x0.86 | 1 |
| p-exponent @110 A/cm2 | +0.03 | +0.57 | +0.45 |
| p-exponent @20 A/cm2 | -0.29 | +0.01 | +0.13 |

Reading. KrF: the pump-independent g0/alpha is now within 4–14 % of
measured at both pressures, while gain and absorption sit low by the SAME
factor (x0.64 / x0.68) — the signature of the anchor's ASSUMED pump
(Kannari's device calibration applied to Klimek's device), not of the
kinetics; the high-current pressure exponent overshoots (+0.57 vs +0.45)
and the low-current one still undershoots (+0.01 vs +0.13). ArF: the
deposition is PRINTED, so the level is firm — x0.63 at the most uniform
point (0.82 atm) is a real, unexplained deficit, while the 1.36-atm point
now reaches x1.01 and the pressure trend runs the wrong way (x1.13 vs
x0.70). Decomposition (quick set: Electra 0.82/1.36 atm + Klimek 110/20
A/cm2 x 1/3 amagat, 0.2 % F2, earlier -3.9 KrF form): baseline Klimek
x0.38 / Electra x0.72; Tg form only x0.48 / x0.80; A2F_AR only x0.49 /
x0.72; both x0.59 / x0.80.

### Measured — the cost: in-loop ASE on MULTI-CELL grids (0-D untouched)

Raising the excimer population raises the gain, and the in-loop ASE term is
exponential in it (engine._ase_rate: phi ~ sum g0 (e^{(g0-alpha)L} - 1) /
(g0-alpha), net gain x chord). On the 24-cell gradient test case
(6x1x4, 3 atm, 0.75 m) that term, not the chemistry, carries the whole
extra solver cost — integrator counters (load-independent), 40 ns:

| | steps | Jacobians |
|---|---|---|
| pre-0.1.40, in-loop ASE ON | 325 | 7 |
| 0.1.40, ASE ON | 503 | 66 |
| pre-0.1.40, ASE OFF | 322 | 5 |
| 0.1.40, ASE OFF | 321 | 5 |

With ASE off the two are identical; 0-D is untouched either way (1 cell,
3 atm, 100 ns: 783 -> 792 steps, 19 -> 20 Jacobians). To 100 ns the same
24-cell case costs x10 steps / x46 Jacobians with the Tg form alone and
x17 / x79 with the A2F_AR change alone (1113 steps / 115 Jacobians before).
The suite's own 100-ns gradient test
(tests/test_engine.py::test_3d_gradient_case_runs_and_closes) is this case
and is affected the same way: each HALF of 0.1.40 integrates it
successfully at 100 ns (success=True, counters above), but the run now
costs of order an hour, so it is DESELECTED from the suite loop and run on
its own. Pre-0.1.40 it is ~5 min. Options for the owner, none taken here:
leave it (and deselect when iterating), shorten its horizon, or run it with
`ase_inloop=False` and keep a shorter ASE-on case beside it.

Observed in passing, benign but latent: in MULTI-CELL mode the 0.1.39
dead-column protection cannot see a structurally inert species, because
`jac_sparsity` is a dense per-cell block, so only the ledger columns are
detected. F2(a3)/F2(A1) (no producers since the 0.1.14 repulsive reroute)
therefore ride their num_jac step factor up to the 1e12 cap in every cell.
Pinning them was measured to change nothing (503 steps / 66 Jacobians
either way) — nothing reads them — so it is recorded, not fixed.

### Measured — decision inputs, NOT adopted (owner's call)

- The MEASURED B/C manifold. Gadomski, Xu, Setser & Nelson, CPL 189, 153
  (1992) fit K_eq = 5.0 +- 1.0 at 300 K ("~4–5", KrF(C) lower, gap ~300
  cm-1) and Xu 1993 confirms tau_eff = 25 +- 2 ns; the engine's K_eq =
  5.0/3.41 = 1.47 (Kannari's 80 cm-1) gives 15.2 ns. Adopting K_eq 3.86
  (tau_eff 25 ns with tau_B 7.0 / tau_C 75 ns) LOWERS gain: full grid Klimek x0.64 -> x0.47, g0/alpha x0.71 /
x0.69, p-exp @110 +0.67.
  Better-sourced than K = 1.47, recorded not shipped: it is a structural
  change to the lasing fraction and the owner decides it.
- The pending flags on top of 0.1.40 (full grid): petrov17_f2 -> Electra x0.79 (flat top x0.74),
Klimek x0.65, absorption x0.73, p-exp @110 +0.43, g0/alpha x0.94 (1, 2
amagat) / x0.81 (3); loki_subfloor -> Electra x0.75, Klimek x0.64,
p-exp @110 +0.42, g0/alpha x0.99 / x0.99 / x0.86. Both still
  OFF (their own entries give the reasons).
- The formation-side W: `spencer_fano_yields` drives with the INTEGRAL W at
  a finite cascade E0 (Ar 29.3 eV), while the scorecard validates the
  DIFFERENTIAL W (Ar +5 % vs ICRU 26.4). For a 500–600 keV beam the
  asymptote is the relevant one, so pair production per joule is likely
  ~5 % low. Not changed (deposition layer, beam-yield path); recorded.

### What moves (0-D benchmark cells, no microwave, pre -> post 0.1.40)

krf_hyper_mw 0.494 -> 0.934 J (eta 1.56 -> 2.96 %); ar_optimized 5.095 ->
5.523 J (5.00 -> 5.42 %); ne_optimized 7.818 -> 7.938 J (6.23 -> 6.33 %);
krf_beam_0d x3.0 and krf_ne_0d x1.04 (both sub-threshold micro-joule
outputs — the knife-edge). The recorded SPECKLE3D baseline (ne_optimized
7.773 J) is therefore +2.1 % stale (merge +0.6 %, 0.1.40 +1.5 %); the
reheat marginals computed against it were NOT recomputed.

### Still open

- Klimek's own pump and mixture (IEEE JQE QE-17, 1847 (1981), DOI
  10.1109/JQE.1981.1071324) — not reachable; they set the KrF level.
- ArF: the level at the most uniform Electra point (0.82 atm) and the
  pressure trend (see the table): the audit's deposition-geometry
  explanation for the trend is unproven; the level is unexplained
  (candidates: the lumped ArF* manifold, Ar* formation per joule).
- AKF_AR energetics and the estimated A2F_KR2 (above) — the KrF level
  from A2F_AR rests on them.
- Shui & Duzy 1980's 300 K theory values (DTIC ADA079033 Table IV:
  KrF*+2Ar 2e-31, KrF*+Ar+Kr -> ArKrF* 1e-30) are theory "within a factor
  of three" of experiment; not used for magnitudes.

### Verified

suite 468 passed + 1 deselected in 8:40 under /usr/bin/python3 (the deselected case is the 100-ns 3-D gradient test — see the ASE cost section: it is run on its own, not in the suite loop); Homebrew smoke (-W error::DeprecationWarning) on krf_beam_0d
and arf_chen_circuit; the shipped network reproduces the exploration
variant it was measured as.

### Imported: PIC-copy [0.1.39] — 2026-09-10 — the stiff blocker was scipy's Jacobian bookkeeping, not the physics


Owner ask: "rebuild the entire numerical Jacobian ... use stiff-source
integration exact per step." The 0.1.35 blocker — `numjac overflow` — had
stood since 2026-08-28 and had killed the 0.8 % F2 lever point and the
multi-cell 3-D runs.

#### Fixed — the runaway, and where it came from

- **Root cause, from a traceback rather than a guess.** BDF built its own
  Jacobian through `scipy.integrate._ivp.common.num_jac`, which keeps a
  per-column step `factor` and multiplies it by ten whenever a column's
  finite difference comes back "too small". **Every ledger column of this
  system is structurally dead** — nothing anywhere reads an accumulated
  ledger value, by design, because the ledger is a pure quadrature of the
  physics. Such a column is "too small" forever, so its factor climbs until
  `new_factor * y_scale[ind]` overflows against a ~1e25 density
  (`_sparse_num_jac`, line 409). Solver bookkeeping meeting a deliberate
  modelling choice; no amount of physics work would have moved it. Supplying
  the sparsity pattern alone does NOT cure it — measured, the pattern went in
  first and the same overflow returned — because an empty column still gets
  its factor inflated.
- **The fix is surgical: `Device.jac` pins the dead columns and leaves every
  live one alone.** A column with no structural nonzeros is the only one that
  can never produce a difference, so it is the only one that climbs on every
  call; live columns self-correct as soon as their step is large enough, and
  their adaptation is what makes the algorithm work. A cap of 1 on the factor
  sits on top as belt-and-braces: a step larger than the state it perturbs is
  not a derivative, and it is ~293 decades below the overflow.
  `atol_vector()` is now one definition shared by the error test and the
  Jacobian, since the two must agree.
- **Nothing about the equations changed.** A Newton matrix sets convergence
  quality, not accuracy; the step controller still governs the error test.

#### The 0-D tolerance floor, and why it is not general

A near-threshold 0-D cavity run needs more: forcing every step to
`h >= atol` takes the 0.8 % F2 point from not finishing in ten minutes to
**8 s, same answer**. That floor is switched on for 0-D only, and the switch
is a measured trade, not a preference. On the 24-cell gradient case the same
floor costs **12x (121 s against 9.8 s)**, because there most species in most
cells are cold and a tolerance-scale probe is a large relative perturbation
that degrades the matrix; in 0-D there is one cell whose species are all
live, so the same probe is informative. Narrower forms were tried and do not
work — flooring only states already below their atol, and flooring only the
photon column, both leave the case unfinished.

#### Measured

| case | before | after |
|---|---|---|
| 0.8 % F2 lever point (0-D, cavity) | never finished | **8 s, 388.1 J** |
| Electra baseline (0-D, cavity) | 480.14 J | **480.05 J, 11 s** |
| 24-cell gradient vs stock scipy | 176 steps / 6 jac | **176 steps / 6 jac** |
| 0-D Electra vs stock scipy | 1162 steps / 37 jac / 2838 rhs | **1164 / 35 / 2814** |

Two warnings about how those numbers were obtained, both learned the hard
way here.

**Wall time on this machine is unusable.** Other sessions held two cores
throughout; the same suite read 5 h and then 1 h 15 for a change that costs
nothing, and one background run reported 266,079 s because it sat idle
across a multi-day gap. Check `ps` before believing any timing.

**But counters are only honest over the window you measure.** The 40 ns
multi-cell row above showed exact parity while the SAME case run to 100 ns
was 12x slower (4259 s against 342 s), because the blow-up lived entirely
in the afterglow a 40 ns window never reaches. That is what forced the cap
from 1 to 1e12. A cheap proxy for a slow test is not a regression check —
time the actual test.

- The F2 lever scan at 40 % Kr is now closed with its upper arm measured
  rather than assumed: 0.1/0.2/0.3/0.5/0.8 % gives 174/406/480/467/388 J, so
  **0.3 % is the optimum** and the fall-off beyond it is real.
- Answers are unchanged: 0.018 % on the baseline, five identical figures on
  peak g0. Every stiff run above was made with `RuntimeWarning` promoted to
  an exception so an overflow could not pass silently.

#### A wrong turn, recorded because it nearly shipped

The first version of this fix replaced scipy's adaptive step with a fixed
one. It removed the overflow and **passed the full suite**, and it was
wrong: on the 24-cell case it needed **296 steps and 50 Jacobians against
stock's 176 and 6**. Killing the adaptation to kill the runaway costs
roughly 7x in work; bounding it costs nothing. A suite that only asserts
correctness will not catch a solver change of this kind, so a solver change
must be gated on the integrator's own counters.

#### Verified

- Full suite under `/usr/bin/python3`: **208 passed in 10:24**, against the
  0.1.38 baseline of 205 in 11:41 on the same machine under the same load.
  Three of those are this entry's new gates, so the fix is slightly FASTER
  than the code it replaces, not merely no slower.

#### Added — gates

- `test_ledger_columns_are_structurally_dead` — pins the exact condition that
  broke `num_jac`, and asserts the delivered Jacobian is zero there and
  finite everywhere. If it ever fails because a rate learned to read a
  ledger, the runaway is cured but the ledger has stopped being a passive
  observer, which is worse.
- `test_jacobian_matches_finite_differences` — the Newton matrix must be the
  real derivative where the derivative is well defined.
- `test_jacobian_costs_no_more_than_stock_num_jac` — the gate this entry
  needed and did not have. Supplying our Jacobian must cost no more steps,
  Jacobians or RHS calls than letting scipy build its own, measured back to
  back in one process. Counters, not seconds: the first version of this fix
  passed every correctness test while doing 8x the solver work, and wall
  time on a loaded machine said nothing useful.

#### Recorded, not done

- **An exponential / quasi-steady-state integrator.** The owner ask named
  "exact per step" integration, the `n/(1 + n k dt)` and `1 - exp(-dt/tau)`
  idioms. Those are exact for a scalar source; the network equivalent is a
  production/loss split integrated as `n exp(-L dt) + (P/L)(1 - exp(-L dt))`,
  which is a DIFFERENT SOLVER — it would move every recorded number and need
  the whole anchor set re-run. Not needed now that the stiff points complete,
  and it should not be bought with a re-validation of the entire scorecard.
  If wanted, it belongs behind a flag with the anchors scored on both.
- A hand-derived analytic Jacobian for the chemistry block. The row structure
  makes it mechanically derivable and it would be faster than differencing,
  but it would also silently omit whatever a hand derivation forgets.

### Imported: PIC-copy [0.1.38] — 2026-09-04 — the krypton-partner rows the 0.1.4 split left out (found by a pure-Kr fill)


Owner ask (PIC/Neon replication thread): "replicate the NRL KrF and ArF,
then test with neon" → "what raises Electra's efficiency" → "is the
krypton optimum real or is the engine wrong?" → "fix the missing rows and
track down the failure." Record: `PIC/Neon/NRL_REPLICATION.md` §6d–6e, §7.

**Version note.** This bundle jumps 0.1.36 → 0.1.38 on purpose. [0.1.37]
(the MAXWELL verification tier, 2026-09-03) was written in the Agent tree
after this bundle was cut and is not present here; taking 0.1.38 keeps the
two histories mergeable rather than colliding two different 0.1.37s. The
patch in `PIC/Neon/patches/` inserts this entry above the Agent tree's
0.1.37, and touches nothing the MAXWELL tier touched.

#### Fixed — a real gap, 32 versions old

- **The KrF network fed the lasing state KrF*(B) from KrF*(C) and from the
  hot manifold KrF*(Bv) through argon and neon partners only.** [0.1.4]
  wrote B<->C mixing (5e-10) and Bv->B relaxation (4e-11) from Kannari's
  abstract, which states them *with argon*, and added the Ne analogs the
  same day; no krypton partner was ever written, although BC_Ne's own
  comment (2026-08-04) quotes Xu 1993 p.2596 that k_BC for He/Ne/Ar "may
  be much smaller than for Kr". Consequence: a buffer-free Kr/F2 fill had
  NO channel into B — KrF*(C) and (Bv) reached 1.6e14 and 1.4e14 cm^-3
  while KrF*(B) sat at its 1e6 floor, g0 < 0, zero output — where Mandl,
  Klimek & Salesky (CLEO 1985) measured 0.82 of the 10 %-Kr optimum and
  Suda 1984 / Razhev 2016 ran buffer-free KrF lasers; and a 40 %-Kr
  Electra-class fill was starved of relaxation (402 J at fixed beam vs
  480 J once the rows exist, Spencer-Fano route; 327 vs 389 J on the
  W-table route the bundle was silently running, see below). Why it survived: every anchor, preset and
  scorecard cell is a 5 %-Kr fill where Ar/Ne do the mixing; the 0.1.24
  ledger's Kr scan stopped at 10 % ("table caveat"); no test touched the
  BC_/CB_/VREL_ rows or ran a Kr-rich fill; `kannari_ps_relaxation` DID
  follow the buffer-partner policy (KPS_Kr exists) but the mixing block
  written the same day did not; EST bracket sweeps price rows that exist
  and cannot find one that is absent. Exposed by an external scan the
  engine had never been asked for (Mandl's 99.6 %-Kr point), then a
  species diagnosis.

#### Added — five Kr-partner rows at the class values (EST, primaries named)

- `krf_bc_mixing`: **BC_Kr** 5.0e-10 (direction: Xu 1993 says the Kr
  partner is LARGER; primary Gadomski/Xu/Setser/Nelson CPL 189,153 (1992)),
  **CB_Kr** 3.41e-10 (detailed balance 0.682), **VREL_Kr** 4e-11 (primary
  Kvaran/Shaw/Simons Appl. Phys. B 46,95 (1988)). `kannari_ps_relaxation`:
  **APS_Kr** 1.0e-10 (Ar** p->s on Kr). `ne_kr_couplings`: **AK8_Kr**
  Kr*+Kr+Kr -> Kr2*+Kr 2.7e-32 (Kr third body; K15 covered Kr2+, nothing
  covered Kr2*). A partner audit (`PIC/Neon/scripts/partner_audit.py`) of
  the 22 catalytic-partner families found only one (Kr** p->s) complete
  before this entry; the remaining Ar-only third-body twins (AKD2/AKI2
  with Kr as M) are recorded, not added — they move the Ar-fill anchors
  and need their own measurement.
- Gates: `test_krf_buffer_partner_families_complete` (the five
  feeding/relaxing families must carry Ar, Ne and Kr with detailed balance
  per partner; the pre-0.1.38 network fails it by construction — teeth) and
  `test_pure_krypton_fill_lases` (Kr/F2 0.997/0.003 at 1 atm must reach
  g0 > 0.5 /m and populate KrF*(B) > 1e19 m^-3). EST census pin 19 -> 24,
  justification in the test.

#### Fixed, bundle only — the PIC bundle carried no LXCat tree

- `cross_sections.lxcat_roots` finds the data through `REBUILD_LXCAT`,
  `MODULES/lxcat/`, or the repo markers (`LoKI-MC/Code/Input`) that exist
  only in the Agent tree. In the PIC bundle none resolved, so
  `Device.__init__` caught the FileNotFoundError and took the W-table
  fallback (`w_table_yields`) for EVERY beam-driven run since the bundle was
  cut (2026-09-01), recording it only in `beam_yields.cite`
  "[SF route unavailable: ...]"; five `test_beam` tests and
  `test_krf_beam_0d_physics_band` (Te peak 2.96 eV against the < 2.0 eV
  pin) failed in the bundle for the same reason, before this entry. The
  five files now ship in `MODULES/lxcat/` (Ar/Kr/F2/He byte-identical to
  `LoKI-MC/Code/Input`, Ne = the Biagi set the root CHANGELOG mandates), so
  the bundle runs the Spencer-Fano route (emergent W 28.5 eV for the
  Electra mixture, 56 processes). The Agent copy is unaffected (it resolves
  the repo tree); the numbers below are on the Spencer-Fano route, with the
  W-table-route values the Neon record quoted up to 2026-09-04 in brackets.
  (The PIC MAIN checkout acquired an untracked `MODULES/lxcat` SYMLINK to the
  Agent tree's `LoKI-MC/Code/Input` at 13:25 on 2026-09-04, from outside this
  job; it resolves the nested paths and so runs the SF route, but it is not
  in git and breaks the moment the bundle moves without the Agent tree. The
  real files are committed on `worktree-neon-ramsauer`. See
  `PIC/Neon/patches/README.md` §3.)
  The routes disagree about MIXTURES, not only about W: the cascade over the
  mixture's cross sections puts 4.3x more quanta on 5 % Kr in neon than in
  argon (the cascade window between Kr's and Ne's thresholds; W-table splits
  by stopping power), which reverses the PIC/Neon record's neon-vs-argon
  verdict at <= 20 % Kr (its Section 7). Nothing in the engine changed there;
  the bundle had simply never run the route.

#### Measured — fixed beam, Electra cell (30x30x100 cm, 500 kV, 20/100/20 ns,
1.36 atm, R_OC 0.4 / L_op 0.10), before -> after (nothing tuned;
`PIC/Neon/scripts/sf_route_rerun.py`, "before" = the five rows removed)

- Electra mixture Ar/Kr/F2 0.597/0.40/0.003: **402 -> 480 J (5.32 -> 6.35 %)**
  [W-table route: 327 -> 389 J]; 60 % Kr 261 -> 420 J (x1.61); 20 % Kr
  473 -> 500 J; 10 % Kr 459 -> 469 J; 5 % Kr 403 -> 407 J (+1 %: Ar does the
  mixing there, so the 5 %-Kr anchors and the recorded SPECKLE/SPECKLE3D
  baselines stand). The 40-vs-10 % laser-energy ratio 0.88 -> 1.02: the
  "10-15 % Kr optimum, +27 %" of the uncorrected engine collapses to a
  20 %-Kr maximum 4 % above the Electra point, inside the x0.5-x2 bracket
  of the Kr2F* three-body rows (K34/K34Ar flip the sign of the 40-vs-10 %
  comparison on the old engine). Electra's mixture is at, or within
  uncertainty of, the engine's optimum; the intrinsic-efficiency optimum
  stays at 10 % Kr (7.65 %) because deposition per atmosphere rises with
  krypton.
- Mandl 1985 Kr scan (350-ns flat top, deposition matched at 23.3 kJ):
  4 %/10 % ratio 0.75 -> 0.74 (measured 0.83); **99.6 %/10 %: 0.00 -> 0.27
  (measured 0.82)** — the pure fill lases now (1.75 %), the remaining x3 is
  the open high-Kr chemistry: K34 (6.0e-31; at 99.6 % Kr and 1.36 atm it
  alone is 4.7x the radiative rate), Kr2*/Kr2F* branching, Kr2+/Kr2F*
  absorbers, and the 5 %-Kr LOKI tables extrapolated to a Kr fill. Named,
  not fixed here.
- Unchanged where it should be: the F2 optimum at 40 % Kr (0.1/0.2/0.3/
  0.5 %: 174/406/480/467 J; the 0.8 % point is stiff on this route —
  numjac overflow, the 0.1.35 blocker in 0-D, killed after 4.5 min; CURED in
  0.1.39 and now measured at 388 J / 5.15 %; W-table
  route 301 J); pulse 200 ns 978 J, 6.35 -> 7.06 % intrinsic (+11 %); L_op
  0.10 -> 0.05 506 J (+5 %); 2.0 atm at 10 % Kr 527 J (5.85 %).

- S5 bracket pricing of the five rows (`PIC/Neon/scripts/est_bracket_0138.py`;
  x3 and x1/3, dg0_rel at the 5 %-Kr anchor / laser energy in the Electra
  40 %-Kr fill): **VREL_Kr** +2.0/-0.7 % / **+16/-10 %** — the one real lever,
  the Kvaran 1988 read is the EST->LIT move that matters; **AK8_Kr**
  -0.3/+0.1 % / -6.6/+2.5 %; **APS_Kr** +0.2/-0.1 % / +0.7/-0.5 %. **BC_Kr and
  CB_Kr one-sided read -11/+4 % and +12/-4 % (anchor), -38/+13 % and
  +23/-18 % (Electra) — that is the detailed-balance ratio being broken,
  not a rate uncertainty**: scaled as a pair (equilibration rate, ratio
  0.682 kept) they move g0 by +0.3/-0.2 % at 40 % Kr and 0.0 % at 5 % Kr,
  because B/C equilibration is 10-100x the radiative rate at every krypton
  fraction. The S5 sweep (`run.py --est-sweep`) brackets rids one at a
  time and will print the one-sided numbers for these two; read them as
  the pair. The active EST set in the Ar fill is now nine
  (`EST_ACTIVE_AR_FILL` in test_network.py), bound set unchanged (15).

#### Verified

- Full suite under `/usr/bin/python3`: **205 passed** (11:15, PIC bundle with
  `MODULES/lxcat/` in place; 199 before this entry + the two new gates + the
  four `test_beam`/physics-band tests that the missing LXCat tree had been
  failing). Before the LXCat files were bundled the same suite ran
  203 passed / 6 failed. Two pins moved and are justified in the tests:
  the EST census 19 -> 24, and `EST_ACTIVE_AR_FILL` (the S5 active set) from
  four rids to nine.

#### Recorded, not done

- **K34/K34Ar look high by ~x2.** Bracketing the Kr2F* three-body pair at
  x0.5/x1/x2 on this route (`PIC/Neon/scripts/kr_bracket.py`, fixed beam,
  Electra cell) gives the Mandl 99.6 %/10 % efficiency ratio 0.75 / 0.29 /
  0.00 against the measured 0.82, and the 40 %/10 % laser ratio 1.30 / 1.02
  / 0.37. Half the class value reproduces the pure-krypton point; the
  Electra point would read 667 J (8.8 %). That is the measurement to chase
  next (a Kr2F* formation rate in a krypton-rich fill), not a buffer study.
- EST -> LIT for the five rows when Gadomski 1992 / Kvaran 1988 are read
  first-hand (the Kr mixing value is expected ABOVE the class value).
- The Ar-only third-body twins AKD2/AKI2 (Ar*+Kr+M, Ar+ +Kr+M with M=Kr).
- The Electra multi-cell ASE check (2x2x4, 6x3x10) did not complete under
  plain `Device.run` (numjac overflow, the 0.1.35 blocker, cured in 0.1.39
  — the multi-cell case now integrates instead of overflowing); it needs the
  SPECKLE3D 30-ns restart driver. The box-geometry ASE report on the 0-D
  state at the gain peak is in `PIC/Neon/NRL_REPLICATION.md` §6c.

## [0.1.39] — 2026-09-04 — LoKI-MC below the floor: the sub-0.5 Td swarm data, measured

Owner ask: "Run the LoKI-MC at the low Td and report the results" — the one
action [0.1.38] named. Seven single-point boltzmannMC runs (2e4 electrons,
smartActivation gas-temperature term, the shipped cross sections), 28-62 min
each, run in parallel on this machine 2026-09-04.

### Measured — the swarm DOES keep cooling with the field (the 0.1.38 caution was wrong)

| fill | E/N (Td) | <eps> (eV) | Te (eV) | v_d (m/s) |
|------|----------|------------|---------|-----------|
| Ar/Kr/F2 | 0.05 | **0.761** | 0.507 | 10600 |
|          | 0.1  | 1.072 | 0.715 | 8579 |
|          | 0.2  | 1.427 | 0.952 | 6621 |
|          | 0.3  | 1.660 | 1.107 | 5665 |
|          | 0.5 (old floor) | 1.980 | 1.320 | 4840 |
| Ne/Kr/F2 | 0.1 | **0.891** | 0.594 | 11861 |
|          | 0.3 | 1.740 | 1.160 | 7669 |
|          | 0.5 | 2.231 | 1.488 | 6708 |
|          | 1.0 (old floor) | 3.063 | 2.042 | — |

[0.1.38] warned the mean energy might saturate (attachment eating the cold
electrons, the Ramsauer window decoupling heating from the field) and that
the sub-floor band might be unreachable by a swarm. It is not: a decade in
E/N buys a factor 2.6 in <eps>, and 0.05 Td lands at 0.761 eV — inside the
beam-mode operating band. The extension works; the caution is withdrawn.

### Measured — the taper error, now against data instead of an estimate

k_MEASURED / k_ENGINE(flat+taper) at the lowest new point of each fill:

| fill | <eps> | F- | vib | elastic |
|------|-------|-----|-----|---------|
| Ar/Kr/F2 @ 0.05 Td | 0.761 eV | **x25.2** | **x14.9** | x0.80-2.17 |
| Ne/Kr/F2 @ 0.1 Td  | 0.891 eV | **x63.8** | **x44.5** | x0.96-6.68 |

Both wrong-way channels are confirmed, in direction and in class: the
engine is LOW by one to two decades where it extrapolates, and the error
GROWS as <eps> falls (Ar F-: x1.18 at 1.66 eV -> x1.47 at 1.43 -> x5.31 at
1.07 -> x25.2 at 0.761). The measured rates rise monotonically toward low
energy exactly as the 0.1.38 sub-floor audit predicted from the cross
sections alone (F-: 1.37 -> 1.70 -> 2.89 -> 5.15e-15 m3/s; vib: 2.75 ->
3.33 -> 4.73 -> 6.02e-15). The atlas's Maxwellian yardstick was the right
sign and about the right size: measured/Maxwellian is 0.43-1.12 across the
ladder, i.e. the true swarm sits within about a factor 2 of a Maxwellian
down here, and nowhere near the tapered table.
Opposite-sign finding: the HIGHER-threshold F2 channels are badly
OVER-predicted by the flat clamp — F2(a3) measured 7.4e-5 of the engine
value at <eps> 1.43 eV, F2(A1) 9.0e-7 at 1.66 eV (small absolute rates,
but the clamp is wrong in both directions).

### Added

- `MODULES/loki_tables/V2point_ArKrF2_subfloor/` (4 points, 0.05-0.3 Td)
  and `V2point_NeKrF2_subfloor/` (3 points, 0.1-0.5 Td), each with
  PROVENANCE.txt naming the setup files, the interpreter run, the wall
  clock, and the fact that a single-E/N LoKI run writes swarmParameters /
  rateCoefficients rather than lookUpTable files, so these two files were
  ASSEMBLED from those outputs value-for-value. Grids are disjoint from the
  base and _lowEN sets by construction. Setups kept in
  LoKI-MC/Code/Input/krf0d/sub_*.in; raw per-point outputs in
  LoKI-MC/Code/Output/krf0d_sub_*.
- `loki_tier.RUN_NAMES_{AR,NE}_FILL_SUBFLOOR` + `rates_for_mode(...,
  subfloor=False)`; `params.loki_subfloor` (default False); engine passes
  it; the three beam scorecard row builders take `subfloor=`.
- `tests/test_loki_subfloor.py` (10). Suite 234 -> 244.

### Fixed

- **Duplicate LOKI run names silently doubled every rate** (measured
  x2.0000). The 0.1.25 overlap guard compares directory names, so the same
  set listed twice passed it. Latent until now; reachable the moment a
  caller concatenates run-name lists, which `subfloor=True` does. Refused
  at construction with the measured factor stated. New Known-bad.

### Measured — what turning it on would do (the flag stays OFF)

Beam presets, ratio to the same case with the flag off:

| case | in-table frac | E_L | n_e peak | F- peak |
|------|---------------|-----|----------|---------|
| krf_beam_0d | 0.331 -> 0.405 | x1.011 | x1.00 | x1.10 |
| krf_ne_0d | 0.025 -> 0.388 | x1.077 | x0.572 | x0.932 |
| ne_optimized | **0.008 -> 0.992** | x0.915 | x0.556 | x0.925 |

The headline is the gauge, not the observable: the ne_optimized lasing cell
— the one behind every SPECKLE3D verdict — goes from 0.8 % of its pulse on
measured electron data to 99.2 %. Its efficiency moves 6.194 -> 5.669 % and
its peak density falls x0.56, so flipping the default would move the
recorded 3-D baselines and the reheat marginals computed against them.
That is the owner's call (the [0.1.36] precedent), and the flag is OFF.

### Measured — what it does NOT fix (an honest negative)

The NRL population anchor is essentially unmoved: F- x0.19 -> x0.20, and
Te x1.15, n_e x1.99, ArF* x0.28, Ar* x0.87, Ar** x0.05, F2 diss x0.43,
g0 x0.54, alpha x0.34 all identical to three digits. That case sits at
Te ~2.3 eV, well inside the tables, so the sub-floor region barely enters
it. Conclusion: the sub-floor extrapolation and the population deficits are
SEPARATE owners. The F- x0.19 and Ar** x0.05 misses survive the fix that
was their leading suspect, which points the F- deficit at the ion-ion
recombination screening (x1.75 form spread, OMISSIONS) and the DA data
itself (SIGLO vs TRINITI x1.58, 0.1.38), and Ar** at the missing
Ar*+e -> Ar** stepwise rows.

### Still open

- The arf_chen (Ne/Ar/F2) fill has no measured extension yet; `subfloor`
  is inert there and says so.
- Even the new floors (0.761 / 0.891 eV) sit above the settled beam-mode Te
  (0.04-0.63 eV), so the taper still carries the deep afterglow. Reaching
  it needs E/N below 0.05 Td, where the run cost and the decaying-swarm
  statistics both worsen; whether a swarm is the right object there at all
  is the open question ([0.1.37]: e-e collisions Maxwellianise that bulk).
- Verified: suite 244/244 under /usr/bin/python3; the assembled tables
  round-trip through LokiRates to the measured values; Homebrew smoke with
  -W error::DeprecationWarning.

## [0.1.38] — 2026-09-03 — the cross-section atlas: why the swarm shape is what it is, and which sigma actually governs the rates

Owner ask, following [0.1.37]: "look at graphing elastic cross section
against eV with the different gases to optimize our rate reaction
coefficients if the EEDF with swarm is better." Built as a diagnostic —
it reads the shipped LXCat sets and the shipped swarm tables and reports.
Nothing it computes enters an engine run, and no cross section was
touched (the D4 rule: a sigma is never adjusted toward an anchor).

### Added

- `tools/sigma_atlas.py` — the atlas. Seven figures into
  `rebuild_out/sigma_atlas/` (elastic sigma per gas with the
  Ramsauer-Townsend minima and table floors marked; momentum-transfer and
  energy-relaxation times at fill density; channels per gas;
  cross-database overlays; the k(Te) spread band; the real fills' cooling
  curve with the measured operating bands; and the loss-channel ranking)
  plus `report.json` with every number. Sections: `ramsauer` (a STATED
  criterion), `relaxation_times`, `mixture_loss_rates`, `database_rates` /
  `spread`, and `subfloor_audit`.
- `tests/test_sigma_atlas.py` (12): the Ramsauer criterion separates Ar/Kr
  (windows) from Ne/He (none); LXCat mass ratios match the atomic masses
  to 2 %; the relaxation identity nu_u = (2 m/M) n sigma v; Ar cannot cool
  through its window within a pulse while Ne can; trace F2 vibration
  outranks buffer elastic at 0.3 eV and the ordering reverses at 3 eV; the
  sub-floor audit finds exactly {F-, vib} and not the threshold buckets;
  every fill's table floor sits above the settled operating Te; the spread
  band is symmetric. Suite 222 -> 234.

### Measured — the elastic cross sections (3 atm, 300 K)

| gas | R-T minimum | depth vs 5-25 eV peak | 2m/M | tau_u at 0.3 eV |
|-----|-------------|-----------------------|------|-----------------|
| Ar  | 0.25 eV     | x165                  | 2.74e-5 | **998 ns**   |
| Kr  | 0.514 eV    | x185                  | 1.30e-5 | 349 ns       |
| Xe  | 0.585 eV    | x97                   | 8.40e-6 | 170 ns       |
| Ne  | (none)      | x6.0                  | 5.43e-5 | 69.8 ns      |
| He  | (none)      | x1.0                  | 2.74e-4 | 2.42 ns      |

Ne is NOT a Ramsauer gas — its shallow dip sits at the bottom of the
tabulation and its zero-energy scattering length is positive (A = +0.21
a0; a MERT zero needs A < 0), so the criterion (interior minimum, >= x20
deep) correctly refuses it. THE POINT: in Ar the elastic energy-relaxation
time in the window is x10 the 100-ns beam pulse. The buffer cannot take
the bulk's energy within a shot, which is exactly why a swarm at 0.5 Td
still reports <eps> = 2 eV rather than the 0.039 eV of a 300 K gas — and
why the beam-mode bulk sits hot and non-thermal. The swarm shape [0.1.37]
preferred is this physics; the atlas is its mechanism.

### Measured — which cross section sets the bulk temperature

Energy-loss rate per electron, inelastic vs elastic, in the real fills:

| <eps> | Ar/Kr/F2 | Ne/Kr/F2 | Ne/Ar/F2 | top channel |
|-------|----------|----------|----------|-------------|
| 0.1 eV | 0 (no channel open) | 0 | 0 | — |
| 0.3 eV | **x187** | x15.8 | x15.6 | F2(V2) 0.219 eV |
| 1 eV   | x6.27    | x2.74  | x2.67  | F2(V3) 0.324 eV |
| 3 eV   | x0.039   | x0.042 | x0.041 | F2(V3) |

At 0.1 % concentration, F2 vibrational excitation is the bulk's dominant
cooling channel wherever the buffer sits in its Ramsauer window; above the
window the ordering reverses. So the sigma that governs the sub-floor rate
coefficients is not the buffer's elastic set at all — it is F2's
vibrational and attachment sets. Below 0.111 eV (the V1 threshold) NO
inelastic channel is open and only the collapsed elastic term remains.

### Measured — the sub-floor audit (the actionable finding)

The engine extrapolates every bucket below the swarm-table floor by one
rule: flat for 0.35 ln-units of <eps>, then an exponential taper — always
DOWN. Against the cross sections (Maxwellian as the common yardstick,
which at these sub-threshold energies is also the defensible shape),
exactly 2 of 14 buckets physically rise there, in all three fills:

| fill | table floor <eps> | F- low by | vib low by |
|------|-------------------|-----------|------------|
| Ar/Kr/F2 | 1.98 eV (Te 1.32) | **x759** | **x142** |
| Ne/Kr/F2 | 3.06 eV (Te 2.04) | x3670 | x728 |
| Ne/Ar/F2 | 2.26 eV (Te 1.51) | x1212 | x229 |

(at <eps> = 0.3 eV, the beam-mode operating point). The true vib rate
peaks at Te 0.414 eV — x1.66-2.45 ABOVE its value at the table floor —
and F- keeps climbing to the bottom of the probe grid (x7.0-11.6). The
tables see only the decreasing branch of each and the taper continues it
downward. This is the SAME signature the scorecard already records for
F- ("the F- DA curve RISES x13 toward low eps where the taper fell",
cured for the KINETIC tier by deepening its measured curves); the atlas
finds vib carries it too, which was not previously recorded, and vib is
the channel that sets Te.

### Measured — what the database choice is worth (an error bar, not a lever)

Maxwellian-averaged k(Te), max/min across published evaluations:

- F2 SIGLO vs TRINITI: elastic x1.99, ionization x1.91, excitation x1.69,
  **attachment x1.58** (worst at Te 0.05 eV).
- Ne Biagi-8.9 vs BSR: ionization x2.14, excitation x1.82, elastic x1.10.
- Kr Biagi-v7.1 vs the lokimc copy: byte-identical, no band.
- Ar: only IST-Lisbon is in-tree, so no band can be formed.

The data spread is a factor ~2. It does NOT explain the F- x0.19 or
Ar** x0.06 anchor deficits, so swapping databases is not the fix — but it
is a real, citable uncertainty the UQ layer should carry per channel.

### VERDICT and the one action

- The elastic sets answer WHY the swarm EEDF is right and where the bulk
  sits; they are not themselves the lever on the rate coefficients.
- The lever is the SUB-FLOOR REGION, and the honest fix is to MEASURE it,
  not to re-shape an extrapolation: LoKI-MC swarm runs below 0.5 Td (the
  present floor) down to ~0.01 Td would carry <eps> from 2 eV to near
  1.5 kT_gas = 0.039 eV and cover the whole operating band. That is a
  data extension of the existing `_lowEN` sets with the same provenance,
  and it removes the extrapolation the engine currently runs on for
  67-99 % of a beam pulse ([0.1.37] te_in_table_frac). Risk to check in
  the run: at low E/N an attaching mixture has no steady swarm, so the
  MC statistics and the decaying-population handling need watching.
- Stated, NOT implemented (an owner decision): the sub-threshold buckets
  (vib, F-, elastic) have a defensible Maxwellian closure at beam
  densities — e-e collisions Maxwellianise the bulk, and below the first
  threshold there is no inelastic depletion to sustain a swarm tail —
  while the threshold buckets keep swarm rates. That is a per-channel
  hybrid, and the [0.1.37] MaxwellShaped machinery already expresses it.
- Verified: 234/234 under /usr/bin/python3; atlas + figures regenerate in
  ~40 s; two new Known-bad entries.

## [0.1.37] — 2026-09-03 — the MAXWELL verification tier: the EEDF-shape assumption scored

Owner ask: "a background verification/approximation module that does a
two-temperature Maxwellian / sum of two Maxwellians" — "I just want a
verification ... do it that way as an approximation." Built as a
side-by-side tier (the S4 class): it scores the LOKI swarm-shape closure,
steers nothing, and is never a default.

### Added

- `MODULES/maxwell_tier.py` (layer 3): `maxwell_rate` (Maxwellian-averaged
  k from sigma(eps), gated against the closed forms), `MaxwellRates`
  (per-bucket k(Te) and `k_two_t` at FIXED mean energy from the LXCat sets
  cross_sections.py resolves — the same files behind loki_tables/ —
  bucketed exactly as loki_tier._bucket), `MaxwellShaped` (the engine
  tier: Maxwellian k inside the base table's <eps> range, the LOKI
  flat+taper rule below its floor so both tiers extrapolate identically
  where neither is measured, E/N lookups delegated, `k_by_en` refused).
- Engine: `electron_tier="maxwell"` (params choices; `run.py --tier
  maxwell`). Beam drive only — the circuit lane is refused at construction
  by capability with the fix stated; the relax bridge is skipped on this
  tier so the comparison isolates the shape.
- scorecard: `kannari_device_rows` / `electra_arf_rows` /
  `electra_population_rows` take `tier=` and `subtable=` (defaults
  unchanged).
- `tools/maxwell_verify.py`: rate level — k_MAXWELL/k_LOKI at the tables'
  own <eps> points, the two-Maxwellian bracket over a STATED tail grid
  (default 1/3/10 % at 2/4/8 eV; `--tail a:T`), and the E/N view for the
  discharge lane (k_by_en against the Maxwellian at <eps>(E/N)); device
  level — `--device` (presets + ne_optimized), `--population`,
  `--kannari` under loki / loki-taper / maxwell, with the engine's
  `te_in_table_frac` gauge printed; JSON report in rebuild_out/.
- `tests/test_maxwell_tier.py` (18): pdf normalisation; constant / step
  sigma closed forms (1e-4 / 3e-3); LOKI-bucket parity for the Ar, Ne and
  Chen fills; two-T reduction + mean-energy conservation; monotonicity;
  table-top order-of-magnitude gate; wrapper delegation + refusal; the
  below-floor rule (flat, then taper, continuous); engine run closes and
  is deterministic on the tier; circuit refused. Suite 203 -> 221.

### Measured — rate level (k_MAXWELL / k_LOKI at the tables' own <eps>)

- Threshold processes at the tables' LOW end (<eps> 2-6 eV): x1e3-1e6.
  Ar fill: Ar* 2.0e5 @3.8 eV, Ar+ 8.9e5 @5.4, Kr* 1.4e5 @3.2, Kr+ 5.4e4
  @5.1; Ne fill: Kr* 509 @3.1, Ne* 4.5e5 @4.5, Ne+ 8.2e4 @5.5; Chen fill:
  Ar* 1.6e3 @3.2, Ar+ 5.1e4 @4.6, Ne+ 8.5e3 @6.5. At the HIGH end
  (8-15 eV): x0.9-2.8 (Ar+ 2.56 @8.4, Kr+ 2.14, Ne** 2.3 @14.9, Ar* 0.92
  @14.5). This is the inelastic wall: a field-driven swarm EEDF is cut
  off above the first excitation threshold, a Maxwellian is not; the
  low-end magnitude is the Townsend class (Ar at ~5 Td: (alpha/N) v_d
  ~2e-21 m3/s against a Maxwellian ~2e-16 at the same <eps>). The lowest
  tabulated points of the high-threshold buckets are also MC-statistics
  limited (2e4 electrons), so the ratios there characterise the table
  floor as much as the shape.
- Sub-threshold processes agree within x2 everywhere: F- (F2 DA)
  1.6-2.3, vib 1.4-1.6, elastic 0.8-1.3, F2(a3) 0.7-2.5. The
  two-Maxwellian bracket at fixed <eps> raises threshold rates further
  (to x1e6) and sub-threshold ones by at most x4.6 (F-): a tail can only
  add high-energy electrons.

### Measured — device level (beam drive; ratios to the LOKI default)

- The LOKI runs' own `te_in_table_frac`: 0.33 (krf_beam_0d), 0.025
  (krf_ne_0d), 0.008 (ne_optimized) — the beam-mode bulk sits BELOW the
  tables' floors nearly all the time, so the standard closure there is
  the relax bridge / flat+taper; the loki-taper column prices the bridge:
  E_L x0.84 and F- x0.37 on krf_ne_0d, n_e x1.29 and F- x0.52 on
  ne_optimized, ~x1.00 on the Ar fill.
- MAXWELL: Te peak x0.15-0.30 (2.2-2.5 -> 0.34-0.74 eV; in-table 0 %):
  the Maxwellian tail drains the electron energy through the thresholds
  1e3-1e5x faster and the fluid Te cannot hold. Then E_L x0.89 (Ar
  sub-threshold cell), x1.10 (Ne), x1.20 (ne_optimized 7.77 -> 9.36 J,
  eta 6.19 -> 7.46 %); n_e peak x0.96 / x1.42 / x2.09; F- x0.03-0.06
  (the below-floor taper of a rate that physically RISES toward low
  energy — the 0.1.13 starvation mechanism, on this tier); Kr**
  x1.2-2.1, Ar* x1.24, Ar** x1.14 (Ar fill).
- NRL population benchmark (petrova2022; engine/NRL): LOKI Te x1.15,
  n_e x1.99, ArF* x0.28, Ar* x0.87, Ar** x0.05, F- x0.19, F2 diss x0.43,
  g0 x0.54; loki-taper identical except F- x0.08; MAXWELL Te x0.20, n_e
  x3.10, ArF* x0.24, Ar* x1.92, Ar** x0.06, F- x0.00 (8e9 cm-3), F2 diss
  x0.44, g0 x0.47, laser out 22.6 J vs 55.9 J.

### VERDICT

- The swarm-shape closure is the better approximation for the beam-mode
  bulk: NRL's EEDF-resolved model reports Te ~2 eV flat, which LOKI meets
  (x1.15) and a Maxwellian cannot (x0.20). A two-temperature Maxwellian
  only adds tail and moves every threshold rate further from the swarm;
  it is not the direction for accuracy. The physics: at n_e ~1e15 cm-3
  e-e collisions do Maxwellianise the sub-threshold bulk, but the tail
  above ~10 eV is depleted by inelastic collisions ~1e3x faster than e-e
  can refill it, so the swarm's inelastic wall is the closer shape.
- The population misses are NOT the bulk shape: Ar** stays x0.05-0.06
  under every tier (owner: the Ar*+e -> Ar** stepwise row class,
  OMISSIONS), and the sub-threshold shape sensitivity of F2 DA is
  <= x2.3 against the x5 F- deficit (owners: ion-ion recombination
  screening x1.75; the sub-floor DA extrapolation, since the DA curve
  rises toward low energy where both tiers taper).
- What the bracket is good for: an EEDF-shape band per bucket for the
  UQ layer, and the measured statement that every beam-mode answer
  depends on the sub-floor rule for 67-99 % of the pulse. What removes
  that dependence is an in-loop Boltzmann solver with the beam source
  term (COMMERCIAL_AUDIT.md section 6, item 6), not a Maxwellian family.

### Numerics (recorded)

- The BARE Maxwellian in the afterglow trips scipy's num_jac factor
  overflow harder than the tables do: krf_beam_0d 100 ns in 7.7 s, the
  next 40 ns in 423 s (28,286 RHS, 9,198 LU factorisations). With the
  LOKI below-floor rule the 200-ns run integrates in 9.5 s (LOKI 8.9 s);
  restart segments were not needed.
- Verified: suite 221/221 under /usr/bin/python3 (534 s, 2 BLAS
  threads); Homebrew 3.14.6 import + 20-ns maxwell run with
  `-W error::DeprecationWarning` (closure 2.1e-12); `np.trapezoid`
  resolved for numpy >= 2. `tools/maxwell_verify.py --device
  --population` runs in ~3 min.

## [0.1.36] — 2026-09-01 — beam-EEDF e+F2 channels (Petrov 2017), built gated and priced

Owner ask: "tackle the open bugs/limitations one by one and see if we can
improve our agreement to the experimental data." Item 1 = the scorecard's
own stated limit of the Electra anchor: the dissociation/vibration/
attachment of F2 by beam secondaries (petrova2022 Table A2; the
Kushner-Moratz 1988 class), read first-hand this session along with their
primary source Petrov 2017 JAP 122 133301 (Table I: k(P, x_F2) fits over
the beam-sustained EEDF, their Fig. 6 — the fits FOLD the secondary tail
into the coefficient, so they REPLACE the <eps>-indexed swarm rates for
those channels, never add to them).

### Added — `petrov17_f2` (default OFF)

- chemistry.petrov17_f2_k: the three Table I fits (vib / att / dis),
  cm^3->m^3, with the paper's own P >= 1e-3 floor. Gate-tested against
  hand-derived literals from the printed table, not against itself.
- CompiledNetwork marks the target rows by PATTERN (reactants e+F2,
  loki-rated, buckets F- / F2(a3) / F2(A1)) — immune to rid drift and
  the 0.1.14 repulsive reroute. The total dissociation coefficient
  branches across the two repulsive rows in proportion to their own
  swarm rates (shape kept, magnitude replaced). The vib SINK swaps at
  its own block. Applied ONLY while the beam envelope is on (engine
  hoists p_dep above the chemistry call; beam-off reverts to the swarm
  rates — one owner per electron population; E/N-direct bookkeeping
  untouched, that path never carries the flag).
- Gates: 5 new tests (formula anchors, limits, marked-row discovery,
  OFF-bit-inert + ON-swaps-only-under-beam, ledger closure ON). Suite
  203/203. One test's own first draft wrongly asserted dis monotonic in
  x_F2 — the printed form is non-monotone; the TEST was fixed.

### Measured — flag ON vs the 2026-09-01 baseline (nothing tuned)

- Population benchmark: F2-burn x0.43 -> x0.53 (0.266 vs NRL 0.5 by
  pulse end — the mechanism moves the right way), alpha x0.34 -> x0.39,
  ArF* x0.28 -> x0.30, Te x1.15 -> x1.22, F- x0.19 flat, Ar** x0.06
  flat, n_e x1.53, g0 x0.56.
- Kannari: x0.83 -> x0.86 (R_OC 80%), x0.20 -> x0.23 (20%, the x9-x11
  threshold amplifier).
- Electra ArF: x0.84 -> x1.06 / x1.41 -> x1.62 / x1.76 -> x2.03; trend
  still NOT reproduced.

### VERDICT — a convergence result, not a fix

With the same class of beam-EEDF F2 channels Orestes uses, this engine
now lands where Orestes lands: matching the 0.82-atm point (x1.06 vs
their x1.02) and over-predicting the high-pressure points — wolford2020's
own words for their model. The Electra high-pressure decline's owner is
therefore SHARED model physics (NRL's stated 3-body suspicion, the
ARF_AA class), not the missing beam channels: that suspect is now
consumed for the trend while remaining real physics for the populations.
Default stays OFF pending the owner's call — flipping it moves every
beam anchor and the recorded SPECKLE/SPECKLE3D baselines; the numbers
above are the decision inputs. (Verification: suite 203 passed under
/usr/bin/python3, 12 min on the loaded box; anchors measured through the
scorecard's own row functions with the flag injected via
engine_parameters.)

## [0.1.35] — 2026-08-28 — the multi-cell cure, the first 3-D verdict, and grinder closeout

### Fixed / corrected

- [0.1.34]'s flagged open item is now DIAGNOSED AND CURED (measured in
  SPECKLE3D, CHANGELOG 2026-08-28): the multi-cell Ne-fill blow-up is
  scipy num_jac FACTOR-state overflow (the 0.1.21 trap class,
  mid-pulse) on BOTH the sparse and dense jac paths — [0.1.34]'s
  "dense is the workaround" attribution was falsified by a dense 2x2
  failure at t = 106 ns. Cure: fine 30-ns checkpoint restarts through
  the pulse (SPECKLE3D driver default); the cured 2x2 sweep completed.
  An engine-side look at num_jac robustness for multi-cell Ne fills
  remains open (workaround-level fix only).

### Measured — first 3-D twin verdict (2x2x1, SPECKLE3D)

- native +10.778% (0-D +10.748% reproduced to 0.03 pp); shell = native
  to 1e-4 (equal-weight consistency check of the [0.1.34] hook,
  passed in flight); **port — all power into half the volume — keeps
  +9.565%, ~89% of the marginal**: the late-extraction gain survives
  strong lateral deposition localization. 3x3/4x4 pending (owner-run).

### Removed — the kinetic-tier grinders (owner ask: terminate unneeded shells)

- Terminated 2026-08-28 without completing: the 2026-08-21 full
  kinetic-tier sweep (7 days, pre-errata, stuck on ne_optimized after
  5/8 cases) and the 2026-08-26 3-case tail (2 days, same case). The
  [0.1.26]/[0.1.28] "still integrating" notes are closed: kinetic-tier
  ne_optimized does NOT converge in practical time under the current
  segmented integrator; **kinetic 5/8 is the standing record**, and
  any future completion needs its own attack (coarser rtol probe or an
  engine-side stiffness fix), not more wall-clock.

## [0.1.34] — 2026-08-27 — stated microwave deposition weights + the SPECKLE3D workspace

Owner ask: "I just made a subdirectory under Kinetic_Modeling called
SPECKLE3D, can we extend SPECKLE from 1D to 3D and see if it still
gives a small fractional percent increase?" The honest instrument is
REBUILD's own multi-cell grid (independent reactors sharing the
cavity/ASE fields — 3-D inhomogeneity + shared extraction, NOT
transport; Nathaniel's 1-D model cannot be taken to 3-D and stays
untouched).

### Added — engine

- `Device.set_mw_cell_weights(w)` (`engine.py`, default-inert): stated
  per-cell microwave DEPOSITION weights — a skin/shell profile, or a
  measured / full-wave (COMSOL-class) deposition map. The window's
  power balance settles against the weighted conductivity and deposits
  per cell as w_i n_e_i; None (default) = the pre-hook conductivity
  weighting, bit-identical. Gate
  (`test_mw_cell_weights_identity_masking_and_conservation`): w = 1 is
  bit-identical to the default; a masked cell receives exactly no
  microwave heat; the drive-ledger rate matches the prescribed power.
  MEASURED BY THE GATE while writing it: delivered power matches the
  prescription to the settle grid's linear-interpolation accuracy
  (~3e-6 relative, weighting-dependent at that scale) — the docstring
  states it; the ledger books the ACTUAL deposit exactly.

### Added — SPECKLE3D (workspace seeded; see SPECKLE3D/CHANGELOG.md)

- `SPECKLE3D/mw_3d_study.py` + README + CHANGELOG: the ne_optimized
  afterglow twin on nx x ny x nz grids with three stated deposition
  geometries (native conductivity / 2-mm shell off all transverse
  walls / single-port skin off one wall), shared off-run per grid,
  dense-Jacobian segmented integration.

### Measured / Observed

- Sanity (1x1x1) reproduces the recorded channel digit-exact
  (off 7.7732 J / 125.49 J; native twin dE +8.5981 mJ = the [0.1.31]
  ensemble nu_hat, marginal +10.75% at the tools' 115-ns window).
- FLAGGED OPEN ITEM: the 4x4x1 Ne-fill off-run under Device.run()'s
  SPARSE numerical Jacobian trips the [0.1.27] non-finite guard at
  beam-on (Ar/Ar+/Ar2+/U_e rows — near-zero-Ar columns under grouped
  finite differences; 1x1x1 dense clean; historical 3-D runs were Ar
  fills). Worked around by dense integration in the study; the
  sparsity pattern on Ne-fill grids needs a look before any Ne-fill
  3-D run uses Device.run().
- The 4x4x1 native/shell/port sweep itself was launched and EXTERNALLY
  STOPPED mid-run (heavy dense-Jacobian compute on a machine busy with
  other builds) — the 3-D verdict is PENDING; relaunch:
  `python3 SPECKLE3D/mw_3d_study.py --skip-sanity`.

### Verified

- Full suite under /usr/bin/python3: **198 passed** (5:38; 197 + the
  new gate). Homebrew python3 -W error::DeprecationWarning smoke: a
  weighted 2-cell device runs clean end-to-end. The gate's first two
  drafts were themselves falsified and corrected by the run (exact-
  conservation and cross-weighting-equality assumptions both fail at
  the settle grid's ~1e-6 interpolation scale — now stated in the
  docstring instead of assumed away).

## [0.1.33] — 2026-08-27 — the derated operating point (owner Q&A before the talk)

Owner asks, verbatim class: what IS the efficiency with/without; what
does "signs not verified" mean; is the 7% bar real; and what operating
point survives realistic coupling losses (IB < 100%, walls, ...).

### Added

- `tools/derated_operating_point.py`: pure arithmetic on the recorded
  artifacts. The exact with/without (Ne-optimized afterglow: eta
  6.194% -> 6.197%, i.e. 7.773 J/125.49 J -> 7.781 J/125.57 J), the
  first-order identity d_eta = (eta_c*m − eta0)*E_RF/E_tot (checked
  against the recorded +0.00244 pp digit-for-digit), the INTRINSIC
  coupling threshold eta_c > eta0/m = 62%, the WALL-PLUG requirement
  eta_src*eta_c > 0.70 at the program-stated 7% bar (magnetron needs
  97% coupling; X-band/gyrotron need >100% — excluded; at an 8% bar
  excluded outright), Drude–Fresnel reflection off the collisional
  plasma face (nu 2e12: 2.45 GHz reflects 5/42/76% at n_e
  1e19/1e20/1e21; 95 GHz reflects 0/0/20% — absorption kinetics is
  frequency-blind but COUPLING favors the high band), and the blended
  gain table over eta_c x E_RF. rebuild_out/mw_benchmarks/
  derated_operating_point.json.

### Measured — the derated recommendation

- DEMONSTRATION POINT (intrinsic accounting): the Ne-optimized cell
  (eta0 6.19%, margin ~19), afterglow window beam-off+10 ns x 150 ns,
  95 GHz launch (Fresnel-favored; the gyrotron's low source
  efficiency is irrelevant to intrinsic bookkeeping), eta_c 0.75–0.9
  plausible vs the 62% threshold. Expected d_eta +0.001 to +0.04 pp
  for E_RF 0.1–2 J: an mJ-class, 16.5-sigma-measurable dE_L — a
  demonstration, not a device lever. No parameter set survives the
  wall-plug criterion at the measured marginal.
- The 7% bar's provenance is now stated ON the slide (program-stated
  current e-beam wall-plug, docs/DE-AR_targets.txt), with sensitivity
  (5%/8% bars) printed by the tool. The backup-slide "signs" note
  reworded: sign-not-tier-stable = the few-percent marginals flip
  between the two independently derived rate sets while the 0.006-pp
  ceiling holds on both — a rate-knowledge statement, not a suspected
  bug. New backup slide carries the blended-gain arithmetic
  (mw_reheat_slides now 13 pages, recompiled clean, inspected).

## [0.1.32] — 2026-08-27 — the group-talk deck

Owner ask: "make this into a few slides in BEAMER to present to my
colleagues tomorrow."

### Added

- `SPECK/PAPER/mw_reheat_slides.tex`: 11-frame Metropolis deck (16:9)
  + one backup — the paper's argument in talk order: question/method
  (twin pairs, three engines, decidability ladder), the answer in
  three lines, the observable (compact benchmark table + the
  don't-mistake-fluorescence-for-efficiency trap), the real channel
  (mechanism + 41/41 bracket strip), the 1-D fingerprints
  (margin/timing panels), the 16.5-sigma placebo-calibrated artifact
  certificate with its scope line, the five walls (+ wall-plug bars),
  the one-rate-wide co-pump sign (KSTEP window + SPECKLE dt series),
  what-to-do-next, summary, and the full sweep table as backup. Title
  slide carries the internal/data-rights note.
- `SPECK/PAPER/make_slide_figs.py`: single-panel, slide-scale
  regenerations of the paper figures (s_margin/s_timing/s_bracket/
  s_wallplug/s_kstep/s_dtseries.pdf) from the SAME recorded JSONs —
  no hand-typed numbers.

### Verified

- pdflatex x2 clean (12 pages); every frame visually inspected (one
  overflow found on the five-walls frame and fixed + re-inspected);
  all numbers cross-read against [0.1.26]-[0.1.31] and the paper
  draft at their recorded precision.

## [0.1.31] — 2026-08-27 — the artifact question, closed statistically: a placebo-calibrated 16.5 sigma

Owner ask: "use a CL_S method or profile likelihood ratio to give me a
3 sigma confidence that this is not a numerical artifact or bug in my
math or code. MAKE NO MISTAKES."

### Added

- `tools/artifact_significance.py`: the [0.1.28]-B ne_optimized
  afterglow twin re-measured under SEVEN solver settings (BDF rtol
  1e-6/1e-8/1e-10, atol x0.1, max_step 1 ns, Radau 1e-6/1e-8;
  segmented), each three ways — off, 500 kW signal, and a 1-W PLACEBO
  window whose twin difference (first-order physical part subtracted)
  calibrates the numerical differencing noise in situ. Statistics:
  exact Gaussian profile likelihood (q0 = (n+m) ln(sigma0^2/sigma^2),
  sigma profiled), EXACT Student-t discovery test (T = nu_hat/(s_p/
  sqrt(n)) ~ t_{n+m-1} under H0), 2e5 toy MC, central-t CLs-allowed
  band, worst case over the 41 [0.1.28] bracket variants, and the
  noise-inflation factor to degrade to 3 sigma. SELF-TESTS RUN FIRST
  and abort on failure (q0 algebra vs direct optimization, t/CLs
  closed forms, synthetic discovery, synthetic null, band-width
  formula); a split-sample NEGATIVE CONTROL on real placebo data is
  gated at |T| < 3.5. rebuild_out/mw_benchmarks/
  artifact_significance.json.

### Measured — the certificate

- Signal nu_hat = +8.5981 mJ over the 7 settings (BDF and Radau agree
  to 2e-8 relative); placebo-calibrated noise s_p = 1.58e-4 mJ.
- **Exact discovery: T = 1.44e5 (13 dof) -> p0 = 1.7e-61 -> Z = 16.5
  sigma** (the honest small-sample conversion — the t tail, not the
  naive normal z). PLR q0 = 296.6 (Z_asym 17.2); conservative
  single-measurement z = 5.4e4; toys 0/200000 exceedances (p0 < 5e-6,
  Z > 4.4 at toy resolution; toy T vs t13 KS p 0.60).
- 95% CLs-allowed signal band [8.5980, 8.5982] mJ; nu = 0 sits 1.4e5
  standard errors below it. Worst single bracket variant of [0.1.28]
  alone: z = 4.2e4. The calibrated noise would have to be
  underestimated x3.9e4 to reduce the discovery to 3 sigma.
- HONESTY NOTES (the audit): the placebo differences cluster at a REAL
  +2.1e-4 mJ — the 0-D small-signal response is much larger per watt
  than at 500 kW (linear-response marginal ~+140% vs +10.7%: the
  response is strongly sublinear in P). Treating that real offset as
  noise only INFLATES s_p (~x40 above the true solver scatter), so
  every Z above is conservative. The negative control (T = +2.83,
  p = 0.01, gate |T| < 3.5) detects exactly that offset — recorded,
  not hidden.
- 1-D (SPECKLE 2026-08-27c placebo ladder): the deterministic
  differencing floor, measured by powers 250 kW -> 250 W, is
  7.8e-4 mJ (plateau reached; the bottom of the ladder is
  non-monotone, the floor signature). The 250-kW signal exceeds it
  x302 at the 300-ns horizon, x~1.9e3 at the 500-ns horizon.
- SCOPE, stated plainly: this rejects a numerics-noise origin of the
  twin difference. No statistic of one model's outputs can certify
  the absence of a shared conceptual error in the physics; that
  burden is carried by the independent-implementation agreement (two
  unrelated codebases, same-sign ~+10% channel), the closed energy
  ledger (1e-13), the 197-test suite, and the bracket program.

### Verified

- Self-tests all pass before any data are touched; toy distribution
  consistent with the exact t13 (KS p 0.60); two independent
  integrators inside the ensemble; run live 2026-08-27 from the MAIN
  checkout on the post-[0.1.27] engine. The paper (SPECK/PAPER) gains
  the corresponding paragraph in Sec. III B and was recompiled clean.

## [0.1.30] — 2026-08-27 — the 1-D verification: the paper's channel survives spatial structure

Owner ask: "can we use SPECKLE, which combines with 1DModelKrF, to
verify the reheating works in the 1D case as well." Measurements live
in SPECKLE (CHANGELOG 2026-08-27b: the --mw-profile deposition knob +
11 new twin pairs); this entry records the verification logic and the
paper update.

### Measured — three predictions of the late-extraction mechanism, all
### confirmed in the independent 1-D model

1. MARGIN: the afterglow marginal rises monotonically with the output
   coupler (+0.01 → +0.20 → +0.78 → +0.89% for r_out 0.3 → 0.9) while
   the ABSOLUTE laser output is non-monotone (peak at r_out 0.6) —
   the marginal follows photon lifetime (drain overlap), not laser
   performance. That dissociation is the mechanism's fingerprint.
2. TIMING: at the high-margin point with the full drain captured
   (500-ns horizon) the marginal is +5.74% (dt 50 ps) → +10.06%
   (dt 25) — robustly positive, the same ~+10% class as the 0-D
   Ne-optimized channel — and decays to ≈0 as the window slides past
   the drain (+1–3% at t0 250 ns once timestep-refined; the dt-50
   −4.4% there was a discretization artifact, resolved by Δt/2 and
   Δt/4 repeats; −0.2% ≈ 0 at t0 330).
3. LOCALITY (1-D-only): re-depositing the same power as a 2-mm edge
   skin vs center-weighted vs uniform moves the marginal ≤ ~0.3 pp
   (edge vs center: 0.04 pp) — WHERE across the gap is second-order,
   WHEN is first-order. The uniform idealization slightly
   under-estimates the marginal; the locality error is bounded at the
   Δt-drift scale.

### Changed — the article draft

- SPECK/PAPER/mw_reheat_draft.tex: new "One-dimensional verification"
  paragraph in Sec. III B + Fig. 2 (fig_1d_verification.pdf, generated
  from the ladder reports by make_figs.py); abstract claim (ii) now
  carries the 1-D reproduction; Limitations rewritten (the 1-D profile
  ladder bounds deposition locality; launch-geometry effects remain
  outside both models); Discussion adds the marginal-vs-delay curve as
  the bench signature. Figure environments reordered so numbering
  follows first citation (robust=1, oned=2, during=3).

### Verified

- pdflatex ×2 clean after the edits: 14 pages, zero errors, zero
  undefined references; fig_1d_verification regenerated from the
  recorded report JSONs (no hand-typed numbers) and inspected (two
  render passes: clipped markers/labels found and fixed). Ladder
  numbers quoted here are read from SPECKLE/runs/*_report.json; the
  refinement pairs (tw_t170_dt25, tw_t250_dt25/dt12) are the
  falsification runs for the two claims that needed them.

## [0.1.29] — 2026-08-27 — the in-silico article draft

Owner ask: "write up a quick draft for an In Silico / Computational
Article describing that it is feasible to reheat as an observable and
that it can be done and is not an artifact but is very difficult."

### Added

- `SPECK/PAPER/`: `mw_reheat_draft.tex` (REVTeX 4.2 aps/pra/preprint —
  "Microwave reheating of electron-beam-pumped rare-gas–halide laser
  plasmas: a robust observable, a real but marginal channel"),
  `make_figs.py` (both figures regenerated from the RECORDED artifacts:
  sign_decidability.json, copump_verify.json, ib_feasibility.json, the
  SPECKLE co-pump reports — no hand-typed numbers), the two figure
  PDFs, the compiled draft PDF, and a README with the
  before-circulation checklist. Structure: the three-part claim
  (observable 8/8 both tiers + SPECKLE; the 41/41 bracket-robust
  late-extraction channel with its mechanism; the five walls), the
  decidability ladder as methods, Table I = the [0.1.26] benchmark
  table + fluorescence column, Fig. 1 = KSTEP window + SPECKLE
  dt-series, Fig. 2 = 41-variant strip + wall-plug bars. Placeholders
  are explicit red [TODO]s: citations, affiliation, funding text
  (program-approved statement to be inserted VERBATIM — the repo's
  frozen-funding-text rule), 1-D-model author permission/citation.
  EXPORT-CONTROL REVIEW REQUIRED before any circulation (header
  comment): program-traceable parameters appear as modeling choices
  only; the controlled document is neither quoted nor cited.

### Verified

- pdflatex ×3 clean: 12 pages, zero errors, zero undefined references
  (final pass grep). Figures generated by script from the JSONs and
  visually inspected page-by-page (title/abstract, Table I, Secs.
  III C–IV) before delivery; the log-axis minor-tick and label-collision
  defects found in the first render were fixed and re-inspected. All
  quantities in the text were cross-read against the recorded tables
  in [0.1.26]–[0.1.28] at their stated precision.

## [0.1.28] — 2026-08-27 — can the model tell? the reheat signs' decidability, measured

Owner ask: "Is there a way to tell if SPECK or SPECKLE indicate a
possible realistic reheating at all?" The answer machinery: demand
invariance across numerics, rate tiers, rate brackets, and model
family — that ladder already DECIDES three of the four reheat
questions (afterglow energy harvest: no, tier-consistent + audited;
fluorescence response: yes, 8/8 both tiers + SPECKLE; magnitude at
feasible coupling: ≤0.006 pp everywhere). This entry measures the
decidability of the fourth (the co-pump sign) and bracket-tests the
one afterglow-positive case.

### Added

- `tools/mw_sign_decidability.py`: two discriminators on recorded
  configurations verbatim. (A) the plateau co-pump marginal vs the
  KSTEP multiplier (the 0.1.23 point: krf_hyper_mw, 2 atm,
  100 A cm⁻², window 35+35 ns, 250 kW, LOKI) — 0.1.23 had shown KSTEP
  (e + Kr* → Kr⁺ + 2e) is the ONLY sign-flipping bracket, negative at
  BOTH ×⅓ and ×3 around a positive nominal, so two zero crossings
  exist; this maps them. (B) the full 0.1.23-style EST ×⅓..×3 +
  LIT-conflict bracket pass on the ne_optimized AFTERGLOW twin (the
  only benchmark whose afterglow marginal is laser-positive),
  segmented integration. rebuild_out/mw_benchmarks/
  sign_decidability.json.

### Measured (A) — the co-pump sign is a ONE-RATE question, and the
### price of deciding it is a factor-2 rate determination

KSTEP map (post-errata engine, from the main checkout; anchors
reproduce 0.1.23: ×⅓ −1.65% vs recorded −1.66%, nominal +8.0% vs
+7.9%, ×3 −5.62% vs −5.62%):

| ×KSTEP | 0.333 | 0.45 | 0.6 | 0.8 | 1.0 | 1.25 | 1.6 | 2.2 | 3.0 |
|---|---|---|---|---|---|---|---|---|---|
| marginal | −1.7% | −1.3% | +2.9% | +5.4% | +8.0% | +15.5% | +6.6% | −3.1% | −5.6% |

- Zero crossings at ×0.49 and ×1.99: **the marginal is positive iff
  KSTEP lies within [×0.49, ×1.99] of the network's nominal** — a
  window of ratio 4.0, geometrically centered almost exactly on ×1
  (√(0.49·1.99) ≈ 0.99). Internal structure: a +15.5% peak at ×1.25
  (an optimum stepwise-ionization rate exists; reported, not
  mechanized).
- DECIDABILITY VERDICT: the EST tier's own ×3 uncertainty convention
  straddles BOTH crossings → the sign is UNDECIDABLE at current rate
  knowledge — but a determination of this single rate to better than
  ~×2 decides it within the LOKI framework. The kinetic tier's
  nominal-point negatives (−4.5% at 0.1.23's window, −2.7% at the
  benchmark window) then say the cross-tier spread in KSTEP-class
  effective physics currently EXCEEDS that window — consistent, and
  it names the deciding measurement: pin the rate. Paths, in cost
  order: a literature adjudication of e + Kr*(4s/4s') ionization
  cross sections; the 0.1.23-proposed field-on relaxation sampler in
  tools/export_kinetic.py (measures the KINETIC tier's own KSTEP
  response); the bench (docs/test_range.txt Tests 1–3).

### Measured (B) — the one afterglow-positive case is NOT rate-hostage:
### 41 of 41 brackets keep it positive

- Every EST-active row at ×⅓ and ×3 plus every applicable LIT-conflict
  alternate on the ne_optimized afterglow twin (window 115+150 ns at
  500 kW — this tool's beam-end estimate; the [0.1.26] benchmark row's
  window started 5 ns later, same afterglow class): **41/41 variants
  POSITIVE, marginal min +8.28% (NK6 ×3) / max +12.85% (NK6 ×⅓)**,
  821 s, one element-bound EST row structurally zero. KSTEP — the rate
  that owns the co-pump sign — moves this marginal only 10.3 → 11.0%:
  the afterglow-class positive does not share the co-pump's hostage.
- Interpretation: the high-margin still-lasing afterglow reheat is the
  models' one BRACKET-ROBUST laser-positive, with cross-framework
  witnesses of the same class (ar_optimized on the KINETIC tier +9.5%;
  SPECKLE's dt-robust +0.4–1.1%). Standing caveats: the SAME-cell tier
  check (ne_optimized on kinetic rates) is still integrating in the
  background sweeps; the blended payoff stays ≈ +0.002 pp at feasible
  coupling — real physics, not an efficiency lever ([0.1.26]'s
  wall-plug arithmetic still caps it); and dE_laser ≫ dFluor here, the
  no-go audit's late-co-pump mechanism (the cavity is still
  extracting), not afterglow harvest.

### Measured (C) — SPECKLE's co-pump at a third timestep: the coarse
### positive does not survive refinement

- New pair cp27k_dt12 (Nathaniel's model at the lasing point, 30 kV ×
  r_out 0.9, window 100–170 ns, 27 kW, dt 12.5 ps × 24000 steps):
  laser marginal **−5.8%**, fluorescence −0.7%
  (SPECKLE/runs/cp27k_dt12_copump_report.json). The dt series is now
  +16.6% (50 ps) → −24.8% (25 ps) → **−5.8% (12.5 ps)**: the two
  finest steps agree in SIGN (negative) and move toward zero, so the
  coarse-step +16.6% was a discretization artifact; the magnitude is
  still unconverged, consistent with a small-or-negative true
  marginal. Within the independent model the plateau co-pump has no
  resolvable benefit; its afterglow marginal (+0.4–1.1%,
  dt-robust) remains that framework's one laser-positive — the same
  high-margin-afterglow class as (B).

### Verified

- All three parts measured live 2026-08-27 from the MAIN checkout
  under /usr/bin/python3, on the post-[0.1.27] engine — that entry's
  changes (swarm-residual sink on the E/N-direct path, circuit filter
  branch, ArF-Chen low-E/N table) are structurally inert on these
  beam-drive, Te-keyed configurations, and part A's exact reproduction
  of the recorded 0.1.23 bracket endpoints and nominal (−1.65 / +8.00 /
  −5.62 vs −1.66 / +7.9 / −5.62) is the in-run check of that. Part B
  reuses mw_benchmarks' own case builder and segmented integrator (the
  0.1.21 cure) and copump_verify's bracket machinery unchanged.
- Version note: numbered 0.1.28 — [0.1.27] was taken by the parallel
  session's Te-closure entry while this work was in flight (the same
  race [0.1.26]/[0.1.27] already record).

## [0.1.27] — 2026-08-26 — the Te closure, the paper's filter branch, and the low-E/N table

Trigger: "begin with the remaining bugs/issues in SPECK" — the three
actionable residuals on the 0.1.25 comparison audit's revised list: the
unattributed x1.55 Te anchor row, the RLC filter branch stated as
unmodeled, and the 70%-below-table breakdown ramp.

### Fixed — the Te x1.55 row: attributed, then closed (the swarm-residual sink)

- **Attribution (measured 2026-08-26).** The LoKI table's own swarm point
  at the anchor's E/N peak is <eps>(14.1 Td) = 6.91 eV -> Te 4.61 eV —
  matching Chen's Fig. 10c (~4.7 eV). The engine's fluid Te sat at
  7.73 eV: x1.67 above ITS OWN tier's equilibrium — an internal
  energy-balance defect, not swarm physics. Mechanism, measured from the
  data: on the E/N-direct path the bucketed inelastic cooling is
  Te-FLAT, and the row inventory books only 88% of e v_d E at 14.1 Td
  (level-lumped bucket thresholds; the unbooked residual runs +5% at
  5 Td -> +45% at 150 Td across the table). The unbooked share has no
  Te-responsive absorber — Ne/Ar elastic alone would need Te ~ 34 eV to
  soak 12% — so the fluid floats until dilution and the stepwise
  channels absorb it.
- **Fix: the swarm-residual sink** (`chemistry.rhs`, E/N-direct path
  only): p_res = max(e n_e v_d(E/N) E − Σ table-row threshold drains −
  elastic@Te_sw − vib − <eps>·(table ionization births), 0), routed
  e -> gas like the vib sink. DERIVED from the same table (v_d, <eps>)
  plus the network's own row inventory — nothing fitted. With it,
  dTe/dt closes EXACTLY at Te = (2/3)<eps>(E/N): probe |dTe/dt| ~
  1e-6 eV/s at te_sw against ±2.6e8 restoring at ±30%, three E/N
  points; the same point read +3e9 eV/s pre-fix. Ionization growth is
  credited as dilution — births, not net: attachment's 1.5 Te capture
  drain cancels in the Te equation (derivation in the sink's comment).
  Gate: `test_swarm_residual_sink_closes_te_at_the_table_point`
  (teeth margin ~1e9).
- **The residual undershoot is attributed too**: the run now settles
  Te 3.4–3.5 eV at the E/N peak (x0.69 of Chen's ~5, from x1.55 OVER).
  At the live 80-ns state dTe/dt = 0 at the settled Te while at
  te_sw 4.61 the −9.0e8 eV/s imbalance equals EXACTLY the dilution from
  beyond-table electron births (nu_total 5.6e8 /s vs table ionization
  3.6e8 — stepwise/Penning are a third of the ionization at breakdown,
  born at 0 fluid energy, Penning's exothermic release correctly
  credited to u_to_e). Legitimate beyond-swarm dynamics, in Chen's own
  post-breakdown 3–4 eV band — recorded, not "fixed".

### Added — the paper's RLC filter branch (circuit 5 -> 7 states)

- `circuit.py`: if/qf states — the Sec II A filter across the gap
  (Rf 70.71 ohm / Lf 100 nH / Cf 10 pF, "to suppress circuit
  oscillations after gap breakdown"); node continuity now
  dE/dt = (i2 − i_cond − i_f)/(eps0 A); `energy_J` carries the Lf/Cf
  stores (Rf dissipation is a circuit-side loss, like R). Params
  rf_ohm/lf_nH/cf_pF (declared, paper defaults); the engine passes them
  and extends the aux atol. This closes 0.1.25's stated simplification
  (b); the remaining one is the fixed Lm = 250 nH.
- Gates: the LC-conservation gate re-pinned at rf = 0 (the branch is a
  lossless LC there — conserved to the same 1e-3); new
  `test_filter_branch_damps_and_only_dissipates` (energy never
  increases within 1e-6 E0; measured 2e-5 of the bank dissipated in the
  no-breakdown case — quasi-static grazing, i_f ~ Cf dV/dt ~ 0.1 A, as
  expected: the branch's real work is post-breakdown ringing).
- Measured on the anchor: totals unchanged at the % level — the branch
  shapes the post-breakdown ringing, exactly the paper's stated role.

### Added — ChenArF_NeArF2_lowEN (the breakdown ramp leaves the stand-in)

- LoKI-MC run 2026-08-26 with spectre's own curated toolchain
  (`Code/lokimc`, arm64; 43.6 min): the ChenArF deck byte-identical
  except reducedElecField [0.5, 1, 2, 3] Td and the output folder — the
  V2point lowEN convention (2e4 electrons/point, boltzmannMC,
  equalSharing). Deck at spectre `Input/Mixture/ChenArF_NeArF2_lowEN.in`;
  job dirs stay with spectre; `lookUpTableRateCoeff.txt` +
  `lookUpTableSwarm.txt` copied byte-verbatim (cmp-checked) to
  `Curated/Output/ChenArF_NeArF2_lowEN/` and to
  `MODULES/loki_tables/ChenArF_NeArF2_lowEN/` (README updated).
  `RUN_NAMES_ARF_CHEN` merges both; the 0.1.25 overlap guard verifies
  grid disjointness at load.
- Physically sensible and smoothly connecting: <eps> 2.26/3.17/4.63/
  5.27 eV at 0.5/1/2/3 Td toward the 5-Td 5.83; BulkV 7.1e3 -> 1.3e4
  m/s; ionization ~0 below 2 Td; F2 attachment RISING toward low E/N
  (1.06e-18 m3/s per-N at 0.5 Td) — the seed-erosion physics the
  stand-in missed.
- Anchor coverage flips 28%-in/72%-below -> **74%-in/26%-below**: the
  breakdown ramp is now on measured rates; the remainder is the deep
  afterglow below 0.5 Td, where every rate is genuinely ~0. Breakdown
  sharpens: E/N peak 14.1 -> 14.6 Td (x0.97 of the 15-Td class,
  toward Chen's ~16), FWHM 19 -> 22 ns.

### Measured — the anchor after all three (live 2026-08-26)

deposits 1.121 J of the 1.62 J bank (x0.692) | n_e 4.49e20 (x0.69) |
Te at E/N peak 3.46 eV (x0.69, was x1.55 OVER) | n_ph 1.97e18 (x0.020) |
P_L 0.092 MW (x0.018) | FWHM 22 ns (x1.47) | E/N pk 14.6 Td (x0.97) |
E_L 2.04 mJ | ledger closes to 5.4e-6. The photon-side x0.02 residual
stands as the honest lane (Chen's printed formation rates by scoring
policy — CHANGELOG [0.1.25]); Te and the ramp are no longer on the
suspect list.

### Verified

- Full suite under `/usr/bin/python3`: **197 passed** (5:27, final run
  with all three changes; 195 + the two new gates; the base-case test's
  5-Td floor pin updated to the new 0.5-Td floor). SPECK/tools coupler
  parity 3 passed (its E/N-supplied path carries the sink on both sides
  of the gate). Homebrew 3.14/numpy 2.5: the anchor preset and every
  probe in this entry ran clean under `/opt/homebrew/bin/python3`.
- Version note: numbered 0.1.27 because a parallel session's
  [0.1.26] (microwave reheat across the benchmarks) landed in the main
  checkout while this work was in flight.

## [0.1.26] — 2026-08-26 — microwave reheat across the benchmarks

Owner ask: "test the reheating on the different benchmarks to see if it
increases efficiency for all of the cases" (+ the frequency BOTE, the
no-go audit, ArF, and the COMSOL decision brief that closes the arc).
Tools built and first measured 2026-08-21; delivery queued behind
[0.1.25], whose CHEN_ERRATA touch Ne-partnered three-body rows (real
flux in Ne fills) — so the LOKI table below was RE-MEASURED 2026-08-26
on the post-errata engine from the main checkout. Ar-fill rows came
back unchanged at table precision; the Ne-optimized row moved (below).
Kinetic-tier rows were measured 2026-08-21/24 pre-errata (5 of 8; both
quoted sign flips are Ar-fill cases, where the errata rows carry ~0
flux).

### Added

- `tools/mw_benchmarks.py`: every beam-drive benchmark configuration the
  package carries, run as twin pairs for the two microwave windows
  (afterglow reheat from beam-off + 10 ns, 150 ns flat; plateau co-pump =
  second half of the beam pulse ending at beam-off) at ONE stated
  coupling rule — the DE-AR booster's 250 kW per litre of cell
  (2.5e8 W m⁻³ × cell volume), 95 GHz, field ceiling 30 Td; the ARPA
  chamber keeps its X-band hardware. Cases: krf_beam_0d (DE-AR class,
  3 atm, Electra-class optics), krf_ne_0d (Ne fill), the Kannari 1985
  anchor cell (R_OC 0.2; j iterated to the recorded 66 J deposit, the
  scorecard's own construction), the NRL Electra ArF anchor at its
  1.09 atm peak point (arf_beam, 193 nm; j iterated to 0.8 MW/cc), the
  HYPER cell, the ARPA column, and the 0.1.24 Ar- and Ne-optimized
  points. Per case/scheme: E_MW, dE_laser, marginal, eta without / with
  microwave (microwave in the denominator), the change in percentage
  points, post-beam fluorescence change. `--tier kinetic` re-runs the
  sweep on the PIC_MCC_3D-measured rates (the 0.1.23 sign-sensitivity
  column); `--cases` selects a subset.

### Added (2) — the frequency question, closed by BOTE

- `tools/ib_feasibility.py`: inverse-bremsstrahlung coupling vs frequency
  (1 GHz → CO₂ 28.3 THz) with ν_m from the engine's OWN swarm transport
  (0.9–8.9×10¹² s⁻¹ over 0.5–30 Td at 2 atm → collisional/IB crossover
  ~0.3 THz). MEASURED CONCLUSIONS: below the crossover the AC factor is
  ~1 AND the absorption length is frequency-independent
  (α ≈ ω_p²/(c ν_m); 1.8 mm at n_e 10²⁰ m⁻³ for EVERY band 1–300 GHz) —
  the plasma cannot tell 2.45 GHz from 95 GHz, so the reheat failure is
  NOT a frequency/coupling failure; it is kinetics + the wall-plug
  inequality η_source × marginal > η_beam(~7%). That inequality is met
  at best MARGINALLY by ~70%-efficient magnetrons (1–2.45 GHz) times the
  fragile +12% plateau co-pump marginal (→ ~8.4%), with the absolute
  gain still capped at P_MW/P_beam ~ 1e-3; gyrotrons (35–40%) and
  everything above the crossover lose outright — the DE-AR's own 95 GHz
  gyrotron reaches η×m ≈ 4.2% < 7%, and the CO₂ laser (ω/ν ≈ 90,
  coupling ∝ ω⁻², source ~12%) reaches ≈ 1.4%: CO₂ reheat cannot raise
  wall-plug efficiency (the owner's hunch, quantified). The one
  EM-positive channel (afterglow reheat of a HIGH-MARGIN cell still
  lasing through the window — Ne-opt +10.0% LOKI / Ar-opt +9.5%
  KINETIC, dE_laser ≫ dFluor there) is the same feed-an-extracting-
  cavity mechanism with the same fragile-sign, wall-plug-capped
  arithmetic. rebuild_out/mw_benchmarks/ib_feasibility.json.

### Added (3) — the no-go audit ("is there any case where reheating cold
### electrons works?")

- `tools/afterglow_nogo_audit.py`: every KrF*-manifold FORMATION row in
  the network, with its numerical response to heating Te 0.5 → 3 eV. All
  22 direct formation rows are heavy-particle; the Te dependence lives
  one step upstream and BOTH supply lines fall when heated: F⁻ (fuel of
  every ion–ion row) via attachment ×0.13, Kr*/Kr** (fuel of every
  harpoon row) via dissociative recombination ×0.37. The one
  heat-positive route (direct e-impact excitation, k(Kr*)
  2.8e-23 → 1.8e-15 m³/s from 0.5 eV to the 30 Td swarm point) needs
  E/N ≳ 20–30 Td while the conductive afterglow pins the settled field
  ≤ 2 Td at any power (measured to 1 GW) — the two walls close every
  route. The measured "+9–10%" exceptions are cells still LASING through
  the window (late co-pump), not afterglow harvest. Contrast stated:
  electron heating pays where the upper level sits at eV energies
  (CO₂/CO vibrational sustainers) — the excimer chain is the wrong
  class. `--arf` runs the same audit on the ArF (193 nm) chain via the
  Electra anchor configuration: identical verdict — all 6 direct ArF*
  formation rows heavy-particle (Ar*/Ar** + F₂ harpoon, Ar₂* rows,
  Ar⁺/Ar₂⁺ + F⁻ ion–ion), F⁻ ×0.13 and Ar* (Ar₂⁺ DR) ×0.37 supply
  falls, k(Ar*) needs the same forbidden 20–30 Td; if anything worse
  positioned than KrF (Ar* threshold 11.6 vs 9.9 eV, no heavier-excimer
  transfer feeds ArF*). Matches the sweep's empirical Electra row
  (afterglow −1.3%, co-pump +5.7% < its own 6.13% baseline).
  rebuild_out/mw_benchmarks/afterglow_nogo_audit{,_arf}.json.

### Added (4) — the COMSOL decision brief

- `Kinetic_Modeling/docs/comsol_case.txt`: the evidence pack for the
  license question ("can we present, with numerical evidence, that a
  COMSOL purchase is warranted to investigate operating parameters and
  horn shape for an efficiency increase?"). The honest answer the
  numbers force: the sweep does NOT support horn-shape-for-efficiency
  (every model already granted the microwave perfect uniform coupling —
  the ceiling a horn can only approach from below — and |Δη| stayed
  ≤ 0.006 pp), but it DOES support the purchase on the DE-AR modeling
  milestones themselves: the award requires EM power deposition into
  the gas cell, <10 ns spatiotemporal deposition maps, and sim↔exp
  deposition agreement within 10%, none of which any in-stack tool
  (0-D/1-D/uniform-field) can compute — while the measured absorption
  length swings ~5 cm → ~2 mm across one pulse (optically-thin →
  skin-dominated, ×25 through the cell size), exactly the regime where
  only a full-wave solver knows where the power lands. Also states the
  one place horn design touches efficiency (the fragile plateau
  co-pump window: power inside the gain volume only, during the 60 ns
  plateau only) and what measurement would promote it.

### Measured — no benchmark gains more than 0.006 percentage points

LOKI tier, 250 kW/L rule, post-[0.1.25] engine
(rebuild_out/mw_benchmarks/mw_benchmarks_loki.json; segmented
integration through the afterglow — the Ne fill empties fastest and
NaN'd the BDF Jacobian in a single solve, the 0.1.21 class):

| case | eta0 | afterglow: marginal / Δeta | co-pump: marginal / Δeta |
|---|---|---|---|
| krf_beam_0d (DE-AR class, 3 atm) | 0.00% — does not lase (the recorded sub-threshold v2 point) | — | — |
| krf_ne_0d (Ne fill, same cell) | 0.00% — does not lase | — | — |
| krf_hyper_mw (DE-AR HYPER) | 1.54% | −2.8% / −0.0055 pp | +12.5% / +0.0029 pp |
| krf_arpa_mw (X-band 1 MW) | 2.16% | −10.4% / −0.0036 pp | −1.9% / −0.0004 pp |
| ar_optimized (0.1.24) | 4.95% | −2.1% / −0.0056 pp | +1.3% / −0.0010 pp |
| ne_optimized (0.1.24) | 6.19% | **+10.0% / +0.0024 pp** | +5.1% / −0.0002 pp |
| kannari_R20 (anchor, 66 J) | 0.70% | −0.0% / −0.0003 pp | +3.0% / +0.0002 pp |
| electra_arf_1p09 (ArF 193 nm anchor) | 6.13% | −1.3% / −0.0031 pp | +5.7% / −0.0001 pp |

- **The [0.1.25] CHEN_ERRATA moved exactly the rows they should**: every
  Ar-fill row reproduced the 2026-08-21 pre-errata measurement at table
  precision (largest drift: HYPER co-pump marginal +12.4 → +12.5%, the
  beam-split F₂ booking fix), while the Ne-optimized cell — where the
  corrected Ne-partnered three-body rows carry real flux — rose
  η₀ 5.76 → 6.19% (+0.43 pp; the 0.1.24 ladder's top rung, now the
  post-errata headline) with afterglow marginal +9.3 → +10.0%. The
  pre-errata table stands in this entry's 2026-08-21 draft history and
  the run logs.
- **Afterglow reheat raises the intrinsic efficiency in 1 of 8 cases**
  (the Ne-optimized cell) and lowers it in the other five that lase;
  the post-beam fluorescence rises in every case (the DE-AR observable
  again).
- **Plateau co-pump has a positive marginal in 5 of 8** (every lasing
  case but the ARPA column) but raises the BLENDED efficiency in only 2
  (HYPER +0.0029 pp, Kannari +0.0002 pp — the marginal must exceed the
  cell's own eta for the microwave joules in the denominator to pay:
  Electra ArF +5.7% marginal vs 6.13% baseline, Ne-optimized +5.1% vs
  6.19% → both DOWN).
- The largest |Δeta| on the whole sweep is **0.0056 percentage points**.
  Across the DE-AR-class cell, the Ne fill, the Kannari and Electra-ArF
  anchors, the HYPER and ARPA cells and the optimized points, neither
  window moves the efficiency by more than the 0.1.24 lever scan's
  smallest step. The answer to "does it increase efficiency for all the
  cases" is: it does not materially change it in ANY case at feasible
  coupling — the efficiency lives in the beam-only levers.
- Both DE-AR-class 0-D benchmarks (Ar and Ne fills at 3 atm,
  Electra-class optics) do not lase at all in this engine — the 0.1.9
  sub-threshold verdict at the v2 point holds for the Ne fill too; the
  cell needs the 0.1.24 levers (F₂, pressure, optics) before any
  microwave question is meaningful there.
- **Kinetic tier (PIC_MCC_3D-measured rates), 5 of 8 measured
  2026-08-21/24 pre-errata** (both flips are Ar-fill cases — errata rows
  ~0 flux there): the sign story of 0.1.23 generalizes. HYPER co-pump
  **+12.5% (LOKI) → −2.7% (KINETIC)**; Ar-optimized afterglow
  **−2.1% (LOKI) → +9.5% (KINETIC)**; ARPA afterglow −10.4 → −3.9%,
  co-pump −1.9 → −3.4%; the two DE-AR-class cells stay non-lasing on
  both tiers; largest kinetic |Δeta| 0.0036 pp. No case has a
  tier-stable positive sign, and no case moves materially on either
  tier — the magnitude verdict is tier-robust, the sign is not. The
  remaining three kinetic cases (ne_optimized, kannari_R20,
  electra_arf_1p09) are still integrating in two nice'd background
  processes (the 2026-08-21 full sweep, pre-errata; a 2026-08-26
  3-case tail, post-errata — day-scale: the kinetic tier's
  relaxation-bridged sub-table curves make the segmented afterglow
  brutally stiff); their JSONs land in rebuild_out/mw_benchmarks/
  whenever they finish, and this table's "5 of 8" is the honest
  current record.

### Verified

- The LOKI sweep was re-run 2026-08-26 from the MAIN checkout on the
  post-[0.1.25] engine (this entry's table) — the tools import and run
  green under the errata'd network, and the before/after diff behaved
  exactly as the errata predicts (Ar rows fixed, Ne rows move; bullet
  above). ib_feasibility and both no-go audits were regenerated from
  the main checkout 2026-08-25 (JSON mtimes) and are errata-unaffected
  (swarm-transport ν_m; Ar-fill audit states).
- Suites: [0.1.25] left REBUILD at 195 green + PIC_MCC_3D 211 +
  SPECK/tools 3; this entry adds documentation and tools only (no
  MODULES/ change), and the three tools were exercised end-to-end by
  the re-runs above. Kinetic-tier completion is PENDING as stated —
  the 5/8 rows are measurements, not projections, and are labeled with
  their engine state.

## [0.1.25] — 2026-08-25 — the Chen lane's transcription defects, measured and fixed

Trigger: the 1DModelKrF-vs-SPECK comparison audit (2026-08-25) re-read
Chen 2025 first-hand, re-ran this anchor live, and attributed the recorded
"x3-5 per-joule formation gap" by controlled swaps through the scorecard's
own rebuild mechanism. Verdict on our books: the gap was mostly OURS — a
circuit transcription defect plus four unflagged impossible printed rates
— not the network physics the 0.1.5 verdict suspected. This entry lands
the fixes; every number below was re-measured after them.

### Fixed

- **Chen circuit constants were the Fig. 5 validation experiment's, not
  the base case's** (`params.py`; first-order in the anchor miss).
  Symptom: anchor n_e x0.17-x0.27, no self-consistent lasing, "1.68 J
  stored" vs the paper's stated E_in = 1.62 J. Root cause: 10.4 nF /
  8.1 nF / 30 nH and the "5 -> 1 uH magnetic switch" are Sec. II D's
  VALIDATION rig (2 atm, 20 mm, -13 kV, mirror-less), grafted onto
  Table II conditions; the base case all results figures use is Sec.
  II A: C0 = 10 nF, Lm = 250 nH, Cp = 10 nF, Ls = 5 nH, V0 = -18 kV
  (re-read first-hand 2026-08-25). Defaults corrected; circuit.py
  docstring now states the two remaining simplifications (fixed Lm at
  250 nH; the paper's RLC filter branch — 70.71 ohm/100 nH/10 pF,
  post-breakdown ringing suppression — is not modeled). New Known-bad.
  Verified: `test_chen_default_circuit_is_the_paper_base_case` pins the
  constants and E_in = 1.62 J exactly.
- **Four physically impossible printed rates ran unflagged in the CHEN
  lane** (`chen2025.py`, `network.py`; the other first-order term).
  R13/R14 print 2.08e-8/5.61e-9 "m3/s" (cm3 misprint, x1e4 above the
  electron gas-kinetic ceiling); R23 prints their own R27 value pasted
  onto the excimer row (x22.7 the Ar*+Ar+M literature rate); R24 prints
  4e-40 m6/s (x1e6; as printed it annihilates Ne* in ~0.4 ps and severs
  the 96%-Ne bath's Penning chain). R36 (x100 low, flagged 0.1.4-era but
  kept verbatim) joins them. All five adjudicated against primary
  sources into `CHEN_ERRATA` (printed value + corrected value + basis
  per row); `chen_rows(errata=True)` is the default lane — the errata
  lane reproduces Chen's MODEL (the printed values cannot describe what
  he ran), the as-printed lane (`chen_errata=False`, declared parameter)
  remains the record of the PAPER. TABLE_III itself is untouched.
  In krf/arf_beam the corrections also land on the Ne-partnered
  three-body rows (R23/R24/R36) — flux ~0 in Ar fills (anchors
  unmoved), real physics in Ne fills. New Known-bad. Verified:
  `test_errata_pins_and_table_stays_as_printed`,
  `test_chen_errata_lane_default_and_as_printed`.
- **ASE window irradiance ~1e7 low — chords marched OUT of the box**
  (`ase.py window_irradiance`). Symptom: 3.37 W/m^2 where the
  uniform-field value is ~5e7. Root cause: the hemisphere filter keeps
  arrival-side directions (d_z < 0) but the march went along -d — from
  a point 1 nm inside the +z wall straight out — so every chord was
  ~1 nm long. The suite's gate only asserted positivity, so it never
  caught it. Fixed (march along d); gate replaced by a quantitative
  closed-form hemisphere pin (`test_window_quantitative_uniform_gate`,
  agreement 1e-9). Diagnostic-only — nothing fed back into the RHS.
- **ASE ray transfer was exact only for uniform fields — segment order
  reversed** (`ase.py ray_intensity`). The outward-from-observer march
  amplified each segment's source by the cells BEHIND it rather than
  the cells between it and the observer; identical for uniform media
  (why the closed-form oracle always passed), ~3% on the 4-pi scalar
  flux and up to 12% on single rays for gain-core/absorbing-skin
  fields. Fixed by applying the recurrence over `reversed(segs)` — the
  transfer is now exact for ANY voxel field; the two-segment hand
  gate re-pins the arrival order and proves the outward order differs
  (`test_two_segment_hand_computation_and_order_physics`).
- **`electron_tier="kinetic"` + circuit drive crashed on the first RHS
  call** (`engine.py`). KineticRates exports k(<eps>) only (no
  E/N-indexed curves, no drift velocity in the npz) while the circuit
  path looks rates and conduction up by E/N — AttributeError at
  `chemistry._tier_k`/`drift_velocity`. Now refused at Device
  construction by CAPABILITY (k_by_en/drift_velocity/en_range_td), with
  the fix stated (`test_kinetic_tier_circuit_drive_refused_at_construction`).
- **Non-finite blow-ups surfaced as an inscrutable `lu_factor`
  LinAlgError** (`engine.py rhs`). One isfinite pass per RHS call now
  raises FloatingPointError naming t and the offending states. Found
  running the as-printed table under the paper's true circuit — that
  lane is measured NOT to integrate (blow-up at ~138 ns in U_e/T_g;
  recorded in the scorecard as a dated measurement, deliberately not a
  live re-run — it grinds the stiff BDF for minutes). The guard
  MECHANISM is suite-pinned cheaply
  (`test_nonfinite_rhs_guard_names_the_state`).
- **Beam split parked F2(C1/H1) excitation as never-decaying F2(A1)**
  (`beam.py spencer_fano_yields`). Only C1/H1 (11.57/13.08 eV) clear
  the noble prompt floor, and they were booked as bound F2(A1) — a
  species the Kannari-chemistry modes never decay (~0.005% of drive at
  0.1% F2, growing with F2 fraction). F2 electronic excitation above
  the floor now books 2 F + the excess over D0 as fragment heat,
  mirroring `reroute_repulsive_f2`; identity exact by construction
  (extended `test_prompt_floor_banks_subexcitation`).
- **The circuit validity gauge was one-sided** (`engine.py _collect`).
  `en_above_table_frac` never counted BELOW-floor excursions — the
  whole breakdown ramp (<5 Td in the Chen set, which ships no low-E/N
  companion) rode the flat+taper stand-in while the scorecard said
  "E/N 100% in-table". Added `en_below_table_frac` +
  `en_table_range_td`; the Chen anchor reports it (72% of the 200 ns
  window at the corrected constants).
- **Jacobian sparsity ignored the MW/ASE global couplings**
  (`engine.py jac_sparsity`). The volume-averaged in-loop ASE field and
  the volume-summed MW power balance couple every cell to every cell;
  the block-diagonal pattern starved the BDF Newton matrix of those
  entries (convergence quality, not answer correctness). Both couplings
  now added when active on a multi-cell grid.
- **Latent LoKI cross-directory double-count** (`loki_tier.py`). Two
  run directories tabulating the SAME E/N would silently double the
  bucket rates (the merge sums). Shipped grids are disjoint, so it
  never fired; a guard now raises on overlap, naming both directories.

### Changed

- `scorecard.chen_figure6_rows`: targets re-read first-hand — deposited
  energy vs the PAPER-STATED E_in = 1.62 J (was "1.68 J", the wrong
  capacitor); n_e 6.5e20 from Fig. 10a's LINEAR axis at 0.1% F2 (the
  old 1e21 was a coarse log-axis read of Fig. 6b — a x1.5 target error
  booked against the engine); added T_e-at-breakdown (Fig. 10c class)
  and pulse-FWHM rows; the as-printed lane is re-run and RECORDED; the
  verdict now carries the measured attribution instead of the 0.1.5
  suspect list. `run_scorecard`'s Chen block and the OMISSIONS circuit
  entry rewritten the same way (the old claims are kept as HISTORY,
  marked measured-wrong).

### Measured — the anchor, before and after (all live 2026-08-25)

Shipped state (Fig-5 circuit + verbatim rates): n_e 1.73e20 (x0.27 of
Fig. 10a's 6.5e20), n_ph 2.5e15, P_L ~0 — no lasing. Corrected lane
(Sec II A circuit + CHEN_ERRATA): deposits 1.056 J of the 1.62 J bank,
n_e 5.11e20 (x0.79), T_e at breakdown 7.7 eV (vs ~5), E/N peak 14.1 Td,
lasing IGNITES — P_L 0.095 MW, FWHM 19 ns, E_L 2.08 mJ (vs 75 mJ:
x0.03; ledger closes to 2e-5). As-printed rates + true circuit: does
not integrate (the two shipped defects masked each other — the softer
wrong circuit was what let the impossible table run). The residual
photon-side gap vs Chen's published output is REAL and now honestly
stated: Chen's printed formation rates (R29/R30/R31) sit below the
measured literature values (Velazco/Mezyk-Flannery class); adopting
those is a beyond-Chen physics decision deliberately NOT taken in the
scoring lane — the lane scores THEIR model.

Verification: full suites green post-fix — REBUILD 195 (incl. the new
gates: circuit-constant pin, errata pins both lanes, quantitative ASE
window closed-form, arrival-order hand arithmetic, non-finite guard,
kinetic+circuit refusal, F2 prompt-dissociation identity, below-table
gauge), PIC_MCC_3D 211, SPECK/tools 3. ASE fixes additionally checked
against independent closed forms out-of-suite (window/scalar/2-cell:
agreement 1e-9 or better). Chen anchor re-run live (numbers above);
status.html, run.py preset text, and the scorecard HISTORY block
updated the same session.

## [0.1.24] — 2026-08-21 — where the efficiency goes, and the levers that raise it

Owner ask: "figure out a way to increase the efficiency or see if it is
just an ionization issue."

### Added

- `CompiledNetwork.row_fluxes(dens, te, tg, ne, ntot, en_td)` +
  `row_stoichiometry()`: per-row reaction flux at a state — the SAME
  k × heavy densities × n_e^pow × T_g/screening arithmetic the rhs loop
  applies, exposed for channel budgets. Gated
  (`tests/test_row_fluxes.py`): Σ(coeff × row flux) reproduces the rhs
  ddens to 1e-10 at n_ph = 0 / ase 0 on a live beam-on state. No
  behavior change.
- `tools/efficiency_budget.py`: the KrF energy/channel budget of one run
  — beam partition (Spencer–Fano yields: ionization / excitation /
  electron heat per pair), KrF* MANIFOLD formation by class (ion–ion,
  neutral harpoon, ArF*/Ar₂F*/Kr₂F* transfer) integrated from the
  per-row fluxes, KrF* loss by class (stimulated, spontaneous, 3-body,
  quench by partner), the fate of the stimulated photons (extracted /
  mirror-window / photo-absorbed BY SPECIES from the cavity absorber
  inventory), ion-pair recovery, eV per KrF*, formation and extraction
  efficiencies, ledger sinks. First cut mis-counted B↔C collisional
  mixing (CB_Ar/BC_Ar at 1e20 each way) as formation/loss — the manifold
  (B, C, v levels) is one reservoir; fixed before any number was used.
  `--mixture ne`, `--f2`, `--set k=v` overrides.
- `tools/efficiency_levers.py`: one-at-a-time lever scan around the
  HYPER reference (pressure, pump density, pulse length, gain length,
  R_OC, L_op, F₂ and Kr fractions — table caveats stated) + the greedy
  combination of the winners.

### Measured — "is it an ionization issue?" NO; here is what it is

Reference (krf_hyper_mw, 2 atm Ar/Kr/F₂ 0.949/0.05/0.001, 100 A cm⁻²,
60 ns, R_OC 0.95 / L_op 0.02): 31.6 J deposited → 0.485 J (1.54%).
- Beam partition: W = 30.4 eV/pair — 51% ionization, 37% excitation,
  12% electron heat. **Formation is healthy: 20.8% of deposited energy
  becomes KrF* quanta (24 eV per KrF*)** — 50% neutral harpoon, 29%
  ion–ion recombination, 18% ArF*/Ar₂F*/Kr₂F* transfer; **37% of all ion
  pairs return as KrF*.** Ionization is the beam's cost, not the loss.
- The losses are downstream: **53% of KrF* are killed by collisions with
  the Ar buffer** (39% three-body → ArKrF*/Kr₂F* trimers, 14% Ar quench);
  13% spontaneous; only 23% emit stimulated photons, and only a third of
  THOSE leave the cavity (the rest die in intracavity absorbers +
  HR/window) → **7.4% of KrF* become laser photons**. Ledger: 75% of
  the deposited energy ends as gas heat, 14% radiated, 8% chemical
  (F atoms), 1.5% laser.
- Levers (one at a time): **F₂ fraction dominates** — 0.05% no lasing,
  0.1% 1.54%, 0.2% 2.77%, 0.3% 2.77% (the reference fill is
  F₂-starved: burn 21% → 13%); **pressure** 1 atm 2.59% (three-body
  39% → 20%, extraction 7.4% → 11.9%), 1.5 atm 2.30%, 3 atm 0.12%;
  **L_op** 0.005 → 2.06%, 0.05 → 0.86%; R_OC optimum AT 0.95 (0.90
  1.09%, 0.99 0.72%); pump 200 A cm⁻² 1.92% (400 → 1.16%, F₂ burn 66%);
  pulse 100 ns 1.84%; cell 0.2 m 1.59%; Kr 5% optimal (2.5% 1.33%, 10%
  1.19%; table caveat).
- **Greedy combination (1 atm, 200 A cm⁻², 100 ns, 0.2 m, L_op 0.005,
  F₂ 0.2%): η = 4.95%** (5.05 J / 102 J; margin g₀/g_th 19). Budget
  there: three-body 7%, stimulated 64% of KrF* — but **64% of the
  stimulated photons are photo-absorbed: Ar* 51%, F₂ 17%, ArKr*/F⁻/Kr₂⁺
  ~7% each**; 21% of KrF* extracted.
- **Ne buffer** (Kannari-class; the Ne-fill LOKI tables): removes the
  Ar* absorber and the ArF* transfer channel alike — un-optimized 2 atm
  2.63%, 1 atm 2.25%; **optimized Ne at 2 atm → η = 5.76%** (7.23 J /
  125 J, ×3.7 the reference; ion-pair recovery 68%; 42% of stimulated
  photons extracted; the residual absorber is **F₂ itself, 60% of what
  is absorbed** — the harpoon fuel is the 248 nm absorber, the trade-off
  that sets the F₂ optimum); 3 atm Ne 4.06%.
- The ladder in the model: Ar 2 atm 1.54% → Ar 1 atm 2.59% → Ne 2 atm
  2.63% → Ar optimized 4.95% → **Ne optimized 5.76%**. Formation runs at
  17–23% of deposited energy everywhere (the KrF quantum ceiling class:
  5 eV photons from ~21–29 eV per KrF*); every gain above comes from
  LOSING FEWER KrF* (buffer quench) and ABSORBING FEWER PHOTONS (Ar*,
  F₂, F⁻, Kr₂F*, Kr₂⁺). The gap to Kannari-class 10–12% intrinsic is
  the engine's recorded formation deficit plus this small-cell optics —
  relative slopes are the robust output.
- Not scannable honestly: mixed Ar/Ne buffers and the F₂ optimum in the
  Ne fill beyond ×2 (new LoKI-MC tables; LoKI-MC is in-repo).

### Verified

- **190 passed** under `/usr/bin/python3` (135 s; 178 + 11 MW/timing
  gates + `test_row_fluxes_reproduce_rhs`), re-run from the main
  checkout after delivery (see below). Homebrew 3.14.6 smoke with
  `-W error::DeprecationWarning` green on `row_fluxes` (169 rows on a
  live beam-on state). The budget's formed/lost manifold accounting
  closes to +0.00% on every run (formed = lost, the stored KrF* at the
  250–400 ns horizons being nil); the stimulated-photon fate sums to
  100% (extracted + HR/window + photo-absorbed) against the ledger's
  laser_out / cavity_loss channels.

## [0.1.23] — 2026-08-21 — is the co-pump gain numerically real? (the verification)

Owner ask: "Is there a way to know that the positive increase of
efficiency is numerically correct? Is there a way to verify statistically
to ensure no false negatives or false positives? If this is working, can
we go forward with testing it using SPECK?"

### Added

- `tools/copump_verify.py`: the verification ladder on the reference
  co-pump point (2 atm, 100 A cm^-2, 250 kW on the 35–74 ns plateau):
  (1) numerical convergence of the twin difference dE_laser across rtol
  1e-6/1e-8/1e-10, max_step 2/1 ns, atol ×0.1, and an INDEPENDENT
  integrator (Radau vs the engine's BDF); (2) the null floor — the spread
  of dE across settings is the differential noise (the twin structure
  cancels common-mode error), z = dE/σ_diff is the false-positive
  statistic, 3σ_diff the false-negative bound (smallest resolvable dE →
  P_min); (3) physics robustness — the marginal under every active EST
  row ×1/3 / ×3 and every LIT-conflict alternate (scorecard's S5 bracket
  machinery, both twins on the same variant network) and on the
  INDEPENDENT rate tier (KINETIC = PIC_MCC_3D-measured vs LOKI); (4) the
  experimental power analysis (interleaved on/off shot pairs, α 0.05 /
  power 0.8, N = (z_a+z_b)² 2σ²/δ²).
- Engine: the window's swarm TRANSPORT (v_d(E/N), <eps>(E/N), E/N range)
  now comes from `Device._mw_swarm` — the LOKI mixture tables whatever
  tier keys the rates (the KINETIC tier exports rates only; transport is
  an E/N-native quantity) — so the same window runs on either tier (the
  cross-tier check; gate `test_mw_window_runs_on_kinetic_tier`, 11 gates).

### Measured — numerically real, physically FRAGILE (the honest answer)

Reference point 2 atm / 100 A cm^-2 / 250 kW on the 35–74 ns plateau;
rebuild_out/copump_scan/copump_verify.json.

- **Numerics PASS.** dE_laser = 0.7308 ± 0.0002 mJ across seven solver
  settings (BDF rtol 1e-6/1e-8/1e-10, max_step 1 ns, atol ×0.1, and the
  independent Radau integrator at 1e-6/1e-8 — all 0.7306–0.7310 mJ,
  marginal 7.90%): z = dE/σ_diff = 3741; baseline spread 0.76 µJ; the
  engine resolves twin differences down to 0.6 µJ (P_min ≈ 200 W at this
  marginal). The +0.73 mJ is the model's genuine output, not integration
  noise — no NUMERICAL false positive.
- **Rate uncertainty: the SIGN is not robust.** 13 variants (4 active
  EST rows ×1/3 and ×3 + 5 LIT-conflict alternates; 15 element-bound EST
  rows structurally zero): 11/13 keep the sign (+2.3% … +8.3%), but
  **KSTEP** — `e + Kr* → Kr+ + 2e`, stepwise ionization of metastable Kr,
  EST tier — flips it BOTH ways: ×1/3 → −1.7%, ×3 → −5.6% (positive only
  near the nominal value: a ridge, not a plateau). APS_Ar ×0.1 → +2.3%,
  R8 ×1/3 / ×3 → +3.2%. The heated electrons compete for the Kr*
  metastables that otherwise harpoon into KrF*; which way the net goes
  depends on an ESTIMATED rate.
- **Independent rate tier: the sign FLIPS.** On the PIC_MCC_3D-measured
  KINETIC tier the same twin pair gives **−0.418 mJ, −4.5%** (LOKI
  +0.731 mJ, +7.9%). The cross-tier check does not confirm the co-pump
  gain — and it sides with the DE-AR line "during the discharge it
  hurts". **CORRECTION to 0.1.21/0.1.22: co-pumping is NOT established as
  laser-positive.** The defensible statement is |dE_laser| ≲ 0.15% of the
  pulse at 250 kW (≲ 3% at 3 MW) with an UNDETERMINED SIGN inside the
  model's uncertainty band (−4.5% … +8.3% marginal); the operating-map
  SHAPE (where the magnitude is large: 1–2 atm, plateau window, moderate
  pump) stands, the sign does not. The efficiency-range statement is
  unaffected (|Δeta| < 0.05 pp either way).
- **The afterglow verdict IS robust** (same two stress tests, scheme
  80–300 ns): laser −2.1% (LOKI) / −1.0% (KINETIC) / −1.9 … −2.3% (KSTEP,
  R8 brackets); fluorescence +1.5% / +2.5% / +1.5 … +1.6% — every variant
  keeps laser DOWN and fluorescence UP. The DE-AR reconciliation, the
  pre-prime null (48 µJ couples — kinematics), the keeper null (pulse 2
  unchanged) and the efficiency null are the robust predictions;
  docs/test_range.txt is corrected accordingly.
- Experimental power analysis (α 0.05, power 0.8, interleaved on/off
  pairs): 250 kW plateau +0.15% → 2788 pairs at 2% shot jitter (697 at
  1%, 17422 at 5%); 1 MW +0.6% → 174 (44 / 1089); 3 MW at 1 atm +3% → 7
  (2 / 44). NB the sign test at 250 kW is now a test BETWEEN the model's
  rate tiers as much as of the model.

### Next SPECK step (the measurement the fragility names)

The co-pump sign hangs on how a beam-shaped EEDF under a weak 95 GHz
field partitions energy between Kr* stepwise ionization, harpoon
formation and attachment — exactly the fluid-Te seam PIC_MCC_3D exists
to price. `tools/export_kinetic.py` already has the uniform-field swarm
driver (`run_point`) and the field-free relaxation sampler
(`relaxation_sampling`, the L3 curves); the next build is a third
sampler — a RELAXING population with a weak AC field on (E/N 0.05–1 Td,
95 GHz, ν/ω ≈ 15) — exporting the KrF-relevant per-bucket rates and
the measured KSTEP-class stepwise rate, field on vs off, as a new npz
tier for the window to consume. That measurement decides the co-pump
sign; nothing in the chemistry network can.

### Verified

- **189 passed** from the main checkout (140 s; 178 + 11 MW/timing
  gates incl. `test_mw_window_runs_on_kinetic_tier`). Homebrew 3.14.6
  smoke with `-W error::DeprecationWarning` green on the kinetic-tier
  window (ledger 2.8e-13, transport from the LOKI tables). The ladder
  and the afterglow robustness check both ran FROM THE MAIN CHECKOUT.

## [0.1.22] — 2026-08-21 — the co-pump operating map

Owner ask: "what range of operating parameters allows for co-pumping?"
— the one laser-positive scheme of 0.1.21 (microwave ON during the beam,
+4.1% marginal at the DE-AR point), mapped over the DE-AR ranges.

### Added

- `tools/copump_scan.py`: twin-pair scan of the co-pump marginal
  (dE_laser / E_MW_absorbed) plus the fluorescence marginal, settled E/N
  and peak Te in the window, and a does-the-baseline-lase flag — axes:
  power x pressure at fixed beam hardware (0.5/1/2/3 atm x 50 kW–3 MW),
  pressure at constant deposited power density (j scaled 1/p), beam pump
  rate (j 50/100/200), window placement inside the pulse (full / first
  half / second half / beam + early afterglow / late beam + afterglow),
  optics margin (R_OC 0.90 vs 0.95). 250 ns horizon (the laser is fully
  extracted by ~150 ns; the 400 ns horizon wedged the 3 atm points in the
  deep-vacuum afterglow — the 0.1.21 BDF-Jacobian regime). Not scanned,
  stated: Kr/F2 fractions (the LOKI tables are the v2 fill's; a
  composition axis needs new LoKI-MC tables — LoKI-MC is in-repo).
- `tools/copump_efficiency.py`: post-processes the map into the
  efficiency range — intrinsic (per deposited J) and per injected beam
  energy, baseline vs co-pumped with the microwave in the denominator,
  plus the coupled power needed for +1 percentage point at each point
  (P = d E_dep / ((m − eta0 − d) t_window)).
- `docs/test_range.txt` (Kinetic_Modeling/docs, beside DE-AR_targets.txt):
  the verification sheet — every predicted range of 0.1.20–0.1.22 as a
  testable claim with observable, number + sign, reproduction command
  and pass/fail tolerance (DE-AR 15% milestone class): the co-pump
  operating window and its interleaved-shot test, the efficiency null
  (< 0.05 pp at any feasible power), the afterglow GATED 248.5 nm
  on/off ratio (×1.0 → ×1.17 → ×1.7 → ×5 across 80–300 ns at 250 kW —
  the fluorescence-vs-laser discriminator, measured for the sheet), and
  the pre-prime / keeper nulls.

### Measured — the co-pump operating map (PREDICTIONS; 26 twin pairs,
### krf_hyper_mw, run from the main checkout;
### rebuild_out/copump_scan/copump_scan_report.json)

Marginal = dE_laser / E_MW_absorbed; the DE-AR point (2 atm, 250 kW,
100 A cm^-2, full-beam window) re-measures at +4.06% on the 250 ns
horizon (4.06% at 400 ns — the metric is horizon-insensitive).

- **Pressure is the dominant knob.** Fixed beam hardware (100 A cm^-2):
  0.5 atm baseline barely lases (0.040 J) → +0.35–0.76%; **1 atm → +6.3%
  (the optimum, 0.408 J baseline)**; 2 atm → +4.1%; **3 atm → ≈0
  (+0.01–0.06%)** — the baseline itself is weak there (0.059 J; the
  high-pressure three-body quench chain) and the field is pinned hardest
  (mean E/N 0.03–0.2 Td). At CONSTANT deposited power density (j ∝ 1/p):
  0.5 atm +1.4%, 1 atm +2.5%, 2 atm +4.1%, 3 atm does not lase. Composite
  window: **1–2 atm**.
- **Power barely moves the marginal** (flat to slightly declining: 1 atm
  6.35% @50 kW → 5.76% @3 MW; 2 atm 4.17% → 3.33%); the ABSOLUTE gain
  scales with power — 3 MW coupled at 1 atm buys +12.4 mJ = **+3.0% of
  the pulse**, the first non-noise number in the hybrid studies. Settled
  window fields stay ≤ 1.7 Td mean even at 3 MW; Te_peak 2.2–3.3 eV —
  no over-ionization anywhere in the DE-AR range.
- **Beam pump rate** (2 atm, 250 kW): j 50 → baseline 0.0026 J, −0.05%
  (nothing to amplify); j 100 → +4.1%; j 200 → +1.4% (denser beam
  plasma pins the field lower, mean E/N 0.061 vs 0.086 Td). Moderate
  pumping — lasing comfortably, not saturated — pays best.
- **Window placement is decisive.** First half of the beam (0–37 ns):
  −0.01% (the laser has not built up; early heating buys nothing).
  **Second half (35–74 ns, the lasing plateau): +7.9%** — the same
  +0.73 mJ as the full window for HALF the microwave energy. Beam + early
  afterglow (0–114): −5.5%; late beam + afterglow (40–114): −8.6% — any
  extension past beam-off flips the sign (the 0.1.20/0.1.21 afterglow
  physics). Rule: fire on the plateau, stop AT beam-off.
- **Optics margin amplifies it**: R_OC 0.90 (baseline 0.345 J) →
  **+8.65%** vs +4.1% at 0.95 — the conditioning-gauge effect
  (dlnE/dlng0 ~ g0/(g0−g_th) grows as the margin shrinks).
- **Mechanism signature, consistent on all 26 points**: extracted energy
  up while fluorescence is DOWN (−0.2 to −2.7%): the microwave nudges
  gain up so a larger share of KrF* decays by stimulated rather than
  spontaneous emission. The DE-AR fluorescence diagnostic would therefore
  read co-pumping as a DECREASE while the laser gains.
- **The ceiling, stated plainly**: marginal ≤ ~9% everywhere, and the
  microwave can only carry P_MW x window — at the DE-AR 250 kW on the
  plateau that is 9 mJ → +0.7 mJ = 0.15% of the pulse; MW-class coupling
  reaches a few %. Co-pumping is real and has a clear operating window,
  but this channel cannot approach the 7% → 14% wall-plug goal on its
  own.
- **The efficiency RANGE with co-pumping** (`tools/copump_efficiency.py`,
  post-processing the map; microwave energy in the denominator per the
  DE-AR metric note): intrinsic efficiency per deposited joule over the
  lasing points is **0.12–3.10% baseline → 0.12–3.10% co-pumped** — the
  largest lift on the whole map is **+0.043 percentage points** (1 atm,
  3 MW: 2.588 → 2.631%); the DE-AR point moves 1.536 → 1.537%. Per
  INJECTED beam energy (the 10 cm cell deposits only 3–17% of a 600 kV
  beam) the numbers are 0.01–0.16% either way. Coupled microwave power
  for **one** percentage point of intrinsic efficiency: 67 MW (R_OC 0.90
  cavity), 80 MW (1 atm), 159 MW (plateau window, 2 atm), 270–550 MW
  (full-beam, 2 atm); unreachable where the marginal is below
  baseline + 1 pp. In the limit where microwave dominates the deposition
  the total efficiency asymptotes to the marginal itself (4–9%) — but
  that is GW-class microwave pumping with a beam seed, outside the
  uniform-field idealization, no longer "co-pumping". Per SOURCE joule:
  at the DE-AR coupling (250 kW of a 2 MW gyrotron, ~12.5%) the marginal
  is 0.5–1.1% against the beam's ~7% wall-plug — the beam wins by
  ×6–14; break-even needs gyrotron→gas coupling ≳ 80–90% even at the
  best marginals.

### Verified

- 188 passed from the main checkout with the tool present (134 s; the
  engine is unchanged since 0.1.21 — this release is a driver + map).
  The 250 ns vs 400 ns horizon check on the reference point (+4.06%
  both) pins the metric's horizon-insensitivity; the 3 atm rows lase
  weakly but integrate cleanly at 250 ns.

## [0.1.21] — 2026-08-20 — HYPER-LASER hybrid-timing study (DE-AR frame)

Owner ask: "Look through Kinetic_Modeling/docs/DE-AR_targets.txt — what
about using it to reheat between lasing? or before, to prime the gas a
little to lower the requirements of the e-beam?" The DE-AR doc reframes
0.1.20: the project IS hybrid pumping (HYPER-LASER, DE-AR0002030) — beam
primary + **95 GHz gyrotron booster (<=250 kW coupled)** (not the ARPA
chamber's X-band WR90 hardware), pumped volume the 5–10 cm cube, and the
recorded prior mapping says "afterglow reheat helps efficiency; during
the discharge it hurts" — whose observable is **248.5 nm fluorescence
(prop. to small-signal gain)**, not extracted laser energy, and whose
efficiency metric (the doc's own note) must subtract microwave input.

### Added

- `beam_t0_ns` (beam start delay — pre-prime timing) and `beam2_t0_ns`
  (second identical pulse; overlap refused; envelope sum shares the
  deposition profile and the drive_in row) — engine + params + gates
  (9 in test_mw_drive.py now): envelope shift exact, two-pulse
  drive_in = 2x single to 1e-5, regime line carries "2-pulse beam".
- Preset `krf_hyper_mw`: the DE-AR frame — 600 kV / 100 A cm^-2 /
  60 ns flat into the 10 cm pumped cube, 2 atm Ar/Kr/F2, 95 GHz,
  250 kW, en_max 30 Td (~x5 field buildup over the 250 kW
  traveling-wave value ~6 Td; sits AT the fill's self-sustainment
  crossing). R_OC 0.95 / L_op 0.02 stated.
- `tools/hybrid_timing_study.py`: every timing scheme against its
  no-microwave twin with identical beam timing — afterglow / copump /
  preprime (delayed beam) / keeper (two pulses, window spanning the
  gap, at both the buildup and traveling-wave ceilings) — scoring BOTH
  metrics: extracted laser energy AND the 248.5 nm B->X fluorescence
  integral (n_B/tau_B x h nu x V from the saved history).
- Physics anchor from the packaged swarm tables (sustain_check):
  **ionization = attachment at ~28 Td** in the Ar/Kr/F2 0.949/0.05/
  0.001 fill (ion/att = 0.58 at 25 Td, 1.55 at 30 Td) — the keeper can
  hold a steady plasma only if the coupled field reaches that; below
  it every window is a decay-slower, not a keeper.

### Measured — the timing verdicts (PREDICTIONS, scored against nothing)

All at 250 kW / 95 GHz on krf_hyper_mw (baseline: 31.6 J deposited →
0.485 J at 248 nm, eta 1.54%), each scheme vs its no-MW twin; signs
reproduced at rtol 1e-8 (identical to the printed digits):

- **afterglow** (window 80–300 ns): E_MW 52.5 mJ absorbed; laser
  −1.1 mJ (marginal **−2.1%**) but 248.5 nm B→X fluorescence
  **+0.8 mJ (+1.5% marginal), all of it post-beam** — THE
  RECONCILIATION with the DE-AR line "afterglow reheat helps
  efficiency": their observable (fluorescence ∝ small-signal gain)
  genuinely RISES while extracted laser energy falls. A gain
  diagnostic brightens; an extraction efficiency does not.
- **copump** (window over the beam): E_MW 18 mJ; laser **+0.7 mJ
  (marginal +4.1%)** — the microwave joules convert BETTER than the
  beam's own 1.54% average when delivered into the fully-driven
  plasma; fluorescence −0.2 mJ. NOT a contradiction of the DE-AR
  "during the discharge it hurts (over-ionizes)": that line is about
  DISCHARGE-pumped operation, where over-ionization is an arc
  instability — a channel an externally-sustained e-beam pulse does
  not have. Absolute effect is small (+0.14% of the pulse energy at
  the full 250 kW).
- **preprime** (window 0–195 ns before a delayed beam): E_MW
  **48 microjoules** — the 1e12 m^-3 seed cannot absorb the feed (the
  settled field rides the ceiling, but heating ∝ n_e), and nothing
  reaches the pulse. Priming CANNOT lower the beam requirement by
  energy delivery; the only lever it could ever pull is held
  ionization (the keeper scheme), not stored energy — KrF has no
  storage (Kr* harpoons into KrF* in ~14 ns, KrF* radiates in ~7 ns).
- **keeper** (two pulses, window spanning the compressed 500 ns gap; the undriven twin is vacuum-class n_e ~1e8 m^-3 by then and past ~550 ns the BDF Jacobian overflows even with restarts):
  at the 30 Td ceiling the window absorbs its full 100 mJ and holds n_e
  at **3.5e16 m^-3 at pulse-2 arrival vs 9.7e7 undriven (x3.6e8)**; at
  the 6 Td traveling-wave ceiling 53 mJ absorbed, n_e 6.3e15 (x6.5e7 —
  heated electrons attach far slower, so even a sub-sustainment field
  holds the decay on this timescale). **Pulse 2 gains NOTHING**: 0.2920 /
  0.2924 J held vs 0.2928 J unheld (−0.3% / −0.1%) — an e-beam pulse
  makes its own ~1e20 m^-3 plasma in nanoseconds; a 1e16 seed is x3e-4
  of it. Total laser marginal −1.9% / −2.8%, fluorescence +2.4% / +1.6%
  (the held plasma glows). The cost scales with the gap: the keeper
  draws the feed for the whole hold — at any real rep rate (ms gaps)
  that is x1e2–1e4 the laser pulse energy. Pulse 2 lases 0.293 J vs
  pulse 1 0.485 J in BOTH arms (F2 burn 21%, no gas flow in the model —
  a real device refreshes gas; the twin difference isolates the
  microwave).

### Fixed

- Deep-afterglow hardening, found by the keeper runs: te_of/rhs divided
  u_e by the quasineutral n_e with no floor — the undriven afterglow
  past ~500 ns reaches n_e → 0 (nine decades down by 400 ns in the
  hyper cell) and the 0/0 NaN'd the numerical Jacobian. A 1 m^-3 floor
  on the te divisor (physically nothing, numerically decisive) +
  a gate (rhs finite with every charged species at zero). Separately,
  a SINGLE BDF solve across the full vacuum decay still overflows its
  numerical-Jacobian factor state (scipy adaptive-perturbation
  internals); measured cure: restart the solver — the keeper schemes
  integrate in checkpoint-restart segments (the checkpoint module's
  stated tolerance-exact semantics), and the probe shows chunked
  restarts walk the same trajectory cleanly.

### Verified

- **188 passed** under `/usr/bin/python3` (135 s; 178 + the 10 MW/timing
  gates incl. the zero-n_e RHS guard), and **188 again FROM THE MAIN
  CHECKOUT after delivery** (140 s). Homebrew 3.14.6 smoke with
  `-W error::DeprecationWarning` green through the new paths (hyper
  preset, 2-pulse beam, MW window; ledger 1.5e-13).
- Every scheme's signs and magnitudes reproduced at rtol 1e-8 (afterglow
  −2.12% / +1.54%, copump +4.07% / −0.90% — identical to the printed
  digits), and the full five-scheme study re-run FROM THE MAIN CHECKOUT
  reproduces the sandbox numbers digit-for-digit
  (rebuild_out/hybrid_timing/hybrid_timing_report.json).
- Ledger closure on every arm: 1e-12 class through windows and pulse 2;
  the segmented keeper's cumulative channels ride the checkpointed state
  (final totals = the sum the single-solve would have produced).
- `verify.py` **ALL CHECKS PASS** (imports-alone, per-file, combined
  forward + reversed at 188 each, composition ladder).

## [0.1.20] — 2026-08-14 — the ARPA-chamber microwave-reheat design study

Owner ask: "set up a simulation of an e-beam lasing a gas in this chamber
geometry [the ARPA-E vacuum chamber measured into
freecad/arpa_chamber_ports.json], and in the afterglow as the gas cools,
reheat the gas with the microwave to attempt to improve efficiency."

### Added — the microwave reheat window (a beam-drive overlay)

- `mw_power_W` / `mw_en_max_td` / `mw_t0_ns` / `mw_flat_ns` / `mw_rise_ns`
  / `mw_freq_GHz` (params.py): a trapezoidal INJECTED-power window on the
  beam drive. Each RHS step the uniform cavity field settles at the
  lowest E/N where the swarm power balance absorbs the feed,
  P_abs = e E v_d(E/N) x nu^2/(nu^2+omega^2) x N_e — capped at
  `mw_en_max_td` (the empty-cavity ceiling the wall Q supports; excess
  power counts as reflected, never deposited). Absorbed power enters
  U_e and drive_in through the same integrated ledger as every drive;
  collision rates stay on the drive's <eps> = (3/2)Te keying (heating
  pushes <eps> toward the swarm point of the settled field, where the
  in-table rates ARE field-driven-EEDF physics). Engine blocks:
  `Device.__init__` (grid precompute), `Device._mw_field_td` (lowest-
  crossing solver — argmax, not searchsorted: the Ar NDC dip makes
  en x v_d non-monotone), the rhs overlay, validity gauges
  (`mw_*` keys), and the regime line.
- WHY power-capped and not a prescribed field: the first cut prescribed
  E/N directly, and 8 Td held across the dense afterglow is an infinite
  generator — it drove a full avalanche to ionization fraction 1.0 and
  deposited 7.2 MJ (ledger closed at 1e-16: the bookkeeping was right,
  the model assumption was wrong). A feed cannot do that; the power
  balance is the physical statement of "three WR90 feeds".
- `r_oc` / `l_op` (params.py): resonator design choices with sentinel
  defaults (0 / negative = mode defaults, gate-tested). The ARPA column
  is 16.9 cm of gain — the Kannari-class L_op 0.088 would put threshold
  at 7.4 /m; the study states R_OC 0.9 / L_op 0.02 research-cell optics
  (g_th 0.37 /m) as construction inputs.
- Presets `krf_arpa_mw` (0-D) and `krf_arpa_mw_3d` (4x6x4, beam axis y
  resolved): pumped column = the O190.5 front-aperture beam
  (area-equivalent 168.8 mm square) through the 264.16 mm chamber depth;
  optical axis z = the top/bottom window axis (internal mirrors stated —
  the physical O76.2/O31.75 windows would vignette an external cavity);
  2 atm (1 atm gauge, vessel-safe) Ar/Kr/F2 0.949/0.05/0.001; 500 kV,
  140 A/cm2, 100 ns flat (~0.9 MW/cm3 peak deposition, Electra-class;
  ~35% of injected beam energy deposits — the 26.4 cm-atm column passes
  the rest to the backplate, stated).
- `tools/mw_reheat_experiment.py`: the study driver — baseline
  (mw_power_W = 0) vs a scanned injected-power ladder, marginal
  efficiency dE_laser/E_MW per rung, absorbed fraction, settled-field
  window profile, Te/n_e response, and the WR90 feed-budget context
  (~1 MW), written to rebuild_out/mw_reheat/mw_reheat_report.json.
- Gates (tests/test_mw_drive.py, 8): envelope shape; beam-drive-only
  refusal; declaration/bounds; cavity sentinel defaults; field-solver
  limits (ceiling for thin plasma, sag for dense, in-balance absorbed
  power reproduces the feed to 2%); never-reached window bit-inert;
  heats-inside-window with S3 closure through the window; ledgered MW
  energy matches the swarm-conduction quadrature to 5%.

### Measured — the design verdict (a PREDICTION, scored against nothing)

Study: `tools/mw_reheat_experiment.py` on `krf_arpa_mw` (0-D), 241 saves.

- BASELINE: 554.1 J deposited (26% of the 2.13 kJ injected beam; the
  26.4 cm-atm column passes the rest to the backplate, stated) →
  **11.95 J at 248 nm, eta_intrinsic 2.16%**, g0 peak 1.95 /m vs g_th
  0.37 /m; n_e 2.7e20 m^-3 at window start (130 ns) collapsing to
  3.8e13 by 300 ns undriven; F2 burn 43.7%.
- **REHEAT NEVER IMPROVES THE LASER.** Marginal efficiency
  dE_laser/E_MW: −6.3% at 0.3 MW, **−6.8% at 1 MW (the stated 3×WR90
  budget)**, −6.9% at 3 MW, −5.9% at 10 MW; starting the window AT
  beam-off (112 ns) is worse, −12.8% at 1 MW. The what-would-it-take
  bracket (en_max 300 Td): STILL negative at 100 MW (−0.96%) and 1 GW
  (−0.13%, 160 J absorbed = 29% of the beam's own deposit).
- WHY (visible in the histories): the window absorbs 100% of the feed
  at a settled field that never exceeds ~2 Td — **the plasma the window
  sustains is too conductive for the field to climb into the
  direct-excitation regime at ANY power** (more power → more n_e → the
  field pins lower: the 1 GW rung ends at 1.5 Td). At the resulting
  Te 2.5–3.3 eV the electrons slow the cold-electron DR / ion–ion
  harvest that feeds KrF* in the tail and quench KrF* — the laser tail
  shrinks. No re-avalanche to lasing conditions exists on this path.
- WHAT THE WINDOW DOES DO: hold the afterglow ionization up 4–7 decades
  (n_e at 300 ns: 3.7e17 at 0.3 MW, 1.1e18 at 1 MW, 3.8e20 at 1 GW, vs
  3.8e13 undriven) at Te ~ 3 eV — an inter-pulse preionization keeper
  for rep-rate operation, not an efficiency booster.
- Cross-checks on every rung: E_MW by ledger subtraction vs the
  swarm-conduction quadrature agree to print precision; ledger residual
  max 9.1e-13 THROUGH the window; nu/omega 148 at the cap (deep in the
  collisional regime, the AC factor is honest); absorbed fraction 100%
  on all rungs (the feed budget, not coupling, is the binding limit).
- Stated limits carried on the preset/tool docstrings: uniform MW field
  (overmoded cavity; cm-class collisional skin at pulse-peak n_e means
  early-window heating is really shell-localized — best late in the
  afterglow); fluid-Te rate keying; internal mirrors R_OC 0.9 /
  L_op 0.02 stated as construction (the physical Ø76.2/Ø31.75 windows
  would vignette an external cavity); no wall cooling of the gas
  (ms-scale). "As the gas cools": the ns-scale cooling the window
  reverses is the ELECTRON gas; T_gas (355 K at end of the baseline
  run) barely moves on these timescales.

### Verified

- **186 passed** under `/usr/bin/python3` (223 s; 178 + the 8 MW gates),
  and 186 again FROM THE MAIN CHECKOUT after delivery (230 s, LXCat via
  the repo markers). `verify.py` **ALL CHECKS PASS** (imports-alone,
  per-file, combined forward + reversed, composition ladder). Homebrew
  3.14.6 smoke with `-W error::DeprecationWarning` green through the MW
  path (0-D beam+window run: success, ledger 1.6e-12, regime line
  carries the window).
- Ledger closure through the window on every scan rung: residual max
  9.1e-13; E_MW by ledger subtraction vs the swarm-conduction
  quadrature agrees to print precision on all 7 rungs (the
  test_mw_drive gate pins 5%).
- Mechanism cross-checked in the histories (not asserted): the photon-
  side deltas at 1 MW are laser −0.0109 J, HR/window −0.0021 J,
  radiated −0.0033 J — the window slightly degrades net gain during
  the 140–200 ns cavity-field drain (F⁻ down ×0.87 at 160 ns, Kr2+ up
  ×1.5 at 200 ns, g0 lower by ~0.5–1% while positive) — while the
  remaining ~0.14 J of the absorbed 0.16 J thermalizes
  (electron → gas heat). Harm is small in absolute terms; benefit is
  zero everywhere scanned.

## [0.1.19] — 2026-08-12 — afterglow visibility: √ tone on the continuum

Owner ask: an mp4 of "gas being excited, forming a plasma, then
cooling in an afterglow." The excite→plasma part was already the glow
movie; the COOLDOWN went black — n_e² falls five decades after
beam-off and the linear brightness map rendered the whole tail dark.

### Changed
- The recombination continuum's frame brightness now uses a √ tone
  (`GLOW_COMPONENTS[..]["tone"] = 0.5`, matching `Math.pow` in the
  viewer's glow mode); the fluorescence stays linear (its collapse IS
  the lasing story). Stated on every caption. The npz keeps exact
  fields; the spread readout is untouched.
- New long-afterglow artifact: `rebuild_out/viz_krf_afterglow/`
  (krf_beam_3d preset run to 400 ns — pump ends ~65 ns, so ~5/6 of the
  movie is the cooldown), rendered as the glow movie.mp4/gif + viewer.

### Fixed
- `glow_cloud` allocation: per-fine-voxel `rint` floored any budget far
  below the fine-voxel count to ZERO everywhere — the deep afterglow
  (6–20 dots over 30k voxels) rendered black even though the physics
  said embers remained. Now largest-remainder cumulative rounding:
  deterministic, total exact (gated, incl. a 7-dot regression case).
  Physics note recorded with it: the visible afterglow of THIS cell is
  genuinely short — F2 attachment + recombination kill n_e within
  ~200 ns (fluor ≥5 dots to ~268 ns; the √-toned continuum holds
  15–30 ember dots to 400 ns). There is no neon-tube-style long glow
  at 1 atm with F2 present; the strip carries the quantitative decay.

## [0.1.18] — 2026-08-12 — the glow view: render what the cell EMITS

Owner feedback on 0.1.17: "I was expecting to see plasma in a container;
I still only see a box that changes color." Root causes, both fixed:
the default-named `movie.gif` was the VOXEL movie (the analysis view —
literally a box changing color), and every view so far spoke data-viz
(voxels, dot lattices, colorbars), not camera. The plasma emits light;
now there is a view that renders exactly that.

### Added
- GLOW style (`--style glow`; part of the new default `--style all`),
  and it is now the HEADLINE `movie.gif` — the voxel/field movie moved
  to `movie_fields.gif` (frames dirs `frames_glow/`, `frames_fields/`,
  `frames_particles/`). The camera view: a dark vessel (edges + optical
  end windows drawn, the 75 cm axis laid horizontally like a bench
  cell) filled with a luminous gas. Brightness ∝ emission, two
  physically-motivated components each tone-mapped to its own run peak
  with the display gains STATED on every frame: KrF*/ArF* B→X
  fluorescence (violet-white; the real line is UV) and the n_e²
  recombination continuum (blue, gain 0.35). Fields are trilinearly
  interpolated between cell centers FOR DISPLAY (stated; the model
  still resolves its coarse grid); luminous point clouds are sampled
  linearly from the refined emission (same deterministic mulberry32
  lattice), drawn halo+core, finished with a PIL bloom pass. Default
  display aspect is √-compressed (stated on the frame; `--aspect
  true|unit` unchanged) so the duct stays a duct while the
  cross-section remains visible.
- viewer.html gains a GLOW mode, now the default (radio glow /
  particles / voxels): additive radial-gradient sprites, brightness
  from the SAME trilinear interpolation implemented in JS over the
  embedded fields, per-component peaks from the stored field limits.
- Gates (tests, 17 → 19): emission/refine/cloud math (deterministic,
  budget-honoring, empty-on-dark, shape-correct trilinear refinement);
  glow frame + bloom render headless; the style→filename contract
  (glow = `movie.gif`, voxels = `movie_fields.gif`); the viewer's glow
  radio is the checked default.

### Measured
- Exposure iterated against real frames (three rounds, inspected):
  45k samples at α 0.05/0.16 blew the true-aspect duct out to a white
  slab; shipped at 30k samples, α 0.028/0.085, bloom r4 g0.8,
  √-compressed aspect — wispy translucent plasma with the vessel and
  the Bragg asymmetry visible, afterglow dimming and granulating.
- Both chamber demos re-rendered with all three styles from their
  saved npz (no re-integration).

## [0.1.17] — 2026-08-12 — uniformity stated as a number; the Bragg preset carries the structure

Owner feedback on 0.1.16: "still appears to be perfectly uniform with
absolutely no deviation or local filamentation/grouping." Measured, the
observation is CORRECT and is the physics of that preset, not a viewer
defect: at the Electra-class point (500 kV through 10 cm of 3 atm) the
per-cell n_e spread is ×1.014 (deposition shape ×1.018) — an e-beam
laser is DESIGNED to pump uniformly. And no filamentation can ever
emerge from this engine: cells are independent reactors and the only
inter-cell couplings (cavity photon field, volume-averaged ASE) are
uniformizing — spatial structure here can only be what the deposition
IMPOSES. Self-generated grouping is kinetic physics and lives in
PIC_MCC_3D (0.8.3's two-stream demo, gated against the closed-form
growth rate).

### Added
- The spread is now STATED on every particle frame ("cell n_e spread
  ×…"), and in the viewer.html particle readout — a uniform-looking
  cloud now says the plasma IS uniform instead of leaving it to the
  eye. The krf_beam_3d preset doc carries the ×1.02 number.
- `krf_beam_3d_bragg` preset (360 cells): 100 kV into 6 cm at 3 atm —
  the recorded mesh-refinement Bragg case. The deposition rises ×2.9
  to the Bragg peak at x ≈ 4 cm and is ZERO past the stop: the
  particle cloud fills the first two-thirds of the chamber and dies,
  with the spread readout in the 1e5+ class. Gated:
  `test_bragg_preset_has_real_spatial_structure` (peak inside the
  domain, zero tail, ×2+ rise — from the built device, no run needed).

### Measured
- Electra-class demo (krf_beam_3d): n_e cell spread ×1.014 at 67 ns,
  ×1.003 by 200 ns; deposition ×1.018 — the uniformity is real.
- Bragg demo (krf_beam_3d_bragg, `--style particles`): run wall
  1453 s (360 cells, stiffer than the uniform case — zero-drive cells
  beside ×2.9-drive cells; ledger residual 1.25e-10, drive 36.4 J).
  Measured n_e cell spread ×2.4e9 at 67 ns, ×3.0e9 at 133 ns — the
  particle cloud fills to the Bragg stop and dies, nine decades of
  real imposed structure vs the Electra case's ×1.014. Artifacts in
  `rebuild_out/viz_krf_beam_3d_bragg/`.
- Suite 176 passed under `/usr/bin/python3` (131 s; 175 + the Bragg
  gate); Homebrew 3.14.6 smoke green (now also covers the PIC
  two-stream path).

## [0.1.16] — 2026-08-12 — particle view: the fields drawn as what they count

Owner feedback on 0.1.15: "I do not see any particles or macro-particles
of the plasma, I instead see the cube just change hues." Correct — the
device grid is a rate-equation model; there ARE no particles in it, and
voxel hue was the entire state. Two answers landed, one per package:

### Added
- PARTICLE VIEW in `tools/visualize_3d.py` (now the default animation
  style; `--style {both,particles,voxels}`, `--pts`): n_e and the upper
  level drawn as dot clouds. Dot COUNT per cell is LINEAR in the real
  particle number (density × cell volume), so every dot carries one
  stated macro-weight — the 240-cell demo renders 1 dot ≈ 7.5e13
  electrons / 6.1e12 KrF* at `--pts 8000`, ignition → plateau →
  afterglow arriving as dots appearing and vanishing. Dot POSITIONS are
  a deterministic per-cell hash lattice (mulberry32, bit-identical in
  Python and the viewer's JS) — a density SAMPLING, labelled "not
  trajectories" on every figure. New outputs `movie_particles.gif/.mp4`
  + `frames_particles/`; `viewer.html` gains a particles/voxels toggle
  (particles default) with per-species macro-weights displayed; the npz
  now carries `cell_vols` (older npz load with a uniform-grid fallback).
- REAL macro-particles, where they actually exist: PIC_MCC_3D 0.8.2's
  `tools/visualize_particles.py` (see that changelog) — every dot IS a
  PIC macro-particle from the leapfrog, in that code's honest
  reduced-density sub-chamber regime.
- Gates (`tests/test_visualize_3d.py`, 13 → 16): sampler determinism +
  [0,1) range; dot counts linear in particle number with the
  weight×dots = real-number identity; particle frame +
  `movie_particles.gif` render headless; the viewer payload carries
  cellvol + per-species pnorm and the mode control; npz round-trip now
  includes cell_vols.

### Measured
- Both JS viewers re-verified headlessly (deno stubbed-DOM → rasterized
  draw calls): REBUILD particle mode = 6000 dots (3000/species at the
  plateau) + 12 chamber edges at yaw −35°; PIC viewer = 1728 dots with
  the v_z slosh stratification visible.
- Suite 175 passed under `/usr/bin/python3` (120 s); Homebrew 3.14.6
  smoke with `-W error::DeprecationWarning` green.

## [0.1.15] — 2026-08-12 — see the run: dense per-cell histories + the 3-D visualizer

Owner ask: visualize the plasma in the chamber in 3-D over the pulse. The
engine already computed per-cell state, but `_collect` volume-averaged the
histories and kept only `y_final`; the `dense_cells` argument in
`Device.run`'s signature was DEAD (accepted, never read), and nothing in
the package could draw anything.

### Added
- `Device.run(dense_cells=...)` HONORED: keeps the full save-time state as
  `res["y_t"]` (`lay.size × n_save`; measured 8.6 MB at the 240-cell
  preset, 121 saves). Default `None` → auto: ON for multi-cell runs, OFF
  for 0-D. `checkpoint.resume` applies the same rule, so a resumed 3-D
  run keeps its cells too.
- `tools/visualize_3d.py` — A VIEWER, NOT AN ORACLE (feeds nothing back):
  per-cell n_e, upper-level density, local g0 (`gain_m1`'s expression
  vectorized over cells), Te, Tg + the static deposition shape, all
  reconstructed from the run state via the engine's own accessors
  (`StateLayout.density_block`/`electron_density`, `Device.te_of`).
  Outputs: `viz_data.npz` (the data of record — exact arrays + meta),
  `contact_sheet.png`, frame PNGs → `movie.gif` (Pillow) / `movie.mp4`
  (`--mp4`, ffmpeg), and `viewer.html` — a self-contained canvas voxel
  viewer (time/yaw sliders, field selector, unit/true-aspect toggle;
  fields QUANTIZED to uint8 with scale + limits stated on the page —
  exact values live in the npz beside it). Every figure carries the
  scope caption (independent reactors, shared cavity/ASE overlay, not a
  PIC/fluid solve). `--npz` re-renders without re-running; `--no-movie`
  is the fast npz → viewer/sheet pass.
- `tools/__init__.py`, and the README's tools/ contract mechanised:
  `test_modules_never_import_tools`.
- `tests/test_visualize_3d.py` (13 gates): dense auto-on/off/explicit;
  volume averages RECOMPUTED from `y_t` reproduce `_collect`'s histories
  (rtol 1e-12); per-cell g0 volume-averages to `_collect`'s gain
  (linearity, rtol 1e-10); npz round-trip exact; quantize/dequantize
  inverts within one 8-bit bin (lin + log); frame, contact sheet and GIF
  render headless (Agg); viewer.html has no external references and its
  embedded data dequantizes to the true field within one bin.

### Fixed
- Viewer camera geometry (pre-ship): the first canvas projection had the
  screen-up/farness components sign-swapped and painted near-to-far.
  Re-derived as an explicit orbit-camera basis (yaw ∘ 28° elevation),
  painter far-first, faces culled by the normal's farness sign — then
  verified HEADLESSLY by executing the shipped HTML's own script under
  deno with a stubbed DOM and rasterizing the recorded draw calls:
  720 fills = 3 visible faces × 240 cells at yaw −35°, and exactly
  480 = 2 × 240 at yaw 0 (x-faces edge-on culled).

### Measured (the 240-cell demo on this checkout)
- `python3 REBUILD/tools/visualize_3d.py --preset krf_beam_3d --mp4`:
  run 230.6 s (regime unchanged: drive 26.76 J, ledger residual max
  3.2e-11), `viz_data.npz` 150 KB, 61 frames 19 MB, `movie.gif` 5.5 MB,
  `movie.mp4` 285 KB, `viewer.html` 284 KB. The fields tell the pulse's
  story: seed → beam-on plateau (n_e ≈ 2e20 m⁻³, g0 ≈ +0.5 /m, mild
  10-cm deposition gradient) → post-beam recombination decay with ⟨g0⟩
  swinging absorber-negative (≈ −0.2 /m).
- Verified: pytest 172 passed under `/usr/bin/python3` (121 s; 159 + 13
  new) and the Homebrew 3.14.6 / numpy 2.5 smoke with
  `-W error::DeprecationWarning` covering a tiny 3-D run + every output
  kind (npz/PNG/GIF/HTML).

## [0.1.14] — 2026-08-05 — the three amplification chains: repulsive-F2 reroute, ion-ion screening, Ar** priced, conditioning gauges

Owner ask: "Proceed in order all three chains" — the amplification-chain
audit of the x0.21/x0.19/x0.06/x0.35 residuals.

### CHAIN 3a — F2's repulsive states routed as physics, not bookkeeping
F2's a3Piu / A1Piu are REPULSIVE: excitation IS dissociation within a
vibrational period. Chen's Table III books them as bound species that
wait for a SECOND electron collision (R13/R14) — at beam-mode ne that
parks excited F2 for ~us. In the Kannari-chemistry modes the measured
LOKI excitation rows now deliver e+F2 -> F+F+e directly (same measured
rate; the electron pays the full 3.16/4.34 eV threshold — D0 = 1.6 eV
stored in the fragments, the excess as fragment heat via
e_cost_override) and R13/R14 are dropped; **arf_chen keeps Chen's
two-step bookkeeping VERBATIM** (anchor fidelity). R11/R12 become the
sole ALLOWED_SAME_CHANNEL exemption (two repulsive states, two measured
rates, same fragments — documented in check_network). Side effect,
honest: Te rises 2.03 -> 2.3 (x1.15 vs NRL) — the two-step bookkeeping
had been double-charging the electron gas (3.16 excitation + 4.69
dissociation per molecule).

### CHAIN 3b — ion-ion Debye screening (chemistry.ion_ion_screening)
Kannari 1985 p.4313 (Morgan [35]): plasma shielding cuts ion-ion
recombination at pulse densities — he prints 0.397 at 1e14 cm^-3 pairs
in 2-atm Ar and applied it as a CONSTANT. The engine now applies the
parameter-free screened-capture form S = 1/(1 + R_c/lambda_D) (R_c the
(3/2)kTg Coulomb radius; lambda_D over all charge carriers at Tg) to
the flannery/lee_ii kinds, computed live per evaluation. At his printed
point the form gives 0.695 vs Morgan's MC 0.397 — a MILDER cut, x1.75
apart, RECORDED as part of the error budget (Morgan's paper is not
in-repo). S -> 1 as the plasma vanishes; gate-tested.

### CHAIN 2 — the Ar** head, priced (not guessed)
The Ar** x0.06 / alpha residuals head at APS_Ar (Ar**+Ar -> Ar*+Ar
1.0e-10): Kannari marks it 'estimated', Lee 1989 Table I prints the
same 1.0e-10 WITH a citation [36] (read first-hand) — two in-repo
prints agree — while NRL's Orestes carries NO neutral 4p->4s quench at
all (petrova2022 A3 relaxes Ar** by electrons only): a STRUCTURAL
inter-model disagreement worth x15 on Ar**. The two-print value is
KEPT; the x0.1 alternate is now SWEPT every est-sweep
(CONFLICT_ALTERNATES; measured at the v2 point: g0 -6.7%). The two
printed electron de-excitation rows are adopted (AR_EDX1/AR_EDX2,
petrova2022 A3 [Vriens-1980/McGeoch-2018], superelastic — energy to
the electron gas).

### CHAIN 1 — conditioning gauges on every anchor point
scorecard._conditioning: each device-anchor row now prints g_th, the
threshold margin, and dlnE/dlng0 ~ g0/(g0-g_th) — an amplified miss
can no longer be mistaken for a raw physics miss. MEASURED: the
Kannari R_OC 20% point reports margin 11%, amplification ~x9 (its
x0.20 is a ~x1.2 gain error through a x9 amplifier); every other
anchor point runs at margin 35-47% (amplification x2-3).

### THE MEASURED STATE (all live, nothing tuned)
- Population benchmark: Te x1.15, Tg x0.95, Ar* x0.87, alpha x0.23 ->
  x0.34 (the screening keeps more Ar2+ absorbing), g0/alpha 18-24 vs
  their gamma ~17.5 (was 27-35); ArF* x0.28, Ar** x0.05, F- x0.19,
  F2-burn x0.43 (now GENUINE dissociation — no parked states).
- Electra: x0.84 / x1.41 / x1.76 (high-pressure overshoot shrinking;
  trend still unreproduced — the 3-body class NRL's own model shares).
- Kannari: x0.83 at R_OC 80% (margin 47%); x0.20 at 20% (the x9
  amplifier, labeled).
- v2: g0 0.0052 /cm, Te 2.0 eV, closure 7e-13; est-sweep re-measured
  (Ar fill worst lever = the APS_Ar x0.1 alternate at -6.7%, the 15
  element-bound unchanged; Ne fill BC/CB_Ne +-68-70% unchanged).

### Verified
- **159 passed** warning-clean (the chain gate
  test_repulsive_f2_reroute_and_screening: reroute + overrides +
  arf_chen verbatim + ALLOWED_SAME_CHANNEL + screening printed-point
  and limits; census CHEN >= 13 re-pinned with history). verify.py
  ALL CHECKS PASS; Homebrew smoke clean; JXA board green.
- Board 58 items: 57 ok / 1 warn / 0 gap (the audit item is ok; the
  warn stays with every residual either shrunk, priced, or labeled).

## [0.1.13] — 2026-08-05 — the Ramsauer stall broken: the export deepened, the F- curve measured

Owner ask: "there is still one open modules, proceed with fixing it" —
the board's one warn, whose dominant quantified residual (F- x0.02 vs
the NRL benchmark) named the export deepening as its strongest fix.

### THE RECORDED FOLLOW-UP, DONE (it was 105 seconds of compute)
`tools/export_kinetic.py --deepen N`: re-runs ONLY the relaxation
sampler — it is independently seeded (seed+5000), so the swarm points
and W-cascade arrays are carried over UNTOUCHED — with max_rounds
raised 25000 -> 250000. Pure extension, mechanically enforced:
- A 25000-round reproduction run first matched the shipped npz
  ARRAY-FOR-ARRAY (NaN-aware) — the replay property proven end-to-end.
- The deepen run's gate froze 99 above-stall bin values bit-identical
  (same seed => rounds 1..25000 replay exactly; the population mean is
  monotone, so no old bin above the stall can gain tallies).
- The stall was max_rounds, NOT attrition: the Ramsauer region is
  cheap per round (small majorant), and the full deepening cost 105 s.
  The 0.1.9 estimate that this needed a compute campaign was wrong.

### WHAT THE NEW BINS MEASURE
- **F-: the dissociative-attachment curve RISES x13 toward low eps**
  (3.5e-16 @ 1.76 eV -> 1.25e-15 @ 1.26 -> 2.90e-15 @ 0.91 ->
  4.60e-15 @ 0.66, stats 4-17%) — exactly where the old taper FELL:
  the F- starvation mechanism, now measured data the bridge consumes.
  (The 0.47-eV bin is a single event, err 1.0 — the >=50%-err filter
  rightly drops it; usable floor 0.66 eV.)
- **el:Ar through the Ramsauer minimum** (1.09e-14 @ 1.26 ->
  2.54e-15 @ 0.47, sub-1% stats down to 0.66): the taper sat x0.42
  LOW there — it over-dragged above the old floor and under-drags
  below it; both directions now pinned in the taper-band gate.
- vib to 0.66 eV (the old 100%-err 1.26 bin now carries 3% stats).

### THE MEASURED CONSEQUENCE
- NRL population benchmark: **F- x0.02 -> x0.19** (1.92e13 vs 1e14
  cm-3); Te x1.01 and the rest hold.
- Kannari x0.81 -> **x0.85** (R_OC 80%); Electra x0.85/1.45/1.87 ->
  **x0.88/1.48/1.89** (trend still unreproduced — the 3-body class,
  shared with NRL's own Orestes).
- Re-pins with history: the vib coverage floor (1.76 -> 0.66), the
  bin-selection in the bridge purity gate now mirrors the bridge's
  own err filter, the taper-band gains the deepened low side
  (0.35 < r and the F- taper-LOW direction).
- Board warn residuals after this: F- x0.19 (consumption-side /
  beam-direct-attachment candidates remain), Ar** x0.06, the
  high-pressure trend, threshold-region points.

### Verified
- **158 passed** (/usr/bin/python3; kinetic+beam gates 22 in
  isolation first) with the three honest re-pins (vib floor, bridge
  purity-gate bin selection, two-direction taper band); v2 unchanged
  (g0 0.0054, Te 1.9 eV) with the coverage gauges now reading <eps>
  in-table 52% / in-table-or-relax 68% (was 0% pre-0.1.12 — the
  gauges reported the gap, now they report the recovery); closure
  1e-13 class. verify.py ALL CHECKS PASS; Homebrew smoke clean; JXA
  board green. The deepened npz ships in
  MODULES/kinetic_tables/pic_mcc_3d_swarm.npz with provenance
  (meta.relaxation_deepened: max_rounds 250000, date, note).
- Board 57 items: 56 ok / 1 warn / 0 gap (the warn's residual list
  shrinks again; the L3 item records the deepening; the OMISSIONS
  Ramsauer follow-up is RESOLVED 0.1.13).

## [0.1.12] — 2026-08-05 — the fluid-Te suspect tested: the prompt-floor fix; Te x0.26 -> x1.01

Owner ask: "Proceed with treatment suspect" — the first-listed
off-network suspect from 0.1.11's cross-device diagnosis.

### THE BENCHMARK — the suspect put on the record
New `scorecard.electra_population_rows`: the NRL Orestes POPULATION/
TEMPERATURE benchmark (petrova2022 Figs. 4/6/8 + text, first-hand),
scored live at 0.5 MW/cc, 12.5 psi, 0.5% F2. MEASURED (before the
fix): Te x0.26 (0.53 vs ~2 eV flat), ne x2.9 HIGH, **F- x0.003**,
ArF* x0.21, g0 x0.40, Tg x0.95, ion ordering reproduced. The deficit
mechanism in one table: cold bulk -> the attachment key sits far
below the measured 1.76-eV floor where the taper FALLS while physical
dissociative attachment RISES -> F- collapses -> the ion-ion channels
(44% of NRL's ArF* formation, their Table 3) are starved -> the
cross-device x0.5 yield deficit.

### THE ROOT CAUSE — derived, not fitted
The Spencer-Fano cascade is TIME-FREE: below the noble-gas inelastic
floor only trace-F2 channels stay open, so the recursion funneled
100% of the sub-excitation flux through them at t=0 — physics that
takes ~ns per event at 0.1-0.5% F2 (tens of events spanning the whole
pulse, NOT prompt) — while the network's own e+F2 rows carry the same
pathway BY RATE: one physical pathway, two owners (the 0.1.11 lesson
at the energy level; the cited w_table convention always banked this
energy with the electron gas). Measured: eps_sec 1.81 eV/pair at
0.5% F2 vs ~4.3 for pure Ar — the trace was draining the Platzman
share.

### THE FIX — `prompt_floor_eV` in the cascade
`spencer_fano_per_process` takes a prompt floor (default 0 = old
behaviour); `spencer_fano_yields` sets it to the minimum noble-gas
inelastic threshold present. Sub-floor trace channels are excluded
from the prompt recursion; their energy stays with the electron gas
(eps_sec 1.81 -> 4.41 eV/pair at 0.5% F2; 3.67 at the v2 krf fill;
subfloor -> 0, F2 vib is bulk chemistry). Pure-noble cascades are
BIT-IDENTICAL — every W anchor and the SF-vs-MC cross-check unmoved
(gate-tested: test_prompt_floor_banks_subexcitation).

### THE MEASURED CONSEQUENCE — nothing tuned, everything re-scored
- Population benchmark: **Te x1.01** (2.03 vs 2 eV), Ar* x0.84,
  Tg x0.95, ne x1.9, g0 x0.54; residuals now sharply localized:
  **F- x0.02** (up from x0.003, still collapsed — the
  export-deepening follow-up's strongest reason), Ar** x0.06, F2
  dissociation x0.41.
- Kannari device: x0.56 -> **x0.81** at R_OC 80% (x0.21 at the
  threshold-sensitive 20% point); g0/alpha 3.0-3.1.
- Electra ArF: x0.07/x0.53/x0.69 -> **x0.85/x1.45/x1.87** — the
  robust low-pressure point nearly closes; the high-pressure points
  now OVERSHOOT with the same unreproduced pressure trend NRL's own
  Orestes has (their stated 3-body suspicion). THE OLD SINGLE x0.5
  CROSS-DEVICE DEFICIT STORY IS BROKEN; the board warn is re-scoped
  to the quantified residuals (F-, Ar**, trend, threshold points).
- v2: g0 0.0034 -> **0.0054 /cm**, Te 0.52 -> 1.88 eV — beam-mode
  <eps> now sits INSIDE the measured relax window (1.26-3.4 eV): the
  L3 bridge consumes measured bins where it used to extrapolate.
  Closure 1e-12 class throughout.
- est-sweep re-measured on the corrected split: Ar fill worst EST
  lever -5.9% (URAD_Bv x3; KR2_RAD conflict falls to +3.4%; the
  stepwise levers R8/KSTEP now register 1-3% on the hotter bulk —
  physically sensible); Ne fill BC/CB_Ne +-66-69% unchanged. The S5
  pricing HOLDS and tightens.

### Verified
- **158 passed** (/usr/bin/python3, forward + verify.py's reversed
  order) incl. the new test_prompt_floor_banks_subexcitation gate
  (pure-Ar cascade bit-identical with the floor; Electra-mixture
  eps_sec > 3.5; subfloor 0; no pre-spent F2(a3); identity < 1e-9);
  physics-band history extended (0.0054, Te 1.9 eV). verify.py ALL
  CHECKS PASS; Homebrew 3.14/numpy 2.5 smoke clean (prompt-floor
  yields + all three modes, DeprecationWarning as error).
- Board **57 items: 56 ok / 1 warn / 0 gap** (the population-benchmark
  item is ok; the warn stays, re-scoped to the quantified residuals),
  JXA-checked.

## [0.1.11] — 2026-08-05 — the Electra ArF anchor lands; building it catches 0.1.10 double-counting

Owner ask: "Proceed with the next milestone" (the Electra ArF anchor
proposed after the anchor-inventory survey). Building the arf_beam mode
forced a Chen-block audit — and the audit FALSIFIED four 0.1.10 claims.

### FIXED — 0.1.10 ran 8 channels double-counted (found by the audit)
SYMPTOM: krf carried two rows for one physical channel on 8 channels
(R29/AF1 harpoon, R30/AIIR1, R31/AIIR2 ion-ion, R27/AAR1, R22/FMN1,
R33/ARF_F2, R34/ARF_AR quench, R37/ARF_RAD radiative — the last three
literally identical values twice). ROOT CAUSE: 0.1.10 established row
absence by GREPPING network.py source, but the CHEN block is GENERATED
from chen2025.TABLE_III — its rows are invisible to source grep. Four
0.1.10 claims were therefore WRONG ('no Ar-atom harpoon', 'the primary
beam ion had no F- path', 'source-less Ar2+', 'sink-less F2+') — the
channels existed; the additions still stand as the better-sourced
first-hand values, and the claims are corrected on the rows, in the
docstrings, on the board, and in OMISSIONS. FIX, with teeth:
- `CHEN_SUPERSEDED` (8 entries): in the Kannari-chemistry modes the
  first-hand LIT rows win and Chen's are dropped (conflicts recorded:
  AF1 7.5e-10 vs R29 4.0e-10 x1.9; Flannery turnover vs R30/R31
  constants; ARF_RAD 2.50e8 vs R37 2.381e8); arf_chen keeps Chen
  VERBATIM — that anchor scores THEIR model. Assembly asserts every
  replacement is present before dropping.
- **check_network channel-uniqueness hard error** (one (reactants,
  products) multiset = one row) — the defect class is now unwritable;
  teeth-tested. New Known-bad entry (walk the BUILT network).
- RE-MEASURED clean: v2 g0 0.0036 -> **0.0034 /cm** (closure 1e-12
  class); Kannari anchor x0.62 -> **x0.56** at R_OC 80% — the 0.1.10
  verdict ('suspect consumed, deficit stands') survives with slightly
  deeper numbers; est-sweep re-measured (Ar fill <=8.0% via URAD_Bv,
  15 structurally zero unchanged; Ne fill +-66-69% unchanged).

### THE MILESTONE — the NRL Electra ArF device anchor, scored live
Sources first-hand 2026-08-05: Wolford 2020 HEDP 36,100801 (the
quadruple-pass oscillator measurements) + the Orestes memo petrova2022
NRL/6720/MR--2022/1 (Table 2 absorbers, App. A4 optics, App. A2-A3
rates for recorded conflicts). New `arf_beam` mode: the SAME
Kannari-chemistry pool as krf, mechanically species-filtered to
ARF_MODE_SPECIES (55 rows; no hand-picked set to drift; krf dedup rule
applied; zero Kr leakage — gate-tested), with NRL's PRINTED 193-nm
optical set (sigma_se effective 1.0e-16 cm2 [their Fig. 10; B(v=0)
1.9e-16 and Lee's 2.75e-16 recorded]; Table 2 absorbers F2/Ar*/Ar**/
Ar2*/Ar2F*/Ar2+/F- incl. the new Ar2* photoionization routing).
`scorecard.electra_arf_rows`: the two folded gain lines unrolled to a
2-m cell (10.5x11 cm2 measured beam; round-trip gain path 4 m = the
real four passes), R_OC 8%, all other passive losses in L_op 0.37,
flat-top power density matched by construction (500 kV, trapezoid
20/100/20 ns). Stated limits: no beam-direct F2 channels (NRL Table
A2; their F2 runs 50% dissociated), no F2*/F* species, Ar-fill LOKI
tables carry 5% Kr.

### THE MEASURED VERDICT — the deficit goes CROSS-DEVICE
Electra ArF: **x0.07 / x0.53 / x0.69** vs measured 90/96/80 J at
0.82/1.09/1.36 atm; the measured pressure trend (peak 1.09, decline
1.36) NOT reproduced — the engine rises, and NRL's own Orestes also
overpredicts the high-pressure points (their stated 3-body suspicion =
the ARF_AA class carried here first-hand). THE HEADLINE: x0.53 at
Electra's robust point is the SAME deficit class as x0.56 at the
Kannari KrF point — a different lab, NO krypton, 193 nm, a different
cavity — and the 0-D construction extracts everything it sees, so the
true deficit is understated if anything. The board's one warn is now a
CROSS-DEVICE shared-layer diagnosis (fluid-Te treatment, 0-D
uniformity, emergent-W +5-8%, beam-partition spreads), not a
Kr-chemistry question. Anchors score, never steer — nothing was tuned.

### Verified
- **157 passed** (~91 s, /usr/bin/python3) incl. the two new gates
  (test_chen_supersede_and_channel_uniqueness,
  test_arf_beam_mode_and_electra_optics); census re-pins CHEN >= 15
  (krf 23 -> 15, the supersede), rows 177 -> 169; physics-band history
  extended (0.0034). verify.py ALL CHECKS PASS; Homebrew smoke clean;
  JXA board check green.
- Board 56 items: 55 ok / 1 warn / 0 gap (the Electra anchor item is
  ok — an anchor that reports its misses is doing its job; the warn
  stays, now cross-device).

## [0.1.10] — 2026-08-05 — the ArKr+ ion manifold lands; the suspect is consumed, the deficit stands

Owner follow-up on 0.1.9: "It says one is still an open problem" — the
board's one open warn (the formation-side deficit, leading suspect the
omitted ArKr+ ion manifold). This release implements the recorded next
milestone and MEASURES the hypothesis.

### THE MILESTONE — the ArKr+ ion manifold, adopted whole
Both primaries re-read first-hand 2026-08-05 at the table level
(Kannari 1985 Table I pp.4311-4313: CT / ion-ion / DR / Penning /
superelastic / absorption blocks; J&H 1980 Table I p.2408 #10,36-40).
`network.ar_kr_ion_manifold()` (18 rows) + the ArKr+ species (EST well
13.5 eV, NeKr+ template):
- **AKI1: Kr+ +2Ar -> ArKr+ +Ar 1.0e-31** [Kannari 16; J&H #37 same] —
  as quantified in 0.1.9, THE dominant Kr+ route in an Ar fill (3.0e8
  s^-1 at the device point vs 8.4e7 for AK7). AKI2: Ar+ +Kr+Ar ->
  ArKr+ 1.0e-31 [66/#36]. AKI3: ArKr+ +Kr -> Kr2+ +Ar 3.2e-10
  [16/#39] (the HK3 analog).
- **Both printed F- branches** (J&H #10 -> KrF*+Ar / #38 -> ArKrF*,
  1.0e-6 each): K108 Flannery form 50/50 (params REPO-mirrored,
  stated; total in-band of the prints). AKI_DR: e+ArKr+ -> Kr**+Ar
  1.2e-6 (Te/300K)^-0.55 [71/#40] on the ADR 300 K-total convention.
  AKI_ION/AKI_SE/AKF_SE close ArKr*'s printed fate set [72/p.4312].
- **The same walk exposed and fixed:** Ar+/Ar2+ +F- recombination
  MISSING (AIIR1/AIIR2 — in an Ar fill the PRIMARY beam ion had no F-
  path; K26/K108 mirrors, existence printed [Kannari 28-35; J&H
  #5/#9]); **Ar2+ SOURCE-LESS since 0.1.0** (AAR1: Ar+ +2Ar 2.5e-31
  [65] — its consumers AK2/ADR and its absorber never fired); **F2+
  SINK-LESS since 0.1.0** (FMN1: F2+ +F- -> 3F 4.0e-8 [91]); the five
  printed Penning rows (PEN_* 5.0e-10 [86,87] — Ne had NR4, Ar/Kr had
  nothing).

### THE RE-WALK RESIDUE — 0.1.9's own misses, caught (`ar_fill_completion`, 9 rows)
The Known-bad rule applied to itself: the fresh row-by-row walk caught
neutral rows the 0.1.9 adoption skipped —
- **The Ar*/Ar** + F2 harpoon quartet (AF1-4)**: krf mode had NO
  Ar-atom harpoon at all (ArF* formed only via Ar2*). AF1 7.5e-10
  [52,59,92; J&H #4 8.5e-10 recorded] is ~40 % of the Ar* fate budget
  at v2. With AK1 (ArF*+Kr -> KrF* 1.6e-9) this is direct formation.
- **The ArF* relaxation block (ARF_AA/F2/AR/RAD)**: AK1 was ArF*'s
  ONLY consumer — 100 % reached KrF*. The printed block (4.0e-31
  three-body [62,63,90], 1.9e-9 F2, 9.0e-12 Ar, 2.50e8 rad [47,48])
  takes ~30 % of the fate budget. Kannari's ArF prints are footnote-a
  reduced by 0.418 ALREADY — consumed as printed (the 0.433-trap rule,
  ArF side).
- **AKD5**: ArKr*+F2 -> ArKrF*+F 3.0e-10 ('analogous to Ar2*') — the
  THIRD ArKr* row; 0.1.9 adopted two of three.
- Absorbers: ArKr* photoionization 1.5e-17 cm2 [85] (needs ArKr+ to
  exist — hence adoptable only now); Ar**/Kr** 6.0e-18 [84] (routing
  existed, sigmas never adopted); F2 EST -> printed 1.2e-20 [81].
  **CORRECTION: the 0.1.9 OMISSIONS 'ArKr+ photoabsorption' was a
  species TYPO — the printed absorber is the NEUTRAL ArKr*; neither
  primary prints a sigma for the ion, so none was invented.**

### Value re-pins by two-primary consensus (K15 precedent)
- **AK7 3.25e-31 -> 2.5e-31, REPO->LIT** (Kannari [43,67] + J&H #43
  [32,40] both print 2.5e-31; Chow's 3.25e-31 recorded in-spread).
- **AK2 7.8e-10 -> 7.5e-10** (Kannari [64] + J&H #35 [31]; Boichenko's
  7.8e-10 x1.04 recorded).
- **R11p re-lumped, REPO->LIT**: products Kr2+ +e -> **Kr*+Kr** (the
  dominant printed branch, was Kr**) at the 300 K two-branch total
  2.24e-13 on the Te^-1/2 template (the ADR convention; the 0.1.0
  'citation debt' is paid; J&H's single Kr**-product row recorded).

### THE MEASURED VERDICT — hypothesis consumed, deficit stands
Kannari device anchor: **x0.60 -> x0.62** at R_OC 80 % (0.37 vs 0.6 J;
20 % still threshold-pinned ~0); g0/alpha 2.6-2.8 vs the empirical
4-12 band. v2 g0 **0.0053 -> 0.0036 /cm** (ne_pk 5.6e20, Te 0.52 eV,
closure 3.7e-13). The formation gains are OFFSET by sinks the same
tables print: the ArF* relaxation block diverts ~30 % of ArF* off
AK1's transfer, AAR1 turns on the previously-dead Ar2+ 248-nm
absorber, AKD5 diverts ArKr*. **The leading-suspect hypothesis is
consumed** — with both primaries' neutral+ion blocks now FULLY
disposed row-by-row, the residual deficit points OFF-network; the
board warn stays OPEN with the suspect list updated (fluid-Te
treatment, 0-D uniformity, emergent-W +5-8 %, the K67 e-quench era
spread). Anchors score, never steer — nothing was tuned to move x0.62
toward 1.
Recorded-not-adopted (quantified in OMISSIONS): e-ionization of
excited/dimer states [72,79] except AKI_ION and the Ar2*/Kr2*/Ar2F*/
ArF* superelastics (each <=1 % of its species' fate budget at beam
points), F+F+M -> F2+M [16] (~1e3 s^-1 on pulse timescales).

### Verified
- **155 passed** (~85 s, /usr/bin/python3 3.9.6/numpy 2.0.2) incl. the
  new `test_ar_ion_manifold_complete_and_sourced` gate (species wiring,
  AKI1 decade-pin, 23 LIT + 5 REPO rids, 50/50-sums-to-K108, absorber
  routing + energetics, no-invented-ion-sigma, the three re-pins);
  census re-pinned LIT >= 100 (measured 103 of 177; EST still 19 — S5
  pricing intact); physics-band floor 0.003 -> 0.002 with history
  (measured 0.0036). `verify.py` ALL CHECKS PASS; Homebrew 3.14/numpy
  2.5 smoke clean (DeprecationWarning as error; same 177 rows, same
  g0 0.00364, closure 3.2e-13 — two-interpreter agreement).
- Full est-sweep RE-MEASURED on the 0.1.10 baseline (the 0.1.9 priced
  numbers were against the old network): Ar fill 15 structurally zero
  (same set), 4 active levers <= 7.4 % g0 (URAD_Bv x3; the KR2_RAD
  LIT-conflict alternate +8.9 %); Ne fill BC/CB_Ne +-66-69 %, NK6
  +-40 % — the S5 pricing holds, board metric refreshed.
- Board **55 items: 54 ok / 1 warn / 0 gap** (the milestone item is ok;
  the warn stays, honestly re-scoped), JXA-checked.

## [0.1.9] — 2026-08-04 — the last two warns closed; the Xu pass finds the missing Ar manifold; the anchor flips

Owner ask: "there are still two outstanding warnings in the code, can you
help me finish and verify the full software structure?"

### L3 CLOSED — beam mode consumes the measured relaxation curves
`kinetic_tier.RelaxationBridged` (engine default for the beam drive;
`beam_subtable="taper"` reproduces pre-0.1.9): below the swarm-table
floor, el:Ar/el:F2/el:Kr (down to 1.26 eV), vib and F- (down to 1.76 eV)
are now the MEASUREMENT itself — log-log interpolation on the shipped
relax curves, a 0.1-ln crossfade for BDF continuity (narrower than the
deepest bin's distance to the floor, so measured bins are consumed
pure), bins with >=50 % stat error excluded (stated rule: the 100 %-err
vib bin at 1.26 eV and every sub-floor Kr+ bin). The x3.3-high F- taper
value at 1.76 eV is gone (junction ratio -> 1.0 exactly at the bin).
Honest residual, now SELF-REPORTED per run (new validity gauge +
regime-line phrase): the v2 point's <eps> (0.3-0.75 eV) sits below even
the measured floor — the exporter's relaxation sampling stalled at the
Ar RAMSAUER MINIMUM — so there the taper still extrapolates, from the
measured anchor instead of the table floor. Deepening the export past
the stall (same seed => pure extension, no pinned band moves) is the
recorded follow-up (OMISSIONS). Measured v2 impact: ne_pk +2.4 %,
Te_end -4 %, g0/eta unchanged. Circuit path exempt (E/N-direct). Gates:
4 bridge unit gates + the Device-level default/off-switch gate.

### THE HEADLINE — the Xu 1993 pass finds the missing Ar heteronuclear manifold
Xu/Gadomski/Setser JCP 99,2591 (1993) read first-hand at the table level
(pp. 2591-2597) to close the recorded K34Ar follow-up — and its p.2595
citation of the KrF*+2Ar three-body ((2-8)e-32 cm6/s "because of ArKrF*
formation") exposed that the ENTIRE ArKr chain was absent here. Kannari
1985 Table I (pp. 4311-4313 re-read) carries every row; J&H 1980 Table I
p.2408 double-sources them (their Fig. 2 literally draws the ArKr*/
ArKrF* boxes). Root cause: the REPO source network (spectre) is
Ne-buffer-oriented, and earlier table passes filtered by species already
tracked — rows of UNTRACKED species were invisible. Added
(`network.ar_kr_manifold`, 30 rows post-split, 2 new species, all LIT,
every value as printed with inter-source conflicts recorded):
- **KAA: KrF*+2Ar -> ArKrF*+Ar 7.0e-32** [Kannari 62,63; J&H #14
  8.0e-32; Xu (2-8)e-32] — ~3.4e8 s^-1 per state at v2, x2.4 the B
  radiative rate: the largest single defect in the Ar-fill answer, x2
  the 0.1.8 K34Ar headline channel.
- **KQAr: KrF*+Ar -> Kr+Ar+F 1.8e-12** (Kannari-estimated; three-source
  record: J&H #13 prints 0.0, Xu measures <=(2.9+-0.5)e-12 as an
  impurity-limited upper bound).
- ArKr*: AKD1/AKD2 formation (1.0e-32 [52]), AKD3 harpoon (6.0e-10
  [52]; J&H 1.0e-9 recorded), AKD4 +Kr -> Kr2* (1.0e-10, dominant),
  AKD_RAD (8.0e7 Kannari-est; J&H 3.0e6 — the x27 triplet-vs-mixed
  spread RECORDED and swept).
- ArKrF* fates: +F2 1.0e-9 [52/#19], +Kr -> Kr2F* 2.0e-11 [43/#16],
  +Ar -> Ar2F* 2.0e-11 [43/#17] (dominant, 1.4e9 s^-1), rad 5.0e6
  (Kr2F*-class; J&H's 5.0e7 recorded as typography-suspect).
- Ar2F* fate block (it was a SINK-LESS ABSORBER since 0.1.0): +F2
  2.0e-10 [61], +Kr -> ArKrF* 1.0e-10 [52,59/#26], +Kr -> KrF* 1.0e-10
  (est), +Ar -> 3Ar+F 1.0e-10 (est; 7.0e9 s^-1 — a flash-through
  species, which is why the omission stayed invisible), rad 5.4e6
  [61,93].
- Kr2*/Ar2* VUV decay (KR2_RAD 8.0e7 / AR2_RAD 6.0e7 [88]) — the 0.1.0
  "spectre parity" scope hole, closed; J&H print 3.3e6/3.8e6 (x24/x16
  spreads recorded and swept). Second harpoon branches K2F_BR/A2_BR
  (-> Kr2F*/Ar2F*, [52] + J&H #20/#23), F-atom partners K2F_F/A2_F,
  AK9 (Ar2*+Kr -> Kr* 4.0e-10 [43/#33] — x20 every previously-present
  Ar2* channel), ADR (Ar2+ DR, absent while Kr2+ had R11p).

### Xu-pass row corrections (dedicated-later-measurement precedent, AK1)
- **K33 2.0e-12 -> 4.0e-12** (Xu Table I/p.2596 recommended
  (4.0+-0.5)e-12 MEASURED; Kannari's 'estimated' 2.0e-12 and J&H #12's
  2.0e-11 recorded).
- **K169 5.6e6 -> 4.5e6 s^-1** (Xu p.2597: tau(Kr2F*) "established as
  220+-20 ns to within 10 %"; Quigley-era 181 ns = Kannari's print,
  Ueda 210, Setser-prev 240+-35, J&H 6.7e6 — all recorded).
- **K34Ar follow-up CLOSED**: Xu Table II mixed Ar+Kr (5.5+-2.0)e-31 —
  6.5e-31 in-band; now triple-sourced. K34 REPO->LIT (Table II 2Kr
  (6.0+-2.0)e-31 verified; Kannari 6.7e-31 in-band), K32 REPO->LIT
  (Table I (3.0+-0.3)e-10 "strongly support"; Kannari's 7.8e-10 = the
  superseded era, recorded), K68Ne/NK5 fetches -> first-hand (NK5 +
  Xu's upper-limit caveat), K125 + Xu's (1.2+-0.2)e-10 in-spread,
  BC_Ne/CB_Ne + Xu's "may be much smaller" bracket direction, AK8
  1.0e-32 -> Kannari's 2.7e-32 [52] (REPO->LIT, conflict recorded).
  johnson1989.pdf checked: XeCl model — not applicable to the KrF EST
  set (recorded).

### THE MEASURED CONSEQUENCE — the anchor flips, and the board says so
v2 point: g0 0.0278 -> **0.0053 /cm** — at 0.42 MW/cm3 deposition the
completed network is SUB-THRESHOLD (Kr2F* peaks 3.2e20 m^-3; its
6.0e-18 cm2 absorption alone is a third of the remaining gain). Sitting
FAR below spectre/phantom is the signature of channels they omit. The
**Kannari device anchor now scores LIVE** (`scorecard.
kannari_device_rows` — the hand-text block is retired): deposited 66 J
by construction, **x0.60 UNDER at the robust R_OC 80 % point,
threshold-pinned ~0 at 20 %** (g0 0.037 vs g_th 0.030 at his L_op
0.64); g0/alpha 2.4-2.9, BELOW the empirical 4-12 band. The old
x1.31-1.34 agreement was CANCELLATION: missing sinks vs a formation-side
deficit — the deficit now shows naked, and the board carries it as the
one open [warn]. Leading suspect, quantified in OMISSIONS: the ArKr+
ion manifold (Kr+ +2Ar -> ArKr+ 1.0e-31 = the dominant Kr+ route,
3.1e8 s^-1 at his fill; F-/e recombination feeds KrF*/Kr**; 1.5e-17
photoabsorption). Corroborating asymmetry: Kannari's model matches his
device carrying a x4.4 LARGER KrF* e-quench — his formation is
correspondingly stronger.

### S5 CLOSED AS PRICED — the x1/3..x3 bracket sweep, implemented + measured
`scorecard.est_bracket_rows` / `run.py --est-sweep` (the engine itself
has NO rate-scaling pathway — the sweep rebuilds a CompiledNetwork per
variant, D4 intact; "stepwise" specs gained an explicit optional
prefactor to be bracketable). Element-bound classification is
structural (a row needing an element the fill lacks runs on 1e12 seed
vs ~1e26 — 1e-8 class). MEASURED (120 ns, 0-D):
- **Ar fill: 15 of 19 EST rows structurally zero; the 4 active levers
  (R8, KSTEP, KPS_Kr, URAD_Bv) move g0 at most -5.6 %** (URAD_Bv x3).
  The only >10 % lever is not EST: KR2_RAD at J&H's alternate print
  (+15.5 %) — a recorded LIT-vs-LIT spread, swept via
  CONFLICT_ALTERNATES.
- **Ne fill: the unmeasured Ne-partner B/C mixing analogs are a REAL
  +-70 % g0 systematic** (BC_Ne/CB_Ne; NK6 +-40 %, NK4 -17 %, NK8
  -12.5 %, VREL_Ne -14 %), direction-informed by Xu p.2596 (toward
  +70 %). Named on the board as the Ne-mode systematic.
Census: 118 -> **149 rows; LIT 40 -> 78; EST 19 (unchanged — every new
row arrived cited); REPO 22 -> 19** (the three upgrades). Literature
exhausted: 7 in-repo papers mined first-hand.

### Also fixed
- Stale scorecard OMISSIONS (in-loop ASE and gated transport were
  still listed as absent — closed 0.1.6), stale flow-board "x100 vib
  unresolved" note (resolved 0.1.4), stale board metrics (v2 g0,
  census, validity numbers).

### Verified
- **154 passed** (~75 s, /usr/bin/python3 3.9.6/numpy 2.0.2) including
  8 new gates (4 bridge + default/off-switch + manifold-complete +
  bracket-machinery + spot-bracket); `verify.py` ALL CHECKS PASS;
  Homebrew 3.14/numpy 2.5 smoke clean (DeprecationWarning as error).
- Ledger closure unchanged (2.9e-13 v2; the manifold adds no closure
  drift — every row balance-checked mechanically, teeth intact).
- Board **54 items: 53 ok / 1 warn / 0 gap** (the two 0.1.8 warns
  CLOSED; the ONE new warn is the exposed formation deficit), JXA-
  checked.

## [0.1.8] — 2026-08-04 — the last three warns: J&H pass, the L3 measurement, W(Ne) decomposed

Owner ask: "Proceed with those last three remaining open issues."

### THE HEADLINE — a missing dominant channel, and the anchor at x1.31
J&H 1980 Table I read first-hand (Johnson.pdf p.2408): **#18 lists
KrF* + Kr + Ar -> Kr2F* + Ar at 6.5e-31 cm6/s** (Kannari's Table I
carries the same value — both first-hand). The Ar-third-body KrF*
quench was MISSING from this network (and from phantom AND spectre —
their v2-point g0 agreement was among models sharing the omission)
while the Kr- and Ne-partner rows existed. At the v2 point it is a
1.7e8 s^-1 sink, comparable to the B radiative rate. Added as K34Ar
(LIT, double-sourced; Xu 1993 Table II re-verification is the recorded
follow-up). Measured impact: v2 point g0 0.0451 -> **0.0278 /cm**, eta
12.7 -> **2.2 %** (now BELOW spectre/phantom, as expected — they lack
the channel); **the Kannari device anchor moves x2.3–3.1 ->
x1.31–1.34** at the working R_OC points (x2.0 near threshold) — the
remaining gap is the size of Kannari's own quantified one-state/v-relax
effect (25–50 %), the known unimplemented remainder.

### S5 — the J&H pass (plus corrections)
- **K33 provenance corrected**: phantom attributed 1.0e-12 to "J&H
  Table 1" — J&H #12 actually prints **2.0e-11** (a x20 mis-cite).
  Adopted Kannari's published 2.0e-12 ('estimated' there) with the
  conflict recorded; value x2.
- **K15 -> 2.5e-31** (LIT): three first-hand sources agree (J&H #41,
  Kannari, Boichenko #18 with the Tg-in-eV conversion); Chow's 2.7e-31
  noted in-spread.
- **AK1 x5.3 spread recorded** (J&H #2: 3.0e-10 vs Kannari [62]
  1.6e-9; the later dedicated measurement kept). KSTEP corroborated
  (J&H #48 discharge-typical 1.5e-8 vs the template's 1.4e-8 at
  1.5 eV). Census: **118 rows, LIT 40, EST 19** — the remaining 19
  have no published measurement in the five papers read.

### L3 — the taper's price tag, MEASURED
The exporter gains ``relaxation_sampling``: field-free relaxing
populations through PIC_MCC_3D's operator, per-process rates binned by
instantaneous <eps> down to 0.5 eV — the beam-mode EEDF class, BELOW
the swarm tables where the engine's flat+taper stand-in operates.
Verdict (shipped in the export, pinned by test, on the scorecard):
**elastic/vib within x0.8–1.9 of measurement; F- attachment
OVERESTIMATED x3.3 by the taper at <eps> = 1.8 eV** (7 % stat). The
circuit path is exempt (E/N-direct). Consuming the relaxation curves
directly in beam mode is the recorded follow-up.

### W(Ne) — decomposed, closed as characterized
Both W routes agree to 0.02 %, so the +16–19 % vs literature is data/
method: decomposition shows our excitation + sub-threshold budget is
22.1 eV vs the 15.0 eV the literature value implies (U_exc/U_ion =
0.62) — carried by the Biagi set's excitation/ionization ratio and the
equal-sharing model (uniform sharing moves W −5.5 %, still +14 %). On
the scorecard as a stated, explained miss.

### Verified
- **147 passed** (~62 s); board **52 items: 50 ok / 2 warn / 0 gap** —
  the two remaining warns are the honest irreducibles (19 EST analogs;
  the taper with its measured band), each carrying its number.

## [0.1.7] — 2026-08-04 — open modules: the Lee/Boichenko EST pass + validity gauges

Owner ask: "Proceed with the open modules" — the remaining [warn] items.

### S5 — the EST retirement pass (Lee 1989 + Boichenko 2000, first-hand)
Both PDFs read at the table level (lee1989.pdf pp. 2055-2056;
Boichenko_2000.pdf pp. 681-683, with the Tg-IN-EV convention decoded:
their Tg^-0.75 at 300 K is x15.5, reconciling Table 1 #16-18 with their
own text and with Kannari/Chow class values).

RETIRED EST -> LIT (direct citations of the actual process):
- **HK4** — Boichenko Table 1 #8: NeKr+ +e -> Kr** +Ne, 1e-7 Te^-1/2.
- **HK1** — Boichenko Table 1 #23: Ne+ +Kr+Ne -> NeKr+ +Ne, 1e-31
  (marked assumed-universal there; Lee [40] analog corroborates).
- **NK2** — Boichenko p.682 bound (<= 5e-13) re-verified; value AT the
  bound. Recorded conflict: Lee [43]'s Ne2+ +Ar 3.0e-11 vs Boichenko's
  Ar bound — the direct Kr bound is kept.
- **APS_Ne** — Lee [37]: Ar** +Ne -> Ar* +Ne 1.0e-10, direct.
Also REPO/EST -> LIT with digits verified: NF1/NF2 (Lee [49] pressure
forms digit-exact), NR3 (Lee [35]), AK2 (Boichenko p.682); APS_Ar gains
Lee [36] as an independent second source; analog digits verified (rows
stay EST — the Kr processes are unmeasured): NK1, NK1b, NK8 (Lee [38]),
HK2 (1.1e-9), HK3 (1.2e-11) (Lee [44]), HK5 (Lee [49] digit-exact).
Census: **EST 23 -> 19, LIT 23 -> 33, network 113 -> 115 rows**.

### FOUND — two channels Lee's table DOES list (contradicting a record)
- **LNF1/LNF2 added**: Ne* +F2 -> NeF* +F (5.3e-10, Lee [32]) and
  Ne** +F2 -> NeF* +F (7.5e-10, [54]) — Lee 1989 Table I p.2056
  "Neutral Formation" lists both, CONTRADICTING the recorded spectre
  2026-07-14 adjudication ("the row does not exist on the page") that
  shaped nef_channels() to ion-ion only. In a Ne fill at 0.1 % F2 this
  is a ~4e7 s^-1 Ne* drain competing with the Kr Penning chain (~14 %
  at 5 % Kr). Flagged for spectre's owner; not edited there.

### FOUND — a third Chen-table typography flag (R36, x100)
Chen R36 (ArF* +Ar+Ne -> Ar2F* +Ne) prints 3.5e-33 cm6/s; Lee [58]/[59]
give the ArF* three-body class at 5.0e-31/1.2e-31 and Kannari Table I
4.0e-31 — both first-hand, both ~x100 above. Recorded in chen2025.py's
transcription-flags block (kept verbatim in the Chen-verbatim scoring
mode by definition); a candidate term in the x3-5 per-joule formation
gap. Corroborations logged: R24 (Lee 4.1e-34), R27 (Lee/Boichenko
2.5e-31 exact), R28 (Lee [39] 4.4e-32 exact), R33/R34/R35/R37 all
confirmed in Lee's table.

### L3/A3 — validity gauges (the ESSim.stability analog)
Every run now reports result["validity"] on its regime line: peak
ionization fraction (the static-N gauge), F2 burn-down, and the fraction
of the pulse with <eps> inside the loaded rate tables (speaks directly
to the <eps>-mapping [warn]). Gated.

### Verified
- **146 passed** (~63 s); board 52 items (49 ok / 3 warn / 0 gap),
  JXA-checked; Homebrew smoke clean. The 3 remaining [warn] are the
  honest opens: 19 EST rows, the <eps> mapping on non-circuit paths,
  the W(Ne) data/method miss.

## [0.1.6] — 2026-08-04 — the four missing modules: board reaches 0 [gap]

Owner ask: "Proceed now with the missing modules" — the four [gap] rows.

### Added
- **In-loop ASE sink (S6 closure)** — `engine._ase_rate` prices the 4-pi
  ASE flux of the CURRENT gain field every RHS call (closed-form chords,
  48 precomputed directions, volume-averaged g0/alpha — approximations
  stated in the docstring; ~1e-5 axial overlap with the cavity mode
  noted); `chemistry.rhs(ase_rate=...)` applies it with stimulated-
  emission routing except the photons ESCAPE. First-principles, no new
  constants. `ase_inloop` param (default ON; off reproduces pre-0.1.6).
  Measured: v2 point g0 0.0458 -> 0.0451 /cm, eta 12.81 -> 12.72 %
  (modest, correct at that gain); direction gate at Kannari-class gain,
  1e-6 on/off equivalence below threshold, every-eV identity intact.
- **Inter-cell transport (gated)** — conservative 6-neighbour species
  diffusion (`transport` + `diff_m2_s` params, EST bracket, Tg/p
  scaling), with neighbour blocks added to the Jacobian sparsity.
  Machine-exact global conservation; sinusoid eigenmode matches the
  discrete -D k^2 to 1e-6. DEFAULT OFF with the measured justification:
  tau_diff ~ 0.1 s at device presets — the former scope note is now a
  number plus an implementation.
- **Mesh refinement, the MESH.txt rung** — `beam.equal_deposition_edges`
  places beam-axis cell edges at the n-quantiles of a fine-march
  cumulative deposition; `Device(axis_edges=...)` runs on the refined
  axis (volumes exact, deposition marched per-width, ledger closes).
  Measured on a fully-stopped 100 kV / 6 cm / 3 atm case: per-cell
  deposit spread x2e11 (uniform grid; empty cells past the stop) ->
  x3.5 (refined, mm cells at the Bragg stop); totals match the fine
  march to 1 %. Static refinement — dynamic re-refinement stays future
  work, stated. Guards: beam-mode only, no transport combination,
  monotone-span validation.
- **Chen-figure scoring (S7, unblocked by the 0.1.5 circuit)** —
  `scorecard.chen_figure6_rows` RUNS the circuit preset (~3 s) and
  scores deposited energy, ne, n_ph, P_L and the E/N class against
  Fig. 6 with ratios and the x3-5 per-joule verdict; a test pins the
  recorded ne miss so it cannot drift silently.

### Verified
- **145 passed** (~62 s, /usr/bin/python3); board **51 items:
  48 ok / 3 warn / 0 gap** — the first zero-[gap] state, JXA-checked.
  Remaining [warn]: EST 23 (S5), the <eps> mapping on non-circuit paths
  (L3), the W(Ne) data/method miss.

## [0.1.5] — 2026-08-04 — the CIRCUIT discharge drive (E/N-direct, self-consistent breakdown)

Owner ask: "Proceed with the next issues" — the top declared item: the
circuit drive that owned the anchored Chen gap.

### Added
- **`MODULES/circuit.py` + the `discharge_mode="circuit"` drive**: Chen
  2025's C0/Cp/L chain (10.4 nF / 8.1 nF / post-switch 1 uH / 30 nH,
  declared parameters) as FIVE engine states (q0, qp, i1, i2, E_gap)
  under the same BDF; conduction j = e n_e v_drift(E/N) from the swarm
  table's BulkV_z; current continuity (conduction + displacement) closes
  the gap field. **Electron-impact rates on this path are looked up by
  E/N DIRECTLY** — the tables' native index; the <eps> = (3/2) Te mapping
  is entirely absent from the circuit-driven path (`LokiRates.k_by_en`,
  `drift_velocity`, `mean_energy_at_en`; `CompiledNetwork` passthrough).
  The ledger's drive_in integrates j.E exactly (the circuit's own stores
  are upstream of the cell budget and reported separately). Prescribed
  mode unchanged and still available. `arf_chen_circuit` preset.
- `tests/test_circuit.py` (6 gates): LC energy conservation of the bare
  circuit (< 1e-3 over 150 ns at R = 0 — the analytic-oscillator gate),
  the 15-Td drift-velocity pin (BulkV_z 5.422e4 m/s), E/N-direct vs
  <eps>-mapped consistency (< x2 — same table, two indexings), the Chen
  base-case class run, 0-D-only enforcement, prescribed-mode isolation.

### Measured — the decisive Chen re-anchor
Chen base case (-18 kV, 3.4 atm, 12 mm x 90 cm2) through the real
circuit: **1.11 J of the 1.68 J stored deposits self-consistently**
(0.39 J left ringing, remainder in the declared loop resistance), U_gap
peaks 14.6 kV and E/N 14.6 Td — their waveform class, **100 % in-table**
— the avalanche runs 1e14 -> 1.7e20 m^-3, ledger closes to 7e-6, wall
3 s. Chen reaches ~1e21 and 5 MW from the same joules. **The drive
pathway is thereby ELIMINATED as the anchored gap's cause: what remains
is a x3-5 PER-JOULE formation gap in the ArF network at ~15 Td**
(suspects, in order: the Ne* production-rate source (our curated LoKI
mixture sweep vs their BOLSIG+/SIGLO), the omitted Xe seed, 0-D vs their
1-D). Scorecard and board carry the sharpened verdict.

### Scope
Kannari-style magnetic switch simplified to its post-switch inductance
(stated); circuit mode is 0-D in v1 (enforced with a clear error); Xe
rows remain omitted.

### Verified
- **137 passed** (~48 s, /usr/bin/python3); board 51 items
  (44 ok / 3 warn / 4 gap), JXA-checked; Homebrew 3.14/numpy 2.5 circuit
  smoke clean (DeprecationWarning as error).

## [0.1.4] — 2026-08-04 — open problems: B/C/Bv split, the vib-bucket BUG, small closures

Owner ask: "Proceed with fixing the open problems." Attacked in anchored-
value order; the circuit drive, in-loop ASE, transport, AMR and the
remaining EST rows stay declared roadmap items.

### Added — the Kannari upper-level split (the x2-3 anchor's leading cause)
- **KrF*(B)/KrF*(C)/KrF*(Bv) replace the lumped KrF*** in the KrF mode:
  a single mechanical transform (`network.split_krf_bc`) rewrites every
  lumped row (consumers -> one row per state at the same k; formation ->
  0.5 Bv + 0.5 C, nascent KrF* born hot per Kannari [41]), plus first-hand
  Table-I rows: URAD_B 1.43e8 / URAD_C 1.33e7 s^-1, B<->C mixing 5e-10
  with detailed balance 0.682 (DE = 80 cm^-1), Bv->B vibrational
  relaxation 4e-11. Only B lases: sigma_B = 2.6e-16 cm2 [ref 46],
  tau_B = 6.99 ns. Network: 89 -> **113 rows** (Ar fill; Orestes: 116).
- **The B/C equilibrium is a real inter-source spread, NOT resolved by the
  anchor** (PLAN (26) forbids picking by it): Kannari's constants give
  f_B = 0.405 (his own tau_eff 15.18 ns derivation reproduces it); the
  retired lumped pair encoded f_B = 0.206 (Xu/Gadomski-era KEQ). The
  Kannari device anchor BRACKETS: x2.1 (0.206) to x2.6 (0.405, shipped)
  at R_OC 20% vs the measured ~2.3 J — with the intra-B vibrational
  ACCESSIBILITY (his 25-50% mechanism, v-resolved B(v<=4) with Boltzmann
  factors) the named unimplemented remainder. v2 point: g0 0.0458 /cm,
  eta 12.82%.

### Fixed — a real bug the S4 comparison caught (the "vib finding")
- **loki_tier._bucket let ANY unmatched e+F2 row fall through to 'vib'**,
  which swallowed the F2 EFFECTIVE momentum-transfer row (~1.8e-13 m3/s,
  elastic-scale, x94–192 the true vib rate by a three-way check against
  the Maxwellian average of the raw sigma). The engine drained 0.2 eV per
  FAKE event in every in-table-Te run (discharge mode; beam mode was
  shielded by the low-Te taper). Effective rows are now dropped for every
  gas (the true Elastic row is tabulated separately) and 'vib' matches
  only F2(V1..V4). Post-fix tier agreement: **vib x0.84 [0.51, 1.15]** —
  every bucket now in the pic1d x1.6 class. **The same defect exists
  upstream in projects/spectre `loki_electron._bucket` — flagged for its
  owner here, deliberately not edited there.** Chen capacitor-matched
  point: ne_pk 3.1e20 -> 4.4e20 (+40%), still sub-threshold — the
  stand-in-drive x3-7 gap stands as the dominant term.
- **Mixture validation** (the PIC_MCC_3D GasMixture rule): a typo'd gas
  name used to vanish silently from the fill; Device now refuses unknown/
  non-ground species and fraction sets that are not mixtures (gated).

### Characterised (report-only)
- **W vs ionization-sharing model** (MC cascade, PIC_MCC_3D operator):
  equal -> uniform sharing moves W(Ne) 44.0 -> 41.6 and W(Ar) 29.6 -> 27.9
  at E0 — toward literature in both, Ne still +14%: the W(Ne) miss is
  confirmed data/method class (Biagi-set coverage + sharing model), not an
  implementation defect. Equal sharing stays the default (matches the SF
  route's assumption; changing it by the anchor would be steering).

### Verified
- Suite **130 passed** (~40 s) under /usr/bin/python3; board 51 items
  (43 ok / 4 warn / 4 gap), JXA-checked. Census: EST 18 -> 23 (documented
  ratchet move: the split's flagged Ne-partner analogs + URAD_Bv + the
  doubled EST formation row), LIT 12 -> 23.
- Analytic trajectory gates re-anchored post-split (URAD_B exponential at
  tau_B; the second-order gate moved to AK6, Kannari-verified first-hand).

## [0.1.3] — 2026-08-04 — validity pass vs primary literature (first-hand)

Owner ask: "test the current validity of the REBUILD against available
data/literature to cross check." Verification policy applied: falsify, not
skim — every claim below is confirmed / wrong / could-not-check against the
in-repo primary PDFs (pages read, not abstracts recalled).

### Verdicts — Kannari 1985 (JAP 57, 4309; Table I read at pp. 4311–4313)

CONFIRMED (14): the PLAN (72) p→s relaxation value — Kr**+Ar → Kr*+Ar AND
Ar**+Ar → Ar*+Ar both 1.0e-10 cm3/s, both marked "estimated" by Kannari
himself (the previously flagged not-yet-verified item; KPS_Ar/APS_Ar now
LIT); K25/K27 harpoon 8.1e-10 [52,59,93]; AK1 1.6e-9 [62]; AK4 6.2e-12;
AK5 1.14e-32 [52]; AK6 3.0e-10 [61]; K_2Fe 1.0e-7 ("estimated" there);
F- 248 nm sigma 5.6e-18 cm2 [82]; sigma_B 2.6e-16 cm2 [46]; W(Ar) 26.2 /
W(Kr) 24.3 / W(F2) 32 eV per pair; beam f_exc = 0.28 x ionization; the
tau_B 6.99 / tau_C 75.19 ns two-state pair (the recorded v2 item).

WRONG, FIXED (4):
- **K169** 5.4e6 → **5.6e6 s^-1**: phantom's 5.4e6 is Kannari's **Ar2F***
  radiative rate — a Kr/Ar TRANSPOSITION in the original extraction.
- **K125** 2.0e-10 → **1.5e-10 cm3/s** [89,94]: same transposition class
  (2.0e-10 is Ar2F*+F2).
- **NK7** (Kr2*+F2 → KrF*+Kr+F) EST 8.1e-10 → **LIT 3.0e-10** [61] (x0.37).
- **AK3** (Ar+ +Kr → Kr+ +Ar) EST 5e-10 → **LIT 3.0e-11** [64] (x0.06 —
  non-resonant atomic CT is slow; the "J&H-class" guess was x17 high).

ABSORBER INVENTORY (Table I, 248 nm, ref 16): **Kr2F* sigma was x27 LOW**
(2.2e-19 → 6.0e-18 cm2); **Ar2F* (1.0e-18), Ar* (4.5e-18, photoionization)
and Kr* (1.0e-19) were MISSING** — all four added with PHOTO_PRODUCTS
routing (energetically valid at 248 nm). Kept, documented: thermal-scaled
Kr2+/Ar2+ sigmas (Kannari uses the hot values — x6 larger; ours carry the
repo's thermal adjudication); Kushner-89 K67 vs Kannari's 2e-7 e-quench
(deliberate later-literature choice). NOT addable: ArKr+ (1.5e-17 cm2) —
species untracked, now named in OMISSIONS via the anchor decomposition.

### Device anchors (scored, nothing tuned)

- **Kannari's own device** (660 kV, 1.88 MW/cm3 x 50 ns, Ar/Kr/F2
  89.8/10/0.2, 2.5 atm, CaF2 windows T = 0.60 → L_op 0.64; deposited
  66 J matched by construction; j required 301 A/cm2 vs diode 188 —
  the volume/backscatter uncertainty band): output OVERPREDICTED
  **x2.1 (R_OC 20%: 4.75 vs ~2.3 J) to x2.8 (80%: 1.69 vs ~0.6 J)**,
  R_OC trend correct. Leading named cause: this engine is a ONE-STATE
  lumped-KrF* model, and Kannari's paper exists to show one-state models
  overpredict extraction by **25–50%** (finite vibrational relaxation
  4e-11 + B/C mixing 5e-10 cm3/s — the recorded two-state v2 item, now
  anchored). Remainder: ArKr+ absorption, hot-ion sigmas, quench-source
  spread, geometry uncertainty. **g0/alpha = 5.6–6.5, INSIDE the
  empirical 4–12 band** (spectre 6.4, phantom11 6.6). Adding the Table-I
  absorbers moved the anchor x2.7–3.7 → x2.1–2.8.
- **Chen 2025** (Phys. Plasmas 32, 113502; pp. 113502-3..6 read):
  Table II parameters and the photon equations (7)/(8) — including the
  ln R_OC + ln(1−L_op) loss split this engine ledgers — CONFIRMED
  verbatim. **Found and fixed: the arf_chen_0d preset geometry read the
  90 cm2 discharge area as the beam cross-section — volume x75 too large**
  (right: 108 cm3, aperture 1.2 cm2, l_g 90 cm). Deposition-matched
  comparison: at their ~1.7 J capacitor energy the engine stays
  sub-threshold (ne x3 low, no pulse); raising the prescribed drive to
  11.5 J reproduces their ENTIRE pulse signature (P_L 5.16 vs 5 MW, FWHM
  15.0 vs 15 ns, n_ph 1.1e20 vs 1e20, ~90 mJ out) — i.e. the
  prescribed-power stand-in carries a **x3–7 formation-efficiency gap**
  (no E/N breakdown kinetics, no Xe seed; their true deposited energy is
  not stated, hence the range). The circuit-drive roadmap item is now a
  measured anchor, on the scorecard.
- Re-confirmed unchanged: W/ESTAR/CSDA anchors, the two-route W
  agreement (0.02–0.19%), and the side finding that BOTH W routes give
  W(Ne) ≈ 43.7 eV at E0 — the +16% literature miss is data/method class.

### Could-not-check (named)

Petrova 2022 / Electra scoring (S7 remainder — the Orestes comparison
proper); Chow 1978 K15 page verification (kept at 2.7e-31, corroborated by
Kannari's 2.5e-31 class and the NK4 adjudication); the hot-vs-thermal
molecular-ion sigma choice (competing treatments both in-repo).

### Changed (code)

network.py (4 rate fixes, 6 cite upgrades EST/REPO→LIT, 4 absorbers +
routing), beam.py (W(F2) 33→32, F_EXC(Kr) 0.5→0.28, Kannari-cited),
run.py (preset geometry fix), scorecard.py (recorded-device-anchors
block), tests (census ratchet: EST ≤ 18, LIT ≥ 12), status.html (51
items; validation leaf + anchored-gap metrics).

### Verified

- Suite **129 passed** (~36 s) after the corrections; ledger closure
  unchanged (1e-11 class on the v2 run; 3e-5 on the stiff Chen runs).
- v2-point impact of the corrections, measured: g0 0.0385 → **0.0378 /cm**
  (spectre +19%, phantom11 +25%), eta 11.16 → **10.62%** (still under the
  hv/W bound 16.4%). Census: **EST 22 → 18, LIT 2 → 12**.

## [0.1.2] — 2026-08-04 — two gaps vs PIC_MCC_3D closed: 2nd W route + analytic trajectory gates

Owner ask: from the REBUILD-vs-PIC_MCC_3D gap list, "proceed with the first
two" (the second independent W route; analytic end-to-end trajectory gates).

### Added
- **The SECOND W route** — `tools/export_kinetic.py` now also runs a
  Monte-Carlo degradation cascade through PIC_MCC_3D's ACTUAL MCC operator
  (`beam3d.cascade_yields`) per pure gas and ships the result in the export
  (`meta["w_cascade"]`: E0, W, U_ion, rounds, seed-derived stat error).
  `kinetic_tier.KineticRates.w_cascade` exposes it;
  `scorecard.two_route_w_rows` compares REBUILD's deterministic
  Spencer–Fano against it AT THE SAME recorded E0 on the same LXCat sets —
  independent code paths, cross-package, no engine import (the coupling
  rule holds; Ne's unmapped file is resolved by REBUILD's loader and handed
  to PIC_MCC_3D as an explicit path rather than editing the sibling's
  GAS_FILES — the ARPAE-Loschmidt rule).
- **Analytic end-to-end trajectory gates**
  (`tests/test_analytic_trajectories.py`, 4 gates) — the
  plasma-oscillation-gate analog this package lacked: REAL network rows
  integrated through the REAL CompiledNetwork (BDF rtol 1e-11, compared at
  1e-5..1e-6 — three decades of slack, per the loose-tolerance Known-bad)
  against closed forms: URAD exponential decay + fragment growth; the K25
  unequal-concentration second-order solution INCLUDING its exact
  release-heating ΔT_g; K15 pseudo-first-order three-body; and the
  Device-level erf integral of the Gaussian drive against the ledger's
  drive_in channel.
- Two new kinetic-tier gates (cascade provenance; two-route agreement at
  5 % ≈ 3–4σ of the MC statistics so a reseeded export cannot flap it).

### Verified (measured 2026-08-04)
- Suite **123 → 129 passed** (~64 s, /usr/bin/python3); the regenerated
  export re-ran clean (same 7 E/N points + 3 cascades).
- **The two W routes agree to Ar 0.19 %, Kr 0.04 %, Ne 0.02 %** (MC stat
  ~1.1–1.4 %, 278–397 rounds, E0 = 531–550 eV) — the cross-package twin of
  PIC_MCC_3D's in-package 0.4 % pairing. Side finding: **both routes give
  W(Ne) ≈ 43.7 eV at this E0**, confirming the +16 % differential-W
  literature miss is data/method class, not an implementation bug.
- The three-body trajectory gate initially FAILED at rtol 2e-6 — because
  the trajectory really does deviate 2.9e-6 from the constant-k
  exponential: the 1.1 eV/event release heats the gas 3.5e-4 K and the
  Tg^-0.75 factor feeds back at the ppm level. The pin is 1e-5 with the
  measurement recorded in the docstring — the coupling resolving ppm
  physics is the finding.

## [0.1.1] — 2026-08-04 — status.html: the system map + data-flow board

Owner ask: "Just like in PIC_MCC_3D let us have an *.html to track the
status of the software."

### Added
- **`status.html` — the system map.** Shell (CSS + render/filter/flow
  machinery) adapted verbatim from `PIC_MCC_3D/status.html` — the point of
  reusing a verified artifact is that its behaviour is known; only the two
  data literals and the framing text are REBUILD's.
  - **System-map view**: inputs → the engine → outputs → gates; **48
    statused items (39 working / 5 open / 4 missing)**, 22 carrying the
    measured number that backs them, every figure naming the dated run it
    comes from (2026-08-04). Status vocabulary is GUI_v2 tab ⑥'s; tags are
    the plan's S1…S7 milestones and D1…D5 defect classes. The five [warn]
    rows are the honest open items (EST count, the vib-bucket finding, the
    ⟨ε⟩-mapping approximation, the W(Ne) miss, the prescribed discharge);
    the four [gap] rows are in-loop ASE, inter-cell transport, AMR, and
    Chen-figure scoring.
  - **Data-flow view**: one RHS EVALUATION as the eight-stage cycle
    (unpack → tier rates → fluxes → energy routing → drives → photon field
    → temperatures → ledger channels) with the BDF integrator as the
    feedback ring wrapping the diagram. Every step-loop edge carries what
    crosses it, the invariant that must hold there, the GATE that enforces
    it (real test names), and the symptom you would see if it broke.
  - A hand-maintained snapshot, not a live probe — the page says so; if it
    and the code disagree, the page is wrong.

### Verified
- The page's real script executed under a stubbed DOM via
  `osascript -l JavaScript` (the PIC_MCC_3D convention; no headless browser
  on this machine): parses and runs end-to-end including the SVG draw path;
  tally 48 items; flow graph 19 nodes / 19 edges, every edge resolving to a
  real node, no duplicate ids, all boxes inside the 1000×726 viewBox,
  **0 box overlaps, 0 dogleg channels passing through boxes**, every
  non-report edge fully annotated. NOT machine-checked (no browser):
  typography/label metrics — same caveat PIC_MCC_3D's board records.

## [0.1.0] — 2026-08-04 — the engine, S1–S7 v1, all gates green

Owner ask: "proceed with the rebuild under REBUILD/" — the
ORESTES_CLASS_REBUILD.txt program, housed in this new package (superseding
that file's §0 recommendation to harden spectre in place; the plan's
architecture, tiers, milestones and policies carry over unchanged).

### Added

- **The engine** (`MODULES/`, 17 modules, import path
  `REBUILD.MODULES.<name>`): tiered reaction network as data
  (CHEN verbatim + REPO/LIT/EST from the audited spectre network + LOKI
  from packaged curated LoKI-MC tables + KINETIC from a PIC_MCC_3D swarm
  export), RHS assembled mechanically, e-beam deposition from Bethe/Møller
  stopping power with EMERGENT W (Spencer–Fano per-process recursion),
  whole-engine energy ledger integrated as state, Chen Eqs. (7)–(18)
  multi-zone cavity with a ledgered loss split, voxel-field 3-D ASE +
  parasitic margins, coarse 3-D device grid under one BDF with
  block-diagonal Jacobian sparsity, declared parameters (flare pattern)
  that REFUSE calibration-knob names, checkpoint/resume, scorecard, CLI
  (`run.py`), verification runner (`verify.py`).
- **Modes**: `krf` (Ar-fill default = the DE-AR device; Ne-fill preset) —
  89 rows Ar-fill / 86 Ne-fill; `arf_chen` (Chen-verbatim scoring mode) —
  40 rows. Kannari p→s relaxation rows are IN unconditionally (KPS_*/APS_*;
  the PLAN (72) closure — no flag exists).
- **Network completeness delta over spectre, recorded**: e+Kr ground
  excitation/ionization rows (LK1–LK3) from the V2point mixture tables —
  spectre's KrF mode had no e+Kr rows at all (beam-only Kr manifolds).
- **`tests/` — 123 falsification gates** (~32 s) + `verify.py` (module
  isolation, per-file == combined == reversed order, 5-rung composition
  ladder). Teeth proven: an unbalanced row, a charge-unbalanced row, an
  over-radiating row, a banned knob name, an unrouted 2 eV/event leak, and
  a source-energy above the cross-section ceiling ALL fail their gates.
- **`tools/export_kinetic.py`** — the PIC_MCC_3D swarm exporter (imports
  PIC_MCC_3D OUTSIDE the engine; 7 E/N points, 800 macros; provenance in
  the npz meta + `MODULES/kinetic_tables/`).
- Data packaged with provenance: `MODULES/loki_tables/` (five curated
  LoKI-MC mixture sweeps, byte-verbatim from projects/spectre, README.txt
  records the copy) and `MODULES/kinetic_tables/pic_mcc_3d_swarm.npz`.

### Verified (all numbers measured on this checkout, 2026-08-04)

- `python3 -m pytest REBUILD/tests -q` → **123 passed in ~32 s**
  (`/usr/bin/python3` 3.9.6 / numpy 2.0.2 / scipy 1.13.1).
- `python3 REBUILD/verify.py` → **ALL CHECKS PASS** (17/17 modules import
  alone silently; 16 test files each green in their own process; combined
  123 = sum of isolated; reversed order green; ladder 5/5).
- Two-interpreter rule (PLAN (68)): all modules import + a 3-D run + ASE
  report complete under Homebrew **3.14.6 / numpy 2.5.0 / scipy 1.18.0**
  with `DeprecationWarning` promoted to error.
- **S3 headline — the energy budget closes as an invariant**: worst
  residual **8.2e-6** of energy moved across the full 200 ns pulse (0-D
  and 96-cell 3-D, beam and discharge drives; rtol 1e-6 class — the
  residual is integrator error by construction). The teeth test injects an
  unrouted 2 eV/event and the gate FAILS as required.
- **Emergent W** (differential window): Ar 27.7 eV (+4.9 % vs 26.4), Kr
  26.2 (+8.1 % vs 24.2) — matching PIC_MCC_3D's independent measurements
  to 3 digits; Ne 42.3 (+15.7 % vs 36.6 — a real miss, stated and pinned).
  Stopping power vs ESTAR: −1.9 % to −5.7 % (the PIC_MCC_3D band). CSDA
  600 keV / 1 atm Ar: 1.938 m (PIC_MCC_3D: 1.94). Conserving per-pair
  identity exact (1e-9 gate) on both yield routes.
- **S4 — KINETIC vs LOKI over the tabulated overlap** (<eps> 3.4–7.7 eV):
  elastic x0.87–0.90, attachment x0.95, F2 electronic x0.77–0.90,
  excitation medians x1.4–1.6 (the pic1d x1.6 class), ionization x1.9–2.8
  with threshold-tail sensitivity at the low-energy edge. **FINDING: the
  F2 'vib' bucket disagrees by ~x100** — reported on the scorecard, engine
  defaults to LOKI.
- Device physics at the v2-class point (Ar 94.9/Kr 5/F2 0.1, 3 atm,
  600 kV, 53 A/cm2, 0-D): ne peak 1.1e15 cm^-3, Te ~0.5 eV, **g0 peak
  0.0385 /cm vs spectre 0.0318 (+21 %) and phantom11 0.0302 (+28 %)** —
  engine-vs-engine at matched conditions, independent architectures;
  eta_intrinsic 11.2 % (< the mechanical hv/W bound 16.4 % — the PLAN (72)
  breach class is impossible here by physics, and a gate asserts it).
  3-D (96 cells) reproduces 0-D in the uniform limit (11.19 % vs 11.16 %).
- ASE: voxel chord integration matches the unified/ase3d closed form to
  1e-15 (axial + oblique); two-segment hand computations exact; the
  pure-gain order-independence identity S(ab−1) verified.

### Scope honesty (also on the scorecard as OMISSIONS)

Xe rows R16/R17 off (trace, stopping only); no Kr2*/Ar2* VUV decay rows
(spectre parity); fluid Te (no e-e operator — both rate tiers are
swarm-kinetic); prescribed discharge waveform (no circuit); ASE is
quasi-steady report + saturation fixed point, not an in-loop sink; no
inter-cell transport (cm cells, 100–200 ns); checkpoint resume is
tolerance-exact, not bitwise (adaptive BDF); MESH.txt adaptive refinement
is a roadmap item — v1 has the density-responsive deposition march and a
stated uniform grid.

### Unverified-by-literature (flagged, per the repo verification policy)

Species formation energies for dimers/excimers are EST-tier estimates
(they move only the heat/potential split, never the closure — see
species.py); the Kannari p→s rate value (1e-16 m^3/s class) is cited via
the repo's PLAN (72) record, not yet re-read from Kannari 1985 first-hand;
Chen-paper output observables are not yet scored against the paper's
figures (the arf_chen drive is a stand-in waveform).
