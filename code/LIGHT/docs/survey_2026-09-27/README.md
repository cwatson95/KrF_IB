# Survey evidence, 2026-09-27

The evidence behind `LIGHT/IMPROVEMENT_PLAN.md`: three survey notes on the
rest of Kinetic_Modeling and the scripts that produced every number the plan
quotes from LIGHT.

## Notes

| File | What it covers |
|---|---|
| `notes_A_reaction_sources.md` | Every reaction set in the tree and in the rate tables the repo holds, with counts. The candidate reactions REBUILD lacks, each with its printed value and page. Value problems in rows REBUILD already has |
| `notes_B_experimental_data.md` | Every experimental dataset in or cited by the tree, what LIGHT would need to reproduce each one, and a ranking for 0-D use |
| `notes_C_missing_components.md` | What beyond the reaction list limits accuracy, with the effects the repo has already measured. The recombination audit. The traps the repo records |

Three read-only survey agents compiled the notes. They cite a file and line,
or a printed page, for each claim. The claims the plan relies on were
re-checked first-hand (list below). Treat the rest as leads to verify before
use.

## Re-checked first-hand for the plan

- REBUILD 0.1.40's krf network, built and dumped: no F⁻ detachment, no
  three-body electron-ion recombination, no F + F + M, one attachment row.
- Xu, Gadomski & Setser 1993, p.2591: k_BC(Kr) = 1.2e-10 cm³/s (their
  ref. 17, Hakuta). The coupled lifetime is 25 ± 2 ns, so the B/C
  equilibrium constant "must be 4–5". The Kr quench constants
  (4.0e-12, 6.0e-31) were assigned under that structure.
- Xu 1993, p.2594: k' = 5.0e7 and 6.2e7 s⁻¹ at 50 and 100 Torr Kr with
  1 Torr F2, and the Table I quench values.
- Kannari 1985, Table I p.4311:
  - Kr* + e ↔ Kr** + e: 9.00(−07) Te^0.50 exp(−1.51/Te)
  - Kr** + e → Kr⁺ + 2e: 2.00(−07) Te^0.62 exp(−2.53/Te)
  - Kr* + e → Kr⁺ + 2e: 7.00(−08) Te^0.68 exp(−4.04/Te)
  - KrF* + F: 7.8(−10), "analogous to F2"
  - F + F + M: 1.0(−33) cm⁶/s
- `SPECK/REBUILD/MODULES/chen2025.py` and `network.py`:
  - Chen's R23 prints 2.5e-43 m⁶/s, and REBUILD's errata correct it to
    1.1e-44.
  - NK6 still carries 2.5e-43 as an "analog of Chen R23".
- `~/coding/PIC/Neon/NRL_REPLICATION.md`: the Mandl 99.6 %/10 % Kr ratio
  comes out 0.27 against 0.82 measured, and halving K34 gives 0.75.
- The phantom11 Stevefelt row
  (`KrF_and_LoKI/KrF Code/factored-pyfiles/phantom11/config.py:657-669`):
  its own formula gives 1.3e-42 cm⁶/s at 0.3 eV, where its comment says
  2.4e-22.

Not re-checked:
- the 500-Torr decay value (2.63e8 s⁻¹), which the survey read off a fitted
  curve in Gadomski 1992;
- the labels in `SPECK/EXPERIMENTS.md`;
- the page citation in `excimer_anchors.py`;
- phantom's r120 row;
- everything else in the notes.

## Scripts

Run them from anywhere with `/usr/bin/python3 <script>`. They find LIGHT by
walking up to `run_v0.py`. They import LIGHT v0 but are not part of it.
Two of them contain prototypes that exist only to measure importance, not to
predict:
- `ebeam_screen.py` adds an e-beam source using REBUILD's Spencer–Fano yields
  for Kr with 0.5 % F2 (27.07 eV per ion pair: Kr⁺ 0.999, Kr** 0.512,
  Kr* 0.263 per pair, and 4.5 eV of slow-electron energy per pair).
- `kstar_prototype.py` patches an explicit Kr** state into LIGHT at run time.

| Script | Measures | Time | Headline |
|---|---|---|---|
| `decay_run.py` | LIGHT against the Xu/Gadomski photolysis decays, pure Kr/F2, seed mode | 20 s | v0 ×2.89 / ×2.41 / ×1.13 at 50 / 100 / 500 Torr Kr. The B/C structure with K_eq 5 gives ×1.01 / ×1.01 / ×1.00 |
| `loo_screen.py` | Leave one row out of v0 + all groups, power drive, 5 points | 10 min | Largest changes in the 248-nm share (points of E_dep): CB_Kr 6.8, BC_Kr 6.1, K67 5.0, K33 4.4, G3 4.1, G5 3.0, KR2_RAD 2.1 |
| `ebeam_screen.py` | Prototype e-beam drive: KrF(B) formation by channel, leave-one-out at 1 atm | 2 min | Ion-ion recombination forms 48 % of KrF(B) at 1 atm and 53 % at 2 atm |
| `recomb_isolate.py` | Recombination (the `ions` + `f2e` groups) on and off, both drives, 0.5 to 2 atm | 3 min | E-beam: 248 nm ×9 to ×14 with it on. Power drive: under 1 point |
| `group_drop.py` | Each optional group removed from the full set, both drives | 5 min | Reports one solver failure at 2 atm under the e-beam drive, the free-n_e floor problem |
| `keq_stakes.py` | B/C equilibrium constant 1.47 against 5, all groups | 1 min | 248 nm 5.94 → 2.93 % (power, 1 atm) and 6.65 → 3.29 % (e-beam) |
| `kstar_prototype.py` | Explicit Kr** with electron up-pumping (Kannari [76], [79]) | 1 min | 248 nm 5.94 → 3.75 % at n_e 1e15 (power drive), 6.65 → 6.54 % (e-beam) |
