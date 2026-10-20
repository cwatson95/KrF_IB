"""LIGHT_Kr v0_5: the lightweight 0-D model of pure krypton (2026-10-16).

v0_5 is v0_4 (the 5s levels resolved, on cross sections trained on
krypton's swarm data) with the limitations of v0_4's accuracy overview
taken up, as the owner chose on 2026-10-16:

* where nothing is measured, a published calculation sets the default
  (bands from 0 to twice it): Vriens & Smeets's ionization of an excited
  level (1980) for the stepwise ionizations out of the 5s levels and the
  5p rung, and for the ionization of Kr2*, which Flannery et al. 1981 put
  at ~1e-15 cm2 (excited.py); the electrons' quenching of Kr2* and the
  mixing of its two states, with no published value, as bands at 0;
* the 5s -> 5p pumps pair by pair, each charged its own threshold (v0_4
  charged every pair the gap to 2p10), their reverses out of the 5p rung's
  levels populated at the gas temperature, its radiation line by line on
  NIST's A-values, its relaxation by krypton with the reverse;
* the 5s levels' superelastic return to the ground state (detailed balance
  on their trained cross sections; Blagoev et al. 1984's afterglow rates a
  score) and their mixing by electrons (Blagoev & Pilosov 1981);
* data: Kholin et al.'s 3P0 rate has argon for its third body (the cite
  corrected), Kr2+ bound by Signorell et al.'s 9267.8 cm-1;
* the two-temperature closure's field tail table with an excited-
  population axis, x_s = [Kr(5s)]/[Kr] (tail_field.py).

See LIGHT_Kr/v0_5/README.md and LIGHT_Kr/CHANGELOG.md.

Canonical imports: ``from v0_5.model import Params, run`` and
``from v0_5.report import summary``. No import-time side effects.
"""

__version__ = "0.5.0-kr"
