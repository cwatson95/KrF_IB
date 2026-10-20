"""LIGHT v0_5: v0_4 with the 5s -> 5p electron excitation on Jung et al.'s
measured cross sections (2026-10-09).

v0_5 starts as a copy of v0_4. K76/K76r, the electron pumping of the 5p
rung from Kr(5s) and its return, run on the electron-impact excitation
cross sections measured out of the two metastable levels of Kr into the
ten 4p5 5p levels (Jung, Stone, Boffard, Anderson & Lin, PRA 73, 022722
(2006); PRL 94, 163202 (2005)), lumped over the 5s manifold
(jung2006.py), under both closures; Kannari's fit, v0_4's default, stays
as Params.k76 = 'kannari'. The band k_5s5p is a scale on the measurement.
docs/krf_data/kleikamp_v0_5.py scores a Kr-buffered e-beam laser's
pressure optima (Kleikamp, Peters & Witteman 1986), held out. v0 to v0_4
stay as they were. See LIGHT/v0_5/README.md.

Canonical imports: ``from v0_5.model import Params, run`` and
``from v0_5.report import summary``. No import-time side effects.
"""

__version__ = "0.1.5"
