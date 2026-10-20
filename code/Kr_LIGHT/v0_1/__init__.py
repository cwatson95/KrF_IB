"""LIGHT_Kr v0_1: the lightweight 0-D model of pure krypton (2026-10-13).

v0_1 starts as a copy of LIGHT's v0_5, the Kr/F2 KrF model, with every
fluorine species, row, band and parameter taken out and the 248-nm optics
with them. What is left is v0_5's krypton unchanged, row for row: the 5s
manifold Kr*, the 5p rung Kr** (relaxed by Kr, radiating the near-IR
lines, pumped from Kr* on Jung et al.'s measured cross sections), Kr+ and
Kr2+ (three-body conversion, dissociative recombination onto the 5p rung),
the excimer Kr2* and its 147-nm second continuum, Penning ionization and
the unmeasured rates as bands. The machinery is v0_5's: the energy ledger,
the kinetic EEDF in the loop (boltzmann.py) and LIGHT v0_2's two-
temperature closure, the walls (gap_cm), the EM wave with Turnbull's
inverse bremsstrahlung (ib.py) and the sensitivity hooks (rates,
row_scale, ladder, k76). With no F2 in the fill v0_5's fluorine rows carry
nothing, so on the same cross-section tables v0_1 gives v0_5's own pure-Kr
answer (tests/test_v0_1.py). Two entries of jung2006's tables are
corrected to the paper (SOURCES.md section 6). See LIGHT_Kr/v0_1/README.md.

Canonical imports: ``from v0_1.model import Params, run`` and
``from v0_1.report import summary``. No import-time side effects.
"""

__version__ = "0.1.0-kr"
