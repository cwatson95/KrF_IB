"""LIGHT_Kr v0_2: the lightweight 0-D model of pure krypton (2026-10-14).

v0_2 is v0_1 with the electron's collisions on ground-state krypton booked
level by level: every one of the 51 excitation levels of Biagi's Magboltz
8.97 cross sections (LXCat's 'Biagi' database) is its own row at its own
threshold, where v0_1 read Biagi-v7.1's four groups at the groups' lowest
thresholds, and the levels at or above Kr2+'s energy carry a
Hornbeck-Molnar branch, Kr(nl) + Kr -> Kr2+ + e, on the band hm_kr (0 by
default: no measurement for krypton). Everything else is v0_1's: LIGHT
v0_5's krypton network and machinery without its fluorine (the energy
ledger, the kinetic EEDF in the loop, the two-temperature closure, the
walls, the EM wave, the bands) and v0_1's two corrections to jung2006's
tables. See LIGHT_Kr/v0_2/README.md.

Canonical imports: ``from v0_2.model import Params, run`` and
``from v0_2.report import summary``. No import-time side effects.
"""

__version__ = "0.2.0-kr"
