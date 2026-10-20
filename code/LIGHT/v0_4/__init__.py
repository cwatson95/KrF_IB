"""LIGHT v0_4: v0_3 with the 5p rung, the 248-nm optics, KrF* + F, ion-ion
recombination's limbs and bands, the ions' drift to the walls and the EM
wave with Turnbull et al.'s inverse bremsstrahlung (2026-10-08).

v0_4 starts as a copy of v0_3. Its electrons are v0_3's kinetic EEDF, now
with the 5s <-> 5p collisions, the ions' Coulomb collisions and the wave's
high-frequency field term in the Boltzmann equation. Its chemistry adds
the 5p rung Kr** (C1 feeds it, K76/K76r pump and return it, KPS relaxes
it, PRAD radiates it, K27 harpoons from it, KSTEP_S/P ionize), KrF(B,C) +
F as a band, the ion-ion yield and B/C bands on a Flannery turnover whose
Langevin limb is the measured mobilities, and the wall loss of a 0-D gap.
optics.py gives g0 and alpha at 248 nm, ib.py the wave's absorbed power.
v0 to v0_3 stay as they were. See LIGHT/v0_4/README.md.

Canonical imports: ``from v0_4.model import Params, run`` and
``from v0_4.report import summary``. No import-time side effects.
"""

__version__ = "0.1.4"
