"""LIGHT_Kr v0_4: the lightweight 0-D model of pure krypton (2026-10-15).

v0_4 is v0_3 (the electron cross sections trained on krypton's swarm data
by ML) with section 7 of v0_3's accuracy overview taken up, as the owner
chose on 2026-10-15:

* the four 5s levels resolved, each its own species at Biagi's energy:
  each forms Kr2* at its own measured rate (Marchal et al. 2018: 3P2
  3.81e-32, 3P1 9.24e-32 cm6/s), is pumped to the 5p rung on its own
  cross sections (Jung et al. 2006, the J = 1 levels on their
  metastable's), and is fed by the 5p lines on NIST's branching; the
  levels mix on krypton at their measured rates (Turner 1967, Audouard
  et al. 1991);
* the redissociation of Kr2* and Kr2+ by detailed balance on their
  measured binding energies (thermo.py);
* bands for what is bounded but not measured at the fill: the fast 0u+
  excimer, Kr2* quenching by Kr, the gas-temperature factor of Kr2+
  recombination, the escaping resonance lines, the 5p relaxation's
  products;
* a field-driven pure-krypton tail table for the two-temperature closure,
  r(n_e/N, Te) from this model's own kinetic EEDF (tail_field.py), the
  afterglow's EEDF time-stepped on request (the bound on its steady
  shape), and the swarm's electroluminescence yield for a new score
  (swarm.el_yield, Mano et al. 2022).

See LIGHT_Kr/v0_4/README.md and LIGHT_Kr/CHANGELOG.md.

Canonical imports: ``from v0_4.model import Params, run`` and
``from v0_4.report import summary``. No import-time side effects.
"""

__version__ = "0.4.0-kr"
