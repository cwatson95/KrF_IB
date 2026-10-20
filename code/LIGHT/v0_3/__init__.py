"""LIGHT v0_3: v0_2 with the kinetic EEDF in the loop.

v0_3 starts as a copy of v0_2 (2026-09-30). Its chemistry, bands and input
are v0_2's. What changes is the electrons. v0_2 integrates every electron
rate over a Maxwellian at Te, with the Kr rows' tail read from a table made
from beam-driven kinetic EEDFs. v0_3 solves the two-term Boltzmann equation
with electron-electron collisions during the run (ML's discretization,
boltzmann.py), with the network's own electron sources and sinks, at
LIGHT's own mean energy (held by a field, the input power in the pulse and
the electrons' own energy after it), and takes every electron rate, the
elastic loss and the captured energies from that EEDF. v0, v0_1 and v0_2
stay as they were. See LIGHT/v0_3/README.md.

Canonical imports: ``from v0_3.model import Params, run`` and
``from v0_3.report import summary``. No import-time side effects.
"""

__version__ = "0.1.3"
