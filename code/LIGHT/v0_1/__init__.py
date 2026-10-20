"""LIGHT v0_1: the lightweight 0-D Kr/F2 model, v0's working copy.

v0_1 starts as a copy of v0 (2026-09-28) and carries the work on the
trimer, ion recombination, fluorine recombination and the KrF C-A band and
decay. v0 stays as it was. The only energy input is a given power,
deposited uniformly and directly into the electron gas. See
LIGHT/v0_1/README.md.

Canonical imports: ``from v0_1.model import Params, run`` and
``from v0_1.report import summary``. No import-time side effects.
"""

__version__ = "0.1.1"
