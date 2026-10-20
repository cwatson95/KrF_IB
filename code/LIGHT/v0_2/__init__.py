"""LIGHT v0_2: the lightweight 0-D Kr/F2 model, v0_1's working copy.

v0_2 starts as a copy of v0_1 (2026-09-29) and carries the work on the Kr2*
route and on the rates that are known only as a band. v0 and v0_1 stay as
they were. The only energy input is a given power, deposited uniformly and
directly into the electron gas. See LIGHT/v0_2/README.md.

Canonical imports: ``from v0_2.model import Params, run`` and
``from v0_2.report import summary``. No import-time side effects.
"""

__version__ = "0.1.2"
