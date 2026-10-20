"""LIGHT_ArF_ebeam v0_2: the lightweight 0-D Ar/F2 (ArF) model.

LIGHT_ArF (2026-09-30): LIGHT's v0_2 converted from Kr/F2 to Ar/F2. The
equations and the ledger are KrF v0_2's; the species, the IST-Lisbon Ar
cross sections, the Ar tail table and the ArF network are new, every rate
cited to a measurement (docs/arf_data/SOURCES_CHECKED.md).
LIGHT_ArF_ebeam (2026-09-30) adds a second way to deposit the given power:
Params.drive = 'beam', an electron beam (beam.py, REBUILD's beam physics on
these cross sections), beside LIGHT's own input into the electron gas. See
LIGHT_ArF_ebeam/v0_2/README.md.

Canonical imports: ``from v0_2.model import Params, run`` and
``from v0_2.report import summary``. No import-time side effects.
"""

__version__ = "0.2.0-arf-ebeam"
