"""LIGHT_ArF_ebeam 1Dv0_2: the 1-D ArF model, LIGHT_ArF_ebeam v0_2's
chemistry across a gap on Fluid1D.

LIGHT_ArF (2026-09-30): LIGHT's 1Dv0_2 converted to Ar/F2. Its chemistry
(model.py, network.py, species.py, xsec.py, report.py and data/) is
LIGHT_ArF v0_2's, byte for byte: the ArF basic network with its Ar2* route
and the rates known only as a band, at their defaults unless Params1D.rates
sets them. model.py still runs it in 0-D. The 1-D model is sim1d.py, on
Nathaniel's Fluid1D (fluid1d.py is the bridge). The energy input is
v0_2's, either way: a given power into the electron gas (drive 'power')
or an electron beam (drive 'beam', beam.py; LIGHT_ArF_ebeam, 2026-09-30).
See LIGHT_ArF_ebeam/1Dv0_2/README.md.

'1Dv0_2' is not a Python identifier, so import it by name:
``importlib.import_module("1Dv0_2.sim1d").run()``. No import-time side
effects: Fluid1D is loaded when a run starts.
"""

__version__ = "0.2.0-arf-ebeam+1d0"
