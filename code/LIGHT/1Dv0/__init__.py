"""LIGHT 1Dv0: the 1-D model, v0_1's chemistry across a gap on Fluid1D.

1Dv0 starts as the owner's copy of v0_1 (2026-09-29). Its chemistry
(model.py, network.py, species.py, xsec.py, two_level.py, report.py and
data/) is v0_1's, unchanged, and model.py still runs it in 0-D. The 1-D
model is sim1d.py, on Nathaniel's Fluid1D (fluid1d.py is the bridge). The
only energy input is v0_1's: a given power into the electron gas. See
LIGHT/1Dv0/README.md.

'1Dv0' is not a Python identifier, so import it by name:
``importlib.import_module("1Dv0.sim1d").run()``. No import-time side
effects: Fluid1D is loaded when a run starts.
"""

__version__ = "0.1.1+1d0"
