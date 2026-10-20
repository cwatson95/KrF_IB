"""LIGHT 1Dv0_2: the 1-D model on v0_2's chemistry, across a gap on Fluid1D.

1Dv0_2 starts as a copy of v0_2 (2026-09-29) with 1Dv0's 1-D layer on it.
Its chemistry (model.py, network.py, species.py, xsec.py, two_level.py,
report.py and data/) is v0_2's, unchanged: the Kr2* route in the basic
network and the rates known only as a band, at their defaults unless
Params1D.rates sets them. model.py still runs it in 0-D. The 1-D model is
sim1d.py, on Nathaniel's Fluid1D (fluid1d.py is the bridge). The only
energy input is v0_2's: a given power into the electron gas. See
LIGHT/1Dv0_2/README.md.

'1Dv0_2' is not a Python identifier, so import it by name:
``importlib.import_module("1Dv0_2.sim1d").run()``. No import-time side
effects: Fluid1D is loaded when a run starts.
"""

__version__ = "0.1.2+1d0"
