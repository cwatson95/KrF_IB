"""LIGHT 1Dv0_3: the 1-D model on v0_2's chemistry, on a mesh graded to the
walls by OSIRIS's load rule.

1Dv0_3 starts as the owner's copy of 1Dv0_2 (2026-10-01). Its chemistry
(model.py, network.py, species.py, xsec.py, two_level.py, report.py and
data/) is v0_2's, unchanged, and model.py still runs it in 0-D. What it
adds is the mesh: mesh.py places the cells as OSIRIS divides its grid
among processors, by equal shares of a load, here dense at the
electrodes, and transport.py runs Fluid1D's transport step on each cell's
own width. Since 2026-10-01 (c) the neutral gas moves too (neutrals.py):
it flows, conducts heat into the electrodes, and its F2 and F diffuse.
Since 2026-10-06 the field's work on the electrons and the ions is kept to
the field's energy (field_work="conserve", transport.py). On a uniform
mesh with the gas still and Fluid1D's field work (mesh="uniform",
neutrals="still", field_work="fluid1d") it runs 1Dv0_2's model to the
bit. The 1-D model is sim1d.py, on Nathaniel's
Fluid1D (fluid1d.py is the bridge). The only energy input is v0_2's: a given power into the
electron gas. See LIGHT/1Dv0_3/README.md.

'1Dv0_3' is not a Python identifier, so import it by name:
``importlib.import_module("1Dv0_3.sim1d").run()``. No import-time side
effects: Fluid1D is loaded when a run starts.
"""

__version__ = "0.1.2+1d3"
