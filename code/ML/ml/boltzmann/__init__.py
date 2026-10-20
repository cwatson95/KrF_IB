"""A two-term Boltzmann solver for the EEDF, vendored from the website's
LXCat simulation build (Agent/projects/simulations/build_lxcat, 2026-09-27):
lx_lib.py, gases_twin.py and boltz_twin.py, copied unchanged.

It is the Hagelaar & Pitchford (2005) scheme that LoKI-B solves, and it reads
LoKI-MC's LXCat files and setups the way LoKI-MC does. Against LoKI-MC's Monte
Carlo it agreed to < 0.5 % in mean energy, < 1 % in drift velocity and 1-2 %
in alpha/N for Kr, Ar, He, Xe and a Kr/F/Ar mixture (build_lxcat/boltz_out.txt,
section D4). Transverse diffusion is the known two-term weakness (Kr D/mu
40-70 % high). Pure F2 below about 30 Td has no positive two-term solution.
md5 at vendoring: lx_lib df49167189c54b93ffebe3b7c7a78217, gases_twin
7ed3a8edc10104caa2ab28ec95bfc41e, boltz_twin 109eeb7e5eedc8e92ce75923a3494473.

ML's check against LoKI-MC for He + 1 % F2 (2026-09-27): 1 % at 6-15 Td,
and a failed (negative mean energy) solution at 3 Td, where attachment
outruns energy relaxation. The same failure explains Kr/F2 99.5/0.5 at
0.3 Td.

ee.py is ML's own addition, not vendored: the electron-electron collision
term and the beam-sustained EEDF.
"""
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
if str(_here) not in sys.path:
    sys.path.insert(0, str(_here))          # the three modules import each other flat
