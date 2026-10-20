"""Launcher for LIGHT 1Dv0_2, the 1-D model on v0_2's chemistry:
python LIGHT/run_1Dv0_2.py --help.

Needs Python >= 3.11 with numba, as Fluid1D does (~/.venvs/f2sweep/bin/python
on this Mac). Puts this file's own directory (LIGHT/) on sys.path; '1Dv0_2'
is not a Python identifier, so the package is imported by name.
"""

import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

main = importlib.import_module("1Dv0_2.run").main

if __name__ == "__main__":
    sys.exit(main())
