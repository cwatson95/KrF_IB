"""Put this tests directory's own LIGHT/ first on sys.path.

pytest finds a test's package by walking up the directories that hold an
__init__.py, and stops at a name that is not a Python identifier. '1Dv0_3'
is not one, so pytest never puts LIGHT/ on sys.path (it does for v0_2),
and importlib.import_module('1Dv0_3.model') would find no package, or the
one in the current directory instead of this one.
"""

import sys
from pathlib import Path

LIGHT = str(Path(__file__).resolve().parents[2])
if sys.path[0] != LIGHT:
    sys.path.insert(0, LIGHT)
