"""Make both halves importable when pytest runs this suite from anywhere.

Mirrors PIC_MCC_3D/tests/conftest.py: put SPECK (the directory containing
REBUILD and PIC_MCC_3D) on sys.path.
"""

import sys
from pathlib import Path

_SPECK = Path(__file__).resolve().parents[2]
if str(_SPECK) not in sys.path:
    sys.path.insert(0, str(_SPECK))
