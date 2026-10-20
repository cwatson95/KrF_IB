"""Make the package importable when pytest is run from anywhere.

Mirrors ``ARPAE/pic3d/tests/conftest.py``: put the directory CONTAINING the
package on sys.path so ``import PIC_MCC_3D`` resolves whether the suite is
started from the repo root or from inside the package.
"""

import sys
from pathlib import Path

_PARENT = Path(__file__).resolve().parents[2]
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))
