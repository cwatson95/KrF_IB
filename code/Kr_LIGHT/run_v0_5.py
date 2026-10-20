"""Launcher for LIGHT_Kr v0_5, the 0-D pure-krypton model: python3
LIGHT_Kr/run_v0_5.py --help (from anywhere).

Puts this file's own directory (LIGHT_Kr/) on sys.path so the v0_5 package
imports; no parent counting, no repo-root assumptions.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v0_5.run import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
