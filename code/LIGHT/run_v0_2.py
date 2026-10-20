"""Launcher for LIGHT v0_2: python3 LIGHT/run_v0_2.py --help (from anywhere).

Puts this file's own directory (LIGHT/) on sys.path so the v0_2 package
imports; no parent counting, no repo-root assumptions.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v0_2.run import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
