"""Launcher for LIGHT v0: python3 LIGHT/run_v0.py --help (from anywhere).

Puts this file's own directory (LIGHT/) on sys.path so the v0 package
imports; no parent counting, no repo-root assumptions.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v0.run import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
