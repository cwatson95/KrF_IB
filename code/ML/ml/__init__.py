"""ML: fitting LIGHT's reaction channels to the measurements.

The package puts LIGHT/ (beside ML/) on sys.path, so the forward models run
LIGHT's own code. That is v0 unless the environment variable LIGHT_VERSION
names another version (v0_1 since 2026-09-29). The version is fixed when ml
is first imported, so scoring another version needs its own process. See
ML/README.md.
"""
import importlib
import os
import sys
from pathlib import Path

_here = Path(__file__).resolve()
# marker walk to Kinetic_Modeling/ (the directory holding LIGHT/run_v0.py)
ROOT = next(p for p in _here.parents if (p / "LIGHT" / "run_v0.py").is_file())
LIGHT_DIR = ROOT / "LIGHT"
if str(LIGHT_DIR) not in sys.path:
    sys.path.insert(0, str(LIGHT_DIR))

# the LIGHT version ML runs: a package in LIGHT/ with its launcher run_<v>.py
LIGHT_VERSION = os.environ.get("LIGHT_VERSION", "v0")
if not (LIGHT_DIR / f"run_{LIGHT_VERSION}.py").is_file():
    raise ImportError(f"LIGHT_VERSION={LIGHT_VERSION!r}: there is no "
                      f"LIGHT/run_{LIGHT_VERSION}.py")


def light(module: str):
    """A module of the LIGHT version ML runs, light('two_level') say."""
    return importlib.import_module(f"{LIGHT_VERSION}.{module}")
