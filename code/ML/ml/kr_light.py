"""LIGHT_Kr, the pure-krypton model, for ML's krypton swarm analysis.

LIGHT_Kr's packages (v0_1, v0_2, v0_3) share their names with LIGHT's,
which ml puts on sys.path, so the LIGHT_Kr version ML trains is loaded
under a name of its own, light_kr_<version> (its relative imports resolve
inside it). LIGHT_KR_VERSION names the version (default v0_3, the first
whose cross sections take trained factors, xsec.Factors); like
LIGHT_VERSION it is read once, when this module is imported.
"""
from __future__ import annotations

import importlib
import importlib.util
import os
import sys

from . import ROOT

LIGHT_KR_DIR = ROOT / "LIGHT_Kr"
LIGHT_KR_VERSION = os.environ.get("LIGHT_KR_VERSION", "v0_3")
ALIAS = f"light_kr_{LIGHT_KR_VERSION}"


def package():
    """The LIGHT_Kr version's package, loaded once under ALIAS."""
    if ALIAS in sys.modules:
        return sys.modules[ALIAS]
    if not (LIGHT_KR_DIR / f"run_{LIGHT_KR_VERSION}.py").is_file():
        raise ImportError(f"LIGHT_KR_VERSION={LIGHT_KR_VERSION!r}: there is "
                          f"no LIGHT_Kr/run_{LIGHT_KR_VERSION}.py")
    pkg_dir = LIGHT_KR_DIR / LIGHT_KR_VERSION
    spec = importlib.util.spec_from_file_location(
        ALIAS, pkg_dir / "__init__.py",
        submodule_search_locations=[str(pkg_dir)])
    mod = importlib.util.module_from_spec(spec)
    sys.modules[ALIAS] = mod
    spec.loader.exec_module(mod)
    return mod


def module(name: str):
    """A module of the LIGHT_Kr version, module('swarm') say."""
    package()
    return importlib.import_module(f"{ALIAS}.{name}")
