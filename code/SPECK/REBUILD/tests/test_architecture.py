"""S2/D5 — the package's SHAPE: layered imports, no mutable module state.

Ported discipline from PIC_MCC_3D/tests/test_architecture.py (the PLAN (74)
rule mechanised): the import graph is acyclic and layered, every module
imports alone, importing the package prints nothing, and module-level state
does not move when the engine runs.
"""

import ast
import subprocess
import sys
from pathlib import Path

import pytest

MODULES = Path(__file__).resolve().parents[1] / "MODULES"
PKG_ROOT = MODULES.parents[1]          # Kinetic_Modeling root

LAYERS = {
    "units": 0,
    "provenance": 1, "params": 1,
    "species": 2, "chen2025": 2, "cross_sections": 2,
    "network": 3, "loki_tier": 3, "kinetic_tier": 3, "maxwell_tier": 3,
    "beam": 3,
    "ledger": 4, "laser": 4, "ase": 4, "chemistry": 4, "circuit": 4,
    "engine": 5, "checkpoint": 5, "scorecard": 5,
}


def _module_names():
    return sorted(p.stem for p in MODULES.glob("*.py")
                  if p.stem != "__init__")


def _internal_imports(name):
    """All intra-package imports of a module (top-level AND nested)."""
    tree = ast.parse((MODULES / f"{name}.py").read_text())
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.level == 1 and node.module:          # from .x import y
                out.add(node.module.split(".")[0])
            elif node.module and node.module.startswith("REBUILD.MODULES."):
                out.add(node.module.split(".")[2])
    return out


def test_every_module_is_layered():
    missing = [m for m in _module_names() if m not in LAYERS]
    assert not missing, f"modules without a declared layer: {missing}"


def test_import_graph_is_acyclic_and_layered():
    graph = {m: _internal_imports(m) for m in _module_names()}
    # layering: an import may only reach the SAME layer or lower — and the
    # same-layer allowance must stay acyclic (checked next)
    for m, deps in graph.items():
        for d in deps:
            assert LAYERS[d] <= LAYERS[m], (
                f"{m} (layer {LAYERS[m]}) imports {d} (layer {LAYERS[d]}) — "
                "the architecture is inverted")
    # acyclicity: Kahn's algorithm
    indeg = {m: 0 for m in graph}
    for m, deps in graph.items():
        for d in deps:
            if d in indeg:
                indeg[m] += 0          # deps counted below
    indeg = {m: sum(1 for d in graph[m] if d in graph) for m in graph}
    order = [m for m, k in indeg.items() if k == 0]
    seen = 0
    from collections import deque
    q = deque(order)
    rdeps = {m: [n for n, ds in graph.items() if m in ds] for m in graph}
    while q:
        m = q.popleft()
        seen += 1
        for n in rdeps[m]:
            indeg[n] -= 1
            if indeg[n] == 0:
                q.append(n)
    assert seen == len(graph), "the import graph has a cycle"


def test_no_sibling_package_imports():
    """The coupling rule: PIC_MCC_3D/phantom11/spectre are reachable only
    as DATA. tools/ may import them; MODULES/ may not."""
    for m in _module_names():
        tree = ast.parse((MODULES / f"{m}.py").read_text())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                mods = [node.module]
            for mod in mods:
                root = mod.split(".")[0]
                assert root not in ("PIC_MCC_3D", "phantom11", "spectre",
                                    "unified", "pic1d", "pic3d"), (m, mod)


@pytest.mark.parametrize("name", sorted(LAYERS))
def test_module_imports_alone_silently(name):
    r = subprocess.run(
        [sys.executable, "-c", f"import REBUILD.MODULES.{name}"],
        cwd=str(PKG_ROOT), capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    assert r.stdout == "" and r.stderr == "", (name, r.stdout, r.stderr)


def _state_snapshot():
    """repr of every module-level data attribute of the engine modules."""
    import importlib
    snap = {}
    for m in _module_names():
        mod = importlib.import_module(f"REBUILD.MODULES.{m}")
        for k, v in vars(mod).items():
            if k.startswith("_") or callable(v) or isinstance(
                    v, type(sys)) or isinstance(v, type):
                continue
            snap[f"{m}.{k}"] = repr(v)
    return snap


def test_module_level_state_does_not_move():
    """PLAN (74): 10 of 17 suite failures traced to leaked globals. Build
    and run a device, then diff every module-level value."""
    before = _state_snapshot()
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "t_end_ns": 20.0,
                "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141,
                "Ly_m": 0.0141, "Lz_m": 0.75, "beam_axis": 0},
               source="problem")
    Device(p).run(n_save=5)
    after = _state_snapshot()
    moved = {k for k in before if before[k] != after.get(k)}
    assert not moved, f"module-level state moved: {sorted(moved)}"


def test_repeat_run_is_deterministic():
    import numpy as np
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters

    def one():
        p = engine_parameters()
        p.set_many({"mode": "krf", "drive": "beam", "t_end_ns": 20.0,
                    "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141,
                    "Ly_m": 0.0141, "Lz_m": 0.75, "beam_axis": 0},
                   source="problem")
        return Device(p).run(n_save=5)["y_final"]

    a, b = one(), one()
    assert np.array_equal(a, b), "same configuration, different answer"


def test_package_import_prints_nothing():
    r = subprocess.run(
        [sys.executable, "-c", "import REBUILD; import REBUILD.MODULES"],
        cwd=str(PKG_ROOT), capture_output=True, text=True, timeout=60)
    assert r.returncode == 0 and r.stdout == "" and r.stderr == ""
