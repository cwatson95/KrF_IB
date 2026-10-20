"""Gates on the SHAPE of the package, not on any one module's physics.

"Does each module work on its own, and do they work together?" splits into
questions the physics gates cannot answer:

  * is the import graph acyclic and layered, or has something started
    reaching upward?
  * does every module import ALONE, in a fresh interpreter, with no
    import-time side effects?
  * does anything mutate module-level state, so that one test silently
    changes the next one's answer?

That last one is not hypothetical in this repo. PLAN (74) traced **10 of 17**
suite failures to exactly this: `krf0d` mutated `phantom11` globals and did
not restore them, and one leaked value (`CHOW_EN_TD_DEFAULT`, 51 -> 90 Td)
was enough on its own to break a parity gate. The rule that came out of it --
diff every module-level scalar before and after -- is mechanised here.

The per-file isolation sweep lives in ``PIC_MCC_3D/verify.py`` instead of
here, because running pytest inside pytest is a poor way to spend a suite.
"""

import ast
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

MODULES_DIR = Path(__file__).resolve().parents[1] / "MODULES"
PKG_PARENT = Path(__file__).resolve().parents[2]

MODULE_NAMES = sorted(p.stem for p in MODULES_DIR.glob("*.py")
                      if p.stem != "__init__")

# The intended layering. A module may import from STRICTLY LOWER layers only.
# Written down so a new import that quietly inverts the architecture fails
# here rather than being discovered later as a circular-import crash.
LAYERS = {
    0: {"units", "params"},
    1: {"gas", "poisson3d", "species3d", "theory"},
    2: {"cross_sections", "es_field", "emission3d"},
    3: {"beam3d", "es_sim", "mcc3d"},
    4: {"checkpoint"},
}
LAYER_OF = {m: k for k, ms in LAYERS.items() for m in ms}

# Third-party imports the package is allowed to have. numpy everywhere; scipy
# only in the Poisson oracle, which is gate-only code.
ALLOWED_EXTERNAL = {"__future__", "numpy", "dataclasses", "pathlib", "json",
                    "hashlib", "os", "sys", "typing", "collections", "math"}
SCIPY_ALLOWED_IN = {"poisson3d"}


def _imports(stem):
    """(intra-package deps, third-party roots) for one module, via AST."""
    tree = ast.parse((MODULES_DIR / f"{stem}.py").read_text())
    deps, ext = set(), set()
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            if n.level == 1:
                if n.module:
                    deps.add(n.module)
                else:
                    deps.update(a.name for a in n.names)
            elif n.module:
                ext.add(n.module.split(".")[0])
        elif isinstance(n, ast.Import):
            for a in n.names:
                ext.add(a.name.split(".")[0])
    return deps, ext


def test_every_module_is_placed_in_the_layering():
    assert set(MODULE_NAMES) == set(LAYER_OF), (
        f"modules not in LAYERS: {set(MODULE_NAMES) - set(LAYER_OF)}; "
        f"LAYERS names that do not exist: {set(LAYER_OF) - set(MODULE_NAMES)}")


@pytest.mark.parametrize("stem", MODULE_NAMES)
def test_imports_only_from_lower_layers(stem):
    """No upward or sideways imports -- the guard against a dependency cycle."""
    deps, _ = _imports(stem)
    for d in deps:
        assert d in LAYER_OF, f"{stem} imports unknown module {d!r}"
        assert LAYER_OF[d] < LAYER_OF[stem], (
            f"{stem} (layer {LAYER_OF[stem]}) imports {d} "
            f"(layer {LAYER_OF[d]}) -- that is not downward")


def test_import_graph_is_acyclic():
    """Independent of the declared layering: a real topological sort."""
    rem = {m: set(_imports(m)[0]) for m in MODULE_NAMES}
    order = []
    while rem:
        ready = sorted(m for m, d in rem.items() if not (d & set(rem)))
        assert ready, f"import cycle among {sorted(rem)}"
        order += ready
        for m in ready:
            rem.pop(m)
    assert len(order) == len(MODULE_NAMES)


@pytest.mark.parametrize("stem", MODULE_NAMES)
def test_third_party_imports_are_declared(stem):
    _, ext = _imports(stem)
    for e in ext:
        if e == "scipy":
            assert stem in SCIPY_ALLOWED_IN, (
                f"{stem} imports scipy; only {SCIPY_ALLOWED_IN} may "
                "(it is gate-only code, and the package must run without it)")
            continue
        assert e in ALLOWED_EXTERNAL, (
            f"{stem} imports {e!r}, which is not in ALLOWED_EXTERNAL -- add "
            "it deliberately or drop the dependency")


def test_each_module_imports_alone_in_a_fresh_interpreter():
    """Every module must stand up by itself, with no import-time side effects.

    Run in ONE subprocess but each import in its own, so a module that only
    works because something else was imported first is caught.
    """
    prog = (
        "import subprocess, sys\n"
        f"mods = {MODULE_NAMES!r}\n"
        "bad = []\n"
        "for m in mods:\n"
        "    r = subprocess.run([sys.executable, '-c',\n"
        "        'import PIC_MCC_3D.MODULES.%s as M; assert M is not None' % m],\n"
        "        capture_output=True, text=True)\n"
        "    if r.returncode: bad.append((m, r.stderr.strip().splitlines()[-1:]))\n"
        "print(bad)\n"
    )
    r = subprocess.run([sys.executable, "-c", prog], capture_output=True,
                       text=True, cwd=str(PKG_PARENT))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "[]", f"modules that fail to import alone: {r.stdout}"


def test_importing_the_package_has_no_side_effects():
    """Importing must not read data files, open sockets, or print."""
    r = subprocess.run(
        [sys.executable, "-c",
         "import PIC_MCC_3D, PIC_MCC_3D.MODULES; print('CLEAN')"],
        capture_output=True, text=True, cwd=str(PKG_PARENT))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "CLEAN", f"import printed: {r.stdout!r}"
    assert r.stderr.strip() == "", f"import wrote to stderr: {r.stderr!r}"


def _module_state():
    """Hashable snapshot of every module-level value that could be mutated."""
    import importlib
    snap = {}
    for m in MODULE_NAMES:
        mod = importlib.import_module(f"PIC_MCC_3D.MODULES.{m}")
        for k, v in vars(mod).items():
            if k.startswith("__"):
                continue
            if isinstance(v, (int, float, str, bool, type(None))):
                snap[f"{m}.{k}"] = repr(v)
            elif isinstance(v, (list, tuple, dict, set)):
                try:
                    snap[f"{m}.{k}"] = repr(sorted(v) if isinstance(v, set)
                                            else v)
                except TypeError:
                    snap[f"{m}.{k}"] = repr(v)
            elif isinstance(v, np.ndarray):
                snap[f"{m}.{k}"] = v.tobytes()
    return snap


def test_running_physics_does_not_mutate_module_globals():
    """The PLAN (74) rule, mechanised: diff every module-level value.

    A leaked global is the worst class of bug in a test suite -- every test
    passes alone and the combination fails, or worse, passes with a wrong
    number. This exercises a representative slice of the package and asserts
    nothing at module scope moved.
    """
    from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
    from PIC_MCC_3D.MODULES.gas import pure
    from PIC_MCC_3D.MODULES.poisson3d import PERIODIC
    from PIC_MCC_3D.MODULES.species3d import (ParticleSpecies, Species,
                                              SpeciesSet)
    from PIC_MCC_3D.MODULES.theory import cold_mode_frequency_ratio

    before = _module_state()

    L, nx = 1e-2, 10
    pop = ParticleSpecies.cold_wave(Species.electron(), 12 ** 3, [0, 0, 0],
                                    [L] * 3, 1e14, amplitude=4e-5, mode=1,
                                    axis=2)
    sim = ESSim(nx, nx, nx, L / nx, SpeciesSet([pop]), 1e-11, bc=PERIODIC,
                background_charge_density=uniform_background_for(pop, L ** 3),
                seed=1)
    sim.run(4)
    cold_mode_frequency_ratio(0.3)
    g = pure("Ar", 1.0)
    g.mean_free_path(1e-20)

    after = _module_state()
    changed = [k for k in before if before[k] != after.get(k)]
    added = sorted(set(after) - set(before))
    assert not changed, f"module-level state mutated: {changed}"
    assert not added, f"module-level names appeared at runtime: {added}"


def test_the_same_seed_gives_the_same_answer_twice_in_one_process():
    """Determinism, and a second net for accumulating hidden state.

    If anything cached across runs (a majorant, a lookup table mutated in
    place, an RNG shared by accident), the second run would differ.
    """
    from PIC_MCC_3D.MODULES import checkpoint as CP
    from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
    from PIC_MCC_3D.MODULES.poisson3d import PERIODIC
    from PIC_MCC_3D.MODULES.species3d import (ParticleSpecies, Species,
                                              SpeciesSet)

    def once():
        L, nx = 1e-2, 10
        pop = ParticleSpecies.cold_wave(Species.electron(), 12 ** 3,
                                        [0, 0, 0], [L] * 3, 1e14,
                                        amplitude=4e-5, mode=1, axis=2)
        sim = ESSim(nx, nx, nx, L / nx, SpeciesSet([pop]), 1e-11, bc=PERIODIC,
                    background_charge_density=uniform_background_for(pop, L**3),
                    seed=5)
        sim.run(5)
        return CP.state_digest(sim)

    assert once() == once()
