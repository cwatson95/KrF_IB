"""Verify the package: each module alone, then all of them together.

    python3 SPECK/PIC_MCC_3D/verify.py            # everything
    python3 SPECK/PIC_MCC_3D/verify.py --quick    # skip the isolation sweep

The ordinary suite (``pytest SPECK/PIC_MCC_3D/tests``) runs everything in ONE
process, which is exactly the condition under which cross-contamination
hides: a module that leaks state, or a test file that only passes because an
earlier one ran first, still shows green. This script separates the two
questions the suite cannot:

  1. INDEPENDENTLY -- each module imports in a fresh interpreter with no
     import-time side effects, and each test file passes ALONE in its own
     pytest process.
  2. TOGETHER -- the combined suite matches the sum of the isolated runs, the
     suite passes again in REVERSED file order, and the composition ladder
     (tests/test_integration.py) builds the full stack a rung at a time.

If (1) passes and (2) fails, the modules are fine and the coupling is not.
If a file passes in the suite but fails alone, it has an undeclared
dependency on another file. Both were real failure modes in this repo:
PLAN (74) traced 10 of 17 suite failures to leaked module globals.

Exit code is 0 only if every check passes, so this is usable as a gate.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKG_PARENT = HERE.parent
MODULES_DIR = HERE / "MODULES"
TESTS_DIR = HERE / "tests"

G = "\033[32m"; R = "\033[31m"; Y = "\033[33m"; D = "\033[2m"; Z = "\033[0m"


def _hdr(n, title):
    print(f"\n{D}{'─' * 72}{Z}\n[{n}] {title}\n{D}{'─' * 72}{Z}")


def _run(args, **kw):
    return subprocess.run(args, capture_output=True, text=True,
                          cwd=str(PKG_PARENT), **kw)


def module_names():
    return sorted(p.stem for p in MODULES_DIR.glob("*.py")
                  if p.stem != "__init__")


def test_files():
    return sorted(p for p in TESTS_DIR.glob("test_*.py"))


# --------------------------------------------------------------- the checks
def check_layers():
    """Topologically sort the intra-package import graph."""
    _hdr(1, "IMPORT GRAPH — acyclic and layered")
    deps = {}
    for m in module_names():
        tree = ast.parse((MODULES_DIR / f"{m}.py").read_text())
        d = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.ImportFrom) and n.level == 1:
                if n.module:
                    d.add(n.module)
                else:
                    d.update(a.name for a in n.names)
        deps[m] = d
    rem, layer, ok = dict(deps), 0, True
    while rem:
        ready = sorted(m for m, d in rem.items() if not (d & set(rem)))
        if not ready:
            print(f"  {R}CYCLE among {sorted(rem)}{Z}")
            ok = False
            break
        print(f"  layer {layer}: {', '.join(ready)}")
        for m in ready:
            rem.pop(m)
        layer += 1
    print(f"  {G}acyclic, {layer} layers{Z}" if ok else f"  {R}FAILED{Z}")
    return ok


def check_module_isolation():
    """Each module imported alone, in its own interpreter."""
    _hdr(2, "MODULE ISOLATION — each module imports alone, no side effects")
    ok = True
    for m in module_names():
        t0 = time.time()
        r = _run([sys.executable, "-c",
                  f"import PIC_MCC_3D.MODULES.{m} as M; print('ok')"])
        good = r.returncode == 0 and r.stdout.strip() == "ok" and not r.stderr
        ok &= good
        mark = f"{G}ok{Z}" if good else f"{R}FAIL{Z}"
        print(f"  {m:<18} {mark}  {D}{time.time() - t0:.2f}s{Z}")
        if not good:
            print(f"    {R}{(r.stderr or r.stdout).strip().splitlines()[-1:]}{Z}")
    return ok


_COUNT = re.compile(r"(\d+) passed")
_FAIL = re.compile(r"(\d+) failed")


def _pytest(args):
    t0 = time.time()
    r = _run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
              *args])
    passed = int(_COUNT.search(r.stdout).group(1)) if _COUNT.search(r.stdout) else 0
    failed = int(_FAIL.search(r.stdout).group(1)) if _FAIL.search(r.stdout) else 0
    return passed, failed, time.time() - t0, r


def check_files_alone():
    """Every test file in its OWN pytest process."""
    _hdr(3, "PER-FILE GATES — each test file alone in its own process")
    total, ok = 0, True
    per_file = {}
    for f in test_files():
        rel = f.relative_to(PKG_PARENT)
        p, fl, dt, r = _pytest([str(rel)])
        total += p
        per_file[f.name] = p
        good = fl == 0 and r.returncode == 0
        ok &= good
        mark = f"{G}{p:>3} passed{Z}" if good else f"{R}{fl} FAILED{Z}"
        print(f"  {f.name:<32} {mark}  {D}{dt:5.1f}s{Z}")
    print(f"  {D}isolated total: {total}{Z}")
    return ok, total, per_file


def check_together(isolated_total):
    """The whole suite in one process, and again in reversed order."""
    _hdr(4, "TOGETHER — one process, then reversed file order")
    p, fl, dt, r = _pytest([str(TESTS_DIR.relative_to(PKG_PARENT))])
    ok = fl == 0 and r.returncode == 0
    print(f"  combined                         "
          f"{(G + str(p) + ' passed' + Z) if ok else (R + str(fl) + ' FAILED' + Z)}"
          f"  {D}{dt:5.1f}s{Z}")

    same = (p == isolated_total)
    print(f"  isolated {isolated_total} vs combined {p}: "
          + (f"{G}equal — no test is order-dependent{Z}" if same
             else f"{R}MISMATCH — a test only runs in one of the two modes{Z}"))

    rev = [str(f.relative_to(PKG_PARENT)) for f in reversed(test_files())]
    p2, fl2, dt2, r2 = _pytest(rev)
    ok2 = fl2 == 0 and r2.returncode == 0
    print(f"  reversed file order              "
          f"{(G + str(p2) + ' passed' + Z) if ok2 else (R + str(fl2) + ' FAILED' + Z)}"
          f"  {D}{dt2:5.1f}s{Z}")
    if not ok2:
        print(f"    {Y}order dependence: a test relies on an earlier file{Z}")
    return ok and same and ok2


def check_ladder():
    """The composition ladder, reported rung by rung."""
    _hdr(5, "COMPOSITION LADDER — modules together, one rung at a time")
    ok = True
    for rung in ("rung1", "rung2", "rung3", "rung4", "rung5"):
        p, fl, dt, r = _pytest([str((TESTS_DIR / "test_integration.py")
                                    .relative_to(PKG_PARENT)), "-k", rung])
        good = fl == 0 and p >= 1
        ok &= good
        name = {"rung1": "fields only",
                "rung2": "+ kinetic particles",
                "rung3": "+ gas, cross sections, collisions",
                "rung4": "+ beam injection & preionization",
                "rung5": "+ persistence (full-stack exact restart)"}[rung]
        mark = f"{G}ok{Z}" if good else (f"{Y}skipped{Z}" if p == 0 and fl == 0
                                         else f"{R}FAIL{Z}")
        print(f"  {rung}  {name:<44} {mark}  {D}{dt:4.1f}s{Z}")
    return ok


def main(argv=None):
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--quick", action="store_true",
                   help="skip the per-file isolation sweep (the slow part)")
    ns = a.parse_args(argv)

    print(f"{D}PIC_MCC_3D verification — {sys.executable}{Z}")
    results = {"import graph": check_layers(),
               "module isolation": check_module_isolation()}

    if ns.quick:
        p, fl, dt, r = _pytest([str(TESTS_DIR.relative_to(PKG_PARENT))])
        _hdr(3, "TOGETHER — combined suite only (--quick)")
        print(f"  {p} passed, {fl} failed  {D}{dt:.1f}s{Z}")
        results["combined suite"] = (fl == 0 and r.returncode == 0)
    else:
        ok_files, total, _ = check_files_alone()
        results["per-file alone"] = ok_files
        results["together"] = check_together(total)
    results["composition ladder"] = check_ladder()

    _hdr(6, "SUMMARY")
    for k, v in results.items():
        print(f"  {k:<24} {(G + 'PASS' + Z) if v else (R + 'FAIL' + Z)}")
    good = all(results.values())
    print(f"\n{(G + 'ALL CHECKS PASS' + Z) if good else (R + 'VERIFICATION FAILED' + Z)}")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
