#!/usr/bin/env python3
"""REBUILD verification runner — each module alone, then together.

The PIC_MCC_3D verify.py discipline: the ordinary suite runs everything in
one process, which is exactly the condition under which cross-contamination
hides. This runner separates the questions:

  1. every engine module imports ALONE in a fresh interpreter, silently;
  2. every test file passes in its OWN pytest process;
  3. the combined suite passes in one process;
  4. the combined suite passes with the file order REVERSED;
  5. the composition ladder builds the stack a rung at a time, so a failure
     names the module that introduced it.

Exit code 0 only if everything passes.

    python3 SPECK/REBUILD/verify.py [--quick]
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PY = sys.executable

MODULE_NAMES = sorted(p.stem for p in (HERE / "MODULES").glob("*.py")
                      if p.stem != "__init__")
TEST_FILES = sorted(p.name for p in (HERE / "tests").glob("test_*.py"))


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=str(ROOT), capture_output=True,
                          text=True, **kw)


def section(title):
    print(f"\n=== {title} " + "=" * max(0, 60 - len(title)))


def main():
    quick = "--quick" in sys.argv
    failures = []

    section("1. module isolation (fresh interpreter each)")
    for m in MODULE_NAMES:
        r = run([PY, "-c", f"import REBUILD.MODULES.{m}"])
        ok = r.returncode == 0 and not r.stdout and not r.stderr
        print(f"  {'ok ' if ok else 'FAIL'} import REBUILD.MODULES.{m}")
        if not ok:
            failures.append(f"import {m}: {r.stderr.strip()[:200]}")

    if not quick:
        section("2. per-file gates (each file its own pytest process)")
        totals = []
        for f in TEST_FILES:
            r = run([PY, "-m", "pytest", f"REBUILD/tests/{f}", "-q",
                     "--no-header", "-p", "no:cacheprovider"])
            last = [ln for ln in r.stdout.splitlines() if ln.strip()][-1:]
            ok = r.returncode == 0
            print(f"  {'ok ' if ok else 'FAIL'} {f:<28s} "
                  f"{last[0] if last else ''}")
            if not ok:
                failures.append(f"per-file {f}")
            totals.append(last[0] if last else "")

    section("3. combined suite, one process")
    r = run([PY, "-m", "pytest", "REBUILD/tests", "-q", "--no-header",
             "-p", "no:cacheprovider"])
    tail = [ln for ln in r.stdout.splitlines() if ln.strip()][-1:]
    print("  " + (tail[0] if tail else "no output"))
    if r.returncode != 0:
        failures.append("combined suite")
        print(r.stdout[-2000:])

    if not quick:
        section("4. combined suite, reversed file order")
        rev = [f"REBUILD/tests/{f}" for f in reversed(TEST_FILES)]
        r = run([PY, "-m", "pytest", *rev, "-q", "--no-header",
                 "-p", "no:cacheprovider"])
        tail = [ln for ln in r.stdout.splitlines() if ln.strip()][-1:]
        print("  " + (tail[0] if tail else "no output"))
        if r.returncode != 0:
            failures.append("reversed order")

    section("5. composition ladder (one module per rung)")
    rungs = [
        ("units+species", "from REBUILD.MODULES.units import LOSCHMIDT_M3;"
         "from REBUILD.MODULES.species import StateLayout, KRF_MODE_SPECIES;"
         "lay = StateLayout(KRF_MODE_SPECIES, 2, 4);"
         "assert abs(LOSCHMIDT_M3/2.6867801118e25 - 1) < 1e-9"),
        ("+network", "from REBUILD.MODULES.network import build_network;"
         "rows = build_network('krf', {'Kr*','Kr**','Kr+','F-','F2(a3)',"
         "'F2(A1)','F2+','Ne*','Ne**','Ne+'}); assert len(rows) >= 80"),
        ("+tiers+chemistry", "import numpy as np;"
         "from REBUILD.MODULES.loki_tier import rates_for_mode;"
         "from REBUILD.MODULES.network import build_network, mode_optics;"
         "from REBUILD.MODULES.chemistry import CompiledNetwork;"
         "from REBUILD.MODULES.species import StateLayout, KRF_MODE_SPECIES;"
         "from REBUILD.MODULES.species import e_form_of;"
         "from REBUILD.MODULES.units import EV, photon_energy_J;"
         "lok = rates_for_mode('krf', {'Ar':0.949,'Kr':0.05,'F2':0.001});"
         "rows = build_network('krf', lok.available_buckets());"
         "lay = StateLayout(KRF_MODE_SPECIES, 2);"
         "cn = CompiledNetwork(rows, lay, lok, mode_optics('krf'));"
         "cn.set_photon_energy(photon_energy_J(248.0)/EV);"
         "rng = np.random.default_rng(1);"
         "dens = np.full((lay.n_species, 2), 1e16);"
         "dens[lay.index['Ar']] = 7e25; dens[lay.index['Kr']] = 3.7e24;"
         "dens[lay.index['F2']] = 7.4e22;"
         "te = np.array([1.0, 2.0]); tg = np.array([300.0, 400.0]);"
         "ne = lay.electron_density(dens);"
         "out = cn.rhs(dens, te, tg, ne, 1e17, 7.4e25);"
         "e_form = np.array([e_form_of(s) for s in lay.species]);"
         "d_pot = (e_form[:,None]*out['ddens']).sum(axis=0)*EV;"
         "hv = cn.hv_eV;"
         "tot = d_pot + out['du_e'] + out['dgas'] + out['rad_escape'] +"
         " out['r_se']*hv*EV + (out['r_ste']-out['r_x'])*hv*EV;"
         "scale = np.abs(d_pot).max() + np.abs(out['dgas']).max();"
         "assert np.abs(tot).max()/scale < 1e-12"),
        ("+engine (0-D run + ledger)", "from REBUILD.MODULES.engine import "
         "Device; from REBUILD.MODULES.params import engine_parameters;"
         "p = engine_parameters();"
         "p.set_many({'mode':'krf','drive':'beam','t_end_ns':30.0,'nx':1,"
         "'ny':1,'nz':1,'Lx_m':0.0141,'Ly_m':0.0141,'Lz_m':0.75,"
         "'beam_axis':0}, source='problem');"
         "r = Device(p).run(n_save=7); assert r['success'];"
         "worst = max(abs(c['residual_rel']) for c in r['closure'][1:]);"
         "assert worst < 5e-4, worst"),
        ("+checkpoint (roundtrip)", "import tempfile, numpy as np;"
         "from pathlib import Path;"
         "from REBUILD.MODULES.engine import Device;"
         "from REBUILD.MODULES.params import engine_parameters;"
         "from REBUILD.MODULES.checkpoint import save, restore;"
         "p = engine_parameters();"
         "p.set_many({'mode':'krf','drive':'beam','t_end_ns':30.0,'nx':1,"
         "'ny':1,'nz':1,'Lx_m':0.0141,'Ly_m':0.0141,'Lz_m':0.75,"
         "'beam_axis':0}, source='problem');"
         "d = Device(p); y = d.initial_state();"
         "f = Path(tempfile.mkdtemp())/'ck.npz'; save(f, d, 0.0, y);"
         "d2, t2, y2 = restore(f); assert np.array_equal(y, y2)"),
    ]
    for name, code in rungs:
        r = run([PY, "-c", code])
        ok = r.returncode == 0
        print(f"  {'ok ' if ok else 'FAIL'} rung: {name}")
        if not ok:
            failures.append(f"ladder {name}: {r.stderr.strip()[-300:]}")

    section("summary")
    if failures:
        for f in failures:
            print(f"  FAIL {f}")
        print(f"\n{len(failures)} failure(s)")
        return 1
    print("  ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
