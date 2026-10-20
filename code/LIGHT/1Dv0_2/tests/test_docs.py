"""1Dv0_2's accuracy overview (LIGHT/docs/accuracy_overview_1Dv0_2.tex) must
show what the model computes today.

The overview takes the chemistry's scores from the v0_2 overview's
generated files, since 1Dv0_2's chemistry is v0_2's file for file, and
v0_2/tests/test_docs.py holds those files current. Here: the chemistry is
still v0_2's, the overview prints only numbers its files define, and, with
Fluid1D (Python >= 3.11 and numba, as test_1d.py needs), the 1-D numbers
against a fresh run, the slow ones with LIGHT_DOCS_SLOW=1 (about 10
minutes)."""

import importlib
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="about 10 minutes (v0_2's reference point "
                          "on four meshes, Tahoma, the time step): set "
                          "LIGHT_DOCS_SLOW=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents
                  if (p / "run_1Dv0_2.py").is_file()), None)
    if light is None or not (light / "docs" /
                             "make_overview_1Dv0_2.py").is_file():
        pytest.skip("LIGHT/docs is not beside 1Dv0_2 (standalone copy)")
    return light / "docs"


def _generator():
    docs = _docs_dir()
    spec = importlib.util.spec_from_file_location(
        "make_overview_1Dv0_2", docs / "make_overview_1Dv0_2.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _needs_fluid1d():
    if sys.version_info < (3, 11):
        pytest.skip("Fluid1D needs Python >= 3.11")
    pytest.importorskip("numba")
    bridge = importlib.import_module("1Dv0_2.fluid1d")
    try:
        bridge.src_dir()
    except ImportError as exc:                   # no Fluid1D beside LIGHT
        pytest.skip(str(exc))


def _check(which):
    _needs_fluid1d()
    script = _docs_dir() / "make_overview_1Dv0_2.py"
    if not (script.parents[2] / "ML" / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script), "--check", which],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_chemistry_is_v0_2s():
    """1Dv0_2's chemistry files have v0_2's md5 sums, so the v0_2
    overview's scores are 1Dv0_2's. A failure means one of the two changed:
    then the 1-D overview can no longer quote them."""
    gen = _generator()
    md5 = gen.check_chemistry()
    sim1d = importlib.import_module("1Dv0_2.sim1d")
    assert set(sim1d.CHEMISTRY_FILES) <= set(md5)
    assert "report.py" in md5


def test_overview_cites_only_numbers_its_files_define():
    """Every \\di, \\wi, \\wm and \\lad number the 1-D overview prints is
    defined in the file it comes from, and so is every table and
    profile."""
    docs = _docs_dir()
    tex = (docs / "accuracy_overview_1Dv0_2.tex").read_text()
    files = {"di": ("overview_1Dv0_2_fast.tex", "overview_1Dv0_2_slow.tex"),
             "wi": ("overview_v0_2_rates.tex", "overview_v0_2_fast.tex",
                    "overview_v0_2_slow.tex"),
             "wm": ("system_map_v0_2_numbers.tex",),
             "lad": ("energy_ladder_v0_2.tex",)}
    generated = ""
    for macro, names in files.items():
        have = set()
        for name in names:
            text = (docs / name).read_text()
            generated += text
            have |= set(re.findall(r"\\" + macro + r"def\{([^}]*)\}", text))
        used = set(re.findall(r"\\" + macro + r"\{([^}]*)\}", tex))
        assert used, macro
        assert used <= have, (macro, sorted(used - have))
    for table in ("diMesh", "diTahoma", "diProfA", "diProfB", "diProfC",
                  "diProfD", "diNe", "diTe", "wirBasic", "wirGroups",
                  "wirXsec", "wirBands", "wirSpecies"):
        assert "\\" + table in tex and "\\def\\" + table + "{" in generated, \
            table


def test_overview_1d_fast_numbers_are_current():
    """The chemistry files and the stage's rates against v0_2's model at
    random cell states, against a fresh run (seconds)."""
    _check("fast")


@SLOW
def test_overview_1d_slow_numbers_are_current():
    """The runs: v0_2's reference point on four meshes, Kr2* + F at the top
    of its band, the time step, and Tahoma from three seeds."""
    _check("slow")
