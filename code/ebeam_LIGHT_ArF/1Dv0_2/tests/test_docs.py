"""The accuracy overview's runs file
(LIGHT_ArF_ebeam/docs/overview_runs.tex) is current:
docs/make_overview.py's fresh runs of v0_2 and 1Dv0_2 give it
again, byte for byte. About 20 minutes, so only with LIGHT_DOCS_SLOW=1; it
needs Fluid1D (Python >= 3.11 with numba), as the 1-D tests do.

    LIGHT_DOCS_SLOW=1 ~/.venvs/f2sweep/bin/python -m pytest \\
        LIGHT_ArF_ebeam/1Dv0_2/tests/test_docs.py -q
"""
import importlib.util
import os
import sys
from pathlib import Path

import pytest

if sys.version_info < (3, 11):
    pytest.skip("Fluid1D needs Python >= 3.11", allow_module_level=True)
pytest.importorskip("numba")

DOCS = Path(__file__).resolve().parents[2] / "docs"


@pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                    reason="about 20 min: set LIGHT_DOCS_SLOW=1")
def test_runs_file_is_current():
    spec = importlib.util.spec_from_file_location(
        "make_overview", DOCS / "make_overview.py")
    mo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mo)
    assert mo.OUT_RUNS.read_text() == mo.runs_tex(write_figure=False), \
        "rerun docs/make_overview.py --runs"
