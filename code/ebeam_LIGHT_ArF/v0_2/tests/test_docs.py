"""The accuracy overview (LIGHT_ArF_ebeam/docs/accuracy_overview.tex) is current.

docs/make_overview.py writes every model number the overview prints. These
tests fail when a generated file no longer matches the network or the
scores, when the overview uses a number or table no file defines, when the
scores were made on other chemistry, or when a short source names someone
its row's citation does not. The energy ladder must match its generator,
and the system map's every node must balance on a fresh run of the model.
The runs file is checked by 1Dv0_2/tests/test_docs.py (it needs Fluid1D
and ~20 min), and the map's numbers file with LIGHT_DOCS_SLOW=1 (its bands
take ~7 min).
"""
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import pytest

LIGHT_ARF = Path(__file__).resolve().parents[2]
DOCS = LIGHT_ARF / "docs"
GENERATED = ("overview_rates.tex", "overview_scores.tex", "overview_runs.tex")
CHEM = ("model.py", "network.py", "species.py", "xsec.py",
        "data/ar_tail_ratio.txt", "beam.py")     # beam.py: v0_2's drive only


def _load(name):
    spec = importlib.util.spec_from_file_location(name, DOCS / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def mo():
    return _load("make_overview")


def test_rates_file_is_current(mo):
    assert mo.OUT_RATES.read_text() == mo.rates_tex(), \
        "rerun docs/make_overview.py"


def test_scores_file_is_current(mo):
    assert mo.OUT_SCORES.read_text() == mo.scores_tex(), \
        "rerun docs/make_overview.py"


@pytest.mark.parametrize("record, pkg", [("score_0d.json", "v0_2"),
                                         ("score_1d.json", "1Dv0_2")])
def test_scores_were_made_on_this_chemistry(record, pkg):
    rec = json.loads((DOCS / "arf_data" / record).read_text())
    now = {f: hashlib.md5((LIGHT_ARF / pkg / f).read_bytes()).hexdigest()
           for f in CHEM if (LIGHT_ARF / pkg / f).is_file()}
    assert rec["chemistry_md5"] == now, \
        f"{record} is from other chemistry: rerun arf_data/score_arf.py"


def test_every_number_and_table_the_overview_uses_is_defined():
    tex = (DOCS / "accuracy_overview.tex").read_text()
    gen = "".join((DOCS / f).read_text() for f in GENERATED)
    used = set(re.findall(r"\\wi\{([^}]+)\}", tex))
    defined = set(re.findall(r"\\widef\{([^}]+)\}", gen))
    assert used and not used - defined, sorted(used - defined)
    tables = set(re.findall(r"\\wiT([A-Za-z]+)", tex))
    tdef = re.findall(r"\\def\\wiT([A-Za-z]+)\{", gen)
    assert tables and not tables - set(tdef), sorted(tables - set(tdef))
    assert len(tdef) == len(set(tdef)), "a table name defined twice"
    # the map's and the ladder's numbers, in the overview and on the map
    wm = set(re.findall(r"\\wmdef\{([^}]+)\}",
                        (DOCS / "system_map_numbers.tex").read_text()))
    lad = set(re.findall(r"\\laddef\{([^}]+)\}",
                         (DOCS / "energy_ladder.tex").read_text()))
    page = (DOCS / "system_map.tex").read_text()
    for text, macro, defined in ((tex, "wm", wm), (page, "wm", wm),
                                 (tex, "lad", lad)):
        used = set(re.findall(r"\\" + macro + r"\{([^}\\]+)\}", text))
        if text is page:                  # \pc{key} prints \wm{key} %
            used |= set(re.findall(r"\\pc\{([^}\\]+)\}", text))
            # \ledger{y}{color}{label}{key}{options} passes its key as #4
            used |= set(re.findall(
                r"\\ledger\{[^}]*\}\{[^}]*\}\{[^}]*\}\{([^}\\]+)\}", text))
        used = {k for k in used if not k.startswith("#")}  # macro parameters
        assert used and not used - defined, sorted(used - defined)


def test_the_figures_exist():
    tex = (DOCS / "accuracy_overview.tex").read_text()
    figs = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", tex)
    figs += re.findall(r"\\includepdf(?:\[[^]]*\])?\{([^}]+)\.pdf\}", tex)
    assert figs and "system_map" in figs
    for f in figs:
        assert (DOCS / f"{f}.pdf").is_file(), f


def test_the_ladder_is_current_and_draws_every_row():
    ml = _load("make_energy_ladder")
    assert ml.OUT.read_text() == ml.render(), "rerun make_energy_ladder.py"


def test_the_system_map_balances_on_a_fresh_run():
    sm = _load("make_system_map")
    f, _ = sm.flows()                 # raises if a node does not balance
    nodes = sm.nodes(f)
    assert len(nodes) >= 30
    assert max(abs(a - b) for a, b in nodes.values()) < 1e-8
    text = (DOCS / "system_map_numbers.tex").read_text()
    assert f"\\wmdef{{nodes}}{{{len(nodes)}}}" in text


@pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                    reason="about 7 min (the bands): set LIGHT_DOCS_SLOW=1")
def test_the_system_map_numbers_are_current():
    sm = _load("make_system_map")
    assert sm.OUT.read_text() == sm.render(sm.numbers()), \
        "rerun make_system_map.py"


def _names(s):
    """The capitalised names in a short source (math and macros dropped)."""
    s = re.sub(r"\$[^$]*\$", " ", s)
    s = re.sub(r"\\[A-Za-z]+(\{\})?", " ", s)
    out = set()
    for t in re.findall(r"[A-Z][A-Za-z'\-]+", s):
        t = t.removesuffix("'s")
        if len(t) >= 3 and t not in ("The", "Table"):
            out.add(t)
    return out


def test_short_sources_name_only_whom_the_citations_name(mo):
    from v0_2 import network as N
    from v0_2.species import SPECIES
    rows = N.build() + N.electron_group_rows()
    for key, short in mo.SHORT.items():
        members = [r for r in rows if mo.short_key(r.rid) == key]
        assert members, key
        for r in members:
            text = r.cite + " ".join(N.BANDS[b].cite for b in N.bands_of(r))
            missing = {t for t in _names(short) if t not in text}
            assert not missing, (r.rid, missing)
    for name, short in mo.BAND_SHORT.items():
        missing = {t for t in _names(short) if t not in N.BANDS[name].cite}
        assert not missing, (name, missing)
    for name, (_, short) in mo.SPECIES_SHORT.items():
        missing = {t for t in _names(short) if t not in SPECIES[name].cite}
        assert not missing, (name, missing)
    assert set(mo.BAND_SHORT) == set(N.BANDS)
    assert set(mo.SPECIES_SHORT) == set(SPECIES)
