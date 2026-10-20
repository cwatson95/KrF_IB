"""v0_1's system map, its accuracy overview and the v0/v0_1 comparison
(LIGHT/docs) must show what the model computes today.

ML fixes its LIGHT version when first imported, and v0's own docs tests
import it bound to v0. So every check that needs ML runs its generator's
--check in a process of its own, which sets LIGHT_VERSION=v0_1."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="minutes (the MCMC, every formation link, "
                          "the kinetic EEDF loop, the grid): set "
                          "LIGHT_DOCS_SLOW=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_1.py").is_file()),
                 None)
    if light is None or not (light / "docs" /
                             "make_system_map_v0_1.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0_1 (standalone copy)")
    return light / "docs"


def _load(docs, name):
    spec = importlib.util.spec_from_file_location(name, docs / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _check(script, *args):
    """A generator's --check, run in a process of its own."""
    if not (script.parents[2] / "ML" / "ml").is_dir() and \
            not (script.parents[3] / "ML" / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    env = dict(os.environ, LIGHT_VERSION="v0_1", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script), "--check", *args],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_system_map_numbers_are_current():
    """The v0_1 map against a fresh run. A failure means the model changed:
    run make_system_map_v0_1.py, then latexmk -pdf system_map_v0_1.tex and
    accuracy_overview_v0_1.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_system_map_v0_1")
    current = (docs / "system_map_v0_1_numbers.tex").read_text()
    assert current == gen.render(gen.numbers()), (
        "system_map_v0_1_numbers.tex is stale: rerun "
        "make_system_map_v0_1.py and rebuild the PDFs")


def test_overview_rate_tables_are_current():
    """Appendix A of accuracy_overview_v0_1.pdf is the network v0_1 builds."""
    docs = _docs_dir()
    gen = _load(docs, "make_overview_v0_1")
    assert (docs / "overview_v0_1_rates.tex").read_text() == gen.rates_tex(), (
        "overview_v0_1_rates.tex is stale: rerun make_overview_v0_1.py")


def test_overview_fast_numbers_are_current():
    """The decays, closure, ML's scorecard with Tahoma and the pure-gas
    checks, against a fresh run on v0_1."""
    _check(_docs_dir() / "make_overview_v0_1.py", "fast")


@SLOW
def test_overview_slow_numbers_are_current():
    _check(_docs_dir() / "make_overview_v0_1.py", "slow")


@SLOW
def test_review_numbers_are_current():
    """docs/REVIEW/v0_vs_v0_1.pdf's own numbers (both versions at Tahoma,
    the reference point one change at a time, the notes' grid)."""
    _check(_docs_dir() / "REVIEW" / "make_review.py")


# a row id, a Biagi group or a state label ("K26", "C1a", "DIS_a3", "F2",
# "A1Pi_u"): its digits are a name, not a value
_ID = re.compile(r"\b[A-Z][A-Za-z]*\d[A-Za-z0-9_]*")
_NUM = re.compile(r"\d+(?:\.\d+)?")


def _quoted_numbers_match(short: str, cite: str, pool: str) -> list:
    """Numbers in a short source missing from its full citation. A year may
    come from another row's citation (the pool)."""
    bad = []
    for tok in _NUM.findall(_ID.sub(" ", short)):
        if tok in cite:
            continue
        if len(tok) == 4 and tok.startswith(("19", "20")) and tok in pool:
            continue
        bad.append(tok)
    return bad


def test_overview_sources_cover_every_row_and_quote_their_citations():
    """Every row v0_1 can run has a short source in Appendix A, and every
    number a short source quotes appears in the row's full citation in
    v0_1/network.py (or, for a year, in another row's). The species notes
    are held to v0_1/species.py the same way."""
    docs = _docs_dir()
    gen = _load(docs, "make_overview_v0_1")
    from v0_1.network import BASIC_GROUPS, build
    from v0_1.species import SPECIES
    basic = build(BASIC_GROUPS, "basic")
    v0 = gen.v0_differences()
    pool = " ".join(r.cite for r in basic + build(("all",), "v0"))
    problems = []
    for r in basic:
        short = gen.SOURCES.get(r.rid)
        if short is None:
            problems.append(f"{r.rid}: no short source")
            continue
        problems += [f"{r.rid}: {t}" for t in
                     _quoted_numbers_match(short, r.cite, pool)]
    both = {**gen.SOURCES, **gen.SOURCES_V0}
    for r in v0:
        short = both.get(r.rid)
        if short is None:
            if not r.rid.endswith("_C"):
                problems.append(f"v0 {r.rid}: no short source")
            continue
        problems += [f"v0 {r.rid}: {t}" for t in
                     _quoted_numbers_match(short, r.cite, pool)]
    for name, note in gen.SPECIES_NOTES.items():
        problems += [f"species {name}: {t}" for t in
                     _quoted_numbers_match(note, SPECIES[name].cite, pool)]
    assert not problems, problems


def test_energy_ladder_is_current():
    """The overview's energy ladder (energy_ladder_v0_1.tex) against v0_1's
    species and network today. A failure means one of them changed: run
    make_energy_ladder.py, then latexmk -pdf accuracy_overview_v0_1.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0_1.tex").read_text() == \
        gen.render("v0_1"), ("energy_ladder_v0_1.tex is stale: rerun "
                             "make_energy_ladder.py and rebuild the PDF")


def test_overview_cites_only_ladder_numbers_the_ladder_defines():
    """Every \\lad number the overview prints beside the ladder is in the
    ladder's file."""
    docs = _docs_dir()
    have = set(re.findall(r"\\laddef\{([^}]*)\}",
                          (docs / "energy_ladder_v0_1.tex").read_text()))
    used = set(re.findall(r"\\lad\{([^}]*)\}",
                          (docs / "accuracy_overview_v0_1.tex").read_text()))
    assert used and used <= have, sorted(used - have)


def test_overview_cites_only_numbers_its_files_define():
    """Every \\vi and \\vm number accuracy_overview_v0_1.tex prints, its
    bottom line included, is defined in the generated file it comes from."""
    docs = _docs_dir()
    tex = (docs / "accuracy_overview_v0_1.tex").read_text()
    files = {"vi": ("overview_v0_1_rates.tex", "overview_v0_1_fast.tex",
                    "overview_v0_1_slow.tex"),
             "vm": ("system_map_v0_1_numbers.tex",)}
    for macro, names in files.items():
        have = set()
        for name in names:
            have |= set(re.findall(r"\\" + macro + r"def\{([^}]*)\}",
                                   (docs / name).read_text()))
        used = {k for k in re.findall(r"\\" + macro + r"\{([^}]*)\}", tex)
                if not k.startswith(("\\", "#"))}
        assert used, macro
        assert used <= have, (macro, sorted(used - have))
