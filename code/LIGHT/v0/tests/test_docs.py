"""The system maps and the accuracy overview (LIGHT/docs) must show what the
model computes today."""

import importlib.util
import os
import re
from pathlib import Path

import pytest


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0.py").is_file()),
                 None)
    if light is None or not (light / "docs" / "make_system_map.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0 (standalone copy)")
    return light / "docs"


def test_system_map_numbers_are_current():
    """Regenerate the map's numbers and compare byte for byte. A failure
    means the model changed and the PDF is stale: run
    LIGHT/docs/make_system_map.py, then latexmk -pdf system_map.tex.
    (Closest printed value to a rounding edge on 2026-09-27: 32.549,
    a 3e-5 relative margin against ~2e-6 cross-interpreter drift.)"""
    docs = _docs_dir()
    spec = importlib.util.spec_from_file_location(
        "make_system_map", docs / "make_system_map.py")
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    current = (docs / "system_map_numbers.tex").read_text()
    assert current == gen.render(gen.numbers()), (
        "system_map_numbers.tex is stale: rerun make_system_map.py and "
        "rebuild system_map.pdf")


def _load(docs, name):
    spec = importlib.util.spec_from_file_location(name, docs / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_basic_system_map_numbers_are_current():
    """The basic network's map (system_map_basic.pdf) against a fresh run.
    A failure means the model changed: run make_system_map_basic.py, then
    latexmk -pdf system_map_basic.tex and accuracy_overview.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_system_map_basic")
    current = (docs / "system_map_basic_numbers.tex").read_text()
    assert current == gen.render(gen.numbers()), (
        "system_map_basic_numbers.tex is stale: rerun "
        "make_system_map_basic.py and rebuild the PDFs")


def _overview(docs):
    if not (docs.parent.parent / "ML" / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    return _load(docs, "make_overview")


def test_overview_rate_tables_are_current():
    """Appendix A of accuracy_overview.pdf is the network LIGHT builds."""
    docs = _docs_dir()
    gen = _overview(docs)
    assert (docs / "overview_rates.tex").read_text() == gen.rates_tex(), (
        "overview_rates.tex is stale: rerun make_overview.py")


def test_overview_fast_numbers_are_current():
    """The accuracy numbers that take seconds (decays, closure, ML's
    scorecard with Tahoma, the pure-gas checks), against a fresh run."""
    docs = _docs_dir()
    gen = _overview(docs)
    assert (docs / "overview_fast.tex").read_text() == gen._render(
        gen.fast_numbers(), "fast"), (
        "overview_fast.tex is stale: rerun make_overview.py")


@pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                    reason="about 5 minutes (the MCMC, every formation "
                    "link, the kinetic EEDF loop): set LIGHT_DOCS_SLOW=1")
def test_overview_slow_numbers_are_current():
    docs = _docs_dir()
    gen = _overview(docs)
    n, tables = gen.slow_numbers()
    assert (docs / "overview_slow.tex").read_text() == gen._render(
        n, "slow", tables), "overview_slow.tex is stale: rerun make_overview.py"


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
    """Every row LIGHT can run has a short source in Appendix A, and every
    number a short source quotes appears in the row's full citation in
    network.py (or, for a year, in another row's). The species notes are
    held to species.py the same way."""
    docs = _docs_dir()
    gen = _overview(docs)
    from v0.network import BASIC_GROUPS, build
    from v0.species import SPECIES
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
    """The overview's energy ladder (energy_ladder_v0.tex) against v0's
    species and network today. A failure means one of them changed: run
    make_energy_ladder.py, then latexmk -pdf accuracy_overview.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0.tex").read_text() == \
        gen.render("v0"), ("energy_ladder_v0.tex is stale: rerun "
                             "make_energy_ladder.py and rebuild the PDF")


def test_overview_cites_only_ladder_numbers_the_ladder_defines():
    """Every \\lad number the overview prints beside the ladder is in the
    ladder's file."""
    docs = _docs_dir()
    have = set(re.findall(r"\\laddef\{([^}]*)\}",
                          (docs / "energy_ladder_v0.tex").read_text()))
    used = set(re.findall(r"\\lad\{([^}]*)\}",
                          (docs / "accuracy_overview.tex").read_text()))
    assert used and used <= have, sorted(used - have)


def test_overview_cites_only_numbers_its_files_define():
    """Every \\ov, \\lm and \\rv number accuracy_overview.tex prints, its
    bottom line included, is defined in the generated file it comes from."""
    docs = _docs_dir()
    tex = (docs / "accuracy_overview.tex").read_text()
    files = {"ov": ("overview_rates.tex", "overview_fast.tex",
                    "overview_slow.tex"),
             "lm": ("system_map_basic_numbers.tex",),
             "rv": ("REVIEW/review_numbers.tex",)}
    for macro, names in files.items():
        have = set()
        for name in names:
            have |= set(re.findall(r"\\" + macro + r"def\{([^}]*)\}",
                                   (docs / name).read_text()))
        used = {k for k in re.findall(r"\\" + macro + r"\{([^}]*)\}", tex)
                if not k.startswith(("\\", "#"))}
        assert used, macro
        assert used <= have, (macro, sorted(used - have))
