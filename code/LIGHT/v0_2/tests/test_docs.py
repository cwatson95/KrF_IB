"""v0_2's documents (LIGHT/docs) must show what the model computes today:
its system map, its accuracy overview, and the review that compares v0,
v0_1 and v0_2 (LIGHT/docs/REVIEW).

ML fixes its LIGHT version when first imported, and v0's own docs tests
import it bound to v0. So every check that needs ML runs its generator's
--check in a process of its own, which sets LIGHT_VERSION=v0_2 (and
make_review.py runs everything ML computes in processes of its own, one
LIGHT version each)."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="minutes (every version at Tahoma, v0_2's "
                          "accuracy tests, the reference point and the "
                          "grid): set LIGHT_DOCS_SLOW=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_2.py").is_file()),
                 None)
    if light is None or not (light / "docs" /
                             "make_overview_v0_2.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0_2 (standalone copy)")
    return light / "docs"


def _load(docs, name):
    spec = importlib.util.spec_from_file_location(name, docs / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _check(script, *args):
    """A generator's --check, run in a process of its own."""
    if not (script.parents[2] / "ML" / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    env = dict(os.environ, LIGHT_VERSION="v0_2", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script), "--check", *args],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def _defined(text, macro):
    return set(re.findall(r"\\" + macro + r"def\{([^}]*)\}", text))


def _used(text, *macros):
    """Keys a document passes to these macros, less \\foreach variables
    (\\wm{\\k}) and macro parameters (\\wm{#1}), whose values the compile
    itself checks."""
    out = set()
    for macro in macros:
        out |= set(re.findall(r"\\" + macro + r"\{([^}]*)\}", text))
    return {k for k in out if not k.startswith(("\\", "#"))}


def test_system_map_numbers_are_current():
    """The v0_2 map against a fresh run, with every band's two ends. A
    failure means the model changed: run make_system_map_v0_2.py, then
    latexmk -pdf system_map_v0_2.tex and the documents that include it."""
    docs = _docs_dir()
    gen = _load(docs, "make_system_map_v0_2")
    current = (docs / "system_map_v0_2_numbers.tex").read_text()
    assert current == gen.render(*gen.numbers()), (
        "system_map_v0_2_numbers.tex is stale: rerun "
        "make_system_map_v0_2.py and rebuild the PDFs")


def test_overview_rate_tables_are_current():
    """Appendix A of accuracy_overview_v0_2.pdf is the network v0_2 builds."""
    docs = _docs_dir()
    gen = _load(docs, "make_overview_v0_2")
    assert (docs / "overview_v0_2_rates.tex").read_text() == gen.rates_tex(), (
        "overview_v0_2_rates.tex is stale: rerun make_overview_v0_2.py")


def test_overview_fast_numbers_are_current():
    """The decays, the closure, ML's scorecard with Tahoma at both ends of
    every band, the Kr2* route and the pure-gas checks, against a fresh
    run on v0_2 (about 2 minutes)."""
    _check(_docs_dir() / "make_overview_v0_2.py", "fast")


@SLOW
def test_overview_slow_numbers_are_current():
    _check(_docs_dir() / "make_overview_v0_2.py", "slow")


def test_map_and_overview_cite_only_numbers_their_files_define():
    """Every \\wm and \\wi number the map and the overview print, and every
    generated table they set, is defined in the file it comes from."""
    docs = _docs_dir()
    wm = _defined((docs / "system_map_v0_2_numbers.tex").read_text(), "wm")
    wi = set()
    for name in ("overview_v0_2_rates.tex", "overview_v0_2_fast.tex",
                 "overview_v0_2_slow.tex"):
        wi |= _defined((docs / name).read_text(), "wi")
    # the map prints its numbers through \wm and \pc (\wm with a percent)
    used = _used((docs / "system_map_v0_2.tex").read_text(), "wm", "pc")
    assert used and used <= wm, sorted(used - wm)
    tex = (docs / "accuracy_overview_v0_2.tex").read_text()
    for macro, have in (("wm", wm), ("wi", wi)):
        used = _used(tex, macro)
        assert used and used <= have, (macro, sorted(used - have))
    generated = "".join((docs / n).read_text() for n in (
        "system_map_v0_2_numbers.tex", "overview_v0_2_rates.tex",
        "overview_v0_2_fast.tex", "overview_v0_2_slow.tex"))
    for table in ("wmBands", "wiTahBands", "wiLinks", "wirBasic",
                  "wirGroups", "wirVzero", "wirXsec", "wirBands",
                  "wirSpecies"):
        assert "\\" + table in tex and "\\def\\" + table + "{" in generated, \
            table


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
    """Every row v0_2 can run has a short source in Appendix A, and every
    number a short source quotes appears in the row's full citation in
    v0_2/network.py (or, for a year, in another row's). The species notes
    are held to v0_2/species.py the same way, and the bands table prints
    each band's own citation, so every band needs one."""
    docs = _docs_dir()
    gen = _load(docs, "make_overview_v0_2")
    from v0_2.network import BANDS, BASIC_GROUPS, build
    from v0_2.species import SPECIES
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
    problems += [f"band {b.name}: no citation" for b in BANDS.values()
                 if not b.cite]
    assert not problems, problems


def _review_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_2.py").is_file()),
                 None)
    if light is None or not (light / "docs" / "REVIEW" /
                             "make_review.py").is_file():
        pytest.skip("LIGHT/docs/REVIEW is not beside v0_2 (standalone copy)")
    return light / "docs" / "REVIEW"


def test_review_cites_only_numbers_its_files_define():
    """Every \\rv, \\ov, \\vi, \\lm and \\vm number the review prints is
    defined in the generated file it comes from, and its two tables too."""
    rev = _review_dir()
    docs = rev.parent
    tex = (rev / "v0_vs_v0_1_vs_v0_2.tex").read_text()
    files = {"rv": ("REVIEW/review_numbers.tex",),
             "ov": ("overview_fast.tex", "overview_slow.tex"),
             "vi": ("overview_v0_1_fast.tex", "overview_v0_1_slow.tex"),
             "lm": ("system_map_basic_numbers.tex",),
             "vm": ("system_map_v0_1_numbers.tex",)}
    for macro, names in files.items():
        have = set()
        for name in names:
            text = (docs / name).read_text()
            have |= set(re.findall(r"\\" + macro + r"def\{([^}]*)\}", text))
        used = set(re.findall(r"\\" + macro + r"\{([^}]*)\}", tex))
        assert used, macro
        assert used <= have, (macro, sorted(used - have))
    numbers = (rev / "review_numbers.tex").read_text()
    for table in ("rvBands", "rvGrid"):
        assert "\\" + table in tex and "\\def\\" + table + "{" in numbers


@SLOW
def test_review_numbers_are_current():
    """docs/REVIEW/v0_vs_v0_1_vs_v0_2.pdf's own numbers against a fresh
    run. A failure means a version changed: rerun make_review.py, then
    latexmk -pdf v0_vs_v0_1_vs_v0_2.tex."""
    script = _review_dir() / "make_review.py"
    if not (script.parents[3] / "ML" / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script), "--check"], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_energy_ladder_is_current():
    """The overview's energy ladder (energy_ladder_v0_2.tex) against v0_2's
    species and network today. A failure means one of them changed: run
    make_energy_ladder.py, then latexmk -pdf accuracy_overview_v0_2.tex and
    accuracy_overview_1Dv0_2.tex, which sets the same ladder."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0_2.tex").read_text() == \
        gen.render("v0_2"), ("energy_ladder_v0_2.tex is stale: rerun "
                             "make_energy_ladder.py and rebuild the PDF")


def test_overview_cites_only_ladder_numbers_the_ladder_defines():
    """Every \\lad number the overview prints beside the ladder is in the
    ladder's file."""
    docs = _docs_dir()
    have = set(re.findall(r"\\laddef\{([^}]*)\}",
                          (docs / "energy_ladder_v0_2.tex").read_text()))
    used = set(re.findall(r"\\lad\{([^}]*)\}",
                          (docs / "accuracy_overview_v0_2.tex").read_text()))
    assert used and used <= have, sorted(used - have)
