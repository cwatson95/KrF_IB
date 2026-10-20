"""v0_4's accuracy overview (LIGHT/docs/accuracy_overview_v0_4.tex) and
system map (system_map_v0_4.tex) must show what the model computes today:
every number they print is defined in the generated file it comes from,
the energy ladder is v0_4's, the generated numbers match a fresh run, and
the overview's words are held to its numbers.

ML fixes its LIGHT version when first imported, so every check that needs
ML runs its generator's --check in a process of its own, which sets
LIGHT_VERSION=v0_4. The fast numbers take about fifteen minutes (Tahoma
and its variants under v0_2's closure, the wave, the walls, the solver's
own checks); the slow ones, the reference point and Tahoma under the
kinetic EEDF, about half an hour on two processes."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="minutes (the reference point and Tahoma "
                          "under the kinetic EEDF, the map's run): set "
                          "LIGHT_DOCS_SLOW=1")
FAST = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_FAST"),
                          reason="about fifteen minutes (Tahoma's variants "
                          "and the wave under v0_2's closure): set "
                          "LIGHT_DOCS_FAST=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_4.py").is_file()),
                 None)
    if light is None or not (light / "docs" /
                             "make_overview_v0_4.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0_4 (standalone copy)")
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
    env = dict(os.environ, LIGHT_VERSION="v0_4", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-W", "ignore", str(script),
                        "--check", *args],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def _defined(text, macro):
    return set(re.findall(r"\\" + macro + r"def\{([^}]*)\}", text))


def _used(text, macro):
    return set(re.findall(r"\\" + macro + r"\{([^}#\\]*)\}", text))


def _values(text, macro):
    return dict(re.findall(r"\\" + macro + r"def\{([^}]*)\}\{(.*)\}$", text,
                           flags=re.M))


def _num(s):
    """A generated value as a float: 5.6, -22.3, $1.2{\\times}10^{-5}$,
    1.2{\\times}10^{-5}, $10^{-12}$, 'below 1e-9'."""
    s = s.replace("$", "").replace("\\,", "").replace("below ", "")
    s = s.replace("\\ensuremath{", "")
    s = s.replace("{\\times}", "x").replace("^{", "^").replace("}", "")
    m = re.fullmatch(r"\s*(-?[\d.]+)?\s*x?\s*(?:10\^(-?\d+))?\s*", s)
    if m and (m.group(1) or m.group(2)):
        mant = float(m.group(1)) if m.group(1) else 1.0
        ex = int(m.group(2)) if m.group(2) else 0
        return mant * 10.0 ** ex
    return float(s)


def _overview(docs):
    return (docs / "accuracy_overview_v0_4.tex").read_text()


def _files(docs):
    v4 = "".join((docs / f).read_text() for f in
                 ("overview_v0_4_fast.tex", "overview_v0_4_slow.tex",
                  "overview_v0_4_unc.tex"))
    v3 = "".join((docs / f).read_text() for f in
                 ("overview_v0_3_fast.tex", "overview_v0_3_slow.tex"))
    v2 = "".join((docs / f).read_text() for f in
                 ("overview_v0_2_fast.tex", "overview_v0_2_slow.tex"))
    wm = (docs / "system_map_v0_4_numbers.tex").read_text()
    lad = (docs / "energy_ladder_v0_4.tex").read_text()
    return {"wf": v4, "wk": v3, "wi": v2, "wm": wm, "lad": lad}


def test_overview_cites_only_numbers_its_files_define():
    """Every \\wf (v0_4), \\wm (the map), \\wk (v0_3), \\wi (v0_2) and \\lad
    key the overview prints is defined in its generated file, every table
    it sets too, and the slow file is not the stub a build may have left."""
    docs = _docs_dir()
    tex = _overview(docs)
    files = _files(docs)
    gen = {"wf": "make_overview_v0_4.py", "wm": "make_system_map_v0_4.py",
           "wk": "make_overview_v0_3.py", "wi": "make_overview_v0_2.py",
           "lad": "make_energy_ladder.py v0_4"}
    for macro, text in files.items():
        missing = _used(tex, macro) - _defined(text, macro)
        assert not missing, f"undefined {macro} keys {sorted(missing)}: " \
            f"rerun {gen[macro]}"
    # the \pc{} shorthand is a map number in percent
    missing = _used(tex, "pc") - _defined(files["wm"], "wm")
    assert not missing, f"undefined map keys {sorted(missing)}"
    for table in re.findall(r"\\(wf[A-Z]\w*)", tex):
        assert f"\\def\\{table}{{" in files["wf"], f"{table}: rerun " \
            "make_overview_v0_4.py"
    for table in re.findall(r"\\(wm[A-Z]\w*)", tex):
        assert f"\\def\\{table}{{" in files["wm"], f"{table}: rerun " \
            "make_system_map_v0_4.py"
    for f, flag in (("overview_v0_4_slow.tex", "--slow"),
                    ("overview_v0_4_unc.tex", "--unc")):
        assert "STUB" not in (docs / f).read_text(), \
            f"{f} is a stub: run make_overview_v0_4.py {flag}"


def test_every_citation_has_its_reference_and_every_reference_is_cited():
    """The overview cites by key; every key has a \\bibitem and every
    \\bibitem is cited (no stray or hand-numbered references)."""
    tex = _overview(_docs_dir())
    keys = {k.strip() for c in re.findall(r"\\cite\{([^}]*)\}", tex)
            for k in c.split(",")}
    items = set(re.findall(r"\\bibitem\{([^}]*)\}", tex))
    assert keys and keys == items, (sorted(keys - items), sorted(items - keys))
    body = tex[:tex.index(r"\begin{thebibliography}")]
    assert not re.search(r"(?<!\d)\[\d+(, ?\d+)*\]\s", body.replace(
        "Table I [81, 82, 62, 84, 16]", "")), "a hand-numbered citation"


def test_the_tests_without_electrons_score_as_v0_2():
    """"The 10 kinetic datasets have no electrons and no Kr**, so v0_4
    scores them exactly as v0_2": every sc.* and xu.* number in v0_4's
    file equals v0_2's, as the generator refuses otherwise."""
    docs = _docs_dir()
    files = _files(docs)
    v4, v2 = _values(files["wf"], "wf"), _values(files["wi"], "wi")
    keys = [k for k in v4 if k.startswith("sc.") and k not in
            ("sc.same", "sc.train.lo", "sc.train.hi")] + ["xu.50", "xu.100"]
    assert len(keys) >= 12
    for k in keys:
        assert v4[k] == v2[k], f"{k}: v0_4 {v4[k]}, v0_2 {v2[k]}"
    assert int(v4["sc.same"]) == len(keys) - 2


def test_overview_prose_matches_its_numbers():
    """What the overview says in words, against the numbers it prints."""
    docs = _docs_dir()
    files = _files(docs)
    v4 = _values(files["wf"], "wf")
    v3 = _values(files["wk"], "wk")
    v2 = _values(files["wi"], "wi")
    wm = _values(files["wm"], "wm")

    def f(k):
        return _num(v4[k])
    # the bottom line: "within ... rms" on the training sets and the
    # held-out decays, the absorption scored, Kr2F* most of it, the
    # compilation value far off
    assert 0.0 < f("sc.train.lo") < f("sc.train.hi") < 10.0
    assert 0.0 < f("sc.xu_figs5_7") < 15.0
    assert 0.5 < f("tah.kin.alpha.x") < 1.0
    assert f("tah.k2f.share") > 50.0
    assert f("tah.comp.alpha.x") > 2.5
    assert 0.9 < f("tah.isat.x") < 1.2 and 0.8 < f("tah.eta.x") < 1.1
    # "the 5s <-> 5p electron cycle ... turns ... of the input into heat
    # and takes the 248-nm share from ... to ...": the share rises with
    # the rung, the cycle is a large share, the relaxation outruns the
    # radiation
    assert _num(wm["out.248"]) > _num(wm["om.v03"])
    assert 10.0 < _num(wm["cyc.e"]) < 60.0
    assert _num(wm["brp.KPS"]) > 90.0
    assert _num(wm["out.IR"]) < 1.0
    # "the kinetic closure moves v0_4's Tahoma numbers by about 1 %"
    for q in ("g0", "alpha", "net"):
        assert abs(f(f"tah.kin.{q}") / f(f"tah.default.{q}") - 1.0) < 0.02, q
    # "the cycle cools the electrons ... and lowers n_e, and the gain
    # rises by a few percent"
    assert f("tah.default.te") < f("tah.v03.te")
    assert f("tah.default.ne") < f("tah.v03.ne")
    assert 1.0 < f("tah.default.g0") / f("tah.v03.g0") < 1.1
    assert f("tah.nocycle.g0.x") < f("tah.default.g0.x")
    # "the rung's own bands move it by under 1 %"; the named bands move it
    for tag in ("kannari", "kpslo"):
        assert abs(f(f"tah.{tag}.g0") / f("tah.default.g0") - 1.0) < 0.01, tag
    for tag in ("tri", "iiyield", "krff", "iifc", "kr2f"):
        assert abs(f(f"tah.{tag}.g0") / f("tah.default.g0") - 1.0) > 0.03, tag
    # the walls: a small loss in a 1-cm gap, more in the thin gap, and
    # 1Dv0_2's line above it
    assert f("wl.share") < 1e-3 and f("wl.thin.share") > f("wl.share")
    assert f("wl.1d") > f("wl.share")
    # the wave: a second input of a few percent, the band moving it
    # little, the ions' share small at 1e15 and rising with n_e, the
    # collisional regime and the two definitions apart
    assert 0.0 < f("wv.share") < 20.0 and f("wv.band") < 1.0
    assert f("wv.mid.ei") < 1.0
    assert f("wv.ei.16.ghz") > f("wv.ei.15.ghz")
    assert f("hf.nuw") > 1.0 and f("hf.share") > 1.1
    assert f("hf.atoms") < 0.02 and f("hf.ions") < 0.02
    # the 5s <-> 5p pair in the equation: the printed law, detailed
    # balance, the field's share and the balance
    assert f("pop.k") < 1e-3 and f("pop.db") < 1e-4
    assert 0.99 < f("pop.short.ratio.lo") <= f("pop.short.ratio.hi") < 1.02
    assert f("pop.short.gap") < 1e-5
    # the kinetic reference point and Tahoma: reached, balanced
    assert v4["tah.kin.reached"] == "every"
    assert f("tah.kin.gap") < 1e-2 and f("ref.gap") < 1e-2
    assert 0.9 < f("ref.ratio.lo") <= f("ref.ratio.hi") < 1.1
    # "the largest misses": every other score of the default model is
    # closer to its data than the named ones
    named = [abs(f("sc.hakuta_fig2")), abs(f("sc.hakuta_fig3")),
             abs(f("sc.bs_tau_b")), 100.0 * abs(_num(v3["pure.vuv.x"]) - 1.0),
             100.0 * abs(f("tah.kin.alpha.x") - 1.0)]
    others = [abs(f(k)) for k in (
        "sc.eden_kr", "sc.eden_f2", "sc.ksu_fig1", "sc.ksu_fig2",
        "sc.xu_figs5_7", "sc.ys_ir300", "sc.ys_ratio_230_300")] + [
        100.0 * abs(f(k) - 1.0) for k in ("xu.50", "xu.100", "tah.kin.g0.x",
                                          "tah.kin.net.x", "tah.isat.x",
                                          "tah.eta.x")]
    assert max(others) < min(named)
    # the seed the text names is v0_2's own (the configuration test holds
    # ML to it)
    assert "tah.seed" in v2
    # the propagated uncertainty: quadrature <= linear, the default inside
    # the envelope, the errors as the factors say, the largest term a band,
    # and the words "of the size of its errors on the gain and the net gain,
    # and the absorption's error is not covered"
    for o in ("ref", "g0", "alpha", "net"):
        assert 0.0 < f(f"unc.{o}.quad") <= f(f"unc.{o}.lin"), o
        assert f(f"unc.{o}.lo") <= f(f"unc.{o}.y0") <= f(f"unc.{o}.hi"), o
    for o, k in (("g0", "tah.default.g0"), ("alpha", "tah.default.alpha"),
                 ("net", "tah.default.net")):
        assert abs(f(f"unc.{o}.y0") - f(k)) < 0.02, o
    assert abs(f("unc.ref.y0") - _num(wm["out.248"])) < 0.06
    assert abs(f("unc.g0.dev.pct") - 100.0 * (f("tah.default.g0.x") - 1.0)) < 1.5
    light = _docs_dir().parent
    if str(light) not in sys.path:
        sys.path.insert(0, str(light))
    gen = _load(_docs_dir(), "make_system_map_v0_4")
    assert v4["unc.ref.top1"] in set(gen.BAND_LABELS.values())
    # "the largest term everywhere but the absorption is <top1>"
    assert v4["unc.g0.top1"] == v4["unc.net.top1"] == v4["unc.ref.top1"]
    # "every Tahoma measurement lies within one quadrature unit"
    for o in ("g0", "alpha", "net"):
        assert f(f"unc.{o}.dev.quad") <= 1.0, o
        assert v4[f"unc.{o}.covers"] == "within", o
        assert f(f"unc.{o}.down") <= 0.0 <= f(f"unc.{o}.up"), o
    # "a third to two fifths of each output"
    for o in ("ref", "g0", "alpha", "net"):
        assert 30.0 <= f(f"unc.{o}.quad.pct") <= 42.0, o
    assert abs(f("tah.kin.g0") / f("tah.g0.meas") - f("tah.kin.g0.x")) < 0.006


def test_tahoma_configuration_is_the_one_the_overview_names():
    """"ML's Tahoma configuration is v0_2's, the default network with every
    band at its default and n_e free from a seed of" v0_2's, the kinetic
    EEDF, no ladder preset, no gap and no wave."""
    docs = _docs_dir()
    ml_dir = docs.parents[1] / "ML"
    if not (ml_dir / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    code = ("import sys; sys.path.insert(0, sys.argv[1]); "
            "from ml import forward, light; p = forward.tahoma_params(); "
            "m = light('model').Model(p); "
            "print(repr((p.network, p.groups, p.rates, p.ne_mode, "
            "p.ne_cm3, m.eedf, p.ladder, p.gap_cm, p.em_E_V_cm)))")
    env = dict(os.environ, LIGHT_VERSION="v0_4", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-c", code, str(ml_dir)], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    seed = _num(_values((docs / "overview_v0_2_fast.tex").read_text(),
                        "wi")["tah.seed"])
    assert r.stdout.strip() == repr(("basic", (), (), "free", seed,
                                     "kinetic", "off", None, 0.0))


def test_energy_ladder_is_current():
    """The overview's energy ladder (energy_ladder_v0_4.tex) against v0_4's
    own species and network. If this fails, rerun make_energy_ladder.py
    v0_4, then latexmk -pdf accuracy_overview_v0_4.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0_4.tex").read_text() == \
        gen.render("v0_4"), ("energy_ladder_v0_4.tex is stale: rerun "
                             "make_energy_ladder.py v0_4 and rebuild the PDF")


def test_system_map_tabulates_every_band_and_names_every_zero_row():
    """The map's bands table lists every band of v0_4's network, and the
    rows it prints as zero are rows of the network (its generator checks
    on its run that every row is drawn and that those rows are 0)."""
    docs = _docs_dir()
    light = docs.parent
    if str(light) not in sys.path:
        sys.path.insert(0, str(light))
    gen = _load(docs, "make_system_map_v0_4")
    from v0_4.network import BANDS, build
    rows = {r.rid for r in build()}
    assert set(gen.ZERO_ROWS) <= rows, sorted(set(gen.ZERO_ROWS) - rows)
    assert set(gen.BAND_ORDER) == set(BANDS), \
        f"bands not tabulated: {sorted(set(BANDS) ^ set(gen.BAND_ORDER))}"
    wm = _values((docs / "system_map_v0_4_numbers.tex").read_text(), "wm")
    assert int(wm["bd.n"]) == len(BANDS)
    assert int(wm["bd.nzero"]) == sum(1 for b in BANDS.values()
                                      if b.default == 0.0)


@FAST
def test_overview_fast_numbers_are_current():
    """ML's scorecard under v0_4 (which the generator holds equal to
    v0_2's), Tahoma and its variants under v0_2's closure, the wave, the
    walls and the solver's checks, against a fresh run."""
    docs = _docs_dir()
    _check(docs / "make_overview_v0_4.py", "fast")


@SLOW
def test_overview_slow_numbers_are_current():
    docs = _docs_dir()
    _check(docs / "make_overview_v0_4.py", "slow")


@SLOW
def test_overview_uncertainty_numbers_are_current():
    """Every band and every constant with a quoted error at its two ends,
    at the reference point and at Tahoma under v0_2's closure."""
    docs = _docs_dir()
    _check(docs / "make_overview_v0_4.py", "unc")


@SLOW
def test_system_map_numbers_are_current():
    """The map's run (the reference point under v0_2's closure, its bands
    and its optics) against a fresh one, in a process of its own."""
    docs = _docs_dir()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-W", "ignore",
                        str(docs / "make_system_map_v0_4.py"), "--check"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
