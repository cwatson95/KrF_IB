"""v0_5's accuracy overview (LIGHT/docs/accuracy_overview_v0_5.tex) and
system map (system_map_v0_5.tex) must show what the model computes today:
every number they print is defined in the generated file it comes from
(v0_5's own numbers are \\wg, v0_4's generated files stay \\wf beside
them), the energy ladder is v0_5's, the generated numbers match a fresh
run, every citation has its reference, and the overview's words are held
to its numbers.

ML fixes its LIGHT version when first imported, so every check that needs
ML runs its generator's --check in a process of its own, which sets
LIGHT_VERSION=v0_5. The fast numbers take about twenty minutes (Tahoma
and its variants, the Kleikamp laser's scans, the wave, the walls, the
solver's checks); the slow ones, the reference point and Tahoma under the
kinetic EEDF, about half an hour on two processes; the propagated
uncertainty about ten minutes."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="minutes (the reference point and Tahoma "
                          "under the kinetic EEDF, the map's run, the "
                          "propagation): set LIGHT_DOCS_SLOW=1")
FAST = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_FAST"),
                          reason="about twenty minutes (Tahoma's variants, "
                          "the Kleikamp scans, the wave under v0_2's "
                          "closure): set LIGHT_DOCS_FAST=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_5.py").is_file()),
                 None)
    if light is None or not (light / "docs" /
                             "make_overview_v0_5.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0_5 (standalone copy)")
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
    env = dict(os.environ, LIGHT_VERSION="v0_5", PYTHONDONTWRITEBYTECODE="1")
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
    1.2{\\times}10^{-5}, $10^{-12}$, 'below 1e-9', \\ensuremath{-22.3}."""
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
    return (docs / "accuracy_overview_v0_5.tex").read_text()


def _files(docs):
    v5 = "".join((docs / f).read_text() for f in
                 ("overview_v0_5_fast.tex", "overview_v0_5_slow.tex",
                  "overview_v0_5_unc.tex"))
    v4 = "".join((docs / f).read_text() for f in
                 ("overview_v0_4_fast.tex", "overview_v0_4_slow.tex",
                  "overview_v0_4_unc.tex"))
    v3 = "".join((docs / f).read_text() for f in
                 ("overview_v0_3_fast.tex", "overview_v0_3_slow.tex"))
    v2 = "".join((docs / f).read_text() for f in
                 ("overview_v0_2_fast.tex", "overview_v0_2_slow.tex"))
    wm = (docs / "system_map_v0_5_numbers.tex").read_text()
    lad = (docs / "energy_ladder_v0_5.tex").read_text()
    return {"wg": v5, "wf": v4, "wk": v3, "wi": v2, "wm": wm, "lad": lad}


def test_overview_cites_only_numbers_its_files_define():
    """Every \\wg (v0_5), \\wf (v0_4), \\wm (the map), \\wk (v0_3), \\wi
    (v0_2) and \\lad key the overview prints is defined in its generated
    file, every table it sets too, and no v0_5 file is a stub a build may
    have left."""
    docs = _docs_dir()
    tex = _overview(docs)
    files = _files(docs)
    gen = {"wg": "make_overview_v0_5.py", "wf": "make_overview_v0_4.py",
           "wm": "make_system_map_v0_5.py", "wk": "make_overview_v0_3.py",
           "wi": "make_overview_v0_2.py", "lad": "make_energy_ladder.py v0_5"}
    for macro, text in files.items():
        missing = _used(tex, macro) - _defined(text, macro)
        assert not missing, f"undefined {macro} keys {sorted(missing)}: " \
            f"rerun {gen[macro]}"
    missing = _used(tex, "pc") - _defined(files["wm"], "wm")
    assert not missing, f"undefined map keys {sorted(missing)}"
    for table in re.findall(r"\\(wg[A-Z]\w*)", tex):
        assert f"\\def\\{table}{{" in files["wg"], f"{table}: rerun " \
            "make_overview_v0_5.py"
    for table in re.findall(r"\\(wm[A-Z]\w*)", tex):
        assert f"\\def\\{table}{{" in files["wm"], f"{table}: rerun " \
            "make_system_map_v0_5.py"
    for f, flag in (("overview_v0_5_fast.tex", "--fast"),
                    ("overview_v0_5_slow.tex", "--slow"),
                    ("overview_v0_5_unc.tex", "--unc")):
        assert "STUB" not in (docs / f).read_text(), \
            f"{f} is a stub: run make_overview_v0_5.py {flag}"


def test_every_citation_has_its_reference_and_every_reference_is_cited():
    tex = _overview(_docs_dir())
    keys = {k.strip() for c in re.findall(r"\\cite\{([^}]*)\}", tex)
            for k in c.split(",")}
    items = set(re.findall(r"\\bibitem\{([^}]*)\}", tex))
    assert keys and keys == items, (sorted(keys - items), sorted(items - keys))
    body = tex[:tex.index(r"\begin{thebibliography}")]
    assert not re.search(r"(?<!\d)\[\d+(, ?\d+)*\]\s", body.replace(
        "Table I [76]", "")), "a hand-numbered citation"


def test_the_tests_without_electrons_score_as_v0_2():
    """Every sc.* and xu.* number in v0_5's file equals v0_2's, as the
    generator refuses otherwise."""
    docs = _docs_dir()
    files = _files(docs)
    v5, v2 = _values(files["wg"], "wg"), _values(files["wi"], "wi")
    keys = [k for k in v5 if k.startswith("sc.") and k not in
            ("sc.same", "sc.train.lo", "sc.train.hi")] + ["xu.50", "xu.100"]
    assert len(keys) >= 12
    for k in keys:
        assert v5[k] == v2[k], f"{k}: v0_5 {v5[k]}, v0_2 {v2[k]}"
    assert int(v5["sc.same"]) == len(keys) - 2


def test_overview_prose_matches_its_numbers():
    """What v0_5's overview says in words, against the numbers it prints."""
    docs = _docs_dir()
    files = _files(docs)
    v5 = _values(files["wg"], "wg")
    v4 = _values(files["wf"], "wf")
    v3 = _values(files["wk"], "wk")
    v2 = _values(files["wi"], "wi")
    wm = _values(files["wm"], "wm")

    def f(k):
        return _num(v5[k])
    # the bottom line's scores
    assert 0.0 < f("sc.train.lo") < f("sc.train.hi") < 10.0
    assert 0.0 < f("sc.xu_figs5_7") < 15.0
    assert 0.5 < f("tah.kin.alpha.x") < 1.0 and f("tah.k2f.share") > 50.0
    assert 0.9 < f("tah.isat.x") < 1.2 and 0.8 < f("tah.eta.x") < 1.1
    # "the measurement confirms the fit": above it below 3 eV, below it
    # at 5.8 eV, within 40 % throughout; the share moves little
    for tag, lo, hi in (("1", 1.2, 1.6), ("116", 1.2, 1.6), ("3", 1.0, 1.25),
                        ("58", 0.8, 1.0)):
        assert lo < f(f"jung.ratio.{tag}") < hi, tag
    assert abs(_num(wm["out.248"]) - _num(wm["om.v04"])) < 0.3
    assert abs(f("ref.248") - _num(v4["ref.248"])) < 0.3
    # "below it at the reference point's Te, above it at Tahoma's": the
    # two_t share falls, the kinetic one rises
    assert _num(wm["out.248"]) < _num(wm["om.v04"])
    assert f("ref.248") > _num(v4["ref.248"])
    assert 0.8 < _num(wm["k.76.ratio"]) < 1.0
    # the switch: the v0_4 column equals v0_4's own two_t Tahoma numbers
    assert abs(f("tah.k76kannari.g0") - f("tah.v04.g0")) < 0.02
    assert abs(f("v04.248") - 5.787) < 0.02
    # "lifts the gain by about 1 %", the kinetic closure within 2 %
    assert 1.0 < f("tah.default.g0") / f("tah.v04.g0") < 1.03
    for q in ("g0", "alpha", "net"):
        assert abs(f(f"tah.kin.{q}") / f(f"tah.default.{q}") - 1.0) < 0.02, q
    # the band's ends bracket the default, the cycle off below them
    assert f("tah.k5s5plo.g0") < f("tah.default.g0") < f("tah.k5s5phi.g0")
    assert f("tah.nocycle.g0") < f("tah.k5s5plo.g0")
    # the propagation: quadrature <= linear, within the excursions, every
    # Tahoma measurement within one unit, the largest term a band, the
    # 5s-5p term no larger than v0_4's
    light = docs.parent
    if str(light) not in sys.path:
        sys.path.insert(0, str(light))
    gen = _load(docs, "make_system_map_v0_5")
    for o in ("ref", "g0", "alpha", "net"):
        assert 0.0 < f(f"unc.{o}.quad") <= f(f"unc.{o}.lin"), o
        assert f(f"unc.{o}.down") <= 0.0 <= f(f"unc.{o}.up"), o
    for o in ("g0", "alpha", "net"):
        assert f(f"unc.{o}.dev.quad") <= 1.0, o
        assert v5[f"unc.{o}.covers"] == "within", o
    assert v5["unc.ref.top1"] in set(gen.BAND_LABELS.values())
    assert v5["unc.g0.top1"] == v5["unc.net.top1"] == v5["unc.ref.top1"]
    assert f("unc.ref.quad") <= _num(v4["unc.ref.quad"]) + 0.05
    # the Kleikamp laser: the optimum below the measured, the saturated
    # proxy above the small-signal one, the pump spread ordered, the F2
    # optimum not reached, the trimer too strong, the sensitivity run near 2
    assert f("kl.kr.sat.lo") < f("kl.kr.sat.mid") < f("kl.kr.sat.hi") < f("kl.meas.kr")
    assert f("kl.kr.ss.mid") < f("kl.kr.sat.mid")
    assert f("kl.f2.sat") >= f("kl.f2.top") and f("kl.f2.ss") >= f("kl.f2.top")
    assert f("kl.k.ratio") > 1.5
    assert abs(f("kl.scaled.k") / f("kl.meas.kR") - 1.0) < 0.05
    assert f("kl.scaled.sat") > f("kl.kr.sat.mid")
    # the largest errors: every other score of the default model is closer
    named = [abs(f("sc.hakuta_fig2")), abs(f("sc.hakuta_fig3")),
             abs(f("sc.bs_tau_b")), 100.0 * abs(_num(v3["pure.vuv.x"]) - 1.0),
             100.0 * abs(f("tah.kin.alpha.x") - 1.0),
             100.0 * abs(f("kl.kr.sat.mid") / f("kl.meas.kr") - 1.0)]
    others = [abs(f(k)) for k in (
        "sc.eden_kr", "sc.eden_f2", "sc.ksu_fig1", "sc.ksu_fig2",
        "sc.xu_figs5_7", "sc.ys_ir300", "sc.ys_ratio_230_300")] + [
        100.0 * abs(f(k) - 1.0) for k in ("xu.50", "xu.100", "tah.kin.g0.x",
                                          "tah.kin.net.x", "tah.isat.x",
                                          "tah.eta.x")]
    assert max(others) < min(named)
    assert "tah.seed" in v2


def test_tahoma_configuration_is_the_one_the_overview_names():
    """"ML's Tahoma configuration is v0_2's, the default network with every
    band at its default and n_e free from a seed of" v0_2's, the kinetic
    EEDF, no ladder preset, no gap, no wave, the measured cross section."""
    docs = _docs_dir()
    ml_dir = docs.parents[1] / "ML"
    if not (ml_dir / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    code = ("import sys; sys.path.insert(0, sys.argv[1]); "
            "from ml import forward, light; p = forward.tahoma_params(); "
            "m = light('model').Model(p); "
            "print(repr((p.network, p.groups, p.rates, p.ne_mode, "
            "p.ne_cm3, m.eedf, p.ladder, p.gap_cm, p.em_E_V_cm, p.k76)))")
    env = dict(os.environ, LIGHT_VERSION="v0_5", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-c", code, str(ml_dir)], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    seed = _num(_values((docs / "overview_v0_2_fast.tex").read_text(),
                        "wi")["tah.seed"])
    assert r.stdout.strip() == repr(("basic", (), (), "free", seed,
                                     "kinetic", "off", None, 0.0, "jung"))


def test_energy_ladder_is_current():
    """The overview's energy ladder (energy_ladder_v0_5.tex) against v0_5's
    own species and network. If this fails, rerun make_energy_ladder.py
    v0_5, then latexmk -pdf accuracy_overview_v0_5.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0_5.tex").read_text() == \
        gen.render("v0_5"), ("energy_ladder_v0_5.tex is stale: rerun "
                             "make_energy_ladder.py v0_5 and rebuild the PDF")


def test_system_map_tabulates_every_band_and_names_every_zero_row():
    docs = _docs_dir()
    light = docs.parent
    if str(light) not in sys.path:
        sys.path.insert(0, str(light))
    gen = _load(docs, "make_system_map_v0_5")
    from v0_5.network import BANDS, build
    rows = {r.rid for r in build()}
    assert set(gen.ZERO_ROWS) <= rows, sorted(set(gen.ZERO_ROWS) - rows)
    assert set(gen.BAND_ORDER) == set(BANDS), \
        f"bands not tabulated: {sorted(set(BANDS) ^ set(gen.BAND_ORDER))}"
    wm = _values((docs / "system_map_v0_5_numbers.tex").read_text(), "wm")
    assert int(wm["bd.n"]) == len(BANDS)
    assert int(wm["bd.nzero"]) == sum(1 for b in BANDS.values()
                                      if b.default == 0.0)


@FAST
def test_overview_fast_numbers_are_current():
    docs = _docs_dir()
    _check(docs / "make_overview_v0_5.py", "fast")


@SLOW
def test_overview_slow_numbers_are_current():
    docs = _docs_dir()
    _check(docs / "make_overview_v0_5.py", "slow")


@SLOW
def test_overview_uncertainty_numbers_are_current():
    docs = _docs_dir()
    _check(docs / "make_overview_v0_5.py", "unc")


@SLOW
def test_system_map_numbers_are_current():
    docs = _docs_dir()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-W", "ignore",
                        str(docs / "make_system_map_v0_5.py"), "--check"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
