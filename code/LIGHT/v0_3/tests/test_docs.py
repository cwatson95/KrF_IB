"""v0_3's accuracy overview (LIGHT/docs/accuracy_overview_v0_3.tex) must
show what the model computes today: every number it prints is defined in
the generated file it comes from, its energy ladder is v0_3's, and its
generated numbers match a fresh run.

ML fixes its LIGHT version when first imported, so every check that needs
ML runs its generator's --check in a process of its own, which sets
LIGHT_VERSION=v0_3. The slow numbers take about twenty-five minutes on
fourteen processes (every Tahoma run under the kinetic EEDF), and the
overview's words are held to its numbers here as well."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_DOCS_SLOW"),
                          reason="minutes (Tahoma, its bands, pure Kr and "
                          "the reference point under the kinetic EEDF): "
                          "set LIGHT_DOCS_SLOW=1")


def _docs_dir():
    here = Path(__file__).resolve()
    light = next((p for p in here.parents if (p / "run_v0_3.py").is_file()),
                 None)
    if light is None or not (light / "docs" /
                             "make_overview_v0_3.py").is_file():
        pytest.skip("LIGHT/docs is not beside v0_3 (standalone copy)")
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
    env = dict(os.environ, LIGHT_VERSION="v0_3", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script), "--check", *args],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def _defined(text, macro):
    return set(re.findall(r"\\" + macro + r"def\{([^}]*)\}", text))


def _used(text, macro):
    return set(re.findall(r"\\" + macro + r"\{([^}#\\]*)\}", text))


def test_overview_cites_only_numbers_its_files_define():
    """Every \\wk (v0_3), \\wi (v0_2) and \\lad key the overview prints is
    defined in its generated file, and every table it sets too."""
    docs = _docs_dir()
    tex = (docs / "accuracy_overview_v0_3.tex").read_text()
    v3 = "".join((docs / f).read_text() for f in
                 ("overview_v0_3_fast.tex", "overview_v0_3_slow.tex"))
    v2 = "".join((docs / f).read_text() for f in
                 ("overview_v0_2_fast.tex", "overview_v0_2_slow.tex"))
    lad = (docs / "energy_ladder_v0_3.tex").read_text()
    for macro, files, gen in (("wk", v3, "make_overview_v0_3.py"),
                              ("wi", v2, "make_overview_v0_2.py"),
                              ("lad", lad, "make_energy_ladder.py")):
        missing = _used(tex, macro) - _defined(files, macro)
        assert not missing, f"undefined {macro} keys {sorted(missing)}: " \
            f"rerun {gen}"
    for table in re.findall(r"\\(wk[A-Z]\w*)", tex):
        assert f"\\def\\{table}{{" in v3, f"{table}: rerun " \
            "make_overview_v0_3.py"


def _values(text, macro):
    return dict(re.findall(r"\\" + macro + r"def\{([^}]*)\}\{(.*)\}$", text,
                           flags=re.M))


def _num(s):
    """A generated number as a float: 4.79, +2.7, \\ensuremath{-}1.8,
    $1.9{\\times}10^{-5}$ or $10^{-14}$."""
    s = s.replace("\\ensuremath{-}", "-").strip("$ ")
    m = re.fullmatch(r"([+-]?[\d.]+)\{\\times\}10\^\{(-?\d+)\}", s)
    if m:
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    m = re.fullmatch(r"10\^\{(-?\d+)\}", s)
    return 10.0 ** int(m.group(1)) if m else float(s)


def test_overview_prose_matches_its_numbers():
    """What the overview says in words, against the numbers it prints: the
    directions of the reference point's shares, the excess field power,
    the agreement with ML's loop, the bands, the pure-Kr gap, the thermal
    electrons and the bottom line's largest misses."""
    docs = _docs_dir()
    v3 = _values("".join((docs / f).read_text() for f in (
        "overview_v0_3_fast.tex", "overview_v0_3_slow.tex")), "wk")
    v2 = _values("".join((docs / f).read_text() for f in (
        "overview_v0_2_fast.tex", "overview_v0_2_slow.tex")), "wi")

    def f(k):
        return _num(v3[k])
    # "spend less on elastic heating and ionization and more on exciting
    # Kr, and more KrF(B) forms"
    assert f("ref.elastic") < f("ref.elastic.v02")
    assert f("ref.krplus") < f("ref.krplus.v02")
    assert f("ref.krstar") > f("ref.krstar.v02")
    assert f("ref.krfb") > f("ref.krfb.v02")
    # "the field brings ... times the input power. The excess is the
    # table's renewal ... Renewed every 5 ns, the excess is at most ..."
    assert f("ref.ratio.lo") > 1.0
    assert f("ref.ratio.dt.hi") < f("ref.ratio.hi")
    # the first nanosecond: "the bands barely notice ... since only ... of
    # the pulse's energy enters then", the steps drop no growth term, and
    # at Tahoma they miss LIGHT's mean energy where steady solutions did
    assert f("first.d.steady") < 1e-3 and f("first.d.maxwell") < 1e-3
    assert f("first.edep") < 1e-2
    assert f("first.mem.goff") == 0 < f("first.steady.goff")
    assert f("first.mem.ratio.lo") > 0.9 and f("first.mem.ratio.hi") < 1.1
    assert f("first.steady.ratio.hi") < 0.5
    assert f("first.tah.mem.unreached") < f("first.tah.steady.unreached")
    # "the kinetic one at every update, held by a field of at most ...,
    # since the beam's secondaries alone heat the electrons almost to it"
    assert v3["tah.reached"] == "every" and f("tah.en") < 1.0
    # "spans the most, as in v0_2"
    assert v3["tah.bd.widest"] == v2["tah.bd.widest"]
    # "the two agree", and the bottom line's "within"
    assert all(f(f"el.abs.{q}") < 1.0 for q in ("te", "ne", "g0",
                                                 "formation"))
    # "do not close v0_2's gap at low pump": still beyond the measurement's
    # error, on v0_2's side
    x3, x2 = f("pure.vuv.x"), _num(v2["pure.vuv.x"])
    assert (x3 - 1.0) * (x2 - 1.0) > 0.0
    assert abs(x3 - 1.0) > f("pure.vuv.err") / f("pure.vuv.meas")
    # "a few percent below the gas temperature" on the default cells
    assert 0.01 < 1.0 - f("und.600.mean") / f("und.32kt") < 0.1
    # the bottom line's largest misses: every other score of the default
    # model is closer to its data
    named = [abs(f("sc.hakuta_fig2")), abs(f("sc.hakuta_fig3")),
             abs(f("sc.bs_tau_b")), 100.0 * abs(x3 - 1.0)]
    others = [abs(f(k)) for k in (
        "sc.eden_kr", "sc.eden_f2", "sc.ksu_fig1", "sc.ksu_fig2",
        "sc.xu_figs5_7", "sc.ys_ir300", "sc.ys_ratio_230_300")] + [
        100.0 * abs(f(k) - 1.0) for k in ("xu.50", "xu.100", "tah.g0.x",
                                          "tah.isat.x", "tah.eta.x")]
    assert max(others) < min(named)


def test_tahoma_configuration_is_the_one_the_overview_names():
    """"ML's Tahoma configuration is v0_2's, the default network with every
    band at its default and n_e free from a seed of" v0_2's, "now with the
    kinetic EEDF in the loop"."""
    docs = _docs_dir()
    ml_dir = docs.parents[1] / "ML"
    if not (ml_dir / "ml").is_dir():
        pytest.skip("ML/ is not beside LIGHT (standalone copy)")
    code = ("import sys; sys.path.insert(0, sys.argv[1]); "
            "from ml import forward, light; p = forward.tahoma_params(); "
            "m = light('model').Model(p); "
            "print(repr((p.network, p.groups, p.rates, p.ne_mode, "
            "p.ne_cm3, m.eedf)))")
    env = dict(os.environ, LIGHT_VERSION="v0_3", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-c", code, str(ml_dir)], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    seed = _num(_values((docs / "overview_v0_2_fast.tex").read_text(),
                        "wi")["tah.seed"])
    assert r.stdout.strip() == repr(("basic", (), (), "free", seed,
                                     "kinetic"))


def test_energy_ladder_is_current():
    """The overview's energy ladder (energy_ladder_v0_3.tex) against v0_3's
    own species and network. If this fails, rerun make_energy_ladder.py
    v0_3, then latexmk -pdf accuracy_overview_v0_3.tex."""
    docs = _docs_dir()
    gen = _load(docs, "make_energy_ladder")
    assert (docs / "energy_ladder_v0_3.tex").read_text() == \
        gen.render("v0_3"), ("energy_ladder_v0_3.tex is stale: rerun "
                             "make_energy_ladder.py v0_3 and rebuild the PDF")


def test_overview_fast_numbers_are_current():
    """ML's scorecard under v0_3 (which the generator holds equal to
    v0_2's) and the kinetic solver's checks, against a fresh run."""
    docs = _docs_dir()
    _check(docs / "make_overview_v0_3.py", "fast")


@SLOW
def test_overview_slow_numbers_are_current():
    docs = _docs_dir()
    _check(docs / "make_overview_v0_3.py", "slow")
