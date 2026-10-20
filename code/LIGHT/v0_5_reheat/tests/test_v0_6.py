"""v0_6: the EM wave's window and the twin pairs (reheat.py). The default
window must be v0_5's wave bit for bit (against v0_5 itself), a window
must cut the run where it opens and closes and book E_em only inside it,
the ledger must close with the wave alone in the afterglow under both
closures, a window past t_end must be no wave, the kinetic EEDF must see
the wave's field only in its window, and the twin pair's numbers must be
the arithmetic reheat.py states."""

import math

import numpy as np
import pytest
from pytest import approx

from v0_6 import reheat
from v0_6.model import Model, Params, run
from v0_6.report import accounting, row_dict, summary

from .test_v0_4 import _check_closures

FAST = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
WAVE = dict(em_E_V_cm=3.0e4, em_f_Hz=6.2e12)
AFTER = dict(em_t_on_s=0.3e-6, em_t_off_s=math.inf)


# ---------------------------------------------------------------------------
# the window
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kw", [
    dict(),
    dict(em_E_V_cm=1500.0, em_f_Hz=95e9),
    dict(em_E_V_cm=3000.0, em_f_Hz=1e12, power_W_cm3=0.0),
])
def test_default_window_is_v0_5_bit_for_bit(kw):
    v05 = pytest.importorskip("v0_5.model")
    a = v05.run(v05.Params(**FAST, **kw))
    b = run(Params(**FAST, **kw))
    assert np.array_equal(a.t, b.t)
    for s in a.n:
        assert np.array_equal(a.n[s], b.n[s]), s
    for k in a.ledger:
        assert np.array_equal(a.ledger[k], b.ledger[k]), k


def test_segments_cut_at_the_window():
    T, F = True, False
    seg = lambda **kw: Model(Params(**FAST, **kw)).segments()  # noqa: E731
    assert seg() == [(0.0, 0.3e-6, T, F), (0.3e-6, 0.6e-6, F, F)]
    assert seg(**WAVE) == [(0.0, 0.3e-6, T, T), (0.3e-6, 0.6e-6, F, F)]
    assert seg(**WAVE, **AFTER) == [(0.0, 0.3e-6, T, F),
                                    (0.3e-6, 0.6e-6, F, T)]
    assert seg(**WAVE, em_t_on_s=0.1e-6, em_t_off_s=0.2e-6) == [
        (0.0, 0.1e-6, T, F), (0.1e-6, 0.2e-6, T, T), (0.2e-6, 0.3e-6, T, F),
        (0.3e-6, 0.6e-6, F, F)]
    assert seg(**WAVE, em_t_on_s=0.3e-6, em_t_off_s=0.4e-6) == [
        (0.0, 0.3e-6, T, F), (0.3e-6, 0.4e-6, F, T), (0.4e-6, 0.6e-6, F, F)]
    # without a field the cuts stay and nothing is on: a twin with the
    # same window integrates the same pieces
    assert seg(em_t_on_s=0.3e-6, em_t_off_s=0.4e-6) == [
        (0.0, 0.3e-6, T, F), (0.3e-6, 0.4e-6, F, F), (0.4e-6, 0.6e-6, F, F)]
    # a decay run (no pulse) is one piece without either drive
    assert Model(Params(eedf="two_t", pulse_s=0.0, t_end_s=1e-7,
                        power_W_cm3=0.0)).segments() == [(0.0, 1e-7, F, F)]


def test_windows_that_cannot_be():
    with pytest.raises(ValueError):
        Model(Params(**FAST, **WAVE, em_t_on_s=-1e-9))
    # opened at the pulse's end with the default close (the pulse's end)
    with pytest.raises(ValueError, match="em_t_off_s None"):
        Model(Params(**FAST, **WAVE, em_t_on_s=0.3e-6))
    with pytest.raises(ValueError, match="empty"):
        Model(Params(**FAST, **WAVE, em_t_on_s=0.4e-6, em_t_off_s=0.4e-6))


def test_a_window_past_t_end_is_no_wave():
    off = run(Params(**FAST))
    late = run(Params(**FAST, **WAVE, em_t_on_s=0.7e-6, em_t_off_s=1e-6))
    assert late.ledger["E_em"][-1] == 0.0
    for s in off.n:
        assert np.array_equal(off.n[s], late.n[s]), s


def test_the_afterglow_wave_books_its_energy_in_its_window_alone():
    off = run(Params(**FAST))
    on = run(Params(**FAST, **WAVE, em_t_on_s=0.3e-6, em_t_off_s=0.45e-6))
    _check_closures(on)
    em = on.ledger["E_em"]
    before = on.t <= 0.3e-6
    assert np.all(em[before] == 0.0) and em[-1] > 0.0
    assert np.ptp(em[on.t >= 0.45e-6]) == 0.0
    # the same pump, and the same run until the wave opens
    assert on.ledger["E_dep"][-1] == approx(off.ledger["E_dep"][-1],
                                            rel=1e-12)
    assert np.array_equal(on.t[before], off.t[off.t <= 0.3e-6])
    for s in off.n:
        assert np.array_equal(on.n[s][before], off.n[s][off.t <= 0.3e-6]), s
    # the wave heats the afterglow's electrons and they live longer
    w = {"w": (0.3e-6, 0.45e-6)}
    a, b = reheat.measure(off, w), reheat.measure(on, w)
    assert b["windows"]["w"]["te_e"] > a["windows"]["w"]["te_e"] + 0.2
    assert b["windows"]["w"]["ne_int"] > a["windows"]["w"]["ne_int"]
    assert b["windows"]["w"]["wave_in"] == approx(b["wave"], rel=1e-12)
    text = summary(on)
    assert "from 0.3 to 0.45 us (the afterglow)" in text
    d = row_dict(on)
    assert (d["em_on_us"], d["em_off_us"]) == approx((0.3, 0.45))
    assert d["E248_J_cm3"] == approx(d["emission_248"] / 100
                                     * d["E_dep_J_cm3"], rel=1e-12)
    assert d["given_J_cm3"] + d["em_J_cm3"] == approx(d["E_dep_J_cm3"])


def test_a_subclass_written_for_v0_5_still_runs():
    """rhs and _kinetic_update are called as v0_5 called them, (t, y, on),
    wherever no window parts the wave from the given power, so a subclass
    or wrapper with v0_5's three arguments (the beam drive of
    docs/survey_2026-09-27/scripts/ebeam_screen.py, test_kinetic.py's
    counters) runs every such run; only a window needs the fourth."""
    class Old(Model):
        def rhs(self, t, y, on):
            return super().rhs(t, y, on)
    for kw in (dict(), WAVE):
        a = Old(Params(**FAST, **kw)).run()
        b = run(Params(**FAST, **kw))
        assert np.array_equal(a.n["KrF(B)"], b.n["KrF(B)"])
        assert np.array_equal(a.ledger["E_em"], b.ledger["E_em"])
    with pytest.raises(TypeError):
        Old(Params(**FAST, **WAVE, **AFTER)).run()


def test_a_window_inside_the_pulse_closes():
    on = run(Params(**FAST, **WAVE, em_t_on_s=0.1e-6, em_t_off_s=0.2e-6))
    _check_closures(on)
    em = on.ledger["E_em"]
    assert np.all(em[on.t <= 0.1e-6] == 0.0)
    assert np.ptp(em[on.t >= 0.2e-6]) == 0.0 and em[-1] > 0.0


def test_the_kinetic_eedf_sees_the_wave_in_its_window_alone():
    """A short kinetic run with the wave alone for 10 ns after a 20-ns
    pulse (about three minutes, most of it the pulse's first nanosecond):
    the ledger closes, the pulse's solutions carry no wave, those after it
    carry the solver's own EM power, and the summary counts them."""
    res = run(Params(eedf="kinetic", pulse_s=20e-9, t_end_s=30e-9,
                     n_save=121, **WAVE, em_t_on_s=20e-9,
                     em_t_off_s=math.inf))
    _check_closures(res, tol=1e-8)
    log = [e for e in res.kinetic if not e.get("skipped")]
    pulse = [e for e in log if e["on"]]
    alone = [e for e in log if e["wave"] and not e["on"]]
    assert pulse and alone
    assert all(not e["wave"] and e["pem"] == 0.0 for e in pulse)
    assert all(e["pem"] > 0.0 for e in alone)
    # the pulse's balance is v0_5's (the field carries the given power)
    late = [e for e in pulse if e["t"] >= 10e-9
            and np.isfinite(e["power_ratio"])]
    assert late and 0.98 < min(e["power_ratio"] for e in late)
    assert max(e["power_ratio"] for e in late) < 1.03
    assert all(np.isfinite(e["dc_share"]) for e in alone)
    assert accounting(res)["e_em"] > 0.0
    assert "with the wave alone" in summary(res)


# ---------------------------------------------------------------------------
# the twin pairs
# ---------------------------------------------------------------------------
def _rec(given, wave, e248, vuv=0.0):
    w = {"all": dict(g0=1.0, alpha=0.5, net=0.5, ne_open=1e14, ne_int=1e7,
                     te_e=1.0, wave_in=wave)}
    return dict(given=given, wave=wave, E248=e248, Evuv=vuv, Eother=0.0,
                heat=given + wave - e248 - vuv, stored=0.0, walls=0.0,
                tg_end=700.0, windows=w, closure=0.0, warnings=[],
                kinetic=None)


def test_twin_arithmetic():
    off = _rec(0.3, 0.0, 0.018)                      # the pump: 6 %
    t = reheat.twin(off, _rec(0.3, 0.1, 0.028), "all")
    assert t["marginal"] == approx(0.10)
    assert t["eff_off"] == approx(0.06)
    assert t["eff_on"] == approx(0.028 / 0.4)
    assert t["d_eff"] == approx(0.07 - 0.06)
    assert t["source"] == approx(0.6)
    assert sum(t["to"].values()) == approx(1.0)
    # the efficiency rises exactly where the marginal beats the pump's
    for e248, sign in ((0.018 + 0.0061, 1), (0.018 + 0.0059, -1)):
        assert np.sign(reheat.twin(off, _rec(0.3, 0.1, e248),
                                   "all")["d_eff"]) == sign
    # a wave that adds no 248 nm never breaks even, and none absorbed
    assert reheat.twin(off, _rec(0.3, 0.1, 0.017), "all")["source"] \
        == math.inf
    assert math.isnan(reheat.twin(off, _rec(0.3, 0.0, 0.018),
                                  "all")["marginal"])
    with pytest.raises(ValueError):
        reheat.twin(off, _rec(0.31, 0.1, 0.02), "all")


def test_unc_stats_is_the_overviews():
    s = reheat.unc_stats(1.0, {"a": (0.9, 1.3), "b": (1.1, 0.8)})
    assert s["h"] == approx({"a": 0.2, "b": -0.15})
    assert s["quad"] == approx(0.25)
    assert (s["lo"], s["hi"]) == approx((0.7, 1.4))


def test_specs_and_their_params():
    p = reheat.params_of(("two_t", "env", "1", "early", "band:kr2_f=1e-09"))
    assert (p.eedf, p.pressure_atm, p.x_f2) == ("two_t", 0.5, 0.02)
    assert p.rates == (("kr2_f", 1e-9),)
    assert p.em_E_V_cm == approx(47.6e3 / math.sqrt(2.0))
    assert p.em_f_Hz == 6.6e12
    assert (p.em_t_on_s, p.em_t_off_s) == approx((1e-6, 1.1e-6))
    p = reheat.params_of(("kinetic", "ref", "off", "", "row:R11p=0.9"))
    assert p.row_scale == (("R11p", 0.9),) and p.em_E_V_cm == 0.0
    p = reheat.params_of(("two_t", "ref", "ctl", "all", ""))
    assert p.em_E_V_cm == reheat.CONTROL_V_CM and p.em_t_off_s == math.inf
    p = reheat.params_of(("two_t", "ref", "1", "all",
                          "band:trimer_tg_pow=-0.5"))
    assert p.trimer_tg_pow == -0.5 and p.rates == ()
    every = reheat.specs()
    assert len({reheat.key(s) for s in every}) == len(every)
    for s in every:                     # every one is a valid model
        Model(reheat.params_of(s))
    # the kinetic bands and the further pairs: valid, unique, each pair
    # with its off twin
    more = reheat.specs(kinetic_bands=("ref",), kinetic_inputs=("kr2_f",),
                        further=True)
    assert len(set(more)) == len(more)
    assert set(reheat.further_specs()) <= set(more)
    for s in more:
        Model(reheat.params_of(s))
    fs = reheat.further_specs()
    assert len(set(fs)) == len(fs)
    assert all((c, p, "off", "", i) in fs
               for c, p, _, _, i in reheat.FURTHER)
    n_inputs = len(reheat.inputs())
    # an end at its default is the base run (kr2_f, kr2_pen, kr2f_f,
    # kr2_form, kr2_f2, kr_pen's hi, kr2f_f2, trimer_tg_pow, krf_f, kps_kr's
    # hi, kstep_s, kstep_p, k5p_f2's hi, ii_yield's hi, ii_f_c)
    ends = sum((lo != d) + (hi != d) for _, _, lo, hi, d in reheat.inputs())
    per_end = 1 + len(reheat.BAND_GRID)          # the off twin and the pairs
    assert len([s for s in every if s[4]]) == 2 * per_end * ends
    assert ends < 2 * n_inputs


def test_the_search_for_the_optimal_reheat():
    """The search's specs (2026-10-20) are valid and unique, each pair
    with its off twin; its windows and points resolve; and the optimum is
    the searched pair with the largest change in efficiency, counting a
    pair whose off twin was saved before the search's windows existed."""
    o = reheat.opt_specs()
    assert len(set(o)) == len(o)
    for s in o:
        Model(reheat.params_of(s))
    assert all(("kinetic", pt, "off", "", inp) in o
               for _, pt, wave, _, inp in o if wave != "off")
    p = reheat.params_of(("kinetic", "ref", "1", "d100", ""))
    assert (p.em_t_on_s, p.em_t_off_s) == approx((1.1e-6, math.inf))
    p = reheat.params_of(("kinetic", "ref", "1", "w30", ""))
    assert p.em_t_off_s == approx(1.03e-6)
    assert reheat.params_of(("kinetic", "ref_f2_0.25", "off", "", "")
                            ).x_f2 == 0.0025
    off = _rec(0.3, 0.0, 0.018)
    cache = {reheat.key(("kinetic", "ref", "off", "", "")): off}
    for wave, w, e248 in (("1", "all", 0.0185), ("1.778", "all", 0.0195),
                          ("1", "d30", 0.0190)):
        r = _rec(0.3, 0.01, e248)
        r["windows"] = {w: r["windows"]["all"]}
        cache[reheat.key(("kinetic", "ref", wave, w, ""))] = r
    best = reheat.optimum(cache, "")
    assert best[:4] == ("field", "ref", "1.778", "all")
    assert best[4]["d_eff"] > 0.0
    t = reheat.pair(cache, "kinetic", "ref", "1", "d30")
    assert math.isnan(t["d_net"]) and np.isfinite(t["d_eff"])
    assert "Optimum (n_e below" in reheat.tables(cache)
    # a wave that drives n_e past the notes' design density is left out
    k = reheat.key(("kinetic", "ref", "1.778", "all", ""))
    cache[k]["windows"]["all"] = dict(cache[k]["windows"]["all"],
                                      ne_max=2.0 * reheat.OPT_NE_LIMIT)
    assert reheat.optimum(cache, "")[:4] == ("window", "ref", "1", "d30")
    assert "(past the limit)" in reheat.tables(cache)


def test_the_reference_pair_under_two_t_is_pinned():
    """The study's two_t headline at the reference point (results/
    reheat.md, 2026-10-19): the notes' wave over the whole afterglow
    sustains the plasma, absorbs three quarters of the pump's energy and
    adds 0.44 % of it as 248 nm, so the efficiency falls 2.25 points; the
    gentle wave (1/100 of the power) takes 248 nm away."""
    rec = {w: reheat.measure(run(reheat.params_of(("two_t", "ref", w,
                                                   "" if w == "off" else
                                                   "all", ""))))
           for w in ("off", "1", "0.1")}
    t = reheat.twin(rec["off"], rec["1"], "all")
    assert t["wave"] == approx(0.22566, rel=2e-3)
    assert t["marginal"] == approx(0.004406, rel=2e-3)
    assert t["eff_off"] == approx(0.056835, rel=2e-4)
    assert t["d_eff"] == approx(-0.022507, rel=2e-3)
    g = reheat.twin(rec["off"], rec["0.1"], "all")
    assert g["wave"] == approx(1.2389e-05, rel=2e-3)
    assert g["marginal"] == approx(-0.014194, rel=1e-2)


def test_cli_twin_and_report(capsys, tmp_path):
    from v0_6.run import main
    assert main(["--pulse", "0.3", "--t-end", "0.6", "--eedf", "two_t",
                 "--em-E", "3e4", "--em-f", "6200", "--em-on", "0.3",
                 "--em-off", "inf", "--twin"]) == 0
    out = capsys.readouterr().out
    assert "TWIN PAIR" in out and "marginal" in out
    assert "from 0.3 to 0.6 us (the afterglow)" in out
    runs, md = tmp_path / "runs.json", tmp_path / "reheat.md"
    assert main(["reheat", "--report", "--runs", str(runs),
                 "--out", str(md)]) == 0
    assert md.read_text().startswith("# LIGHT v0_6: the afterglow reheat")
