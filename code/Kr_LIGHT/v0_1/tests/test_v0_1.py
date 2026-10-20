"""LIGHT_Kr v0_1 (2026-10-13): the pure-Kr model is LIGHT v0_5 with the
fluorine taken out, and nothing else changed but two corrections of
jung2006's tables to the paper. Held here:

* the species, rows and bands are v0_5's krypton ones, field for field and
  in v0_5's order, and no fluorine is left anywhere;
* jung2006 differs from v0_5's in the two corrected entries alone;
* from one state, v0_1 and v0_5 (no F2 in its fill, the same tables) build
  the same kinetic EEDF table and right-hand side to round-off;
* over whole runs they agree to the solver's tolerance (the state vectors
  differ in size, so the stiff solver's steps differ), v0_5's fluorine
  densities and events stay exactly 0, and the difference shrinks with
  rtol. Under the kinetic closure the runs differ by the closure's own
  spread between two tolerances (3e-4 on the 147-nm share on the short
  pulse, 2026-10-13); that run is slow (LIGHT_KINETIC_SLOW=1).

Every comparison with v0_5 skips when v0_5 is not beside v0_1."""

import copy
import importlib
import os
from contextlib import contextmanager
from dataclasses import astuple, fields

import numpy as np
import pytest

from v0_1 import boltzmann as KB1
from v0_1 import ib as IB1
from v0_1 import jung2006 as J1
from v0_1 import model as M1
from v0_1 import network as N1
from v0_1 import report as R1
from v0_1 import species as S1
from v0_1 import xsec as X1

try:
    M5 = importlib.import_module("v0_5.model")
    N5 = importlib.import_module("v0_5.network")
    S5 = importlib.import_module("v0_5.species")
    J5 = importlib.import_module("v0_5.jung2006")
except ImportError:                      # v0_5 not beside v0_1
    M5 = N5 = S5 = J5 = None

needs_v0_5 = pytest.mark.skipif(M5 is None, reason="LIGHT v0_5 is not "
                                "beside v0_1 (LIGHT_Kr/v0_5)")
SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_KINETIC_SLOW"),
                          reason="minutes: set LIGHT_KINETIC_SLOW=1")
KR = ("Kr", "Kr*", "Kr**", "Kr+", "Kr2+", "Kr2*")


def _fluorine(row):
    return any("F" in s for s in row.reactants + row.products)


def _clear(J):
    for f in (J.process, J.rate_table, J.superelastic_table,
              J.capture_table):
        f.cache_clear()


@contextmanager
def v0_5_on_v0_1s_tables():
    """v0_5 with jung2006's tables as v0_1 carries them (the two
    corrections), restored afterwards. TABLES holds the two dicts, so they
    change in place."""
    saved = copy.deepcopy(J5.TABLE_II), copy.deepcopy(J5.TABLE_III)
    try:
        for mine, theirs in ((J5.TABLE_II, J1.TABLE_II),
                             (J5.TABLE_III, J1.TABLE_III)):
            mine.clear()
            mine.update(copy.deepcopy(theirs))
        _clear(J5)
        yield
    finally:
        for mine, old in zip((J5.TABLE_II, J5.TABLE_III), saved):
            mine.clear()
            mine.update(old)
        _clear(J5)


# ---------------------------------------------------------------------------
# the model as data
# ---------------------------------------------------------------------------
@needs_v0_5
@pytest.mark.parametrize("kw", [{}, {"groups": ("all",)}, {"k76": "kannari"},
                                {"groups": ("electron",), "k76": "kannari"}])
def test_rows_are_v0_5s_krypton_rows(kw):
    r1 = N1.build(**kw)
    r5 = [r for r in N5.build(**kw) if not _fluorine(r)]
    assert [r.rid for r in r1] == [r.rid for r in r5]
    for a, b in zip(r1, r5):
        diff = [f.name for f in fields(a) if getattr(a, f.name)
                != getattr(b, f.name)]
        assert not diff, (a.rid, diff)


@needs_v0_5
def test_species_bands_and_constants_are_v0_5s():
    assert list(S1.SPECIES) == [s for s in S5.SPECIES if "F" not in s]
    for name, s in S1.SPECIES.items():
        assert astuple(s) == astuple(S5.SPECIES[name]), name
    read5 = N5.bands_read([r for r in N5.build(("all",))
                           if not _fluorine(r)])
    assert set(N1.BANDS) == read5 | set(N1.FREE_BANDS)
    for name, b in N1.BANDS.items():
        assert astuple(b) == astuple(N5.BANDS[name]), name
    assert N1.ION_K0 == {k: v for k, v in N5.ION_K0.items() if "F" not in k}
    assert N1.LADDER_KANNARI == {k: v for k, v in N5.LADDER_KANNARI.items()
                                 if k != "k5p_f2"}
    for k in ("KR_5P_TAU_S", "G_5S", "G_5P", "KANNARI_K76", "KANNARI_KSTEP_S",
              "KANNARI_KSTEP_P", "KR_CONV_300"):
        assert getattr(N1, k) == getattr(N5, k), k


def test_no_fluorine_is_left():
    """No species, row, band, cross-section file, gas, parameter or report
    line of fluorine."""
    assert all("F" not in s for s in S1.SPECIES)
    assert not any(_fluorine(r) for r in N1.build(("all",)))
    assert list(X1.FILES) == ["Kr"] and KB1.GASES == ("Kr",)
    assert list(IB1.MT_BLOCK) == ["Kr"]
    names = {f.name for f in fields(M1.Params)}
    assert not names & {"x_f2", "ne_mode", "tau_sp_s", "trimer_tg_pow",
                        "network"}
    text = R1.summary(M1.run(M1.Params(eedf="two_t", pulse_s=0.1e-6,
                                       t_end_s=0.2e-6, n_save=101)))
    for word in ("KrF", "F2", "248", "Kr2F", "fluorine"):
        assert word not in text, word


@needs_v0_5
def test_jung2006_differs_from_v0_5s_by_the_two_corrections_alone():
    """Table II's 1s5 -> 2p6 at 4 eV: 23 as printed (v0_5: 21). Table
    III: no 1s3 -> 2p5 (v0_5 booked Table IV's 1s5 -> 2p5 there)."""
    t2 = copy.deepcopy(J5.TABLE_II)
    assert t2[4]["2p6"] == 21.0 and J1.TABLE_II[4]["2p6"] == 23.0
    t2[4]["2p6"] = 23.0
    assert t2 == J1.TABLE_II
    t3 = copy.deepcopy(J5.TABLE_III)
    assert t3.pop(3.5) == {"2p5": 0.4} and t3[6].pop("2p5") == 0.11
    assert t3 == J1.TABLE_III
    for k in ("LEVELS", "WEIGHTS_5S", "PROXY", "UNCERTAINTY", "SHAPE_REF",
              "E_MAX_TABLE", "CM2"):
        assert getattr(J1, k) == getattr(J5, k), k


# ---------------------------------------------------------------------------
# the kinetic machinery from one state
# ---------------------------------------------------------------------------
@needs_v0_5
@pytest.mark.parametrize("memory", [0.0, np.inf])
def test_kinetic_table_and_rhs_from_one_state_are_v0_5s(memory):
    """A mid-pulse state of the reference point handed to both: the EEDF
    table (every krypton row's rates, the captured energies, the elastic
    and wave powers) is v0_5's bit for bit, and the right-hand side to
    round-off. v0_5's table carries its idle fluorine rows besides."""
    ref = M1.run(M1.Params(eedf="two_t"))
    k = int(np.searchsorted(ref.t, 0.5e-6))
    with v0_5_on_v0_1s_tables():
        m1 = M1.Model(M1.Params(eedf="kinetic", kin_memory_s=memory))
        m5 = M5.Model(M5.Params(x_f2=0.0, eedf="kinetic",
                                kin_memory_s=memory))
        ys = []
        for m in (m1, m5):
            y = np.zeros(m.size)
            for s in KR:
                y[m.idx[s]] = ref.n[s][k] / m.N0
            y[m.i_ue] = ref.ue[k] / m.N0
            y[m.i_tg] = ref.tg[k]
            ys.append(y)
            m._kinetic_update(ref.t[k], y, True)
            if m._tab is None:                  # a time step's start
                m._kinetic_update(ref.t[k] + 1e-12, y, True)
        t1, t5 = m1._tab, m5._tab
        assert np.array_equal(t1["lnte"], t5["lnte"])
        by5 = {m5.rows[j].rid: j for j in t5["lnk"]}
        for j, lnk in t1["lnk"].items():
            assert np.array_equal(lnk, t5["lnk"][by5[m1.rows[j].rid]])
        for j, cap in t1["cap"].items():
            assert np.array_equal(cap, t5["cap"][by5[m1.rows[j].rid]])
        assert np.array_equal(t1["pel"], t5["pel"])
        assert np.array_equal(t1["pem"], t5["pem"])
        assert all(_fluorine(m5.rows[j]) for j in t5["lnk"]
                   if m5.rows[j].rid not in {r.rid for r in m1.rows})
        d1 = m1.rhs(ref.t[k], ys[0], True)
        d5 = m5.rhs(ref.t[k], ys[1], True)
    for s in KR:
        assert d1[m1.idx[s]] == pytest.approx(d5[m5.idx[s]], rel=1e-14,
                                              abs=0.0), s
    for i1, i5 in ((m1.i_ue, m5.i_ue), (m1.i_tg, m5.i_tg)):
        assert d1[i1] == pytest.approx(d5[i5], rel=1e-14, abs=0.0)


# ---------------------------------------------------------------------------
# whole runs
# ---------------------------------------------------------------------------
def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))
                 / max(float(np.max(np.abs(b))), 1e-300))


def _parity(kw, rtol):
    """Max relative difference of the krypton densities, U_e, Tg, the
    ledger and the krypton events, and of v0_1's row_dict on both runs;
    v0_5's fluorine stays exactly 0."""
    r1 = M1.Model(M1.Params(rtol=rtol, **kw)).run()
    with v0_5_on_v0_1s_tables():
        r5 = M5.Model(M5.Params(x_f2=0.0, rtol=rtol, **kw)).run()
    assert np.array_equal(r1.t, r5.t)
    worst = max([_rel(r1.n[s], r5.n[s]) for s in KR]
                + [_rel(r1.ue, r5.ue), _rel(r1.tg, r5.tg)]
                + [_rel(r1.ledger[k], r5.ledger[k]) for k in r1.ledger
                   if np.abs(r5.ledger[k]).max() > 0]
                + [_rel(r1.events[r], r5.events[r]) for r in r1.events
                   if np.abs(r5.events[r]).max() > 0])
    for k in r1.ledger:
        if not np.abs(r5.ledger[k]).max() > 0:
            assert not np.abs(r1.ledger[k]).max() > 0, k
    for s in r5.n:
        if s not in KR:
            assert not np.any(r5.n[s]), s
    for r in r5.events:
        if r not in r1.events:
            assert not np.any(r5.events[r]), r
    d1, d5 = R1.row_dict(r1), R1.row_dict(r5)
    rows = max(abs(d1[k] - d5[k]) / abs(d5[k]) for k in d1
               if isinstance(d1[k], float) and k != "closure"
               and np.isfinite(d5[k]) and d5[k] != 0.0)
    return worst, rows


@needs_v0_5
@pytest.mark.parametrize("kw", [
    dict(eedf="two_t"),
    dict(eedf="two_t", groups=("all",)),
    dict(eedf="two_t", k76="kannari"),
    dict(eedf="two_t", ladder="kannari"),
    dict(eedf="two_t", gap_cm=0.05, pressure_atm=0.1),
    dict(eedf="two_t", em_E_V_cm=1500.0, em_f_Hz=95e9, pulse_s=0.3e-6,
         t_end_s=0.6e-6),
    dict(eedf="two_t", pressure_atm=2.0, power_W_cm3=1e5, ne_cm3=1e16),
    dict(eedf="two_t", rates=(("kr2_pen", 5e-10), ("kr2_form", 5.5e-32),
                              ("kstep_s", 1.5e-7), ("kstep_p", 2e-7),
                              ("kps_kr", 3e-11))),
    dict(eedf="maxwellian", capture="three_halves"),
])
def test_runs_are_v0_5s_without_f2(kw):
    """Measured at the default rtol (2026-10-13): 3e-7 at most on the
    trajectories, 6e-7 on the row_dict numbers."""
    worst, rows = _parity(kw, 1e-8)
    assert worst < 1e-5 and rows < 1e-5


@needs_v0_5
def test_the_difference_is_the_solvers_path():
    """At rtol 1e-10 the reference point's difference falls a hundredfold
    (3e-9 measured), as a difference of solver steps does and a difference
    of physics would not."""
    kw = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
    coarse, _ = _parity(kw, 1e-8)
    fine, _ = _parity(kw, 1e-10)
    assert fine < 1e-7 and fine < 0.1 * coarse


@SLOW
@needs_v0_5
def test_kinetic_runs_agree_within_the_closures_own_spread():
    """Under the kinetic EEDF the table is renewed where the drift triggers
    fire, which the solver's steps place. v0_1 against v0_5 differs by what
    v0_1 at rtol 1e-8 differs from itself at rtol 1e-9 (the same 1.2e-3 on
    the elastic share and 3e-4 on the 147-nm share, 2026-10-13)."""
    kw = dict(eedf="kinetic", pulse_s=0.3e-6, t_end_s=0.6e-6)
    a = R1.row_dict(M1.Model(M1.Params(rtol=1e-8, **kw)).run())
    b = R1.row_dict(M1.Model(M1.Params(rtol=1e-9, **kw)).run())
    with v0_5_on_v0_1s_tables():
        c = R1.row_dict(M5.Model(M5.Params(x_f2=0.0, rtol=1e-8, **kw)).run())

    def spread(x, y):
        return max(abs(x[k] - y[k]) / abs(y[k]) for k in x
                   if isinstance(x[k], float) and k != "closure"
                   and np.isfinite(y[k]) and y[k] != 0.0)
    own, cross = spread(a, b), spread(a, c)
    assert cross < 3.0 * own + 1e-6
    assert abs(a["emission_vuv"] / c["emission_vuv"] - 1.0) < 2e-3
