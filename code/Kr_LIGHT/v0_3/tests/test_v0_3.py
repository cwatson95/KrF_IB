"""LIGHT_Kr v0_3 (2026-10-15): v0_2 with its electron cross sections
trained on krypton's swarm data (ML/ml/kr_swarm.py). The gate:

* everything but the cross sections' factors is v0_2's: rows, bands,
  species, jung2006, the tail table and Biagi's file, byte for byte;
* with unit factors (xsec.UNIT) v0_3 is v0_2: the cross sections, the
  swarm and a model run number for number;
* the trained factors are the pinned file's, inside the box the training
  searched, and their record names the train set;
* the swarm's transverse diffusion obeys the Einstein relation in the
  thermal limit;
* v0_2's Hornbeck-Molnar and associative-ionization checks, and the swarm
  against LoKI-MC's Monte Carlo on Biagi's own file (UNIT);
* the kinetic solver's one change, a failure path v0_2's runs never took:
  where no EEDF is found from a warm start it starts again cold, once,
  and counts it (Solver.restarts), with the cold start's answer.

Comparisons with v0_2 skip when v0_2 is not beside v0_3."""

import importlib
import json
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

from v0_3 import boltzmann as KB
from v0_3 import network, swarm, xsec
from v0_3.model import Model, Params
from v0_3.report import accounting, row_dict
from v0_3.species import SPECIES, e_form

HERE = Path(__file__).resolve()
V0_2 = HERE.parents[2] / "v0_2"
needs_v0_2 = pytest.mark.skipif(not (V0_2 / "network.py").is_file(),
                                reason="v0_2 is not beside v0_3")
FAST = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
KT300 = 1.380649e-23 * 300.0 / 1.602176634e-19


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


# ------------------------------------------------------------- vs v0_2
@needs_v0_2
def test_everything_but_the_factors_is_v0_2s():
    old_net = importlib.import_module("v0_2.network")
    old_sp = importlib.import_module("v0_2.species")
    old = old_net.build(("all",))
    new = network.build(("all",))
    assert [r.rid for r in new] == [r.rid for r in old]
    for a, b in zip(old, new):
        assert astuple(a) == astuple(b), a.rid
    assert {k: astuple(v) for k, v in network.BANDS.items()} == {
        k: astuple(v) for k, v in old_net.BANDS.items()}
    assert {k: astuple(v) for k, v in SPECIES.items()} == {
        k: astuple(v) for k, v in old_sp.SPECIES.items()}
    j2 = importlib.import_module("v0_2.jung2006")
    j3 = importlib.import_module("v0_3.jung2006")
    for name in ("LEVELS", "PROXY", "SHAPE_REF", "TABLES", "TABLE_II",
                 "TABLE_III", "UNCERTAINTY", "WEIGHTS_5S", "E_MAX_TABLE",
                 "CM2"):
        assert getattr(j3, name) == getattr(j2, name), name
    for f in ("data/kr_tail_ratio.txt", "data/Kr_Biagi897_LXCat.txt"):
        assert (HERE.parents[1] / f).read_bytes() == (V0_2 / f).read_bytes()


@needs_v0_2
def test_unit_factors_are_v0_2():
    """xsec.UNIT: the file's own arrays, the swarm and a model run number
    for number v0_2's."""
    x2 = importlib.import_module("v0_2.xsec")
    s2 = importlib.import_module("v0_2.swarm")
    m2 = importlib.import_module("v0_2.model")
    r2 = importlib.import_module("v0_2.report")
    with xsec.using(xsec.UNIT):
        a, b = xsec.processes("Kr"), x2.processes("Kr")
        assert len(a) == len(b)
        for p, q in zip(a, b):
            assert (p.product, p.kind, p.threshold) == (q.product, q.kind,
                                                       q.threshold)
            assert np.array_equal(p.sigma, q.sigma)
        for en in (0.01, 1.0, 30.0, 200.0):
            d, d2 = swarm.point(en), s2.point(en)
            assert (d["W"], d["alphaN"], d["mean"]) == (
                d2["W"], d2["alphaN"], d2["mean"])
        res = Model(Params(**FAST)).run()
    res2 = m2.Model(m2.Params(**FAST)).run()
    assert row_dict(res) == r2.row_dict(res2)


# ------------------------------------------------------------- the factors
def test_trained_factors_are_the_pinned_files():
    name, md5 = xsec.TRAINED_FILE
    assert md5 is not None, "v0_3 runs on trained factors"
    path = xsec.DATA_DIR / name
    assert xsec.md5_of(path) == md5
    f = xsec.trained()
    assert f == xsec.read_factors(path) == xsec.factors()
    rec = json.loads(path.read_text())
    assert rec["train"] and all(n.startswith(("w_", "a_", "dt_"))
                                for n in rec["train"])
    # inside the box the emulator was trained on (ML/ml/kr_swarm.BOX = 0.6)
    z = np.log(np.concatenate([f.m, f.exc]))
    assert np.all(np.abs(z) < 0.6)
    assert f.ion == 1.0
    assert len(f.m) == len(f.m_knots) == 6


def test_factors_scale_the_named_groups_only():
    f = xsec.Factors((0.1, 10.0), (1.0, 2.0), (0.5, 0.7, 0.9), 1.0)
    with xsec.using(f):
        new = {p.product: p for p in xsec.processes("Kr")}
    raw = {p.product: p for p in xsec.raw_processes("Kr")}
    for prod, p in new.items():
        q = raw[prod]
        if p.kind == "excitation":
            x = {"1s": 0.5, "2p": 0.7, "hi": 0.9}[xsec.exc_group(prod)]
            assert np.allclose(p.sigma, x * q.sigma, rtol=1e-15, atol=0)
        elif p.kind == "elastic":
            r = f.m_at(p.energy)
            assert np.allclose(p.sigma, r * q.sigma, rtol=1e-15, atol=0)
            assert r[0] == 1.0 and r[-1] == 2.0
        else:
            assert np.array_equal(p.sigma, q.sigma)
    assert sum(xsec.exc_group(p) == "1s" for p in raw
               if raw[p].kind == "excitation") == 4
    assert sum(xsec.exc_group(p) == "2p" for p in raw
               if raw[p].kind == "excitation") == 10


# ------------------------------------------------------------- the swarm
def test_transverse_diffusion_is_einsteins_in_the_thermal_limit():
    """D_T/mu = kTg/e when the field barely heats the electrons (Einstein's
    relation; 2.6e-5 off at 1e-4 Td, 2026-10-14)."""
    for tg in (195.0, 300.0):
        d = swarm.point(1e-4, tg)
        assert d["DTmu"] == approx(KT300 * tg / 300.0, rel=2e-4)
        assert d["mean"] == approx(1.5 * KT300 * tg / 300.0, rel=1e-3)


# LoKI-MC (the vendored Monte Carlo solver) on Biagi's 8.97 file, 300 K,
# 1 Torr, 1e5 electrons, equal sharing (2026-10-14, v0_2's run): the SST
# reduced Townsend coefficient [m^2], the mean energy [eV] and the
# ionizations counted (a statistical error of 19, 4.3, 1.8 and 1.0 %)
LOKI_MC_897 = {20.0: (3.6057e-25, 4.5844, 29), 30.0: (7.1147e-24, 4.7536, 543),
               40.0: (3.3102e-23, 4.9022, 3127),
               50.0: (8.7447e-23, 5.0360, 9521)}


@pytest.mark.parametrize("en", sorted(LOKI_MC_897))
def test_swarm_against_monte_carlo_on_biagis_file(en):
    with xsec.using(xsec.UNIT):
        d = swarm.point(en)
    a_mc, m_mc, counted = LOKI_MC_897[en]
    assert d["alphaN"] * 1e-4 == approx(a_mc, rel=0.02 + 2.0 / counted ** 0.5)
    assert d["mean"] == approx(m_mc, rel=2e-3)


# --------------------------------------------- v0_2's Hornbeck-Molnar checks
def _state(m, te=1.8):
    y = m.initial_state()
    y[m.idx["Kr*"]] = 1e-5
    y[m.idx["Kr2+"]] = 2e-6
    y[m.idx["Kr+"]] = 1e-6
    ne = m.electron_density(np.maximum(y[:m.ns], 0.0) * m.N0)
    y[m.i_ue] = 1.5 * ne * te / m.N0
    return y


def test_hm_rows_at_zero_change_nothing():
    m = Model(Params(eedf="two_t"))
    bare = Model(Params(eedf="two_t"))
    bare.compiled = [c for c in bare.compiled
                     if not c.row.rid.startswith("HM_")]
    y = _state(m)
    d, d0 = m.rhs(0.0, y, True), bare.rhs(0.0, y, True)
    live = list(range(m.ns)) + [m.i_ue, m.i_tg]
    assert np.array_equal(d[live], d0[live])
    k = Model(Params())
    rows = [pr.rows[0] for pr in k.kin.procs]
    assert len(rows) == 52 and not any(r.startswith("HM_") for r in rows)
    sig = {pr.rows[0]: pr.p.sigma for pr in k.kin.procs}
    for c in k.compiled:
        if c.row.rid.startswith("C1c_"):
            p = xsec.find(*c.row.rate[1:4])
            assert np.array_equal(sig[c.row.rid], p.sigma)


def test_hm_at_one_ionizes_every_such_excitation():
    res0 = Model(Params(**FAST)).run()
    res1 = Model(Params(**FAST, rates=(("hm_kr", 1.0),))).run()
    hm = [r for r in res1.events if r.startswith("HM_")]
    assert len(hm) == 20
    made = sum(res1.events[r][-1] for r in hm)
    assert made > 0.0
    for r in hm:
        assert res1.events["C1c_" + r[3:]][-1] == 0.0
        assert res0.events[r][-1] == 0.0
    lev0 = sum(res0.events["C1c_" + r[3:]][-1] for r in hm)
    assert made == approx(lev0, rel=0.05)
    a1 = accounting(res1)
    assert abs(a1["closure"]) < 1e-12
    assert a1["flows"]["   of which Hornbeck-Molnar (HM)"] == approx(
        made * e_form("Kr2+"), rel=1e-12)
    d0, d1 = row_dict(res0), row_dict(res1)
    assert d1["to_Krplus"] > d0["to_Krplus"]
    assert d1["emission_vuv"] == approx(d0["emission_vuv"], rel=2e-3)


def test_associative_ionization_in_the_boltzmann_equation():
    """v0_2's check on v0_3's cross sections: the associative electrons
    count in the growth term and the energy balance closes, less what they
    bring, born in the first cell's centre (du/2 = 0.025 eV here: 7.3e-5 of
    the field's power on the trained set, 8.6e-5 on Biagi's file as
    published)."""
    procs = swarm.processes(hm=1.0)
    assert sum(q.kind == "assoc" for q in procs) == 20
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    st = KB.State(x={"Kr": 1.0}, tg=450.0, n_col=2.45e25, ne=7.4e19)
    sol = s.solve(st, 60.0)
    g, F, N = sol.grid, sol.F, st.n_col
    k_dir = sum(g.rate(q.p, F) for q in procs if q.kind == "ionization")
    k_as = sum(g.rate(q.p, F) for q in procs if q.kind == "assoc")
    assert k_as > 0.5 * k_dir
    assert sol.nu == approx((k_dir + k_as) * N, rel=1e-6)
    inel = N * sum(g.rate(q.p, F) * q.threshold for q in procs)
    el = N * s.elastic_power(st, sol, "Kr")
    gain = s.field_power(st, sol)
    born = k_as * N * g.ec[0]
    loss = inel + el + sol.nu * sol.mean - born
    assert abs(gain - loss) / gain < 2e-5
    assert born / gain == approx(7.32e-5, rel=0.05)


def test_a_failed_warm_start_restarts_cold(monkeypatch):
    """Solver.eedf: where a grid's solve finds no solution from the warm
    start, it starts again cold, once, and counts the restart; the mean
    match has one solution from any start, so the answer is the cold
    start's. A cold start that fails still raises."""
    procs = swarm.processes()
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    st = KB.State(x={"Kr": 1.0}, tg=450.0, n_col=2.45e25, ne=7.4e19)
    cold = s.eedf(st, 3.0)
    assert s.restarts == 0
    real = KB.Solver._eedf_on_grid

    def no_warm_start(self, state, target_mean, warm, growth):
        if warm is not None:
            raise RuntimeError("kinetic EEDF: no solution (test)")
        return real(self, state, target_mean, warm, growth)

    monkeypatch.setattr(KB.Solver, "_eedf_on_grid", no_warm_start)
    s.grid = cold.grid
    again = s.eedf(st, 3.0, warm=cold)
    assert s.restarts == 1
    assert again.mean == approx(cold.mean, rel=1e-6)

    def never(self, state, target_mean, warm, growth):
        raise RuntimeError("kinetic EEDF: no solution (test)")

    monkeypatch.setattr(KB.Solver, "_eedf_on_grid", never)
    with pytest.raises(RuntimeError):
        s.eedf(st, 3.0, warm=cold)
    assert s.restarts == 2
