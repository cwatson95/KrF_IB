"""v0_4's gate (LIGHT_Kr v0_4, 2026-10-15), the checks v0_5 still has to
pass (v0_5's own gate is test_v0_5.py; its comparison of the rows with
v0_3's is gone). v0_4's gate was:

* everything v0_4 does not change is v0_3's: the data files (Biagi's
  file, the trained factors, the beam tail table), the cross sections, the
  swarm, jung2006's measured tables, every row that books no 5s level
  and no excimer;
* the 5s levels resolved: four species at Biagi's thresholds, each made by
  its own G1 row; the per-level 5s -> 5p rates average to v0_3's lumped
  one by statistical weight, exactly; the near-IR branching is the NIST
  A-values', and its mean decay the measured lifetimes';
* detailed balance: the 1s5 <-> 1s4 transfer's reverse against Turner's
  measured one, the equilibrium constant (thermo.py) against the classical
  configuration integral, and the reverses at equilibrium;
* the ledger closes on both closures, and the new bands at their ends;
* the numerics: the two-temperature closure reads the field-driven tail
  table at the run's (n_e/N, Te), v0_3's beam table stays as
  Params.tail = 'beam'; the afterglow's EEDF can be time-stepped
  (Params.kin_afterglow_memory); the swarm's electroluminescence yield
  against Mano et al.'s measurement, a score the training never saw;
* v0_3's checks that still hold: the trained factors, the swarm (Einstein,
  LoKI-MC on Biagi's own file), the Hornbeck-Molnar band, the associative
  electrons in the Boltzmann equation, the solver's cold restart.

Comparisons with v0_3 skip when v0_3 is not beside v0_4."""

import importlib
import json
import math
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

from v0_5 import boltzmann as KB
from v0_5 import jung2006, network, swarm, thermo, xsec
from v0_5.model import Model, Params
from v0_5.report import accounting, row_dict
from v0_5.species import SPECIES, e_form

HERE = Path(__file__).resolve()
V0_3 = HERE.parents[2] / "v0_3"
needs_v0_3 = pytest.mark.skipif(not (V0_3 / "network.py").is_file(),
                                reason="v0_3 is not beside v0_4")
FAST = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
KT300 = 1.380649e-23 * 300.0 / 1.602176634e-19
LEVELS = network.LEVELS_5S


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


# ------------------------------------------------------------- vs v0_3
# the rows v0_4 rewrote or added: every row on a 5s level or an excimer
# state, the 5p rung's rows, the ion's reverse and the recombination's band
V0_4_ROWS = ("G1_", "G2_", "AK8", "KR2", "QKR2", "MIX_", "RES_", "KPS",
             "PRAD", "K76", "KSTEP", "PEN_KR2", "K15", "R11p")   # (v0_4's)


@needs_v0_3
def test_cross_sections_and_swarm_are_v0_3s():
    x3 = importlib.import_module("v0_3.xsec")
    s3 = importlib.import_module("v0_3.swarm")
    a, b = xsec.processes("Kr"), x3.processes("Kr")
    assert len(a) == len(b)
    for p, q in zip(a, b):
        assert (p.product, p.kind, p.threshold) == (q.product, q.kind,
                                                   q.threshold)
        assert np.array_equal(p.sigma, q.sigma)
    for en in (0.01, 1.0, 30.0):
        d, d3 = swarm.point(en), s3.point(en)
        assert (d["W"], d["alphaN"], d["mean"]) == (d3["W"], d3["alphaN"],
                                                    d3["mean"])


# ------------------------------------------------------------- the 5s levels
def test_the_5s_levels_are_biagis_and_each_g1_row_makes_its_own():
    rows = {r.rid: r for r in network.build()}
    for p in xsec.raw_processes("Kr"):
        if p.kind == "excitation" and p.threshold < e_form("Kr**"):
            lv = network.level_tag(p.product)
            assert lv in LEVELS
            assert e_form(network.kr5s(lv)) == p.threshold
            assert rows[f"G1_{lv}"].products == ("e", network.kr5s(lv))
    assert [network.kr5s(lv) for lv in LEVELS] == list(network.KR_5S)
    assert sum(network.G_LEVEL_5S.values()) == network.G_5S == 12


@needs_v0_3
def test_per_level_k76_averages_to_v0_3s_lump():
    """The four levels' 5s -> 5p rates, weighted 2J + 1, are v0_3's lumped
    rate, and so are their reverses' weights: one cross section split, not
    changed. To 1e-6: each table is read log-log in Te on its own, so the
    weighted sum of the readings is the reading of the sum only at the
    tables' nodes."""
    j3 = importlib.import_module("v0_3.jung2006")
    for te in (0.5, 1.0, 2.0, 4.0):
        mean = sum(network.G_LEVEL_5S[lv] / 12.0
                   * float(jung2006.level_rate_table(lv)(te))
                   for lv in LEVELS)
        assert mean == approx(float(j3.rate_table()(te)), rel=1e-6)
    assert sum(network.G_LEVEL_5S[lv] / network.G_5P for lv in LEVELS) == \
        approx(network.G_5S / network.G_5P, rel=1e-15)


def test_near_ir_branching_is_nists():
    """v0_4's statistical branching of NIST's A-values, which v0_5 keeps as
    a reference: v0_5 radiates line by line from the rung's levels at Tg
    (test_v0_5)."""
    b = network.KR_5P_BRANCH
    assert sum(b.values()) == approx(1.0, rel=1e-12)
    assert b == approx({"1s5": 0.4229, "1s4": 0.2351, "1s3": 0.0892,
                        "1s2": 0.2529}, rel=2e-3)
    a_mean = sum(network.KR_5P_G[u] * a for (u, _), a
                 in network.KR_5P_5S_A.items()) / sum(network.KR_5P_G.values())
    # NIST's A-values against the measured lifetimes the rung keeps
    assert 1.0 / a_mean == approx(network.KR_5P_TAU_S, rel=0.02)


# ------------------------------------------------------------- detailed balance
def test_the_1s4_1s5_transfer_reverse_is_turners():
    """MIX_1s5r by detailed balance from Turner's 1s4 -> 1s5 (1.3e-13 at
    295 K) against his measured 1s5 -> 1s4, 7.6e-16 (Phys. Rev. 158, 121
    (1967)). The forward passes through both of his points, 1.3e-13 at
    295 K and 22e-14 at 196 K."""
    m = Model(Params(eedf="two_t"))
    c = {x.row.rid: x for x in m.compiled}
    assert m._k(c["MIX_1s4"], 1.0, 295.0, None, None) == approx(1.3e-13,
                                                               rel=1e-12)
    assert m._k(c["MIX_1s4"], 1.0, 196.0, None, None) == approx(22e-14,
                                                               rel=0.01)
    k = m._k(c["MIX_1s5r"], 1.0, 295.0, None, None)
    assert k == approx(7.6e-16, rel=0.03)


def test_equilibrium_constant_against_the_classical_integral():
    """thermo.k_eq: in the classical limit of a narrow harmonic well it is
    exactly the configuration integral over sigma = 2; against the integral
    of the real quantum at 3000 K to 6 % (the rigid rotor and harmonic
    oscillator holding to kT/D0 = 0.36)."""
    from scipy.integrate import quad
    mu = 0.5 * thermo.M_KR * thermo.AMU
    g_m, we, re = thermo.MOLECULES["Kr2*"]
    nu = thermo.C_CM * we
    k_force = mu * (2.0 * math.pi * nu) ** 2
    d0 = 0.709
    de = d0 + 0.5 * we / thermo.CM_PER_EV
    re_m = re * 1e-10
    for tg, tol in ((3000.0, 0.06), (2000.0, 0.04)):
        kt = thermo.KB * tg
        half = math.sqrt(2.0 * de * thermo.EV / k_force)
        f = lambda r: 4.0 * math.pi * r * r * math.exp(
            (de * thermo.EV - 0.5 * k_force * (r - re_m) ** 2) / kt)
        integ, _ = quad(f, re_m - half, re_m + half, points=[re_m],
                        limit=400)
        k_cl = g_m / 5.0 * 0.5 * integ * 1e6
        assert thermo.k_eq("Kr2*", d0, 5, tg) == approx(k_cl, rel=tol)


def test_reverses_balance_their_forwards_at_equilibrium():
    """At equilibrium densities (thermo.k_eq) every association and its
    reverse carry the same flux, at any gas temperature."""
    m = Model(Params(eedf="two_t"))
    c = {x.row.rid: x for x in m.compiled}
    n_kr = 2.45e19
    for tg in (300.0, 800.0, 1500.0):
        for fwd, rev, mol, lv, g_a in (
                ("AK8_1s5", "AK8r_1s5", "Kr2*", "1s5", 5),
                ("AK8_1s4", "AK8r_1s4", "Kr2*", "1s4", 3),
                ("AK8f_1s4", "AK8fr_1s4", "Kr2*(0u+)", "1s4", 3)):
            d0 = e_form(network.kr5s(lv)) - e_form(mol)
            k_eq = thermo.k_eq(mol, d0, g_a, tg)
            kf = m._k(c[fwd], 1.0, tg, None, None) / max(c[fwd].factor, 1e-300)
            kr = m._k(c[rev], 1.0, tg, None, None) / max(c[rev].factor, 1e-300)
            n_a = 1e12
            n_m = k_eq * n_a * n_kr
            assert kf * n_a * n_kr ** 2 == approx(kr * n_m * n_kr, rel=1e-12)
        k15, k15r = c["K15"], c["K15r"]
        kf = m._k(k15, 1.0, tg, None, None)
        kr = m._k(k15r, 1.0, tg, None, None)
        k_eq = thermo.k_eq("Kr2+", e_form("Kr+") - e_form("Kr2+"), 4, tg)
        assert kf / kr == approx(k_eq, rel=1e-12)


def test_redissociation_negligible_cold_dominant_hot():
    m = Model(Params(eedf="two_t"))
    c = next(x for x in m.compiled if x.row.rid == "AK8r_1s5")
    n_kr = 2.45e19
    rad = 3.7e6
    assert m._k(c, 1.0, 300.0, None, None) * n_kr < 1e-6 * rad
    assert m._k(c, 1.0, 1500.0, None, None) * n_kr > rad
    lo = Model(Params(eedf="two_t", rates=(("kr2_d0", 0.62),)))
    c_lo = next(x for x in lo.compiled if x.row.rid == "AK8r_1s5")
    assert lo._k(c_lo, 1.0, 800.0, None, None) > 2.0 * m._k(c, 1.0, 800.0,
                                                             None, None)


# ------------------------------------------------------------- the runs
def test_ledger_closes_and_the_new_bands_act():
    base = Model(Params(**FAST)).run()
    a = accounting(base)
    assert abs(a["closure"]) < 1e-12
    d0 = row_dict(base)
    # the fast state takes its share of the 1s4 association and radiates
    fast = Model(Params(**FAST, rates=(("kr2_fast", 0.17),))).run()
    assert fast.events["KR2F_RAD"][-1] > 0.0
    assert base.events["KR2F_RAD"][-1] == 0.0
    assert abs(accounting(fast)["closure"]) < 1e-12
    # quenching by Kr at its top costs light
    q = Model(Params(**FAST, rates=(("kr2_quench", 8.6e-15),))).run()
    assert row_dict(q)["emission_vuv"] < d0["emission_vuv"]
    # the resonance lines radiate only when the band lets them escape
    assert a["bands"].get("123.6 nm (resonance)", 0.0) == 0.0
    esc = Model(Params(**FAST, rates=(("res_escape", 1.7e5),))).run()
    assert accounting(esc)["bands"]["123.6 nm (resonance)"] > 0.0
    # all relaxation into 1s5: in v0_4 it slowed the excimer's formation
    # and cost light; in v0_5 it moves the share the other way (the
    # overview's band ends and its reading of them)
    low = Model(Params(**FAST, rates=(("kps_5s", 0.0),))).run()
    assert row_dict(low)["emission_vuv"] > d0["emission_vuv"]


def test_recombination_band_is_one_at_300_k():
    m0 = Model(Params(eedf="two_t"))
    m1 = Model(Params(eedf="two_t", rates=(("dr_tg", 1.0),)))
    c0 = next(x for x in m0.compiled if x.row.rid == "R11p")
    c1 = next(x for x in m1.compiled if x.row.rid == "R11p")
    assert m1._k(c1, 1.0, 300.0, None, None) == approx(
        m0._k(c0, 1.0, 300.0, None, None), rel=1e-12)
    assert m1._k(c1, 1.0, 900.0, None, None) == approx(
        4.80 * m0._k(c0, 1.0, 900.0, None, None), rel=1e-9)


# ------------------------------------------------------------- the factors
def test_trained_factors_are_the_pinned_files():
    name, md5 = xsec.TRAINED_FILE
    assert md5 is not None, "v0_4 runs on trained factors"
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
    y[m.idx["Kr*(1s5)"]] = 1e-5
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


# ------------------------------------------------------------ the numerics
def test_two_t_reads_the_field_table_at_the_runs_state():
    """Each tail row's rate under the two-temperature closure is its
    two-temperature table at T_t = r Te, r the field table's at the run's
    (n_e/N, Te); Params.tail = 'beam' reads v0_3's table at its F2 = 0
    column, whatever Te. Anything else is refused."""
    n_tot = 2.45e19
    for tail, r_of in (("field", lambda x, te: xsec.tail_field_table()(x, te)),
                       ("beam", lambda x, te: xsec.tail_ratio_table()(x, 0.0))):
        m = Model(Params(eedf="two_t", tail=tail))
        c = {x.row.rid: x for x in m.compiled}
        n = np.zeros(m.ns)
        n[m.idx["Kr"]] = n_tot
        for ne, te in ((3e14, 2.3), (1e13, 3.0), (1e17, 1.2)):
            r = r_of(ne / n_tot, te)
            for rid, prod in (("G1_1s5", "Kr(1S5)"), ("C2", None)):
                kind = "ionization" if rid == "C2" else "excitation"
                want = xsec.two_t_table("Kr", kind, prod)(te, r * te)
                got = m._k(c[rid], te, 300.0, n, ne)
                assert got == approx(want, rel=1e-12), (tail, rid, ne, te)
    a = Model(Params(eedf="two_t"))._k
    with pytest.raises(ValueError):
        Model(Params(eedf="two_t", tail="hot"))


def test_the_field_table_brings_the_closure_toward_the_kinetic_eedf():
    """v0_3's beam table put the two-temperature tail too hot (r 0.18-0.20
    against the kinetic EEDF's 0.12-0.16 along the reference run); the
    field table's colder tail lowers the closure's 147-nm share and
    raises Te (the electrons need a hotter bulk to lose the same energy)."""
    f = row_dict(Model(Params(eedf="two_t")).run())
    b = row_dict(Model(Params(eedf="two_t", tail="beam")).run())
    assert f["emission_vuv"] < b["emission_vuv"]
    assert f["te_mean_eV"] > b["te_mean_eV"]


def test_the_afterglow_eedf_can_be_time_stepped():
    """Params.kin_afterglow_memory keeps the kinetic EEDF's time step after
    the pulse (the log's 'memory'), where by default it is steady; on a
    short pulse the step after it moves the bands little and the ledger
    closes either way."""
    kw = dict(eedf="kinetic", pulse_s=40e-9, t_end_s=80e-9, n_save=81)
    out = {}
    for flag in (False, True):
        r = Model(Params(kin_afterglow_memory=flag, **kw)).run()
        assert abs(accounting(r)["closure"]) < 1e-9
        after = [e for e in r.kinetic if not e["on"] and not e.get("skipped")]
        assert after and all(e["memory"] is flag for e in after)
        assert all(e["memory"] for e in r.kinetic
                   if e["on"] and not e.get("skipped"))
        out[flag] = row_dict(r)
    for k in ("emission_vuv", "emission_nir"):
        assert out[True][k] == approx(out[False][k], rel=1e-2), k


def test_electroluminescence_yield_against_mano():
    """The swarm's excitations per electron and per cm of drift, one photon
    each, against Mano et al.'s krypton yield at 1.1 bar (Phys. Lett. B
    824, 136782 (2022), the fit 0.113 E/N - 0.30 in 1e-17 cm^2): within 5
    % of it from 6 to 13.5 Td, nearly every excitation a 5s level's, and
    below the fit toward its 2.7-Td intercept (2026-10-15: 0.99 to 1.04,
    and 0.30 at 3 Td). Nothing ionizes below 14 Td."""
    for en in (6.0, 8.0, 10.0, 13.5):
        y = swarm.el_yield(en)
        assert 0.95 < y["ratio"] < 1.05, (en, y["ratio"])
        assert y["share_5s"] > 0.97
        assert y["alphaN"] < 1e-5 * y["model"]
    assert swarm.el_yield(3.0)["ratio"] < 0.5
    # photons per excitation scale it
    y1, y2 = swarm.el_yield(10.0), swarm.el_yield(10.0, photons={"1s": 0.5})
    assert y2["model"] == approx(y1["model"] * (0.5 * y1["share_5s"]
                                                + 1 - y1["share_5s"]),
                                 rel=1e-9)


# the state where the scan's hottest corner (0.25 atm, 500 kW/cm3) stopped
# on 2026-10-15: 52 ns into the afterglow of 3970-K gas, a table point
# below the no-field mean, told to leave the growth term out (its centre
# had none), on the centre's 16-eV grid. The Penning frequencies are the
# run's; their energies are the rows' own (a change in the fifth decimal
# moves the no-growth solve across its convergence edge)
HOT_PEN = (585953.5883000462, 511468.321801811, 157960.7814724454,
           397233.28362642805, 111612.868932891, 68940.55896526361,
           173368.88534949513, 10645.727315952887, 53542.875375622294,
           67323.71162638624, 0.0)
HOT_STATE = dict(
    x={"Kr": 1.0}, tg=3967.5796566349404, n_col=6.104764900811287e+24,
    ne=5.913018902430488e+21, rec=462733.06272305077, rec_pow=-0.55,
    n_ion=5.913018902430488e+21,
    pops={"Kr*(1s5)": 0.00043120351956801915,
          "Kr*(1s4)": 0.0001881948885647534,
          "Kr*(1s3)": 5.812170646675417e-05,
          "Kr*(1s2)": 0.00014616207956522165,
          "Kr**": 1.1467680361319614e-05})
HOT_TARGET, HOT_EMAX = 0.3604272572559662, 15.972000000000005


def _hot():
    m = Model(Params(pressure_atm=0.25, power_W_cm3=5e5))
    pen = tuple((f, c.route["to_e"]) for f, c in zip(HOT_PEN, m._kin_pen))
    s = KB.State(pen=pen, **HOT_STATE)

    def solver():
        sv = KB.Solver(m.kin.procs, n_cells=m.p.kin_cells,
                       emax=m.p.kin_emax_eV, tol=m.p.kin_tol,
                       elastic=("Kr",), pops=m.kin.pops)
        sv.grid = sv._grid(HOT_EMAX)
        return sv
    return s, solver


def test_a_table_point_without_a_steady_solution_falls_back(monkeypatch):
    """In 4000-K gas's afterglow Penning ionization makes electrons four
    times faster than recombination takes them, so a table point told to
    leave the growth term out may have no steady solution: v0_3's solver,
    and v0_4's before 2026-10-15, stopped the run there. The point now
    falls back to the growth term on the same grid, where it is the
    no-field solution, below which no field takes the mean; failing that,
    to a grid of its own (each counted in Solver.fallbacks)."""
    s, solver = _hot()
    assert sum(v for v, _ in s.pen) > 4.0 * s.rec
    # the run's own state: never a stop
    sv = solver()
    sol = sv.eedf(s, target_mean=HOT_TARGET, warm=None, growth=False,
                  adapt=False)
    assert not sol.reached and sol.en_td == 0.0
    assert sol.mean > 1.1 * HOT_TARGET
    ref = solver().eedf(s, target_mean=HOT_TARGET, warm=None, growth=True,
                        adapt=False)
    assert not ref.reached and not ref.growth_off
    # no solution without the growth term: the growth term's, on its grid
    real = KB.Solver._eedf_on_grid

    def no_steady(self, st, target, warm, growth):
        if not growth:
            raise RuntimeError("no steady solution (test)")
        return real(self, st, target, warm, growth)
    monkeypatch.setattr(KB.Solver, "_eedf_on_grid", no_steady)
    sv = solver()
    sol = sv.eedf(s, target_mean=HOT_TARGET, warm=None, growth=False,
                  adapt=False)
    assert (sv.fallbacks, sv.restarts) == (1, 0)
    assert sol.mean == approx(ref.mean, rel=1e-9)
    assert sol.grid.emax == approx(HOT_EMAX, rel=1e-12)
    # nothing on the centre's grid: a grid of its own

    def not_here(self, st, target, warm, growth):
        if abs(self.grid.emax - HOT_EMAX) < 1e-9:
            raise RuntimeError("no solution on this grid (test)")
        return real(self, st, target, warm, growth)
    monkeypatch.setattr(KB.Solver, "_eedf_on_grid", not_here)
    sv = solver()
    sol = sv.eedf(s, target_mean=HOT_TARGET, warm=None, growth=False,
                  adapt=False)
    assert sv.fallbacks == 2
    assert sol.grid.emax != approx(HOT_EMAX, rel=1e-9)
    assert not sol.reached
    # a state with no solution anywhere still stops, with its message
    monkeypatch.setattr(KB.Solver, "_eedf_on_grid",
                        lambda *a: (_ for _ in ()).throw(
                            RuntimeError("nowhere (test)")))
    with pytest.raises(RuntimeError, match="nowhere"):
        solver().eedf(s, target_mean=HOT_TARGET, warm=None, growth=True,
                      adapt=True)
