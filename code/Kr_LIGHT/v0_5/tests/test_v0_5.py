"""LIGHT_Kr v0_5 (2026-10-16): v0_4 with the limitations of its accuracy
overview taken up. The gate:

* everything v0_5 does not change is v0_4's: the data files but the field
  tail table, the cross sections and the swarm, jung2006's tables, every
  row v0_5 did not write or re-cite, every species but Kr2+ (Signorell et
  al.'s D0); the field table's x_s = 0 block is v0_4's table;
* the pumps pair by pair: out of each level they sum to v0_4's level rate,
  the energy an electron pays per pump is the measured cross sections'
  (v0_4 undercharged it), and each pair balances its reverse out of the
  rung's level at Te = Tg;
* Vriens & Smeets's cross section integrates to their rate formula, and its
  peak for Kr2* is Flannery et al.'s order; the pseudo cross sections give
  their rate laws exactly;
* in the kinetic EEDF every population operator gives or takes rate x u
  of the electrons' energy, and the population ionizations feed the growth
  term;
* the bands v0_5 adds act on the answer in their direction, the 5p rung's
  reverse in hot gas alone;
* the superelastic returns against Blagoev et al.'s afterglow rates (a
  score, pinned: 5-10 times theirs for 3P2 and 3P1, theirs for 3P0).

Comparisons with v0_4 skip when v0_4 is not beside v0_5."""

import importlib
import math
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

from v0_5 import boltzmann as KB
from v0_5 import excited, jung2006, network, swarm, xsec
from v0_5.model import KB_EV, Model, Params
from v0_5.report import accounting, row_dict
from v0_5.species import SPECIES, e_form

HERE = Path(__file__).resolve()
V0_4 = HERE.parents[2] / "v0_4"
needs_v0_4 = pytest.mark.skipif(not (V0_4 / "network.py").is_file(),
                                reason="v0_4 is not beside v0_5")
FAST = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
LEVELS = network.LEVELS_5S
# the rows v0_5 wrote or added: the 5s -> 5p pairs, the lines, the
# stepwise ionizations, the 5p rung's reverse, the superelastic returns,
# the 5s mixing, the electrons on Kr2*
V0_5_ROWS = ("K76", "PRAD", "KSTEP", "KPSr", "SUPER", "EMIX", "EI_", "EQ_",
             "EM_")
# rows re-cited, their rates and species unchanged: Kholin's third body
# (AK8_1s3) and the Hornbeck-Molnar rows' height above Kr2+ (its D0)
RECITED = ("AK8_1s3",)


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


# ------------------------------------------------------------- vs v0_4
@needs_v0_4
def test_everything_else_is_v0_4s():
    old_net = importlib.import_module("v0_4.network")
    old = {r.rid: r for r in old_net.build(("all",))}
    new = {r.rid: r for r in network.build(("all",))}
    kept = [rid for rid in new if not rid.startswith(V0_5_ROWS)]
    assert kept and all(rid in old for rid in kept)
    gone = [rid for rid in old if rid not in new]
    assert sorted(gone) == sorted(
        [f"PRAD_{lv}" for lv in LEVELS] + ["KSTEP_P"]
        + [x for lv in LEVELS for x in (f"K76_{lv}", f"K76r_{lv}")])
    for rid in kept:
        a, b = astuple(old[rid]), astuple(new[rid])
        n = len(a)
        if rid in RECITED or rid.startswith("HM_"):
            assert (a[:5], a[6:]) == (b[:5], b[6:n]), rid    # the cite
        else:
            assert a == b[:n], rid
        assert new[rid].rung_level is None and new[rid].e_gain_eV is None
    for name in SPECIES:
        if name != "Kr2+":
            assert astuple(SPECIES[name]) == astuple(importlib.import_module(
                "v0_4.species").SPECIES[name]), name
    j4 = importlib.import_module("v0_4.jung2006")
    for name in ("LEVELS", "PROXY", "SHAPE_REF", "TABLES", "TABLE_II",
                 "TABLE_III", "UNCERTAINTY", "WEIGHTS_5S", "E_MAX_TABLE",
                 "CM2"):
        assert getattr(jung2006, name) == getattr(j4, name), name
    for f in ("data/kr_tail_ratio.txt", "data/Kr_Biagi897_LXCat.txt",
              "data/kr_swarm_trained.json"):
        assert (HERE.parents[1] / f).read_bytes() == (V0_4 / f).read_bytes()
    t4 = importlib.import_module("v0_4.xsec").tail_field_table()
    assert np.array_equal(xsec.tail_field_table().r[:, :, 0], t4.r)


@needs_v0_4
def test_cross_sections_and_swarm_are_v0_4s():
    x4 = importlib.import_module("v0_4.xsec")
    s4 = importlib.import_module("v0_4.swarm")
    a, b = xsec.processes("Kr"), x4.processes("Kr")
    assert len(a) == len(b)
    for p, q in zip(a, b):
        assert (p.product, p.kind, p.threshold) == (q.product, q.kind,
                                                   q.threshold)
        assert np.array_equal(p.sigma, q.sigma)
    for en in (0.01, 1.0, 30.0):
        d, d4 = swarm.point(en), s4.point(en)
        assert (d["W"], d["alphaN"], d["mean"]) == (d4["W"], d4["alphaN"],
                                                    d4["mean"])


def test_kr2_ion_is_signorells():
    assert e_form("Kr2+") == approx(13.9996 - 9267.8 / 8065.543937,
                                    rel=1e-12)
    assert "Signorell" in SPECIES["Kr2+"].cite
    by = {r.rid: r for r in network.build()}
    assert "argon is the third body" in by["AK8_1s3"].cite


# ------------------------------------------------------------- the pumps
def test_pairs_sum_to_the_level_and_charge_the_measured_gap():
    """Out of each level the pairs' Maxwellian rates sum to v0_4's level
    rate (the onsets at Biagi's energies, within 1 %), and the energy an
    electron pays per pump, the pairs' thresholds weighted by their rates,
    is 1.536, 1.419, 1.560 and 1.478 eV at 2.6 eV (on Biagi's level
    energies; 1.540, 1.420, 1.561 and 1.481 on Table I's two decimals, the
    groundwork's), where v0_4 charged the gap to 2p10: 1.388, 1.271, 0.741
    and 0.660 eV."""
    want = {"1s5": 1.5363, "1s4": 1.4191, "1s3": 1.5597, "1s2": 1.4784}
    for lv in LEVELS:
        finals = jung2006.pair_finals(lv)
        for te in (1.0, 2.6, 5.0):
            k = {f: float(jung2006.pair_rate_table(lv, f)(te))
                 for f in finals}
            assert sum(k.values()) == approx(
                float(jung2006.level_rate_table(lv)(te)), rel=0.01), lv
            if te == 2.6:
                mean = sum(k[f] * jung2006.pair_threshold(lv, f)
                           for f in finals) / sum(k.values())
                assert mean == approx(want[lv], abs=2e-4), lv
                assert mean > e_form("Kr**") - e_form(f"Kr*({lv})")


def test_each_pair_balances_its_reverse_out_of_the_rungs_level():
    """At Te = Tg = T with the rung's levels in Boltzmann against the 5s
    level, every pair's forward and reverse carry the same flux (the
    Klein-Rosseland tables to their quadrature)."""
    for te in (1.0, 2.5):
        for lv in LEVELS:
            for f in jung2006.pair_finals(lv):
                kf = float(jung2006.pair_rate_table(lv, f)(te))
                kr = float(jung2006.pair_superelastic_table(lv, f)(te))
                thr = jung2006.pair_threshold(lv, f)
                n_ratio = (excited.RUNG[f][1] / network.G_LEVEL_5S[lv]
                           * math.exp(-thr / te))
                assert kf == approx(kr * n_ratio, rel=3e-3), (lv, f, te)


def test_the_rung_at_the_gas_temperature():
    """Z(Tg) runs from 3 (2p10 alone) to 36 (every level by weight); the
    shares sum to 1; the hot gas puts the 5p' levels within reach."""
    assert excited.rung_partition(50.0) == approx(3.0, rel=1e-9)
    assert excited.rung_partition(1e9) == approx(36.0, rel=1e-5)
    for tg in (300.0, 1213.0, 4000.0):
        assert sum(excited.rung_fraction(x, tg) for x in excited.RUNG) == \
            approx(1.0, rel=1e-12)
    assert excited.rung_fraction("2p10", 300.0) > 0.98
    assert excited.rung_fraction("2p2", 300.0) < 1e-12
    assert excited.rung_fraction("2p2", 4000.0) > 0.01
    assert excited.rung_mean_excess(700.0) == approx(0.04699, rel=1e-3)


# ------------------------------------------------------------- calculations
def test_vriens_smeets_cross_section_integrates_to_their_rate():
    """Their cross section and their rate formula, printed separately,
    agree: to 0.2 % at Te 2-5 eV, 0.8 % at 1 eV and 3.5 % at 0.5 eV (the
    1s5 level's, the worst), for the 5s levels,
    the rung's bottom and Kr2*. Its peak for Kr2*, 1.0e-15 cm2 near 10 eV,
    is Flannery, McCann & Winter's order (about 1e-15 cm2)."""
    for e_n in (13.9996 - 9.9152, 13.9996 - 10.6436, 13.9996 - 11.3035,
                e_form("Kr2+") - e_form("Kr2*")):
        t = excited.vs_rate_table(e_n)
        for te, tol in ((0.5, 0.035), (1.0, 0.01), (2.0, 2e-3),
                        (3.0, 2e-3), (5.0, 2e-3)):
            assert float(t(te)) == approx(excited.vs_rate(e_n, te),
                                          rel=tol), (e_n, te)
    e_n = e_form("Kr2+") - e_form("Kr2*")
    e = np.linspace(e_n, 60.0, 20000)
    s = excited.vs_sigma(e, e_n)
    assert s.max() == approx(1.0e-15, rel=0.05)
    assert 9.0 < e[np.argmax(s)] < 11.0
    assert excited.vs_sigma(np.array([e_n]), e_n)[0] == 0.0


def test_pseudo_cross_sections_give_their_rate_laws():
    """A Te^B exp(-C/Te) is the Maxwellian of sigma0 (eps - C)^(B + 1/2) /
    eps exactly, and a constant superelastic rate the Klein-Rosseland
    reverse of excited.gain_process."""
    for a, b, c in ((4.3e-8, 0.0, 0.1172), (9e-7, 0.5, 1.51),
                    (1e-7, 0.0, 0.142)):
        p = excited.arrhenius_process(a, b, c, c)
        k = xsec._mean(p.sigma_at(xsec.EPS) * xsec.V_CM_S
                       * np.sqrt(xsec.EPS))
        tab = xsec.RateTable(k, "")
        for te in (0.3, 1.0, 2.5, 5.0):
            assert float(tab(te)) == approx(a * te ** b * math.exp(-c / te),
                                            rel=1e-3), (a, b, c, te)
    k0, u = 1e-8, 8.43
    p = excited.gain_process(k0, u, 1.0)
    e = xsec.EPS
    sig_r = (e + u) / e * p.sigma_at(e + u)
    kr = xsec.RateTable(xsec._mean(sig_r * xsec.V_CM_S * np.sqrt(e)), "")
    for te in (0.5, 1.5, 5.0):
        assert float(kr(te)) == approx(k0, rel=1e-2), te


# ------------------------------------------------------------- kinetic EEDF
def _state_mid_pulse():
    two = Model(Params(eedf="two_t"))
    res = two.run()
    k = int(np.searchsorted(res.t, 0.3e-6))
    y = np.zeros(two.size)
    for s, i in two.idx.items():
        y[i] = res.n[s][k] / two.N0
    y[two.i_ue] = res.ue[k] / two.N0
    y[two.i_tg] = res.tg[k]
    return y


def test_every_population_operator_exchanges_rate_times_u():
    """On a solved EEDF at a mid-pulse state (5s, 5p and Kr2* populated),
    each population operator changes the electrons' energy by its own rate
    times its threshold: the pairs, the reverses, the superelastic returns
    (10 eV up, on the run's grid), the ionizations and the mixing; and the
    ionizations make electrons one for one."""
    m = Model(Params())
    st, te, ne = m._kinetic_state(_state_mid_pulse())
    assert st.pops["Kr2*"] > 0.0 and st.pops["Kr**:2p10"] > 0.0
    sol = m.kin.eedf(st, 1.5 * te)
    g, F = sol.grid, sol.F
    seen = set()
    for name, q in g.pops.items():
        op_f, op_r = g.op_pop[name]
        kf, kr = g.pop_rates(name, F)
        for op, rate, sign in ((op_f, kf, -1.0), (op_r, kr, 1.0)):
            if op is None:
                continue
            dF = -(op @ F)
            assert float(np.sum(g.ec * dF)) == approx(
                sign * q.p.threshold * rate, rel=2e-3), name
            n_new = float(np.sum(dF))
            if q.kind == "ionization":
                assert n_new == approx(rate, rel=1e-9), name
            else:
                assert abs(n_new) < 1e-9 * rate, name
            seen.add(name.split("_")[0])
    assert {"K76", "K76r", "SUPER", "KSTEP", "EI", "EMIX", "EMIXr"} <= seen
    # the ionizations out of the excited populations feed the growth term
    ion0 = m.kin._frequencies(g, st, F)[0]
    st0 = KB.State(**{**st.__dict__, "pops": {}})
    ion1 = m.kin._frequencies(g, st0, F)[0]
    extra = sum(g.rate(q.p, F) * st.pops.get(q.lower_pop, 0.0) * st.n_col
                for q in g.pop_ion)
    assert ion0 - ion1 == approx(extra, rel=1e-9) and extra > 0.0


# ------------------------------------------------------------- the bands
def test_the_new_bands_act_in_their_direction():
    """The electrons' ionization of Kr2* costs light (at 0 the share more
    than doubles), so do its quench and the stepwise ionizations; the
    mixing of its states toward the fast one gives some back; the measured
    5s mixing at its top moves the share by under 2 %."""
    d = row_dict(Model(Params(**FAST)).run())["emission_vuv"]

    def at(**rates):
        return row_dict(Model(Params(**FAST, rates=tuple(rates.items())))
                        .run())["emission_vuv"]
    assert at(ei_kr2=0.0) > 1.5 * d > d > at(ei_kr2=2.0)
    assert at(eq_kr2=1e-7) < 0.5 * d
    assert at(em_kr2=1e-7) > d
    assert at(kstep_s=0.0) > d > at(kstep_s=2.0)
    assert at(kstep_p=0.0) > d > at(kstep_p=2.0)
    assert at(emix_5s=2.3e-7) == approx(d, rel=0.02)


def test_the_rungs_reverse_acts_in_hot_gas_alone():
    """KPSr: under 1e-3 of the relaxation at the reference point (700-860
    K); at 0.5 atm and 500 kW/cm3 (the gas at 2200 K) about 2 %, its gap
    drawn from the gas, which ends the cooler for it."""
    def run_(**kw):
        return Model(Params(eedf="two_t", **kw)).run()

    def share(r):
        return sum(r.events[f"KPSr_{lv}"][-1] for lv in LEVELS) / sum(
            r.events[f"KPS_{lv}"][-1] for lv in LEVELS)
    off = tuple((f"KPSr_{lv}", 0.0) for lv in LEVELS) + (("KPSr_low", 0.0),)
    assert share(run_()) < 1e-3
    hot = run_(pressure_atm=0.5, power_W_cm3=5e5)
    cold = run_(pressure_atm=0.5, power_W_cm3=5e5, row_scale=off)
    assert hot.tg.max() > 2000.0 and 0.01 < share(hot) < 0.05
    assert share(cold) == 0.0 and hot.tg[-1] < cold.tg[-1]
    assert abs(accounting(hot)["closure"]) < 1e-9


# ------------------------------------------------------------- the score
# Blagoev, Mishonov & Popov, J. Phys. B 17, 435 (1984), abstract: thermal
# electrons' superelastic rates in a krypton afterglow [cm3/s]
BLAGOEV = {"1s5": 1.5e-10, "1s4": 5.2e-10, "1s3": 9.8e-10}


def test_superelastic_returns_against_blagoev():
    """The Klein-Rosseland returns on the trained G1 cross sections at
    thermal electrons (300-2000 K): 5.3-6.6 times Blagoev et al.'s for
    3P2, 5.5-9.5 for 3P1, 0.55-1.24 for 3P0 (pinned). The returns move the
    reference point by 0.4 % (the 2026-10-16 groundwork), so the
    disagreement, in the cross sections' first 0.2 eV above threshold, does
    not reach the answer."""
    m = Model(Params(eedf="two_t"))
    c = {x.row.rid: x for x in m.compiled}
    ratio = {lv: [m._k(c[f"SUPER_{lv}"], KB_EV * t, 300.0, None, None)
                  / v for t in (300.0, 500.0, 1000.0, 2000.0)]
             for lv, v in BLAGOEV.items()}
    assert ratio["1s5"] == approx([5.87, 6.60, 6.42, 5.30], rel=0.01)
    assert ratio["1s4"] == approx([9.52, 8.95, 7.72, 5.53], rel=0.01)
    assert ratio["1s3"] == approx([1.24, 1.03, 0.74, 0.55], rel=0.01)


def test_the_table_survives_a_point_or_a_centre_with_no_eedf():
    """v0_5: a table point with no kinetic EEDF, after every restart and
    fallback, is left out of the table (counted 'dropped'), and an update
    whose centre has none leaves LIGHT's Maxwellian in charge until the
    next renewal (logged 'failed'). Before, either stopped the run: the
    scan's 0.25-atm, 500-kW corner, a side point in its afterglow
    (2026-10-16). Steady table and time step alike."""
    for memory in (0.0, math.inf):
        m = Model(Params(kin_memory_s=memory))
        y = m.initial_state()
        y[m.idx["Kr2+"]] = y[m.idx["Kr+"]]        # recombination at work
        n = np.maximum(y[:m.ns], 0.0) * m.N0
        ne = m.electron_density(n)
        y[m.i_ue] = 1.5 * ne * 2.0 / m.N0         # the table at Te = 2 eV
        if memory:
            m._kinetic_update(0.0, y, True)       # the start: Maxwellian
        name = "step" if memory else "eedf"
        solve = getattr(m.kin, name)
        low = 3.0 * math.exp(-0.12)    # the lowest side point, 3 eV e^-0.16

        def failing(s, *a, **kw):
            # step takes the mean by keyword, eedf as its second argument
            mean = kw.get("target_mean", a[0] if a and not memory else None)
            if mean is not None and mean < low:
                raise RuntimeError("no solution (test)")
            return solve(s, *a, **kw)
        setattr(m.kin, name, failing)
        m._kinetic_update(1e-12, y, True)
        log = m.kin_log[-1]
        assert log["dropped"] == 1 and not log.get("failed"), memory
        lnte = m._tab["lnte"]
        assert len(lnte) == 4
        assert lnte[0] > math.log(2.0) - 0.12     # the point left out
        # a centre with no solution: the update leaves no table, and every
        # row reads LIGHT's Maxwellian until the next renewal

        def none(*a, **kw):
            raise RuntimeError("no solution (test)")
        setattr(m.kin, name, none)
        m._kinetic_update(2e-12, y, True)
        log = m.kin_log[-1]
        assert log["failed"] and log["skipped"] and m._tab is None
        tg = y[m.i_tg]
        for c in m.compiled:
            if c.row.rid in ("C2", "G1_1s5", "KSTEP_S_1s5", "EI_KR2"):
                assert m._rate(c, 2.0, tg, n, ne) == m._k(c, 2.0, tg, n,
                                                          ne)
