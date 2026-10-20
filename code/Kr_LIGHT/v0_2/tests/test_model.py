"""The 0-D pure-Kr model: ledger closures, analytic limits, convergence,
CLI.

These are LIGHT v0_5's model tests carried to krypton. They run on LIGHT
v0_2's electron closure (eedf='two_t', the Maxwellian with the tabulated
tail) unless a test names another: Params and run below default to it. The
kinetic EEDF, the model's own default, has its tests in test_kinetic.py."""

import math

import numpy as np
import pytest

from v0_2.model import KB_EV, QE, Model
from v0_2.model import Params as _Params
from v0_2.model import run as _run
from v0_2.report import accounting, reservoir, row_dict, summary


def Params(**kw):
    kw.setdefault("eedf", "two_t")
    return _Params(**kw)


def run(params=None, **kw):
    if params is None:
        kw.setdefault("eedf", "two_t")
    return _run(params, **kw)


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


FAST = dict(pulse_s=0.3e-6, t_end_s=0.6e-6)


def _heat_booked(res):
    m = res.model
    return (res.ledger["heat_elastic"] + res.ledger["heat_capture"]
            + sum(res.events[c.row.rid] * c.route["heat"] for c in m.compiled))


def _check_closures(res, tol=1e-9):
    a = accounting(res)
    assert abs(a["closure"]) < tol                       # chemical ledger
    e_dep = a["e_dep"]
    # the same closure recomputed from the raw arrays, not from report.py
    rad = sum(res.events[c.row.rid][-1] * c.route["rad"]
              for c in res.model.compiled)
    raw = (res.ledger["E_dep"][-1] - (res.stored[-1] - res.stored[0]) - rad
           - _heat_booked(res)[-1]) / e_dep
    assert abs(raw) < tol and raw == approx(a["closure"], abs=1e-12)
    gap = np.abs((res.u_thermal - res.u_thermal[0]) - _heat_booked(res))
    assert gap.max() / e_dep < tol                       # gas thermal energy
    for name in ("Kr*", "Kr**", "Kr2*"):
        rv = reservoir(res, name)
        assert abs(rv["balance"]) <= 1e-6 * max(rv["formed"], 1.0)


def test_reference_point_closes():
    res = run()
    _check_closures(res)
    e_dep_j = res.ledger["E_dep"][-1] * QE
    assert e_dep_j == approx(0.3, rel=1e-9)       # 300 kW/cm^3 x 1 us


@pytest.mark.parametrize("kw", [
    dict(groups=("all",)),
    dict(pressure_atm=0.1, power_W_cm3=5e5, ne_cm3=1e14),
    dict(pressure_atm=2.0, power_W_cm3=1e5, ne_cm3=1e16),
    dict(eedf="maxwellian"),
    dict(eedf="maxwellian", capture="three_halves"),
])
def test_closures_across_modes(kw):
    _check_closures(run(**FAST, **kw))


def test_n_e_is_solved_and_quasineutral():
    """No electrons are injected: n_e starts at the seed and is the
    cations' sum. The parameters only fluorine needed are gone."""
    res = run(**FAST, groups=("all",))
    assert res.model.ne_mode == "free"
    assert res.ne[0] == approx(1e15, rel=1e-12)   # the seed
    # the model's own rule: densities clipped at zero (the solver's noise
    # near zero is within its absolute tolerance, 1e-18 N0 = 25 cm^-3),
    # the cations' sum, floored at 1 cm^-3
    cat = sum(np.maximum(res.n[s], 0.0) for s in ("Kr+", "Kr2+"))
    assert np.allclose(res.ne, np.maximum(cat, 1.0), rtol=1e-12)
    for gone in ("x_f2", "ne_mode", "tau_sp_s", "trimer_tg_pow", "network"):
        with pytest.raises(TypeError):
            _Params(**{gone: 1})


def test_dilute_kr_star_seed_ends_as_147_nm_photons():
    """No drive and no electrons to speak of (a 1 cm^-3 seed): a dilute
    Kr* seed can only form Kr2* (Kr* + Kr* is second order in the seed) and
    Kr2* only radiates (Kr2* + Kr2* is 0 by default). Every Kr* becomes one
    147-nm photon, on the two-step chain's analytic time course with
    Ueda's 1.7e-32 cm6/s and the 270-ns lifetime at 300 K."""
    p = Params(power_W_cm3=0.0, init_extra_cm3=(("Kr*", 1e10),),
               pulse_s=0.0, t_end_s=3e-6, ne_cm3=1.0, n_save=301)
    res = Model(p).run()
    ev = res.events
    assert ev["AK8_Kr"][-1] / 1e10 == approx(1.0, rel=1e-5)
    k3 = 1.7e-32 * res.n["Kr"][0] ** 2
    a = 3.7e6
    t = res.t
    chain = 1.0 - (a * np.exp(-k3 * t) - k3 * np.exp(-a * t)) / (a - k3)
    assert np.max(np.abs(ev["KR2_RAD"] / 1e10 - chain)) < 1e-5
    rad = accounting(res)["bands"]["147 nm (Kr2* VUV)"]
    assert rad == approx(8.43 * ev["KR2_RAD"][-1], rel=1e-12)


def test_a_5p_seed_relaxes_or_radiates_in_the_analytic_ratio():
    """No drive, no electrons: a Kr** seed relaxes by Kr (KPS) or radiates
    the near-IR lines (PRAD). The radiated share is A / (A + k_ps N) with
    the 29-ns lifetime and Gazeli et al.'s 1.7e-10 cm3/s, and every one
    ends in Kr*."""
    from v0_2.network import KR_5P_TAU_S
    p = Params(power_W_cm3=0.0, init_extra_cm3=(("Kr**", 1e10),),
               pulse_s=0.0, t_end_s=50e-9, ne_cm3=1.0)
    res = Model(p).run()
    a = 1.0 / KR_5P_TAU_S
    kps = 1.7e-10 * res.n["Kr"][0]
    ev = res.events
    assert ev["PRAD"][-1] / 1e10 == approx(a / (a + kps), rel=1e-5)
    assert (ev["PRAD"][-1] + ev["KPS"][-1]) / 1e10 == approx(1.0, rel=1e-6)


def test_electron_energy_is_quasi_steady_mid_pulse():
    """Te relaxes in ~ns; mid-pulse the heating equals the losses: the
    electron energy changes by 2e-5 of the input power at 0.5 us."""
    m = Model(Params())
    res = m.run()
    k = np.searchsorted(res.t, 0.5e-6)
    y = np.concatenate([[res.n[s][k] / m.N0 for s in m.species],
                        [res.ue[k] / m.N0, res.tg[k]],
                        np.zeros(m.size - m.ns - 2)])
    dy = m.rhs(res.t[k], y, True)
    p_in = m.p.power_W_cm3 / QE / m.N0
    assert abs(dy[m.i_ue]) / p_in < 1e-3


def test_rtol_converged():
    a = row_dict(run(**FAST, rtol=1e-7))
    b = row_dict(run(**FAST, rtol=1e-10))
    for key in ("to_Kr2star", "emission_vuv", "emission_nir",
                "heat_reactions", "elastic", "to_Krplus"):
        assert a[key] == approx(b[key], rel=1e-4), key


def test_te_scales_the_right_way_with_ne():
    """More seed electrons share the same power early in the pulse: the
    mean Te falls."""
    te = [np.mean(run(**FAST, ne_cm3=n).te) for n in (1e14, 1e15, 1e16)]
    assert te[0] > te[1] > te[2]


def test_summary_and_cli(capsys):
    from v0_2.run import main
    assert main(["--pulse", "0.2", "--t-end", "0.3", "--eedf", "two_t"]) == 0
    out = capsys.readouterr().out
    assert "LIGHT_Kr v0_2: 0-D pure-krypton model" in out
    assert "FATE AT t_end" in out and "closure residual" in out
    assert "no measurement, 0: hm_kr, kr2_pen, kstep_p, kstep_s" in out
    assert "electrons: two_t EEDF, captured energy exact" in out
    assert "Kr2* fate" in out and "KrF" not in out
    assert main(["--pulse", "0.2", "--t-end", "0.3", "--eedf", "two_t",
                 "--rate", "kr2_pen=3e-10"]) == 0
    out = capsys.readouterr().out
    assert ("bands: kr2_pen 3e-10;  no measurement, 0: hm_kr, kstep_p, "
            "kstep_s") in out
    with pytest.raises(SystemExit):
        main(["--rate", "kr2_f=3e-10"])          # a fluorine band
    with pytest.raises(SystemExit):
        main(["--f2-pct", "0.5"])                # no F2 to set
    text = summary(run(**FAST, groups=("electron",)))
    assert "optional groups: electron" in text and "87 rows" in text


def test_tg_rise_matches_monatomic_heat_capacity_for_pure_heat():
    """All deposited energy that leaves the electrons is heat, and the gas
    is monatomic but for its small Kr2* and Kr2+ populations, so
    1.5 N k dTg = heat."""
    res = run(**FAST)
    heat = _heat_booked(res)[-1]
    n_tot = sum(res.n[s][-1] for s in res.model.species)
    d_tg = heat / (1.5 * n_tot * KB_EV)
    assert res.tg[-1] - res.tg[0] == approx(d_tg, rel=2e-3)
    assert math.isfinite(res.tg[-1])


def test_thermal_electrons_stay_thermal_without_drive():
    """Te = Tg and no drive: the elastic term (1 - Tg/Te) must transfer
    nothing, so Te stays at Tg (a missing gas-temperature factor cools the
    electrons below the gas). Every row is taken out, so the elastic
    exchange alone acts (LIGHT used its v0 network, which had no
    recombination)."""
    m = Model(Params(power_W_cm3=0.0, pulse_s=0.0, t_end_s=200e-9))
    m.compiled = []
    res = m.run()
    tg_ev = res.tg * KB_EV
    assert np.max(np.abs(res.te / tg_ev - 1.0)) < 1e-6


def _recombination_only(**kw):
    """Kr2+ + e alone (no drive), the seed put in Kr2+ and the electrons
    started hot, at 20 x 300 K = 0.52 eV. k = 1.6e-6 (300 K/Te)^0.55 = 3.1e-7
    cm3/s on n_e = 1e15: a 3.2-ns recombination time, so a 3-ns window
    takes n_e to about half."""
    m = Model(Params(power_W_cm3=0.0, pulse_s=0.0, t_end_s=3e-9,
                     tg0_K=300.0, **kw))
    m.compiled = [c for c in m.compiled if c.row.rid == "R11p"]
    y0 = m.initial_state()
    seed = y0[m.idx["Kr+"]]
    y0[m.idx["Kr+"]] = 0.0
    y0[m.idx["Kr2+"]] = seed
    y0[m.idx["Kr"]] -= seed
    y0[m.i_ue] *= 20.0
    m.initial_state = lambda: y0
    return m, y0


def test_three_halves_capture_keeps_te_while_n_e_decays():
    """The notes' closure: each captured electron takes (3/2)Te with it, so
    Te stays put while n_e decays (elastic still relaxes it toward Tg).
    Without the capture term the energy per electron would RISE. (LIGHT
    held it for F2 attachment; here it is krypton's recombination.)"""
    m, _ = _recombination_only(capture="three_halves")
    res = m.run()
    assert 0.2 * res.ne[0] < res.ne[-1] < 0.8 * res.ne[0]
    assert np.all(np.diff(res.te) <= 1e-7 * res.te[0])
    assert res.te[-1] < res.te[0]


def test_exact_capture_books_the_recombined_electrons_energy():
    """The default closure: a recombining electron takes (3/2 + n) Te =
    0.95 Te with Shiu & Biondi's Te^-0.55, less than the mean (3/2) Te, so
    the rest heat up while n_e decays (recombination heating), and
    dU_e/dt is exactly -0.95 Te k [Kr2+] n_e - P_el. What the electrons
    lose is the booked capture and elastic heat."""
    m, y0 = _recombination_only()
    assert m.capture == "exact"
    n = y0[:m.ns] * m.N0
    ne = m.electron_density(n)
    te = y0[m.i_ue] * m.N0 / (1.5 * ne)
    c = next(c for c in m.compiled if c.row.rid == "R11p")
    k = 1.6e-6 * (te / (300.0 * 8.617333262e-5)) ** -0.55
    assert m._k(c, te, 300.0, n, ne) == approx(k, rel=1e-12)
    assert c.cap(te) == approx(0.95 * te, rel=1e-12)
    p_el = ne * n[m.idx["Kr"]] * float(m.kel(te)) * (1 - 300.0 * KB_EV / te)
    dy = m.rhs(0.0, y0, False)
    assert dy[m.i_ue] * m.N0 == approx(
        -0.95 * te * k * n[m.idx["Kr2+"]] * ne - p_el, rel=1e-9)
    res = m.run()
    assert 0.2 * res.ne[0] < res.ne[-1] < 0.8 * res.ne[0]
    assert res.te[-1] > 1.05 * res.te[0]
    lost = res.ledger["heat_capture"][-1] + res.ledger["heat_elastic"][-1]
    assert res.ue[0] - res.ue[-1] == approx(lost, rel=1e-6)


def test_closures_auto_two_t_and_the_pure_kr_tail():
    """'auto' is the kinetic EEDF with the exact captured energy. LIGHT
    v0_2's closure stays as eedf='two_t': its Kr electron-impact rows (since
    LIGHT_Kr v0_2 every level, its Hornbeck-Molnar twins and the
    ionization, 72 rows) read the tail table, at its F2 = 0 column (pure
    Kr), and at a laser plasma's n_e/N the tail row falls well below the
    Maxwellian."""
    from v0_2 import xsec
    auto = Model(_Params(groups=("all",)))
    assert (auto.eedf, auto.capture) == ("kinetic", "exact")
    b = Model(Params(groups=("all",)))
    assert (b.eedf, b.capture) == ("two_t", "exact")
    from v0_2.network import group_rows
    level_ids = {r.rid for g in ("G1", "C1a", "C1b", "C1c", "HM")
                 for r in group_rows(b.rows, g)}
    assert {c.row.rid for c in b.compiled if c.kind == "xsec_tail"} \
        == level_ids | {"C2"}
    assert len(level_ids) == 71
    cb = {c.row.rid: c for c in b.compiled}
    assert cb["R11p"].cap(1.3) == approx(0.95 * 1.3, rel=1e-12)  # Te^-0.55
    n = np.zeros(b.ns)
    n[b.idx["Kr"]] = 2.0e19
    te, ne = 1.6, 2.0e19 * 2.9e-5
    tail = xsec.tail_ratio_table()
    r = tail(2.9e-5, 0.0)
    assert r != tail(2.9e-5, 0.005)               # the column matters
    g1 = ("Kr", "excitation", "Kr(1S5)")
    assert b._k(cb["G1_1s5"], te, 300.0, n, ne) == approx(
        xsec.two_t_table(*g1)(te, r * te), rel=1e-12)
    assert b._k(cb["G1_1s5"], te, 300.0, n, ne) < \
        0.3 * xsec.rate_table(*g1)(te)
    # a Hornbeck-Molnar row reads its level's table times hm_kr (0 here)
    assert b._k(cb["HM_high"], te, 300.0, n, ne) == 0.0
    bh = Model(Params(rates=(("hm_kr", 0.25),)))
    ch = {c.row.rid: c for c in bh.compiled}
    hi = xsec.two_t_table("Kr", "excitation", "Kr(HIGH)")(te, r * te)
    assert bh._k(ch["HM_high"], te, 300.0, n, ne) == approx(0.25 * hi,
                                                            rel=1e-12)
    assert bh._k(ch["C1c_high"], te, 300.0, n, ne) == approx(0.75 * hi,
                                                             rel=1e-12)
    with pytest.raises(ValueError):
        Model(Params(eedf="boltzmann"))
    with pytest.raises(ValueError):
        Model(Params(capture="half"))
    with pytest.raises(ValueError):             # the kinetic EEDF books its
        Model(_Params(capture="three_halves"))  # own captured energies


def test_rate_laws_reproduce_the_printed_values():
    """Engine-side evaluation of the non-constant rate laws, each against a
    number printed in its own citation."""
    m = Model(Params(groups=("all",)))
    c = {x.row.rid: x for x in m.compiled}
    n = np.full(m.ns, 0.0)
    n[m.idx["Kr"]] = m.N0
    # KSTEP: 'the template evaluates 1.4e-8 at Te = 1.5 eV' (REBUILD cite)
    assert m._k(c["KSTEP"], 1.5, 300.0, n, 1e15) == approx(1.4e-8,
                                                           rel=0.03)
    # R11p: Shiu & Biondi, 1.6e-6 (300/Te)^0.55 cm3/s; Mikus 1978 measured
    # 1.48e-6 (Te/300 K)^-0.53 to 21300 K, within 8 %
    te_300 = 300.0 * KB_EV
    assert m._k(c["R11p"], te_300, 300.0, n, 1e15) == approx(1.6e-6,
                                                             rel=1e-6)
    for te in (0.3, 1.0, 1.8):
        mikus = 1.48e-6 * (te / te_300) ** -0.53
        assert m._k(c["R11p"], te, 300.0, n, 1e15) == approx(mikus, rel=0.08)
    # K15: the mean of three 300 K measurements, on Perez-Rios & Greene's
    # T^-3/4
    assert m._k(c["K15"], 1.0, 300.0, n, 1e15) == approx(2.4e-31, rel=1e-12)
    assert m._k(c["K15"], 1.0, 600.0, n, 1e15) == approx(2.4e-31 * 2 ** -0.75,
                                                         rel=1e-12)
    # AK8_Kr: Ueda's 1.7e-32 on Wieme & Lenaerts' (Tg/300)^-0.64
    assert m._k(c["AK8_Kr"], 1.0, 600.0, n, 1e15) == approx(
        1.7e-32 * 2.0 ** -0.64, rel=1e-12)
    # LIGHT v0_5: K76 is the Maxwellian table of Jung et al.'s lumped cross
    # section and K76r its Klein-Rosseland reverse (detailed balance over
    # the ledger's gap, to the table's quadrature)
    from v0_2 import jung2006
    assert m._k(c["K76"], 1.3, 300.0, n, 1e15) == approx(
        float(jung2006.rate_table()(1.3)), rel=1e-12)
    ratio = (m._k(c["K76r"], 1.3, 300.0, n, 1e15)
             / m._k(c["K76"], 1.3, 300.0, n, 1e15))
    gap = 11.3035 - 9.9152
    assert ratio == approx(12 / 36 * math.exp(gap / 1.3), rel=2e-2)
    # the band scales both; Params.k76 = 'kannari' is LIGHT v0_4's fit
    m5 = Model(Params(rates=(("k_5s5p", 0.5),)))
    c5 = {x.row.rid: x for x in m5.compiled}
    assert m5._k(c5["K76"], 1.3, 300.0, n, 1e15) == approx(
        0.5 * m._k(c["K76"], 1.3, 300.0, n, 1e15), rel=1e-12)
    mj = Model(Params(k76="kannari"))
    cj = {x.row.rid: x for x in mj.compiled}
    assert mj._k(cj["K76"], 1.3, 300.0, n, 1e15) == approx(
        9.0e-7 * 1.3 ** 0.5 * math.exp(-1.51 / 1.3), rel=1e-12)
    ratio = (mj._k(cj["K76r"], 1.3, 300.0, n, 1e15)
             / mj._k(cj["K76"], 1.3, 300.0, n, 1e15))
    assert ratio == approx(12 / 36 * math.exp(gap / 1.3), rel=1e-9)
    assert m._k(c["KSTEP_S"], 1.3, 300.0, n, 1e15) == 0.0     # band at 0
    mk = Model(Params(ladder="kannari"))
    ck = {x.row.rid: x for x in mk.compiled}
    assert mk._k(ck["KSTEP_S"], 1.3, 300.0, n, 1e15) == approx(
        7.0e-8 * 1.3 ** 0.68 * math.exp(-4.04 / 1.3), rel=1e-12)
    assert mk._k(ck["KSTEP_P"], 1.3, 300.0, n, 1e15) == approx(
        2.0e-7 * 1.3 ** 0.62 * math.exp(-2.53 / 1.3), rel=1e-12)
    assert mk._k(ck["KPS"], 1.3, 300.0, n, 1e15) == 1.0e-10


# ---------------------------------------------------------------------------
# the reference point, pinned (2026-10-14): LIGHT's drive (300 kW/cm^3 for
# 1 us at 1 atm, a 1e15 cm^-3 seed, 2 us) in pure Kr, under LIGHT v0_2's
# closure, on Biagi 8.97's levels and jung2006's tables as corrected to the
# paper. A regression pin, not a fit. v0_1 on Biagi-v7.1's groups gave
# 15.051, 0.6340, 802.38 K, 2.2555 eV, 81.812 and 3.0270
V0_2_REF = dict(emission_vuv=15.9976, emission_nir=0.62578, tg_end_K=796.76,
                te_mean_eV=2.0036, heat_reactions=81.398, elastic=2.4915)


def test_reference_point_regression():
    from v0_2.report import warnings
    res = run()
    d = row_dict(res)
    for k, v in V0_2_REF.items():
        assert d[k] == approx(v, rel=2e-3), k
    assert abs(d["closure"]) < 1e-12
    assert any("Kr* + 2Kr's (Tg/300)^-0.64" in w for w in warnings(res))
    # the seed only sets the first nanoseconds
    seed = row_dict(run(ne_cm3=1e12))
    assert seed["emission_vuv"] == approx(d["emission_vuv"], rel=0.02)


# ---------------------------------------------------------------------------
# LIGHT v0_2: the bands
# ---------------------------------------------------------------------------
def test_param_rows_read_their_band():
    """A ("param", name[, frac]) row takes frac x Params.rates[name], or the
    band's default. A band no row reads, a fluorine band, or a negative
    rate, is refused."""
    from v0_2.network import BANDS
    c = {x.row.rid: x for x in Model(Params()).compiled}
    for rid, name in (("G2", "kr_pen"), ("AK8_Kr", "kr2_form"),
                      ("PEN_KR2", "kr2_pen"), ("KPS", "kps_kr")):
        assert c[rid].kind == "param" and c[rid].cap is None, rid
        assert c[rid].rate == approx(BANDS[name].default, rel=1e-12), rid
    assert c["PEN_KR2"].rate == 0.0
    m = Model(Params(rates=(("kr2_pen", 3e-10), ("kr2_form", 4.4e-32))))
    c = {x.row.rid: x for x in m.compiled}
    assert c["PEN_KR2"].rate == approx(3e-10, rel=1e-12)
    n = np.zeros(m.ns)
    n[m.idx["Kr"]] = m.N0
    assert m._k(c["AK8_Kr"], 1.0, 600.0, n, 1e15) == approx(
        4.4e-32 * 2.0 ** -0.64, rel=1e-12)
    for bad in ((("k2f", 1e-10),), (("kr2_pen", -1e-10),),
                (("kr2_f", 3e-10),), (("trimer_tg_pow", -0.5),)):
        with pytest.raises(ValueError):
            Model(Params(rates=bad))


def test_bands_run_each_end_once():
    """--bands runs the network as given, then each band's two ends one at a
    time. An end at the band's current value is that run again."""
    from v0_2.network import BANDS
    from v0_2.run import band_table, do_bands, with_band
    base = Params(pulse_s=0.1e-6, t_end_s=0.2e-6, n_save=101)
    assert with_band(base.replace(rates=(("kr2_pen", 1e-10),)), "kr2_pen",
                     2e-10).rates == (("kr2_pen", 2e-10),)
    rows = do_bands(base)
    assert len(rows) == 1 + 2 * len(BANDS) and rows[0]["band"] == "as given"
    ref = rows[0]
    by = {(r["band"], r["end"]): r for r in rows[1:]}
    for b in BANDS.values():
        for end, v in (("low", b.lo), ("high", b.hi)):
            assert by[(b.name, end)]["value"] == v
            if v == b.default:
                assert by[(b.name, end)]["emission_vuv"] == ref["emission_vuv"]
    # slower 5p relaxation leaves more of the 5p rung to radiate
    assert by[("kps_kr", "low")]["emission_nir"] > 2 * ref["emission_nir"]
    assert "Kr2* + Kr2* -> Kr2+ + 2Kr + e" in band_table(rows)


def test_kr2_warning_needs_the_route_and_hot_gas():
    """The Kr* + 2Kr temperature law is measured to 451 K. The report warns
    above it, and not when the gas stays cool or the route is off."""
    from v0_2.report import warnings
    cool = run(**FAST, power_W_cm3=3e4)
    assert cool.tg.max() < 451.0
    assert not any("Kr* + 2Kr" in w for w in warnings(cool))
    m = Model(Params(**FAST, power_W_cm3=1.5e6, rates=(("kr2_form", 0.0),)))
    hot = m.run()
    assert hot.tg.max() > 451.0
    assert not any("Kr* + 2Kr" in w for w in warnings(hot))
