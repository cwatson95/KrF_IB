"""The 0-D ArF model: ledger closures, analytic limits, convergence, CLI."""

import math

import numpy as np
import pytest

from v0_2 import network as N
from v0_2.model import AMAGAT_CM3, KB_EV, QE, Model, Params, run
from v0_2.report import accounting, reservoir, row_dict, summary


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
            + res.ledger["heat_beam"]
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
    for name in ("Ar*", "ArF*", "Ar2F*", "Ar2*"):
        rv = reservoir(res, name)
        assert abs(rv["balance"]) <= 1e-6 * max(rv["formed"], 1.0), name


def test_reference_point_closes():
    res = run()
    _check_closures(res)
    assert res.ledger["E_dep"][-1] * QE == approx(0.3, rel=1e-9)


@pytest.mark.parametrize("kw", [
    dict(groups=("all",)),
    dict(eedf="maxwellian", capture="three_halves"),
    dict(pressure_atm=0.5, power_W_cm3=5e5, ne_cm3=1e14),
    dict(pressure_atm=3.0, power_W_cm3=1e5, ne_cm3=1e16),
    dict(rates=(("ar2_f", 1e-9), ("ar2_pen", 5e-10), ("ar_pen", 6.2e-10),
                ("ar2f_ar", 2.2e-14), ("arf_yield", 1.0),
                ("ar2_rad", 6.0e7)), trimer_tg_pow=-3.9),
])
def test_closures_across_modes(kw):
    _check_closures(run(**FAST, **kw))


def test_free_ne_is_quasineutral_and_clamped_is_refused():
    res = run(**FAST)
    assert res.model.ne_mode == "free"
    assert res.ne[0] == approx(1e15, rel=1e-12)          # the seed
    cat = sum(np.maximum(res.n[s], 0.0) for s in ("Ar+", "Ar2+"))
    assert np.allclose(res.ne, np.maximum(cat - np.maximum(res.n["F-"], 0.0),
                                          1.0), rtol=1e-12)
    with pytest.raises(ValueError):
        Model(Params(ne_mode="clamped"))


def _arf_seed(**kw):
    """No power, no electrons: an ArF(B) seed in Ar with a little F2."""
    p = Params(power_W_cm3=0.0, init_extra_cm3=(("ArF(B)", 1e10),),
               pulse_s=0.0, t_end_s=200e-9, ne_cm3=1.0, **kw)
    return Model(p).run()


@pytest.mark.parametrize("tau", [3.0e-9, 4.2e-9, 6.0e-9])
def test_arf_branching_is_the_measured_stern_volmer(tau):
    """The emitted share of an ArF(B) seed is 1 / (1 + (k tau)_F2 [F2] +
    (k tau)_Ar [Ar] + (k tau)_3 [Ar]^2) with the measured products,
    whatever the calculated tau_B: the rows carry k tau (network.py)."""
    res = _arf_seed(tau_sp_s=tau, pressure_atm=0.5, x_f2=0.002)
    ar, f2 = res.n["Ar"][0], res.n["F2"][0]
    kt3 = N.BANDS["arf_trimer"].default
    expect = 1.0 / (1.0 + N.KT_ARF_F2 * f2 + N.KT_ARF_AR * ar
                    + kt3 * ar * ar)
    assert res.events["G4"][-1] / 1e10 == approx(expect, rel=1e-5)
    assert res.events["G5"][-1] / 1e10 == approx(kt3 * ar * ar * expect,
                                                 rel=1e-5)


def test_harpoon_takes_a_dilute_ar_star_seed():
    """With Ar* + 2Ar off, every Ar* reacts with F2: half of it makes
    ArF(B) (Setser's yield), half dissociates (G3X). With the Ar2* route
    on, Ar* + 2Ar takes its analytic share at 300 K."""
    p = Params(power_W_cm3=0.0, init_extra_cm3=(("Ar*", 1e10),),
               pulse_s=0.0, t_end_s=200e-9, ne_cm3=1.0)   # 20 lifetimes
    ev = Model(p.replace(rates=(("ar2_form", 0.0),))).run().events
    assert (ev["G3"][-1] + ev["G3X"][-1]) / 1e10 == approx(1.0, rel=1e-5)
    assert ev["G3"][-1] / 1e10 == approx(0.5, rel=1e-5)
    res = Model(p).run()
    harpoon = N.BANDS["ar_f2"].default * res.n["F2"][0]
    three = 6.0e-33 * res.n["Ar"][0] ** 2
    share = harpoon / (harpoon + three)
    ev = res.events
    assert (ev["G3"][-1] + ev["G3X"][-1]) / 1e10 == approx(share, rel=1e-5)
    assert ev["AK8_Ar"][-1] / 1e10 == approx(1.0 - share, rel=1e-4)


def test_tau_sp_scales_every_stern_volmer_row():
    """Changing the calculated lifetime keeps each measured k tau product."""
    for tau in (3e-9, 4.2e-9, 6e-9):
        c = {x.row.rid: x for x in Model(Params(tau_sp_s=tau)).compiled}
        assert c["G4"].rate * tau == approx(1.0, rel=1e-12)
        assert c["QF2"].rate * tau == approx(N.KT_ARF_F2, rel=1e-12)
        assert c["QAR"].rate * tau == approx(N.KT_ARF_AR, rel=1e-12)
        assert c["G5"].rate * tau == approx(N.BANDS["arf_trimer"].default,
                                            rel=1e-12)
        assert c["QE"].rate * tau == approx(N.BANDS["arf_e"].default,
                                            rel=1e-12)


def test_rtol_converged():
    a = row_dict(run(**FAST, rtol=1e-7))
    b = row_dict(run(**FAST, rtol=1e-10))
    for key in ("to_ArFB", "emission_193", "heat_reactions", "elastic",
                "to_Arplus"):
        assert a[key] == approx(b[key], rel=1e-4), key


def test_summary_and_cli(capsys):
    from v0_2.run import main
    assert main(["--pulse", "0.2", "--t-end", "0.3"]) == 0
    out = capsys.readouterr().out
    assert "FATE AT t_end" in out and "closure residual" in out
    assert "no measurement, 0: ar2_f, ar2_pen, ar2f_ar, ar_pen" in out
    assert main(["--pulse", "0.2", "--t-end", "0.3", "--rate",
                 "ar2_f=3e-10"]) == 0
    out = capsys.readouterr().out
    assert "bands: ar2_f 3e-10;  no measurement, 0: ar2_pen, ar2f_ar, " \
        "ar_pen" in out
    for bad in (["--rate", "a2f=3e-10"], ["--rate", "trimer_tg_pow=-1"]):
        with pytest.raises(SystemExit):
            main(bad)


def test_tg_rise_matches_monatomic_heat_capacity_for_pure_heat():
    """A pure-Ar fill: the energy that leaves the electrons as heat raises
    Tg by heat / (1.5 N k) (the composition barely changes)."""
    res = run(**FAST, x_f2=0.0)
    heat = _heat_booked(res)[-1]
    n_tot = sum(res.n[s][-1] for s in res.model.species)
    d_tg = heat / (1.5 * n_tot * KB_EV)
    assert res.tg[-1] - res.tg[0] == approx(d_tg, rel=5e-3)
    assert math.isfinite(res.tg[-1])


def test_thermal_electrons_stay_thermal_without_drive():
    """Te = Tg and no drive: the elastic term (1 - Tg/Te) must transfer
    nothing (a missing gas-temperature factor cools the electrons below the
    gas). Every reaction row is taken out."""
    m = Model(Params(power_W_cm3=0.0, x_f2=0.0, pulse_s=0.0,
                     t_end_s=200e-9))
    m.compiled = []
    res = m.run()
    assert np.max(np.abs(res.te / (res.tg * KB_EV) - 1.0)) < 1e-6


def _attachment_only(**kw):
    """Attachment alone (free n_e, no drive), electrons started hot at
    0.5 eV: a ns-scale n_e decay that does not reach the floor."""
    m = Model(Params(power_W_cm3=0.0, pulse_s=0.0, t_end_s=1.5e-9,
                     tg0_K=300.0, **kw))
    m.compiled = [c for c in m.compiled if c.row.rid == "ATT"]
    y0 = m.initial_state()
    y0[m.i_ue] *= 20.0
    m.initial_state = lambda: y0
    return m, y0


def test_three_halves_capture_keeps_te_while_n_e_decays():
    m, _ = _attachment_only(capture="three_halves")
    res = m.run()
    assert 1e-3 * res.ne[0] < res.ne[-1] < 0.5 * res.ne[0]
    assert np.all(np.diff(res.te) <= 1e-7 * res.te[0])


def test_exact_capture_books_the_attached_electrons_energy():
    """dU_e/dt = -E_cap(Te) k_att [F2] n_e - P_el exactly, and since F2
    attaches the slow electrons the rest heat up."""
    from v0_2 import xsec
    m, y0 = _attachment_only()
    assert m.capture == "exact"
    n = y0[:m.ns] * m.N0
    ne = m.electron_density(n)
    te = y0[m.i_ue] * m.N0 / (1.5 * ne)
    k = float(xsec.rate_table("F2", "attachment")(te))
    e_cap = xsec.capture_energy_table("F2", "attachment")(te)
    p_el = ne * n[m.idx["Ar"]] * float(m.kel(te)) * (1 - 300.0 * KB_EV / te)
    dy = m.rhs(0.0, y0, False)
    assert dy[m.i_ue] * m.N0 == approx(
        -e_cap * k * n[m.idx["F2"]] * ne - p_el, rel=1e-9)
    res = m.run()
    assert res.te[-1] > 1.05 * res.te[0]
    lost = res.ledger["heat_capture"][-1] + res.ledger["heat_elastic"][-1]
    assert res.ue[0] - res.ue[-1] == approx(lost, rel=1e-6)


def test_basic_network_runs_the_two_temperature_closure():
    """'auto' gives the network the two-temperature tail on all 38 Ar rows
    and the exact captured energy; a tail row reads the kinetic table's r."""
    from v0_2 import xsec
    b = Model(Params())
    assert (b.eedf, b.capture) == ("two_t", "exact")
    tail = sorted(c.row.rid for c in b.compiled if c.kind == "xsec_tail")
    assert len(tail) == 38 and "C2" in tail and "G1_3P2" in tail
    cb = {c.row.rid: c for c in b.compiled}
    assert cb["R11p"].cap(1.3) == approx((1.5 - 0.67) * 1.3, rel=1e-12)
    n = np.zeros(b.ns)
    n[b.idx["Ar"]], n[b.idx["F2"]] = 2.0e19, 1.0e17
    te, ne = 1.6, 2.01e19 * 2.0e-5
    r = xsec.tail_ratio_table()(2.0e-5, 1.0e17 / 2.01e19)
    g1 = ("Ar", "excitation", "Ar(3P2)")
    assert b._k(cb["G1_3P2"], te, 300.0, n, ne) == approx(
        xsec.two_t_table(*g1)(te, r * te), rel=1e-12)
    m = Model(Params(eedf="maxwellian"))
    assert not any(c.kind == "xsec_tail" for c in m.compiled)
    with pytest.raises(ValueError):
        Model(Params(eedf="kinetic"))


def test_rate_laws_reproduce_the_printed_values():
    """Engine-side evaluation of the non-constant rate laws, each against
    the form its citation prints."""
    m = Model(Params())
    c = {x.row.rid: x for x in m.compiled}
    n = np.zeros(m.ns)
    n[m.idx["Ar"]] = AMAGAT_CM3                     # 1 amagat
    te300 = 300.0 * KB_EV
    # R11p: 7.5e-7 (Te/300 K)^-0.67
    assert m._k(c["R11p"], te300, 300.0, n, 1e15) == approx(7.5e-7)
    assert m._k(c["R11p"], 1.0, 300.0, n, 1e15) == approx(
        7.5e-7 * (1.0 / te300) ** -0.67)
    # K15 on the T^-3/4 theory
    assert m._k(c["K15"], 1.0, 600.0, n, 1e15) == approx(2.2e-31 * 2 ** -0.75)
    # K26 and K108 at 1 amagat, 300 K: Lee's printed Hoffman & Moreno forms
    # times the Debye screening; at 600 K the three-body limit falls x8
    s300 = m._screening(1e15, n, 300.0)
    assert 0.0 < s300 < 1.0
    assert m._k(c["K26"], 1.0, 300.0, n, 1e15) == approx(
        9.32e-6 / (1 + 6.21e-2 + 1.69) * s300, rel=1e-9)
    assert m._k(c["K108"], 1.0, 300.0, n, 1e15) == approx(
        1.51e-5 / (1 + 1.81 + 2.25) * s300, rel=1e-9)
    s600 = m._screening(1e15, n, 600.0)
    assert m._k(c["K26"], 1.0, 600.0, n, 1e15) == approx(
        9.32e-6 / (1 + 6.21e-2 + 8 * 1.69) * s600, rel=1e-9)
    # the trimer: no Tg law unless a band run sets one
    assert m._k(c["G5"], 1.0, 700.0, n, 1e15) == approx(
        m._k(c["G5"], 1.0, 300.0, n, 1e15), rel=1e-12)
    hot = Model(Params(trimer_tg_pow=-3.9))
    ch = {x.row.rid: x for x in hot.compiled}
    assert hot._k(ch["G5"], 1.0, 700.0, n, 1e15) == approx(
        m._k(c["G5"], 1.0, 300.0, n, 1e15) * (700 / 300) ** -3.9, rel=1e-9)
    # F + F + Ar: the 295 K measurement, held
    assert m._k(c["FF"], 1.0, 1300.0, n, 1e15) == approx(8.0e-35)


def test_fluorine_recombination_is_too_slow_to_refill_f2():
    res = run()
    assert 0.0 < res.events["FF"][-1] < 0.01 * res.n["F2"][0]


def test_trimer_switch_moves_only_the_trimer():
    base = Model(Params())
    alt = Model(Params(trimer_tg_pow=-3.9))
    n = np.zeros(base.ns)
    n[base.idx["Ar"]] = base.N0
    cb = {c.row.rid: c for c in base.compiled}
    ca = {c.row.rid: c for c in alt.compiled}
    for rid in cb:
        if cb[rid].kind == "elastic":
            continue
        kb = base._k(cb[rid], 1.5, 700.0, n, 1e14)
        ka = alt._k(ca[rid], 1.5, 700.0, n, 1e14)
        if rid == "G5":
            assert ka / kb == approx((700 / 300) ** -3.9, rel=1e-9)
        else:
            assert ka == approx(kb, rel=1e-12), rid


def test_param_rows_read_their_band():
    """A band row takes Params.rates or the band's default; a band no row
    reads, a negative rate or a fraction above 1 is refused."""
    B = N.BANDS
    c = {x.row.rid: x for x in Model(Params()).compiled}
    assert c["G3"].rate == approx(B["ar_f2"].default * 0.5, rel=1e-12)
    assert c["NK7"].rate + c["K2F_BR"].rate == approx(B["ar2_f2"].default,
                                                      rel=1e-12)
    assert c["G2"].rate == c["K2F_F"].rate == c["PEN_AR2"].rate == 0.0
    assert c["Q2AR"].rate == 0.0 and c["AR2_RAD"].rate == approx(1 / 3.2e-6)
    m = Model(Params(rates=(("arf_yield", 1.0), ("ar2_f", 3e-10))))
    c = {x.row.rid: x for x in m.compiled}
    assert c["G3X"].rate == 0.0
    assert c["G3"].rate == approx(B["ar_f2"].default, rel=1e-12)
    assert c["K2F_F"].rate == approx(3e-10, rel=1e-12)
    for bad in ((("a2f", 1e-10),), (("ar2_f", -1e-10),),
                (("trimer_tg_pow", -0.5),), (("arf_yield", 1.5),)):
        with pytest.raises(ValueError):
            Model(Params(rates=bad))


def test_bands_run_each_end_once():
    from v0_2.run import band_table, do_bands, with_band
    base = Params(pulse_s=0.1e-6, t_end_s=0.2e-6, n_save=101)
    assert with_band(base, "trimer_tg_pow", -3.9).trimer_tg_pow == -3.9
    rows = do_bands(base)
    assert len(rows) == 1 + 2 * len(N.BANDS) and rows[0]["band"] == "as given"
    ref = rows[0]["emission_193"]
    by = {(r["band"], r["end"]): r for r in rows[1:]}
    for b in N.BANDS.values():
        for end, v in (("low", b.lo), ("high", b.hi)):
            assert by[(b.name, end)]["value"] == v
            if v == b.default:
                assert by[(b.name, end)]["emission_193"] == ref
    assert by[("arf_yield", "high")]["emission_193"] > ref
    assert "Ar2* + F -> ArF(B) + Ar" in band_table(rows)


def test_warnings_name_the_unmeasured_tg_laws():
    from v0_2.report import warnings
    res = run()
    assert res.tg.max() > 400.0
    w = warnings(res)
    assert any("ArF* + 2Ar is measured near 300 K" in x for x in w)
    assert any("F2 burn-out" in x for x in w)
    cool = run(**FAST, power_W_cm3=2e4)
    assert cool.tg.max() < 400.0
    assert not any("measured near 300 K" in x for x in warnings(cool))


def test_reference_point_regression():
    """The reference point (1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, 2 us) on
    the defaults (2026-09-30), and the three bands that move 193 nm most:
    Ar2* + F at its top (the Ar2* left when F2 burns out makes ArF), the
    ArF* yield at 1.0 (every published model's assumption) and NRL's
    Tg^-3.9 trimer law. The seed only sets the first nanoseconds."""
    d = row_dict(run())
    assert d["emission_193"] == approx(5.052, rel=2e-3)
    assert d["emission_vuv"] == approx(4.604, rel=2e-3)
    assert d["to_ArFB"] == approx(18.69, rel=2e-3)
    assert d["tg_end_K"] == approx(675.2, rel=2e-3)
    assert abs(d["closure"]) < 1e-12
    assert row_dict(run(rates=(("ar2_f", 1e-9),)))["emission_193"] == \
        approx(8.558, rel=2e-3)
    assert row_dict(run(rates=(("arf_yield", 1.0),)))["emission_193"] == \
        approx(8.893, rel=2e-3)
    assert row_dict(run(trimer_tg_pow=-3.9))["emission_193"] == \
        approx(7.010, rel=2e-3)
    assert row_dict(run(ne_cm3=1e13))["emission_193"] == approx(
        d["emission_193"], rel=0.01)
