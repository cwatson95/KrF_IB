"""The 0-D model: ledger closures, analytic limits, convergence, CLI.

These are v0_2's tests. v0_4's chemistry is v0_2's, so they run on v0_2's
electron closure (eedf='two_t', the Maxwellian with the tabulated tail) unless
a test names another: Params and run below default to it. The kinetic EEDF,
v0_4's own default, has its tests in test_kinetic.py."""

import math
from pathlib import Path

import numpy as np
import pytest

from v0_4.model import KB_EV, QE, Model
from v0_4.model import Params as _Params
from v0_4.model import run as _run
from v0_4.report import accounting, reservoir, row_dict, summary


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
    for name in ("Kr*", "KrF*"):
        rv = reservoir(res, name)
        assert abs(rv["balance"]) <= 1e-6 * max(rv["formed"], 1.0)


def test_reference_point_closes():
    res = run()
    _check_closures(res)
    e_dep_j = res.ledger["E_dep"][-1] * QE
    assert e_dep_j == approx(0.3, rel=1e-9)       # 300 kW/cm^3 x 1 us


@pytest.mark.parametrize("kw", [
    dict(groups=("all",)),
    dict(network="v0", groups=("all",)),
    dict(ne_mode="free"),
    dict(ne_mode="free", groups=("all",)),
    dict(pressure_atm=0.1, power_W_cm3=5e5, ne_cm3=1e14),
    dict(pressure_atm=2.0, power_W_cm3=1e5, ne_cm3=1e16),
])
def test_closures_across_modes(kw):
    _check_closures(run(**FAST, **kw))


def test_clamped_ne_is_constant_and_free_ne_is_quasineutral():
    res = run(**FAST, ne_mode="clamped")
    assert np.all(res.ne == 1e15)
    assert Model(Params()).ne_mode == "free"           # the basic default
    assert Model(Params(network="v0")).ne_mode == "clamped" 
    res = run(**FAST, ne_mode="free", groups=("all",))
    assert res.ne[0] == approx(1e15, rel=1e-12)   # the injected seed
    # the model's own rule: densities clipped at zero (the solver's noise
    # near zero is within its absolute tolerance, 1e-18 N0 = 25 cm^-3),
    # cations less anions, floored at 1 cm^-3
    cat = sum(np.maximum(res.n[s], 0.0) for s in ("Kr+", "Kr2+"))
    assert np.allclose(res.ne, np.maximum(cat - np.maximum(res.n["F-"], 0.0),
                                          1.0), rtol=1e-12)


def test_radiative_vs_trimer_branching_is_analytic():
    """No drive, a small KrF(B) seed: the emitted share must be
    A / (A + k5 [Kr]^2) at 300 K."""
    p = Params(network="v0", power_W_cm3=0.0,
               init_extra_cm3=(("KrF(B)", 1e10),),
               pulse_s=0.0, t_end_s=100e-9, x_f2=0.0)
    res = Model(p).run()
    kr = res.n["Kr"][0]
    a_rad = 1.43e8
    expect = a_rad / (a_rad + 6.0e-31 * kr * kr)
    got = res.events["G4"][-1] / 1e10
    assert got == approx(expect, rel=1e-6)


def test_harpoon_takes_a_dilute_Kr_star_seed():
    """With Kr* + 2Kr off, every Kr* is harpooned. In the basic network 90 %
    of the reactions make KrF* (Velazco, Kolts & Setser 1976: 65e-11 of
    72e-11), 62 % of that in KrF(B) and 38 % in KrF(C) (Yu & Setser's
    k_B/k_C 1.6), and 10 % dissociate (G3X). In v0 all go into B. With the
    Kr2* route on (v0_2), Kr* + 2Kr takes its analytic share at 300 K."""
    p = Params(power_W_cm3=0.0, init_extra_cm3=(("Kr*", 1e10),),
               pulse_s=0.0, t_end_s=200e-9, ne_cm3=1.0)   # no recombination
    ev = Model(p.replace(rates=(("kr2_form", 0.0),))).run().events
    assert (ev["G3"][-1] + ev["G3C"][-1] + ev["G3X"][-1]) / 1e10 == \
        approx(1.0, rel=1e-5)
    assert ev["G3C"][-1] / 1e10 == approx(0.90 * 0.38, rel=1e-5)
    assert ev["G3X"][-1] / 1e10 == approx(0.10, rel=1e-5)
    v0 = Model(p.replace(network="v0")).run().events
    assert v0["G3"][-1] / 1e10 == approx(1.0, rel=1e-5)
    res = Model(p).run()
    ev = res.events
    harpoon = 8.47e-10 * res.n["F2"][0]
    share = harpoon / (harpoon + 1.7e-32 * res.n["Kr"][0] ** 2)
    assert (ev["G3"][-1] + ev["G3C"][-1] + ev["G3X"][-1]) / 1e10 == \
        approx(share, rel=1e-5)
    assert ev["AK8_Kr"][-1] / 1e10 == approx(1.0 - share, rel=1e-4)


def test_electron_energy_is_quasi_steady_mid_pulse():
    """Te relaxes in ~ns; mid-pulse, at a held n_e, the heating must equal
    the losses."""
    m = Model(Params(ne_mode="clamped"))
    res = m.run()
    k = np.searchsorted(res.t, 0.5e-6)
    y = np.concatenate([[res.n[s][k] / m.N0 for s in m.species],
                        [res.ue[k] / m.N0, res.tg[k]],
                        np.zeros(m.size - m.ns - 2)])
    dy = m.rhs(res.t[k], y, True)
    p_in = m.p.power_W_cm3 / QE / m.N0
    assert abs(dy[m.i_ue]) / p_in < 1e-3


def test_tau_sp_override():
    """tau_sp is KrF(B)'s whole radiative lifetime: B-X (G4) and B-A
    (G4A, A(B-A)/A(B-X) = 0.035) share it."""
    m = Model(Params(tau_sp_s=4e-9))
    c = {x.row.rid: x for x in m.compiled}
    assert c["G4"].rate + c["G4A"].rate == approx(2.5e8, rel=1e-12)
    assert c["G4A"].rate / c["G4"].rate == approx(0.035, rel=1e-12)


def test_rtol_converged():
    a = row_dict(run(**FAST, rtol=1e-7))
    b = row_dict(run(**FAST, rtol=1e-10))
    for key in ("to_KrFB", "emission_248", "heat_reactions", "elastic",
                "to_Krplus"):
        assert a[key] == approx(b[key], rel=1e-4), key


def test_te_scales_the_right_way_with_ne():
    """More electrons share the same power: Te must fall."""
    te = [np.mean(run(**FAST, ne_cm3=n).te) for n in (1e14, 1e15, 1e16)]
    assert te[0] > te[1] > te[2]


def test_summary_and_cli(capsys):
    from v0_4.run import main
    assert main(["--pulse", "0.2", "--t-end", "0.3", "--eedf", "two_t"]) == 0
    out = capsys.readouterr().out
    assert "FATE AT t_end" in out and "closure residual" in out
    assert "no measurement, 0: kr2_f, kr2_pen, kr2f_f, krf_f" in out
    assert "electrons: two_t EEDF, captured energy exact" in out
    assert main(["--pulse", "0.2", "--t-end", "0.3", "--eedf", "two_t",
                 "--rate", "kr2_f=3e-10"]) == 0
    out = capsys.readouterr().out
    assert ("bands: kr2_f 3e-10;  no measurement, 0: kr2_pen, kr2f_f, krf_f"
            in out)
    with pytest.raises(SystemExit):
        main(["--rate", "k2f=3e-10"])
    text = summary(run(**FAST, network="v0", groups=("quench",)))
    assert "quench" in text and "bands:" not in text


def test_tg_rise_matches_monatomic_heat_capacity_for_pure_heat():
    """A pure-Kr fill, no reactions but elastic: all deposited energy that
    leaves the electrons is heat, and 1.5 N k dTg = heat."""
    res = run(**FAST, x_f2=0.0)
    heat = _heat_booked(res)[-1]
    n_tot = sum(res.n[s][-1] for s in res.model.species)
    d_tg = heat / (1.5 * n_tot * KB_EV)
    # composition barely changes (no F2), so the monatomic estimate holds
    assert res.tg[-1] - res.tg[0] == approx(d_tg, rel=2e-3)
    assert math.isfinite(res.tg[-1])


def test_thermal_electrons_stay_thermal_without_drive():
    """Te = Tg and no drive: the elastic term (1 - Tg/Te) must transfer
    nothing, so Te stays at Tg (a missing gas-temperature factor cools the
    electrons below the gas)."""
    res = run(power_W_cm3=0.0, x_f2=0.0, pulse_s=0.0, t_end_s=200e-9,
              network="v0")        # elastic only: no recombination there
    tg_ev = res.tg * KB_EV
    assert np.max(np.abs(res.te / tg_ev - 1.0)) < 1e-6


def _attachment_only(**kw):
    """Attachment alone (free n_e, no drive), electrons started hot at
    0.5 eV. k_att(0.5 eV) = 6.6e-9 cm3/s on 1.2e17 F2: a 1.25 ns lifetime,
    so a 1.5 ns window takes n_e to ~1/4 without reaching the floor."""
    m = Model(Params(power_W_cm3=0.0, pulse_s=0.0, t_end_s=1.5e-9,
                     ne_mode="free", tg0_K=300.0, **kw))
    m.compiled = [c for c in m.compiled if c.row.rid == "ATT"]
    y0 = m.initial_state()
    y0[m.i_ue] *= 20.0
    m.initial_state = lambda: y0
    return m, y0


def test_three_halves_capture_keeps_te_while_n_e_decays():
    """v0's closure: each captured electron takes (3/2)Te with it, so Te
    stays put while n_e decays (elastic still relaxes it toward Tg).
    Without the capture term the energy per electron would RISE."""
    m, _ = _attachment_only(capture="three_halves")
    res = m.run()
    assert 1e-3 * res.ne[0] < res.ne[-1] < 0.5 * res.ne[0]
    assert np.all(np.diff(res.te) <= 1e-7 * res.te[0])
    assert res.te[-1] < res.te[0]


def test_exact_capture_books_the_attached_electrons_energy():
    """The basic network's closure: an attached electron takes <eps sigma
    v>/<sigma v> from the electron gas. dU_e/dt is exactly
    -E_cap(Te) k_att [F2] n_e - P_el, and since F2 attaches the slow
    electrons the rest heat up (attachment heating): Te rises while n_e
    decays. What the electrons lose is the booked capture and elastic
    heat."""
    from v0_4 import xsec
    m, y0 = _attachment_only()
    assert m.capture == "exact"
    n = y0[:m.ns] * m.N0
    ne = m.electron_density(n)
    te = y0[m.i_ue] * m.N0 / (1.5 * ne)
    k = float(xsec.rate_table("F2", "attachment")(te))
    e_cap = xsec.capture_energy_table("F2", "attachment")(te)
    assert e_cap < 0.5 * 1.5 * te
    p_el = ne * n[m.idx["Kr"]] * float(m.kel(te)) * (1 - 300.0 * KB_EV / te)
    dy = m.rhs(0.0, y0, False)
    assert dy[m.i_ue] * m.N0 == approx(
        -e_cap * k * n[m.idx["F2"]] * ne - p_el, rel=1e-9)
    res = m.run()
    assert 1e-3 * res.ne[0] < res.ne[-1] < 0.5 * res.ne[0]
    assert res.te[-1] > 1.1 * res.te[0]
    lost = res.ledger["heat_capture"][-1] + res.ledger["heat_elastic"][-1]
    assert res.ue[0] - res.ue[-1] == approx(lost, rel=1e-6)


def test_basic_network_runs_the_two_temperature_closure():
    """v0_4's 'auto' gives the basic network the kinetic EEDF and the exact
    captured energy (v0_2's gave it the two-temperature tail), and keeps
    v0 on the notes' Maxwellian with (3/2)Te. v0_2's closure stays
    available as eedf='two_t': its five Kr rows read the tail table, and at
    Tahoma's n_e/N and F2 the tail row falls below the Maxwellian."""
    from v0_4 import xsec
    auto_b = Model(_Params(groups=("all",)))
    auto_v = Model(_Params(network="v0", groups=("all",)))
    assert (auto_b.eedf, auto_b.capture) == ("kinetic", "exact")
    assert (auto_v.eedf, auto_v.capture) == ("maxwellian", "three_halves")
    b = Model(Params(groups=("all",)))
    v = Model(Params(network="v0", groups=("all",), eedf="auto"))
    assert (b.eedf, b.capture) == ("two_t", "exact")
    assert (v.eedf, v.capture) == ("maxwellian", "three_halves")
    assert sorted(c.row.rid for c in b.compiled if c.kind == "xsec_tail") \
        == ["C1a", "C1b", "C1c", "C2", "G1"]
    assert not any(c.kind == "xsec_tail" for c in v.compiled)
    assert all(c.cap is None for c in v.compiled)
    cb = {c.row.rid: c for c in b.compiled}
    att = xsec.capture_energy_table("F2", "attachment")
    assert [cb["ATT"].cap(t) for t in (0.3, 1.3)] == [att(0.3), att(1.3)]
    assert cb["R11p"].cap(1.3) == approx(0.95 * 1.3, rel=1e-12)  # Te^-0.55
    n = np.zeros(b.ns)
    n[b.idx["Kr"]], n[b.idx["F2"]] = 2.0e19, 1.0e17
    te, ne = 1.6, 2.01e19 * 2.9e-5
    r = xsec.tail_ratio_table()(2.9e-5, 1.0e17 / 2.01e19)
    g1 = ("Kr", "excitation", "Kr(s)(9.915eV)")
    assert b._k(cb["G1"], te, 300.0, n, ne) == approx(
        xsec.two_t_table(*g1)(te, r * te), rel=1e-12)
    assert b._k(cb["G1"], te, 300.0, n, ne) < 0.3 * xsec.rate_table(*g1)(te)
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
    # K15: the mean of three 300 K measurements, on Perez-Rios & Greene's T^-3/4
    assert m._k(c["K15"], 1.0, 300.0, n, 1e15) == approx(2.4e-31, rel=1e-12)
    assert m._k(c["K15"], 1.0, 600.0, n, 1e15) == approx(2.4e-31 * 2 ** -0.75,
                                                         rel=1e-12)
    # G5: Shui, 'a change in temperature from 300 K to 600 K lowers kr(Kr)
    # by a factor of 10'
    k300 = m._k(c["G5"], 1.0, 300.0, n, 1e15)
    k600 = m._k(c["G5"], 1.0, 600.0, n, 1e15)
    assert k300 == approx(1.90e-30, rel=1e-12)
    assert k600 / k300 == approx(0.1, rel=1e-9)
    # K26 at 1 amagat (v0_4): the harmonic mean of REBUILD's three-body
    # limb k3 = 2.4 x 3.417e-6 / 1.7 and the Langevin limb of the Ellis
    # mobilities kL = 1.8095e-6 (0.84 + 2.29), times Debye S. v0_3's
    # Flannery form gave 2.4/2.7 x 3.417e-6 S = 3.04e-6 S; this is 2.60e-6 S
    n_am = n.copy()
    n_am[m.idx["Kr"]] = 2.6867811e19
    s_scr = m._screening(1e15, n_am, 300.0)
    assert 0.0 < s_scr < 1.0
    k3, kl = 2.4 * 3.417e-6 / 1.7, 1.8095e-6 * (0.84 + 2.29)
    assert m._k(c["K26"], 1.0, 300.0, n_am, 1e15) == approx(
        k3 * kl / (k3 + kl) * s_scr, rel=1e-9)
    assert k3 * kl / (k3 + kl) == approx(2.60e-6, rel=5e-3)
    # the Thomson (Tg/300)^-3 acts on the three-body limb alone, at fixed
    # N: at 600 K that limb is k3/8
    s600 = m._screening(1e15, n_am, 600.0)
    a3 = k3 * 0.125
    assert m._k(c["K26"], 1.0, 600.0, n_am, 1e15) == approx(
        a3 * kl / (a3 + kl) * s600, rel=1e-9)
    # the limbs: at low density the three-body one, k3 N; at high density
    # the Langevin one, kL / N (both x S)
    lo = n_am * 1e-3
    assert m._k(c["K26"], 1.0, 300.0, lo, 1e15) == approx(
        k3 * 1e-3 * m._screening(1e15, lo, 300.0), rel=2e-3)
    hi = n_am * 1e3
    assert m._k(c["K26"], 1.0, 300.0, hi, 1e15) == approx(
        kl / 1e3 * m._screening(1e15, hi, 300.0), rel=2e-3)
    # the yield and the C share multiply the three product rows
    tot = sum(m._k(c[r], 1.0, 300.0, n_am, 1e15)
              for r in ("K26", "K26C", "K26X"))
    assert tot == approx(k3 * kl / (k3 + kl) * s_scr, rel=1e-9)
    m2 = Model(Params(rates=(("ii_yield", 0.8), ("ii_f_c", 0.25))))
    c2 = {x.row.rid: x for x in m2.compiled}
    k26 = m2._k(c2["K26"], 1.0, 300.0, n_am, 1e15)
    assert k26 == approx(0.8 * 0.75 * tot, rel=1e-9)
    assert m2._k(c2["K26C"], 1.0, 300.0, n_am, 1e15) == approx(
        0.8 * 0.25 * tot, rel=1e-9)
    assert m2._k(c2["K26X"], 1.0, 300.0, n_am, 1e15) == approx(0.2 * tot,
                                                               rel=1e-9)
    # the te_fit rows: Kannari's forms, and K76r by detailed balance
    assert m._k(c["K76"], 1.3, 300.0, n, 1e15) == approx(
        9.0e-7 * 1.3 ** 0.5 * math.exp(-1.51 / 1.3), rel=1e-12)
    ratio = (m._k(c["K76r"], 1.3, 300.0, n, 1e15)
             / m._k(c["K76"], 1.3, 300.0, n, 1e15))
    assert ratio == approx(12 / 36 * math.exp(1.385 / 1.3), rel=1e-9)
    assert m._k(c["KSTEP_S"], 1.3, 300.0, n, 1e15) == 0.0     # band at 0
    mk = Model(Params(ladder="kannari"))
    ck = {x.row.rid: x for x in mk.compiled}
    assert mk._k(ck["KSTEP_S"], 1.3, 300.0, n, 1e15) == approx(
        7.0e-8 * 1.3 ** 0.68 * math.exp(-4.04 / 1.3), rel=1e-12)
    assert mk._k(ck["KSTEP_P"], 1.3, 300.0, n, 1e15) == approx(
        2.0e-7 * 1.3 ** 0.62 * math.exp(-2.53 / 1.3), rel=1e-12)
    assert mk._k(ck["KPS"], 1.3, 300.0, n, 1e15) == 1.0e-10
    assert mk._k(ck["K27"], 1.3, 300.0, n, 1e15) == approx(
        8.1e-10 * 0.9 * 0.62, rel=1e-12)
    # F + F + Kr: the Cohen & Bott evaluation, held at its 300 K value
    assert m._k(c["FF"], 1.0, 1300.0, n_am, 1e15) == approx(2.8e-34,
                                                            rel=1e-12)


TORR_CM3 = 133.322368 / (1.380649e-23 * 300.0) * 1e-6


def _krf_decay(kr_torr, f2_torr, net="basic"):
    """Seed KrF(B) in pure Kr/F2 at 300 K (no power, no electrons) and fit
    the slow decay of KrF(B) + KrF(C), as in the photolysis experiments."""
    tot = kr_torr + f2_torr
    p = Params(network=net, pressure_atm=tot / 760.0, x_f2=f2_torr / tot,
               power_W_cm3=0.0, pulse_s=0.0, t_end_s=400e-9,
               ne_mode="free", ne_cm3=1.0, n_save=2001,
               init_extra_cm3=(("KrF(B)", 1e12),))
    r = Model(p).run()
    n = r.n["KrF(B)"] + r.n.get("KrF(C)", 0.0)
    rel = n / n[0]
    sel = (r.t > 20e-9) & (rel < 1e-3) & (rel > 1e-8)
    return -np.polyfit(r.t[sel], np.log(n[sel]), 1)[0]


def _gadomski_eq5(kr_torr, f2_torr):
    """Gadomski et al. CPL 189, 153 (1992) eq. 5, the slow eigenvalue of the
    two-level model, with the basic network's own constants at 300 K (tau_B
    from Params, the rest read from network.basic_rows). B and C each lose
    to their own quench constants."""
    k, rad_b = _basic_consts()
    kr, f2 = kr_torr * TORR_CM3, f2_torr * TORR_CM3
    qb = k["QKR"] * kr + k["G5"] * kr * kr + k["QF2"] * f2
    qc = k["QKRC"] * kr + k["G5C"] * kr * kr + k["QF2C"] * f2
    k1, k3 = rad_b + qb, k["G4C"] + qc
    k2, km2 = k["BC"] * kr, k["CB"] * kr
    return 0.5 * ((k1 + k2 + km2 + k3)
                  - math.sqrt((k1 + k2 - km2 - k3) ** 2 + 4 * k2 * km2))


@pytest.mark.parametrize("kr,f2", [(50, 1), (100, 1), (500, 1), (100, 2),
                                   (300, 2), (50, 6)])
def test_basic_network_is_gadomskis_two_level_model(kr, f2):
    """Implementation check: the basic network's KrF(B,C) decay equals the
    analytic slow eigenvalue of the model the constants were fitted with."""
    assert _krf_decay(kr, f2) == approx(_gadomski_eq5(kr, f2),
                                               rel=1e-3)


def test_basic_network_against_the_measured_krf_decays():
    """Score, not fit: Xu, Gadomski & Setser 1993 p.2594 print k' = 5.0e7
    and 6.2e7 1/s at 50 and 100 Torr Kr with 1 Torr F2. The basic network
    (ML's joint fit, adopted 2026-09-28) gives x0.91 and x0.97. Before, with
    Hakuta's tau_C and k_BC, it gave x0.97 and x1.01, and Gadomski's own
    constants x0.89 and x0.95. v0's B-only KrF misses them x2.9 and x2.4.
    The Kansas State figures scatter by up to 30 % at one condition (4.1e7
    in Gadomski Fig. 2a and about 5.4e7 in Fig. 1 at 50 Torr), so the match
    is to that level. The constants were fitted to the Kansas State decays
    (with Eden's and Yu & Setser's), so this checks structure."""
    for kr, meas, lo, hi in ((50, 5.0e7, 0.86, 0.96), (100, 6.2e7, 0.92, 1.02)):
        assert lo < _krf_decay(kr, 1) / meas < hi
    assert _krf_decay(50, 1, "v0") / 5.0e7 > 2.5


def _basic_consts():
    """The basic network's B/C constants at 300 K (C -> B is k_BC/K_eq)."""
    from v0_4.network import basic_rows
    rows = {r.rid: r for r in basic_rows()}
    k = {rid: rows[rid].rate[1] for rid in
         ("BC", "G4C", "QKR", "QKRC", "G5", "G5C", "QF2", "QF2C", "QE",
          "QEC")}
    kind, k_fwd, k_ref = rows["CB"].rate[:3]
    assert kind == "balance"
    k["CB"] = k_fwd / k_ref
    return k, 1.0 / Params().tau_sp_s


def _sv_tau_b(kr):
    """Stern-Volmer lifetime of B-X emission, KrF formed in B, no F2 and no
    electrons, 300 K: d(I0/I)/d(k_e n_e) at n_e = 0, from the steady state of
    the B/C block (equal electron quench on B and C, as the network has)."""
    k, rad_b = _basic_consts()
    qb = k["QKR"] * kr + k["G5"] * kr * kr
    qc = k["QKRC"] * kr + k["G5C"] * kr * kr
    assert k["QE"] == k["QEC"]
    a = rad_b + qb + k["BC"] * kr                # B out
    d = k["G4C"] + qc + k["CB"] * kr             # C out
    b2c, c2b = k["BC"] * kr, k["CB"] * kr
    det = a * d - b2c * c2b
    b0 = d / det                                 # steady B per unit source
    m_b = (d * b0 + c2b * (b2c / det)) / det     # (M^-1 M^-1 R)_B
    return m_b / b0


def test_electron_quench_is_trainor_jacob_under_this_b_c_model():
    """k_e is tied to the B/C constants. Trainor & Jacob measured the product
    k_e tau / k_att = 2.9e-7 s at 190 Torr Kr (Duzy et al., AVCO 1980, p.21).
    In this network tau is the two-level Stern-Volmer lifetime of B-X
    emission, 5.8 ns (8.3 ns with the constants before 2026-09-28, 6.5 ns in
    their one-state reading). With LIGHT's own F2 attachment rate at the
    1.58 eV Kushner, Hanson & Schneider computed for the experiment, and
    their recombination factor 0.62, the product gives 6.6e-8, the QE row.
    If the B/C constants change, re-derive k_e (ML forward.
    with_rederived_k_e, docs/krf_data/electron_quench.py)."""
    from v0_4 import xsec
    tau = _sv_tau_b(190.0 * 3.2188e16)
    assert tau == approx(5.80e-9, rel=0.01)
    k_att = float(xsec.rate_table("F2", "attachment")(1.58))
    k_e = 2.9e-7 * k_att * 0.62 / tau
    k, _ = _basic_consts()
    assert k_e == approx(k["QE"], rel=0.10)


def test_basic_network_against_eden_b_only_decays():
    """Score (a known failure until 2026-09-28). Eden, Waynant, Searles &
    Burnham,
    APL 32, 733 (1978) photolysed KrF2 at 193 nm to make vibrationally cold
    KrF(B) alone in 2-200 Torr Kr and printed, from single-exponential fits
    with the pump deconvolved, Sigma = 1/tau + 8.6e-12 [Kr] + 9.7e-31 [Kr]^2
    (tau 6.8 ns, plus 0.5 Torr KrF2 at 1.4e-10). His signals decayed
    "exponentially over several e-foldings". With Hakuta's measured k_BC =
    1.2e-10 the network's B(t) after a B-only start fell to 21-51 % within
    2 ns (the transfer to C), and the fitted rate was x1.35, x0.95 and x0.82
    of Eden's at 50, 100 and 200 Torr. With ML's joint fit (slow mixing,
    k_BC 1.05e-11) it is x0.96, x0.95 and x0.98. Hakuta is now the known
    conflict (test_basic_network_against_hakuta_c_decays)."""
    n_torr = 3.22e16                             # Eden's figure density axis
    for p_torr in (50.0, 100.0, 200.0):
        kr = p_torr * n_torr
        eden = (1 / 6.8e-9 + 1.4e-10 * 0.5 * n_torr + 8.6e-12 * kr
                + 9.7e-31 * kr * kr)
        prm = Params(network="basic", pressure_atm=p_torr / 760.0,
                     x_f2=1e-9, power_W_cm3=0.0, pulse_s=0.0,
                     t_end_s=150e-9, ne_mode="free", ne_cm3=1.0,
                     n_save=3001, init_extra_cm3=(("KrF(B)", 1e12),))
        r = Model(prm).run()
        b = r.n["KrF(B)"] / r.n["KrF(B)"][0]
        sel = b > 0.05
        sigma = -np.polyfit(r.t[sel], np.log(b[sel]), 1)[0]
        assert sigma / eden == approx(1.0, abs=0.10)


@pytest.mark.xfail(strict=True, reason="known conflict: Hakuta 1990 is "
                   "inconsistent with Eden, Kansas State and Yu & Setser "
                   "(ML consistency check, chi2 12 per point)")
def test_basic_network_against_hakuta_c_decays():
    """Score, recorded as a known failure. Hakuta, Fujino, Nakayama &
    Takuma, JJAP 29, 1194 (1990) Figs. 2-3: the C-A decay rate after 248-nm
    photoassociation against F2 (9 Torr Kr) and against Kr (1 Torr F2),
    digitized in docs/krf_data/digitized/. The basic network's slow
    eigenvalue runs 17 % (Fig. 2) and 31 % (Fig. 3) below them, 19 % and
    34 % rms. Their rates are about 1.8x the Kansas State rates at the same
    conditions. Strict: if new data or constants make it pass, remove the
    xfail."""
    import csv
    from v0_4 import two_level as L
    c = L.Consts.from_light()
    base = Path(__file__).resolve().parents[2] / "docs" / "krf_data" / "digitized"
    for fn, fixed in (("hakuta1990_fig3_Kr.csv", "f2"),
                      ("hakuta1990_fig2_F2.csv", "kr")):
        with open(base / fn) as fh:
            pts = [(float(a), float(b)) for a, b in csv.reader(
                ln for ln in fh if ln[0].isdigit())]
        rel = []
        for p_torr, v in pts:
            kr, f2 = (p_torr, 1.0) if fixed == "f2" else (9.0, p_torr)
            slow = L.eigen(c, kr=kr * 3.22e16, f2=f2 * 3.22e16)[1]
            rel.append(slow / (v * 1e7) - 1.0)
        assert math.sqrt(np.mean(np.square(rel))) < 0.15, fn


def test_c_to_b_follows_detailed_balance():
    """K_eq = k_BC/k_CB (1.10 at 300 K since 2026-09-28) follows the gas
    temperature with the B-C gap the energy ledger books (species.py), so
    the kinetics and the heat routed by B <-> C exchange use one number."""
    from v0_4.network import KRF_K_BC, KRF_K_EQ
    from v0_4.species import e_form
    m = Model(Params())
    rows = {c.row.rid: c for c in m.compiled}
    theta = (e_form("KrF(B)") - e_form("KrF(C)")) / KB_EV
    assert theta == approx(120.0 * 1.23984193e-4 / KB_EV, rel=1e-12)
    zero = np.zeros(m.ns)
    for tg in (300.0, 460.0, 800.0):
        keq = (m._k(rows["BC"], 1.0, tg, zero, 0.0)
               / m._k(rows["CB"], 1.0, tg, zero, 0.0))
        assert keq == approx(
            KRF_K_EQ * math.exp(theta * (1 / tg - 1 / 300)), rel=1e-12)
    assert m._k(rows["BC"], 1.0, 300.0, zero, 0.0) == approx(KRF_K_BC)


# ---------------------------------------------------------------------------
# v0_1: recombination, fluorine recombination, the trimer switch
# ---------------------------------------------------------------------------
# v0_4's reference point under v0_2's closure (two_t), pinned once its
# new rows carried their literature values (2026-10-08); see the README
V0_4_REF = dict(emission_248=5.787, emission_vuv=4.53, tg_end_K=742.0)


def test_reference_point_regression():
    """v0_3's reference point (v0_2's network and numbers) is recovered
    when the 5p rung is neutralised: the 5s <-> 5p electron cycle off
    (k_5s5p = 0) and the 5p relaxed to 5s at once (kps_kr large), so
    that the higher excitations end in Kr* as v0_3's prompt cascade had
    them. The ion-ion limb change (v0_4) moves the 248-nm share by less
    than 0.1 %. With Kannari's Kr2* estimates on top it is v0_1's kr2
    group (9.17 %), and with the Penning one alone 4.65 % (CHANGELOG
    2026-09-29 (b)). The seed only sets the first nanoseconds."""
    from v0_4.report import warnings
    off = (("k_5s5p", 0.0), ("kps_kr", 1e-10))
    res = run(rates=off)
    d = row_dict(res)
    assert res.model.ne_mode == "free"
    assert d["emission_248"] == approx(4.667, rel=2e-3)
    assert d["emission_vuv"] == approx(20.0, rel=5e-3)
    assert d["tg_end_K"] == approx(649.4, rel=2e-3)
    assert abs(d["closure"]) < 1e-12
    assert any("Kr* + 2Kr's (Tg/300)^-0.64" in w for w in warnings(res))
    kannari = row_dict(run(rates=off + (("kr2_f", 3.0e-10),
                                        ("kr2_pen", 5.0e-10))))
    assert kannari["emission_248"] == approx(9.173, rel=4e-3)
    assert kannari["tg_end_K"] == approx(723.6, rel=2e-3)
    pen = row_dict(run(rates=off + (("kr2_pen", 5.0e-10),)))
    assert pen["emission_248"] == approx(4.652, rel=2e-3)
    seed = row_dict(run(ne_cm3=1e12, rates=off))
    assert seed["emission_248"] == approx(d["emission_248"], rel=0.02)
    # v0_4's own defaults, pinned
    if V0_4_REF["emission_248"] is not None:
        d4 = row_dict(run())
        for k, v in V0_4_REF.items():
            assert d4[k] == approx(v, rel=3e-3), k


def test_fluorine_recombination_is_too_slow_to_refill_f2():
    """F + F + Kr at the Cohen & Bott evaluation re-forms under 1 % of the
    F2 fill within the run (the F atoms' recombination time is
    milliseconds at 1 atm)."""
    res = run()
    assert 0.0 < res.events["FF"][-1] < 0.01 * res.n["F2"][0]


def test_trimer_switch_moves_only_the_trimer():
    """trimer_tg_pow spans the unmeasured high-temperature law for a
    sensitivity run. It changes G5 and G5C and nothing else."""
    base = Model(Params())
    alt = Model(Params(trimer_tg_pow=-0.5))
    n = np.zeros(base.ns)
    n[base.idx["Kr"]] = base.N0
    cb = {c.row.rid: c for c in base.compiled}
    ca = {c.row.rid: c for c in alt.compiled}
    for rid in cb:
        if cb[rid].kind == "elastic":        # an energy loss, not a rate
            continue
        kb = base._k(cb[rid], 1.5, 700.0, n, 1e14)
        ka = alt._k(ca[rid], 1.5, 700.0, n, 1e14)
        if rid in ("G5", "G5C"):
            assert ka / kb == approx((700 / 300) ** (-0.5 + 3.321928094887362),
                                     rel=1e-9)
        else:
            assert ka == approx(kb, rel=1e-12), rid


def test_held_n_e_with_attachment_is_flagged():
    """A held n_e with attachment in the network makes charge from nothing
    (CHANGELOG Known-bad). The report says so."""
    from v0_4.report import warnings
    res = run(**FAST, ne_mode="clamped")
    assert any("charge is not conserved" in w for w in warnings(res))
    assert not any("charge is not conserved" in w for w in warnings(run(**FAST)))



# ---------------------------------------------------------------------------
# v0_2: the Kr2* route and the bands
# ---------------------------------------------------------------------------
def test_param_rows_read_their_band():
    """A ("param", name[, frac]) row takes frac x Params.rates[name], or the
    band's default. A band no row reads, or a negative rate, is refused."""
    from v0_4.network import BANDS
    c = {x.row.rid: x for x in Model(Params()).compiled}
    for rid, name, frac in (("G2", "kr_pen", 1.0), ("Q2F", "kr2f_f2", 1.0),
                            ("Q2FF", "kr2f_f", 1.0),
                            ("AK8_Kr", "kr2_form", 1.0),
                            ("NK7", "kr2_f2", 0.5), ("K2F_BR", "kr2_f2", 0.5),
                            ("K2F_F", "kr2_f", 1.0),
                            ("PEN_KR2", "kr2_pen", 1.0)):
        assert c[rid].kind == "param" and c[rid].cap is None, rid
        assert c[rid].rate == approx(frac * BANDS[name].default,
                                     rel=1e-12), rid
    assert c["K2F_F"].rate == c["PEN_KR2"].rate == c["Q2FF"].rate == 0.0
    m = Model(Params(rates=(("kr2_f", 3e-10), ("kr2_f2", 5e-10))))
    c = {x.row.rid: x for x in m.compiled}
    assert c["K2F_F"].rate == approx(3e-10, rel=1e-12)
    assert c["NK7"].rate == c["K2F_BR"].rate == approx(2.5e-10, rel=1e-12)
    n = np.zeros(m.ns)
    n[m.idx["Kr"]] = m.N0
    assert m._k(c["K2F_F"], 1.0, 600.0, n, 1e15) == approx(3e-10, rel=1e-12)
    assert m._k(c["AK8_Kr"], 1.0, 600.0, n, 1e15) == approx(
        1.7e-32 * 2.0 ** -0.64, rel=1e-12)
    for bad in ((("k2f", 1e-10),), (("kr2_f", -1e-10),),
                (("trimer_tg_pow", -0.5),)):
        with pytest.raises(ValueError):
            Model(Params(rates=bad))
    with pytest.raises(ValueError):                    # v0 reads no band
        Model(Params(network="v0", rates=(("kr2_f", 3e-10),)))


def test_bands_run_each_end_once():
    """--bands runs the network as given, then each band's two ends one at a
    time. An end at the band's current value is that run again."""
    from v0_4.network import BANDS
    from v0_4.run import band_table, do_bands, with_band
    base = Params(pulse_s=0.1e-6, t_end_s=0.2e-6, n_save=101)
    assert with_band(base, "trimer_tg_pow", -0.5).trimer_tg_pow == -0.5
    assert with_band(base.replace(rates=(("kr2_f", 1e-10),)), "kr2_f",
                     2e-10).rates == (("kr2_f", 2e-10),)
    rows = do_bands(base)
    assert len(rows) == 1 + 2 * len(BANDS) and rows[0]["band"] == "as given"
    ref = rows[0]["emission_248"]
    by = {(r["band"], r["end"]): r for r in rows[1:]}
    for b in BANDS.values():
        for end, v in (("low", b.lo), ("high", b.hi)):
            assert by[(b.name, end)]["value"] == v
            if v == b.default:
                assert by[(b.name, end)]["emission_248"] == ref
    assert by[("kr2_f", "high")]["emission_248"] > ref
    assert "Kr2* + F -> KrF(B) + Kr" in band_table(rows)
    with pytest.raises(SystemExit):
        do_bands(base.replace(network="v0"))


def test_kr2_warning_needs_the_route_and_hot_gas():
    """The Kr* + 2Kr temperature law is measured to 451 K. The report warns
    above it, and not when the gas stays cool or the route is off."""
    from v0_4.report import warnings
    cool = run(**FAST, power_W_cm3=3e4)
    assert cool.tg.max() < 451.0
    assert not any("Kr* + 2Kr" in w for w in warnings(cool))
    m = Model(Params(**FAST, power_W_cm3=1.5e6, rates=(("kr2_form", 0.0),)))
    hot = m.run()
    assert hot.tg.max() > 451.0
    assert not any("Kr* + 2Kr" in w for w in warnings(hot))
