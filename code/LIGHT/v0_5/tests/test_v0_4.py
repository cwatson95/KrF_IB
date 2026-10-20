"""v0_4 (2026-10-08): the 5p rung, KrF* + F, the ion-ion bands, the ions'
drift to the walls, the 248-nm optics and the EM wave with Turnbull's
inverse bremsstrahlung. The ledger must close in every new mode, each new
term must reproduce its own closed form, and the EM wave's absorbed power
under the kinetic closure must agree with the Maxwellian kernel where the
two should agree.

Fast tests run on v0_2's closure (eedf='two_t'); the kinetic EM test is
short (a 20-ns pulse)."""

import math

import numpy as np
import pytest

from v0_5 import ib, optics
from v0_5.model import KB_EV, QE, Model, WALL_E_AMBIPOLAR_TE, \
    WALL_E_ELECTRON_TE
from v0_5.model import Params as _Params
from v0_5.model import run as _run
from v0_5 import network
from v0_5.network import ION_K0, LADDER_KANNARI
from v0_5.report import accounting, reservoir, row_dict, summary
from v0_5.species import e_form

AMAGAT = 2.6867811e19


def Params(**kw):
    kw.setdefault("eedf", "two_t")
    return _Params(**kw)


def run(params=None, **kw):
    if params is None:
        kw.setdefault("eedf", "two_t")
    return _run(params, **kw)


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


FAST = dict(pulse_s=0.3e-6, t_end_s=0.6e-6)


def _heat_booked(res):
    m = res.model
    return (res.ledger["heat_elastic"] + res.ledger["heat_capture"]
            + res.ledger["E_ion_heat"]
            + sum(res.events[c.row.rid] * c.route["heat"] for c in m.compiled))


def _check_closures(res, tol=1e-9):
    a = accounting(res)
    assert abs(a["closure"]) < tol
    e_dep = a["e_dep"]
    rad = sum(res.events[c.row.rid][-1] * c.route["rad"]
              for c in res.model.compiled)
    raw = (res.ledger["E_dep"][-1] + res.ledger["E_em"][-1]
           - (res.stored[-1] - res.stored[0]) - rad - _heat_booked(res)[-1]
           - res.ledger["E_wall"][-1]) / e_dep
    assert abs(raw) < tol
    # the gas's thermal energy: the heat booked, less what the particles
    # lost to the walls carried out at Tg
    gap = np.abs((res.u_thermal - res.u_thermal[0])
                 - (_heat_booked(res) - res.ledger["E_wall_th"]))
    assert gap.max() / e_dep < tol
    for name in ("Kr*", "Kr**", "KrF*"):
        rv = reservoir(res, name)
        assert abs(rv["balance"]) <= 1e-6 * max(rv["formed"], 1.0)


# ---------------------------------------------------------------------------
# closures in every new mode
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kw", [
    dict(),
    dict(ladder="kannari"),
    dict(rates=(("krf_f", 7.8e-10),)),
    dict(rates=(("ii_yield", 0.5), ("ii_f_c", 0.5))),
    dict(gap_cm=1.0),
    dict(gap_cm=0.1, pressure_atm=0.1),
    dict(em_E_V_cm=1500.0, em_f_Hz=95e9),
    dict(em_E_V_cm=1500.0, em_f_Hz=95e9, power_W_cm3=0.0),
    dict(em_E_V_cm=3000.0, em_f_Hz=1e12, rates=(("ib_scale", 1.1),)),
    dict(ladder="kannari", gap_cm=1.0, em_E_V_cm=1000.0, em_f_Hz=3e11,
         rates=(("krf_f", 7.8e-10),)),
])
def test_closures_in_the_new_modes(kw):
    _check_closures(run(**FAST, **kw))


def test_neutralised_ladder_is_v0_3():
    """With the 5s <-> 5p cycle off and the 5p relaxed at once, the higher
    excitations end in Kr* as in v0_3, and the reference point is v0_3's
    (4.667 % at 248 nm, VUV 20.0 %) to the ion-ion limb change."""
    d = row_dict(run(rates=(("k_5s5p", 0.0), ("kps_kr", 1e-10))))
    assert d["emission_248"] == approx(4.667, rel=2e-3)
    assert d["emission_vuv"] == approx(20.0, rel=5e-3)


# ---------------------------------------------------------------------------
# the 5p rung
# ---------------------------------------------------------------------------
def test_ladder_detailed_balance_and_energy():
    """A seed of Kr* and Kr** in pure Kr with hot electrons and no drive,
    with the collisional and radiative losses off: the electrons alone
    move the two levels to their Boltzmann ratio g_p/g_s exp(-dE/Te), and
    the energy they take and return nets to the population change times
    the gap."""
    te = 2.0
    # the gas at the electrons' temperature, so the elastic exchange (not
    # a row, always on) moves no energy and Te stays put; with the gas at
    # 300 K the electrons cooled to 1.1 eV in 20 ns and the populations
    # lagged Te by 1 %
    m = Model(Params(network="basic", pressure_atm=1.0, x_f2=1e-12,
                     power_W_cm3=0.0, pulse_s=0.0, t_end_s=1e-9,
                     ne_mode="clamped", ne_cm3=1e16, tg0_K=te / KB_EV,
                     rates=(("kps_kr", 0.0), ("k_5s5p", 1.0)),
                     init_extra_cm3=(("Kr*", 1e12), ("Kr**", 1e10))))
    # take the radiative row and every other Kr* loss out: the cycle alone
    # (the event slots of the dropped rows stay, at zero)
    m.compiled = [c for c in m.compiled if c.row.rid in ("K76", "K76r")]
    y0 = m.initial_state()
    y0[m.i_ue] = 1.5 * 1e16 * te / m.N0
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, y: m.rhs(t, y, False), (0.0, 20e-9), y0,
                    method="BDF", rtol=1e-10, atol=1e-30)
    y = sol.y[:, -1]
    n_s, n_p = y[m.idx["Kr*"]], y[m.idx["Kr**"]]
    gap = e_form("Kr**") - e_form("Kr*")
    te_end = y[m.i_ue] * m.N0 / (1.5 * 1e16)
    assert n_p / n_s == approx(36 / 12 * math.exp(-gap / te_end), rel=1e-4)
    # the electrons paid the gap for every net Kr* -> Kr**, less the
    # elastic exchange with the gas (always on, and not zero once the
    # cycle has moved Te a hair off Tg), which the ledger books
    d_p = (y[m.idx["Kr**"]] - y0[m.idx["Kr**"]]) * m.N0
    d_ue = (y[m.i_ue] - y0[m.i_ue]) * m.N0
    e_el = (y[m.i_led + 1] - y0[m.i_led + 1]) * m.N0
    assert d_ue == approx(-gap * d_p - e_el, rel=1e-6)
    assert abs(e_el) < 0.05 * abs(gap * d_p)


def test_5p_rung_routes_the_higher_excitations():
    """At the reference point the C1 rows make Kr**, which the electrons
    return or the radiation drains (the default bands), and Kr** balances;
    with the ladder preset the 5p relaxes to 5s and the electrons ionize
    it, and the near-IR band shrinks."""
    res = run(**FAST)
    rv = reservoir(res, "Kr**")
    formed = rv["formed"]
    assert formed > 0.0
    c1 = sum(res.events[r][-1] for r in ("C1a", "C1b", "C1c"))
    assert c1 + res.events["R11p"][-1] + res.events["K76"][-1] == approx(
        formed, rel=1e-9)
    assert set(rv["losses"]) >= {"K76r", "PRAD", "KPS", "K27"}
    assert rv["losses"]["KPS"] > rv["losses"]["PRAD"]      # 1.7e-10 x N >> 1/29 ns
    assert rv["losses"]["KSTEP_P"] == 0.0                   # no measurement
    a = accounting(res)
    ir = a["bands"]["5p-5s (near IR)"] / a["e_dep"]
    assert ir < 0.02
    # with the relaxation and the 5p harpoon off the 5p can only radiate
    # or be returned by the electrons, and the near-IR band takes over a
    # tenth of the input (the first build's 17 %, CHANGELOG 2026-10-08)
    a0 = accounting(run(**FAST, rates=(("kps_kr", 0.0), ("k5p_f2", 0.0))))
    ir0 = a0["bands"]["5p-5s (near IR)"] / a0["e_dep"]
    assert ir0 > 0.05 and ir0 > 10 * ir          # 7.6 % on the short pulse
    k = row_dict(run(**FAST, ladder="kannari"))
    assert k["heat_reactions"] > 0.9 * row_dict(res)["heat_reactions"]


def test_ladder_preset_sets_every_ladder_band():
    m = Model(Params(ladder="kannari"))
    for name, v in LADDER_KANNARI.items():
        assert m.rates[name] == v
    m2 = Model(Params(ladder="kannari", rates=(("kps_kr", 3e-11),)))
    assert m2.rates["kps_kr"] == 3e-11          # a named band wins
    with pytest.raises(ValueError):
        Model(Params(ladder="all_on"))


def test_krf_f_quench_outruns_radiation_after_burnout():
    """At the band's top the F atoms left by the burnt F2 quench KrF(B)
    faster than it radiates late in the pulse, and the 248-nm share falls
    by a quarter."""
    base = row_dict(run())
    top = run(rates=(("krf_f", 7.8e-10),))
    d = row_dict(top)
    assert d["emission_248"] < 0.8 * base["emission_248"]
    k = {c.row.rid: c for c in top.model.compiled}
    f_end = top.n["F"][top.t <= top.p.pulse_s][-1]
    assert 7.8e-10 * f_end > 1.43e8
    assert reservoir(top, "KrF*")["losses"]["QF"] > 0.0


# ---------------------------------------------------------------------------
# ion-ion recombination
# ---------------------------------------------------------------------------
def test_ion_ion_products_follow_the_bands():
    res = run(**FAST, rates=(("ii_yield", 0.8), ("ii_f_c", 0.25)))
    ev = {k: v[-1] for k, v in res.events.items()}
    tot = ev["K26"] + ev["K26C"] + ev["K26X"]
    assert ev["K26"] == approx(0.8 * 0.75 * tot, rel=1e-9)
    assert ev["K26C"] == approx(0.8 * 0.25 * tot, rel=1e-9)
    assert ev["K26X"] == approx(0.2 * tot, rel=1e-9)
    base = run(**FAST)
    eb = {k: v[-1] for k, v in base.events.items()}
    assert eb["K26C"] == 0.0 and eb["K26X"] == 0.0 and eb["K108X"] == 0.0
    with pytest.raises(ValueError):
        Model(Params(rates=(("ii_yield", 1.5),)))


# ---------------------------------------------------------------------------
# the ions' drift to the walls
# ---------------------------------------------------------------------------
def test_wall_loss_is_the_ambipolar_rate_of_the_measured_mobilities():
    """A seeded, undriven plasma in a gap loses its cations at D_a /
    Lambda^2 with D_a = mu_i k (Te + Tg)/e, mu_i = K0 N0/N, Lambda =
    gap/pi, and F- at its free diffusion."""
    m = Model(Params(gap_cm=0.5, pressure_atm=0.5, tg0_K=300.0))
    n = np.zeros(m.ns)
    n[m.idx["Kr"]] = m.N0
    n[m.idx["Kr+"]] = 1e12
    n[m.idx["Kr2+"]] = 2e12
    n[m.idx["F-"]] = 5e11
    te, tg = 1.2, 400.0
    ue = 1.5 * te * 3e12                 # U_e over the cations: Te_w = te
    assert m.wall_te(ue, n) == approx(te, rel=1e-12)
    loss = {m.species[i]: (q, f) for i, q, f in m.wall_loss(ue, tg, n)}
    lam2 = (0.5 / math.pi) ** 2
    n_am = n.sum() / AMAGAT
    for s, nd in (("Kr+", 1e12), ("Kr2+", 2e12)):
        mu = ION_K0[s] / n_am
        assert loss[s][1] == approx(mu * (te + tg * KB_EV) / lam2 * nd,
                                    rel=1e-9), s
    mu_f = ION_K0["F-"] / n_am
    assert loss["F-"][1] == approx(mu_f * tg * KB_EV / lam2 * 5e11, rel=1e-9)
    assert Model(Params()).wall_loss(ue, tg, n) == ()
    # the drive is smooth and bounded where the electrons are gone
    assert m.wall_te(1e30, n) == 30.0 and m.wall_te(0.0, n) == 0.0


def test_wall_loss_books_its_energy_and_is_small_in_a_cm_gap():
    """The walls take the lost charges' formation energy and 2 Te per
    pair, the gas Te per pair, and the ledger closes. In a 1-cm gap at the
    reference point the walls take well under a percent of the input,
    as 1Dv0_3's resolved sheaths found (0.008 %)."""
    res = run(**FAST, gap_cm=1.0)
    _check_closures(res)
    a = accounting(res)
    assert 0.0 < a["walls"] / a["e_dep"] < 1e-3
    assert res.ledger["E_ion_heat"][-1] > 0.0
    # a thin gap at low pressure loses a visible share
    thin = run(**FAST, gap_cm=0.05, pressure_atm=0.1)
    _check_closures(thin)
    at = accounting(thin)
    assert at["walls"] / at["e_dep"] > 10 * a["walls"] / a["e_dep"]
    assert "walls" in summary(thin)
    # the electron's and the ion's shares per lost pair
    assert (WALL_E_ELECTRON_TE, WALL_E_AMBIPOLAR_TE) == (2.0, 1.0)


# ---------------------------------------------------------------------------
# optics
# ---------------------------------------------------------------------------
def test_optics_is_the_sum_over_the_absorbers():
    res = run(**FAST)
    g = optics.gain_cm(res)
    assert np.allclose(g, 2.6e-16 * res.n["KrF(B)"])
    al = optics.absorption_cm(res)
    tot = sum(sig * res.n[s] for s, (sig, _) in optics.ABSORBERS.items()
              if s in res.n)
    assert np.allclose(al["total"], tot)
    w = optics.window(res)
    assert w["net"] == approx(w["g0"] - w["alpha"], rel=1e-9)
    assert sum(w["by_species"].values()) == approx(w["alpha"], rel=1e-9)
    d = row_dict(res)
    assert d["g0_pct_cm"] == approx(w["g0"]) and \
        d["alpha_pct_cm"] == approx(w["alpha"])
    assert "OPTICS at 248 nm" in summary(res)
    for s in ("F-", "Kr2+", "Kr**", "F2"):
        assert s in optics.ABSORBERS


# ---------------------------------------------------------------------------
# row_scale: a source's quoted error through the model
# ---------------------------------------------------------------------------
def test_row_scale_multiplies_one_rows_rate_and_nothing_else():
    """Params.row_scale scales the named row's compiled rate by its factor
    (every rate kind), leaves the others alone, and refuses an unknown
    id or a negative factor. It is the sensitivity hook of the accuracy
    overview's propagated uncertainty."""
    m0 = Model(Params(eedf="two_t"))
    m1 = Model(Params(eedf="two_t", row_scale=(("R11p", 1.1), ("C3", 0.9))))
    f0 = {c.row.rid: c.factor for c in m0.compiled}
    f1 = {c.row.rid: c.factor for c in m1.compiled}
    assert f1["R11p"] == approx(1.1 * f0["R11p"])
    assert f1["C3"] == approx(0.9 * f0["C3"])
    assert all(f1[k] == f0[k] for k in f0 if k not in ("R11p", "C3"))
    n = np.zeros(m0.ns)
    n[m0.idx["Kr"]] = 0.995 * m0.N0
    n[m0.idx["Kr2+"]] = 1e14
    assert m1._k(m1.compiled[[c.row.rid for c in m1.compiled].index("R11p")],
                 1.5, 300.0, n, 1e14) == approx(
        1.1 * m0._k(m0.compiled[[c.row.rid for c in m0.compiled].index("R11p")],
                    1.5, 300.0, n, 1e14))
    with pytest.raises(ValueError):
        Model(Params(row_scale=(("NO_SUCH_ROW", 1.1),)))
    with pytest.raises(ValueError):
        Model(Params(row_scale=(("R11p", -0.5),)))


def _hot_initial(m, te=2.0):
    """The model's initial state with its electrons at te [eV]."""
    y = m.initial_state()
    ne = m.electron_density(np.maximum(y[:m.ns], 0.0) * m.N0)
    y[m.i_ue] = 1.5 * ne * te / m.N0
    return y


def test_row_scale_reaches_the_kinetic_table():
    """Under the kinetic EEDF the rows whose rates the kinetic table carries
    (the cross-section rows, the superelastic row and R11p's recombination
    template) carry their row_scale in the processes the Boltzmann equation
    solves, as K76's band always did, so the rate a run uses is the scaled
    one. Found in LIGHT_Kr v0_1 on 2026-10-13 and fixed here on 2026-10-14:
    the table's rate was read without the row's factor, so a row_scale on
    R11p, G1 or C2 moved nothing on the default closure. The 5s <-> 5p pair
    is one cross section and its detailed-balance reverse, so K76 and K76r
    take one factor; the elastic row cannot be scaled (its exchange is the
    momentum transfer of the whole EEDF)."""
    from v0_5 import xsec
    g1 = "G1"
    base = Params(eedf="kinetic")
    m0 = Model(base)
    m1 = Model(base.replace(row_scale=(("R11p", 1.1),)))
    m2 = Model(base.replace(row_scale=((g1, 0.8),)))
    m3 = Model(base.replace(row_scale=(("K76", 1.1), ("K76r", 1.1))))
    # the processes the equation solves
    s0 = {pr.rows[0]: pr.p.sigma for pr in m0.kin.procs}
    s2 = {pr.rows[0]: pr.p.sigma for pr in m2.kin.procs}
    assert np.array_equal(s2[g1], 0.8 * s0[g1])
    assert all(np.array_equal(s2[k], s0[k]) for k in s0 if k != g1)
    rec0 = {j: cons for j, cons, _, _ in m0._kin_rec}
    rec1 = {j: cons for j, cons, _, _ in m1._kin_rec}
    assert rec1 == {j: approx(1.1 * v, rel=1e-15, abs=0.0)
                    for j, v in rec0.items()}
    # the rates a run reads at one state: no molecular ion and no excited
    # atom yet, so R11p's sink and the 5s <-> 5p collisions leave the EEDF
    # alone and their rates scale exactly
    rates = []
    for m in (m0, m1, m2, m3):
        y = _hot_initial(m)
        m._kinetic_update(0.0, y, False)
        n = np.maximum(y[:m.ns], 0.0) * m.N0
        ne = m.electron_density(n)
        te = y[m.i_ue] * m.N0 / (1.5 * ne)
        rates.append({c.row.rid: m._rate(c, te, y[m.i_tg], n, ne)
                      for c in m.compiled if c.kind != "elastic"})
    r0, r1, r2, r3 = rates
    assert r1["R11p"] == approx(1.1 * r0["R11p"], rel=1e-9, abs=0.0)
    assert all(r1[k] == approx(r0[k], rel=1e-9, abs=0.0)
               for k in r0 if k != "R11p")
    for rid in ("K76", "K76r"):
        assert r3[rid] == approx(1.1 * r0[rid], rel=1e-9, abs=0.0)
    # the 5s group's row at 0.8: the rate read is 0.8 of its measured cross
    # section's over the run's own EEDF, whose tail rises as the row takes
    # less
    sol = m2._kin_sol
    row = next(c.row for c in m2.compiled if c.row.rid == g1)
    p = xsec.find(*row.rate[1:4])
    assert r2[g1] == approx(0.8 * sol.grid.rate(p, sol.F) * 1e6, rel=1e-6,
                            abs=0.0)  # the table's centre sits at the
    #                                   solver's mean energy
    assert 0.8 < r2[g1] / r0[g1] < 1.0
    for bad in ((("G6", 1.1),), (("K76", 1.1),), (("K76r", 1.1),),
                (("K76", 1.1), ("K76r", 1.2))):
        with pytest.raises(ValueError):
            Model(Params(row_scale=bad))


# ---------------------------------------------------------------------------
# the EM wave
# ---------------------------------------------------------------------------
def test_em_wave_is_a_second_input_with_its_band():
    """Its absorbed energy is booked apart (E_em), enters the fate's
    denominator, and the ib_scale band moves the electron-ion share
    alone: at 1 THz in 1 atm Kr that share is a few percent at most."""
    kw = dict(power_W_cm3=0.0, em_E_V_cm=1500.0, em_f_Hz=95e9, **FAST)
    mid = run(**kw)
    a = accounting(mid)
    assert a["e_em"] > 0.0 and a["e_given"] == 0.0
    assert a["e_dep"] == approx(a["e_em"])
    assert a["e_em"] * QE > 1e-4                 # J/cm3: a real input
    lo = accounting(run(**kw, rates=(("ib_scale", 0.9),)))["e_em"]
    hi = accounting(run(**kw, rates=(("ib_scale", 1.1),)))["e_em"]
    assert lo < a["e_em"] < hi
    assert (hi - lo) / a["e_em"] < 0.2 * 0.6     # the e-i share is < 60 %
    assert "EM wave" in summary(mid) and "absorbed" in summary(mid)
    with pytest.raises(ValueError):
        Model(Params(em_E_V_cm=100.0, em_f_Hz=0.0))


def test_em_kernel_is_v0s_lorentz_gas_without_ions():
    """With no ions the Maxwellian kernel is v0's Lorentz-gas absorption
    (ib.py test 6 holds it to 1e-6); here the model's own call."""
    m = Model(Params(em_E_V_cm=1500.0, em_f_Hz=95e9))
    n = np.zeros(m.ns)
    n[m.idx["Kr"]] = 0.995 * m.N0
    n[m.idx["F2"]] = 0.005 * m.N0
    per_e = m._em_maxwell(1.3, n, 1e15)
    w = 2 * math.pi * 95e9
    ref = ib.em_power_per_electron(1.3, w, 1.5e5, m.N0 * 1e6,
                                   {"Kr": 0.995, "F2": 0.005}, 1e21, 0.0,
                                   ei_scale=0.0)
    assert per_e == approx(ref, rel=1e-9)       # no cations in n
    assert m.em_power(1.3, n, 1e15) == approx(per_e * 1e15, rel=1e-9)


def test_kinetic_em_power_matches_the_maxwellian_kernel_where_it_should():
    """Under the kinetic closure the wave's power is the solver's own HF
    field term (Hagelaar & Pitchford). For a Maxwellian EEDF handed to
    the solver it must equal the Maxwellian kernel without the Langdon
    factor, and the DC solver is unchanged at omega = 0."""
    from v0_5 import boltzmann as KB
    m = Model(Params(eedf="kinetic", em_E_V_cm=1500.0, em_f_Hz=95e9))
    n = np.zeros(m.ns)
    n[m.idx["Kr"]] = 0.995 * m.N0
    n[m.idx["F2"]] = 0.005 * m.N0
    n[m.idx["Kr+"]] = 1e15
    y = m.initial_state()
    y[:m.ns] = n / m.N0
    te = 1.5
    y[m.i_ue] = 1.5 * 1e15 * te / m.N0
    s, _, ne = m._kinetic_state(y, 0.0, True)
    assert s.omega == approx(2 * math.pi * 95e9) and s.en_em_td > 0.0
    assert s.n_ion == approx(1e21)
    g = m.kin._grid(12.0)
    m.kin.grid = g
    F = g.maxwellian(te)
    sol = KB.Solution(F, 0.0, 0.0, g.mean(F), False, 0, g)
    p_kin = m.kin.em_power(s, sol)
    p_kin_atoms = m.kin.em_power(s, sol, ions=False)
    w = 2 * math.pi * 95e9
    parts = ib.em_power_per_electron(te, w, 1.5e5, s.n_col,
                                     {"Kr": 0.995, "F2": 0.005}, 1e21, 1e21,
                                     parts=True, langdon=False)
    atoms = ib.em_power_per_electron(te, w, 1.5e5, s.n_col,
                                     {"Kr": 0.995, "F2": 0.005}, 1e21, 0.0,
                                     parts=True, langdon=False)
    # the atoms' part: the same kernel on the same cross sections (the
    # solver's grid against a fine quadrature, a few percent)
    assert p_kin_atoms == approx(parts["en"], rel=0.05)
    assert p_kin > p_kin_atoms
    # the ions' part is the MARGINAL one, the power with the ions in the
    # momentum transfer less the power without them, to 1 %. The kernel's
    # printed "ei" is the nu_ei share of the numerator, which is not the
    # same thing where nu_en is not << omega: here nu_en = 4.4 omega and
    # the share is 1.5x the marginal part
    marginal = parts["total"] - atoms["total"]
    assert (p_kin - p_kin_atoms) == approx(marginal, rel=0.01)
    assert parts["ei"] / marginal > 1.3
    assert parts["nu_en_mean"] / w > 3.0
    # omega = 0: the DC assembly is the old one (no HF part)
    s0 = KB.State(x=s.x, tg=s.tg, n_col=s.n_col, ne=s.ne, n_ion=0.0)
    parts0 = m.kin._parts(g, s0, 0.0, F)
    assert parts0[4][0] is None
    assert m.kin.em_power(s0, sol) == 0.0


def test_kinetic_run_with_the_wave_closes():
    """A short kinetic run with the wave on: the ledger closes, the EM
    energy is booked, and the solver's EM power is read from the table."""
    res = _run(_Params(eedf="kinetic", em_E_V_cm=1500.0, em_f_Hz=95e9,
                       pulse_s=20e-9, t_end_s=25e-9, n_save=101))
    _check_closures(res)
    a = accounting(res)
    assert a["e_em"] > 0.0
    # in the pulse the table carried the solver's own EM power
    on = [e for e in res.kinetic if e["on"] and not e.get("skipped")]
    assert on and all(e["pem"] > 0.0 for e in on[1:])
