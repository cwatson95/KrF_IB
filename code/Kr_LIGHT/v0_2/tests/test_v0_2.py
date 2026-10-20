"""LIGHT_Kr v0_2 (2026-10-14): v0_1 with the electron's collisions on
ground-state krypton booked level by level on Biagi's Magboltz 8.97
cross sections, and the Hornbeck-Molnar branch of the levels at or above
E(Kr2+) on the band hm_kr (0 by default). The gate:

* every other row, band and table is v0_1's, the species at 8.97's
  thresholds;
* at hm_kr = 0 the Hornbeck-Molnar rows carry nothing: the right-hand side
  is the one without them and the kinetic equation leaves them out;
* at hm_kr = 1 every such excitation ionizes associatively: Kr2+ is made
  from the levels, their cascade carries nothing, and the ledger closes on
  both closures;
* the kinetic equation's associative ionization: its electrons count in
  the growth term and the solution's energy balance holds;
* the swarm (swarm.py): LoKI-MC's Monte Carlo on the same file gives
  alpha/N within 4 % of the two-term solution at 20 to 50 Td, and the band
  raises alpha/N 17-fold at 14 Td, less as E/N rises.

Comparisons with v0_1 skip when v0_1 is not beside v0_2."""

import importlib
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

from v0_2 import boltzmann as KB
from v0_2 import network, swarm, xsec
from v0_2.model import Model, Params
from v0_2.report import accounting, row_dict
from v0_2.species import SPECIES, e_form

HERE = Path(__file__).resolve()
V0_1 = HERE.parents[2] / "v0_1"
needs_v0_1 = pytest.mark.skipif(not (V0_1 / "network.py").is_file(),
                                reason="v0_1 is not beside v0_2")
FAST = dict(eedf="two_t", pulse_s=0.3e-6, t_end_s=0.6e-6)
LEVEL_GROUPS = ("G1", "C1a", "C1b", "C1c", "HM")


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


def _is_level(rid):
    return rid.startswith(tuple(g + "_" for g in LEVEL_GROUPS))


# ------------------------------------------------------------- vs v0_1
@needs_v0_1
def test_v0_1s_rows_bands_and_tables_carry_over():
    """Every row that is not a ground-state level is v0_1's, field for
    field but for the text, in v0_1's order; v0_1's bands are unchanged and
    hm_kr is the one more; jung2006 and the two-temperature tail table are
    v0_1's byte for byte; the species are v0_1's at Biagi 8.97's
    thresholds (NIST's levels)."""
    old_net = importlib.import_module("v0_1.network")
    old_sp = importlib.import_module("v0_1.species")
    old = [r for r in old_net.build(("all",))
           if r.rid not in ("G1", "C1a", "C1b", "C1c")]
    new = [r for r in network.build(("all",)) if not _is_level(r.rid)]
    assert [r.rid for r in new] == [r.rid for r in old]
    for a, b in zip(old, new):
        for f in ("reactants", "products", "group", "tg_pow", "band",
                  "energy_to", "rebuild", "frac"):
            assert getattr(a, f) == getattr(b, f), (a.rid, f)
        if a.rid == "SUPER":            # v0_2: the lumped 5s at 8.97's
            assert b.rate == ("superelastic", "Kr", xsec.LUMPED_5S, 1 / 12)
        elif a.rid in ("G6", "C2"):     # the same blocks of the new file
            assert b.rate[:3] == a.rate[:3]
        else:
            assert b.rate == a.rate, a.rid
    by = {r.rid: r for r in new}
    assert by["PRAD"].hv_eV == approx(e_form("Kr**") - e_form("Kr*"),
                                      rel=1e-12)
    assert by["C2"].e_cost_eV == 13.9996
    assert {k: astuple(v) for k, v in network.BANDS.items()
            if k != "hm_kr"} == {k: astuple(v) for k, v in
                                 old_net.BANDS.items()}
    assert list(network.BANDS)[-1] == "hm_kr"
    # jung2006's data are v0_1's (its docstring gives v0_2's energies), and
    # the two-temperature tail table is v0_1's file byte for byte
    j1 = importlib.import_module("v0_1.jung2006")
    j2 = importlib.import_module("v0_2.jung2006")
    for name in ("LEVELS", "PROXY", "SHAPE_REF", "TABLES", "TABLE_II",
                 "TABLE_III", "UNCERTAINTY", "WEIGHTS_5S", "E_MAX_TABLE",
                 "CM2"):
        assert getattr(j2, name) == getattr(j1, name), name
    for e in (1.5, 2.0, 4.0, 6.0, 20.0, 100.0):
        assert float(j2.lumped_sigma(e)) == float(j1.lumped_sigma(e))
    f = "data/kr_tail_ratio.txt"
    assert (HERE.parents[1] / f).read_bytes() == (V0_1 / f).read_bytes()
    assert set(SPECIES) == set(old_sp.SPECIES)
    shift = {"Kr*": 9.9152, "Kr**": 11.3035, "Kr+": 13.9996,
             "Kr2+": 13.9996 - 1.15}
    for s, sp in SPECIES.items():
        o = old_sp.SPECIES[s]
        assert (sp.atoms, sp.charge, sp.cv) == (o.atoms, o.charge, o.cv)
        assert sp.e_form_eV == approx(shift.get(s, o.e_form_eV), rel=1e-12)
        assert abs(sp.e_form_eV - o.e_form_eV) < 4e-3, s


# --------------------------------------------- the Hornbeck-Molnar rows
def _state(m, te=1.8):
    y = m.initial_state()
    y[m.idx["Kr*"]] = 1e-5
    y[m.idx["Kr2+"]] = 2e-6
    y[m.idx["Kr+"]] = 1e-6
    ne = m.electron_density(np.maximum(y[:m.ns], 0.0) * m.N0)
    y[m.i_ue] = 1.5 * ne * te / m.N0
    return y


def test_hm_rows_at_zero_change_nothing():
    """At hm_kr = 0 (the default) the right-hand side of the model is that
    of the model without its Hornbeck-Molnar rows, number for number, and
    the kinetic equation leaves them out (the cascade carries the levels'
    whole cross sections)."""
    m = Model(Params(eedf="two_t"))
    bare = Model(Params(eedf="two_t"))
    bare.compiled = [c for c in bare.compiled
                     if not c.row.rid.startswith("HM_")]
    y = _state(m)
    d, d0 = m.rhs(0.0, y, True), bare.rhs(0.0, y, True)
    live = list(range(m.ns)) + [m.i_ue, m.i_tg]
    assert np.array_equal(d[live], d0[live])
    for c in m.compiled:
        if c.row.rid.startswith("HM_"):
            assert c.factor == 0.0 and d[m.i_ev + c.j] == 0.0
    k = Model(Params())
    rows = [pr.rows[0] for pr in k.kin.procs]
    assert len(rows) == 52 and not any(r.startswith("HM_") for r in rows)
    sig = {pr.rows[0]: pr.p.sigma for pr in k.kin.procs}
    for c in k.compiled:
        if c.row.rid.startswith("C1c_"):
            p = xsec.find(*c.row.rate[1:4])
            assert np.array_equal(sig[c.row.rid], p.sigma)


def test_hm_at_one_ionizes_every_such_excitation():
    """hm_kr = 1: the twenty levels at or above E(Kr2+) make Kr2+ and an
    electron, their cascade carries nothing, the ledger closes, and the
    ions' formation energy the report books includes them. At the
    reference drive the 147-nm share moves by a few 1e-4 of itself: direct
    and associative ionization make well under a percent of the ions next
    to Penning's (2026-10-14)."""
    res0 = Model(Params(**FAST)).run()
    res1 = Model(Params(**FAST, rates=(("hm_kr", 1.0),))).run()
    hm = [r for r in res1.events if r.startswith("HM_")]
    assert len(hm) == 20
    made = sum(res1.events[r][-1] for r in hm)
    assert made > 0.0
    for r in hm:
        assert res1.events["C1c_" + r[3:]][-1] == 0.0
        assert res0.events[r][-1] == 0.0
    # the levels' total excitation is about the same either way
    lev0 = sum(res0.events["C1c_" + r[3:]][-1] for r in hm)
    assert made == approx(lev0, rel=0.05)
    a1 = accounting(res1)
    assert abs(a1["closure"]) < 1e-12
    assert a1["flows"]["   of which Hornbeck-Molnar (HM)"] == approx(
        made * e_form("Kr2+"), rel=1e-12)
    d0, d1 = row_dict(res0), row_dict(res1)
    assert d1["to_Krplus"] > d0["to_Krplus"]
    assert d1["emission_vuv"] == approx(d0["emission_vuv"], rel=2e-3)


def test_kinetic_run_with_hm_closes():
    """A short kinetic pulse at hm_kr = 1: the associative ionizations are
    in the Boltzmann equation, the ledger closes, the field carries the
    input and the two energy balances agree as at hm_kr = 0
    (test_jung2006)."""
    r = Model(Params(eedf="kinetic", pulse_s=40e-9, t_end_s=45e-9,
                     n_save=91, rates=(("hm_kr", 1.0),))).run()
    assert abs(accounting(r)["closure"]) < 1e-8
    assert sum(v[-1] for k, v in r.events.items() if k.startswith("HM_")) > 0
    kinds = {pr.kind for pr in r.model.kin.procs}
    assert kinds == {"excitation", "ionization", "assoc"}
    on = [e for e in r.kinetic if e["on"] and not e.get("skipped")
          and e["t"] >= 20e-9]
    ratio = [e["power_ratio"] for e in on if np.isfinite(e["power_ratio"])]
    assert 0.98 < min(ratio) <= max(ratio) < 1.03
    assert max(abs(e["balance"]) for e in on) < 2e-5


def test_associative_ionization_in_the_boltzmann_equation():
    """A driven state with every level at or above E(Kr2+) ionizing
    associatively: the growth frequency counts those electrons with the
    direct ones, and the solution's energy balance holds to the solver's
    tolerance: field power = thresholds paid + elastic loss + nu <eps>,
    less what the associative electrons bring, born in the first cell
    (its centre, du/2 = 0.025 eV here: 8.6e-5 of the field's power, which
    the ledger books at rest)."""
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
    assert born / gain == approx(8.6e-5, rel=0.05)


# ------------------------------------------------------------- the swarm
# LoKI-MC (the vendored Monte Carlo solver, LoKI-MC/Code/lokimc) on the
# same Biagi 8.97 file, rewritten in LoKI-MC's block form, 300 K, 1 Torr,
# 1e5 electrons, equal sharing (2026-10-14): the SST reduced Townsend
# coefficient [m^2], the mean energy [eV] and the ionizations the run
# counted after its steady state (29, 543, 3127 and 9521: a statistical
# error of 19, 4.3, 1.8 and 1.0 %)
LOKI_MC_897 = {20.0: (3.6057e-25, 4.5844, 29), 30.0: (7.1147e-24, 4.7536, 543),
               40.0: (3.3102e-23, 4.9022, 3127),
               50.0: (8.7447e-23, 5.0360, 9521)}


@pytest.mark.parametrize("en", sorted(LOKI_MC_897))
def test_swarm_against_monte_carlo(en):
    """The two-term solution on v0_2's cross sections against a Monte
    Carlo solution of the same: the ionization coefficient within 2 % (the
    two-term approximation) and two of the counting's standard errors
    (found: +3.7, -4.0, -3.5 and -2.5 % at 20, 30, 40 and 50 Td), the mean
    energy within 0.2 %. Bordage et al. (J. Phys. D 46, 334003 (2013))
    found the two-term alpha/N accurate in Kr against Monte Carlo too."""
    d = swarm.point(en)
    a_mc, m_mc, counted = LOKI_MC_897[en]
    assert d["alphaN"] * 1e-4 == approx(a_mc, rel=0.02 + 2.0 / counted ** 0.5)
    assert d["mean"] == approx(m_mc, rel=2e-3)


def test_swarm_hornbeck_molnar_band():
    """The band hm_kr in the swarm: at its top every excitation of the
    twenty levels ionizes, which raises the ionization coefficient 17-fold
    at 14 Td, where the tail falls by orders of magnitude between the
    levels (12.86-13.6 eV) and the ionization threshold (14.0 eV), and
    1.7-fold at 100 Td; the direct part stays within 3.4 % (the growth
    term). hm = 0 is v0_1's swarm, direct ionization alone."""
    lo, hi = swarm.point(14.12), swarm.point(14.12, hm=1.0)
    assert lo["alphaN"] == lo["alphaN_direct"]
    assert hi["alphaN"] / lo["alphaN"] == approx(16.84, rel=0.01)
    assert hi["alphaN_direct"] == approx(lo["alphaN_direct"], rel=1e-3)
    r100 = swarm.point(100.0, hm=1.0)["alphaN"] / swarm.point(100.0)["alphaN"]
    assert 1.6 < r100 < 1.8
    t100 = swarm.point(100.0, hm=1.0)
    assert t100["alphaN_direct"] == approx(swarm.point(100.0)["alphaN"],
                                           rel=0.034)
    with pytest.raises(ValueError):
        swarm.processes(1.5)
