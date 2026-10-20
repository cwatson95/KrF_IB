"""The electron side: the e-e operator, the solver against LoKI-MC, the
beam-sustained EEDF and LIGHT with it."""
import math

import numpy as np
import pytest

from ml import boltzmann  # noqa: F401
from ml import eedf, electron, forward, score, swarm
import boltz_twin as BT
import ee
import lx_lib as L


def _grid(N=400, emax=20.0):
    du = emax / N
    return du, (np.arange(N) + 0.5) * du, np.arange(1, N) * du


@pytest.mark.parametrize("te", [0.5, 1.13, 3.0])
def test_ee_flux_vanishes_for_a_maxwellian(te):
    """Coulomb collisions leave a Maxwellian of any temperature alone."""
    du, ec, eb = _grid()
    F = ee.maxwellian(ec, du, te)
    We, De = ee.ee_coeffs(F, ec, eb, du, 1.0)
    z = We * du / De
    flux = De / du * (BT.bern(-z) * F[:-1] - BT.bern(z) * F[1:])
    assert np.max(np.abs(flux[:200]) / (De / du * F[:-1])[:200]) < 2e-3


def test_ee_slows_a_fast_electron_at_the_spitzer_rate():
    """dE/dt = -n_e e^4 lnL / (4 pi eps0^2 m v) for v >> v_thermal."""
    ne, N, te = 4.8e20, 2.14e25, 1.13
    du, ec, eb = _grid(2000, 60.0)
    F = ee.maxwellian(ec, du, te)
    lnl = ee.coulomb_log(ne, 1.5 * te)
    alpha = ne / N * ee.C_EE * lnl * L.GAMMA
    We, _ = ee.ee_coeffs(F, ec, eb, du, alpha)
    i = int(30.0 / du)
    v = math.sqrt(2 * L.QE * eb[i] / L.ME)
    spitzer = -ne * L.QE ** 3 * lnl / (4 * math.pi * ee.EPS0 ** 2 * L.ME * v)
    assert N * We[i] / math.sqrt(eb[i]) == pytest.approx(spitzer, rel=2e-3)


def test_ee_makes_a_field_driven_eedf_maxwellian():
    """Kr/F2 at 3 Td: far from Maxwellian as a swarm, Maxwellian at a high
    ionization degree, and the e-e term moves (almost) no energy."""
    kr, f2 = eedf.gases(300.0)
    mix = [(kr, 0.995), (f2, 0.005)]

    def dev(neN):
        sol = ee.solve_field(mix, 3.0, 300.0, 400, 20.0, ne_N=neN,
                             ne_m3=4.8e20)
        sw = BT.swarm(mix, sol)
        Fm = ee.maxwellian(sol["ec"], sol["du"], 2 / 3 * sw["mean"])
        d = float(np.sum(np.sqrt(sol["ec"]) * np.abs(sol["F"] - Fm))
                  * sol["du"])
        return d, sol, sw
    assert dev(0.0)[0] > 0.3
    d, sol, sw = dev(1e-3)
    assert d < 0.03
    eb = np.arange(1, 400) * sol["du"]
    alpha = 1e-3 * ee.C_EE * ee.coulomb_log(4.8e20, sw["mean"]) * L.GAMMA
    p_ee = ee.ee_power(sol["F"], sol["ec"], eb, sol["du"], alpha)
    assert abs(p_ee) < 0.05 * 3.0 * BT.TD * sw["w"]


def test_two_term_matches_lokimc_in_helium_fluorine():
    """He + 1 % F2 at 10 Td, 300 K. LoKI-MC 2026-09-27 (2e4 electrons):
    <eps> 5.033 eV, flux W 2.091e4 m/s, k_att per F2 5.630e-16 m3/s."""
    he, f2 = swarm.gas("He", 300.0), swarm.gas("F2", 300.0)
    mix = [(he, 0.99), (f2, 0.01)]
    sol = BT.auto_solve(mix, 10.0, 300.0, 400, 20.0)
    sw = BT.swarm(mix, sol)
    k = [kk for g, d, kk in sw["rates"] if "Attachment" in d][0]
    assert sw["mean"] == pytest.approx(5.033, rel=0.01)
    assert sw["w"] == pytest.approx(2.091e4, rel=0.01)
    assert k == pytest.approx(5.630e-16, rel=0.02)


@pytest.fixture(scope="module")
def tahoma():
    """The fit set and LIGHT's Tahoma state on its Maxwellian closure, the
    state the kinetic loop starts from."""
    c = score.candidate_sets()["fit (train set)"]
    return c, forward.run_device(c, closure=forward.MAXWELL)


@pytest.fixture(scope="module")
def loop(tahoma):
    return eedf.self_consistent(tahoma[0])


def test_beam_eedf_bulk_is_maxwellian_and_tail_is_empty(tahoma):
    """At LIGHT's Tahoma state every rate below 9 eV is within 5 % of the
    Maxwellian with the same mean energy, and the Kr rates are cut."""
    b = eedf.beam_eedf(tahoma[1])
    for r in ("ATT", "VIB1", "VIB2", "VIB3", "VIB4", "DIS_a3", "DIS_A1"):
        assert b["f"][r] == pytest.approx(1.0, abs=0.05), r
    assert b["f"]["G1"] < 0.2
    assert b["f"]["C1a"] < 0.02
    assert b["f"]["C2"] < 1e-4
    # electrons balance: pairs = attachment + recombination
    n = tahoma[1]["n"]
    att = b["k"]["ATT"] * b["ne"] * n["F2"]
    rec = tahoma[1]["events"]["R11p"] / tahoma[1]["ne"] * b["ne"]
    assert att + rec == pytest.approx(tahoma[1]["pairs"], rel=0.01)


def test_beam_eedf_without_ee_is_not_maxwellian(tahoma):
    """Drop e-e (LoKI-MC's physics) and the bulk departs: attachment falls
    below half its Maxwellian value."""
    b = eedf.beam_eedf(tahoma[1], ee_on=False)
    assert b["f"]["ATT"] < 0.5


def test_light_and_boltzmann_agree_at_the_fixed_point(loop):
    """The loop's regression values (2026-09-28, today's network) and its
    self-check: with the exact capture energies LIGHT's Te and n_e land
    within 2 % and 8 % of the Boltzmann solution's."""
    dev, b, hist = loop
    assert dev["g0"] == pytest.approx(6.07, abs=0.05)
    assert dev["formation"] == pytest.approx(0.317, abs=0.003)
    assert dev["te"] == pytest.approx(2 / 3 * b["mean"], rel=0.02)
    assert dev["ne"] == pytest.approx(b["ne"], rel=0.08)


def test_an_attaching_electron_carries_off_far_less_than_1p5_te(tahoma):
    """F2 attaches below 1 eV: the energy an attached electron removes is
    about 0.5 eV at Tahoma, a quarter of the (3/2)Te LIGHT booked until
    2026-09-28, and the Boltzmann bulk gives the Maxwellian's value, which
    LIGHT's capture table now books."""
    b = eedf.beam_eedf(tahoma[1])
    te = 2 / 3 * b["mean"]
    kr, f2 = eedf.gases(tahoma[1]["tg"])
    a = eedf.collision(kr, f2, "F2", 0.0)
    du = 0.002
    ec = (np.arange(15000) + 0.5) * du
    w = ec * L.interp(a, ec) * ee.maxwellian(ec, du, te)
    e_max = float(np.sum(ec * w) / np.sum(w))
    assert b["e_att"] == pytest.approx(e_max, rel=0.05)
    assert b["e_att"] < 0.35 * 1.5 * te
    from v0 import xsec
    assert xsec.capture_energy_table("F2", "attachment")(te) == \
        pytest.approx(e_max, rel=0.01)


def test_att_table_unscaled_is_lights_own():
    from v0 import xsec
    ref = xsec.rate_table("F2", "attachment")
    new = electron.att_table(1.0, 1.0)
    assert np.allclose(new.values, ref.values, rtol=1e-12)


def test_swarm_registry_roles_and_linearity():
    """Every swarm set has a role, no device data is among them, and for
    trace F2 in N2 the attachment rate is linear in the cross section."""
    sets = swarm.load()
    assert {s.role for s in sets} <= {"train", "held_out", "check"}
    assert all(np.all(s.y > 0) and np.all(s.sigma > 0) for s in sets)
    n2, f2 = swarm.gas("N2", 300.0), swarm.gas("F2", 300.0)
    f2m = swarm.f2_modified(300.0, 0.7, 1.5)
    x = swarm.X_F2_TJ
    sol, _ = swarm._solve([(n2, 1 - x), (f2, x)], 5.0, 300.0)
    solm, _ = swarm._solve([(n2, 1 - x), (f2m, x)], 5.0, 300.0)
    c = next(q for q in f2m["colls"] if q["type"] == "Attachment")
    lin = np.sum(sol["ec"] * L.interp(c, sol["ec"]) * sol["F"]) * sol["du"]
    exact = (np.sum(solm["ec"] * L.interp(c, solm["ec"]) * solm["F"])
             * solm["du"])
    assert lin == pytest.approx(exact, rel=0.01)


def test_light_closure_lands_on_the_kinetic_fixed_point(tahoma, loop):
    """LIGHT on its own closure (the two-temperature tail from the kinetic
    table and the exact captured energies, no factors from outside) against
    the loop's fixed point (its Maxwellian closure with every electron row
    scaled by the Boltzmann factors): Te within 2 %, n_e within 8 %, g0 and
    the formation efficiency within 2 %."""
    ref = loop[0]
    dev = forward.run_device(tahoma[0])
    assert dev["te"] == pytest.approx(ref["te"], rel=0.02)
    assert dev["ne"] == pytest.approx(ref["ne"], rel=0.08)
    assert dev["g0"] == pytest.approx(ref["g0"], rel=0.02)
    assert dev["formation"] == pytest.approx(ref["formation"], rel=0.02)
