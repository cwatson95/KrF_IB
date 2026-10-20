"""The closed form of the basic network's B/C block against LIGHT's own run."""

import importlib
import pytest

L = importlib.import_module("1Dv0.two_level")


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


@pytest.mark.parametrize("kr,f2", [(50.0, 1.0), (100.0, 1.0), (500.0, 1.0),
                                   (190.0, 2.0)])
def test_closed_form_equals_light_at_300_k(kr, f2):
    assert L.check_against_light(kr, f2) < 1e-5


def test_closed_form_equals_light_at_500_k():
    """K_eq(Tg) and the trimer's Tg law, through a seed run started hot."""
    import numpy as np
    _mod = importlib.import_module("1Dv0.model")
    Model, Params = _mod.Model, _mod.Params
    tot = 101.0
    p = Params(network="basic", pressure_atm=tot / 760, x_f2=1 / tot,
               tg0_K=500.0, power_W_cm3=0.0, pulse_s=0.0, t_end_s=300e-9,
               ne_mode="free", ne_cm3=1.0, n_save=3001,
               init_extra_cm3=(("KrF(B)", 1e12),))
    r = Model(p).run()
    n_tot = p.pressure_atm * 101325 / (1.380649e-23 * 500.0) * 1e-6
    Ba, Ca = L.impulse(L.Consts.from_light(), r.t, B0=r.n["KrF(B)"][0],
                       kr=n_tot * 100 / tot, f2=n_tot / tot, T=500.0)
    sel = r.t < 100e-9
    b0 = r.n["KrF(B)"][0]
    assert np.max(abs(r.n["KrF(B)"][sel] - Ba[sel])) / b0 < 1e-5
    assert np.max(abs(r.n["KrF(C)"][sel] - Ca[sel])) / b0 < 1e-5


def test_from_light_reads_the_network():
    c = L.Consts.from_light()
    assert c.k_BC == 1.05e-11 and c.K_eq == 1.10
    assert c.tau_C == approx(50e-9)
    assert (c.k3_Kr, c.k3_Kr_C) == (1.90e-30, 3.79e-31)
    assert (c.k_F2, c.k_F2_C) == (4.75e-10, 3.04e-10)
    assert c.gap_cm1 == approx(120.0, rel=1e-6)          # Yu & Setser 1990
