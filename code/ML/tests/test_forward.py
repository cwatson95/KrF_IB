"""Forward models against numbers the sources print."""

import pytest

from ml import data, forward, params


def test_yu_setser_eq18_returns_their_ratio():
    """Their constants (tau_B 7.0, tau_C 75.2 ns, K_eq 1.66) give their
    measured I_r(300 K) = 0.16 (J. Phys. Chem. 94, 2934, p.2941)."""
    c = params.light_now().with_(tau_B=7.0e-9, tau_C=75.2e-9, K_eq=1.66)
    ys = data.by_name("ys_ir300")[0]
    assert forward.predict(ys, c)[0] == pytest.approx(0.160, abs=0.001)


def test_gadomski_eq5_with_his_constants():
    """Gadomski's published set gives 4.45e7 1/s at 50 Torr Kr, 1 Torr F2
    (the value his fit leaves against Xu's quoted 5.0e7)."""
    gad = params.gadomski_1992()
    ds = data.Dataset("x", "slow_decay", "held_out", "test",
                      (({"kr_torr": 50.0, "f2_torr": 1.0}, 5.0e7),), 0.1)
    assert forward.predict(ds, gad)[0] == pytest.approx(4.45e7, rel=0.01)


def test_vector_round_trip():
    now = params.light_now()
    back = params.consts(params.vector(now))
    for n in params.NAMES:
        want = getattr(now, n) if getattr(now, n) is not None \
            else getattr(now, n[:-2])
        assert getattr(back, n) == pytest.approx(want, rel=1e-12)


def test_rederived_electron_quench_matches_light():
    """LIGHT's QE (6.6e-8 since 2026-09-28, 4.5e-8 before) is Trainor &
    Jacob re-read with its own B/C constants."""
    now = params.light_now()
    assert forward.with_rederived_k_e(now).k_e == pytest.approx(now.k_e,
                                                                rel=0.02)
    old = params.light_to_0927()
    assert forward.with_rederived_k_e(old).k_e == pytest.approx(4.5e-8,
                                                                rel=0.10)


def test_device_forward_runs_light_with_the_set():
    """Tahoma through the prototype drive on today's network and LIGHT's
    own closure (since 2026-09-28: the measured formation steps, the Kr2*
    route, the two-temperature tail, the exact captured energies and the
    absolute harpoon yield 0.90). LIGHT's adopted constants give g0 6.11
    %/cm (measured 6.39) with B share 0.53, and the constants it had until
    2026-09-27 give 3.17 (3.67 on that day's network and closure) with B
    share 0.22."""
    old = forward.run_device(forward.with_rederived_k_e(params.light_to_0927()))
    assert old["g0"] == pytest.approx(3.17, rel=0.03)
    assert 0.2 < old["b_share"] < 0.3
    now = forward.run_device(forward.with_rederived_k_e(params.light_now()))
    assert now["g0"] == pytest.approx(6.11, rel=0.03)
    assert 0.45 < now["b_share"] < 0.60
