"""The Maxwellian integrator against closed forms, and the data pins."""

import math
from pathlib import Path

import numpy as np
import pytest

from v0 import xsec

TE_CHECK = np.geomspace(0.05, 20.0, 25)


def mean_speed_cm_s(te):
    """<v> of a Maxwellian = sqrt(8 e Te / (pi m)) [cm/s]."""
    return math.sqrt(8.0 * xsec.QE * te / (math.pi * xsec.ME)) * 100.0


def _k_of(sigma_eps, te):
    """Run an arbitrary sigma(EPS) array through the module's integrator."""
    table = xsec.RateTable(xsec._mean(sigma_eps * xsec.V_CM_S
                                      * np.sqrt(xsec.EPS)), "test")
    return table(te)


def test_constant_sigma_gives_sigma_times_mean_speed():
    sig = 1e-16
    for te in TE_CHECK:
        k = _k_of(np.full_like(xsec.EPS, sig), te)
        assert k == pytest.approx(sig * mean_speed_cm_s(te), rel=2e-4)


def test_step_sigma_closed_form():
    """sigma = s0 above e0: k = s0 <v> (1 + e0/Te) exp(-e0/Te)."""
    s0, e0 = 1e-16, 10.0
    sig = np.where(xsec.EPS >= e0, s0, 0.0)
    for te in (0.8, 1.5, 3.0, 8.0):
        exact = s0 * mean_speed_cm_s(te) * (1 + e0 / te) * math.exp(-e0 / te)
        assert _k_of(sig, te) == pytest.approx(exact, rel=2e-3)


def test_energy_weighted_average_closed_form():
    """<sigma v eps> for constant sigma = 2 Te sigma <v> (the elastic
    energy-loss integral's kernel)."""
    sig = 1e-16
    for te in TE_CHECK:
        val = xsec.RateTable(xsec._mean(sig * xsec.V_CM_S * np.sqrt(xsec.EPS)
                                        * xsec.EPS), "t")(te)
        assert val == pytest.approx(2 * te * sig * mean_speed_cm_s(te),
                                    rel=2e-4)


def test_superelastic_is_detailed_balance():
    g = 1.0 / 12.0
    fwd = xsec.rate_table("Kr", "excitation", "Kr(s)(9.915eV)")
    sup = xsec.superelastic_table("Kr", "Kr(s)(9.915eV)", g)
    for te in (0.7, 1.0, 2.0, 5.0, 10.0):
        assert sup(te) == pytest.approx(g * fwd(te) * math.exp(9.915 / te),
                                        rel=5e-3)


def test_em_power_is_drude_for_constant_collision_frequency():
    """sigma_m = nu0 / (N v) makes nu(eps) = nu0: the Lorentz average must
    reduce to Drude, e^2 E^2 nu0 / (m (nu0^2 + w^2)) per electron."""
    n_gas, nu0, f, e_rms = 2.5e19, 3.0e11, 95e9, 1000.0
    w = 2 * math.pi * f

    def sigma_m(eps):
        return nu0 / (n_gas * xsec.V_CM_S * np.sqrt(eps))

    for te in (0.5, 1.5, 4.0):
        got = xsec.em_power_per_electron(te, n_gas, 0.0, 1e15, f, e_rms,
                                         sigma_m=sigma_m)
        drude = (xsec.QE ** 2 * (e_rms * 100) ** 2 / xsec.ME) * nu0 \
            / (nu0 ** 2 + w ** 2) / xsec.QE
        assert got == pytest.approx(drude, rel=2e-3)


def test_coulomb_log_nrl_form():
    assert xsec.coulomb_log(1e15, 2.0) == pytest.approx(
        23.0 - math.log(math.sqrt(1e15) * 2.0 ** -1.5))
    assert xsec.coulomb_log(1e16, 0.03) == 1.0      # strongly coupled floor


def test_biagi_blocks_and_thresholds():
    procs = xsec.processes("Kr")
    kinds = [(p.kind, p.threshold) for p in procs]
    assert kinds == [("elastic", 0.0), ("excitation", 9.915),
                     ("excitation", 11.30), ("excitation", 11.998),
                     ("excitation", 12.75), ("ionization", 13.996)]


def test_vendored_files_pinned_and_match_the_repo_copies():
    for gas, (name, md5, origin) in xsec.FILES.items():
        path = xsec.DATA_DIR / name
        assert xsec.md5_of(path) == md5
        rel = Path(origin.split(" ")[0])
        repo = next((b / rel for b in path.resolve().parents
                     if (b / rel).is_file()), None)
        if repo is None:
            pytest.skip(f"repo copy {rel} not found above LIGHT "
                        "(standalone LIGHT); pin checked, origin not")
        assert repo.read_bytes() == path.read_bytes(), repo


def test_elastic_energy_loss_fraction_is_2m_over_M():
    printed = xsec.find("Kr", "elastic").mass_ratio      # 6.5e-6 in the file
    assert printed == 6.5e-6
    assert xsec.two_m_over_M("Kr") == pytest.approx(2 * printed, rel=0.01)


# --- the captured energy and the two-temperature tail (2026-09-28) ---------
KR_TAIL = (("excitation", "Kr(s)(9.915eV)"), ("excitation", "Kr(p)(11.30eV)"),
           ("excitation", "Kr(d+p)(11.998eV)"), ("excitation", "Kr*(12.75eV)"),
           ("ionization", "Kr^+"))


def test_capture_energy_is_the_maxwellian_identity():
    """<eps sigma v>/<sigma v> = (3/2 + d ln k/d ln Te) Te, the identity
    model.py books the analytic rate laws by, here for F2 attachment on
    the table grid. F2 attaches below 1 eV, so at Tahoma's Te an attached
    electron carries off far less than (3/2)Te."""
    cap = xsec.capture_energy_table("F2", "attachment")
    k = xsec.rate_table("F2", "attachment")
    s = np.gradient(np.log(k.values), np.log(xsec.TE))
    sel = slice(100, 560)                       # 0.035 to 19 eV
    assert np.allclose(cap.values[sel], ((1.5 + s) * xsec.TE)[sel],
                       rtol=2e-3)
    assert cap(1.13) < 0.4 * 1.5 * 1.13


def test_only_kr_excitation_and_ionization_are_in_the_tail():
    assert xsec.e_break() == 9.915
    for kind, product in KR_TAIL:
        assert xsec.in_tail("Kr", kind, product)
    assert not xsec.in_tail("F2", "attachment")
    assert not xsec.in_tail("F2", "excitation", "F2(A1)(4.34eV)")


@pytest.mark.parametrize("kind,product", KR_TAIL)
def test_two_t_at_equal_temperatures_is_the_maxwellian(kind, product):
    """T_t = Te: the same discrete sums as the Maxwellian table."""
    two = xsec.two_t_table("Kr", kind, product)
    mx = xsec.rate_table("Kr", kind, product)
    for i in range(60, 600, 60):
        te = float(xsec.TE[i])
        assert two(te, te) == pytest.approx(mx(te), rel=1e-9)


def test_two_t_linear_ramp_closed_form():
    """sigma = s1 (eps - E_b) above the break:
    k = s1 V (E_b T_t^2 + 2 T_t^3) exp(-E_b/Te) / Z,
    Z = Te^1.5 gamma(3/2, E_b/Te) + exp(-E_b/Te) T_t^1.5 exp(E_b/T_t)
    Gamma(3/2, E_b/T_t), V = sqrt(2e/m) in cm/s per sqrt(eV)."""
    from scipy.special import gamma, gammainc, gammaincc
    eb, s1 = xsec.e_break(), 1e-17
    two = xsec.two_t_from_sigma(np.where(xsec.EPS >= eb,
                                         s1 * (xsec.EPS - eb), 0.0), "ramp")
    g32 = gamma(1.5)
    for te, tt in ((1.5, 0.4), (1.5, 1.5), (0.8, 0.3), (3.0, 1.0)):
        num = s1 * xsec.V_CM_S * (eb * tt ** 2 + 2 * tt ** 3) \
            * math.exp(-eb / te)
        z = (te ** 1.5 * g32 * gammainc(1.5, eb / te) + math.exp(-eb / te)
             * tt ** 1.5 * math.exp(eb / tt) * g32 * gammaincc(1.5, eb / tt))
        assert two(te, tt) == pytest.approx(num / z, rel=5e-3), (te, tt)


def test_a_cooler_tail_lowers_every_tail_rate():
    for kind, product in KR_TAIL:
        two = xsec.two_t_table("Kr", kind, product)
        k = [two(1.6, r * 1.6) for r in (1.0, 0.5, 0.25, 0.1)]
        assert all(a > b for a, b in zip(k, k[1:])), product


def test_tail_ratio_table_is_pinned_and_physical():
    """The kinetic table loads only at its pinned md5 and r lies in (0, 1].
    With F2 present, more electron-electron collisions (higher n_e/N) fill
    the tail. Pure Kr is not monotone at its lowest pumps: the injected
    electrons, slowed by nothing but elastic and e-e collisions, reach the
    break on their own (ML/ml/tail.py)."""
    t = xsec.tail_ratio_table()
    assert np.all(t.r > 0.0) and np.all(t.r <= 1.0)
    assert t.f2[0] == 0.0 and np.all(np.diff(t.f2) > 0)
    assert np.all(np.diff(t.r[:, 1:], axis=0) >= 0.0)
    assert "ML/ml/tail.py" in t.header
    # bilinear, held at the edges
    assert t(10.0 ** t.x[0] / 10.0, -1.0) == pytest.approx(t.r[0, 0])
    assert t(1.0, 1.0) == pytest.approx(t.r[-1, -1])
    xm, fm = 0.5 * (t.x[3] + t.x[4]), 0.5 * (t.f2[2] + t.f2[3])
    assert t(10.0 ** xm, t.f2[2]) == pytest.approx(0.5 * (t.r[3, 2]
                                                          + t.r[4, 2]))
    assert t(10.0 ** t.x[5], fm) == pytest.approx(0.5 * (t.r[5, 2]
                                                         + t.r[5, 3]))
