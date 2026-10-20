"""The Maxwellian integrator against closed forms, and the data pins."""

import math
from pathlib import Path

import numpy as np
import pytest

from v0_5 import xsec

HERE = Path(__file__).resolve()

def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)

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
        assert k == approx(sig * mean_speed_cm_s(te), rel=2e-4)


def test_step_sigma_closed_form():
    """sigma = s0 above e0: k = s0 <v> (1 + e0/Te) exp(-e0/Te). A step
    cannot sit between grid points, so the quadrature is good to about half
    a cell of the Boltzmann factor, d eps / (2 Te) (0.7 % at 0.8 eV, 0.04 %
    at 8 eV for a step at 10 eV). v0's version asked 0.2 % but its
    pytest.approx also passed anything below 1e-12, so it never checked."""
    s0 = 1e-16
    for e0 in (10.0, 9.915):
        sig = np.where(xsec.EPS >= e0, s0, 0.0)
        i = int(np.searchsorted(xsec.EPS, e0))
        cell = xsec.EPS[i] - xsec.EPS[i - 1]
        for te in (0.8, 1.5, 3.0, 8.0):
            exact = (s0 * mean_speed_cm_s(te) * (1 + e0 / te)
                     * math.exp(-e0 / te))
            assert _k_of(sig, te) == approx(exact, rel=cell / (2 * te))


def test_energy_weighted_average_closed_form():
    """<sigma v eps> for constant sigma = 2 Te sigma <v> (the elastic
    energy-loss integral's kernel)."""
    sig = 1e-16
    for te in TE_CHECK:
        val = xsec.RateTable(xsec._mean(sig * xsec.V_CM_S * np.sqrt(xsec.EPS)
                                        * xsec.EPS), "t")(te)
        assert val == approx(2 * te * sig * mean_speed_cm_s(te),
                                    rel=2e-4)


def test_superelastic_is_detailed_balance():
    """On one level (1s2) and on the 5s manifold's summed cross section at
    its 1s5 energy (v0_2's SUPER): the Klein-Rosseland table is detailed
    balance over the cross section's own threshold."""
    g = 1.0 / 12.0
    for product, u in (("Kr(1S2)", 10.6436), (xsec.LUMPED_5S, 9.9152)):
        fwd = xsec.rate_table("Kr", "excitation", product)
        sup = xsec.superelastic_table("Kr", product, g)
        for te in (0.7, 1.0, 2.0, 5.0, 10.0):
            assert sup(te) == approx(g * fwd(te) * math.exp(u / te),
                                     rel=5e-3), (product, te)


def test_the_lumped_5s_manifold_is_the_four_levels():
    """v0_2: SUPER's cross section is the sum of the 1s5, 1s4, 1s3 and 1s2
    blocks at the 1s5 threshold, and it is not one of the file's processes
    (the kinetic equation's momentum transfer must not count it twice)."""
    lump = xsec.find("Kr", "excitation", xsec.LUMPED_5S)
    assert lump.threshold == 9.9152
    assert lump not in xsec.processes("Kr")
    five_s = [xsec.find("Kr", "excitation", f"Kr(1S{i})") for i in (5, 4, 3, 2)]
    for e in (9.95, 10.3, 10.6, 10.7, 12.0, 20.0, 100.0):
        assert float(lump.sigma_at(e)) == approx(
            sum(float(q.sigma_at(e)) for q in five_s), rel=1e-12)


def test_biagi_blocks_and_thresholds():
    """Biagi 8.97: the elastic momentum transfer, 51 excitation levels in
    rising threshold from 1s5 (9.9152 eV) to 'HIGH' (13.6 eV), and the
    ionization at 13.9996 eV. v0_1's Biagi-v7.1 had four groups."""
    procs = xsec.processes("Kr")
    kinds = [p.kind for p in procs]
    assert kinds == ["elastic"] + ["excitation"] * 51 + ["ionization"]
    exc = [p for p in procs if p.kind == "excitation"]
    thr = [p.threshold for p in exc]
    assert thr == sorted(thr)
    assert [(p.product, p.threshold) for p in exc[:5]] == [
        ("Kr(1S5)", 9.9152), ("Kr(1S4)", 10.0324), ("Kr(1S3)", 10.5624),
        ("Kr(1S2)", 10.6436), ("Kr(2P10)", 11.3035)]
    assert (exc[-1].product, exc[-1].threshold) == ("Kr(HIGH)", 13.6)
    assert procs[-1].threshold == 13.9996
    assert len({p.product for p in exc}) == 51


def test_vendored_files_pinned_and_match_the_repo_copies():
    """Krypton's file alone is vendored (LIGHT's F2 file is not carried);
    v0_3 adds the trained factors, pinned as the file is."""
    assert list(xsec.FILES) == ["Kr"] and list(xsec.MASS_AMU) == ["Kr"]
    assert sorted(p.name for p in xsec.DATA_DIR.iterdir()) == [
        "Kr_Biagi897_LXCat.txt", "kr_swarm_trained.json",
        "kr_tail_field.txt", "kr_tail_ratio.txt"]
    name, md5 = xsec.TRAINED_FILE
    assert xsec.md5_of(xsec.DATA_DIR / name) == md5
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
    assert xsec.two_m_over_M("Kr") == approx(2 * printed, rel=0.01)


# --- the captured energy and the two-temperature tail (2026-09-28) ---------
KR_TAIL = tuple([(p.kind, p.product) for p in xsec.processes("Kr")
                 if p.kind != "elastic"] + [("excitation", xsec.LUMPED_5S)])


def test_capture_energy_is_the_maxwellian_identity():
    """<eps sigma v>/<sigma v> = (3/2 + d ln k/d ln Te) Te, the identity
    model.py books the analytic rate laws by (R11p's (3/2 + n) Te), here on
    the table grid for krypton's elastic and 5s excitation blocks. (LIGHT
    held it for F2 attachment, which no pure-Kr row has.)"""
    for kind, product, sel in (("elastic", None, slice(100, 560)),
                               ("excitation", "Kr(1S5)",
                                slice(300, 560))):
        cap = xsec.capture_energy_table("Kr", kind, product)
        k = xsec.rate_table("Kr", kind, product)
        s = np.gradient(np.log(k.values), np.log(xsec.TE))
        assert np.allclose(cap.values[sel], ((1.5 + s) * xsec.TE)[sel],
                           rtol=2e-3), kind
    # the 5s excitation takes electrons above its threshold
    assert xsec.capture_energy_table("Kr", "excitation",
                                     "Kr(1S5)")(1.13) > 9.9152


def test_only_kr_excitation_and_ionization_are_in_the_tail():
    assert xsec.e_break() == 9.9152
    for kind, product in KR_TAIL:
        assert xsec.in_tail("Kr", kind, product)
    assert not xsec.in_tail("Kr", "elastic")


@pytest.mark.parametrize("kind,product", KR_TAIL)
def test_two_t_at_equal_temperatures_is_the_maxwellian(kind, product):
    """T_t = Te: the same discrete sums as the Maxwellian table."""
    two = xsec.two_t_table("Kr", kind, product)
    mx = xsec.rate_table("Kr", kind, product)
    for i in range(60, 600, 60):
        te = float(xsec.TE[i])
        assert two(te, te) == approx(mx(te), rel=1e-9)


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
        assert two(te, tt) == approx(num / z, rel=5e-3), (te, tt)


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
    assert t(10.0 ** t.x[0] / 10.0, -1.0) == approx(t.r[0, 0])
    assert t(1.0, 1.0) == approx(t.r[-1, -1])
    xm, fm = 0.5 * (t.x[3] + t.x[4]), 0.5 * (t.f2[2] + t.f2[3])
    assert t(10.0 ** xm, t.f2[2]) == approx(0.5 * (t.r[3, 2]
                                                          + t.r[4, 2]))
    assert t(10.0 ** t.x[5], fm) == approx(0.5 * (t.r[5, 2]
                                                         + t.r[5, 3]))


def test_field_tail_table_is_pinned_and_physical():
    """The field-driven table (tail_field.py) loads only at its pinned
    md5. v0_5: r(n_e/N, Te, x_s), its x_s = 0 block v0_4's table to the
    printed digit (v0_5's solver with no excited atom is v0_4's); at the
    reference point's state without excited atoms (n_e/N 1.2e-5, Te 2.3
    eV) it gives 0.10 to 0.13, where v0_3's beam table gave 0.19; there
    more electron-electron collisions (higher n_e/N) fill the tail, and
    above 4 eV the tail is close to the Maxwellian's whatever n_e/N and
    x_s. With excited atoms the superelastic electrons heat the tail: at
    x_s = 1e-3 and 2.5 eV it is more than twice the atoms-free one, and
    some low-Te points are hotter than the form follows (held at the top,
    3). The reader is trilinear in (log10 n_e/N, ln Te, log10(x_s +
    1e-7)) on ln r, held at its edges."""
    t = xsec.tail_field_table()
    name, md5 = xsec.TAIL_FIELD_FILE
    assert xsec.md5_of(xsec.DATA_DIR / name) == md5
    assert t.r.shape == (14, 17, 7) and t.xs[0] == 0.0
    assert np.all(t.r >= 1e-3) and np.all(t.r <= 3.0)
    assert "tail_field.py" in t.header
    v04 = HERE.parents[2] / "v0_4" / "data" / name
    if v04.is_file():
        assert np.array_equal(xsec.read_tail_field(v04).r[:, :, 0],
                              t.r[:, :, 0])
    z = t.r[:, :, 0]
    assert 0.10 < t(1.2e-5, 2.3) < 0.13 and t(1.2e-5, 2.3, 0.0) == t(
        1.2e-5, 2.3)
    assert xsec.tail_ratio_table()(1.2e-5, 0.0) > 0.18
    assert np.all(np.diff(z[-6:, :], axis=0) > 0.0)     # n_e/N >= 1e-5
    assert np.all(t.r[:, t.te >= 4.0, :] > 0.55)
    assert t(1e-5, 2.5, 1e-3) > 2.0 * t(1e-5, 2.5, 0.0)
    assert np.sum(t.r[:, :, 0] == 3.0) == 0 and np.sum(t.r == 3.0) > 100
    assert t(1e-12, 0.1) == approx(t.r[0, 0, 0])
    assert t(1.0, 100.0, 1.0) == approx(t.r[-1, -1, -1])
    i, j, k = 9, 7, 3
    xm = 0.5 * (t.x[i] + t.x[i + 1])
    assert t(10.0 ** xm, t.te[j], t.xs[k]) == approx(math.sqrt(
        t.r[i, j, k] * t.r[i + 1, j, k]))
    tm = math.sqrt(t.te[j] * t.te[j + 1])
    assert t(10.0 ** t.x[i], tm, t.xs[k]) == approx(math.sqrt(
        t.r[i, j, k] * t.r[i, j + 1, k]))
    zm = 0.5 * (t.zs[k] + t.zs[k + 1])
    assert t(10.0 ** t.x[i], t.te[j], 10.0 ** zm - t.xs0) == approx(
        math.sqrt(t.r[i, j, k] * t.r[i, j, k + 1]))


def test_field_tail_table_regenerates():
    """Points of the table solved again by its generator (the run's
    kinetic solver with its population processes, cells and tolerance)
    give the file's values to its five decimals, with and without excited
    atoms, and a point hotter than the form follows is held at the top."""
    from v0_5 import tail_field
    t = xsec.tail_field_table()
    for i, j, k in ((8, 8, 0), (9, 9, 0), (11, 10, 4), (9, 8, 6)):
        r, held, _ = tail_field.build(x_grid=(float(t.x[i]),),
                                      te_grid=(float(t.te[j]),),
                                      xs_grid=(float(t.xs[k]),))
        assert not held[0, 0, 0]
        assert r[0, 0, 0] == approx(t.r[i, j, k], rel=0.0, abs=6e-6), \
            (i, j, k)
    r, held, info = tail_field.build(x_grid=(-6.0,), te_grid=(0.5,),
                                     xs_grid=(1e-5,))
    assert held[0, 0, 0] and r[0, 0, 0] == tail_field.R_RANGE[1]
    assert info[0]["loss_kin"] > 100.0 * info[0]["loss_maxw"]
