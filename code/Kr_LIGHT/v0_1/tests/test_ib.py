"""LIGHT_Kr v0_1's ib.py (LIGHT v0_4's module, on krypton alone) against
the IB project's Turnbull oracles.

Run from the Kinetic_Modeling root:
    /usr/bin/python3 -m pytest LIGHT_Kr/v0_1/tests/test_ib.py -q

Oracles (read-only, test-time only; the module re-implements every formula):
  projects/IB/quant_IB/turnbull.py      Turnbull et al. [T24] Eqs. (2)-(11)
  projects/IB/quant_IB/qib_kinetics.py  per-energy L(eps) and screening
  projects/IB/nonlinear/nlib_kinetics.py  sigma_c = pi (e/(4 pi eps0 eps))^2 L
  LIGHT/v0/xsec.py                      the Lorentz-gas kernel this replaces
                                        (LIGHT/ beside LIGHT_Kr/)
Tests needing an absent oracle skip with the import error as reason.
The IB project is looked for at <Agent>/projects/IB; override with the
environment variable IB_PROJECT_DIR.
"""

import importlib
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve()
PKG_DIR = HERE.parents[1]                      # LIGHT_Kr/v0_1
LIGHT = HERE.parents[2]                        # LIGHT_Kr
IB_DIR = Path(os.environ.get("IB_PROJECT_DIR", LIGHT.parents[1] / "IB"))
for p in (LIGHT, IB_DIR / "quant_IB", IB_DIR / "nonlinear"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
# the v0 oracle lives in LIGHT/ beside LIGHT_Kr/ (LIGHT_Kr carries no v0);
# appended, so it never shadows LIGHT_Kr's own packages
SIBLING = LIGHT.parent / "LIGHT"
if (SIBLING / "v0" / "xsec.py").is_file() and str(SIBLING) not in sys.path:
    sys.path.append(str(SIBLING))

ib = importlib.import_module(f"{PKG_DIR.name}.ib")


def _try_import(name):
    try:
        return importlib.import_module(name), None
    except Exception as exc:                   # noqa: BLE001  any failure = oracle absent
        return None, f"oracle {name} not importable: {exc!r}"


tb, TB_REASON = _try_import("turnbull")
qk, QK_REASON = _try_import("qib_kinetics")    # pulls in nlib_kinetics (NEON scripts)
nk, NK_REASON = _try_import("nlib_kinetics")
v0x, V0_REASON = _try_import("v0.xsec")

needs_tb = pytest.mark.skipif(tb is None, reason=TB_REASON or "")
needs_qk = pytest.mark.skipif(qk is None or tb is None, reason=(QK_REASON or TB_REASON or ""))
needs_nk = pytest.mark.skipif(nk is None, reason=NK_REASON or "")
needs_v0 = pytest.mark.skipif(v0x is None, reason=V0_REASON or "")

THZ = 1e12
TWO_PI = 2 * math.pi
N_1ATM_300K = 101325.0 / (ib.KB * 300.0)       # m^-3


def rel(a, b):
    return abs(a / b - 1.0)


# --------------------------------------------------------------------------
# 1. Turnbull's Maxwellian nu_ei with the Sommerfeld/Pradler logarithm
# --------------------------------------------------------------------------
@needs_tb
@pytest.mark.parametrize("te", [0.5, 1.5, 3.0])
@pytest.mark.parametrize("f_thz", [0.1, 1.0, 10.0, 30.0])
@pytest.mark.parametrize("n", [1e20, 1e21, 1e22])
@pytest.mark.parametrize("average", ["eq8", "lower", "net"])
def test_nu_ei_maxwell_matches_turnbull(te, f_thz, n, average):
    w = TWO_PI * f_thz * THZ
    ours = ib.nu_ei_maxwell(te, w, n, n, Z=1.0, average=average)
    oracle = tb.nu_ei(te, [(1.0, n, tb.lnL_sommerfeld(te, w, 1.0, average))])
    assert rel(ours, oracle) < 1e-12


@needs_tb
def test_nu_ei_prefactor_is_turnbulls_2p91e12():
    assert rel(ib.NU_EI_PREF, tb.NU_EI_PREF) < 1e-14
    assert rel(ib.NU_EI_PREF, 2.91e-12) < 2e-3


@needs_tb
def test_thermal_gff_averages_differ_only_at_order_xi():
    te, w = 1.5, TWO_PI * 10 * THZ
    xi = ib.HBAR * w / (ib.QE * te)
    eq8 = ib.thermal_gff(te, w, average="eq8")
    assert rel(ib.thermal_gff(te, w, average="lower"), eq8 * math.exp(xi)) < 1e-14
    assert rel(ib.thermal_gff(te, w, average="net"), eq8 * math.expm1(xi) / xi) < 1e-14
    assert rel(eq8, tb.thermal_gff(te, w)) < 1e-12


# --------------------------------------------------------------------------
# 2. screening and Langdon factors
# --------------------------------------------------------------------------
@needs_tb
@pytest.mark.parametrize("te", [0.5, 1.5, 3.0])
@pytest.mark.parametrize("f_thz", [0.095, 1.0, 10.0])
@pytest.mark.parametrize("ni", [1e20, 1e21, 1e22])
def test_f_screen_matches_turnbull_at_ion_sphere_radius(te, f_thz, ni):
    w = TWO_PI * f_thz * THZ
    assert rel(ib.a_ws(ni), tb.a_ws(ni)) < 1e-14
    assert rel(ib.f_screen(te, w, ni), tb.f_sc(te, w, tb.a_ws(ni))) < 1e-13
    assert 0.0 < ib.f_screen(te, w, ni) <= 1.0


@needs_tb
@pytest.mark.parametrize("te", [0.5, 1.3, 3.0])
@pytest.mark.parametrize("f_thz", [0.095, 1.0, 10.0])
@pytest.mark.parametrize("e_rms", [1e4, 1.5e5, 2e6])
@pytest.mark.parametrize("z_eff", [1.0, 2.0])
def test_f_langdon_matches_turnbull(te, f_thz, e_rms, z_eff):
    w = TWO_PI * f_thz * THZ
    e0 = math.sqrt(2.0) * e_rms                 # amplitude from rms
    a_l = tb.alpha_L(e0, w, te, z_eff)
    assert rel(ib.alpha_langdon(e0, w, te, z_eff), a_l) < 1e-14
    assert rel(ib.matte_m(a_l), tb.matte_m(a_l)) < 1e-14
    new = ib.f_langdon(e0, w, te, Z_eff=z_eff, version="new")
    orig = ib.f_langdon(e0, w, te, Z_eff=z_eff, version="orig")
    assert rel(new, tb.fL_new(tb.matte_m(a_l), te, w)) < 1e-12
    assert rel(orig, tb.fL_orig(a_l)) < 1e-14
    assert 0.0 < new <= 1.0 + 1e-12 and 0.0 < orig <= 1.0


def test_f_langdon_is_one_without_a_field():
    assert ib.matte_m(0.0) == 2.0
    assert ib.fL_orig(0.0) == 1.0
    assert rel(ib.f_langdon(0.0, TWO_PI * THZ, 1.5), 1.0) < 1e-12
    with pytest.raises(ValueError):
        ib.f_langdon(1e5, TWO_PI * THZ, 1.5, version="other")


# --------------------------------------------------------------------------
# 3. per-energy logarithm and its Maxwellian average
# --------------------------------------------------------------------------
EPS_GRID = np.geomspace(1e-4, 100.0, 300)


@needs_qk
@pytest.mark.parametrize("f_thz", [0.1, 1.0, 10.0, 30.0])
@pytest.mark.parametrize("ni", [1e20, 1e22])
def test_lnL_per_energy_matches_qib_kinetics(f_thz, ni):
    w = TWO_PI * f_thz * THZ
    hw = tb.HBAR * w / tb.QE
    above = EPS_GRID >= 0.5 * hw
    l_sc = tb.a_ws(ni)
    oracle = np.where(above, qk.L_sommerfeld(EPS_GRID, w) + qk.screen_add(EPS_GRID, w, l_sc), 0.0)
    ours = ib.lnL_per_energy(EPS_GRID, w, ni, ni, screening="turnbull", floor=False)
    assert ours.shape == EPS_GRID.shape
    np.testing.assert_allclose(ours, oracle, rtol=1e-10, atol=0.0)
    # zeros below the threshold, and only there
    assert np.all(ours[~above] == 0.0)
    if f_thz >= 10.0:
        assert np.sum(~above) > 0, "the grid should reach below hbar w/2 at this frequency"
    assert np.all(ours[above] != 0.0)
    # the floor clips the negative values the screening term makes for slow electrons
    floored = ib.lnL_per_energy(EPS_GRID, w, ni, ni, screening="turnbull")
    np.testing.assert_array_equal(floored, np.maximum(ours, 0.0))
    # unscreened and Debye-screened variants
    np.testing.assert_allclose(ib.lnL_per_energy(EPS_GRID, w, ni, ni, screening="none"),
                               np.where(above, qk.L_sommerfeld(EPS_GRID, w), 0.0), rtol=1e-10, atol=0.0)
    te = 1.5
    l_d = qk.debye_e(ni, te)
    assert rel(ib.debye_length(ni, te), l_d) < 1e-14
    deb = np.where(above, qk.L_sommerfeld(EPS_GRID, w) + qk.screen_add(EPS_GRID, w, l_d), 0.0)
    np.testing.assert_allclose(ib.lnL_per_energy(EPS_GRID, w, ni, ni, screening="debye", te_eV=te, floor=False),
                               deb, rtol=1e-10, atol=0.0)


@needs_qk
def test_sigma_ei_matches_qib_coulomb_model_turnbull():
    """qib's coulomb_model("turnbull") screens at a_i of CAT_PER_E n_e cations and
    floors the logarithm; sigma_ei with ni = CAT_PER_E ne is the same function."""
    w, ne = TWO_PI * THZ, 1e21
    oracle = qk.coulomb_model("turnbull")(EPS_GRID, w, ne)
    ours = ib.sigma_ei(EPS_GRID, w, ne, qk.CAT_PER_E * ne)
    np.testing.assert_allclose(ours, oracle, rtol=1e-10, atol=0.0)


def test_lnL_per_energy_argument_errors():
    with pytest.raises(ValueError):
        ib.lnL_per_energy(1.0, 0.0, 1e21, 1e21)
    with pytest.raises(ValueError):
        ib.lnL_per_energy(1.0, TWO_PI * THZ, 1e21, 1e21, screening="debye")
    with pytest.raises(ValueError):
        ib.lnL_per_energy(1.0, TWO_PI * THZ, 1e21, 1e21, screening="other")


@needs_tb
def test_lnL_maxwell_relation_to_turnbull_thermal_average():
    """lnL_M = int L(eps) e^(-eps/T) deps/T (qib's absorption_parts) against
    Turnbull's thermal logarithm at Te = 1.5 eV, f = 10 THz.

    Exact relation (derived in ib.lnL_maxwell): L(eps) sits at the upper
    energy E_i = eps + hbar w/2, so e^(-eps/T) = e^(xi/2) e^(-E_i/T) and
        lnL_M = e^(xi/2) lnL_eq8 = e^(-xi/2) lnL_lower,
    neither e^xi nor (1 - e^-xi)/xi. Against the "net" average the ratio is
    e^(xi/2) xi/(e^xi - 1) = 1 - xi^2/24 + O(xi^4)."""
    te, w, n = 1.5, TWO_PI * 10 * THZ, 1e21
    xi = ib.HBAR * w / (ib.QE * te)
    ours = ib.lnL_maxwell(te, w, n, n, screening="none")
    lower = tb.lnL_sommerfeld(te, w, 1.0, "lower")
    net = tb.lnL_sommerfeld(te, w, 1.0, "net")
    eq8 = tb.lnL_sommerfeld(te, w, 1.0, "eq8")
    assert rel(ours, math.exp(-0.5 * xi) * lower) < 1e-12
    assert rel(ours, math.exp(0.5 * xi) * eq8) < 1e-12
    assert rel(ours, lower) < 0.02                      # the 2 % of the specification
    assert rel(ours, net) < xi * xi / 12                # twice the leading xi^2/24
    assert ours < net < lower
    # the module's quadrature against an independent integrator
    from scipy.integrate import quad
    hw = ib.HBAR * w / ib.QE

    def integrand(e):
        return float(ib.lnL_per_energy(e, w, n, n, screening="none")) * math.exp(-e / te) / te

    val, _ = quad(integrand, 0.5 * hw, 60.0 * te, limit=400, epsabs=0.0, epsrel=1e-12)
    assert rel(val, ours) < 1e-9


# --------------------------------------------------------------------------
# 4. the cross-section prefactor against nlib_kinetics
# --------------------------------------------------------------------------
@needs_nk
def test_sigma_ei_prefactor_matches_nlib_sigma_coulomb():
    eps = np.geomspace(10.0, 1000.0, 13)
    w, ne = TWO_PI * 10 * THZ, 1e21
    same_l = nk.L_ib(eps, w, ne)
    np.testing.assert_allclose(ib.sigma_ei_prefactor(eps) * same_l, nk.sigma_coulomb(eps, w, ne),
                               rtol=1e-12, atol=0.0)
    assert rel(ib.E4PE, nk.E4PE) < 1e-14
    # sigma_ei with the Sommerfeld L against nlib's classical/quantum b-ratio L: the
    # same prefactor, logarithms of the same size (a sanity bracket, not an identity)
    ratio = ib.sigma_ei(eps, w, ne, ne) / nk.sigma_coulomb(eps, w, ne)
    assert np.all((ratio > 0.5) & (ratio < 1.5)), ratio


def test_sigma_ei_prefactor_is_pi_b0_squared():
    eps = 2.0
    b0 = ib.QE**2 / (4 * math.pi * ib.EPS0 * (eps * ib.QE))        # e^2/(4 pi eps0 eps_J) [m]
    assert rel(float(ib.sigma_ei_prefactor(eps)), math.pi * b0 * b0) < 1e-14


# --------------------------------------------------------------------------
# 5. Drude limits of the Lorentz kernel
# --------------------------------------------------------------------------
def test_drude_limits_with_constant_collision_frequency(monkeypatch):
    nu0 = 3e11
    monkeypatch.setattr(ib, "nu_en_energy",
                        lambda eps, n_gas, x_gas: np.full_like(np.asarray(eps, dtype=float), nu0))
    te, e_rms = 1.3, 1e5
    kw = dict(n_gas_m3=1e25, x_gas={"Kr": 1.0}, ne_m3=0.0, ni_m3=0.0, ei_scale=0.0)
    p0 = ib.em_power_per_electron(te, 0.0, e_rms, **kw)
    drude = ib.QE**2 * e_rms**2 / (ib.ME * nu0) / ib.QE            # eV/s
    assert rel(p0, drude) < 1e-6
    w = 100.0 * nu0
    pw = ib.em_power_per_electron(te, w, e_rms, **kw)
    assert rel(pw / p0, nu0**2 / w**2) < 1e-2                       # exact: nu^2/(nu^2 + w^2)
    assert rel(pw / p0, nu0**2 / (nu0**2 + w**2)) < 1e-6
    parts = ib.em_power_per_electron(te, w, e_rms, parts=True, **kw)
    assert parts["ei"] == 0.0 and parts["f_L"] == 1.0
    assert rel(parts["nu_en_mean"], nu0) < 1e-6
    assert rel(parts["total"], parts["en"]) < 1e-14


# --------------------------------------------------------------------------
# 6. the electron-neutral part against LIGHT v0
# --------------------------------------------------------------------------
V0_POINT = dict(te=1.3, n_gas=2.45e25, e_rms=1.5e5, f_hz=95e9, x_gas={"Kr": 1.0})


@needs_v0
def test_en_part_matches_v0_em_power_per_electron():
    """v0's kernel with n_ion = 0 (cm^-3 and V/cm at its interface) on
    krypton, v0's default sigma_m (LIGHT checked its Kr/F2 mixture too)."""
    p = V0_POINT
    v0_kr = v0x.em_power_per_electron(p["te"], p["n_gas"] * 1e-6, 0.0, 1e15, p["f_hz"], p["e_rms"] / 100.0)
    ours_kr = ib.em_power_per_electron(p["te"], TWO_PI * p["f_hz"], p["e_rms"], p["n_gas"], p["x_gas"],
                                       1e21, 1e21, ei_scale=0.0)
    assert rel(ours_kr, v0_kr) < 0.02, (ours_kr, v0_kr)


def test_nu_en_energy_is_the_lxcat_momentum_transfer_sum():
    eps = np.array([0.5, 1.5, 4.0])
    n, x = 2.45e25, {"Kr": 1.0}
    s_kr = ib.xsec.find("Kr", "elastic").sigma_at(eps) * 1e-4
    v = np.sqrt(2 * ib.QE * eps / ib.ME)
    np.testing.assert_allclose(ib.nu_en_energy(eps, n, x), n * s_kr * v, rtol=1e-14)
    assert ib.sigma_mt("Kr").kind == "elastic" and list(ib.MT_BLOCK) == ["Kr"]


# --------------------------------------------------------------------------
# 7. the +-10 % model band
# --------------------------------------------------------------------------
BAND_POINT = dict(te=1.5, e_rms=1.5e5, n_gas=N_1ATM_300K, x_gas={"Kr": 1.0}, ne=1e22)


@pytest.mark.parametrize("f_thz", [0.095, 1.0, 10.0, 30.0])
def test_band_is_ordered(f_thz):
    p = BAND_POINT
    lo, mid, hi = ib.em_power_per_electron_band(p["te"], TWO_PI * f_thz * THZ, p["e_rms"], p["n_gas"],
                                                p["x_gas"], p["ne"], p["ne"])
    assert lo < mid < hi
    assert rel(mid, ib.em_power_per_electron(p["te"], TWO_PI * f_thz * THZ, p["e_rms"], p["n_gas"],
                                             p["x_gas"], p["ne"], p["ne"])) < 1e-14


@pytest.mark.parametrize("f_thz, tol", [(30.0, 1e-3), (10.0, 5e-3)])
def test_band_width_is_model_error_times_ei_share(f_thz, tol):
    """ei_scale scales nu_ei inside nu/(nu^2 + w^2), so (high - mid)/mid =
    TURNBULL_MODEL_ERROR x (e-i share) holds where the kernel is linear in
    nu, w >> nu: the deviation is O((nu/w)^2) (2e-4 at 30 THz, 3e-3 at
    10 THz; at 1 THz nu/w = 0.4 at 1 atm Kr and it is 10 % off)."""
    p = BAND_POINT
    w = TWO_PI * f_thz * THZ
    args = (p["te"], w, p["e_rms"], p["n_gas"], p["x_gas"], p["ne"], p["ne"])
    lo, mid, hi = ib.em_power_per_electron_band(*args)
    d = ib.em_power_per_electron(*args, parts=True)
    share = d["ei"] / d["total"]
    assert (d["nu_en_mean"] + d["nu_ei_mean"]) / w < 0.05
    assert rel((hi - mid) / mid, ib.TURNBULL_MODEL_ERROR * share) < tol
    assert rel((mid - lo) / mid, ib.TURNBULL_MODEL_ERROR * share) < tol
    assert ib.TURNBULL_MODEL_ERROR == 0.10


def test_parts_are_consistent():
    p = BAND_POINT
    w = TWO_PI * THZ
    d = ib.em_power_per_electron(p["te"], w, p["e_rms"], p["n_gas"], p["x_gas"], p["ne"], p["ne"], parts=True)
    assert set(d) == {"total", "en", "ei", "f_L", "nu_en_mean", "nu_ei_mean"}
    assert rel(d["total"], d["en"] + d["ei"]) < 1e-14
    # the Langdon factor multiplies the e-i part only
    d0 = ib.em_power_per_electron(p["te"], w, p["e_rms"], p["n_gas"], p["x_gas"], p["ne"], p["ne"],
                                  parts=True, langdon=False)
    assert d0["f_L"] == 1.0 and d0["en"] == d["en"]
    assert rel(d["ei"], d0["ei"] * ib.f_langdon(math.sqrt(2.0) * p["e_rms"], w, p["te"])) < 1e-12
    assert rel(d["f_L"], ib.f_langdon(math.sqrt(2.0) * p["e_rms"], w, p["te"])) < 1e-14
    # nu_ei_mean = (2/3T)<eps nu_ei> is Turnbull's rate with the Maxwellian average of L(eps)
    nu_turnbull_form = ib.NU_EI_PREF * p["te"]**-1.5 * p["ne"] * ib.lnL_maxwell(p["te"], w, p["ne"], p["ne"])
    assert rel(d["nu_ei_mean"], nu_turnbull_form) < 1e-3
    # ei_scale = 0 and ni = 0 both switch the ion terms off. The e-n part then moves a
    # little, because nu_ei leaves the denominator nu^2 + w^2 as well (0.2 % here).
    for kw in (dict(ei_scale=0.0), dict(ni_m3=0.0)):
        a = dict(te_eV=p["te"], w=w, E_rms_V_m=p["e_rms"], n_gas_m3=p["n_gas"], x_gas=p["x_gas"],
                 ne_m3=p["ne"], ni_m3=p["ne"], parts=True)
        a.update(kw)
        dd = ib.em_power_per_electron(**a)
        assert dd["ei"] == 0.0 and dd["nu_ei_mean"] == 0.0 and dd["f_L"] == 1.0
        assert dd["total"] == dd["en"]
        assert dd["en"] > d["en"] and rel(dd["en"], d["en"]) < 0.01


# --------------------------------------------------------------------------
# 8. the e-i share at the reference points
# --------------------------------------------------------------------------
def _share(ne, f_thz, te=1.5, e_rms=1.5e5):
    d = ib.em_power_per_electron(te, TWO_PI * f_thz * THZ, e_rms, N_1ATM_300K, {"Kr": 1.0}, ne, ne, parts=True)
    return d["ei"] / d["total"], d


def test_ei_share_sanity_bracket_1e21():
    """ne = ni = 1e21 m^-3, Te = 1.5 eV, 1 atm Kr (300 K), f = 1 THz.

    The specification's bracket was 1-60 %. The physics gives 0.35 %: at
    1 atm the Kr momentum-transfer frequency (2/3T)<eps nu_en> is 2.6e12/s
    against Turnbull's nu_ei of 6.4e9/s (lnL_M = 4.1, f_sc-like screening
    included), and v0's Spitzer kernel (NRL lnL = 6.3) gives 0.52 % at the
    same point. The share scales with n_e and crosses 1 % near 3e21 m^-3,
    so the floor here is 0.1 %; the 1-60 % bracket is tested at 1e22."""
    share, d = _share(1e21, 1.0)
    assert 1e-3 < share < 0.6, f"e-i share {share:.4%} (f_L {d['f_L']:.4f}, nu_en {d['nu_en_mean']:.3e}, nu_ei {d['nu_ei_mean']:.3e})"


def test_ei_share_sanity_bracket_1e22():
    share, d = _share(1e22, 1.0)
    assert 0.01 < share < 0.6, f"e-i share {share:.4%} (f_L {d['f_L']:.4f}, nu_en {d['nu_en_mean']:.3e}, nu_ei {d['nu_ei_mean']:.3e})"
    share10, _ = _share(1e22, 10.0)
    assert 0.01 < share10 < 0.6, f"e-i share at 10 THz {share10:.4%}"
    # the share grows with n_i (nu_ei ~ n_i) and shrinks with w at w >> nu through the
    # smaller logarithm (b_max = v/w)
    share21, _ = _share(1e21, 1.0)
    assert share21 < share and share10 < share


# --------------------------------------------------------------------------
# propagation hints
# --------------------------------------------------------------------------
def test_plasma_frequency_and_critical_density():
    w = TWO_PI * 95e9
    assert rel(ib.plasma_frequency(ib.n_crit(w)), w) < 1e-14
    # v0's report: n_crit = 1.24e-8 f^2 cm^-3
    assert rel(ib.n_crit(w) * 1e-6, 1.24e-8 * 95e9**2) < 5e-3
    h = ib.skin_hints(1.5, w, 1e21, 3e11)
    assert rel(h["ne_over_n_crit"], 1e21 / ib.n_crit(w)) < 1e-14
    assert rel(h["nu_over_w"], 3e11 / w) < 1e-14
    assert rel(h["w_p_over_w"], ib.plasma_frequency(1e21) / w) < 1e-14
    assert "0-D" in h["note"]


def test_module_uses_light_xsec_reader():
    assert hasattr(ib.xsec, "find") and hasattr(ib.xsec, "processes")
    assert ib.xsec.__name__.endswith(".xsec")
