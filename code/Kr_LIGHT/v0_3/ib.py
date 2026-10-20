"""Inverse-bremsstrahlung (IB) heating of the electrons of a weakly ionized
krypton gas by an EM wave of angular frequency w (LIGHT v0_4's module,
written for its Kr/F2 laser gas; LIGHT_Kr carries krypton alone).

The electron-neutral part is the Lorentz-gas kernel of LIGHT v0
(v0/xsec.py em_power_per_electron); the electron-ion part follows

  [T23] D. Turnbull et al., Phys. Rev. Lett. 130, 145103 (2023)
  [T24] D. Turnbull et al., Phys. Plasmas 31, 063304 (2024)

with the pieces of their composite model:

  * the Sommerfeld Coulomb logarithm lnL_IB = (pi/sqrt 3) <g_ff>
    [T24 Eq. (5)] from Pradler & Semmelrock's fit to the velocity-
    dependent Gaunt factor [T24 Eq. (7)], thermally averaged by
    [T24 Eq. (8)];
  * the ion-screening factor f_sc [T24 Eq. (11)] at the ion-sphere radius
    l_sc = a_i = (3/(4 pi n_i))^(1/3);
  * the Langdon factor f_L: the original [T24 Eq. (9)] or the new
    [T24 Eq. (10)] over Matte's super-Gaussian [T24 Eqs. (1)-(2)] at the
    Langdon parameter alpha_L = Z_eff v_o^2/v_t^2 [T24 Eq. (3)].

Two forms are provided. (i) Turnbull's own Maxwellian-averaged rate
nu_ei [T24 Eq. (4)-(5)] with f_sc and f_L as multiplicative corrections
(nu_ei_maxwell, f_screen, f_langdon). (ii) A per-electron-energy
logarithm L(eps) (lnL_per_energy) and cross section sigma_ei(eps) for a
kinetic (Boltzmann) solver: the Pradler Gaunt factor at the upper energy
of the absorption jump plus the Born-Yukawa screening term per electron,
zero below the threshold eps < hbar w/2. For a Maxwellian, (2/3T)<eps
nu_ei(eps)> with nu_ei(eps) = n_i sigma_ei v equals form (i) with
lnL = int L(u T) e^-u du, so the two forms are the same physics
(see lnL_maxwell for the exact relation to [T24 Eq. (8)]).

Units: SI inside. At the interface energies and temperatures are in eV,
densities in m^-3, fields in V/m, w in rad/s, powers per electron in eV/s.
Cross sections read from the LXCat files through xsec are cm^2 and are
converted here.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gamma as Gamma

try:                                   # inside the v0_3 package
    from . import xsec
except ImportError:                    # script use: LIGHT_Kr/ on sys.path
    from v0_3 import xsec

# --------------------------------------------------------------------------
# constants (CODATA 2018, as turnbull.py and xsec.py)
# --------------------------------------------------------------------------
ME = 9.1093837015e-31          # kg
QE = 1.602176634e-19           # C (also J per eV)
EPS0 = 8.8541878128e-12        # F/m
HBAR = 1.054571817e-34         # J s
C = 299792458.0                # m/s
KB = 1.380649e-23              # J/K
ALPHA = QE**2 / (4 * math.pi * EPS0 * HBAR * C)       # fine-structure constant
EULER = 0.5772156649015329                             # Euler's gamma
RY = 0.5 * ME * (ALPHA * C) ** 2 / QE                  # Rydberg energy [eV]
SQ3PI = math.sqrt(3.0) / math.pi
E4PE = QE / (4 * math.pi * EPS0)                       # e/(4 pi eps0) [V m]
# nu_ei = NU_EI_PREF T^-3/2 sum Z^2 n_i lnL  (the 2.91e-12 of [T24], SI, T in eV)
NU_EI_PREF = 4 * math.sqrt(2 * math.pi) / 3 * (QE * E4PE) ** 2 / (math.sqrt(ME) * QE**1.5)

# Turnbull et al.'s composite model (Sommerfeld lnL with the Pradler Gaunt
# factor, ion screening at a_i, the Langdon factor) reproduces their
# measured absorption to about 10 %, so the electron-ion absorption carries
# a +-10 % model band. Their words (read from the local PDFs, 2026-10-08):
# [T24] Sec. IV, p. 063304-8: "It is particularly satisfying to see all of
# the mixture shots with measured absorption > 1.5% (which are therefore
# the least affected by the absolute measurement uncertainty) all predicted
# to better than 10% and well within the uncertainty."; Sec. I: "a
# screening correction factor adapted from Rozsnyai was added because the
# models still overestimated absorption by about 10% unless it was
# included." [T23] p. 145103-4: with screening "Oster being on average
# just 2% below the data, and Sommerfeld just 2% above" over 22 shots.
# The 10 % is the per-shot scatter of the best data, not the 2 % average.
TURNBULL_MODEL_ERROR = 0.10

# Maxwellian quadrature for em_power_per_electron: geometric grid in eps,
# N_EPS points from EPS_LO_OVER_TE*Te to EPS_HI_OVER_TE*Te (trapezoid
# weights on a geometric grid are the uniform-in-ln(eps) trapezoid rule up
# to a constant factor that the normalization removes, so the average is
# spectrally accurate; the lower end leaves < 1e-9 of the Maxwellian out).
N_EPS = 4000
EPS_LO_OVER_TE = 1e-6
EPS_HI_OVER_TE = 50.0


# --------------------------------------------------------------------------
# Gaunt factor and its thermal average   [T24 Eqs. (5), (7), (8)]
# --------------------------------------------------------------------------
def gff_pradler(nu, x):
    """[T24 Eq. (7)]: Pradler & Semmelrock's fit to Sommerfeld's free-free
    Gaunt factor. nu = alpha Z c/v_i is the Sommerfeld parameter of the
    initial (upper) velocity, x = hbar w/E_i in (0, 1]. The last term
    stands outside the bracket, as printed. (1+sqrt(1-x))/(1-sqrt(1-x)) is
    written (1+sqrt(1-x))^2/x and (pi nu)^2/sinh^2(pi nu) with e^(-pi nu),
    which keep their digits at small x and large nu."""
    nu = np.asarray(nu, dtype=float)
    x = np.asarray(x, dtype=float)
    s = np.sqrt(np.clip(1.0 - x, 0.0, None))
    ratio = (1.0 + s) ** 2 / x
    eg = math.exp(EULER)
    t1 = np.log1p(ratio / (nu * eg))
    q = np.exp(-math.pi * nu)
    sinh_term = (2 * math.pi * nu * q / (-np.expm1(-2 * math.pi * nu))) ** 2
    t2 = sinh_term * np.log(nu * eg)
    t3 = (10.0 / 3.0) * nu**2 * np.exp(-1.5 * math.pi * nu) * (
        5.0 / (1.0 + np.exp(10.0 * (x - 6.0 / 7.0))) - 2.0)
    return SQ3PI * (t1 + t2 - t3) + x * nu / (math.pi ** (-2.0 / 3.0) + x * nu)


def _s_nodes(xi, s_max=8.0, panels=72, order=10):
    """Gauss-Legendre nodes and weights for int_0^s_max ds, panels geometric
    in s near 0 (scale sqrt xi) and uniform beyond 1. Used with
    u = xi + s^2 (du = 2 s ds), which removes the sqrt(1-x) of the
    threshold from the integrand."""
    g_x, g_w = np.polynomial.legendre.leggauss(order)
    s0 = max(1e-3 * math.sqrt(xi), 1e-12)
    n_geo = panels // 2
    edges = np.concatenate([[0.0], np.geomspace(s0, 1.0, n_geo),
                            np.linspace(1.0, s_max, panels - n_geo + 1)[1:]])
    a, b = edges[:-1], edges[1:]
    s = (0.5 * (b - a)[:, None] * (g_x[None, :] + 1.0) + a[:, None]).ravel()
    w = (0.5 * (b - a)[:, None] * g_w[None, :]).ravel()
    return s, w


def thermal_gff(te_eV, w, Z=1.0, average="net"):
    """Thermal average <g_ff> at photon angular frequency w [rad/s] and
    Te [eV]; xi = hbar w/(e Te), nu_p = Z sqrt(Ry/Te).

    average="eq8": [T24 Eq. (8)] as printed, int_xi^inf g(xi/u, nu_p/sqrt u)
    e^-u du over the initial (upper) energy u;
    "lower": the usual average over the lower energy, e^xi times eq8;
    "net": absorption minus stimulated emission of a Maxwellian,
    (1 - e^-xi)/xi times "lower" = (e^xi - 1)/xi times eq8.
    The three agree to O(xi); xi < 0.1 for Te > 0.5 eV below 30 THz."""
    xi = HBAR * w / (QE * te_eV)
    nu_p = Z * math.sqrt(RY / te_eV)
    s, ws = _s_nodes(xi)
    u = xi + s**2
    val = float(np.sum(ws * 2 * s * gff_pradler(nu_p / np.sqrt(u), xi / u) * np.exp(-u)))
    if average == "eq8":
        return val
    if average == "lower":
        return val * math.exp(xi)
    if average == "net":
        return val * math.expm1(xi) / xi
    raise ValueError(f"average {average!r}: use 'eq8', 'lower' or 'net'")


def lnL_thermal(te_eV, w, Z=1.0, average="net"):
    """[T24 Eq. (5)]: lnL_IB = (pi/sqrt 3) <g_ff>, Maxwellian-averaged."""
    return thermal_gff(te_eV, w, Z, average) / SQ3PI


# --------------------------------------------------------------------------
# per-electron-energy logarithm and cross section (kinetic form)
# --------------------------------------------------------------------------
def a_ws(ni_m3):
    """Ion-sphere (Wigner-Seitz) radius a_i = (3/(4 pi n_i))^(1/3) [m], the
    screening length l_sc of [T24 Eq. (11)]."""
    return (3.0 / (4 * math.pi * ni_m3)) ** (1.0 / 3.0)


def debye_length(ne_m3, te_eV):
    """Electron Debye length sqrt(eps0 Te/(n_e e)) [m] (Te in eV)."""
    return math.sqrt(EPS0 * te_eV / (ne_m3 * QE))


def _upper_energy(eps_eV, w):
    """E_i [eV] of the absorption jump eps - hbar w/2 -> eps + hbar w/2
    across the face eps (clipped to >= hbar w so the fit's arguments stay
    valid where the threshold cut applies), and hbar w [eV]."""
    hw = HBAR * w / QE
    eps = np.asarray(eps_eV, dtype=float)
    return np.maximum(eps + 0.5 * hw, hw), hw, eps


def L_sommerfeld_energy(eps_eV, w, Z=1.0):
    """(pi/sqrt 3) g_ff(nu_i, x) of [T24 Eq. (7)] at the upper energy
    E_i = eps + hbar w/2 of the jump across the face eps (nu_i =
    Z sqrt(Ry/E_i), x = hbar w/E_i); zero for eps < hbar w/2, where the
    lower energy would be negative. Its Maxwellian average is e^(xi/2)
    times [T24 Eq. (8)] (see lnL_maxwell)."""
    Ei, hw, eps = _upper_energy(eps_eV, w)
    nu = Z * math.sqrt(RY) / np.sqrt(Ei)
    return np.where(eps >= 0.5 * hw, gff_pradler(nu, hw / Ei) / SQ3PI, 0.0)


def screen_add_energy(eps_eV, w, l_sc):
    """Born-Yukawa change of the logarithm per electron at the upper
    energy: ln(y/sqrt(1+y^2)) - 1/(2(1+y^2)), y = l_sc w/v_i. This is the
    numerator of [T24 Eq. (11)] minus its denominator, so that averaged
    over the Maxwellian with lnL it reproduces f_sc lnL at v_i = v_t."""
    Ei, hw, eps = _upper_energy(eps_eV, w)
    v = np.sqrt(2 * QE * Ei / ME)
    y = l_sc * w / v
    return np.log(y / np.sqrt(1 + y * y)) - 0.5 / (1 + y * y)


def lnL_per_energy(eps_eV, w, ne_m3, ni_m3, screening="turnbull", te_eV=None, floor=True):
    """Per-electron-energy Coulomb logarithm L(eps) [dimensionless] for e-i
    free-free absorption at angular frequency w [rad/s], for a kinetic
    solver: sigma_ei(eps) = pi b0^2 L(eps).

    L = L_Sommerfeld(eps) + screening term, where
      * L_Sommerfeld is (pi/sqrt 3) g_ff of [T24 Eq. (7)] at the upper
        energy of the jump ([T24 Eq. (5)] before the thermal average of
        [T24 Eq. (8)]), zero below the threshold eps < hbar w/2;
      * screening="turnbull": the Born-Yukawa term of [T24 Eq. (11)] per
        electron at l_sc = a_i = (3/(4 pi n_i))^(1/3) (ion sphere);
        "debye": the same at the electron Debye length at te_eV (required);
        "none": no screening term.
    floor=True clips L at zero (a negative cross section is meaningless;
    the screened logarithm goes negative for slow electrons at low w).
    Returns an array of the shape of eps_eV. w must be > 0."""
    if w <= 0:
        raise ValueError("lnL_per_energy needs w > 0 (the IB logarithm diverges at w = 0)")
    hw = HBAR * w / QE
    eps = np.asarray(eps_eV, dtype=float)
    L = L_sommerfeld_energy(eps, w)
    if screening == "turnbull":
        L = L + screen_add_energy(eps, w, a_ws(ni_m3))
    elif screening == "debye":
        if te_eV is None:
            raise ValueError("screening='debye' needs te_eV")
        L = L + screen_add_energy(eps, w, debye_length(ne_m3, te_eV))
    elif screening != "none":
        raise ValueError(f"screening {screening!r}: use 'turnbull', 'debye' or 'none'")
    L = np.where(eps >= 0.5 * hw, L, 0.0)
    if floor:
        L = np.maximum(L, 0.0)
    return L


def sigma_ei_prefactor(eps_eV):
    """pi b0^2 [m^2], b0 = e^2/(4 pi eps0 eps_J) = (e/(4 pi eps0))/eps_eV the
    classical distance of closest approach at energy eps."""
    return math.pi * (E4PE / np.asarray(eps_eV, dtype=float)) ** 2


def sigma_ei(eps_eV, w, ne_m3, ni_m3, **kw):
    """Momentum-transfer-like e-i cross section per singly charged ion for
    the kinetic equation [m^2]: sigma = pi b0^2 L(eps), b0 = e^2/(4 pi eps0
    eps_J), L from lnL_per_energy (kw: screening, te_eV, floor). For a
    Maxwellian the absorbed power of nu_ei(eps) = n_i sigma v is that of
    Turnbull's nu_ei [T24 Eq. (5)] with lnL = int L(u Te) e^-u du."""
    return sigma_ei_prefactor(eps_eV) * lnL_per_energy(eps_eV, w, ne_m3, ni_m3, **kw)


def lnL_maxwell(te_eV, w, ne_m3, ni_m3, **kw):
    """Maxwellian average of the per-energy logarithm, lnL_M = int_0^inf
    L(eps) e^(-eps/Te) d eps/Te, on the quadrature of thermal_gff (u =
    xi/2 + s^2). This is the lnL that makes NU_EI_PREF Te^-3/2 n_i lnL_M
    equal to (2/3Te)<eps n_i sigma_ei v>_Maxwell exactly.

    Relation to [T24 Eq. (8)]: with screening="none", L(eps) is the Gaunt
    factor at the upper energy E_i = eps + hbar w/2, and the Maxwellian
    weight is e^(-eps/Te) = e^(xi/2) e^(-E_i/Te), so lnL_M = e^(xi/2) lnL_eq8
    = e^(-xi/2) lnL_lower exactly; it agrees with the "net" average
    (e^xi - 1)/xi lnL_eq8 to O(xi^2) (xi^2/24)."""
    xi = HBAR * w / (QE * te_eV)
    s, ws = _s_nodes(xi)
    u = 0.5 * xi + s**2
    L = lnL_per_energy(u * te_eV, w, ne_m3, ni_m3, **kw)
    return float(np.sum(ws * 2 * s * L * np.exp(-u)))


# --------------------------------------------------------------------------
# Turnbull's Maxwellian-averaged rate and its corrections
# --------------------------------------------------------------------------
def nu_ei_maxwell(te_eV, w, ne_m3, ni_m3, Z=1.0, average="net", **kw):
    """Turnbull's electron-ion collision frequency [1/s], [T24 Eq. (5) in
    Eq. (4)]: nu_ei = NU_EI_PREF Te^-3/2 Z^2 n_i lnL_IB with lnL_IB =
    (pi/sqrt 3)<g_ff> from [T24 Eqs. (7)-(8)], the thermal average as in
    thermal_gff (default "net": absorption minus stimulated emission;
    "eq8" is [T24 Eq. (8)] as printed). ne_m3 and kw are accepted for
    interface symmetry with the per-energy functions and not used: n_e
    enters Turnbull's rate only through n_e/n_c in kappa [T24 Eq. (4)]."""
    return NU_EI_PREF * te_eV**-1.5 * Z * Z * ni_m3 * lnL_thermal(te_eV, w, Z, average)


def f_screen(te_eV, w, ni_m3):
    """Ion-screening factor [T24 Eq. (11)] at l_sc = a_i(n_i):
    [ln((2 Te e/hbar w) y/sqrt(1+y^2)) - 1/(2(1+y^2))] / ln(2 Te e/hbar w),
    y = l_sc w/v_t, v_t = sqrt(e Te/m_e). Multiplies nu_ei_maxwell."""
    vt = math.sqrt(QE * te_eV / ME)
    y = a_ws(ni_m3) * w / vt
    L0 = math.log(2 * QE * te_eV / (HBAR * w))
    return (math.log(2 * QE * te_eV / (HBAR * w) * y / math.sqrt(1 + y * y)) - 0.5 / (1 + y * y)) / L0


def alpha_langdon(E0_V_m, w, te_eV, Z_eff=1.0):
    """Langdon parameter [T24 Eq. (3)]: alpha_L = Z_eff v_o^2/v_t^2 with the
    quiver velocity v_o = e E0/(m_e w) of the field AMPLITUDE E0 =
    sqrt(2) E_rms and v_t = sqrt(e Te/m_e)."""
    vo = QE * E0_V_m / (ME * w)
    return Z_eff * vo**2 / (QE * te_eV / ME)


def matte_m(alpha_L):
    """Matte's super-Gaussian exponent [T24 Eq. (2)]:
    m = 2 + 3/(1 + 1.66/alpha_L^0.724)."""
    return 2.0 + 3.0 / (1.0 + 1.66 / alpha_L**0.724) if alpha_L > 0 else 2.0


def fL_orig(alpha_L):
    """Langdon's original factor [T24 Eq. (9)]: 1 - 0.553/[1 + (0.27/alpha_L)^0.75]."""
    return 1.0 - 0.553 / (1.0 + (0.27 / alpha_L) ** 0.75) if alpha_L > 0 else 1.0


def _a_m(m):
    """m_e v_m^2/(2 e Te) = (3/2) Gamma(3/m)/Gamma(5/m) of [T24 Eq. (1)] (1 at m = 2)."""
    return 1.5 * Gamma(3.0 / m) / Gamma(5.0 / m)


def _cm_over_c2(m):
    """f_SG(0)/f_Maxwell(0) at equal density and temperature, c_m/c_2 of [T24 Eq. (1)]."""
    return m * Gamma(1.5) / (2.0 * Gamma(3.0 / m)) * _a_m(m) ** -1.5


def fL_new(m, te_eV, w, Z=1.0):
    """The new Langdon factor [T24 Eq. (10)] at exponent m: the absorption
    of the super-Gaussian [T24 Eq. (1)] over that of the Maxwellian with the
    velocity-dependent Gaunt factor,
      (m/2) a_m^(-m/2) (c_m/c_2) int g e^-(u/a_m)^(m/2) u^(m/2-1) du / int g e^-u du,
    u = E_i/Te, g zero below the threshold u = xi as in [T24 Eq. (8)]."""
    xi = HBAR * w / (QE * te_eV)
    nu_p = Z * math.sqrt(RY / te_eV)
    s, ws = _s_nodes(xi)
    u = xi + s**2
    g = gff_pradler(nu_p / np.sqrt(u), xi / u)
    am = _a_m(m)
    num = float(np.sum(ws * 2 * s * g * np.exp(-((u / am) ** (m / 2))) * u ** (m / 2 - 1)))
    den = float(np.sum(ws * 2 * s * g * np.exp(-u)))
    return float(0.5 * m * am ** (-m / 2) * _cm_over_c2(m) * num / den)


def f_langdon(E0_V_m, w, te_eV, Z_eff=1.0, version="new"):
    """Langdon factor f_L on the Maxwellian-averaged e-i absorption:
    alpha_L [T24 Eq. (3)] -> m [T24 Eq. (2)] -> f_L, version "new" =
    [T24 Eq. (10)] (fL_new(m, te, w)), "orig" = [T24 Eq. (9)]. E0_V_m is
    the field AMPLITUDE, E0 = sqrt(2) E_rms. Z_eff enters alpha_L only;
    the Gaunt factor inside fL_new is taken at Z = 1 (singly charged
    Kr+, Kr2+, F-). A solver that computes the EEDF must not apply f_L:
    the distortion it stands for is then already in the distribution."""
    aL = alpha_langdon(E0_V_m, w, te_eV, Z_eff)
    if version == "new":
        return fL_new(matte_m(aL), te_eV, w)
    if version == "orig":
        return fL_orig(aL)
    raise ValueError(f"version {version!r}: use 'new' or 'orig'")


# --------------------------------------------------------------------------
# electron-neutral momentum transfer from the LXCat files
# --------------------------------------------------------------------------
# gas -> (LXCat block kind, product) that holds its momentum-transfer total:
# Kr: the Biagi ELASTIC block (momentum transfer), what a Lorentz-gas
# kernel wants. (LIGHT's Kr/F2 also read F2's SIGLO EFFECTIVE block.)
MT_BLOCK = {"Kr": ("elastic", None)}
_MT_CACHE = {}


def sigma_mt(gas):
    """The momentum-transfer Process of a gas (cached)."""
    if gas not in _MT_CACHE:
        kind, product = MT_BLOCK[gas]
        _MT_CACHE[gas] = xsec.find(gas, kind, product)
    return _MT_CACHE[gas]


def speed(eps_eV):
    """v = sqrt(2 e eps/m_e) [m/s]."""
    return np.sqrt(2.0 * QE * np.asarray(eps_eV, dtype=float) / ME)


def nu_en_energy(eps_eV, n_gas_m3, x_gas):
    """Electron-neutral momentum-transfer collision frequency [1/s] at each
    energy: n_gas sum_g x_g sigma_m,g(eps) v(eps), sigma_m,g from the
    LXCat blocks of MT_BLOCK (cm^2 in xsec, converted to m^2)."""
    eps = np.asarray(eps_eV, dtype=float)
    sig = np.zeros_like(eps)
    for gas, x in x_gas.items():
        if x:
            sig = sig + x * sigma_mt(gas).sigma_at(eps) * 1e-4
    return n_gas_m3 * sig * speed(eps)


# --------------------------------------------------------------------------
# absorbed power of a Maxwellian electron gas (Lorentz-gas kernel)
# --------------------------------------------------------------------------
def maxwell_grid(te_eV):
    """(eps [eV], normalized Maxwellian weights) for <.>_Maxwell: geometric
    grid, trapezoid weights, g = (2/sqrt pi) Te^-3/2 sqrt(eps) e^(-eps/Te),
    the weights summing to 1."""
    eps = np.geomspace(EPS_LO_OVER_TE * te_eV, EPS_HI_OVER_TE * te_eV, N_EPS)
    wq = np.zeros_like(eps)
    d = np.diff(eps)
    wq[:-1] += 0.5 * d
    wq[1:] += 0.5 * d
    g = np.sqrt(eps) * np.exp(-eps / te_eV) * wq
    return eps, g / g.sum()


def em_power_per_electron(te_eV, w, E_rms_V_m, n_gas_m3, x_gas, ne_m3, ni_m3,
                          parts=False, langdon=True, screening="turnbull", ei_scale=1.0):
    """Absorbed power per electron [eV/s] of a MAXWELLIAN electron gas at
    Te in the field E_rms cos(w t), Lorentz-gas kernel with the
    energy-dependent collision frequency nu(eps) = nu_en(eps) + ei_scale
    nu_ei(eps), nu_ei(eps) = n_i sigma_ei(eps) v(eps) (the Turnbull
    per-energy logarithm, lnL_per_energy with `screening`):

        P = (e^2 E_rms^2/m) < (2 eps/(3 Te)) nu(eps)/(nu(eps)^2 + w^2) >_Maxwell

    (v0's em_power_per_electron with its Spitzer nu_ei replaced). For
    constant nu this is Drude: e^2 E_rms^2/(m nu) at w = 0 and e^2 E_rms^2
    nu/(m w^2) for w >> nu. E_rms is the field INSIDE the plasma (no
    propagation, reflection or skin depth here: see skin_hints).

    The e-i part is the nu_ei share of the numerator, P_ei = (e^2 E^2/m)
    <(2 eps/3Te) ei_scale nu_ei/(nu^2 + w^2)>; langdon=True multiplies it
    by f_langdon(sqrt 2 E_rms, w, Te) ([T24 Eq. (10)]): Turnbull's Langdon
    factor corrects the Maxwellian-averaged e-i absorption for the
    super-Gaussian distortion of the field-heated electrons. A kinetic
    solver that solves the EEDF must NOT apply it (double counting).
    ei_scale scales nu_ei inside the kernel (the +-TURNBULL_MODEL_ERROR
    band of em_power_per_electron_band); with ei_scale = 0 the ion
    terms are not evaluated (w = 0 is then allowed).

    parts=True returns {"total", "en", "ei" [eV/s], "f_L", "nu_en_mean",
    "nu_ei_mean" [1/s]} with nu_X_mean = (2/3Te)<eps nu_X(eps)>, the
    effective collision frequencies of the w >> nu limit (P_X -> (e^2
    E^2/m) nu_X_mean/w^2); nu_ei_mean includes ei_scale and equals
    NU_EI_PREF Te^-3/2 n_i lnL_maxwell, not f_L."""
    eps, g = maxwell_grid(te_eV)
    nu_en = np.asarray(nu_en_energy(eps, n_gas_m3, x_gas), dtype=float) * np.ones_like(eps)
    if ei_scale != 0.0 and ni_m3 > 0.0:
        nu_ei = ei_scale * ni_m3 * sigma_ei(eps, w, ne_m3, ni_m3, screening=screening,
                                            te_eV=te_eV) * speed(eps)
    else:
        nu_ei = np.zeros_like(eps)
    nu = nu_en + nu_ei
    kern = (2.0 * eps / (3.0 * te_eV)) / (nu * nu + w * w)
    pref = QE**2 * E_rms_V_m**2 / ME / QE          # eV/s per (1/s)
    p_en = pref * float(np.sum(g * kern * nu_en))
    p_ei = pref * float(np.sum(g * kern * nu_ei))
    f_L = f_langdon(math.sqrt(2.0) * E_rms_V_m, w, te_eV) if (langdon and p_ei != 0.0) else 1.0
    p_ei *= f_L
    total = p_en + p_ei
    if not parts:
        return total
    return {"total": total, "en": p_en, "ei": p_ei, "f_L": f_L,
            "nu_en_mean": float(np.sum(g * (2.0 * eps / (3.0 * te_eV)) * nu_en)),
            "nu_ei_mean": float(np.sum(g * (2.0 * eps / (3.0 * te_eV)) * nu_ei))}


def em_power_per_electron_band(te_eV, w, E_rms_V_m, n_gas_m3, x_gas, ne_m3, ni_m3,
                               langdon=True, screening="turnbull"):
    """(low, mid, high) [eV/s]: em_power_per_electron at ei_scale =
    1 - TURNBULL_MODEL_ERROR, 1, 1 + TURNBULL_MODEL_ERROR, the +-10 %
    model band of the electron-ion absorption."""
    return tuple(em_power_per_electron(te_eV, w, E_rms_V_m, n_gas_m3, x_gas, ne_m3, ni_m3,
                                       langdon=langdon, screening=screening, ei_scale=s)
                 for s in (1.0 - TURNBULL_MODEL_ERROR, 1.0, 1.0 + TURNBULL_MODEL_ERROR))


# --------------------------------------------------------------------------
# propagation hints (0-D: printed, not applied)
# --------------------------------------------------------------------------
def plasma_frequency(ne_m3):
    """w_p = sqrt(n_e e^2/(eps0 m_e)) [rad/s]."""
    return math.sqrt(ne_m3 * QE**2 / (EPS0 * ME))


def n_crit(w):
    """Critical density eps0 m_e w^2/e^2 [m^-3] (w_p = w)."""
    return EPS0 * ME * w**2 / QE**2


def skin_hints(te_eV, w, ne_m3, nu_mean):
    """What v0's report prints beside an EM run: n_e/n_crit, nu/w, w_p/w
    (and the inputs), for judging whether the wave at w reaches the
    electrons as a uniform field. te_eV is carried for the caller's
    bookkeeping (the hints do not depend on it)."""
    nc = n_crit(w)
    wp = plasma_frequency(ne_m3)
    return {"te_eV": te_eV, "w": w, "ne_m3": ne_m3, "n_crit_m3": nc,
            "ne_over_n_crit": ne_m3 / nc, "w_p": wp, "w_p_over_w": wp / w,
            "nu_mean": nu_mean, "nu_over_w": nu_mean / w,
            "note": "E_rms is the field inside the plasma; no propagation, "
                    "reflection or skin depth in 0-D"}
