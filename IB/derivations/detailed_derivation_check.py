#!/usr/bin/env python3
"""Numerical verification of detailed_derivation.tex (projects/IB, quantum inverse bremsstrahlung).

Every boxed result of the notes is checked here, in the order of the notes, and
every check carries the LaTeX label of the equation (or table, or appendix
section) it verifies, so that Appendix C of the notes can be read against this
output line by line.  The script also writes Appendix C's table,
detailed_derivation_checks.tex, which the notes \\input.

  2  classical baseline   Drude conductivity (and the printed sign of Atzeni's
                          Eq. 11.8), the absorbed power by time integration of
                          Newton's equation, the absorption coefficient from the
                          complex wave number
  3  Hamiltonian          quiver motion by ODE integration, the Yukawa Fourier
                          transform by radial quadrature, the KH cancellation
  4  Volkov states        the phase by finite differences, Jacobi-Anger, the
                          three sum rules, the generating function of J^2
  5  rates                the Volkov matrix element against its Bessel series,
                          the golden-rule limit at finite times, the ion average
                          by Monte Carlo, the delta-function integral, Rutherford,
                          the one-photon cross section at fixed orientation by 2-D
                          quadrature, Bunkin-Fedorov Eq. (10) as printed, the
                          arcsine limit
  6  kinetic equation     the operator equation of motion in a fermionic Fock
                          space (Jordan-Wigner), the Markov limit, no blocking,
                          conservation, the two forms of P, positivity
  7  Lindhard             Im chi0 by direct 3-D quadrature (Maxwell, Fermi-Dirac,
                          an arbitrary f), the classical limit
  8  closed forms         Lambda(u) by the cross-section route, limits, Table 2,
                          Fermi-Dirac with explicit Pauli factors, T = 0, Ziman,
                          the universal tail for three distributions
  9  Gaunt / Kirchhoff    the mode sum, the Born cross section, the Gaunt factor,
                          its thermal and power averages, Kirchhoff's law
 10  many-body / RPA      the two-frame identity, chi0 by principal-value
                          quadrature against the Z form, the dielectric identity,
                          the classical limit = Atzeni Eq. (11.19), Table 4
 11  strong fields        the B_l recursion, the exact sum against the 4-D
                          cross-section route for l = 1, 2, Table 5, the
                          quasi-static limit, g, G = 3Q/r^3 (Silin), the limits
                          and the constant c_S
 12  Langdon              detailed balance of the rates, the Fokker-Planck limit
                          against the exact averaged master equation, F(0) law,
                          the super-Gaussian factors
  A  Atzeni               every recomputed number of Table A.1 and Table 1
  B  identities           the erf integral, ln 2, the Gaussian phase average,
                          the large-field average

Units: Gaussian.  Dimensionless checks use m = hbar = k_B T = 1 unless stated.
Run:  /usr/bin/python3 detailed_derivation_check.py      (numpy + scipy)
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import quad, solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.special import jv, k0, kve, erf, wofz, ive, gamma as Gamma

T0 = time.time()
EG = 0.5772156649015329
FAILS = []
RESULTS = {}          # label -> list of (name, err, ok)
ORDER = []
HERE = os.path.dirname(os.path.abspath(__file__))


def _log(label, name, ok, err, kind='quant'):
    """kind: 'quant' (a relative or absolute deviation), 'printed' (a number printed in the notes,
    reproduced to its printed digits), 'logical' (a condition)."""
    if label not in RESULTS:
        RESULTS[label] = []
        ORDER.append(label)
    RESULTS[label].append((name, err, ok, kind))
    if not ok:
        FAILS.append(f"{label}: {name}")


def report(label, name, got, want, tol):
    """relative error |got - want| / max(|want|, tiny) below tol."""
    got = complex(got) if np.iscomplexobj(got) else float(got)
    want = complex(want) if np.iscomplexobj(want) else float(want)
    err = float(abs(got - want) / max(abs(want), 1e-300))
    ok = err < tol
    print(f"  {'PASS' if ok else 'FAIL':4s} [{label}] {name:66s} err={err:.1e}")
    _log(label, name, ok, err)
    return ok


def absreport(label, name, got, want, tol):
    """absolute error below tol."""
    err = float(abs(got - want))
    ok = err < tol
    print(f"  {'PASS' if ok else 'FAIL':4s} [{label}] {name:66s} abs={err:.1e}")
    _log(label, name, ok, err)
    return ok


def printed(label, name, value, shown, decimals):
    """a number printed in the notes: |value - shown| <= half a unit of the last digit."""
    err = abs(value - shown)
    ok = err <= 0.5 * 10.0**(-decimals) * (1 + 1e-9)
    print(f"  {'PASS' if ok else 'FAIL':4s} [{label}] {name:66s} {value:.6g} -> {shown}")
    _log(label, name, ok, None, kind='printed')
    return ok


def check(label, name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL':4s} [{label}] {name:66s} {detail}")
    _log(label, name, bool(cond), None, kind='logical')
    return cond


def rich(G, h):
    """central difference at h and h/2, Richardson-extrapolated (error O(h^4))."""
    d1 = (G(h) - G(-h)) / (2 * h)
    d2 = (G(h / 2) - G(-h / 2)) / h
    return (4 * d2 - d1) / 3


def Lam(u):
    return -np.expm1(-u) / u * kve(0, u / 2)          # kve(0, x) = e^x K0(x), safe at large u


def gavg(u):
    return np.sqrt(3) / np.pi * kve(0, u / 2)


def Zf(z):
    return 1j * np.sqrt(np.pi) * wofz(z)


# CGS constants (CODATA 2018)
C_LIGHT = 2.99792458e10
E_CH = 4.80320471e-10
M_E = 9.1093837015e-28
M_P = 1.67262192369e-24
HBAR = 1.054571817e-27
H_PL = 2 * np.pi * HBAR
KEV = 1.602176634e-9
EV = KEV / 1e3
ALPHA = E_CH**2 / (HBAR * C_LIGHT)
A_B = HBAR**2 / (M_E * E_CH**2)


# =====================================================================
print("\n== 2  the classical baseline")
# =====================================================================
e, m, ne = 1.3, 0.7, 2.1
w, nu = 3.0, 0.11
wp2 = 4 * np.pi * e**2 * ne / m
v_over_E = -e / (m * (nu - 1j * w))                     # from -i w m v = -e E - nu m v
sigma = -e * ne * v_over_E
report('eq:drude-sigma', 'sigma = i e^2 n/(m(w+i nu)) from Newton (11.7)', sigma, 1j * e**2 * ne / (m * (w + 1j * nu)), 1e-14)
sig_atz = e**2 * ne / (1j * m * (w + 1j * nu))           # Atzeni Eq. (11.8) as printed
report('eq:drude-sigma', 'Atzeni (11.8) as printed equals -sigma (sign slip)', sig_atz, -sigma, 1e-14)
eps_T = 1 + 4j * np.pi * sigma / w
report('eq:drude-eps', 'eps_T = 1 + 4 pi i sigma/w = 1 - wp^2/(w(w+i nu))', eps_T, 1 - wp2 / (w * (w + 1j * nu)), 1e-14)
eps_atz = 1 + 4j * np.pi * sig_atz / w
report('eq:drude-eps', 'with (11.8) as printed the plasma term flips sign', eps_atz, 1 + wp2 / (w * (w + 1j * nu)), 1e-14)
report('eq:drude-imeps', 'Im eps_T = wp^2 nu/(w(w^2+nu^2))', eps_T.imag, wp2 * nu / (w * (w**2 + nu**2)), 1e-14)
# absorbed power by integrating Newton's equation in real time to the steady state
E0 = 0.8
sol = solve_ivp(lambda t, y: [-(e / m) * E0 * np.cos(w * t) - nu * y[0]], [0, 400 * 2 * np.pi / w], [0.0],
                rtol=1e-11, atol=1e-13, dense_output=True)
Tp = 2 * np.pi / w
tt = np.linspace(398 * Tp, 399 * Tp, 4001)
pw = -e * ne * sol.sol(tt)[0] * E0 * np.cos(w * tt)
P_num = np.trapezoid(pw, tt) / Tp
v0 = e * E0 / (m * w)
report('eq:P-drude', 'cycle-averaged j.E (ODE) = (w/8pi) Im eps E0^2', P_num, w / (8 * np.pi) * eps_T.imag * E0**2, 1e-7)
report('eq:P-drude', '(w/8pi) Im eps E0^2 = ne nu m v0^2/2 /(1+nu^2/w^2)', w / (8 * np.pi) * eps_T.imag * E0**2,
       ne * nu * m * v0**2 / 2 / (1 + nu**2 / w**2), 1e-14)
# absorption coefficient from the complex wave number
cl = 50.0
for (wp_, nu_) in [(0.5, 1e-3), (2.0, 1e-3), (2.9, 1e-4)]:
    epsx = 1 - wp_**2 / (w * (w + 1j * nu_))
    kc = w / cl * np.sqrt(epsx)
    K_exact = 2 * kc.imag
    K_form = nu_ / cl * wp_**2 / w**2 / np.sqrt(1 - wp_**2 / w**2)
    report('eq:drude-K', f'K = 2 Im k vs (nu/c)(wp^2/w^2)/n_r, wp/w={wp_/w:.2f}', K_exact, K_form, 30 * (nu_ / w)**2 / (1 - wp_**2 / w**2)**2 + 1e-12)
# nu_T coefficient of Atzeni (10.133)
nuT_1keV = 4 * np.sqrt(2 * np.pi) / 3 * 1e21 * E_CH**4 / (np.sqrt(M_E) * KEV**1.5)
printed('eq:nuT', 'Atzeni (10.133): tau_e = 1.09e-11 s T^1.5/(n_i21 Z^2 lnL)', 1 / nuT_1keV * 1e11, 1.09, 2)

# =====================================================================
print("\n== 3  the Hamiltonian")
# =====================================================================
e, m, w, E0, c = 1.1, 0.9, 2.3, 0.7, 40.0
v0 = e * E0 / (m * w)
a0 = v0 / w
sol = solve_ivp(lambda t, y: [y[1], -(e / m) * E0 * np.cos(w * t)], [0, 7.0], [a0, 0.0], rtol=1e-12, atol=1e-14, dense_output=True)
tt = np.linspace(0, 7, 701)
xx, vv = sol.sol(tt)
report('eq:quiver', 'pure quiver x(t) = a0 cos wt from x(0)=a0, v(0)=0', np.max(np.abs(xx - a0 * np.cos(w * tt))) + a0, a0, 1e-9)
report('eq:quiver', 'v(t) = -v0 sin wt', np.max(np.abs(vv + v0 * np.sin(w * tt))) + v0, v0, 1e-9)
Up_num = np.trapezoid(0.5 * m * (v0 * np.sin(w * np.linspace(0, 2 * np.pi / w, 20001)))**2, np.linspace(0, 2 * np.pi / w, 20001)) / (2 * np.pi / w)
report('eq:quiver', 'U_p = <m v_osc^2/2> = m v0^2/4 = e^2E0^2/(4 m w^2)', Up_num, e**2 * E0**2 / (4 * m * w**2), 1e-10)
for t in [0.3, 1.7]:
    A = -(c * E0 / w) * np.sin(w * t)
    report('eq:alpha-A', f'alpha_dot = (e/mc) A at t={t}', -v0 * np.sin(w * t), e * A / (m * c), 1e-14)
    p = np.array([0.4, -1.2, 0.7])
    ehat = np.array([0.0, 0.6, 0.8])
    Avec = A * ehat
    adot = (e / (m * c)) * Avec
    lhs = np.sum((p + e * Avec / c)**2) / (2 * m) - adot @ p - e**2 * (Avec @ Avec) / (2 * m * c**2)
    report('eq:H-KH1', f'(p+eA/c)^2/2m - alpha_dot.p - e^2A^2/2mc^2 = p^2/2m, t={t}', lhs, p @ p / (2 * m), 1e-13)
Ze2 = 1.0
for (kk, kap) in [(0.7, 0.3), (2.5, 1.1), (0.4, 0.0)]:
    if kap > 0:
        rad = quad(lambda r: np.sin(kk * r) * np.exp(-kap * r), 0, np.inf, limit=400)[0]
    else:
        rad = 1 / kk     # Abel limit of int sin(kr) dr
    ft = -Ze2 * 4 * np.pi / kk * rad
    report('eq:vk', f'Fourier transform of -Ze^2 e^(-kap r)/r, k={kk}, kap={kap}', ft, -4 * np.pi * Ze2 / (kk**2 + kap**2), 1e-10)

# =====================================================================
print("\n== 4  the Volkov states")
# =====================================================================
hb = 1.0
p = np.array([0.6, -0.3, 1.1])
ehat = np.array([0.0, 0.6, 0.8])
Up = m * v0**2 / 4
Ep = p @ p / (2 * m)


def hPhi(t, pv=p):
    return (pv @ pv / (2 * m) + Up) * t - (pv @ (a0 * ehat)) * (1 - np.cos(w * t)) - Up / (2 * w) * np.sin(2 * w * t)


for t in [0.37, 2.2]:
    der = rich(lambda h: hPhi(t + h), 1e-3)
    A = -(c * E0 / w) * np.sin(w * t)
    rhs = np.sum((p + e * A * ehat / c)**2) / (2 * m)
    report('eq:volkov', f'hbar dPhi/dt = (p+eA/c)^2/2m at t={t}', der, rhs, 1e-11)
for (z, th) in [(0.7, 0.4), (5.3, 2.1), (19.0, -1.0)]:
    ell = np.arange(-120, 121)
    ser = np.sum((1j)**(ell % 4) * jv(ell, z) * np.exp(1j * ell * th))
    report('eq:jacobi-anger', f'sum i^l J_l(z) e^(il th) = e^(iz cos th), z={z}', ser, np.exp(1j * z * np.cos(th)), 1e-13)
    J2 = jv(ell, z)**2
    absreport('eq:bessel-sums', f'sum J_l^2 = 1, z={z}', J2.sum(), 1.0, 1e-13)
    absreport('eq:bessel-sums', f'sum l J_l^2 = 0, z={z}', (ell * J2).sum(), 0.0, 1e-13)
    report('eq:bessel-sums', f'sum l^2 J_l^2 = z^2/2, z={z}', (ell**2 * J2).sum(), z**2 / 2, 1e-12)
    for s in [0.3, 2.0]:
        report('eq:bessel-genfun', f'sum J_l^2 e^(-ils) = J0(2z sin(s/2)), z={z}, s={s}', np.sum(J2 * np.exp(-1j * ell * s)),
               jv(0, 2 * z * np.sin(s / 2)), 1e-12)
z = 1e-3
report('eq:bessel-small', 'J_1^2(z) = z^2/4 (1 + O(z^2))', jv(1, z)**2, z**2 / 4, 1e-6)
report('eq:bessel-small', 'J_0^2(z) = 1 - z^2/2 + O(z^4)', jv(0, z)**2, 1 - z**2 / 2, 1e-12)

# =====================================================================
print("\n== 5  multiphoton rates in Born approximation")
# =====================================================================
pp = np.array([-0.2, 0.9, 0.4])
kvec = (pp - p) / hb
zz = kvec @ (a0 * ehat)
ell = np.arange(-80, 81)
for t in [0.5, 3.1]:
    direct = np.exp(1j * (hPhi(t, pp) - hPhi(t, p)) / hb)
    Dl = (pp @ pp - p @ p) / (2 * m)
    series = np.exp(-1j * zz) * np.sum((1j)**(ell % 4) * jv(ell, zz) * np.exp(1j * (Dl + ell * hb * w) * t / hb))
    report('eq:Vpp', f'Volkov phase factor = e^(-iz) sum i^l J_l e^(i(dE+l hw)t), t={t}', series, direct, 1e-12)
# golden-rule limit at finite times, extrapolated in 1/T
wq, zq = 1.0, 1.7
phi = lambda d: np.exp(-((d - 0.3)**2) / 2.0) * (1 + 0.2 * d)
lq = np.arange(-30, 31)
rhs = 2 * np.pi * np.sum(jv(lq, zq)**2 * phi(-lq * wq))


def lhs_T(T, dt=0.02):
    """(1/T) int dD phi(D) |int_0^T e^{iDt} g(t) dt|^2 = (1/T) int int g(t) g*(t') phihat(t - t'),
    phihat(tau) = int phi(D) e^{iD tau} dD = e^{0.3 i tau} sqrt(2pi) e^{-tau^2/2} (1.06 + 0.2 i tau)."""
    n = int(round(T / dt))
    t = (np.arange(n) + 0.5) * dt
    g = np.exp(1j * zq * np.cos(wq * t))
    K = int(12 / dt)
    tot = 0.0 + 0.0j
    for lag in range(-K, K + 1):
        tau = lag * dt
        ph = np.exp(0.3j * tau) * np.sqrt(2 * np.pi) * np.exp(-tau * tau / 2) * (1.06 + 0.2j * tau)
        if lag >= 0:
            tot += ph * np.sum(g[lag:] * np.conj(g[:n - lag]))
        else:
            tot += ph * np.sum(g[:n + lag] * np.conj(g[-lag:]))
    return (tot * dt * dt / T).real


l1, l2 = lhs_T(2 * np.pi * 64), lhs_T(2 * np.pi * 128)        # whole periods, so the cross terms do not beat
report('eq:rate', '(1/T) int phi |I|^2 -> 2 pi sum J_l^2 phi(-l w)  (T = 128 periods)', l2, rhs, 2e-3)
report('eq:rate', 'same, extrapolated 2 L(2T) - L(T) (1/T error removed)', 2 * l2 - l1, rhs, 2e-5)
rng = np.random.default_rng(7)
Ni, L = 60, 1.0
kk = 2 * np.pi * np.array([3, -1, 2]) / L
vals = []
for _ in range(20000):
    R = rng.random((Ni, 3)) * L
    vals.append(abs(np.exp(-1j * R @ kk).sum())**2 / Ni)
vals = np.array(vals)
absreport('eq:rate-ions', '<|rho_i(k)|^2>/N_i = 1 for independent ions (Monte Carlo)', vals.mean(), 1.0, 4 * vals.std() / np.sqrt(len(vals)))
# the delta-function integral of the cross section: int p'^2 dp' delta(E'-E-l hw) = m p'_l
mm, hw0, pin, lq_ = 1.3, 0.37, 1.1, 2
pl = np.sqrt(pin**2 + 2 * mm * lq_ * hw0)
for s in [1e-3, 3e-4]:
    sq = s * mm / pl
    val = quad(lambda q: q**2 * np.exp(-((q * q - pin * pin) / (2 * mm) - lq_ * hw0)**2 / (2 * s * s)) / np.sqrt(2 * np.pi) / s,
               pl - 40 * sq, pl + 40 * sq, points=[pl], limit=400)[0]
    report('eq:bf', f"int p'^2 dp' delta(E'-E-l hw) = m p'_l (Gaussian delta, width {s})", val, mm * pl, 3e-6)
for (pv, th) in [(1.0, 0.4), (2.3, 1.9)]:
    kk_ = 2 * pv * np.sin(th / 2) / hb
    born = (mm / (2 * np.pi * hb**2))**2 * (4 * np.pi * 1.0 / kk_**2)**2
    ruth = (1.0 / (2 * mm * (pv / mm)**2))**2 / np.sin(th / 2)**4
    report('eq:rutherford', f'Born Coulomb = Rutherford, p={pv}, theta={th}', born, ruth, 1e-13)


# one-photon cross section at fixed orientation
def sig_theta0_closed(pv, ppv, th0):
    Lg = np.log((pv + ppv) / abs(pv - ppv))
    br = np.sin(th0)**2 * ((pv**2 + ppv**2) * Lg / (4 * pv * ppv) - 0.5) + np.cos(th0)**2 * (1 + (pv**2 - ppv**2) * Lg / (2 * pv * ppv))
    return 2 * np.pi / pv**2 * (ppv / pv) * br       # times m^2 Z^2 e^4 a0^2 / hbar^2


def sig_theta0_quad(pv, ppv, th0):
    eh = np.array([np.sin(th0), 0.0, np.cos(th0)])
    def f(c_, ph):
        s_ = np.sqrt(1 - c_ * c_)
        d = np.array([ppv * s_ * np.cos(ph), ppv * s_ * np.sin(ph), ppv * c_ - pv])
        return (d @ eh)**2 / (d @ d)**2
    inner = lambda c_: quad(lambda ph: f(c_, ph), 0, 2 * np.pi, limit=200)[0]
    return (ppv / pv) * quad(inner, -1, 1, limit=200, points=[1.0 - 1e-6])[0]


for (pv, ppv, th0) in [(1.0, 1.25, 0.3), (1.0, 0.7, 1.2), (2.0, 2.2, 2.6)]:
    report('eq:sigma-theta0', f"closed form vs 2-D quadrature, p={pv}, p'={ppv}, th0={th0}", sig_theta0_closed(pv, ppv, th0),
           sig_theta0_quad(pv, ppv, th0), 1e-8)
    cos0, wts0 = np.polynomial.legendre.leggauss(40)
    avg = 0.5 * np.sum(wts0 * np.array([sig_theta0_closed(pv, ppv, np.arccos(cc)) for cc in cos0]))
    report('eq:sigma-iso', f"orientation average = (2pi/3) L/p^2 (p'/p)(p/p'), p={pv}, p'={ppv}", avg,
           2 * np.pi / 3 * np.log((pv + ppv) / abs(pv - ppv)) / pv**2, 1e-13)
# Bunkin-Fedorov Eq. (10) as printed against the exact one-photon cross sections
for xi in [1e-3, 1e-4]:
    for th0 in [0.0, 0.7, np.pi / 2]:
        pv = 1.0
        for sgn, name in [(+1, 'a'), (-1, 'e')]:
            ppv = pv * np.sqrt(1 + 2 * sgn * xi)
            exact = sig_theta0_closed(pv, ppv, th0) * pv**2 / np.pi          # = 2 (p'/p) {curly bracket}
            l2x = np.log(2 / xi)
            bfp = (3 * np.cos(th0)**2 - 1) + (1 - np.cos(th0)**2) * l2x - sgn * (3 * np.cos(th0)**2 - 1) * xi * l2x + sgn * 2 * xi * np.cos(th0)**2
            report('eq:bf10', f'BF Eq.(10) sigma_{name}, xi={xi:g}, th0={th0:.2f} (residual O(xi^2 ln))', exact, bfp, 30 * xi**2 * l2x + 1e-12)
# BF prefactor conversion: n_i v hw (8 pi / c E0^2) * 2 pi m^2 Z^2 e^4 a0^2/(hbar^2 p^2) = 16 pi^2 n_i Z^2 e^6/(m^2 v c hbar w^3)
ni_, Zc, e_, m_, c_, hb_, w_, v_, E0_ = 1.7, 2.0, 0.9, 1.3, 30.0, 0.8, 2.2, 0.6, 0.45
a0_ = e_ * E0_ / (m_ * w_**2)
lhs = ni_ * v_ * hb_ * w_ * 8 * np.pi / (c_ * E0_**2) * 2 * np.pi * m_**2 * Zc**2 * e_**4 * a0_**2 / (hb_**2 * (m_ * v_)**2)
report('eq:bf10', 'prefactor in BF normalization = 16 pi^2 n_i Z^2 e^6/(m^2 v c hbar w^3)', lhs,
       16 * np.pi**2 * ni_ * Zc**2 * e_**6 / (m_**2 * v_ * c_ * hb_ * w_**3), 1e-14)
# net emission positive along the field, isotropic average -2/3
xi = 1e-4
l2x = np.log(2 / xi)
check('eq:bf10', 'net emission (3cos^2-1)ln(2/xi)-2cos^2 > 0 at th0 = 0 (Marcuse)', 2 * l2x - 2 > 0, f'{2*l2x-2:.3f}')
avg = 0.5 * quad(lambda cc: (3 * cc * cc - 1) * l2x - 2 * cc * cc, -1, 1)[0]
report('eq:bf10', 'isotropic average of the net-emission bracket = -2/3', avg, -2.0 / 3, 1e-13)
# the arcsine (classical) limit
ZE = 1.3
phi_t = lambda x: np.exp(-(x - 0.4)**2) * (1 + 0.3 * x)
target = quad(lambda th: phi_t(ZE * np.cos(th)), 0, 2 * np.pi, limit=200)[0] / (2 * np.pi)
errs = []
for hwq in [0.05, 0.02, 0.01]:
    zq = ZE / hwq
    lq = np.arange(-int(zq) - 60, int(zq) + 61)
    errs.append(abs(np.sum(jv(lq, zq)**2 * phi_t(lq * hwq)) - target) / target)
report('eq:arcsine', 'sum J_l^2(Z_E/hw) phi(l hw) -> (1/2pi) int phi(Z_E cos th), hw = 0.01', target * (1 + errs[-1]), target, 2e-4)
check('eq:arcsine', 'the error falls as hw -> 0', errs[0] > errs[1] > errs[2], ', '.join(f'{x:.1e}' for x in errs))

# =====================================================================
print("\n== 6  the kinetic equation")
# =====================================================================
M = 3
a_ = np.array([[0, 1], [0, 0]], dtype=complex)
Zp = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def kron_list(ms):
    out = ms[0]
    for x in ms[1:]:
        out = np.kron(out, x)
    return out


b = [kron_list([Zp] * j + [a_] + [I2] * (M - j - 1)) for j in range(M)]
bd = [x.conj().T for x in b]
ok_car = all(np.allclose(b[i] @ bd[j] + bd[j] @ b[i], (i == j) * np.eye(2**M)) for i in range(M) for j in range(M))
check('eq:bilinear', 'Jordan-Wigner operators obey the anticommutation relations', ok_car)
mx = 0.0
for a1 in range(M):
    for b1 in range(M):
        for c1 in range(M):
            for d1 in range(M):
                lhs = (bd[a1] @ b[b1]) @ (bd[c1] @ b[d1]) - (bd[c1] @ b[d1]) @ (bd[a1] @ b[b1])
                rhs = (b1 == c1) * (bd[a1] @ b[d1]) - (a1 == d1) * (bd[c1] @ b[b1])
                mx = max(mx, np.max(np.abs(lhs - rhs)))
absreport('eq:bilinear', '[b+a b_b, b+c b_d] = d_bc b+a b_d - d_ad b+c b_b (all 81 cases)', mx, 0.0, 1e-14)
rng = np.random.default_rng(11)
Vh = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
Vh = Vh + Vh.conj().T
H = sum(Vh[c1, d1] * bd[c1] @ b[d1] for c1 in range(M) for d1 in range(M))
X = rng.normal(size=(2**M, 2**M)) + 1j * rng.normal(size=(2**M, 2**M))
rho_mb = X @ X.conj().T
rho_mb /= np.trace(rho_mb)
rho1 = np.array([[np.trace(rho_mb @ bd[a1] @ b[b1]) for b1 in range(M)] for a1 in range(M)])
hbq = 0.7
mx = 0.0
for a1 in range(M):
    for b1 in range(M):
        O = bd[a1] @ b[b1]
        lhs = (1j / hbq) * np.trace(rho_mb @ (H @ O - O @ H))
        rhs = (1j / hbq) * sum(Vh[c1, a1] * rho1[c1, b1] - Vh[b1, c1] * rho1[a1, c1] for c1 in range(M))
        mx = max(mx, abs(lhs - rhs))
absreport('eq:eom', 'd rho_ab/dt from the Fock-space Heisenberg equation = Eq. (eom)', mx, 0.0, 1e-13)
# Markov step
phm = lambda O: np.exp(-O * O) * (1 + 0.5 * O)
for t in [200.0, 800.0]:
    val = quad(lambda O: phm(O) * np.sin(O * t) / O, -12, 12, limit=4000, points=[0.0])[0]
    report('eq:markov', f'int phi(W) sin(Wt)/W dW -> pi phi(0), t={t:.0f}', val, np.pi * phm(0.0), 2e-3 / t * 200)
# (2/hbar^2) Re[V*(t) int_0^t V] = d/dt |int_0^t V/hbar|^2
Vt = lambda t: (0.3 + 0.2j) * np.exp(1.7j * t) + 0.1 * np.exp(-0.4j * t)
intV = lambda t: quad(lambda s: Vt(s).real, 0, t)[0] + 1j * quad(lambda s: Vt(s).imag, 0, t)[0]
t = 2.3
lhs = 2 / hbq**2 * (np.conj(Vt(t)) * intV(t)).real
rhs = rich(lambda h: abs(intV(t + h) / hbq)**2, 1e-3)
report('eq:second-order', '(2/hbar^2) Re[V* int V] = d/dt |int V/hbar|^2', lhs, rhs, 1e-9)
for eta in [+1, -1]:
    f1, f2 = 0.31, 0.77
    report('eq:no-blocking', f'f2(1+eta f1) - f1(1+eta f2) = f2 - f1, eta={eta:+d}', f2 * (1 + eta * f1) - f1 * (1 + eta * f2), f2 - f1, 1e-15)
# conservation and the two forms of P on a discrete model with energies on a ladder
rng = np.random.default_rng(3)
nst, hwl = 40, 0.25
nlev = rng.integers(0, 12, nst)
En = nlev * hwl + 0.0
W = rng.random((nst, nst))
W = W + W.T
np.fill_diagonal(W, 0.0)
fr = rng.random(nst)
dfdt = W.T @ fr - W.sum(axis=0) * fr          # sum_p' W_{p'p}(f_p' - f_p)
absreport('eq:master', 'sum_p df_p/dt = 0 for symmetric W', dfdt.sum(), 0.0, 1e-13)
P1 = 2 * np.sum(En * dfdt)
P2 = np.sum((En[:, None] - En[None, :]) * W.T * (fr[None, :] - fr[:, None]))
report('eq:P-general', '(2/V) sum E_p df_p/dt = (1/V) sum (E_p - E_p\') W (f_p\' - f_p)', P1, P2, 1e-12)
fdec = np.exp(-En / 0.6)
dfd = W.T @ fdec - W.sum(axis=0) * fdec
P3 = 2 * np.sum(En * dfd)
check('eq:P-general', 'P >= 0 for f decreasing in energy (random symmetric W)', P3 >= 0, f'P={P3:.3e}')

# =====================================================================
print("\n== 7  the Lindhard function")
# =====================================================================


def imchi_direct(kk, Om, fofE):
    """-2 pi int d3q/(2 pi)^3 [f(q) - f(q+k)] delta(Om + E_q - E_{q+k}), m = hbar = 1, by cylinder quadrature."""
    def term(qpar):
        return quad(lambda qt: 2 * np.pi * qt * fofE((qpar**2 + qt**2) / 2), 0, np.inf, limit=300)[0]
    return -2 * np.pi / (2 * np.pi)**3 * (1.0 / kk) * (term(Om / kk - kk / 2) - term(Om / kk + kk / 2))


def imchi_general(kk, Om, fofE):
    Em = 0.5 * (Om / kk - kk / 2)**2
    Ep = 0.5 * (Om / kk + kk / 2)**2
    return -1 / (2 * np.pi * kk) * quad(fofE, Em, Ep, limit=300)[0]


farb = lambda E: np.exp(-E * E / 3) * (1 + 0.4 * np.cos(1.7 * E))      # neither Maxwell nor Fermi-Dirac, not monotone
for (kk, Om) in [(0.7, 0.4), (2.0, 1.3), (0.3, 2.0)]:
    report('eq:imchi-general', f'-(m^2/2pi hbar^4 k) int_E-^E+ f dE vs 3-D quadrature, arbitrary f, k={kk}', imchi_general(kk, Om, farb),
           imchi_direct(kk, Om, farb), 1e-8)
mu_nd = -3.0
fnd = lambda E: np.exp(mu_nd - E)
ne_nd = 2 * (2 * np.pi)**-3 * (2 * np.pi)**1.5 * np.exp(mu_nd)
for (kk, Om) in [(0.7, 0.4), (2.0, 1.3), (0.2, 0.05)]:
    zm, zp = (Om / kk - kk / 2) / np.sqrt(2), (Om / kk + kk / 2) / np.sqrt(2)
    cf = -ne_nd * np.sqrt(np.pi / 2) / kk * (np.exp(-zm**2) - np.exp(-zp**2))
    report('eq:imchi-mb', f'Maxwell closed form vs 3-D quadrature, k={kk}, W={Om}', cf, imchi_direct(kk, Om, fnd), 1e-10)
for (muq, Tq) in [(2.0, 0.7), (8.0, 0.7), (-1.0, 1.3)]:
    ffd = lambda E, muq=muq, Tq=Tq: 1.0 / (np.exp(np.clip((E - muq) / Tq, -700, 700)) + 1.0)
    for (kk, Om) in [(0.7, 0.4), (3.0, 5.0)]:
        Em, Ep = 0.5 * (Om / kk - kk / 2)**2, 0.5 * (Om / kk + kk / 2)**2
        cf = -(Tq / (2 * np.pi * kk)) * np.log((1 + np.exp((muq - Em) / Tq)) / (1 + np.exp((muq - Ep) / Tq)))
        report('eq:imchi-fd', f'Fermi-Dirac closed form vs 3-D quadrature, mu={muq}, T={Tq}, k={kk}', cf, imchi_direct(kk, Om, ffd), 1e-10)
# classical limit: units m = kT = 1, hbar -> 0 at fixed k, W
for (kk, Om) in [(0.8, 0.5), (0.3, 0.6)]:
    zeta = Om / (np.sqrt(2) * kk)
    cl = -1.0 * np.sqrt(np.pi) * zeta * np.exp(-zeta**2)          # n_e = kT = 1
    errs = []
    for hq in [1e-1, 1e-2, 1e-3]:
        zm, zp = (Om / kk - hq * kk / 2) / np.sqrt(2), (Om / kk + hq * kk / 2) / np.sqrt(2)
        qv = -np.sqrt(np.pi / 2) / (hq * kk) * (np.exp(-zm**2) - np.exp(-zp**2))
        errs.append(abs(qv - cl) / abs(cl))
    report('eq:imchi-classical', f'hbar -> 0 limit (hbar = 1e-3) = -(n/kT) sqrt(pi) z e^-z^2, k={kk}', cl * (1 + errs[-1]), cl, 1e-5)
    check('eq:imchi-classical', f'error falls as hbar^2, k={kk}', errs[0] / errs[1] > 50 and errs[1] / errs[2] > 50, ', '.join(f'{x:.1e}' for x in errs))

# =====================================================================
print("\n== 8  closed forms for weak fields")
# =====================================================================
for (aa, bb) in [(0.3, 2.0), (1.7, 0.01)]:
    num = quad(lambda k: np.exp(-aa / k**2 - bb * k**2) / k, 0, np.inf, limit=400)[0]
    report('eq:K0-integral', f'int dk/k e^(-a/k^2 - b k^2) = K0(2 sqrt(ab)), a={aa}, b={bb}', num, k0(2 * np.sqrt(aa * bb)), 1e-9)


def Lpm(pv, s, hwq):
    pp2 = pv * pv + 2 * s * hwq
    if pp2 <= 0:
        return 0.0
    ppv = np.sqrt(pp2)
    return np.log((pv + ppv) / abs(pv - ppv))


def nu_crosssection_maxwell(u):
    """cross-section route, units m = kT = hbar = 1, n_i Z^2 e^4 = 1; returns nu_eff/nu_T."""
    fM = lambda pv: (2 * np.pi)**-1.5 * np.exp(-pv * pv / 2)
    integ = lambda pv: 4 * np.pi * pv * pv * fM(pv) * (1 / pv) * (Lpm(pv, +1, u) - Lpm(pv, -1, u))
    val = quad(integ, 0, np.sqrt(2 * u), limit=400)[0] + quad(integ, np.sqrt(2 * u), np.inf, limit=400)[0]
    nu = (2 * u / u**2) * (2 * np.pi / 3) * val
    return nu / (4 * np.sqrt(2 * np.pi) / 3)


for u in [0.003, 0.03, 0.3, 1.0, 5.0]:
    report('eq:Lambda-mb', f'cross-section route (Bunkin-Fedorov, Maxwell) = Lambda(u), u={u}', nu_crosssection_maxwell(u), Lam(u), 1e-9)
    report('eq:Lambda-mb', f'two forms of Lambda agree, u={u}', 2 * np.sinh(u / 2) / u * k0(u / 2), Lam(u), 1e-13)
for u in [1e-3, 1e-2]:
    d = abs(Lam(u) - np.log(4 / (np.exp(EG) * u)))
    check('eq:Lambda-limits', f'|Lambda - ln(4/(e^g u))| < u^2 ln(1/u), u={u}', d < u**2 * np.log(1 / u), f'{d:.2e}')
for u in [100.0, 1000.0]:
    report('eq:Lambda-limits', f'Lambda -> sqrt(pi) u^-3/2 (error O(1/u)), u={u:g}', Lam(u), np.sqrt(np.pi) * u**-1.5, 1.0 / u)
printed('eq:cutoffs', 'b_min = e^g/4 hbar/(m v_t) = 0.445', np.exp(EG) / 4, 0.445, 3)
printed('eq:cutoffs', 'de Broglie hbar/(2 m v), m v^2 = 3kT: 0.289 hbar/(m v_t)', 1 / (2 * np.sqrt(3)), 0.289, 3)
printed('eq:cutoffs', 'ln(0.445/0.289) = 0.43', np.log(np.exp(EG) / 4 * 2 * np.sqrt(3)), 0.43, 2)
printed('eq:cutoffs', 'ratio 0.289/0.445 = 0.65', 1 / (2 * np.sqrt(3)) / (np.exp(EG) / 4), 0.65, 2)
tabL = [(0.001, 7.717, 7.717, 56050, 4.257), (0.003, 6.618, 6.618, 10787, 3.654), (0.01, 5.414, 5.414, 1772, 3.000),
        (0.03, 4.316, 4.316, 341.1, 2.415), (0.1, 3.116, 3.112, 56.05, 1.805), (0.3, 2.038, 2.013, 10.79, 1.300),
        (1, 0.963, 0.809, 1.772, 0.840), (3, 0.304, -0.290, 0.341, 0.528), (10, 0.0548, -1.494, 0.0560, 0.302)]
for (u, L1, L2, L3, g4) in tabL:
    printed('tab:Lambda', f'Lambda({u})', Lam(u), L1, 4 if u == 10 else 3)
    printed('tab:Lambda', f'ln(4/e^g u) at u={u}', np.log(4 / (np.exp(EG) * u)), L2, 3)
    dec = {0.001: 0, 0.003: 0, 0.01: 0, 0.03: 1, 0.1: 2, 0.3: 2, 1: 3, 3: 3, 10: 4}[u]
    printed('tab:Lambda', f'sqrt(pi) u^-1.5 at u={u}', np.sqrt(np.pi) * u**-1.5, L3, dec)
    printed('tab:Lambda', f'<g>({u})', gavg(u), g4, 3)


# Fermi-Dirac: closed form vs the cross-section route with explicit Pauli factors
def nu_fd_closed(u, muq):
    """units m = hbar = kT = 1, n_i Z^2 e^4 = 1 per (Z n_e = Z^2 n_i): returns nu_eff * (n_e/(Z n_i Z))..."""
    def Lfd(kk):
        Em, Ep = 0.5 * (u / kk - kk / 2)**2, 0.5 * (u / kk + kk / 2)**2
        return np.logaddexp(0, muq - Em) - np.logaddexp(0, muq - Ep)
    kw = np.sqrt(2 * u)
    integ = sum(quad(lambda kk: Lfd(kk) / kk, a, b_, limit=600)[0] for a, b_ in [(0, kw), (kw, 4 * kw + 10)]) \
        + quad(lambda kk: Lfd(kk) / kk, 4 * kw + 10, np.inf, limit=600)[0]
    return 4 / (3 * np.pi * u) * integ          # (4 Z e^4 m kT/(3 pi hbar^4 w)) with Z e^4 -> 1


def nu_fd_crosssection(u, muq):
    """2/(2pi)^3 int d3p hw v sigma_+(p) [f(1-f') - f'(1-f)] / (n_e m v0^2/2), n_i Z^2 e^4 = 1, then x n_e/(Z n_i)."""
    f = lambda E: 1 / (np.exp(np.clip(E - muq, -700, 700)) + 1)
    def integ(pv):
        E = pv * pv / 2
        f1, f2 = f(E), f(E + u)
        net = f1 * (1 - f2) - f2 * (1 - f1)
        sig = 2 * np.pi / 3 * Lpm(pv, +1, u) / pv**2           # times a0^2 m^2 Z^2 e^4/hbar^2
        return 4 * np.pi * pv * pv * pv * sig * net
    pmax = np.sqrt(2 * max(muq, 0) + 80)
    val = 2 / (2 * np.pi)**3 * quad(integ, 0, pmax, limit=600)[0]
    ne_ = 2 / (2 * np.pi)**3 * quad(lambda pv: 4 * np.pi * pv * pv * f(pv * pv / 2), 0, pmax, limit=600)[0]
    # P = hw * val * a0^2 ; nu = P/(n_e w^2 a0^2/2) = 2 val/(u n_e); Lindhard form is per n_e/(Z n_i): multiply by n_e
    return 2 * val / u


for (u, muq) in [(0.1, -4.0), (0.1, 2.0), (0.5, 10.0), (0.05, 30.0)]:
    report('eq:nu-fd', f'Fermi-Dirac closed form = cross-section route with Pauli factors, u={u}, mu/kT={muq}',
           nu_fd_closed(u, muq), nu_fd_crosssection(u, muq), 1e-7)
u = 0.2
muq = -12.0
ne_nd = 2 * (2 * np.pi)**-3 * (2 * np.pi)**1.5 * np.exp(muq)
report('eq:nu-fd', 'nondegenerate limit of the Fermi-Dirac form = Lambda(u)', nu_fd_closed(u, muq) / ne_nd / (4 * np.sqrt(2 * np.pi) / 3), Lam(u), 1e-5)


def I0w(w):
    s1 = np.sqrt(1 + 2 * w)
    Im_ = (1 + w) * np.log((s1 + 1) / (s1 - 1)) - s1
    Ip_ = 0.0
    if w < 0.5:
        s2 = np.sqrt(1 - 2 * w)
        Ip_ = (1 - w) * np.log((1 + s2) / (1 - s2)) - s2
    return Im_ - Ip_


for w in [0.01, 0.1, 0.3, 0.6, 2.0, 10.0]:
    f = lambda x: (max(1 - (w / x - x / 2)**2, 0) - max(1 - (w / x + x / 2)**2, 0)) / x
    pts = sorted([np.sqrt(1 + 2 * w) - 1, np.sqrt(1 + 2 * w) + 1] + ([1 - np.sqrt(1 - 2 * w), 1 + np.sqrt(1 - 2 * w)] if w < 0.5 else []))
    num = sum(quad(f, a, b_, limit=200)[0] for a, b_ in zip([1e-9] + pts, pts + [12.0]))
    report('eq:nu-T0', f'I_-(w) - I_+(w) vs quadrature of the T=0 integral, w={w}', I0w(w), num, 1e-10)
# T = 0 against the finite-temperature form at kT = E_F/2000 (units m = hbar = 1, E_F = 1/2 -> v_F = k_F = 1)
EF = 0.5
for w in [0.1, 0.8]:
    hw_ = 2 * EF * w
    Tq = EF / 2000
    def Lfd(kk):
        Em, Ep = 0.5 * (hw_ / kk - kk / 2)**2, 0.5 * (hw_ / kk + kk / 2)**2
        return Tq * (np.logaddexp(0, (EF - Em) / Tq) - np.logaddexp(0, (EF - Ep) / Tq))
    pts = sorted(set([np.sqrt(1 + 2 * w) - 1, np.sqrt(1 + 2 * w) + 1] + ([1 - np.sqrt(1 - 2 * w), 1 + np.sqrt(1 - 2 * w)] if w < 0.5 else [])))
    num = sum(quad(lambda kk: Lfd(kk) / kk, a, b_, limit=400)[0] for a, b_ in zip([1e-9] + pts, pts + [12.0]))
    # nu = (4 Z e^4 m/(3 pi hbar^4 w)) int dk/k [..] with E_F factor:  L0 = I0/(2w) = int/(EF * 2w)  (x = k/k_F, k_F = 1)
    report('eq:nu-T0', f'T=0 closed form = Fermi-Dirac form at kT = E_F/2000, w={w}', I0w(w) / (2 * w), num / (EF * 2 * w), 2e-4)
for w in [1e-4, 1e-3]:
    d = abs(I0w(w) / (2 * w) - np.log(2 / w))
    check('eq:L0-limits', f'L0 -> ln(4E_F/hw) = ln(2/w), error < 2 w ln(1/w), w={w}', d < 2 * w * np.log(1 / w), f'{d:.1e}')
for w in [1e2, 1e3]:
    report('eq:L0-limits', f'L0 -> (4/3)(2w)^(-3/2), w={w:g}', I0w(w) / (2 * w), 4 / 3 * (2 * w)**-1.5, 3.0 / w)
# Ziman
kF, kmin = 1.0, 0.013
val = quad(lambda q: q**3 * (4 * np.pi / q**2)**2, kmin, 2 * kF)[0]
ni_ = 1.0 / (3 * np.pi**2)              # Z = 1, n_e = k_F^3/(3 pi^2)
report('eq:L0-limits', 'Ziman 1/tau = (n_i m/4 pi hbar^3 k_F^3) int q^3 |v|^2 = (4 m/3 pi) ln(2k_F/k_min)',
       ni_ / (4 * np.pi * kF**3) * val, 4 / (3 * np.pi) * np.log(2 * kF / kmin), 1e-12)
report('eq:L0-limits', 'ln(2k_F/k_min) with k_min = w/v_F equals ln(4E_F/hw)', np.log(2 * 1.0 / (0.02 / 1.0)), np.log(4 * 0.5 / 0.02), 1e-14)
# universal tail: Maxwell, T = 0 and a shell distribution
uu = 1e5
report('eq:tail', 'Maxwell nu_T sqrt(pi) u^-3/2 = (4 sqrt2 pi/3) Z n_e e^4/(m^1/2 (hw)^3/2)', 4 * np.sqrt(2 * np.pi) / 3 * np.sqrt(np.pi) * uu**-1.5,
       4 * np.sqrt(2) * np.pi / 3 * uu**-1.5, 1e-14)
report('eq:tail', 'Maxwell Lambda(u) -> sqrt(pi) u^-3/2 at u = 1e5', Lam(uu), np.sqrt(np.pi) * uu**-1.5, 1e-4)
w = 1e5
nuT0 = 4 / (3 * np.pi) * I0w(w) / (2 * w)          # (4 Z e^4 m/(3 pi hbar^3)) L0 with Z e^4 = m = hbar = 1, E_F = 1/2, n_e = 1/(3 pi^2)
report('eq:tail', 'T=0 form -> universal tail at hw = 2e5 E_F', nuT0, 4 * np.sqrt(2) * np.pi / 3 * (1 / (3 * np.pi**2)) * (2 * EF * w)**-1.5, 2e-4)
# a shell distribution f = exp(-(E-1)^2/0.02) (m = hbar = 1), its n_e, and the tail via Eq. (imchi-general) at large w
fsh = lambda E: np.exp(-(E - 1.0)**2 / 0.02)
ne_sh = 2 / (2 * np.pi)**3 * quad(lambda pv: 4 * np.pi * pv * pv * fsh(pv * pv / 2), 0, 3, points=[np.sqrt(2.0)])[0]


def imchi_shell(kk, Om):
    Em, Ep = 0.5 * (Om / kk - kk / 2)**2, 0.5 * (Om / kk + kk / 2)**2
    integral = np.sqrt(0.02 * np.pi) / 2 * (erf((Ep - 1.0) / np.sqrt(0.02)) - erf((Em - 1.0) / np.sqrt(0.02)))
    return -1 / (2 * np.pi * kk) * integral


for wl in [2e3, 2e4]:
    kw = np.sqrt(2 * wl)
    intk = quad(lambda kk: imchi_shell(kk, wl), kw - 4.0, kw + 4.0, points=[kw - 1.5, kw, kw + 1.5], limit=400)[0]
    nu_sh = -8 / (3 * wl) * intk                     # -(8 Z e^4/(3 m w)) int Im chi0 dk, Z e^4 = 1
    report('eq:tail', f'shell distribution -> universal tail, w={wl:g} (error O(E/hw))', nu_sh,
           4 * np.sqrt(2) * np.pi / 3 * ne_sh * wl**-1.5, 10.0 / wl)

# =====================================================================
print("\n== 9  spontaneous bremsstrahlung, the Gaunt factor and Kirchhoff's law")
# =====================================================================
hbq, wq, Vq, Nq = 1.0, 1.3, 2.0, 3.0
E0q = 2 * np.sqrt(2 * np.pi * hbq * wq * Nq / Vq)                      # E0/2 = (2 pi hbar w N/V)^1/2
report('eq:spont-stim', 'E0/2 = (2 pi hbar w N/V)^1/2 gives E0^2/8pi = N hbar w/V', E0q**2 / (8 * np.pi), Nq * hbq * wq / Vq, 1e-14)
# polarization sum and photon-direction integral, by sampling
kv = np.array([0.3, -0.8, 0.5])
cosl, wl_ = np.polynomial.legendre.leggauss(60)
phil = np.linspace(0, 2 * np.pi, 121)[:-1]
tot = 0.0
for cl_, wc in zip(cosl, wl_):
    sl = np.sqrt(1 - cl_ * cl_)
    for ph in phil:
        kL = np.array([sl * np.cos(ph), sl * np.sin(ph), cl_])
        e1 = np.cross(kL, [0.3, 0.1, 0.9])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(kL, e1)
        tot += wc * (2 * np.pi / len(phil)) * ((kv @ e1)**2 + (kv @ e2)**2)
report('eq:spont-diff', 'int dOmega_L sum_lambda (k.e_lambda)^2 = (8 pi/3) k^2', tot, 8 * np.pi / 3 * (kv @ kv), 1e-12)
# Born cross section vs Kramers: e2, c, hbar, m, Z arbitrary units
e2, cc_, hbq, mq, Zq = 1.0, 7.0, 1.0, 1.0, 1.0
for (pv, hwq) in [(1.0, 0.3), (2.0, 1.9), (3.0, 0.01)]:
    ppv = np.sqrt(pv * pv - 2 * mq * hwq)
    wf = hwq / hbq
    integ = lambda c_: 2 * np.pi * (e2 / (6 * np.pi**3 * cc_**3 * hbq**3 * wf)) * (ppv / pv) * ((pv * pv + ppv * ppv - 2 * pv * ppv * c_) / hbq**2) \
        * (4 * np.pi * Zq * e2 / ((pv * pv + ppv * ppv - 2 * pv * ppv * c_) / hbq**2))**2
    ds = quad(integ, -1, 1, limit=200)[0]
    born = 16 / 3 * Zq**2 * e2**3 / (mq**2 * cc_**3 * hbq * (pv / mq)**2 * wf) * np.log((pv + ppv) / (pv - ppv))
    report('eq:born-brems', f'mode-summed Eq. (spont-diff) integrated over Omega\' = Eq. (born-brems), p={pv}', ds, born, 1e-9)
    kram = (16 * np.pi / (3 * np.sqrt(3))) * Zq**2 * e2**3 / (mq**2 * cc_**3 * (pv / mq)**2 * hbq) / wf
    report('eq:gaunt-born', f'Born/Kramers = (sqrt3/pi) ln((p+p\')/(p-p\')), p={pv}', born / kram, np.sqrt(3) / np.pi * np.log((pv + ppv) / (pv - ppv)), 1e-13)
for u in [0.01, 0.3, 1.0, 4.0]:
    gB = lambda eps: np.sqrt(3) / np.pi * np.log((np.sqrt(eps) + np.sqrt(eps - u)) / (np.sqrt(eps) - np.sqrt(eps - u)))
    num = np.exp(u) * quad(lambda eps: np.exp(-eps) * gB(eps), u, np.inf, limit=300)[0]
    report('eq:gaunt-avg', f'thermal average of g_B = (sqrt3/pi) e^(u/2) K0(u/2), u={u}', num, gavg(u), 1e-9)
    report('eq:Lambda-gaunt', f'Lambda = (pi/sqrt3) <g> (1-e^-u)/u, u={u}', np.pi / np.sqrt(3) * gavg(u) * (1 - np.exp(-u)) / u, Lam(u), 1e-13)
gbar = quad(lambda u: gavg(u) * np.exp(-u), 0, np.inf, limit=400)[0]
report('eq:gaunt-power', 'int <g>(u) e^-u du = 2 sqrt3/pi', gbar, 2 * np.sqrt(3) / np.pi, 1e-10)
report('eq:gaunt-power', 'int dt/(1+cosh t) = 1 (so int e^-x K0 = 1)', quad(lambda t: 1 / (1 + np.cosh(t)), 0, np.inf)[0], 1.0, 1e-12)
printed('eq:gaunt-power', '2 sqrt3/pi = 1.103', 2 * np.sqrt(3) / np.pi, 1.103, 3)
# Kirchhoff, CGS, at T = 300 eV and hw = 30 eV, n_e = n_i = 1e20, Z = 1 (vacuum-like, n_r = 1)
kT, hw = 300 * EV, 30 * EV
u = hw / kT
ne_, ni_, Zk = 1e20, 1e20, 1.0
vt = np.sqrt(kT / M_E)
wq = hw / HBAR


def avg_L_over_v():
    """Maxwell average of L/v over the speed, with v = sqrt(v_min^2 + y^2) to smooth the threshold."""
    fM = lambda v: 4 * np.pi * v * v * (M_E / (2 * np.pi * kT))**1.5 * np.exp(-M_E * v * v / (2 * kT))
    vmin = np.sqrt(2 * hw / M_E)
    def integ(y):
        v = np.sqrt(vmin * vmin + y * y)
        p_ = M_E * v
        pp_ = M_E * y                                   # sqrt(p^2 - 2 m hw) = m y
        return fM(v) * np.log((p_ + pp_) / (p_ - pp_)) / v * (y / v)     # dv = (y/v) dy
    return quad(integ, 0, 30 * vt, limit=400, epsabs=0, epsrel=1e-12)[0]


jw = 1 / (4 * np.pi) * ne_ * ni_ * hw * 16 / 3 * Zk**2 * E_CH**6 / (M_E**2 * C_LIGHT**3 * HBAR * wq) * avg_L_over_v()
nuT_k = 4 * np.sqrt(2 * np.pi) / 3 * ni_ * Zk**2 * E_CH**4 / (np.sqrt(M_E) * kT**1.5)
wp2 = 4 * np.pi * E_CH**2 * ne_ / M_E
Kk = nuT_k * Lam(u) / C_LIGHT * wp2 / wq**2
Bw = HBAR * wq**3 / (4 * np.pi**3 * C_LIGHT**2) / (np.exp(u) - 1)
report('eq:kirchhoff', 'j_w (Born emission, Maxwell quadrature) = K B_w (master-equation absorption)', jw, Kk * Bw, 1e-8)
report('eq:kirchhoff', '<L/v> = (2/pi)^1/2 (pi/sqrt3) <g> e^-u / v_t', avg_L_over_v(), np.sqrt(2 / np.pi) * np.pi / np.sqrt(3) * gavg(u) * np.exp(-u) / vt, 1e-8)

# =====================================================================
print("\n== 10  the interacting electron gas and dynamic screening")
# =====================================================================
# two-frame identity with a random periodic P(t)
wq = 1.3
tt = np.linspace(0, 2 * np.pi / wq, 40001)
coef = [(0.4, 0.2), (0.1, -0.3), (0.05, 0.07)]
Pt = sum(a * np.cos((n + 1) * wq * tt) + b_ * np.sin((n + 1) * wq * tt) for n, (a, b_) in enumerate(coef)) + 0.3
dPt = sum(-(n + 1) * wq * a * np.sin((n + 1) * wq * tt) + (n + 1) * wq * b_ * np.cos((n + 1) * wq * tt) for n, (a, b_) in enumerate(coef))
eq_, mq_, E0q, Nq = 1.0, 1.0, 0.7, 5.0
Et = E0q * np.cos(wq * tt)
adot = -(eq_ * E0q / (mq_ * wq)) * np.sin(wq * tt)
Plab = -eq_ * Et * (Pt / mq_ + Nq * adot)
dH = -adot * dPt
report('eq:frames', '<P_lab> = <d<H\'>/dt> for a periodic total momentum', np.trapezoid(Plab, tt), np.trapezoid(dH, tt), 1e-9)
# chi0 for Maxwell electrons: Z form vs principal-value quadrature (units m = kT = 1, hbar = h)
for (kk, Om, hq) in [(0.7, 0.9, 0.3), (2.0, 1.1, 1.0), (0.4, 0.2, 0.05)]:
    cm, cp = Om / kk - hq * kk / 2, Om / kk + hq * kk / 2
    g1 = lambda v: np.exp(-v * v / 2) / np.sqrt(2 * np.pi)
    pv_m = quad(g1, -40, 40, weight='cauchy', wvar=cm, limit=400)[0]        # PV int g1/(v - c)
    pv_p = quad(g1, -40, 40, weight='cauchy', wvar=cp, limit=400)[0]
    # 1/(c - v + i0) = -PV 1/(v - c) - i pi delta(v - c)
    I_m = -pv_m - 1j * np.pi * g1(cm)
    I_p = -pv_p - 1j * np.pi * g1(cp)
    direct = (1.0 / (hq * kk)) * (I_m - I_p)          # n_e = 1
    zform = 1.0 / (np.sqrt(2) * hq * kk) * (Zf(cp / np.sqrt(2)) - Zf(cm / np.sqrt(2)))
    report('eq:chi0-Z', f'chi0 = n/(sqrt2 hbar k v_t)[Z(z+) - Z(z-)] vs principal-value quadrature, k={kk}', zform, direct, 1e-9)
# the dielectric identity with RPA
for chi0v in [(-0.3 - 0.12j), (0.8 - 0.05j)]:
    kq, eq2, Zq = 0.9, 1.0, 2.0
    vc = 4 * np.pi * eq2 / kq**2
    epsv = 1 - vc * chi0v
    chi = chi0v / epsv
    report('eq:P-silin', '|v|^2 Im chi = (4 pi Z^2 e^2/k^2) Im(1/eps), RPA', (4 * np.pi * Zq * eq2 / kq**2)**2 * chi.imag,
           4 * np.pi * Zq**2 * eq2 / kq**2 * (1 / epsv).imag, 1e-13)
    report('eq:P-silin', '1/eps = 1 + v_c chi in RPA', 1 / epsv, 1 + vc * chi, 1e-13)


def eps_x(x, r, u):
    zp = (1 / x + u * x / 2) / np.sqrt(2)
    zm = (1 / x - u * x / 2) / np.sqrt(2)
    return 1 - r**2 * (Zf(zp) - Zf(zm)) / (np.sqrt(2) * u * x**3)


def Lam_rpa(r, u, screen=True, analytic_pole=None):
    """Eq. (Lambda-rpa).  Near the Langmuir zero of Re eps the integrand is a Lorentzian of width
    Im eps/|d Re eps/dx|, which falls below double precision as r -> 1.  Inside the window
    x_pl (1 -+ 1e-3) the Lorentzian is integrated analytically, -x_pl^2 (2/|s|) arctan(|s| d/g0),
    when its width is below 1e-4 of the window, and numerically otherwise."""
    def f(x):
        ee = eps_x(x, r, u)
        return x**2 * ((1 / ee).imag if screen else -ee.imag)
    xmax = 60.0 / u + 60.0
    edges = [1e-8, 0.05, 0.2, 0.5, 1, 2, 5, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4]
    edges = [q for q in edges if q < xmax] + [xmax]
    inner = 0.0
    if screen and r < 1:
        g = lambda x: eps_x(x, r, u).real
        x0 = np.sqrt((1 / r**2 - 1) / 3)
        lo, hi = 0.3 * x0, min(1.6 * x0, 0.9)
        if g(lo) * g(hi) < 0:
            pole = brentq(g, lo, hi, xtol=1e-15)
            s_ = (g(pole * (1 + 1e-6)) - g(pole * (1 - 1e-6))) / (2e-6 * pole)
            g0 = eps_x(pole, r, u).imag
            width = abs(g0 / s_)
            d = 1e-3 * pole
            use_analytic = (width < 1e-4 * d) if analytic_pole is None else analytic_pole
            if use_analytic:
                inner = -pole**2 * (2 / abs(s_)) * np.arctan(abs(s_) * d / g0)
            else:
                pts = sorted(set([pole] + [pole + sg * kq * width for sg in (-1, 1) for kq in (1, 3, 10, 30, 100, 1e3, 1e4)
                                           if kq * width < d]))
                ie = [pole - d] + pts + [pole + d]
                inner = sum(quad(f, a, b_, limit=500, epsabs=0, epsrel=1e-12)[0] for a, b_ in zip(ie[:-1], ie[1:]))
            edges = sorted(set([q for q in edges if not (pole - d < q < pole + d)] + [pole - d, pole + d]))
            edges = [q for q in edges]
            tot = 0.0
            for a, b_ in zip(edges[:-1], edges[1:]):
                if a == pole - d and b_ == pole + d:
                    continue
                tot += quad(f, a, b_, limit=500, epsabs=0, epsrel=1e-12)[0]
            return -(1 / r**2) * np.sqrt(2 / np.pi) * (tot + inner)
    tot = 0.0
    for a, b_ in zip(edges[:-1], edges[1:]):
        tot += quad(f, a, b_, limit=500, epsabs=0, epsrel=1e-12)[0]
    return -(1 / r**2) * np.sqrt(2 / np.pi) * tot


for u in [0.03, 0.3]:
    report('eq:Lambda-rpa', f'unscreened limit of Eq. (Lambda-rpa) = Lambda(u), u={u}', Lam_rpa(0.5, u, screen=False), Lam(u), 1e-9)
# dimensional Dawson-Oberman = dimensionless form (CGS check at one point)
kT, ne_, Zq = 500 * EV, 3e21, 1.0
vt = np.sqrt(kT / M_E)
wp = np.sqrt(4 * np.pi * E_CH**2 * ne_ / M_E)
wq = wp / 0.6
u = HBAR * wq / kT
xg = lambda k: k * vt / wq
nu_dim = -(2 * Zq * E_CH**2 / (3 * np.pi * M_E * wq)) * (wq / vt)**3 * quad(lambda x: x**2 * (1 / eps_x(x, 0.6, u)).imag, 1e-8, 60 / u + 60,
                                                                             limit=2000, epsabs=0, epsrel=1e-11, points=[0.5, 1, 2, 5, 20])[0]
nuT_d = 4 * np.sqrt(2 * np.pi) / 3 * (ne_ / Zq) * Zq**2 * E_CH**4 / (np.sqrt(M_E) * kT**1.5)
report('eq:dawson-oberman', 'Dawson-Oberman in CGS / nu_T = Lambda_RPA(r = 0.6, u)', nu_dim / nuT_d, Lam_rpa(0.6, u), 1e-7)
# classical limit of eps = Atzeni Eq. (11.19)
for (x, r) in [(0.8, 0.7), (2.0, 0.4)]:
    u = 1e-5
    zeta = 1 / (np.sqrt(2) * x)
    Wat = 1 + (np.sqrt(2) * zeta / np.sqrt(2)) * Zf(zeta)          # W(X) = 1 + (X/sqrt2) Z(X/sqrt2), X = w/(k v_t) = sqrt2 zeta
    eps_cl = 1 + r**2 / x**2 * Wat
    report('eq:chi0-Z', f'u -> 0 limit of eps = 1 + (k_D^2/k^2) W(w/k v_t), Atzeni (11.19), x={x}', eps_x(x, r, u), eps_cl, 1e-6)
tabR = {0.003: [-0.440, -0.375, -0.278, -0.204, -0.099, -0.034, -0.004],
        0.03: [-0.439, -0.374, -0.277, -0.203, -0.098, -0.034, -0.004],
        0.3: [-0.393, -0.333, -0.244, -0.178, -0.085, -0.029, -0.003]}
rlist = [0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1]
RPA_CACHE = {}
for u, row in tabR.items():
    for r, shown in zip(rlist, row):
        d = Lam_rpa(r, u) - Lam(u)
        RPA_CACHE[(r, u)] = d
        printed('tab:rpa', f'Lambda_RPA - Lambda at u={u}, r={r}', d, shown, 3)
d97 = Lam_rpa(0.97, 0.003) - Lam(0.003)
printed('tab:rpa', 'correction near the critical density, r = 0.97: about -0.47', d97, -0.47, 2)
printed('tab:rpa', 'r = 0.97 is n_e = 0.94 n_c', 0.97**2, 0.94, 2)
printed('tab:rpa', 'correction at n_e = n_c/4 (r = 0.5, u = 0.003): -0.10', RPA_CACHE[(0.5, 0.003)], -0.10, 2)
printed('tab:rpa', 'limit w -> w_p (r = 0.99999, u = 0.003): about -0.56', Lam_rpa(0.99999, 0.003) - Lam(0.003), -0.56, 2)
printed('tab:rpa', 'limit w -> w_p at u = 0.3: about -0.51', Lam_rpa(0.99999, 0.3) - Lam(0.3), -0.51, 2)
# the analytic treatment of the Langmuir pole against the resolved numerical one where both work
for (r, u) in [(0.95, 0.003), (0.95, 0.3), (0.93, 0.03)]:
    report('eq:Lambda-rpa', f'Langmuir pole: analytic Lorentzian = resolved quadrature, r={r}, u={u}',
           Lam_rpa(r, u, analytic_pole=True), Lam_rpa(r, u, analytic_pole=False), 1e-6)
for r in [0.97, 0.98, 0.99]:
    a_, b_ = Lam_rpa(r, 0.003) - Lam(0.003), Lam_rpa(r, 0.03) - Lam(0.03)
    absreport('eq:Lambda-rpa', f'near the critical density, r={r}: u = 0.003 and 0.03 agree as at r <= 0.95', a_, b_, 3e-3)

# =====================================================================
print("\n== 11  strong fields")
# =====================================================================


def I_all(lmax, z):
    """I_l(z) = int_0^z J_l^2, l = 0..lmax, by the downward recursion of Eq. (Bl-recursion)."""
    z = np.atleast_1d(z)
    Lst = int(max(np.max(z) + 40 + 10 * np.max(z)**(1 / 3), lmax + 20))
    out = np.zeros((lmax + 1, z.size))
    Il = np.zeros(z.size)
    Jn1 = jv(Lst + 1, z)
    for l in range(Lst, -1, -1):
        Jl = jv(l, z)
        Il = Il + z * (Jl**2 + Jn1**2) / (2 * l + 1)
        if l <= lmax:
            out[l] = Il
        Jn1 = Jl
    return out


for z in [0.5, 7.3, 40.0]:
    Ia = I_all(12, z)
    for l in [0, 3, 11]:
        q_ = quad(lambda t: jv(l, t)**2, 0, z, limit=400)[0]
        report('eq:Bl-recursion', f'downward recursion = quadrature of int_0^z J_l^2, z={z}, l={l}', Ia[l, 0], q_, 1e-10)
    Ib = I_all(int(z) + 60, z)
    report('eq:Bl-recursion', f'I_0 + 2 sum_l>=1 I_l = z (from sum J_l^2 = 1), z={z}', Ib[0, 0] + 2 * Ib[1:, 0].sum(), z, 1e-12)


def S_exact(beta, u, nx=5000):
    lmax = int(2.5 * 2 * beta**2 / u + 60)
    xmax = max(1.6 * np.sqrt(2 * lmax / u) + 12 / u, 40 / u) + 20
    x = np.geomspace(1e-2, xmax, nx)
    I = I_all(lmax, x * beta)
    tot = 0.0
    for l in range(1, lmax + 1):
        D = np.exp(-0.5 * (l / x - u * x / 2)**2) - np.exp(-0.5 * (l / x + u * x / 2)**2)
        B = 2 / (x * beta) * I[l]
        tot += 2 * l * np.trapezoid(D * B / x**2, np.log(x))      # dx/x^3 = dlnx/x^2
    return tot


def R_exact(beta, u, nx=5000):
    return 3 * S_exact(beta, u, nx) / (u * beta**2 * Lam(u))


u = 0.3
report('eq:S-exact', 'weak-field limit of S: R(beta = 0.01) -> 1', R_exact(0.01, u), 1.0, 1e-4)
S1 = R_exact(2.0, u, 5000)
S2 = R_exact(2.0, u, 10000)
report('eq:S-exact', 'x-grid convergence of R (5000 vs 10000 points), beta=2', S1, S2, 1e-6)
tabS = [(0.05, 0.9996, 0.9994), (0.25, 0.9909, 0.9861), (1, 0.8673, 0.8050), (2, 0.5909, 0.4604), (4, 0.2102, 0.1242), (8, 0.0440, 0.0224)]


def g_fun(s):
    return np.sqrt(np.pi) * erf(s) / s - 2 * np.exp(-s * s) if s > 1e-4 else 4 / 3 * s * s - 4 / 5 * s**4


def avg_g(s0):
    return 2 / np.pi * quad(lambda th: g_fun(s0 * np.sin(th)), 0, np.pi / 2, limit=400, epsabs=0, epsrel=1e-13)[0]


def Gf(s0):
    return 3 * avg_g(s0) / (2 * s0**2)


def Qf(r):
    return quad(lambda z: z * z * (ive(0, z * z) - ive(1, z * z)), 0, r, limit=400, epsabs=0, epsrel=1e-13)[0]


for (beta, Rs, Gs) in tabS:
    printed('tab:strong', f'R(beta={beta}, u=0.3)', R_exact(beta, 0.3), Rs, 4)
    printed('tab:strong', f'G(beta/sqrt2), beta={beta}', Gf(beta / np.sqrt(2)), Gs, 4)
# the l-th term of the Lindhard route against the cross-section route (4-D quadrature), beta = 1, u = 0.5


def lindhard_term(l, beta, u):
    x = np.geomspace(1e-3, 60, 20000)
    I = I_all(l, x * beta)
    D = np.exp(-0.5 * (l / x - u * x / 2)**2) - np.exp(-0.5 * (l / x + u * x / 2)**2)
    B = 2 / (x * beta) * I[l]
    return 3 * 2 * l * np.trapezoid(D * B / x**2, np.log(x)) / (u * beta**2)        # contribution to nu/nu_T


def crosssection_term(l, beta, u, ny=48, nmu=24, nk=96, nphi=48):
    """2 l u [<v sigma_l> - <v sigma_-l>]/beta^2 / nu_T, units m = hbar = kT = 1, n_i Z^2 e^4 = 1, a0 = beta/u.

    The Maxwell average int_0^inf e^-t (...) dt, t = p^2/2, is done with t = t0 + y^2 (t0 = 0 for absorption,
    t0 = l u for emission, where sigma switches on like sqrt), so that the integrand is smooth in y."""
    a0q = beta / u
    ys, yw = np.polynomial.legendre.leggauss(ny)
    ys = 0.5 * 6.5 * (ys + 1)
    yw = 0.5 * 6.5 * yw
    mu0, mw0 = np.polynomial.legendre.leggauss(nmu)       # cos(theta_0) in (0,1), by the symmetry mu0 -> -mu0, phi -> phi + pi
    mu0 = 0.5 * (mu0 + 1)
    mw0 = 0.5 * mw0
    xs, xw = np.polynomial.legendre.leggauss(nk)
    phis = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    res = 0.0
    for sgn in (+1, -1):
        acc = 0.0
        t0 = 0.0 if sgn > 0 else l * u
        for y, wy in zip(ys, yw):
            t = t0 + y * y
            wt = wy * 2 * y * np.exp(-t)
            pv = np.sqrt(2 * t)
            pp2 = pv * pv + 2 * sgn * l * u
            if pp2 <= 0 or pv == 0:
                continue
            ppv = np.sqrt(pp2)
            k2lo, k2hi = (pv - ppv)**2, (pv + ppv)**2
            lk = 0.5 * (np.log(k2hi) - np.log(k2lo)) * xs + 0.5 * (np.log(k2hi) + np.log(k2lo))
            k2 = np.exp(lk)
            jac = 0.5 * (np.log(k2hi) - np.log(k2lo)) * xw * k2 / (2 * pv * ppv)      # dc = dk^2/(2 p p')
            cth = (pv * pv + ppv * ppv - k2) / (2 * pv * ppv)
            sth = np.sqrt(np.clip(1 - cth * cth, 0, None))
            sig = 0.0
            for m0, mw in zip(mu0, mw0):
                s0 = np.sqrt(1 - m0 * m0)
                kdote = (ppv * sth[:, None] * np.cos(phis)[None, :] * s0 + (ppv * cth[:, None] - pv) * m0) * a0q
                val = jv(l, kdote)**2 * 16 * np.pi**2 / k2[:, None]**2 / (4 * np.pi**2)
                sig += mw * np.sum(jac[:, None] * val) * (2 * np.pi / len(phis)) * (ppv / pv) * 2      # mu0 symmetric: x2 on (0,1) half
            sig *= 0.5     # average over mu0 in (-1,1)
            # <v sigma> with Maxwell: int 4 pi p^2 dp (2pi)^-1.5 e^-p^2/2 p sigma ; p dp = dt  ->  4 pi (2pi)^-1.5 int dt e^-t p^2 sigma
            acc += wt * 4 * np.pi * (2 * np.pi)**-1.5 * pv * pv * sig
        res += sgn * acc
    nu = 2 * l * u * res / beta**2
    return nu / (4 * np.sqrt(2 * np.pi) / 3)


for l in [1, 2]:
    report('eq:S-exact', f'term l={l}: Lindhard route (Eq. S-exact) = 4-D cross-section route, beta=0.8, u=0.8',
           lindhard_term(l, 0.8, 0.8), crosssection_term(l, 0.8, 0.8), 2e-5)
# the quasi-static (arcsine) limit pointwise in x
beta = 3.0
s0 = beta / np.sqrt(2)
errs = []
for x in [3.0, 6.0, 12.0, 24.0]:
    Nl = int(6 * x * beta + 50)
    I = I_all(Nl, np.array([x * beta]))[:, 0]
    tot = 0.0
    for l in range(1, Nl + 1):
        tot += 2 * l * l * np.exp(-l * l / (2 * x * x)) * 2 / (x * beta) * I[l]
    errs.append(abs(tot / (x**2 * avg_g(s0)) - 1))
report('eq:quasi-static', 'classical l-sum -> x^2 <g(s0 cos th)> at x = 24 (beta = 3)', 1 + errs[-1], 1.0, 1e-4)
check('eq:quasi-static', 'error falls as 1/x^2 (ratios ~4 per doubling)', all(3.5 < errs[i] / errs[i + 1] < 4.5 for i in range(3)),
      ', '.join(f'{q:.2e}' for q in errs))
# g, G, Silin
for s in [0.3, 1.0, 4.0]:
    report('eq:g-function', f'g(s) = 4 s^2 int_0^1 t^2 e^(-s^2 t^2) = sqrt(pi) erf(s)/s - 2 e^-s^2, s={s}',
           4 * s * s * quad(lambda t: t * t * np.exp(-s * s * t * t), 0, 1)[0], g_fun(s), 1e-12)
    report('eq:g-function', f'angular integral 2 pi int mu^2 s^2 e^(-s^2 mu^2) dmu = pi g(s), s={s}',
           2 * np.pi * quad(lambda mu_: s * s * mu_ * mu_ * np.exp(-s * s * mu_ * mu_), -1, 1)[0], np.pi * g_fun(s), 1e-12)
s = 0.05
report('eq:g-function', 'series g = 4/3 s^2 - 4/5 s^4 + 2/7 s^6', 4 / 3 * s**2 - 4 / 5 * s**4 + 2 / 7 * s**6, g_fun(s), 1e-8)
for r in [0.1, 0.5, 1, 2, 5, 20, 100]:
    report('eq:silin-Q', f'G(sqrt2 r) = 3 Q(r)/r^3, r={r}', Gf(np.sqrt(2) * r), 3 * Qf(r) / r**3, 1e-11)
r = 0.7
report('eq:silin-Q', "Silin's prefactor 32 sqrt(2pi)/(8 r^3) = (4 sqrt(2pi)/3)(3/r^3)", 32 * np.sqrt(2 * np.pi) / (8 * r**3),
       4 * np.sqrt(2 * np.pi) / 3 * 3 / r**3, 1e-14)
report('eq:G', 'G(s0 -> 0) = 1', Gf(1e-3), 1.0, 1e-6)
s0 = 0.08
report('eq:G-limits', 'G = 1 - 9/20 s0^2 + 15/112 s0^4 + O(s0^6), s0 = 0.08', Gf(s0), 1 - 9 / 20 * s0**2 + 15 / 112 * s0**4, 5 * s0**6)
for s0 in [30.0, 300.0, 3000.0]:
    asy = 3 / (np.sqrt(np.pi) * s0**3) * (np.log(4 * s0) + EG / 2 - 1)
    report('eq:G-limits', f'G -> 3/(sqrt(pi) s0^3)[ln(4 s0) + g/2 - 1], s0={s0:g}', Gf(s0), asy, 0.5 / s0)
cS = np.log(4 * np.sqrt(2)) + EG / 2 - 1
printed('eq:G-limits', 'c_S = ln(4 sqrt2) + g/2 - 1 = 1.0215', cS, 1.0215, 4)
report('eq:G-limits', '2 sqrt(2pi) Q(r) - ln r -> c_S at r = 1e3', 2 * np.sqrt(2 * np.pi) * Qf(1e3) - np.log(1e3), cS, 1e-4)

# =====================================================================
print("\n== 12  the kinetic equation of an isotropic distribution")
# =====================================================================
for (pv, hwq) in [(1.0, 0.2), (0.4, 0.05), (3.0, 1.0)]:
    ppv = np.sqrt(pv * pv + 2 * hwq)
    gR_up = pv * (1 / pv) * np.log((pv + ppv) / abs(pv - ppv))         # g ~ p, R_+ ~ L/p
    gR_dn = ppv * (1 / ppv) * np.log((ppv + pv) / abs(ppv - pv))       # at E + hw going down to p
    report('eq:db-rates', f'g(E) R_+(E) = g(E+hw) R_-(E+hw), p={pv}', gR_up, gR_dn, 1e-14)
    report('eq:Rpm', f'R_+ = n_i v sigma_bar_+ (Eq. sigma-iso), p={pv}', (pv / 1.0) * 2 * np.pi / 3 * np.log((pv + ppv) / abs(pv - ppv)) / pv**2,
           2 * np.pi / 3 * np.log((pv + ppv) / abs(pv - ppv)) / pv, 1e-14)
# Fokker-Planck limit: exact averaged master equation vs the Langdon operator (m = 1, n_i Z^2 e^4 = 1, a0 = v0/w)
Fv = lambda v: np.exp(-(v / 1.3)**4)                  # a smooth flat-topped test distribution
dF = lambda v: -4 * v**3 / 1.3**4 * Fv(v)
def langdon_rhs(v, hwq, v0):
    Lv = np.log(2 * v * v / hwq)
    def flux(vv):
        return np.log(2 * vv * vv / hwq) / vv * dF(vv)
    dflux = rich(lambda h: flux(v + h), 1e-4 * v)
    return 2 * np.pi * v0**2 / 3 / v**2 * dflux


def exact_rhs(v, hwq, v0):
    E = v * v / 2
    a0q = v0 / hwq            # hbar = 1, w = hw
    pv = v
    Rp = 2 * np.pi / 3 * a0q**2 / pv * np.log((pv + np.sqrt(pv * pv + 2 * hwq)) / abs(pv - np.sqrt(pv * pv + 2 * hwq)))
    Rm = 0.0
    if E > hwq:
        pm = np.sqrt(pv * pv - 2 * hwq)
        Rm = 2 * np.pi / 3 * a0q**2 / pv * np.log((pv + pm) / abs(pv - pm))
    Fe = lambda EE: Fv(np.sqrt(2 * EE))
    return Rp * (Fe(E + hwq) - Fe(E)) + Rm * (Fe(E - hwq) - Fe(E))


for v in [0.6, 1.1]:
    e1 = abs(exact_rhs(v, 1e-3, 0.05) / langdon_rhs(v, 1e-3, 0.05) - 1)
    e2 = abs(exact_rhs(v, 1e-4, 0.05) / langdon_rhs(v, 1e-4, 0.05) - 1)
    report('eq:langdon', f'exact averaged master equation / Langdon operator -> 1 at hw = 1e-4, v={v}', 1 + e2, 1.0, 2e-3)
    check('eq:langdon', f'the difference falls with hw (v={v})', e2 < e1, f'{e1:.1e} -> {e2:.1e}')
# heating rate of an arbitrary isotropic F with fixed L: P = 4 pi m C L F(0)
Cq, Lq, mq = 0.37, 5.0, 1.0
Ft = lambda v: 2.0 * np.exp(-(v / 1.1)**3)
dFt = lambda v: -3 * v * v / 1.1**3 * Ft(v)
op = lambda v: Cq * Lq / v**2 * rich(lambda h: (1 / (v + h)) * dFt(v + h), 1e-5 * v)
Pnum = quad(lambda v: 4 * np.pi * v * v * (mq * v * v / 2) * op(v), 1e-4, 12, limit=400)[0]
report('eq:F0', 'int 4 pi v^2 (m v^2/2) (Langdon operator, fixed L) = 4 pi m C L F(0)', Pnum, 4 * np.pi * mq * Cq * Lq * Ft(0.0), 1e-6)
# Maxwell with the velocity-dependent logarithm: average of ln(4y/u) under y ~ Exp(1) = ln(4/u) - gamma
u = 0.01
avg = quad(lambda y: np.log(4 * y / u) * np.exp(-y), 0, np.inf)[0]
report('eq:F0', 'Maxwell average of L(v) = ln(4/u) - gamma (= Lambda to O(u))', avg, np.log(4 / u) - EG, 1e-10)


def superG(Mq):
    return (2 * np.pi)**1.5 * Mq / (4 * np.pi * Gamma(3 / Mq)) * (Gamma(5 / Mq) / (3 * Gamma(3 / Mq)))**1.5


for Mq, shown in [(2, 1.0), (3, 0.621), (4, 0.501), (5, 0.446)]:
    printed('eq:superG', f'P_M/P_2 for M={Mq}', superG(Mq), shown, 3)
printed('eq:superG', 'P_M/P_2 for M -> infinity', superG(2000.0), 0.336, 3)
# F(0) of the super-Gaussian by quadrature
Mq = 4
vm = np.sqrt(3 * Gamma(3 / Mq) / Gamma(5 / Mq))
Fm = lambda v: np.exp(-(v / vm)**Mq)
nrm = quad(lambda v: 4 * np.pi * v * v * Fm(v), 0, 20)[0]
ke = quad(lambda v: 4 * np.pi * v**4 * Fm(v), 0, 20)[0] / nrm
report('eq:superG', 'super-Gaussian v_m gives <v^2> = 3 v_t^2', ke, 3.0, 1e-10)
report('eq:superG', 'F(0)/n_e of the super-Gaussian = M/(4 pi v_m^3 Gamma(3/M))', 1 / nrm, Mq / (4 * np.pi * vm**3 * Gamma(3 / Mq)), 1e-10)

# =====================================================================
print("\n== A  reduction to Atzeni and Meyer-ter-Vehn")
# =====================================================================
T = 1.0
tau = 3 * np.sqrt(M_E) * (T * KEV)**1.5 / (4 * np.sqrt(2 * np.pi) * 1e21 * E_CH**4)
printed('tab:atzeni', '(10.133) coefficient 1.088e-11 s', tau * 1e11, 1.088, 3)
lam = 1e-4
nc1 = np.pi * M_E * C_LIGHT**2 / (E_CH**2 * lam**2)
printed('tab:atzeni', '(11.30) n_c = 1.115e21/lambda^2', nc1 / 1e21, 1.115, 3)
nuec = 1 / (3 * np.sqrt(M_E) * (T * KEV)**1.5 / (4 * np.sqrt(2 * np.pi) * nc1 * E_CH**4 * 8.0))
printed('tab:atzeni', '(11.33) nu_ec/c = 27.3 (printed 25)', nuec / C_LIGHT, 27.3, 1)
kT = KEV
lD = np.sqrt(kT / (4 * np.pi * E_CH**2 * 1e21))
v3 = np.sqrt(3 * kT / M_E)
ldB = HBAR / (2 * M_E * v3)
printed('app:atzeni-log', 'lambda_D(1 keV, 1e21) = 7.43e-7 cm', lD * 1e7, 7.43, 2)
printed('app:atzeni-log', 'lambda_dB = 2.52e-10 cm', ldB * 1e10, 2.52, 2)
printed('app:atzeni-log', 'ln(lambda_D/lambda_dB) = 7.99', np.log(lD / ldB), 7.99, 2)
printed('app:atzeni-log', 'printed constant 7.1 is 0.89 low', np.log(lD / ldB) - 7.1, 0.89, 2)
lam = 0.5e-4
nc5 = np.pi * M_E * C_LIGHT**2 / (E_CH**2 * lam**2)
printed('tab:atzeni', 'n_c(0.5 um) = 4.46e21', nc5 / 1e21, 4.46, 2)
printed('tab:atzeni', '(10.136) at 1 keV, n_c(0.5 um) = 6.35', 7.1 - 0.5 * np.log(nc5 / 1e21), 6.35, 2)
lD5 = np.sqrt(kT / (4 * np.pi * E_CH**2 * nc5))
printed('tab:atzeni', 'stated cutoffs at 1 keV, n_c(0.5 um) = 7.24', np.log(lD5 / ldB), 7.24, 2)
u05 = H_PL * C_LIGHT / lam / kT
printed('tab:atzeni', 'u = hw/kT at 0.5 um, 1 keV = 0.00248', u05, 0.00248, 5)
printed('tab:atzeni', 'quantum Lambda(u) = 6.81 there', Lam(u05), 6.81, 2)
printed('tab:atzeni', 'Born <g> = 3.758 there', gavg(u05), 3.758, 3)
# Kramers opacity (10.92) and power (10.93)
A_, rho_, hnu_ = 1.0, 1.0, 1.0
nu_ = hnu_ * KEV / H_PL
kap = 4 / (3 * np.sqrt(6 * np.pi)) * ALPHA**5 * A_B**2 * rho_ * C_LIGHT**3 / ((A_ * M_P)**2 * nu_**3) * np.sqrt(M_E * C_LIGHT**2 / (T * KEV))
printed('tab:atzeni', '(10.92) coefficient 2.74 (printed 2.78)', kap, 2.74, 2)
Pbr = 32 * np.pi / (3 * np.sqrt(6 * np.pi)) * E_CH**6 / (M_E * HBAR * C_LIGHT**2) * np.sqrt(T * KEV / (M_E * C_LIGHT**2)) / M_P**2 * 1e-7
printed('tab:atzeni', '(10.93) coefficient 1.74e17 W/g (printed 1.76e17)', Pbr / 1e17, 1.74, 2)
printed('tab:atzeni', '(10.92) and (10.93) both 1.4 percent below print', (2.78 / kap - 1) * 100, 1.4, 1)
printed('tab:atzeni', '(10.93): 1.4 percent', (1.76e17 / Pbr - 1) * 100, 1.4, 1)
# (10.91) = Maxwell average of (10.89), symbolically (CGS at T = 2 keV, h nu = 1 keV)
kT2, hnu2 = 2 * KEV, 1 * KEV
nu2 = hnu2 / H_PL
ne2, A2, Z2 = 1e22, 2.5, 1.0
vmin = np.sqrt(2 * hnu2 / M_E)
fv = lambda v: ne2 * 4 * np.pi * v * v * (M_E / (2 * np.pi * kT2))**1.5 * np.exp(-M_E * v * v / (2 * kT2))
dsdnu = lambda v: 32 * np.pi**2 * Z2**2 * E_CH**6 / (3 * np.sqrt(3) * M_E**2 * C_LIGHT**3 * v * v * hnu2)
eta_num = hnu2 / (4 * np.pi * A2 * M_P) * quad(lambda v: fv(v) * v * dsdnu(v), vmin, 40 * np.sqrt(kT2 / M_E), limit=400)[0]
eta_1091 = 16 * np.pi / (3 * np.sqrt(6 * np.pi)) * E_CH**6 * Z2**2 * ne2 / (M_E**2 * C_LIGHT**3 * np.sqrt(kT2 / M_E) * A2 * M_P) * np.exp(-hnu2 / kT2)
report('tab:atzeni', '(10.91) symbolic = Maxwell average of Kramers (10.89)', eta_num, eta_1091, 1e-7)
kap_kirch = eta_1091 * C_LIGHT**2 * np.exp(hnu2 / kT2) / (2 * H_PL * nu2**3)
rho2 = A2 * M_P * ne2 / Z2
kap_1092 = 4 / (3 * np.sqrt(6 * np.pi)) * ALPHA**5 * A_B**2 * Z2**3 * rho2 * C_LIGHT**3 / ((A2 * M_P)**2 * nu2**3) * np.sqrt(M_E * C_LIGHT**2 / kT2)
report('tab:atzeni', '(10.92) symbolic = Kirchhoff of (10.91) (Wien B for true absorption)', kap_kirch, kap_1092, 1e-12)
# (11.31) vs (11.34): numbers at lambda = 0.5 um, T = 1 keV, n_e = n_c, Z = 1
for Zi, Ai in [(1.0, 1.0), (3.5, 6.5)]:
    ne = nc5
    rho = Ai * M_P * ne / Zi
    nu5 = C_LIGHT / lam
    kap5 = 4 / (3 * np.sqrt(6 * np.pi)) * ALPHA**5 * A_B**2 * Zi**3 * rho * C_LIGHT**3 / ((Ai * M_P)**2 * nu5**3) * np.sqrt(M_E * C_LIGHT**2 / kT)
    K1134 = rho * kap5 * (1 - np.exp(-u05))
    nuT5 = 4 * np.sqrt(2 * np.pi) / 3 * (ne / Zi) * Zi**2 * E_CH**4 / (np.sqrt(M_E) * kT**1.5)
    K1131 = nuT5 * 8.0 / C_LIGHT
    Kours = nuT5 * (np.pi / np.sqrt(3)) * (1 - np.exp(-u05)) / u05 / C_LIGHT
    report('app:atzeni-ff', f'(11.34) with g=1 = (nu_T/c)(pi/sqrt3)(1-e^-u)/u exactly, Z={Zi}', K1134, Kours, 1e-12)
    printed('app:atzeni-ff', f'ratio (11.31)/(11.34) = 4.42, Z={Zi}', K1131 / K1134, 4.42, 2)
    printed('app:atzeni-ff', f'(11.34) per Z/(lambda^2 T^1.5) = 6.19, Z={Zi}', K1134 * 0.25 / Zi, 6.19, 2)
printed('app:atzeni-ff', 'g needed for agreement = (sqrt3/pi) 8 = 4.4', np.sqrt(3) / np.pi * 8, 4.4, 1)
# Table 1 (tab:scales)
I15 = 1e15 * 1e7
rows = [(0.351, 0.1, 0.0353, 0.52, 0.68, 4.15, 4.05, 4.39), (0.351, 1, 0.0035, 0.16, 0.21, 6.45, 6.36, 6.69),
        (0.351, 3, 0.0012, 0.10, 0.12, 7.55, 7.45, 7.79), (0.248, 0.1, 0.0500, 0.52, 0.48, 3.81, 3.71, 4.04),
        (0.248, 1, 0.0050, 0.16, 0.15, 6.11, 6.01, 6.34), (0.248, 3, 0.0017, 0.10, 0.09, 7.21, 7.11, 7.44)]
for (lum, TkeV, us, eta_s, b_s, L_s, Lr_s, La_s) in rows:
    lam = lum * 1e-4
    wq = 2 * np.pi * C_LIGHT / lam
    kT = TkeV * KEV
    u = HBAR * wq / kT
    nc = M_E * wq**2 / (4 * np.pi * E_CH**2)
    ne = nc / 4
    vt = np.sqrt(kT / M_E)
    E0 = np.sqrt(8 * np.pi * I15 / C_LIGHT)
    v0 = E_CH * E0 / (M_E * wq)
    tag = f'{lum} um, {TkeV} keV'
    printed('tab:scales', f'u at {tag}', u, us, 4)
    printed('tab:scales', f'eta/Z at {tag}', E_CH**2 / (HBAR * vt), eta_s, 2)
    printed('tab:scales', f'v0/vt at {tag}', v0 / vt, b_s, 2)
    printed('tab:scales', f'Lambda(u) at {tag}', Lam(u), L_s, 2)
    printed('tab:scales', f'Lambda_RPA(r=0.5) at {tag}', Lam_rpa(0.5, u), Lr_s, 2)
    printed('tab:scales', f'Atzeni (10.136) at {tag}', 7.1 - 0.5 * np.log(ne / 1e21) + np.log(TkeV), La_s, 2)
lam = 0.351e-4
printed('sec:dipole', 'v0/c = 0.0095 at 1e15 W/cm^2, 0.351 um', E_CH * np.sqrt(8 * np.pi * I15 / C_LIGHT) / (M_E * 2 * np.pi * C_LIGHT / lam) / C_LIGHT, 0.0095, 4)
printed('sec:dipole', 'v_t/c = 0.044 at 1 keV', np.sqrt(KEV / (M_E * C_LIGHT**2)), 0.044, 3)
nc351 = M_E * (2 * np.pi * C_LIGHT / lam)**2 / (4 * np.pi * E_CH**2)
lD = np.sqrt(KEV / (4 * np.pi * E_CH**2 * nc351 / 4))
printed('sec:closed', 'stated cutoffs of (10.136) at 1 keV, n_c(0.351)/4 = 7.58', np.log(lD / ldB), 7.58, 2)
printed('sec:closed', 'eta at thermal velocity = Z (2 Ry/kT)^1/2, 2 Ry = 27.2 eV', 2 * 13.605693122994, 27.2, 1)

# =====================================================================
print("\n== B  mathematical identities")
# =====================================================================
S = 50.0
val = quad(lambda s: erf(s) / s, 0, S, limit=400)[0]
report('eq:erf-integral', 'int_0^S erf(s)/s = ln S + g/2 + ln 2 + o(1), S = 50', val, np.log(S) + EG / 2 + np.log(2), 1e-12)
report('eq:erf-integral', 'int e^-s^2 ln s = -(sqrt(pi)/4)(g + 2 ln 2)', quad(lambda s: np.exp(-s * s) * np.log(s), 0, np.inf)[0],
       -np.sqrt(np.pi) / 4 * (EG + 2 * np.log(2)), 1e-10)
report('eq:ln2-integral', 'int_0^1 [(1-x^2)^-1/2 - 1] dx/x = ln 2', quad(lambda x: (1 / np.sqrt(1 - x * x) - 1) / x, 0, 1, limit=200)[0], np.log(2), 1e-9)
for a in [0.3, 5.0]:
    num = quad(lambda th: np.exp(-a * np.cos(th)**2), 0, 2 * np.pi)[0] / (2 * np.pi)
    report('app:math', f'<e^(-a cos^2 th)> = e^(-a/2) I0(a/2), a={a}', num, ive(0, a / 2), 1e-12)
    numc = quad(lambda th: np.cos(th)**2 * np.exp(-a * np.cos(th)**2), 0, 2 * np.pi)[0] / (2 * np.pi)
    report('app:math', f'<cos^2 e^(-a cos^2)> = (1/2) e^(-a/2)[I0 - I1](a/2), a={a}', numc, 0.5 * (ive(0, a / 2) - ive(1, a / 2)), 1e-12)
for s0 in [100.0, 1000.0]:
    report('app:math', f'<g(s0 cos)> -> (2/(sqrt(pi) s0))[ln(4 s0) + g/2 - 1], s0={s0:g}', avg_g(s0),
           2 / (np.sqrt(np.pi) * s0) * (np.log(4 * s0) + EG / 2 - 1), 0.3 / s0)

# =====================================================================
# summary and Appendix C table
# =====================================================================
ntot = sum(len(v) for v in RESULTS.values())
print(f"\n{ntot} checks, {len(ORDER)} labels, {len(FAILS)} failures, {time.time() - T0:.1f} s")


def ref_of(label):
    """Plain (link-free) references: hyperref's \\ref* makes no link annotation.  Link annotations
    inside table cells crash Preview on macOS 26.5 (its Vision table detector traps in
    CoreGraphics PageLayout::GetBoundsForRangeWithinLine while inserting the links into the
    structure tree); see the CHANGELOG Known-bad list."""
    if label.startswith('eq:'):
        return f'Eq.~(\\ref*{{{label}}})'
    if label.startswith('tab:'):
        return f'Table~\\ref*{{{label}}}'
    if label.startswith('fig:'):
        return f'Fig.~\\ref*{{{label}}}'
    return f'Sec.~\\ref*{{{label}}}'


nprinted = sum(1 for v in RESULTS.values() for r_ in v if r_[3] == 'printed')
nlogical = sum(1 for v in RESULTS.values() for r_ in v if r_[3] == 'logical')
lines = [r'\begin{longtable}{L{4.2cm}ccc}',
         r'\caption{The numerical checks of \texttt{detailed\_derivation\_check.py}, by the equation, table or section they verify. '
         r'The second column counts all checks of the row. The third gives the largest deviation among its quantitative checks, '
         r'relative, or absolute for quantities near zero, and is empty when the row has none. The fourth counts the numbers '
         r'printed in these notes that the script reproduces to their last printed digit, and the fifth the logical conditions, '
         r'such as a sign or a monotone decrease.}\label{tab:checks}\\',
         r'\toprule', r'Equation & Checks & Largest deviation & Printed & Logical\\', r'\midrule', r'\endfirsthead',
         r'\toprule', r'Equation & Checks & Largest deviation & Printed & Logical\\', r'\midrule', r'\endhead',
         r'\bottomrule', r'\endfoot']
lines[0] = r'\begin{longtable}{L{4.0cm}cccc}'
for lab in ORDER:
    rows_ = RESULTS[lab]
    errs = [r_[1] for r_ in rows_ if r_[3] == 'quant' and r_[1] is not None]
    npr = sum(1 for r_ in rows_ if r_[3] == 'printed')
    nlg = sum(1 for r_ in rows_ if r_[3] == 'logical')
    dev = ''
    if errs:
        mval = max(errs)
        if not np.isfinite(mval):
            dev = 'not finite'
        elif mval == 0:
            dev = '$0$'
        else:
            mant, ex = f'{mval:.1e}'.split('e')
            dev = f'${mant}\\times10^{{{int(ex)}}}$'
    lines.append(f'{ref_of(lab)} & {len(rows_)} & {dev} & {npr if npr else ""} & {nlg if nlg else ""}\\\\')
lines.append(f'\\midrule\nTotal & {ntot} & & {nprinted} & {nlogical}\\\\')
lines.append(r'\end{longtable}')
with open(os.path.join(HERE, 'detailed_derivation_checks.tex'), 'w') as fh:
    fh.write('% written by detailed_derivation_check.py -- do not edit by hand\n')
    fh.write('\n'.join(lines) + '\n')
if FAILS:
    print("FAILED:")
    for f_ in FAILS:
        print("   ", f_)
    sys.exit(1)
print("ALL CHECKS PASS")
