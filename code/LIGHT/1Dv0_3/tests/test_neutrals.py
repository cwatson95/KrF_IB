"""1Dv0_3's neutral gas (neutrals.py): its transport coefficients against
their sources and a direct quadrature, and its flow against exact
solutions. NumPy and SciPy only, so these run under any Python."""

import importlib
import math

import numpy as np
import pytest

G = importlib.import_module("1Dv0_3.neutrals")
M = importlib.import_module("1Dv0_3.mesh")
KB = G.KB
M_KR = G.mass_kg((("Kr", 1),))
M_F2 = G.mass_kg((("F", 2),))
M_F = G.mass_kg((("F", 1),))
M_E = 9.1093837139e-31


def approx(x, rel):
    return pytest.approx(x, rel=rel, abs=0)


# ---------------------------------------------------------------------------
# The coefficients
# ---------------------------------------------------------------------------
def _lj_quadrature(T_stars):
    """Omega(1,1)* and Omega(2,2)* of the (12-6) potential by quadrature of
    the Chapman-Enskog definitions (reduced units): the deflection angle by
    Gauss-Legendre in t (u = 1 - t^2, u = r_m/r), the cross sections by the
    trapezoid in b, the thermal average by the trapezoid in E/T*."""
    tn, tw = np.polynomial.legendre.leggauss(120)
    tn, tw = 0.5 * (tn + 1.0), 0.5 * tw
    B = np.concatenate([np.linspace(0.0, 2.5, 1501)[1:],
                        np.linspace(2.5, 7.0, 451)[1:]])
    E_grid = np.geomspace(2e-3, 2e3, 150)

    def r_outer(b, E):
        def g(r):
            return r * r - b * b - (4.0 / E) * (r ** -10 - r ** -4)
        r = np.maximum(b, 1.0) + 3.0
        hi = r.copy()
        active = np.ones_like(b, dtype=bool)
        for _ in range(3000):
            neg = (g(r) < 0.0) & active
            active &= ~neg
            if not active.any():
                break
            hi = np.where(active, r, hi)
            r = np.where(active, r - 0.005, r)
        lo = r
        for _ in range(55):
            mid = 0.5 * (lo + hi)
            below = g(mid) < 0.0
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        return 0.5 * (lo + hi)

    q1, q2 = [], []
    for E in E_grid:
        rm = r_outer(B, E)[:, None]
        beta = B[:, None] / rm
        u = 1.0 - tn[None, :] ** 2
        r = rm / u
        F = np.maximum(1.0 - beta ** 2 * u ** 2
                       - 4.0 * (r ** -12 - r ** -6) / E, 1e-300)
        chi = math.pi - 2.0 * beta[:, 0] * np.sum(
            tw * 2.0 * tn / np.sqrt(F), axis=1)
        bb = np.concatenate([[0.0], B])
        q1.append(2.0 * np.trapezoid(
            np.concatenate([[0.0], (1.0 - np.cos(chi)) * B]), bb))
        q2.append(3.0 * np.trapezoid(
            np.concatenate([[0.0], (1.0 - np.cos(chi) ** 2) * B]), bb))
    q1, q2 = np.array(q1), np.array(q2)
    out = []
    x = np.geomspace(1e-5, 70.0, 4000)
    for T in T_stars:
        res = []
        for q, s in ((q1, 1), (q2, 2)):
            qx = np.exp(np.interp(np.log(x * T), np.log(E_grid), np.log(q)))
            res.append(np.trapezoid(np.exp(-x) * x ** (s + 1) * qx, x)
                       / math.factorial(s + 1))
        out.append(res)
    return np.array(out)


def test_collision_integral_fits_reproduce_a_direct_quadrature():
    """The Neufeld, Janzen & Aziz (1972) fits against a quadrature of the
    (12-6) collision integrals from their definitions, over the T* the gap
    reaches (Kr at 300 to 3000 K is T* 1.7 to 17, F2-Kr from 2.1)."""
    T_stars = np.array([1.0, 1.7, 2.1, 5.0, 11.0, 20.0])
    quad = _lj_quadrature(T_stars)
    np.testing.assert_allclose(G.omega11(T_stars), quad[:, 0], rtol=3e-3)
    np.testing.assert_allclose(G.omega22(T_stars), quad[:, 1], rtol=3e-3)


# Svehla (1962) Table III, eta x 1e6 poise, which used HCB's tabulated
# Omega(2,2)*: Kr (3.655, 178.9) p.80 and F2 (3.357, 112.6) p.70
SVEHLA_ETA = {"Kr": [(300, 252.8), (1000, 641.0), (2000, 1015.4),
                     (3000, 1321.6)],
              "F2": [(300, 236.2), (1000, 554.0), (2000, 869.0)]}


def test_viscosities_follow_svehlas_tables():
    """Svehla's Kr and F2 viscosities from the same force constants, to
    0.9 % (his tables follow HCB's 1954 Omega(2,2)*, which the fit and the
    quadrature both put 0.1 to 0.9 % lower), and Kr's conductivity is
    (15/4)(k/m) eta, Svehla's 22.5e-6 cal/(cm s K) at 300 K."""
    for name, rows in SVEHLA_ETA.items():
        s, e, _ = G.LJ[name]
        m = M_KR if name == "Kr" else M_F2
        for T, ref in rows:
            eta = G.viscosity(np.array([T], float), m, s, e)[0]
            assert eta * 10.0 * 1e6 == approx(ref, 9e-3), (name, T)
    tr = G.BathTransport(("Kr",), np.array([M_KR]))
    kappa = tr.conductivity(np.array([300.0]))[0]
    assert kappa / 418.4 * 1e6 == approx(22.5, 1e-2)


def test_diffusion_follows_brokaws_examples_and_their_measurements():
    """Brokaw (1961) eqs. (11) and (22)-(24): argon's self-diffusion at
    22 C and 1 atm (his chart 0.182 cm^2/s, measured 0.180, p.11) and H2-N2
    at 85 C (his chart 1.058, measured 1.052, p.22), from his Table I
    constants."""
    amu = G.AMU
    p = 101325.0

    def d_cm2(T, mA, mB, ljA, ljB):
        n = p / (KB * T)
        return G.n_diffusion(np.array([T]), mA * amu, mB * amu, ljA,
                             ljB)[0] / n * 1e4
    ar = d_cm2(295.2, 39.94, 39.94, (3.421, 119.5), (3.421, 119.5))
    h2n2 = d_cm2(358.2, 2.016, 28.02, (2.915, 38.0), (3.681, 91.5))
    assert ar == approx(0.180, 5e-3) and ar == approx(0.182, 2e-2)
    assert h2n2 == approx(1.052, 1.5e-2) and h2n2 == approx(1.058, 2e-2)


def test_f2_in_kr():
    """F2 and F diffuse in the Kr bath, each pair's sigma the mean and
    epsilon the geometric mean (Brokaw eqs. (23), (24)). F2-Kr at 300 K and
    1 atm is 0.1456 cm^2/s by hand from his eq. (22) (sigma_12 3.506 A,
    epsilon_12/k 141.9 K, T* 2.114, Omega(1,1)* 1.0553), and F, lighter and
    smaller, diffuses faster."""
    names = ("Kr", "F2", "F")
    tr = G.BathTransport(names, np.array([M_KR, M_F2, M_F]))
    assert tr.diffusing == [1, 2]
    T = np.array([300.0])
    n = 101325.0 / (KB * 300.0)
    d = {i: v[0] / n * 1e4 for i, v in tr.n_diffusion(T).items()}
    assert d[1] == approx(0.1456, 5e-3)
    assert d[2] > d[1]                               # lighter, smaller


# ---------------------------------------------------------------------------
# The flow
# ---------------------------------------------------------------------------
def _kr_gas(mesh, transport=None, t_wall=None):
    return G.NeutralGas(mesh, ("Kr",), [M_KR], [1.5], [True], bath="Kr",
                        transport=transport, t_wall=t_wall)


def _mix(mesh, transport="bath", t_wall=300.0):
    names = ("e-", "Kr", "F2", "F")
    masses = np.array([M_E, M_KR, M_F2, M_F])
    tr = G.BathTransport(names, masses) if transport == "bath" else transport
    return G.NeutralGas(mesh, names, masses, [0.0, 1.5, 2.5, 1.5],
                        [False, True, True, True], e_col=0, transport=tr,
                        t_wall=t_wall)


def _run(gas, dens, T, u, t_end, eps=None, dt_max=None, record=None):
    t = 0.0
    while t < t_end * (1 - 1e-12):
        dt = min(gas.stable_dt(dens, T, u), t_end - t)
        if dt_max:
            dt = min(dt, dt_max)
        r = gas.step(dens, T, u, eps, dt)
        dens, T, u, eps = r.dens, r.T, r.u, r.eps
        t += dt
        if record is not None:
            record(t, dens, T, u)
    return dens, T, u, eps


def _sod_exact(x, t, x0, rL, pL, rR, pR, g):
    """Toro's exact solution of Sod's problem (a rarefaction left, a shock
    right, both states at rest): the density at x, and u* (the contact's
    speed)."""
    cL, cR = math.sqrt(g * pL / rL), math.sqrt(g * pR / rR)

    def f(p, r, pk, c):
        if p > pk:
            A, B = 2 / ((g + 1) * r), (g - 1) / (g + 1) * pk
            return (p - pk) * math.sqrt(A / (p + B))
        return 2 * c / (g - 1) * ((p / pk) ** ((g - 1) / (2 * g)) - 1)
    lo, hi = pR, pL
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = ((lo, mid) if f(mid, rL, pL, cL) + f(mid, rR, pR, cR) > 0
                  else (mid, hi))
    ps = 0.5 * (lo + hi)
    us = 0.5 * (f(ps, rR, pR, cR) - f(ps, rL, pL, cL))
    k = (g - 1) / (g + 1)
    shock = cR * math.sqrt((g + 1) / (2 * g) * ps / pR + (g - 1) / (2 * g))
    c_star = cL * (ps / pL) ** ((g - 1) / (2 * g))
    rho = np.empty_like(x)
    for i, xx in enumerate(x):
        s = (xx - x0) / t
        if s < -cL:
            rho[i] = rL
        elif s < us - c_star:
            c = 2 / (g + 1) * (cL - (g - 1) / 2 * s)
            rho[i] = rL * (c / cL) ** (2 / (g - 1))
        elif s < us:
            rho[i] = rL * (ps / pL) ** (1 / g)
        elif s < shock:
            rho[i] = rR * (ps / pR + k) / (k * ps / pR + 1)
        else:
            rho[i] = rR
    return rho, us


@pytest.mark.parametrize("kind", ["uniform", "graded"])
def test_sod_shock_tube(kind):
    """Sod's problem in krypton (densities 8:1, pressures 10:1) against the
    exact solution: the density's L1 error is below 0.6 % on 200 uniform
    cells (1 % on the graded mesh, whose middle cells, where the waves
    are, are twice as wide) and falls as the cells shrink, and the mass is
    kept exactly."""
    L = 1e-2
    errs = []
    for N in (200, 400):
        mesh = (M.uniform(N, L) if kind == "uniform"
                else M.graded(N, L, 5e-6 * 200 / N))
        gas = _kr_gas(mesh)
        x = mesh.centres
        nL, TL, nR, TR = 2.45e25, 300.0, 2.45e25 / 8, 240.0
        dens = np.where(x < L / 2, nL, nR)[:, None].astype(float)
        T = np.where(x < L / 2, TL, TR).astype(float)
        cL = math.sqrt(5 / 3 * KB * TL / M_KR)
        t_end = 0.2 * L / cL
        mass0 = dens[:, 0] @ mesh.widths
        dens, T, u, _ = _run(gas, dens, T, np.zeros(N), t_end)
        rho = dens[:, 0] * M_KR
        exact, _ = _sod_exact(x, t_end, L / 2, nL * M_KR, nL * KB * TL,
                              nR * M_KR, nR * KB * TR, 5 / 3)
        errs.append(np.sum(np.abs(rho - exact) * mesh.widths)
                    / (L * nL * M_KR))
        assert dens[:, 0] @ mesh.widths == approx(mass0, 1e-14)
    limit = 6e-3 if kind == "uniform" else 1e-2
    assert errs[0] < limit and errs[1] < 0.75 * errs[0], errs


def test_a_species_front_rides_the_contact():
    """Sod's problem again with F atoms (half the particles) on the right,
    one gamma throughout: the F front sits where the exact contact is, u*
    times t from the diaphragm, to a third of a cell, the F stays between
    its two starting fractions, and none is lost."""
    L, N = 1e-2, 200
    mesh = M.uniform(N, L)
    gas = G.NeutralGas(mesh, ("Kr", "F"), [M_KR, M_F], [1.5, 1.5],
                       [True, True], bath="Kr")
    x = mesh.centres
    nL, TL = 2.45e25, 300.0
    m_right = 0.5 * (M_KR + M_F)
    nR = nL * M_KR / 8 / m_right                 # 8:1 in mass density
    pL = nL * KB * TL
    TR = pL / 10 / (nR * KB)                     # 10:1 in pressure
    left = x < L / 2
    dens = np.stack([np.where(left, nL, 0.5 * nR),
                     np.where(left, 0.0, 0.5 * nR)], axis=1)
    T = np.where(left, TL, TR)
    t_end = 0.2 * L / math.sqrt(5 / 3 * KB * TL / M_KR)
    f0 = dens[:, 1] @ mesh.widths
    dens, T, u, _ = _run(gas, dens, T, np.zeros(N), t_end)
    _, us = _sod_exact(x, t_end, L / 2, nL * M_KR, pL, nR * m_right,
                       pL / 10, 5 / 3)
    xf = dens[:, 1] / dens.sum(axis=1)
    i = int(np.argmax(xf > 0.25))
    front = x[i - 1] + (0.25 - xf[i - 1]) / (xf[i] - xf[i - 1]) * (
        x[i] - x[i - 1])
    assert abs(front - (L / 2 + us * t_end)) < 0.34 * L / N
    assert xf.min() >= 0.0 and xf.max() <= 0.5 * (1 + 1e-12)
    assert dens[:, 1] @ mesh.widths == approx(f0, 1e-14)


def test_a_contact_at_rest_stays_at_rest():
    """Equal pressures, a temperature jump of 4 and another composition on
    either side (F2 on the cold side, F atoms on the hot), no transport:
    HLLC holds a contact, so nothing moves and nothing mixes."""
    mesh = M.graded(64, 1e-4, 0.2e-6)
    gas = _mix(mesh, transport=None, t_wall=None)
    x = mesh.centres / mesh.length
    hot = x > 0.4
    n = np.where(hot, 2.45e25 / 4, 2.45e25)
    xf2 = np.where(hot, 0.0, 0.005)
    xf = np.where(hot, 0.01, 0.0)
    dens = np.stack([np.full(64, 1e18), n * (1 - xf2 - xf), n * xf2,
                     n * xf], axis=1)
    T = np.where(hot, 1200.0, 300.0)
    p = (dens[:, 1:] @ np.ones(3)) * KB * T
    d1, T1, u1, _ = _run(gas, dens, T, np.zeros(64), 2e-7)
    p1 = (d1[:, 1:] @ np.ones(3)) * KB * T1
    assert np.max(np.abs(u1)) < 1e-9
    assert np.max(np.abs(p1 / p - 1)) < 1e-12
    assert np.array_equal(d1[:, 2] > 0, dens[:, 2] > 0)


def test_sound_rings_at_c_and_decays_at_the_classical_rate():
    """A standing acoustic wave, cos(pi x/L) between adiabatic walls 20 um
    apart in Kr at 300 K: it rings at c pi/L, c = sqrt(5 k T/3 m), and its
    amplitude decays at Stokes and Kirchhoff's (k^2/2)((4/3) mu/rho +
    (gamma - 1) kappa/(rho c_p)), the viscosity and conductivity the gas
    carries."""
    L, N = 20e-6, 100
    mesh = M.uniform(N, L)
    tr = G.BathTransport(("Kr",), np.array([M_KR]))
    gas = _kr_gas(mesh, tr, t_wall=None)
    x = mesh.centres
    eps = 1e-4
    n0, T0 = 2.45e25, 300.0
    mode = np.cos(math.pi * x / L)
    dens = (n0 * (1 + eps * mode))[:, None]
    T = T0 * (1 + (2 / 3) * eps * mode)
    ts, amps = [], []

    def rec(t, d, TT, uu):
        p = d[:, 0] * KB * TT
        ts.append(t)
        amps.append(2.0 / L * np.sum((p - np.mean(p)) * mode * mesh.widths))
    _run(gas, dens, T, np.zeros(N), 1.2e-6, dt_max=2e-10, record=rec)
    ts, amps = np.array(ts), np.array(amps)
    from scipy.optimize import curve_fit
    c = math.sqrt(5 / 3 * KB * T0 / M_KR)
    k = math.pi / L
    rho = n0 * M_KR
    mu = tr.viscosity(np.array([T0]))[0]
    kap = tr.conductivity(np.array([T0]))[0]
    cp = 2.5 * KB / M_KR
    alpha = 0.5 * k * k * (4 / 3 * mu / rho + (2 / 3) * kap / (rho * cp))
    p0 = n0 * KB * T0 * 5 / 3 * eps

    def model(t, a, al, om, ph):
        return a * np.exp(-al * t) * np.cos(om * t + ph)
    fit, _ = curve_fit(model, ts, amps, p0=[p0, alpha, c * k, 0.0])
    assert fit[2] == approx(c * k, 2e-3)
    assert fit[1] == approx(alpha, 5e-2)


def test_heat_spreads_at_constant_pressure():
    """An isobaric temperature ripple, cos(pi x/L) at 1 % between adiabatic
    walls, decays at kappa k^2/(rho c_p), the gas expanding where it heats
    and keeping its pressure."""
    L, N = 50e-6, 100
    mesh = M.uniform(N, L)
    tr = G.BathTransport(("Kr",), np.array([M_KR]))
    gas = _kr_gas(mesh, tr, t_wall=None)
    x = mesh.centres
    mode = np.cos(math.pi * x / L)
    T = 300.0 * (1 + 0.01 * mode)
    dens = (2.45e25 * 300.0 / T)[:, None]
    ts, amps = [], []

    def rec(t, d, TT, uu):
        ts.append(t)
        amps.append(2.0 / L * np.sum((TT - np.mean(TT)) * mode * mesh.widths))
    _run(gas, dens, T, np.zeros(N), 2e-6, dt_max=2e-10, record=rec)
    ts, amps = np.array(ts), np.array(amps)
    late = ts > 0.5e-6
    rate = -np.polyfit(ts[late], np.log(amps[late]), 1)[0]
    rho_cp = 2.45e25 * 2.5 * KB
    kap = tr.conductivity(np.array([300.0]))[0]
    assert rate == approx(kap * (math.pi / L) ** 2 / rho_cp, 3e-2)


def _interdiffusing(N=100, L=50e-6):
    """Kr and F2 at rest, 300 K and one pressure, the F2 fraction stepping
    from 0 to 0.2 over 5 um (tanh) mid-gap, adiabatic walls."""
    names = ("Kr", "F2")
    m = np.array([M_KR, M_F2])
    mesh = M.uniform(N, L)
    gas = G.NeutralGas(mesh, names, m, [1.5, 2.5], [True, True], bath="Kr",
                       transport=G.BathTransport(names, m), t_wall=None)
    xf2 = 0.1 + 0.1 * np.tanh((mesh.centres - L / 2) / 5e-6)
    n0 = 2.45e25
    dens = np.stack([n0 * (1 - xf2), n0 * xf2], axis=1)
    return gas, dens, m


def test_diffusion_carries_no_mass():
    """One substep from rest at one pressure: F2 and Kr interdiffuse, so
    the particles move (the masses differ), but Kr's flux balances F2's
    mass flux and no cell's mass changes beyond round-off."""
    gas, dens, m = _interdiffusing()
    T, u = np.full(len(dens), 300.0), np.zeros(len(dens))
    d1 = gas._substep(dens, T, u, None, gas.stable_dt(dens, T, u))[0]
    assert np.max(np.abs(d1 @ m / (dens @ m) - 1)) < 1e-14
    assert np.max(np.abs(d1.sum(axis=1) / dens.sum(axis=1) - 1)) > 1e-6


def test_interdiffusion_releases_no_heat():
    """Ideal gases mix at one pressure and temperature without heat: over
    0.2 us of Kr-F2 interdiffusion the temperature moves by the sound the
    first instants make (0.06 K), where an enthalpy flux without the pV
    term, h_s = cv_s k T, heats and cools the front by 0.4 K."""
    gas, dens, m = _interdiffusing()
    T = np.full(len(dens), 300.0)
    _, T1, u1, _ = _run(gas, dens, T, np.zeros(len(dens)), 2e-7)
    assert np.max(np.abs(T1 - 300.0)) < 0.15
    assert np.max(np.abs(u1)) > 0.01                 # the mixture does move


def test_conduction_into_isothermal_electrodes():
    """A temperature ripple sin(2 pi x/L) at 1 % between electrodes held
    at 300 K, 20 cells: it decays at kappa k^2/(rho c_p) (its mean is zero,
    so the pressure holds), the electrode a half cell from the first
    centre. A full cell there slows it by 9 %."""
    L, N = 50e-6, 20
    mesh = M.uniform(N, L)
    tr = G.BathTransport(("Kr",), np.array([M_KR]))
    gas = _kr_gas(mesh, tr, t_wall=300.0)
    mode = np.sin(2 * math.pi * mesh.centres / L)
    T = 300.0 * (1 + 0.01 * mode)
    n0 = 2.45e25
    ts, amps = [], []

    def rec(t, d, TT, uu):
        ts.append(t)
        amps.append(2.0 / L * np.sum((TT - 300.0) * mode * mesh.widths))
    _run(gas, (n0 * 300.0 / T)[:, None], T, np.zeros(N), 2e-6, dt_max=2e-10,
         record=rec)
    ts, amps = np.array(ts), np.array(amps)
    late = ts > 0.5e-6
    rate = -np.polyfit(ts[late], np.log(amps[late]), 1)[0]
    kap = tr.conductivity(np.array([300.0]))[0]
    expect = kap * (2 * math.pi / L) ** 2 / (n0 * 2.5 * KB)
    assert rate == approx(expect, 3e-2)


def test_a_trace_species_spreads_as_2_d_t():
    """F2 at 1e-4 in a Gaussian 5 um wide, in Kr at rest: its variance grows
    by 2 D t, D the n D the gas carries over n, and its number is kept."""
    L, N = 200e-6, 400
    mesh = M.uniform(N, L)
    n0 = 2.45e25
    nd = 0.1456e-4 * n0
    gas = _mix(mesh, transport=G.ConstantTransport(n_d={2: nd}),
               t_wall=None)
    x = mesh.centres - L / 2
    s0 = 5e-6
    xf2 = 1e-4 * np.exp(-0.5 * (x / s0) ** 2)
    dens = np.stack([np.zeros(N), n0 * (1 - xf2), n0 * xf2, np.zeros(N)],
                    axis=1)
    t_end = 1e-6
    w = mesh.widths

    def var(d):
        f = d[:, 2] * w
        return np.sum(f * x * x) / np.sum(f)
    v0, tot0 = var(dens), dens[:, 2] @ w
    d1, *_ = _run(gas, dens, np.full(N, 300.0), np.zeros(N), t_end)
    assert var(d1) - v0 == approx(2 * nd / n0 * t_end, 1e-2)
    assert d1[:, 2] @ w == approx(tot0, 1e-13)


def test_a_closed_box_keeps_every_species_and_its_energy():
    """A random smooth state on a graded mesh, every transport on,
    adiabatic walls: each species' number, the energy and the electrons'
    energy are kept to round-off, and the electrons' mean energy stays
    within its starting range."""
    mesh = M.graded(48, 1e-4, 0.5e-6)
    gas = _mix(mesh, t_wall=None)
    x = mesh.centres / mesh.length
    dens = np.stack([1e20 * (1 + 0.5 * np.sin(3 * x)),
                     2.4e25 * (1 + 0.2 * np.cos(2 * math.pi * x)),
                     1e23 * (1 + 0.8 * np.sin(5 * x)),
                     5e22 * (1 + 0.5 * np.cos(7 * x))], axis=1)
    T = 300 + 900 * np.exp(-((x - 0.3) / 0.1) ** 2)
    u = 20 * np.sin(math.pi * x)
    eps = 1 + np.exp(-((x - 0.6) / 0.2) ** 2)
    h = mesh.widths
    tot0, e0 = dens.T @ h, gas.energy(dens, T, u) @ h
    ee0 = (dens[:, 0] * eps) @ h
    d1, T1, u1, eps1 = _run(gas, dens, T, u, 1e-8, eps=eps)
    np.testing.assert_allclose(d1.T @ h, tot0, rtol=1e-14)
    assert gas.energy(d1, T1, u1) @ h == approx(e0, 1e-14)
    assert (d1[:, 0] * eps1) @ h == approx(ee0, 1e-14)
    assert eps1.min() >= eps.min() * (1 - 1e-12)
    assert eps1.max() <= eps.max() * (1 + 1e-12)


def test_isothermal_walls_take_the_heat_the_books_say():
    """Hot gas between electrodes at 300 K: the energy the gas loses is the
    heat the step books into the two walls, to round-off, and a symmetric
    state stays symmetric."""
    mesh = M.graded(40, 1e-4, 0.5e-6)
    gas = _mix(mesh, t_wall=300.0)
    x = mesh.centres / mesh.length
    T = 300 + 600 * np.sin(math.pi * x)
    n = 2.45e25 * 300 / T
    dens = np.stack([np.full(40, 1e19), 0.995 * n, 0.005 * n,
                     np.zeros(40)], axis=1)
    h = mesh.widths
    e0 = gas.energy(dens, T, np.zeros(40)) @ h
    lost = 0.0
    u = np.zeros(40)
    eps = np.full(40, 1.0)
    for _ in range(300):
        r = gas.step(dens, T, u, eps, 20e-12)
        dens, T, u, eps = r.dens, r.T, r.u, r.eps
        lost += float(np.sum(r.wall_heat))
    de = gas.energy(dens, T, u) @ h - e0
    assert lost > 0.0
    assert -de == approx(lost, 1e-9)
    assert np.max(np.abs(T - T[::-1])) < 1e-9 * T.max()
    assert np.max(np.abs(u + u[::-1])) < 1e-9 * np.max(np.abs(u))


def test_the_substep_is_the_tighter_of_sound_and_diffusion():
    """stable_dt is 0.5 h/(|u| + c) on the narrowest cell, or 0.4 h^2/(2
    chi), chi the largest of (4/3) nu, kappa/(n c_v k) and D, where that is
    shorter: by hand, Kr at 1 atm and 300 K on 10 um cells (sound binds),
    and at 0.01 atm and 3000 K on 0.5 um cells (conduction binds)."""
    tr = G.BathTransport(("Kr",), np.array([M_KR]))
    for N, n0, T0, u0, binds in ((100, 2.45e25, 300.0, 5.0, "sound"),
                                 (200, 2.45e23, 3000.0, 0.0, "heat")):
        mesh = M.uniform(N, N * (10e-6 if binds == "sound" else 0.5e-6))
        gas = _kr_gas(mesh, tr, t_wall=None)
        h = mesh.widths[0]
        T = np.full(N, T0)
        c = math.sqrt(5 / 3 * KB * T0 / M_KR)
        mu = G.viscosity(np.array([T0]), M_KR, *G.LJ["Kr"][:2])[0]
        kap = 3.75 * KB / M_KR * mu
        chi = max(4 / 3 * mu / (n0 * M_KR), kap / (1.5 * n0 * KB))
        sound, heat = 0.5 * h / (abs(u0) + c), 0.4 * h * h / (2 * chi)
        dt = gas.stable_dt(np.full((N, 1), n0), T, np.full(N, u0))
        assert dt == approx(min(sound, heat), 1e-12)
        assert (sound < heat) == (binds == "sound")


def test_a_uniform_gas_stays_put_to_the_bit():
    """A uniform gas at rest, every transport on, walls at its temperature,
    on a graded mesh: two hundred steps leave it exactly as it was."""
    mesh = M.graded(64, 1e-2, 2e-6)
    gas = _mix(mesh, t_wall=300.0)
    dens = np.tile([1e20, 2.45e25 * 0.995, 2.45e25 * 0.005, 0.0], (64, 1))
    T, u, eps = np.full(64, 300.0), np.zeros(64), np.full(64, 2.0)
    d1, T1, u1, eps1 = dens, T, u, eps
    for _ in range(200):
        r = gas.step(d1, T1, u1, eps1, 50e-12)
        d1, T1, u1, eps1 = r.dens, r.T, r.u, r.eps
    assert np.array_equal(d1, dens) and np.array_equal(T1, T)
    assert np.array_equal(u1, u) and np.array_equal(eps1, eps)
