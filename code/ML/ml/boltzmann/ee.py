"""Electron-electron collisions and the beam-sustained EEDF, on top of the
vendored two-term solver (boltz_twin.py stays byte-identical).

LoKI-MC and the vendored solver describe a swarm: electrons too few to
collide with one another. In an e-beam-pumped laser the ionization degree
is about 1e-5, and Coulomb collisions between electrons then act on the EEDF
faster than F2 does. This module adds them.

The e-e term. The isotropic Landau operator for like particles, written in
energy for the EEDF normalised as int sqrt(eps) F0 d eps = 1, adds to the
energy-space flux G = W F0 - D dF0/d eps of boltz_twin
    W_ee = -3 alpha A1,   D_ee = 2 alpha (A2 + eps^(3/2) A3)
    A1 = int_0^eps u^(1/2) F0 du,  A2 = int_0^eps u^(3/2) F0 du,
    A3 = int_eps^inf F0 du
    alpha = (n_e/N) (e^2 / (24 pi eps0^2)) ln(Lambda) gamma
    ln(Lambda) = ln[12 pi (eps0 kTe/e^2)^(3/2) / n_e^(1/2)],  kTe = 2/3 <eps>
(Hagelaar & Pitchford 2005, PSST 14 722, sec. 2.5, re-derived here from
the Rosenbluth coefficients: a Maxwellian of any temperature makes the e-e
flux vanish, and a fast electron's friction and diffusion come out as the
Spitzer values; tests/test_electron.py checks both).
The coefficients depend on F0, so the solution is iterated to
self-consistency.

Two drives:
  solve_field   the swarm drive (a DC field), plus e-e at a given n_e/N.
                With n_e/N = 0 it is the vendored solver's answer.
  solve_beam    no field. Electrons appear at the rate of ion-pair
                production with the sub-excitation spectrum of the
                degradation cascade, are heated by superelastic collisions
                with Kr*, lose energy to elastic, vibrational and electronic
                collisions, and leave by attachment and by recombination
                (a loss frequency the same for every electron). The solution
                is the absolute EEDF, so n_e and <eps> come out of it.
"""
from __future__ import annotations

import math

import numpy as np

import boltz_twin as BT
import lx_lib as L

EPS0 = 8.8541878128e-12
C_EE = L.QE ** 2 / (24.0 * math.pi * EPS0 ** 2)     # V^2 m^2


def coulomb_log(ne_m3, mean_eV):
    kte = 2.0 / 3.0 * max(mean_eV, 1e-3) * L.QE     # J
    lam = 12.0 * math.pi * (EPS0 * kte / L.QE ** 2) ** 1.5 / math.sqrt(ne_m3)
    return math.log(max(lam, 2.0))


def ee_coeffs(F, ec, eb, du, alpha):
    """W_ee, D_ee at the N-1 inner boundaries from the cell values of F0.
    The powers of u are integrated exactly over each cell (F constant in a
    cell): the midpoint rule is 6 % off in the first cell, where sqrt(u)
    is not smooth."""
    edges = np.arange(len(ec) + 1) * du
    w1 = 2.0 / 3.0 * np.diff(edges ** 1.5)
    w2 = 2.0 / 5.0 * np.diff(edges ** 2.5)
    a1 = np.cumsum(w1 * F)[:-1]
    a2 = np.cumsum(w2 * F)[:-1]
    a3 = (np.sum(F) - np.cumsum(F))[:-1] * du
    return -3.0 * alpha * a1, 2.0 * alpha * (a2 + eb ** 1.5 * a3)


def flux_matrix(W, D, du, N):
    """The Scharfetter-Gummel flux part of boltz_twin.assemble."""
    z = W * du / D
    a = D / du * BT.bern(-z)
    b = -D / du * BT.bern(z)
    A = np.zeros((N, N))
    i = np.arange(N - 1)
    A[i, i] += a
    A[i, i + 1] += b
    A[i + 1, i] -= a
    A[i + 1, i + 1] -= b
    return A


def base_WD(EN, Tg, eb, sig_t, sig_eps_b):
    """W and D exactly as boltz_twin.assemble builds them."""
    G = L.GAMMA
    W = -G * eb ** 2 * sig_eps_b
    D = (G / 3 * (EN * BT.TD) ** 2 * eb / sig_t
         + G * (L.KB * Tg / L.QE) * eb ** 2 * sig_eps_b)
    return W, D


def _with_ee(mix, EN, Tg, N, emax, nubar, F, ne_N, ne_m3, ion):
    A, ec, eb, du, sig_t, sig_eps = BT.assemble(mix, EN, Tg, N, emax, nubar,
                                                ion)
    if ne_N > 0:
        W0, D0 = base_WD(EN, Tg, eb, sig_t, sig_eps)
        mean = float(np.sum(ec ** 1.5 * F) * du)
        alpha = ne_N * C_EE * coulomb_log(ne_m3, mean) * L.GAMMA
        We, De = ee_coeffs(F, ec, eb, du, alpha)
        A = A - flux_matrix(W0, D0, du, N) + flux_matrix(W0 + We, D0 + De,
                                                         du, N)
    return A, ec, eb, du


def ee_power(F, ec, eb, du, alpha):
    """int G_ee d eps [per N per n_e, eV m^3/s]: zero for the exact operator,
    so its size against the field power measures the discretization's
    energy error."""
    We, De = ee_coeffs(F, ec, eb, du, alpha)
    dF = (F[1:] - F[:-1]) / du
    Fb = 0.5 * (F[1:] + F[:-1])
    return float(np.sum(We * Fb - De * dF) * du)


def solve_field(mix, EN, Tg, N, emax, ne_N=0.0, ne_m3=1e20, ion="equal",
                growth=True, tol=1e-9, maxit=4000):
    """Two-term EEDF in a DC field with e-e collisions at ionization degree
    ne_N (ne_m3 sets the Coulomb logarithm). Returns boltz_twin's solution
    dict (usable by BT.swarm) plus 'ee_iters'."""
    sol = BT.solve_fixed(mix, EN, Tg, N, emax, ion, growth=growth)
    if ne_N <= 0:
        sol["ee_iters"] = 0
        return sol
    F, nubar = sol["F"], sol["nubar"]
    for it in range(1, maxit + 1):
        A, ec, eb, du = _with_ee(mix, EN, Tg, N, emax, nubar, F, ne_N, ne_m3,
                                 ion)
        A[0, :] = 0.0
        A[0, 0] = 1.0
        rhs = np.zeros(N)
        rhs[0] = 1.0
        Fn = np.linalg.solve(A, rhs)
        Fn /= float(np.sum(np.sqrt(ec) * Fn) * du)
        Fn = 0.5 * (Fn + F)
        Fn /= float(np.sum(np.sqrt(ec) * Fn) * du)
        change = float(np.max(np.abs(Fn - F)) / np.max(F))
        F = Fn
        if growth:
            nubar = 0.5 * nubar + 0.5 * BT.net_ion(mix, F, ec, du)
        if change < tol:
            break
    else:
        raise RuntimeError(f"e-e iteration: no convergence ({change:.1e})")
    sol.update(F=F, nubar=nubar, ee_iters=it)
    return sol


def injection(ec, e_top, mean):
    """A linear sub-excitation spectrum on [0, e_top] with the given mean
    [per eV, integrates to 1]. mean = e_top/2 is flat and mean = e_top/3 is
    a triangle falling to zero at e_top."""
    b = (e_top / 2.0 - mean) / (e_top / 6.0)
    if not -1.0 <= b <= 1.0:
        raise ValueError("mean outside [e_top/3, 2 e_top/3]")
    g = np.where(ec < e_top, 1.0 + b * (1.0 - 2.0 * ec / e_top), 0.0)
    return g / float(np.sum(g) * (ec[1] - ec[0]))


def solve_beam(mix, Tg, N_gas, N, emax, pairs_m3s, spectrum, nu_rec,
               ne_guess, ee_on=True, rec_te=None, tol=1e-9, maxit=4000,
               rec_pow=-0.5):
    """The steady EEDF of a beam-sustained plasma (no field).
    mix: [(gas, x)] with any superelastic populations already set on the
    gas's collisions (c['wsup'], c['gratio']); N_gas [m^-3]; pairs_m3s the
    ion-pair production rate; spectrum(ec) the injection spectrum [1/eV];
    nu_rec a loss frequency for every electron (recombination) [1/s].
    ee_on=False drops the e-e term (the swarm codes' physics).
    rec_te: None makes recombination the same for every electron. A
    temperature makes it dissociative recombination, sigma ~ 1/eps
    (k ~ Te^-1/2, LIGHT's R11p template): nu_rec is then that template's
    rate at rec_te, and the sink falls as eps^-1/2 across the EEDF.
    rec_pow: the template's exponent n in k ~ Te^n, sigma v ~ eps^n
    (-1/2 is v0's R11p, -0.55 is v0_1's Shiu & Biondi law).
    Cell i then holds sqrt(e_i) Phi_i du electrons per m^3 with
    Phi = n_e F0, and N_gas (A Phi)_i = Q_i du.
    Returns dict(F (normalised), ne, mean, ec, du, emax, N, iters)."""
    du = emax / N
    ec = (np.arange(N) + 0.5) * du
    q = pairs_m3s * spectrum(ec) * du
    te = 1.0
    F = np.exp(-ec / te)
    F /= float(np.sum(np.sqrt(ec) * F) * du)
    ne = ne_guess
    # sigma v = C eps^-1/2 with <sigma v>_Maxwellian(rec_te) = nu_rec:
    # the loss from cell i is C (sqrt(e_i) F_i du) / sqrt(e_i) = C F_i du.
    # For sigma v = C eps^n, <sigma v> = C Te^n Gamma(n + 3/2)/Gamma(3/2),
    # and the loss from cell i is C eps_i^(n + 1/2) F_i du
    w_dr = None
    if rec_te is None:
        c_dr = None
    elif rec_pow == -0.5:
        c_dr = nu_rec * math.sqrt(math.pi * rec_te) / 2.0
    else:
        c_dr = (nu_rec * math.gamma(1.5)
                / (math.gamma(1.5 + rec_pow) * rec_te ** rec_pow))
        w_dr = ec ** (rec_pow + 0.5)
    for it in range(1, maxit + 1):
        A, ec, eb, du = _with_ee(mix, 0.0, Tg, N, emax,
                                 nu_rec / N_gas if c_dr is None else 0.0, F,
                                 ne / N_gas if ee_on else 0.0, ne, "equal")
        if c_dr is not None and w_dr is None:
            A[np.arange(N), np.arange(N)] += c_dr / N_gas * du
        elif c_dr is not None:
            A[np.arange(N), np.arange(N)] += c_dr / N_gas * du * w_dr
        phi = np.linalg.solve(A * N_gas, q)
        ne_new = float(np.sum(np.sqrt(ec) * phi) * du)
        Fn = phi / ne_new
        Fn = 0.5 * (Fn + F)
        Fn /= float(np.sum(np.sqrt(ec) * Fn) * du)
        change = max(float(np.max(np.abs(Fn - F)) / np.max(F)),
                     abs(ne_new - ne) / ne)
        F, ne = Fn, 0.5 * (ne + ne_new)
        if change < tol:
            break
    else:
        raise RuntimeError(f"beam EEDF: no convergence ({change:.1e})")
    mean = float(np.sum(ec ** 1.5 * F) * du)
    return dict(F=F, ne=ne, mean=mean, ec=ec, du=du, emax=emax, N=N,
                iters=it, Tg=Tg, EN=0.0, nubar=nu_rec / N_gas)


def rate(c, sol):
    """<sigma v> of collision c over the solution's F0 [m^3/s]."""
    ec, F, du = sol["ec"], sol["F"], sol["du"]
    return L.GAMMA * float(np.sum(ec * L.interp(c, ec) * F)) * du


def maxwellian(ec, du, te):
    F = np.exp(-ec / te)
    return F / float(np.sum(np.sqrt(ec) * F) * du)
