#!/usr/bin/env python3
"""Figures 2 to 4 of detailed_derivation.tex (projects/IB).

  fig_lambda.pdf  Lambda(u) of Maxwell electrons with its asymptotes, and L0(w) of the
                  degenerate gas at T = 0 with its asymptotes            (Fig. 2)
  fig_rpa.pdf     the dynamic screening correction Lambda_RPA - Lambda    (Fig. 3)
  fig_strong.pdf  the exact quantum strong-field ratio R(beta, u) against Silin's G,
                  and the x-resolved integrand, exact and quasi-static     (Fig. 4)

The functions are the ones of detailed_derivation_check.py, which verifies the values
plotted here at the points of Tables 2, 4 and 5.  Colors: categorical slots 1-3 of the
dataviz reference palette (validated light, all pairs); asymptotes and reference curves
in the secondary ink.  PDFs are written without creation dates, so the bytes are
reproducible.

Run:  /usr/bin/python3 detailed_derivation_fig.py
"""
import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import jv, kve, erf, wofz

HERE = os.path.dirname(os.path.abspath(__file__))
EG = 0.5772156649015329
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e2dc'
META = {'CreationDate': None, 'ModDate': None}

plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9, 'legend.fontsize': 8,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'axes.linewidth': 0.6,
    'grid.color': GRID, 'grid.linewidth': 0.5, 'grid.linestyle': '-', 'lines.linewidth': 1.6,
    'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round', 'legend.frameon': False,
    'pdf.fonttype': 42, 'mathtext.fontset': 'cm',
})


def style(ax):
    ax.grid(True, which='major')
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def Lam(u):
    return -np.expm1(-u) / u * kve(0, u / 2)


def I0w(w):
    s1 = np.sqrt(1 + 2 * w)
    Im_ = (1 + w) * np.log((s1 + 1) / (s1 - 1)) - s1
    Ip_ = 0.0
    if w < 0.5:
        s2 = np.sqrt(1 - 2 * w)
        Ip_ = (1 - w) * np.log((1 + s2) / (1 - s2)) - s2
    return Im_ - Ip_


def Zf(z):
    return 1j * np.sqrt(np.pi) * wofz(z)


def eps_x(x, r, u):
    zp = (1 / x + u * x / 2) / np.sqrt(2)
    zm = (1 / x - u * x / 2) / np.sqrt(2)
    return 1 - r**2 * (Zf(zp) - Zf(zm)) / (np.sqrt(2) * u * x**3)


def Lam_rpa(r, u):
    """Eq. (Lambda-rpa), with the Langmuir pole as in detailed_derivation_check.py: inside x_pl (1 -+ 1e-3)
    the Lorentzian is integrated analytically when its width is below 1e-4 of that window."""
    f = lambda x: x**2 * (1 / eps_x(x, r, u)).imag
    xmax = 60.0 / u + 60.0
    edges = [1e-8, 0.05, 0.2, 0.5, 1, 2, 5, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4]
    edges = [q for q in edges if q < xmax] + [xmax]
    inner = 0.0
    if r < 1:
        g = lambda x: eps_x(x, r, u).real
        x0 = np.sqrt((1 / r**2 - 1) / 3)
        lo, hi = 0.3 * x0, min(1.6 * x0, 0.9)
        if g(lo) * g(hi) < 0:
            pole = brentq(g, lo, hi, xtol=1e-15)
            s_ = (g(pole * (1 + 1e-6)) - g(pole * (1 - 1e-6))) / (2e-6 * pole)
            g0 = eps_x(pole, r, u).imag
            width = abs(g0 / s_)
            d = 1e-3 * pole
            if width < 1e-4 * d:
                inner = -pole**2 * (2 / abs(s_)) * np.arctan(abs(s_) * d / g0)
            else:
                pts = sorted(set([pole] + [pole + sg * k * width for sg in (-1, 1) for k in (1, 3, 10, 30, 100, 1e3, 1e4)
                                           if k * width < d]))
                ie = [pole - d] + pts + [pole + d]
                inner = sum(quad(f, a, b, limit=500, epsabs=0, epsrel=1e-11)[0] for a, b in zip(ie[:-1], ie[1:]))
            edges = sorted(set([q for q in edges if not (pole - d < q < pole + d)] + [pole - d, pole + d]))
            tot = sum(quad(f, a, b, limit=500, epsabs=0, epsrel=1e-11)[0] for a, b in zip(edges[:-1], edges[1:])
                      if not (a == pole - d and b == pole + d))
            return -(1 / r**2) * np.sqrt(2 / np.pi) * (tot + inner)
    tot = sum(quad(f, a, b, limit=500, epsabs=0, epsrel=1e-11)[0] for a, b in zip(edges[:-1], edges[1:]))
    return -(1 / r**2) * np.sqrt(2 / np.pi) * tot


def I_all(lmax, z):
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


def S_exact(beta, u, nx=600):
    lmax = int(2.5 * 2 * beta**2 / u + 60)
    xmax = max(1.6 * np.sqrt(2 * lmax / u) + 12 / u, 40 / u) + 20
    x = np.geomspace(1e-2, xmax, nx)
    I = I_all(lmax, x * beta)
    tot = 0.0
    for l in range(1, lmax + 1):
        D = np.exp(-0.5 * (l / x - u * x / 2)**2) - np.exp(-0.5 * (l / x + u * x / 2)**2)
        tot += 2 * l * np.trapezoid(D * 2 / (x * beta) * I[l] / x**2, np.log(x))
    return tot


def R_exact(beta, u, nx=600):
    return 3 * S_exact(beta, u, nx) / (u * beta**2 * Lam(u))


def g_fun(s):
    return np.sqrt(np.pi) * erf(s) / s - 2 * np.exp(-s * s) if s > 1e-4 else 4 / 3 * s * s


def Gf(s0):
    avg = 2 / np.pi * quad(lambda th: g_fun(s0 * np.sin(th)), 0, np.pi / 2, limit=400)[0]
    return 3 * avg / (2 * s0**2)


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), metadata=META, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print('wrote', name)


T0 = time.time()
# ---------------------------------------------------------------- Fig. 2
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.6))
u = np.geomspace(1e-3, 30, 300)
a1.loglog(u, Lam(u), color=BLUE, label=r'$\Lambda(u)$')
us = u[u < 1.2]
a1.loglog(us, np.log(4 / (np.exp(EG) * us)), color=INK2, lw=1.0, ls='--')
ul = u[u > 0.3]
a1.loglog(ul, np.sqrt(np.pi) * ul**-1.5, color=INK2, lw=1.0, ls=':')
a1.text(1.5e-3, 3.3, r'$\ln(4/e^{\gamma}u)$', color=INK2, fontsize=8)
a1.text(4.0, 0.25, r'$\sqrt{\pi}\,u^{-3/2}$', color=INK2, fontsize=8)
a1.text(2.0e-2, 6.2, r'$\Lambda(u)$', color=INK, fontsize=8)
a1.set_xlabel(r'$u=\hbar\omega/k_BT$')
a1.set_ylabel(r'Coulomb logarithm $\Lambda$')
a1.set_ylim(5e-3, 12)
style(a1)
w = np.geomspace(1e-3, 1e2, 300)
L0 = np.array([I0w(x) / (2 * x) for x in w])
a2.loglog(w, L0, color=ORANGE)
ws = w[w < 1.0]
a2.loglog(ws, np.log(2 / ws), color=INK2, lw=1.0, ls='--')
wl = w[w > 0.2]
a2.loglog(wl, 4 / 3 * (2 * wl)**-1.5, color=INK2, lw=1.0, ls=':')
a2.text(1.5e-3, 3.0, r'$\ln(4E_F/\hbar\omega)$', color=INK2, fontsize=8)
a2.text(3.0, 0.12, r'$\frac{4}{3}(E_F/\hbar\omega)^{3/2}$', color=INK2, fontsize=8)
a2.text(1.5e-2, 6.0, r'$L_0(w)$, $T=0$', color=INK, fontsize=8)
a2.set_xlabel(r'$w=\hbar\omega/2E_F$')
a2.set_ylabel(r'$L_0$')
a2.set_ylim(1e-3, 12)
style(a2)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig_lambda.pdf')

# ---------------------------------------------------------------- Fig. 3
fig, ax = plt.subplots(figsize=(4.0, 2.7))
r = np.concatenate([np.linspace(0.05, 0.9, 18), [0.93, 0.95, 0.97, 0.98, 0.99, 0.995]])
for uu, col, mk, ms, lab in [(0.003, BLUE, 'o', 5.0, r'$u=0.003$'), (0.03, ORANGE, 's', 2.8, r'$u=0.03$'),
                             (0.3, AQUA, '^', 3.4, r'$u=0.3$')]:
    d = np.array([Lam_rpa(rr, uu) - Lam(uu) for rr in r])
    ax.plot(r, d, color=col, marker=mk, ms=ms, mew=0, label=lab)
ax.axhline(0, color=INK2, lw=0.6)
ax.set_xlabel(r'$r=\omega_p/\omega=(n_e/n_c)^{1/2}$')
ax.set_ylabel(r'$\Lambda_{\rm RPA}-\Lambda$')
ax.legend(loc='lower left')
style(ax)
fig.tight_layout()
save(fig, 'fig_rpa.pdf')

# ---------------------------------------------------------------- Fig. 4
fig, (b1, b2) = plt.subplots(1, 2, figsize=(6.5, 2.75))
bG = np.geomspace(0.05, 10, 60)
b1.loglog(bG, [Gf(b / np.sqrt(2)) for b in bG], color=INK2, lw=1.2, ls='--', label='Silin, fixed logarithm')
for uu, bmax, col, mk in [(0.3, 10.0, AQUA, '^'), (0.1, 6.0, ORANGE, 's'), (0.03, 3.0, BLUE, 'o')]:
    bb = np.geomspace(0.05, bmax, 16)
    R = [R_exact(b, uu) for b in bb]
    b1.loglog(bb, R, color=col, marker=mk, ms=3.2, mew=0, label=rf'exact, $u={uu}$')
b1.set_xlabel(r'$\beta=v_0/v_t$')
b1.set_ylabel(r'$\nu_{\rm eff}(v_0)/\nu_{\rm eff}(0)$')
b1.set_ylim(5e-3, 1.5)
b1.legend(loc='lower left')
style(b1)
# x-resolved integrand, beta = 8, u = 0.3
beta, uu = 8.0, 0.3
x = np.geomspace(0.05, 200, 220)
lmax = 1500
I = I_all(lmax, x * beta)
Iex = np.zeros_like(x)
for l in range(1, lmax + 1):
    D = np.exp(-0.5 * (l / x - uu * x / 2)**2) - np.exp(-0.5 * (l / x + uu * x / 2)**2)
    Iex += 2 * l * D * 2 / (x * beta) * I[l] / x**3
th = (np.arange(1200) + 0.5) / 1200 * (np.pi / 2)
mu, wmu = np.polynomial.legendre.leggauss(300)
mu = 0.5 * (mu + 1)
wmu = 0.5 * wmu
Iqs = np.zeros_like(x)
for m_, wm in zip(mu, wmu):
    nu = x[:, None] * beta * m_ * np.cos(th)[None, :]
    Dq = np.exp(-0.5 * (nu / x[:, None] - uu * x[:, None] / 2)**2) - np.exp(-0.5 * (nu / x[:, None] + uu * x[:, None] / 2)**2)
    Iqs += 2 * wm * (nu * Dq).mean(axis=1)
Iqs /= x**3
b2.semilogx(x, x * Iqs, color=INK2, lw=1.2, ls='--', label='quasi-static')
b2.semilogx(x, x * Iex, color=AQUA, label='exact sum')
b2.set_xlabel(r'$x=kv_t/\omega$')
b2.set_ylabel(r'$x\times$ integrand of $\mathcal{S}$')
b2.legend(loc='upper right')
b2.text(0.06, 0.07, r'$\beta=8,\ u=0.3$', color=INK2, fontsize=8)
style(b2)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig_strong.pdf')
print(f'done in {time.time() - T0:.1f} s')
