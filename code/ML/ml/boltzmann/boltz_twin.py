#!/usr/bin/env python3
"""boltz_twin.py -- derive.py's two-term Boltzmann solver: the discretization lxcat.html uses, solved densely
(numpy.linalg.solve on the whole matrix, superelastic and growth terms implicit), written apart from the page.

The equation (Hagelaar & Pitchford 2005, PSST 14 722, the form LoKI-B solves), steady state, DC field, the isotropic
part F0(eps) of the EEDF normalised as int sqrt(eps) F0 d eps = 1:

    d/d eps [ W F0 - D dF0/d eps ] = S
    W = -gamma eps^2 sigma_eps,        sigma_eps = sum_g x_g (2 m/M_g) sigma_el,g
    D = gamma/3 (E/N)^2 eps / sigma~ + gamma (k Tg/e) eps^2 sigma_eps
    sigma~ = sum_g x_g [sigma_el,g + sum_c p_c sigma_c + sum_c,rev p'_c sigma_c,super] + (nu/N)/(gamma sqrt eps)
    S = the inelastic terms (loss at eps, the scattered electrons appearing at eps - u) - (nu/N) sqrt(eps) F0

Discretization (the page's, stated in its calculation panel):
  * cells i = 0..N-1 of width du, centres e_i = (i + 1/2) du, F_i = F0(e_i); zero flux through 0 and through e_max;
  * the flux through the boundary b = (i+1) du by the exponential (Scharfetter-Gummel) scheme:
        Gamma_b = (D_b/du) [ B(-z) F_i - B(z) F_{i+1} ],  z = W_b du / D_b,  B(x) = x/(e^x - 1);
  * a collision c of weight w_c (= fraction x of its gas times the population of its target) removes
        L_cj = gamma w_c e_j sigma_c(e_j) F_j du
    from cell j and puts the scattered electron(s) at e_j - u (excitation; with equal sharing both electrons at
    (e_j - u)/2, with one-takes-all the second at 0) split linearly between the two cells whose centres bracket that
    energy (into cell 0 below e_0): electrons are conserved exactly, and the mean energy of the deposit is exact;
  * a reversible collision whose product is populated (p') adds the de-excitation: from cell j at the rate
        gamma x p' (g_t/g_p) (e_j + u) sigma_c(e_j + u) F_j du, deposited at e_j + u (into the top cell above e_max);
  * attachment removes; temporal growth removes (nu/N) sqrt(e_i) F_i du from every cell, with nu/N the net
    ionization (ionization less attachment) of the current F, iterated to self-consistency;
  * the bottom cell's balance is replaced by F_0 = 1 (the system is singular: its columns sum to zero without
    growth, and the dropped balance holds at the self-consistent nu), then F is normalised.
Cross sections are evaluated as LoKI-MC evaluates them (lx_lib.interp): linear, zero outside the table and at or
below a threshold.
"""
import math
import numpy as np
import lx_lib as L

TD = 1e-21   # V m^2


def bern(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-3
    xs = x[small]
    out[small] = 1 - xs / 2 + xs * xs / 12 - xs ** 4 / 720
    xl = x[~small]
    with np.errstate(over='ignore'):
        out[~small] = xl / np.expm1(xl)
    return out


def deposit_weights(pos):
    """position in cell-index units (centre of cell i at i): [(cell, weight), ...]"""
    if pos <= 0:
        return [(0, 1.0)]
    i0 = int(math.floor(pos))
    w = pos - i0
    return [(i0, 1.0 - w), (i0 + 1, w)]


def assemble(mix, EN, Tg, N, emax, nubar, ion='equal'):
    """the matrix of the page's discretization (without the F_0 = 1 row); mix = [(gas, x), ...]"""
    du = emax / N
    ec = (np.arange(N) + 0.5) * du
    eb = (np.arange(1, N)) * du          # the N-1 inner boundaries
    G = L.GAMMA
    EN_V = EN * TD
    sig_eps_b = np.zeros(N - 1)
    sig_m_b = np.zeros(N - 1)
    A = np.zeros((N, N))
    for gas, x in mix:
        sig_el_b = gas['sigma_el'](eb)
        sig_eps_b += x * 2 * L.ME / gas['M'] * sig_el_b
        sig_m_b += x * sig_el_b
        for c in gas['colls']:
            if c['w'] > 0:
                sig_m_b += x * c['w'] * L.interp(c, eb)
            if c.get('wsup', 0) > 0:
                sig_m_b += x * c['wsup'] * L.superelastic(c, eb, c['gratio'])
    sig_t = sig_m_b + nubar / (G * np.sqrt(eb))
    W = -G * eb ** 2 * sig_eps_b
    D = G / 3 * EN_V ** 2 * eb / sig_t + G * (L.KB * Tg / L.QE) * eb ** 2 * sig_eps_b
    z = W * du / D
    a = D / du * bern(-z)
    b = -D / du * bern(z)
    for i in range(N - 1):
        A[i, i] += a[i]
        A[i, i + 1] += b[i]
        A[i + 1, i] -= a[i]
        A[i + 1, i + 1] -= b[i]
    for gas, x in mix:
        for c in gas['colls']:
            wt = x * c['w']
            u = c['thr']
            if wt > 0:
                sc = L.interp(c, ec)
                for j in np.nonzero(sc)[0]:
                    Lj = G * wt * ec[j] * sc[j] * du
                    A[j, j] += Lj
                    if c['type'] == 'Attachment':
                        continue
                    if c['type'] == 'Ionization' and ion == 'equal':
                        for k, wk in deposit_weights((ec[j] - u) / 2 / du - 0.5):
                            A[k, j] -= 2 * wk * Lj
                    else:
                        for k, wk in deposit_weights((ec[j] - u) / du - 0.5):
                            A[k, j] -= wk * Lj
                        if c['type'] == 'Ionization':    # one-takes-all: the second electron at zero energy
                            A[0, j] -= Lj
            ws = x * c.get('wsup', 0)
            if ws > 0:
                sup = c['gratio'] * (ec + u) * L.interp(c, ec + u)    # e * sigma_super(e)
                for j in np.nonzero(sup)[0]:
                    Lj = G * ws * sup[j] * du
                    A[j, j] += Lj
                    for k, wk in deposit_weights((ec[j] + u) / du - 0.5):
                        A[min(k, N - 1), j] -= wk * Lj
    if nubar != 0:
        A[np.arange(N), np.arange(N)] += nubar * np.sqrt(ec) * du
    return A, ec, eb, du, sig_t, sig_eps_b


def net_ion(mix, F, ec, du):
    """nu/N = sum of ionization rates less attachment rates of F (the solver's own quadrature)"""
    s = 0.0
    for gas, x in mix:
        for c in gas['colls']:
            if c['type'] in ('Ionization', 'Attachment') and c['w'] > 0:
                k = L.GAMMA * float(np.sum(ec * L.interp(c, ec) * F)) * du
                s += (1 if c['type'] == 'Ionization' else -1) * x * c['w'] * k
    return s


def solve_at(mix, EN, Tg, N, emax, nubar, ion):
    A, ec, eb, du, sig_t, sig_eps = assemble(mix, EN, Tg, N, emax, nubar, ion)
    A[0, :] = 0
    A[0, 0] = 1
    rhs = np.zeros(N)
    rhs[0] = 1
    F = np.linalg.solve(A, rhs)
    F = F / float(np.sum(np.sqrt(ec) * F) * du)
    return F, ec, du, float(sig_t.min())


def solve_fixed(mix, EN, Tg, N, emax, ion='equal', growth=True, tol=1e-13, maxit=200):
    """nu/N = netIon(F(nu/N)) by Brent's method on h(nu) = netIon(F(nu)) - nu (the page uses secant steps: a different
    route to the same root); plain iteration diverges where attachment dominates (F2)"""
    from scipy.optimize import brentq
    F, ec, du, smin = solve_at(mix, EN, Tg, N, emax, 0.0, ion)
    g0 = net_ion(mix, F, ec, du) if growth else 0.0
    its = [1]
    if growth and g0 != 0.0:
        def h(nu):
            Fh, ech, duh, sm = solve_at(mix, EN, Tg, N, emax, nu, ion)
            its[0] += 1
            return net_ion(mix, Fh, ech, duh) - nu
        a, ha = 0.0, g0
        b = g0
        hb = h(b)
        k = 0
        while ha * hb > 0:
            k += 1
            if k > 80:
                raise RuntimeError('growth: no bracket')
            a, ha = b, hb
            b = b * 2
            hb = h(b)
        nubar = brentq(h, min(a, b), max(a, b), xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=200)
        F, ec, du, smin = solve_at(mix, EN, Tg, N, emax, nubar, ion)
        nubar = net_ion(mix, F, ec, du)
    else:
        nubar = g0
    return dict(F=F, ec=ec, du=du, emax=emax, N=N, nubar=nubar, iters=its[0], EN=EN, Tg=Tg, sigt_min=smin)


def swarm(mix, sol):
    """the page's swarm parameters from a solution"""
    F, ec, du, N = sol['F'], sol['ec'], sol['du'], sol['N']
    EN_V = sol['EN'] * TD
    G = L.GAMMA
    A, ec, eb, du, sig_t, sig_eps = assemble(mix, sol['EN'], sol['Tg'], N, sol['emax'], sol['nubar'])
    mean = float(np.sum(ec ** 1.5 * F) * du)
    muN = -G / 3 * float(np.sum(eb / sig_t * (F[1:] - F[:-1])))
    # D N at the cell centres: sigma~ there
    sig_c = np.zeros(N)
    for gas, x in mix:
        sig_c += x * gas['sigma_el'](ec)
        for c in gas['colls']:
            if c['w'] > 0:
                sig_c += x * c['w'] * L.interp(c, ec)
            if c.get('wsup', 0) > 0:
                sig_c += x * c['wsup'] * L.superelastic(c, ec, c['gratio'])
    sig_c = sig_c + sol['nubar'] / (G * np.sqrt(ec))
    DN = G / 3 * float(np.sum(ec / sig_c * F) * du)
    w = muN * EN_V
    rates = []
    kion = katt = 0.0
    for gas, x in mix:
        for c in gas['colls']:
            k = G * float(np.sum(ec * L.interp(c, ec) * F)) * du
            rates.append((gas['id'], c['desc'], k))
            if c['type'] == 'Ionization':
                kion += x * c['w'] * k
            if c['type'] == 'Attachment':
                katt += x * c['w'] * k
    return dict(mean=mean, muN=muN, w=w, DN=DN, char=DN / muN, kion=kion, katt=katt,
                alphaN=kion / w if w > 0 else 0.0, etaN=katt / w if w > 0 else 0.0, rates=rates)


def decades_of(F):
    return math.log10(F.max() / max(abs(F[-1]), 1e-300))


def auto_solve(mix, EN, Tg, N, emax0, ion='equal', emax_cap=900.0):
    """the page's grid rule: e_max so that F falls between 11 and 15 decades from its peak to the top cell; steps of
    13/decades clamped to [1/2, 2] until a short (<11) and a long (>15) grid are both known, then bisection in log
    e_max; within 2 % the longer is kept; at most 20 tries; never past emax_cap"""
    emax = min(emax0, emax_cap)
    hist = []
    lo = hi = sol_hi = None
    at_cap = False
    for t in range(20):
        try:
            sol = solve_fixed(mix, EN, Tg, N, emax, ion)
        except RuntimeError:
            sol = solve_fixed(mix, EN, Tg, N, emax, ion, growth=False)
            sol['nogrowth'] = True
        dec = decades_of(sol['F'])
        hist.append((emax, dec))
        if 11 <= dec <= 15:
            break
        if dec < 11:
            if emax >= emax_cap:
                at_cap = True
                break
            lo = emax if lo is None else max(lo, emax)
        elif hi is None or emax < hi:
            hi, sol_hi = emax, sol
        if lo is not None and hi is not None:
            if hi / lo < 1.02:
                sol = sol_hi
                dec = decades_of(sol['F'])
                break
            emax = math.sqrt(lo * hi)
        else:
            emax = min(emax * min(max(13 / max(dec, 0.5), 0.5), 2.0), emax_cap)
    if not (11 <= dec <= 15) and not at_cap and sol_hi is not None:
        sol = sol_hi
        dec = decades_of(sol['F'])
    sol['decades'] = dec
    sol['hist'] = hist
    sol['at_cap'] = at_cap
    sol['invalid'] = not (sol['F'].min() > -1e-8 * sol['F'].max())
    return sol
