"""1Dv0_3's neutral gas: it flows, conducts heat and diffuses.

1Dv0_2 held every neutral still, since 1 us at 1 atm spreads the neutrals
by microns. On 1Dv0_3's graded mesh the cells beside the electrodes are a
few microns wide, and there the gas must move. The F2 burns out a band at
a time, and a burned band's gas heats toward 2000 K at constant volume
beside cells at 600 K. Here the gas obeys the 1-D Navier-Stokes equations
of the mixture along the gap,

    dn_s/dt + d/dx (n_s u + J_s) = 0                        every species
    d(rho u)/dt + d/dx (rho u^2 + p - tau) = 0
    dE/dt + d/dx ((E + p) u - tau u + q + sum_s J_s h_s) = 0,

with E = U_th + rho u^2/2, U_th = sum_s n_s cv_s k T over the heavy
species (LIGHT's own heat capacities, species.py), p = sum_s n_s k T,
tau = (4/3) mu du/dx, q = -kappa dT/dx, h_s = (cv_s + 1) k T. Each
species' formation energy moves with it, so the chemical energy needs no
flux of its own. The electrons ride with the gas, their mean energy with
them, and drift and diffuse relative to it in Fluid1D's transport.

The electrodes are impermeable and fixed (u = 0, no slip), and isothermal
at t_wall, or adiabatic with t_wall None. The heat they take is returned
per step, so the energy books close on it.

Transport, Chapman-Enskog in the first approximation, Lennard-Jones (12-6):
  viscosity  eta = (5/16) sqrt(pi m k T) / (pi sigma^2 Omega(2,2)*)
             (HCB eq. 8.2-18, Svehla 1962 eq. (1), Brokaw 1961 eq. (4))
  conductivity, monatomic: kappa = (15/4)(k/m) eta  (Brokaw eq. (7))
  binary diffusion  n D_12 = (3/16) sqrt(2 pi k T/m_12) / (pi sigma_12^2
             Omega(1,1)*), m_12 the reduced mass (HCB eq. 8.2-44, Brokaw
             eq. (22)), sigma_12 the mean and epsilon_12 the geometric
             mean of the pair's (Brokaw eqs. (23), (24))
  collision integrals: the fits of Neufeld, Janzen & Aziz, J. Chem.
             Phys. 57, 1100 (1972), as commonly quoted (the paper is not
             held here). A direct quadrature of the (12-6) definitions
             reproduces them to 0.2 % over T* = 0.5-50, and the tests hold
             them to 0.3 % over T* = 1-20 (tests/test_neutrals.py).
The gas is 99.5 % Kr, so its viscosity and conductivity are krypton's
(the bath). F2 and F diffuse in it, and Kr's own flux balances theirs so
diffusion carries no net mass. The other species have no Lennard-Jones
parameters (excited states, the excimers) and only ride with the gas.

Numerics: finite volumes on the mesh (mesh.py), MUSCL-Hancock with a van
Leer limiter on the primitive variables and the HLLC Riemann solver, the
viscous, conductive and diffusive fluxes centred, all explicit and
conservative, with substeps under the acoustic and diffusive limits.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np

from .model import KB
from .xsec import AMU

ATOM_AMU = {"Kr": 83.798, "F": 18.998403}

# Lennard-Jones (12-6) parameters, Svehla (1962) NASA TR R-132, Table I(a),
# pp. 22-24, with Table I(b)'s methods: sigma [Angstrom], epsilon/k [K]
LJ = {
    "Kr": (3.655, 178.9, "Svehla Table I(a) p.22: least-squares fit of "
           "measured viscosity (methods 1, 1; data refs 103-105)"),
    "F2": (3.357, 112.6, "Svehla Table I(a) p.22: least-squares fit of "
           "Franck & Stober's measured viscosity (methods 1, 1; data ref "
           "22; Aoyama & Kanda's set discarded, p.14)"),
    "F": (2.968, 112.6, "Svehla Table I(a) p.22: empirical rules 13 "
          "(sigma = 3/4 sigma(F2) + 0.45) and 15 (epsilon equal to F2's), an "
          "estimate"),
}
BATH = "Kr"


def omega11(t):
    """Omega(1,1)*(T*), Neufeld, Janzen & Aziz (1972)."""
    t = np.asarray(t, dtype=float)
    return (1.06036 * t ** -0.15610 + 0.19300 * np.exp(-0.47635 * t)
            + 1.03587 * np.exp(-1.52996 * t)
            + 1.76474 * np.exp(-3.89411 * t))


def omega22(t):
    """Omega(2,2)*(T*), Neufeld, Janzen & Aziz (1972)."""
    t = np.asarray(t, dtype=float)
    return (1.16145 * t ** -0.14874 + 0.52487 * np.exp(-0.77320 * t)
            + 2.16178 * np.exp(-2.43787 * t))


def mass_kg(atoms) -> float:
    return sum(ATOM_AMU[a] * k for a, k in atoms) * AMU


def viscosity(T, m, sigma_A, eps_k):
    """Pure-gas viscosity [Pa s] (HCB eq. 8.2-18)."""
    s = sigma_A * 1e-10
    return (5.0 / 16.0) * np.sqrt(math.pi * m * KB * T) / (
        math.pi * s * s * omega22(T / eps_k))


def n_diffusion(T, m1, m2, lj1, lj2):
    """n D_12 [m^-1 s^-1] for a binary pair, the pair's sigma the mean and
    epsilon the geometric mean (Brokaw eqs. (22)-(24), HCB eq. 8.2-44),
    independent of the density."""
    s12 = 0.5 * (lj1[0] + lj2[0]) * 1e-10
    e12 = math.sqrt(lj1[1] * lj2[1])
    mu = m1 * m2 / (m1 + m2)
    return (3.0 / 16.0) * np.sqrt(2.0 * math.pi * KB * T / mu) / (
        math.pi * s12 * s12 * omega11(T / e12))


class BathTransport:
    """The mixture's transport as the bath's (Kr), with F2 and F diffusing
    in it: viscosity, conductivity and n D_s per cell from T alone."""

    def __init__(self, names, masses):
        self.names = list(names)
        i_bath = self.names.index(BATH)
        self.m_bath = masses[i_bath]
        sb, eb, _ = LJ[BATH]
        self.lj_bath = (sb, eb)
        self.diffusing = [i for i, n in enumerate(self.names)
                          if n in LJ and n != BATH]
        self.lj = {i: LJ[self.names[i]][:2] for i in self.diffusing}
        self.masses = masses

    def viscosity(self, T):
        return viscosity(T, self.m_bath, *self.lj_bath)

    def conductivity(self, T):
        return 3.75 * (KB / self.m_bath) * self.viscosity(T)

    def n_diffusion(self, T):
        """{species column: n D in the bath} for the diffusing species."""
        return {i: n_diffusion(T, self.masses[i], self.m_bath, self.lj[i],
                               self.lj_bath)
                for i in self.diffusing}


class ConstantTransport:
    """Constant mu [Pa s], kappa [W/m/K] and n D [1/m/s] (the tests')."""

    def __init__(self, mu=0.0, kappa=0.0, n_d=None):
        self.mu, self.kappa, self.n_d = mu, kappa, dict(n_d or {})
        self.diffusing = sorted(self.n_d)

    def viscosity(self, T):
        return np.full_like(T, self.mu)

    def conductivity(self, T):
        return np.full_like(T, self.kappa)

    def n_diffusion(self, T):
        return {i: np.full_like(T, v) for i, v in self.n_d.items()}


def _dot(A, w):
    """A @ w as multiply-and-sum: numpy 2's matmul on macOS Accelerate
    warns spuriously over operands that hold zeros (LIGHT's Known-bad
    list), and the zeros are common here (empty species)."""
    return (A * w).sum(axis=-1)


def _vanleer(a, b):
    """van Leer's harmonic slope of two one-sided slopes (0 at extrema)."""
    prod = a * b
    out = np.zeros_like(a)
    ok = prod > 0.0
    out[ok] = 2.0 * prod[ok] / (a[ok] + b[ok])
    return out


class NeutralGas:
    """The gas on a mesh.

    names, masses [kg], cv [k per particle], heavy (mask: in the gas's
    mass and pressure; the electrons are not), e_col: the electrons'
    column (their mean energy rides with them). bath: the species whose
    diffusive flux balances the others' (Kr). transport: BathTransport or
    ConstantTransport. t_wall: the electrodes' temperature [K], or None
    for adiabatic walls."""

    CFL = 0.5
    DIFFUSION_SAFETY = 0.4

    def __init__(self, mesh, names, masses, cv, heavy, e_col=None,
                 bath=BATH, transport=None, t_wall=None):
        self.mesh = mesh
        self.names = tuple(names)
        self.m = np.asarray(masses, dtype=float)
        self.cv = np.asarray(cv, dtype=float)
        self.heavy = np.asarray(heavy, dtype=bool)
        self.e_col = e_col
        self.i_bath = self.names.index(bath) if bath in self.names else None
        self.transport = transport
        self.t_wall = t_wall
        h = np.asarray(mesh.widths, dtype=float)
        self.h = h
        self.d = np.asarray(mesh.spacing, dtype=float)       # N - 1 faces
        wl, wr = mesh.weights
        self.wl, self.wr = np.asarray(wl), np.asarray(wr)
        hmin = h.copy()
        hmin[:-1] = np.minimum(hmin[:-1], self.d)
        hmin[1:] = np.minimum(hmin[1:], self.d)
        self.hmin = hmin
        self.m_heavy = np.where(self.heavy, self.m, 0.0)
        self.cv_heavy = np.where(self.heavy, self.cv, 0.0)
        self.one_heavy = self.heavy.astype(float)

    # ----------------------------------------------------------- state
    def primitives(self, dens, T, u):
        rho = _dot(dens, self.m_heavy)
        n = _dot(dens, self.one_heavy)
        C = _dot(dens, self.cv_heavy)              # sum n cv (units of k)
        p = n * KB * T
        gamma = 1.0 + n / C
        return rho, n, C, p, gamma

    def energy(self, dens, T, u):
        """Thermal + kinetic energy per cell [J/m^3]."""
        rho, n, C, p, gamma = self.primitives(dens, T, u)
        return C * KB * T + 0.5 * rho * u * u

    def stable_dt(self, dens, T, u):
        rho, n, C, p, gamma = self.primitives(dens, T, u)
        c = np.sqrt(gamma * p / rho)
        dt = self.CFL * np.min(self.h / (np.abs(u) + c))
        if self.transport is not None:
            mu = self.transport.viscosity(T)
            kap = self.transport.conductivity(T)
            chi = np.maximum(4.0 / 3.0 * mu / rho, kap / (C * KB))
            for i, nd in self.transport.n_diffusion(T).items():
                chi = np.maximum(chi, nd / n)
            if np.any(chi > 0.0):
                dt = min(dt, self.DIFFUSION_SAFETY * float(
                    np.min(self.hmin ** 2 / (2.0 * np.maximum(chi, 1e-300)))))
        return dt

    # ------------------------------------------------------------ step
    def step(self, dens, T, u, eps=None, dt=0.0):
        """Advance by dt [s]: dens (cells, species) [m^-3], T [K], u [m/s],
        eps the electrons' mean energy [eV] (or None). Returns a namespace:
        dens, T, u, eps, wall_heat (J/m^2 conducted into the left and right
        electrodes), substeps."""
        dens = np.array(dens, dtype=float)
        T = np.array(T, dtype=float)
        u = np.array(u, dtype=float)
        eps = None if eps is None else np.array(eps, dtype=float)
        wall = np.zeros(2)
        if dt <= 0.0:
            return SimpleNamespace(dens=dens, T=T, u=u, eps=eps,
                                   wall_heat=wall, substeps=0)
        n_sub = max(1, int(math.ceil(dt / self.stable_dt(dens, T, u))))
        h = dt / n_sub
        count = [0]
        for _ in range(n_sub):
            dens, T, u, eps, w = self._safe_substep(dens, T, u, eps, h, count)
            wall += w
        return SimpleNamespace(dens=dens, T=T, u=u, eps=eps, wall_heat=wall,
                               substeps=count[0])

    def _safe_substep(self, dens, T, u, eps, dt, count, depth=0):
        """A substep, split in two where it would leave a density or an
        energy non-positive (a sharp front the limits did not foresee)."""
        try:
            out = self._substep(dens, T, u, eps, dt)
            count[0] += 1
            return out
        except FloatingPointError:
            if depth >= 8:
                raise
        wall = np.zeros(2)
        for _ in range(2):
            dens, T, u, eps, w = self._safe_substep(dens, T, u, eps, 0.5 * dt,
                                                    count, depth + 1)
            wall += w
        return dens, T, u, eps, wall

    def _faces(self, W, odd):
        """Limited face values (minus side, plus side) of the cell values
        W (cells, k): one-sided slopes on the mesh, mirror ghosts at the
        walls (odd columns, the velocity, change sign there), van Leer, and
        the faces held within their cell's neighbours."""
        h, d = self.h, self.d
        sign = np.where(odd, -1.0, 1.0)
        ghost_l = W[0] * sign
        ghost_r = W[-1] * sign
        ext = np.vstack([ghost_l, W, ghost_r])
        spacing = np.concatenate([[h[0]], d, [h[-1]]])[:, None]
        delta = (ext[1:] - ext[:-1]) / spacing          # N + 1 faces
        slope = _vanleer(delta[:-1], delta[1:])
        half = 0.5 * h[:, None] * slope
        lo = np.minimum(np.minimum(ext[:-2], ext[1:-1]), ext[2:])
        hi = np.maximum(np.maximum(ext[:-2], ext[1:-1]), ext[2:])
        minus = np.clip(W - half, lo, hi)
        plus = np.clip(W + half, lo, hi)
        return minus, plus, slope

    def _hllc(self, rL, uL, pL, gL, rR, uR, pR, gR):
        """HLLC fluxes (mass, momentum, energy) and which side is upwind of
        the contact (True: left)."""
        EL = pL / (gL - 1.0) + 0.5 * rL * uL * uL
        ER = pR / (gR - 1.0) + 0.5 * rR * uR * uR
        cL = np.sqrt(gL * pL / rL)
        cR = np.sqrt(gR * pR / rR)
        SL = np.minimum(uL - cL, uR - cR)
        SR = np.maximum(uL + cL, uR + cR)
        den = rL * (SL - uL) - rR * (SR - uR)
        Ss = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / den
        FL = np.stack([rL * uL, rL * uL * uL + pL, (EL + pL) * uL])
        FR = np.stack([rR * uR, rR * uR * uR + pR, (ER + pR) * uR])
        UL = np.stack([rL, rL * uL, EL])
        UR = np.stack([rR, rR * uR, ER])

        def star(r, uu, p, E, S):
            gap = S - Ss
            gap = np.where(np.abs(gap) > 0.0, gap, 1e-300)
            f = r * (S - uu) / gap
            return np.stack([f, f * Ss,
                             f * (E / r + (Ss - uu) * (Ss + p / (r * (S - uu))))])
        UsL = star(rL, uL, pL, EL, SL)
        UsR = star(rR, uR, pR, ER, SR)
        F = np.where(SL >= 0.0, FL,
                     np.where(Ss >= 0.0, FL + SL * (UsL - UL),
                              np.where(SR >= 0.0, FR + SR * (UsR - UR), FR)))
        return F, Ss >= 0.0

    def _substep(self, dens, T, u, eps, dt):
        h, d, wl, wr = self.h, self.d, self.wl, self.wr
        N = h.size
        rho, n, C, p, gamma = self.primitives(dens, T, u)
        phi = dens / rho[:, None]                         # per kg
        has_e = eps is not None and self.e_col is not None
        cols = [rho[:, None], u[:, None], p[:, None], phi]
        if has_e:
            cols.append(eps[:, None])
        W = np.hstack(cols)
        odd = np.zeros(W.shape[1], dtype=bool)
        odd[1] = True
        minus, plus, slope = self._faces(W, odd)

        # Hancock: both faces of a cell move half a step with its slopes
        S = phi.shape[1]
        s_r, s_u, s_p = slope[:, 0], slope[:, 1], slope[:, 2]
        dW = np.zeros_like(W)
        dW[:, 0] = -(u * s_r + rho * s_u)
        dW[:, 1] = -(u * s_u + s_p / rho)
        dW[:, 2] = -(u * s_p + gamma * p * s_u)
        dW[:, 3:] = -u[:, None] * slope[:, 3:]
        minus = minus + 0.5 * dt * dW
        plus = plus + 0.5 * dt * dW
        for A in (minus, plus):
            # a face the predictor took to a non-positive density or
            # pressure falls back to its cell's value (first order there)
            bad = (A[:, 0] <= 0.0) | (A[:, 2] <= 0.0)
            if np.any(bad):
                A[bad] = W[bad]
            # no species below zero at a face (the predictor can step past
            # the limiter's range by round-off), and the composition
            # renormalised so its heavy mass fractions sum to one: the
            # species fluxes then carry exactly the mass flux
            np.maximum(A[:, 3:3 + S], 0.0, out=A[:, 3:3 + S])
            A[:, 3:3 + S] /= _dot(A[:, 3:3 + S], self.m_heavy)[:, None]

        def gamma_of(Wf):
            ph = Wf[:, 3:3 + S]
            return 1.0 + _dot(ph, self.one_heavy) / _dot(ph, self.cv_heavy)

        # interior faces: left state = plus side of cell i, right = minus
        # side of cell i + 1
        L, R = plus[:-1], minus[1:]
        F, left_up = self._hllc(L[:, 0], L[:, 1], L[:, 2], gamma_of(L),
                                R[:, 0], R[:, 1], R[:, 2], gamma_of(R))
        up = np.where(left_up[:, None], L, R)
        f_rho = F[0]
        flux_n = f_rho[:, None] * up[:, 3:3 + S]        # every species
        flux_m = F[1].copy()
        flux_E = F[2].copy()
        flux_eps = (flux_n[:, self.e_col] * up[:, 3 + S]) if has_e else None

        # the walls: no mass through them, the mirror state's pressure
        r0, u0, p0 = minus[0, 0], minus[0, 1], minus[0, 2]
        c0 = math.sqrt(gamma_of(minus[:1])[0] * p0 / r0)
        rN, uN, pN = plus[-1, 0], plus[-1, 1], plus[-1, 2]
        cN = math.sqrt(gamma_of(plus[-1:])[0] * pN / rN)
        pw_l = p0 - r0 * u0 * (abs(u0) + c0 - u0)
        pw_r = pN + rN * uN * (abs(uN) + cN + uN)

        wall_m = np.array([pw_l, pw_r])
        wall_E = np.zeros(2)

        # viscous, conductive and diffusive fluxes, centred, from the
        # cell values at the start of the step
        if self.transport is not None:
            mu = self.transport.viscosity(T)
            kap = self.transport.conductivity(T)

            def series(c):
                return d / (0.5 * h[:-1] / c[:-1] + 0.5 * h[1:] / c[1:])
            mu_f = series(np.maximum(mu, 1e-300))
            kap_f = series(np.maximum(kap, 1e-300))
            tau = (4.0 / 3.0) * mu_f * (u[1:] - u[:-1]) / d
            q = -kap_f * (T[1:] - T[:-1]) / d
            u_f = wl * u[:-1] + wr * u[1:]
            T_f = wl * T[:-1] + wr * T[1:]
            flux_m -= tau
            flux_E += q - tau * u_f
            nd = self.transport.n_diffusion(T)
            if nd:
                x = dens / n[:, None]
                J = np.zeros((N - 1, S))
                for i, ndi in nd.items():
                    J[:, i] = -series(np.maximum(ndi, 1e-300)) * (
                        x[1:, i] - x[:-1, i]) / d
                if self.i_bath is not None:
                    J[:, self.i_bath] = -_dot(J, self.m_heavy) / self.m[
                        self.i_bath]
                flux_n += J
                flux_E += _dot(J, (self.cv_heavy + 1.0) * self.one_heavy) \
                    * KB * T_f
            # no slip at the walls, the stress over half a cell
            wall_m -= (4.0 / 3.0) * np.array([
                mu[0] * u[0] / (0.5 * h[0]),
                mu[-1] * (0.0 - u[-1]) / (0.5 * h[-1])])
            if self.t_wall is not None:
                wall_E = np.array([
                    -kap[0] * (T[0] - self.t_wall) / (0.5 * h[0]),
                    -kap[-1] * (self.t_wall - T[-1]) / (0.5 * h[-1])])

        # the update: every face's flux in and out of its two cells
        def divergence(interior, left, right):
            faces = np.concatenate([np.atleast_1d(left), interior,
                                    np.atleast_1d(right)], axis=0)
            return (faces[1:] - faces[:-1]) / (h if faces.ndim == 1
                                               else h[:, None])
        zero = np.zeros(S)
        dens_new = dens - dt * divergence(flux_n, zero[None, :],
                                          zero[None, :])
        mom = rho * u - dt * divergence(flux_m, wall_m[0], wall_m[1])
        dE = -dt * divergence(flux_E, wall_E[0], wall_E[1])
        rho_new = _dot(dens_new, self.m_heavy)
        if np.any(rho_new <= 0.0) or np.any(dens_new[:, self.heavy] < 0.0):
            i, s = np.unravel_index(np.argmin(np.where(
                self.heavy, dens_new, np.inf)), dens_new.shape)
            raise FloatingPointError(
                f"a gas density went negative: {self.names[s]} in cell {i}, "
                f"{dens_new[i, s]:.3e} m^-3 (was {dens[i, s]:.3e})")
        u_new = mom / rho_new
        C_new = _dot(dens_new, self.cv_heavy)
        dK = 0.5 * (rho_new * u_new * u_new - rho * u * u)
        # U_th' = U_th + dE - dK, written so a cell where nothing happens
        # keeps its temperature to the bit
        T_new = T * (C / C_new) + (dE - dK) / (C_new * KB)
        if np.any(T_new <= 0.0):
            raise FloatingPointError("the gas energy went non-positive")
        eps_new = eps
        if has_e:
            ne0 = dens[:, self.e_col]
            ne1 = dens_new[:, self.e_col]
            N_eps = ne0 * eps - dt * divergence(flux_eps, 0.0, 0.0)
            live = ne1 > 0.0
            eps_new = eps.copy()
            eps_new[live] = N_eps[live] / ne1[live]
        # heat into the electrodes [J/m^2] over the step (positive: out of
        # the gas)
        wall_heat = dt * np.array([-wall_E[0], wall_E[1]])
        return dens_new, T_new, u_new, eps_new, wall_heat
