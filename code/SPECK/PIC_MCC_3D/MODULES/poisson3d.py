"""MISSING MODULE 1 -- the 3-D electrostatic Poisson solve.

    laplacian(phi) = -rho / eps0

This is the blocker named in ``../PIC_MCC_3D.txt`` section 2.1: without it
there is no sheath, no cathode fall, no ambipolar field and no breakdown, and
the electromagnetic (Yee/FDTD) path cannot substitute because its timestep is
Courant-limited by the speed of light on a cell that must resolve lambda_D
(124x-601x more steps over Te = 4 eV .. 0.17 eV; see ``units.es_over_em_dt_ratio``).

THE METHOD (why this is affordable in 3-D)
------------------------------------------
The 7-point stencil on a uniform grid is SEPARABLE.  In the periodic (x, y)
plane the discrete Laplacian is diagonalised exactly by the DFT: a plane wave
exp(2 pi i (p x_i / nx + q y_j / ny)) is an eigenvector with eigenvalue

    lam_p = (2 cos(2 pi p / nx) - 2) / h^2        (and likewise for q)

So after ONE 2-D FFT of rho over (x, y), each (p, q) mode is an independent
TRIDIAGONAL system in z, which a Thomas algorithm solves in O(nz).  Vectorised
over all nx*ny modes at once this is O(N log N) total, exact to round-off, with
no iteration, no preconditioner and no dependency beyond numpy.

Two boundary modes:

  * ``"periodic_xy_dirichlet_z"`` -- THE DEVICE GEOMETRY.  Electrodes at
    z = 0 and z = (nz-1) h hold prescribed potentials (scalar or an (nx, ny)
    profile); periodic transversely.  Non-singular, so no neutrality
    requirement.
  * ``"periodic"`` -- THE VERIFICATION GEOMETRY.  Fully spectral in all three
    directions.  The k = 0 mode of the discrete Laplacian is EXACTLY singular,
    so a net charge has no periodic solution: :meth:`solve` raises rather than
    silently zeroing the mode and returning a wrong answer (Known-bad list).

``reference_solve_lu`` assembles the same operator explicitly with
scipy.sparse and factorises it.  It is O(N^2)-ish in memory and is used ONLY
by the gates, where it is the oracle: the fast path must reproduce it to
round-off.

Node layout: phi and rho live on nodes, shape (nx, ny, nz), spacing h.  In the
periodic directions node nx is node 0 (not stored).  In the Dirichlet
direction nodes k = 0 and k = nz-1 ARE the electrodes and their rho is
ignored.
"""

from __future__ import annotations

import numpy as np

from .units import EPS0

PERIODIC_XY_DIRICHLET_Z = "periodic_xy_dirichlet_z"
PERIODIC = "periodic"
_BCS = (PERIODIC_XY_DIRICHLET_Z, PERIODIC)


def _dft_eigenvalues(n: int, h: float) -> np.ndarray:
    """Eigenvalues of the periodic 1-D second difference, DFT-ordered."""
    k = np.arange(n)
    return (2.0 * np.cos(2.0 * np.pi * k / n) - 2.0) / (h * h)


def _thomas_batch(diag: np.ndarray, off: float, rhs: np.ndarray) -> np.ndarray:
    """Solve a batch of tridiagonal systems sharing constant off-diagonals.

    ``diag``  (nx, ny)      -- the (mode-dependent) main diagonal
    ``off``   scalar        -- sub- and super-diagonal (symmetric)
    ``rhs``   (nx, ny, m)   -- right-hand sides
    returns   (nx, ny, m)

    Forward sweep then back substitution, vectorised across all modes.  The
    systems are weakly diagonally dominant (|diag| >= 2|off|, with equality
    only at the (0,0) mode), where the Thomas algorithm is stable without
    pivoting -- the LU gate is what proves that claim on real numbers.
    """
    m = rhs.shape[-1]
    cp = np.empty(rhs.shape, dtype=rhs.dtype)
    dp = np.empty(rhs.shape, dtype=rhs.dtype)
    cp[..., 0] = off / diag
    dp[..., 0] = rhs[..., 0] / diag
    for k in range(1, m):
        denom = diag - off * cp[..., k - 1]
        cp[..., k] = off / denom
        dp[..., k] = (rhs[..., k] - off * dp[..., k - 1]) / denom
    x = np.empty(rhs.shape, dtype=rhs.dtype)
    x[..., m - 1] = dp[..., m - 1]
    for k in range(m - 2, -1, -1):
        x[..., k] = dp[..., k] - cp[..., k] * x[..., k + 1]
    return x


class Poisson3D:
    """Separable 3-D electrostatic Poisson solver (see the module docstring)."""

    def __init__(self, nx: int, ny: int, nz: int, h: float,
                 bc: str = PERIODIC_XY_DIRICHLET_Z,
                 neutrality_tol: float = 1e-9) -> None:
        if bc not in _BCS:
            raise ValueError(f"unknown bc {bc!r}; choose from {_BCS}")
        if min(nx, ny, nz) < 3:
            raise ValueError(f"grid {nx}x{ny}x{nz}: need at least 3 nodes/axis")
        if h <= 0:
            raise ValueError(f"h = {h}: cell size must be positive")
        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
        self.h = float(h)
        self.bc = bc
        self.neutrality_tol = float(neutrality_tol)
        self.shape = (self.nx, self.ny, self.nz)

        lam_x = _dft_eigenvalues(self.nx, self.h)
        lam_y = _dft_eigenvalues(self.ny, self.h)
        self._lam_xy = lam_x[:, None] + lam_y[None, :]          # (nx, ny)
        if bc == PERIODIC:
            lam_z = _dft_eigenvalues(self.nz, self.h)
            self._lam_xyz = self._lam_xy[:, :, None] + lam_z[None, None, :]
        else:
            self._off = 1.0 / (self.h * self.h)
            # main diagonal per mode: -2/h^2 + lam_xy
            self._diag = self._lam_xy - 2.0 / (self.h * self.h)

    # ------------------------------------------------------------------ solve
    def solve(self, rho: np.ndarray, phi_lo=0.0, phi_hi=0.0) -> np.ndarray:
        """Return phi [V] for a charge density ``rho`` [C/m^3] on the nodes.

        ``phi_lo`` / ``phi_hi`` are the z = 0 / z = Lz electrode potentials
        (scalar or (nx, ny)); they are ignored in the fully periodic mode.
        """
        rho = np.asarray(rho, dtype=float)
        if rho.shape != self.shape:
            raise ValueError(f"rho shape {rho.shape} != grid {self.shape}")
        if self.bc == PERIODIC:
            return self._solve_periodic(rho)
        return self._solve_dirichlet_z(rho, phi_lo, phi_hi)

    def _solve_periodic(self, rho: np.ndarray) -> np.ndarray:
        # The k=0 eigenvalue is exactly zero: a net charge has NO periodic
        # solution.  Refuse rather than return a wrong answer (Known-bad).
        cell = self.h ** 3
        net = float(np.sum(rho)) * cell
        scale = float(np.sum(np.abs(rho))) * cell
        if scale > 0.0 and abs(net) > self.neutrality_tol * scale:
            raise ValueError(
                f"fully periodic Poisson needs a net-neutral rho: net charge "
                f"{net:.3e} C is {abs(net) / scale:.3e} of the total |charge| "
                f"(tol {self.neutrality_tol:.1e}).  Add a neutralising "
                "background, or use the Dirichlet-z boundary mode.")
        rho_h = np.fft.fftn(rho)
        lam = self._lam_xyz.copy()
        lam[0, 0, 0] = 1.0                      # placeholder; mode zeroed below
        phi_h = -rho_h / (EPS0 * lam)
        phi_h[0, 0, 0] = 0.0                    # gauge: <phi> = 0
        return np.real(np.fft.ifftn(phi_h))

    def _solve_dirichlet_z(self, rho, phi_lo, phi_hi) -> np.ndarray:
        nz = self.nz
        lo = np.broadcast_to(np.asarray(phi_lo, dtype=float),
                             (self.nx, self.ny))
        hi = np.broadcast_to(np.asarray(phi_hi, dtype=float),
                             (self.nx, self.ny))
        # transform the transverse plane: rho and both electrode profiles
        rho_h = np.fft.fft2(rho, axes=(0, 1))
        lo_h = np.fft.fft2(lo)
        hi_h = np.fft.fft2(hi)

        # unknowns are the interior nodes k = 1 .. nz-2
        b = -rho_h[:, :, 1:nz - 1] / EPS0
        b[:, :, 0] -= self._off * lo_h
        b[:, :, -1] -= self._off * hi_h
        phi_int = _thomas_batch(self._diag, self._off, b)

        phi_h = np.empty_like(rho_h)
        phi_h[:, :, 0] = lo_h
        phi_h[:, :, nz - 1] = hi_h
        phi_h[:, :, 1:nz - 1] = phi_int
        return np.real(np.fft.ifft2(phi_h, axes=(0, 1)))

    # ------------------------------------------------------------- residual
    def residual(self, phi: np.ndarray, rho: np.ndarray) -> float:
        """max |lap(phi) + rho/eps0| over the INTERIOR, normalised.

        An independent check of a returned solution that does not reuse the
        solver: it re-applies the stencil directly.
        """
        lap = laplacian(phi, self.h, self.bc)
        want = -np.asarray(rho, dtype=float) / EPS0
        if self.bc == PERIODIC:
            r, w = lap, want
        else:
            r, w = lap[:, :, 1:-1], want[:, :, 1:-1]
        den = float(np.max(np.abs(w))) + 1e-300
        return float(np.max(np.abs(r - w))) / den


def laplacian(phi: np.ndarray, h: float, bc: str = PERIODIC_XY_DIRICHLET_Z
              ) -> np.ndarray:
    """The same 7-point stencil the solver inverts, applied directly.

    In the Dirichlet mode the z-boundary planes are returned as zero (they are
    not interior equations); callers slice them off.
    """
    phi = np.asarray(phi, dtype=float)
    lap = (np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0)
           + np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1)
           - 4.0 * phi) / (h * h)
    if bc == PERIODIC:
        lap += (np.roll(phi, 1, axis=2) + np.roll(phi, -1, axis=2)
                - 2.0 * phi) / (h * h)
        return lap
    out = np.zeros_like(phi)
    out[:, :, 1:-1] = (lap[:, :, 1:-1]
                       + (phi[:, :, :-2] + phi[:, :, 2:]
                          - 2.0 * phi[:, :, 1:-1]) / (h * h))
    return out


# ------------------------------------------------------------------ oracle
def reference_solve_lu(rho: np.ndarray, h: float,
                       bc: str = PERIODIC_XY_DIRICHLET_Z,
                       phi_lo=0.0, phi_hi=0.0) -> np.ndarray:
    """Assemble the operator explicitly and solve it with a sparse LU.

    THE GATE ORACLE, not a production path: it builds an (N x N) sparse matrix
    for N = nx*ny*nz, so it is only usable on small grids.  The fast separable
    solver must reproduce this to round-off; that equivalence is what licenses
    using the fast path everywhere else.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import splu

    rho = np.asarray(rho, dtype=float)
    nx, ny, nz = rho.shape
    n = nx * ny * nz
    inv_h2 = 1.0 / (h * h)

    def idx(i, j, k):
        return (i * ny + j) * nz + k

    rows, cols, vals = [], [], []
    b = np.zeros(n)
    lo = np.broadcast_to(np.asarray(phi_lo, dtype=float), (nx, ny))
    hi = np.broadcast_to(np.asarray(phi_hi, dtype=float), (nx, ny))

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                r = idx(i, j, k)
                dirichlet_node = (bc != PERIODIC and (k == 0 or k == nz - 1))
                if dirichlet_node:
                    # ROW EQUILIBRATION (do not "simplify" this to 1.0).  The
                    # interior rows carry entries ~1/h^2; at h = 1e-3 that is
                    # 1e6, so an unscaled boundary row (phi = value, entry 1)
                    # leaves the matrix badly row-scaled and the sparse LU
                    # loses ~6 digits -- which reads as the FAST solver being
                    # wrong when it is the oracle that is noisy.  Scaling the
                    # row by 1/h^2 is the same equation and restores
                    # round-off agreement.  Measured 2026-08-03: unscaled
                    # 2.4e-10 vs scaled 1e-16.
                    rows.append(r); cols.append(r); vals.append(inv_h2)
                    b[r] = inv_h2 * (lo[i, j] if k == 0 else hi[i, j])
                    continue
                diag = -6.0 * inv_h2 if bc == PERIODIC else -6.0 * inv_h2
                rows.append(r); cols.append(r); vals.append(diag)
                for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                   (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                    ii, jj, kk = (i + di) % nx, (j + dj) % ny, k + dk
                    if bc == PERIODIC:
                        kk %= nz
                    rows.append(r); cols.append(idx(ii, jj, kk))
                    vals.append(inv_h2)
                b[r] = -rho[i, j, k] / EPS0

    A = csr_matrix((vals, (rows, cols)), shape=(n, n))
    if bc == PERIODIC:
        # pin the singular constant mode so the LU is defined, then re-gauge.
        # Same row-equilibration reasoning as the Dirichlet rows above.
        A = A.tolil()
        A[0, :] = 0.0
        A[0, 0] = inv_h2
        b[0] = 0.0
        A = A.tocsc()
    phi = splu(A.tocsc()).solve(b).reshape(nx, ny, nz)
    if bc == PERIODIC:
        phi -= phi.mean()
    return phi
