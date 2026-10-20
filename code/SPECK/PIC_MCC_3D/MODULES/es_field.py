"""Electrostatic field <-> particle coupling: deposit, gradient, gather.

Three operations, and the ONE invariant that ties them together:

    deposit_rho     particles -> rho on the nodes        (trilinear / CIC)
    efield_from_phi phi -> E = -grad phi on the nodes     (finite difference)
    gather_E        E on the nodes -> E at the particles  (trilinear / CIC)

THE INVARIANT: ``deposit_rho`` and ``gather_E`` use the SAME trilinear shape
function.  That is the momentum-conserving pairing -- the internal forces of a
closed system cancel to round-off, and an isolated particle exerts no net
force on itself.  Mixing shapes (a CIC deposit with a nearest-grid-point
gather, say) breaks both properties and accelerates particles out of nothing.
This is on the Known-bad list; a gate measures the residual self-force and the
sum of internal forces.

Geometry conventions (shared with :mod:`poisson3d`): fields live on NODES,
shape ``(nx, ny, nz)``, uniform spacing ``h``.  x and y are periodic, so node
``nx`` is node ``0`` and a particle's x wraps into ``[0, nx h)``.  z is either
periodic too, or bounded by electrode planes at ``k = 0`` and ``k = nz-1``, in
which case particles live in ``[0, (nz-1) h]`` and leaving that range is a
wall event the caller handles.
"""

from __future__ import annotations

import numpy as np

from .poisson3d import PERIODIC, PERIODIC_XY_DIRICHLET_Z


def _weights_periodic(coord: np.ndarray, h: float, n: int):
    """Base index, upper index and the linear weight, periodic axis."""
    s = coord / h
    i0 = np.floor(s).astype(np.intp)
    w = s - i0
    i0 = np.mod(i0, n)
    i1 = np.mod(i0 + 1, n)
    return i0, i1, w


def _weights_bounded(coord: np.ndarray, h: float, n: int):
    """Base index, upper index and the linear weight, bounded (z) axis.

    The base cell is clamped to ``[0, n-2]`` so both node indices are always
    valid; a particle exactly on the far plane lands with weight 1 on it.
    """
    s = coord / h
    i0 = np.floor(s).astype(np.intp)
    i0 = np.clip(i0, 0, n - 2)
    w = s - i0
    return i0, i0 + 1, w


def _cic_indices(shape, h, x, y, z, bc):
    nx, ny, nz = shape
    i0, i1, wx = _weights_periodic(np.asarray(x, dtype=float), h, nx)
    j0, j1, wy = _weights_periodic(np.asarray(y, dtype=float), h, ny)
    if bc == PERIODIC:
        k0, k1, wz = _weights_periodic(np.asarray(z, dtype=float), h, nz)
    else:
        k0, k1, wz = _weights_bounded(np.asarray(z, dtype=float), h, nz)
    return (i0, i1, wx), (j0, j1, wy), (k0, k1, wz)


def _corner_terms(ix, iy, iz):
    """The eight (i, j, k, weight-factors) corners of the trilinear stencil."""
    (i0, i1, wx), (j0, j1, wy), (k0, k1, wz) = ix, iy, iz
    return (
        (i0, j0, k0, (1 - wx) * (1 - wy) * (1 - wz)),
        (i1, j0, k0, wx * (1 - wy) * (1 - wz)),
        (i0, j1, k0, (1 - wx) * wy * (1 - wz)),
        (i1, j1, k0, wx * wy * (1 - wz)),
        (i0, j0, k1, (1 - wx) * (1 - wy) * wz),
        (i1, j0, k1, wx * (1 - wy) * wz),
        (i0, j1, k1, (1 - wx) * wy * wz),
        (i1, j1, k1, wx * wy * wz),
    )


def deposit_rho(shape, h: float, x, y, z, charge_weight,
                bc: str = PERIODIC_XY_DIRICHLET_Z, out=None) -> np.ndarray:
    """Trilinear (CIC) deposit of ``charge_weight`` [C] onto the nodes.

    ``charge_weight`` is q*w per macro-particle (signed).  The result is a
    charge DENSITY [C/m^3]: the deposited charge is divided by the cell
    volume h^3.  Total charge is conserved exactly (the eight corner weights
    sum to 1 identically), which a gate checks.
    """
    nx, ny, nz = shape
    rho = np.zeros(shape) if out is None else out
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return rho
    q = np.asarray(charge_weight, dtype=float)
    q = np.broadcast_to(q, x.shape)
    ix, iy, iz = _cic_indices(shape, h, x, y, z, bc)
    inv_vol = 1.0 / (h * h * h)
    # bincount over FLAT indices, not np.add.at.  np.add.at is an unbuffered
    # ufunc call and runs ~20-50x slower here; at 1e5 particles x 8 corners x
    # 1e3 steps that is the difference between a gate that runs in seconds and
    # one that times out.  Concatenating the eight corners into a single
    # bincount also beats eight separate calls.
    flat = np.empty(8 * x.size, dtype=np.intp)
    vals = np.empty(8 * x.size, dtype=float)
    for c, (i, j, k, wgt) in enumerate(_corner_terms(ix, iy, iz)):
        s = slice(c * x.size, (c + 1) * x.size)
        flat[s] = (i * ny + j) * nz + k
        vals[s] = q * wgt * inv_vol
    rho += np.bincount(flat, weights=vals,
                       minlength=nx * ny * nz).reshape(shape)
    return rho


def efield_from_phi(phi: np.ndarray, h: float,
                    bc: str = PERIODIC_XY_DIRICHLET_Z):
    """E = -grad phi on the nodes; returns (Ex, Ey, Ez).

    Centred differences in the periodic directions.  In the bounded direction
    the interior is centred and the two electrode planes use one-sided
    (first-order) differences -- the field AT a perfectly conducting plane is
    not a centred quantity, and particles are absorbed there anyway.
    """
    phi = np.asarray(phi, dtype=float)
    inv_2h = 1.0 / (2.0 * h)
    Ex = -(np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) * inv_2h
    Ey = -(np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) * inv_2h
    if bc == PERIODIC:
        Ez = -(np.roll(phi, -1, axis=2) - np.roll(phi, 1, axis=2)) * inv_2h
        return Ex, Ey, Ez
    Ez = np.empty_like(phi)
    Ez[:, :, 1:-1] = -(phi[:, :, 2:] - phi[:, :, :-2]) * inv_2h
    Ez[:, :, 0] = -(phi[:, :, 1] - phi[:, :, 0]) / h
    Ez[:, :, -1] = -(phi[:, :, -1] - phi[:, :, -2]) / h
    return Ex, Ey, Ez


def gather_E(Ex, Ey, Ez, h: float, x, y, z,
             bc: str = PERIODIC_XY_DIRICHLET_Z):
    """Trilinear gather of the node fields to particle positions.

    The SAME shape function as :func:`deposit_rho` -- see the module
    docstring.  Exact for a field that is linear in space (gated).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    ex = np.zeros(n)
    ey = np.zeros(n)
    ez = np.zeros(n)
    if n == 0:
        return ex, ey, ez
    ix, iy, iz = _cic_indices(Ex.shape, h, x, y, z, bc)
    for i, j, k, wgt in _corner_terms(ix, iy, iz):
        ex += Ex[i, j, k] * wgt
        ey += Ey[i, j, k] * wgt
        ez += Ez[i, j, k] * wgt
    return ex, ey, ez


def field_energy_J(Ex, Ey, Ez, h: float) -> float:
    """(eps0/2) integral |E|^2 dV over the grid [J]."""
    from .units import EPS0
    e2 = Ex * Ex + Ey * Ey + Ez * Ez
    return float(0.5 * EPS0 * np.sum(e2) * h ** 3)
