"""Scharfetter–Gummel drift-diffusion flux on the 1D mesh.

The charged-species continuity equation (Chen et al. 2025, Eq. 1)

    ∂n/∂t + ∂ₓΓ = S,   with the drift-diffusion flux (Eq. 2)   Γ = −D ∂ₓn + v n,

is discretised finite-volume on the cell-centred mesh. The face fluxes use the
**Scharfetter–Gummel exponential-fitting** scheme — the standard for drift-diffusion because
it is stable and non-oscillatory when drift dominates diffusion (the high-field excimer
regime), degrading gracefully to upwinding as the cell Péclet number ``Δ = v·dx/D`` grows and
to centred differencing as ``Δ → 0``. This is the same exponential scheme Hagelaar & Pitchford
(2005) use inside the Boltzmann solve itself.

For a face between cells ``i`` (left) and ``i+1`` (right), with face drift velocity ``v`` and
diffusion ``D``:

    J_{i+½} = (D/dx)·[ B(−Δ)·n_i − B(Δ)·n_{i+1} ],   Δ = v·dx/D,   B(t) = t/(eᵗ−1).

``v`` is the signed drift velocity ``v = z·μ·E`` (electrons z=−1 drift against E). The module
provides the interior face fluxes and the flux-divergence update; the caller supplies the two
boundary-face fluxes (electrode BCs), which are physics, not geometry.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# Floor on |D| at a face so Δ = v·dx/D stays finite (pure-drift limit handled by SG anyway).
_D_FLOOR = 1.0e-30


def bernoulli(x: NDArray[np.floating]) -> NDArray[np.float64]:
    """Bernoulli function B(x) = x/(eˣ−1), numerically stable (B(0)=1, B(x→+∞)→0)."""
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    small = np.abs(x) < 1.0e-10
    # Series near 0: B(x) ≈ 1 − x/2 (+ x²/12); avoids 0/0.
    out[small] = 1.0 - 0.5 * x[small]
    xl = x[~small]
    # expm1 is accurate; xl/inf → 0 for large +x, xl/(−1) → −x for large −x (correct linear growth).
    with np.errstate(over="ignore", invalid="ignore"):
        out[~small] = xl / np.expm1(xl)
    return out


def sg_interior_flux(
    n: NDArray[np.float64],
    v_face: NDArray[np.float64],
    d_face: NDArray[np.float64],
    dx: float,
) -> NDArray[np.float64]:
    """Scharfetter–Gummel fluxes at the ``N−1`` interior faces (signed, +x direction).

    ``n`` is the ``(N,)`` cell-centred density; ``v_face``/``d_face`` are the ``(N−1,)`` face
    drift velocity (m/s) and diffusion (m²/s). Returns the ``(N−1,)`` face fluxes (m⁻²s⁻¹).
    """
    n = np.asarray(n, dtype=np.float64)
    v = np.asarray(v_face, dtype=np.float64)
    d = np.asarray(d_face, dtype=np.float64)
    d_safe = np.where(np.abs(d) < _D_FLOOR, _D_FLOOR, d)
    delta = v * dx / d_safe
    return (d_safe / dx) * (bernoulli(-delta) * n[:-1] - bernoulli(delta) * n[1:])


def flux_divergence(
    interior_flux: NDArray[np.float64],
    flux_left: float,
    flux_right: float,
    dx: float,
) -> NDArray[np.float64]:
    """Transport rate ∂n/∂t|_transport = −(J_{i+½} − J_{i−½})/dx for all ``N`` cells.

    ``interior_flux`` is the ``(N−1,)`` SG face fluxes; ``flux_left`` / ``flux_right`` are the
    +x-directed fluxes at the domain boundary faces (x=0 and x=L) — the electrode BCs.
    """
    interior = np.asarray(interior_flux, dtype=np.float64)
    n_cells = interior.size + 1
    faces = np.empty(n_cells + 1, dtype=np.float64)  # flux at every face, incl. the two walls
    faces[0] = flux_left
    faces[-1] = flux_right
    faces[1:-1] = interior
    return -(faces[1:] - faces[:-1]) / dx


def face_average(cell_values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Arithmetic-mean a cell-centred field onto the ``N−1`` interior faces."""
    a = np.asarray(cell_values, dtype=np.float64)
    return 0.5 * (a[:-1] + a[1:])
