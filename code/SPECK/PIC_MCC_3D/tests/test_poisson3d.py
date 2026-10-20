"""Gates for MISSING MODULE 1 -- the 3-D electrostatic Poisson solve.

Every gate is falsifiable against something the solver cannot influence: an
explicitly assembled sparse-LU factorisation of the same operator, a
manufactured solution with a known analytic source, or a configuration whose
exact answer is a polynomial the stencil represents exactly.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.poisson3d import (PERIODIC, PERIODIC_XY_DIRICHLET_Z,
                                          Poisson3D, laplacian,
                                          reference_solve_lu)
from PIC_MCC_3D.MODULES.units import EPS0


def _rel(a, b):
    return float(np.max(np.abs(a - b))) / (float(np.max(np.abs(b))) + 1e-300)


# ------------------------------------------------------- the LU oracle gate
@pytest.mark.parametrize("h", [1.0, 1e-3, 5e-5])
def test_matches_lu_oracle_dirichlet_z(h):
    """The fast separable solve must reproduce an explicit LU to round-off.

    Run across three decades of h: the equilibration of the oracle's boundary
    rows (see poisson3d.reference_solve_lu) is what makes this hold at small
    h, and an unequilibrated oracle regresses to ~1e-10 here.
    """
    rng = np.random.default_rng(3)
    nx, ny, nz = 8, 6, 10
    rho = rng.normal(0.0, 1e-6, (nx, ny, nz))
    lo, hi = 0.3, -1.7
    fast = Poisson3D(nx, ny, nz, h, bc=PERIODIC_XY_DIRICHLET_Z).solve(rho, lo, hi)
    ref = reference_solve_lu(rho, h, bc=PERIODIC_XY_DIRICHLET_Z,
                             phi_lo=lo, phi_hi=hi)
    assert _rel(fast, ref) < 1e-12


@pytest.mark.parametrize("h", [1.0, 1e-3])
def test_matches_lu_oracle_periodic(h):
    rng = np.random.default_rng(5)
    nx, ny, nz = 8, 6, 10
    rho = rng.normal(0.0, 1e-6, (nx, ny, nz))
    rho -= rho.mean()
    fast = Poisson3D(nx, ny, nz, h, bc=PERIODIC).solve(rho)
    ref = reference_solve_lu(rho, h, bc=PERIODIC)
    assert _rel(fast, ref) < 1e-12


def test_boundary_planes_are_honoured_exactly():
    """Electrode planes carry their prescribed potential, not an approximation."""
    rng = np.random.default_rng(11)
    nx, ny, nz, h = 8, 8, 12, 1e-3
    rho = rng.normal(0.0, 1e-6, (nx, ny, nz))
    lo = rng.normal(0.0, 1.0, (nx, ny))          # a PROFILE, not a scalar
    hi = 2.5
    phi = Poisson3D(nx, ny, nz, h).solve(rho, lo, hi)
    assert np.max(np.abs(phi[:, :, 0] - lo)) < 1e-12
    assert np.max(np.abs(phi[:, :, -1] - hi)) < 1e-12


# ------------------------------------------------------------- MMS gates
def _mms_error(nx, bc):
    """L2 error of the solve against a manufactured analytic solution."""
    L = 0.02
    h = L / nx
    ax = np.arange(nx) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    if bc == PERIODIC:
        kx = ky = kz = 2.0 * np.pi / L
        phi_exact = np.sin(kx * X) * np.sin(ky * Y) * np.sin(kz * Z)
        lo = hi = 0.0
    else:
        # z runs 0..(nz-1)h = L - h; use a sine that vanishes on both planes
        Lz = (nx - 1) * h
        kx = ky = 2.0 * np.pi / L
        kz = np.pi / Lz
        phi_exact = np.sin(kx * X) * np.cos(ky * Y) * np.sin(kz * Z)
        lo = hi = 0.0
    c = kx ** 2 + ky ** 2 + kz ** 2
    rho = EPS0 * c * phi_exact                 # rho = -eps0 lap(phi)
    phi = Poisson3D(nx, nx, nx, h, bc=bc).solve(rho, lo, hi)
    if bc == PERIODIC:
        phi = phi - phi.mean()
        phi_exact = phi_exact - phi_exact.mean()
    err = np.sqrt(np.mean((phi - phi_exact) ** 2))
    return h, float(err)


@pytest.mark.parametrize("bc", [PERIODIC_XY_DIRICHLET_Z, PERIODIC])
def test_mms_second_order_accuracy(bc):
    """Method of manufactured solutions: the scheme must converge at O(h^2).

    This is rung 1 of the verification ladder (Riva 2017's code-verification
    step, ../PIC_MCC_3D.txt section 5).  A solver can match an LU oracle
    perfectly and still discretise the wrong equation; only a convergence
    order test catches that.
    """
    hs, errs = [], []
    for nx in (8, 12, 16, 24):
        h, e = _mms_error(nx, bc)
        hs.append(h)
        errs.append(e)
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    assert 1.9 < slope < 2.1, f"observed order {slope:.3f}, expected ~2"


# --------------------------------------------- exact-polynomial configurations
def test_vacuum_capacitor_is_a_linear_ramp():
    """No charge, unequal plates -> phi linear in z, exact to round-off."""
    nx, ny, nz, h = 6, 6, 17, 1e-3
    V = 12.0
    phi = Poisson3D(nx, ny, nz, h).solve(np.zeros((nx, ny, nz)), 0.0, V)
    z = np.arange(nz) * h
    want = np.broadcast_to(V * z / ((nz - 1) * h), (nx, ny, nz))
    assert _rel(phi, want) < 1e-12


def test_parallel_plate_uniform_charge_is_exact():
    """Uniform rho between grounded plates -> a parabola the stencil is EXACT on.

    The 3-point second difference reproduces a quadratic with no truncation
    error, so this is a round-off gate, not a convergence gate.
    """
    nx, ny, nz, h = 6, 6, 21, 1e-3
    rho0 = 1e-6
    d = (nz - 1) * h
    rho = np.full((nx, ny, nz), rho0)
    phi = Poisson3D(nx, ny, nz, h).solve(rho, 0.0, 0.0)
    z = np.arange(nz) * h
    want = np.broadcast_to(rho0 * z * (d - z) / (2.0 * EPS0), (nx, ny, nz))
    assert _rel(phi, want) < 1e-12


# --------------------------------------------------------- discrete identity
def test_solution_satisfies_the_discrete_operator():
    """lap(phi) = -rho/eps0 on the 7-point stencil -- the discrete Gauss law.

    Applied by re-running the stencil directly, so it does not reuse any of
    the solver's own machinery.
    """
    rng = np.random.default_rng(17)
    nx, ny, nz, h = 10, 8, 12, 1e-3
    rho = rng.normal(0.0, 1e-6, (nx, ny, nz))
    s = Poisson3D(nx, ny, nz, h)
    phi = s.solve(rho, 0.4, -0.9)
    assert s.residual(phi, rho) < 1e-10
    # and independently, for the periodic mode
    rho -= rho.mean()
    sp = Poisson3D(nx, ny, nz, h, bc=PERIODIC)
    phip = sp.solve(rho)
    lap = laplacian(phip, h, PERIODIC)
    assert _rel(lap, -rho / EPS0) < 1e-10


# ------------------------------------------------------------------- guards
def test_periodic_rejects_net_charge():
    """The k=0 mode is singular: a net charge has NO periodic solution."""
    nx, h = 8, 1e-3
    rho = np.full((nx, nx, nx), 1e-6)          # emphatically not neutral
    with pytest.raises(ValueError, match="net-neutral"):
        Poisson3D(nx, nx, nx, h, bc=PERIODIC).solve(rho)


def test_periodic_accepts_neutral_charge():
    rng = np.random.default_rng(2)
    nx, h = 8, 1e-3
    rho = rng.normal(0.0, 1e-6, (nx, nx, nx))
    rho -= rho.mean()
    phi = Poisson3D(nx, nx, nx, h, bc=PERIODIC).solve(rho)
    assert np.isfinite(phi).all()
    assert abs(float(phi.mean())) < 1e-12      # gauge: <phi> = 0


def test_constructor_rejects_bad_inputs():
    with pytest.raises(ValueError, match="unknown bc"):
        Poisson3D(8, 8, 8, 1e-3, bc="nonsense")
    with pytest.raises(ValueError, match="at least 3 nodes"):
        Poisson3D(2, 8, 8, 1e-3)
    with pytest.raises(ValueError, match="must be positive"):
        Poisson3D(8, 8, 8, 0.0)


def test_solve_rejects_wrong_shape():
    s = Poisson3D(8, 8, 8, 1e-3)
    with pytest.raises(ValueError, match="shape"):
        s.solve(np.zeros((8, 8, 7)))
