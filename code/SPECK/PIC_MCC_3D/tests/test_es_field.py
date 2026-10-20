"""Gates for the deposit / gradient / gather trio.

The load-bearing property is the one in the module docstring: deposit and
gather share a shape function, which makes the internal forces of a closed
system cancel and the self-force vanish.  Both are measured here, because a
scheme that quietly violates them accelerates particles out of nothing.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.es_field import (deposit_rho, efield_from_phi,
                                         field_energy_J, gather_E)
from PIC_MCC_3D.MODULES.poisson3d import PERIODIC, PERIODIC_XY_DIRICHLET_Z, Poisson3D
from PIC_MCC_3D.MODULES.units import E_CHARGE, EPS0


def test_deposit_conserves_total_charge():
    """The eight corner weights sum to 1, so charge is conserved exactly."""
    rng = np.random.default_rng(1)
    shape, h = (12, 10, 14), 1e-3
    n = 20000
    x = rng.random(n) * shape[0] * h
    y = rng.random(n) * shape[1] * h
    z = rng.random(n) * (shape[2] - 1) * h
    q = rng.normal(0.0, 1e-15, n)
    rho = deposit_rho(shape, h, x, y, z, q, PERIODIC_XY_DIRICHLET_Z)
    deposited = float(np.sum(rho)) * h ** 3
    assert abs(deposited - q.sum()) / np.abs(q).sum() < 1e-14


def test_deposit_of_a_uniform_load_is_uniform():
    """A regular lattice load must produce a flat density (the quiet start)."""
    shape, h = (8, 8, 8), 1e-3
    per = 16
    g = (np.arange(per) + 0.5) / per
    gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
    L = shape[0] * h
    x, y, z = (gx.ravel() * L, gy.ravel() * L, gz.ravel() * L)
    q = np.full(x.size, -E_CHARGE)
    rho = deposit_rho(shape, h, x, y, z, q, PERIODIC)
    assert np.std(rho) / abs(np.mean(rho)) < 1e-12


def test_gather_is_exact_for_a_linear_field():
    """Trilinear interpolation reproduces a linear field exactly."""
    rng = np.random.default_rng(4)
    shape, h = (10, 10, 12), 1e-3
    ii, jj, kk = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    # linear only along z, so the periodic x/y wrap seam cannot pollute it
    Ez = 3.0 + 7.0 * kk * h
    Ex = np.zeros(shape)
    Ey = np.zeros(shape)
    n = 5000
    x = rng.random(n) * shape[0] * h
    y = rng.random(n) * shape[1] * h
    z = rng.random(n) * (shape[2] - 1) * h
    _, _, ez = gather_E(Ex, Ey, Ez, h, x, y, z, PERIODIC_XY_DIRICHLET_Z)
    assert np.max(np.abs(ez - (3.0 + 7.0 * z))) < 1e-12


def test_efield_of_a_linear_potential_is_constant():
    shape, h = (8, 8, 10), 1e-3
    _, _, kk = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    phi = -5.0 * kk * h
    _, _, Ez = efield_from_phi(phi, h, PERIODIC_XY_DIRICHLET_Z)
    assert np.max(np.abs(Ez - 5.0)) < 1e-9


def test_internal_forces_cancel():
    """Sum of forces over a closed system = 0 (momentum conservation).

    Guaranteed by deposit and gather sharing a shape function together with a
    symmetric field solve.  Measured against the largest single-particle
    force, so the bound is meaningful rather than absolute.
    """
    rng = np.random.default_rng(7)
    nx, h = 16, 1e-3
    shape, L = (nx, nx, nx), nx * h
    n = 500
    x, y, z = (rng.random(n) * L, rng.random(n) * L, rng.random(n) * L)
    q = np.full(n, -E_CHARGE * 1e6)
    rho = deposit_rho(shape, h, x, y, z, q, PERIODIC)
    rho -= rho.mean()                          # neutralising background
    phi = Poisson3D(nx, nx, nx, h, bc=PERIODIC).solve(rho)
    E = efield_from_phi(phi, h, PERIODIC)
    ex, ey, ez = gather_E(*E, h, x, y, z, PERIODIC)
    forces = q[:, None] * np.column_stack((ex, ey, ez))
    scale = float(np.max(np.abs(forces)))
    assert float(np.max(np.abs(forces.sum(axis=0)))) / scale < 1e-12


def test_self_force_is_zero():
    """One particle in a neutralising background feels no force from itself."""
    nx, h = 16, 1e-3
    shape, L = (nx, nx, nx), nx * h
    x = np.array([0.37 * L])
    y = np.array([0.51 * L])
    z = np.array([0.62 * L])
    q = np.array([-E_CHARGE * 1e6])
    rho = deposit_rho(shape, h, x, y, z, q, PERIODIC)
    rho -= rho.mean()
    phi = Poisson3D(nx, nx, nx, h, bc=PERIODIC).solve(rho)
    E = efield_from_phi(phi, h, PERIODIC)
    ex, ey, ez = gather_E(*E, h, x, y, z, PERIODIC)
    self_force = abs(float(q[0] * ex[0])) + abs(float(q[0] * ey[0])) \
        + abs(float(q[0] * ez[0]))
    # reference scale: the force this charge would feel from a like charge one
    # cell away, i.e. the smallest force the scheme is meant to resolve
    ref = float(q[0] ** 2 / (4.0 * np.pi * EPS0 * h ** 2))
    assert self_force / ref < 1e-12


def test_field_energy_of_a_uniform_field():
    shape, h = (8, 8, 8), 1e-3
    Ex = np.full(shape, 100.0)
    Ey = np.zeros(shape)
    Ez = np.zeros(shape)
    want = 0.5 * EPS0 * 100.0 ** 2 * (8 * 8 * 8) * h ** 3
    assert abs(field_energy_J(Ex, Ey, Ez, h) - want) / want < 1e-12


def test_empty_population_is_a_no_op():
    shape, h = (6, 6, 6), 1e-3
    empty = np.zeros(0)
    rho = deposit_rho(shape, h, empty, empty, empty, empty, PERIODIC)
    assert np.all(rho == 0.0)
    ex, ey, ez = gather_E(rho, rho, rho, h, empty, empty, empty, PERIODIC)
    assert ex.size == ey.size == ez.size == 0


@pytest.mark.parametrize("bc", [PERIODIC, PERIODIC_XY_DIRICHLET_Z])
def test_deposit_keeps_particles_on_the_grid(bc):
    """Indices stay in range for particles at and beyond the domain edges."""
    shape, h = (6, 6, 8), 1e-3
    edge = np.array([0.0, shape[0] * h, -1e-9, shape[0] * h + 1e-9])
    z = np.array([0.0, (shape[2] - 1) * h, 1e-9, (shape[2] - 1) * h - 1e-9])
    q = np.full(4, -E_CHARGE)
    rho = deposit_rho(shape, h, edge, edge, z, q, bc)
    assert np.isfinite(rho).all()
