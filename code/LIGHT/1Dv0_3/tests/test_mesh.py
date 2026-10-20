"""1Dv0_3's mesh (mesh.py): OSIRIS's equal-load rule on the cell edges.

NumPy only, so these run under LIGHT's usual python3 as well as with
Fluid1D's Python."""

import importlib
import math

import numpy as np
import pytest

M = importlib.import_module("1Dv0_3.mesh")

L = 1e-2                                  # the reference point's 1-cm gap


def test_uniform_mesh_is_fluid1d_s():
    """Mesh1D's arithmetic: dx = L/N, edges arange(N+1) dx, centres
    (i + 1/2) dx, and face weights of exactly 1/2 (so a face average is
    Fluid1D's 0.5 (a + b) to the bit)."""
    for n in (8, 100, 128, 140):
        m = M.uniform(n, L)
        dx = L / n
        assert np.array_equal(m.widths, np.full(n, dx))
        assert np.array_equal(m.edges, np.arange(n + 1) * dx)
        assert np.array_equal(m.centres, (np.arange(n) + 0.5) * dx)
        wl, wr = m.weights
        assert np.all(wl == 0.5) and np.all(wr == 0.5)
        assert np.array_equal(m.spacing, np.full(n - 1, dx))
    with pytest.raises(ValueError):
        M.Mesh(np.array([0.0, 1e-3, 3e-3]))      # not arange(N+1) dx


def test_graded_cells_grow_from_the_wall_cell():
    """The first cell is the wall cell asked for, the cells grow by one
    ratio e^s up to the middle, the two walls are mirror images (to the
    rounding of L - x, since the far wall's edges are L less a distance),
    and N cells fill the gap."""
    for n, wall in ((64, 2e-6), (140, 0.5e-6), (99, 1e-6)):
        m = M.graded(n, L, wall)
        h = m.widths
        assert m.n_cells == n and m.edges[0] == 0.0 and m.edges[-1] == L
        assert h[0] == pytest.approx(wall, rel=1e-12, abs=0)
        _, s, _ = m.load
        half = h[: n // 2 - 1]
        np.testing.assert_allclose(half[1:] / half[:-1], math.exp(s),
                                   rtol=1e-9, atol=0)
        np.testing.assert_allclose(h, h[::-1], rtol=0,
                                   atol=4 * np.spacing(L))
        assert m.ratio <= M.RATIO_MAX


def test_face_weights_are_the_half_cells_in_series():
    """A face average weighs each side by its half-width: the conductance
    of the two half-cells in series, d / (h_l/2a_l + h_r/2a_r), is that
    average to first order in the contrast a_r/a_l - 1 (the field's
    coefficient and the diffusion are conductances)."""
    m = M.graded(64, L, 2e-6)
    h = m.widths
    wl, wr = m.weights
    np.testing.assert_allclose(wl, h[:-1] / (h[:-1] + h[1:]), rtol=1e-15,
                               atol=0)
    np.testing.assert_allclose(wl + wr, 1.0, rtol=1e-15, atol=0)
    a = 1.0 + 1e-4 * np.random.default_rng(0).standard_normal(64)
    series = m.spacing / (h[:-1] / (2 * a[:-1]) + h[1:] / (2 * a[1:]))
    mean = wl * a[:-1] + wr * a[1:]
    assert np.max(np.abs(series / mean - 1.0)) < 1e-7


def test_every_cell_holds_an_equal_share_of_the_load():
    """OSIRIS's rule: the edges sit at the k/N quantiles of the load's
    integral, so every cell holds 1/N of it."""
    m = M.graded(128, L, 1e-6)
    h0, s, h_max = m.load
    ell, _ = M.wall_load(h0, s, h_max, L / 2)
    x = m.edges
    cum = np.where(x <= L / 2, ell(x), 2 * ell(L / 2) - ell(L - x))
    shares = np.diff(cum)
    np.testing.assert_allclose(shares, shares.mean(), rtol=1e-10, atol=0)
    assert shares.mean() == pytest.approx(1.0, rel=1e-12, abs=0)


def test_osiris_recursive_bisection_gives_the_same_edges():
    """rec_tree_dist (OSIRIS, os-grid-parallel.f90:236) cuts a range of n
    parts where the load reaches its left end's plus nl/n of its own, nl =
    n // 2, then cuts each side. On a continuous load that recursion lands
    on the quantiles, which mesh.py places directly."""
    m = M.graded(77, L, 1e-6)
    h0, s, h_max = m.load
    ell, inverse = M.wall_load(h0, s, h_max, L / 2)
    total = 2 * float(ell(L / 2))

    def cum(x):
        return float(ell(x)) if x <= L / 2 else total - float(ell(L - x))

    def place(v):                        # the point where the load is v
        return float(inverse(v)) if v <= total / 2 \
            else L - float(inverse(total - v))

    edges = {0: 0.0, 77: L}

    def split(lo, hi):
        n = hi - lo
        if n < 2:
            return
        nl = n // 2
        a, b = cum(edges[lo]), cum(edges[hi])
        edges[lo + nl] = place(a + (b - a) * nl / n)
        split(lo, lo + nl)
        split(lo + nl, hi)
    split(0, 77)
    rec = np.array([edges[k] for k in range(78)])
    np.testing.assert_allclose(rec, m.edges, rtol=0, atol=1e-15)


def test_general_rule_on_the_wall_load_gives_the_mesh():
    """equidistribute, the rule for any load density (OSIRIS's
    spatial_loaddensity), reproduces the analytic edges on w = 1/h(d)."""
    m = M.graded(70, L, 1e-6)
    h0, s, _ = m.load
    edges = M.equidistribute(lambda x: 1.0 / (h0 + s * np.minimum(x, L - x)),
                             70, L)
    assert np.max(np.abs(edges - m.edges)) < 1e-4 * m.wall_cell


def test_refine_is_nested_and_rebuilds_from_its_parameters():
    """refine doubles the cells on the same load: every edge stays, the
    neighbour ratio is its square root, and the refined mesh is the graded
    mesh of its own cell count and wall cell."""
    m = M.graded(64, L, 2e-6)
    for levels in (1, 2):
        r = M.refine(m, levels)
        assert r.n_cells == 64 * 2 ** levels
        np.testing.assert_allclose(r.edges[:: 2 ** levels], m.edges,
                                   rtol=0, atol=1e-18)
        assert math.exp(r.load[1]) == pytest.approx(
            math.exp(m.load[1] / 2 ** levels), rel=1e-12, abs=0)
        again = M.graded(r.n_cells, L, r.wall_cell)
        np.testing.assert_allclose(again.edges, r.edges, rtol=0, atol=1e-18)
    u = M.refine(M.uniform(32, L))
    assert u.uniform and u.n_cells == 64


def test_a_cap_holds_the_largest_cell():
    m = M.graded(100, L, 1e-6, max_cell_m=200e-6)
    assert m.widths.max() <= 200e-6 * (1 + 1e-9)
    middle = m.widths[40:60]
    np.testing.assert_allclose(middle, 200e-6, rtol=1e-9, atol=0)
    assert m.widths[0] == pytest.approx(1e-6, rel=1e-12, abs=0)


def test_impossible_meshes_are_refused():
    with pytest.raises(ValueError, match="grow by"):
        M.graded(40, L, 5e-6)                 # neighbours 1.34 apart
    with pytest.raises(ValueError, match="fill"):
        M.graded(100, L, 100e-6)              # 100 x 100 um is the gap
    with pytest.raises(ValueError, match="cannot"):
        M.graded(40, L, 1e-6, max_cell_m=100e-6)
    with pytest.raises(ValueError):
        M.build("stretched", 64, L)
