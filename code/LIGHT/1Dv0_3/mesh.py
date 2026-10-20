"""1Dv0_3's mesh: cells placed by OSIRIS's load-balancing rule.

OSIRIS (``Kinetic_Modeling/osiris``) keeps one cell size, dx = (xmax -
xmin)/nx_p (``source/os-simulation.f03``, line 330), and scales the
computation instead: its load balancing gives each processor an equal
share of a load, the particles in its cells or a load density given as an
expression (``grid{ lb_type, spatial_loaddensity }``). The cuts come from
recursive bisection of the running load total (``rec_tree_dist``,
``source/grid/os-grid-parallel.f90``, line 236): a range holding n parts
is cut where the load reaches its left end's value plus nl/n of its own,
nl = n/2, and each side is cut again. For a continuous load this puts the
k-th cut at the k/N quantile of the total.

1Dv0_3 applies the rule to the cell edges themselves. The N + 1 edges sit
at the quantiles of a load density w(x), so every cell holds an equal
share of the load and the cells are small where the load is dense. The
load is the number of cells wanted per unit length, w = 1/h(d), with

    h(d) = min(h_max, h0 + s d),    d = the distance to the nearer wall,

which makes the cells grow geometrically from each electrode: the k-th is
h_wall r^k with r = e^s and h_wall = h0 (r - 1)/s, until they reach h_max
(no cap by default). The important interactions in 1Dv0_2's gap are at
the walls. The sheath is a few Debye lengths thick (about 1 um in the
pulse at the reference point), and a uniform mesh of 39-um cells does not
resolve it (1Dv0_2's overview, 2026-09-30). OSIRIS sets no limit on how
fast its ranks' widths change, but a mesh needs one, and here neighbours
differ by the ratio r.

``graded(cells, length, wall_cell)`` solves for s so that the load totals
N cells with the first cell h_wall wide. ``refine`` doubles N on the same
load: h0, s and h_max halve, every edge stays, and each cell splits at its
load midpoint, so a refinement study runs on nested meshes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

RATIO_MAX = 1.3            # a mesh grading faster is refused: add cells


@dataclass(frozen=True)
class Mesh:
    """N cells between two walls, edges x_0 = 0 .. x_N = L."""

    edges: np.ndarray
    area_m2: float = 1.0
    kind: str = "uniform"
    # the graded mesh's load (refine() keeps it): h0 [m], s, h_max [m]
    load: tuple = None
    widths: np.ndarray = field(init=False, repr=False)
    centres: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        e = np.asarray(self.edges, dtype=np.float64)
        if e.ndim != 1 or e.size < 3 or not np.all(np.diff(e) > 0.0):
            raise ValueError("a mesh needs at least two cells of positive "
                             "width")
        object.__setattr__(self, "edges", e)
        n = e.size - 1
        if self.kind == "uniform":
            # Fluid1D's arithmetic exactly (Mesh1D): one dx, centres (i+1/2)dx.
            # The edges are arange(N+1) dx, so edges[1] is dx to the bit
            # (edges[-1]/N need not be)
            dx = float(e[1] - e[0])
            if not np.array_equal(e, np.arange(n + 1, dtype=np.float64) * dx):
                raise ValueError("uniform edges must be arange(N+1) dx")
            object.__setattr__(self, "widths", np.full(n, dx))
            object.__setattr__(self, "centres",
                               (np.arange(n, dtype=np.float64) + 0.5) * dx)
        else:
            object.__setattr__(self, "widths", np.diff(e))
            object.__setattr__(self, "centres", 0.5 * (e[:-1] + e[1:]))

    @property
    def n_cells(self) -> int:
        return self.widths.size

    @property
    def cell_centers(self) -> np.ndarray:
        """Fluid1D's name for the centres (its orchestrator reads it)."""
        return self.centres

    @property
    def length(self) -> float:
        return float(self.edges[-1] - self.edges[0])

    @property
    def uniform(self) -> bool:
        return self.kind == "uniform"

    @property
    def spacing(self) -> np.ndarray:
        """The centre-to-centre distance at the N - 1 interior faces."""
        if self.uniform:
            return self.widths[:-1].copy()
        return 0.5 * (self.widths[:-1] + self.widths[1:])

    @property
    def weights(self) -> tuple:
        """(left, right) weights of a face average at the interior faces,
        the shares of the half-widths on either side: the two half-cells'
        conductances in series, to first order in their contrast (the
        field's coefficient and the diffusion are conductances). They are
        exactly 0.5 and 0.5 on a uniform mesh, where w a + w b is
        Fluid1D's 0.5 (a + b) to the bit."""
        h = self.widths
        tot = h[:-1] + h[1:]
        return h[:-1] / tot, h[1:] / tot

    @property
    def volumes(self) -> np.ndarray:
        return self.widths * self.area_m2

    @property
    def ratio(self) -> float:
        """The largest width ratio of two neighbouring cells."""
        h = self.widths
        r = h[1:] / h[:-1]
        return float(np.max(np.maximum(r, 1.0 / r)))

    @property
    def wall_cell(self) -> float:
        return float(self.widths[0])

    def describe(self) -> str:
        h = self.widths * 1e6
        if self.uniform:
            return f"{self.n_cells} uniform cells of {h[0]:.4g} um"
        return (f"{self.n_cells} graded cells, {h[0]:.4g} um at the walls to "
                f"{h.max():.4g} um, neighbours within {self.ratio:.4f}")


def uniform(cells: int, length_m: float, area_m2: float = 1.0) -> Mesh:
    """Fluid1D's mesh: N cells of length/N (Mesh1D's dx, edges and centres
    to the bit)."""
    n = int(cells)
    dx = length_m / n
    return Mesh(np.arange(n + 1, dtype=np.float64) * dx, area_m2)


def build(kind: str, cells: int, length_m: float, area_m2: float = 1.0,
          wall_cell_m: float = None, max_cell_m: float = None) -> Mesh:
    """The mesh a run asks for: 'uniform' (Fluid1D's) or 'graded'."""
    if kind == "uniform":
        return uniform(cells, length_m, area_m2)
    if kind == "graded":
        if wall_cell_m is None or wall_cell_m <= 0.0:
            raise ValueError("a graded mesh needs a positive wall cell")
        return graded(cells, length_m, wall_cell_m, max_cell_m, area_m2)
    raise ValueError(f"mesh is 'uniform' or 'graded', not {kind!r}")


# ---------------------------------------------------------------------------
# The load and its quantiles (OSIRIS's rule)
# ---------------------------------------------------------------------------
def wall_load(h0: float, s: float, h_max, half: float):
    """The wall load's integral from a wall, ell(d) = int_0^d dx/h(x), with
    h = min(h_max, h0 + s x), and its inverse, for d in [0, half]."""
    d_cap = math.inf if h_max is None else max((h_max - h0) / s, 0.0)
    ell_cap = math.log1p(s * d_cap / h0) / s if math.isfinite(d_cap) \
        else math.inf

    def ell(d):
        d = np.asarray(d, dtype=np.float64)
        v = np.log1p(s * np.minimum(d, d_cap) / h0) / s
        if math.isfinite(d_cap):
            v = v + np.maximum(d - d_cap, 0.0) / h_max
        return v

    def inverse(v):
        v = np.asarray(v, dtype=np.float64)
        d = h0 * np.expm1(s * np.minimum(v, ell_cap)) / s
        if math.isfinite(d_cap):
            d = d + np.maximum(v - ell_cap, 0.0) * h_max
        return np.minimum(d, half)
    return ell, inverse


def quantile_edges(cells: int, length: float, ell, inverse) -> np.ndarray:
    """N + 1 edges at the k/N quantiles of a load symmetric about the
    gap's middle (ell and its inverse over a half-gap, from either wall):
    OSIRIS's equal shares, each half placed from its own wall, so the two
    walls get the same cells to the rounding of L - x (1e-18 m on a 1-cm
    gap)."""
    n = int(cells)
    half = 0.5 * length
    total = 2.0 * float(ell(half))
    k = np.arange(n + 1, dtype=np.float64)
    left = inverse(total * k / n)
    right = length - inverse(total * (n - k) / n)
    edges = np.where(k <= n / 2, left, right)
    edges[0], edges[-1] = 0.0, length
    if n % 2 == 0:
        edges[n // 2] = half
    return edges


def equidistribute(load, cells: int, length: float,
                   samples: int = 400001) -> np.ndarray:
    """The rule for any load density w(x) > 0 on [0, length] (OSIRIS's
    spatial_loaddensity): the N + 1 quantiles of its integral, from a
    trapezoid sum on points that crowd the walls."""
    u = np.linspace(0.0, 1.0, samples)
    x = 0.5 * length * (1.0 - np.cos(np.pi * u))
    w = np.asarray(load(x), dtype=np.float64)
    if np.any(~np.isfinite(w)) or np.any(w <= 0.0):
        raise ValueError("the load density must be finite and positive")
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1])
                                           * np.diff(x))])
    edges = np.interp(cum[-1] * np.arange(cells + 1) / cells, cum, x)
    edges[0], edges[-1] = 0.0, length
    return edges


def _h0(wall_cell: float, s: float) -> float:
    """The load's h0 whose first cell is wall_cell wide."""
    return wall_cell * s / math.expm1(s)


def solve_rate(cells: int, length: float, wall_cell: float, h_max=None):
    """The growth rate s at which N cells, the first wall_cell wide, fill
    the gap: 2 ell(L/2) = N with h0 = wall_cell s/(e^s - 1)."""
    n, half = int(cells), 0.5 * length
    if n * wall_cell >= length:
        raise ValueError(f"{n} cells of {wall_cell * 1e6:.3g} um already fill "
                         f"the {length * 100:.3g}-cm gap: use a uniform mesh")
    if h_max is not None and n * h_max <= length:
        raise ValueError(f"{n} cells of at most {h_max * 1e6:.3g} um cannot "
                         f"fill the {length * 100:.3g}-cm gap")

    def excess(s):
        ell, _ = wall_load(_h0(wall_cell, s), s, h_max, half)
        return 2.0 * float(ell(half)) - n

    lo, hi = 1e-12, 1.0
    while excess(hi) > 0.0:
        hi *= 2.0
        if hi > 50.0:
            raise ValueError("no grading gives this mesh")
    for _ in range(300):                 # the load falls as s grows
        mid = 0.5 * (lo + hi)
        if excess(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 4e-16 * hi:
            break
    return 0.5 * (lo + hi)


def _from_load(cells, length_m, h0, s, h_max, area_m2):
    ell, inverse = wall_load(h0, s, h_max, 0.5 * length_m)
    edges = quantile_edges(cells, length_m, ell, inverse)
    return Mesh(edges, area_m2, kind="graded", load=(h0, s, h_max))


def graded(cells: int, length_m: float, wall_cell_m: float,
           max_cell_m: float = None, area_m2: float = 1.0) -> Mesh:
    """N cells, wall_cell_m at each electrode, growing by a ratio e^s solved
    to fill the gap (and capped at max_cell_m if given)."""
    s = solve_rate(cells, length_m, wall_cell_m, max_cell_m)
    if math.exp(s) > RATIO_MAX:
        raise ValueError(
            f"{cells} cells from {wall_cell_m * 1e6:.3g} um at the walls grow "
            f"by {math.exp(s):.3f} a cell (at most {RATIO_MAX}): use more "
            "cells or a larger wall cell")
    return _from_load(cells, length_m, _h0(wall_cell_m, s), s, max_cell_m,
                      area_m2)


def refine(mesh: Mesh, levels: int = 1) -> Mesh:
    """The mesh with 2^levels times the cells on the same load: every edge
    kept, each cell split at its load midpoints (nested meshes)."""
    if mesh.uniform:
        return uniform(mesh.n_cells * 2 ** levels, mesh.length, mesh.area_m2)
    h0, s, h_max = mesh.load
    f = 2.0 ** levels
    return _from_load(mesh.n_cells * 2 ** levels, mesh.length, h0 / f, s / f,
                      None if h_max is None else h_max / f, mesh.area_m2)
