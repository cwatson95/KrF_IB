"""ase.py gates: the uniform-field oracle, quadrature, traversal, margins."""

import math

import numpy as np

from REBUILD.MODULES.ase import (_segments, parasitic_margins, planck_source,
                                 point_fluxes, ray_intensity,
                                 sphere_quadrature, sum_rule_bandwidth_hz,
                                 window_irradiance)

EXTENTS = np.array([0.02, 0.02, 0.75])
SHAPE = (4, 4, 16)
N = 4 * 4 * 16
LAM, SIG, TAU = 248.0, 5.354e-21, 25e-9


def _uniform(g0, al):
    return np.full(N, g0), np.full(N, al)


def _closed_form(g0, al, L):
    dnu = sum_rule_bandwidth_hz(LAM, SIG, TAU)
    sbd = planck_source(LAM) * dnu
    gn = g0 - al
    if abs(gn * L) < 1e-12:
        return g0 * L * sbd
    return (g0 / gn) * sbd * (math.exp(gn * L) - 1.0)


def test_uniform_field_reproduces_unified_closed_form():
    """The oracle: a uniform voxel field must telescope to unified/ase3d's
    exact uniform-medium solution."""
    g0v, av = 2.0, 0.35
    g0, al = _uniform(g0v, av)
    dnu = sum_rule_bandwidth_hz(LAM, SIG, TAU)
    sbd = planck_source(LAM) * dnu
    I = ray_intensity(np.array([0.01, 0.01, 0.0]), np.array([0.0, 0.0, 1.0]),
                      g0, al, EXTENTS, SHAPE, sbd)
    assert abs(I / _closed_form(g0v, av, 0.75) - 1.0) < 1e-12


def test_traversal_total_length_is_exact():
    p0 = np.array([0.005, 0.013, 0.2])
    d = np.array([0.3, -0.2, 0.93])
    d = d / np.linalg.norm(d)
    segs = _segments(p0, d, EXTENTS, SHAPE)
    # exact chord length to the box wall
    with np.errstate(divide="ignore"):
        t_hi = np.where(d > 0, (EXTENTS - p0) / d, np.inf)
        t_lo = np.where(d < 0, -p0 / d, np.inf)
    L_exact = float(np.min(np.minimum(t_hi, t_lo)))
    assert abs(sum(ds for _, ds in segs) - L_exact) < 1e-12


def test_quadrature_weights_sum_to_4pi():
    _, w = sphere_quadrature(24, 48)
    assert abs(w.sum() - 4.0 * math.pi) < 1e-10


def test_sum_rule_bandwidth_no_measured_linewidth():
    dnu = sum_rule_bandwidth_hz(LAM, SIG, TAU)
    dlam_nm = dnu * (LAM * 1e-9) ** 2 / 2.99792458e8 * 1e9
    assert 2.0 < dlam_nm < 6.0        # the KrF sum-rule class (~3.8 nm)


def test_two_segment_hand_computation_and_order_physics():
    """The voxel treatment against exact hand arithmetic on a two-segment
    chord, in the ARRIVAL order (0.1.25 fix — the pre-fix march applied
    the recurrence outward from the observer, amplifying each source by
    the cells BEHIND it; correct transfer amplifies by the cells BETWEEN
    source and observer): radiation observed at p0 from the d side starts
    at the FAR wall and traverses toward p0, so the far segment's
    contribution is amplified through the near segment. For PURE gain the
    order provably cancels (both orderings give S(ab-1), the fully-
    inverted-amplifier identity — the reason the outward march passed the
    uniform oracle for a year); nonuniform absorption breaks it through
    the g0/g_net source weighting. All hand-computed."""
    dnu = sum_rule_bandwidth_hz(LAM, SIG, TAU)
    S = planck_source(LAM) * dnu
    shape2 = (1, 1, 2)
    ext2 = np.array([0.02, 0.02, 0.30])
    L = 0.15
    p0 = np.array([0.01, 0.01, 0.0])
    dz = np.array([0.0, 0.0, 1.0])
    # pure gain, order swap: identical (S(ab-1))
    a, b = math.exp(2.0 * L), math.exp(0.7 * L)
    I_fwd = ray_intensity(p0, dz, np.array([2.0, 0.7]), np.zeros(2),
                          ext2, shape2, S)
    I_rev = ray_intensity(p0, dz, np.array([0.7, 2.0]), np.zeros(2),
                          ext2, shape2, S)
    assert abs(I_fwd / (S * (a * b - 1.0)) - 1.0) < 1e-12
    assert abs(I_rev / I_fwd - 1.0) < 1e-12
    # gain cell at the observer (index 0), absorbing cell at the far wall
    # (index 1, g_net = 0.5): the FAR cell's source crosses the near gain
    # cell — amplified by e^{1 L} — and the near cell's own source arrives
    # unattenuated. Hand-composed in arrival order:
    g0 = np.array([1.0, 1.0])
    al = np.array([0.0, 0.5])
    gn2 = 0.5
    I_far = S * (1.0 / gn2) * (math.exp(gn2 * L) - 1.0)
    I_hand = I_far * math.exp(1.0 * L) + S * (math.exp(1.0 * L) - 1.0)
    I_num = ray_intensity(p0, dz, g0, al, ext2, shape2, S)
    assert abs(I_num / I_hand - 1.0) < 1e-12
    # and reversing the absorber's position genuinely changes the answer
    I_swap = ray_intensity(p0, dz, g0, al[::-1].copy(), ext2, shape2, S)
    assert abs(I_swap / I_num - 1.0) > 0.02
    # the old outward-order value is provably NOT the arrival value here
    I_outward = (S * (math.exp(1.0 * L) - 1.0)) * math.exp(gn2 * L) \
        + S * (1.0 / gn2) * (math.exp(gn2 * L) - 1.0)
    assert abs(I_outward / I_hand - 1.0) > 0.005


def test_window_quantitative_uniform_gate():
    """The +z-face irradiance against an INDEPENDENT closed-form hemisphere
    integral on a uniform field (same quadrature, chord length from the
    face geometry, uniform-medium intensity per direction). Quantitative
    on purpose: the pre-0.1.25 gate asserted only positivity, and the
    march-direction bug it never caught reported window ASE ~1e7 LOW
    (chords exited the +z wall after ~1 nm — CHANGELOG [0.1.25])."""
    g0v, av = 1.5, 0.2
    g0, al = _uniform(g0v, av)
    dnu = sum_rule_bandwidth_hz(LAM, SIG, TAU)
    sbd = planck_source(LAM) * dnu
    dirs, wts = sphere_quadrature(8, 16)
    p0 = np.array([0.5 * EXTENTS[0], 0.5 * EXTENTS[1],
                   EXTENTS[2] * (1.0 - 1e-9)])
    gn = g0v - av
    ref = 0.0
    for k in range(dirs.shape[0]):
        d = dirs[k]
        if d[2] >= 0.0:
            continue
        with np.errstate(divide="ignore"):
            t_hi = np.where(d > 0, (EXTENTS - p0) / d, np.inf)
            t_lo = np.where(d < 0, -p0 / d, np.inf)
        L = float(np.min(np.minimum(t_hi, t_lo)))
        ref += wts[k] * abs(d[2]) * (g0v / gn) * sbd * (
            math.exp(gn * L) - 1.0)
    win = window_irradiance(g0, al, EXTENTS, SHAPE, LAM, SIG, TAU,
                            n_mu=8, n_phi=16)
    assert abs(win / ref - 1.0) < 1e-9
    fl = point_fluxes(0.5 * EXTENTS, g0, al, EXTENTS, SHAPE, LAM, SIG, TAU,
                      n_mu=8, n_phi=16)
    assert 0.0 < win < fl["scalar_W_m2"]


def test_parasitic_margins_uniform_analytic():
    g0v, av = 2.0, 0.35
    g0, al = _uniform(g0v, av)
    m = parasitic_margins(g0, al, EXTENTS, SHAPE, 0.1, 0.99, 0.0, 0.04)
    ax = m["axial (laser cavity)"]
    ref = 0.1 * 0.99 * math.exp(2.0 * (g0v - av) * 0.75)
    assert abs(ax["round_trip_gain"] / ref - 1.0) < 1e-9
    assert m["transverse x (walls)"]["round_trip_gain"] == 0.0
