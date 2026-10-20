"""3-D ASE / parasitics over the VOXEL gain field (S6, C5).

PHYSICS ADAPTED from ``unified/ase3d.py`` (PLAN (27) — the Lehmberg &
Giuliani-class ray treatment Orestes uses), UPGRADED from the uniform box to
the device grid's non-uniform gain field:

* Radiative transfer along a chord through a fully-inverted medium (the
  excimer lower state is dissociative, so complete inversion is exact):
      dI/ds = g0(s) * S_bb * dnu_eff + (g0(s) - alpha(s)) * I
  integrated cell-by-cell with the EXACT uniform-segment update
      I <- I e^{g_net ds} + (g0/g_net) S_bb dnu_eff (e^{g_net ds} - 1)
  (series-safe as g_net -> 0). For a uniform field this telescopes to
  unified/ase3d's closed form — the oracle gate.

* The band integral uses the sigma<->A sum rule (no measured linewidth
  enters): dnu_eff = lambda^2 A_band / (8 pi sigma_line), with the SAME
  lumped (sigma_eff, A = 1/tau_eff) pair the engine's URAD row radiates —
  so the ASE report and the ledger's radiated_escape describe one medium.

* Ray traversal is Amanatides–Woo DDA over the grid; quadrature is
  Gauss–Legendre in mu x uniform phi (sum of weights = 4 pi).

Scope (v1, honest, same convention as unified/ase3d): quasi-steady ASE
evaluated ON a state of the run (peak gain by default) — a report and a
saturation fixed point, not an in-loop RHS sink. The in-loop coupling is the
recorded follow-up.
"""

from __future__ import annotations

import math

import numpy as np

from .units import C_LIGHT, H_PLANCK


def sum_rule_bandwidth_hz(wavelength_nm: float, sigma_se_m2: float,
                          tau_upper_s: float) -> float:
    lam = wavelength_nm * 1e-9
    A_band = 1.0 / tau_upper_s
    return lam * lam * A_band / (8.0 * math.pi * sigma_se_m2)


def planck_source(wavelength_nm: float) -> float:
    """S_bb = 2 h nu^3 / c^2 [W / m^2 / sr / Hz] — the one-noise-photon-
    per-mode source function of an inverted two-level medium."""
    nu = C_LIGHT / (wavelength_nm * 1e-9)
    return 2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2


def sphere_quadrature(n_mu: int = 24, n_phi: int = 48):
    """(directions (N,3), weights (N,)), sum(weights) = 4 pi."""
    mu, w_mu = np.polynomial.legendre.leggauss(int(n_mu))
    phi = (np.arange(int(n_phi)) + 0.5) * (2.0 * math.pi / int(n_phi))
    s = np.sqrt(np.clip(1.0 - mu ** 2, 0.0, None))
    dirs = np.stack([
        np.outer(s, np.cos(phi)).ravel(),
        np.outer(s, np.sin(phi)).ravel(),
        np.repeat(mu, int(n_phi)),
    ], axis=1)
    wts = np.repeat(w_mu, int(n_phi)) * (2.0 * math.pi / int(n_phi))
    return dirs, wts


def _segments(p0, d, extents, shape):
    """Amanatides–Woo voxel traversal: yields (cell_flat_index, ds) from p0
    along unit direction d to the box wall."""
    p0 = np.asarray(p0, float)
    d = np.asarray(d, float)
    extents = np.asarray(extents, float)
    shape = np.asarray(shape, int)
    h = extents / shape
    idx = np.clip((p0 / h).astype(int), 0, shape - 1)
    step = np.where(d > 0, 1, -1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_delta = np.abs(h / d)
        nxt = np.where(d > 0, (idx + 1) * h, idx * h)
        t_max = np.where(d != 0, (nxt - p0) / d, np.inf)
    t_max = np.where(np.isnan(t_max), np.inf, t_max)
    t = 0.0
    out = []
    ny, nz = int(shape[1]), int(shape[2])
    while True:
        ax = int(np.argmin(t_max))
        t_next = float(t_max[ax])
        ds = t_next - t
        if ds > 0.0:
            flat = (idx[0] * ny + idx[1]) * nz + idx[2]
            out.append((int(flat), ds))
        idx[ax] += step[ax]
        if idx[ax] < 0 or idx[ax] >= shape[ax]:
            break
        t = t_next
        t_max[ax] += t_delta[ax]
    return out


def ray_intensity(p0, d, g0_flat, alpha_flat, extents, shape, s_bb_dnu):
    """ASE intensity [W/m^2/sr] arriving AT p0 from the d side: radiation
    that propagates along -d, born and amplified on the chord
    {p0 + s d, s in [0, L]} and observed at p0.

    EXACT for any voxel field (0.1.25): the chord's segments are collected
    outward from p0 but the transfer recurrence is applied in REVERSED
    order — far-wall segment first, the segment nearest p0 last — so each
    segment's source is amplified by exactly the cells BETWEEN it and the
    observer. The pre-0.1.25 outward-order march amplified each source by
    the cells BEHIND it instead, which is the same number only for a
    uniform medium (the two orderings coincide there because g0/g_net is
    constant); on gain-core/absorbing-skin fields it biased the 4-pi
    scalar flux ~3% and single rays up to 12% (found 2026-08-25). The
    uniform-field closed form — the oracle gate — is unchanged."""
    segs = _segments(p0, d, extents, shape)
    I = 0.0
    for flat, ds in reversed(segs):
        g0 = g0_flat[flat]
        gn = g0 - alpha_flat[flat]
        x = gn * ds
        if abs(x) < 1e-8:
            e1 = ds * (1.0 + 0.5 * x)
            I = I * math.exp(x) + g0 * s_bb_dnu * e1
        else:
            ex = math.exp(min(x, 700.0))
            I = I * ex + (g0 / gn) * s_bb_dnu * (ex - 1.0)
    return I


def point_fluxes(point, g0_flat, alpha_flat, extents, shape, wavelength_nm,
                 sigma_se_m2, tau_upper_s, n_mu: int = 24, n_phi: int = 48):
    """Angular-integrated ASE at an interior point: the 4-pi scalar flux
    [W/m^2] (what saturates a homogeneously broadened line) and the mean
    axial intensity [W/m^2/sr]."""
    dnu = sum_rule_bandwidth_hz(wavelength_nm, sigma_se_m2, tau_upper_s)
    s_bb_dnu = planck_source(wavelength_nm) * dnu
    dirs, wts = sphere_quadrature(n_mu, n_phi)
    total = 0.0
    for k in range(dirs.shape[0]):
        total += wts[k] * ray_intensity(point, dirs[k], g0_flat, alpha_flat,
                                        extents, shape, s_bb_dnu)
    Iz = 0.5 * (ray_intensity(point, np.array([0.0, 0.0, 1.0]), g0_flat,
                              alpha_flat, extents, shape, s_bb_dnu)
                + ray_intensity(point, np.array([0.0, 0.0, -1.0]), g0_flat,
                                alpha_flat, extents, shape, s_bb_dnu))
    return {"scalar_W_m2": total, "along_z_W_m2_sr": Iz}


def window_irradiance(g0_flat, alpha_flat, extents, shape, wavelength_nm,
                      sigma_se_m2, tau_upper_s, n_mu: int = 24,
                      n_phi: int = 48) -> float:
    """ASE irradiance at the centre of the +z face [W/m^2]: hemisphere
    integral of I |cos theta| dOmega over rays arriving from inside.

    Convention (matches ray_intensity): a kept quadrature direction d has
    d_z < 0 and is the ARRIVAL side — the chord from the face point along
    d runs INTO the box; the radiation observed propagates along -d
    (d_z<0 -> u_z>0, hitting the face from inside). BUG FIXED 2026-08-25:
    the pre-0.1.25 code marched along -d — from a point 1 nm inside the
    +z wall straight OUT of the box — so every chord was ~1 nm long and
    the reported irradiance sat ~1e7 LOW (measured 3.37 W/m^2 where the
    uniform-field value is 4.6e7). The suite's gate only asserted
    positivity; test_ase now pins the uniform-field hemisphere value
    quantitatively."""
    dnu = sum_rule_bandwidth_hz(wavelength_nm, sigma_se_m2, tau_upper_s)
    s_bb_dnu = planck_source(wavelength_nm) * dnu
    dirs, wts = sphere_quadrature(n_mu, n_phi)
    point = np.array([0.5 * extents[0], 0.5 * extents[1],
                      extents[2] * (1.0 - 1e-9)])
    total = 0.0
    for k in range(dirs.shape[0]):
        if dirs[k, 2] >= 0.0:
            continue                       # not an arrival side of +z
        I = ray_intensity(point, dirs[k], g0_flat, alpha_flat, extents,
                          shape, s_bb_dnu)
        total += wts[k] * abs(dirs[k, 2]) * I
    return total


def path_gain(p0, d, g0_flat, alpha_flat, extents, shape) -> float:
    """integral of (g0 - alpha) ds along the chord from p0 in direction d."""
    return sum((g0_flat[f] - alpha_flat[f]) * ds
               for f, ds in _segments(p0, d, extents, shape))


def parasitic_margins(g0_flat, alpha_flat, extents, shape, R_mirror1: float,
                      R_mirror2: float, R_wall: float = 0.0,
                      R_window: float = 0.04) -> dict:
    """Round-trip gain of the closed parasitic path families through the
    box centre; oscillation at >= 1 (unified/ase3d's table, evaluated on
    the ACTUAL field via chord integrals)."""
    c = 0.5 * np.asarray(extents, float)
    def chord(dvec):
        d = np.asarray(dvec, float)
        d = d / np.linalg.norm(d)
        return (path_gain(c, d, g0_flat, alpha_flat, extents, shape)
                + path_gain(c, -d, g0_flat, alpha_flat, extents, shape))
    gl_axial = chord([0, 0, 1.0])
    gl_x = chord([1.0, 0, 0])
    gl_y = chord([0, 1.0, 0])
    gl_diag = chord(list(extents / np.linalg.norm(extents)))
    paths = {
        "axial (laser cavity)": (R_mirror1 * R_mirror2, gl_axial),
        "transverse x (walls)": (R_wall ** 2, gl_x),
        "transverse y (walls)": (R_wall ** 2, gl_y),
        "body diagonal (walls)": (R_wall ** 2, gl_diag),
        "windows axial (Fresnel)": (R_window ** 2, gl_axial),
    }
    out = {}
    for name, (Rprod, gl) in paths.items():
        G = Rprod * math.exp(min(2.0 * gl, 700.0))
        out[name] = {"round_trip_gain": G, "gl_one_way": gl,
                     "R_product": Rprod, "oscillates": bool(G >= 1.0)}
    return out


def ase_report(result: dict, at: str = "peak_gain", n_mu: int = 24,
               n_phi: int = 48) -> dict:
    """Quasi-steady ASE/parasitics report for a finished engine run.

    ``at``: 'peak_gain' evaluates on the saved state nearest the g0 peak
    (v1 keeps the final full state only, so peak-gain uses the volume-
    average trace to pick the time and scales the final FIELD SHAPE to it —
    stated; 'final' uses the final state exactly)."""
    dev = result["device"]
    lay = dev.lay
    y = result["y_final"]
    dens = np.maximum(lay.density_block(y), 0.0)
    opt = dev.optics
    i_up = lay.index[opt["upper"]]
    g0 = dev.cav.sigma_se_m2 * dens[i_up].astype(float)
    alpha = np.zeros_like(g0)
    for sp, sig in dev.cav.absorbers.items():
        i_sp = lay.index.get(sp)
        if i_sp is not None:
            alpha += sig * dens[i_sp]
    scale_note = "final state"
    if at == "peak_gain":
        g0_t = np.asarray(result["g0_m1"], float)
        j_pk = int(np.argmax(g0_t))
        g_fin = float(g0_t[-1])
        if g_fin > 0.0 and g0_t[j_pk] > 0.0:
            s = g0_t[j_pk] / g_fin
            g0 = g0 * s
            alpha = alpha * s          # same shape assumption, stated
            scale_note = (f"final field shape scaled x{s:.3g} to the "
                          f"t = {result['t_ns'][j_pk]:.0f} ns gain peak")
    extents = dev.extents
    shape = dev.shape
    lam = dev.cav.wavelength_nm
    sig_se = dev.cav.sigma_se_m2
    tau = opt["tau_upper_s"]
    centre = 0.5 * extents
    fl = point_fluxes(centre, g0, alpha, extents, shape, lam, sig_se, tau,
                      n_mu, n_phi)
    win = window_irradiance(g0, alpha, extents, shape, lam, sig_se, tau,
                            n_mu, n_phi)
    margins = parasitic_margins(g0, alpha, extents, shape,
                                dev.cav.R_OC, 0.99)
    dnu = sum_rule_bandwidth_hz(lam, sig_se, tau)
    return {
        "ok": True, "state": scale_note,
        "dnu_eff_hz": dnu,
        "dlam_eff_nm": dnu * (lam * 1e-9) ** 2 / C_LIGHT * 1e9,
        "g0_max_m1": float(g0.max()), "g0_mean_m1": float(g0.mean()),
        "ase_scalar_W_m2": fl["scalar_W_m2"],
        "ase_window_W_m2": win,
        "parasitic_margins": margins,
        "worst_parasitic": max(margins.items(),
                               key=lambda kv: kv[1]["round_trip_gain"])[0],
    }
