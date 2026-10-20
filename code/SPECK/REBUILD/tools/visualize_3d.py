#!/usr/bin/env python3
"""3-D device-grid visualization — the per-cell plasma state over time.

A VIEWER, NOT AN ORACLE: this script renders what ``Device.run`` recorded
(``res["y_t"]``, kept for every multi-cell run since 0.1.15) and feeds
nothing back into the engine. The per-cell fields are reconstructed with
the engine's own accessors (``StateLayout.density_block`` /
``electron_density``, ``Device.te_of``) and the gain expression is
``laser.gain_m1``'s, vectorized over cells — the numbers here ARE the
run's numbers, only reshaped to (nx, ny, nz).

What it renders, per save time:
  n_e           quasineutral electron density [m^-3]        (log scale)
  <upper>       the lasing level's density, KrF*/ArF* [m^-3] (log scale)
  g0            local small-signal gain [1/m]     (diverging, 0-centred)
  Te            electron temperature [eV]
  Tg            gas temperature [K]
  deposition    the static beam-deposition shape [peak-normalized]
plus the cavity time strip (P_out, <g0>) with a cursor at the frame time.

Outputs (all under --out, default ``<root>/rebuild_out/viz_<preset>/``):
  viz_data.npz     the exact fields + meta (the data of record here;
                   re-render from it with --npz, skipping the run)
  contact_sheet.png     n_e snapshots on one page
  frames/frame_*.png -> movie.gif (Pillow) and optionally movie.mp4
                   (ffmpeg, --mp4)
  viewer.html      self-contained interactive viewer (canvas voxels,
                   time/yaw sliders, field selector; embeds the fields
                   QUANTIZED to uint8 — read exact values from the npz)

Honesty, restated on every figure: the grid is the engine's coarse device
grid — cells are independent kinetic reactors coupled by the shared cavity
photon field and the volume-averaged ASE overlay; there is no inter-cell
transport in the default presets, and this is not a PIC/fluid solve.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/visualize_3d.py --preset krf_beam_3d
    python3 SPECK/REBUILD/tools/visualize_3d.py --npz SPECK/rebuild_out/viz_krf_beam_3d/viz_data.npz
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")


# --------------------------------------------------------------------- run
def build_device(preset: str, t_end_ns=None, tier=None):
    """Build the Device exactly the way run.py does (same presets dict)."""
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    from REBUILD.run import PRESETS
    spec = PRESETS[preset]
    P = engine_parameters()
    P.set_many(spec["params"], source="problem")
    if tier:
        P.set("electron_tier", tier, source="cli")
    if t_end_ns is not None:
        P.set("t_end_ns", float(t_end_ns), source="cli")
    return Device(P, mixture=spec.get("mixture"))


# ------------------------------------------------------------------ fields
def extract_fields(device, res) -> dict:
    """res (with y_t) -> {name: (n_t, nx, ny, nz)} + geometry/meta arrays.

    Pure reads; every number comes from the run state via the engine's own
    accessors. Raises if the run was not collected densely.
    """
    if "y_t" not in res:
        raise ValueError(
            "res has no 'y_t' — run with dense_cells=True (multi-cell "
            "runs keep it by default since 0.1.15)")
    lay = device.lay
    y_t = res["y_t"]
    n_t = y_t.shape[1]
    nx, ny, nz = device.shape
    upper = device.optics["upper"]
    cav = device.cav
    # gain_m1's expression, per cell: n_upper*sigma_se - sum absorbers
    abs_rows = [(lay.index[sp], sig) for sp, sig in cav.absorbers.items()
                if sp in lay.index]

    def grid(v):                       # flat cell vector -> (nx, ny, nz)
        return np.asarray(v).reshape(nx, ny, nz)

    ne = np.empty((n_t, nx, ny, nz))
    nup = np.empty_like(ne)
    g0 = np.empty_like(ne)
    te = np.empty_like(ne)
    tg = np.empty_like(ne)
    for j in range(n_t):
        yj = y_t[:, j]
        dens = np.maximum(lay.density_block(yj), 0.0)
        ne[j] = grid(lay.electron_density(dens))
        nup[j] = grid(dens[lay.index[upper]])
        g = dens[lay.index[upper]] * cav.sigma_se_m2
        for i, sig in abs_rows:
            g = g - dens[i] * sig
        g0[j] = grid(g)
        te[j] = grid(device.te_of(yj))
        tg[j] = grid(lay.cells_view(yj)[lay.i_tg])
    out = {
        "t_ns": np.asarray(res["t_ns"], float),
        "shape": np.array([nx, ny, nz]),
        "extents_cm": np.asarray(device.extents, float) * 100.0,
        "cell_vols": np.asarray(device.cell_vols, float).copy(),
        "P_out_W": np.asarray(res["P_out_W"], float),
        "g0avg_m1": np.asarray(res["g0_m1"], float),
        "fields": {
            "n_e": ne, upper: nup, "g0": g0, "Te": te, "Tg": tg,
        },
        "meta": {
            "upper": upper,
            "regime": res.get("regime", ""),
            "mode": device.mode,
            "stamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        },
    }
    if device._dep_frac is not None:
        dep = device._dep_frac / device.cell_vols          # ~W/m^3 shape
        out["dep_shape"] = grid(dep / dep.max())
    return out


FIELD_UNITS = {"n_e": "m^-3", "g0": "1/m", "Te": "eV", "Tg": "K"}


def field_specs(F) -> dict:
    """Per-field global (whole-run) scale decisions, shared by every
    renderer so the PNGs and the HTML viewer agree by construction."""
    specs = {}
    for name, vals in F["fields"].items():
        vmax = float(vals.max())
        vmin = float(vals.min())
        unit = FIELD_UNITS.get(name, "m^-3")
        if name in ("n_e",) or name == F["meta"]["upper"]:
            vmax = vmax if vmax > 0 else 1.0
            specs[name] = {"scale": "log", "vmin": vmax * 1e-5,
                           "vmax": vmax, "cmap": "viridis", "unit": unit}
        elif name == "g0":
            a = max(abs(vmin), abs(vmax), 1e-30)
            if vmin < 0.0 < vmax:
                specs[name] = {"scale": "lin", "vmin": -a, "vmax": a,
                               "cmap": "diverging", "unit": unit}
            else:
                specs[name] = {"scale": "lin", "vmin": vmin, "vmax": vmax,
                               "cmap": "diverging", "unit": unit}
        else:
            pad = 0.05 * max(vmax - vmin, 1e-30)
            specs[name] = {"scale": "lin", "vmin": vmin - pad,
                           "vmax": vmax + pad, "cmap": "viridis",
                           "unit": unit}
    if "dep_shape" in F:
        specs["deposition"] = {"scale": "lin", "vmin": 0.0, "vmax": 1.0,
                               "cmap": "magma", "unit": "peak-norm."}
    return specs


# ------------------------------------------------------------------- npz
def save_npz(path, F) -> None:
    arrays = {"t_ns": F["t_ns"], "shape": F["shape"],
              "extents_cm": F["extents_cm"], "P_out_W": F["P_out_W"],
              "g0avg_m1": F["g0avg_m1"],
              "meta_json": np.array(json.dumps(F["meta"]))}
    if "cell_vols" in F:
        arrays["cell_vols"] = F["cell_vols"]
    for name, vals in F["fields"].items():
        arrays[f"field_{name}"] = vals
    if "dep_shape" in F:
        arrays["dep_shape"] = F["dep_shape"]
    np.savez_compressed(Path(path), **arrays)


def load_npz(path) -> dict:
    z = np.load(Path(path), allow_pickle=False)
    F = {"t_ns": z["t_ns"], "shape": z["shape"],
         "extents_cm": z["extents_cm"], "P_out_W": z["P_out_W"],
         "g0avg_m1": z["g0avg_m1"],
         "meta": json.loads(str(z["meta_json"])), "fields": {}}
    for key in z.files:
        if key.startswith("field_"):
            F["fields"][key[len("field_"):]] = z[key]
    if "dep_shape" in z.files:
        F["dep_shape"] = z["dep_shape"]
    n_cells = int(np.prod(F["shape"]))
    if "cell_vols" in z.files:
        F["cell_vols"] = z["cell_vols"]
    else:                       # pre-0.1.16 npz: presets used uniform grids
        v = float(np.prod(F["extents_cm"] / 100.0)) / n_cells
        F["cell_vols"] = np.full(n_cells, v)
    return F


# ------------------------------------------------------- particle sampling
# The particle view draws each field as a cloud of dots: dot COUNT per cell
# is proportional to the number of real particles there (density x cell
# volume), so one dot carries a fixed, stated macro-weight. Dot POSITIONS
# come from a deterministic per-cell hash lattice — a SAMPLING of the
# density field, not trajectories (the engine is a rate-equation model;
# real macro-particles live in PIC_MCC_3D). The same hash runs in the HTML
# viewer's JS (mulberry32), so every renderer shows the same cloud.

def _mulberry32(seeds):
    """Vectorized mulberry32 -> uniforms in [0, 1). uint32 wraparound
    matches the JS Math.imul implementation bit for bit."""
    M = np.uint64(0xFFFFFFFF)
    a = (np.asarray(seeds, np.uint64) + np.uint64(0x6D2B79F5)) & M
    t = a
    t = (t ^ (t >> np.uint64(15))) * (t | np.uint64(1)) & M
    t = (t ^ (t + ((t ^ (t >> np.uint64(7))) * (t | np.uint64(61)) & M))) & M
    return ((t ^ (t >> np.uint64(14))) & M).astype(float) / 2.0 ** 32


def sample_positions(n_cells, k_max, species_idx):
    """(n_cells, k_max, 3) unit-cell offsets, deterministic (hash lattice)."""
    c = np.arange(n_cells, dtype=np.uint64)[:, None, None]
    k = np.arange(k_max, dtype=np.uint64)[None, :, None]
    ax = np.arange(3, dtype=np.uint64)[None, None, :]
    seed = (c * np.uint64(9781) + k * np.uint64(2621)
            + np.uint64(species_idx) * np.uint64(131071)
            + ax * np.uint64(40503))
    return _mulberry32(seed)


def particle_counts(vals_t, cell_vols, n_target):
    """Dots per (frame, cell) linear in real particle number, plus the
    macro-weight (real particles per dot). Peak-frame total ~ n_target."""
    per = np.asarray(vals_t, float).reshape(len(vals_t), -1) * cell_vols
    norm = float(per.sum(axis=1).max())
    if norm <= 0.0:
        return np.zeros(per.shape, int), 0.0
    counts = np.rint(per / norm * n_target).astype(int)
    return counts, norm / n_target


# ------------------------------------------------------------- quantization
def quantize(vals, vmin, vmax, scale):
    """Field values -> uint8 for the HTML viewer (dequantize inverts to
    within one bin; the npz keeps the exact numbers)."""
    v = np.asarray(vals, float)
    if scale == "log":
        v = np.log10(np.clip(v, vmin, vmax))
        lo, hi = np.log10(vmin), np.log10(vmax)
    else:
        v = np.clip(v, vmin, vmax)
        lo, hi = vmin, vmax
    q = np.rint((v - lo) / max(hi - lo, 1e-300) * 255.0)
    return q.astype(np.uint8)


def dequantize(q, vmin, vmax, scale):
    f = np.asarray(q, float) / 255.0
    if scale == "log":
        lo, hi = np.log10(vmin), np.log10(vmax)
        return 10.0 ** (lo + f * (hi - lo))
    return vmin + f * (vmax - vmin)


# --------------------------------------------------------------- rendering
_CAPTION = ("coarse device grid — cells are independent kinetic reactors "
            "coupled by the shared cavity photon field and volume-averaged "
            "ASE; no inter-cell transport; not a PIC/fluid solve")


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return matplotlib, plt


def _norm_of(mpl, spec):
    from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
    if spec["scale"] == "log":
        return LogNorm(vmin=spec["vmin"], vmax=spec["vmax"])
    if spec["cmap"] == "diverging" and spec["vmin"] < 0.0 < spec["vmax"]:
        return TwoSlopeNorm(vmin=spec["vmin"], vcenter=0.0,
                            vmax=spec["vmax"])
    return Normalize(vmin=spec["vmin"], vmax=spec["vmax"])


def _cmap_of(mpl, spec):
    name = {"viridis": "viridis", "diverging": "RdBu_r",
            "magma": "magma"}[spec["cmap"]]
    return mpl.colormaps[name]


def _normed(norm, vals3d):
    """norm() over a 3-D block — some norms' scale transforms only take
    1-D input, so normalize raveled and reshape."""
    flat = np.ma.filled(norm(vals3d.ravel()), 0.0)
    return np.clip(flat, 0.0, 1.0).reshape(vals3d.shape)


def _voxel_panel(mpl, ax, vals3d, spec, extents_cm, aspect, title):
    norm, cmap = _norm_of(mpl, spec), _cmap_of(mpl, spec)
    nvals = _normed(norm, vals3d)
    rgba = cmap(nvals)
    rgba[..., 3] = 0.30 + 0.65 * np.abs(2.0 * (nvals - 0.5)) \
        if spec["cmap"] == "diverging" else 0.30 + 0.65 * nvals
    filled = np.ones(vals3d.shape, bool)
    ax.voxels(filled, facecolors=rgba, edgecolor=(0, 0, 0, 0.18),
              linewidth=0.2)
    nx, ny, nz = vals3d.shape
    if aspect == "true":
        ax.set_box_aspect(tuple(extents_cm))
    else:
        ax.set_box_aspect((1.0, 1.0, 1.0))
    for setter, n, l_cm in ((ax.set_xticks, nx, extents_cm[0]),
                            (ax.set_yticks, ny, extents_cm[1]),
                            (ax.set_zticks, nz, extents_cm[2])):
        setter([0, n])
    ax.set_xticklabels(["0", f"{extents_cm[0]:g}"], fontsize=6)
    ax.set_yticklabels(["0", f"{extents_cm[1]:g}"], fontsize=6)
    ax.set_zticklabels(["0", f"{extents_cm[2]:g}"], fontsize=6)
    ax.set_xlabel("x cm", fontsize=7, labelpad=-6)
    ax.set_ylabel("y cm", fontsize=7, labelpad=-6)
    ax.set_zlabel("z cm", fontsize=7, labelpad=-6)
    ax.set_title(title, fontsize=9)
    ax.tick_params(pad=-2)
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cb = ax.figure.colorbar(sm, ax=ax, shrink=0.60, pad=0.12,
                            fraction=0.05)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(spec["unit"], fontsize=7)


def _time_strip(fig, gs_cell, F, j):
    axt = fig.add_subplot(gs_cell)
    t, p_out, g0avg = F["t_ns"], F["P_out_W"], F["g0avg_m1"]
    if np.any(p_out > 0):
        axt.semilogy(t, np.clip(p_out, p_out[p_out > 0].min(), None),
                     color="#c96a2b", lw=1.4, label="P_out")
    else:
        axt.plot(t, p_out, color="#c96a2b", lw=1.4, label="P_out")
    axt.set_ylabel("P_out W", fontsize=8)
    axt.set_xlabel("t ns", fontsize=8)
    axg = axt.twinx()
    axg.plot(t, g0avg, color="#3f7ec0", lw=1.2, label="<g0>")
    axg.set_ylabel("<g0> 1/m", fontsize=8)
    for a in (axt, axg):
        a.tick_params(labelsize=7)
    axt.axvline(t[j], color="#888", lw=1.0, ls="--")
    axt.set_xlim(t[0], t[-1])
    return axt


def render_frame(F, specs, j, aspect="unit", fig=None):
    """One full figure at save index j: field panels + cavity time strip."""
    mpl, plt = _mpl()
    if fig is None:
        fig = plt.figure(figsize=(13.0, 8.2), dpi=110)
    fig.clf()
    panels = [("n_e", "n_e"), (F["meta"]["upper"], F["meta"]["upper"]),
              ("g0", "g0 (local gain)"), ("Te", "Te"), ("Tg", "Tg")]
    if "dep_shape" in F:
        panels.append(("deposition", "deposition shape (static)"))
    gs = fig.add_gridspec(3, 3, height_ratios=(1.0, 1.0, 0.34),
                          hspace=0.34, wspace=0.30,
                          left=0.03, right=0.96, top=0.90, bottom=0.07)
    for k, (name, title) in enumerate(panels[:6]):
        ax = fig.add_subplot(gs[divmod(k, 3)], projection="3d")
        vals = (F["dep_shape"] if name == "deposition"
                else F["fields"][name][j])
        _voxel_panel(mpl, ax, vals, specs[name], F["extents_cm"],
                     aspect, title)
    _time_strip(fig, gs[2, :], F, j)
    asp = ("true aspect" if aspect == "true"
           else "axes independently normalized")
    nx, ny, nz = (int(v) for v in F["shape"])
    fig.suptitle(
        f"REBUILD {F['meta'].get('mode', '')} 3-D device grid "
        f"{nx}x{ny}x{nz} — t = {F['t_ns'][j]:.1f} ns — {asp}",
        fontsize=12)
    fig.text(0.5, 0.005, _CAPTION, ha="center", fontsize=7,
             color="#777")
    return fig


# ------------------------------------------------------------- glow view
# The camera view: brightness = what the plasma actually EMITS. Two
# physically-motivated components, each tone-mapped to its own run peak
# (relative display gains stated on every figure):
#   fluorescence ∝ n_upper   — the KrF*/ArF* B→X spontaneous emission,
#                              overwhelmingly the brightest light in the
#                              cell (rendered violet-white; the real line
#                              is 248/193 nm UV);
#   afterglow    ∝ n_e²      — the recombination continuum, orders
#                              dimmer (rendered deep blue).
# Fields are trilinearly interpolated between cell centers FOR DISPLAY
# (stated); the model still resolves only its coarse device grid.

# tone: exponent on the per-frame brightness fraction (1.0 = linear).
# The continuum gets √ so the AFTERGLOW stays visible — n_e² falls five
# decades after beam-off and a linear map renders the cooldown black.
# Stated on every caption; the npz keeps the exact fields.
GLOW_COMPONENTS = [
    {"key": "fluor", "label": "KrF* fluorescence", "color": "#d9c8ff",
     "core": "#ffffff", "gain": 1.0, "tone": 1.0},
    {"key": "recomb", "label": "n_e² recombination continuum",
     "color": "#3b6fd4", "core": "#7fa8ff", "gain": 0.35, "tone": 0.5},
]


def emission_fields(F):
    """(n_t, nx, ny, nz) emission proxies per glow component."""
    upper = F["meta"]["upper"]
    return {"fluor": F["fields"][upper],
            "recomb": F["fields"]["n_e"] ** 2}


def refine_field(vals3d, factor=4):
    """Trilinear display interpolation to a factor-finer grid (order-1
    zoom = exactly linear between cell centers)."""
    from scipy.ndimage import zoom
    return np.maximum(zoom(vals3d, factor, order=1, mode="nearest",
                           grid_mode=True), 0.0)


def glow_cloud(vals3d, n_target, seed_offset, factor=4):
    """Sample a luminous point cloud from the refined field: per-fine-
    voxel counts linear in emission, hash-lattice offsets (deterministic,
    same mulberry32 as the particle view). Returns (n, 3) positions in
    CELL-index coordinates."""
    fine = refine_field(vals3d, factor)
    flat = fine.reshape(-1)
    total = float(flat.sum())
    if total <= 0.0 or n_target <= 0:
        return np.zeros((0, 3))
    # largest-remainder (cumulative-rounding) allocation: per-voxel rint
    # floors a small budget spread over many voxels to zero everywhere
    # (the deep afterglow rendered black); this keeps the total EXACT
    # and stays deterministic
    cum = np.cumsum(flat) / total * n_target
    counts = np.diff(np.rint(np.concatenate(([0.0], cum)))).astype(int)
    k_max = max(int(counts.max()), 1)
    pos = sample_positions(flat.size, k_max, seed_offset)
    nx, ny, nz = fine.shape
    idx = np.indices((nx, ny, nz)).reshape(3, -1).T
    pts = []
    for c in np.nonzero(counts)[0]:
        pts.append((idx[c] + pos[c, : counts[c]]) / factor)
    return np.concatenate(pts) if pts else np.zeros((0, 3))


PARTICLE_SPECIES = [("n_e", "#6ec6ff", "electrons"),
                    (None, "#ffab4d", "upper level")]   # None -> meta upper


def prepare_particles(F, pts=8000):
    """Per-species dot bookkeeping for the particle view: counts per
    (frame, cell), the constant macro-weight, and the fixed hash-lattice
    positions. Same numbers the HTML viewer computes from its payload."""
    n_cells = int(np.prod(F["shape"]))
    out = []
    for s_idx, (name, color, label) in enumerate(PARTICLE_SPECIES):
        fname = F["meta"]["upper"] if name is None else name
        vals = F["fields"][fname].reshape(len(F["t_ns"]), n_cells)
        counts, weight = particle_counts(vals, F["cell_vols"], pts)
        k_max = max(int(counts.max()), 1)
        pos = sample_positions(n_cells, k_max, s_idx)
        out.append({"field": fname, "color": color, "label": label,
                    "counts": counts, "weight": weight, "pos": pos})
    return out


def render_particle_frame(F, pcloud, j, aspect="unit", fig=None):
    """One particle-view figure at save index j: density-sampled dot
    clouds for electrons + the upper level, in the chamber outline."""
    mpl, plt = _mpl()
    if fig is None:
        fig = plt.figure(figsize=(11.0, 8.2), dpi=110)
    fig.clf()
    fig.patch.set_facecolor("#101318")
    gs = fig.add_gridspec(2, 1, height_ratios=(1.0, 0.22), hspace=0.16,
                          left=0.02, right=0.98, top=0.92, bottom=0.07)
    ax = fig.add_subplot(gs[0], projection="3d")
    ax.set_facecolor("#101318")
    nx, ny, nz = (int(v) for v in F["shape"])
    shape = np.array([nx, ny, nz], float)
    idx = np.indices((nx, ny, nz)).reshape(3, -1).T          # (ncells, 3)
    legend = []
    for P in pcloud:
        cnt = P["counts"][j]
        total = int(cnt.sum())
        if total:
            xs, ys, zs = [], [], []
            for c in np.nonzero(cnt)[0]:
                p = idx[c] + P["pos"][c, : cnt[c]]           # index coords
                xs.append(p[:, 0]); ys.append(p[:, 1]); zs.append(p[:, 2])
            ax.scatter(np.concatenate(xs), np.concatenate(ys),
                       np.concatenate(zs), s=2.6, c=P["color"],
                       alpha=0.65, linewidths=0, depthshade=False)
        legend.append(f"{P['label']} ({P['field']}): {total} dots x "
                      f"{P['weight']:.2e} real")
    # chamber outline
    for s, e in [((0, 0, 0), (nx, 0, 0)), ((0, 0, 0), (0, ny, 0)),
                 ((0, 0, 0), (0, 0, nz)), ((nx, ny, nz), (0, ny, nz)),
                 ((nx, ny, nz), (nx, 0, nz)), ((nx, ny, nz), (nx, ny, 0)),
                 ((nx, 0, 0), (nx, ny, 0)), ((nx, 0, 0), (nx, 0, nz)),
                 ((0, ny, 0), (nx, ny, 0)), ((0, ny, 0), (0, ny, nz)),
                 ((0, 0, nz), (nx, 0, nz)), ((0, 0, nz), (0, ny, nz))]:
        ax.plot(*zip(s, e), color="#3a4250", lw=0.8)
    ax.set_xlim(0, nx); ax.set_ylim(0, ny); ax.set_zlim(0, nz)
    ax.set_box_aspect(tuple(F["extents_cm"]) if aspect == "true"
                      else (1.0, 1.0, 1.0))
    ax.set_axis_off()
    ex = F["extents_cm"]
    ax.text2D(0.01, 0.02, f"x {ex[0]:g} cm   y {ex[1]:g} cm   "
              f"z {ex[2]:g} cm", transform=ax.transAxes, fontsize=8,
              color="#8a919c")
    for k, txt in enumerate(legend):
        ax.text2D(0.01, 0.965 - 0.045 * k, txt, transform=ax.transAxes,
                  fontsize=9, color=PARTICLE_SPECIES[k][1])
    # spatial-structure honesty: the measured cell-to-cell spread, so a
    # uniform-looking cloud states that the PLASMA is uniform (x1.01 at
    # the Electra-class preset; x1e5+ in the Bragg-stopped preset)
    ne_j = F["fields"]["n_e"][j].ravel()
    spread = float(ne_j.max() / max(ne_j.min(), 1e-300))
    label = f"×{spread:.2f}" if spread < 1e4 else f"×{spread:.2e}"
    ax.text2D(0.01, 0.965 - 0.045 * len(legend),
              f"cell n_e spread {label}", transform=ax.transAxes,
              fontsize=8, color="#8a919c")
    st = _time_strip(fig, gs[1], F, j)
    st.set_facecolor("#171b22")
    for a in fig.axes[-2:]:                     # strip + its twin
        a.tick_params(colors="#aab0b8", labelsize=7)
        a.xaxis.label.set_color("#aab0b8")
        a.yaxis.label.set_color("#aab0b8")
        for sp in a.spines.values():
            sp.set_color("#3a4250")
    fig.suptitle(
        f"REBUILD {F['meta'].get('mode', '')} — particle view — "
        f"t = {F['t_ns'][j]:.1f} ns — grid {nx}x{ny}x{nz}",
        fontsize=12, color="#d8dbe0")
    fig.text(0.5, 0.005,
             "dot count ∝ particles per cell (macro-weight above); "
             "positions = density-sampling lattice, NOT trajectories — "
             "independent-reactor device grid, not a PIC solve",
             ha="center", fontsize=7, color="#767c85")
    return fig


def prepare_glow(F, pts=30000):
    """Per-component emission fields + run-peak frame sums, so the dot
    budget of each frame is linear in how much light the cell emits."""
    em = emission_fields(F)
    n_t = len(F["t_ns"])
    prep = {"em": em, "pts": int(pts), "sums": {}, "peaks": {}}
    for comp in GLOW_COMPONENTS:
        s = em[comp["key"]].reshape(n_t, -1).sum(axis=1)
        prep["sums"][comp["key"]] = s
        prep["peaks"][comp["key"]] = float(s.max()) if s.max() > 0 else 1.0
    return prep


def _dark_axes(ax_list):
    for a in ax_list:
        a.tick_params(colors="#aab0b8", labelsize=7)
        a.xaxis.label.set_color("#aab0b8")
        a.yaxis.label.set_color("#aab0b8")
        for sp in a.spines.values():
            sp.set_color("#3a4250")


def render_glow_frame(F, prep, j, aspect="true", fig=None):
    """The camera view: a dark vessel with the light the plasma emits.

    Luminous point clouds sampled linearly from the trilinearly refined
    emission fields, layered halo+core scatter for glow (plus a PIL
    bloom pass in render_animation); true aspect by default — the cell
    looks like the long thin duct it is."""
    mpl, plt = _mpl()
    if fig is None:
        fig = plt.figure(figsize=(12.5, 7.6), dpi=110)
    fig.clf()
    fig.patch.set_facecolor("#07090c")
    gs = fig.add_gridspec(2, 1, height_ratios=(1.0, 0.20), hspace=0.10,
                          left=0.02, right=0.98, top=0.93, bottom=0.07)
    ax = fig.add_subplot(gs[0], projection="3d")
    ax.set_facecolor("#07090c")
    nx, ny, nz = (int(v) for v in F["shape"])
    # bench frame: plot (z, x, y) so the 75 cm optical axis lies
    # horizontally — the cell reads as the long thin duct it is
    legend = []
    for ci, comp in enumerate(GLOW_COMPONENTS):
        key = comp["key"]
        frac = (prep["sums"][key][j] / prep["peaks"][key]) \
            ** comp.get("tone", 1.0)
        n_dots = int(round(prep["pts"] * comp["gain"] * frac))
        cloud = glow_cloud(prep["em"][key][j], n_dots,
                           seed_offset=11 + ci)
        if len(cloud):
            cz, cx, cyy = cloud[:, 2], cloud[:, 0], cloud[:, 1]
            ax.scatter(cz, cx, cyy, s=26, c=comp["color"], alpha=0.028,
                       linewidths=0, depthshade=False)
            ax.scatter(cz, cx, cyy, s=5, c=comp["core"], alpha=0.085,
                       linewidths=0, depthshade=False)
        legend.append(f"{comp['label']} (display gain "
                      f"{comp['gain']:g}): {len(cloud)} samples")
    # the vessel in the bench frame: edges + the optical end windows
    bx, by, bz = nz, nx, ny
    for s, e in [((0, 0, 0), (bx, 0, 0)), ((0, 0, 0), (0, by, 0)),
                 ((0, 0, 0), (0, 0, bz)), ((bx, by, bz), (0, by, bz)),
                 ((bx, by, bz), (bx, 0, bz)), ((bx, by, bz), (bx, by, 0)),
                 ((bx, 0, 0), (bx, by, 0)), ((bx, 0, 0), (bx, 0, bz)),
                 ((0, by, 0), (bx, by, 0)), ((0, by, 0), (0, by, bz)),
                 ((0, 0, bz), (bx, 0, bz)), ((0, 0, bz), (0, by, bz))]:
        ax.plot(*zip(s, e), color="#55607a", lw=1.1)
    for xc in (0, bx):
        ax.plot([xc] * 5, [0, by, by, 0, 0], [0, 0, bz, bz, 0],
                color="#7c89a0", lw=1.4, alpha=0.8)
    ax.set_xlim(0, bx); ax.set_ylim(0, by); ax.set_zlim(0, bz)
    ex = F["extents_cm"]
    if aspect == "true":
        dims = (ex[2], ex[0], ex[1])
    elif aspect == "auto":          # √-compressed: still clearly a long
        dims = tuple(np.sqrt((ex[2], ex[0], ex[1])))   # duct, but the
    else:                           # cross-section stays visible
        dims = (1.0, 1.0, 1.0)
    ax.set_box_aspect(dims)
    ax.view_init(elev=15, azim=-80)
    ax.set_axis_off()
    ex = F["extents_cm"]
    for kk, txt in enumerate(legend):
        ax.text2D(0.01, 0.97 - 0.05 * kk, txt, transform=ax.transAxes,
                  fontsize=8, color=GLOW_COMPONENTS[kk]["color"])
    ne_j = F["fields"]["n_e"][j].ravel()
    spread = float(ne_j.max() / max(ne_j.min(), 1e-300))
    slabel = f"×{spread:.2f}" if spread < 1e4 else f"×{spread:.2e}"
    ax.text2D(0.01, 0.97 - 0.05 * len(legend),
              f"cell n_e spread {slabel}", transform=ax.transAxes,
              fontsize=8, color="#8a919c")
    ax.text2D(0.01, 0.015,
              f"cell {ex[0]:g} × {ex[1]:g} × {ex[2]:g} cm — optical "
              f"axis horizontal (end windows drawn); beam pumps "
              f"across the thin x side",
              transform=ax.transAxes, fontsize=8, color="#8a919c")
    st = _time_strip(fig, gs[1], F, j)
    st.set_facecolor("#0d1117")
    _dark_axes(fig.axes[-2:])
    asp = {"true": "true aspect",
           "unit": "axes independently normalized"}.get(
        aspect, f"√-compressed aspect (true {ex[0]:g}×{ex[1]:g}"
                f"×{ex[2]:g} cm)")
    fig.suptitle(
        f"REBUILD {F['meta'].get('mode', '')} — what the cell looks "
        f"like — t = {F['t_ns'][j]:.1f} ns — {asp}",
        fontsize=12, color="#d8dbe0")
    fig.text(0.5, 0.004,
             "brightness ∝ emission: KrF* fluorescence violet-white "
             "(real line is UV, linear tone) + n_e² continuum blue "
             "(gain 0.35, √ tone so the afterglow stays visible), each "
             "vs its run peak; trilinear DISPLAY interpolation between "
             "model cells — independent-reactor grid, not a PIC solve",
             ha="center", fontsize=6.8, color="#767c85")
    return fig


def bloom_png(path, radius=4, gain=0.8):
    """Camera-style bloom: max-composite a blurred copy over the frame."""
    from PIL import Image, ImageFilter
    img = Image.open(path).convert("RGB")
    blur = img.filter(ImageFilter.GaussianBlur(radius))
    arr = np.maximum(np.asarray(img, dtype=np.int16),
                     (np.asarray(blur, dtype=np.float32) * gain)
                     .astype(np.int16)).astype(np.uint8)
    Image.fromarray(arr).save(path)


def _frame_indices(n_t, frames):
    return np.unique(np.linspace(0, n_t - 1, min(frames, n_t)).round()
                     .astype(int))


def render_animation(F, specs, out_dir, frames=61, fps=12, aspect="unit",
                     gif=True, mp4=False, log=print, style="voxels",
                     pts=8000):
    """Frame PNGs -> a movie per style: 'glow' (the camera view — the
    headline movie.gif, true aspect, bloomed), 'particles'
    (movie_particles), 'voxels' (the field panels, movie_fields)."""
    _, plt = _mpl()
    out_dir = Path(out_dir)
    fdir = out_dir / {"glow": "frames_glow", "particles":
                      "frames_particles", "voxels": "frames_fields"}[style]
    stem = {"glow": "movie", "particles": "movie_particles",
            "voxels": "movie_fields"}[style]
    fdir.mkdir(parents=True, exist_ok=True)
    pcloud = prepare_particles(F, pts) if style == "particles" else None
    gprep = prepare_glow(F) if style == "glow" else None
    g_aspect = "auto" if style == "glow" and aspect == "unit" else aspect
    idx = _frame_indices(F["t_ns"].size, frames)
    fig = None
    paths = []
    for k, j in enumerate(idx):
        if style == "glow":
            fig = render_glow_frame(F, gprep, int(j), aspect=g_aspect,
                                    fig=fig)
        elif style == "particles":
            fig = render_particle_frame(F, pcloud, int(j), aspect=aspect,
                                        fig=fig)
        else:
            fig = render_frame(F, specs, int(j), aspect=aspect, fig=fig)
        p = fdir / f"frame_{k:04d}.png"
        fig.savefig(p, facecolor=fig.get_facecolor())
        if style == "glow":
            bloom_png(p)
        paths.append(p)
        if k % 10 == 0 or k == len(idx) - 1:
            log(f"  {style} frame {k + 1}/{len(idx)} "
                f"(t = {F['t_ns'][j]:.1f} ns)")
    plt.close(fig)
    made = {"frames": fdir}
    if gif:
        from PIL import Image
        imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE)
                for p in paths]
        gif_path = out_dir / f"{stem}.gif"
        imgs[0].save(gif_path, save_all=True, append_images=imgs[1:],
                     duration=int(1000 / fps), loop=0, optimize=True)
        made["gif"] = gif_path
    if mp4:
        ff = shutil.which("ffmpeg")
        if ff is None:
            log("  ffmpeg not found — skipping mp4")
        else:
            mp4_path = out_dir / f"{stem}.mp4"
            r = subprocess.run(
                [ff, "-y", "-loglevel", "error", "-framerate", str(fps),
                 "-i", str(fdir / "frame_%04d.png"),
                 "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                 "-pix_fmt", "yuv420p", str(mp4_path)],
                capture_output=True, text=True)
            if r.returncode == 0:
                made["mp4"] = mp4_path
            else:
                log(f"  ffmpeg failed ({r.stderr.strip()[:200]}) — no mp4")
    return made


def render_contact_sheet(F, specs, path, n=12, aspect="unit"):
    """n_e snapshots (the single most telling field) on one page."""
    mpl, plt = _mpl()
    idx = _frame_indices(F["t_ns"].size, n)
    rows, cols = 3, 4
    fig = plt.figure(figsize=(13.0, 9.5), dpi=110)
    spec = specs["n_e"]
    norm, cmap = _norm_of(mpl, spec), _cmap_of(mpl, spec)
    for k, j in enumerate(idx[: rows * cols]):
        ax = fig.add_subplot(rows, cols, k + 1, projection="3d")
        vals = F["fields"]["n_e"][j]
        nvals = _normed(norm, vals)
        rgba = cmap(nvals)
        rgba[..., 3] = 0.30 + 0.65 * nvals
        ax.voxels(np.ones(vals.shape, bool), facecolors=rgba,
                  edgecolor=(0, 0, 0, 0.15), linewidth=0.15)
        ax.set_box_aspect((1, 1, 1) if aspect == "unit"
                          else tuple(F["extents_cm"]))
        ax.set_axis_off()
        ax.set_title(f"{F['t_ns'][j]:.0f} ns", fontsize=8)
    from matplotlib.cm import ScalarMappable
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                      ax=fig.axes, shrink=0.5, pad=0.02)
    cb.set_label("n_e m^-3", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("n_e over the pulse — REBUILD 3-D device grid "
                 "(axes normalized)", fontsize=12)
    fig.text(0.5, 0.01, _CAPTION, ha="center", fontsize=7, color="#777")
    fig.savefig(path)
    plt.close(fig)
    return path


# ------------------------------------------------------------------- html
def build_viewer_payload(F, specs, frames=61) -> dict:
    idx = _frame_indices(F["t_ns"].size, frames)
    fields = {}
    for name, vals in F["fields"].items():
        s = specs[name]
        fields[name] = {
            "scale": s["scale"], "vmin": s["vmin"], "vmax": s["vmax"],
            "cmap": s["cmap"], "unit": s["unit"],
            "frames": [quantize(vals[j].ravel(), s["vmin"], s["vmax"],
                                s["scale"]).tolist() for j in idx],
        }
    if "dep_shape" in F:
        s = specs["deposition"]
        fields["deposition"] = {
            "scale": s["scale"], "vmin": 0.0, "vmax": 1.0,
            "cmap": s["cmap"], "unit": s["unit"], "static": True,
            "frames": [quantize(F["dep_shape"].ravel(), 0.0, 1.0,
                                "lin").tolist()],
        }
    n_cells = int(np.prod(F["shape"]))
    species = []
    for s_idx, (name, color, label) in enumerate(PARTICLE_SPECIES):
        fname = F["meta"]["upper"] if name is None else name
        vals = F["fields"][fname].reshape(len(F["t_ns"]), n_cells)
        per = vals * F["cell_vols"]
        species.append({"field": fname, "color": color, "label": label,
                        "sidx": s_idx,
                        "pnorm": float(per.sum(axis=1).max())})
    return {
        "shape": [int(v) for v in F["shape"]],
        "extents_cm": [float(v) for v in F["extents_cm"]],
        "t_ns": [float(F["t_ns"][j]) for j in idx],
        "pout_W": [float(F["P_out_W"][j]) for j in idx],
        "g0avg_m1": [float(F["g0avg_m1"][j]) for j in idx],
        "meta": F["meta"], "fields": fields,
        "particles": {"npts": 3000,
                      "cellvol": [float(v) for v in F["cell_vols"]],
                      "species": species},
    }


def write_html(F, specs, path, frames=61) -> Path:
    payload = build_viewer_payload(F, specs, frames)
    data_json = json.dumps(payload, separators=(",", ":"))
    title = (f"REBUILD {payload['meta'].get('mode', '')} 3-D grid "
             f"{'x'.join(str(v) for v in payload['shape'])}")
    html = _HTML_TEMPLATE.replace("__TITLE__", title) \
                         .replace("__CAPTION__", _CAPTION) \
                         .replace("__DATA__", data_json)
    Path(path).write_text(html)
    return Path(path)


_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  body{margin:0;background:#14161a;color:#d8dbe0;
       font:13px/1.45 -apple-system,'Segoe UI',sans-serif}
  .wrap{max-width:1080px;margin:0 auto;padding:14px}
  h1{font-size:17px;margin:4px 0 2px}
  .sub{color:#8a8f98;font-size:12px;margin:0 0 10px}
  .row{display:flex;gap:14px;flex-wrap:wrap;align-items:flex-start}
  canvas{background:#0d0f12;border:1px solid #2a2e35;border-radius:6px}
  .ctl{display:flex;flex-direction:column;gap:9px;min-width:230px;flex:1}
  .ctl label{display:flex;align-items:center;gap:8px;font-size:12px;
             color:#aab0b8}
  .ctl input[type=range]{flex:1}
  select,button{background:#1e2228;color:#d8dbe0;border:1px solid #343a43;
                border-radius:4px;padding:4px 8px;font-size:12px}
  button{cursor:pointer}
  .readout{font-family:ui-monospace,Menlo,monospace;font-size:11px;
           color:#9aa0a8;white-space:pre-wrap}
  .cap{color:#767c85;font-size:11px;margin-top:10px}
</style>
</head>
<body>
<div class="wrap">
  <h1>__TITLE__</h1>
  <p class="sub" id="regime"></p>
  <div class="row">
    <canvas id="view" width="720" height="520"></canvas>
    <div class="ctl">
      <label>view
        <label><input type="radio" name="mode" value="glow" checked>
          glow</label>
        <label><input type="radio" name="mode" value="particles">
          particles</label>
        <label><input type="radio" name="mode" value="voxels">
          voxels</label>
      </label>
      <label>field (voxel view)
        <select id="field"></select>
      </label>
      <label>t <input type="range" id="time" min="0" value="0" step="1">
        <span id="tlabel" class="readout"></span></label>
      <label>yaw <input type="range" id="yaw" min="-180" max="180"
        value="-35" step="1"></label>
      <label><input type="checkbox" id="aspect"> true aspect
        (device is long and thin)</label>
      <button id="play">play</button>
      <canvas id="cbar" width="230" height="46"></canvas>
      <div class="readout" id="range"></div>
      <div class="readout" id="strip"></div>
    </div>
  </div>
  <p class="cap">__CAPTION__.
  Colors are quantized to 8 bits for this viewer — exact values live in
  viz_data.npz alongside this file.</p>
</div>
<script id="vizdata" type="application/json">__DATA__</script>
<script>
"use strict";
const D = JSON.parse(document.getElementById("vizdata").textContent);
const [NX, NY, NZ] = D.shape, NC = NX * NY * NZ;
const PD = D.particles;   // dot-cloud bookkeeping (same math as the tool)
const view = document.getElementById("view"), ctx = view.getContext("2d");
const fieldSel = document.getElementById("field");
const timeSl = document.getElementById("time");
const yawSl = document.getElementById("yaw");
const aspectCb = document.getElementById("aspect");
const playBt = document.getElementById("play");
document.getElementById("regime").textContent =
  (D.meta && D.meta.regime) ? D.meta.regime : "";
for (const name of Object.keys(D.fields)) {
  const o = document.createElement("option");
  o.value = name; o.textContent = name; fieldSel.appendChild(o);
}
timeSl.max = String(D.t_ns.length - 1);
// --- colormaps (compact LUTs, lerped) --------------------------------
const LUTS = {
  viridis: [[68,1,84],[71,44,122],[59,81,139],[44,113,142],[33,144,141],
            [39,173,129],[92,200,99],[170,220,50],[253,231,37]],
  magma:   [[0,0,4],[40,11,84],[101,21,110],[159,42,99],[212,72,66],
            [245,125,21],[250,193,39],[252,253,191]],
  diverging:[[33,102,172],[103,169,207],[209,229,240],[247,247,247],
             [253,219,199],[239,138,98],[178,24,43]],
};
function lut(name, f) {
  const L = LUTS[name] || LUTS.viridis;
  const x = Math.min(Math.max(f, 0), 1) * (L.length - 1);
  const i = Math.min(Math.floor(x), L.length - 2), u = x - i;
  const a = L[i], b = L[i + 1];
  return [a[0]+(b[0]-a[0])*u, a[1]+(b[1]-a[1])*u, a[2]+(b[2]-a[2])*u];
}
// --- geometry --------------------------------------------------------
// Orbit camera: yaw about world z, then elevation ELEV above horizontal.
// Camera basis (yaw'd frame): view-to-viewer v = (0,-ce,se),
// up u = (0,se,ce), right r = x. project() returns
// [screen-x, screen-up, farness]; faces are visible when the outward
// normal's farness component is negative (it points at the viewer),
// and the painter draws descending farness (far first).
const ELEV = 28 * Math.PI / 180;
function boxDims() {
  if (aspectCb.checked) {
    const m = Math.max(...D.extents_cm);
    return D.extents_cm.map(e => e / m);
  }
  return [1, 1, 1];
}
function project(p, cy, sy, ce, se) {
  const x = p[0]*cy - p[1]*sy, y = p[0]*sy + p[1]*cy, z = p[2];
  return [x, y*se + z*ce, y*ce - z*se];   // [sx, up, farness]
}
function mulberry32(a) {
  a = (a + 0x6D2B79F5) >>> 0;
  let t = a;
  t = Math.imul(t ^ (t >>> 15), t | 1) >>> 0;
  t = (t ^ ((t + (Math.imul(t ^ (t >>> 7), t | 61) >>> 0)) >>> 0)) >>> 0;
  return (((t ^ (t >>> 14)) >>> 0) / 4294967296);
}
function samplePos(c, k, sidx, axis) {
  return mulberry32((c * 9781 + k * 2621 + sidx * 131071
                     + axis * 40503) >>> 0);
}
function mode() {
  return document.querySelector('input[name="mode"]:checked').value;
}
function geom() {
  const yaw = +yawSl.value * Math.PI / 180;
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  const ct = Math.cos(ELEV), st = Math.sin(ELEV);
  const [bx, by, bz] = boxDims();
  const R = Math.hypot(bx, by, bz) / 2;
  const S = Math.min(view.width, view.height) / (2.35 * R);
  const cxp = view.width / 2, cyp = view.height / 2 + 12;
  const P = (x, y, z) => {
    const pr = project([x, y, z], cy, sy, ct, st);
    return [cxp + S * pr[0], cyp - S * pr[1]];
  };
  return {cy, sy, ct, st, bx, by, bz, dx: bx / NX, dy: by / NY,
          dz: bz / NZ, P};
}
function cellOrder(g) {
  const order = [];
  for (let i = 0; i < NX; i++) for (let j = 0; j < NY; j++)
    for (let k = 0; k < NZ; k++) {
      const c = [(i + 0.5) * g.dx - g.bx / 2, (j + 0.5) * g.dy - g.by / 2,
                 (k + 0.5) * g.dz - g.bz / 2];
      order.push([project(c, g.cy, g.sy, g.ct, g.st)[2], i, j, k]);
    }
  order.sort((a, b) => b[0] - a[0]);      // far first (painter)
  return order;
}
function boxEdges(g) {
  const [bx, by, bz] = [g.bx, g.by, g.bz];
  const C = [[0,0,0],[bx,0,0],[0,by,0],[0,0,bz],[bx,by,0],[bx,0,bz],
             [0,by,bz],[bx,by,bz]].map(p =>
               g.P(p[0]-bx/2, p[1]-by/2, p[2]-bz/2));
  const E = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],
             [4,7],[5,7],[6,7]];
  ctx.strokeStyle = "#3a4250"; ctx.lineWidth = 1;
  for (const [a, b] of E) {
    ctx.beginPath(); ctx.moveTo(C[a][0], C[a][1]);
    ctx.lineTo(C[b][0], C[b][1]); ctx.stroke();
  }
}
function dequant(f, spec) {
  if (spec.scale === "log") {
    const lo = Math.log10(spec.vmin), hi = Math.log10(spec.vmax);
    return Math.pow(10, lo + f * (hi - lo));
  }
  return spec.vmin + f * (spec.vmax - spec.vmin);
}
function drawFrame() {
  const m = mode();
  if (m === "glow") { drawGlow(); return; }
  if (m === "particles") { drawParticles(); return; }
  drawVoxels();
}
// trilinear sample of a field at cell-index coords (values live at
// cell centers; clamped at the walls)
function trilerp(q, fld, x, y, z) {
  const fx = Math.min(Math.max(x - 0.5, 0), NX - 1),
        fy = Math.min(Math.max(y - 0.5, 0), NY - 1),
        fz = Math.min(Math.max(z - 0.5, 0), NZ - 1);
  const i0 = Math.min(Math.floor(fx), NX - 1), i1 = Math.min(i0 + 1, NX - 1),
        j0 = Math.min(Math.floor(fy), NY - 1), j1 = Math.min(j0 + 1, NY - 1),
        k0 = Math.min(Math.floor(fz), NZ - 1), k1 = Math.min(k0 + 1, NZ - 1);
  const tx = fx - i0, ty = fy - j0, tz = fz - k0;
  const V = (i, j, k) => dequant(q[(i * NY + j) * NZ + k] / 255, fld);
  const c00 = V(i0,j0,k0)*(1-tx) + V(i1,j0,k0)*tx;
  const c01 = V(i0,j0,k1)*(1-tx) + V(i1,j0,k1)*tx;
  const c10 = V(i0,j1,k0)*(1-tx) + V(i1,j1,k0)*tx;
  const c11 = V(i0,j1,k1)*(1-tx) + V(i1,j1,k1)*tx;
  return (c00*(1-ty) + c10*ty)*(1-tz) + (c01*(1-ty) + c11*ty)*tz;
}
const GLOWC = [   // tone: sqrt on the continuum keeps the afterglow lit
  {field: null, sq: false, rgb: [217,200,255], gain: 1.0, tone: 1.0},
  {field: "n_e", sq: true, rgb: [59,111,212], gain: 0.35, tone: 0.5},
];
GLOWC[0].field = (D.meta && D.meta.upper) ? D.meta.upper : "n_e";
const SPRITES_PER_CELL = 26;
function drawGlow() {
  const ti = +timeSl.value;
  const g = geom();
  ctx.clearRect(0, 0, view.width, view.height);
  boxEdges(g);
  const R = Math.hypot(g.bx, g.by, g.bz) / 2;
  const S = Math.min(view.width, view.height) / (2.35 * R);
  const rp = Math.max(5, 1.6 * S * Math.cbrt(g.dx * g.dy * g.dz));
  ctx.globalCompositeOperation = "lighter";
  for (let ci = 0; ci < GLOWC.length; ci++) {
    const comp = GLOWC[ci];
    const fld = D.fields[comp.field];
    if (!fld) continue;
    const q = fld.frames[fld.static ? 0 : ti];
    const peak = comp.sq ? fld.vmax * fld.vmax : fld.vmax;
    for (let c = 0; c < NC; c++) {
      const i = Math.floor(c / (NY * NZ)),
            jj = Math.floor(c / NZ) % NY, kk = c % NZ;
      for (let s = 0; s < SPRITES_PER_CELL; s++) {
        const ux = samplePos(c, s, 21 + ci, 0),
              uy = samplePos(c, s, 21 + ci, 1),
              uz = samplePos(c, s, 21 + ci, 2);
        let val = trilerp(q, fld, i + ux, jj + uy, kk + uz);
        if (comp.sq) val = val * val;
        const a = comp.gain * Math.pow(val / peak, comp.tone);
        if (a < 0.006) continue;
        const [px, py] = g.P((i + ux) * g.dx - g.bx / 2,
                             (jj + uy) * g.dy - g.by / 2,
                             (kk + uz) * g.dz - g.bz / 2);
        const gr = ctx.createRadialGradient(px, py, 0, px, py, rp);
        const al = Math.min(0.55 * a, 0.5);
        gr.addColorStop(0, `rgba(${comp.rgb[0]},${comp.rgb[1]},` +
                           `${comp.rgb[2]},${al})`);
        gr.addColorStop(1, `rgba(${comp.rgb[0]},${comp.rgb[1]},` +
                           `${comp.rgb[2]},0)`);
        ctx.fillStyle = gr;
        ctx.fillRect(px - rp, py - rp, 2 * rp, 2 * rp);
      }
    }
  }
  ctx.globalCompositeOperation = "source-over";
  ctx.fillStyle = "#8a8f98"; ctx.font = "11px sans-serif";
  ctx.fillText("brightness ∝ emission: upper-level fluorescence " +
    "(violet-white) + n_e² continuum (blue, gain 0.35); trilinear " +
    "display interpolation", 10, view.height - 26);
  ctx.fillText(`x ${D.extents_cm[0].toFixed(2)} cm  y ` +
    `${D.extents_cm[1].toFixed(2)} cm  z ${D.extents_cm[2].toFixed(1)} cm`,
    10, view.height - 10);
  updateSide("glow", null, ti);
}
function drawParticles() {
  const ti = +timeSl.value;
  const g = geom();
  ctx.clearRect(0, 0, view.width, view.height);
  boxEdges(g);
  const order = cellOrder(g);
  const totals = [];
  for (const sp of PD.species) {
    const fld = D.fields[sp.field];
    const q = fld.frames[ti];
    let total = 0;
    ctx.fillStyle = sp.color;
    ctx.globalAlpha = 0.75;
    for (const [, i, j, k] of order) {
      const c = (i * NY + j) * NZ + k;
      const val = dequant(q[c] / 255, fld);
      const cnt = Math.round(PD.npts * val * PD.cellvol[c] / sp.pnorm);
      total += cnt;
      const x0 = i * g.dx - g.bx / 2, y0 = j * g.dy - g.by / 2,
            z0 = k * g.dz - g.bz / 2;
      for (let m = 0; m < cnt; m++) {
        const [px, py] = g.P(x0 + samplePos(c, m, sp.sidx, 0) * g.dx,
                             y0 + samplePos(c, m, sp.sidx, 1) * g.dy,
                             z0 + samplePos(c, m, sp.sidx, 2) * g.dz);
        ctx.fillRect(px - 1, py - 1, 2, 2);
      }
    }
    totals.push(total);
  }
  ctx.globalAlpha = 1.0;
  ctx.font = "12px sans-serif";
  PD.species.forEach((sp, s) => {
    ctx.fillStyle = sp.color;
    const w = sp.pnorm / PD.npts;
    ctx.fillText(`${sp.label} (${sp.field}): ${totals[s]} dots x ` +
                 `${w.toExponential(2)} real`, 10, 20 + 16 * s);
  });
  ctx.fillStyle = "#8a8f98"; ctx.font = "11px sans-serif";
  ctx.fillText("dot count ∝ particles per cell; positions = fixed " +
               "sampling lattice, not trajectories", 10,
               view.height - 26);
  ctx.fillText(`x ${D.extents_cm[0].toFixed(2)} cm  y ` +
    `${D.extents_cm[1].toFixed(2)} cm  z ${D.extents_cm[2].toFixed(1)} cm`,
    10, view.height - 10);
  updateSide("particles", null, ti);
}
function drawVoxels() {
  const name = fieldSel.value, fld = D.fields[name];
  const ti = fld.static ? 0 : +timeSl.value;
  const q = fld.frames[ti];
  const g = geom();
  const cy = g.cy, sy = g.sy, ct = g.ct, st = g.st;
  const dx = g.dx, dy = g.dy, dz = g.dz;
  const bx = g.bx, by = g.by, bz = g.bz;
  ctx.clearRect(0, 0, view.width, view.height);
  const order = cellOrder(g);
  const P = (x, y, z) => g.P(x, y, z);
  const FACES = [   // [corner offsets x4, brightness, outward normal]
    [[[0,0,1],[1,0,1],[1,1,1],[0,1,1]], 1.00, [0,0,1]],   // top (z+)
    [[[0,0,0],[1,0,0],[1,0,1],[0,0,1]], 0.80, [0,-1,0]],  // y-
    [[[1,0,0],[1,1,0],[1,1,1],[1,0,1]], 0.66, [1,0,0]],   // x+
    [[[0,1,0],[0,0,0],[0,0,1],[0,1,1]], 0.66, [-1,0,0]],  // x-
    [[[1,1,0],[0,1,0],[0,1,1],[1,1,1]], 0.80, [0,1,0]],   // y+
  ];
  for (const [, i, j, k] of order) {
    const f = q[(i*NY + j)*NZ + k] / 255;
    const rgb = lut(fld.cmap, f);
    const x0 = i*dx - bx/2, y0 = j*dy - by/2, z0 = k*dz - bz/2;
    for (const [corners, br, n] of FACES) {
      // visible if the outward normal faces the viewer (depth decreases)
      const nr = project(n, cy, sy, ct, st);
      if (nr[2] >= 0) continue;
      ctx.beginPath();
      for (let c = 0; c < 4; c++) {
        const [px, py] = P(x0 + corners[c][0]*dx, y0 + corners[c][1]*dy,
                           z0 + corners[c][2]*dz);
        if (c === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.closePath();
      const al = 0.35 + 0.6 * f;
      ctx.fillStyle = `rgba(${rgb[0]*br|0},${rgb[1]*br|0},` +
                      `${rgb[2]*br|0},${al})`;
      ctx.fill();
      ctx.strokeStyle = "rgba(0,0,0,0.28)"; ctx.lineWidth = 0.5;
      ctx.stroke();
    }
  }
  // labels
  ctx.fillStyle = "#8a8f98"; ctx.font = "11px sans-serif";
  ctx.fillText(`x ${D.extents_cm[0].toFixed(2)} cm  y ` +
    `${D.extents_cm[1].toFixed(2)} cm  z ${D.extents_cm[2].toFixed(1)} cm`,
    10, view.height - 10);
  updateSide(name, fld, ti);
}
function fmt(v, scale) {
  return (scale === "log" || Math.abs(v) >= 1e4 ||
          (v !== 0 && Math.abs(v) < 1e-2)) ? v.toExponential(2)
                                           : v.toPrecision(4);
}
function updateSide(name, fld, ti) {
  const t = D.t_ns[+timeSl.value];
  document.getElementById("tlabel").textContent = t.toFixed(1) + " ns";
  const cb = document.getElementById("cbar"), c2 = cb.getContext("2d");
  c2.clearRect(0, 0, cb.width, cb.height);
  if (!fld) {              // glow / particle views: no colorbar
    const ne = D.fields["n_e"];
    const q = ne.frames[+timeSl.value];
    let qlo = 255, qhi = 0;
    for (const v of q) { if (v < qlo) qlo = v; if (v > qhi) qhi = v; }
    const spread = dequant(qhi / 255, ne) / dequant(qlo / 255, ne);
    document.getElementById("range").textContent =
      (name === "glow"
       ? "camera view — what the cell emits (fluorescence + faint " +
         "n_e² continuum); exact fields in the npz"
       : PD.species.map(sp => `${sp.label}: 1 dot ≈ ` +
           `${(sp.pnorm / PD.npts).toExponential(2)} real`).join("   ")
      ) + `   cell n_e spread ×${spread.toPrecision(3)}`;
  } else {
    for (let x = 0; x < cb.width; x++) {
      const rgb = lut(fld.cmap, x / (cb.width - 1));
      c2.fillStyle = `rgb(${rgb[0]|0},${rgb[1]|0},${rgb[2]|0})`;
      c2.fillRect(x, 6, 1, 16);
    }
    c2.fillStyle = "#aab0b8"; c2.font = "10px sans-serif";
    c2.fillText(fmt(fld.vmin, fld.scale), 0, 36);
    const hi = fmt(fld.vmax, fld.scale);
    c2.fillText(hi, cb.width - c2.measureText(hi).width, 36);
    document.getElementById("range").textContent =
      `${name} ${fld.unit} — ${fld.scale} scale` +
      (fld.static ? " — static (no time dependence)" : "");
  }
  document.getElementById("strip").textContent =
    `P_out ${fmt(D.pout_W[+timeSl.value], "log")} W    ` +
    `<g0> ${fmt(D.g0avg_m1[+timeSl.value], "lin")} 1/m`;
}
let timer = null;
playBt.addEventListener("click", () => {
  if (timer) { clearInterval(timer); timer = null;
               playBt.textContent = "play"; return; }
  playBt.textContent = "pause";
  timer = setInterval(() => {
    timeSl.value = String((+timeSl.value + 1) % D.t_ns.length);
    drawFrame();
  }, 90);
});
for (const el of [fieldSel, timeSl, yawSl, aspectCb])
  el.addEventListener("input", drawFrame);
for (const el of document.querySelectorAll('input[name="mode"]'))
  el.addEventListener("input", drawFrame);
drawFrame();
</script>
</body>
</html>
"""


# ------------------------------------------------------------------- main
def main(argv=None):
    from REBUILD.run import PRESETS
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--preset", choices=sorted(PRESETS))
    src.add_argument("--npz", metavar="PATH",
                     help="re-render from a saved viz_data.npz (no run)")
    ap.add_argument("--out", metavar="DIR", default=None,
                    help="output dir (default rebuild_out/viz_<preset>/)")
    ap.add_argument("--t-end", type=float, default=None, metavar="NS")
    ap.add_argument("--n-save", type=int, default=121)
    ap.add_argument("--tier", choices=("loki", "kinetic"), default=None)
    ap.add_argument("--frames", type=int, default=61,
                    help="animation/viewer frame count (default 61)")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--aspect", choices=("unit", "true"), default="unit")
    ap.add_argument("--style",
                    choices=("all", "glow", "particles", "voxels",
                             "both"),
                    default="all",
                    help="animation style(s): glow = the camera view "
                         "(movie.gif), particles, voxels (field "
                         "panels, movie_fields.gif); all = every one "
                         "(default); both = particles+voxels")
    ap.add_argument("--pts", type=int, default=8000,
                    help="particle-view dots per species at peak "
                         "(macro-weight adjusts; default 8000)")
    ap.add_argument("--mp4", action="store_true",
                    help="also encode movie.mp4 (needs ffmpeg)")
    ap.add_argument("--no-gif", action="store_true")
    ap.add_argument("--no-movie", action="store_true",
                    help="skip frames/gif/mp4 (fast npz -> viewer/sheet)")
    ap.add_argument("--no-html", action="store_true")
    args = ap.parse_args(argv)

    if args.npz:
        F = load_npz(args.npz)
        out = Path(args.out) if args.out else Path(args.npz).parent
        print(f"[viz] loaded {args.npz}")
    else:
        dev = build_device(args.preset, t_end_ns=args.t_end,
                           tier=args.tier)
        print(f"[viz] running preset {args.preset} "
              f"(grid {dev.shape[0]}x{dev.shape[1]}x{dev.shape[2]}, "
              f"{args.n_save} saves) ...")
        res = dev.run(n_save=args.n_save, dense_cells=True)
        print(res["regime"])
        F = extract_fields(dev, res)
        F["meta"]["preset"] = args.preset
        out = Path(args.out) if args.out else \
            ROOT / "rebuild_out" / f"viz_{args.preset}"
        out.mkdir(parents=True, exist_ok=True)
        save_npz(out / "viz_data.npz", F)
        print(f"[viz] data of record -> {out / 'viz_data.npz'}")
    out.mkdir(parents=True, exist_ok=True)
    specs = field_specs(F)

    sheet = render_contact_sheet(F, specs, out / "contact_sheet.png",
                                 aspect=args.aspect)
    print(f"[viz] contact sheet -> {sheet}")
    if not args.no_movie:
        styles = {"all": ("glow", "particles", "voxels"),
                  "both": ("particles", "voxels")}.get(args.style,
                                                      (args.style,))
        for style in styles:
            made = render_animation(F, specs, out, frames=args.frames,
                                    fps=args.fps, aspect=args.aspect,
                                    gif=not args.no_gif, mp4=args.mp4,
                                    style=style, pts=args.pts)
            for k in ("gif", "mp4"):
                if k in made:
                    print(f"[viz] {style} {k} -> {made[k]}")
    if not args.no_html:
        viewer = write_html(F, specs, out / "viewer.html",
                            frames=args.frames)
        print(f"[viz] interactive viewer -> {viewer} (open in a browser)")
    print("[viz] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
