#!/usr/bin/env python3
"""REAL macro-particles, visible: snapshot a PIC run's particle state
over time and render it.

Unlike REBUILD's density-sampled particle view (a rate-equation model has
no particles to show), every dot here IS a macro-particle of the ES
PIC loop — position from the leapfrog, color from its v_z — captured by
reading the public particle arrays between steps. Nothing in MODULES/ is
modified or imported FROM here by the engine.

The stock demo is the package's own headline-gate configuration made
visible: a cold electron plasma in a neutralising jelly, quiet-started
with a sinusoidal displacement (``ParticleSpecies.cold_wave``), sloshing
at the cold-plasma frequency in a fully periodic box. The mode amplitude
trace below the dots is the same signal the frequency gate measures.

Scale honesty: this is the code's honest regime — a reduced-density
(1e14 m^-3), centimetre, sub-relativistic box — NOT the laser chamber;
the resolution budget in ../PIC_MCC_3D.txt section 3 is why.

Outputs (default ``<root>/picmcc_out/viz_langmuir/``): particles_data.npz
(positions/velocities per snapshot — the data of record), frame PNGs →
movie.gif / movie.mp4 (--mp4), and viewer.html (self-contained canvas
viewer; positions quantized to uint8 of the box — exact values in the
npz).

Usage (from the Kinetic_Modeling root):
    python3 SPECK/PIC_MCC_3D/tools/visualize_particles.py
    python3 SPECK/PIC_MCC_3D/tools/visualize_particles.py --npz SPECK/picmcc_out/viz_langmuir/particles_data.npz
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "PIC_MCC_3D" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("PIC_MCC_3D/MODULES not found above this script")


# --------------------------------------------------------------------- run
def run_langmuir(nx=16, per_axis=12, periods=3.0, wpe_dt=0.05,
                 amp_frac=0.05, every=4, n_e_m3=1.0e14, box_m=1.0e-2):
    """The cold Langmuir slosh, snapshotting REAL macro-particle state.

    amp_frac is larger than the frequency gate's 0.004 so the motion is
    visible; the gate's linear-regime number lives in tests/test_es_sim.py.
    """
    from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
    from PIC_MCC_3D.MODULES.poisson3d import PERIODIC
    from PIC_MCC_3D.MODULES.species3d import ParticleSpecies, Species, \
        SpeciesSet
    from PIC_MCC_3D.MODULES.units import plasma_frequency

    L = float(box_m)
    h = L / nx
    wpe = plasma_frequency(n_e_m3)
    dt = wpe_dt / wpe
    pop = ParticleSpecies.cold_wave(Species.electron(), per_axis ** 3,
                                    [0, 0, 0], [L, L, L], n_e_m3,
                                    amplitude=amp_frac * L, mode=1, axis=2)
    sim = ESSim(nx, nx, nx, h, SpeciesSet([pop]), dt, bc=PERIODIC,
                background_charge_density=uniform_background_for(pop,
                                                                 L ** 3))
    sim.set_mode(axis=2, mode=1)
    sim.prepare()
    n_steps = int(round(periods * 2.0 * np.pi / wpe / dt))
    t, pos, vz, ke, mode_re = [], [], [], [], []

    def snap():
        n = pop.n
        t.append(sim.time)
        pos.append(np.stack([pop.x[:n], pop.y[:n], pop.z[:n]],
                            axis=1).astype(np.float32))
        vz.append(pop.vz[:n].astype(np.float32))
        ke.append(pop.kinetic_energy_J())
        sim.solve_fields()
        mode_re.append(sim.density_mode().real)

    for k in range(n_steps):
        if k % every == 0:
            snap()
        sim.step()
    snap()
    weight = pop.real_count() / max(pop.n, 1)
    return {
        "t_s": np.asarray(t), "pos": np.stack(pos), "vz": np.stack(vz),
        "ke_J": np.asarray(ke), "mode_re": np.asarray(mode_re),
        "box_m": L,
        "meta": {
            "demo": "langmuir", "n_macro": int(pop.n),
            "weight_real": float(weight), "n_e_m3": n_e_m3,
            "wpe_rad_s": float(wpe), "dt_s": float(dt),
            "amp_frac": float(amp_frac), "grid": [nx, nx, nx],
            "periods": float(periods),
            "note": ("REAL macro-particles from the ES PIC loop "
                     "(cold Langmuir slosh, periodic, neutralising "
                     "jelly); reduced-density cm box — the code's "
                     "honest regime, not the laser chamber"),
        },
    }


def twostream_growth_rate(k, v0, wb):
    """Closed-form cold symmetric two-stream growth rate [rad/s].

    Two equal cold beams (each with beam plasma frequency wb) drifting
    at ±v0: the dispersion (u²−κ²)² = 2(u²+κ²) with u = ω/wb,
    κ = k·v0/wb gives u² = (κ²+1) − sqrt(4κ²+1) on the unstable branch
    (κ < √2); γ = wb·sqrt(−u²). Maximum γ = wb/2 at κ = √3/2.
    """
    kappa = k * v0 / wb
    s = (kappa ** 2 + 1.0) - np.sqrt(4.0 * kappa ** 2 + 1.0)
    return wb * np.sqrt(max(-s, 0.0))


def run_twostream(nx=16, per_axis=10, periods=10.0, wpe_dt=0.05,
                  amp_frac=0.01, every=10, n_e_m3=1.0e14, box_m=1.0e-2,
                  vth_frac=0.02, seed=20260812):
    """Two-stream instability: REAL self-generated bunching.

    Two equal electron beams counter-stream at ±v0 along z, with v0
    chosen so the box's mode-1 wavenumber sits at the fastest-growing
    κ = k·v0/wb = √3/2 (γ = wb/2, closed form above — a falsifiable
    target, nothing tuned). Opposite-phase displacement seeds put the
    start on the growing eigenmode; the macro-particles then bunch,
    trap, and roll into phase-space vortices — filamentation-class
    structure from the dynamics, not from an imposed profile.

    ``vth_frac``: a small Gaussian v_z spread (default 2 % of v0, from a
    SEEDED generator — runs stay deterministic). Without it the regular
    transverse lattice is exactly degenerate in (z, v_z) — the mode-1
    field is x,y-independent, so all particles in a z-plane share one
    phase-space point forever and the trapping vortex cannot fill in.
    The spread is far below v0, so the cold growth rate still applies
    (and is still gated).
    """
    from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
    from PIC_MCC_3D.MODULES.poisson3d import PERIODIC
    from PIC_MCC_3D.MODULES.species3d import ParticleSpecies, Species, \
        SpeciesSet
    from PIC_MCC_3D.MODULES.units import plasma_frequency

    L = float(box_m)
    h = L / nx
    wpe = plasma_frequency(n_e_m3)              # total (both beams)
    wb = wpe / np.sqrt(2.0)                     # per beam (n/2 each)
    k = 2.0 * np.pi / L                         # seeded mode m=1
    v0 = (np.sqrt(3.0) / 2.0) * wb / k          # fastest-growing kappa
    dt = wpe_dt / wpe
    rng = np.random.default_rng(seed)
    beams = []
    for sgn in (+1.0, -1.0):
        pop = ParticleSpecies.cold_wave(
            Species.electron(), per_axis ** 3, [0, 0, 0], [L, L, L],
            n_e_m3 / 2.0, amplitude=sgn * amp_frac * L, mode=1, axis=2)
        pop.vz[:pop.n] += sgn * v0 + \
            vth_frac * v0 * rng.standard_normal(pop.n)
        beams.append(pop)
    bg = sum(uniform_background_for(p, L ** 3) for p in beams)
    sim = ESSim(nx, nx, nx, h, SpeciesSet(beams), dt, bc=PERIODIC,
                background_charge_density=bg)
    sim.set_mode(axis=2, mode=1)
    sim.prepare()
    n_steps = int(round(periods * 2.0 * np.pi / wpe / dt))
    n1 = beams[0].n
    beam_id = np.concatenate([np.zeros(n1, np.int8),
                              np.ones(beams[1].n, np.int8)])
    t, pos, vz, ke, mode_re, mode_im, bunch = [], [], [], [], [], [], []

    def snap():
        t.append(sim.time)
        pos.append(np.concatenate(
            [np.stack([p.x[:p.n], p.y[:p.n], p.z[:p.n]], axis=1)
             for p in beams]).astype(np.float32))
        vz.append(np.concatenate([p.vz[:p.n]
                                  for p in beams]).astype(np.float32))
        ke.append(sum(p.kinetic_energy_J() for p in beams))
        # the textbook bunching parameter per beam, |<exp(-ikz)>|: the
        # clean monotone growth signal (the TOTAL-charge mode starts
        # cancelled — the seeds are anti-phase — and beats at k*v0)
        bunch.append([float(np.abs(np.mean(
            np.exp(-1j * k * p.z[:p.n])))) for p in beams])
        sim.solve_fields()
        mm = sim.density_mode()
        mode_re.append(mm.real)
        mode_im.append(mm.imag)

    for kk in range(n_steps):
        if kk % every == 0:
            snap()
        sim.step()
    snap()
    weight = beams[0].real_count() / max(beams[0].n, 1)
    return {
        "t_s": np.asarray(t), "pos": np.stack(pos), "vz": np.stack(vz),
        "ke_J": np.asarray(ke), "mode_re": np.asarray(mode_re),
        "mode_im": np.asarray(mode_im), "bunch": np.asarray(bunch),
        "beam": beam_id, "box_m": L,
        "meta": {
            "demo": "twostream", "n_macro": int(beam_id.size),
            "weight_real": float(weight), "n_e_m3": n_e_m3,
            "wpe_rad_s": float(wpe), "wb_rad_s": float(wb),
            "v0_m_s": float(v0), "dt_s": float(dt),
            "gamma_theory_rad_s": float(twostream_growth_rate(k, v0, wb)),
            "amp_frac": float(amp_frac), "vth_frac": float(vth_frac),
            "seed": int(seed), "grid": [nx, nx, nx],
            "periods": float(periods),
            "note": ("REAL macro-particles: cold symmetric two-stream "
                     "at the fastest-growing mode (kappa = sqrt(3)/2, "
                     "gamma = wb/2 closed-form) — self-generated "
                     "bunching and phase-space vortices; reduced-"
                     "density cm box (plan §3), not the laser chamber"),
        },
    }


# ------------------------------------------------------------------- npz
def save_npz(path, F):
    arrays = {"t_s": F["t_s"], "pos": F["pos"], "vz": F["vz"],
              "ke_J": F["ke_J"], "mode_re": F["mode_re"],
              "box_m": np.float64(F["box_m"]),
              "meta_json": np.array(json.dumps(F["meta"]))}
    for k in ("mode_im", "beam", "bunch"):
        if k in F:
            arrays[k] = F[k]
    np.savez_compressed(Path(path), **arrays)


def load_npz(path):
    z = np.load(Path(path), allow_pickle=False)
    out = {"t_s": z["t_s"], "pos": z["pos"], "vz": z["vz"],
           "ke_J": z["ke_J"], "mode_re": z["mode_re"],
           "box_m": float(z["box_m"]),
           "meta": json.loads(str(z["meta_json"]))}
    for k in ("mode_im", "beam", "bunch"):
        if k in z.files:
            out[k] = z[k]
    return out


def mode_abs(F):
    if "mode_im" in F:
        return np.hypot(F["mode_re"], F["mode_im"])
    return np.abs(F["mode_re"])


# --------------------------------------------------------------- rendering
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return matplotlib, plt


BEAM_COLORS = ("#6ec6ff", "#ff8a5e")


def _box_edges(ax, L):
    for s, e in [((0, 0, 0), (L, 0, 0)), ((0, 0, 0), (0, L, 0)),
                 ((0, 0, 0), (0, 0, L)), ((L, L, L), (0, L, L)),
                 ((L, L, L), (L, 0, L)), ((L, L, L), (L, L, 0)),
                 ((L, 0, 0), (L, L, 0)), ((L, 0, 0), (L, 0, L)),
                 ((0, L, 0), (L, L, 0)), ((0, L, 0), (0, L, L)),
                 ((0, 0, L), (L, 0, L)), ((0, 0, L), (0, L, L))]:
        ax.plot(*zip(s, e), color="#3a4250", lw=0.8)


def render_twostream_frame(F, j, fig=None):
    """Two-stream frame: 3-D beams + the (z, v_z) phase space where the
    bunching and trapping vortices live + the |mode| growth trace with
    the closed-form gamma as a slope guide."""
    mpl, plt = _mpl()
    if fig is None:
        fig = plt.figure(figsize=(12.5, 8.0), dpi=110)
    fig.clf()
    fig.patch.set_facecolor("#101318")
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.30),
                          width_ratios=(1.0, 1.05), hspace=0.18,
                          wspace=0.14, left=0.05, right=0.97, top=0.90,
                          bottom=0.08)
    m = F["meta"]
    L = F["box_m"]
    p, v, b = F["pos"][j], F["vz"][j], F["beam"]
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax.set_facecolor("#101318")
    for bi, color in enumerate(BEAM_COLORS):
        s = b == bi
        ax.scatter(p[s, 0], p[s, 1], p[s, 2], s=2.6, c=color,
                   alpha=0.75, linewidths=0, depthshade=False)
    _box_edges(ax, L)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.text2D(0.01, 0.97, f"{m['n_macro']} macro-particles x "
              f"{m['weight_real']:.2e} real e⁻ each",
              transform=ax.transAxes, fontsize=9, color="#9fd0ff")
    ax.text2D(0.01, 0.015, f"box {L * 100:.1f} cm — beams ±v0 = "
              f"{m['v0_m_s']:.2e} m/s along z", transform=ax.transAxes,
              fontsize=8, color="#8a919c")
    axp = fig.add_subplot(gs[0, 1])
    axp.set_facecolor("#171b22")
    vlim = 3.2 * m["v0_m_s"]
    for bi, color in enumerate(BEAM_COLORS):
        s = b == bi
        axp.scatter(p[s, 2] * 100.0, v[s], s=2.0, c=color, alpha=0.6,
                    linewidths=0)
    axp.set_xlim(0, L * 100.0)
    axp.set_ylim(-vlim, vlim)
    axp.set_xlabel("z cm", fontsize=8)
    axp.set_ylabel("v_z m/s", fontsize=8)
    axp.set_title("phase space (z, v_z) — bunching → trapping vortices",
                  fontsize=9, color="#aab0b8")
    axp.tick_params(colors="#aab0b8", labelsize=7)
    axp.xaxis.label.set_color("#aab0b8")
    axp.yaxis.label.set_color("#aab0b8")
    for sp in axp.spines.values():
        sp.set_color("#3a4250")
    axt = fig.add_subplot(gs[1, :])
    axt.set_facecolor("#171b22")
    tt = F["t_s"] * m["wpe_rad_s"] / (2.0 * np.pi)
    amp = F["bunch"][:, 0]
    floor = amp[amp > 0].min() if np.any(amp > 0) else 1e-30
    axt.semilogy(tt, np.clip(amp, floor, None), color="#6ec6ff", lw=1.3,
                 label="beam-0 bunching |⟨e^{-ikz}⟩|")
    g = m["gamma_theory_rad_s"]
    guide = amp[0] * np.exp(g * F["t_s"])
    axt.semilogy(tt, np.clip(guide, floor, amp.max() * 3.0), "--",
                 color="#c9a34a", lw=1.0,
                 label="closed-form γ = ω_b/2 slope")
    axt.axvline(tt[j], color="#888", lw=1.0, ls="--")
    axt.set_xlim(tt[0], tt[-1])
    axt.set_xlabel("t in plasma periods", fontsize=8)
    axt.set_ylabel("bunching", fontsize=8)
    axt.legend(fontsize=7, loc="lower right", framealpha=0.2,
               labelcolor="#aab0b8")
    axt.tick_params(colors="#aab0b8", labelsize=7)
    axt.xaxis.label.set_color("#aab0b8")
    axt.yaxis.label.set_color("#aab0b8")
    for sp in axt.spines.values():
        sp.set_color("#3a4250")
    tau = F["t_s"][j] * m["wpe_rad_s"] / (2.0 * np.pi)
    fig.suptitle("PIC_MCC_3D — two-stream instability, real "
                 f"macro-particles — t = {tau:.2f} plasma periods",
                 fontsize=12, color="#d8dbe0")
    fig.text(0.5, 0.005,
             "self-generated bunching: every dot IS a macro-particle; "
             "the growth follows the closed-form cold-beam rate (nothing "
             "tuned); reduced-density cm box — plan §3 regime",
             ha="center", fontsize=7, color="#767c85")
    return fig


def render_frame(F, j, fig=None):
    """One frame: the macro-particles (colored by v_z) + the mode trace."""
    if "beam" in F:
        return render_twostream_frame(F, j, fig=fig)
    mpl, plt = _mpl()
    if fig is None:
        fig = plt.figure(figsize=(9.5, 8.0), dpi=110)
    fig.clf()
    fig.patch.set_facecolor("#101318")
    gs = fig.add_gridspec(2, 1, height_ratios=(1.0, 0.24), hspace=0.15,
                          left=0.03, right=0.97, top=0.92, bottom=0.08)
    ax = fig.add_subplot(gs[0], projection="3d")
    ax.set_facecolor("#101318")
    L = F["box_m"]
    p = F["pos"][j]
    v = F["vz"][j]
    vmax = float(np.abs(F["vz"]).max())
    vmax = vmax if vmax > 0 else 1.0
    sc = ax.scatter(p[:, 0], p[:, 1], p[:, 2], c=v, cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax, s=3.2, alpha=0.85,
                    linewidths=0, depthshade=False)
    for s, e in [((0, 0, 0), (L, 0, 0)), ((0, 0, 0), (0, L, 0)),
                 ((0, 0, 0), (0, 0, L)), ((L, L, L), (0, L, L)),
                 ((L, L, L), (L, 0, L)), ((L, L, L), (L, L, 0)),
                 ((L, 0, 0), (L, L, 0)), ((L, 0, 0), (L, 0, L)),
                 ((0, L, 0), (L, L, 0)), ((0, L, 0), (0, L, L)),
                 ((0, 0, L), (L, 0, L)), ((0, 0, L), (0, L, L))]:
        ax.plot(*zip(s, e), color="#3a4250", lw=0.8)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    m = F["meta"]
    tau = F["t_s"][j] * m["wpe_rad_s"] / (2.0 * np.pi)
    ax.text2D(0.01, 0.965,
              f"{m['n_macro']} macro-particles x "
              f"{m['weight_real']:.2e} real electrons each",
              transform=ax.transAxes, fontsize=9, color="#9fd0ff")
    ax.text2D(0.01, 0.02, f"box {L * 100:.1f} cm   n_e {m['n_e_m3']:.0e} "
              f"m^-3   color = v_z", transform=ax.transAxes, fontsize=8,
              color="#8a919c")
    cb = fig.colorbar(sc, ax=ax, shrink=0.45, pad=0.02)
    cb.set_label("v_z m/s", fontsize=8, color="#aab0b8")
    cb.ax.tick_params(labelsize=6, colors="#aab0b8")
    axt = fig.add_subplot(gs[1])
    axt.set_facecolor("#171b22")
    tt = F["t_s"] * m["wpe_rad_s"] / (2.0 * np.pi)
    axt.plot(tt, F["mode_re"], color="#6ec6ff", lw=1.3)
    axt.axvline(tau, color="#888", lw=1.0, ls="--")
    axt.set_xlabel("t in plasma periods", fontsize=8)
    axt.set_ylabel("mode amplitude", fontsize=8)
    axt.tick_params(colors="#aab0b8", labelsize=7)
    axt.xaxis.label.set_color("#aab0b8")
    axt.yaxis.label.set_color("#aab0b8")
    for sp in axt.spines.values():
        sp.set_color("#3a4250")
    fig.suptitle(f"PIC_MCC_3D — real macro-particles — cold Langmuir "
                 f"slosh — t = {tau:.2f} plasma periods",
                 fontsize=12, color="#d8dbe0")
    fig.text(0.5, 0.005,
             "every dot IS a PIC macro-particle (leapfrog positions, "
             "v_z color); reduced-density cm box — the honest 3-D PIC "
             "regime, not the laser chamber (plan §3)",
             ha="center", fontsize=7, color="#767c85")
    return fig


def render_animation(F, out_dir, fps=14, gif=True, mp4=False, log=print):
    _, plt = _mpl()
    out_dir = Path(out_dir)
    fdir = out_dir / "frames"
    fdir.mkdir(parents=True, exist_ok=True)
    n_f = len(F["t_s"])
    fig = None
    paths = []
    for k in range(n_f):
        fig = render_frame(F, k, fig=fig)
        p = fdir / f"frame_{k:04d}.png"
        fig.savefig(p, facecolor=fig.get_facecolor())
        paths.append(p)
        if k % 15 == 0 or k == n_f - 1:
            log(f"  frame {k + 1}/{n_f}")
    plt.close(fig)
    made = {"frames": fdir}
    if gif:
        from PIL import Image
        imgs = [Image.open(q).convert("P", palette=Image.ADAPTIVE)
                for q in paths]
        gif_path = out_dir / "movie.gif"
        imgs[0].save(gif_path, save_all=True, append_images=imgs[1:],
                     duration=int(1000 / fps), loop=0, optimize=True)
        made["gif"] = gif_path
    if mp4:
        ff = shutil.which("ffmpeg")
        if ff is None:
            log("  ffmpeg not found — skipping mp4")
        else:
            mp4_path = out_dir / "movie.mp4"
            r = subprocess.run(
                [ff, "-y", "-loglevel", "error", "-framerate", str(fps),
                 "-i", str(fdir / "frame_%04d.png"),
                 "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                 "-pix_fmt", "yuv420p", str(mp4_path)],
                capture_output=True, text=True)
            if r.returncode == 0:
                made["mp4"] = mp4_path
            else:
                log(f"  ffmpeg failed ({r.stderr.strip()[:200]})")
    return made


# ------------------------------------------------------------------- html
def quantize01(vals):
    """[0,1)-normalized values -> uint8 (invert with /255)."""
    return np.clip(np.rint(np.asarray(vals, float) * 255.0), 0,
                   255).astype(np.uint8)


def write_html(F, path, max_frames=80):
    idx = np.unique(np.linspace(0, len(F["t_s"]) - 1,
                                min(max_frames, len(F["t_s"]))).round()
                    .astype(int))
    L = F["box_m"]
    vmax = float(np.abs(F["vz"]).max())
    vmax = vmax if vmax > 0 else 1.0
    m = F["meta"]
    amp = mode_abs(F)
    payload = {
        "meta": m, "box_cm": L * 100.0, "vmax": vmax,
        "tau": [float(F["t_s"][j] * m["wpe_rad_s"] / (2 * np.pi))
                for j in idx],
        "mode_abs": [float(amp[j]) for j in idx],
        "frames": [{
            "p": quantize01(F["pos"][j] / L).ravel().tolist(),
            "v": quantize01((F["vz"][j] / vmax + 1.0) / 2.0).tolist(),
        } for j in idx],
    }
    if "beam" in F:
        payload["beam"] = [int(x) for x in F["beam"]]
    if "bunch" in F:
        payload["bunch"] = [float(F["bunch"][j, 0]) for j in idx]
    title = ("PIC_MCC_3D real macro-particles — two-stream instability"
             if "beam" in F else
             "PIC_MCC_3D real macro-particles — Langmuir slosh")
    html = (_HTML_TEMPLATE
            .replace("__TITLE__", title)
            .replace("__DATA__", json.dumps(payload,
                                            separators=(",", ":"))))
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
  .wrap{max-width:1000px;margin:0 auto;padding:14px}
  h1{font-size:17px;margin:4px 0 2px}
  .row{display:flex;gap:14px;flex-wrap:wrap;align-items:flex-start}
  canvas{background:#0d0f12;border:1px solid #2a2e35;border-radius:6px}
  .ctl{display:flex;flex-direction:column;gap:9px;min-width:220px;flex:1}
  .ctl label{display:flex;align-items:center;gap:8px;font-size:12px;
             color:#aab0b8}
  .ctl input[type=range]{flex:1}
  button{background:#1e2228;color:#d8dbe0;border:1px solid #343a43;
         border-radius:4px;padding:4px 8px;font-size:12px;cursor:pointer}
  .readout{font-family:ui-monospace,Menlo,monospace;font-size:11px;
           color:#9aa0a8;white-space:pre-wrap}
  .cap{color:#767c85;font-size:11px;margin-top:10px}
</style>
</head>
<body>
<div class="wrap">
  <h1>__TITLE__</h1>
  <p class="readout" id="head"></p>
  <div class="row">
    <canvas id="view" width="680" height="540"></canvas>
    <div class="ctl">
      <label>t <input type="range" id="time" min="0" value="0" step="1">
        <span id="tlabel" class="readout"></span></label>
      <label>yaw <input type="range" id="yaw" min="-180" max="180"
        value="-35" step="1"></label>
      <button id="play">play</button>
      <div class="readout" id="strip"></div>
    </div>
  </div>
  <p class="cap">Every dot IS a PIC macro-particle (leapfrog positions;
  color = v_z, blue−/red+). Cold Langmuir slosh in a periodic
  neutralising jelly — the package's frequency-gate configuration made
  visible. Reduced-density cm box: the honest 3-D PIC regime, not the
  laser chamber (plan §3). Positions are uint8-quantized here — exact
  values live in particles_data.npz.</p>
</div>
<script id="vizdata" type="application/json">__DATA__</script>
<script>
"use strict";
const D = JSON.parse(document.getElementById("vizdata").textContent);
const view = document.getElementById("view"), ctx = view.getContext("2d");
const timeSl = document.getElementById("time");
const yawSl = document.getElementById("yaw");
document.getElementById("head").textContent =
  `${D.meta.n_macro} macro-particles x ` +
  `${D.meta.weight_real.toExponential(2)} real electrons each — ` +
  `box ${D.box_cm.toFixed(1)} cm, n_e ${D.meta.n_e_m3.toExponential(0)} m^-3`;
timeSl.max = String(D.tau.length - 1);
const ELEV = 28 * Math.PI / 180;
function project(p, cy, sy, ce, se) {
  const x = p[0]*cy - p[1]*sy, y = p[0]*sy + p[1]*cy, z = p[2];
  return [x, y*se + z*ce, y*ce - z*se];   // [sx, up, farness]
}
function vcolor(f) {           // 0..1 -> blue..white..red
  const r = f < 0.5 ? 90 + 300 * f : 240;
  const b = f > 0.5 ? 240 - 300 * (f - 0.5) * 2 + 60 : 240;
  const g = 120 + 160 * (0.5 - Math.abs(f - 0.5)) * 2 - 40;
  return `rgba(${r|0},${g|0},${Math.min(b,255)|0},0.85)`;
}
function drawFrame() {
  const ti = +timeSl.value;
  const fr = D.frames[ti];
  const yaw = +yawSl.value * Math.PI / 180;
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  const ce = Math.cos(ELEV), se = Math.sin(ELEV);
  const S = Math.min(view.width, view.height) / (2.35 * Math.sqrt(3) / 2);
  const cxp = view.width / 2, cyp = view.height / 2 + 10;
  ctx.clearRect(0, 0, view.width, view.height);
  const C = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],
             [1,1,1]].map(p => {
    const pr = project([p[0]-0.5, p[1]-0.5, p[2]-0.5], cy, sy, ce, se);
    return [cxp + S*pr[0], cyp - S*pr[1]];
  });
  ctx.strokeStyle = "#3a4250"; ctx.lineWidth = 1;
  for (const [a,b] of [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],
                       [3,6],[4,7],[5,7],[6,7]]) {
    ctx.beginPath(); ctx.moveTo(C[a][0], C[a][1]);
    ctx.lineTo(C[b][0], C[b][1]); ctx.stroke();
  }
  const n = fr.v.length;
  const BEAM = D.beam || null;
  const BCOL = ["rgba(110,198,255,0.8)", "rgba(255,138,94,0.8)"];
  // painter: precompute farness, draw far first
  const ord = new Array(n);
  for (let i = 0; i < n; i++) {
    const px = fr.p[3*i] / 255 - 0.5, py = fr.p[3*i+1] / 255 - 0.5,
          pz = fr.p[3*i+2] / 255 - 0.5;
    ord[i] = [project([px, py, pz], cy, sy, ce, se)[2], i];
  }
  ord.sort((a, b) => b[0] - a[0]);
  for (const [, i] of ord) {
    const px = fr.p[3*i] / 255 - 0.5, py = fr.p[3*i+1] / 255 - 0.5,
          pz = fr.p[3*i+2] / 255 - 0.5;
    const pr = project([px, py, pz], cy, sy, ce, se);
    ctx.fillStyle = BEAM ? BCOL[BEAM[i]] : vcolor(fr.v[i] / 255);
    ctx.fillRect(cxp + S*pr[0] - 1.5, cyp - S*pr[1] - 1.5, 3, 3);
  }
  document.getElementById("tlabel").textContent =
    D.tau[ti].toFixed(2) + " T_pe";
  document.getElementById("strip").textContent =
    (D.bunch ? `beam-0 bunching ${D.bunch[ti].toExponential(2)}`
             : `|density mode| ${D.mode_abs[ti].toExponential(2)}`) +
    (BEAM ? "   beams: blue +v0 / orange −v0" : "");
}
let timer = null;
document.getElementById("play").addEventListener("click", (e) => {
  if (timer) { clearInterval(timer); timer = null;
               e.target.textContent = "play"; return; }
  e.target.textContent = "pause";
  timer = setInterval(() => {
    timeSl.value = String((+timeSl.value + 1) % D.tau.length);
    drawFrame();
  }, 70);
});
for (const el of [timeSl, yawSl]) el.addEventListener("input", drawFrame);
drawFrame();
</script>
</body>
</html>
"""


# ------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--demo", choices=("langmuir", "twostream"),
                     default="langmuir")
    src.add_argument("--npz", metavar="PATH",
                     help="re-render from saved particles_data.npz")
    ap.add_argument("--out", metavar="DIR", default=None)
    ap.add_argument("--nx", type=int, default=16)
    ap.add_argument("--per-axis", type=int, default=12,
                    help="macro-particles per axis (cube per beam; "
                         "default 12 -> 1728)")
    ap.add_argument("--periods", type=float, default=None,
                    help="plasma periods (default: 3 langmuir, "
                         "10 twostream)")
    ap.add_argument("--amp", type=float, default=None,
                    help="seed displacement fraction of the box "
                         "(default: 0.05 langmuir, 0.01 twostream)")
    ap.add_argument("--every", type=int, default=None,
                    help="snapshot every N steps (default: 4 langmuir, "
                         "10 twostream)")
    ap.add_argument("--fps", type=int, default=14)
    ap.add_argument("--mp4", action="store_true")
    ap.add_argument("--no-gif", action="store_true")
    ap.add_argument("--no-html", action="store_true")
    ap.add_argument("--no-movie", action="store_true")
    args = ap.parse_args(argv)

    if args.npz:
        F = load_npz(args.npz)
        out = Path(args.out) if args.out else Path(args.npz).parent
        print(f"[picviz] loaded {args.npz}")
    else:
        two = args.demo == "twostream"
        periods = args.periods if args.periods is not None else \
            (10.0 if two else 3.0)
        amp = args.amp if args.amp is not None else \
            (0.01 if two else 0.05)
        n_total = args.per_axis ** 3 * (2 if two else 1)
        print(f"[picviz] running {args.demo} demo (nx={args.nx}, "
              f"{n_total} macro-particles, "
              f"{periods:g} plasma periods) ...")
        every = args.every if args.every is not None else \
            (10 if two else 4)
        runner = run_twostream if two else run_langmuir
        F = runner(nx=args.nx, per_axis=args.per_axis,
                   periods=periods, amp_frac=amp, every=every)
        out = Path(args.out) if args.out else \
            ROOT / "picmcc_out" / f"viz_{args.demo}"
        out.mkdir(parents=True, exist_ok=True)
        save_npz(out / "particles_data.npz", F)
        print(f"[picviz] {len(F['t_s'])} snapshots of "
              f"{F['meta']['n_macro']} macro-particles "
              f"-> {out / 'particles_data.npz'}")
    out.mkdir(parents=True, exist_ok=True)
    if not args.no_movie:
        made = render_animation(F, out, fps=args.fps,
                                gif=not args.no_gif, mp4=args.mp4)
        for k in ("gif", "mp4"):
            if k in made:
                print(f"[picviz] {k} -> {made[k]}")
    if not args.no_html:
        v = write_html(F, out / "viewer.html")
        print(f"[picviz] interactive viewer -> {v}")
    print("[picviz] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
