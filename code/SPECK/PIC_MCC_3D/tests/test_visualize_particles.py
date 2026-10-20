"""gates for tools/visualize_particles.py — real macro-particles, visible.

The contract: snapshots ARE the PIC state (shapes/bounds/energies match
the loop), runs are deterministic, the npz round-trips, renders work
headless, and the HTML viewer is self-contained with positions that
dequantize to the true coordinates within one 8-bit bin.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from PIC_MCC_3D.tools import visualize_particles as vp

TOOLS = Path(__file__).resolve().parents[1] / "tools"


@pytest.fixture(scope="module")
def tiny_run():
    """8^3 grid, 216 macro-particles, 1.2 plasma periods (~seconds)."""
    return vp.run_langmuir(nx=8, per_axis=6, periods=1.2, every=5)


def test_snapshots_are_the_pic_state(tiny_run):
    F = tiny_run
    n_f = len(F["t_s"])
    assert F["pos"].shape == (n_f, 216, 3)
    assert F["vz"].shape == (n_f, 216)
    L = F["box_m"]
    assert (F["pos"] >= 0.0).all() and (F["pos"] <= L).all()
    # cold quiet start: snapshot 0 carries only the leapfrog's retarded
    # half-step kick, far below the slosh peak
    assert F["ke_J"][0] < 1e-2 * F["ke_J"].max()
    assert F["ke_J"].max() > 0.0          # ... and it does slosh
    assert F["meta"]["n_macro"] == 216
    assert F["meta"]["weight_real"] > 0.0


def test_mode_amplitude_oscillates(tiny_run):
    """The visible signal is the frequency gate's signal: it must swing
    through zero repeatedly over 1.2 periods."""
    s = np.sign(tiny_run["mode_re"])
    assert (np.diff(s[s != 0]) != 0).sum() >= 2


def test_run_is_deterministic():
    a = vp.run_langmuir(nx=8, per_axis=4, periods=0.3, every=6)
    b = vp.run_langmuir(nx=8, per_axis=4, periods=0.3, every=6)
    assert np.array_equal(a["pos"], b["pos"])
    assert np.array_equal(a["vz"], b["vz"])


def test_npz_roundtrip(tiny_run, tmp_path):
    p = tmp_path / "particles_data.npz"
    vp.save_npz(p, tiny_run)
    G = vp.load_npz(p)
    assert G["meta"] == tiny_run["meta"]
    assert np.array_equal(G["pos"], tiny_run["pos"])
    assert np.array_equal(G["vz"], tiny_run["vz"])
    assert G["box_m"] == tiny_run["box_m"]


def test_frame_and_gif_render_headless(tiny_run, tmp_path):
    fig = vp.render_frame(tiny_run, len(tiny_run["t_s"]) // 2)
    out = tmp_path / "f.png"
    fig.savefig(out, facecolor=fig.get_facecolor())
    _, plt = vp._mpl()
    plt.close(fig)
    assert out.stat().st_size > 20_000
    small = {k: (v[:2] if isinstance(v, np.ndarray) and v.ndim >= 1
                 and len(v) == len(tiny_run["t_s"]) else v)
             for k, v in tiny_run.items()}
    made = vp.render_animation(small, tmp_path, gif=True, mp4=False,
                               log=lambda *a, **k: None)
    assert made["gif"].stat().st_size > 2_000


def test_html_selfcontained_and_invertible(tiny_run, tmp_path):
    out = vp.write_html(tiny_run, tmp_path / "viewer.html", max_frames=5)
    text = out.read_text()
    assert "http" not in text.lower()
    payload = json.loads(
        text.split('<script id="vizdata" type="application/json">')[1]
        .split("</script>")[0])
    assert payload["meta"]["n_macro"] == 216
    assert len(payload["frames"]) == len(payload["tau"]) <= 5
    q = np.array(payload["frames"][0]["p"]).reshape(-1, 3)
    back = q / 255.0 * tiny_run["box_m"]
    err = np.abs(back - tiny_run["pos"][0])
    assert err.max() <= tiny_run["box_m"] / 255.0 * 0.51


# ------------------------------------------------------------ two-stream
@pytest.fixture(scope="module")
def twostream_run():
    """8^3 grid, 2x216 macro-particles, 8 plasma periods (~seconds)."""
    return vp.run_twostream(nx=8, per_axis=6, periods=8.0, every=8)


def test_twostream_growth_matches_closed_form(twostream_run):
    """The bunching is the cold-beam instability, not an artefact: the
    seeded mode grows exponentially at the closed-form rate (kappa =
    sqrt(3)/2 by construction -> gamma = wb/2), then saturates by
    particle trapping. Nothing here is tunable."""
    F = twostream_run
    amp = F["bunch"][:, 0]          # per-beam bunching: the clean,
    g_th = F["meta"]["gamma_theory_rad_s"]   # cancellation-free signal
    assert np.isclose(g_th, F["meta"]["wb_rad_s"] / 2.0, rtol=1e-9)
    # fit the log-linear window: clear of the seed and of trapping
    lo, hi = 3.0 * amp[0], 0.3 * amp.max()
    win = (amp > lo) & (amp < hi) & \
        (np.arange(len(amp)) < int(np.argmax(amp)))
    assert win.sum() >= 4, "no clean exponential window resolved"
    t, y = F["t_s"][win], np.log(amp[win])
    slope = np.polyfit(t, y, 1)[0]
    assert 0.8 * g_th <= slope <= 1.2 * g_th, (slope, g_th)
    # saturation by trapping: the peak is reached well before the end
    assert int(np.argmax(amp)) < 0.8 * len(amp)


def test_twostream_beams_bunch_in_space(twostream_run):
    """Visible grouping, quantified: the z-histogram of one beam starts
    near-uniform (quiet start) and develops strong modulation."""
    F = twostream_run
    L = F["box_m"]
    b0 = F["beam"] == 0
    # bins = lattice planes (6 for per_axis=6), else the histogram
    # aliases the quiet-start lattice into fake modulation
    edges = np.linspace(0.0, L, 7)

    def modulation(j):
        h, _ = np.histogram(F["pos"][j][b0, 2], bins=edges)
        return h.std() / h.mean()
    early = modulation(0)
    late = max(modulation(j) for j in range(len(F["t_s"])))
    assert early < 0.05                     # quiet start: near-uniform
    assert late > 0.5                       # strong real grouping forms


def test_twostream_schema_roundtrip_and_html(twostream_run, tmp_path):
    F = twostream_run
    p = tmp_path / "ts.npz"
    vp.save_npz(p, F)
    G = vp.load_npz(p)
    assert np.array_equal(G["beam"], F["beam"])
    assert np.array_equal(G["mode_im"], F["mode_im"])
    out = vp.write_html(F, tmp_path / "v.html", max_frames=4)
    text = out.read_text()
    assert "http" not in text.lower()
    payload = json.loads(
        text.split('<script id="vizdata" type="application/json">')[1]
        .split("</script>")[0])
    assert payload["beam"] == [int(x) for x in F["beam"]]
    assert len(payload["mode_abs"]) == len(payload["tau"])
    fig = vp.render_frame(F, len(F["t_s"]) - 1)
    q = tmp_path / "ts.png"
    fig.savefig(q, facecolor=fig.get_facecolor())
    _, plt = vp._mpl()
    plt.close(fig)
    assert q.stat().st_size > 20_000


def test_modules_never_import_tools():
    import ast
    mods = Path(__file__).resolve().parents[1] / "MODULES"
    for py in sorted(mods.glob("*.py")):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            mods_named = []
            if isinstance(node, ast.Import):
                mods_named = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    mods_named = [node.module]
                elif node.level:
                    mods_named = [a.name for a in node.names]
            for mod in mods_named:
                assert "tools" not in mod.split("."), (py.name, mod)
