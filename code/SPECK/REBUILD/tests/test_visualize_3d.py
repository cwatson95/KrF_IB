"""gates for dense per-cell histories + tools/visualize_3d.py.

The viewer contract: every rendered number IS the run's number (engine
accessors, gain_m1's expression), volume averages of the per-cell fields
reproduce the run's own averaged histories to machine precision, the npz
round-trips exactly, and the HTML viewer's quantized values invert to
within one 8-bit bin. Renders are gated headless (Agg).
"""

import ast
import json
from pathlib import Path

import numpy as np
import pytest

from REBUILD.tools import visualize_3d as viz

MODULES = Path(__file__).resolve().parents[1] / "MODULES"


@pytest.fixture(scope="module")
def tiny3d():
    """A 2x1x2-cell KrF beam device, 40 ns — the smallest honest 3-D run."""
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": 40.0, "nx": 2, "ny": 1, "nz": 2,
                "Lx_m": 0.04, "Ly_m": 0.0141, "Lz_m": 0.10,
                "beam_axis": 0}, source="problem")
    dev = Device(p)
    res = dev.run(n_save=5)              # dense auto-on: n_cells == 4
    return dev, res


@pytest.fixture(scope="module")
def tiny_fields(tiny3d):
    dev, res = tiny3d
    F = viz.extract_fields(dev, res)
    return F, viz.field_specs(F)


# ------------------------------------------------------- engine: dense runs
def test_dense_auto_on_multicell(tiny3d):
    dev, res = tiny3d
    assert "y_t" in res
    assert res["y_t"].shape == (dev.lay.size, 5)


def test_dense_auto_off_0d(krf_result):
    assert "y_t" not in krf_result


def test_dense_explicit_on_0d(krf_device):
    res = krf_device.run(n_save=3, dense_cells=True)
    assert res["y_t"].shape == (krf_device.lay.size, 3)


def test_volume_averages_reproduce_collect(tiny3d):
    """The averaged histories and y_t describe the SAME run: recomputing
    the volume-weighted averages from the dense state matches _collect."""
    dev, res = tiny3d
    lay, w = dev.lay, dev.cell_vols
    v_g = float(w.sum())
    for sp in (dev.optics["upper"], "F2"):
        i = lay.index[sp]
        for j in range(res["y_t"].shape[1]):
            dens = np.maximum(lay.density_block(res["y_t"][:, j]), 0.0)
            avg = float(dens[i] @ w) / v_g
            assert np.isclose(avg, res["dens_avg"][sp][j], rtol=1e-12,
                              atol=0.0)
    for j in range(res["y_t"].shape[1]):
        dens = np.maximum(lay.density_block(res["y_t"][:, j]), 0.0)
        ne = float(lay.electron_density(dens) @ w) / v_g
        assert np.isclose(ne, res["ne_m3"][j], rtol=1e-12, atol=0.0)


# ------------------------------------------------------------- field extract
def test_extract_fields_shapes(tiny3d, tiny_fields):
    dev, res = tiny3d
    F, _ = tiny_fields
    n_t = res["y_t"].shape[1]
    for name, vals in F["fields"].items():
        assert vals.shape == (n_t, 2, 1, 2), name
    assert (F["fields"]["n_e"] > 0).all()
    assert "dep_shape" in F                     # beam drive
    assert F["dep_shape"].shape == (2, 1, 2)
    assert np.isclose(F["dep_shape"].max(), 1.0)
    assert list(F["shape"]) == [2, 1, 2]


def test_extract_requires_dense(tiny3d, krf_result):
    dev, _ = tiny3d
    with pytest.raises(ValueError, match="dense_cells"):
        viz.extract_fields(dev, {k: v for k, v in krf_result.items()
                                 if k != "y_t"})


def test_local_gain_averages_to_collect_gain(tiny3d, tiny_fields):
    """gain_m1 is linear in densities, so the volume average of the
    per-cell g0 field must equal _collect's averaged-density gain."""
    dev, res = tiny3d
    F, _ = tiny_fields
    w = dev.cell_vols
    v_g = float(w.sum())
    for j in range(F["fields"]["g0"].shape[0]):
        avg = float(F["fields"]["g0"][j].ravel() @ w) / v_g
        assert np.isclose(avg, res["g0_m1"][j], rtol=1e-10, atol=1e-12)


# --------------------------------------------------------------- npz + bits
def test_npz_roundtrip(tiny_fields, tmp_path):
    F, _ = tiny_fields
    p = tmp_path / "viz_data.npz"
    viz.save_npz(p, F)
    G = viz.load_npz(p)
    assert G["meta"] == F["meta"]
    assert np.array_equal(G["t_ns"], F["t_ns"])
    assert np.array_equal(G["cell_vols"], F["cell_vols"])
    assert np.array_equal(G["dep_shape"], F["dep_shape"])
    assert set(G["fields"]) == set(F["fields"])
    for name in F["fields"]:
        assert np.array_equal(G["fields"][name], F["fields"][name]), name


def test_quantize_inverts_within_one_bin():
    v = np.linspace(2.0, 350.0, 977)
    q = viz.quantize(v, 2.0, 350.0, "lin")
    back = viz.dequantize(q, 2.0, 350.0, "lin")
    assert np.abs(back - v).max() <= (350.0 - 2.0) / 255.0 * 0.51
    vl = np.logspace(np.log10(3e8), np.log10(3e13), 977)
    ql = viz.quantize(vl, 3e8, 3e13, "log")
    backl = viz.dequantize(ql, 3e8, 3e13, "log")
    bin_f = np.log(3e13 / 3e8) / 255.0
    assert np.abs(np.log(backl / vl)).max() <= bin_f * 0.51


# ---------------------------------------------------------------- rendering
def test_frame_renders_headless(tiny_fields, tmp_path):
    F, specs = tiny_fields
    fig = viz.render_frame(F, specs, j=F["t_ns"].size - 1)
    out = tmp_path / "frame.png"
    fig.savefig(out)
    _, plt = viz._mpl()
    plt.close(fig)
    assert out.stat().st_size > 20_000


def test_contact_sheet_and_gif(tiny_fields, tmp_path):
    F, specs = tiny_fields
    sheet = viz.render_contact_sheet(F, specs, tmp_path / "sheet.png", n=4)
    assert Path(sheet).stat().st_size > 20_000
    made = viz.render_animation(F, specs, tmp_path, frames=3,
                                gif=True, mp4=False,
                                log=lambda *a, **k: None)
    assert sorted(made["frames"].glob("frame_*.png"))
    assert made["gif"].stat().st_size > 5_000


# ---------------------------------------------------------- particle view
def test_sampler_deterministic_unit_interval():
    a = viz.sample_positions(7, 5, 0)
    b = viz.sample_positions(7, 5, 0)
    assert np.array_equal(a, b)
    assert a.shape == (7, 5, 3)
    assert (a >= 0.0).all() and (a < 1.0).all()
    assert not np.array_equal(a, viz.sample_positions(7, 5, 1))


def test_particle_counts_linear_in_particle_number():
    vols = np.array([1.0, 2.0])                            # m^3
    vals = np.array([[0.0, 0.0], [2.0, 3.0], [1.0, 1.5]])  # m^-3
    counts, w = viz.particle_counts(vals, vols, 100)
    assert counts[0].sum() == 0            # below one macro-weight: no dots
    assert counts[1].sum() == 100          # peak frame hits the target
    assert counts[2].sum() == 50           # half the real particle number
    assert np.isclose(w * counts[1].sum(), float((vals[1] * vols).sum()))


def test_particle_frame_and_gif(tiny_fields, tmp_path):
    F, specs = tiny_fields
    pc = viz.prepare_particles(F, pts=500)
    assert {p["field"] for p in pc} == {"n_e", F["meta"]["upper"]}
    fig = viz.render_particle_frame(F, pc, j=F["t_ns"].size - 1)
    out = tmp_path / "pframe.png"
    fig.savefig(out, facecolor=fig.get_facecolor())
    _, plt = viz._mpl()
    plt.close(fig)
    assert out.stat().st_size > 20_000
    made = viz.render_animation(F, specs, tmp_path, frames=2, gif=True,
                                mp4=False, style="particles", pts=300,
                                log=lambda *a, **k: None)
    assert made["gif"].name == "movie_particles.gif"
    assert made["gif"].stat().st_size > 2_000


# --------------------------------------------------------------- glow view
def test_glow_cloud_math(tiny_fields):
    F, _ = tiny_fields
    em = viz.emission_fields(F)
    assert set(em) == {"fluor", "recomb"}
    fine = viz.refine_field(em["fluor"][-1], factor=3)
    assert fine.shape == (6, 3, 6)
    assert (fine >= 0.0).all()
    a = viz.glow_cloud(em["fluor"][-1], 500, seed_offset=11)
    b = viz.glow_cloud(em["fluor"][-1], 500, seed_offset=11)
    assert np.array_equal(a, b)                     # deterministic
    assert len(a) == 500                            # budget EXACT
    assert (a >= 0.0).all()
    assert (a[:, 0] <= 2.0).all() and (a[:, 2] <= 2.0).all()
    assert len(viz.glow_cloud(np.zeros((2, 1, 2)), 500, 1)) == 0
    # the deep-afterglow regression: a budget far below the fine-voxel
    # count must still land its dots (per-voxel rint floored it to zero)
    assert len(viz.glow_cloud(em["fluor"][-1], 7, seed_offset=11)) == 7


def test_glow_frame_and_headline_movie_names(tiny_fields, tmp_path):
    F, specs = tiny_fields
    prep = viz.prepare_glow(F, pts=800)
    fig = viz.render_glow_frame(F, prep, j=F["t_ns"].size - 1)
    out = tmp_path / "g.png"
    fig.savefig(out, facecolor=fig.get_facecolor())
    _, plt = viz._mpl()
    plt.close(fig)
    viz.bloom_png(out)
    assert out.stat().st_size > 20_000
    made = viz.render_animation(F, specs, tmp_path, frames=2, gif=True,
                                mp4=False, style="glow",
                                log=lambda *a, **k: None)
    assert made["gif"].name == "movie.gif"          # the headline file
    made2 = viz.render_animation(F, specs, tmp_path, frames=2, gif=True,
                                 mp4=False, style="voxels",
                                 log=lambda *a, **k: None)
    assert made2["gif"].name == "movie_fields.gif"


def test_bragg_preset_has_real_spatial_structure():
    """The structure demo is physics, not decoration: at 100 kV into
    6 cm of 3 atm the deposition rises to a Bragg peak inside the cell
    and is ZERO past the stop."""
    dev = viz.build_device("krf_beam_3d_bragg")
    dep = dev._dep_frac.reshape(*dev.shape)[:, 0, 0]
    dep = dep / dep.max()
    assert 3 <= int(np.argmax(dep)) <= 10       # peak inside the domain
    assert dep[-1] < 1e-6                       # fully stopped
    assert dep.max() / dep[0] > 2.0             # visible Bragg rise


# --------------------------------------------------------------------- html
def test_html_viewer_selfcontained_and_invertible(tiny_fields, tmp_path):
    F, specs = tiny_fields
    out = viz.write_html(F, specs, tmp_path / "viewer.html", frames=4)
    text = out.read_text()
    assert "http" not in text.lower()           # no external requests
    payload = json.loads(
        text.split('<script id="vizdata" type="application/json">')[1]
        .split("</script>")[0])
    assert payload["shape"] == [2, 1, 2]
    assert len(payload["t_ns"]) == len(payload["fields"]["n_e"]["frames"])
    spec = payload["fields"]["n_e"]
    j = viz._frame_indices(F["t_ns"].size, 4)[-1]
    true = np.clip(F["fields"]["n_e"][j].ravel(), spec["vmin"],
                   spec["vmax"])
    back = viz.dequantize(np.array(spec["frames"][-1]), spec["vmin"],
                          spec["vmax"], spec["scale"])
    bin_f = np.log(spec["vmax"] / spec["vmin"]) / 255.0
    assert np.abs(np.log(back / true)).max() <= bin_f * 0.51
    assert payload["fields"]["deposition"]["static"] is True
    # particle-view payload: volumes + per-species normalization + controls
    pd = payload["particles"]
    assert len(pd["cellvol"]) == 4
    assert {s["field"] for s in pd["species"]} == \
        {"n_e", F["meta"]["upper"]}
    assert all(s["pnorm"] > 0 for s in pd["species"])
    assert 'name="mode"' in text
    assert 'value="glow" checked' in text           # camera view default


# ----------------------------------------------------------- architecture
def test_modules_never_import_tools():
    """README's tools/ contract, mechanised: the engine never imports a
    tool — couplings flow tool -> MODULES, and only as file outputs."""
    for py in sorted(MODULES.glob("*.py")):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    mods = [node.module]
                elif node.level:                 # from . import x
                    mods = [a.name for a in node.names]
            for mod in mods:
                assert "tools" not in mod.split("."), (py.name, mod)
