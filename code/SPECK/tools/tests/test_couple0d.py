"""gates for SPECK/tools/couple0d.py -- the coupler's reproduce-limit.

THE gate: clamp the coupler to a REBUILD engine run's own recorded
(n_e, T_e, T_g, n_ph) trajectory across an afterglow window (drive checked
to be exactly zero there) and require every major heavy species to land on
the engine's answer.  Same compiled network, same clamps => same chemistry;
if this fails, the coupler is advancing different physics than the engine
and has no business inside a PIC loop.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_SPECK = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "couple0d", _SPECK / "tools" / "couple0d.py")
couple0d = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("couple0d", couple0d)
_spec.loader.exec_module(couple0d)


@pytest.fixture(scope="module")
def engine_run():
    """A short 0-D KrF beam shot with dense saves, afterglow included."""
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "nx": 1, "ny": 1, "nz": 1, "t_end_ns": 140.0},
               source="api")
    # n_save sets the CLAMP resolution, and clamp resolution is the coupler's
    # error knob (measured 2026-08-14: 2 ns clamps -> 16% worst species (F-),
    # 1 ns clamps -> 4% -- interpolation error, not structural drift).  1 ns
    # here; the same number is the macro-step guidance for a real PIC loop.
    dev = Device(p)
    res = dev.run(n_save=141, dense_cells=True)  # 0-D keeps y_t only if asked
    assert res["success"]
    return dev, res


def _clamps_at(dev, res, j):
    lay = dev.lay
    y = res["y_t"][:, j]
    cells = lay.cells_view(y)
    dens = np.maximum(cells[: lay.n_species], 0.0)
    return {
        "ne": float(lay.electron_density(dens)[0]),
        "te": float(dev.te_of(y)[0]),
        "tg": float(cells[lay.i_tg][0]),
        "n_ph": float(y[lay.i_ph]),
    }, dens.copy()


def test_reproduce_limit_afterglow_parity(engine_run):
    dev, res = engine_run
    lay = dev.lay
    t_s = res["t_ns"] * 1e-9
    # the afterglow window: drive strictly zero at every save used
    drive = np.array([dev.spec.power_areal_W_m2(t) for t in t_s])
    idx = np.nonzero(drive == 0.0)[0]
    assert idx.size >= 10, "no afterglow saves -- lengthen t_end_ns"
    j0 = int(idx[1])            # skip one save clear of the beam edge
    j1 = int(len(t_s) - 1)
    chem = couple0d.ClampedChemistry(device=dev)
    clamp, dens = _clamps_at(dev, res, j0)
    for j in range(j0 + 1, j1 + 1):
        clamp1, _ = _clamps_at(dev, res, j)
        dt = float(t_s[j] - t_s[j - 1])
        dens = chem.advance(dens, clamp, clamp1, dt)["dens"]
        clamp = clamp1
    _, dens_engine = _clamps_at(dev, res, j1)
    compared = 0
    worst = ("", 0.0)
    for s, i in lay.index.items():
        eng = float(dens_engine[i, 0])
        got = float(dens[i, 0])
        if eng < 1e14:
            continue
        compared += 1
        rel = abs(got - eng) / eng
        if rel > worst[1]:
            worst = (s, rel)
        assert rel < 0.15, (f"{s}: coupler {got:.4e} vs engine {eng:.4e} "
                            f"({rel:.1%}) -- reproduce-limit broken")
    assert compared >= 8, f"only {compared} species above the floor"
    print(f"parity over {compared} species, worst {worst[0]} {worst[1]:.2%}")


def test_advance_is_deterministic_and_shaped(engine_run):
    dev, _ = engine_run
    chem = couple0d.ClampedChemistry(device=dev)
    d0 = chem.initial_dens()
    assert d0.shape == (chem.lay.n_species, 1)
    clamp = {"ne": 1e20, "te": 1.2, "tg": 400.0, "n_ph": 0.0}
    a = chem.advance(d0, clamp, clamp, 5e-9)["dens"]
    b = chem.advance(d0, clamp, clamp, 5e-9)["dens"]
    assert np.array_equal(a, b)
    assert np.isfinite(a).all() and (a >= 0.0).all()


def test_contract_rejects_missing_clamps(engine_run):
    dev, _ = engine_run
    chem = couple0d.ClampedChemistry(device=dev)
    d0 = chem.initial_dens()
    with pytest.raises(TypeError):
        chem.advance(d0, {"te": 1.0}, {"te": 1.0}, 1e-9)   # ne missing
