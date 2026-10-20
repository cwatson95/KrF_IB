"""checkpoint gates: roundtrip fidelity, resume within tolerance."""

import numpy as np
import pytest

from REBUILD.MODULES.checkpoint import restore, resume, save
from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters


def _small_params(t_end=80.0):
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": t_end, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    return p


def test_roundtrip_preserves_state_and_provenance(tmp_path):
    dev = Device(_small_params())
    y = dev.initial_state()
    y[dev.lay.index["Kr+"]] = 3.14159e17
    f = tmp_path / "ck.npz"
    save(f, dev, 12.5e-9, y)
    dev2, t2, y2 = restore(f)
    assert t2 == 12.5e-9
    assert np.array_equal(y, y2)                     # the ARRAY is bitwise
    assert dev2.p["pressure_atm"] == 3.0
    assert dev2.p.source("pressure_atm") == "file"
    assert dev2.mode == dev.mode and dev2.mixture == dev.mixture


def test_resume_matches_uninterrupted_within_tolerance(tmp_path):
    """Adaptive BDF: resumed-vs-uninterrupted agree to integration
    tolerance, NOT bitwise — the honestly-stated semantics (module doc)."""
    full = Device(_small_params(80.0)).run(n_save=17)
    half_dev = Device(_small_params(40.0))
    r_half = half_dev.run(n_save=9)
    f = tmp_path / "mid.npz"
    save(f, half_dev, r_half["t_ns"][-1] * 1e-9, r_half["y_final"])
    r_res = resume(f, t_end_ns=80.0, n_save=9)
    ne_full = full["ne_m3"][-1]
    ne_res = r_res["ne_m3"][-1]
    assert abs(ne_res / ne_full - 1.0) < 5e-3
    assert abs(r_res["dens_avg"]["KrF*(B)"][-1]
               / full["dens_avg"]["KrF*(B)"][-1] - 1.0) < 2e-2


def test_wrong_format_raises(tmp_path):
    import json
    f = tmp_path / "junk.npz"
    np.savez(f, y=np.zeros(3), t_s=np.array([0.0]),
             meta=json.dumps({"format": "other"}))
    with pytest.raises(ValueError, match="rebuild-checkpoint"):
        restore(f)


def test_species_set_change_detected(tmp_path):
    dev = Device(_small_params())
    f = tmp_path / "ck.npz"
    save(f, dev, 0.0, dev.initial_state()[:-1])      # truncated state
    with pytest.raises(ValueError, match="species set changed|state size"):
        restore(f)
