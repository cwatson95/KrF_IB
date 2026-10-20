"""Save a run, resume it — with honestly-stated semantics.

ADAPTED from ``PIC_MCC_3D/MODULES/checkpoint.py`` (itself from
projects/MHD flare/io.py, FLASH ch. 9): full state to ``.npz`` with the
parameter set AND its provenance, so a checkpoint answers "what exactly was
this run?".

SEMANTICS DIFFERENCE, stated up front: PIC_MCC_3D's restart gate is
BITWISE-exact because its loop is fixed-step explicit with an owned RNG.
This engine integrates with ADAPTIVE BDF, whose internal step history is not
(and should not be) part of the physics state — so a resumed run re-enters
the same trajectory to integrator tolerance, not bit-for-bit. The gate for
this module is therefore: resume(t1) -> t2 matches an uninterrupted run at
t2 within the run's own rtol-class error. Claiming bitwise here would be
claiming a property the method does not have.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FORMAT = "rebuild-checkpoint-1"


def save(path, device, t_s: float, y: np.ndarray) -> None:
    p = device.p
    meta = {
        "format": FORMAT,
        "mode": device.mode,
        "mixture": device.mixture,
        "params": {k: v for k, v in p.as_dict().items()},
        "params_provenance": p.provenance(),
        "beam_yields_mode": (device.beam_yields.mode
                             if device.beam_yields else None),
        "tier_runs": getattr(device.tiers, "run_names", None),
    }
    np.savez(Path(path), y=np.asarray(y, dtype=float),
             t_s=np.array([float(t_s)]), meta=json.dumps(meta))


def restore(path):
    """-> (device, t_s, y). The Device is rebuilt from the saved parameters
    (declared-parameter validation applies again on the way in)."""
    from .engine import Device
    from .params import engine_parameters
    z = np.load(Path(path), allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    if meta.get("format") != FORMAT:
        raise ValueError(f"not a {FORMAT} file: {path}")
    P = engine_parameters()
    P.set_many(meta["params"], source="file")
    dev = Device(P, mixture=meta["mixture"])
    y = np.asarray(z["y"], dtype=float)
    if y.size != dev.lay.size:
        raise ValueError(
            f"state size {y.size} != rebuilt layout {dev.lay.size} — the "
            "code's species set changed since this checkpoint was written")
    return dev, float(z["t_s"][0]), y


def resume(path, t_end_ns: float = None, n_save: int = 121):
    """Continue a checkpointed run to ``t_end_ns`` (default: the saved
    configuration's own t_end)."""
    from scipy.integrate import solve_ivp
    dev, t0, y0 = restore(path)
    t_end = (dev.p["t_end_ns"] if t_end_ns is None else float(t_end_ns)) * 1e-9
    if t_end <= t0:
        raise ValueError(f"t_end {t_end * 1e9:.1f} ns is not after the "
                         f"checkpoint time {t0 * 1e9:.1f} ns")
    t_eval = np.linspace(t0, t_end, int(n_save))
    atol = np.full(dev.lay.size, 1e3)
    cells_atol = dev.lay.cells_view(atol)
    cells_atol[dev.lay.i_ue] = 1e-10
    cells_atol[dev.lay.i_tg] = 1e-6
    atol[dev.lay.i_ph] = 1e-3
    atol[dev.lay.i_ledger0:] = 1e-12
    if dev.lay.n_aux:
        atol[dev.lay.i_aux0:dev.lay.i_aux0 + 5] = [1e-12, 1e-12, 1e-6,
                                                   1e-6, 1e-2]
    kw = {}
    if dev.n_cells > 1:
        kw["jac_sparsity"] = dev.jac_sparsity()
    import time as _time
    w0 = _time.time()
    sol = solve_ivp(dev.rhs, (t0, t_end), y0, method="BDF",
                    rtol=dev.p["rtol"], atol=atol, t_eval=t_eval,
                    max_step=2e-9, **kw)
    # same dense-cells rule as Device.run: a multi-cell resume keeps y_t
    return dev._collect(sol, y0, _time.time() - w0, dev.n_cells > 1)
