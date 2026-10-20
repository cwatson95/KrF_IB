"""The scorecard: candidate constant sets against every dataset.

Kinetic datasets report rms % deviation, or the signed % for single values.
Device datasets report the model value and its ratio to the measurement.
Each set's electron quench is re-derived first (forward.with_rederived_k_e),
so every set is self-consistent.
"""
from __future__ import annotations

import math

import numpy as np

from . import data, fit, forward, params


def candidate_sets(train=None) -> dict:
    """LIGHT's constants, the two published sets it replaced, and the live
    fit. Since 2026-09-28 'LIGHT now' is the adopted joint fit."""
    sets = {"LIGHT now": params.light_now(),
            "LIGHT to 09-27": params.light_to_0927(),
            "Gadomski 1992": params.gadomski_1992()}
    train = data.by_role("train") if train is None else train
    sets["fit (train set)"] = params.consts(fit.fit(train))
    return {k: forward.with_rederived_k_e(v) for k, v in sets.items()}


def kinetic_score(ds, c) -> float:
    pred = forward.predict(ds, c)
    meas = np.array(ds.values)
    if len(meas) == 1:
        return 100.0 * (pred[0] / meas[0] - 1.0)
    return 100.0 * math.sqrt(np.mean((pred / meas - 1.0) ** 2))


def scorecard(sets=None) -> list:
    sets = sets or candidate_sets()
    rows = []
    for ds in data.all_datasets():
        if ds.kind in data.DEVICE_KINDS:
            continue
        rows.append((ds, {k: kinetic_score(ds, c) for k, c in sets.items()}))
    from . import formation
    dev = {k: formation.observables(c) for k, c in sets.items()}
    key = {"device_g0": "g0", "device_formation": "eta",
           "device_isat": "isat"}
    for ds in data.by_role("device"):
        rows.append((ds, {k: dev[k][key[ds.kind]] for k in sets}))
    return rows, dev


def device_posterior(chain, n_draws=40, seed=5) -> tuple:
    """Tahoma g0 over posterior draws: 16th, 50th and 84th percentiles."""
    rng = np.random.default_rng(seed)
    g = []
    for x in chain[rng.choice(len(chain), n_draws, replace=False)]:
        c = forward.with_rederived_k_e(params.consts(x))
        g.append(forward.run_device(c)["g0"])
    return tuple(np.percentile(g, [16, 50, 84]))
