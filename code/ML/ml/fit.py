"""Fitting the channels to kinetic datasets: least squares, the consistency
check, and an ensemble MCMC (numpy only).

Residuals follow each dataset's error model. Rates and lifetimes use
log(model / measured) / sigma (a fractional error). Ratios use
(model - measured) / sigma (an absolute error). The tau_B prior adds one
term. chi2 per point near 1 means a dataset agrees with its own stated
error.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

from . import forward, params
from .data import DEVICE_KINDS


def _check(datasets):
    for d in datasets:
        if d.kind in DEVICE_KINDS:
            raise ValueError(f"{d.name}: device data score and never train "
                             "(predictive-only policy, ml/data.py)")


def residuals(x, datasets):
    c = params.consts(x)
    out = []
    for d in datasets:
        pred = forward.predict(d, c)
        meas = np.array(d.values)
        if d.sigma_is_fraction:
            out.append(np.log(pred / meas) / d.sigma)
        else:
            out.append((pred - meas) / d.sigma)
    mu, sd = params.TAU_B_PRIOR
    out.append(np.array([(math.exp(x[3]) - mu) / sd]))
    return np.concatenate(out)


def chi2_per_point(x, datasets) -> dict:
    c = params.consts(x)
    res = {}
    for d in datasets:
        pred = forward.predict(d, c)
        meas = np.array(d.values)
        r = (np.log(pred / meas) / d.sigma if d.sigma_is_fraction
             else (pred - meas) / d.sigma)
        res[d.name] = float(np.mean(r ** 2))
    return res


def fit(datasets, starts=(1e-11, 3e-11, 1e-10, 3e-10)):
    """Least squares from several k_BC starts; returns the best vector."""
    _check(datasets)
    best = None
    for kbc in starts:
        x0 = np.log([kbc, 2.0, 60e-9, 6.9e-9, 3e-12, 3e-12, 8e-31, 4e-31,
                     4e-10, 3e-10])
        sol = least_squares(lambda x: residuals(x, datasets), x0,
                            bounds=(params.LO, params.HI), diff_step=1e-3,
                            max_nfev=1500)
        if best is None or sol.cost < best.cost:
            best = sol
    return best.x


def consistency(groups: dict) -> list:
    """Fit all groups together, then with each group left out. Every fit
    reports chi2 per point for every group, including the one left out.
    `groups` maps a label to a list of datasets."""
    labels = list(groups)
    rows = []
    for leave in [None] + labels:
        use = [d for k in labels if k != leave for d in groups[k]]
        x = fit(use)
        chi = chi2_per_point(x, [d for k in labels for d in groups[k]])
        per_group = {k: float(np.mean([chi[d.name] for d in groups[k]]))
                     for k in labels}
        rows.append({"left_out": leave, "x": x, "chi2": per_group,
                     "consts": params.consts(x)})
    return rows


def log_post(x, datasets):
    if np.any(x < params.LO) or np.any(x > params.HI):
        return -np.inf
    r = residuals(x, datasets)
    return -0.5 * float(r @ r)


def mcmc(datasets, x0, n_walkers=40, n_steps=2500, burn=800, seed=7):
    """Goodman & Weare (2010) affine-invariant ensemble sampler (stretch
    move, a = 2). Returns the post-burn chain (walkers stacked) and the
    acceptance fraction."""
    _check(datasets)
    rng = np.random.default_rng(seed)
    d = len(x0)
    w = np.clip(x0 + 1e-3 * rng.standard_normal((n_walkers, d)),
                params.LO + 1e-9, params.HI - 1e-9)
    lp = np.array([log_post(v, datasets) for v in w])
    chain, acc, half = [], 0, n_walkers // 2
    for step in range(n_steps):
        for s in (0, 1):
            idx = np.arange(s * half, (s + 1) * half)
            other = np.arange((1 - s) * half, (2 - s) * half)
            z = ((rng.random(half) * (math.sqrt(2.0) - math.sqrt(0.5))
                  + math.sqrt(0.5)) ** 2)
            partners = w[rng.choice(other, half)]
            prop = partners + z[:, None] * (w[idx] - partners)
            lp_new = np.array([log_post(v, datasets) for v in prop])
            ok = np.log(rng.random(half)) < (d - 1) * np.log(z) + lp_new \
                - lp[idx]
            w[idx[ok]] = prop[ok]
            lp[idx[ok]] = lp_new[ok]
            acc += int(ok.sum())
        if step >= burn:
            chain.append(w.copy())
    return np.concatenate(chain), acc / (n_steps * n_walkers)


def summarize(chain, n_share=2000) -> dict:
    """Median and 16th/84th percentiles of every channel, and of the lasing
    B share at 1 atm (300 K)."""
    q = np.percentile(np.exp(chain), [16, 50, 84], axis=0)
    out = {n: tuple(q[:, i]) for i, n in enumerate(params.NAMES)}
    sub = chain[:: max(1, len(chain) // n_share)]
    shares = [forward.b_share_1atm(params.consts(x)) for x in sub]
    out["b_share_1atm"] = tuple(np.percentile(shares, [16, 50, 84]))
    return out
