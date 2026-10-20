"""Training LIGHT_Kr's cross sections on krypton's swarm data.

The owner (2026-10-14): "Let us try to improve like we did for LIGHT/, we
can also use ML to train it to improve the accuracy". LIGHT's channels
were fitted to their own measurements (ml/fit.py, ml/nn.py); LIGHT_Kr's
largest scored error is the electron swarm, and swarm data measure the
cross sections themselves, so they may train them (ml/kr_data.py).

What is trained (xsec.Factors of LIGHT_Kr v0_3, natural logs z):
  m_<e>   the elastic momentum transfer's factor at six knots, 0.03 to 20
          eV, log-linear in energy between them (held beyond the ends);
  x_1s, x_2p, x_hi  one factor on each group of Biagi's excitation levels
          (the four 5s levels, the ten 5p levels, the 37 others).
The ionization cross section is not trained: beam experiments measure it,
and in a swarm its scale trades against the excitation's. The prior is
broad, ln-factor sd PRIOR_SD (a factor 1.65 at one sigma), so the data
decide; where the posterior lands is reported against it.

The forward model is LIGHT_Kr's own two-term solver as a swarm (swarm.py:
W, alpha/N = k_ion/W, D_T/mu), run on LIGHT_Kr's processes with the
factors (xsec.using). The error model is each dataset's fractional sigma
(kr_data) and an E/N uncertainty EN_SIGMA through the model's own slope,
    sigma_eff^2 = sigma^2 + (EN_SIGMA d ln y / d ln(E/N))^2,
which matters where W climbs out of the Ramsauer minimum and for alpha/N
at low field. Residuals are ln(model/measured)/sigma_eff, chi2 per point
near 1 means agreement within the stated errors.

The true solver costs 0.01 to 0.2 s a point, so the fits run on a neural
emulator of it (ml/kr_nn.py) and every answer is re-checked on the true
solver. Fits: least squares (scipy), the consistency check (leave one out,
pruning the worst until every train dataset agrees within its errors),
an ensemble MCMC (Goodman & Weare, as ml/fit.py) and neural posterior
estimation (ml/kr_nn.py).
"""
from __future__ import annotations

import json
import math
import os
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

from . import kr_data, kr_light

KNOTS_EV = (0.03, 0.2, 0.6, 2.0, 6.0, 20.0)
PNAMES = tuple(f"m_{k:g}" for k in KNOTS_EV) + ("x_1s", "x_2p", "x_hi")
NP = len(PNAMES)
# Which factors train (2026-10-15). With all nine free the knots at 0.03
# and 0.2 eV fit Hunter 1988's drift velocity below its stated error (chi2
# per point 0.79 -> 0.13) and move the momentum transfer away from every
# other low-field dataset: Koizumi's characteristic energy 2.5 -> 8.1,
# Pack's 195-K drift velocities 3.1 -> 5.5 and 6.1 -> 8.7. With those two
# knots held at Biagi's value the train sets fit as well (all within the
# threshold, Hunter 0.70) and the same datasets score 0.84, 2.8 and 5.3.
# So the momentum transfer below 0.2 eV stays Biagi's; Koizumi's and the
# dropped drift velocities informed this choice (validation), and the
# beam partition and the low-field alpha/N, which no choice saw, test it.
FREE = (False, False, True, True, True, True, True, True, True)
PRIOR_SD = 0.5
BOX = 0.6           # |z| the emulator is trained over (kr_nn); fits stay in
EN_SIGMA = 0.02
CELLS = 1000
WORKERS = int(os.environ.get("KR_WORKERS", "6"))
KINDS = ("W", "alphaN", "DTmu")


def factors(z):
    """xsec.Factors from the ln-factor vector z (UNIT at z = 0)."""
    xs = kr_light.module("xsec")
    z = np.asarray(z, dtype=float)
    return xs.Factors(KNOTS_EV, tuple(float(v) for v in np.exp(z[:6])),
                      tuple(float(v) for v in np.exp(z[6:9])), 1.0)


# ------------------------------------------------- what the two-term leaves out
# The two-term swarm gives the flux drift velocity, k_i/W and its own
# transverse diffusion; a drift tube measures the bulk drift velocity, a
# Townsend experiment the SST coefficient, and the two-term D_T is far off
# above a few Td. LIGHT_Kr's swarm_mc.py holds LoKI-MC's ratios of each to
# the two-term's (on Biagi's file; the trained set's within 1.3 %), and the
# forward model multiplies the two-term swarm by them.
def mc_ratio(kind, en) -> np.ndarray:
    """The Monte Carlo's value over the two-term's for one kind at E/N."""
    return kr_light.module("swarm_mc").ratio(kind, en)


# ------------------------------------------------------------ true solver
def _solve(args):
    """One swarm point on the true solver: (W, alpha/N, D_T/mu)."""
    z, en, tg = args
    xs = kr_light.module("xsec")
    sw = kr_light.module("swarm")
    with xs.using(factors(z)):
        d = sw.point(en, tg, cells=CELLS)
    return d["W"], d["alphaN"], d["DTmu"]


def points(datasets) -> list:
    """The distinct (E/N, Tg) of the datasets."""
    return sorted({(float(x), float(d.tg)) for d in datasets for x in d.x})


def solve_many(zs, pts, workers=WORKERS) -> np.ndarray:
    """(len(zs), len(pts), 3) on the true solver, in parallel."""
    jobs = [(tuple(map(float, z)), en, tg) for z in zs for en, tg in pts]
    if workers > 1:
        with Pool(workers) as pool:
            out = pool.map(_solve, jobs, chunksize=4)
    else:
        out = [_solve(j) for j in jobs]
    return np.array(out, dtype=float).reshape(len(zs), len(pts), 3)


class TrueModel:
    """The forward model on the true solver, cached by z."""

    def __init__(self, datasets, workers=WORKERS):
        self.datasets = list(datasets)
        self.pts = points(self.datasets)
        self.index = {p: i for i, p in enumerate(self.pts)}
        self.workers = workers
        self._cache = {}

    def table(self, z):
        key = tuple(np.round(np.asarray(z, dtype=float), 12))
        if key not in self._cache:
            self._cache[key] = solve_many([z], self.pts, self.workers)[0]
        return self._cache[key]

    def __call__(self, z, d, corrected=True):
        t = self.table(z)
        col = KINDS.index(d.kind)
        y = np.array([t[self.index[(float(x), float(d.tg))], col]
                      for x in d.x])
        return y * mc_ratio(d.kind, d.x) if corrected else y


# ------------------------------------------------------------ error model
def slopes(model, datasets, h=0.01) -> dict:
    """d ln y / d ln(E/N) of the model at z = 0, by a central difference on
    the true solver (one solve each side of every point)."""
    out = {}
    pts = points(datasets)
    lo = [(en * math.exp(-h), tg) for en, tg in pts]
    hi = [(en * math.exp(h), tg) for en, tg in pts]
    z0 = np.zeros(NP)
    a = solve_many([z0], lo, model.workers if hasattr(model, "workers")
                   else WORKERS)[0]
    b = solve_many([z0], hi, model.workers if hasattr(model, "workers")
                   else WORKERS)[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        s = (np.log(b) - np.log(a)) / (2 * h)
    idx = {p: i for i, p in enumerate(pts)}
    for d in datasets:
        col = KINDS.index(d.kind)
        out[d.name] = np.array([s[idx[(float(x), float(d.tg))], col]
                                for x in d.x])
    return out


def effective_sigma(d, slope) -> np.ndarray:
    return np.sqrt(d.sigma ** 2 + (EN_SIGMA * np.nan_to_num(slope)) ** 2)


# ------------------------------------------------------------ fits
def residuals(z, datasets, model, sig) -> np.ndarray:
    r = [np.log(model(z, d) / d.y) / sig[d.name] for d in datasets]
    r.append(np.asarray(z) / PRIOR_SD)
    return np.concatenate(r)


def chi2_per_point(z, datasets, model, sig) -> dict:
    return {d.name: float(np.mean((np.log(model(z, d) / d.y)
                                   / sig[d.name]) ** 2)) for d in datasets}


def embed(zs, free=FREE) -> np.ndarray:
    """The full z with the free components zs and the others at 0."""
    z = np.zeros(NP)
    z[np.asarray(free)] = zs
    return z


def fit(datasets, model, sig, z0=None, bound=BOX, free=FREE):
    """Least squares in the free components of z with the prior, inside
    the emulator's box; returns (z, covariance), the fixed components 0
    with no variance."""
    for d in datasets:
        if d.role == "check":
            raise ValueError(f"{d.name}: a check dataset never trains")
    f = np.asarray(free)
    zs0 = (np.zeros(NP) if z0 is None else np.asarray(z0, dtype=float))[f]
    sol = least_squares(lambda zs: residuals(embed(zs, f), datasets, model,
                                             sig), zs0,
                        bounds=(-bound, bound), x_scale=0.1,
                        diff_step=1e-4, max_nfev=400)
    J = sol.jac
    cov = np.zeros((NP, NP))
    cov[np.ix_(f, f)] = np.linalg.pinv(J.T @ J)
    return embed(sol.x, f), cov


def consistency(candidates, model, sig, threshold=3.0, log=print) -> dict:
    """The train set: fit every candidate together and score each against
    that joint fit (its chi2 per point there); drop the worst while it
    exceeds the threshold, refit, and repeat. Every round also records each
    kept set's leave-one-out chi2 per point (the set against the fit of the
    others), shown but not used: a set that alone pins some factor (the
    momentum transfer at a few eV, say) predicts poorly when left out
    however well it agrees, so leaving one out reads 'informative' as
    'inconsistent' (2026-10-14: it dropped Hunter 1988's drift velocity)."""
    keep = list(candidates)
    rounds, dropped = [], []
    while True:
        z_all, _ = fit(keep, model, sig)
        joint = chi2_per_point(z_all, keep, model, sig)
        loo = {}
        for d in keep:
            others = [e for e in keep if e is not d]
            z_o, _ = fit(others, model, sig, z0=z_all)
            loo[d.name] = chi2_per_point(z_o, [d], model, sig)[d.name]
        rounds.append({"kept": [d.name for d in keep], "joint": joint,
                       "loo": loo, "z": z_all.tolist()})
        worst = max(keep, key=lambda d: joint[d.name])
        if log:
            log(f"  round {len(rounds)}: {len(keep)} sets, worst "
                f"{worst.name} chi2/pt {joint[worst.name]:.2f} at the "
                f"joint fit")
        if joint[worst.name] <= threshold or len(keep) <= 2:
            break
        keep.remove(worst)
        dropped.append((worst.name, joint[worst.name]))
    return {"kept": [d.name for d in keep], "dropped": dropped,
            "rounds": rounds}


def log_post(z, datasets, model, sig, bound=BOX):
    if np.any(np.abs(z) > bound):
        return -np.inf
    r = residuals(z, datasets, model, sig)
    return -0.5 * float(r @ r)


def mcmc(datasets, model, sig, z0, n_walkers=40, n_steps=3000, burn=1000,
         seed=11, free=FREE):
    """Goodman & Weare (2010) stretch-move ensemble sampler (a = 2), as
    ml/fit.py, over the free components. Returns the post-burn chain (full
    z, the fixed components 0) and the acceptance fraction."""
    f = np.asarray(free)
    nf = int(f.sum())
    rng = np.random.default_rng(seed)
    w = np.asarray(z0, dtype=float)[f] + 1e-3 * rng.standard_normal(
        (n_walkers, nf))

    def lpost(v):
        return log_post(embed(v, f), datasets, model, sig)
    lp = np.array([lpost(v) for v in w])
    chain, acc, half = [], 0, n_walkers // 2
    for step in range(n_steps):
        for s in (0, 1):
            idx = np.arange(s * half, (s + 1) * half)
            other = np.arange((1 - s) * half, (2 - s) * half)
            zz = ((rng.random(half) * (math.sqrt(2.0) - math.sqrt(0.5))
                   + math.sqrt(0.5)) ** 2)
            partners = w[rng.choice(other, half)]
            prop = partners + zz[:, None] * (w[idx] - partners)
            lp_new = np.array([lpost(v) for v in prop])
            ok = np.log(rng.random(half)) < (nf - 1) * np.log(zz) + lp_new \
                - lp[idx]
            w[idx[ok]] = prop[ok]
            lp[idx[ok]] = lp_new[ok]
            acc += int(ok.sum())
        if step >= burn:
            chain.append(np.array([embed(v, f) for v in w]))
    return np.concatenate(chain), acc / (n_steps * n_walkers)


def summarize(chain) -> dict:
    """Median and 16th/84th percentiles of every factor (not its log)."""
    q = np.percentile(np.exp(chain), [16, 50, 84], axis=0)
    return {n: tuple(float(v) for v in q[:, i]) for i, n in enumerate(PNAMES)}


def factors_record(z, cov, chain_summary, train, chi2, note) -> dict:
    """The trained file's content (xsec.read_factors reads the first four
    keys; the rest is the fit's record)."""
    f = factors(z)
    return {"m_knots_eV": list(KNOTS_EV), "m": list(f.m),
            "exc": dict(zip(("1s", "2p", "hi"), f.exc)), "ion": 1.0,
            "z": [float(v) for v in z],
            "z_sd": [float(v) for v in np.sqrt(np.diag(cov))],
            "z_cov": [[float(v) for v in row] for row in cov],
            "posterior": chain_summary, "train": list(train),
            "chi2_per_point": chi2, "note": note}


def write_trained(record, path) -> str:
    """Writes the trained factors file; returns its md5 (to pin in LIGHT_Kr's
    xsec.TRAINED_FILE)."""
    text = json.dumps(record, indent=1, sort_keys=True) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return kr_data.md5_of(path)
