"""A neural emulator of LIGHT_Kr's swarm solver, and neural posterior
estimation of the cross-section factors (numpy only, as ml/nn.py).

The emulator. The true solver takes 0.01 to 0.2 s a point and a fit
evaluates hundreds of points per trial, so the fits run on a network
trained to reproduce it:
  inputs   the nine ln-factors z, ln(E/N), Tg/300 K;
  outputs  the change from Biagi's set (z = 0) of ln W, ln D_T/mu and
           ln alpha/N at that E/N and Tg (alpha/N from ALPHA_FROM_TD up);
  data     the true solver at NODES for a Latin hypercube of z in the box
           [-BOX, BOX]^9 (simulate), each design's base taken at z = 0;
  check    held-out designs: the rms of its error per output (the
           emulator's error must sit well inside the data's).
EmuModel puts it on a dataset's points: the true solver's base value at
each measured E/N times exp(the emulated change).

Neural posterior estimation (as ml/nn.py): z drawn from the prior box,
the train datasets simulated through the emulator with their noise, and a
network with a Gaussian head reads z off the observations. Like ml/nn.py
it cannot tell conflicting datasets apart, so it runs on the consistency
check's train set only.
"""
from __future__ import annotations

import hashlib
import inspect
import math
from pathlib import Path

import numpy as np

from . import kr_light, kr_swarm
from .nn import MLP as GaussMLP

BOX = kr_swarm.BOX
ALPHA_FROM_TD = 12.0
RESULTS = Path(__file__).resolve().parent.parent / "results"


def nodes() -> list:
    """(E/N [Td], Tg [K]) where the true solver is sampled: 300 K from
    0.0007 to 330 Td, 195 K and 368 K where Pack's data are."""
    n300 = np.geomspace(7e-4, 330.0, 66)
    n195 = np.geomspace(1.1e-3, 0.1, 10)
    n368 = np.geomspace(1.4e-3, 7e-3, 5)
    return ([(float(e), 300.0) for e in n300] +
            [(float(e), 195.0) for e in n195] +
            [(float(e), 368.0) for e in n368])


def lhs(n, rng) -> np.ndarray:
    """A Latin hypercube in [-BOX, BOX]^NP."""
    u = (np.argsort(rng.random((kr_swarm.NP, n)), axis=1).T
         + rng.random((n, kr_swarm.NP))) / n
    return -BOX + 2.0 * BOX * u


def code_key() -> str:
    """What the simulations depend on: LIGHT_Kr's solver, swarm and factor
    code, Biagi's file, the knots, the nodes and the grid."""
    xs = kr_light.module("xsec")
    h = hashlib.md5()
    for name in ("boltzmann", "swarm"):
        h.update(Path(kr_light.module(name).__file__).read_bytes())
    for obj in (xs.Factors, xs.apply_factors, xs.exc_group, xs.parse_lxcat):
        h.update(inspect.getsource(obj).encode())
    h.update(xs.md5_of(xs.DATA_DIR / xs.FILES["Kr"][0]).encode())
    h.update(repr((kr_swarm.KNOTS_EV, kr_swarm.CELLS, nodes(), BOX,
                   kr_light.LIGHT_KR_VERSION)).encode())
    return h.hexdigest()[:12]


def simulate(n_design, seed, log=print) -> dict:
    """The true solver at every node for n_design z (plus z = 0 first),
    cached under ML/results/ by code_key and the seed."""
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / f"kr_sims_{code_key()}_{n_design}_{seed}.npz"
    if path.is_file():
        d = np.load(path)
        return {k: d[k] for k in d.files}
    rng = np.random.default_rng(seed)
    zs = np.vstack([np.zeros(kr_swarm.NP), lhs(n_design, rng)])
    pts = nodes()
    log(f"  simulating {len(zs)} designs x {len(pts)} nodes on the true "
        f"solver ({kr_swarm.WORKERS} workers) ...")
    y = kr_swarm.solve_many(zs, pts)
    out = {"z": zs, "en": np.array([p[0] for p in pts]),
           "tg": np.array([p[1] for p in pts]), "y": y}
    np.savez(path, **out)
    return out


def _xy(sims):
    """Training pairs: inputs (z, ln E/N, Tg/300), targets the change of
    ln (W, D_T/mu, alpha/N) from z = 0, and the alpha mask."""
    z, en, tg, y = sims["z"], sims["en"], sims["tg"], sims["y"]
    base = np.log(y[0])                         # (nodes, 3)
    X, Y, M = [], [], []
    for i in range(1, len(z)):
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.log(y[i]) - base
        X.append(np.column_stack([np.repeat(z[i][None, :], len(en), 0),
                                  np.log(en), tg / 300.0]))
        Y.append(d[:, [0, 2, 1]])               # W, D_T/mu, alpha/N
        M.append(np.column_stack([np.ones(len(en)), np.ones(len(en)),
                                  (en >= ALPHA_FROM_TD).astype(float)]))
    Y, M = np.vstack(Y), np.vstack(M)
    Y = np.where(M > 0, np.nan_to_num(Y), 0.0)
    return np.vstack(X), Y, M


class Regressor:
    """Two (or three) hidden tanh layers, linear outputs, Adam on a masked
    mean squared error."""

    def __init__(self, n_in, n_out, width, depth, rng):
        sizes = [n_in] + [width] * depth + [n_out]
        self.W = [rng.standard_normal((a, b)) * math.sqrt(1.0 / a)
                  for a, b in zip(sizes[:-1], sizes[1:])]
        self.b = [np.zeros(b) for b in sizes[1:]]
        self.m = [np.zeros_like(w) for w in self.W + self.b]
        self.v = [np.zeros_like(w) for w in self.W + self.b]
        self.t = 0

    def forward(self, X):
        with np.errstate(all="ignore"):
            hs = [X]
            for W, b in zip(self.W[:-1], self.b[:-1]):
                hs.append(np.tanh(hs[-1] @ W + b))
            return hs[-1] @ self.W[-1] + self.b[-1], hs

    def step(self, X, Y, M, lr):
        with np.errstate(all="ignore"):
            out, hs = self.forward(X)
            g = 2.0 * (out - Y) * M / max(M.sum(), 1.0)
            gW, gb = [], []
            for k in range(len(self.W) - 1, -1, -1):
                gW.insert(0, hs[k].T @ g)
                gb.insert(0, g.sum(0))
                if k:
                    g = (g @ self.W[k].T) * (1.0 - hs[k] ** 2)
            self.t += 1
            params = self.W + self.b
            for i, (p, gr) in enumerate(zip(params, gW + gb)):
                self.m[i] = 0.9 * self.m[i] + 0.1 * gr
                self.v[i] = 0.999 * self.v[i] + 0.001 * gr * gr
                mh = self.m[i] / (1 - 0.9 ** self.t)
                vh = self.v[i] / (1 - 0.999 ** self.t)
                p -= lr * mh / (np.sqrt(vh) + 1e-8)


class Emulator:
    def __init__(self, net, xm, xs, ys):
        self.net, self.xm, self.xs, self.ys = net, xm, xs, ys

    def delta(self, X) -> np.ndarray:
        out, _ = self.net.forward((X - self.xm) / self.xs)
        return out * self.ys


def train_emulator(sims, n_val=100, width=160, depth=3, epochs=300,
                   seed=5, log=print) -> tuple:
    """The emulator and its held-out rms error per output (ln units)."""
    rng = np.random.default_rng(seed)
    n = len(sims["z"]) - 1
    val = set(range(n - n_val, n))
    X, Y, M = _xy(sims)
    per = len(sims["en"])
    rows = np.arange(len(X))
    is_val = np.isin(rows // per, list(val))
    Xt, Yt, Mt = X[~is_val], Y[~is_val], M[~is_val]
    xm, xs = Xt.mean(0), Xt.std(0) + 1e-9
    ys = np.sqrt((Yt ** 2 * Mt).sum(0) / Mt.sum(0)) + 1e-6
    net = Regressor(X.shape[1], 3, width, depth, rng)
    Xn, Yn = (Xt - xm) / xs, Yt / ys
    for ep in range(epochs):
        lr = 2e-3 if ep < epochs * 0.5 else (6e-4 if ep < epochs * 0.8
                                              else 1.5e-4)
        order = rng.permutation(len(Xn))
        for i in range(0, len(Xn), 512):
            b = order[i:i + 512]
            net.step(Xn[b], Yn[b], Mt[b], lr)
        if log and (ep + 1) % 50 == 0:
            emu = Emulator(net, xm, xs, ys)
            err = _rms(emu, X[is_val], Y[is_val], M[is_val])
            log(f"    epoch {ep + 1}: held-out rms (W, D_T/mu, alpha/N) "
                + ", ".join(f"{e:.4f}" for e in err))
    assert all(np.all(np.isfinite(w)) for w in net.W), "training diverged"
    emu = Emulator(net, xm, xs, ys)
    return emu, _rms(emu, X[is_val], Y[is_val], M[is_val])


def _rms(emu, X, Y, M):
    d = (emu.delta(X) - Y) * M
    return np.sqrt((d ** 2).sum(0) / M.sum(0))


class EmuModel:
    """The forward model on the emulator: the true solver's z = 0 value at
    each measured point, with the Monte Carlo's correction
    (kr_swarm.mc_ratio), times exp(the emulated change)."""

    def __init__(self, emu, datasets, base_model):
        self.emu = emu
        self.base = {d.name: base_model(np.zeros(kr_swarm.NP), d,
                                        corrected=False)
                     * kr_swarm.mc_ratio(d.kind, d.x) for d in datasets}
        self.X = {d.name: np.column_stack([np.log(d.x),
                                           np.full(len(d.x), d.tg / 300.0)])
                  for d in datasets}

    def __call__(self, z, d):
        xe = self.X[d.name]
        X = np.column_stack([np.repeat(np.asarray(z, dtype=float)[None, :],
                                       len(xe), 0), xe])
        col = {"W": 0, "DTmu": 1, "alphaN": 2}[d.kind]
        return self.base[d.name] * np.exp(self.emu.delta(X)[:, col])


# ------------------------------------------------ neural posterior estimation
def npe(datasets, model, sig, n_train=40000, n_val=4000, epochs=60,
        width=128, seed=3, free=None) -> dict:
    """ml/nn.py's estimator on the swarm: the free components of z uniform
    in the box (the others 0), every train dataset simulated through
    `model` with its noise (sigma_eff), and a Gaussian-head network
    trained to read them off. Returns its held-out R^2 and +-1 sigma
    coverage per free factor and its answer on the measured data (mean and
    sd of z; the fixed components 0, their R^2 and coverage NaN)."""
    f = np.asarray(kr_swarm.FREE if free is None else free)
    rng = np.random.default_rng(seed)
    sg = np.concatenate([sig[d.name] for d in datasets])
    obs = np.concatenate([np.log(d.y) for d in datasets])

    def sim(n):
        Z = np.zeros((n, kr_swarm.NP))
        Z[:, f] = -BOX + 2 * BOX * rng.random((n, int(f.sum())))
        O = np.array([np.concatenate([np.log(model(z, d)) for d in datasets])
                      for z in Z])
        return O + sg * rng.standard_normal(O.shape), Z[:, f]

    Xtr, Ytr = sim(n_train)
    Xva, Yva = sim(n_val)
    xm, xs = Xtr.mean(0), Xtr.std(0) + 1e-9
    ym, ys = Ytr.mean(0), Ytr.std(0)
    net = GaussMLP(Xtr.shape[1], Ytr.shape[1], width, rng)
    Xn, Yn = (Xtr - xm) / xs, (Ytr - ym) / ys
    for ep in range(epochs):
        lr = 2e-3 if ep < epochs // 2 else 5e-4
        order = rng.permutation(len(Xn))
        for i in range(0, len(Xn), 256):
            b = order[i:i + 256]
            net.step(Xn[b], Yn[b], lr)
    assert all(np.all(np.isfinite(w)) for w in net.p), "training diverged"
    mu, lv, _ = net.forward((Xva - xm) / xs)
    pred, sd = mu * ys + ym, np.exp(0.5 * lv) * ys
    mu0, lv0, _ = net.forward(((obs - xm) / xs)[None, :])

    def full(v, fill):
        out = np.full(kr_swarm.NP, fill, dtype=float)
        out[f] = v
        return out
    return {"r2": full(1 - ((pred - Yva) ** 2).mean(0) / Yva.var(0), np.nan),
            "coverage": full((np.abs(pred - Yva) < sd).mean(0), np.nan),
            "z": full(mu0[0] * ys + ym, 0.0),
            "z_sd": full(np.exp(0.5 * lv0[0]) * ys, 0.0)}
