"""Neural posterior estimation: a network that reads the channels off the
data (numpy only).

Method:
  1. Draw channel vectors from the training priors below.
  2. Simulate each dataset's observables with the forward models, and add
     that dataset's noise.
  3. Train a two-layer tanh MLP with a Gaussian head (a mean and a variance
     per channel) to map observables to channels.
  4. Check it on held-out simulations (R^2 and +-1 sigma coverage), then
     apply it to the measured data.
A network reads constants off consistent data as well as the MCMC does. It
cannot tell which of several conflicting datasets is wrong, and fed them all
it returns a confident compromise. Check consistency first (fit.consistency).
"""
from __future__ import annotations

import math

import numpy as np

from . import forward, params

# training priors: wide, but inside the region where every observable is a
# measurable decay (a 1e12 s^-1 decay constant is not an experiment)
P_LO = np.log([1e-12, 0.3, 25e-9, 6.3e-9, 1e-14, 1e-14, 1e-32, 1e-32, 1e-10,
               1e-10])
P_HI = np.log([1e-9, 10.0, 120e-9, 7.5e-9, 1e-11, 1e-11, 3e-30, 3e-30, 1e-9,
               1e-9])


def _observe(x, datasets):
    c = params.consts(x)
    return np.concatenate([np.log(forward.predict(d, c)) for d in datasets])


def _sigma(datasets):
    s = []
    for d in datasets:
        frac = d.sigma if d.sigma_is_fraction else d.sigma / np.mean(d.values)
        s += [frac] * len(d.points)
    return np.array(s)


def simulate(n, datasets, rng):
    X, Y = [], []
    sig = _sigma(datasets)
    while len(X) < n:
        x = P_LO + (P_HI - P_LO) * rng.random(len(P_LO))
        o = _observe(x, datasets)
        if np.all(np.isfinite(o)):
            X.append(o + sig * rng.standard_normal(len(o)))
            Y.append(x)
    return np.array(X), np.array(Y)


class MLP:
    """Two hidden tanh layers and a Gaussian head, trained with Adam on the
    Gaussian negative log-likelihood."""

    def __init__(self, n_in, n_out, width, rng):
        def init(a, b):
            return rng.standard_normal((a, b)) * math.sqrt(1.0 / a)
        self.p = [init(n_in, width), np.zeros(width), init(width, width),
                  np.zeros(width), init(width, 2 * n_out), np.zeros(2 * n_out)]
        self.m = [np.zeros_like(w) for w in self.p]
        self.v = [np.zeros_like(w) for w in self.p]
        self.t, self.n_out = 0, n_out

    # numpy 2 on macOS Accelerate raises spurious divide/overflow warnings in
    # matmuls whose operands hold exact zeros (LIGHT Known-bad list). They are
    # silenced here, and train() asserts every number is finite.
    def forward(self, X):
        with np.errstate(all="ignore"):
            return self._forward(X)

    def _forward(self, X):
        W1, b1, W2, b2, W3, b3 = self.p
        h1 = np.tanh(X @ W1 + b1)
        h2 = np.tanh(h1 @ W2 + b2)
        out = h2 @ W3 + b3
        return (out[:, :self.n_out], np.clip(out[:, self.n_out:], -12, 6),
                (X, h1, h2))

    def step(self, X, Y, lr):
        with np.errstate(all="ignore"):
            W1, b1, W2, b2, W3, b3 = self.p
            mu, lv, (x0, h1, h2) = self._forward(X)
            inv, n = np.exp(-lv), len(X)
            g_out = np.concatenate([-(Y - mu) * inv / n,
                                    0.5 * (1.0 - (Y - mu) ** 2 * inv) / n],
                                   axis=1)
            g_h2 = (g_out @ W3.T) * (1 - h2 ** 2)
            g_h1 = (g_h2 @ W2.T) * (1 - h1 ** 2)
            grads = [x0.T @ g_h1, g_h1.sum(0), h1.T @ g_h2, g_h2.sum(0),
                     h2.T @ g_out, g_out.sum(0)]
            self.t += 1
            for i, g in enumerate(grads):
                self.m[i] = 0.9 * self.m[i] + 0.1 * g
                self.v[i] = 0.999 * self.v[i] + 0.001 * g * g
                mh = self.m[i] / (1 - 0.9 ** self.t)
                vh = self.v[i] / (1 - 0.999 ** self.t)
                self.p[i] -= lr * mh / (np.sqrt(vh) + 1e-8)


def train(datasets, n_train=30000, n_val=3000, epochs=80, width=128, seed=3):
    rng = np.random.default_rng(seed)
    Xtr, Ytr = simulate(n_train, datasets, rng)
    Xva, Yva = simulate(n_val, datasets, rng)
    xm, xs = Xtr.mean(0), Xtr.std(0) + 1e-9
    ym, ys = Ytr.mean(0), Ytr.std(0)
    net = MLP(Xtr.shape[1], Ytr.shape[1], width, rng)
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
    return {"net": net, "xm": xm, "xs": xs, "ym": ym, "ys": ys,
            "r2": 1 - ((pred - Yva) ** 2).mean(0) / Yva.var(0),
            "coverage": (np.abs(pred - Yva) < sd).mean(0),
            "datasets": list(datasets)}


def apply(model):
    """The network's answer (natural-log mean and sd) for the measured data."""
    x_obs = np.concatenate([np.log(d.values) for d in model["datasets"]])
    mu, lv, _ = model["net"].forward(((x_obs - model["xm"]) / model["xs"])
                                     [None, :])
    return mu[0] * model["ys"] + model["ym"], np.exp(0.5 * lv[0]) * model["ys"]
