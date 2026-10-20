"""The krypton swarm training, end to end (`python3 ML/run.py kr-train`).

1. The datasets (kr_data) and the true solver's base and slope at every
   measured point (kr_swarm: the error model's sigma_eff).
2. The emulator (kr_nn), trained on the true solver over the prior box and
   checked on held-out designs.
3. The consistency check on the candidates (kr_swarm.consistency): every
   dataset the forward model can claim (no 'check' set), pruned until each
   kept one agrees with the fit of the others within its errors.
4. Least squares, the MCMC and neural posterior estimation on that train
   set; the answer re-checked on the true solver at every dataset.
5. The record LIGHT_Kr pins (kr_swarm.factors_record).
Everything is cached under ML/results/ by the simulations' code key.
"""
from __future__ import annotations

import json
import pickle
import time

import numpy as np

from . import kr_data, kr_nn, kr_swarm

N_DESIGN, SIM_SEED = 2000, 1
RESULTS = kr_nn.RESULTS


def _cached(name, make):
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / f"kr_{name}_{kr_nn.code_key()}.pkl"
    if path.is_file():
        with open(path, "rb") as fh:
            return pickle.load(fh)
    obj = make()
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)
    return obj


def setup(log=print) -> dict:
    """Datasets, true model, sigma_eff, simulations and the emulator."""
    t0 = time.time()
    ds = kr_data.load()
    true = kr_swarm.TrueModel(ds)
    base = _cached("base", lambda: true.table(np.zeros(kr_swarm.NP)))
    true._cache[tuple(np.zeros(kr_swarm.NP))] = base
    slope = _cached("slopes", lambda: kr_swarm.slopes(true, ds))
    sig = {d.name: kr_swarm.effective_sigma(d, slope[d.name]) for d in ds}
    sims = kr_nn.simulate(N_DESIGN, SIM_SEED, log=log)

    def make_emu():
        emu, err = kr_nn.train_emulator(sims, log=log)
        return {"emu": emu, "err": err}
    e = _cached("emulator", make_emu)
    model = kr_nn.EmuModel(e["emu"], ds, true)
    log(f"  setup {time.time() - t0:.0f} s; emulator held-out rms (ln W, "
        "ln D_T/mu, ln alpha/N): " + ", ".join(f"{v:.4f}" for v in e["err"]))
    return {"ds": ds, "true": true, "sig": sig, "model": model,
            "emu_err": e["err"], "sims": sims}


def candidates(ds, low_field=False) -> list:
    """The consistency check's candidates: every dataset but the check
    sets, D_T/mu (held out: a test of the trained momentum transfer that
    the training never sees) and, unless low_field (a diagnostic), the
    low-field alpha/N."""
    return [d for d in ds if d.name not in kr_data.CHECK and
            d.kind != "DTmu" and
            (low_field or d.name not in kr_data.HELD_BY_RULE)]


def run_consistency(s, threshold=3.0, low_field=False, log=print) -> dict:
    cand = candidates(s["ds"], low_field)
    for d in cand:
        d.role = "train"            # the check fits them; roles follow it
    tag = (f"consistency_joint_{threshold:g}_"
           + "".join("1" if x else "0" for x in kr_swarm.FREE)
           + ("_low" if low_field else ""))
    return _cached(tag, lambda: kr_swarm.consistency(
        cand, s["model"], s["sig"], threshold=threshold, log=log))


def run_fit(s, train_names, n_steps=3000, log=print) -> dict:
    by = {d.name: d for d in s["ds"]}
    train = [by[n] for n in train_names]
    for d in s["ds"]:
        d.role = ("train" if d.name in train_names else
                  "check" if d.name in kr_data.CHECK else "held_out")
    z, cov = kr_swarm.fit(train, s["model"], s["sig"])
    chain, acc = kr_swarm.mcmc(train, s["model"], s["sig"], z,
                               n_steps=n_steps)
    post = kr_swarm.summarize(chain)
    zmed = np.median(chain, axis=0)
    nn = kr_nn.npe(train, s["model"], s["sig"])
    # the answer on the true solver, every dataset
    chi_true = kr_swarm.chi2_per_point(z, s["ds"], s["true"], s["sig"])
    chi_base = kr_swarm.chi2_per_point(np.zeros(kr_swarm.NP), s["ds"],
                                       s["true"], s["sig"])
    chi_emu = kr_swarm.chi2_per_point(z, s["ds"], s["model"], s["sig"])
    # 40 draws, evenly spaced through the post-burn chain: what a document
    # propagates (the chain is walkers stacked per step, so a stride that is
    # a multiple of the walker count would take one walker's path)
    idx = np.linspace(0, len(chain) - 1, 40).round().astype(int)
    return {"z": z, "cov": cov, "chain_acc": acc, "post": post,
            "z_median": zmed, "z_chain_sd": chain.std(0), "nn": nn,
            "draws": chain[idx], "chain_cov": np.cov(chain.T),
            "chi_true": chi_true, "chi_base": chi_base, "chi_emu": chi_emu,
            "train": list(train_names)}


VARIANTS = {
    "all nine": (True,) * 9,
    "knots below 0.6 eV fixed": (False, False) + (True,) * 7,
    "knots below 1 eV fixed": (False, False, False) + (True,) * 6,
    "excitation only": (False,) * 6 + (True,) * 3,
}


def free_variants(s, train_names) -> dict:
    """Which factors to free: each variant's least squares on the train
    set, and every W and D_T/mu dataset's chi2 per point on the true
    solver (the evidence for kr_swarm.FREE), with Biagi's file beside."""
    by = {d.name: d for d in s["ds"]}
    train = [by[n] for n in train_names]
    for d in train:
        d.role = "train"
    look = [d for d in s["ds"] if d.kind in ("W", "DTmu")
            and d.name not in kr_data.CHECK] + [
        by[n] for n in train_names if by[n].kind == "alphaN"]
    out = {}
    for label, free in list(VARIANTS.items()) + [("Biagi's file", None)]:
        if free is None:
            z = np.zeros(kr_swarm.NP)
        else:
            z, _ = kr_swarm.fit(train, s["model"], s["sig"], free=free)
        chi = kr_swarm.chi2_per_point(z, look, s["true"], s["sig"])
        out[label] = {"free": None if free is None else list(free),
                      "factors": [float(v) for v in np.exp(z)],
                      "chi2": {k: round(v, 3) for k, v in chi.items()}}
    return out


def threshold_scan(s, thresholds=(2.0, 3.0, 4.0, 6.0)) -> list:
    """The train set and least-squares factors at other consistency
    thresholds (how much the choice of 3 decides)."""
    out = []
    for th in thresholds:
        c = kr_swarm.consistency(candidates(s["ds"]), s["model"], s["sig"],
                                 threshold=th, log=None)
        by = {d.name: d for d in s["ds"]}
        train = [by[n] for n in c["kept"]]
        for d in train:
            d.role = "train"
        z, cov = kr_swarm.fit(train, s["model"], s["sig"])
        out.append({"threshold": th, "train": c["kept"],
                    "factors": [float(v) for v in np.exp(z)],
                    "z_sd": [float(v) for v in np.sqrt(np.diag(cov))]})
    return out


def record(fit, note, s=None, consistency=None) -> dict:
    # W and alpha/N only: D_T/mu is held out and scored by the documents on
    # the Monte Carlo's low-field correction (LIGHT_Kr's swarm_mc.py)
    keep = {d.name for d in kr_data.load() if d.kind != "DTmu"}
    rec = kr_swarm.factors_record(
        fit["z"], fit["cov"], {k: list(v) for k, v in fit["post"].items()},
        fit["train"], {k: round(v, 4) for k, v in fit["chi_true"].items()
                       if k in keep}, note)
    rec["chi2_per_point_biagi"] = {k: round(v, 4)
                                   for k, v in fit["chi_base"].items()
                                   if k in keep}
    rec["mcmc"] = {"acceptance": round(float(fit["chain_acc"]), 3),
                   "median": [float(v) for v in np.exp(fit["z_median"])],
                   "z_sd": [float(v) for v in fit["z_chain_sd"]],
                   "z_cov": [[float(v) for v in row]
                             for row in fit["chain_cov"]],
                   "draws": [[float(v) for v in row]
                             for row in fit["draws"]]}
    def num(v):
        return None if not np.isfinite(v) else round(float(v), 3)
    rec["nn"] = {"factors": [float(v) for v in np.exp(fit["nn"]["z"])],
                 "z_sd": [float(v) for v in fit["nn"]["z_sd"]],
                 "r2": [num(v) for v in fit["nn"]["r2"]],
                 "coverage": [num(v) for v in fit["nn"]["coverage"]]}
    rec["free"] = [bool(v) for v in kr_swarm.FREE]
    if s is not None:
        rec["emulator"] = {"n_design": N_DESIGN, "box": kr_swarm.BOX,
                           "nodes": len(kr_nn.nodes()),
                           "rms_ln": [round(float(v), 5)
                                      for v in s["emu_err"]],
                           "code_key": kr_nn.code_key()}
        rec["threshold_scan"] = threshold_scan(s)
    if s is not None:
        rec["free_variants"] = free_variants(s, fit["train"])
    if consistency is not None:
        rec["consistency"] = {
            "threshold": 3.0, "kept": consistency["kept"],
            "dropped": [[n, round(v, 3)] for n, v in consistency["dropped"]],
            "rounds": len(consistency["rounds"])}
    rec["names"] = list(kr_swarm.PNAMES)
    return rec


def report(s, fit, log=print):
    """The fit's factors and every dataset's chi2 per point."""
    log("\nTrained factors (least squares; MCMC median [16, 84 %]; neural "
        "estimator mean +- sd):")
    for i, n in enumerate(kr_swarm.PNAMES):
        lo, mid, hi = fit["post"][n]
        log(f"  {n:7s} {np.exp(fit['z'][i]):6.3f}   {mid:6.3f} [{lo:.3f}, "
            f"{hi:.3f}]   {np.exp(fit['nn']['z'][i]):6.3f} +- "
            f"{fit['nn']['z_sd'][i]:.3f} (ln)")
    log(f"  MCMC acceptance {fit['chain_acc']:.2f}; network held-out R^2 "
        + ", ".join(f"{v:.2f}" for v in fit["nn"]["r2"]) + "; coverage "
        + ", ".join(f"{v:.2f}" for v in fit["nn"]["coverage"]))
    log("\nchi2 per point, Biagi's set -> trained (true solver; emulator):")
    for role in ("train", "held_out", "check"):
        log(f" {role}")
        for d in s["ds"]:
            if d.role != role:
                continue
            log(f"  {d.name:20s} {len(d.x):3d} pts  {fit['chi_base'][d.name]:8.2f}"
                f" -> {fit['chi_true'][d.name]:8.2f}   ({fit['chi_emu'][d.name]:.2f})")


def dumps(obj) -> str:
    return json.dumps(obj, indent=1, sort_keys=True, default=float)
