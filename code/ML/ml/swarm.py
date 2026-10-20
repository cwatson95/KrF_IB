"""Electron swarm measurements: the data that may train LIGHT's cross sections.

Swarm data (drift velocity, ionization and attachment coefficients, rate
constants in a buffer gas) are measurements of the cross sections
themselves, so under the predictive-only policy they may fit them. The
forward model is the two-term Boltzmann solver vendored in ml/boltzmann
(the website LXCat sim's solver, held to LoKI-MC within 1 % where both
apply). ML checked it against LoKI-MC for He + 1 % F2 (2026-09-27):
  - 6-15 Td: the two agree to 1 % in <eps>, W and k_att;
  - 3 Td: the two-term solution fails (attachment outruns energy
    relaxation), so He/F2 points below 6 Td are not modelled.

The files in ML/data/swarm/ carry their provenance in their headers.

Datasets (name: gas, what, role):
  kr_w_pack       Kr, W(E/N) 0.0008-3 Td, Pack 1962 via Dutton 1975   check
  kr_w_nakamura   Kr, W(E/N) 0.34-49 Td, Nakamura 1991 (figure)        check
  kr_alpha_kruithof, kr_alpha_heylen  Kr, alpha/N, via Dutton 1975     check
  f2_katt_nygaard He + 1 % F2, k_att = (eta/N_F2) x W, both measured
                  (Nygaard's own route to his rate constant), 6-15 Td  train
  f2_katt_tj300, f2_katt_tj500  N2 + F2 (1 Torr in 1 atm), k_att vs mean
                  energy, Trainor & Jacob 1979 (Avco 1980 figure)      train
  f2_katt_chen    N2 + F2, k_att vs mean energy, Chen 1977 (figure)    train
  f2_katt_sb      N2 or Ar + F2, Schneider & Brau 1978 (figure; the
                  buffer of each point is not given)                   held_out
  f2_katt_thermal thermal k_att: Sides 1976 (350, 600 K), McCorkle 1986
                  (298 K), which disagree by x5                        held_out
The N2 data give only the mean energy, which the authors computed with
their own N2 Boltzmann code. The model compares at equal mean energy from
this solver.

Training target: the F2 attachment cross section as two resonances. The
main one (0-2.5 eV) sets the laser's attachment, the 6 eV one sets
Nygaard's upper E/N. Each has a scale factor:
    sigma'(e) = sigma(e) [s_low + (s_high - s_low) / (1 + exp(-(e - 2.5)/0.2))]
For trace F2 the EEDF does not depend on the attachment cross section, so
k_att is a linear functional of it over EEDFs solved once. For 1 % F2 in He
this holds to within `linear_check`, which re-solves.
"""
from __future__ import annotations

import copy
import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import ROOT
from . import boltzmann  # noqa: F401  (puts the vendored solver on sys.path)
import boltz_twin as BT  # noqa: E402
import gases_twin as GT  # noqa: E402
import lx_lib as L  # noqa: E402

DATA = Path(__file__).resolve().parent.parent / "data" / "swarm"
INP = str(ROOT / "LoKI-MC" / "Code" / "Input")
GT.load_tables(INP)
GT.SETS.setdefault("F2_TRINITI", (["Fluorine2/F2_LXCat_TRINITI.txt"], "F2"))
E_SPLIT, E_WIDTH = 2.5, 0.2      # eV: the minimum between the resonances
X_F2_TJ = 1.0 / 760.0            # 1 Torr F2 in 1 atm N2 (Trainor & Jacob)
X_F2_TRACE = X_F2_TJ             # Chen, Schneider & Brau: trace F2 (the
                                 # EEDF of N2 is the same at 1e-3 or 1.3e-3)


@dataclass
class Swarm:
    name: str
    kind: str            # W | alphaN | katt_EN | katt_mean | katt_thermal
    buffer: str          # Kr | He | N2
    x_f2: float
    T: float
    x: np.ndarray        # E/N [Td], mean energy [eV] or gas temperature [K]
    y: np.ndarray        # W [cm/s], alpha/N [cm2] or k [cm3/s]
    sigma: np.ndarray    # fractional 1-sigma error per point
    role: str            # train | held_out | check
    source: str
    note: str = ""
    cache: dict = field(default_factory=dict, repr=False)


def _rows(fname):
    with open(DATA / fname) as fh:
        return list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))


def _interp_log(x, xp, fp):
    return np.exp(np.interp(x, xp, np.log(fp)))


def load() -> list:
    out = []
    r = [q for q in _rows("Kr_W_Pack1962_via_Dutton1975.csv")
         if q["T_gas_K"] == "300"]
    out.append(Swarm(
        "kr_w_pack", "W", "Kr", 0.0, 300.0,
        np.array([float(q["E_N_Td"]) for q in r]),
        np.array([float(q["W_cm_per_s"]) for q in r]), np.full(len(r), 0.05),
        "check", "Pack, Voshall & Phelps 1962, via Dutton 1975 Table 1.9(a)",
        "no error printed; 5 % assumed. Two W values repeat (Dutton's "
        "typos?)"))
    r = [q for q in _rows("Kr_W_NDL_Nakamura1991_JPSabstract_digitized.csv")
         if q["quantity"] == "W"]
    out.append(Swarm(
        "kr_w_nakamura", "W", "Kr", 0.0, 300.0,
        np.array([float(q["E_N_Td"]) for q in r]),
        np.array([float(q["value"]) for q in r]), np.full(len(r), 0.05),
        "check", "Nakamura 1991, JPS meeting abstract Fig. 1 (digitized)",
        "error not stated; 3 % reading + 4 % assumed"))
    r = _rows("Kr_alphaN_Kruithof1940_Heylen1971_via_Dutton1975.csv")
    for src, nm in (("Kruithof1940", "kr_alpha_kruithof"),
                    ("Heylen1971", "kr_alpha_heylen")):
        q = [a for a in r if a["source"] == src]
        out.append(Swarm(
            nm, "alphaN", "Kr", 0.0, 300.0,
            np.array([float(a["E_N_Td"]) for a in q]),
            np.array([float(a["alpha_over_N_cm2"]) for a in q]),
            np.full(len(q), 0.10), "check",
            f"{src[:-4]} {src[-4:]}, via Dutton 1975 Table 7.4(a)",
            "no error printed; 10 % assumed. Dutton: Heylen falls x2 below "
            "Kruithof near 30 Td"))
    # Nygaard: k_att = (eta/N_F2) W, W interpolated in the measured W(E/N)
    eta = _rows("F2_He_etaN_Nygaard1977_digitized.csv")
    wr = [q for q in _rows("F2_He_W_Nygaard1977_digitized.csv")
          if q["F2_percent"] == "1.0"]
    wx = np.array([float(q["E_N_Td"]) for q in wr])
    wy = np.array([float(q["W_cm_per_s"]) for q in wr])
    o = np.argsort(wx)
    wx, wy = wx[o], wy[o]
    en = np.array([float(q["E_N_Td"]) for q in eta])
    ev = np.array([float(q["eta_over_N_F2_cm2"]) for q in eta])
    sel = (en >= wx[0]) & (en <= wx[-1])
    # W(E/N) is noisy (7.57 and 7.97 Td cross); a straight line through
    # the six points stands in for it
    fitw = np.polyfit(wx, wy, 1)
    wmod = np.polyval(fitw, en[sel])
    out.append(Swarm(
        "f2_katt_nygaard", "katt_EN", "He", 0.01, 300.0, en[sel],
        ev[sel] * wmod, np.full(int(sel.sum()), math.hypot(0.12, 0.13)),
        "train", "Nygaard, Hunter, Fletcher & Foltyn 1977/78 (DTIC "
        "AD-A048764 Figs. 3 and 2, digitized)",
        "k = (eta/N_F2) x W as the authors did; W is a line through their "
        "six 1 % F2 points (errors 12 % and 13 % added in quadrature)"))
    av = _rows("F2_katt_vs_mean_energy_Avco1980_digitized.csv")
    for ds, T, nm, xf2, sg, role in (
            ("TrainorJacob1979", "300", "f2_katt_tj300", X_F2_TJ, 0.30,
             "train"),
            ("TrainorJacob1979", "500", "f2_katt_tj500", X_F2_TJ, 0.30,
             "train"),
            ("ChenEtAl1977", "300", "f2_katt_chen", X_F2_TRACE, 0.30,
             "train"),
            ("SchneiderBrau1978", "300", "f2_katt_sb", X_F2_TRACE, 0.10,
             "held_out")):
        q = [a for a in av if a["dataset"] == ds and a["T_gas_K"] == T]
        out.append(Swarm(
            nm, "katt_mean", "N2", xf2, float(T),
            np.array([float(a["mean_energy_eV"]) for a in q]),
            np.array([float(a["k_att_cm3_per_s"]) for a in q]),
            np.full(len(q), sg), role,
            f"{ds[:-4]} {ds[-4:]} ({T} K), via Avco DOE/DP/40103-4 App. B",
            "mean energy from the authors' N2 Boltzmann code, not measured"))
    th = _rows("F2_attachment_rate_constants_printed.csv")
    pick = [(350.0, 3.1e-9, 1.2 / 3.1), (600.0, 4.6e-9, 1.2 / 4.6),
            (298.0, 1.7e-8, 0.5)]
    out.append(Swarm(
        "f2_katt_thermal", "katt_thermal", "N2", X_F2_TRACE, 300.0,
        np.array([p[0] for p in pick]), np.array([p[1] for p in pick]),
        np.array([p[2] for p in pick]), "held_out",
        "Sides, Tiernan & Hanrahan 1976 (350, 600 K); McCorkle et al. 1986 "
        "(298 K), from Nygaard's report and the OSTI abstract",
        "Maxwellian at the stated temperature; McCorkle's F2 density was "
        "uncertain (50 % assumed)"))
    assert len(th) >= 6
    return out


def by_name(*names):
    d = {s.name: s for s in load()}
    return [d[n] for n in names]


# ------------------------------------------------------------ the forward model
_GAS: dict = {}


def gas(gid, T):
    key = (gid, round(T, 1))
    if key not in _GAS:
        _GAS[key] = GT.build(gid, T, INP)
    return _GAS[key]


def att_factor(e, s_low=1.0, s_high=1.0):
    return s_low + (s_high - s_low) / (1.0 + np.exp(-(e - E_SPLIT) / E_WIDTH))


def f2_modified(T, s_low=1.0, s_high=1.0, f2_set="F2"):
    """F2 with the attachment cross section rescaled (a copy)."""
    g = copy.deepcopy(gas(f2_set, T))
    for c in g["colls"]:
        if c["type"] == "Attachment":
            c["s"] = c["s"] * att_factor(c["e"], s_low, s_high)
    return g


def _att(g):
    return next(c for c in g["colls"] if c["type"] == "Attachment")


def _mix(ds, f2g):
    return [(gas(ds.buffer, ds.T), 1.0 - ds.x_f2), (f2g, ds.x_f2)]


_SOL: dict = {}


def _solve(mix, EN, T, key=None):
    """Two-term solution; key (hashable) caches it."""
    if key is not None and (key, round(EN, 9)) in _SOL:
        return _SOL[(key, round(EN, 9))]
    sol = BT.auto_solve(mix, EN, T, 400, 20.0)
    out = (sol, BT.swarm(mix, sol))
    if key is not None:
        _SOL[(key, round(EN, 9))] = out
    return out


def eedfs(ds, f2_set="F2"):
    """The EEDF of each point with the unmodified F2 set (cached)."""
    if f2_set in ds.cache:
        return ds.cache[f2_set]
    f2g = gas(f2_set, ds.T)
    out = []
    if ds.kind == "katt_EN":
        for en in ds.x:
            sol, sw = _solve(_mix(ds, f2g), en, ds.T)
            out.append(dict(sol=sol, mean=sw["mean"], EN=en))
    elif ds.kind == "katt_mean":
        from scipy.optimize import brentq
        key = (ds.buffer, ds.T, ds.x_f2, f2_set)

        def mean_at(lg):
            return _solve(_mix(ds, f2g), 10.0 ** lg, ds.T, key)[1]["mean"]
        grid = np.linspace(-1.0, 2.5, 15)          # log10 E/N [Td]
        means = np.array([mean_at(z) for z in grid])
        for m in ds.x:
            i = int(np.searchsorted(means, m))
            lg = brentq(lambda z: mean_at(z) - m, grid[i - 1], grid[i],
                        xtol=1e-4)
            sol, sw = _solve(_mix(ds, f2g), 10.0 ** lg, ds.T, key)
            out.append(dict(sol=sol, mean=sw["mean"], EN=10.0 ** lg))
    elif ds.kind == "katt_thermal":
        for T in ds.x:
            du = 0.0005
            ec = (np.arange(4000) + 0.5) * du
            kte = L.KB_EV * T
            F = np.exp(-ec / kte)
            F /= float(np.sum(np.sqrt(ec) * F) * du)
            out.append(dict(sol=dict(F=F, ec=ec, du=du), mean=1.5 * kte))
    else:
        raise ValueError(ds.kind)
    ds.cache[f2_set] = out
    return out


def predict(ds, s_low=1.0, s_high=1.0, f2_set="F2"):
    """Model values at every point of a dataset."""
    if ds.kind in ("W", "alphaN"):
        key = ("kr", ds.kind)
        if key not in ds.cache:
            vals = []
            for en in ds.x:
                sol, sw = _solve([(gas("Kr", ds.T), 1.0)], en, ds.T)
                vals.append(sw["w"] * 100.0 if ds.kind == "W"
                            else sw["alphaN"] * 1e4)
            ds.cache[key] = np.array(vals)
        return ds.cache[key]
    f2g = f2_modified(ds.T, s_low, s_high, f2_set)
    c = _att(f2g)
    return np.array([L.GAMMA * float(np.sum(
        p["sol"]["ec"] * L.interp(c, p["sol"]["ec"]) * p["sol"]["F"]))
        * p["sol"]["du"] * 1e6 for p in eedfs(ds, f2_set)])


def linear_check(ds, s_low, s_high, f2_set="F2"):
    """Re-solve the Boltzmann equation with the rescaled cross section:
    (exact k, linear k) per point."""
    f2g = f2_modified(ds.T, s_low, s_high, f2_set)
    lin = predict(ds, s_low, s_high, f2_set)
    exact = []
    for p in eedfs(ds, f2_set):
        sol, sw = _solve(_mix(ds, f2g), p["EN"], ds.T)
        exact.append(L.GAMMA * float(np.sum(
            sol["ec"] * L.interp(_att(f2g), sol["ec"]) * sol["F"]))
            * sol["du"] * 1e6)
    return np.array(exact), lin


def chi2(ds, s_low=1.0, s_high=1.0, f2_set="F2"):
    """chi2 per point in log space (errors are fractional)."""
    m = predict(ds, s_low, s_high, f2_set)
    r = np.log(m / ds.y) / ds.sigma
    return float(np.mean(r ** 2))


def fit(sets, f2_set="F2"):
    """Least squares for (s_low, s_high) in log space."""
    from scipy.optimize import least_squares

    def res(z):
        return np.concatenate([np.log(predict(d, *np.exp(z), f2_set) / d.y)
                               / d.sigma for d in sets])
    r = least_squares(res, np.zeros(2), bounds=(np.log(0.05), np.log(20.0)))
    J = r.jac
    cov = np.linalg.pinv(J.T @ J)
    return np.exp(r.x), np.sqrt(np.diag(cov)), r


# The F2 attachment train set. Leave-one-out (2026-09-27): Nygaard and
# Trainor & Jacob at 300 K agree with each other (chi2 per point 0.56, 0.86
# at the joint fit). Chen conflicts with Trainor & Jacob (x2 at 1 eV, the
# known group-to-group spread). Trainor & Jacob's 500 K points carry a
# +25 % gas-temperature rise that no 300 K cross section has.
TRAIN_F2 = ("f2_katt_nygaard", "f2_katt_tj300")


def fit_low(sets, f2_set="F2"):
    """One scale on the main attachment resonance (the 6 eV one is left as
    the set has it: it carries under 10 % of any dataset's rate and 0.3 %
    of the laser's). Returns (s_low, 1-sigma of ln s_low)."""
    from scipy.optimize import minimize_scalar

    def tot(z):
        return sum(float(np.sum((np.log(predict(d, math.exp(z), 1.0, f2_set)
                                        / d.y) / d.sigma) ** 2))
                   for d in sets)
    r = minimize_scalar(tot, bounds=(-3.0, 3.0), method="bounded")
    h = 1e-3
    curv = (tot(r.x + h) - 2 * tot(r.x) + tot(r.x - h)) / h ** 2
    return math.exp(r.x), math.sqrt(2.0 / curv)
