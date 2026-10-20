"""Forward models: what each dataset's experiment would measure, for a given
set of KrF(B,C) constants.

Kinetic kinds use LIGHT's closed-form B/C block (v0/two_level.py, held to
LIGHT's own run within 1e-5):
  b_only_decay  single exponential fitted to B(t) after a B-only start, down
                to 5 % (what Eden's deconvolution recovers). KrF2 at the
                stated pressure quenches at 1.4e-10 (Eden Fig. 2).
  slow_decay    the slow eigenvalue lam- (Gadomski eq. 5)
  eq_ratio      Yu & Setser eq. 18: I_r = 1.035 (tau_B/tau_C) K_eq(T)
  eq_ratio_T    I_r(T1)/I_r(T2) = K_eq(T1)/K_eq(T2)
  tau_B         the B radiative lifetime
Device kinds run LIGHT itself (basic network + ions + f2e) under the
prototype beam drive of LIGHT/docs/survey_2026-09-27/scripts/ebeam_screen.py.
The set's constants are imposed on LIGHT's rows, and its k_e is re-derived.
Under LIGHT_VERSION=v0_1 the ion and F2-electron rows are basic rows, B's
radiative rate is shared between B-X and B-A, and the drive's script is
loaded against v0_1's model (ml/__init__.py).
"""
from __future__ import annotations

import sys
import types

import numpy as np

from . import LIGHT_DIR, LIGHT_VERSION, light
from .data import TORR

L = light("two_level")
xsec = light("xsec")

KRF2_QUENCH = 1.4e-10          # Eden Fig. 2
QE_SI = 1.602176634e-19        # C (eV -> J)
_T_GRID = np.arange(0.0, 150e-9, 0.1e-9)
SIGMA_B = 2.6e-16              # cm2, Tellinghuisen 1976 via Kannari Table I
HV_EV = 5.00
# Trainor & Jacob's measured product and its re-reading
# (LIGHT/docs/krf_data/electron_quench.py)
TJ_PRODUCT, TJ_KR_TORR, TE_KHS, F_REC_KHS = 2.9e-7, 190.0, 1.58, 0.62


def b_only_sigma(c, kr, f2, krf2=0.0):
    Bn, _ = L.impulse(c, _T_GRID, B0=1.0, C0=0.0, kr=kr, f2=f2,
                      extra=KRF2_QUENCH * krf2)
    sel = Bn > 0.05
    return -np.polyfit(_T_GRID[sel], np.log(Bn[sel]), 1)[0]


def predict(ds, c) -> np.ndarray:
    """Model values for every point of a kinetic dataset (300 K data)."""
    c = c.with_(k_e=0.0)
    out = []
    for cond, _ in ds.points:
        if ds.kind == "b_only_decay":
            out.append(b_only_sigma(c, cond["kr_torr"] * TORR,
                                    cond["f2_torr"] * TORR,
                                    cond.get("krf2_torr", 0.0) * TORR))
        elif ds.kind == "slow_decay":
            out.append(L.eigen(c, kr=cond["kr_torr"] * TORR,
                               f2=cond["f2_torr"] * TORR)[1])
        elif ds.kind == "eq_ratio":
            out.append(1.035 * c.tau_B / c.tau_C * c.keq(cond["T"]))
        elif ds.kind == "eq_ratio_T":
            out.append(c.keq(cond["T1"]) / c.keq(cond["T2"]))
        elif ds.kind == "tau_B":
            out.append(c.tau_B)
        else:
            raise ValueError(f"{ds.name}: {ds.kind} is not a kinetic kind")
    return np.array(out)


def with_rederived_k_e(c):
    """k_e from Trainor & Jacob's product k_e tau / k_att = 2.9e-7 s, read
    with this set's own Stern-Volmer lifetime at 190 Torr. The lifetime does
    not depend on k_e's size, so a unit value stands in while computing it."""
    ka = float(xsec.rate_table("F2", "attachment")(TE_KHS))
    tau = L.sv_tau(c.with_(k_e=1.0, k_e_C=None), kr=TJ_KR_TORR * L.TORR_300K)
    return c.with_(k_e=TJ_PRODUCT * ka * F_REC_KHS / tau)


def impose(model, c):
    """Put a constant set onto a LIGHT model's compiled rows."""
    def pick(n):
        v = getattr(c, n)
        return getattr(c, n[:-2]) if v is None else v
    table = {"G4C": 1.0 / c.tau_C, "BC": c.k_BC,
             "QKR": c.k_Kr, "QKRC": pick("k_Kr_C"), "G5": c.k3_Kr,
             "G5C": pick("k3_Kr_C"), "QF2": c.k_F2, "QF2C": pick("k_F2_C"),
             "QE": c.k_e, "QEC": c.k_e}
    seen = set()
    for cc in model.compiled:
        rid = cc.row.rid
        if cc.row.rate[:2] == ("inv_tau", "tau_sp_s"):
            # B's radiative rows take their share of 1/tau_B (v0: G4 alone;
            # v0_1: B-X and B-A, G4 and G4A)
            frac = cc.row.rate[2] if len(cc.row.rate) > 2 else 1.0
            cc.rate = frac / c.tau_B
            seen.add(rid)
        elif rid in table:
            cc.rate = table[rid]
            seen.add(rid)
        elif rid == "CB":                 # ("balance"): (k_BC/K_eq, theta)
            cc.rate = (c.k_BC / c.K_eq, cc.rate[1])
            seen.add(rid)
    missing = set(table) | {"CB", "G4"}
    assert missing <= seen, f"rows not found: {missing - seen}"
    return model


TAHOMA_COND = {"pressure_torr": 665.0, "f2_torr": 2.8, "power_W_cm3": 3.3e5,
               "pulse_s": 450e-9}
WINDOW = (200e-9, 400e-9)        # s, the quasi-steady window scored


def _drive():
    """The survey's beam-drive script. For another LIGHT version the same
    file is loaded with its two v0 imports pointed at that version, so the
    drive is identical and only the model under it changes."""
    scripts = LIGHT_DIR / "docs" / "survey_2026-09-27" / "scripts"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    if LIGHT_VERSION == "v0":
        import ebeam_screen
        return ebeam_screen
    name = f"ebeam_screen_{LIGHT_VERSION}"
    if name not in sys.modules:
        path = scripts / "ebeam_screen.py"
        src = path.read_text()
        if src.count("from v0.") != 2:
            raise RuntimeError(f"{path.name}: expected two 'from v0.' "
                               "imports to point at " + LIGHT_VERSION)
        mod = types.ModuleType(name)
        mod.__file__ = str(path)
        exec(compile(src.replace("from v0.", f"from {LIGHT_VERSION}."),
                     str(path), "exec"), mod.__dict__)
        sys.modules[name] = mod
    return sys.modules[name]


# LIGHT's electron closure until 2026-09-28: every rate over the Maxwellian
# at Te, (3/2)Te per captured electron. The kinetic loop (ml/eedf.py) and
# the tail-table states (ml/tail.py) run on it, so they never depend on the
# closure they test.
MAXWELL = {"eedf": "maxwellian", "capture": "three_halves"}


def tahoma_params(cond=None, kr2=True, closure=None, afterglow_s=150e-9):
    """LIGHT's Params for Tahoma: basic network + ions + f2e + kr2 (the
    Kr2* route with its measured constants since 2026-09-28). A version
    that carries those rows in its basic network (v0_1 the ions and f2e,
    v0_2 the Kr2* route too) takes them from there, so kr2=False cannot
    remove v0_2's route (formation.build zeroes its formation instead).
    closure: Params fields of the electron closure (default LIGHT's own)."""
    Params = light("model").Params
    cond = cond or TAHOMA_COND
    # v0_1 carries the ion and F2-electron rows in its basic network
    basic = light("network").BASIC_GROUPS
    groups = (tuple(g for g in ("ions", "f2e") if g in basic)
              + (("kr2",) if kr2 and "kr2" in basic else ()))
    return Params(network="basic", groups=groups,
                  pressure_atm=cond["pressure_torr"] / 760.0,
                  x_f2=cond["f2_torr"] / cond["pressure_torr"],
                  power_W_cm3=cond["power_W_cm3"], pulse_s=cond["pulse_s"],
                  t_end_s=cond["pulse_s"] + afterglow_s, ne_mode="free",
                  ne_cm3=1e10, n_save=1201, **(closure or {}))


def pair_rate(p):
    """Ion pairs cm^-3 s^-1 of the beam drive (Spencer-Fano W_eff)."""
    return p.power_W_cm3 / QE_SI / _drive().W_EFF


def run_device(c, cond=None, modify=None, closure=None, afterglow_s=150e-9):
    """Tahoma under the prototype beam drive: g0 (%/cm), formation
    efficiency, the B share of the B/C manifold, Tg, n_e and Te over
    200-400 ns. `modify(model)` may change the electron side (cross sections
    or EEDF, ml/electron.py) after the B/C constants are imposed. closure:
    as tahoma_params (MAXWELL for LIGHT's closure until 2026-09-28)."""
    BeamModel, W_EFF = _drive().BeamModel, _drive().W_EFF
    QE = light("model").QE
    p = tahoma_params(cond, closure=closure, afterglow_s=afterglow_s)
    m = impose(BeamModel(p), c)
    if modify is not None:
        modify(m)
    r = m.run()
    t = r.t
    man = ("KrF(B)", "KrF(C)")
    formed = np.zeros_like(t)
    for cc in m.compiled:
        net = (sum(cc.row.products.count(s) for s in man)
               - sum(cc.row.reactants.count(s) for s in man))
        if net > 0:
            formed += net * r.events[cc.row.rid]
    win = (t >= WINDOW[0]) & (t <= WINDOW[1])
    rate = np.polyfit(t[win], formed[win], 1)[0]
    nb = float(np.mean(r.n["KrF(B)"][win]))
    nc = float(np.mean(r.n["KrF(C)"][win]))
    # the window's mean densities and event rates [cm^-3, cm^-3 s^-1],
    # what the electron side (ml/eedf.py) needs
    dens = {s: float(np.mean(r.n[s][win])) for s in r.n}
    events = {cc.row.rid: float(np.polyfit(t[win], r.events[cc.row.rid][win],
                                           1)[0]) for cc in m.compiled}
    return {"g0": SIGMA_B * nb * 100.0,
            "formation": rate * HV_EV * QE / p.power_W_cm3,
            "b_share": nb / (nb + nc) if nb + nc > 0 else float("nan"),
            "tg": float(np.mean(r.tg[win])),
            "ne": float(np.mean(r.ne[win])), "te": float(np.mean(r.te[win])),
            "n": dens, "events": events,
            "pairs": p.power_W_cm3 / QE / W_EFF}


def b_share_1atm(c, T_gas=300.0):
    """Steady [B]/([B]+[C]) at 1 atm, 0.5 % F2, KrF formed in B."""
    n = 760.0 * TORR * 300.0 / T_gas
    return L.steady_fraction_B(c, kr=0.995 * n, f2=0.005 * n, T=T_gas)
