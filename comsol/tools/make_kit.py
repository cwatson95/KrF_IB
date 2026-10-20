"""Generate the COMSOL kit: the afterglow-reheat slab for Kr/F2 (KrF) and
pure Kr (Kr), from LIGHT's own networks and kinetic EEDF.

    python3 tools/make_kit.py runs      # LIGHT's kinetic reference runs (cached)
    python3 tools/make_kit.py tables    # the closure tables (closure.py, cached)
    python3 tools/make_kit.py export    # the COMSOL files, both gases
    python3 tools/make_kit.py all       # the three in turn

(from KrF_IB/comsol; --gas KrF or --gas Kr for one gas). The reference
runs take 10 to 20 minutes each on the kinetic EEDF and are cached in
tools/cache/ (delete a file there to run it again).

What it writes, per gas (KrF/ or Kr/):
  kinetics_<gas>.txt   the electron closure as ONE interpolation function of
                       three arguments (COMSOL's Grid format): ln Te [ln eV],
                       ln n_e [ln cm^-3], and the quantity's index
  heavy_<gas>.txt      (Kr only) the heavy-species rates that depend on Tg
                       other than as a power law, as ln k in (ln Tg, index)
  params_<gas>.txt     COMSOL parameters (Load from File)
  variables_<gas>_plasma.txt   COMSOL variables for the Plasma Module model
  variables_<gas>_ode.txt      the same for the Domain-ODE model
  chemistry_<gas>.txt  species and reactions for Plasma Chemistry Import
  initial_<gas>.txt    the initial values (the end of the pulse), listed
  reference_<gas>.csv  LIGHT's afterglow, 0-D, wave off and on (to compare)
  kit_<gas>.json       everything above as data (check_kit.py and
                       slab_reference.py read it)
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np

import closure as C

HERE = Path(__file__).resolve().parent
KIT = HERE.parent
CACHE = HERE / "cache"
AMU = {"Kr": 83.798, "F": 18.998403}
# COMSOL names (Plasma Module: an ion's tag is name_1p / name_1m)
NAMES = {
    "KrF": {"Kr": "Kr", "F2": "F2", "F": "F", "Kr*": "Krs", "Kr**": "Krp",
            "Kr+": "Kr+", "KrF(B)": "KrFB", "Kr2F*": "Kr2Fs",
            "KrF(C)": "KrFC", "Kr2+": "Kr2+", "F-": "F-", "Kr2*": "Kr2s"},
    "Kr": {"Kr": "Kr", "Kr*(1s5)": "Kr1s5", "Kr*(1s4)": "Kr1s4",
           "Kr*(1s3)": "Kr1s3", "Kr*(1s2)": "Kr1s2", "Kr**": "Krp",
           "Kr+": "Kr+", "Kr2+": "Kr2+", "Kr2*": "Kr2s",
           "Kr2*(0u+)": "Kr2s0u"},
}
# Lennard-Jones parameters for the user-defined species [Angstrom, K]: they
# set the heavy species' diffusion only, which moves nothing in 1 us at 1 atm
LJ = {"Kr": (3.655, 178.9), "F2": (3.357, 112.6), "F": (2.968, 112.6),
      "KrF": (3.9, 200.0), "Kr2F": (4.6, 250.0), "Kr2": (4.4, 220.0)}


def tag(name: str) -> str:
    """COMSOL's tag of a species name: Kr+ -> Kr_1p, F- -> F_1m."""
    if name.endswith("+"):
        return name[:-1] + "_1p"
    if name.endswith("-"):
        return name[:-1] + "_1m"
    return name


# ---------------------------------------------------------------------------
# LIGHT's reference runs
# ---------------------------------------------------------------------------
def record(res, kin_log=None) -> dict:
    out = {"t": np.asarray(res.t), "n": {k: np.asarray(v) for k, v in res.n.items()},
           "ne": np.asarray(res.ne), "ue": np.asarray(res.ue),
           "te": np.asarray(res.te), "tg": np.asarray(res.tg),
           "ledger": {k: np.asarray(v) for k, v in res.ledger.items()},
           "events": {k: np.asarray(v) for k, v in res.events.items()}}
    return out


class Rec:
    def __init__(self, d):
        self.__dict__.update(d)


def cached(name: str, fn):
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"{name}.pkl"
    if path.is_file():
        return pickle.loads(path.read_bytes())
    t0 = time.time()
    print(f"  running {name} ...", flush=True)
    out = fn()
    out["seconds"] = time.time() - t0
    path.write_bytes(pickle.dumps(out))
    print(f"  {name}: {out['seconds']:.0f} s", flush=True)
    return out


def at(rec, t: float) -> int:
    return int(np.argmin(abs(np.asarray(rec["t"]) - t)))


def runs(gas: str, mods=None) -> dict:
    """The kinetic reference runs: the pump and the afterglow without the
    wave ('full'), and the afterglow twins from the end of the pulse with
    LIGHT's renewed kinetic EEDF ('ag_off', 'ag_on'; for KrF the reheat
    study's own 'on' twin instead, 'full_on')."""
    mods = mods or C.load(gas)
    model = mods["model"]
    out = {}
    out["full"] = cached(f"{gas}_full", lambda: record(
        model.run(C.reference_params(mods, gas, wave=False))))
    full = Rec(out["full"])
    p_full = C.reference_params(mods, gas, wave=False)
    k1 = at(out["full"], C.PULSE_S)
    Start = C.start_class(mods)
    M = Start(p_full)
    y1 = C.state_vector(M, full, k1)
    if gas == "KrF":
        out["full_on"] = cached(f"{gas}_full_on", lambda: record(
            model.run(C.reference_params(mods, gas, wave=True))))
    for name, e in (("ag_off", 0.0), ("ag_on", C.E_RMS_V_CM)):
        def go(e=e):
            p = C.afterglow_params(mods, p_full, e, C.F_BUILD_HZ)
            return record(Start(p, y0=y1).run())
        out[name] = cached(f"{gas}_{name}", go)
    out["y1"] = y1
    return out


# ---------------------------------------------------------------------------
# the closure tables
# ---------------------------------------------------------------------------
def tables(gas: str, mods=None, R=None) -> list:
    """N_SLICES tables at states of the 'on' afterglow whose n_e is
    log-spaced from 2 ns after the pulse to 0.8 us after it."""
    mods = mods or C.load(gas)
    R = R or runs(gas, mods)

    def build():
        on = R["full_on"] if gas == "KrF" else R["ag_on"]
        t0 = C.PULSE_S if gas == "KrF" else 0.0
        t = np.asarray(on["t"])
        ne = np.asarray(on["ne"])
        m = (t >= t0 + 2e-9) & (t <= t0 + 0.8e-6)
        targets = np.geomspace(ne[m][0], ne[m][-1], C.N_SLICES)
        picks = []
        for v in targets:
            k = np.flatnonzero(m)[int(np.argmin(abs(np.log(ne[m] / v))))]
            if k not in picks:
                picks.append(int(k))
        p_full = C.reference_params(mods, gas, wave=False)
        p = C.afterglow_params(mods, p_full, C.E_RMS_V_CM, C.F_BUILD_HZ)
        M = C.start_class(mods)(p)
        tabs = []
        for k in picks:
            y = C.state_vector(M, Rec(on), k)
            fam = C.hf_family(mods, M, y, t[k])
            tab = C.table_on_grid(mods, M, y, fam)
            tab["t_us"] = float(t[k] * 1e6)
            tabs.append(tab)
            print(f"    {gas} table at t = {t[k]*1e6:.4f} us, n_e "
                  f"{tab['state']['ne_cm3']:.3g} cm-3, family Te "
                  f"{math.exp(tab['family_lnte'][0]):.3f}-"
                  f"{math.exp(tab['family_lnte'][1]):.3f} eV", flush=True)
        return {"tables": tabs}
    return cached(f"{gas}_tables", build)["tables"]


# ---------------------------------------------------------------------------
# the network as COMSOL reactions
# ---------------------------------------------------------------------------
def species_rows(mods, gas, M) -> list:
    sp = mods["species"]
    out = []
    for s in M.species:
        d = sp.SPECIES[s]
        atoms = dict(d.atoms)
        mass = sum(AMU[a] * k for a, k in atoms.items())
        nkr, nf = atoms.get("Kr", 0), atoms.get("F", 0)
        lj = (LJ["Kr"] if (nkr, nf) == (1, 0) else LJ["F2"] if (nkr, nf) == (0, 2)
              else LJ["F"] if (nkr, nf) == (0, 1) else LJ["KrF"] if (nkr, nf) == (1, 1)
              else LJ["Kr2F"] if (nkr, nf) == (2, 1) else LJ["Kr2"])
        name = NAMES[gas][s]
        out.append({"light": s, "name": name, "tag": tag(name),
                    "charge": d.charge, "mass_kg_mol": mass * 1e-3,
                    "cv": d.cv, "e_form_eV": d.e_form_eV, "lj": lj})
    return out


def formula(names, side) -> str:
    return "+".join(names.get(s, s) for s in side)


def group_rows(M) -> tuple:
    """Active rows (factor > 0) lumped by (reactants, products): one COMSOL
    reaction each. The elastic row is not a reaction (its energy is the
    elastic power, pel)."""
    el, hv = OrderedDict(), OrderedDict()
    for c in M.compiled:
        if c.factor == 0.0:
            continue
        if c.kind == "elastic":
            continue
        key = (tuple(sorted(c.row.reactants)), tuple(sorted(c.row.products)))
        (el if c.row.n_electrons_in else hv).setdefault(key, []).append(c)
    return list(el.values()), list(hv.values())


def gname(cs) -> str:
    """A group's name: its row id, or the first id + '_grp' when lumped."""
    if len(cs) == 1:
        return cs[0].row.rid.replace("-", "_")
    return cs[0].row.rid.replace("-", "_") + f"_x{len(cs)}"


def heavy_law(mods, M, c) -> dict:
    """A heavy row's rate law, as COMSOL can carry it [LIGHT's cm units]."""
    kind = c.kind
    tgp = c.row.tg_pow
    if M.p.__dict__.get("trimer_tg_pow") is not None and c.row.rid in ("G5", "G5C"):
        tgp = M.p.trimer_tg_pow
    tgfac = getattr(M, "_tg_fac", {})
    if kind in ("const", "inv_tau", "param") and c.j not in tgfac:
        return {"type": "power", "k300": float(c.rate) * c.factor, "n": float(tgp)}
    if kind == "balance":
        return {"type": "balance", "a": float(c.rate[0]) * c.factor,
                "theta": float(c.rate[1]), "n": float(tgp)}
    if kind == "ionion":
        k3, kl, n3 = c.rate
        return {"type": "ionion", "k3": float(k3), "kL": float(kl),
                "n3": float(n3), "f": float(c.factor), "n": float(tgp)}
    return {"type": "tgtab"}


def molecularity(c) -> int:
    return len(c.row.reactants)


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------
def build_kit(gas: str) -> dict:
    mods = C.load(gas)
    R = runs(gas, mods)
    tabs = tables(gas, mods, R)
    p_full = C.reference_params(mods, gas, wave=False)
    p_on = C.afterglow_params(mods, p_full, C.E_RMS_V_CM, C.F_BUILD_HZ)
    Blend = C.blend_class(mods)
    M = Blend(p_on, tabs, y0=R["y1"])
    spec = species_rows(mods, gas, M)
    names = {s["light"]: s["name"] for s in spec}
    egroups, hgroups = group_rows(M)
    byrid = {c.row.rid: c for c in M.compiled}
    elastic = [c for c in M.compiled if c.kind == "elastic" and c.factor != 0]
    if [c.row.rate[1] for c in elastic] != ["Kr"]:
        raise ValueError("the kit books Kr's elastic power only")
    lnte = tabs[0]["lnte"]
    tabs_sorted = sorted(tabs, key=lambda t: t["ln_ne"])
    lnne = np.array([t["ln_ne"] for t in tabs_sorted])
    # ---- the electron quantities on the grid ------------------------------
    quants = []             # (key, description, values[slice, lnte])

    def add(key, desc, fn):
        vals = np.empty((len(tabs_sorted), len(lnte)))
        for s, tab in enumerate(tabs_sorted):
            M._tab = M.blend(math.exp(tab["ln_ne"]))
            st = tab["state"]
            for i, lt in enumerate(lnte):
                vals[s, i] = fn(M, math.exp(lt), st, lt)
        quants.append((key, desc, vals))
        return len(quants)                       # COMSOL's index (1-based)

    n_ref = np.maximum(R["y1"][:M.ns], 0.0) * M.N0

    def rate(c, te, st):
        return M._rate(c, te, st["tg_K"], n_ref, st["ne_cm3"])

    eidx = []
    for cs in egroups:
        name = gname(cs)
        i_k = add(f"lnk_{name}", "ln of the rate coefficient [cm3/s]",
                  lambda M, te, st, lt, cs=cs: math.log(max(sum(
                      rate(c, te, st) for c in cs), 1e-300)))
        de = [c.route["to_e"] - c.route["from_e"] for c in cs]
        heat = [c.route["heat"] for c in cs]
        rec = {"name": name, "rows": [c.row.rid for c in cs],
               "reactants": list(cs[0].row.reactants),
               "products": list(cs[0].row.products), "i_k": i_k}
        if len(set(np.round(de, 9))) > 1:
            rec["i_de"] = add(f"de_{name}", "electron energy change per event [eV]",
                              lambda M, te, st, lt, cs=cs: sum(
                                  rate(c, te, st) * (c.route["to_e"] - c.route["from_e"])
                                  for c in cs) / max(sum(rate(c, te, st) for c in cs), 1e-300))
        else:
            rec["de"] = float(de[0])
        if len(set(np.round(heat, 9))) > 1:
            rec["i_heat"] = add(f"heat_{name}", "gas heat per event [eV]",
                                lambda M, te, st, lt, cs=cs: sum(
                                    rate(c, te, st) * c.route["heat"] for c in cs)
                                / max(sum(rate(c, te, st) for c in cs), 1e-300))
        else:
            rec["heat"] = float(heat[0])
        caps = [c for c in cs if c.route["capture"]]
        if caps:
            if len(cs) != 1:
                raise ValueError(f"{name}: a lumped capture group")
            c = caps[0]

            def capfn(M, te, st, lt, c=c):
                tab = M._tab
                if c.j in tab["cap"]:
                    return float(np.interp(lt, tab["lnte"], tab["cap"][c.j]))
                return 1.5 * te if c.cap is None else float(c.cap(te))
            rec["i_cap"] = add(f"cap_{name}", "energy of the captured electron [eV]", capfn)
            rec["capture"] = int(c.route["capture"])
        rads = {}
        for c in cs:
            if c.route["rad"]:
                rads.setdefault(c.row.band, []).append(c)
        if rads:
            raise ValueError(f"{name}: a radiative electron row")
        eidx.append(rec)
    i_pel = add("pel", "Kr elastic power per electron per Kr atom [eV cm3/s]",
                lambda M, te, st, lt: float(np.interp(lt, M._tab["lnte"], M._tab["pel"])))
    i_nu = add("lnnu", "ln of the effective collision frequency for the wave [1/s]",
               lambda M, te, st, lt: math.log(max(float(C.nu_eff(
                   np.interp(lt, M._tab["lnte"], M._tab["A"]))), 1.0)))
    # ---- heavy reactions -------------------------------------------------
    hidx, tg_tab = [], []
    lntg = np.linspace(math.log(250.0), math.log(5000.0), 80)
    for cs in hgroups:
        name = gname(cs)
        laws = [heavy_law(mods, M, c) for c in cs]
        rec = {"name": name, "rows": [c.row.rid for c in cs],
               "reactants": list(cs[0].row.reactants),
               "products": list(cs[0].row.products), "order": molecularity(cs[0])}
        if len(cs) == 1 and laws[0]["type"] != "tgtab":
            rec["law"] = laws[0]
        elif all(l["type"] == "power" and l["n"] == laws[0]["n"] for l in laws):
            rec["law"] = {"type": "power", "k300": sum(l["k300"] for l in laws),
                          "n": laws[0]["n"]}
        else:
            vals = []
            for lt in lntg:
                tg = math.exp(lt)
                k = sum(M._k(c, 1.0, tg, n_ref, float(M.electron_density(n_ref))) for c in cs)
                # the law must not move with n or n_e
                k2 = sum(M._k(c, 1.0, tg, n_ref * 1.3, 3.0 * float(M.electron_density(n_ref))) for c in cs)
                if not math.isclose(k, k2, rel_tol=1e-9, abs_tol=1e-300):
                    raise ValueError(f"{name}: depends on more than Tg")
                vals.append(math.log(max(k, 1e-300)))
            tg_tab.append(vals)
            rec["law"] = {"type": "tgtab", "index": len(tg_tab)}
        de = [c.route["to_e"] - c.route["from_e"] for c in cs]
        heat = [c.route["heat"] for c in cs]
        ks = [M._k(c, 1.0, 300.0, n_ref, 1e12) for c in cs]
        ksum = sum(ks) or 1.0
        rec["de"] = float(sum(k * d for k, d in zip(ks, de)) / ksum) if len(cs) > 1 else float(de[0])
        rec["heat"] = float(sum(k * h for k, h in zip(ks, heat)) / ksum) if len(cs) > 1 else float(heat[0])
        if len(cs) > 1 and (len(set(np.round(de, 9))) > 1 or len(set(np.round(heat, 9))) > 1):
            if any(l["type"] == "tgtab" for l in laws):
                raise ValueError(f"{name}: lumped rows of different energies and laws")
        rads = {}
        for c, k in zip(cs, ks):
            if c.route["rad"]:
                rads[c.row.band] = rads.get(c.row.band, 0.0) + k * c.route["rad"] / ksum
        rec["rad"] = rads
        rec["capture"] = int(max(c.route["capture"] for c in cs))
        if rec["capture"]:
            raise ValueError(f"{name}: a heavy capture row")
        hidx.append(rec)
    # ---- the initial state: the end of the pulse ---------------------------
    full = R["full"]
    k1 = at(full, C.PULSE_S)
    init = {s["light"]: float(full["n"][s["light"]][k1]) for s in spec}
    ne1 = float(full["ne"][k1])
    te1 = float(full["te"][k1])
    tg1 = float(full["tg"][k1])
    # what the pulse left behind (energies per cm3), common to both twins
    rows = {c.row.rid: c for c in M.compiled}
    pump_bands = {}
    for rid, c in rows.items():
        if c.route["rad"]:
            pump_bands[c.row.band] = pump_bands.get(c.row.band, 0.0) + \
                float(full["events"][rid][k1]) * c.route["rad"] * C.QE
    kit = {
        "gas": gas, "model": C.GASES[gas]["label"],
        "point": {"pressure_atm": M.p.pressure_atm,
                  "x_f2": getattr(M.p, "x_f2", 0.0),
                  "power_W_cm3": p_full.power_W_cm3, "pulse_s": p_full.pulse_s,
                  "tg0_K": p_full.tg0_K},
        "N0_cm3": float(M.N0), "species": spec,
        "initial": {"t_s": float(full["t"][k1]), "n_cm3": init, "ne_cm3": ne1,
                    "te_eV": te1, "tg_K": tg1},
        "pump_phase": {"E_given_J_cm3": float(full["ledger"]["E_dep"][k1]) * C.QE,
                       "bands_J_cm3": pump_bands},
        "grid": {"lnte": lnte.tolist(), "lnne": lnne.tolist(),
                 "quantities": [{"key": k, "desc": d, "values": v.tolist()} for k, d, v in quants],
                 "i_pel": i_pel, "i_nu": i_nu,
                 "states": [{"t_us": t["t_us"], **t["state"]} for t in tabs_sorted]},
        "heavy_grid": {"lntg": lntg.tolist(), "values": tg_tab},
        "electron_reactions": eidx, "heavy_reactions": hidx,
        "f_build_Hz": C.F_BUILD_HZ, "E_build_V_cm": C.E_RMS_V_CM,
    }
    return kit


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=("runs", "tables", "export", "all"))
    ap.add_argument("--gas", choices=("KrF", "Kr"), default=None)
    a = ap.parse_args()
    gases = [a.gas] if a.gas else ["KrF", "Kr"]
    for gas in gases:
        if a.step in ("runs", "all"):
            runs(gas)
        if a.step in ("tables", "all"):
            tables(gas)
        if a.step in ("export", "all"):
            import export
            export.write_all(build_kit(gas))
