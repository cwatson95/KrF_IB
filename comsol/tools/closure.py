"""The electron closure the COMSOL kit carries, and the LIGHT runs it is
built from and tested against.

LIGHT's kinetic EEDF (two-term Boltzmann with e-e collisions) is solved
again and again during a run, at the local state. COMSOL cannot call it, so
the kit carries it as lookup tables: every electron rate, captured energy,
the elastic power and the wave's absorption, tabulated in ln Te at several
states of LIGHT's own afterglow with the wave on, and read at the local
electron density between those states (linear in ln n_e). The EEDF along a
table is the one the wave alone holds at each mean energy: the field of the
wave scanned from nothing to 30 times the compendium's 3-MW/cm2 field, no DC
field. Below the no-field mean and above the strongest field's, a table
point is LIGHT's own extension (its end value carried by the Maxwellian
change, model._tab_end). The tables are frozen: nothing is solved during a
run, and the same tables give the 0-D runs here, the 1-D slab
(slab_reference.py) and the COMSOL model.

Two LIGHT models carry the two gases:
  KrF  LIGHT v0_6 (code/LIGHT/v0_5_reheat): Kr + 0.5 % F2 at 1 atm, the
       afterglow reheat study's twins (its reheat.params_of);
  Kr   LIGHT_Kr v0_5 (code/Kr_LIGHT/v0_5): pure Kr at 1 atm. LIGHT_Kr
       carries the EM wave only with the given power ('on'), so its
       afterglow runs are one 'on' segment with the given power 0 and the
       steady EEDF (kin_memory_s 0), as LIGHT v0_6 treats the wave alone.

Every run of the kit is an afterglow run: 1 us of afterglow from LIGHT's
kinetic state at the end of the pulse (300 kW/cm3 for 1 us), which both
twins share. The code here never writes into code/: LIGHT is imported
read-only, without byte-code files.
"""
from __future__ import annotations

import dataclasses
import importlib
import inspect
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

CODE = Path(__file__).resolve().parents[2] / "code"
GASES = {
    "KrF": {"root": CODE / "LIGHT", "pkg": "v0_5_reheat",
            "label": "LIGHT v0_6 (code/LIGHT/v0_5_reheat)",
            "bands": {"248 nm (KrF B-X)": "P248", "147 nm (Kr2* VUV)": "P147"}},
    "Kr": {"root": CODE / "Kr_LIGHT", "pkg": "v0_5",
           "label": "LIGHT_Kr v0_5 (code/Kr_LIGHT/v0_5)",
           "bands": {"147 nm (Kr2* VUV)": "P147"}},
}
# the compendium's wave: 3 MW/cm2, E0 = 47.6 kV/cm (rms 33.66 kV/cm), at the
# envelope notes' pure-Kr frequency, which LIGHT v0_6's reheat study takes
E_RMS_V_CM = 47.6e3 / math.sqrt(2.0)
F_BUILD_HZ = 6.2e12
PULSE_S = 1.0e-6
AFTERGLOW_S = 1.0e-6
# the wave's field factors the tables are built along, and their ln Te grid
FACTORS = np.concatenate([[1e-3], np.geomspace(0.12, 30.0, 60)])
LNTE = np.linspace(math.log(0.10), math.log(20.0), 96)
N_SLICES = 8
QE = 1.602176634e-19
ME = 9.1093837015e-31


def load(gas: str) -> dict:
    """The gas's LIGHT modules (model, network, species, report, boltzmann,
    and reheat for KrF), imported from code/ as the launchers do."""
    g = GASES[gas]
    pkg = g["pkg"]
    for name in [m for m in sys.modules if m == pkg or m.startswith(pkg + ".")]:
        del sys.modules[name]
    sys.path.insert(0, str(g["root"]))
    try:
        names = ["model", "network", "species", "report", "boltzmann"]
        if gas == "KrF":
            names.append("reheat")
        mods = {m: importlib.import_module(f"{pkg}.{m}") for m in names}
    finally:
        sys.path.remove(str(g["root"]))
    return mods


def _wave_arg(model_cls) -> bool:
    """LIGHT v0_6's rhs and _kinetic_state take the wave's flag apart from
    the given power's; LIGHT_Kr v0_5 ties them."""
    return "wave" in inspect.signature(model_cls.rhs).parameters


def reference_params(mods, gas, wave: bool):
    """The full twin's Params (pulse + afterglow) on the kinetic EEDF: the
    reheat study's reference point for KrF, LIGHT_Kr's for Kr. With wave,
    the compendium's field at F_BUILD_HZ over the afterglow (KrF only: the
    pure-Kr wave runs from the end-of-pulse state, afterglow_params)."""
    if gas == "KrF":
        spec = ("kinetic", "ref", "1" if wave else "off",
                "all" if wave else "", "")
        return mods["reheat"].params_of(spec)
    if wave:
        raise ValueError("LIGHT_Kr's wave runs from the end of the pulse")
    return mods["model"].Params(eedf="kinetic")


def afterglow_params(mods, p_full, e_v_cm: float, f_hz: float,
                     duration_s: float = AFTERGLOW_S):
    """The afterglow alone, starting at the end of the pulse: no given
    power, the wave (rms field e_v_cm, 0 = the off twin) over all of it."""
    if _wave_arg(mods["model"].Model):
        kw = dict(power_W_cm3=0.0, pulse_s=0.0, t_end_s=duration_s,
                  em_E_V_cm=e_v_cm, em_f_Hz=f_hz)
        if e_v_cm > 0.0:
            kw.update(em_t_on_s=0.0, em_t_off_s=math.inf)
        return p_full.replace(**kw)
    return p_full.replace(power_W_cm3=0.0, pulse_s=duration_s,
                          t_end_s=duration_s, em_E_V_cm=e_v_cm, em_f_Hz=f_hz,
                          kin_memory_s=0.0)


def state_vector(M, res, k) -> np.ndarray:
    """LIGHT's scaled state at output point k of a run (densities, U_e,
    Tg); its ledger and event counters start at 0."""
    y = np.zeros(M.size)
    for s, i in M.idx.items():
        y[i] = res.n[s][k] / M.N0
    y[M.i_ue] = res.ue[k] / M.N0
    y[M.i_tg] = res.tg[k]
    return y


def start_class(mods):
    """LIGHT's Model, started from a given state."""
    base = mods["model"].Model

    class Start(base):
        def __init__(self, p, y0=None):
            super().__init__(p)
            self._y0 = y0

        def initial_state(self):
            return super().initial_state() if self._y0 is None \
                else self._y0.copy()
    return Start


def kinetic_state(M, y, t):
    """The kinetic solver's State at y with the wave in it."""
    if _wave_arg(type(M)):
        s, te, ne = M._kinetic_state(y, t, on=False, wave=True)
    else:
        s, te, ne = M._kinetic_state(y, t, on=True)
    return dataclasses.replace(s, omega=2.0 * math.pi * M.p.em_f_Hz)


def hf_family(mods, M, y, t, factors=FACTORS):
    """The kinetic EEDF the wave alone holds at state y, for the wave's rms
    field factor x M.p.em_E_V_cm, factor by factor (each warm-started from
    the last). Returns [(Te, k{j}, cap{j}, pel, A, factor)], A the absorbed
    power per electron per (V/cm)^2 of rms field [eV/s]."""
    s = kinetic_state(M, y, t)
    e_ref = M.p.em_E_V_cm
    out, warm = [], None
    for f in factors:
        e = f * e_ref
        sE = dataclasses.replace(s, en_em_td=e * 100.0 / s.n_col
                                 / mods["boltzmann"].TD)
        sol = M.kin.eedf(sE, None, warm=warm)
        warm = sol
        k, cap, pel, pem = M._kin_quantities(sE, sol)
        out.append((sol.te, k, cap, pel, pem / e ** 2, f))
    return out


def table_on_grid(mods, M, y, family, lnte=LNTE) -> dict:
    """One state's table on the kit's ln Te grid: LIGHT's interpolation
    (linear in ln Te; ln k for rates) inside the family's range, LIGHT's
    own extension outside it (model._tab_end: the end value carried by the
    Maxwellian change), and the Maxwellian values at the grid's ends, which
    carry the table beyond the grid as LIGHT would."""
    kb = mods["model"].KB_EV
    pts = {}
    for te, k, cap, pel, a, f in family:
        pts[round(math.log(te), 9)] = (k, cap, pel, a)
    lt_f, rows = [], []
    for key in sorted(pts):
        if lt_f and key - lt_f[-1] < 1e-4:
            continue
        lt_f.append(key)
        rows.append(pts[key])
    lt_f = np.array(lt_f)
    n = np.maximum(y[:M.ns], 0.0) * M.N0
    ne = M.electron_density(n)
    tg = y[M.i_tg]
    byj = {c.j: c for c in M.compiled}
    e2 = M.p.em_E_V_cm ** 2
    raw = {"lnk": {j: np.log(np.maximum([r[0][j] for r in rows], 1e-300))
                   for j in rows[0][0]},
           "cap": {j: np.array([r[1][j] for r in rows]) for j in rows[0][1]},
           "pel": np.array([r[2] for r in rows]),
           "A": np.array([r[3] for r in rows])}
    tab = {"lnte": lnte.copy(),
           "lnk": {j: np.interp(lnte, lt_f, v) for j, v in raw["lnk"].items()},
           "cap": {j: np.interp(lnte, lt_f, v) for j, v in raw["cap"].items()},
           "pel": np.interp(lnte, lt_f, raw["pel"]),
           "A": np.interp(lnte, lt_f, raw["A"])}
    for i, lt in enumerate(lnte):
        if lt_f[0] <= lt <= lt_f[-1]:
            continue
        end = 0 if lt < lt_f[0] else -1
        te, te_e = math.exp(lt), math.exp(lt_f[end])
        for j in tab["lnk"]:
            km = max(M._k(byj[j], te, tg, n, ne), 1e-300)
            kme = max(M._k(byj[j], te_e, tg, n, ne), 1e-300)
            tab["lnk"][j][i] = raw["lnk"][j][end] + math.log(km / kme)
        for j in tab["cap"]:
            if byj[j].cap is not None:
                tab["cap"][j][i] = (raw["cap"][j][end] * byj[j].cap(te)
                                    / byj[j].cap(te_e))
        tab["pel"][i] = (raw["pel"][end] + M._pel_maxwell(te, tg * kb)
                         - M._pel_maxwell(te_e, tg * kb))
        tab["A"][i] = raw["A"][end] + (M._em_maxwell(te, n, ne)
                                       - M._em_maxwell(te_e, n, ne)) / e2
    ends = (math.exp(lnte[0]), math.exp(lnte[-1]))
    tab["lnkM"] = {j: np.log(np.maximum([M._k(byj[j], e, tg, n, ne)
                                         for e in ends], 1e-300))
                   for j in tab["lnk"]}
    tab["capM"] = {j: np.array([byj[j].cap(e) for e in ends])
                   for j in tab["cap"] if byj[j].cap is not None}
    tab["pelM"] = np.array([M._pel_maxwell(e, tg * kb) for e in ends])
    tab["AM"] = np.array([M._em_maxwell(e, n, ne) / e2 for e in ends])
    tab["ln_ne"] = math.log(ne)
    tab["family_lnte"] = (float(lt_f[0]), float(lt_f[-1]))
    tab["state"] = {"ne_cm3": float(ne), "tg_K": float(tg),
                    "n_ion_cm3": float(n[M.cations].sum())}
    return tab


def blend_class(mods):
    """LIGHT's Model with the frozen tables, read at the local n_e (linear
    in ln n_e between the states, held beyond the first and the last), and
    never renewed. The wave's power per electron is A(Te) E^2 at LIGHT's
    field, or, in a slab cell, the cell's own absorbed power (p_em_cell)."""
    base = start_class(mods)

    class Blend(base):
        def __init__(self, p, tables, y0=None):
            self.tables = sorted(tables, key=lambda t: t["ln_ne"])
            self.ln_ne = np.array([t["ln_ne"] for t in self.tables])
            self.p_em_cell = None
            super().__init__(p, y0=y0)

        def blend(self, ne):
            x = math.log(max(ne, 1.0))
            tabs, xs = self.tables, self.ln_ne
            if x <= xs[0]:
                a = b = tabs[0]
                w = 0.0
            elif x >= xs[-1]:
                a = b = tabs[-1]
                w = 0.0
            else:
                i = int(np.searchsorted(xs, x)) - 1
                w = (x - xs[i]) / (xs[i + 1] - xs[i])
                a, b = tabs[i], tabs[i + 1]

            def mix(u, v):
                return (1.0 - w) * u + w * v
            e = self.p.em_E_V_cm
            out = {"lnte": a["lnte"],
                   "lnk": {j: mix(a["lnk"][j], b["lnk"][j]) for j in a["lnk"]},
                   "cap": {j: mix(a["cap"][j], b["cap"][j]) for j in a["cap"]},
                   "pel": mix(a["pel"], b["pel"]), "A": mix(a["A"], b["A"]),
                   "lnkM": {j: mix(a["lnkM"][j], b["lnkM"][j])
                            for j in a["lnkM"]},
                   "capM": {j: mix(a["capM"][j], b["capM"][j])
                            for j in a["capM"]},
                   "pelM": mix(a["pelM"], b["pelM"]),
                   "AM": mix(a["AM"], b["AM"])}
            out["pem"] = out["A"] * e * e
            out["pemM"] = out["AM"] * e * e
            return out

        def _set_table(self, y):
            n = np.maximum(y[:self.ns], 0.0) * self.N0
            self._tab = self.blend(self.electron_density(n))

        def rhs(self, t, y, *args):
            self._set_table(y)
            return super().rhs(t, y, *args)

        def _kinetic_update(self, t, y, *args):
            self._set_table(y)
            self._kin_ref = self._kin_observed(y)

        def _drift_events(self):
            return []

        def em_power(self, te, n, ne):
            if self.p_em_cell is not None:
                return self.p_em_cell
            return super().em_power(te, n, ne)
    return Blend


def nu_eff(A, omega_build=2.0 * math.pi * F_BUILD_HZ):
    """The effective collision frequency [1/s] that gives the table's
    absorbed power per electron in the Drude form, A E_rms^2 = (e^2 E_rms^2
    / m) nu / (nu^2 + w^2) at the build frequency (the smaller root; nu <<
    w at THz). The tables are built at w >> nu, where the EEDF the wave
    holds depends on E/w alone, so nu carries over to other frequencies
    while nu stays << w."""
    g = np.asarray(A, float) * ME / (QE * 1e4)
    disc = np.maximum(1.0 - 4.0 * g * g * omega_build ** 2, 0.0)
    return np.where(g > 0.0, (1.0 - np.sqrt(disc))
                    / (2.0 * np.maximum(g, 1e-300)), 0.0)
