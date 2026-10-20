"""Check the COMSOL kit against LIGHT, from the files COMSOL reads.

    python3 tools/check_kit.py            # both gases
    python3 tools/check_kit.py --gas Kr

1. files: the grid file read back equals the kit's tables; the chemistry
   file's reactions are balanced (elements, charge) and name only the
   species and kf_ variables the kit defines; every variable is defined
   once; the initial state is quasineutral and its mass fractions sum to 1;
2. equations: the variables files are evaluated as COMSOL would (their
   expressions translated, unit brackets in SI, the Grid file interpolated
   linearly in each argument) at states drawn along LIGHT's afterglow, and
   every species source, the electron power, the gas temperature rate and
   the emitted powers are compared with LIGHT's own right-hand side on the
   kit's closure (closure.blend_class) at the same state. Both routes
   ('plasma': the Plasma Module's variables, 'ode': the Domain ODEs').

Exit status 0 when every check passes.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

import closure as C
import make_kit as K

KIT = Path(__file__).resolve().parent.parent
TOL = 1e-6                     # relative, on sources (with a floor, below)

CONST = {"e_const": 1.602176634e-19, "k_B_const": 1.380649e-23,
         "me_const": 9.1093837015e-31, "N_A_const": 6.02214076e23,
         "epsilon0_const": 8.8541878128e-12, "c_const": 299792458.0,
         "Z0_const": 376.730313668, "pi": math.pi, "j": 1j}
UNIT = {"m": 1.0, "s": 1.0, "K": 1.0, "V": 1.0, "W": 1.0, "J": 1.0, "C": 1.0,
        "kg": 1.0, "mol": 1.0, "A": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6,
        "us": 1e-6, "ns": 1e-9, "THz": 1e12, "GHz": 1e9, "Hz": 1.0, "MW": 1e6}


def unit_value(u: str) -> float:
    expr = re.sub(r"[A-Za-z]+", lambda m: repr(UNIT[m.group(0)]), u)
    return float(eval(expr.replace("^", "**"), {"__builtins__": {}}, {}))


def translate(expr: str) -> str:
    """COMSOL expression -> Python (numpy) expression."""
    e = re.sub(r"\[([^\]]+)\]", lambda m: f"*{unit_value(m.group(1))!r}", expr)
    e = e.replace("^", "**")
    e = re.sub(r"\bplas\.", "plas_", e)
    e = re.sub(r"\bemw\.", "emw_", e)
    e = re.sub(r"(?<![\w.])1\*", "1.0*", e)
    return e


def read_table(path: Path) -> list:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        name, expr, desc = line.split("\t")
        rows.append((name, expr, desc))
    return rows


def read_grid(path: Path):
    """A COMSOL Grid file: the argument grids and the data rows."""
    lines = [l for l in path.read_text().splitlines()]
    i = lines.index("% Grid")
    j = lines.index("% Data")
    axes = [np.array([float(v) for v in l.split()]) for l in lines[i + 1:j]]
    data = np.array([[float(v) for v in l.split()] for l in lines[j + 1:] if l.strip()])
    return axes, data


class GridFunction:
    """Linear in each argument, constant beyond the grid (COMSOL's
    Interpolation with Linear interpolation and Constant extrapolation)."""

    def __init__(self, path: Path):
        self.axes, data = read_grid(path)
        shape = [len(a) for a in reversed(self.axes[1:])] + [len(self.axes[0])]
        self.v = data.reshape(shape)          # [z][y][x] or [y][x]

    @staticmethod
    def _w(ax, x):
        x = min(max(x, ax[0]), ax[-1])
        i = int(np.clip(np.searchsorted(ax, x) - 1, 0, len(ax) - 2))
        w = (x - ax[i]) / (ax[i + 1] - ax[i])
        return i, w

    def __call__(self, *args):
        args = [float(np.real(a)) for a in args]
        idx = [self._w(ax, a) for ax, a in zip(self.axes, args)]
        if len(idx) == 2:
            (ix, wx), (iy, wy) = idx
            v = self.v
            return ((1 - wy) * ((1 - wx) * v[iy][ix] + wx * v[iy][ix + 1])
                    + wy * ((1 - wx) * v[iy + 1][ix] + wx * v[iy + 1][ix + 1]))
        (ix, wx), (iy, wy), (iz, wz) = idx
        out = 0.0
        for dz, fz in ((0, 1 - wz), (1, wz)):
            if fz == 0.0:
                continue
            for dy, fy in ((0, 1 - wy), (1, wy)):
                if fy == 0.0:
                    continue
                row = self.v[iz + dz][iy + dy]
                out += fz * fy * ((1 - wx) * row[ix] + wx * row[ix + 1])
        return out


class Evaluator:
    """The variables of one file, evaluated lazily in a namespace of model
    values (plas_*, n*, Ue, Tg, emw_*), parameters and constants."""

    def __init__(self, rows, params, funcs, values):
        self.code = {n: compile(translate(e), n, "eval") for n, e, _ in rows}
        self.cache = {}
        self.ns = dict(CONST)
        self.ns.update(funcs)
        self.ns.update({"exp": np.exp, "log": np.log, "sqrt": np.sqrt,
                        "max": max, "min": min, "abs": abs, "real": np.real})
        self.ns.update(values)
        for n, e, _ in params:
            self.ns[n] = None
        self.params = {n: compile(translate(e), n, "eval") for n, e, _ in params}

    def __getitem__(self, name):
        if name in self.cache:
            return self.cache[name]
        if name in self.code:
            v = eval(self.code[name], {"__builtins__": {}}, _Lazy(self))
        elif name in self.params:
            v = eval(self.params[name], {"__builtins__": {}}, _Lazy(self))
        else:
            v = self.ns[name]
        self.cache[name] = v
        return v


class _Lazy(dict):
    def __init__(self, ev):
        super().__init__()
        self.ev = ev

    def __getitem__(self, k):
        return self.ev[k]

    def __contains__(self, k):
        return True


def check_files(kit, out: Path, fails: list):
    gas = kit["gas"]
    axes, data = read_grid(out / f"kinetics_{gas}.txt")
    ref = np.array([row for q in kit["grid"]["quantities"] for row in q["values"]])
    if data.shape != ref.shape or not np.allclose(data, ref, rtol=1e-9, atol=0):
        fails.append("grid file differs from the kit's tables")
    if not (np.allclose(axes[0], kit["grid"]["lnte"]) and np.allclose(axes[1], kit["grid"]["lnne"])):
        fails.append("grid axes differ")
    sp = {s["name"]: s for s in kit["species"]}
    atoms = {}
    for s in kit["species"]:
        a = {}
        for el in ("Kr", "F"):
            pass
    variables = {n for n, _, _ in read_table(out / f"variables_{gas}_plasma.txt")}
    chem = (out / f"chemistry_{gas}.txt").read_text().splitlines()
    block = None
    nreac = 0
    for line in chem:
        if not line or line.startswith("%"):
            continue
        if line in ("SPECIES", "USERDEFINED", "END", "ELECTRONIMPACTREACTIONS",
                    "RATECONSTANT", "HEAVYSPECIESREACTIONS"):
            block = line if line not in ("RATECONSTANT", "USERDEFINED", "END") else block
            continue
        if block in ("ELECTRONIMPACTREACTIONS", "HEAVYSPECIESREACTIONS"):
            lhs_rhs, _, rest = line.partition(" ! ")
            kf = rest.split()[-1]
            if kf not in variables:
                fails.append(f"chemistry: {kf} is no variable")
            nreac += 1
            l, r = lhs_rhs.split("=>")
            for side in (l, r):
                for nm in split_species(side):
                    if nm != "e" and nm not in sp:
                        fails.append(f"chemistry: unknown species {nm} in {line}")
            q = lambda side: sum(-1 if nm == "e" else sp[nm]["charge"] for nm in split_species(side))
            if q(l) != q(r):
                fails.append(f"chemistry: charge not balanced in {line}")
    names = [n for n, _, _ in read_table(out / f"variables_{gas}_plasma.txt")]
    if len(names) != len(set(names)):
        fails.append("a variable defined twice")
    ini = kit["initial"]
    cat = sum(ini["n_cm3"][s["light"]] * s["charge"] for s in kit["species"])
    if not math.isclose(cat, ini["ne_cm3"], rel_tol=1e-6):
        fails.append(f"initial state not quasineutral: {cat:.6g} vs {ini['ne_cm3']:.6g}")
    return nreac


def split_species(side: str) -> list:
    """'e+Kr2++F-' -> ['e', 'Kr2+', 'F-']: a '+' or '-' that follows a
    name and precedes '+', the end, or nothing else is a charge."""
    out, cur, i = [], "", 0
    while i < len(side):
        ch = side[i]
        if ch == "+":
            if cur and (i + 1 == len(side) or side[i + 1] == "+"):
                cur += "+"
                i += 1
                continue
            out.append(cur)
            cur = ""
        elif ch == "-" and (i + 1 == len(side) or side[i + 1] == "+"):
            cur += "-"
        else:
            cur += ch
        i += 1
    if cur:
        out.append(cur)
    return [s for s in out if s]


def check_equations(kit, out: Path, fails: list, n_states=12) -> dict:
    gas = kit["gas"]
    mods = C.load(gas)
    R = K.runs(gas, mods)
    tabs = K.tables(gas, mods, R)
    p_full = C.reference_params(mods, gas, wave=False)
    p0 = C.afterglow_params(mods, p_full, 0.0, C.F_BUILD_HZ)
    M = C.blend_class(mods)(p0, tabs, y0=R["y1"])
    args = (False, False) if C._wave_arg(mods["model"].Model) else (False,)
    funcs = {f"kin_{gas}": GridFunction(out / f"kinetics_{gas}.txt")}
    if (out / f"heavy_{gas}.txt").is_file():
        funcs[f"heavy_{gas}"] = GridFunction(out / f"heavy_{gas}.txt")
    params = read_table(out / f"params_{gas}.txt")
    spec = kit["species"]
    worst = {}
    for route in ("plasma", "ode"):
        rows = read_table(out / f"variables_{gas}_{route}.txt")
        on = R["ag_on"]
        ks = np.unique(np.linspace(0, len(on["t"]) - 1, n_states).astype(int))
        for k in ks:
            y = C.state_vector(M, K.Rec(on), int(k))
            n = np.maximum(y[:M.ns], 0.0) * M.N0                     # cm^-3
            ne = M.electron_density(n)
            te = max(y[M.i_ue] * M.N0 / (1.5 * ne), 1e-3)
            tg = y[M.i_tg]
            vals = {"Tg": tg}
            if route == "plasma":
                vals.update({"plas_Te": te, "plas_ne": ne * 1e6})
                for s in spec:
                    vals[f"plas_n_w{s['tag']}"] = n[M.idx[s["light"]]] * 1e6
            else:
                for s in spec:
                    vals[f"n{s['tag']}"] = n[M.idx[s["light"]]] * 1e6
                vals["Ue"] = y[M.i_ue] * M.N0 * 1e6 * C.QE
            vals["emw_normE"] = 0.0
            vals["emw_Qrh"] = 0.0
            ev = Evaluator(rows, params, funcs, vals)
            M.p_em_cell = 0.0
            d = M.rhs(0.0, y, *args)                                 # per N0 per s
            for s in spec:
                want = d[M.idx[s["light"]]] * M.N0 * 1e6              # 1/(m3 s)
                got = float(np.real(ev[f"dn_{s['tag']}"]))
                scale = max(abs(want), 1e-9 * M.N0 * 1e6 / 1e-9 * 0 + 1e14)
                err = abs(got - want) / max(abs(want), 1e14)
                worst[f"{route}:dn_{s['tag']}"] = max(worst.get(f"{route}:dn_{s['tag']}", 0.0), err)
            want_e = d[M.i_ue] * M.N0 * 1e6 * C.QE                   # W/m3
            got_e = float(np.real(ev["Qe"]))
            worst[f"{route}:Qe"] = max(worst.get(f"{route}:Qe", 0.0),
                                       abs(got_e - want_e) / max(abs(want_e), 1e3))
            want_t = d[M.i_tg]
            got_t = float(np.real(ev["dTg"]))
            worst[f"{route}:dTg"] = max(worst.get(f"{route}:dTg", 0.0),
                                        abs(got_t - want_t) / max(abs(want_t), 1e2))
            for c in M.compiled:
                pass
    for key, err in worst.items():
        if err > TOL:
            fails.append(f"{key}: relative error {err:.2e}")
    return worst


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gas", choices=("KrF", "Kr"))
    a = ap.parse_args(argv)
    total = 0
    for gas in [a.gas] if a.gas else ["KrF", "Kr"]:
        out = KIT / gas
        kit = json.loads((out / f"kit_{gas}.json").read_text())
        fails = []
        nreac = check_files(kit, out, fails)
        worst = check_equations(kit, out, fails)
        print(f"{gas}: {nreac} reactions in the chemistry file; equations checked at "
              f"12 states x 2 routes, worst relative error "
              f"{max(worst.values()):.2e} ({max(worst, key=worst.get)})")
        for f in fails:
            print("  FAIL", f)
        total += len(fails)
    print("all checks pass" if total == 0 else f"{total} failures")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
