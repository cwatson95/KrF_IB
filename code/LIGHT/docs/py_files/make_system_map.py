"""Numbers for the LIGHT v0 system map (system_map.tex), taken from the model.

The map carries no hand-typed number. This script runs the v0 reference
point (and the --omissions set at it), computes every energy flow drawn on
the map as a percent of the deposited energy, checks that each node of the
map conserves them, and writes system_map_numbers.tex as one
\\lmdef{key}{value} per number. test_docs.py fails when that file is stale.

    python3 LIGHT/docs/make_system_map.py              # rewrite the numbers
    cd LIGHT/docs && latexmk -pdf system_map.tex       # rebuild the PDF

Energies are booked at the ledger's formation energies (Kr* 9.915, Kr+
13.996, KrF(B) 6.15, Kr2F* 4.10, F 0.80 eV above ground Kr + F2), so
what enters a node equals what leaves it.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# marker walk to LIGHT/, the directory holding run_v0.py (no parent count)
LIGHT = next(p for p in Path(__file__).resolve().parents
             if (p / "run_v0.py").is_file())
sys.path.insert(0, str(LIGHT))

from v0 import __version__                                   # noqa: E402
from v0.model import Params, run                              # noqa: E402
from v0.network import GROUPS                                 # noqa: E402
from v0.report import (accounting, f2_burnout_time, reservoir,  # noqa: E402
                       row_dict)
from v0.species import e_form                                 # noqa: E402

OUT = Path(__file__).resolve().parent / "system_map_numbers.tex"
C1_ROWS = ("C1a", "C1b", "C1c")


def _sci(x: float) -> str:
    """1.64e17 -> '1.64\\times10^{17}' (three significant figures)."""
    e = int(math.floor(math.log10(abs(x))))
    return f"{x / 10 ** e:.2f}{{\\times}}10^{{{e}}}"


def _pow10(x: float) -> str:
    e = round(math.log10(x))
    if abs(x / 10 ** e - 1.0) > 1e-12:
        return _sci(x)
    return f"10^{{{e}}}"


def flows(params: Params = None) -> tuple:
    """(flows, context): every energy flow on the map as percent of E_dep,
    unrounded, plus the run it came from."""
    # the map draws the notes' list as written (network v0)
    p = params or Params(network="v0")
    res = run(p)
    a = accounting(res)
    e_dep = a["e_dep"]
    ev = a["events"]
    route = {c.row.rid: c.route for c in res.model.compiled}
    krs, krf, k2f = (reservoir(res, n) for n in ("Kr*", "KrF*", "Kr2F*"))
    e_krs, e_krp, e_krf, e_k2f, e_f = (
        e_form(s) for s in ("Kr*", "Kr+", "KrF(B)", "Kr2F*", "F"))

    def pct(x):
        return 100.0 * x / e_dep

    f = {
        # electron gas: in = power + Penning return, out = G1 C1 C2 G6
        "e.P": pct(e_dep),
        "e.pen": pct(ev["G2"] * route["G2"]["to_e"]),
        "e.G1": pct(ev["G1"] * route["G1"]["from_e"]),
        "e.C1": pct(sum(ev[r] * route[r]["from_e"] for r in C1_ROWS)),
        "e.C2": pct(ev["C2"] * route["C2"]["from_e"]),
        "e.G6": pct(a["heat_elastic"]),
        "e.dU": pct(a["stores"]["electron gas"]),
        # C1: the cascade delivers Kr*(5s) and heats the gas with the excess
        "C1.kr": pct(sum(ev[r] for r in C1_ROWS) * e_krs),
        "C1.heat": pct(sum(ev[r] * route[r]["heat"] for r in C1_ROWS)),
        # the Kr* reservoir
        "krs.in": pct(krs["formed"] * e_krs),
        "krs.G3": pct(ev["G3"] * e_krs),
        "krs.G2": pct(2 * ev["G2"] * e_krs),
        "krs.left": pct(krs["left"] * e_krs),
        "brs.G3": 100.0 * krs["losses"]["G3"] / krs["formed"],
        "brs.G2": 100.0 * krs["losses"]["G2"] / krs["formed"],
        # G3 harpoon and G2 Penning
        "G3.krf": pct(ev["G3"] * e_krf),
        "G3.F": pct(ev["G3"] * e_f),
        "G3.heat": pct(ev["G3"] * route["G3"]["heat"]),
        "G2.ion": pct(ev["G2"] * e_krp),
        "C2.ion": pct(ev["C2"] * e_krp),
        # the KrF(B) reservoir
        "krf.G4": pct(ev["G4"] * e_krf),
        "krf.G5": pct(ev["G5"] * e_krf),
        "krf.left": pct(krf["left"] * e_krf),
        "brf.G4": 100.0 * krf["losses"]["G4"] / krf["formed"],
        "brf.G5": 100.0 * krf["losses"]["G5"] / krf["formed"],
        "krf.eff": pct(ev["G3"] * route["G4"]["rad"]),
        # G4 248 nm, G5 trimer, C3 415 nm
        "G4.hv": pct(ev["G4"] * route["G4"]["rad"]),
        "G4.F": pct(ev["G4"] * e_f),
        "G4.heat": pct(ev["G4"] * route["G4"]["heat"]),
        "G5.k2f": pct(ev["G5"] * e_k2f),
        "G5.heat": pct(ev["G5"] * route["G5"]["heat"]),
        "k2f.C3": pct(ev["C3"] * e_k2f),
        "k2f.left": pct(k2f["left"] * e_k2f),
        "C3.hv": pct(ev["C3"] * route["C3"]["rad"]),
        "C3.F": pct(ev["C3"] * e_f),
        "C3.heat": pct(ev["C3"] * route["C3"]["heat"]),
        # the fate ledger at t_end
        "out.248": pct(a["bands"]["248 nm (KrF B-X)"]),
        "out.415": pct(a["bands"]["~415 nm (Kr2F*)"]),
        "out.heat": pct(a["heat_reactions"]),
        "out.elastic": pct(a["heat_elastic"]),
        "out.ions": pct(a["stores"]["ions (Kr+, Kr2+, F-)"]),
        "out.F": pct(a["stores"]["F atoms (F2 broken)"]),
        "out.other": pct(a["stores"]["Kr*"] + a["stores"]["KrF(B)+KrF(C)"]
                         + a["stores"]["Kr2F*"] + a["stores"]["electron gas"]),
    }
    _check_conservation(f, a["closure"])
    extra = {"res": res, "params": p, "krs_formed": krs["formed"],
             "krf_formed": krf["formed"]}
    return f, extra


def _check_conservation(f: dict, closure: float, tol: float = 1e-8) -> None:
    """Each node of the map must balance before a number is printed."""
    nodes = {
        "electron gas": (f["e.P"] + f["e.pen"],
                         f["e.G1"] + f["e.C1"] + f["e.C2"] + f["e.G6"]
                         + f["e.dU"]),
        "C1": (f["e.C1"], f["C1.kr"] + f["C1.heat"]),
        "Kr*": (f["e.G1"] + f["C1.kr"],
                f["krs.G3"] + f["krs.G2"] + f["krs.left"]),
        "G3": (f["krs.G3"], f["G3.krf"] + f["G3.F"] + f["G3.heat"]),
        "G2": (f["krs.G2"], f["G2.ion"] + f["e.pen"]),
        "KrF(B)": (f["G3.krf"], f["krf.G4"] + f["krf.G5"] + f["krf.left"]),
        "G4": (f["krf.G4"], f["G4.hv"] + f["G4.F"] + f["G4.heat"]),
        "G5": (f["krf.G5"], f["G5.k2f"] + f["G5.heat"]),
        "Kr2F*": (f["G5.k2f"], f["k2f.C3"] + f["k2f.left"]),
        "C3": (f["k2f.C3"], f["C3.hv"] + f["C3.F"] + f["C3.heat"]),
        "ions": (f["C2.ion"] + f["G2.ion"], f["out.ions"]),
        "F atoms": (f["G3.F"] + f["G4.F"] + f["C3.F"], f["out.F"]),
        "heat": (f["C1.heat"] + f["G3.heat"] + f["G4.heat"] + f["G5.heat"]
                 + f["C3.heat"], f["out.heat"]),
        "ledger": (100.0, f["out.248"] + f["out.415"] + f["out.heat"]
                   + f["out.elastic"] + f["out.ions"] + f["out.F"]
                   + f["out.other"]),
    }
    for name, (into, out) in nodes.items():
        if abs(into - out) > tol:
            raise AssertionError(f"map node {name!r} does not balance: "
                                 f"{into:.12f} in, {out:.12f} out")
    if abs(closure) > 1e-12:
        raise AssertionError(f"ledger closure {closure:.2e} above 1e-12")


def numbers(params: Params = None) -> dict:
    """The formatted strings the map prints, keyed as in the .tex."""
    f, x = flows(params)
    p, res = x["params"], x["res"]
    n = {k: f"{v:.1f}" for k, v in f.items()}
    on = res.t <= p.pulse_s
    tb = f2_burnout_time(res)
    n.update({
        "version": __version__,
        "p": f"{p.pressure_atm:g}",
        "xF2": f"{100 * p.x_f2:g}",
        "Tg0": f"{p.tg0_K:.0f}",
        "P": f"{p.power_W_cm3 / 1e3:.0f}",
        "pulse": f"{p.pulse_s * 1e6:g}",
        "tend": f"{p.t_end_s * 1e6:g}",
        "ne": _pow10(p.ne_cm3),
        "N": _sci(res.model.N0),
        "tau": f"{p.tau_sp_s * 1e9:.2f}",
        "Edep": f"{res.ledger['E_dep'][-1] * 1.602176634e-19:.2f}",
        "Te": f"{float(res.te[on].mean()):.2f}",
        "TgEnd": f"{res.tg[-1]:.0f}",
        "pEnd": f"{res.pressure_atm[-1]:.2f}",
        "F2burn": f"{tb * 1e9:.0f}" if tb is not None else "none",
        "F2left": f"{100 * res.n['F2'][-1] / res.n['F2'][0]:.1f}",
        "krsFormed": _sci(x["krs_formed"]),
        "krfFormed": _sci(x["krf_formed"]),
        # the notes' output line (flows, they overlap)
        "notes.krs": f"{f['e.G1'] + f['e.C1']:.1f}",
        "notes.krp": f"{f['C2.ion'] + f['G2.ion']:.1f}",
        "notes.emit": f"{f['out.248'] + f['out.415']:.1f}",
    })
    # what v0 leaves out: the 248-nm share with each optional group on
    for g in GROUPS + ("all",):
        d = row_dict(run(p.replace(groups=(g,))))
        n[f"om.{g}"] = f"{d['emission_248']:.1f}"
    n["om.v0"] = n["out.248"]
    return n


def render(n: dict) -> str:
    lines = [
        "% LIGHT system-map numbers. GENERATED by make_system_map.py from the",
        f"% v0 model ({n['version']}). Do not edit. Rerun the script instead.",
        "\\providecommand\\lmdef[2]{\\expandafter\\def"
        "\\csname lm@#1\\endcsname{#2}}",
    ]
    lines += [f"\\lmdef{{{k}}}{{{v}}}" for k, v in n.items()]
    return "\n".join(lines) + "\n"


def main() -> int:
    OUT.write_text(render(numbers()))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
