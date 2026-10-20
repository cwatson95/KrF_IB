"""Numbers for LIGHT_ArF_ebeam v0_2's system map (system_map.tex).

LIGHT's make_system_map_v0_2.py for the ArF network. The map draws v0_2's
basic network at the reference point with every band at its default: every
energy flow as a percent of the deposited energy, each node checked to
conserve them (to 1e-8 of E_dep) before a number is printed, every row of
the network drawn. It writes system_map_numbers.tex as one
\\wmdef{key}{value} per number. v0_2/tests/test_docs.py checks the balance
and, with LIGHT_DOCS_SLOW=1, the file.

    python3 LIGHT_ArF_ebeam/docs/make_system_map.py            # about 7 minutes
    cd LIGHT_ArF_ebeam/docs && /opt/homebrew/bin/latexmk -pdf system_map.tex

What differs from KrF's map. ArF(B) is the one excimer state, so there is no
C node and no mixing. The 37 argon excitation rows are drawn as two
families, the four 4s rows (G1, whose three upper levels leave their excess
over the lumped Ar* as heat) and the 33 higher rows (C1, which cascade to
Ar*). Four rows have no measurement and are 0 by default: Ar* + Ar* (G2),
Ar2* + F (K2F_F), Ar2* + Ar2* (PEN_AR2) and Ar2F* + Ar (Q2AR). The map
prints them and draws no band, and this script refuses to run when one of
them carries energy, since a band set to a nonzero value needs its own
drawing. The box at the bottom gives the 193-nm share at the two ends of
every band (run_v0_2.py --bands).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIGHT_ARF = HERE.parent
if str(LIGHT_ARF) not in sys.path:
    sys.path.insert(0, str(LIGHT_ARF))

from v0_2 import __version__                                  # noqa: E402
from v0_2.model import Params, run                            # noqa: E402
from v0_2.network import BANDS, BASIC_GROUPS                  # noqa: E402
from v0_2.report import (accounting, f2_burnout_time,         # noqa: E402
                         reservoir, row_dict)
from v0_2.run import do_bands                                 # noqa: E402
from v0_2.species import e_form                               # noqa: E402

OUT = HERE / "system_map_numbers.tex"
HARPOON = ("G3", "G3X")
QUENCH_B = ("QF2", "QAR", "QE")
F2E = ("DIS_a3", "DIS_A1", "VIB1", "VIB2", "VIB3", "VIB4")
ION_ION = ("K26", "K108")
AR2_F2 = ("NK7", "K2F_BR")
CATIONS = ("Ar+", "Ar2+")
# the rows with no measurement, 0 unless a band run sets them: the map
# prints them and draws no band
ZERO_ROWS = ("G2", "K2F_F", "PEN_AR2", "Q2AR")
BANDS_OUT = {"193": "193 nm (ArF B-X)", "284": "284 nm (Ar2F*)",
             "127": "127 nm (Ar2* VUV)"}
# the order the map's box lists the bands in: the levers first
BAND_ORDER = ("ar2_f", "arf_yield", "arf_trimer", "trimer_tg_pow",
              "ar2_form", "ar_f2", "ar2_rad", "ar2_f2", "ar2_f2_arf",
              "arf_e", "ar2f_f", "ar2f_ar", "ar_pen", "ar2_pen")


def _sci(x: float) -> str:
    """1.64e17 -> '1.64\\times10^{17}' (three significant figures)."""
    e = int(math.floor(math.log10(abs(x))))
    return f"{x / 10 ** e:.2f}{{\\times}}10^{{{e}}}"


def _pow10(x: float) -> str:
    e = round(math.log10(x))
    if abs(x / 10 ** e - 1.0) > 1e-12:
        return _sci(x)
    return f"10^{{{e}}}"


def _k(x: float) -> str:
    """A rate constant as printed: at most three significant figures and at
    least two, trailing zeros beyond two dropped."""
    mant, exp = f"{x:.2e}".split("e")
    mant = mant.rstrip("0")
    if mant.endswith("."):
        mant += "0"
    return f"{mant}{{\\times}}10^{{{int(exp)}}}"


def _g(x: float, sig: int) -> str:
    """sig significant figures without an exponent (203.3 -> 203)."""
    return f"{float(f'{x:.{sig}g}'):g}"


def _f1(x: float) -> str:
    """One decimal, with no negative zero."""
    s = f"{x:.1f}"
    return "0.0" if s == "-0.0" else s


def _net_f(row) -> int:
    """F atoms a row makes, net of those it consumes."""
    return row.products.count("F") - row.reactants.count("F")


def _check_defaults(p: Params) -> None:
    if p.network != "basic" or p.groups:
        raise ValueError("the map draws the basic network with no optional "
                         "group")
    if p.rates or p.trimer_tg_pow is not None:
        raise ValueError("the map draws every band at its default; a band "
                         "run needs its own drawing")


def flows(params: Params = None) -> tuple:
    """(flows, context): every energy flow on the map as percent of E_dep,
    unrounded, plus the run it came from."""
    p = params or Params()
    _check_defaults(p)
    res = run(p)
    if res.model.ne_mode != "free":
        raise ValueError("the map draws free n_e (ne_mode 'auto')")
    a = accounting(res)
    e_dep = a["e_dep"]
    ev = a["events"]
    rows = {c.row.rid: c.row for c in res.model.compiled}
    route = {c.row.rid: c.route for c in res.model.compiled}
    E = {s: e_form(s) for s in ("Ar*", "Ar+", "Ar2+", "F-", "ArF(B)",
                                "Ar2F*", "F", "Ar2*")}
    g1 = tuple(r for r in rows if r.startswith("G1_"))
    c1 = tuple(r for r in rows if r.startswith("C1_"))
    ars = reservoir(res, "Ar*")
    a2s = reservoir(res, "Ar2*")
    a2f = reservoir(res, "Ar2F*")
    used = set()
    nonzero = [r for r in ZERO_ROWS if ev[r] != 0.0]
    if nonzero:
        raise AssertionError(f"rows the map prints as 0 carry energy: "
                             f"{nonzero}. Redraw them in system_map.tex")

    def pct(x):
        return 100.0 * x / e_dep

    def ev_(*rids):
        used.update(rids)
        return sum(ev[r] for r in rids)

    def from_e(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["from_e"] for r in rids))

    def heat(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["heat"] for r in rids))

    def f_atoms(*rids):
        used.update(rids)
        return pct(sum(ev[r] * _net_f(rows[r]) for r in rids) * E["F"])

    def rad(rid):
        used.add(rid)
        return pct(ev[rid] * route[rid]["rad"])

    def left(*species):
        return pct(sum((res.n[s][-1] - res.n[s][0]) * E[s] for s in species))

    def ion_in(rid):
        """Positive-ion energy a row consumes (per the ledger)."""
        return sum(rows[rid].reactants.count(s) * E[s] for s in CATIONS)

    f = {
        # electron gas: in = power + the two Penning returns (0 here), out =
        # the e-impact rows, elastic loss, the captured electrons' energy
        # and the store
        "e.P": pct(e_dep),
        "e.pen": pct(ev_("G2") * route["G2"]["to_e"]),
        "e.pen2": pct(ev_("PEN_AR2") * route["PEN_AR2"]["to_e"]),
        "e.G1": from_e(*g1),
        "e.C1": from_e(*c1),
        "e.C2": from_e("C2"),
        "e.F2": from_e(*F2E),
        "e.G6": pct(a["heat_elastic"]),
        "e.dU": pct(a["stores"]["electron gas"]),
        "e.cap": pct(float(res.ledger["heat_capture"][-1])),
        # G1: the 4s levels make Ar*, the three above 3P2 heat the gas
        "G1.ar": pct(ev_(*g1) * E["Ar*"]),
        "G1.heat": heat(*g1),
        # C1: the cascade delivers Ar* and heats the gas with the excess
        "C1.ar": pct(ev_(*c1) * E["Ar*"]),
        "C1.heat": heat(*c1),
        # F2 under electron impact: dissociation, vibration and attachment
        "F2e.F": f_atoms(*F2E, "ATT"),
        "F2e.heat": heat(*F2E, "ATT"),
        "F2e.anion": pct(ev_("ATT") * E["F-"]),
        # the Ar* reservoir, fed by excitation and by recombination
        "ars.in": pct(ars["formed"] * E["Ar*"]),
        "ars.harp": pct(ev_(*HARPOON) * E["Ar*"]),
        "ars.G2": pct(2 * ev_("G2") * E["Ar*"]),
        "ars.AK8": pct(ev_("AK8_Ar") * E["Ar*"]),
        "ars.left": left("Ar*"),
        "brs.harp": 100.0 * sum(ars["losses"].get(r, 0.0) for r in HARPOON)
        / ars["formed"],
        "brs.AK8": 100.0 * ars["losses"]["AK8_Ar"] / ars["formed"],
        # Ar* + F2: ArF(B), or no ArF at all
        "harp.B": pct(ev_("G3") * E["ArF(B)"]),
        "harp.F": f_atoms(*HARPOON),
        "harp.heat": heat(*HARPOON),
        # Penning (0) and direct ionization into the positive ions
        "G2.ion": pct(ev_("G2") * E["Ar+"]),
        "C2.ion": pct(ev_("C2") * E["Ar+"]),
        # the Ar2* node: Ar* + 2Ar (AK8_Ar) inside it, then 127-nm light,
        # Ar2* + F2, and the two rows at 0
        "AK8.heat": heat("AK8_Ar"),
        "a2s.in": pct(a2s["formed"] * E["Ar2*"]),
        "a2s.RAD": pct(ev_("AR2_RAD") * E["Ar2*"]),
        "a2s.F2": pct(ev_(*AR2_F2) * E["Ar2*"]),
        "a2s.F": pct(ev_("K2F_F") * E["Ar2*"]),
        "a2s.pen": pct(2 * ev_("PEN_AR2") * E["Ar2*"]),
        "a2s.left": left("Ar2*"),
        "bra2.RAD": 100.0 * a2s["losses"]["AR2_RAD"] / a2s["formed"],
        "bra2.F2": 100.0 * sum(a2s["losses"][r] for r in AR2_F2)
        / a2s["formed"],
        "RAD.hv": rad("AR2_RAD"),
        "RAD.heat": heat("AR2_RAD"),
        # Ar2* + F2: ArF(B) + Ar + F (NK7) or Ar2F* + F (K2F_BR)
        "NK7.B": pct(ev["NK7"] * E["ArF(B)"]),
        "BR.a2f": pct(ev["K2F_BR"] * E["Ar2F*"]),
        "A2F2.F": f_atoms(*AR2_F2),
        "A2F2.heat": heat(*AR2_F2),
        # Ar2* + F -> ArF(B) + Ar (0): its F atom comes in from the F store
        "a2F.B": pct(ev["K2F_F"] * E["ArF(B)"]),
        "a2F.Fin": -f_atoms("K2F_F"),
        "a2F.heat": heat("K2F_F"),
        # Ar2* + Ar2* -> Ar2+ + 2Ar + e (0)
        "pen2.ion": pct(ev["PEN_AR2"] * E["Ar2+"]),
        # the positive ions: Ar+ -> Ar2+ (K15), Ar2+ + e -> Ar* (R11p),
        # ion-ion recombination with F- into ArF(B) (K26, K108)
        "K15.heat": heat("K15"),
        "R11p.in": pct(ev_("R11p") * ion_in("R11p")),
        "R11p.ar": pct(ev["R11p"] * E["Ar*"]),
        "R11p.heat": heat("R11p"),
        "ii.in": pct(sum(ev_(r) * ion_in(r) for r in ION_ION)),
        "ii.anion": pct(ev_(*ION_ION) * E["F-"]),
        "ii.B": pct(ev_(*ION_ION) * E["ArF(B)"]),
        "ii.heat": heat(*ION_ION),
        "ions.left": left(*CATIONS),
        "anion.left": left("F-"),
        # the ArF(B) store
        "B.G4": pct(ev_("G4") * E["ArF(B)"]),
        "B.G5": pct(ev_("G5") * E["ArF(B)"]),
        "B.Q": pct(ev_(*QUENCH_B) * E["ArF(B)"]),
        "B.QF2": pct(ev["QF2"] * E["ArF(B)"]),
        "B.QAR": pct(ev["QAR"] * E["ArF(B)"]),
        "B.QE": pct(ev["QE"] * E["ArF(B)"]),
        "B.left": left("ArF(B)"),
        # 193-nm emission, and the quenching of B by F2, Ar and electrons
        "G4.hv": rad("G4"),
        "G4.F": f_atoms("G4"),
        "G4.heat": heat("G4"),
        "QB.F": f_atoms(*QUENCH_B),
        "QB.heat": heat(*QUENCH_B),
        # the trimer
        "tri.a2f": pct(ev_("G5") * E["Ar2F*"]),
        "tri.heat": heat("G5"),
        # the Ar2F* store, its light and its quenchers
        "a2f.C3": pct(ev_("C3") * E["Ar2F*"]),
        "a2f.Q2F": pct(ev_("Q2F") * E["Ar2F*"]),
        "a2f.Q2FF": pct(ev_("Q2FF") * E["Ar2F*"]),
        "a2f.Q2AR": pct(ev_("Q2AR") * E["Ar2F*"]),
        "a2f.left": left("Ar2F*"),
        "brt.C3": 100.0 * a2f["losses"]["C3"] / a2f["formed"],
        "brt.Q2F": 100.0 * a2f["losses"]["Q2F"] / a2f["formed"],
        "brt.Q2FF": 100.0 * a2f["losses"]["Q2FF"] / a2f["formed"],
        "C3.hv": rad("C3"),
        "C3.F": f_atoms("C3"),
        "C3.heat": heat("C3"),
        "Q2F.F": f_atoms("Q2F"),
        "Q2F.heat": heat("Q2F"),
        "Q2FF.F": f_atoms("Q2FF"),
        "Q2FF.heat": heat("Q2FF"),
        "Q2AR.F": f_atoms("Q2AR"),
        "Q2AR.heat": heat("Q2AR"),
        # F + F + Ar -> F2 + Ar: F atoms back to F2, their energy to heat
        "FF.in": -f_atoms("FF"),
        "FF.heat": heat("FF"),
        # the fate ledger at t_end
        "out.193": pct(a["bands"][BANDS_OUT["193"]]),
        "out.284": pct(a["bands"][BANDS_OUT["284"]]),
        "out.127": pct(a["bands"][BANDS_OUT["127"]]),
        "out.heat": pct(a["heat_reactions"]),
        "out.elastic": pct(a["heat_elastic"]),
        "out.ions": pct(a["stores"]["ions (Ar+, Ar2+, F-)"]),
        "out.F": pct(a["stores"]["F atoms (F2 broken)"]),
        "out.other": pct(a["stores"]["Ar*"] + a["stores"]["Ar2*"]
                         + a["stores"]["ArF(B)"] + a["stores"]["Ar2F*"]
                         + a["stores"]["electron gas"]),
    }
    used.add("G6")                        # G6 is the elastic loss, e.G6
    missing = set(rows) - used
    if missing:
        raise AssertionError(f"rows not drawn on the map: {sorted(missing)}")
    if set(a["bands"]) != set(BANDS_OUT.values()):
        raise AssertionError(f"bands changed: {sorted(a['bands'])}")
    _check_conservation(f, a["closure"])
    arf = reservoir(res, "ArF*")
    extra = {"res": res, "params": p, "ars": ars, "arf": arf, "a2f": a2f,
             "a2s": a2s, "rows": rows, "events": ev}
    return f, extra


def _check_conservation(f: dict, closure: float, tol: float = 1e-8) -> None:
    """Each node of the map must balance before a number is printed."""
    for name, (into, out) in nodes(f).items():
        if abs(into - out) > tol:
            raise AssertionError(f"map node {name!r} does not balance: "
                                 f"{into:.12f} in, {out:.12f} out")
    if abs(closure) > 1e-12:
        raise AssertionError(f"ledger closure {closure:.2e} above 1e-12")


def nodes(f: dict) -> dict:
    """Every node of the map as (energy in, energy out), percent of E_dep."""
    return {
        "electron gas": (f["e.P"] + f["e.pen"] + f["e.pen2"],
                         f["e.G1"] + f["e.C1"] + f["e.C2"] + f["e.F2"]
                         + f["e.G6"] + f["e.dU"] + f["e.cap"]),
        "G1": (f["e.G1"], f["G1.ar"] + f["G1.heat"]),
        "C1": (f["e.C1"], f["C1.ar"] + f["C1.heat"]),
        "F2 electron rows": (f["e.F2"], f["F2e.F"] + f["F2e.heat"]
                             + f["F2e.anion"]),
        "Ar*": (f["G1.ar"] + f["C1.ar"] + f["R11p.ar"],
                f["ars.harp"] + f["ars.G2"] + f["ars.AK8"] + f["ars.left"]),
        "Ar* in": (f["ars.in"], f["G1.ar"] + f["C1.ar"] + f["R11p.ar"]),
        "harpoon": (f["ars.harp"], f["harp.B"] + f["harp.F"]
                    + f["harp.heat"]),
        "G2": (f["ars.G2"], f["G2.ion"] + f["e.pen"]),
        "Ar2* in": (f["a2s.in"], f["ars.AK8"] - f["AK8.heat"]),
        "Ar2* (with AK8_Ar)": (f["ars.AK8"],
                               f["AK8.heat"] + f["a2s.RAD"] + f["a2s.F2"]
                               + f["a2s.F"] + f["a2s.pen"] + f["a2s.left"]),
        "AR2_RAD": (f["a2s.RAD"], f["RAD.hv"] + f["RAD.heat"]),
        "Ar2* + F2": (f["a2s.F2"], f["NK7.B"] + f["BR.a2f"] + f["A2F2.F"]
                      + f["A2F2.heat"]),
        "Ar2* + F": (f["a2s.F"] + f["a2F.Fin"], f["a2F.B"] + f["a2F.heat"]),
        "Ar2* + Ar2*": (f["a2s.pen"], f["pen2.ion"] + f["e.pen2"]),
        "positive ions": (f["C2.ion"] + f["G2.ion"] + f["pen2.ion"],
                          f["K15.heat"] + f["R11p.in"] + f["ii.in"]
                          + f["ions.left"]),
        "R11p": (f["R11p.in"], f["R11p.ar"] + f["R11p.heat"]),
        "F-": (f["F2e.anion"], f["ii.anion"] + f["anion.left"]),
        "ion-ion": (f["ii.in"] + f["ii.anion"], f["ii.B"] + f["ii.heat"]),
        "ArF(B)": (f["harp.B"] + f["ii.B"] + f["NK7.B"] + f["a2F.B"],
                   f["B.G4"] + f["B.G5"] + f["B.Q"] + f["B.left"]),
        "B quench split": (f["B.Q"], f["B.QF2"] + f["B.QAR"] + f["B.QE"]),
        "G4": (f["B.G4"], f["G4.hv"] + f["G4.F"] + f["G4.heat"]),
        "B quench": (f["B.Q"], f["QB.F"] + f["QB.heat"]),
        "trimer": (f["B.G5"], f["tri.a2f"] + f["tri.heat"]),
        "Ar2F*": (f["tri.a2f"] + f["BR.a2f"],
                  f["a2f.C3"] + f["a2f.Q2F"] + f["a2f.Q2FF"] + f["a2f.Q2AR"]
                  + f["a2f.left"]),
        "C3": (f["a2f.C3"], f["C3.hv"] + f["C3.F"] + f["C3.heat"]),
        "Q2F": (f["a2f.Q2F"], f["Q2F.F"] + f["Q2F.heat"]),
        "Q2FF": (f["a2f.Q2FF"], f["Q2FF.F"] + f["Q2FF.heat"]),
        "Q2AR": (f["a2f.Q2AR"], f["Q2AR.F"] + f["Q2AR.heat"]),
        "FF": (f["FF.in"], f["FF.heat"]),
        "ions": (f["ions.left"] + f["anion.left"], f["out.ions"]),
        "F atoms": (f["harp.F"] + f["G4.F"] + f["QB.F"] + f["C3.F"]
                    + f["Q2F.F"] + f["Q2FF.F"] + f["Q2AR.F"] + f["F2e.F"]
                    + f["A2F2.F"],
                    f["FF.in"] + f["a2F.Fin"] + f["out.F"]),
        "heat": (f["G1.heat"] + f["C1.heat"] + f["harp.heat"] + f["G4.heat"]
                 + f["QB.heat"] + f["tri.heat"] + f["C3.heat"]
                 + f["Q2F.heat"] + f["Q2FF.heat"] + f["Q2AR.heat"]
                 + f["e.cap"] + f["F2e.heat"] + f["K15.heat"]
                 + f["R11p.heat"] + f["ii.heat"] + f["FF.heat"]
                 + f["AK8.heat"] + f["RAD.heat"] + f["A2F2.heat"]
                 + f["a2F.heat"],
                 f["out.heat"]),
        "193 nm": (f["G4.hv"], f["out.193"]),
        "284 nm": (f["C3.hv"], f["out.284"]),
        "127 nm": (f["RAD.hv"], f["out.127"]),
        "other stores": (f["ars.left"] + f["a2s.left"] + f["B.left"]
                         + f["a2f.left"] + f["e.dU"], f["out.other"]),
        "ledger": (100.0, f["out.193"] + f["out.284"] + f["out.127"]
                   + f["out.heat"] + f["out.elastic"] + f["out.ions"]
                   + f["out.F"] + f["out.other"]),
    }


def _rate(res, rid) -> float:
    """The rate one row ran with, as compiled (a band row's resolved
    value, a radiative row's share of its lifetime)."""
    c = next(c for c in res.model.compiled if c.row.rid == rid)
    if c.kind not in ("const", "param", "inv_tau"):
        raise ValueError(f"{rid}: rate kind {c.kind!r} has no single value")
    return float(c.rate)


def _band_value(b, v) -> str:
    """A band end as the map prints it (math)."""
    if b.unit in ("exponent", "fraction"):
        return f"{v:.2f}".replace("-", "{-}")
    return "0" if v == 0.0 else _k(v)


def bands_numbers(p: Params) -> dict:
    """The 193-nm share of E_dep at the two ends of every band, one band at
    a time (run_v0_2.py --bands), for the map's box."""
    rows = do_bands(p)
    by = {(r["band"], r["end"]): r for r in rows[1:]}
    if set(BAND_ORDER) != set(BANDS):
        raise AssertionError("the bands changed: update BAND_ORDER and the "
                             "map's box")
    n = {}
    for name in BAND_ORDER:
        lo, hi = by[(name, "low")], by[(name, "high")]
        n[f"bd.{name}.lo193"] = f"{lo['emission_193']:.2f}"
        n[f"bd.{name}.hi193"] = f"{hi['emission_193']:.2f}"
    n["bd.n"] = str(len(BANDS))
    return n


def numbers(params: Params = None, with_bands: bool = True) -> dict:
    """The formatted strings the map prints, keyed as in the .tex."""
    f, x = flows(params)
    p, res, rows, ev = x["params"], x["res"], x["rows"], x["events"]
    n = {k: _f1(v) for k, v in f.items()}
    # the Ar2* fates box draws AR2_RAD and Ar2* + F2 together: its input
    # band and its heat tag add the two
    n["a2s.fates"] = _f1(f["a2s.RAD"] + f["a2s.F2"])
    n["ARx.heat"] = _f1(f["RAD.heat"] + f["A2F2.heat"])
    # the ion-ion band carries the ArF(B) of Ar2* + F2 above the point
    # where the two meet
    n["iiNK7.B"] = _f1(f["ii.B"] + f["NK7.B"])
    # the Q2 box draws Ar2F* + F2 and Ar2F* + F together (Q2AR is 0)
    n["a2f.Q2x"] = _f1(f["a2f.Q2F"] + f["a2f.Q2FF"])
    n["Q2x.F"] = _f1(f["Q2F.F"] + f["Q2FF.F"])
    n["Q2x.heat"] = _f1(f["Q2F.heat"] + f["Q2FF.heat"])
    on = res.t <= p.pulse_s
    tb = f2_burnout_time(res)
    m = res.model
    # the electron closure at mid-pulse: the tail temperature T_t = r Te
    k = int(res.t.searchsorted(0.5 * p.pulse_s))
    n_gas = float(sum(res.n[s][k] for s in m.species))
    r_tail = m.tail(res.ne[k] / n_gas, res.n["F2"][k] / n_gas)
    arf = x["arf"]
    lost = arf["losses"]

    def brf(*rids):
        return _f1(100.0 * sum(lost.get(r, 0.0) for r in rids)
                   / arf["formed"])
    k_h = _rate(res, "G3") + _rate(res, "G3X")
    y_h = _rate(res, "G3") / k_h
    formed_ii = sum(ev[r] for r in ION_ION)
    formed_a2 = ev["NK7"] + ev["K2F_F"]
    formed_all = ev["G3"] + formed_ii + formed_a2
    if abs(formed_all - arf["formed"]) > 1e-9 * arf["formed"]:
        raise AssertionError("the ArF sources on the map miss one")
    kind, k0, te0, n_pow = rows["R11p"].rate
    assert kind == "te_pow", rows["R11p"].rate
    k_f2 = _rate(res, "NK7") + _rate(res, "K2F_BR")
    burn = tb if tb is not None else p.pulse_s
    before = res.t <= burn
    i_b = int(res.t.searchsorted(burn))
    ar_in = {rid: float(res.events[rid][i_b]) for rid in
             ("AK8_Ar",) + HARPOON}
    ar_in["G2"] = 2.0 * float(res.events["G2"][i_b])
    g1 = sorted(r.e_cost_eV for r in rows.values() if r.rid.startswith("G1_"))
    c1 = sorted(r.e_cost_eV for r in rows.values() if r.rid.startswith("C1_"))
    n.update({
        "version": __version__.replace("_", r"\_"),
        "p": f"{p.pressure_atm:g}",
        "xF2": f"{100 * p.x_f2:g}",
        "Tg0": f"{p.tg0_K:.0f}",
        "P": f"{p.power_W_cm3 / 1e3:.0f}",
        "pulse": f"{p.pulse_s * 1e6:g}",
        "tend": f"{p.t_end_s * 1e6:g}",
        "ne": _pow10(p.ne_cm3),
        # n_e while F2 lasts (time mean up to the burnout) and its peak
        "neF2": _sci(float(res.ne[before].mean())),
        "nePeak": _sci(float(res.ne.max())),
        "N": _sci(m.N0),
        "Edep": f"{res.ledger['E_dep'][-1] * 1.602176634e-19:.2f}",
        "Te": f"{float(res.te[on].mean()):.2f}",
        "TeMid": f"{float(res.te[k]):.2f}",
        "rTail": f"{r_tail:.2f}",
        "TtMid": f"{r_tail * float(res.te[k]):.2f}",
        "TgEnd": f"{res.tg[-1]:.0f}",
        "pEnd": f"{res.pressure_atm[-1]:.2f}",
        "F2burn": f"{tb * 1e9:.0f}" if tb is not None else "none",
        "arsFormed": _sci(x["ars"]["formed"]),
        "a2sFormed": _sci(x["a2s"]["formed"]),
        "arfFormed": _sci(arf["formed"]),
        "a2fFormed": _sci(x["a2f"]["formed"]),
        "arf.eff": _f1(100.0 * arf["formed"] * rows["G4"].hv_eV
                       / float(res.ledger["E_dep"][-1])),
        # where ArF(B) comes from, percent of the ArF formed
        "frm.harp": _f1(100.0 * ev["G3"] / formed_all),
        "frm.ii": _f1(100.0 * formed_ii / formed_all),
        "frm.ar2": _f1(100.0 * formed_a2 / formed_all),
        # where ArF(B) goes, percent of the ArF formed
        "brf.G4": brf("G4"), "brf.tri": brf("G5"), "brf.e": brf("QE"),
        "brf.F2": brf("QF2"), "brf.Ar": brf("QAR"),
        "brf.QB": brf(*QUENCH_B),
        # how the ions recombine: with an electron (R11p) or with F-
        "ion.R11p": _f1(100.0 * ev["R11p"] / (ev["R11p"] + formed_ii)),
        "ion.ii": _f1(100.0 * formed_ii / (ev["R11p"] + formed_ii)),
        # Ar* into Ar2* before the F2 burns out, percent of the Ar* used
        # up to then
        "brs.AK8.F2": _f1(100.0 * ar_in["AK8_Ar"] / sum(ar_in.values())),
        # the constants the boxes print, read from the rows the model ran
        "k.harp": _k(k_h),
        "pc.harpY": f"{100 * y_h:.0f}",
        "pc.harpX": f"{100 * (1 - y_h):.0f}",
        "tauB": f"{p.tau_sp_s * 1e9:.1f}",
        "tauC3": _g(1e9 / _rate(res, "C3"), 3),
        "k.G5": _k(_rate(res, "G5")),
        "k.Q2F": _k(_rate(res, "Q2F")),
        "k.K15": _k(_rate(res, "K15")),
        "tpow.K15": f"{rows['K15'].tg_pow:.2f}".replace("-", "{-}"),
        "k.R11p": _k(k0),
        "pow.R11p": f"{-n_pow:.2f}",
        "k.FF": _k(_rate(res, "FF")),
        "n3.ii": f"{rows['K26'].rate[4]:g}".replace("-", "{-}"),
        "k.AK8": _k(_rate(res, "AK8_Ar")),
        "tauRAD": _g(1e6 / _rate(res, "AR2_RAD"), 2),
        "k.ar2F2": _k(k_f2),
        "pc.NK7": f"{100 * _rate(res, 'NK7') / k_f2:.0f}",
        "E.ar2s": f"{e_form('Ar2*'):.3f}",
        "E.ar2s.bind": f"{e_form('Ar*') - e_form('Ar2*'):.2f}",
        "hv.RAD": f"{rows['AR2_RAD'].hv_eV:.2f}",
        # how far F- lies below Ar + F2 on the ledger's scale (a magnitude)
        "E.anion": f"{-e_form('F-'):.1f}",
        "hv.G4": f"{rows['G4'].hv_eV:.2f}",
        "hv.C3": f"{rows['C3'].hv_eV:.2f}",
        "E.B": f"{e_form('ArF(B)'):.3f}",
        "E.X": f"{e_form('Ar*'):.3f}",
        "E.P": f"{e_form('Ar+'):.2f}",
        "g1.lo": f"{g1[0]:.3f}", "g1.hi": f"{g1[-1]:.3f}", "nG1": str(len(g1)),
        "c1.lo": f"{c1[0]:.3f}", "c1.hi": f"{c1[-1]:.3f}", "nC1": str(len(c1)),
        # the notes' output line (flows, they overlap)
        "notes.ars": _f1(f["e.G1"] + f["e.C1"]),
        "notes.arp": _f1(f["C2.ion"] + f["G2.ion"] + f["pen2.ion"]),
        "notes.emit": _f1(f["out.193"] + f["out.284"] + f["out.127"]),
        "nodes": str(len(nodes(f))),
    })
    # what the basic network leaves out: the 193-nm share with the electron
    # group on (the only group v0_2 has)
    for g in BASIC_GROUPS:
        d = row_dict(run(p.replace(groups=(g,))))
        n[f"om.{g}"] = f"{d['emission_193']:.2f}"
    n["om.basic"] = f"{row_dict(res)['emission_193']:.2f}"
    if with_bands:
        n.update(bands_numbers(p))
    return n


def render(n: dict) -> str:
    lines = [
        "% LIGHT_ArF_ebeam v0_2 map numbers. GENERATED by make_system_map.py",
        f"% from the model ({n['version']}). Do not edit. Rerun the script "
        "instead.",
        "\\providecommand\\wmdef[2]{\\expandafter\\def"
        "\\csname wm@#1\\endcsname{#2}}",
    ]
    lines += [f"\\wmdef{{{k}}}{{{v}}}" for k, v in n.items()]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if system_map_numbers.tex is stale")
    a = ap.parse_args(argv)
    text = render(numbers())
    if a.check:
        if not OUT.is_file() or OUT.read_text() != text:
            print(f"{OUT.name} is stale. Rerun make_system_map.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(text)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
