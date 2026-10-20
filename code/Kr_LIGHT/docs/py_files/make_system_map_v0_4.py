"""Numbers for LIGHT_Kr v0_4's system map (system_map_v0_4.tex).

The map draws pure krypton at LIGHT's reference point (1 atm, 300 kW/cm3
for 1 us into the electrons, run to 2 us) under the default electron
closure, the kinetic EEDF, with every band at its default. The rules are
LIGHT's maps': every energy flow is a percent of the deposited energy, each
node is checked to conserve them before a number is written, and every row
of the network is drawn or printed. v0_4 resolves the 5s levels, so the
map draws the four as one node, the 5s manifold Kr*, whose energy is each
level's at its own energy: G1 makes the levels (no excess to heat), the
transfers between them (MIX_*) move energy between the manifold and the
gas, and every row's flow is its energy change in each node, summed over
the node's rows (dE). Reverses are drawn net of their forwards: the
excimer's formation net of its redissociation (AK8_* less AK8r_*), the
ion's conversion net of its (K15 less K15r). The 5s <-> 5p electron cycle
(K76_* up, K76r_* back) circulates many times the input, so it is drawn as
its net, as v0_3's map drew it: the electrons' net contribution, each
level's (K76 - K76r) x gap, enters Kr** through the K76 box, and the 5s
energy the rung hands back, the manifold's gain from KPS, PRAD, K76 and
K76r, runs from Kr** to Kr*. KPS leaves the gaps as heat, the largest
heat, so the map draws it as a band to the ledger; the other heats are
tags.

    python3 LIGHT_Kr/docs/make_system_map_v0_4.py          # the numbers
    cd LIGHT_Kr/docs && latexmk -pdf system_map_v0_4.tex   # the PDF

The runs come from kr_docs_v0_4.run_many (cached per version of the model's
code); a first run of the kinetic closure takes minutes. --check compares
the file with a fresh computation and writes nothing.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kr_docs_v0_4 as K                                   # noqa: E402
from v0_4 import jung2006                                  # noqa: E402
from v0_4.model import Model                               # noqa: E402
from v0_4.network import (BANDS, E_ASSOC_EV, KR_5P_TAU_S,  # noqa: E402
                          KR_5S, LEVELS_5S, OFF_BY_DEFAULT, build,
                          group_rows, kr5s)
from v0_4.species import e_form                            # noqa: E402

OUT = HERE / "system_map_v0_4_numbers.tex"
QE = 1.602176634e-19
CATIONS = ("Kr+", "Kr2+")
KR2 = ("Kr2*", "Kr2*(0u+)")
# the level groups (since v0_2): their rows' ids
G1_ROWS = tuple(r.rid for r in group_rows(build(), "G1"))
C1_ROWS = tuple(r.rid for g in ("C1a", "C1b", "C1c")
                for r in group_rows(build(), g))
HM_ROWS = tuple(r.rid for r in group_rows(build(), "HM"))
# v0_4: the rows by box
G2_ROWS = tuple(r.rid for r in build() if r.rid.startswith("G2_"))
K76F = tuple(f"K76_{lv}" for lv in LEVELS_5S)
K76R = tuple(f"K76r_{lv}" for lv in LEVELS_5S)
KPS_ROWS = tuple(f"KPS_{lv}" for lv in LEVELS_5S) + ("KPS_low",)
PRAD_ROWS = tuple(f"PRAD_{lv}" for lv in LEVELS_5S)
AK8F = ("AK8_1s5", "AK8_1s4", "AK8f_1s4", "AK8_1s3")
AK8R = ("AK8r_1s5", "AK8r_1s4", "AK8fr_1s4", "AK8r_1s3")
MIX_ROWS = ("MIX_1s4", "MIX_1s5r", "MIX_1s2", "MIX_1s4r2")
RES_ROWS = ("RES_1s4", "RES_1s2")
STEPS_ROWS = tuple(f"KSTEP_S_{lv}" for lv in LEVELS_5S)
RAD_ROWS = ("KR2_RAD", "KR2F_RAD")
# the rows 0 by default (no measurement, or off by default), which the map
# prints and draws no band for; v0_4 adds the fast excimer's rows, the
# quenching, the resonance lines and the 5p relaxation's 1s5-only share
ZERO_ROWS = (("PEN_KR2", "KSTEP_P", "AK8f_1s4", "AK8fr_1s4", "KR2F_RAD",
              "QKR2", "KPS_low") + STEPS_ROWS + RES_ROWS + HM_ROWS)
BAND_ORDER = ("k_5s5p", "k_5s5p_j1", "kps_kr", "kps_5s", "kr_pen",
              "kr2_form", "kr2_form_1s4", "kr2_form_1s3", "kr2_fast",
              "kr2_quench", "kr2_d0", "kr2_pen", "dr_tg", "res_escape",
              "kstep_s", "kstep_p", "hm_kr")


def _f1(x: float) -> str:
    """A flow as the map prints it: 0.1 resolution, no '-0.0'."""
    s = f"{x:.1f}"
    return "0.0" if s == "-0.0" else s


def _sci(x: float, d: int = 2) -> str:
    """x as m\\times10^{e} for math mode."""
    if x == 0.0:
        return "0"
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


def _k(x: float) -> str:
    """A rate constant as the map prints it (math)."""
    return _sci(x, 1)


def flows(s: dict, p) -> dict:
    """Every energy flow on the map as percent of E_dep (unrounded), from
    one run's summary (kr_docs_v0_4.summarize) and the model's rows."""
    m = Model(p)
    route = {c.row.rid: c.route for c in m.compiled}
    rows = {c.row.rid: c.row for c in m.compiled}
    ev = s["events"]
    e_dep = s["e_dep_eV"]
    E = {x: e_form(x) for x in KR_5S + ("Kr**", "Kr+", "Kr2+") + KR2}
    used = set()

    def de(rid, members):
        """The energy one event of the row moves into the node's species
        (members) [eV]: its products' there less its reactants'."""
        r = rows[rid]
        return (sum(E[x] for x in r.products if x in members)
                - sum(E[x] for x in r.reactants if x in members))
    # a row printed as 0 that carries energy (the solver's noise on a
    # species no row makes is far below 1e-9 of the input)
    def scale(r):
        return max([abs(v) for v in route[r].values()
                    if isinstance(v, float)]
                   + [E[x] for x in rows[r].reactants if x in E] + [1.0])
    loud = [r for r in ZERO_ROWS if abs(ev[r]) * scale(r) > 1e-9 * e_dep]
    if loud:
        raise AssertionError(f"rows the map prints as 0 carry energy: "
                             f"{loud}. Redraw them in system_map_v0_4.tex")

    def pct(x):
        return 100.0 * x / e_dep

    def ev_(*rids):
        used.update(rids)
        return sum(ev[r] for r in rids)

    def into(members, *rids):
        used.update(rids)
        return pct(sum(ev[r] * de(r, members) for r in rids))

    def heat(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["heat"] for r in rids))

    def rad(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["rad"] for r in rids))

    def e_from(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["from_e"] for r in rids))

    def e_to(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["to_e"] for r in rids))

    def left(*species):
        return pct(sum((s["n_end"][x] - s["n0"][x]) * E[x] for x in species))

    pp = ("Kr**",)
    f = {
        # the electron gas: in, the power and the Penning returns; out,
        # the electron-impact rows, the net 5s -> 5p pumping, the stepwise
        # rows, elastic loss, the captured electrons' energy and the store
        "e.P": pct(s["ledger"]["E_dep"]),
        "e.pen": e_to(*G2_ROWS),
        "e.pen2": e_to("PEN_KR2"),
        "e.G1": e_from(*G1_ROWS),
        "e.C1": e_from(*C1_ROWS),
        "e.HM": e_from(*HM_ROWS),
        "e.C2": e_from("C2"),
        "e.76": e_from(*K76F) - e_to(*K76R),
        "e.stepS": e_from(*STEPS_ROWS),
        "e.stepP": e_from("KSTEP_P"),
        "e.G6": pct(s["heat_elastic_eV"]),
        "e.dU": pct(s["ue_end"] - s["ue0"]),
        "e.cap": pct(s["ledger"]["heat_capture"]),
        # G1: the four 5s levels, each its own species (v0_4: no excess);
        # C1: the higher levels deliver Kr** and heat the gas with the
        # excess over the 2p10 energy; HM: the Hornbeck-Molnar rows
        # deliver Kr2+ (0 at the default band)
        "G1.krs": into(KR_5S, *G1_ROWS),
        "G1.heat": heat(*G1_ROWS),
        "C1.kr": into(pp, *C1_ROWS),
        "C1.heat": heat(*C1_ROWS),
        "HM.ion": into(CATIONS, *HM_ROWS),
        "HM.heat": heat(*HM_ROWS),
        # the Kr** reservoir: the 5s energy the rung hands back, net of the
        # electrons' pumping (the manifold's gain from these rows)
        "kpp.back": into(KR_5S, *(K76F + K76R + KPS_ROWS + PRAD_ROWS)),
        "KPS.heat": heat(*KPS_ROWS),
        "PRAD.hv": rad(*PRAD_ROWS),
        "PRAD.heat": heat(*PRAD_ROWS),
        "kpp.stepP": -into(pp, "KSTEP_P"),
        "kpp.left": left("Kr**"),
        # the 5s manifold, fed by G1 and the rung's net return; out to
        # Penning, to the excimer (net of its redissociation), the
        # stepwise ionization, the resonance lines and, through the
        # transfers between the levels, the gas
        "krs.G2": -into(KR_5S, *G2_ROWS),
        "krs.AK8": -into(KR_5S, *(AK8F + AK8R)),
        "krs.AK8f": -into(KR_5S, *AK8F),
        "krs.AK8r": into(KR_5S, *AK8R),
        "krs.stepS": -into(KR_5S, *STEPS_ROWS),
        "krs.res": -into(KR_5S, *RES_ROWS),
        "krs.mix": -into(KR_5S, *MIX_ROWS),
        "MIX.heat": heat(*MIX_ROWS),
        "RES.hv": rad(*RES_ROWS),
        "krs.left": left(*KR_5S),
        # Penning and stepwise ionization into the ions
        "G2.ion": into(CATIONS, *G2_ROWS),
        "C2.ion": into(CATIONS, "C2"),
        "stepS.ion": into(CATIONS, *STEPS_ROWS),
        "stepP.ion": into(CATIONS, "KSTEP_P"),
        # the Kr2* node (both states, the formation net of the
        # redissociation inside it), its light, its quenching, the row at 0
        "AK8.heat": heat(*(AK8F + AK8R)),
        "k2s.RAD": -into(KR2, *RAD_ROWS),
        "k2s.Q": -into(KR2, "QKR2"),
        "k2s.pen": -into(KR2, "PEN_KR2"),
        "k2s.left": left(*KR2),
        "RAD.hv": rad(*RAD_ROWS),
        "RAD.heat": heat(*RAD_ROWS),
        "Q.heat": heat("QKR2"),
        "pen2.ion": into(CATIONS, "PEN_KR2"),
        # the ions: Kr+ -> Kr2+ net of its redissociation (K15, K15r),
        # Kr2+ + e -> Kr** (R11p)
        "K15.heat": heat("K15", "K15r"),
        "R11p.in": -into(CATIONS, "R11p"),
        "R11p.kr": into(pp, "R11p"),
        "R11p.heat": heat("R11p"),
        "ions.left": left(*CATIONS),
        # the fate at t_end
        "out.147": pct(s["bands_eV"]["147 nm (Kr2* VUV)"]),
        "out.IR": pct(s["bands_eV"]["5p-5s (near IR)"]),
        "out.res": pct(s["bands_eV"]["123.6 nm (resonance)"]
                       + s["bands_eV"]["116.5 nm (resonance)"]),
        "out.heat": pct(s["heat_reactions_eV"]),
        "out.elastic": pct(s["heat_elastic_eV"]),
        "out.ions": pct(s["stores_eV"]["ions (Kr+, Kr2+)"]),
        "out.other": pct(s["stores_eV"]["Kr*"] + s["stores_eV"]["Kr** (5p)"]
                         + s["stores_eV"]["Kr2*"]
                         + s["stores_eV"]["electron gas"]),
    }
    used.add("G6")                         # G6 is the elastic loss, e.G6
    missing = set(rows) - used
    if missing:
        raise AssertionError(f"rows not drawn on the map: {sorted(missing)}")
    if f["e.76"] <= 0.0 or f["kpp.back"] <= 0.0:
        raise AssertionError("the 5s <-> 5p cycle runs the other way: "
                             "redraw its two net bands")
    if f["krs.AK8"] <= 0.0:
        raise AssertionError("the excimer redissociates more than it forms: "
                             "redraw the formation band")
    if f["ions.left"] > 0.0:
        raise AssertionError("the ions gain energy: the map's ledger note "
                             "says the seed's ions recombine")
    if set(s["bands_eV"]) != {"147 nm (Kr2* VUV)", "5p-5s (near IR)",
                              "123.6 nm (resonance)", "116.5 nm (resonance)"}:
        raise AssertionError(f"bands changed: {sorted(s['bands_eV'])}")
    if s["ledger"]["E_wall"] != 0.0 or s["ledger"]["E_em"] != 0.0:
        raise AssertionError("the map's run has walls or a wave")
    _check(f, s["closure"])
    return f


def _check(f: dict, closure: float, tol: float = 1e-8) -> None:
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
                         f["e.G1"] + f["e.C1"] + f["e.HM"] + f["e.C2"]
                         + f["e.76"] + f["e.stepS"] + f["e.stepP"]
                         + f["e.G6"] + f["e.dU"] + f["e.cap"]),
        "G1": (f["e.G1"], f["G1.krs"] + f["G1.heat"]),
        "C1": (f["e.C1"], f["C1.kr"] + f["C1.heat"]),
        "HM": (f["e.HM"], f["HM.ion"] + f["HM.heat"]),
        "Kr**": (f["C1.kr"] + f["R11p.kr"] + f["e.76"],
                 f["kpp.back"] + f["KPS.heat"] + f["PRAD.hv"]
                 + f["PRAD.heat"] + f["kpp.stepP"] + f["kpp.left"]),
        "KSTEP_P": (f["kpp.stepP"] + f["e.stepP"], f["stepP.ion"]),
        "Kr* (the 5s levels)": (
            f["G1.krs"] + f["kpp.back"],
            f["krs.G2"] + f["krs.AK8"] + f["krs.stepS"] + f["krs.res"]
            + f["krs.mix"] + f["krs.left"]),
        "the formation, net": (f["krs.AK8f"] - f["krs.AK8r"],
                               f["krs.AK8"]),
        "MIX": (f["krs.mix"], f["MIX.heat"]),
        "RES": (f["krs.res"], f["RES.hv"]),
        "KSTEP_S": (f["krs.stepS"] + f["e.stepS"], f["stepS.ion"]),
        "G2": (f["krs.G2"], f["G2.ion"] + f["e.pen"]),
        "Kr2* (with AK8, net)": (f["krs.AK8"],
                                 f["AK8.heat"] + f["k2s.RAD"] + f["k2s.Q"]
                                 + f["k2s.pen"] + f["k2s.left"]),
        "KR2_RAD": (f["k2s.RAD"], f["RAD.hv"] + f["RAD.heat"]),
        "QKR2": (f["k2s.Q"], f["Q.heat"]),
        "Kr2* + Kr2*": (f["k2s.pen"], f["pen2.ion"] + f["e.pen2"]),
        "ions": (f["C2.ion"] + f["G2.ion"] + f["pen2.ion"] + f["stepS.ion"]
                 + f["stepP.ion"] + f["HM.ion"],
                 f["K15.heat"] + f["R11p.in"] + f["ions.left"]),
        "R11p": (f["R11p.in"], f["R11p.kr"] + f["R11p.heat"]),
        "heat": (f["G1.heat"] + f["C1.heat"] + f["HM.heat"]
                 + f["KPS.heat"] + f["PRAD.heat"] + f["MIX.heat"]
                 + f["AK8.heat"] + f["RAD.heat"] + f["Q.heat"]
                 + f["K15.heat"] + f["R11p.heat"] + f["e.cap"],
                 f["out.heat"]),
        "147 nm": (f["RAD.hv"], f["out.147"]),
        "near IR": (f["PRAD.hv"], f["out.IR"]),
        "resonance lines": (f["RES.hv"], f["out.res"]),
        "stored ions": (f["ions.left"], f["out.ions"]),
        "other stores": (f["krs.left"] + f["kpp.left"] + f["k2s.left"]
                         + f["e.dU"], f["out.other"]),
        "ledger": (100.0, f["out.147"] + f["out.IR"] + f["out.res"]
                   + f["out.heat"] + f["out.elastic"] + f["out.ions"]
                   + f["out.other"]),
    }


def _rate(m, rid) -> float:
    c = next(c for c in m.compiled if c.row.rid == rid)
    if c.kind not in ("const", "param"):
        raise ValueError(f"{rid}: rate kind {c.kind!r} has no single value")
    return float(c.rate) * c.factor


def _band_value(b, v) -> str:
    if b.unit in ("factor", "fraction", "eV"):
        return f"{v:g}"
    return "0" if v == 0.0 else _k(v)


# the bands in LaTeX (the documents' macros)
BAND_LABELS = {
    "kr2_pen": r"\Krts{} + \Krts{} \ra \Krtp{} + 2Kr + \el{}",
    "kr2_form": r"\Krs(1s$_5$) + 2Kr \ra \Krts{} + Kr",
    "kr2_form_1s4": r"\Krs(1s$_4$) + 2Kr \ra \Krts{} + Kr",
    "kr2_form_1s3": r"\Krs(1s$_3$) + 2Kr \ra \Krts{} + Kr",
    "kr2_fast": r"the share of the 1s$_4$ formation into the fast 0$_u^+$ "
                r"state",
    "kr2_quench": r"\Krts{} + Kr \ra 3Kr",
    "kr2_d0": r"the excimer's well, for its redissociation",
    "dr_tg": r"the gas-temperature factor on \Krtp{} + \el{}",
    "res_escape": r"the escaping resonance lines, per 1s$_4$ and 1s$_2$",
    "kps_5s": r"the share of \Krpp{} + Kr's products over the 5s levels",
    "kr_pen": r"\Krs{} + \Krs{} \ra \Krp{} + Kr + \el{}",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "k_5s5p": r"\el{} + \Krs{} \ra \el{} + \Krpp{}, a scale on the "
              r"measured cross sections",
    "k_5s5p_j1": r"the same, a scale on the $J=1$ levels' proxies",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}, pre-exponential",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}, pre-exponential",
    "hm_kr": r"Kr(nl) + Kr \ra \Krtp{} + \el{}, the share of the levels' "
             r"excitations",
}
UNITS = {"cm3/s": r"\unit{cm^3\,s^{-1}}", "cm6/s": r"\unit{cm^6\,s^{-1}}",
         "factor": "", "fraction": "", "1/s": r"\unit{s^{-1}}",
         "eV": r"\unit{eV}"}


def bands_numbers(ref: dict, ref2: dict) -> tuple:
    """The 147-nm share at each end of every band the rows read, one band
    at a time, on the kinetic EEDF (the map's box) and the two-temperature
    closure, and the overview's table of both (\\wmBands)."""
    jobs = {}
    for eedf in K.CLOSURES:
        for name, end, v, p in K.band_params(eedf):
            jobs[(eedf, name, end)] = ("field", p)
    keys = list(jobs)
    got = dict(zip(keys, K.run_many([jobs[k] for k in keys])))
    names = [b.name for b in K.row_bands()]
    if set(names) != set(BAND_ORDER):
        raise AssertionError(f"the bands the rows read changed: {names}")
    n, lines = {}, []
    for name in BAND_ORDER:
        b = BANDS[name]
        ends = {}
        for eedf, r0 in (("auto", ref), ("two_t", ref2)):
            for end, v in (("low", b.lo), ("high", b.hi)):
                d = r0 if v == b.default else got[(eedf, name, end)]
                ends[(eedf, end)] = d
        lo, hi = ends[("auto", "low")], ends[("auto", "high")]
        n[f"bd.{name}.def"] = _band_value(b, b.default)
        n[f"bd.{name}.lo"] = _band_value(b, b.lo)
        n[f"bd.{name}.hi"] = _band_value(b, b.hi)
        n[f"bd.{name}.lo147"] = f"{lo['emission_vuv']:.2f}"
        n[f"bd.{name}.hi147"] = f"{hi['emission_vuv']:.2f}"
        lo2, hi2 = ends[("two_t", "low")], ends[("two_t", "high")]
        lines.append(
            f"{BAND_LABELS[name]} & ${_band_value(b, b.default)}$ & "
            f"${_band_value(b, b.lo)}$ to ${_band_value(b, b.hi)}$"
            f"{UNITS[b.unit]} & "
            f"{lo['emission_vuv']:.2f} to {hi['emission_vuv']:.2f} & "
            f"{lo2['emission_vuv']:.2f} to {hi2['emission_vuv']:.2f} & "
            f"{lo['emission_nir']:.2f} to {hi['emission_nir']:.2f} & "
            f"{lo['tg_end_K']:.0f} to {hi['tg_end_K']:.0f}\\\\")
    n["bd.n"] = str(len(BAND_ORDER))
    # v0_4: the 0 defaults with no measurement, and those off by default
    # though a measurement bounds them (network.OFF_BY_DEFAULT)
    n["bd.nzero"] = str(sum(1 for x in BAND_ORDER if BANDS[x].default == 0.0
                            and x not in OFF_BY_DEFAULT))
    n["bd.noff"] = str(sum(1 for x in BAND_ORDER if x in OFF_BY_DEFAULT))
    return n, "\\def\\wmBands{" + "\n".join(lines) + "}"


def numbers() -> tuple:
    """(numbers, tables): the strings the map prints, keyed as in the
    .tex, and the overview's bands table."""
    p = K.ref_params()
    p2 = K.ref_params("two_t")
    sw = K.switch_params()
    ref, ref2, kann, nocyc = K.run_many(
        [("field", p), ("field", p2), ("field", sw["kannari"]),
         ("field", sw["nocycle"])])
    f = flows(ref, p)
    n = {k: _f1(v) for k, v in f.items()}
    m = Model(p)
    rows = {c.row.rid: c.row for c in m.compiled}
    ev = ref["events"]
    kin = ref["kin"]
    krs, kpp, k2s = ref["res"]["Kr*"], ref["res"]["Kr**"], ref["res"]["Kr2*"]
    kind, k0, te0, n_pow = ("te_pow",) + tuple(rows["R11p"].rate[1:])

    def br(rv, *rids):
        return f"{100.0 * sum(rv['losses'].get(r, 0.0) for r in rids) / rv['formed']:.1f}"

    def gain(rv, *rids):
        """A share of a reservoir's formation, percent: the rows that make
        it (reservoir() books them as formation)."""
        tot = 0.0
        for r in rids:
            row = rows[r]
            net = (sum(row.products.count(x) for x in RES_MEMBERS[rv])
                   - sum(row.reactants.count(x) for x in RES_MEMBERS[rv]))
            tot += max(net, 0) * ev[r]
        return f"{100.0 * tot / ref['res'][rv]['formed']:.1f}"

    k_ion = sum(ev[r] for r in ("C2", "PEN_KR2", "KSTEP_P") + G2_ROWS
                + STEPS_ROWS + HM_ROWS)
    form = sum(ev[r] for r in AK8F)
    redis = sum(ev[r] for r in AK8R)
    # v0_3: the trained factors the cross sections carry (xsec.trained())
    fac = K.xsec.trained()
    n.update({"fac.x1s": f"{fac.exc[0]:.2f}", "fac.x2p": f"{fac.exc[1]:.2f}",
              "fac.xhi": f"{fac.exc[2]:.2f}", "fac.mlo": f"{min(fac.m):.2f}",
              "fac.mhi": f"{max(fac.m):.2f}"})
    n.update({
        "version": K.V,
        "md5": K.full_md5()[:12],
        "p": f"{p.pressure_atm:g}",
        "Tg0": f"{p.tg0_K:.0f}",
        "P": f"{p.power_W_cm3 / 1e3:.0f}",
        "pulse": f"{p.pulse_s * 1e6:g}",
        "tend": f"{p.t_end_s * 1e6:g}",
        "ne": _sci(p.ne_cm3, 0).replace("1\\times", ""),
        "N": _sci(ref["n0"]["Kr"]),
        "Edep": f"{ref['ledger']['E_dep'] * QE:.2f}",
        "Te": f"{ref['te_mean_eV']:.2f}",
        "TeMax": f"{ref['te_max']:.2f}",
        "enLo": f"{kin['en_min']:.2g}", "enHi": f"{kin['en_max']:.3g}",
        "kinN": str(kin["solutions"]), "kinOn": str(kin["on"]),
        "kinRatioLo": f"{kin['ratio_min']:.4f}",
        "kinRatioHi": f"{kin['ratio_max']:.4f}",
        "kinGap": _sci(kin["balance_max"], 0).replace("1\\times", ""),
        "closure": _sci(abs(ref["closure"]), 0).replace("1\\times", ""),
        "TgEnd": f"{ref['tg_end_K']:.0f}",
        "TgMax": f"{ref['tg_max']:.0f}",
        "pEnd": f"{ref['p_end_atm']:.2f}",
        "nePeak": _sci(ref["ne_peak"]),
        "neEnd": _sci(ref["ne_end_cm3"]),
        "krsFormed": _sci(krs["formed"]),
        "kppFormed": _sci(kpp["formed"]),
        "k2sFormed": _sci(k2s["formed"]),
        # where Kr*, Kr** and Kr2* go, percent of each formed (gross)
        "brs.76": br(krs, *K76F), "brs.AK8": br(krs, *AK8F),
        "brs.G2": br(krs, *G2_ROWS),
        "brp.KPS": br(kpp, *KPS_ROWS), "brp.76r": br(kpp, *K76R),
        "brp.PRAD": br(kpp, *PRAD_ROWS),
        "brk2.RAD": br(k2s, *RAD_ROWS), "brk2.redis": br(k2s, *AK8R),
        "k2s.leftpc": f"{100.0 * k2s['left'] / k2s['formed']:.1f}",
        # v0_4: what makes the 5s manifold, percent of its formation
        "gs.G1": gain("Kr*", *G1_ROWS), "gs.KPS": gain("Kr*", *KPS_ROWS),
        "gs.PRAD": gain("Kr*", *PRAD_ROWS), "gs.76r": gain("Kr*", *K76R),
        "gs.redis": gain("Kr*", *AK8R),
        # the excimer's formation by level, percent of the gross formation
        **{f"form.{lv}": f"{100.0 * sum(ev[r] for r in AK8F if r.endswith(lv)) / form:.1f}"
           for lv in ("1s5", "1s4", "1s3")},
        "redis.pc": f"{100.0 * redis / form:.1f}",
        # the cycle: K76 events per Kr2* formed net
        "cyc.perK2": f"{sum(ev[r] for r in K76F) / (form - redis):.0f}",
        "cyc.e": _f1(f["e.76"]),
        "cyc.net": _sci(sum(ev[r] for r in K76F) - sum(ev[r] for r in K76R)),
        # the ions: who makes them
        "ion.G2": f"{100.0 * sum(ev[r] for r in G2_ROWS) / k_ion:.0f}"
                  if k_ion else "0",
        "ion.C2": f"{100.0 * ev['C2'] / k_ion:.2g}" if k_ion else "0",
        # the constants the boxes print, read from the rows the model runs
        "k.G2": _k(_rate(m, "G2_1s5_1s5")),
        "k.AK8": _k(_rate(m, "AK8_1s5")),
        "k.AK8.1s4": _k(_rate(m, "AK8_1s4")),
        "tpow.AK8": f"{rows['AK8_1s5'].tg_pow:.2f}",
        "tauRAD": f"{1e9 / _rate(m, 'KR2_RAD'):.0f}",
        "tauRADF": f"{1e9 / _rate(m, 'KR2F_RAD'):.1f}",
        "hv.RAD": f"{rows['KR2_RAD'].hv_eV:.2f}",
        "E.kr2s": f"{e_form('Kr2*'):.2f}",
        "E.kr2p": f"{e_form('Kr2+'):.3f}",
        "k.K15": _k(_rate(m, "K15")),
        "tpow.K15": f"{rows['K15'].tg_pow:.2f}",
        "k.R11p": _k(k0),
        "pow.R11p": f"{-n_pow:.2f}",
        "k.KPS": _k(sum(_rate(m, r) for r in KPS_ROWS)),
        "k.MIX": _k(m._k(next(c for c in m.compiled
                               if c.row.rid == "MIX_1s4"), 1.0, 300.0,
                          None, None)),
        "tauPRAD": f"{KR_5P_TAU_S * 1e9:.1f}",
        "E.PP": f"{e_form('Kr**'):.2f}",
        "E.S": f"{e_form(kr5s('1s5')):.4f}",
        **{f"E.{lv}": f"{e_form(kr5s(lv)):.4f}" for lv in LEVELS_5S},
        "E.ion": f"{e_form('Kr+'):.4f}",
        "E.ai": f"{E_ASSOC_EV:.2f}",
        "nG1": str(len(G1_ROWS)), "nC1": str(len(C1_ROWS)),
        "nHM": str(len(HM_ROWS)),
        "C1.lo": f"{min(rows[r].e_cost_eV for r in C1_ROWS):.2f}",
        "C1.hi": f"{max(rows[r].e_cost_eV for r in C1_ROWS):.1f}",
        "G1.hi": f"{max(rows[r].e_cost_eV for r in G1_ROWS):.2f}",
        "gap5p": f"{e_form('Kr**') - e_form(kr5s('1s5')):.3f}",
        "gap5p.lo": f"{e_form('Kr**') - e_form(kr5s('1s2')):.3f}",
        "pen.e": f"{2 * e_form(kr5s('1s5')) - e_form('Kr+'):.3f}",
        "sig.76": _k(float(jung2006.lumped_sigma(6.0))),
        "nodes": str(len(nodes(f))),
        "nrows": str(len(rows)),
        # the notes' output line (flows, they overlap)
        "notes.krs": f"{ref['to_Krstar']:.1f}",
        "notes.krp": f"{ref['to_Krplus']:.1f}",
        "notes.k2s": f"{ref['to_Kr2star']:.1f}",
        # the comparisons the box prints
        "om.kannari": f"{kann['emission_vuv']:.2f}",
        "om.nocycle": f"{nocyc['emission_vuv']:.2f}",
        "om.twot": f"{ref2['emission_vuv']:.2f}",
        "om.twot.te": f"{ref2['te_mean_eV']:.2f}",
    })
    bn, table = bands_numbers(ref, ref2)
    n.update(bn)
    return n, [table]


RES_MEMBERS = {"Kr*": KR_5S, "Kr**": ("Kr**",), "Kr2*": KR2}


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_4 map numbers. GENERATED by make_system_map_v0_4.py",
        f"% from the model ({n['version']}). Do not edit. Rerun the script "
        "instead.",
        "\\providecommand\\wmdef[2]{\\expandafter\\def"
        "\\csname wm@#1\\endcsname{#2}}",
    ]
    lines += [f"\\wmdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_4's system map "
                                 "numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_system_map_v0_4.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
