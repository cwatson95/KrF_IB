"""Numbers for LIGHT_Kr v0_2's system map (system_map_v0_2.tex).

The map draws pure krypton at LIGHT's reference point (1 atm, 300 kW/cm3
for 1 us into the electrons, run to 2 us) under the default electron
closure, the kinetic EEDF, with every band at its default. The rules are
LIGHT's maps': every energy flow is a percent of the deposited energy, each
node is checked to conserve them before a number is written, and every row
of the network is drawn or printed. v0_2 books the ground-state
excitations level by level, so the map draws their groups: G1 (the four
5s levels, into Kr*, the excess over the 1s5 energy to heat), C1 (every
higher level, into the 5p rung Kr**, the excess over the 2p10 energy to
heat) and HM (the Hornbeck-Molnar rows, 0 at the default band, printed).
The 5s <-> 5p electron cycle (K76 up, K76r back) circulates many times the
input, so it is drawn as its net, as LIGHT v0_4's map drew it: the
electrons' net contribution, (K76 - K76r) x E(Kr**) - E(Kr*), enters Kr**
through the K76 box, and the Kr* energy the rung hands back, (KPS + PRAD
+ K76r - K76) x E(Kr*), runs from Kr** to Kr*. KPS leaves the gap as heat,
three quarters of the input here, so the map draws that heat as a band to
the ledger; the other heats are tags.

    python3 LIGHT_Kr/docs/make_system_map_v0_2.py          # the numbers
    cd LIGHT_Kr/docs && latexmk -pdf system_map_v0_2.tex   # the PDF

The runs come from kr_docs_v0_2.run_many (cached per version of the model's
code); a first run of the kinetic closure takes minutes. --check compares
the file with a fresh computation and writes nothing.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kr_docs_v0_2 as K                                   # noqa: E402
from v0_2 import jung2006                                  # noqa: E402
from v0_2.model import Model                               # noqa: E402
from v0_2.network import (BANDS, E_ASSOC_EV, KR_5P_TAU_S,  # noqa: E402
                          build, group_rows)
from v0_2.species import e_form                            # noqa: E402

OUT = HERE / "system_map_v0_2_numbers.tex"
QE = 1.602176634e-19
CATIONS = ("Kr+", "Kr2+")
# the level groups (v0_2): their rows' ids
G1_ROWS = tuple(r.rid for r in group_rows(build(), "G1"))
C1_ROWS = tuple(r.rid for g in ("C1a", "C1b", "C1c")
                for r in group_rows(build(), g))
HM_ROWS = tuple(r.rid for r in group_rows(build(), "HM"))
# the rows with no measurement, 0 unless a band run sets them: the map
# prints them and draws no band
ZERO_ROWS = ("PEN_KR2", "KSTEP_S", "KSTEP_P") + HM_ROWS
BAND_ORDER = ("k_5s5p", "kps_kr", "kr_pen", "kr2_form", "kr2_pen",
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
    one run's summary (kr_docs.summarize) and the model's routes."""
    m = Model(p)
    route = {c.row.rid: c.route for c in m.compiled}
    rows = {c.row.rid: c.row for c in m.compiled}
    ev = s["events"]
    e_dep = s["e_dep_eV"]
    E = {x: e_form(x) for x in ("Kr*", "Kr**", "Kr+", "Kr2+", "Kr2*")}
    used = set()
    nonzero = [r for r in ZERO_ROWS if ev[r] != 0.0]
    if nonzero:
        raise AssertionError(f"rows the map prints as 0 carry energy: "
                             f"{nonzero}. Redraw them in system_map_v0_2.tex")

    def pct(x):
        return 100.0 * x / e_dep

    def ev_(*rids):
        used.update(rids)
        return sum(ev[r] for r in rids)

    def heat(*rids):
        used.update(rids)
        return pct(sum(ev[r] * route[r]["heat"] for r in rids))

    def rad(rid):
        used.add(rid)
        return pct(ev[rid] * route[rid]["rad"])

    def left(*species):
        return pct(sum((s["n_end"][x] - s["n0"][x]) * E[x] for x in species))

    gap = E["Kr**"] - E["Kr*"]
    net76 = ev_("K76") - ev_("K76r")
    back = (ev_("KPS") + ev_("PRAD") - net76) * E["Kr*"]
    f = {
        # the electron gas: in, the power and the Penning returns; out,
        # the electron-impact rows, the net 5s -> 5p pumping, the stepwise
        # rows, elastic loss, the captured electrons' energy and the store
        "e.P": pct(s["ledger"]["E_dep"]),
        "e.pen": pct(ev_("G2") * route["G2"]["to_e"]),
        "e.pen2": pct(ev_("PEN_KR2") * route["PEN_KR2"]["to_e"]),
        "e.G1": pct(sum(ev_(r) * route[r]["from_e"] for r in G1_ROWS)),
        "e.C1": pct(sum(ev_(r) * route[r]["from_e"] for r in C1_ROWS)),
        "e.HM": pct(sum(ev_(r) * route[r]["from_e"] for r in HM_ROWS)),
        "e.C2": pct(ev_("C2") * route["C2"]["from_e"]),
        "e.76": pct(net76 * gap),
        "e.stepS": pct(ev_("KSTEP_S") * route["KSTEP_S"]["from_e"]),
        "e.stepP": pct(ev_("KSTEP_P") * route["KSTEP_P"]["from_e"]),
        "e.G6": pct(s["heat_elastic_eV"]),
        "e.dU": pct(s["ue_end"] - s["ue0"]),
        "e.cap": pct(s["ledger"]["heat_capture"]),
        # G1: the four 5s levels deliver Kr* and heat the gas with the
        # excess over the 1s5 energy; C1: the higher levels deliver Kr**
        # and heat the gas with the excess over the 2p10 energy; HM: the
        # Hornbeck-Molnar rows deliver Kr2+ (0 at the default band)
        "G1.krs": pct(ev_(*G1_ROWS) * E["Kr*"]),
        "G1.heat": heat(*G1_ROWS),
        "C1.kr": pct(ev_(*C1_ROWS) * E["Kr**"]),
        "C1.heat": heat(*C1_ROWS),
        "HM.ion": pct(ev_(*HM_ROWS) * E["Kr2+"]),
        "HM.heat": heat(*HM_ROWS),
        # the Kr** reservoir
        "kpp.back": pct(back),
        "KPS.heat": heat("KPS"),
        "PRAD.hv": rad("PRAD"),
        "PRAD.heat": heat("PRAD"),
        "kpp.stepP": pct(ev_("KSTEP_P") * E["Kr**"]),
        "kpp.left": left("Kr**"),
        # the Kr* reservoir, fed by G1 and by the rung's net return
        "krs.G2": pct(2 * ev_("G2") * E["Kr*"]),
        "krs.AK8": pct(ev_("AK8_Kr") * E["Kr*"]),
        "krs.stepS": pct(ev_("KSTEP_S") * E["Kr*"]),
        "krs.left": left("Kr*"),
        # Penning and stepwise ionization into the ions
        "G2.ion": pct(ev["G2"] * E["Kr+"]),
        "C2.ion": pct(ev["C2"] * E["Kr+"]),
        "stepS.ion": pct(ev["KSTEP_S"] * E["Kr+"]),
        "stepP.ion": pct(ev["KSTEP_P"] * E["Kr+"]),
        # the Kr2* node (AK8_Kr inside it), its light and the row at 0
        "AK8.heat": heat("AK8_Kr"),
        "k2s.RAD": pct(ev_("KR2_RAD") * E["Kr2*"]),
        "k2s.pen": pct(2 * ev["PEN_KR2"] * E["Kr2*"]),
        "k2s.left": left("Kr2*"),
        "RAD.hv": rad("KR2_RAD"),
        "RAD.heat": heat("KR2_RAD"),
        "pen2.ion": pct(ev["PEN_KR2"] * E["Kr2+"]),
        # the ions: Kr+ -> Kr2+ (K15), Kr2+ + e -> Kr** (R11p)
        "K15.heat": heat("K15"),
        "R11p.in": pct(ev_("R11p") * E["Kr2+"]),
        "R11p.kr": pct(ev["R11p"] * E["Kr**"]),
        "R11p.heat": heat("R11p"),
        "ions.left": left(*CATIONS),
        # the fate at t_end
        "out.147": pct(s["bands_eV"]["147 nm (Kr2* VUV)"]),
        "out.IR": pct(s["bands_eV"]["5p-5s (near IR)"]),
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
    if net76 <= 0.0 or back <= 0.0:
        raise AssertionError("the 5s <-> 5p cycle runs the other way: "
                             "redraw its two net bands")
    if f["ions.left"] > 0.0:
        raise AssertionError("the ions gain energy: the map's ledger note "
                             "says the seed's ions recombine")
    if set(s["bands_eV"]) != {"147 nm (Kr2* VUV)", "5p-5s (near IR)"}:
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
        "Kr*": (f["G1.krs"] + f["kpp.back"],
                f["krs.G2"] + f["krs.AK8"] + f["krs.stepS"] + f["krs.left"]),
        "KSTEP_S": (f["krs.stepS"] + f["e.stepS"], f["stepS.ion"]),
        "G2": (f["krs.G2"], f["G2.ion"] + f["e.pen"]),
        "Kr2* (with AK8_Kr)": (f["krs.AK8"], f["AK8.heat"] + f["k2s.RAD"]
                               + f["k2s.pen"] + f["k2s.left"]),
        "KR2_RAD": (f["k2s.RAD"], f["RAD.hv"] + f["RAD.heat"]),
        "Kr2* + Kr2*": (f["k2s.pen"], f["pen2.ion"] + f["e.pen2"]),
        "ions": (f["C2.ion"] + f["G2.ion"] + f["pen2.ion"] + f["stepS.ion"]
                 + f["stepP.ion"] + f["HM.ion"],
                 f["K15.heat"] + f["R11p.in"] + f["ions.left"]),
        "R11p": (f["R11p.in"], f["R11p.kr"] + f["R11p.heat"]),
        "heat": (f["G1.heat"] + f["C1.heat"] + f["HM.heat"]
                 + f["KPS.heat"] + f["PRAD.heat"]
                 + f["AK8.heat"] + f["RAD.heat"] + f["K15.heat"]
                 + f["R11p.heat"] + f["e.cap"], f["out.heat"]),
        "147 nm": (f["RAD.hv"], f["out.147"]),
        "near IR": (f["PRAD.hv"], f["out.IR"]),
        "stored ions": (f["ions.left"], f["out.ions"]),
        "other stores": (f["krs.left"] + f["kpp.left"] + f["k2s.left"]
                         + f["e.dU"], f["out.other"]),
        "ledger": (100.0, f["out.147"] + f["out.IR"] + f["out.heat"]
                   + f["out.elastic"] + f["out.ions"] + f["out.other"]),
    }


def _rate(m, rid) -> float:
    c = next(c for c in m.compiled if c.row.rid == rid)
    if c.kind not in ("const", "param"):
        raise ValueError(f"{rid}: rate kind {c.kind!r} has no single value")
    return float(c.rate) * c.factor


def _band_value(b, v) -> str:
    if b.unit in ("factor", "fraction"):
        return f"{v:g}"
    return "0" if v == 0.0 else _k(v)


# the bands in LaTeX (the documents' macros)
BAND_LABELS = {
    "kr2_pen": r"\Krts{} + \Krts{} \ra \Krtp{} + 2Kr + \el{}",
    "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{} \ra \Krp{} + Kr + \el{}",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "k_5s5p": r"\el{} + \Krs{} \ra \el{} + \Krpp{}, a scale on the "
              r"measured cross section",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}, pre-exponential",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}, pre-exponential",
    "hm_kr": r"Kr(nl) + Kr \ra \Krtp{} + \el{}, the share of the levels' "
             r"excitations",
}
UNITS = {"cm3/s": r"\unit{cm^3\,s^{-1}}", "cm6/s": r"\unit{cm^6\,s^{-1}}",
         "factor": "", "fraction": ""}


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
    n["bd.nzero"] = str(sum(1 for x in BAND_ORDER
                            if BANDS[x].default == 0.0))
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

    k_ion = sum(ev[r] for r in ("C2", "G2", "KSTEP_S", "KSTEP_P", "PEN_KR2")
                + HM_ROWS)
    n.update({
        "version": K.V,
        "md5": K.code_md5()[:12],
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
        # where Kr*, Kr** and Kr2* go, percent of each formed
        "brs.76": br(krs, "K76"), "brs.AK8": br(krs, "AK8_Kr"),
        "brs.G2": br(krs, "G2"),
        "brp.KPS": br(kpp, "KPS"), "brp.76r": br(kpp, "K76r"),
        "brp.PRAD": br(kpp, "PRAD"),
        "brk2.RAD": br(k2s, "KR2_RAD"),
        "k2s.leftpc": f"{100.0 * k2s['left'] / k2s['formed']:.1f}",
        # the cycle: K76 events per Kr2* formed, and per Kr* made by the
        # electrons (G1) or handed down by the rung's other sources
        "cyc.perK2": f"{ev['K76'] / ev['AK8_Kr']:.0f}",
        "cyc.e": _f1(f["e.76"]),
        "cyc.net": _sci(ev["K76"] - ev["K76r"]),
        # the ions: who makes them
        "ion.G2": f"{100.0 * ev['G2'] / k_ion:.0f}" if k_ion else "0",
        "ion.C2": f"{100.0 * ev['C2'] / k_ion:.2g}" if k_ion else "0",
        # the constants the boxes print, read from the rows the model runs
        "k.G2": _k(_rate(m, "G2")),
        "k.AK8": _k(_rate(m, "AK8_Kr")),
        "tpow.AK8": f"{rows['AK8_Kr'].tg_pow:.2f}",
        "tauRAD": f"{1e9 / _rate(m, 'KR2_RAD'):.0f}",
        "hv.RAD": f"{rows['KR2_RAD'].hv_eV:.2f}",
        "E.kr2s": f"{e_form('Kr2*'):.2f}",
        "E.kr2p": f"{e_form('Kr2+'):.3f}",
        "k.K15": _k(_rate(m, "K15")),
        "tpow.K15": f"{rows['K15'].tg_pow:.2f}",
        "k.R11p": _k(k0),
        "pow.R11p": f"{-n_pow:.2f}",
        "k.KPS": _k(_rate(m, "KPS")),
        "tauPRAD": f"{KR_5P_TAU_S * 1e9:.1f}",
        "E.PP": f"{e_form('Kr**'):.2f}",
        "E.S": f"{e_form('Kr*'):.4f}",
        "E.ion": f"{e_form('Kr+'):.4f}",
        "E.ai": f"{E_ASSOC_EV:.2f}",
        "nG1": str(len(G1_ROWS)), "nC1": str(len(C1_ROWS)),
        "nHM": str(len(HM_ROWS)),
        "C1.lo": f"{min(rows[r].e_cost_eV for r in C1_ROWS):.2f}",
        "C1.hi": f"{max(rows[r].e_cost_eV for r in C1_ROWS):.1f}",
        "G1.hi": f"{max(rows[r].e_cost_eV for r in G1_ROWS):.2f}",
        "gap5p": f"{e_form('Kr**') - e_form('Kr*'):.3f}",
        "pen.e": f"{2 * e_form('Kr*') - e_form('Kr+'):.3f}",
        "sig.76": _k(float(jung2006.lumped_sigma(6.0))),
        "nodes": str(len(nodes(f))),
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


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_2 map numbers. GENERATED by make_system_map_v0_2.py",
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
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_2's system map "
                                 "numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_system_map_v0_2.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
