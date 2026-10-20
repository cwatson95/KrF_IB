"""Numbers for v0_1's system map (system_map_v0_1.tex).

v0_1's basic network is v0's with the charge chemistry and fluorine
recombination in it (v0_1/README.md). The map draws it at the reference
point with the same rules as the basic map (make_system_map_basic.py): every
energy flow as a percent of the deposited energy, each node checked to
conserve them, every row of the network drawn. It writes
system_map_v0_1_numbers.tex as one \\vmdef{key}{value} per number, a
macro of its own so that one document can hold v0's and v0_1's numbers.
test_docs.py (v0_1) fails when the file is stale.

    python3 LIGHT/docs/make_system_map_v0_1.py        # rewrite the numbers
    cd LIGHT/docs && latexmk -pdf system_map_v0_1.tex # rebuild the PDF

What v0_1 adds to the map. Ionization no longer ends in a store. Kr+ turns
into Kr2+ (K15), which recombines with an electron into Kr* (R11p), or
with F- into KrF(B) (K26, K108). F- comes from F2 attachment (ATT), and
F2 also loses electrons' energy to dissociation and vibration (DIS, VIB).
F + F + Kr gives F atoms back to F2 (FF). KrF(B) radiates at 248 nm and,
3.5 % as often, at 272 nm (G4A). F- lies 2.6 eV below Kr + F2 on the
ledger's scale, so its flows are negative. The map prints them and draws
no band for them.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_1.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_1.py").is_file())
for _p in (str(LIGHT), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from make_system_map import _pow10, _sci                     # noqa: E402
from make_system_map_basic import _f1, _k, _net_f, _sig2      # noqa: E402
from v0_1 import __version__                                  # noqa: E402
from v0_1.model import Model, Params, run                     # noqa: E402
from v0_1.network import BASIC_GROUPS, KRF_BA_OVER_BX         # noqa: E402
from v0_1.report import (accounting, f2_burnout_time,         # noqa: E402
                         reservoir, row_dict)
from v0_1.species import e_form                               # noqa: E402

OUT = HERE / "system_map_v0_1_numbers.tex"
C1_ROWS = ("C1a", "C1b", "C1c")
HARPOON = ("G3", "G3C", "G3X")
QUENCH_B = ("QF2", "QKR", "QE")
QUENCH_C = ("QF2C", "QKRC", "QEC")
F2E = ("DIS_a3", "DIS_A1", "VIB1", "VIB2", "VIB3", "VIB4")
ION_ION = ("K26", "K108")
CATIONS = ("Kr+", "Kr2+")
BANDS = {"248": "248 nm (KrF B-X)", "272": "272 nm (KrF B-A)",
         "275": "275 nm (KrF C-A)", "420": "420 nm (Kr2F*)"}


def flows(params: Params = None) -> tuple:
    """(flows, context): every energy flow on the v0_1 map as percent of
    E_dep, unrounded, plus the run it came from."""
    p = params or Params()
    if p.network != "basic" or p.groups:
        raise ValueError("the v0_1 map draws the basic network with no "
                         "optional group")
    res = run(p)
    if res.model.ne_mode != "free":
        raise ValueError("the v0_1 map draws free n_e (ne_mode 'auto')")
    a = accounting(res)
    e_dep = a["e_dep"]
    ev = a["events"]
    rows = {c.row.rid: c.row for c in res.model.compiled}
    route = {c.row.rid: c.route for c in res.model.compiled}
    E = {s: e_form(s) for s in ("Kr*", "Kr+", "Kr2+", "F-", "KrF(B)",
                                "KrF(C)", "Kr2F*", "F")}
    krs = reservoir(res, "Kr*")
    k2f = reservoir(res, "Kr2F*")
    used = set()

    def pct(x):
        return 100.0 * x / e_dep

    def ev_(*rids):
        used.update(rids)
        return sum(ev[r] for r in rids)

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

    mix = ev_("CB") - ev_("BC")          # net C -> B events
    f = {
        # electron gas: in = power + Penning return, out = the e-impact rows,
        # elastic loss, the captured electrons' energy and the store
        "e.P": pct(e_dep),
        "e.pen": pct(ev_("G2") * route["G2"]["to_e"]),
        "e.G1": pct(ev_("G1") * route["G1"]["from_e"]),
        "e.C1": pct(sum(ev_(r) * route[r]["from_e"] for r in C1_ROWS)),
        "e.C2": pct(ev_("C2") * route["C2"]["from_e"]),
        "e.F2": pct(sum(ev_(r) * route[r]["from_e"] for r in F2E)),
        "e.G6": pct(a["heat_elastic"]),
        "e.dU": pct(a["stores"]["electron gas"]),
        "e.cap": pct(float(res.ledger["heat_capture"][-1])),
        # C1: the cascade delivers Kr*(5s) and heats the gas with the excess
        "C1.kr": pct(ev_(*C1_ROWS) * E["Kr*"]),
        "C1.heat": heat(*C1_ROWS),
        # F2 under electron impact: dissociation, vibration and attachment
        "F2e.F": f_atoms(*F2E, "ATT"),
        "F2e.heat": heat(*F2E, "ATT"),
        "F2e.anion": pct(ev_("ATT") * E["F-"]),
        # the Kr* reservoir, fed by excitation and by recombination
        "krs.in": pct(krs["formed"] * E["Kr*"]),
        "krs.harp": pct(ev_(*HARPOON) * E["Kr*"]),
        "krs.G2": pct(2 * ev_("G2") * E["Kr*"]),
        "krs.left": left("Kr*"),
        "brs.harp": 100.0 * sum(krs["losses"][r] for r in HARPOON)
        / krs["formed"],
        "brs.G2": 100.0 * krs["losses"]["G2"] / krs["formed"],
        # Kr* + F2: KrF(B), KrF(C), or no KrF at all
        "harp.B": pct(ev_("G3") * E["KrF(B)"]),
        "harp.C": pct(ev_("G3C") * E["KrF(C)"]),
        "harp.F": f_atoms(*HARPOON),
        "harp.heat": heat(*HARPOON),
        # Penning and direct ionization into the positive ions
        "G2.ion": pct(ev_("G2") * E["Kr+"]),
        "C2.ion": pct(ev_("C2") * E["Kr+"]),
        # the positive ions: Kr+ -> Kr2+ (K15), Kr2+ + e -> Kr* (R11p),
        # ion-ion recombination with F- into KrF(B) (K26, K108)
        "K15.heat": heat("K15"),
        "R11p.in": pct(ev_("R11p") * ion_in("R11p")),
        "R11p.kr": pct(ev["R11p"] * E["Kr*"]),
        "R11p.heat": heat("R11p"),
        "ii.in": pct(sum(ev_(r) * ion_in(r) for r in ION_ION)),
        "ii.anion": pct(ev_(*ION_ION) * E["F-"]),
        "ii.B": pct(ev_(*ION_ION) * E["KrF(B)"]),
        "ii.heat": heat(*ION_ION),
        "ions.left": left(*CATIONS),
        "anion.left": left("F-"),
        # B-C mixing, net C -> B: out of C at C's energy, into B at B's
        "mix.net": pct(mix * E["KrF(B)"]),
        "mix.netC": pct(mix * E["KrF(C)"]),
        "mix.heat": heat("BC", "CB"),
        "mix.BC": pct(ev["BC"] * E["KrF(B)"]),
        "mix.CB": pct(ev["CB"] * E["KrF(C)"]),
        # the KrF(B) store
        "B.G4": pct(ev_("G4") * E["KrF(B)"]),
        "B.G4A": pct(ev_("G4A") * E["KrF(B)"]),
        "B.G5": pct(ev_("G5") * E["KrF(B)"]),
        "B.Q": pct(ev_(*QUENCH_B) * E["KrF(B)"]),
        "B.QF2": pct(ev["QF2"] * E["KrF(B)"]),
        "B.QKR": pct(ev["QKR"] * E["KrF(B)"]),
        "B.QE": pct(ev["QE"] * E["KrF(B)"]),
        "B.left": left("KrF(B)"),
        # the KrF(C) store
        "C.G4C": pct(ev_("G4C") * E["KrF(C)"]),
        "C.G5C": pct(ev_("G5C") * E["KrF(C)"]),
        "C.Q": pct(ev_(*QUENCH_C) * E["KrF(C)"]),
        "C.QF2": pct(ev["QF2C"] * E["KrF(C)"]),
        "C.QKR": pct(ev["QKRC"] * E["KrF(C)"]),
        "C.QE": pct(ev["QEC"] * E["KrF(C)"]),
        "C.left": left("KrF(C)"),
        # emission from B (248 nm and 272 nm) and C (the C-A band)
        "G4.hv": rad("G4"),
        "G4.F": f_atoms("G4"),
        "G4.heat": heat("G4"),
        "G4A.hv": rad("G4A"),
        "G4A.F": f_atoms("G4A"),
        "G4A.heat": heat("G4A"),
        "G4C.hv": rad("G4C"),
        "G4C.F": f_atoms("G4C"),
        "G4C.heat": heat("G4C"),
        # quenching of B and of C by F2, Kr and electrons
        "QB.F": f_atoms(*QUENCH_B),
        "QB.heat": heat(*QUENCH_B),
        "QC.F": f_atoms(*QUENCH_C),
        "QC.heat": heat(*QUENCH_C),
        # the trimer, from B and from C
        "tri.k2f": pct(ev_("G5", "G5C") * E["Kr2F*"]),
        "tri.heat": heat("G5", "G5C"),
        # the Kr2F* store, its light and its F2 quench
        "k2f.C3": pct(ev_("C3") * E["Kr2F*"]),
        "k2f.Q2F": pct(ev_("Q2F") * E["Kr2F*"]),
        "k2f.left": left("Kr2F*"),
        "brk.C3": 100.0 * k2f["losses"]["C3"] / k2f["formed"],
        "brk.Q2F": 100.0 * k2f["losses"]["Q2F"] / k2f["formed"],
        "C3.hv": rad("C3"),
        "C3.F": f_atoms("C3"),
        "C3.heat": heat("C3"),
        "Q2F.F": f_atoms("Q2F"),
        "Q2F.heat": heat("Q2F"),
        # F + F + Kr -> F2 + Kr: F atoms back to F2, their energy to heat
        "FF.in": -f_atoms("FF"),
        "FF.heat": heat("FF"),
        # the fate ledger at t_end
        "out.248": pct(a["bands"][BANDS["248"]]),
        "out.272": pct(a["bands"][BANDS["272"]]),
        "out.275": pct(a["bands"][BANDS["275"]]),
        "out.420": pct(a["bands"][BANDS["420"]]),
        "out.heat": pct(a["heat_reactions"]),
        "out.elastic": pct(a["heat_elastic"]),
        "out.ions": pct(a["stores"]["ions (Kr+, Kr2+, F-)"]),
        "out.F": pct(a["stores"]["F atoms (F2 broken)"]),
        "out.other": pct(a["stores"]["Kr*"] + a["stores"]["KrF(B)+KrF(C)"]
                         + a["stores"]["Kr2F*"]
                         + a["stores"]["electron gas"]),
    }
    used.add("G6")                        # G6 is the elastic loss, e.G6
    missing = set(rows) - used
    if missing:
        raise AssertionError(f"rows not drawn on the map: {sorted(missing)}")
    if mix <= 0.0:
        raise AssertionError("net B-C mixing now runs B -> C: redraw the "
                             "mixing band in system_map_v0_1.tex")
    if set(a["bands"]) != set(BANDS.values()):
        raise AssertionError(f"bands changed: {sorted(a['bands'])}")
    _check_conservation(f, a["closure"])
    krf = reservoir(res, "KrF*")
    extra = {"res": res, "params": p, "krs": krs, "krf": krf, "k2f": k2f,
             "rows": rows, "events": ev}
    return f, extra


def _check_conservation(f: dict, closure: float, tol: float = 1e-8) -> None:
    """Each node of the map must balance before a number is printed."""
    nodes = {
        "electron gas": (f["e.P"] + f["e.pen"],
                         f["e.G1"] + f["e.C1"] + f["e.C2"] + f["e.F2"]
                         + f["e.G6"] + f["e.dU"] + f["e.cap"]),
        "C1": (f["e.C1"], f["C1.kr"] + f["C1.heat"]),
        "F2 electron rows": (f["e.F2"], f["F2e.F"] + f["F2e.heat"]
                             + f["F2e.anion"]),
        "Kr*": (f["e.G1"] + f["C1.kr"] + f["R11p.kr"],
                f["krs.harp"] + f["krs.G2"] + f["krs.left"]),
        "Kr* in": (f["krs.in"], f["e.G1"] + f["C1.kr"] + f["R11p.kr"]),
        "harpoon": (f["krs.harp"], f["harp.B"] + f["harp.C"] + f["harp.F"]
                    + f["harp.heat"]),
        "G2": (f["krs.G2"], f["G2.ion"] + f["e.pen"]),
        "positive ions": (f["C2.ion"] + f["G2.ion"],
                          f["K15.heat"] + f["R11p.in"] + f["ii.in"]
                          + f["ions.left"]),
        "R11p": (f["R11p.in"], f["R11p.kr"] + f["R11p.heat"]),
        "F-": (f["F2e.anion"], f["ii.anion"] + f["anion.left"]),
        "ion-ion": (f["ii.in"] + f["ii.anion"], f["ii.B"] + f["ii.heat"]),
        "mixing": (f["mix.netC"], f["mix.net"] + f["mix.heat"]),
        "KrF(B)": (f["harp.B"] + f["mix.net"] + f["ii.B"],
                   f["B.G4"] + f["B.G4A"] + f["B.G5"] + f["B.Q"]
                   + f["B.left"]),
        "KrF(C)": (f["harp.C"], f["mix.netC"] + f["C.G4C"] + f["C.G5C"]
                   + f["C.Q"] + f["C.left"]),
        "B quench split": (f["B.Q"], f["B.QF2"] + f["B.QKR"] + f["B.QE"]),
        "C quench split": (f["C.Q"], f["C.QF2"] + f["C.QKR"] + f["C.QE"]),
        "G4": (f["B.G4"], f["G4.hv"] + f["G4.F"] + f["G4.heat"]),
        "G4A": (f["B.G4A"], f["G4A.hv"] + f["G4A.F"] + f["G4A.heat"]),
        "G4C": (f["C.G4C"], f["G4C.hv"] + f["G4C.F"] + f["G4C.heat"]),
        "B quench": (f["B.Q"], f["QB.F"] + f["QB.heat"]),
        "C quench": (f["C.Q"], f["QC.F"] + f["QC.heat"]),
        "trimer": (f["B.G5"] + f["C.G5C"], f["tri.k2f"] + f["tri.heat"]),
        "Kr2F*": (f["tri.k2f"], f["k2f.C3"] + f["k2f.Q2F"] + f["k2f.left"]),
        "C3": (f["k2f.C3"], f["C3.hv"] + f["C3.F"] + f["C3.heat"]),
        "Q2F": (f["k2f.Q2F"], f["Q2F.F"] + f["Q2F.heat"]),
        "FF": (f["FF.in"], f["FF.heat"]),
        "ions": (f["ions.left"] + f["anion.left"], f["out.ions"]),
        "F atoms": (f["harp.F"] + f["G4.F"] + f["G4A.F"] + f["G4C.F"]
                    + f["QB.F"] + f["QC.F"] + f["C3.F"] + f["Q2F.F"]
                    + f["F2e.F"], f["FF.in"] + f["out.F"]),
        "heat": (f["C1.heat"] + f["harp.heat"] + f["mix.heat"] + f["G4.heat"]
                 + f["G4A.heat"] + f["G4C.heat"] + f["QB.heat"]
                 + f["QC.heat"] + f["tri.heat"] + f["C3.heat"]
                 + f["Q2F.heat"] + f["e.cap"] + f["F2e.heat"]
                 + f["K15.heat"] + f["R11p.heat"] + f["ii.heat"]
                 + f["FF.heat"], f["out.heat"]),
        "248 nm": (f["G4.hv"], f["out.248"]),
        "272 nm": (f["G4A.hv"], f["out.272"]),
        "C-A band": (f["G4C.hv"], f["out.275"]),
        "420 nm": (f["C3.hv"], f["out.420"]),
        "other stores": (f["krs.left"] + f["B.left"] + f["C.left"]
                         + f["k2f.left"] + f["e.dU"], f["out.other"]),
        "ledger": (100.0, f["out.248"] + f["out.272"] + f["out.275"]
                   + f["out.420"] + f["out.heat"] + f["out.elastic"]
                   + f["out.ions"] + f["out.F"] + f["out.other"]),
    }
    for name, (into, out) in nodes.items():
        if abs(into - out) > tol:
            raise AssertionError(f"map node {name!r} does not balance: "
                                 f"{into:.12f} in, {out:.12f} out")
    if abs(closure) > 1e-12:
        raise AssertionError(f"ledger closure {closure:.2e} above 1e-12")


def _rate(rows, rid) -> float:
    kind, *args = rows[rid].rate
    if kind != "const":
        raise ValueError(f"{rid}: rate kind {kind!r} has no single value")
    return float(args[0])


def numbers(params: Params = None) -> dict:
    """The formatted strings the map prints, keyed as in the .tex."""
    f, x = flows(params)
    p, res, rows, ev = x["params"], x["res"], x["rows"], x["events"]
    n = {k: _f1(v) for k, v in f.items()}
    # the G4 box draws B -> X and B -> A together: its tags add the two
    n["G4x.F"] = _f1(f["G4.F"] + f["G4A.F"])
    n["G4x.heat"] = _f1(f["G4.heat"] + f["G4A.heat"])
    on = res.t <= p.pulse_s
    tb = f2_burnout_time(res)
    m = res.model
    # the electron closure at mid-pulse: the tail temperature T_t = r Te
    k = int(res.t.searchsorted(0.5 * p.pulse_s))
    n_gas = float(sum(res.n[s][k] for s in m.species))
    r_tail = m.tail(res.ne[k] / n_gas, res.n["F2"][k] / n_gas)
    krf = x["krf"]
    lost = krf["losses"]

    def brf(*rids):
        return f"{100.0 * sum(lost.get(r, 0.0) for r in rids) / krf['formed']:.1f}"
    k_h = _rate(rows, "G3") + _rate(rows, "G3C") + _rate(rows, "G3X")
    y_h = (_rate(rows, "G3") + _rate(rows, "G3C")) / k_h
    f_c = _rate(rows, "G3C") / (_rate(rows, "G3") + _rate(rows, "G3C"))
    k_eq = rows["CB"].rate[2]
    formed_ii = sum(ev[r] for r in ("K26", "K108"))
    formed_all = ev["G3"] + ev["G3C"] + formed_ii
    kind, k0, te0, n_pow = rows["R11p"].rate
    assert kind == "te_pow", rows["R11p"].rate
    burn = tb if tb is not None else p.pulse_s
    before = res.t <= burn
    n.update({
        "version": __version__,
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
        "F2left": f"{100 * res.n['F2'][-1] / res.n['F2'][0]:.2f}",
        "krsFormed": _sci(x["krs"]["formed"]),
        "krfFormed": _sci(krf["formed"]),
        "k2fFormed": _sci(x["k2f"]["formed"]),
        "krf.eff": _f1(100.0 * krf["formed"] * rows["G4"].hv_eV
                       / float(res.ledger["E_dep"][-1])),
        # where KrF(B,C) comes from, percent of the KrF formed
        "frm.ii": _f1(100.0 * formed_ii / formed_all),
        # where KrF(B,C) goes, percent of the KrF formed (B and C together)
        "brf.G4": brf("G4"), "brf.G4A": brf("G4A"), "brf.G4C": brf("G4C"),
        "brf.tri": brf("G5", "G5C"), "brf.e": brf("QE", "QEC"),
        "brf.F2": brf("QF2", "QF2C"), "brf.Kr": brf("QKR", "QKRC"),
        "brf.QB": brf(*QUENCH_B), "brf.QC": brf(*QUENCH_C),
        # how the ions recombine: with an electron (R11p) or with F-
        "ion.R11p": _f1(100.0 * ev["R11p"] / (ev["R11p"] + formed_ii)),
        "ion.ii": _f1(100.0 * formed_ii / (ev["R11p"] + formed_ii)),
        # the constants the boxes print, read from the rows the model ran
        "k.harp": _k(k_h),
        "pc.harpY": f"{100 * y_h:.0f}",
        "pc.harpB": f"{100 * y_h * (1 - f_c):.0f}",
        "pc.harpC": f"{100 * y_h * f_c:.0f}",
        "pc.harpX": f"{100 * (1 - y_h):.0f}",
        "pc.fC": f"{100 * f_c:.0f}",
        "k.G2": _k(_rate(rows, "G2")),
        "tauB": f"{p.tau_sp_s * 1e9:.2f}",
        "pc.BA": f"{100 * KRF_BA_OVER_BX:.1f}",
        "tauC": _sig2(1e9 / _rate(rows, "G4C")),
        "tauC3": _sig2(1e9 / _rate(rows, "C3")),
        "gapBC": f"{(e_form('KrF(B)') - e_form('KrF(C)')) / 1.23984e-4:.0f}",
        "k.BC": _k(_rate(rows, "BC")),
        "Keq": f"{k_eq:.2f}",
        "k.G5": _k(_rate(rows, "G5")),
        "k.G5C": _k(_rate(rows, "G5C")),
        "tpow.G5": f"{rows['G5'].tg_pow:.2f}",
        "k.Q2F": _k(_rate(rows, "Q2F")),
        "k.K15": _k(_rate(rows, "K15")),
        "tpow.K15": f"{rows['K15'].tg_pow:.2f}",
        "k.R11p": _k(k0),
        "pow.R11p": f"{-n_pow:.2f}",
        "cap.R11p": f"{1.5 + n_pow:.2f}",
        "k.FF": _k(_rate(rows, "FF")),
        "n3.ii": f"{rows['K26'].rate[3]:g}",
        # how far F- lies below Kr + F2 on the ledger's scale (a magnitude)
        "E.anion": f"{-e_form('F-'):.1f}",
        "hv.G4C": f"{rows['G4C'].hv_eV:.2f}",
        "hv.G4": f"{rows['G4'].hv_eV:.2f}",
        "hv.G4A": f"{rows['G4A'].hv_eV:.2f}",
        "hv.C3": f"{rows['C3'].hv_eV:.2f}",
        "E.B": f"{e_form('KrF(B)'):.3f}",
        # the notes' output line (flows, they overlap)
        "notes.krs": _f1(f["e.G1"] + f["e.C1"]),
        "notes.krp": _f1(f["C2.ion"] + f["G2.ion"]),
        "notes.emit": _f1(f["out.248"] + f["out.272"] + f["out.275"]
                          + f["out.420"]),
    })
    # what the basic network leaves out: the 248-nm share with each
    # optional group on, and the notes' list as written (v0 network)
    for g in BASIC_GROUPS + ("all",):
        d = row_dict(run(p.replace(groups=(g,))))
        n[f"om.{g}"] = f"{d['emission_248']:.1f}"
    n["om.basic"] = n["out.248"]
    # the Kr2* group without Kr2* + F -> KrF(B) + Kr, Kannari's estimate
    m = Model(p.replace(groups=("kr2",)))
    m.compiled = [c for c in m.compiled if c.row.rid != "K2F_F"]
    n["om.kr2noF"] = f"{row_dict(m.run())['emission_248']:.1f}"
    n["om.v0"] = f"{row_dict(run(p.replace(network='v0')))['emission_248']:.1f}"
    # the trimer band: the same point with the pure-Kr T^-1/2 law
    d = row_dict(run(p.replace(trimer_tg_pow=-0.5)))
    n["om.trimer"] = f"{d['emission_248']:.1f}"
    return n


def render(n: dict) -> str:
    lines = [
        "% LIGHT v0_1 map numbers. GENERATED by make_system_map_v0_1.py",
        f"% from the model ({n['version']}). Do not edit. Rerun the script "
        "instead.",
        "\\providecommand\\vmdef[2]{\\expandafter\\def"
        "\\csname vm@#1\\endcsname{#2}}",
    ]
    lines += [f"\\vmdef{{{k}}}{{{v}}}" for k, v in n.items()]
    return "\n".join(lines) + "\n"


def main() -> int:
    OUT.write_text(render(numbers()))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
