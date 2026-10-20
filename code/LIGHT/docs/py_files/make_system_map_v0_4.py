"""Numbers for v0_4's system map (system_map_v0_4.tex).

v0_4's basic network is v0_3's with the 5p rung of krypton (Kr**, fed by
the higher excitations C1a-c and by Kr2+ + e, pumped from Kr* by the
electrons and returned, relaxed by Kr, radiating to Kr*, harpooning F2),
KrF(B,C) + F as a band, ion-ion recombination with its yield and B/C
bands, and the 248-nm optics (v0_4/README.md). The map draws the network
at the reference point with every band at its default, under v0_2's
electron closure (the kinetic closure moves the reference point's bands
by about 2 % of themselves, v0_4/README.md), with the rules of the basic
map (make_system_map_basic.py): every energy flow as a percent of the
deposited energy, each node checked to conserve them, every row of the
network drawn. It writes system_map_v0_4_numbers.tex as one
\\wmdef{key}{value} per number. test_docs.py (v0_4) fails when the file is
stale.

    python3 LIGHT/docs/make_system_map_v0_4.py        # about five minutes
    cd LIGHT/docs && latexmk -pdf system_map_v0_4.tex # rebuild the PDF

What v0_4 adds to the map. A Kr** node between C1 and Kr*. The 5s <-> 5p
electron cycle (K76 up, K76r back, each event moving 9.915 eV of Kr*
energy and 1.385 eV of electron energy) circulates more than twice the
input, so the map draws it as its net: the electrons' net contribution,
(K76 - K76r) x 1.385 eV, enters the Kr** node from the electron gas, and
the Kr* energy the rung hands back, (KPS + PRAD + K76r - K76) x 9.915 eV,
runs from Kr** to Kr* as one band. KPS leaves the gap as heat and PRAD as
a near-IR photon (a new ledger line). K27, K27C and K27X, the harpoon
from 5p, join the harpoon box. Rows with no measurement are 0 by default
and are printed, not drawn: KrF(B,C) + F (QF, QFC), the stepwise
ionizations (KSTEP_S, KSTEP_P), the C and dissociating branches of
ion-ion recombination (K26C, K26X, K108C, K108X) and v0_2's three. The
bands box gives the 248-nm share at the two ends of every band
(run_v0_4.py --bands), and the footer the 248-nm gain and absorption.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_4.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_4.py").is_file())
for _p in (str(LIGHT), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from make_system_map import _pow10, _sci                     # noqa: E402
from make_system_map_basic import _f1, _k, _net_f, _sig2      # noqa: E402
from v0_4 import __version__, optics                          # noqa: E402
from v0_4.model import Params, run                            # noqa: E402
from v0_4.network import (BANDS, BASIC_GROUPS, KRF_BA_OVER_BX,  # noqa: E402
                          KR_5P_TAU_S)
from v0_4.report import (accounting, f2_burnout_time,         # noqa: E402
                         reservoir, row_dict)
from v0_4.run import do_bands                                 # noqa: E402
from v0_4.species import e_form                               # noqa: E402

OUT = HERE / "system_map_v0_4_numbers.tex"
C1_ROWS = ("C1a", "C1b", "C1c")
HARPOON = ("G3", "G3C", "G3X")
HARPOON_P = ("K27", "K27C", "K27X")
QUENCH_B = ("QF2", "QKR", "QE", "QF")
QUENCH_C = ("QF2C", "QKRC", "QEC", "QFC")
F2E = ("DIS_a3", "DIS_A1", "VIB1", "VIB2", "VIB3", "VIB4")
ION_ION = ("K26", "K108", "K26C", "K108C", "K26X", "K108X")
ION_ION_B = ("K26", "K108")
ION_ION_C = ("K26C", "K108C")
ION_ION_X = ("K26X", "K108X")
KR2_F2 = ("NK7", "K2F_BR")
CATIONS = ("Kr+", "Kr2+")
# the rows with no measurement, 0 unless a band run sets them: the map
# prints them and draws no band
ZERO_ROWS = ("K2F_F", "PEN_KR2", "Q2FF", "QF", "QFC", "KSTEP_S", "KSTEP_P",
             "K26C", "K108C", "K26X", "K108X")
BANDS_OUT = {"248": "248 nm (KrF B-X)", "272": "272 nm (KrF B-A)",
             "275": "275 nm (KrF C-A)", "420": "420 nm (Kr2F*)",
             "147": "147 nm (Kr2* VUV)", "IR": "5p-5s (near IR)"}
# the order the map's box and the overview's table list the bands in: the
# wide ones first, v0_2's eight then v0_4's nine
BAND_ORDER = ("trimer_tg_pow", "krf_f", "kr2_f", "k_5s5p", "kstep_s",
              "kr2_form", "kstep_p", "kr_pen", "kr2f_f2", "kr2f_f",
              "ii_yield", "ii_f_c", "kps_kr", "k5p_f2", "kr2_f2", "kr2_pen",
              "ib_scale")


def _check_defaults(p: Params) -> None:
    if p.network != "basic" or p.groups:
        raise ValueError("the v0_4 map draws the basic network with no "
                         "optional group")
    if p.rates or p.trimer_tg_pow is not None or p.ladder != "off":
        raise ValueError("the v0_4 map draws every band at its default; a "
                         "band run needs its own drawing")
    if p.gap_cm is not None or p.em_E_V_cm > 0.0:
        raise ValueError("the v0_4 map draws the box of the notes: no walls, "
                         "no wave")


def flows(params: Params = None) -> tuple:
    """(flows, context): every energy flow on the v0_4 map as percent of
    E_dep, unrounded, plus the run it came from."""
    p = params or Params(eedf="two_t")
    _check_defaults(p)
    res = run(p)
    if res.model.ne_mode != "free":
        raise ValueError("the v0_4 map draws free n_e (ne_mode 'auto')")
    a = accounting(res)
    e_dep = a["e_dep"]
    ev = a["events"]
    rows = {c.row.rid: c.row for c in res.model.compiled}
    route = {c.row.rid: c.route for c in res.model.compiled}
    E = {s: e_form(s) for s in ("Kr*", "Kr**", "Kr+", "Kr2+", "F-", "KrF(B)",
                                "KrF(C)", "Kr2F*", "F", "Kr2*")}
    krs = reservoir(res, "Kr*")
    kpp = reservoir(res, "Kr**")
    k2s = reservoir(res, "Kr2*")
    k2f = reservoir(res, "Kr2F*")
    used = set()
    nonzero = [r for r in ZERO_ROWS if ev[r] != 0.0]
    if nonzero:
        raise AssertionError(f"rows the map prints as 0 carry energy: "
                             f"{nonzero}. Redraw them in system_map_v0_4.tex")

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
    gap = E["Kr**"] - E["Kr*"]
    net76 = ev_("K76") - ev_("K76r")     # net 5s -> 5p events
    # the Kr* energy the rung hands back to Kr*, net: KPS and PRAD return
    # it, K76 takes it up and K76r brings it down
    back = (ev_("KPS") + ev_("PRAD") - net76) * E["Kr*"]
    f = {
        # electron gas: in = power + the two Penning returns, out = the
        # e-impact rows, the net 5s -> 5p pumping, the stepwise rows,
        # elastic loss, the captured electrons' energy and the store
        "e.P": pct(e_dep),
        "e.pen": pct(ev_("G2") * route["G2"]["to_e"]),
        "e.pen2": pct(ev_("PEN_KR2") * route["PEN_KR2"]["to_e"]),
        "e.G1": pct(ev_("G1") * route["G1"]["from_e"]),
        "e.C1": pct(sum(ev_(r) * route[r]["from_e"] for r in C1_ROWS)),
        "e.C2": pct(ev_("C2") * route["C2"]["from_e"]),
        "e.F2": pct(sum(ev_(r) * route[r]["from_e"] for r in F2E)),
        "e.76": pct(net76 * gap),
        "e.stepS": pct(ev_("KSTEP_S") * route["KSTEP_S"]["from_e"]),
        "e.stepP": pct(ev_("KSTEP_P") * route["KSTEP_P"]["from_e"]),
        "e.G6": pct(a["heat_elastic"]),
        "e.dU": pct(a["stores"]["electron gas"]),
        "e.cap": pct(float(res.ledger["heat_capture"][-1])),
        # C1: the higher excitations deliver Kr** (5p) and heat the gas
        # with the excess over 11.30 eV
        "C1.kr": pct(ev_(*C1_ROWS) * E["Kr**"]),
        "C1.heat": heat(*C1_ROWS),
        # F2 under electron impact: dissociation, vibration and attachment
        "F2e.F": f_atoms(*F2E, "ATT"),
        "F2e.heat": heat(*F2E, "ATT"),
        "F2e.anion": pct(ev_("ATT") * E["F-"]),
        # the Kr** reservoir: in from C1, R11p and the net pumping (Kr*
        # energy and electron energy); out as the net return to Kr*, the
        # relaxation's heat, the near-IR photon, the 5p harpoon and the
        # stepwise ionization
        "kpp.in": pct(kpp["formed"] * E["Kr**"]),
        "kpp.back": pct(back),
        "KPS.heat": heat("KPS"),
        "PRAD.hv": rad("PRAD"),
        "PRAD.heat": heat("PRAD"),
        "kpp.harp": pct(ev_(*HARPOON_P) * E["Kr**"]),
        "kpp.stepP": pct(ev_("KSTEP_P") * E["Kr**"]),
        "kpp.left": left("Kr**"),
        "brp.KPS": 100.0 * kpp["losses"]["KPS"] / kpp["formed"],
        "brp.PRAD": 100.0 * kpp["losses"]["PRAD"] / kpp["formed"],
        "brp.76r": 100.0 * kpp["losses"]["K76r"] / kpp["formed"],
        "brp.harp": 100.0 * sum(kpp["losses"][r] for r in HARPOON_P)
        / kpp["formed"],
        # Kr(5p) + F2: KrF(B), KrF(C), or no KrF at all
        "harpP.B": pct(ev_("K27") * E["KrF(B)"]),
        "harpP.C": pct(ev_("K27C") * E["KrF(C)"]),
        "harpP.F": f_atoms(*HARPOON_P),
        "harpP.heat": heat(*HARPOON_P),
        # the Kr* reservoir, fed by excitation and by the rung's return
        "krs.in": pct(krs["formed"] * E["Kr*"]),
        "krs.harp": pct(ev_(*HARPOON) * E["Kr*"]),
        "krs.G2": pct(2 * ev_("G2") * E["Kr*"]),
        "krs.AK8": pct(ev_("AK8_Kr") * E["Kr*"]),
        "krs.stepS": pct(ev_("KSTEP_S") * E["Kr*"]),
        "krs.left": left("Kr*"),
        "brs.harp": 100.0 * sum(krs["losses"][r] for r in HARPOON)
        / krs["formed"],
        "brs.G2": 100.0 * krs["losses"]["G2"] / krs["formed"],
        "brs.AK8": 100.0 * krs["losses"]["AK8_Kr"] / krs["formed"],
        "brs.76": 100.0 * krs["losses"]["K76"] / krs["formed"],
        # Kr* + F2: KrF(B), KrF(C), or no KrF at all
        "harp.B": pct(ev_("G3") * E["KrF(B)"]),
        "harp.C": pct(ev_("G3C") * E["KrF(C)"]),
        "harp.F": f_atoms(*HARPOON),
        "harp.heat": heat(*HARPOON),
        # Penning, direct and stepwise ionization into the positive ions
        "G2.ion": pct(ev_("G2") * E["Kr+"]),
        "C2.ion": pct(ev_("C2") * E["Kr+"]),
        "stepS.ion": pct(ev["KSTEP_S"] * E["Kr+"]),
        "stepP.ion": pct(ev["KSTEP_P"] * E["Kr+"]),
        # the Kr2* node: Kr* + 2Kr (AK8_Kr) inside it, then 147-nm light,
        # Kr2* + F2, and the two rows at 0
        "AK8.heat": heat("AK8_Kr"),
        "k2s.in": pct(k2s["formed"] * E["Kr2*"]),
        "k2s.RAD": pct(ev_("KR2_RAD") * E["Kr2*"]),
        "k2s.F2": pct(ev_(*KR2_F2) * E["Kr2*"]),
        "k2s.F": pct(ev_("K2F_F") * E["Kr2*"]),
        "k2s.pen": pct(2 * ev_("PEN_KR2") * E["Kr2*"]),
        "k2s.left": left("Kr2*"),
        "brk2.RAD": 100.0 * k2s["losses"]["KR2_RAD"] / k2s["formed"],
        "brk2.F2": 100.0 * sum(k2s["losses"][r] for r in KR2_F2)
        / k2s["formed"],
        "RAD.hv": rad("KR2_RAD"),
        "RAD.heat": heat("KR2_RAD"),
        # Kr2* + F2: KrF(B) + Kr + F (NK7) or Kr2F* + F (K2F_BR)
        "NK7.B": pct(ev["NK7"] * E["KrF(B)"]),
        "BR.k2f": pct(ev["K2F_BR"] * E["Kr2F*"]),
        "K2F2.F": f_atoms(*KR2_F2),
        "K2F2.heat": heat(*KR2_F2),
        # Kr2* + F -> KrF(B) + Kr: its F atom comes in from the F store
        "k2F.B": pct(ev["K2F_F"] * E["KrF(B)"]),
        "k2F.Fin": -f_atoms("K2F_F"),
        "k2F.heat": heat("K2F_F"),
        # Kr2* + Kr2* -> Kr2+ + 2Kr + e
        "pen2.ion": pct(ev["PEN_KR2"] * E["Kr2+"]),
        # the positive ions: Kr+ -> Kr2+ (K15), Kr2+ + e -> Kr** (R11p),
        # ion-ion recombination with F- (K26, K108 and their bands)
        "K15.heat": heat("K15"),
        "R11p.in": pct(ev_("R11p") * ion_in("R11p")),
        "R11p.kr": pct(ev["R11p"] * E["Kr**"]),
        "R11p.heat": heat("R11p"),
        "ii.in": pct(sum(ev_(r) * ion_in(r) for r in ION_ION)),
        "ii.anion": pct(ev_(*ION_ION) * E["F-"]),
        "ii.B": pct(ev_(*ION_ION_B) * E["KrF(B)"]),
        "ii.C": pct(ev_(*ION_ION_C) * E["KrF(C)"]),
        "ii.X": f_atoms(*ION_ION_X),
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
        "B.QF": pct(ev["QF"] * E["KrF(B)"]),
        "B.left": left("KrF(B)"),
        # the KrF(C) store
        "C.G4C": pct(ev_("G4C") * E["KrF(C)"]),
        "C.G5C": pct(ev_("G5C") * E["KrF(C)"]),
        "C.Q": pct(ev_(*QUENCH_C) * E["KrF(C)"]),
        "C.QF2": pct(ev["QF2C"] * E["KrF(C)"]),
        "C.QKR": pct(ev["QKRC"] * E["KrF(C)"]),
        "C.QE": pct(ev["QEC"] * E["KrF(C)"]),
        "C.QF": pct(ev["QFC"] * E["KrF(C)"]),
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
        # quenching of B and of C by F2, Kr, electrons and F
        "QB.F": f_atoms(*QUENCH_B),
        "QB.heat": heat(*QUENCH_B),
        "QC.F": f_atoms(*QUENCH_C),
        "QC.heat": heat(*QUENCH_C),
        # the trimer, from B and from C
        "tri.k2f": pct(ev_("G5", "G5C") * E["Kr2F*"]),
        "tri.heat": heat("G5", "G5C"),
        # the Kr2F* store, its light, its F2 quench and its F quench (0)
        "k2f.C3": pct(ev_("C3") * E["Kr2F*"]),
        "k2f.Q2F": pct(ev_("Q2F") * E["Kr2F*"]),
        "k2f.Q2FF": pct(ev_("Q2FF") * E["Kr2F*"]),
        "k2f.left": left("Kr2F*"),
        "brk.C3": 100.0 * k2f["losses"]["C3"] / k2f["formed"],
        "brk.Q2F": 100.0 * k2f["losses"]["Q2F"] / k2f["formed"],
        "C3.hv": rad("C3"),
        "C3.F": f_atoms("C3"),
        "C3.heat": heat("C3"),
        "Q2F.F": f_atoms("Q2F"),
        "Q2F.heat": heat("Q2F"),
        "Q2FF.F": f_atoms("Q2FF"),
        "Q2FF.heat": heat("Q2FF"),
        # F + F + Kr -> F2 + Kr: F atoms back to F2, their energy to heat
        "FF.in": -f_atoms("FF"),
        "FF.heat": heat("FF"),
        # the fate ledger at t_end
        "out.248": pct(a["bands"][BANDS_OUT["248"]]),
        "out.272": pct(a["bands"][BANDS_OUT["272"]]),
        "out.275": pct(a["bands"][BANDS_OUT["275"]]),
        "out.420": pct(a["bands"][BANDS_OUT["420"]]),
        "out.147": pct(a["bands"][BANDS_OUT["147"]]),
        "out.IR": pct(a["bands"][BANDS_OUT["IR"]]),
        "out.heat": pct(a["heat_reactions"]),
        "out.elastic": pct(a["heat_elastic"]),
        "out.ions": pct(a["stores"]["ions (Kr+, Kr2+, F-)"]),
        "out.F": pct(a["stores"]["F atoms (F2 broken)"]),
        "out.other": pct(a["stores"]["Kr*"] + a["stores"]["Kr** (5p)"]
                         + a["stores"]["Kr2*"]
                         + a["stores"]["KrF(B)+KrF(C)"]
                         + a["stores"]["Kr2F*"]
                         + a["stores"]["electron gas"]),
    }
    used.add("G6")                        # G6 is the elastic loss, e.G6
    missing = set(rows) - used
    if missing:
        raise AssertionError(f"rows not drawn on the map: {sorted(missing)}")
    if mix <= 0.0:
        raise AssertionError("net B-C mixing now runs B -> C: redraw the "
                             "mixing band in system_map_v0_4.tex")
    if net76 <= 0.0 or back <= 0.0:
        raise AssertionError("the 5s <-> 5p cycle runs the other way: "
                             "redraw its two net bands in system_map_v0_4.tex")
    if set(a["bands"]) != set(BANDS_OUT.values()):
        raise AssertionError(f"bands changed: {sorted(a['bands'])}")
    if a["walls"] != 0.0 or a["e_em"] != 0.0:
        raise AssertionError("the map's run has walls or a wave")
    _check_conservation(f, a["closure"])
    krf = reservoir(res, "KrF*")
    extra = {"res": res, "params": p, "krs": krs, "kpp": kpp, "krf": krf,
             "k2f": k2f, "k2s": k2s, "rows": rows, "events": ev,
             "net76": net76}
    return f, extra


def _check_conservation(f: dict, closure: float, tol: float = 1e-8) -> None:
    """Each node of the map must balance before a number is printed."""
    for name, (into, out) in _nodes(f).items():
        if abs(into - out) > tol:
            raise AssertionError(f"map node {name!r} does not balance: "
                                 f"{into:.12f} in, {out:.12f} out")
    if abs(closure) > 1e-12:
        raise AssertionError(f"ledger closure {closure:.2e} above 1e-12")


def _nodes(f: dict) -> dict:
    """Every node of the map as (energy in, energy out), percent of E_dep."""
    return {
        "electron gas": (f["e.P"] + f["e.pen"] + f["e.pen2"],
                         f["e.G1"] + f["e.C1"] + f["e.C2"] + f["e.F2"]
                         + f["e.76"] + f["e.stepS"] + f["e.stepP"]
                         + f["e.G6"] + f["e.dU"] + f["e.cap"]),
        "C1": (f["e.C1"], f["C1.kr"] + f["C1.heat"]),
        "F2 electron rows": (f["e.F2"], f["F2e.F"] + f["F2e.heat"]
                             + f["F2e.anion"]),
        # the rung: C1 and R11p bring Kr** in, the electrons' net pumping
        # brings the gap's energy; out go the net return to Kr*, the heat
        # and the photon of the gap, the 5p harpoon and the stepwise row
        "Kr**": (f["C1.kr"] + f["R11p.kr"] + f["e.76"],
                 f["kpp.back"] + f["KPS.heat"] + f["PRAD.hv"]
                 + f["PRAD.heat"] + f["kpp.harp"] + f["kpp.stepP"]
                 + f["kpp.left"]),
        "5p harpoon": (f["kpp.harp"], f["harpP.B"] + f["harpP.C"]
                       + f["harpP.F"] + f["harpP.heat"]),
        "KSTEP_P": (f["kpp.stepP"] + f["e.stepP"], f["stepP.ion"]),
        "Kr*": (f["e.G1"] + f["kpp.back"],
                f["krs.harp"] + f["krs.G2"] + f["krs.AK8"] + f["krs.stepS"]
                + f["krs.left"]),
        "KSTEP_S": (f["krs.stepS"] + f["e.stepS"], f["stepS.ion"]),
        "harpoon": (f["krs.harp"], f["harp.B"] + f["harp.C"] + f["harp.F"]
                    + f["harp.heat"]),
        "G2": (f["krs.G2"], f["G2.ion"] + f["e.pen"]),
        "Kr2* in": (f["k2s.in"], f["krs.AK8"] - f["AK8.heat"]),
        "Kr2* (with AK8_Kr)": (f["krs.AK8"],
                               f["AK8.heat"] + f["k2s.RAD"] + f["k2s.F2"]
                               + f["k2s.F"] + f["k2s.pen"] + f["k2s.left"]),
        "KR2_RAD": (f["k2s.RAD"], f["RAD.hv"] + f["RAD.heat"]),
        "Kr2* + F2": (f["k2s.F2"], f["NK7.B"] + f["BR.k2f"] + f["K2F2.F"]
                      + f["K2F2.heat"]),
        "Kr2* + F": (f["k2s.F"] + f["k2F.Fin"], f["k2F.B"] + f["k2F.heat"]),
        "Kr2* + Kr2*": (f["k2s.pen"], f["pen2.ion"] + f["e.pen2"]),
        "positive ions": (f["C2.ion"] + f["G2.ion"] + f["pen2.ion"]
                          + f["stepS.ion"] + f["stepP.ion"],
                          f["K15.heat"] + f["R11p.in"] + f["ii.in"]
                          + f["ions.left"]),
        "R11p": (f["R11p.in"], f["R11p.kr"] + f["R11p.heat"]),
        "F-": (f["F2e.anion"], f["ii.anion"] + f["anion.left"]),
        "ion-ion": (f["ii.in"] + f["ii.anion"], f["ii.B"] + f["ii.C"]
                    + f["ii.X"] + f["ii.heat"]),
        "mixing": (f["mix.netC"], f["mix.net"] + f["mix.heat"]),
        "KrF(B)": (f["harp.B"] + f["harpP.B"] + f["mix.net"] + f["ii.B"]
                   + f["NK7.B"] + f["k2F.B"],
                   f["B.G4"] + f["B.G4A"] + f["B.G5"] + f["B.Q"]
                   + f["B.left"]),
        "KrF(C)": (f["harp.C"] + f["harpP.C"] + f["ii.C"],
                   f["mix.netC"] + f["C.G4C"] + f["C.G5C"] + f["C.Q"]
                   + f["C.left"]),
        "B quench split": (f["B.Q"], f["B.QF2"] + f["B.QKR"] + f["B.QE"]
                           + f["B.QF"]),
        "C quench split": (f["C.Q"], f["C.QF2"] + f["C.QKR"] + f["C.QE"]
                           + f["C.QF"]),
        "G4": (f["B.G4"], f["G4.hv"] + f["G4.F"] + f["G4.heat"]),
        "G4A": (f["B.G4A"], f["G4A.hv"] + f["G4A.F"] + f["G4A.heat"]),
        "G4C": (f["C.G4C"], f["G4C.hv"] + f["G4C.F"] + f["G4C.heat"]),
        "B quench": (f["B.Q"], f["QB.F"] + f["QB.heat"]),
        "C quench": (f["C.Q"], f["QC.F"] + f["QC.heat"]),
        "trimer": (f["B.G5"] + f["C.G5C"], f["tri.k2f"] + f["tri.heat"]),
        "Kr2F*": (f["tri.k2f"] + f["BR.k2f"],
                  f["k2f.C3"] + f["k2f.Q2F"] + f["k2f.Q2FF"]
                  + f["k2f.left"]),
        "C3": (f["k2f.C3"], f["C3.hv"] + f["C3.F"] + f["C3.heat"]),
        "Q2F": (f["k2f.Q2F"], f["Q2F.F"] + f["Q2F.heat"]),
        "Q2FF": (f["k2f.Q2FF"], f["Q2FF.F"] + f["Q2FF.heat"]),
        "FF": (f["FF.in"], f["FF.heat"]),
        "ions": (f["ions.left"] + f["anion.left"], f["out.ions"]),
        "F atoms": (f["harp.F"] + f["harpP.F"] + f["G4.F"] + f["G4A.F"]
                    + f["G4C.F"] + f["QB.F"] + f["QC.F"] + f["C3.F"]
                    + f["Q2F.F"] + f["Q2FF.F"] + f["F2e.F"] + f["K2F2.F"]
                    + f["ii.X"],
                    f["FF.in"] + f["k2F.Fin"] + f["out.F"]),
        "heat": (f["C1.heat"] + f["harp.heat"] + f["harpP.heat"]
                 + f["mix.heat"] + f["G4.heat"]
                 + f["G4A.heat"] + f["G4C.heat"] + f["QB.heat"]
                 + f["QC.heat"] + f["tri.heat"] + f["C3.heat"]
                 + f["Q2F.heat"] + f["Q2FF.heat"] + f["e.cap"]
                 + f["F2e.heat"] + f["K15.heat"] + f["R11p.heat"]
                 + f["ii.heat"] + f["FF.heat"] + f["AK8.heat"]
                 + f["RAD.heat"] + f["K2F2.heat"] + f["k2F.heat"]
                 + f["KPS.heat"] + f["PRAD.heat"],
                 f["out.heat"]),
        "248 nm": (f["G4.hv"], f["out.248"]),
        "272 nm": (f["G4A.hv"], f["out.272"]),
        "C-A band": (f["G4C.hv"], f["out.275"]),
        "420 nm": (f["C3.hv"], f["out.420"]),
        "147 nm": (f["RAD.hv"], f["out.147"]),
        "near IR": (f["PRAD.hv"], f["out.IR"]),
        "other stores": (f["krs.left"] + f["kpp.left"] + f["k2s.left"]
                         + f["B.left"] + f["C.left"] + f["k2f.left"]
                         + f["e.dU"],
                         f["out.other"]),
        "ledger": (100.0, f["out.248"] + f["out.272"] + f["out.275"]
                   + f["out.420"] + f["out.147"] + f["out.IR"]
                   + f["out.heat"] + f["out.elastic"] + f["out.ions"]
                   + f["out.F"] + f["out.other"]),
    }


def _rate(res, rid) -> float:
    """The rate one row ran with, as compiled (a band row's resolved
    value, a radiative row's share of its lifetime)."""
    c = next(c for c in res.model.compiled if c.row.rid == rid)
    if c.kind not in ("const", "param", "inv_tau"):
        raise ValueError(f"{rid}: rate kind {c.kind!r} has no single value")
    return float(c.rate) * c.factor


def _band_value(b, v) -> str:
    """A band end as the map prints it (math)."""
    if b.unit in ("exponent", "fraction", "factor"):
        return f"{v:.2f}" if b.unit == "exponent" else f"{v:g}"
    return "0" if v == 0.0 else _k(v)


# the bands in LaTeX (the documents' macros), as the review labels them
BAND_LABELS = {
    "kr2_f": r"\Krts{} + F \ra KrF(B) + Kr",
    "kr2_pen": r"\Krts{} + \Krts{} \ra \Krtp{} + 2Kr + \el{}",
    "kr2f_f": r"\KtF{} + F \ra 2Kr + 2F",
    "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr, 300\unit{K}",
    "kr2_f2": r"\Krts{} + \Ft{}, both channels",
    "kr_pen": r"\Krs{} + \Krs{} \ra \Krp{} + Kr + \el{}",
    "kr2f_f2": r"\KtF{} + \Ft{}",
    "trimer_tg_pow": r"the trimer's $(T_g/300)^n$",
    "krf_f": r"KrF(B,C) + F \ra Kr + 2F",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "k_5s5p": r"\el{} + \Krs{} \ra \el{} + \Krpp{}, pre-exponential",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}, pre-exponential",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}, pre-exponential",
    "k5p_f2": r"\Krpp{} + \Ft{}, all channels",
    "ii_yield": r"KrF\ensuremath{^*} per ion-ion recombination",
    "ii_f_c": r"share of the ion-ion KrF\ensuremath{^*} in C",
    "ib_scale": r"the wave's electron-ion absorption, a factor",
}
UNITS = {"cm3/s": r"\unit{cm^3\,s^{-1}}", "cm6/s": r"\unit{cm^6\,s^{-1}}",
         "exponent": "", "fraction": "", "factor": ""}


def bands_numbers(p: Params) -> tuple:
    """(numbers, table): the 248-nm share of E_dep at the two ends of every
    band, one band at a time (run_v0_4.py --bands), for the map's box, and
    the same with the 147-nm share, the gain, the absorption and the final
    gas temperature as the rows of the overview's table (\\wmBands)."""
    rows = do_bands(p)
    ref = rows[0]
    by = {(r["band"], r["end"]): r for r in rows[1:]}
    if set(BAND_ORDER) != set(BANDS):
        raise AssertionError("the bands changed: update BAND_ORDER and the "
                             "map's box")
    n, lines, minor = {}, [], 0.0
    for name in BAND_ORDER:
        b = BANDS[name]
        lo, hi = by[(name, "low")], by[(name, "high")]
        n[f"bd.{name}.def"] = _band_value(b, b.default)
        n[f"bd.{name}.lo"] = _band_value(b, b.lo)
        n[f"bd.{name}.hi"] = _band_value(b, b.hi)
        n[f"bd.{name}.lo248"] = f"{lo['emission_248']:.2f}"
        n[f"bd.{name}.hi248"] = f"{hi['emission_248']:.2f}"
        span = max(abs(lo["emission_248"] - ref["emission_248"]),
                   abs(hi["emission_248"] - ref["emission_248"]))
        if name not in ("kr2_f", "trimer_tg_pow", "krf_f", "k_5s5p"):
            minor = max(minor, span)
        lines.append(
            f"{BAND_LABELS[name]} & ${_band_value(b, b.default)}$ & "
            f"${_band_value(b, b.lo)}$ to ${_band_value(b, b.hi)}$"
            f"{UNITS[b.unit]} & "
            f"{lo['emission_248']:.2f} to {hi['emission_248']:.2f} & "
            f"{lo['emission_vuv']:.1f} to {hi['emission_vuv']:.1f} & "
            f"{lo['g0_pct_cm']:.2f} to {hi['g0_pct_cm']:.2f} & "
            f"{lo['alpha_pct_cm']:.2f} to {hi['alpha_pct_cm']:.2f} & "
            f"{lo['tg_end_K']:.0f} to {hi['tg_end_K']:.0f}\\\\")
    n["bd.minor"] = _f2(minor)
    n["bd.n"] = str(len(BANDS))
    n["bd.nzero"] = str(sum(1 for b in BANDS.values() if b.default == 0.0))
    return n, "\\def\\wmBands{" + "\n".join(lines) + "}"


def _f2(x: float) -> str:
    return f"{x:.2f}"


def numbers(params: Params = None) -> tuple:
    """(numbers, tables): the formatted strings the map prints, keyed as
    in the .tex, and the overview's bands table."""
    f, x = flows(params)
    p, res, rows, ev = x["params"], x["res"], x["rows"], x["events"]
    n = {k: _f1(v) for k, v in f.items()}
    # the G4 box draws B -> X and B -> A together: its tags add the two
    n["G4x.F"] = _f1(f["G4.F"] + f["G4A.F"])
    n["G4x.heat"] = _f1(f["G4.heat"] + f["G4A.heat"])
    # the ion-ion band carries the KrF(B) of Kr2* + F2 above the point
    # where the two meet
    n["iiNK7.B"] = _f1(f["ii.B"] + f["NK7.B"])
    # the Kr2* fates box draws KR2_RAD and Kr2* + F2 together: its input
    # band and its heat tag add the two
    n["k2s.fates"] = _f1(f["k2s.RAD"] + f["k2s.F2"])
    n["KRx.heat"] = _f1(f["RAD.heat"] + f["K2F2.heat"])
    # the harpoon box takes Kr* and Kr(5p): its tags add the two
    n["harpx.F"] = _f1(f["harp.F"] + f["harpP.F"])
    n["harpx.heat"] = _f1(f["harp.heat"] + f["harpP.heat"])
    n["harpx.B"] = _f1(f["harp.B"] + f["harpP.B"])
    n["harpx.C"] = _f1(f["harp.C"] + f["harpP.C"])
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
    k_h = _rate(res, "G3") + _rate(res, "G3C") + _rate(res, "G3X")
    y_h = (_rate(res, "G3") + _rate(res, "G3C")) / k_h
    f_c = _rate(res, "G3C") / (_rate(res, "G3") + _rate(res, "G3C"))
    k_hp = _rate(res, "K27") + _rate(res, "K27C") + _rate(res, "K27X")
    k_eq = rows["CB"].rate[2]
    formed_ii = sum(ev[r] for r in ION_ION_B + ION_ION_C)
    formed_k2 = ev["NK7"] + ev["K2F_F"]
    formed_p = ev["K27"] + ev["K27C"]
    formed_all = ev["G3"] + ev["G3C"] + formed_ii + formed_k2 + formed_p
    if abs(formed_all - krf["formed"]) > 1e-9 * krf["formed"]:
        raise AssertionError("the KrF sources on the map miss one")
    kind, k0, te0, n_pow = rows["R11p"].rate
    assert kind == "te_pow", rows["R11p"].rate
    k_f2 = _rate(res, "NK7") + _rate(res, "K2F_BR")
    burn = tb if tb is not None else p.pulse_s
    before = res.t <= burn
    i_b = int(res.t.searchsorted(burn))
    kr_in = {rid: float(res.events[rid][i_b]) for rid in
             ("AK8_Kr",) + HARPOON}
    kr_in["G2"] = 2.0 * float(res.events["G2"][i_b])
    kr_in["K76"] = float(res.events["K76"][i_b]) - float(
        res.events["K76r"][i_b])
    a76 = {c.row.rid: c for c in m.compiled}["K76"].rate
    kstep_s = {c.row.rid: c for c in m.compiled}["KSTEP_S"].rate
    w = optics.window(res)
    byab = sorted(w["by_species"].items(), key=lambda kv: -kv[1])
    n.update({
        "version": __version__,
        "p": f"{p.pressure_atm:g}",
        "xF2": f"{100 * p.x_f2:g}",
        "Tg0": f"{p.tg0_K:.0f}",
        "P": f"{p.power_W_cm3 / 1e3:.0f}",
        "pulse": f"{p.pulse_s * 1e6:g}",
        "tend": f"{p.t_end_s * 1e6:g}",
        "ne": _pow10(p.ne_cm3),
        "neF2": _sci(float(res.ne[before].mean())),
        "nePeak": _sci(float(res.ne.max())),
        "N": _sci(m.N0),
        "Edep": f"{res.ledger['E_dep'][-1] * 1.602176634e-19:.2f}",
        "Te": f"{float(res.te[on].mean()):.2f}",
        "TeMid": f"{float(res.te[k]):.2f}",
        "rTail": f"{r_tail:.2f}",
        "TtMid": f"{r_tail * float(res.te[k]):.2f}",
        "TgEnd": f"{res.tg[-1]:.0f}",
        "TgMax": f"{res.tg.max():.0f}",
        "pEnd": f"{res.pressure_atm[-1]:.2f}",
        "F2burn": f"{tb * 1e9:.0f}" if tb is not None else "none",
        "F2left": f"{100 * res.n['F2'][-1] / res.n['F2'][0]:.2f}",
        "krsFormed": _sci(x["krs"]["formed"]),
        "kppFormed": _sci(x["kpp"]["formed"]),
        "k2sFormed": _sci(x["k2s"]["formed"]),
        "krfFormed": _sci(krf["formed"]),
        "k2fFormed": _sci(x["k2f"]["formed"]),
        "krf.eff": _f1(100.0 * krf["formed"] * rows["G4"].hv_eV
                       / float(res.ledger["E_dep"][-1])),
        # where KrF(B,C) comes from, percent of the KrF formed
        "frm.ii": _f1(100.0 * formed_ii / formed_all),
        "frm.kr2": _f1(100.0 * formed_k2 / formed_all),
        "frm.5p": _f1(100.0 * formed_p / formed_all),
        # where KrF(B,C) goes, percent of the KrF formed (B and C together)
        "brf.G4": brf("G4"), "brf.G4A": brf("G4A"), "brf.G4C": brf("G4C"),
        "brf.tri": brf("G5", "G5C"), "brf.e": brf("QE", "QEC"),
        "brf.F2": brf("QF2", "QF2C"), "brf.Kr": brf("QKR", "QKRC"),
        "brf.F": brf("QF", "QFC"),
        "brf.QB": brf(*QUENCH_B), "brf.QC": brf(*QUENCH_C),
        # how the ions recombine: with an electron (R11p) or with F-
        "ion.R11p": _f1(100.0 * ev["R11p"] / (ev["R11p"] + formed_ii)),
        "ion.ii": _f1(100.0 * formed_ii / (ev["R11p"] + formed_ii)),
        # Kr* into Kr2* and up to 5p before the F2 burns out, percent of
        # the Kr* used up to then
        "brs.AK8.F2": _f1(100.0 * kr_in["AK8_Kr"] / sum(kr_in.values())),
        "brs.76.F2": _f1(100.0 * kr_in["K76"] / sum(kr_in.values())),
        # the 5p rung: gross cycle events over the Kr* formed, and the
        # cycle's share of the input
        "cyc.gross": _f1(ev["K76"] / x["krs"]["formed"]),
        "cyc.e": _f1(f["e.76"]),
        # the constants the boxes print, read from the rows the model ran
        "k.harp": _k(k_h),
        "k.harpP": _k(k_hp),
        "pc.harpY": f"{100 * y_h:.0f}",
        "pc.harpB": f"{100 * y_h * (1 - f_c):.0f}",
        "pc.harpC": f"{100 * y_h * f_c:.0f}",
        "pc.harpX": f"{100 * (1 - y_h):.0f}",
        "pc.fC": f"{100 * f_c:.0f}",
        "k.G2": _k(_rate(res, "G2")),
        "tauB": f"{p.tau_sp_s * 1e9:.2f}",
        "pc.BA": f"{100 * KRF_BA_OVER_BX:.1f}",
        "tauC": _sig2(1e9 / _rate(res, "G4C")),
        "tauC3": _sig2(1e9 / _rate(res, "C3")),
        "gapBC": f"{(e_form('KrF(B)') - e_form('KrF(C)')) / 1.23984e-4:.0f}",
        "k.BC": _k(_rate(res, "BC")),
        "Keq": f"{k_eq:.2f}",
        "k.G5": _k(_rate(res, "G5")),
        "k.G5C": _k(_rate(res, "G5C")),
        "tpow.G5": f"{rows['G5'].tg_pow:.2f}",
        "k.Q2F": _k(_rate(res, "Q2F")),
        "k.K15": _k(_rate(res, "K15")),
        "tpow.K15": f"{rows['K15'].tg_pow:.2f}",
        "k.R11p": _k(k0),
        "pow.R11p": f"{-n_pow:.2f}",
        "cap.R11p": f"{1.5 + n_pow:.2f}",
        "k.FF": _k(_rate(res, "FF")),
        "n3.ii": f"{rows['K26'].rate[3]:g}",
        "k3.ii": _k(rows["K26"].rate[1]),
        "kL.ii": _k(rows["K26"].rate[2]),
        # the Kr2* route
        "k.AK8": _k(_rate(res, "AK8_Kr")),
        "tpow.AK8": f"{rows['AK8_Kr'].tg_pow:.2f}",
        "tauRAD": _sig2(1e9 / _rate(res, "KR2_RAD")),
        "k.kr2F2": _k(k_f2),
        "pc.NK7": f"{100 * _rate(res, 'NK7') / k_f2:.0f}",
        "E.kr2s": f"{e_form('Kr2*'):.2f}",
        "E.kr2s.bind": f"{e_form('Kr*') - e_form('Kr2*'):.3f}",
        "hv.RAD": f"{rows['KR2_RAD'].hv_eV:.2f}",
        # the 5p rung's constants
        "E.PP": f"{e_form('Kr**'):.2f}",
        "gap5p": f"{e_form('Kr**') - e_form('Kr*'):.3f}",
        "k.KPS": _k(_rate(res, "KPS")),
        "tauPRAD": _sig2(KR_5P_TAU_S * 1e9),
        "A.76": _k(a76[0]),
        "B.76": f"{a76[1]:.2f}",
        "C.76": f"{a76[2]:.2f}",
        "C.stepS": f"{kstep_s[2]:.2f}",
        "hv.PRAD": f"{rows['PRAD'].hv_eV:.3f}",
        # how far F- lies below Kr + F2 on the ledger's scale (a magnitude)
        "E.anion": f"{-e_form('F-'):.1f}",
        "hv.G4C": f"{rows['G4C'].hv_eV:.2f}",
        "hv.G4": f"{rows['G4'].hv_eV:.2f}",
        "hv.G4A": f"{rows['G4A'].hv_eV:.2f}",
        "hv.C3": f"{rows['C3'].hv_eV:.2f}",
        "E.B": f"{e_form('KrF(B)'):.3f}",
        # the notes' output line (flows, they overlap)
        "notes.krs": _f1(f["e.G1"] + f["e.C1"]),
        "notes.krp": _f1(f["C2.ion"] + f["G2.ion"] + f["pen2.ion"]
                         + f["stepS.ion"] + f["stepP.ion"]),
        "notes.emit": _f1(f["out.248"] + f["out.272"] + f["out.275"]
                          + f["out.420"] + f["out.147"] + f["out.IR"]),
        "nodes": str(len(_nodes(f))),
        # the 248-nm optics over the pulse's second half (optics.py)
        "op.t0": f"{w['t0'] * 1e9:.0f}", "op.t1": f"{w['t1'] * 1e9:.0f}",
        "op.g0": f"{w['g0']:.2f}", "op.alpha": f"{w['alpha']:.2f}",
        "op.net": f"{w['net']:.2f}",
        "op.g0peak": f"{w['g0_peak']:.2f}",
        "op.tpeak": f"{w['t_peak'] * 1e9:.0f}",
        "op.by": ", ".join(f"{_sp(s)} {v:.2f}" for s, v in byab if v > 0.0),
        "op.sigB": _k(optics.SIGMA_SE_CM2),
        "op.sigK2F": _k(optics.ABSORBERS["Kr2F*"][0]),
    })
    # what the basic network leaves out: the 248-nm share with the electron
    # group on (the only group v0_4 has), the notes' list as written, and
    # the 5p rung neutralised (the 5s <-> 5p cycle off, the rung relaxed
    # at once: v0_3's reference point)
    for g in BASIC_GROUPS:
        d = row_dict(run(p.replace(groups=(g,))))
        n[f"om.{g}"] = f"{d['emission_248']:.1f}"
    n["om.basic"] = n["out.248"]
    n["om.v0"] = f"{row_dict(run(p.replace(network='v0')))['emission_248']:.1f}"
    n["om.v03"] = f"{row_dict(run(p.replace(rates=(('k_5s5p', 0.0), ('kps_kr', 1e-10)))))['emission_248']:.2f}"
    n["om.kannari"] = f"{row_dict(run(p.replace(ladder='kannari')))['emission_248']:.2f}"
    # the bands, each end alone
    bn, table = bands_numbers(p)
    n.update(bn)
    return n, [table]


def _sp(s: str) -> str:
    """A species in the map's LaTeX macros."""
    return {"F-": r"\Fm{}", "Kr2+": r"\Krtp{}", "Kr**": r"\Krpp{}",
            "F2": r"\Ft{}", "Kr2F*": r"\KtF{}", "Kr*": r"\Krs{}"}[s]


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT v0_4 map numbers. GENERATED by make_system_map_v0_4.py",
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
    ap = argparse.ArgumentParser(description="v0_4's system map numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh run, write nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_system_map_v0_4.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
