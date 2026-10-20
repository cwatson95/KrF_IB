"""The energy ladder of LIGHT_ArF_ebeam v0_2 (energy_ladder.tex), for its overview.

LIGHT's make_energy_ladder.py drawn for the ArF network: every species of
v0_2 stands at its formation energy, the energy above ground-state Ar and F2
at rest (species.py), and the rows that move energy between the species are
arrows labelled with their row ids (network.py). 1Dv0_2 runs v0_2's
chemistry unchanged, so the one ladder is both models'.

Every number on the ladder comes from the code. A level's height and label
are its E_form, and a pair's (dashed) is the sum of its two. An
electron-impact threshold is its row's e_cost_eV (the cross section's
threshold), a photon is its row's hv_eV and band, and the heat a photon
leaves is its row's energy_routing. The layout, which column a level sits
in, the x of every arrow and where each label goes, is this script's. An
arrow's style says where its row stands: solid in the default network,
dashed in the optional group, dotted for a band whose default is 0. Every
row the network can run is drawn, or named in NOT_DRAWN, and the script
refuses a row it has no drawing for. The 37 argon excitation rows are drawn
as two family arrows, the four 4s rows (G1) and the 33 higher ones (C1),
each with its lowest and highest threshold.

ArF(B) is the one excimer state (no C state is carried), and it sits at the
spectroscopic lower bound, exactly one 193-nm quantum above Ar + F, so the
193-nm photon leaves no lower-state heat.

    python3 LIGHT_ArF_ebeam/docs/make_energy_ladder.py          # writes the file
    python3 LIGHT_ArF_ebeam/docs/make_energy_ladder.py --check  # exit 1 when stale

The file defines its numbers (\\laddef{key}{value}) and the picture
(\\ladpicture). The overview inputs it in its preamble and sets
\\ladpicture in a figure.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIGHT_ARF = HERE.parent
if str(LIGHT_ARF) not in sys.path:
    sys.path.insert(0, str(LIGHT_ARF))

OUT = HERE / "energy_ladder.tex"
VER = "v0_2"

# rows that leave every level where it is: elastic scattering and F2
# vibration (the electron's loss heats the gas). The caption names them.
NOT_DRAWN = ("G6", "VIB1", "VIB2", "VIB3", "VIB4")
G1_FAMILY = "G1_"          # the four 4s rows
C1_FAMILY = "C1_"          # the 33 higher rows, cascading to 4s
SUPER_FAMILY = "SUPER_"    # the optional group's superelastics

# the vertical scale: cm per eV and the energy range drawn
CM_PER_EV = 0.9
E_LO, E_HI = -3.0, 16.3
HEAD_E = E_HI + 0.6

TEX = {"Ar": "Ar", "Ar*": r"Ar\ensuremath{^*}", "Ar+": r"Ar\ensuremath{^+}",
       "Ar2*": r"Ar\ensuremath{_2^*}", "Ar2+": r"Ar\ensuremath{_2^+}",
       "ArF(B)": "ArF(B)", "Ar2F*": r"Ar\ensuremath{_2}F\ensuremath{^*}",
       "F2": r"F\ensuremath{_2}", "F": "F", "F-": r"F\ensuremath{^-}"}
# F2's two dissociative states, as the cross sections name them (sansmath
# has no upright capital Greek, so the Pi is set in the serif math font)
F2_STATE = {"a3": r"{\unsansmath\ensuremath{a\,^3\Pi_u}}",
            "A1": r"{\unsansmath\ensuremath{A\,^1\Pi_u}}"}

# --- the layout (cm) -------------------------------------------------------
X_AXIS = 0.95
AR = (2.65, 6.30)          # Ar atoms, the full width
AR2 = (5.00, 6.30)         # Ar2 dimers, the right of the Ar column
ARF = (7.00, 12.30)        # the ArF column
PAIRS = (10.20, 12.30)     # the ion pairs K26 and K108 start from
B_BAR = (8.40, 11.05)      # ArF(B)
T_BAR = (7.20, 8.75)       # Ar2F*, below ArF(B)'s left end
F2 = (13.90, 15.10)


def _neg(s: str) -> str:
    """A proper minus sign for LaTeX text."""
    return s.replace("-", r"\ensuremath{-}")


def _e(v: float) -> str:
    """An energy as the overview's formation-energy table prints it, to the
    meV and no further (ArF(B) is 0.80 + hc/193.3 nm, 7.21408...)."""
    return _neg(f"{round(v, 3):g}")


def _h(v: float) -> str:
    """A derived energy: two decimals, the precision of the inputs, rounded
    half up on its decimal value (binary floating point stores some ties
    just below them)."""
    d = Decimal(repr(round(v, 9)))
    return _neg(str(d.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)))


def _names(rids) -> str:
    return ", ".join(r.replace("_", r"\_") for r in rids)


def _band(label: str) -> tuple:
    """'193 nm (ArF B-X)' -> ('193 nm', 'X'); '284 nm (Ar2F*)' -> (.., None)."""
    m = re.match(r"(~?\d+ nm)(?: \((.*)\))?$", label)
    if m is None:
        raise ValueError(f"band label {label!r} not understood")
    wl = m.group(1).replace("~", r"\ensuremath{\sim}")
    low = re.search(r"\b[A-Z]-([A-Z])\b", m.group(2) or "")
    return wl, (low.group(1) if low else None)


def _threshold(r) -> str:
    """An F2 electron-impact row's threshold as its cross section prints it."""
    m = re.search(r"\(([\d.]+)eV\)", r.rate[3] or "")
    if m is None:
        raise ValueError(f"{r.rid}: no threshold in {r.rate!r}")
    if abs(float(m.group(1)) - r.e_cost_eV) > 1e-9:
        raise ValueError(f"{r.rid}: label {m.group(1)} is not e_cost "
                         f"{r.e_cost_eV}")
    return m.group(1)


class Version:
    """What the code says: energies, rows and their standing."""

    def __init__(self, ver: str = VER):
        self.ver = ver
        sp = importlib.import_module(f"{ver}.species")
        net = importlib.import_module(f"{ver}.network")
        self.version = importlib.import_module(ver).__version__
        self.E = {n: sp.e_form(n) for n in sp.SPECIES}
        self.default = {r.rid: r for r in net.build((), "basic")}
        every = net.build(net.BASIC_GROUPS, "basic")
        self.rows = {r.rid: r for r in every}
        self.bands = net.BANDS
        self.route = {r.rid: net.energy_routing(r) for r in every}
        made = {p for r in self.default.values() for p in r.products}
        self.optional_only = {n for n in self.E
                              if n not in made and self.E[n] != 0.0}
        self.groups = tuple(net.BASIC_GROUPS)

    def family(self, prefix: str) -> tuple:
        return tuple(r for r in self.rows if r.startswith(prefix))

    def style(self, rid: str) -> str:
        r = self.rows[rid]
        if rid not in self.default:
            return "opt"
        if r.rate[0] == "param" and self.bands[r.rate[1]].default == 0.0:
            return "zero"
        return "on"

    def has(self, rid: str) -> bool:
        return rid in self.rows

    def pair(self, *names) -> float:
        return sum(self.E[n] for n in names)


def y(e: float) -> float:
    return (e - E_LO) * CM_PER_EV


def _c(x: float, e: float) -> str:
    return f"({x:.3f},{y(e):.3f})"


class Picture:
    def __init__(self, v: Version):
        self.v = v
        self.cmds = {"grid": [], "lev": [], "arr": [], "lab": []}
        self.drawn = set()

    # levels -------------------------------------------------------------
    def level(self, x0, x1, e, color, kind="lev", pale=False, dy=0.0):
        st = {"lev": "lad lev", "pair": "lad pair", "tick": "lad tick",
              "ground": "lad lev"}[kind]
        col = f"{color}!40" if pale else color
        self.cmds["lev"].append(
            rf"\draw[{st},draw={col}] ({x0:.3f},{y(e) + dy:.3f}) -- "
            rf"({x1:.3f},{y(e) + dy:.3f});")

    def label(self, x, e, text, anchor="east", dy=0.0, style="lad lab"):
        self.cmds["lab"].append(
            rf"\node[{style},anchor={anchor}] at ({x:.3f},{y(e) + dy:.3f}) "
            rf"{{{text}}};")

    # arrows ---------------------------------------------------------------
    def _rids(self, rid) -> tuple:
        rids = (rid,) if isinstance(rid, str) else tuple(rid)
        return tuple(r for r in rids if self.v.has(r))

    def _style(self, rids) -> str:
        styles = {self.v.style(r) for r in rids}
        if len(styles) != 1:
            raise ValueError(f"rows {rids} drawn as one arrow differ in "
                             f"standing: {styles}")
        return {"on": "lad on", "opt": "lad opt", "zero": "lad zero"}[
            styles.pop()]

    def arrow(self, rid, p0, p1, color, t=0.5, side=None, lab=True,
              kind="lad arr", lstyle="lad alab"):
        """An arrow for row(s) rid (a str, or a tuple drawn as one) from
        p0=(x, E) to p1=(x, E). Its label (the row ids, or lab) sits on the
        line at fraction t, along it, or beside it when side=(anchor, dx,
        dy) gives where."""
        rids = self._rids(rid)
        if not rids:
            return None
        self.drawn |= set(rids)
        text = _names(rids) if lab is True else lab
        st = self._style(rids)
        node = ""
        if text and side is None:
            turn = "rotate=90" if p0[0] == p1[0] else "sloped"
            node = rf"node[pos={t:.2f},{turn},lad alab] {{{text}}} "
        self.cmds["arr"].append(
            rf"\draw[{kind},{st},draw={color}] {_c(*p0)} -- {node}"
            rf"{_c(*p1)};")
        if text and side is not None:
            anchor, dx, dy = side
            x = p0[0] + t * (p1[0] - p0[0]) + dx
            e = p0[1] + t * (p1[1] - p0[1])
            self.label(x, e, text, anchor=anchor, dy=dy, style=lstyle)
        return rids

    def photon(self, rid, x, e_up, color, heat_at=("west", 0.08, 0.0),
               t=0.5):
        """A radiative row: the quantum as a straight arrow labelled with the
        row and its band, then the heat its lower state leaves as a squiggle
        down to where it falls apart, labelled at heat_at=(anchor, dx, dE)
        from the photon's end."""
        if not self.v.has(rid):
            return None
        r = self.v.rows[rid]
        wl, lower = _band(r.band)
        e_low = e_up - r.hv_eV
        self.arrow(rid, (x, e_up), (x, e_low), color, t=t, kind="lad hv",
                   lab=rf"{_names((rid,))}, {wl}")
        heat = self.v.route[rid]["heat"]
        self.level(x - 0.14, x + 0.14, e_low, "ladMuted", kind="tick")
        if heat * CM_PER_EV > 0.12:
            self.squiggle(x, e_low, e_low - heat)
        low = f"{lower}, " if lower else ""
        anchor, dx, de = heat_at
        self.label(x + dx, e_low + de, rf"{low}{_h(heat)}", anchor=anchor,
                   style="lad hlab")
        return e_low

    def stub(self, rids, x, e, below=None, length=0.8):
        """Rows that return a level's energy to the gas: a short red arrow,
        labelled below its tip."""
        rids = self._rids(rids)
        if not rids:
            return
        self.drawn |= set(rids)
        st = self._style(rids)
        self.cmds["arr"].append(
            rf"\draw[lad arr,{st},draw=ladH] {_c(x, e)} -- "
            rf"{_c(x, e - length)};")
        if below:
            self.label(x, e - length - 0.05, below, anchor="north",
                       style="lad alab")

    def squiggle(self, x, e0, e1):
        self.cmds["arr"].append(
            rf"\draw[lad heat] {_c(x, e0)} -- {_c(x, e1)};")


def _family_arrow(P, v, prefix, x, t, lab, tick_x1, top_label):
    """One arrow for a family of electron-impact rows from the ground state
    to its highest threshold, ticks at its lowest and highest, and the
    cascade's heat as a squiggle down to Ar*."""
    fam = v.family(prefix)
    thr = sorted(v.rows[r].e_cost_eV for r in fam)
    lo, hi = thr[0], thr[-1]
    P.arrow(fam, (x, 0.0), (x, hi), "ladE", t=t, lab=lab)
    for e in {lo, hi}:
        P.level(x - 0.07, tick_x1, e, "ladMuted", kind="tick")
    if top_label:
        P.label(tick_x1 + 0.06, hi, top_label.format(lo=_e(lo), hi=_e(hi),
                                                     n=len(fam)),
                anchor="west", style="lad tlab")
    return lo, hi


def build(v: Version) -> Picture:
    """The ladder."""
    E = v.E
    P = Picture(v)
    pale = v.optional_only.__contains__

    # grid, axis, heads
    for k in range(-2, 17, 2):
        P.cmds["grid"].append(
            rf"\draw[lad grid] ({X_AXIS:.3f},{y(k):.3f}) -- "
            rf"({F2[1] + 1.1:.3f},{y(k):.3f});")
        P.label(X_AXIS - 0.08, k, _neg(str(k)), anchor="east",
                style="lad axis")
    P.cmds["grid"].append(
        rf"\draw[lad axisline] ({X_AXIS:.3f},{y(E_LO):.3f}) -- "
        rf"({X_AXIS:.3f},{y(E_HI):.3f});")
    P.label(X_AXIS - 0.08, HEAD_E, "eV", anchor="east", style="lad axis")
    for (x0, x1), head in ((AR, r"Ar and Ar\ensuremath{_2}"),
                           (ARF, r"ArF and Ar\ensuremath{_2}F"),
                           (F2, r"F\ensuremath{_2}")):
        P.label(0.5 * (x0 + x1), HEAD_E, head, anchor="center",
                style="lad head")

    # --- Ar and Ar2 --------------------------------------------------------
    eX, eP, eV, eI = E["Ar*"], E["Ar+"], E["Ar2*"], E["Ar2+"]
    P.level(*AR, 0.0, "ladInk2", kind="ground")
    P.label(AR[0] - 0.08, 0.0, "Ar")
    P.level(*AR, eX, "ladX")
    P.label(AR[0] - 0.08, eX, rf"{TEX['Ar*']} {_e(eX)}")
    P.level(*AR, eP, "ladI")
    P.label(AR[0] - 0.08, eP, rf"{TEX['Ar+']} {_e(eP)}")
    P.level(*AR2, eV, "ladV", pale=pale("Ar2*"))
    P.label(AR2[0] + 0.03, eV, rf"{TEX['Ar2*']}\\{_e(eV)}",
            anchor="north west", dy=-0.05, style="lad lab2l")
    P.level(*AR2, eI, "ladI", pale=pale("Ar2+"))
    P.label(AR2[1], eI, rf"{TEX['Ar2+']}\\{_e(eI)}", anchor="south east",
            dy=0.05, style="lad lab2")

    # electron impact from the ground state, and Ar* up and down
    P.arrow("C2", (2.80, 0.0), (2.80, eP), "ladE", t=0.52)
    lo4, hi4 = _family_arrow(P, v, G1_FAMILY, 3.08, 0.40, "G1 (4s)",
                             3.30, "")
    P.squiggle(3.30, hi4, eX)
    P.arrow(v.family(SUPER_FAMILY), (3.36, eX), (3.36, 0.0), "ladE",
            t=0.78, lab="SUPER")
    P.arrow("G2", (3.08, eX), (3.08, eP), "ladX", t=0.36)
    P.arrow("KSTEP", (3.36, eX), (3.36, eP), "ladE", t=0.62)
    lo, hi = _family_arrow(P, v, C1_FAMILY, 3.70, 0.55, "C1", 4.40,
                           r"{hi}")
    P.label(4.46, lo, _e(lo), anchor="west", style="lad tlab")
    P.squiggle(4.40, hi, eX)
    # G3X: Ar* + F2 -> Ar + 2F, a quench, leaves from Ar*'s left end
    if v.has("G3X"):
        P.drawn.add("G3X")
        P.cmds["arr"].append(
            rf"\draw[lad arr,{P._style(('G3X',))},draw=ladH] "
            rf"{_c(AR[0] + 0.05, eX)} -- {_c(AR[0] - 0.30, eX - 0.85)};")
        P.label(AR[0] - 0.34, eX - 0.72, "G3X", anchor="east",
                style="lad alab")
    # the dimers
    P.arrow("AK8_Ar", (5.10, eX), (5.10, eV), "ladX", t=0.5,
            side=("east", -0.06, 0.0))
    P.arrow("K15", (5.25, eP), (5.25, eI), "ladI", t=0.5)
    P.arrow("R11p", (5.80, eI), (5.80, eX), "ladI", t=0.45)
    P.arrow("PEN_AR2", (6.12, eV), (6.12, eI), "ladV", t=0.62)
    P.photon("AR2_RAD", 5.60, eV, "ladV", heat_at=("west", 0.10, -0.55))

    # --- ArF and Ar2F ------------------------------------------------------
    eB, eT = E["ArF(B)"], E["Ar2F*"]
    e_ip, e_ip2 = v.pair("Ar+", "F-"), v.pair("Ar2+", "F-")
    e_af = v.pair("Ar", "F")
    P.level(*PAIRS, e_ip, "ladInk2", kind="pair", pale=pale("F-"))
    P.label(PAIRS[1] + 0.08, e_ip,
            rf"{TEX['Ar+']} + {TEX['F-']} {_e(e_ip)}", anchor="west")
    P.level(*PAIRS, e_ip2, "ladInk2", kind="pair",
            pale=pale("Ar2+") or pale("F-"))
    P.label(PAIRS[1] + 0.08, e_ip2,
            rf"{TEX['Ar2+']} + {TEX['F-']} {_e(e_ip2)}", anchor="west")
    P.level(ARF[0], PAIRS[1], e_af, "ladInk2", kind="pair")
    P.label(PAIRS[1], e_af, rf"Ar + F {_e(e_af)}", anchor="south east",
            dy=0.04)
    P.level(*B_BAR, eB, "ladB")
    P.label(B_BAR[1] + 0.08, eB, rf"{TEX['ArF(B)']} {_e(eB)}\\"
            r"one state, no C carried", anchor="west", style="lad lab2l")
    P.level(*T_BAR, eT, "ladT")
    P.label(T_BAR[0] - 0.08, eT, rf"{TEX['Ar2F*']}\\{_e(eT)}",
            style="lad lab2")

    # formation: the lower source (Ar2*) lands left of the upper (Ar*),
    # and the ion pairs drop to the right of both, so nothing crosses
    P.arrow("K2F_BR", (AR2[1], eV), (7.52, eT), "ladV", t=0.42)
    P.arrow("NK7", (AR2[1], eV), (8.66, eB), "ladV", t=0.42)
    P.arrow("K2F_F", (AR2[1], eV), (8.90, eB), "ladV", t=0.62)
    P.arrow("G3", (AR[1], eX), (9.14, eB), "ladX", t=0.40)
    P.arrow("K108", (10.45, e_ip2), (10.45, eB), "ladI", t=0.5)
    P.arrow("K26", (10.72, e_ip), (10.72, eB), "ladI", t=0.62)
    # the trimer, the 193-nm light and B's quenchers
    P.arrow("G5", (8.62, eB), (8.62, eT), "ladB", t=0.5,
            side=("west", 0.05, 0.0), lstyle="lad slab")
    P.photon("G4", 9.55, eB, "ladB", heat_at=("north", 0.0, -0.10))
    P.stub(("QF2", "QAR", "QE"), 10.93, eB, length=0.8)
    P.label(10.62, eB - 0.85, _names(("QF2", "QAR", "QE")),
            anchor="north west", style="lad alab")
    # Ar2F*: its light and its quenchers (each photon leaves between two
    # arrivals, so no arrow seems to run through the bar)
    P.photon("C3", 7.78, eT, "ladT", heat_at=("west", 0.10, -0.55), t=0.64)
    P.stub(("Q2F",), 7.32, eT, below="Q2F")
    P.stub(("Q2FF",), 8.14, eT, below="Q2FF")
    P.stub(("Q2AR",), 8.46, eT, length=0.5)
    P.label(8.52, eT - 0.55, "Q2AR", anchor="north west", style="lad alab")

    # --- F2 ------------------------------------------------------------------
    eF, eFm = E["F"], E["F-"]
    e_ff, e_fmf = v.pair("F", "F"), v.pair("F-", "F")
    P.level(*F2, 0.0, "ladInk2", kind="ground")
    P.label(F2[0] - 0.08, 0.0, TEX["F2"])
    P.level(*F2, eF, "ladF")
    P.label(F2[0] - 0.08, eF, rf"F {_e(eF)}")
    P.level(*F2, e_ff, "ladInk2", kind="pair")
    P.label(F2[0] - 0.08, e_ff, rf"F + F {_e(e_ff)}")
    P.level(*F2, e_fmf, "ladInk2", kind="pair", pale=pale("F-"))
    P.label(F2[0] - 0.08, e_fmf, rf"{TEX['F-']} + F {_e(e_fmf)}")
    P.level(*F2, eFm, "ladI", pale=pale("F-"))
    P.label(F2[0] - 0.08, eFm, rf"{TEX['F-']} {_e(eFm)}")
    dis = [r for r in ("DIS_a3", "DIS_A1") if v.has(r)]
    top = e_ff
    for r, x, t in zip(dis, (14.36, 14.64), (0.62, 0.66)):
        e = v.rows[r].e_cost_eV
        top = max(top, e)
        state = re.search(r"F2\((\w+)\)", v.rows[r].rate[3]).group(1)
        P.arrow(r, (x, 0.0), (x, e), "ladE", t=t)
        P.level(x - 0.07, F2[1], e, "ladMuted", kind="tick")
        P.label(F2[1] + 0.05, e, rf"{F2_STATE[state]} {_threshold(v.rows[r])}",
                anchor="west", style="lad tlab")
    if dis:
        P.squiggle(F2[1], top, e_ff)
    P.arrow("ATT", (14.08, 0.0), (14.08, e_fmf), "ladE", t=0.5)
    P.arrow("FF", (14.88, e_ff), (14.88, 0.0), "ladF", t=0.25)

    missing = set(v.rows) - P.drawn - set(NOT_DRAWN)
    if missing:
        raise ValueError(f"rows with no drawing: {sorted(missing)}")
    return P


def _and(items) -> str:
    """'A', 'A and B', 'A, B and C'."""
    items = list(items)
    if len(items) < 2:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def numbers(v: Version) -> dict:
    """The numbers and row lists the overview cites next to the ladder."""
    E, R = v.E, v.route

    def ids(rows):
        return _and(r.replace("_", r"\_") for r in rows)

    def wl(rid):
        return _band(v.rows[rid].band)[0].replace(" nm", "")

    g1 = sorted(v.rows[r].e_cost_eV for r in v.family(G1_FAMILY))
    c1 = sorted(v.rows[r].e_cost_eV for r in v.family(C1_FAMILY))
    return {
        "version": v.version.replace("_", r"\_"),
        "E.X": _e(E["Ar*"]), "E.P": _e(E["Ar+"]),
        "E.B": _e(E["ArF(B)"]), "E.T": _e(E["Ar2F*"]),
        "E.V": _e(E["Ar2*"]), "E.I": _e(E["Ar2+"]),
        "DeB": _h(v.pair("Ar+", "F-") - E["ArF(B)"]),
        "BoverArF": _h(E["ArF(B)"] - v.pair("Ar", "F")),
        "wl.G4": wl("G4"), "hv.G4": _h(v.rows["G4"].hv_eV),
        "heat.G4": _h(R["G4"]["heat"]),
        "heat.G3": _h(R["G3"]["heat"]),
        "heat.NK7": _h(R["NK7"]["heat"]),
        "heat.K26": _h(R["K26"]["heat"]),
        "heat.K108": _h(R["K108"]["heat"]),
        "heat.G5": _h(R["G5"]["heat"]),
        "wl.C3": wl("C3"), "hv.C3": _h(v.rows["C3"].hv_eV),
        "heat.C3": _h(R["C3"]["heat"]),
        "wl.V": wl("AR2_RAD"), "hv.V": _h(v.rows["AR2_RAD"].hv_eV),
        "heat.V": _h(R["AR2_RAD"]["heat"]),
        "nG1": str(len(g1)), "nC1": str(len(c1)),
        "g1.lo": _e(g1[0]), "g1.hi": _e(g1[-1]),
        "c1.lo": _e(c1[0]), "c1.hi": _e(c1[-1]),
        "notdrawn": ids(r for r in NOT_DRAWN if v.has(r)),
        "toe": ids(sorted({("SUPER" if r.startswith(SUPER_FAMILY) else r)
                           for r in v.rows if R[r]["to_e"] > 0})),
        "zero": ids(r for r in v.rows if v.style(r) == "zero"),
        "groups": _and(rf"\texttt{{{g}}}" for g in v.groups),
    }


PREAMBLE = (r"\providecommand\laddef[2]{\expandafter\def"
            r"\csname lad@#1\endcsname{#2}}")

# LIGHT's ladder palette (checked with the dataviz validator for its map);
# the ArF ladder uses the KrF(B) hue for ArF(B) and draws no C state
STYLES = r"""\definecolor{ladE}{HTML}{2A78D6}% electron gas
\definecolor{ladX}{HTML}{EDA100}% Ar*
\definecolor{ladB}{HTML}{1BAF7A}% ArF(B), 193 nm
\definecolor{ladT}{HTML}{4A3AA7}% Ar2F*, 284 nm
\definecolor{ladV}{HTML}{A63A83}% Ar2*, 127 nm
\definecolor{ladH}{HTML}{E34948}% heat
\definecolor{ladI}{HTML}{52514E}% ions
\definecolor{ladF}{HTML}{A8A69E}% F atoms
\definecolor{ladInk}{HTML}{0B0B0B}%
\definecolor{ladInk2}{HTML}{52514E}%
\definecolor{ladMuted}{HTML}{898781}%
\tikzset{
  lad grid/.style={draw=ladMuted!22, line width=0.3pt},
  lad axisline/.style={draw=ladMuted, line width=0.5pt},
  lad lev/.style={line width=1.4pt},
  lad pair/.style={line width=0.8pt, dash pattern=on 3pt off 2pt},
  lad tick/.style={line width=0.5pt},
  lad arr/.style={line width=0.7pt, -{Triangle[width=3.2pt,length=3.4pt]}},
  lad hv/.style={line width=1.3pt, -{Triangle[width=4.2pt,length=4.4pt]}},
  lad on/.style={},
  lad opt/.style={dash pattern=on 2.4pt off 1.6pt},
  lad zero/.style={dash pattern=on 0.6pt off 1.4pt, line cap=round},
  lad heat/.style={draw=ladH, line width=0.6pt, decorate,
    decoration={snake, amplitude=0.5mm, segment length=1.4mm,
    post length=0.5mm}},
  lad lab/.style={font=\scriptsize, text=ladInk, inner sep=1pt},
  lad lab2/.style={font=\scriptsize, text=ladInk, inner sep=1pt,
    align=right},
  lad lab2l/.style={font=\scriptsize, text=ladInk, inner sep=1pt,
    align=left},
  lad alab/.style={font=\scriptsize, text=ladInk, inner sep=0.8pt,
    fill=white, fill opacity=0.85, text opacity=1, align=center},
  lad slab/.style={font=\tiny, text=ladInk, inner sep=0.6pt,
    fill=white, fill opacity=0.85, text opacity=1, align=left},
  lad hlab/.style={font=\tiny, text=ladInk2, inner sep=0.8pt,
    fill=white, fill opacity=0.85, text opacity=1},
  lad tlab/.style={font=\tiny, text=ladInk2, inner sep=0.8pt},
  lad axis/.style={font=\scriptsize, text=ladInk2, inner sep=1pt},
  lad head/.style={font=\footnotesize\bfseries, text=ladInk2},
}"""


def legend(v: Version) -> list:
    """The key, two columns below the ladder."""
    items = [("lad arr,lad on,draw=ladInk2", "a row of the default network")]
    if any(v.style(r) == "opt" for r in v.rows):
        items.append(("lad arr,lad opt,draw=ladInk2",
                      "a row of the optional group, off by default"))
    if any(v.style(r) == "zero" for r in v.rows):
        items.append(("lad arr,lad zero,draw=ladInk2",
                      "a band whose default is 0, no measurement"))
    items += [("lad hv,lad on,draw=ladB", "a photon, with its row and band"),
              ("lad arr,lad on,draw=ladH",
               "a quench, the level's energy to heat"),
              ("lad heat",
               "heat as a lower state falls apart (eV) or cascades"),
              ("lad pair,draw=ladInk2",
               "a pair, the sum of two formation energies"),
              ("lad tick,draw=ladMuted",
               "a threshold, or where a photon leaves its lower state")]
    if v.optional_only:
        items.append(("lad lev,draw=ladInk2!40",
                      "a species that only the optional group makes"))
    out = []
    x = X_AXIS + 0.3
    col_w = 7.6
    y0 = y(E_LO) - 0.55
    for i, (st, text) in enumerate(items):
        xi = x + (i % 2) * col_w
        yi = y0 - (i // 2) * 0.42
        out.append(rf"\draw[{st}] ({xi:.3f},{yi:.3f}) -- "
                   rf"({xi + 0.7:.3f},{yi:.3f});")
        out.append(rf"\node[lad lab,anchor=west] at ({xi + 0.8:.3f},"
                   rf"{yi:.3f}) {{{text}}};")
    return out


def render() -> str:
    v = Version()
    P = build(v)
    head = ("% LIGHT_ArF_ebeam v0_2 energy ladder. GENERATED by make_energy_ladder.py\n"
            f"% from LIGHT_ArF_ebeam v0_2 ({v.version}). Do not edit. Rerun the "
            "script instead.\n")
    lines = [head + PREAMBLE]
    for k, val in numbers(v).items():
        lines.append(rf"\laddef{{{k}}}{{{val}}}")
    lines.append(r"\def\ladpicture{%")
    lines.append(STYLES)
    lines.append(r"\begin{tikzpicture}")
    for part in ("grid", "lev", "arr", "lab"):
        lines += P.cmds[part]
    lines += legend(v)
    lines.append(r"\end{tikzpicture}}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if energy_ladder.tex is stale")
    a = ap.parse_args(argv)
    text = render()
    if a.check:
        if not OUT.is_file() or OUT.read_text() != text:
            print(f"{OUT.name} is stale. Rerun make_energy_ladder.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(text)
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
