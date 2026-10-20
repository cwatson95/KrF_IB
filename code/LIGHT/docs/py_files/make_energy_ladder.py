"""Energy ladders for the accuracy overviews (energy_ladder_<v>.tex).

The ladder is the LADDER view of the website's excimer simulation
(website/simulations/excimer.html) drawn for LIGHT with the model's own
numbers. Every species stands at its formation energy, the energy above
ground-state Kr and F2 at rest (species.py), and the rows that move energy
between the species are arrows labelled with their row ids (network.py).
There is one ladder per version, v0, v0_1 and v0_2. 1Dv0_2 runs v0_2's
chemistry unchanged (its species.py and network.py are v0_2's, which its
test_docs.py checks), so its overview sets v0_2's ladder.

Every number on a ladder comes from the version's code. A level's height and
label are its E_form, and a pair's (dashed) is the sum of its two. An
electron-impact threshold is its cross section's label, a photon is its
row's hv_eV and band, and the heat a photon leaves is its row's
energy_routing. The layout, which column a level sits in, the x of every
arrow and where each label goes, is this script's. An arrow's style says
where its row stands: solid in the default network, dashed in an optional
group, dotted for a band whose default is 0. Every row the version can run
is drawn, or named in NOT_DRAWN, and the script refuses a row it has no
drawing for.

KrF(B) and KrF(C) lie 80 or 120 cm-1 apart, a hundredth of the scale's
smallest gap, so the ladder draws them as a pair of bars, as excimer.html
does, and every B and C row meets the pair.

    python3 LIGHT/docs/make_energy_ladder.py          # writes all three
    python3 LIGHT/docs/make_energy_ladder.py --check  # exit 1 when stale

Each file defines its numbers (\\laddef{key}{value}) and the picture
(\\ladpicture). An overview inputs it in its preamble and sets \\ladpicture
in a figure.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_2.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_2.py").is_file())
if str(LIGHT) not in sys.path:
    sys.path.insert(0, str(LIGHT))

VERSIONS = ("v0", "v0_1", "v0_2", "v0_3", "v0_4", "v0_5")

# rows that leave every level where it is: elastic scattering and F2
# vibration (the electron's loss heats the gas). The caption names them.
NOT_DRAWN = ("G6", "VIB1", "VIB2", "VIB3", "VIB4")

# the vertical scale: cm per eV and the energy range drawn
CM_PER_EV = 1.0
E_LO, E_HI = -3.0, 14.45
HEAD_E = E_HI + 0.55

TEX = {"Kr": "Kr", "Kr*": r"Kr\ensuremath{^*}", "Kr+": r"Kr\ensuremath{^+}",
       "Kr**": r"Kr\ensuremath{^{**}}",
       "Kr2*": r"Kr\ensuremath{_2^*}", "Kr2+": r"Kr\ensuremath{_2^+}",
       "KrF(B)": "KrF(B)", "KrF(C)": "KrF(C)",
       "Kr2F*": r"Kr\ensuremath{_2}F\ensuremath{^*}",
       "F2": r"F\ensuremath{_2}", "F": "F", "F-": r"F\ensuremath{^-}"}
# F2's two dissociative states, as the cross sections name them (sansmath
# has no upright capital Greek, so the Pi is set in the serif math font)
F2_STATE = {"a3": r"{\unsansmath\ensuremath{a\,^3\Pi_u}}",
            "A1": r"{\unsansmath\ensuremath{A\,^1\Pi_u}}"}

# --- the layout (cm) -------------------------------------------------------
X_AXIS = 0.95
KR = (2.65, 6.30)          # Kr atoms, the full width
KR2 = (5.00, 6.30)         # Kr2 dimers, the right of the Kr column
KRF = (7.00, 12.30)        # the KrF column
PAIRS = (10.20, 12.30)     # the ion pairs K26 and K108 start from
BC = (8.40, 11.05)         # the KrF(B) and KrF(C) bars
BC_SEP = 0.03              # half the drawn B-C separation (cm)
T_BAR = (7.30, 8.65)       # Kr2F*, below the pair's left end
F2 = (13.90, 15.10)


def _neg(s: str) -> str:
    """A proper minus sign for LaTeX text."""
    return s.replace("-", r"\ensuremath{-}")


def _e(v: float) -> str:
    """An energy as the overview's formation-energy table prints it."""
    return _neg(f"{v:.6g}")


def _h(v: float) -> str:
    """A derived energy: two decimals, the precision of the inputs, rounded
    half up on its decimal value (9.915 - 6.95 is 2.965, which binary
    floating point stores just below the tie)."""
    d = Decimal(repr(round(v, 9)))
    return str(d.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def _names(rids) -> str:
    return ", ".join(r.replace("_", r"\_") for r in rids)


def _band(label: str) -> tuple:
    """'248 nm (KrF B-X)' -> ('248 nm', 'X'); '420 nm (Kr2F*)' -> (.., None).
    A '~' is the approximate sign, not TeX's tie."""
    m = re.match(r"(~?\d+ nm)(?: \((.*)\))?$", label)
    if m is None:
        # v0_4's 5p -> 5s lines, "5p-5s (near IR)": the parenthetical names
        # the band, and no lower state falls apart
        m2 = re.match(r"[\w-]+ \((.*)\)$", label)
        if m2 is None:
            raise ValueError(f"band label {label!r} not understood")
        return m2.group(1), None
    wl = m.group(1).replace("~", r"\ensuremath{\sim}")
    low = re.search(r"\b[A-Z]-([A-Z])\b", m.group(2) or "")
    return wl, (low.group(1) if low else None)


def _threshold(r) -> str:
    """An electron-impact row's threshold as its cross section prints it."""
    m = re.search(r"\(([\d.]+)eV\)", r.rate[3] or "")
    if m is None:
        raise ValueError(f"{r.rid}: no threshold in {r.rate!r}")
    if abs(float(m.group(1)) - r.e_cost_eV) > 1e-9:
        raise ValueError(f"{r.rid}: label {m.group(1)} is not e_cost "
                         f"{r.e_cost_eV}")
    return m.group(1)


class Version:
    """What one version's code says: energies, rows and their standing."""

    def __init__(self, ver: str):
        self.ver = ver
        sp = importlib.import_module(f"{ver}.species")
        net = importlib.import_module(f"{ver}.network")
        self.version = importlib.import_module(ver).__version__
        self.E = {n: sp.e_form(n) for n in sp.SPECIES}
        self.default = {r.rid: r for r in net.build((), "basic")}
        every = net.build(net.BASIC_GROUPS, "basic")
        self.rows = {r.rid: r for r in every}
        self.bands = getattr(net, "BANDS", {})
        self.route = {r.rid: net.energy_routing(r) for r in every}
        # species only an optional row makes, drawn pale
        made = {p for r in self.default.values() for p in r.products}
        self.optional_only = {n for n in self.E
                              if n not in made and self.E[n] != 0.0}
        self.groups = tuple(net.BASIC_GROUPS)
        from v0_2.species import CM1_EV          # a unit, not a model number
        self.cm1 = CM1_EV

    def style(self, rid: str) -> str:
        r = self.rows[rid]
        if rid not in self.default:
            return "opt"
        if r.rate[0] == "param" and self.bands[r.rate[1]].default == 0.0:
            return "zero"
        # v0_4: a te_fit row whose pre-exponential is a band at 0, and a
        # row whose band fractions (Row.frac) make it 0 by default
        if r.rate[0] == "te_fit" and isinstance(r.rate[1], str) \
                and self.bands[r.rate[1]].default == 0.0:
            return "zero"
        for name, how in getattr(r, "frac", ()):
            d = self.bands[name].default
            if (d if how == "x" else 1.0 - d) == 0.0:
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
            # v0_4: one process whose branches are bands (K26 and K26C, the
            # B and C shares of ion-ion recombination) is drawn once, in
            # the standing of its default branch
            if styles == {"on", "zero"}:
                styles = {"on"}
            else:
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
        down to where it falls apart. The heat's label sits at
        heat_at=(anchor, dx, dE) from the photon's end, or with
        heat_at=("row", k) in row k below the dissociation limit, under the
        photon (neighbouring photons take alternate rows)."""
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
        if heat_at[0] == "row":
            self.label(x, e_low - heat - 0.08 - 0.30 * heat_at[1],
                       rf"{low}{_h(heat)}", anchor="north", style="lad hlab")
        else:
            anchor, dx, de = heat_at
            self.label(x + dx, e_low + de, rf"{low}{_h(heat)}",
                       anchor=anchor, style="lad hlab")
        return e_low

    def stub(self, rids, x, e, below=None, side=None, length=0.8):
        """Rows that return a level's energy to the gas: a short red arrow,
        labelled below its tip (below=text) or beside it (side=(anchor,
        text))."""
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
        if side:
            anchor, text = side
            dx = 0.06 if anchor == "west" else -0.06
            self.label(x + dx, e - 0.5 * length, text, anchor=anchor,
                       style="lad alab")

    def squiggle(self, x, e0, e1):
        self.cmds["arr"].append(
            rf"\draw[lad heat] {_c(x, e0)} -- {_c(x, e1)};")


def build(v: Version) -> Picture:
    """The ladder of one version."""
    E = v.E
    P = Picture(v)
    pale = v.optional_only.__contains__

    # grid, axis, heads
    for k in range(-2, 15, 2):
        P.cmds["grid"].append(
            rf"\draw[lad grid] ({X_AXIS:.3f},{y(k):.3f}) -- "
            rf"({F2[1] + 1.1:.3f},{y(k):.3f});")
        P.label(X_AXIS - 0.08, k, _neg(str(k)), anchor="east",
                style="lad axis")
    P.cmds["grid"].append(
        rf"\draw[lad axisline] ({X_AXIS:.3f},{y(E_LO):.3f}) -- "
        rf"({X_AXIS:.3f},{y(E_HI):.3f});")
    P.label(X_AXIS - 0.08, HEAD_E, "eV", anchor="east", style="lad axis")
    for (x0, x1), head in ((KR, r"Kr and Kr\ensuremath{_2}"),
                           (KRF, r"KrF and Kr\ensuremath{_2}F"),
                           (F2, r"F\ensuremath{_2}")):
        P.label(0.5 * (x0 + x1), HEAD_E, head, anchor="center",
                style="lad head")

    # --- Kr and Kr2 ------------------------------------------------------
    eX, eP, eV, eI = E["Kr*"], E["Kr+"], E["Kr2*"], E["Kr2+"]
    P.level(*KR, 0.0, "ladInk2", kind="ground")
    P.label(KR[0] - 0.08, 0.0, "Kr")
    P.level(*KR, eX, "ladX")
    P.label(KR[0] - 0.08, eX, rf"{TEX['Kr*']} {_e(eX)}")
    P.level(*KR, eP, "ladI")
    P.label(KR[0] - 0.08, eP, rf"{TEX['Kr+']} {_e(eP)}")
    ePP = E.get("Kr**")                 # v0_4: the 5p rung
    if ePP is not None:
        P.level(*KR, ePP, "ladX")
        P.label(KR[0] - 0.08, ePP, rf"{TEX['Kr**']} {_e(ePP)}")
    P.level(*KR2, eV, "ladV", pale=pale("Kr2*"))
    P.label(KR2[0] + 0.03, eV, rf"{TEX['Kr2*']}\\{_e(eV)}",
            anchor="north west", dy=-0.05, style="lad lab2l")
    P.level(*KR2, eI, "ladI", pale=pale("Kr2+"))
    P.label(KR2[1], eI, rf"{TEX['Kr2+']}\\{_e(eI)}", anchor="south east",
            dy=0.05, style="lad lab2")

    # electron impact from the ground state, and Kr* up and down
    P.arrow("C2", (2.80, 0.0), (2.80, eP), "ladE", t=0.52)
    P.arrow("G1", (3.08, 0.0), (3.08, eX), "ladE", t=0.40)
    P.arrow("SUPER", (3.36, eX), (3.36, 0.0), "ladE", t=0.78)
    P.arrow("G2", (3.08, eX), (3.08, eP), "ladX", t=0.36)
    P.arrow("KSTEP", (3.36, eX), (3.36, eP), "ladE", t=0.62)
    c1 = [r for r in ("C1a", "C1b", "C1c") if v.has(r)]
    top = eX
    # the cascade ends on Kr* (v0 to v0_3) or on the 5p rung (v0_4)
    e_end = ePP if ePP is not None else eX
    for i, (r, x) in enumerate(zip(c1, (3.62, 3.80, 3.98))):
        e = v.rows[r].e_cost_eV
        top = max(top, e)
        lab = r"C1a--c" if i == 1 else False
        P.arrow(r, (x, 0.0), (x, e), "ladE", t=0.55, lab=lab)
        if abs(e - e_end) > 1e-9:       # a threshold on the rung's own
            P.level(x - 0.07, 4.40, e, "ladMuted", kind="tick")   # bar needs
            P.label(4.46, e, _threshold(v.rows[r]), anchor="west",  # no tick
                    style="lad tlab")
    P.squiggle(4.40, top, e_end)
    if ePP is not None:
        # the rung's rows: the electrons pump Kr* up (K76) and return it
        # (K76r), Kr relaxes it (KPS, a quench), it radiates to Kr* (PRAD),
        # the electrons ionize both (KSTEP_S, KSTEP_P), and the harpoon
        # from 5p (K27, K27C) with its dissociating branch (K27X)
        P.photon("PRAD", 4.24, ePP, "ladX", heat_at=("west", 0.10, -0.30),
                 t=0.5)
        P.arrow("K76", (4.58, eX), (4.58, ePP), "ladE", t=0.32)
        P.arrow("K76r", (4.82, ePP), (4.82, eX), "ladE", t=0.30)
        P.stub(("KPS",), 5.06, ePP, length=ePP - eX,
               side=("west", "KPS"))
        # the stepwise ionizations beside C2 and the optional KSTEP, their
        # labels high up where nothing else sits
        P.arrow("KSTEP_S", (2.94, eX), (2.94, eP), "ladE", t=0.85)
        P.arrow("KSTEP_P", (3.50, ePP), (3.50, eP), "ladE", t=0.85)
        if v.has("K27X"):
            P.drawn.add("K27X")
            P.cmds["arr"].append(
                rf"\draw[lad arr,{P._style(('K27X',))},draw=ladH] "
                rf"{_c(KR[0] + 0.05, ePP)} -- {_c(KR[0] - 0.30, ePP - 0.85)};")
            P.label(KR[0] - 0.34, ePP - 0.72, "K27X", anchor="east",
                    style="lad alab")
    # G3X: Kr* + F2 -> Kr + 2F, a quench, leaves from Kr*'s left end
    if v.has("G3X"):
        P.drawn.add("G3X")
        P.cmds["arr"].append(
            rf"\draw[lad arr,{P._style(('G3X',))},draw=ladH] "
            rf"{_c(KR[0] + 0.05, eX)} -- {_c(KR[0] - 0.30, eX - 0.85)};")
        P.label(KR[0] - 0.34, eX - 0.72, "G3X", anchor="east",
                style="lad alab")
    # the dimers
    P.arrow("AK8_Kr", (5.10, eX), (5.10, eV), "ladX", t=0.5,
            side=("east", -0.06, 0.0))
    P.arrow("K15", (5.25, eP), (5.25, eI), "ladI", t=0.5)
    # Kr2+ + e makes Kr* (v0 to v0_3) or the 5p rung (v0_4)
    P.arrow("R11p", (5.80, eI), (5.80, e_end), "ladI", t=0.45)
    P.arrow("PEN_KR2", (6.12, eV), (6.12, eI), "ladV", t=0.62)
    P.photon("KR2_RAD", 5.88, eV, "ladV", heat_at=("west", 0.10, -0.40))

    # --- KrF and Kr2F ----------------------------------------------------
    eB, eC, eT = E["KrF(B)"], E["KrF(C)"], E["Kr2F*"]
    e_ip, e_ip2 = v.pair("Kr+", "F-"), v.pair("Kr2+", "F-")
    e_kf = v.pair("Kr", "F")
    P.level(*PAIRS, e_ip, "ladInk2", kind="pair", pale=pale("F-"))
    P.label(PAIRS[1] + 0.08, e_ip,
            rf"{TEX['Kr+']} + {TEX['F-']} {_e(e_ip)}", anchor="west")
    P.level(*PAIRS, e_ip2, "ladInk2", kind="pair",
            pale=pale("Kr2+") or pale("F-"))
    P.label(PAIRS[1] + 0.08, e_ip2,
            rf"{TEX['Kr2+']} + {TEX['F-']} {_e(e_ip2)}", anchor="west")
    P.level(KRF[0], PAIRS[1], e_kf, "ladInk2", kind="pair")
    P.label(PAIRS[1], e_kf, rf"Kr + F {_e(e_kf)}", anchor="south east",
            dy=0.04)
    # B above C, drawn BC_SEP either side of their mean
    P.level(*BC, eB, "ladB", dy=BC_SEP)
    P.level(*BC, eC, "ladC", dy=-BC_SEP)
    gap = (eB - eC) / v.cm1
    bc_rows = _names([r for r in ("BC", "CB") if v.has(r)])
    P.drawn |= {r for r in ("BC", "CB") if v.has(r)}
    P.label(BC[1] + 0.08, eB, rf"{TEX['KrF(B)']} {_e(eB)}\\"
            rf"{TEX['KrF(C)']} {gap:.0f}\,cm\ensuremath{{^{{-1}}}} "
            rf"below\\{bc_rows} mix them", anchor="west", style="lad lab2l")
    P.level(*T_BAR, eT, "ladT")
    P.label(T_BAR[0] - 0.08, eT, rf"{TEX['Kr2F*']}\\{_e(eT)}",
            style="lad lab2")

    # formation: the lower source (Kr2*) lands left of the upper (Kr*),
    # and the ion pairs drop to the right of both, so nothing crosses
    P.arrow("K2F_BR", (KR2[1], eV), (7.62, eT), "ladV", t=0.42)
    P.arrow("NK7", (KR2[1], eV), (8.66, eB), "ladV", t=0.42)
    P.arrow("K2F_F", (KR2[1], eV), (8.90, eB), "ladV", t=0.62)
    P.arrow(("G3", "G3C"), (KR[1], eX), (9.14, eB), "ladX", t=0.40)
    if ePP is not None:
        P.arrow(("K27", "K27C"), (KR[1], ePP), (9.26, eB), "ladX", t=0.30)
    # ion-ion recombination: the B and C branches as one arrow (v0_4 adds
    # K26C and K108C, bands at 0), and the dissociating branches as stubs
    P.arrow(("K108", "K108C"), (10.45, e_ip2), (10.45, eB), "ladI", t=0.5)
    P.arrow(("K26", "K26C"), (10.72, e_ip), (10.72, eB), "ladI", t=0.62)
    if v.has("K26X"):
        P.stub(("K108X",), 10.30, e_ip2, side=("east", "K108X"))
        P.stub(("K26X",), 10.30, e_ip, side=("east", "K26X"), length=0.6)
    # the trimer, and light and loss (each photon leaves between two
    # arrivals, so no arrow seems to run through the bars)
    P.arrow(("G5", "G5C"), (8.52, eB), (8.52, eT), "ladB", t=0.5,
            side=("east", -0.05, 0.0), lab=r"G5,\\G5C", lstyle="lad slab")
    P.photon("G4", 9.35, eB, "ladB", heat_at=("row", 1))
    P.photon("G4A", 9.75, eB, "ladB", heat_at=("row", 0))
    P.photon("G4C", 10.15, eC, "ladC", heat_at=("row", 1))
    b_q = [r for r in ("QF2", "QKR", "QE", "QF") if v.has(r)]
    c_q = [r for r in ("QF2C", "QKRC", "QEC", "QFC") if v.has(r)]
    P.stub(tuple(b_q + c_q), 10.93, eC, length=0.8)
    P.label(10.62, eC - 0.85, rf"{_names(b_q)}\\{_names(c_q)}",
            anchor="north west", style="lad alab")
    P.photon("C3", 7.98, eT, "ladT", heat_at=("row", 0), t=0.64)
    P.stub(("Q2F",), 7.52, eT, below="Q2F")
    P.stub(("Q2FF",), 8.40, eT, below="Q2FF")

    # --- F2 ----------------------------------------------------------------
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
        raise ValueError(f"{v.ver}: rows with no drawing: {sorted(missing)}")
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

    n = {"version": v.version,
         "E.B": _e(E["KrF(B)"]),
         "DeB": _h(v.pair("Kr+", "F-") - E["KrF(B)"]),
         "BoverKrF": _h(E["KrF(B)"] - v.pair("Kr", "F")),
         "wl.G4": wl("G4"), "hv.G4": _h(v.rows["G4"].hv_eV),
         "heat.G4": _h(R["G4"]["heat"]),
         "heat.G3": _h(R["G3"]["heat"]),
         "heat.K26": _h(R["K26"]["heat"]),
         "heat.G5": _h(R["G5"]["heat"]),
         "wl.C3": wl("C3"), "hv.C3": _h(v.rows["C3"].hv_eV),
         "heat.C3": _h(R["C3"]["heat"]),
         "gapBC": f"{(E['KrF(B)'] - E['KrF(C)']) / v.cm1:.0f}",
         "notdrawn": ids(r for r in NOT_DRAWN if v.has(r)),
         "toe": ids(r for r in v.rows if R[r]["to_e"] > 0),
         "zero": ids(r for r in v.rows if v.style(r) == "zero"),
         "groups": _and(rf"\texttt{{{g}}}" for g in v.groups),
         "pale": _and(TEX[s] for s in v.E if s in v.optional_only)}
    if v.has("KR2_RAD"):
        n["wl.V"] = wl("KR2_RAD")
        n["hv.V"] = _h(v.rows["KR2_RAD"].hv_eV)
        n["heat.V"] = _h(R["KR2_RAD"]["heat"])
    if "Kr**" in E:                      # v0_4: the 5p rung
        n["E.PP"] = _e(E["Kr**"])
        n["gap5p"] = _h(E["Kr**"] - E["Kr*"])
        n["hv.PRAD"] = _h(v.rows["PRAD"].hv_eV)
        n["heat.KPS"] = _h(R["KPS"]["heat"])
        n["heat.C1b"] = _h(R["C1b"]["heat"])
        n["heat.R11p"] = _h(R["R11p"]["heat"])
    return n


PREAMBLE = (r"\providecommand\laddef[2]{\expandafter\def"
            r"\csname lad@#1\endcsname{#2}}")

STYLES = r"""\definecolor{ladE}{HTML}{2A78D6}% electron gas
\definecolor{ladX}{HTML}{EDA100}% Kr*
\definecolor{ladB}{HTML}{1BAF7A}% KrF(B), 248 and 272 nm
\definecolor{ladC}{HTML}{008300}% KrF(C), 275 nm
\definecolor{ladT}{HTML}{4A3AA7}% Kr2F*, 420 nm
\definecolor{ladV}{HTML}{A63A83}% Kr2*, 147 nm
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
    fill=white, fill opacity=0.85, text opacity=1, align=right},
  lad hlab/.style={font=\tiny, text=ladInk2, inner sep=0.8pt,
    fill=white, fill opacity=0.85, text opacity=1},
  lad tlab/.style={font=\tiny, text=ladInk2, inner sep=0.8pt},
  lad axis/.style={font=\scriptsize, text=ladInk2, inner sep=1pt},
  lad head/.style={font=\footnotesize\bfseries, text=ladInk2},
}"""


def legend(v: Version) -> list:
    """The key, two columns below the ladder."""
    items = [("lad arr,lad on,draw=ladInk2", "arrow",
              "a row of the default network")]
    if any(v.style(r) == "opt" for r in v.rows):
        items.append(("lad arr,lad opt,draw=ladInk2", "arrow",
                      "a row of an optional group, off by default"))
    if any(v.style(r) == "zero" for r in v.rows):
        items.append(("lad arr,lad zero,draw=ladInk2", "arrow",
                      "a band whose default is 0, no measurement"))
    items += [("lad hv,lad on,draw=ladB", "arrow",
               "a photon, with its row and band"),
              ("lad arr,lad on,draw=ladH", "arrow",
               "a quench, the level's energy to heat"),
              ("lad heat", "arrow",
               "heat as a lower state falls apart (eV) or cascades"),
              ("lad pair,draw=ladInk2", "bar",
               "a pair, the sum of two formation energies"),
              ("lad tick,draw=ladMuted", "bar",
               "a threshold, or where a photon leaves its lower state")]
    if v.optional_only:
        items.append(("lad lev,draw=ladInk2!40", "bar",
                      "a species that only an optional group makes"))
    out = []
    x = X_AXIS + 0.3
    col_w = 7.6
    y0 = y(E_LO) - 0.55
    for i, (st, _kind, text) in enumerate(items):
        xi = x + (i % 2) * col_w
        yi = y0 - (i // 2) * 0.42
        out.append(rf"\draw[{st}] ({xi:.3f},{yi:.3f}) -- "
                   rf"({xi + 0.7:.3f},{yi:.3f});")
        out.append(rf"\node[lad lab,anchor=west] at ({xi + 0.8:.3f},"
                   rf"{yi:.3f}) {{{text}}};")
    return out


def render(ver: str) -> str:
    v = Version(ver)
    P = build(v)
    head = (f"% LIGHT {ver} energy ladder. GENERATED by make_energy_ladder.py\n"
            f"% from LIGHT {ver} ({v.version}). Do not edit. Rerun the "
            f"script instead.\n")
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
                    help="exit 1 if a ladder file is stale")
    ap.add_argument("versions", nargs="*", default=list(VERSIONS))
    a = ap.parse_args(argv)
    stale = []
    for ver in a.versions:
        out = HERE / f"energy_ladder_{ver}.tex"
        text = render(ver)
        if a.check:
            if not out.is_file() or out.read_text() != text:
                stale.append(out.name)
        else:
            out.write_text(text)
            print(f"wrote {out.name}")
    if stale:
        print("stale: " + ", ".join(stale) + ". Rerun make_energy_ladder.py")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
