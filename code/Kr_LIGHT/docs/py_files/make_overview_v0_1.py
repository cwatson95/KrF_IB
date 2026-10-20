"""Numbers for 'How accurate LIGHT_Kr v0_1 is' (accuracy_overview_v0_1.tex).

Writes overview_v0_1_numbers.tex, one \\wgdef{key}{value} per number and a
few \\def tables, from the model:

  the swarm      v0_1's two-term solver (boltzmann.Solver, the kinetic
                 EEDF's own) on its Biagi v7.1 cross sections, without
                 electron-electron collisions, against krypton's measured
                 drift velocity (Pack, Voshall & Phelps 1962) and
                 ionization coefficient (Kruithof 1940, Heylen 1971), as
                 Dutton's 1975 survey prints them. The data are typed below
                 from Dutton's tables (DUTTON_W, DUTTON_ALPHA), printed
                 values, no figure readings.
  the e-beam     v0_1 under the beam drive (kr_docs.BeamModel) against
                 Morozov et al. 2008 (0.42 +- 0.05 above 500 hPa) and
                 Eckstrom et al. 1988 (near 46 % up to 3e5 W/cm3), and the
                 drive's partition against Saito et al. 2003 (W_s 15.4 eV,
                 N_ex/N_i 0.55).
  uncertainty    every band the rows read at its two ends and the three
                 measured constants with a quoted error at theirs, one at a
                 time, on both closures (kr_docs.band_params,
                 const_params).
  the reference  the kinetic EEDF and the two-temperature closure, the
                 5s -> 5p scale, Kannari's fit, the scan, and the kinetic
                 closure's convergence (kr_docs.conv_params).
  the conversion v0_5 (no F2, its own tables and v0_1's) against v0_1.

    python3 LIGHT_Kr/docs/make_overview_v0_1.py        # minutes (cached)
    cd LIGHT_Kr/docs && latexmk -pdf accuracy_overview_v0_1.tex

The runs come from kr_docs.run_many (cached per version of the model's
code). --check compares the file with a fresh computation and writes
nothing.
"""

from __future__ import annotations

import copy
import importlib
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kr_docs as K                                        # noqa: E402
from v0_1 import boltzmann as KB                           # noqa: E402
from v0_1 import jung2006, xsec                            # noqa: E402
from v0_1.model import Model, Params, run                  # noqa: E402
from v0_1.network import BANDS, KR_5P_TAU_S, build         # noqa: E402
from v0_1.report import row_dict                           # noqa: E402
from v0_1.species import e_form                            # noqa: E402

OUT = HERE / "overview_v0_1_numbers.tex"

# ------------------------------------------------------------- the data
# Pure Kr, the electron drift velocity W [cm/s] against E/N [Td], as
# Dutton (J. Phys. Chem. Ref. Data 4, 577 (1975)) prints Pack, Voshall &
# Phelps (Phys. Rev. 127, 2084 (1962)): Table 1.9(a), 300 K (key a) and
# 368 K (key b), and Table 1.9(b), 195 K. Printed values (the table's unit
# 1e5 cm/s multiplied out). The table prints 1.04e5 at 0.02772 and 0.0306
# Td and 2.51e5 at 1.786 and 2.440 Td; copied as printed.
DUTTON_W = (
    (0.000791, 1.08e3, 300), (0.000924, 1.28e3, 300), (0.001525, 2.26e3, 368),
    (0.002144, 3.1e3, 368), (0.00306, 4.2e3, 368), (0.00308, 3.6e3, 300),
    (0.00458, 5.8e3, 300), (0.00461, 6.6e3, 368), (0.00617, 8.5e3, 300),
    (0.00622, 9.5e3, 368), (0.00774, 1.06e4, 300), (0.00918, 1.33e4, 300),
    (0.01254, 2.10e4, 300), (0.01550, 3.3e4, 300), (0.01712, 4.2e4, 300),
    (0.01838, 5.1e4, 300), (0.02118, 6.2e4, 300), (0.02440, 8.6e4, 300),
    (0.02772, 1.04e5, 300), (0.0306, 1.04e5, 300), (0.0758, 1.27e5, 300),
    (0.1496, 1.34e5, 300), (0.304, 1.57e5, 300), (0.609, 1.78e5, 300),
    (0.892, 2.00e5, 300), (1.53, 2.30e5, 300), (1.786, 2.51e5, 300),
    (2.440, 2.51e5, 300), (3.06, 3.2e5, 300),
    (0.001210, 1.58e3, 195), (0.001561, 1.74e3, 195), (0.001851, 2.38e3, 195),
    (0.00308, 3.8e3, 195), (0.00465, 5.2e3, 195), (0.00617, 6.9e3, 195),
    (0.00918, 1.09e4, 195), (0.01236, 1.62e4, 195), (0.01891, 4.0e4, 195),
    (0.02475, 8.0e4, 195), (0.0442, 1.19e5, 195), (0.0924, 1.27e5, 195),
)
# The first Townsend coefficient alpha/N [cm2] against E/N [Td], Dutton's
# Table 7.4(a): Kruithof (Physica 7, 519 (1940), 'K') and Heylen (Int. J.
# Electron. 31, 19 (1971), 'H'), printed values (unit 1e-17 cm2 multiplied
# out); E/N with N at 0 C as Dutton converts E/p0. To 1017 Td: above it the
# mean energy passes 13 eV and the two-term equation is not the tool.
DUTTON_ALPHA = (
    (14.12, 4.4e-21, "K"), (15.54, 6.8e-21, "K"), (16.95, 1.03e-20, "K"),
    (19.77, 2.18e-20, "K"), (22.60, 4.1e-20, "K"), (25.42, 7.1e-20, "K"),
    (27.51, 4.6e-20, "H"), (28.25, 1.17e-19, "K"), (30.5, 1.18e-19, "H"),
    (30.8, 1.14e-19, "H"), (31.1, 1.82e-19, "K"), (33.7, 1.73e-19, "H"),
    (33.9, 2.68e-19, "K"), (33.9, 1.97e-19, "H"), (36.4, 2.61e-19, "H"),
    (36.9, 3.0e-19, "H"), (39.6, 5.2e-19, "K"), (42.3, 4.6e-19, "H"),
    (43.1, 6.4e-19, "H"), (45.2, 8.9e-19, "K"), (49.4, 8.7e-19, "H"),
    (50.9, 1.37e-18, "K"), (54.9, 1.35e-18, "H"), (56.5, 1.95e-18, "K"),
    (61.5, 1.98e-18, "H"), (62.2, 2.62e-18, "K"), (62.4, 2.33e-18, "H"),
    (68.4, 2.82e-18, "H"), (68.4, 3.16e-18, "H"), (70.6, 3.8e-18, "K"),
    (77.2, 3.8e-18, "H"), (78.3, 4.5e-18, "H"), (79.1, 5.1e-18, "K"),
    (87.1, 5.6e-18, "H"), (87.1, 5.9e-18, "H"), (90.4, 7.1e-18, "K"),
    (99.0, 7.9e-18, "H"), (101.7, 9.3e-18, "K"), (110.9, 1.08e-17, "H"),
    (113.0, 1.18e-17, "K"), (124.2, 1.42e-17, "H"), (127.1, 1.51e-17, "K"),
    (130.6, 1.48e-17, "H"), (139.2, 1.78e-17, "H"), (140.2, 1.71e-17, "H"),
    (141.2, 1.87e-17, "K"), (157.1, 2.22e-17, "H"), (166.4, 2.39e-17, "H"),
    (169.5, 2.62e-17, "K"), (193.2, 3.2e-17, "H"), (197.7, 3.4e-17, "K"),
    (203.1, 3.4e-17, "H"), (226.0, 4.2e-17, "K"), (246.1, 4.8e-17, "H"),
    (249.7, 4.7e-17, "H"), (254.2, 4.9e-17, "K"), (282.5, 5.7e-17, "K"),
    (339, 7.4e-17, "K"), (396, 9.0e-17, "K"), (452, 1.06e-16, "K"),
    (509, 1.22e-16, "K"), (565, 1.37e-16, "K"), (678, 1.64e-16, "K"),
    (791, 1.89e-16, "K"), (904, 2.12e-16, "K"), (1017, 2.33e-16, "K"),
)
ALPHA_RANGES = ((14.0, 50.0), (50.0, 100.0), (100.0, 300.0), (300.0, 1100.0))
W_RANGES = ((0.0, 0.01), (0.01, 0.05), (0.05, 3.1))
# the e-beam measurements
MOROZOV = (0.42, 0.05)        # J. Appl. Phys. 103, 103301 (2008), > 500 hPa
ECKSTROM = 0.46               # J. Appl. Phys. 64, 1691 (1988), abstract:
#                               'near the theoretical limit of ~46%'
SAITO_WS, SAITO_NEX = 15.4, 0.55   # IEEE TNS 50, 2452 (2003), abstract
# the stated uncertainties of Jung et al.'s cross sections (their Sec. II F
# and Tables II, III): out of 1s5 +-35 %, out of 1s3 +-40 %
JUNG_UNC = (35, 40)
CELLS_SWARM = 1000


def _f(x, d=1):
    s = f"{x:.{d}f}"
    return s.replace("-0." + "0" * d, "0." + "0" * d) if d else s


def _pm(x, d=0):
    """A signed percent as text: +31, -4."""
    s = f"{x:+.{d}f}"
    return "0" if s in ("+0", "-0") else s.replace("-", "$-$")


def _sci(x, d=1):
    if x == 0.0:
        return "0"
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


# ------------------------------------------------------------- the swarm
_PROCS = None


def _procs():
    global _PROCS
    if _PROCS is None:
        _PROCS = [KB.Process("Kr", p.kind, p.threshold, p)
                  for p in xsec.processes("Kr")
                  if p.kind in ("excitation", "ionization")]
    return _PROCS


def swarm_point(en_td: float, tg: float, cells: int = CELLS_SWARM) -> dict:
    """A swarm in pure Kr at E/N: the two-term solution with the temporal
    growth term, no e-e (n_e/N = 4e-29), on a grid whose EEDF falls 11 to
    25 decades (the solver's own window; bracketed in emax). Returns the
    drift velocity [cm/s] (Joule power over eE, the solver's own flux),
    alpha/N = k_ion / W [cm2] and the mean energy [eV]."""
    n_col = 2.5e25
    p_ion = next(q.p for q in _procs() if q.kind == "ionization")
    emax = 1.5 if en_td < 0.01 else 4.0 * max(en_td, 0.1) ** 0.35 + 8.0
    lo = hi = None
    for _ in range(30):
        s = KB.Solver(_procs(), n_cells=cells, emax=emax, elastic=("Kr",),
                      tol=1e-10, maxit=4000)
        st = KB.State(x={"Kr": 1.0}, tg=tg, n_col=n_col, ne=1e-3)
        sol = s.solve(st, en_td)
        dec = KB._decades(sol.F)
        if 11.0 <= dec <= 25.0:
            break
        if dec < 11.0:
            lo = emax
            emax = emax * 1.6 if hi is None else math.sqrt(emax * hi)
        else:
            hi = emax
            emax = emax / 1.6 if lo is None else math.sqrt(lo * emax)
    else:
        raise RuntimeError(f"no grid in the window at {en_td} Td")
    w = s.field_power(st, sol) / (en_td * KB.TD * n_col)      # m/s
    k = sol.grid.rate(p_ion, sol.F)                            # m3/s
    return {"W": w * 100.0, "alphaN": k / w * 1e4, "mean": sol.mean,
            "emax": emax, "decades": dec}


def _rms(e):
    e = np.asarray(e)
    return float(np.sqrt(np.mean(e ** 2)))


def swarm_numbers() -> tuple:
    n = {}
    ew, ea = [], []
    for en, wm, tg in DUTTON_W:
        d = swarm_point(en, tg)
        ew.append((en, tg, d["W"] / wm - 1.0, wm, d["W"]))
    for en, am, src in DUTTON_ALPHA:
        d = swarm_point(en, 300.0)
        ea.append((en, src, d["alphaN"] / am - 1.0, am, d["alphaN"]))
    e = [x[2] for x in ew]
    n["sw.W.n"] = str(len(ew))
    n["sw.W.rms"] = _f(100 * _rms(e))
    n["sw.W.mean"] = _pm(100 * np.mean(e), 1)
    n["sw.W.enlo"] = f"{min(x[0] for x in ew):.4f}"
    n["sw.W.enhi"] = f"{max(x[0] for x in ew):.2f}"
    worst = max(ew, key=lambda x: abs(x[2]))
    n["sw.W.worst"] = _pm(100 * worst[2])
    n["sw.W.worst.en"] = f"{worst[0]:.3g}"
    for i, (a, b) in enumerate(W_RANGES):
        sel = [x[2] for x in ew if a <= x[0] < b]
        n[f"sw.W.r{i}.n"] = str(len(sel))
        n[f"sw.W.r{i}.rms"] = _f(100 * _rms(sel))
        n[f"sw.W.r{i}.mean"] = _pm(100 * np.mean(sel), 1)
    for i, (a, b) in enumerate(ALPHA_RANGES):
        sel = [x[2] for x in ea if a <= x[0] < b]
        n[f"sw.a.r{i}.n"] = str(len(sel))
        n[f"sw.a.r{i}.rms"] = _f(100 * _rms(sel))
        n[f"sw.a.r{i}.mean"] = _pm(100 * np.mean(sel))
        for src in ("K", "H"):
            ss = [x[2] for x in ea if a <= x[0] < b and x[1] == src]
            n[f"sw.a.r{i}.{src}"] = _pm(100 * np.mean(ss)) if ss else "--"
    run_pts = sorted((en, wm) for en, wm, tg in DUTTON_W
                     if tg == 300 and W_RANGES[1][0] <= en < W_RANGES[1][1])
    n["sw.W.r1.wx"] = f"{run_pts[-1][1] / run_pts[0][1]:.0f}"
    n["sw.W.r1.ex"] = f"{run_pts[-1][0] / run_pts[0][0]:.1f}"
    first300 = min(en for en, _, _ in DUTTON_ALPHA if en >= 300.0)
    n["sw.a.r3.eps"] = f"{swarm_point(first300, 300.0)['mean']:.0f}"
    raw = {"r0": np.mean([x[2] for x in ea if x[0] < 50.0]),
           "r3": np.mean([x[2] for x in ea if x[0] >= 300.0])}
    if not (raw["r0"] < 0.0 < raw["r3"]):
        raise AssertionError("the text says the model falls short below 50 "
                             "Td and lies above past 300 Td")
    n["sw.a.n"] = str(len(ea))
    lo14 = min(ea, key=lambda x: x[0])
    n["sw.a.lo.en"] = f"{lo14[0]:.1f}"
    n["sw.a.lo.x"] = f"{lo14[4] / lo14[3]:.2f}"
    # the grid: the same scores on half and twice the cells
    alt = {}
    for cells in (CELLS_SWARM // 2, CELLS_SWARM * 2):
        w_e = [swarm_point(en, tg, cells)["W"] / wm - 1.0
               for en, wm, tg in DUTTON_W]
        a_e = [swarm_point(en, 300.0, cells)["alphaN"] / am - 1.0
               for en, am, _ in DUTTON_ALPHA]
        alt[cells] = (100 * _rms(w_e), [100 * _rms(
            [x for (en, _, _), x in zip(DUTTON_ALPHA, a_e) if a <= en < b])
            for a, b in ALPHA_RANGES])
    base = (100 * _rms(e), [100 * _rms([x[2] for x in ea if a <= x[0] < b])
                            for a, b in ALPHA_RANGES])
    dmax = max(max(abs(alt[c][0] - base[0]),
                   max(abs(u - v) for u, v in zip(alt[c][1], base[1])))
               for c in alt)
    n["sw.cells"] = str(CELLS_SWARM)
    n["sw.grid"] = f"{dmax:.2f}"
    # the figure: the model on a fine grid and the measured points
    lines = []
    ens = np.geomspace(5e-4, 4.0, 60)
    for tg, key, top in ((300.0, "Wa", 4.0), (368.0, "Wb", 0.012),
                         (195.0, "Wc", 0.2)):
        pts = " ".join(f"({en:.5g},{swarm_point(en, tg)['W']:.5g})"
                       for en in ens if en <= top)
        lines.append(f"\\def\\wgSw{key}model{{{pts}}}")
    for tg, key in ((300, "Wa"), (368, "Wb"), (195, "Wc")):
        pts = " ".join(f"({en:.5g},{wm:.5g})" for en, wm, t in DUTTON_W
                       if t == tg)
        lines.append(f"\\def\\wgSw{key}data{{{pts}}}")
    ens = np.geomspace(12.0, 1100.0, 50)
    pts = " ".join(f"({en:.5g},{swarm_point(en, 300.0)['alphaN']:.5g})"
                   for en in ens)
    lines.append(f"\\def\\wgSwAmodel{{{pts}}}")
    for src in ("K", "H"):
        pts = " ".join(f"({en:.5g},{am:.5g})" for en, am, s in DUTTON_ALPHA
                       if s == src)
        lines.append(f"\\def\\wgSwA{src}data{{{pts}}}")
    return n, lines


# ------------------------------------------------------------- e-beam
def beam_numbers() -> dict:
    n = {}
    Y = K.beam_yields()
    hv = next(r.hv_eV for r in build() if r.rid == "KR2_RAD")
    ceil = (1.0 + Y["nex"]) * hv / Y["W"]
    n.update({
        "beam.W": f"{Y['W']:.1f}", "beam.nex": f"{Y['nex']:.2f}",
        "beam.pp": f"{Y['pp']:.3f}", "beam.s": f"{Y['s']:.3f}",
        "beam.ws": f"{Y['ws']:.2f}", "beam.e0": f"{Y['e0']:.0f}",
        "beam.ceiling.eV": f"{Y['ceiling']:.0f}",
        "beam.eps": f"{Y['eps_sec']:.2f}",
        "beam.cascade": f"{Y['cascade']:.2f}",
        "beam.ceil": f"{ceil:.3f}",
        "beam.php": f"{1.0 + Y['nex']:.3f}",
        "hv147": f"{hv:.2f}",
        "saito.ws": f"{SAITO_WS:.1f}", "saito.nex": f"{SAITO_NEX:.2f}",
        "saito.ideal": f"{hv / SAITO_WS:.3f}",
        "saito.W": f"{SAITO_WS * (1 + SAITO_NEX):.1f}",
        "moz.short": f"{100 * (1 - MOROZOV[0] / (hv / SAITO_WS)):.0f}",
        "beam.ws.dev": _pm(100 * (Y["ws"] / SAITO_WS - 1), 1),
        "beam.nex.dev": _pm(100 * (Y["nex"] / SAITO_NEX - 1)),
        "beam.W.dev": _pm(100 * (Y["W"] / (SAITO_WS * (1 + SAITO_NEX)) - 1)),
        "moz.meas": f"{MOROZOV[0]:.2f}", "moz.err": f"{MOROZOV[1]:.2f}",
        "eck.meas": f"{ECKSTROM:.2f}",
        "moz.atm": f"{K.MOROZOV_ATM:.3f}", "moz.hpa": "1000",
        "moz.w.lo": f"{K.MOROZOV_W[0]:.0f}", "moz.w.hi": f"{K.MOROZOV_W[1]:.0f}",
        "eck.w": f"{K.ECKSTROM_W / 1e3:.0f}",
        "beam.pulse": f"{K.BEAM_PULSE * 1e6:g}",
        "beam.tail.lo": f"{K.BEAM_TAIL_LO * 1e6:g}",
        "beam.tail.hi": f"{K.BEAM_TAIL_HI * 1e6:g}",
    })
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        bp = K.beam_params(eedf)
        keys = list(bp)
        got = dict(zip(keys, K.run_many([("beam", bp[k]) for k in keys])))
        for k, s in got.items():
            pairs = s["e_dep_eV"] / s["beam"]["W"]
            ev = s["events"]
            vuv = s["emission_vuv"] / 100.0
            name = f"{k[0]}.{k[1]:g}".replace("eckstrom", "eck").replace(
                "morozov", "moz")
            n[f"{name}.{tag}"] = f"{vuv:.3f}"
            n[f"{name}.{tag}.G2"] = f"{ev['G2'] / pairs:.3f}"
            n[f"{name}.{tag}.G1"] = f"{(ev['G1'] + ev['C1a'] + ev['C1b'] + ev['C1c']) / pairs:.3f}"
            n[f"{name}.{tag}.te"] = f"{s['te_mean_eV']:.2f}"
            n[f"{name}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
            n[f"{name}.{tag}.ph"] = f"{ev['KR2_RAD'] / pairs:.3f}"
            n[f"{name}.{tag}.rec"] = f"{ev['R11p'] / pairs:.3f}"
            n[f"{name}.{tag}.krs"] = _sci(s["n_max"]["Kr*"])
            n[f"{name}.{tag}.ne"] = _sci(s["ne_peak"])
            n[f"{name}.{tag}.closure"] = _sci(abs(s["closure"]), 0)
            # the loss below v0_1's ceiling (one photon per excitation and
            # per ion), the quantity a report 'near the limit' tests
            n[f"{name}.{tag}.loss"] = f"{100 * (1 - vuv / ceil):.0f}"
        m = got[("morozov", K.MOROZOV_W[1])]
        v = m["emission_vuv"] / 100.0
        n[f"moz.{tag}.dev"] = _pm(100 * (v / MOROZOV[0] - 1))
        n[f"moz.{tag}.sig"] = f"{(v - MOROZOV[0]) / MOROZOV[1]:.1f}"
        e1 = got[("eckstrom", 1.0)]["emission_vuv"] / 100.0
        n[f"eck.{tag}.dev"] = _pm(100 * (e1 / ECKSTROM - 1))
        # what the text says of these runs
        if not v > MOROZOV[0]:
            raise AssertionError("the text says v0_1 lies above Morozov")
        if not 1 - got[("eckstrom", 2.0)]["emission_vuv"] / 100 / ceil < 0.10:
            raise AssertionError("the text says v0_1 loses little at 2 atm")
        if not (got[("eckstrom", 0.5)]["emission_vuv"] < 100 * e1
                < got[("eckstrom", 2.0)]["emission_vuv"]
                and e1 < ECKSTROM
                and got[("eckstrom_pen", 1.0)]["emission_vuv"] > 100 * e1
                and got[("eckstrom_form", 1.0)]["emission_vuv"] > 100 * e1):
            raise AssertionError("the text's account of the high-pump runs "
                                 "no longer holds")
        n[f"moz.{tag}.flat"] = _pm(100 * (
            got[("morozov", K.MOROZOV_W[0])]["emission_vuv"]
            / got[("morozov", K.MOROZOV_W[1])]["emission_vuv"] - 1), 1)
    return n


# ------------------------------------------------------ the uncertainty
UNC_LABELS = {
    "k_5s5p": r"\el{} + \Krs{} \ra \Krpp{}, scale",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{}, Penning",
    "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr",
    "kr2_pen": r"\Krts{} + \Krts{}, Penning",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}",
    "R11p": r"\Krtp{} + \el{} \ra \Krpp{} + Kr",
    "K15": r"\Krp{} + 2Kr \ra \Krtp{} + Kr",
    "KR2_RAD": r"\Krts{} lifetime",
}
UNC_SHORT = {
    "k_5s5p": r"the 5s \ra 5p cross section", "kps_kr": r"\Krpp{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{}", "kr2_form": r"\Krs{} + 2Kr",
    "kr2_pen": r"\Krts{} + \Krts{}", "kstep_s": r"\el{} + \Krs{} \ra \Krp",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp", "R11p": r"\Krtp{} + \el",
    "K15": r"\Krp{} + 2Kr", "KR2_RAD": r"the \Krts{} lifetime",
}


def _ends_text(name, lo, hi):
    if name in BANDS:
        b = BANDS[name]
        if b.unit == "factor":
            return f"{lo:g} to {hi:g}"
        unit = {"cm3/s": r"\unit{cm^3\,s^{-1}}",
                "cm6/s": r"\unit{cm^6\,s^{-1}}"}[b.unit]
        lo_t = "0" if lo == 0 else _sci(lo)
        hi_t = "0" if hi == 0 else _sci(hi)
        return f"${lo_t}$ to ${hi_t}${unit}"
    if name == "R11p":
        return f"$\\pm{100 * K.R11P_ERR:.0f}\\,\\%$"
    if name == "K15":
        return f"${_sci(K.K15_LO)}$ to ${_sci(K.K15_HI)}\\unit{{cm^6\\,s^{{-1}}}}$"
    if name == "KR2_RAD":
        return (f"{K.TAU_RAD - K.TAU_RAD_ERR:.0f} to "
                f"{K.TAU_RAD + K.TAU_RAD_ERR:.0f}\\unit{{ns}}")
    raise KeyError(name)


def uncertainty_numbers() -> tuple:
    """Each input alone at its two ends, on both closures: the 147-nm
    share, the near IR and the final gas temperature; half the difference
    of its two ends is its term, summed in quadrature, and the excursions
    from the reference summed by sign."""
    n, lines = {}, []
    out = {}
    for eedf in K.CLOSURES:
        ref = K.run_many([("field", K.ref_params(eedf))])[0]
        items = []
        for name, end, v, p in K.band_params(eedf):
            items.append((name, end, v, p))
        for rid, end, f, p in K.const_params(eedf):
            items.append((rid, end, f, p))
        got = K.run_many([("field", p) for *_, p in items])
        res = {}
        for (name, end, v, _), s in zip(items, got):
            res[(name, end)] = s
        out[eedf] = (ref, res)
    names = ([b.name for b in K.row_bands()]
             + [rid for rid, _, _ in K.CONSTANTS])
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        ref, res = out[eedf]
        y0 = ref["emission_vuv"]
        terms, up, down = {}, 0.0, 0.0
        for name in names:
            vals = []
            for end in ("low", "high"):
                s = res.get((name, end), ref)
                vals.append(s["emission_vuv"])
                dv = s["emission_vuv"] - y0
                up += max(dv, 0.0)
                down += min(dv, 0.0)
            terms[name] = 0.5 * abs(vals[1] - vals[0])
        quad = math.sqrt(sum(t * t for t in terms.values()))
        top = sorted(terms.items(), key=lambda kv: -kv[1])
        n[f"unc.{tag}.y0"] = f"{y0:.2f}"
        n[f"unc.{tag}.quad"] = f"{quad:.1f}"
        n[f"unc.{tag}.quad.pct"] = f"{100 * quad / y0:.0f}"
        n[f"unc.{tag}.up"] = f"+{up:.1f}"
        n[f"unc.{tag}.down"] = f"$-${-down:.1f}"
        n[f"unc.{tag}.top1"] = UNC_SHORT[top[0][0]]
        n[f"unc.{tag}.top1.t"] = f"{top[0][1]:.2f}"
        n[f"unc.{tag}.top2"] = UNC_SHORT[top[1][0]]
        n[f"unc.{tag}.top2.t"] = f"{top[1][1]:.2f}"
        n[f"unc.{tag}.top3"] = UNC_SHORT[top[2][0]]
        n[f"unc.{tag}.top3.t"] = f"{top[2][1]:.2f}"
        n[f"unc.{tag}.top1.share"] = f"{100 * top[0][1] ** 2 / quad ** 2:.0f}"
    ref_k, res_k = out["auto"]
    top_k = max(names, key=lambda x: 0.5 * abs(
        res_k.get((x, "high"), ref_k)["emission_vuv"]
        - res_k.get((x, "low"), ref_k)["emission_vuv"]))
    if top_k != "kr2_form":
        raise AssertionError(f"the text says Kr* + 2Kr leads; {top_k} does")
    for name in ("kr2_pen", "kstep_s", "kstep_p"):
        cut = 1 - res_k[(name, "high")]["emission_vuv"] / ref_k["emission_vuv"]
        if not 0.15 < cut < 0.55:
            raise AssertionError(f"the text says {name} at its top cuts the "
                                 f"share by a fifth to a half: {cut:.2f}")
    b = BANDS["kr2_form"]
    n["form.ratio"] = f"{b.hi / b.lo:.1f}"
    n["unc.n"] = str(len(names))
    n["unc.nband"] = str(len(K.row_bands()))
    n["unc.nrow"] = str(len(K.CONSTANTS))
    # the table: each input, its ends, and the outputs at its two ends,
    # the largest term (kinetic) first
    ref_k, res_k = out["auto"]
    ref_t, res_t = out["two_t"]

    def term(x):
        return abs(res_k.get((x, "high"), ref_k)["emission_vuv"]
                   - res_k.get((x, "low"), ref_k)["emission_vuv"])
    for name in sorted(names, key=lambda x: -term(x)):
        if name in BANDS:
            b = BANDS[name]
            ends = _ends_text(name, b.lo, b.hi)
        else:
            ends = _ends_text(name, None, None)
        sk = [res_k.get((name, e), ref_k) for e in ("low", "high")]
        st = [res_t.get((name, e), ref_t) for e in ("low", "high")]
        lines.append(
            f"{UNC_LABELS[name]} & {ends} & "
            f"{sk[0]['emission_vuv']:.2f} to {sk[1]['emission_vuv']:.2f} & "
            f"{st[0]['emission_vuv']:.2f} to {st[1]['emission_vuv']:.2f} & "
            f"{sk[0]['emission_nir']:.2f} to {sk[1]['emission_nir']:.2f} & "
            f"{sk[0]['tg_end_K']:.0f} to {sk[1]['tg_end_K']:.0f}\\\\")
    return n, ["\\def\\wgUncTable{" + "\n".join(lines) + "}"]


# ------------------------------------------------------ the reference
def reference_numbers() -> tuple:
    n, lines = {}, []
    ref, ref2 = K.run_many([("field", K.ref_params()),
                            ("field", K.ref_params("two_t"))])
    for s, tag in ((ref, "kin"), (ref2, "2t")):
        ev = s["events"]
        krs = s["res"]["Kr*"]
        n.update({
            f"ref.{tag}.147": f"{s['emission_vuv']:.2f}",
            f"ref.{tag}.IR": f"{s['emission_nir']:.2f}",
            f"ref.{tag}.heat": f"{s['heat_reactions']:.2f}",
            f"ref.{tag}.elastic": f"{s['elastic']:.2f}",
            f"ref.{tag}.krstar": f"{s['to_Krstar']:.2f}",
            f"ref.{tag}.krplus": f"{s['to_Krplus']:.2f}",
            f"ref.{tag}.kr2": f"{s['to_Kr2star']:.2f}",
            f"ref.{tag}.ions": f"{s['stored_ions']:.2f}",
            f"ref.{tag}.te": f"{s['te_mean_eV']:.2f}",
            f"ref.{tag}.tg": f"{s['tg_end_K']:.0f}",
            f"ref.{tag}.ne": _sci(s["ne_end_cm3"]),
            f"ref.{tag}.cyc": f"{ev['K76'] / ev['AK8_Kr']:.0f}",
            f"ref.{tag}.brs76": f"{100 * krs['losses']['K76'] / krs['formed']:.1f}",
            f"ref.{tag}.ak8": f"{100 * krs['losses']['AK8_Kr'] / krs['formed']:.1f}",
            f"ref.{tag}.closure": _sci(abs(s["closure"]), 0),
            f"ref.{tag}.min": f"{s['seconds'] / 60:.0f}",
        })
        # the cycle's heat: the net pumping's 1.385 eV per event
        gap = e_form("Kr**") - e_form("Kr*")
        n[f"ref.{tag}.cycheat"] = f"{100 * (ev['K76'] - ev['K76r']) * gap / s['e_dep_eV']:.1f}"
        n[f"ref.{tag}.perphoton"] = f"{(ev['K76'] - ev['K76r']) * gap / ev['KR2_RAD']:.0f}"
    kin = ref["kin"]
    k15 = next(r for r in build() if r.rid == "K15").rate[1]
    n_kr = ref["n0"]["Kr"]
    n["k15.ns"] = f"{1e9 / (k15 * n_kr ** 2):.0f}"
    ev = ref["events"]
    c2 = ev["C2"] / sum(ev[r] for r in ("C2", "G2", "KSTEP_S", "KSTEP_P",
                                         "PEN_KR2"))
    if not c2 < 0.01:
        raise AssertionError("the text says Penning ionization makes the "
                             "reference point's ions")
    n.update({
        "ref.kin.n": str(kin["solutions"]), "ref.kin.on": str(kin["on"]),
        "ref.kin.en.lo": f"{kin['en_min']:.2f}",
        "ref.kin.en.hi": f"{kin['en_max']:.1f}",
        "ref.kin.ratio.lo": f"{kin['ratio_min']:.4f}",
        "ref.kin.ratio.hi": f"{kin['ratio_max']:.4f}",
        "ref.kin.gap": _sci(kin["balance_max"], 0),
        "ref.kin.reach": _sci(kin["reach_max"], 0),
        "ref.kin.emax.lo": f"{kin['emax_min']:.0f}",
        "ref.kin.emax.hi": f"{kin['emax_max']:.0f}",
        "ref.kin.vs2t": _pm(100 * (ref["emission_vuv"] / ref2["emission_vuv"]
                                   - 1)),
    })
    # the 5s -> 5p scale and Kannari's fit, both closures
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        jobs = [("field", K.ref_params(eedf).replace(
            rates=(("k_5s5p", k),))) if k != 1.0 else
            ("field", K.ref_params(eedf)) for k in K.K5S5P_TABLE]
        got = K.run_many(jobs)
        for k, s in zip(K.K5S5P_TABLE, got):
            n[f"k5.{k:g}.{tag}"] = f"{s['emission_vuv']:.1f}"
            n[f"k5.{k:g}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
        sw = K.switch_params(eedf)
        kann = K.run_many([("field", sw["kannari"])])[0]
        n[f"kann.{tag}"] = f"{kann['emission_vuv']:.1f}"
    # the scan, two-temperature closure
    sp = K.scan_params()
    keys = list(sp)
    got = dict(zip(keys, K.run_many([("field", sp[k]) for k in keys])))
    for w in K.SCAN_KW:
        cells = [f"{got[(a, w)]['emission_vuv']:.1f}" for a in K.SCAN_ATM]
        lines.append(f"{w:.0f} & " + " & ".join(cells) + "\\\\")
    tg01 = [got[(0.1, w)]["tg_end_K"] for w in K.SCAN_KW]
    n["scan.tg01.lo"] = f"{min(tg01):.0f}"
    n["scan.tg01.hi"] = f"{max(tg01):.0f}"
    n["scan.best"] = f"{max(s['emission_vuv'] for s in got.values()):.0f}"
    return n, ["\\def\\wgScan{" + "\n".join(lines) + "}"]


# ------------------------------------------------- numerical accuracy
def convergence_numbers() -> dict:
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    got = K.run_many([("field", p) for *_, p in K.conv_params()])
    worst = 0.0
    for (k, v, _), s in zip(K.conv_params(), got):
        d147 = abs(s["emission_vuv"] / ref["emission_vuv"] - 1)
        dir_ = abs(s["emission_nir"] / ref["emission_nir"] - 1)
        n[f"conv.{k}.147"] = _sci(d147, 0)
        n[f"conv.{k}.IR"] = _sci(dir_, 0)
        if k != "rtol":
            worst = max(worst, d147, dir_)
        else:
            n["conv.rtol.147.pct"] = f"{100 * d147:.2f}"
    n["conv.worst"] = _sci(worst, 0)
    return n


def conversion_numbers() -> dict:
    """v0_5 with no F2 against v0_1 on the two-temperature closure: with
    its own tables (the two misreadings of Jung et al. 2006) and with
    v0_1's."""
    n = {}
    M5 = importlib.import_module("v0_5.model")
    R5 = importlib.import_module("v0_5.report")
    J5 = importlib.import_module("v0_5.jung2006")
    p5 = M5.Params(x_f2=0.0, eedf="two_t")
    own = R5.row_dict(M5.run(p5))
    saved = copy.deepcopy(J5.TABLE_II), copy.deepcopy(J5.TABLE_III)

    def clear():
        for f in (J5.process, J5.rate_table, J5.superelastic_table,
                  J5.capture_table):
            f.cache_clear()
    try:
        for mine, theirs in ((J5.TABLE_II, jung2006.TABLE_II),
                             (J5.TABLE_III, jung2006.TABLE_III)):
            mine.clear()
            mine.update(copy.deepcopy(theirs))
        clear()
        same = R5.row_dict(M5.run(p5))
    finally:
        for mine, old in zip((J5.TABLE_II, J5.TABLE_III), saved):
            mine.clear()
            mine.update(old)
        clear()
    v1 = row_dict(run(Params(eedf="two_t")))
    # the keys the two reports share (v0_5's has no near-IR line of its
    # own; with no F2 its all-band emission is the 147 nm and the near IR)
    keys = ("emission_vuv", "emission_all", "heat_reactions", "elastic",
            "tg_end_K", "te_mean_eV", "to_Krstar", "to_Krplus")
    dev = max(abs(v1[k] / same[k] - 1) for k in keys)
    n["par.own"] = f"{own['emission_vuv']:.3f}"
    n["par.v01"] = f"{v1['emission_vuv']:.3f}"
    n["par.dev"] = _sci(dev, 0)
    n["par.rows"] = str(len(build()))
    return n


def jung_numbers() -> dict:
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    te_k = ref["te_mean_eV"]
    table = jung2006.rate_table()
    for te, key in ((1.0, "1"), (2.0, "2"), (3.0, "3"), (te_k, "ref")):
        r = float(table(te)) / float(jung2006.kannari_rate(te))
        n[f"jung.ratio.{key}"] = f"{r:.2f}"
    n["jung.te.ref"] = f"{te_k:.2f}"
    n["jung.sig6"] = _sci(float(jung2006.lumped_sigma(6.0)))
    n["jung.sig6.meta"] = _sci(float(jung2006.lumped_sigma(6.0,
                                                           proxy_j1=False)))
    n["jung.unc.1s5"], n["jung.unc.1s3"] = (str(x) for x in JUNG_UNC)
    n["jung.k.ref"] = _sci(float(table(te_k)))
    return n


def ladder_numbers() -> dict:
    """The energy ladder's levels and the energies its rows move [eV]."""
    rows = {r.rid: r for r in build()}
    E = {x: e_form(x) for x in ("Kr*", "Kr**", "Kr+", "Kr2*", "Kr2+")}
    hv = rows["KR2_RAD"].hv_eV
    return {
        "lad.S": f"{E['Kr*']:.3f}", "lad.P": f"{E['Kr**']:.2f}",
        "lad.I": f"{E['Kr+']:.3f}", "lad.X2": f"{E['Kr2*']:.2f}",
        "lad.I2": f"{E['Kr2+']:.3f}",
        "lad.hv": f"{hv:.2f}", "lad.low": f"{E['Kr2*'] - hv:.2f}",
        "lad.gap": f"{E['Kr**'] - E['Kr*']:.3f}",
        "lad.ak8": f"{E['Kr*'] - E['Kr2*']:.3f}",
        "lad.k15": f"{E['Kr+'] - E['Kr2+']:.2f}",
        "lad.r11": f"{E['Kr2+'] - E['Kr**']:.3f}",
        "lad.pen": f"{2 * E['Kr*'] - E['Kr+']:.3f}",
        "lad.c1hi": f"{rows['C1c'].e_cost_eV:.2f}",
    }


def numbers() -> tuple:
    # every run at once, in parallel, so a stale cache is one batch
    K.run_many(K.all_jobs())
    n = {"version": K.V, "md5": K.code_md5()[:12]}
    tables = []
    for part in (reference_numbers, uncertainty_numbers, swarm_numbers):
        a, t = part()
        n.update(a)
        tables += t
    for part in (beam_numbers, convergence_numbers, conversion_numbers,
                 jung_numbers, ladder_numbers):
        n.update(part())
    n["tau5p"] = f"{KR_5P_TAU_S * 1e9:.1f}"
    return n, tables


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_1 overview numbers. GENERATED by make_overview_v0_1.py",
        f"% from the model ({n['version']}). Do not edit. Rerun the script "
        "instead.",
        "\\providecommand\\wgdef[2]{\\expandafter\\def"
        "\\csname wg@#1\\endcsname{#2}}",
    ]
    lines += [f"\\wgdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_1's accuracy "
                                 "overview numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_overview_v0_1.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
