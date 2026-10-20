"""Numbers for 'How accurate LIGHT_Kr v0_2 is' (accuracy_overview_v0_2.tex).

Writes overview_v0_2_numbers.tex, one \\wgdef{key}{value} per number and a
few \\def tables, from the model:

  the swarm      v0_2's two-term solver run as a swarm (v0_2/swarm.py, the
                 kinetic EEDF's own solver) on its Biagi 8.97 cross
                 sections, without electron-electron collisions, against
                 krypton's measured drift velocity (Pack, Voshall & Phelps
                 1962) and ionization coefficient (Kruithof 1940, Heylen
                 1971), as Dutton's 1975 survey prints them, with direct
                 ionization alone (the default) and with every excitation
                 of the levels at or above E(Kr2+) ionizing associatively
                 (the band hm_kr at its top). v0_1's swarm (its Biagi-v7.1
                 groups, make_overview_v0_1.swarm_point) beside it, and
                 LoKI-MC's Monte Carlo solution of v0_2's cross sections
                 as a check on the two-term solution (LOKI_MC, typed below
                 with the runs' setup). The data are typed from Dutton's
                 tables (DUTTON_W, DUTTON_ALPHA), printed values, no
                 figure readings.
  the e-beam     v0_2 under the beam drive (kr_docs_v0_2.BeamModel) against
                 Morozov et al. 2008 (0.42 +- 0.05 above 500 hPa) and
                 Eckstrom et al. 1988 (near 46 % up to 3e5 W/cm3), and the
                 drive's partition against Saito et al. 2003 (W_s 15.4 eV,
                 N_ex/N_i 0.55).
  uncertainty    every band the rows read at its two ends and the three
                 measured constants with a quoted error at theirs, one at a
                 time, on both closures (kr_docs_v0_2.band_params,
                 const_params).
  the reference  the kinetic EEDF and the two-temperature closure, the
                 5s -> 5p scale, Kannari's fit, the scan, the kinetic
                 closure's convergence (kr_docs_v0_2.conv_params), and v0_1
                 at the same points (kr_docs, its own cache).

    python3 LIGHT_Kr/docs/make_overview_v0_2.py        # minutes (cached)
    cd LIGHT_Kr/docs && latexmk -pdf accuracy_overview_v0_2.tex

The runs come from kr_docs_v0_2.run_many (cached per version of the
model's code). --check compares the file with a fresh computation and
writes nothing.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kr_docs as K1                                       # noqa: E402
import kr_docs_v0_2 as K                                   # noqa: E402
import make_overview_v0_1 as OV1                           # noqa: E402
from v0_2 import jung2006, swarm                           # noqa: E402
from v0_2.network import (BANDS, E_ASSOC_EV, KR_5P_TAU_S,  # noqa: E402
                          build, group_rows, levels)
from v0_2.species import e_form                            # noqa: E402

OUT = HERE / "overview_v0_2_numbers.tex"

# the measurements, the scores' ranges and the e-beam's: v0_1's overview's,
# typed there from their sources (make_overview_v0_1.py)
DUTTON_W, DUTTON_ALPHA = OV1.DUTTON_W, OV1.DUTTON_ALPHA
ALPHA_RANGES, W_RANGES = OV1.ALPHA_RANGES, OV1.W_RANGES
MOROZOV, ECKSTROM = OV1.MOROZOV, OV1.ECKSTROM
SAITO_WS, SAITO_NEX = OV1.SAITO_WS, OV1.SAITO_NEX
JUNG_UNC = OV1.JUNG_UNC
CELLS_SWARM = OV1.CELLS_SWARM
# LoKI-MC (Dias et al., Comput. Phys. Commun. 282, 108554 (2023); the
# vendored LoKI-MC/Code/lokimc) on v0_2's Biagi 8.97 file rewritten in
# LoKI-MC's block form (as v0_2/tests/test_kinetic.py's _loki_input does),
# pure Kr at 300 K and 1 Torr, 1e5 electrons, equal sharing, 2026-10-14:
# E/N [Td] -> (the SST reduced Townsend coefficient [m^2], the mean energy
# [eV], the ionizations counted after the steady state)
LOKI_MC = {20.0: (3.6057e-25, 4.5844, 29), 30.0: (7.1147e-24, 4.7536, 543),
           40.0: (3.3102e-23, 4.9022, 3127), 50.0: (8.7447e-23, 5.0360, 9521)}
# Hornbeck & Molnar, Phys. Rev. 84, 621 (1951), abstract: the molecular
# ions' appearance potentials lie below the atomic ions' by "0.7-0.3+0.7
# volt in Kr" (the energy below Kr+'s at which Kr2+ appears, V)
HORNBECK_KR = (0.7, 0.3, 0.7)
_f, _pm, _sci, _rms = OV1._f, OV1._pm, OV1._sci, OV1._rms


# ------------------------------------------------------------- the swarm
def swarm_point(en_td: float, tg: float, cells: int = CELLS_SWARM,
                hm: float = 0.0) -> dict:
    """v0_2's swarm (v0_2/swarm.py): W [cm/s], alpha/N [cm2], the mean
    energy [eV]."""
    return swarm.point(en_td, tg, hm=hm, cells=cells)


def _alpha_scores(points, key="alphaN"):
    """((E/N, source, error, measured, model), ...) for Dutton's alpha/N."""
    return [(en, src, d[key] / am - 1.0, am, d[key])
            for (en, am, src), d in zip(DUTTON_ALPHA, points)]


def swarm_numbers() -> tuple:
    n = {}
    ew = []
    for en, wm, tg in DUTTON_W:
        d = swarm_point(en, tg)
        ew.append((en, tg, d["W"] / wm - 1.0, wm, d["W"]))
    direct = [swarm_point(en, 300.0) for en, _, _ in DUTTON_ALPHA]
    top = [swarm_point(en, 300.0, hm=1.0) for en, _, _ in DUTTON_ALPHA]
    old = [OV1.swarm_point(en, 300.0) for en, _, _ in DUTTON_ALPHA]
    ea = _alpha_scores(direct)
    eh = _alpha_scores(top)
    e1 = _alpha_scores(old)
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
    w1 = [OV1.swarm_point(en, tg)["W"] / wm - 1.0 for en, wm, tg in DUTTON_W]
    n["sw1.W.rms"] = _f(100 * _rms(w1))
    for tag, rows in (("sw", ea), ("swhm", eh), ("sw1", e1)):
        for i, (a, b) in enumerate(ALPHA_RANGES):
            sel = [x[2] for x in rows if a <= x[0] < b]
            n[f"{tag}.a.r{i}.n"] = str(len(sel))
            n[f"{tag}.a.r{i}.rms"] = _f(100 * _rms(sel))
            n[f"{tag}.a.r{i}.mean"] = _pm(100 * np.mean(sel))
            for src in ("K", "H"):
                ss = [x[2] for x in rows if a <= x[0] < b and x[1] == src]
                n[f"{tag}.a.r{i}.{src}"] = (_pm(100 * np.mean(ss)) if ss
                                            else "--")
    run_pts = sorted((en, wm) for en, wm, tg in DUTTON_W
                     if tg == 300 and W_RANGES[1][0] <= en < W_RANGES[1][1])
    n["sw.W.r1.wx"] = f"{run_pts[-1][1] / run_pts[0][1]:.0f}"
    n["sw.W.r1.ex"] = f"{run_pts[-1][0] / run_pts[0][0]:.1f}"
    first300 = min(en for en, _, _ in DUTTON_ALPHA if en >= 300.0)
    n["sw.a.r3.eps"] = f"{swarm_point(first300, 300.0)['mean']:.0f}"
    low = [x for x in ea if x[0] < 50.0]
    raw = {"r0": np.mean([x[2] for x in low]),
           "r3": np.mean([x[2] for x in ea if x[0] >= 300.0]),
           "v01": np.mean([x[2] for x in e1 if x[0] < 50.0])}
    if not (raw["r0"] < 0.0 < raw["r3"]):
        raise AssertionError("the text says the model falls short below 50 "
                             "Td and lies above past 300 Td")
    if not raw["v01"] < raw["r0"]:
        raise AssertionError("the text says v0_2 is closer than v0_1 below "
                             "50 Td")
    n["sw.a.n"] = str(len(ea))
    # how far the two measurements part below 50 Td: Kruithof's, read in
    # log-log between his points, over Heylen's at Heylen's E/N
    kx = np.log([en for en, _, s_ in DUTTON_ALPHA if s_ == "K"])
    ky = np.log([am for _, am, s_ in DUTTON_ALPHA if s_ == "K"])
    kh = [math.exp(float(np.interp(math.log(en), kx, ky))) / am
          for en, am, s_ in DUTTON_ALPHA if s_ == "H" and en < 50.0]
    n["kh.max"] = f"{max(kh):.1f}"
    lo14 = min(ea, key=lambda x: x[0])
    n["sw.a.lo.en"] = f"{lo14[0]:.1f}"
    n["sw.a.lo.x"] = f"{lo14[4] / lo14[3]:.3f}"
    # the band hm_kr: where the measurements lie between its two ends
    pts = [(en, src, am, d["alphaN"], t["alphaN"])
           for (en, am, src), d, t in zip(DUTTON_ALPHA, direct, top)
           if en < 50.0]
    inside = [p for p in pts if p[3] <= p[2] <= p[4]]
    above = [p for p in pts if p[2] > p[4]]
    below = [p for p in pts if p[2] < p[3]]
    if below or any(p[1] != "K" for p in above) or \
            not all(p[0] < 18.0 for p in above):
        raise AssertionError("the text says every measurement below 50 Td "
                             "lies above the direct ionization, and only "
                             "Kruithof's below 18 Td above the band's top")
    n["hm.n"] = str(len(pts))
    n["hm.inside"] = str(len(inside))
    n["hm.above"] = str(len(above))
    n["hm.above.en"] = f"{max(p[0] for p in above):.2f}" if above else "--"
    n["hm.above.lo"] = (f"{min(p[4] / p[2] for p in above):.2f}"
                        if above else "--")
    n["hm.above.hi"] = (f"{max(p[4] / p[2] for p in above):.2f}"
                        if above else "--")
    for en, key in ((14.12, "14"), (30.0, "30"), (100.0, "100")):
        r = swarm_point(en, 300.0, hm=1.0)["alphaN"] / \
            swarm_point(en, 300.0)["alphaN"]
        n[f"hm.x{key}"] = f"{r:.1f}"
    n["hm.levels"] = str(sum(1 for p in levels() if p.threshold >= E_ASSOC_EV))
    n["hm.elo"] = f"{min(p.threshold for p in levels() if p.threshold >= E_ASSOC_EV):.2f}"
    n["hm.ehi"] = f"{max(p.threshold for p in levels()):.1f}"
    n["hm.eai"] = f"{E_ASSOC_EV:.2f}"
    v, lo, hi = HORNBECK_KR
    n["hm.ap"] = f"${v:g}^{{+{hi:g}}}_{{-{lo:g}}}$"
    # the Monte Carlo check of the two-term solution, 20 to 50 Td
    dev, mdev = [], []
    for en, (a_mc, m_mc, counted) in sorted(LOKI_MC.items()):
        d = swarm_point(en, 300.0)
        dev.append((en, d["alphaN"] * 1e-4 / a_mc - 1.0, counted))
        mdev.append(abs(d["mean"] / m_mc - 1.0))
    n["mc.max"] = f"{100 * max(abs(x[1]) for x in dev):.1f}"
    n["mc.mean.max"] = f"{100 * max(mdev):.2f}"
    n["mc.list"] = ", ".join(_pm(100 * x[1], 1) for x in dev)
    n["mc.stat20"] = f"{100 / math.sqrt(LOKI_MC[20.0][2]):.0f}"
    n["mc.stat50"] = f"{100 / math.sqrt(LOKI_MC[50.0][2]):.0f}"
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
    for tg, key, topw in ((300.0, "Wa", 4.0), (368.0, "Wb", 0.012),
                          (195.0, "Wc", 0.2)):
        p_ = " ".join(f"({en:.5g},{swarm_point(en, tg)['W']:.5g})"
                      for en in ens if en <= topw)
        lines.append(f"\\def\\wgSw{key}model{{{p_}}}")
    for tg, key in ((300, "Wa"), (368, "Wb"), (195, "Wc")):
        p_ = " ".join(f"({en:.5g},{wm:.5g})" for en, wm, t in DUTTON_W
                      if t == tg)
        lines.append(f"\\def\\wgSw{key}data{{{p_}}}")
    ens = np.geomspace(12.0, 1100.0, 50)
    for key, hm in (("Amodel", 0.0), ("AHMmodel", 1.0)):
        p_ = " ".join(f"({en:.5g},{swarm_point(en, 300.0, hm=hm)['alphaN']:.5g})"
                      for en in ens)
        lines.append(f"\\def\\wgSw{key}{{{p_}}}")
    p_ = " ".join(f"({en:.5g},{OV1.swarm_point(en, 300.0)['alphaN']:.5g})"
                  for en in ens)
    lines.append(f"\\def\\wgSwAoldmodel{{{p_}}}")
    for src in ("K", "H"):
        p_ = " ".join(f"({en:.5g},{am:.5g})" for en, am, s in DUTTON_ALPHA
                      if s == src)
        lines.append(f"\\def\\wgSwA{src}data{{{p_}}}")
    p_ = " ".join(f"({en:.5g},{a * 1e4:.5g})"
                  for en, (a, _, _) in sorted(LOKI_MC.items()))
    lines.append(f"\\def\\wgSwAMCdata{{{p_}}}")
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
            exc = sum(v for r, v in ev.items()
                      if r.startswith(("G1_", "C1")))
            n[f"{name}.{tag}.G1"] = f"{exc / pairs:.3f}"
            n[f"{name}.{tag}.te"] = f"{s['te_mean_eV']:.2f}"
            n[f"{name}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
            n[f"{name}.{tag}.ph"] = f"{ev['KR2_RAD'] / pairs:.3f}"
            n[f"{name}.{tag}.rec"] = f"{ev['R11p'] / pairs:.3f}"
            n[f"{name}.{tag}.krs"] = _sci(s["n_max"]["Kr*"])
            n[f"{name}.{tag}.ne"] = _sci(s["ne_peak"])
            n[f"{name}.{tag}.closure"] = _sci(abs(s["closure"]), 0)
            # the loss below v0_2's ceiling (one photon per excitation and
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
            raise AssertionError("the text says v0_2 lies above Morozov")
        if not 1 - got[("eckstrom", 2.0)]["emission_vuv"] / 100 / ceil < 0.10:
            raise AssertionError("the text says v0_2 loses little at 2 atm")
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
    "hm_kr": r"Kr(nl) + Kr \ra \Krtp{} + \el{}, the share",
    "R11p": r"\Krtp{} + \el{} \ra \Krpp{} + Kr",
    "K15": r"\Krp{} + 2Kr \ra \Krtp{} + Kr",
    "KR2_RAD": r"\Krts{} lifetime",
}
UNC_SHORT = {
    "k_5s5p": r"the 5s \ra 5p cross section", "kps_kr": r"\Krpp{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{}", "kr2_form": r"\Krs{} + 2Kr",
    "kr2_pen": r"\Krts{} + \Krts{}", "kstep_s": r"\el{} + \Krs{} \ra \Krp",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp", "R11p": r"\Krtp{} + \el",
    "hm_kr": r"the Hornbeck--Molnar share",
    "K15": r"\Krp{} + 2Kr", "KR2_RAD": r"the \Krts{} lifetime",
}


def _ends_text(name, lo, hi):
    if name in BANDS:
        b = BANDS[name]
        if b.unit in ("factor", "fraction"):
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
    # v0_2: the Hornbeck-Molnar share at its top barely moves the
    # reference point (direct and associative ionization make well under
    # a percent of its ions next to Penning's)
    hm_cut = res_k[("hm_kr", "high")]["emission_vuv"] / ref_k["emission_vuv"] - 1
    if not abs(hm_cut) < 0.01:
        raise AssertionError(f"the text says hm_kr at its top moves the "
                             f"share by under 1 %: {hm_cut:+.4f}")
    n["unc.hm.cut"] = _pm(100 * hm_cut, 2)
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
        # the cycle's heat: the net pumping's 1.388 eV per event
        gap = e_form("Kr**") - e_form("Kr*")
        n[f"ref.{tag}.cycheat"] = f"{100 * (ev['K76'] - ev['K76r']) * gap / s['e_dep_eV']:.1f}"
        n[f"ref.{tag}.perphoton"] = f"{(ev['K76'] - ev['K76r']) * gap / ev['KR2_RAD']:.0f}"
    kin = ref["kin"]
    k15 = next(r for r in build() if r.rid == "K15").rate[1]
    n_kr = ref["n0"]["Kr"]
    n["k15.ns"] = f"{1e9 / (k15 * n_kr ** 2):.0f}"
    ev = ref["events"]
    hm_ev = sum(v for r, v in ev.items() if r.startswith("HM_"))
    c2 = (ev["C2"] + hm_ev) / (sum(ev[r] for r in ("C2", "G2", "KSTEP_S",
                                                   "KSTEP_P", "PEN_KR2"))
                               + hm_ev)
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
            # v0_2: the renewals of the tighter run, against the reference's
            n["conv.rtol.n"] = str(s["kin"]["solutions"])
    n["conv.worst"] = _sci(worst, 0)
    return n


def change_numbers() -> dict:
    """v0_1 (Biagi-v7.1's groups, its own cache, kr_docs) against v0_2 at
    the reference point on both closures, and the network's size."""
    n = {}
    old = K1.run_many([("field", K1.ref_params()),
                       ("field", K1.ref_params("two_t"))])
    new = K.run_many([("field", K.ref_params()),
                      ("field", K.ref_params("two_t"))])
    for (o, s), tag in zip(zip(old, new), ("kin", "2t")):
        n[f"chg.{tag}.v01"] = f"{o['emission_vuv']:.2f}"
        n[f"chg.{tag}.v02"] = f"{s['emission_vuv']:.2f}"
        n[f"chg.{tag}.dev"] = _pm(100 * (s["emission_vuv"]
                                         / o["emission_vuv"] - 1), 1)
        n[f"chg.{tag}.te.v01"] = f"{o['te_mean_eV']:.2f}"
        n[f"chg.{tag}.te.v02"] = f"{s['te_mean_eV']:.2f}"
        n[f"chg.{tag}.el.v01"] = f"{o['elastic']:.2f}"
        n[f"chg.{tag}.el.v02"] = f"{s['elastic']:.2f}"
    n["par.rows"] = str(len(build()))
    n["par.levels"] = str(len(levels()))
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
        "lad.c1hi": f"{max(r.e_cost_eV for r in group_rows(build(), 'C1c')):.2f}",
        "lad.c1lo": f"{min(r.e_cost_eV for r in group_rows(build(), 'C1a')):.2f}",
        "lad.s1": f"{min(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
        "lad.s2": f"{max(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
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
    for part in (beam_numbers, convergence_numbers, change_numbers,
                 jung_numbers, ladder_numbers):
        n.update(part())
    n["tau5p"] = f"{KR_5P_TAU_S * 1e9:.1f}"
    return n, tables


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_2 overview numbers. GENERATED by make_overview_v0_2.py",
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
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_2's accuracy "
                                 "overview numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_overview_v0_2.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
