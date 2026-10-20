"""Numbers, tables and figures for accuracy_overview.tex, LIGHT_ArF_ebeam's accuracy.

No number in the overview is typed by hand. This script writes three files,
each one \\widef{key}{value} per number and one \\def\\wiT<Name>{...} per
finished table, and the figures:

  overview_rates.tex   Appendix A: v0_2's rows, bands and species, from
                       v0_2/network.py and v0_2/species.py
  overview_scores.tex  the scores: arf_data/score_0d.json and score_1d.json
                       through arf_data/make_score_report.py's statistics,
                       so SCORE.md and the overview cannot disagree
  overview_runs.tex    fresh runs: the reference point with its Maxwellian
                       and seed checks, each band at both ends, Morozov's
                       point across the Ar* and Ar2* bands, R11p held at
                       30 eV at Rokni's pressures, and 1Dv0_2's reference
                       point on 64 and 128 cells with its Rokni 2-atm trace
  figures/*.pdf        from the scores, and the trace (with --runs)

From Kinetic_Modeling/:
    python3 LIGHT_ArF_ebeam/docs/make_overview.py                    # rates, scores, figures
    ~/.venvs/f2sweep/bin/python LIGHT_ArF_ebeam/docs/make_overview.py --runs
                                               # + the runs, ~20 min, Fluid1D
    python3 LIGHT_ArF_ebeam/docs/make_overview.py --check scores     # compare only
    cd LIGHT_ArF_ebeam/docs && /opt/homebrew/bin/latexmk -pdf accuracy_overview.tex

The runs import 1Dv0_2, which finds Fluid1D through FLUID1D_SRC or the
nathaniel/ copy (Python >= 3.11 with numba). v0_2/tests/test_docs.py
compares overview_rates.tex and overview_scores.tex with fresh ones, and
1Dv0_2/tests/test_docs.py overview_runs.tex (LIGHT_DOCS_SLOW=1).
"""
from __future__ import annotations

import argparse
import importlib
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

DOCS = Path(__file__).resolve().parent
LIGHT_ARF = DOCS.parent
ARF_DATA = DOCS / "arf_data"
FIGS = DOCS / "figures"
OUT_RATES = DOCS / "overview_rates.tex"
OUT_SCORES = DOCS / "overview_scores.tex"
OUT_RUNS = DOCS / "overview_runs.tex"
for _p in (str(ARF_DATA), str(LIGHT_ARF)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import make_score_report as MSR  # noqa: E402
from v0_2 import __version__ as VERSION  # noqa: E402
from v0_2 import network as N  # noqa: E402
from v0_2.species import SPECIES  # noqa: E402

HDR = r"\providecommand\widef[2]{\expandafter\def\csname wi@#1\endcsname{#2}}"


# ---------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------
def sci(v, sig=2):
    """v for LaTeX: plain for 1e-3 <= |v| < 1e4, else m x 10^e."""
    v = float(v)
    if v == 0.0:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    m = round(v / 10 ** e, sig - 1)
    if abs(m) >= 10.0:
        e += 1
        m = round(v / 10 ** e, sig - 1)
    if -3 <= e < 4:
        return f"{v:.{max(sig - 1 - e, 0)}f}"
    return f"${m:.{sig - 1}f}{{\\times}}10^{{{e}}}$"


def fx(v, nd=2):
    """Fixed decimals; a negative number gets a math minus (text mode)."""
    v = float(v)
    s = f"{abs(v):.{nd}f}"
    return ("$-$" + s) if v < 0 and float(s) != 0.0 else s


def rng(lo, hi, nd=1):
    """'lo to hi' at nd decimals (the prose rule: no dash ranges)."""
    return f"{fx(lo, nd)} to {fx(hi, nd)}"


SP_TEX = {"e": r"\el{}", "Ar": "Ar", "F2": r"\Ft{}", "F": "F", "F-": r"\Fm{}",
          "Ar*": r"\Ars{}", "Ar+": r"\Arp{}", "Ar2+": r"\Artp{}",
          "Ar2*": r"\Arts{}", "ArF(B)": "ArF(B)", "Ar2F*": r"\AtF{}"}


def side(names) -> str:
    counts = {}
    for s in names:
        counts[s] = counts.get(s, 0) + 1
    return " + ".join((f"{k}" if k > 1 else "") + SP_TEX[s]
                      for s, k in counts.items())


def reaction(r) -> str:
    rhs = side(r.products) + (r" + \hnu{}" if r.is_rad else "")
    if r.rate[0] == "elastic":
        return side(r.reactants) + ", elastic"
    return side(r.reactants) + r" \ra{} " + rhs


def tex_name(s) -> str:
    return s.replace("_", r"\_")


def tex_free(s) -> str:
    """Free text (a figure reading's header) as LaTeX text."""
    s = s.replace("\\", "/").replace(" & ", " and ").replace("&", " and ")
    s = s.replace("; ", ", ")
    for a, b in (("%", r"\%"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}"), ("$", r"\$")):
        s = s.replace(a, b)
    return s.replace("+-", r"$\pm$")


# text of Band.what, plain to LaTeX (longest names first)
_TEXT = ((" & ", " and "), ("Ar2F*", r"\AtF{}"), ("ArF*", r"\ArFs{}"), ("Ar2*", r"\Arts{}"),
         ("Ar2+", r"\Artp{}"), ("Ar*", r"\Ars{}"), ("Ar+", r"\Arp{}"),
         ("F2", r"\Ft{}"), ("->", r"\ra{}"), ("k_e tau_B", r"$k_e\tau_B$"),
         ("k tau_B", r"$k\tau_B$"), ("(Tg/300)^n", r"$(T_g/300)^n$"),
         ("hv", r"\hnu{}"))


def tex_text(s) -> str:
    for a, b in _TEXT:
        s = s.replace(a, b)
    return re.sub(r"\^(\d+)", r"$^{\1}$", s)       # Ar^2 -> Ar$^{2}$


UNIT_K = {1: r"s$^{-1}$", 2: r"cm$^3$s$^{-1}$", 3: r"cm$^6$s$^{-1}$"}
UNIT_KT = {2: r"cm$^3$", 3: r"cm$^6$"}
UNIT_BAND = {"cm3/s": r"cm$^3$s$^{-1}$", "cm6/s": r"cm$^6$s$^{-1}$",
             "cm6": r"cm$^6$", "cm3": r"cm$^3$", "1/s": r"s$^{-1}$",
             "fraction": "", "exponent": ""}


def render(title, n: dict, tables: dict) -> str:
    lines = [f"% {title}. GENERATED by make_overview.py from LIGHT_ArF_ebeam",
             f"% v0_2 {VERSION}. Do not edit. Rerun the script instead.", HDR]
    lines += [f"\\widef{{{k}}}{{{v}}}" for k, v in n.items()]
    for name, rows in tables.items():
        assert name.isalpha(), name
        lines.append(f"\\def\\wiT{name}{{%")
        lines += [r + "%" if not r.endswith("\\\\") else r for r in rows]
        lines.append("}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Appendix A: the rates
# ---------------------------------------------------------------------------
# short sources, held to each row's (and its bands') citation by
# v0_2/tests/test_docs.py: every name here appears there
SHORT = {
    "G1": "IST-Lisbon set, Yanguas-Gil, Cotrino and Alves 2005",
    "C1": "IST-Lisbon set, prompt cascade to 4s after Setser 1977",
    "C2": "IST-Lisbon set, Rapp and Englander-Golden",
    "G6": "IST-Lisbon momentum transfer, after Phelps",
    "ATT": "SIGLO", "DIS_a3": "SIGLO", "DIS_A1": "SIGLO", "VIB": "SIGLO",
    "G2": "no measurement, band top from Lee and Lieberman",
    "G3": "total the mean of Velazco, Chen and Ueda, ArF$^*$ yield 0.5 "
          "measured by Setser 1977",
    "G3X": "the rest of G3, products assigned",
    "AK8_Ar": "Ueda 1984, band Keto 1974 to Kolts and Setser 1977",
    "G4": r"$\tau_B$ calculated by Dunning and Hay, band peak Chen and Payne",
    "QF2": "mean of Rokni 1977 and Chen and Payne 1979",
    "QAR": "Rokni 1977",
    "G5": "mean of Rokni, Chen and Payne, and Jancaitis, top Marowsky",
    "QE": "Duzy 1980",
    "C3": "mean of Chen and Payne, Jancaitis, and Ueda",
    "Q2F": "mean of Chen and Payne, Jancaitis, and Ueda",
    "Q2FF": "0.5 to 1 times Q2F, Jancaitis 1983",
    "Q2AR": "below $10^{-14}$, Jancaitis, top Boewering",
    "K15": "median of eight, Brau's Table 4.8, $T_g$ law Perez-Rios and "
           "Greene",
    "R11p": "Oskam and Mittelstadt, Mehr and Biondi, via Brau",
    "K26": "Hoffman and Moreno's calculation as Lee 1989 prints it",
    "K108": "the same form as K26",
    "FF": "Ganguli and Kaufman via Baulch",
    "AR2_RAD": "Keto, Gleason and Walters 1974, top from Kannari",
    "NK7": "Chen et al. via Chen and Payne, share Setser over Chen and Payne",
    "K2F_BR": "the rest of NK7",
    "K2F_F": "no measurement, band top from SRI",
    "PEN_AR2": "no measurement, band top from Lorents and Olson",
    "KSTEP": "REBUILD's template, an estimate",
    "SUPER": "detailed balance of the G1 rows, an estimate",
}
BAND_SHORT = {
    "ar_f2": "mean of Velazco, Chen and Ueda, ends Ueda 1984 and 1983",
    "arf_yield": "Setser 1977, bottom Chen and Payne, top the published "
                 "models",
    "ar2_form": "Ueda 1984, ends Keto 1974 and Kolts and Setser 1977",
    "ar2_rad": "Keto 1974, top Kannari's mixing",
    "ar2_f2": "Chen et al., bottom SRI",
    "ar2_f2_arf": "Setser over Chen and Payne",
    "ar2_f": "no measurement, top SRI",
    "ar2_pen": "no measurement, top Kannari",
    "ar_pen": "no measurement, top Lee and Lieberman",
    "arf_trimer": "mean of Rokni, Chen and Payne, and Jancaitis, top "
                  "Marowsky",
    "trimer_tg_pow": "no measurement, bottom NRL's form after Shui",
    "arf_e": "Duzy 1980, two measurements",
    "ar2f_f": "Jancaitis 1983",
    "ar2f_ar": "below $10^{-14}$, Jancaitis and Marowsky, top Boewering",
}
SPECIES_SHORT = {
    "Ar": ("", "ground state"),
    "F2": ("", "ground state"),
    "F": ("measured", r"half the F$_2$ bond, REBUILD"),
    "F-": ("measured", "F less the electron affinity, REBUILD"),
    "Ar*": ("threshold", "IST-Lisbon 3P2 threshold"),
    "Ar+": ("threshold", "IST-Lisbon ionization threshold"),
    "Ar2+": ("measured", "well depth by photoionization, via Krauss and Mies"),
    "Ar2*": ("calculated", "well depth of Saxon and Liu"),
    "ArF(B)": ("lower bound", "F plus the 1933-\\AA{} band peak of Chen and "
                              "Payne"),
    "Ar2F*": ("calculated", "well depth of Huestis, Marowsky and Tittel"),
}


def short_key(rid) -> str:
    if rid.startswith(("G1_", "C1_")):
        return rid[:2]
    if rid.startswith("VIB"):
        return "VIB"
    if rid.startswith("SUPER_"):
        return "SUPER"
    return rid


def rate_tex(r) -> str:
    kind = r.rate[0]
    nr = len(r.reactants)
    if kind == "xsec":
        return r"$\langle\sigma v\rangle$"
    if kind == "elastic":
        return r"energy loss"
    if kind == "const":
        return f"{sci(r.rate[1], 3)} {UNIT_K[nr]}"
    if kind == "inv_tau":
        if len(r.rate) > 2:
            return (f"$k\\tau_B$ {sci(r.rate[2], 3)} {UNIT_KT[nr]}, "
                    r"over $\tau_B$")
        return rf"$1/\tau_B$, $\tau_B$ {N.ARF_TAU_B * 1e9:.1f} ns"
    if kind == "param":
        name = r.rate[1]
        frac = float(r.rate[2]) if len(r.rate) > 2 else 1.0
        mods = r.rate[3] if len(r.rate) > 3 else ()
        how = "band " + tex_name(name)
        for m in mods:
            if m[0] == "share":
                how += r" $\times$ " + tex_name(m[1])
            elif m[0] == "rest":
                how += r" $\times$ (1 $-$ " + tex_name(m[1]) + ")"
        if any(m[0] == "per_tau" for m in mods):
            kt = frac * N.BANDS[name].default
            return (f"$k\\tau_B$ {sci(kt, 3)} {UNIT_KT[nr]}, {how}, "
                    r"over $\tau_B$")
        v = N.param_rate(r.rate, lambda b: N.BANDS[b].default, N.ARF_TAU_B)
        return ("0, " + how) if v == 0 else f"{sci(v, 2)} {UNIT_K[nr]}, {how}"
    if kind == "te_pow":
        k0, _, n = r.rate[1:4]
        return (f"{sci(k0, 2)} $(T_e/300\\,\\mathrm{{K}})^{{{n:g}}}$ "
                f"{UNIT_K[nr]}")
    if kind == "hm":
        a, b, c = r.rate[1:4]
        return (f"{sci(a, 3)} $N$/($N^2$ + {b:g}$N$ + {c:g}) {UNIT_K[nr]}, "
                "$N$ in amagat")
    if kind == "stepwise":
        return (r"$1.56{\times}10^{-7}$ $T_e^{0.71}$ "
                f"$e^{{-{r.rate[1]:g}/T_e}}$ {UNIT_K[nr]}")
    if kind == "superelastic":
        return f"detailed balance, $g$ ratio 1/{round(1 / r.rate[3])}"
    raise ValueError(f"{r.rid}: no text for rate kind {kind!r}")


def tdep_tex(r) -> str:
    kind = r.rate[0]
    if kind == "xsec":
        return (r"$T_e$, $T_\mathrm{tail}$" if r.rate[1] == "Ar"
                else r"$T_e$")
    if kind == "elastic":
        return r"$T_e$, $T_g$"
    if kind == "te_pow":
        return f"$T_e^{{{r.rate[3]:g}}}$"
    if kind == "hm":
        return r"$N$, screening, $(T_g/300)^{-3}$ on the 3-body limit"
    if kind in ("stepwise", "superelastic"):
        return r"$T_e$"
    if r.tg_pow:
        return f"$(T_g/300)^{{{r.tg_pow:g}}}$"
    return ""


def _row_line(rid_tex, reac, rate, tdep, src) -> str:
    return f"{rid_tex} & {reac} & {rate} & {tdep} & {src}\\\\"


def rate_rows(rows) -> list:
    """One line per row, the Ar excitation families and F2's vibration
    folded into one line each."""
    out, done = [], set()
    for r in rows:
        key = short_key(r.rid)
        if key in ("G1", "C1", "VIB"):
            if key in done:
                continue
            done.add(key)
            fam = [q for q in rows if short_key(q.rid) == key]
            thr = [q.e_cost_eV for q in fam]
            if key == "VIB":
                reac = r"\el{} + \Ft{} \ra{} \el{} + \Ft{}($v$ = 1 to 4)"
            elif key == "G1":
                reac = (r"\el{} + Ar \ra{} \el{} + \Ars{}, the four 4s "
                        f"levels, {thr[0]:.3f} to {thr[-1]:.3f} eV")
            else:
                reac = (r"\el{} + Ar \ra{} \el{} + Ar(4p, 3d, 5p, 5s, 4d, "
                        r"6s) \ra{} \Ars{}, " f"{len(fam)} levels, "
                        f"{thr[0]:.3f} to {thr[-1]:.3f} eV")
            ids = f"{tex_name(fam[0].rid)} to {tex_name(fam[-1].rid)}"
            out.append(_row_line(ids, reac, rate_tex(r), tdep_tex(r),
                                 SHORT[key]))
            continue
        out.append(_row_line(tex_name(r.rid), reaction(r), rate_tex(r),
                             tdep_tex(r), SHORT[key]))
    return out


def band_rows() -> list:
    out = []
    for b in N.BANDS.values():
        u = UNIT_BAND[b.unit]
        f = (lambda v: fx(v, 2)) if b.unit in ("fraction", "exponent") \
            else (lambda v: sci(v, 3))
        out.append(f"{tex_name(b.name)} & {tex_text(b.what)} & "
                   f"{f(b.default)} & {f(b.lo)} to {f(b.hi)} & {u} & "
                   f"{BAND_SHORT[b.name]}\\\\")
    return out


def species_rows() -> list:
    out = []
    for s in SPECIES.values():
        status, src = SPECIES_SHORT[s.name]
        out.append(f"{SP_TEX[s.name]} & {fx(s.e_form_eV, 3)} & {status} & "
                   f"{src}\\\\")
    return out


def rates_parts():
    basic = N.build()
    group = N.electron_group_rows()
    n = {"version": VERSION.replace("_", r"\_"),
         "rates.nBasic": len(basic), "rates.nGroup": len(group),
         "rates.nBands": len(N.BANDS), "rates.nAr": len(N.AR_LEVELS),
         "rates.nAr4s": sum(1 for r in basic if r.rid.startswith("G1_")),
         "rates.nArHi": sum(1 for r in basic if r.rid.startswith("C1_")),
         "rates.nBandRows": sum(1 for r in basic if N.bands_of(r)),
         "rates.nZero": sum(1 for b in N.BANDS.values() if b.default == 0),
         "rates.tauB": fx(N.ARF_TAU_B * 1e9, 1)}
    from v0_2 import xsec
    n["rates.ebreak"] = fx(xsec.e_break(), 3)
    n["rates.temax"] = fx(xsec.TE_MAX, 0)
    # the electron-beam drive (v0_2/beam.py) at the reference fill and in
    # pure Ar, its REBUILD variant, and the literature W it is scored on
    from v0_2 import beam as B
    y, y0 = B.yields(0.005), B.yields(0.0)
    yr = B.yields(0.005, "integral", "equal")
    n["beam.W"] = fx(y.W, 1)
    n["beam.exc"] = fx(y.ar_star, 2)
    n["beam.eps"] = fx(y.eps_sec, 1)
    n["beam.W0"] = fx(y0.W, 1)
    n["beam.ceil"] = fx(y0.one_photon_ceiling(N.HV_AR2), 3)
    n["beam.Wlit"] = fx(B.W_LIT_AR, 1)
    n["beam.Wmiss"] = fx(100 * (1 - y0.W / B.W_LIT_AR), 1)
    n["beam.rebW"] = fx(yr.W, 2)
    n["beam.rebexc"] = fx(yr.ar_star, 2)
    n["beam.ebar"] = fx(B.OPB_EBAR["Ar"], 1)
    tables = {"Rows": rate_rows(basic), "Group": rate_rows(group),
              "BandList": band_rows(), "Species": species_rows()}
    return n, tables


def rates_tex() -> str:
    n, tables = rates_parts()
    return render("LIGHT_ArF_ebeam reaction-rate tables", n, tables)


# ---------------------------------------------------------------------------
# the scores (from the recorded runs)
# ---------------------------------------------------------------------------
ALIAS = {"electra_formation": "elf", "electra_yield_band": "ely",
         "electra_gain_trend": "elt", "electra_gain_pump": "elp",
         "keio_formation": "keio", "sandia_formation": "san",
         "rokni_vs_ar": "rkp", "rokni_current_ratio": "rkc",
         "jancaitis_scan": "jan", "sri_vs_ar": "sria", "sri_vs_f2": "srif",
         "morozov_vuv": "moz", "kinetic_decays": "dec",
         "stern_volmer": "sv"}
MAXW = "_maxwellian"
MOROZOV_PRINTED = (0.33, 0.04)   # Morozov et al. 2008, 500-1400 hPa, printed


def alias(name) -> str:
    if name.endswith(MAXW):
        return ALIAS[name[:-len(MAXW)]] + "M"
    return ALIAS[name]


def _load():
    d0, d1 = MSR.load("score_0d.json"), MSR.load("score_1d.json")
    if d0 is None or d1 is None:
        raise SystemExit("arf_data/score_0d.json and score_1d.json are "
                         "needed: run arf_data/score_arf.py first")
    return d0["tests"], d1["tests"], d0, d1


def _pts(t, key="own"):
    return t.get(key) or []


def _stat(points, key="model"):
    r = MSR.ratios(points, key)
    return r, (MSR.gm(r) if r else float("nan")), \
        (MSR.rms_log(r) if r else float("nan"))


def _at(points, x, key="model", **match):
    """The point nearest x (with matching extra fields)."""
    c = [p for p in points if all(p.get(k) == v for k, v in match.items())]
    return min(c, key=lambda p: abs(p["x"] - x))


def score_parts():
    T0, T1, d0, d1 = _load()
    n, tables = {}, {}
    order = [k for k in MSR.TITLES if k in T0]
    for name in order:
        t, a = T0[name], alias(name)
        r, g, s = _stat(_pts(t))
        n[f"{a}.n"] = len(r)
        n[f"{a}.own"] = fx(g)
        n[f"{a}.ownrms"] = fx(s)
        if "beam" in t:
            _, g, s = _stat(t["beam"])
            n[f"{a}.beam"] = fx(g)
            n[f"{a}.beamrms"] = fx(s)
        if name in T1:
            _, g, s = _stat(T1[name], "centre")
            n[f"{a}.cone"] = fx(g)
            n[f"{a}.conerms"] = fx(s)
            _, g, _ = _stat(T1[name], "gap")
            n[f"{a}.gap"] = fx(g)

    # --- Electra, Keio, Sandia: the formation efficiency -----------------
    ef, eb = T0["electra_formation"]["own"], T0["electra_formation"]["beam"]
    for i, p in enumerate(ef, 1):
        b = eb[i - 1]
        n[f"elf.x{i}p"] = fx(p["x"], 2)
        n[f"elf.meas{i}"] = fx(p["measured"], 2)
        n[f"elf.mod{i}"] = fx(p["model"], 3)
        n[f"elf.r{i}"] = fx(p["model"] / p["measured"])
        n[f"elf.b{i}"] = fx(b["model"], 3)
        n[f"elf.br{i}"] = fx(b["model"] / p["measured"])
    for name, a in (("electra_formation", "elf"), ("keio_formation", "keio"),
                    ("sandia_formation", "san")):
        c1 = {round(q["x"], 2): q for q in T1.get(name, [])}
        for i, p in enumerate(T0[name]["own"], 1):
            q = c1.get(round(p["x"], 2))
            if q:
                n[f"{a}.c{i}"] = fx(q["centre"] / p["measured"])
    st = ef[0]["state"]
    n["elf.te"] = fx(st["te"], 1)
    n["elf.ne"] = sci(st["ne"], 2)
    n["elf.fleft"] = fx(100 * st["f2_left"], 0)
    y, yb = T0["electra_yield_band"]["own"][0], T0["electra_yield_band"]["beam"][0]
    n["ely.mod"] = fx(y["model"], 3)
    n["ely.r"] = fx(y["model"] / y["measured"])
    n["ely.b"] = fx(yb["model"], 3)
    n["ely.br"] = fx(yb["model"] / y["measured"])
    for name, a in (("keio_formation", "keio"), ("sandia_formation", "san")):
        own, beam = T0[name]["own"], T0[name]["beam"]
        for i, p in enumerate(own, 1):
            n[f"{a}.x{i}p"] = fx(p["x"], 2)
            n[f"{a}.meas{i}"] = fx(p["measured"], 2)
            n[f"{a}.mod{i}"] = fx(p["model"], 3)
            n[f"{a}.r{i}"] = fx(p["model"] / p["measured"])
            n[f"{a}.br{i}"] = fx(beam[i - 1]["model"] / p["measured"])
    # the Maxwellian closure at the same points
    for name, a in (("electra_formation", "elfM"), ("keio_formation", "keioM"),
                    ("sandia_formation", "sanM")):
        if name + MAXW in T0:
            for i, p in enumerate(T0[name + MAXW]["own"], 1):
                n[f"{a}.r{i}"] = fx(p["model"] / p["measured"])
                n[f"{a}.te{i}"] = fx(p["state"]["te"], 1)
    # Electra's single pass against pressure
    et, etb = T0["electra_gain_trend"]["own"], T0["electra_gain_trend"]["beam"]
    e1 = {round(p["x"], 2): p for p in T1.get("electra_gain_trend", [])}
    rows = []
    for p, b in zip(et, etb):
        q = e1.get(round(p["x"], 2))
        rows.append(f"{fx(p['x'], 2)} & {fx(p['g0_meas_pct_cm'], 2)} & "
                    f"{fx(p['measured'], 2)} & {fx(p['model'], 2)} & "
                    f"{fx(b['model'], 2)} & "
                    f"{fx(q['centre'], 2) if q else ''} & "
                    f"{fx(p['g0_calc_pct_cm'], 2)}\\\\")
    tables["ElectraTrend"] = rows
    hi = et[-1]
    n["elt.meas"] = fx(hi["measured"], 2)
    n["elt.mod"] = fx(hi["model"], 2)
    n["elt.hiatm"] = fx(hi["x"], 2)
    n["elt.loatm"] = fx(et[0]["x"], 2)
    n["elt.x"] = fx(hi["model"] / hi["measured"])
    n["elt.bx"] = fx(etb[-1]["model"] / hi["measured"])
    q = e1.get(round(hi["x"], 2))
    if q:
        n["elt.cx"] = fx(q["centre"] / hi["measured"])
    ep, epb = T0["electra_gain_pump"]["own"][0], T0["electra_gain_pump"]["beam"][0]
    n["elp.meas"] = fx(ep["measured"], 2)
    n["elp.mod"] = fx(ep["model"], 2)
    n["elp.b"] = fx(epb["model"], 2)
    n["elp.x"] = fx(ep["x"], 2)
    n["elp.r"] = fx(ep["model"] / ep["measured"])
    n["elp.br"] = fx(epb["model"] / ep["measured"])

    # --- Rokni ---------------------------------------------------------
    rk, rkb = T0["rokni_vs_ar"]["own"], T0["rokni_vs_ar"]["beam"]
    n["rkp.te"] = rng(min(p["te"] for p in rk), max(p["te"] for p in rk), 0)
    n["rkp.bte"] = rng(min(p["te"] for p in rkb), max(p["te"] for p in rkb), 1)
    one = {round(p["x"], 2): p for p in T1["rokni_vs_ar"]}
    own = {round(p["x"], 2): p for p in rk}
    beam0 = {round(p["x"], 2): p for p in rkb}
    for atm, key in ((2.0, "dep2"), (4.0, "dep4")):
        q = one[atm]
        n[f"rkp.{key}"] = fx(q["centre"] / own[atm]["model"])
        n[f"rkp.b{key}"] = fx(q["centre_beam"] / beam0[atm]["model"])
    n["rkp.te2"] = fx(own[2.0]["te"], 0)
    n["rkp.te4"] = fx(own[4.0]["te"], 0)
    S = importlib.import_module("score_arf")
    n["rkp.pw"] = fx(S.ROKNI_P_PER_ATM / 1e3, 0)
    n["rkp.xmin"] = fx(min(p["x"] for p in rk), 2)
    n["rkp.xmax"] = fx(max(p["x"] for p in rk), 0)
    if "rokni_vs_ar" + MAXW in T0:
        rm = T0["rokni_vs_ar" + MAXW]["own"]
        n["rkpM.te"] = rng(min(p["te"] for p in rm), max(p["te"] for p in rm), 1)
        mm = {round(p["x"], 2): p for p in rm}
        for atm, k in ((0.24, "lo"), (2.0, "two"), (4.0, "four")):
            n[f"rkp.p{k}"] = fx(atm, 2) if atm < 1 else fx(atm, 0)
            n[f"rkpM.r{k}"] = fx(mm[atm]["model"] / mm[atm]["measured"])
            n[f"rkp.r{k}"] = fx(own[atm]["model"] / own[atm]["measured"])
            b = {round(q["x"], 2): q for q in rkb}[atm]
            n[f"rkpB.r{k}"] = fx(b["model"] / b["measured"])
    rc, rcb = T0["rokni_current_ratio"]["own"], T0["rokni_current_ratio"]["beam"]
    n["rkc.te"] = rng(min(p["te"] for p in rc), max(p["te"] for p in rc), 0)
    n["rkc.telow"] = rng(min(p["te_low"] for p in rc),
                         max(p["te_low"] for p in rc), 0)
    n["rkc.meas"] = rng(min(p["measured"] for p in rc),
                        max(p["measured"] for p in rc), 0)
    rows = []
    rcm = {round(p["x"], 3): p for p in T0.get("rokni_current_ratio" + MAXW,
                                               {}).get("own", [])}
    for p, b in zip(rc, rcb):
        m = rcm.get(round(p["x"], 3))
        rows.append(f"{fx(p['x'], 2)} & {fx(p['measured'], 1)} & "
                    f"{fx(p['model'], 1)} & {fx(p['te'], 0)} & "
                    f"{fx(p['te_low'], 0)} & "
                    f"{fx(m['model'], 1) if m else ''} & "
                    f"{fx(m['te'], 1) if m else ''} & "
                    f"{fx(b['model'], 1)} & {fx(b['te'], 1)}\\\\")
    tables["RokniRatio"] = rows

    # --- Jancaitis -------------------------------------------------------
    jn, jb = T0["jancaitis_scan"]["own"], T0["jancaitis_scan"]["beam"]
    a2 = [p for p in jn if p["species"] == "Ar2F"]
    top2 = sorted(a2, key=lambda p: p["x"])[-2:]
    n["jan.axone"] = fx(top2[0]["x"], 0)
    n["jan.axtwo"] = fx(top2[1]["x"], 0)
    n["jan.aone"] = fx(top2[0]["model"] / top2[0]["measured"])
    n["jan.atwo"] = fx(top2[1]["model"] / top2[1]["measured"])
    for sp, k in (("ArF", "frms"), ("Ar2F", "arms")):
        n[f"jan.{k}"] = fx(_stat([p for p in jn if p["species"] == sp])[2])
    af = sorted((p for p in jn if p["species"] == "ArF"), key=lambda p: p["x"])
    n["jan.ftop"] = fx(af[-1]["model"] / af[-1]["measured"])
    n["jan.fxtop"] = fx(af[-1]["x"], 0)
    n["jan.jmin"] = fx(min(p["x"] for p in jn), 1)
    n["jan.jmax"] = fx(max(p["x"] for p in jn), 0)
    n["jan.te"] = rng(min(p["te"] for p in jn), max(p["te"] for p in jn), 1)
    n["jan.bte"] = rng(min(p["te"] for p in jb), max(p["te"] for p in jb), 1)
    if "jancaitis_scan" + MAXW in T0:
        jm = T0["jancaitis_scan" + MAXW]["own"]
        n["janM.te"] = rng(min(p["te"] for p in jm), max(p["te"] for p in jm), 1)

    # --- SRI ----------------------------------------------------------------
    sa = T0["sri_vs_ar"]["own"]
    for band, k in (("193", "a"), ("290", "b")):
        pts = [p for p in sa if str(p["band"]) == band]
        r, g, s = _stat(pts)
        n[f"sria.{k}n"] = len(r)
        n[f"sria.{k}gm"] = fx(g)
        n[f"sria.{k}rms"] = fx(s)
        n[f"sria.{k}rng"] = rng(min(r), max(r), 2)
    b290 = sorted((p for p in sa if str(p["band"]) == "290"),
                  key=lambda p: p["x"])
    n["sria.blo"] = fx(b290[0]["model"] / b290[0]["measured"])
    n["sria.bhi"] = fx(b290[-1]["model"] / b290[-1]["measured"])
    n["sria.bloatm"] = fx(b290[0]["x"], 0)
    n["sria.bhiatm"] = fx(b290[-1]["x"], 0)
    r, _, _ = _stat(T0["sri_vs_f2"]["own"])
    n["srif.rng"] = rng(min(r), max(r), 2)

    # --- Morozov ---------------------------------------------------------------
    mo, mb = T0["morozov_vuv"]["own"], T0["morozov_vuv"]["beam"]
    plat = [p for p in mb if p["x"] >= 500]
    n["moz.bplat"] = fx(np.mean([p["model"] for p in plat]), 2)
    n["moz.oplat"] = fx(np.mean([p["model"] for p in mo if p["x"] >= 500]), 2)
    n["moz.meas"] = fx(MOROZOV_PRINTED[0], 2)
    n["moz.err"] = fx(MOROZOV_PRINTED[1], 2)
    n["moz.frac"] = fx(MOROZOV_PRINTED[0]
                       / np.mean([p["model"] for p in plat]), 2)

    # --- consistency sets -------------------------------------------------------
    for name, tab in (("kinetic_decays", "Decays"), ("stern_volmer", "SV")):
        sets = {}
        for p in T0[name]["own"]:
            sets.setdefault(p["set"], []).append(p)
        rows = []
        for sname, pts in sets.items():
            r, g, s = _stat(pts)
            rows.append(f"{tex_text(sname)} & {len(r)} & {fx(g)} & "
                        f"{fx(s)}\\\\")
        tables[tab] = rows
    svr = [p for p in T0["stern_volmer"]["own"] if p["set"].startswith("Rokni")]
    n["sv.rokni"] = fx(_stat(svr)[1])

    # --- 1Dv0_2's centre against v0_2 at the device points, on both drives -----
    dev, gap, bdev, bgap = [], [], [], []
    for name in ("electra_formation", "keio_formation", "sandia_formation",
                 "electra_gain_trend"):
        own = {round(p["x"], 2): p for p in T0[name]["own"]}
        beam = {round(p["x"], 2): p for p in T0[name]["beam"]}
        for q in T1[name]:
            p = own[round(q["x"], 2)]
            dev.append((abs(q["centre"] / p["model"] - 1), name, q["x"]))
            b = beam[round(q["x"], 2)]
            bdev.append((abs(q["centre_beam"] / b["model"] - 1), name, q["x"]))
            if name != "electra_gain_trend":
                gap.append(q["gap"] / q["centre"])
                bgap.append(q["gap_beam"] / q["centre_beam"])
    test_name = {"sandia_formation": "Sandia", "keio_formation": "Keio"}
    worst = max(dev)
    n["od.dev"] = sci(worst[0], 1)
    n["od.devtest"] = test_name.get(worst[1], "Electra")
    n["od.gap"] = fx(100 * (1 - np.mean(gap)), 1)
    worst = max(bdev)
    n["od.bdev"] = sci(worst[0], 1)
    n["od.bdevtest"] = test_name.get(worst[1], "Electra")
    n["od.bgap"] = fx(100 * (1 - np.mean(bgap)), 1)
    n["od.cells"] = d1["cells"]

    # --- the Te regimes ----------------------------------------------------------
    def te_rng(points, key="te", nd=1):
        return rng(min(p[key] for p in points), max(p[key] for p in points), nd)
    rows = []
    rm = T0.get("rokni_vs_ar" + MAXW, {}).get("own")
    rcm_l = T0.get("rokni_current_ratio" + MAXW, {}).get("own")
    jm = T0.get("jancaitis_scan" + MAXW, {}).get("own")
    rows.append(r"Rokni Fig.~4, 20\unit{kW/cm^3} per atm & "
                f"{te_rng(rk)} & {te_rng(rm) if rm else ''} & "
                f"{te_rng(rkb)}\\\\")
    rows.append(r"Rokni Fig.~5, full current & "
                f"{te_rng(rc)} & {te_rng(rcm_l) if rcm_l else ''} & "
                f"{te_rng(rcb)}\\\\")
    rows.append(r"Rokni Fig.~5, 1/25 of the current & "
                f"{te_rng(rc, 'te_low')} & "
                f"{te_rng(rcm_l, 'te_low') if rcm_l else ''} & "
                f"{te_rng(rcb, 'te_low')}\\\\")
    rows.append(r"Jancaitis Fig.~3, 20 to 50\unit{ns} & "
                f"{te_rng(jn)} & {te_rng(jm) if jm else ''} & "
                f"{te_rng(jb)}\\\\")
    ref = T0["electra_formation"]["own"]
    rows.append(r"Electra, Keio, Sandia gain points & "
                + rng(min(p["state"]["te"] for t in ("electra_formation",
                                                    "keio_formation",
                                                    "sandia_formation")
                          for p in T0[t]["own"]),
                      max(p["state"]["te"] for t in ("electra_formation",
                                                    "keio_formation",
                                                    "sandia_formation")
                          for p in T0[t]["own"]), 1) + " & "
                + (rng(min(p["state"]["te"] for t in ("electra_formation",
                                                     "keio_formation",
                                                     "sandia_formation")
                           for p in T0[t + MAXW]["own"]),
                       max(p["state"]["te"] for t in ("electra_formation",
                                                     "keio_formation",
                                                     "sandia_formation")
                           for p in T0[t + MAXW]["own"]), 1)
                   if "electra_formation" + MAXW in T0 else "")
                + " & " + rng(min(p["state"]["te"] for t in (
                    "electra_formation", "keio_formation", "sandia_formation")
                    for p in T0[t]["beam"]), max(p["state"]["te"] for t in (
                        "electra_formation", "keio_formation",
                        "sandia_formation") for p in T0[t]["beam"]), 1)
                + "\\\\")
    del ref
    tables["TeRegimes"] = rows
    rows = []
    for f in sorted((ARF_DATA / "fig_readings").glob("*.csv")):
        head = [ln[1:].strip() for ln in f.read_text().splitlines()
                if ln.startswith("#")]
        method = head[1] if len(head) > 1 and " " in head[1] else ""
        rows.append(r"\texttt{"
                    + tex_name(f.stem).replace(r"\_", r"\_\allowbreak{}")
                    + "} & "
                    + tex_free(head[0]) + " & " + tex_free(method) + r"\\")
    tables["Readings"] = rows
    n["readings.n"] = len(rows)
    n["chem.md5"] = d0["chemistry_md5"]["network.py"][:8]
    return n, tables


def scores_tex() -> str:
    n, tables = score_parts()
    return render("LIGHT_ArF_ebeam accuracy numbers (the scores)", n, tables)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
COL = {"meas": "#000000", "own": "#0072B2", "maxw": "#CC79A7",
       "beam": "#E69F00", "cone": "#009E73", "y1": "#D55E00"}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
        "lines.markersize": 4, "legend.frameon": False,
        "legend.fontsize": 7, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False})
    return plt


def _save(fig, name):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight",
                metadata={"CreationDate": None, "ModDate": None,
                          "Creator": None, "Producer": None})
    import matplotlib.pyplot as plt
    plt.close(fig)


def figures_scores():
    plt = _plt()
    T0, T1, _, _ = _load()

    # --- the scorecard -------------------------------------------------------
    tests = [("electra_formation", "Electra formation efficiency"),
             ("keio_formation", "Keio formation efficiency"),
             ("sandia_formation", "Sandia formation efficiency"),
             ("electra_gain_trend", "Electra gain against pressure"),
             ("rokni_vs_ar", "Rokni 193 nm against Ar"),
             ("rokni_current_ratio", "Rokni current ratio"),
             ("jancaitis_scan", "Jancaitis against current"),
             ("sri_vs_ar", "SRI yields against Ar"),
             ("sri_vs_f2", "SRI 193-nm yield against F$_2$"),
             ("morozov_vuv", "Morozov pure-Ar VUV"),
             ("kinetic_decays", "decays (consistency)"),
             ("stern_volmer", "Stern-Volmer (consistency)")]
    fig, ax = plt.subplots(figsize=(6.2, 3.3))
    ax.axvspan(0.8, 1.25, color="0.92", lw=0)
    ax.axvline(1.0, color="0.5", lw=0.6)
    for i, (name, label) in enumerate(tests):
        y = len(tests) - i
        for key, src, dy, mk, col in (
                ("own", T0[name].get("own"), 0.2, "s", COL["own"]),
                ("maxw", T0.get(name + MAXW, {}).get("own"), 0.0, "o",
                 COL["maxw"]),
                ("beam", T0[name].get("beam"), -0.2, "^", COL["beam"])):
            if not src:
                continue
            r, g, s = _stat(src)
            if len(r) > 2:
                ax.plot([g * math.exp(-s), g * math.exp(s)], [y + dy] * 2,
                        color=col, lw=0.8)
            ax.plot([g], [y + dy], mk, color=col, ms=4)
        if name in T1:
            r, g, s = _stat(T1[name], "centre")
            ax.plot([g], [y - 0.38], "D", color=COL["cone"], ms=3)
    ax.set_yticks(range(len(tests), 0, -1))
    ax.set_yticklabels([lab for _, lab in tests])
    ax.set_xscale("log")
    ax.set_xlim(0.25, 4)
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.set_xlabel("model / measured, geometric mean (bar: $\\times e^{\\pm"
                  "\\mathrm{rms}}$)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="s", ls="", color=COL["own"],
                label="own input, two-temperature tail (the model)"),
         Line2D([], [], marker="o", ls="", color=COL["maxw"],
                label="own input, Maxwellian EEDF (sensitivity)"),
         Line2D([], [], marker="^", ls="", color=COL["beam"],
                label="e-beam drive"),
         Line2D([], [], marker="D", ls="", color=COL["cone"], ms=3,
                label="1Dv0_2 centre cell")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, "fig_scorecard")

    # --- formation efficiency ----------------------------------------------
    fig, ax = plt.subplots(figsize=(6.2, 2.5))
    cats = []
    for name, lab in (("electra_formation", "Electra"),
                      ("keio_formation", "Keio"),
                      ("sandia_formation", "Sandia")):
        for i, p in enumerate(T0[name]["own"]):
            b = T0[name]["beam"][i]
            m = T0.get(name + MAXW, {}).get("own")
            cats.append((f"{lab}\n{p['x']:.2f} MW/cm$^3$", p, b,
                         m[i] if m else None))
    for j, (lab, p, b, m) in enumerate(cats):
        ax.plot([j], [p["measured"]], "o", color=COL["meas"], ms=5)
        ax.plot([j - 0.12], [p["model"]], "s", color=COL["own"])
        if m:
            ax.plot([j], [m["model"]], "o", color=COL["maxw"], mfc="none")
        ax.plot([j + 0.12], [b["model"]], "^", color=COL["beam"])
    y = T0["electra_yield_band"]
    ax.plot([-0.12], [y["own"][0]["model"]], "s", color=COL["y1"], mfc="none")
    ax.plot([0.12], [y["beam"][0]["model"]], "^", color=COL["y1"], mfc="none")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c[0] for c in cats])
    ax.set_ylabel("$g_0 I_\\mathrm{sat}/P$")
    ax.set_ylim(0, None)
    h = [Line2D([], [], marker="o", ls="", color=COL["meas"], label="measured"),
         Line2D([], [], marker="s", ls="", color=COL["own"], label="own input"),
         Line2D([], [], marker="o", ls="", color=COL["maxw"], mfc="none",
                label="own input, Maxwellian"),
         Line2D([], [], marker="^", ls="", color=COL["beam"], label="beam drive"),
         Line2D([], [], marker="s", ls="", color=COL["y1"], mfc="none",
                label="ArF$^*$ yield 1.0 (Electra)")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, "fig_formation")

    # --- Rokni against Ar -------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.4))
    rk, rkb = T0["rokni_vs_ar"]["own"], T0["rokni_vs_ar"]["beam"]
    rm = T0.get("rokni_vs_ar" + MAXW, {}).get("own")
    x = [p["x"] for p in rk]
    a1.plot(x, [p["measured"] for p in rk], "o", color=COL["meas"],
            label="measured", ms=4)
    a1.plot(x, [p["model"] for p in rk], "s-", color=COL["own"],
            label="own input")
    if rm:
        a1.plot(x, [p["model"] for p in rm], "o--", color=COL["maxw"],
                mfc="none", label="own, Maxwellian")
    a1.plot(x, [p["model"] for p in rkb], "^-", color=COL["beam"],
            label="beam drive")
    c1 = sorted(T1["rokni_vs_ar"], key=lambda p: p["x"])
    a1.plot([p["x"] for p in c1], [p["centre"] for p in c1], "D",
            color=COL["cone"], ms=3, label="1Dv0_2 centre")
    a1.set_xlabel("Ar (atm)")
    a1.set_ylabel("193 nm, relative to 1 atm")
    a2.semilogy(x, [p["te"] for p in rk], "s-", color=COL["own"])
    if rm:
        a2.semilogy(x, [p["te"] for p in rm], "o--", color=COL["maxw"],
                    mfc="none")
    a2.semilogy(x, [p["te"] for p in rkb], "^-", color=COL["beam"])
    a2.axhline(40, color="0.5", lw=0.6, ls=":")
    a2.text(0.1, 43, "Fluid1D's 40 eV", fontsize=6.5, color="0.4")
    a2.set_xlabel("Ar (atm)")
    a2.set_ylabel("$T_e$ (eV), 250 to 300 ns")
    h, lab = a1.get_legend_handles_labels()
    a2.legend(h, lab, loc="center", fontsize=6.5)
    _save(fig, "fig_rokni")

    # --- Jancaitis -----------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.4))
    jn, jb = T0["jancaitis_scan"]["own"], T0["jancaitis_scan"]["beam"]
    jm = T0.get("jancaitis_scan" + MAXW, {}).get("own")
    for ax, sp, lab in ((axs[0], "ArF", "ArF$^*$ 193 nm"),
                        (axs[1], "Ar2F", "Ar$_2$F$^*$ 290 nm")):
        pts = sorted((p for p in jn if p["species"] == sp), key=lambda p: p["x"])
        bts = sorted((p for p in jb if p["species"] == sp), key=lambda p: p["x"])
        xs = [p["x"] for p in pts]
        ax.loglog(xs, [p["measured"] for p in pts], "o", color=COL["meas"],
                  label="measured", ms=4)
        ax.loglog(xs, [p["model"] for p in pts], "s-", color=COL["own"],
                  label="own input")
        if jm:
            mts = sorted((p for p in jm if p["species"] == sp),
                         key=lambda p: p["x"])
            ax.loglog([p["x"] for p in mts], [p["model"] for p in mts], "o--",
                      color=COL["maxw"], mfc="none", label="own, Maxwellian")
        ax.loglog([p["x"] for p in bts], [p["model"] for p in bts], "^-",
                  color=COL["beam"], label="beam drive")
        ax.set_title(lab, fontsize=8)
        ax.set_xlabel("current density (A/cm$^2$)")
    axs[0].set_ylabel("peak, relative to 74 A/cm$^2$")
    axs[0].legend(loc="upper left")
    _save(fig, "fig_jancaitis")

    # --- SRI -----------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.4))
    sa = T0["sri_vs_ar"]["own"]
    for band, mk, ls in (("193", "o", "-"), ("290", "v", "--")):
        pts = sorted((p for p in sa if str(p["band"]) == band),
                     key=lambda p: p["x"])
        a1.plot([p["x"] for p in pts], [p["measured"] for p in pts], mk,
                color=COL["meas"], mfc="none" if band == "290" else None,
                label=f"{band} nm measured")
        xs = sorted({p["x"] for p in pts})
        mod = [next(p["model"] for p in pts if p["x"] == x) for x in xs]
        a1.plot(xs, mod, ls, color=COL["own"], label=f"{band} nm model")
    a1.set_xlabel("Ar (atm), 1 Torr F$_2$")
    a1.set_ylabel("photons per initial excitation")
    a1.legend(loc="upper right", fontsize=6.5)
    sf = T0["sri_vs_f2"]["own"]
    for atm, mk in ((1.0, "o"), (2.0, "s"), (5.0, "^")):
        pts = sorted((p for p in sf if abs(p["ar_atm"] - atm) < 1e-6),
                     key=lambda p: p["x"])
        a2.plot([p["x"] for p in pts], [p["measured"] for p in pts], mk,
                color=COL["meas"], mfc="none", label=f"{atm:g} atm measured")
        a2.plot([p["x"] for p in pts], [p["model"] for p in pts], "-",
                color=COL["own"])
    a2.set_yscale("log")
    a2.set_xlabel("F$_2$ (Torr)")
    a2.set_ylabel("193-nm photons per excitation")
    a2.legend(loc="upper right", fontsize=6.5)
    _save(fig, "fig_sri")

    # --- Morozov --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    mo, mb = T0["morozov_vuv"]["own"], T0["morozov_vuv"]["beam"]
    ax.axhspan(MOROZOV_PRINTED[0] - MOROZOV_PRINTED[1],
               MOROZOV_PRINTED[0] + MOROZOV_PRINTED[1], xmin=0.3, color="0.9",
               lw=0)
    ax.plot([p["x"] for p in mo], [p["measured"] for p in mo], "o",
            color=COL["meas"], label="measured", ms=4)
    ax.plot([p["x"] for p in mo], [p["model"] for p in mo], "s-",
            color=COL["own"], label="own input")
    ax.plot([p["x"] for p in mb], [p["model"] for p in mb], "^-",
            color=COL["beam"], label="beam drive")
    ax.set_ylim(0, 0.8)
    ax.set_xlabel("pure Ar (hPa)")
    ax.set_ylabel("VUV power / deposited power")
    ax.legend(loc="lower right")
    _save(fig, "fig_morozov")


# ---------------------------------------------------------------------------
# the runs
# ---------------------------------------------------------------------------
def _score_arf():
    import score_arf as S
    return S


def run_parts(tmp, quick=False):
    """(numbers, tables, trace) from fresh runs. quick: one band and coarse
    meshes, a smoke test of this code only (its numbers are not the
    overview's)."""
    from v0_2 import report as R
    from v0_2 import run as RUN
    from v0_2.model import Model, Params, run
    S = _score_arf()
    n, tables = {}, {}

    # --- the reference point ---------------------------------------------------
    p0 = Params()
    res = run(p0)
    a = R.accounting(res)
    d = R.row_dict(res)
    e = a["e_dep"]
    n["ref.p"] = fx(p0.pressure_atm, 0)
    n["ref.f2"] = fx(100 * p0.x_f2, 1)
    n["ref.kw"] = fx(p0.power_W_cm3 / 1e3, 0)
    n["ref.pulse"] = fx(p0.pulse_s * 1e6, 0)
    n["ref.tend"] = fx(p0.t_end_s * 1e6, 0)
    n["ref.e193"] = fx(d["emission_193"])
    n["ref.e193x"] = fx(d["emission_193"], 4)
    n["ref.e127"] = fx(d["emission_vuv"])
    n["ref.e127x"] = fx(d["emission_vuv"], 4)
    n["ref.e284"] = fx(100 * a["bands"].get(N.TRIMER_BAND, 0.0) / e)
    n["ref.formed"] = fx(d["to_ArFB"], 1)
    n["ref.toars"] = fx(d["to_Arstar"], 1)
    n["ref.toarp"] = fx(d["to_Arplus"], 2)
    n["ref.tg"] = fx(d["tg_end_K"], 0)
    n["ref.tgx"] = fx(d["tg_end_K"], 2)
    n["ref.te"] = fx(d["te_mean_eV"], 1)
    n["ref.temax"] = fx(float(np.max(res.te[res.t <= p0.pulse_s])), 1)
    n["ref.closure"] = sci(abs(a["closure"]), 1)
    burn = R.f2_burnout_time(res)
    n["ref.burn"] = fx(burn * 1e9, 0) if burn else "none"
    rows = [f"193\\unit{{nm}} (ArF B--X) & {fx(d['emission_193'])}\\\\",
            f"284\\unit{{nm}} (\\AtF{{}}) & "
            f"{fx(100 * a['bands'].get(N.TRIMER_BAND, 0.0) / e)}\\\\",
            f"127\\unit{{nm}} (\\Arts{{}}) & {fx(d['emission_vuv'])}\\\\",
            f"heat, elastic & {fx(100 * a['heat_elastic'] / e)}\\\\",
            f"heat, reactions and captured electrons & "
            f"{fx(100 * a['heat_reactions'] / e)}\\\\"]
    for label, v in a["stores"].items():
        rows.append(f"stored, {tex_text(label)} & {fx(100 * v / e)}\\\\")
    tables["Budget"] = rows
    arf = a["arf"]
    ev = a["events"]
    formed_by = {}
    for c in res.model.compiled:
        r = c.row
        k = r.products.count("ArF(B)") - r.reactants.count("ArF(B)")
        if k > 0:
            formed_by[r.rid] = k * ev[r.rid]
    f = arf["formed"]
    n["ref.chharp"] = fx(100 * formed_by.get("G3", 0.0) / f, 0)
    n["ref.chdimer"] = fx(100 * (formed_by.get("NK7", 0.0)
                                 + formed_by.get("K2F_F", 0.0)) / f, 0)
    n["ref.chion"] = fx(100 * (formed_by.get("K26", 0.0)
                               + formed_by.get("K108", 0.0)) / f, 0)
    for rid, k in (("G4", "rad"), ("G5", "tri"), ("QAR", "qar"),
                   ("QF2", "qf"), ("QE", "qe")):
        n[f"ref.arf.{k}"] = fx(100 * arf["losses"].get(rid, 0.0) / f, 0)
    ars = a["arstar"]
    for rid, k in (("G3", "harp"), ("G3X", "rest"), ("AK8_Ar", "dimer")):
        n[f"ref.ars.{k}"] = fx(100 * ars["losses"].get(rid, 0.0)
                               / ars["formed"], 0)

    # --- the Maxwellian and the seeds ----------------------------------------------
    dm = R.row_dict(run(p0.replace(eedf="maxwellian")))
    n["ref.m193"] = fx(dm["emission_193"])
    n["ref.m127"] = fx(dm["emission_vuv"])
    n["ref.mformed"] = fx(dm["to_ArFB"], 1)
    n["ref.mte"] = fx(dm["te_mean_eV"], 1)
    n["ref.mtg"] = fx(dm["tg_end_K"], 0)
    n["ref.seedlo"] = sci(1e10)
    n["ref.seedhi"] = sci(1e15)
    seeds = {1e10: None, 1e13: None, 1e15: d["emission_193"]}
    for s in (1e10, 1e13):
        seeds[s] = R.row_dict(run(p0.replace(ne_cm3=s)))["emission_193"]
    n["ref.seeds"] = rng(min(seeds.values()), max(seeds.values()), 2)

    # --- the bands at both ends -------------------------------------------------------
    saved = RUN.BANDS
    if quick:
        RUN.BANDS = {"ar_f2": N.BANDS["ar_f2"]}
    try:
        brs = RUN.do_bands(p0)
    finally:
        RUN.BANDS = saved
    base = brs[0]["emission_193"]
    rows, moves = [], {}
    for b in N.BANDS.values():
        if not any(q["band"] == b.name for q in brs):
            continue
        lo = next(q for q in brs if q["band"] == b.name and q["end"] == "low")
        hi = next(q for q in brs if q["band"] == b.name and q["end"] == "high")
        f_ = (lambda v: fx(v, 2)) if b.unit in ("fraction", "exponent") \
            else (lambda v: sci(v, 2))
        rows.append(f"{tex_name(b.name)} & {f_(b.lo)} & "
                    f"{fx(lo['emission_193'])} & {fx(lo['emission_vuv'], 1)} & "
                    f"{f_(b.hi)} & {fx(hi['emission_193'])} & "
                    f"{fx(hi['emission_vuv'], 1)}\\\\")
        moves[b.name] = max(abs(lo["emission_193"] - base),
                            abs(hi["emission_193"] - base))
        n[f"band.{b.name.replace('_', '')}.lo"] = fx(lo["emission_193"])
        n[f"band.{b.name.replace('_', '')}.hi"] = fx(hi["emission_193"])
    tables["BandEnds"] = rows
    levers = ("ar2_f", "arf_yield", "arf_trimer", "trimer_tg_pow", "ar2_form")
    n["band.othermax"] = fx(max([v for k, v in moves.items()
                                 if k not in levers] or [0.0]))

    # --- Morozov's point across the Ar* and Ar2* bands ---------------------------------
    pa = 1000.0 / 1013.25

    def moz(beam, rates=()):
        p = S.params(pa, 0.0, 330.0, 450e-9, 450e-9 + 20e-6, ne_cm3=1e8,
                     n_save=2001, rates=rates)
        m, r = S.run0(p, beam)
        return float(r.events["AR2_RAD"][-1]) * N.HV_AR2 \
            / float(r.ledger["E_dep"][-1])
    cases = {"default": (), "ar_pen": (("ar_pen", N.BANDS["ar_pen"].hi),),
             "ar2_pen": (("ar2_pen", N.BANDS["ar2_pen"].hi),),
             "ar2_form lo": (("ar2_form", N.BANDS["ar2_form"].lo),),
             "ar2_form hi": (("ar2_form", N.BANDS["ar2_form"].hi),),
             "ar2_rad": (("ar2_rad", N.BANDS["ar2_rad"].hi),)}
    vals = {k: moz(True, v) for k, v in cases.items()}
    n["mozb.def"] = fx(vals["default"], 3)
    n["mozb.min"] = fx(min(vals.values()), 3)
    n["mozb.max"] = fx(max(vals.values()), 3)
    n["mozb.minband"] = tex_name(min(vals, key=vals.get).split()[0])

    # --- R11p held at 30 eV at Rokni's pressures -------------------------------------------
    from v0_2 import xsec

    class ClipR11p(Model):
        def _k(self, c, te, tg, nn, ne):
            if c.kind == "te_pow":
                te = min(te, xsec.TE_MAX)
            return super()._k(c, te, tg, nn, ne)

    def rokni_em(cls, pa_):
        tot = pa_ + 2 / 760
        p = S.params(tot, (2 / 760) / tot, S.ROKNI_P_PER_ATM * pa_, 300e-9)
        m = cls(p)
        r = m.run()
        sel = (r.t >= 250e-9) & (r.t <= 300e-9)
        return float(np.mean(S.emission(m, r, "G4")[sel]))
    dev = max(abs(rokni_em(ClipR11p, x) / rokni_em(Model, x) - 1)
              for x in (0.243, 0.494, 1.0, 2.0, 4.003))
    n["rclip.dev"] = fx(100 * dev, 2)

    # --- 1Dv0_2: the reference point on 64 and 128 cells ----------------------------------
    sim1d = importlib.import_module("1Dv0_2.sim1d")
    rows = []
    mesh = {}
    meshes = (16, 32) if quick else (64, 128)
    for cells in meshes:
        p1 = sim1d.Params1D(cells=cells, verbose=False, name=f"ref{cells}",
                            output_dir=str(tmp), n_profiles=3)
        r1 = sim1d.run(p1)
        last = r1.data["records"][-1]
        led, ledc = last["ledger"], last["ledger_centre"]
        vuv = f"rad:{N.VUV_BAND}"
        fate = r1.fate("ledger")
        mesh[cells] = {
            "c193": r1.share_193("ledger_centre"), "g193": r1.share_193(),
            "c127": 100 * ledc[vuv] / ledc["input"],
            "g127": 100 * led[vuv] / led["input"],
            "ctg": last["tg_centre"], "gtg": last["tg_mean"],
            "cut": fate["cut by the mean-energy bounds"],
            "walls": fate["lost to the walls"],
            "closure": fate["closure residual"],
            "phantom": -100 * r1.data["field_dropped"] / led["input"]}
    for key, lab, nd, zero in (
            ("193", r"193\unit{nm}, \% of the input", 4, d["emission_193"]),
            ("127", r"127\unit{nm}, \% of the input", 4, d["emission_vuv"]),
            ("tg", r"$T_g$ at the end, K", 1, d["tg_end_K"])):
        c0, c1_ = meshes
        rows.append(f"{lab} & {fx(zero, nd)} & "
                    f"{fx(mesh[c0]['c' + key], nd)} & "
                    f"{fx(mesh[c0]['g' + key], nd)} & "
                    f"{fx(mesh[c1_]['c' + key], nd)} & "
                    f"{fx(mesh[c1_]['g' + key], nd)}\\\\")
    tables["Mesh"] = rows
    for cells, tag in zip(meshes, ("a", "b")):
        m_ = mesh[cells]
        n[f"m{tag}.cells"] = cells
        n[f"m{tag}.cut"] = fx(m_["cut"])
        n[f"m{tag}.walls"] = fx(m_["walls"], 3)
        n[f"m{tag}.phantom"] = fx(m_["phantom"], 1)
        n[f"m{tag}.closure"] = sci(abs(m_["closure"]), 1) \
            if m_["closure"] else "0"
        n[f"m{tag}.centredev"] = sci(abs(m_["c193"] / d["emission_193"] - 1),
                                     1)
        n[f"m{tag}.g193"] = fx(m_["g193"], 3)

    tn, trace = trace_parts(tmp, quick)
    n.update(tn)
    return n, tables, trace


def trace_parts(tmp, quick=False):
    """1Dv0_2 at Rokni's 2 atm against v0_2: (numbers, trace)."""
    S = _score_arf()
    sim1d = importlib.import_module("1Dv0_2.sim1d")
    n = {}
    pa_ = 2.0
    tot = pa_ + 2 / 760
    x_f2 = (2 / 760) / tot
    pw = S.ROKNI_P_PER_ATM * pa_
    p1 = sim1d.Params1D(pressure_atm=tot, x_f2=x_f2, power_W_cm3=pw,
                        pulse_s=300e-9, t_end_s=300e-9, ne_cm3=1e14,
                        cells=16 if quick else 64,
                        dt_s=50e-12, n_save=401, n_profiles=3, verbose=False,
                        name="rokni2", output_dir=str(tmp))
    r1 = sim1d.run(p1)
    m0, r0 = S.run0(S.params(tot, x_f2, pw, 300e-9))
    t1 = r1.t
    te1 = np.array([q["te_centre"] for q in r1.data["records"]])
    ceiling = 2.0 / 3.0 * r1.data["eps_bounds_eV"][1]
    hit = np.nonzero(te1 >= 0.995 * ceiling)[0]
    t_hit = float(t1[hit[0]]) if hit.size else float("nan")
    # the two Te part (by PART) when the 1-D one meets the ceiling; until
    # then every species agrees but Ar+, whose production (the tail rows'
    # ionization) is so steep in Te that it moves first
    PART = 0.01
    te0 = np.interp(t1, r0.t, r0.te)
    t_part = float(t1[np.nonzero((t1 > 50e-9)
                                 & (np.abs(te1 / te0 - 1) > PART))[0][0]])
    track = ("e-", "Ar*", "Ar2*", "Ar2+", "F-", "ArF(B)", "Ar2F*")
    sel = (t1 >= 100e-9) & (t1 <= t_part)
    worst = 0.0
    for sp in track:
        c1 = r1.series(sp, "centre")[sel]
        c0 = np.interp(t1[sel], r0.t, r0.ne if sp == "e-" else r0.n[sp])
        worst = max(worst, float(np.max(np.abs(c1 / c0 - 1))))
    c1 = r1.series("Ar+", "centre")[sel]
    arp = float(np.max(np.abs(c1 / np.interp(t1[sel], r0.t, r0.n["Ar+"]) - 1)))
    n["tr.tpart"] = fx(t_part * 1e9, 0)
    n["tr.partpct"] = fx(100 * PART, 0)
    n["tr.arpagree"] = fx(100 * arp, 0)
    n["tr.ceil"] = fx(ceiling, 0)
    n["tr.eps"] = fx(r1.data["eps_bounds_eV"][1], 0)
    n["tr.from"] = fx(100, 0)
    n["tr.thit"] = fx(t_hit * 1e9, 0)
    n["tr.agree"] = fx(100 * worst, 1)
    n["tr.te0"] = fx(float(r0.te[-1]), 0)
    n["tr.arp"] = fx(float(r0.n["Ar+"][-1]) / float(r1.series("Ar+", "centre")[-1]), 1)
    n["tr.arf"] = fx(float(r0.n["ArF(B)"][-1])
                     / float(r1.series("ArF(B)", "centre")[-1]), 2)
    # the ionization rate at the 0-D state at 300 ns, Te 30, 40 and the 0-D's
    j = len(r0.t) - 1
    dens = np.array([r0.n[s][j] for s in m0.species])
    ne = float(r0.ne[j])
    ng = float(dens.sum())
    n["tr.nen"] = sci(ne / ng, 1)
    n["tr.tail"] = fx(float(m0.tail(ne / ng, r0.n["F2"][j] / ng)), 3)
    c2 = next(c for c in m0.compiled if c.row.rid == "C2")
    for te, k in ((30.0, "thirty"), (40.0, "forty"), (float(r0.te[j]), "end")):
        n[f"tr.ctwo.{k}"] = sci(m0._k(c2, te, float(r0.tg[j]), dens, ne), 2)
    trace = {"t1": t1, "te1": te1, "t0": r0.t, "te0": r0.te,
             "arf1": r1.series("ArF(B)", "centre"), "arf0": r0.n["ArF(B)"],
             "arp1": r1.series("Ar+", "centre"), "arp0": r0.n["Ar+"],
             "ceiling": ceiling}
    return n, trace


def figure_trace(tr):
    plt = _plt()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.4))
    a1.plot(tr["t0"] * 1e9, tr["te0"], "-", color=COL["own"], label="v0_2")
    a1.plot(tr["t1"] * 1e9, tr["te1"], "--", color=COL["cone"],
            label="1Dv0_2 centre")
    a1.axhline(tr["ceiling"], color="0.5", lw=0.6, ls=":")
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("$T_e$ (eV)")
    a1.legend(loc="lower right")
    a2.semilogy(tr["t0"] * 1e9, tr["arf0"], "-", color=COL["own"],
                label="ArF(B), v0_2")
    a2.semilogy(tr["t1"] * 1e9, tr["arf1"], "--", color=COL["cone"],
                label="ArF(B), 1Dv0_2")
    a2.semilogy(tr["t0"] * 1e9, tr["arp0"], "-", color=COL["beam"],
                label="Ar$^+$, v0_2")
    a2.semilogy(tr["t1"] * 1e9, tr["arp1"], "--", color=COL["y1"],
                label="Ar$^+$, 1Dv0_2")
    a2.set_xlabel("time (ns)")
    a2.set_ylabel("density (cm$^{-3}$)")
    a2.set_ylim(1e8, None)
    a2.legend(loc="center right", fontsize=6.5)
    _save(fig, "fig_trace")


def runs_tex(write_figure=True) -> str:
    tmp = Path(tempfile.mkdtemp(prefix="light_arf_overview_"))
    try:
        n, tables, trace = run_parts(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if write_figure:
        figure_trace(trace)
    return render("LIGHT_ArF_ebeam accuracy numbers (the runs)", n, tables)


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--runs", action="store_true",
                    help="also make the runs (~20 min, needs Fluid1D)")
    ap.add_argument("--check", choices=("rates", "scores", "runs"),
                    help="compare one file with a fresh one, write nothing")
    ap.add_argument("--trace-figure", action="store_true",
                    help="redraw figures/fig_trace only (runs the trace)")
    ap.add_argument("--quick-runs", action="store_true",
                    help="smoke-test the runs code (one band, coarse "
                         "meshes), print its numbers, write nothing")
    a = ap.parse_args(argv)
    if a.trace_figure:
        tmp = Path(tempfile.mkdtemp(prefix="light_arf_overview_"))
        try:
            _, trace = trace_parts(tmp)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        figure_trace(trace)
        print(f"wrote {FIGS}/fig_trace")
        return 0
    if a.quick_runs:
        tmp = Path(tempfile.mkdtemp(prefix="light_arf_overview_"))
        try:
            n, tables, _ = run_parts(tmp, quick=True)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        for k, v in n.items():
            print(f"{k:22s} {v}")
        for k, rows in tables.items():
            print(k, len(rows), rows[0][:150])
        return 0
    if a.check:
        fresh = {"rates": rates_tex, "scores": scores_tex,
                 "runs": lambda: runs_tex(write_figure=False)}[a.check]()
        path = {"rates": OUT_RATES, "scores": OUT_SCORES,
                "runs": OUT_RUNS}[a.check]
        if not path.is_file() or path.read_text() != fresh:
            print(f"{path.name} is stale: rerun make_overview.py")
            return 1
        print(f"{path.name} is current")
        return 0
    OUT_RATES.write_text(rates_tex())
    print(f"wrote {OUT_RATES}")
    OUT_SCORES.write_text(scores_tex())
    print(f"wrote {OUT_SCORES}")
    figures_scores()
    print(f"wrote {FIGS}/fig_scorecard, fig_formation, fig_rokni, "
          "fig_jancaitis, fig_sri, fig_morozov")
    if a.runs:
        OUT_RUNS.write_text(runs_tex())
        print(f"wrote {OUT_RUNS} and {FIGS}/fig_trace")
    return 0


if __name__ == "__main__":
    sys.exit(main())
