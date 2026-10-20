"""Numbers and tables for accuracy_overview_v0_1.tex, v0_1's accuracy.

The v0 overview's generator (make_overview.py) run on v0_1, with the checks
v0_1 needs besides. No number in the overview is typed by hand. This script
writes three files, each one \\videf{key}{value} per number or a finished
table (a macro of its own, so one document can hold v0's numbers and
v0_1's):

  overview_v0_1_rates.tex  Appendix A, the reaction rates of v0_1's rows
  overview_v0_1_fast.tex   the checks that take seconds: the photolysis
                           decays, the ledger closure, ML's scorecard with
                           Tahoma (a score, never a fit target), Tahoma
                           with v0_1's default network, with the trimer's
                           weak temperature law and with v0_1's own power
                           input, and the pure-gas checks
  overview_v0_1_slow.tex   the checks that take minutes: the posterior's
                           Tahoma gain, each unmeasured formation link, the
                           electron closure against the kinetic EEDF loop,
                           and the cross sections against swarm data

The accuracy numbers come from ML/, which runs v0_1 when LIGHT_VERSION
names it. ML fixes its version when first imported, so the numbers are
made in a process of their own: main() sets LIGHT_VERSION=v0_1 before ML
is imported, and every ML-backed function refuses to run on another
version. test_docs_v0_1.py compares the files through --check.

    python3 LIGHT/docs/make_overview_v0_1.py               # all, ~6 min
    python3 LIGHT/docs/make_overview_v0_1.py --fast        # rates + fast
    python3 LIGHT/docs/make_overview_v0_1.py --check fast  # compare only
    cd LIGHT/docs
    latexmk -pdf system_map_v0_1.tex
    latexmk -pdf accuracy_overview_v0_1.tex
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_1.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_1.py").is_file())
ML_DIR = LIGHT.parent / "ML"
for _p in (str(LIGHT), str(HERE), str(ML_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                            # noqa: E402

import make_overview as MO                                    # noqa: E402
from make_system_map_basic import _k                          # noqa: E402
from v0_1 import __version__, xsec                            # noqa: E402
from v0_1.model import AMAGAT_CM3, Model, Params, run         # noqa: E402
from v0_1.network import BASIC_GROUPS, GROUPS, build          # noqa: E402
from v0_1.report import accounting                            # noqa: E402
from v0_1.species import SPECIES, e_form                      # noqa: E402

VERSION = "v0_1"
OUT_RATES = HERE / "overview_v0_1_rates.tex"
OUT_FAST = HERE / "overview_v0_1_fast.tex"
OUT_SLOW = HERE / "overview_v0_1_slow.tex"

# ---------------------------------------------------------------------------
# Appendix A: the reaction rates
# ---------------------------------------------------------------------------
# the order the tables list the rows in (v0's, with v0_1's new rows where
# they belong). A row missing here is appended at the end, never dropped.
ORDER = ("G1", "C1a", "C1b", "C1c", "C2", "G6", "ATT", "DIS_a3", "DIS_A1",
         "VIB1", "VIB2", "VIB3", "VIB4", "G2", "G3", "G3C", "G3X",
         "K15", "R11p", "K26", "K108",
         "G4", "G4A", "G4C", "BC", "CB", "G5", "G5C", "QE", "QEC", "QF2",
         "QF2C", "QKR", "QKRC", "C3", "Q2F", "FF",
         "AK8_Kr", "KR2_RAD", "NK7", "K2F_BR", "K2F_F", "PEN_KR2", "KSTEP",
         "SUPER", "K33", "K32", "K67", "K125", "K_2Fe", "BC_Kr", "CB_Kr",
         "URAD_C")


def ordered(rows) -> list:
    """Rows in ORDER. A C-state copy (bc group, 'G5_C') follows its B row."""
    rank = {rid: i for i, rid in enumerate(ORDER)}

    def key(r):
        if r.rid not in rank and r.rid.endswith("_C"):
            return (rank.get(r.rid[:-2], len(ORDER)), 1, r.rid)
        return (rank.get(r.rid, len(ORDER)), 0, r.rid)
    return sorted(rows, key=key)


# Short sources, one per row, condensed from each row's full citation in
# v0_1/network.py. Every number here must appear in that citation
# (test_docs_v0_1.py checks it), so the table cannot drift from the code.
SOURCES = {
    **{k: v for k, v in MO.SOURCES.items()},
    "G4": "Kannari 1985 Table I, the 96.6 % B-X share of tau_B = 6.99 ns. "
          "Eden 1978 measured 6.8 +- 0.2 ns",
    "G4A": "A(B-A)/A(B-X) = 0.035, Dunning and Hay's calculation via Yu "
           "and Setser 1990. 272 nm is a theoretical position",
    "G4C": "ML joint fit. Hakuta 1990 measured 60 +- 10 ns. Band at 275 "
           "nm, printed",
    "C3": "Xu, Gadomski and Setser 1993, 220 +- 20 ns. Band at 420 nm, "
          "printed",
    "CB": "detailed balance of BC, with the joint fit's K_eq and Yu and "
          "Setser's measured 120 +- 50 cm-1 gap",
    "G5": "ML joint fit. T_g law from Shui's calculation for 300-800 K, "
          "measured only to 425 K",
    "K15": "the mean of three measurements, 2.4e-31. T_g law from Perez-Rios "
           "and Greene 2018",
    "R11p": "Shiu and Biondi 1977, 1.6e-6 (300 K/T_e)^0.55, +- 10 %. Mikus "
            "1978 within 8 %",
    "K26": "Kannari 1985, Flannery turnover with Debye screening, and "
           "(T_g/300)^-3 on its three-body part from Thomson theory. "
           "Measured with F- only below 300 Torr, Jacob et al. 1978",
    "K108": "as K26",
    "FF": "Cohen and Bott 1976 evaluation, 2.8e-34, every rare gas given "
          "the Ar value. Held at its 300 K value",
    "PEN_KR2": "Kannari 1985 Table I, from an unpublished 1973 SRI report "
               "that gives Ar2* the same value. No measurement held",
}
# v0_1's v0 network (the notes' list as written, REBUILD's copies) where it
# differs from the basic network
SOURCES_V0 = {
    **{k: v for k, v in MO.SOURCES_V0.items()},
    "G4": "Kannari 1985 Table I, all of tau_B as B-X",
    "K15": "Kannari 1985, Boichenko 2000",
    "R11p": "Kannari 1985 Table I, 1.39e-6 at 300 K on a T_e^-1/2 template",
    "K26": "Kannari 1985, Flannery turnover with Debye screening, no T_g "
           "law",
    "K108": "as K26",
}
SPECIES_NOTES = {
    **{k: v for k, v in MO.SPECIES_NOTES.items()},
    "KrF(B)": "D_e(B) = 5.54 eV below the ion pair, Tellinghuisen 1976",
    "KrF(C)": "120 cm-1 below B, Yu and Setser 1990",
    "Kr2F*": "0.58 eV below KrF* + Kr, calculated by Huestis et al. 1984",
}

# LaTeX for the laws the short sources print (before the generic rules,
# which would split their T_g and T_e)
_LAWS = {
    "(300 K/T_e)^0.55": r"$(300\,\mathrm{K}/T_e)^{0.55}$",
    "(T_g/300)^-3": r"$(T_g/300)^{-3}$",
    "tau_B": r"$\tau_B$",
    "A(B-A)/A(B-X)": r"$A_{B\text{-}A}/A_{B\text{-}X}$",
    "D_e(B)": r"$D_e(B)$",
    "Ar2*": r"Ar$_2^*$",
}


def _tex_text(s: str) -> str:
    """A short source as LaTeX text: v0's rules (make_overview._tex_text)
    after v0_1's laws are set as math."""
    keep = {}
    for i, (plain, tex) in enumerate(_LAWS.items()):
        if plain in s:
            mark = f"@@{i}@@"
            s = s.replace(plain, mark)
            keep[mark] = tex
    # a bare % would comment out the rest of the table row, its \\ too
    s = re.sub(r" ?%", r"\\,\\%", s)
    s = MO._tex_text(s)
    for mark, tex in keep.items():
        s = s.replace(mark, tex)
    return s


def rate_text(r, p: Params) -> str:
    """The rate coefficient at 300 K as printed in the table."""
    kind = r.rate[0]
    if kind == "inv_tau":
        frac = r.rate[2] if len(r.rate) > 2 else 1.0
        return f"${_k(frac / getattr(p, r.rate[1]))}$ {MO._unit(r)}"
    if kind == "te_pow":
        k0, te0, n_pow = r.rate[1:]
        return (f"${_k(k0)}$ $(T_e/{te0 * 1.160451812e4:.0f}"
                rf"\,\mathrm{{K}})^{{{n_pow:g}}}$ {MO._unit(r)}")
    if kind == "flannery":
        n_am = 101325.0 / (1.380649e-23 * 300.0) * 1e-6 / AMAGAT_CM3
        k = r.rate[1] * n_am / (n_am * n_am + 1.7) * r.rate[2]
        return f"${_k(k)}$ {MO._unit(r)} at 1\\,atm, 300\\,K, unscreened"
    return MO.rate_text(r, p)


def tdep_text(r) -> str:
    kind = r.rate[0]
    if kind == "te_pow":
        parts = [f"$T_e^{{{r.rate[3]:g}}}$"]
        if r.tg_pow:
            parts.append(f"$(T_g/300)^{{{r.tg_pow:.2f}}}$")
        return ", ".join(parts)
    if kind == "flannery" and len(r.rate) > 3:
        return (f"$N$, screening, $(T_g/300)^{{{r.rate[3]:g}}}$ on the "
                "3-body part")
    return MO.tdep_text(r)


def _rate_rows(rows, sources, p) -> str:
    out = []
    for r in rows:
        src = sources.get(r.rid)
        if src is None and r.rid.endswith("_C"):
            src = f"as {r.rid[:-2]}, applied to KrF(C)"
        if src is None:
            raise KeyError(f"no short source for row {r.rid}")
        rid = r.rid.replace("_", "\\_")
        out.append(f"{rid} & {MO.reaction(r)} & {rate_text(r, p)} & "
                   f"{tdep_text(r)} & {_tex_text(src)}\\\\")
    return "\n".join(out)


def v0_differences() -> list:
    """v0_1's v0 network rows (every group) that its basic network lacks or
    runs with another value."""
    basic = {r.rid: r for r in build(("all",), "basic")}
    out = []
    for r in build(("all",), "v0"):
        b = basic.get(r.rid)
        if b is None or b.rate != r.rate or b.tg_pow != r.tg_pow:
            out.append(r)
    return out


def xsec_table(p: Params, te_mid: float, r_mid: float) -> str:
    """Table A.3 as in v0's overview, from v0_1's tables."""
    rows = {r.rid: r for r in build(("all",), "basic")}
    lines = []
    for rid in MO.XSEC_ROWS:
        r = rows[rid]
        gas, kind, prod = r.rate[1:4]
        tab = xsec.rate_table(gas, kind, prod)
        cells = [f"${_k(float(tab(te)))}$" for te in MO.TE_GRID]
        if xsec.in_tail(gas, kind, prod):
            two = xsec.two_t_table(gas, kind, prod)
            f = float(two(te_mid, r_mid * te_mid)) / float(tab(te_mid))
            cells.append(MO._fac(f))
        else:
            cells.append("1")
        tid = rid.replace("_", "\\_")
        label = f"{MO.XSEC_LABELS[rid]}, {r.e_cost_eV:g}\\,eV" \
            if r.e_cost_eV is not None else MO.XSEC_LABELS[rid]
        lines.append(f"{tid} & {label} & " + " & ".join(cells) + "\\\\")
    cap = xsec.capture_energy_table("F2", "attachment")
    cells = [f"{float(cap(te)):.2f}" for te in MO.TE_GRID]
    lines.append(r"\midrule" "\n" r"ATT & energy an attached \el{} takes, eV & "
                 + " & ".join(cells) + r" & \\")
    cells = [f"{1.5 * te:.2f}" for te in MO.TE_GRID]
    lines.append(r" & $\tfrac32 T_e$, for comparison & "
                 + " & ".join(cells) + r" & \\")
    return "\n".join(lines)


def species_table() -> str:
    names = ("Kr", "F2", "F", "Kr*", "Kr+", "KrF(B)", "KrF(C)", "Kr2F*",
             "Kr2*", "Kr2+", "F-")
    lines = []
    for n in names:
        s = SPECIES[n]
        lines.append(f"{MO.TEX[n]} & {MO._neg(f'{e_form(n):.6g}')} & "
                     f"{s.cv:g} & {_tex_text(SPECIES_NOTES[n])}\\\\")
    return "\n".join(lines)


def rates_tex() -> str:
    p = Params()
    res = run(p)
    k = int(res.t.searchsorted(0.5 * p.pulse_s))
    m = res.model
    n_gas = float(sum(res.n[s][k] for s in m.species))
    te_mid = float(res.te[k])
    r_mid = float(m.tail(res.ne[k] / n_gas, res.n["F2"][k] / n_gas))
    basic = build((), "basic")
    groups = [r for r in build(BASIC_GROUPS, "basic")
              if r.rid not in {b.rid for b in basic}]
    parts = [
        "% LIGHT v0_1 reaction-rate tables. GENERATED by make_overview_v0_1.py",
        f"% from LIGHT v0_1 ({__version__}). Do not edit. Rerun the script "
        "instead.",
        r"\providecommand\videf[2]{\expandafter\def\csname vi@#1\endcsname{#2}}",
        r"\videf{rates.teMid}{" + f"{te_mid:.2f}" + "}",
        r"\videf{rates.rMid}{" + f"{r_mid:.2f}" + "}",
        r"\videf{rates.nBasic}{" + str(len(basic)) + "}",
        r"\videf{rates.nGroups}{" + str(len(groups)) + "}",
        r"\videf{rates.groups}{" + ", ".join(BASIC_GROUPS) + "}",
        r"\videf{rates.v0groups}{" + ", ".join(GROUPS) + "}",
        r"\def\virBasic{" + _rate_rows(ordered(basic), SOURCES, p) + "}",
        r"\def\virGroups{" + _rate_rows(ordered(groups), SOURCES, p) + "}",
        r"\def\virVzero{" + _rate_rows(ordered(v0_differences()),
                                       {**SOURCES, **SOURCES_V0}, p) + "}",
        r"\def\virXsec{" + xsec_table(p, te_mid, r_mid) + "}",
        r"\def\virTe{" + " & ".join(f"{t:g}" for t in MO.TE_GRID) + "}",
        r"\def\virSpecies{" + species_table() + "}",
    ]
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# the accuracy numbers
# ---------------------------------------------------------------------------
def _render(n: dict, what: str, tables=()) -> str:
    lines = [
        f"% LIGHT v0_1 accuracy numbers ({what}). GENERATED by "
        "make_overview_v0_1.py",
        f"% from LIGHT v0_1 ({__version__}) and ML/. Do not edit. Rerun the "
        "script instead.",
        r"\providecommand\videf[2]{\expandafter\def\csname vi@#1\endcsname{#2}}",
    ]
    lines += [f"\\videf{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def _ml():
    """ML, bound to v0_1. It fixes its LIGHT version when first imported,
    so the accuracy numbers need a process that set LIGHT_VERSION=v0_1
    first (main() does)."""
    import ml
    if ml.LIGHT_VERSION != VERSION:
        raise RuntimeError(f"ML runs LIGHT {ml.LIGHT_VERSION} in this "
                           "process. Run make_overview_v0_1.py as a program, "
                           "which sets LIGHT_VERSION=v0_1 before ML loads")
    return ml


def _krf_decay(kr_torr, f2_torr, net="basic"):
    """make_overview._krf_decay on v0_1: KrF(B) seeded in Kr/F2 at 300 K,
    no power, and the slow decay of KrF(B) + KrF(C) fitted."""
    tot = kr_torr + f2_torr
    p = Params(network=net, pressure_atm=tot / 760.0, x_f2=f2_torr / tot,
               power_W_cm3=0.0, pulse_s=0.0, t_end_s=400e-9,
               ne_mode="free", ne_cm3=1.0, n_save=2001,
               init_extra_cm3=(("KrF(B)", 1e12),))
    r = Model(p).run()
    n = r.n["KrF(B)"] + r.n.get("KrF(C)", 0.0)
    rel = n / n[0]
    sel = (r.t > 20e-9) & (rel < 1e-3) & (rel > 1e-8)
    return -np.polyfit(r.t[sel], np.log(n[sel]), 1)[0]


def _kr2_split() -> dict:
    """What the optional Kr2* group does at the reference point, taken
    apart: the 248-nm share without it, with it, and with it but without
    Kr2* + F -> KrF(B) + Kr (K2F_F, Kannari's estimate), and the part of the
    gain that arrives after F2 burns out."""
    from v0_1.report import f2_burnout_time

    def run_share(groups=(), drop=None, k_f=None):
        m = Model(Params(groups=groups))
        if drop:
            m.compiled = [c for c in m.compiled if c.row.rid != drop]
        for c in m.compiled:
            if k_f is not None and c.row.rid == "K2F_F":
                c.rate = k_f
        res = m.run()
        a = accounting(res)
        tb = f2_burnout_time(res)
        tb = res.t[-1] if tb is None else tb
        hv = next(c.row.hv_eV for c in m.compiled if c.row.rid == "G4")
        g4 = res.events["G4"]
        after = (g4[-1] - float(np.interp(tb, res.t, g4))) * hv
        return (100.0 * g4[-1] * hv / a["e_dep"],
                100.0 * after / a["e_dep"], res, tb)
    off, off_after, _, _ = run_share()
    on, on_after, res, tb = run_share(("kr2",))
    no_f = run_share(("kr2",), drop="K2F_F")[0]
    # the rate is no bound, and the share saturates: a slower and a faster
    # rate than Kannari's (1e-9 is the top of SRI's range for Ar2* + F)
    k_lo, k_hi = 1.0e-11, 1.0e-9
    lo = run_share(("kr2",), k_f=k_lo)[0]
    hi = run_share(("kr2",), k_f=k_hi)[0]
    rows = {c.row.rid: c.row for c in res.model.compiled}
    k_f = rows["K2F_F"].rate[1]
    k_rad = rows["KR2_RAD"].rate[1]
    f_atoms = float(np.interp(tb, res.t, res.n["F"]))
    gain, gain_after = on - off, on_after - off_after
    return {"kr2.off": f"{off:.2f}", "kr2.on": f"{on:.2f}",
            "kr2.noF": f"{no_f:.2f}", "kr2.gain": f"{gain:.2f}",
            "kr2.after": f"{gain_after:.2f}",
            "kr2.afterpc": f"{100.0 * gain_after / gain:.0f}",
            "kr2.burn": f"{tb * 1e9:.0f}",
            "kr2.F": f"${_k(f_atoms)}$", "kr2.k": f"${_k(k_f)}$",
            "kr2.ratio": f"{k_f * f_atoms / k_rad:.0f}",
            "kr2.lo.k": f"${_k(k_lo)}$", "kr2.lo": f"{lo:.2f}",
            "kr2.hi.k": f"${_k(k_hi)}$", "kr2.hi": f"{hi:.2f}"}


def _tahoma_window(m):
    """Run a Tahoma model; the 200-400 ns window's g0, Tg, n_e and Te."""
    from ml import forward
    r = m.run()
    win = (r.t >= forward.WINDOW[0]) & (r.t <= forward.WINDOW[1])
    nb = float(np.mean(r.n["KrF(B)"][win]))
    return {"g0": forward.SIGMA_B * nb * 100.0,
            "tg": float(np.mean(r.tg[win])), "ne": float(np.mean(r.ne[win])),
            "te": float(np.mean(r.te[win]))}


def fast_numbers() -> dict:
    _ml()
    from ml import data, forward, formation, params, pure, pure_data, score
    from v0_1 import two_level as L
    n = {}
    # the photolysis decays against Xu 1993's printed constants
    for kr, meas in MO.XU_1993:
        n[f"xu.{kr:.0f}"] = MO._x(_krf_decay(kr, 1.0), meas)
        n[f"xu.{kr:.0f}.v0"] = MO._x(_krf_decay(kr, 1.0, "v0"), meas)
    for kr, meas in MO.XU_1993:
        n[f"xu.{kr:.0f}.meas"] = f"${_k(meas)}$"
    # the basic network's B/C decay against Gadomski's analytic eq. 5
    c_light = L.Consts.from_light()
    dev = 0.0
    for kr, f2 in MO.EQ5_POINTS:
        slow = L.eigen(c_light, kr=kr * L.TORR_300K, f2=f2 * L.TORR_300K)[1]
        dev = max(dev, abs(_krf_decay(kr, f2) / slow - 1.0))
    mant, ex = f"{dev:.1e}".split("e")
    n["eq5.dev"] = f"${mant}{{\\times}}10^{{{int(ex)}}}$"
    n["gap"] = f"{c_light.gap_cm1:.0f}"
    # the ledger closure in every mode
    base = Params()
    runs = [base.replace(groups=(g,)) for g in BASIC_GROUPS]
    runs += [base, base.replace(groups=("all",)),
             base.replace(ne_mode="clamped"),
             base.replace(groups=("all",), ne_mode="clamped"),
             base.replace(trimer_tg_pow=-0.5),
             base.replace(network="v0"),
             base.replace(network="v0", groups=("all",))]
    worst = max(abs(accounting(run(p))["closure"]) for p in runs)
    n["closure.max"] = f"$10^{{{math.ceil(math.log10(worst))}}}$"
    n["closure.runs"] = str(len(runs))
    # the hottest corner of the notes' ranges: recombination returns the
    # ionization energy as heat
    hot = base.replace(pressure_atm=0.1, power_W_cm3=5e5)
    n["hot.p"] = f"{hot.pressure_atm:g}"
    n["hot.P"] = f"{hot.power_W_cm3 / 1e3:.0f}"
    n["hot.tg"] = f"{run(hot).tg[-1]:.0f}"
    # ML's scorecard for v0_1's constants (k_e re-derived as ML does)
    c = forward.with_rederived_k_e(params.light_now())
    for ds in data.all_datasets():
        if ds.kind in data.DEVICE_KINDS:
            continue
        v = score.kinetic_score(ds, c)
        n[f"sc.{ds.name}"] = (MO._pct_signed(v) if len(ds.values) == 1
                              else f"{v:.1f}")
    for key, name in (("ys", "ys_ir300"), ("ysr", "ys_ratio_230_300")):
        d = data.by_name(name)[0]
        n[f"{key}.meas"] = f"{d.values[0]:.2f}"
        n[f"{key}.err"] = f"{d.sigma:.2f}"
    rms = [score.kinetic_score(d, c) for d in data.by_role("train")
           if len(d.values) > 1]
    n["sc.train.lo"], n["sc.train.hi"] = f"{min(rms):.1f}", f"{max(rms):.1f}"
    n["sc.bs.meas"] = f"{data.by_name('bs_tau_b')[0].values[0] * 1e9:.1f}"
    n["sc.bs.err"] = f"{data.by_name('bs_tau_b')[0].sigma * 1e9:.1f}"
    # Tahoma under the beam drive, as ML scores it (basic + kr2), a score
    o = formation.observables(c)
    meas = {d.name: d for d in data.by_role("device")}
    g0, eta, isat = (meas[k].values[0] for k in ("tahoma_g0",
                                                  "tahoma_formation",
                                                  "tahoma_isat"))
    n.update({
        "tah.g0.meas": f"{g0:.2f}", "tah.g0": f"{o['g0']:.2f}",
        "tah.g0.x": MO._x(o["g0"], g0),
        "tah.isat.meas": f"{isat:.1f}", "tah.isat": f"{o['isat']:.2f}",
        "tah.isat.err": f"{meas['tahoma_isat'].sigma * isat:.1f}",
        "tah.isat.x": MO._x(o["isat"], isat),
        "tah.eta.meas": f"{eta:.2f}", "tah.eta": f"{o['eta']:.3f}",
        "tah.eta.x": MO._x(o["eta"], eta),
        "tah.form": f"{o['formation']:.3f}",
        "tah.tg": f"{o['tg']:.0f}",
        "tah.pin": f"{forward.TAHOMA_COND['power_W_cm3'] / 1e3:.0f}",
        "tah.torr": f"{forward.TAHOMA_COND['pressure_torr']:.0f}",
        "tah.f2torr": f"{forward.TAHOMA_COND['f2_torr']:.1f}",
        "tah.groups": ", ".join(forward.tahoma_params().groups),
    })
    # the same beam with v0_1's default network (no Kr2* group)
    drive = forward._drive()
    p_def = forward.tahoma_params(kr2=False)
    d = _tahoma_window(forward.impose(drive.BeamModel(p_def), c))
    n["tah.def.g0"] = f"{d['g0']:.2f}"
    n["tah.def.x"] = MO._x(d["g0"], g0)
    # the Kr2* group without its unmeasured Kr2* + F -> KrF(B) + Kr
    m = forward.impose(drive.BeamModel(forward.tahoma_params()), c)
    m.compiled = [cc for cc in m.compiled if cc.row.rid != "K2F_F"]
    d = _tahoma_window(m)
    n["tah.kr2noF.g0"] = f"{d['g0']:.2f}"
    n["tah.kr2noF.x"] = MO._x(d["g0"], g0)
    # the same at the reference point, with the part after F2 burnout
    n.update(_kr2_split())
    # the trimer's weak temperature law (pure-Kr T^-1/2) at Tahoma
    p_tri = forward.tahoma_params().replace(trimer_tg_pow=-0.5)
    d = _tahoma_window(forward.impose(drive.BeamModel(p_tri), c))
    n["tah.tri.g0"] = f"{d['g0']:.2f}"
    n["tah.tri.x"] = MO._x(d["g0"], g0)
    n["tri.lo.pow"] = "0.5"
    n["tri.klimek.hi"] = "425"
    # v0_1's own input at Tahoma's conditions: the power into the
    # electrons, no beam
    own = _tahoma_window(forward.impose(Model(forward.tahoma_params()), c))
    beam = _tahoma_window(forward.impose(drive.BeamModel(
        forward.tahoma_params()), c))
    for tag, d in (("own", own), ("beam", beam)):
        n[f"tah.{tag}.g0"] = f"{d['g0']:.2f}"
        n[f"tah.{tag}.te"] = f"{d['te']:.1f}"
        n[f"tah.{tag}.ne"] = f"${_k(d['ne'])}$"
        n[f"tah.{tag}.tg"] = f"{d['tg']:.0f}"
    n["tah.own.x"] = MO._x(own["g0"], g0)
    # the pure-gas checks (the F2-free half of the chain). v0_1 sets R11p
    # and K15 from the very measurements the second and third compare with
    for g in ("Kr", "Ar"):
        x = pure.supply(g)
        n[f"pure.{g}.x"] = MO._x(x["per_eV_sf"], x["per_eV_meas"])
    rows, k15, conv = pure.kr_rates()
    ratios = [light / m for _, light, m in rows]
    n["pure.dr.lo"], n["pure.dr.hi"] = (f"{min(ratios):.2f}",
                                        f"{max(ratios):.2f}")
    n["pure.dr.te.lo"] = f"{rows[0][0]:.2g}"
    n["pure.dr.te.hi"] = f"{rows[-1][0]:.2g}"
    n["pure.k15"] = f"${_k(k15)}$"
    n["pure.conv.lo"], n["pure.conv.hi"] = (f"${_k(conv[1])}$",
                                            f"${_k(conv[0])}$")
    v = pure.kr_vuv(3.3e2)
    vk = next(x for x in pure_data.VUV if x.link == "vuv_Kr")
    n.update({"pure.vuv": f"{v['vuv']:.2f}", "pure.vuv.meas": f"{vk.value:.2f}",
              "pure.vuv.err": f"{vk.err:.2f}",
              "pure.vuv.x": MO._x(v["vuv"], vk.value),
              "pure.vuv.ideal": f"{pure.ideal_vuv('Kr'):.2f}",
              "pure.vuv.nothermal":
                  f"{pure.kr_vuv(3.3e2, thermal=False)['vuv']:.2f}",
              "pure.thermal": f"{v['thermal']:.3f}"})
    return n


def slow_numbers() -> tuple:
    """(numbers, tables): the \\videf numbers and the finished table rows."""
    _ml()
    from ml import (data, eedf, fit, formation, formation_data, forward,
                    params, score, swarm)
    n, tables = {}, []
    # the posterior of the B/C block: its constants are v0's (the train
    # sets are at 300 K, where the B-C gap does not enter), and its Tahoma
    # gain is v0_1's
    train = data.by_role("train")
    x0 = fit.fit(train)
    chain, _ = fit.mcmc(train, x0, n_steps=2500)
    s = fit.summarize(chain)
    for name, scale, fmt in (("k_BC", 1.0, "k"), ("K_eq", 1.0, ".2f"),
                             ("tau_C", 1e9, ".0f"), ("tau_B", 1e9, ".2f"),
                             ("k_F2", 1.0, "k"), ("k_F2_C", 1.0, "k"),
                             ("k3_Kr", 1.0, "k"), ("k3_Kr_C", 1.0, "k")):
        for part, v in zip(("lo", "mid", "hi"), s[name]):
            v = v * scale
            n[f"fit.{name}.{part}"] = f"${_k(v)}$" if fmt == "k" else \
                f"{v:{fmt}}"
    for part, v in zip(("lo", "mid", "hi"), s["b_share_1atm"]):
        n[f"fit.bshare.{part}"] = f"{v:.2f}"
    n["fit.bshare.light"] = f"{forward.b_share_1atm(params.light_now()):.2f}"
    n["fit.n"] = str(len(params.NAMES))
    n_draws = 40
    lo, mid, hi = score.device_posterior(chain, n_draws=n_draws)
    n.update({"post.g0.lo": f"{lo:.2f}", "post.g0.mid": f"{mid:.2f}",
              "post.g0.hi": f"{hi:.2f}", "post.n": str(n_draws)})
    # each unmeasured link at Tahoma, one at a time from v0_1's chain
    c = forward.with_rederived_k_e(params.light_now())
    have = [kw for _, kw in formation_data.UNMEASURED]
    for kw, key, _, shown, _ in MO.LINKS:
        if kw not in have:
            raise KeyError(f"{kw} is no longer in formation_data.UNMEASURED")
    if len(have) != len(MO.LINKS):
        raise KeyError("formation_data.UNMEASURED has links this overview "
                       "does not show; add them to make_overview.LINKS")
    chain_m = formation.LIGHT_CHAIN.with_(**formation_data.MEASURED)
    base = formation.observables(c, chain_m)["eta"]
    n["link.base"] = f"{base:.3f}"
    minor = 0.0
    rows = []
    for kw, key, label, shown, what in MO.LINKS:
        eta = formation.observables(c, chain_m.with_(**kw))["eta"]
        n[f"link.{key}"] = f"{eta:.3f}"
        n[f"link.{key}.d"] = MO._neg(f"{eta - base:+.3f}")
        if label is None:
            minor = max(minor, abs(eta - base))
            continue
        name = list(kw)[0]
        light = MO._light_link(name, getattr(chain_m, name)) if label else ""
        rows.append(f"{label} & {light} & {shown} & {what} & {eta:.3f} & "
                    f"{n[f'link.{key}.d']}\\\\")
    n["link.minor"] = f"{math.ceil(minor * 1000) / 1000:.3f}"
    tables.append("\\def\\viLinks{" + "\n".join(rows) + "}")
    # the electron closure against the kinetic loop at Tahoma
    loop, b, _ = eedf.self_consistent(c)
    own = forward.run_device(c)
    for tag, d in (("loop", loop), ("light", own)):
        n[f"el.{tag}.te"] = f"{d['te']:.2f}"
        n[f"el.{tag}.ne"] = f"${_k(d['ne'])}$"
        n[f"el.{tag}.g0"] = f"{d['g0']:.2f}"
        n[f"el.{tag}.form"] = f"{d['formation']:.3f}"
    for q in ("te", "ne", "g0", "formation"):
        dev = 100.0 * abs(own[q] / loop[q] - 1.0)
        n[f"el.dev.{q}"] = f"{dev:.1f}"
    n["el.eatt"] = f"{b['e_att']:.2f}"
    n["el.eatt.32"] = f"{b['mean']:.2f}"
    n["el.rpow"] = f"{-eedf.r11p_pow():.2f}"
    n["el.rcap"] = f"{1.5 + eedf.r11p_pow():.2f}"
    # the cross sections against swarm data (the same files as v0)
    by = {ds.name: ds for ds in swarm.load()}
    r = swarm.predict(by["kr_w_nakamura"]) / by["kr_w_nakamura"].y
    n["sw.krw.lo"], n["sw.krw.hi"] = f"{r.min():.2f}", f"{r.max():.2f}"
    lo, hi = [], []
    for name in ("f2_katt_nygaard", "f2_katt_tj300", "f2_katt_tj500",
                 "f2_katt_chen", "f2_katt_sb"):
        q = swarm.predict(by[name]) / by[name].y
        lo.append(float(q.min()))
        hi.append(float(q.max()))
    n["sw.f2.lo"], n["sw.f2.hi"] = f"{min(lo):.2f}", f"{max(hi):.2f}"
    return n, tables


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--fast", action="store_true",
                    help="write the rate tables and the fast numbers only")
    ap.add_argument("--check", choices=("rates", "fast", "slow"),
                    help="compare one file with a fresh run, write nothing")
    a = ap.parse_args(argv)
    # before ML is imported anywhere in this process (see _ml)
    os.environ["LIGHT_VERSION"] = VERSION
    if a.check:
        fresh = {"rates": rates_tex,
                 "fast": lambda: _render(fast_numbers(), "fast"),
                 "slow": lambda: _render(*_slow_parts(), )}[a.check]()
        path = {"rates": OUT_RATES, "fast": OUT_FAST, "slow": OUT_SLOW}[a.check]
        if path.read_text() != fresh:
            print(f"{path.name} is stale: rerun make_overview_v0_1.py")
            return 1
        print(f"{path.name} is current")
        return 0
    OUT_RATES.write_text(rates_tex())
    print(f"wrote {OUT_RATES}")
    OUT_FAST.write_text(_render(fast_numbers(), "fast"))
    print(f"wrote {OUT_FAST}")
    if not a.fast:
        OUT_SLOW.write_text(_render(*_slow_parts()))
        print(f"wrote {OUT_SLOW}")
    return 0


def _slow_parts():
    n, tables = slow_numbers()
    return n, "slow", tables


if __name__ == "__main__":
    sys.exit(main())
