"""Numbers and tables for accuracy_overview_v0_2.tex, v0_2's accuracy.

The v0_1 overview's generator (make_overview_v0_1.py) run on v0_2, with
the checks v0_2 needs besides: its bands, each end at Tahoma, the Kr2*
route at the reference point and in pure krypton. No number in the
overview is typed by hand. This script writes three files, each one
\\widef{key}{value} per number or a finished table (a macro of its own, so
one document can hold v0_1's numbers and v0_2's):

  overview_v0_2_rates.tex  Appendix A, the reaction rates of v0_2's rows
                           and its bands
  overview_v0_2_fast.tex   the checks that take seconds: the photolysis
                           decays, the ledger closure, ML's scorecard with
                           Tahoma (a score, never a fit target), Tahoma at
                           both ends of every band and on v0_2's own power
                           input, the Kr2* route at the reference point,
                           and the pure-gas checks
  overview_v0_2_slow.tex   the checks that take minutes: the posterior's
                           Tahoma gain, each unmeasured formation link, the
                           electron closure against the kinetic EEDF loop,
                           and the cross sections against swarm data

The accuracy numbers come from ML/, which runs v0_2 when LIGHT_VERSION
names it. ML fixes its version when first imported, so the numbers are
made in a process of their own: main() sets LIGHT_VERSION=v0_2 before ML
is imported, and every ML-backed function refuses to run on another
version. v0_2/tests/test_docs.py compares the files through --check. The
map's numbers (system_map_v0_2_numbers.tex, make_system_map_v0_2.py) hold
the bands at the reference point, which the overview reuses.

    python3 LIGHT/docs/make_system_map_v0_2.py             # about 1 min
    python3 LIGHT/docs/make_overview_v0_2.py               # all, ~7 min
    python3 LIGHT/docs/make_overview_v0_2.py --fast        # rates + fast
    python3 LIGHT/docs/make_overview_v0_2.py --check fast  # compare only
    cd LIGHT/docs
    latexmk -pdf system_map_v0_2.tex
    latexmk -pdf accuracy_overview_v0_2.tex
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_2.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_2.py").is_file())
ML_DIR = LIGHT.parent / "ML"
for _p in (str(LIGHT), str(HERE), str(ML_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                            # noqa: E402

import make_overview as MO                                    # noqa: E402
import make_overview_v0_1 as MO1                              # noqa: E402
from make_system_map_basic import _k                          # noqa: E402
from make_system_map_v0_2 import (BAND_LABELS, BAND_ORDER,    # noqa: E402
                                  UNITS)
from v0_2 import __version__, xsec                            # noqa: E402
from v0_2.model import Model, Params, run                     # noqa: E402
from v0_2.network import (BANDS, BASIC_GROUPS, GROUPS,        # noqa: E402
                          build)
from v0_2.report import accounting, f2_burnout_time           # noqa: E402
from v0_2.species import SPECIES, e_form                      # noqa: E402

VERSION = "v0_2"
OUT_RATES = HERE / "overview_v0_2_rates.tex"
OUT_FAST = HERE / "overview_v0_2_fast.tex"
OUT_SLOW = HERE / "overview_v0_2_slow.tex"
TRIMER_WEAK = -0.5                        # the pure-Kr T^-1/2 law's exponent

# ---------------------------------------------------------------------------
# Appendix A: the reaction rates
# ---------------------------------------------------------------------------
# the order the tables list the rows in (v0_1's, with Q2FF after Q2F). A
# row missing here is appended at the end, never dropped.
ORDER = tuple(r for rid in MO1.ORDER for r in
              ((rid, "Q2FF") if rid == "Q2F" else (rid,)))


def ordered(rows) -> list:
    """Rows in ORDER. A C-state copy (bc group, 'G5_C') follows its B row."""
    rank = {rid: i for i, rid in enumerate(ORDER)}

    def key(r):
        if r.rid not in rank and r.rid.endswith("_C"):
            return (rank.get(r.rid[:-2], len(ORDER)), 1, r.rid)
        return (rank.get(r.rid, len(ORDER)), 0, r.rid)
    return sorted(rows, key=key)


# Short sources, one per row, condensed from each row's full citation in
# v0_2/network.py. Every number here must appear in that citation
# (test_docs.py checks it), so the table cannot drift from the code. v0_1's
# where the row is v0_1's, and the rows v0_2 changed or added.
SOURCES = {
    **MO1.SOURCES,
    "G2": "Kannari 1985 Table I Penning block, from Leichner and Ericson "
          "1974. Band down to Werner's 1.5e-10, via Jancaitis 1983",
    "Q2F": "Xu 1993 Table I, 1.2 +- 0.2e-10. Band up to Quigley and "
           "Hughes' 4.3e-10",
    "Q2FF": "no Kr measurement, so 0 unless a band run sets it, 0 to "
            "1.8e-10. The top is Jancaitis 1983's Ar2F* analog",
    "AK8_Kr": "Ueda 1984 Table V, 1.7 +- 0.2e-32 at 0.5-1 atm. T_g law "
              "from Wieme and Lenaerts 1981, 300-451 K. Band up to Hughes "
              "1977, 5.5e-32",
    "KR2_RAD": "Morikawa 1989, 270 +- 5 ns, and five other measurements",
    "NK7": "half of Quigley and Hughes 1978, 2.1 +- 0.2e-10. The even split "
           "is Kannari's. Band up to Chen and Payne's 5.0e-10",
    "K2F_BR": "the other half of Quigley and Hughes 1978",
    "K2F_F": "no measurement, so 0 unless a band run sets it, 0 to 1e-9. "
             "Kannari 1985's 3.0e-10 is an estimate",
    "PEN_KR2": "no measurement, so 0 unless a band run sets it, 0 to "
               "5.0e-10. Kannari 1985 took 5.0e-10 from an unpublished 1973 "
               "SRI report",
    "K26": "Kannari 1985, Flannery turnover with Debye screening, and "
           "(T_g/300)^-3 on its three-body part from Thomson theory. "
           "Measured with F- only below 300 Torr, Jacob et al. 1978. No "
           "measurement near 1 atm",
}
# v0_2's v0 network (the notes' list as written, REBUILD's copies) where it
# differs from the basic network
SOURCES_V0 = {
    **MO1.SOURCES_V0,
    "G2": MO.SOURCES["G2"],
    "K2F_F": MO.SOURCES["K2F_F"],
    "PEN_KR2": MO1.SOURCES["PEN_KR2"],
}
SPECIES_NOTES = {
    **MO1.SPECIES_NOTES,
    "Kr2+": "IP(Kr) less D0(Kr2+) = 1.15 eV, measured by photoionization, "
            "Ng et al. 1977",
}


# the excimers the band citations name, set as math before v0_1's rules
# (which would read their digits as numbers)
_CHEM = (("Kr2F*", "Kr$_2$F$^*$"), ("Ar2F*", "Ar$_2$F$^*$"),
         ("Kr2*", "Kr$_2^*$"), ("Ar2*", "Ar$_2^*$"), ("Kr*", "Kr$^*$"))


def _tex_text(s: str) -> str:
    """A short source as LaTeX text (v0_1's rules, the excimers first). A
    band's citation names authors with '&', a column break in a table, so
    it reads 'and' here, and a number's leading minus is set as math."""
    s = s.replace(" & ", " and ")
    s = re.sub(r"(?<=\s)-(\d)", r"\\ensuremath{-}\1", s)
    for plain, tex in _CHEM:
        s = s.replace(plain, tex)
    return MO1._tex_text(s)


def rate_text(r, p: Params) -> str:
    """The rate coefficient at 300 K as printed in the table. A band row
    prints the value v0_2 runs by default and the band's name."""
    if r.rate[0] == "param":
        name = r.rate[1]
        frac = float(r.rate[2]) if len(r.rate) > 2 else 1.0
        v = frac * BANDS[name].default
        half = r", half" if frac != 1.0 else ""
        tex = name.replace("_", r"\_")
        if v == 0.0:
            return f"0, band {tex}"
        return f"${_k(v)}$ {MO._unit(r)}, band {tex}{half}"
    return MO1.rate_text(r, p)


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
                   f"{MO1.tdep_text(r)} & {_tex_text(src)}\\\\")
    return "\n".join(out)


def v0_differences() -> list:
    """v0_2's v0 network rows (every group) that its basic network lacks or
    runs with another value."""
    basic = {r.rid: r for r in build(("all",), "basic")}
    out = []
    for r in build(("all",), "v0"):
        b = basic.get(r.rid)
        if b is None or b.rate != r.rate or b.tg_pow != r.tg_pow:
            out.append(r)
    return out


def xsec_table(te_mid: float, r_mid: float) -> str:
    """Table A.3 as in v0's overview, from v0_2's tables."""
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


def _band_num(b, v) -> str:
    if b.unit == "exponent":
        return f"${v:.2f}$"
    return "0" if v == 0.0 else f"${_k(v)}$"


def bands_table() -> str:
    """Table A.4, every band: its process, default, range and source."""
    lines = []
    for name in BAND_ORDER:
        b = BANDS[name]
        unit = {"cm3/s": r"cm$^3$s$^{-1}$", "cm6/s": r"cm$^6$s$^{-1}$",
                "exponent": "exponent"}[b.unit]
        lines.append(f"{name.replace('_', chr(92) + '_')} & "
                     f"{BAND_LABELS[name]} & {_band_num(b, b.default)} & "
                     f"{_band_num(b, b.lo)} to {_band_num(b, b.hi)} {unit} & "
                     f"{_tex_text(b.cite)}\\\\")
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
    n_band_rows = sum(1 for r in basic if r.rate[0] == "param")
    parts = [
        "% LIGHT v0_2 reaction-rate tables. GENERATED by make_overview_v0_2.py",
        f"% from LIGHT v0_2 ({__version__}). Do not edit. Rerun the script "
        "instead.",
        r"\providecommand\widef[2]{\expandafter\def\csname wi@#1\endcsname{#2}}",
        r"\widef{rates.teMid}{" + f"{te_mid:.2f}" + "}",
        r"\widef{rates.rMid}{" + f"{r_mid:.2f}" + "}",
        r"\widef{rates.nBasic}{" + str(len(basic)) + "}",
        r"\widef{rates.nBandRows}{" + str(n_band_rows) + "}",
        r"\widef{rates.nGroups}{" + str(len(groups)) + "}",
        r"\widef{rates.groups}{" + ", ".join(BASIC_GROUPS) + "}",
        r"\widef{rates.v0groups}{" + ", ".join(GROUPS) + "}",
        r"\def\wirBasic{" + _rate_rows(ordered(basic), SOURCES, p) + "}",
        r"\def\wirGroups{" + _rate_rows(ordered(groups), SOURCES, p) + "}",
        r"\def\wirVzero{" + _rate_rows(ordered(v0_differences()),
                                       {**SOURCES, **SOURCES_V0}, p) + "}",
        r"\def\wirXsec{" + xsec_table(te_mid, r_mid) + "}",
        r"\def\wirTe{" + " & ".join(f"{t:g}" for t in MO.TE_GRID) + "}",
        r"\def\wirBands{" + bands_table() + "}",
        r"\def\wirSpecies{" + species_table() + "}",
    ]
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# the accuracy numbers
# ---------------------------------------------------------------------------
def _render(n: dict, what: str, tables=()) -> str:
    lines = [
        f"% LIGHT v0_2 accuracy numbers ({what}). GENERATED by "
        "make_overview_v0_2.py",
        f"% from LIGHT v0_2 ({__version__}) and ML/. Do not edit. Rerun the "
        "script instead.",
        r"\providecommand\widef[2]{\expandafter\def\csname wi@#1\endcsname{#2}}",
    ]
    lines += [f"\\widef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def _ml():
    """ML, bound to v0_2. It fixes its LIGHT version when first imported,
    so the accuracy numbers need a process that set LIGHT_VERSION=v0_2
    first (main() does)."""
    import ml
    if ml.LIGHT_VERSION != VERSION:
        raise RuntimeError(f"ML runs LIGHT {ml.LIGHT_VERSION} in this "
                           "process. Run make_overview_v0_2.py as a program, "
                           "which sets LIGHT_VERSION=v0_2 before ML loads")
    return ml


def _krf_decay(kr_torr, f2_torr, net="basic"):
    """make_overview._krf_decay on v0_2: KrF(B) seeded in Kr/F2 at 300 K,
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


def _kannari(name) -> float:
    """Kannari's estimate among a band's published values."""
    return next(v for v, who in BANDS[name].marks if who.startswith("Kannari"))


def _kr2_route() -> dict:
    """The Kr2* route at the reference point: what Kr2* + F -> KrF(B) + Kr
    (band kr2_f, no measurement, 0 by default) would do across its band,
    the part of that gain that comes after F2 burns out, and Kannari's two
    estimates together (v0_1's kr2 group)."""
    def run_share(rates=()):
        res = run(Params(rates=rates))
        a = accounting(res)
        tb = f2_burnout_time(res)
        tb = res.t[-1] if tb is None else tb
        hv = next(c.row.hv_eV for c in res.model.compiled
                  if c.row.rid == "G4")
        g4 = res.events["G4"]
        after = (g4[-1] - float(np.interp(tb, res.t, g4))) * hv
        return {"s248": 100.0 * g4[-1] * hv / a["e_dep"],
                "after": 100.0 * after / a["e_dep"],
                "vuv": 100.0 * a["bands"]["147 nm (Kr2* VUV)"] / a["e_dep"],
                "res": res, "tb": tb}
    base = run_share()
    k_kann = _kannari("kr2_f")
    k_lo, k_hi = 1.0e-11, BANDS["kr2_f"].hi
    lo = run_share((("kr2_f", k_lo),))
    kann = run_share((("kr2_f", k_kann),))
    hi = run_share((("kr2_f", k_hi),))
    both = run_share((("kr2_f", k_kann), ("kr2_pen", _kannari("kr2_pen"))))
    res = base["res"]
    rows = {c.row.rid: c for c in res.model.compiled}
    k_rad = float(rows["KR2_RAD"].rate)
    f_atoms = float(np.interp(base["tb"], res.t, res.n["F"]))
    gain = hi["s248"] - base["s248"]
    gain_after = hi["after"] - base["after"]
    return {"kr2.def": f"{base['s248']:.2f}", "kr2.vuv": f"{base['vuv']:.1f}",
            "kr2.burn": f"{base['tb'] * 1e9:.0f}",
            "kr2.F": f"${_k(f_atoms)}$",
            "kr2.lo.k": f"${_k(k_lo)}$", "kr2.lo": f"{lo['s248']:.2f}",
            "kr2.kann.k": f"${_k(k_kann)}$", "kr2.kann": f"{kann['s248']:.2f}",
            "kr2.kann.ratio": f"{k_kann * f_atoms / k_rad:.0f}",
            "kr2.hi.k": f"${_k(k_hi)}$", "kr2.hi": f"{hi['s248']:.2f}",
            "kr2.hi.vuv": f"{hi['vuv']:.1f}",
            "kr2.hi.ratio": f"{k_hi * f_atoms / k_rad:.0f}",
            "kr2.hi.afterpc": f"{100.0 * gain_after / gain:.0f}",
            "kr2.both": f"{both['s248']:.2f}",
            "kr2.pen.k": f"${_k(_kannari('kr2_pen'))}$"}


def _tahoma_window(m):
    """Run a Tahoma model; the 200-400 ns window's g0, Tg, n_e and Te, and
    the share of the F2 fill left when the window closes."""
    from ml import forward
    r = m.run()
    win = (r.t >= forward.WINDOW[0]) & (r.t <= forward.WINDOW[1])
    nb = float(np.mean(r.n["KrF(B)"][win]))
    f2 = float(np.interp(forward.WINDOW[1], r.t, r.n["F2"])) / r.n["F2"][0]
    return {"g0": forward.SIGMA_B * nb * 100.0,
            "tg": float(np.mean(r.tg[win])), "ne": float(np.mean(r.ne[win])),
            "te": float(np.mean(r.te[win])), "f2left": 100.0 * f2}


def _with_band(p, name, v):
    """Params with one band at v (the trimer law is its own field)."""
    if name == "trimer_tg_pow":
        return p.replace(trimer_tg_pow=v)
    return p.replace(rates=((name, v),))


# pure krypton at Morozov's low pump and at the notes' power: the conditions
# the v0_2 README quotes (665 Torr; a 1-us pulse at 300 kW/cm3)
VUV_HIGH = {"power_W_cm3": 3.0e5, "pulse": 1.0e-6, "torr": 665.0}


def _pure_bands() -> dict:
    """ML's pure-Kr VUV yield with each band at each end: its largest move
    at Morozov's pump, and the span at 300 kW/cm3. kr_vuv sets Kr* + 2Kr
    from its k3 argument, so that band goes through k3, not rates."""
    from ml import pure
    base_lo = pure.kr_vuv(3.3e2)["vuv"]
    base_hi = pure.kr_vuv(**VUV_HIGH)["vuv"]
    moves, hi_vals = [], {}
    for name in BAND_ORDER:
        if name == "trimer_tg_pow":
            continue            # no KrF in pure Kr: the trimer never forms
        b = BANDS[name]
        for v in (b.lo, b.hi):
            kw = ({"k3": v} if name == "kr2_form"
                  else {"rates": ((name, v),)})
            lo = pure.kr_vuv(3.3e2, **kw)["vuv"]
            hi = pure.kr_vuv(**VUV_HIGH, **kw)["vuv"]
            moves.append(abs(lo - base_lo))
            hi_vals[(name, v)] = hi
    low_end = min(hi_vals, key=hi_vals.get)
    high_end = max(hi_vals, key=hi_vals.get)

    def what(end):
        name, v = end
        b = BANDS[name]
        return f"{TEXT_LABELS[name]} at {_band_num(b, v)}{UNITS[b.unit]}"
    return {"pure.band.move": f"{max(moves):.3f}",
            "pure.hi.P": f"{VUV_HIGH['power_W_cm3'] / 1e3:.0f}",
            "pure.hi.pulse": f"{VUV_HIGH['pulse'] * 1e6:g}",
            "pure.hi.torr": f"{VUV_HIGH['torr']:.0f}",
            "pure.hi.def": f"{base_hi:.2f}",
            "pure.hi.lo": f"{hi_vals[low_end]:.2f}",
            "pure.hi.lo.what": what(low_end),
            "pure.hi.hi": f"{hi_vals[high_end]:.2f}",
            "pure.hi.hi.what": what(high_end)}


# the bands in running text: the table labels, less the 300-K note
TEXT_LABELS = {**BAND_LABELS, "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr"}


def fast_numbers() -> tuple:
    """(numbers, tables)."""
    _ml()
    from ml import data, forward, formation, params, pure, pure_data, score
    from v0_2 import two_level as L
    n, tables = {}, []
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
    # the ledger closure in every mode, with the three unmeasured rows at
    # the tops of their bands among them
    base = Params()
    tops = tuple((b.name, b.hi) for b in BANDS.values() if b.default == 0.0)
    runs = [base.replace(groups=(g,)) for g in BASIC_GROUPS]
    runs += [base, base.replace(groups=("all",)),
             base.replace(ne_mode="clamped"),
             base.replace(groups=("all",), ne_mode="clamped"),
             base.replace(trimer_tg_pow=TRIMER_WEAK),
             base.replace(rates=tops),
             base.replace(network="v0"),
             base.replace(network="v0", groups=("all",))]
    worst = max(abs(accounting(run(p))["closure"]) for p in runs)
    n["closure.max"] = f"$10^{{{math.ceil(math.log10(worst))}}}$"
    n["closure.runs"] = str(len(runs))
    # the hottest corner of the notes' ranges
    hot = base.replace(pressure_atm=0.1, power_W_cm3=5e5)
    n["hot.p"] = f"{hot.pressure_atm:g}"
    n["hot.P"] = f"{hot.power_W_cm3 / 1e3:.0f}"
    n["hot.tg"] = f"{run(hot).tg[-1]:.0f}"
    # ML's scorecard for v0_2's constants (k_e re-derived as ML does)
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
    # Tahoma under the beam drive, as ML scores it. v0_2's default network
    # carries the Kr2* route, so ML's configuration is the default network
    o = formation.observables(c)
    meas = {d.name: d for d in data.by_role("device")}
    g0, eta, isat = (meas[k].values[0] for k in ("tahoma_g0",
                                                  "tahoma_formation",
                                                  "tahoma_isat"))
    p_t = forward.tahoma_params()
    if p_t.groups or p_t.rates or p_t.trimer_tg_pow is not None:
        raise RuntimeError("ML's Tahoma configuration is no longer v0_2's "
                           "default network: say which one it is")
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
        "tah.seed": f"${_k(p_t.ne_cm3)}$",
    })
    # every band at both ends at Tahoma, one band at a time, the rest at
    # their defaults (scores: device data set no band)
    drive = forward._drive()
    ref = _tahoma_window(forward.impose(drive.BeamModel(p_t), c))
    n["tah.def.g0"] = f"{ref['g0']:.2f}"
    n["tah.win.lo"] = f"{forward.WINDOW[0] * 1e9:.0f}"
    n["tah.win.hi"] = f"{forward.WINDOW[1] * 1e9:.0f}"
    n["tah.f2left"] = f"{ref['f2left']:.0f}"
    rows, widest = [], (0.0, "")
    for name in BAND_ORDER:
        b = BANDS[name]
        ends = [_tahoma_window(forward.impose(drive.BeamModel(
            _with_band(p_t, name, v)), c))["g0"] for v in (b.lo, b.hi)]
        n[f"tah.bd.{name}.lo"] = f"{ends[0]:.2f}"
        n[f"tah.bd.{name}.hi"] = f"{ends[1]:.2f}"
        span = abs(ends[1] - ends[0])
        if span > widest[0]:
            widest = (span, name)
        rows.append(f"{BAND_LABELS[name]} & {ends[0]:.2f} to {ends[1]:.2f} & "
                    f"\\x{{{MO._x(ends[0], g0)}}} to "
                    f"\\x{{{MO._x(ends[1], g0)}}}\\\\")
    tables.append("\\def\\wiTahBands{" + "\n".join(rows) + "}")
    minor = max(abs(float(n[f"tah.bd.{nm}.{e}"]) - ref["g0"])
                for nm in BAND_ORDER if nm not in ("kr2_f", "trimer_tg_pow")
                for e in ("lo", "hi"))
    n["tah.bd.minor"] = f"{minor:.2f}"
    n["tah.bd.widest"] = BAND_LABELS[widest[1]]
    for name, key in (("kr2_f", "hiF"), ("trimer_tg_pow", "tri")):
        n[f"tah.{key}.g0"] = n[f"tah.bd.{name}.hi"]
        n[f"tah.{key}.x"] = MO._x(float(n[f"tah.bd.{name}.hi"]), g0)
    n["tri.lo.pow"] = f"{-TRIMER_WEAK:g}"
    n["tri.klimek.hi"] = "425"
    # Kannari's two Kr2* estimates together (v0_1's kr2 group) at Tahoma
    kann = (("kr2_f", _kannari("kr2_f")), ("kr2_pen", _kannari("kr2_pen")))
    d = _tahoma_window(forward.impose(drive.BeamModel(p_t.replace(
        rates=kann)), c))
    n["tah.kannari.g0"] = f"{d['g0']:.2f}"
    n["tah.kannari.x"] = MO._x(d["g0"], g0)
    # v0_2's own input at Tahoma's conditions: the power into the
    # electrons, no beam
    own = _tahoma_window(forward.impose(Model(p_t), c))
    for tag, dd in (("own", own), ("beam", ref)):
        n[f"tah.{tag}.g0"] = f"{dd['g0']:.2f}"
        n[f"tah.{tag}.te"] = f"{dd['te']:.1f}"
        n[f"tah.{tag}.ne"] = f"${_k(dd['ne'])}$"
        n[f"tah.{tag}.tg"] = f"{dd['tg']:.0f}"
    n["tah.own.x"] = MO._x(own["g0"], g0)
    # the Kr2* route at the reference point
    n.update(_kr2_route())
    # the pure-gas checks (the F2-free half of the chain). v0_2 sets R11p
    # and K15 from the very measurements the second and third compare with
    for g in ("Kr", "Ar"):
        x = pure.supply(g)
        n[f"pure.{g}.x"] = MO._x(x["per_eV_sf"], x["per_eV_meas"])
    rows_dr, k15, conv = pure.kr_rates()
    ratios = [light / m for _, light, m in rows_dr]
    n["pure.dr.lo"], n["pure.dr.hi"] = (f"{min(ratios):.2f}",
                                        f"{max(ratios):.2f}")
    n["pure.dr.te.lo"] = f"{rows_dr[0][0]:.2g}"
    n["pure.dr.te.hi"] = f"{rows_dr[-1][0]:.2g}"
    n["pure.k15"] = f"${_k(k15)}$"
    n["pure.conv.lo"], n["pure.conv.hi"] = (f"${_k(conv[1])}$",
                                            f"${_k(conv[0])}$")
    v = pure.kr_vuv(3.3e2)
    vk = next(x for x in pure_data.VUV if x.link == "vuv_Kr")
    n.update({"pure.vuv": f"{v['vuv']:.2f}",
              "pure.vuv.meas": f"{vk.value:.2f}",
              "pure.vuv.err": f"{vk.err:.2f}",
              "pure.vuv.x": MO._x(v["vuv"], vk.value),
              "pure.vuv.ideal": f"{pure.ideal_vuv('Kr'):.2f}",
              "pure.vuv.nothermal":
                  f"{pure.kr_vuv(3.3e2, thermal=False)['vuv']:.2f}",
              "pure.thermal": f"{v['thermal']:.3f}"})
    n.update(_pure_bands())
    return n, tables


# the unmeasured formation links shown in the overview, keyed to
# formation_data.UNMEASURED as in make_overview.LINKS. One differs on v0_2:
# ML's "no Kr2* + F" is v0_2's own default (the band's 0), so v0_2's table
# tries Kannari's estimate and the top of the band instead
LINKS = tuple(
    x for link in MO.LINKS for x in (
        ((dict(k_kr2_f=3.0e-10), "kr2F.kann",
          r"\Krts{} + F \ra KrF(B) + Kr", r"$3.0{\times}10^{-10}$",
          "Kannari's estimate"),
         (dict(k_kr2_f=1.0e-9), "kr2F.hi", "", r"$1.0{\times}10^{-9}$",
          "the top of its band"))
        if link[1] == "kr2F" else (link,)))


def _light_link(name, v) -> str:
    if name.startswith("k_") and v == 0.0:
        return "0"
    return MO._light_link(name, v)


def slow_numbers() -> tuple:
    """(numbers, tables): the \\widef numbers and the finished table rows."""
    _ml()
    from ml import (data, eedf, fit, formation, formation_data, forward,
                    params, score, swarm)
    n, tables = {}, []
    # the posterior of the B/C block: its constants are v0's (the train
    # sets are at 300 K, where the B-C gap does not enter), and its Tahoma
    # gain is v0_2's
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
    # each unmeasured link at Tahoma, one at a time from v0_2's chain
    c = forward.with_rederived_k_e(params.light_now())
    chain_m = formation.LIGHT_CHAIN.with_(**formation_data.MEASURED)
    if chain_m.k_kr2_f != 0.0:
        raise RuntimeError("v0_2's chain no longer runs Kr2* + F at 0: "
                           "update LINKS")
    shown = [kw for kw, *_ in LINKS]
    for _, kw in formation_data.UNMEASURED:
        if kw == {"k_kr2_f": 0.0}:
            continue                     # v0_2's own value, the table's base
        if kw not in shown:
            raise KeyError(f"{kw} is in formation_data.UNMEASURED but not "
                           "in this overview's LINKS")
    base = formation.observables(c, chain_m)["eta"]
    n["link.base"] = f"{base:.3f}"
    minor = 0.0
    rows = []
    for kw, key, label, shown_v, what in LINKS:
        eta = formation.observables(c, chain_m.with_(**kw))["eta"]
        n[f"link.{key}"] = f"{eta:.3f}"
        n[f"link.{key}.d"] = MO._neg(f"{eta - base:+.3f}")
        if label is None:
            minor = max(minor, abs(eta - base))
            continue
        name = list(kw)[0]
        light = _light_link(name, getattr(chain_m, name)) if label else ""
        rows.append(f"{label} & {light} & {shown_v} & {what} & {eta:.3f} & "
                    f"{n[f'link.{key}.d']}\\\\")
    n["link.minor"] = f"{math.ceil(minor * 1000) / 1000:.3f}"
    tables.append("\\def\\wiLinks{" + "\n".join(rows) + "}")
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
                 "fast": lambda: _render(*_fast_parts()),
                 "slow": lambda: _render(*_slow_parts())}[a.check]()
        path = {"rates": OUT_RATES, "fast": OUT_FAST, "slow": OUT_SLOW}[a.check]
        if path.read_text() != fresh:
            print(f"{path.name} is stale: rerun make_overview_v0_2.py")
            return 1
        print(f"{path.name} is current")
        return 0
    OUT_RATES.write_text(rates_tex())
    print(f"wrote {OUT_RATES}")
    OUT_FAST.write_text(_render(*_fast_parts()))
    print(f"wrote {OUT_FAST}")
    if not a.fast:
        OUT_SLOW.write_text(_render(*_slow_parts()))
        print(f"wrote {OUT_SLOW}")
    return 0


def _fast_parts():
    n, tables = fast_numbers()
    return n, "fast", tables


def _slow_parts():
    n, tables = slow_numbers()
    return n, "slow", tables


if __name__ == "__main__":
    sys.exit(main())
