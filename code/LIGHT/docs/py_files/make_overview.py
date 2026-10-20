"""Numbers and tables for accuracy_overview.tex, LIGHT's accuracy today.

No number in the overview is typed by hand. This script writes three
files, each one \\ovdef{key}{value} per number or a finished table:

  overview_rates.tex  Appendix A, the reaction rates, read from the rows
                      LIGHT builds (network.build) and from species.py.
                      Takes a second.
  overview_fast.tex   the checks that take seconds: the photolysis decays
                      Xu 1993 printed, the ledger closure in every mode,
                      ML's scorecard of the kinetic datasets and of Tahoma
                      (a score, never a fit target), and the pure-gas
                      checks.
  overview_slow.tex   the checks that take minutes: the MCMC posterior of
                      the B/C constants, each unmeasured formation link at
                      Tahoma, the electron closure against the kinetic
                      EEDF loop, and the cross sections against swarm data.

test_docs.py regenerates the first two and compares them byte for byte.
It compares the third only with LIGHT_DOCS_SLOW=1 (about 5 minutes).
The accuracy numbers come from ML/ (beside LIGHT/), which runs LIGHT itself.

    python3 LIGHT/docs/make_overview.py           # all three, about 6 min
    python3 LIGHT/docs/make_overview.py --fast    # the first two, 30 s
    cd LIGHT/docs
    latexmk -pdf system_map_basic.tex             # the map page it includes
    latexmk -pdf accuracy_overview.tex
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0.py").is_file())
ML_DIR = LIGHT.parent / "ML"
for _p in (str(LIGHT), str(HERE), str(ML_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                            # noqa: E402

from make_system_map_basic import _k                          # noqa: E402
from v0 import __version__, xsec                              # noqa: E402
from v0.model import AMAGAT_CM3, Model, Params, run           # noqa: E402
from v0.network import BASIC_GROUPS, GROUPS, build            # noqa: E402
from v0.report import accounting                              # noqa: E402
from v0.species import SPECIES, e_form                        # noqa: E402

OUT_RATES = HERE / "overview_rates.tex"
OUT_FAST = HERE / "overview_fast.tex"
OUT_SLOW = HERE / "overview_slow.tex"

# ---------------------------------------------------------------------------
# Appendix A: the reaction rates
# ---------------------------------------------------------------------------
TEX = {"e": r"\el{}", "Kr": "Kr", "Kr*": r"\Krs{}", "Kr+": r"\Krp{}",
       "Kr2+": r"\Krtp{}", "Kr2*": r"\Krts{}", "Kr2F*": r"\KtF{}",
       "F2": r"\Ft{}", "F": "F", "F-": r"F\ensuremath{^-}",
       "KrF(B)": "KrF(B)", "KrF(C)": "KrF(C)"}
# the order the tables list the rows in: the electrons, Kr*, KrF formation,
# emission, B-C mixing, the trimer, quenching, Kr2F*. A row missing here is
# appended at the end, never dropped.
ORDER = ("G1", "C1a", "C1b", "C1c", "C2", "G6", "G2", "G3", "G3C", "G3X",
         "G4", "G4C", "BC", "CB", "G5", "G5C", "QE", "QEC", "QF2", "QF2C",
         "QKR", "QKRC", "C3", "Q2F",
         "AK8_Kr", "KR2_RAD", "NK7", "K2F_BR", "K2F_F", "PEN_KR2", "K15",
         "R11p", "KSTEP", "SUPER", "ATT", "DIS_a3", "DIS_A1", "VIB1", "VIB2",
         "VIB3", "VIB4", "K26", "K108",
         "K33", "K32", "K67", "K125", "K_2Fe", "BC_Kr", "CB_Kr", "URAD_C")
# Table A.3's labels for the electron-impact rows (the threshold is read
# from the row)
XSEC_LABELS = {"G1": "Kr, 5s manifold", "C1a": "Kr, 5p group",
               "C1b": "Kr, d+p group", "C1c": "Kr, higher group",
               "C2": "Kr, ionization", "ATT": "F$_2$, attachment",
               "DIS_a3": r"F$_2$, {\unsansmath$a\,^3\Pi_u$}",
               "DIS_A1": r"F$_2$, {\unsansmath$A\,^1\Pi_u$}",
               "VIB1": "F$_2$, $v=1$", "VIB2": "F$_2$, $v=2$",
               "VIB3": "F$_2$, $v=3$", "VIB4": "F$_2$, $v=4$"}


def ordered(rows) -> list:
    """Rows in ORDER. A C-state copy (bc group, 'G5_C') follows its B row."""
    rank = {rid: i for i, rid in enumerate(ORDER)}

    def key(r):
        if r.rid not in rank and r.rid.endswith("_C"):
            return (rank.get(r.rid[:-2], len(ORDER)), 1, r.rid)
        return (rank.get(r.rid, len(ORDER)), 0, r.rid)
    return sorted(rows, key=key)

# Short sources, one per row, condensed from each row's full citation in
# v0/network.py. Every number here must appear in that citation
# (test_docs.py checks it), so the table cannot drift from the code.
SOURCES = {
    "G1": "Biagi v7.1, the whole 5s manifold at 9.915 eV",
    "G2": "Kannari 1985 Table I, Penning block",
    "G3": "rate from Velazco, Kolts and Setser 1976 and Ueda 1984, KrF* "
          "yield 0.90 from Velazco, Kolts and Setser 1976, B share from Yu "
          "and Setser 1990",
    "G3C": "as G3, with Yu and Setser's k_B/k_C of 1.6",
    "G3X": "as G3, the KrF*-free rest, assigned to Kr + 2F",
    "G4": "Kannari 1985 Table I. The notes write about 4 ns",
    "G5": "ML joint fit. T_g law from Shui's calculation",
    "G6": "Biagi v7.1 momentum transfer",
    "C1a": "Biagi v7.1 5p group at 11.30 eV, prompt cascade to 5s",
    "C1b": "Biagi v7.1 d+p group at 11.998 eV, as C1a",
    "C1c": "Biagi v7.1 higher group at 12.75 eV, as C1a",
    "C2": "Biagi v7.1 ionization at 13.996 eV",
    "C3": "Xu, Gadomski and Setser 1993, 220 +- 20 ns",
    "G4C": "ML joint fit. Hakuta 1990 measured 60 +- 10 ns",
    "BC": "ML joint fit. Hakuta 1990's 1.2e-10 is held out",
    "CB": "detailed balance of BC, with the joint fit's K_eq",
    "G5C": "ML joint fit",
    "QF2": "ML joint fit. Xu 1993 printed 3.0 +- 0.3e-10 with fast mixing",
    "QF2C": "ML joint fit",
    "QKR": "ML joint fit, at its lower bound",
    "QKRC": "ML joint fit",
    "QE": "Trainor and Jacob's measured product, re-read with this "
          "network's lifetime",
    "QEC": "as QE",
    "Q2F": "Xu 1993 Table I, 1.2 +- 0.2e-10",
    # optional groups, as the basic network runs them
    "AK8_Kr": "Ueda 1984 Table V, 1.7 +- 0.2e-32. T_g law from Wieme and "
              "Lenaerts 1981",
    "KR2_RAD": "Morikawa 1989, 270 +- 5 ns, and five other measurements",
    "NK7": "half of Quigley and Hughes 1978, 2.1 +- 0.2e-10. The even split "
           "is Kannari's",
    "K2F_BR": "the other half of Quigley and Hughes 1978",
    "K2F_F": "Kannari 1985, an estimate. No measurement exists",
    "PEN_KR2": "Kannari 1985 Table I",
    "K15": "Kannari 1985, Boichenko 2000",
    "R11p": "Kannari 1985 Table I, on a T_e^-1/2 template",
    "KSTEP": "Chen R9 template, an estimate",
    "SUPER": "detailed balance of G1, statistical weight ratio 1/12",
    "ATT": "SIGLO F2 set",
    "DIS_a3": "SIGLO a3Pi_u, repulsive, so excitation dissociates",
    "DIS_A1": "SIGLO A1Pi_u, as DIS_a3",
    "VIB1": "SIGLO, v = 1", "VIB2": "SIGLO, v = 2", "VIB3": "SIGLO, v = 3",
    "VIB4": "SIGLO, v = 4",
    "K26": "Kannari 1985, Flannery turnover with Debye screening. No "
           "measurement with F- exists",
    "K108": "as K26",
}
# v0's rows where they differ from the basic network (the notes' list as
# written and REBUILD's copies)
SOURCES_V0 = {
    "G3": "Kannari 1985 Table I, all to B. The report it credits gives "
          "7.2e-10",
    "G5": "Xu, Gadomski and Setser 1993 Table II, 6.0 +- 2.0e-31",
    "AK8_Kr": "REBUILD's estimate, Kannari's Ar-partner rate given to Kr",
    "KR2_RAD": "Kannari 1985 Table I, 8.0e7 per second",
    "NK7": "Kannari 1985 Table I",
    "K2F_BR": "Kannari 1985 Table I",
    "K33": "Xu, Gadomski and Setser 1993, 4.0 +- 0.5e-12",
    "K32": "Xu, Gadomski and Setser 1993, 3.0 +- 0.3e-10",
    "K67": "Kushner, Hanson and Schneider 1989",
    "K125": "Kannari 1985 Table I",
    "K_2Fe": "Kannari 1985, estimated there",
    "BC_Kr": "Kannari's B-C mixing with a Kr partner, an estimate",
    "CB_Kr": "detailed balance of BC_Kr at 300 K",
    "URAD_C": "Kannari 1985 Table I, 1.33e7 per second",
}
SPECIES_NOTES = {
    "Kr": "ground state", "F2": "ground state",
    "F": "half the F2 bond, D0 = 1.60 eV",
    "Kr*": "Biagi's 5s threshold, the 5s[3/2]2 level",
    "Kr+": "Biagi's ionization threshold",
    "KrF(B)": "F + the 248-nm quantum + 0.35 eV of X-state repulsion",
    "KrF(C)": "80 cm-1 below B",
    "Kr2F*": "F + about 3.0 eV of 415-nm band + about 0.3 eV",
    "Kr2*": "the 147-nm excimer, an estimated well",
    "Kr2+": "IP(Kr) less D0(Kr2+) = 1.15 eV",
    "F-": "E_form(F) less the electron affinity, 3.40 eV",
}
TE_GRID = (0.5, 1.0, 1.5, 2.0, 3.0)                # eV, Table A.3
XSEC_ROWS = ("G1", "C1a", "C1b", "C1c", "C2", "ATT", "DIS_a3", "DIS_A1",
             "VIB1", "VIB2", "VIB3", "VIB4")


def _tex_text(s: str) -> str:
    """A short source as LaTeX text."""
    s = s.replace("+-", r"$\pm$")
    s = re.sub(r"(\d)e(-?\d+)", lambda m: m.group(1) + r"$\times10^{"
               + str(int(m.group(2))) + "}$", s)
    s = s.replace("k_B/k_C", "$k_B/k_C$").replace("T_g", "$T_g$")
    s = s.replace("T_e^-1/2", "$T_e^{-1/2}$").replace("K_eq", r"$K_\mathrm{eq}$")
    # sansmath has no upright-sans uppercase Greek (LIGHT's Known-bad list)
    s = s.replace("a3Pi_u", r"{\unsansmath$a\,^3\Pi_u$}")
    s = s.replace("A1Pi_u", r"{\unsansmath$A\,^1\Pi_u$}")
    s = s.replace("F2", r"F$_2$").replace("F-", r"F$^-$")
    s = s.replace("KrF*", r"KrF$^*$").replace("cm-1", r"cm$^{-1}$")
    s = s.replace("Kr2+", r"Kr$_2^+$").replace("D0", r"$D_0$")
    s = s.replace("E_form", r"$E_\mathrm{form}$")
    # escape the underscores left in text (row ids), not those in math
    parts = s.split("$")
    parts[::2] = [t.replace("_", r"\_") for t in parts[::2]]
    return "$".join(parts)


def _side(species) -> str:
    out, seen = [], []
    for s in species:
        if s in seen:
            continue
        seen.append(s)
        n = species.count(s)
        out.append((str(n) if n > 1 else "") + TEX[s])
    return " + ".join(out)


def reaction(r) -> str:
    rhs = _side(r.products) if r.products else ""
    if r.is_rad:
        rhs = (rhs + r" + \hnu") if rhs else r"\hnu"
    if r.rate[0] == "elastic":
        return _side(r.reactants) + ", elastic"
    return _side(r.reactants) + r" \ra " + rhs


def _order(r) -> int:
    return len(r.heavy) + r.n_electrons_in


def _unit(r) -> str:
    return {1: r"s$^{-1}$", 2: r"cm$^3$s$^{-1}$", 3: r"cm$^6$s$^{-1}$"}[
        _order(r)]


def rate_text(r, p: Params) -> str:
    """The rate coefficient at 300 K as printed in the table."""
    kind = r.rate[0]
    if kind == "const":
        return f"${_k(r.rate[1])}$ {_unit(r)}"
    if kind == "inv_tau":
        return f"${_k(1.0 / getattr(p, r.rate[1]))}$ {_unit(r)}"
    if kind == "balance":
        k_fwd, k_ref = r.rate[1], r.rate[2]
        return f"${_k(k_fwd / k_ref)}$ {_unit(r)}"
    if kind == "xsec":
        return r"$\langle\sigma v\rangle$, Table A.3"
    if kind == "elastic":
        return r"energy loss $(2m/M)\langle\sigma_m v\varepsilon\rangle$"
    if kind == "stepwise":
        return (r"$1.56{\times}10^{-7}\,T_e^{0.71}\times$ $e^{-"
                + f"{r.rate[1]:g}" + r"/T_e}$ " + _unit(r))
    if kind == "dr_sqrt":
        return f"${_k(r.rate[1])}\\,T_e^{{-1/2}}$ {_unit(r)}"
    if kind == "flannery":
        n_am = 101325.0 / (1.380649e-23 * 300.0) * 1e-6 / AMAGAT_CM3
        k = r.rate[1] * n_am / (n_am * n_am + 1.7) * r.rate[2]
        return f"${_k(k)}$ {_unit(r)} at 1\\,atm, 300\\,K, unscreened"
    if kind == "superelastic":
        return r"from G1 by detailed balance"
    raise ValueError(f"{r.rid}: rate kind {kind!r}")


def tdep_text(r) -> str:
    kind = r.rate[0]
    parts = []
    if kind in ("xsec", "stepwise", "dr_sqrt", "superelastic"):
        parts.append("$T_e$")
    if kind == "elastic":
        parts.append("$T_e$, $T_g$")
    if kind == "balance":
        parts.append(r"$1/K_\mathrm{eq}(T_g)$")
    if kind == "flannery":
        parts.append("$N$, screening")
    if r.tg_pow:
        parts.append(f"$(T_g/300)^{{{r.tg_pow:.2f}}}$")
    return ", ".join(parts) if parts else ""


def _rate_rows(rows, sources, p) -> str:
    out = []
    for r in rows:
        src = sources.get(r.rid)
        if src is None and r.rid.endswith("_C"):
            src = f"as {r.rid[:-2]}, applied to KrF(C)"
        if src is None:
            raise KeyError(f"no short source for row {r.rid}")
        rid = r.rid.replace("_", "\\_")
        out.append(f"{rid} & {reaction(r)} & {rate_text(r, p)} & "
                   f"{tdep_text(r)} & {_tex_text(src)}\\\\")
    return "\n".join(out)


def v0_differences() -> list:
    """v0's rows (every group) that the basic network lacks or runs with
    another value."""
    basic = {r.rid: r for r in build(("all",), "basic")}
    out = []
    for r in build(("all",), "v0"):
        b = basic.get(r.rid)
        if b is None or b.rate != r.rate or b.tg_pow != r.tg_pow:
            out.append(r)
    return out


def _fac(f: float) -> str:
    """A factor: two significant figures, powers of ten below 0.01."""
    return f"{f:#.2g}" if f >= 0.01 else f"${_k(f)}$"


def xsec_table(p: Params, te_mid: float, r_mid: float) -> str:
    """Table A.3: Maxwellian rate coefficients of the electron-impact rows,
    the captured energy of attachment, and the two-temperature factor at
    the reference point's mid-pulse state."""
    rows = {r.rid: r for r in build(("all",), "basic")}
    lines = []
    for rid in XSEC_ROWS:
        r = rows[rid]
        gas, kind, prod = r.rate[1:4]
        tab = xsec.rate_table(gas, kind, prod)
        cells = [f"${_k(float(tab(te)))}$" for te in TE_GRID]
        if xsec.in_tail(gas, kind, prod):
            two = xsec.two_t_table(gas, kind, prod)
            f = float(two(te_mid, r_mid * te_mid)) / float(tab(te_mid))
            cells.append(_fac(f))
        else:
            cells.append("1")
        tid = rid.replace("_", "\\_")
        label = f"{XSEC_LABELS[rid]}, {r.e_cost_eV:g}\\,eV" \
            if r.e_cost_eV is not None else XSEC_LABELS[rid]
        lines.append(f"{tid} & {label} & " + " & ".join(cells) + "\\\\")
    cap = xsec.capture_energy_table("F2", "attachment")
    cells = [f"{float(cap(te)):.2f}" for te in TE_GRID]
    lines.append(r"\midrule" "\n" r"ATT & energy an attached \el{} takes, eV & "
                 + " & ".join(cells) + r" & \\")
    cells = [f"{1.5 * te:.2f}" for te in TE_GRID]
    lines.append(r" & $\tfrac32 T_e$, for comparison & "
                 + " & ".join(cells) + r" & \\")
    return "\n".join(lines)


def species_table() -> str:
    names = ("Kr", "F2", "F", "Kr*", "Kr+", "KrF(B)", "KrF(C)", "Kr2F*",
             "Kr2*", "Kr2+", "F-")
    lines = []
    for n in names:
        s = SPECIES[n]
        lines.append(f"{TEX[n]} & {_neg(f'{e_form(n):.6g}')} & {s.cv:g} & "
                     f"{_tex_text(SPECIES_NOTES[n])}\\\\")
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
        "% LIGHT reaction-rate tables. GENERATED by make_overview.py from",
        f"% LIGHT v0 ({__version__}). Do not edit. Rerun the script instead.",
        r"\providecommand\ovdef[2]{\expandafter\def\csname ov@#1\endcsname{#2}}",
        r"\ovdef{rates.teMid}{" + f"{te_mid:.2f}" + "}",
        r"\ovdef{rates.rMid}{" + f"{r_mid:.2f}" + "}",
        r"\ovdef{rates.nBasic}{" + str(len(basic)) + "}",
        r"\ovdef{rates.nGroups}{" + str(len(groups)) + "}",
        r"\ovdef{rates.groups}{" + ", ".join(BASIC_GROUPS) + "}",
        r"\ovdef{rates.v0groups}{" + ", ".join(GROUPS) + "}",
        r"\def\ratesBasic{" + _rate_rows(ordered(basic), SOURCES, p) + "}",
        r"\def\ratesGroups{" + _rate_rows(ordered(groups), SOURCES, p)
        + "}",
        r"\def\ratesVzero{" + _rate_rows(ordered(v0_differences()),
                                         {**SOURCES, **SOURCES_V0}, p) + "}",
        r"\def\ratesXsec{" + xsec_table(p, te_mid, r_mid) + "}",
        r"\def\ratesTe{" + " & ".join(f"{t:g}" for t in TE_GRID) + "}",
        r"\def\speciesTable{" + species_table() + "}",
    ]
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# formatting helpers for the accuracy numbers
# ---------------------------------------------------------------------------
def _neg(s: str) -> str:
    """A proper minus sign for LaTeX text."""
    return s.replace("-", r"\ensuremath{-}")


def _x(model: float, meas: float) -> str:
    return f"{model / meas:.2f}"


def _pct_signed(v: float) -> str:
    return _neg(f"{v:+.1f}")


def _render(n: dict, what: str, tables=()) -> str:
    lines = [
        f"% LIGHT accuracy numbers ({what}). GENERATED by make_overview.py",
        f"% from LIGHT v0 ({__version__}) and ML/. Do not edit. Rerun the "
        "script instead.",
        r"\providecommand\ovdef[2]{\expandafter\def\csname ov@#1\endcsname{#2}}",
    ]
    lines += [f"\\ovdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# the fast accuracy numbers
# ---------------------------------------------------------------------------
# Xu, Gadomski & Setser, JCP 99, 2591 (1993) p.2594 print the KrF(B,C)
# decay constant k' after photolysis in Kr with 1 Torr F2 at 300 K
XU_1993 = ((50.0, 5.0e7), (100.0, 6.2e7))
# Torr of Kr and F2 where the eq. 5 agreement is checked (test_model.py's)
EQ5_POINTS = ((50, 1), (100, 1), (500, 1), (100, 2), (300, 2), (50, 6))


def _krf_decay(kr_torr, f2_torr, net="basic"):
    """Seed KrF(B) in pure Kr/F2 at 300 K (no power, no electrons) and fit
    the slow decay of KrF(B) + KrF(C), as in the photolysis experiments
    (the same run as test_model.py's check)."""
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


def fast_numbers() -> dict:
    from ml import data, forward, formation, params, pure, pure_data, score
    n = {}
    # the photolysis decays against Xu 1993's printed constants
    for kr, meas in XU_1993:
        n[f"xu.{kr:.0f}"] = _x(_krf_decay(kr, 1.0), meas)
        n[f"xu.{kr:.0f}.v0"] = _x(_krf_decay(kr, 1.0, "v0"), meas)
    for kr, meas in XU_1993:
        n[f"xu.{kr:.0f}.meas"] = f"${_k(meas)}$"
    # the basic network's B/C decay against Gadomski's analytic eq. 5
    from v0 import two_level as L
    c_light = L.Consts.from_light()
    dev = 0.0
    for kr, f2 in EQ5_POINTS:
        slow = L.eigen(c_light, kr=kr * L.TORR_300K, f2=f2 * L.TORR_300K)[1]
        dev = max(dev, abs(_krf_decay(kr, f2) / slow - 1.0))
    mant, ex = f"{dev:.1e}".split("e")
    n["eq5.dev"] = f"${mant}{{\\times}}10^{{{int(ex)}}}$"
    # the ledger closure in every mode
    base = Params()
    runs = [base.replace(groups=(g,)) for g in BASIC_GROUPS]
    runs += [base, base.replace(groups=("all",)),
             base.replace(groups=("all",), ne_mode="free"),
             base.replace(network="v0"),
             base.replace(network="v0", groups=("all",))]
    worst = max(abs(accounting(run(p))["closure"]) for p in runs)
    n["closure.max"] = f"$10^{{{math.ceil(math.log10(worst))}}}$"
    n["closure.runs"] = str(len(runs))
    # ML's scorecard for LIGHT's own constants (k_e re-derived as ML does)
    c = forward.with_rederived_k_e(params.light_now())
    for ds in data.all_datasets():
        if ds.kind in data.DEVICE_KINDS:
            continue
        v = score.kinetic_score(ds, c)
        n[f"sc.{ds.name}"] = (_pct_signed(v) if len(ds.values) == 1
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
    # Tahoma, a score only
    o = formation.observables(c)
    meas = {d.name: d for d in data.by_role("device")}
    g0, eta, isat = (meas[k].values[0] for k in ("tahoma_g0",
                                                  "tahoma_formation",
                                                  "tahoma_isat"))
    n.update({
        "tah.g0.meas": f"{g0:.2f}", "tah.g0": f"{o['g0']:.2f}",
        "tah.g0.x": _x(o["g0"], g0),
        "tah.isat.meas": f"{isat:.1f}", "tah.isat": f"{o['isat']:.2f}",
        "tah.isat.err": f"{meas['tahoma_isat'].sigma * isat:.1f}",
        "tah.isat.x": _x(o["isat"], isat),
        "tah.eta.meas": f"{eta:.2f}", "tah.eta": f"{o['eta']:.3f}",
        "tah.eta.x": _x(o["eta"], eta),
        "tah.form": f"{o['formation']:.3f}",
        "tah.pin": f"{forward.TAHOMA_COND['power_W_cm3'] / 1e3:.0f}",
        "tah.torr": f"{forward.TAHOMA_COND['pressure_torr']:.0f}",
        "tah.f2torr": f"{forward.TAHOMA_COND['f2_torr']:.1f}",
    })
    # the pure-gas checks (the F2-free half of the chain)
    for g in ("Kr", "Ar"):
        x = pure.supply(g)
        n[f"pure.{g}.x"] = _x(x["per_eV_sf"], x["per_eV_meas"])
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
              "pure.vuv.x": _x(v["vuv"], vk.value),
              "pure.vuv.ideal": f"{pure.ideal_vuv('Kr'):.2f}",
              "pure.vuv.nothermal": f"{pure.kr_vuv(3.3e2, thermal=False)['vuv']:.2f}",
              "pure.thermal": f"{v['thermal']:.3f}"})
    return n


# ---------------------------------------------------------------------------
# the slow accuracy numbers
# ---------------------------------------------------------------------------
# the unmeasured formation links shown in the overview, keyed to
# formation_data.UNMEASURED (the script fails if ML's list changes)
# (the change, key, table label or None for the minor Kr2* links, the value
# tried as printed, what that value is)
LINKS = (
    (dict(gamma_H=0.70), "yield.lo", r"KrF$^*$ yield of \Krs{} + \Ft",
     "0.70", "the measured 0.90 less its stated 0.20"),
    (dict(gamma_H=1.0), "yield.hi", "", "1.0", "the published assignment"),
    (dict(gamma_H=0.5), "yield.arf", "", "0.5", r"as for Ar$^*$ + F$_2$"),
    (dict(f_C_H=0.0), "nascentC", r"share of that KrF formed in C", "0",
     "if nascent C relaxes to B at 1 atm"),
    (dict(y_ion=0.64), "ionY", r"KrF$^*$ yield of ion-ion recombination",
     "0.64", r"as for Cl$^-$"),
    (dict(f_C_ion=0.38), "ionC", r"share of that KrF formed in C", "0.38",
     r"as for SF$_6^-$"),
    (dict(ion_scale=0.38), "ionRate", r"ion-ion recombination rate",
     r"\x{0.38}", "Avco's model input"),
    (dict(k_kr2_f=0.0), "kr2F", r"\Krts{} + F \ra KrF(B) + Kr", "0",
     "no measurement exists"),
    (dict(kr2_krf=0.0), "kr2split.lo", None, "", ""),
    (dict(kr2_krf=1.0), "kr2split.hi", None, "", ""),
    (dict(k3=4.4e-32), "k3", None, "", ""),
)


def _light_link(name, v) -> str:
    """LIGHT's own value of a formation link, as the links table prints it."""
    if name == "ion_scale":
        return "Flannery"
    if name.startswith("k_"):
        return f"${_k(v)}$"
    if name in ("y_ion", "f_C_ion"):
        return f"{v:g}"
    return f"{v:.2f}"


def slow_numbers() -> tuple:
    """(numbers, tables): the \\ovdef numbers and the finished table rows."""
    from ml import (data, eedf, fit, formation, formation_data, forward,
                    params, score, swarm)
    n, tables = {}, []
    # the MCMC posterior of the B/C block on the train set
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
    n["fit.names"] = ", ".join(d.name.replace("_", r"\_") for d in train)
    n_draws = 40
    lo, mid, hi = score.device_posterior(chain, n_draws=n_draws)
    n.update({"post.g0.lo": f"{lo:.2f}", "post.g0.mid": f"{mid:.2f}",
              "post.g0.hi": f"{hi:.2f}", "post.n": str(n_draws),
              "fit.n": str(len(params.NAMES))})
    # each unmeasured link at Tahoma, one at a time from LIGHT's chain
    c = forward.with_rederived_k_e(params.light_now())
    have = [kw for _, kw in formation_data.UNMEASURED]
    for kw, key, _, shown, _ in LINKS:
        if shown and not shown.startswith("\\") and \
                float(shown) != list(kw.values())[0]:
            raise ValueError(f"{key}: printed {shown} is not {kw}")
        if kw not in have:
            raise KeyError(f"{kw} is no longer in formation_data.UNMEASURED")
    if len(have) != len(LINKS):
        raise KeyError("formation_data.UNMEASURED has links this overview "
                       "does not show; add them to LINKS")
    chain_m = formation.LIGHT_CHAIN.with_(**formation_data.MEASURED)
    base = formation.observables(c, chain_m)["eta"]
    n["link.base"] = f"{base:.3f}"
    minor = 0.0
    rows = []
    for kw, key, label, shown, what in LINKS:
        eta = formation.observables(c, chain_m.with_(**kw))["eta"]
        n[f"link.{key}"] = f"{eta:.3f}"
        n[f"link.{key}.d"] = _neg(f"{eta - base:+.3f}")
        if label is None:
            minor = max(minor, abs(eta - base))
            continue
        name = list(kw)[0]
        light = _light_link(name, getattr(chain_m, name)) if label else ""
        rows.append(f"{label} & {light} & {shown} & {what} & {eta:.3f} & "
                    f"{n[f'link.{key}.d']}\\\\")
    n["link.minor"] = f"{math.ceil(minor * 1000) / 1000:.3f}"
    tables.append("\\def\\linksTable{" + "\n".join(rows) + "}")
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
    n["el.eatt.te"] = f"{2.0 / 3.0 * b['mean']:.2f}"
    n["el.eatt.32"] = f"{b['mean']:.2f}"
    # the cross sections against swarm data
    by = {ds.name: ds for ds in swarm.load()}
    r = swarm.predict(by["kr_w_nakamura"]) / by["kr_w_nakamura"].y
    n["sw.krw.lo"], n["sw.krw.hi"] = f"{r.min():.2f}", f"{r.max():.2f}"
    s_low, err = swarm.fit_low([by[k] for k in swarm.TRAIN_F2])
    n["sw.slow"] = f"{s_low:.2f}"
    n["sw.slow.err"] = f"{s_low * err:.2f}"
    n["sw.slow.sigma"] = f"{abs(math.log(s_low)) / err:.1f}"
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
    a = ap.parse_args(argv)
    OUT_RATES.write_text(rates_tex())
    print(f"wrote {OUT_RATES}")
    OUT_FAST.write_text(_render(fast_numbers(), "fast"))
    print(f"wrote {OUT_FAST}")
    if not a.fast:
        n, tables = slow_numbers()
        OUT_SLOW.write_text(_render(n, "slow", tables))
        print(f"wrote {OUT_SLOW}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
