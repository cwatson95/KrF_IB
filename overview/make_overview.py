#!/usr/bin/env python3
"""Write the generated tables of KrF_IB/overview.pdf from the models' code.

For LIGHT v0_5 (Kr/F2, KrF) and LIGHT_Kr v0_5 (pure Kr, code/Kr_LIGHT/) it
reads each model's default network (network.build(), the rows a default
run integrates), its bands (network.BANDS) and its species' formation
energies (species.e_form), and evaluates every row's rate constant with the
model's own Model._k at one reference state:

  Tg = 300 K, N = the 1-atm fill, n_e = 1e14 cm-3 (the ion-ion rows' Debye
  screening), electron rows on a Maxwellian at Te = 1 eV and at 3 eV,
  every band at its default.

The models' default closure is the kinetic EEDF, solved during a run; the
Maxwellian values here are a reference, not what a run uses. Nothing is
typed by hand: rerun this file after a change to either network.

    /usr/bin/python3 overview/make_overview.py      # from KrF_IB/
    cd overview && latexmk -pdf overview.tex && cp build/overview.pdf ../

Writes overview/generated/{rates,bands,species}_{light,kr}.tex and
counts.tex.
"""

from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code"
OUT = HERE / "generated"

TE_REF = (1.0, 3.0)        # eV, Maxwellian
TG_REF = 300.0             # K
NE_REF = 1e14              # cm-3, for the ion-ion rows' screening only
SOURCE_CHARS = 170         # the source column's excerpt

# (tag, directory under code/, package, label)
MODELS = (("light", "LIGHT", "v0_5", r"LIGHT v0\_5"),
          ("kr", "Kr_LIGHT", "v0_5", r"LIGHT\_Kr v0\_5"))


# ------------------------------------------------------------ LaTeX helpers
_TEX = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
        "<": r"\textless{}", ">": r"\textgreater{}", "—": "---"}


def tex(s: str) -> str:
    return "".join(_TEX.get(ch, ch) for ch in str(s))


def sci(x: float, sig: int = 3) -> str:
    """x in math mode, a x 10^b; 0 as 0."""
    if x is None:
        return ""
    if x == 0:
        return "$0$"
    m, e = f"{x:.{sig - 1}e}".split("e")
    e = int(e)
    m = f"{float(m):.{sig}g}"
    if e == 0:
        return f"${m}$"
    if m == "1":
        return f"$10^{{{e}}}$"
    return f"${m}\\times10^{{{e}}}$"


def num(x: float, sig: int = 3) -> str:
    """A plain number where it reads better than sci()."""
    if x != 0 and (abs(x) < 1e-2 or abs(x) >= 1e4):
        return sci(x, sig)
    return f"${float(f'{x:.{sig}g}'):g}$"


_SUFFIX = {"(0u+)": r"(0$_u^+$)"}


def species(s: str) -> str:
    """'Kr2F*' -> Kr$_2$F$^*$, 'Kr*(1s5)' -> Kr$^*$(1s$_5$), 'e' -> e$^-$."""
    if s == "e":
        return r"e$^-$"
    m = re.match(r"^([^()]*)(\(.*\))?$", s)
    core, suffix = m.group(1), m.group(2) or ""
    body = re.match(r"^((?:[A-Z][a-z]?\d*)+)([*+\-]*)$", core)
    if not body:
        return tex(s)
    out = ""
    for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", body.group(1)):
        out += el + (f"$_{cnt}$" if cnt else "")
    charge = body.group(2)
    if charge:
        if out.endswith("$"):          # merge with the subscript: Kr$_2^+$
            out = out[:-1] + f"^{{{charge}}}$"
        else:
            out += f"$^{{{charge}}}$"
    if suffix in _SUFFIX:
        suffix = _SUFFIX[suffix]
    else:
        suffix = re.sub(r"\((1s|2p)(\d+)\)", r"(\1$_{\2}$)", suffix)
    return out + suffix


def side(names) -> str:
    """Reactants or products, repeated species counted: Kr + Kr -> 2Kr."""
    seen = []
    for s in names:
        if s not in seen:
            seen.append(s)
    return " + ".join((f"{names.count(s)}" if names.count(s) > 1 else "")
                      + species(s) for s in seen)


def excerpt(cite: str, n: int = SOURCE_CHARS) -> str:
    cite = " ".join(cite.split())
    if len(cite) <= n:
        return tex(cite)
    cut = cite[:n].rsplit(" ", 1)[0].rstrip(",;:(")
    return tex(cut) + r"\,\ldots"


# ------------------------------------------------------------ the models
def load(root: str, pkg: str):
    """Import code/<root>/<pkg> as the launchers do (root on sys.path).
    Both models' packages are named v0_5, so the previous one's modules
    are dropped first."""
    for name in [m for m in sys.modules if m == pkg or m.startswith(pkg + ".")]:
        del sys.modules[name]
    path = str(CODE / root)
    sys.path.insert(0, path)
    try:
        mods = {m: importlib.import_module(f"{pkg}.{m}")
                for m in ("model", "network", "species")}
    finally:
        sys.path.remove(path)
    return mods


def row_bands(net, r) -> list:
    if hasattr(net, "_row_bands"):          # LIGHT_Kr: every band a row reads
        return sorted(net._row_bands(r))
    return sorted(net.bands_read([r]))


def law(r, model) -> str:
    """The row's rate law in words and symbols, from its rate spec."""
    spec = r.rate
    kind = spec[0]

    def inner(s):
        if s[0] == "const":
            return sci(s[1])
        if s[0] == "param":
            f = s[2] if len(s) > 2 else 1.0
            return (f"{num(f)}$\\,\\times$\\," if f != 1 else "") + \
                f"\\textsf{{{tex(s[1])}}}"
        return tex(str(s))

    if kind == "const":
        out = sci(spec[1])
    elif kind == "inv_tau":
        f = spec[2] if len(spec) > 2 else 1.0
        tau = float(getattr(model.p, spec[1]))
        out = (f"{num(f)}$\\,/\\tau$, $\\tau$ = {num(tau * 1e9)}\\,ns "
               f"(\\textsf{{{tex(spec[1])}}})")
    elif kind == "param":
        out = "band " + inner(spec)
    elif kind == "xsec":
        out = (f"Maxwellian $\\langle\\sigma v\\rangle$, {tex(spec[1])} "
               f"{tex(spec[2])} \\texttt{{{tex(spec[3] or '')}}}")
    elif kind == "elastic":
        out = (f"elastic energy loss on the {tex(spec[1])} momentum-transfer "
               "cross section; no density change")
    elif kind == "stepwise":
        out = (f"$1.56\\times10^{{-7}}T_e^{{0.71}}e^{{-E/T_e}}$, "
               f"$E$ = {num(spec[1])}\\,eV")
    elif kind == "dr_sqrt":
        out = f"{sci(spec[1])}$\\,/\\sqrt{{T_e}}$"
    elif kind == "te_pow":
        k0, te0, n = spec[1:]
        ref = ("300\\,K" if abs(te0 - 0.025852) < 1e-4
               else f"{num(te0)}\\,eV")
        out = f"{sci(k0)}$\\,(T_e/{ref})^{{{n:g}}}$"
    elif kind == "flannery":
        out = (f"Flannery: $aN/(N^2+1.7)$, $a$ = {sci(spec[1])}, "
               f"$\\times${num(spec[2])}, $N$ in amagat")
    elif kind == "ionion":
        k3, kl, n3 = spec[1:]
        out = (f"harmonic mean of $k_3N(T_g/300)^{{{n3:g}}}$ and $k_L/N$, "
               f"$k_3$ = {sci(k3)}, $k_L$ = {sci(kl)} (1 amagat), "
               "$\\times$ Debye screening")
    elif kind == "superelastic":
        out = (f"Klein--Rosseland reverse of \\texttt{{{tex(spec[2])}}}, "
               f"$g$ ratio {num(spec[3])}")
    elif kind == "balance":
        k_f, k_ref, up, lo = spec[1:]
        out = (f"detailed balance: $k_f/K(T_g)$, $k_f$ = {sci(k_f)}, "
               f"$K_{{300}}$ = {num(k_ref)} over the "
               f"{species(up)}--{species(lo)} gap")
    elif kind == "te_fit":
        a = (f"\\textsf{{{tex(spec[1])}}}" if isinstance(spec[1], str)
             else sci(spec[1]))
        out = f"{a}$\\,T_e^{{{spec[2]:.3g}}}e^{{-{spec[3]:.3g}/T_e}}$"
        if len(spec) > 4:
            out = (f"reverse by detailed balance ($g$ {num(spec[5])}) of "
                   + out)
    elif kind == "pop_xsec":
        if len(spec) > 3 and spec[3] in ("pair", "pair_rev"):
            pair = f"{spec[2]}$\\to${spec[4]}"
            out = ("Jung et al.\\ measured $\\sigma$, " + pair
                   if spec[3] == "pair"
                   else "Klein--Rosseland reverse of the measured " + pair)
        elif len(spec) > 2:
            out = (f"Klein--Rosseland reverse of the lumped measured "
                   f"$\\sigma$, $g$ ratio {num(spec[3])}")
        else:
            out = ("Maxwellian rate of the lumped measured 5s$\\to$5p "
                   "$\\sigma$ (Jung et al.)")
    elif kind == "assoc":
        out = (f"Hornbeck--Molnar branch of \\texttt{{{tex(spec[3])}}}'s "
               "$\\langle\\sigma v\\rangle$")
    elif kind == "vs_ion":
        out = (f"Vriens--Smeets $\\sigma$, $E_n$ = {num(spec[1])}\\,eV, "
               f"$\\times$\\,\\textsf{{{tex(spec[2])}}}")
    elif kind == "rung_reverse":
        out = (f"reverse of {inner(spec[1])} by detailed balance at $T_g$ "
               f"(5p levels at $T_g$, $g_{{5s}}$ = {spec[2]})")
    elif kind == "eq_reverse":
        out = (f"reverse of {inner(spec[1])}"
               f"$(T_g/300)^{{{spec[2]:g}}}$ by detailed balance on "
               f"{species(spec[3])}'s well")
        if len(spec) > 5:
            out += f" (\\textsf{{{tex(spec[5])}}})"
    elif kind == "level_reverse":
        out = (f"reverse of {inner(spec[1])}"
               f"$(T_g/300)^{{{spec[2]:g}}}$ by detailed balance, "
               f"$g$ ratio {num(spec[3])}")
    else:
        out = tex(str(spec))
    if r.tg_pow:
        out += f" $\\times(T_g/300)^{{{r.tg_pow:.3g}}}$"
    for name, how in r.frac:
        out += (f" $\\times$\\,\\textsf{{{tex(name)}}}" if how == "x"
                else f" $\\times(1-$\\textsf{{{tex(name)}}}$)$")
    tg_fac = getattr(r, "tg_fac", ())
    if tg_fac:
        out += f" $\\times F(T_g)^x$, $x$ = \\textsf{{{tex(tg_fac[0])}}}"
    level = getattr(r, "rung_level", None)
    if level:
        out += f"; per atom in {level}$^\\dagger$"
    return out


def unit(r) -> str:
    if r.rate[0] == "assoc":        # first order in Kr (network.py)
        return r"cm$^3$/s"
    order = len(r.reactants)
    return {1: r"s$^{-1}$", 2: r"cm$^3$/s", 3: r"cm$^6$/s"}.get(order, "")


_LEVEL = re.compile(r"^(1s\d|2p\d+|low|high|\d[a-z]\d+|\d[a-z]|\d)$")


def family(rid: str) -> str:
    """A row id less its level tags: K76_1s5_2p2 -> K76, KSTEP_S_1s5 ->
    KSTEP_S, HM_3s1_3 -> HM."""
    parts = rid.split("_")
    while len(parts) > 1 and _LEVEL.match(parts[-1].lower()):
        parts.pop()
    return "_".join(parts)


def rates_table(tag, label, m) -> tuple:
    model, net, spc = m["model"], m["network"], m["species"]
    M = model.Model(model.Params(eedf="maxwellian"))
    n = M.initial_state()[:len(M.species)] * M.N0
    lines = []
    last = None
    zero = 0
    for c in M.compiled:
        r = c.row
        if c.kind == "elastic":
            k = (None, None)
        elif r.n_electrons_in:
            k = tuple(M._k(c, te, TG_REF, n, NE_REF) for te in TE_REF)
        else:
            k = (M._k(c, TE_REF[0], TG_REF, n, NE_REF), None)
        if k[0] == 0:
            zero += 1
        fam = family(r.rid)
        if last is not None and fam != last:
            lines.append(r"\addlinespace[2pt]")
        last = fam
        bands = ", ".join(f"\\textsf{{{tex(b)}}}" for b in row_bands(net, r))
        k1 = sci(k[0]) if k[0] is not None else "--"
        k3 = (sci(k[1]) if k[1] is not None
              else ("--" if c.kind == "elastic" else r"\textcolor{ink2}{same}"))
        lines.append(" & ".join([
            f"\\texttt{{{tex(r.rid)}}}",
            f"{side(r.reactants)} $\\to$ {side(r.products)}",
            law(r, M), k1, k3, unit(r), bands,
            f"\\textcolor{{ink2}}{{{excerpt(r.cite)}}}"]) + r" \\")
    body = "\n".join(lines)
    (OUT / f"rates_{tag}.tex").write_text(
        f"% generated by make_overview.py from code/{label}: do not edit\n"
        + body + "\n")
    return len(M.compiled), zero, M


def optional_table(tag, m, M) -> int:
    """LIGHT's optional rows a default run leaves out (groups=())."""
    net = m["network"]
    allowed = getattr(net, "BASIC_GROUPS", ())
    rows = [r for r in net.optional_rows() if r.group in allowed]
    lines = []
    for r in rows:
        lines.append(" & ".join([
            f"\\texttt{{{tex(r.rid)}}}", tex(r.group),
            f"{side(r.reactants)} $\\to$ {side(r.products)}",
            law(r, M), f"\\textcolor{{ink2}}{{{excerpt(r.cite)}}}"])
            + r" \\")
    (OUT / f"optional_{tag}.tex").write_text(
        "% generated by make_overview.py: do not edit\n"
        + "\n".join(lines) + "\n")
    return len(rows)


def bands_table(tag, m, M) -> int:
    net = m["network"]
    readers = {}
    for c in M.compiled:
        for b in row_bands(net, c.row):
            readers.setdefault(b, []).append(c.row.rid)
    lines = []
    for b in net.BANDS.values():
        who = readers.get(b.name, [])
        fams = []
        for rid in who:
            f = family(rid)
            if f not in fams:
                fams.append(f)
        rows = (f"\\texttt{{{tex(', '.join(fams))}}}" if fams
                else r"\textcolor{ink2}{none (not a row rate)}")
        lines.append(" & ".join([
            f"\\textsf{{{tex(b.name)}}}", tex(b.what), num(b.default),
            f"{num(b.lo)} to {num(b.hi)}", tex(b.unit), rows]) + r" \\")
    (OUT / f"bands_{tag}.tex").write_text(
        "% generated by make_overview.py: do not edit\n"
        + "\n".join(lines) + "\n")
    return len(net.BANDS)


def species_table(tag, m, M) -> int:
    spc = m["species"]
    names = sorted(M.species, key=lambda s: -spc.e_form(s))
    lines = [f"{species(s)} & {spc.e_form(s):.3f} \\\\" for s in names]
    (OUT / f"species_{tag}.tex").write_text(
        "% generated by make_overview.py: do not edit\n"
        + "\n".join(lines) + "\n")
    return len(names)


def main() -> int:
    OUT.mkdir(exist_ok=True)
    counts = []
    for tag, root, pkg, label in MODELS:
        m = load(root, pkg)
        nrows, zero, M = rates_table(tag, f"{root}/{pkg}", m)
        nopt = optional_table(tag, m, M)
        nb = bands_table(tag, m, M)
        ns = species_table(tag, m, M)
        name = {"light": "Light", "kr": "Kr"}[tag]
        counts += [f"\\newcommand\\nrows{name}{{{nrows}}}",
                   f"\\newcommand\\nzero{name}{{{zero}}}",
                   f"\\newcommand\\nopt{name}{{{nopt}}}",
                   f"\\newcommand\\nbands{name}{{{nb}}}",
                   f"\\newcommand\\nspecies{name}{{{ns}}}"]
        print(f"{label}: {nrows} rows ({zero} at 0 by default), {nopt} "
              f"optional, {nb} bands, {ns} species")
    (OUT / "counts.tex").write_text(
        "% generated by make_overview.py: do not edit\n"
        + "\n".join(counts) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
