"""Write a kit (make_kit.build_kit) as COMSOL files.

Conventions of every file:
  * SI with COMSOL's unit brackets; energies per event in eV are numbers
    times e_const*1[V];
  * the closure function kin_<gas>(lnTe, lnne, i): lnTe = ln(Te/1 eV),
    lnne = ln(n_e/1 cm^-3), i the quantity's index (an exact integer, so the
    linear interpolation in i returns that quantity alone); linear in the
    other two, constant beyond the grid;
  * every electron energy term is in the variable Qe (the Plasma Module's
    General Power Deposition), so every Electron Impact Reaction carries an
    energy loss of 0: the ledger is LIGHT's own (energy_routing), with
    nothing counted twice;
  * rate constants for the Plasma Module's reactions are the variables
    kf_<reaction>, molar (k N_A^(order-1)), as its typed rate constants are.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
KIT = HERE.parent
UNITS = {1: "1[1/s]", 2: "1e-6[m^3/s]", 3: "1e-12[m^6/s]"}


def fmt(x: float) -> str:
    return f"{x:.10g}"


def nvar(sp) -> str:
    return "n_" + sp["tag"]


# ---------------------------------------------------------------------------
def write_grid(path: Path, x, y, data, z=None, header=()):
    """COMSOL Grid format: '% Grid', one line per argument, '% Data', then
    one row per (y, z) with y varying fastest."""
    with open(path, "w") as f:
        for h in header:
            f.write(f"% {h}\n")
        f.write("% Grid\n")
        f.write(" ".join(fmt(v) for v in x) + "\n")
        f.write(" ".join(fmt(v) for v in y) + "\n")
        if z is not None:
            f.write(" ".join(fmt(v) for v in z) + "\n")
        f.write("% Data\n")
        for row in data:
            f.write(" ".join(fmt(v) for v in row) + "\n")


def kinetics_rows(kit):
    """The 3-D grid's data rows: for each index, for each n_e slice, the
    values along ln Te."""
    rows = []
    for q in kit["grid"]["quantities"]:
        for s in q["values"]:
            rows.append(s)
    return rows


# ---------------------------------------------------------------------------
def species_by_light(kit):
    return {s["light"]: s for s in kit["species"]}


def term(kit, side):
    sp = species_by_light(kit)
    return [s for s in side if s != "e"], [sp[s] for s in side if s != "e"]


def rate_expr(kit, rec, electron: bool) -> str:
    gas = kit["gas"]
    if electron:
        return f"exp(kin_{gas}(lnTe,lnne,{rec['i_k']}))*1e-6[m^3/s]"
    law = rec["law"]
    unit = UNITS[rec["order"]]
    if law["type"] == "power":
        k = fmt(law["k300"])
        if law["n"] == 0.0:
            return f"{k}*{unit}"
        return f"{k}*(Tg_/300[K])^({fmt(law['n'])})*{unit}"
    if law["type"] == "balance":
        e = f"{fmt(law['a'])}*exp(-{fmt(law['theta'])}[K]*(1/Tg_-1/300[K]))"
        if law["n"]:
            e += f"*(Tg_/300[K])^({fmt(law['n'])})"
        return f"{e}*{unit}"
    if law["type"] == "ionion":
        a3 = f"{fmt(law['k3'])}*Nam*(Tg_/300[K])^({fmt(law['n3'])})"
        al = f"{fmt(law['kL'])}/Nam"
        e = f"({a3})*({al})/(({a3})+({al}))*Sscr*{fmt(law['f'])}"
        if law["n"]:
            e += f"*(Tg_/300[K])^({fmt(law['n'])})"
        return f"{e}*{unit}"
    if law["type"] == "tgtab":
        return f"exp(heavy_{gas}(log(Tg_/1[K]),{law['index']}))*{unit}"
    raise ValueError(law)


def flux_expr(kit, rec, electron: bool) -> str:
    sp = species_by_light(kit)
    parts = [f"k_{rec['name']}"]
    if electron:
        parts.append("ne_")
    for s in rec["reactants"]:
        if s == "e":
            continue
        parts.append(nvar(sp[s]))
    return "*".join(parts)


def variables(kit, route: str) -> list:
    """[(name, expression, description)] for route 'plasma' (the Plasma
    Module model) or 'ode' (the Domain-ODE model)."""
    gas = kit["gas"]
    sp = kit["species"]
    V = []

    def v(name, expr, desc):
        V.append((name, expr, desc))
    cations = [s for s in sp if s["charge"] > 0]
    anions = [s for s in sp if s["charge"] < 0]
    if route == "plasma":
        v("Te_eV", "plas.Te/1[V]", "electron temperature [eV], from the Plasma interface")
        v("ne_", "plas.ne", "electron density")
        for s in sp:
            v(nvar(s), f"plas.n_w{s['tag']}", f"number density of {s['light']}")
    else:
        for s in sp:
            v(nvar(s), f"n{s['tag']}", f"number density of {s['light']} (Domain ODE)")
        cat = "+".join(nvar(s) for s in cations)
        an = "-".join([""] + [nvar(s) for s in anions]) if anions else ""
        v("ne_", f"max({cat}{an},1[1/m^3])", "electron density: cations - anions (quasineutral, LIGHT's)")
        v("Te_eV", "max(Ue/(1.5*ne_*e_const*1[V]),1e-3)", "electron temperature [eV]: Ue = (3/2) n_e Te")
    v("Tg_", "Tg", "gas temperature (Domain ODE)")
    v("lnTe", "log(max(Te_eV,0.01))", "ln(Te / 1 eV): the closure's first argument")
    v("lnne", "log(max(ne_*1e-6[m^3],1))", "ln(n_e / 1 cm^-3): the closure's second argument")
    v("N_", "+".join(nvar(s) for s in sp), "heavy-particle density")
    if any(r["law"]["type"] == "ionion" for r in kit["heavy_reactions"]):
        an = "+".join(nvar(s) for s in anions)
        v("Nam", "N_/2.6867811e25[1/m^3]", "density in amagat")
        v("n_scr", f"max(2*(ne_+{an}),1[1/m^3])", "screening density (LIGHT's _screening)")
        v("lam_D", "sqrt(epsilon0_const*k_B_const*Tg_/(n_scr*e_const^2))", "Debye length")
        v("r_cap", "2*e_const^2/(3*4*pi*epsilon0_const*k_B_const*Tg_)", "capture radius")
        v("Sscr", "1/(1+r_cap/lam_D)", "Debye screening of ion-ion recombination")
    for rec in kit["electron_reactions"]:
        v(f"k_{rec['name']}", rate_expr(kit, rec, True),
          f"{' + '.join(rec['reactants'])} -> {' + '.join(rec['products'])}"
          + (f" ({len(rec['rows'])} rows lumped)" if len(rec['rows']) > 1 else ""))
    for rec in kit["heavy_reactions"]:
        v(f"k_{rec['name']}", rate_expr(kit, rec, False),
          f"{' + '.join(rec['reactants'])} -> {' + '.join(rec['products'])}"
          + (f" ({len(rec['rows'])} rows lumped)" if len(rec['rows']) > 1 else ""))
    for rec in kit["electron_reactions"]:
        v(f"R_{rec['name']}", flux_expr(kit, rec, True), "event rate")
    for rec in kit["heavy_reactions"]:
        v(f"R_{rec['name']}", flux_expr(kit, rec, False), "event rate")
    if route == "plasma":
        for rec in kit["electron_reactions"]:
            v(f"kf_{rec['name']}", f"k_{rec['name']}*N_A_const",
              "molar rate constant for the Electron Impact Reaction")
        for rec in kit["heavy_reactions"]:
            o = rec["order"]
            f = "" if o == 1 else ("*N_A_const" if o == 2 else "*N_A_const^2")
            v(f"kf_{rec['name']}", f"k_{rec['name']}{f}", "molar rate constant for the Reaction")
    # species sources
    names = {s["light"]: s for s in sp}
    for s in sp:
        terms = []
        for rec, el in [(r, True) for r in kit["electron_reactions"]] + \
                       [(r, False) for r in kit["heavy_reactions"]]:
            nu = rec["products"].count(s["light"]) - rec["reactants"].count(s["light"])
            if nu:
                terms.append(("+" if nu > 0 else "-")
                             + (f"{abs(nu)}*" if abs(nu) != 1 else "") + f"R_{rec['name']}")
        expr = "".join(terms).lstrip("+") or "0[1/(m^3*s)]"
        v(f"dn_{s['tag']}", expr, f"net production of {s['light']}")
    # energies
    gi = kit["grid"]
    v("pel_", f"kin_{gas}(lnTe,lnne,{gi['i_pel']})*1e-6[m^3/s]",
      "Kr elastic power per electron per Kr atom [eV m^3/s as V*m^3/s numbers]")
    kr = nvar(names["Kr"])
    v("Pel", f"e_const*1[V]*pel_*ne_*{kr}", "electron power lost in elastic collisions with Kr")
    de, heat, cap = [], [], []
    for rec in kit["electron_reactions"]:
        n = rec["name"]
        if "i_de" in rec:
            de.append(f"kin_{gas}(lnTe,lnne,{rec['i_de']})*R_{n}")
        elif rec["de"]:
            de.append(f"({fmt(rec['de'])})*R_{n}")
        if "i_heat" in rec:
            heat.append(f"kin_{gas}(lnTe,lnne,{rec['i_heat']})*R_{n}")
        elif rec["heat"]:
            heat.append(f"({fmt(rec['heat'])})*R_{n}")
        if "i_cap" in rec:
            c = f"{rec['capture']}*" if rec["capture"] != 1 else ""
            cap.append(f"{c}kin_{gas}(lnTe,lnne,{rec['i_cap']})*R_{n}")
    for rec in kit["heavy_reactions"]:
        n = rec["name"]
        if rec["de"]:
            de.append(f"({fmt(rec['de'])})*R_{n}")
        if rec["heat"]:
            heat.append(f"({fmt(rec['heat'])})*R_{n}")
    join = lambda xs: "+".join(xs) if xs else "0"
    v("Qcap", f"e_const*1[V]*({join(cap)})" if cap else "0[W/m^3]",
      "energy the captured electrons take from the electron gas to the gas")
    v("Qe_rx", f"e_const*1[V]*({join(de)})", "electron energy gained (+) or paid (-) in reactions")
    v("Qe", "Qe_rx-Qcap-Pel",
      "electron power from everything but the wave (the General Power Deposition)")
    v("Qg", f"e_const*1[V]*({join(heat)})+Qcap+Pel", "power into the gas (heat)")
    v("cvn", "k_B_const*(" + "+".join(f"{fmt(s['cv'])}*{nvar(s)}" for s in sp) + ")",
      "constant-volume heat capacity of the gas")
    v("dTg", "(Qg-k_B_const*Tg_*(" + "+".join(f"{fmt(s['cv'])}*dn_{s['tag']}" for s in sp)
      + "))/cvn", "gas temperature rate (LIGHT's: created particles take cv k Tg)")
    bands = {}
    for rec in kit["heavy_reactions"]:
        for b, e in rec["rad"].items():
            bands.setdefault(b, []).append(f"({fmt(e)})*R_{rec['name']}")
    for b, terms in bands.items():
        v(band_var(b), f"e_const*1[V]*({join(terms)})", f"emitted power, {b}")
    v("nu_e", f"exp(kin_{gas}(lnTe,lnne,{gi['i_nu']}))*1[1/s]",
      "effective collision frequency for the wave (the kinetic absorption)")
    if route == "plasma":
        v("muN_e", "e_const*N0/(me_const*nu_e)",
          "reduced electron mobility for the Plasma Model: nu_m = e/(m mu) = nu_e")
        v("Qw", "emw.Qrh", "absorbed wave power")
    else:
        v("sig_r", "ne_*e_const^2*nu_e/(me_const*(nu_e^2+omega0^2))", "Re(sigma), Drude")
        v("sig_i", "-ne_*e_const^2*omega0/(me_const*(nu_e^2+omega0^2))", "Im(sigma), Drude")
        v("epsr_D", "1+sig_i/(omega0*epsilon0_const)-j*sig_r/(omega0*epsilon0_const)",
          "relative permittivity of the plasma (exp(+j w t))")
        v("Qw", "0.5*sig_r*emw.normE^2", "absorbed wave power")
    return V


def band_var(b: str) -> str:
    return "P" + b.split()[0].replace(".", "p").replace("-", "_")


# ---------------------------------------------------------------------------
def params(kit) -> list:
    gas = kit["gas"]
    ini = kit["initial"]
    sp = kit["species"]
    P = []

    def p(name, expr, desc):
        P.append((name, expr, desc))
    p("L_gas", "10[cm]", "gas depth along the wave")
    p("h_strip", "0.1[mm]", "height of the quasi-1-D strip (PMC/insulating sides)")
    p("L_vac", "0.5[mm]", "vacuum in front of and behind the gas")
    p("f0", "1.3[THz]", "wave frequency (1.3 THz: absorbed over ~10 cm in the afterglow)")
    p("omega0", "2*pi*f0", "angular frequency")
    p("lambda0", "c_const/f0", "vacuum wavelength")
    p("hmax", "lambda0/6", "largest element along the wave")
    p("I0_MW", "3", "incident intensity [MW/cm^2]; 0 for the off twin")
    p("I0", "I0_MW*1e6[W/cm^2]", "incident intensity")
    p("E0", "sqrt(2*Z0_const*I0)", "peak amplitude of the incident field")
    p("t_end", "1[us]", "afterglow simulated (from the end of the pulse)")
    p("N0", f"{fmt(kit['N0_cm3'] * 1e6)}[1/m^3]", "heavy-particle density (constant volume)")
    p("Tg0", f"{fmt(ini['tg_K'])}[K]", "gas temperature at the end of the pulse")
    p("ne0", f"{fmt(ini['ne_cm3'] * 1e6)}[1/m^3]", "electron density at the end of the pulse")
    p("Te0", f"{fmt(ini['te_eV'])}[V]", "electron temperature at the end of the pulse")
    p("ebar0", "1.5*Te0", "mean electron energy at the end of the pulse")
    p("Ue0", f"1.5*ne0*Te0*e_const", "electron energy density at the end of the pulse")
    tot = sum(ini["n_cm3"][s["light"]] * s["mass_kg_mol"] for s in sp)
    for s in sp:
        n = ini["n_cm3"][s["light"]]
        p(f"n0_{s['tag']}", f"{fmt(n * 1e6)}[1/m^3]", f"{s['light']} at the end of the pulse")
    for s in sp:
        n = ini["n_cm3"][s["light"]]
        p(f"w0_{s['tag']}", fmt(n * s["mass_kg_mol"] / tot), f"mass fraction of {s['light']}")
    pp = kit["pump_phase"]
    p("E_given", f"{fmt(pp['E_given_J_cm3'])}[J/cm^3]*L_gas", "pump energy (300 kW/cm3 x 1 us) per area")
    for b, e in pp["bands_J_cm3"].items():
        p(f"Epump_{band_var(b)[1:]}", f"{fmt(e)}[J/cm^3]*L_gas", f"{b} emitted during the pulse, per area")
    return P


def write_table(path: Path, rows, header=()):
    """name expression description, tab-separated (COMSOL Load from File)."""
    with open(path, "w") as f:
        for name, expr, desc in rows:
            f.write(f"{name}\t{expr}\t{desc}\n")


def chemistry(kit) -> str:
    gas = kit["gas"]
    sp = kit["species"]
    names = {s["light"]: s["name"] for s in sp}
    L = [f"% {gas}: species and reactions of {kit['model']} for COMSOL's",
         "% Plasma Chemistry Import (6.4, Global) or the Plasma Chemistry add-in.",
         "% Rate constants are the variables kf_<reaction> (molar, units carried);",
         "% load the variables file first. Energy losses are 0: every electron",
         "% energy term is in the General Power Deposition Qe.",
         "SPECIES", "USERDEFINED"]
    for s in sp:
        sig, eps = s["lj"]
        L.append(f"{s['tag']} ! {s['mass_kg_mol']:.6f} {sig} {eps} 0")
    L += ["END", "ELECTRONIMPACTREACTIONS", "RATECONSTANT"]
    skipped = []
    for rec in kit["electron_reactions"]:
        r, p = rec["reactants"], rec["products"]
        if sorted(r) == sorted(p):
            skipped.append(rec["name"])
            continue
        dn = p.count("e") - r.count("e")
        kind = "Ionization" if dn > 0 else "Attachment" if dn < 0 else "Excitation"
        lhs = "+".join(["e"] + [names[s] for s in r if s != "e"])
        rhs = "+".join(["e"] * p.count("e") + [names[s] for s in p if s != "e"])
        if kind == "Attachment":
            L.append(f"{lhs}=>{rhs} ! Attachment kf_{rec['name']}")
        else:
            L.append(f"{lhs}=>{rhs} ! {kind} 0 kf_{rec['name']}")
    L += ["END", "HEAVYSPECIESREACTIONS", "RATECONSTANT"]
    for rec in kit["heavy_reactions"]:
        r, p = rec["reactants"], rec["products"]
        lhs = "+".join(names[s] for s in r)
        rhs = "+".join(["e"] * p.count("e") + [names[s] for s in p if s != "e"])
        L.append(f"{lhs}=>{rhs} ! kf_{rec['name']}")
    L += ["END"]
    if skipped:
        L.insert(5, "% Left out (no species change; their energy is in Qe): "
                 + ", ".join(skipped))
    return "\n".join(L) + "\n"


def write_all(kit, ref=None):
    gas = kit["gas"]
    out = KIT / gas
    out.mkdir(exist_ok=True)
    g = kit["grid"]
    nq = len(g["quantities"])
    head = [f"{gas}: the kinetic electron closure of {kit['model']}, frozen in tables",
            "arguments: lnTe = ln(Te/eV), lnne = ln(ne/cm^-3), index (1.." + str(nq) + ")"]
    head += [f"index {i + 1}: {q['key']} - {q['desc']}" for i, q in enumerate(g["quantities"])]
    write_grid(out / f"kinetics_{gas}.txt", g["lnte"], g["lnne"], kinetics_rows(kit),
               z=list(range(1, nq + 1)), header=head)
    if kit["heavy_grid"]["values"]:
        hv = kit["heavy_grid"]
        idx = [r["law"]["index"] for r in kit["heavy_reactions"] if r["law"]["type"] == "tgtab"]
        write_grid(out / f"heavy_{gas}.txt", hv["lntg"], list(range(1, len(hv["values"]) + 1)),
                   hv["values"], header=[f"{gas}: ln k of the heavy reactions whose rate "
                                         "depends on Tg other than as a power law",
                                         "arguments: ln(Tg/K), index; values ln k in LIGHT's "
                                         "cm units (1/s, cm3/s, cm6/s by order)"])
    write_table(out / f"params_{gas}.txt", params(kit))
    write_table(out / f"variables_{gas}_plasma.txt", variables(kit, "plasma"))
    write_table(out / f"variables_{gas}_ode.txt", variables(kit, "ode"))
    (out / f"chemistry_{gas}.txt").write_text(chemistry(kit))
    (out / f"kit_{gas}.json").write_text(json.dumps(kit, indent=1))
    return out
