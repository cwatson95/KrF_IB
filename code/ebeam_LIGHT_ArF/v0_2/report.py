"""Energy accounting and the text report (the notes' 'Output' line).

LIGHT_ArF (2026-09-30): the Ar/F2 version of LIGHT/v0_2/report.py. The
ledger and its closure test are unchanged; the species, bands and warnings
are ArF's. LIGHT_ArF_ebeam: under the beam drive the beam makes Ar+, Ar*
and F without a row, so the flows add its split and the Ar* reservoir
counts its excited atoms as formed (pairs = E_dep / W exactly); its prompt
heat is in the reaction heat (ledger 'heat_beam').

Two views of the same ledger, both in percent of the deposited energy:

FLOWS — what passed through each channel ("percent to Ar+, ArF(B), Ar*,
heat, elastic losses and emission losses"). Flows overlap: energy that
reaches Ar* goes on to ArF(B) and then to emission, so they do not sum.

FATE AT t_end — where the energy sits when the run stops: radiated by band,
heat (elastic and reaction, separately), or still stored in a species or
the electron gas. These sum to 100 % up to the printed closure residual,
which is the model's bookkeeping test (events x per-event energies vs the
densities the ODE produced).

Reservoir accounting treats the ArF excimer states (species.MANIFOLD) as
ONE manifold (REBUILD 0.1.24 lesson: state-mixing rows swamp per-state
bookkeeping).
"""

from __future__ import annotations

import numpy as np

from . import network as N
from .model import QE
from .network import BANDS
from .species import (EXCIMER, MANIFOLD, RG, RG_ION, RG_STAR, TRIMER,
                      e_form)

DIMER = "Ar2*"
DIMER_ION = "Ar2+"
RESERVOIRS = {
    RG_STAR: (RG_STAR,),
    "ArF*": MANIFOLD,
    TRIMER: (TRIMER,),
    DIMER: (DIMER,),
}
STORE_GROUPS = (
    ("ions (Ar+, Ar2+, F-)", (RG_ION, DIMER_ION, "F-")),
    (RG_STAR, (RG_STAR,)),
    (DIMER, (DIMER,)),
    ("+".join(MANIFOLD), MANIFOLD),
    (TRIMER, (TRIMER,)),
    ("F atoms (F2 broken)", ("F",)),
)
HOT_K = 1000.0            # redissociation of Ar2+ and Ar2F* not modeled
BX_BAND = N.BX_BAND
VUV_BAND = N.VUV_BAND


def _pct(x, e_dep):
    return 100.0 * x / e_dep if e_dep > 0 else float("nan")


def accounting(res) -> dict:
    """All energies in eV/cm^3 at t_end, plus reservoir branchings."""
    m = res.model
    e_dep = float(res.ledger["E_dep"][-1])
    ev = {rid: float(v[-1]) for rid, v in res.events.items()}
    rows = {c.row.rid: c for c in m.compiled}

    bands = {}
    heat_rx = (float(res.ledger["heat_capture"][-1])
               + float(res.ledger["heat_beam"][-1]))
    for rid, c in rows.items():
        if c.route["rad"]:
            b = c.row.band
            bands[b] = bands.get(b, 0.0) + ev[rid] * c.route["rad"]
        heat_rx += ev[rid] * c.route["heat"]
    heat_el = float(res.ledger["heat_elastic"][-1])

    stores = {}
    for label, members in STORE_GROUPS:
        d = sum((res.n[s][-1] - res.n[s][0]) * e_form(s)
                for s in members if s in res.n)
        if any(s in res.n for s in members):
            stores[label] = float(d)
    stores["electron gas"] = float(res.ue[-1] - res.ue[0])

    fate_total = sum(bands.values()) + heat_el + heat_rx + sum(stores.values())
    closure = (e_dep - fate_total) / e_dep if e_dep > 0 else float("nan")

    # --- flows: the notes' output list -----------------------------------
    exc_rids = [r for r in rows if r.startswith(("G1", "C1"))]
    e_exc = sum(ev[r] * rows[r].route["from_e"] for r in exc_rids)
    e_exc_4s = sum(ev[r] * rows[r].route["from_e"] for r in exc_rids
                   if r.startswith("G1"))
    # every row that frees an electron makes an ion; book the cation's
    # formation energy (Ar+ 15.76 eV, Ar2+ 14.53 eV) per event
    e_ion_all = sum(ev[r] * sum(e_form(s) for s in c.row.products
                                if s.endswith("+"))
                    for r, c in rows.items() if c.row.ne_delta > 0)
    e_ion_e = ev.get("C2", 0.0) * e_form(RG_ION)
    arf = reservoir(res, "ArF*")
    ars = reservoir(res, RG_STAR)
    hv_bx = rows["G4"].route["rad"]
    e_to_e = sum(ev[r] * c.route["to_e"] for r, c in rows.items())
    flows = {
        "to Ar* (e-impact excitation, all levels)": e_exc,
        "   of which 4s direct (G1)": e_exc_4s,
        "   of which 4p/3d/5s/higher -> 4s (C1)": e_exc - e_exc_4s,
        "to Ar+ (every ionization x ion energy)": e_ion_all,
        "   of which e-impact on Ar (C2)": e_ion_e,
        f"to ArF(B) (formed x {hv_bx:.2f} eV quantum)": arf["formed"] * hv_bx,
        "heat (reactions + relaxation)": heat_rx,
        "elastic losses (G6)": heat_el,
        "emission losses (all bands)": sum(bands.values()),
        "returned to electrons (Penning etc.)": e_to_e,
    }
    if m.beam is not None:
        y, pairs = m.beam, e_dep / m.beam.W
        flows.update({
            "beam: to Ar+ (every pair x 15.76 eV)": pairs * e_form(RG_ION),
            "beam: to Ar* (excitations x 11.548 eV)":
                pairs * y.ar_star * e_form(RG_STAR),
            "beam: prompt heat (above 4s, F2 fragments)": pairs * y.heat,
            "beam: to the electron gas (secondaries)": pairs * y.eps_sec,
            "beam: to F atoms (prompt F2 dissociation)":
                pairs * 2.0 * y.f2_diss * e_form("F"),
        })
    return {"e_dep": e_dep, "bands": bands, "heat_elastic": heat_el,
            "heat_reactions": heat_rx, "stores": stores, "closure": closure,
            "flows": flows, "arf": arf, "arstar": ars, "events": ev,
            "hv_bx": hv_bx}


def reservoir(res, name) -> dict:
    """Formation and per-row losses of a reservoir (species set)."""
    members = RESERVOIRS[name]
    formed, losses = 0.0, {}
    for c in res.model.compiled:
        r = c.row
        net = (sum(r.products.count(s) for s in members)
               - sum(r.reactants.count(s) for s in members))
        e = float(res.events[r.rid][-1])
        if net > 0:
            formed += net * e
        elif net < 0:
            losses[r.rid] = -net * e
    y = res.model.beam
    if y is not None:                      # made by the beam, not by a row
        net = sum(y.per_pair().get(s, 0.0) for s in members)
        made = net * float(res.ledger["E_dep"][-1]) / y.W
        if made > 0:
            formed += made
        elif made < 0:
            losses["beam"] = -made
    left = sum(res.n[s][-1] - res.n[s][0] for s in members if s in res.n)
    return {"formed": formed, "losses": losses, "left": float(left),
            "balance": formed - sum(losses.values()) - left}


def f2_burnout_time(res, frac=0.1):
    """First time [s] F2 falls below ``frac`` of its initial density."""
    f2 = res.n["F2"]
    below = np.nonzero(f2 < frac * f2[0])[0]
    return float(res.t[below[0]]) if below.size else None


def warnings(res) -> list:
    out = []
    from . import xsec
    if res.te.max() > xsec.TE_MAX:
        out.append(f"Te reached {res.te.max():.1f} eV, above the rate-table "
                   f"edge {xsec.TE_MAX:.0f} eV: rows on Maxwellian tables were "
                   "held at the edge, two-temperature tail rows follow T_tail "
                   "= r Te up to it, and power-law rows follow their laws")
    tmax = float(res.tg.max())
    for rid, t_valid, text in N.TG_VALIDITY:
        rows = [c for c in res.model.compiled if c.row.rid == rid]
        if rows and tmax > t_valid:
            out.append(f"Tg reached {tmax:.0f} K. {text}")
    if tmax > HOT_K:
        out.append(f"Tg above {HOT_K:.0f} K: Ar2+ (bound by 1.23 eV) and "
                   "Ar2F* (0.65 eV) redissociate on these time scales, "
                   "which the model does not include, and the ion "
                   "conversion and ion-ion laws are extrapolated beyond "
                   "their data")
    tb = f2_burnout_time(res)
    if tb is not None:
        out.append(f"F2 burn-out: below 10 % of the fill at t = "
                   f"{tb * 1e9:.0f} ns. F + F + Ar re-forms F2 on a "
                   "millisecond scale, far too slowly to refill it")
    rg_loss = 1.0 - res.n[RG][-1] / res.n[RG][0]
    if rg_loss > 0.05:
        out.append(f"ground-state Ar depleted by {100 * rg_loss:.1f} %")
    return out


def _fmt_pct(x, e_dep):
    return f"{_pct(x, e_dep):7.2f} %"


def summary(res) -> str:
    p, a = res.p, accounting(res)
    e = a["e_dep"]
    j = e * QE
    lines = []
    y = res.model.beam
    drive = (f"P = {p.power_W_cm3 / 1e3:.0f} kW/cm^3 into the electron gas"
             if y is None else
             f"e-beam P = {p.power_W_cm3 / 1e3:.0f} kW/cm^3 (W {y.W:.2f} eV "
             f"per pair, Ar* {y.ar_star:.3f} per pair, secondaries "
             f"{y.eps_sec:.2f} eV, {y.route}/{y.sharing})")
    groups = ",".join(p.groups) if p.groups else "none"
    lines.append("LIGHT_ArF_ebeam: 0-D Ar/F2 ArF model (v0_2)")
    lines.append(f"  fill {p.pressure_atm:.3g} atm, F2 {100 * p.x_f2:.3g} % "
                 f"in Ar, Tg0 {p.tg0_K:.0f} K;  {drive} for "
                 f"{p.pulse_s * 1e6:.3g} us;  n_e seed {p.ne_cm3:.2g} cm^-3 "
                 f"({res.model.ne_mode});  t_end {p.t_end_s * 1e6:.3g} us")
    lines.append(f"  network {p.network};  tau_B {p.tau_sp_s * 1e9:.2f} ns;"
                 f"  EEDF {res.model.eedf};  optional groups: {groups}")
    read = {b for c in res.model.compiled for b in N.bands_of(c.row)}
    if read:
        set_ = [f"{k} {v:.3g}" for k, v in p.rates]
        if p.trimer_tg_pow is not None:
            set_.append(f"trimer_tg_pow {p.trimer_tg_pow:.3g}")
        off = [k for k in sorted(read) if BANDS[k].default == 0.0
               and k not in res.model.rates]
        lines.append("  bands: " + (", ".join(set_) if set_ else
                                    "every band at its default")
                     + (";  no measurement, 0: " + ", ".join(off)
                        if off else "") + "  (--bands runs each end)")
    lines.append(f"  deposited energy E_dep = {j:.4g} J/cm^3")
    lines.append("")
    lines.append("FLOWS (the notes' outputs; percent of E_dep, they overlap)")
    for k, v in a["flows"].items():
        lines.append(f"  {k:<44s}{_fmt_pct(v, e)}")
    lines.append("")
    lines.append("FATE AT t_end (sums to 100 %)")
    for b, v in a["bands"].items():
        lines.append(f"  emission  {b:<34s}{_fmt_pct(v, e)}")
    lines.append(f"  heat      {'elastic collisions':<34s}"
                 f"{_fmt_pct(a['heat_elastic'], e)}")
    lines.append(f"  heat      {'reactions + relaxation':<34s}"
                 f"{_fmt_pct(a['heat_reactions'], e)}")
    for k, v in a["stores"].items():
        lines.append(f"  stored    {k:<34s}{_fmt_pct(v, e)}")
    lines.append(f"  closure residual (bookkeeping test)   "
                 f"{a['closure']:.1e}")
    lines.append("")
    for name, title in ((RG_STAR, "Ar* (4s) fate"),
                        ("ArF*", "ArF* fate")):
        rv = reservoir(res, name)
        f = rv["formed"]
        if f <= 0:
            continue
        parts = ", ".join(f"{rid} {100 * v / f:.1f} %" for rid, v in
                          sorted(rv["losses"].items(), key=lambda kv: -kv[1]))
        lines.append(f"{title}: formed {f:.3g} cm^-3 -> {parts}, "
                     f"left {100 * rv['left'] / f:.1f} %")
    t = res.t
    on = t <= p.pulse_s
    lines.append("")
    lines.append(f"state: Te pulse-mean {np.mean(res.te[on]):.2f} eV "
                 f"(max {res.te.max():.2f});  Tg end {res.tg[-1]:.0f} K;  "
                 f"p end {res.pressure_atm[-1]:.3g} atm;  F2 left "
                 f"{100 * res.n['F2'][-1] / res.n['F2'][0]:.2g} %")
    lines.append(f"       n_e end {res.ne[-1]:.3g} cm^-3 "
                 f"(peak {res.ne.max():.3g})"
                 f";  Ar ground left "
                 f"{100 * res.n[RG][-1] / res.n[RG][0]:.3g} %")
    lines.append(f"       peaks: Ar* {res.n[RG_STAR].max():.3g}, ArF(B) "
                 f"{res.n[EXCIMER].max():.3g}, Ar2F* "
                 f"{res.n[TRIMER].max():.3g} cm^-3;  RHS evals {res.nfev}")
    w = warnings(res)
    if w:
        lines.append("")
        lines.append("WARNINGS")
        lines.extend("  - " + s for s in w)
    return "\n".join(lines)


def row_dict(res) -> dict:
    """Flat numbers for scans / CSV (percent of E_dep unless stated)."""
    p, a = res.p, accounting(res)
    e = a["e_dep"]
    arf = a["arf"]
    ar = a["arstar"]
    band = a["bands"]
    flows = a["flows"]
    d = {
        "pressure_atm": p.pressure_atm, "x_f2": p.x_f2,
        "power_kW_cm3": p.power_W_cm3 / 1e3, "drive": p.drive,
        "ne_cm3": p.ne_cm3,
        "ne_mode": res.model.ne_mode,
        "groups": "+".join(p.groups) or p.network,
        "E_dep_J_cm3": e * QE,
        "to_Arstar": _pct(flows["to Ar* (e-impact excitation, all "
                                "levels)"], e),
        "to_Arplus": _pct(flows["to Ar+ (every ionization x ion energy)"],
                          e),
        "to_ArFB": _pct(arf["formed"] * a["hv_bx"], e),
        "heat_reactions": _pct(a["heat_reactions"], e),
        "elastic": _pct(a["heat_elastic"], e),
        "emission_193": _pct(band.get(BX_BAND, 0.0), e),
        "emission_vuv": _pct(band.get(VUV_BAND, 0.0), e),
        "emission_all": _pct(sum(band.values()), e),
        "stored_ions": _pct(a["stores"].get("ions (Ar+, Ar2+, F-)", 0.0), e),
        "stored_F": _pct(a["stores"].get("F atoms (F2 broken)", 0.0), e),
        "arf_to_193_pct": 100 * arf["losses"].get("G4", 0.0) / arf["formed"]
        if arf["formed"] else float("nan"),
        "arstar_to_arf_pct": 100 * sum(
            v for rid, v in ar["losses"].items() if rid in N.HARPOON_RIDS)
        / ar["formed"] if ar["formed"] else float("nan"),
        "te_mean_eV": float(np.mean(res.te[res.t <= p.pulse_s])),
        "tg_end_K": float(res.tg[-1]),
        "f2_left_pct": 100 * res.n["F2"][-1] / res.n["F2"][0],
        "closure": a["closure"],
    }
    return d
