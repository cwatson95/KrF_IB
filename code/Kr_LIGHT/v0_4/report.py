"""Energy accounting and the text report (the notes' 'Output' line), for a
fill of pure krypton (LIGHT_Kr v0_4).

Two views of the same ledger, both in percent of the deposited energy:

FLOWS — what passed through each channel (the notes' "percent to Kr+,
Kr*, heat, elastic losses and emission losses", with Kr2*, the emitter of a
krypton fill, in KrF(B)'s place). Flows overlap: energy that reaches Kr*
goes on to Kr2* and then to emission, so they do not sum.

FATE AT t_end — where the energy sits when the run stops: radiated by band,
heat (elastic and reaction, separately), or still stored in a species or
the electron gas. These sum to 100 % up to the printed closure residual,
which is the model's bookkeeping test (events x per-event energies vs the
densities the ODE produced).

LIGHT v0_5's report, with its fluorine lines (248 nm, KrF*, Kr2F*, F
atoms, F2 burn-out, the trimer's temperature law) and its 248-nm optics
taken out. v0_2 sums the level rows by their groups (G1_*, C1a_* to
C1c_*, HM_*; network.group_rows) where v0_1 read one row per group.
"""

from __future__ import annotations

import numpy as np

from .model import QE
from .network import BANDS, KR_5S, OFF_BY_DEFAULT, bands_read
from .species import e_form

# v0_4: Kr* is the 5s manifold, the four levels' sum (transfers between
# them move nothing in or out), and Kr2* both of the excimer's states
RESERVOIRS = {
    "Kr*": KR_5S,
    "Kr**": ("Kr**",),
    "Kr2*": ("Kr2*", "Kr2*(0u+)"),
}
STORE_GROUPS = (
    ("ions (Kr+, Kr2+)", ("Kr+", "Kr2+")),
    ("Kr*", KR_5S),
    ("Kr** (5p)", ("Kr**",)),
    ("Kr2*", ("Kr2*", "Kr2*(0u+)")),
)
# Kr* + 2Kr's Tg law, Wieme & Lenaerts: measured over 200-400 K (their
# title and abstract; v0_1 to v0_3 said 300-451 K, the range of the
# power-law conversion, not the measurement's, 2026-10-15)
KR2_TG_VALID_K = 400.0
HOT_K = 1000.0            # above it Kr2+ redissociates (K15r, v0_4)
VUV_BAND = "147 nm (Kr2* VUV)"
NIR_BAND = "5p-5s (near IR)"


def _pct(x, e_dep):
    return 100.0 * x / e_dep if e_dep > 0 else float("nan")


def _kr2_label(hv):
    return f"to Kr2* (formed x {hv:.2f} eV quantum)"


def accounting(res) -> dict:
    """All energies in eV/cm^3 at t_end, plus reservoir branchings. e_dep
    is everything put in: the given power and the EM wave's absorbed
    energy (E_em)."""
    m = res.model
    e_given = float(res.ledger["E_dep"][-1])
    e_em = float(res.ledger["E_em"][-1])
    e_dep = e_given + e_em
    ev = {rid: float(v[-1]) for rid, v in res.events.items()}
    rows = {c.row.rid: c for c in m.compiled}

    bands = {}
    # the ambipolar field's work on the ions (E_ion_heat) heats the gas and
    # belongs to no row
    heat_rx = (float(res.ledger["heat_capture"][-1])
               + float(res.ledger["E_ion_heat"][-1]))
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
    # the walls take the lost charges' formation energy and the
    # electrons' 2 Te (E_wall). The thermal energy the lost particles
    # carry out (E_wall_th) is part of the heat already booked, so it is
    # not a second fate line (it would be counted twice)
    walls = float(res.ledger["E_wall"][-1])

    fate_total = (sum(bands.values()) + heat_el + heat_rx
                  + sum(stores.values()) + walls)
    closure = (e_dep - fate_total) / e_dep if e_dep > 0 else float("nan")

    # --- flows: the notes' output list -----------------------------------
    exc_rids = [r for r in rows if r.startswith(("G1_", "C1"))]
    e_exc = sum(ev[r] * rows[r].route["from_e"] for r in exc_rids)
    e_exc_5s = sum(ev[r] * rows[r].route["from_e"] for r in exc_rids
                   if r.startswith("G1_"))
    # every row that frees an electron makes an ion; book the cation's
    # formation energy (Kr+ 13.9996 eV, Kr2+ 12.8496 eV) per event
    e_ion_all = sum(ev[r] * sum(e_form(s) for s in c.row.products
                                if s.endswith("+"))
                    for r, c in rows.items() if c.row.ne_delta > 0)
    e_ion_e = ev.get("C2", 0.0) * e_form("Kr+")
    # v0_2: the Hornbeck-Molnar rows' ions (Kr2+), at the band hm_kr
    e_ion_hm = sum(ev[r] * e_form("Kr2+") for r in rows
                   if r.startswith("HM_"))
    kr2 = reservoir(res, "Kr2*")
    krs = reservoir(res, "Kr*")
    hv147 = rows["KR2_RAD"].route["rad"]
    e_to_e = sum(ev[r] * c.route["to_e"] for r, c in rows.items())
    flows = {
        "to Kr* (e-impact excitation, all levels)": e_exc,
        "   of which 5s direct (G1)": e_exc_5s,
        "   of which 5p/4d/higher -> Kr** (C1)": e_exc - e_exc_5s,
        "to Kr+ (every ionization x ion energy)": e_ion_all,
        "   of which e-impact on Kr (C2)": e_ion_e,
        "   of which Hornbeck-Molnar (HM)": e_ion_hm,
        _kr2_label(hv147): kr2["formed"] * hv147,
        "heat (reactions + relaxation)": heat_rx,
        "elastic losses (G6)": heat_el,
        "emission losses (all bands)": sum(bands.values()),
        "returned to electrons (Penning etc.)": e_to_e,
    }
    return {"e_dep": e_dep, "e_given": e_given, "e_em": e_em,
            "bands": bands, "heat_elastic": heat_el,
            "heat_reactions": heat_rx, "stores": stores, "walls": walls,
            "closure": closure, "flows": flows, "kr2": kr2, "krstar": krs,
            "hv147": hv147, "events": ev}


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
    left = sum(res.n[s][-1] - res.n[s][0] for s in members if s in res.n)
    return {"formed": formed, "losses": losses, "left": float(left),
            "balance": formed - sum(losses.values()) - left}


def warnings(res) -> list:
    out = []
    from . import xsec
    if res.te.max() > xsec.TE_MAX:
        out.append(f"Te reached {res.te.max():.1f} eV, above the rate-table "
                   f"edge {xsec.TE_MAX:.0f} eV: rates were held at the edge")
    ak8 = [c for c in res.model.compiled if c.row.rid.startswith("AK8_")]
    if any(c.rate > 0 for c in ak8) and res.tg.max() > KR2_TG_VALID_K:
        out.append(f"Tg reached {res.tg.max():.0f} K. Kr* + 2Kr's "
                   "(Tg/300)^-0.64 is measured over 200-400 K (Wieme & "
                   "Lenaerts 1981) and extrapolated above it")
    if res.tg.max() > HOT_K:
        out.append(f"Tg above {HOT_K:.0f} K: Kr2+ (bound by 1.15 eV) "
                   "redissociates there (K15r, by detailed balance), and "
                   "the ion conversion law is extrapolated beyond its data")
    kr_loss = 1.0 - res.n["Kr"][-1] / res.n["Kr"][0]
    if kr_loss > 0.05:
        out.append(f"ground-state Kr depleted by {100 * kr_loss:.1f} %")
    return out


def _fmt_pct(x, e_dep):
    return f"{_pct(x, e_dep):7.2f} %"


def electrons_line(res) -> str:
    """How the run's electron rates were made. For the kinetic EEDF: how
    many solutions, the pulse's E/N range, the field power over the part
    of the drive a field carries (all of the given power, none of a beam's,
    whose electrons the kinetic equation injects) in the pulse's second
    half, and there the largest gap between the kinetic energy balance and
    LIGHT's (the field power less that part, plus LIGHT's n_e d<eps>/dt,
    which the steady-shape kinetic balance leaves out, over the drive's
    electron power: 0 when the two agree)."""
    m = res.model
    if m.eedf == "two_t":                 # v0_4: which tail table
        return (f"electrons: two_t EEDF ({m.p.tail} tail table), captured "
                f"energy {m.capture}")
    if m.eedf != "kinetic":
        return f"electrons: {m.eedf} EEDF, captured energy {m.capture}"
    log = [e for e in getattr(res, "kinetic", []) if not e.get("skipped")]
    on = [e for e in log if e["on"]]
    stepped = sum(1 for e in log if not e["on"] and e.get("memory"))
    txt = (f"electrons: kinetic EEDF (two-term Boltzmann with e-e, "
           f"boltzmann.py), {len(log)} solutions, {len(on)} in the pulse")
    if stepped:
        txt += f" ({stepped} after it time-stepped)"
    if on:
        en = [e["en_td"] for e in on if e["en_td"] > 0.0]
        late = [e for e in on if e["t"] >= 0.5 * res.p.pulse_s]
        if en:
            txt += f";  pulse E/N {min(en):.3g}-{max(en):.3g} Td"
        # the ratio exists where a field carries the drive (the given
        # power), not where the drive injects its electrons (a beam)
        r = [e["power_ratio"] for e in late
             if e["en_td"] > 0.0 and np.isfinite(e["power_ratio"])]
        b = [abs(e["balance"]) for e in late if np.isfinite(e["balance"])]
        if r:
            txt += (f";  field/input power {np.mean(r):.3f} "
                    f"(min {min(r):.3f}, max {max(r):.3f})")
        if b:
            txt += f"{', ' if r else ';  '}energy balance gap {max(b):.1e}"
    return txt


def bands_line(res) -> str:
    """The bands the run's rows read: the values set, and those with no
    measurement left at their default of 0. (LIGHT v0_5's line listed the
    'param' bands alone and so left out the 0-default stepwise
    ionizations, kstep_s and kstep_p, which its te_fit rows read.)"""
    p, m = res.p, res.model
    read = bands_read([c.row for c in m.compiled])
    set_ = [f"{k} {v:.3g}" for k, v in p.rates]
    zero = [k for k in sorted(read) if BANDS[k].default == 0.0
            and k not in m.rates]
    # v0_4: the bands a measurement bounds, at 0 by choice, apart
    off = [k for k in zero if k not in OFF_BY_DEFAULT]
    chosen = [k for k in zero if k in OFF_BY_DEFAULT]
    return ("  bands: " + (", ".join(set_) if set_ else
                           "every band at its default")
            + (";  no measurement, 0: " + ", ".join(off) if off else "")
            + (";  off by default: " + ", ".join(chosen) if chosen else "")
            + "  (--bands runs each end)")


def summary(res) -> str:
    p, a = res.p, accounting(res)
    m = res.model
    e = a["e_dep"]
    j = e * QE
    lines = []
    drive = f"P = {p.power_W_cm3 / 1e3:.0f} kW/cm^3"
    groups = ",".join(p.groups) if p.groups else "none"
    lines.append("LIGHT_Kr v0_4: 0-D pure-krypton model")
    lines.append(f"  fill {p.pressure_atm:.3g} atm Kr, Tg0 {p.tg0_K:.0f} K;  "
                 f"{drive} for {p.pulse_s * 1e6:.3g} us;  n_e seed "
                 f"{p.ne_cm3:.2g} cm^-3 (then solved);  t_end "
                 f"{p.t_end_s * 1e6:.3g} us")
    lines.append(f"  network: pure Kr, {len(m.rows)} rows (the ground-state "
                 f"excitations level by level);  optional groups: {groups}")
    lines.append("  " + electrons_line(res))
    lines.append(bands_line(res))
    if p.em_E_V_cm > 0.0:
        lines.append(f"  EM wave: {p.em_E_V_cm:.4g} V/cm rms inside the "
                     f"plasma at {p.em_f_Hz / 1e9:.4g} GHz, absorbed "
                     f"{a['e_em'] * QE:.4g} J/cm^3 "
                     f"({100 * a['e_em'] / e:.2f} % of all input; "
                     f"ib_scale {m.ib_scale:.2f})")
    if p.gap_cm is not None:
        lines.append(f"  gap {p.gap_cm:.3g} cm: the charges drift to the "
                     f"electrodes (walls line in the fate)")
    if p.ladder != "off":
        lines.append(f"  ladder preset: {p.ladder} (every ladder band at "
                     f"Kannari's printed value)")
    lines.append(f"  deposited energy E_dep = {j:.4g} J/cm^3"
                 + (f" (given {a['e_given'] * QE:.4g} + EM "
                    f"{a['e_em'] * QE:.4g})" if a["e_em"] > 0.0 else ""))
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
    if p.gap_cm is not None:
        lines.append(f"  walls     {'charges lost to the electrodes':<34s}"
                     f"{_fmt_pct(a['walls'], e)}")
    lines.append(f"  closure residual (bookkeeping test)   "
                 f"{a['closure']:.1e}")
    lines.append("")
    for name, title in (("Kr*", "Kr* (5s) fate"), ("Kr**", "Kr** (5p) fate"),
                        ("Kr2*", "Kr2* fate")):
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
                 f"p end {res.pressure_atm[-1]:.3g} atm")
    lines.append(f"       n_e end {res.ne[-1]:.3g} cm^-3 "
                 f"(peak {res.ne.max():.3g})"
                 f";  Kr ground left "
                 f"{100 * res.n['Kr'][-1] / res.n['Kr'][0]:.3g} %")
    kr5 = sum(np.asarray(res.n[s]) for s in KR_5S)
    lines.append(f"       peaks: Kr* {kr5.max():.3g}, Kr** "
                 f"{res.n['Kr**'].max():.3g}, Kr2* "
                 f"{res.n['Kr2*'].max():.3g} cm^-3;  RHS evals {res.nfev}")
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
    kr = a["krstar"]
    kr2 = a["kr2"]
    band = a["bands"]
    return {
        "pressure_atm": p.pressure_atm,
        "power_kW_cm3": p.power_W_cm3 / 1e3, "ne_cm3": p.ne_cm3,
        "groups": "+".join(p.groups) or "basic",
        "E_dep_J_cm3": e * QE,
        "to_Krstar": _pct(a["flows"]["to Kr* (e-impact excitation, all "
                                "levels)"], e),
        "to_Krplus": _pct(a["flows"]["to Kr+ (every ionization x ion energy)"],
                          e),
        "to_Kr2star": _pct(kr2["formed"] * a["hv147"], e),
        "heat_reactions": _pct(a["heat_reactions"], e),
        "elastic": _pct(a["heat_elastic"], e),
        "emission_vuv": _pct(band.get(VUV_BAND, 0.0), e),
        "emission_nir": _pct(band.get(NIR_BAND, 0.0), e),
        "emission_all": _pct(sum(band.values()), e),
        "stored_ions": _pct(a["stores"].get("ions (Kr+, Kr2+)", 0.0), e),
        "stored_Krpp": _pct(a["stores"].get("Kr** (5p)", 0.0), e),
        "krstar_to_kr2_pct": 100 * sum(
            v for r, v in kr["losses"].items() if r.startswith("AK8"))
        / kr["formed"] if kr["formed"] else float("nan"),
        "kr2_to_vuv_pct": 100 * (kr2["losses"].get("KR2_RAD", 0.0)
                                 + kr2["losses"].get("KR2F_RAD", 0.0))
        / kr2["formed"] if kr2["formed"] else float("nan"),
        "te_mean_eV": float(np.mean(res.te[res.t <= p.pulse_s])),
        "tg_end_K": float(res.tg[-1]),
        "ne_end_cm3": float(res.ne[-1]),
        "closure": a["closure"],
        "walls": _pct(a["walls"], e), "em_J_cm3": a["e_em"] * QE,
    }
