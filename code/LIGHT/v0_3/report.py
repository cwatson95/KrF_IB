"""Energy accounting and the text report (the notes' 'Output' line).

Two views of the same ledger, both in percent of the deposited energy:

FLOWS — what passed through each channel ("percent to Kr+, KrF(B), Kr*,
heat, elastic losses and emission losses"). Flows overlap: energy that
reaches Kr* goes on to KrF(B) and then to emission, so they do not sum.

FATE AT t_end — where the energy sits when the run stops: radiated by band,
heat (elastic and reaction, separately), or still stored in a species or
the electron gas. These sum to 100 % up to the printed closure residual,
which is the model's bookkeeping test (events x per-event energies vs the
densities the ODE produced).

Reservoir accounting treats KrF(B) + KrF(C) as ONE manifold (REBUILD 0.1.24
lesson: B<->C mixing rows swamp per-state bookkeeping).
"""

from __future__ import annotations

import numpy as np

from .model import QE
from .network import BANDS
from .species import e_form

RESERVOIRS = {
    "Kr*": ("Kr*",),
    "KrF*": ("KrF(B)", "KrF(C)"),
    "Kr2F*": ("Kr2F*",),
    "Kr2*": ("Kr2*",),
}
STORE_GROUPS = (
    ("ions (Kr+, Kr2+, F-)", ("Kr+", "Kr2+", "F-")),
    ("Kr*", ("Kr*",)),
    ("Kr2*", ("Kr2*",)),
    ("KrF(B)+KrF(C)", ("KrF(B)", "KrF(C)")),
    ("Kr2F*", ("Kr2F*",)),
    ("F atoms (F2 broken)", ("F",)),
)
G5_TG_VALID_K = 800.0     # Shui's calc 300-800 K; measured 294-425 K
KR2_TG_VALID_K = 451.0    # Kr* + 2Kr's Tg law, Wieme & Lenaerts 300-451 K
HOT_K = 1000.0            # redissociation of Kr2+ and Kr2F* not modeled
VUV_BAND = "147 nm (Kr2* VUV)"


def _pct(x, e_dep):
    return 100.0 * x / e_dep if e_dep > 0 else float("nan")


def accounting(res) -> dict:
    """All energies in eV/cm^3 at t_end, plus reservoir branchings."""
    m = res.model
    e_dep = float(res.ledger["E_dep"][-1])
    ev = {rid: float(v[-1]) for rid, v in res.events.items()}
    rows = {c.row.rid: c for c in m.compiled}

    bands = {}
    heat_rx = float(res.ledger["heat_capture"][-1])
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
    def evs(*rids):
        return sum(ev.get(r, 0.0) for r in rids)

    exc_rids = [r for r in rows if r == "G1" or r.startswith("C1")]
    e_exc = sum(ev[r] * rows[r].route["from_e"] for r in exc_rids)
    e_exc_5s = ev["G1"] * rows["G1"].route["from_e"]
    # every row that frees an electron makes an ion; book the cation's
    # formation energy (Kr+ 13.996 eV, Kr2+ 12.846 eV) per event
    e_ion_all = sum(ev[r] * sum(e_form(s) for s in c.row.products
                                if s.endswith("+"))
                    for r, c in rows.items() if c.row.ne_delta > 0)
    e_ion_e = ev.get("C2", 0.0) * e_form("Kr+")
    krf = reservoir(res, "KrF*")
    krs = reservoir(res, "Kr*")
    hv248 = rows["G4"].route["rad"]
    e_to_e = sum(ev[r] * c.route["to_e"] for r, c in rows.items())
    flows = {
        "to Kr* (e-impact excitation, all levels)": e_exc,
        "   of which 5s direct (G1)": e_exc_5s,
        "   of which 5p/4d/higher -> 5s (C1)": e_exc - e_exc_5s,
        "to Kr+ (every ionization x ion energy)": e_ion_all,
        "   of which e-impact on Kr (C2)": e_ion_e,
        "to KrF(B) (formed x 5.00 eV quantum)": krf["formed"] * hv248,
        "heat (reactions + relaxation)": heat_rx,
        "elastic losses (G6)": heat_el,
        "emission losses (all bands)": sum(bands.values()),
        "returned to electrons (Penning etc.)": e_to_e,
    }
    return {"e_dep": e_dep, "bands": bands, "heat_elastic": heat_el,
            "heat_reactions": heat_rx, "stores": stores, "closure": closure,
            "flows": flows, "krf": krf, "krstar": krs, "events": ev}


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


def f2_burnout_time(res, frac=0.1):
    """First time [s] F2 falls below ``frac`` of its initial density."""
    f2 = res.n["F2"]
    below = np.nonzero(f2 < frac * f2[0])[0]
    return float(res.t[below[0]]) if below.size else None


def warnings(res) -> list:
    p, out = res.p, []
    from . import xsec
    if res.te.max() > xsec.TE_MAX:
        out.append(f"Te reached {res.te.max():.1f} eV, above the rate-table "
                   f"edge {xsec.TE_MAX:.0f} eV: rates were held at the edge")
    if res.tg.max() > G5_TG_VALID_K:
        out.append(f"Tg reached {res.tg.max():.0f} K. The trimer's "
                   "(Tg/300)^-3.32 follows Shui's calculation for 300-800 K "
                   "(measured 294-425 K). Above that it is an "
                   "extrapolation, and Kr2F*, bound by 0.58 eV, may "
                   "dissociate back to KrF* + Kr, for which no rate is "
                   "published. Pure-Kr efficiencies suggest a weaker "
                   "T^-1/2 (try --trimer-tg-pow -0.5 for the band)")
    ak8 = [c for c in res.model.compiled if c.row.rid == "AK8_Kr"]
    if ak8 and ak8[0].rate > 0 and res.tg.max() > KR2_TG_VALID_K:
        out.append(f"Tg reached {res.tg.max():.0f} K. Kr* + 2Kr's "
                   "(Tg/300)^-0.64 is measured over 300-451 K (Wieme & "
                   "Lenaerts 1981) and extrapolated above it")
    if res.tg.max() > HOT_K:
        out.append(f"Tg above {HOT_K:.0f} K: Kr2+ (bound by 1.15 eV) and "
                   "Kr2F* (0.58 eV) redissociate on these time scales, "
                   "which the model does not include, and the ion "
                   "conversion and ion-ion laws are extrapolated beyond "
                   "their data")
    tb = f2_burnout_time(res)
    if tb is not None:
        out.append(f"F2 burn-out: below 10 % of the fill at t = "
                   f"{tb * 1e9:.0f} ns. F + F + Kr re-forms F2 on a "
                   "millisecond scale, far too slowly to refill it")
    kr_loss = 1.0 - res.n["Kr"][-1] / res.n["Kr"][0]
    if kr_loss > 0.05:
        out.append(f"ground-state Kr depleted by {100 * kr_loss:.1f} %")
    if res.model.ne_mode == "clamped" and any(
            c.row.ne_delta < 0 for c in res.model.compiled):
        out.append("n_e is held while attachment or recombination removes "
                   "electrons, so charge is not conserved (CHANGELOG "
                   "Known-bad). Use the default free n_e for results")
    if res.model.ne_mode == "clamped":
        ions = sum(res.n[s][-1] for s in ("Kr+", "Kr2+") if s in res.n)
        ratio = ions / p.ne_cm3
        if ratio > 1.5:
            out.append(f"clamped n_e: the ion inventory is {ratio:.1f} x n_e "
                       "at t_end. Without the 'ions' group nothing "
                       "recombines, so ions are bookkeeping only")
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
    if m.eedf != "kinetic":
        return f"electrons: {m.eedf} EEDF, captured energy {m.capture}"
    log = [e for e in getattr(res, "kinetic", []) if not e.get("skipped")]
    on = [e for e in log if e["on"]]
    txt = (f"electrons: kinetic EEDF (two-term Boltzmann with e-e, "
           f"boltzmann.py), {len(log)} solutions, {len(on)} in the pulse")
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


def summary(res) -> str:
    p, a = res.p, accounting(res)
    e = a["e_dep"]
    j = e * QE
    lines = []
    drive = f"P = {p.power_W_cm3 / 1e3:.0f} kW/cm^3"
    groups = ",".join(p.groups) if p.groups else "none"
    lines.append("LIGHT: 0-D Kr/F2 KrF model (KrF_model_notes_with_Gary)")
    lines.append(f"  fill {p.pressure_atm:.3g} atm, F2 {100 * p.x_f2:.3g} % "
                 f"in Kr, Tg0 {p.tg0_K:.0f} K;  {drive} for "
                 f"{p.pulse_s * 1e6:.3g} us;  n_e {p.ne_cm3:.2g} cm^-3 "
                 f"({res.model.ne_mode});  t_end {p.t_end_s * 1e6:.3g} us")
    lines.append(f"  network {p.network};  tau_sp {p.tau_sp_s * 1e9:.2f} ns;"
                 f"  optional groups: "
                 f"{groups}")
    lines.append("  " + electrons_line(res))
    read = {c.row.rate[1] for c in res.model.compiled
            if c.row.rate[0] == "param"}
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
    for name, title in (("Kr*", "Kr* (5s) fate"), ("KrF*", "KrF(B) fate")):
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
                 f";  Kr ground left "
                 f"{100 * res.n['Kr'][-1] / res.n['Kr'][0]:.3g} %")
    lines.append(f"       peaks: Kr* {res.n['Kr*'].max():.3g}, KrF(B) "
                 f"{res.n['KrF(B)'].max():.3g}, Kr2F* "
                 f"{res.n['Kr2F*'].max():.3g} cm^-3;  RHS evals {res.nfev}")
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
    krf = a["krf"]
    kr = a["krstar"]
    band = a["bands"]
    d = {
        "pressure_atm": p.pressure_atm, "x_f2": p.x_f2,
        "power_kW_cm3": p.power_W_cm3 / 1e3, "ne_cm3": p.ne_cm3,
        "ne_mode": res.model.ne_mode,
        "groups": "+".join(p.groups) or p.network,
        "E_dep_J_cm3": e * QE,
        "to_Krstar": _pct(a["flows"]["to Kr* (e-impact excitation, all "
                                "levels)"], e),
        "to_Krplus": _pct(a["flows"]["to Kr+ (every ionization x ion energy)"],
                          e),
        "to_KrFB": _pct(krf["formed"] * 5.0, e),
        "heat_reactions": _pct(a["heat_reactions"], e),
        "elastic": _pct(a["heat_elastic"], e),
        "emission_248": _pct(band.get("248 nm (KrF B-X)", 0.0), e),
        "emission_vuv": _pct(band.get(VUV_BAND, 0.0), e),
        "emission_all": _pct(sum(band.values()), e),
        "stored_ions": _pct(a["stores"].get("ions (Kr+, Kr2+, F-)", 0.0), e),
        "stored_F": _pct(a["stores"].get("F atoms (F2 broken)", 0.0), e),
        "krf_to_248_pct": 100 * krf["losses"].get("G4", 0.0) / krf["formed"]
        if krf["formed"] else float("nan"),
        "krstar_to_krf_pct": 100 * kr["losses"].get("G3", 0.0) / kr["formed"]
        if kr["formed"] else float("nan"),
        "te_mean_eV": float(np.mean(res.te[res.t <= p.pulse_s])),
        "tg_end_K": float(res.tg[-1]),
        "f2_left_pct": 100 * res.n["F2"][-1] / res.n["F2"][0],
        "closure": a["closure"],
    }
    return d
