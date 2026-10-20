"""Paschen curve for argon -- Townsend breakdown, measured, vs Raizer.

THE CRITERION, exactly as the source states it (Raizer, *Gas Discharge
Physics*, Springer 1991, Eq. (7.1), read first-hand 2026-08-14):

    gamma [exp(alpha d) - 1] = 1

The gas-phase half, M0 = exp(alpha d) - 1 = ions created per seed electron,
is MEASURED here by direct PIC-MCC avalanche counting: seed electrons leave
the cathode, multiply through the real LXCat argon chain, and every
ionization event is counted (ions stay as immobile charge -- their return
trip does not change how many there are).  The surface half, gamma, is a
STATED PARAMETER (it belongs to the cathode, not the gas): the breakdown
voltage for a given gamma is the bisection root of gamma * M0(V; pd) = 1.

WHAT THE COMPARISON IS (all read from the book, not recalled):
  * the printed experimental anchor, argon with an IRON cathode:
      (pd)_min = 1.5 Torr cm, V_min = 265 V, (E/p)_min = 176 V/(cm Torr)
      (Raizer Sect. 7.2.4, p. 135)
  * the similarity-law curve Eq. (7.3): V_t = B pd / (C + ln pd),
      C = ln A - ln ln(1/gamma + 1), with Table 4.1 argon constants
      A = 12 /(cm Torr), B = 180 V/(cm Torr), stated valid E/p 100-600;
      that fit + the printed anchor imply gamma_iron ~ 1.5e-3
  * the alpha/p check: alpha d = ln(1 + M0) against BOTH Raizer forms,
      (4.5) alpha/p = A exp(-B/(E/p)) and the inert-gas form
      (4.6) alpha/p = C' exp(-D sqrt(p/E)), C' = 29.2, D = 26.6.

STATED LIMITS: gamma here is energy-independent (ion-neutral collisions,
plan section 2.1 item 3, are open -- ions free-fall, so an ion-energy-
resolved gamma would be dishonest); the M0 measurement is space-charge-free
by construction (unit weights), which IS the Townsend regime; points with
E/p outside Table 4.1's stated 100-600 window are flagged in the output.

Run from the Kinetic_Modeling root:

    python3 SPECK/PIC_MCC_3D/tools/paschen.py --curve          # ~minutes
    python3 SPECK/PIC_MCC_3D/tools/paschen.py --demo           # ~minutes
    python3 SPECK/PIC_MCC_3D/tools/paschen.py --point 1.5 265  # one M0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "PIC_MCC_3D" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("PIC_MCC_3D/MODULES not found above this script")

from PIC_MCC_3D.MODULES import emission3d as EM                # noqa: E402
from PIC_MCC_3D.MODULES import species3d as SP                 # noqa: E402
from PIC_MCC_3D.MODULES.cross_sections import mixture_set      # noqa: E402
from PIC_MCC_3D.MODULES.es_sim import ESSim                    # noqa: E402
from PIC_MCC_3D.MODULES.gas import pure                        # noqa: E402
from PIC_MCC_3D.MODULES.mcc3d import MCC3D                     # noqa: E402
from PIC_MCC_3D.MODULES.units import E_CHARGE, M_E             # noqa: E402

OUT_DEFAULT = ROOT / "picmcc_out" / "paschen"

# Raizer, read first-hand 2026-08-14 (see module docstring for pages)
RAIZER_AR = {"A": 12.0, "B": 180.0, "ep_lo": 100.0, "ep_hi": 600.0,
             "C_inert": 29.2, "D_inert": 26.6,
             "anchor": {"pd_torr_cm": 1.5, "V_min": 265.0,
                        "ep_min": 176.0, "cathode": "iron"}}
GAMMA_IRON = 1.5e-3          # implied by (7.4) + the printed anchor
TORR_PA = 133.322


def raizer_vt(pd_torr_cm, gamma: float) -> np.ndarray:
    """Similarity law Eq. (7.3) with Table 4.1 argon constants."""
    pd = np.asarray(pd_torr_cm, dtype=float)
    C = np.log(RAIZER_AR["A"]) - np.log(np.log(1.0 / gamma + 1.0))
    den = C + np.log(pd)
    v = RAIZER_AR["B"] * pd / den
    return np.where(den > 0, v, np.inf)


def _mcc_for(p_torr: float, seed_extra: int = 0):
    mix = pure("Ar", pressure_atm=p_torr / 760.0, T_gas_K=300.0)
    cs = mixture_set(mix)
    return mix, cs


def measure_M0(pd_torr_cm: float, V: float, d_m: float = 1.0e-2,
               n_seed: int = 32, nz: int = 33, seed: int = 20260814,
               T_seed_eV: float = 2.0, max_steps: int = 120000,
               sharing: str = "equal", max_macro: int = 150000) -> dict:
    """Ions created per cathode seed electron at (pd, V) -- gamma-free.

    ``max_macro`` caps the avalanche: once the live electron count passes
    it the run stops and ``capped`` is set -- M0 is then a LOWER BOUND.
    With n_seed ~ 50 that bound is ~3000, above every breakdown target this
    tool uses (1/gamma <= 1000), so a capped point reads as "above
    threshold" to the bisection.  (Known-bad 2026-08-14: without the cap a
    fixed high bracket voltage at large pd asked for e^37 particles and ran
    13 CPU-hours producing nothing.)"""
    p_torr = pd_torr_cm / (d_m * 100.0)
    mix, cs = _mcc_for(p_torr)
    e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=4096)
    mcc = MCC3D(cs, mixture=mix,          # no ion species: counting mode
                ionization_sharing=sharing)
    h = d_m / (nz - 1)
    v_cap = np.sqrt(2.0 * E_CHARGE * max(V, 50.0) * 1.5 / M_E)
    dt = 0.35 * h / v_cap
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop]), dt,
                phi_lo=0.0, phi_hi=float(V), mcc=mcc, seed=seed)
    rng = np.random.default_rng(seed + 101)
    vx, vy, vz = EM.flux_maxwellian(M_E, T_seed_eV, n_seed, rng, +1.0)
    e_pop.spawn(rng.random(n_seed) * sim.L[0], rng.random(n_seed) * sim.L[1],
                rng.random(n_seed) * np.abs(vz) * dt, vx, vy, vz,
                np.ones(n_seed))
    sim.prepare()
    steps = 0
    peak = e_pop.n
    capped = False
    while e_pop.n and steps < max_steps:
        sim.step()
        steps += 1
        peak = max(peak, e_pop.n)
        if e_pop.n > max_macro:
            capped = True
            break
    ion_k = [k for k in range(cs.n_proc) if cs.kinds[k] == "ionization"]
    n_ion = float(sim.collision_tally[ion_k].sum()) if \
        sim.collision_tally is not None else 0.0
    eps_ceiling = cs.energy_ceiling_eV      # @property
    return {
        "pd_torr_cm": pd_torr_cm, "V": V, "d_m": d_m, "p_torr": p_torr,
        "ep_V_cm_torr": V / (d_m * 100.0) / p_torr,
        "M0": n_ion / n_seed, "n_ion": n_ion, "n_seed": n_seed,
        "alpha_d": float(np.log1p(n_ion / n_seed)),
        "back_absorbed": sim.wall_flux_w.get(("e-", "lo"), 0.0),
        "anode_absorbed": sim.wall_flux_w.get(("e-", "hi"), 0.0),
        "steps": steps, "step_cap_hit": bool(e_pop.n) and not capped,
        "capped": capped,
        "peak_macro": peak, "eps_ceiling_eV": float(eps_ceiling),
    }


def breakdown_voltage(pd_torr_cm: float, gamma: float, v_lo: float = 60.0,
                      v_max: float = 3000.0, grow: float = 1.3,
                      n_bisect: int = 7, **kw) -> dict:
    """Bisect gamma * M0(V) = 1.  Deterministic (fixed seeds per eval).

    The upper bracket is FOUND by a geometric climb (x``grow`` per step)
    from ``v_lo`` rather than assumed: M0 grows like exp(alpha d), so a
    fixed high bracket is an astronomically large avalanche at large pd
    (see measure_M0's Known-bad).  The climb keeps every evaluation within
    a bounded factor of the threshold avalanche."""
    target = 1.0 / gamma
    trace = []

    def excess(v):
        r = measure_M0(pd_torr_cm, v, **kw)
        trace.append((v, r["M0"], r["capped"]))
        if r["capped"]:
            return np.inf, r
        return r["M0"] - target, r

    e_lo, r_lo = excess(v_lo)
    if e_lo > 0:
        return {"pd_torr_cm": pd_torr_cm, "gamma": gamma, "V_b": v_lo,
                "bracket": "below", "trace": trace, "last": r_lo}
    v_a, v_b = v_lo, v_lo * grow
    last = r_lo
    while True:
        e_b, last = excess(v_b)
        if e_b >= 0:
            break
        v_a, v_b = v_b, v_b * grow
        if v_b > v_max:
            return {"pd_torr_cm": pd_torr_cm, "gamma": gamma, "V_b": np.inf,
                    "bracket": "above", "trace": trace, "last": last}
    v_lo, v_hi = v_a, v_b
    for _ in range(n_bisect):
        v_mid = 0.5 * (v_lo + v_hi)
        e_mid, last = excess(v_mid)
        if e_mid < 0:
            v_lo = v_mid
        else:
            v_hi = v_mid
    vb = 0.5 * (v_lo + v_hi)
    return {"pd_torr_cm": pd_torr_cm, "gamma": gamma, "V_b": vb,
            "bracket": "ok", "trace": trace, "last": last,
            "ep_V_cm_torr": vb / pd_torr_cm,
            "in_table41_window": bool(RAIZER_AR["ep_lo"] <= vb / pd_torr_cm
                                      <= RAIZER_AR["ep_hi"])}


def paschen_curve(pds=(0.7, 0.9, 1.2, 1.5, 2.5, 4.0, 8.0),
                  gammas=(GAMMA_IRON, 1.0e-2), **kw) -> dict:
    out = {}
    for g in gammas:
        rows = []
        for pd in pds:
            r = breakdown_voltage(pd, g, **kw)
            rows.append(r)
            print(f"  gamma={g:g} pd={pd:5.2f} -> V_b={r['V_b']:7.1f} V "
                  f"({r['bracket']}, E/p={r.get('ep_V_cm_torr', np.nan):6.1f}"
                  f"{'' if r.get('in_table41_window', True) else ' OUTSIDE 100-600'})")
        out[g] = rows
    return out


def alpha_two_length(ep_V_cm_torr: float, p_torr: float = 1.5,
                     d1_cm: float = 0.6, d2_cm: float = 1.2,
                     n_seed: int = 64, shrink: float = 0.6,
                     max_tries: int = 5, **kw) -> dict:
    """alpha/p at one E/p by THE SOURCE'S OWN METHOD (Raizer Sect. 4.1.4):
    two gap lengths at fixed E/p, alpha = slope of ln M0 vs d.  The
    'const' -- seed spectrum, back-diffusion, the equilibration length --
    cancels in the difference; a single-length ln(1+M0)/(pd) does NOT
    cancel it and under-reads alpha (measured x0.68 at the anchor point
    before this fix, 2026-08-14).

    CAP-AWARE: if either gap's avalanche hits measure_M0's macro cap the
    slope would be a lower bound, so both gaps shrink by ``shrink`` and the
    pair is re-run (the lengths actually used are reported).  Without this
    the E/p = 220 and 300 points read x0.81 and x0.47 -- artifacts of a
    truncated long gap, not physics (found 2026-08-14)."""
    d1, d2 = float(d1_cm), float(d2_cm)
    rows = None
    for _ in range(max_tries):
        rows = []
        for d_cm in (d1, d2):
            pd = p_torr * d_cm
            V = ep_V_cm_torr * pd
            rows.append(measure_M0(pd, V, d_m=d_cm / 100.0,
                                   n_seed=n_seed, **kw))
        if not any(r["capped"] for r in rows):
            break
        d1, d2 = d1 * shrink, d2 * shrink
    m1, m2 = rows[0]["M0"], rows[1]["M0"]
    alpha = (np.log(m2) - np.log(m1)) / (d2 - d1)
    return {"ep": ep_V_cm_torr, "alpha_p_pic": alpha / p_torr,
            "M0_short": m1, "M0_long": m2, "p_torr": p_torr,
            "d1_cm": d1, "d2_cm": d2,
            "capped": bool(any(r["capped"] for r in rows)),
            "alpha_p_45": RAIZER_AR["A"]
            * np.exp(-RAIZER_AR["B"] / ep_V_cm_torr),
            "alpha_p_46": RAIZER_AR["C_inert"]
            * np.exp(-RAIZER_AR["D_inert"] / np.sqrt(ep_V_cm_torr))}


def alpha_scan(eps=(100.0, 120.0, 150.0, 176.0, 220.0, 300.0), **kw) -> list:
    """alpha/p vs E/p by the two-length method, against both Raizer fits."""
    rows = []
    for ep in eps:
        r = alpha_two_length(ep, **kw)
        rows.append(r)
        print(f"  E/p={ep:6.1f}  alpha/p PIC={r['alpha_p_pic']:.3f}  "
              f"(4.5)={r['alpha_p_45']:.3f}  (4.6)={r['alpha_p_46']:.3f}  "
              f"ratio vs (4.5) x{r['alpha_p_pic']/r['alpha_p_45']:.2f}  "
              f"[gaps {r['d1_cm']:.2f}/{r['d2_cm']:.2f} cm"
              f"{', CAPPED' if r['capped'] else ''}]")
    return rows


def selfsustain_demo(gamma: float = 0.1, pd_torr_cm: float = 1.0,
                     d_m: float = 5.0e-3, nz: int = 33,
                     v_over: float = 1.3, v_under: float = 0.7,
                     cycles: float = 3.0, n_seed: int = 30,
                     seed: int = 20260814) -> dict:
    """The gamma loop CLOSED dynamically: real Ar+ macro-ions free-fall back
    to the cathode and knock out real secondaries (M5).  Above threshold the
    electron population must grow cycle over cycle; below, die.

    gamma=0.1 is a stated demo value (a poor/dirty cathode class), chosen so
    one multiplication cycle is small and the run stays short; the physics
    statement is growth vs decay, not the iron-anchor voltage.
    """
    vb = breakdown_voltage(pd_torr_cm, gamma, d_m=d_m, nz=nz,
                           n_seed=24, seed=seed)["V_b"]
    p_torr = pd_torr_cm / (d_m * 100.0)
    runs = {}
    for tag, fac in (("over", v_over), ("under", v_under)):
        V = fac * vb
        mix, cs = _mcc_for(p_torr)
        e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=65536)
        ar_pop = SP.ParticleSpecies(SP.Species.ion("Ar", 39.948),
                                    capacity=65536)
        ar_pop.subcycle = 20
        mcc = MCC3D(cs, mixture=mix, ion_species={"Ar": ar_pop})
        h = d_m / (nz - 1)
        v_cap = np.sqrt(2.0 * E_CHARGE * V * 1.5 / M_E)
        dt = 0.35 * h / v_cap
        sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop, ar_pop]), dt,
                    phi_lo=0.0, phi_hi=float(V), mcc=mcc, seed=seed)
        see = EM.SecondaryEmitter(e_pop, ar_pop, gamma=gamma, wall="lo",
                                  T_emit_eV=2.0, w_macro=1.0)
        sim.wall_hooks.append(see.wall_hit)
        rng = np.random.default_rng(seed + 7)
        vx, vy, vz = EM.flux_maxwellian(M_E, 2.0, n_seed, rng, +1.0)
        e_pop.spawn(rng.random(n_seed) * sim.L[0],
                    rng.random(n_seed) * sim.L[1],
                    rng.random(n_seed) * np.abs(vz) * dt,
                    vx, vy, vz, np.ones(n_seed))
        # one gamma cycle ~ the ion free-fall transit
        E_field = V / d_m
        a_ion = E_CHARGE * E_field / (39.948 * 1.66053906660e-27)
        t_cycle = np.sqrt(2.0 * d_m / a_ion)
        n_steps = int(cycles * t_cycle / dt)
        t_hist, ne_hist, ni_hist, sec_hist = [], [], [], []
        sim.prepare()
        for k in range(n_steps):
            sim.step()
            if k % 200 == 0:
                t_hist.append(sim.time)
                ne_hist.append(e_pop.real_count())
                ni_hist.append(ar_pop.real_count())
                sec_hist.append(see.acc.emitted_real)
        runs[tag] = {"V": V, "t_s": np.array(t_hist),
                     "ne": np.array(ne_hist), "ni": np.array(ni_hist),
                     "secondaries": np.array(sec_hist),
                     "t_cycle_s": t_cycle, "dt": dt, "steps": n_steps}
    over, under = runs["over"], runs["under"]
    return {"gamma": gamma, "pd_torr_cm": pd_torr_cm, "V_b": vb,
            "over": over, "under": under,
            "grew": bool(over["ne"][-1] > 3.0 * n_seed),
            "died": bool(under["ne"][-1] < 0.5 * n_seed)}


# ----------------------------------------------------------------- artifacts
def _npz_path(out_dir: Path) -> Path:
    return out_dir / "paschen_data.npz"


def _load_npz(out_dir: Path) -> dict:
    pth = _npz_path(out_dir)
    if not pth.exists():
        return {}
    with np.load(pth) as z:
        return {k: z[k] for k in z.files}


def save_npz(out_dir: Path, curves: dict = None, alphas: list = None,
             demo: dict = None) -> Path:
    """Merge into the existing npz (never drops keys another run wrote)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    save = _load_npz(out_dir)
    for g, rows in (curves or {}).items():
        ok = [(r["pd_torr_cm"], r["V_b"]) for r in rows
              if np.isfinite(r["V_b"]) and r["bracket"] == "ok"]
        if ok:
            save[f"pd_gamma_{g:g}"] = np.array([x for x, _ in ok])
            save[f"vb_gamma_{g:g}"] = np.array([y for _, y in ok])
    if alphas:
        save["alpha_ep"] = np.array([r["ep"] for r in alphas])
        save["alpha_p_pic"] = np.array([r["alpha_p_pic"] for r in alphas])
        save["alpha_d1_cm"] = np.array([r["d1_cm"] for r in alphas])
        save["alpha_d2_cm"] = np.array([r["d2_cm"] for r in alphas])
    if demo:
        for tag in ("over", "under"):
            save[f"demo_{tag}_t"] = demo[tag]["t_s"]
            save[f"demo_{tag}_ne"] = demo[tag]["ne"]
            save[f"demo_{tag}_V"] = np.array([demo[tag]["V"]])
            save[f"demo_{tag}_tcycle"] = np.array([demo[tag]["t_cycle_s"]])
        save["demo_meta"] = np.array([demo["gamma"], demo["pd_torr_cm"],
                                      demo["V_b"]])
    np.savez_compressed(_npz_path(out_dir), **save)
    return _npz_path(out_dir)


def curves_from_npz(out_dir: Path) -> dict:
    z = _load_npz(out_dir)
    curves = {}
    for k in z:
        if k.startswith("pd_gamma_"):
            g = float(k[len("pd_gamma_"):])
            vb = z["vb_gamma_" + k[len("pd_gamma_"):]]
            curves[g] = [{"pd_torr_cm": float(a), "V_b": float(b),
                          "bracket": "ok"} for a, b in zip(z[k], vb)]
    return curves


def alphas_from_npz(out_dir: Path) -> list:
    z = _load_npz(out_dir)
    if "alpha_ep" not in z:
        return []
    out = []
    for i, ep in enumerate(z["alpha_ep"]):
        ep = float(ep)
        out.append({"ep": ep, "alpha_p_pic": float(z["alpha_p_pic"][i]),
                    "alpha_p_45": RAIZER_AR["A"]
                    * np.exp(-RAIZER_AR["B"] / ep),
                    "alpha_p_46": RAIZER_AR["C_inert"]
                    * np.exp(-RAIZER_AR["D_inert"] / np.sqrt(ep)),
                    "d1_cm": float(z["alpha_d1_cm"][i])
                    if "alpha_d1_cm" in z else np.nan,
                    "d2_cm": float(z["alpha_d2_cm"][i])
                    if "alpha_d2_cm" in z else np.nan})
    return out


def render_curve(curves: dict, out_dir: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    pp = np.geomspace(0.55, 12.0, 200)
    marks = ["o", "s", "^"]
    for i, (g, rows) in enumerate(sorted(curves.items())):
        pd_v = [(r["pd_torr_cm"], r["V_b"]) for r in rows
                if np.isfinite(r["V_b"]) and r["bracket"] == "ok"]
        if pd_v:
            x, y = zip(*pd_v)
            ax.plot(x, y, marks[i % 3], ms=6, label=f"PIC-MCC, gamma={g:g}")
        ax.plot(pp, raizer_vt(pp, g), "--", lw=1.0,
                label=f"Raizer (7.3), A=12 B=180, gamma={g:g}")
    a = RAIZER_AR["anchor"]
    ax.plot([a["pd_torr_cm"]], [a["V_min"]], "k*", ms=13,
            label=f"measured, Ar/{a['cathode']} (Raizer p.135)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pd [Torr cm]")
    ax.set_ylabel("V_b [V]")
    ax.set_title("Paschen curve, argon -- gamma is the cathode's, stated")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    pth = out_dir / "paschen_ar.png"
    fig.savefig(pth, dpi=140)
    plt.close(fig)
    return pth


def render_alpha(alphas: list, out_dir: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ep = np.array([r["ep"] for r in alphas])
    ax.semilogy(ep, [r["alpha_p_pic"] for r in alphas], "o", ms=6,
                label="PIC-MCC, two-length slope of ln M0 (Raizer 4.1.4)")
    ee = np.linspace(min(ep) * 0.9, max(ep) * 1.08, 120)
    ax.semilogy(ee, RAIZER_AR["A"] * np.exp(-RAIZER_AR["B"] / ee), "--",
                lw=1.0, label="Raizer (4.5): 12 exp(-180 p/E)")
    ax.semilogy(ee, RAIZER_AR["C_inert"]
                * np.exp(-RAIZER_AR["D_inert"] / np.sqrt(ee)), ":",
                lw=1.2, label="Raizer (4.6): 29.2 exp(-26.6 sqrt(p/E))")
    ax.set_xlabel("E/p [V/(cm Torr)]")
    ax.set_ylabel("alpha/p [1/(cm Torr)]")
    ax.set_title("first Townsend coefficient, argon")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    pth = out_dir / "alpha_ar.png"
    fig.savefig(pth, dpi=140)
    plt.close(fig)
    return pth


def render_demo(demo: dict, out_dir: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for tag, style, lab in (("over", "-", "above threshold"),
                            ("under", "--", "below threshold")):
        r = demo[tag]
        ax.semilogy(r["t_s"] / r["t_cycle_s"], r["ne"], style, lw=1.6,
                    label=f"{lab} (V={r['V']:.0f} V)")
    ax.set_xlabel("t / ion-transit cycles")
    ax.set_ylabel("real electrons in gap")
    ax.set_title(f"gamma loop closed dynamically (gamma={demo['gamma']}, "
                 f"V_b={demo['V_b']:.0f} V, pd={demo['pd_torr_cm']})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    pth = out_dir / "selfsustain_demo.png"
    fig.savefig(pth, dpi=140)
    plt.close(fig)
    return pth


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--curve", action="store_true")
    ap.add_argument("--alpha", action="store_true")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="re-render curve/alpha figures from the saved npz")
    ap.add_argument("--point", nargs=2, type=float, metavar=("PD", "V"))
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    if args.point:
        r = measure_M0(args.point[0], args.point[1])
        print(f"pd={r['pd_torr_cm']} V={r['V']}  M0={r['M0']:.2f} "
              f"(alpha d={r['alpha_d']:.3f}, E/p={r['ep_V_cm_torr']:.1f}, "
              f"back={r['back_absorbed']:.0f}/{r['n_seed']}, "
              f"capped={r['capped']})")
    made = []
    if args.curve:
        curves = paschen_curve()
        save_npz(args.out, curves=curves)
        made.append(render_curve(curves, args.out))
    if args.alpha:
        alphas = alpha_scan()
        save_npz(args.out, alphas=alphas)
        made.append(render_alpha(alphas, args.out))
    if args.demo:
        demo = selfsustain_demo()
        print(f"demo: V_b={demo['V_b']:.0f} V, grew={demo['grew']}, "
              f"died={demo['died']}")
        save_npz(args.out, demo=demo)
        made.append(render_demo(demo, args.out))
    if args.replot:
        c = curves_from_npz(args.out)
        if c:
            made.append(render_curve(c, args.out))
        al = alphas_from_npz(args.out)
        if al:
            made.append(render_alpha(al, args.out))
    for pth in made:
        print("wrote", pth)
    if made:
        print("data of record:", _npz_path(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
