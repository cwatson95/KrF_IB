"""The cross-section atlas — why the swarm EEDF has the shape it has, per
gas, and what the choice of published database is worth in rate
coefficients (0.1.38).

Owner ask 2026-09-03, after the MAXWELL verification found the swarm shape
the better beam-mode closure: "look at graphing elastic cross section
against eV with the different gases to optimize our rate reaction
coefficients".

POLICY LINE, stated once. A cross section is never adjusted to improve an
anchor — that is the D4 calibration path this package forbids
structurally. What this tool does instead:

  (a) DIAGNOSE. sigma_el(eps) sets both how a swarm electron heats (the
      momentum-transfer frequency) and how it cools (the (2m/M) elastic
      energy-loss rate). The Ramsauer-Townsend window in Ar and Kr is why
      those fills run a hot mean energy at a tiny field, and it sits
      exactly where the beam-mode bulk lives.
  (b) BRACKET. Where LXCat ships more than one evaluation of the same gas
      (F2: SIGLO vs TRINITI; Ne: Biagi vs BSR; Kr, Xe), the spread in
      Maxwellian-averaged k(Te) is a real, citable uncertainty band on the
      engine's rates. Choosing a set is a PROVENANCE decision recorded
      with its citation; it is not a fit.
  (c) LOCATE THE HOLE. The swarm tables stop at 0.5-1 Td; the engine's
      beam-mode Te spends most of a pulse below their floor (the [0.1.37]
      te_in_table_frac gauge: 0.8-33%). This tool marks the floor on
      every axis so the extrapolated region is visible, and prices what
      the cross sections say happens there.

Figures (PNG, into rebuild_out/sigma_atlas/):
  fig1_elastic_sigma      sigma_el(eps) for every gas, Ramsauer minima and
                          first thresholds marked, table floors shaded
  fig2_relaxation         momentum-transfer time and ENERGY-relaxation time
                          (2m/M weighted) at the device fill density
  fig3_channels           per gas: elastic vs total inelastic vs attachment
  fig4_databases          same gas, different published evaluations
  fig5_rate_band          Maxwellian-averaged k(Te) per database and the
                          spread band on each engine bucket
  fig6_mixture            the real fills: mixture cooling rate vs eps with
                          the measured operating points overlaid

    python3 SPECK/REBUILD/tools/sigma_atlas.py
    python3 SPECK/REBUILD/tools/sigma_atlas.py --gases Ar Kr --no-figs

Numbers land in report.json. Scores the data; steers nothing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.MODULES.cross_sections import (gas_file, lxcat_roots,  # noqa: E402
                                            parse_lxcat)
from REBUILD.MODULES.loki_tier import (RUN_NAMES_AR_FILL,  # noqa: E402
                                       RUN_NAMES_ARF_CHEN,
                                       RUN_NAMES_NE_FILL, LokiRates)
from REBUILD.MODULES.maxwell_tier import (EPS_GRID, MaxwellRates,  # noqa: E402
                                          maxwell_rate)
from REBUILD.MODULES.units import EV, M_E, gas_density_m3  # noqa: E402

GASES = ("Ne", "Ar", "Kr", "F2", "He", "Xe")
COLOR = {"Ne": "#e6550d", "Ar": "#3182bd", "Kr": "#31a354", "F2": "#756bb1",
         "He": "#636363", "Xe": "#d62728"}

#: alternate published evaluations of the SAME gas, relative to the LoKI-MC
#: input tree (the vendored sets the swarm tables were computed against).
#: Each entry: label -> path fragment. The engine's own choice is resolved
#: by cross_sections.gas_file and always appears as "engine".
ALTERNATES = {
    "F2": {"SIGLO": "Fluorine2/F2_LXCat_SIGLO.txt",
           "TRINITI": "Fluorine2/F2_LXCat_TRINITI.txt"},
    "Ne": {"Biagi-8.9": "Neon/Ne_LXCat.txt",
           "BSR": "Neon/Ne_LXCat_mt.txt"},
    "Kr": {"Biagi-7.1": "Krypton/Kr_LXCat.txt",
           "lokimc": "Krypton/Kr_LXCat_lokimc.txt"},
    "Ar": {"IST-Lisbon": "Argon/Ar_LXCat.txt"},
    "Xe": {"BSR": "Xenon/Xe_LXCat.txt"},
}

#: the fills the engine actually runs (mole fractions) -> its LOKI table set
FILLS = {
    "Ar/Kr/F2 (krf_beam_0d)": ({"Ar": 0.949, "Kr": 0.05, "F2": 0.001},
                               RUN_NAMES_AR_FILL, 3.0),
    "Ne/Kr/F2 (krf_ne_0d)": ({"Ne": 0.949, "Kr": 0.05, "F2": 0.001},
                             RUN_NAMES_NE_FILL, 3.0),
    "Ne/Ar/F2 (arf_chen)": ({"Ne": 0.964, "Ar": 0.035, "F2": 0.001},
                            RUN_NAMES_ARF_CHEN, 3.4),
}

#: measured beam-mode operating points ([0.1.37] maxwell_verify, LOKI tier):
#: (label, Te peak [eV], Te settled [eV])
OPERATING = [("krf_beam_0d", 2.205, 0.043), ("krf_ne_0d", 2.460, 0.043),
             ("ne_optimized", 2.278, 0.630)]


# --------------------------------------------------------------- data access
def _loki_input_root() -> Path | None:
    for r in lxcat_roots():
        if r.name == "Input" or (r / "Fluorine2").is_dir():
            return r
        cand = r.parent
        if (cand / "Fluorine2").is_dir():
            return cand
    for base in [ROOT, *ROOT.parents]:
        cand = base / "LoKI-MC" / "Code" / "Input"
        if cand.is_dir():
            return cand
    return None


#: gases the ENGINE has no mapping for, readable from the vendored tree so
#: the atlas can still show them (the engine's own skip is unchanged)
ATLAS_ONLY = {"Xe": "Xenon/Xe_LXCat.txt"}


def load_set(gas: str, path: Path = None) -> dict:
    """-> {elastic: Process, inelastic: [Process], attachment: [Process],
    ionization: [Process], path, database}."""
    if path is not None:
        p = Path(path)
    else:
        try:
            p = gas_file(gas)
        except KeyError:
            root = _loki_input_root()
            frag = ATLAS_ONLY.get(gas)
            if root is None or frag is None:
                raise
            p = root / frag
    procs = parse_lxcat(p)
    for x in procs:
        x.gas = gas
    el = [x for x in procs if x.kind == "elastic"]
    eff = [x for x in procs if x.kind == "effective"]
    db = "unstated"
    for line in p.read_text().splitlines()[:80]:
        s = line.strip()
        if s.startswith("DATABASE:"):
            db = s.split(":", 1)[1].strip()
            break
        if s.startswith("Source:"):
            db = s.split(":", 1)[1].strip()[:60]
            break
    return {"path": str(p), "database": db,
            "elastic": (el or eff)[0] if (el or eff) else None,
            "elastic_is_effective": not el and bool(eff),
            "excitation": [x for x in procs if x.kind == "excitation"],
            "attachment": [x for x in procs if x.kind == "attachment"],
            "ionization": [x for x in procs if x.kind == "ionization"]}


def sigma_on(proc, eps):
    return np.interp(eps, proc.energy, proc.sigma, left=0.0, right=0.0)


# ------------------------------------------------------------ derived physics
def ramsauer(proc, eps_max: float = 5.0) -> dict:
    """Locate the elastic minimum below eps_max and its depth vs the 5-25 eV
    peak. A 'window' is a real Ramsauer-Townsend minimum: the depth ratio
    is the number that matters, not the presence of a local dip."""
    m = (proc.energy > 1e-3) & (proc.energy < eps_max)
    if not m.any():
        return {}
    i = int(np.argmin(proc.sigma[m]))
    e_min, s_min = float(proc.energy[m][i]), float(proc.sigma[m][i])
    hi = proc.sigma[(proc.energy > 5.0) & (proc.energy < 25.0)]
    s_hi = float(hi.max()) if hi.size else float("nan")
    depth = (s_hi / s_min) if s_min > 0 else float("nan")
    # STATED criterion for calling it a Ramsauer-Townsend window: a genuine
    # interior minimum at least x20 below the 5-25 eV peak. Ne fails it
    # (depth ~6, minimum at the bottom of the tabulation): neon has no RT
    # minimum — its zero-energy scattering length is POSITIVE (A = +0.21 a0)
    # and a MERT zero needs A < 0 (PIC/Neon note, 2026-09-02).
    return {"eps_min_eV": e_min, "sigma_min_m2": s_min,
            "sigma_peak_m2": s_hi, "depth": depth,
            "at_table_edge": bool(i == 0),
            "is_rt_window": bool(np.isfinite(depth) and depth >= 20.0
                                 and i != 0)}


def relaxation_times(st: dict, n_gas_m3: float, eps=EPS_GRID) -> dict:
    """Momentum-transfer and ENERGY-relaxation times [s] vs eps for one gas
    at number density n_gas_m3.

    nu_m = n sigma_el v; the elastic ENERGY-loss rate carries the mass
    ratio: nu_u = (2 m_e / M) nu_m. tau_u = 1/nu_u is how long a bulk
    electron takes to give its energy to the gas — the quantity the
    Ramsauer window makes enormous, and the reason a low-field swarm mean
    energy stays far above the gas temperature."""
    el = st["elastic"]
    v = np.sqrt(2.0 * eps * EV / M_E)
    nu_m = n_gas_m3 * sigma_on(el, eps) * v
    ratio = float(el.mass_ratio) if el.mass_ratio > 0 else float("nan")
    nu_u = 2.0 * ratio * nu_m
    with np.errstate(divide="ignore"):
        return {"eps_eV": eps, "nu_m": nu_m, "nu_u": nu_u,
                "tau_m_s": 1.0 / np.where(nu_m > 0, nu_m, np.nan),
                "tau_u_s": 1.0 / np.where(nu_u > 0, nu_u, np.nan),
                "mass_ratio": ratio}


def inelastic_total(st: dict, eps=EPS_GRID) -> dict:
    out = {}
    for key in ("excitation", "attachment", "ionization"):
        tot = np.zeros_like(eps)
        for p in st[key]:
            tot += sigma_on(p, eps)
        out[key] = tot
    thr = sorted(p.threshold for p in st["excitation"] if p.threshold > 0)
    out["first_threshold_eV"] = float(thr[0]) if thr else float("nan")
    return out


def table_floor_te(run_names) -> dict:
    """Te [eV] of the lowest tabulated <eps> per bucket for a table set, and
    the set's overall floor (the engine maps <eps> = 1.5 Te)."""
    lok = LokiRates(run_names)
    per = {b: float(np.exp(lok.curves[b][0][0])) / 1.5 for b in lok.buckets}
    return {"per_bucket_te_eV": per,
            "floor_te_eV": min(per.values()) if per else float("nan"),
            "en_range_td": lok.en_range_td}


# ------------------------------------------------------- cross-database bands
BUCKET_KEYS = {"excitation": "exc", "ionization": "ion", "attachment": "att"}


def database_rates(gas: str, te_grid, root: Path) -> dict:
    """Maxwellian-averaged k(Te) per available evaluation of one gas, for
    the lumped channels (elastic / excitation / ionization / attachment).
    The Maxwellian is a COMMON YARDSTICK here, not a physics claim: the
    same weight applied to every database isolates the sigma difference."""
    out = {}
    sets = {"engine": None}
    sets.update({k: root / v for k, v in ALTERNATES.get(gas, {}).items()}
                if root else {})
    for label, path in sets.items():
        try:
            st = load_set(gas, path) if path is not None else load_set(gas)
        except Exception as e:
            out[label] = {"error": repr(e)[:200]}
            continue
        if path is not None and not Path(path).exists():
            continue
        if label != "engine" and out.get("engine", {}).get("path"):
            a, b = Path(out["engine"]["path"]), Path(st["path"])
            if a.exists() and b.exists() and a.read_bytes() == b.read_bytes():
                continue          # same data under another name: no band
        chans = inelastic_total(st)
        rec = {"database": st["database"], "path": st["path"],
               "elastic_is_effective": st["elastic_is_effective"]}
        for key, short in BUCKET_KEYS.items():
            sig = chans[key]
            if not np.any(sig > 0):
                continue
            rec[short] = [float(maxwell_rate(EPS_GRID, sig, float(t)))
                          for t in te_grid]
        rec["el"] = [float(maxwell_rate(EPS_GRID,
                                        sigma_on(st["elastic"], EPS_GRID),
                                        float(t))) for t in te_grid]
        out[label] = rec
    return out


def spread(rates: dict, te_grid) -> dict:
    """max/min across databases per channel per Te — the band the engine's
    rate for that channel carries from the data alone."""
    labels = [k for k, v in rates.items() if "error" not in v]
    if len(labels) < 2:
        return {}
    out = {}
    for short in ("el", "exc", "ion", "att"):
        cols = [np.array(rates[l][short]) for l in labels if short in rates[l]]
        if len(cols) < 2:
            continue
        stack = np.vstack(cols)
        with np.errstate(divide="ignore", invalid="ignore"):
            band = np.where(stack.min(axis=0) > 0,
                            stack.max(axis=0) / stack.min(axis=0), np.nan)
        out[short] = {"te_eV": list(map(float, te_grid)),
                      "band": [float(x) for x in band],
                      "worst": float(np.nanmax(band)),
                      "worst_te_eV": float(te_grid[int(np.nanargmax(band))]),
                      "n_sets": len(cols)}
    return out


# ------------------------------------------------------------------- figures
def make_figures(atlas: dict, out_dir: Path, te_grid, root: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gases = [g for g in atlas["gases"] if "error" not in atlas["gas"][g]]
    floors = [v["floor_te_eV"] * 1.5 for v in atlas["table_floors"].values()]
    floor_lo, floor_hi = min(floors), max(floors)

    # fig 1 — the elastic cross sections
    fig, ax = plt.subplots(figsize=(9, 6))
    _stagger = iter([(8, -34), (10, 16), (8, -62), (8, 30)])
    for g in gases:
        st = load_set(g)
        ax.loglog(EPS_GRID, sigma_on(st["elastic"], EPS_GRID), lw=1.8,
                  color=COLOR[g],
                  label=f"{g}" + (" (effective)" if st["elastic_is_effective"] else ""))
        r = atlas["gas"][g]["ramsauer"]
        if r and r.get("is_rt_window"):
            ax.plot([r["eps_min_eV"]], [r["sigma_min_m2"]], "v", ms=9,
                    color=COLOR[g], mec="k", mew=0.6, zorder=5)
            ax.annotate(f"{g} R-T min: {r['eps_min_eV']:.2g} eV, "
                        f"x{r['depth']:.0f} deep",
                        (r["eps_min_eV"], r["sigma_min_m2"]),
                        textcoords="offset points",
                        xytext=next(_stagger, (8, -34)), fontsize=8,
                        color=COLOR[g],
                        arrowprops=dict(arrowstyle="-", lw=0.6,
                                        color=COLOR[g], alpha=0.7))
        t = atlas["gas"][g]["channels"]["first_threshold_eV"]
        if np.isfinite(t):
            ax.axvline(t, color=COLOR[g], ls=":", lw=0.8, alpha=0.55)
    ax.axvspan(floor_lo, floor_hi, color="0.85", alpha=0.5, zorder=0)
    ax.annotate("swarm-table floors\n(<eps> at the lowest tabulated E/N)",
                (floor_lo, 4e-24), fontsize=8, color="0.35")
    ax.axvspan(EPS_GRID[0], floor_lo, color="#fdd", alpha=0.35, zorder=0)
    ax.annotate("the engine's beam-mode bulk lives here\n"
                "(left of the grey band the tables are extrapolated)",
                (1.3e-3, 1.4e-24), fontsize=8.5, color="#a33")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(1e-24, 1e-18)
    ax.set_xlabel("electron energy  eps  [eV]")
    ax.set_ylabel(r"elastic momentum-transfer cross section  [m$^2$]")
    ax.set_title("Elastic cross sections: the Ramsauer window sits where the "
                 "beam-mode bulk sits\n(dotted = first excitation threshold "
                 "per gas)", fontsize=10.5)
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(out_dir / "fig1_elastic_sigma.png", dpi=150)
    plt.close(fig)

    # fig 2 — relaxation times at the device density
    n_gas = gas_density_m3(3.0, 300.0)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for g in gases:
        st = load_set(g)
        rel = relaxation_times(st, n_gas)
        axs[0].loglog(EPS_GRID, rel["tau_m_s"] * 1e12, color=COLOR[g], lw=1.7,
                      label=f"{g}")
        axs[1].loglog(EPS_GRID, rel["tau_u_s"] * 1e9, color=COLOR[g], lw=1.7,
                      label=f"{g}  (2m/M = {2 * rel['mass_ratio']:.1e})")
    for a in axs:
        a.axvspan(floor_lo, floor_hi, color="0.85", alpha=0.5, zorder=0)
        a.axvspan(EPS_GRID[0], floor_lo, color="#fdd", alpha=0.3, zorder=0)
        a.set_xlim(1e-2, 1e2)
        a.grid(alpha=0.25, which="both")
        a.set_xlabel("electron energy  eps  [eV]")
    axs[1].axhline(100.0, color="k", ls="--", lw=1.0)
    axs[1].annotate("100 ns beam pulse", (1.2e-2, 115), fontsize=8)
    axs[0].set_ylabel("momentum-transfer time  [ps]")
    axs[1].set_ylabel("elastic ENERGY-relaxation time  [ns]")
    axs[0].set_title("How often a bulk electron collides", fontsize=10)
    axs[1].set_title("How long it takes to give that energy to the gas\n"
                     "(above the dashed line the bulk cannot cool within a "
                     "pulse)", fontsize=10)
    axs[0].legend(fontsize=8.5)
    axs[1].legend(fontsize=8)
    fig.suptitle(f"Relaxation at 3 atm, 300 K (n = {n_gas:.2e} m^-3)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "fig2_relaxation.png", dpi=150)
    plt.close(fig)

    # fig 3 — channels per gas
    n = len(gases)
    fig, axs = plt.subplots(2, (n + 1) // 2, figsize=(4.2 * ((n + 1) // 2), 7.5),
                            squeeze=False)
    for k, g in enumerate(gases):
        a = axs[k // ((n + 1) // 2)][k % ((n + 1) // 2)]
        st = load_set(g)
        ch = inelastic_total(st)
        a.loglog(EPS_GRID, sigma_on(st["elastic"], EPS_GRID), color=COLOR[g],
                 lw=1.8, label="elastic")
        if np.any(ch["excitation"] > 0):
            a.loglog(EPS_GRID, ch["excitation"], "k--", lw=1.2, label="excitation")
        if np.any(ch["ionization"] > 0):
            a.loglog(EPS_GRID, ch["ionization"], "k-.", lw=1.2, label="ionization")
        if np.any(ch["attachment"] > 0):
            a.loglog(EPS_GRID, ch["attachment"], color="#d62728", lw=1.4,
                     label="attachment")
        a.axvspan(EPS_GRID[0], floor_lo, color="#fdd", alpha=0.3, zorder=0)
        a.set_xlim(1e-2, 1e3)
        a.set_ylim(1e-25, 1e-18)
        a.set_title(g, fontsize=11)
        a.grid(alpha=0.25, which="both")
        a.legend(fontsize=7.5, loc="lower left")
        a.set_xlabel("eps [eV]")
        a.set_ylabel(r"sigma [m$^2$]")
    for k in range(len(gases), axs.size):
        axs[k // ((n + 1) // 2)][k % ((n + 1) // 2)].axis("off")
    fig.suptitle("Channels per gas: what is available to an electron at each "
                 "energy (pink = below the swarm tables)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "fig3_channels.png", dpi=150)
    plt.close(fig)

    # fig 4 + 5 — cross-database
    multi = [g for g in gases if len(atlas["databases"].get(g, {})) > 1]
    if multi:
        fig, axs = plt.subplots(1, len(multi), figsize=(5.2 * len(multi), 4.6),
                                squeeze=False)
        for k, g in enumerate(multi):
            a = axs[0][k]
            for label, rec in atlas["databases"][g].items():
                if "error" in rec:
                    continue
                st = load_set(g, rec["path"])
                a.loglog(EPS_GRID, sigma_on(st["elastic"], EPS_GRID), lw=1.6,
                         label=f"{label}: {rec['database'][:26]}")
                att = inelastic_total(st)["attachment"]
                if np.any(att > 0):
                    a.loglog(EPS_GRID, att, ls="--", lw=1.3,
                             color=a.lines[-1].get_color())
            a.axvspan(EPS_GRID[0], floor_lo, color="#fdd", alpha=0.3, zorder=0)
            a.set_xlim(1e-2, 1e3)
            a.set_ylim(1e-24, 1e-18)
            a.set_title(f"{g} — published evaluations\n(solid elastic, dashed "
                        "attachment)", fontsize=10)
            a.set_xlabel("eps [eV]")
            a.set_ylabel(r"sigma [m$^2$]")
            a.grid(alpha=0.25, which="both")
            a.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "fig4_databases.png", dpi=150)
        plt.close(fig)

        fig, axs = plt.subplots(1, len(multi), figsize=(5.2 * len(multi), 4.6),
                                squeeze=False)
        for k, g in enumerate(multi):
            a = axs[0][k]
            sp = atlas["spread"].get(g, {})
            for short, style in (("el", "-"), ("exc", "--"), ("ion", "-."),
                                 ("att", ":")):
                if short not in sp:
                    continue
                a.semilogx(sp[short]["te_eV"], sp[short]["band"], style, lw=1.7,
                           label=f"{short}  worst x{sp[short]['worst']:.2f} "
                                 f"@ {sp[short]['worst_te_eV']:.2g} eV")
            a.axhline(1.0, color="k", lw=0.8)
            a.set_title(f"{g} — database spread in k(Te)", fontsize=10)
            a.set_xlabel("Te [eV]")
            a.set_ylabel("max / min across evaluations")
            a.grid(alpha=0.25, which="both")
            a.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "fig5_rate_band.png", dpi=150)
        plt.close(fig)

    # fig 7 — which channel cools the bulk (built before fig 6 uses ax)
    fig, axs = plt.subplots(1, len(FILLS), figsize=(5.4 * len(FILLS), 4.8),
                            squeeze=False)
    for k, (name, (mix, _rn, p_atm)) in enumerate(FILLS.items()):
        a = axs[0][k]
        lr = mixture_loss_rates(mix, p_atm)
        a.loglog(EPS_GRID, lr["elastic"], "k-", lw=2.0, label="elastic (2m/M)")
        tot = np.zeros_like(EPS_GRID)
        for v in lr["by_channel"].values():
            tot += v
        a.loglog(EPS_GRID, tot, color="#d62728", lw=1.8, label="all inelastic")
        ranked = sorted(lr["by_channel"].items(),
                        key=lambda kv: -np.nanmax(kv[1][(EPS_GRID > 0.05)
                                                        & (EPS_GRID < 3.0)]))
        for nm, v in ranked[:3]:
            a.loglog(EPS_GRID, v, lw=1.1, ls="--", label=nm.split(":")[0] + " "
                     + nm.split("@")[-1])
        fl = 1.5 * atlas["table_floors"][name]["floor_te_eV"]
        a.axvspan(EPS_GRID[0], fl, color="#fdd", alpha=0.3, zorder=0)
        a.set_xlim(2e-2, 3e1)
        a.set_ylim(1e2, 1e12)
        a.set_title(f"{name}\n{p_atm:g} atm", fontsize=9.5)
        a.set_xlabel("eps [eV]")
        a.set_ylabel("energy-loss rate per electron [eV/s]")
        a.grid(alpha=0.25, which="both")
        a.legend(fontsize=7)
    fig.suptitle("Which cross section sets the bulk temperature: in the "
                 "Ramsauer window the elastic term collapses and trace F2 "
                 "vibrational excitation takes over", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out_dir / "fig7_loss_channels.png", dpi=150)
    plt.close(fig)

    # fig 6 — the real fills
    fig, ax = plt.subplots(figsize=(9.5, 6))
    for name, (mix, run_names, p_atm) in FILLS.items():
        n_tot = gas_density_m3(p_atm, 300.0)
        nu_u = np.zeros_like(EPS_GRID)
        for g, x in mix.items():
            try:
                st = load_set(g)
            except Exception:
                continue
            rel = relaxation_times(st, n_tot * x)
            nu_u += rel["nu_u"]
        with np.errstate(divide="ignore"):
            ax.loglog(EPS_GRID, 1e9 / np.where(nu_u > 0, nu_u, np.nan), lw=1.9,
                      label=f"{name}, {p_atm:g} atm")
        fl = atlas["table_floors"][name]["floor_te_eV"] * 1.5
        ax.axvline(fl, ls=":", lw=1.1, color=ax.lines[-1].get_color())
    for label, te_pk, te_end in OPERATING:
        ax.axvspan(1.5 * te_end, 1.5 * te_pk, color="#3182bd", alpha=0.07)
    ax.axhline(100.0, color="k", ls="--", lw=1.0)
    ax.annotate("100 ns beam pulse", (1.1e-2, 115), fontsize=8)
    ax.annotate("measured operating <eps> bands\n(settled -> peak, [0.1.37])",
                (0.07, 1.6e4), fontsize=8.5, color="#26557f")
    ax.set_xlim(1e-2, 1e2)
    ax.set_xlabel(r"mean energy  <eps>  [eV]   (the engine maps <eps> = 1.5 Te)")
    ax.set_ylabel("elastic energy-relaxation time of the mixture  [ns]")
    ax.set_title("The real fills: dotted = that fill's swarm-table floor.\n"
                 "Left of it the engine extrapolates; the Ramsauer window "
                 "makes cooling there slower than the pulse.", fontsize=10.5)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(out_dir / "fig6_mixture.png", dpi=150)
    plt.close(fig)


def mixture_loss_rates(mix: dict, p_atm: float, eps=EPS_GRID) -> dict:
    """Energy-loss RATE [eV/s per electron] by channel for a real fill:
    elastic (2m/M weighted) vs each inelastic channel (its threshold per
    event). This is the answer to 'which cross section sets the bulk
    temperature': in the Ramsauer window the elastic term collapses and a
    trace inelastic channel with a low threshold can dominate."""
    n_tot = gas_density_m3(p_atm, 300.0)
    v = np.sqrt(2.0 * eps * EV / M_E)
    out = {"eps_eV": eps, "elastic": np.zeros_like(eps), "by_channel": {}}
    for g, x in mix.items():
        try:
            st = load_set(g)
        except Exception:
            continue
        n = n_tot * float(x)
        out["elastic"] += (2.0 * float(st["elastic"].mass_ratio) * n
                           * sigma_on(st["elastic"], eps) * v * eps)
        for key in ("excitation", "attachment", "ionization"):
            for pr in st[key]:
                loss = float(pr.threshold) if pr.threshold > 0 else (
                    float(eps.mean()) if key == "attachment" else 0.0)
                if key == "attachment":
                    loss = 0.0        # attachment removes the electron, not its energy
                r = n * sigma_on(pr, eps) * v * loss
                if not np.any(r > 0):
                    continue
                name = f"{g}:{key}:{pr.product or pr.target}@{pr.threshold:g}eV"
                out["by_channel"][name] = out["by_channel"].get(
                    name, np.zeros_like(eps)) + r
    return out


def subfloor_audit(mixture: dict, run_names, eps_probe=(0.15, 0.3, 0.6, 1.2)) -> dict:
    """THE ACTIONABLE PART. Below the swarm table's floor the engine
    extrapolates every bucket by one rule: flat for 0.35 ln-units of <eps>,
    then an exponential taper — i.e. monotonically DOWN. This audit asks,
    per bucket, whether the underlying cross sections say the rate really
    goes down there, using a Maxwellian as the common yardstick (at these
    sub-threshold energies it is also the defensible shape: e-e collisions
    Maxwellianise the bulk at beam densities, and below the first threshold
    there is no inelastic depletion to sustain a swarm tail).

    A bucket whose true rate RISES toward low energy is extrapolated the
    WRONG WAY. That is a known defect class here: the scorecard's OMISSIONS
    record it for F- ('the F- DA curve RISES x13 toward low eps where the
    taper fell'), cured for the KINETIC tier by deepening its measured
    curves. This finds every bucket with the same signature."""
    lok = LokiRates(run_names)
    mx = MaxwellRates(mixture)
    out = {"eps_probe_eV": list(eps_probe), "buckets": {}}
    for b in sorted(set(lok.buckets) & set(mx.buckets)):
        xs, ys = lok.curves[b]
        floor_eps = float(np.exp(xs[0]))
        te_grid = np.geomspace(0.02, max(floor_eps, 4.0) / 1.5 * 3.0, 90)
        k_m = np.array([float(mx.k(b, float(t_))) for t_ in te_grid])
        i_pk = int(np.argmax(k_m))
        te_pk = float(te_grid[i_pk])
        # is the TRUE rate still rising as energy falls to the table floor?
        k_at_floor = float(mx.k(b, floor_eps / 1.5))
        rises = bool(te_pk < floor_eps / 1.5 and k_m[i_pk] > k_at_floor * 1.05)
        rec = {"floor_eps_eV": floor_eps, "floor_te_eV": floor_eps / 1.5,
               "maxwell_peak_te_eV": te_pk,
               "maxwell_peak_over_floor": float(k_m[i_pk] / k_at_floor)
               if k_at_floor > 0 else float("inf"),
               "true_rate_rises_below_floor": rises, "probe": {}}
        for e in eps_probe:
            te = e / 1.5
            kl = float(lok.k(b, te))
            km = float(mx.k(b, te))
            rec["probe"][f"{e:g} eV"] = {
                "k_loki_tapered": kl, "k_maxwell": km,
                "ratio": (km / kl) if kl > 0 else float("inf")}
        out["buckets"][b] = rec
    return out


# ---------------------------------------------------------------------- main
def build_atlas(gases, te_grid) -> dict:
    root = _loki_input_root()
    atlas = {"gases": list(gases), "gas": {}, "databases": {}, "spread": {},
             "table_floors": {}, "loki_input_root": str(root) if root else None,
             "operating_points": OPERATING}
    for g in gases:
        try:
            st = load_set(g)
        except Exception as e:
            atlas["gas"][g] = {"error": repr(e)[:200]}
            continue
        n_gas = gas_density_m3(3.0, 300.0)
        rel = relaxation_times(st, n_gas)
        ch = inelastic_total(st)
        r = ramsauer(st["elastic"])
        i_ops = [int(np.argmin(np.abs(EPS_GRID - e))) for e in (0.1, 0.3, 1.0, 3.0)]
        atlas["gas"][g] = {
            "database": st["database"], "path": st["path"],
            "elastic_is_effective": st["elastic_is_effective"],
            "mass_ratio": rel["mass_ratio"],
            "ramsauer": r,
            "channels": {"first_threshold_eV": ch["first_threshold_eV"],
                         "has_attachment": bool(np.any(ch["attachment"] > 0))},
            "tau_u_ns_at": {f"{EPS_GRID[i]:.3g} eV": float(rel["tau_u_s"][i] * 1e9)
                            for i in i_ops},
            "tau_m_ps_at": {f"{EPS_GRID[i]:.3g} eV": float(rel["tau_m_s"][i] * 1e12)
                            for i in i_ops},
        }
        atlas["databases"][g] = database_rates(g, te_grid, root)
        sp = spread(atlas["databases"][g], te_grid)
        if sp:
            atlas["spread"][g] = sp
    for name, (mix, run_names, p_atm) in FILLS.items():
        atlas["table_floors"][name] = table_floor_te(run_names)
        lr = mixture_loss_rates(mix, p_atm)
        rec = {}
        for e_probe in (0.1, 0.3, 1.0, 3.0):
            i = int(np.argmin(np.abs(EPS_GRID - e_probe)))
            tot_in = {k: float(v[i]) for k, v in lr["by_channel"].items()
                      if v[i] > 0}
            top = sorted(tot_in.items(), key=lambda kv: -kv[1])[:3]
            el = float(lr["elastic"][i])
            s_in = sum(tot_in.values())
            rec[f"{EPS_GRID[i]:.3g} eV"] = {
                "elastic_eV_per_s": el,
                "inelastic_eV_per_s": s_in,
                "inelastic_over_elastic": (s_in / el) if el > 0 else float("inf"),
                "top_channels": [{"channel": k, "eV_per_s": v} for k, v in top]}
        atlas.setdefault("loss_channels", {})[name] = rec
        atlas.setdefault("subfloor", {})[name] = subfloor_audit(mix, run_names)
    return atlas


def format_atlas(atlas: dict) -> str:
    out = []
    out.append("== elastic cross sections: the Ramsauer window per gas "
               "(3 atm, 300 K for the times)")
    out.append(f"  {'gas':4s} {'database':26s} {'2m/M':>9s} {'R-T min':>10s} "
               f"{'depth':>7s} {'exc thr':>8s} {'tau_u(0.3eV)':>13s} "
               f"{'tau_u(3eV)':>11s}")
    for g in atlas["gases"]:
        d = atlas["gas"].get(g, {})
        if "error" in d or not d:
            out.append(f"  {g:4s} -- {d.get('error', 'missing')}")
            continue
        r = d["ramsauer"]
        tu = d["tau_u_ns_at"]
        k03 = next((v for k, v in tu.items() if k.startswith("0.3")), float("nan"))
        k3 = next((v for k, v in tu.items() if k.startswith("3")), float("nan"))
        out.append(f"  {g:4s} {d['database'][:26]:26s} {2 * d['mass_ratio']:9.2e} "
                   f"{r.get('eps_min_eV', float('nan')):8.3g}eV "
                   f"{r.get('depth', float('nan')):7.1f} "
                   f"{d['channels']['first_threshold_eV']:8.3g} "
                   f"{k03:11.4g}ns {k3:9.4g}ns")
    out.append("")
    out.append("== swarm-table floors (below these the engine extrapolates by "
               "the flat+taper rule)")
    for name, d in atlas["table_floors"].items():
        out.append(f"  {name:24s} E/N {d['en_range_td'][0]:g}-{d['en_range_td'][1]:g} Td, "
                   f"floor <eps> {1.5 * d['floor_te_eV']:.3g} eV "
                   f"(Te {d['floor_te_eV']:.3g} eV)")
    if atlas.get("loss_channels"):
        out.append("")
        out.append("== what actually cools the bulk, per fill "
                   "(energy-loss rate per electron; elastic vs inelastic)")
        for name, rec in atlas["loss_channels"].items():
            out.append(f"  {name}")
            for e, d in rec.items():
                top = d["top_channels"][0]["channel"] if d["top_channels"] else "-"
                out.append(f"      <eps> {e:>8s}: inelastic/elastic "
                           f"x{d['inelastic_over_elastic']:.3g}   top: {top}")
    if atlas.get("subfloor"):
        out.append("")
        out.append("== sub-floor audit: buckets the flat+taper rule "
                   "extrapolates the WRONG WAY")
        out.append("   (the engine's rule always falls below the floor; a "
                   "'rises' bucket physically climbs toward low energy)")
        for name, aud in atlas["subfloor"].items():
            wrong = {b: r for b, r in aud["buckets"].items()
                     if r["true_rate_rises_below_floor"]}
            out.append(f"  {name}: {len(wrong)} of {len(aud['buckets'])} "
                       f"buckets rise below the floor")
            for b, r in sorted(wrong.items(),
                               key=lambda kv: -kv[1]["probe"]["0.3 eV"]["ratio"]):
                pr = r["probe"]["0.3 eV"]
                out.append(f"      {b:8s} floor <eps> {r['floor_eps_eV']:.3g} eV, "
                           f"true peak at Te {r['maxwell_peak_te_eV']:.3g} eV "
                           f"(x{r['maxwell_peak_over_floor']:.2f} the floor value); "
                           f"at <eps> 0.3 eV engine is LOW by x{pr['ratio']:.4g}")
    if atlas["spread"]:
        out.append("")
        out.append("== published-evaluation spread in Maxwellian-averaged k(Te) "
                   "(max/min; a data uncertainty, not a knob)")
        for g, sp in atlas["spread"].items():
            labels = ", ".join(
                f"{k}={v.get('database', '?')[:22]}"
                for k, v in atlas["databases"][g].items() if "error" not in v)
            out.append(f"  {g}: {labels}")
            for short, rec in sp.items():
                out.append(f"      {short:4s} worst x{rec['worst']:.3g} at Te "
                           f"{rec['worst_te_eV']:.3g} eV over {rec['n_sets']} sets")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gases", nargs="+", default=list(GASES))
    ap.add_argument("--no-figs", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "rebuild_out" / "sigma_atlas"))
    args = ap.parse_args(argv)

    te_grid = np.geomspace(0.05, 20.0, 60)
    atlas = build_atlas(args.gases, te_grid)
    print(format_atlas(atlas))

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(atlas, indent=1, default=str))
    if not args.no_figs:
        make_figures(atlas, out_dir, te_grid, _loki_input_root())
        print(f"\nfigures + report -> {out_dir}")
    else:
        print(f"\nreport -> {out_dir / 'report.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
