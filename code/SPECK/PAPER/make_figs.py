"""Generate the two figures for the in-silico microwave-reheat draft from the
RECORDED run JSONs (no numbers typed by hand). Run from Kinetic_Modeling root:

    python3 SPECK/PAPER/make_figs.py

Writes fig_during_pulse.pdf and fig_robust_and_wallplug.pdf next to itself.
Data of record:
  SPECK/rebuild_out/mw_benchmarks/sign_decidability.json   ([0.1.28])
  SPECK/rebuild_out/copump_scan/copump_verify.json         ([0.1.23])
  SPECK/rebuild_out/mw_benchmarks/ib_feasibility.json      ([0.1.26])
  SPECKLE/runs/cp27k{,_dt25,_dt12}_copump_report.json      (SPECKLE 08-26/27)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent          # Kinetic_Modeling root (SPECK/PAPER/..)
SPECK = ROOT / "SPECK"

INK = "#23202C"
MUTED = "#6B6E7C"
ACCENT = "#5B2D8F"      # deep violet (the 248 nm identity)
BAND = "#5B2D8F"
GRID = "#DDDAD4"

plt.rcParams.update({
    "font.size": 8.0, "axes.labelsize": 8.0, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
})


def _panel(ax):
    ax.grid(True, color=GRID, linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ data
sd = json.loads((SPECK / "rebuild_out/mw_benchmarks/sign_decidability.json").read_text())
cv = json.loads((SPECK / "rebuild_out/copump_scan/copump_verify.json").read_text())
ib = json.loads((SPECK / "rebuild_out/mw_benchmarks/ib_feasibility.json").read_text())

A = sd["A_kstep_map"]
B = sd["B_neopt_afterglow_brackets"]
f_lo, f_hi = A["f_lo"], A["f_hi"]
xs = [r["factor"] for r in A["rows"]]
ms = [100 * r["marginal"] for r in A["rows"]]
kin_nom = 100 * cv["tier"]["kinetic"]["marginal"]

dts, sp_marg = [], []
for name in ("cp27k", "cp27k_dt25", "cp27k_dt12"):
    r = json.loads((ROOT / "SPECKLE" / "runs" / f"{name}_copump_report.json").read_text())
    dts.append(1e12 * r["config"]["dt"])
    sp_marg.append(100 * r["marginal"])

# ------------------------------------------- Fig 1: the during-pulse channel
fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.6, 4.3))
fig.subplots_adjust(left=0.16, right=0.97, top=0.94, bottom=0.10, hspace=0.42)

_panel(a1)
a1.axvspan(f_lo, f_hi, color=BAND, alpha=0.10, lw=0)
a1.axhline(0, color=MUTED, lw=0.8)
a1.plot(xs, ms, "-o", color=ACCENT, lw=1.4, ms=3.4, zorder=3)
a1.plot([1.0], [kin_nom], "D", mfc="none", mec=MUTED, ms=5, mew=1.1, zorder=3)
a1.annotate("kinetic tier,\nnominal rate", (1.0, kin_nom), xytext=(0.36, -8.0),
            fontsize=7, color=MUTED, ha="left",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
a1.annotate(f"sign-positive window\n$\\times${f_lo:.2f} – $\\times${f_hi:.2f}",
            (1.0, 16.4), fontsize=7, color=ACCENT, ha="center")
a1.set_xscale("log")
a1.minorticks_off()
a1.set_xticks([1 / 3, 0.5, 1, 2, 3])
a1.set_xticklabels([r"$\times\frac{1}{3}$", r"$\times$0.5", r"$\times$1",
                    r"$\times$2", r"$\times$3"])
a1.set_xlim(0.3, 3.3)
a1.set_ylim(-9, 19)
a1.set_xlabel("KSTEP multiplier  (e + Kr$^*$ $\\to$ Kr$^+$ + 2e)")
a1.set_ylabel("co-pump marginal  $\\Delta E_L/E_{\\rm MW}$  (%)")
a1.set_title("(a) the during-pulse sign is a one-rate question", loc="left")

_panel(a2)
a2.axhline(0, color=MUTED, lw=0.8)
a2.plot(dts, sp_marg, "-s", color=ACCENT, lw=1.4, ms=3.8, zorder=3)
for x, y in zip(dts, sp_marg):
    a2.annotate(f"{y:+.1f}%", (x, y), xytext=(0, 6 if y > -20 else -11),
                textcoords="offset points", fontsize=7, ha="center", color=INK)
a2.set_xscale("log")
a2.minorticks_off()
a2.set_xticks(dts)
a2.set_xticklabels([f"{d:g}" for d in dts])
a2.set_xlim(9, 70)
a2.set_ylim(-32, 26)
a2.invert_xaxis()
a2.set_xlabel("timestep  $\\Delta t$  (ps)   [independent 1-D discharge model]")
a2.set_ylabel("co-pump marginal  (%)")
a2.set_title("(b) the coarse-$\\Delta t$ positive dies under refinement",
             loc="left")

fig.savefig(HERE / "fig_during_pulse.pdf")
plt.close(fig)

# ----------------------------- Fig 2: robustness strip + wall-plug inequality
fig, (b1, b2) = plt.subplots(2, 1, figsize=(3.6, 4.3))
fig.subplots_adjust(left=0.16, right=0.97, top=0.94, bottom=0.13, hspace=0.52)

_panel(b1)
bm = sorted(100 * r["marginal"] for r in B["rows"])
b1.axhline(0, color=MUTED, lw=0.8)
b1.plot(range(1, len(bm) + 1), bm, "o", color=ACCENT, ms=3.0, zorder=3)
b1.annotate(f"all {len(bm)} rate variants positive\n"
            f"({bm[0]:+.1f}% … {bm[-1]:+.1f}%)",
            (len(bm) / 2 + 0.5, 3.6), fontsize=7, color=ACCENT, ha="center")
b1.set_xlim(0, len(bm) + 1)
b1.set_ylim(-2, 15)
b1.set_xlabel("bracket variant, ranked")
b1.set_ylabel("afterglow marginal  (%)")
b1.set_title("(a) the afterglow gain survives every bracket", loc="left")

_panel(b2)
rows = [r for r in ib["rows"]
        if r["label"].startswith(("1 GHz", "2.45 GHz", "9.5 GHz", "95 GHz"))
        or "CO2" in r["label"]]
labels = ["1 GHz\nmagnetron", "2.45 GHz\nmagnetron", "9.5 GHz\nX-band",
          "95 GHz\ngyrotron", "CO$_2$ laser\n28 THz"]
vals = [100 * r["eta_times_m_best"] for r in rows]
xpos = range(len(rows))
b2.bar(xpos, vals, width=0.62, color=ACCENT, alpha=0.85, zorder=3)
bench = 100 * ib["eta_beam_wallplug"]
b2.axhline(bench, color=MUTED, lw=1.0, ls=(0, (4, 2)))
b2.annotate(f"e-beam benchmark {bench:.0f}%", (4.45, bench + 0.3),
            fontsize=7, color=MUTED, ha="right")
for x, v in zip(xpos, vals):
    b2.annotate(f"{v:.1f}", (x, v + 0.18), fontsize=7, ha="center", color=INK)
b2.set_xticks(list(xpos))
b2.set_xticklabels(labels, fontsize=6.8)
b2.set_ylim(0, 10.4)
b2.set_ylabel("$\\eta_{\\rm source}\\times m_{\\rm best}$  (%)")
b2.set_title("(b) the wall-plug inequality is barely reachable", loc="left")

fig.savefig(HERE / "fig_robust_and_wallplug.pdf")
plt.close(fig)

print("wrote", HERE / "fig_during_pulse.pdf")
print("wrote", HERE / "fig_robust_and_wallplug.pdf")

# --------------------- Fig 3: the 1-D verification (SPECKLE ladders) ------
RUNS = ROOT / "SPECKLE" / "runs"


def _rep(name):
    p = RUNS / f"{name}_afterglow_report.json"
    return json.loads(p.read_text()) if p.exists() else None


margin_pts = []
for name in ("mg_r03", "mg_r06", "mg_r08", "ag250k"):
    r = _rep(name)
    if r:
        margin_pts.append((r["config"]["r_out"], 100 * r["marginal"],
                           1e3 * r["baseline"]["E_laser_J"]))

timing_pts = []
for name in ("tw_t170", "tw_t250", "tw_t330"):
    r = _rep(name)
    if r:
        timing_pts.append((r["config"]["mw_t0"], 100 * r["marginal"]))
refine = []                      # (t0, marginal %, label) open markers
for name, lab in (("tw_t170_dt25", "$\\Delta t/2$"),
                  ("tw_t250_dt25", "$\\Delta t/2$"),
                  ("tw_t250_dt12", "$\\Delta t/4$")):
    r = _rep(name)
    if r:
        refine.append((r["config"]["mw_t0"], 100 * r["marginal"], lab))

prof_pts = {}
for name, lab in (("pf_edge", "edge"), ("pf_center", "center"),
                  ("ag250k", "uniform")):
    r = _rep(name)
    if r:
        prof_pts[lab] = 100 * r["marginal"]

if margin_pts and timing_pts:
    fig, (c1, c2) = plt.subplots(2, 1, figsize=(3.6, 4.3))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.94, bottom=0.11,
                        hspace=0.46)

    _panel(c1)
    c1.axhline(0, color=MUTED, lw=0.8)
    xs_m, ys_m, e_off = zip(*sorted(margin_pts))
    c1.plot(xs_m, ys_m, "-o", color=ACCENT, lw=1.4, ms=3.6, zorder=3)
    for x, y, e in zip(xs_m, ys_m, e_off):
        c1.annotate(f"{e:.1f} mJ", (x, y), xytext=(0, -11),
                    textcoords="offset points", fontsize=6.3, ha="center",
                    color=MUTED)
    off = {"edge": -0.028, "uniform": 0.0, "center": 0.028}
    mk = {"edge": "s", "uniform": "o", "center": "^"}
    for lab, y in prof_pts.items():
        if lab == "uniform":
            continue
        c1.plot([0.9 + off[lab]], [y], mk[lab], mfc="none", mec=MUTED,
                ms=4.6, mew=1.0, zorder=3)
        c1.annotate(lab, (0.9 + off[lab], y), xytext=(-2, -11),
                    textcoords="offset points", fontsize=6.3, color=MUTED,
                    ha="center")
    c1.set_xlim(0.24, 0.99)
    c1.set_ylim(-0.15, 1.65)
    c1.set_xlabel("output-coupler reflectivity $r_{\\rm out}$ "
                  "(optical margin knob)")
    c1.set_ylabel("afterglow marginal  (%)")
    c1.set_title("(a) the gain follows the optical margin", loc="left")

    _panel(c2)
    c2.axhline(0, color=MUTED, lw=0.8)
    xs_t, ys_t = zip(*sorted(timing_pts))
    c2.plot(xs_t, ys_t, "-o", color=ACCENT, lw=1.4, ms=3.6, zorder=3)
    for x, y in zip(xs_t, ys_t):
        c2.annotate(f"{y:+.2f}%", (x, y), xytext=(0, 6),
                    textcoords="offset points", fontsize=6.6, ha="center",
                    color=INK)
    for x, y, lab in refine:
        c2.plot([x], [y], "D", mfc="none", mec=MUTED, ms=5, mew=1.1,
                zorder=3)
        c2.annotate(lab, (x, y), xytext=(6, -9),
                    textcoords="offset points", fontsize=6.6, color=MUTED)
    c2.set_xlim(150, 350)
    c2.set_ylim(-7.5, 13.5)
    c2.set_xlabel("afterglow window start  $t_0$  (ns)")
    c2.set_ylabel("afterglow marginal  (%)")
    c2.set_title("(b) the gain dies past the cavity drain", loc="left")

    fig.savefig(HERE / "fig_1d_verification.pdf")
    plt.close(fig)
    print("wrote", HERE / "fig_1d_verification.pdf")
else:
    print("fig_1d_verification: ladder reports not found, skipped")
