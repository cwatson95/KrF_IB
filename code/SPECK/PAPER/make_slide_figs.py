"""Single-panel, slide-scale versions of the paper figures, for the
Beamer deck. Same recorded JSON sources as make_figs.py — no hand-typed
numbers. Run from the Kinetic_Modeling root:

    python3 SPECK/PAPER/make_slide_figs.py

Writes s_margin/s_timing/s_bracket/s_wallplug/s_kstep/s_dtseries.pdf
next to itself.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SPECK = ROOT / "SPECK"
RUNS = ROOT / "SPECKLE" / "runs"

INK = "#23202C"
MUTED = "#6B6E7C"
ACCENT = "#5B2D8F"
GRID = "#DDDAD4"

plt.rcParams.update({
    "font.size": 11.5, "axes.labelsize": 11.5, "axes.titlesize": 12,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})
SIZE = (5.4, 3.15)


def panel():
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.17)
    ax.grid(True, color=GRID, linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return fig, ax


def save(fig, name):
    fig.savefig(HERE / name)
    plt.close(fig)
    print("wrote", HERE / name)


def rep(name):
    p = RUNS / f"{name}_afterglow_report.json"
    return json.loads(p.read_text()) if p.exists() else None


sd = json.loads((SPECK / "rebuild_out/mw_benchmarks/sign_decidability.json").read_text())
cv = json.loads((SPECK / "rebuild_out/copump_scan/copump_verify.json").read_text())
ib = json.loads((SPECK / "rebuild_out/mw_benchmarks/ib_feasibility.json").read_text())

# ---------------------------------------------------------- s_margin
pts = []
for name in ("mg_r03", "mg_r06", "mg_r08", "ag250k"):
    r = rep(name)
    if r:
        pts.append((r["config"]["r_out"], 100 * r["marginal"],
                    1e3 * r["baseline"]["E_laser_J"]))
prof = {}
for name, lab in (("pf_edge", "edge"), ("pf_center", "center")):
    r = rep(name)
    if r:
        prof[lab] = 100 * r["marginal"]
fig, ax = panel()
ax.axhline(0, color=MUTED, lw=0.9)
xs, ys, eo = zip(*sorted(pts))
ax.plot(xs, ys, "-o", color=ACCENT, lw=1.8, ms=5, zorder=3)
for x, y, e in zip(xs, ys, eo):
    ax.annotate(f"{e:.0f} mJ", (x, y), xytext=(0, -15),
                textcoords="offset points", fontsize=9, ha="center",
                color=MUTED)
off = {"edge": -0.03, "center": 0.03}
mk = {"edge": "s", "center": "^"}
for lab, y in prof.items():
    ax.plot([0.9 + off[lab]], [y], mk[lab], mfc="none", mec=MUTED,
            ms=6.5, mew=1.2, zorder=3)
    ax.annotate(lab, (0.9 + off[lab], y), xytext=(-4, -14 if lab == "edge" else 8),
                textcoords="offset points", fontsize=9, color=MUTED,
                ha="center")
ax.set_xlim(0.24, 0.99)
ax.set_ylim(-0.2, 1.7)
ax.set_xlabel("output-coupler reflectivity $r_{\\rm out}$")
ax.set_ylabel("afterglow marginal (%)")
save(fig, "s_margin.pdf")

# ---------------------------------------------------------- s_timing
tp = []
for name in ("tw_t170", "tw_t250", "tw_t330"):
    r = rep(name)
    if r:
        tp.append((r["config"]["mw_t0"], 100 * r["marginal"]))
refine = []
for name, lab in (("tw_t170_dt25", "$\\Delta t/2$"),
                  ("tw_t250_dt25", "$\\Delta t/2$"),
                  ("tw_t250_dt12", "$\\Delta t/4$")):
    r = rep(name)
    if r:
        refine.append((r["config"]["mw_t0"], 100 * r["marginal"], lab))
fig, ax = panel()
ax.axhline(0, color=MUTED, lw=0.9)
xs, ys = zip(*sorted(tp))
ax.plot(xs, ys, "-o", color=ACCENT, lw=1.8, ms=5, zorder=3)
for x, y in zip(xs, ys):
    ax.annotate(f"{y:+.1f}%", (x, y), xytext=(0, 8),
                textcoords="offset points", fontsize=10, ha="center",
                color=INK)
for x, y, lab in refine:
    ax.plot([x], [y], "D", mfc="none", mec=MUTED, ms=6.5, mew=1.2,
            zorder=3)
    ax.annotate(lab, (x, y), xytext=(8, -10), textcoords="offset points",
                fontsize=9, color=MUTED)
ax.set_xlim(150, 350)
ax.set_ylim(-8, 14)
ax.set_xlabel("afterglow window start $t_0$ (ns)")
ax.set_ylabel("afterglow marginal (%)")
save(fig, "s_timing.pdf")

# --------------------------------------------------------- s_bracket
bm = sorted(100 * r["marginal"]
            for r in sd["B_neopt_afterglow_brackets"]["rows"])
fig, ax = panel()
ax.axhline(0, color=MUTED, lw=0.9)
ax.plot(range(1, len(bm) + 1), bm, "o", color=ACCENT, ms=4.2, zorder=3)
ax.annotate(f"all {len(bm)} rate variants positive "
            f"({bm[0]:+.1f}% … {bm[-1]:+.1f}%)",
            (len(bm) / 2 + 0.5, 3.6), fontsize=10.5, color=ACCENT,
            ha="center")
ax.set_xlim(0, len(bm) + 1)
ax.set_ylim(-2, 15)
ax.set_xlabel("bracket variant (EST $\\times\\frac{1}{3}..\\times$3 + "
              "literature alternates), ranked")
ax.set_ylabel("afterglow marginal (%)")
save(fig, "s_bracket.pdf")

# -------------------------------------------------------- s_wallplug
rows = [r for r in ib["rows"]
        if r["label"].startswith(("1 GHz", "2.45 GHz", "9.5 GHz", "95 GHz"))
        or "CO2" in r["label"]]
labels = ["1 GHz\nmagnetron", "2.45 GHz\nmagnetron", "9.5 GHz\nX-band",
          "95 GHz\ngyrotron", "CO$_2$ laser\n28 THz"]
vals = [100 * r["eta_times_m_best"] for r in rows]
bench = 100 * ib["eta_beam_wallplug"]
fig, ax = panel()
ax.bar(range(len(rows)), vals, width=0.62, color=ACCENT, alpha=0.88,
       zorder=3)
ax.axhline(bench, color=MUTED, lw=1.2, ls=(0, (4, 2)))
ax.annotate(f"e-beam benchmark {bench:.0f}%", (4.45, bench + 0.32),
            fontsize=10, color=MUTED, ha="right")
for x, v in zip(range(len(rows)), vals):
    ax.annotate(f"{v:.1f}", (x, v + 0.2), fontsize=10, ha="center",
                color=INK)
ax.set_xticks(range(len(rows)))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylim(0, 10.6)
ax.set_ylabel("$\\eta_{\\rm source}\\times m_{\\rm best}$ (%)")
save(fig, "s_wallplug.pdf")

# ----------------------------------------------------------- s_kstep
A = sd["A_kstep_map"]
f_lo, f_hi = A["f_lo"], A["f_hi"]
xs = [r["factor"] for r in A["rows"]]
ms = [100 * r["marginal"] for r in A["rows"]]
kin = 100 * cv["tier"]["kinetic"]["marginal"]
fig, ax = panel()
ax.axvspan(f_lo, f_hi, color=ACCENT, alpha=0.10, lw=0)
ax.axhline(0, color=MUTED, lw=0.9)
ax.plot(xs, ms, "-o", color=ACCENT, lw=1.8, ms=5, zorder=3)
ax.plot([1.0], [kin], "D", mfc="none", mec=MUTED, ms=7, mew=1.3,
        zorder=3)
ax.annotate("kinetic tier,\nnominal rate", (1.0, kin), xytext=(0.36, -7.5),
            fontsize=9.5, color=MUTED, ha="left",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
ax.annotate(f"sign-positive window  $\\times${f_lo:.2f} – $\\times${f_hi:.2f}",
            (1.0, 16.6), fontsize=10.5, color=ACCENT, ha="center")
ax.set_xscale("log")
ax.minorticks_off()
ax.set_xticks([1 / 3, 0.5, 1, 2, 3])
ax.set_xticklabels([r"$\times\frac{1}{3}$", r"$\times$0.5", r"$\times$1",
                    r"$\times$2", r"$\times$3"])
ax.set_xlim(0.3, 3.3)
ax.set_ylim(-9, 20)
ax.set_xlabel("KSTEP multiplier  (e + Kr$^*$ $\\to$ Kr$^+$ + 2e)")
ax.set_ylabel("co-pump marginal (%)")
save(fig, "s_kstep.pdf")

# --------------------------------------------------------- s_dtseries
dts, sm = [], []
for name in ("cp27k", "cp27k_dt25", "cp27k_dt12"):
    r = json.loads((RUNS / f"{name}_copump_report.json").read_text())
    dts.append(1e12 * r["config"]["dt"])
    sm.append(100 * r["marginal"])
fig, ax = panel()
ax.axhline(0, color=MUTED, lw=0.9)
ax.plot(dts, sm, "-s", color=ACCENT, lw=1.8, ms=5.5, zorder=3)
for x, y in zip(dts, sm):
    ax.annotate(f"{y:+.1f}%", (x, y), xytext=(0, 9 if y > -20 else -15),
                textcoords="offset points", fontsize=10, ha="center",
                color=INK)
ax.set_xscale("log")
ax.minorticks_off()
ax.set_xticks(dts)
ax.set_xticklabels([f"{d:g}" for d in dts])
ax.set_xlim(9, 70)
ax.set_ylim(-32, 26)
ax.invert_xaxis()
ax.set_xlabel("timestep $\\Delta t$ (ps)  [independent 1-D model]")
ax.set_ylabel("co-pump marginal (%)")
save(fig, "s_dtseries.pdf")
