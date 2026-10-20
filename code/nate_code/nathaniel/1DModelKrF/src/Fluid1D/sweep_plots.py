"""Sweep-overlay figures: one figure per diagnostic, one color-coded curve per F₂ value.

Consumes the per-point ``diagnostics.h5`` archives a :func:`Fluid1D.sweep.sweep_f2_parallel`
run writes to ``runs/<sweep>/f2_<value>/`` (via the :class:`~Fluid1D.diagnostics
.PlasmaScalarTraces` / ``LaserOutput`` diagnostics), plus mixture-level Boltzmann tables for
the α/η figure. Experiment files declare every output filename, x/y dataset, scale, and
optional fixed range. Unpinned axes fit the finite visible data without decade rounding.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from diagnostics import DataRef, data_limits
from particles.constants import ATM_PA, K_B

from .convergence import pulse_fwhm


@dataclass(frozen=True)
class OverlayFigure:
    """One run-file-owned sweep overlay, one curve per F₂ value."""

    filename: str | Path
    x: DataRef
    y: DataRef
    x_scale: str = "linear"
    y_scale: str = "linear"
    x_limits: tuple[float, float] | None = None
    y_limits: tuple[float, float] | None = None
    x_include_zero: bool | None = None
    y_include_zero: bool | None = None
    title: str | None = None


@dataclass(frozen=True)
class OverlayPanel:
    """One axes of a multi-panel sweep overlay (a panel of Chen Fig. 15)."""

    x: DataRef
    y: DataRef
    tag: str | None = None  # corner annotation, e.g. "(a)"
    x_scale: str = "linear"
    y_scale: str = "linear"
    x_limits: tuple[float, float] | None = None
    y_limits: tuple[float, float] | None = None
    x_include_zero: bool | None = None
    y_include_zero: bool | None = None


@dataclass(frozen=True)
class PanelOverlayFigure:
    """Run-file-owned grid of sweep-overlay panels sharing one F₂ legend (Chen Fig. 15 style)."""

    filename: str | Path
    panels: tuple[OverlayPanel, ...]
    ncols: int = 2
    legend_title: str = r"F$_2$ [%]"


def fit_y_range(arrays, *, pad: float = 1.2) -> tuple[float, float]:
    """Return finite-data-aware limits for in-memory sweep arrays.

    ``pad`` uses a multiplicative API (``1.2`` means 20% total padding), while
    signed data now fit their actual asymmetric range instead of wasting space symmetrically.
    """
    return data_limits(
        *(np.asarray(a, dtype=np.float64) for a in arrays),
        pad=max(float(pad) - 1.0, 0.0),
    )


def sweep_point_dirs(sweep_dir: Path | str) -> list[tuple[float, Path]]:
    """(f2_fraction, run_dir) for every ``f2_*`` point under ``sweep_dir``, sorted by F₂."""
    out = []
    for p in Path(sweep_dir).glob("f2_*"):
        if (p / "diagnostics.h5").exists():
            out.append((float(p.name[3:]), p))
    return sorted(out)


def _f2_colors(n: int):
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    return [cmap(0.1 + 0.8 * i / max(n - 1, 1)) for i in range(n)]


def _read_point_curves(points, colors, x_ref: DataRef, y_ref: DataRef):
    """(f2, x, y, color) per sweep point holding both datasets; shape-checked."""
    import h5py

    plotted: list[tuple[float, NDArray[np.float64], NDArray[np.float64], object]] = []
    for (f2, run_dir), color in zip(points, colors):
        with h5py.File(run_dir / "diagnostics.h5", "r") as h5:
            if x_ref.dataset not in h5 or y_ref.dataset not in h5:
                continue
            x = x_ref.read(h5)
            y = y_ref.read(h5)
        if x.shape != y.shape or x.ndim != 1:
            raise ValueError(
                f"sweep overlay x/y shapes differ for {run_dir}: {x.shape} != {y.shape}"
            )
        plotted.append((f2, x, y, color))
    return plotted


def _draw_overlay_axis(ax, spec, plotted, *, labelled: bool) -> None:
    """Draw the per-F₂ curves + fitted/pinned limits of one overlay axes."""
    fitted_x = spec.x_limits or data_limits(
        *(x for _, x, _, _ in plotted),
        scale=spec.x_scale, pad=0.02, include_zero=spec.x_include_zero,
    )
    visible_y: list[NDArray[np.float64]] = []
    for f2, x, y, color in plotted:
        mask = (x >= fitted_x[0]) & (x <= fitted_x[1])
        visible_y.append(y[mask])
        label = f"{100 * f2:.3g}% F$_2$" if labelled else None
        ax.plot(x, y, color=color, lw=1.4, label=label)
    ax.set_xlabel(spec.x.label)
    ax.set_ylabel(spec.y.label)
    ax.set_xscale(spec.x_scale)
    ax.set_yscale(spec.y_scale)
    ax.set_xlim(*fitted_x)
    ax.set_ylim(*(spec.y_limits or data_limits(
        *visible_y, scale=spec.y_scale, pad=0.05,
        include_zero=spec.y_include_zero,
    )))


def plot_sweep_overlays(
    sweep_dir: Path | str,
    figures: Sequence[OverlayFigure],
    *,
    out_dir: Path | str | None = None,
) -> list[Path]:
    """Render the run-file-declared overlays from per-point HDF5 archives."""
    import matplotlib.pyplot as plt

    points = sweep_point_dirs(sweep_dir)
    if not points:
        raise FileNotFoundError(f"no f2_*/diagnostics.h5 under {sweep_dir}")
    out = Path(out_dir) if out_dir is not None else Path(sweep_dir)
    out.mkdir(parents=True, exist_ok=True)
    colors = _f2_colors(len(points))
    paths: list[Path] = []
    for spec in figures:
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        plotted = _read_point_curves(points, colors, spec.x, spec.y)
        if not plotted:
            plt.close(fig)
            raise KeyError(
                f"no sweep point contains {spec.x.dataset!r} and {spec.y.dataset!r}"
            )
        _draw_overlay_axis(ax, spec, plotted, labelled=True)
        if spec.title:
            ax.set_title(spec.title)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = out / spec.filename
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_sweep_panels(
    sweep_dir: Path | str,
    figure: PanelOverlayFigure,
    *,
    out_dir: Path | str | None = None,
) -> Path:
    """Render one multi-panel overlay figure (Chen Fig. 15 style, shared F₂ legend)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    points = sweep_point_dirs(sweep_dir)
    if not points:
        raise FileNotFoundError(f"no f2_*/diagnostics.h5 under {sweep_dir}")
    out = Path(out_dir) if out_dir is not None else Path(sweep_dir)
    out.mkdir(parents=True, exist_ok=True)
    colors = _f2_colors(len(points))
    n = len(figure.panels)
    ncols = max(1, int(figure.ncols))
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.2 * nrows),
                             squeeze=False)
    for idx, panel in enumerate(figure.panels):
        ax = axes[idx // ncols][idx % ncols]
        plotted = _read_point_curves(points, colors, panel.x, panel.y)
        if not plotted:
            plt.close(fig)
            raise KeyError(
                f"no sweep point contains {panel.x.dataset!r} and {panel.y.dataset!r}"
            )
        _draw_overlay_axis(ax, panel, plotted, labelled=False)
        if panel.tag:
            ax.annotate(panel.tag, xy=(0.05, 0.88), xycoords="axes fraction",
                        fontsize=11, fontweight="bold")
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    handles = [Line2D([0.0], [0.0], color=color, lw=1.8) for color in colors]
    labels = [f"{100 * f2:.3g}" for f2, _ in points]
    fig.legend(handles, labels, title=figure.legend_title, loc="upper center",
               ncol=min(len(points), 6), frameon=False, fontsize=9, title_fontsize=9)
    legend_rows = -(-len(points) // 6)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0 - 0.055 * (legend_rows + 1)))
    path = out / figure.filename
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ per-point scalars
def pulse_fwhm_ns(t_ns: NDArray[np.float64], p_w: NDArray[np.float64]) -> float:
    """FWHM (ns) of a pulse trace; 0 when the trace never rises above noise."""
    return pulse_fwhm(t_ns, p_w)


def plot_fwhm_ce(
    sweep_dir: Path | str,
    *,
    filename: str | Path,
    x_range_pct: tuple[float, float] | None = None,
    fwhm_range_ns: tuple[float, float] | None = None,
    ce_range_pct: tuple[float, float] | None = None,
    ce_denominator: str = "deposited",  # or "bank" / "stored"
    out_dir: Path | str | None = None,
) -> Path:
    """Laser-pulse FWHM (ns) and conversion efficiency CE (%) vs F₂ concentration.

    CE = E_laser / E_deposited by default. ``ce_denominator="bank"`` uses Chen's
    ½C₀V₀² denominator; ``"stored"`` uses total initial circuit energy.
    """
    import h5py
    import matplotlib.pyplot as plt

    if ce_denominator not in {"deposited", "bank", "stored"}:
        raise ValueError("ce_denominator must be 'deposited', 'bank', or 'stored'")

    points = sweep_point_dirs(sweep_dir)
    f2_pct, fwhm, ce = [], [], []
    for f2, run_dir in points:
        with h5py.File(run_dir / "diagnostics.h5", "r") as h5:
            if "laser_output" in h5:
                t_ns = h5["laser_output/times_ps"][:] * 1e-3
                fwhm.append(pulse_fwhm_ns(t_ns, h5["laser_output/power_out_w"][:]))
            else:
                fwhm.append(0.0)
            e_out = e_den = 0.0
            if "plasma_scalars" in h5:
                a = h5["plasma_scalars"].attrs
                e_out = float(a.get("e_out_j", 0.0))
                denominator_key = {
                    "deposited": "e_dep_j",
                    "bank": "e_bank0_j",
                    "stored": "e_stored0_j",
                }[ce_denominator]
                e_den = float(a.get(denominator_key, 0.0))
        ce.append(100.0 * e_out / e_den if e_den > 0 else 0.0)
        f2_pct.append(100.0 * f2)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(f2_pct, fwhm, "o-", color="tab:blue", label="FWHM")
    ax.set_xlabel(r"F$_2$ concentration (%)")
    ax.set_ylabel("laser pulse FWHM (ns)")
    ax.set_xlim(*(x_range_pct or data_limits(np.asarray(f2_pct), pad=0.02)))
    ax.set_ylim(*(fwhm_range_ns if fwhm_range_ns is not None else data_limits(
        np.asarray(fwhm), include_zero=True,
    )))
    ax2 = ax.twinx()
    ax2.plot(f2_pct, ce, "s--", color="tab:red", label="CE")
    denominator_label = {
        "deposited": "deposited-energy efficiency",
        "bank": r"bank efficiency $E_L/(\frac{1}{2}C_0V_0^2)$",
        "stored": "initial stored-circuit efficiency",
    }[ce_denominator]
    ax2.set_ylabel(f"{denominator_label} (%)")
    ax2.set_ylim(*(ce_range_pct if ce_range_pct is not None else data_limits(
        np.asarray(ce), include_zero=True,
    )))
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="upper right", fontsize=8)
    fig.tight_layout()
    out = Path(out_dir) if out_dir is not None else Path(sweep_dir)
    path = out / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ α/η vs E/N (Chen Fig. 4)
def plot_alpha_eta(
    base_config,
    f2_values,
    *,
    out_path: Path | str,
    solver=None,
    en_range_td: tuple[float, float] = (1.0, 20.0),
    en_points: int = 40,
    y_range: tuple[float, float] | None = None,
    rebalance: str | None = None,
) -> Path:
    """Townsend ionization α and attachment η coefficients (1/m) vs E/N, F₂-swept.

    Mixture-level (no time integration): for each F₂ value a Boltzmann table is built on a
    linear E/N grid; α = Σ_ion k_r(E/N)·Πn_partners / W and η likewise over net-electron-
    consuming reactions, with W the drift speed and partner densities from the fill at the
    configured pressure/temperature. Log y (a log axis cannot start at 0 — the range floor
    is the first decade), linear x.
    """
    import matplotlib.pyplot as plt

    from particles import excimer_system_registry
    from reactions import load_reaction_network

    from .boltzmann import AnalyticStubSolver, build_lut

    registry = excimer_system_registry()
    compiled = load_reaction_network(base_config.reactions_path, registry).compile()
    solver = solver or AnalyticStubSolver()
    if rebalance is None:
        fr = base_config.gas.normalized_fractions
        rebalance = "Ne" if "Ne" in fr else ("Kr" if "Kr" in fr else max(fr, key=fr.get))
    # net electron production per reaction (from the full network)
    e_net: dict[str, float] = {}
    heavies: dict[str, list[tuple[str, int]]] = {}
    for rxn in tuple(compiled.kinetic) + tuple(compiled.fluid):
        prod_e = sum(c for c, n in rxn.products if n == "e-")
        reac_e = sum(c for c, n in rxn.reactants if n == "e-")
        if reac_e == 0:
            continue  # not electron-impact
        e_net[rxn.reaction_id] = prod_e - reac_e
        heavies[rxn.reaction_id] = [(n, c) for c, n in rxn.reactants if n != "e-"]
    en = np.linspace(*en_range_td, en_points)
    gas0 = base_config.gas
    n_tot = gas0.pressure_atm * ATM_PA / (K_B * gas0.temperature_k)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    colors = _f2_colors(len(f2_values))
    plotted_coefficients: list[NDArray[np.float64]] = []
    for f2, color in zip(f2_values, colors):
        gas = gas0.with_fraction("F2", float(f2), rebalance=rebalance)
        table = build_lut(gas.normalized_fractions, solver=solver,
                          en_grid_td=np.geomspace(max(en_range_td[0], 0.1),
                                                  en_range_td[1], en_points),
                          gas_temp_k=gas.temperature_k, reactions=compiled)
        # Real solvers return mean-energy-keyed tables regardless of the requested keying;
        # interpolating E/N against an <eps> axis silently produces invalid values.
        # Map the E/N axis through the table's local-field-equilibrium <eps>(E/N) first.
        if table.key_kind == "mean_energy_eV":
            key = np.asarray(table.equilibrium_mean_energy(en), dtype=np.float64)
        else:
            key = en
        w = np.maximum(np.interp(key, table.key_grid, table.drift_velocity), 1e-30)
        alpha = np.zeros_like(en)
        eta = np.zeros_like(en)
        for rid, dn_e in e_net.items():
            if rid not in table.rate_coeff or dn_e == 0:
                continue
            k = np.interp(key, table.key_grid, table.rate_coeff[rid])
            n_part = 1.0
            for name, count in heavies[rid]:
                n_part *= (gas.normalized_fractions.get(name, 0.0) * n_tot) ** count
            nu = k * n_part
            if dn_e > 0:
                alpha += dn_e * nu / w
            else:
                eta += (-dn_e) * nu / w
        label = f"{100 * float(f2):.3g}% F$_2$"
        plotted_coefficients.extend((alpha, eta))
        ax.plot(en, np.ma.masked_less_equal(alpha, 0.0), color=color, lw=1.4, label=label)
        ax.plot(en, np.ma.masked_less_equal(eta, 0.0), color=color, lw=1.4, ls="--")
    ax.set_yscale("log")
    ax.set_xlim(*en_range_td)
    ax.set_ylim(*(y_range or data_limits(*plotted_coefficients, scale="log")))
    ax.set_xlabel("E/N (Td)")
    ax.set_ylabel(r"$\alpha$ (solid), $\eta$ (dashed)  (m$^{-1}$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


__all__ = [
    "fit_y_range",
    "OverlayFigure",
    "OverlayPanel",
    "PanelOverlayFigure",
    "plot_alpha_eta",
    "plot_fwhm_ce",
    "plot_sweep_overlays",
    "plot_sweep_panels",
    "pulse_fwhm_ns",
    "sweep_point_dirs",
]
