"""Shared diagnostics and figure declarations for ArF and KrF F₂ sweeps."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from diagnostics import DataRef, LinePlot, Plot

from .diagnostics import (
    CircuitWaveforms,
    ExcimerRates,
    LaserOutput,
    PlasmaScalarTraces,
    SpeciesDensities,
)
from .simulation import Simulation
from .sweep import SweepPoint
from .sweep_plots import (
    OverlayFigure,
    PanelOverlayFigure,
    plot_alpha_eta,
    plot_fwhm_ce,
    plot_sweep_overlays,
    plot_sweep_panels,
)


STANDARD_OVERLAY_FIGURES = (
    OverlayFigure(
        filename="sweep_electron_density.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/n_e_m3", r"$n_e$ (m$^{-3}$)"),
        y_scale="log",
    ),
    OverlayFigure(
        filename="sweep_reduced_field.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/en_td", "E/N (Td)"),
    ),
    OverlayFigure(
        filename="sweep_electron_temperature.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/te_ev", r"$T_e$ (eV)"),
    ),
    OverlayFigure(
        filename="sweep_discharge_power.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/p_dis_wm3", r"$P_{dis}$ (W m$^{-3}$)"),
    ),
    OverlayFigure(
        filename="sweep_laser_power.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/p_l_wm3", r"$P_L$ (W m$^{-3}$)"),
    ),
)

STANDARD_FWHM_CE_OPTIONS = {
    "filename": "sweep_pulse_width_efficiency.png",
    "x_range_pct": (0.0, 0.54),
    "ce_denominator": "bank",
}

STANDARD_ALPHA_ETA_RANGES = {
    "en_range_td": (1.0, 20.0),
    "y_range": (1.0, 1.0e5),
}


def standard_point_diagnostics(
    f2: float,
    *,
    excimer: str,
    circuit_backend: str,
    panel_species: tuple[str, ...],
) -> list[object]:
    """Return the numeric collectors shared by both production sweep scripts."""
    del f2
    diagnostics: list[object] = [
        PlasmaScalarTraces(),
        ExcimerRates(excimer=excimer),
        LaserOutput(),
        SpeciesDensities(track=panel_species),
    ]
    if circuit_backend == "mpc":
        diagnostics.insert(1, CircuitWaveforms())
    return diagnostics


def standard_point_plots(
    f2: float,
    *,
    system: str,
    formation_ids: tuple[str, ...],
    quenching_ids: tuple[str, ...],
    circuit_backend: str,
) -> list[Plot]:
    """Return per-point rate, laser-power, and optional MPC waveform figures."""
    del f2
    plots: list[Plot] = [
        LinePlot(
            filename=f"{system}_formation_rates.png",
            x=DataRef("excimer_rates_formation/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_formation/{reaction}", reaction)
                for reaction in formation_ids
            ),
            ylabel="formation rate (mol m⁻³ s⁻¹)",
        ),
        LinePlot(
            filename=f"{system}_quenching_rates.png",
            x=DataRef("excimer_rates_quenching/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_quenching/{channel}", channel)
                for channel in quenching_ids
            ),
            ylabel="quenching rate (mol m⁻³ s⁻¹)",
        ),
        LinePlot(
            filename="laser_output_power.png",
            x=DataRef("laser_output/times_ps", "time (ns)", scale=1.0e-3),
            y=(
                DataRef("laser_output/power_out_spont_w", "spontaneous"),
                DataRef("laser_output/power_out_stim_w", "stimulated"),
                DataRef("laser_output/power_out_w", "total"),
            ),
            ylabel="output power (W)",
            y_scale="log",
        ),
    ]
    if circuit_backend == "mpc":
        plots.append(
            LinePlot(
                filename="mpc_circuit_waveforms.png",
                x=DataRef("circuit_waveforms/times_ns", "time (ns)"),
                y=(
                    DataRef("circuit_waveforms/u_gap_kv", "U_gap"),
                    DataRef("circuit_waveforms/u_c0_kv", "U_C0"),
                ),
                ylabel="voltage (kV)",
                secondary_y=(DataRef("circuit_waveforms/i_d_ka", "I_d"),),
                secondary_ylabel="discharge current (kA)",
                x_limits=(0.0, 200.0),
            )
        )
    return plots


def render_standard_sweep_figures(
    base: Simulation,
    points: Sequence[SweepPoint],
    sweep_dir: Path,
    species_panels: PanelOverlayFigure,
) -> tuple[Path, ...]:
    """Render the aggregate figures shared by the ArF and KrF sweep scripts."""
    figures = list(plot_sweep_overlays(sweep_dir, STANDARD_OVERLAY_FIGURES))
    figures.append(plot_sweep_panels(sweep_dir, species_panels))
    figures.append(plot_fwhm_ce(sweep_dir, **STANDARD_FWHM_CE_OPTIONS))
    figures.append(
        plot_alpha_eta(
            base,
            [point.f2_fraction for point in points],
            **STANDARD_ALPHA_ETA_RANGES,
            out_path=sweep_dir / "alpha_eta_en.png",
        )
    )
    return tuple(figures)


__all__ = (
    "STANDARD_ALPHA_ETA_RANGES",
    "STANDARD_FWHM_CE_OPTIONS",
    "STANDARD_OVERLAY_FIGURES",
    "render_standard_sweep_figures",
    "standard_point_diagnostics",
    "standard_point_plots",
)
