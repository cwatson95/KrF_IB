#!/usr/bin/env python3
"""ArF* F₂-concentration sweep + Chen-figure overlays — the project's headline experiment.

Runs `sweep_f2_parallel` (one process per F₂ point) over the Chen et al. (2025) ArF
configuration, writing per-point ``runs/<name>/f2_<value>/diagnostics.h5``, then renders the
sweep figures. Every filename, x/y dataset, and scale is declared here. Explicit ranges win;
otherwise axes tightly fit all finite visible data with a small margin:

- per point (inside each f2_* dir): circuit waveforms U_gap/U_C0/I_d (Chen Fig. 6),
  ArF* formation + quenching reaction rates in mol m⁻³ s⁻¹ including the optics-owned
  SE/STE radiative sinks (Chen Fig. 7);
- sweep overlays (one curve per F₂): n_e, E/N, T_e vs t (Chen Fig. 10) and P_dis, P_L
  power densities vs t (Chen Fig. 12), plus the Ar⁺/Ar*/ArF*/F⁻ density panel figure
  (Chen Fig. 15);
- per-point scalars: laser FWHM (ns) + conversion efficiency (%) vs F₂ (Chen Figs. 16/17);
- mixture-level: Townsend α/η vs E/N, F₂-swept, log y (Chen Fig. 9).

Dev testing: ``--tiny`` shrinks to a stub-backend 24-cell few-ns configuration (full-scale
sweeps are production runs — see CLAUDE.md).
"""

from __future__ import annotations

import argparse
from functools import partial
from pathlib import Path

from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    DataRef,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    Preionization,
    Simulation,
)
from Fluid1D.sweep import sweep_f2_parallel
from Fluid1D.sweep_plots import OverlayPanel, PanelOverlayFigure
from Fluid1D.sweep_presets import (
    render_standard_sweep_figures,
    standard_point_diagnostics,
    standard_point_plots,
)

_ROOT = Path(__file__).resolve().parents[1]
_PANEL_SPECIES = ("Ar+", "Ar*", "ArF*", "F-")
_PANEL_TIME = DataRef("species_densities/times_ps", "t (ns)", scale=1.0e-3)
# Chen Fig. 15: 2x2 density panels, one curve per F₂ value, shared legend.
SPECIES_PANEL_FIGURE = PanelOverlayFigure(
    filename="sweep_species_densities.png",
    panels=(
        OverlayPanel(tag="(a)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/Ar+",
                               r"$n_{\mathrm{Ar}^+}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
        OverlayPanel(tag="(b)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/Ar*",
                               r"$n_{\mathrm{Ar}^*}$ ($\times 10^{21}$ m$^{-3}$)",
                               scale=1.0e-21)),
        OverlayPanel(tag="(c)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/ArF*",
                               r"$n_{\mathrm{ArF}^*}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
        OverlayPanel(tag="(d)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/F-",
                               r"$n_{\mathrm{F}^-}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
    ),
)

_AR_FRAC, _XE_FRAC = 0.035, 0.0001


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--name", default="arf_f2_sweep")
    p.add_argument("--start", type=float, default=0.0005, help="F2 fraction sweep start")
    p.add_argument("--stop", type=float, default=0.005, help="F2 fraction sweep stop")
    p.add_argument("--step", type=float, default=0.0001, help="F2 fraction sweep step")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--excimer", default="ArF*", help="excimer species for the rate figures")
    p.add_argument("--reactions", default=str(_ROOT / "data" / "reactions" / "arf.toml"))
    p.add_argument("--tiny", action="store_true",
                   help="dev-test config: stub backend+circuit, 24 cells, ~2 ns, 3 points")
    args = p.parse_args()

    mixture = {"Ne": 1.0 - _AR_FRAC - _XE_FRAC - 0.001, "Ar": _AR_FRAC, "F2": 0.001,
               "Xe": _XE_FRAC}
    if args.tiny:
        grid = Grid1D(cells=24, length=12e-3, area_m2=90e-4)
        boltzmann = BoltzmannSpec(backend="stub")
        circuit = CircuitSpec(backend="stub", phi_left=0.0, phi_right=-5.0e3)
        numerics = Numerics(steps=40, dt=50e-12, fluid_every=10)
        selection = {"values": (0.0005, 0.001, 0.002)}
    else:
        grid = Grid1D(cells=512, length=12e-3, area_m2=90e-4)
        boltzmann = BoltzmannSpec(backend="two_term", key_kind="mean_energy_eV")
        circuit = CircuitSpec(backend="mpc", v0_v=-18.0e3)
        numerics = Numerics(steps=4000, dt=50e-12, fluid_every=20)
        selection = {"start": args.start, "stop": args.stop, "step": args.step}

    base = Simulation(
        name=args.name,
        grid=grid,
        gas=GasMixture(mixture, pressure_atm=3.4, temperature_k=300.0),
        plasma=Plasma(density=1.0e16, ion="Ne+", electron_mean_energy_ev=2.0),
        preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0),
        cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),  # Chen: 100% HR
        boltzmann=boltzmann,
        circuit=circuit,
        reactions=args.reactions,
        output_dir=_ROOT / "runs",
        diagnostics=[],
    )
    factory = partial(
        standard_point_diagnostics,
        excimer=args.excimer,
        circuit_backend=circuit.backend,
        panel_species=_PANEL_SPECIES,
    )
    plot_factory = partial(
        standard_point_plots,
        system="arf",
        formation_ids=("R29", "R29c", "R30", "R31"),
        quenching_ids=("R32", "R33", "R34", "R35", "R36", "SE", "STE"),
        circuit_backend=circuit.backend,
    )
    result = sweep_f2_parallel(
        base, **selection,
        numerics=numerics, rebalance="Ne", max_workers=args.workers,
        diagnostics_factory=factory,
        plots_factory=plot_factory,
    )

    sweep_dir = _ROOT / "runs" / args.name
    figs = render_standard_sweep_figures(
        base, result.points, sweep_dir, SPECIES_PANEL_FIGURE
    )
    print(f"sweep of {len(result.points)} F2 points -> {sweep_dir}")
    for pt in result.points:
        print(f"  F2={100 * pt.f2_fraction:.3g}%  E_out={pt.output_energy_j:.3e} J  "
              f"eff={pt.intrinsic_efficiency:.3e}")
    for f in figs:
        print(f"  fig: {f}")


if __name__ == "__main__":
    main()
