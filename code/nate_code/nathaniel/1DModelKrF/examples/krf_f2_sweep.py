#!/usr/bin/env python3
"""Ne/Kr/F2 KrF* concentration sweep and Chen-style diagnostic overlays.

This is the 248 nm mirror of ``arf_f2_sweep.py``. It runs one process per F2 point with
the 45-reaction ``krf.toml`` network, writes each point to
``runs/<name>/f2_<fraction>/diagnostics.h5``, and renders:

- circuit waveforms and KrF* formation/quenching rates (including the optics-owned SE/STE
  radiative sinks, Chen Fig. 7 style) for every mixture;
- electron density, E/N, electron temperature, discharge power, and laser power overlays;
- the Kr⁺/Kr*/KrF*/F⁻ density panel figure (Chen Fig. 15 analog), one curve per F2 value;
- laser FWHM and conversion efficiency versus F2 concentration; and
- Townsend ionization/attachment coefficients versus E/N for every mixture.

The production default scans 0.05--0.50% F2 in a 96%-class Ne buffer with 3.5% Kr. Use
``--tiny`` for a three-point, stub-backed development smoke test.
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
_PANEL_SPECIES = ("Kr+", "Kr*", "KrF*", "F-")
_PANEL_TIME = DataRef("species_densities/times_ps", "t (ns)", scale=1.0e-3)
# Chen Fig. 15 analog: 2x2 density panels, one curve per F2 value, shared legend.
SPECIES_PANEL_FIGURE = PanelOverlayFigure(
    filename="sweep_species_densities.png",
    panels=(
        OverlayPanel(tag="(a)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/Kr+",
                               r"$n_{\mathrm{Kr}^+}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
        OverlayPanel(tag="(b)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/Kr*",
                               r"$n_{\mathrm{Kr}^*}$ ($\times 10^{21}$ m$^{-3}$)",
                               scale=1.0e-21)),
        OverlayPanel(tag="(c)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/KrF*",
                               r"$n_{\mathrm{KrF}^*}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
        OverlayPanel(tag="(d)", x=_PANEL_TIME, x_include_zero=True, y_include_zero=True,
                     y=DataRef("species_densities/F-",
                               r"$n_{\mathrm{F}^-}$ ($\times 10^{20}$ m$^{-3}$)",
                               scale=1.0e-20)),
    ),
)

_KR_FRAC = 0.035


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--name", default="krf_f2_sweep")
    parser.add_argument("--start", type=float, default=0.0005, help="F2 fraction sweep start")
    parser.add_argument("--stop", type=float, default=0.005, help="F2 fraction sweep stop")
    parser.add_argument("--step", type=float, default=0.0001, help="F2 fraction sweep step")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--kr", type=float, default=_KR_FRAC,
                        help="fixed Kr mole fraction; Ne absorbs the F2 sweep delta")
    parser.add_argument("--excimer", default="KrF*", help="excimer species for rate figures")
    parser.add_argument(
        "--reactions", default=str(_ROOT / "data" / "reactions" / "krf.toml")
    )
    parser.add_argument(
        "--tiny", action="store_true",
        help="dev-test config: stub backend+circuit, 24 cells, about 2 ns, 3 points",
    )
    args = parser.parse_args()

    if not 0.0 < args.kr < 1.0:
        parser.error("--kr must be between zero and one")
    baseline_f2 = 0.001
    ne_fraction = 1.0 - args.kr - baseline_f2
    if ne_fraction <= 0.0:
        parser.error("--kr plus the 0.1% baseline F2 fraction must be less than one")
    mixture = {"Ne": ne_fraction, "Kr": args.kr, "F2": baseline_f2}

    if args.tiny:
        grid = Grid1D(cells=24, length=12.0e-3, area_m2=90.0e-4)
        boltzmann = BoltzmannSpec(backend="stub")
        circuit = CircuitSpec(backend="stub", phi_left=0.0, phi_right=-5.0e3)
        numerics = Numerics(steps=40, dt=50.0e-12, fluid_every=10)
        selection = {"values": (0.0005, 0.001, 0.002)}
    else:
        grid = Grid1D(cells=512, length=12.0e-3, area_m2=90.0e-4)
        boltzmann = BoltzmannSpec(backend="two_term", key_kind="mean_energy_eV")
        circuit = CircuitSpec(backend="mpc", v0_v=-18.0e3)
        numerics = Numerics(steps=4000, dt=50.0e-12, fluid_every=20)
        selection = {"start": args.start, "stop": args.stop, "step": args.step}

    base = Simulation(
        name=args.name,
        grid=grid,
        gas=GasMixture(mixture, pressure_atm=3.4, temperature_k=300.0),
        plasma=Plasma(density=1.0e16, ion="Ne+", electron_mean_energy_ev=2.0),
        preionization=Preionization(
            mode="tanh", amplitude=1.0e14, floor=1.0e13, steepness=100.0
        ),
        # F2, F-, and Kr** are already explicit PHOTO absorbers in krf.toml. Setting k_c=0
        # prevents Chen's empirical F2 loss term from attenuating the same molecules twice.
        cavity=Cavity(wavelength_nm=248.0, r_out=0.1, r_back=1.0, loss_coeff_per_f2=0.0),
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
        system="krf",
        formation_ids=("R28", "R29", "R30", "R31", "R32"),
        quenching_ids=("R33", "R34", "R35", "R36", "R37", "SE", "STE"),
        circuit_backend=circuit.backend,
    )
    result = sweep_f2_parallel(
        base,
        **selection,
        numerics=numerics,
        rebalance="Ne",
        max_workers=args.workers,
        diagnostics_factory=factory,
        plots_factory=plot_factory,
    )

    sweep_dir = _ROOT / "runs" / args.name
    figures = render_standard_sweep_figures(
        base, result.points, sweep_dir, SPECIES_PANEL_FIGURE
    )

    print(f"sweep of {len(result.points)} F2 points -> {sweep_dir}")
    for point in result.points:
        print(
            f"  F2={100 * point.f2_fraction:.3g}%  "
            f"E_out={point.output_energy_j:.3e} J  "
            f"eff={point.intrinsic_efficiency:.3e}"
        )
    for figure in figures:
        print(f"  fig: {figure}")


if __name__ == "__main__":
    main()
