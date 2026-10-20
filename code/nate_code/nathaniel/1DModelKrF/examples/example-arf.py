#!/usr/bin/env python3
"""
ArF* (193 nm) excimer discharge — replication setup of Chen et al. 2025.

A declarative end-to-end run of the 1D charged-species fluid model configured for the
**ArF** system of Jiayi Chen, Liyang Zhang, Zhihang Zhao & Yangyang Fu, "Fluid modeling of
the discharge in an ArF excimer laser with varying F2 concentration," *Phys. Plasmas* **32**,
113502 (2025). It exercises the real physics engine end-to-end on the ArF gas mixture:

- **Gas fill** Ne/Ar/F₂/Xe at the paper's 96.4/3.5/0.1/0.01 ratio, 3.4 atm, 300 K, in a
  d = 12 mm interelectrode gap (Chen §II.B). ``--f2`` varies the F₂ mole fraction — the
  paper's central study — rebalancing the Ne buffer.
- **Chemistry** a Chen Table III-derived network with documented corrections and cited additions
  (``data/reactions/arf.toml``), with measured electron-impact cross sections.
- **Charged transport** the real 1D drift-diffusion-reaction model (``DriftDiffusionPlasmaModel``):
  Scharfetter–Gummel drift-diffusion + the electron energy PDE, coupled bidirectionally to the
  **MPC discharge circuit** (C₀ = 10 nF pre-charged to −18 kV → L_m → C_p → gap) via the 1D
  Poisson field. Electron transport + rate coefficients come from the Boltzmann LUT keyed on the
  local mean energy (LMEA); ion mobilities from the cited swarm data in ``data/ion_transport/``.
- **Pre-ionization** Chen Eq. 6, ``n_{e,0}(x) = 1e14·tanh(100·x) + 1e13`` m⁻³.
- **Optics** the 0D excimer gain/output model (Chen §II.C, Eqs. 7/10/14/18): the 193 nm cavity
  photon density is evolved by :class:`optics.photon_density.PhotonDensity`, with gain from ArF*
  stimulated emission minus Ar**/F⁻ UV absorption (from ``arf.toml``); the laser power follows
  from the output-coupling efficiency (R_OC = 8 %). So this run produces a **real 193 nm laser
  pulse** — its output energy / efficiency track ArF* build-up during the discharge.

The Boltzmann backend is the project's own **two-term solver** (``TwoTermBoltzmannSolver``,
Hagelaar & Pitchford 2005) over the measured SIGLO cross-sections. It writes static figures
and ``diagnostics.h5`` to ``runs/example_arf/``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    ChargedDensityXT,
    CircuitSpec,
    CircuitWaveforms,
    DataRef,
    EEDF,
    Efficiency,
    ExcimerRates,
    FieldLine,
    GasMixture,
    Grid1D,
    HeatmapPanel,
    HeatmapPlot,
    LaserOutput,
    LinePlot,
    Numerics,
    OpticalGain,
    Plasma,
    PlasmaScalarTraces,
    Plot,
    Preionization,
    ProgressReport,
    Simulation,
    SpeciesDensities,
    Units,
)

# Chen et al. (2025) fill: Ne/Ar/F2/Xe = 96.4 / 3.5 / 0.1 / 0.01 (mole %), 3.4 atm.
_AR_FRAC = 0.035
_XE_FRAC = 0.0001


def _mixture(f2: float) -> dict[str, float]:
    """Ne/Ar/F2/Xe fractions with F2 = ``f2``; the Ne buffer absorbs the balance."""
    ne = 1.0 - _AR_FRAC - _XE_FRAC - f2
    return {"Ne": ne, "Ar": _AR_FRAC, "F2": f2, "Xe": _XE_FRAC}


def _plots(units: Units, *, circuit_backend: str) -> list[Plot]:
    """Every figure choice for this run: names, datasets, scales, and optional limits."""
    density_unit = units.density_label
    length_unit = units.length_label
    time_ps = DataRef("species_densities/times_ps", "time (ns)", scale=1.0e-3)
    plots: list[Plot] = [
        LinePlot(
            filename="electron_distribution.png",
            x=DataRef("eedf/energy_ev", "energy (eV)"),
            y=(DataRef("eedf/f_e", "Maxwellian reconstruction"),),
            ylabel="f(E) (normalized)",
            title="electron EEDF reconstructed from mean energy",
        ),
        LinePlot(
            filename="tracked_species_history.png",
            x=time_ps,
            y=tuple(
                DataRef(f"species_densities/{name}", name)
                for name in ("ArF*", "Ar", "F-")
            ),
            ylabel=f"number density ({density_unit})",
            y_scale="log",
            title="tracked species densities",
        ),
        LinePlot(
            filename="small_signal_gain.png",
            x=DataRef("optical_gain/times_ps", "time (ns)", scale=1.0e-3),
            y=(DataRef(f"optical_gain/gain_per_{length_unit}", "net gain"),),
            ylabel=f"optical gain (1/{length_unit})",
        ),
        LinePlot(
            filename="cavity_output_power.png",
            x=DataRef("laser_output/times_ps", "time (ns)", scale=1.0e-3),
            y=(
                DataRef("laser_output/power_out_spont_w", "spontaneous"),
                DataRef("laser_output/power_out_stim_w", "stimulated"),
                DataRef("laser_output/power_out_w", "total"),
            ),
            ylabel="output power (W)",
            y_scale="log",
        ),
        LinePlot(
            filename="cumulative_efficiency.png",
            x=DataRef("efficiency/times_ps", "time (ns)", scale=1.0e-3),
            y=(DataRef("efficiency/cumulative_efficiency", "intrinsic efficiency"),),
            ylabel="cumulative intrinsic efficiency",
        ),
        LinePlot(
            filename="final_axial_field.png",
            x=DataRef(f"field_line/x_{length_unit}", f"x ({length_unit})"),
            y=(DataRef("field_line/e_mag_v_m", "|E|"),),
            ylabel="electric field magnitude (V/m)",
        ),
        LinePlot(
            filename="final_net_charge.png",
            x=DataRef(f"field_line/x_{length_unit}", f"x ({length_unit})"),
            y=(DataRef("field_line/rho_c_m3", "net charge"),),
            ylabel="net charge density (C/m³)",
        ),
        LinePlot(
            filename="final_spatial_densities.png",
            x=DataRef(f"field_line/x_{length_unit}", f"x ({length_unit})"),
            y=tuple(
                DataRef(f"field_line/n_{name}_per_{length_unit}3", name)
                for name in ("ArF*", "e-", "Ar+")
            ),
            ylabel=f"number density ({density_unit})",
            y_scale="log",
        ),
        HeatmapPlot(
            filename="charged_species_space_time.png",
            x=DataRef(f"charged_density_xt/x_{length_unit}", f"x ({length_unit})"),
            y=DataRef("charged_density_xt/times_ps", "time (ns)", scale=1.0e-3),
            panels=tuple(
                HeatmapPanel(
                    DataRef(f"charged_density_xt/n_{name}", f"n ({density_unit})"),
                    f"n({name})(x, t)",
                )
                for name in ("e-", "Ne+", "Ar+", "F-")
            ),
            color_scale="log",
        ),
        LinePlot(
            filename="mean_electron_density.png",
            x=DataRef("plasma_scalars/times_ns", "time (ns)"),
            y=(DataRef("plasma_scalars/n_e_m3", "mean electron density"),),
            ylabel="electron density (m⁻³)",
            y_scale="log",
        ),
        LinePlot(
            filename="mean_reduced_field.png",
            x=DataRef("plasma_scalars/times_ns", "time (ns)"),
            y=(DataRef("plasma_scalars/en_td", "mean E/N"),),
            ylabel="E/N (Td)",
        ),
        LinePlot(
            filename="mean_electron_temperature.png",
            x=DataRef("plasma_scalars/times_ns", "time (ns)"),
            y=(DataRef("plasma_scalars/te_ev", "mean electron temperature"),),
            ylabel="electron temperature (eV)",
        ),
        LinePlot(
            filename="discharge_power_density.png",
            x=DataRef("plasma_scalars/times_ns", "time (ns)"),
            y=(DataRef("plasma_scalars/p_dis_wm3", "discharge power"),),
            ylabel="power density (W m⁻³)",
        ),
        LinePlot(
            filename="laser_power_density.png",
            x=DataRef("plasma_scalars/times_ns", "time (ns)"),
            y=(DataRef("plasma_scalars/p_l_wm3", "laser power"),),
            ylabel="power density (W m⁻³)",
        ),
        LinePlot(
            filename="arf_formation_rates.png",
            x=DataRef("excimer_rates_formation/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_formation/{reaction}", reaction)
                for reaction in ("R29", "R29c", "R30", "R31")
            ),
            ylabel="formation rate (mol m⁻³ s⁻¹)",
        ),
        LinePlot(
            filename="arf_quenching_rates.png",
            x=DataRef("excimer_rates_quenching/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_quenching/{channel}", channel)
                # SE/STE: optics-owned spontaneous/stimulated sinks (Chen Fig. 7(b)).
                for channel in ("R32", "R33", "R34", "R35", "R36", "SE", "STE")
            ),
            ylabel="quenching rate (mol m⁻³ s⁻¹)",
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--name", type=str, default="example_arf")
    parser.add_argument("--circuit", choices=("mpc", "stub"), default="mpc",
                        help="discharge driver: MPC circuit ODE (default) or a fixed-voltage stub")
    parser.add_argument("--volume-units", choices=("cm", "m"), default="cm")
    parser.add_argument("--cells", type=int, default=256, help="1D mesh cells along the axis")
    parser.add_argument("--f2", type=float, default=0.001,
                        help="F2 mole fraction (Chen's study knob; paper baseline 0.1%%)")
    parser.add_argument("--pressure", type=float, default=3.4, help="gas pressure (atm)")
    parser.add_argument("--voltage", type=float, default=18.0e3, help="C0 charge / electrode voltage (V)")
    parser.add_argument("--dt", type=float, default=None, help="timestep (s); default per circuit")
    parser.add_argument("--steps", type=int, default=None, help="step count; default per circuit")
    parser.add_argument("--fluid-every", type=int, default=20, help="chemistry cadence (steps)")
    parser.add_argument("--report-every", type=int, default=0, help="progress cadence (0 = every 10%%)")
    args = parser.parse_args()

    # Explicit ArF-tuned fields (the default "excimer" field auto-resolves to ArF* here;
    # these add e-); "photons" is the live 193 nm cavity photon count from the 0D
    # gain/output model (stimulated + spontaneous).
    pr_kwargs: dict = {"fields": ("total:ArF*", "total:e-", "temp:e-", "photons")}
    if args.report_every:
        pr_kwargs["every"] = args.report_every
    progress = ProgressReport(**pr_kwargs)

    if args.circuit == "mpc":
        circuit = CircuitSpec(backend="mpc", v0_v=-args.voltage)
        dt = args.dt if args.dt is not None else 50.0e-12
        steps = args.steps if args.steps is not None else 4000  # ~200 ns
    else:
        circuit = CircuitSpec(backend="stub", phi_left=0.0, phi_right=-args.voltage)
        dt = args.dt if args.dt is not None else 1.0e-12
        steps = args.steps if args.steps is not None else 1000

    # Work in the chosen volume units; the Chen values are SI (m, m^-3) expressed in --volume-units.
    units = Units(args.volume_units)
    gap = units.length_from_si(12.0e-3)     # 12 mm interelectrode gap (Chen Table II)
    area = units.area_from_si(90.0e-4)       # 90 cm^2 discharge area (Chen Table II)

    # t=0 electron pre-ionization (Chen Eq. 6), SI values expressed in --volume-units.
    preion = Preionization(mode="tanh",
                           amplitude=units.density_from_si(1.0e14),
                           floor=units.density_from_si(1.0e13),
                           steepness=units.inv_length_from_si(100.0))

    diagnostics = [
        EEDF(energy_range=(0.1, 40.0), bins=80),
        SpeciesDensities(("ArF*", "Ar", "F-")),
        OpticalGain(),
        LaserOutput(),
        Efficiency(),
        FieldLine(("ArF*", "e-", "Ar+")),
        ChargedDensityXT(("e-", "Ne+", "Ar+", "F-")),
        PlasmaScalarTraces(),
        ExcimerRates(excimer="ArF*"),
    ]
    if args.circuit == "mpc":
        diagnostics.append(CircuitWaveforms())

    sim = Simulation(
        name=args.name,
        volume_units=args.volume_units,
        grid=Grid1D(cells=args.cells, length=gap, area_m2=area),
        gas=GasMixture(_mixture(args.f2), pressure_atm=args.pressure, temperature_k=300.0),
        # Majority-gas ion balances the pre-ionization electrons quasi-neutrally at t=0.
        # NOTE: with `preionization` set, the electron profile comes from Chen Eq. 6 alone —
        # `density` here is only the fallback used when no Preionization spec is given.
        plasma=Plasma(density=units.density_from_si(1.0e16), ion="Ne+", electron_mean_energy_ev=2.0),
        preionization=preion,
        # ArF*; R_OC=8% (Chen Table II); HR mirror 100% ("assuming 100% reflectivity", §II.C)
        cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),
        # Electron coefficients: the two-term Boltzmann solve of the measured SIGLO
        # cross-sections, including the non-Maxwellian tail response.
        boltzmann=BoltzmannSpec(backend="two_term", key_kind="mean_energy_eV"),
        circuit=circuit,
        reactions=_ROOT / "data" / "reactions" / "arf.toml",
        output_dir=_ROOT / "runs",
        diagnostics=diagnostics,
        plots=_plots(units, circuit_backend=args.circuit),
    )
    sim.run(
        Numerics(steps=steps, dt=dt, fluid_every=args.fluid_every),
        progress=progress,
    )


if __name__ == "__main__":
    main()
