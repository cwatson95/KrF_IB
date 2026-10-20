#!/usr/bin/env python3
"""Single Ne/Kr/F2 KrF* discharge with the 45-reaction reduced chemistry.

A declarative end-to-end run of the complete 1D fluid model. The
experiment is composed from component specs (1D grid / gas / plasma / cavity / Boltzmann
source / discharge circuit) plus a list of diagnostics — change any of them here without
touching the engine:

- **Charged species** (e⁻, Kr⁺, F⁻, …) are a 1D fluid; electron transport + electron-impact
  rate coefficients come from the project's two-term Boltzmann solve of the cited
  Ne/Kr/F₂ cross sections. Use the sweep script's ``--tiny`` mode for a stub-backed smoke test.
- **Neutral species** (Ne/Kr/F₂ fill, plus Kr*, Ne*, KrF*, … produced during the run) are a
  per-cell fluid; the combined charged+neutral chemistry is integrated on the fluid cadence.
- **Photons** (248 nm) follow the 0D excimer gain/output model (GainCavity/PhotonDensity):
  seeded by KrF* → Kr + F + hν decay, amplified by stimulated gain, and attenuated by the
  F₂ / F⁻ / Kr** photoabsorption channels declared in ``krf.toml``.
- **Field** comes from a 1D Poisson solve with Dirichlet electrode potentials from a
  discharge circuit — by default the ``MpcCircuit`` magnetic-pulse-compression driver
  (Chen et al. 2025) whose 7-state ODE ramps the cathode over ~100 ns; ``--circuit stub``
  applies a fixed voltage instead.

Writes static figures + ``diagnostics.h5`` to ``runs/01_krf_excimer/``. The matching
production F₂ sweep is ``examples/krf_f2_sweep.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    Absorption,
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

# Diagnostic mirror of krf.toml R43--R45. The cavity itself reads these channels directly from
# the network, so Cavity.absorption stays empty and does not double-count them.
PHOTOABSORPTION_248 = {"Kr**": 6.0e-22, "F-": 5.6e-22, "F2": 1.48e-24}
_KR_FRAC = 0.035  # match the active-rare-gas fraction in the Chen ArF comparison


def _mixture(f2: float, kr: float) -> dict[str, float]:
    """Ne-buffered three-gas fill, with F2 and Kr explicit and Ne taking the balance."""
    ne = 1.0 - kr - f2
    if ne <= 0.0:
        raise ValueError("Kr and F2 fractions must sum to less than one")
    return {"Ne": ne, "Kr": kr, "F2": f2}


def _plots(units: Units, *, circuit_backend: str) -> list[Plot]:
    """Every figure choice for this run: names, datasets, scales, and optional limits."""
    density_unit = units.density_label
    length_unit = units.length_label
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
            x=DataRef("species_densities/times_ps", "time (ns)", scale=1.0e-3),
            y=tuple(
                DataRef(f"species_densities/{name}", name)
                for name in ("F", "F-", "Ne*", "Kr*", "Kr**", "KrF*", "Kr2F*", "e-")
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
            filename="photoabsorption.png",
            x=DataRef("absorption/times_ps", "time (ns)", scale=1.0e-3),
            y=(DataRef(f"absorption/alpha_per_{length_unit}", "absorption"),),
            ylabel=f"absorption coefficient (1/{length_unit})",
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
                for name in ("KrF*", "e-", "Kr*")
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
                for name in ("e-", "Ne+", "Kr+", "F-")
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
            filename="krf_formation_rates.png",
            x=DataRef("excimer_rates_formation/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_formation/{reaction}", reaction)
                for reaction in ("R28", "R29", "R30", "R31", "R32")
            ),
            ylabel="formation rate (mol m⁻³ s⁻¹)",
        ),
        LinePlot(
            filename="krf_quenching_rates.png",
            x=DataRef("excimer_rates_quenching/times_ns", "time (ns)"),
            y=tuple(
                DataRef(f"excimer_rates_quenching/{channel}", channel)
                # SE/STE: optics-owned spontaneous/stimulated sinks (Chen Fig. 7(b)).
                for channel in ("R33", "R34", "R35", "R36", "R37", "SE", "STE")
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
    parser.add_argument("--name", type=str, default="01_krf_excimer")
    parser.add_argument("--circuit", choices=("mpc", "stub"), default="mpc",
                        help="discharge driver: MPC circuit ODE (default) or a fixed-voltage stub")
    parser.add_argument("--volume-units", choices=("cm", "m"), default="cm",
                        help="volume/length unit system for inputs + reported densities/lengths")
    parser.add_argument("--cells", type=int, default=128, help="1D mesh cells along the axis")
    parser.add_argument("--length", type=float, default=None,
                        help="electrode gap (in --volume-units; default 12 mm)")
    parser.add_argument("--area", type=float, default=None,
                        help="discharge area (in --volume-units²; default 90 cm²)")
    parser.add_argument("--pressure", type=float, default=3.4, help="gas pressure (atm)")
    parser.add_argument("--f2", type=float, default=0.001, help="initial F2 mole fraction")
    parser.add_argument("--kr", type=float, default=_KR_FRAC,
                        help="Kr mole fraction (default 3.5%%; Ne takes the balance)")
    parser.add_argument("--preionization", choices=("tanh", "uniform"), default="tanh",
                        help="t=0 electron pre-ionization profile n_e0(x) (default: Chen et al. tanh)")
    parser.add_argument("--voltage", type=float, default=18.0e3, help="C0 charge / electrode voltage (V)")
    parser.add_argument("--dt", type=float, default=None, help="timestep (s); default per circuit")
    parser.add_argument("--steps", type=int, default=None, help="step count; default per circuit")
    parser.add_argument("--fluid-every", type=int, default=20, help="chemistry+photon cadence (steps)")
    parser.add_argument("--report-every", type=int, default=0, help="progress cadence (0 = every 10%%)")
    args = parser.parse_args()

    progress_kwargs: dict = {"fields": ("total:KrF*", "total:e-", "temp:e-", "photons")}
    if args.report_every:
        progress_kwargs["every"] = args.report_every
    progress = ProgressReport(**progress_kwargs)

    # The MPC driver ramps the gap over ~100 ns; the stub applies a fixed voltage immediately.
    if args.circuit == "mpc":
        circuit = CircuitSpec(backend="mpc", v0_v=-args.voltage)
        dt = args.dt if args.dt is not None else 50.0e-12
        steps = args.steps if args.steps is not None else 4000  # ~200 ns
    else:
        circuit = CircuitSpec(backend="stub", phi_left=0.0, phi_right=-args.voltage)
        dt = args.dt if args.dt is not None else 1.0e-12
        steps = args.steps if args.steps is not None else 1000

    # Work in the chosen volume units: the Chen et al. paper values are SI (m^-3, m); express them
    # in --volume-units so the specs + all reported densities/lengths use it. Physics is identical
    # (Simulation converts back to SI internally) — this just sets the user's working unit.
    units = Units(args.volume_units)
    gap = args.length if args.length is not None else units.length_from_si(12.0e-3)   # 12 mm
    area = args.area if args.area is not None else units.area_from_si(90.0e-4)         # 90 cm²

    # t=0 electron pre-ionization profile. Chen et al. Eq. 6 (SI): n_e0(x)=1e14*tanh(100x)+1e13 m^-3,
    # rising from the anode (x=0) toward the cathode (x=length); expressed here in --volume-units.
    if args.preionization == "tanh":
        preion = Preionization(mode="tanh",
                               amplitude=units.density_from_si(1.0e14),
                               floor=units.density_from_si(1.0e13),
                               steepness=units.inv_length_from_si(100.0))
    else:
        preion = Preionization(mode="uniform", density=units.density_from_si(1.0e16))

    diagnostics = [
        EEDF(energy_range=(0.1, 40.0), bins=80),
        SpeciesDensities(("F", "F-", "Ne*", "Kr*", "Kr**", "KrF*", "Kr2F*", "e-")),
        OpticalGain(),
        Absorption(PHOTOABSORPTION_248),
        LaserOutput(),
        Efficiency(),
        FieldLine(("KrF*", "e-", "Kr*")),
        ChargedDensityXT(("e-", "Ne+", "Ne2+", "Kr+", "Kr2+", "F-", "F2+")),
        PlasmaScalarTraces(),
        ExcimerRates(excimer="KrF*"),
    ]
    if args.circuit == "mpc":
        diagnostics.append(CircuitWaveforms())

    sim = Simulation(
        name=args.name,
        volume_units=args.volume_units,
        grid=Grid1D(cells=args.cells, length=gap, area_m2=area),
        gas=GasMixture(_mixture(args.f2, args.kr), pressure_atm=args.pressure,
                       temperature_k=300.0),
        # Ne is the majority gas, so Ne+ balances the pre-ionization seed.
        plasma=Plasma(density=units.density_from_si(1.0e16), ion="Ne+",
                      electron_mean_energy_ev=2.0),
        preionization=preion,
        # R43--R45 already contribute explicit species-resolved absorption. Disable Chen's
        # empirical k_c*x_F2 loss term here so F2 is not attenuated twice.
        # r_back=1.0 aligns with Chen's stated 100% HR (kept consistent with the sweeps).
        cavity=Cavity(wavelength_nm=248.0, r_out=0.1, r_back=1.0,
                      loss_coeff_per_f2=0.0),
        boltzmann=BoltzmannSpec(backend="two_term", key_kind="mean_energy_eV"),
        circuit=circuit,
        reactions=_ROOT / "data" / "reactions" / "krf.toml",
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
