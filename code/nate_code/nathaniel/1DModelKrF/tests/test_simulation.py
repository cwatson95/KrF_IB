"""End-to-end smoke test of the declarative 1D fluid Simulation + diagnostics (stub-backed)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless plotting

_ROOT = Path(__file__).resolve().parents[1]
import h5py
import numpy as np

from Fluid1D import (
    Absorption,
    BoltzmannSpec,
    Cavity,
    ChargedDensityXT,
    CircuitSpec,
    DataRef,
    EEDF,
    Efficiency,
    FieldLine,
    GasMixture,
    Grid1D,
    LaserOutput,
    LinePlot,
    Numerics,
    OpticalGain,
    Plasma,
    Preionization,
    Simulation,
    SpeciesDensities,
)

KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"


def _sim(diagnostics, name="smoke", output_dir=None, plots=None):
    return Simulation(
        name=name,
        grid=Grid1D(cells=4, length=2.0e-2, area_m2=1.0e-4),
        gas=GasMixture({"Kr": 0.92, "F2": 0.08}),
        plasma=Plasma(density=1.0e18, ion="Kr+"),
        cavity=Cavity(),
        boltzmann=BoltzmannSpec(backend="stub"),
        circuit=CircuitSpec(),
        reactions=KRF_TOML,
        diagnostics=diagnostics,
        plots=plots,
        output_dir=Path(output_dir or tempfile.gettempdir()),
    )


class TestSimulation(unittest.TestCase):
    def test_field_heating_cap_control_validation(self) -> None:
        self.assertIsNone(Numerics(field_heating_cap_factor=None).field_heating_cap_factor)
        for value in (0.0, -1.0, np.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Numerics(field_heating_cap_factor=value)

    def test_short_run_flushes_partial_fluid_cadence(self) -> None:
        sim = _sim([], name="partial-cadence")
        sim.run(Numerics(steps=1, dt=1e-12, fluid_every=10), verbose=False)
        self.assertEqual(len(sim.energetics.times_s), 1)
        self.assertEqual(sim.energetics.interval_s, [1e-12])
        self.assertGreaterEqual(sim.energetics.input_energy_j(), 0.0)

    def test_initial_ion_must_be_positive_and_known(self) -> None:
        for ion in ("F-", "not-a-species"):
            sim = Simulation(
                name="bad-ion", grid=Grid1D(cells=4, length=2e-2, area_m2=1e-4),
                gas=GasMixture({"Kr": 1.0}), plasma=Plasma(density=1e18, ion=ion),
                cavity=Cavity(), reactions=KRF_TOML, diagnostics=[],
            )
            with self.subTest(ion=ion), self.assertRaises(ValueError):
                sim.run(Numerics(steps=1, dt=1e-12), verbose=False)

    def test_declarative_diagnostics_produce_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            plots = [
                LinePlot("custom_eedf.png", DataRef("eedf/energy_ev", "energy (eV)"),
                         (DataRef("eedf/f_e", "f(E)"),), "f(E)"),
                LinePlot("custom_species.png", DataRef("species_densities/times_ps", "time (ps)"),
                         tuple(DataRef(f"species_densities/{name}", name)
                               for name in ("Kr", "KrF*", "e-")), "density", y_scale="log"),
                LinePlot("custom_gain.png", DataRef("optical_gain/times_ps", "time (ps)"),
                         (DataRef("optical_gain/gain_per_m", "gain"),), "gain (1/m)"),
                LinePlot("custom_absorption.png", DataRef("absorption/times_ps", "time (ps)"),
                         (DataRef("absorption/alpha_per_m", "absorption"),), "absorption (1/m)"),
                LinePlot("custom_output.png", DataRef("laser_output/times_ps", "time (ps)"),
                         (DataRef("laser_output/power_out_w", "output"),), "power (W)"),
                LinePlot("custom_efficiency.png", DataRef("efficiency/times_ps", "time (ps)"),
                         (DataRef("efficiency/cumulative_efficiency", "efficiency"),),
                         "efficiency"),
                LinePlot("custom_field.png", DataRef("field_line/x_m", "x (m)"),
                         (DataRef("field_line/e_mag_v_m", "|E|"),), "field (V/m)"),
            ]
            sim = _sim(
                [EEDF(at_step=16), SpeciesDensities(("Kr", "KrF*", "e-")),
                 OpticalGain(), Absorption({"F2": 4.5e-24}), LaserOutput(),
                 Efficiency(), FieldLine(("KrF*", "e-"))],
                output_dir=d, plots=plots,
            )
            out = sim.run(Numerics(steps=16, dt=1e-12, fluid_every=4), verbose=False)
            produced = {p.name for p in out.iterdir()}
            for expected in ("custom_eedf.png", "custom_species.png", "custom_gain.png",
                             "custom_absorption.png", "custom_output.png",
                             "custom_efficiency.png", "custom_field.png", "diagnostics.h5"):
                self.assertIn(expected, produced, msg=f"missing {expected}")
            for source_owned_name in ("eedf.png", "species_densities.png", "optical_gain.png"):
                self.assertNotIn(source_owned_name, produced)
            # No animation / 3D-volumetric artifacts in the 1D model.
            for forbidden in ("fields.h5", "fields.xmf", "field_slices.png",
                              "charged_and_photons_3d.mp4", "animation_frames"):
                self.assertNotIn(forbidden, produced, msg=f"unexpected {forbidden}")

    def test_numerics_mirrored_to_archive(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            sim = _sim(
                [EEDF(at_step=8), SpeciesDensities(("Kr", "KrF*")), OpticalGain(),
                 Absorption({"F2": 4.5e-24}), LaserOutput(), Efficiency()],
                name="archive", output_dir=d,
            )
            out = sim.run(Numerics(steps=8, dt=1e-12, fluid_every=2), verbose=False)
            with h5py.File(out / "diagnostics.h5", "r") as f:
                for grp in ("eedf", "species_densities", "optical_gain",
                            "absorption", "laser_output", "efficiency"):
                    self.assertIn(grp, f, msg=f"missing archive group {grp}")
                self.assertIn("Kr", f["species_densities"])
                # EEDF is now a smooth Maxwellian curve (energy_ev, f_e), not a histogram.
                self.assertEqual(f["eedf/energy_ev"].shape, f["eedf/f_e"].shape)
                self.assertGreater(f["eedf"].attrs["te_ev"], 0.0)
                for col in ("power_out_w", "power_out_spont_w", "power_out_stim_w",
                            "power_in_w", "cumulative_efficiency"):
                    self.assertIn(col, f["efficiency"], msg=f"missing efficiency column {col}")

    def test_photons_reachable_via_species_total(self) -> None:
        sim = _sim([SpeciesDensities(("hv248",))], name="photons")
        sim.run(Numerics(steps=60, dt=1e-12, fluid_every=10), verbose=False)
        # The 0D gain model seeds n_ph from spontaneous emission, so the cavity count is > 0.
        self.assertGreater(sim.gain.total_photons(), 0.0)
        self.assertEqual(sim.species_total("hv248"), sim.gain.total_photons())
        self.assertEqual(sim.species_total("hv193"), 0.0)  # mode not modeled -> 0, not a crash

    def _units_sim(self, volume_units, density):
        return Simulation(
            name="units", grid=Grid1D(cells=4, length=2.0e-2, area_m2=1.0e-4),
            gas=GasMixture({"Kr": 0.92, "F2": 0.08}),
            plasma=Plasma(density=density, ion="Kr+"),
            cavity=Cavity(), boltzmann=BoltzmannSpec(), circuit=CircuitSpec(),
            volume_units=volume_units, reactions=KRF_TOML,
            diagnostics=[SpeciesDensities(("e-",), every=1)],
            output_dir=tempfile.gettempdir(),
        )

    def test_density_reported_in_volume_units(self) -> None:
        # Same physics (1e18 m^-3 electrons), reported in the chosen volume units.
        sim_m = self._units_sim("m", 1.0e18)           # 1e18 m^-3 in
        sim_m.run(Numerics(steps=1, dt=1e-12, fluid_every=1), verbose=False)
        self.assertAlmostEqual(sim_m.diagnostics[0]._series["e-"][-1], 1.0e18, delta=1.0e17)

        sim_cm = self._units_sim("cm", 1.0e12)         # 1e12 cm^-3 in (= 1e18 m^-3)
        sim_cm.run(Numerics(steps=1, dt=1e-12, fluid_every=1), verbose=False)
        self.assertAlmostEqual(sim_cm.diagnostics[0]._series["e-"][-1], 1.0e12, delta=1.0e11)
        # internal state is identical SI (1e18 m^-3) regardless of the display unit
        self.assertAlmostEqual(float(sim_cm.charged.electron_density().mean()), 1.0e18, delta=1.0e17)

    def test_composition_state(self) -> None:
        sim = _sim([SpeciesDensities(("Kr",))])
        sim.run(Numerics(steps=1, dt=1e-12, fluid_every=1), verbose=False)
        # quasi-neutral seed: electrons balance the singly-charged ions
        self.assertAlmostEqual(sim.charged.total_number("e-"),
                               sim.charged.total_number("Kr+"), delta=1.0)
        # gas mixture seeded at the 92/8 partial-density ratio
        ratio = sim.species_total("Kr") / sim.species_total("F2")
        self.assertAlmostEqual(ratio, 0.92 / 0.08, delta=0.5)

    def test_tanh_preionization_seeds_1d_quasineutral_profile(self) -> None:
        sim = Simulation(
            name="preion", grid=Grid1D(cells=32, length=12e-3, area_m2=90e-4),
            gas=GasMixture({"Kr": 0.99, "F2": 0.01}), plasma=Plasma(density=1e16, ion="Kr+"),
            cavity=Cavity(), boltzmann=BoltzmannSpec(), circuit=CircuitSpec(),
            preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0),
            reactions=KRF_TOML, diagnostics=[], output_dir=tempfile.gettempdir(),
        )
        num = Numerics(steps=1, dt=1e-12)
        sim.numerics = num
        sim._build(num)  # seed the fluids without stepping
        x = sim.mesh.cell_centers
        n_e = sim.charged.density_of("e-")
        np.testing.assert_allclose(n_e, 1e14 * np.tanh(100.0 * x) + 1e13, rtol=1e-9)
        np.testing.assert_allclose(sim.charged.density_of("Kr+"), n_e)  # quasi-neutral seed
        self.assertGreater(n_e[-1], 3.0 * n_e[0])  # rises anode->cathode (1D, not flat)

    def test_real_transport_evolves_charged_density(self) -> None:
        # With the MPC circuit, the real drift-diffusion transport model is selected and the
        # charged densities evolve from the quasi-neutral seed.
        with tempfile.TemporaryDirectory() as d:
            sim = Simulation(
                name="ddt", grid=Grid1D(cells=24, length=12e-3, area_m2=90e-4),
                gas=GasMixture({"Kr": 0.999, "F2": 0.001}, pressure_atm=3.4),
                plasma=Plasma(density=1e15, ion="Kr+", electron_mean_energy_ev=2.0),
                preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0),
                cavity=Cavity(), boltzmann=BoltzmannSpec(key_kind="mean_energy_eV"),
                circuit=CircuitSpec(backend="mpc", v0_v=-18e3),
                reactions=KRF_TOML, diagnostics=[ChargedDensityXT(("e-", "Kr+", "F-"), every=10)],
                output_dir=Path(d),
            )
            out = sim.run(Numerics(steps=60, dt=5e-11, fluid_every=10), verbose=False)
            from Fluid1D.plasma_model import DriftDiffusionPlasmaModel
            self.assertIsInstance(sim.plasma, DriftDiffusionPlasmaModel)
            # densities changed from the seed (electrons transported/attached)
            self.assertFalse(np.allclose(sim.charged.density, sim._charged_seed))
            # spatio-temporal archive present with (n_times, n_cells) blocks
            with h5py.File(out / "diagnostics.h5", "r") as f:
                g = f["charged_density_xt"]
                self.assertIn("times_ps", g)
                self.assertEqual(g["n_e-"].ndim, 2)
                self.assertEqual(g["n_e-"].shape[1], 24)  # n_cells
                self.assertIn(g.attrs["density_units"], ("m^-3", "cm^-3"))

    def test_diagnostics_are_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            sim = Simulation(
                name="bare", grid=Grid1D(cells=4, length=2.0e-2, area_m2=1.0e-4),
                gas=GasMixture({"Kr": 1.0}), plasma=Plasma(density=1e18, ion="Kr+"),
                cavity=Cavity(), boltzmann=BoltzmannSpec(), circuit=CircuitSpec(),
                reactions=KRF_TOML, diagnostics=[], output_dir=Path(d),
            )
            out = sim.run(Numerics(steps=4, dt=1e-12, fluid_every=2), verbose=False)
            self.assertEqual(list(out.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
