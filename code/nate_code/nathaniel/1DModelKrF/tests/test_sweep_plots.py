"""Phase 7: parallel F₂ sweep + overlay-plot module (tiny dev configs only)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    DataRef,
    GasMixture,
    Grid1D,
    LaserOutput,
    Numerics,
    PlasmaScalarTraces,
    Plasma,
    Simulation,
    SpeciesDensities,
)
from Fluid1D.sweep import sweep_f2_parallel
from Fluid1D.sweep_plots import (
    OverlayFigure,
    OverlayPanel,
    PanelOverlayFigure,
    plot_alpha_eta,
    plot_fwhm_ce,
    plot_sweep_overlays,
    plot_sweep_panels,
    pulse_fwhm_ns,
)

_ARF = ROOT / "data" / "reactions" / "arf.toml"
TEST_OVERLAYS = (
    OverlayFigure(
        filename="chosen_density_overlay.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/n_e_m3", "electron density (m⁻³)"),
        y_scale="log",
    ),
    OverlayFigure(
        filename="chosen_temperature_overlay.png",
        x=DataRef("plasma_scalars/times_ns", "time (ns)"),
        y=DataRef("plasma_scalars/te_ev", "electron temperature (eV)"),
    ),
)


def _tiny_base(tmp: Path) -> Simulation:
    mix = {"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001}
    return Simulation(
        name="tiny_sweep",
        grid=Grid1D(cells=16, length=12e-3, area_m2=90e-4),
        gas=GasMixture(mix, pressure_atm=3.4, temperature_k=300.0),
        plasma=Plasma(density=1e15, ion="Ne+", electron_mean_energy_ev=3.0),
        cavity=Cavity(wavelength_nm=193.0, r_out=0.08),
        boltzmann=BoltzmannSpec(backend="stub"),
        circuit=CircuitSpec(backend="stub", phi_left=0.0, phi_right=-4.0e3),
        reactions=_ARF,
        output_dir=tmp,
        diagnostics=[],
    )


TEST_PANELS = PanelOverlayFigure(
    filename="chosen_species_panels.png",
    panels=tuple(
        OverlayPanel(
            tag=tag,
            x=DataRef("species_densities/times_ps", "t (ns)", scale=1.0e-3),
            y=DataRef(f"species_densities/{name}", f"n_{name} (m⁻³)"),
            x_include_zero=True, y_include_zero=True,
        )
        for tag, name in (("(a)", "Ar+"), ("(b)", "Ar*"), ("(c)", "ArF*"), ("(d)", "F-"))
    ),
)


def _point_diags(f2: float) -> list:
    return [PlasmaScalarTraces(), LaserOutput(),
            SpeciesDensities(track=("Ar+", "Ar*", "ArF*", "F-"))]


class TestFitYRange(unittest.TestCase):
    def test_nonnegative_data_fits_120_percent(self) -> None:
        from Fluid1D.sweep_plots import fit_y_range

        lo, hi = fit_y_range([np.array([0.0, 2.0, 5.0]), np.array([1.0])])
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 6.0)

    def test_signed_data_fits_asymmetric_extent(self) -> None:
        from Fluid1D.sweep_plots import fit_y_range

        lo, hi = fit_y_range([np.array([-4.0, 3.0])])
        self.assertAlmostEqual(lo, -5.4)
        self.assertAlmostEqual(hi, 4.4)

    def test_empty_and_flat_zero(self) -> None:
        from Fluid1D.sweep_plots import fit_y_range

        self.assertEqual(fit_y_range([]), (0.0, 1.0))
        self.assertEqual(fit_y_range([np.zeros(4)]), (0.0, 1.0))


class TestPulseFwhm(unittest.TestCase):
    def test_triangle_pulse(self) -> None:
        t = np.linspace(0.0, 20.0, 401)
        p = np.maximum(0.0, 1.0 - np.abs(t - 10.0) / 5.0)  # FWHM = 5 ns
        self.assertAlmostEqual(pulse_fwhm_ns(t, p), 5.0, delta=0.2)

    def test_flat_zero_trace(self) -> None:
        self.assertEqual(pulse_fwhm_ns(np.arange(10.0), np.zeros(10)), 0.0)


class TestSweepEndToEnd(unittest.TestCase):
    """Tiny stub-backed sweep (16 cells, 30 steps): engine + h5 + every figure renders."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        base = _tiny_base(tmp)
        cls.result = sweep_f2_parallel(
            base, values=(0.001, 0.002), numerics=Numerics(steps=30, dt=1e-12, fluid_every=10),
            rebalance="Ne", max_workers=2, diagnostics_factory=_point_diags,
        )
        cls.sweep_dir = tmp / "tiny_sweep"
        cls.base = base

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_points_and_archives(self) -> None:
        self.assertEqual(len(self.result.points), 2)
        for pt in self.result.points:
            self.assertTrue((pt.out_dir / "diagnostics.h5").exists())

    def test_plasma_scalars_series(self) -> None:
        import h5py

        with h5py.File(self.result.points[0].out_dir / "diagnostics.h5", "r") as h5:
            grp = h5["plasma_scalars"]
            for key in ("times_ns", "n_e_m3", "en_td", "te_ev", "p_dis_wm3", "p_l_wm3"):
                self.assertIn(key, grp)
            self.assertGreater(grp["n_e_m3"][:].max(), 0.0)
            for attr in ("e_out_j", "e_dep_j", "e_bank0_j", "e_stored0_j"):
                self.assertIn(attr, grp.attrs)
            self.assertAlmostEqual(float(grp.attrs["e_bank0_j"]), 1.62, places=12)

    def test_overlay_figures_render(self) -> None:
        paths = plot_sweep_overlays(self.sweep_dir, TEST_OVERLAYS)
        self.assertEqual(len(paths), len(TEST_OVERLAYS))
        for p in paths:
            self.assertTrue(p.exists(), p)
        self.assertEqual(
            {path.name for path in paths},
            {"chosen_density_overlay.png", "chosen_temperature_overlay.png"},
        )

    def test_fwhm_ce_figure(self) -> None:
        p = plot_fwhm_ce(self.sweep_dir, filename="chosen_fwhm_ce.png")
        self.assertTrue(p.exists())

    def test_species_panel_figure_renders(self) -> None:
        p = plot_sweep_panels(self.sweep_dir, TEST_PANELS)
        self.assertTrue(p.exists())
        self.assertEqual(p.name, "chosen_species_panels.png")

    def test_alpha_eta_figure(self) -> None:
        p = plot_alpha_eta(self.base, (0.0005, 0.005),
                           out_path=self.sweep_dir / "alpha_eta_en.png")
        self.assertTrue(p.exists())


if __name__ == "__main__":
    unittest.main()
