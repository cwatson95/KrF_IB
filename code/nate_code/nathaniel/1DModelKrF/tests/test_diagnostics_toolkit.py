"""Diagnostics toolkit: headless plotting, archive, tallies, efficiency."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless

import numpy as np

import h5py

from diagnostics import (
    DataRef,
    DiagnosticsArchive,
    HeatmapPanel,
    HeatmapPlot,
    LinePlot,
    data_limits,
    export_figure,
    plot_time_series,
)
from Fluid1D.diagnostics import AVOGADRO, ExcimerRates
from kinetics.efficiency import LaserEnergetics


class _Namespace:
    def __init__(self, **kw) -> None:
        self.__dict__.update(kw)


class TestExcimerRatesSeries(unittest.TestCase):
    """Formation/quenching split by ν sign + the optics-owned SE/STE quench series."""

    def _fake_sim(self, with_gain: bool) -> _Namespace:
        chem = _Namespace(
            species_names=["ArF*"],
            nu=np.array([[+1.0, -1.0]]),  # column 0 forms ArF*, column 1 quenches it
        )
        gain = _Namespace(
            last_dt_s=1.0e-9,
            last_se_sink_m3={"ArF*": 5.0 * AVOGADRO * 1.0e-9},   # -> 5 mol m^-3 s^-1
            last_ste_sink_m3={"ArF*": 2.0 * AVOGADRO * 1.0e-9},  # -> 2 mol m^-3 s^-1
        ) if with_gain else None
        sim = _Namespace(
            numerics=_Namespace(fluid_every=1),
            plasma_chem=chem,
            reaction_rates_mean=lambda: (
                ("R_f", "R_q"), np.array([3.0, 4.0]) * AVOGADRO,
            ),
        )
        if gain is not None:
            sim.gain = gain
        return sim

    def test_split_includes_se_ste_from_the_cavity_ledger(self) -> None:
        diag = ExcimerRates(excimer="ArF*", every=1)
        diag.record(self._fake_sim(with_gain=True), 1, 1.0e-9)
        form, quench = diag._split()
        self.assertAlmostEqual(float(form["R_f"][0]), 3.0, places=9)
        self.assertAlmostEqual(float(quench["R_q"][0]), 4.0, places=9)
        self.assertAlmostEqual(float(quench["SE"][0]), 5.0, places=9)
        self.assertAlmostEqual(float(quench["STE"][0]), 2.0, places=9)

    def test_no_cavity_records_zero_radiative_sinks(self) -> None:
        diag = ExcimerRates(excimer="ArF*", every=1)
        diag.record(self._fake_sim(with_gain=False), 1, 1.0e-9)
        _form, quench = diag._split()
        self.assertEqual(float(quench["SE"][0]), 0.0)
        self.assertEqual(float(quench["STE"][0]), 0.0)


class TestPlots(unittest.TestCase):
    """Every plot renders + writes a non-empty PNG headlessly."""

    def _assert_png(self, out: Path) -> None:
        self.assertTrue(out.is_file() and out.stat().st_size > 0)

    def test_time_series(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            t = np.linspace(0, 50, 40)
            plot_time_series(t, {"n_e": np.exp(t / 20), "KrF*": np.sin(t) ** 2 + 0.1}, log_y=True)
            self._assert_png(export_figure(Path(d) / "series.png"))

    def test_time_series_distinguishes_coincident_lines(self) -> None:
        # KrF* and F track 1:1 and have identical values; distinct styles/markers must
        # keep both visible (one plain line would sit exactly on top of the other).
        import matplotlib.pyplot as plt

        t = np.linspace(0, 10, 11)
        same = t**2
        ax = plot_time_series(t, {"KrF*": same, "F": same.copy()})
        lines = ax.get_lines()
        self.assertNotEqual(
            (lines[0].get_linestyle(), lines[0].get_marker()),
            (lines[1].get_linestyle(), lines[1].get_marker()),
        )
        plt.close(ax.figure)

    def test_data_limits_fit_without_decade_whitespace(self) -> None:
        lo, hi = data_limits(np.array([0.0, 1.5e-3]))
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 1.575e-3)
        lo, hi = data_limits(np.array([-3.0, 2.0]))
        self.assertAlmostEqual(lo, -3.25)
        self.assertAlmostEqual(hi, 2.25)

    def test_data_limits_ignore_nonfinite_and_pad_logarithmically(self) -> None:
        lo, hi = data_limits(np.array([np.nan, 1.0e2, 1.0e4, np.inf]), scale="log")
        self.assertLess(lo, 1.0e2)
        self.assertGreater(hi, 1.0e4)
        self.assertAlmostEqual(np.log10(hi / lo), 2.2)

    def test_time_series_linear_ylim(self) -> None:
        import matplotlib.pyplot as plt

        t = np.linspace(0, 10, 5)
        ax = plot_time_series(t, {"a": np.array([0.0, 2.0e-3, 1.0e-3, 0.0, 0.0])},
                              ylim=(0.0, 0.01))
        self.assertEqual(ax.get_ylim(), (0.0, 0.01))
        plt.close(ax.figure)

    def test_archive_backed_line_plot_uses_run_owned_name_and_visible_data(self) -> None:
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            archive_path = root / "diagnostics.h5"
            with DiagnosticsArchive(archive_path) as archive:
                archive.series(
                    "trace", np.array([0.0, 1.0, 2.0, 3.0]),
                    {"value": np.array([0.0, 2.0, 4.0, 1.0e6])},
                )
            plot = LinePlot(
                filename="chosen_by_test.py.png",
                x=DataRef("trace/times_ps", "time (ps)"),
                y=(DataRef("trace/value", "value"),),
                ylabel="value",
                x_limits=(0.0, 2.0),
            )
            path = plot.render(archive_path, root)[0]
            self._assert_png(path)
            self.assertEqual(path.name, "chosen_by_test.py.png")
            self.assertFalse((root / "trace.png").exists())
            self.assertEqual(plt.get_fignums(), [])

    def test_archive_backed_heatmap_renders(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            archive_path = root / "diagnostics.h5"
            with DiagnosticsArchive(archive_path) as archive:
                group = archive.group("xt")
                group.create_dataset("x", data=np.array([0.0, 1.0, 2.0]))
                group.create_dataset("t", data=np.array([0.0, 1.0]))
                group.create_dataset("n", data=np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 8.0]]))
            plot = HeatmapPlot(
                filename="chosen_heatmap.png",
                x=DataRef("xt/x", "x (m)"),
                y=DataRef("xt/t", "time (s)"),
                panels=(HeatmapPanel(DataRef("xt/n", "n (m⁻³)"), "density"),),
                color_scale="log",
            )
            self._assert_png(plot.render(archive_path, root)[0])



class TestDiagnosticsArchive(unittest.TestCase):
    def test_series_and_array_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "diagnostics.h5"
            t = np.linspace(0, 10, 11)
            with DiagnosticsArchive(path) as arc:
                arc.series("optical_gain", t, {"gain_per_m": t**2})
                arc.array("eedf", "bin_edges_ev", np.logspace(-1, 2, 31))
            with h5py.File(path, "r") as f:
                np.testing.assert_allclose(f["optical_gain/times_ps"][:], t)
                np.testing.assert_allclose(f["optical_gain/gain_per_m"][:], t**2)
                self.assertEqual(f["eedf/bin_edges_ev"].shape, (31,))



class TestEfficiency(unittest.TestCase):
    def test_interval_average_integrates_first_and_unequal_intervals(self) -> None:
        e = LaserEnergetics()
        e.record(2.0, 0.0, 3.0, 5.0, interval_s=2.0)
        self.assertEqual(e.output_energy_j(), 6.0)
        self.assertEqual(e.input_energy_j(), 10.0)
        e.record(2.5, 0.0, 7.0, 11.0, interval_s=0.5)
        self.assertEqual(e.output_energy_j(), 9.5)
        self.assertEqual(e.input_energy_j(), 15.5)

    def test_energies_and_efficiencies(self) -> None:
        e = LaserEnergetics()
        # constant 1.5 W spont + 0.5 W stim out, 10 W in over 1e-6 s
        # -> 2e-6 J out (1.5e-6 spont / 5e-7 stim), 1e-5 J in
        for k in range(11):
            e.record(k * 1e-7, 1.5, 0.5, 10.0)
        self.assertAlmostEqual(e.output_energy_j(), 2.0e-6, places=12)
        self.assertAlmostEqual(e.output_energy_spont_j(), 1.5e-6, places=12)
        self.assertAlmostEqual(e.output_energy_stim_j(), 5.0e-7, places=12)
        self.assertAlmostEqual(e.input_energy_j(), 1.0e-5, places=12)
        # no pre-pulse stored energy -> denominator is the drive work alone
        self.assertAlmostEqual(e.intrinsic_efficiency(), 0.2, places=9)
        self.assertAlmostEqual(e.wall_plug_efficiency(4.0e-6), 0.5, places=9)

    def test_initial_energy_enters_the_denominator(self) -> None:
        # a pre-ionized plasma's stored energy is deposited energy too:
        # eta = out / (E_init + integral P_in)
        e = LaserEnergetics(initial_energy_j=1.0e-5)
        for k in range(11):
            e.record(k * 1e-7, 2.0, 0.0, 10.0)
        self.assertAlmostEqual(e.deposited_energy_j(), 2.0e-5, places=12)
        self.assertAlmostEqual(e.intrinsic_efficiency(), 0.1, places=9)

    def test_input_power_is_signed(self) -> None:
        # intervals where the particles do net work on the drive integrate negative;
        # the ledger must not rectify them away
        e = LaserEnergetics(initial_energy_j=4.0e-6)
        for k in range(11):
            e.record(k * 1e-7, 1.0, 0.0, -2.0)  # particles push back on the drive
        self.assertAlmostEqual(e.input_energy_j(), -2.0e-6, places=12)
        self.assertAlmostEqual(e.deposited_energy_j(), 2.0e-6, places=12)
        self.assertAlmostEqual(e.intrinsic_efficiency(), 0.5, places=9)

    def test_cumulative_efficiency_stays_finite_with_returned_drive_work(self) -> None:
        from Fluid1D.diagnostics import Efficiency

        e = LaserEnergetics(initial_energy_j=1.0e-9)
        e.record(1e-9, 0.0, 1.0e-6, 1.0e3)
        e.record(2e-9, 0.0, 1.0e-6, -2.0e4)  # net return drives signed deposit through ~0
        cum = Efficiency._cumulative(e)
        self.assertTrue(np.all(np.isfinite(cum)))
        self.assertLess(cum.max(), 1.0)

        # a non-positive total deposit yields 0, not a nonsense ratio
        drained = LaserEnergetics(initial_energy_j=1.0e-6)
        for k in range(11):
            drained.record(k * 1e-7, 1.0, 0.0, -2.0)
        self.assertEqual(drained.intrinsic_efficiency(), 0.0)


if __name__ == "__main__":
    unittest.main()
