"""Configurable terminal progress reporting (Fluid1D.progress / ProgressReport)."""

from __future__ import annotations

import contextlib
import io
import json
import re
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    ProgressReport,
    Simulation,
    SpeciesDensities,
)
from Fluid1D.progress import ProgressReporter, format_duration

STEPS = 20


def _simulation(tmpdir: str) -> Simulation:
    return Simulation(
        name="progress_smoke",
        grid=Grid1D(cells=4, length=2.0e-2, area_m2=1.0e-4),
        gas=GasMixture({"Kr": 0.92, "F2": 0.08}, pressure_atm=1.0, temperature_k=300.0),
        plasma=Plasma(density=1.0e18, ion="Kr+"),
        cavity=Cavity(wavelength_nm=248.0),
        boltzmann=BoltzmannSpec(backend="stub"),
        circuit=CircuitSpec(),
        reactions=_ROOT / "data" / "reactions" / "krf.toml",
        output_dir=tmpdir,
        diagnostics=[SpeciesDensities(("e-",))],
    )


def _run(
    tmpdir: str, progress: ProgressReport | None, *, verbose: bool = True
) -> tuple[list[str], list[str]]:
    """Run a tiny simulation; return (progress lines, other non-empty lines) from stdout."""
    sim = _simulation(tmpdir)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sim.run(Numerics(steps=STEPS, dt=1.0e-12, fluid_every=10),
                progress=progress, verbose=verbose)
    lines = [ln for ln in buf.getvalue().splitlines() if ln.strip()]
    is_progress = [bool(re.match(r"\s+(\d+%|\d+/\d+)\s", ln)) for ln in lines]
    return ([ln for ln, hit in zip(lines, is_progress) if hit],
            [ln for ln, hit in zip(lines, is_progress) if not hit])


class TestProgressLines(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        cls.tmpdir = cls._tmp.name

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_default_format(self) -> None:
        lines, _ = _run(self.tmpdir, progress=None)
        self.assertEqual(len(lines), 10)  # every 10% of 20 steps -> every 2 steps
        # The default ``excimer`` field resolves to this network's lasing upper state (KrF*).
        pattern = (r"^\s+\d+%\s\s"
                   r"KrF\*=\d\.\d{2}e[+-]\d{2}\s\s"
                   r"photons=\d\.\d{2}e[+-]\d{2}\s\s"
                   r"T\(e-\)=[\d.eE+-]+ eV\s\s"
                   r"η=\d\.\d{2}e[+-]\d{2}\s\s"
                   r"elapsed=\S+\s\s"
                   r"eta=\S+$")
        for ln in lines:
            self.assertRegex(ln, pattern)

    def test_timing_disabled_uses_untimed_format(self) -> None:
        lines, _ = _run(self.tmpdir, ProgressReport(timing=False))
        self.assertEqual(len(lines), 10)
        untimed = (r"^\s+\d+%\s\s"
                   r"KrF\*=\d\.\d{2}e[+-]\d{2}\s\s"
                   r"photons=\d\.\d{2}e[+-]\d{2}\s\s"
                   r"T\(e-\)=[\d.eE+-]+ eV\s\s"
                   r"η=\d\.\d{2}e[+-]\d{2}$")
        for ln in lines:
            self.assertRegex(ln, untimed)

    def test_steps_style_and_every(self) -> None:
        lines, _ = _run(self.tmpdir, ProgressReport(every=5, style="steps"))
        self.assertEqual(len(lines), 4)
        self.assertTrue(lines[0].lstrip().startswith("5/20"))
        self.assertIn("20/20", lines[-1])

    def test_final_step_reported_with_uneven_cadence(self) -> None:
        # 20 % 3 != 0: the run must still close on the 20/20 line (7 lines, not 6).
        lines, _ = _run(self.tmpdir, ProgressReport(every=3, style="steps"))
        self.assertEqual(len(lines), 7)
        self.assertIn("20/20", lines[-1])

    def test_explicit_progress_reports_when_quiet(self) -> None:
        # An explicit spec reports even under verbose=False; banner/summary stay silent.
        lines, other = _run(self.tmpdir, ProgressReport(every=10), verbose=False)
        self.assertEqual(len(lines), 2)
        self.assertEqual(other, [])

    def test_quiet_run_without_spec_prints_nothing(self) -> None:
        lines, other = _run(self.tmpdir, progress=None, verbose=False)
        self.assertEqual(lines, [])
        self.assertEqual(other, [])

    def test_custom_fields_render(self) -> None:
        fields = ("total:e-", "density:e-", "temp:e-", "photons",
                  "p_in", "p_out", "efficiency", "gain")
        lines, _ = _run(self.tmpdir, ProgressReport(every=10, fields=fields))
        self.assertEqual(len(lines), 2)
        last = lines[-1]
        for token in ("e-=", "n(e-)=", "T(e-)=", "photons=", "P_in=", "P_out=",
                      "η=", "gain="):
            self.assertIn(token, last)
        m = re.search(r"n\(e-\)=(\S+)", last)
        self.assertIsNotNone(m)
        self.assertTrue(np.isfinite(float(m.group(1))))

    def test_invalid_field_fails_before_stepping(self) -> None:
        sim = _simulation(self.tmpdir)
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                sim.run(Numerics(steps=STEPS, dt=1.0e-12),
                        progress=ProgressReport(fields=("bogus_field",)))

    def test_unknown_species_fails_before_stepping(self) -> None:
        sim = _simulation(self.tmpdir)
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                sim.run(Numerics(steps=STEPS, dt=1.0e-12),
                        progress=ProgressReport(fields=("density:Unobtanium",)))

    def test_bad_every_and_style_rejected(self) -> None:
        sim = _simulation(self.tmpdir)
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                sim.run(Numerics(steps=STEPS, dt=1.0e-12),
                        progress=ProgressReport(every=0))
            with self.assertRaises(ValueError):
                sim.run(Numerics(steps=STEPS, dt=1.0e-12),
                        progress=ProgressReport(style="fraction"))

    def test_terminal_and_heartbeat_both_off_rejected(self) -> None:
        sim = _simulation(self.tmpdir)
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                sim.run(Numerics(steps=STEPS, dt=1.0e-12),
                        progress=ProgressReport(terminal=False))

    def test_heartbeat_single_run(self) -> None:
        sim = _simulation(self.tmpdir)
        before = time.time()
        with contextlib.redirect_stdout(io.StringIO()):
            out = sim.run(Numerics(steps=STEPS, dt=1.0e-12, fluid_every=10),
                          progress=ProgressReport(every=5, heartbeat=True))
        after = time.time()
        record = json.loads((out / "progress.json").read_text(encoding="utf-8"))
        self.assertEqual(record["schema"], 1)
        self.assertEqual(record["status"], "done")
        self.assertEqual(record["step"], STEPS)
        self.assertEqual(record["steps"], STEPS)
        self.assertEqual(record["pct"], 100.0)
        self.assertAlmostEqual(record["t_sim_s"], STEPS * 1.0e-12)
        self.assertGreaterEqual(record["wall_s"], 0.0)
        self.assertTrue(before <= record["updated_unix"] <= after)
        self.assertEqual(set(record["fields"]), {"KrF*", "photons", "T(e-)", "η"})
        self.assertTrue(all(isinstance(v, str) for v in record["fields"].values()))
        # Atomic-write hygiene: no tmp residue beside the heartbeat.
        self.assertFalse((out / "progress.json.tmp").exists())

    def test_heartbeat_without_terminal_is_silent(self) -> None:
        sim = _simulation(self.tmpdir)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            out = sim.run(Numerics(steps=STEPS, dt=1.0e-12, fluid_every=10),
                          progress=ProgressReport(every=5, terminal=False, heartbeat=True),
                          verbose=False)
        self.assertEqual([ln for ln in buf.getvalue().splitlines() if ln.strip()], [])
        record = json.loads((out / "progress.json").read_text(encoding="utf-8"))
        self.assertEqual(record["status"], "done")


class TestFormatDuration(unittest.TestCase):
    def test_ranges(self) -> None:
        self.assertEqual(format_duration(0.0), "0.0s")
        self.assertEqual(format_duration(12.34), "12.3s")
        self.assertEqual(format_duration(75.0), "1m15s")
        self.assertEqual(format_duration(3725.0), "1h02m")


class TestReporterValues(unittest.TestCase):
    """Numeric spot-checks against directly computed values on a built sim."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        cls.sim = _simulation(cls._tmp.name)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.sim.run(Numerics(steps=2, dt=1.0e-12, fluid_every=10),
                        progress=ProgressReport(every=1))

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_density_matches_species_total(self) -> None:
        reporter = ProgressReporter(
            ProgressReport(fields=("density:e-",), timing=False), self.sim
        )
        rendered = reporter.line(2).split("  ")[-1]
        value = float(rendered.split("=")[1])
        vol_m3 = self.sim.domain_length * self.sim.area_si
        expected = self.sim.units.density_from_si(self.sim.species_total("e-") / vol_m3)
        self.assertAlmostEqual(value, expected, delta=0.01 * expected)

    def test_temp_matches_electron_temperature(self) -> None:
        reporter = ProgressReporter(
            ProgressReport(fields=("temp:e-",), timing=False), self.sim
        )
        rendered = reporter.line(2).split("  ")[-1]
        value = float(rendered.split("=")[1].split()[0])
        self.assertAlmostEqual(value, self.sim.electron_temperature(), delta=0.05)

    def test_excimer_field_renders_lasing_species(self) -> None:
        reporter = ProgressReporter(
            ProgressReport(fields=("excimer",), timing=False), self.sim
        )
        rendered = reporter.line(2).split("  ")[-1]
        self.assertTrue(rendered.startswith("KrF*="))
        expected = self.sim.species_total("KrF*")
        self.assertAlmostEqual(
            float(rendered.split("=")[1]), expected, delta=0.01 * expected + 1e-30
        )
        # A network with no STIMULATED channel reports n/a instead of failing.
        original = self.sim.gain.lasing
        try:
            self.sim.gain.lasing = None
            reporter = ProgressReporter(
                ProgressReport(fields=("excimer",), timing=False), self.sim
            )
            self.assertEqual(reporter.line(2).split("  ")[-1], "excimer=n/a")
        finally:
            self.sim.gain.lasing = original


if __name__ == "__main__":
    unittest.main()
