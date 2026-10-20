"""Parallel-sweep live reporting: parent lines, failure isolation, worker heartbeats."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import tempfile
import time
import unittest
from pathlib import Path

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
)
from Fluid1D.progress import write_heartbeat
from Fluid1D.sweep import (
    _point_dir_name,
    _rollup_line,
    _running_points,
    _worker_progress_spec,
    sweep_f2_parallel,
)

_ARF = _ROOT / "data" / "reactions" / "arf.toml"
_NUMERICS = Numerics(steps=30, dt=1e-12, fluid_every=10)


def _tiny_base(tmp: Path, name: str = "progress_sweep") -> Simulation:
    mix = {"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001}
    return Simulation(
        name=name,
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


def _exit_factory(f2: float) -> list:  # pragma: no cover - runs in a child that dies
    """Diagnostics factory that kills the worker outright (BrokenProcessPool probe)."""
    os._exit(13)


def _fail_second_factory(f2: float) -> list:  # pragma: no cover - runs in a worker
    """Raise for one valid point while allowing the other point to finish."""
    if f2 == 0.002:
        raise RuntimeError("intentional per-point failure")
    return []


class TestSweepFailureIsolation(unittest.TestCase):
    """One real 2-worker sweep with one success and one factory failure."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cls.result = sweep_f2_parallel(
                _tiny_base(tmp), values=(0.001, 0.002), numerics=_NUMERICS,
                rebalance="Ne", max_workers=2, diagnostics_factory=_fail_second_factory,
            )
        cls.stdout = buf.getvalue()
        cls.sweep_dir = tmp / "progress_sweep"

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_success_recorded_with_wall_time(self) -> None:
        self.assertEqual(len(self.result.points), 1)
        point = self.result.points[0]
        self.assertEqual(point.f2_fraction, 0.001)
        self.assertGreater(point.wall_s, 0.0)
        self.assertEqual(self.result.f2().tolist(), [0.001])
        self.assertIs(self.result.best(), point)

    def test_failure_recorded_with_traceback(self) -> None:
        self.assertEqual(len(self.result.failures), 1)
        failure = self.result.failures[0]
        self.assertEqual(failure.f2_fraction, 0.002)
        self.assertIn("RuntimeError", failure.error)
        self.assertIn("Traceback", failure.traceback)
        self.assertEqual(failure.out_dir.name, "f2_0.002")

    def test_success_heartbeat_done(self) -> None:
        record = json.loads(
            (self.sweep_dir / "f2_0.001" / "progress.json").read_text(encoding="utf-8"))
        self.assertEqual(record["schema"], 1)
        self.assertEqual(record["status"], "done")
        self.assertEqual(record["step"], 30)
        self.assertEqual(record["pct"], 100.0)

    def test_failed_heartbeat_records_error(self) -> None:
        record = json.loads(
            (self.sweep_dir / "f2_0.002" / "progress.json").read_text(encoding="utf-8"))
        self.assertEqual(record["status"], "failed")
        self.assertIn("RuntimeError", record["error"])

    def test_parent_lines(self) -> None:
        self.assertRegex(self.stdout,
                         r"\[sweep\] progress_sweep: 2 points, .*2 workers, 30 steps")
        completion = [ln for ln in self.stdout.splitlines()
                      if re.match(r"\[sweep\] [12]/2\s", ln)]
        self.assertEqual(len(completion), 2)
        self.assertEqual(sum("FAILED" in ln for ln in completion), 1)
        self.assertIn("1 of 2 points FAILED", self.stdout)
        self.assertRegex(self.stdout, r"\[sweep\] complete: 1 ok, 1 failed")

    def test_parent_owns_the_terminal(self) -> None:
        # Workers are heartbeat-only by construction; every captured line is a parent line.
        for ln in self.stdout.splitlines():
            if ln.strip():
                self.assertTrue(ln.startswith("[sweep]"), ln)


class TestSweepProgressModes(unittest.TestCase):
    def test_progress_false_is_fully_silent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = _tiny_base(Path(tmp), name="silent_sweep")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                result = sweep_f2_parallel(
                    base, values=(0.001,), numerics=_NUMERICS,
                    rebalance="Ne", max_workers=1, progress=False,
                )
            self.assertEqual(buf.getvalue(), "")
            self.assertEqual(len(result.points), 1)
            self.assertFalse(
                (Path(tmp) / "silent_sweep" / "f2_0.001" / "progress.json").exists())

    def test_custom_spec_forced_to_heartbeat_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = _tiny_base(Path(tmp), name="custom_sweep")
            with contextlib.redirect_stdout(io.StringIO()):
                result = sweep_f2_parallel(
                    base, values=(0.001,), numerics=_NUMERICS, rebalance="Ne",
                    max_workers=1,
                    progress=ProgressReport(every=10, fields=("photons",)),
                )
            self.assertEqual(len(result.points), 1)
            record = json.loads((Path(tmp) / "custom_sweep" / "f2_0.001" /
                                 "progress.json").read_text(encoding="utf-8"))
            self.assertEqual(record["status"], "done")
            self.assertEqual(list(record["fields"]), ["photons"])

    def test_broken_pool_recorded_as_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = _tiny_base(Path(tmp), name="broken_sweep")
            with contextlib.redirect_stdout(io.StringIO()):
                result = sweep_f2_parallel(
                    base, values=(0.001,), numerics=_NUMERICS, rebalance="Ne",
                    max_workers=1, diagnostics_factory=_exit_factory,
                )
            self.assertEqual(result.points, [])
            self.assertEqual(len(result.failures), 1)
            self.assertIn("BrokenProcessPool", result.failures[0].error)


class TestProgressHelpers(unittest.TestCase):
    """Pure-function units: no process pool, no timing dependence."""

    def test_worker_progress_spec_normalization(self) -> None:
        self.assertIsNone(_worker_progress_spec(False))
        default = _worker_progress_spec(True)
        self.assertFalse(default.terminal)
        self.assertTrue(default.heartbeat)
        custom = _worker_progress_spec(ProgressReport(every=7, fields=("photons",)))
        self.assertEqual(custom.every, 7)
        self.assertEqual(custom.fields, ("photons",))
        self.assertFalse(custom.terminal)
        self.assertTrue(custom.heartbeat)
        with self.assertRaises(TypeError):
            _worker_progress_spec("yes")

    def test_point_dir_name_matches_run_layout(self) -> None:
        self.assertEqual(_point_dir_name(0.0015), "f2_0.0015")
        self.assertEqual(_point_dir_name(2.0), "f2_2")

    def test_running_points_filters_stale_heartbeats(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sweep_dir = Path(tmp)
            now = time.time()
            fresh = sweep_dir / _point_dir_name(0.001)
            stale = sweep_dir / _point_dir_name(0.002)
            fresh.mkdir()
            stale.mkdir()
            write_heartbeat(fresh / "progress.json",
                            {"schema": 1, "status": "running", "pct": 40.0,
                             "updated_unix": now})
            write_heartbeat(stale / "progress.json",
                            {"schema": 1, "status": "running", "pct": 80.0,
                             "updated_unix": now - 3600.0})
            pcts = _running_points(sweep_dir, [0.001, 0.002, 0.003], t0_wall=now - 5.0)
            self.assertEqual(pcts, [40.0])

    def test_rollup_line_text(self) -> None:
        line = _rollup_line(3, 10, [40.0, 44.0], 5, elapsed_s=60.0)
        self.assertIn("done 3/10", line)
        self.assertIn("running 2 (avg 42%)", line)
        self.assertIn("queued 5", line)
        self.assertIn("eta~", line)
        self.assertIn("eta n/a", _rollup_line(0, 10, [], 10, elapsed_s=5.0))


if __name__ == "__main__":
    unittest.main()
