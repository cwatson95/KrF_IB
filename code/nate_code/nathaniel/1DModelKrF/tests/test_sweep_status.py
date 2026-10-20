"""sweep_status heartbeat table: scanning, staleness, rendering, CLI return codes."""

from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from Fluid1D.sweep_status import main, render_table, scan_sweep

_NOW = 1_754_000_000.0


def _record(**overrides) -> dict:
    record = {
        "schema": 1, "status": "running", "step": 12, "steps": 30, "pct": 40.0,
        "t_sim_s": 1.2e-11, "wall_s": 3.2, "updated_unix": _NOW,
        "fields": {"photons": "1.00e+05"},
    }
    record.update(overrides)
    return record


def _build_sweep(tmp: Path) -> Path:
    sweep = tmp / "sweep"
    cases = {
        "f2_0.001": _record(status="done", step=30, pct=100.0),
        "f2_0.002": _record(),                                    # running, fresh
        "f2_0.003": _record(updated_unix=_NOW - 1000.0),          # running, stale
        "f2_0.004": _record(status="failed",
                            error="Traceback ...\nValueError: boom"),
    }
    for name, record in cases.items():
        (sweep / name).mkdir(parents=True)
        (sweep / name / "progress.json").write_text(json.dumps(record), encoding="utf-8")
    (sweep / "f2_0.005").mkdir()                                  # no heartbeat yet
    (sweep / "f2_readme").mkdir()                                 # non-numeric: skipped
    (sweep / "f2_0.002" / "progress.json.tmp").write_text("{", encoding="utf-8")
    return sweep


class TestScanAndRender(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.sweep = _build_sweep(Path(self._tmp.name))

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_scan_orders_and_skips_non_numeric(self) -> None:
        statuses = scan_sweep(self.sweep)
        self.assertEqual([status.f2 for status in statuses],
                         [0.001, 0.002, 0.003, 0.004, 0.005])
        self.assertIsNone(statuses[-1].record)

    def test_render_statuses_and_footer(self) -> None:
        table = render_table(scan_sweep(self.sweep), now_unix=_NOW)
        rows = table.splitlines()
        self.assertIn("done", rows[1])
        self.assertIn("running", rows[2])
        self.assertIn("40%", rows[2])
        self.assertIn("stalled", rows[3])
        self.assertIn("failed", rows[4])
        self.assertIn("ValueError: boom", rows[4])
        self.assertIn("0.005", rows[5])
        footer = rows[-1]
        self.assertIn("[5 points]", footer)
        for token in ("done 1", "running 1", "stalled 1", "failed 1", "no-data 1"):
            self.assertIn(token, footer)

    def test_stale_threshold_honored(self) -> None:
        table = render_table(scan_sweep(self.sweep), stale_after_s=1e9, now_unix=_NOW)
        self.assertNotIn("stalled", table)

    def test_main_return_codes(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main([str(self.sweep)]), 0)
            self.assertEqual(main([str(self.sweep / "missing")]), 2)
            with tempfile.TemporaryDirectory() as empty:
                self.assertEqual(main([empty]), 1)


if __name__ == "__main__":
    unittest.main()
