"""CLI selection regressions for the ArF and KrF sweep programs."""

from __future__ import annotations

import importlib.util
import io
import sys
import unittest
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def _load_example(stem: str) -> ModuleType:
    path = REPOSITORY_ROOT / "examples" / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"test_{stem}_entrypoint", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestSweepEntrypoints(unittest.TestCase):
    def test_tiny_presets_pass_only_explicit_values(self) -> None:
        for stem in ("arf_f2_sweep", "krf_f2_sweep"):
            with self.subTest(stem=stem):
                module = _load_example(stem)
                call: dict[str, object] = {}

                def fake_sweep(base: object, **kwargs: object) -> SimpleNamespace:
                    call.update(kwargs)
                    return SimpleNamespace(points=[])

                argv = [str(module.__file__), "--tiny", "--name", f"test_{stem}"]
                with (
                    patch.object(module, "Simulation", return_value=object()),
                    patch.object(module, "sweep_f2_parallel", side_effect=fake_sweep),
                    patch.object(module, "render_standard_sweep_figures", return_value=()),
                    patch.object(sys, "argv", argv),
                    patch("sys.stdout", new=io.StringIO()),
                ):
                    module.main()

                self.assertEqual(call["values"], (0.0005, 0.001, 0.002))
                self.assertNotIn("start", call)
                self.assertNotIn("stop", call)
                self.assertNotIn("step", call)

    def test_production_preset_passes_only_range(self) -> None:
        module = _load_example("arf_f2_sweep")
        call: dict[str, object] = {}

        def fake_sweep(base: object, **kwargs: object) -> SimpleNamespace:
            call.update(kwargs)
            return SimpleNamespace(points=[])

        argv = [
            str(module.__file__),
            "--name", "test_arf_range",
            "--start", "0.001",
            "--stop", "0.003",
            "--step", "0.0005",
        ]
        with (
            patch.object(module, "Simulation", return_value=object()),
            patch.object(module, "sweep_f2_parallel", side_effect=fake_sweep),
            patch.object(module, "render_standard_sweep_figures", return_value=()),
            patch.object(sys, "argv", argv),
            patch("sys.stdout", new=io.StringIO()),
        ):
            module.main()

        self.assertEqual(call["start"], 0.001)
        self.assertEqual(call["stop"], 0.003)
        self.assertEqual(call["step"], 0.0005)
        self.assertNotIn("values", call)


if __name__ == "__main__":
    unittest.main()
