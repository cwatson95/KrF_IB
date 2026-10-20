"""Tests for the runtime compute backend (CPU topology + thread/layer selection).

The suite exercises the topology detection and thread-count heuristic that drive Numba/FFT
parallelism, plus a BLAS-agnostic linear-algebra sanity check.
"""

from __future__ import annotations

import unittest

import numpy as np

import backend
from backend.threading import (
    available_threading_layers,
    optimal_thread_count,
    select_threading_layer,
)
from backend.topology import CPUTopology, detect_topology


def _topo(logical, physical, perf, eff, *, alloc=None, system="linux", arch="x86_64"):
    return CPUTopology(
        logical_cores=logical,
        physical_cores=physical,
        performance_cores=perf,
        efficiency_cores=eff,
        platform=system,
        arch=arch,
        allocation_limit=alloc,
    )


class TestTopologyDetection(unittest.TestCase):
    def test_detect_returns_sane_values(self) -> None:
        t = detect_topology()
        self.assertGreaterEqual(t.logical_cores, 1)
        self.assertGreaterEqual(t.physical_cores, 1)
        self.assertGreaterEqual(t.performance_cores, 0)
        self.assertGreaterEqual(t.efficiency_cores, 0)
        self.assertIsInstance(t.describe(), str)

    def test_heterogeneous_flag(self) -> None:
        self.assertTrue(_topo(8, 8, 4, 4).is_heterogeneous)
        self.assertFalse(_topo(8, 8, 0, 0).is_heterogeneous)


class TestThreadingLayer(unittest.TestCase):
    def test_workqueue_always_available(self) -> None:
        layers = available_threading_layers()
        self.assertIn("workqueue", layers)

    def test_select_prefers_tbb_then_omp_then_workqueue(self) -> None:
        self.assertEqual(select_threading_layer(["workqueue"]), "workqueue")
        self.assertEqual(select_threading_layer(["workqueue", "omp"]), "omp")
        self.assertEqual(select_threading_layer(["workqueue", "omp", "tbb"]), "tbb")

    def test_selected_layer_is_available(self) -> None:
        self.assertIn(select_threading_layer(), available_threading_layers())


class TestThreadHeuristic(unittest.TestCase):
    def test_explicit_override_wins(self) -> None:
        self.assertEqual(optimal_thread_count(_topo(8, 8, 4, 4), override=3), 3)

    def test_asymmetric_uses_performance_cores(self) -> None:
        # 4P + 4E -> use the 4 performance cores (efficiency cores straggle).
        self.assertEqual(optimal_thread_count(_topo(8, 8, 4, 4)), 4)

    def test_homogeneous_uses_physical_cores(self) -> None:
        self.assertEqual(optimal_thread_count(_topo(16, 16, 0, 0)), 16)

    def test_allocation_limit_clamps(self) -> None:
        self.assertEqual(optimal_thread_count(_topo(64, 64, 0, 0, alloc=8)), 8)

    def test_floor_at_one(self) -> None:
        self.assertGreaterEqual(optimal_thread_count(_topo(1, 1, 0, 0)), 1)


class TestConfigure(unittest.TestCase):
    def test_configure_returns_runtime_config(self) -> None:
        cfg = backend.configure(threads=2, force=True)
        self.assertEqual(cfg.n_threads, 2)
        self.assertEqual(cfg.fft_workers, 2)
        self.assertIn(cfg.layer, (*available_threading_layers(), "none"))
        self.assertIs(backend.get_config(), cfg)

    def test_default_configure_is_positive(self) -> None:
        cfg = backend.configure(force=True)
        self.assertGreaterEqual(cfg.n_threads, 1)


class TestLinearAlgebraSanity(unittest.TestCase):
    """BLAS-agnostic: correctness only, no backend-name assertion."""

    def test_matrix_multiply(self) -> None:
        rng = np.random.default_rng(0)
        a = rng.standard_normal((128, 128))
        b = rng.standard_normal((128, 128))
        np.testing.assert_allclose(a @ b, np.dot(a, b), rtol=1e-10, atol=1e-10)

    def test_linear_solve(self) -> None:
        rng = np.random.default_rng(1)
        a = rng.standard_normal((64, 64))
        a = a @ a.T + 64.0 * np.eye(64)
        x = rng.standard_normal(64)
        np.testing.assert_allclose(np.linalg.solve(a, a @ x), x, rtol=1e-9, atol=1e-9)


if __name__ == "__main__":
    unittest.main()
