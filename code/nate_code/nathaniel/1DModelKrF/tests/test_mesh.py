"""1D mesh geometry + boundary helpers."""

from __future__ import annotations

import unittest

import numpy as np

from Fluid1D.mesh import Mesh1D


class TestMesh1D(unittest.TestCase):
    def test_geometry(self) -> None:
        m = Mesh1D(10, 1.0, area_m2=2.0)
        self.assertAlmostEqual(m.dx, 0.1)
        self.assertAlmostEqual(m.cell_volume, 0.2)
        self.assertEqual(m.cell_centers.shape, (10,))
        self.assertAlmostEqual(m.cell_centers[0], 0.05)
        self.assertEqual(m.edges.shape, (11,))
        self.assertAlmostEqual(m.edges[-1], 1.0)

    def test_cell_of_periodic(self) -> None:
        m = Mesh1D(10, 1.0)
        idx = m.cell_of(np.array([0.0, 0.15, 1.05, -0.05]))
        # 0.0->0, 0.15->1, 1.05 wraps to 0.05->0, -0.05 wraps to 0.95->9
        np.testing.assert_array_equal(idx, [0, 1, 0, 9])

