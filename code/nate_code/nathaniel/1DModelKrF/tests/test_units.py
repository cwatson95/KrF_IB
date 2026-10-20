"""Volume/length unit system (cm|m) — conversions, labels, round-trips."""

from __future__ import annotations

import unittest

from Fluid1D.components import Units


class TestUnits(unittest.TestCase):
    def test_length_scale(self) -> None:
        self.assertEqual(Units("m").length_m, 1.0)
        self.assertEqual(Units("cm").length_m, 1.0e-2)

    def test_invalid_unit_raises(self) -> None:
        with self.assertRaises(ValueError):
            Units("mm")

    def test_length_and_area_to_si(self) -> None:
        cm = Units("cm")
        self.assertAlmostEqual(cm.length_to_si(1.2), 1.2e-2)      # 1.2 cm -> 0.012 m
        self.assertAlmostEqual(cm.area_to_si(90.0), 90.0e-4)      # 90 cm^2 -> 9e-3 m^2

    def test_density_conversions(self) -> None:
        cm = Units("cm")
        self.assertAlmostEqual(cm.density_to_si(1.0e12) / 1.0e18, 1.0, places=9)  # cm^-3 -> m^-3
        self.assertAlmostEqual(cm.density_from_si(1.0e18) / 1.0e12, 1.0, places=9)  # m^-3 -> cm^-3
        # round-trip is exact
        self.assertAlmostEqual(cm.density_from_si(cm.density_to_si(3.3e11)) / 3.3e11, 1.0, places=12)

    def test_inv_length_conversions(self) -> None:
        cm = Units("cm")
        self.assertAlmostEqual(cm.inv_length_to_si(1.0), 100.0)   # 1/cm -> 1/m (x100)
        self.assertAlmostEqual(cm.inv_length_from_si(100.0), 1.0)  # 1/m -> 1/cm

    def test_m_is_identity(self) -> None:
        m = Units("m")
        for f, x in [(m.density_to_si, 1e18), (m.density_from_si, 1e18),
                     (m.length_to_si, 0.012), (m.area_to_si, 9e-3),
                     (m.inv_length_to_si, 100.0)]:
            self.assertEqual(f(x), x)

    def test_labels(self) -> None:
        self.assertEqual(Units("cm").density_unit, "cm^-3")
        self.assertEqual(Units("m").density_unit, "m^-3")
        self.assertEqual(Units("cm").length_label, "cm")
        self.assertEqual(Units("cm").inv_length_label, "1/cm")


if __name__ == "__main__":
    unittest.main()
