"""Ion mobility/diffusion loader (reactions.ion_transport)."""

from __future__ import annotations

import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
from particles.constants import E_CHARGE, K_B
from reactions import load_ion_mobility
from reactions.ion_transport import LOSCHMIDT_M3, IonMobilityRecord, IonMobilityTable

TOML = _ROOT / "data" / "ion_transport" / "ion_mobility.toml"


class TestIonMobilityLoader(unittest.TestCase):
    def test_loads_real_table(self) -> None:
        t = load_ion_mobility(TOML)
        self.assertTrue(t.has_ion("Kr+"))
        self.assertTrue(t.has_ion("F-"))
        # μ·N invariant = K0 * 1e-4 * N0 for Kr+ (K0 = 0.84 cm^2/Vs).
        muN = t.reduced_mobility("Kr+", {"Kr": 1.0})
        self.assertAlmostEqual(muN / (0.84e-4 * LOSCHMIDT_M3), 1.0, places=9)

    def test_reduced_mobility_scales_with_density(self) -> None:
        t = load_ion_mobility(TOML)
        mix = {"Kr": 1.0}
        # μ = K0 * N0 / N ; at N = N0, μ = K0 (in SI).
        mu = float(t.mobility("Kr+", LOSCHMIDT_M3, mix))
        self.assertAlmostEqual(mu / (0.84e-4), 1.0, places=9)

    def test_einstein_diffusion(self) -> None:
        t = load_ion_mobility(TOML)
        n, tg = 8.3e25, 300.0
        mu = float(t.mobility("Kr+", n, {"Kr": 1.0}))
        d = float(t.diffusion("Kr+", n, tg, {"Kr": 1.0}))
        self.assertAlmostEqual(d, mu * K_B * tg / E_CHARGE, places=18)

    def test_blancs_law_single_gas(self) -> None:
        # F- in a Kr/F2 mix: only Kr has an F- record -> uses F-/Kr.
        t = load_ion_mobility(TOML)
        muN_mix = t.reduced_mobility("F-", {"Kr": 0.999, "F2": 0.001})
        muN_kr = t.reduced_mobility("F-", {"Kr": 1.0})
        self.assertAlmostEqual(muN_mix, muN_kr, places=6)

    def test_blancs_law_two_gases(self) -> None:
        # Blanc's law: 1/μ = Σ x_g/μ_g. Construct a synthetic two-gas ion.
        recs = (
            IonMobilityRecord("X+", "A", 2.0, 2.0e-4 * LOSCHMIDT_M3, "s"),
            IonMobilityRecord("X+", "B", 4.0, 4.0e-4 * LOSCHMIDT_M3, "s"),
        )
        t = IonMobilityTable(recs)
        muN = t.reduced_mobility("X+", {"A": 0.5, "B": 0.5})
        inv = 0.5 / (2.0e-4 * LOSCHMIDT_M3) + 0.5 / (4.0e-4 * LOSCHMIDT_M3)
        self.assertAlmostEqual(muN, 1.0 / inv, places=6)

    def test_unknown_ion_raises(self) -> None:
        t = load_ion_mobility(TOML)
        with self.assertRaises(KeyError):
            t.reduced_mobility("Unobtanium+", {"Kr": 1.0})


if __name__ == "__main__":
    unittest.main()
