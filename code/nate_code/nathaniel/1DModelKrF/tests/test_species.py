"""Tests for the species representation (electrons, ions, excited states, excimers)."""

from __future__ import annotations

import unittest

import numpy as np

from particles import (
    ELECTRON,
    F2,
    KR,
    KR_2STAR,
    KR_PLUS,
    KR_STAR,
    KRF_EXCIMER,
    M_F,
    M_KR,
    IonSpecies,
    SpeciesCategory,
    krf_system_registry,
    molecular_mass,
)
from particles.constants import E_CHARGE


class TestSpeciesModel(unittest.TestCase):
    def test_charge_and_field_coupling(self) -> None:
        self.assertEqual(ELECTRON.charge_number, -1)
        self.assertAlmostEqual(KR_PLUS.charge, E_CHARGE)
        self.assertTrue(KR_PLUS.is_field_coupled)
        self.assertFalse(KR.is_field_coupled)
        self.assertFalse(KRF_EXCIMER.is_charged)

    def test_excitation_flags(self) -> None:
        self.assertFalse(KR.is_excited)
        self.assertTrue(KR_STAR.is_excited)
        self.assertGreater(KR_2STAR.excitation_energy_ev, KR_STAR.excitation_energy_ev)
        self.assertTrue(KR_STAR.metastable)
        self.assertEqual(KR_STAR.relaxes_to, "Kr")

    def test_molecule_and_excimer_composition(self) -> None:
        self.assertTrue(F2.is_molecule)
        self.assertEqual(F2.constituents, ("F", "F"))
        self.assertFalse(KR.is_molecule)
        self.assertTrue(KRF_EXCIMER.is_molecule)
        self.assertEqual(KRF_EXCIMER.category, SpeciesCategory.EXCIMER)
        self.assertTrue(KRF_EXCIMER.is_excited)

    def test_molecular_mass_matches_constituents(self) -> None:
        self.assertAlmostEqual(F2.mass, 2.0 * M_F)
        self.assertAlmostEqual(KRF_EXCIMER.mass, M_KR + M_F)
        self.assertAlmostEqual(molecular_mass("Kr", "F"), M_KR + M_F)


class TestIonSpeciesBackwardCompat(unittest.TestCase):
    def test_positional_constructor(self) -> None:
        ion = IonSpecies("Kr+", M_KR)
        self.assertEqual(ion.name, "Kr+")
        self.assertEqual(ion.mass, M_KR)
        self.assertEqual(ion.charge_number, 1)
        self.assertEqual(ion.category, SpeciesCategory.ION)

    def test_explicit_charge_number(self) -> None:
        ion = IonSpecies("Kr2+", M_KR, 2)
        self.assertEqual(ion.charge_number, 2)


class TestSpeciesRegistry(unittest.TestCase):
    def test_krf_system_lookup_and_grouping(self) -> None:
        registry = krf_system_registry()
        self.assertIn("KrF*", registry)
        self.assertIn("Ne*", registry)
        self.assertIn("F2+", registry)
        self.assertIn("hv248", registry)
        self.assertIs(registry.get("Kr*"), KR_STAR)
        excited_names = {s.name for s in registry.excited()}
        self.assertLessEqual(
            {"Kr*", "Kr**", "KrF*", "Kr2*", "Kr2F*", "Ne*", "Ne**", "Ne2*"},
            excited_names,
        )
        charged_names = {s.name for s in registry.charged()}
        self.assertEqual(
            charged_names,
            {"e-", "Kr+", "Kr2+", "F-", "F2+", "Ne+", "Ne2+"},
        )
        self.assertEqual(
            {s.name for s in registry.by_category(SpeciesCategory.EXCIMER)},
            {"F2(a3)", "KrF*", "Kr2*", "Kr2F*", "Ne2*"},
        )

    def test_duplicate_rejected(self) -> None:
        registry = krf_system_registry()
        with self.assertRaises(ValueError):
            registry.add(KR_STAR)

    def test_unknown_lookup_raises(self) -> None:
        with self.assertRaises(KeyError):
            krf_system_registry().get("Xe*")


if __name__ == "__main__":
    unittest.main()
