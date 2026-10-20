"""Tests for Cross-Section Browser import."""

from __future__ import annotations

import math
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from reactions import (
    ReactionKind,
    import_cross_sections,
    interpolate_sigma_linear,
    load_cross_section,
    load_cross_section_directory,
    resolve_cross_section_path,
)
from reactions.cross_sections import (
    DifferentialCrossSectionInterpolator,
    TabulatedSigma,
    TotalCrossSectionInterpolator,
)
from reactions.network import ReactionNetwork

# test fixtures only -- never point these tests at data/cross_sections
# (experiment inputs there may be regenerated at any time; see data/tests/readme.md)
_DATA = _ROOT / "data" / "tests"
_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "cross_sections"
_KR_JSON = _DATA / "example_e_Kr_ionization.json"


class TestCrossSections(unittest.TestCase):
    def test_load_amdis_kr_ionization(self) -> None:
        record = load_cross_section(_KR_JSON)
        self.assertEqual(record.data_kind, "total")
        self.assertEqual(record.reaction_kind, ReactionKind.IONIZATION)
        self.assertEqual(record.threshold_ev, 14.0)
        self.assertEqual(record.energy_ev.size, 220)
        self.assertEqual(record.projectile, "e")
        self.assertEqual(record.target, "Kr")
        self.assertEqual(record.process_id, "example_e_Kr_ionization")
        # the CSV column is cm^2 (units metadata); values must land in SI m^2
        self.assertAlmostEqual(float(record.sigma_m2.max()), 4.1467e-20, delta=1e-23)

    def test_total_interpolation_positive_above_threshold(self) -> None:
        record = load_cross_section(_KR_JSON)
        interp = TotalCrossSectionInterpolator.from_record(record)
        below = float(interp(10.0))
        above = float(interp(20.0))
        self.assertEqual(below, 0.0)
        self.assertGreater(above, 0.0)

    def test_interpolate_sigma_linear_bracketing(self) -> None:
        energy = np.array([10.0, 16.0, 18.0, 20.0], dtype=np.float64)
        sigma = np.array([0.0, 1.0e-20, 3.0e-20, 5.0e-20], dtype=np.float64)
        mid = float(interpolate_sigma_linear(17.0, energy, sigma, threshold_ev=10.0))
        self.assertAlmostEqual(mid, 2.0e-20)

    def test_interpolate_below_threshold_returns_zero(self) -> None:
        energy = np.array([14.0, 20.0], dtype=np.float64)
        sigma = np.array([0.0, 1.0e-20], dtype=np.float64)
        self.assertEqual(
            float(interpolate_sigma_linear(10.0, energy, sigma, threshold_ev=14.0)),
            0.0,
        )

    def test_interpolation_clamps_beyond_measured_range(self) -> None:
        energy = np.array([1.0, 2.0])
        sigma = np.array([3.0e-21, 5.0e-21])
        self.assertEqual(interpolate_sigma_linear(20.0, energy, sigma), sigma[-1])

    def test_interpolate_at_tabulated_energy(self) -> None:
        record = load_cross_section(_KR_JSON)
        interp = TotalCrossSectionInterpolator.from_record(record)
        idx = 5
        e_tab = float(record.energy_ev[idx])
        s_tab = float(record.sigma_m2[idx])
        self.assertTrue(math.isclose(float(interp(e_tab)), s_tab, rel_tol=1e-12))

    def test_resolve_cross_section_by_basename_and_process_id(self) -> None:
        by_basename = resolve_cross_section_path("example_e_Kr_ionization", _DATA)
        self.assertEqual(by_basename, _KR_JSON)
        by_process_id = resolve_cross_section_path("example_e_Kr_ionization", _DATA)
        self.assertEqual(by_process_id, _KR_JSON)

    def test_load_differential_fixture(self) -> None:
        record = load_cross_section(_FIXTURES / "TEST_Differential_00.json")
        self.assertEqual(record.data_kind, "differential")
        self.assertIsNotNone(record.secondary_axis)
        self.assertEqual(record.sigma_m2.shape, (2, 2))
        interp = DifferentialCrossSectionInterpolator.from_record(record)
        self.assertEqual(interp(1.0e30, 1.0e30), 0.0)

    def test_catalog_from_directory(self) -> None:
        catalog = load_cross_section_directory(_DATA, projectile="e", target="Kr")
        self.assertGreaterEqual(len(catalog.records), 1)
        sigma = catalog.total_cross_section(20.0)
        self.assertGreater(float(sigma), 0.0)

    def test_kr_catalog_from_fixture_directory(self) -> None:
        # e/Kr filter picks the ionization + 1s5 excitation fixtures, not F2
        catalog = load_cross_section_directory(_DATA, projectile="e", target="Kr")
        self.assertEqual(catalog.target, "Kr")
        self.assertEqual(len(catalog.records), 2)

    def test_vibrational_metadata_mapping(self) -> None:
        record = load_cross_section(_FIXTURES / "TEST_F2_Vibrational_00.json")
        self.assertEqual(record.reaction_kind, ReactionKind.VIBRATIONAL)
        self.assertEqual(record.threshold_ev, 0.109)


def _kr_ionization_network(
    sigma_key: str | None = "example_e_Kr_ionization",
) -> "ReactionNetwork":
    params: dict = {"threshold_eV": 14.0}
    net = ReactionNetwork()
    net.add("Kr + e", "Kr+ + 2 e", sigma=sigma_key,
            rate_form="electron_impact", params=params, id="5")
    return net


class TestTabulatedSigma(unittest.TestCase):
    def setUp(self) -> None:
        self.record = load_cross_section(_KR_JSON)
        self.tab = TabulatedSigma.from_record(self.record)
        self.ref = TotalCrossSectionInterpolator.from_record(self.record)

    def test_matches_bracketing_interpolation(self) -> None:
        # log-grid LUT must reproduce exact bracketing interpolation to well
        # under the dataset's stated 10% uncertainty (inside the tabulated
        # range only: beyond it the LUT clamps while the reference extrapolates)
        e_lo = float(self.record.energy_ev.min()) * 1.01
        e_hi = float(self.record.energy_ev.max()) * 0.999
        energies = np.logspace(math.log10(e_lo), math.log10(e_hi), 2000)
        ref = np.asarray(self.ref(energies), dtype=np.float64)
        tab = np.asarray(self.tab(energies), dtype=np.float64)
        mask = ref > 0.0
        rel = np.abs(tab[mask] - ref[mask]) / ref[mask]
        self.assertLess(float(rel.max()), 5.0e-3)

    def test_zero_below_threshold_and_clamped_above_table(self) -> None:
        self.assertEqual(float(self.tab(10.0)), 0.0)
        self.assertEqual(float(self.tab(14.0)), 0.0)  # gate is strict: sigma(E0) = 0
        last = float(self.record.sigma_m2[-1])
        self.assertAlmostEqual(float(self.tab(1.0e6)), last, delta=1e-12 * last)

    def test_scalar_and_array_calls_agree(self) -> None:
        energies = np.array([12.0, 20.0, 72.0, 500.0])
        vec = np.asarray(self.tab(energies), dtype=np.float64)
        for e, v in zip(energies, vec, strict=True):
            self.assertEqual(float(self.tab(float(e))), float(v))


class TestCrossSectionImport(unittest.TestCase):
    def test_successful_import_reports_and_overrides(self) -> None:
        report = import_cross_sections(_kr_ionization_network().compile(), _DATA)
        self.assertEqual(len(report.imports), 1)
        imp = report.imports[0]
        self.assertEqual(imp.status, "imported")
        self.assertEqual(imp.key, "example_e_Kr_ionization")
        self.assertIn("successfully imported", imp.describe())
        self.assertIn("[example_e_Kr_ionization]", imp.describe())
        record = load_cross_section(_KR_JSON)
        self.assertIn(record.reaction, imp.detail)  # dataset reaction string, for reference
        self.assertIsInstance(imp.sigma, TabulatedSigma)
        self.assertEqual(report.n_failed, 0)

    def test_missing_dataset_is_reported(self) -> None:
        report = import_cross_sections(
            _kr_ionization_network(sigma_key="does_not_exist").compile(), _DATA
        )
        imp = report.imports[0]
        self.assertEqual(imp.status, "failed")
        self.assertIn("FAILED", imp.detail)
        self.assertIsNone(imp.sigma)

    def test_mismatched_reaction_metadata_is_reported(self) -> None:
        # link the Kr ionization dataset to an F2 attachment reaction: the dataset's
        # projectile/target/kind metadata must reject the pairing at import time
        net = ReactionNetwork()
        net.add("F2 + e", "F + F-", sigma="example_e_Kr_ionization",
                rate_form="electron_impact",
                params={"threshold_eV": 0.1}, id="10")
        report = import_cross_sections(net.compile(), _DATA)
        imp = report.imports[0]
        self.assertEqual(imp.status, "failed")
        self.assertIn("not a reactant", imp.detail)
        self.assertIsNone(imp.sigma)

    def test_kind_mismatch_is_reported(self) -> None:
        # excitation-shaped reaction (no extra electron) linked to an ionization dataset
        net = ReactionNetwork()
        net.add("Kr + e", "Kr* + e", sigma="example_e_Kr_ionization",
                rate_form="electron_impact",
                params={"threshold_eV": 9.95}, id="1")
        report = import_cross_sections(net.compile(), _DATA)
        imp = report.imports[0]
        self.assertEqual(imp.status, "failed")
        self.assertIn("does not match", imp.detail)

    def test_resolve_accepts_csv_filename(self) -> None:
        # requirement: the cross-section csv/json *filename* works as the TOML key
        path = resolve_cross_section_path("example_e_Kr_ionization.csv", _DATA)
        self.assertTrue(path.name.startswith("example_e_Kr_ionization"))
        record = load_cross_section(path)
        self.assertEqual(record.process_id, "example_e_Kr_ionization")

    def test_attachment_dataset_cm2_conversion(self) -> None:
        # Morgan F2 dissociative attachment: cm^2 CSV -> SI m^2 via units metadata
        record = load_cross_section(_DATA / "example_e_F2_attachment.json")
        self.assertEqual(record.reaction_kind, ReactionKind.ATTACHMENT)
        self.assertEqual(record.target, "F2")
        self.assertAlmostEqual(float(record.sigma_m2.max()), 8.0e-19, delta=1e-25)
        tab = TabulatedSigma.from_record(record)
        self.assertAlmostEqual(float(tab(0.01)), 4.44e-19, delta=2e-21)

    def test_same_stem_csv_fallback_for_stale_data_file(self) -> None:
        # a renamed export whose JSON still names the old CSV must load via the
        # same-stem CSV sitting next to the JSON
        import json
        import shutil
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            shutil.copy(_DATA / "example_e_F2_attachment.csv", tmp_path / "renamed.csv")
            meta = json.loads((_DATA / "example_e_F2_attachment.json").read_text())
            meta["data_file"] = "no_longer_exists.csv"
            (tmp_path / "renamed.json").write_text(json.dumps(meta))
            record = load_cross_section(tmp_path / "renamed.json")
            self.assertEqual(record.reaction_kind, ReactionKind.ATTACHMENT)
            self.assertEqual(record.energy_ev.size, 61)


if __name__ == "__main__":
    unittest.main()
