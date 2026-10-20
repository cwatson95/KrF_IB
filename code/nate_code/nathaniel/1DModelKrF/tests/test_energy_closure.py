"""Electron-energy closure interfaces and signed diagnostic reconciliation."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    EnergyClosure,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    Preionization,
    Simulation,
    reconcile_continuous_energy,
)
from diagnostics.archive import DiagnosticsArchive
from particles.constants import E_CHARGE


class TestClosureInterface(unittest.TestCase):
    def test_production_default_and_diagnostic_alternative(self) -> None:
        self.assertEqual(Numerics().electron_energy_closure, "source_coupled")
        self.assertEqual(
            Numerics(electron_energy_closure="continuous_reconciled").electron_energy_closure,
            "continuous_reconciled",
        )
        with self.assertRaisesRegex(ValueError, "electron_energy_closure"):
            Numerics(electron_energy_closure="invalid")

    def test_diagnostic_archive_marks_continuous_mode_and_totals(self) -> None:
        collector = EnergyClosure(every=1)

        class FakeSimulation(SimpleNamespace):
            def energy_closure_snapshot(self):
                return {
                    "closure_mode": "continuous_reconciled",
                    "field_work_j": 1.0,
                    "limiter_activity": 2.0,
                    "max_ampere_residual": 0.01,
                }

        sim = FakeSimulation(_energy_closure_totals={"field_work_j": 1.0})
        collector.record(sim, 1, 1.0e-9)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "closure.h5"
            with DiagnosticsArchive(path) as archive:
                collector.export(sim, archive)
            with h5py.File(path, "r") as h5:
                attrs = h5["energy_closure"].attrs
                self.assertEqual(attrs["closure_mode"], "continuous_reconciled")
                self.assertFalse(attrs["production_closure"])
                self.assertIn("diagnostic only", attrs["status"])
                self.assertEqual(attrs["total_limiter_activity"], 2.0)
                self.assertEqual(attrs["total_field_work_j"], 1.0)

    def test_chen_bank_target_is_derived_from_published_denominator(self) -> None:
        fixture = json.loads(
            (_ROOT / "data/reference/chen_2025_arf.json").read_text(encoding="utf-8")
        )
        case = fixture["case"]
        expected = case["bank_energy_j"] * case["bank_efficiency"]["value"]
        self.assertAlmostEqual(expected, 0.09396, places=12)
        self.assertAlmostEqual(case["laser_output_energy_j"]["value"], expected, places=12)


class TestStageSourceDelivery(unittest.TestCase):
    """End-to-end source_coupled delivery: Σ delivered ≡ field_work − elastic exactly."""

    def test_delivered_source_matches_the_transport_ledgers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sim = Simulation(
                name="stage_delivery",
                grid=Grid1D(cells=16, length=12e-3, area_m2=90e-4),
                gas=GasMixture({"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
                               pressure_atm=3.4, temperature_k=300.0),
                plasma=Plasma(density=1e15, ion="Ne+", electron_mean_energy_ev=3.0),
                preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13,
                                            steepness=100.0),
                cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),
                boltzmann=BoltzmannSpec(backend="stub", key_kind="mean_energy_eV"),
                circuit=CircuitSpec(backend="mpc", v0_v=-6.0e3),
                reactions=str(_ROOT / "data" / "reactions" / "arf.toml"),
                output_dir=Path(tmp),
                diagnostics=[],
            )
            sim.run(Numerics(steps=24, dt=5.0e-11, fluid_every=8), verbose=False)
            self.assertTrue(sim._real_transport)
            totals = sim._energy_closure_totals
            # Per-interval residual delivery makes the run-integrated identity exact:
            # delivered ≡ rectified (field_work − elastic); the discarded phantom-negative
            # part is quarantined in undelivered_field_work_j (≤ 0), and the three ledgers
            # reproduce the signed books.
            signed = totals["field_work_j"] - totals["elastic_transfer_j"]
            recovered = totals["delivered_source_j"] + totals["undelivered_field_work_j"]
            self.assertAlmostEqual(
                recovered, signed, delta=max(abs(signed), 1e-30) * 1e-9,
            )
            self.assertLessEqual(totals["undelivered_field_work_j"], 0.0)
            self.assertGreaterEqual(totals["delivered_source_j"], 0.0)
            live = sim.last_energy_closure_report
            np.testing.assert_allclose(
                live.stage_delivered_j_m3,
                np.maximum(
                    live.field_work_j_m3 - live.elastic_transfer_j_m3, 0.0,
                ),
                rtol=1e-9, atol=1e-25,
            )


class TestCrossModeAgreement(unittest.TestCase):
    """source_coupled vs continuous_reconciled must describe the SAME plasma.

    Pre-repair the reconciled mode collapsed x13 against production (its lump-sum refund
    was destroyed above the ceiling by the stage pre-clamp).  With the refund delivered as
    an in-stage rate, the two closures are alternative splittings of one physics and a
    short drive must land within a modest factor.
    """

    def _run(self, closure: str) -> tuple[float, float]:
        with tempfile.TemporaryDirectory() as tmp:
            sim = Simulation(
                name=f"xmode_{closure}",
                grid=Grid1D(cells=16, length=12e-3, area_m2=90e-4),
                gas=GasMixture({"Ne": 0.9639, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
                               pressure_atm=3.4, temperature_k=300.0),
                plasma=Plasma(density=1e15, ion="Ne+", electron_mean_energy_ev=3.0),
                preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13,
                                            steepness=100.0),
                cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),
                boltzmann=BoltzmannSpec(backend="stub", key_kind="mean_energy_eV"),
                circuit=CircuitSpec(backend="mpc", v0_v=-8.0e3),
                reactions=str(_ROOT / "data" / "reactions" / "arf.toml"),
                output_dir=Path(tmp),
                diagnostics=[],
            )
            sim.run(Numerics(steps=60, dt=5.0e-11, fluid_every=10,
                             electron_energy_closure=closure), verbose=False)
            return (
                float(sim.charged.electron_density().mean()),
                float(sim.charged.electron_temperature_ev().mean()),
            )

    def test_modes_agree_within_a_modest_factor(self) -> None:
        ne_src, te_src = self._run("source_coupled")
        ne_rec, te_rec = self._run("continuous_reconciled")
        self.assertGreater(ne_src, 0.0)
        self.assertGreater(ne_rec, 0.0)
        ratio = ne_rec / ne_src
        self.assertGreater(ratio, 0.5, f"reconciled n_e collapsed: ratio={ratio:.3g}")
        self.assertLess(ratio, 2.0, f"reconciled n_e inflated: ratio={ratio:.3g}")
        # The two splittings apply heating on different cadences, so a short -8 kV kick
        # leaves an O(1 eV) transient T_e skew (measured ~1.3 eV); the failure class this
        # guards against is the pre-repair multi-eV starvation collapse.
        self.assertLess(abs(te_rec - te_src), 2.0)


class TestSignedReconciliation(unittest.TestCase):
    def test_refund_and_additional_debit_keep_their_sign(self) -> None:
        result = reconcile_continuous_energy(
            transport_j_m3=np.array([10.0, 10.0]),
            prepaid_sink_j_m3=np.array([4.0, 1.0]),
            accepted_cost_j_m3=np.array([1.0, 3.0]),
            electron_density_m3=np.full(2, 1.0 / E_CHARGE),
            energy_bounds_ev=(0.0, 1.0e30),
            ev_to_j=E_CHARGE,
        )
        np.testing.assert_array_equal(result.signed_reconciliation_j_m3, [3.0, -2.0])
        np.testing.assert_array_equal(result.final_j_m3, [13.0, 8.0])
        np.testing.assert_array_equal(result.bound_adjustment_j_m3, [0.0, 0.0])

    def test_floor_intervention_is_a_separate_ledger_term(self) -> None:
        result = reconcile_continuous_energy(
            transport_j_m3=np.array([1.0]),
            prepaid_sink_j_m3=np.array([0.0]),
            accepted_cost_j_m3=np.array([2.0]),
            electron_density_m3=np.array([1.0 / E_CHARGE]),
            energy_bounds_ev=(2.0, 20.0),
            ev_to_j=E_CHARGE,
        )
        self.assertEqual(result.signed_reconciliation_j_m3[0], -2.0)
        self.assertEqual(result.final_j_m3[0], 2.0)
        self.assertEqual(result.bound_adjustment_j_m3[0], 3.0)


if __name__ == "__main__":
    unittest.main()
