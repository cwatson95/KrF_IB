#!/usr/bin/env python3
"""Reproduce ArF electron-energy closure costs for stub and MPC drives.

The JSON output records the provisional continuous sink, accepted source-stage chemistry cost,
their measured ratio (when defined), electron growth, equilibrium-temperature residual, limiter
activity, and warmed runtime.  This is causal evidence, not a Chen validation run.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitWaveforms,
    CircuitSpec,
    EnergyClosure,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    PlasmaScalarTraces,
    Preionization,
    Simulation,
)
from Fluid1D.plasma_model import DriftDiffusionPlasmaModel
from particles import excimer_system_registry
from reactions import load_ion_mobility, load_reaction_network


def _simulation(
    name: str, drive: str, closure: str, *, cells: int, coefficient_backend: str,
) -> Simulation:
    fractions = {"Ne": 0.9634, "Ar": 0.035, "F2": 0.0015, "Xe": 0.0001}
    plasma_model = None
    if drive == "stub":
        registry = excimer_system_registry()
        compiled = load_reaction_network(
            _ROOT / "data/reactions/arf.toml", registry,
        ).compile()
        plasma_model = DriftDiffusionPlasmaModel(
            reactions=compiled,
            ion_mobility=load_ion_mobility(_ROOT / "data/ion_transport/ion_mobility.toml"),
            mixture=fractions,
            gas_temp_k=300.0,
            area_m2=90e-4,
            gap_m=12e-3,
            electron_energy_closure=closure,
        )
    return Simulation(
        name=name,
        grid=Grid1D(cells=cells, length=12e-3, area_m2=90e-4),
        gas=GasMixture(fractions, pressure_atm=3.4),
        plasma=Plasma(density=1e16, ion="Ne+", electron_mean_energy_ev=2.0),
        preionization=Preionization(
            mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0,
        ),
        cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),
        boltzmann=BoltzmannSpec(
            backend=coefficient_backend,
            key_kind="mean_energy_eV",
            cache_dir=str(_ROOT / "runs/lut_cache"),
        ),
        circuit=(CircuitSpec(backend="mpc", v0_v=-18e3) if drive == "mpc"
                 else CircuitSpec(backend="stub", phi_left=0.0, phi_right=-18e3)),
        reactions=_ROOT / "data/reactions/arf.toml",
        output_dir=_ROOT / "runs",
        diagnostics=[EnergyClosure(every=1), PlasmaScalarTraces(), CircuitWaveforms()],
        plasma_model=plasma_model,
    )


def _run_case(
    drive: str,
    closure: str,
    *,
    cells: int,
    steps: int,
    dt_s: float,
    coefficient_backend: str = "stub",
) -> dict:
    tag = f"{drive}_{closure}"
    sim = _simulation(
        f"analysis/arf_energy_closure/{coefficient_backend}_{tag}",
        drive,
        closure,
        cells=cells,
        coefficient_backend=coefficient_backend,
    )
    n0 = float(sim.preionization_spec.floor + sim.preionization_spec.amplitude)
    numerics = Numerics(
        steps=steps,
        dt=dt_s,
        fluid_every=max(1, min(10, steps)),
        fluid_substeps=4,
        electron_energy_closure=closure,
    )
    started = time.perf_counter()
    sim.run(numerics, verbose=False)
    runtime = time.perf_counter() - started
    plasma_trace = next(d for d in sim.diagnostics if isinstance(d, PlasmaScalarTraces))
    circuit_trace = next(d for d in sim.diagnostics if isinstance(d, CircuitWaveforms))
    closure_trace = next(d for d in sim.diagnostics if isinstance(d, EnergyClosure))
    totals = closure_trace.totals
    accepted = float(totals["accepted_reaction_cost_j"])
    provisional = float(totals["continuous_inelastic_sink_j"])
    snapshot = sim.energy_closure_snapshot()
    electron_density_m3 = plasma_trace.series("n_e_m3")
    current_a = circuit_trace.current_a
    return {
        "drive": drive,
        "closure": closure,
        "coefficient_backend": coefficient_backend,
        "diagnostic_only": closure == "continuous_reconciled",
        "cells": cells,
        "steps": steps,
        "dt_s": dt_s,
        "runtime_s": runtime,
        "provisional_continuous_sink_j": provisional,
        "accepted_reaction_cost_j": accepted,
        "signed_reconciliation_j": float(totals["signed_reconciliation_j"]),
        "coefficient_bound_adjustment_j": float(totals["coefficient_bound_adjustment_j"]),
        "sink_to_accepted_cost_ratio": provisional / accepted if accepted > 0.0 else None,
        "electron_growth_ratio": float(np.max(sim.charged.electron_density())) / max(n0, 1.0),
        "peak_electron_density_m3": (
            float(np.max(electron_density_m3)) if electron_density_m3.size else 0.0
        ),
        "peak_current_a": float(np.max(np.abs(current_a))) if current_a.size else 0.0,
        "deposited_energy_j": sim.energetics.deposited_energy_j(),
        "output_energy_j": sim.energetics.output_energy_j(),
        "peak_output_power_w": max(sim.energetics.power_out_w, default=0.0),
        "mean_electron_temperature_ev": sim.electron_temperature(),
        "equilibrium_te_residual_ev": float(snapshot["equilibrium_te_residual_ev"]),
        "bulk_equilibrium_te_residual_ev": float(
            snapshot.get("bulk_equilibrium_te_residual_ev", 0.0)
        ),
        "delivered_source_j": float(totals.get("delivered_source_j", 0.0)),
        "undelivered_field_work_j": float(totals.get("undelivered_field_work_j", 0.0)),
        "limiter_activity": int(np.sum(closure_trace.series("limiter_activity"))),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=24)
    parser.add_argument("--steps", type=int, default=40)
    parser.add_argument("--dt-ps", type=float, default=50.0)
    parser.add_argument(
        "--coefficient-backend", choices=("stub", "eedf", "two_term"), default="stub",
    )
    parser.add_argument("--drive", choices=("all", "stub", "mpc"), default="all")
    parser.add_argument(
        "--closure", choices=("all", "source_coupled", "continuous_reconciled"),
        default="all",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_ROOT / "runs/analysis/arf_energy_closure/closure_probe.json",
    )
    args = parser.parse_args(argv)
    # Compile both transport variants before timing.  Stub cases have no JIT transport loop.
    closures = (
        ("source_coupled", "continuous_reconciled")
        if args.closure == "all" else (args.closure,)
    )
    drives = ("stub", "mpc") if args.drive == "all" else (args.drive,)
    if "mpc" in drives:
        for closure in closures:
            _run_case(
                "mpc", closure, cells=args.cells, steps=1, dt_s=args.dt_ps * 1e-12,
                coefficient_backend=args.coefficient_backend,
            )
    cases = [
        _run_case(
            drive,
            closure,
            cells=args.cells,
            steps=args.steps,
            dt_s=args.dt_ps * 1e-12,
            coefficient_backend=args.coefficient_backend,
        )
        for drive in drives
        for closure in closures
    ]
    payload = {
        "schema": "arf-energy-closure-probe/1",
        "coefficient_backend": args.coefficient_backend,
        "note": "Ratios are calculated from the energy ledgers produced by these cases.",
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
