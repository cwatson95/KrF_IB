#!/usr/bin/env python3
"""Run the ArF 0.15% F2 spatial and timestep convergence gates.

The production matrix is 128/256/512 cells at 50 ps plus 256 cells at 100/50/25 ps,
all over the same 200 ns physical window. The process exits nonzero when either successive
refinement misses a declared gate. ``--quick`` is a wiring smoke test, not a physics gate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from Fluid1D import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    CircuitWaveforms,
    EnergyClosure,
    ExcimerRates,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    PlasmaScalarTraces,
    Preionization,
    Simulation,
    SpeciesDensities,
)
from Fluid1D.convergence import (
    ConvergenceMetrics,
    compare_refinement,
    pulse_fwhm,
)


def _run(
    cells: int,
    dt_s: float,
    duration_s: float,
    f2: float,
    root: Path,
    *,
    field_heating_cap_factor: float | None = 2.0,
    tag_suffix: str = "",
    electron_energy_closure: str = "source_coupled",
) -> ConvergenceMetrics:
    trace = PlasmaScalarTraces()
    circuit_trace = CircuitWaveforms()
    excimer_rates = ExcimerRates(excimer="ArF*", every=1)
    excimer_density = SpeciesDensities(track=("ArF*",), every=1)
    closure_trace = EnergyClosure(every=1)
    tag = f"c{cells}_dt{dt_s * 1e12:g}ps{tag_suffix}"
    sim = Simulation(
        name=f"analysis/arf_convergence/{tag}",
        grid=Grid1D(cells=cells, length=12e-3, area_m2=90e-4),
        gas=GasMixture(
            {"Ne": 1.0 - 0.035 - 0.0001 - f2, "Ar": 0.035, "F2": f2, "Xe": 0.0001},
            pressure_atm=3.4,
        ),
        plasma=Plasma(density=1e16, ion="Ne+", electron_mean_energy_ev=2.0),
        preionization=Preionization(mode="tanh", amplitude=1e14, floor=1e13, steepness=100.0),
        cavity=Cavity(wavelength_nm=193.0, r_out=0.08, r_back=1.0),
        boltzmann=BoltzmannSpec(
            backend="two_term", key_kind="mean_energy_eV",
            cache_dir=str(root / "runs/lut_cache"),
        ),
        circuit=CircuitSpec(backend="mpc", v0_v=-18e3),
        reactions=root / "data/reactions/arf.toml",
        output_dir=root / "runs",
        diagnostics=[trace, circuit_trace, excimer_rates, excimer_density, closure_trace],
    )
    steps = int(round(duration_s / dt_s))
    sim.run(
        Numerics(
            steps=steps,
            dt=dt_s,
            fluid_every=max(1, round(1e-9 / dt_s)),
            field_heating_cap_factor=field_heating_cap_factor,
            electron_energy_closure=electron_energy_closure,
        ),
        verbose=False,
    )
    t = np.asarray(sim.energetics.times_s, dtype=np.float64)
    p = np.asarray(sim.energetics.power_out_w, dtype=np.float64)
    t_ne = trace.times_ns * 1.0e-9
    ne = trace.series("n_e_m3")
    breakdown = np.flatnonzero(ne >= 1.0e20)
    formation = excimer_rates.formation_rates()
    quenching = excimer_rates.quenching_rates()
    rate_t_s = excimer_rates.times_ns * 1.0e-9
    formation_total = sum((values for values in formation.values()), start=np.zeros(rate_t_s.size))
    quenching_total = sum((values for values in quenching.values()), start=np.zeros(rate_t_s.size))
    closure_residual = closure_trace.series("max_energy_closure_residual_j_m3")
    current_a = circuit_trace.current_a
    te_ev = trace.series("te_ev")
    excimer_m3 = excimer_density.series("ArF*")
    return ConvergenceMetrics(
        cells=cells,
        dt_s=dt_s,
        duration_s=duration_s,
        peak_electron_density_m3=float(np.max(ne)) if ne.size else 0.0,
        deposited_energy_j=sim.energetics.deposited_energy_j(),
        output_energy_j=sim.energetics.output_energy_j(),
        peak_output_power_w=float(np.max(p)) if p.size else 0.0,
        pulse_fwhm_s=pulse_fwhm(t, p),
        chemistry_energy_cost_j=closure_trace.chemistry_energy_cost_j,
        chemistry_bound_adjustment_j=closure_trace.chemistry_bound_adjustment_j,
        f2_mole_fraction=f2,
        breakdown_time_s=float(t_ne[breakdown[0]]) if breakdown.size else 0.0,
        peak_current_a=float(np.max(np.abs(current_a))) if current_a.size else 0.0,
        peak_electron_temperature_ev=float(np.max(te_ev)) if te_ev.size else 0.0,
        peak_excimer_density_m3=float(np.max(excimer_m3)) if excimer_m3.size else 0.0,
        excimer_formation_integral_mol_m3=(
            float(np.trapezoid(formation_total, rate_t_s)) if rate_t_s.size > 1 else 0.0
        ),
        excimer_quenching_integral_mol_m3=(
            float(np.trapezoid(quenching_total, rate_t_s)) if rate_t_s.size > 1 else 0.0
        ),
        max_energy_closure_residual_j_m3=max(closure_residual, default=0.0),
        electron_energy_closure=electron_energy_closure,
    )


def _checkpoint(metric: ConvergenceMetrics) -> None:
    """Persist each expensive member immediately; the aggregate report is written last."""
    tag = f"c{metric.cells}_dt{metric.dt_s * 1e12:g}ps"
    path = _ROOT / f"runs/analysis/arf_convergence/{tag}/convergence_metrics.json"
    path.write_text(json.dumps(metric.as_dict(), indent=2) + "\n", encoding="utf-8")
    print(f"completed {tag}: {json.dumps(metric.as_dict(), sort_keys=True)}", flush=True)


def _load_checkpoint(
    cells: int,
    dt_s: float,
    duration_s: float,
    f2: float,
    electron_energy_closure: str = "source_coupled",
) -> ConvergenceMetrics | None:
    """Load an exact completed member, rejecting stale or differently configured output."""
    tag = f"c{cells}_dt{dt_s * 1e12:g}ps"
    path = _ROOT / f"runs/analysis/arf_convergence/{tag}/convergence_metrics.json"
    if not path.exists():
        return None
    metric = ConvergenceMetrics(**json.loads(path.read_text(encoding="utf-8")))
    expected = (cells, dt_s, duration_s, f2)
    actual = (metric.cells, metric.dt_s, metric.duration_s, metric.f2_mole_fraction)
    if (not all(np.isclose(a, b, rtol=0.0, atol=1e-18) for a, b in zip(actual, expected))
            or metric.electron_energy_closure != electron_energy_closure):
        return None
    print(f"resumed {tag} from {path}", flush=True)
    return metric


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--f2", type=float, default=0.0015)
    parser.add_argument("--duration-ns", type=float, default=200.0)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--closure", choices=("source_coupled", "continuous_reconciled"),
        default="source_coupled",
    )
    args = parser.parse_args()
    duration_s = args.duration_ns * 1e-9
    cells = (32, 64, 128) if args.quick else (128, 256, 512)
    dts = (200e-12, 100e-12, 50e-12) if args.quick else (100e-12, 50e-12, 25e-12)
    if args.quick and args.duration_ns == 200.0:
        duration_s = 5e-9

    runs: dict[tuple[int, float], ConvergenceMetrics] = {}
    for c in cells:
        metric = _load_checkpoint(
            c, 50e-12, duration_s, args.f2, args.closure,
        ) if args.resume else None
        runs[(c, 50e-12)] = metric or _run(
            c, 50e-12, duration_s, args.f2, _ROOT,
            electron_energy_closure=args.closure,
        )
        if metric is None:
            _checkpoint(runs[(c, 50e-12)])
    for dt in dts:
        if (cells[1], dt) not in runs:
            metric = _load_checkpoint(
                cells[1], dt, duration_s, args.f2, args.closure,
            ) if args.resume else None
            runs[(cells[1], dt)] = metric or _run(
                cells[1], dt, duration_s, args.f2, _ROOT,
                electron_energy_closure=args.closure,
            )
            if metric is None:
                _checkpoint(runs[(cells[1], dt)])

    comparisons = {
        f"space_{cells[0]}_{cells[1]}": compare_refinement(
            runs[(cells[0], 50e-12)], runs[(cells[1], 50e-12)],
        ),
        f"space_{cells[1]}_{cells[2]}": compare_refinement(
            runs[(cells[1], 50e-12)], runs[(cells[2], 50e-12)],
        ),
        "time_coarse_medium": compare_refinement(
            runs[(cells[1], dts[0])], runs[(cells[1], dts[1])],
        ),
        "time_medium_fine": compare_refinement(
            runs[(cells[1], dts[1])], runs[(cells[1], dts[2])],
        ),
    }
    payload = {
        "runs": {f"c{c}_dt{dt * 1e12:g}ps": metric.as_dict() for (c, dt), metric in runs.items()},
        "comparisons": {
            name: {"passed": result.passed, "changes": result.changes, "failures": result.failures}
            for name, result in comparisons.items()
        },
    }
    out = _ROOT / "runs/analysis/arf_convergence/convergence.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not args.quick and not all(result.passed for result in comparisons.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
