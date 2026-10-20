"""F2-concentration parameter sweep — the top-level use case.

The program is built around scanning the initial F₂ mole fraction of the gas fill. For each
F₂ value :func:`sweep_f2` rebuilds the :class:`~Fluid1D.components.GasMixture`, builds (and
caches) one Boltzmann :class:`~Fluid1D.boltzmann.CoefficientTable` for that mixture, runs a
fresh :class:`~Fluid1D.simulation.Simulation`, and collects per-point scalars (intrinsic
efficiency, laser output energy) into a :class:`SweepResult`.
"""

from __future__ import annotations

import math
import os
import pickle
import time
import traceback
from collections.abc import Callable, Sequence
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field, replace
from numbers import Integral
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from particles import excimer_system_registry
from reactions import load_reaction_network

from .boltzmann import BoltzmannSolver, CoefficientTable, build_lut, en_grid_from_spec
from .components import Numerics, ProgressReport
from .plasma_model import PlasmaModel
from .progress import format_duration, read_heartbeat, write_heartbeat
from .simulation import Simulation

# Parent-side rollup-line cadence while sweep points are in flight (seconds).
_ROLLUP_PERIOD_S = 10.0


@dataclass
class SweepPoint:
    """One F₂-fraction sweep sample."""

    f2_fraction: float
    output_energy_j: float
    intrinsic_efficiency: float
    out_dir: Path
    wall_s: float = 0.0


@dataclass
class SweepFailure:
    """One F₂ point that raised instead of completing.

    ``error`` is the one-line ``ExcType: message``; ``traceback`` the full worker
    ``traceback.format_exc()`` (empty for pool-level failures where no worker frame
    exists). ``out_dir`` is the point's run directory, which may not exist when the
    failure predates the build.
    """

    f2_fraction: float
    error: str
    traceback: str
    out_dir: Path
    wall_s: float = 0.0


@dataclass
class SweepResult:
    """Collected per-F₂ scalars for the parameter sweep (successes in ``points``,
    recorded per-point errors in ``failures``)."""

    points: list[SweepPoint] = field(default_factory=list)
    failures: list[SweepFailure] = field(default_factory=list)

    def f2(self) -> NDArray[np.float64]:
        return np.array([p.f2_fraction for p in self.points], dtype=np.float64)

    def efficiency(self) -> NDArray[np.float64]:
        return np.array([p.intrinsic_efficiency for p in self.points], dtype=np.float64)

    def output_energy_j(self) -> NDArray[np.float64]:
        return np.array([p.output_energy_j for p in self.points], dtype=np.float64)

    def best(self) -> SweepPoint:
        """The sweep point of maximum intrinsic efficiency."""
        return max(self.points, key=lambda p: p.intrinsic_efficiency)


def f2_range(start: float, stop: float, step: float) -> tuple[float, ...]:
    """Inclusive F₂ sweep values from ``start`` to ``stop`` in increments of ``step``."""
    if not all(math.isfinite(value) for value in (start, stop, step)):
        raise ValueError("start, stop, and step must be finite")
    if step <= 0.0:
        raise ValueError("step must be finite and positive")
    if start > stop:
        raise ValueError("start must be less than or equal to stop")
    if start < 0.0 or stop > 1.0:
        raise ValueError("F2 range bounds must lie in the closed interval [0, 1]")
    n = int(round((stop - start) / step))
    vals = [start + i * step for i in range(n + 1)]
    if vals and vals[-1] > stop + 1e-12:
        vals.pop()
    return tuple(round(v, 12) for v in vals)


def _validated_values(
    *,
    values: Sequence[float] | None,
    start: float | None = None,
    stop: float | None = None,
    step: float | None = None,
) -> tuple[float, ...]:
    """Validate and normalize either an explicit value list or an inclusive range."""
    range_args = (start, stop, step)
    if values is not None and any(value is not None for value in range_args):
        raise ValueError("values=... is mutually exclusive with start/stop/step")
    if values is None:
        if any(value is None for value in range_args):
            raise ValueError("give either values=... or all of start/stop/step")
        normalized = f2_range(float(start), float(stop), float(step))
    else:
        normalized = tuple(float(value) for value in values)

    if not normalized:
        raise ValueError("the F2 sweep must contain at least one value")
    if not all(math.isfinite(value) for value in normalized):
        raise ValueError("F2 sweep values must be finite")
    if not all(0.0 <= value <= 1.0 for value in normalized):
        raise ValueError("F2 sweep values must lie in the closed interval [0, 1]")
    if len(set(normalized)) != len(normalized):
        raise ValueError("F2 sweep values must be unique")
    return normalized


def _validate_factory(name: str, factory: Callable | None, *, picklable: bool) -> None:
    if factory is None:
        return
    if not callable(factory):
        raise TypeError(f"{name} must be callable")
    if picklable:
        try:
            pickle.dumps(factory)
        except Exception as exc:
            raise TypeError(
                f"{name} must be picklable for a process-parallel sweep"
            ) from exc


def _validate_injected_runtime(
    base_config: Simulation,
    *,
    plasma_model_factory: Callable[[float], PlasmaModel] | None,
    coeff_table_factory: Callable[[float], CoefficientTable] | None,
    picklable: bool,
) -> None:
    """Require per-point construction for mutable/injected runtime objects."""
    _validate_factory("plasma_model_factory", plasma_model_factory, picklable=picklable)
    _validate_factory("coeff_table_factory", coeff_table_factory, picklable=picklable)
    if base_config._plasma_model is not None and plasma_model_factory is None:
        raise ValueError(
            "base_config has an injected plasma_model; provide plasma_model_factory(f2) "
            "so every sweep point receives a fresh model"
        )
    if base_config._coeff_table is not None and coeff_table_factory is None:
        raise ValueError(
            "base_config has an injected coeff_table; provide coeff_table_factory(f2) "
            "so each mixture receives an explicit table"
        )


def _base_point_payload(
    base_config: Simulation,
    *,
    rebalance: str,
    diagnostics_factory: Callable | None,
    plots_factory: Callable | None,
    plasma_model_factory: Callable[[float], PlasmaModel] | None,
    coeff_table_factory: Callable[[float], CoefficientTable] | None,
) -> dict:
    """Physical configuration shared by serial and parallel point construction."""
    return dict(
        sweep_name=base_config.name,
        gas=base_config.gas,
        grid=base_config.grid_spec,
        plasma=base_config.plasma_spec,
        preionization=base_config.preionization_spec,
        cavity=base_config.cavity,
        reactions=base_config.reactions_path,
        boltzmann=base_config.boltzmann_spec,
        circuit=base_config.circuit_spec,
        volume_units=base_config.units.volume,
        output_dir=base_config.output_dir,
        cross_sections_dir=base_config.cross_sections_dir,
        ion_mobility_path=base_config.ion_mobility_path,
        drive=base_config.drive,
        rebalance=rebalance,
        diagnostics_factory=diagnostics_factory,
        plots_factory=plots_factory,
        plasma_model_factory=plasma_model_factory,
        coeff_table_factory=coeff_table_factory,
    )


def _simulation_for_point(payload: dict, f2: float, *, name: str) -> Simulation:
    """Build one point while preserving every representable base configuration field."""
    gas = payload["gas"].with_fraction("F2", f2, rebalance=payload["rebalance"])
    diagnostics_factory = payload["diagnostics_factory"]
    plots_factory = payload["plots_factory"]
    plasma_model_factory = payload["plasma_model_factory"]
    coeff_table_factory = payload["coeff_table_factory"]
    coeff_table = payload.get("coeff_table")
    if coeff_table_factory is not None:
        coeff_table = coeff_table_factory(f2)
    return Simulation(
        name=name,
        grid=payload["grid"],
        gas=gas,
        plasma=payload["plasma"],
        preionization=payload["preionization"],
        cavity=payload["cavity"],
        reactions=payload["reactions"],
        boltzmann=payload["boltzmann"],
        circuit=payload["circuit"],
        volume_units=payload["volume_units"],
        drive=payload["drive"],
        diagnostics=(list(diagnostics_factory(f2))
                     if diagnostics_factory is not None else []),
        plots=list(plots_factory(f2)) if plots_factory is not None else [],
        plasma_model=(plasma_model_factory(f2)
                      if plasma_model_factory is not None else None),
        coeff_table=coeff_table,
        output_dir=payload["output_dir"],
        cross_sections_dir=payload["cross_sections_dir"],
        ion_mobility_path=payload["ion_mobility_path"],
    )


def _point_dir_name(f2: float) -> str:
    """Directory name of one sweep point (must match the worker's Simulation name)."""
    return f"f2_{float(f2):.6g}"


def _worker_progress_spec(progress: bool | ProgressReport) -> ProgressReport | None:
    """Normalize the sweep ``progress`` argument into the per-worker reporter spec.

    Workers never own the terminal (the parent is the single writer), so an enabled spec
    is forced to heartbeat-only.
    """
    if progress is False:
        return None
    if progress is True:
        return ProgressReport(terminal=False, heartbeat=True)
    if isinstance(progress, ProgressReport):
        return replace(progress, terminal=False, heartbeat=True)
    raise TypeError(f"progress must be bool or ProgressReport, got {type(progress).__name__}")


def _sweep_worker_init() -> None:  # pragma: no cover - runs in child processes
    """Pin each sweep worker to one numba thread (the sweep is the parallel axis)."""
    try:
        import numba

        numba.set_num_threads(1)
    except Exception:
        pass
    # The parent process already surfaced any orphan-charge warnings once; silence the
    # per-worker repeats (same network, same finding, N workers).
    try:
        import warnings

        from reactions.network import OrphanChargeWarning

        warnings.simplefilter("ignore", OrphanChargeWarning)
    except Exception:
        pass


def _run_sweep_point(payload: dict) -> "SweepPoint | SweepFailure":  # pragma: no cover - child process
    """One F₂ point in a worker process: build the mixture's LUT, run, return scalars.

    Failures come back as a :class:`SweepFailure` sentinel rather than a raised exception:
    raised exceptions must pickle through the executor's result queue (and can poison the
    pool when they don't), while a primitives-only dataclass always transports and carries
    the full worker traceback explicitly.
    """
    t0 = time.monotonic()
    f2 = float(payload["f2"])
    run_dir = Path(payload["output_dir"]) / payload["sweep_name"] / _point_dir_name(f2)
    spec = payload["progress"]
    numerics = payload["numerics"] or Numerics()
    if spec is not None:
        # Entry heartbeat before the (possibly minutes-long) LUT build, so the parent's
        # rollup sees the point as running at 0% instead of queued.
        run_dir.mkdir(parents=True, exist_ok=True)
        write_heartbeat(run_dir / "progress.json", {
            "schema": 1, "status": "running", "step": 0, "steps": numerics.steps,
            "pct": 0.0, "t_sim_s": 0.0, "wall_s": 0.0,
            "updated_unix": time.time(), "fields": {},
        })
    try:
        sim = _simulation_for_point(
            payload,
            f2,
            name=f"{payload['sweep_name']}/{_point_dir_name(f2)}",
        )
        out = sim.run(numerics, progress=spec, verbose=False)
        return SweepPoint(
            f2_fraction=f2,
            output_energy_j=sim.energetics.output_energy_j(),
            intrinsic_efficiency=sim.energetics.intrinsic_efficiency(),
            out_dir=out,
            wall_s=time.monotonic() - t0,
        )
    except Exception as exc:
        wall = time.monotonic() - t0
        tb = traceback.format_exc()
        if spec is not None:
            # Preserve the last running step/pct/fields under the failed status.
            record = read_heartbeat(run_dir / "progress.json") or {
                "schema": 1, "step": 0, "steps": numerics.steps,
                "pct": 0.0, "t_sim_s": 0.0, "fields": {},
            }
            record.update(status="failed", error=tb,
                          updated_unix=time.time(), wall_s=round(wall, 3))
            write_heartbeat(run_dir / "progress.json", record)
        return SweepFailure(f2, f"{type(exc).__name__}: {exc}", tb, run_dir, wall_s=wall)


def _future_outcome(future, f2: float, sweep_dir: Path) -> "SweepPoint | SweepFailure":
    """The worker's sentinel result, with pool-level errors folded into failures."""
    try:
        return future.result()
    except BrokenProcessPool as exc:
        return SweepFailure(
            float(f2),
            f"BrokenProcessPool: {exc} (a worker died abruptly; this point may be the "
            f"victim or collateral)",
            "",
            sweep_dir / _point_dir_name(f2),
        )
    except Exception as exc:  # e.g. result transport errors — never abort the sweep
        return SweepFailure(float(f2), f"{type(exc).__name__}: {exc}",
                            traceback.format_exc(), sweep_dir / _point_dir_name(f2))


def _completion_line(done_n: int, n: int, point: SweepPoint, elapsed_s: float) -> str:
    eta = elapsed_s / done_n * (n - done_n)
    return (f"[sweep] {done_n}/{n}  f2={point.f2_fraction:.6g}  "
            f"E_out={point.output_energy_j:.3e} J  eff={point.intrinsic_efficiency:.3e}  "
            f"wall={format_duration(point.wall_s)}  eta~{format_duration(eta)}")


def _failure_line(done_n: int, n: int, failure: SweepFailure) -> str:
    first = failure.error.splitlines()[0] if failure.error else "unknown error"
    return (f"[sweep] {done_n}/{n}  f2={failure.f2_fraction:.6g}  "
            f"FAILED after {format_duration(failure.wall_s)}: {first}")


def _running_points(sweep_dir: Path, remaining_f2: Sequence[float],
                    t0_wall: float) -> list[float]:
    """pct values of the not-yet-collected points with fresh ``running`` heartbeats.

    Heartbeats older than the sweep start (1 s slack) are leftovers from a previous
    invocation of the same sweep name and count as queued.
    """
    pcts: list[float] = []
    for f2 in remaining_f2:
        record = read_heartbeat(Path(sweep_dir) / _point_dir_name(f2) / "progress.json")
        if (record is not None and record.get("schema") == 1
                and record.get("status") == "running"
                and float(record.get("updated_unix", 0.0)) >= t0_wall - 1.0):
            pcts.append(float(record.get("pct", 0.0)))
    return pcts


def _rollup_line(done_n: int, n: int, running_pcts: Sequence[float], queued: int,
                 elapsed_s: float) -> str:
    effective = done_n + sum(running_pcts) / 100.0
    eta = (f"eta~{format_duration(elapsed_s * (n - effective) / effective)}"
           if effective > 0.0 else "eta n/a")
    running = (f"running {len(running_pcts)} "
               f"(avg {sum(running_pcts) / len(running_pcts):.0f}%)"
               if running_pcts else "running 0")
    return f"[sweep] done {done_n}/{n}  {running}  queued {queued}  {eta}"


def _print_failure_summary(result: SweepResult, n: int, elapsed_s: float,
                           sweep_dir: Path) -> None:
    if result.failures:
        print(f"[sweep] {len(result.failures)} of {n} points FAILED:", flush=True)
        for failure in result.failures:
            first = failure.error.splitlines()[0] if failure.error else "unknown error"
            print(f"[sweep]   f2={failure.f2_fraction:.6g}  {first}", flush=True)
        print("[sweep] tracebacks: SweepResult.failures; "
              "per-point progress.json has status='failed'", flush=True)
    print(f"[sweep] complete: {len(result.points)} ok, {len(result.failures)} failed "
          f"in {format_duration(elapsed_s)} -> {sweep_dir}", flush=True)


def sweep_f2_parallel(
    base_config: Simulation,
    *,
    start: float | None = None,
    stop: float | None = None,
    step: float | None = None,
    values: Sequence[float] | None = None,
    numerics: Numerics | None = None,
    rebalance: str = "Kr",
    max_workers: int | None = None,
    diagnostics_factory=None,
    plots_factory=None,
    plasma_model_factory: Callable[[float], PlasmaModel] | None = None,
    coeff_table_factory: Callable[[float], CoefficientTable] | None = None,
    progress: bool | ProgressReport = True,
) -> SweepResult:
    """Embarrassingly-parallel F₂ sweep in worker-process batches.

    Select the sweep either by ``(start, stop, step)`` (inclusive range) or an explicit
    ``values`` list. Each point runs a fresh :class:`Simulation` in its own process (workers
    pinned to one numba thread each — the sweep itself is the parallel axis) and writes its
    outputs to ``runs/<name>/f2_<value>/`` (including ``diagnostics.h5`` when a
    ``diagnostics_factory`` — a picklable ``f2 -> [Diagnostic, ...]`` callable — is given,
    which the sweep-overlay plots consume). A separate picklable ``plots_factory`` lets the
    experiment file choose per-point filenames, datasets, and scales without putting those
    choices in the simulation runtime. Injected runtime objects are deliberately not reused:
    pass picklable ``plasma_model_factory(f2)`` or ``coeff_table_factory(f2)`` callables to
    construct fresh per-point instances. A base configuration containing either injection is
    rejected when its corresponding factory is omitted.

    ``progress`` drives the live reporting. The default ``True`` prints a start banner,
    one completion line per finished point, a periodic in-flight rollup, and an end
    summary, while each worker rewrites ``f2_<value>/progress.json`` heartbeats on the
    default cadence (watch them with ``python -m Fluid1D.sweep_status runs/<name>``). Pass
    a :class:`~Fluid1D.components.ProgressReport` to choose the worker cadence/fields
    (terminal output stays parent-owned: the spec is forced to heartbeat-only), or
    ``False`` for a fully silent sweep.

    A point that raises is recorded on ``SweepResult.failures`` (message + worker
    traceback) and the remaining points complete; failed points leave no
    ``diagnostics.h5``, so the sweep-overlay figures render from the survivors. If a
    worker process dies outright (``BrokenProcessPool``), every not-yet-collected point is
    recorded as failed — the true victim is the one whose heartbeat is stuck ``running``.
    """
    values = _validated_values(values=values, start=start, stop=stop, step=step)
    if max_workers is not None:
        if isinstance(max_workers, bool) or not isinstance(max_workers, Integral):
            raise TypeError("max_workers must be an integer or None")
        if max_workers <= 0:
            raise ValueError("max_workers must be positive")
        max_workers = int(max_workers)
    _validate_factory("diagnostics_factory", diagnostics_factory, picklable=True)
    _validate_factory("plots_factory", plots_factory, picklable=True)
    _validate_injected_runtime(
        base_config,
        plasma_model_factory=plasma_model_factory,
        coeff_table_factory=coeff_table_factory,
        picklable=True,
    )
    worker_spec = _worker_progress_spec(progress)
    reporting = progress is not False
    numerics_eff = numerics or Numerics()
    base_payload = _base_point_payload(
        base_config,
        rebalance=rebalance,
        diagnostics_factory=diagnostics_factory,
        plots_factory=plots_factory,
        plasma_model_factory=plasma_model_factory,
        coeff_table_factory=coeff_table_factory,
    )
    payloads = [dict(base_payload, f2=f2, numerics=numerics, progress=worker_spec)
                for f2 in values]
    workers = max_workers if max_workers is not None else max(1, (os.cpu_count() or 2) - 1)
    sweep_dir = Path(base_config.output_dir) / base_config.name
    t0_mono, t0_wall = time.monotonic(), time.time()
    if reporting and values:
        print(f"[sweep] {base_config.name}: {len(values)} points, "
              f"f2 {min(values):.6g}..{max(values):.6g}, {workers} workers, "
              f"{numerics_eff.steps} steps @ dt={numerics_eff.dt:.1e}s -> {sweep_dir}",
              flush=True)
    result = SweepResult()
    n = len(payloads)
    done_n = 0
    with ProcessPoolExecutor(max_workers=workers, initializer=_sweep_worker_init) as pool:
        fut_to_f2 = {pool.submit(_run_sweep_point, p): float(p["f2"]) for p in payloads}
        pending = set(fut_to_f2)
        last_rollup = t0_mono
        try:
            while pending:
                done, pending = wait(pending, timeout=_ROLLUP_PERIOD_S,
                                     return_when=FIRST_COMPLETED)
                for future in done:
                    done_n += 1
                    outcome = _future_outcome(future, fut_to_f2[future], sweep_dir)
                    if isinstance(outcome, SweepFailure):
                        result.failures.append(outcome)
                        if reporting:
                            print(_failure_line(done_n, n, outcome), flush=True)
                    else:
                        result.points.append(outcome)
                        if reporting:
                            print(_completion_line(done_n, n, outcome,
                                                   time.monotonic() - t0_mono), flush=True)
                now = time.monotonic()
                if reporting and pending and now - last_rollup >= _ROLLUP_PERIOD_S:
                    remaining = sorted(fut_to_f2[f] for f in pending)
                    pcts = _running_points(sweep_dir, remaining, t0_wall)
                    print(_rollup_line(done_n, n, pcts, len(remaining) - len(pcts),
                                       now - t0_mono), flush=True)
                    last_rollup = now
        except KeyboardInterrupt:
            pool.shutdown(wait=False, cancel_futures=True)
            if reporting:
                print(f"[sweep] interrupted: {done_n}/{n} finished, "
                      f"{len(pending)} cancelled or aborted", flush=True)
            raise
    result.points.sort(key=lambda p: p.f2_fraction)
    result.failures.sort(key=lambda p: p.f2_fraction)
    if reporting:
        _print_failure_summary(result, n, time.monotonic() - t0_mono, sweep_dir)
    return result


def sweep_f2(
    f2_values: Sequence[float],
    base_config: Simulation,
    *,
    numerics: Numerics | None = None,
    solver: BoltzmannSolver | None = None,
    cache_dir: Path | str | None = None,
    rebalance: str = "Kr",
    verbose: bool = False,
    diagnostics_factory: Callable | None = None,
    plots_factory: Callable | None = None,
    plasma_model_factory: Callable[[float], PlasmaModel] | None = None,
    coeff_table_factory: Callable[[float], CoefficientTable] | None = None,
) -> SweepResult:
    """Sweep the F₂ mole fraction, running one :class:`Simulation` per value.

    ``base_config`` supplies the fixed specs (grid / plasma / cavity / circuit / boltzmann /
    reactions); only the gas mixture's F₂ fraction varies (the delta is absorbed into the
    ``rebalance`` species). One Boltzmann LUT is built per mixture (cached on disk when
    ``cache_dir`` is set). Diagnostics and plots default to empty; optional factories use the
    same per-point construction contract as :func:`sweep_f2_parallel`. The scalars come from
    ``sim.energetics``.
    """
    values = _validated_values(values=f2_values)
    _validate_factory("diagnostics_factory", diagnostics_factory, picklable=False)
    _validate_factory("plots_factory", plots_factory, picklable=False)
    _validate_injected_runtime(
        base_config,
        plasma_model_factory=plasma_model_factory,
        coeff_table_factory=coeff_table_factory,
        picklable=False,
    )
    if solver is not None and coeff_table_factory is not None:
        raise ValueError("solver and coeff_table_factory are mutually exclusive")

    registry = excimer_system_registry()
    compiled = load_reaction_network(base_config.reactions_path, registry).compile()
    en_grid = en_grid_from_spec(base_config.boltzmann_spec.en_grid_td)
    boltzmann = base_config.boltzmann_spec
    if cache_dir is not None:
        boltzmann = replace(boltzmann, cache_dir=str(cache_dir))

    payload = _base_point_payload(
        base_config,
        rebalance=rebalance,
        diagnostics_factory=diagnostics_factory,
        plots_factory=plots_factory,
        plasma_model_factory=plasma_model_factory,
        coeff_table_factory=coeff_table_factory,
    )
    payload["boltzmann"] = boltzmann

    result = SweepResult()
    for f2 in values:
        table = None
        if solver is not None:
            gas = base_config.gas.with_fraction("F2", f2, rebalance=rebalance)
            table = build_lut(
                gas.normalized_fractions, solver=solver, en_grid_td=en_grid,
                gas_temp_k=gas.temperature_k, reactions=compiled,
                cache_dir=boltzmann.cache_dir, key_kind=boltzmann.key_kind,
            )
        point_payload = dict(payload, coeff_table=table)
        sim = _simulation_for_point(
            point_payload,
            float(f2),
            name=f"{base_config.name}_f2_{float(f2):.4f}",
        )
        out = sim.run(numerics, verbose=verbose)
        result.points.append(SweepPoint(
            f2_fraction=float(f2),
            output_energy_j=sim.energetics.output_energy_j(),
            intrinsic_efficiency=sim.energetics.intrinsic_efficiency(),
            out_dir=out,
        ))
    return result
