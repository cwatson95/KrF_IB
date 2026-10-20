"""Terminal progress-line and heartbeat-file reporting for :meth:`Simulation.run`.

A :class:`ProgressReporter` binds a :class:`~Fluid1D.components.ProgressReport` spec to a
built simulation and delivers one report per due step to the enabled sinks: a terminal
line (``spec.terminal``) and/or an atomically rewritten ``<run_dir>/progress.json``
heartbeat (``spec.heartbeat``). Every field reads existing engine state. Field and species
names are validated up front so a typo fails before the first timestep.

Vocabulary (``<sp>`` is any registered species name):

===============  ==============================================================
``total:<sp>``   real particles (charged/neutral fluid total, or photon mode)
``density:<sp>`` number density in the simulation's declared volume units
``temp:<sp>``    electron temperature Tₑ = ⅔⟨ε⟩ (eV); only meaningful for ``e-``
``excimer``      total of the network's lasing upper state (KrF*/ArF*); n/a
                 when the network declares no STIMULATED channel
``photons``      total real photons in the cavity mode
``p_in``         drive power deposited, W (last fluid-cadence record)
``p_out``        laser output power, W (last fluid-cadence record)
``efficiency``   cumulative intrinsic efficiency
``gain``         small-signal gain of the declared lasing upper state, 1/m
===============  ==============================================================

With ``spec.timing`` the terminal line carries trailing ``elapsed=``/``eta=`` wall-clock
columns. Heartbeat records hold the schema-1 dict (status/step/pct/wall/fields) with field
values as the same formatted strings the terminal shows — strict-JSON safe even when a
diverging run renders ``nan``/``inf``. Heartbeat write failures warn once and disable the
sink: progress reporting must never abort the physics.
"""

from __future__ import annotations

import json
import os
import time
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .components import ProgressReport

_SPECIES_FIELDS = ("total", "density", "temp")
_PLAIN_FIELDS = ("excimer", "photons", "p_in", "p_out", "efficiency", "gain")
_STYLES = ("percent", "steps")


def _last_power(series: list[float]) -> str:
    return f"{series[-1]:.3e}" if series else "n/a"


def format_duration(seconds: float) -> str:
    """``12.3s`` below one minute, ``4m05s`` below one hour, ``1h02m`` above."""
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    minutes, secs = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m{secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m"


def write_heartbeat(path: Path | str, record: dict) -> None:
    """Atomically replace ``path`` with ``record`` as JSON (same-dir tmp + ``os.replace``)."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(record), encoding="utf-8")
    os.replace(tmp, path)


def read_heartbeat(path: Path | str) -> dict | None:
    """Parsed heartbeat dict, or ``None`` when missing, partial, or not a JSON object."""
    try:
        record = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return record if isinstance(record, dict) else None


@dataclass(frozen=True)
class _Field:
    """One resolved progress field: display label plus value renderer."""

    label: str
    render: Callable[[], str]


class ProgressReporter:
    """Deliver progress reports for one built :class:`~Fluid1D.simulation.Simulation`."""

    def __init__(self, spec: ProgressReport, sim) -> None:
        if spec.style not in _STYLES:
            raise ValueError(f"ProgressReport.style must be one of {_STYLES}, got {spec.style!r}")
        if spec.every is not None and spec.every < 1:
            raise ValueError(f"ProgressReport.every must be >= 1 (or None), got {spec.every}")
        if not (spec.terminal or spec.heartbeat):
            raise ValueError(
                "ProgressReport reports nowhere: terminal and heartbeat are both disabled"
            )
        self.spec = spec
        self.sim = sim
        self.steps = sim.numerics.steps
        self.every = spec.every if spec.every is not None else max(self.steps // 10, 1)
        self._t0 = time.monotonic()
        self._fields: list[_Field] = [self._resolve(field) for field in spec.fields]
        self._heartbeat_path: Path | None = None

    # ------------------------------------------------------------ lifecycle
    def start(self, run_dir: Path | str) -> None:
        """Bind the heartbeat sink to ``run_dir`` and write the initial running record."""
        if self.spec.heartbeat:
            self._heartbeat_path = Path(run_dir) / "progress.json"
            self._write(self._record(0))

    def due(self, step: int) -> bool:
        # The final step always reports so every run closes on a 100% / N-of-N line.
        return step % self.every == 0 or step == self.steps

    def emit(self, step: int) -> None:
        """Deliver one due report to every enabled sink."""
        if self.spec.terminal:
            print(f"  {self.line(step)}", flush=True)
        if self._heartbeat_path is not None:
            self._write(self._record(step))

    def finalize(self, status: str = "done", error: str | None = None) -> None:
        """Write the closing heartbeat record (no-op without a bound heartbeat sink)."""
        if self._heartbeat_path is not None:
            self._write(self._record(self.steps, status=status, error=error))

    # ------------------------------------------------------------ rendering
    def line(self, step: int) -> str:
        cols = [self._progress_column(step)]
        cols += [f"{field.label}={field.render()}" for field in self._fields]
        if self.spec.timing:
            elapsed = time.monotonic() - self._t0
            cols.append(f"elapsed={format_duration(elapsed)}")
            cols.append(f"eta={format_duration(elapsed * (self.steps - step) / step)}")
        return "  ".join(cols)

    def _record(self, step: int, status: str = "running", error: str | None = None) -> dict:
        record = {
            "schema": 1,
            "status": status,
            "step": step,
            "steps": self.steps,
            "pct": round(100.0 * step / self.steps, 1),
            "t_sim_s": step * self.sim.numerics.dt,
            "wall_s": round(time.monotonic() - self._t0, 3),
            "updated_unix": time.time(),
            "fields": {field.label: field.render() for field in self._fields},
        }
        if error is not None:
            record["error"] = error
        return record

    def _write(self, record: dict) -> None:
        try:
            write_heartbeat(self._heartbeat_path, record)
        except OSError as exc:
            # Progress must never abort the physics: warn once, drop the sink.
            warnings.warn(f"progress heartbeat write failed ({exc}); heartbeat disabled")
            self._heartbeat_path = None

    # ------------------------------------------------------------ resolution
    def _resolve(self, field: str) -> _Field:
        sim = self.sim
        kind, _, species = field.partition(":")
        if species:
            if kind not in _SPECIES_FIELDS:
                raise ValueError(
                    f"unknown progress field {field!r}; species-prefixed fields are "
                    f"{[k + ':<species>' for k in _SPECIES_FIELDS]}"
                )
            if species not in sim.registry.names:
                raise ValueError(
                    f"progress field {field!r} names unknown species {species!r}; "
                    f"registered: {', '.join(sim.registry.names)}"
                )
            if kind == "total":
                return _Field(species, lambda: f"{sim.species_total(species):.2e}")
            if kind == "density":
                vol_m3 = sim.domain_length * sim.area_si  # SI domain volume
                return _Field(f"n({species})", lambda: (
                    f"{sim.units.density_from_si(sim.species_total(species) / vol_m3):.2e}"
                ))
            # temp: electron temperature only
            if species == "e-":
                return _Field("T(e-)", lambda: f"{sim.electron_temperature():.3g} eV")
            return _Field(f"T({species})", lambda: "n/a")

        if field == "excimer":
            lasing = sim.gain.lasing
            if lasing is None:
                return _Field("excimer", lambda: "n/a")
            upper = lasing[0]
            return _Field(upper, lambda: f"{sim.species_total(upper):.2e}")
        if field == "photons":
            return _Field("photons", lambda: f"{sim.gain.total_photons():.2e}")
        if field == "p_in":
            return _Field("P_in", lambda: f"{_last_power(sim.energetics.power_in_w)} W")
        if field == "p_out":
            return _Field("P_out", lambda: f"{_last_power(sim.energetics.power_out_w)} W")
        if field == "efficiency":
            return _Field("η", lambda: f"{sim.energetics.intrinsic_efficiency():.2e}")
        if field == "gain":
            return _Field("gain", lambda: f"{sim.gain_coefficient():.3e} /m")
        raise ValueError(
            f"unknown progress field {field!r}; plain fields are {_PLAIN_FIELDS}, "
            f"species fields are {[k + ':<species>' for k in _SPECIES_FIELDS]}"
        )

    def _progress_column(self, step: int) -> str:
        if self.spec.style == "steps":
            width = len(str(self.steps))
            return f"{step:{width}d}/{self.steps}"
        return f"{100 * step // self.steps:3d}%"
