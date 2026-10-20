"""
Compute backend: runtime CPU-topology detection and thread configuration.

One entry point, :func:`configure`, detects the host's cores, selects a Numba threading
layer, picks a thread count, and exposes the result via :func:`get_config`. The simulation
configures this backend before warming the transport and chemistry kernels.

Typical use (once, at startup, *before* any grid is built / kernels warm up)::

    import backend
    cfg = backend.configure(threads=args.threads)   # None -> heuristic
    ...
    workers = backend.get_config().fft_workers

Backend selection is explicit and process-global rather than threaded through every
simulation component.
"""

from __future__ import annotations

from dataclasses import dataclass

from .threading import (
    available_threading_layers,
    configure_numba,
    optimal_thread_count,
    select_threading_layer,
)
from .topology import CPUTopology, detect_topology

__all__ = (
    "CPUTopology",
    "RuntimeConfig",
    "available_threading_layers",
    "configure",
    "detect_topology",
    "get_config",
    "optimal_thread_count",
    "select_threading_layer",
)


@dataclass(frozen=True)
class RuntimeConfig:
    """Resolved compute configuration for this process."""

    n_threads: int
    fft_workers: int
    layer: str
    topology: CPUTopology

    def describe(self) -> str:
        return (
            f"backend: {self.n_threads} threads, layer={self.layer}, "
            f"fft_workers={self.fft_workers} [{self.topology.describe()}]"
        )


_CONFIG: RuntimeConfig | None = None


def configure(
    threads: int | None = None,
    layer: str | None = None,
    *,
    force: bool = False,
) -> RuntimeConfig:
    """Detect topology, select threading layer + thread count, apply to Numba.

    Idempotent: the first call wins (the Numba threading layer cannot change once
    a parallel kernel has run). Pass ``force=True`` only in tests. ``threads``
    overrides the heuristic; ``layer`` overrides layer auto-selection.
    """
    global _CONFIG
    if _CONFIG is not None and not force:
        return _CONFIG

    topology = detect_topology()
    n_threads = optimal_thread_count(topology, override=threads)
    resolved_layer = configure_numba(n_threads, layer=layer)
    _CONFIG = RuntimeConfig(
        n_threads=n_threads,
        fft_workers=n_threads,
        layer=resolved_layer,
        topology=topology,
    )
    return _CONFIG


def get_config() -> RuntimeConfig:
    """Return the active config, lazily configuring with defaults if needed."""
    if _CONFIG is None:
        return configure()
    return _CONFIG
