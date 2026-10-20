"""
Threading-layer selection and Numba thread configuration.

Numba's parallel ``prange`` kernels run under a threading layer. ``workqueue``
(static chunking, no work-stealing) is always available but scales poorly on
asymmetric chips; ``tbb`` (work-stealing) and ``omp`` are preferred where their
runtimes are installed (typically Linux/x86 HPC). This module probes what
actually loads on the host, picks the best layer, and chooses a thread count from
the :class:`~backend.topology.CPUTopology` heuristic.

The threading layer must be selected *before* the first parallel kernel executes,
so :func:`configure_numba` is called once at startup (before grid warmup).
"""

from __future__ import annotations

import os

from .topology import CPUTopology

# Preference order: work-stealing schedulers first, static workqueue last.
_LAYER_PREFERENCE = ("tbb", "omp", "workqueue")


def available_threading_layers() -> list[str]:
    """Threading layers whose backing runtime is importable on this host.

    ``workqueue`` ships with Numba and is always available. ``tbb``/``omp`` are
    reported only if their pool extension imports (the runtime is present).
    """
    layers = ["workqueue"]
    try:
        import numba.np.ufunc.omppool  # noqa: F401

        layers.append("omp")
    except ImportError:
        pass
    try:
        import numba.np.ufunc.tbbpool  # noqa: F401

        layers.append("tbb")
    except ImportError:
        pass
    return layers


def select_threading_layer(available: list[str] | None = None) -> str:
    """Best available layer by preference (tbb > omp > workqueue)."""
    avail = set(available if available is not None else available_threading_layers())
    for layer in _LAYER_PREFERENCE:
        if layer in avail:
            return layer
    return "workqueue"


def optimal_thread_count(topology: CPUTopology, override: int | None = None) -> int:
    """Pick a Numba thread count for this host.

    Precedence: explicit ``override`` > ``NUMBA_NUM_THREADS`` env > heuristic.
    Heuristic: on an asymmetric chip use only the performance cores (efficiency
    cores straggle under static chunking and stall the reduction barrier);
    otherwise use physical cores. Always clamped to any HPC allocation limit and
    floored at 1.
    """
    if override is not None and override > 0:
        chosen = override
    else:
        env = os.environ.get("NUMBA_NUM_THREADS")
        if env and env.isdigit() and int(env) > 0:
            chosen = int(env)
        elif topology.is_heterogeneous:
            chosen = topology.performance_cores
        else:
            chosen = topology.physical_cores or topology.logical_cores

    if topology.allocation_limit is not None:
        chosen = min(chosen, topology.allocation_limit)
    return max(1, chosen)


def configure_numba(threads: int, layer: str | None = None) -> str:
    """Set the Numba threading layer and thread count for this process.

    Must run before the first parallel kernel executes. Returns the layer name
    actually configured. Also pins ``OPENBLAS_NUM_THREADS=1`` (BLAS is off the fluid
    hot path; this prevents a BLAS pool from oversubscribing against Numba).
    Numba being absent is not an error — the call becomes a no-op.
    """
    if os.environ.get("OPENBLAS_NUM_THREADS") is None:
        os.environ["OPENBLAS_NUM_THREADS"] = "1"

    try:
        import numba
    except ImportError:
        return "none"

    resolved = layer or select_threading_layer()
    # THREADING_LAYER must be set before the first parallel dispatch; setting both
    # the env var and the config object covers import-order edge cases.
    os.environ.setdefault("NUMBA_THREADING_LAYER", resolved)
    try:
        numba.config.THREADING_LAYER = resolved
    except Exception:  # pragma: no cover - defensive; never block startup
        pass
    numba.set_num_threads(threads)
    return resolved
