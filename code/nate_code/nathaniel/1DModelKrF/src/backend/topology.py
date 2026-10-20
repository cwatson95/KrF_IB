"""
CPU topology detection for the compute backend.

Detects logical/physical core counts, the performance vs. efficiency core split
on asymmetric chips (Apple Silicon, ARM big.LITTLE), the platform/architecture,
and any HPC scheduler allocation limit (SLURM/PBS). The thread-count heuristic in
:mod:`backend.threading` consumes this to pick an optimal Numba thread count per
machine, so the same code expands sensibly across a laptop and heterogeneous HPC
compute units.
"""

from __future__ import annotations

import os
import platform
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class CPUTopology:
    """Resolved CPU layout for one host.

    ``performance_cores`` / ``efficiency_cores`` are ``0`` when the split is
    unknown or the chip is homogeneous. ``allocation_limit`` is the core budget
    imposed by an HPC scheduler (SLURM/PBS) or CPU affinity, or ``None`` if
    unconstrained.
    """

    logical_cores: int
    physical_cores: int
    performance_cores: int
    efficiency_cores: int
    platform: str
    arch: str
    allocation_limit: int | None

    @property
    def is_heterogeneous(self) -> bool:
        """True when the chip has both performance and efficiency cores."""
        return self.performance_cores > 0 and self.efficiency_cores > 0

    def describe(self) -> str:
        parts = [f"{self.platform}/{self.arch}", f"{self.logical_cores} logical"]
        if self.physical_cores:
            parts.append(f"{self.physical_cores} physical")
        if self.is_heterogeneous:
            parts.append(f"{self.performance_cores}P+{self.efficiency_cores}E")
        if self.allocation_limit is not None:
            parts.append(f"alloc<={self.allocation_limit}")
        return ", ".join(parts)


def _int_or_zero(value: str | None) -> int:
    try:
        return int(value) if value is not None else 0
    except (TypeError, ValueError):
        return 0


def _sysctl(*keys: str) -> dict[str, int]:
    """Read integer sysctl keys on macOS; missing keys are omitted."""
    out: dict[str, int] = {}
    for key in keys:
        try:
            res = subprocess.run(
                ["sysctl", "-n", key],
                capture_output=True,
                text=True,
                timeout=2.0,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if res.returncode == 0:
            val = _int_or_zero(res.stdout.strip())
            if val > 0:
                out[key] = val
    return out


def _allocation_limit(logical: int) -> int | None:
    """Core budget from an HPC scheduler or CPU affinity, if smaller than ``logical``."""
    candidates: list[int] = []
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NP", "NSLOTS"):
        candidates.append(_int_or_zero(os.environ.get(var)))
    # CPU affinity mask (cgroups / taskset) caps usable cores on Linux.
    if hasattr(os, "sched_getaffinity"):
        try:
            candidates.append(len(os.sched_getaffinity(0)))
        except OSError:
            pass
    limits = [c for c in candidates if c > 0 and c < logical]
    return min(limits) if limits else None


def _macos_cores() -> tuple[int, int, int]:
    """(physical, performance, efficiency) cores on macOS; zeros where unknown."""
    keys = _sysctl(
        "hw.physicalcpu",
        "hw.perflevel0.logicalcpu",  # performance cores
        "hw.perflevel1.logicalcpu",  # efficiency cores
    )
    return (
        keys.get("hw.physicalcpu", 0),
        keys.get("hw.perflevel0.logicalcpu", 0),
        keys.get("hw.perflevel1.logicalcpu", 0),
    )


def _linux_cores(logical: int) -> tuple[int, int, int]:
    """(physical, performance, efficiency) on Linux via sysfs; best-effort."""
    physical = 0
    try:
        # Count distinct (physical_id, core_id) pairs from /proc/cpuinfo.
        with open("/proc/cpuinfo", encoding="ascii", errors="ignore") as fh:
            pkgs: set[tuple[str, str]] = set()
            phys = core = None
            for line in fh:
                if line.startswith("physical id"):
                    phys = line.split(":")[-1].strip()
                elif line.startswith("core id"):
                    core = line.split(":")[-1].strip()
                elif line.strip() == "" and phys is not None and core is not None:
                    pkgs.add((phys, core))
                    phys = core = None
            physical = len(pkgs)
    except OSError:
        physical = 0

    # big.LITTLE detection: distinct cpu_capacity values across CPUs.
    perf = eff = 0
    try:
        caps: dict[int, int] = {}
        base = "/sys/devices/system/cpu"
        for entry in os.listdir(base):
            if not entry.startswith("cpu") or not entry[3:].isdigit():
                continue
            cap_path = os.path.join(base, entry, "cpu_capacity")
            try:
                with open(cap_path, encoding="ascii") as fh:
                    caps[int(entry[3:])] = int(fh.read().strip())
            except (OSError, ValueError):
                continue
        if caps and len(set(caps.values())) > 1:
            top = max(caps.values())
            perf = sum(1 for v in caps.values() if v == top)
            eff = len(caps) - perf
    except OSError:
        pass

    return (physical or logical, perf, eff)


def detect_topology() -> CPUTopology:
    """Probe the running host's CPU layout (best-effort, never raises)."""
    system = platform.system().lower()  # darwin / linux / windows
    arch = platform.machine().lower()  # arm64 / x86_64 / ...
    logical = os.cpu_count() or 1

    physical = perf = eff = 0
    if system == "darwin":
        physical, perf, eff = _macos_cores()
    elif system == "linux":
        physical, perf, eff = _linux_cores(logical)
    else:  # windows / other: only logical count is portable.
        physical = logical

    return CPUTopology(
        logical_cores=logical,
        physical_cores=physical or logical,
        performance_cores=perf,
        efficiency_cores=eff,
        platform=system,
        arch=arch,
        allocation_limit=_allocation_limit(logical),
    )
