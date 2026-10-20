"""Heartbeat-based sweep status table (``python -m Fluid1D.sweep_status``).

Reads the per-point ``progress.json`` heartbeats that parallel-sweep workers rewrite on
their progress cadence (see :mod:`Fluid1D.progress`) and renders a plain-text table — one
row per ``f2_*`` point directory. Run it from a second terminal against a live sweep::

    python -m Fluid1D.sweep_status runs/krf_f2_sweep --watch 5

``running`` rows whose heartbeat has not been rewritten within ``--stale`` seconds are
flagged ``stalled`` (dead or interrupted workers never write a final status). The reader
is invocation-agnostic: it shows whatever is on disk, including heartbeats left by an
earlier sweep of the same name.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .progress import format_duration, read_heartbeat

_COLUMNS = ("f2", "status", "%", "step", "t_sim", "wall", "updated", "fields")


@dataclass(frozen=True)
class PointStatus:
    """One sweep point's heartbeat state (``record`` is None when progress.json is absent)."""

    f2: float
    path: Path
    record: dict | None


def scan_sweep(sweep_dir: Path | str) -> list[PointStatus]:
    """All ``f2_*`` point dirs under ``sweep_dir`` with their parsed heartbeats, by f2."""
    statuses: list[PointStatus] = []
    for path in Path(sweep_dir).glob("f2_*"):
        if not path.is_dir():
            continue
        try:
            f2 = float(path.name[3:])
        except ValueError:
            continue  # non-numeric f2_* sibling (never written by the sweep itself)
        statuses.append(PointStatus(f2, path, read_heartbeat(path / "progress.json")))
    statuses.sort(key=lambda status: status.f2)
    return statuses


def _row_status(record: dict | None, stale_after_s: float, now_unix: float) -> str:
    if record is None or record.get("schema") != 1:
        return "-"
    status = str(record.get("status", "-"))
    if status == "running":
        age = now_unix - float(record.get("updated_unix", 0.0))
        if age > stale_after_s:
            return "stalled"
    return status


def _row(status: PointStatus, label: str, now_unix: float) -> tuple[str, ...]:
    record = status.record
    if record is None or record.get("schema") != 1:
        return (f"{status.f2:.6g}", label, "-", "-", "-", "-", "-", "-")
    step, steps = record.get("step", 0), record.get("steps", 0)
    t_sim = float(record.get("t_sim_s", 0.0))
    age = now_unix - float(record.get("updated_unix", now_unix))
    if label == "failed":
        error = str(record.get("error", ""))
        fields = error.strip().splitlines()[-1] if error.strip() else "-"
    else:
        fields = "  ".join(f"{k}={v}" for k, v in record.get("fields", {}).items()) or "-"
    return (
        f"{status.f2:.6g}",
        label,
        f"{float(record.get('pct', 0.0)):.0f}%",
        f"{step}/{steps}",
        f"{t_sim * 1e9:.3g} ns",
        format_duration(float(record.get("wall_s", 0.0))),
        f"{format_duration(max(age, 0.0))} ago",
        fields,
    )


def render_table(statuses: Sequence[PointStatus], *, stale_after_s: float = 120.0,
                 now_unix: float | None = None) -> str:
    """Fixed-width status table plus a one-line footer of per-status counts."""
    now = time.time() if now_unix is None else now_unix
    labels = [_row_status(status.record, stale_after_s, now) for status in statuses]
    rows = [_row(status, label, now) for status, label in zip(statuses, labels)]
    widths = [max(len(_COLUMNS[i]), *(len(row[i]) for row in rows)) if rows else
              len(_COLUMNS[i]) for i in range(len(_COLUMNS))]
    lines = ["  ".join(name.ljust(widths[i]) for i, name in enumerate(_COLUMNS)).rstrip()]
    lines += ["  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)).rstrip()
              for row in rows]
    counts = {label: labels.count(label)
              for label in ("done", "running", "stalled", "failed", "-")}
    footer = "  ".join(f"{name if name != '-' else 'no-data'} {count}"
                       for name, count in counts.items() if count)
    lines.append(f"[{len(statuses)} points]  {footer}" if footer
                 else f"[{len(statuses)} points]")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m Fluid1D.sweep_status",
        description="Render the per-point heartbeat status of a parallel F2 sweep.",
    )
    parser.add_argument("sweep_dir", help="sweep output directory (runs/<name>)")
    parser.add_argument("--watch", type=float, default=None, metavar="SECONDS",
                        help="re-render every SECONDS until interrupted")
    parser.add_argument("--stale", type=float, default=120.0, metavar="SECONDS",
                        help="flag 'running' heartbeats older than this as stalled"
                             " (default: 120)")
    args = parser.parse_args(argv)
    sweep_dir = Path(args.sweep_dir)
    if not sweep_dir.is_dir():
        print(f"sweep_status: no such directory: {sweep_dir}", file=sys.stderr)
        return 2
    if args.watch is None:
        statuses = scan_sweep(sweep_dir)
        if not statuses:
            print(f"sweep_status: no f2_* point dirs under {sweep_dir}", file=sys.stderr)
            return 1
        print(render_table(statuses, stale_after_s=args.stale))
        return 0
    try:
        while True:  # pragma: no branch - loop exits only via KeyboardInterrupt
            statuses = scan_sweep(sweep_dir)
            print(f"--- {time.strftime('%H:%M:%S')}  {sweep_dir}", flush=True)
            if statuses:
                print(render_table(statuses, stale_after_s=args.stale), flush=True)
            else:
                print("(no f2_* point dirs yet)", flush=True)
            time.sleep(args.watch)
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
