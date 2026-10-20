"""Data-aware, archive-backed static plots for fluid diagnostics.

The physics diagnostics only collect and export numeric arrays.  Run files compose the
figures from those arrays with :class:`LinePlot` and :class:`HeatmapPlot`, including every
output filename, x/y dataset selection, scale, and optional fixed range.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Protocol, runtime_checkable

import h5py
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

AxisKind = Literal["linear", "log", "symlog"]
ColorKind = Literal["linear", "log"]


@dataclass(frozen=True)
class DataRef:
    """One HDF5 dataset selected by a run-file plot declaration."""

    dataset: str
    label: str
    scale: float = 1.0

    def read(self, archive: h5py.File) -> NDArray[np.float64]:
        if self.dataset not in archive:
            raise KeyError(f"diagnostic dataset {self.dataset!r} is not in {archive.filename}")
        return np.asarray(archive[self.dataset], dtype=np.float64) * self.scale


@runtime_checkable
class Plot(Protocol):
    """A run-file-owned figure rendered from ``diagnostics.h5``."""

    def render(self, archive_path: Path, out_dir: Path) -> list[Path]: ...


def _finite_values(arrays: tuple[NDArray[np.float64], ...]) -> NDArray[np.float64]:
    finite = [values[np.isfinite(values)] for values in arrays if values.size]
    if not finite:
        return np.empty(0, dtype=np.float64)
    return np.concatenate(finite)


def data_limits(
    *arrays: NDArray[np.float64],
    scale: AxisKind = "linear",
    pad: float = 0.05,
    include_zero: bool | None = None,
) -> tuple[float, float]:
    """Return tight limits that retain all finite plottable data.

    Linear limits use the actual minimum and maximum rather than rounding to a decade.
    ``include_zero=None`` includes zero only when the data already approach or cross it;
    this keeps pulse baselines meaningful without wasting most of an axis on a remote zero.
    Log limits use multiplicative padding and ignore non-positive values, which Matplotlib
    cannot display on a logarithmic axis. An all-zero log series gets a harmless default
    range so an inactive diagnostic does not abort an otherwise valid simulation run.
    """
    if pad < 0.0:
        raise ValueError("plot padding must be non-negative")
    values = _finite_values(tuple(np.asarray(a, dtype=np.float64).ravel() for a in arrays))
    if scale == "log":
        values = values[values > 0.0]
        if values.size == 0:
            return (1.0, 10.0)
        lo = float(values.min())
        hi = float(values.max())
        if lo == hi:
            factor = 10.0 ** max(pad, 0.05)
            return (lo / factor, hi * factor)
        decades = np.log10(hi) - np.log10(lo)
        factor = 10.0 ** (max(decades, 1.0) * pad)
        return (lo / factor, hi * factor)
    if values.size == 0:
        return (0.0, 1.0)

    lo = float(values.min())
    hi = float(values.max())
    if lo == hi:
        if lo == 0.0:
            return (0.0, 1.0)
        delta = max(abs(lo) * max(pad, 0.05), np.finfo(np.float64).eps)
        return (lo - delta, hi + delta)

    if include_zero is None:
        include_zero = (
            (lo >= 0.0 and lo <= 0.05 * hi)
            or (hi <= 0.0 and abs(hi) <= 0.05 * abs(lo))
        )
    if include_zero:
        lo = min(lo, 0.0)
        hi = max(hi, 0.0)
    delta = (hi - lo) * pad
    lower = lo - delta
    upper = hi + delta
    if include_zero and lo == 0.0:
        lower = 0.0
    if include_zero and hi == 0.0:
        upper = 0.0
    return (lower, upper)


def _visible_values(
    x: NDArray[np.float64],
    values: tuple[NDArray[np.float64], ...],
    x_limits: tuple[float, float],
) -> tuple[NDArray[np.float64], ...]:
    visible: list[NDArray[np.float64]] = []
    mask = (x >= x_limits[0]) & (x <= x_limits[1])
    for y in values:
        if y.ndim != 1:
            raise ValueError("LinePlot datasets must be one-dimensional")
        if y.shape != x.shape:
            raise ValueError(f"LinePlot x/y shapes differ: {x.shape} != {y.shape}")
        visible.append(y[mask])
    return tuple(visible)


def _coordinate_edges(centers: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return cell edges that fully bracket pcolormesh coordinate centres."""
    if centers.size == 0:
        raise ValueError("a heatmap coordinate axis cannot be empty")
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=np.float64)
    if np.any(np.diff(centers) <= 0.0):
        raise ValueError("heatmap coordinates must be strictly increasing")
    middle = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - (middle[0] - centers[0])
    last = centers[-1] + (centers[-1] - middle[-1])
    return np.concatenate(([first], middle, [last]))


def _output_path(out_dir: Path, filename: str | Path) -> Path:
    path = out_dir / Path(filename)
    if path.suffix == "":
        raise ValueError(f"plot filename needs an extension: {filename!s}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


@dataclass(frozen=True)
class LinePlot:
    """One archive-backed line figure, optionally with a secondary y-axis."""

    filename: str | Path
    x: DataRef
    y: tuple[DataRef, ...]
    ylabel: str
    secondary_y: tuple[DataRef, ...] = ()
    secondary_ylabel: str | None = None
    title: str | None = None
    x_scale: AxisKind = "linear"
    y_scale: AxisKind = "linear"
    secondary_y_scale: AxisKind = "linear"
    x_limits: tuple[float, float] | None = None
    y_limits: tuple[float, float] | None = None
    secondary_y_limits: tuple[float, float] | None = None
    x_include_zero: bool | None = None
    y_include_zero: bool | None = None
    secondary_y_include_zero: bool | None = None
    x_padding: float = 0.02
    y_padding: float = 0.05
    markers: bool = False
    legend_columns: int = 2
    figure_size: tuple[float, float] = (7.0, 4.2)
    dpi: int = 150

    def __post_init__(self) -> None:
        if not self.y:
            raise ValueError("LinePlot needs at least one primary y dataset")
        if self.secondary_y and self.secondary_ylabel is None:
            raise ValueError("secondary_ylabel is required with secondary_y datasets")

    @staticmethod
    def _draw(
        ax: plt.Axes,
        x: NDArray[np.float64],
        refs: tuple[DataRef, ...],
        values: tuple[NDArray[np.float64], ...],
        *,
        markers: bool,
        color_offset: int = 0,
    ) -> None:
        styles = ("-", "--", "-.", ":")
        marker_kinds = ("o", "s", "^", "D", "v", "P", "X", "*")
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        for index, (ref, y) in enumerate(zip(refs, values)):
            kwargs: dict[str, object] = {
                "label": ref.label,
                "linestyle": styles[index % len(styles)],
                "linewidth": 1.5,
                "color": colors[(index + color_offset) % len(colors)],
            }
            if markers:
                kwargs.update(
                    marker=marker_kinds[index % len(marker_kinds)],
                    markersize=4,
                    markevery=max(1, len(y) // 20),
                )
            ax.plot(x, y, **kwargs)

    def render(self, archive_path: Path, out_dir: Path) -> list[Path]:
        with h5py.File(archive_path, "r") as archive:
            x = self.x.read(archive)
            primary = tuple(ref.read(archive) for ref in self.y)
            secondary = tuple(ref.read(archive) for ref in self.secondary_y)
        if x.ndim != 1:
            raise ValueError("LinePlot x dataset must be one-dimensional")

        fitted_x = self.x_limits or data_limits(
            x, scale=self.x_scale, pad=self.x_padding, include_zero=self.x_include_zero
        )
        primary_visible = _visible_values(x, primary, fitted_x)
        secondary_visible = _visible_values(x, secondary, fitted_x)

        fig, ax = plt.subplots(figsize=self.figure_size)
        self._draw(ax, x, self.y, primary, markers=self.markers)
        ax.set_xscale(self.x_scale)
        ax.set_yscale(self.y_scale)
        ax.set_xlim(*fitted_x)
        ax.set_ylim(*(self.y_limits or data_limits(
            *primary_visible,
            scale=self.y_scale,
            pad=self.y_padding,
            include_zero=self.y_include_zero,
        )))
        ax.set_xlabel(self.x.label)
        ax.set_ylabel(self.ylabel)
        axes = [ax]

        if secondary:
            ax2 = ax.twinx()
            self._draw(
                ax2, x, self.secondary_y, secondary,
                markers=self.markers, color_offset=len(primary),
            )
            ax2.set_yscale(self.secondary_y_scale)
            ax2.set_ylim(*(self.secondary_y_limits or data_limits(
                *secondary_visible,
                scale=self.secondary_y_scale,
                pad=self.y_padding,
                include_zero=self.secondary_y_include_zero,
            )))
            ax2.set_ylabel(self.secondary_ylabel)
            axes.append(ax2)
        if self.title:
            ax.set_title(self.title)
        lines = [line for axis in axes for line in axis.get_lines()]
        if lines:
            ax.legend(
                lines,
                [line.get_label() for line in lines],
                fontsize="small",
                ncol=self.legend_columns,
            )
        fig.tight_layout()
        path = _output_path(out_dir, self.filename)
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        return [path]


@dataclass(frozen=True)
class HeatmapPanel:
    """One z dataset and title inside a :class:`HeatmapPlot`."""

    z: DataRef
    title: str


@dataclass(frozen=True)
class HeatmapPlot:
    """One or more x-y heatmaps with independent data-aware color ranges."""

    filename: str | Path
    x: DataRef
    y: DataRef
    panels: tuple[HeatmapPanel, ...]
    color_scale: ColorKind = "linear"
    x_scale: AxisKind = "linear"
    y_scale: AxisKind = "linear"
    x_limits: tuple[float, float] | None = None
    y_limits: tuple[float, float] | None = None
    color_limits: tuple[float, float] | None = None
    cmap: str = "viridis"
    figure_height: float = 3.8
    dpi: int = 150

    def __post_init__(self) -> None:
        if not self.panels:
            raise ValueError("HeatmapPlot needs at least one panel")

    def render(self, archive_path: Path, out_dir: Path) -> list[Path]:
        with h5py.File(archive_path, "r") as archive:
            x = self.x.read(archive)
            y = self.y.read(archive)
            blocks = tuple(panel.z.read(archive) for panel in self.panels)
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("HeatmapPlot x and y datasets must be one-dimensional")
        x_edges = _coordinate_edges(x)
        y_edges = _coordinate_edges(y)

        fig, axes = plt.subplots(
            1, len(self.panels),
            figsize=(4.2 * len(self.panels), self.figure_height),
            squeeze=False,
        )
        for ax, panel, block in zip(axes[0], self.panels, blocks):
            if block.shape != (y.size, x.size):
                raise ValueError(
                    f"heatmap {panel.z.dataset!r} shape {block.shape} does not match "
                    f"y/x sizes {(y.size, x.size)}"
                )
            if self.color_limits is not None:
                vmin, vmax = self.color_limits
            else:
                color_kind: AxisKind = "log" if self.color_scale == "log" else "linear"
                vmin, vmax = data_limits(block, scale=color_kind, pad=0.0)
            norm: mcolors.Normalize
            if self.color_scale == "log":
                norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
                plotted = np.ma.masked_less_equal(block, 0.0)
            else:
                norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
                plotted = block
            mesh = ax.pcolormesh(
                x_edges, y_edges, plotted, cmap=self.cmap, norm=norm, shading="flat",
            )
            fig.colorbar(mesh, ax=ax, label=panel.z.label)
            ax.set_title(panel.title)
            ax.set_xlabel(self.x.label)
            ax.set_ylabel(self.y.label)
            ax.set_xscale(self.x_scale)
            ax.set_yscale(self.y_scale)
            ax.set_xlim(*(self.x_limits or data_limits(x_edges, scale=self.x_scale, pad=0.0)))
            ax.set_ylim(*(self.y_limits or data_limits(y_edges, scale=self.y_scale, pad=0.0)))
        fig.tight_layout()
        path = _output_path(out_dir, self.filename)
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        return [path]


def plot_time_series(
    times: NDArray[np.float64],
    series: dict[str, NDArray[np.float64]],
    *,
    ax: plt.Axes | None = None,
    xlabel: str = "time (ps)",
    ylabel: str = "value",
    log_y: bool = False,
    ylim: tuple[float, float] | None = None,
) -> plt.Axes:
    """Draw named in-memory series against a shared time axis."""
    if ax is None:
        _, ax = plt.subplots()
    refs = tuple(DataRef(label, label) for label in series)
    values = tuple(np.asarray(value, dtype=np.float64) for value in series.values())
    LinePlot._draw(ax, np.asarray(times), refs, values, markers=True)
    y_scale: AxisKind = "log" if log_y else "linear"
    ax.set_yscale(y_scale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif values:
        ax.set_ylim(*data_limits(*values, scale=y_scale))
    ax.legend(fontsize="small", ncol=2)
    return ax


def export_figure(path: str | Path, *, dpi: int = 150) -> Path:
    """Save and close the current Matplotlib figure."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out
