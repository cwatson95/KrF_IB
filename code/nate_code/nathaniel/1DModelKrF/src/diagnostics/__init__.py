"""Diagnostics toolkit for the 1D fluid engine: HDF5 archive + static figures.

Numeric collectors live in :mod:`Fluid1D.diagnostics`; run files use the archive-backed plot
declarations here to select filenames, axes, datasets, and scales.
"""

from .archive import DiagnosticsArchive
from .plot import (
    DataRef,
    HeatmapPanel,
    HeatmapPlot,
    LinePlot,
    Plot,
    data_limits,
    export_figure,
    plot_time_series,
)

__all__ = (
    "DataRef",
    "DiagnosticsArchive",
    "HeatmapPanel",
    "HeatmapPlot",
    "LinePlot",
    "Plot",
    "data_limits",
    "export_figure",
    "plot_time_series",
)
