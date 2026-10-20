"""A single HDF5 archive of the numeric data behind the rendered diagnostics.

Numeric collectors keep light time-series (gain, densities, efficiency, EEDF curve) in
memory and write them here, one group per diagnostic. Run-file plot declarations render
selected datasets after the archive closes, so figures are reproducible without re-running.

``DiagnosticsArchive`` is a thin, dependency-light wrapper over :mod:`h5py`: each
diagnostic's :meth:`export` hook is handed the archive and writes its own group via
:meth:`series` (a time axis + named line series) or :meth:`array` (a bare dataset).
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType

import h5py
import numpy as np
from numpy.typing import NDArray


class DiagnosticsArchive:
    """Write the numerics behind every in-memory diagnostic to one HDF5 file."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._file = h5py.File(self.path, "w")

    def group(self, name: str) -> h5py.Group:
        """Create (or fetch) a named group for a diagnostic to write into."""
        return self._file.require_group(name)

    def series(
        self,
        name: str,
        times: NDArray[np.float64],
        columns: dict[str, NDArray[np.float64]],
        *,
        time_key: str = "times_ps",
    ) -> None:
        """Write a time axis + one dataset per named line under group ``name``."""
        grp = self.group(name)
        grp.create_dataset(time_key, data=np.asarray(times, dtype=np.float64))
        for key, values in columns.items():
            grp.create_dataset(key, data=np.asarray(values, dtype=np.float64))

    def array(self, name: str, key: str, data: NDArray[np.float64]) -> None:
        """Write a single named dataset under group ``name``."""
        self.group(name).create_dataset(key, data=np.asarray(data, dtype=np.float64))

    def close(self) -> None:
        if self._file.id.valid:
            self._file.close()

    def __enter__(self) -> DiagnosticsArchive:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()
