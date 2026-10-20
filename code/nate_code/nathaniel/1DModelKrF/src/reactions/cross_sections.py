"""Load Cross-Section Browser JSON+CSV exports and build interpolation tables."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import RegularGridInterpolator

from reactions.network import CompiledNetwork, Reaction

from .types import ReactionKind, reaction_kind_from_metadata

SCHEMA_VERSION = "cross-section-metadata/1.1"
DataKind = Literal["total", "differential"]
Metadata = dict[str, object]


def _parse_metadata_line(line: str) -> tuple[str, str] | None:
    stripped = line.strip()
    if not stripped.startswith("#"):
        return None
    body = stripped.lstrip("#").strip()
    if ":" not in body:
        return None
    key, value = body.split(":", 1)
    return key.strip().lower(), value.strip()


def _resolve_paired_paths(path: Path) -> tuple[Path, Path]:
    """Return ``(json_path, csv_path)`` for a Cross-Section Browser export pair.

    The JSON's ``data_file`` names the CSV; when that file is missing but a CSV
    with the JSON's own stem exists (a renamed export whose metadata went
    stale), the same-stem CSV is used instead.
    """
    path = Path(path)
    if path.suffix.lower() == ".json":
        json_path = path
        meta = json.loads(json_path.read_text(encoding="utf-8"))
        csv_name = str(meta.get("data_file", ""))
        if not csv_name:
            raise ValueError(f"JSON cross-section metadata missing data_file: {json_path}")
        csv_path = json_path.parent / csv_name
        if not csv_path.is_file() and json_path.with_suffix(".csv").is_file():
            csv_path = json_path.with_suffix(".csv")
    elif path.suffix.lower() == ".csv":
        csv_path = path
        json_path = csv_path.with_suffix(".json")
    else:
        raise ValueError(f"unsupported cross-section path extension: {path.suffix!r}")

    if not json_path.is_file():
        raise FileNotFoundError(f"missing JSON metadata for cross-section: {json_path}")
    if not csv_path.is_file():
        raise FileNotFoundError(f"missing CSV data for cross-section: {csv_path}")
    return json_path, csv_path


def _load_csv_table(csv_path: Path) -> tuple[Metadata, list[str], NDArray[np.float64]]:
    """Parse comment metadata, header row, and numeric data from a cross-section CSV."""
    comment_meta: Metadata = {}
    header: list[str] = []
    rows: list[list[float]] = []

    with csv_path.open(newline="", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                parsed = _parse_metadata_line(line)
                if parsed is not None:
                    comment_meta[parsed[0]] = parsed[1]
                continue
            header = [part.strip() for part in line.strip().split(",")]
            break

        if not header:
            raise ValueError(f"cross-section CSV has no header row: {csv_path}")

        reader = csv.reader(handle)
        for row in reader:
            if not row or not row[0].strip():
                continue
            rows.append([float(value) for value in row])

    if not rows:
        raise ValueError(f"cross-section CSV contains no data rows: {csv_path}")

    data = np.asarray(rows, dtype=np.float64)
    return comment_meta, header, data


def _apply_axis_conversion(values: NDArray[np.float64], factor: float) -> NDArray[np.float64]:
    return np.asarray(values, dtype=np.float64) * float(factor)


# CSV column unit -> SI multiplier (energies in eV, cross sections in m^2).
_CSV_UNIT_TO_SI: dict[str, float] = {
    "": 1.0,  # unit not declared -> assume SI
    "ev": 1.0,
    "m^2": 1.0,
    "m2": 1.0,
    "cm^2": 1e-4,
    "cm2": 1e-4,
}


def _csv_unit_to_si(unit: str, json_path: Path) -> float:
    """SI multiplier for a CSV column unit; unknown units fail loudly."""
    key = unit.strip().lower()
    factor = _CSV_UNIT_TO_SI.get(key)
    if factor is None:
        raise ValueError(
            f"unknown CSV column unit {unit!r} in {json_path}; "
            f"known: {sorted(k for k in _CSV_UNIT_TO_SI if k)}"
        )
    return factor


def interpolate_sigma_linear(
    energy_ev: float | NDArray[np.float64],
    energy_grid_ev: NDArray[np.float64],
    sigma_m2: NDArray[np.float64],
    *,
    threshold_ev: float = 0.0,
) -> float | NDArray[np.float64]:
    """Linearly interpolate σ(E) between the two nearest tabulated energies (m²)."""
    scalar = np.isscalar(energy_ev)
    e_arr = np.atleast_1d(np.asarray(energy_ev, dtype=np.float64))
    grid = np.asarray(energy_grid_ev, dtype=np.float64)
    sigma = np.asarray(sigma_m2, dtype=np.float64)
    out = np.zeros_like(e_arr, dtype=np.float64)

    if grid.size == 0:
        return float(out[0]) if scalar else out

    above = e_arr >= threshold_ev
    if not np.any(above):
        return float(out[0]) if scalar else out

    e_use = e_arr[above]
    if grid.size == 1:
        out[above] = sigma[0]
        return float(out[0]) if scalar else out

    idx = np.searchsorted(grid, e_use, side="right") - 1
    idx = np.clip(idx, 0, grid.size - 2)
    e0 = grid[idx]
    e1 = grid[idx + 1]
    s0 = sigma[idx]
    s1 = sigma[idx + 1]
    denom = e1 - e0
    # Endpoint clamp, never linear extrapolation. Extrapolation made a falling final segment
    # negative and a rising one grow without evidence beyond the measured energy range.
    frac = np.clip(np.where(denom > 0.0, (e_use - e0) / denom, 0.0), 0.0, 1.0)
    out[above] = s0 + frac * (s1 - s0)
    return float(out[0]) if scalar else out


@dataclass
class CrossSectionRecord:
    """One Cross-Section Browser process (total or differential)."""

    process_id: str
    data_kind: DataKind
    reaction: str
    reaction_kind: ReactionKind
    projectile: str
    target: str
    threshold_ev: float | None
    metadata: Metadata
    energy_ev: NDArray[np.float64]
    sigma_m2: NDArray[np.float64]
    secondary_axis: NDArray[np.float64] | None = None
    secondary_quantity: str | None = None
    source_path: Path | None = None

    @property
    def is_total(self) -> bool:
        return self.data_kind == "total"

    @property
    def is_differential(self) -> bool:
        return self.data_kind == "differential"


@dataclass
class TotalCrossSectionInterpolator:
    """Linear bracketing interpolator for a total cross section σ(E)."""

    record: CrossSectionRecord

    @classmethod
    def from_record(cls, record: CrossSectionRecord) -> TotalCrossSectionInterpolator:
        if not record.is_total:
            raise ValueError("TotalCrossSectionInterpolator requires data_kind='total'")
        return cls(record=record)

    def __call__(self, energy_ev: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
        return interpolate_sigma_linear(
            energy_ev,
            self.record.energy_ev,
            self.record.sigma_m2,
            threshold_ev=self.record.threshold_ev or 0.0,
        )


_LUT_POINTS_DEFAULT = 1024


@dataclass
class TabulatedSigma:
    """σ(ε) resampled onto a uniform log-energy grid for O(1) evaluation.

    One log, one multiply-add index, and one linear blend per lookup provide
    constant-time evaluation with no per-call binary search.
    Coefficient builders evaluate the same table arithmetic. Sigma is zero at or below
    ``threshold_ev``; below the tabulated range sigma
    clamps to the first value and above it to the last (no extrapolation).
    """

    log_e0: float
    inv_dlog: float
    values: NDArray[np.float64]
    threshold_ev: float
    e_min_ev: float
    e_max_ev: float
    record: CrossSectionRecord = field(repr=False)

    @classmethod
    def from_record(
        cls, record: CrossSectionRecord, *, n_points: int = _LUT_POINTS_DEFAULT
    ) -> TabulatedSigma:
        if not record.is_total:
            raise ValueError(f"tabulated σ requires data_kind='total', got {record.data_kind!r}")
        energy = np.asarray(record.energy_ev, dtype=np.float64)
        if energy.size == 0 or not np.all(np.isfinite(energy)) or energy.max() <= 0.0:
            raise ValueError("cross-section table has no usable energy grid")
        threshold = max(float(record.threshold_ev or 0.0), 0.0)
        e_min = float(max(energy.min(), 1e-4))
        e_max = float(max(energy.max(), e_min * (1.0 + 1e-9)))
        grid = np.logspace(math.log10(e_min), math.log10(e_max), n_points)
        vals = np.asarray(
            interpolate_sigma_linear(grid, energy, record.sigma_m2, threshold_ev=threshold),
            dtype=np.float64,
        )
        np.maximum(vals, 0.0, out=vals)
        log_e0 = math.log(e_min)
        dlog = (math.log(e_max) - log_e0) / (n_points - 1)
        return cls(
            log_e0=log_e0,
            inv_dlog=1.0 / dlog,
            values=vals,
            threshold_ev=threshold,
            e_min_ev=e_min,
            e_max_ev=e_max,
            record=record,
        )

    def __call__(self, energy_ev: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
        """Evaluate σ (m²) for scalar or array ε (eV); mirrors the numba kernel."""
        scalar = np.isscalar(energy_ev)
        e_arr = np.atleast_1d(np.asarray(energy_ev, dtype=np.float64))
        out = np.zeros_like(e_arr)
        above = e_arr > self.threshold_ev if self.threshold_ev > 0.0 else e_arr > 0.0
        if np.any(above):
            x = (np.log(e_arr[above]) - self.log_e0) * self.inv_dlog
            idx = np.clip(np.floor(x).astype(np.int64), 0, self.values.size - 2)
            frac = np.clip(x - idx, 0.0, 1.0)
            v = self.values
            out[above] = v[idx] + frac * (v[idx + 1] - v[idx])
        return float(out[0]) if scalar else out


@dataclass
class DifferentialCrossSectionInterpolator:
    """2D log-log interpolator for differential cross sections."""

    record: CrossSectionRecord
    _interp: RegularGridInterpolator = field(repr=False)

    @classmethod
    def from_record(cls, record: CrossSectionRecord) -> DifferentialCrossSectionInterpolator:
        if not record.is_differential:
            raise ValueError("DifferentialCrossSectionInterpolator requires data_kind='differential'")
        if record.secondary_axis is None:
            raise ValueError("differential record missing secondary_axis")

        energy = record.energy_ev
        secondary = record.secondary_axis
        table = record.sigma_m2
        if table.ndim != 2:
            raise ValueError(f"expected 2D sigma table, got shape {table.shape}")

        log_e = np.log10(np.maximum(energy, 1e-30))
        log_s = np.log10(np.maximum(secondary, 1e-30))
        safe = np.maximum(table, 1e-40)
        log_table = np.log10(safe)

        interp = RegularGridInterpolator(
            (log_e, log_s),
            log_table,
            bounds_error=False,
            fill_value=np.nan,
        )
        return cls(record=record, _interp=interp)

    def __call__(self, energy_ev: float, secondary_value: float) -> float:
        log_e = math.log10(max(energy_ev, 1e-30))
        log_s = math.log10(max(secondary_value, 1e-30))
        log_val = float(self._interp((log_e, log_s)))
        if not math.isfinite(log_val):
            return 0.0
        return math.pow(10.0, log_val)

def load_cross_section(path: str | Path) -> CrossSectionRecord:
    """Load a Cross-Section Browser JSON+CSV pair."""
    json_path, csv_path = _resolve_paired_paths(Path(path))
    meta: Metadata = json.loads(json_path.read_text(encoding="utf-8"))

    schema = str(meta.get("schema", ""))
    if schema != SCHEMA_VERSION:
        raise ValueError(f"unsupported cross-section schema {schema!r} in {json_path}")

    comment_meta, header, data = _load_csv_table(csv_path)
    data_kind = str(meta.get("data_kind", comment_meta.get("data_kind", "total"))).lower()
    if data_kind not in ("total", "differential"):
        raise ValueError(f"unsupported data_kind {data_kind!r}")

    # CSV -> SI factors come from the `units` dict describing the CSV columns.
    # `conversion_from_source` documents the source->CSV provenance and must NOT
    # be applied here (doing so re-scaled cm^2 columns by 1e4 instead of 1e-4).
    units = meta.get("units", {})
    units = units if isinstance(units, dict) else {}
    x_factor = _csv_unit_to_si(str(units.get("x", "")), json_path)
    y_factor = _csv_unit_to_si(str(units.get("y", "")), json_path)

    axes = meta.get("axes", {})
    x_quantity = str(axes.get("x_quantity", "Energy")) if isinstance(axes, dict) else "Energy"
    y_quantity = str(axes.get("y_quantity", "Cross section")) if isinstance(axes, dict) else "Cross section"

    threshold_raw = meta.get("threshold", comment_meta.get("threshold"))
    threshold_ev: float | None
    if threshold_raw is None:
        threshold_ev = None
    else:
        threshold_ev = float(str(threshold_raw).split()[0])

    type_name = str(meta.get("type", ""))
    category = str(meta.get("category", ""))
    reaction_kind = reaction_kind_from_metadata(type_name, category)

    if data_kind == "total":
        if data.shape[1] < 2:
            raise ValueError(f"total cross-section CSV needs 2 columns: {csv_path}")
        energy_ev = _apply_axis_conversion(data[:, 0], x_factor)
        sigma_m2 = _apply_axis_conversion(data[:, 1], y_factor)
        if (not np.all(np.isfinite(energy_ev)) or not np.all(np.isfinite(sigma_m2))
                or np.any(energy_ev < 0.0) or np.any(sigma_m2 < 0.0)):
            raise ValueError(f"total cross-section data must be finite and non-negative: {csv_path}")
        order = np.argsort(energy_ev, kind="stable")
        energy_ev = energy_ev[order]
        sigma_m2 = sigma_m2[order]
        unique, inverse, counts = np.unique(energy_ev, return_inverse=True, return_counts=True)
        if unique.size != energy_ev.size:
            summed = np.zeros(unique.size, dtype=np.float64)
            np.add.at(summed, inverse, sigma_m2)
            energy_ev = unique
            sigma_m2 = summed / counts
        secondary_axis = None
        secondary_quantity = None
    else:
        # differential exports: `units.y` describes the SECONDARY AXIS column
        # (e.g. secondary energy in eV), not the tabulated values, which stay
        # in SI as exported.
        if data.shape[1] < 3:
            raise ValueError(f"differential cross-section CSV needs >=3 columns: {csv_path}")
        energy_ev = np.unique(_apply_axis_conversion(data[:, 0], x_factor))
        secondary_axis = np.unique(_apply_axis_conversion(data[:, 1], y_factor))
        n_e = energy_ev.size
        n_s = secondary_axis.size
        table = np.zeros((n_e, n_s), dtype=np.float64)
        e_to_i = {float(v): i for i, v in enumerate(energy_ev)}
        s_to_i = {float(v): i for i, v in enumerate(secondary_axis)}
        for row in data:
            e_val = float(row[0]) * x_factor
            s_val = float(row[1]) * y_factor
            ei = e_to_i[e_val]
            si = s_to_i[s_val]
            table[ei, si] = float(row[2])
        if (not np.all(np.isfinite(energy_ev)) or not np.all(np.isfinite(secondary_axis))
                or not np.all(np.isfinite(table)) or np.any(energy_ev < 0.0)
                or np.any(table < 0.0)):
            raise ValueError(f"differential cross-section data must be finite and non-negative: {csv_path}")
        sigma_m2 = table
        secondary_quantity = y_quantity if y_quantity.lower() != "cross section" else header[1]

    return CrossSectionRecord(
        process_id=str(meta.get("process_id", json_path.stem)),
        data_kind=data_kind,  # type: ignore[arg-type]
        reaction=str(meta.get("reaction", comment_meta.get("reaction", ""))),
        reaction_kind=reaction_kind,
        projectile=str(meta.get("projectile", "")),
        target=str(meta.get("target", "")),
        threshold_ev=threshold_ev,
        metadata=meta,
        energy_ev=energy_ev,
        sigma_m2=sigma_m2,
        secondary_axis=secondary_axis,
        secondary_quantity=secondary_quantity,
        source_path=json_path,
    )


@dataclass
class CrossSectionCatalog:
    """Index of Cross-Section Browser exports for one projectile–target pair."""

    projectile: str
    target: str
    records: list[CrossSectionRecord] = field(default_factory=list)
    _total_interpolators: dict[str, TotalCrossSectionInterpolator] = field(
        default_factory=dict, repr=False
    )
    def add(self, record: CrossSectionRecord) -> None:
        self.records.append(record)
        key = record.process_id
        if record.is_total:
            self._total_interpolators[key] = TotalCrossSectionInterpolator.from_record(record)

    def total_cross_section(
        self,
        energy_ev: float | NDArray[np.float64],
    ) -> float | NDArray[np.float64]:
        """Sum enabled total cross sections at ``energy_ev`` (m²)."""
        scalar = np.isscalar(energy_ev)
        e_arr = np.atleast_1d(np.asarray(energy_ev, dtype=np.float64))
        total = np.zeros_like(e_arr, dtype=np.float64)
        for interp in self._total_interpolators.values():
            total += np.asarray(interp(e_arr), dtype=np.float64)
        return float(total[0]) if scalar else total

def load_cross_section_directory(
    directory: str | Path,
    *,
    projectile: str | None = None,
    target: str | None = None,
) -> CrossSectionCatalog:
    """Load all ``*.json`` exports in ``directory`` into a catalog."""
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"cross-section directory not found: {directory}")

    records: list[CrossSectionRecord] = []
    for json_path in sorted(directory.glob("*.json")):
        record = load_cross_section(json_path)
        if projectile is not None and record.projectile != projectile:
            continue
        if target is not None and record.target != target:
            continue
        records.append(record)

    if not records:
        raise ValueError(f"no cross-section records found in {directory}")

    cat = CrossSectionCatalog(
        projectile=records[0].projectile,
        target=records[0].target,
    )
    for record in records:
        cat.add(record)
    return cat


def build_cross_section_index(base_dir: Path) -> dict[str, Path]:
    """Map ``process_id`` from JSON metadata to each export's ``.json`` path."""
    base_dir = Path(base_dir)
    index: dict[str, Path] = {}
    if not base_dir.is_dir():
        return index
    for json_path in sorted(base_dir.glob("*.json")):
        meta = json.loads(json_path.read_text(encoding="utf-8"))
        process_id = str(meta.get("process_id", json_path.stem))
        index[process_id] = json_path
    return index


def resolve_cross_section_path(key: str, base_dir: Path) -> Path:
    """Resolve a TOML ``cross_section`` key to a JSON+CSV export path.

    ``key`` is normally the export's basename (``"e_Kr_ionization"`` →
    ``base_dir/e_Kr_ionization.json``, paired to its CSV via ``data_file``).
    A ``.json``/``.csv`` filename (looked up inside ``base_dir``), an explicit
    path, or a ``process_id`` declared in any JSON metadata under ``base_dir``
    also resolve.
    """
    base_dir = Path(base_dir)
    key_path = Path(key)
    candidates: list[Path] = []
    if key_path.suffix.lower() in (".json", ".csv"):
        candidates.append(base_dir / key_path.name)
        candidates.append(key_path)
    else:
        candidates.append(base_dir / f"{key}.json")
        candidates.append(base_dir / f"{key}.csv")
    for cand in candidates:
        if cand.is_file():
            return cand

    indexed = build_cross_section_index(base_dir).get(key)
    if indexed is not None and indexed.is_file():
        return indexed

    searched = [str(c) for c in candidates] + [f"{base_dir}/*.json (process_id={key!r})"]
    raise FileNotFoundError(
        f"cross-section {key!r} not found; searched: {', '.join(searched)}"
    )


_SPECIES_ALIASES = {"e": "e-", "E": "e-"}


def _normalize_species(name: str) -> str:
    name = name.strip()
    return _SPECIES_ALIASES.get(name, name)


def _expected_reaction_kind(rxn: Reaction) -> ReactionKind:
    """Reaction kind implied by TOML electron stoichiometry."""
    e_in = sum(c for c, n in rxn.reactants if n == "e-")
    e_out = sum(c for c, n in rxn.products if n == "e-")
    if e_out > e_in:
        return ReactionKind.IONIZATION
    if e_out < e_in:
        return ReactionKind.ATTACHMENT
    threshold = rxn.params.get("threshold_eV") or 0.0
    return ReactionKind.EXCITATION if threshold > 0.0 else ReactionKind.ELASTIC


_COMPATIBLE_KINDS: dict[ReactionKind, frozenset[ReactionKind]] = {
    ReactionKind.IONIZATION: frozenset({ReactionKind.IONIZATION}),
    ReactionKind.ATTACHMENT: frozenset({ReactionKind.ATTACHMENT, ReactionKind.RECOMBINATION}),
    ReactionKind.EXCITATION: frozenset({ReactionKind.EXCITATION, ReactionKind.VIBRATIONAL}),
    ReactionKind.ELASTIC: frozenset({ReactionKind.ELASTIC}),
}


def match_record_to_reaction(record: CrossSectionRecord, rxn: Reaction) -> None:
    """Raise ``ValueError`` when a dataset's metadata does not describe ``rxn``.

    Checks the export's projectile/target species and reaction kind against the
    TOML reaction, so a mislinked ``cross_section`` key is reported instead of silently
    colliding with the wrong physics. The dataset's own reaction string is included in
    every error for reference.
    """
    reactant_names = {n for _, n in rxn.reactants}
    projectile = _normalize_species(record.projectile)
    target = _normalize_species(record.target)
    if projectile and projectile not in reactant_names:
        raise ValueError(
            f"dataset projectile {record.projectile!r} is not a reactant of "
            f"{rxn.describe()!r} (dataset reaction: {record.reaction!r})"
        )
    if target and target not in reactant_names:
        raise ValueError(
            f"dataset target {record.target!r} is not a reactant of "
            f"{rxn.describe()!r} (dataset reaction: {record.reaction!r})"
        )
    expected = _expected_reaction_kind(rxn)
    allowed = _COMPATIBLE_KINDS.get(expected, frozenset({expected}))
    if record.reaction_kind not in allowed:
        raise ValueError(
            f"dataset kind {record.reaction_kind.value!r} does not match "
            f"{expected.value!r} implied by {rxn.describe()!r} "
            f"(dataset reaction: {record.reaction!r})"
        )


@dataclass(frozen=True)
class CrossSectionImport:
    """Outcome of importing one reaction's ``cross_section`` dataset."""

    reaction_id: str
    reaction: str  # TOML reaction, rxn.describe()
    key: str  # the TOML cross_section value
    status: Literal["imported", "failed"]
    detail: str
    sigma: TabulatedSigma | None = None

    @property
    def ok(self) -> bool:
        return self.status == "imported"

    def describe(self) -> str:
        return f"cross-section [{self.key}] {self.detail}"


@dataclass
class CrossSectionImportReport:
    """Per-reaction ``cross_section`` import outcomes for one compiled network."""

    imports: list[CrossSectionImport] = field(default_factory=list)

    @property
    def n_failed(self) -> int:
        return sum(1 for imp in self.imports if not imp.ok)

    def lines(self) -> list[str]:
        return [imp.describe() for imp in self.imports]


def import_cross_sections(
    compiled: CompiledNetwork,
    base_dir: str | Path,
    *,
    n_points: int = _LUT_POINTS_DEFAULT,
) -> CrossSectionImportReport:
    """Import every kinetic reaction's ``cross_section`` dataset; never raises.

    Success builds a :class:`TabulatedSigma` lookup. Any unresolvable key,
    unreadable file, invalid dataset, or reaction-metadata mismatch is recorded
    as a failed validation. Coefficient solvers load datasets directly and
    assign zero to missing electron-impact channels. One import result per
    reaction is retained for the simulation startup report.
    """
    report = CrossSectionImportReport()
    for rxn in compiled.kinetic:
        if not rxn.cross_section:
            continue
        key = str(rxn.cross_section)
        rid = rxn.reaction_id or rxn.describe()
        sigma: TabulatedSigma | None = None
        try:
            if not rxn.reaction_id:
                raise ValueError("reaction has no id= (overrides are keyed by reaction id)")
            path = resolve_cross_section_path(key, Path(base_dir))
            record = load_cross_section(path)
            if not record.is_total:
                raise ValueError(f"data_kind {record.data_kind!r} is not 'total'")
            match_record_to_reaction(record, rxn)
            sigma = TabulatedSigma.from_record(record, n_points=n_points)
        except Exception as exc:  # noqa: BLE001 — report every validation failure
            status = "failed"
            detail = f"validation FAILED for reaction {rid} ({exc})"
        else:
            status = "imported"
            toml_thr = float(rxn.params.get("threshold_eV") or 0.0)
            thr_note = ""
            if (record.threshold_ev is not None and toml_thr > 0.0
                    and abs(record.threshold_ev - toml_thr) > max(0.25 * toml_thr, 0.5)):
                thr_note = (f"; NOTE dataset threshold {record.threshold_ev:g} eV "
                            f"vs TOML {toml_thr:g} eV")
            detail = (
                f"successfully imported for reaction {rid} ({rxn.describe()}): "
                f"{record.reaction or record.process_id}, {record.energy_ev.size} points, "
                f"{record.energy_ev.min():g}-{record.energy_ev.max():g} eV{thr_note}"
            )
        report.imports.append(CrossSectionImport(
            reaction_id=rxn.reaction_id, reaction=rxn.describe(), key=key,
            status=status, detail=detail, sigma=sigma,
        ))
    return report
