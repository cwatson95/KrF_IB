"""Numeric diagnostics for the 1D fluid simulation — listed in the experiment file.

Each diagnostic is a small collector with two hooks the engine calls:

- ``record(sim, step, t)`` — every step; the diagnostic decides (from its own cadence)
  whether to sample the live simulation state;
- ``export(sim, archive)`` — after the run; write the diagnostic's raw numerics into the
  shared ``diagnostics.h5``.

Collectors do not choose figures or filenames. Experiment files compose archive-backed
``LinePlot``/``HeatmapPlot`` declarations so the run author owns every output name, x/y
dataset, scale, and optional fixed range.

All quantities are derived from the **1D fluid representation** — per-cell species densities,
the electron mean-energy field, scalar cavity photon state, and energetics ledger. Outputs are
numeric arrays in ``diagnostics.h5``; run-owned declarations can turn selected arrays into
static 1D/time-series plots. There are deliberately **no animations and no 3D/ParaView
(XDMF) volumetric exports**.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from particles.constants import AVOGADRO, TD


def _copy_float_array(values: object) -> NDArray[np.float64]:
    """Return diagnostic values as an independent one-dimensional float array."""
    return np.array(values, dtype=np.float64, copy=True)


def _named_series(
    series: dict[str, list[float]],
    name: str,
) -> NDArray[np.float64]:
    """Return one copied series with a useful error for an unknown name."""
    try:
        values = series[name]
    except KeyError:
        available = ", ".join(sorted(series)) or "none"
        raise KeyError(f"unknown diagnostic series {name!r}; available: {available}") from None
    return _copy_float_array(values)


@runtime_checkable
class Diagnostic(Protocol):
    def record(self, sim: object, step: int, t: float) -> None: ...
    def export(self, sim: object, archive: object) -> None: ...


@dataclass
class EnergyClosure:
    """Electron-energy closure terms, residuals, and limiter activity versus time.

    Energy columns are interval-integrated joules.  ``continuous_reconciled`` archives are
    explicitly marked diagnostic because their accepted chemistry cost remains an interval
    reconciliation rather than the production source stage.
    """

    every: int = 1
    _t: list[float] = field(default_factory=list, repr=False)
    _series: dict[str, list[float]] = field(default_factory=dict, repr=False)
    _mode: str = field(default="", repr=False)
    _totals: dict[str, float] = field(default_factory=dict, repr=False)

    @property
    def times_ps(self) -> NDArray[np.float64]:
        """Sample times in picoseconds as an independent array."""
        return _copy_float_array(self._t)

    @property
    def series_names(self) -> tuple[str, ...]:
        """Names accepted by :meth:`series`."""
        return tuple(sorted(self._series))

    def series(self, name: str) -> NDArray[np.float64]:
        """Return a copied time series by its archived dataset name."""
        return _named_series(self._series, name)

    @property
    def totals(self) -> dict[str, float]:
        """Run-integrated energy-closure terms in joules as a copied mapping."""
        return dict(self._totals)

    @property
    def chemistry_energy_cost_j(self) -> float:
        """Accepted inelastic reaction cost accumulated over the run (J)."""
        return float(self._totals.get("accepted_reaction_cost_j", 0.0))

    @property
    def chemistry_bound_adjustment_j(self) -> float:
        """Coefficient-table-bound electron-energy adjustment over the run (J)."""
        return float(self._totals.get("coefficient_bound_adjustment_j", 0.0))

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        snapshot = sim.energy_closure_snapshot()  # type: ignore[attr-defined]
        self._mode = str(snapshot.pop("closure_mode"))
        self._t.append(t * 1e12)
        for name, value in snapshot.items():
            self._series.setdefault(name, []).append(float(value))

    def export(self, sim: object, archive: object) -> None:
        self._totals = {
            str(name): float(value)
            for name, value in sim._energy_closure_totals.items()  # type: ignore[attr-defined]
        }
        if not self._t:
            return
        archive.series(  # type: ignore[attr-defined]
            "energy_closure", np.asarray(self._t),
            {name: np.asarray(values) for name, values in self._series.items()},
        )
        group = archive.group("energy_closure")  # type: ignore[attr-defined]
        group.attrs["closure_mode"] = self._mode
        group.attrs["production_closure"] = self._mode == "source_coupled"
        if self._mode == "continuous_reconciled":
            group.attrs["status"] = (
                "diagnostic only; continuous cooling reconciled to accepted interval chemistry cost"
            )
        for name, value in self._totals.items():
            group.attrs[f"total_{name}"] = float(value)
        group.attrs["total_limiter_activity"] = float(sum(
            self._series.get("limiter_activity", ()),
        ))
        for name in (
            "equilibrium_te_residual_ev",
            "maximum_source_step_energy_jump_ev",
            "max_energy_closure_residual_j_m3",
            "max_ampere_residual",
        ):
            group.attrs[f"maximum_{name}"] = float(max(self._series.get(name, (0.0,))))


def _due(every: int | None, sim: object, step: int) -> bool:
    """True when a fluid-cadence diagnostic should sample at ``step``."""
    cadence = every if every is not None else sim.numerics.fluid_every  # type: ignore[attr-defined]
    return step > 0 and step % cadence == 0


@dataclass
class EEDF:
    """Electron energy distribution at one step: a smooth (Maxwellian) f(E) at Tₑ.

    The runtime fluid represents electrons by a per-cell mean energy, not a resolved EEDF field.
    This diagnostic exports the analytic Maxwellian EEDF at the domain-mean electron
    temperature as a *qualitative* f(E).
    """

    at_step: int | None = None  # default: final step
    energy_range: tuple[float, float] = (0.1, 100.0)
    bins: int = 60
    _te_ev: float | None = field(default=None, repr=False)

    def record(self, sim: object, step: int, t: float) -> None:
        target = self.at_step if self.at_step is not None else sim.numerics.steps  # type: ignore[attr-defined]
        if step == target:
            self._te_ev = float(sim.electron_temperature())  # type: ignore[attr-defined]

    def _curve(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        lo, hi = self.energy_range
        energy = np.linspace(lo, hi, self.bins)
        te = max(self._te_ev or 0.0, 1e-6)
        # Maxwellian EEDF normalized so ∫ f(E) dE = 1.
        f = 2.0 * np.sqrt(energy / np.pi) * te**-1.5 * np.exp(-energy / te)
        return energy, f

    def export(self, sim: object, archive: object) -> None:
        if self._te_ev is None:
            return
        energy, f = self._curve()
        grp = archive.group("eedf")  # type: ignore[attr-defined]
        grp.create_dataset("energy_ev", data=energy)
        grp.create_dataset("f_e", data=f)
        grp.attrs["te_ev"] = self._te_ev


@dataclass
class SpeciesDensities:
    """Number density of each tracked species vs time, in the experiment's volume units.

    ``species_total`` returns a *count* of real particles in the simulated volume; this divides
    by the (SI) domain volume and reports the number density in ``sim.units`` (cm⁻³ or m⁻³).
    """

    track: tuple[str, ...] = ()
    every: int | None = None
    _t: list[float] = field(default_factory=list, repr=False)
    _series: dict[str, list[float]] = field(default_factory=dict, repr=False)

    @property
    def times_ps(self) -> NDArray[np.float64]:
        """Sample times in picoseconds as an independent array."""
        return _copy_float_array(self._t)

    @property
    def series_names(self) -> tuple[str, ...]:
        """Recorded species names accepted by :meth:`series`."""
        return tuple(sorted(self._series))

    def series(self, name: str) -> NDArray[np.float64]:
        """Return one recorded species-density series as an independent array."""
        return _named_series(self._series, name)

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        vol_m3 = sim.domain_length * sim.area_si  # type: ignore[attr-defined]  # SI domain volume
        self._t.append(t * 1e12)
        for name in self.track:
            n_si = sim.species_total(name) / vol_m3  # type: ignore[attr-defined]  # m^-3
            self._series.setdefault(name, []).append(sim.units.density_from_si(n_si))  # type: ignore[attr-defined]

    def export(self, sim: object, archive: object) -> None:
        if not self._t:
            return
        archive.series("species_densities", np.asarray(self._t),  # type: ignore[attr-defined]
                       {k: np.asarray(v) for k, v in self._series.items()})
        archive.group("species_densities").attrs["units"] = sim.units.density_unit  # type: ignore[attr-defined]


@dataclass
class OpticalGain:
    """Volume-averaged small-signal gain of the declared lasing channel vs time."""

    every: int | None = None
    _t: list[float] = field(default_factory=list, repr=False)
    _g: list[float] = field(default_factory=list, repr=False)

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        self._t.append(t * 1e12)
        self._g.append(sim.gain_coefficient())  # type: ignore[attr-defined]

    def export(self, sim: object, archive: object) -> None:
        if not self._t:
            return
        u = sim.units  # type: ignore[attr-defined]
        archive.series("optical_gain", np.asarray(self._t),  # type: ignore[attr-defined]
                       {f"gain_per_{u.length_label}": u.inv_length_from_si(np.asarray(self._g))})


@dataclass
class Absorption:
    """Cavity-mode absorption coefficient α = Σ σ_i·⟨n_i⟩ over user channels.

    Channels resolve against the neutral fluid (e.g. F₂, Kr*) or the charged fluid (e.g. F⁻).
    """

    channels: dict[str, float] = field(default_factory=dict)  # species -> σ_abs (m^2)
    every: int | None = None
    _t: list[float] = field(default_factory=list, repr=False)
    _a: list[float] = field(default_factory=list, repr=False)

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        fluid = sim.fluid  # type: ignore[attr-defined]
        charged = sim.charged  # type: ignore[attr-defined]
        alpha = 0.0
        for name, sigma in self.channels.items():
            if name in fluid.column:
                alpha += sigma * float(fluid.density[:, fluid.column[name]].mean())
            elif name in charged.column:
                alpha += sigma * float(charged.density_of(name).mean())
        self._t.append(t * 1e12)
        self._a.append(alpha)

    def export(self, sim: object, archive: object) -> None:
        if not self._t:
            return
        u = sim.units  # type: ignore[attr-defined]
        archive.series("absorption", np.asarray(self._t),  # type: ignore[attr-defined]
                       {f"alpha_per_{u.length_label}": u.inv_length_from_si(np.asarray(self._a))})


def _cumulative_trapezoid(y: NDArray[np.float64], t: NDArray[np.float64]) -> NDArray[np.float64]:
    """Cumulative ∫y dt (trapezoid, starts at 0) — matches LaserEnergetics' integrals."""
    if t.size < 2:
        return np.zeros_like(t)
    increments = 0.5 * (y[1:] + y[:-1]) * np.diff(t)
    return np.concatenate([[0.0], np.cumsum(increments)])


@dataclass
class Efficiency:
    """Cumulative intrinsic efficiency ∫P_out / (E_init + ∫P_in) from the energetics ledger.

    The **pulse-integrated** scalar ``LaserEnergetics.intrinsic_efficiency()`` uses signed
    ``∫P_in`` (net drive work). The **cumulative curve** uses ``∫max(P_in, 0)`` in the
    denominator so energy returned to the field mid-pulse cannot drive the ratio through zero.
    """

    def record(self, sim: object, step: int, t: float) -> None:
        return  # power bookkeeping is filled by the engine into sim.energetics

    @staticmethod
    def _cumulative(e: object) -> NDArray[np.float64]:
        intervals = getattr(e, "interval_s", [])
        interval_mode = (len(intervals) == len(e.times_s)  # type: ignore[attr-defined]
                         and all(dt is not None for dt in intervals))
        if interval_mode:
            dt = np.asarray(intervals, dtype=np.float64)
            cum_out = np.cumsum(np.asarray(e.power_out_w) * dt)  # type: ignore[attr-defined]
        else:
            cum_out = _cumulative_trapezoid(np.asarray(e.power_out_w), e.t)  # type: ignore[attr-defined]
        # Signed P_in is net drive work — it can go negative when the plasma returns
        # energy to the field, driving E_init + ∫P_in through zero and blowing up η(t).
        # For the cumulative curve use one-sided deposition: energy put *into* the plasma.
        p_in_pos = np.maximum(np.asarray(e.power_in_w, dtype=np.float64), 0.0)  # type: ignore[attr-defined]
        if interval_mode:
            cum_in = np.cumsum(p_in_pos * dt)
        else:
            cum_in = _cumulative_trapezoid(p_in_pos, e.t)  # type: ignore[attr-defined]
        deposited = e.initial_energy_j + cum_in
        return cum_out / np.maximum(deposited, 1e-300)

    def export(self, sim: object, archive: object) -> None:
        e = sim.energetics  # type: ignore[attr-defined]
        if e.t.size == 0:
            return
        archive.series("efficiency", e.t * 1e12, {  # type: ignore[attr-defined]
            "power_out_w": np.asarray(e.power_out_w),
            "power_out_spont_w": np.asarray(e.power_out_spont_w),
            "power_out_stim_w": np.asarray(e.power_out_stim_w),
            "power_in_w": np.asarray(e.power_in_w),
            "cumulative_efficiency": self._cumulative(e),
        })
        archive.group("efficiency").attrs["initial_energy_j"] = e.initial_energy_j  # type: ignore[attr-defined]


@dataclass
class LaserOutput:
    """Output-coupler power split by photon origin: spontaneous, stimulated, total.

    The photon layer tags every photon by its creation event; the engine tallies the coupler
    flux per origin into ``sim.energetics`` on the fluid cadence. Below threshold the output
    is spontaneous-origin fluorescence; lasing shows up as the stimulated fraction → 1.
    """

    def record(self, sim: object, step: int, t: float) -> None:
        return  # power bookkeeping is filled by the engine into sim.energetics

    def export(self, sim: object, archive: object) -> None:
        e = sim.energetics  # type: ignore[attr-defined]
        if e.t.size == 0:
            return
        archive.series("laser_output", e.t * 1e12, {  # type: ignore[attr-defined]
            "power_out_spont_w": np.asarray(e.power_out_spont_w),
            "power_out_stim_w": np.asarray(e.power_out_stim_w),
            "power_out_w": np.asarray(e.power_out_w),
        })


@dataclass
class FieldLine:
    """Final-step 1D axial profiles: |E|(x), net charge ρ(x), and species densities n_s(x).

    The diagnostic records the arrays once at the final step.
    """

    species: tuple[str, ...] = ("KrF*", "e-")
    _x: NDArray[np.float64] | None = field(default=None, repr=False)
    _emag: NDArray[np.float64] | None = field(default=None, repr=False)
    _rho: NDArray[np.float64] | None = field(default=None, repr=False)
    _dens: dict[str, NDArray[np.float64]] = field(default_factory=dict, repr=False)

    def record(self, sim: object, step: int, t: float) -> None:
        if step != sim.numerics.steps:  # type: ignore[attr-defined]
            return
        # Stored in SI internally; converted to the experiment's volume units at export.
        self._x = np.asarray(sim.cell_centers, dtype=np.float64)  # type: ignore[attr-defined]  # m
        self._emag = np.linalg.norm(sim._last_efield, axis=1)  # type: ignore[attr-defined]  # V/m
        self._rho = sim.charged.net_charge_density().copy()  # type: ignore[attr-defined]  # C/m^3
        for name in self.species:
            fluid = sim.fluid  # type: ignore[attr-defined]
            charged = sim.charged  # type: ignore[attr-defined]
            if name in fluid.column:
                self._dens[name] = fluid.density[:, fluid.column[name]].copy()  # m^-3
            elif name in charged.column:
                self._dens[name] = charged.density_of(name).copy()  # m^-3

    def export(self, sim: object, archive: object) -> None:
        if self._x is None:
            return
        u = sim.units  # type: ignore[attr-defined]
        # Unit-in-name convention: densities in the volume unit, field V/m, charge C/m^3.
        cols: dict[str, NDArray[np.float64]] = {"e_mag_v_m": self._emag, "rho_c_m3": self._rho}
        for name, dens in self._dens.items():
            cols[f"n_{name}_per_{u.length_label}3"] = u.density_from_si(dens)
        archive.series(  # type: ignore[attr-defined]
            "field_line", u.length_from_si(self._x), cols,
            time_key=f"x_{u.length_label}",
        )


@dataclass
class ChargedDensityXT:
    """Spatio-temporal charged-species densities ``n_s(x, t)`` — one snapshot per fluid step.

    The temporal evolution of the spatial charge profiles: for each tracked charged species a
    ``(n_times, n_cells)`` block is accumulated (every ``every`` steps, default the fluid
    cadence) and written to ``diagnostics.h5`` under ``charged_density_xt`` with the time (ps)
    and position (volume-unit) axes. Run files may select any of them for x–t heatmaps. This
    is the export the real drift-diffusion transport makes meaningful (the stub held these flat).
    """

    species: tuple[str, ...] = ("e-", "Kr+", "Kr2+", "F-")
    every: int | None = None
    _x: NDArray[np.float64] | None = field(default=None, repr=False)
    _t: list[float] = field(default_factory=list, repr=False)
    _rows: dict[str, list[NDArray[np.float64]]] = field(default_factory=dict, repr=False)

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        charged = sim.charged  # type: ignore[attr-defined]
        if self._x is None:
            self._x = np.asarray(sim.cell_centers, dtype=np.float64)  # type: ignore[attr-defined]  # m
            self._rows = {name: [] for name in self.species if name in charged.column}
        self._t.append(t * 1.0e12)  # ps
        for name in self._rows:
            self._rows[name].append(charged.density_of(name).copy())  # m^-3 (view -> copy!)

    def _stacked(self) -> dict[str, NDArray[np.float64]]:
        return {name: np.asarray(rows, dtype=np.float64) for name, rows in self._rows.items() if rows}

    def export(self, sim: object, archive: object) -> None:
        if not self._t or self._x is None:
            return
        u = sim.units  # type: ignore[attr-defined]
        grp = archive.group("charged_density_xt")  # type: ignore[attr-defined]
        grp.create_dataset("times_ps", data=np.asarray(self._t, dtype=np.float64))
        grp.create_dataset(f"x_{u.length_label}", data=np.asarray(u.length_from_si(self._x), dtype=np.float64))
        for name, block in self._stacked().items():
            grp.create_dataset(f"n_{name}", data=np.asarray(u.density_from_si(block), dtype=np.float64))
        grp.attrs["density_units"] = u.density_unit
        grp.attrs["length_units"] = u.length_label


@dataclass
class PlasmaScalarTraces:
    """Per-step scalar traces feeding the sweep-overlay figures (Chen Figs. 7/12 style).

    Records mean electron density n_e [m^-3], mean reduced field E/N [Td], mean electron
    temperature T_e [eV], average discharge power density P_dis = I_cond*U_gap/V [W/m^3],
    and laser power density P_L = P_out/V [W/m^3] (V = gap x discharge area). Exports the
    series to ``plasma_scalars`` in diagnostics.h5 plus the per-run energy scalars
    (``e_out_j``, ``e_dep_j``, ``e_bank0_j``, ``e_stored0_j``) the FWHM/CE sweep plot consumes.
    """

    every: int = 1
    en_statistic: str = "mean"  # "mean" (Chen Fig. 7 bulk scale) or "max" (sheath-inclusive)
    _t: list[float] = field(default_factory=list, repr=False)
    _ne: list[float] = field(default_factory=list, repr=False)
    _en: list[float] = field(default_factory=list, repr=False)
    _te: list[float] = field(default_factory=list, repr=False)
    _pdis: list[float] = field(default_factory=list, repr=False)
    _pl: list[float] = field(default_factory=list, repr=False)
    _e0: float | None = field(default=None, repr=False)

    @property
    def times_ns(self) -> NDArray[np.float64]:
        """Sample times in nanoseconds as an independent array."""
        return _copy_float_array(self._t)

    @property
    def series_names(self) -> tuple[str, ...]:
        """Names accepted by :meth:`series`."""
        return ("en_td", "n_e_m3", "p_dis_wm3", "p_l_wm3", "te_ev")

    def series(self, name: str) -> NDArray[np.float64]:
        """Return a copied scalar series by its archived dataset name."""
        values = {
            "n_e_m3": self._ne,
            "en_td": self._en,
            "te_ev": self._te,
            "p_dis_wm3": self._pdis,
            "p_l_wm3": self._pl,
        }
        try:
            return _copy_float_array(values[name])
        except KeyError:
            available = ", ".join(self.series_names)
            raise KeyError(
                f"unknown diagnostic series {name!r}; available: {available}"
            ) from None

    def record(self, sim: object, step: int, t: float) -> None:
        if step % self.every:
            return
        vol = sim.domain_length * sim.area_si  # type: ignore[attr-defined]
        self._t.append(t * 1e9)
        self._ne.append(float(sim.charged.electron_density().mean()))  # type: ignore[attr-defined]
        n_gas = float(sim.fluid.density.sum(axis=1).mean())  # type: ignore[attr-defined]
        e_ax = np.abs(np.asarray(sim._last_efield)[:, 0])  # type: ignore[attr-defined]
        e_stat = float(e_ax.max() if self.en_statistic == "max" else e_ax.mean())
        self._en.append(e_stat / max(n_gas, 1.0) / TD)
        self._te.append(float(sim.charged.electron_temperature_ev().mean()))  # type: ignore[attr-defined]
        circ = sim.circuit  # type: ignore[attr-defined]
        if self._e0 is None:
            se = getattr(circ, "stored_energy_j", None)
            self._e0 = float(se()) if callable(se) else 0.0
        i_cond = float(getattr(sim, "_conduction_current", 0.0))
        v_b = float(getattr(circ, "v_b", 0.0) or 0.0)
        self._pdis.append(abs(i_cond * v_b) / vol)
        e = sim.energetics  # type: ignore[attr-defined]
        p_out = float(e.power_out_w[-1]) if len(e.power_out_w) else 0.0
        self._pl.append(p_out / vol)

    def export(self, sim: object, archive: object) -> None:
        if not self._t:
            return
        t = np.asarray(self._t)
        archive.series("plasma_scalars", t, {  # type: ignore[attr-defined]
            "n_e_m3": np.asarray(self._ne), "en_td": np.asarray(self._en),
            "te_ev": np.asarray(self._te), "p_dis_wm3": np.asarray(self._pdis),
            "p_l_wm3": np.asarray(self._pl),
        }, time_key="times_ns")
        grp = archive.group("plasma_scalars")  # type: ignore[attr-defined]
        grp.attrs["times_unit"] = "ns"
        vol = sim.domain_length * sim.area_si  # type: ignore[attr-defined]
        grp.attrs["volume_m3"] = vol
        grp.attrs["e_out_j"] = float(sim.energetics.output_energy_j())  # type: ignore[attr-defined]
        grp.attrs["e_dep_j"] = float(np.trapezoid(np.asarray(self._pdis) * vol, t * 1e-9))
        grp.attrs["e_stored0_j"] = float(self._e0 or 0.0)
        circuit_spec = sim.circuit_spec  # type: ignore[attr-defined]
        grp.attrs["e_bank0_j"] = float(0.5 * circuit_spec.c0_f * circuit_spec.v0_v**2)


@dataclass
class CircuitWaveforms:
    """Collect U_gap, U_C0 (kV), and I_d (kA) versus time for MPC circuits."""

    every: int = 1
    _t: list[float] = field(default_factory=list, repr=False)
    _ugap: list[float] = field(default_factory=list, repr=False)
    _uc0: list[float] = field(default_factory=list, repr=False)
    _id: list[float] = field(default_factory=list, repr=False)

    @property
    def times_ns(self) -> NDArray[np.float64]:
        """Sample times in nanoseconds as an independent array."""
        return _copy_float_array(self._t)

    @property
    def series_names(self) -> tuple[str, ...]:
        """Names accepted by :meth:`series`."""
        return ("i_d_ka", "u_c0_kv", "u_gap_kv")

    def series(self, name: str) -> NDArray[np.float64]:
        """Return a copied waveform by its archived dataset name."""
        values = {"u_gap_kv": self._ugap, "u_c0_kv": self._uc0, "i_d_ka": self._id}
        try:
            return _copy_float_array(values[name])
        except KeyError:
            available = ", ".join(self.series_names)
            raise KeyError(
                f"unknown diagnostic series {name!r}; available: {available}"
            ) from None

    @property
    def current_a(self) -> NDArray[np.float64]:
        """Discharge-current samples in amperes as an independent array."""
        return 1.0e3 * self.series("i_d_ka")

    def record(self, sim: object, step: int, t: float) -> None:
        circ = sim.circuit  # type: ignore[attr-defined]
        if step % self.every or not hasattr(circ, "v_c0"):
            return  # stub circuit: nothing to record
        self._t.append(t * 1e9)
        self._ugap.append(circ.v_b / 1e3)
        self._uc0.append(circ.v_c0 / 1e3)
        self._id.append(circ.i_s / 1e3)

    def export(self, sim: object, archive: object) -> None:
        if not self._t:
            return
        archive.series("circuit_waveforms", np.asarray(self._t), {  # type: ignore[attr-defined]
            "u_gap_kv": np.asarray(self._ugap), "u_c0_kv": np.asarray(self._uc0),
            "i_d_ka": np.asarray(self._id),
        }, time_key="times_ns")
        archive.group("circuit_waveforms").attrs["times_unit"] = "ns"  # type: ignore[attr-defined]


@dataclass
class ExcimerRates:
    """Temporal evolution of the excimer's formation and quenching reaction rates.

    Chen Fig. 7 style, generalized to any excimer (``ArF*`` / ``KrF*``): every sampled step,
    the volume-mean rate of progress of each reaction with a net stoichiometric source
    (formation) or sink (quenching) of ``excimer`` is recorded, in **mol m^-3 s^-1**
    (particle rates / N_A). Formation and quenching series are exported separately.

    The two radiative sinks are owned by the optics layer (excluded from the chemistry
    network), so they are sampled from the cavity's per-advance ledger instead: ``SE``
    (spontaneous) and ``STE`` (stimulated) join the quenching series as interval-averaged
    rates over the last optics advance, on the same domain-mean volumetric basis.
    """

    excimer: str = "ArF*"
    every: int | None = None  # default: the chemistry (fluid_every) cadence
    _t: list[float] = field(default_factory=list, repr=False)
    _rows: list[np.ndarray] = field(default_factory=list, repr=False)
    _se: list[float] = field(default_factory=list, repr=False)
    _ste: list[float] = field(default_factory=list, repr=False)
    _ids: tuple[str, ...] = field(default=(), repr=False)
    _nu_x: np.ndarray | None = field(default=None, repr=False)

    @property
    def times_ns(self) -> NDArray[np.float64]:
        """Sample times in nanoseconds as an independent array."""
        return _copy_float_array(self._t)

    def formation_rates(self) -> dict[str, NDArray[np.float64]]:
        """Formation rates in mol m^-3 s^-1 as an independent mapping of arrays."""
        formation, _ = self._split()
        return {name: _copy_float_array(values) for name, values in formation.items()}

    def quenching_rates(self) -> dict[str, NDArray[np.float64]]:
        """Quenching rates in mol m^-3 s^-1 as an independent mapping of arrays."""
        _, quenching = self._split()
        return {name: _copy_float_array(values) for name, values in quenching.items()}

    def record(self, sim: object, step: int, t: float) -> None:
        if not _due(self.every, sim, step):
            return
        chem = sim.plasma_chem  # type: ignore[attr-defined]
        if self._nu_x is None:
            if self.excimer not in chem.species_names:
                return
            self._nu_x = np.asarray(chem.nu[chem.species_names.index(self.excimer), :])
        ids, rates = sim.reaction_rates_mean()  # type: ignore[attr-defined]
        self._ids = ids
        self._t.append(t * 1e9)
        self._rows.append(rates / AVOGADRO)
        gain = getattr(sim, "gain", None)
        dt_last = float(getattr(gain, "last_dt_s", 0.0)) if gain is not None else 0.0
        if dt_last > 0.0:
            self._se.append(gain.last_se_sink_m3.get(self.excimer, 0.0) / dt_last / AVOGADRO)
            self._ste.append(gain.last_ste_sink_m3.get(self.excimer, 0.0) / dt_last / AVOGADRO)
        else:
            self._se.append(0.0)
            self._ste.append(0.0)

    def _split(self) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        rows = np.asarray(self._rows)  # (T, R) mol m^-3 s^-1
        form: dict[str, np.ndarray] = {}
        quench: dict[str, np.ndarray] = {}
        for r, rid in enumerate(self._ids):
            nu = float(self._nu_x[r])
            if nu > 0:
                form[rid] = rows[:, r] * nu
            elif nu < 0:
                quench[rid] = rows[:, r] * (-nu)
        if self._se:
            quench["SE"] = np.asarray(self._se)
            quench["STE"] = np.asarray(self._ste)
        return form, quench

    def export(self, sim: object, archive: object) -> None:
        if not self._t or self._nu_x is None:
            return
        t = np.asarray(self._t)
        form, quench = self._split()
        for name, series in (("formation", form), ("quenching", quench)):
            if series:
                archive.series(  # type: ignore[attr-defined]
                    f"excimer_rates_{name}", t, series, time_key="times_ns",
                )
                grp = archive.group(f"excimer_rates_{name}")  # type: ignore[attr-defined]
                grp.attrs["excimer"] = self.excimer
                grp.attrs["units"] = "mol m^-3 s^-1"
                grp.attrs["times_unit"] = "ns"
