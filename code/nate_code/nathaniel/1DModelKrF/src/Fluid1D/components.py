"""Composable specification objects for a 1D charged-species fluid excimer sim (Keras-like).

Each spec is a small, frozen dataclass describing one facet of the experiment — 1D mesh,
gas fill, initial plasma, optical cavity, Boltzmann-coefficient source, discharge circuit,
and numerics — so an experiment file reads as a declarative composition rather than a flat
bag of parameters. ``Simulation`` consumes these; diagnostics are supplied separately as a
list (see ``Fluid1D.diagnostics``).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from numpy.typing import NDArray

from particles.constants import ATM_PA, K_B
_VOLUME_LENGTH_M = {"cm": 1.0e-2, "m": 1.0}  # metres per user length unit


@dataclass(frozen=True)
class Units:
    """User-facing volume/length unit system for the experiment (``"cm"`` or ``"m"``).

    The engine always computes in SI (m, m⁻³) and rate constants / cross-sections carry their own
    units in the data files (converted to SI at load) — those are NEVER touched by this. ``Units``
    only rescales the **volumetric/spatial** quantities the user enters (mesh geometry, seed
    densities) on the way in, and the ones reported back (densities, lengths, gain/absorption
    coefficients) on the way out, so a user can work entirely in cm or entirely in m. Electrical
    (V, A, F, H), thermal (K, eV, J), and reduced-field (Td) quantities stay SI.

    Conversions use a single length scale ``length_m`` (metres per user unit): lengths ×L,
    areas ×L², volumes ×L³, number/charge densities ÷L³, inverse lengths ÷L (user→SI).
    """

    volume: str = "m"

    def __post_init__(self) -> None:
        if self.volume not in _VOLUME_LENGTH_M:
            raise ValueError(f"volume units must be 'cm' or 'm', got {self.volume!r}")

    @property
    def length_m(self) -> float:
        return _VOLUME_LENGTH_M[self.volume]

    # -- user -> SI (inputs) --
    def length_to_si(self, x: float) -> float:
        return x * self.length_m

    def area_to_si(self, a: float) -> float:
        return a * self.length_m**2

    def density_to_si(self, n: float) -> float:
        return n / self.length_m**3  # unit⁻³ -> m⁻³

    def inv_length_to_si(self, k: float) -> float:
        return k / self.length_m  # unit⁻¹ -> m⁻¹

    # -- SI -> user (outputs) --
    def length_from_si(self, x):
        return x / self.length_m

    def area_from_si(self, a):
        return a / self.length_m**2

    def density_from_si(self, n):
        return n * self.length_m**3  # m⁻³ -> unit⁻³

    def inv_length_from_si(self, k):
        return k * self.length_m  # m⁻¹ -> unit⁻¹

    # -- labels --
    @property
    def length_label(self) -> str:
        return self.volume

    @property
    def density_label(self) -> str:
        return f"{self.volume}$^{{-3}}$"  # matplotlib

    @property
    def density_unit(self) -> str:
        return f"{self.volume}^-3"  # plain text (archive attrs)

    @property
    def inv_length_label(self) -> str:
        return f"1/{self.volume}"


@dataclass(frozen=True)
class Grid1D:
    """1D mesh in the experiment's declared length units.

    ``length`` and ``area_m2`` are converted by :class:`Units` when the simulation builds;
    despite its name, ``area_m2`` follows the square of the declared length unit.
    """

    cells: int = 128
    length: float = 2.0e-2
    area_m2: float = 1.0

    def __post_init__(self) -> None:
        if not isinstance(self.cells, int) or self.cells < 2:
            raise ValueError("Grid1D.cells must be an integer >= 2")
        if not np.isfinite(self.length) or self.length <= 0.0:
            raise ValueError("Grid1D.length must be finite and positive")
        if not np.isfinite(self.area_m2) or self.area_m2 <= 0.0:
            raise ValueError("Grid1D.area_m2 must be finite and positive")


@dataclass(frozen=True)
class GasMixture:
    """Neutral background: mole fractions per species + pressure/temperature.

    Fractions are by number (partial pressure); e.g. ``{"Kr": 0.99, "F2": 0.01}`` or a
    buffered fill ``{"He": 0.90, "Kr": 0.095, "F2": 0.005}``. Every named species must be a
    neutral in the registry (validated at build). With ``normalize=True`` (default) the
    fractions are renormalized to sum to 1; with ``normalize=False`` they must already sum
    to ~1. ``with_fraction`` is the F2-sweep knob.
    """

    fractions: dict[str, float]
    pressure_atm: float = 1.0
    temperature_k: float = 300.0
    normalize: bool = True

    def __post_init__(self) -> None:
        if not self.fractions:
            raise ValueError("GasMixture requires at least one species fraction")
        for name, frac in self.fractions.items():
            if not np.isfinite(frac) or frac <= 0.0:
                raise ValueError(f"GasMixture fraction for {name!r} must be > 0 (got {frac})")
        total = sum(self.fractions.values())
        if not np.isfinite(self.pressure_atm) or self.pressure_atm <= 0.0:
            raise ValueError("GasMixture.pressure_atm must be finite and positive")
        if not np.isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise ValueError("GasMixture.temperature_k must be finite and positive")
        if not self.normalize and abs(total - 1.0) > 1e-6:
            raise ValueError(
                f"GasMixture fractions must sum to 1 when normalize=False (sum={total:.6f})"
            )

    @property
    def normalized_fractions(self) -> dict[str, float]:
        total = sum(self.fractions.values())
        return {name: frac / total for name, frac in self.fractions.items()}

    @property
    def number_density(self) -> float:
        """Total neutral number density n = P / (k_B T) (m^-3)."""
        return self.pressure_atm * ATM_PA / (K_B * self.temperature_k)

    def partial_density(self, name: str) -> float:
        return self.normalized_fractions[name] * self.number_density

    def with_fraction(self, name: str, value: float, *, rebalance: str = "Kr") -> "GasMixture":
        """Copy with ``name`` set to ``value``, absorbing the delta into ``rebalance``.

        Used by the F2 sweep: ``mix.with_fraction("F2", 0.01)`` shifts the change into the
        Kr partial pressure so the total (and any buffer gas) is preserved.
        """
        new = dict(self.fractions)
        delta = value - new.get(name, 0.0)
        new[name] = value
        new[rebalance] = new.get(rebalance, 0.0) - delta
        return replace(self, fractions=new)


@dataclass(frozen=True)
class Plasma:
    """Initial charged fill for the 1D fluid model (quasi-neutral electrons + one ion).

    The charged species are per-cell density fields. The *spatial* profile of the seed
    electron/ion density is set by :class:`Preionization` (defaults to a
    uniform ``density``); ``ion`` names the balancing positive species and
    ``electron_mean_energy_ev`` the seed ⟨ε⟩.
    """

    density: float = 1.0e18  # seed density in the declared volume unit^-3
    ion: str = "Kr+"
    electron_mean_energy_ev: float = 4.0  # seed <eps>; Te = (2/3)<eps>
    # Cathode secondary-emission coefficient γ: electrons injected per incident positive
    # ion at the local mean energy. Chen (2025) does not state this input; it is exposed as
    # a breakdown-timing sensitivity parameter.
    secondary_emission_gamma: float = 0.05

    def __post_init__(self) -> None:
        if not np.isfinite(self.density) or self.density < 0.0:
            raise ValueError("Plasma.density must be finite and non-negative")
        if not np.isfinite(self.electron_mean_energy_ev) or self.electron_mean_energy_ev < 0.0:
            raise ValueError("Plasma.electron_mean_energy_ev must be finite and non-negative")
        if not np.isfinite(self.secondary_emission_gamma) or not (
            0.0 <= self.secondary_emission_gamma < 1.0
        ):
            raise ValueError("Plasma.secondary_emission_gamma must be finite in [0, 1)")


@dataclass(frozen=True)
class Preionization:
    """Initial electron pre-ionization profile in the declared experiment units.

    A necessary condition for uniform breakdown at high pressure — corona/spark pre-ionizers
    near the cathode seed a spatially non-uniform electron density. Modes:

    - ``"uniform"``: ``n_{e,0}(x) = density`` everywhere (the flat seed).
    - ``"tanh"``: ``n_{e,0}(x) = amplitude·tanh(steepness·x) + floor`` — the standard profile
      used by Chen et al. (Eq. 6), rising from ``floor`` at ``x=0`` toward the far electrode.
      Defaults are interpreted through :class:`Units` when the simulation builds.

    ``density``, ``amplitude``, and ``floor`` use the declared volume unit⁻³;
    ``steepness`` uses the inverse declared length unit. The built SI profile rises from the
    grounded anode at x=0 toward the powered cathode at x=length.
    """

    mode: str = "uniform"
    density: float = 1.0e18
    amplitude: float = 1.0e14
    floor: float = 1.0e13
    steepness: float = 100.0

    def __post_init__(self) -> None:
        if self.mode not in {"uniform", "tanh"}:
            raise ValueError(f"unknown preionization mode {self.mode!r}")
        values = (self.density, self.amplitude, self.floor, self.steepness)
        if not all(np.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("preionization densities and steepness must be finite and non-negative")

    def profile(self, cell_centers: NDArray[np.floating]) -> NDArray[np.float64]:
        """Evaluate n_{e,0}(x) (m⁻³) at the given axial cell-centre coordinates (m)."""
        x = np.asarray(cell_centers, dtype=np.float64)
        if self.mode == "uniform":
            return np.full_like(x, self.density)
        if self.mode == "tanh":
            return self.amplitude * np.tanh(self.steepness * x) + self.floor
        raise ValueError(f"unknown preionization mode {self.mode!r} (use 'uniform' or 'tanh')")


@dataclass(frozen=True)
class Cavity:
    """Optical cavity for the 0D excimer gain/output model (Chen et al. 2025 §II.C).

    The laser beam runs along the **gain axis** (``gain_length_m``, l_g) — transverse to the 1D
    discharge axis — bounced between a rear mirror and the output coupler (reflectivity ``r_out``
    = R_OC). The cavity photon density is evolved by Chen Eq. 7 and the output power by Eqs. 10/18
    (see :mod:`optics.gain` / :mod:`optics.photon_density`).

    ``absorption`` maps extra UV-absorber species to photoabsorption cross-sections (m^2); it is
    **unioned** with the network's PHOTO reactions (:func:`reactions.photo_absorbers`). ArF sources
    its absorbers (Ar**, F⁻) from ``arf.toml``, so its ``absorption`` is usually empty; a KrF fill
    typically passes ``{"F2": 4.5e-24, "F-": 1.7e-21, "Kr*": 3.9e-23}``.

    ``cavity_length_m`` (l_c) is the full optical cavity, ``gain_length_m`` (l_g) the pumped gain
    region, ``beam_area_m2`` (S) the output-beam cross-section (``None`` → choose S so the pumped
    gain volume ``l_g*S`` equals the simulated discharge volume).
    ``axial_se_fraction`` (K) is the fraction of isotropic spontaneous emission launched along the
    axis; ``loss_coeff_per_f2`` (k_c) sets the single-pass optical loss ``L_op = k_c·x_F2`` (Chen
    Eq. 9). Defaults are Chen's Table II ArF values.
    """

    wavelength_nm: float = 248.0
    r_out: float = 0.1  # output-coupler reflectivity R_OC (1 - r_out transmits as laser output)
    r_back: float = 0.99  # rear-mirror reflectivity
    absorption: dict[str, float] = field(default_factory=dict)  # extra absorbers: species -> sigma (m^2)
    gain_length_m: float = 0.90  # l_g — pumped gain length along the optical axis
    cavity_length_m: float = 1.05  # l_c — full optical cavity length
    beam_area_m2: float | None = None  # S; None -> l_g*S equals the simulated discharge volume
    axial_se_fraction: float = 1.0e-5  # K — axial fraction of isotropic spontaneous emission
    loss_coeff_per_f2: float = 80.0  # k_c — single-pass loss L_op = k_c · x_F2 (Chen Eq. 9)

    def __post_init__(self) -> None:
        if not np.isfinite(self.wavelength_nm) or self.wavelength_nm <= 0.0:
            raise ValueError("Cavity.wavelength_nm must be finite and positive")
        if not (0.0 < self.r_out < 1.0 and 0.0 < self.r_back <= 1.0):
            raise ValueError("cavity mirror reflectivities must lie in (0, 1]")
        if (not np.isfinite(self.gain_length_m) or self.gain_length_m <= 0.0
                or not np.isfinite(self.cavity_length_m)
                or self.cavity_length_m < self.gain_length_m):
            raise ValueError("cavity lengths require 0 < gain_length_m <= cavity_length_m")
        if self.beam_area_m2 is not None and (
                not np.isfinite(self.beam_area_m2) or self.beam_area_m2 <= 0.0):
            raise ValueError("Cavity.beam_area_m2 must be finite and positive")
        if not (0.0 <= self.axial_se_fraction <= 1.0):
            raise ValueError("Cavity.axial_se_fraction must lie in [0, 1]")
        if not np.isfinite(self.loss_coeff_per_f2) or self.loss_coeff_per_f2 < 0.0:
            raise ValueError("Cavity.loss_coeff_per_f2 must be finite and non-negative")
        if any(not np.isfinite(value) or value < 0.0 for value in self.absorption.values()):
            raise ValueError("Cavity absorption cross sections must be finite and non-negative")


@dataclass(frozen=True)
class BoltzmannSpec:
    """How to build/look up the electron Boltzmann LUT (transport + rate coefficients).

    ``backend``: ``"two_term"`` (:class:`Fluid1D.two_term.TwoTermBoltzmannSolver` — an in-tree
    two-term Boltzmann solve with non-Maxwellian rates),
    ``"eedf"`` (:class:`Fluid1D.eedf_solver.MaxwellianEedfSolver` — Maxwellian-EEDF moments of
    the same cross-sections; faster build, over-populated tail), or ``"stub"``
    (AnalyticStubSolver, placeholder coefficients for fast tests). ``en_grid_td`` is
    ``(min, max, points)`` for the reduced-field (E/N, townsend) sweep, log-spaced.
    ``key_kind`` selects the runtime lookup axis (``"EN_Td"`` or ``"mean_energy_eV"``);
    ``cache_dir`` caches one LUT per mixture.
    """

    backend: str = "stub"
    en_grid_td: tuple[float, float, int] = (1.0, 1000.0, 64)
    key_kind: str = "EN_Td"
    cache_dir: str | None = None

    def __post_init__(self) -> None:
        if self.backend not in {"stub", "eedf", "two_term"}:
            raise ValueError("BoltzmannSpec.backend must be 'stub', 'eedf', or 'two_term'")
        lo, hi, count = self.en_grid_td
        if not (np.isfinite(lo) and np.isfinite(hi) and 0.0 < lo < hi and int(count) >= 2):
            raise ValueError("BoltzmannSpec.en_grid_td requires 0 < min < max and count >= 2")
        if self.key_kind not in {"EN_Td", "mean_energy_eV"}:
            raise ValueError("BoltzmannSpec.key_kind must be 'EN_Td' or 'mean_energy_eV'")


@dataclass(frozen=True)
class CircuitSpec:
    """External discharge-circuit configuration feeding the electrode Dirichlet BCs.

    ``backend``: ``"stub"`` (fixed/prescribed electrode potentials via ``phi_left``/
    ``phi_right``) or ``"mpc"`` (the magnetic-pulse-compression driver of Chen et al.,
    Phys. Plasmas 32, 113502 (2025) — see :class:`Fluid1D.field.MpcCircuit`). The ``mpc``
    component values default to the paper's; all are modifiable. The gap capacitance is
    computed from the mesh geometry (``ε₀·A/d``) plus the lumped parasitic ``c_head_f``.
    """

    backend: str = "stub"
    # --- stub backend ---
    phi_left: float = 0.0
    phi_right: float = -20.0e3
    # --- mpc backend (Chen et al. 2025 defaults) ---
    c0_f: float = 10.0e-9  # energy-storage capacitor (F)
    lm_h: float = 250.0e-9  # magnetic-switch inductor (H)
    cp_f: float = 10.0e-9  # peaking capacitor (F)
    ls_h: float = 5.0e-9  # discharge stray inductor (H)
    rf_ohm: float = 70.71  # snubber resistance (Ω) — optimally-damped R_f = sqrt(L_f/2C_f)
    lf_h: float = 100.0e-9  # snubber inductor (H)
    cf_f: float = 10.0e-12  # snubber capacitor (F)
    v0_v: float = -18.0e3  # C₀ pre-charge voltage (V, negative)
    r1_ohm: float = 0.02  # MPC-loop bus resistance (Ω)
    rs_ohm: float = 0.02  # discharge-loop bus resistance (Ω)
    c_head_f: float = 20.0e-12  # lumped electrode/feedthrough parasitic capacitance (F)

    def __post_init__(self) -> None:
        if self.backend not in {"stub", "mpc"}:
            raise ValueError("CircuitSpec.backend must be 'stub' or 'mpc'")
        if not all(np.isfinite(value) for value in (self.phi_left, self.phi_right, self.v0_v)):
            raise ValueError("circuit voltages must be finite")
        positive = (self.c0_f, self.lm_h, self.cp_f, self.ls_h, self.lf_h, self.cf_f)
        if not all(np.isfinite(value) and value > 0.0 for value in positive):
            raise ValueError("circuit capacitances and inductances must be positive")
        nonnegative = (self.rf_ohm, self.r1_ohm, self.rs_ohm, self.c_head_f)
        if not all(np.isfinite(value) and value >= 0.0 for value in nonnegative):
            raise ValueError("circuit resistances and parasitic capacitance must be non-negative")


@dataclass(frozen=True)
class Numerics:
    """Time-stepping controls passed to :meth:`Simulation.run`.

    ``field_heating_cap_factor`` limits resolved electron ``J.E`` work to that multiple of
    terminal gap power (plus the small bootstrap allowance used by the transport closure).
    ``None`` disables the audited supply cap for sensitivity studies; it is not a claim that
    an uncapped under-resolved sheath is energy-conservative.

    ``electron_energy_closure="source_coupled"`` is the production closure.  The alternative
    ``"continuous_reconciled"`` is a Chen-equivalence diagnostic: transport applies continuous
    inelastic cooling, then the chemistry interval records a signed reconciliation against its
    accepted reaction cost.
    """

    steps: int = 2000
    dt: float = 1.0e-12  # s
    fluid_every: int = 10  # neutral-chemistry + photon update cadence (steps)
    fluid_substeps: int = 20  # implicit substeps per fluid-chemistry call
    field_heating_cap_factor: float | None = 2.0
    electron_energy_closure: str = "source_coupled"

    def __post_init__(self) -> None:
        if not isinstance(self.steps, int) or self.steps <= 0:
            raise ValueError("Numerics.steps must be a positive integer")
        if not np.isfinite(self.dt) or self.dt <= 0.0:
            raise ValueError("Numerics.dt must be finite and positive")
        if not isinstance(self.fluid_every, int) or self.fluid_every <= 0:
            raise ValueError("Numerics.fluid_every must be a positive integer")
        if not isinstance(self.fluid_substeps, int) or self.fluid_substeps <= 0:
            raise ValueError("Numerics.fluid_substeps must be a positive integer")
        if self.field_heating_cap_factor is not None and (
            not np.isfinite(self.field_heating_cap_factor)
            or self.field_heating_cap_factor <= 0.0
        ):
            raise ValueError("field_heating_cap_factor must be finite and positive or None")
        if self.electron_energy_closure not in {"source_coupled", "continuous_reconciled"}:
            raise ValueError(
                "electron_energy_closure must be 'source_coupled' or 'continuous_reconciled'"
            )


# Default progress-line fields for the fluid engine. ``excimer`` resolves to the built
# network's lasing upper state (KrF*/ArF*), so the default line is system-aware.
DEFAULT_PROGRESS_FIELDS: tuple[str, ...] = (
    "excimer", "photons", "temp:e-", "efficiency",
)


@dataclass(frozen=True)
class ProgressReport:
    """Progress reporting configuration for :meth:`Simulation.run`.

    ``every`` is the report cadence in timesteps (``None`` = one line per 10% of the run);
    the final step always reports. ``style`` selects the leading column: ``"percent"`` or
    ``"steps"``. ``fields`` picks reported quantities from the vocabulary in
    :mod:`Fluid1D.progress` — per-species ``total:/density:/temp:`` plus
    ``excimer, photons, efficiency, gain, p_in, p_out``. ``timing`` appends
    ``elapsed=``/``eta=`` wall-clock columns to the terminal line. Sinks: ``terminal``
    prints the per-cadence lines; ``heartbeat`` atomically rewrites
    ``<run_dir>/progress.json`` on the same cadence plus a final status record. At least
    one sink must stay enabled. Passing a spec to :meth:`Simulation.run` enables reporting
    even under ``verbose=False`` (banner and summary stay quiet).
    """

    every: int | None = None
    style: str = "percent"
    fields: tuple[str, ...] = DEFAULT_PROGRESS_FIELDS
    timing: bool = True
    terminal: bool = True
    heartbeat: bool = False
