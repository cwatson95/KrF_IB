"""Declared runtime parameters: validation, precedence, provenance.

ADAPTED, NOT RE-DERIVED, from ``PIC_MCC_3D/MODULES/params.py`` (itself
adapted from ``projects/MHD`` flare/params.py, FLASH guide ch. 10). The
mechanism is deliberately identical — declared parameters with doc/type/range,
fixed precedence ``default < problem < file < cli < api``, every value
remembering its source — because the three failures it closes are the same
ones: a misspelled override silently ignored, an out-of-range value failing
late and obscurely, and "why is this value what it is?" being unanswerable.

What differs here is only the declaration set at the bottom: these are the
knobs of a chemistry/radiation device engine, not of a PIC loop.

THE NO-CALIBRATION RULE IS ENFORCED AT DECLARATION TIME: a parameter whose
name or doc matches a fit/scale/calibration pattern is refused outright
(``_BANNED_KNOB``). phantom11's ``DIS_SCALE`` — auto-tuned to hit a target
density (defect class D4) — could not be declared in this store.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

SOURCES = ("default", "problem", "file", "cli", "api")

# Patterns that describe a fittable knob. Deliberately aggressive: a false
# positive costs a renamed parameter; a false negative costs defect class D4.
_BANNED_KNOB = re.compile(
    r"(calibrat|fudge|tune[d_]|_scale$|^scale_|scale_factor|fit_to|"
    r"target_match|anchor_adjust)", re.IGNORECASE)


class ParamError(ValueError):
    """Unknown name, bad value, banned knob, or constant reassignment."""


@dataclass
class Param:
    """One declared parameter."""

    default: object
    doc: str = ""
    type: type = None
    choices: tuple = None
    low: float = None
    high: float = None
    constant: bool = False
    units: str = ""

    def __post_init__(self):
        if self.type is None and self.default is not None:
            self.type = type(self.default)

    def coerce(self, name: str, value):
        if self.type is bool and isinstance(value, str):
            v = value.strip().lower()
            if v in ("true", ".true.", "yes", "1"):
                return True
            if v in ("false", ".false.", "no", "0"):
                return False
            raise ParamError(f"{name}: {value!r} is not a boolean")
        if self.type is not None and not isinstance(value, self.type):
            try:
                value = self.type(value)
            except (TypeError, ValueError):
                raise ParamError(f"{name}: cannot convert {value!r} to "
                                 f"{self.type.__name__}") from None
        if self.choices is not None and value not in self.choices:
            raise ParamError(f"{name}: {value!r} is not one of "
                             f"{list(self.choices)}")
        if self.low is not None and value < self.low:
            raise ParamError(f"{name}: {value} is below the minimum {self.low}"
                             + (f" ({self.doc})" if self.doc else ""))
        if self.high is not None and value > self.high:
            raise ParamError(f"{name}: {value} is above the maximum {self.high}"
                             + (f" ({self.doc})" if self.doc else ""))
        return value


@dataclass
class _Entry:
    param: Param
    value: object
    source: str
    owner: str = ""


class ParameterSet:
    """A validated, provenance-tracking parameter store."""

    def __init__(self) -> None:
        self._e = {}

    # ----------------------------------------------------------- declaration
    def declare(self, name: str, param: Param, owner: str = "") -> None:
        if name in self._e:
            raise ParamError(f"{name}: declared twice "
                             f"(by {self._e[name].owner!r} and {owner!r})")
        if _BANNED_KNOB.search(name) or _BANNED_KNOB.search(param.doc or ""):
            raise ParamError(
                f"{name}: name/doc matches a calibration-knob pattern — "
                "fittable scale factors are banned in live code (defect "
                "class D4, ORESTES_CLASS_REBUILD.txt §2)")
        self._e[name] = _Entry(param, param.coerce(name, param.default)
                               if param.default is not None else None,
                               "default", owner)

    def declare_many(self, mapping: dict, owner: str = "") -> "ParameterSet":
        for n, p in mapping.items():
            self.declare(n, p, owner)
        return self

    # -------------------------------------------------------------- access
    def __contains__(self, name: str) -> bool:
        return name in self._e

    def __getitem__(self, name: str):
        return self.get(name)

    def get(self, name: str):
        if name not in self._e:
            raise ParamError(f"unknown parameter {name!r}; "
                             f"declared: {sorted(self._e)}")
        return self._e[name].value

    def source(self, name: str) -> str:
        if name not in self._e:
            raise ParamError(f"unknown parameter {name!r}")
        return self._e[name].source

    def doc(self, name: str) -> str:
        return self._e[name].param.doc

    def changed(self, name: str) -> bool:
        return self._e[name].source != "default"

    def names(self) -> list:
        return sorted(self._e)

    def as_dict(self) -> dict:
        return {n: e.value for n, e in self._e.items()}

    def provenance(self) -> dict:
        return {n: e.source for n, e in self._e.items()}

    # --------------------------------------------------------------- setting
    def set(self, name: str, value, source: str = "api"):
        """Assign, honouring precedence. Unknown names are an ERROR."""
        if source not in SOURCES:
            raise ParamError(f"unknown source {source!r}; use {list(SOURCES)}")
        if name not in self._e:
            raise ParamError(
                f"unknown parameter {name!r} — a typo here would otherwise be "
                f"silently ignored. Declared: {sorted(self._e)}")
        e = self._e[name]
        if e.param.constant and e.source != "default":
            raise ParamError(f"{name} is constant and already set")
        if SOURCES.index(source) < SOURCES.index(e.source):
            return e.value
        e.value = e.param.coerce(name, value)
        e.source = source
        return e.value

    def set_many(self, mapping: dict, source: str = "api") -> "ParameterSet":
        for n, v in (mapping or {}).items():
            self.set(n, v, source)
        return self

    def copy(self) -> "ParameterSet":
        out = ParameterSet()
        for n, e in self._e.items():
            out._e[n] = _Entry(e.param, e.value, e.source, e.owner)
        return out

    # -------------------------------------------------------------- reporting
    def format_table(self, only_changed: bool = False) -> str:
        lines = []
        for n in self.names():
            e = self._e[n]
            if only_changed and e.source == "default":
                continue
            mark = "" if e.source == "default" else f"  [CHANGED: {e.source}]"
            unit = f" {e.param.units}" if e.param.units else ""
            lines.append(f"{n:<24} = {e.value!r}{unit}{mark}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# The engine's declared parameters. Ranges encode this package's actual
# validity limits so a run cannot silently sit outside them.
# ---------------------------------------------------------------------------
def engine_parameters() -> ParameterSet:
    """A fresh ParameterSet with the device engine's declarations."""
    P = ParameterSet()
    P.declare_many({
        # mode
        "mode": Param("krf", "network/species mode", str,
                      choices=("krf", "arf_chen", "arf_beam")),
        "drive": Param("beam", "energy source", str,
                       choices=("beam", "discharge")),
        "electron_tier": Param("loki", "electron-impact rate tier "
                               "(maxwell = the Maxwellian-EEDF VERIFICATION "
                               "tier, beam drive only, never a default)", str,
                               choices=("loki", "kinetic", "maxwell")),
        # CHEN-tier errata lane (0.1.25). True applies the adjudicated
        # corrections to the five NUMERIC transcription errors of Chen's
        # printed Table III (chen2025.CHEN_ERRATA: R13/R14 cm3-misprint,
        # R23 R27-paste, R24 x1e6 slip, R36 exponent typo — each
        # physically impossible or source-contradicted AS PRINTED); the
        # errata lane reproduces Chen's MODEL. False runs the table as
        # printed — the record of the PAPER, which (measured, 0.1.25)
        # cannot integrate under the paper's own base-case circuit.
        "chen_errata": Param(True, "apply chen2025.CHEN_ERRATA to the "
                             "CHEN-tier rows (False = as printed)", bool),
        # gas
        "pressure_atm": Param(3.0, "fill pressure", float, low=0.05, high=10.0,
                              units="[atm]"),
        "T_gas_K": Param(300.0, "initial gas temperature", float, low=77.0,
                         high=1000.0, units="[K]"),
        # time
        "t_end_ns": Param(200.0, "pulse window", float, low=1.0, high=5000.0,
                          units="[ns]"),
        "rtol": Param(1e-6, "BDF relative tolerance — the root Known-bad "
                      "forbids loose-tolerance tight comparisons", float,
                      low=1e-12, high=1e-3),
        # beam drive
        "beam_kV": Param(600.0, "beam voltage at the gas", float, low=50.0,
                         high=2000.0, units="[kV]"),
        "beam_j_A_cm2": Param(53.0, "beam current density at the gas", float,
                              low=0.0, high=1e4, units="[A/cm2]"),
        "beam_flat_ns": Param(60.0, "flat-top duration", float, low=0.0,
                              units="[ns]"),
        "beam_rise_ns": Param(5.0, "linear rise/fall", float, low=0.0,
                              units="[ns]"),
        "beam_axis": Param(0, "injection axis (0/1/2 = x/y/z; laser axis is "
                           "z)", int, low=0, high=2),
        "beam_t0_ns": Param(0.0, "beam pulse start (delay lets a microwave "
                            "window run BEFORE the beam — the pre-prime "
                            "timing study)", float, low=0.0, units="[ns]"),
        "beam2_t0_ns": Param(0.0, "second identical beam pulse start "
                             "(0 = single pulse; must clear the first "
                             "pulse — the inter-pulse/rep-rate timing "
                             "study)", float, low=0.0, units="[ns]"),
        # sub-table electron rates in beam mode (the L3 closure): 'relax'
        # consumes the MEASURED PIC_MCC_3D relaxing-population curves below
        # the swarm-table floor (kinetic_tier.RelaxationBridged — the
        # beam-mode EEDF class the tables cannot represent); 'taper' keeps
        # the flat+taper extrapolation (pre-0.1.9 behaviour, and the
        # measured x3.3 F- overestimate that goes with it).
        # measured sub-0.5 Td swarm extension (0.1.39). OFF by default:
        # turning it on moves every beam anchor and the recorded SPECKLE/
        # SPECKLE3D baselines, so it is the owner's call (the 0.1.36
        # petrov17_f2 precedent).
        "loki_subfloor": Param(False, "append the measured sub-0.5 Td LOKI "
                               "table extension", bool),
        "beam_subtable": Param("relax", "sub-table rate source for the beam "
                               "drive", str, choices=("relax", "taper")),
        # discharge drive (prescribed deposited power; circuit model is a
        # roadmap item, as in spectre v1)
        "dis_peak_W_m3": Param(5e11, "discharge power peak", float, low=0.0,
                               units="[W/m3]"),
        "dis_t0_ns": Param(100.0, "discharge pulse centre", float, units="[ns]"),
        "dis_fwhm_ns": Param(30.0, "discharge pulse FWHM", float, low=0.1,
                             units="[ns]"),
        "preion_m3": Param(1e14, "preionization seed density (an initial "
                           "condition, not a simulated flashover — the "
                           "PIC_MCC_3D M9 scope call)", float, low=0.0,
                           units="[m^-3]"),
        # circuit discharge drive (Chen 2025 Sec II A MPC topology;
        # E/N-direct rates, swarm-drift conduction — see MODULES/circuit.py).
        # Defaults are the paper's Sec II A / Table II BASE CASE, corrected
        # 0.1.25: the 0.1.5..0.1.24 values (10.4 nF / 8.1 nF / 30 nH and a
        # "5 -> 1 uH switch") were the Fig. 5 VALIDATION experiment's
        # components (2 atm, 20 mm, -13 kV, mirror-less) grafted onto the
        # base-case conditions — see CHANGELOG [0.1.25] and the Known-bad.
        "discharge_mode": Param("prescribed", "discharge energy pathway",
                                str, choices=("prescribed", "circuit")),
        "c0_nF": Param(10.0, "storage capacitor C0 (Chen Sec II A base "
                       "case: 10; E_in = C0 V0^2/2 = 1.62 J as the paper "
                       "states. 10.4 was the Fig. 5 experiment's)", float,
                       low=0.1, high=1e4, units="[nF]"),
        "cp_nF": Param(10.0, "peaking capacitor Cp (Chen Sec II A: 10; "
                       "8.1 was the Fig. 5 experiment's)", float,
                       low=0.1, high=1e4, units="[nF]"),
        "lm_uH": Param(0.25, "C0->Cp magnetic-switch inductance (Chen "
                       "Sec II A: Lm = 250 nH, 'a saturable inductor with "
                       "a low inductance'. The 5 -> 1 uH switching curve "
                       "belongs to the Fig. 5 validation circuit only)",
                       float, low=0.01, high=100.0, units="[uH]"),
        "ld_nH": Param(5.0, "Cp->gap stray inductance (Chen Sec II A: "
                       "Ls = 5; 30 was the Fig. 5 experiment's Ld)", float,
                       low=1.0, high=1e4, units="[nH]"),
        "v0_kV": Param(18.0, "C0 charge voltage magnitude (Chen base "
                       "case: 18)", float, low=0.1, high=100.0,
                       units="[kV]"),
        "r_loop_ohm": Param(0.2, "series loop resistance (copper/ESR "
                            "class, a physical circuit parameter)", float,
                            low=0.0, high=10.0, units="[ohm]"),
        # the paper's RLC filter branch across the gap (Sec. II A prints
        # Rf = 70.71 ohm, Lf = 100 nH, Cf = 10 pF, "to suppress circuit
        # oscillations after gap breakdown") — modeled as of 0.1.27
        "rf_ohm": Param(70.71, "gap filter-branch resistance (Chen Sec "
                        "II A: 70.71; 0 disables its damping)", float,
                        low=0.0, high=1e4, units="[ohm]"),
        "lf_nH": Param(100.0, "gap filter-branch inductance (Chen Sec "
                       "II A: 100)", float, low=0.1, high=1e5,
                       units="[nH]"),
        "cf_pF": Param(10.0, "gap filter-branch capacitance (Chen Sec "
                       "II A: 10)", float, low=0.01, high=1e6,
                       units="[pF]"),
        # microwave reheat window (an afterglow overlay on the BEAM drive):
        # a prescribed INJECTED power, trapezoidal in time. Each RHS step
        # the uniform cavity field settles where the swarm power balance
        # absorbs it — e E v_d(E/N) x nu^2/(nu^2+omega^2) summed over the
        # plasma electrons — capped at mw_en_max_td (the empty-cavity
        # field the wall Q supports). Dense plasma -> the field sags and
        # per-electron heating vanishes; dilute plasma -> the field rides
        # the ceiling and the absorbed fraction dies with n_e. Absorbed
        # (not injected) power is what enters drive_in. engine._mw block.
        "mw_power_W": Param(0.0, "microwave window injected power "
                            "(0 = off; beam drive only)", float,
                            low=0.0, high=1e9, units="[W]"),
        "mw_en_max_td": Param(30.0, "cavity-field ceiling as E/N (the "
                              "empty-cavity buildup at the feed power and "
                              "wall Q; 30 Td ~ 1 MW into the 16.7 L "
                              "chamber at Q ~ 1e4, 2 atm)", float,
                              low=0.1, high=300.0, units="[Td]"),
        "mw_t0_ns": Param(130.0, "microwave window start", float, low=0.0,
                          units="[ns]"),
        "mw_flat_ns": Param(150.0, "microwave window flat-top", float,
                            low=0.0, units="[ns]"),
        "mw_rise_ns": Param(10.0, "microwave window linear rise/fall",
                            float, low=0.0, units="[ns]"),
        "mw_freq_GHz": Param(9.5, "microwave frequency (WR90 X-band "
                             "mid-band); sets the collisional AC factor "
                             "nu^2/(nu^2+omega^2) and the nu/omega "
                             "validity gauge", float, low=0.5, high=100.0,
                             units="[GHz]"),
        # resonator design (physical device choices; 0 / negative keep the
        # mode's default optics — krf: R_OC 0.1, L_op 0.088 fixed;
        # arf: R_OC 0.08, L_op = k_c x C_F2)
        "r_oc": Param(0.0, "output-coupler reflectivity (0 = mode "
                      "default)", float, low=0.0, high=0.999),
        "l_op": Param(-1.0, "fixed single-pass optics loss (negative = "
                      "mode default)", float, low=-1.0, high=0.9),
        # geometry (device grid; z is the optical axis)
        "nx": Param(1, "device-grid cells along x", int, low=1, high=256),
        "ny": Param(1, "device-grid cells along y", int, low=1, high=256),
        "nz": Param(1, "device-grid cells along z", int, low=1, high=1024),
        "Lx_m": Param(0.014, "cell extent x", float, low=1e-4, units="[m]"),
        "Ly_m": Param(0.014, "cell extent y", float, low=1e-4, units="[m]"),
        "Lz_m": Param(0.75, "cell extent z (gain length)", float, low=1e-3,
                      units="[m]"),
        # in-loop ASE (S6): the 4-pi amplified-spontaneous-emission flux
        # of the current gain field depopulates the upper level in the RHS
        # (first-principles: sum-rule bandwidth, closed-form chords on the
        # volume-averaged field — the approximations are stated in
        # engine._ase_rate). Off reproduces the pre-0.1.6 behaviour.
        "ase_inloop": Param(True, "couple ASE depopulation into the RHS",
                            bool),
        # beam-sustained-EEDF rates for the low-threshold e+F2 channels
        # (vibrational sink, repulsive dissociation, dissociative
        # attachment): while the beam envelope is ON, the <eps>-indexed
        # swarm coefficients under-drive the 3-4 eV repulsive channels
        # the secondary tail pumps (the F2-burn x0.43 / F- residuals).
        # Petrov 2017 JAP 122 133301 Table I fits k(P, x_F2) over the
        # beam-sustained EEDF (their Fig. 6; NRL-computed, nothing tuned
        # here) replace those coefficients for the marked rows; beam-off
        # phases revert to the swarm rates (one owner per population).
        # Default OFF: flipping it moves the beam anchors — measured
        # before/after lives in CHANGELOG [0.1.36].
        "petrov17_f2": Param(False, "beam-EEDF e+F2 low-threshold rates "
                             "(Petrov 2017 Table I) while the beam is on",
                             bool),
        # inter-cell transport (conservative 6-neighbour diffusion).
        # MEASURED negligible at device presets (cm cells, <= 200 ns:
        # tau_diff ~ 0.1 s) — available and gated, default off.
        "transport": Param(False, "species diffusion between cells", bool),
        "diff_m2_s": Param(2e-5, "diffusion coefficient at 1 atm, 300 K "
                           "(neutral-neutral class; EST bracket x1/3..x3; "
                           "scaled by (T_g/300)^1.5 / p)", float, low=0.0,
                           high=1e-2, units="[m^2/s]"),
        # seeds
        "te_seed_eV": Param(0.3, "initial electron temperature", float,
                            low=0.026, high=20.0, units="[eV]"),
        "seed_density_m3": Param(1e12, "seed density for minor species",
                                 float, low=1.0, units="[m^-3]"),
    }, owner="Device")
    return P
