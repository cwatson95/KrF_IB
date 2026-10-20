"""Declared runtime parameters: validation, precedence, and PROVENANCE.

ADAPTED FROM ``projects/MHD`` (``flare/params.py``), which recreates FLASH's
runtime-parameter system (guide ch. 10): every unit declares its parameters
with a doc string and a validity range, overrides are applied in a fixed
precedence order, and the resolved set remembers *where each value came from*.

WHY THIS PACKAGE NEEDED IT.  Until now every knob was a bare keyword argument
threaded through constructors, with validation scattered across them and no
record of who set what.  Three concrete failures that buys you:

  * **A misspelled parameter was silently ignored.**  ``ESSim(..., ppc=100)``
    would raise, but a typo in a dict of overrides would not -- the value just
    never took effect and the run quietly used the default.  :meth:`set`
    rejects unknown names outright.
  * **Out-of-range values failed late and obscurely.**  ``ppw = 4`` surfaced
    as a strange dispersion result, not as an error.  Ranges are declared once
    and checked at assignment.
  * **"Why is dt this value?" was unanswerable.**  Now
    ``ps.source("dt")`` says ``default`` / ``problem`` / ``file`` / ``cli`` /
    ``api``, and :meth:`format_table` prints the whole resolved set with
    ``[CHANGED]`` markers -- the same convention FLASH's logfile uses.

Precedence, lowest to highest: ``default < problem < file < cli < api``.  A
lower-precedence source cannot overwrite a value already set by a higher one,
so a problem definition cannot silently clobber what the caller passed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# lowest .. highest.  A set() from a weaker source than the current one is
# ignored (not an error): that is what makes layering defaults safe.
SOURCES = ("default", "problem", "file", "cli", "api")


class ParamError(ValueError):
    """Raised for an unknown name, a bad value, or a constant reassignment."""


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
        """Type-convert and validate, or raise :class:`ParamError`."""
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
                raise ParamError(
                    f"{name}: cannot convert {value!r} to "
                    f"{self.type.__name__}") from None
        if self.choices is not None and value not in self.choices:
            raise ParamError(
                f"{name}: {value!r} is not one of {list(self.choices)}")
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
        self._e: dict = {}

    # ----------------------------------------------------------- declaration
    def declare(self, name: str, param: Param, owner: str = "") -> None:
        if name in self._e:
            raise ParamError(f"{name}: declared twice "
                             f"(by {self._e[name].owner!r} and {owner!r})")
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
        """Where the current value came from -- the debugging question."""
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
        """Assign, honouring precedence. Unknown names are an ERROR, not a no-op."""
        if source not in SOURCES:
            raise ParamError(f"unknown source {source!r}; use {list(SOURCES)}")
        if name not in self._e:
            raise ParamError(
                f"unknown parameter {name!r} -- a typo here would otherwise "
                f"be silently ignored. Declared: {sorted(self._e)}")
        e = self._e[name]
        if e.param.constant and e.source != "default":
            raise ParamError(f"{name} is constant and already set")
        if SOURCES.index(source) < SOURCES.index(e.source):
            return e.value                    # weaker source: keep what we have
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
        """The resolved set, FLASH-style, with [CHANGED] markers."""
        lines = []
        for n in self.names():
            e = self._e[n]
            if only_changed and e.source == "default":
                continue
            mark = "" if e.source == "default" else f"  [CHANGED: {e.source}]"
            unit = f" {e.param.units}" if e.param.units else ""
            lines.append(f"{n:<24} = {e.value!r}{unit}{mark}")
        return "\n".join(lines)


# --------------------------------------------------------------------------
# The declared parameters of an electrostatic PIC-MCC run.  Ranges encode the
# physics limits this package actually has, so a run cannot silently sit
# outside them -- e.g. ppw below 8 makes a dispersion comparison meaningless,
# and wpe_dt above 2 is unconditionally unstable for the explicit mover.
# --------------------------------------------------------------------------
def sim_parameters() -> ParameterSet:
    """A fresh ParameterSet with the simulation's declarations."""
    P = ParameterSet()
    P.declare_many({
        # grid
        "nx": Param(32, "cells along x", int, low=3, high=1024),
        "ny": Param(32, "cells along y", int, low=3, high=1024),
        "nz": Param(32, "cells along z", int, low=3, high=1024),
        "h": Param(1.0e-3, "cell size", float, low=0.0, units="[m]"),
        "bc": Param("periodic_xy_dirichlet_z", "boundary geometry", str,
                    choices=("periodic", "periodic_xy_dirichlet_z")),
        # time
        "dt": Param(1.0e-12, "field timestep", float, low=0.0, units="[s]"),
        "n_steps": Param(100, "steps to run", int, low=0),
        "wpe_dt_max": Param(0.4, "explicit-mover stability ceiling on "
                            "omega_pe*dt (unstable above 2)", float,
                            low=0.0, high=2.0),
        # plasma
        "n_e": Param(1.0e14, "seed electron density", float, low=0.0,
                     units="[m^-3]"),
        "Te": Param(1.0, "seed electron temperature", float, low=0.0,
                    units="[eV]"),
        "ppc": Param(50, "macro-particles per cell", int, low=1),
        "seed": Param(0, "RNG seed -- part of the run's identity", int),
        # electrodes
        "phi_lo": Param(0.0, "potential at z = 0", float, units="[V]"),
        "phi_hi": Param(0.0, "potential at z = Lz", float, units="[V]"),
        # gas / collisions
        "gas": Param("Ar", "dominant neutral species", str),
        "pressure_atm": Param(1.0, "fill pressure", float, low=0.0,
                              units="[atm]"),
        "T_gas": Param(300.0, "neutral temperature", float, low=1.0,
                       units="[K]"),
        "mcc_on": Param(True, "run the electron-neutral MCC", bool),
        "ion_particles": Param(True, "ionization makes kinetic ions "
                               "(False = immobile grid charge)", bool),
    }, owner="ESSim")
    return P
