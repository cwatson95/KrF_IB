"""MAXWELL tier — Maxwellian and two-Maxwellian EEDF rate coefficients
integrated from the LXCat cross sections, as a VERIFICATION of the LOKI
swarm-shape tables (S4 class: reported side by side, never averaged).

WHY (owner ask 2026-09-03). The LOKI tier reads its rates from LoKI-MC
swarm tables — field-driven, Druyvesteyn-class EEDF shapes computed
WITHOUT an electron-electron operator (KrF_and_LoKI/docs/two-region-
eedf.md §1) — and looks them up by <eps> = 1.5 Te. Where electron-electron
collisions dominate the energy exchange (a beam-pumped bulk at
n_e ~1e14-1e15 cm^-3 and Te ~1 eV: nu_ee ~1e9-1e10 s^-1 against
~1e7-1e8 s^-1 of elastic + vibrational cooling) the sub-threshold bulk is
Maxwellian, not swarm-shaped. This tier integrates every electron-impact
bucket over a Maxwellian at the engine's Te from the SAME cross-section
files the tables were computed from, so the LOKI/MAXWELL difference
isolates the EEDF-shape assumption and nothing else.

WHAT IT IS NOT. Not a better closure by construction, and never a
default. It is Te-indexed only: a Maxwellian carries no field -> EEDF
relation, so the circuit / E/N-direct path (the discharge lane) cannot
run on it — the engine refuses that combination at construction, and
``k_by_en`` raises. For the discharge the comparison is rate-level:
tools/maxwell_verify.py indexes the same ratios by E/N through the swarm
table's own <eps>(E/N). The two-Maxwellian form (bulk + tail at FIXED
mean energy) is a sweep coordinate for the shape bracket — the tail
fraction and temperature are values stated by the caller of a tool,
never parameters of the engine (the no-knob rule, params.py).

Provenance. Cross sections = the sets cross_sections.py resolves
(IST-Lisbon Ar/He, Biagi-v7.1 Kr, SIGLO F2, Biagi/Magboltz-8.9 Ne) — the
files behind MODULES/loki_tables/ (its README.txt). Bucket names and
rules mirror loki_tier._bucket on the LXCat labels/thresholds:
  Ar: 3P2/3P1/3P0/1P1 (4s, 11.5-11.8 eV) -> Ar*, every higher level ->
      Ar**; Kr: threshold < 11 eV (the 5s manifold, 9.9 eV) -> Kr*, else
      Kr**; Ne: threshold < 17.5 eV (3s, 16.6-16.85 eV) -> Ne*, else
      Ne**; F2: V1..V4 -> vib, a3 -> F2(a3), A1/C1/H1 -> F2(A1),
      attachment -> F-, ionization -> F2+; 'Elastic' -> el:<gas>.
F2 ships only an Effective set; CrossSectionSet unfolds it into
elastic = effective - sum(inelastic), clipped (the pic1d remedy; LoKI-MC
performs the same class of unfolding, root CHANGELOG Known-bad on the
Phelps Ne set). Xe has no sigma file and is skipped (recorded on the
instance) exactly as the beam module skips it.

Numerics. k(T) = int sigma(eps) v(eps) f_M(eps; T) d eps with
f_M = (2/sqrt(pi)) T^-3/2 sqrt(eps) exp(-eps/T) on the 8192-point
log grid 1e-5..2e3 eV (trapezoid in eps); gated against the closed forms
for a constant and a step cross section (tests/test_maxwell_tier.py).
"""

from __future__ import annotations

import numpy as np

from .cross_sections import CrossSectionSet, gas_file
from .units import EV, M_E

#: the energy grid shared by every integral here (eV): 1e-5 so the
#: truncated low-energy tail is < 1e-6 of a 0.026 eV (Te-floor) integral
EPS_GRID = np.geomspace(1e-5, 2e3, 8192)
# numpy >= 2.0 renamed trapz; the two-interpreter rule runs the Homebrew
# smoke with -W error::DeprecationWarning, so resolve the name once
_trapz = getattr(np, "trapezoid", None) or np.trapz

_ION = {"Ar": "Ar+", "Kr": "Kr+", "Ne": "Ne+", "F2": "F2+", "He": "He+"}
_ELASTIC_GASES = ("Ne", "Ar", "Kr", "F2")       # the gases loki_tier drains


def bucket_of(proc) -> str | None:
    """Engine bucket for an LXCat process (mirrors loki_tier._bucket)."""
    g = getattr(proc, "gas", "")
    kind = proc.kind
    prod = proc.product
    thr = float(proc.threshold)
    if kind == "effective":
        return None
    if kind == "elastic":
        return f"el:{g}" if g in _ELASTIC_GASES else None
    if kind == "ionization":
        return _ION.get(g)
    if kind == "attachment":
        return "F-" if g == "F2" else None
    # excitation
    if g == "Ar":
        return "Ar*" if any(s in prod for s in ("Ar(3P", "Ar(1P1)")) else "Ar**"
    if g == "Kr":
        return "Kr*" if thr < 11.0 else "Kr**"
    if g == "Ne":
        return "Ne*" if thr < 17.5 else "Ne**"
    if g == "He":
        return "He*"
    if g == "F2":
        if "F2(V" in prod:
            return "vib"
        if "F2(a3)" in prod:
            return "F2(a3)"
        if any(s in prod for s in ("F2(A1)", "F2(C1)", "F2(H1)")):
            return "F2(A1)"
        return None
    return None


def maxwell_pdf(eps_eV, T_eV):
    """Maxwellian energy distribution f(eps) [1/eV], int f d eps = 1.
    ``T_eV`` scalar or (n,) -> (n_eps,) or (n, n_eps)."""
    eps = np.asarray(eps_eV, dtype=float)
    T = np.asarray(T_eV, dtype=float)
    if T.ndim == 0:
        return (2.0 / np.sqrt(np.pi)) * T ** -1.5 * np.sqrt(eps) * np.exp(-eps / T)
    return ((2.0 / np.sqrt(np.pi)) * T[:, None] ** -1.5 * np.sqrt(eps)[None, :]
            * np.exp(-eps[None, :] / T[:, None]))


def maxwell_rate(eps_eV, sigma_m2, T_eV) -> float:
    """[m^3/s] k(T) = int sigma v f_M d eps on the given grid (trapezoid).

    Closed forms (the gates): constant sigma -> sigma sqrt(8 e T / pi m_e);
    step sigma_0 above eps_th -> sigma_0 sqrt(8 e T / pi m_e)
    (1 + eps_th/T) exp(-eps_th/T)."""
    eps = np.asarray(eps_eV, dtype=float)
    v = np.sqrt(2.0 * eps * EV / M_E)
    g = np.asarray(sigma_m2, dtype=float) * v * maxwell_pdf(eps, float(T_eV))
    return float(_trapz(g, eps))


def bulk_temperature(te_eV, tail_fraction: float, tail_T_eV: float):
    """Bulk temperature of a two-Maxwellian with mean energy 1.5 Te:
    (1 - a) T_b + a T_t = Te. Raises if the stated tail already carries
    more than the whole mean energy (T_b would be <= 0)."""
    a = float(tail_fraction)
    if not 0.0 <= a < 1.0:
        raise ValueError(f"tail_fraction must be in [0, 1): {a}")
    te = np.asarray(te_eV, dtype=float)
    tb = (te - a * float(tail_T_eV)) / (1.0 - a)
    if np.any(tb <= 0.0):
        raise ValueError(
            f"two-Maxwellian tail ({a}, {tail_T_eV} eV) exceeds the mean "
            f"energy at Te = {np.min(te):.4g} eV (bulk temperature <= 0)")
    return tb


class MaxwellRates:
    """k(bucket, Te) integrated over a Maxwellian from the LXCat sets of a
    mixture; ``k_two_t`` adds a stated tail at fixed mean energy."""

    def __init__(self, mixture: dict):
        dens = {}
        self.skipped = []
        for g, x in mixture.items():
            if float(x) <= 0.0:
                continue
            try:
                gas_file(g)
            except KeyError:
                self.skipped.append(g)          # Xe: no sigma file (as beam.py)
                continue
            dens[g] = float(x)
        if not dens:
            raise ValueError("MaxwellRates: no gas of the mixture has a "
                             f"cross-section file (mixture {mixture})")
        self.cs = CrossSectionSet(dens)
        self.eps = EPS_GRID
        sig = self.cs.sigma_matrix(self.eps)             # (n_proc, n_eps)
        v = np.sqrt(2.0 * self.eps * EV / M_E)
        self._sv = {}
        self.processes_in = {}
        for i, proc in enumerate(self.cs.processes):
            b = bucket_of(proc)
            if b is None:
                continue
            self._sv[b] = self._sv.get(b, 0.0) + sig[i] * v
            self.processes_in.setdefault(b, []).append(
                f"{proc.gas}:{proc.kind}:{proc.product or '-'}@{proc.threshold:g}eV")
        self.buckets = sorted(self._sv)
        self.mixture = dict(mixture)

    def available_buckets(self):
        return set(self.buckets)

    def k(self, bucket: str, te_eV):
        """[m^3/s] Maxwellian at temperature Te (mean energy 1.5 Te);
        scalar or array Te, vectorized."""
        te = np.asarray(te_eV, dtype=float)
        scalar = te.ndim == 0
        te = np.atleast_1d(te)
        sv = self._sv.get(bucket)
        if sv is None:
            out = np.zeros_like(te)
            return float(out[0]) if scalar else out
        f = maxwell_pdf(self.eps, np.maximum(te, 1e-3))
        out = _trapz(sv[None, :] * f, self.eps, axis=1)
        return float(out[0]) if scalar else out

    def k_two_t(self, bucket: str, te_eV, tail_fraction: float,
                tail_T_eV: float):
        """[m^3/s] (1 - a) k(T_b) + a k(T_t) with the bulk T_b chosen so the
        mean energy stays 1.5 Te — the shape bracket at fixed <eps>."""
        a = float(tail_fraction)
        if a == 0.0:
            return self.k(bucket, te_eV)
        te = np.asarray(te_eV, dtype=float)
        scalar = te.ndim == 0
        tb = bulk_temperature(np.atleast_1d(te), a, tail_T_eV)
        kb = np.atleast_1d(self.k(bucket, tb))
        kt = float(self.k(bucket, float(tail_T_eV)))
        out = (1.0 - a) * kb + a * kt
        return float(out[0]) if scalar else out


class MaxwellShaped:
    """The engine tier: Maxwellian k(Te) for every bucket that has a sigma
    set, delegating everything else (E/N-native lookups, drift velocity,
    the table curves the validity gauge reads, run names) to the base
    LOKI tier — the same wrapper pattern as kinetic_tier.RelaxationBridged.

    ``tail=(fraction, T_eV)`` switches every Maxwell bucket to the
    two-Maxwellian at fixed mean energy (verification-tool use only).
    Below the base table's floor the tier follows the LOKI flat+taper rule
    (see ``k``); ``MaxwellRates`` itself is the bare integral at any Te.
    Buckets the base tabulates but no sigma set covers fall back to the
    base rate and are listed in ``fallback_buckets``."""

    def __init__(self, base, mixture: dict, tail=None):
        self.base = base
        self.mx = MaxwellRates(mixture)
        self.tail = tuple(tail) if tail is not None else None
        base_b = set(base.available_buckets())
        self.maxwell_buckets = sorted(base_b & self.mx.available_buckets())
        self.fallback_buckets = sorted(base_b - self.mx.available_buckets())

    def __getattr__(self, name):
        return getattr(self.base, name)

    def available_buckets(self):
        return self.base.available_buckets()

    #: below the base table's floor: the LOKI tier's own rule (flat for
    #: 0.35 ln-units of <eps>, then an exponential taper), natural-log units
    _TAPER = 0.35

    def floor_te(self, bucket: str):
        """Te [eV] of the base table's lowest tabulated <eps> for the bucket
        (None when the base has no curve for it)."""
        cur = self.base.curves.get(bucket) if hasattr(self.base, "curves") else None
        if cur is None:
            return None
        return float(np.exp(cur[0][0])) / 1.5

    def k(self, bucket: str, te_eV):
        """Maxwellian k(Te) inside the base table's <eps> range; BELOW the
        table floor the SAME flat+taper rule as loki_tier.LokiRates.k, so
        the LOKI/MAXWELL comparison isolates the in-table shape and both
        tiers extrapolate identically where neither is measured.

        Why the rule is needed (measured 2026-09-03, krf_beam_0d 0-D): the
        bare Maxwellian keeps steepening as exp(-eps_th/Te) into the
        afterglow, and scipy's finite-difference Jacobian then grinds —
        100 ns integrated in 7.7 s, the next 40 ns took 423 s with 9,198
        LU factorizations; the LOKI docstring's warning that a hard zero
        collapses the BDF step applies to a near-hard exponential too."""
        if bucket not in self.mx._sv:
            return self.base.k(bucket, te_eV)
        te = np.asarray(te_eV, dtype=float)
        scalar = te.ndim == 0
        te = np.atleast_1d(te)
        floor = self.floor_te(bucket)
        te_eval = te if floor is None else np.maximum(te, floor)
        if self.tail is None:
            out = np.atleast_1d(self.mx.k(bucket, te_eval)).astype(float)
        else:
            out = np.atleast_1d(self.mx.k_two_t(bucket, te_eval,
                                                self.tail[0], self.tail[1])).astype(float)
        if floor is not None:
            x = np.log(np.maximum(1.5 * te, 1e-3))
            x0 = np.log(1.5 * floor)
            lo = x <= x0
            if lo.any():
                dx = (x0 - self._TAPER) - x[lo]
                out[lo] *= np.where(dx <= 0.0, 1.0, np.exp(-dx / self._TAPER))
        return float(out[0]) if scalar else out

    def k_by_en(self, bucket: str, en_td: float) -> float:
        raise NotImplementedError(
            "the MAXWELL tier is Te-indexed only (a Maxwellian has no "
            "field -> EEDF relation): E/N-direct lookups belong to the LOKI "
            "tier; for the discharge lane use the rate-level comparison in "
            "tools/maxwell_verify.py (indexed by E/N through <eps>(E/N))")

    def describe(self) -> str:
        tail = ("single Maxwellian" if self.tail is None
                else f"two-Maxwellian tail {self.tail[0]:g} @ {self.tail[1]:g} eV")
        return (f"MAXWELL tier ({tail}) over base {type(self.base).__name__}: "
                f"{len(self.maxwell_buckets)} buckets from sigma, "
                f"{len(self.fallback_buckets)} fall back to base "
                f"{self.fallback_buckets}; skipped gases {self.mx.skipped}")
