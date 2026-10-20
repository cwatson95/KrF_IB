"""LOKI tier — electron-impact rates from the curated LoKI-MC mixture tables.

ADAPTED from spectre's ``loki_electron.py`` (the audited PLAN (45)/(48)
reader) with two changes, both stated:

* the data ships WITH this package (``MODULES/loki_tables/``, provenance in
  its README.txt) so the engine runs standalone — the raw sweep directories
  stay with spectre;
* the mean electron energy per E/N comes from ``lookUpTableSwarm.txt``'s
  ``MeanE(eV)`` column instead of scraping per-job ``swarmParameters.txt``
  files (same numbers, one authoritative table), and the k(Te) lookup is
  vectorized over cell arrays.

THE DOCUMENTED APPROXIMATION, carried on the class: rates are tabulated
against the swarm mean energy <eps>(E/N) and looked up via <eps> = (3/2) Te
— a Maxwellian mapping of a non-Maxwellian EEDF, the standard fluid-model
compromise (same class as Chen's BOLSIG+ precomputation). This is exactly
the assumption the KINETIC tier (kinetic_tier.py) exists to SCORE (S4);
the difference is reported, never averaged away.

Low-energy behaviour (spectre's rule, kept with its reasoning): below the
table the rate holds flat for 0.35 decades then tapers exponentially — a
hard zero collapses the BDF step size, a bare flat clamp lets edge rates
cool Te without bound after the pump ends.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_NUM = r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"

# run-name sets per operating point (see loki_tables/README.txt)
RUN_NAMES_NE_FILL = ["V2point_NeKrF2", "V2point_NeKrF2_lowEN"]
RUN_NAMES_AR_FILL = ["V2point_ArKrF2", "V2point_ArKrF2_lowEN"]
RUN_NAMES_ARF_CHEN = ["ChenArF_NeArF2", "ChenArF_NeArF2_lowEN"]

#: sub-0.5 Td extensions measured 2026-09-04 (0.1.39). Appending one drops
#: the set's <eps> floor from 1.98 -> 0.761 eV (Ar fill) and 3.06 -> 0.891 eV
#: (Ne fill), i.e. into the beam-mode operating band, so the flat+taper rule
#: stops carrying the two channels it extrapolates the WRONG way (F-, vib —
#: the 0.1.38 sub-floor audit). OPT-IN: appending changes every beam anchor,
#: so the default table sets are unchanged pending the owner's call.
RUN_NAMES_AR_FILL_SUBFLOOR = ["V2point_ArKrF2_subfloor"]
RUN_NAMES_NE_FILL_SUBFLOOR = ["V2point_NeKrF2_subfloor"]


def _bucket(description: str):
    """Map a LoKI process description to an engine bucket.

    Vibrational rows are a Te-equation energy sink only ('vib'); elastic
    rows become per-gas momentum-transfer drains ('el:<gas>').

    BUG FIXED 2026-08-04 (found by the S4 KINETIC-vs-LOKI comparison, then
    pinned by a three-way check against the Maxwellian average of the raw
    F2 vib cross sections): spectre's mapper — inherited here — let ANY
    unmatched e+F2 row fall through to 'vib', which swallowed the
    F2 EFFECTIVE momentum-transfer row (~1.8e-13 m3/s, elastic-scale,
    x94-192 the true vib rate). The engine then drained 0.2 eV per FAKE
    event from the electron gas in every in-table-Te run. Effective rows
    are now DROPPED for every gas (the true Elastic row is tabulated
    separately — remapping Effective to el:<gas> would double-count), and
    'vib' matches only the explicit F2(V1..V4) excitations. The same
    defect exists upstream in projects/spectre loki_electron._bucket —
    flagged for its owner in this package's CHANGELOG, not edited there."""
    d = description
    if ",Effective" in d:
        return None
    if ",Elastic" in d:
        for g in ("Kr", "Ne", "F2", "Ar"):
            if d.startswith(f"e+{g}("):
                return f"el:{g}"
        return None
    if d.startswith("e+Kr"):
        if ",Ionization" in d:
            return "Kr+"
        return "Kr*" if "Kr(5s" in d else "Kr**"
    if d.startswith("e+Ar"):
        if ",Ionization" in d:
            return "Ar+"
        return "Ar*" if any(s in d for s in ("Ar(3P", "Ar(1P1)")) else "Ar**"
    if d.startswith("e+Ne"):
        if ",Ionization" in d:
            return "Ne+"
        return "Ne*" if "Ne(3s" in d else "Ne**"
    if d.startswith("e+F2"):
        if ",Attachment" in d:
            return "F-"
        if ",Ionization" in d:
            return "F2+"
        if "F2(a3)" in d:
            return "F2(a3)"
        if "F2(A1)" in d or "F2(C1)" in d or "F2(H1)" in d:
            return "F2(A1)"
        if "F2(V" in d.split("->", 1)[-1]:
            return "vib"               # the V1..V4 excitations, explicitly
        return None
    return None


def _parse_rate_table(path: Path):
    """-> (id -> description, {en_Td: {id: k_m3s}})."""
    desc, rows, header_ids = {}, {}, []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("#"):
            m = re.match(r"#\s*(\d+)\s+(\S+)", s)
            if m:
                desc[int(m.group(1))] = m.group(2)
            continue
        if not s:
            continue
        toks = s.split()
        if not re.match(rf"^{_NUM}$", toks[0]):
            header_ids = [int(m.group(1)) for t in toks[1:]
                          if (m := re.match(r"R(\d+)_ine", t))]
            continue
        vals = [float(t) for t in toks]
        if header_ids and vals[0] not in rows:
            rows[vals[0]] = {pid: v for pid, v in zip(header_ids, vals[1:])}
    return desc, rows


def _parse_swarm(path: Path):
    """-> {en_Td: (<eps>_eV, |v_drift|_m_s)} from lookUpTableSwarm.txt.

    Mean energy from MeanE(eV); the drift speed from BulkV_z(m/s) — the
    bulk drift is the transport quantity that carries conduction current
    (j = e n_e v_d). Magnitude taken: LoKI's sign convention follows its
    field axis; the engine handles polarity."""
    out = {}
    names = None
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        toks = s.split()
        if not re.match(rf"^{_NUM}$", toks[0]):
            names = toks
            continue
        if names is None:
            continue
        try:
            i_mean = names.index("MeanE(eV)")
            i_vd = names.index("BulkV_z(m/s)")
        except ValueError:
            raise ValueError(f"{path}: missing MeanE(eV)/BulkV_z in {names}")
        vals = [float(t) for t in toks]
        out.setdefault(vals[0], (vals[i_mean], abs(vals[i_vd])))
    return out


def _parse_mean_energies(path: Path):
    return {en: me for en, (me, _vd) in _parse_swarm(path).items()}


class LokiRates:
    """k(bucket, Te) from the packaged mixture tables, vectorized.

    Curves are log(k) vs log(<eps>), linearly interpolated, with the flat+
    taper extension below the table and a flat clamp above it."""

    #: taper width/offset in log-energy decades... (natural log units below)
    _TAPER = 0.35

    def __init__(self, run_names, root: Path = None):
        root = Path(root) if root is not None else _HERE / "loki_tables"
        # DUPLICATE guard (0.1.39). The cross-directory guard below catches
        # two DIFFERENT dirs tabulating one E/N, but the same dir listed
        # twice passed it (prev == d.name) and silently DOUBLED every
        # bucket — measured x2.0000 exactly. Latent since the guard landed
        # in 0.1.25; reachable now that callers concatenate run-name lists
        # (rates_for_mode's subfloor option).
        seen = [n for i, n in enumerate(run_names) if n in list(run_names)[:i]]
        if seen:
            raise ValueError(
                f"duplicate LOKI run names {sorted(set(seen))} — each "
                "directory's rates would be summed into the buckets twice "
                "(measured: exactly x2). Pass each table set once.")
        dirs = [root / n for n in run_names]
        missing = [str(d) for d in dirs if not d.is_dir()]
        if missing:
            raise FileNotFoundError(
                "LOKI tables not found: " + ", ".join(missing)
                + " (see MODULES/loki_tables/README.txt)")
        samples = {}
        en_samples = {}
        swarm = {}
        # cross-directory overlap guard (0.1.25, latent through 0.1.24):
        # summing process rates into a bucket is correct WITHIN one table
        # (a bucket is a sum over processes at one E/N), but a second
        # directory tabulating the SAME E/N would silently DOUBLE the
        # bucket. The shipped grids are disjoint (base + _lowEN), so this
        # never fired — the guard keeps it that way loudly.
        en_owner = {}
        for d in dirs:
            desc, rows = _parse_rate_table(d / "lookUpTableRateCoeff.txt")
            sw = _parse_swarm(d / "lookUpTableSwarm.txt")
            for en, rec in sw.items():
                swarm.setdefault(en, rec)
            for en, kk in rows.items():
                if en not in sw:
                    continue
                prev = en_owner.get(en)
                if prev is not None and prev != d.name:
                    raise ValueError(
                        f"E/N = {en:g} Td tabulated by BOTH {prev!r} and "
                        f"{d.name!r} — summing them would double-count "
                        "the bucket rates; the run-name sets must use "
                        "disjoint E/N grids")
                en_owner[en] = d.name
                eps = sw[en][0]
                for pid, k in kk.items():
                    b = _bucket(desc.get(pid, ""))
                    if b is None:
                        continue
                    samples.setdefault(b, {})
                    en_samples.setdefault(b, {})
                    if k > 0.0:
                        samples[b][eps] = samples[b].get(eps, 0.0) + k
                        en_samples[b][en] = en_samples[b].get(en, 0.0) + k
        self.curves = {}
        for b, agg in samples.items():
            xy = sorted(agg.items())
            if xy:
                self.curves[b] = (np.log([x for x, _ in xy]),
                                  np.log([y for _, y in xy]))
        # E/N-indexed curves for the circuit-driven discharge path — the
        # tables are NATIVELY E/N-indexed, so this lookup has no mapping
        # approximation at all (2026-08-04, the circuit-drive closure)
        self.en_curves = {}
        for b, agg in en_samples.items():
            xy = sorted(agg.items())
            if xy:
                self.en_curves[b] = (np.log([x for x, _ in xy]),
                                     np.log([y for _, y in xy]))
        sw_xy = sorted(swarm.items())
        self._sw_en = np.log([x for x, _ in sw_xy])
        self._sw_eps = np.array([r[0] for _, r in sw_xy])
        self._sw_vd = np.array([r[1] for _, r in sw_xy])
        self.en_range_td = (float(np.exp(self._sw_en[0])),
                            float(np.exp(self._sw_en[-1])))
        self.buckets = sorted(self.curves)
        self.run_names = list(run_names)

    def available_buckets(self):
        return set(self.buckets)

    def k(self, bucket: str, te_eV):
        """Rate coefficient [m^3/s] at Te (scalar or array), via <eps>=1.5 Te."""
        te = np.asarray(te_eV, dtype=float)
        scalar = te.ndim == 0
        te = np.atleast_1d(te)
        cur = self.curves.get(bucket)
        if cur is None:
            out = np.zeros_like(te)
            return float(out[0]) if scalar else out
        xs, ys = cur
        x = np.log(np.maximum(1.5 * te, 1e-3))
        out = np.exp(np.interp(x, xs, ys))
        # below-table: flat to xs[0]-TAPER, then exponential taper (spectre)
        lo = x <= xs[0]
        if lo.any():
            dx = (xs[0] - self._TAPER) - x[lo]
            out[lo] = np.exp(ys[0]) * np.where(dx <= 0.0, 1.0,
                                               np.exp(-dx / self._TAPER))
        return float(out[0]) if scalar else out

    # ------------------------------------------------ E/N-direct lookups
    def k_by_en(self, bucket: str, en_td: float) -> float:
        """[m^3/s] looked up DIRECTLY by reduced field [Td] — the native
        index of the tables; no <eps> = (3/2)Te mapping is involved. Below
        the table: the same flat+taper rule; above: flat clamp (the
        out-of-table fraction is a reported validity gauge)."""
        cur = self.en_curves.get(bucket)
        if cur is None:
            return 0.0
        xs, ys = cur
        x = float(np.log(max(en_td, 1e-3)))
        if x <= xs[0]:
            dx = (xs[0] - self._TAPER) - x
            return float(np.exp(ys[0]) * (1.0 if dx <= 0.0
                                          else np.exp(-dx / self._TAPER)))
        if x >= xs[-1]:
            return float(np.exp(ys[-1]))
        return float(np.exp(np.interp(x, xs, ys)))

    def drift_velocity(self, en_td: float) -> float:
        """|v_drift| [m/s] at E/N [Td] (BulkV from the swarm table),
        linear in log(E/N), linearly extrapolated to zero below the table
        (v_d -> 0 with the field), clamped above it."""
        if en_td <= 0.0:
            return 0.0
        lo_td = float(np.exp(self._sw_en[0]))
        if en_td < lo_td:
            return float(self._sw_vd[0]) * en_td / lo_td
        x = float(np.log(en_td))
        if x >= self._sw_en[-1]:
            return float(self._sw_vd[-1])
        return float(np.interp(x, self._sw_en, self._sw_vd))

    def mean_energy_at_en(self, en_td: float) -> float:
        """<eps> [eV] at E/N [Td] from the swarm table (consistency gauge
        for the fluid Te alongside the E/N-driven rates)."""
        x = float(np.log(max(en_td, 1e-3)))
        x = min(max(x, self._sw_en[0]), self._sw_en[-1])
        return float(np.interp(x, self._sw_en, self._sw_eps))


def rates_for_mode(mode: str, mixture: dict, root: Path = None,
                   subfloor: bool = False) -> LokiRates:
    """Pick the table set for a mode/mixture (Ar- vs Ne-buffered fill).

    ``subfloor`` appends the measured sub-0.5 Td extension (0.1.39) where
    one exists for the fill; the Chen/ArF set has none yet, so the flag is
    inert there and says so by returning the base set unchanged."""
    if mode == "arf_chen":
        return LokiRates(RUN_NAMES_ARF_CHEN, root)
    ar = float(mixture.get("Ar", 0.0))
    ne = float(mixture.get("Ne", 0.0))
    if ar >= ne:
        names = list(RUN_NAMES_AR_FILL)
        extra = RUN_NAMES_AR_FILL_SUBFLOOR
    else:
        names = list(RUN_NAMES_NE_FILL)
        extra = RUN_NAMES_NE_FILL_SUBFLOOR
    if subfloor:
        names += extra
    return LokiRates(names, root)
