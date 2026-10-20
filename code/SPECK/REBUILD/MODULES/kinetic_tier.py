"""KINETIC tier — PIC_MCC_3D swarm exports, read as DATA (S4).

THE COUPLING RULE (ORESTES_CLASS_REBUILD.txt §5): PIC_MCC_3D exports data
with provenance; this engine reads it as one more tier. No import
dependency, no shared state — the architecture gate fails if any module
under ``MODULES/`` imports PIC_MCC_3D. The exporter that PRODUCES the file
is ``REBUILD/tools/export_kinetic.py``, which does import PIC_MCC_3D and is
therefore quarantined outside the engine package.

WHAT THE COMPARISON MEANS, precisely: both this tier and the LOKI tier are
kinetic Monte-Carlo swarm sources looked up by mean electron energy
(<eps> = (3/2) Te — the documented fluid-model mapping, which NEITHER tier
removes). The per-bucket ratio ``compare_tiers`` reports is therefore an
independent-code agreement band (PIC_MCC_3D's MCC lineage vs LoKI-MC), the
same class of statement as pic1d's x1.6 swarm band — a checkable
cross-validation, not a Maxwellian-vs-kinetic delta. It is REPORTED, never
averaged away, and the engine can run on either tier by configuration.

File format (.npz): ``eps_eV`` — mean-energy grid; ``k_<bucket>`` — rate
coefficient [m^3/s] per bucket on that grid (NaN where unmeasured);
``err_<bucket>`` — relative statistical error; ``meta`` — JSON string with
provenance (source package + commit, mixture, E/N points, particle count,
seed, date).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
DEFAULT_EXPORT = _HERE / "kinetic_tables" / "pic_mcc_3d_swarm.npz"


class KineticRates:
    """Same interface as loki_tier.LokiRates, backed by an export file."""

    def __init__(self, path):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(
                f"KINETIC export not found: {path}. Generate it with "
                "'python3 REBUILD/tools/export_kinetic.py' (imports "
                "PIC_MCC_3D — deliberately not part of the engine).")
        z = np.load(path, allow_pickle=False)
        self.meta = json.loads(str(z["meta"]))
        self.eps = np.asarray(z["eps_eV"], dtype=float)
        self.curves = {}
        self.errors = {}
        for key in z.files:
            if key.startswith("k_"):
                b = key[2:]
                if b.startswith("el_"):        # npz-safe name round-trip
                    b = "el:" + b[3:]
                k = np.asarray(z[key], dtype=float)
                good = np.isfinite(k) & (k > 0.0)
                if good.sum() >= 2:
                    self.curves[b] = (np.log(self.eps[good]),
                                      np.log(k[good]))
        for key in z.files:
            if key.startswith("err_"):
                b = key[4:]
                if b.startswith("el_"):
                    b = "el:" + b[3:]
                self.errors[b] = np.asarray(z[key], dtype=float)
        self.buckets = sorted(self.curves)
        # the SECOND W route: per-gas MC-cascade results (integral W at the
        # recorded E0, run through PIC_MCC_3D's actual MCC operator by the
        # exporter). Empty dict on older exports.
        self.w_cascade = dict(self.meta.get("w_cascade", {}))
        # the L3 measurement: field-free RELAXING-population rates binned
        # by instantaneous <eps> — the beam-mode EEDF class, reaching
        # BELOW the swarm tables (empty on older exports)
        self.relax = {}
        if "relax_eps_eV" in z.files:
            r_eps = np.asarray(z["relax_eps_eV"], dtype=float)
            for key in z.files:
                if key.startswith("relaxk_"):
                    b = key[7:]
                    if b.startswith("el_"):
                        b = "el:" + b[3:]
                    k = np.asarray(z[key], dtype=float)
                    e = np.asarray(z.get("relaxerr_" + key[7:],
                                         np.full_like(k, np.nan))
                                   if hasattr(z, "get") else
                                   z["relaxerr_" + key[7:]], dtype=float)
                    self.relax[b] = (r_eps, k, e)
        self.path = str(path)

    @classmethod
    def load_default(cls) -> "KineticRates":
        return cls(DEFAULT_EXPORT)

    def available_buckets(self):
        return set(self.buckets)

    def k(self, bucket: str, te_eV):
        """[m^3/s] at Te via <eps> = 1.5 Te (the shared mapping — see
        module docstring), log-log interpolated, clamped at the ends."""
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
        lo = x <= xs[0] - 0.35
        if lo.any():
            dx = (xs[0] - 0.35) - x[lo]
            out[lo] *= np.exp(-dx / 0.35)
        return float(out[0]) if scalar else out


class RelaxationBridged:
    """A tier whose SUB-TABLE region consumes the MEASURED relaxation
    curves — the L3 closure (2026-08-04).

    Beam-mode plasmas live at <eps> ~ 0.5–2 eV, BELOW every swarm table,
    where the base tiers extrapolate with the flat+taper stand-in. That
    stand-in was priced against PIC_MCC_3D's field-free relaxing-population
    measurement (0.1.8): elastic/vib within x0.8–1.9, F- attachment x3.3
    HIGH at 1.8 eV. This wrapper replaces the stand-in with the measurement
    itself wherever the measurement reaches:

    * at or above the base tier's table floor: the base lookup, untouched
      (in-table swarm data stays authoritative; the shared <eps> = (3/2) Te
      mapping there remains the documented S4-scored approximation);
    * below the floor, buckets with measured sub-floor bins: log-log
      interpolation on the relaxation curve, with a log-linear crossfade
      over the last ``BLEND`` (natural-log) energy units below the floor so
      k(Te) stays continuous for the BDF step controller — the band is
      narrower than the deepest measured bin's distance to the floor, so
      every measured bin is consumed pure;
    * below the measured range, and for buckets with no usable sub-floor
      bins: the base tier's own flat+taper rule (stated fallback — a hard
      zero collapses the BDF step size).

    Bins whose relative statistical error exceeds ``MAX_REL_ERR`` carry no
    information (a +-1 sigma band spanning more than x2.7) and are excluded
    — with the shipped export this drops the 100 %-error vib bin at 1.26 eV
    and every sub-floor Kr+ bin, leaving el:Ar/el:F2/el:Kr (down to
    1.26 eV), vib and F- (down to 1.76 eV): exactly the buckets whose price
    tags were measured. The circuit drive is exempt by construction — it
    looks rates up by E/N (`k_by_en`), which this wrapper delegates.

    Everything except ``k`` delegates to the base tier, so the engine's
    network build, regime line, validity gauges and circuit path see the
    base tier unchanged.
    """

    #: exclusion threshold on a bin's relative statistical error (stated
    #: rule, not a knob: err >= 0.5 means the +-1 sigma band spans > x2.7)
    MAX_REL_ERR = 0.5
    #: crossfade half-width at the table floor, natural-log energy units
    #: (same convention as the 0.35 flat/taper constant). 0.1 is narrower
    #: than the deepest measured bin's ln-distance to any floor (0.118 for
    #: the 1.76 eV bin under the 1.98 eV floor), so measured bins are
    #: consumed pure; the only blended region is unmeasured either way.
    BLEND = 0.1
    _TAPER = 0.35

    def __init__(self, base, relax_source):
        self.base = base
        self.relax_source_path = getattr(relax_source, "path", "")
        self._sub = {}
        base_curves = getattr(base, "curves", {})
        for b, (eps, k, err) in getattr(relax_source, "relax", {}).items():
            cur = base_curves.get(b)
            if cur is None:
                continue
            x_floor = float(cur[0][0])
            eps = np.asarray(eps, dtype=float)
            k = np.asarray(k, dtype=float)
            err = np.asarray(err, dtype=float)
            good = (np.isfinite(k) & (k > 0.0)
                    & np.isfinite(err) & (err < self.MAX_REL_ERR))
            if good.sum() < 2:
                continue
            xs = np.log(eps[good])
            ys = np.log(k[good])
            if xs[0] >= x_floor:        # nothing measured below the floor
                continue
            self._sub[b] = (xs, ys, x_floor)
        self.bridged_buckets = sorted(self._sub)

    def __getattr__(self, name):
        return getattr(self.base, name)

    def sub_table_coverage(self):
        """bucket -> (lowest measured <eps> [eV], base table floor [eV])."""
        return {b: (float(np.exp(xs[0])), float(np.exp(xf)))
                for b, (xs, _ys, xf) in self._sub.items()}

    def k(self, bucket: str, te_eV):
        cur = self._sub.get(bucket)
        if cur is None:
            return self.base.k(bucket, te_eV)
        te = np.asarray(te_eV, dtype=float)
        scalar = te.ndim == 0
        te = np.atleast_1d(te)
        base_out = np.asarray(self.base.k(bucket, te), dtype=float)
        xs, ys, x_floor = cur
        x = np.log(np.maximum(1.5 * te, 1e-3))
        below = x < x_floor
        if below.any():
            xb = x[below]
            r = np.interp(xb, xs, ys)
            lo = xb <= xs[0]
            if lo.any():
                dx = (xs[0] - self._TAPER) - xb[lo]
                r[lo] = ys[0] - np.where(dx <= 0.0, 0.0, dx / self._TAPER)
            w = np.clip((xb - (x_floor - self.BLEND)) / self.BLEND, 0.0, 1.0)
            ln_base = np.log(np.maximum(base_out[below], 1e-300))
            out = base_out.copy()
            out[below] = np.exp(w * ln_base + (1.0 - w) * r)
        else:
            out = base_out
        return float(out[0]) if scalar else out


def _eps_range(tier, bucket):
    cur = tier.curves.get(bucket)
    if cur is None:
        return None
    return float(np.exp(cur[0][0])), float(np.exp(cur[0][-1]))


def taper_vs_relaxation(kinetic, loki) -> list:
    """The L3 verdict rows: at beam-mode mean energies (BELOW the swarm
    tables) the engine's flat+taper extrapolation vs the MEASURED
    relaxing-population rates. Reported per bucket at each relax bin under
    the LOKI table floor; a ratio far from 1 quantifies exactly what the
    <eps>-mapping/taper stand-in costs in beam mode."""
    rows = []
    for b in sorted(kinetic.relax):
        cur = loki.curves.get(b)
        if cur is None:
            continue
        lo_eps = float(np.exp(cur[0][0]))
        eps, k, err = kinetic.relax[b]
        below = (eps < lo_eps) & (k > 0)
        if not below.any():
            continue
        for e_val, k_val, r_err in zip(eps[below], k[below],
                                       np.asarray(err)[below]):
            k_tap = float(loki.k(b, e_val / 1.5))
            ratio = k_tap / k_val if k_val > 0 else float("inf")
            rows.append((b, float(e_val), float(k_val), k_tap,
                         float(ratio), float(r_err)))
    return rows


def compare_tiers(kinetic, loki, te_eV=None) -> dict:
    """Per-bucket KINETIC/LOKI ratio — the S4 report.

    Compared ONLY over the overlap of the two tiers' tabulated mean-energy
    ranges: outside its table a tier returns a clamped extrapolation, and a
    ratio of clamp against table characterises the clamp, not the physics
    (found 2026-08-04: naive full-range ratios reached x1e7 on excitation
    buckets whose true-overlap agreement is x0.8–1.1). The discrepancy is
    characterised, never averaged into either tier."""
    out = {"buckets": {}}
    for b in sorted(set(kinetic.buckets) & set(loki.buckets)):
        rk = _eps_range(kinetic, b)
        rl = _eps_range(loki, b)
        if rk is None or rl is None:
            continue
        lo, hi = max(rk[0], rl[0]), min(rk[1], rl[1])
        if hi <= lo:
            continue
        eps = (np.geomspace(lo, hi, 15) if te_eV is None
               else 1.5 * np.asarray(te_eV))
        te = eps / 1.5
        kk = kinetic.k(b, te)
        kl = loki.k(b, te)
        good = (kk > 0) & (kl > 0)
        if not good.any():
            continue
        ratio = kk[good] / kl[good]
        out["buckets"][b] = {
            "ratio_median": float(np.median(ratio)),
            "ratio_min": float(ratio.min()),
            "ratio_max": float(ratio.max()),
            "n_points": int(good.sum()),
            "overlap_eps_eV": (lo, hi),
        }
    only_k = sorted(set(kinetic.buckets) - set(loki.buckets))
    only_l = sorted(set(loki.buckets) - set(kinetic.buckets))
    out["only_kinetic"] = only_k
    out["only_loki"] = only_l
    return out


def format_comparison(rep: dict) -> str:
    lines = ["KINETIC vs LOKI per-bucket ratio over the tabulated overlap "
             "(scored, never merged):"]
    for b, r in rep["buckets"].items():
        lo, hi = r["overlap_eps_eV"]
        lines.append(f"  {b:<8s} median x{r['ratio_median']:.2f}  "
                     f"range [{r['ratio_min']:.2f}, {r['ratio_max']:.2f}]  "
                     f"<eps> {lo:.1f}-{hi:.1f} eV")
    if rep["only_loki"]:
        lines.append("  LOKI-only buckets: " + ", ".join(rep["only_loki"]))
    if rep["only_kinetic"]:
        lines.append("  KINETIC-only buckets: "
                     + ", ".join(rep["only_kinetic"]))
    return "\n".join(lines)
