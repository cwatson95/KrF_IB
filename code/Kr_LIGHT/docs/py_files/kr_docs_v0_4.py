"""What LIGHT_Kr v0_4's two documents share: the runs, cached and in
parallel, and the e-beam drive the overview scores the model under. v0_1's
to v0_3's documents keep their own (kr_docs.py, kr_docs_v0_2.py,
kr_docs_v0_3.py).

    from kr_docs_v0_4 import V, Params, run_many, beam_yields

run_many(jobs) runs each job once per version of the model's code: the
cache (docs/cache/runs_v0_4.json) is keyed by the md5 of every file under
v0_4/ that a run reads (the .py files and data/) and by the job's repr,
so a change to the model reruns everything and nothing else does. A job
is ("field", Params) for the model as it is, or ("beam", Params) for the
e-beam drive below, or ("xs", (Params, factors)) for the model on other
cross-section factors (a draw of the training's posterior). Each result is
summarize(res): row_dict's numbers and what the documents need besides
(the kinetic EEDF's log, the reservoirs, the events, the ledger, the
densities at the two ends and, v0_4, the kinetic EEDF's own tail ratio
along the run against the two-temperature closure's tables).

The e-beam drive is a test fixture of the documents, not part of the
model: LIGHT's prototype beam drive (the BeamFixture of the tests) with
pure-krypton yields. A fast electron's energy goes to ion pairs, to the
four 5s levels (v0_4: each its own species, at its own energy), to the
higher levels (the 5p rung Kr**, the network's C1 rows, the excess over
the 2p10 energy to heat) and to the secondaries' sub-excitation energy,
which the kinetic equation injects (linear on [0, E(Kr*(1s5))] with
their mean energy). The numbers per pair come from the
Fowler (Spencer-Fano) degradation of a fast electron on the same Biagi
8.97 cross sections the model uses, level by level (spencer_fano),
REBUILD's recursion: an excitation continues the electron at E - I, an
ionization shares E - I equally between the two electrons, elastic losses
are left out. The drive books no Hornbeck-Molnar ionization (the band
hm_kr at 0); BeamModel refuses a run that sets it.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIGHT_KR = HERE.parent
if str(LIGHT_KR) not in sys.path:
    sys.path.insert(0, str(LIGHT_KR))

from v0_4 import __version__ as V                      # noqa: E402
from v0_4 import xsec                                   # noqa: E402
from v0_4.model import QE, Model, Params, run          # noqa: E402
from v0_4.network import (KR_5S, LEVELS_5S, kr5s,      # noqa: E402
                          level_group, level_tag)
from v0_4.report import accounting, reservoir, row_dict  # noqa: E402
from v0_4.species import e_form                         # noqa: E402

CACHE = HERE / "cache" / "runs_v0_4.json"
RESERVOIRS = ("Kr*", "Kr**", "Kr2*")


# ---------------------------------------------------------------- the code
# v0_3: the swarm modules are the documents' own computation, read live by
# the generators; no run imports them, so the runs' cache does not key on
# them (the numbers files are stamped with full_md5, which does). v0_4:
# tail_field.py generates a data file the runs read (and key on), and no
# run imports it
NOT_RUN = ("swarm.py", "swarm_mc.py", "tail_field.py")


def code_md5() -> str:
    """md5 over every file of v0_4/ a run reads: its modules and data."""
    root = LIGHT_KR / "v0_4"
    h = hashlib.md5()
    files = sorted([f for f in root.glob("*.py") if f.name not in NOT_RUN]
                   + list((root / "data").glob("*")))
    for f in files:
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def full_md5() -> str:
    """md5 over every module and data file of v0_4/ (the numbers files'
    stamp)."""
    root = LIGHT_KR / "v0_4"
    h = hashlib.md5()
    files = sorted(list(root.glob("*.py")) + list((root / "data").glob("*")))
    for f in files:
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


# ----------------------------------------------------- the e-beam's yields
def spencer_fano(e0_eV: float = None, n_grid: int = 1200) -> dict:
    """The Fowler (Spencer-Fano) degradation of one fast electron in pure
    Kr on the vendored Biagi 8.97 cross sections: the expected number of
    each inelastic collision, on a uniform grid from 0 to e0 (55 % of the
    tables' end, REBUILD's rule, where the tables running out does not yet
    bend the count). Returns the counts per process and the energy left
    below the lowest threshold."""
    procs = [p for p in xsec.processes("Kr")
             if p.kind in ("excitation", "ionization")]
    ceiling = min(float(p.energy[-1]) for p in procs)
    e0 = 0.55 * ceiling if e0_eV is None else float(e0_eV)
    grid = np.linspace(0.0, e0, n_grid)
    du = grid[1] - grid[0]
    thr = np.array([p.threshold for p in procs])
    ion = np.array([p.kind == "ionization" for p in procs])
    U = np.zeros((n_grid, len(procs)))

    def at(e):
        if e <= 0.0:
            return np.zeros(len(procs))
        x = min(e / du, n_grid - 1.000001)
        i = int(x)
        f = x - i
        return U[i] * (1.0 - f) + U[i + 1] * f

    for i in range(1, n_grid):
        e = grid[i]
        sig = np.array([float(p.sigma_at(e)) for p in procs])
        sig[thr > e] = 0.0
        tot = float(sig.sum())
        if tot <= 0.0:
            continue
        acc = np.zeros(len(procs))
        for j in np.nonzero(sig)[0]:
            unit = np.zeros(len(procs))
            unit[j] = 1.0
            if ion[j]:
                acc += sig[j] / tot * (unit + 2.0 * at(0.5 * (e - thr[j])))
            else:
                acc += sig[j] / tot * (unit + at(e - thr[j]))
        U[i] = acc
    counts = {p.product: float(U[-1][j]) for j, p in enumerate(procs)}
    return {"e0": e0, "ceiling": ceiling, "counts": counts,
            "below": e0 - float(U[-1] @ thr),
            "thresholds": {p.product: p.threshold for p in procs}}


_YIELDS = None


def beam_yields() -> dict:
    """Per ion pair: Kr+ (1), Kr** (the levels above the 5s manifold, the
    network's C1 products), each 5s level (the G1 rows; v0_4: 's_<level>',
    and 's' their sum), the cascade heat (each level's excess over the
    energy of the species it makes: none for a 5s level since v0_4), the
    secondaries' sub-excitation energy and W."""
    global _YIELDS
    if _YIELDS is None:
        sf = spencer_fano()
        c, thr = sf["counts"], sf["thresholds"]
        levels = {p.product: level_group(p) for p in xsec.processes("Kr")
                  if p.kind == "excitation"}
        n_ion = c["Kr^+"]
        hi = [k for k, g in levels.items() if g != "G1"]
        lo = {level_tag(k): k for k, g in levels.items() if g == "G1"}
        if sorted(lo) != sorted(LEVELS_5S):
            raise AssertionError(f"the 5s levels: {sorted(lo)}")
        y_pp = sum(c[k] for k in hi) / n_ion
        y_lv = {lv: c[lo[lv]] / n_ion for lv in LEVELS_5S}
        y_s = sum(y_lv.values())
        cascade = (sum(c[k] * (thr[k] - e_form("Kr**")) for k in hi)
                   + sum(c[lo[lv]] * (thr[lo[lv]] - e_form(kr5s(lv)))
                         for lv in LEVELS_5S)) / n_ion
        w = sf["e0"] / n_ion
        eps = w - (e_form("Kr+") + y_pp * e_form("Kr**")
                   + sum(y * e_form(kr5s(lv)) for lv, y in y_lv.items())
                   + cascade)
        if eps <= 0.0:
            raise AssertionError(f"the beam's yields overspend: {eps}")
        _YIELDS = {"W": w, "ion": 1.0, "pp": y_pp, "s": y_s,
                   "cascade": cascade, "eps_sec": eps, "e0": sf["e0"],
                   "ceiling": sf["ceiling"], "nex": y_pp + y_s,
                   "ws": w / (1.0 + y_pp + y_s),
                   **{f"s_{lv}": y for lv, y in y_lv.items()}}
    return _YIELDS


class BeamModel(Model):
    """v0_4 under a uniform e-beam: the given power deposited as ion pairs
    with the yields above, the secondaries' energy into the electron gas
    and their spectrum into the kinetic equation. No field carries it (the
    kinetic log's field ratio is undefined). The yields book no
    Hornbeck-Molnar ionization, so the band hm_kr must stay at 0."""

    def __init__(self, p):
        if dict(p.rates).get("hm_kr", 0.0) != 0.0:
            raise ValueError("the beam's yields book no Hornbeck-Molnar "
                             "ionization: hm_kr must be 0")
        super().__init__(p)

    def rhs(self, t, y, on):
        if not on:
            return super().rhs(t, y, False)
        Y = beam_yields()
        dy = super().rhs(t, y, False)          # the chemistry, no field
        N0 = self.N0
        pairs = self.p.power_W_cm3 / QE / Y["W"]
        src = {"Kr": -(Y["ion"] + Y["pp"] + Y["s"]), "Kr+": Y["ion"],
               "Kr**": Y["pp"],
               **{kr5s(lv): Y[f"s_{lv}"] for lv in LEVELS_5S}}
        dn = np.zeros(self.ns)
        for s, v in src.items():
            dy[self.idx[s]] += pairs * v / N0
            dn[self.idx[s]] = pairs * v
        dy[self.i_ue] += pairs * Y["eps_sec"] / N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        tg_ev = y[self.i_tg] * 8.617333262e-5
        c_v = float((n * self.cv).sum()) * 8.617333262e-5
        dy[self.i_tg] += (pairs * Y["cascade"]
                          - tg_ev * float((dn * self.cv).sum())) / c_v
        dy[self.i_led + 0] += pairs * Y["W"] / N0          # E_dep
        dy[self.i_led + 2] += pairs * Y["cascade"] / N0    # cascade heat
        return dy

    def kinetic_injection(self, t, y, on):
        if not on:
            return ()
        Y = beam_yields()
        return ((self.p.power_W_cm3 / QE / Y["W"], e_form(kr5s("1s5")),
                 Y["eps_sec"]),)


# ------------------------------------------------------------- one run
_TAIL_ROWS = {}


def _tail_rows(m):
    """(rid, threshold, two-temperature table, factor) of the tail rows of
    a model like m, on the two-temperature closure (the field table's
    definition, v0_4/tail_field.py)."""
    key = (m.p.k76, tuple(m.p.rates), tuple(m.p.row_scale),
           repr(xsec.factors()))
    if key not in _TAIL_ROWS:
        two = Model(m.p.replace(eedf="two_t", capture="auto"))
        out = []
        for c in two.compiled:
            if c.kind == "xsec_tail" and c.factor != 0.0:
                p = xsec.find(*c.row.rate[1:4])
                out.append((c.row.rid, p.threshold, c.rate, c.factor))
        _TAIL_ROWS[key] = out
    return _TAIL_ROWS[key]


def tail_series(res) -> list:
    """v0_4: the kinetic EEDF's own tail ratio at each of its updates (the
    r at which the two-temperature EEDF at the run's Te removes the energy
    the kinetic one does through the tail rows), with the field table's
    and v0_3's beam table's r at the update's (n_e/N, Te): ((t, on, Te,
    n_e/N, r_kin, r_field, r_beam), ...). r_kin is nan where no r from
    1e-3 to 3 matches."""
    from scipy.optimize import brentq
    m = res.model
    rows = _tail_rows(m)
    field, beam = xsec.tail_field_table(), xsec.tail_ratio_table()
    n_tot = sum(res.n[s] for s in res.n)
    out = []
    for e in res.kinetic:
        if e.get("skipped") or "k" not in e or e["te"] < 0.3:
            continue
        te = e["te"]
        x = float(np.interp(e["t"], res.t, res.ne)
                  / np.interp(e["t"], res.t, n_tot))
        target = sum(e["k"].get(rid, 0.0) * thr for rid, thr, _, _ in rows)
        r = float("nan")
        if target > 0.0:
            def fn(lr):
                return math.log(sum(f * tab(te, math.exp(lr) * te) * thr
                                    for _, thr, tab, f in rows) / target)
            lo, hi = math.log(1e-3), math.log(3.0)
            if fn(lo) < 0.0 < fn(hi):
                r = math.exp(brentq(fn, lo, hi, xtol=1e-8))
        out.append((float(e["t"]), bool(e["on"]), float(te), x, r,
                    field(x, te), beam(x, 0.0)))
    return out


def summarize(res) -> dict:
    """The numbers the documents read from one run (JSON-safe)."""
    p, m = res.p, res.model
    a = accounting(res)
    out = {k: (v if isinstance(v, str) else float(v))
           for k, v in row_dict(res).items()}
    on = res.t <= p.pulse_s
    out.update({
        "e_dep_eV": float(a["e_dep"]),
        "heat_reactions_eV": float(a["heat_reactions"]),
        "heat_elastic_eV": float(a["heat_elastic"]),
        "bands_eV": {k: float(v) for k, v in a["bands"].items()},
        "stores_eV": {k: float(v) for k, v in a["stores"].items()},
        "events": {k: float(v[-1]) for k, v in res.events.items()},
        "ledger": {k: float(v[-1]) for k, v in res.ledger.items()},
        "n0": {k: float(v[0]) for k, v in res.n.items()},
        "n_end": {k: float(v[-1]) for k, v in res.n.items()},
        "n_max": {k: float(v.max()) for k, v in res.n.items()},
        # v0_4: the 5s levels' sum, and the excimer's two states'
        "n_max_krstar": float(sum(res.n[s] for s in KR_5S).max()),
        "n_max_kr2": float((res.n["Kr2*"] + res.n["Kr2*(0u+)"]).max()),
        "ue0": float(res.ue[0]), "ue_end": float(res.ue[-1]),
        "te_mid": float(res.te[int(res.t.searchsorted(0.5 * p.pulse_s))]),
        "te_max": float(res.te.max()), "tg_max": float(res.tg.max()),
        "ne_peak": float(res.ne.max()), "ne_mean_on": float(res.ne[on].mean()),
        "p_end_atm": float(res.pressure_atm[-1]), "nfev": int(res.nfev),
        "res": {name: reservoir(res, name) for name in RESERVOIRS},
        "eedf": m.eedf,
    })
    log = [e for e in getattr(res, "kinetic", []) if not e.get("skipped")]
    if log:
        on_log = [e for e in log if e["on"]]
        late = [e for e in on_log if e["t"] >= 0.5 * p.pulse_s]
        en = [e["en_td"] for e in on_log if e["en_td"] > 0.0]
        ratio = [e["power_ratio"] for e in late
                 if e["en_td"] > 0.0 and np.isfinite(e["power_ratio"])]
        bal = [abs(e["balance"]) for e in late if np.isfinite(e["balance"])]
        reach = [abs(e["te_kin"] / e["te"] - 1.0) for e in late]
        out["kin"] = {
            "solutions": len(log), "on": len(on_log),
            "en_min": min(en) if en else 0.0, "en_max": max(en) if en else 0.0,
            "ratio_min": min(ratio) if ratio else float("nan"),
            "ratio_max": max(ratio) if ratio else float("nan"),
            "ratio_mean": float(np.mean(ratio)) if ratio else float("nan"),
            "balance_max": max(bal) if bal else float("nan"),
            "reach_max": max(reach) if reach else float("nan"),
            "emax_min": min(e["emax"] for e in log),
            "emax_max": max(e["emax"] for e in log),
            "growth_off": sum(1 for e in log if e["growth_off"]),
            "unreached": sum(1 for e in on_log if not e["reached"]),
            # v0_3: the solver's cold restarts (boltzmann.Solver.eedf)
            "restarts": sum(e.get("restarts", 0) for e in log),
            # v0_4: its fallbacks to the growth term or a grid of its own
            "fallbacks": sum(e.get("fallbacks", 0) for e in log),
            # v0_4: the time-stepped solutions after the pulse
            # (Params.kin_afterglow_memory)
            "stepped_after": sum(1 for e in log if not e["on"]
                                 and e.get("memory")),
        }
        # v0_4: the kinetic EEDF's own tail ratio against the closure's
        # tables, over the pulse body (50 ns on) and as a series
        ser = tail_series(res)
        body = [s for s in ser if s[1] and s[0] >= 50e-9
                and np.isfinite(s[4])]
        if body:
            qf = np.array([s[5] / s[4] for s in body])
            qb = np.array([s[6] / s[4] for s in body])
            rk = np.array([s[4] for s in body])
            out["tail"] = {
                "n": len(body), "r_kin_lo": float(rk.min()),
                "r_kin_hi": float(rk.max()),
                "field_lo": float(qf.min()), "field_hi": float(qf.max()),
                "field_mean": float(qf.mean()),
                "beam_lo": float(qb.min()), "beam_hi": float(qb.max()),
                "beam_mean": float(qb.mean()),
                "series": ser[::max(1, len(ser) // 60)]}
    return out


def _job_key(job) -> str:
    kind, p = job
    return f"{kind}|{p!r}"


def _work(job):
    kind, p = job
    t0 = time.time()
    if kind == "field":
        res = run(p)
    elif kind == "beam":
        res = BeamModel(p).run()
    elif kind == "xs":
        # v0_3: the model on other cross-section factors (a posterior draw
        # of the training, or Biagi's file as published): p = (Params,
        # (m_knots, m, exc, ion)); summarized on the same factors
        params, f = p
        with xsec.using(xsec.Factors(*f)):
            res = run(params)
            out = summarize(res)
    else:
        raise ValueError(f"job kind {kind!r}")
    if kind != "xs":
        out = summarize(res)
    out["seconds"] = time.time() - t0
    if kind == "beam":
        out["beam"] = dict(beam_yields())
    return out


def _load() -> dict:
    try:
        return json.loads(CACHE.read_text())
    except (OSError, ValueError):
        return {}


def _save(cache: dict) -> None:
    CACHE.parent.mkdir(exist_ok=True)
    tmp = CACHE.with_suffix(".tmp")
    tmp.write_text(json.dumps(cache, indent=0, sort_keys=True,
                              allow_nan=True))
    os.replace(tmp, CACHE)


def run_many(jobs, workers: int = None, log=print) -> list:
    """Each job's summary, from the cache where this code ran it, else run
    (in parallel, `workers` processes) and cached as each finishes."""
    md5 = code_md5()
    cache = _load()
    if cache.get("md5") != md5:
        cache = {"md5": md5, "version": V, "runs": {}}
    keys = [_job_key(j) for j in jobs]
    todo = {}
    for k, j in zip(keys, jobs):
        if k not in cache["runs"] and k not in todo:
            todo[k] = j
    if todo:
        n = workers or min(len(todo), max(1, (os.cpu_count() or 2) - 6))
        log(f"running {len(todo)} of {len(jobs)} jobs on {n} processes")
        t0 = time.time()
        if n == 1:
            for k, j in todo.items():
                cache["runs"][k] = _work(j)
                _save(cache)
        else:
            with ProcessPoolExecutor(max_workers=n) as ex:
                futs = {ex.submit(_work, j): k for k, j in todo.items()}
                from concurrent.futures import as_completed
                for i, f in enumerate(as_completed(futs), 1):
                    k = futs[f]
                    try:
                        cache["runs"][k] = f.result()
                    except Exception:
                        # v0_3: name the failed job at once and cancel the
                        # ones not started; the pool still finishes the
                        # jobs in flight before the error leaves this block
                        # (2026-10-15: 20 minutes with nothing in the log)
                        log(f"  FAILED ({time.time() - t0:.0f} s): {k[:300]}")
                        for g in futs:
                            g.cancel()
                        raise
                    _save(cache)
                    log(f"  {i}/{len(todo)} done ({time.time() - t0:.0f} s,"
                        f" this one {cache['runs'][k]['seconds']:.0f} s): "
                        f"{k[:90]}")
    return [cache["runs"][k] for k in keys]


def cache_md5() -> str:
    """The md5 the cache was made with ('' when there is none)."""
    return _load().get("md5", "")


# ------------------------------------------------- the documents' runs
from v0_4.network import BANDS, bands_read, build   # noqa: E402

CLOSURES = ("auto", "two_t")      # the kinetic EEDF (default) and LIGHT v0_2's
# the measured constants with a quoted error that no band covers, as
# ((row, low factor, high factor), why): the overview propagates each
R11P_ERR = 0.10           # Shiu & Biondi 1977, +-10 %
K15_LO, K15_HI = 2.2e-31, 3.30e-31   # Grossl 1981 (2.2 +- 0.2) to Smith
K15_DEF = 2.4e-31                    # 1972 read at 300 K; the mean is 2.4
TAU_RAD, TAU_RAD_ERR = 270.0, 5.0    # Morikawa 1989, ns
CONSTANTS = (
    ("R11p", 1.0 - R11P_ERR, 1.0 + R11P_ERR),
    ("K15", K15_LO / K15_DEF, K15_HI / K15_DEF),
    ("KR2_RAD", TAU_RAD / (TAU_RAD + TAU_RAD_ERR),
     TAU_RAD / (TAU_RAD - TAU_RAD_ERR)),
)
# the convergence of the kinetic closure at the reference point: half the
# renewal step, a renewal at least every 5 ns, twice the energy cells,
# and the solver's tolerance 10x tighter
CONVERGENCE = (("kin_step", 0.05), ("kin_max_dt_s", 5e-9),
               ("kin_cells", 1200), ("rtol", 1e-9))
# the e-beam scores. Morozov et al. 2008: a continuous 8-15 keV beam of
# 5-50 mW into 500-1400 hPa, the efficiency flat over the beam power, so
# two low pump rates a decade apart at 1000 hPa. Eckstrom et al. 1988:
# pump rates up to 3e5 W/cm3 'at modest pressures', so that rate at 0.5,
# 1 and 2 atm, and at 1 atm with Kr* + Kr* at its band's low end
# and with Kr* + 2Kr at its band's top
MOROZOV_ATM = 1000.0 / 1013.25
MOROZOV_W = (30.0, 300.0)
ECKSTROM_W = 3.0e5
ECKSTROM_ATM = (0.5, 1.0, 2.0)
BEAM_PULSE, BEAM_TAIL_LO, BEAM_TAIL_HI = 1e-6, 50e-6, 20e-6


def ref_params(eedf="auto") -> Params:
    return Params(eedf=eedf)


def row_bands() -> list:
    """The bands the rows read (ib_scale acts on the EM wave alone and the
    reference point has none), in the order of network.BANDS."""
    read = bands_read(build())
    return [b for b in BANDS.values() if b.name in read]


def band_params(eedf="auto") -> list:
    """((band, end, value, Params), ...) for every end that is not the
    band's default."""
    out = []
    for b in row_bands():
        for end, v in (("low", b.lo), ("high", b.hi)):
            if v != b.default:
                out.append((b.name, end, v,
                            ref_params(eedf).replace(rates=((b.name, v),))))
    return out


def const_params(eedf="auto") -> list:
    out = []
    for rid, lo, hi in CONSTANTS:
        for end, f in (("low", lo), ("high", hi)):
            out.append((rid, end, f, ref_params(eedf).replace(
                row_scale=((rid, f),))))
    return out


def conv_params() -> list:
    return [(k, v, ref_params().replace(**{k: v})) for k, v in CONVERGENCE]


def switch_params(eedf="auto") -> dict:
    """Kannari's fit, the cycle off; v0_4: the two-temperature closure on
    v0_3's beam table, and the kinetic EEDF time-stepped after the pulse
    (the afterglow's bound)."""
    p = ref_params(eedf)
    out = {"kannari": p.replace(k76="kannari"),
           "nocycle": p.replace(rates=(("k_5s5p", 0.0),))}
    if eedf == "two_t":
        out["beam_tail"] = p.replace(tail="beam")
    else:
        out["afterglow"] = p.replace(kin_afterglow_memory=True)
    return out


# v0_4: what each of its changes does at the reference point, each undone
# alone from v0_4 (the 'undo' table): the formation at Ueda's 1.7e-32 for
# every level (v0_1 to v0_3's lumped default), the 3P1 forming at the 3P2's
# rate, the 5p relaxation's products all to 1s5, no redissociation of the
# excimer, of the ion, no transfer between the 5s levels, the J = 1 levels
# pumped by no electron
UEDA = 1.7e-32
REDISSOC = ("AK8r_1s5", "AK8r_1s4", "AK8r_1s3", "AK8fr_1s4")
MIXING = ("MIX_1s4", "MIX_1s5r", "MIX_1s2", "MIX_1s4r2")


def undo_params(eedf="auto") -> dict:
    p = ref_params(eedf)
    return {
        "ueda": p.replace(rates=(("kr2_form", UEDA), ("kr2_form_1s4", UEDA),
                                 ("kr2_form_1s3", UEDA))),
        "form_1s4": p.replace(rates=(("kr2_form_1s4",
                                      BANDS["kr2_form"].default),)),
        "kps_low": p.replace(rates=(("kps_5s", 0.0),)),
        "no_redissoc": p.replace(row_scale=tuple((r, 0.0)
                                                 for r in REDISSOC)),
        "no_k15r": p.replace(row_scale=(("K15r", 0.0),)),
        "no_mixing": p.replace(row_scale=tuple((r, 0.0) for r in MIXING)),
        "no_j1": p.replace(rates=(("k_5s5p_j1", 0.0),)),
    }


def beam_params(eedf="auto") -> dict:
    base = Params(eedf=eedf, ne_cm3=1e9, pulse_s=BEAM_PULSE, n_save=2001)
    out = {}
    for w in MOROZOV_W:
        out[("morozov", w)] = base.replace(
            pressure_atm=MOROZOV_ATM, power_W_cm3=w,
            t_end_s=BEAM_PULSE + BEAM_TAIL_LO)
    for atm in ECKSTROM_ATM:
        out[("eckstrom", atm)] = base.replace(
            pressure_atm=atm, power_W_cm3=ECKSTROM_W,
            t_end_s=BEAM_PULSE + BEAM_TAIL_HI)
    out[("eckstrom_pen", 1.0)] = out[("eckstrom", 1.0)].replace(
        rates=(("kr_pen", BANDS["kr_pen"].lo),))
    out[("eckstrom_form", 1.0)] = out[("eckstrom", 1.0)].replace(
        rates=(("kr2_form", BANDS["kr2_form"].hi),))
    return out


SCAN_ATM = (0.1, 0.25, 0.5, 1.0, 2.0)
SCAN_KW = (100.0, 300.0, 500.0)
K5S5P_TABLE = (0.0, 0.6, 1.0, 1.4)
# v0_4: the scan's corners (and its centre, the reference point) on the
# kinetic EEDF and on the closure with v0_3's beam table, which score the
# closure the scan uses
SCAN_CORNERS = ((0.25, 100.0), (0.25, 500.0), (2.0, 100.0), (2.0, 500.0))


def scan_params(eedf="two_t", tail="field", points=None) -> dict:
    pts = points or [(a, w) for a in SCAN_ATM for w in SCAN_KW]
    return {(a, w): ref_params(eedf).replace(pressure_atm=a,
                                             power_W_cm3=w * 1e3, tail=tail)
            for a, w in pts}


# v0_4: the photons each excitation makes in the swarm regime, for the
# electroluminescence score (Mano et al. at 1.1 bar and 296 K): a dilute
# seed of each 5s level and of the 5p rung, no drive and no electrons to
# speak of, to the end of its cascade; at the bands' defaults and with
# Kr2* + Kr at its band's top
EL_SEEDS = KR_5S + ("Kr**",)
EL_BAR, EL_TG, EL_SEED = 1.1, 296.0, 1e10


def el_seed_params(seed, rates=()) -> Params:
    return Params(eedf="two_t", power_W_cm3=0.0, pulse_s=0.0, t_end_s=40e-6,
                  ne_cm3=1.0, n_save=401, pressure_atm=EL_BAR / 1.01325,
                  tg0_K=EL_TG, init_extra_cm3=((seed, EL_SEED),),
                  rates=rates)


def el_params() -> dict:
    out = {}
    for s in EL_SEEDS:
        out[(s, "default")] = el_seed_params(s)
        out[(s, "quench")] = el_seed_params(
            s, (("kr2_quench", BANDS["kr2_quench"].hi),))
    return out


# v0_3: the trained cross sections' posterior, as draws of the training's
# MCMC (xsec.TRAINED_FILE's record), each a run of the reference point
XS_DRAWS = {"auto": 12, "two_t": 40}


def xs_factors(z) -> tuple:
    rec = json.loads((xsec.DATA_DIR / xsec.TRAINED_FILE[0]).read_text())
    return (tuple(rec["m_knots_eV"]), tuple(math.exp(v) for v in z[:6]),
            tuple(math.exp(v) for v in z[6:9]), 1.0)


def xs_jobs(eedf="auto") -> list:
    rec = json.loads((xsec.DATA_DIR / xsec.TRAINED_FILE[0]).read_text())
    draws = rec["mcmc"]["draws"][:XS_DRAWS[eedf]]
    return [("xs", (ref_params(eedf), xs_factors(z))) for z in draws]


def unit_factors() -> tuple:
    """Biagi's file as published (xsec.UNIT) as an 'xs' job's factors."""
    u = xsec.UNIT
    return (tuple(u.m_knots), tuple(u.m), tuple(u.exc), float(u.ion))


def biagi_jobs() -> list:
    """v0_4 at the reference point on Biagi's file as published: what the
    training does to v0_4's answer (v0_3's overview compared v0_2)."""
    return [("xs", (ref_params(eedf), unit_factors())) for eedf in CLOSURES]


def all_jobs() -> list:
    """Every run the two documents read, for one batch."""
    jobs = []
    for eedf in CLOSURES:
        jobs.append(("field", ref_params(eedf)))
        jobs += [("field", p) for *_, p in band_params(eedf)]
        jobs += [("field", p) for *_, p in const_params(eedf)]
        jobs += [("field", p) for p in switch_params(eedf).values()]
        jobs += [("beam", p) for p in beam_params(eedf).values()]
        jobs += [("field", ref_params(eedf).replace(rates=(("k_5s5p", k),)))
                 for k in K5S5P_TABLE if k not in (0.0, 1.0)
                 and ("k_5s5p", k) not in [(b, v) for b, _, v, _ in
                                           band_params(eedf)]]
        jobs += [("field", p) for p in undo_params(eedf).values()]
    jobs += [("field", p) for *_, p in conv_params()]
    jobs += [("field", p) for p in scan_params().values()]
    jobs += [("field", p) for p in scan_params(
        "auto", points=SCAN_CORNERS).values()]
    jobs += [("field", p) for p in scan_params(
        "two_t", "beam", points=SCAN_CORNERS + ((1.0, 300.0),)).values()]
    jobs += [("field", p) for p in el_params().values()]
    for eedf in CLOSURES:
        jobs += xs_jobs(eedf)
    jobs += biagi_jobs()
    seen, out = set(), []
    for j in jobs:
        k = _job_key(j)
        if k not in seen:
            seen.add(k)
            out.append(j)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="run (and cache) every run "
                                 "the two documents read")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--two-t-only", action="store_true")
    ap.add_argument("--list", action="store_true",
                    help="print the jobs and how many are cached")
    a = ap.parse_args()
    jobs = all_jobs()
    if a.list:
        cache = _load()
        have = cache.get("runs", {}) if cache.get("md5") == code_md5() else {}
        kin = sum(1 for j in jobs if (j[1][0] if j[0] == "xs" else j[1]).eedf
                  != "two_t")
        print(f"{len(jobs)} jobs ({kin} kinetic), "
              f"{sum(_job_key(j) in have for j in jobs)} cached")
        sys.exit(0)
    if a.two_t_only:
        jobs = [j for j in jobs
                if (j[1][0] if j[0] == "xs" else j[1]).eedf == "two_t"]
    t0 = time.time()
    run_many(jobs, workers=a.workers)
    print(f"{len(jobs)} jobs current ({time.time() - t0:.0f} s)")
