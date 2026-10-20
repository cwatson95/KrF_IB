"""v0_6: does an EM wave that reheats the afterglow raise the 248-nm
efficiency? The twin pairs (2026-10-19).

A twin pair is one pump, LIGHT's given power for pulse_s, run twice:
without the wave (off) and with the wave on in a window after the pulse
(on), every other input the same. From the pair:

  E_wave    the wave's absorbed energy, the ledger's E_em [J/cm3]
  marginal  (E_248,on - E_248,off) / E_wave: the 248-nm energy each
            absorbed joule of the wave adds
  eff_off   E_248,off / E_given, the pump's own 248-nm efficiency
  eff_on    E_248,on / (E_given + E_wave), the wave counted as input
  d_eff     eff_on - eff_off, positive exactly where marginal > eff_off:
            the wave must turn its energy into 248 nm better than the
            pump does
  source    eff_off / marginal, how efficient (electricity to absorbed
            power) the wave's source must be, relative to the pump's, for
            the pair to break even at the wall plug (none where marginal
            <= 0)

and where the wave's energy went (each fate's change over E_wave, summing
to 1 where both ledgers close), the optics in the window (g0, alpha,
g0 - alpha: LIGHT's small-signal diagnostic, no cavity), the gas
temperature and the electrons.

The points (POINTS): LIGHT's reference point (1 atm, 0.5 % F2), and the
envelope notes' worked example with fluorine (IB/envelope_kr/
envelope_kr_v6.tex, 'A worked example', item 10: 0.5 atm, 2 % F2, n_e
1e15 cm-3, 300 kW/cm3, L = 10 cm), where the notes put the wave at
6.6 +- 2.5 THz with E0 = 47.6 kV/cm, the amplitude of a 3 MW/cm2 wave
(33.7 kV/cm rms, the field LIGHT takes). The reference point takes the
notes' pure-Kr frequency, 6.2 +- 0.6 THz. At both, LIGHT's pump: 300
kW/cm3 for 1 us, run to 2 us. The field runs at 1, 1/sqrt(10) and 1/10 of
the notes' (1, 1/10 and 1/100 of the power), the window over the whole
afterglow and over its first 100 ns.

The bands (network.BANDS) and the three constants with a quoted error
that the accuracy overview propagates (UNC_ROWS) run at their two ends,
one at a time, under v0_2's closure (two_t, seconds a run), as the
overview's propagated uncertainty does; each moves both twins. The
half-differences combine in quadrature and as an envelope. The kinetic
EEDF (minutes a run) gives the headline pairs and its difference from
two_t is a line of its own; --kinetic-bands runs the inputs' ends on the
kinetic headline pair too, where the two closures put the wave in
different regimes. A control twin with a 1 V/cm wave measures each
closure's numerical floor: the wave's power is nothing there, but the run
is cut and, under the kinetic EEDF, renewed as often as with the wave.

    python3 LIGHT/run_v0_6.py reheat --kinetic-bands --further   # the study as reported
    python3 LIGHT/run_v0_6.py reheat                   # the pairs and the two_t bands
    python3 LIGHT/run_v0_6.py reheat --fast            # two_t only, bands included
    python3 LIGHT/run_v0_6.py reheat --report          # the tables again from the saved runs

--kinetic-bands alone runs the kinetic bands at the reference point (about
50 runs of 10-20 min); --further adds FURTHER's 13.

The runs are saved as they finish (results/reheat_runs.json), so a study
that stops resumes where it stopped (a failed run is tried again), and
results/reheat.md holds the tables.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from . import optics
from .model import QE, Params, run
from .network import BANDS
from .report import VUV_BAND, accounting, row_dict, warnings

HERE = Path(__file__).resolve().parent
RUNS_JSON = HERE / "results" / "reheat_runs.json"
TABLES_MD = HERE / "results" / "reheat.md"
B248 = "248 nm (KrF B-X)"

# the envelope notes' wave (envelope_kr_v6.tex, 'A worked example', items
# 6, 9 and 10): E0 = 47.6 kV/cm, the amplitude of 3 MW/cm2. LIGHT's field
# is the rms one
E0_KV_CM = 47.6
E_RMS_V_CM = E0_KV_CM * 1e3 / math.sqrt(2.0)
POINTS = {
    "ref": {"label": "LIGHT's reference point (1 atm, 0.5 % F2, 6.2 THz)",
            "params": {}, "f_Hz": 6.2e12},
    "env": {"label": "the envelope's worked example (0.5 atm, 2 % F2, "
                     "6.6 THz)",
            "params": {"pressure_atm": 0.5, "x_f2": 0.02}, "f_Hz": 6.6e12},
}
FIELDS = (1.0, 1.0 / math.sqrt(10.0), 0.1)     # x E_RMS_V_CM
CONTROL_V_CM = 1.0
# the windows after the pulse, (open, close) [s] from pulse_s
WINDOWS = {"all": (0.0, math.inf), "early": (0.0, 100e-9)}
HEADLINE = (1.0, "all")
# the pairs the inputs' two ends run on: the notes' wave and a gentle one
# (1/100 of the power), both over the whole afterglow
BAND_GRID = (HEADLINE, (FIELDS[2], "all"))
# the kinetic EEDF's pairs (minutes each): the headline, the early window
# and the weaker fields over the whole afterglow
KINETIC_GRID = ((1.0, "all"), (1.0, "early"), (FIELDS[1], "all"),
                (FIELDS[2], "all"))
CLOSURES = ("kinetic", "two_t")
# further pairs, each a question the bands raised (2026-10-19): under the
# kinetic EEDF Kr2* + F at the top of its band (1e-9, no measurement)
# turned the reference point's reheat positive, so Kannari's value inside
# the band (3e-10) and five below it (where the sign turns), the gentle
# wave at the top, and the envelope point at the top and at the three
# inputs that move its two_t marginal most
FURTHER = (
    ("kinetic", "ref", "1", "all", "band:kr2_f=3e-10"),
    ("kinetic", "ref", "1", "all", "band:kr2_f=1e-10"),
    ("kinetic", "ref", "1", "all", "band:kr2_f=3e-11"),
    ("kinetic", "ref", "1", "all", "band:kr2_f=1e-11"),
    ("kinetic", "ref", "1", "all", "band:kr2_f=3e-12"),
    ("kinetic", "ref", "1", "all", "band:kr2_f=1e-12"),
    ("kinetic", "ref", "0.1", "all", "band:kr2_f=1e-09"),
    ("kinetic", "env", "1", "all", "band:kr2_f=1e-09"),
    ("kinetic", "env", "1", "all", "band:k_5s5p=0.5"),
    ("kinetic", "env", "1", "all", "band:k_5s5p=1.4"),
    ("kinetic", "env", "1", "all", "band:krf_f=7.8e-10"),
    ("kinetic", "env", "1", "all", "band:trimer_tg_pow=-0.5"),
)
# the search for the optimal reheat (2026-10-20), under the kinetic EEDF
# and at two values of the rate that decides the reference point: LIGHT's
# default Kr2* + F (0) and Kannari's 3e-10. One axis at a time from the
# reference point and the notes' wave over the whole afterglow: the field,
# the window, then the operating point (the notes' ranges, 0.1-2 atm and
# 100-500 kW/cm3, and F2 fills either side of the reference's burn-out).
# At omega >> nu the wave's heating depends on E / f alone (the field term
# goes as (E/omega)^2), so the frequency is not an axis of its own
OPT_KR2F = ("", "band:kr2_f=3e-10")
OPT_FIELDS = (FIELDS[1], 0.5623, 1.778, 3.162)        # x E_RMS_V_CM, and 1
OPT_WINDOWS = {"w30": (0.0, 30e-9), "w300": (0.0, 300e-9),
               "d30": (30e-9, math.inf), "d100": (100e-9, math.inf),
               "d300": (300e-9, math.inf)}
OPT_POINTS = {
    "ref_f2_0.25": {"label": "the reference point with 0.25 % F2",
                    "params": {"x_f2": 0.0025}, "f_Hz": 6.2e12},
    "ref_f2_1": {"label": "the reference point with 1 % F2",
                 "params": {"x_f2": 0.01}, "f_Hz": 6.2e12},
    "ref_p0.5": {"label": "the reference point at 0.5 atm",
                 "params": {"pressure_atm": 0.5}, "f_Hz": 6.2e12},
    "ref_p2": {"label": "the reference point at 2 atm",
               "params": {"pressure_atm": 2.0}, "f_Hz": 6.2e12},
    "ref_P100": {"label": "the reference point at 100 kW/cm3",
                 "params": {"power_W_cm3": 1e5}, "f_Hz": 6.2e12},
    "ref_P500": {"label": "the reference point at 500 kW/cm3",
                 "params": {"power_W_cm3": 5e5}, "f_Hz": 6.2e12},
}
# the combinations of the axes' best settings, added once the axes have run
OPT_COMBINED = ()
# every point's afterglow (t_end_s - pulse_s) [s]
AFTERGLOW_S = 1e-6
# the fixed field inside the plasma holds while the electrons stay below
# the envelope notes' design density (alpha L = 1 over L = 10 cm at 1e15
# cm-3): a wave that drives n_e past it re-ignites the afterglow (the
# notes' field x 1.78 grew n_e 12-fold in 300 ns, and the stiff solver
# stopped at 1.31 us), which a 0-D model with a given field cannot follow.
# The optimum leaves such pairs out
OPT_NE_LIMIT = 1e15
# the constants with a quoted error that the accuracy overview propagates
# beside the bands (docs/make_overview_v0_5.py UNC_ROWS): row id ->
# relative error (Shiu & Biondi 1977, Grossl et al. 1981, Xu et al. 1993)
UNC_ROWS = {"R11p": 0.10, "K15": 0.09, "C3": 0.09}
FATES = ("E248", "Evuv", "Eother", "heat", "stored", "walls")
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")


# ---------------------------------------------------------------------------
# the runs
# ---------------------------------------------------------------------------
def wave_tag(factor: float) -> str:
    return f"{factor:.4g}"


def inputs() -> list:
    """(kind, name, lo, hi, default): every band, then every row of
    UNC_ROWS at 1 -+ its quoted error."""
    out = [("band", n, b.lo, b.hi, b.default) for n, b in BANDS.items()]
    out += [("row", r, 1.0 - e, 1.0 + e, 1.0) for r, e in UNC_ROWS.items()]
    return out


def input_tag(kind, name, v) -> str:
    return f"{kind}:{name}={v!r}"


def _input_kw(tag: str) -> dict:
    """The Params fields that set one input: a band through `rates`, the
    trimer's temperature law through its own field, a row through
    `row_scale` (as the accuracy overview's _unc_kw)."""
    kind, _, rest = tag.partition(":")
    name, _, value = rest.partition("=")
    v = float(value)
    if kind == "row":
        return {"row_scale": ((name, v),)}
    if name == "trimer_tg_pow":
        return {"trimer_tg_pow": v}
    return {"rates": ((name, v),)}


def point_of(name) -> dict:
    """A point of POINTS or of the search's OPT_POINTS."""
    return POINTS[name] if name in POINTS else OPT_POINTS[name]


def window_of(name) -> tuple:
    """A window of WINDOWS or of the search's OPT_WINDOWS."""
    return WINDOWS[name] if name in WINDOWS else OPT_WINDOWS[name]


def params_of(spec) -> Params:
    """spec = (closure, point, wave, window, input): point a name of POINTS
    or OPT_POINTS, wave 'off', 'ctl' (the 1 V/cm control) or a factor on
    E_RMS_V_CM (wave_tag), window a name of WINDOWS or OPT_WINDOWS ('' when
    off), input '' or input_tag()."""
    closure, point, wave, window, inp = spec
    pt = point_of(point)
    kw = dict(eedf=closure, **pt["params"])
    if inp:
        kw.update(_input_kw(inp))
    p = Params(**kw)
    if wave != "off":
        a, b = window_of(window)
        e = CONTROL_V_CM if wave == "ctl" else float(wave) * E_RMS_V_CM
        p = p.replace(em_E_V_cm=e, em_f_Hz=pt["f_Hz"],
                      em_t_on_s=p.pulse_s + a, em_t_off_s=p.pulse_s + b)
    return p


def key(spec) -> str:
    return "|".join(spec)


def specs(closures=CLOSURES, points=tuple(POINTS), bands=True,
          kinetic_bands=(), kinetic_inputs=None, further=False) -> list:
    """Every run of the study: per closure and point the off twin, the
    pairs (every field and window under two_t, KINETIC_GRID under the
    kinetic EEDF) and the control; then, under two_t, the off twin and
    the pairs of BAND_GRID at each input's two ends (an end at the default
    is the base run and is not repeated); at the points named in
    kinetic_bands, the same for the headline pair under the kinetic EEDF
    (minutes a run), for the inputs named in kinetic_inputs (None: every
    one); with further, FURTHER's pairs."""
    out = []
    for closure in closures:
        for pt in points:
            out.append((closure, pt, "off", "", ""))
            grid = (KINETIC_GRID if closure == "kinetic" else
                    tuple((f, w) for f in FIELDS for w in WINDOWS))
            out += [(closure, pt, wave_tag(f), w, "") for f, w in grid]
            out.append((closure, pt, "ctl", "all", ""))
    if bands:
        for pt in points:
            for kind, name, lo, hi, default in inputs():
                for v in (lo, hi):
                    if v == default:
                        continue
                    tag = input_tag(kind, name, v)
                    out.append(("two_t", pt, "off", "", tag))
                    out += [("two_t", pt, wave_tag(f), w, tag)
                            for f, w in BAND_GRID]
    f, w = HEADLINE
    for pt in (p for p in kinetic_bands if p in points):
        for kind, name, lo, hi, default in inputs():
            if kinetic_inputs is not None and name not in kinetic_inputs:
                continue
            for v in (lo, hi):
                if v == default:
                    continue
                tag = input_tag(kind, name, v)
                out += [("kinetic", pt, "off", "", tag),
                        ("kinetic", pt, wave_tag(f), w, tag)]
    if further:
        out += [s for s in further_specs() if s[1] in points
                and s not in out]
    return out


def windows_of(p: Params) -> dict:
    """WINDOWS and the search's OPT_WINDOWS as times [s] for these Params,
    closed at t_end_s (records saved before 2026-10-20 hold WINDOWS only)."""
    return {name: (p.pulse_s + a, min(p.pulse_s + b, p.t_end_s))
            for name, (a, b) in {**WINDOWS, **OPT_WINDOWS}.items()}


def opt_specs() -> list:
    """The search for the optimal reheat (OPT_KR2F, OPT_FIELDS, OPT_WINDOWS,
    OPT_POINTS, OPT_COMBINED): every pair with its off twin, under the
    kinetic EEDF."""
    out = []
    for inp in OPT_KR2F:
        out.append(("kinetic", "ref", "off", "", inp))
        out += [("kinetic", "ref", wave_tag(f), "all", inp)
                for f in OPT_FIELDS + (1.0,)]
        out += [("kinetic", "ref", "1", w, inp)
                for w in ("early",) + tuple(OPT_WINDOWS)]
        for pt in OPT_POINTS:
            out += [("kinetic", pt, "off", "", inp),
                    ("kinetic", pt, "1", "all", inp)]
    for pt, wave, w, inp in OPT_COMBINED:
        out += [("kinetic", pt, "off", "", inp), ("kinetic", pt, wave, w, inp)]
    seen, uniq = set(), []
    for s in out:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq


def further_specs() -> list:
    """FURTHER's pairs, each with its off twin."""
    out = []
    for closure, pt, wave, w, inp in FURTHER:
        off = (closure, pt, "off", "", inp)
        if off not in out:
            out.append(off)
        out.append((closure, pt, wave, w, inp))
    return out


def _window_stats(res, t0, t1) -> dict:
    """The window's optics (optics.window), its electrons (n_e at the
    opening, int n_e dt, and the electron-weighted Te, (2/3) int U_e dt /
    int n_e dt, which the n_e floor cannot spoil) and the wave's energy
    absorbed in it [J/cm3]."""
    m = (res.t >= t0) & (res.t <= t1)
    if m.sum() < 2:
        raise ValueError(f"the window [{t0:.4g}, {t1:.4g}] s holds fewer "
                         "than two saved points")
    w = optics.window(res, t0, t1)
    t, ne, ue = res.t[m], res.ne[m], res.ue[m]
    ne_int = float(_trapz(ne, t))
    em = res.ledger["E_em"][m]
    return {"g0": w["g0"], "alpha": w["alpha"], "net": w["net"],
            "ne_open": float(ne[0]), "ne_max": float(ne.max()),
            "ne_int": ne_int,
            "te_e": (2.0 / 3.0) * float(_trapz(ue, t)) / ne_int
            if ne_int > 0.0 else float("nan"),
            "wave_in": float(em[-1] - em[0]) * QE}


def _kinetic_stats(res):
    """The kinetic EEDF with the wave alone (report.electrons_line): its
    solutions, those at LIGHT's mean, and the holding field's share of the
    wave's power (None for the other closures)."""
    if res.model.eedf != "kinetic":
        return None
    log = [e for e in getattr(res, "kinetic", []) if not e.get("skipped")]
    alone = [e for e in log if e.get("wave") and not e["on"]]
    sh = [e["dc_share"] for e in alone
          if np.isfinite(e.get("dc_share", float("nan")))]
    nan = float("nan")
    return {"solutions": len(log), "alone": len(alone),
            "reached": sum(1 for e in alone if e.get("reached")),
            "stepped": sum(1 for e in log if e.get("stepped")),
            "dc_median": float(np.median(sh)) if sh else nan,
            "dc_mean": float(np.mean(sh)) if sh else nan,
            "dc_max": float(max(sh)) if sh else nan}


def measure(res, windows=None) -> dict:
    """One run's record: its energies [J/cm3], every fate, the windows
    ({name: (t0, t1)} [s], by default WINDOWS after this run's pulse),
    the electrons and the warnings."""
    if windows is None:
        windows = windows_of(res.p)
    a = accounting(res)
    d = row_dict(res)
    bands = a["bands"]
    e248, evuv = bands.get(B248, 0.0), bands.get(VUV_BAND, 0.0)
    return {
        "given": a["e_given"] * QE, "wave": a["e_em"] * QE,
        "E248": e248 * QE, "Evuv": evuv * QE,
        "Eother": (sum(bands.values()) - e248 - evuv) * QE,
        "heat": (a["heat_elastic"] + a["heat_reactions"]) * QE,
        "stored": sum(a["stores"].values()) * QE,
        "walls": a["walls"] * QE, "closure": a["closure"],
        "tg_end": float(res.tg[-1]), "f2_left": d["f2_left_pct"],
        "krf_formed": a["krf"]["formed"], "row": d,
        "windows": {name: _window_stats(res, t0, t1)
                    for name, (t0, t1) in windows.items()},
        "kinetic": _kinetic_stats(res), "warnings": warnings(res)}


def _task(spec):
    t0 = time.time()
    rec = measure(run(params_of(spec)))
    rec["seconds"] = time.time() - t0
    return spec, rec


def load(path=RUNS_JSON) -> dict:
    """The saved runs, without the failed ones (they are tried again)."""
    path = Path(path)
    if not path.is_file():
        return {}
    cache = json.loads(path.read_text())
    return {k: v for k, v in cache.items() if "error" not in v}


def save(cache: dict, path=RUNS_JSON):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(cache, indent=0, sort_keys=True,
                              default=lambda o: o.item()
                              if hasattr(o, "item") else str(o)))
    tmp.replace(path)


def run_all(todo_specs, workers: int, path=RUNS_JSON) -> dict:
    """Run every spec not yet saved, the kinetic ones first (they are the
    long ones), each in a fresh process, saving each as it finishes."""
    cache = load(path)
    todo = [s for s in todo_specs if key(s) not in cache]
    # the longest first: kinetic before two_t, and the envelope point's
    # attaching afterglow (35-75 min a kinetic run) before the reference's
    todo.sort(key=lambda s: (s[0] != "kinetic", s[1] != "env"))
    if not todo:
        return cache
    # one thread a process: the processes are the parallelism
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    t0 = time.time()
    print(f"{len(todo)} runs to do ({sum(s[0] == 'kinetic' for s in todo)} "
          f"kinetic) on {workers} processes; saving to {path}", flush=True)
    with ProcessPoolExecutor(max_workers=workers,
                             mp_context=mp.get_context("spawn")) as ex:
        futs = {ex.submit(_task, s): s for s in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            spec = futs[fut]
            try:
                _, rec = fut.result()
                note = f"{rec['seconds']:.0f} s"
            except Exception as exc:          # a solver failure is a result
                rec = {"error": f"{type(exc).__name__}: {exc}"}
                note = rec["error"][:120]
            cache[key(spec)] = rec
            save(cache, path)
            print(f"  [{i}/{len(todo)}, {time.time() - t0:.0f} s] "
                  f"{key(spec)}: {note}", flush=True)
    return cache


# ---------------------------------------------------------------------------
# the pairs
# ---------------------------------------------------------------------------
def twin(off: dict, on: dict, window: str) -> dict:
    """The pair's numbers (the module's docstring), as fractions."""
    if not math.isclose(off["given"], on["given"], rel_tol=1e-9):
        raise ValueError("the twins must share their pump: given "
                         f"{off['given']} and {on['given']} J/cm3")
    ew = on["wave"]
    nan = float("nan")
    d248 = on["E248"] - off["E248"]
    eff_off = off["E248"] / off["given"]
    eff_on = on["E248"] / (on["given"] + ew)
    marginal = d248 / ew if ew > 0.0 else nan
    # a window the off twin's record does not hold (one saved before the
    # search's windows) leaves the window's numbers unknown
    blank = dict.fromkeys(("g0", "alpha", "net", "te_e", "ne_int"), nan)
    wo = off["windows"].get(window, blank)
    wn = on["windows"].get(window, blank)
    # the electrons the wave keeps in its window: the mean (int n_e dt over
    # the window, every point's afterglow 1 us) and the peak (records from
    # 2026-10-20 on)
    a, b = window_of(window) if window in WINDOWS or window in OPT_WINDOWS \
        else (nan, nan)
    span = min(b, AFTERGLOW_S) - a
    return {
        "ne_mean_on": wn["ne_int"] / span if span > 0.0 else nan,
        "ne_max_on": wn.get("ne_max", nan),
        "wave": ew, "wave_share": ew / on["given"], "d248": d248,
        "marginal": marginal, "eff_off": eff_off, "eff_on": eff_on,
        "d_eff": eff_on - eff_off,
        "source": (eff_off / marginal if marginal > 0.0 else math.inf)
        if ew > 0.0 else nan,
        "to": {f: (on[f] - off[f]) / ew if ew > 0.0 else nan
               for f in FATES},
        "d_g0": wn["g0"] - wo["g0"], "d_alpha": wn["alpha"] - wo["alpha"],
        "d_net": wn["net"] - wo["net"], "net_off": wo["net"],
        "d_tg_end": on["tg_end"] - off["tg_end"],
        "te_off": wo["te_e"], "te_on": wn["te_e"],
        "ne_int_off": wo["ne_int"], "ne_int_on": wn["ne_int"]}


def unc_stats(y0: float, pairs: dict) -> dict:
    """Each input's signed half-difference h = (y_hi - y_lo)/2, their
    quadrature sum, and the envelope: y0 plus the sum of every input's
    negative excursions, and y0 plus the sum of the positive ones (the
    accuracy overview's _unc_stats)."""
    h = {k: 0.5 * (b - a) for k, (a, b) in pairs.items()}
    quad = math.sqrt(sum(v * v for v in h.values()))
    lo = y0 + sum(min(a - y0, b - y0, 0.0) for a, b in pairs.values())
    hi = y0 + sum(max(a - y0, b - y0, 0.0) for a, b in pairs.values())
    return {"h": h, "quad": quad, "lo": lo, "hi": hi}


def _rec(cache, spec):
    """A run's record, None where it is missing or failed (run_all keeps a
    failed run's error in the cache it returns)."""
    r = cache.get(key(spec))
    return None if r is None or "error" in r else r


def pair(cache, closure, point, wave, window, inp="") -> dict:
    off = _rec(cache, (closure, point, "off", "", inp))
    on = _rec(cache, (closure, point, wave, window, inp))
    if off is None or on is None:
        return None
    return twin(off, on, window)


def band_pairs(cache, point, metric, f=HEADLINE[0], w=HEADLINE[1],
               closure="two_t") -> tuple:
    """The pair's `metric` (field f, window w) at its default and at each
    input's two ends: (y0, {name: (y_lo, y_hi)}). Inputs with a run
    missing are left out."""
    base = pair(cache, closure, point, wave_tag(f), w)
    if base is None:
        return None, {}
    out = {}
    for kind, name, lo, hi, default in inputs():
        ys = []
        for v in (lo, hi):
            tag = "" if v == default else input_tag(kind, name, v)
            t = pair(cache, closure, point, wave_tag(f), w, tag)
            ys.append(None if t is None else t[metric])
        if None not in ys:
            out[name] = tuple(ys)
    return base[metric], out


# ---------------------------------------------------------------------------
# the tables
# ---------------------------------------------------------------------------
def _pc(x, digits=2) -> str:
    if x is None or not np.isfinite(x):
        return "inf" if x == math.inf else "n/a"
    return f"{100.0 * x:.{digits}f}"


def _g(x, fmt=".3g") -> str:
    if x is None or not np.isfinite(x):
        return "n/a"
    return format(x, fmt)


def _src(x) -> str:
    """The break-even source efficiency: never where the wave adds no
    248 nm."""
    return "never" if x == math.inf else _g(x)


def _new_warnings(on: dict, off: dict) -> list:
    """The first sentences of the on twin's warnings that the off twin
    does not raise (compared with their numbers taken out)."""
    def norm(s):
        return re.sub(r"[-+]?\d[\d.,e+-]*", "#", s)
    seen = {norm(s) for s in off["warnings"]}
    return [s.split(". ")[0] for s in on["warnings"] if norm(s) not in seen]


def _wave_label(wave) -> str:
    if wave == "ctl":
        return f"control, {CONTROL_V_CM:g} V/cm"
    return f"{float(wave) * E_RMS_V_CM / 1e3:.3g} kV/cm"


def tables(cache) -> str:
    """The study's tables (Markdown)."""
    out = ["# LIGHT v0_6: the afterglow reheat, twin pairs",
           "",
           "Generated by `v0_6/reheat.py` from `results/reheat_runs.json`. "
           "Do not edit; rerun `python3 LIGHT/run_v0_6.py reheat --report`.",
           "",
           f"Wave: E0 = {E0_KV_CM} kV/cm amplitude ({E_RMS_V_CM / 1e3:.2f} "
           "kV/cm rms) at the envelope notes' frequency, on from the end "
           "of the pulse (1 us) to the end of the run (`all`, 1-2 us) or "
           "for 100 ns (`early`). Pump: 300 kW/cm3 for 1 us at both "
           "points. Percent of the stated input; d_eff in percentage "
           "points; source = eff_off / marginal, the wave source's "
           "electricity-to-absorbed efficiency, relative to the pump's, "
           "at which the pair breaks even.",
           ""]
    # 1. every pair
    out += ["## 1. The pairs", "",
            "| closure | point | field (rms) | window | E_wave, J/cm3 | "
            "E_wave / E_given, % | dE_248, J/cm3 | marginal, % | "
            "eff_off, % | eff_on, % | d_eff, pp | source | "
            "d(g0 - alpha) in the window, %/cm | dTg_end, K |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for closure in CLOSURES:
        for pt in POINTS:
            grid = (KINETIC_GRID if closure == "kinetic" else
                    tuple((f, w) for f in FIELDS for w in WINDOWS))
            rows = [(wave_tag(f), w) for f, w in grid] + [("ctl", "all")]
            for wave, w in rows:
                t = pair(cache, closure, pt, wave, w)
                if t is None:
                    raw = cache.get(key((closure, pt, wave, w, "")))
                    why = "failed" if raw and "error" in raw else "missing"
                    out.append(f"| {closure} | {pt} | {_wave_label(wave)} | "
                               f"{w} | {why} | | | | | | | | | |")
                    continue
                # the control's wave absorbs nothing: its ratios are noise
                ctl = wave == "ctl"
                out.append(
                    f"| {closure} | {pt} | {_wave_label(wave)} | {w} | "
                    f"{_g(t['wave'])} | {_pc(t['wave_share'], 3)} | "
                    f"{_g(t['d248'], '.3g')} | "
                    f"{'-' if ctl else _pc(t['marginal'])} | "
                    f"{_pc(t['eff_off'], 3)} | {_pc(t['eff_on'], 3)} | "
                    f"{100 * t['d_eff']:+.4f} | "
                    f"{'-' if ctl else _src(t['source'])} | "
                    f"{t['d_net']:+.4f} | {t['d_tg_end']:+.1f} |")
    out.append("")
    # 2. where the wave's energy went
    out += ["## 2. Where the absorbed wave went (window `all`)", "",
            "Change of each fate at t_end over E_wave, percent (sums to 100 "
            "where both ledgers close). Stored: species and electrons left "
            "at t_end (F atoms, ions, Kr*, Kr2*). Te: the electron-weighted "
            "temperature in the window.", "",
            "| closure | point | field (rms) | 248 nm | 147 nm VUV | "
            "other bands | heat | stored | walls | sum | Te off -> on, eV | "
            "int n_e dt off -> on, cm-3 s |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for closure in CLOSURES:
        for pt in POINTS:
            for f, w in BAND_GRID:
                t = pair(cache, closure, pt, wave_tag(f), w)
                if t is None:
                    continue
                to = t["to"]
                out.append(
                    f"| {closure} | {pt} | {_wave_label(wave_tag(f))} | "
                    + " | ".join(_pc(to[k], 1) for k in FATES)
                    + f" | {_pc(sum(to.values()), 1)} | "
                    f"{_g(t['te_off'], '.2f')} -> {_g(t['te_on'], '.2f')} | "
                    f"{_g(t['ne_int_off'])} -> {_g(t['ne_int_on'])} |")
    out.append("")
    # 3. the bands
    out += ["## 3. The bands", "",
            "Each input at its two ends, one at a time, moving both twins; "
            "the half-differences in quadrature and the envelope (the sum "
            "of every input's excursions each way), as the accuracy "
            "overview propagates them.", ""]
    sweeps = (("two_t", BAND_GRID), ("kinetic", (HEADLINE,)))
    for pt in POINTS:
        for closure, f, w in ((c, f, w) for c, g in sweeps for f, w in g):
            y0, mp_ = band_pairs(cache, pt, "marginal", f, w, closure)
            e0, ep = band_pairs(cache, pt, "d_eff", f, w, closure)
            if y0 is None or not mp_:
                continue
            sm, se = unc_stats(y0, mp_), unc_stats(e0, ep)
            out += [f"### {POINTS[pt]['label']}, "
                    f"{_wave_label(wave_tag(f))} rms, window `{w}`, "
                    f"{closure}", "",
                    "| input | lo -> hi | marginal at lo, hi, % | "
                    "d_eff at lo, hi, pp |", "|---|---|---|---|"]
            for kind, name, lo, hi, default in inputs():
                if name not in mp_ or name not in ep:
                    continue
                a, b = mp_[name]
                c, d = ep[name]
                out.append(f"| {name} | {lo:.3g} -> {hi:.3g} | {_pc(a)}, "
                           f"{_pc(b)} | {100 * c:+.4f}, {100 * d:+.4f} |")
            out += ["",
                    f"marginal {_pc(y0)} % +- {_pc(sm['quad'])} "
                    f"(quadrature), envelope {_pc(sm['lo'])} to "
                    f"{_pc(sm['hi'])} %; d_eff {100 * e0:+.4f} pp +- "
                    f"{100 * se['quad']:.4f}, envelope "
                    f"{100 * se['lo']:+.4f} to {100 * se['hi']:+.4f} pp "
                    f"({len(mp_)} of {len(inputs())} inputs).", ""]
            ck = pair(cache, "kinetic", pt, wave_tag(f), w)
            if closure == "two_t" and ck is not None:
                out += [f"The kinetic EEDF against two_t: marginal "
                        f"{_pc(ck['marginal'])} against {_pc(y0)} %, d_eff "
                        f"{100 * ck['d_eff']:+.4f} against "
                        f"{100 * e0:+.4f} pp.", ""]
        for closure in CLOSURES:
            c = pair(cache, closure, pt, "ctl", "all")
            if c is not None:
                out += [f"Control ({closure}, {POINTS[pt]['label']}): "
                        f"dE_248 {_g(c['d248'])} J/cm3, d_eff "
                        f"{100 * c['d_eff']:+.2e} pp.", ""]
    # 4. further pairs (FURTHER)
    rows = []
    for closure, pt, wave, w, inp in FURTHER:
        t = pair(cache, closure, pt, wave, w, inp)
        if t is not None:
            rows.append(f"| {closure} | {pt} | {inp.split(':')[1]} | "
                        f"{_wave_label(wave)} | {w} | {_g(t['wave'])} | "
                        f"{_pc(t['marginal'])} | {_pc(t['eff_off'], 3)} | "
                        f"{100 * t['d_eff']:+.4f} | {_src(t['source'])} | "
                        f"{_pc(t['to']['E248'], 1)} / "
                        f"{_pc(t['to']['Evuv'], 1)} / "
                        f"{_pc(t['to']['heat'], 1)} |")
    if rows:
        out += ["## 4. Further pairs", "",
                "Questions the bands raised: Kr2* + F -> KrF(B) (kr2_f, no "
                "measurement, 0 by default, band 0 to 1e-9 cm3/s) at its "
                "top turned the kinetic reference point's reheat positive. "
                "Kannari's 3e-10 lies inside the band; the values "
                "below it find where the sign turns.", "",
                "| closure | point | input | field (rms) | window | E_wave, "
                "J/cm3 | marginal, % | eff_off, % | d_eff, pp | source | "
                "wave to 248 / VUV / heat, % |",
                "|---|---|---|---|---|---|---|---|---|---|---|"] + rows + [""]
    # 5. the kinetic EEDF with the wave alone, and the warnings
    out += ["## 5. Checks", "",
            "| closure | point | field | window | closure residual off, on |"
            " kinetic solutions with the wave alone (at LIGHT's mean; "
            "time-stepped in the run) | holding field's share, median / "
            "mean / max | warnings only with the wave |",
            "|---|---|---|---|---|---|---|---|"]
    for closure in CLOSURES:
        for pt in POINTS:
            grid = (KINETIC_GRID if closure == "kinetic" else
                    tuple((f_, w_) for f_ in FIELDS for w_ in WINDOWS))
            off = _rec(cache, (closure, pt, "off", "", ""))
            for f_, w_ in grid:
                on = _rec(cache, (closure, pt, wave_tag(f_), w_, ""))
                if off is None or on is None:
                    continue
                k = on.get("kinetic") or {}
                new = _new_warnings(on, off)
                out.append(
                    f"| {closure} | {pt} | {_wave_label(wave_tag(f_))} | "
                    f"{w_} | {off['closure']:.1e}, {on['closure']:.1e} | "
                    + (f"{k['alone']} ({k['reached']}; "
                       f"{k.get('stepped', 0)}) | "
                       f"{_g(k['dc_median'], '.2g')} / "
                       f"{_g(k['dc_mean'], '.2g')} / {_g(k['dc_max'], '.2g')}"
                       if k else "- | -")
                    + f" | {'; '.join(new) or 'none'} |")
    out.append("")
    out += opt_tables(cache)
    return "\n".join(out)


def opt_pairs(cache, inp="") -> list:
    """The search's pairs found in the cache, at one Kr2* + F input:
    [(axis, point, wave, window, twin), ...], the reference pair first."""
    rows = [("reference", "ref", "1", "all")]
    rows += [("field", "ref", wave_tag(f), "all") for f in OPT_FIELDS]
    rows += [("window", "ref", "1", w)
             for w in ("early",) + tuple(OPT_WINDOWS)]
    rows += [("point", pt, "1", "all") for pt in OPT_POINTS]
    rows += [("combined", pt, wave, w) for pt, wave, w, i in OPT_COMBINED
             if i == inp]
    out = []
    for axis, pt, wave, w in rows:
        t = pair(cache, "kinetic", pt, wave, w, inp)
        if t is not None:
            out.append((axis, pt, wave, w, t))
    return out


def _within(t) -> bool:
    """The pair keeps its electrons below OPT_NE_LIMIT in the window (the
    peak where the record holds it, else the mean)."""
    ne = t["ne_max_on"] if np.isfinite(t["ne_max_on"]) else t["ne_mean_on"]
    return not (ne > OPT_NE_LIMIT)


def optimum(cache, inp=""):
    """The searched pair with the largest change in efficiency (the wave
    counted as input) at one Kr2* + F input, among those that keep the
    electrons below OPT_NE_LIMIT, or None."""
    rows = [r for r in opt_pairs(cache, inp)
            if np.isfinite(r[4]["d_eff"]) and _within(r[4])]
    return max(rows, key=lambda r: r[4]["d_eff"]) if rows else None


def _window_label(w) -> str:
    a, b = window_of(w)
    a_ns, b_ns = a * 1e9, b * 1e9
    end = "the end" if b == math.inf else f"+{b_ns:.0f} ns"
    return f"`{w}`: +{a_ns:.0f} ns to {end}"


def opt_tables(cache) -> list:
    """The search's section of the tables."""
    out = ["## 6. The search for the optimal reheat (kinetic EEDF)", "",
           "One axis at a time from the reference point and the notes' "
           "wave over the whole afterglow: the field, the window (from "
           "the pulse's end), then the operating point; then the axes' "
           "best settings combined. The best change in efficiency, the "
           "wave counted as input, is the optimum.", ""]
    for inp in OPT_KR2F:
        rows = opt_pairs(cache, inp)
        if not rows:
            continue
        title = ("LIGHT's default Kr2* + F (0)" if not inp else
                 f"Kr2* + F at {inp.split('=')[1]} cm3/s (Kannari)")
        out += [f"### {title}", "",
                "| axis | point | field (rms) | window | E_wave, J/cm3 | "
                "marginal, % | eff_off, % | d_eff, pp | source | "
                "dTg_end, K | Te in the window off -> on, eV | n_e in the "
                "window with the wave, mean / peak, cm-3 |",
                "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for axis, pt, wave, w, t in rows:
            out.append(
                f"| {axis} | {pt} | {_wave_label(wave)} | "
                f"{_window_label(w)} | {_g(t['wave'])} | "
                f"{_pc(t['marginal'])} | {_pc(t['eff_off'], 3)} | "
                f"{100 * t['d_eff']:+.4f} | {_src(t['source'])} | "
                f"{t['d_tg_end']:+.0f} | {_g(t['te_off'], '.2f')} -> "
                f"{_g(t['te_on'], '.2f')} | {_g(t['ne_mean_on'])} / "
                f"{_g(t['ne_max_on'])}"
                f"{'' if _within(t) else ' (past the limit)'} |")
        done = {(pt, wave, w) for _, pt, wave, w, _ in rows}
        lost = [key(s) for s in opt_specs() if s[4] == inp
                and s[2] != "off" and (s[1], s[2], s[3]) not in done]
        if lost:
            out += ["", "Not in the table (not run, or failed): "
                    + ", ".join(f"`{k}`" for k in lost) + "."]
        best = optimum(cache, inp)
        if best is None:
            out.append("")
            continue
        axis, pt, wave, w, t = best
        out += ["", f"Optimum (n_e below {OPT_NE_LIMIT:.0e} cm-3): "
                f"{point_of(pt)['label']}, "
                f"{_wave_label(wave)} rms at "
                f"{point_of(pt)['f_Hz'] / 1e12:g} THz, window "
                f"{_window_label(w)}: d_eff {100 * t['d_eff']:+.4f} pp "
                f"(marginal {_pc(t['marginal'])} % against the pump's "
                f"{_pc(t['eff_off'], 3)} %; source "
                f"{_src(t['source'])}).", ""]
    return out


def twin_text(off: dict, on: dict, window: str) -> str:
    """One pair in words, for run.py --twin."""
    t = twin(off, on, window)
    return "\n".join([
        "TWIN PAIR (the same pump without and with the wave)",
        f"  E_given {off['given']:.4g} J/cm3;  E_wave {t['wave']:.4g} J/cm3 "
        f"({_pc(t['wave_share'], 3)} % of E_given)",
        f"  E_248 off {off['E248']:.5g}, on {on['E248']:.5g} J/cm3;  "
        f"marginal {_pc(t['marginal'])} % of the wave",
        f"  efficiency off {_pc(t['eff_off'], 3)} %, on {_pc(t['eff_on'], 3)}"
        f" % of all input;  change {100 * t['d_eff']:+.4f} pp",
        f"  the wave's source must be {_g(t['source'])} x as efficient as "
        "the pump's to break even" if np.isfinite(t["source"]) else
        "  no source efficiency breaks even (the wave adds no 248 nm)"
        if t["source"] == math.inf else "  the wave absorbed nothing",
        "  where the wave went: " + ", ".join(
            f"{k} {_pc(v, 1)} %" for k, v in t["to"].items()),
        f"  in the window: g0 - alpha {t['net_off']:.4f} -> "
        f"{t['net_off'] + t['d_net']:.4f} %/cm;  Te {_g(t['te_off'], '.3f')}"
        f" -> {_g(t['te_on'], '.3f')} eV (electron-weighted)"])


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="run_v0_6.py reheat",
                                 description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--fast", action="store_true",
                    help="two_t only (the bands included): minutes")
    ap.add_argument("--no-bands", action="store_true",
                    help="leave out the inputs' two ends")
    ap.add_argument("--kinetic-bands", nargs="?", const="ref", default="",
                    metavar="POINTS",
                    help="the inputs' two ends on the headline pair under "
                    "the kinetic EEDF too, at these points (comma list; "
                    "alone: ref). About 50 runs a point, 10-60 min each")
    ap.add_argument("--kinetic-inputs", default=None,
                    help="with --kinetic-bands, only these inputs (comma "
                    "list of band or UNC_ROWS names)")
    ap.add_argument("--points", default=",".join(POINTS),
                    help="comma list of " + ", ".join(POINTS))
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--further", action="store_true",
                    help="and FURTHER's pairs (the questions the bands "
                    "raised; 13 kinetic runs)")
    ap.add_argument("--optimize", action="store_true",
                    help="and the search for the optimal reheat "
                    "(opt_specs: the field, the window, the operating "
                    "point; about 45 kinetic runs)")
    ap.add_argument("--report", action="store_true",
                    help="no runs: the tables from the saved runs")
    ap.add_argument("--runs", default=str(RUNS_JSON),
                    help="the saved runs (JSON)")
    ap.add_argument("--out", default=str(TABLES_MD), help="the tables")
    a = ap.parse_args(argv)
    points = tuple(x for x in a.points.split(",") if x)
    bad = set(points) - set(POINTS)
    if bad:
        raise SystemExit(f"--points: no point {sorted(bad)}")
    k_inputs = None
    if a.kinetic_inputs:
        k_inputs = tuple(x for x in a.kinetic_inputs.split(",") if x)
        bad = set(k_inputs) - {name for _, name, *_ in inputs()}
        if bad:
            raise SystemExit(f"--kinetic-inputs: no input {sorted(bad)}")
    k_points = tuple(x for x in a.kinetic_bands.split(",") if x)
    bad = set(k_points) - set(POINTS)
    if bad:
        raise SystemExit(f"--kinetic-bands: no point {sorted(bad)}")
    if a.report:
        cache = load(a.runs)
    else:
        todo = specs(closures=("two_t",) if a.fast else CLOSURES,
                     points=points, bands=not a.no_bands,
                     kinetic_bands=() if a.fast else k_points,
                     kinetic_inputs=k_inputs,
                     further=a.further and not a.fast)
        if a.optimize and not a.fast:
            todo += [s for s in opt_specs() if s not in todo]
        cache = run_all(todo, a.workers, a.runs)
    text = tables(cache)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text)
    print(text)
    print(f"\n(tables in {a.out})")
    return 0
