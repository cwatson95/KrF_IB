"""Numbers for v0_vs_v0_1_vs_v0_2.tex, LIGHT v0, v0_1 and v0_2 on the same
tests.

v0's and v0_1's accuracy numbers come from their own overviews' generated
files (docs/overview_*.tex for v0, docs/overview_v0_1_*.tex for v0_1). v0_2
has no overview, so this script makes v0_2's column of the same tests with
the same ML calls (scores). It writes review_numbers.tex, one
\\rvdef{key}{value} per number:
  - Tahoma under the same beam drive with each version's default network
    (v0's holds n_e at 1e15 cm^-3 and has no charge chemistry), with the
    Kr2* group added (v0, v0_1), with v0_2's unmeasured Kr2* rates at
    Kannari's estimates and at the top of the Kr2* + F band, and with each
    version's own power input
  - v0_2's column of the accuracy tests, the kinetic datasets, the pure-gas
    checks, Tahoma's g0, Isat and formation efficiency, the posterior of the
    B/C constants at Tahoma and the electron closure against the kinetic
    loop
  - the reference point for the three defaults, what each v0_1 change moves
    there, one change undone at a time, and v0_2's bands, each end alone
  - the notes' pressure x power grid for the three defaults, with v0_2 at
    the top of the Kr2* + F band and with the trimer's weak temperature law
ML fixes its LIGHT version when first imported, so everything ML computes
runs in a process of its own (--tahoma VERSION, --scores v0_2). Those
processes run side by side while this one runs the reference point and the
grid, which need no ML.

    python3 LIGHT/docs/REVIEW/make_review.py           # about 6 minutes
    python3 LIGHT/docs/REVIEW/make_review.py --check   # compare, write nothing
    cd LIGHT/docs/REVIEW && latexmk -pdf v0_vs_v0_1_vs_v0_2.tex
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_1.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_1.py").is_file())
ML_DIR = LIGHT.parent / "ML"
for _p in (str(LIGHT), str(ML_DIR), str(HERE.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

OUT = HERE / "review_numbers.tex"
VERSIONS = ("v0", "v0_1", "v0_2")
TAG = {"v0": "v0", "v0_1": "v01", "v0_2": "v02"}
PRESSURES = (0.1, 0.25, 0.5, 1.0, 2.0)          # atm, the notes' range
POWERS = (100.0, 300.0, 500.0)                   # kW/cm3
TRIMER_WEAK = -0.5                               # the pure-Kr T^-1/2 law
# v0_2's bands in the review's table, in LaTeX (the document's macros)
BAND_LABELS = {
    "kr2_f": r"\Krts{} + F \ra KrF(B) + Kr",
    "kr2_pen": r"\Krts{} + \Krts{} \ra \Krtp{} + 2Kr + \el{}",
    "kr2f_f": r"\KtF{} + F \ra 2Kr + 2F",
    "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr, 300\unit{K}",
    "kr2_f2": r"\Krts{} + \Ft{}, both channels",
    "kr_pen": r"\Krs{} + \Krs{} \ra \Krp{} + Kr + \el{}",
    "kr2f_f2": r"\KtF{} + \Ft{}",
    "trimer_tg_pow": r"the trimer's $(T_g/300)^n$",
}
UNITS = {"cm3/s": r"\unit{cm^3\,s^{-1}}", "cm6/s": r"\unit{cm^6\,s^{-1}}",
         "exponent": ""}


def _f2(x: float) -> str:
    return f"{x:.2f}"


def _kannari(band) -> float:
    """Kannari's estimate among a band's published values."""
    return next(v for v, who in band.marks if who.startswith("Kannari"))


# ---------------------------------------------------------------------------
# Tahoma, one version per process
# ---------------------------------------------------------------------------
def tahoma(version: str) -> dict:
    """The window's g0 (and T_e, n_e) at Tahoma for one LIGHT version:
    its default network under the beam, the same with the Kr2* group (v0,
    v0_1) or with v0_2's unmeasured Kr2* rates set (Kannari's estimates, the
    top of the Kr2* + F band, the trimer's weak law), and ML's
    configuration run on the version's own power input."""
    import ml
    if ml.LIGHT_VERSION != version:
        raise RuntimeError(f"ML runs {ml.LIGHT_VERSION}, asked for {version}")
    import numpy as np
    from ml import data, forward, params
    model = ml.light("model")
    bands = getattr(ml.light("network"), "BANDS", None)
    c = forward.with_rederived_k_e(params.light_now())
    drive = forward._drive()
    defaults = model.Params()

    def window(m):
        r = m.run()
        win = (r.t >= forward.WINDOW[0]) & (r.t <= forward.WINDOW[1])
        nb = float(np.mean(r.n["KrF(B)"][win]))
        return {"g0": forward.SIGMA_B * nb * 100.0,
                "te": float(np.mean(r.te[win])),
                "ne": float(np.mean(r.ne[win])),
                "tg": float(np.mean(r.tg[win]))}

    out = {}
    # v0_2's default network already holds the Kr2* route
    configs = (("def", False),) if bands else (("def", False),
                                                ("defkr2", True))
    for tag, kr2 in configs:
        p = forward.tahoma_params(kr2=kr2)
        p = p.replace(groups=tuple(g for g in p.groups if g == "kr2"),
                      ne_mode=defaults.ne_mode)
        if defaults.ne_mode == "clamped":
            # a held n_e is the model's assumption (v0's 1e15). A free n_e
            # keeps ML's seed, which only starts the run
            p = p.replace(ne_cm3=defaults.ne_cm3)
        out[tag] = window(forward.impose(drive.BeamModel(p), c))
        out[tag]["ne_mode"] = defaults.ne_mode
    out["own"] = window(forward.impose(model.Model(forward.tahoma_params()),
                                       c))
    if bands is None:
        # ML's configuration without Kr2* + F -> KrF(B) + Kr, Kannari's
        # estimate
        m = forward.impose(drive.BeamModel(forward.tahoma_params()), c)
        m.compiled = [cc for cc in m.compiled if cc.row.rid != "K2F_F"]
        out["noF"] = window(m)
    else:
        p = forward.tahoma_params()
        kannari = (("kr2_f", _kannari(bands["kr2_f"])),
                   ("kr2_pen", _kannari(bands["kr2_pen"])))
        for tag, kw in (("kannari", dict(rates=kannari)),
                        ("hiF", dict(rates=(("kr2_f", bands["kr2_f"].hi),))),
                        ("tri", dict(trimer_tg_pow=TRIMER_WEAK))):
            out[tag] = window(forward.impose(drive.BeamModel(p.replace(**kw)),
                                             c))
    meas = data.by_name("tahoma_g0")[0].values[0]
    for d in out.values():
        d["x"] = d["g0"] / meas
    return out


# ---------------------------------------------------------------------------
# v0_2's column of the accuracy tests, as the v0_1 overview makes v0_1's
# ---------------------------------------------------------------------------
def _krf_decay(model, kr_torr, f2_torr):
    """make_overview._krf_decay on a version's basic network: KrF(B) seeded
    in Kr/F2 at 300 K, no power, and the slow decay of KrF(B) + KrF(C)
    fitted."""
    import numpy as np
    tot = kr_torr + f2_torr
    p = model.Params(pressure_atm=tot / 760.0, x_f2=f2_torr / tot,
                     power_W_cm3=0.0, pulse_s=0.0, t_end_s=400e-9,
                     ne_mode="free", ne_cm3=1.0, n_save=2001,
                     init_extra_cm3=(("KrF(B)", 1e12),))
    r = model.Model(p).run()
    n = r.n["KrF(B)"] + r.n.get("KrF(C)", 0.0)
    rel = n / n[0]
    sel = (r.t > 20e-9) & (rel < 1e-3) & (rel > 1e-8)
    return -np.polyfit(r.t[sel], np.log(n[sel]), 1)[0]


def scores(version: str) -> dict:
    """The review's Table 1 for one version, made with the calls of
    make_overview_v0_1.py (fast_numbers and slow_numbers) and ML bound to
    the version. v0_2 has no overview, so its column comes from here. The
    keys are the overviews' own."""
    import ml
    if ml.LIGHT_VERSION != version:
        raise RuntimeError(f"ML runs {ml.LIGHT_VERSION}, asked for {version}")
    import make_overview as MO
    from make_system_map_basic import _k
    from ml import (data, eedf, fit, formation, forward, params, pure,
                    score)
    model = ml.light("model")
    n = {}
    # the photolysis decays against Xu 1993's printed constants
    for kr, meas in MO.XU_1993:
        n[f"xu.{kr:.0f}"] = MO._x(_krf_decay(model, kr, 1.0), meas)
    # ML's scorecard for the version's constants (k_e re-derived as ML does)
    c = forward.with_rederived_k_e(params.light_now())
    for ds in data.all_datasets():
        if ds.kind in data.DEVICE_KINDS:
            continue
        v = score.kinetic_score(ds, c)
        n[f"sc.{ds.name}"] = (MO._pct_signed(v) if len(ds.values) == 1
                              else f"{v:.1f}")
    # Tahoma under the beam drive, as ML scores the version
    o = formation.observables(c)
    meas = {d.name: d for d in data.by_role("device")}
    g0, eta, isat = (meas[k].values[0] for k in ("tahoma_g0",
                                                  "tahoma_formation",
                                                  "tahoma_isat"))
    n.update({"tah.g0": f"{o['g0']:.2f}", "tah.g0.x": MO._x(o["g0"], g0),
              "tah.isat": f"{o['isat']:.2f}",
              "tah.isat.x": MO._x(o["isat"], isat),
              "tah.eta": f"{o['eta']:.3f}", "tah.eta.x": MO._x(o["eta"], eta),
              "tah.tg": f"{o['tg']:.0f}"})
    # the pure-gas checks
    for g in ("Kr", "Ar"):
        x = pure.supply(g)
        n[f"pure.{g}.x"] = MO._x(x["per_eV_sf"], x["per_eV_meas"])
    rows, k15, _ = pure.kr_rates()
    ratios = [light / m for _, light, m in rows]
    n["pure.dr.lo"], n["pure.dr.hi"] = (f"{min(ratios):.2f}",
                                        f"{max(ratios):.2f}")
    n["pure.k15"] = f"${_k(k15)}$"
    n["pure.vuv"] = f"{pure.kr_vuv(3.3e2)['vuv']:.2f}"
    # the posterior of the B/C block, its Tahoma gain on this version (the
    # MCMC and the draws are seeded, as in the overviews)
    train = data.by_role("train")
    chain, _ = fit.mcmc(train, fit.fit(train), n_steps=2500)
    lo, mid, hi = score.device_posterior(chain, n_draws=40)
    n.update({"post.g0.lo": f"{lo:.2f}", "post.g0.mid": f"{mid:.2f}",
              "post.g0.hi": f"{hi:.2f}"})
    # the electron closure against the kinetic loop at Tahoma
    loop, _, _ = eedf.self_consistent(c)
    own = forward.run_device(c)
    for q in ("te", "ne"):
        n[f"el.dev.{q}"] = f"{100.0 * abs(own[q] / loop[q] - 1.0):.1f}"
    return n


def _start(flag: str, version: str):
    env = dict(os.environ, LIGHT_VERSION=version,
               PYTHONDONTWRITEBYTECODE="1")
    return subprocess.Popen([sys.executable, str(Path(__file__)), flag,
                             version], env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)


def _finish(proc) -> dict:
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"{' '.join(proc.args[-2:])} failed:\n{err}")
    return json.loads(out.strip().splitlines()[-1])


# ---------------------------------------------------------------------------
# the reference point and the grid (no ML, the three versions in this
# process)
# ---------------------------------------------------------------------------
def _v01_run(edit=None, energies=None, **kw):
    """row_dict of a v0_1 run with an edit to its compiled rows, or with
    formation energies swapped in for the run (restored after)."""
    import v0_1.model as M1
    import v0_1.species as S1
    from v0_1.report import row_dict
    saved = {}
    for name, e in (energies or {}).items():
        saved[name] = S1.SPECIES[name]
        S1.SPECIES[name] = replace(S1.SPECIES[name], e_form_eV=e)
    try:
        m = M1.Model(M1.Params(**kw))
        if edit is not None:
            edit(m)
        return row_dict(m.run())
    finally:
        S1.SPECIES.update(saved)


def _no_ba(m):
    """B's radiative decay all B-X, as in v0."""
    for c in m.compiled:
        if c.row.rid == "G4A":
            c.rate = 0.0
        elif c.row.rid == "G4":
            c.rate = 1.0 / m.p.tau_sp_s


def _no_ff(m):
    m.compiled = [c for c in m.compiled if c.row.rid != "FF"]


def _v0_ion_values(m):
    """v0's Kr+ conversion and Kr2+ recombination (REBUILD's copies)."""
    from v0.network import optional_rows
    own = {r.rid: r for r in optional_rows()}
    for c in m.compiled:
        if c.row.rid == "K15":
            c.rate = float(own["K15"].rate[1])
        elif c.row.rid == "R11p":
            kind, k = own["R11p"].rate
            c.kind, c.rate = kind, (k,)
            c.cap = lambda te: te if te > 0.02 else 1.5 * te


def _no_ion_ion_tg(m):
    for c in m.compiled:
        if c.row.rid in ("K26", "K108"):
            c.rate = tuple(c.rate[:2])


def _value(v: float, unit: str) -> str:
    from make_system_map_basic import _k
    if unit == "exponent":
        return f"${v:.2f}$"
    return "0" if v == 0.0 else f"${_k(v)}$"


def reference() -> tuple:
    """(numbers, table rows): the three defaults at the reference point,
    v0_1's changes undone one at a time, v0_2 with Kannari's two estimates,
    and v0_2's bands, each end alone (run_v0_2.py --bands)."""
    import v0.model as M0
    import v0_2.model as M2
    from v0.report import row_dict as rd0
    from v0.species import e_form as e0
    from v0_2.network import BANDS
    from v0_2.report import row_dict as rd2
    from v0_2.run import do_bands
    rows = {
        "v0": rd0(M0.run(M0.Params())),
        "v0groups": rd0(M0.run(M0.Params(ne_mode="free",
                                          groups=("ions", "f2e")))),
        "v01": _v01_run(),
        "noBA": _v01_run(_no_ba),
        "noFF": _v01_run(_no_ff),
        "v0ions": _v01_run(_v0_ion_values),
        "noIIT": _v01_run(_no_ion_ion_tg),
        "v0E": _v01_run(energies={s: e0(s) for s in ("KrF(B)", "KrF(C)",
                                                     "Kr2F*")}),
    }
    bands = do_bands(M2.Params())
    kannari = (("kr2_f", _kannari(BANDS["kr2_f"])),
               ("kr2_pen", _kannari(BANDS["kr2_pen"])))
    rows2 = {"v02": bands[0],
             "v02kannari": rd2(M2.run(M2.Params(rates=kannari))),
             "v02pen": rd2(M2.run(M2.Params(rates=kannari[1:])))}
    n = {}
    # v0 and v0_1 by default have no Kr2* route and so no 147-nm band
    for key, d in list(rows.items()) + list(rows2.items()):
        n[f"ref.{key}"] = _f2(d["emission_248"])
        n[f"ref.{key}.tg"] = f"{d['tg_end_K']:.0f}"
        n[f"ref.{key}.vuv"] = f"{d.get('emission_vuv', 0.0):.1f}"
    # the gap between the defaults, v0 minus v0_1, in two parts: v0's
    # default network against v0 with the charge chemistry, then the values
    # v0_1 changed. And v0_1 minus v0_2, the Kr2* route's measured rows
    e = {k: d["emission_248"] for k, d in rows.items()}
    n["ref.gap"] = _f2(e["v0"] - e["v01"])
    n["ref.d.network"] = _f2(e["v0"] - e["v0groups"])
    n["ref.d.values"] = _f2(e["v0groups"] - e["v01"])
    n["ref.gap12"] = _f2(e["v01"] - rows2["v02"]["emission_248"])
    n["ref.cool12"] = f"{rows['v01']['tg_end_K'] - rows2['v02']['tg_end_K']:.0f}"
    # v0_2's bands: each one's two ends, the rest at their defaults
    ref = bands[0]["emission_248"]
    by = {(r["band"], r["end"]): r for r in bands[1:]}
    lines, minor = [], 0.0
    for name, b in BANDS.items():
        lo, hi = by[(name, "low")], by[(name, "high")]
        span = max(abs(lo["emission_248"] - ref), abs(hi["emission_248"] - ref))
        n[f"band.{name}.lo"] = _f2(lo["emission_248"])
        n[f"band.{name}.hi"] = _f2(hi["emission_248"])
        if name not in ("kr2_f", "trimer_tg_pow"):
            minor = max(minor, span)
        unit = UNITS[b.unit]
        lines.append(
            f"{BAND_LABELS[name]} & {_value(b.default, b.unit)} & "
            f"{_value(b.lo, b.unit)} to {_value(b.hi, b.unit)}{unit} & "
            f"{lo['emission_248']:.2f} to {hi['emission_248']:.2f} & "
            f"{lo['emission_vuv']:.1f} to {hi['emission_vuv']:.1f} & "
            f"{lo['tg_end_K']:.0f} to {hi['tg_end_K']:.0f}\\\\")
    n["band.minor"] = _f2(minor)
    n["band.n"] = str(len(BANDS))
    return n, ["\\def\\rvBands{" + "\n".join(lines) + "}"]


def grid() -> tuple:
    """(numbers, table rows): the 248-nm share over the notes' ranges for
    the three defaults, and v0_2's two wide bands."""
    import v0.model as M0
    import v0_2.model as M2
    from v0.report import row_dict as rd0
    from v0_2.network import BANDS
    from v0_2.report import row_dict as rd2
    top = BANDS["kr2_f"].hi
    n, lines = {}, []
    d01, d12, spread_max, kr2_max, tri_max = 0.0, 0.0, 0.0, 1.0, 0.0
    n_band, hot = 0, []
    for p in PRESSURES:
        for pw in POWERS:
            kw = dict(pressure_atm=p, power_W_cm3=pw * 1e3)
            a = rd0(M0.run(M0.Params(**kw)))["emission_248"]
            b = _v01_run(**kw)["emission_248"]
            d = rd2(M2.run(M2.Params(**kw)))
            e, tg = d["emission_248"], d["tg_end_K"]
            f = rd2(M2.run(M2.Params(rates=(("kr2_f", top),),
                                     **kw)))["emission_248"]
            w = rd2(M2.run(M2.Params(trimer_tg_pow=TRIMER_WEAK,
                                     **kw)))["emission_248"]
            d01, d12 = max(d01, abs(a / b - 1.0)), max(d12, abs(e / b - 1.0))
            # how far the versions spread, against how far v0_2's two wide
            # bands reach, each relative to v0_2
            spread = (max(a, b, e) - min(a, b, e)) / e
            reach = max(f / e - 1.0, 1.0 - w / e)
            spread_max = max(spread_max, spread)
            kr2_max, tri_max = max(kr2_max, f / e), max(tri_max, 1.0 - w / e)
            if reach > spread:
                n_band += 1
            else:
                hot.append(tg)
            lines.append(f"{p:g} & {pw:.0f} & {a:.2f} & {b:.2f} & {e:.2f} & "
                         f"{f:.2f} & {w:.2f} & {tg:.0f}\\\\")
    n["grid.diff"] = f"{100 * d01:.0f}"
    n["grid.diff12"] = f"{100 * d12:.0f}"
    n["grid.spread"] = f"{100 * spread_max:.0f}"
    n["grid.kr2"] = f"{kr2_max:.1f}"
    n["grid.tri"] = f"{100 * tri_max:.0f}"
    n["grid.n"] = str(len(PRESSURES) * len(POWERS))
    n["grid.nband"] = str(n_band)
    n["grid.hotmin"] = f"{min(hot):.0f}" if hot else "none"
    return n, ["\\def\\rvGrid{" + "\n".join(lines) + "}"]


# ---------------------------------------------------------------------------
def numbers() -> tuple:
    import v0
    import v0_1
    import v0_2
    kids = {("tahoma", v): _start("--tahoma", v) for v in VERSIONS}
    kids[("scores", "v0_2")] = _start("--scores", "v0_2")
    ref, ref_tables = reference()
    g, grid_tables = grid()
    n = {"v0.version": v0.__version__, "v01.version": v0_1.__version__,
         "v02.version": v0_2.__version__}
    for v in VERSIONS:
        for cfg, d in _finish(kids[("tahoma", v)]).items():
            tag = f"tah.{TAG[v]}.{cfg}"
            n[f"{tag}.g0"] = _f2(d["g0"])
            n[f"{tag}.x"] = _f2(d["x"])
            n[f"{tag}.te"] = f"{d['te']:.1f}"
            n[f"{tag}.tg"] = f"{d['tg']:.0f}"
    for k, v in _finish(kids[("scores", "v0_2")]).items():
        n[f"v02.{k}"] = v
    n.update(ref)
    n.update(g)
    return n, ref_tables + grid_tables


def render(n: dict, tables) -> str:
    lines = [
        "% LIGHT v0, v0_1 and v0_2 comparison numbers. GENERATED by "
        "make_review.py",
        f"% from v0 ({n['v0.version']}), v0_1 ({n['v01.version']}) and v0_2 "
        f"({n['v02.version']}). Do not edit. Rerun the script instead.",
        r"\providecommand\rvdef[2]{\expandafter\def\csname rv@#1\endcsname{#2}}",
    ]
    lines += [f"\\rvdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--tahoma", choices=VERSIONS,
                    help="(internal) one version's Tahoma runs, as JSON")
    ap.add_argument("--scores", choices=VERSIONS,
                    help="(internal) one version's accuracy tests, as JSON")
    ap.add_argument("--check", action="store_true",
                    help="compare review_numbers.tex with a fresh run")
    a = ap.parse_args(argv)
    if a.tahoma:
        print(json.dumps(tahoma(a.tahoma)))
        return 0
    if a.scores:
        print(json.dumps(scores(a.scores)))
        return 0
    text = render(*numbers())
    if a.check:
        if OUT.read_text() != text:
            print(f"{OUT.name} is stale: rerun make_review.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(text)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
