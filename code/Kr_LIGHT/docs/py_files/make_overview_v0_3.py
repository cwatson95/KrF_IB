"""Numbers for 'How accurate LIGHT_Kr v0_3 is' (accuracy_overview_v0_3.tex).

Writes overview_v0_3_numbers.tex, one \\wgdef{key}{value} per number and a
few \\def tables, from the model:

  the swarm      v0_3's two-term solver run as a swarm (v0_3/swarm.py) on its
                 trained cross sections and on Biagi's file as published
                 (xsec.UNIT: v0_2's), as the experiments measure it
                 (v0_3/swarm_mc.py: LoKI-MC's bulk drift velocity, SST
                 ionization coefficient and transverse diffusion over the
                 two-term's), against every krypton dataset of ML's registry
                 (ML/ml/kr_data.py: Dutton's printed tables and LXCat's
                 LAPLACE database), each in the role ML's consistency check
                 gave it, with the training's record (the pinned
                 xsec.TRAINED_FILE: the factors, their MCMC and neural
                 estimates, the consistency check, the emulator's accuracy)
                 and the Monte Carlo's checks of the two-term solution.
  the e-beam     v0_3 under the beam drive (kr_docs_v0_3.BeamModel) against
                 Morozov et al. 2008 and Eckstrom et al. 1988, and the
                 drive's partition against Saito et al. 2003, on the trained
                 cross sections and (the partition) on v0_2's.
  uncertainty    every band the rows read at its two ends and the three
                 measured constants with a quoted error at theirs, one at a
                 time, on both closures (kr_docs_v0_3.band_params,
                 const_params), and the trained cross sections as draws of
                 the training's posterior (kr_docs_v0_3.xs_jobs).
  the reference  the kinetic EEDF and the two-temperature closure, the
                 5s -> 5p scale, Kannari's fit, the scan, the kinetic
                 closure's convergence (kr_docs_v0_3.conv_params), and v0_2
                 at the same points (kr_docs_v0_2, its own cache).

    python3 LIGHT_Kr/docs/make_overview_v0_3.py        # minutes (cached)
    cd LIGHT_Kr/docs && latexmk -pdf accuracy_overview_v0_3.tex

The runs come from kr_docs_v0_3.run_many (cached per version of the
model's code). --check compares the file with a fresh computation and
writes nothing.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kr_docs_v0_2 as K2                                  # noqa: E402
import kr_docs_v0_3 as K                                   # noqa: E402
import make_overview_v0_1 as OV1                           # noqa: E402
import make_overview_v0_2 as OV2                           # noqa: E402
from v0_3 import jung2006, swarm, swarm_mc, xsec           # noqa: E402
from v0_3.network import (BANDS, E_ASSOC_EV, KR_5P_TAU_S,  # noqa: E402
                          build, group_rows, levels)
from v0_3.species import e_form                            # noqa: E402

OUT = HERE / "overview_v0_3_numbers.tex"

# the e-beam's measurements and the scores' grid: v0_1's overview's, typed
# there from their sources (make_overview_v0_1.py)
MOROZOV, ECKSTROM = OV1.MOROZOV, OV1.ECKSTROM
SAITO_WS, SAITO_NEX = OV1.SAITO_WS, OV1.SAITO_NEX
JUNG_UNC = OV1.JUNG_UNC
CELLS_SWARM = OV1.CELLS_SWARM
# LoKI-MC on Biagi's file as published, 20 to 50 Td (v0_2's overview's
# run, typed there with its setup): E/N [Td] -> (the SST reduced Townsend
# coefficient [m^2], the mean energy [eV], the ionizations counted)
LOKI_MC = OV2.LOKI_MC
# Hornbeck & Molnar, Phys. Rev. 84, 621 (1951), abstract: the molecular
# ions' appearance potentials lie below the atomic ions' by "0.7-0.3+0.7
# volt in Kr" (the energy below Kr+'s at which Kr2+ appears, V)
HORNBECK_KR = (0.7, 0.3, 0.7)
_f, _pm, _sci, _rms = OV1._f, OV1._pm, OV1._sci, OV1._rms
WORDS = ("no", "one", "two", "three", "four", "five", "six", "seven",
         "eight", "nine")


# ------------------------------------------------------------- the swarm
# ML's registry of krypton's swarm data (ML/ml/kr_data.py), loaded on its
# own: it imports nothing of ML's, and importing the ml package would put
# LIGHT/ on sys.path, whose v0_2 and v0_3 are the KrF models
def _load_kr_data():
    path = HERE.parent.parent / "ML" / "ml" / "kr_data.py"
    spec = importlib.util.spec_from_file_location("ml_kr_data", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ml_kr_data"] = mod          # its dataclasses look it up
    spec.loader.exec_module(mod)
    return mod


KD = _load_kr_data()
TRAINED = json.loads((xsec.DATA_DIR / xsec.TRAINED_FILE[0]).read_text())
EN_SIGMA = 0.02          # E/N's uncertainty in the error model, as ML's
SETS = {"tr": None, "bi": xsec.UNIT}     # v0_3's trained set, Biagi's file
_SW: dict = {}


def swarm_point(en_td: float, tg: float, which: str = "tr",
                hm: float = 0.0, cells: int = CELLS_SWARM) -> dict:
    """v0_3's swarm (v0_3/swarm.py) on the trained set ('tr') or Biagi's
    file as published ('bi', v0_2's cross sections), with the measured
    quantities (swarm_mc: the bulk W, the SST alpha/N, the Monte Carlo's
    D_T/mu) beside the two-term ones."""
    key = (round(en_td, 10), tg, which, hm, cells)
    if key not in _SW:
        f = SETS[which]
        if f is None:
            d = swarm.point(en_td, tg, hm=hm, cells=cells)
        else:
            with xsec.using(f):
                d = swarm.point(en_td, tg, hm=hm, cells=cells)
        d = dict(d)
        d.update({k + "_m": v for k, v in swarm_mc.as_measured(d, en_td).items()})
        _SW[key] = d
    return _SW[key]


def _model(d, which):
    return np.array([swarm_point(float(x), d.tg, which)[d.kind + "_m"]
                     for x in d.x])


def _sigma_eff(d):
    """ML's error model: the dataset's sigma and E/N's 2 % through the
    model's slope (Biagi's file, a central difference in ln E/N)."""
    h = 0.01
    sl = []
    for x in d.x:
        lo = swarm_point(float(x) * math.exp(-h), d.tg, "bi")[d.kind]
        hi = swarm_point(float(x) * math.exp(h), d.tg, "bi")[d.kind]
        sl.append((math.log(hi) - math.log(lo)) / (2 * h))
    return np.sqrt(d.sigma ** 2 + (EN_SIGMA * np.array(sl)) ** 2)


def dataset_scores() -> list:
    """Every dataset: role, the model's error on both sets (rms and mean
    percent, chi2 per point with ML's error model)."""
    out = []
    train = set(TRAINED["train"])
    for d in KD.load():
        role = ("train" if d.name in train else
                "check" if d.name in KD.CHECK else
                "validation" if d.name in VALIDATION else "held_out")
        s = _sigma_eff(d)
        row = {"d": d, "role": role, "sig": s}
        for which in ("bi", "tr"):
            m = _model(d, which)
            r = m / d.y - 1.0
            row[which] = {"rms": 100 * _rms(r), "mean": 100 * float(np.mean(r)),
                          "chi": float(np.mean((np.log(m / d.y) / s) ** 2)),
                          "ratio": m / d.y}
        out.append(row)
    return out


def _pooled(rows, which, key="rms"):
    r = np.concatenate([x[which]["ratio"] for x in rows]) - 1.0
    return 100 * (_rms(r) if key == "rms" else float(np.mean(r)))


LABEL = {"w_pack62_300": "Pack et al.\\ 1962, 300\\unit{K}",
         "w_pack62_368": "Pack et al.\\ 1962, 368\\unit{K}",
         "w_pack62_195": "Pack et al.\\ 1962, 195\\unit{K}",
         "w_bowe60": "Bowe 1960", "w_english53": "English and Hanna 1953",
         "w_hunter88": "Hunter et al.\\ 1988",
         "w_kucuk81": "Kucukarpaci and Lucas 1981",
         "w_kucuk81_high": "the same, above 300\\unit{Td}",
         "w_nakamura90": "Nakamura 1990",
         "w_pack92_195": "Pack et al.\\ 1992, 195\\unit{K}",
         "w_pack92": "Pack et al.\\ 1992, 300\\unit{K}",
         "w_dallarmi92": "Dall'Armi et al.\\ 1992",
         "dt_alamin87": "Al-Amin and Lucas 1987",
         "dt_alamin87_high": "the same, above 300\\unit{Td}",
         "dt_koizumi86": "Koizumi et al.\\ 1986",
         "a_specht80_direct": "Specht et al.\\ 1980, direct",
         "a_specht80_total": "Specht et al.\\ 1980, total",
         "a_bhatta79": "Bhattacharya 1979",
         "a_bhatta79_high": "the same, above 300\\unit{Td}",
         "a_dallarmi92_low": "Dall'Armi et al.\\ 1992, below 50\\unit{Td}",
         "a_dallarmi92": "Dall'Armi et al.\\ 1992",
         "a_jacques86_low": "Jacques et al.\\ 1986, 40\\unit{Td}",
         "a_jacques86": "Jacques et al.\\ 1986",
         "a_kruithof40_low": "Kruithof 1940, below 50\\unit{Td}",
         "a_kruithof40": "Kruithof 1940",
         "a_kruithof40_high": "the same, above 300\\unit{Td}",
         "a_heylen71_low": "Heylen 1971, below 50\\unit{Td}",
         "a_heylen71": "Heylen 1971"}
CITE = {"Hunter 1988": "hunter88", "Pack 1992": "pack92",
        "Pack 1962": "pack62,dutton75", "Nakamura 1990": "nakamura90",
        "Dall'Armi 1992": "dallarmi92", "Bowe 1960": "bowe60",
        "English 1953": "english53", "Kucukarpaci 1981": "kucukarpaci81",
        "Koizumi 1986": "koizumi86", "Al-Amin 1987": "alamin87",
        "Jacques 1986": "jacques86", "Bhattacharya 1979": "bhattacharya79",
        "Specht 1980": "specht80", "Kruithof 1940": "kruithof40,dutton75",
        "Heylen 1971": "heylen71,dutton75"}
KIND_TEXT = {"W": "$W$", "alphaN": "$\\alpha/N$", "DTmu": "$D_T/\\mu$"}
ROLE_TEXT = {"train": "train", "held_out": "held out", "check": "check",
             "validation": "validation"}
# held out, but they chose which factors train (the variant table's
# held-out columns, checked there)
VALIDATION = ("dt_koizumi86", "w_pack62_195", "w_pack92_195")


def swarm_numbers() -> tuple:
    n, lines = {}, []
    rows = dataset_scores()
    by = {r["d"].name: r for r in rows}
    # the dataset table, by kind and role
    order = {"W": 0, "DTmu": 1, "alphaN": 2}
    rank = {"train": 0, "validation": 1, "held_out": 1, "check": 2}
    tab = []
    for r in sorted(rows, key=lambda r: (order[r["d"].kind],
                                         rank[r["role"]], r["d"].x.min())):
        d = r["d"]
        tab.append(
            f"{LABEL[d.name]}~\\cite{{{CITE[d.group]}}} & {KIND_TEXT[d.kind]} & "
            f"{_rng(d.x.min())}--{_rng(d.x.max())} & {len(d.x)} & "
            f"{ROLE_TEXT[r['role']]} & {d.sigma * 100:.0f} & "
            f"{_chi(r['bi']['chi'])} & {_chi(r['tr']['chi'])} & "
            f"{r['bi']['rms']:.0f} & {r['tr']['rms']:.0f}\\\\")
    lines.append("\\def\\wgSwTable{" + "\n".join(tab) + "}")
    n["sw.nsets"] = str(len(rows))
    n["sw.npts"] = str(sum(len(r["d"].x) for r in rows))
    for role in ("train", "held_out", "check", "validation"):
        n[f"sw.n.{role}"] = str(sum(r["role"] == role for r in rows))
    # the pooled scores the score table and the bottom line quote
    groups = {
        "wtr": [r for r in rows if r["d"].kind == "W" and r["role"] == "train"],
        "atr": [r for r in rows if r["d"].kind == "alphaN" and r["role"] == "train"],
        "alow": [r for r in rows if r["d"].name in KD.HELD_BY_RULE],
        "dtk": [by["dt_koizumi86"]],
    }
    for g, rs in groups.items():
        for which in ("bi", "tr"):
            n[f"sw.{g}.{which}.rms"] = f"{_pooled(rs, which):.1f}"
            n[f"sw.{g}.{which}.mean"] = _pm(_pooled(rs, which, 'mean'), 0)
        n[f"sw.{g}.n"] = str(sum(len(r["d"].x) for r in rs))
        n[f"sw.{g}.sets"] = str(len(rs))
    for name in ("w_hunter88", "w_nakamura90", "a_jacques86", "a_kruithof40",
                 "a_heylen71", "a_bhatta79", "a_dallarmi92", "dt_koizumi86",
                 "a_specht80_direct", "a_kruithof40_low", "a_heylen71_low",
                 "a_dallarmi92_low", "w_pack92", "w_dallarmi92", "w_bowe60",
                 "w_english53", "w_kucuk81", "w_pack62_300", "w_pack62_195",
                 "w_pack92_195", "dt_alamin87"):
        r = by[name]
        for which in ("bi", "tr"):
            n[f"sw.{name}.{which}.rms"] = f"{r[which]['rms']:.1f}"
            n[f"sw.{name}.{which}.mean"] = _pm(r[which]["mean"], 0)
            n[f"sw.{name}.{which}.chi"] = _chi(r[which]["chi"])
    # the 1960 drift velocities against the 1988 ones where both measured
    bo, hu = by["w_bowe60"]["d"], by["w_hunter88"]["d"]
    o = np.argsort(hu.x)
    rat = [y / math.exp(float(np.interp(math.log(x), np.log(hu.x[o]),
                                        np.log(hu.y[o]))))
           for x, y in zip(bo.x, bo.y) if hu.x.min() <= x <= hu.x.max()]
    n["bowe.hunter"] = _pm(100 * (float(np.mean(rat)) - 1), 0)
    # Pack et al. 1992 at 300 K: neighbouring points (within 5 % in E/N)
    # against the model
    pk = by["w_pack92"]
    o = np.argsort(pk["d"].x)
    lx, lr = np.log(pk["d"].x[o]), np.log(pk["tr"]["ratio"][o])
    jumps = [abs(lr[i + 1] - lr[i]) for i in range(len(lx) - 1)
             if lx[i + 1] - lx[i] < 0.05]
    n["pack92.jump"] = f"{100 * (math.exp(max(jumps)) - 1):.0f}"
    # the excitation against the ionization below 50 Td (trained set)
    for en in (30.0,):
        st = swarm.KB.State(x={"Kr": 1.0}, tg=300.0, n_col=swarm.N_COL,
                            ne=swarm.NE_FRAC * swarm.N_COL)
        procs = swarm.processes(0.0)
        d = swarm_point(en, 300.0)
        s = swarm.KB.Solver(procs, n_cells=CELLS_SWARM, emax=d["emax"],
                            elastic=("Kr",), tol=1e-10, maxit=4000)
        sol = s.solve(st, en)
        kx = sum(sol.grid.rate(q.p, sol.F) for q in procs
                 if q.kind == "excitation")
        ki = sum(sol.grid.rate(q.p, sol.F) for q in procs
                 if q.kind == "ionization")
        n["low.exc_ion"] = f"{kx / ki:.0f}"
    # how far the low-field datasets part from one another: each pair at the
    # E/N where both measured (one read log-log between its own points)
    low = [by[nm]["d"] for nm in KD.HELD_BY_RULE if len(by[nm]["d"].x) > 1]
    part = 1.0
    for a in low:
        for b in low:
            if a is b:
                continue
            o = np.argsort(b.x)
            bx, by_ = np.log(b.x[o]), np.log(b.y[o])
            for x, y in zip(a.x, a.y):
                if b.x.min() <= x <= b.x.max():
                    yb = math.exp(float(np.interp(math.log(x), bx, by_)))
                    part = max(part, y / yb, yb / y)
    n["sw.low.part"] = f"{part:.1f}"
    n["sw.low.sets"] = str(len(KD.HELD_BY_RULE))
    # the grid: the train sets' pooled scores on half and twice the cells
    dmax = 0.0
    for g in ("wtr", "atr"):
        base = _pooled(groups[g], "tr")
        for cells in (CELLS_SWARM // 2, CELLS_SWARM * 2):
            r = np.concatenate([np.array([swarm_point(float(x), d.tg, "tr",
                                                      cells=cells)
                                          [d.kind + "_m"] for x in d.x])
                                / d.y for d in (q["d"] for q in groups[g])])
            dmax = max(dmax, abs(100 * _rms(r - 1.0) - base))
    n["sw.cells"] = str(CELLS_SWARM)
    n["sw.grid"] = f"{dmax:.2f}"
    # held-out drift velocities: the range of their rms errors
    wh = [r for r in rows if r["d"].kind == "W"
          and r["role"] in ("held_out", "validation")]
    n["sw.wh.sets"] = str(len(wh))
    n["sw.wh.tr.lo"] = f"{min(r['tr']['rms'] for r in wh):.0f}"
    n["sw.wh.tr.hi"] = f"{max(r['tr']['rms'] for r in wh):.0f}"
    worse = sorted(r["d"].name for r in wh if r["tr"]["chi"] > r["bi"]["chi"])
    if worse != ["w_dallarmi92", "w_kucuk81", "w_pack92"]:
        raise AssertionError(f"the text names the held-out drift velocities "
                             f"the training moves away: {worse}")
    n["sw.wh.better"] = str(len(wh) - len(worse))
    # the claims the text makes of these scores
    for g in ("wtr", "atr", "alow"):
        if not float(n[f"sw.{g}.tr.rms"]) < float(n[f"sw.{g}.bi.rms"]):
            raise AssertionError(f"the text says training improves {g}")
    if not all(r["tr"]["chi"] <= 3.0 for r in groups["wtr"] + groups["atr"]):
        raise AssertionError("the text says every train set is within the "
                             "consistency threshold")
    # the training's record
    names = TRAINED["names"]
    z = np.array(TRAINED["z"])
    post = TRAINED["posterior"]
    nn = TRAINED["nn"]
    ftab = []
    what = {"m_0.03": "momentum transfer at 0.03\\unit{eV}",
            "m_0.2": "the same at 0.2\\unit{eV}",
            "m_0.6": "at 0.6\\unit{eV} (the Ramsauer minimum)",
            "m_2": "at 2\\unit{eV}", "m_6": "at 6\\unit{eV}",
            "m_20": "at 20\\unit{eV} and above",
            "x_1s": "the four 5s levels (1s)",
            "x_2p": "the ten 5p levels (2p)",
            "x_hi": "the 37 higher levels"}
    free = TRAINED["free"]
    for i, nm in enumerate(names):
        if not free[i]:
            # held at Biagi's value: the variants below chose which train
            if z[i] != 0.0 or post[nm] != [1.0, 1.0, 1.0]:
                raise AssertionError(f"{nm} is not trained but moved")
            ftab.append(f"{what[nm]} & \\multicolumn{{5}}{{l}}{{1, not "
                        f"trained (Table~\\ref{{tab:var}})}}\\\\")
            continue
        lo, mid, hi = post[nm]
        ftab.append(f"{what[nm]} & {math.exp(z[i]):.3f} & {mid:.3f} & "
                    f"[{lo:.3f}, {hi:.3f}] & {nn['factors'][i]:.3f} "
                    f"$\\pm$ {nn['z_sd'][i]:.2f} & {nn['r2'][i]:.2f}\\\\")
        n[f"fac.{nm}"] = f"{math.exp(z[i]):.2f}"
        n[f"fac.{nm}.lo"] = f"{lo:.2f}"
        n[f"fac.{nm}.hi"] = f"{hi:.2f}"
    lines.append("\\def\\wgFacTable{" + "\n".join(ftab) + "}")
    n["fac.n"] = str(sum(free))
    knots = TRAINED["m_knots_eV"]
    held = [f"{e:g}" for e, f in zip(knots, free) if not f]
    moved = [e for e, f in zip(knots, free) if f]
    if held != ["0.03", "0.2"] or moved != knots[2:] or not all(free[6:]):
        raise AssertionError("the text says the momentum transfer is held "
                             "at 0.03 and 0.2 eV and the rest train")
    n["fac.held"] = " and ".join(held)
    n["fac.free.lo"] = f"{moved[0]:g}"
    n["fac.free.hi"] = f"{moved[-1]:g}"
    n["fac.free.n"] = WORDS[len(moved)]
    n["fac.box"] = f"{math.exp(-TRAINED['emulator']['box']):.2f}"
    n["fac.box.hi"] = f"{math.exp(TRAINED['emulator']['box']):.2f}"
    emu = TRAINED["emulator"]
    n["emu.n"] = str(emu["n_design"])
    n["emu.nodes"] = str(emu["nodes"])
    n["emu.W"] = f"{100 * emu['rms_ln'][0]:.2f}"
    n["emu.DT"] = f"{100 * emu['rms_ln'][1]:.2f}"
    n["emu.a"] = f"{100 * emu['rms_ln'][2]:.2f}"
    n["mcmc.acc"] = f"{TRAINED['mcmc']['acceptance']:.2f}"
    # the network reads the trained factors only (None at the others)
    cov = [c for c, f in zip(nn["coverage"], free) if f]
    r2 = sorted((r, nm) for r, nm, f in zip(nn["r2"], names, free) if f)
    n["nn.cov.lo"] = f"{min(cov):.2f}"
    n["nn.cov.hi"] = f"{max(cov):.2f}"
    n["nn.r2.lo"] = f"{r2[0][0]:.2f}"
    n["nn.r2.next"] = f"{r2[1][0]:.2f}"
    if r2[0][1] != "x_2p":
        raise AssertionError("the text says the data pin the 5p excitation "
                             "least")
    con = TRAINED["consistency"]
    n["cons.dropped"] = ", ".join(
        f"{LABEL[nm]} ({v:.1f})" for nm, v in con["dropped"])
    n["cons.ndropped"] = str(len(con["dropped"]))
    n["cons.cands"] = str(len(con["kept"]) + len(con["dropped"]))
    scan = {s["threshold"]: s for s in TRAINED["threshold_scan"]}
    if not scan[2.0]["train"] == scan[3.0]["train"] == TRAINED["train"]:
        raise AssertionError("the text says thresholds 2 and 3 give one set")
    extra = [x for x in scan[4.0]["train"] if x not in scan[3.0]["train"]]
    if sorted(extra) != ["w_bowe60", "w_pack62_195"]:
        raise AssertionError("the text says the 1960 drift velocities and "
                             "one more set join at 4")
    n["cons.at4"] = " and ".join(LABEL[x] for x in extra)
    n["cons.at4.n"] = WORDS[len(extra)]
    i6 = names.index("m_6")
    n["cons.at4.m6"] = f"{scan[4.0]['factors'][i6]:.2f}"
    i1 = names.index("x_1s")
    n["cons.at4.x1s"] = f"{scan[4.0]['factors'][i1]:.2f}"
    # which factors train: each variant's least squares on the train set,
    # every dataset scored on the true solver (the factors: the record's
    # free_variants)
    var = TRAINED["free_variants"]
    rows_v = (("Biagi's file", "bi", None,
               "none (Biagi's file as published)"),
              ("excitation only", "exc", (False,) * 6 + (True,) * 3,
               "the three excitations"),
              ("knots below 1 eV fixed", "b1", (False,) * 3 + (True,) * 6,
               "these, and the momentum transfer at 2 to 20\\unit{eV}"),
              ("knots below 0.6 eV fixed", "v3", (False,) * 2 + (True,) * 7,
               "these, and at 0.6\\unit{eV} too (v0\\_3)"),
              ("all nine", "all", (True,) * 9, "all nine"))
    cols = (("w_hunter88", "hunter"), ("dt_koizumi86", "dtk"),
            ("w_pack62_195", "p62"), ("w_pack92_195", "p92"))
    if tuple(c for c, _ in cols[1:]) != VALIDATION:
        raise AssertionError("the validation sets are the variant table's "
                             "held-out columns")
    if set(var) != {r[0] for r in rows_v}:
        raise AssertionError(f"the record's variants: {sorted(var)}")
    npts = {nm: len(by[nm]["d"].x) for nm in TRAINED["train"]}
    n["var.npts"] = str(sum(npts.values()))
    if tuple(var["knots below 0.6 eV fixed"]["free"]) != tuple(free):
        raise AssertionError("v0_3 trains another variant than the text says")
    # the variants' factors from the record, scored here on the solver as
    # every other dataset is (so a change to swarm_mc's ratios reaches
    # them); the record's own chi2 must agree where those ratios did not
    # move (W and alpha/N)
    vt, pooled, ch = [], {}, {}
    for key, tag, mask, text in rows_v:
        v = var[key]
        if (None if v["free"] is None else tuple(v["free"])) != mask:
            raise AssertionError(f"variant {key}: {v['free']}")
        which = {"bi": "bi", "v3": "tr"}.get(tag, "var_" + tag)
        if which.startswith("var_"):
            fv = v["factors"]
            SETS[which] = xsec.Factors(m_knots=tuple(knots), m=tuple(fv[:6]),
                                       exc=tuple(fv[6:9]), ion=1.0)
        elif max(abs(a - b) for a, b in zip(v["factors"], (
                [1.0] * 9 if tag == "bi" else TRAINED["m"] + [
                    TRAINED["exc"][g] for g in xsec.EXC_GROUPS]))) > 1e-12:
            raise AssertionError(f"variant {key} is not the {which} set")
        chi = {}
        for k in list(npts) + [c for c, _ in cols]:
            r = by[k]
            m = _model(r["d"], which)
            chi[k] = float(np.mean((np.log(m / r["d"].y) / r["sig"]) ** 2))
            if (r["d"].kind != "DTmu"
                    and abs(chi[k] - v["chi2"][k]) > 0.01 * max(1.0, chi[k])):
                raise AssertionError(f"{key} {k}: the record's chi2 "
                                     f"{v['chi2'][k]} against {chi[k]:.3f}")
        ch[tag] = chi
        pooled[tag] = sum(npts[k] * chi[k] for k in npts) / sum(npts.values())
        vt.append(f"{text} & {pooled[tag]:.2f} & "
                  f"{max(chi[k] for k in npts):.2f} & "
                  + " & ".join(f"{chi[k]:.2f}" for k, _ in cols) + "\\\\")
        n[f"var.{tag}.train"] = f"{pooled[tag]:.2f}"
        for k, short in cols:
            n[f"var.{tag}.{short}"] = f"{chi[k]:.1f}"
    lines.append("\\def\\wgVarTable{" + "\n".join(vt) + "}")
    for k, _ in cols[1:]:
        if not ch["all"][k] > ch["bi"][k] > ch["v3"][k]:
            raise AssertionError(f"the text says all nine worsen {k} and "
                                 f"v0_3's improve it")
        if any(ch[t][k] < ch["v3"][k] for t in ch):
            raise AssertionError(f"the text says v0_3's variant fits {k} best")
    if not (pooled["all"] < pooled["v3"] < pooled["b1"] < pooled["exc"]
            < pooled["bi"]):
        raise AssertionError("the text's order of the train sets' fits")
    fa = var["all nine"]["factors"]
    n["var.all.m003"] = f"{fa[0]:.2f}"
    n["var.all.m06"] = f"{fa[2]:.2f}"
    # the Monte Carlo: what the two-term leaves out, Biagi's file and the
    # trained set
    mc = []
    chk = {r[0]: r for r in swarm_mc.TRAINED_CHECK}
    for en, rw, ra, rd in swarm_mc.RATIOS:
        c = chk.get(en)
        mc.append(f"{en:g} & {rw:.3f} & {ra:.3f} & {rd:.2f} & "
                  + (f"{c[1]:.3f} & {c[2]:.3f}" if c else "-- & --") + "\\\\")
    lines.append("\\def\\wgMcTable{" + "\n".join(mc) + "}")
    worst = max(max(abs(c[1] / r[1] - 1), abs(c[2] / r[2] - 1))
                for r in swarm_mc.RATIOS for c in swarm_mc.TRAINED_CHECK
                if c[0] == r[0])
    n["mc.chk.max"] = f"{100 * worst:.1f}"
    dts = [r[3] for r in swarm_mc.RATIOS] + [r[1] for r in swarm_mc.DT_LOW
                                            if r[0] >= 10.0]
    n["mc.dt.lo"] = f"{100 * (1 / max(dts) - 1):.0f}"
    n["mc.dt.hi"] = f"{100 * (1 / min(dts) - 1):.0f}"
    # how close the two-term solution is to the Monte Carlo's flux values
    n["mc.flux.w"] = f"{100 * max(abs(r[1] - 1) for r in swarm_mc.FLUX if r[0] >= 10):.1f}"
    low = [r for r in swarm_mc.FLUX if r[0] < 10]
    n["mc.flux.wlow"] = f"{100 * max(abs(r[1] - 1) for r in low):.1f}"
    n["mc.flux.wlow.en"] = f"{min(r[0] for r in low):g}"
    fa = [r[2] - 1 for r in swarm_mc.FLUX if r[2] is not None and r[0] >= 40]
    n["mc.flux.a.lo"] = f"{100 * min(fa):.1f}"
    n["mc.flux.a.hi"] = f"{100 * max(fa):.1f}"
    n["mc.w300"] = f"{100 * (swarm_mc.RATIOS[-1][1] - 1):.0f}"
    n["mc.a300"] = f"{100 * (1 - swarm_mc.RATIOS[-1][2]):.0f}"
    # v0_2's Monte Carlo check of the two-term alpha/N, 20 to 50 Td
    dev, mdev = [], []
    for en, (a_mc, m_mc, counted) in sorted(LOKI_MC.items()):
        d = swarm_point(en, 300.0, "bi")
        dev.append((en, d["alphaN"] * 1e-4 / a_mc - 1.0, counted))
        mdev.append(abs(d["mean"] / m_mc - 1.0))
    n["mc.max"] = f"{100 * max(abs(x[1]) for x in dev):.1f}"
    n["mc.mean.max"] = f"{100 * max(mdev):.2f}"
    n["mc.list"] = ", ".join(_pm(100 * x[1], 1) for x in dev)
    # the Hornbeck-Molnar band on the trained set
    for en, key in ((14.12, "14"), (30.0, "30"), (100.0, "100")):
        r = swarm_point(en, 300.0, hm=1.0)["alphaN"] / swarm_point(en, 300.0)["alphaN"]
        n[f"hm.x{key}"] = f"{r:.1f}"
    n["hm.levels"] = str(sum(1 for p in levels() if p.threshold >= E_ASSOC_EV))
    n["hm.elo"] = f"{min(p.threshold for p in levels() if p.threshold >= E_ASSOC_EV):.2f}"
    n["hm.ehi"] = f"{max(p.threshold for p in levels()):.1f}"
    n["hm.eai"] = f"{E_ASSOC_EV:.2f}"
    v, lo, hi = HORNBECK_KR
    n["hm.ap"] = f"${v:g}^{{+{hi:g}}}_{{-{lo:g}}}$"
    atr_hm = [r for r in groups["atr"]]
    rat = np.concatenate([np.array([swarm_point(float(x), 300.0, hm=1.0)["alphaN_m"]
                                    for x in r["d"].x]) / r["d"].y
                          for r in atr_hm]) - 1.0
    n["swhm.atr.mean"] = _pm(100 * float(np.mean(rat)), 0)
    # the figures
    lines += _swarm_figures(rows)
    return n, lines


def _rng(x):
    return f"{x:.3g}" if x < 1 else f"{x:.0f}" if x >= 10 else f"{x:.1f}"


def _chi(c):
    return f"{c:.1f}" if c < 99.5 else f"{c:.0f}"


def _swarm_figures(rows) -> list:
    lines = []
    ens = np.geomspace(1e-3, 60.0, 70)
    for which, key in (("tr", "Wtr"), ("bi", "Wbi")):
        p_ = " ".join(f"({en:.5g},{swarm_point(en, 300.0, which)['W_m']:.5g})"
                      for en in ens)
        lines.append(f"\\def\\wgSw{key}{{{p_}}}")
    ens = np.geomspace(12.0, 300.0, 40)
    for which, key, hm in (("tr", "Atr", 0.0), ("bi", "Abi", 0.0),
                           ("tr", "AHM", 1.0)):
        p_ = " ".join(f"({en:.5g},{swarm_point(en, 300.0, which, hm)['alphaN_m']:.5g})"
                      for en in ens)
        lines.append(f"\\def\\wgSw{key}{{{p_}}}")
    ens = np.geomspace(3e-3, 1.0, 40)
    for which, key in (("tr", "Dtr"), ("bi", "Dbi")):
        p_ = " ".join(f"({en:.5g},{swarm_point(en, 300.0, which)['DTmu_m']:.5g})"
                      for en in ens)
        lines.append(f"\\def\\wgSw{key}{{{p_}}}")
    for kind, role, key in (("W", "train", "WT"), ("W", "held_out", "WH"),
                            ("alphaN", "train", "AT"),
                            ("alphaN", "held_out", "AH"),
                            ("alphaN", "check", "AC"),
                            ("DTmu", "held_out", "DH")):
        pts = []
        for r in rows:
            d = r["d"]
            # the figure draws a validation set as held out (open marks)
            shown = "held_out" if r["role"] == "validation" else r["role"]
            if d.kind != kind or shown != role or d.tg != 300.0 \
                    and not (d.tg == 301.0):
                continue
            if kind == "alphaN" and role == "check" and \
                    d.name != "a_specht80_total":
                continue
            pts += [f"({x:.5g},{y:.5g})" for x, y in zip(d.x, d.y)
                    if (kind != "W" or x <= 60.0) and
                    (kind != "alphaN" or 12.0 <= x <= 300.0)]
        lines.append(f"\\def\\wgSw{key}data{{{' '.join(pts)}}}")
    return lines


# ------------------------------------------------------------- e-beam
def beam_numbers() -> dict:
    n = {}
    Y = K.beam_yields()
    # v0_3: the same degradation on v0_2's cross sections (Biagi's file as
    # published): the partition the trained factors move
    Y2 = K2.beam_yields()
    for key in ("W", "nex", "ws"):
        meas = {"W": SAITO_WS * (1 + SAITO_NEX), "nex": SAITO_NEX,
                "ws": SAITO_WS}[key]
        fmt = {"W": ".1f", "nex": ".2f", "ws": ".2f"}[key]
        n[f"beam.{key}.v02"] = f"{Y2[key]:{fmt}}"
        n[f"beam.{key}.dev.v02"] = _pm(100 * (Y2[key] / meas - 1),
                                       1 if key == "ws" else 0)
    if not (abs(Y["W"] / (SAITO_WS * (1 + SAITO_NEX)) - 1)
            < abs(Y2["W"] / (SAITO_WS * (1 + SAITO_NEX)) - 1)
            and abs(Y["nex"] / SAITO_NEX - 1) < abs(Y2["nex"] / SAITO_NEX - 1)):
        raise AssertionError("the text says the trained set moves the "
                             "beam's W and N_ex/N_i toward Saito's")
    hv = next(r.hv_eV for r in build() if r.rid == "KR2_RAD")
    ceil = (1.0 + Y["nex"]) * hv / Y["W"]
    n.update({
        "beam.W": f"{Y['W']:.1f}", "beam.nex": f"{Y['nex']:.2f}",
        "beam.pp": f"{Y['pp']:.3f}", "beam.s": f"{Y['s']:.3f}",
        "beam.ws": f"{Y['ws']:.2f}", "beam.e0": f"{Y['e0']:.0f}",
        "beam.ceiling.eV": f"{Y['ceiling']:.0f}",
        "beam.eps": f"{Y['eps_sec']:.2f}",
        "beam.cascade": f"{Y['cascade']:.2f}",
        "beam.ceil": f"{ceil:.3f}",
        "beam.php": f"{1.0 + Y['nex']:.3f}",
        "hv147": f"{hv:.2f}",
        "saito.ws": f"{SAITO_WS:.1f}", "saito.nex": f"{SAITO_NEX:.2f}",
        "saito.ideal": f"{hv / SAITO_WS:.3f}",
        "saito.W": f"{SAITO_WS * (1 + SAITO_NEX):.1f}",
        "moz.short": f"{100 * (1 - MOROZOV[0] / (hv / SAITO_WS)):.0f}",
        "beam.ws.dev": _pm(100 * (Y["ws"] / SAITO_WS - 1), 1),
        "beam.nex.dev": _pm(100 * (Y["nex"] / SAITO_NEX - 1)),
        "beam.W.dev": _pm(100 * (Y["W"] / (SAITO_WS * (1 + SAITO_NEX)) - 1)),
        "moz.meas": f"{MOROZOV[0]:.2f}", "moz.err": f"{MOROZOV[1]:.2f}",
        "eck.meas": f"{ECKSTROM:.2f}",
        "moz.atm": f"{K.MOROZOV_ATM:.3f}", "moz.hpa": "1000",
        "moz.w.lo": f"{K.MOROZOV_W[0]:.0f}", "moz.w.hi": f"{K.MOROZOV_W[1]:.0f}",
        "eck.w": f"{K.ECKSTROM_W / 1e3:.0f}",
        "beam.pulse": f"{K.BEAM_PULSE * 1e6:g}",
        "beam.tail.lo": f"{K.BEAM_TAIL_LO * 1e6:g}",
        "beam.tail.hi": f"{K.BEAM_TAIL_HI * 1e6:g}",
    })
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        bp = K.beam_params(eedf)
        keys = list(bp)
        got = dict(zip(keys, K.run_many([("beam", bp[k]) for k in keys])))
        for k, s in got.items():
            pairs = s["e_dep_eV"] / s["beam"]["W"]
            ev = s["events"]
            vuv = s["emission_vuv"] / 100.0
            name = f"{k[0]}.{k[1]:g}".replace("eckstrom", "eck").replace(
                "morozov", "moz")
            n[f"{name}.{tag}"] = f"{vuv:.3f}"
            n[f"{name}.{tag}.G2"] = f"{ev['G2'] / pairs:.3f}"
            exc = sum(v for r, v in ev.items()
                      if r.startswith(("G1_", "C1")))
            n[f"{name}.{tag}.G1"] = f"{exc / pairs:.3f}"
            n[f"{name}.{tag}.te"] = f"{s['te_mean_eV']:.2f}"
            n[f"{name}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
            n[f"{name}.{tag}.ph"] = f"{ev['KR2_RAD'] / pairs:.3f}"
            n[f"{name}.{tag}.rec"] = f"{ev['R11p'] / pairs:.3f}"
            n[f"{name}.{tag}.krs"] = _sci(s["n_max"]["Kr*"])
            n[f"{name}.{tag}.ne"] = _sci(s["ne_peak"])
            n[f"{name}.{tag}.closure"] = _sci(abs(s["closure"]), 0)
            # the loss below v0_3's ceiling (one photon per excitation and
            # per ion), the quantity a report 'near the limit' tests
            n[f"{name}.{tag}.loss"] = f"{100 * (1 - vuv / ceil):.0f}"
        m = got[("morozov", K.MOROZOV_W[1])]
        v = m["emission_vuv"] / 100.0
        n[f"moz.{tag}.dev"] = _pm(100 * (v / MOROZOV[0] - 1))
        n[f"moz.{tag}.sig"] = f"{(v - MOROZOV[0]) / MOROZOV[1]:.1f}"
        e1 = got[("eckstrom", 1.0)]["emission_vuv"] / 100.0
        n[f"eck.{tag}.dev"] = _pm(100 * (e1 / ECKSTROM - 1))
        # what the text says of these runs
        if not v > MOROZOV[0]:
            raise AssertionError("the text says v0_3 lies above Morozov")
        if not 1 - got[("eckstrom", 2.0)]["emission_vuv"] / 100 / ceil < 0.10:
            raise AssertionError("the text says v0_3 loses little at 2 atm")
        if not (got[("eckstrom", 0.5)]["emission_vuv"] < 100 * e1
                < got[("eckstrom", 2.0)]["emission_vuv"]
                and e1 < ECKSTROM
                and got[("eckstrom_pen", 1.0)]["emission_vuv"] > 100 * e1
                and got[("eckstrom_form", 1.0)]["emission_vuv"] > 100 * e1):
            raise AssertionError("the text's account of the high-pump runs "
                                 "no longer holds")
        n[f"moz.{tag}.flat"] = _pm(100 * (
            got[("morozov", K.MOROZOV_W[0])]["emission_vuv"]
            / got[("morozov", K.MOROZOV_W[1])]["emission_vuv"] - 1), 1)
    # v0_2 at the same two points, from its own cache (kr_docs_v0_2)
    bp2 = K2.beam_params("auto")
    m2, e2 = K2.run_many([("beam", bp2[("morozov", K2.MOROZOV_W[1])]),
                          ("beam", bp2[("eckstrom", 1.0)])])
    n["moz.kin.dev.v02"] = _pm(100 * (m2["emission_vuv"] / 100 / MOROZOV[0]
                                      - 1))
    n["eck.kin.dev.v02"] = _pm(100 * (e2["emission_vuv"] / 100 / ECKSTROM
                                      - 1))
    return n


# ------------------------------------------------------ the uncertainty
UNC_LABELS = {
    "k_5s5p": r"\el{} + \Krs{} \ra \Krpp{}, scale",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{}, Penning",
    "kr2_form": r"\Krs{} + 2Kr \ra \Krts{} + Kr",
    "kr2_pen": r"\Krts{} + \Krts{}, Penning",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}",
    "hm_kr": r"Kr(nl) + Kr \ra \Krtp{} + \el{}, the share",
    "R11p": r"\Krtp{} + \el{} \ra \Krpp{} + Kr",
    "K15": r"\Krp{} + 2Kr \ra \Krtp{} + Kr",
    "KR2_RAD": r"\Krts{} lifetime",
}
UNC_SHORT = {
    "k_5s5p": r"the 5s \ra 5p cross section", "kps_kr": r"\Krpp{} + Kr",
    "kr_pen": r"\Krs{} + \Krs{}", "kr2_form": r"\Krs{} + 2Kr",
    "kr2_pen": r"\Krts{} + \Krts{}", "kstep_s": r"\el{} + \Krs{} \ra \Krp",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp", "R11p": r"\Krtp{} + \el",
    "hm_kr": r"the Hornbeck--Molnar share",
    "K15": r"\Krp{} + 2Kr", "KR2_RAD": r"the \Krts{} lifetime",
    "xs": r"the trained cross sections",
}
UNC_LABELS["xs"] = r"the electron cross sections (trained)"


def _ends_text(name, lo, hi):
    if name in BANDS:
        b = BANDS[name]
        if b.unit in ("factor", "fraction"):
            return f"{lo:g} to {hi:g}"
        unit = {"cm3/s": r"\unit{cm^3\,s^{-1}}",
                "cm6/s": r"\unit{cm^6\,s^{-1}}"}[b.unit]
        lo_t = "0" if lo == 0 else _sci(lo)
        hi_t = "0" if hi == 0 else _sci(hi)
        return f"${lo_t}$ to ${hi_t}${unit}"
    if name == "R11p":
        return f"$\\pm{100 * K.R11P_ERR:.0f}\\,\\%$"
    if name == "K15":
        return f"${_sci(K.K15_LO)}$ to ${_sci(K.K15_HI)}\\unit{{cm^6\\,s^{{-1}}}}$"
    if name == "KR2_RAD":
        return (f"{K.TAU_RAD - K.TAU_RAD_ERR:.0f} to "
                f"{K.TAU_RAD + K.TAU_RAD_ERR:.0f}\\unit{{ns}}")
    raise KeyError(name)


def uncertainty_numbers() -> tuple:
    """Each input alone at its two ends, on both closures: the 147-nm
    share, the near IR and the final gas temperature; half the difference
    of its two ends is its term, summed in quadrature, and the excursions
    from the reference summed by sign."""
    n, lines = {}, []
    out = {}
    for eedf in K.CLOSURES:
        ref = K.run_many([("field", K.ref_params(eedf))])[0]
        items = []
        for name, end, v, p in K.band_params(eedf):
            items.append((name, end, v, p))
        for rid, end, f, p in K.const_params(eedf):
            items.append((rid, end, f, p))
        got = K.run_many([("field", p) for *_, p in items])
        res = {}
        for (name, end, v, _), s in zip(items, got):
            res[(name, end)] = s
        out[eedf] = (ref, res)
    names = ([b.name for b in K.row_bands()]
             + [rid for rid, _, _ in K.CONSTANTS])
    # v0_3: the trained cross sections, as the reference point over draws
    # of the training's posterior (kr_docs_v0_3.xs_jobs): their standard
    # deviation is the term, one standard deviation like the others
    xs = {eedf: K.run_many(K.xs_jobs(eedf)) for eedf in K.CLOSURES}

    def xs_stat(eedf, key):
        v = np.array([s[key] for s in xs[eedf]])
        return float(v.mean()), float(v.std(ddof=1))
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        ref, res = out[eedf]
        y0 = ref["emission_vuv"]
        terms, up, down = {}, 0.0, 0.0
        for name in names:
            vals = []
            for end in ("low", "high"):
                s = res.get((name, end), ref)
                vals.append(s["emission_vuv"])
                dv = s["emission_vuv"] - y0
                up += max(dv, 0.0)
                down += min(dv, 0.0)
            terms[name] = 0.5 * abs(vals[1] - vals[0])
        m, sd = xs_stat(eedf, "emission_vuv")
        terms["xs"] = sd
        up += sd
        down -= sd
        n[f"unc.{tag}.xs.sd"] = f"{sd:.2f}"
        n[f"unc.{tag}.xs.mean"] = f"{m:.2f}"
        n[f"unc.{tag}.xs.n"] = str(len(xs[eedf]))
        n[f"unc.{tag}.xs.pct"] = f"{100 * sd / y0:.1f}"
        quad = math.sqrt(sum(t * t for t in terms.values()))
        top = sorted(terms.items(), key=lambda kv: -kv[1])
        n[f"unc.{tag}.y0"] = f"{y0:.2f}"
        n[f"unc.{tag}.quad"] = f"{quad:.1f}"
        n[f"unc.{tag}.quad.pct"] = f"{100 * quad / y0:.0f}"
        n[f"unc.{tag}.up"] = f"+{up:.1f}"
        n[f"unc.{tag}.down"] = f"$-${-down:.1f}"
        n[f"unc.{tag}.top1"] = UNC_SHORT[top[0][0]]
        n[f"unc.{tag}.top1.t"] = f"{top[0][1]:.2f}"
        n[f"unc.{tag}.top2"] = UNC_SHORT[top[1][0]]
        n[f"unc.{tag}.top2.t"] = f"{top[1][1]:.2f}"
        n[f"unc.{tag}.top3"] = UNC_SHORT[top[2][0]]
        n[f"unc.{tag}.top3.t"] = f"{top[2][1]:.2f}"
        n[f"unc.{tag}.top1.share"] = f"{100 * top[0][1] ** 2 / quad ** 2:.0f}"
    ref_k, res_k = out["auto"]
    top_k = max(names, key=lambda x: 0.5 * abs(
        res_k.get((x, "high"), ref_k)["emission_vuv"]
        - res_k.get((x, "low"), ref_k)["emission_vuv"]))
    if top_k != "kr2_form":
        raise AssertionError(f"the text says Kr* + 2Kr leads; {top_k} does")
    for name in ("kr2_pen", "kstep_s", "kstep_p"):
        cut = 1 - res_k[(name, "high")]["emission_vuv"] / ref_k["emission_vuv"]
        if not 0.15 < cut < 0.55:
            raise AssertionError(f"the text says {name} at its top cuts the "
                                 f"share by a fifth to a half: {cut:.2f}")
    # v0_2 and v0_3: the Hornbeck-Molnar share at its top barely moves the
    # reference point (direct and associative ionization make well under
    # a percent of its ions next to Penning's)
    hm_cut = res_k[("hm_kr", "high")]["emission_vuv"] / ref_k["emission_vuv"] - 1
    if not abs(hm_cut) < 0.01:
        raise AssertionError(f"the text says hm_kr at its top moves the "
                             f"share by under 1 %: {hm_cut:+.4f}")
    n["unc.hm.cut"] = _pm(100 * hm_cut, 2)
    b = BANDS["kr2_form"]
    n["form.ratio"] = f"{b.hi / b.lo:.1f}"
    n["unc.n"] = str(len(names))
    n["unc.nband"] = str(len(K.row_bands()))
    n["unc.nrow"] = str(len(K.CONSTANTS))
    # the table: each input, its ends, and the outputs at its two ends,
    # the largest term (kinetic) first
    ref_k, res_k = out["auto"]
    ref_t, res_t = out["two_t"]

    def term(x):
        return abs(res_k.get((x, "high"), ref_k)["emission_vuv"]
                   - res_k.get((x, "low"), ref_k)["emission_vuv"])
    xs_k = xs_stat("auto", "emission_vuv")[1]
    for name in sorted(names + ["xs"],
                       key=lambda x: -(2 * xs_k if x == "xs" else term(x))):
        if name == "xs":
            cells = []
            for eedf, key, fmt in (("auto", "emission_vuv", ".2f"),
                                   ("two_t", "emission_vuv", ".2f"),
                                   ("auto", "emission_nir", ".2f"),
                                   ("auto", "tg_end_K", ".0f")):
                m, sd = xs_stat(eedf, key)
                cells.append(f"{m:{fmt}} $\\pm$ {sd:{fmt}}")
            lines.append(
                f"{UNC_LABELS['xs']} & {len(xs['auto'])} and "
                f"{len(xs['two_t'])} draws of the posterior & "
                + " & ".join(cells) + "\\\\")
            continue
        if name in BANDS:
            b = BANDS[name]
            ends = _ends_text(name, b.lo, b.hi)
        else:
            ends = _ends_text(name, None, None)
        sk = [res_k.get((name, e), ref_k) for e in ("low", "high")]
        st = [res_t.get((name, e), ref_t) for e in ("low", "high")]
        lines.append(
            f"{UNC_LABELS[name]} & {ends} & "
            f"{sk[0]['emission_vuv']:.2f} to {sk[1]['emission_vuv']:.2f} & "
            f"{st[0]['emission_vuv']:.2f} to {st[1]['emission_vuv']:.2f} & "
            f"{sk[0]['emission_nir']:.2f} to {sk[1]['emission_nir']:.2f} & "
            f"{sk[0]['tg_end_K']:.0f} to {sk[1]['tg_end_K']:.0f}\\\\")
    return n, ["\\def\\wgUncTable{" + "\n".join(lines) + "}"]


# ------------------------------------------------------ the reference
def reference_numbers() -> tuple:
    n, lines = {}, []
    ref, ref2 = K.run_many([("field", K.ref_params()),
                            ("field", K.ref_params("two_t"))])
    for s, tag in ((ref, "kin"), (ref2, "2t")):
        ev = s["events"]
        krs = s["res"]["Kr*"]
        n.update({
            f"ref.{tag}.147": f"{s['emission_vuv']:.2f}",
            f"ref.{tag}.IR": f"{s['emission_nir']:.2f}",
            f"ref.{tag}.heat": f"{s['heat_reactions']:.2f}",
            f"ref.{tag}.elastic": f"{s['elastic']:.2f}",
            f"ref.{tag}.krstar": f"{s['to_Krstar']:.2f}",
            f"ref.{tag}.krplus": f"{s['to_Krplus']:.2f}",
            f"ref.{tag}.kr2": f"{s['to_Kr2star']:.2f}",
            f"ref.{tag}.ions": f"{s['stored_ions']:.2f}",
            f"ref.{tag}.te": f"{s['te_mean_eV']:.2f}",
            f"ref.{tag}.tg": f"{s['tg_end_K']:.0f}",
            f"ref.{tag}.ne": _sci(s["ne_end_cm3"]),
            f"ref.{tag}.cyc": f"{ev['K76'] / ev['AK8_Kr']:.0f}",
            f"ref.{tag}.brs76": f"{100 * krs['losses']['K76'] / krs['formed']:.1f}",
            f"ref.{tag}.ak8": f"{100 * krs['losses']['AK8_Kr'] / krs['formed']:.1f}",
            f"ref.{tag}.closure": _sci(abs(s["closure"]), 0),
            f"ref.{tag}.min": f"{s['seconds'] / 60:.0f}",
        })
        # the cycle's heat: the net pumping's 1.388 eV per event
        gap = e_form("Kr**") - e_form("Kr*")
        n[f"ref.{tag}.cycheat"] = f"{100 * (ev['K76'] - ev['K76r']) * gap / s['e_dep_eV']:.1f}"
        n[f"ref.{tag}.perphoton"] = f"{(ev['K76'] - ev['K76r']) * gap / ev['KR2_RAD']:.0f}"
    kin = ref["kin"]
    k15 = next(r for r in build() if r.rid == "K15").rate[1]
    n_kr = ref["n0"]["Kr"]
    n["k15.ns"] = f"{1e9 / (k15 * n_kr ** 2):.0f}"
    ev = ref["events"]
    hm_ev = sum(v for r, v in ev.items() if r.startswith("HM_"))
    c2 = (ev["C2"] + hm_ev) / (sum(ev[r] for r in ("C2", "G2", "KSTEP_S",
                                                   "KSTEP_P", "PEN_KR2"))
                               + hm_ev)
    if not c2 < 0.01:
        raise AssertionError("the text says Penning ionization makes the "
                             "reference point's ions")
    n.update({
        "ref.kin.n": str(kin["solutions"]), "ref.kin.on": str(kin["on"]),
        "ref.kin.en.lo": f"{kin['en_min']:.2f}",
        "ref.kin.en.hi": f"{kin['en_max']:.1f}",
        "ref.kin.ratio.lo": f"{kin['ratio_min']:.4f}",
        "ref.kin.ratio.hi": f"{kin['ratio_max']:.4f}",
        "ref.kin.gap": _sci(kin["balance_max"], 0),
        "ref.kin.reach": _sci(kin["reach_max"], 0),
        "ref.kin.emax.lo": f"{kin['emax_min']:.0f}",
        "ref.kin.emax.hi": f"{kin['emax_max']:.0f}",
        "ref.kin.vs2t": _pm(100 * (ref["emission_vuv"] / ref2["emission_vuv"]
                                   - 1)),
    })
    # the 5s -> 5p scale and Kannari's fit, both closures
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        jobs = [("field", K.ref_params(eedf).replace(
            rates=(("k_5s5p", k),))) if k != 1.0 else
            ("field", K.ref_params(eedf)) for k in K.K5S5P_TABLE]
        got = K.run_many(jobs)
        for k, s in zip(K.K5S5P_TABLE, got):
            n[f"k5.{k:g}.{tag}"] = f"{s['emission_vuv']:.1f}"
            n[f"k5.{k:g}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
        sw = K.switch_params(eedf)
        kann = K.run_many([("field", sw["kannari"])])[0]
        n[f"kann.{tag}"] = f"{kann['emission_vuv']:.1f}"
    # the scan, two-temperature closure
    sp = K.scan_params()
    keys = list(sp)
    got = dict(zip(keys, K.run_many([("field", sp[k]) for k in keys])))
    for w in K.SCAN_KW:
        cells = [f"{got[(a, w)]['emission_vuv']:.1f}" for a in K.SCAN_ATM]
        lines.append(f"{w:.0f} & " + " & ".join(cells) + "\\\\")
    tg01 = [got[(0.1, w)]["tg_end_K"] for w in K.SCAN_KW]
    n["scan.tg01.lo"] = f"{min(tg01):.0f}"
    n["scan.tg01.hi"] = f"{max(tg01):.0f}"
    n["scan.best"] = f"{max(s['emission_vuv'] for s in got.values()):.0f}"
    return n, ["\\def\\wgScan{" + "\n".join(lines) + "}"]


# ------------------------------------------------- numerical accuracy
def convergence_numbers() -> dict:
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    got = K.run_many([("field", p) for *_, p in K.conv_params()])
    worst = 0.0
    for (k, v, _), s in zip(K.conv_params(), got):
        d147 = abs(s["emission_vuv"] / ref["emission_vuv"] - 1)
        dir_ = abs(s["emission_nir"] / ref["emission_nir"] - 1)
        n[f"conv.{k}.147"] = _sci(d147, 0)
        n[f"conv.{k}.IR"] = _sci(dir_, 0)
        if k != "rtol":
            worst = max(worst, d147, dir_)
        else:
            n["conv.rtol.147.pct"] = f"{100 * d147:.2f}"
            # v0_2: the renewals of the tighter run, against the reference's
            n["conv.rtol.n"] = str(s["kin"]["solutions"])
    n["conv.worst"] = _sci(worst, 0)
    # v0_3: the kinetic solver's cold restarts over every run of this note
    # (boltzmann.Solver.eedf: a warm start from which no EEDF was found)
    kin = [s["kin"] for s in K.run_many(K.all_jobs()) if "kin" in s]
    n["kin.runs"] = str(len(kin))
    times = sum(k["restarts"] for k in kin)
    runs = sum(1 for k in kin if k["restarts"])
    if not times:
        raise AssertionError("the text says the solver restarted")
    n["kin.restarts"] = ({1: "once", 2: "twice"}.get(times, f"{times} times")
                         + f", in {WORDS[runs] if runs < 10 else runs} "
                         + ("run" if runs == 1 else "runs"))
    return n


def change_numbers() -> dict:
    """v0_2 (Biagi's file as published, its own cache, kr_docs_v0_2)
    against v0_3 at the reference point on both closures, and the
    network's size."""
    n = {}
    old = K2.run_many([("field", K2.ref_params()),
                       ("field", K2.ref_params("two_t"))])
    new = K.run_many([("field", K.ref_params()),
                      ("field", K.ref_params("two_t"))])
    for (o, s), tag in zip(zip(old, new), ("kin", "2t")):
        n[f"chg.{tag}.v02"] = f"{o['emission_vuv']:.2f}"
        n[f"chg.{tag}.v03"] = f"{s['emission_vuv']:.2f}"
        n[f"chg.{tag}.dev"] = _pm(100 * (s["emission_vuv"]
                                         / o["emission_vuv"] - 1), 1)
        n[f"chg.{tag}.te.v02"] = f"{o['te_mean_eV']:.2f}"
        n[f"chg.{tag}.te.v03"] = f"{s['te_mean_eV']:.2f}"
        n[f"chg.{tag}.el.v02"] = f"{o['elastic']:.2f}"
        n[f"chg.{tag}.el.v03"] = f"{s['elastic']:.2f}"
        n[f"chg.{tag}.heat.v02"] = f"{o['heat_reactions']:.2f}"
        n[f"chg.{tag}.heat.v03"] = f"{s['heat_reactions']:.2f}"
    n["par.rows"] = str(len(build()))
    n["par.levels"] = str(len(levels()))
    return n


def jung_numbers() -> dict:
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    te_k = ref["te_mean_eV"]
    table = jung2006.rate_table()
    for te, key in ((1.0, "1"), (2.0, "2"), (3.0, "3"), (te_k, "ref")):
        r = float(table(te)) / float(jung2006.kannari_rate(te))
        n[f"jung.ratio.{key}"] = f"{r:.2f}"
    n["jung.te.ref"] = f"{te_k:.2f}"
    n["jung.sig6"] = _sci(float(jung2006.lumped_sigma(6.0)))
    n["jung.sig6.meta"] = _sci(float(jung2006.lumped_sigma(6.0,
                                                           proxy_j1=False)))
    n["jung.unc.1s5"], n["jung.unc.1s3"] = (str(x) for x in JUNG_UNC)
    n["jung.k.ref"] = _sci(float(table(te_k)))
    return n


def ladder_numbers() -> dict:
    """The energy ladder's levels and the energies its rows move [eV]."""
    rows = {r.rid: r for r in build()}
    E = {x: e_form(x) for x in ("Kr*", "Kr**", "Kr+", "Kr2*", "Kr2+")}
    hv = rows["KR2_RAD"].hv_eV
    return {
        "lad.S": f"{E['Kr*']:.3f}", "lad.P": f"{E['Kr**']:.2f}",
        "lad.I": f"{E['Kr+']:.3f}", "lad.X2": f"{E['Kr2*']:.2f}",
        "lad.I2": f"{E['Kr2+']:.3f}",
        "lad.hv": f"{hv:.2f}", "lad.low": f"{E['Kr2*'] - hv:.2f}",
        "lad.gap": f"{E['Kr**'] - E['Kr*']:.3f}",
        "lad.ak8": f"{E['Kr*'] - E['Kr2*']:.3f}",
        "lad.k15": f"{E['Kr+'] - E['Kr2+']:.2f}",
        "lad.r11": f"{E['Kr2+'] - E['Kr**']:.3f}",
        "lad.pen": f"{2 * E['Kr*'] - E['Kr+']:.3f}",
        "lad.c1hi": f"{max(r.e_cost_eV for r in group_rows(build(), 'C1c')):.2f}",
        "lad.c1lo": f"{min(r.e_cost_eV for r in group_rows(build(), 'C1a')):.2f}",
        "lad.s1": f"{min(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
        "lad.s2": f"{max(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
    }


def numbers() -> tuple:
    # every run at once, in parallel, so a stale cache is one batch
    K.run_many(K.all_jobs())
    n = {"version": K.V, "md5": K.full_md5()[:12]}
    tables = []
    for part in (reference_numbers, uncertainty_numbers, swarm_numbers):
        a, t = part()
        n.update(a)
        tables += t
    for part in (beam_numbers, convergence_numbers, change_numbers,
                 jung_numbers, ladder_numbers):
        n.update(part())
    n["tau5p"] = f"{KR_5P_TAU_S * 1e9:.1f}"
    return n, tables


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_3 overview numbers. GENERATED by make_overview_v0_3.py",
        f"% from the model ({n['version']}). Do not edit. Rerun the script "
        "instead.",
        "\\providecommand\\wgdef[2]{\\expandafter\\def"
        "\\csname wg@#1\\endcsname{#2}}",
    ]
    lines += [f"\\wgdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_3's accuracy "
                                 "overview numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_overview_v0_3.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
