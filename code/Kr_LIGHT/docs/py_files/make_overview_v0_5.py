"""Numbers for 'How accurate LIGHT_Kr v0_5 is' (accuracy_overview_v0_5.tex).

Writes overview_v0_5_numbers.tex, one \\wgdef{key}{value} per number and a
few \\def tables, from the model. v0_5's cross sections and swarm are
v0_3's and v0_4's, so the swarm part is theirs, computed again here:

  the swarm      v0_5's two-term solver run as a swarm (v0_5/swarm.py) on its
                 trained cross sections and on Biagi's file as published,
                 as the experiments measure it (v0_5/swarm_mc.py), against
                 every krypton dataset of ML's registry (ML/ml/kr_data.py),
                 each in the role ML's consistency check gave it, with the
                 training's record (xsec.TRAINED_FILE).
  the e-beam     v0_5 under the beam drive (kr_docs_v0_5.BeamModel) against
                 Morozov et al. 2008 and Eckstrom et al. 1988, the drive's
                 partition against Saito et al. 2003, and v0_4 at the same
                 points (kr_docs_v0_4, its own cache).
  uncertainty    every band the rows read at its two ends and the three
                 measured constants with a quoted error at theirs, one at a
                 time, on both closures, and the trained cross sections as
                 draws of the training's posterior.
  the reference  the kinetic EEDF and the two-temperature closure, where the
                 energy goes (the 5s <-> 5p pumps at their own thresholds,
                 the electrons' ionizations of the excited atoms and of
                 Kr2*), the 5s -> 5p scale, Kannari's fit, the scan, the
                 kinetic closure's convergence, and v0_4 at the same points.
  v0_5's own     each v0_5 change undone alone (kr_docs_v0_5.undo_params);
                 the calculated rows against their calculation (Vriens &
                 Smeets's cross section against their rate formula,
                 Kannari's fits, Flannery et al.'s order); the superelastic
                 returns against Blagoev et al.'s afterglow rates; the
                 field tail table's excited axis against the kinetic EEDF's
                 own tail along its runs and the closure at the scan's
                 corners; Mano et al.'s electroluminescence yield.

    python3 LIGHT_Kr/docs/make_overview_v0_5.py        # minutes (cached)
    cd LIGHT_Kr/docs && latexmk -pdf accuracy_overview_v0_5.tex

The runs come from kr_docs_v0_5.run_many (cached per version of the
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

import kr_docs_v0_4 as K4                                  # noqa: E402
import kr_docs_v0_5 as K                                   # noqa: E402
import make_overview_v0_1 as OV1                           # noqa: E402
import make_overview_v0_2 as OV2                           # noqa: E402
from v0_5 import excited, jung2006, swarm, swarm_mc, thermo, xsec  # noqa
from v0_5.network import (BANDS, CALCULATED, E_ASSOC_EV,   # noqa: E402
                          G_LEVEL_5S, KR_5P_LINES, KR_5P_TAU_S, KR_5S,
                          LEVELS_5S, OFF_BY_DEFAULT, build, group_rows,
                          kr5s, levels)
from v0_5.species import e_form                            # noqa: E402
from v0_5.report import KR2_TG_VALID_K                    # noqa: E402

OUT = HERE / "overview_v0_5_numbers.tex"

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
# what the parts computed, for the text's claims (claims())
_CLAIMS: dict = {}


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
G2_ROWS = tuple(r.rid for r in build() if r.rid.startswith("G2_"))
AK8F = ("AK8_1s5", "AK8_1s4", "AK8f_1s4", "AK8_1s3")
AK8R = ("AK8r_1s5", "AK8r_1s4", "AK8fr_1s4", "AK8r_1s3")
# v0_5: the 5s -> 5p pumps pair by pair, and their reverses
K76F = tuple(r.rid for r in build() if r.rid.startswith("K76_"))
K76R = tuple(r.rid for r in build() if r.rid.startswith("K76r_"))
RAD_ROWS = ("KR2_RAD", "KR2F_RAD")
KSTEP_S = tuple(f"KSTEP_S_{lv}" for lv in LEVELS_5S)
KSTEP_P = tuple(r.rid for r in build() if r.rid.startswith("KSTEP_P_"))
EI_ROWS = ("EI_KR2", "EI_KR2F")
SUPER_ROWS = tuple(f"SUPER_{lv}" for lv in LEVELS_5S)
KPSR_ROWS = tuple(r.rid for r in build() if r.rid.startswith("KPSr_"))
PRAD_ROWS = tuple(r.rid for r in build() if r.rid.startswith("PRAD_"))


def _sum(ev, rids):
    return sum(ev[r] for r in rids)


def beam_numbers() -> dict:
    n = {}
    Y = K.beam_yields()
    # v0_3's overview compared the trained set's degradation with Biagi's
    # file as published (v0_2's): the same numbers, from v0_2's documents
    import kr_docs_v0_2 as K2
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
    # v0_4: the 5s yields level by level
    for lv in LEVELS_5S:
        n[f"beam.s.{lv}"] = f"{Y[f's_{lv}']:.3f}"
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
            n[f"{name}.{tag}.G2"] = f"{_sum(ev, G2_ROWS) / pairs:.3f}"
            exc = sum(v for r, v in ev.items()
                      if r.startswith(("G1_", "C1")))
            n[f"{name}.{tag}.G1"] = f"{exc / pairs:.3f}"
            n[f"{name}.{tag}.te"] = f"{s['te_mean_eV']:.2f}"
            n[f"{name}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
            n[f"{name}.{tag}.ph"] = f"{_sum(ev, RAD_ROWS) / pairs:.3f}"
            n[f"{name}.{tag}.rec"] = f"{ev['R11p'] / pairs:.3f}"
            n[f"{name}.{tag}.krs"] = _sci(s["n_max_krstar"])
            # v0_4: the excimers redissociated, per pair
            n[f"{name}.{tag}.redis"] = f"{_sum(ev, AK8R) / pairs:.3f}"
            n[f"{name}.{tag}.ne"] = _sci(s["ne_peak"])
            n[f"{name}.{tag}.closure"] = _sci(abs(s["closure"]), 0)
            # the loss below v0_3's ceiling (one photon per excitation and
            # per ion), the quantity a report 'near the limit' tests
            loss = f"{100 * (1 - vuv / ceil):.0f}"
            n[f"{name}.{tag}.loss"] = "0" if loss == "-0" else loss
        m = got[("morozov", K.MOROZOV_W[1])]
        v = m["emission_vuv"] / 100.0
        n[f"moz.{tag}.dev"] = _pm(100 * (v / MOROZOV[0] - 1))
        n[f"moz.{tag}.sig"] = f"{(v - MOROZOV[0]) / MOROZOV[1]:.1f}"
        e1 = got[("eckstrom", 1.0)]["emission_vuv"] / 100.0
        n[f"eck.{tag}.dev"] = _pm(100 * (e1 / ECKSTROM - 1))
        # what the text says of these runs
        if not v > MOROZOV[0]:
            raise AssertionError("the text says v0_5 lies above Morozov")
        if not 1 - got[("eckstrom", 2.0)]["emission_vuv"] / 100 / ceil < 0.10:
            raise AssertionError("the text says v0_5 loses little at 2 atm")
        if not abs(e1 / ECKSTROM - 1) < 0.05:
            raise AssertionError(f"the text says v0_5 meets Eckstrom at 1 atm:"
                                 f" {e1:.3f}")
        if not (got[("eckstrom", 0.5)]["tg_end_K"] > 1000.0):
            raise AssertionError("the text says the 0.5-atm run heats past "
                                 "1000 K")
        # v0_4: the high-pump account is written from these facts
        n[f"eck.{tag}.rising"] = str(int(
            got[("eckstrom", 0.5)]["emission_vuv"] < 100 * e1
            < got[("eckstrom", 2.0)]["emission_vuv"]))
        n[f"eck.{tag}.below"] = str(int(e1 < ECKSTROM))
        if not (got[("eckstrom_pen", 1.0)]["emission_vuv"] > 100 * e1):
            raise AssertionError("the text says a weaker Kr* + Kr* raises "
                                 "the high-pump efficiency")
        n[f"moz.{tag}.flat"] = _pm(100 * (
            got[("morozov", K.MOROZOV_W[0])]["emission_vuv"]
            / got[("morozov", K.MOROZOV_W[1])]["emission_vuv"] - 1), 1)
    # v0_4 at the same points, from its own cache (kr_docs_v0_4)
    bp4 = K4.beam_params("auto")
    keys4 = [("morozov", K4.MOROZOV_W[1])] + [("eckstrom", a)
                                              for a in K4.ECKSTROM_ATM]
    got4 = dict(zip(keys4, K4.run_many([("beam", bp4[k]) for k in keys4])))
    m4 = got4[("morozov", K4.MOROZOV_W[1])]
    n["moz.kin.dev.v04"] = _pm(100 * (m4["emission_vuv"] / 100 / MOROZOV[0]
                                      - 1))
    n["moz.kin.v04"] = f"{m4['emission_vuv'] / 100:.3f}"
    for a in K4.ECKSTROM_ATM:
        e4 = got4[("eckstrom", a)]["emission_vuv"] / 100
        n[f"eck.{a:g}.kin.v04"] = f"{e4:.3f}"
    n["eck.kin.dev.v04"] = _pm(100 * (got4[("eckstrom", 1.0)]["emission_vuv"]
                                      / 100 / ECKSTROM - 1))
    # v0_5 against v0_4 under the beam: at low pump and at each of
    # Eckstrom's pressures
    bp5 = K.beam_params("auto")
    m5 = K.run_many([("beam", bp5[("morozov", K.MOROZOV_W[1])])])[0]
    _CLAIMS["moz_same"] = (round(m5["emission_vuv"] / 100, 3)
                           == round(m4["emission_vuv"] / 100, 3))
    pairs5 = m5["e_dep_eV"] / m5["beam"]["W"]
    _CLAIMS["moz_ei"] = _sum(m5["events"], EI_ROWS) / pairs5
    _CLAIMS["moz_ph"] = _sum(m5["events"], RAD_ROWS) / pairs5
    n["moz.300.kin.ei"] = f"{_CLAIMS['moz_ei']:.3f}"
    for a in K4.ECKSTROM_ATM:
        e4 = got4[("eckstrom", a)]["emission_vuv"]
        e5 = K.run_many([("beam", bp5[("eckstrom", a)])])[0]["emission_vuv"]
        n[f"eck.{a:g}.vs.v04"] = _pm(100 * (e5 / e4 - 1))
        _CLAIMS[f"eck_vs_v04_{a:g}"] = e5 / e4 - 1
    # v0_5: what its electron rows do under the beam, whose electrons are
    # cold: the stepwise and the excimer's ionizations per pair at 1 atm
    s1 = None
    for eedf in ("auto",):
        bp = K.beam_params(eedf)
        s1 = K.run_many([("beam", bp[("eckstrom", 1.0)])])[0]
    pairs = s1["e_dep_eV"] / s1["beam"]["W"]
    ev = s1["events"]
    n["eck.1.kin.step"] = f"{(_sum(ev, KSTEP_S) + _sum(ev, KSTEP_P)) / pairs:.3f}"
    n["eck.1.kin.ei"] = f"{_sum(ev, EI_ROWS) / pairs:.3f}"
    return n


# ------------------------------------------------------ the uncertainty
UNC_LABELS = {
    "k_5s5p": r"\el{} + \Krs{} \ra \Krpp{}, scale",
    "k_5s5p_j1": r"\el{} + \Krs($J=1$) \ra \Krpp{}, scale",
    "kps_kr": r"\Krpp{} + Kr \ra \Krs{} + Kr",
    "kps_5s": r"its products over the 5s levels, share",
    "kr_pen": r"\Krs{} + \Krs{}, Penning",
    "kr2_form": r"\Krs(1s$_5$) + 2Kr \ra \Krts{} + Kr",
    "kr2_form_1s4": r"\Krs(1s$_4$) + 2Kr \ra \Krts{} + Kr",
    "kr2_form_1s3": r"\Krs(1s$_3$) + 2Kr \ra \Krts{} + Kr",
    "kr2_fast": r"the 1s$_4$ formation's share into 0$_u^+$",
    "kr2_quench": r"\Krts{} + Kr \ra 3Kr",
    "kr2_d0": r"the excimer's well (its redissociation)",
    "dr_tg": r"\Krtp{} + \el{}, the gas-temperature factor",
    "res_escape": r"the resonance lines' escape",
    "kr2_pen": r"\Krts{} + \Krts{}, Penning",
    "kstep_s": r"\el{} + \Krs{} \ra \Krp{} + 2\el{}, scale",
    "kstep_p": r"\el{} + \Krpp{} \ra \Krp{} + 2\el{}, scale",
    "ei_kr2": r"\el{} + \Krts{} \ra \Krtp{} + 2\el{}, scale",
    "eq_kr2": r"\el{} + \Krts{} \ra \el{} + 2Kr",
    "em_kr2": r"\el{} + \Krts{} \ra \el{} + \Krts(0$_u^+$)",
    "emix_5s": r"\el{} + \Krs(1s$_5$) \ra \el{} + \Krs(1s$_4$)",
    "hm_kr": r"Kr(nl) + Kr \ra \Krtp{} + \el{}, the share",
    "R11p": r"\Krtp{} + \el{} \ra \Krpp{} + Kr",
    "K15": r"\Krp{} + 2Kr \ra \Krtp{} + Kr",
    "KR2_RAD": r"\Krts{} lifetime",
}
UNC_SHORT = {
    "k_5s5p": r"the 5s \ra 5p cross sections", "kps_kr": r"\Krpp{} + Kr",
    "k_5s5p_j1": r"the $J=1$ levels' pumping",
    "kps_5s": r"the 5p relaxation's products",
    "kr2_form_1s4": r"the 1s$_4$ formation",
    "kr2_form_1s3": r"the 1s$_3$ formation",
    "kr2_fast": r"the fast excimer's share",
    "kr2_quench": r"\Krts{} + Kr", "kr2_d0": r"the excimer's well",
    "dr_tg": r"the recombination's gas-temperature factor",
    "res_escape": r"the resonance lines",
    "kr_pen": r"\Krs{} + \Krs{}", "kr2_form": r"the 1s$_5$ formation",
    "kr2_pen": r"\Krts{} + \Krts{}",
    "kstep_s": r"the stepwise ionization of the 5s levels",
    "kstep_p": r"the stepwise ionization of \Krpp",
    "ei_kr2": r"the ionization of \Krts{} by electrons",
    "eq_kr2": r"the quenching of \Krts{} by electrons",
    "em_kr2": r"the electrons' mixing of the excimer's states",
    "emix_5s": r"the electrons' 1s$_5$ \ra 1s$_4$",
    "R11p": r"\Krtp{} + \el",
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
        if b.unit == "eV":
            return f"{lo:g} to {hi:g}\\unit{{eV}}"
        unit = {"cm3/s": r"\unit{cm^3\,s^{-1}}",
                "cm6/s": r"\unit{cm^6\,s^{-1}}",
                "1/s": r"\unit{s^{-1}}"}[b.unit]
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
    # the rows with no measurement or set by a calculation, and those off
    # by default, at their band's top and (the calculated) bottom, as
    # changes of the share
    for name in ("kr2_pen", "kstep_s", "kstep_p", "ei_kr2", "eq_kr2",
                 "em_kr2", "kr2_quench", "kr2_fast", "dr_tg", "res_escape",
                 "emix_5s", "hm_kr"):
        for end in ("low", "high"):
            s = res_k.get((name, end), ref_k)
            d = 100 * (s["emission_vuv"] / ref_k["emission_vuv"] - 1)
            n[f"unc.{end}.{name}"] = _pm(d)
            n[f"unc.{end}.{name}.v"] = f"{s['emission_vuv']:.2f}"
    # every band's own ends, kinetic
    for name in [b.name for b in K.row_bands()]:
        for end in ("low", "high"):
            s = res_k.get((name, end), ref_k)
            n[f"unc.end.{name}.{end}"] = f"{s['emission_vuv']:.2f}"
    # v0_2 to v0_4: the Hornbeck-Molnar share at its top barely moves the
    # reference point
    hm_cut = res_k[("hm_kr", "high")]["emission_vuv"] / ref_k["emission_vuv"] - 1
    n["unc.hm.cut"] = _pm(100 * hm_cut, 2)

    def term_k(x):
        return 0.5 * abs(res_k.get((x, "high"), ref_k)["emission_vuv"]
                         - res_k.get((x, "low"), ref_k)["emission_vuv"])
    order = sorted(names, key=lambda x: -term_k(x))
    _CLAIMS["unc_order"] = order
    _CLAIMS["unc_terms"] = {x: term_k(x) for x in names}
    _CLAIMS["unc_ref"] = ref_k["emission_vuv"]
    _CLAIMS["unc_res"] = res_k
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
    rows = {r.rid: r for r in build()}
    for s, tag in ((ref, "kin"), (ref2, "2t")):
        ev = s["events"]
        e_dep = s["e_dep_eV"]
        krs = s["res"]["Kr*"]
        k2s = s["res"]["Kr2*"]
        kpp = s["res"]["Kr**"]
        form, redis = _sum(ev, AK8F), _sum(ev, AK8R)
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
            f"ref.{tag}.ne.on": _sci(s["ne_mean_on"]),
            f"ref.{tag}.cyc": f"{_sum(ev, K76F) / (form - redis):.0f}",
            f"ref.{tag}.ak8": f"{100 * sum(krs['losses'].get(r, 0.0) for r in AK8F) / krs['formed']:.1f}",
            f"ref.{tag}.redis": f"{100 * redis / form:.1f}",
            f"ref.{tag}.k2rad": f"{100 * sum(k2s['losses'].get(r, 0.0) for r in RAD_ROWS) / k2s['formed']:.1f}",
            # v0_5: the excimers the electrons ionize, and those formed
            f"ref.{tag}.k2ei": f"{100 * sum(k2s['losses'].get(r, 0.0) for r in EI_ROWS) / k2s['formed']:.1f}",
            f"ref.{tag}.k2redis": f"{100 * sum(k2s['losses'].get(r, 0.0) for r in AK8R) / k2s['formed']:.1f}",
            f"ref.{tag}.closure": _sci(abs(s["closure"]), 0),
            f"ref.{tag}.min": f"{s['seconds'] / 60:.0f}",
            # v0_5: what the 5s atoms' electrons do: their stepwise
            # ionization and their return to the ground state, per 5s atom
            # formed
            f"ref.{tag}.krs.step": f"{100 * sum(krs['losses'].get(r, 0.0) for r in KSTEP_S) / krs['formed']:.1f}",
            f"ref.{tag}.krs.super": f"{100 * sum(krs['losses'].get(r, 0.0) for r in SUPER_ROWS) / krs['formed']:.2f}",
            f"ref.{tag}.kpp.step": f"{100 * sum(kpp['losses'].get(r, 0.0) for r in KSTEP_P) / kpp['formed']:.1f}",
            # the pumps' energy above the rung's 2p10, which v0_4 did not
            # charge, as a share of the input
            f"ref.{tag}.pumpx": f"{100 * s['pump_excess_eV'] / e_dep:.1f}",
        })
        for lv in ("1s5", "1s4", "1s3"):
            f_lv = sum(ev[r] for r in AK8F if r.endswith(lv))
            r_lv = sum(ev[r] for r in AK8R if r.endswith(lv))
            n[f"ref.{tag}.form.{lv}"] = f"{100 * f_lv / form:.1f}"
        # the pumping by level, percent of all 5s -> 5p
        k76 = _sum(ev, K76F)
        for lv in LEVELS_5S:
            k_lv = sum(v for r, v in ev.items() if r.startswith(f"K76_{lv}_"))
            n[f"ref.{tag}.pump.{lv}"] = f"{100 * k_lv / k76:.1f}"
        # the cycle's heat: each net pump's threshold, pair by pair
        heat = sum(ev[r] * rows[r].e_cost_eV for r in K76F) - sum(
            ev[r] * rows[r].e_gain_eV for r in K76R)
        n[f"ref.{tag}.cycheat"] = f"{100 * heat / e_dep:.1f}"
        n[f"ref.{tag}.perphoton"] = f"{heat / _sum(ev, RAD_ROWS):.0f}"
        # the ions' sources: direct and associative ionization, Penning
        # (5s on 5s and Kr2* on Kr2*), stepwise out of 5s and 5p, Kr2*'s
        src = {"C2": ev["C2"] + sum(v for r, v in ev.items()
                                    if r.startswith("HM_")),
               "pen": _sum(ev, G2_ROWS) + ev["PEN_KR2"],
               "step_s": _sum(ev, KSTEP_S), "step_p": _sum(ev, KSTEP_P),
               "ei": _sum(ev, EI_ROWS)}
        tot = sum(src.values())
        for k, v in src.items():
            n[f"ref.{tag}.ion.{k}"] = f"{100 * v / tot:.1f}"
        _CLAIMS[f"ion_{tag}"] = {k: v / tot for k, v in src.items()}
        _CLAIMS[f"ref_{tag}"] = s
    kin = ref["kin"]
    k15 = next(r for r in build() if r.rid == "K15").rate[1]
    n_kr = ref["n0"]["Kr"]
    n["k15.ns"] = f"{1e9 / (k15 * n_kr ** 2):.0f}"
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
    # the excimer's bound share at the fill's density (thermo.k_eq) and
    # the gas temperature where it falls to one
    n_fill = ref["n0"]["Kr"]
    d0 = e_form(kr5s("1s5")) - e_form("Kr2*")
    for tg in (300.0, 700.0, 1000.0, 1500.0):
        ratio = thermo.k_eq("Kr2*", d0, G_LEVEL_5S["1s5"], tg) * n_fill
        n[f"keq.{tg:.0f}"] = (f"{ratio:.2g}" if ratio < 1e4 else
                              _sci(ratio, 0))
    lo, hi = 300.0, 3000.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if thermo.k_eq("Kr2*", d0, 5, mid) * n_fill > 1.0:
            lo = mid
        else:
            hi = mid
    n["keq.half.K"] = f"{lo:.0f}"
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
        cells = [f"{got[(a, w)]['emission_vuv']:.2f}" for a in K.SCAN_ATM]
        lines.append(f"{w:.0f} & " + " & ".join(cells) + "\\\\")
    tg01 = [got[(0.1, w)]["tg_end_K"] for w in K.SCAN_KW]
    n["scan.tg01.lo"] = f"{min(tg01):.0f}"
    n["scan.tg01.hi"] = f"{max(tg01):.0f}"
    n["scan.best"] = f"{max(s['emission_vuv'] for s in got.values()):.0f}"
    best = max(got, key=lambda k: got[k]["emission_vuv"])
    n["scan.best.atm"] = f"{best[0]:g}"
    n["scan.best.kw"] = f"{best[1]:.0f}"
    ions01 = got[(0.1, 500.0)]["stored_ions"]
    n["scan.ions01"] = f"{ions01:.0f}"
    # v0_5: the 5p rung's reverse at the scan's points (its share of the
    # relaxation), where the gas is hottest
    sh = {}
    for k, s in got.items():
        kps = sum(v for r, v in s["events"].items()
                  if r.startswith("KPS_"))
        kpsr = _sum(s["events"], KPSR_ROWS)
        sh[k] = kpsr / kps if kps > 0 else 0.0
    kmax = max(sh, key=sh.get)
    n["kpsr.max"] = f"{100 * sh[kmax]:.1f}"
    n["kpsr.max.atm"] = f"{kmax[0]:g}"
    n["kpsr.max.kw"] = f"{kmax[1]:.0f}"
    n["kpsr.ref"] = _sci(sh[(1.0, 300.0)], 0)
    _CLAIMS["scan"] = got
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
    # v0_4: the afterglow's EEDF time-stepped, against its steady shape
    # (Params.kin_afterglow_memory): the bound the groundwork measured
    ag = K.run_many([("field", K.switch_params()["afterglow"])])[0]
    n["ag.147"] = _sci(abs(ag["emission_vuv"] / ref["emission_vuv"] - 1), 0)
    n["ag.IR"] = _sci(abs(ag["emission_nir"] / ref["emission_nir"] - 1), 0)
    n["ag.heat"] = _sci(abs(ag["heat_reactions"] / ref["heat_reactions"]
                            - 1), 0)
    n["ag.n"] = str(ag["kin"]["stepped_after"])
    if not ag["kin"]["stepped_after"] > 0:
        raise AssertionError("the afterglow run stepped no EEDF")
    if not abs(ag["emission_vuv"] / ref["emission_vuv"] - 1) < 1e-3:
        raise AssertionError("the text says the afterglow's shape moves the "
                             "share by under 1e-3")
    # v0_3: the kinetic solver's cold restarts over every run of this note
    # (boltzmann.Solver.eedf: a warm start from which no EEDF was found)
    kin = [s["kin"] for s in K.run_many(K.all_jobs()) if "kin" in s]
    n["kin.runs"] = str(len(kin))
    times = sum(k["restarts"] for k in kin)
    runs = sum(1 for k in kin if k["restarts"])
    n["kin.restarts"] = ("never" if not times else
                         {1: "once", 2: "twice"}.get(times, f"{times} times")
                         + f", in {WORDS[runs] if runs < 10 else runs} "
                         + ("run" if runs == 1 else "runs"))
    # v0_4: the fallbacks (boltzmann.Solver.eedf), and in which runs
    fb = sum(k.get("fallbacks", 0) for k in kin)
    fb_runs = sum(1 for k in kin if k.get("fallbacks", 0))
    n["kin.fallbacks"] = ("never" if not fb else
                          {1: "once", 2: "twice"}.get(fb, f"{fb} times")
                          + f", in {WORDS[fb_runs] if fb_runs < 10 else fb_runs} "
                          + ("run" if fb_runs == 1 else "runs"))
    # v0_5: table points left out, and updates LIGHT's Maxwellian took
    # (model._kinetic_update), and in which runs
    for key in ("dropped", "failed"):
        tot = sum(k.get(key, 0) for k in kin)
        nr = sum(1 for k in kin if k.get(key, 0))
        n[f"kin.{key}"] = ("never" if not tot else
                           {1: "once", 2: "twice"}.get(tot, f"{tot} times")
                           + f", in {WORDS[nr] if nr < 10 else nr} "
                           + ("run" if nr == 1 else "runs"))
    return n


def change_numbers() -> dict:
    """v0_4 (its own cache, kr_docs_v0_4) against v0_5 at the reference
    point on both closures, and the network's size."""
    n = {}
    old = K4.run_many([("field", K4.ref_params()),
                       ("field", K4.ref_params("two_t"))])
    new = K.run_many([("field", K.ref_params()),
                      ("field", K.ref_params("two_t"))])
    for (o, s), tag in zip(zip(old, new), ("kin", "2t")):
        n[f"chg.{tag}.v04"] = f"{o['emission_vuv']:.2f}"
        n[f"chg.{tag}.v05"] = f"{s['emission_vuv']:.2f}"
        n[f"chg.{tag}.dev"] = _pm(100 * (s["emission_vuv"]
                                         / o["emission_vuv"] - 1), 0)
        n[f"chg.{tag}.ratio"] = f"{o['emission_vuv'] / s['emission_vuv']:.1f}"
        for key, k in (("te", "te_mean_eV"), ("el", "elastic"),
                       ("heat", "heat_reactions"), ("tg", "tg_end_K"),
                       ("IR", "emission_nir")):
            fmt = ".0f" if key == "tg" else ".2f"
            n[f"chg.{tag}.{key}.v04"] = f"{o[k]:{fmt}}"
            n[f"chg.{tag}.{key}.v05"] = f"{s[k]:{fmt}}"
        n[f"chg.{tag}.ne.v04"] = _sci(o["ne_mean_on"])
        n[f"chg.{tag}.ne.v05"] = _sci(s["ne_mean_on"])
    # v0_5 on Biagi's file as published: what the training does to it
    bi = K.run_many(K.biagi_jobs())
    for s, tag, ref in zip(bi, ("kin", "2t"), new):
        n[f"chg.{tag}.bi"] = f"{s['emission_vuv']:.2f}"
        n[f"chg.{tag}.bi.dev"] = _pm(100 * (ref["emission_vuv"]
                                            / s["emission_vuv"] - 1), 1)
    n["par.rows"] = str(len(build()))
    n["par.rows.v04"] = str(len(K4.build()))
    n["par.levels"] = str(len(levels()))
    n["par.species"] = str(len(K.Model(K.ref_params("two_t")).species))
    n["par.pops"] = str(len(K.Model(K.ref_params()).kin.pops))
    n["par.lines"] = str(len(KR_5P_LINES))
    n["par.pairs"] = str(len(K76F))
    return n


def jung_numbers() -> dict:
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    te_k = ref["te_mean_eV"]
    table = jung2006.rate_table()
    for te, key in ((1.0, "1"), (2.0, "2"), (3.0, "3"), (te_k, "ref")):
        r = float(table(te)) / float(jung2006.kannari_rate(te))
        n[f"jung.ratio.{key}"] = f"{r:.2f}"
    # each level's rate at the reference point's Te, the pairs' sum against
    # v0_4's level table, and the rate-weighted energy per pump (v0_5)
    worst = 0.0
    for lv in LEVELS_5S:
        fs = jung2006.pair_finals(lv)
        k = {f: float(jung2006.pair_rate_table(lv, f)(te_k)) for f in fs}
        tot = sum(k.values())
        n[f"jung.k.{lv}"] = _sci(tot)
        worst = max(worst, abs(tot / float(jung2006.level_rate_table(lv)(te_k))
                               - 1))
        mean = sum(k[f] * jung2006.pair_threshold(lv, f) for f in fs) / tot
        n[f"jung.gap.{lv}"] = f"{mean:.3f}"
        n[f"jung.gap0.{lv}"] = f"{e_form('Kr**') - e_form(kr5s(lv)):.3f}"
        n[f"jung.npairs.{lv}"] = str(len(fs))
        top = max(k, key=k.get)
        n[f"jung.top.{lv}"] = top
        n[f"jung.top.{lv}.pct"] = f"{100 * k[top] / tot:.0f}"
    n["jung.sum.dev"] = _sci(worst, 0)
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
    E = {x: e_form(x) for x in ("Kr**", "Kr+", "Kr2*", "Kr2+", "Kr2*(0u+)")}
    E["Kr*"] = e_form(kr5s("1s5"))
    hv = rows["KR2_RAD"].hv_eV
    out = {f"lad.{lv}": f"{e_form(kr5s(lv)):.3f}" for lv in LEVELS_5S}
    out.update({f"lad.gap.{lv}": f"{E['Kr**'] - e_form(kr5s(lv)):.3f}"
                for lv in LEVELS_5S})
    # v0_5: the rung's levels, the 5p' (2P1/2-core) ones 0.8 eV up
    e0 = excited.RUNG["2p10"][0]
    out["lad.rung.top"] = f"{excited.RUNG['2p1'][0]:.3f}"
    out["lad.rung.5p.top"] = f"{excited.RUNG['2p5'][0] - e0:.2f}"
    out["lad.rung.5pp.lo"] = f"{excited.RUNG['2p4'][0] - e0:.2f}"
    out["lad.rung.5pp.hi"] = f"{excited.RUNG['2p1'][0] - e0:.2f}"
    for tg in (300.0, 700.0, 1000.0, 3000.0):
        out[f"lad.rung.f10.{tg:.0f}"] = f"{100 * excited.rung_fraction('2p10', tg):.0f}"
        out[f"lad.rung.z.{tg:.0f}"] = f"{excited.rung_partition(tg):.2f}"
    out["lad.rung.dE.700"] = f"{excited.rung_mean_excess(700.0):.3f}"
    a300 = sum(excited.rung_fraction(u, 300.0) * v
               for (u, _), v in KR_5P_LINES.items())
    out["lad.rung.tau300"] = f"{1e9 / a300:.1f}"
    a_inf = sum(excited.RUNG[u][1] * v for (u, _), v in KR_5P_LINES.items()) / 36.0
    out["lad.rung.tauinf"] = f"{1e9 / a_inf:.1f}"
    thr = {p.product: p.threshold for p in levels()}
    core = [thr[f"Kr(2P{i})"] - e_form(kr5s(lv)) for i in (1, 2, 3, 4)
            for lv in ("1s3", "1s2")]
    out["lad.core.lo"] = f"{min(core):.2f}"
    out["lad.core.hi"] = f"{max(core):.2f}"
    for lv in LEVELS_5S:
        out[f"lad.ion.{lv}"] = f"{E['Kr+'] - e_form(kr5s(lv)):.3f}"
    return {
        **out,
        "lad.F2": f"{E['Kr2*(0u+)']:.2f}",
        "lad.ak8f": f"{e_form(kr5s('1s4')) - E['Kr2*(0u+)']:.3f}",
        "lad.S": f"{E['Kr*']:.3f}", "lad.P": f"{E['Kr**']:.2f}",
        "lad.I": f"{E['Kr+']:.3f}", "lad.X2": f"{E['Kr2*']:.2f}",
        "lad.I2": f"{E['Kr2+']:.3f}",
        "lad.hv": f"{hv:.2f}", "lad.low": f"{E['Kr2*'] - hv:.2f}",
        "lad.gap": f"{E['Kr**'] - E['Kr*']:.3f}",
        "lad.ak8": f"{E['Kr*'] - E['Kr2*']:.3f}",
        "lad.k15": f"{E['Kr+'] - E['Kr2+']:.4f}",
        "lad.r11": f"{E['Kr2+'] - E['Kr**']:.3f}",
        "lad.pen": f"{2 * E['Kr*'] - E['Kr+']:.3f}",
        "lad.ei": f"{E['Kr2+'] - E['Kr2*']:.3f}",
        "lad.eiF": f"{E['Kr2+'] - E['Kr2*(0u+)']:.3f}",
        "lad.stepP": f"{E['Kr+'] - E['Kr**']:.3f}",
        "lad.mix2": f"{E['Kr2*(0u+)'] - E['Kr2*']:.3f}",
        "lad.mix15": f"{e_form(kr5s('1s4')) - e_form(kr5s('1s5')):.3f}",
        "lad.c1hi": f"{max(r.e_cost_eV for r in group_rows(build(), 'C1c')):.2f}",
        "lad.c1lo": f"{min(r.e_cost_eV for r in group_rows(build(), 'C1a')):.2f}",
        "lad.s1": f"{min(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
        "lad.s2": f"{max(r.e_cost_eV for r in group_rows(build(), 'G1')):.3f}",
    }


# ------------------------------------------------------ v0_5's changes
UNDO_TEXT = {
    "no_calc": r"the three calculated rows at 0 (v0\_4's rule)",
    "no_kstep": r"no stepwise ionization (5s and 5p)",
    "no_ei": r"no ionization of \Krts{} by electrons",
    "no_super": r"no superelastic return of the 5s levels",
    "no_emix": r"no 1s$_5$ \ra 1s$_4$ by electrons",
    "no_kpsr": r"no reverse of the 5p rung's relaxation",
    "ueda": r"every level forms the excimer at Ueda's "
            r"$1.7\times10^{-32}$ (v0\_3's)",
    "no_redissoc": r"no redissociation of \Krts",
    "no_j1": r"the $J=1$ levels pumped by no electron",
}


def undo_numbers() -> tuple:
    """Each v0_5 change undone alone from v0_5 (kr_docs_v0_5.undo_params)
    on both closures: the 147-nm share and the final gas temperature."""
    n, lines = {}, []
    got = {}
    refs = {}
    for eedf, tag in (("auto", "kin"), ("two_t", "2t")):
        ref = K.run_many([("field", K.ref_params(eedf))])[0]
        refs[tag] = ref
        up = K.undo_params(eedf)
        got[tag] = dict(zip(up, K.run_many([("field", p)
                                             for p in up.values()])))
        for k, s in got[tag].items():
            dev = 100 * (s["emission_vuv"] / ref["emission_vuv"] - 1)
            back = 100 * (ref["emission_vuv"] / s["emission_vuv"] - 1)
            n[f"undo.{k}.{tag}"] = f"{s['emission_vuv']:.2f}"
            # v0_5: a tenth of a percent where the change is under 10 %
            n[f"undo.{k}.{tag}.dev"] = _pm(dev, 1 if abs(dev) < 10 else 0)
            n[f"undo.{k}.{tag}.tg"] = f"{s['tg_end_K']:.0f}"
            n[f"undo.{k}.{tag}.ratio"] = f"{s['emission_vuv'] / ref['emission_vuv']:.2f}"
            # what the change itself does: v0_5 against it undone
            n[f"undo.{k}.{tag}.back"] = _pm(back, 1 if abs(back) < 10 else 0)
    for k, text in UNDO_TEXT.items():
        s, s2 = got["kin"][k], got["2t"][k]
        lines.append(f"{text} & {s['emission_vuv']:.2f} & "
                     f"{n[f'undo.{k}.kin.dev']} & {s2['emission_vuv']:.2f} & "
                     f"{s['tg_end_K']:.0f}\\\\")
    _CLAIMS["undo"] = got
    _CLAIMS["undo_ref"] = refs
    return n, ["\\def\\wgUndo{" + "\n".join(lines) + "}"]


def tail_numbers() -> tuple:
    """v0_5's field tail table with its excited axis against the kinetic
    EEDF's own tail along the kinetic runs (and its x_s = 0 block, v0_4's
    table), at the reference point and the scan's corners; and the
    two-temperature closure on the table, on its zero block and on v0_3's
    beam table against the kinetic EEDF at those points."""
    n, lines = {}, []
    ref = K.run_many([("field", K.ref_params())])[0]
    tl = ref["tail"]
    n.update({"tail.rk.lo": f"{tl['r_kin_lo']:.3f}",
              "tail.rk.hi": f"{tl['r_kin_hi']:.3f}",
              "tail.f.lo": f"{tl['field_lo']:.2f}",
              "tail.f.hi": f"{tl['field_hi']:.2f}",
              "tail.f.mean": f"{tl['field_mean']:.2f}",
              "tail.f0.lo": f"{tl['field0_lo']:.2f}",
              "tail.f0.hi": f"{tl['field0_hi']:.2f}",
              "tail.f0.mean": f"{tl['field0_mean']:.2f}",
              "tail.b.lo": f"{tl['beam_lo']:.2f}",
              "tail.b.hi": f"{tl['beam_hi']:.2f}",
              "tail.b.mean": f"{tl['beam_mean']:.2f}",
              "tail.xs.lo": _sci(tl["xs_lo"], 1),
              "tail.xs.hi": _sci(tl["xs_hi"], 1),
              "tail.n": str(tl["n"])})
    pts = list(K.SCAN_CORNERS) + [(1.0, 300.0)]
    kin = K.run_many([("field", K.scan_params("auto", points=[p])[p])
                      if p != (1.0, 300.0) else ("field", K.ref_params())
                      for p in pts])
    two_f = K.run_many([("field", K.scan_params("two_t", points=[p])[p])
                        for p in pts])
    two_b = K.run_many([("field", K.scan_params("two_t", "beam",
                                                points=[p])[p])
                        for p in pts])
    two_0 = K.run_many([("field", K.switch_params("two_t")["field0_tail"])])
    worst_f = worst_b = 0.0
    for (a, w), sk, sf, sb in zip(pts, kin, two_f, two_b):
        k = sk["emission_vuv"]
        df = 100 * (sf["emission_vuv"] / k - 1)
        db = 100 * (sb["emission_vuv"] / k - 1)
        if (a, w) != (1.0, 300.0):     # v0_5: the corners alone
            worst_f, worst_b = max(worst_f, abs(df)), max(worst_b, abs(db))
        tl = sk.get("tail")
        rr = (f"{tl['field_mean']:.2f} & {tl['field0_mean']:.2f} & "
              f"{tl['beam_mean']:.2f}" if tl else "-- & -- & --")
        lines.append(f"{a:g} & {w:.0f} & {k:.2f} & {sf['emission_vuv']:.2f} "
                     f"({_pm(df)}) & {sb['emission_vuv']:.2f} ({_pm(db)}) & "
                     f"{rr} & {sk['tg_end_K']:.0f}\\\\")
        key = f"{a:g}.{w:.0f}"
        n[f"tail.pt.{key}.kin"] = f"{k:.2f}"
        n[f"tail.pt.{key}.f"] = _pm(df)
        n[f"tail.pt.{key}.b"] = _pm(db)
    n["tail.worst.f"] = f"{worst_f:.0f}"
    n["tail.worst.b"] = f"{worst_b:.0f}"
    k0 = kin[-1]["emission_vuv"]
    n["tail.ref.f0"] = f"{two_0[0]['emission_vuv']:.2f}"
    n["tail.ref.f0.dev"] = _pm(100 * (two_0[0]["emission_vuv"] / k0 - 1))
    n["tail.ref.f.dev"] = _pm(100 * (two_f[-1]["emission_vuv"] / k0 - 1))
    n["tail.ref.b.dev"] = _pm(100 * (two_b[-1]["emission_vuv"] / k0 - 1))
    # v0_5: the corners alone (the reference point is the last)
    corners = [sk for sk in kin[:-1] if sk.get("tail")]
    qf = [sk["tail"]["field_mean"] for sk in corners]
    q0 = [sk["tail"]["field0_mean"] for sk in corners]
    qb = [sk["tail"]["beam_mean"] for sk in corners]
    n["tail.cr.f.lo"], n["tail.cr.f.hi"] = f"{min(qf):.2f}", f"{max(qf):.2f}"
    n["tail.cr.f0.lo"], n["tail.cr.f0.hi"] = f"{min(q0):.2f}", f"{max(q0):.2f}"
    n["tail.cr.b.lo"], n["tail.cr.b.hi"] = f"{min(qb):.2f}", f"{max(qb):.2f}"
    n["tail.ref.kin"] = f"{k0:.2f}"
    n["tail.ref.f"] = f"{two_f[-1]['emission_vuv']:.2f}"
    n["tail.ref.b"] = f"{two_b[-1]['emission_vuv']:.2f}"
    _CLAIMS["tail"] = dict(kin=kin, two_f=two_f, two_b=two_b, two_0=two_0,
                           pts=pts)
    tab = xsec.tail_field_table()
    n["tail.npts"] = str(tab.r.size)
    n["tail.nfloor"] = str(int((tab.r <= 1.0001e-3).sum()))
    n["tail.ntop"] = str(int((tab.r >= 2.9999).sum()))
    n["tail.nxs"] = str(len(tab.xs))
    n["tail.xs.top"] = _sci(float(tab.xs[-1]), 0)
    n["tail.x.lo"] = f"{tab.x[0]:g}"
    n["tail.x.hi"] = f"{tab.x[-1]:g}"
    n["tail.te.lo"] = f"{tab.te[0]:g}"
    n["tail.te.hi"] = f"{tab.te[-1]:g}"
    n["tail.ref.r"] = f"{tab(1.2e-5, 2.5):.3f}"
    n["tail.ref.r3"] = f"{tab(1e-5, 2.5, 1e-3):.3f}"
    n["tail.ref.r0"] = f"{tab(1e-5, 2.5, 0.0):.3f}"
    n["tail.ref.rb"] = f"{xsec.tail_ratio_table()(1.2e-5, 0.0):.3f}"
    # the table's own sensitivity at the states the runs visit: half and
    # twice the pressure, and 700 K gas (tail_field.build at those points),
    # with and without excited atoms
    from v0_5 import tail_field
    pts = [(-5.0, 2.0), (-5.0, 1.414214), (-4.5, 1.681793),
           (-4.5, 2.0), (-4.0, 1.414214), (-5.5, 2.378414)]
    dev = {}
    for key, kw in (("p", dict(p_atm=0.5)), ("p2", dict(p_atm=2.0)),
                    ("tg", dict(tg=700.0))):
        worst = 0.0
        for x, te in pts:
            for xs in (0.0, 1e-4):
                r, held, _ = tail_field.build(x_grid=(x,), te_grid=(te,),
                                              xs_grid=(xs,), **kw)
                if not held[0, 0, 0]:
                    worst = max(worst, abs(r[0, 0, 0]
                                           / tab(10.0 ** x, te, xs) - 1))
        dev[key] = worst
    n["tail.p.dev"] = f"{100 * max(dev['p'], dev['p2']):.0f}"
    n["tail.tg.dev"] = f"{100 * dev['tg']:.0f}"
    return n, ["\\def\\wgTail{" + "\n".join(lines) + "}"]


def el_numbers() -> tuple:
    """Mano et al.'s krypton electroluminescence yield against the swarm,
    each excitation weighted by the photons the network's own cascade makes
    of it at 1.1 bar and 296 K (kr_docs_v0_4.el_params: a dilute seed of
    each 5s level and of the 5p rung, no drive), at the bands' defaults and
    with Kr2* + Kr at its band's top."""
    n, lines = {}, []
    ep = K.el_params()
    got = dict(zip(ep, K.run_many([("field", p) for p in ep.values()])))
    ph = {}
    for (seed, case), s in got.items():
        ev = s["events"]
        vuv = (sum(ev[r] for r in RAD_ROWS)
               + sum(ev[r] for r in ("RES_1s4", "RES_1s2")))
        ph[(seed, case)] = vuv / K.EL_SEED
        name = "5p" if seed == "Kr**" else seed[4:-1]       # 'Kr*(1s5)'
        n[f"el.ph.{name}.{case}"] = f"{vuv / K.EL_SEED:.3f}"
    # the swarm's groups: the 5s levels by Biagi's level cross sections
    # ('1s'), the rest to the 5p rung ('2p', 'hi'); a 5s excitation's
    # photons are its level's, weighted by the swarm's own split
    out = {}
    for case in ("default", "quench"):
        for en in (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 13.5):
            d = swarm.point(en, tg=swarm.MANO_EL_TG, cells=CELLS_SWARM)
            procs = [q for q in swarm.processes(0.0) if q.kind == "excitation"]
            st = swarm.KB.State(x={"Kr": 1.0}, tg=swarm.MANO_EL_TG,
                                n_col=swarm.N_COL,
                                ne=swarm.NE_FRAC * swarm.N_COL)
            s5 = swarm.KB.Solver(procs, n_cells=CELLS_SWARM,
                                 emax=d["emax"], elastic=("Kr",), tol=1e-10,
                                 maxit=4000)
            sol = s5.solve(st, en)
            k5 = {lv: 0.0 for lv in LEVELS_5S}
            for q in procs:
                tag = q.p.product.strip()[3:].rstrip(")").lower()
                if tag in k5:
                    k5[tag] += sol.grid.rate(q.p, sol.F)
            k5t = sum(k5.values())
            p5 = (sum(k5[lv] * ph[(kr5s(lv), case)] for lv in LEVELS_5S)
                  / k5t if k5t else 1.0)
            y = swarm.el_yield(en, photons={"1s": p5,
                                            "2p": ph[("Kr**", case)],
                                            "hi": ph[("Kr**", case)]},
                               cells=CELLS_SWARM)
            out[(case, en)] = y
    for en in (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 13.5):
        y, yq = out[("default", en)], out[("quench", en)]
        lines.append(f"{en:g} & {y['mano']:.3f} & {y['model']:.3f} & "
                     f"{y['ratio']:.2f} & {yq['ratio']:.2f} & "
                     f"{100 * y['share_5s']:.1f}\\\\")
    hi = [out[("default", en)]["ratio"] for en in (6.0, 8.0, 10.0, 12.0, 13.5)]
    n["el.lo"] = f"{min(hi):.2f}"
    n["el.hi"] = f"{max(hi):.2f}"
    n["el.rms"] = f"{100 * _rms(np.array(hi) - 1):.1f}"
    hq = [out[("quench", en)]["ratio"] for en in (6.0, 8.0, 10.0, 12.0, 13.5)]
    n["el.q.lo"] = f"{min(hq):.2f}"
    n["el.q.hi"] = f"{max(hq):.2f}"
    n["el.4"] = f"{out[('default', 4.0)]['ratio']:.2f}"
    means = [swarm.point(en, tg=swarm.MANO_EL_TG, cells=CELLS_SWARM)["mean"]
             for en in (6.0, 13.5)]
    n["el.mean.lo"] = f"{min(means):.1f}"
    n["el.mean.hi"] = f"{max(means):.1f}"
    n["el.5"] = f"{out[('default', 5.0)]['ratio']:.2f}"
    def fit(en):
        return swarm.MANO_EL[0] * en + swarm.MANO_EL[1]
    n["el.frac4"] = f"{fit(4.0) / fit(13.5):.2f}"
    n["el.frac5"] = f"{fit(5.0) / fit(13.5):.2f}"
    n["el.mano.a"] = f"{swarm.MANO_EL[0]:g}"
    n["el.mano.b"] = f"{-swarm.MANO_EL[1]:g}"
    n["el.mano.thr"] = f"{swarm.MANO_EL_RANGE[0]:g}"
    n["el.mano.top"] = f"{swarm.MANO_EL_RANGE[1]:g}"
    if not all(abs(r - 1) < 0.06 for r in hi):
        raise AssertionError("the text says the swarm meets Mano's fit to "
                             "6 % from 6 to 13.5 Td")
    return n, ["\\def\\wgEl{" + "\n".join(lines) + "}"]


# Blagoev, Mishonov & Popov, J. Phys. B 17, 435 (1984), abstract: the
# superelastic rates of thermal electrons on the 5s levels [cm3/s]
BLAGOEV = {"1s5": 1.5e-10, "1s4": 5.2e-10, "1s3": 9.8e-10}
BLAGOEV_PILOSOV = 4.3e-8      # e + 3P2 -> 3P1, Phys. Lett. A 83, 161 (1981)
FLANNERY = 1.0e-15            # cm2, Kr2*'s ionization, the order (1981)


def calc_numbers() -> dict:
    """v0_5's calculated rows against their calculation and their checks:
    Vriens & Smeets's cross section integrated against their rate formula,
    Kannari's printed fits against it at the reference point's Te, its
    peak for Kr2* against Flannery et al.'s order; the superelastic
    returns against Blagoev et al.'s thermal rates (the score)."""
    n = {}
    ref = K.run_many([("field", K.ref_params())])[0]
    te = ref["te_mean_eV"]
    worst = {}
    for tag, e_n in (("1s5", e_form("Kr+") - e_form(kr5s("1s5"))),
                     ("2p10", e_form("Kr+") - e_form("Kr**")),
                     ("kr2", e_form("Kr2+") - e_form("Kr2*"))):
        t = excited.vs_rate_table(e_n)
        worst[tag] = max(abs(float(t(x)) / excited.vs_rate(e_n, x) - 1)
                         for x in (1.0, 2.0, 3.0, 5.0))
        e = np.linspace(e_n, 80.0, 40000)
        sg = excited.vs_sigma(e, e_n)
        n[f"vs.{tag}.En"] = f"{e_n:.2f}"
        n[f"vs.{tag}.peak"] = _sci(float(sg.max()), 1)
        n[f"vs.{tag}.peakE"] = f"{e[np.argmax(sg)]:.1f}"
        n[f"vs.{tag}.k.ref"] = _sci(excited.vs_rate(e_n, te), 1)
    n["vs.int.dev"] = _sci(max(worst.values()), 0)
    # the excimer's ionization time in the pulse, the formula on a
    # Maxwellian at the pulse's mean Te and n_e (an estimate the text
    # calls one), against its radiative lifetime
    n["ei.tau"] = f"{1e9 / (excited.vs_rate(e_form('Kr2+') - e_form('Kr2*'), te) * ref['ne_mean_on']):.0f}"
    # the formula's Kr2* peak against 1s5's (Flannery et al.: 40-50 % above
    # the metastable atom's)
    e = np.linspace(1.0, 80.0, 40000)
    n["vs.kr2.over"] = f"{100 * (float(excited.vs_sigma(e, e_form('Kr2+') - e_form('Kr2*')).max()) / float(excited.vs_sigma(e, e_form('Kr+') - e_form(kr5s('1s5'))).max()) - 1):.0f}"
    n["vs.int.dev05"] = f"{100 * abs(float(excited.vs_rate_table(4.0844)(0.5)) / excited.vs_rate(4.0844, 0.5) - 1):.0f}"
    # Kannari's printed fits (Table I [79], [72]) over the formula at Te,
    # and their range over 1 to 3 eV, all three rows
    kr = []
    for tag, (a, b, c), e_n in (
            ("s", (7.0e-8, 0.68, 4.04), e_form("Kr+") - e_form(kr5s("1s5"))),
            ("p", (2.0e-7, 0.62, 2.53), e_form("Kr+") - e_form("Kr**")),
            ("kr2", (9.0e-8, 0.70, 3.53), e_form("Kr2+") - e_form("Kr2*"))):
        for x, key in ((1.0, "1"), (te, "ref"), (3.0, "3")):
            r = a * x ** b * math.exp(-c / x) / excited.vs_rate(e_n, x)
            n[f"kann.{tag}.{key}"] = f"{r:.2f}"
        kr += [a * x ** b * math.exp(-c / x) / excited.vs_rate(e_n, x)
               for x in np.linspace(1.0, 3.0, 21)]
    n["kann.lo"], n["kann.hi"] = f"{min(kr):.2f}", f"{max(kr):.2f}"
    flan = float(excited.vs_sigma(np.linspace(1, 80, 40000),
                                  e_form("Kr2+") - e_form("Kr2*")).max())
    n["vs.flannery"] = _sci(FLANNERY, 0)
    n["vs.flannery.ratio"] = f"{flan / FLANNERY:.1f}"
    # the superelastic returns against Blagoev et al. (thermal electrons,
    # 300 to 2000 K), the model's own Maxwellian rates
    m = K.Model(K.ref_params("two_t"))
    c = {x.row.rid: x for x in m.compiled}
    kb = 8.617333262e-5
    for lv, v in BLAGOEV.items():
        rr = [m._k(c[f"SUPER_{lv}"], kb * t, 300.0, None, None) / v
              for t in (300.0, 500.0, 1000.0, 2000.0)]
        n[f"blag.{lv}.lo"] = f"{min(rr):.1f}"
        n[f"blag.{lv}.hi"] = f"{max(rr):.1f}"
        n[f"blag.{lv}.meas"] = _sci(v, 1)
        n[f"blag.{lv}.k300"] = _sci(m._k(c[f"SUPER_{lv}"], kb * 300.0,
                                        300.0, None, None), 1)
    n["blag.pil"] = _sci(BLAGOEV_PILOSOV, 1)
    n["blag.pil.top"] = _sci(BANDS["emix_5s"].hi, 1)
    # the 5s mixing by electrons at the reference point's Te, per 1s5 atom,
    # against the excimer's formation out of 1s5
    n_kr = ref["n0"]["Kr"]
    k_up = m._k(c["EMIX_1s5"], te, 300.0, None, None) * ref["ne_mean_on"]
    k3 = BANDS["kr2_form"].default * n_kr ** 2
    n["emix.vs.form"] = f"{k_up / k3:.2f}"
    _CLAIMS["calc"] = dict(worst=worst)
    return n


def claims(n: dict) -> None:
    """What the text says, held against what the parts computed: raises
    when a claim fails (written with the text, 2026-10-16)."""
    c = _CLAIMS

    def need(ok, says):
        if not ok:
            raise AssertionError(f"the text says {says}")
    ref = c["ref_kin"]
    ion = c["ion_kin"]
    # the bottom line
    need(float(n["chg.kin.ratio"]) > 1.0, "the answer falls from v0_4's")
    need(ion["step_s"] + ion["step_p"] + ion["ei"] > 0.5,
         "the stepwise ionizations and the excimer's make most of the ions")
    need(ref["emission_nir"] < 1.0, "the rung's light is under 1 % of the "
         "input (so its lines are not propagated)")
    need(0.5 < float(n["vs.flannery.ratio"]) < 2.0,
         "the formula gives Kr2* Flannery et al.'s 1e-15 cm2")
    # the calculated bands, 0 to twice the calculation; the electrons' rows
    # on Kr2* at 0, to Kannari's 1e-7
    for name in CALCULATED:
        b = BANDS[name]
        need((b.lo, b.default, b.hi) == (0.0, 1.0, 2.0),
             f"{name} runs from 0 to twice the calculation")
    for name in ("eq_kr2", "em_kr2"):
        b = BANDS[name]
        need((b.lo, b.default, b.hi) == (0.0, 0.0, 1.0e-7),
             f"{name} is 0, to Kannari's 1e-7")
    # the electrons' 1s5 -> 1s4: the band's top is the measurement taken
    # at 0.07 eV
    b = BANDS["emix_5s"]
    gap = e_form(kr5s("1s4")) - e_form(kr5s("1s5"))
    need(b.default == BLAGOEV_PILOSOV
         and abs(b.hi / (BLAGOEV_PILOSOV * math.exp(gap / 0.07)) - 1) < 0.02,
         "emix_5s runs from Blagoev and Pilosov's value to it read at "
         "0.07 eV")
    # the rung at 300 K is nearly all 2p10
    need(excited.rung_fraction("2p10", 300.0) > 0.95,
         "the rung's 300-K lifetime is nearly 2p10's alone")
    # the hot corner lies outside the formation's measured range
    need(float(n["scan.tg01.lo"]) > KR2_TG_VALID_K,
         "the 0.1-atm column ends above the formation's measured range")
    # what dominates the uncertainty: the electrons on Kr2*, the ionization
    # first; the quenching cuts, the mixing raises (the fast state radiating
    # before the electrons ionize it); the J = 1 proxies now the other way
    order, res, y0 = c["unc_order"], c["unc_res"], c["unc_ref"]
    need(order[0] == "ei_kr2" and set(order[:3]) == {"ei_kr2", "em_kr2",
                                                     "eq_kr2"},
         f"the three largest terms are the electrons on Kr2*: {order[:4]}")

    def end(name, e):
        return res.get((name, e), c["ref_kin"])["emission_vuv"]
    need(end("eq_kr2", "high") < y0 < end("em_kr2", "high"),
         "the quenching at its top cuts the share and the mixing raises it")
    rows = {r.rid: r for r in build()}
    need(1e9 / rows["KR2F_RAD"].rate[1] < float(n["ei.tau"]),
         "the fast state radiates before the electrons can ionize it")
    need(end("k_5s5p_j1", "low") < y0,
         "with no J = 1 pumping the share falls (v0_4: it rose)")
    # the changes undone: the calculated rows at 0 lie below v0_4; the J = 1
    # pumping raises the share in v0_5 and lowered it in v0_4
    u, uref = c["undo"]["kin"], c["undo_ref"]["kin"]["emission_vuv"]
    need(u["no_calc"]["emission_vuv"] < float(n["chg.kin.v04"]),
         "with the calculated rows at 0 the share is below v0_4's")
    need(u["no_j1"]["emission_vuv"] < uref, "the J = 1 pumping raises it")
    # what lies between the calculated rows at 0 and v0_4 is chiefly what
    # cannot be undone alone: the changes that can (the superelastic
    # return, the electrons' 1s5 -> 1s4, the rung's reverse) move the share
    # by far less than that gap
    small = sum(abs(uref / u[k]["emission_vuv"] - 1)
                for k in ("no_super", "no_emix", "no_kpsr"))
    gap = float(n["chg.kin.v04"]) / u["no_calc"]["emission_vuv"] - 1
    need(small < gap / 3, "the rest of v0_5's change is chiefly what cannot "
         "be undone alone")
    j4 = K4.run_many([("field", K4.undo_params()["no_j1"])])[0]
    r4 = K4.run_many([("field", K4.ref_params())])[0]
    need(j4["emission_vuv"] > r4["emission_vuv"],
         "in v0_4 the J = 1 pumping lowered the share")
    # the reference point: the excimers mostly do not radiate
    need(ref["to_Kr2star"] > 0 and float(n["ref.kin.k2rad"]) < 50.0
         and float(n["ref.kin.k2ei"]) > float(n["ref.kin.k2rad"]),
         "the excimers mostly do not radiate; the electrons ionize more")
    # the scan: rises with pressure, falls with power, nothing at 0.1 atm
    sc = c["scan"]
    for w in K.SCAN_KW:
        v = [sc[(a, w)]["emission_vuv"] for a in K.SCAN_ATM]
        need(all(x < y for x, y in zip(v, v[1:])),
             f"the share rises with pressure ({w:.0f} kW)")
    for a in K.SCAN_ATM:
        v = [sc[(a, w)]["emission_vuv"] for w in K.SCAN_KW]
        need(all(x > y for x, y in zip(v, v[1:])),
             f"the share falls with power ({a:g} atm)")
    need(max(sc[(0.1, w)]["emission_vuv"] for w in K.SCAN_KW) < 0.05,
         "nothing at 0.1 atm")
    # the closure: its tail follows the kinetic one along the reference
    # run, the zero block's is hotter, and the zero block's share overshoots
    tl = ref["tail"]
    t = c["tail"]
    k0 = t["kin"][-1]["emission_vuv"]
    need(0.9 < tl["field_lo"] and tl["field_hi"] < 1.1
         and tl["field0_lo"] > 1.1,
         "the table's tail follows the kinetic one; the zero block's is "
         "hotter")
    need(t["two_f"][-1]["emission_vuv"] < k0 < t["two_0"][0]["emission_vuv"],
         "the bulk's error is below, the zero block's tail more than "
         "offsets it")
    # the beam: low pump as v0_4, no excimer ionized there; at 0.5 atm and
    # high pump the new rows act
    need(c["moz_same"], "the low-pump efficiency is v0_4's")
    need(c["moz_ei"] < 0.01 * c["moz_ph"],
         "nearly every excimer radiates at low pump")
    need(c["eck_vs_v04_0.5"] < -0.05,
         "the new rows act at 0.5 atm under the high pump")


def numbers() -> tuple:
    # every run at once, in parallel, so a stale cache is one batch
    K.run_many(K.all_jobs())
    n = {"version": K.V, "md5": K.full_md5()[:12]}
    tables = []
    for part in (reference_numbers, uncertainty_numbers, swarm_numbers,
                 undo_numbers, tail_numbers, el_numbers):
        a, t = part()
        n.update(a)
        tables += t
    for part in (beam_numbers, convergence_numbers, change_numbers,
                 jung_numbers, ladder_numbers, calc_numbers):
        n.update(part())
    n["tau5p"] = f"{KR_5P_TAU_S * 1e9:.1f}"
    claims(n)
    return n, tables


def render(n: dict, tables=()) -> str:
    lines = [
        "% LIGHT_Kr v0_5 overview numbers. GENERATED by make_overview_v0_5.py",
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
    ap = argparse.ArgumentParser(description="LIGHT_Kr v0_5's accuracy "
                                 "overview numbers")
    ap.add_argument("--check", action="store_true",
                    help="compare the file with a fresh computation, write "
                    "nothing")
    a = ap.parse_args(argv)
    fresh = render(*numbers())
    if a.check:
        if not OUT.exists() or OUT.read_text() != fresh:
            print(f"{OUT.name} is stale: rerun make_overview_v0_5.py")
            return 1
        print(f"{OUT.name} is current")
        return 0
    OUT.write_text(fresh)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
