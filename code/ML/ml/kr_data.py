"""Krypton's electron swarm measurements: the data that train and score
LIGHT_Kr's cross sections (ml/kr_swarm.py).

Swarm data measure the cross sections themselves, so under the
predictive-only policy they may train them (ML/README.md, "The rule for
data"); laser and lamp data never do. Every dataset here has a role:
  train     fits the cross sections (the consistency check's group)
  held_out  scores the fit and never enters it
  check     outside what the forward model can claim (above 300 Td the
            two-term approximation and the growth model are not the tool;
            Specht's 'total' coefficient includes the surface photoemission
            his 'direct' one removes)
The roles of the candidates are set by the leave-one-out consistency check
(kr_swarm.consistency, `python3 ML/run.py kr-consistency`, 2026-10-14);
see ROLES below and ML/README.md.

Where the numbers come from:
  Dutton's survey (J. Phys. Chem. Ref. Data 4, 577 (1975)), its printed
    tables, transcribed in ML/data/swarm/*_via_Dutton1975.csv: Pack 1962's
    drift velocity, Kruithof 1940's and Heylen 1971's alpha/N;
  LXCat's LAPLACE database (measurements after 1975), the owner's
    approval to accept LXCat's terms and read it (2026-10-14): the four
    Kr_*_LAPLACE_LXCat.txt files, LXCat's own text output, unedited and
    md5-pinned below. LAPLACE digitized them from the papers' tables
    (Hunter 1988, Jacques 1986, Pack 1992) or figures (the rest).
Errors: none of these files carries one. A dataset's sigma is a class
value, stated as assumed in its note: drift velocity 3 % (after 1985) or
5 % (before), characteristic energy 7 %, ionization coefficient 10 %.
An E/N uncertainty of 2 % enters through the model's own slope
(kr_swarm.effective_sigma), which matters where W climbs steeply out of
the Ramsauer minimum (0.01 to 0.03 Td).
"""
from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data" / "swarm"
TD_M2 = 1e-21                     # 1 Td in V m^2

# LXCat's text output, as served (CRLF), pinned
LAPLACE = {
    "alphaN": ("Kr_alphaN_LAPLACE_LXCat.txt", "ff721f19cf4aaa982044951b9ed880cf"),
    "muN": ("Kr_muN_LAPLACE_LXCat.txt", "4735b6c834e3fd99b07f8c198aa09194"),
    "energy": ("Kr_char_energy_LAPLACE_LXCat.txt",
               "b8be54610dcd4430784ae3f125c56812"),
    "NDL": ("Kr_NDL_LAPLACE_LXCat.txt", "32a8add103c49ac2c755813ee6e920f3"),
}
DUTTON_W = "Kr_W_Pack1962_via_Dutton1975.csv"
DUTTON_ALPHA = "Kr_alphaN_Kruithof1940_Heylen1971_via_Dutton1975.csv"

# the forward model's reach for alpha/N (two-term, temporal growth)
ALPHA_TOP_TD = 300.0
# below this the low-field alpha/N is its own dataset (Bordage et al. 2013:
# impurity Penning ionization; LIGHT_Kr's Known-bad list)
ALPHA_LOW_TD = 50.0
# no point above this is loaded: the mean energy passes 13 eV and the
# solver's grid runs past the end of Biagi's tables (LIGHT_Kr's overview
# stopped at 1017 Td for the same reason)
MODEL_TOP_TD = 1100.0

SOURCES = {
    "Hunter 1988": "S. R. Hunter, J. G. Carter and L. G. Christophorou, "
                   "Phys. Rev. A 38, 5539 (1988), doi:10.1103/PhysRevA.38.5539",
    "Pack 1992": "J. L. Pack, R. E. Voshall, A. V. Phelps and L. E. Kline, "
                 "J. Appl. Phys. 71, 5363 (1992), doi:10.1063/1.350555",
    "Pack 1962": "J. L. Pack, R. E. Voshall and A. V. Phelps, Phys. Rev. "
                 "127, 2084 (1962), doi:10.1103/PhysRev.127.2084",
    "Nakamura 1990": "Y. Nakamura, in Nonequilibrium Effects in Ion and "
                     "Electron Transport (Plenum, 1990), p. 363, "
                     "doi:10.1007/978-1-4613-0661-0_33",
    "Dall'Armi 1992": "G. Dall'Armi, K. L. Brown, P. H. Purdie and J. "
                      "Fletcher, Aust. J. Phys. 45, 185 (1992), "
                      "doi:10.1071/PH920185",
    "Bowe 1960": "J. C. Bowe, Phys. Rev. 117, 1411 (1960), "
                 "doi:10.1103/PhysRev.117.1411",
    "English 1953": "W. N. English and G. C. Hanna, Can. J. Phys. 31, 768 "
                    "(1953), doi:10.1139/p53-070",
    "Kucukarpaci 1981": "H. N. Kucukarpaci and J. Lucas, J. Phys. D 14, 2001 "
                        "(1981), doi:10.1088/0022-3727/14/11/008",
    "Koizumi 1986": "T. Koizumi, E. Shirakawa and I. Ogawa, J. Phys. B 19, "
                    "2331 (1986), doi:10.1088/0022-3700/19/15/014",
    "Al-Amin 1987": "S. A. J. Al-Amin and J. Lucas, J. Phys. D 20, 1590 "
                    "(1987), doi:10.1088/0022-3727/20/12/007",
    "Jacques 1986": "L. Jacques, W. Bruynooghe, R. Boucique and W. Wieme, "
                    "J. Phys. D 19, 1731 (1986), "
                    "doi:10.1088/0022-3727/19/9/017",
    "Bhattacharya 1979": "A. K. Bhattacharya, J. Appl. Phys. 50, 6207 "
                         "(1979), doi:10.1063/1.325754",
    "Specht 1980": "L. T. Specht, S. A. Lawton and T. A. DeTemple, J. Appl. "
                   "Phys. 51, 166 (1980), doi:10.1063/1.327395",
    "Kruithof 1940": "A. A. Kruithof, Physica 7, 519 (1940), "
                     "doi:10.1016/S0031-8914(40)90043-X",
    "Heylen 1971": "A. E. D. Heylen, Int. J. Electron. 31, 19 (1971), "
                   "doi:10.1080/00207217108938196",
}
VIA_DUTTON = ("J. Dutton, J. Phys. Chem. Ref. Data 4, 577 (1975), "
              "doi:10.1063/1.555525")
VIA_LXCAT = ("LXCat's LAPLACE database (L. C. Pitchford et al., Plasma "
             "Process. Polym. 14, 1600098 (2017), doi:10.1002/ppap.201600098)")


@dataclass
class KrSwarm:
    name: str
    kind: str              # W [cm/s] | alphaN [cm^2] | DTmu [V]
    tg: float              # gas temperature [K]
    x: np.ndarray          # E/N [Td]
    y: np.ndarray
    sigma: float           # fractional 1-sigma error (assumed; see note)
    group: str             # the experiment (SOURCES key)
    via: str
    note: str = ""
    role: str = ""         # set from ROLES
    cache: dict = field(default_factory=dict, repr=False)

    @property
    def source(self) -> str:
        return SOURCES[self.group]


def md5_of(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def _pinned(key) -> Path:
    name, md5 = LAPLACE[key]
    path = DATA / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}")
    return path


def lxcat_blocks(path) -> list:
    """The e / Kr blocks of an LXCat swarm text file: process, comment,
    gas temperature (PARAM 'Tgas = ...', else None) and the table."""
    txt = Path(path).read_text(encoding="utf-8").replace("\r\n", "\n")
    out = []
    for blk in txt.split("SPECIES: e / Kr")[1:]:
        proc = re.search(r"PROCESS: (.*)", blk).group(1).strip()
        com = re.search(r"COMMENT: (.*)", blk).group(1).strip()
        par = re.search(r"PARAM\.:\s*Tgas = (\d+)", blk)
        tab = blk.split("-----------------------------")[1]
        xy = np.array([[float(v) for v in ln.split()]
                       for ln in tab.strip().splitlines()])
        out.append({"process": proc, "comment": com,
                    "tg": float(par.group(1)) if par else None,
                    "x": xy[:, 0], "y": xy[:, 1]})
    return out


def _csv(name):
    with open(DATA / name) as fh:
        return list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))


def _cut(name, d: KrSwarm, lo, hi, note) -> KrSwarm:
    m = (d.x >= lo) & (d.x < hi)
    return KrSwarm(name, d.kind, d.tg, d.x[m], d.y[m], d.sigma, d.group,
                   d.via, note)


def load() -> list:
    """Every dataset, its role from ROLES."""
    out = []
    # --- drift velocity
    rows = _csv(DUTTON_W)
    for tg in ("300", "368", "195"):
        q = [r for r in rows if r["T_gas_K"] == tg]
        out.append(KrSwarm(
            f"w_pack62_{tg}", "W", float(tg),
            np.array([float(r["E_N_Td"]) for r in q]),
            np.array([float(r["W_cm_per_s"]) for r in q]), 0.05, "Pack 1962",
            VIA_DUTTON, "printed; no error printed, 5 % assumed. Two W "
            "values repeat in Dutton's table"))
    muN = {
        "Hunter et al 1988.": ("w_hunter88", "Hunter 1988", 301.0, 0.03,
                               "LAPLACE from the paper's table; 301 K"),
        "Pack et al 1992.": ("w_pack92", "Pack 1992", None, 0.03,
                             "LAPLACE from the paper's table"),
        "Nakamura 1990.": ("w_nakamura90", "Nakamura 1990", 300.0, 0.03,
                           "LAPLACE digitized from a figure; room "
                           "temperature taken as 300 K"),
        "DallArmi et al 1992.": ("w_dallarmi92", "Dall'Armi 1992", 300.0,
                                 0.03, "LAPLACE digitized from a figure"),
        "Bowe 1960.": ("w_bowe60", "Bowe 1960", 300.0, 0.05,
                       "LAPLACE digitized from a figure"),
        "English and Hanna 1953.": ("w_english53", "English 1953", 300.0,
                                    0.05, "LAPLACE digitized from a "
                                    "figure"),
        "kucukarpaci et al 1981.": ("w_kucuk81", "Kucukarpaci 1981", 300.0,
                                    0.05, "LAPLACE digitized from a "
                                    "figure"),
    }
    for b in lxcat_blocks(_pinned("muN")):
        nm, grp, tg, sg, note = muN[b["comment"]]
        tg = b["tg"] if tg is None else tg
        if b["tg"] == 195.0:
            nm += "_195"
        out.append(KrSwarm(nm, "W", tg, b["x"], b["y"] * b["x"] * TD_M2
                           * 100.0, sg, grp, VIA_LXCAT,
                           note + f"; mobility x N times E/N; {sg:.0%} "
                           "assumed"))
    # --- characteristic energy D_T/mu
    en = {"Koizumi et al 1986 (DT/mu).": ("dt_koizumi86", "Koizumi 1986"),
          "Al-Amin et al 1987 (DT/mu).": ("dt_alamin87", "Al-Amin 1987")}
    for b in lxcat_blocks(_pinned("energy")):
        if b["comment"] not in en:
            continue                      # D_L/mu: no longitudinal model
        nm, grp = en[b["comment"]]
        out.append(KrSwarm(nm, "DTmu", 300.0, b["x"], b["y"], 0.07, grp,
                           VIA_LXCAT, "LAPLACE digitized from a figure; "
                           "7 % assumed"))
    # --- ionization coefficient
    al = {"(Direct) Specht et al 1980.": ("a_specht80_direct", "Specht 1980"),
          "(Total) Specht et al 1980.": ("a_specht80_total", "Specht 1980"),
          "Bhattacharya 1979.": ("a_bhatta79", "Bhattacharya 1979"),
          "DallArmi et al 1992.": ("a_dallarmi92", "Dall'Armi 1992"),
          "Jacques et al 1986.": ("a_jacques86", "Jacques 1986")}
    sets = []
    for b in lxcat_blocks(_pinned("alphaN")):
        nm, grp = al[b["comment"]]
        sets.append(KrSwarm(nm, "alphaN", 300.0, b["x"], b["y"] * 1e4, 0.10,
                            grp, VIA_LXCAT, "LAPLACE "
                            + ("from the paper's table" if "Jacques" in grp
                               else "digitized from a figure")
                            + "; 10 % assumed"))
    rows = _csv(DUTTON_ALPHA)
    for src, nm, grp in (("Kruithof1940", "a_kruithof40", "Kruithof 1940"),
                         ("Heylen1971", "a_heylen71", "Heylen 1971")):
        q = [r for r in rows if r["source"] == src]
        sets.append(KrSwarm(
            nm, "alphaN", 300.0, np.array([float(r["E_N_Td"]) for r in q]),
            np.array([float(r["alpha_over_N_cm2"]) for r in q]), 0.10, grp,
            VIA_DUTTON, "printed; no error printed, 10 % assumed"))
    for d in sets:
        if d.name in ("a_specht80_total", "a_specht80_direct"):
            out.append(d)
            continue
        if d.x.min() < ALPHA_LOW_TD:
            out.append(_cut(d.name + "_low", d, 0.0, ALPHA_LOW_TD,
                            d.note + f"; below {ALPHA_LOW_TD:g} Td"))
        out.append(_cut(d.name, d, ALPHA_LOW_TD, 1e9,
                        d.note + f"; from {ALPHA_LOW_TD:g} Td"))
    # every kind: the points above ALPHA_TOP_TD are their own dataset (check)
    split = []
    for d in out:
        if d.x.max() > ALPHA_TOP_TD:
            split.append(_cut(d.name, d, 0.0, ALPHA_TOP_TD + 1e-9, d.note))
            split.append(_cut(d.name + "_high", d, ALPHA_TOP_TD + 1e-9, 1e9,
                              d.note + f"; above {ALPHA_TOP_TD:g} Td"))
        else:
            split.append(d)
    out = []
    for d in split:
        keep = d.x <= MODEL_TOP_TD
        d.x, d.y = d.x[keep], d.y[keep]
        if len(d.x):
            out.append(d)
    for d in out:
        d.role = ROLES.get(d.name, "held_out")
    return out


# The roles. 'check' is set by what the forward model can claim: above
# ALPHA_TOP_TD, Specht's total coefficient, and D_T/mu above 1 Td, where
# LoKI-MC's Monte Carlo puts the two-term transverse diffusion far off
# (+70 % at 50 Td on Biagi's set, 2026-10-14). Below ALPHA_LOW_TD the
# ionization coefficient never trains (LIGHT_Kr's Known-bad list:
# Bordage et al. 2013 ascribe the measurements' excess there to impurity
# Penning ionization, and the measurements part by x2 to x30): those sets
# are held out by that rule. The train set is the consistency check's
# among the rest (kr_swarm.consistency): filled in from its result,
# every other dataset held out.
CHECK = ("a_specht80_total", "a_kruithof40_high", "a_bhatta79_high",
         "w_kucuk81_high", "dt_alamin87", "dt_alamin87_high")
HELD_BY_RULE = ("a_specht80_direct", "a_kruithof40_low", "a_heylen71_low",
                "a_dallarmi92_low", "a_jacques86_low")
# the consistency check's train set (kr_train.run_consistency, threshold
# 3, 2026-10-14; the same at threshold 2)
TRAIN = ("w_pack62_368", "w_hunter88", "w_nakamura90", "a_bhatta79",
         "a_dallarmi92", "a_jacques86", "a_kruithof40", "a_heylen71")
ROLES = {**{n: "check" for n in CHECK},
         **{n: "held_out" for n in HELD_BY_RULE},
         **{n: "train" for n in TRAIN}}


def by_name(*names) -> list:
    d = {s.name: s for s in load()}
    return [d[n] for n in names]


def by_role(role) -> list:
    return [d for d in load() if d.role == role]
