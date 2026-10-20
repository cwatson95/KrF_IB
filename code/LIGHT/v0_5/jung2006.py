"""The measured electron-impact excitation cross sections out of the
metastable levels of krypton into the ten 4p5 5p levels, and the lumped
Kr(5s) -> Kr(5p) cross section LIGHT v0_5 builds from them.

Source: Jung, Stone, Boffard, Anderson and Lin, Phys. Rev. A 73, 022722
(2006), Tables I (level energies), II (apparent cross sections out of the
J = 2 metastable 1s5, +-35 %, 2p3 and 2p4 +-50 %), III (out of the J = 0
metastable 1s3, +-40 %) and IV (the 3.5-eV and 6-eV values including the
core-changing 1s3 -> 2p5), read from the paper in LIGHT's new_articles/.
The same group's Letter, Phys. Rev. Lett. 94, 163202 (2005), Table I,
prints the 3.5-eV column. Apparent cross sections include cascades (under
10 %, their Sec. II F); they are used as printed.

What the model needs is one cross section for its lumped Kr* (the Biagi
'Kr(s)' manifold, the four 4p5 5s levels at 9.915 eV) into its lumped Kr**
(the 4p5 5p manifold at 11.30 eV). Jung et al. measured out of the two
metastables, 1s5 (J = 2, statistical weight 5) and 1s3 (J = 0, weight 1);
nothing is measured out of the two J = 1 resonance levels 1s4 and 1s2
(weight 3 each), which at 1 atm are radiation-trapped and populated with
the metastables. The lumped cross section is the statistical-weight mean
over the manifold, each level's cross section summed over the ten 5p
levels, with the J = 1 levels given the measured set of the metastable of
the same ion core (1s4 the 1s5 set, 1s2 the 1s3 set; the thresholds
shifted to their own energies). That proxy is the one assumption here and
the band k_5s5p runs down to the metastables alone (the J = 1 levels
contributing nothing, half the manifold's weight) and up to the
measurement's +40 %.

Energy dependence: linear between the tabulated energies, linear from 0
at each pair's own threshold (E_f - E_i, Table I) to the first tabulated
point, and above the last tabulated point each column follows the shape
of the 1s5 -> 2p9 column, which is tabulated to 300 eV (the dipole-allowed
fall-off), then ln(E)/E. The energy an electron loses in the lumped
collision is the ledger's gap E(Kr**) - E(Kr*) = 1.385 eV (the Process
threshold), as in v0_4, so the ledger closes; the per-level thresholds
(1.39 to 1.76 eV) only shape the onset.

No number here is fitted to anything. k_5s5p is a scale on this cross
section, 1 by default.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import xsec

# Table I: level energies [eV] above the ground state, Paschen notation
LEVELS = {
    "2p1": 12.26, "2p2": 12.14, "2p3": 12.14, "2p4": 12.10, "2p5": 11.67,
    "2p6": 11.55, "2p7": 11.53, "2p8": 11.44, "2p9": 11.44, "2p10": 11.30,
    "1s2": 10.64, "1s3": 10.56, "1s4": 10.03, "1s5": 9.91,
}
# statistical weights 2J + 1 of the 5s levels (J: 1s5 2, 1s4 1, 1s3 0, 1s2 1)
WEIGHTS_5S = {"1s5": 5, "1s4": 3, "1s3": 1, "1s2": 3}
# the J = 1 levels take the measured set of the metastable with the same core
PROXY = {"1s4": "1s5", "1s2": "1s3"}
UNCERTAINTY = {"1s5": 0.35, "1s3": 0.40}       # the papers' totals

# Table II: out of 1s5 (J = 2), 10^-16 cm2, energy [eV] -> {level: sigma}
TABLE_II = {
    2: {"2p2": 0.0, "2p3": 0.0, "2p4": 0.0, "2p5": 0.32, "2p6": 8.6,
        "2p7": 4.7, "2p8": 4.8, "2p9": 27.0, "2p10": 9.8},
    3: {"2p2": 0.47, "2p3": 0.5, "2p4": 0.5, "2p5": 0.76, "2p6": 18.0,
        "2p7": 3.8, "2p8": 7.6, "2p9": 47.0, "2p10": 15.0},
    4: {"2p2": 0.66, "2p3": 0.4, "2p4": 0.3, "2p5": 0.34, "2p6": 21.0,
        "2p7": 2.8, "2p8": 10.0, "2p9": 52.0, "2p10": 15.0},
    6: {"2p2": 0.30, "2p3": 0.2, "2p4": 0.1, "2p5": 0.11, "2p6": 25.0,
        "2p7": 2.6, "2p8": 9.6, "2p9": 52.0, "2p10": 17.0},
    8: {"2p2": 0.10, "2p3": 0.2, "2p4": 0.05, "2p5": 0.04, "2p6": 24.0,
        "2p7": 2.2, "2p8": 8.6, "2p9": 51.0, "2p10": 17.0},
    10: {"2p6": 22.0, "2p8": 8.2, "2p9": 49.0},
    15: {"2p6": 17.0, "2p9": 41.0},
    25: {"2p6": 14.0, "2p9": 37.0},
    50: {"2p6": 8.8, "2p9": 23.0},
    100: {"2p6": 4.8, "2p9": 16.0},
    200: {"2p6": 3.1, "2p9": 8.9},
    300: {"2p6": 2.5, "2p9": 6.0},
}
# Table III: out of 1s3 (J = 0), with the core-changing 1s3 -> 2p5 of
# Table IV (0.4 at 3.5 eV, 0.11 at 6 eV)
TABLE_III = {
    2: {"2p1": 0.90, "2p2": 7.0, "2p3": 17.0, "2p4": 18.0},
    3: {"2p1": 1.3, "2p2": 18.0, "2p3": 43.0, "2p4": 48.0},
    3.5: {"2p5": 0.4},
    4: {"2p1": 1.1, "2p2": 14.0, "2p3": 55.0, "2p4": 65.0},
    6: {"2p1": 0.63, "2p2": 5.5, "2p3": 59.0, "2p4": 62.0, "2p5": 0.11},
    8: {"2p1": 0.23, "2p2": 3.6, "2p3": 60.0, "2p4": 63.0},
    10: {"2p3": 56.0, "2p4": 67.0},
}
TABLES = {"1s5": TABLE_II, "1s3": TABLE_III}
SHAPE_REF = ("1s5", "2p9")                     # the column followed above the tables
E_MAX_TABLE = 300.0
CM2 = 1e-16


def _column(initial: str, final: str):
    """(energies, cross sections [cm2]) of one measured pair, ascending."""
    pts = sorted((e, row[final]) for e, row in TABLES[initial].items()
                 if final in row)
    if not pts:
        return None
    e = np.array([p[0] for p in pts], dtype=float)
    s = np.array([p[1] for p in pts], dtype=float) * CM2
    return e, s


@lru_cache(maxsize=None)
def _shape_ref():
    e, s = _column(*SHAPE_REF)
    return e, s


def sigma_pair(initial: str, final: str, eps) -> np.ndarray:
    """sigma(eps) [cm2] of one measured pair at the energies eps [eV]: 0
    below the pair's threshold E_f - E_i, linear to the first tabulated
    point, linear between points, the 1s5 -> 2p9 shape above the last
    point (ln E / E above 300 eV)."""
    scalar = np.ndim(eps) == 0
    eps = np.atleast_1d(np.asarray(eps, dtype=float))
    col = _column(initial, final)
    if col is None:
        return np.zeros_like(eps) if not scalar else 0.0
    e_tab, s_tab = col
    thr = LEVELS[final] - LEVELS[initial]
    e_pts = np.concatenate([[thr], e_tab])
    s_pts = np.concatenate([[0.0], s_tab])
    out = np.interp(eps, e_pts, s_pts, left=0.0, right=np.nan)
    above = eps > e_tab[-1]
    if above.any():
        e_ref, s_ref = _shape_ref()
        ref_last = np.interp(e_tab[-1], e_ref, s_ref)
        ea = eps[above]
        ref = np.interp(ea, e_ref, s_ref, right=s_ref[-1])
        tail = s_ref[-1] * (np.log(ea) / ea) / (math.log(E_MAX_TABLE) / E_MAX_TABLE)
        ref = np.where(ea > E_MAX_TABLE, tail, ref)
        out[above] = s_tab[-1] * ref / ref_last
    out = np.where(eps < thr, 0.0, out)
    return float(out[0]) if scalar else out


def level_total(initial: str, eps, as_level: str = None) -> np.ndarray:
    """The sum over the ten 5p levels of the cross section out of one 5s
    level [cm2]. as_level: book the thresholds at another 5s level's
    energy (the J = 1 proxies)."""
    eps = np.asarray(eps, dtype=float)
    shift = 0.0 if as_level is None else LEVELS[as_level] - LEVELS[initial]
    # a higher-lying initial level sees the same final levels at a
    # smaller energy defect: shift the energy axis by the level spacing
    tot = np.zeros_like(eps, dtype=float)
    for final in [k for k in LEVELS if k.startswith("2p")]:
        tot += sigma_pair(initial, final, eps + shift)
    return tot


def lumped_sigma(eps, proxy_j1: bool = True) -> np.ndarray:
    """The statistical-weight mean over the 5s manifold of each level's
    total cross section into the 5p manifold [cm2]."""
    eps = np.asarray(eps, dtype=float)
    tot = np.zeros_like(eps, dtype=float)
    wsum = sum(WEIGHTS_5S.values())
    for lvl, w in WEIGHTS_5S.items():
        if lvl in TABLES:
            tot += w / wsum * level_total(lvl, eps)
        elif proxy_j1:
            tot += w / wsum * level_total(PROXY[lvl], eps, as_level=lvl)
    return tot


def energy_grid() -> np.ndarray:
    thr = min(LEVELS[f] - LEVELS[i] for i in WEIGHTS_5S for f in LEVELS
              if f.startswith("2p"))
    return np.concatenate([[0.0, 0.999 * thr], np.geomspace(thr, 1000.0, 600)])


@lru_cache(maxsize=None)
def process(scale: float = 1.0, proxy_j1: bool = True, threshold: float = None):
    """The lumped Kr* -> Kr** excitation as an xsec.Process for the kinetic
    solver. threshold: the energy the electron loses (the ledger's gap by
    default, E(Kr**) - E(Kr*)); the cross section's onset is the measured
    one regardless."""
    from .species import e_form
    gap = (e_form("Kr**") - e_form("Kr*")) if threshold is None else threshold
    e = energy_grid()
    return xsec.Process("excitation", "Kr*", "Kr**(5p)", gap, e,
                        scale * lumped_sigma(e, proxy_j1))


@lru_cache(maxsize=None)
def rate_table(scale: float = 1.0, proxy_j1: bool = True) -> xsec.RateTable:
    """k(Te) = <sigma v> over a Maxwellian [cm3/s]."""
    p = process(scale, proxy_j1)
    k = xsec._mean(p.sigma_at(xsec.EPS) * xsec.V_CM_S * np.sqrt(xsec.EPS))
    return xsec.RateTable(k, f"Kr* -> Kr** (Jung 2006, x{scale:g})")


@lru_cache(maxsize=None)
def superelastic_table(g_ratio: float, scale: float = 1.0,
                       proxy_j1: bool = True) -> xsec.RateTable:
    """The reverse Kr** + e -> Kr* + e by the Klein-Rosseland relation on
    the lumped cross section at the ledger's gap, integrated directly (as
    xsec.superelastic_table does for the LXCat rows)."""
    p = process(scale, proxy_j1)
    e_thr = p.threshold
    sig = g_ratio * (xsec.EPS + e_thr) / xsec.EPS * p.sigma_at(xsec.EPS + e_thr)
    k = xsec._mean(sig * xsec.V_CM_S * np.sqrt(xsec.EPS))
    return xsec.RateTable(k, f"Kr** -> Kr* superelastic (Jung 2006, x{scale:g})")


@lru_cache(maxsize=None)
def capture_table(scale: float = 1.0, proxy_j1: bool = True):
    """<eps sigma v> / <sigma v> [eV]: the mean energy of the electrons the
    forward collision takes from a Maxwellian at Te (they lose the gap and
    keep the rest; model.py books the gap through the row's routing and
    the Boltzmann solver through the Process threshold)."""
    p = process(scale, proxy_j1)
    sv = p.sigma_at(xsec.EPS) * xsec.V_CM_S * np.sqrt(xsec.EPS)
    k = xsec._mean(sv)
    ke = xsec._mean(sv * xsec.EPS)
    ok = k > xsec._K_FLOOR
    e = np.where(ok, ke / np.where(ok, k, 1.0), 1.5 * xsec.TE)
    return xsec.EnergyTable(e, "Kr* -> Kr** (Jung 2006): captured energy")


def kannari_rate(te, a=9.0e-7, b=0.5, c=1.51):
    """Kannari's Table I [76] fit (Jacob and Mangano), for comparison."""
    te = np.asarray(te, dtype=float)
    return a * te ** b * np.exp(-c / te)


def summary(te=(1.0, 1.16, 2.0, 3.0, 5.8)) -> str:
    """The lumped rate against Kannari's fit at a few Te, for the README."""
    k = rate_table()
    lines = [f"{'Te':>6} {'Jung k':>10} {'Kannari k':>10} {'ratio':>6}"]
    for t in te:
        kj, kk = float(k(t)), float(kannari_rate(t))
        lines.append(f"{t:6.2f} {kj:10.3e} {kk:10.3e} {kj / kk:6.2f}")
    e6 = lumped_sigma(6.0)
    lines.append(f"lumped sigma at 6 eV: {e6:.3e} cm2 (metastables only: "
                 f"{lumped_sigma(6.0, proxy_j1=False):.3e})")
    return "\n".join(lines)
