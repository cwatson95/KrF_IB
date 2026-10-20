"""The electron-beam drive: REBUILD's beam physics on LIGHT's own data.

LIGHT_ArF_ebeam (2026-09-30). LIGHT's own input is a given power into the
electron gas. An electron beam does not heat the electron gas: its fast
electrons make ions and excited atoms directly and leave each secondary
electron a few eV. This module is that drive, adapted from
SPECK/REBUILD/MODULES/beam.py (itself from PIC_MCC_3D/beam3d and
spectre/ebeam), with its provenance:

* Stopping power. The relativistic Bethe/Moller collisional formula with
  the ICRU 37 mean excitation energies, REBUILD's stopping_power line for
  line (checked there against NIST ESTAR to 2-6 % at 600 keV). It turns a
  printed beam voltage and current density into a deposited power
  (deposited_power_W_cm3: a straight ray through a thin cell, no detour or
  backscatter, the O(10 %) corrections REBUILD states).

* Yields per ion pair. REBUILD's Fowler/Spencer-Fano recursion
  (spencer_fano_per_process) run on the LXCat files the model's electron
  rows read (xsec.FILES: the IST-Lisbon Ar and SIGLO F2 sets, md5-identical
  to the ARPAE and LoKI-MC copies REBUILD reads), with REBUILD's prompt
  floor at the lowest Ar threshold, 11.548 eV. Below it, F2's vibrational
  and a3/A1 channels are rate-limited bulk processes that the network's
  electron rows already carry (REBUILD 0.1.12).

* The secondary's energy (found 2026-09-30). REBUILD splits what an
  ionization leaves, E - I, equally between the two electrons. On these
  sets that rule makes the yields per pair oscillate with the source
  energy: between 150 and 900 eV, E/U_ion runs 28 to 35 eV and the Ar
  excitations per ion 0.31 to 1.28, and some points move by x2 with the
  grid. REBUILD's single source energy (55 % of the tables' ceiling, 516
  eV with F2) is one point on that oscillation. The default here gives the
  secondary the measured Opal-Peterson-Beaty spectrum, p(eps) ~ 1 / (1 +
  (eps/Ebar)^2) on [0, (E - I)/2] (J. Chem. Phys. 55, 4100 (1971); Ebar(Ar)
  = 10.0 eV as LoKI-MC's Databases/OPBParameter.txt carries it, the same
  value PIC_MCC_3D/beam3d samples). The cascade then approaches its
  asymptote smoothly, and W and the yields are the same over every window
  from 150 to 900 eV to 0.1 %. sharing='equal' is REBUILD's rule.

* The route. A beam of 12 to 600 keV cascades at the asymptote, so the
  default reads the yields per pair as the least-squares slopes of the
  cumulative counts over 20-90 % of the tables' ceiling, from one sweep
  ('asymptotic'). route='integral' is REBUILD's integral ratio at 55 % of
  the ceiling. yields(x, 'integral', 'equal') reproduces REBUILD's
  spencer_fano_yields to round-off (tests/test_beam.py).

* The split onto LIGHT's species. Every ionization, of Ar or of F2, is one
  Ar+ and one electron: LIGHT carries no F2+, and the 0.07 eV between the
  two ionization energies comes out of the secondary's share (F2 makes 0.2
  % of the pairs at 0.5 % F2). Every Ar excitation is one Ar*, the 4s
  manifold at 11.548 eV, and the rest of its threshold is prompt heat: the
  cascade LIGHT's G1 and C1 rows book the same way. A prompt F2 C1 or H1
  excitation is dissociation into 2 F, with the energy above D0 = 1.6 eV
  as heat (REBUILD's reroute of the repulsive states). The secondaries'
  share eps_sec closes the identity

      W = sum_s y_s E_form(s) + heat + eps_sec

  exactly and goes to the electron gas.

W is an output. It is scored against the measured 26.4 eV per ion pair for
electrons in Ar (ICRU 31, REBUILD's anchor) and never set. With the
Opal-Peterson-Beaty spectrum these sets give 25.2 eV, 4.5 % low.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from . import xsec
from .species import RG, RG_ION, RG_STAR, e_form

# --- stopping power (REBUILD MODULES/beam.py) ------------------------------
MEC2_EV = 510998.95            # electron rest energy [eV]
R_E_M = 2.8179403262e-15       # classical electron radius [m]
KB = 1.380649e-23              # J/K
ATM_PA = 101325.0
# ICRU 37 mean excitation energy I [eV] and electrons per molecule; F2 as 2
# F atoms (I_F = 115 eV), as REBUILD
ICRU = {"Ar": (188.0, 18), "F2": (115.0, 18)}

# --- the cascade -------------------------------------------------------------
W_LIT_AR = 26.4                # eV per ion pair, electrons in Ar (ICRU 31)
SPLIT_4S = 12.9                # Ar 4s / 4p split [eV] (REBUILD _LO_HI_SPLIT)
D0_F2 = 2.0 * e_form("F")      # 1.6 eV
F2_ELEC_MIN = 2.5              # below: F2 vibrations (REBUILD)
# Opal-Peterson-Beaty Ebar [eV], LoKI-MC Databases/OPBParameter.txt: Ar 10.0
# (from Opal, Peterson & Beaty 1971). F2 is not in that table, whose rule
# for such a gas is the ionization potential, 15.69 eV (F2 ionizes 0.2 % of
# the pairs here)
OPB_EBAR = {"Ar": 10.0, "F2": 15.69}
SHARINGS = ("opb", "equal")
ROUTES = ("asymptotic", "integral")
N_GRID = {"opb": 4000, "equal": 1200}   # opb: W to 0.01 %, Ar* per pair to
                                        # 1 % (2000/4000/8000); equal:
                                        # REBUILD's grid
FIT = (0.20, 0.90)             # the asymptotic route's window, of the ceiling
INTEGRAL_AT = 0.55             # REBUILD's source energy, of the ceiling

# REBUILD's collision-frequency table: a fixed log grid read at the nearest
# point (cross_sections.CrossSectionSet), kept so the recursion is REBUILD's
_TAB_EPS = np.geomspace(1e-3, 2e3, 8192)
_TAB_LOG0 = math.log(_TAB_EPS[0])
_TAB_DLOG = math.log(_TAB_EPS[1] / _TAB_EPS[0])
_QE, _ME = 1.602176634e-19, 9.1093837015e-31
_V_TAB = np.sqrt(2.0 * _TAB_EPS * _QE / _ME)


# ------------------------------------------------------------ stopping power
def gas_densities_m3(pressure_atm: float, x_f2: float,
                     tg_K: float = 300.0) -> dict:
    n = pressure_atm * ATM_PA / (KB * tg_K)
    return {"Ar": (1.0 - x_f2) * n, "F2": x_f2 * n}


def stopping_power(densities_m3: dict, E_eV: float) -> float:
    """Collisional stopping power -dE/dx [eV/m] of a {gas: n_m3} mixture:
    relativistic Bethe/Moller, density effect 0, bremsstrahlung excluded
    (REBUILD beam.stopping_power)."""
    E_eV = float(E_eV)
    if E_eV <= 0.0:
        return 0.0
    tau = E_eV / MEC2_EV
    beta2 = 1.0 - 1.0 / (1.0 + tau) ** 2
    if beta2 <= 0.0:
        return 0.0
    total = 0.0
    for gas, n_gas in densities_m3.items():
        if n_gas <= 0.0:
            continue
        i_ev, z = ICRU[gas]
        f = (1.0 - beta2 + (tau ** 2 / 8.0 - (2.0 * tau + 1.0) * math.log(2.0))
             / (1.0 + tau) ** 2)
        big_l = math.log(tau ** 2 * (tau + 2.0)
                         / (2.0 * (i_ev / MEC2_EV) ** 2)) + f
        if big_l <= 0.0:
            continue
        total += (2.0 * math.pi * R_E_M ** 2 * MEC2_EV * n_gas * z / beta2) \
            * big_l
    return total


def csda_range_cm(pressure_atm: float, x_f2: float, E0_eV: float,
                  tg_K: float = 300.0, e_min_eV: float = 1.0e3,
                  n_steps: int = 400) -> float:
    """Continuous-slowing-down range [cm] (REBUILD beam.csda_range)."""
    if E0_eV <= e_min_eV:
        return 0.0
    dens = gas_densities_m3(pressure_atm, x_f2, tg_K)
    e = np.linspace(e_min_eV, E0_eV, int(n_steps))
    s = np.array([stopping_power(dens, x) for x in e])
    good = s > 0.0
    trap = getattr(np, "trapezoid", None) or np.trapz
    return 100.0 * float(trap(1.0 / s[good], e[good]))


def deposited_power_W_cm3(voltage_kV: float, j_A_cm2: float,
                          pressure_atm: float, x_f2: float,
                          tg_K: float = 300.0) -> float:
    """The power a beam of this voltage and current density (at the gas,
    after the foil) deposits per cm^3 while it crosses a cell thin next to
    its range: j S(E). A straight ray, so a lower bound by the detour."""
    s_ev_m = stopping_power(gas_densities_m3(pressure_atm, x_f2, tg_K),
                            voltage_kV * 1e3)
    return j_A_cm2 * s_ev_m * 1e-2          # A/cm2 x V/m -> W/cm3


# ------------------------------------------------------------- the cascade
class Cascade:
    """The inelastic processes of Ar + x_f2 F2 that the recursion needs,
    with REBUILD's tabulated collision frequencies. Only the composition
    matters: the recursion reads frequency ratios."""

    def __init__(self, x_f2: float):
        self.x_f2 = float(x_f2)
        frac = {"Ar": 1.0 - self.x_f2, "F2": self.x_f2}
        self.gases, self.kinds, self.products, thr, sig = [], [], [], [], []
        ends = []
        for gas in ("Ar", "F2"):
            if frac[gas] <= 0.0:
                continue
            procs = xsec.processes(gas)
            has_el = any(p.kind == "elastic" for p in procs)
            for p in procs:
                # REBUILD's ceiling: the last point of the elastic (or the
                # effective an elastic is synthesised from) and ionization
                if p.kind == "ionization" or p.kind == "elastic" or (
                        p.kind == "effective" and not has_el):
                    ends.append(float(p.energy[-1]))
                if p.kind not in ("excitation", "ionization"):
                    continue
                self.gases.append(gas)
                self.kinds.append(p.kind)
                self.products.append(p.product)
                thr.append(float(p.threshold))
                sig.append(frac[gas] * p.sigma_at(_TAB_EPS) * _V_TAB)
        self.thr = np.asarray(thr)
        self.nu_tab = np.asarray(sig)                  # (n_proc, n_tab)
        self.ceiling = min(ends)
        self.is_ion = np.array([k == "ionization" for k in self.kinds])
        self.ebar = np.array([OPB_EBAR[g] for g in self.gases])
        self.floor = float(self.thr[[g == "Ar" for g in self.gases]].min())

    def nu(self, e):
        i = int((math.log(max(e, _TAB_EPS[0])) - _TAB_LOG0) / _TAB_DLOG + 0.5)
        return self.nu_tab[:, min(max(i, 0), _TAB_EPS.size - 1)]


def sweep(cas: Cascade, e_top: float, sharing: str = "opb",
          n_grid: int = None):
    """Expected count of every prompt process per source electron at every
    energy of a uniform grid up to e_top: REBUILD's spencer_fano_per_process
    loop, which fills the grid from the bottom. An excitation of threshold
    I continues the electron at E - I. An ionization counts one pair and
    shares T = E - I: 'equal' gives each electron T/2 (REBUILD), 'opb' gives
    the secondary eps on [0, T/2] with the Opal-Peterson-Beaty weight (exact
    bin masses of its arctangent CDF) and the primary T - eps. Elastic loss
    is neglected (2m/M per event). Processes below the prompt floor are
    left to the network.

    Returns (grid, U (n_grid, n_proc)): zero columns for the left-out ones.
    """
    if sharing not in SHARINGS:
        raise ValueError(f"sharing {sharing!r}; one of {SHARINGS}")
    if e_top > cas.ceiling:
        raise ValueError(f"{e_top:.1f} eV is above the tables' ceiling "
                         f"{cas.ceiling:.1f} eV, where every sigma reads 0")
    n = int(n_grid or N_GRID[sharing])
    use = np.nonzero(cas.thr >= cas.floor)[0]
    thr, is_ion, ebar = cas.thr[use], cas.is_ion[use], cas.ebar[use]
    grid = np.linspace(0.0, float(e_top), n)
    de = grid[1] - grid[0]
    u = np.zeros((n, use.size))
    eye = np.eye(use.size)

    def at(e):
        pos = np.minimum(e / de, n - 1.000001)
        i0 = pos.astype(int)
        f = (pos - i0)[:, None]
        out = u[i0] * (1.0 - f) + u[np.minimum(i0 + 1, n - 1)] * f
        out[e <= 0.0] = 0.0
        return out

    for i in range(1, n):
        e = grid[i]
        nu = cas.nu(e)[use]
        tot = float(nu[thr <= e].sum())
        open_ = (thr <= e) & (nu > 0.0)
        if tot <= 0.0 or not open_.any():
            continue
        j = np.nonzero(open_)[0]
        w = nu[j] / tot
        rest = e - thr[j]
        acc = np.zeros(use.size)
        ex = ~is_ion[j]
        if ex.any():                                   # excitations
            acc += (w[ex][:, None] * (eye[j[ex]] + at(rest[ex]))).sum(axis=0)
        for k in np.nonzero(~ex)[0]:                   # ionizations
            t = rest[k]
            if sharing == "equal":
                share = 2.0 * at(np.array([0.5 * t]))[0]
            else:
                m = max(int(0.5 * t / de), 1)
                edges = np.linspace(0.0, 0.5 * t, m + 1)
                b = ebar[j[k]]
                cdf = np.arctan(edges / b) / math.atan(0.5 * t / b)
                eps = 0.5 * (edges[:-1] + edges[1:])
                share = (np.diff(cdf)[:, None] * (at(eps) + at(t - eps))
                         ).sum(axis=0)
            acc += w[k] * (eye[j[k]] + share)
        u[i] = acc
    full = np.zeros((n, cas.thr.size))
    full[:, use] = u
    return grid, full


def spencer_fano(cas: Cascade, E0_eV: float, sharing: str = "equal",
                 n_grid: int = None):
    """Counts per source electron at E0 and the sub-threshold terminal
    energy (REBUILD's spencer_fano_per_process with its prompt floor)."""
    grid, u = sweep(cas, E0_eV, sharing, n_grid)
    return u[-1], float(E0_eV) - float(np.dot(u[-1], cas.thr))


@dataclass(frozen=True)
class BeamYields:
    """Production per ion pair, with the identity attached."""
    x_f2: float
    route: str
    sharing: str
    window_eV: tuple        # asymptotic: the fit window; integral: (0, E0)
    W: float                # eV per ion pair
    ar_star: float          # Ar* per pair (all Ar excitations)
    ar_star_4s: float       # of which the 4s levels directly
    f2_diss: float          # prompt F2 dissociations per pair (C1, H1)
    f2_ion: float           # F2 ionizations per pair (booked as Ar+)
    heat: float             # prompt heat to the gas per pair [eV]
    eps_sec: float          # to the electron gas per pair [eV]

    def per_pair(self) -> dict:
        """Species made (+) and used (-) per ion pair."""
        return {RG_ION: 1.0, RG_STAR: self.ar_star, "F": 2.0 * self.f2_diss,
                "F2": -self.f2_diss, RG: -(1.0 + self.ar_star)}

    def identity_residual(self) -> float:
        spend = sum(n * e_form(s) for s, n in self.per_pair().items())
        return self.W - spend - self.heat - self.eps_sec

    def one_photon_ceiling(self, hv_eV: float) -> float:
        """The pure-Ar VUV efficiency if every ion and every excitation
        ends as one photon of hv_eV (the network's accounting)."""
        return hv_eV * (1.0 + self.ar_star) / self.W


@lru_cache(maxsize=256)
def yields(x_f2: float, route: str = "asymptotic",
           sharing: str = "opb") -> BeamYields:
    """Per-pair yields for Ar + x_f2 F2 (see the module docstring)."""
    if route not in ROUTES:
        raise ValueError(f"route {route!r}; one of {ROUTES}")
    cas = Cascade(x_f2)
    if route == "integral":
        e0 = INTEGRAL_AT * cas.ceiling
        grid, u = sweep(cas, e0, sharing)
        d, d_e, window = u[-1], e0, (0.0, e0)
    else:
        lo, hi = FIT[0] * cas.ceiling, FIT[1] * cas.ceiling
        grid, u = sweep(cas, hi, sharing)
        sel = grid >= lo
        d = np.polyfit(grid[sel], u[sel], 1)[0]          # counts per eV
        d_e, window = 1.0, (lo, hi)
    n_ion = float(d[cas.is_ion].sum())
    if n_ion <= 0.0:
        raise ValueError("the cascade made no ionization")
    y = d / n_ion
    w = d_e / n_ion
    ar_star = ar_4s = f2_diss = f2_ion = heat = 0.0
    e_star = e_form(RG_STAR)
    for k in range(cas.thr.size):
        if y[k] <= 0.0 or cas.thr[k] < cas.floor:
            continue
        thr = float(cas.thr[k])
        if cas.kinds[k] == "ionization":
            if cas.gases[k] == "F2":
                f2_ion += y[k]
        elif cas.gases[k] == "Ar":
            ar_star += y[k]
            if thr < SPLIT_4S:
                ar_4s += y[k]
            heat += y[k] * (thr - e_star)
        elif thr >= F2_ELEC_MIN:            # F2 C1 / H1: prompt dissociation
            f2_diss += y[k]
            heat += y[k] * (thr - D0_F2)
    spend = (e_form(RG_ION) + ar_star * e_star
             + 2.0 * f2_diss * e_form("F"))
    eps_sec = w - spend - heat
    if eps_sec < 0.0:
        raise ValueError(f"the beam split over-spends: eps_sec {eps_sec:.3f} eV")
    return BeamYields(float(x_f2), route, sharing, tuple(window), w, ar_star,
                      ar_4s, f2_diss, f2_ion, heat, eps_sec)
