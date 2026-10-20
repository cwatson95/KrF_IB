"""Electron collisions with excited krypton that nothing measures, and the
5p rung's levels (LIGHT_Kr v0_5, 2026-10-16).

Three things the network's v0_5 rows read:

* **Vriens & Smeets's ionization of an excited level** (Phys. Rev. A 22,
  940 (1980)), the owner's choice of 2026-10-16 for the rows with no
  measurement (the stepwise ionizations out of the 5s levels and the 5p
  rung, and the ionization of Kr2*, which Flannery, McCann & Winter, J.
  Phys. B 14, 3789 (1981) put at about 1e-15 cm2 by binary-encounter
  theory). The cross section of a level bound by E_n,

      sigma(E) = 4 pi a0^2 R^2 / (E + 3.25 E_n)
                 x (5 / (3 E_n) - 1/E - 2 E_n / (3 E^2)),     E > E_n,

  and their rate formula, its Maxwellian average,

      k(Te) = 9.56e-6 Te^-3/2 exp(-e) / (e^2.33 + 4.38 e^1.72 + 1.32 e)
              cm3/s,  e = E_n / Te  (energies in eV).

  The two agree to 0.1 % for Te = 1-5 eV and 3 % at 0.5 eV (checked when
  v0_5 was built, CHANGELOG 2026-10-16 (b), and by the tests). The kinetic
  EEDF reads the cross section, the two-temperature closure its
  Maxwellian table. Kannari et al.'s printed fits of the same rows (JAP
  57, 4309 (1985) Table I [72], [79]) lie 10-40 % from it over 1-3 eV.

* **Pseudo cross sections** that carry a rate law into the kinetic EEDF:
  A Te^B exp(-C/Te) is exactly the Maxwellian average of

      sigma(eps) = sigma0 (eps - C)^(B + 1/2) / eps,   eps > C,
      sigma0 = A / (sqrt(2 e / m) Gamma(B + 3/2) (2 / sqrt pi)),

  (v0_4's form for B = 1/2 is the case sigma0 (1 - C/eps)), and a constant
  rate k of a superelastic collision that gives the electron u is the
  reverse, by the Klein-Rosseland relation at g, of a forward cross
  section k (E - u) / (g E v(E - u)).

* **The 5p rung's levels.** The rung Kr** is the ten 4p5 5p levels,
  booked at the lowest (2p10). v0_5 takes them populated at the gas
  temperature: krypton mixes them faster than they radiate (Gazeli et al.'s
  1.7e-10 cm3/s total loss of 2p2 is 4e9 /s at 1 atm, against 3e7 /s of
  radiation; Whitehead et al., JCP 102, 1965 (1995) measured the 2p8 <->
  2p9 mixing in krypton). A row that leaves from one level reads the
  rung times that level's share f_x(Tg) = g_x exp(-(E_x - E_2p10)/kTg) / Z;
  the level's excitation above 2p10, which the ledger does not book in the
  rung, comes from the gas when the level de-excites (network.py).
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import xsec

A0_CM = 0.529177210903e-8          # Bohr radius
RYD_EV = 13.605693122994           # Rydberg energy
KB_EV = 1.380649e-23 / 1.602176634e-19


# --------------------------------------------------------- Vriens & Smeets
def vs_sigma(e, e_n: float) -> np.ndarray:
    """Vriens & Smeets's ionization cross section [cm2] of a level bound
    by e_n [eV], at electron energies e [eV]; 0 at and below e_n."""
    e = np.asarray(e, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = (4.0 * math.pi * A0_CM ** 2 * RYD_EV ** 2 / (e + 3.25 * e_n)
             * (5.0 / (3.0 * e_n) - 1.0 / e - 2.0 * e_n / (3.0 * e * e)))
    return np.where(e > e_n, np.maximum(s, 0.0), 0.0)


def vs_rate(e_n: float, te: float) -> float:
    """Their rate formula [cm3/s] at Te [eV]."""
    x = e_n / te
    return (9.56e-6 * te ** -1.5 * math.exp(-x)
            / (x ** 2.33 + 4.38 * x ** 1.72 + 1.32 * x))


def _grid(lo: float) -> np.ndarray:
    return np.concatenate([[0.0, lo], np.geomspace(lo * (1.0 + 1e-9), 1000.0,
                                                   700)])


@lru_cache(maxsize=None)
def vs_process(e_n: float, label: str = "excited") -> xsec.Process:
    """The cross section as an ionization Process (threshold e_n)."""
    e = _grid(e_n)
    return xsec.Process("ionization", label, "Kr+ (Vriens-Smeets)", e_n, e,
                        vs_sigma(e, e_n))


@lru_cache(maxsize=None)
def vs_rate_table(e_n: float) -> xsec.RateTable:
    """Its Maxwellian <sigma v>(Te) [cm3/s] on xsec's grid."""
    p = vs_process(e_n)
    k = xsec._mean(p.sigma_at(xsec.EPS) * xsec.V_CM_S * np.sqrt(xsec.EPS))
    return xsec.RateTable(k, f"Vriens-Smeets ionization, E_n {e_n:.4f} eV")


@lru_cache(maxsize=None)
def vs_capture_table(e_n: float) -> xsec.EnergyTable:
    """<eps sigma v> / <sigma v> [eV] of the ionizing electrons."""
    p = vs_process(e_n)
    sv = p.sigma_at(xsec.EPS) * xsec.V_CM_S * np.sqrt(xsec.EPS)
    k = xsec._mean(sv)
    ke = xsec._mean(sv * xsec.EPS)
    ok = k > xsec._K_FLOOR
    return xsec.EnergyTable(np.where(ok, ke / np.where(ok, k, 1.0),
                                     1.5 * xsec.TE),
                            f"Vriens-Smeets E_n {e_n:.4f}: ionizing energy")


# ------------------------------------------------------- pseudo processes
def arrhenius_sigma0(a_cm3: float, b: float) -> float:
    """sigma0 [cm2] of (eps - C)^(B + 1/2) / eps for A Te^B exp(-C/Te)."""
    return a_cm3 / (xsec.V_CM_S * math.gamma(b + 1.5) * 2.0
                    / math.sqrt(math.pi))


@lru_cache(maxsize=None)
def arrhenius_process(a_cm3: float, b: float, c_ev: float, loss_ev: float,
                      target: str = "excited", product: str = "excited"):
    """An excitation Process whose Maxwellian rate is exactly A Te^B
    exp(-C/Te); the electron loses loss_ev (the energy gap, C itself for a
    transfer between two levels)."""
    e = _grid(c_ev)
    s0 = arrhenius_sigma0(a_cm3, b)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(e > c_ev, s0 * np.maximum(e - c_ev, 0.0) ** (b + 0.5)
                     / np.where(e > 0, e, 1.0), 0.0)
    return xsec.Process("excitation", target, product, loss_ev, e, s)


@lru_cache(maxsize=None)
def gain_process(k_cm3: float, u_ev: float, g_ratio: float = 1.0,
                 target: str = "excited"):
    """The forward Process whose Klein-Rosseland reverse at g_ratio is a
    superelastic collision of constant rate k [cm3/s] that gives the
    electron u_ev: sigma_f(E) = k (E - u) / (g E v(E - u)), E > u."""
    e = _grid(u_ev)
    x = np.maximum(e - u_ev, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(e > u_ev, k_cm3 * x / (g_ratio * e * xsec.V_CM_S
                                            * np.sqrt(np.where(x > 0, x, 1.0))),
                     0.0)
    return xsec.Process("excitation", target, "ground", u_ev, e, s)


# --------------------------------------------------------- the 5p rung
# the ten 4p5 5p levels: Paschen label -> (Biagi's threshold [eV], which is
# NIST's level, and 2J + 1). 2p1-2p4 have the 2P1/2 core, 2p5-2p10 the
# 2P3/2 one
RUNG = {
    "2p10": (11.3035, 3), "2p9": (11.4430, 7), "2p8": (11.4447, 5),
    "2p7": (11.5261, 3), "2p6": (11.5458, 5), "2p5": (11.6660, 1),
    "2p4": (12.1004, 3), "2p3": (12.1404, 3), "2p2": (12.1437, 5),
    "2p1": (12.2565, 1),
}
RUNG_BOTTOM = "2p10"


@lru_cache(maxsize=65536)
def _rung(tg_r: float):
    kt = KB_EV * tg_r
    e0 = RUNG[RUNG_BOTTOM][0]
    w = {x: g * math.exp(-(e - e0) / kt) for x, (e, g) in RUNG.items()}
    z = sum(w.values())
    return z, {x: v / z for x, v in w.items()}


def rung_partition(tg: float) -> float:
    """Z(Tg) = sum g_x exp(-(E_x - E_2p10)/kTg), the rung's levels against
    its bottom (3 at 0 K, 36 at infinite Tg)."""
    return _rung(round(max(float(tg), 50.0), 3))[0]


def rung_fraction(level: str, tg: float) -> float:
    """f_x(Tg): the share of the rung in one level."""
    return _rung(round(max(float(tg), 50.0), 3))[1][level]


def rung_mean_excess(tg: float) -> float:
    """<E_x - E_2p10> over the rung's levels at Tg [eV]."""
    z, f = _rung(round(max(float(tg), 50.0), 3))
    e0 = RUNG[RUNG_BOTTOM][0]
    return sum(f[x] * (RUNG[x][0] - e0) for x in RUNG)
