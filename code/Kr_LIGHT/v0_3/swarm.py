"""The kinetic EEDF's solver run as a swarm in pure krypton (LIGHT_Kr v0_3).

No electron-electron collisions, the temporal growth term on, and every
grid sized so that the EEDF falls 11 to 25 decades across it (the
solver's own window, bracketed in its top energy). The drift velocity is
the Joule power over eE (the solver's own flux) and the ionization
coefficient alpha/N = k_ion / W. v0_3 adds the transverse diffusion,
N D_T = (gamma/3) int eps F / sigma~ d eps over the same sigma~ (momentum
transfer and growth) the field term uses (Hagelaar & Pitchford 2005), and
the characteristic energy D_T/mu. ML's swarm analysis trains v0_3's cross
sections on these against krypton's measured swarm data (xsec.Factors),
and the accuracy overview scores them.

hm is the network's band hm_kr: the share of the excitations of the
levels at or above E(Kr2+) that ionize associatively, Kr(nl) + Kr ->
Kr2+ + e (Hornbeck-Molnar), as the HM rows book it. A steady-state
Townsend experiment counts those electrons with the direct ones, so the
ionization coefficient at hm is (k_ion + k_assoc) / W; hm = 0 is direct
ionization alone, v0_1's swarm. Their electrons enter the equation as the
network's do (boltzmann 'assoc': born at rest, counted in the growth).
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from . import boltzmann as KB
from . import xsec
from .network import E_ASSOC_EV

N_COL = 2.5e25          # m^-3; a swarm's results are per molecule
NE_FRAC = 4e-29         # n_e/N: no electron-electron collisions to speak of


def processes(hm: float = 0.0) -> list:
    """Krypton's excitations and ionization as the solver's processes, the
    levels at or above E(Kr2+) split between their cascade (1 - hm) and
    their associative ionization (hm, kind 'assoc')."""
    if not 0.0 <= hm <= 1.0:
        raise ValueError(f"hm {hm}: a share, 0 to 1")
    out = []
    for p in xsec.processes("Kr"):
        if p.kind == "ionization":
            out.append(KB.Process("Kr", "ionization", p.threshold, p))
        elif p.kind == "excitation":
            if p.threshold >= E_ASSOC_EV and hm > 0.0:
                if hm < 1.0:
                    out.append(KB.Process("Kr", "excitation", p.threshold,
                                          replace(p, sigma=p.sigma
                                                  * (1.0 - hm))))
                out.append(KB.Process("Kr", "assoc", p.threshold,
                                      replace(p, sigma=p.sigma * hm)))
            else:
                out.append(KB.Process("Kr", "excitation", p.threshold, p))
    return out


def point(en_td: float, tg: float = 300.0, hm: float = 0.0,
          cells: int = 1000) -> dict:
    """The swarm at E/N [Td] and gas temperature tg [K]: drift velocity W
    [cm/s], alpha/N [cm^2] (direct and associative), its direct part, the
    transverse diffusion N D_T [1/(cm s)] and the characteristic energy
    D_T/mu [V], the mean energy [eV], and the grid's top [eV] and
    decades."""
    procs = processes(hm)
    emax = 1.5 if en_td < 0.01 else 4.0 * max(en_td, 0.1) ** 0.35 + 8.0
    lo = hi = None
    for _ in range(30):
        s = KB.Solver(procs, n_cells=cells, emax=emax, elastic=("Kr",),
                      tol=1e-10, maxit=4000)
        st = KB.State(x={"Kr": 1.0}, tg=tg, n_col=N_COL,
                      ne=NE_FRAC * N_COL)
        sol = s.solve(st, en_td)
        dec = KB._decades(sol.F)
        if 11.0 <= dec <= 25.0:
            break
        if dec < 11.0:
            lo = emax
            emax = emax * 1.6 if hi is None else math.sqrt(emax * hi)
        else:
            hi = emax
            emax = emax / 1.6 if lo is None else math.sqrt(lo * emax)
    else:
        raise RuntimeError(f"no grid in the window at {en_td} Td")
    w = s.field_power(st, sol) / (en_td * KB.TD * N_COL)          # m/s
    k_dir = sum(sol.grid.rate(q.p, sol.F) for q in procs
                if q.kind == "ionization")
    k_as = sum(sol.grid.rate(q.p, sol.F) for q in procs if q.kind == "assoc")
    # v0_3: N D_T = (gamma/3) int eps/sigma~ F d eps at the inner
    # boundaries, with the field term's own gamma eps / (3 sigma~)
    dfield, _, _ = s._split(st, sol)
    fb = 0.5 * (sol.F[1:] + sol.F[:-1])
    dt_n = float(np.sum(dfield * fb)) * sol.grid.du                # 1/(m s)
    mu_n = w / (en_td * KB.TD)                                      # 1/(V m s)
    return {"W": w * 100.0, "alphaN": (k_dir + k_as) / w * 1e4,
            "alphaN_direct": k_dir / w * 1e4, "DT_N": dt_n / 100.0,
            "DTmu": dt_n / mu_n, "mean": sol.mean,
            "emax": emax, "decades": dec}
