"""The kinetic EEDF's solver run as a swarm in pure krypton (LIGHT_Kr v0_3;
v0_4 adds the electroluminescence yield, el_yield).

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
    k_exc = {}
    for q in procs:
        if q.kind == "excitation":
            g = xsec.exc_group(q.p.product)
            k_exc[g] = k_exc.get(g, 0.0) + sol.grid.rate(q.p, sol.F)
    return {"W": w * 100.0, "alphaN": (k_dir + k_as) / w * 1e4,
            "alphaN_direct": k_dir / w * 1e4, "DT_N": dt_n / 100.0,
            "DTmu": dt_n / mu_n, "mean": sol.mean,
            "emax": emax, "decades": dec,
            # v0_4: the excitations per electron and per cm of drift, per
            # atom [cm^2], by xsec.exc_group ('1s' the 5s levels, '2p' the
            # 5p rung, 'hi' the rest): the electroluminescence yield's
            # count where each excitation ends as one photon (el_yield)
            "excN": {g: v / w * 1e4 for g, v in k_exc.items()}}


# v0_4: Mano et al., Phys. Lett. B 824, 136782 (2022): krypton's reduced
# electroluminescence yield at 1.1 bar and room temperature, the linear fit
# Y/N = 0.113 (E/N) - 0.30 [1e-17 photons per electron per cm per atom,
# cm^2] with E/N in Td, measured below 14 Td from its threshold, 2.7 Td.
# The training never saw it: a swarm-regime score of the excitation cross
# sections, with no electron-electron collisions and no excited atoms
MANO_EL = (0.113, -0.30)
MANO_EL_UNIT = 1e-17           # cm^2
MANO_EL_RANGE = (2.7, 14.0)    # Td
MANO_EL_BAR, MANO_EL_TG = 1.1, 296.0


def el_yield(en_td: float, photons=None, tg: float = MANO_EL_TG,
             cells: int = 1000) -> dict:
    """The swarm's reduced electroluminescence yield at E/N [Td]: the
    excitations per electron and per cm of drift per atom, each group's
    times its photons per excitation (photons: {group: yield}, by default
    one each), against Mano et al.'s fit, in units of 1e-17 cm^2. Below 20
    Td nothing ionizes and the flux drift velocity is the bulk one
    (swarm_mc)."""
    pt = point(en_td, tg=tg, cells=cells)
    photons = photons or {}
    y = sum(v * photons.get(g, 1.0) for g, v in pt["excN"].items())
    mano = MANO_EL[0] * en_td + MANO_EL[1]
    tot = sum(pt["excN"].values())
    return {"model": y / MANO_EL_UNIT, "mano": mano,
            "ratio": y / MANO_EL_UNIT / mano if mano > 0.0 else float("nan"),
            "share_5s": pt["excN"].get("1s", 0.0) / tot if tot else 0.0,
            "alphaN": pt["alphaN"] / MANO_EL_UNIT, "W": pt["W"]}
