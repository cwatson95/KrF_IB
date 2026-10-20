"""The channels being fitted: the KrF(B,C) block of LIGHT's basic network.

Ten channels. Rates are log-uniform inside wide bounds. tau_B carries a
Gaussian prior of 6.9 +- 0.3 ns (Eden measured 6.8 +- 0.2, theory 7.0).
Held fixed at LIGHT's values:
  - the B-C gap in the energy ledger, which sets K_eq's temperature
    dependence (80 cm^-1 in v0, 120 in v0_1, species.py);
  - the trimer's (Tg/300)^-3.32 law;
  - the electron quench, which is re-derived for each set
    (forward.with_rederived_k_e).
"""
from __future__ import annotations

import numpy as np

from . import light

L = light("two_level")

NAMES = ("k_BC", "K_eq", "tau_C", "tau_B", "k_Kr", "k_Kr_C", "k3_Kr",
         "k3_Kr_C", "k_F2", "k_F2_C")
LO = np.log([1e-13, 0.2, 20e-9, 6.0e-9, 1e-14, 1e-14, 1e-33, 1e-33, 5e-11,
             5e-11])
HI = np.log([1e-9, 20.0, 150e-9, 7.8e-9, 5e-11, 5e-11, 3e-30, 3e-30, 2e-9,
             2e-9])
TAU_B_PRIOR = (6.9e-9, 0.3e-9)


def light_now() -> L.Consts:
    return L.Consts.from_light()


def gadomski_1992() -> L.Consts:
    """Gadomski, Xu, Setser & Nelson 1992's published two-level set: tau_C
    75 ns, k_BC 1e-10, K_eq 5.0, one quench set for B and C (k_Kr 4e-12,
    k3 6e-31, k_F2 3e-10), k_e 4.5e-8, with LIGHT's ledger gap and trimer
    temperature law."""
    return light_now().with_(tau_C=75e-9, k_BC=1.0e-10, K_eq=5.0,
                             k_Kr=4.0e-12, k3_Kr=6.0e-31, k_F2=3.0e-10,
                             k_e=4.5e-8, k_Kr_C=None, k3_Kr_C=None,
                             k_F2_C=None, k_e_C=None)


def light_to_0927() -> L.Consts:
    """LIGHT's basic-network constants until 2026-09-28: Gadomski's quench
    set and K_eq with Hakuta's tau_C (60 ns) and k_BC (1.2e-10)."""
    return gadomski_1992().with_(tau_C=60e-9, k_BC=1.2e-10)


def consts(x) -> L.Consts:
    """A channel vector (natural logs) as two-level constants, with LIGHT's
    ledger gap and trimer temperature law."""
    now = light_now()
    return L.Consts(gap_cm1=now.gap_cm1, k3_tpow=now.k3_tpow, k_e=now.k_e,
                    **dict(zip(NAMES, np.exp(np.asarray(x, dtype=float)))))


def vector(c: L.Consts) -> np.ndarray:
    """The channel vector of a constant set (C values default to B's)."""
    vals = []
    for n in NAMES:
        v = getattr(c, n)
        if v is None:
            v = getattr(c, n[:-2])
        vals.append(v)
    return np.log(vals)
