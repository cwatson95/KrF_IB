"""The KrF(B,C) two-level block of the basic network in closed form.

LIGHT's basic network holds KrF in two coupled states, B and C (Gadomski,
Xu, Setser & Nelson, CPL 189, 153 (1992) eq. 5). With no power and no
electron production the B/C block is linear, so any seeded experiment
(photolysis, photoassociation, pulse radiolysis) has a closed-form answer:
two exponentials with rates lam+ and lam-. This module gives that answer
for any constants, so a published decay can be reproduced the way it was
measured (excitation pulse, detector response, fit window) and re-read
under the two-level model.

`Consts.from_light()` takes the constants from LIGHT's basic network, and
`light_seed_run()` runs LIGHT itself. `check_against_light()` measures how
far the two differ, and `tests/test_two_level.py` holds them to 1e-5, so a
result here is a result for LIGHT. The evidence scripts in docs/krf_data/
and the fitting project ML/ (beside LIGHT) use this module.

Densities in cm^-3, rates in s^-1, times in s. 1 Torr at 300 K is
3.2188e16 cm^-3.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

TORR_300K = 3.2188e16            # cm^-3 per Torr at 300 K
CM1_K = 1.4387769                # K per cm^-1


@dataclass(frozen=True)
class Consts:
    """Two-level constants. B and C share the quench constants unless the
    *_C fields are set (Gadomski's model has them equal)."""
    tau_B: float = 6.993e-9      # s
    tau_C: float = 75.0e-9       # s
    k_BC: float = 1.0e-10        # B -> C by Kr, cm^3/s
    K_eq: float = 5.0            # [C]/[B] at equilibrium, at T_ref
    gap_cm1: float = 0.0         # B-C gap for K_eq(T) by detailed balance
    T_ref: float = 300.0
    k_Kr: float = 4.0e-12        # KrF + Kr -> quench, cm^3/s
    k3_Kr: float = 6.0e-31       # KrF + 2Kr -> Kr2F* + Kr, cm^6/s (300 K)
    k3_tpow: float = 0.0         # T exponent of k3_Kr (LIGHT: -3.32)
    k_F2: float = 3.0e-10        # KrF + F2, cm^3/s
    k_e: float = 4.5e-8          # KrF + e, cm^3/s
    k_Kr_C: float | None = None
    k3_Kr_C: float | None = None
    k_F2_C: float | None = None
    k_e_C: float | None = None

    def keq(self, T):
        return self.K_eq * math.exp(self.gap_cm1 * CM1_K
                                    * (1.0 / T - 1.0 / self.T_ref))

    @classmethod
    def from_light(cls, **kw) -> "Consts":
        """The constants of LIGHT's basic network (network.basic_rows and
        Params.tau_sp_s), read from the rows so they cannot drift."""
        from .model import Params
        from .network import basic_rows
        from .species import e_form
        rows = {r.rid: r for r in basic_rows()}

        def const(rid):
            kind, val = rows[rid].rate
            assert kind == "const", rid
            return val
        k_bc = const("BC")
        kind, k_fwd, k_ref, upper, lower = rows["CB"].rate
        assert kind == "balance" and k_fwd == k_bc, rows["CB"].rate
        gap_cm1 = (e_form(upper) - e_form(lower)) / 1.239841984e-4
        base = dict(tau_B=Params().tau_sp_s, tau_C=1.0 / const("G4C"),
                    k_BC=k_bc, K_eq=k_ref, gap_cm1=gap_cm1,
                    k_Kr=const("QKR"), k3_Kr=const("G5"),
                    k3_tpow=rows["G5"].tg_pow or 0.0,
                    k_F2=const("QF2"), k_e=const("QE"),
                    k_Kr_C=const("QKRC"), k3_Kr_C=const("G5C"),
                    k_F2_C=const("QF2C"), k_e_C=const("QEC"))
        assert rows["G5C"].tg_pow == rows["G5"].tg_pow
        base.update(kw)
        return cls(**base)

    def with_(self, **kw) -> "Consts":
        return replace(self, **kw)


def rates(c: Consts, kr, f2=0.0, ne=0.0, T=300.0, extra=0.0):
    """K1 (B loss other than to C), K2 (B -> C), Km2 (C -> B), K3 (C loss
    other than to B). ``extra`` is an added loss on both states (a
    quencher the model does not carry, such as KrF2)."""
    k3 = c.k3_Kr * (T / 300.0) ** c.k3_tpow
    k3c = (c.k3_Kr_C if c.k3_Kr_C is not None else c.k3_Kr) \
        * (T / 300.0) ** c.k3_tpow
    qB = c.k_Kr * kr + k3 * kr * kr + c.k_F2 * f2 + c.k_e * ne + extra
    qC = ((c.k_Kr_C if c.k_Kr_C is not None else c.k_Kr) * kr
          + k3c * kr * kr
          + (c.k_F2_C if c.k_F2_C is not None else c.k_F2) * f2
          + (c.k_e_C if c.k_e_C is not None else c.k_e) * ne + extra)
    K2 = c.k_BC * kr
    return 1.0 / c.tau_B + qB, K2, K2 / c.keq(T), 1.0 / c.tau_C + qC


def eigen(c: Consts, kr, f2=0.0, ne=0.0, T=300.0, extra=0.0):
    """lam+ and lam- (s^-1), Gadomski eq. 5 and its partner."""
    K1, K2, Km2, K3 = rates(c, kr, f2, ne, T, extra)
    a, d = K1 + K2, Km2 + K3
    s = math.sqrt((a - d) ** 2 + 4.0 * K2 * Km2)
    return 0.5 * (a + d + s), 0.5 * (a + d - s)


def impulse(c: Consts, t, B0=1.0, C0=0.0, **cond):
    """B(t), C(t) after an instantaneous seed (B0, C0)."""
    K1, K2, Km2, K3 = rates(c, **cond)
    lp, lm = eigen(c, **cond)
    a, d = K1 + K2, Km2 + K3
    # B = bp e^-lp t + bm e^-lm t with B(0) = B0, B'(0) = -a B0 + Km2 C0
    dB0 = -a * B0 + Km2 * C0
    bp = (-dB0 - lm * B0) / (lp - lm)
    bm = B0 - bp
    dC0 = K2 * B0 - d * C0
    cp = (-dC0 - lm * C0) / (lp - lm)
    cm = C0 - cp
    ep, em = np.exp(-lp * t), np.exp(-lm * t)
    return bp * ep + bm * em, cp * ep + cm * em


def steady_fraction_B(c: Consts, **cond) -> float:
    """[B]/([B]+[C]) under a steady source that forms B only."""
    K1, K2, Km2, K3 = rates(c, **cond)
    c_over_b = K2 / (Km2 + K3)
    return 1.0 / (1.0 + c_over_b)


def manifold_rate(c: Consts, source_B=1.0, source_C=0.0, **cond):
    """Steady state under a constant source: returns (B, C, total loss
    rate of the B/C manifold per molecule, excluding B<->C exchange)."""
    K1, K2, Km2, K3 = rates(c, **cond)
    M = np.array([[K1 + K2, -Km2], [-K2, Km2 + K3]])
    B, C = np.linalg.solve(M, np.array([source_B, source_C]))
    return B, C, (source_B + source_C) / (B + C)


def sv_tau(c: Consts, source_B=1.0, source_C=0.0, **cond) -> float:
    """The lifetime a Stern-Volmer electron-quench analysis of steady B-X
    emission measures: d(I0/I)/d(n_e) at n_e = 0, divided by k_e. For one
    state it is that state's lifetime. For fast B/C mixing with equal
    quench constants it is the manifold lifetime."""
    cond = dict(cond, ne=0.0)
    K1, K2, Km2, K3 = rates(c, **cond)
    kC = c.k_e_C if c.k_e_C is not None else c.k_e
    M = np.array([[K1 + K2, -Km2], [-K2, Km2 + K3]])
    R = np.array([source_B, source_C])
    x = np.linalg.solve(M, R)                  # steady B, C with no e
    D = np.diag([c.k_e, kC])
    dx = -np.linalg.solve(M, D @ x)            # d(B, C)/d(n_e) at 0
    return float(-dx[0] / x[0] / c.k_e)        # d(B0/B)/dn = -dB/B0


def instrument(t, y, fwhm_s):
    """Convolve a uniformly sampled impulse response with a Gaussian of
    the given FWHM (excitation pulse and detector folded together)."""
    if fwhm_s <= 0.0:
        return y.copy()
    dt = t[1] - t[0]
    sig = fwhm_s / 2.354820045
    tt = np.arange(-6.0 * sig, 6.0 * sig + 0.5 * dt, dt)
    g = np.exp(-0.5 * (tt / sig) ** 2)
    g /= g.sum()
    pad = np.concatenate([np.zeros(len(tt)), y, np.zeros(len(tt))])
    out = np.convolve(pad, g, mode="same")
    return out[len(tt):len(tt) + len(y)]


def fit_exp(t, y, hi=0.5, lo=0.05, t_min=None, t_max=None):
    """Single-exponential rate fitted to y on its falling side, between the
    fractions hi and lo of the peak (or between t_min and t_max)."""
    i0 = int(np.argmax(y))
    sel = np.arange(len(y)) > i0
    if t_min is not None or t_max is not None:
        sel &= (t >= (t_min if t_min is not None else -np.inf))
        sel &= (t <= (t_max if t_max is not None else np.inf))
    else:
        sel &= (y < hi * y[i0]) & (y > lo * y[i0])
    if sel.sum() < 3:
        raise ValueError("fit window holds fewer than 3 points")
    return -np.polyfit(t[sel], np.log(y[sel]), 1)[0]


def light_seed_run(kr_torr, f2_torr=0.0, B0=1e12, C0=0.0, t_end=300e-9,
                   n_save=3001, overrides=None):
    """LIGHT's basic network in seed mode (no power, free n_e = 1 cm^-3,
    300 K). ``overrides`` maps row id -> new constant, for sensitivity
    runs only. Returns t, KrF(B), KrF(C)."""
    from .model import Model, Params
    tot = kr_torr + f2_torr
    seeds = [("KrF(B)", B0)] + ([("KrF(C)", C0)] if C0 else [])
    p = Params(network="basic", pressure_atm=tot / 760.0,
               x_f2=f2_torr / tot if tot else 0.0, power_W_cm3=0.0,
               pulse_s=0.0, t_end_s=t_end, ne_mode="free", ne_cm3=1.0,
               n_save=n_save, init_extra_cm3=tuple(seeds))
    m = Model(p)
    for rid, val in (overrides or {}).items():
        hit = [cc for cc in m.compiled if cc.row.rid == rid]
        assert hit and hit[0].kind == "const", rid
        hit[0].rate = val               # const rows only (not CB, "balance")
    r = m.run()
    return r.t, r.n["KrF(B)"], r.n["KrF(C)"]


def check_against_light(kr_torr=100.0, f2_torr=1.0):
    """Max relative difference between LIGHT's run and the closed form
    for a KrF(B) seed, over the first 100 ns."""
    c = Consts.from_light()
    t, B, C = light_seed_run(kr_torr, f2_torr)
    tot = kr_torr + f2_torr
    # LIGHT's gas density is set by p/(k Tg) at 300 K, as TORR_300K is
    kr = kr_torr * TORR_300K
    f2 = f2_torr * TORR_300K
    Ba, Ca = impulse(c, t, B0=B[0], C0=0.0, kr=kr, f2=f2)
    sel = t <= 100e-9
    return float(max(np.max(np.abs(B[sel] - Ba[sel]) / B[0]),
                     np.max(np.abs(C[sel] - Ca[sel]) / B[0])))
