"""Analytic END-TO-END trajectory gates — the plasma-oscillation-gate analog.

The PIC_MCC_3D gap-closure (2026-08-04): the suite gated the RHS *algebra*
per evaluation but never integrated a network against a closed-form
solution. These gates integrate REAL network rows (pulled by rid from
network.py, not synthetic stand-ins) through the REAL CompiledNetwork with
scipy BDF and compare whole trajectories — densities AND the gas-heating
they release — against exact solutions.

Tolerance discipline (root Known-bad: never integrate loose and compare
tight): rtol = 1e-11 here, comparisons at 1e-6 — three decades of slack
between the integrator's accuracy and the gate.

Electron temperature is held at 0.026 eV (no electron-impact rows are in
these sub-networks; the only Te-dependent terms are the tapered elastic/vib
drains, negligible at the quasineutral floor densities involved and checked
harmless by the Tg gates)."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from REBUILD.MODULES.chemistry import CompiledNetwork
from REBUILD.MODULES.network import (ar_kr_couplings, build_network,
                                     kr_core, mode_optics,
                                     upper_level_radiation)
from REBUILD.MODULES.species import KRF_MODE_SPECIES, StateLayout
from REBUILD.MODULES.units import EV, K_B, photon_energy_J

NTOT = 2.4e25          # heat-capacity reference density for the Tg equation
TG0 = 300.0


def _row(rid):
    for r in (kr_core() + upper_level_radiation("krf")
              + ar_kr_couplings()):
        if r.rid == rid:
            return r
    raise KeyError(rid)


def _integrate(rows, init_dens, t_end_s, tiers, n_eval=9, rtol=1e-11):
    """Integrate [dens..., U_e, T_g] for one cell through the real
    CompiledNetwork; returns (t, dens(t) dict, Tg(t))."""
    lay = StateLayout(KRF_MODE_SPECIES, 1)
    cn = CompiledNetwork(rows, lay, tiers, mode_optics("krf"))
    cn.set_photon_energy(photon_energy_J(248.0) / EV)
    S = lay.n_species
    y0 = np.zeros(S + 2)
    for s, v in init_dens.items():
        y0[lay.index[s]] = v
    y0[S + 1] = TG0

    def rhs(t, y):
        dens = np.maximum(y[:S], 0.0).reshape(S, 1)
        tg = np.clip(y[S + 1], 200.0, 3000.0)
        ne = lay.electron_density(dens)
        out = cn.rhs(dens, np.array([0.026]), np.array([tg]), ne, 0.0, NTOT)
        dy = np.empty_like(y)
        dy[:S] = out["ddens"][:, 0]
        dy[S] = out["du_e"][0]
        dy[S + 1] = out["dgas"][0] / (1.5 * NTOT * K_B)
        return dy

    t_eval = np.linspace(0.0, t_end_s, n_eval)
    sol = solve_ivp(rhs, (0.0, t_end_s), y0, method="BDF", rtol=rtol,
                    atol=1e-2, t_eval=t_eval)
    assert sol.success
    dens_t = {s: sol.y[lay.index[s]] for s in lay.species}
    return sol.t, dens_t, sol.y[S + 1]


def test_pure_radiative_decay_is_exponential(loki_ar_fill):
    """One row (URAD_B: KrF*(B) -> Kr + F at 1.43e8 s^-1, the Kannari
    B-state rate): n(t) = n0 e^{-t/tau}, fragments grow as
    n0 (1 - e^{-t/tau}). Compared over four lifetimes."""
    n0, tau = 1.0e19, 1.0 / 1.43e8
    urad_b = _row("URAD_B")
    t, d, _ = _integrate([urad_b], {"KrF*(B)": n0},
                         4.0 * tau, loki_ar_fill)
    ref = n0 * np.exp(-t / tau)
    assert np.allclose(d["KrF*(B)"], ref, rtol=1e-7)
    assert np.allclose(d["Kr"][1:], n0 * (1.0 - np.exp(-t[1:] / tau)),
                       rtol=1e-7)
    assert np.allclose(d["F"], d["Kr"], rtol=1e-9, atol=1e-3)


def test_second_order_harpoon_matches_closed_form(loki_ar_fill):
    """One row (AK6: Ar2* + F2 -> ArF* + Ar + F, k = 3.0e-16 m^3/s,
    Kannari [61] first-hand), unequal initial concentrations. Closed form
    for A + B -> products:
        b(t) = b0 (a0 - b0) / (a0 e^{(a0-b0) k t} - b0),  a = b + (a0 - b0)
    and the released heat is release_eV x (number of events), so
        dTg = release * e * (a0 - a(t)) / (3/2 n k_B)  — the trajectory-level
    twin of the per-evaluation every-eV-routed identity. (Was K25 — its
    KrF* product is now B/C/Bv-split; AK6 keeps the single-row form.)"""
    row = _row("AK6")
    k = row.rate[1]
    a0, b0 = 2.0e22, 1.0e22
    t_end = 400e-9
    t, d, tg = _integrate([row], {"Ar2*": a0, "F2": b0}, t_end, loki_ar_fill)
    dcon = a0 - b0
    b_ref = b0 * dcon / (a0 * np.exp(dcon * k * t) - b0)
    a_ref = b_ref + dcon
    assert np.allclose(d["Ar2*"], a_ref, rtol=1e-6)
    assert np.allclose(d["F2"], b_ref, rtol=1e-6)
    assert np.allclose(d["ArF*"][1:], (a0 - a_ref)[1:], rtol=1e-6)
    release = row.release_eV()
    assert release > 0.0
    tg_ref = TG0 + release * EV * (a0 - a_ref) / (1.5 * NTOT * K_B)
    assert np.allclose(tg, tg_ref, rtol=1e-6)


def test_three_body_dimerization_pseudo_first_order(loki_ar_fill):
    """One row (K15: Kr+ + 2Kr -> Kr2+ + Kr, 2.7e-43 m^6/s, Tg^-0.75).
    With the third body 2.4e7x the ion density, Kr is constant to 8e-8 and
    the ion decays exponentially at k n_Kr^2. MEASURED (2026-08-04, the
    reason for the 1e-5 pin): the trajectory deviates from the constant-k
    exponential by 2.9e-6 — which is the REAL Tg^-0.75 back-reaction (the
    1.1 eV/event release heats the gas 3.5e-4 K, lowering k by ~1 ppm).
    The coupling working at the ppm level is the finding; a stoichiometry
    or rate bug sits orders of magnitude above this pin."""
    row = _row("K15")
    k3 = row.rate[1]
    n_kr, n_ion = 2.4e25, 1.0e18
    lam = k3 * n_kr * n_kr
    t_end = 4.0 / lam
    t, d, tg = _integrate([row], {"Kr": n_kr, "Kr+": n_ion}, t_end,
                          loki_ar_fill)
    assert np.allclose(d["Kr+"], n_ion * np.exp(-lam * t), rtol=1e-5)
    assert np.allclose(d["Kr2+"][1:],
                       n_ion * (1.0 - np.exp(-lam * t[1:])), rtol=1e-5)
    # the release heating itself, against its own closed form
    assert abs((tg[-1] - TG0)
               / (row.release_eV() * EV * n_ion * (1.0 - math.exp(-4.0))
                  / (1.5 * NTOT * K_B)) - 1.0) < 1e-3


def test_device_drive_integral_matches_erf():
    """Device-level closed form: the ledger's drive_in channel for the
    Gaussian discharge waveform must equal P0 V sqrt(pi/a)/fwhm-scaled erf
    over the finite window — gates the channel integration and the drive
    wiring end to end through the real engine."""
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    P0, t0_ns, fwhm = 2.0e11, 30.0, 20.0
    p.set_many({"mode": "krf", "drive": "discharge", "t_end_ns": 60.0,
                "dis_peak_W_m3": P0, "dis_t0_ns": t0_ns,
                "dis_fwhm_ns": fwhm, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75},
               source="problem")
    dev = Device(p)
    res = dev.run(n_save=7)
    assert res["success"]
    a = 2.772588722 / (fwhm * 1e-9) ** 2          # exponent: -a (t - t0)^2
    sa = math.sqrt(a)
    t0 = t0_ns * 1e-9
    integral = (P0 * math.sqrt(math.pi / a) / 2.0
                * (math.erf(sa * (60e-9 - t0)) + math.erf(sa * t0)))
    vol = float(dev.cell_vols.sum())
    drive = res["final_channels_J"]["drive_in"]
    assert abs(drive / (integral * vol) - 1.0) < 1e-4
