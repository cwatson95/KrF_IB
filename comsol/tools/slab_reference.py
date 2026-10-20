"""The answer the COMSOL slab should reproduce: LIGHT's afterglow in every
cell of a 1-D slab, coupled only through the wave.

Each cell runs LIGHT's own right-hand side with the kit's frozen closure
(closure.blend_class): the same tables COMSOL reads. The wave is a plane
wave launched from vacuum at x = 0 into gas of depth L, solved exactly
through the stack of cells for piecewise-constant permittivity (transfer
matrices) with each cell's Drude permittivity, exp(+j w t):

    sigma = n_e e^2 / (m (nu + j w)),   eps = 1 - j sigma / (w eps0),

nu the closure's effective collision frequency (closure.nu_eff). A cell
absorbs the drop of the Poynting flux across it. Behind the gas: vacuum
(the transmitted wave is lost) or a mirror (back = 'pec'). Cells are much
thicker than a wavelength, so a mirror's standing wave is averaged over
each cell (the electrons' energy relaxation length is of order lambda/4).

The off twin is the same slab without the wave, the 0-D run times L (the
afterglow is uniform without it). Both twins start from LIGHT's kinetic
state at the end of the pulse and share its pump.

    python3 tools/slab_reference.py KrF 1.3 0.3            # f [THz], I0 [MW/cm2]
    python3 tools/slab_reference.py Kr 6.2 3 --L 10 --cells 30 --back pec
    python3 tools/slab_reference.py KrF --zero-d           # the 0-D twins
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import block_diag, csr_matrix

import closure as C
import make_kit as K

QE = C.QE
ME = C.ME
EPS0 = 8.8541878128e-12
C0 = 299792458.0
Z0 = 376.730313668
MAIN_BAND = {"KrF": "248 nm (KrF B-X)", "Kr": "147 nm (Kr2* VUV)"}
OUT = Path(__file__).resolve().parent.parent / "reference" / "runs"


def tmm(eps, dx, k0, back="vacuum"):
    """Fields through layers of permittivity eps (exp(+jwt)) and thickness
    dx [m] for a unit incident wave from vacuum: the Poynting flux at every
    interface over the incident flux (S[0] = 1 - R, S[-1] = T), the share
    each layer absorbs, and R."""
    n = np.sqrt(eps.astype(complex))
    n = np.where(n.imag > 0, -n, n)              # Im n < 0: decays forward
    E, H = (0.0 + 0j, 1.0 + 0j) if back == "pec" else (1.0 + 0j, 1.0 + 0j)
    fields = [(E, H)]
    for i in range(len(eps) - 1, -1, -1):
        d = k0 * n[i] * dx[i]
        c, s = np.cos(d), np.sin(d)
        E, H = c * E + 1j * s / n[i] * H, 1j * n[i] * s * E + c * H
        scale = max(abs(E), abs(H))
        E, H = E / scale, H / scale
        fields = [(f[0] / scale, f[1] / scale) for f in fields]
        fields.insert(0, (E, H))
    E0, H0 = fields[0]
    inc, ref = 0.5 * (E0 + H0), 0.5 * (E0 - H0)
    S = np.array([0.5 * (e * np.conj(h)).real for e, h in fields])
    S = S / (0.5 * abs(inc) ** 2)
    return S, S[:-1] - S[1:], abs(ref / inc) ** 2


def cells_for(L_m, ncell, first=20e-6):
    """Geometric cells from `first` at the front, summing to L."""
    if ncell * first >= L_m:
        return np.full(ncell, L_m / ncell)
    lo, hi = 1.0, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        tot = first * (r ** ncell - 1.0) / (r - 1.0)
        lo, hi = (r, hi) if tot < L_m else (lo, r)
    c = first * r ** np.arange(ncell)
    return c * L_m / c.sum()


class Setup:
    """The gas's LIGHT modules, runs, tables and start state."""

    def __init__(self, gas):
        self.gas = gas
        self.mods = C.load(gas)
        self.R = K.runs(gas, self.mods)
        self.tabs = K.tables(gas, self.mods, self.R)
        self.p_full = C.reference_params(self.mods, gas, wave=False)
        self.Blend = C.blend_class(self.mods)
        self.wave_arg = C._wave_arg(self.mods["model"].Model)
        self.y1 = self.R["y1"]
        self.LEDGER = self.mods["model"].LEDGER

    def model(self, e_v_cm, f_hz):
        p = C.afterglow_params(self.mods, self.p_full, e_v_cm, f_hz)
        return self.Blend(p, self.tabs, y0=self.y1)

    def rhs_args(self, wave_on: bool):
        return (False, wave_on) if self.wave_arg else (wave_on,)


def band_energy(M, y, band):
    """Energy [eV/cm3] emitted in band, from a state's event counters."""
    return sum(y[M.i_ev + c.j] * M.N0 * c.route["rad"] for c in M.compiled
               if c.route["rad"] and c.row.band == band)


class Slab:
    def __init__(self, S: Setup, f_hz, i0_w_cm2, cells, back="vacuum"):
        self.S = S
        self.e_inc = math.sqrt(i0_w_cm2 * 1e4 * Z0) / 100.0  # rms V/cm
        self.models = [S.model(self.e_inc, f_hz) for _ in cells]
        self.m0 = self.models[0]
        self.dx = np.asarray(cells, float)
        self.nc = len(cells)
        self.size = self.m0.size
        self.omega = 2.0 * math.pi * f_hz
        self.k0 = self.omega / C0
        self.I0 = i0_w_cm2 * 1e4                              # W/m2
        self.back = back
        self.args = S.rhs_args(True)
        self.y0 = np.tile(S.y1, self.nc)
        for m in self.models:
            m._kinetic_update(0.0, S.y1, *self.args)

    def wave(self, y):
        m = self.m0
        Y = y.reshape(self.nc, self.size)
        n = np.maximum(Y[:, :m.ns], 0.0) * m.N0
        ne = np.array([m.electron_density(row) for row in n])
        te = np.maximum(Y[:, m.i_ue] * m.N0 / (1.5 * ne), 1e-3)
        nu = np.empty(self.nc)
        for c in range(self.nc):
            tab = self.models[c].blend(ne[c])
            lt = min(max(math.log(te[c]), tab["lnte"][0]), tab["lnte"][-1])
            nu[c] = C.nu_eff(np.interp(lt, tab["lnte"], tab["A"]))
        sigma = ne * 1e6 * QE ** 2 / (ME * (nu + 1j * self.omega))
        eps = 1.0 - 1j * sigma / (self.omega * EPS0)
        S, absorbed, R = tmm(eps, self.dx, self.k0, self.back)
        return absorbed * self.I0 / self.dx * 1e-6 / QE, R, S[-1]   # eV/cm3/s

    def rhs(self, t, y):
        p_em, _, _ = self.wave(y)
        out = np.empty_like(y)
        Y = y.reshape(self.nc, self.size)
        for c, m in enumerate(self.models):
            m.p_em_cell = p_em[c]
            out[c * self.size:(c + 1) * self.size] = m.rhs(t, Y[c], *self.args)
        return out

    def run(self, t_end=C.AFTERGLOW_S, n_out=201, rtol=1e-6):
        m = self.m0
        atol = np.tile(np.where(np.arange(self.size) == m.i_tg, 1e-6, 1e-12), self.nc)
        blk = csr_matrix(np.ones((self.size, self.size)))
        sol = solve_ivp(self.rhs, (0.0, t_end), self.y0, method="BDF",
                        t_eval=np.linspace(0.0, t_end, n_out), rtol=rtol,
                        atol=atol, jac_sparsity=block_diag([blk] * self.nc, format="csr"),
                        first_step=1e-15)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol


def zero_d(S: Setup, e_v_cm, f_hz):
    M = S.model(e_v_cm, f_hz)
    res = M.run()
    return M, res


def pump_numbers(S: Setup):
    full = S.R["full"]
    k1 = K.at(full, C.PULSE_S)
    M = S.model(0.0, C.F_BUILD_HZ)
    y = np.zeros(M.size)
    for name in ("events",):
        pass
    bands = {}
    for c in M.compiled:
        if c.route["rad"]:
            bands[c.row.band] = bands.get(c.row.band, 0.0) + \
                float(full["events"][c.row.rid][k1]) * c.route["rad"]
    return float(full["ledger"]["E_dep"][k1]), bands          # eV/cm3


def twin_numbers(given, pump_band, off_band, on_band, wave_abs, wave_inc=None):
    """Efficiencies of the main band (fractions) and the pair's numbers;
    energies in any one unit (per cm3 or per cm2)."""
    eff_off = (pump_band + off_band) / given
    eff_on = (pump_band + on_band) / (given + wave_abs)
    d = on_band - off_band
    out = {"eff_off": eff_off, "eff_on": eff_on, "d_eff_pp": (eff_on - eff_off) * 100,
           "marginal": d / wave_abs if wave_abs > 0 else float("nan"),
           "source_breakeven": (eff_off / (d / wave_abs)) if wave_abs > 0 and d > 0 else float("inf")}
    if wave_inc:
        out["eff_on_incident"] = (pump_band + on_band) / (given + wave_inc)
        out["d_eff_incident_pp"] = (out["eff_on_incident"] - eff_off) * 100
    return out


def run_zero_d(gas):
    """LIGHT's kinetic afterglow twins against the kit's closure, 0-D."""
    S = Setup(gas)
    band = MAIN_BAND[gas]
    given, pump = pump_numbers(S)
    M_off, r_off = zero_d(S, 0.0, C.F_BUILD_HZ)
    M_on, r_on = zero_d(S, C.E_RMS_V_CM, C.F_BUILD_HZ)
    i_em = S.LEDGER.index("E_em")

    def light(rec):
        return float(sum(rec["events"][c.row.rid][-1] * c.route["rad"]
                         for c in M_off.compiled if c.route["rad"] and c.row.band == band))
    out = {"gas": gas, "band": band, "E_given_eV_cm3": given, "pump_band_eV_cm3": pump.get(band, 0.0)}
    lk_off, lk_on = light(S.R["ag_off"]), light(S.R["ag_on"])
    lk_w = float(S.R["ag_on"]["ledger"]["E_em"][-1])
    out["light_kinetic"] = {"off_band": lk_off, "on_band": lk_on, "wave_abs": lk_w,
                            **twin_numbers(given, pump.get(band, 0.0), lk_off, lk_on, lk_w)}
    bk_off = float(sum(r_off.events[c.row.rid][-1] * c.route["rad"] for c in M_off.compiled
                       if c.route["rad"] and c.row.band == band))
    bk_on = float(sum(r_on.events[c.row.rid][-1] * c.route["rad"] for c in M_on.compiled
                      if c.route["rad"] and c.row.band == band))
    bk_w = float(r_on.ledger["E_em"][-1])
    out["kit_closure"] = {"off_band": bk_off, "on_band": bk_on, "wave_abs": bk_w,
                          **twin_numbers(given, pump.get(band, 0.0), bk_off, bk_on, bk_w)}
    return out


def run_slab(gas, f_thz, i_mw, L_cm=10.0, ncell=30, back="vacuum"):
    S = Setup(gas)
    band = MAIN_BAND[gas]
    given, pump = pump_numbers(S)
    t0 = time.time()
    M_off, r_off = zero_d(S, 0.0, f_thz * 1e12)
    off_band = float(sum(r_off.events[c.row.rid][-1] * c.route["rad"] for c in M_off.compiled
                         if c.route["rad"] and c.row.band == band))
    cells = cells_for(L_cm * 1e-2, ncell)
    sl = Slab(S, f_thz * 1e12, i_mw * 1e6, cells, back)
    sol = sl.run()
    m = sl.m0
    dx_cm = cells * 100.0
    Yend = sol.y[:, -1].reshape(sl.nc, sl.size)
    i_em = m.i_led + S.LEDGER.index("E_em")
    on_band = float(sum(band_energy(m, Yend[c], band) * dx_cm[c] for c in range(sl.nc)))
    wave_abs = float((Yend[:, i_em] * m.N0 * dx_cm).sum())
    Rs, Ts = [], []
    for k in range(sol.y.shape[1]):
        _, R, T = sl.wave(sol.y[:, k])
        Rs.append(R)
        Ts.append(T)
    I0_eV = i_mw * 1e6 / QE                                    # eV/cm2/s
    inc = I0_eV * C.AFTERGLOW_S
    out = {"gas": gas, "band": band, "f_THz": f_thz, "I0_MW_cm2": i_mw, "L_cm": L_cm,
           "cells": ncell, "back": back, "per_cm2_eV": True,
           "E_given": given * L_cm, "pump_band": pump.get(band, 0.0) * L_cm,
           "off_band": off_band * L_cm, "on_band": on_band, "wave_abs": wave_abs,
           "wave_inc": inc, "wave_refl": float(np.trapezoid(np.array(Rs), sol.t)) * I0_eV,
           "wave_trans": float(np.trapezoid(np.array(Ts), sol.t)) * I0_eV,
           "x_cm": (np.cumsum(dx_cm) - dx_cm / 2).tolist(),
           "absorbed_eV_cm3": (Yend[:, i_em] * m.N0).tolist(),
           "tg_end_K": Yend[:, m.i_tg].tolist(),
           "seconds": time.time() - t0}
    out.update(twin_numbers(out["E_given"], out["pump_band"], out["off_band"], on_band,
                            wave_abs, inc))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gas", choices=("KrF", "Kr"))
    ap.add_argument("f_THz", type=float, nargs="?", default=1.3)
    ap.add_argument("I0_MW_cm2", type=float, nargs="?", default=0.3)
    ap.add_argument("--L", type=float, default=10.0)
    ap.add_argument("--cells", type=int, default=30)
    ap.add_argument("--back", choices=("vacuum", "pec"), default="vacuum")
    ap.add_argument("--zero-d", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.zero_d:
        res = run_zero_d(a.gas)
        name = f"{a.gas}_zero_d.json"
    else:
        res = run_slab(a.gas, a.f_THz, a.I0_MW_cm2, a.L, a.cells, a.back)
        name = f"{a.gas}_f{a.f_THz:g}_I{a.I0_MW_cm2:g}_L{a.L:g}_{a.back}.json"
    (OUT / name).write_text(json.dumps(res, indent=1))
    keys = [k for k in ("eff_off", "eff_on", "d_eff_pp", "marginal", "eff_on_incident",
                        "d_eff_incident_pp") if k in res]
    print(name, {k: res[k] for k in keys} if keys else
          {k: {kk: vv for kk, vv in v.items() if kk in ("eff_off", "eff_on", "d_eff_pp", "marginal", "wave_abs")}
           for k, v in res.items() if isinstance(v, dict) and "eff_off" in v})
