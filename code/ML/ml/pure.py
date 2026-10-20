"""Verification on pure rare gases: the half of the KrF chain without F2.

KrF formation runs through three segments:
  1. beam energy -> Kr+ and Kr* (the energy per ion pair W and the
     excitations per ion N_ex/N_i),
  2. Kr+ and Kr* -> Kr2* (ion conversion, recombination, three-body
     excimer formation, VUV emission),
  3. the F2 steps (harpooning, attachment, ion-ion recombination).
Pure gases test the first two alone, on their own measurements
(ml/pure_data.py).
  supply(gas)   segment 1: REBUILD's Spencer-Fano (LIGHT's beam drive
                takes its Kr yields from it) against the measured W and
                N_ex/N_i. What formation uses is their combination, the
                Kr* + Kr+ made per eV, (1 + N_ex/N_i) / W.
  kr_rates()    segment 2, link by link: LIGHT's Kr2+ recombination and Kr+
                conversion against pure-Kr afterglow measurements. v0_1
                sets both from these measurements, so there they are the
                inputs and no longer a check.
  kr_vuv(...)   segments 1 and 2 together: LIGHT with no F2 under the same
                beam drive. Every Kr* and every recombined ion should end as
                one Kr2* photon (Saito 2003: 'one scintillation photon is
                emitted from one recombination process'). The light is
                counted to the end of the afterglow, as a VUV yield
                measurement counts it. At a low pump rate this is Morozov's
                experiment (0.42 +- 0.05 for Kr).
"""
from __future__ import annotations

import dataclasses
import os
import sys

import numpy as np

from . import ROOT, forward, light
from . import pure_data as P

HV_KR2 = P.HV_VUV["Kr"]


def _rebuild():
    speck = os.environ.get("LIGHT_SPECK_DIR", str(ROOT / "SPECK"))
    if speck not in sys.path:
        sys.path.insert(0, speck)
    from REBUILD.MODULES import beam
    return beam


def supply(gas, n_m3=2.45e25):
    """Segment 1 for one gas: Spencer-Fano and measured W, N_ex/N_i and
    Kr* + Kr+ per eV (with Kannari's ledger for comparison)."""
    B = _rebuild()
    y = B.spencer_fano_yields({gas: n_m3})
    nex_sf = sum(v for s, v in y.per_pair.items()
                 if s.startswith(gas) and s.endswith("*"))
    out = {"W_sf": y.W_eV, "nex_sf": nex_sf,
           "per_eV_sf": (1.0 + nex_sf) / y.W_eV,
           "W_kannari": B.W_EV.get(gas), "nex_kannari": B.F_EXC.get(gas)}
    if out["W_kannari"]:
        out["per_eV_kannari"] = (1.0 + out["nex_kannari"]) / out["W_kannari"]
    try:
        w = P.value("W_" + gas)
        out["W_meas"] = w
        nex = P.value("nex_" + gas)
        out["nex_meas"] = nex
        out["per_eV_meas"] = (1.0 + nex) / w
        out["per_eV_ws"] = 1.0 / P.value("ws_" + gas)
    except IndexError:
        pass
    return out


def kr_rates(te_eV=(0.026, 0.3, 1.0, 2.0)):
    """Segment 2 links: LIGHT's R11p and K15, read from its rows, against
    Shiu & Biondi and Smith et al. v0 carries Kannari's R11p (1.39e-6 at
    300 K on Te^-1/2) and K15 2.5e-31. v0_1 carries Shiu & Biondi's R11p
    and a K15 of 2.4e-31 that averages Smith et al. with two others."""
    network = light("network")
    own = {r.rid: r for r in network.build(("all",), "basic")}
    kind, *arg = own["R11p"].rate
    rows = []
    for te in te_eV:
        if kind == "dr_sqrt":
            k = arg[0] / np.sqrt(max(te, 0.02))
        elif kind == "te_pow":
            k = arg[0] * (max(te, 0.02) / arg[1]) ** arg[2]
        else:
            raise ValueError(f"R11p: rate kind {kind!r}")
        meas = P.KR_DR.value * (0.025852 / te) ** P.KR_DR_POW
        rows.append((te, k, meas))
    beta = P.KR_CONV.value
    conv = (beta / (3.22e16) ** 2, beta / (3.536e16) ** 2)   # 300 K, 273 K
    k15_kind, k15 = own["K15"].rate
    assert k15_kind == "const", own["K15"].rate
    return rows, k15, conv


def kr_vuv(power_W_cm3=3.3e2, torr=665.0, pulse=450e-9, tail=20e-6,
           k3=1.7e-32, tpow=-0.64, kr2_rad=3.7e6, thermal=True, rates=()):
    """LIGHT on pure Kr (basic network + ions, kr2 and electron groups) under
    the prototype beam drive. Returns the VUV yield (Kr2* photons x 8.43 eV
    over the deposited energy, counted to `tail` after the pulse) and the
    per-pair budget. rates: band values for a version that carries bands
    (v0_2's Params.rates, e.g. (("kr2_pen", 5e-10),))."""
    model = light("model")
    Params, QE = model.Params, model.QE
    E = forward._drive()
    basic = light("network").BASIC_GROUPS
    p = Params(network="basic",
               groups=tuple(g for g in ("ions", "kr2", "electron")
                            if g in basic),
               pressure_atm=torr / 760.0, x_f2=0.0,
               power_W_cm3=power_W_cm3, pulse_s=pulse,
               t_end_s=pulse + tail, ne_mode="free", ne_cm3=1e10,
               n_save=2001, **({"rates": tuple(rates)} if rates else {}))
    m = E.BeamModel(p)
    by = {cc.row.rid: cc for cc in m.compiled}
    by["AK8_Kr"].rate = k3
    by["AK8_Kr"].row = dataclasses.replace(by["AK8_Kr"].row, tg_pow=tpow)
    by["KR2_RAD"].rate = kr2_rad
    if not thermal:
        m.compiled = [cc for cc in m.compiled if cc.row.rid
                      not in ("G1", "C1a", "C1b", "C1c", "C2")]
    r = m.run()
    pairs = p.power_W_cm3 / QE / E.W_EFF * pulse
    tot = {cc.row.rid: float(r.events[cc.row.rid][-1]) / pairs
           for cc in m.compiled}
    thermal_kr = sum(tot.get(k, 0.0) for k in ("G1", "C1a", "C1b", "C1c"))
    return {"vuv": tot["KR2_RAD"] * HV_KR2 / E.W_EFF, "per_pair": tot,
            "thermal": thermal_kr, "W_eff": E.W_EFF,
            "ceiling": (E.Y_ION + E.Y_PP + E.Y_S) * HV_KR2 / E.W_EFF}


def ideal_vuv(gas="Kr"):
    """One photon per excitation and per recombined ion, with the measured
    W and N_ex/N_i: (1 + N_ex/N_i) h nu / W."""
    s = supply(gas)
    return s["per_eV_meas"] * P.HV_VUV[gas]
