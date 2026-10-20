"""The whole-engine energy budget as a GATED INVARIANT (S3 — the headline).

Defect class D3 (ORESTES_CLASS_REBUILD.txt §2): phantom11 has no enforced
energy budget, and PLAN (70)'s ~18x electron-heating → light inconsistency
is the open consequence. The prevention is structural: every joule entering
the cell is accounted to a named store or sink, the flow terms are computed
INSIDE the same RHS the physics integrates (so the identity holds
analytically and any residual is integrator error), and the suite fails if
the budget does not close.

STORES (recomputed from the state, never integrated — no drift possible):
  chem_potential    sum over cells of V_c * sum_s n_s * E_form(s)
  electron_thermal  sum V_c * U_e
  gas_thermal       sum V_c * (3/2) n_tot k_B (T_g - T_g0)
  photon_field      n_ph * h nu * V_cavity

CUMULATIVE CHANNELS (integrated as state rows by the same BDF step):
  drive_in          beam energy deposited + discharge energy delivered [J]
  laser_out         extracted through the output coupler [J]
  cavity_loss       photon losses at HR/windows (the 1 - L_op ln-term) [J]
  radiated_escape   spontaneous/VUV emission not captured by the cavity [J]

THE INVARIANT, at every saved time:

  drive_in = (laser_out + cavity_loss + radiated_escape)
             + delta(chem_potential + electron_thermal + gas_thermal
                     + photon_field)
             + residual,   |residual| <= tol * max(drive_in, floor)

Anything that would break this — an unrouted exothermic release, a photon
absorbed without transforming its absorber, a beam split that over-spends W
— shows up as residual, not as a silently wrong efficiency.
"""

from __future__ import annotations

import numpy as np

from .species import e_form_of
from .units import EV, K_B

CHANNELS = ("drive_in", "laser_out", "cavity_loss", "radiated_escape")
N_CHANNELS = len(CHANNELS)


def stores_J(layout, y, cell_vol_m3, n_tot_m3, tg0_K, cavity_vol_m3,
             hv_J) -> dict:
    """The four stores, evaluated from one state vector."""
    dens = layout.density_block(y)
    cells = layout.cells_view(y)
    u_e = cells[layout.i_ue]
    tg = cells[layout.i_tg]
    e_form = np.array([e_form_of(s) for s in layout.species])
    chem = float((e_form[:, None] * dens).sum(axis=0) @ cell_vol_m3) * EV
    elec = float(u_e @ cell_vol_m3)
    gas = float((1.5 * n_tot_m3 * K_B * (tg - tg0_K)) @ cell_vol_m3)
    ph = float(y[layout.i_ph]) * hv_J * cavity_vol_m3
    return {"chem_potential": chem, "electron_thermal": elec,
            "gas_thermal": gas, "photon_field": ph}


def channels_J(layout, y) -> dict:
    off = layout.i_ledger0
    return {name: float(y[off + i]) for i, name in enumerate(CHANNELS)}


def closure_report(layout, y0, y1, cell_vol_m3, n_tot_m3, tg0_K,
                   cavity_vol_m3, hv_J) -> dict:
    """Budget between two states; the gate calls this at every saved step."""
    s0 = stores_J(layout, y0, cell_vol_m3, n_tot_m3, tg0_K, cavity_vol_m3, hv_J)
    s1 = stores_J(layout, y1, cell_vol_m3, n_tot_m3, tg0_K, cavity_vol_m3, hv_J)
    c0 = channels_J(layout, y0)
    c1 = channels_J(layout, y1)
    d_stores = {k: s1[k] - s0[k] for k in s0}
    d_chan = {k: c1[k] - c0[k] for k in c0}
    sinks = (d_chan["laser_out"] + d_chan["cavity_loss"]
             + d_chan["radiated_escape"])
    residual = d_chan["drive_in"] - sinks - sum(d_stores.values())
    # normalize by the total energy MOVEMENT, not drive alone: before the
    # pulse the drive is ~0 while seed species rearrange picojoules, and
    # dividing by ~nothing would report closure noise as failure
    scale = max(abs(d_chan["drive_in"]),
                sum(abs(v) for v in d_stores.values()) + abs(sinks), 1e-9)
    return {"stores_delta_J": d_stores, "channels_delta_J": d_chan,
            "residual_J": residual, "residual_rel": residual / scale,
            "drive_in_J": d_chan["drive_in"]}


def format_budget(rep: dict) -> str:
    d = rep["channels_delta_J"]
    s = rep["stores_delta_J"]
    lines = [f"drive in            {d['drive_in']:12.5e} J"]
    for k in ("laser_out", "cavity_loss", "radiated_escape"):
        lines.append(f"  sink {k:<14s} {d[k]:12.5e} J")
    for k, v in s.items():
        lines.append(f"  store {k:<13s} {v:+12.5e} J")
    lines.append(f"  RESIDUAL           {rep['residual_J']:+12.5e} J "
                 f"({rep['residual_rel']:+.3e} of drive)")
    return "\n".join(lines)
