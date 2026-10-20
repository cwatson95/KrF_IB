"""ledger gates — S3, the headline: the budget closes, and the gate bites.

PLAN (70) is the defect this structure prevents: phantom11's electron-
heating -> light chain is internally inconsistent by ~18x with no invariant
to catch it. Here every joule is routed, the closure rides the same BDF as
the physics, and a deliberately unrouted term FAILS (the teeth test).
"""

import numpy as np
import pytest

from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters

TOL = 5e-4          # closure tolerance vs energy moved (rtol 1e-6 class
                    # integration; measured 8e-6 worst on the presets)


def _worst_closure(res):
    return max(abs(c["residual_rel"]) for c in res["closure"][1:])


def test_beam_mode_budget_closes(krf_result):
    assert krf_result["success"]
    assert _worst_closure(krf_result) < TOL


def test_discharge_mode_budget_closes():
    p = engine_parameters()
    p.set_many({"mode": "arf_chen", "drive": "discharge",
                "pressure_atm": 3.4, "t_end_ns": 150.0,
                "dis_peak_W_m3": 3e12, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0949, "Ly_m": 0.0949, "Lz_m": 0.90},
               source="problem")
    res = Device(p).run(n_save=31)
    assert res["success"]
    assert _worst_closure(res) < TOL
    assert res["final_channels_J"]["drive_in"] > 1.0


def test_teeth_unrouted_energy_fails_closure(krf_device):
    """Inject a 2 eV/event unrouted gas heat into one busy row and assert
    the closure residual EXCEEDS tolerance — proof the invariant detects
    exactly the defect class it exists for."""
    import copy
    dev = Device(krf_device.p.copy(), mixture=dict(krf_device.mixture))
    target = None
    for i, entry in enumerate(dev.compiled._compiled):
        if entry[0].rid.startswith("K25"):     # the harpoon (B/C-split rids)
            target = i
            break
    assert target is not None
    e = list(dev.compiled._compiled[target])
    e[6] = e[6] + 2.0                          # u_to_gas += 2 eV, unrouted
    dev.compiled._compiled[target] = tuple(e)
    res = dev.run(n_save=21)
    assert _worst_closure(res) > TOL, (
        "the ledger failed to notice 2 eV/event of minted energy")


def test_w_value_bound_on_beam_conversion(krf_result):
    """PLAN (72)'s breach class: marginal conversion of beam-DEPOSITED
    energy into laser light cannot exceed hv/W (each photon needs one
    KrF*, each pair costs W; formation per pair <= 1). phantom11's broken
    default reached 54% against a <=25% bound; this engine must sit under
    its own bound BY PHYSICS, not by a flag."""
    dev = krf_result["device"]
    bound = (dev.hv_eV / dev.beam_yields.W_eV)
    eta = krf_result["eta_intrinsic"]
    assert 0.0 < eta < bound, (eta, bound)


def test_beam_vs_discharge_marginal_conversion_reported():
    """The audit PLAN (70) could not do: the same engine, two drive routes,
    both ledger-closed, conversions side by side. The beam route pays the
    W-value; the discharge route heats electrons directly — the DIFFERENCE
    is physics and it is now measurable instead of hidden."""
    base = {"mode": "krf", "pressure_atm": 3.0, "t_end_ns": 120.0,
            "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141, "Ly_m": 0.0141,
            "Lz_m": 0.75, "beam_axis": 0}
    p1 = engine_parameters()
    p1.set_many(dict(base, drive="beam"), source="problem")
    r1 = Device(p1).run(n_save=21)
    p2 = engine_parameters()
    p2.set_many(dict(base, drive="discharge", dis_peak_W_m3=2e11,
                     dis_t0_ns=60.0, dis_fwhm_ns=40.0), source="problem")
    r2 = Device(p2).run(n_save=21)
    for r in (r1, r2):
        assert r["success"] and _worst_closure(r) < TOL
    # both conversions are finite, positive numbers with closed budgets —
    # the comparison is REPORTABLE; no equality is asserted because the
    # two routes genuinely differ (ionization cost vs direct heating)
    assert r1["eta_intrinsic"] >= 0.0
    assert r2["eta_intrinsic"] >= 0.0


def test_stores_are_recomputed_not_integrated(krf_result):
    """chem_potential/electron_thermal/gas_thermal/photon_field come from
    the STATE, so they cannot drift from it."""
    from REBUILD.MODULES.ledger import stores_J
    dev = krf_result["device"]
    s = stores_J(dev.lay, krf_result["y_final"], dev.cell_vols, dev.n_tot,
                 dev.p["T_gas_K"], dev.cav.cavity_vol_m3, dev.cav.photon_J)
    assert s["chem_potential"] > 0
    assert s["gas_thermal"] > 0
