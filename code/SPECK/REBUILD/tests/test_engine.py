"""engine gates: the composed device — 0-D, 3-D, both modes, honest lines."""

import numpy as np
import pytest

from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters


def test_krf_beam_0d_physics_band(krf_result):
    """The v2-class operating point: ne ~1e15 cm^-3, Te sub-eV — and,
    since the 0.1.9/0.1.10 manifold completions, g0 ~0.004 /cm: at
    0.42 MW/cm3 deposition the completed network is SUB-THRESHOLD.
    History of this pin: 0.0385 (0.1.3) -> 0.0278 (0.1.8, +K34Ar) ->
    0.0053 (0.1.9, the Ar heteronuclear manifold) -> 0.0036 (0.1.10:
    the ArKr+ ion manifold + Ar-fill completion) -> 0.0034 (0.1.11:
    the Chen-overlap dedup — 0.1.10's number ran 8 channels
    double-counted) -> 0.0054 with Te 1.9 eV (0.1.12: the prompt-floor
    beam split banks the Platzman share in the electron gas — hotter
    bulk, more attachment/ion-ion formation; <eps> now sits INSIDE the
    measured relax window). Sitting far below the sibling engines
    (spectre 0.0318, phantom 0.0302) is the expected signature of
    channels they omit — engine-vs-engine context lives on the
    scorecard."""
    r = krf_result
    assert r["success"]
    ne_pk = r["ne_m3"].max()
    assert 1e20 < ne_pk < 1e22
    assert 0.2 < r["Te_eV"].max() < 2.0
    g0_cm = r["g0_m1"].max() / 100.0
    assert 0.002 < g0_cm < 0.10
    assert 0.0 <= r["eta_intrinsic"] < 0.25


def test_3d_uniform_limit_matches_0d(krf_result):
    """A thin transverse cell deposits near-uniformly, so a small 3-D grid
    must reproduce the 0-D answer — the composition gate."""
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": 120.0, "nx": 2, "ny": 1, "nz": 3,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    r3 = Device(p).run(n_save=31)
    assert r3["success"]
    assert abs(r3["eta_intrinsic"] / krf_result["eta_intrinsic"] - 1) < 0.05
    worst = max(abs(c["residual_rel"]) for c in r3["closure"][1:])
    assert worst < 5e-4


def test_3d_gradient_case_runs_and_closes():
    """10 cm beam path at 3 atm: a real deposition gradient; the ledger
    must still close and the entry side must outproduce the exit side."""
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": 100.0, "nx": 6, "ny": 1, "nz": 4,
                "Lx_m": 0.30, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    dev = Device(p)
    frac = dev._dep_frac.reshape(dev.shape).sum(axis=(1, 2))
    assert frac[-1] > frac[0] * 1.05      # Bragg rise along 30 cm at 3 atm
    r = dev.run(n_save=21)
    assert r["success"]
    assert max(abs(c["residual_rel"]) for c in r["closure"][1:]) < 5e-4


def test_regime_line_states_scope(krf_result):
    line = krf_result["regime"]
    for token in ("grid=", "reactions", "ledger residual", "wall",
                  "not a PIC grid"):
        assert token in line


def test_geometry_volume_consistency_enforced():
    from REBUILD.MODULES.laser import Cavity
    p = engine_parameters()
    p.set_many({"mode": "krf", "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75},
               source="problem")
    bad = Cavity(l_c_m=0.75, l_g_m=0.60, area_m2=2e-4,
                 sigma_se_m2=5.354e-21)
    with pytest.raises(ValueError, match="gain volume"):
        Device(p, cavity=bad)


def test_w_route_recorded_and_never_defaulted(krf_device):
    y = krf_device.beam_yields
    assert y.mode in ("spencer_fano", "w_table")
    assert y.cite
    assert abs(y.identity_residual_eV()) < 1e-9


def test_mixture_validation_refuses_typos_and_nonmixtures():
    """The PIC_MCC_3D GasMixture rule, closed here 2026-08-04: a typo'd
    gas name used to vanish silently from the fill."""
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "t_end_ns": 20.0,
                "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141,
                "Ly_m": 0.0141, "Lz_m": 0.75, "beam_axis": 0},
               source="problem")
    with pytest.raises(KeyError, match="ground-state"):
        Device(p, mixture={"Argon": 0.95, "Kr": 0.05})
    with pytest.raises(KeyError, match="ground-state"):
        Device(p, mixture={"Kr*": 0.95, "F2": 0.05})
    with pytest.raises(ValueError, match="not a mixture"):
        Device(p, mixture={"Ar": 94.9, "Kr": 5.0, "F2": 0.1})


def test_beam_subtable_bridge_default_and_off_switch(krf_device, krf_result):
    """The L3 closure at the Device level: a beam device runs on the
    relaxation-bridged tier by DEFAULT (regime line + validity gauge say
    so), beam_subtable='taper' reproduces the pre-0.1.9 stand-in, and the
    two runs really differ — the measured x3.3 F- overestimate was
    load-bearing."""
    from REBUILD.MODULES.kinetic_tier import RelaxationBridged
    assert isinstance(krf_device.tiers, RelaxationBridged)
    assert krf_device.p.source("beam_subtable") == "default"
    assert "MEASURED relax curves" in krf_result["regime"]
    v = krf_result["validity"]
    assert "F-" in v["relax_bridged_buckets"]
    assert v["te_in_covered_frac"] >= v["te_in_table_frac"]
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": 120.0, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0, "beam_subtable": "taper"}, source="problem")
    dev_t = Device(p)
    assert not isinstance(dev_t.tiers, RelaxationBridged)
    r_t = dev_t.run(n_save=31)
    assert r_t["success"]
    assert "flat+taper stand-in" in r_t["regime"]
    # the bridge must matter: with F- attachment x3.3 lower at afterglow
    # <eps>, the electron density history cannot be identical
    ne_rel = abs(r_t["ne_m3"][-1] / krf_result["ne_m3"][-1] - 1.0)
    assert ne_rel > 0.01, "taper vs relax made no measurable difference"


def test_kinetic_tier_circuit_drive_refused_at_construction():
    """0.1.25: electron_tier='kinetic' + discharge_mode='circuit' was a
    CLI-reachable combination that died with AttributeError on the FIRST
    RHS call — KineticRates exports k(<eps>) only (the npz carries no
    E/N-indexed curves and no drift velocity), while the circuit path
    looks rates and conduction up by E/N. Refused at Device construction,
    by capability rather than by tier name, with the fix stated."""
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "discharge",
                "discharge_mode": "circuit", "electron_tier": "kinetic",
                "nx": 1, "ny": 1, "nz": 1}, source="problem")
    with pytest.raises(ValueError, match="electron_tier='loki'"):
        Device(p)


def test_electron_tier_selectable():
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "electron_tier": "kinetic",
                "t_end_ns": 40.0, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    dev = Device(p)
    from REBUILD.MODULES.kinetic_tier import KineticRates, RelaxationBridged
    # beam mode bridges the selected tier (the L3 closure); the ACTIVE
    # swarm source underneath must be the KINETIC export, and the bridge's
    # relax data comes from that same export (one file, one provenance)
    assert isinstance(dev.tiers, RelaxationBridged)
    assert isinstance(dev.tiers.base, KineticRates)
    assert dev.tiers.relax_source_path == dev.tiers.base.path
    r = dev.run(n_save=11)
    assert r["success"]


def test_pure_krypton_fill_lases():
    """0.1.38: a buffer-free Kr/F2 fill must reach a positive gain (Mandl,
    Klimek & Salesky 1985 measured 9.2 % efficiency at 99.6 % Kr; Suda 1984
    and Razhev 2016 ran buffer-free KrF lasers). Before 0.1.38 the lasing
    state KrF*(B) stayed at its floor (1e6 cm^-3) because nothing fed it
    from C or Bv without an Ar/Ne partner. Short beam run, 0-D."""
    from REBUILD.MODULES.laser import Cavity
    from REBUILD.MODULES.network import mode_optics
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 1.0,
                "beam_kV": 500.0, "beam_j_A_cm2": 150.0, "beam_flat_ns": 60.0,
                "beam_rise_ns": 10.0, "t_end_ns": 80.0, "nx": 1, "ny": 1,
                "nz": 1, "Lx_m": 0.30, "Ly_m": 0.30, "Lz_m": 1.0,
                "beam_axis": 0}, source="problem")
    opt = mode_optics("krf")
    cav = Cavity(l_c_m=1.0, l_g_m=1.0, K_axial=1e-5, R_OC=0.4, k_c_loss=0.0,
                 L_op_fixed=0.10, area_m2=0.09, wavelength_nm=248.0,
                 sigma_se_m2=opt["sigma_se_m2"], absorbers=dict(opt["absorbers"]))
    res = Device(p, mixture={"Kr": 0.997, "F2": 0.003}, cavity=cav).run(n_save=11)
    assert res["success"]
    g0 = np.asarray(res["g0_m1"], float)
    assert g0.max() > 0.5, f"pure Kr/F2 fill shows no gain: g0 max {g0.max():.3e} /m"
    b_peak = float(np.max(np.asarray(res["dens_avg"]["KrF*(B)"], float)))
    assert b_peak > 1e19, f"KrF*(B) never populated: {b_peak:.3e} m^-3"


def test_ledger_columns_are_structurally_dead():
    """0.1.39 teeth, the condition that broke scipy's num_jac.

    Nothing in this system reads an accumulated ledger value: the ledger
    rows are pure quadratures of the physics, so every ledger COLUMN of the
    Jacobian is exactly zero. scipy's `num_jac` cannot cope with that. It
    keeps a per-column step `factor`, multiplies it by ten every time a
    column's difference comes back "too small", and a structurally dead
    column is too small forever -- the factor runs away until
    `new_factor * y_scale` overflows against a density scale of ~1e25.
    That was the 0.1.35 stiff blocker. If this test ever fails because a
    rate learned to read a ledger, the runaway is cured but the ledger has
    stopped being a passive observer, which is a bigger problem.
    """
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 1.36,
                "t_end_ns": 40.0, "nx": 1, "ny": 1, "nz": 1},
               source="problem")
    dev = Device(p, mixture={"Ar": 0.597, "Kr": 0.40, "F2": 0.003})
    struct = dev.jac_sparsity().tocsc()
    lay = dev.lay
    for j in range(lay.i_ledger0, lay.i_ledger0 + lay.n_ledger):
        assert struct[:, j].nnz == 0, f"ledger column {j} claims a coupling"
    # and the delivered Jacobian really is zero there, finite everywhere
    y = dev.initial_state()
    J = dev.jac(0.0, y)
    A = np.asarray(J.todense()) if hasattr(J, "todense") else np.asarray(J)
    assert np.isfinite(A).all()
    assert A.shape == (lay.size, lay.size)
    led = A[:, lay.i_ledger0:lay.i_ledger0 + lay.n_ledger]
    assert not led.any(), "a ledger column carries a nonzero derivative"


def test_jacobian_matches_finite_differences():
    """The Newton matrix must be the real derivative where the derivative
    is well defined. Compared against CENTRAL differences on a mid-pulse
    state, over the columns that are not sitting on a floor: species pinned
    at the density floor sit on the `maximum(dens, 0)` and `maximum(ne, 1)`
    kinks, where a one-sided and a two-sided difference legitimately
    disagree and neither is 'the' derivative. Accuracy is governed by the
    step controller either way; this pins that the matrix is not garbage.
    """
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 1.36,
                "beam_kV": 500.0, "beam_j_A_cm2": 160.0,
                "beam_flat_ns": 100.0, "beam_rise_ns": 20.0,
                "t_end_ns": 60.0, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.30, "Ly_m": 0.30, "Lz_m": 1.00, "beam_axis": 0},
               source="problem")
    dev = Device(p, mixture={"Ar": 0.597, "Kr": 0.40, "F2": 0.003})
    res = dev.run(n_save=5, dense_cells=True)
    y = res["y_t"][:, -1]
    t = 60e-9
    J = dev.jac(t, y)
    A = np.asarray(J.todense() if hasattr(J, "todense") else J)
    lay = dev.lay
    dens = y[:lay.n_species]
    well_scaled = np.where(dens > 1e15)[0]      # off the floor by decades
    assert well_scaled.size > 8, "no well-scaled species to compare"
    devs = []
    for j in well_scaled:
        h = abs(y[j]) * 1e-6
        yp = y.copy()
        yp[j] += h
        ym = y.copy()
        ym[j] -= h
        col = (dev.rhs(t, yp) - dev.rhs(t, ym)) / (2.0 * h)
        s = max(np.abs(col).max(), 1e-30)
        devs.append(np.abs(A[:, j] - col).max() / s)
    devs = np.array(devs)
    assert np.median(devs) < 1e-4, \
        f"median column deviation {np.median(devs):.2e}"
    assert devs.max() < 5e-2, f"worst column deviation {devs.max():.2e}"


def test_jacobian_costs_no_more_than_stock_num_jac():
    """0.1.39 gate, on the integrator's own counters rather than the clock.

    A solver change is invisible to a correctness suite: the first version of
    this fix passed everything and still needed 296 steps and 50 Jacobians on
    a 24-cell case where stock scipy needs 176 and 6, because it had thrown
    away num_jac's step adaptation to stop the runaway. Wall time cannot be
    the gate either — this machine produced a 5-hour and a 75-minute reading
    of the same suite while other work held two cores. Steps, Jacobian
    evaluations and RHS calls are load-independent, so they are what is
    pinned: supplying our Jacobian must cost no more than letting scipy build
    its own, on the same problem, back to back in one process.
    """
    from scipy.integrate import solve_ivp

    def mk():
        p = engine_parameters()
        p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 1.36,
                    "beam_kV": 500.0, "beam_j_A_cm2": 160.0,
                    "beam_flat_ns": 100.0, "beam_rise_ns": 20.0,
                    "t_end_ns": 60.0, "nx": 1, "ny": 1, "nz": 1,
                    "Lx_m": 0.30, "Ly_m": 0.30, "Lz_m": 1.00,
                    "beam_axis": 0}, source="problem")
        return Device(p, mixture={"Ar": 0.597, "Kr": 0.40, "F2": 0.003})

    d = mk()
    y0 = d.initial_state()
    atol = d.atol_vector()
    rtol = d.p["rtol"]
    kw = dict(method="BDF", rtol=rtol, atol=atol, max_step=2e-9)
    mine = solve_ivp(d.rhs, (0.0, 60e-9), y0, jac=d.jac, **kw)
    stock = solve_ivp(mk().rhs, (0.0, 60e-9), y0, **kw)
    assert mine.success and stock.success
    # generous bands: the point is to catch a 3-8x blow-up, not to freeze
    # scipy's exact step sequence, which may shift with a scipy release
    assert mine.njev <= 2 * stock.njev + 3, \
        f"Jacobians {mine.njev} vs stock {stock.njev}"
    assert mine.t.size <= 1.5 * stock.t.size + 10, \
        f"steps {mine.t.size} vs stock {stock.t.size}"
    assert mine.nfev <= 1.5 * stock.nfev + 50, \
        f"rhs calls {mine.nfev} vs stock {stock.nfev}"
    # and the two must agree on the answer, not merely on the effort
    scale = np.maximum(np.abs(stock.y[:, -1]), atol)
    assert np.max(np.abs(mine.y[:, -1] - stock.y[:, -1]) / scale) < 1e-3
