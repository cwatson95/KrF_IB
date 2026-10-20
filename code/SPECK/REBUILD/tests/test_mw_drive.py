"""Gates for the microwave reheat window (the beam-drive power overlay).

The window's contract: OFF it is bit-inert; ON it deposits only inside its
trapezoid, ABSORBED power never exceeds injected (the field solver caps at
mw_en_max_td and dense plasma pulls the field down — no infinite-generator
pathology), its joules enter drive_in through the same integrated ledger
as every other drive (closure is the S3 invariant, not a new bookkeeping
path), and the deposition follows e n_e v_d(E/N) E — dying with the
plasma. Cavity r_oc/l_op params keep their mode defaults on sentinels.
"""

import numpy as np
import pytest

from REBUILD.MODULES.engine import Device, _mw_envelope
from REBUILD.MODULES.params import ParamError, engine_parameters
from REBUILD.MODULES.units import E_CHARGE, M_E

MW_TEST_POWER_W = 3.0e6


def _params(**over):
    p = engine_parameters()
    base = {"mode": "krf", "drive": "beam", "pressure_atm": 2.0,
            "t_end_ns": 60.0, "beam_kV": 500.0, "beam_j_A_cm2": 140.0,
            "beam_flat_ns": 10.0, "beam_rise_ns": 2.0, "beam_axis": 1,
            "nx": 1, "ny": 1, "nz": 1,
            "Lx_m": 0.1688, "Ly_m": 0.26416, "Lz_m": 0.1688,
            "r_oc": 0.9, "l_op": 0.02,
            "mw_t0_ns": 20.0, "mw_flat_ns": 20.0, "mw_rise_ns": 2.0}
    base.update(over)
    p.set_many(base, source="problem")
    return p


@pytest.fixture(scope="module")
def mw_baseline():
    dev = Device(_params(mw_power_W=0.0))
    return dev, dev.run(n_save=121)


@pytest.fixture(scope="module")
def mw_reheat():
    dev = Device(_params(mw_power_W=MW_TEST_POWER_W))
    return dev, dev.run(n_save=121)


def test_envelope_shape():
    t0, r, flat = 20.0, 2.0, 20.0
    ns = 1e-9

    def env(t_ns):
        return _mw_envelope(t_ns * ns, t0, r, flat)

    assert env(0.0) == 0.0
    assert env(19.99) == 0.0
    assert env(21.0) == pytest.approx(0.5)       # mid-rise
    assert env(25.0) == 1.0                      # flat
    assert env(41.0) == 1.0                      # end of flat (t0+r+flat)
    assert env(43.0) == pytest.approx(0.5)       # mid-fall
    assert env(44.0) == 0.0
    assert env(100.0) == 0.0
    # zero-rise degenerate window is a clean rectangle
    assert _mw_envelope(20.5 * ns, 20.0, 0.0, 1.0) == 1.0
    assert _mw_envelope(21.5 * ns, 20.0, 0.0, 1.0) == 0.0


def test_mw_requires_beam_drive():
    p = _params(mw_power_W=1e5)
    p.set("drive", "discharge", source="cli")
    with pytest.raises(ValueError, match="beam-drive overlay"):
        Device(p)


def test_mw_params_declared_and_bounded():
    p = engine_parameters()
    assert p["mw_power_W"] == 0.0                # off by default
    with pytest.raises(ParamError):
        p.set("mw_power_W", -1.0)
    with pytest.raises(ParamError):
        p.set("mw_en_max_td", 0.0)               # below the 0.1 floor
    with pytest.raises(ParamError):
        p.set("mw_no_such_param", 1.0)


def test_cavity_design_params():
    dev = Device(_params(mw_power_W=0.0))
    assert dev.cav.R_OC == pytest.approx(0.9)
    assert dev.cav.L_op_fixed == pytest.approx(0.02)
    d_def = Device(_params(mw_power_W=0.0, r_oc=0.0, l_op=-1.0))
    assert d_def.cav.R_OC == pytest.approx(0.1)      # krf mode default
    assert d_def.cav.L_op_fixed == pytest.approx(0.088)


def test_mw_field_solver_limits():
    """The settled field: ceiling for a plasma too thin to absorb the
    feed, tiny for a dense one, and monotone between (power balance)."""
    dev = Device(_params(mw_power_W=MW_TEST_POWER_W))
    cap = float(dev.p["mw_en_max_td"])
    assert dev._mw_field_td(0.0, MW_TEST_POWER_W) == pytest.approx(cap)
    assert dev._mw_field_td(1e6, MW_TEST_POWER_W) == pytest.approx(cap)
    en_dense = dev._mw_field_td(1e24, MW_TEST_POWER_W)   # ~1e26 m^-3 0-D
    assert en_dense < 0.1
    en_mid = dev._mw_field_td(1e20, MW_TEST_POWER_W)
    assert en_dense < en_mid <= cap
    assert dev._mw_field_td(1e20, 0.0) == 0.0
    # absorbed power at the settled field reproduces the feed (in-balance
    # branch): P_abs = e (E/N) N v_d fac Ne
    en = en_mid
    v_d = dev._mw_swarm.drift_velocity(en)
    e_f = en * 1e-21 * dev.n_tot
    nu = E_CHARGE * e_f / (M_E * v_d)
    fac = nu * nu / (nu * nu + dev._mw_omega ** 2)
    p_abs = E_CHARGE * 1e20 * v_d * e_f * fac
    assert p_abs == pytest.approx(MW_TEST_POWER_W, rel=0.02)


def test_beam_delay_and_second_pulse():
    # start delay shifts the envelope exactly
    d = Device(_params(mw_power_W=0.0, beam_t0_ns=50.0))
    d0 = Device(_params(mw_power_W=0.0))
    assert d.spec.power_areal_W_m2(49e-9) == 0.0
    assert d.spec.power_areal_W_m2(55e-9) == d0.spec.power_areal_W_m2(5e-9)
    # an overlapping second pulse is refused
    with pytest.raises(ValueError, match="overlaps pulse 1"):
        Device(_params(mw_power_W=0.0, beam2_t0_ns=5.0))
    # two identical pulses deposit exactly twice the energy (deposition
    # depends on gas density only, and the density is fixed)
    p1 = _params(mw_power_W=0.0, beam_flat_ns=6.0, beam_rise_ns=2.0,
                 t_end_ns=50.0)
    r1 = Device(p1).run(n_save=51)
    p2 = _params(mw_power_W=0.0, beam_flat_ns=6.0, beam_rise_ns=2.0,
                 t_end_ns=50.0, beam2_t0_ns=25.0)
    r2 = Device(p2).run(n_save=51)
    assert r2["final_channels_J"]["drive_in"] == pytest.approx(
        2.0 * r1["final_channels_J"]["drive_in"], rel=1e-5)
    assert "2-pulse beam" in r2["regime"]


def test_rhs_finite_at_zero_electron_density():
    """The deep-afterglow guard: with every charged species at zero (the
    n_e -> 0 endpoint the undriven afterglow reaches past ~500 ns), the
    RHS must stay finite — the bare u_e/(n_e EV) 0/0 NaN'd the numerical
    Jacobian before the 1 m^-3 floor (found by the keeper runs)."""
    from REBUILD.MODULES.species import charge_of
    dev = Device(_params(mw_power_W=0.0))
    y = dev.initial_state()
    lay = dev.lay
    cells = lay.cells_view(y)
    for s, i in lay.index.items():
        if charge_of(s) != 0:
            cells[i] = 0.0
    cells[lay.i_ue] = 0.0
    dy = dev.rhs(50e-9, y)
    assert np.all(np.isfinite(dy))


def test_mw_window_beyond_t_end_is_inert(mw_baseline):
    """A scheduled-but-never-reached window must not perturb the run."""
    dev0, res0 = mw_baseline
    dev = Device(_params(mw_power_W=MW_TEST_POWER_W, mw_t0_ns=1000.0))
    res = dev.run(n_save=121)
    assert res["final_channels_J"]["drive_in"] == pytest.approx(
        res0["final_channels_J"]["drive_in"], rel=1e-10)
    assert res["final_channels_J"]["laser_out"] == pytest.approx(
        res0["final_channels_J"]["laser_out"], rel=1e-8, abs=1e-12)


def test_mw_heats_inside_window_and_ledger_closes(mw_baseline, mw_reheat):
    dev0, res0 = mw_baseline
    dev1, res1 = mw_reheat
    # the window added drive energy, and no more than injected
    e_mw = (res1["final_channels_J"]["drive_in"]
            - res0["final_channels_J"]["drive_in"])
    assert e_mw > 0.0
    e_inj_max = MW_TEST_POWER_W * (2 * 2.0 + 20.0) * 1e-9
    assert e_mw <= e_inj_max * (1.0 + 1e-9)
    # Te inside the flat top exceeds the un-driven afterglow's
    t = res1["t_ns"]
    inside = (t >= 30.0) & (t <= 40.0)
    assert np.max(res1["Te_eV"][inside]) > np.max(res0["Te_eV"][inside])
    # before the window the trajectories agree (same beam by construction)
    before = t <= 18.0
    assert np.allclose(res1["ne_m3"][before], res0["ne_m3"][before],
                       rtol=1e-6)
    # the S3 invariant holds THROUGH the window
    worst = max(abs(c["residual_rel"]) for c in res1["closure"][1:])
    assert worst < 1e-6
    # and the gauges are reported
    v = res1["validity"]
    assert v["mw_power_W"] == pytest.approx(MW_TEST_POWER_W)
    assert v["mw_nu_over_omega_at_cap"] > 3.0    # collisional regime here
    assert v["mw_cap_in_table"]
    assert "MW window" in res1["regime"]


def test_mw_energy_matches_swarm_conduction_quadrature(mw_baseline,
                                                       mw_reheat):
    """drive_in's MW share must be the integral of e n_e v_d E fac V at
    the SETTLED field — the same solver + expression the rhs uses,
    reconstructed from saved history."""
    dev, res = mw_reheat
    p = dev.p
    v_tot = float(dev.cell_vols.sum())
    t_s = res["t_ns"] * 1e-9
    p_mw = np.zeros_like(t_s)
    for j, ts in enumerate(t_s):
        env = _mw_envelope(ts, p["mw_t0_ns"], p["mw_rise_ns"],
                           p["mw_flat_ns"])
        if env <= 0.0:
            continue
        ne = res["ne_m3"][j]
        en = dev._mw_field_td(ne * v_tot, p["mw_power_W"] * env)
        v_d = dev._mw_swarm.drift_velocity(en)
        if v_d <= 0.0:
            continue
        e_f = en * 1e-21 * dev.n_tot
        nu = E_CHARGE * e_f / (M_E * v_d)
        fac = nu * nu / (nu * nu + dev._mw_omega ** 2)
        p_mw[j] = E_CHARGE * ne * v_d * e_f * fac * v_tot
    trap = getattr(np, "trapezoid", None) or np.trapz
    e_quad = float(trap(p_mw, t_s))
    # the beam share is identical by construction (same waveform, same gas
    # density) — the baseline run's ledger carries it exactly
    _dev0, res0 = mw_baseline
    e_mw_ledger = (res["final_channels_J"]["drive_in"]
                   - res0["final_channels_J"]["drive_in"])
    assert e_quad == pytest.approx(e_mw_ledger, rel=0.05)


def test_mw_window_runs_on_kinetic_tier():
    """Cross-tier check enabler: the KINETIC tier exports rates only; the
    window's swarm transport must come from the LOKI tables so the same
    window runs on either rate tier (the 0.1.23 cross-tier verification)."""
    p = _params(mw_power_W=MW_TEST_POWER_W, t_end_ns=30.0,
                mw_t0_ns=12.0, mw_flat_ns=10.0)
    p.set("electron_tier", "kinetic", source="cli")
    dev = Device(p)
    assert dev._mw_swarm is not dev.tiers          # LOKI transport object
    assert hasattr(dev._mw_swarm, "drift_velocity")
    res = dev.run(n_save=16)
    assert res["success"]
    assert "MW window" in res["regime"]


def test_mw_cell_weights_identity_masking_and_conservation():
    """[0.1.34] stated-deposition weights: w=1 per cell is bit-identical
    to the default conductivity weighting; a masked cell receives
    exactly no microwave heat; and the drive-ledger rate stays at the
    prescribed injected power either way (the settled field rescales)."""
    kw = dict(nx=2, t_end_ns=50.0)
    y = Device(_params(mw_power_W=0.0, **kw)).run(n_save=31)["y_final"]
    t_mid = 30.0e-9                       # beam off, window flat
    d_def = Device(_params(mw_power_W=MW_TEST_POWER_W, **kw))
    d_w11 = Device(_params(mw_power_W=MW_TEST_POWER_W, **kw))
    d_w11.set_mw_cell_weights([1.0, 1.0])
    d_w10 = Device(_params(mw_power_W=MW_TEST_POWER_W, **kw))
    d_w10.set_mw_cell_weights([1.0, 0.0])
    d_off = Device(_params(mw_power_W=0.0, **kw))
    f_def = d_def.rhs(t_mid, y)
    f_w11 = d_w11.rhs(t_mid, y)
    f_w10 = d_w10.rhs(t_mid, y)
    f_off = d_off.rhs(t_mid, y)
    lay = d_def.lay
    # (1) w = 1 is the default, bit for bit
    assert np.array_equal(f_w11, f_def)
    # (2) the masked cell's electron-energy derivative is exactly the
    #     no-microwave one (its deposit is exactly 0.0)
    ue_w10 = lay.cells_view(f_w10)[lay.i_ue]
    ue_off = lay.cells_view(f_off)[lay.i_ue]
    ue_def = lay.cells_view(f_def)[lay.i_ue]
    assert ue_w10[1] == ue_off[1]
    assert ue_w10[1] < ue_def[1]
    # (3) prescribed power conserved: drive-ledger rate (channel 0)
    #     equals the injected power under BOTH weightings
    drv_def = f_def[lay.i_ledger0] - f_off[lay.i_ledger0]
    drv_w10 = f_w10[lay.i_ledger0] - f_off[lay.i_ledger0]
    # each weighting settles at its own grid interval, so delivered
    # power matches the prescription at the settle grid's interpolation
    # accuracy (~1e-6 relative), for ANY weighting; the ledger books
    # the actual deposit exactly
    assert drv_def == pytest.approx(MW_TEST_POWER_W, rel=1e-4)
    assert drv_w10 == pytest.approx(MW_TEST_POWER_W, rel=1e-4)
    # (4) validation: bad shapes / negatives refused
    with pytest.raises(ValueError):
        d_def.set_mw_cell_weights([1.0])
    with pytest.raises(ValueError):
        d_def.set_mw_cell_weights([1.0, -0.5])
    with pytest.raises(ValueError):
        d_def.set_mw_cell_weights([0.0, 0.0])
