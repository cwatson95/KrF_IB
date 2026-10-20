"""circuit.py + the E/N-driven discharge drive — gates.

The drive that replaces the prescribed-power stand-in: Chen's C0/Cp/L
chain as engine state, conduction from the swarm drift velocity, and
electron-impact rates looked up by E/N DIRECTLY (the tables' native index
— no <eps> = (3/2)Te mapping on this path at all)."""

import numpy as np
import pytest

from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters


def _circuit_params(**over):
    p = engine_parameters()
    d = {"mode": "arf_chen", "drive": "discharge",
         "discharge_mode": "circuit", "pressure_atm": 3.4,
         "t_end_ns": 200.0, "nx": 1, "ny": 1, "nz": 1,
         "Lx_m": 0.012, "Ly_m": 0.010, "Lz_m": 0.90}
    d.update(over)
    p.set_many(d, source="problem")
    return p


def test_lc_energy_conserved_without_plasma():
    """R = 0, Rf = 0, no preionization, 1 kV (E/N ~1 Td — every rate
    tapered to nothing): the C0-Lm-Cp-Ld-gap chain PLUS the (lossless at
    Rf = 0) Lf-Cf filter branch must conserve its own energy. The
    analytic-oscillator gate for the circuit alone."""
    p = _circuit_params(v0_kV=1.0, r_loop_ohm=0.0, rf_ohm=0.0,
                        preion_m3=0.0, t_end_ns=150.0)
    res = Device(p).run(n_save=31)
    e = res["circuit"]["circuit_energy_J"]
    assert e[0] > 0.0
    assert np.abs(e / e[0] - 1.0).max() < 1e-3
    # and essentially nothing entered the (empty) plasma
    assert res["final_channels_J"]["drive_in"] < 1e-4 * e[0]


def test_filter_branch_damps_and_only_dissipates():
    """0.1.27: the paper's RLC filter branch (Rf 70.71 / Lf 100 nH /
    Cf 10 pF, Sec. II A) is modeled. With R = 0 and no plasma, Rf is the
    ONLY loss: circuit energy must be non-increasing (successive saves,
    up to integrator tolerance) and must measurably fall — the branch is
    alive, and it can only dissipate, never source. Without breakdown
    the branch is only GRAZED (quasi-static: i_f ~ Cf dV/dt ~ 0.1 A ->
    ~2e-5 of the bank over 150 ns, measured); its real work is
    post-breakdown ringing, exercised by the base-case run below."""
    p = _circuit_params(v0_kV=1.0, r_loop_ohm=0.0, preion_m3=0.0,
                        t_end_ns=150.0)
    res = Device(p).run(n_save=31)
    e = res["circuit"]["circuit_energy_J"]
    assert e[0] > 0.0
    assert float(np.diff(e).max()) < 1e-6 * e[0]   # never sources
    assert e[-1] < (1.0 - 3e-6) * e[0]             # alive: it damps


def test_drift_velocity_pin(loki_arf):
    """BulkV_z at 15 Td from the packaged ChenArF swarm table."""
    assert abs(loki_arf.drift_velocity(15.0) / 5.422e4 - 1.0) < 0.01
    assert loki_arf.drift_velocity(0.0) == 0.0
    # linear-to-zero below the table, monotone in the table
    assert loki_arf.drift_velocity(2.5) < loki_arf.drift_velocity(5.0)


def test_en_direct_vs_eps_mapped_consistency(loki_arf):
    """The same table indexed two ways (E/N direct vs <eps>-mapped Te)
    must agree — interpolation differences only, well under x2."""
    for b in ("Ar+", "Ne*", "F-"):
        for en in (10.0, 20.0, 50.0):
            k_en = loki_arf.k_by_en(b, en)
            te = loki_arf.mean_energy_at_en(en) / 1.5
            k_te = float(loki_arf.k(b, te))
            if k_en > 0 and k_te > 0:
                assert 0.5 < k_en / k_te < 2.0, (b, en)


def test_chen_default_circuit_is_the_paper_base_case():
    """0.1.25 regression pin: the defaults are Chen Sec II A / Table II
    (10 nF / 250 nH / 10 nF / 5 nH, 18 kV; E_in = C0 V0^2/2 = 1.62 J as
    the paper states). The 0.1.5..0.1.24 defaults were the Fig. 5
    VALIDATION experiment's components (10.4 / 8.1 / 30 nH, '5 -> 1 uH
    switch') grafted onto base-case conditions — the Known-bad class."""
    p = engine_parameters()
    assert p["c0_nF"] == 10.0 and p["cp_nF"] == 10.0
    assert p["lm_uH"] == 0.25 and p["ld_nH"] == 5.0
    assert p["v0_kV"] == 18.0
    e_in = 0.5 * p["c0_nF"] * 1e-9 * (p["v0_kV"] * 1e3) ** 2
    assert abs(e_in - 1.62) < 1e-12


def test_chen_base_case_class():
    """The Chen 2025 base case through the real circuit (0.1.25: the
    paper's Sec II A constants + CHEN_ERRATA): self-consistent deposition
    ~1 J class from the 1.62 J bank, E/N never exceeds the table, the
    breakdown avalanche runs, lasing ignites, ledger closes. The
    BELOW-table fraction is REPORTED (the gauge was one-sided through
    0.1.24). Residual misses vs the paper are the scorecard's business
    (chen_figure6_rows), not asserted away here."""
    res = Device(_circuit_params()).run(n_save=41)
    assert res["success"]
    worst = max(abs(c["residual_rel"]) for c in res["closure"][1:])
    assert worst < 5e-4
    cz = res["circuit"]
    assert 0.5 < res["final_channels_J"]["drive_in"] < 1.62
    assert cz["en_above_table_frac"] == 0.0
    assert 0.0 < cz["en_below_table_frac"] < 1.0   # afterglow reported
    # the Chen set's floor is 0.5 Td as of 0.1.27 (ChenArF_NeArF2_lowEN
    # takes the breakdown ramp off the flat+taper stand-in; only the
    # deep afterglow below 0.5 Td remains on it)
    assert abs(cz["en_table_range_td"][0] - 0.5) < 1e-9
    assert cz["en_below_table_frac"] < 0.5         # the ramp is IN-table now
    assert 8.0 < cz["en_td"].max() < 30.0
    assert 1e20 < res["ne_m3"].max() < 1e21        # x0.79-class avalanche
    assert res["n_ph"].max() > 1e17                # lasing ignites (0.1.25)
    # the drive_in channel equals the circuit's losses to the plasma:
    # stored - (left in circuit) - (R dissipation) = deposited; R > 0
    # makes this an inequality
    e0 = cz["circuit_energy_J"][0]
    assert res["final_channels_J"]["drive_in"] < e0


def test_nonfinite_rhs_guard_names_the_state():
    """The 0.1.25 non-finite guard: an inf/NaN derivative raises a
    diagnostic that NAMES the state instead of surfacing minutes later as
    an inscrutable LinAlgError in the BDF lu_factor. (Found running the
    AS-PRINTED Chen table under the paper's own circuit — that lane's
    non-integrability is a recorded scorecard measurement, 2026-08-25,
    deliberately not re-run here: it grinds the stiff solver for minutes
    on its way to the blow-up. This test pins the MECHANISM cheaply.)"""
    import numpy as np
    dev = Device(_circuit_params())
    y = dev.initial_state()
    y[dev.lay.i_ue] = np.nan                     # cell 0's U_e
    with pytest.raises(FloatingPointError, match="non-finite RHS"):
        with np.errstate(all="ignore"):
            dev.rhs(0.0, y)


def test_circuit_mode_is_0d_only():
    with pytest.raises(ValueError, match="0-D"):
        Device(_circuit_params(nx=2))


def test_prescribed_mode_unchanged():
    """The circuit is opt-in: discharge_mode='prescribed' keeps the old
    Gaussian path bit-for-bit (no aux states, no E/N lookup)."""
    p = _circuit_params()
    p.set("discharge_mode", "prescribed", source="api")
    d = Device(p)
    assert d.circuit is None and d.lay.n_aux == 0
