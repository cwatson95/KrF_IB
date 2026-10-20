"""The cross-section atlas (0.1.38): the Ramsauer criterion, the elastic
relaxation identity, the mixture loss-channel ranking, and the sub-floor
audit that finds buckets the flat+taper rule extrapolates the wrong way.

These gate the DIAGNOSIS, not a physics change: the tool reads the shipped
LXCat sets and the shipped swarm tables and reports. Nothing it computes
enters an engine run.
"""

import numpy as np
import pytest

from REBUILD.MODULES.loki_tier import RUN_NAMES_AR_FILL, RUN_NAMES_NE_FILL
from REBUILD.MODULES.maxwell_tier import EPS_GRID
from REBUILD.MODULES.units import EV, M_E, gas_density_m3
from REBUILD.tools.sigma_atlas import (FILLS, load_set, mixture_loss_rates,
                                       ramsauer, relaxation_times, sigma_on,
                                       spread, subfloor_audit, table_floor_te)

AR_FILL = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}

#: atomic/molecular masses [amu] for the mass-ratio identity
AMU = {"Ne": 20.18, "Ar": 39.95, "Kr": 83.80, "He": 4.0026, "F2": 37.997}


@pytest.mark.parametrize("gas,expect_window", [
    ("Ar", True), ("Kr", True), ("Ne", False), ("He", False)])
def test_ramsauer_criterion_separates_the_real_windows(gas, expect_window):
    """Ar and Kr have genuine Ramsauer-Townsend minima; Ne and He do not.
    Neon's shallow low-energy dip must NOT be reported as a window: its
    zero-energy scattering length is positive, so no MERT zero exists."""
    r = ramsauer(load_set(gas)["elastic"])
    assert r["is_rt_window"] is expect_window, (gas, r)
    if expect_window:
        assert r["depth"] >= 20.0
        assert 0.1 < r["eps_min_eV"] < 1.0          # the window sits sub-eV


def test_mass_ratios_match_the_atomic_masses():
    """The LXCat mass_ratio field is m_e/M — the factor the elastic energy
    loss carries. A wrong one silently rescales every cooling rate."""
    for gas, amu in AMU.items():
        st = load_set(gas)
        ratio = st["elastic"].mass_ratio
        assert ratio == pytest.approx(1.0 / (amu * 1822.888), rel=0.02), gas


def test_relaxation_time_matches_its_closed_form():
    st = load_set("Ar")
    n = gas_density_m3(3.0, 300.0)
    rel = relaxation_times(st, n)
    i = int(np.argmin(np.abs(EPS_GRID - 3.0)))
    eps = EPS_GRID[i]
    v = np.sqrt(2.0 * eps * EV / M_E)
    nu_m = n * sigma_on(st["elastic"], eps) * v
    assert rel["nu_m"][i] == pytest.approx(nu_m, rel=1e-9)
    assert rel["nu_u"][i] == pytest.approx(2.0 * st["elastic"].mass_ratio * nu_m,
                                           rel=1e-9)
    assert rel["tau_u_s"][i] == pytest.approx(1.0 / rel["nu_u"][i], rel=1e-9)


def test_argon_cannot_cool_through_its_window_within_a_pulse():
    """The physical claim the atlas exists to make: in the Ramsauer window
    the elastic energy-relaxation time exceeds a 100 ns beam pulse by an
    order of magnitude, while neon (no window) stays under it."""
    n = gas_density_m3(3.0, 300.0)
    i = int(np.argmin(np.abs(EPS_GRID - 0.3)))
    tau_ar = relaxation_times(load_set("Ar"), n)["tau_u_s"][i]
    tau_ne = relaxation_times(load_set("Ne"), n)["tau_u_s"][i]
    assert tau_ar > 5e-7, tau_ar                    # ~1 us, x10 the pulse
    assert tau_ne < 1e-7, tau_ne                    # under the pulse
    assert tau_ar / tau_ne > 5.0


def test_trace_f2_vibration_outranks_buffer_elastic_in_the_window():
    """At 0.3 eV the buffer's elastic loss has collapsed into the Ramsauer
    minimum and 0.1% F2 vibrational excitation dominates the cooling."""
    lr = mixture_loss_rates(AR_FILL, 3.0)
    i = int(np.argmin(np.abs(lr["eps_eV"] - 0.3)))
    inelastic = sum(v[i] for v in lr["by_channel"].values())
    assert inelastic > 10.0 * lr["elastic"][i]
    top = max(lr["by_channel"].items(), key=lambda kv: kv[1][i])[0]
    assert top.startswith("F2:excitation") and "V" in top, top
    # and at 3 eV, above the window, the ordering reverses
    j = int(np.argmin(np.abs(lr["eps_eV"] - 3.0)))
    assert sum(v[j] for v in lr["by_channel"].values()) < lr["elastic"][j]


@pytest.mark.parametrize("mixture,run_names", [
    (AR_FILL, RUN_NAMES_AR_FILL),
    ({"Ne": 0.949, "Kr": 0.05, "F2": 0.001}, RUN_NAMES_NE_FILL)])
def test_subfloor_audit_finds_exactly_the_two_f2_channels(mixture, run_names):
    """The flat+taper rule always falls below the table floor. Exactly two
    buckets physically RISE there, and both are F2: dissociative attachment
    (the F- source, the recorded x13 starvation class) and vibrational
    excitation (the bulk's cooling channel). Every other bucket is
    threshold-limited and genuinely falls."""
    aud = subfloor_audit(mixture, run_names)
    rising = {b for b, r in aud["buckets"].items()
              if r["true_rate_rises_below_floor"]}
    assert rising == {"F-", "vib"}, rising
    for b in rising:
        r = aud["buckets"][b]
        assert r["maxwell_peak_te_eV"] < r["floor_te_eV"]
        assert r["maxwell_peak_over_floor"] > 1.05
        assert r["probe"]["0.3 eV"]["ratio"] > 50.0
    # a threshold bucket must NOT be flagged: its rate really does collapse
    assert not aud["buckets"]["Kr+"]["true_rate_rises_below_floor"]


def test_table_floors_are_above_the_measured_operating_band():
    """Why the audit matters: every fill's table floor sits ABOVE the
    settled beam-mode Te, so the engine runs on the extrapolation."""
    for name, (_mix, run_names, _p) in FILLS.items():
        d = table_floor_te(run_names)
        assert d["floor_te_eV"] > 0.6, name        # settled Te 0.04-0.63 eV
        assert d["en_range_td"][0] <= 1.0


def test_database_spread_is_a_symmetric_band():
    te = np.geomspace(0.1, 10.0, 6)
    rates = {"a": {"el": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]},
             "b": {"el": [2.0, 2.0, 6.0, 4.0, 5.0, 3.0]}}
    sp = spread(rates, te)
    assert sp["el"]["band"][0] == pytest.approx(2.0)
    assert sp["el"]["band"][1] == pytest.approx(1.0)
    assert sp["el"]["worst"] == pytest.approx(2.0)
    assert min(sp["el"]["band"]) >= 1.0
    assert not spread({"only": {"el": [1.0]}}, te[:1])   # one set, no band
