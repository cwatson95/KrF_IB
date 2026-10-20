"""MAXWELL verification tier (0.1.37): closed-form gates on the integrator,
bucket parity with the LOKI tier, the fixed-mean-energy two-Maxwellian,
wrapper delegation, and the engine path (beam drive runs and closes;
the circuit lane is refused by capability)."""

import numpy as np
import pytest

from REBUILD.MODULES.loki_tier import (RUN_NAMES_AR_FILL, RUN_NAMES_NE_FILL,
                                       RUN_NAMES_ARF_CHEN, LokiRates)
from REBUILD.MODULES.maxwell_tier import (EPS_GRID, MaxwellRates,
                                          MaxwellShaped, bulk_temperature,
                                          maxwell_pdf, maxwell_rate)
from REBUILD.MODULES.units import EV, M_E

_trapz = getattr(np, "trapezoid", None) or np.trapz

AR_FILL = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
NE_FILL = {"Ne": 0.949, "Kr": 0.05, "F2": 0.001}
CHEN = {"Ne": 0.964, "Ar": 0.035, "F2": 0.001}


def _vbar(T):
    return np.sqrt(8.0 * EV * T / (np.pi * M_E))


def test_pdf_is_normalised_with_mean_three_halves_T():
    f = maxwell_pdf(EPS_GRID, 2.0)
    assert abs(_trapz(f, EPS_GRID) - 1.0) < 1e-6
    assert abs(_trapz(EPS_GRID * f, EPS_GRID) / 3.0 - 1.0) < 1e-6


@pytest.mark.parametrize("T", [0.05, 0.5, 3.0, 20.0])
def test_constant_sigma_closed_form(T):
    sigma = np.full(EPS_GRID.size, 2.5e-20)
    k = maxwell_rate(EPS_GRID, sigma, T)
    assert abs(k / (2.5e-20 * _vbar(T)) - 1.0) < 1e-4


@pytest.mark.parametrize("T,thr", [(1.0, 3.0), (3.0, 11.5), (8.0, 16.0)])
def test_step_sigma_closed_form(T, thr):
    sigma = np.where(EPS_GRID >= thr, 1.0e-20, 0.0)
    k = maxwell_rate(EPS_GRID, sigma, T)
    ref = 1.0e-20 * _vbar(T) * (1.0 + thr / T) * np.exp(-thr / T)
    assert abs(k / ref - 1.0) < 3e-3          # the step sits between grid points


@pytest.mark.parametrize("mixture,run_names", [
    (AR_FILL, RUN_NAMES_AR_FILL), (NE_FILL, RUN_NAMES_NE_FILL),
    (CHEN, RUN_NAMES_ARF_CHEN)])
def test_every_loki_bucket_has_a_maxwell_counterpart(mixture, run_names):
    lok = LokiRates(run_names)
    mx = MaxwellRates(mixture)
    missing = set(lok.available_buckets()) - mx.available_buckets()
    assert not missing, f"LOKI buckets without a sigma-set counterpart: {missing}"


def test_two_maxwellian_reduces_and_conserves_mean_energy():
    mx = MaxwellRates(AR_FILL)
    te = np.array([0.5, 1.0, 2.0])
    assert np.allclose(mx.k_two_t("Ar*", te, 0.0, 5.0), mx.k("Ar*", te))
    a, tt = 0.05, 5.0
    tb = bulk_temperature(te, a, tt)
    assert np.allclose((1.0 - a) * tb + a * tt, te)
    with pytest.raises(ValueError):
        bulk_temperature(0.2, 0.05, 5.0)          # tail exceeds the mean energy
    # a stated tail at fixed mean energy raises the threshold rates
    assert mx.k_two_t("Ar*", 1.0, a, tt) > mx.k("Ar*", 1.0)
    assert mx.k_two_t("Ar+", 1.0, a, tt) > mx.k("Ar+", 1.0)


def test_rates_are_monotone_in_te_and_vectorised():
    mx = MaxwellRates(AR_FILL)
    te = np.array([0.3, 1.0, 3.0, 10.0])
    for b in ("Ar*", "Ar**", "Ar+", "Kr*", "Kr+", "F2(a3)"):
        k = mx.k(b, te)
        assert k.shape == te.shape and np.all(np.diff(k) > 0.0), b
        assert k[1] == pytest.approx(mx.k(b, 1.0))


def test_maxwell_vs_loki_same_order_at_the_table_top():
    """Units/mapping smoke gate: at the highest tabulated <eps> the two
    shapes agree to well within a decade on ionization (loose on purpose —
    the shapes differ; a x100 miss would be a units or mapping defect)."""
    lok = LokiRates(RUN_NAMES_AR_FILL)
    mx = MaxwellRates(AR_FILL)
    xs, ys = lok.curves["Ar+"]
    eps_top = float(np.exp(xs[-1]))
    ratio = mx.k("Ar+", eps_top / 1.5) / float(np.exp(ys[-1]))
    assert 0.05 < ratio < 20.0, ratio


def test_shaped_wrapper_delegates_and_refuses_en_lookups():
    lok = LokiRates(RUN_NAMES_AR_FILL)
    tier = MaxwellShaped(lok, AR_FILL)
    assert tier.available_buckets() == lok.available_buckets()
    assert tier.drift_velocity(50.0) == lok.drift_velocity(50.0)
    assert tier.curves is lok.curves
    assert not tier.fallback_buckets
    assert tier.k("Ar+", 2.0) == pytest.approx(tier.mx.k("Ar+", 2.0))
    with pytest.raises(NotImplementedError):
        tier.k_by_en("Ar+", 50.0)
    tail = MaxwellShaped(lok, AR_FILL, tail=(0.05, 5.0))
    assert tail.k("Ar*", 1.0) == pytest.approx(tier.mx.k_two_t("Ar*", 1.0, 0.05, 5.0))
    assert "two-Maxwellian" in tail.describe()


def test_shaped_wrapper_follows_the_loki_rule_below_the_table_floor():
    """In-table: the bare Maxwellian. Below the floor: flat for 0.35
    ln-units of <eps>, then the exponential taper — continuous at the
    floor, the same rule LokiRates.k applies."""
    lok = LokiRates(RUN_NAMES_AR_FILL)
    tier = MaxwellShaped(lok, AR_FILL)
    b = "Ar+"
    f = tier.floor_te(b)
    assert f is not None and f > 0.0
    assert tier.k(b, 2.0 * f) == pytest.approx(tier.mx.k(b, 2.0 * f))
    assert tier.k(b, f) == pytest.approx(tier.mx.k(b, f))
    # flat just below the floor (within 0.35 ln-units of <eps>)
    assert tier.k(b, f * np.exp(-0.2)) == pytest.approx(tier.mx.k(b, f))
    # tapering further down, exactly as the LOKI rule
    te = f * np.exp(-1.0)
    dx = (np.log(1.5 * f) - 0.35) - np.log(1.5 * te)
    assert tier.k(b, te) == pytest.approx(tier.mx.k(b, f) * np.exp(-dx / 0.35))
    # vectorised path agrees with the scalar path across the floor
    tes = np.array([te, f, 3.0 * f])
    assert np.allclose(tier.k(b, tes), [tier.k(b, x) for x in tes])


def _beam_params(tier, mode="krf", t_end=20.0):
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": mode, "drive": "beam", "t_end_ns": t_end,
                "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141, "Ly_m": 0.0141,
                "Lz_m": 0.75, "beam_axis": 0}, source="problem")
    p.set("electron_tier", tier, source="cli")
    return p


def test_engine_runs_on_the_maxwell_tier_closes_and_is_deterministic():
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.maxwell_tier import MaxwellShaped
    dev = Device(_beam_params("maxwell"))
    assert isinstance(dev.tiers, MaxwellShaped), type(dev.tiers)
    res = dev.run(n_save=9)
    assert res["success"]
    assert max(abs(c["residual_rel"]) for c in res["closure"][1:]) < 5e-4
    again = Device(_beam_params("maxwell")).run(n_save=9)
    assert np.array_equal(res["y_final"], again["y_final"])
    ref = Device(_beam_params("loki")).run(n_save=9)
    # the tiers differ (that is the point) but both are finite and positive
    assert np.isfinite(res["ne_m3"]).all() and res["ne_m3"][-1] > 0.0
    assert res["ne_m3"][-1] != pytest.approx(ref["ne_m3"][-1], rel=1e-6)


def test_engine_refuses_the_maxwell_tier_on_the_circuit_lane():
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "arf_chen", "drive": "discharge",
                "discharge_mode": "circuit", "t_end_ns": 5.0,
                "nx": 1, "ny": 1, "nz": 1}, source="problem")
    p.set("electron_tier", "maxwell", source="cli")
    with pytest.raises(ValueError, match="MAXWELL"):
        Device(p)
