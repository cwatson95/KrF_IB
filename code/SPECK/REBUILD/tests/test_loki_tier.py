"""loki_tier.py gates: all packaged table sets load and behave."""

import numpy as np

from REBUILD.MODULES.loki_tier import (RUN_NAMES_AR_FILL, RUN_NAMES_ARF_CHEN,
                                       RUN_NAMES_NE_FILL, LokiRates,
                                       rates_for_mode)


def test_all_three_sets_load():
    for names, need in (
            (RUN_NAMES_ARF_CHEN, {"Ne+", "Ar+", "F-"}),
            (RUN_NAMES_NE_FILL, {"Ne+", "Kr+", "F-"}),
            (RUN_NAMES_AR_FILL, {"Ar+", "Kr+", "F-"})):
        lr = LokiRates(names)
        assert need <= lr.available_buckets(), names


def test_mixture_routing():
    assert "Ar+" in rates_for_mode(
        "krf", {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}).buckets
    ne = rates_for_mode("krf", {"Ne": 0.949, "Kr": 0.05, "F2": 0.001})
    assert "Ne+" in ne.buckets and "Ar+" not in ne.buckets


def test_ionization_rises_with_te(loki_ar_fill):
    k1 = loki_ar_fill.k("Kr+", 1.0)
    k5 = loki_ar_fill.k("Kr+", 5.0)
    assert 0.0 < k1 < k5


def test_low_end_taper_not_hard_zero(loki_ar_fill):
    """spectre's rule: a hard zero collapses the BDF step; a bare flat clamp
    cools Te without bound. Below the table the rate tapers smoothly."""
    edge = loki_ar_fill.k("Kr+", 2.0)
    low = loki_ar_fill.k("Kr+", 0.05)
    assert 0.0 < low < edge


def test_vectorized_matches_scalar(loki_ar_fill):
    te = np.array([0.3, 1.0, 4.0])
    vec = loki_ar_fill.k("F-", te)
    for i, t in enumerate(te):
        assert vec[i] == loki_ar_fill.k("F-", float(t))


def test_unknown_bucket_returns_zero(loki_ar_fill):
    assert loki_ar_fill.k("Xx+", 1.0) == 0.0
