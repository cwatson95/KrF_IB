"""The measured excimer-gain anchors (0.1.40, tools/excimer_anchors.py) and
the gas-temperature rate form they motivated (Reaction.tg_floor).

These gate the MEASUREMENT APPARATUS: that the printed deposition is matched
exactly, that the scoring refuses the silent W-table fallback, that the
measured targets are the printed ones, and that the a + b Tg^n form
evaluates as documented. They do not assert any anchor ratio — anchors
score, they never become test targets (the no-calibration Known-bad).
"""
import math
from dataclasses import replace

import numpy as np
import pytest

from REBUILD.MODULES import beam as B
from REBUILD.MODULES.chemistry import CompiledNetwork
from REBUILD.MODULES.network import NetworkError, check_reaction
from REBUILD.tools import excimer_anchors as XA

KLIMEK_MIX = {"Ar": 0.948, "Kr": 0.05, "F2": 0.002}


def _klimek_dev(j=20.0, t_end=160.0):
    return XA.build("krf", KLIMEK_MIX, 1.0 * XA.AMAGAT_TO_ATM_300K, 600.0, j,
                    150.0, 5.0, t_end, XA.KLIMEK_BOX, 248.0)


def test_analytic_deposition_equals_the_ledger_drive():
    """deposited_J integrates the beam drive WITHOUT the plasma; the ledger's
    drive_in integrates the same term inside the run. They must agree to the
    quadrature error, or deposition matching would be matching the wrong
    number."""
    dev, _, _ = _klimek_dev()
    analytic = XA.deposited_J(dev, 160.0)
    res = dev.run(n_save=21)
    ledger = res["final_channels_J"]["drive_in"]
    assert analytic == pytest.approx(ledger, rel=1e-6)


def test_deposition_is_linear_in_current_and_matching_is_exact():
    """The premise of one-step matching: the drive is J x V x geometry, so
    doubling j doubles the deposition exactly; matched() then lands on the
    printed target to 1e-9 (it raises otherwise)."""
    d1, _, _ = _klimek_dev(j=20.0)
    d2, _, _ = _klimek_dev(j=40.0)
    assert XA.deposited_J(d2, 160.0) == pytest.approx(
        2.0 * XA.deposited_J(d1, 160.0), rel=1e-12)
    target = 3.21
    dev, _, _, j = XA.matched("krf", KLIMEK_MIX, 1.0, 600.0, 150.0, 5.0,
                              160.0, XA.KLIMEK_BOX, 248.0, target)
    assert XA.deposited_J(dev, 160.0) == pytest.approx(target, rel=1e-9)
    assert j > 0.0


def test_refuses_the_silent_w_table_fallback(monkeypatch):
    """2026-09-04: without LXCat data the engine silently swaps to fixed-W
    beam yields and every anchor moves x2-7. The anchor tool must refuse to
    score in that state rather than print numbers."""
    def missing(*_a, **_k):
        raise FileNotFoundError("LXCat data absent (simulated)")
    monkeypatch.setattr(B, "spencer_fano_yields", missing)
    with pytest.raises(RuntimeError, match="Spencer-Fano"):
        _klimek_dev()


def test_measured_targets_are_the_printed_ones():
    """Pins of the printed targets (Kannari 1985 Tables II/III, NRL
    MR-2023/1 Sect. 3.3, Kannari p.4314 calibration): a transcription edit
    must fail loudly, not move the score quietly."""
    assert sum(len(v) for v in XA.KLIMEK_GAIN.values()) == 18
    assert sum(len(v) for v in XA.KLIMEK_ABS.values()) == 18
    assert XA.KLIMEK_GAIN[110.0] == [9.1, 8.3, 14.0, 11.0, 15.0, 11.0]
    assert XA.KLIMEK_ABS[110.0] == [0.86, 0.73, 1.8, 1.3, 2.2, 1.6]
    # measured pressure exponent at 110 A/cm2, 0.2 % F2, 1 -> 3 amagat
    assert math.log(15.0 / 9.1) / math.log(3.0) == pytest.approx(0.4549, abs=1e-4)
    assert [g for (_p, _m, g, _w, _l) in XA.ELECTRA_GAIN] == [1.84, 1.77, 1.29, 3.21]
    assert XA.ELECTRA_TREND_MEASURED == pytest.approx(0.701, abs=1e-3)
    assert XA.KANNARI_MW_CM3_PER_ATM_A == pytest.approx(0.004, rel=1e-12)


def _flux_ratio(dev, rid, tg_hot, tg_cold=300.0):
    """row flux at tg_hot / at tg_cold, same densities, one cell."""
    lay = dev.lay
    rng = np.random.default_rng(7)
    dens = 10.0 ** rng.uniform(14, 22, size=(len(lay.species), 1))
    te = np.array([1.0])
    ne = np.array([1e20])
    hot = dev.compiled.row_fluxes(dens, te, np.array([tg_hot]), ne, dev.n_tot)
    cold = dev.compiled.row_fluxes(dens, te, np.array([tg_cold]), ne, dev.n_tot)
    return float(hot[rid][0] / cold[rid][0])


def _recompile(dev, rows):
    dev.net = rows
    dev.compiled = CompiledNetwork(rows, dev.lay, dev.tiers, dev.optics,
                                   petrov17=False)
    dev.compiled.set_photon_energy(dev.hv_eV)


def test_tg_floor_evaluates_the_normalised_sum_form():
    """k *= f + (1 - f)(Tg/300)^n: exactly 1 at 300 K, the documented value
    at 400 K, and f = 0 bit-identical to the pure power law (the pre-0.1.40
    arithmetic every existing tg_pow row still runs)."""
    dev, _, _ = _klimek_dev()
    base = {r.rid: r for r in dev.net}
    rid = "K33_B"                                    # a plain 2-body const row
    assert base[rid].tg_pow == 0.0
    f, n = 0.032, -3.9
    rows = [replace(r, tg_pow=n, tg_floor=f) if r.rid == rid else r
            for r in dev.net]
    _recompile(dev, rows)
    assert _flux_ratio(dev, rid, 300.0) == pytest.approx(1.0, rel=1e-13)
    assert _flux_ratio(dev, rid, 400.0) == pytest.approx(
        f + (1.0 - f) * (400.0 / 300.0) ** n, rel=1e-12)
    # f = 0 takes the untouched `tg_fac ** tg_pow` branch (the pre-0.1.40
    # arithmetic, so existing tg_pow rows are bit-identical); the ratio of
    # two flux products agrees with the power law to rounding
    rows0 = [replace(r, tg_pow=n, tg_floor=0.0) if r.rid == rid else r
             for r in dev.net]
    _recompile(dev, rows0)
    assert _flux_ratio(dev, rid, 400.0) == pytest.approx(
        (400.0 / 300.0) ** n, rel=1e-14)


def test_tg_floor_misuse_is_refused():
    dev, _, _ = _klimek_dev()
    r = next(r for r in dev.net if r.rid == "K33_B")
    with pytest.raises(NetworkError, match="without tg_pow"):
        check_reaction(replace(r, tg_floor=0.1))
    with pytest.raises(NetworkError, match="outside"):
        check_reaction(replace(r, tg_pow=-1.0, tg_floor=1.0))
    check_reaction(replace(r, tg_pow=-3.9, tg_floor=0.032))   # accepted
