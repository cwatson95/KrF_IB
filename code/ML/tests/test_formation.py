"""The formation model: LIGHT's chain reproduced, the probe observable
checked against the two-level identity, and the chain's bookkeeping."""
import pytest

from ml import formation as F
from ml import forward, score


@pytest.fixture(scope="module")
def fit_set():
    return score.candidate_sets()["fit (train set)"]


def test_default_chain_is_light(fit_set):
    """With LIGHT's values the extra rows carry nothing: the run is
    run_device's (to the solver's 1e-9)."""
    dev = forward.run_device(fit_set)
    r = F.run(fit_set)
    assert r["g"] == pytest.approx(dev["g0"], rel=1e-6)
    assert r["formation"] == pytest.approx(dev["formation"], rel=1e-6)


def test_probe_books_the_extracted_power(fit_set):
    """The power the probe takes out, g(I) I, is the stimulated events
    times h nu (the units of sigma_B, I and h nu agree). The 0.4 % left is
    the window average: events are least-squares slopes of the cumulative
    count, the gain is the mean of KrF(B), and the medium drifts."""
    r = F.run(fit_set, I_MW=2.0)
    extracted = r["g"] / 100.0 * 2.0e6                       # W/cm3
    stim = r["events"]["STIM"] * r["pairs"] * F.HV_PROBE * forward.QE_SI
    assert stim == pytest.approx(extracted, rel=1e-2)


def test_lasing_cools_the_gas(fit_set):
    """Light that leaves as laser output is not heat: the gas runs cooler
    under a 4 MW/cm2 flux, which is why Isat grows slowly with I."""
    assert F.run(fit_set, I_MW=4.0)["tg"] < F.run(fit_set)["tg"] - 5.0


def test_branching_keeps_the_total_harpoon_rate(fit_set):
    """gamma_H and f_C_H split Kr* + F2 without changing how fast Kr*
    reacts, and the run's events carry the shares."""
    ch = F.LIGHT_CHAIN.with_(gamma_H=0.7, f_C_H=0.3)
    k = {cc.row.rid: cc.rate for cc in F.build(fit_set, ch).compiled}
    assert k["G3"] + k["G3C"] + k["G3X"] == pytest.approx(8.47e-10)
    ev = F.run(fit_set, ch)["events"]
    krf = ev["G3"] + ev["G3C"]
    assert krf / (krf + ev["G3X"]) == pytest.approx(0.7, rel=1e-6)
    assert ev["G3C"] / krf == pytest.approx(0.3, rel=1e-6)


def test_isat_is_a_lifetime_not_a_formation_number(fit_set):
    """Changing only how much KrF forms (the harpoon yield from 0.90 to
    0.63) leaves Isat within 3 %."""
    a = F.observables(fit_set)
    b = F.observables(fit_set, F.LIGHT_CHAIN.with_(gamma_H=0.63))
    assert b["isat"] == pytest.approx(a["isat"], rel=0.03)
    assert b["eta"] < 0.9 * a["eta"]


def test_measured_chain_regression(fit_set):
    """Every measured link at its measured value, which since step 2 of
    2026-09-28 is LIGHT's own chain, with LIGHT's electron closure (step 3)
    and the absolute harpoon yield 0.90 (Velazco 1977): formation x1.05 of
    Tahoma with the fitted B/C set, g0 x0.96, Isat x1.09. With the yield at
    1 they were x1.10, x1.00 and x1.10, and on the Maxwellian closure
    x1.13, x1.07 and x1.06."""
    assert F.LIGHT_CHAIN == F.LIGHT_CHAIN.with_(**__import__(
        "ml.formation_data", fromlist=["MEASURED"]).MEASURED)
    from ml import formation_data as D
    o = F.observables(fit_set, F.LIGHT_CHAIN.with_(**D.MEASURED))
    assert o["eta"] == pytest.approx(0.284, abs=0.003)
    assert o["g0"] == pytest.approx(6.13, abs=0.05)
    assert o["isat"] == pytest.approx(1.53, abs=0.02)


def test_kr2_route_is_not_a_loss_with_its_measured_lifetime(fit_set):
    """With the measured 270 ns Kr2* lifetime most Kr2* reacts before it
    radiates; REBUILD's 12.5 ns made it look like a 4 % formation loss."""
    ch = F.LIGHT_CHAIN.with_(k3=1.7e-32, k3_tpow=-0.64, k_kr2_f2=2.1e-10)
    slow = F.run(fit_set, ch.with_(kr2_rad=3.7e6))
    fast = F.run(fit_set, ch.with_(kr2_rad=8.0e7))
    assert slow["events"]["KR2_RAD"] < 0.15 * slow["events"]["AK8_Kr"]
    assert fast["events"]["KR2_RAD"] > 0.6 * fast["events"]["AK8_Kr"]
    assert slow["per_pair"] > fast["per_pair"]


def test_formation_registry_is_complete():
    """Every measured value names a Chain link, a source, a place in it and
    how it was read."""
    from dataclasses import fields
    from ml import formation_data as D
    links = {f.name for f in fields(F.Chain)}
    for v in D.VALUES + D.ION_NOTES + D.ANALOGS + D.ABSENT:
        assert v.link in links, v
        assert v.value > 0 or v.link.endswith("tpow"), v
        assert v.err >= 0 and v.source and v.where
        assert v.status in ("checked", "secondhand", "abstract"), v
