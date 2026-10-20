"""Verification on pure rare gases (ml/pure.py): the F2-free half of the
KrF chain against its own measurements."""
import pytest

from ml import pure
from ml import pure_data as P


def test_pure_registry_is_complete():
    for v in P.ALL:
        assert v.value > 0 and v.err >= 0 and v.source and v.where, v
        assert v.status in ("checked", "secondhand", "abstract"), v


@pytest.mark.parametrize("gas", ["Ar", "Kr"])
def test_beam_supply_per_ev_matches_the_measured(gas):
    """Spencer-Fano's W is 11 % high and its N_ex/N_i 40 % high, but the
    Kr* + Kr+ it makes per eV, which is what formation uses, is within 5 %
    of the measured (1 + N_ex/N_i)/W (Saito 2003 and the measured W)."""
    x = pure.supply(gas)
    assert x["W_sf"] / x["W_meas"] == pytest.approx(1.11, abs=0.04)
    assert x["per_eV_sf"] == pytest.approx(x["per_eV_meas"], rel=0.05)
    assert x["per_eV_kannari"] < 0.88 * x["per_eV_meas"]


def test_kr_ion_chemistry_matches_pure_kr_afterglows():
    rows, k15, conv = pure.kr_rates(te_eV=(0.3, 1.0, 2.0))
    for te, light, meas in rows:
        assert light == pytest.approx(meas, rel=0.15), te
    assert 0.75 * conv[1] < k15 < 1.25 * conv[0]


@pytest.fixture(scope="module")
def low_pump():
    return pure.kr_vuv(3.3e2)


def test_light_pure_kr_vuv_regression(low_pump):
    """LIGHT at 0.33 kW/cm3 on its two-temperature closure (2026-09-28,
    step 3): VUV yield 0.570, still above its own lossless ceiling (0.553)
    because thermal electrons excite Kr, 0.060 per pair. The Maxwellian
    closure gave 0.595 and 0.14. The kinetic loop (ml/eedf.py) puts the
    Kr(5s) excitation at about half of LIGHT's, and it is not zero: with no
    F2 to slow them, e-e lifts some injected electrons over the break. At
    this density the bulk is not Maxwellian (LIGHT's Te 0.88 eV, the
    kinetic 0.78-0.82), which no tail closure can represent."""
    assert low_pump["vuv"] == pytest.approx(0.570, abs=0.01)
    assert low_pump["vuv"] > low_pump["ceiling"]
    assert 0.03 < low_pump["thermal"] < 0.09


@pytest.mark.xfail(strict=True, reason="LIGHT's pure-Kr VUV yield is x1.4 "
                   "Morozov 2008's 0.42 +- 0.05 (ML/CHANGELOG 2026-09-28 b)")
def test_light_pure_kr_vuv_matches_morozov(low_pump):
    assert abs(low_pump["vuv"] - 0.42) < 2 * 0.05
