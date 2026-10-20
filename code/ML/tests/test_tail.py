"""The two-temperature tail (ml/tail.py): LIGHT's table against the kinetic
EEDF it stands for."""
import numpy as np
import pytest

from ml import score, tail
from v0 import xsec

TAHOMA_X = 2.8 / 665.0


@pytest.fixture(scope="module")
def light_now():
    return score.candidate_sets()["LIGHT now"]


@pytest.mark.parametrize("pump", [3.3e5, 3.3e4])
def test_table_gives_the_kinetic_tail_at_tahoma(light_now, pump):
    """Tahoma's state (665 Torr, 2.8 Torr F2) is not one of the table's
    states. At its pump and a tenth of it the table's r is the kinetic one
    within 10 %, and the Kr(5s) rate it gives is the kinetic rate within
    15 %."""
    o = tail.point(light_now, TAHOMA_X, pump, 665.0)
    r = xsec.tail_ratio_table()(o["ne_N"], o["x_win"])
    assert r == pytest.approx(o["r"], rel=0.10)
    assert tail.rows_at(o["beam"], o["tg"], r)["G1"] == pytest.approx(
        1.0, rel=0.15)


def test_matching_the_tail_energy_matches_kr_5s(light_now):
    """The criterion matches the energy the five Kr rows take, and Kr(5s)
    carries it: at Tahoma G1 is then the kinetic rate within 3 %."""
    o = tail.point(light_now, TAHOMA_X, 3.3e5, 665.0)
    assert o["rows"]["G1"] == pytest.approx(1.0, abs=0.03)
    assert 0.1 < o["r"] < 0.5


def test_generator_reads_the_table_as_light_does():
    """tail.lookup and LIGHT's TailRatio agree on the pinned table."""
    t = xsec.tail_ratio_table()
    assert np.allclose(t.x, tail.GRID_X) and np.allclose(t.f2, tail.X_F2)
    for ne_n, f2 in ((3e-5, 0.004), (2e-7, 0.0), (1e-6, 0.0002),
                     (5e-9, 0.05), (1e-2, 0.01)):
        assert tail.lookup(t.r, ne_n, f2) == pytest.approx(t(ne_n, f2),
                                                           rel=1e-12)
