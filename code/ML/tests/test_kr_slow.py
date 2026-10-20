"""ML's krypton swarm training end to end, on its cached simulations: the
emulator against the true solver, the consistency check's train set, and
the factors LIGHT_Kr pins against a fresh true-solver run. Slow (the true
solver at every measured point, and the emulator's training when its
cache is cold): run with ML_KR_SLOW=1."""
import json
import os

import numpy as np
import pytest

from ml import kr_data, kr_light, kr_nn, kr_swarm, kr_train

pytestmark = pytest.mark.skipif(not os.environ.get("ML_KR_SLOW"),
                                reason="slow: set ML_KR_SLOW=1")


@pytest.fixture(scope="module")
def s():
    return kr_train.setup(log=lambda *a: None)


def test_emulator_is_well_inside_the_datas_errors(s):
    """Held-out designs: the emulator's rms error (ln units) under 1 % for
    W and D_T/mu and 2 % for alpha/N, against data errors of 3 to 10 %."""
    err = s["emu_err"]
    assert err[0] < 0.01 and err[1] < 0.01 and err[2] < 0.02


def test_emulator_matches_the_true_solver_at_a_random_trial(s):
    rng = np.random.default_rng(42)
    z = rng.uniform(-0.5, 0.5, kr_swarm.NP) * kr_swarm.BOX
    for name in ("w_hunter88", "a_jacques86", "w_nakamura90"):
        d = next(e for e in s["ds"] if e.name == name)
        emu = s["model"](z, d)
        true = s["true"](z, d)
        assert np.max(np.abs(np.log(emu / true))) < 0.03, name


def test_the_train_set_is_the_consistency_checks(s):
    c = kr_train.run_consistency(s, log=None)
    final = c["rounds"][-1]
    assert max(final["joint"].values()) <= 3.0
    for name, chi in c["dropped"]:
        assert chi > 3.0
    # no check set and no low-field alpha/N is ever a candidate
    assert not set(c["kept"]) & set(kr_data.CHECK + kr_data.HELD_BY_RULE)


def test_the_pinned_factors_are_this_fits(s):
    xs = kr_light.module("xsec")
    rec = json.loads((xs.DATA_DIR / xs.TRAINED_FILE[0]).read_text())
    c = kr_train.run_consistency(s, log=None)
    assert rec["train"] == c["kept"] == list(kr_data.TRAIN)
    by = {d.name: d for d in s["ds"]}
    z = np.array(rec["z"])
    # the true solver, now, at every train dataset: within the threshold
    chi = kr_swarm.chi2_per_point(z, [by[n] for n in rec["train"]],
                                  s["true"], s["sig"])
    for n, v in chi.items():
        assert v <= 3.0, n
        assert v == pytest.approx(rec["chi2_per_point"][n], rel=1e-3)
