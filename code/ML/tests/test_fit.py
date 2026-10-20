"""The fitter recovers known constants, and the consistency check finds
the dataset that disagrees."""

import dataclasses

import numpy as np
import pytest

from ml import data, fit, forward, params


def _synthetic(truth):
    out = []
    for d in data.by_role("train"):
        pred = forward.predict(d, truth)
        pts = tuple((cond, float(v)) for (cond, _), v in zip(d.points, pred))
        out.append(dataclasses.replace(d, points=pts))
    return out


def test_fit_recovers_known_constants():
    truth = params.consts(np.log([2e-11, 1.5, 55e-9, 6.9e-9, 1e-12, 2e-12,
                                  8e-31, 4e-31, 4.5e-10, 3e-10]))
    x = fit.fit(_synthetic(truth))
    got = params.consts(x)
    for n in ("k_BC", "K_eq", "tau_C", "k_F2", "k_F2_C"):
        assert getattr(got, n) == pytest.approx(getattr(truth, n), rel=0.05), n


def test_consistency_singles_out_hakuta():
    groups = {"Eden": data.by_name("eden_kr", "eden_f2"),
              "Kansas State": data.by_name("ksu_fig1", "ksu_fig2"),
              "Hakuta": data.by_name("hakuta_fig2", "hakuta_fig3"),
              "Yu & Setser": data.by_name("ys_ir300")}
    rows = {r["left_out"]: r for r in fit.consistency(groups)}
    without = rows["Hakuta"]["chi2"]
    assert max(without["Eden"], without["Kansas State"],
               without["Yu & Setser"]) < 0.5
    assert without["Hakuta"] > 5.0
    for other in ("Eden", "Kansas State", "Yu & Setser"):
        assert max(rows[other]["chi2"].values()) > 5.0, other
