"""gates for tools/paschen.py -- Townsend counting + the Raizer transcription.

Fast gates only (the full curve is a tool run, cached to picmcc_out): the
similarity-law transcription reproduces the book's own printed argon anchor,
the avalanche counter conserves every electron it made, counting is
deterministic and monotone in V, and the two-length alpha sits in a stated
band of the Table 4.1 fit (the x0.8 deficit is a recorded finding, not noise).
"""

import numpy as np

from PIC_MCC_3D.tools import paschen as P


def test_raizer_transcription_reproduces_the_printed_anchor():
    # Eq. (7.3)+(7.4) with A=12, B=180, gamma_iron=1.5e-3 must land on the
    # printed Ar/iron minimum (265 V at 1.5 Torr cm) -- checks OUR reading
    # of Table 4.1 / p.135, nothing about the PIC.
    v = float(P.raizer_vt(1.5, P.GAMMA_IRON))
    assert abs(v - 265.0) < 4.0
    pd = np.geomspace(0.6, 8.0, 60)
    vt = P.raizer_vt(pd, P.GAMMA_IRON)
    i_min = int(np.nanargmin(vt))
    assert 1.2 < pd[i_min] < 1.9          # (pd)_min ~ 1.5


def test_avalanche_counter_ledger_and_determinism():
    a = P.measure_M0(1.5, 180.0, seed=3)
    b = P.measure_M0(1.5, 180.0, seed=3)
    assert a["M0"] == b["M0"]
    assert a["M0"] > 0.5                     # multiplication happened
    assert not a["step_cap_hit"]
    # every electron that ever existed is at a wall now: seeds + births
    assert np.isclose(a["back_absorbed"] + a["anode_absorbed"],
                      a["n_seed"] + a["n_ion"], rtol=1e-12)
    # nobody sampled past the sigma-table ceiling energy
    assert a["V"] < a["eps_ceiling_eV"]


def test_M0_monotone_in_V():
    lo = P.measure_M0(1.5, 170.0, n_seed=32, seed=9)
    hi = P.measure_M0(1.5, 230.0, n_seed=32, seed=9)
    assert hi["M0"] > 2.0 * lo["M0"]


def test_two_length_alpha_in_stated_band():
    r = P.alpha_two_length(150.0, p_torr=1.5, d1_cm=0.5, d2_cm=0.9,
                           n_seed=48, seed=17)
    # measured 2026-08-14: x0.80 of Raizer (4.5) at 176 V/(cm Torr); the
    # band says the deficit stays bounded, not that it is zero
    ratio = r["alpha_p_pic"] / r["alpha_p_45"]
    assert 0.5 < ratio < 1.2
    assert r["M0_long"] > r["M0_short"]


def test_breakdown_bisection_brackets():
    # a huge gamma makes the target M0 tiny -> V_b must come out low and
    # bracketed, with a monotone-consistent trace
    res = P.breakdown_voltage(1.5, gamma=0.5, v_lo=80.0, v_max=400.0,
                              n_bisect=5, n_seed=24, seed=21)
    assert res["bracket"] == "ok"
    assert 80.0 < res["V_b"] < 400.0
