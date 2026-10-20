"""gate for tools/ebsd.py -- the beam-sustained (non-self-sustaining) gap.

The structural physics of Raizer Fig. 7.1, cheaply: with an external beam
source and NO gamma process, collected current must be finite, grow with
bias (collection -> amplification), and report a small ion-space-charge
ratio.  The full i-V artifact lives in picmcc_out/ebsd (tool run; measured
2026-08-14: 5.6 -> 49.8 A/m^2 over 20 -> 240 V, x5.8 amplification 60 -> 240 V,
consistent with exp(alpha d) at that E/p).
"""

import numpy as np

from PIC_MCC_3D.tools import ebsd


def test_beam_sustained_current_grows_with_bias_and_stays_bounded():
    lo = ebsd.run_point(30.0, nz=9, window_ns=6.0, beam_j_A_cm2=1e-5,
                        seed=6)
    hi = ebsd.run_point(150.0, nz=9, window_ns=6.0, beam_j_A_cm2=1e-5,
                        seed=6)
    assert lo["j_anode_A_m2"] > 0.0
    assert hi["j_anode_A_m2"] > 1.3 * lo["j_anode_A_m2"]
    # the smallness bound the tool reports must actually be small here
    assert hi["spacecharge_over_applied"] < 0.2
    # beam energy ledger: something entered and something stopped
    assert hi["beam_in_J"] > 0.0
    assert 0.0 < hi["beam_dep_J"] <= hi["beam_in_J"]


def test_ebsd_deterministic():
    a = ebsd.run_point(80.0, nz=9, window_ns=3.0, seed=11)
    b = ebsd.run_point(80.0, nz=9, window_ns=3.0, seed=11)
    assert a["j_anode_A_m2"] == b["j_anode_A_m2"]
    assert a["n_ionizations"] == b["n_ionizations"]
