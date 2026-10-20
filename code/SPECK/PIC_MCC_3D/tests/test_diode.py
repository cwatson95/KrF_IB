"""gate G4.1 -- the Child-Langmuir diode (tools/diode.py).

The plan's own wording: "a DC sheath in the collisionless limit follows
j = (4/9) eps0 sqrt(2q/M) V^{3/2}/d^2 within 10 %".  Run the experiment the
way the tool runs it (emit 3x more than the gap carries) and require the
transmitted current to land on the law, the virtual cathode to exist, and
the limiting to be REAL (most of the offered current turned back).
"""

import numpy as np

from PIC_MCC_3D.tools import diode


def test_child_langmuir_constant_value():
    # the classic constant: 2.334e-6 A V^-3/2 m^... x V^1.5/d^2
    assert np.isclose(diode.child_langmuir_j(250.0, 5.0e-3),
                      2.334e-6 * 250.0 ** 1.5 / 25.0e-6, rtol=1e-3)


def test_g41_transmitted_current_on_the_law():
    res = diode.run_diode(V=250.0, d_m=5.0e-3, nz=65, transits=8.0,
                          seed=20260814)
    # G4.1: within 10 % (measured 2026-08-14: 1.012 at this config; the
    # small positive residual is the stated O(sqrt(T/V)) thermal bias)
    assert abs(res["ratio"] - 1.0) < 0.10
    # the virtual cathode formed, close to the emitter
    assert res["phi_min_V"] < -0.05
    assert res["phi_min_z_m"] < 0.1 * res["d_m"]
    # limiting is real: 3x was offered, ~1x crossed
    assert res["j_meas_A_m2"] < 0.5 * res["j_emitted_A_m2"]


def test_diode_deterministic():
    a = diode.run_diode(V=150.0, d_m=3.0e-3, nz=33, transits=3.0, seed=5)
    b = diode.run_diode(V=150.0, d_m=3.0e-3, nz=33, transits=3.0, seed=5)
    assert a["j_meas_A_m2"] == b["j_meas_A_m2"]
    assert np.array_equal(a["phi_axis_V"], b["phi_axis_V"])
