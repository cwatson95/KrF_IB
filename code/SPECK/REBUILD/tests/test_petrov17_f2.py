"""Gates for the petrov17 beam-EEDF e+F2 channels ([0.1.36], flag-gated).

The contract: OFF the engine is bit-inert (default); ON, the marked
low-threshold e+F2 rows (dissociative attachment, the two repulsive
dissociation rows) and the vibrational sink evaluate Petrov 2017 JAP
122 133301 Table I coefficients k(P, x_F2) WHILE the beam envelope is
on, and revert to the swarm rates when it is off (one owner per
electron population). The formulas are pinned against hand-derived
literals from the printed table — not against this module's own code.
"""

import numpy as np

from REBUILD.MODULES.chemistry import petrov17_f2_k
from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters


def _params(petrov17: bool, **over):
    p = engine_parameters()
    base = {"mode": "krf", "drive": "beam", "pressure_atm": 2.0,
            "t_end_ns": 30.0, "beam_kV": 500.0, "beam_j_A_cm2": 140.0,
            "beam_flat_ns": 10.0, "beam_rise_ns": 2.0, "beam_axis": 1,
            "nx": 1, "ny": 1, "nz": 1,
            "Lx_m": 0.1, "Ly_m": 0.1, "Lz_m": 0.1,
            "petrov17_f2": petrov17}
    base.update(over)
    p.set_many(base, source="problem")
    return p


def test_formula_anchors_hand_derived():
    # Table I at P = 1, x_F2 = 0 (all literals derived by hand from the
    # printed formulas; cm^3/s -> m^3/s is 1e-6):
    #   vib: a + b = 1e-9 + (2 - 0.45 ln 1)e-9 = 3.0e-9 cm^3/s
    #   att: b = (2.3 - 0.5 ln 1)e-9 = 2.3e-9 cm^3/s
    #   dis: a/(1 + 1/b) with a = 2.4e-9, b = 0.35
    #        -> 2.4e-9/3.857142857... = 6.2222e-10 cm^3/s
    assert np.isclose(petrov17_f2_k("vib", 1.0, 0.0), 3.0e-15,
                      rtol=1e-12)
    assert np.isclose(petrov17_f2_k("att", 1.0, 0.0), 2.3e-15,
                      rtol=1e-12)
    assert np.isclose(petrov17_f2_k("dis", 1.0, 0.0),
                      2.4e-15 / (1.0 + 1.0 / 0.35), rtol=1e-12)
    # P = 1.5, x_F2 = 0.003 (the Electra class), attachment by hand:
    #   b = (2.3 - 0.5 ln 1.5)e-9 = 2.0973e-9; c = 50/7 = 7.1429
    #   k = b/(1 + c sqrt(0.003)) = b/1.39124 = 1.5075e-9 cm^3/s
    b = (2.3 - 0.5 * np.log(1.5)) * 1.0e-9
    k_hand = b / (1.0 + (50.0 / 7.0) * np.sqrt(0.003)) * 1e-6
    assert np.isclose(petrov17_f2_k("att", 1.5, 0.003), k_hand,
                      rtol=1e-12)


def test_formula_limits_and_shapes():
    # finite at the paper's own P floor; F2 loading suppresses att/vib.
    # Dissociation is NON-monotone by the printed form (the (b+cx)^-1
    # suppression weakens before the d*sqrt(x) term wins), so it gets a
    # positivity check only.
    for ch in ("vib", "att", "dis"):
        assert np.isfinite(petrov17_f2_k(ch, 1e-3, 0.0))
        assert petrov17_f2_k(ch, 1.0, 0.01) > 0.0
    for ch in ("vib", "att"):
        assert petrov17_f2_k(ch, 1.0, 0.01) < petrov17_f2_k(ch, 1.0, 0.0)
    # vectorized over cells
    P = np.array([0.5, 1.0, 2.0])
    x = np.array([0.002, 0.003, 0.005])
    out = petrov17_f2_k("dis", P, x)
    assert out.shape == (3,) and np.all(out > 0.0)


def test_marked_rows_discovered_in_beam_modes():
    dev = Device(_params(True), mixture={"Ar": 0.948, "Kr": 0.05,
                                         "F2": 0.002})
    comp = dev.compiled
    assert comp._petrov17
    assert comp._p17_att, "no dissociative-attachment row marked"
    assert set(comp._p17_dis.values()) == {"F2(a3)", "F2(A1)"}
    # OFF: nothing marked, flag false
    dev0 = Device(_params(False), mixture={"Ar": 0.948, "Kr": 0.05,
                                           "F2": 0.002})
    assert not dev0.compiled._petrov17
    assert not dev0.compiled._p17_att and not dev0.compiled._p17_dis


def test_off_is_inert_and_on_swaps_only_under_beam():
    mix = {"Ar": 0.948, "Kr": 0.05, "F2": 0.002}
    dev0 = Device(_params(False), mixture=mix)
    dev1 = Device(_params(True), mixture=mix)
    y0 = dev0.initial_state()
    # mid-beam (t = 5 ns, envelope on): the swap must CHANGE the rhs
    d_on0 = dev0.rhs(5e-9, y0)
    d_on1 = dev1.rhs(5e-9, y0.copy())
    assert not np.allclose(d_on0, d_on1), \
        "petrov17 ON changed nothing under the beam"
    # beam off (t = 25 ns > rise+flat+fall): identical evaluations —
    # the swarm rates own the afterglow on both paths
    d_off0 = dev0.rhs(25e-9, y0)
    d_off1 = dev1.rhs(25e-9, y0.copy())
    assert np.array_equal(d_off0, d_off1), \
        "petrov17 leaked outside the beam envelope"


def test_ledger_closes_with_flag_on():
    dev = Device(_params(True, t_end_ns=25.0),
                 mixture={"Ar": 0.948, "Kr": 0.05, "F2": 0.002})
    res = dev.run(n_save=11)
    assert res["success"]
    worst = max(abs(c["residual_rel"]) for c in res["closure"])
    assert worst < 1e-5, f"ledger residual {worst:.2e} with petrov17 on"
