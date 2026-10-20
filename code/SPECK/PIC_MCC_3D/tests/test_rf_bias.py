"""gate for the C3/RF slice of M4: a phi(t) electrode waveform drives real
dynamics.

A single electron in a vacuum gap under phi_hi(t) = V0 cos(w t) feels the
uniform field E_z = -V0 cos(w t)/L and, starting from rest, follows exactly

    v(t) = u sin(w t),  z(t) = z0 + (u/w)(1 - cos(w t)),  u = e V0/(m L w)

(the cos drive has no secular drift term; a sin drive from rest carries the
physical ponderomotive-free DC drift u*t and walks out of any small gap).
The leapfrog must land on that closed form to O(dt^2) over several RF
cycles.  This is the piece that makes an RF CCP drivable at all (the Turner
G5.1 benchmark needs it); the external circuit and self-bias remain open M4
work and are NOT claimed here.
"""

import numpy as np

from PIC_MCC_3D.MODULES import species3d as SP
from PIC_MCC_3D.MODULES.es_sim import ESSim
from PIC_MCC_3D.MODULES.units import E_CHARGE, M_E


def test_rf_gap_single_electron_matches_closed_form():
    V0 = 10.0
    omega = 2.0 * np.pi * 500.0e6          # 500 MHz
    nz = 33
    h = 1.0e-4
    L = (nz - 1) * h                       # 3.2 mm gap
    e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=8)
    z0 = 0.5 * L
    e_pop.spawn(np.array([2e-4]), np.array([2e-4]), np.array([z0]),
                np.zeros(1), np.zeros(1), np.zeros(1), np.ones(1))
    period = 2.0 * np.pi / omega
    dt = period / 400.0
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop]), dt,
                phi_lo=0.0, phi_hi=lambda t: V0 * np.cos(omega * t), seed=2)
    u = E_CHARGE * V0 / (M_E * L * omega)
    amp = u / omega                        # ~5.6e-5 m << L
    sim.prepare()
    worst = 0.0
    for _ in range(int(3 * 400)):
        sim.step()
        t = sim.time
        z_exact = z0 + amp * (1.0 - np.cos(omega * t))
        worst = max(worst, abs(e_pop.z[0] - z_exact))
    assert e_pop.n == 1                     # never left the gap
    assert worst < 0.01 * amp


def test_rf_ccp_case_table_and_build():
    """tools/rf_ccp: the transcribed Turner case table is shaped right and
    a scaled build runs a few RF steps with collisions live (G5.1 itself is
    NOT claimed -- see the tool's docstring)."""
    from PIC_MCC_3D.tools import rf_ccp

    cases = rf_ccp.turner_cases()
    assert sorted(cases) == [1, 2, 3, 4]
    for c in cases.values():
        assert c["gas"] == "He" and c["gap_m"] == 6.7e-2
        assert c["f_Hz"] == 13.56e6
    assert cases[1]["V0"] == 450.0 and cases[4]["n_neutral_m3"] == 321.0e20
    b = rf_ccp.build_rf_gap(1, nz=33, ppc=4, seed=8)
    assert b["built"]["scaled"]
    b["sim"].run(20)
    assert b["e"].n > 0 and b["ion"].n > 0


def test_rf_bias_number_and_callable_agree_frozen():
    """A callable frozen at a constant equals the plain-number bias."""
    def build(bias):
        p = SP.ParticleSpecies(SP.Species.electron(), capacity=8)
        p.spawn(np.array([1e-4]), np.array([1e-4]), np.array([1.6e-3]),
                np.zeros(1), np.zeros(1), np.array([1e5]), np.ones(1))
        s = ESSim(4, 4, 33, 1e-4, SP.SpeciesSet([p]), 1e-12,
                  phi_lo=0.0, phi_hi=bias, seed=4)
        s.run(200)
        return p.z[0] if p.n else np.nan
    za = build(25.0)
    zb = build(lambda t: 25.0)
    assert np.isclose(za, zb, rtol=0, atol=0)
