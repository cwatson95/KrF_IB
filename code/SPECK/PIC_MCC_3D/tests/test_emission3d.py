"""gates for MODULES/emission3d.py (M5 partial) + the es_sim hook points.

The contract: the stated emission energy distribution is the one actually
sampled; the accumulator emits whole super-particles and loses nothing; the
flux emitter delivers exactly j*A*t/e real electrons; the gamma process emits
gamma x incident real weight, attributed to the right electrode; the bias
accepts phi(t); and none of it disturbs the old bookkeeping.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES import emission3d as EM
from PIC_MCC_3D.MODULES import species3d as SP
from PIC_MCC_3D.MODULES.es_sim import ESSim
from PIC_MCC_3D.MODULES.units import E_CHARGE, EV, M_E


def _empty_electrons(cap=4096):
    return SP.ParticleSpecies(SP.Species.electron(), capacity=cap)


def _sim(pops, nz=17, h=1e-4, dt=5e-12, phi_hi=0.0, seed=7):
    return ESSim(4, 4, nz, h, SP.SpeciesSet(list(pops)), dt,
                 phi_lo=0.0, phi_hi=phi_hi, seed=seed)


# ------------------------------------------------------- the stated sampling
def test_flux_maxwellian_moments_and_direction():
    rng = np.random.default_rng(3)
    T = 2.0
    vth = np.sqrt(2.0 * T * EV / M_E)
    vx, vy, vz = EM.flux_maxwellian(M_E, T, 200000, rng, +1.0)
    assert (vz > 0).all()
    # flux-weighted half-Maxwellian: <v_n> = vth * sqrt(pi)/2
    assert np.isclose(vz.mean(), vth * np.sqrt(np.pi) / 2.0, rtol=2e-2)
    # tangentials: plain thermals at eT/m
    assert abs(vx.mean()) < 0.01 * vth
    assert np.isclose(vy.var(), T * EV / M_E, rtol=3e-2)
    _, _, vz_dn = EM.flux_maxwellian(M_E, T, 100, rng, -1.0)
    assert (vz_dn < 0).all()


def test_accumulator_whole_particles_nothing_lost():
    acc = EM._Accumulator(w_macro=5.0)
    # 2.5/5 = 0.5 is exactly representable, so the count is float-safe
    got = [acc.take(2.5) for _ in range(10)]        # 25 real / 5 = 5 macros
    assert sum(got) == 5 and all(isinstance(g, int) for g in got)
    assert acc.emitted_macro == 5
    assert np.isclose(acc.emitted_real + acc.frac * 5.0, 25.0)


# ------------------------------------------------------------- flux emitter
def test_flux_emitter_delivers_j_a_t_over_e():
    e_pop = _empty_electrons()
    sim = _sim([e_pop], phi_hi=100.0)
    emitter = EM.FluxEmitter(e_pop, j_A_m2=3.0, wall="lo",
                             T_emit_eV=0.5, w_macro=5.0)
    sim.injectors.append(emitter)
    n_steps = 2000            # 10 ns; transit ~0.5 ns, so >90 % has crossed
    sim.run(n_steps, record_every=500)
    area = sim.L[0] * sim.L[1]
    expected_real = 3.0 * area * (n_steps * sim.dt) / E_CHARGE
    assert abs(emitter.acc.emitted_real + emitter.acc.frac * 5.0
               - expected_real) < 1e-9 * expected_real
    # conservation: emitted = collected at the anode + still in flight
    collected = sim.wall_flux_w.get(("e-", "hi"), 0.0)
    in_flight = e_pop.real_count()
    back = sim.wall_flux_w.get(("e-", "lo"), 0.0)
    assert np.isclose(collected + in_flight + back,
                      emitter.acc.emitted_real, rtol=1e-12)
    # with 100 V pulling, essentially everything crossed
    assert collected > 0.9 * emitter.acc.emitted_real
    # old global bookkeeping still runs (electrons carry -e)
    assert np.isclose(sim.wall_charge,
                      -E_CHARGE * (collected + back), rtol=1e-12)


def test_flux_emitter_deterministic():
    outs = []
    for _ in range(2):
        e_pop = _empty_electrons()
        sim = _sim([e_pop], phi_hi=50.0, seed=11)
        sim.injectors.append(EM.FluxEmitter(e_pop, 8.0, w_macro=2.0))
        sim.run(60)
        outs.append((e_pop.n, e_pop.x[:e_pop.n].copy(),
                     e_pop.vz[:e_pop.n].copy()))
    assert outs[0][0] == outs[1][0]
    assert np.array_equal(outs[0][1], outs[1][1])
    assert np.array_equal(outs[0][2], outs[1][2])


# ------------------------------------------------------------ gamma process
def test_secondary_emitter_gamma_bookkeeping():
    e_pop = _empty_electrons()
    ar_ion = SP.ParticleSpecies(SP.Species.ion("Ar", 39.948), capacity=64)
    n_i = 7
    rng = np.random.default_rng(1)
    ar_ion.spawn(x=rng.random(n_i) * 2e-4, y=rng.random(n_i) * 2e-4,
                 z=np.full(n_i, 3e-4), vx=np.zeros(n_i), vy=np.zeros(n_i),
                 vz=np.full(n_i, -2.0e4), w=np.full(n_i, 3.0))
    sim = _sim([e_pop, ar_ion])
    see = EM.SecondaryEmitter(e_pop, ar_ion, gamma=0.4, wall="lo",
                              T_emit_eV=2.0, w_macro=1.0)
    sim.wall_hooks.append(see.wall_hit)
    # ions reach the cathode at 3e-4 / 2e4 = 15 ns = step 3000; stop 0.3 ns
    # later so the 2 eV secondaries (<~3e6 m/s) are still inside the gap
    sim.run(3060)
    assert ar_ion.n == 0
    assert np.isclose(see.incident_real, 21.0)
    assert np.isclose(sim.wall_flux_w[("Ar", "lo")], 21.0)
    # gamma x 21 = 8.4 real -> 8 whole macros + 0.4 carried
    assert see.acc.emitted_macro == 8
    assert np.isclose(see.acc.frac, 0.4)
    assert e_pop.n == 8           # all 8 secondaries, none lost this fast
    assert (e_pop.vz[:8] > 0).all()
    assert (e_pop.z[:8] < 1.0e-3).all()


def test_secondary_emitter_rejects_self_watch():
    e_pop = _empty_electrons()
    with pytest.raises(ValueError):
        EM.SecondaryEmitter(e_pop, e_pop, gamma=0.1)


def test_secondary_emitter_ignores_other_wall_and_species():
    e_pop = _empty_electrons()
    ar_ion = SP.ParticleSpecies(SP.Species.ion("Ar", 39.948), capacity=64)
    ar_ion.spawn(np.array([1e-4]), np.array([1e-4]), np.array([1.3e-3]),
                 np.zeros(1), np.zeros(1), np.array([+2e4]), np.array([2.0]))
    sim = _sim([e_pop, ar_ion])
    see = EM.SecondaryEmitter(e_pop, ar_ion, gamma=0.5, wall="lo")
    sim.wall_hooks.append(see.wall_hit)
    sim.run(int(20e-9 / sim.dt))
    assert ar_ion.n == 0                      # absorbed at the ANODE
    assert see.acc.emitted_macro == 0         # cathode hook stayed silent
    assert sim.wall_flux_w[("Ar", "hi")] == 2.0


# ------------------------------------------------------------- phi(t) bias
def test_callable_bias_is_dirichlet_exact():
    e_pop = _empty_electrons()
    sim = _sim([e_pop], phi_hi=lambda t: 40.0 * np.sin(2e9 * t))
    for t_probe in (0.0, 3e-10, 7e-10):
        sim.time = t_probe
        sim.solve_fields()
        want = 40.0 * np.sin(2e9 * t_probe)
        assert np.allclose(sim.phi[:, :, -1], want, atol=1e-12)
        assert np.allclose(sim.phi[:, :, 0], 0.0, atol=1e-12)
