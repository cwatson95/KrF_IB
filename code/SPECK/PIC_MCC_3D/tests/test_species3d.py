"""Gates for MISSING MODULE 2 -- mobile, multi-species kinetic particles.

The point of this module is that ions are PARTICLES, not grid charge, so the
gates check that the same movers handle both and that the mass ratio actually
shows up in the dynamics.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species, SpeciesSet,
                                          ion_subcycle_ratio, move,
                                          push_boris, push_es)
from PIC_MCC_3D.MODULES.units import AMU, E_CHARGE, EV, M_E


def test_species_records():
    e = Species.electron()
    assert e.charge < 0 and e.mass == M_E
    ar = Species.ion("Ar+", ION_MASSES_AMU["Ar+"])
    assert ar.charge > 0
    assert abs(ar.mass / AMU - 39.948) < 1e-12
    assert abs(ar.q_over_m / (E_CHARGE / ar.mass) - 1.0) < 1e-15
    f_minus = Species("F-", -E_CHARGE, ION_MASSES_AMU["F-"] * AMU)
    assert f_minus.q_over_m < 0


def test_maxwellian_load_recovers_its_temperature():
    rng = np.random.default_rng(0)
    e = Species.electron()
    T = 2.5
    p = ParticleSpecies.maxwellian_box(e, 200000, [0, 0, 0], [1e-2] * 3,
                                       1e14, T, rng)
    assert p.n == 200000
    # (2/3)<eps> over an isotropic Maxwellian is the temperature
    assert abs(p.temperature_eV() / T - 1.0) < 0.02
    assert abs(p.real_count() - 1e14 * 1e-6) / (1e14 * 1e-6) < 1e-12


def test_cold_wave_seeds_one_mode_and_is_cold():
    L, amp0 = 1e-2, 1e-4
    p = ParticleSpecies.cold_wave(Species.electron(), 16 ** 3, [0, 0, 0],
                                  [L] * 3, 1e14, amplitude=amp0,
                                  mode=1, axis=2)
    assert p.n == 16 ** 3
    assert p.kinetic_energy_J() == 0.0          # exactly cold
    # Measure the harmonics DIRECTLY from the particle positions.  A binned
    # histogram is the wrong instrument here: the quiet start is a regular
    # lattice of 16 planes, so any histogram with more bins than planes shows
    # the lattice granularity (a huge peak at the plane spacing) rather than
    # the seeded perturbation.
    k = 2.0 * np.pi / L
    z = p.z[:p.n]
    harm = [abs(complex(np.sum(np.exp(-1j * m * k * z)) / z.size))
            for m in (1, 2, 3, 4, 5)]
    assert np.argmax(harm) == 0, f"mode 1 should dominate, got {harm}"
    # and its amplitude is the linear-theory value A k / 2
    assert abs(harm[0] / (amp0 * k / 2.0) - 1.0) < 0.05


def test_cold_wave_requires_a_cube():
    with pytest.raises(ValueError, match="perfect cube"):
        ParticleSpecies.cold_wave(Species.electron(), 1000 + 1, [0, 0, 0],
                                  [1e-2] * 3, 1e14)


def test_push_es_applies_qE_over_m():
    e = Species.electron()
    p = ParticleSpecies(e, capacity=4)
    p.spawn([0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [1.0])
    E, dt = 1000.0, 1e-12
    push_es(p, np.array([E]), np.array([0.0]), np.array([0.0]), dt)
    assert abs(p.vx[0] - e.q_over_m * E * dt) < 1e-9 * abs(p.vx[0])


def test_ion_and_electron_accelerations_differ_by_the_mass_ratio():
    """The same field, the same mover: only m changes. This is the whole point."""
    e = Species.electron()
    ar = Species.ion("Ar+", ION_MASSES_AMU["Ar+"])
    E, dt = 500.0, 1e-12
    out = {}
    for sp in (e, ar):
        p = ParticleSpecies(sp, capacity=2)
        p.spawn([0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [1.0])
        push_es(p, np.array([E]), np.array([0.0]), np.array([0.0]), dt)
        out[sp.name] = p.vx[0]
    ratio = abs(out["e-"] / out["Ar+"])
    assert abs(ratio - ar.mass / M_E) / (ar.mass / M_E) < 1e-12


def test_boris_conserves_speed_in_pure_B():
    rng = np.random.default_rng(2)
    e = Species.electron()
    n = 500
    p = ParticleSpecies(e, capacity=n)
    v = rng.normal(0.0, 1e6, (n, 3))
    p.spawn(np.zeros(n), np.zeros(n), np.zeros(n), v[:, 0], v[:, 1], v[:, 2],
            np.ones(n))
    s0 = np.sqrt(p.vx[:n] ** 2 + p.vy[:n] ** 2 + p.vz[:n] ** 2)
    z = np.zeros(n)
    B = np.full(n, 0.05)
    for _ in range(20):
        push_boris(p, z, z, z, z, z, B, 1e-12)
    s1 = np.sqrt(p.vx[:n] ** 2 + p.vy[:n] ** 2 + p.vz[:n] ** 2)
    assert np.max(np.abs(s1 - s0) / s0) < 1e-13


def test_boris_rotation_angle_matches_theory():
    """One step rotates by exactly 2 atan(omega_c dt / 2), not omega_c dt."""
    e = Species.electron()
    p = ParticleSpecies(e, capacity=1)
    p.spawn([0.0], [0.0], [0.0], [1e6], [0.0], [0.0], [1.0])
    Bz, dt = 0.05, 2e-11
    z = np.zeros(1)
    push_boris(p, z, z, z, z, z, np.array([Bz]), dt)
    got = np.arctan2(p.vy[0], p.vx[0])
    wc = e.q_over_m * Bz                        # signed; negative for e-
    # SIGN: the Boris half-rotation vector is t = (q/2m) B dt, and
    # v' = v + v x t sends +x to -t_z in y.  So the in-plane rotation is
    # -2 atan(omega_c dt / 2): an ELECTRON (q/m < 0) turns the positive way.
    want = -2.0 * np.arctan(wc * dt / 2.0)
    assert want > 0.0
    assert abs(got - want) < 1e-9 * abs(want)


def test_move_advances_positions():
    p = ParticleSpecies(Species.electron(), capacity=2)
    p.spawn([1.0], [2.0], [3.0], [10.0], [-5.0], [0.5], [1.0])
    move(p, 0.5)
    assert (p.x[0], p.y[0], p.z[0]) == (6.0, -0.5, 3.25)


def test_spawn_kill_compact():
    p = ParticleSpecies(Species.electron(), capacity=2)
    n = 10
    p.spawn(np.arange(n) * 1.0, np.zeros(n), np.zeros(n), np.zeros(n),
            np.zeros(n), np.zeros(n), np.ones(n))
    assert p.n == n and p.capacity >= n
    p.kill(np.arange(n) % 2 == 0)
    p.compact()
    assert p.n == n // 2
    assert np.allclose(np.sort(p.x[:p.n]), np.arange(1, n, 2))


def test_diagnostics_are_weight_aware():
    e = Species.electron()
    p = ParticleSpecies(e, capacity=4)
    p.spawn([0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1e6, 0.0], [0.0, 0.0],
            [0.0, 0.0], [3.0, 1.0])
    assert abs(p.real_count() - 4.0) < 1e-15
    want_ke = 0.5 * M_E * 3.0 * 1e12
    assert abs(p.kinetic_energy_J() - want_ke) / want_ke < 1e-15
    want_px = M_E * 3.0 * 1e6
    assert abs(p.momentum()[0] - want_px) / want_px < 1e-15
    want_T = (2.0 / 3.0) * (3.0 * 0.5 * M_E * 1e12 / EV) / 4.0
    assert abs(p.temperature_eV() - want_T) / want_T < 1e-12


def test_subcycle_ratio_follows_the_mass_ratio_and_cfl():
    ar = Species.ion("Ar+", ION_MASSES_AMU["Ar+"])
    k_mass = ion_subcycle_ratio(ar, dt=1e-12, h=1e-3)
    assert k_mass == int(round(np.sqrt(ar.mass / M_E)))
    # a CFL cap must be able to override the mass-ratio suggestion
    k_cfl = ion_subcycle_ratio(ar, dt=1e-12, h=1e-6, v_max=1e5)
    assert 1 <= k_cfl < k_mass


def test_species_set_aggregates():
    e = Species.electron()
    ar = Species.ion("Ar+", ION_MASSES_AMU["Ar+"])
    pe = ParticleSpecies(e, capacity=2)
    pe.spawn([0.0], [0.0], [0.0], [1e6], [0.0], [0.0], [5.0])
    pi = ParticleSpecies(ar, capacity=2)
    pi.spawn([0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [5.0])
    s = SpeciesSet([pe, pi])
    assert len(s) == 2
    assert s.by_name("Ar+") is pi
    assert abs(s.net_charge()) < 1e-30          # 5 electrons + 5 Ar+ = neutral
    assert s.kinetic_energy_J() == pe.kinetic_energy_J()
    with pytest.raises(KeyError):
        s.by_name("Kr+")
