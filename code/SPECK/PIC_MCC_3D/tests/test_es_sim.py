"""Integration gates: the closed deposit -> solve -> gather -> push loop.

The headline gate is the cold plasma oscillation.  A Poisson solver can match
an LU oracle exactly and a pusher can be exact in isolation while the two are
wired together wrongly (a sign, a factor of h, a stale field); only a closed
loop with a known analytic frequency catches that.

The target is the DISCRETE dispersion relation of this scheme
(``theory.cold_mode_frequency_ratio``), not a loose band around omega_pe --
the scheme's frequency deficit is computable, so scoring against the
continuum value with a sloppy tolerance would hide real error.  A separate
gate checks that the measured frequency converges to omega_pe as h -> 0.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
from PIC_MCC_3D.MODULES.poisson3d import PERIODIC, PERIODIC_XY_DIRICHLET_Z
from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species, SpeciesSet)
from PIC_MCC_3D.MODULES.theory import (cold_mode_frequency_ratio,
                                       zero_crossing_frequency)
from PIC_MCC_3D.MODULES.units import E_CHARGE, plasma_frequency

N_E = 1.0e14        # m^-3
L = 1.0e-2          # m


def _oscillation_run(nx, per_axis, n_periods=8, wpe_dt=0.05, amp_frac=0.004):
    """A cold, quiet-start electron plasma in a neutralising jelly."""
    h = L / nx
    wpe = plasma_frequency(N_E)
    dt = wpe_dt / wpe
    pop = ParticleSpecies.cold_wave(Species.electron(), per_axis ** 3,
                                    [0, 0, 0], [L, L, L], N_E,
                                    amplitude=amp_frac * L, mode=1, axis=2)
    sim = ESSim(nx, nx, nx, h, SpeciesSet([pop]), dt, bc=PERIODIC,
                background_charge_density=uniform_background_for(pop, L ** 3))
    sim.set_mode(axis=2, mode=1)
    hist = sim.run(int(round(n_periods * 2.0 * np.pi / wpe / dt)))
    f = zero_crossing_frequency(hist["t"], hist["mode_re"])
    return f / (wpe / (2.0 * np.pi)), hist, sim


def test_cold_plasma_oscillates_at_the_discrete_frequency():
    """THE headline gate: measured f/f_pe matches this scheme's own theory."""
    nx = 24
    got, _, _ = _oscillation_run(nx, per_axis=32)
    predicted = cold_mode_frequency_ratio(2.0 * np.pi / L * (L / nx))
    assert abs(got / predicted - 1.0) < 3e-3, (
        f"measured {got:.5f} vs discrete prediction {predicted:.5f}")


def test_oscillation_converges_to_omega_pe_with_resolution():
    """Refining h must drive the measured frequency toward the continuum value."""
    coarse, _, _ = _oscillation_run(16, per_axis=24, n_periods=6)
    fine, _, _ = _oscillation_run(32, per_axis=40, n_periods=6)
    assert abs(fine - 1.0) < abs(coarse - 1.0)
    assert abs(fine - 1.0) < 0.01


def test_oscillation_frequency_is_insensitive_to_dt():
    """The frequency error is SPATIAL; halving dt must not move it much.

    Distinguishes a discretisation error (expected, computable) from a
    time-integration error (would indicate a broken mover).
    """
    a, _, _ = _oscillation_run(16, per_axis=24, n_periods=6, wpe_dt=0.05)
    b, _, _ = _oscillation_run(16, per_axis=24, n_periods=6, wpe_dt=0.025)
    assert abs(a - b) < 2e-3


def test_energy_is_exchanged_between_field_and_particles():
    """A cold start converts field energy to kinetic energy and back."""
    _, hist, _ = _oscillation_run(16, per_axis=24, n_periods=4)
    ke, ue = hist["ke_J"], hist["ue_J"]
    # NOT ke[0] == 0: run() calls prepare(), which retards v by half a step,
    # so the first sample already carries the half-step velocity.  The load
    # itself is exactly cold -- that is gated in test_species3d.
    assert ke[0] < 0.05 * ke.max()
    assert ke.max() > 0.0
    # the total stays bounded: no runaway grid heating over a few periods
    total = ke + ue
    assert total.max() / total[0] < 1.5


def test_periodic_run_conserves_charge_and_particle_count():
    _, hist, sim = _oscillation_run(16, per_axis=24, n_periods=3)
    assert hist["n_macro"][0] == hist["n_macro"][-1]
    assert sim.wall_charge == 0.0           # nothing is lost in a periodic box


# --------------------------------------------------------------- geometry
def test_periodic_boundaries_wrap_particles():
    nx, h = 8, 1e-3
    pop = ParticleSpecies(Species.electron(), capacity=4)
    box = nx * h
    pop.spawn([0.5 * box], [0.5 * box], [0.99 * box], [0.0], [0.0], [1e7],
              [1.0])
    # a periodic box needs the neutralising jelly or the solve refuses
    sim = ESSim(nx, nx, nx, h, SpeciesSet([pop]), 1e-9, bc=PERIODIC,
                background_charge_density=uniform_background_for(pop, box ** 3))
    sim.step()
    assert 0.0 <= pop.z[0] < box            # wrapped, not lost


def test_electrodes_absorb_particles_and_tally_the_charge():
    """Particles leaving through an electrode are removed and BOOKKEPT."""
    nx, nz, h = 8, 10, 1e-3
    pop = ParticleSpecies(Species.electron(), capacity=4)
    Lz = (nz - 1) * h
    pop.spawn([4e-3], [4e-3], [0.9 * Lz], [0.0], [0.0], [1e8], [7.0])
    sim = ESSim(nx, nx, nz, h, SpeciesSet([pop]), 1e-10,
                bc=PERIODIC_XY_DIRICHLET_Z)
    sim.step()
    assert pop.n == 0
    assert abs(sim.wall_charge - (-E_CHARGE * 7.0)) < 1e-30
    assert sim.wall_energy_J > 0.0


def test_applied_electrode_voltage_drives_a_field():
    """A vacuum gap at fixed potentials produces the expected uniform E_z."""
    nx, nz, h, V = 6, 21, 1e-3, 100.0
    sim = ESSim(nx, nx, nz, h, SpeciesSet([]), 1e-12,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=V)
    sim.solve_fields()
    want = -V / ((nz - 1) * h)
    assert abs(float(sim.E[2][3, 3, nz // 2]) - want) / abs(want) < 1e-9


# -------------------------------------------------------------- subcycling
def test_subcycled_heavy_species_matches_an_unsubcycled_run():
    """Subcycling must be a cost optimisation, not a physics change.

    A uniform applied field is integrated exactly by the leapfrog, so the two
    trajectories should agree to round-off; any drift here would mean the
    subcycle bookkeeping (which steps happen, with which dt) is wrong.
    """
    nx, nz, h, V = 6, 21, 1e-3, 500.0
    dt = 1e-11
    out = {}
    for k in (1, 8):
        pop = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                              capacity=4)
        pop.spawn([3e-3], [3e-3], [0.5 * (nz - 1) * h], [0.0], [0.0], [0.0],
                  [1.0])
        pop.subcycle = k
        sim = ESSim(nx, nx, nz, h, SpeciesSet([pop]), dt,
                    bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=V)
        sim.prepare()
        for _ in range(64):
            sim.step()
        out[k] = (pop.z[0], pop.vz[0])
    assert abs(out[1][1] - out[8][1]) / abs(out[1][1]) < 1e-9
    assert abs(out[1][0] - out[8][0]) / abs(out[1][0]) < 1e-3


# --------------------------------------------------------------- stability
def test_stability_reports_both_explicit_walls():
    """The dx<lambda_D and omega_pe dt<2 walls must be reported, not hidden."""
    nx, h = 16, 1e-4
    sim = ESSim(nx, nx, nx, h, SpeciesSet([]), 1e-11, bc=PERIODIC)
    s = sim.stability(n_e_m3=1e16, Te_eV=1.0)
    assert s["dx_over_debye"] > 0 and np.isfinite(s["dx_over_debye"])
    assert s["wpe_dt"] > 0
    assert s["resolved"] == (s["dx_over_debye"] < 1.0 and s["wpe_dt"] < 2.0)
    # a deliberately under-resolved configuration must report itself as such
    bad = ESSim(nx, nx, nx, 1e-2, SpeciesSet([]), 1e-8, bc=PERIODIC)
    assert bad.stability(1e16, 1.0)["resolved"] is False


@pytest.mark.parametrize("bc", [PERIODIC, PERIODIC_XY_DIRICHLET_Z])
def test_empty_simulation_runs(bc):
    """No particles is a valid (vacuum) run, not a crash."""
    sim = ESSim(6, 6, 8, 1e-3, SpeciesSet([]), 1e-12, bc=bc)
    hist = sim.run(3)
    assert np.all(np.asarray(hist["ke_J"]) == 0.0)


# ----------------------------------------------------- collisions in the loop
def _mcc_bits():
    """(gas, cross-sections) or None if this checkout has no LXCat tree."""
    from PIC_MCC_3D.MODULES.cross_sections import lxcat_roots, mixture_set
    from PIC_MCC_3D.MODULES.gas import pure
    if not any(r.is_dir() for r in lxcat_roots()):
        return None
    g = pure("Ar", 0.01, 300.0)
    return g, mixture_set(g)


def test_collisions_run_inside_the_loop_and_make_ions():
    """The PIC and the MCC actually composed: a driven gap ionizes.

    Not a physics claim about a discharge -- it is the wiring gate. Electrons
    accelerated by the applied field must produce ionization events, those
    events must appear in the run's tally, and every one must leave a real
    ion macro-particle behind.
    """
    bits = _mcc_bits()
    if bits is None:
        pytest.skip("no LXCat data tree reachable")
    from PIC_MCC_3D.MODULES.mcc3d import MCC3D
    g, cs = bits

    nx, nz, h = 8, 21, 2e-4
    rng = np.random.default_rng(3)
    Lz = (nz - 1) * h
    electrons = ParticleSpecies.maxwellian_box(
        Species.electron(), 4000, [0, 0, 0.25 * Lz],
        [nx * h, nx * h, 0.75 * Lz], 1e15, 5.0, rng)
    ions = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                           capacity=4096)
    ions.subcycle = 16
    mcc = MCC3D(cs, mixture=g, ion_species={"Ar": ions})
    sim = ESSim(nx, nx, nz, h, SpeciesSet([electrons, ions]), 2e-12,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=800.0,
                mcc=mcc, rng=rng)
    sim.run(60)

    kinds = cs.kinds
    n_ion = int(sim.collision_tally[[i for i, k in enumerate(kinds)
                                     if k == "ionization"]].sum())
    assert int(sim.collision_tally.sum()) > 0, "no collisions occurred"
    assert n_ion > 0, "no ionization -- raise the field or the density"
    # every ionization left a real ion; some may have hit a wall since
    assert ions.n + _wall_ions(sim) == n_ion
    assert np.all(sim.rho_heavy == 0.0), "ions were registered; none should be frozen"


def _wall_ions(sim):
    """Ions collected by the electrodes, from the wall-charge tally."""
    from PIC_MCC_3D.MODULES.units import E_CHARGE
    # electrons carry -e and ions +e; the tally is signed, so recover the ion
    # count only when it is unambiguous (this helper is used where the ion
    # weight is 1 per macro-particle)
    return 0 if sim.wall_charge <= 0 else int(round(sim.wall_charge / E_CHARGE))


def test_unregistered_ions_become_immobile_grid_charge():
    """Without an ion species the charge is still accounted for, on the grid."""
    bits = _mcc_bits()
    if bits is None:
        pytest.skip("no LXCat data tree reachable")
    from PIC_MCC_3D.MODULES.mcc3d import MCC3D
    g, cs = bits
    nx, nz, h = 8, 21, 2e-4
    rng = np.random.default_rng(3)
    Lz = (nz - 1) * h
    electrons = ParticleSpecies.maxwellian_box(
        Species.electron(), 4000, [0, 0, 0.25 * Lz],
        [nx * h, nx * h, 0.75 * Lz], 1e15, 5.0, rng)
    sim = ESSim(nx, nx, nz, h, SpeciesSet([electrons]), 2e-12,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=800.0,
                mcc=MCC3D(cs, mixture=g), rng=rng)
    sim.run(60)
    n_ion = int(sim.collision_tally[[i for i, k in enumerate(cs.kinds)
                                     if k == "ionization"]].sum())
    assert n_ion > 0
    deposited = float(np.sum(sim.rho_heavy)) * h ** 3
    from PIC_MCC_3D.MODULES.units import E_CHARGE
    assert abs(deposited - n_ion * E_CHARGE * electrons.w[0]) \
        < 1e-9 * abs(deposited)


def test_mcc_target_inference_refuses_to_guess():
    bits = _mcc_bits()
    if bits is None:
        pytest.skip("no LXCat data tree reachable")
    from PIC_MCC_3D.MODULES.mcc3d import MCC3D
    g, cs = bits
    two = SpeciesSet([ParticleSpecies(Species.electron()),
                      ParticleSpecies(Species.electron())])
    with pytest.raises(ValueError, match="cannot infer the MCC target"):
        ESSim(6, 6, 8, 1e-3, two, 1e-12, mcc=MCC3D(cs, mixture=g))
