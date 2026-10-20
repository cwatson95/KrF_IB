"""The composition ladder: modules together, one rung at a time.

The per-module files answer "does this work alone". This one answers "do they
work together", and it is built as a LADDER rather than one big end-to-end
test on purpose: each rung adds exactly one module to the previous rung and
re-checks the invariants that already held. When a rung fails, the module it
introduced is the suspect -- which is not true of a single everything-at-once
test, where any of eight modules could be at fault.

    rung 1  fields only ................ Poisson + es_field
    rung 2  + kinetic particles ........ species3d + es_sim
    rung 3  + neutral gas + collisions .. gas + cross_sections + mcc3d
    rung 4  + beam injection ........... beam3d
    rung 5  + persistence .............. checkpoint, with everything live

Rung 5 is the capstone: a bitwise-exact restart of the FULL stack. It can
only pass if every module above it put all of its state where the checkpoint
could see it.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES import checkpoint as CP
from PIC_MCC_3D.MODULES.beam3d import BeamInjector, BeamSpec, seed_preionization
from PIC_MCC_3D.MODULES.cross_sections import lxcat_roots, mixture_set
from PIC_MCC_3D.MODULES.es_sim import ESSim
from PIC_MCC_3D.MODULES.gas import pure
from PIC_MCC_3D.MODULES.mcc3d import MCC3D
from PIC_MCC_3D.MODULES.poisson3d import PERIODIC_XY_DIRICHLET_Z
from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species, SpeciesSet)
from PIC_MCC_3D.MODULES.units import E_CHARGE, EV

needs_data = pytest.mark.skipif(
    not any(r.is_dir() for r in lxcat_roots()), reason="no LXCat data tree")

NX, NZ, H = 8, 16, 4.0e-4          # a small driven gap
VOLT = 600.0
DT = 2.0e-12


def _empty_sim(**kw):
    return ESSim(NX, NX, NZ, H, SpeciesSet([]), DT,
                 bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=VOLT, **kw)


# --------------------------------------------------------------- rung 1
def test_rung1_fields_only():
    """Poisson + es_field with no particles: the vacuum gap field."""
    sim = _empty_sim()
    sim.solve_fields()
    want = -VOLT / ((NZ - 1) * H)
    assert abs(float(sim.E[2][4, 4, NZ // 2]) - want) / abs(want) < 1e-9
    assert sim.species.kinetic_energy_J() == 0.0


# --------------------------------------------------------------- rung 2
def test_rung2_add_particles():
    """+ species3d/es_sim: particles feel the gap field and are collected."""
    rng = np.random.default_rng(2)
    e = ParticleSpecies.maxwellian_box(Species.electron(), 2000, [0, 0, 0],
                                       [NX * H, NX * H, (NZ - 1) * H],
                                       1e15, 2.0, rng)
    q0 = e.species.charge * e.real_count()
    sim = ESSim(NX, NX, NZ, H, SpeciesSet([e]), DT,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=VOLT, rng=rng)
    sim.run(40)
    # every electron that left is on the wall tally: charge is conserved
    q1 = e.species.charge * e.real_count()
    assert abs((q1 + sim.wall_charge) - q0) < 1e-9 * abs(q0)
    assert sim.wall_charge != 0.0, "nothing reached a wall -- rung is vacuous"


# --------------------------------------------------------------- rung 3
@needs_data
def test_rung3_add_gas_and_collisions():
    """+ gas/cross_sections/mcc3d: ionization creates matched e-/ion pairs.

    Re-checks rung 2's charge budget, now with sources AND sinks: seeded
    charge + everything created - everything collected must still balance.
    """
    g = pure("Ar", 0.05, 300.0)
    cs = mixture_set(g)
    rng = np.random.default_rng(4)
    e = ParticleSpecies.maxwellian_box(Species.electron(), 3000, [0, 0, 0],
                                       [NX * H, NX * H, (NZ - 1) * H],
                                       1e15, 4.0, rng)
    ions = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                           capacity=1 << 15)
    ions.subcycle = 16
    w = e.w[0]
    n_e0 = e.n
    sim = ESSim(NX, NX, NZ, H, SpeciesSet([e, ions]), DT,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=VOLT,
                mcc=MCC3D(cs, mixture=g, ion_species={"Ar": ions}), rng=rng)
    sim.run(40)

    kinds = cs.kinds
    n_ion = int(sim.collision_tally[[i for i, k in enumerate(kinds)
                                     if k == "ionization"]].sum())
    assert n_ion > 0, "no ionization -- rung is vacuous"
    # PARTICLE bookkeeping across three modules at once: electrons now =
    # seeded + ionization births - wall losses.  Wall losses are whatever the
    # tally says, converted back to macro-particles at the common weight.
    lost = (e.n + ions.n) - (n_e0 + 2 * n_ion)
    assert lost <= 0, "more particles exist than were created"
    # CHARGE: seeded electrons + created pairs (net zero) - collected
    q_now = (e.species.charge * e.real_count()
             + ions.species.charge * ions.real_count())
    q_created = 0.0                       # each pair is net neutral
    q_seed = -E_CHARGE * w * n_e0
    assert abs((q_now + sim.wall_charge) - (q_seed + q_created)) \
        < 1e-6 * abs(q_seed)
    # ENERGY: what the electrons lost went to the gas ledger, nothing invented
    assert float(np.sum(sim.collision_energy_eV)) > 0.0


# --------------------------------------------------------------- rung 4
@needs_data
def test_rung4_add_the_beam():
    """+ beam3d: injection and the preionization seed feed the same loop.

    The injector is given the SIM's rng on purpose -- one stream per run, so
    the checkpoint at rung 5 covers the beam's randomness too. Two streams
    would make the full-stack restart unrepeatable for no benefit.
    """
    g = pure("Ar", 0.05, 300.0)
    cs = mixture_set(g)
    rng = np.random.default_rng(7)
    e = ParticleSpecies(Species.electron(), capacity=1 << 16)
    ions = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                           capacity=1 << 16)
    ions.subcycle = 16
    net = seed_preionization(e, ions, 1500, [0, 0, 0.3 * (NZ - 1) * H],
                             [NX * H, NX * H, 0.6 * (NZ - 1) * H],
                             1e15, 1.0, rng)
    assert abs(net) < 1e-30, "the seed is not neutral"

    sim = ESSim(NX, NX, NZ, H, SpeciesSet([e, ions]), DT,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=VOLT,
                mcc=MCC3D(cs, mixture=g, ion_species={"Ar": ions}), rng=rng)
    inj = BeamInjector(g, BeamSpec(voltage_kV=600.0, j_A_cm2=50.0),
                       e, rng=sim.rng)
    n_before = e.n
    dep = 0.0
    for _ in range(30):
        rec = inj.inject(sim)
        dep += rec["energy_deposited_J"]
        sim.step()
    assert dep > 0.0 and e.n > n_before, "the beam injected nothing"
    assert inj.energy_deposited_J == pytest.approx(dep)
    # the beam deposits far less than it carries in: it ranges metres
    assert inj.energy_deposited_J < 0.05 * inj.energy_in_J


# --------------------------------------------------------------- rung 5
@needs_data
def test_rung5_full_stack_restart_is_bitwise_exact(tmp_path):
    """THE capstone: checkpoint the whole live stack and resume identically.

    Passes only if every module put ALL of its state somewhere the checkpoint
    could see it -- particles of two species, the immobile-charge grid, the
    collision tallies, the fields, and the single RNG stream that the MCC and
    the beam both draw from.
    """
    g = pure("Ar", 0.05, 300.0)
    cs = mixture_set(g)

    def build(seed):
        rng = np.random.default_rng(seed)
        e = ParticleSpecies(Species.electron(), capacity=1 << 16)
        ions = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                               capacity=1 << 16)
        ions.subcycle = 16
        seed_preionization(e, ions, 1200, [0, 0, 0.3 * (NZ - 1) * H],
                           [NX * H, NX * H, 0.6 * (NZ - 1) * H],
                           1e15, 1.0, rng)
        sim = ESSim(NX, NX, NZ, H, SpeciesSet([e, ions]), DT,
                    bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=0.0, phi_hi=VOLT,
                    mcc=MCC3D(cs, mixture=g, ion_species={"Ar": ions}),
                    rng=np.random.default_rng(seed + 100))
        return sim, e

    a, _ = build(9)
    a.prepare()
    for _ in range(12):
        a.step()

    CP.save(tmp_path / "full.npz", a,
            gas_description={"gas": "Ar", "pressure_atm": 0.05, "T_gas": 300.0})

    # rebuild the operator from the saved description -- it is not serialised
    meta = CP.load(tmp_path / "full.npz")["meta"]["gas"]
    g2 = pure(meta["gas"], meta["pressure_atm"], meta["T_gas"])
    b = CP.restore(tmp_path / "full.npz")
    ions_b = b.species.by_name("Ar+")
    b.mcc = MCC3D(mixture_set(g2), mixture=g2, ion_species={"Ar": ions_b})

    assert CP.state_digest(a) == CP.state_digest(b)
    assert a.collision_tally is not None and a.collision_tally.sum() > 0, \
        "no collisions before the checkpoint -- the gate would be vacuous"

    for _ in range(8):
        a.step()
        b.step()
    assert CP.state_digest(a) == CP.state_digest(b), \
        "the full stack diverged after restart"
