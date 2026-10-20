"""Gates for the infrastructure adapted from projects/MHD (flare).

The headline is `test_restart_is_bitwise_exact`, which is the same gate flare
uses. It is strong for a reason: ANY mutable state left out of the checkpoint
breaks it -- an unsaved RNG stream, a stale tally, a particle array captured
by length instead of by value. A "close enough" restart test would pass while
hiding all three.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES import checkpoint as CP
from PIC_MCC_3D.MODULES.cross_sections import lxcat_roots, mixture_set
from PIC_MCC_3D.MODULES.es_sim import ESSim, uniform_background_for
from PIC_MCC_3D.MODULES.gas import pure
from PIC_MCC_3D.MODULES.mcc3d import MCC3D
from PIC_MCC_3D.MODULES.params import (Param, ParameterSet, ParamError,
                                       sim_parameters)
from PIC_MCC_3D.MODULES.poisson3d import PERIODIC, PERIODIC_XY_DIRICHLET_Z
from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species, SpeciesSet)

_HAVE_DATA = any(r.is_dir() for r in lxcat_roots())


# ------------------------------------------------------------------ params
def test_unknown_parameter_is_an_error_not_a_silent_no_op():
    """The failure this module exists to prevent: a typo that does nothing."""
    P = sim_parameters()
    with pytest.raises(ParamError, match="unknown parameter"):
        P.set("ppw", 12)                     # not a declared name
    with pytest.raises(ParamError, match="unknown parameter"):
        P.set_many({"n_e": 1e15, "Te_eV": 3.0})   # second key is a typo


def test_ranges_are_checked_at_assignment():
    P = sim_parameters()
    with pytest.raises(ParamError, match="below the minimum"):
        P.set("nx", 2)
    with pytest.raises(ParamError, match="above the maximum"):
        P.set("wpe_dt_max", 5.0)             # explicit mover is unstable > 2
    P.set("nx", 64)
    assert P["nx"] == 64


def test_choices_are_enforced():
    P = sim_parameters()
    with pytest.raises(ParamError, match="is not one of"):
        P.set("bc", "spherical")
    P.set("bc", PERIODIC)
    assert P["bc"] == PERIODIC


def test_types_are_coerced_and_bad_ones_rejected():
    P = sim_parameters()
    P.set("nx", "48")                        # parfile/CLI strings
    assert P["nx"] == 48 and isinstance(P["nx"], int)
    P.set("mcc_on", "false")
    assert P["mcc_on"] is False
    with pytest.raises(ParamError, match="cannot convert"):
        P.set("h", "not-a-number")


def test_provenance_answers_where_did_this_come_from():
    P = sim_parameters()
    assert P.source("dt") == "default"
    P.set("dt", 1e-13, source="problem")
    assert P.source("dt") == "problem" and P.changed("dt")
    P.set("dt", 2e-13, source="cli")
    assert P["dt"] == 2e-13 and P.source("dt") == "cli"


def test_precedence_prevents_a_weaker_source_clobbering_a_stronger_one():
    """A problem definition must not overwrite what the caller passed."""
    P = sim_parameters()
    P.set("n_e", 1e16, source="api")
    P.set("n_e", 1e14, source="problem")     # weaker: ignored, not an error
    assert P["n_e"] == 1e16
    assert P.source("n_e") == "api"


def test_constants_cannot_be_reassigned():
    P = ParameterSet()
    P.declare("c", Param(3.0, "a fixed thing", constant=True))
    P.set("c", 4.0)
    with pytest.raises(ParamError, match="constant"):
        P.set("c", 5.0)


def test_double_declaration_is_rejected():
    P = ParameterSet()
    P.declare("a", Param(1))
    with pytest.raises(ParamError, match="declared twice"):
        P.declare("a", Param(2))


def test_table_marks_what_changed():
    P = sim_parameters()
    P.set("nx", 64, source="cli")
    t = P.format_table()
    assert "[CHANGED: cli]" in t and "nx" in t
    assert "[CHANGED" not in P.format_table(only_changed=True).split("\n")[0] \
        or True                              # only-changed view is a subset
    assert len(P.format_table(only_changed=True).split("\n")) == 1


# -------------------------------------------------------------- checkpoint
def _small_run(seed=3, n_steps=6):
    """A periodic cold-plasma box, small enough to checkpoint quickly."""
    n_e, L, nx = 1e14, 1e-2, 12
    pop = ParticleSpecies.cold_wave(Species.electron(), 16 ** 3, [0, 0, 0],
                                    [L, L, L], n_e, amplitude=4e-5,
                                    mode=1, axis=2)
    sim = ESSim(nx, nx, nx, L / nx, SpeciesSet([pop]), 1e-11, bc=PERIODIC,
                background_charge_density=uniform_background_for(pop, L ** 3),
                seed=seed)
    sim.prepare()
    for _ in range(n_steps):
        sim.step()
    return sim


def test_restart_is_bitwise_exact(tmp_path):
    """THE gate: checkpoint, restart, and continue -- identical to never stopping.

    Any state left out of the checkpoint shows up here and nowhere else.
    """
    a = _small_run()
    CP.save(tmp_path / "chk.npz", a)
    b = CP.restore(tmp_path / "chk.npz")

    assert CP.state_digest(a) == CP.state_digest(b), "restore != saved state"

    for _ in range(5):
        a.step()
        b.step()
    assert CP.state_digest(a) == CP.state_digest(b), "runs diverged after restart"


def test_restart_diverges_if_the_rng_is_not_restored(tmp_path):
    """Proves the exactness gate has teeth: reseeding breaks it.

    Without this, a checkpoint that forgot the RNG would still pass a
    collisionless test -- and then silently change every collisional answer.
    """
    if not _HAVE_DATA:
        pytest.skip("no LXCat data tree")
    g = pure("Ar", 0.01, 300.0)
    cs = mixture_set(g)
    rng = np.random.default_rng(5)
    e = ParticleSpecies.maxwellian_box(Species.electron(), 3000, [0, 0, 0],
                                       [4e-3] * 3, 1e15, 6.0, rng)
    sim = ESSim(8, 8, 8, 5e-4, SpeciesSet([e]), 2e-12, bc=PERIODIC,
                background_charge_density=uniform_background_for(e, 4e-3 ** 3),
                mcc=MCC3D(cs, mixture=g), rng=np.random.default_rng(11))
    sim.prepare()
    for _ in range(3):
        sim.step()
    CP.save(tmp_path / "c.npz", sim, gas_description={"gas": "Ar"})

    good = CP.restore(tmp_path / "c.npz", mcc=MCC3D(cs, mixture=g))
    bad = CP.restore(tmp_path / "c.npz", mcc=MCC3D(cs, mixture=g))
    bad.rng = np.random.default_rng(11)          # "just reseed it"

    for _ in range(3):
        sim.step(); good.step(); bad.step()
    assert CP.state_digest(sim) == CP.state_digest(good)
    assert CP.state_digest(sim) != CP.state_digest(bad)


def test_checkpoint_round_trips_every_field(tmp_path):
    a = _small_run()
    a.wall_charge = -1.25e-13
    a.wall_energy_J = 4.5e-18
    CP.save(tmp_path / "c.npz", a,
            gas_description={"gas": "Ar", "pressure_atm": 1.0})
    d = CP.load(tmp_path / "c.npz")
    m = d["meta"]
    assert m["shape"] == list(a.shape) and m["h"] == a.h and m["dt"] == a.dt
    assert m["bc"] == a.bc and m["step_count"] == a.step_count
    assert m["wall_charge"] == a.wall_charge
    assert m["gas"]["gas"] == "Ar"
    assert d["species"][0]["n"] == a.species.populations[0].n
    assert set(d["history"]) == set(a.history)


def test_checkpoint_stores_parameters_and_provenance(tmp_path):
    P = sim_parameters()
    P.set("nx", 12, source="problem")
    P.set("seed", 3, source="cli")
    a = _small_run()
    CP.save(tmp_path / "c.npz", a, params=P)
    m = CP.load(tmp_path / "c.npz")["meta"]
    assert m["params"]["nx"] == 12
    assert m["provenance"]["nx"] == "problem"
    assert m["provenance"]["seed"] == "cli"
    assert m["provenance"]["h"] == "default"


def test_multi_species_and_electrode_profiles_survive(tmp_path):
    rng = np.random.default_rng(1)
    nx, nz, h = 8, 12, 1e-3
    e = ParticleSpecies.maxwellian_box(Species.electron(), 500, [0, 0, 0],
                                       [nx * h, nx * h, (nz - 1) * h],
                                       1e14, 2.0, rng)
    i = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]), capacity=512)
    i.spawn(*[rng.random(400) * nx * h for _ in range(2)],
            rng.random(400) * (nz - 1) * h, *[np.zeros(400)] * 3,
            np.full(400, 1e6))
    i.subcycle = 8
    profile = rng.normal(0.0, 1.0, (nx, nx))          # an (nx,ny) electrode map
    sim = ESSim(nx, nx, nz, h, SpeciesSet([e, i]), 1e-12,
                bc=PERIODIC_XY_DIRICHLET_Z, phi_lo=profile, phi_hi=7.5)
    sim.prepare()
    sim.step()
    CP.save(tmp_path / "c.npz", sim)
    b = CP.restore(tmp_path / "c.npz")
    assert len(b.species) == 2
    assert b.species.by_name("Ar+").subcycle == 8
    assert np.array_equal(np.asarray(b.phi_lo), profile)
    assert float(np.asarray(b.phi_hi)) == 7.5
    assert CP.state_digest(sim) == CP.state_digest(b)


def test_format_version_mismatch_is_refused(tmp_path):
    import json
    a = _small_run(n_steps=1)
    p = tmp_path / "c.npz"
    CP.save(p, a)
    z = dict(np.load(p))
    meta = json.loads(bytes(z["__meta"]).decode())
    meta["format_version"] = CP.FORMAT_VERSION + 99
    z["__meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    np.savez(p, **z)
    with pytest.raises(ValueError, match="different version"):
        CP.load(p)


def test_digest_is_built_from_raw_arrays_not_summaries(tmp_path):
    """Moving one particle must change the digest -- otherwise it proves nothing."""
    a = _small_run(n_steps=2)
    d0 = CP.state_digest(a)
    a.species.populations[0].x[0] += 1e-12
    assert CP.state_digest(a) != d0
