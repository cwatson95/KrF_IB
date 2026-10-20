"""Gates for the electron-neutral MCC (missing module 3, electron half).

Three classes of gate:

  * BOOKKEEPING -- the electron energy ledger closes exactly, and every
    ionization produces exactly one electron and one ion. These are
    invariants, so they are held at round-off, not at a tolerance.
  * PHYSICS -- the measured collision rate reproduces n sigma v for a
    monoenergetic beam, and thresholds are respected.
  * PROVENANCE -- the adapted operator reproduces ``ARPAE/pic3d/mcc3d.py``
    EXACTLY on a shared RNG stream. That module is itself oracle-matched
    against ``pic1d.mcc_ref``, so this gate is what transfers the lineage's
    validation to this copy rather than re-asserting it.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.cross_sections import lxcat_roots, mixture_set
from PIC_MCC_3D.MODULES.gas import pure
from PIC_MCC_3D.MODULES.mcc3d import MCC3D
from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species)
from PIC_MCC_3D.MODULES.units import EV, M_E

pytestmark = pytest.mark.skipif(
    not any(r.is_dir() for r in lxcat_roots()),
    reason="no LXCat data tree reachable from this checkout")

BOX = [0.0, 0.0, 0.0], [1e-2, 1e-2, 1e-2]


def _electrons(n, Te, rng, seed_box=BOX):
    return ParticleSpecies.maxwellian_box(Species.electron(), n, seed_box[0],
                                          seed_box[1], 1e14, Te, rng)


def _ion_pop(name="Ar+"):
    return ParticleSpecies(Species.ion(name, ION_MASSES_AMU[name]),
                           capacity=1024)


# ------------------------------------------------------------- bookkeeping
def test_electron_energy_ledger_closes_exactly():
    """KE(before) - KE(after, incl. secondaries) = sum(ledger)."""
    rng = np.random.default_rng(4)
    g = pure("Ar", 0.01, 300.0)
    mcc = MCC3D(mixture_set(g), mixture=g, ion_species={"Ar": _ion_pop()})
    p = _electrons(20000, 8.0, rng)
    ke0 = p.kinetic_energy_J()
    ev = mcc.collide_adaptive(p, 2e-10, rng)
    ke1 = p.kinetic_energy_J()
    ledger = float(np.sum(ev["energy_eV"])) * EV
    assert int(ev["tally"].sum()) > 1000, "no collisions -- gate is vacuous"
    assert abs((ke0 - ke1) - ledger) / ke0 < 1e-10


def test_every_ionization_makes_one_electron_and_one_ion():
    rng = np.random.default_rng(4)
    g = pure("Ar", 0.01, 300.0)
    cs = mixture_set(g)
    ions = _ion_pop()
    mcc = MCC3D(cs, mixture=g, ion_species={"Ar": ions})
    p = _electrons(20000, 8.0, rng)
    n0 = p.n
    ev = mcc.collide_adaptive(p, 2e-10, rng)
    n_ion = int(ev["tally"][[i for i, k in enumerate(cs.kinds)
                             if k == "ionization"]].sum())
    assert n_ion > 100
    assert p.n - n0 == n_ion            # one new electron per event
    assert ions.n == n_ion              # one ion per event
    assert ev["ion_w"].size == 0        # all births went into the species


def test_unregistered_targets_fall_back_to_immobile_charge():
    """Without an ion species the events are still reported, never dropped."""
    rng = np.random.default_rng(4)
    g = pure("Ar", 0.01, 300.0)
    cs = mixture_set(g)
    mcc = MCC3D(cs, mixture=g)                    # no ion_species
    p = _electrons(20000, 8.0, rng)
    ev = mcc.collide_adaptive(p, 2e-10, rng)
    n_ion = int(ev["tally"][[i for i, k in enumerate(cs.kinds)
                             if k == "ionization"]].sum())
    assert ev["ion_w"].size == n_ion
    assert ev["ion_pos"][0].size == n_ion
    assert ev["gas_reservoir_eV"] == 0.0


def test_born_ions_carry_the_neutral_thermal_energy_and_it_is_booked():
    rng = np.random.default_rng(11)
    g = pure("Ar", 0.01, 300.0)
    ions = _ion_pop()
    mcc = MCC3D(mixture_set(g), mixture=g, ion_species={"Ar": ions})
    ev = mcc.collide_adaptive(_electrons(20000, 8.0, rng), 2e-10, rng)
    assert ions.n > 0
    assert ev["gas_reservoir_eV"] > 0.0
    # the booked reservoir energy IS the ions' kinetic energy
    assert abs(ev["gas_reservoir_eV"] * EV - ions.kinetic_energy_J()) \
        / ions.kinetic_energy_J() < 1e-12
    # and it is thermal: (2/3)<eps> ~ the gas temperature in eV
    kT_eV = 300.0 * 1.380649e-23 / EV
    assert abs(ions.temperature_eV() / kT_eV - 1.0) < 0.15


def test_ions_can_be_born_at_rest_instead():
    rng = np.random.default_rng(11)
    g = pure("Ar", 0.01, 300.0)
    ions = _ion_pop()
    mcc = MCC3D(mixture_set(g), mixture=g, ion_species={"Ar": ions},
                born_ion_thermal=False)
    ev = mcc.collide_adaptive(_electrons(20000, 8.0, rng), 2e-10, rng)
    assert ions.n > 0
    assert ions.kinetic_energy_J() == 0.0
    assert ev["gas_reservoir_eV"] == 0.0


def test_attachment_removes_electrons():
    rng = np.random.default_rng(5)
    g = pure("F2", 0.01, 300.0)
    cs = mixture_set(g)
    assert "attachment" in cs.kinds
    mcc = MCC3D(cs, mixture=g)
    p = _electrons(20000, 2.0, rng)
    n0 = p.n
    ev = mcc.collide_adaptive(p, 2e-10, rng)
    n_att = int(ev["tally"][[i for i, k in enumerate(cs.kinds)
                             if k == "attachment"]].sum())
    n_ion = int(ev["tally"][[i for i, k in enumerate(cs.kinds)
                             if k == "ionization"]].sum())
    assert n_att > 0
    assert p.n - n0 == n_ion - n_att
    assert ev["att_w"].size == n_att


# ---------------------------------------------------------------- physics
def test_measured_rate_reproduces_n_sigma_v():
    """A monoenergetic beam must collide at exactly n sigma(eps) v.

    The strongest physics check available without a swarm experiment: the
    expected count is computed from the SAME tables, so any disagreement is
    in the sampling, not the data.
    """
    rng = np.random.default_rng(9)
    g = pure("Ar", 0.001, 300.0)
    cs = mixture_set(g)
    mcc = MCC3D(cs, mixture=g)
    eps0 = 1.0                                   # eV, below every threshold
    v0 = np.sqrt(2.0 * eps0 * EV / M_E)
    n = 200000
    p = ParticleSpecies(Species.electron(), capacity=n + 16)
    z = np.zeros(n)
    p.spawn(rng.random(n) * 1e-2, rng.random(n) * 1e-2, rng.random(n) * 1e-2,
            np.full(n, v0), z, z, np.ones(n))
    dt = 1e-12
    nu_expected = float(cs.nu_total(np.array([eps0]))[0])
    assert nu_expected * dt < 0.05, "dt too large for a single-collision count"
    ev = mcc.collide_adaptive(p, dt, rng)
    got = int(ev["tally"].sum())
    want = n * (1.0 - np.exp(-nu_expected * dt))
    # Poisson counting error on ~`want` events
    assert abs(got - want) < 5.0 * np.sqrt(want), (got, want)


def test_no_inelastic_events_below_threshold():
    """Cold electrons in argon can only scatter elastically."""
    rng = np.random.default_rng(3)
    g = pure("Ar", 0.1, 300.0)
    cs = mixture_set(g)
    mcc = MCC3D(cs, mixture=g)
    p = _electrons(50000, 0.05, rng)             # ~0.05 eV, far below 11.5 eV
    ev = mcc.collide_adaptive(p, 1e-11, rng)
    for k, kind in enumerate(cs.kinds):
        if kind in ("excitation", "ionization"):
            assert ev["tally"][k] == 0, f"{kind} fired below threshold"
    assert ev["tally"].sum() > 0


def test_elastic_collisions_cool_electrons_slowly():
    """Elastic loss per collision is ~2m/M -- a very small fraction."""
    rng = np.random.default_rng(6)
    g = pure("Ar", 0.1, 300.0)
    mcc = MCC3D(mixture_set(g), mixture=g)
    p = _electrons(50000, 0.5, rng)
    ke0 = p.kinetic_energy_J()
    for _ in range(5):
        mcc.collide_adaptive(p, 1e-11, rng)
    frac = (ke0 - p.kinetic_energy_J()) / ke0
    assert 0.0 < frac < 0.05, frac


def test_zero_timestep_produces_no_collisions():
    rng = np.random.default_rng(2)
    g = pure("Ar", 1.0, 300.0)
    mcc = MCC3D(mixture_set(g), mixture=g)
    p = _electrons(1000, 5.0, rng)
    ev = mcc.collide_adaptive(p, 0.0, rng)
    assert ev["tally"].sum() == 0


def test_empty_population_is_a_no_op():
    rng = np.random.default_rng(2)
    g = pure("Ar", 1.0, 300.0)
    mcc = MCC3D(mixture_set(g), mixture=g)
    p = ParticleSpecies(Species.electron(), capacity=4)
    ev = mcc.collide_adaptive(p, 1e-12, rng)
    assert ev["tally"].sum() == 0
    assert ev["ion_w"].size == 0


def test_rate_coefficients_are_physical():
    """<sigma v> for Ar at 4 eV: elastic dominant, ionization tiny."""
    g = pure("Ar", 1.0, 300.0)
    mcc = MCC3D(mixture_set(g), mixture=g)
    r = mcc.rate_coefficients(4.0)
    assert r["elastic"] > 1e-14
    assert r["elastic"] > 100.0 * r.get("ionization", 0.0)


def test_thermal_birth_requires_a_mixture():
    g = pure("Ar", 1.0, 300.0)
    with pytest.raises(ValueError, match="needs a GasMixture"):
        MCC3D(mixture_set(g), ion_species={"Ar": _ion_pop()})


# ------------------------------------------------------------- provenance
def _arpae_pic3d():
    """Import ARPAE's pic3d if this checkout has it."""
    here = Path(__file__).resolve()
    for base in here.parents:
        cand = base / "ARPAE"
        if (cand / "pic3d" / "mcc3d.py").is_file():
            if str(cand) not in sys.path:
                sys.path.insert(0, str(cand))
            import pic3d.cross_sections as acs
            import pic3d.mcc3d as amcc
            import pic3d.particles3d as apart
            return acs, amcc, apart
    return None


def test_matches_the_arpae_operator_on_a_shared_rng_stream():
    """The adapted MCC must reproduce its source EXACTLY, event for event.

    ARPAE/pic3d/mcc3d.py is oracle-matched against pic1d.mcc_ref; reproducing
    it bit-for-bit is what carries that validation across to this copy. Run
    with ion births OFF so both operators consume the RNG identically.
    """
    mods = _arpae_pic3d()
    if mods is None:
        pytest.skip("ARPAE/pic3d not present in this checkout")
    acs, amcc, apart = mods

    # ONE shared density, handed to both.  Letting each package compute its
    # own from (p, T) makes them differ by 3.7e-7 -- ARPAE's Loschmidt
    # constant carries a rounding slip in its 7th digit (see units.py) -- and
    # that is enough to move a few borderline null-collision draws and
    # desynchronise the streams. This gate is about the OPERATOR, so the
    # constant is taken out of it.
    from PIC_MCC_3D.MODULES.cross_sections import (CrossSectionSet as MyCS,
                                                   gas_file)
    n_ar = 2.4463132910655e23
    g = pure("Ar", 0.01, 300.0)
    mine = MCC3D(MyCS({"Ar": (gas_file("Ar"), n_ar)}), mixture=g)
    theirs = amcc.MCC3D(acs.CrossSectionSet({"Ar": (acs.gas_file("Ar"), n_ar)}))
    assert mine.cs.n_proc == theirs.cs.n_proc, "different process sets"
    assert np.array_equal(mine.cs.density, theirs.cs.density)

    # identical initial states
    rng0 = np.random.default_rng(21)
    n = 20000
    pos = rng0.random((n, 3)) * 1e-2
    vel = rng0.normal(0.0, np.sqrt(6.0 * EV / M_E), (n, 3))
    w = np.full(n, 3.0)

    p_mine = ParticleSpecies(Species.electron(), capacity=2 * n)
    p_mine.spawn(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1],
                 vel[:, 2], w)
    p_theirs = apart.Particles3D(capacity=2 * n)
    p_theirs.spawn(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1],
                   vel[:, 2], w)

    a = mine.collide_adaptive(p_mine, 2e-10, np.random.default_rng(77))
    b = theirs.collide_adaptive(p_theirs, 2e-10, np.random.default_rng(77))

    assert int(a["tally"].sum()) > 1000
    assert np.array_equal(a["tally"], b["tally"]), "process tallies diverged"
    assert np.allclose(a["energy_eV"], b["energy_eV"], rtol=0, atol=0), \
        "energy ledgers diverged"
    assert p_mine.n == p_theirs.n
    assert np.allclose(p_mine.vx[:p_mine.n], p_theirs.vx[:p_theirs.n],
                       rtol=0, atol=0), "velocities diverged"


# ------------------------------------------------- scattering-angle sampling
# Gap identified 2026-08-03 from the PIC-MCC verification literature (Parodi &
# Petronio, Phys. Plasmas 32, 013902 (2025); arXiv:2602.15130): comparing the
# EMPIRICAL angular CDF against the analytic one catches "coding errors that do
# not modify the expected outcome from the collision algorithm" -- a wrong
# angular distribution still conserves energy, still tallies correctly, and
# still passes every gate this suite had. These close that hole.
def _ks(sample, cdf) -> float:
    """Kolmogorov-Smirnov statistic of ``sample`` against an analytic CDF."""
    x = np.sort(np.asarray(sample, dtype=float))
    n = x.size
    f = cdf(x)
    return float(max(np.max(np.arange(1, n + 1) / n - f),
                     np.max(f - np.arange(0, n) / n)))


def test_isotropic_scattering_is_uniform_in_cos_chi():
    """Isotropic means uniform in cos(chi), NOT uniform in chi."""
    from PIC_MCC_3D.MODULES.mcc3d import sample_cos_chi
    rng = np.random.default_rng(0)
    n = 200000
    c = sample_cos_chi("isotropic", np.full(n, 5.0), rng)
    assert -1.0 <= c.min() and c.max() <= 1.0
    d = _ks(c, lambda x: (x + 1.0) / 2.0)
    assert d < 1.36 / np.sqrt(n) * 1.5, f"KS = {d:.4f}"


def test_screened_model_matches_its_analytic_cdf():
    """The screened-Coulomb sampler must invert its own distribution.

    Sampled by inverting  P(cos chi) with eta = e_screen/eps; the analytic CDF
    follows from that inversion, so a sign or algebra slip shows up here.
    """
    from PIC_MCC_3D.MODULES.mcc3d import sample_cos_chi
    rng = np.random.default_rng(1)
    n, eps, e_scr = 200000, 20.0, 27.21
    c = sample_cos_chi("screened", np.full(n, eps), rng, e_scr)
    eta = e_scr / eps
    # cos = (eta + 1 - (2 eta + 1) p) / (eta + 1 - p) DECREASES with p:
    # p = 0 -> cos = +1, p = 1 -> cos = -1.  So P(cos <= x) = 1 - p(x), not
    # p(x); using p(x) directly gives KS = 1 and looks like a broken sampler.
    def cdf(x):
        return 1.0 - (eta + 1.0) * (1.0 - x) / (2.0 * eta + 1.0 - x)
    d = _ks(c, cdf)
    assert d < 1.36 / np.sqrt(n) * 1.5, f"KS = {d:.4f}"
    assert c.mean() > 0.0, "screened scattering must be forward-peaked"


def test_surendra_model_matches_its_analytic_cdf():
    from PIC_MCC_3D.MODULES.mcc3d import sample_cos_chi
    rng = np.random.default_rng(2)
    n, eps = 200000, 15.0
    c = sample_cos_chi("surendra", np.full(n, eps), rng)
    # same inversion direction as the screened model: cos falls from +1 to -1
    # as p runs 0 -> 1, so the CDF is 1 - p(x).
    def cdf(x):
        return 1.0 - np.log((2.0 + eps - eps * x) / 2.0) / np.log(1.0 + eps)
    d = _ks(c, cdf)
    assert d < 1.36 / np.sqrt(n) * 1.5, f"KS = {d:.4f}"
    assert c.mean() > 0.0, "Surendra scattering must be forward-peaked"


def test_isotropic_directions_are_uniform_on_the_sphere():
    from PIC_MCC_3D.MODULES.mcc3d import _isotropic_dirs
    rng = np.random.default_rng(3)
    n = 200000
    d = _isotropic_dirs(rng, n)
    assert np.max(np.abs(np.linalg.norm(d, axis=1) - 1.0)) < 1e-12
    # <cos> = 0 and <cos^2> = 1/3 on a sphere, in every axis
    for k in range(3):
        assert abs(float(d[:, k].mean())) < 4.0 / np.sqrt(n)
        assert abs(float((d[:, k] ** 2).mean()) - 1.0 / 3.0) < 0.01
    phi = np.arctan2(d[:, 1], d[:, 0])
    assert _ks(phi, lambda x: (x + np.pi) / (2 * np.pi)) < 1.36 / np.sqrt(n) * 1.5


def test_scatter_rotation_produces_exactly_the_requested_angle():
    """The deflection applied must EQUAL the angle sampled.

    Independent of any distribution: it checks the basis construction in
    _scatter_dirs, including its near-x special case. A bug there would
    silently replace the intended angular model with a different one while
    every energy and tally gate kept passing.
    """
    from PIC_MCC_3D.MODULES.mcc3d import _scatter_dirs
    rng = np.random.default_rng(4)
    n = 50000
    v = rng.normal(size=(n, 3))
    vhat = v / np.linalg.norm(v, axis=1)[:, None]
    cos_chi = 1.0 - 2.0 * rng.random(n)
    eta = 2.0 * np.pi * rng.random(n)
    out = _scatter_dirs(vhat, cos_chi, eta)
    assert np.max(np.abs(np.linalg.norm(out, axis=1) - 1.0)) < 1e-12
    got = np.sum(out * vhat, axis=1)
    assert np.max(np.abs(got - cos_chi)) < 1e-12, "deflection != sampled angle"


def test_scatter_rotation_is_azimuthally_uniform_about_the_incident_axis():
    """Rotating by eta must sweep the full cone, not a preferred side."""
    from PIC_MCC_3D.MODULES.mcc3d import _scatter_dirs
    rng = np.random.default_rng(5)
    n = 100000
    vhat = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    cos_chi = np.full(n, 0.3)
    eta = 2.0 * np.pi * rng.random(n)
    out = _scatter_dirs(vhat, cos_chi, eta)
    phi = np.arctan2(out[:, 1], out[:, 0])
    assert _ks(phi, lambda x: (x + np.pi) / (2 * np.pi)) < 1.36 / np.sqrt(n) * 1.5
