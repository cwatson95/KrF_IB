"""Gates for the e-beam module (missing module 9) and the preionization seed.

The scoring gates here are the strongest first-principles checks in the
package, and two of them are SCORED ANCHORS, not pass/fail physics: the
W-value comparison against the literature is reported with its miss (PLAN
(26) -- anchors score, they never steer). Nothing in this module is tuned to
reproduce 26.4 eV, and a gate that demanded it would be a tuning instruction
in disguise.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.beam3d import (LITERATURE_W_EV, BeamInjector, BeamSpec,
                                       cascade_yields, csda_range,
                                       default_cascade_energy, deposition_profile,
                                       differential_W, sample_opal_peterson,
                                       seed_preionization, spencer_fano_yields,
                                       stopping_power)
from PIC_MCC_3D.MODULES.cross_sections import lxcat_roots, mixture_set
from PIC_MCC_3D.MODULES.gas import GasMixture, pure
from PIC_MCC_3D.MODULES.mcc3d import MCC3D
from PIC_MCC_3D.MODULES.species3d import (ION_MASSES_AMU, ParticleSpecies,
                                          Species)
from PIC_MCC_3D.MODULES.units import EV

_HAVE_DATA = any(r.is_dir() for r in lxcat_roots())
needs_data = pytest.mark.skipif(not _HAVE_DATA, reason="no LXCat data tree")

# NIST ESTAR collisional mass stopping powers at ~600 keV [MeV cm^2/g], as
# carried by unified/ebeam.py. An INDEPENDENT tabulation: our Bethe/Moller
# evaluation must reproduce them without having seen them.
ESTAR_600KEV = {"He": 1.94, "Ar": 1.52, "Kr": 1.34, "F2": 1.63}


# ------------------------------------------------------------ stopping power
@pytest.mark.parametrize("gas,estar", sorted(ESTAR_600KEV.items()))
def test_bethe_stopping_power_matches_nist_estar(gas, estar):
    """First-principles Bethe/Moller vs the NIST ESTAR tabulation, at 600 keV."""
    m = pure(gas, 1.0, 300.0)
    rho_g_cm3 = m.total_density_m3 * m.mass_of(gas) * 1e3 / 1e6
    s_mass = stopping_power(m, 600e3) / 1e6 / 100.0 / rho_g_cm3
    assert abs(s_mass / estar - 1.0) < 0.10, f"{gas}: {s_mass:.3f} vs {estar}"


def test_stopping_power_has_the_right_energy_dependence():
    """S falls steeply from keV to the ~1 MeV minimum -- the shape sets deposition."""
    m = pure("Ar", 1.0)
    s = [stopping_power(m, e) for e in (1e4, 1e5, 6e5)]
    assert s[0] > 4.0 * s[1] > 0.0
    assert s[1] > s[2] > 0.0


def test_stopping_power_scales_with_density():
    a = stopping_power(pure("Ar", 0.5), 600e3)
    b = stopping_power(pure("Ar", 1.0), 600e3)
    assert abs(b / a - 2.0) < 1e-9


def test_stopping_power_of_a_mixture_is_additive():
    mix = GasMixture({"Ar": 0.5, "Kr": 0.5}, 1.0)
    want = (stopping_power(pure("Ar", 0.5), 600e3)
            + stopping_power(pure("Kr", 0.5), 600e3))
    assert abs(stopping_power(mix, 600e3) / want - 1.0) < 1e-12


def test_csda_range_is_physical():
    """A 600 keV electron ranges out over ~1-2 m of 1 atm argon."""
    r = csda_range(pure("Ar", 1.0), 600e3)
    assert 0.5 < r < 4.0, r
    # and the range shortens at higher pressure, in proportion
    r2 = csda_range(pure("Ar", 2.0), 600e3)
    assert abs(r2 / r - 0.5) < 0.02


def test_unknown_gas_has_no_stopping_data():
    m = pure("Ar", 1.0)
    m._densities["Xe"] = 1e24                      # force an unknown component
    with pytest.raises(KeyError, match="ICRU"):
        stopping_power(m, 6e5)


# ------------------------------------------------------- deposition profile
def test_deposition_profile_conserves_and_orders_energy():
    m = pure("Ar", 1.0)
    spec = BeamSpec(voltage_kV=600.0)
    prof = deposition_profile(m, spec, h=1e-3, n_cells=200)
    assert prof.sum() <= 1.0 + 1e-12
    assert np.all(prof >= 0.0)
    # a 600 keV beam barely attenuates over 20 cm, so deposition is near-flat
    # but must still be monotonically non-decreasing as the electron slows
    nz = prof[prof > 0]
    assert nz[-1] >= nz[0]


def test_deposition_direction_flips():
    m = pure("Ar", 1.0)
    a = deposition_profile(m, BeamSpec(from_low=True), 1e-2, 40)
    b = deposition_profile(m, BeamSpec(from_low=False), 1e-2, 40)
    assert np.allclose(a, b[::-1])


def test_beam_envelope_shape():
    spec = BeamSpec(pulse_ns=50.0, rise_ns=5.0)
    assert spec.envelope(0.0) == 0.0
    assert abs(spec.envelope(2.5e-9) - 0.5) < 1e-12       # mid-rise
    assert spec.envelope(30e-9) == 1.0                     # flat top
    assert spec.envelope(200e-9) == 0.0                    # after the pulse
    assert spec.power_areal_W_m2(30e-9) > 0.0


# --------------------------------------------------------- the energy ceiling
@needs_data
def test_cross_section_sets_report_their_validity_ceiling():
    """The ceiling must be real and BELOW the lookup table's span.

    This is the trap the module exists to avoid: sigma_at returns ZERO above
    the last tabulated point, so an electron launched past the ceiling
    collides with nothing and the cascade reports zero ionizations rather
    than an error.
    """
    for gas, lo, hi in (("Ar", 900.0, 1100.0), ("Kr", 900.0, 1000.0),
                        ("He", 900.0, 1100.0)):
        c = mixture_set(pure(gas, 1.0)).energy_ceiling_eV
        assert lo < c < hi, f"{gas} ceiling {c}"


@needs_data
def test_cascades_refuse_to_start_above_the_ceiling():
    cs = mixture_set(pure("Ar", 1.0))
    too_high = cs.energy_ceiling_eV * 1.5
    with pytest.raises(ValueError, match="validity ceiling"):
        spencer_fano_yields(cs, too_high)
    with pytest.raises(ValueError, match="validity ceiling"):
        cascade_yields(MCC3D(cs, mixture=pure("Ar", 1.0)), too_high)


@needs_data
def test_default_cascade_energy_is_inside_the_clean_window():
    cs = mixture_set(pure("Ar", 1.0))
    e = default_cascade_energy(cs)
    assert 0.2 * cs.energy_ceiling_eV < e < 0.7 * cs.energy_ceiling_eV


# ------------------------------------------------------------- the W-value
@needs_data
def test_the_two_W_routes_agree():
    """THE gate: deterministic recursion vs Monte-Carlo through the real MCC.

    Both are priced from the same LXCat data, so agreement tests the SAMPLING
    -- that the MCC selects processes with the probabilities the cross
    sections imply. A disagreement here is a bug; a common offset from the
    literature is not (see the scored anchor below).
    """
    m = pure("Ar", 1.0)
    cs = mixture_set(m)
    E0 = default_cascade_energy(cs)
    sf = spencer_fano_yields(cs, E0)
    mc = cascade_yields(MCC3D(cs, mixture=m), E0, n_electrons=120,
                        rng=np.random.default_rng(2))
    assert mc["converged"]
    assert abs(mc["W_eV"] / sf["W_eV"] - 1.0) < 0.05, (mc["W_eV"], sf["W_eV"])
    assert abs(mc["U_ion"] / sf["U_ion"] - 1.0) < 0.05


@needs_data
def test_W_falls_toward_its_asymptote_then_is_spoiled_by_the_table_edge():
    """The energy dependence that makes differential_W the right quantity."""
    cs = mixture_set(pure("Ar", 1.0))
    w = [spencer_fano_yields(cs, e)["W_eV"] for e in (100, 200, 400, 600)]
    assert w[0] > w[1] > w[2] > w[3]            # falling toward the asymptote
    near_edge = spencer_fano_yields(cs, 0.95 * cs.energy_ceiling_eV)["W_eV"]
    assert near_edge > w[3]                      # and bending back up at the edge


@needs_data
@pytest.mark.parametrize("gas,tol", [("Ar", 0.15), ("Kr", 0.15)])
def test_W_value_scored_against_the_literature(gas, tol):
    """SCORED ANCHOR, not a tuning target.

    Nothing in this module is fitted to the literature W. Argon and krypton
    come out ~5 % and ~8 % high, which is the expected direction for a
    Spencer-Fano treatment with equal ionization sharing on a
    momentum-transfer LXCat basis. The band is deliberately loose enough that
    passing it is not evidence of agreement to better than it states.
    """
    W = differential_W(mixture_set(pure(gas, 1.0)))
    assert abs(W / LITERATURE_W_EV[gas] - 1.0) < tol, W


@needs_data
def test_helium_W_is_a_reported_miss():
    """He comes out ~33 % HIGH. Recorded as a miss, not hidden or tuned away.

    Pinned so the number cannot drift silently: if a future change moves it,
    that is a result worth noticing in either direction. Candidate causes are
    in the changelog -- equal-sharing ionization, and the He set's channel
    coverage.
    """
    W = differential_W(mixture_set(pure("He", 1.0)))
    ratio = W / LITERATURE_W_EV["He"]
    assert 1.2 < ratio < 1.5, f"He W ratio moved to {ratio:.3f}"


# -------------------------------------------------- secondary-electron spectrum
def test_opal_peterson_spectrum_is_sane():
    rng = np.random.default_rng(1)
    eps = sample_opal_peterson(rng, 200000, "Ar", 200.0)
    assert np.all(eps >= 0.0) and np.all(eps <= 200.0 + 1e-9)
    # median of 1/(1+(e/Ebar)^2) truncated at eps_max, computed from the CDF
    e_bar = 10.0
    want = e_bar * np.tan(0.5 * np.arctan(200.0 / e_bar))
    assert abs(np.median(eps) / want - 1.0) < 0.05


# ------------------------------------------------------------- the injector
@needs_data
def test_injector_conserves_the_energy_it_deposits():
    """Kinetic mode: the secondaries carry exactly the deposited energy."""
    from PIC_MCC_3D.MODULES.es_sim import ESSim
    from PIC_MCC_3D.MODULES.poisson3d import PERIODIC_XY_DIRICHLET_Z
    from PIC_MCC_3D.MODULES.species3d import SpeciesSet
    m = pure("Ar", 1.0)
    electrons = ParticleSpecies(Species.electron(), capacity=1 << 16)
    sim = ESSim(6, 6, 12, 1e-3, SpeciesSet([electrons]), 1e-12,
                bc=PERIODIC_XY_DIRICHLET_Z)
    sim.time = 20e-9                                  # inside the flat top
    inj = BeamInjector(m, BeamSpec(), electrons, rng=np.random.default_rng(0))
    rec = inj.inject(sim)
    assert rec["n_macro"] > 0
    assert electrons.n == rec["n_macro"]
    ke = electrons.kinetic_energy_J()
    # against the DEPOSITED energy, not what entered: a 600 keV beam ranges
    # metres, so a 12 mm cell stops a tiny fraction of the incident power
    assert abs(ke / rec["energy_deposited_J"] - 1.0) < 1e-9
    assert rec["energy_deposited_J"] < 0.01 * rec["energy_in_J"]


@needs_data
def test_injector_is_silent_outside_the_pulse():
    from PIC_MCC_3D.MODULES.es_sim import ESSim
    from PIC_MCC_3D.MODULES.species3d import SpeciesSet
    m = pure("Ar", 1.0)
    electrons = ParticleSpecies(Species.electron(), capacity=256)
    sim = ESSim(6, 6, 12, 1e-3, SpeciesSet([electrons]), 1e-12)
    sim.time = 500e-9                                 # long after the pulse
    inj = BeamInjector(m, BeamSpec(), electrons)
    assert inj.inject(sim)["n_macro"] == 0
    assert electrons.n == 0


def test_w_mode_refuses_to_default_the_W_value():
    """pair_mode='w' makes W an INPUT, so it must be stated explicitly."""
    m = pure("Ar", 1.0)
    e = ParticleSpecies(Species.electron())
    with pytest.raises(ValueError, match="W-value an INPUT"):
        BeamInjector(m, BeamSpec(), e, pair_mode="w")
    with pytest.raises(ValueError, match="unknown pair_mode"):
        BeamInjector(m, BeamSpec(), e, pair_mode="magic")


# --------------------------------------------------- the preionization seed
def test_preionization_seed_is_quasi_neutral_and_co_located():
    """The surface-discharge stand-in: matched, co-located electron/ion pairs.

    Co-location matters. Seeding the two populations independently would
    leave a sampling-noise charge imbalance in every cell and launch a
    spurious space-charge field the moment the run starts.
    """
    rng = np.random.default_rng(5)
    e = ParticleSpecies(Species.electron(), capacity=1 << 14)
    i = ParticleSpecies(Species.ion("Ar+", ION_MASSES_AMU["Ar+"]),
                        capacity=1 << 14)
    net = seed_preionization(e, i, 8000, [0, 0, 0], [1e-3] * 3, 1e16, 1.0, rng)
    assert e.n == i.n == 8000
    assert abs(net) < 1e-30
    assert np.array_equal(e.x[:e.n], i.x[:i.n])       # exactly co-located
    assert np.array_equal(e.z[:e.n], i.z[:i.n])
    # the seeded density is what was asked for
    assert abs(e.real_count() / (1e16 * 1e-9) - 1.0) < 1e-12
    # electrons are hot, ions are cold -- and the ions are far slower
    assert e.temperature_eV() > 10.0 * i.temperature_eV()
