"""Gates for the LXCat loader.

These check the DATA PATH and the tabulation, not the physics of the cross
sections themselves — that provenance belongs to the LXCat files and to the
pic1d/LoKI-MC lineage this module was adapted from.
"""

import numpy as np
import pytest

from PIC_MCC_3D.MODULES.cross_sections import (CrossSectionSet, gas_file,
                                               lxcat_roots, mixture_set,
                                               parse_lxcat)
from PIC_MCC_3D.MODULES.gas import GasMixture, pure

pytestmark = pytest.mark.skipif(
    not any(r.is_dir() for r in lxcat_roots()),
    reason="no LXCat data tree reachable from this checkout")


def test_data_root_is_found_by_searching_upward():
    """The resolver must work from a worktree, not just the main checkout.

    Counting parent directories (``parents[2]``) silently resolved to a
    non-existent path inside a git worktree, where the package sits two levels
    deeper. This gate pins the searching behaviour that replaced it.
    """
    roots = [r for r in lxcat_roots() if r.is_dir()]
    assert roots, "no data root resolved"
    assert gas_file("Ar").is_file()


def test_every_mapped_gas_resolves():
    for name in ("Ar", "Kr", "F2", "He"):
        assert gas_file(name).is_file()


def test_unknown_gas_raises():
    with pytest.raises(KeyError, match="no LXCat file mapped"):
        gas_file("Xe")


def test_parse_returns_physical_processes():
    procs = parse_lxcat(gas_file("Ar"))
    assert procs
    kinds = {p.kind for p in procs}
    assert "ionization" in kinds
    for p in procs:
        assert p.energy.size == p.sigma.size > 1
        assert np.all(np.diff(p.energy) > 0), f"{p.kind} energies not increasing"
        assert np.all(p.sigma >= 0.0)
        if p.kind in ("excitation", "ionization"):
            assert p.threshold > 0.0
        if p.kind in ("elastic", "effective"):
            assert 0.0 < p.mass_ratio < 1e-2


def test_argon_ionization_threshold_is_physical():
    """Ar first ionization sits at 15.76 eV; anything far off means a bad parse."""
    procs = parse_lxcat(gas_file("Ar"))
    ion = [p for p in procs if p.kind == "ionization"]
    assert ion
    assert 15.0 < min(p.threshold for p in ion) < 16.5


def test_effective_is_dropped_and_elastic_exists():
    """'effective' would double-count with elastic + inelastic."""
    cs = mixture_set(pure("Ar", 1.0))
    assert "effective" not in cs.kinds
    assert "elastic" in cs.kinds


def test_collision_frequency_is_in_the_expected_band():
    """nu_m(1 eV, 1 atm Ar) -- the same sanity band ARPAE's system map probes."""
    cs = mixture_set(pure("Ar", 1.0, 300.0))
    nu = float(cs.nu_momentum(np.array([1.0]))[0])
    assert 5e10 < nu < 5e11, nu


def test_frequency_scales_linearly_with_density():
    a = mixture_set(pure("Ar", 0.5)).nu_total(np.array([5.0]))[0]
    b = mixture_set(pure("Ar", 1.0)).nu_total(np.array([5.0]))[0]
    assert abs(b / a - 2.0) < 1e-9


def test_majorant_bounds_the_total_rate():
    """The null-collision majorant must never be exceeded, or the scheme is wrong."""
    cs = mixture_set(pure("Ar", 1.0))
    eps = np.geomspace(1e-3, 200.0, 4000)
    nu = cs.nu_total(eps)
    assert cs.nu_majorant(200.0) >= nu.max() * (1.0 - 1e-12)
    # and it is monotone non-decreasing in the energy ceiling
    ceilings = [1.0, 10.0, 100.0, 500.0]
    vals = [cs.nu_majorant(c) for c in ceilings]
    assert all(b >= a - 1e-30 for a, b in zip(vals, vals[1:]))


def test_mixture_carries_every_component():
    mix = GasMixture({"Ar": 0.9495, "Kr": 0.05, "F2": 0.0005}, 1.0)
    cs = mixture_set(mix)
    assert set(cs.species_density) == {"Ar", "Kr", "F2"}
    assert "attachment" in cs.kinds          # F2 brings dissociative attachment
    for name, dens in cs.species_density.items():
        assert abs(dens - mix.density_of(name)) < 1e-6 * dens


def test_sigma_is_zero_outside_the_tabulated_range():
    procs = parse_lxcat(gas_file("Ar"))
    p = procs[0]
    # NOT `energy[0] * 1e-6`: Ar's elastic table starts at exactly 0 eV, so
    # that scales to 0 -- which is ON the first point, not below it, and
    # correctly returns sigma[0].
    assert p.sigma_at(np.array([-1.0]))[0] == 0.0
    assert p.sigma_at(np.array([p.energy[-1] * 1e3]))[0] == 0.0


def test_empty_mixture_rejected():
    class _Empty:
        densities = {}
    with pytest.raises(ValueError, match="no components"):
        mixture_set(_Empty())


def test_describe_is_informative():
    cs = mixture_set(pure("Ar", 1.0))
    d = cs.describe()
    assert "processes" in d and "ionization" in d


def test_set_can_be_built_directly_from_paths():
    cs = CrossSectionSet({"Ar": (gas_file("Ar"), 1.0e24)})
    assert cs.n_proc > 0
    assert abs(cs.species_density["Ar"] - 1.0e24) < 1e9
