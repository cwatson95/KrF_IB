"""species.py gates: registry completeness, energy-reference consistency."""

import numpy as np
import pytest

from REBUILD.MODULES.provenance import TIER_ORDER
from REBUILD.MODULES.species import (ARF_MODE_SPECIES, KRF_MODE_SPECIES,
                                     REGISTRY, StateLayout, atoms_of,
                                     charge_of, e_form_of)


def test_every_species_carries_full_data():
    for name, sp in REGISTRY.items():
        assert sp.atoms, name
        assert sp.mass_amu > 0, name
        assert sp.tier in TIER_ORDER, name
        assert sp.cite, name


def test_mode_lists_resolve():
    for lst in (KRF_MODE_SPECIES, ARF_MODE_SPECIES):
        for s in lst:
            assert s in REGISTRY, s


def test_charges_match_names():
    for s in REGISTRY:
        q = charge_of(s)
        if s.endswith("+"):
            assert q == +1, s
        elif s.endswith("-"):
            assert q == -1, s
        else:
            assert q == 0, s
    assert charge_of("e") == -1


def test_excimer_energy_reference_consistency():
    """E_form(excimer) must cover fragments + emitted quantum, else the
    stimulated/spontaneous bookkeeping would mint energy (the rule stated
    in the registry)."""
    hv = {"KrF*": 5.00, "ArF*": 6.42, "Kr2F*": 3.0, "Ar2F*": 4.35}
    frags = {"KrF*": ("Kr", "F"), "ArF*": ("Ar", "F"),
             "Kr2F*": ("Kr", "Kr", "F"), "Ar2F*": ("Ar", "Ar", "F")}
    for exc, quantum in hv.items():
        need = quantum + sum(e_form_of(f) for f in frags[exc])
        assert e_form_of(exc) >= need - 1e-9, exc


def test_f_minus_formation_energy_is_negative():
    # F- = E_form(F) - EA(F) = 0.80 - 3.40
    assert abs(e_form_of("F-") - (-2.60)) < 1e-12


def test_layout_indexing_and_quasineutrality():
    lay = StateLayout(KRF_MODE_SPECIES, n_cells=3, n_ledger=4)
    assert lay.size == 3 * (len(KRF_MODE_SPECIES) + 2) + 1 + 4
    y = np.zeros(lay.size)
    dens = lay.density_block(y)
    dens[lay.index["Kr+"], :] = 5e18
    dens[lay.index["Ar+"], :] = 2e18
    dens[lay.index["F-"], :] = 3e18
    ne = lay.electron_density(dens)
    assert np.allclose(ne, 4e18)
    dens[lay.index["F-"], 1] = 9e18          # anion excess -> floor
    ne = lay.electron_density(dens)
    assert ne[1] == 1.0


def test_atoms_of_photon_and_electron_are_empty():
    assert atoms_of("e") == {}
    assert atoms_of("hv193") == {}
    assert atoms_of("Kr2F*") == {"Kr": 2, "F": 1}


def test_unknown_species_raises():
    with pytest.raises(KeyError):
        atoms_of("Xy*")
