"""Tests for the reaction network: parsing, unit conversion, layer routing, Σ/O."""

from __future__ import annotations

import unittest
import warnings
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
from reactions import Layer, RateForm, ReactionNetwork, load_reaction_network, parse_side
from reactions.units import to_si
from particles import excimer_system_registry

KRF_TOML = _ROOT / "data" / "reactions" / "krf.toml"
ARF_TOML = _ROOT / "data" / "reactions" / "arf.toml"
ARF_CHEN_TOML = _ROOT / "data" / "reactions" / "arf_chen_2025.toml"


class TestArfNetwork(unittest.TestCase):
    """The Chen-derived ArF network loads and compiles against the registry."""

    def test_arf_compiles(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # orphan-charge warnings (Xe+/Ar2F* have no loss in Table III)
            net = load_reaction_network(ARF_TOML, excimer_system_registry()).compile()
        total = len(net.fluid) + len(net.kinetic) + len(net.photon)
        # Chen R1-R41 plus active cited additions R21b and R29c; R29b is disabled.
        self.assertEqual(total, 43)
        # the resolved ArF species are registered
        for name in ("Ne**", "F2(a3)", "F2(A1)", "F2+", "ArF*", "Ar2F*"):
            self.assertIn(name, net.index)

    def test_arf_photon_accessors(self) -> None:
        from reactions import photo_absorbers, stimulated_channel
        reg = excimer_system_registry()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            net = load_reaction_network(ARF_TOML, reg).compile()
        # R38 stimulated emission: ArF* -> Ar + F @ hv193, σ = 4e-20 m².
        upper, sigma, lower, mode = stimulated_channel(net, reg)
        self.assertEqual(upper, "ArF*")
        self.assertAlmostEqual(sigma, 4.0e-20, places=24)
        self.assertEqual(set(lower), {(1, "Ar"), (1, "F")})
        self.assertEqual(mode, "hv193")
        # R40/R41 UV absorbers: Ar** (1e-22) and F- (1e-21) at hv193.
        absorbers = {sp: sig for sp, sig, m in photo_absorbers(net, reg) if m == "hv193"}
        self.assertAlmostEqual(absorbers["Ar**"], 1.0e-22, places=26)
        self.assertAlmostEqual(absorbers["F-"], 1.0e-21, places=25)

    def test_te_dependent_energy_metadata_survives_load_and_compile(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            compiled = load_reaction_network(ARF_TOML, excimer_system_registry()).compile()
        by_id = {
            reaction.reaction_id: reaction
            for reaction in (*compiled.fluid, *compiled.kinetic, *compiled.photon)
        }
        self.assertEqual(by_id["R9"].params["energy_loss_eV"], 2.63)
        self.assertEqual(by_id["R18"].params["energy_loss_eV"], 3.4)

    def test_chen_reproduction_overlay_is_separate_from_production(self) -> None:
        reg = excimer_system_registry()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            production = load_reaction_network(ARF_TOML, reg).compile()
            reproduction = load_reaction_network(ARF_CHEN_TOML, reg).compile()
        production_ids = {
            reaction.reaction_id
            for reaction in (*production.fluid, *production.kinetic, *production.photon)
        }
        reproduction_rows = {
            reaction.reaction_id: reaction
            for reaction in (*reproduction.fluid, *reproduction.kinetic, *reproduction.photon)
        }
        self.assertEqual(len(production_ids), 43)
        self.assertEqual(len(reproduction_rows), 41)
        self.assertNotIn("R21b", reproduction_rows)
        self.assertNotIn("R29c", reproduction_rows)
        self.assertEqual(reproduction_rows["R26"].product_names, ("Ne+", "Ar"))


class TestParsing(unittest.TestCase):
    def test_string_with_multiplicity_and_charge(self) -> None:
        self.assertEqual(parse_side("Kr+ + 2e"), ((1, "Kr+"), (2, "e-")))

    def test_tuple_form(self) -> None:
        self.assertEqual(parse_side(("Kr+", "2e")), ((1, "Kr+"), (2, "e-")))

    def test_charge_plus_not_a_separator(self) -> None:
        # The '+' in Kr+ must survive; only the whitespace-delimited '+' splits.
        self.assertEqual(parse_side("Kr+ + F-"), ((1, "Kr+"), (1, "F-")))

    def test_aliases_and_spaced_coefficient(self) -> None:
        self.assertEqual(parse_side("2 Kr + hv"), ((2, "Kr"), (1, "hv248")))


class TestUnits(unittest.TestCase):
    def test_two_body(self) -> None:
        self.assertAlmostEqual(to_si(1.0, "cm3/s"), 1e-6)

    def test_whitespace_is_normalized(self) -> None:
        self.assertAlmostEqual(to_si(1.0, "  cm3   eV^0.5/s  "), 1e-6)

    def test_three_body(self) -> None:
        self.assertAlmostEqual(to_si(1.0, "cm6/s"), 1e-12)

    def test_cross_section(self) -> None:
        self.assertAlmostEqual(to_si(1.0, "cm2"), 1e-4)

    def test_unknown_unit_raises(self) -> None:
        with self.assertRaises(KeyError):
            to_si(1.0, "furlong/fortnight")


class TestRoutingAndMatrices(unittest.TestCase):
    def test_layer_routing(self) -> None:
        net = ReactionNetwork()
        net.add("Kr* + F2", "KrF* + F", r=8.1e-10)  # all neutral -> FLUID
        net.add("Kr + e", "Kr+ + 2e", sigma="x")  # charged reactant -> KINETIC
        net.add("KrF* + hv248", "Kr + F + 2 hv248", rate_form=RateForm.STIMULATED,
                params={"sigma_eff_cm2": 5.354e-17})  # photon reactant -> PHOTON
        c = net.compile()
        self.assertEqual(len(c.fluid), 1)
        self.assertEqual(len(c.kinetic), 1)
        self.assertEqual(len(c.photon), 1)

    def test_fluid_stoichiometry_and_orders(self) -> None:
        net = ReactionNetwork()
        net.add("Kr* + F2", "KrF* + F", r=8.1e-10)
        c = net.compile()
        i = c.index
        # single fluid reaction, column 0
        self.assertEqual(c.orders[i["Kr*"], 0], 1)
        self.assertEqual(c.orders[i["F2"], 0], 1)
        self.assertEqual(c.orders[i["KrF*"], 0], 0)
        self.assertEqual(c.stoichiometry[i["Kr*"], 0], -1)
        self.assertEqual(c.stoichiometry[i["F2"], 0], -1)
        self.assertEqual(c.stoichiometry[i["KrF*"], 0], +1)
        self.assertEqual(c.stoichiometry[i["F"], 0], +1)
        # SI conversion applied
        self.assertAlmostEqual(c.rate_si[0], 8.1e-10 * 1e-6)

    def test_three_body_order(self) -> None:
        net = ReactionNetwork()
        net.add("Kr* + 2 Kr", "Kr2* + Kr", r=4.5e-32, r_units="cm6/s")
        c = net.compile()
        i = c.index
        self.assertEqual(c.orders[i["Kr"], 0], 2)  # third-body order 2
        self.assertEqual(c.stoichiometry[i["Kr"], 0], -1)  # net: 2 consumed, 1 produced
        self.assertEqual(c.stoichiometry[i["Kr2*"], 0], +1)

    def test_radiative_rate_is_inverse_tau(self) -> None:
        net = ReactionNetwork()
        net.add("KrF*", "Kr + F + hv248", tau_s=25e-9, photon_mode="hv248")
        c = net.compile()
        self.assertAlmostEqual(c.rate_si[0], 1.0 / 25e-9)

    def test_orphan_charge_warns(self) -> None:
        from reactions import network as _netmod
        _netmod._warned_orphans.clear()  # per-process dedup would swallow the re-warn
        net = ReactionNetwork()
        net.add("Ar + e", "Ar+ + 2e", sigma="x")  # Ar+ produced, never consumed
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            net.compile()
        self.assertTrue(any("Ar+" in str(w.message) for w in caught))


class TestTomlDataset(unittest.TestCase):
    def test_photo_row_requires_cross_section(self) -> None:
        from reactions.data import _add_row

        row = {
            "id": "photo-with-rate",
            "reactants": ["Kr**", "hv248"],
            "products": ["Kr+", "e-"],
            "rate_form": "photo",
            "k": 1.0,
            "k_units": "1/s",
        }
        with self.assertRaisesRegex(ValueError, "requires sigma_m2"):
            _add_row(ReactionNetwork(), row)

    def test_load_and_compile_krf(self) -> None:
        net = load_reaction_network(KRF_TOML)
        c = net.compile()
        self.assertEqual(len(c.fluid) + len(c.kinetic) + len(c.photon), 45)
        # Pin representative channels in every execution layer.
        kinetic_ids = {r.reaction_id for r in c.kinetic}
        self.assertLessEqual({"R1", "R5", "R9", "R12", "R15", "R30", "R33"}, kinetic_ids)
        photon_ids = {r.reaction_id for r in c.photon}
        self.assertEqual(photon_ids, {"R39", "R43", "R44", "R45"})
        fluid_ids = {r.reaction_id for r in c.fluid}
        self.assertLessEqual({"R16", "R22", "R28", "R34", "R38", "R40", "R42"}, fluid_ids)
        n_fluid = len(c.fluid)
        # matrices are (S, R_fluid) and integer-valued
        self.assertEqual(c.stoichiometry.shape, (c.n_species, n_fluid))
        self.assertEqual(c.orders.shape, (c.n_species, n_fluid))
        self.assertTrue(np.all(c.orders >= 0))
        # every fluid rate constant converted to SI (finite, > 0)
        self.assertTrue(np.all(np.isfinite(c.rate_si)))

    def test_krf_photon_channels(self) -> None:
        from reactions import photo_absorbers, stimulated_channel

        registry = excimer_system_registry()
        compiled = load_reaction_network(KRF_TOML, registry).compile()
        upper, sigma, lower, mode = stimulated_channel(compiled, registry)
        self.assertEqual(upper, "KrF*")
        self.assertAlmostEqual(sigma, 2.6e-20, places=24)
        self.assertEqual(set(lower), {(1, "Kr"), (1, "F")})
        self.assertEqual(mode, "hv248")
        absorbers = {
            species: cross_section
            for species, cross_section, photon_mode in photo_absorbers(compiled, registry)
            if photon_mode == "hv248"
        }
        self.assertEqual(set(absorbers), {"Kr**", "F-", "F2"})
        self.assertAlmostEqual(absorbers["Kr**"], 6.0e-22, places=26)
        self.assertAlmostEqual(absorbers["F-"], 5.6e-22, places=26)
        self.assertAlmostEqual(absorbers["F2"], 1.48e-24, places=28)


class TestNetworkQueries(unittest.TestCase):
    def test_radiative_channels_query(self) -> None:
        from particles.species import excimer_system_registry
        from reactions import radiative_channels

        registry = excimer_system_registry()
        net = ReactionNetwork()
        net.add("KrF*", "Kr + F + hv248", tau_s=25e-9, photon_mode="hv248")
        net.add("Kr* + F2", "KrF* + F", r=8.1e-10)  # non-radiative: excluded
        chans = radiative_channels(net.compile(), registry)
        self.assertEqual(len(chans), 1)
        upper, tau, lower = chans[0]
        self.assertEqual(upper, "KrF*")
        self.assertAlmostEqual(tau, 25e-9)
        self.assertEqual(lower, ((1, "Kr"), (1, "F")))  # photon product excluded

    def test_stimulated_sigma_query(self) -> None:
        from reactions import RateForm, stimulated_sigma

        net = ReactionNetwork()
        net.add("KrF* + hv248", "Kr + F + 2 hv248", rate_form=RateForm.STIMULATED,
                params={"sigma_eff_cm2": 5.354e-17})
        compiled = net.compile()
        self.assertAlmostEqual(stimulated_sigma(compiled), 5.354e-17 * 1e-4)
        empty = ReactionNetwork()
        empty.add("Kr* + F2", "KrF* + F", r=8.1e-10)
        self.assertEqual(stimulated_sigma(empty.compile()), 0.0)


if __name__ == "__main__":
    unittest.main()
