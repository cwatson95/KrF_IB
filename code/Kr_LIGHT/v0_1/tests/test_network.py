"""The pure-Kr network as data: its rows, balance, routing, bands, and the
drift gates against SPECK/REBUILD (LIGHT v0_5's tests on krypton)."""

import importlib
import os
import sys
from pathlib import Path

import pytest

from v0_1 import network
from v0_1.species import SPECIES, e_form


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


PURE_KR = ["G1", "G2", "G6", "C1a", "C1b", "C1c", "C2", "K15", "R11p",
           "AK8_Kr", "KR2_RAD", "PEN_KR2", "KPS", "PRAD", "K76", "K76r",
           "KSTEP_S", "KSTEP_P"]


def test_the_network_is_krypton_and_electrons_alone():
    """18 rows, v0_5's krypton rows in v0_5's order, every species krypton
    or the electron; the optional group adds KSTEP and SUPER."""
    rows = network.build()
    assert [r.rid for r in rows] == PURE_KR
    assert set(SPECIES) == {"Kr", "Kr*", "Kr**", "Kr+", "Kr2*", "Kr2+"}
    for r in rows + network.optional_rows():
        for s in r.reactants + r.products:
            assert s == "e" or s.startswith("Kr"), (r.rid, s)
    every = network.build(("all",))
    assert [r.rid for r in every] == PURE_KR + ["KSTEP", "SUPER"]
    assert network.GROUPS == network.BASIC_GROUPS == ("electron",)


def test_every_row_balances_and_is_cited():
    rows = network.build(("all",))          # build() runs check_row
    assert len(rows) == 20 and all(r.cite for r in rows)


def test_unknown_or_fluorine_group_refused():
    for g in ("nonsense", "f2e", "kr2", "ions", "quench", "bc"):
        with pytest.raises(network.NetworkError):
            network.build((g,))


def test_energy_routing_examples():
    by = {r.rid: r for r in network.build(("all",))}
    route = network.energy_routing
    assert route(by["G2"])["to_e"] == approx(2 * 9.915 - 13.996)
    assert route(by["C1a"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["C2"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["SUPER"])["to_e"] == approx(9.915)
    # the 147-nm photon leaves 0.77 eV of ground-state repulsion with the
    # ledger's Kr2* at 9.2 eV
    kr2 = route(by["KR2_RAD"])
    assert (kr2["rad"], kr2["heat"]) == (8.43, approx(9.2 - 8.43, rel=1e-9))
    # three-body association heats the gas by the excimer's well
    assert route(by["AK8_Kr"])["heat"] == approx(9.915 - 9.2, rel=1e-9)
    # Kr2* + Kr2* Penning gives its excess to the new electron
    assert route(by["PEN_KR2"])["to_e"] == approx(2 * 9.2 - (13.996 - 1.15),
                                                  rel=1e-12)
    # dissociative recombination captures an electron
    assert route(by["R11p"])["capture"] == 1


# ---------------------------------------------------------------------------
# drift gate: every value copied from SPECK/REBUILD must still equal it
# ---------------------------------------------------------------------------
def _find_speck():
    env = os.environ.get("LIGHT_SPECK_DIR")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for base in here.parents:
        cand = base / "SPECK" / "REBUILD" / "MODULES" / "network.py"
        if cand.is_file():
            return base / "SPECK"
    return None


@pytest.fixture(scope="module")
def rebuild_rows():
    speck = _find_speck()
    if speck is None:
        pytest.skip("SPECK/REBUILD not found above LIGHT_Kr; set "
                    "LIGHT_SPECK_DIR to run the drift gate")
    sys.path.insert(0, str(speck))
    try:
        rn = importlib.import_module("REBUILD.MODULES.network")
        rows = rn.kr_core()
    finally:
        sys.path.remove(str(speck))
    return {r.rid: r for r in rows}


def test_values_copied_from_rebuild_have_not_drifted(rebuild_rows):
    """The basic network carries its own measured values (no rebuild=
    tag). The one row still copied from REBUILD is the optional KSTEP."""
    tagged = [r for r in network.build(("all",)) if r.rebuild is not None]
    assert [r.rid for r in tagged] == ["KSTEP"]
    for r in tagged:
        rb = rebuild_rows[r.rebuild]
        assert rb.rate[0] == "stepwise" and r.rate[1] == rb.rate[1]
        assert len(rb.rate) == 2          # REBUILD default A = 1.56e-13
        assert r.tg_pow == approx(rb.tg_pow, abs=1e-12), r.rid


def test_species_energies_match_rebuild():
    """E_form drift gate. Kr*, Kr+ (and Kr2+ through Kr+) deliberately use
    the Biagi thresholds (9.915 / 13.996 eV) instead of REBUILD's rounded
    9.92 / 14.00 so collision cost and stored energy are one number."""
    speck = _find_speck()
    if speck is None:
        pytest.skip("SPECK/REBUILD not found above LIGHT_Kr")
    sys.path.insert(0, str(speck))
    try:
        rs = importlib.import_module("REBUILD.MODULES.species")
    finally:
        sys.path.remove(str(speck))
    assert e_form("Kr2*") == approx(rs.REGISTRY["Kr2*"].e_form_eV, abs=1e-12)
    for light in ("Kr*", "Kr+", "Kr2+"):
        assert e_form(light) == approx(rs.REGISTRY[light].e_form_eV,
                                       abs=0.06), light


# ---------------------------------------------------------------------------
# the rows and their bands (LIGHT v0_1, v0_2, v0_4, v0_5)
# ---------------------------------------------------------------------------
def test_charge_and_kr2_rows():
    """Kr+ + 2Kr at the mean of three measurements on Perez-Rios & Greene's
    T^-3/4; Kr2+ + e on Shiu & Biondi onto the 5p rung; the Kr2* route on
    Ueda's three-body rate with Wieme & Lenaerts' Tg form and the 270-ns
    lifetime; Kr2* + Kr2* on a 0-default band."""
    b = {r.rid: r for r in network.build()}
    assert b["K15"].rate == ("const", 2.4e-31) and b["K15"].tg_pow == -0.75
    assert b["R11p"].rate[0] == "te_pow" and b["R11p"].rate[1] == 1.6e-6
    assert b["R11p"].products == ("Kr**", "Kr")
    assert b["AK8_Kr"].rate == ("param", "kr2_form")
    assert b["AK8_Kr"].tg_pow == -0.64
    assert b["KR2_RAD"].rate == ("const", 3.7e6)
    assert b["KR2_RAD"].band == "147 nm (Kr2* VUV)"
    assert b["PEN_KR2"].rate == ("param", "kr2_pen")
    assert b["G2"].rate == ("param", "kr_pen")
    for k in ("AK8_Kr", "KR2_RAD", "PEN_KR2", "K15", "R11p", "G2"):
        assert b[k].rebuild is None and b[k].group == "basic", k


def test_bands_hold_the_published_values():
    """Each band keeps a measured rate at its measurement and spans the
    other published values. A rate with no measurement is 0 by default and
    runs to the largest published value (LIGHT 2026-09-29). Only the
    krypton bands are here."""
    bands = network.BANDS
    assert list(bands) == ["kr2_pen", "kr2_form", "kr_pen", "kps_kr",
                           "k_5s5p", "kstep_s", "kstep_p", "ib_scale"]
    ends = {k: (b.default, b.lo, b.hi) for k, b in bands.items()}
    assert ends["kr2_pen"] == (0.0, 0.0, 5.0e-10)   # Kannari's estimate
    assert ends["kr2_form"] == (1.7e-32, 1.7e-32, 5.5e-32)   # Ueda .. Hughes
    assert ends["kr_pen"] == (5.0e-10, 1.5e-10, 5.0e-10)     # Werner .. v0
    rows = network.build()
    read = network.bands_read(rows)
    assert read == set(bands) - set(network.FREE_BANDS)
    for b in bands.values():
        assert b.lo <= b.default <= b.hi and b.cite and b.what, b.name
        assert all(b.lo <= v <= b.hi for v, _ in b.marks), b.name
        if b.unit not in ("exponent", "factor", "fraction"):
            assert b.lo >= 0.0, b.name
            if b.default == 0.0:         # no measurement: 0 up to the top
                assert b.lo == 0.0 and b.hi > 0.0, b.name


def test_ladder_rows_and_their_routing():
    """The higher excitations make Kr** (5p, 11.30 eV, the Biagi group's
    own threshold); Kr** relaxes, radiates its 1.385 eV untrapped, is
    pumped from Kr* and returned by detailed balance; the stepwise rows
    take the ledger's gaps from the electron gas."""
    by = {r.rid: r for r in network.build()}
    route = network.energy_routing
    for rid in ("C1a", "C1b", "C1c"):
        assert by[rid].products == ("e", "Kr**")
    assert route(by["C1a"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["C1b"])["heat"] == approx(11.998 - 11.30)
    assert route(by["C1c"])["heat"] == approx(12.75 - 11.30)
    gap = e_form("Kr**") - e_form("Kr*")
    assert gap == approx(1.385)
    assert route(by["KPS"])["heat"] == approx(gap)
    prad = route(by["PRAD"])
    assert prad["rad"] == approx(gap) and prad["heat"] == approx(0.0,
                                                                  abs=1e-12)
    assert route(by["K76"])["from_e"] == approx(gap)
    assert route(by["K76r"])["to_e"] == approx(gap)
    # LIGHT v0_5: the pair on Jung et al.'s measured cross sections, the
    # band a scale on both through frac
    assert by["K76"].rate == ("pop_xsec", "jung2006")
    assert by["K76r"].rate == ("pop_xsec", "jung2006", "reverse", 12 / 36,
                               "Kr**", "Kr*")
    assert by["K76"].frac == by["K76r"].frac == (("k_5s5p", "x"),)
    # Params.k76 = 'kannari' gives LIGHT v0_4's fit rows, the band scaling
    byk = {r.rid: r for r in network.build(k76="kannari")}
    assert byk["K76"].rate == ("te_fit", 9.0e-7, 0.50, 1.51)
    assert byk["K76r"].rate[4:] == ("reverse", 12 / 36, "Kr**", "Kr*")
    assert byk["K76"].frac == (("k_5s5p", "x"),)
    assert [r.rid for r in network.build()] == \
        [r.rid for r in network.build(k76="kannari")]
    with pytest.raises(network.NetworkError):
        network.build(k76="jacob")
    assert route(by["KSTEP_S"])["from_e"] == approx(13.996 - 9.915)
    assert route(by["KSTEP_P"])["from_e"] == approx(13.996 - 11.30)
    # Kr2+ + e makes the 5p rung (Shiu & Biondi saw the 5p lines)
    assert route(by["R11p"])["heat"] == approx(13.996 - 1.15 - 11.30)
    # the unmeasured ladder bands are 0, the measured ones carry their
    # measurements, and the preset sets them together
    B = network.BANDS
    for name in ("kstep_s", "kstep_p"):
        assert B[name].default == 0.0, name
    assert B["k_5s5p"].default == 1.0 and (B["k_5s5p"].lo,
                                           B["k_5s5p"].hi) == (0.5, 1.4)
    assert B["kps_kr"].default == 1.7e-10        # Gazeli 2021, one 5p' level
    assert network.KR_5P_TAU_S == approx(29.0e-9, rel=5e-3)
    assert set(network.LADDER_KANNARI) == {"kstep_s", "kstep_p", "kps_kr"}
    for name, v in network.LADDER_KANNARI.items():
        assert B[name].lo <= v <= B[name].hi, name


def test_ion_mobilities_are_ellis():
    """The walls' reduced mobilities (Ellis et al. 1976), krypton's ions
    alone."""
    assert network.ION_K0 == {"Kr+": 0.84, "Kr2+": 1.22}
