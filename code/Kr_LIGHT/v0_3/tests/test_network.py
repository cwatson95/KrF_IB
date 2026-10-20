"""The pure-Kr network as data: its rows, balance, routing, bands, and the
drift gates against SPECK/REBUILD (LIGHT v0_5's tests on krypton)."""

import importlib
import os
import sys
from pathlib import Path

import pytest

from v0_3 import network
from v0_3.species import SPECIES, e_form


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


# v0_2: the rows that are not a ground-state level, v0_1's in v0_1's order
OTHER = ["G2", "G6", "C2", "K15", "R11p", "AK8_Kr", "KR2_RAD", "PEN_KR2",
         "KPS", "PRAD", "K76", "K76r", "KSTEP_S", "KSTEP_P"]
GROUPS = ("G1", "C1a", "C1b", "C1c", "HM")


def test_the_network_is_krypton_and_electrons_alone():
    """85 rows: every one of Biagi 8.97's 51 levels once at its own
    threshold in its group, twenty Hornbeck-Molnar twins, and v0_1's other
    rows in v0_1's order; every species krypton or the electron; the
    optional group adds KSTEP and SUPER."""
    rows = network.build()
    ids = [r.rid for r in rows]
    assert len(rows) == 85 and len(set(ids)) == 85
    grouped = [r.rid for g in GROUPS for r in network.group_rows(rows, g)]
    assert [i for i in ids if i not in grouped] == OTHER
    # the order: the 5s levels, G2, G6, the higher levels, HM, C2, the rest
    assert ids[:4] == [f"G1_1s{i}" for i in (5, 4, 3, 2)]
    assert ids[4:6] == ["G2", "G6"] and ids[73] == "C2"
    assert [len(network.group_rows(rows, g)) for g in GROUPS] == \
        [4, 10, 10, 27, 20]
    # each level of the file is one excitation row, at its own threshold,
    # the groups in turn and each group in the file's order (C1a holds the
    # 5p' levels at 12.10-12.26 eV, which the file lists after 3d6)
    exc = [r for r in rows if r.rate[0] == "xsec" and r.rate[2] ==
           "excitation"]
    levels = network.levels()
    assert [r.rate[3] for r in exc] == [
        p.product for g in GROUPS[:4] for p in levels
        if network.level_group(p) == g]
    by_product = {p.product: p for p in levels}
    assert len(by_product) == len(exc) == 51
    for r in exc:
        p = by_product[r.rate[3]]
        assert r.e_cost_eV == p.threshold, r.rid
        assert r.rid == f"{network.level_group(p)}_" \
                        f"{network.level_tag(p.product)}"
    assert set(SPECIES) == {"Kr", "Kr*", "Kr**", "Kr+", "Kr2*", "Kr2+"}
    for r in rows + network.optional_rows():
        for s in r.reactants + r.products:
            assert s == "e" or s.startswith("Kr"), (r.rid, s)
    every = network.build(("all",))
    assert [r.rid for r in every] == ids + ["KSTEP", "SUPER"]
    assert network.GROUPS == network.BASIC_GROUPS == ("electron",)


def test_level_groups_and_tags():
    """The groups by the levels' Paschen labels and thresholds, and the
    row ids from Biagi's labels (his '!' marks levels sharing a label)."""
    tag = network.level_tag
    assert tag("Kr(3D1!!)") == "3d1_2" and tag("Kr(3S1!!!!") == "3s1_4"
    assert tag("Kr(4F SUM") == "4f" and tag("Kr(HIGH)") == "high"
    assert tag("Kr(1S5)") == "1s5" and tag("Kr(2P10)") == "2p10"
    g = {network.level_tag(p.product): (network.level_group(p), p.threshold)
         for p in network.levels()}
    assert {k for k, (grp, _) in g.items() if grp == "G1"} == \
        {"1s5", "1s4", "1s3", "1s2"}
    assert {k for k, (grp, _) in g.items() if grp == "C1a"} == \
        {f"2p{i}" for i in range(1, 11)}
    c1b = {k for k, (grp, _) in g.items() if grp == "C1b"}
    assert c1b == {"3d6", "3d5", "3d3", "3d4_1", "3d4", "3d1_2", "3d1_1",
                   "2s5", "3d2", "2s4"}
    assert max(t for grp, t in g.values() if grp == "C1b") == 12.3853
    assert min(t for grp, t in g.values() if grp == "C1c") == 12.7564
    # the Hornbeck-Molnar levels: at or above E(Kr2+) = IP - D0
    assert network.E_ASSOC_EV == approx(13.9996 - 1.15, rel=1e-12)
    hm = sorted(t for k, (grp, t) in g.items()
                if grp == "C1c" and t >= network.E_ASSOC_EV)
    assert len(hm) == 20 and hm[0] == 12.8573 and hm[-1] == 13.6
    below = max(t for grp, t in g.values()
                if grp == "C1c" and t < network.E_ASSOC_EV)
    assert below == 12.8253


def test_hornbeck_molnar_twins():
    """Each level at or above E(Kr2+) is two rows: the cascade at the
    share 1 - hm_kr and the associative ionization at hm_kr, on the same
    cross section and threshold; the ionization's rate is first order in
    Kr (model.py), its electron pays the threshold, the excess over Kr2+
    heats the gas and the new electron is born at rest."""
    rows = network.build()
    by = {r.rid: r for r in rows}
    route = network.energy_routing
    for r in network.group_rows(rows, "HM"):
        twin = by["C1c_" + r.rid[3:]]
        assert r.rate == ("assoc",) + twin.rate[1:]
        assert r.frac == (("hm_kr", "x"),)
        assert twin.frac == (("hm_kr", "1-x"),)
        assert r.e_cost_eV == twin.e_cost_eV >= network.E_ASSOC_EV
        assert (r.reactants, r.products) == (("e", "Kr", "Kr"),
                                             ("e", "e", "Kr2+"))
        assert r.ne_delta == 1
        rt = route(r)
        assert rt["from_e"] == r.e_cost_eV and rt["to_e"] == 0.0
        assert rt["heat"] == approx(r.e_cost_eV - e_form("Kr2+"), rel=1e-12)
        assert rt["heat"] >= 0.0
    # every other level has no band
    for g in ("G1", "C1a", "C1b"):
        assert all(r.frac == () for r in network.group_rows(rows, g)), g
    assert network.BANDS["hm_kr"].default == 0.0


def test_every_row_balances_and_is_cited():
    rows = network.build(("all",))          # build() runs check_row
    assert len(rows) == 87 and all(r.cite for r in rows)


def test_unknown_or_fluorine_group_refused():
    for g in ("nonsense", "f2e", "kr2", "ions", "quench", "bc"):
        with pytest.raises(network.NetworkError):
            network.build((g,))


def test_energy_routing_examples():
    by = {r.rid: r for r in network.build(("all",))}
    route = network.energy_routing
    assert route(by["G2"])["to_e"] == approx(2 * 9.9152 - 13.9996)
    assert route(by["C1a_2p10"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["G1_1s5"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["G1_1s2"])["heat"] == approx(10.6436 - 9.9152)
    assert route(by["C2"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["SUPER"])["to_e"] == approx(9.9152)
    # the 147-nm photon leaves 0.77 eV of ground-state repulsion with the
    # ledger's Kr2* at 9.2 eV
    kr2 = route(by["KR2_RAD"])
    assert (kr2["rad"], kr2["heat"]) == (8.43, approx(9.2 - 8.43, rel=1e-9))
    # three-body association heats the gas by the excimer's well
    assert route(by["AK8_Kr"])["heat"] == approx(9.9152 - 9.2, rel=1e-9)
    # Kr2* + Kr2* Penning gives its excess to the new electron
    assert route(by["PEN_KR2"])["to_e"] == approx(2 * 9.2 - (13.9996 - 1.15),
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
    the Biagi thresholds (Magboltz 8.97: 9.9152 / 13.9996 eV) instead of
    REBUILD's rounded 9.92 / 14.00 so collision cost and stored energy are
    one number."""
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
                           "k_5s5p", "kstep_s", "kstep_p", "ib_scale",
                           "hm_kr"]
    ends = {k: (b.default, b.lo, b.hi) for k, b in bands.items()}
    assert ends["kr2_pen"] == (0.0, 0.0, 5.0e-10)   # Kannari's estimate
    assert ends["kr2_form"] == (1.7e-32, 1.7e-32, 5.5e-32)   # Ueda .. Hughes
    assert ends["kr_pen"] == (5.0e-10, 1.5e-10, 5.0e-10)     # Werner .. v0
    # v0_2: no measurement, 0; every such excitation ionizing, the top
    assert ends["hm_kr"] == (0.0, 0.0, 1.0)
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
    """The higher excitations make Kr** (5p, 11.3035 eV, the 2p10 level),
    each paying its own threshold with the excess to heat; Kr** relaxes,
    radiates its 1.388 eV untrapped, is pumped from Kr* and returned by
    detailed balance; the stepwise rows take the ledger's gaps from the
    electron gas."""
    rows = network.build()
    by = {r.rid: r for r in rows}
    route = network.energy_routing
    for g in ("C1a", "C1b", "C1c"):
        for r in network.group_rows(rows, g):
            assert r.products == ("e", "Kr**"), r.rid
            assert route(r)["heat"] == approx(r.e_cost_eV - 11.3035,
                                              abs=1e-12), r.rid
    for r in network.group_rows(rows, "G1"):
        assert r.products == ("e", "Kr*")
        assert route(r)["heat"] == approx(r.e_cost_eV - 9.9152, abs=1e-12)
    assert route(by["C1b_3d6"])["heat"] == approx(11.9981 - 11.3035)
    assert route(by["C1c_high"])["heat"] == approx(13.6 - 11.3035)
    gap = e_form("Kr**") - e_form("Kr*")
    assert gap == approx(11.3035 - 9.9152, rel=1e-12)
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
    assert route(by["KSTEP_S"])["from_e"] == approx(13.9996 - 9.9152)
    assert route(by["KSTEP_P"])["from_e"] == approx(13.9996 - 11.3035)
    # Kr2+ + e makes the 5p rung (Shiu & Biondi saw the 5p lines)
    assert route(by["R11p"])["heat"] == approx(13.9996 - 1.15 - 11.3035)
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
