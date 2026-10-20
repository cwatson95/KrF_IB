"""The pure-Kr network as data: its rows, balance, routing, bands, and the
drift gates against SPECK/REBUILD (LIGHT v0_5's tests on krypton)."""

import importlib
import os
import sys
from pathlib import Path

import pytest

from v0_4 import network
from v0_4.species import SPECIES, e_form


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


# the rows that are not a ground-state level, in order: v0_4's resolved 5s
# levels (the Penning pairs, the excimer's formation from each level and
# their reverses, the mixing, the resonance lines, the 5p rung's rows per
# level) and the ion's reverse
LV = ("1s5", "1s4", "1s3", "1s2")
PAIRS = [f"G2_{a}_{b}" for i, a in enumerate(LV) for b in LV[i:]]
OTHER = (PAIRS + ["G6", "C2", "K15", "K15r", "R11p",
                  "AK8_1s5", "AK8_1s4", "AK8f_1s4", "AK8_1s3",
                  "AK8r_1s5", "AK8r_1s4", "AK8fr_1s4", "AK8r_1s3",
                  "KR2_RAD", "KR2F_RAD", "QKR2",
                  "MIX_1s4", "MIX_1s5r", "MIX_1s2", "MIX_1s4r2",
                  "RES_1s4", "RES_1s2", "PEN_KR2"]
         + [f"KPS_{lv}" for lv in LV] + ["KPS_low"]
         + [f"PRAD_{lv}" for lv in LV]
         + [x for lv in LV for x in (f"K76_{lv}", f"K76r_{lv}")]
         + [f"KSTEP_S_{lv}" for lv in LV] + ["KSTEP_P"])
GROUPS = ("G1", "C1a", "C1b", "C1c", "HM")


def test_the_network_is_krypton_and_electrons_alone():
    """126 rows: every one of Biagi 8.97's 51 levels once at its own
    threshold in its group, twenty Hornbeck-Molnar twins, and the other
    rows in order (v0_4: the 5s levels resolved); every species krypton or
    the electron; no optional group."""
    rows = network.build()
    ids = [r.rid for r in rows]
    assert len(rows) == 126 and len(set(ids)) == 126
    grouped = [r.rid for g in GROUPS for r in network.group_rows(rows, g)]
    assert [i for i in ids if i not in grouped] == OTHER
    # the order: the 5s levels, the Penning pairs, G6, the higher levels,
    # HM, C2, the rest
    assert ids[:4] == [f"G1_1s{i}" for i in (5, 4, 3, 2)]
    assert ids[4:14] == PAIRS and ids[14] == "G6" and ids[82] == "C2"
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
    assert set(SPECIES) == {"Kr", "Kr*(1s5)", "Kr*(1s4)", "Kr*(1s3)",
                            "Kr*(1s2)", "Kr**", "Kr+", "Kr2*", "Kr2*(0u+)",
                            "Kr2+"}
    for r in rows + network.optional_rows():
        for s in r.reactants + r.products:
            assert s == "e" or s.startswith("Kr"), (r.rid, s)
    every = network.build(("all",))
    assert [r.rid for r in every] == ids
    assert network.GROUPS == network.BASIC_GROUPS == ()
    assert network.optional_rows() == []


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
    assert len(rows) == 126 and all(r.cite for r in rows)


def test_unknown_or_fluorine_group_refused():
    for g in ("nonsense", "f2e", "kr2", "ions", "quench", "bc", "electron"):
        with pytest.raises(network.NetworkError):
            network.build((g,))


def test_energy_routing_examples():
    by = {r.rid: r for r in network.build(("all",))}
    route = network.energy_routing
    e_kr2 = 9.9152 - 0.709                       # v0_4: the measured well
    assert route(by["G2_1s5_1s5"])["to_e"] == approx(2 * 9.9152 - 13.9996)
    assert route(by["G2_1s4_1s2"])["to_e"] == approx(10.0324 + 10.6436
                                                     - 13.9996)
    assert route(by["C1a_2p10"])["heat"] == approx(0.0, abs=1e-12)
    # v0_4: each 5s level stores its own energy, no heat at excitation
    for lv in LV:
        assert route(by[f"G1_{lv}"])["heat"] == approx(0.0, abs=1e-12)
    assert route(by["C2"])["heat"] == approx(0.0, abs=1e-12)
    # the 147-nm photon leaves the ground state's repulsion as heat
    kr2 = route(by["KR2_RAD"])
    assert (kr2["rad"], kr2["heat"]) == (8.43, approx(e_kr2 - 8.43,
                                                      rel=1e-9))
    fast = route(by["KR2F_RAD"])
    assert fast["rad"] == 8.43 and fast["heat"] == approx(10.0324 - 0.684
                                                           - 8.43, rel=1e-9)
    # three-body association heats the gas by the excimer's well, the
    # reverse draws it back from the gas
    assert route(by["AK8_1s5"])["heat"] == approx(0.709, rel=1e-9)
    assert route(by["AK8r_1s5"])["heat"] == approx(-0.709, rel=1e-9)
    assert route(by["AK8_1s4"])["heat"] == approx(10.0324 - e_kr2, rel=1e-9)
    assert route(by["K15r"])["heat"] == approx(-1.15, rel=1e-9)
    # Kr2* + Kr2* Penning gives its excess to the new electron
    assert route(by["PEN_KR2"])["to_e"] == approx(2 * e_kr2 - (13.9996
                                                               - 1.15),
                                                  rel=1e-12)
    # the escaping resonance photon carries the level's whole energy
    assert route(by["RES_1s4"])["rad"] == approx(10.0324, rel=1e-12)
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
    """The network carries its own measured values: v0_4 has no row
    copied from REBUILD (v0_3's optional KSTEP, the last, read the lumped
    Kr*)."""
    tagged = [r for r in network.build(("all",)) if r.rebuild is not None]
    assert tagged == []


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
    # v0_4: Kr2* at its measured well, 0.006 eV above REBUILD's estimate
    assert e_form("Kr2*") == approx(rs.REGISTRY["Kr2*"].e_form_eV, abs=0.01)
    assert e_form("Kr*(1s5)") == approx(rs.REGISTRY["Kr*"].e_form_eV,
                                        abs=0.06)
    for light in ("Kr+", "Kr2+"):
        assert e_form(light) == approx(rs.REGISTRY[light].e_form_eV,
                                       abs=0.06), light


# ---------------------------------------------------------------------------
# the rows and their bands (LIGHT v0_1, v0_2, v0_4, v0_5)
# ---------------------------------------------------------------------------
def test_charge_and_kr2_rows():
    """Kr+ + 2Kr at the mean of three measurements on Perez-Rios & Greene's
    T^-3/4, and its reverse by detailed balance; Kr2+ + e on Shiu & Biondi
    onto the 5p rung, its gas-temperature band; the Kr2* route from each 5s
    level on its own band with Wieme & Lenaerts' Tg form, each with its
    reverse; the 270-ns and 3.4-ns lifetimes; Kr2* + Kr2* and Kr2* + Kr
    on 0-default bands; the Penning pairs on kr_pen, 2k across levels."""
    b = {r.rid: r for r in network.build()}
    assert b["K15"].rate == ("const", 2.4e-31) and b["K15"].tg_pow == -0.75
    assert b["K15r"].rate == ("eq_reverse", ("const", 2.4e-31), -0.75,
                              "Kr2+", 4)
    assert b["R11p"].rate[0] == "te_pow" and b["R11p"].rate[1] == 1.6e-6
    assert b["R11p"].products == ("Kr**", "Kr")
    assert b["R11p"].tg_fac[0] == "dr_tg"
    for rid, band in (("AK8_1s5", "kr2_form"), ("AK8_1s4", "kr2_form_1s4"),
                      ("AK8f_1s4", "kr2_form_1s4"),
                      ("AK8_1s3", "kr2_form_1s3")):
        assert b[rid].rate == ("param", band) and b[rid].tg_pow == -0.64
        rev = b[rid.replace("AK8", "AK8r").replace("AK8rf", "AK8fr")]
        assert rev.rate[:3] == ("eq_reverse", ("param", band), -0.64)
        assert rev.tg_pow == 0.0 and rev.frac == b[rid].frac
    assert b["AK8_1s4"].frac == (("kr2_fast", "1-x"),)
    assert b["AK8f_1s4"].frac == (("kr2_fast", "x"),)
    assert b["AK8f_1s4"].products == ("Kr2*(0u+)", "Kr")
    assert b["KR2_RAD"].rate == ("const", 3.7e6)
    assert b["KR2_RAD"].band == b["KR2F_RAD"].band == "147 nm (Kr2* VUV)"
    assert b["KR2F_RAD"].rate == ("const", approx(1.0 / 3.4e-9, rel=1e-12))
    assert b["PEN_KR2"].rate == ("param", "kr2_pen")
    assert b["QKR2"].rate == ("param", "kr2_quench")
    for pair in PAIRS:
        a, c = pair.split("_")[1:]
        assert b[pair].rate == ("param", "kr_pen", 1.0 if a == c else 2.0)
    for k in ("AK8_1s5", "KR2_RAD", "PEN_KR2", "K15", "R11p", "G2_1s5_1s5"):
        assert b[k].rebuild is None and b[k].group == "basic", k


def test_bands_hold_the_published_values():
    """Each band keeps a measured rate at its measurement and spans the
    other published values. A rate with no measurement is 0 by default and
    runs to the largest published value (LIGHT 2026-09-29). Only the
    krypton bands are here."""
    bands = network.BANDS
    assert list(bands) == ["kr2_pen", "kr2_form", "kr2_form_1s4",
                           "kr2_form_1s3", "kr2_fast", "kr2_quench",
                           "kr2_d0", "dr_tg", "res_escape", "kps_5s",
                           "kr_pen", "kps_kr", "k_5s5p", "k_5s5p_j1",
                           "kstep_s", "kstep_p", "ib_scale", "hm_kr"]
    ends = {k: (b.default, b.lo, b.hi) for k, b in bands.items()}
    assert ends["kr2_pen"] == (0.0, 0.0, 5.0e-10)   # Kannari's estimate
    # v0_4: Marchal's state-selective 1s5 value, Ueda .. Smith & Turner
    assert ends["kr2_form"] == (3.81e-32, 1.7e-32, 4.2e-32)
    # Marchal's 3P1, down to Audouard's ratio 3/1
    assert ends["kr2_form_1s4"] == (9.24e-32, 1.27e-32, 9.24e-32)
    assert ends["kr2_form_1s3"] == (3.81e-32, 0.0, 9.24e-32)
    assert ends["kr2_fast"] == (0.0, 0.0, 0.17)
    assert ends["kr2_quench"] == (0.0, 0.0, 8.6e-15)
    assert ends["kr2_d0"] == (0.709, 0.62, 0.709)
    assert ends["dr_tg"] == (0.0, 0.0, 1.0)
    assert ends["res_escape"] == (0.0, 0.0, 1.7e5)
    assert ends["kps_5s"] == (1.0, 0.0, 1.0)
    assert ends["k_5s5p"] == (1.0, 0.6, 1.4)
    assert ends["k_5s5p_j1"] == (1.0, 0.0, 1.0)
    assert ends["kr_pen"] == (5.0e-10, 1.5e-10, 5.0e-10)     # Werner .. v0
    # v0_2: no measurement, 0; every such excitation ionizing, the top
    assert ends["hm_kr"] == (0.0, 0.0, 1.0)
    rows = network.build()
    read = network.bands_read(rows)
    assert read == set(bands) - set(network.FREE_BANDS)
    for b in bands.values():
        assert b.lo <= b.default <= b.hi and b.cite and b.what, b.name
        assert all(b.lo <= v <= b.hi for v, _ in b.marks), b.name
        if b.unit not in ("exponent", "factor", "fraction", "eV"):
            assert b.lo >= 0.0, b.name
            if b.default == 0.0:         # no measurement: 0 up to the top
                assert b.lo == 0.0 and b.hi > 0.0, b.name


def test_ladder_rows_and_their_routing():
    """The higher excitations make Kr** (5p, 11.3035 eV, the 2p10 level),
    each paying its own threshold with the excess to heat; v0_4: each 5s
    level is its own species; Kr** relaxes to the levels (KPS_*, the
    statistical share or all to 1s5), radiates to them on NIST's branching
    (PRAD_*), is pumped from each and returned by detailed balance (K76_*,
    K76r_*), the J = 1 levels on their proxies' band as well; the stepwise
    rows take each level's own gap from the electron gas."""
    rows = network.build()
    by = {r.rid: r for r in rows}
    route = network.energy_routing
    for g in ("C1a", "C1b", "C1c"):
        for r in network.group_rows(rows, g):
            assert r.products == ("e", "Kr**"), r.rid
            assert route(r)["heat"] == approx(r.e_cost_eV - 11.3035,
                                              abs=1e-12), r.rid
    for r in network.group_rows(rows, "G1"):
        lv = r.rid[3:]
        assert r.products == ("e", f"Kr*({lv})")
        assert route(r)["heat"] == approx(0.0, abs=1e-12)
    assert route(by["C1b_3d6"])["heat"] == approx(11.9981 - 11.3035)
    assert route(by["C1c_high"])["heat"] == approx(13.6 - 11.3035)
    G = network.G_LEVEL_5S
    for lv in LV:
        gap = e_form("Kr**") - e_form(f"Kr*({lv})")
        kps = by[f"KPS_{lv}"]
        assert kps.rate == ("param", "kps_kr", G[lv] / 12)
        assert kps.frac == (("kps_5s", "x"),)
        assert route(kps)["heat"] == approx(gap)
        prad = route(by[f"PRAD_{lv}"])
        assert prad["rad"] == approx(gap) and prad["heat"] == approx(
            0.0, abs=1e-12)
        assert by[f"PRAD_{lv}"].rate[1] == approx(
            network.KR_5P_BRANCH[lv] / network.KR_5P_TAU_S, rel=1e-12)
        assert route(by[f"K76_{lv}"])["from_e"] == approx(gap)
        assert route(by[f"K76r_{lv}"])["to_e"] == approx(gap)
        assert by[f"K76_{lv}"].rate == ("pop_xsec", "jung2006", lv)
        assert by[f"K76r_{lv}"].rate == ("pop_xsec", "jung2006", lv,
                                         "reverse", G[lv] / 36, "Kr**",
                                         f"Kr*({lv})")
        frac = (("k_5s5p", "x"),) + ((("k_5s5p_j1", "x"),)
                                     if lv in ("1s4", "1s2") else ())
        assert by[f"K76_{lv}"].frac == by[f"K76r_{lv}"].frac == frac
        assert route(by[f"KSTEP_S_{lv}"])["from_e"] == approx(
            13.9996 - e_form(f"Kr*({lv})"))
        assert by[f"KSTEP_S_{lv}"].rate[3] == approx(
            4.04 - (e_form(f"Kr*({lv})") - 9.9152), rel=1e-12)
    assert by["KPS_low"].products == ("Kr*(1s5)", "Kr")
    assert by["KPS_low"].frac == (("kps_5s", "1-x"),)
    # Params.k76 = 'kannari' gives LIGHT v0_4's fit out of every level
    byk = {r.rid: r for r in network.build(k76="kannari")}
    for lv in LV:
        assert byk[f"K76_{lv}"].rate == ("te_fit", 9.0e-7, 0.50, 1.51)
        assert byk[f"K76r_{lv}"].rate[4:] == ("reverse", G[lv] / 36, "Kr**",
                                              f"Kr*({lv})")
        assert byk[f"K76_{lv}"].frac == (("k_5s5p", "x"),)
    assert [r.rid for r in network.build()] == \
        [r.rid for r in network.build(k76="kannari")]
    with pytest.raises(network.NetworkError):
        network.build(k76="jacob")
    assert route(by["KSTEP_P"])["from_e"] == approx(13.9996 - 11.3035)
    # Kr2+ + e makes the 5p rung (Shiu & Biondi saw the 5p lines)
    assert route(by["R11p"])["heat"] == approx(13.9996 - 1.15 - 11.3035)
    # the unmeasured ladder bands are 0, the measured ones carry their
    # measurements, and the preset sets them together
    B = network.BANDS
    for name in ("kstep_s", "kstep_p"):
        assert B[name].default == 0.0, name
    assert B["kps_kr"].default == 1.7e-10        # Gazeli 2021, one 5p' level
    assert network.KR_5P_TAU_S == approx(29.0e-9, rel=5e-3)
    assert set(network.LADDER_KANNARI) == {"kstep_s", "kstep_p", "kps_kr"}
    for name, v in network.LADDER_KANNARI.items():
        assert B[name].lo <= v <= B[name].hi, name


def test_ion_mobilities_are_ellis():
    """The walls' reduced mobilities (Ellis et al. 1976), krypton's ions
    alone."""
    assert network.ION_K0 == {"Kr+": 0.84, "Kr2+": 1.22}
