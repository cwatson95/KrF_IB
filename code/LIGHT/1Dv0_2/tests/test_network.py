"""The network as data: the notes' list, balance, routing, drift vs SPECK."""

import importlib
import os
import sys
from pathlib import Path

import pytest

network = importlib.import_module("1Dv0_2.network")
_mod = importlib.import_module("1Dv0_2.model")
Params = _mod.Params
_mod = importlib.import_module("1Dv0_2.species")
e_form = _mod.e_form


def approx(expected, rel=1e-6, abs=0.0):
    """pytest.approx without its default absolute tolerance of 1e-12, which
    lets any two rate constants below 1e-12 compare equal (a three-body
    rate of 1.9e-30 'equals' 6.0e-31). Tolerances are relative unless a
    test names an absolute one."""
    return pytest.approx(expected, rel=rel, abs=abs)


def test_v0_is_the_notes_list_plus_three_closures():
    rows = network.build(net="v0")
    assert [r.rid for r in rows] == ["G1", "G2", "G3", "G4", "G5", "G6",
                                     "C1a", "C1b", "C1c", "C2", "C3"]
    assert {r.group for r in rows[:6]} == {"gary"}
    assert {r.group for r in rows[6:]} == {"closure"}


def test_every_row_balances_and_is_cited():
    rows = network.build(("all",), "v0")     # build() runs check_row
    assert len(rows) == 42
    assert all(r.cite for r in rows)
    basic = network.build(("all",))
    assert len(basic) == 46 and all(r.cite for r in basic)


def test_unknown_group_refused():
    with pytest.raises(network.NetworkError):
        network.build(("nonsense",))


def test_energy_routing_examples():
    by = {r.rid: r for r in network.build(("all",), "v0")}
    route = network.energy_routing
    assert route(by["G2"])["to_e"] == approx(2 * 9.915 - 13.996)
    assert route(by["G3"])["heat"] == approx(9.915 - e_form("KrF(B)") - 0.80,
                                             rel=1e-12)
    g4 = route(by["G4"])
    # the 248-nm photon leaves 0.056 eV of X-state repulsion (Tellinghuisen)
    assert (g4["rad"], g4["heat"]) == (5.00, approx(0.056, rel=1e-9))
    assert route(by["C1a"])["heat"] == approx(11.30 - 9.915)
    assert route(by["C2"])["heat"] == approx(0.0, abs=1e-12)
    att = route(by["ATT"])
    assert att["capture"] == 1
    assert att["heat"] == approx(-(e_form("F-") + e_form("F")))
    assert route(by["SUPER"])["to_e"] == approx(9.915)
    assert route(by["CB_Kr"])["heat"] < 0          # uphill, drawn from gas


def test_bc_group_copies_every_KrF_B_consumer_to_C():
    rids = {r.rid for r in network.build(("bc", "quench"), "v0")}
    assert {"G5_C", "K33_C", "K32_C", "K67_C"} <= rids
    assert "G4_C" not in rids                    # C radiates via URAD_C


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
        pytest.skip("SPECK/REBUILD not found above LIGHT (standalone LIGHT); "
                    "set LIGHT_SPECK_DIR to run the drift gate")
    sys.path.insert(0, str(speck))
    try:
        rn = importlib.import_module("REBUILD.MODULES.network")
        rows = (rn.kr_core() + rn.ne_kr_couplings() + rn.ar_kr_manifold()
                + rn.ar_kr_ion_manifold() + rn.upper_level_radiation("krf")
                + rn.krf_bc_mixing())
    finally:
        sys.path.remove(str(speck))
    return {r.rid: r for r in rows}


def test_values_copied_from_rebuild_have_not_drifted(rebuild_rows):
    checked = set()
    rows = network.build(("all",), "v0") + network.build(("all",))
    for r in rows:      # KrF(C) copies re-check their row
        if r.rebuild is None:
            continue
        rb = rebuild_rows[r.rebuild]
        kind = r.rate[0]
        order = len(r.reactants)                  # electrons included
        si_to_cgs = {1: 1.0, 2: 1e6, 3: 1e12}[order]
        if kind == "const":
            assert rb.rate[0] == "const"
            assert r.rate[1] == approx(rb.rate[1] * si_to_cgs,
                                              rel=1e-12), r.rid
        elif kind == "inv_tau":
            assert 1.0 / Params().tau_sp_s == approx(rb.rate[1],
                                                            rel=1e-12)
        elif kind == "stepwise":
            assert rb.rate[0] == "stepwise" and r.rate[1] == rb.rate[1]
            assert len(rb.rate) == 2          # REBUILD default A = 1.56e-13
        elif kind == "dr_sqrt":
            assert r.rate[1] == approx(rb.rate[1] * 1e6, rel=1e-12)
        elif kind == "flannery":
            assert tuple(r.rate[1:]) == tuple(rb.rate[1:])
        else:
            raise AssertionError(f"{r.rid}: no drift rule for {kind}")
        assert r.tg_pow == approx(rb.tg_pow, abs=1e-12), r.rid
        if r.is_rad:
            assert r.hv_eV == approx(rb.hv_eV), r.rid
        checked.add(r.rebuild)
    assert len(checked) == 24


def test_kr2_route_is_basic_in_v0_2_and_rebuilds_in_v0():
    """v0_2 puts the Kr2* route in the basic network (2026-09-29). Ueda's
    three-body rate on Wieme & Lenaerts' Tg form, the 270 ns lifetime and
    Quigley & Hughes' Kr2* + F2 split evenly, each through its band. Kr2* +
    F and Kr2* + Kr2* have no measurement and read 0-default bands. v0
    keeps REBUILD's copies in its optional group, which the drift gate
    checks."""
    b = {r.rid: r for r in network.build()}
    v = {r.rid: r for r in network.build(("kr2",), "v0")}
    assert b["AK8_Kr"].rate == ("param", "kr2_form")
    assert b["AK8_Kr"].tg_pow == -0.64
    assert b["KR2_RAD"].rate == ("const", 3.7e6)
    assert b["NK7"].rate == b["K2F_BR"].rate == ("param", "kr2_f2", 0.5)
    assert b["K2F_F"].rate == ("param", "kr2_f")
    assert b["PEN_KR2"].rate == ("param", "kr2_pen")
    for k in ("AK8_Kr", "KR2_RAD", "NK7", "K2F_BR", "K2F_F", "PEN_KR2"):
        assert b[k].rebuild is None and v[k].rebuild == k, k
        assert b[k].group == "basic" and v[k].group == "kr2", k
        assert b[k].reactants == v[k].reactants, k
        assert b[k].products == v[k].products, k
    assert v["KR2_RAD"].rate == ("const", 8.0e7)
    assert v["K2F_F"].rate == ("const", 3.0e-10)
    assert v["PEN_KR2"].rate == ("const", 5.0e-10)
    with pytest.raises(network.NetworkError, match="Kr2"):
        network.build(("kr2",))


def test_bands_hold_the_published_values():
    """Each band keeps a measured rate at its measurement and spans the
    other published values. A rate with no measurement is 0 by default and
    runs to the largest published value (2026-09-29)."""
    bands = network.BANDS
    assert list(bands) == ["kr2_f", "kr2_pen", "kr2f_f", "kr2_form",
                           "kr2_f2", "kr_pen", "kr2f_f2", "trimer_tg_pow"]
    ends = {k: (b.default, b.lo, b.hi) for k, b in bands.items()}
    assert ends["kr2_f"] == (0.0, 0.0, 1.0e-9)      # SRI's 3-10e-10, Ar2*
    assert ends["kr2_pen"] == (0.0, 0.0, 5.0e-10)   # Kannari's estimate
    assert ends["kr2f_f"] == (0.0, 0.0, 1.8e-10)    # the Ar2F* analog
    assert ends["kr2_form"] == (1.7e-32, 1.7e-32, 5.5e-32)   # Ueda .. Hughes
    assert ends["kr2_f2"] == (2.1e-10, 2.1e-10, 5.0e-10)     # Q&H .. C&P
    assert ends["kr_pen"] == (5.0e-10, 1.5e-10, 5.0e-10)     # Werner .. v0
    assert ends["kr2f_f2"] == (1.2e-10, 1.2e-10, 4.3e-10)    # Xu .. Q&H
    rows = network.build()
    by = {r.rid: r for r in rows}
    assert ends["trimer_tg_pow"] == (by["G5"].tg_pow, by["G5"].tg_pow, -0.5)
    assert by["G5C"].tg_pow == by["G5"].tg_pow
    read = {r.rate[1] for r in rows if r.rate[0] == "param"}
    assert read == set(bands) - {"trimer_tg_pow"}
    for b in bands.values():
        assert b.lo <= b.default <= b.hi and b.cite and b.what, b.name
        assert all(b.lo <= v <= b.hi for v, _ in b.marks), b.name
        if b.unit != "exponent":
            assert b.lo >= 0.0, b.name
            if b.default == 0.0:         # no measurement: 0 up to the top
                assert b.lo == 0.0 and b.hi > 0.0, b.name
    assert (by["G2"].rate, by["Q2F"].rate, by["Q2FF"].rate) == (
        ("param", "kr_pen"), ("param", "kr2f_f2"), ("param", "kr2f_f"))
    assert by["Q2FF"].reactants == ("Kr2F*", "F")
    assert network.energy_routing(by["Q2FF"])["heat"] == approx(
        e_form("Kr2F*") - e_form("F"), rel=1e-12)
    # v0 keeps the notes' G2 and REBUILD's quench values, untouched
    v0 = {r.rid: r for r in network.build(("quench",), "v0")}
    assert v0["G2"].rate == ("const", 5.0e-10) and v0["G2"].rebuild == "PEN_KR"
    assert "Q2FF" not in v0


def test_species_energies_match_rebuild():
    """E_form drift gate. Kr*, Kr+ (and Kr2+ through Kr+) deliberately use
    the Biagi thresholds (9.915 / 13.996 eV) instead of REBUILD's rounded
    9.92 / 14.00 so collision cost and stored energy are one number."""
    speck = _find_speck()
    if speck is None:
        pytest.skip("SPECK/REBUILD not found above LIGHT")
    sys.path.insert(0, str(speck))
    try:
        rs = importlib.import_module("REBUILD.MODULES.species")
    finally:
        sys.path.remove(str(speck))
    names = {"F": "F", "F-": "F-", "Kr2*": "Kr2*"}
    # v0_1 takes KrF(B), KrF(C) and Kr2F* from measurement and theory
    # instead of REBUILD's estimated wells: D_e(B) = 5.54 eV below the ion
    # pair (Tellinghuisen 1976), C 120 cm^-1 below B (Yu & Setser 1990),
    # Kr2F* 0.58 eV below KrF* + Kr (Huestis et al. 1984)
    assert e_form("KrF(B)") == approx(0.80 + 13.996 - 3.40 - 5.54, rel=1e-12)
    assert e_form("KrF(B)") - e_form("KrF(C)") == approx(
        120.0 * 1.23984193e-4, rel=1e-12)
    assert e_form("KrF(B)") - e_form("Kr2F*") == approx(0.58, rel=1e-12)
    for light, rb in names.items():
        assert e_form(light) == approx(rs.REGISTRY[rb].e_form_eV,
                                              abs=1e-12), light
    for light in ("Kr*", "Kr+", "Kr2+"):
        assert e_form(light) == approx(rs.REGISTRY[light].e_form_eV,
                                              abs=0.06), light


def test_basic_is_the_measured_krf_set():
    rids = [r.rid for r in network.build()]
    assert rids == ["G1", "G2", "G3", "G4", "G5", "G6", "C1a", "C1b", "C1c",
                    "C2", "C3", "G3C", "G3X", "G4C", "G4A", "BC", "CB", "G5C",
                    "QF2", "QF2C", "QKR", "QKRC", "QE", "QEC", "Q2F",
                    "Q2FF", "K15", "R11p", "ATT", "DIS_a3", "DIS_A1", "VIB1",
                    "VIB2", "VIB3", "VIB4", "K26", "K108", "FF", "AK8_Kr",
                    "KR2_RAD", "NK7", "K2F_BR", "K2F_F", "PEN_KR2"]
    by = {r.rid: r for r in network.build()}
    # Gadomski 1992 eq. 5 structure with ML's joint-fit constants (adopted
    # 2026-09-28) and K_eq by detailed balance over the ledger's B-C gap
    assert by["BC"].rate == ("const", network.KRF_K_BC) == ("const", 1.05e-11)
    assert by["CB"].rate == ("balance", 1.05e-11, 1.10, "KrF(B)", "KrF(C)")
    assert by["G4C"].rate[1] == approx(1 / 50e-9)
    # B and C each carry their own quench constants, and the basic network's
    # G5 is its own (v0 keeps REBUILD's 6.0e-31)
    assert by["G5"].rate[1] == 1.90e-30 and by["G5C"].rate[1] == 3.79e-31
    assert by["QKR"].rate[1] < by["QKRC"].rate[1]
    assert by["QE"].rate == by["QEC"].rate == ("const", 6.6e-8)
    v0 = {r.rid: r for r in network.build(net="v0")}
    assert v0["G5"].rate == ("const", 6.0e-31) and v0["G5"].rebuild == "K34"
    for grp in ("quench", "bc", "ions", "f2e", "kr2"):
        with pytest.raises(network.NetworkError):
            network.build((grp,))
    # v0_1's charge chemistry and fluorine recombination, with their sources
    assert by["K15"].rate == ("const", 2.4e-31) and by["K15"].tg_pow == -0.75
    assert by["R11p"].rate[0] == "te_pow" and by["R11p"].rate[1] == 1.6e-6
    assert by["K26"].rate == ("flannery", 2.4, 3.417e-6, -3.0)
    assert by["K108"].rate == ("flannery", 2.65, 3.096e-6, -3.0)
    assert by["FF"].rate == ("const", 2.8e-34)
    assert by["FF"].reactants == ("F", "F", "Kr")
    assert by["G4C"].hv_eV == 4.51 and by["G4A"].hv_eV == 4.56
    assert by["C3"].hv_eV == 2.95
