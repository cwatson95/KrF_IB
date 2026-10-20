"""The network as data: the notes' list, balance, routing, drift vs SPECK."""

import importlib
import os
import sys
from pathlib import Path

import pytest

from v0 import network
from v0.model import Params
from v0.species import e_form


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
    assert len(basic) == 43 and all(r.cite for r in basic)


def test_unknown_group_refused():
    with pytest.raises(network.NetworkError):
        network.build(("nonsense",))


def test_energy_routing_examples():
    by = {r.rid: r for r in network.build(("all",), "v0")}
    route = network.energy_routing
    assert route(by["G2"])["to_e"] == pytest.approx(2 * 9.915 - 13.996)
    assert route(by["G3"])["heat"] == pytest.approx(9.915 - 6.15 - 0.80)
    g4 = route(by["G4"])
    assert (g4["rad"], g4["heat"]) == (5.00, pytest.approx(0.35))
    assert route(by["C1a"])["heat"] == pytest.approx(11.30 - 9.915)
    assert route(by["C2"])["heat"] == pytest.approx(0.0, abs=1e-12)
    att = route(by["ATT"])
    assert att["capture"] == 1
    assert att["heat"] == pytest.approx(-(e_form("F-") + e_form("F")))
    assert route(by["SUPER"])["to_e"] == pytest.approx(9.915)
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
            assert r.rate[1] == pytest.approx(rb.rate[1] * si_to_cgs,
                                              rel=1e-12), r.rid
        elif kind == "inv_tau":
            assert 1.0 / Params().tau_sp_s == pytest.approx(rb.rate[1],
                                                            rel=1e-12)
        elif kind == "stepwise":
            assert rb.rate[0] == "stepwise" and r.rate[1] == rb.rate[1]
            assert len(rb.rate) == 2          # REBUILD default A = 1.56e-13
        elif kind == "dr_sqrt":
            assert r.rate[1] == pytest.approx(rb.rate[1] * 1e6, rel=1e-12)
        elif kind == "flannery":
            assert tuple(r.rate[1:]) == tuple(rb.rate[1:])
        else:
            raise AssertionError(f"{r.rid}: no drift rule for {kind}")
        assert r.tg_pow == pytest.approx(rb.tg_pow, abs=1e-12), r.rid
        if r.is_rad:
            assert r.hv_eV == pytest.approx(rb.hv_eV), r.rid
        checked.add(r.rebuild)
    assert len(checked) == 24


def test_kr2_route_is_measured_in_basic_and_rebuilds_in_v0():
    """The basic network's Kr2* route carries the measured constants (step
    2, 2026-09-28): Ueda's three-body rate on Wieme & Lenaerts' Tg form, the
    270 ns lifetime and Quigley & Hughes' Kr2* + F2 split evenly. v0 keeps
    REBUILD's copies, which the drift gate checks."""
    b = {r.rid: r for r in network.build(("kr2",))}
    v = {r.rid: r for r in network.build(("kr2",), "v0")}
    assert b["AK8_Kr"].rate == ("const", 1.7e-32)
    assert b["AK8_Kr"].tg_pow == -0.64
    assert b["KR2_RAD"].rate == ("const", 3.7e6)
    assert b["NK7"].rate[1] + b["K2F_BR"].rate[1] == pytest.approx(2.1e-10)
    for k in ("AK8_Kr", "KR2_RAD", "NK7", "K2F_BR"):
        assert b[k].rebuild is None and v[k].rebuild == k, k
        assert b[k].group == v[k].group == "kr2"
    assert v["KR2_RAD"].rate == ("const", 8.0e7)
    assert b["K2F_F"] == v["K2F_F"] and b["PEN_KR2"] == v["PEN_KR2"]


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
    names = {"F": "F", "F-": "F-", "KrF(B)": "KrF*(B)", "KrF(C)": "KrF*(C)",
             "Kr2F*": "Kr2F*", "Kr2*": "Kr2*"}
    for light, rb in names.items():
        assert e_form(light) == pytest.approx(rs.REGISTRY[rb].e_form_eV,
                                              abs=1e-12), light
    for light in ("Kr*", "Kr+", "Kr2+"):
        assert e_form(light) == pytest.approx(rs.REGISTRY[light].e_form_eV,
                                              abs=0.06), light


def test_basic_is_the_measured_krf_set():
    rids = [r.rid for r in network.build()]
    assert rids == ["G1", "G2", "G3", "G4", "G5", "G6", "C1a", "C1b", "C1c",
                    "C2", "C3", "G3C", "G3X", "G4C", "BC", "CB", "G5C", "QF2",
                    "QF2C", "QKR", "QKRC", "QE", "QEC", "Q2F"]
    by = {r.rid: r for r in network.build()}
    # Gadomski 1992 eq. 5 structure with ML's joint-fit constants (adopted
    # 2026-09-28) and K_eq by detailed balance over the ledger's B-C gap
    assert by["BC"].rate == ("const", network.KRF_K_BC) == ("const", 1.05e-11)
    assert by["CB"].rate == ("balance", 1.05e-11, 1.10, "KrF(B)", "KrF(C)")
    assert by["G4C"].rate[1] == pytest.approx(1 / 50e-9)
    # B and C each carry their own quench constants, and the basic network's
    # G5 is its own (v0 keeps REBUILD's 6.0e-31)
    assert by["G5"].rate[1] == 1.90e-30 and by["G5C"].rate[1] == 3.79e-31
    assert by["QKR"].rate[1] < by["QKRC"].rate[1]
    assert by["QE"].rate == by["QEC"].rate == ("const", 6.6e-8)
    v0 = {r.rid: r for r in network.build(net="v0")}
    assert v0["G5"].rate == ("const", 6.0e-31) and v0["G5"].rebuild == "K34"
    for grp in ("quench", "bc"):
        with pytest.raises(network.NetworkError):
            network.build((grp,))
