"""LIGHT v0_5's 5s -> 5p electron excitation on Jung et al.'s measured
cross sections (jung2006.py), under both closures, with Kannari's fit kept
as Params.k76 = 'kannari', on pure krypton. Every check is a
falsification: the tables as printed, the lumped cross section's sum, its
Maxwellian rate against a direct quadrature, the reverse by detailed
balance, the kinetic solver's population process against the tables, the
fit switch, and the pure-Kr reference point pinned.

LIGHT's version of this file compared through pytest's default approx,
whose absolute tolerance of 1e-12 passes any two cross sections of 1e-16
cm2 (the Known-bad list); here every comparison is relative (abs=0)."""

import math

import numpy as np
import pytest

from v0_4 import boltzmann as KB
from v0_4 import jung2006 as J
from v0_4 import xsec
from v0_4.model import Model, Params, run
from v0_4.report import accounting, row_dict
from v0_4.species import e_form


def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)


_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# the lumped process of v0_1 to v0_3 (kept for comparison) at the 1s5
# energy; v0_4's levels each pay their own gap
GAP = e_form("Kr**") - e_form("Kr*(1s5)")
LEVELS = ("1s5", "1s4", "1s3", "1s2")
G5 = {"1s5": 5, "1s4": 3, "1s3": 1, "1s2": 3}
# the pure-Kr reference point under LIGHT v0_2's closure (two_t): v0_4 on
# the resolved 5s levels and the field tail table (2026-10-15; on v0_3's
# beam table 33.838, 0.44541, 691.16 K); v0_3 gave 15.7759, 0.62625, 798.08 K
# and, with Kannari's fit, 17.175, 0.60692, 789.81 K; v0_2 (the file as
# published) 15.9976, 0.62578, 796.76 K and 17.452, 0.60584, 788.16 K;
# v0_1 on Biagi-v7.1's groups 15.051, 0.6340, 802.38 K and 16.278, 0.6172,
# 795.12 K
V0_4_REF = dict(emission_vuv=32.9227, emission_nir=0.45166, tg_end_K=696.59)
KANNARI_REF = dict(emission_vuv=33.5856, emission_nir=0.44161,
                   tg_end_K=692.64)
FAST = dict(eedf="two_t")


# Jung et al., PRA 73, 022722 (2006), as printed (p.022722-6 and -7),
# 10^-16 cm2: Table II out of 1s5, energies 2-8 eV for 2p2-2p10 and the
# columns that go on to 300 eV; Table III out of 1s3, 2p1-2p4
PAPER_II = {
    2: (0, 0, 0, 0.32, 8.6, 4.7, 4.8, 27, 9.8),
    3: (0.47, 0.5, 0.5, 0.76, 18, 3.8, 7.6, 47, 15),
    4: (0.66, 0.4, 0.3, 0.34, 23, 2.8, 10, 52, 15),
    6: (0.30, 0.2, 0.1, 0.11, 25, 2.6, 9.6, 52, 17),
    8: (0.10, 0.2, 0.05, 0.04, 24, 2.2, 8.6, 51, 17),
}
PAPER_II_LOWER = ("2p2", "2p3", "2p4", "2p5", "2p6", "2p7", "2p8", "2p9",
                  "2p10")
PAPER_II_HIGH = {10: {"2p6": 22, "2p8": 8.2, "2p9": 49},
                 15: {"2p6": 17, "2p9": 41}, 25: {"2p6": 14, "2p9": 37},
                 50: {"2p6": 8.8, "2p9": 23}, 100: {"2p6": 4.8, "2p9": 16},
                 200: {"2p6": 3.1, "2p9": 8.9}, 300: {"2p6": 2.5, "2p9": 6.0}}
PAPER_III = {
    2: (0.90, 7.0, 17, 18), 3: (1.3, 18, 43, 48), 4: (1.1, 14, 55, 65),
    6: (0.63, 5.5, 59, 62), 8: (0.23, 3.6, 60, 63),
}
PAPER_III_HIGH = {10: {"2p3": 56, "2p4": 67}}


def test_tables_are_the_papers():
    """Every tabulated entry of Tables II and III equals the paper's, typed
    again here from the PDF (2026-10-13). LIGHT v0_5 carried 21 for 1s5 ->
    2p6 at 4 eV, where the paper prints 23, and booked Table IV's last row,
    1s5 -> 2p5, as a core-changing 1s3 -> 2p5; LIGHT's test of the 4-eV
    row passed through pytest's default absolute tolerance."""
    for e, row in PAPER_II.items():
        assert J.TABLE_II[e] == dict(zip(PAPER_II_LOWER, map(float, row))), e
    for e, row in PAPER_II_HIGH.items():
        assert J.TABLE_II[e] == {k: float(v) for k, v in row.items()}, e
    assert set(J.TABLE_II) == set(PAPER_II) | set(PAPER_II_HIGH)
    for e, row in PAPER_III.items():
        assert J.TABLE_III[e] == dict(zip(("2p1", "2p2", "2p3", "2p4"),
                                          map(float, row))), e
    for e, row in PAPER_III_HIGH.items():
        assert J.TABLE_III[e] == {k: float(v) for k, v in row.items()}, e
    assert set(J.TABLE_III) == set(PAPER_III) | set(PAPER_III_HIGH)


def test_tables_as_printed():
    """Spot values of Jung et al. 2006 Tables II and III through the
    interpolation, and its shape outside the tables."""
    assert J.sigma_pair("1s5", "2p9", 6.0) == approx(52e-16)
    assert J.sigma_pair("1s5", "2p6", 2.0) == approx(8.6e-16)
    assert J.sigma_pair("1s5", "2p6", 4.0) == approx(23e-16)
    assert J.sigma_pair("1s5", "2p9", 300.0) == approx(6.0e-16)
    assert J.sigma_pair("1s3", "2p4", 4.0) == approx(65e-16)
    assert J.sigma_pair("1s3", "2p1", 8.0) == approx(0.23e-16)
    # 1s5 -> 2p5 (Table II, and Table IV's last row at 6 eV)
    assert J.sigma_pair("1s5", "2p5", 6.0) == approx(0.11e-16)
    # nothing is measured out of 1s3 into the lower tier
    for f in ("2p5", "2p6", "2p7", "2p8", "2p9", "2p10"):
        assert J.sigma_pair("1s3", f, 3.5) == 0.0, f
    # below each pair's threshold the cross section is 0; the lowest
    # threshold is 1s5 -> 2p10 at 1.39 eV
    assert J.sigma_pair("1s5", "2p10", 1.3) == 0.0
    assert J.sigma_pair("1s5", "2p10", 1.5) > 0.0
    assert J.sigma_pair("1s3", "2p4", 1.5) == 0.0          # threshold 1.54
    # nothing measured out of 1s5 into 2p1
    assert J.sigma_pair("1s5", "2p1", 6.0) == 0.0
    # above its last point (10 eV, 8.2) a column follows the 1s5 -> 2p9
    # shape (49 there, 6.0 at 300 eV), then ln E / E. LIGHT's version of
    # this check took the 8-eV point (8.6 x 6.0/51) and passed only through
    # pytest's default absolute tolerance
    s300, s600 = J.sigma_pair("1s5", "2p8", [300.0, 600.0])
    assert s300 == approx(8.2e-16 * 6.0 / 49.0, rel=1e-9)
    assert s600 / s300 == approx(math.log(600.0) / 600.0
                                 / (math.log(300.0) / 300.0), rel=1e-9)


def test_lumped_cross_section_is_the_manifold_mean():
    """At 6 eV the sums of Tables II (1s5) and III (1s3) and the
    statistical mean with the J = 1 proxies; without the proxies the
    metastables alone carry 6 of 12."""
    assert J.level_total("1s5", 6.0) == approx(106.91e-16, rel=1e-6)
    assert J.level_total("1s3", 6.0) == approx(127.13e-16, rel=1e-6)
    want = (8 * 106.91e-16 + 4 * 127.13e-16) / 12
    # the proxies see the same final levels at their own thresholds, so
    # at 6 eV the shifted evaluation differs slightly from the plain sum
    assert J.lumped_sigma(6.0) == approx(want, rel=0.03)
    assert J.lumped_sigma(6.0, proxy_j1=False) == approx(
        (5 * 106.91e-16 + 1 * 127.13e-16) / 12, rel=1e-6)
    assert J.lumped_sigma(6.0, proxy_j1=False) / J.lumped_sigma(6.0) == \
        approx(0.5, rel=0.0, abs=0.03)
    assert J.lumped_sigma(1.0) == 0.0


def test_maxwellian_rate_against_a_direct_quadrature():
    """The table's k(Te) against a trapezoid of sigma v over a Maxwellian
    on a fine grid of its own."""
    p = J.process()
    for te in (1.0, 2.0, 5.8):
        e = np.linspace(0.0, 60.0, 120001)
        f = 2.0 / math.sqrt(math.pi) * te ** -1.5 * np.sqrt(e) * np.exp(-e / te)
        k = _trapezoid(f * p.sigma_at(e) * xsec.V_CM_S * np.sqrt(e), e)
        assert float(J.rate_table()(te)) == approx(k, rel=2e-3), te


def test_reverse_is_detailed_balance_over_the_gap():
    for te in (0.5, 1.3, 3.0):
        ratio = (float(J.superelastic_table(12 / 36)(te))
                 / float(J.rate_table()(te)))
        assert ratio == approx(12 / 36 * math.exp(GAP / te), rel=2e-3), te


def test_each_level_against_a_direct_quadrature_and_its_own_balance():
    """v0_4: each 5s level's rate against a trapezoid of its own cross
    section (the J = 1 levels on their proxy's, at their own threshold),
    and its reverse by detailed balance over its own gap with its own
    weight, g_level / 36."""
    e = np.linspace(0.0, 60.0, 120001)
    for lv in LEVELS:
        p = J.level_process(lv)
        gap = e_form("Kr**") - e_form(f"Kr*({lv})")
        assert p.threshold == approx(gap), lv
        for te in (1.0, 2.0, 5.8):
            f = (2.0 / math.sqrt(math.pi) * te ** -1.5 * np.sqrt(e)
                 * np.exp(-e / te))
            k = _trapezoid(f * p.sigma_at(e) * xsec.V_CM_S * np.sqrt(e), e)
            assert float(J.level_rate_table(lv)(te)) == approx(
                k, rel=2e-3), (lv, te)
        for te in (0.5, 1.3, 3.0):
            ratio = (float(J.level_superelastic_table(lv, G5[lv] / 36)(te))
                     / float(J.level_rate_table(lv)(te)))
            assert ratio == approx(G5[lv] / 36 * math.exp(gap / te),
                                   rel=2e-3), (lv, te)


def test_measured_rate_against_kannaris_fit():
    """LIGHT v0_5's finding: the measured set lies within 40 % of Jacob and
    Mangano's fit from 1 to 6 eV, above it below 3 eV and below it
    above."""
    k = J.rate_table()
    for te, lo, hi in ((1.0, 1.2, 1.6), (2.0, 1.1, 1.4), (3.0, 1.0, 1.25),
                       (5.8, 0.8, 1.0)):
        r = float(k(te)) / float(J.kannari_rate(te))
        assert lo < r < hi, (te, r)


def test_kinetic_population_process_matches_the_tables():
    """The solver's forward and reverse population rates for a Maxwellian
    handed to it equal the tables (its grid against the quadrature), level
    by level (v0_4)."""
    m = Model(Params(eedf="kinetic"))
    g = KB.Grid(2400, 40.0, m.kin.procs, m.kin.superelastic, m.kin.pops)
    assert sorted(g.pops) == sorted(f"Kr*({lv})" for lv in LEVELS)
    for lv in LEVELS:
        for te in (1.0, 2.5, 5.8):
            F = g.maxwellian(te)
            kf, kr = g.pop_rates(f"Kr*({lv})", F)
            assert kf * 1e6 == approx(float(J.level_rate_table(lv)(te)),
                                      rel=3e-3), (lv, te)
            assert kr * 1e6 == approx(float(J.level_superelastic_table(
                lv, G5[lv] / 36)(te)), rel=3e-3), (lv, te)
    # the bands scale the process the solver sees: k_5s5p every level,
    # k_5s5p_j1 the J = 1 levels alone, and a level at 0 leaves the solver
    m2 = Model(Params(eedf="kinetic", rates=(("k_5s5p", 0.6),
                                             ("k_5s5p_j1", 0.5))))
    g2 = KB.Grid(2400, 40.0, m2.kin.procs, m2.kin.superelastic, m2.kin.pops)
    F = g2.maxwellian(2.5)
    for lv in LEVELS:
        x = 0.6 * (0.5 if lv in ("1s4", "1s2") else 1.0)
        assert g2.pop_rates(f"Kr*({lv})", F)[0] == approx(
            x * g.pop_rates(f"Kr*({lv})", F)[0], rel=1e-9), lv
    m3 = Model(Params(eedf="kinetic", rates=(("k_5s5p_j1", 0.0),)))
    assert sorted(q.name for q in m3.kin.pops) == ["Kr*(1s3)", "Kr*(1s5)"]


def test_the_pair_pays_and_returns_the_gap_and_captures_no_electron():
    """The exciting electron loses the ledger's gap and stays in the gas
    (the row's routing books from_e = gap, no capture), the reverse gives
    it back; the mean energy of the electrons that excite, above the gap
    and rising with Te, is a diagnostic of jung2006 only."""
    m = Model(Params(**FAST))
    byrid = {x.row.rid: x for x in m.compiled}
    for lv in LEVELS:
        gap = e_form("Kr**") - e_form(f"Kr*({lv})")
        c = J.level_capture_table(lv)
        assert float(c(1.16)) > gap and float(c(5.8)) > float(c(1.16))
        f, r = byrid[f"K76_{lv}"], byrid[f"K76r_{lv}"]
        assert f.cap is None and r.cap is None
        assert f.kind == r.kind == "pop_xsec"
        assert f.route["from_e"] == approx(gap)
        assert r.route["to_e"] == approx(gap)
        assert f.route["capture"] == 0


def test_kannari_switch_on_pure_krypton():
    """Params.k76 = 'kannari' runs LIGHT v0_4's fit out of every 5s level;
    in pure Kr it moves the 147-nm share by 2 % of itself (v0_3: 9 %), the
    measurement lying above the fit at the reference point's Te."""
    d = row_dict(run(Params(**FAST, k76="kannari")))
    for k, v in KANNARI_REF.items():
        assert d[k] == approx(v, rel=2e-3), k
    with pytest.raises(ValueError):
        Model(Params(k76="jacob"))


def test_reference_point_on_the_measured_cross_sections():
    """The pure-Kr reference point, pinned, and what the 5s <-> 5p cycle
    does to it: the electrons pump Kr* to Kr**, which the gas relaxes back
    with the gap as heat, so the stronger the pumping the less of the
    input reaches Kr2* and the 147-nm band. The band's ends bracket the
    default (v0_4: 39.5 % at 0.6, 28.8 % at 1.4), and with the cycle off
    69.4 % of the input is 147-nm light (v0_3: 51.8 %)."""
    r = run(Params(**FAST))
    d = row_dict(r)
    for k, v in V0_4_REF.items():
        assert d[k] == approx(v, rel=2e-3), k
    assert abs(accounting(r)["closure"]) < 1e-9
    lo = row_dict(run(Params(**FAST, rates=(("k_5s5p", 0.6),))))
    hi = row_dict(run(Params(**FAST, rates=(("k_5s5p", 1.4),))))
    assert hi["emission_vuv"] < d["emission_vuv"] < lo["emission_vuv"]
    assert hi["heat_reactions"] > d["heat_reactions"] > lo["heat_reactions"]
    off = row_dict(run(Params(**FAST, rates=(("k_5s5p", 0.0),))))
    assert off["emission_vuv"] == approx(69.397, rel=2e-3)
    assert off["emission_vuv"] > 2 * d["emission_vuv"]


def test_ledger_closes_with_the_measured_cycle_under_the_kinetic_eedf():
    """A short kinetic pulse with the measured pair inside the Boltzmann
    equation: the field carries the input and the two energy balances
    agree, as in LIGHT v0_5 (0.994 to 0.999 of the input power and 8e-6 in
    pure Kr, 2026-10-13). v0_2's 51 thresholds on the 600-cell grid give
    1.1e-5, the grid's own error: 3.2e-6 on 1200 cells (2026-10-14)."""
    r = Model(Params(eedf="kinetic", pulse_s=40e-9, t_end_s=45e-9,
                     n_save=91)).run()
    assert abs(accounting(r)["closure"]) < 1e-8
    on = [e for e in r.kinetic if e["on"] and not e.get("skipped")
          and e["t"] >= 20e-9]
    ratio = [e["power_ratio"] for e in on if np.isfinite(e["power_ratio"])]
    assert 0.98 < min(ratio) <= max(ratio) < 1.03
    assert max(abs(e["balance"]) for e in on) < 2e-5
