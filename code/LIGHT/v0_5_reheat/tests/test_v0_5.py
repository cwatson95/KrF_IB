"""v0_5: the 5s -> 5p electron excitation on Jung et al.'s measured cross
sections (jung2006.py) in place of Kannari's fit, under both closures,
with the fit kept as Params.k76 = 'kannari'. Every check here is a
falsification: the tables as printed, the lumped cross section's sum, its
Maxwellian rate against a direct quadrature, the reverse by detailed
balance, the kinetic solver's population process against the tables, the
fit switch reproducing v0_4, and the reference point pinned."""

import math

import numpy as np
import pytest
from pytest import approx

from v0_6 import boltzmann as KB
from v0_6 import jung2006 as J
from v0_6 import xsec
from v0_6.model import Model, Params, run
from v0_6.report import accounting, row_dict
from v0_6.species import e_form

GAP = e_form("Kr**") - e_form("Kr*")
# v0_4's reference point (two_t), which Params.k76 = 'kannari' must give
V0_4_REF = dict(emission_248=5.787, emission_vuv=4.53, tg_end_K=742.0)
# v0_5's reference point (two_t) on Jung's cross sections, 2026-10-09
V0_5_REF = dict(emission_248=5.684, emission_vuv=4.56, tg_end_K=741.8)
FAST = dict(eedf="two_t")


def test_tables_as_printed():
    """Spot values of Jung et al. 2006 Tables II and III and the Letter's
    Table I (the 3.5-eV column of Table IV)."""
    assert J.sigma_pair("1s5", "2p9", 6.0) == approx(52e-16)
    assert J.sigma_pair("1s5", "2p6", 2.0) == approx(8.6e-16)
    assert J.sigma_pair("1s5", "2p9", 300.0) == approx(6.0e-16)
    assert J.sigma_pair("1s3", "2p4", 4.0) == approx(65e-16)
    assert J.sigma_pair("1s3", "2p1", 8.0) == approx(0.23e-16)
    assert J.sigma_pair("1s3", "2p5", 3.5) == approx(0.4e-16)
    # below each pair's threshold the cross section is 0; the lowest
    # threshold is 1s5 -> 2p10 at 1.39 eV
    assert J.sigma_pair("1s5", "2p10", 1.3) == 0.0
    assert J.sigma_pair("1s5", "2p10", 1.5) > 0.0
    assert J.sigma_pair("1s3", "2p4", 1.5) == 0.0          # threshold 1.54
    # nothing measured out of 1s5 into 2p1
    assert J.sigma_pair("1s5", "2p1", 6.0) == 0.0
    # above the tables the 1s5 -> 2p9 shape, then ln E / E
    s300, s600 = J.sigma_pair("1s5", "2p8", [300.0, 600.0])
    assert s300 == approx(8.6e-16 * 6.0 / 51.0, rel=1e-9)
    assert 0.3 < s600 / s300 < 0.7


def test_lumped_cross_section_is_the_manifold_mean():
    """At 6 eV the sums of Tables II (1s5) and III (1s3) and the
    statistical mean with the J = 1 proxies; without the proxies the
    metastables alone carry 6 of 12."""
    assert J.level_total("1s5", 6.0) == approx(106.91e-16, rel=1e-6)
    assert J.level_total("1s3", 6.0) == approx(127.24e-16, rel=1e-6)
    want = (8 * 106.91e-16 + 4 * 127.24e-16) / 12
    # the proxies see the same final levels at their own thresholds, so
    # at 6 eV the shifted evaluation differs slightly from the plain sum
    assert J.lumped_sigma(6.0) == approx(want, rel=0.03)
    assert J.lumped_sigma(6.0, proxy_j1=False) == approx(
        (5 * 106.91e-16 + 1 * 127.24e-16) / 12, rel=1e-6)
    assert J.lumped_sigma(6.0, proxy_j1=False) / J.lumped_sigma(6.0) == approx(
        0.5, abs=0.03)
    assert J.lumped_sigma(1.0) == 0.0


def test_maxwellian_rate_against_a_direct_quadrature():
    """The table's k(Te) against a trapezoid of sigma v over a Maxwellian
    on a fine grid of its own."""
    p = J.process()
    for te in (1.0, 2.0, 5.8):
        e = np.linspace(0.0, 60.0, 120001)
        f = 2.0 / math.sqrt(math.pi) * te ** -1.5 * np.sqrt(e) * np.exp(-e / te)
        k = np.trapz(f * p.sigma_at(e) * xsec.V_CM_S * np.sqrt(e), e)
        assert float(J.rate_table()(te)) == approx(k, rel=2e-3), te


def test_reverse_is_detailed_balance_over_the_gap():
    for te in (0.5, 1.3, 3.0):
        ratio = float(J.superelastic_table(12 / 36)(te)) / float(J.rate_table()(te))
        assert ratio == approx(12 / 36 * math.exp(GAP / te), rel=2e-3), te


def test_measured_rate_against_kannaris_fit():
    """The finding: the measured set lies within 40 % of Jacob and
    Mangano's fit from 1 to 6 eV, above it below 3 eV and below it above."""
    k = J.rate_table()
    for te, lo, hi in ((1.0, 1.2, 1.6), (2.0, 1.1, 1.4), (3.0, 1.0, 1.25),
                       (5.8, 0.8, 1.0)):
        r = float(k(te)) / float(J.kannari_rate(te))
        assert lo < r < hi, (te, r)


def test_kinetic_population_process_matches_the_tables():
    """The solver's forward and reverse population rates for a Maxwellian
    handed to it equal the tables (its grid against the quadrature)."""
    m = Model(Params(eedf="kinetic"))
    g = KB.Grid(2400, 40.0, m.kin.procs, m.kin.superelastic, m.kin.pops)
    for te in (1.0, 2.5, 5.8):
        F = g.maxwellian(te)
        kf, kr = g.pop_rates("Kr*", F)
        assert kf * 1e6 == approx(float(J.rate_table()(te)), rel=3e-3), te
        assert kr * 1e6 == approx(float(J.superelastic_table(12 / 36)(te)),
                                  rel=3e-3), te
    # the band scales the process the solver sees
    m2 = Model(Params(eedf="kinetic", rates=(("k_5s5p", 0.5),)))
    g2 = KB.Grid(2400, 40.0, m2.kin.procs, m2.kin.superelastic, m2.kin.pops)
    F = g2.maxwellian(2.5)
    assert g2.pop_rates("Kr*", F)[0] == approx(0.5 * g.pop_rates("Kr*", F)[0],
                                               rel=1e-9)


def test_the_pair_pays_and_returns_the_gap_and_captures_no_electron():
    """The exciting electron loses the ledger's gap and stays in the gas
    (the row's routing books from_e = gap, no capture), the reverse gives
    it back; the mean energy of the electrons that excite, above the gap
    and rising with Te, is a diagnostic of jung2006 only."""
    c = J.capture_table()
    assert float(c(1.16)) > GAP and float(c(5.8)) > float(c(1.16))
    m = Model(Params(**FAST))
    byrid = {x.row.rid: x for x in m.compiled}
    assert byrid["K76"].cap is None and byrid["K76r"].cap is None
    assert byrid["K76"].kind == byrid["K76r"].kind == "pop_xsec"
    assert byrid["K76"].route["from_e"] == approx(GAP)
    assert byrid["K76r"].route["to_e"] == approx(GAP)
    assert byrid["K76"].route["capture"] == 0


def test_kannari_switch_reproduces_v0_4():
    d = row_dict(run(Params(**FAST, k76="kannari")))
    for k, v in V0_4_REF.items():
        assert d[k] == approx(v, rel=2e-3), k
    with pytest.raises(ValueError):
        Model(Params(k76="jacob"))


def test_reference_point_on_the_measured_cross_sections():
    r = run(Params(**FAST))
    d = row_dict(r)
    for k, v in V0_5_REF.items():
        assert d[k] == approx(v, rel=2e-3), k
    assert abs(accounting(r)["closure"]) < 1e-9
    # the band's two ends bracket the default on the 248-nm share
    lo = row_dict(run(Params(**FAST, rates=(("k_5s5p", 0.5),))))
    hi = row_dict(run(Params(**FAST, rates=(("k_5s5p", 1.4),))))
    assert lo["emission_248"] < d["emission_248"] < hi["emission_248"]
    # the cycle off (0) is v0_3's share, as in v0_4
    off = row_dict(run(Params(**FAST, rates=(("k_5s5p", 0.0),))))
    assert off["emission_248"] < lo["emission_248"]


def test_ledger_closes_with_the_measured_cycle_under_the_kinetic_eedf():
    """A short kinetic pulse with the measured pair inside the Boltzmann
    equation: the field carries the input and the two energy balances
    agree, as v0_4's pseudo cross section did."""
    r = Model(Params(eedf="kinetic", pulse_s=40e-9, t_end_s=45e-9,
                     n_save=91)).run()
    assert abs(accounting(r)["closure"]) < 1e-8
    on = [e for e in r.kinetic if e["on"] and not e.get("skipped")
          and e["t"] >= 20e-9]
    ratio = [e["power_ratio"] for e in on if np.isfinite(e["power_ratio"])]
    assert 0.98 < min(ratio) <= max(ratio) < 1.03
    assert max(abs(e["balance"]) for e in on) < 1e-5
