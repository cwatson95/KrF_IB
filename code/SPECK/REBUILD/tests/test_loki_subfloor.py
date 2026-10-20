"""The measured sub-0.5 Td LOKI extension (0.1.39): the tables load, their
E/N grid is disjoint from the shipped sets, merging lowers the <eps> floor
into the beam-mode operating band, the two wrong-way channels are now
MEASURED rather than tapered, and the opt-in is inert by default.
"""

import numpy as np
import pytest

from REBUILD.MODULES.loki_tier import (RUN_NAMES_AR_FILL,
                                       RUN_NAMES_AR_FILL_SUBFLOOR,
                                       RUN_NAMES_NE_FILL,
                                       RUN_NAMES_NE_FILL_SUBFLOOR,
                                       LokiRates, rates_for_mode)

AR = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
NE = {"Ne": 0.949, "Kr": 0.05, "F2": 0.001}

#: measured 2026-09-04, LoKI-MC boltzmannMC 2e4 electrons (PROVENANCE.txt
#: in each table dir). These literals are the run's own output — a gate on
#: the shipped data, not a tuned number.
AR_POINTS = {0.05: 0.7612, 0.1: 1.072, 0.2: 1.427, 0.3: 1.660}
NE_POINTS = {0.1: 0.8906, 0.3: 1.740, 0.5: 2.231}


@pytest.mark.parametrize("names,pts", [(RUN_NAMES_AR_FILL_SUBFLOOR, AR_POINTS),
                                       (RUN_NAMES_NE_FILL_SUBFLOOR, NE_POINTS)])
def test_subfloor_tables_carry_the_measured_points(names, pts):
    lok = LokiRates(names)
    assert lok.en_range_td[0] == pytest.approx(min(pts), rel=1e-9)
    assert lok.en_range_td[1] == pytest.approx(max(pts), rel=1e-9)
    for en, eps in pts.items():
        assert lok.mean_energy_at_en(en) == pytest.approx(eps, rel=2e-3), en
    for b in ("vib", "F-"):
        assert b in lok.curves, b


@pytest.mark.parametrize("base,extra", [(RUN_NAMES_AR_FILL, RUN_NAMES_AR_FILL_SUBFLOOR),
                                        (RUN_NAMES_NE_FILL, RUN_NAMES_NE_FILL_SUBFLOOR)])
def test_grids_are_disjoint_so_the_overlap_guard_stays_quiet(base, extra):
    """LokiRates sums processes into a bucket WITHIN one E/N; a second
    directory tabulating the same E/N would silently double it. The guard
    raises on overlap — merging must not trip it, and re-adding the same
    set must."""
    merged = LokiRates(list(base) + list(extra))
    assert merged.en_range_td[0] == LokiRates(extra).en_range_td[0]
    # the SAME set listed twice used to pass the cross-directory guard and
    # silently double every bucket (x2.0000, measured 0.1.39) — now refused
    with pytest.raises(ValueError, match="duplicate"):
        LokiRates(list(extra) + list(extra))
    with pytest.raises(ValueError, match="duplicate"):
        LokiRates(list(base) + list(base))


@pytest.mark.parametrize("base,extra,mix,old,new", [
    (RUN_NAMES_AR_FILL, RUN_NAMES_AR_FILL_SUBFLOOR, AR, 1.980, 0.7612),
    (RUN_NAMES_NE_FILL, RUN_NAMES_NE_FILL_SUBFLOOR, NE, 3.063, 0.8906)])
def test_merging_lowers_the_floor_into_the_operating_band(base, extra, mix,
                                                          old, new):
    def floor(lok):
        return min(float(np.exp(lok.curves[b][0][0])) for b in lok.buckets)
    assert floor(LokiRates(list(base))) == pytest.approx(old, rel=2e-3)
    assert floor(LokiRates(list(base) + list(extra))) == pytest.approx(new, rel=2e-3)
    # the settled beam-mode Te (0.04-0.63 eV) is still below even the new
    # floor, so the taper is reduced in reach, not eliminated
    assert new / 1.5 > 0.4


@pytest.mark.parametrize("mix,extra", [(AR, RUN_NAMES_AR_FILL_SUBFLOOR),
                                       (NE, RUN_NAMES_NE_FILL_SUBFLOOR)])
def test_the_two_wrong_way_channels_are_now_measured_not_tapered(mix, extra):
    """The 0.1.38 audit found F- and vib rise below the old floor while the
    taper falls. At the new floor they must come from the table, and be
    HIGHER than the old set's extrapolation there."""
    old = rates_for_mode("krf", mix)
    new = rates_for_mode("krf", mix, subfloor=True)
    te = float(np.exp(min(new.curves[b][0][0] for b in new.buckets))) / 1.5
    for b in ("vib", "F-"):
        k_old, k_new = float(old.k(b, te)), float(new.k(b, te))
        assert k_new > 5.0 * k_old, (b, k_new, k_old)
    # a threshold bucket is not inflated by the extension
    assert float(new.k("Kr+", te)) <= float(old.k("Kr+", te)) * 1.05


def test_the_flag_is_off_by_default_and_inert():
    from REBUILD.MODULES.params import engine_parameters
    assert engine_parameters()["loki_subfloor"] is False
    a, b = rates_for_mode("krf", AR), rates_for_mode("krf", AR, subfloor=False)
    assert a.run_names == b.run_names
    assert a.en_range_td == b.en_range_td
    # arf_chen has no measured extension yet: the flag must be inert there
    c = rates_for_mode("arf_chen", {"Ne": 0.964, "Ar": 0.035, "F2": 0.001})
    d = rates_for_mode("arf_chen", {"Ne": 0.964, "Ar": 0.035, "F2": 0.001},
                       subfloor=True)
    assert c.run_names == d.run_names


def test_engine_runs_on_the_extension_and_closes():
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "t_end_ns": 20.0, "nx": 1,
                "ny": 1, "nz": 1, "Lx_m": 0.0141, "Ly_m": 0.0141,
                "Lz_m": 0.75, "beam_axis": 0}, source="problem")
    p.set("loki_subfloor", True, source="cli")
    res = Device(p).run(n_save=9)
    assert res["success"]
    assert max(abs(c["residual_rel"]) for c in res["closure"][1:]) < 5e-4
    # the reported in-table fraction must RISE against the same run without it
    p2 = engine_parameters()
    p2.set_many({"mode": "krf", "drive": "beam", "t_end_ns": 20.0, "nx": 1,
                 "ny": 1, "nz": 1, "Lx_m": 0.0141, "Ly_m": 0.0141,
                 "Lz_m": 0.75, "beam_axis": 0}, source="problem")
    ref = Device(p2).run(n_save=9)
    assert (res["validity"]["te_in_table_frac"]
            >= ref["validity"]["te_in_table_frac"])
