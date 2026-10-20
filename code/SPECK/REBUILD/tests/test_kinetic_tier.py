"""kinetic_tier.py gates — the S4 coupling: data in, comparison out.

The band pins below are REGRESSION PINS on the shipped export
(kinetic_tables/pic_mcc_3d_swarm.npz, provenance in its meta), not physics
claims: they hold the recorded PIC_MCC_3D-vs-LoKI-MC agreement fixed so it
cannot drift silently. Regenerating the export legitimately moves them.
"""

import numpy as np
import pytest

from REBUILD.MODULES.kinetic_tier import (KineticRates, compare_tiers,
                                          format_comparison)


@pytest.fixture(scope="module")
def kin():
    return KineticRates.load_default()


def test_export_loads_with_provenance(kin):
    assert kin.meta["source"].startswith("PIC_MCC_3D")
    assert "mixture" in kin.meta and "seed" in kin.meta
    assert {"Ar+", "Kr+", "F-", "el:Ar"} <= set(kin.buckets)


def test_lookup_vectorized(kin):
    te = np.array([3.0, 4.0, 5.0])
    v = kin.k("Ar+", te)
    assert v.shape == te.shape and (v > 0).all()
    assert kin.k("Ar+", 4.0) == v[1]


def test_comparison_uses_overlap_only(kin, loki_ar_fill):
    rep = compare_tiers(kin, loki_ar_fill)
    for b, r in rep["buckets"].items():
        lo, hi = r["overlap_eps_eV"]
        assert lo < hi, b
        assert lo >= 3.0, "comparison reached below the export's table"


def test_recorded_agreement_bands(kin, loki_ar_fill):
    """Elastic x0.87-0.90 and attachment x0.95 medians measured 2026-08-04;
    excitation medians sit in the pic1d x1.6 class. Pins are generous so
    statistics can move them, but a x10 silent drift fails."""
    rep = compare_tiers(kin, loki_ar_fill)["buckets"]
    assert 0.5 < rep["el:Ar"]["ratio_median"] < 2.0
    assert 0.5 < rep["F-"]["ratio_median"] < 2.0
    assert 0.5 < rep["Kr*"]["ratio_median"] < 5.0
    # RESOLVED 2026-08-04: the x100 vib "finding" was a LOKI-side bucket
    # bug (the F2 Effective row fell through to vib — loki_tier._bucket
    # docstring). Post-fix the bucket sits in-band: measured x0.84
    # [0.51, 1.15]. This pin now guards the FIX.
    assert 0.3 < rep["vib"]["ratio_median"] < 2.5


def test_report_formats(kin, loki_ar_fill):
    txt = format_comparison(compare_tiers(kin, loki_ar_fill))
    assert "never merged" in txt and "el:Ar" in txt


def test_missing_export_raises_with_instructions(tmp_path):
    with pytest.raises(FileNotFoundError, match="export_kinetic"):
        KineticRates(tmp_path / "nope.npz")


def test_w_cascade_shipped_with_provenance(kin):
    """The SECOND W route (the PIC_MCC_3D gap-closure, 2026-08-04): per-gas
    MC-cascade results ride the export with full provenance."""
    wc = kin.w_cascade
    assert {"Ar", "Kr", "Ne"} <= set(wc)
    for gas, rec in wc.items():
        assert rec["converged"], gas
        assert 15.0 < rec["W_eV"] < 60.0, gas
        assert rec["E0_eV"] > 400.0, gas
        assert rec["n_electrons"] >= 150, gas
        assert 0.0 < rec["stat_err_rel"] < 0.05, gas


def test_relaxation_data_shipped_and_taper_band(kin, loki_ar_fill):
    """L3 (2026-08-04, deepened 0.1.13): the export carries
    relaxing-population rates below the swarm tables; vs them the engine
    taper sat x0.8-1.9 (elastic/vib) and x3.3 HIGH (F- at 1.8 eV) on the
    stalled export. The 0.1.13 deepening (same seed, 10x rounds, floors
    1.26/1.76 -> 0.47/0.66 eV) extends the band down through the Ar
    Ramsauer minimum, where the taper now sits x0.42 LOW for el:Ar —
    both directions pinned so the measured band cannot drift silently
    (a regenerated export moves it legitimately)."""
    from REBUILD.MODULES.kinetic_tier import taper_vs_relaxation
    assert {"F-", "el:Ar", "el:Kr"} <= set(kin.relax)
    rows = taper_vs_relaxation(kin, loki_ar_fill)
    by = {}
    for b, e, km, kt, ratio, err in rows:
        by.setdefault(b, []).append(ratio)
    assert 2.0 < max(by["F-"]) < 6.0          # the x3.3 finding (1.76 eV)
    assert min(by["F-"]) < 0.5     # the deepened bins: taper LOW where
    #                                measured DA rises toward low eps
    for b in ("el:Ar", "el:Kr"):
        for r in by[b]:
            assert 0.35 < r < 3.0  # was 0.5 floor; 0.42 at the 0.47-eV
            #                        Ramsauer bin (deepened 0.1.13)


def test_bridge_consumes_measured_bins_pure(kin, loki_ar_fill):
    """The L3 closure: below the table floor, at a measured bin outside the
    blend band, the bridged rate IS the measurement (the x3.3-high taper
    value is gone)."""
    from REBUILD.MODULES.kinetic_tier import RelaxationBridged
    br = RelaxationBridged(loki_ar_fill, kin)
    assert {"F-", "el:Ar", "el:F2", "el:Kr", "vib"} <= set(br.bridged_buckets)
    eps, k, err = kin.relax["F-"]
    # select with the bridge's OWN filter (err < MAX_REL_ERR) — post the
    # 0.1.13 deepened export the raw curve reaches 0.47 eV but that bin
    # is a single event (err 1.0) the bridge rightly drops; the first
    # USABLE bin is 0.66 eV (was 1.76 at the Ramsauer-stalled export)
    good = (np.isfinite(k) & (k > 0) & np.isfinite(err)
            & (err < RelaxationBridged.MAX_REL_ERR))
    e0, k0 = float(eps[good][0]), float(k[good][0])
    assert e0 < 1.9                                       # truly sub-floor
    assert abs(br.k("F-", e0 / 1.5) / k0 - 1.0) < 1e-9
    taper = float(loki_ar_fill.k("F-", e0 / 1.5))
    # direction-agnostic x2 gate: at the old 1.76-eV floor the stand-in
    # was x3.3 HIGH; at the deepened 0.66-eV floor it is x10-class LOW
    # (measured DA RISES toward low eps while the taper falls) — either
    # way the measurement, not the stand-in, must be what the bridge
    # returns (asserted exactly above)
    assert taper / k0 > 2.0 or k0 / taper > 2.0


def test_bridge_error_filter_drops_uninformative_bins(kin, loki_ar_fill):
    """err >= 0.5 bins carry no information and are excluded; sub-floor
    Kr+ is not bridged at all (no sub-threshold events exist). History of
    the vib floor: 1.76 eV at the Ramsauer-stalled export (its 1.26-eV
    bin was a 100 %-err single event) -> 0.66 eV after the 0.1.13
    deepening (same seed, 10x rounds: that bin now carries 3 % stats and
    the floor moved two bins down). The F- single-event 0.47-eV bin is
    the CURRENT filter demonstration: raw coverage reaches it, bridged
    coverage rightly does not."""
    from REBUILD.MODULES.kinetic_tier import RelaxationBridged
    br = RelaxationBridged(loki_ar_fill, kin)
    cov = br.sub_table_coverage()
    assert abs(cov["vib"][0] - 0.66) < 0.05
    assert cov["F-"][0] > 0.5          # the 0.47 single-event bin dropped
    assert "Kr+" not in cov


def test_bridge_identity_above_floor_and_delegation(kin, loki_ar_fill):
    """At or above the table floor the bridge returns the base tier's value
    exactly (scalar and array), and every non-k attribute delegates —
    the circuit path (k_by_en) and the engine's gauges see the base tier."""
    from REBUILD.MODULES.kinetic_tier import RelaxationBridged
    br = RelaxationBridged(loki_ar_fill, kin)
    te = np.array([1.5, 2.0, 3.0, 5.0])                  # eps 2.25-7.5, in
    for b in ("F-", "el:Ar", "Kr+", "Ar*"):
        assert np.array_equal(br.k(b, te), loki_ar_fill.k(b, te)), b
    assert br.k("F-", 2.0) == loki_ar_fill.k("F-", 2.0)
    assert br.k_by_en("F-", 15.0) == loki_ar_fill.k_by_en("F-", 15.0)
    assert br.available_buckets() == loki_ar_fill.available_buckets()
    assert br.en_range_td == loki_ar_fill.en_range_td


def test_bridge_continuity_for_bdf(kin, loki_ar_fill):
    """k(Te) crossing the floor and the whole sub-table region must be
    continuous (C0) — a jump would stall the BDF error controller. Gate:
    on a dense log grid, adjacent samples never differ by more than the
    grid's worst-case smooth slope."""
    from REBUILD.MODULES.kinetic_tier import RelaxationBridged
    br = RelaxationBridged(loki_ar_fill, kin)
    te = np.geomspace(0.05, 3.0, 800)
    for b in br.bridged_buckets:
        k = br.k(b, te)
        assert (k > 0).all(), b
        step = np.abs(np.diff(np.log(k)))
        assert step.max() < 0.12, (b, float(step.max()))


def test_w_two_routes_agree(kin):
    """REBUILD's deterministic Spencer–Fano vs PIC_MCC_3D's Monte-Carlo
    null-collision cascade — independent code paths, same E0, same LXCat
    data. Measured 2026-08-04: Ar 0.19 %, Kr 0.04 %, Ne 0.02 % with MC
    statistics ~1.2 %. Gate at 5 % (~3–4 sigma of the MC error) so
    regenerated exports with fresh seeds cannot flap it."""
    from REBUILD.MODULES.beam import spencer_fano_per_process
    from REBUILD.MODULES.cross_sections import CrossSectionSet
    from REBUILD.MODULES.units import gas_density_m3
    n = gas_density_m3(1.0, 300.0)
    for gas, rec in kin.w_cascade.items():
        cs = CrossSectionSet({gas: n})
        U, _ = spencer_fano_per_process(cs, rec["E0_eV"])
        ion = [k for k, kind in enumerate(cs.kinds) if kind == "ionization"]
        W_sf = rec["E0_eV"] / float(U[ion].sum())
        assert abs(W_sf / rec["W_eV"] - 1.0) < 0.05, (gas, W_sf, rec["W_eV"])
