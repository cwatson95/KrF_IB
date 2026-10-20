"""beam.py gates: stopping power, emergent W, the conserving identity."""

import numpy as np
import pytest

from REBUILD.MODULES import beam as B
from REBUILD.MODULES.cross_sections import CrossSectionSet
from REBUILD.MODULES.units import N_AVOGADRO, gas_density_m3

N1 = gas_density_m3(1.0, 300.0)
_MOLAR = {"He": 4.0026, "Ne": 20.18, "Ar": 39.948, "Kr": 83.80,
          "F2": 37.9968}


def _mass_stopping(gas):
    S = B.stopping_power({gas: N1}, 6.0e5)
    rho = N1 * _MOLAR[gas] / N_AVOGADRO / 1e6
    return S * 1e-8 / rho


def test_stopping_vs_estar_band():
    """PIC_MCC_3D's measured Bethe/Moller-vs-ESTAR band is 2-6 % at
    600 keV; the same formula must land in the same band here."""
    for gas, ref in (("He", 1.94), ("Ne", 1.71), ("Ar", 1.52),
                     ("Kr", 1.34), ("F2", 1.63)):
        assert abs(_mass_stopping(gas) / ref - 1.0) < 0.08, gas


def test_csda_range_matches_pic_mcc_3d():
    # PIC_MCC_3D CHANGELOG [0.4.0]: 1.94 m at 600 keV, 1 atm Ar
    r = B.csda_range({"Ar": N1}, 6.0e5)
    assert abs(r - 1.94) < 0.05
    # exact 1/p scaling
    r3 = B.csda_range({"Ar": 3 * N1}, 6.0e5)
    assert abs(r3 * 3.0 / r - 1.0) < 1e-9


def test_emergent_w_scored_not_tuned():
    """W comes OUT; the misses vs literature are pinned so they cannot
    drift silently in either direction (the PIC_MCC_3D convention: its Ar
    +5 % / Kr +8 % are the expected equal-sharing Spencer-Fano class)."""
    W_ar = B.differential_W(CrossSectionSet({"Ar": N1}))
    W_kr = B.differential_W(CrossSectionSet({"Kr": N1}))
    assert 26.0 < W_ar < 29.5          # lit 26.4: +5% class
    assert 24.5 < W_kr < 28.0          # lit 24.2: +8% class
    W_ne = B.differential_W(CrossSectionSet({"Ne": N1}))
    assert 38.0 < W_ne < 47.0          # lit 36.6: +16% MISS, stated


def test_prompt_floor_banks_subexcitation():
    """The 0.1.12 fluid-Te fix: the prompt cascade ends at the noble-gas
    inelastic floor. The recursion is TIME-FREE, so below that floor it
    funneled 100% of the flux through the trace-F2 channels however slow
    their absolute collision rate (~ns per event at 0.1-0.5% F2 — the
    whole pulse, not prompt), pre-spending the Platzman sub-excitation
    share at t=0 while the network's own e+F2 rows re-treat the same
    pathway by rate. MEASURED before the fix: eps_sec 1.81 eV/pair at
    0.5% F2 (vs ~4.3 pure Ar), bulk Te x0.26 and F- x0.003 vs the NRL
    benchmark (scorecard.electra_population_rows). Pure-noble cascades
    must be bit-identical — the W anchors and the SF-vs-MC cross-check
    never move."""
    n = 2.1e25
    cs = CrossSectionSet({"Ar": n})
    U0, t0 = B.spencer_fano_per_process(cs, 550.0)
    U1, t1 = B.spencer_fano_per_process(cs, 550.0, prompt_floor_eV=11.0)
    assert np.allclose(U0, U1) and abs(t0 - t1) < 1e-9
    y = B.spencer_fano_yields({"Ar": n * 0.995, "F2": n * 0.005})
    assert y.eps_sec_eV > 3.5           # was 1.81 with the prompt drain
    assert "F2(a3)" not in y.per_pair   # the sub-floor state, un-pre-spent
    assert abs(y.identity_residual_eV()) < 1e-9
    # 0.1.25: the only F2 channels above the noble floor (C1/H1 at
    # 11.57/13.08 eV) book PROMPT dissociation — 2 F + fragment heat —
    # never a bound F2(A1) that nothing in the Kannari-chemistry modes
    # decays. subfloor now carries exactly that fragment heat.
    assert "F2(A1)" not in y.per_pair
    cs2 = CrossSectionSet({"Ar": n * 0.995, "F2": n * 0.005})
    floor = min(float(cs2.thresholds[k]) for k, kind in enumerate(cs2.kinds)
                if kind in ("excitation", "ionization")
                and cs2.gases[k] == "Ar")
    U, _term = B.spencer_fano_per_process(cs2, B.default_cascade_energy(cs2),
                                          prompt_floor_eV=floor)
    u_ion = sum(float(U[k]) for k, kind in enumerate(cs2.kinds)
                if kind == "ionization")
    heat = sum(float(U[k]) / u_ion * (float(cs2.thresholds[k]) - 1.6)
               for k, kind in enumerate(cs2.kinds)
               if kind == "excitation" and cs2.gases[k] == "F2"
               and float(cs2.thresholds[k]) >= 2.5)
    assert abs(y.subfloor_eV - heat) < 1e-12


def test_conserving_identity_both_routes():
    dens = {g: f * gas_density_m3(3.0, 300.0)
            for g, f in (("Ar", 0.949), ("Kr", 0.05), ("F2", 0.001))}
    y = B.spencer_fano_yields(dens)
    assert abs(y.identity_residual_eV()) < 1e-9
    assert y.eps_sec_eV > 0 and y.subfloor_eV >= 0
    yw = B.w_table_yields({"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
    assert abs(yw.identity_residual_eV()) < 1e-9
    assert y.mode == "spencer_fano" and yw.mode == "w_table"


def test_ceiling_refusal_has_teeth():
    """PIC_MCC_3D Known-bad: above the tables every sigma reads zero and a
    cascade silently returns nothing. The refusal must raise."""
    cs = CrossSectionSet({"Ar": N1})
    with pytest.raises(ValueError, match="ceiling"):
        B.spencer_fano_per_process(cs, cs.energy_ceiling_eV * 2.0)


def test_deposition_profile_responds_to_density():
    """The MESH.txt density-adaptivity rung: the march uses LOCAL density,
    so a denser slab absorbs more."""
    spec = B.BeamSpec()
    lo = {"Ar": N1}
    hi = {"Ar": 3 * N1}
    prof = B.deposition_profile(spec, [lo, hi, lo, hi], 0.05)
    assert prof[1] > 2.5 * prof[0]
    assert prof.sum() <= 1.0 + 1e-12


def test_beam_sources_exhaust_p_dep():
    dens = {"Ar": 3 * N1}
    y = B.spencer_fano_yields(dens)
    p = np.array([1e10, 5e10])
    prod, p_e, p_gas = B.beam_sources(y, p)
    from REBUILD.MODULES.species import e_form_of
    from REBUILD.MODULES.units import EV
    pairs = p / (y.W_eV * EV)
    spend = sum(n * e_form_of(s) for s, n in y.per_pair.items())
    total = pairs * spend * EV + p_e + p_gas
    assert np.allclose(total, p, rtol=1e-12)


def test_pulse_envelope():
    s = B.BeamSpec(t_flat_ns=60.0, t_rise_ns=5.0, j_A_cm2=10.0,
                   foil_transmission=0.5)
    assert s.current_density_A_m2(-1e-9) == 0.0
    assert abs(s.current_density_A_m2(2.5e-9) / (10e4 * 0.5) - 0.5) < 1e-12
    assert s.current_density_A_m2(30e-9) == 10e4 * 0.5
    assert s.current_density_A_m2(80e-9) == 0.0
