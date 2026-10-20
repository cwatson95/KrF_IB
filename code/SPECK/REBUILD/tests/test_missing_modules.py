"""The four board [gap] closures (2026-08-04): in-loop ASE, transport,
the MESH.txt refinement rung, Chen-figure scoring."""

import numpy as np
import pytest

from REBUILD.MODULES.engine import Device
from REBUILD.MODULES.params import engine_parameters
from REBUILD.MODULES.units import gas_density_m3


def _krf(**over):
    p = engine_parameters()
    d = {"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
         "t_end_ns": 120.0, "nx": 1, "ny": 1, "nz": 1, "Lx_m": 0.0141,
         "Ly_m": 0.0141, "Lz_m": 0.75, "beam_axis": 0}
    d.update(over)
    p.set_many(d, source="problem")
    return p


# --- in-loop ASE (S6) ------------------------------------------------------
def test_ase_identity_holds_with_sink(loki_ar_fill):
    """The every-eV-routed identity with the ASE sink active: same
    mechanical routing as stimulated emission, photons to escape."""
    from REBUILD.MODULES.chemistry import CompiledNetwork
    from REBUILD.MODULES.network import build_network, mode_optics
    from REBUILD.MODULES.species import (KRF_MODE_SPECIES, StateLayout,
                                         e_form_of)
    from REBUILD.MODULES.units import EV, photon_energy_J
    rows = build_network("krf", loki_ar_fill.available_buckets())
    lay = StateLayout(KRF_MODE_SPECIES, 2)
    cn = CompiledNetwork(rows, lay, loki_ar_fill, mode_optics("krf"))
    cn.set_photon_energy(photon_energy_J(248.0) / EV)
    rng = np.random.default_rng(3)
    dens = np.full((lay.n_species, 2), 1e16)
    dens[lay.index["Ar"]] = 7e25
    dens[lay.index["Kr"]] = 3.7e24
    dens[lay.index["F2"]] = 7.4e22
    dens[lay.index["KrF*(B)"]] = 1e19
    te = np.array([0.5, 1.0]); tg = np.array([300.0, 400.0])
    ne = lay.electron_density(dens)
    out = cn.rhs(dens, te, tg, ne, 1e17, 7.4e25, ase_rate=3e7)
    e_form = np.array([e_form_of(s) for s in lay.species])
    d_pot = (e_form[:, None] * out["ddens"]).sum(axis=0) * EV
    hv = cn.hv_eV
    tot = (d_pot + out["du_e"] + out["dgas"] + out["rad_escape"]
           + out["r_se"] * hv * EV + (out["r_ste"] - out["r_x"]) * hv * EV)
    scale = np.abs(d_pot).max() + np.abs(out["dgas"]).max()
    assert np.abs(tot).max() / scale < 1e-12


def test_ase_clamps_gain_at_high_gain_and_vanishes_at_low():
    """Direction gate: at Kannari-class gain the ASE sink lowers peak g0;
    at low gain on/off agree to integrator noise."""
    hi = {"pressure_atm": 2.5, "beam_kV": 660.0, "beam_j_A_cm2": 300.0,
          "beam_flat_ns": 50.0, "Lx_m": 0.04, "Ly_m": 0.04, "Lz_m": 0.40,
          "t_end_ns": 90.0}
    g = {}
    for on in (True, False):
        r = Device(_krf(ase_inloop=on, **hi)).run(n_save=19)
        g[on] = r["g0_m1"].max()
    assert g[True] < g[False] * 0.995
    lo = {"beam_j_A_cm2": 0.5, "t_end_ns": 60.0}
    e = {}
    for on in (True, False):
        r = Device(_krf(ase_inloop=on, **lo)).run(n_save=13)
        e[on] = r["ne_m3"].max()
    assert abs(e[True] / e[False] - 1.0) < 1e-6


# --- inter-cell transport --------------------------------------------------
def test_transport_conserves_exactly_and_relaxes_gradients():
    p = _krf(transport=True, nx=6, ny=1, nz=2, t_end_ns=40.0,
             beam_j_A_cm2=1.0)
    dev = Device(p)
    rng = np.random.default_rng(5)
    dens = np.abs(rng.normal(1e18, 5e17, (dev.lay.n_species,
                                          dev.n_cells)))
    term = dev._transport_term(dens)
    # conservation: machine-exact (uniform volumes here)
    assert np.abs(term.sum(axis=1)).max() < 1e-6 * np.abs(term).max()
    # a uniform field does not move
    assert np.abs(dev._transport_term(np.full_like(dens, 3e18))).max() == 0
    # sinusoid eigenmode decays at ~ -D k^2 (O(h^2) discretization)
    nx = 6
    k = 2 * np.pi / dev.p["Lx_m"]
    prof = 1e18 * (1.0 + 0.5 * np.sin(k * (np.arange(nx) + 0.5)
                                      * dev.d_cell[0]))
    d1 = np.tile(prof.repeat(2), (dev.lay.n_species, 1))
    t1 = dev._transport_term(d1)
    i = dev.lay.index["Ar"]
    num = t1[i] / (d1[i] - 1e18)
    lam = -np.median(num[np.abs(d1[i] - 1e18) > 1e16])
    lam_ref = dev._diff0 * (2.0 / dev.d_cell[0] ** 2
                            * (1.0 - np.cos(k * dev.d_cell[0])))
    assert abs(lam / lam_ref - 1.0) < 1e-6


def test_transport_off_is_default_and_gauged():
    p = _krf()
    dev = Device(p)
    assert dev._transport is False


# --- MESH.txt refinement rung ---------------------------------------------
def test_equal_deposition_refinement_resolves_the_bragg_stop():
    """A 150 kV beam STOPS inside 3 cm of 3 atm: the refined mesh
    equidistributes the deposit (max/min ~1) where the uniform grid piles
    x3+ into the Bragg cells; both totals match the fine march to 1 %."""
    from REBUILD.MODULES.beam import (BeamSpec, deposition_profile,
                                      equal_deposition_edges)
    dens = {g: f * gas_density_m3(3.0, 300.0)
            for g, f in (("Ar", 0.949), ("Kr", 0.05), ("F2", 0.001))}
    spec = BeamSpec(voltage_kV=100.0, axis=0)     # CSDA 4.05 cm < L: stops
    L, n = 0.06, 12
    edges = equal_deposition_edges(spec, dens, L, n)
    assert edges[0] == 0.0 and abs(edges[-1] - L) < 1e-12
    w = np.diff(edges)
    assert (w > 0).all() and abs(w.sum() - L) < 1e-12
    fr = deposition_profile(spec, [dens] * n, w)
    fu = deposition_profile(spec, [dens] * n, L / n)
    fine = deposition_profile(spec, [dens] * 600, L / 600).sum()
    # measured 2026-08-04: refined 3.5 (the flat-F tail cell), uniform 2e11
    assert fr.max() / fr.min() < 5.0
    assert fu.max() / max(fu.min(), 1e-12) > 100.0
    assert abs(fr.sum() / fine - 1.0) < 0.01
    assert abs(fu.sum() / fine - 1.0) < 0.01


def test_engine_runs_on_a_refined_axis_and_closes():
    from REBUILD.MODULES.beam import BeamSpec, equal_deposition_edges
    dens = {g: f * gas_density_m3(3.0, 300.0)
            for g, f in (("Ar", 0.949), ("Kr", 0.05), ("F2", 0.001))}
    p = _krf(nx=6, beam_kV=100.0, Lx_m=0.06, t_end_ns=60.0)
    edges = equal_deposition_edges(BeamSpec(voltage_kV=100.0, axis=0),
                                   dens, 0.06, 6)
    dev = Device(p, axis_edges=edges)
    assert abs(dev.cell_vols.sum()
               - 0.06 * 0.0141 * 0.75) < 1e-12
    r = dev.run(n_save=13)
    assert r["success"]
    worst = max(abs(c["residual_rel"]) for c in r["closure"][1:])
    assert worst < 5e-4


def test_refined_axis_guards():
    from REBUILD.MODULES.beam import BeamSpec, equal_deposition_edges
    dens = {"Ar": gas_density_m3(3.0, 300.0)}
    edges = equal_deposition_edges(BeamSpec(voltage_kV=100.0, axis=0),
                                   dens, 0.06, 6)
    with pytest.raises(ValueError, match="BEAM axis"):
        Device(_krf(drive="discharge", nx=6, Lx_m=0.06),
               axis_edges=edges)
    with pytest.raises(ValueError, match="monotone|span"):
        Device(_krf(nx=6, Lx_m=0.06), axis_edges=edges[:-1])


# --- Chen-figure scoring (S7) ---------------------------------------------
def test_chen_figure6_scoring_renders_and_pins_the_gap():
    """0.1.25: the anchor renders the MEASURED attribution (the 0.1.5..24
    'x3-5 formation gap' suspect list is kept as history, marked wrong),
    the as-printed lane's recorded non-integrability, the below-table E/N
    fraction, and the corrected-lane n_e ratio (x0.79 class vs Fig. 10a's
    linear-axis 6.5e20) — each pinned against silent drift."""
    from REBUILD.MODULES.scorecard import chen_figure6_rows
    rows = chen_figure6_rows(t_end_ns=200.0)
    txt = "\n".join(rows)
    assert "VERDICT" in txt and "WRONG RANKING" in txt
    assert "as-printed lane" in txt and "does NOT integrate" in txt
    assert "BELOW the floor" in txt
    ne_row = [r for r in rows if "n_e" in r][0]
    ratio = float(ne_row.split("(x")[1].split(")")[0])
    # the corrected-lane miss class: visible, pinned against silent drift
    assert 0.6 < ratio < 1.0, ne_row


# --- validity gauges (the ESSim.stability analog) -------------------------
def test_validity_gauges_present_and_sane(krf_result):
    v = krf_result["validity"]
    assert 0.0 < v["ionization_frac_max"] < 0.01     # static-gas assumption
    assert 0.0 <= v["f2_burn_frac"] < 0.5
    assert 0.0 <= v["te_in_table_frac"] <= 1.0
    lo, hi = v["tier_eps_range_eV"]
    assert 0.0 < lo < hi
    assert "validity:" in krf_result["regime"]
