"""chemistry.py gates — the assembled RHS conserves what it must.

The two mechanical identities here are the structural answer to phantom11's
hand-written-RHS defect class: EVERY evaluation, on random states, must
conserve each chemical element exactly and route every eV to a named
destination. A hidden term, a flipped sign (the r115 class), or an unrouted
release breaks one of them at machine precision.
"""

import numpy as np
import pytest

from REBUILD.MODULES.chemistry import CompiledNetwork
from REBUILD.MODULES.network import build_network, mode_optics
from REBUILD.MODULES.species import (KRF_MODE_SPECIES, REGISTRY, StateLayout,
                                     atoms_of, e_form_of)
from REBUILD.MODULES.units import EV, photon_energy_J


@pytest.fixture(scope="module")
def compiled(loki_ar_fill):
    rows = build_network("krf", loki_ar_fill.available_buckets())
    lay = StateLayout(KRF_MODE_SPECIES, n_cells=4)
    cn = CompiledNetwork(rows, lay, loki_ar_fill, mode_optics("krf"))
    cn.set_photon_energy(photon_energy_J(248.0) / EV)
    return cn, lay


def _random_state(lay, seed=7, n_tot=7.4e25):
    rng = np.random.default_rng(seed)
    C = lay.n_cells
    dens = np.zeros((lay.n_species, C))
    dens[lay.index["Ar"]] = 0.949 * n_tot
    dens[lay.index["Kr"]] = 0.05 * n_tot
    dens[lay.index["F2"]] = 0.001 * n_tot
    for s in lay.species:
        if s not in ("Ar", "Kr", "F2"):
            dens[lay.index[s]] = rng.uniform(1e14, 1e19, C)
    te = rng.uniform(0.3, 3.0, C)
    tg = rng.uniform(300.0, 500.0, C)
    ne = lay.electron_density(dens)
    return dens, te, tg, ne


def test_elemental_mass_conserved_exactly(compiled):
    cn, lay = compiled
    dens, te, tg, ne = _random_state(lay)
    out = cn.rhs(dens, te, tg, ne, n_ph=1e18, ntot=7.4e25)
    elements = sorted({el for s in lay.species for el in atoms_of(s)})
    A = np.zeros((len(elements), lay.n_species))
    for j, s in enumerate(lay.species):
        for el, n in atoms_of(s).items():
            A[elements.index(el), j] = n
    imbalance = A @ out["ddens"]                      # atoms m^-3 s^-1
    scale = np.abs(A @ np.abs(out["ddens"])).max() + 1e-300
    assert np.abs(imbalance).max() / scale < 1e-12


def test_every_ev_routed_exactly(compiled):
    """Per-cell, per-evaluation energy identity (see chemistry.py):
    d(potential) + du_e + dgas + rad_escape + r_se*hv + (r_ste - r_x)*hv = 0
    """
    cn, lay = compiled
    hv = cn.hv_eV
    e_form = np.array([e_form_of(s) for s in lay.species])
    for n_ph in (0.0, 1e18):
        dens, te, tg, ne = _random_state(lay, seed=int(3 + n_ph % 7))
        out = cn.rhs(dens, te, tg, ne, n_ph=n_ph, ntot=7.4e25)
        d_pot = (e_form[:, None] * out["ddens"]).sum(axis=0) * EV
        photon = (out["r_ste"] - out["r_x"]) * hv * EV
        upper_rad = out["r_se"] * hv * EV
        total = (d_pot + out["du_e"] + out["dgas"] + out["rad_escape"]
                 + upper_rad + photon)
        scale = (np.abs(d_pot) + np.abs(out["du_e"]) + np.abs(out["dgas"])
                 + np.abs(out["rad_escape"]) + np.abs(photon)).max()
        assert np.abs(total).max() / scale < 1e-12


def test_vectorization_matches_per_cell(compiled, loki_ar_fill):
    cn, lay = compiled
    dens, te, tg, ne = _random_state(lay, seed=11)
    out = cn.rhs(dens, te, tg, ne, n_ph=3e17, ntot=7.4e25)
    lay1 = StateLayout(KRF_MODE_SPECIES, n_cells=1)
    rows = build_network("krf", loki_ar_fill.available_buckets())
    cn1 = CompiledNetwork(rows, lay1, loki_ar_fill, mode_optics("krf"))
    cn1.set_photon_energy(cn.hv_eV)
    for c in range(lay.n_cells):
        o1 = cn1.rhs(dens[:, c:c + 1], te[c:c + 1], tg[c:c + 1],
                     ne[c:c + 1], n_ph=3e17, ntot=7.4e25)
        assert np.allclose(o1["ddens"][:, 0], out["ddens"][:, c],
                           rtol=1e-13, atol=0)
        assert np.allclose(o1["du_e"][0], out["du_e"][c], rtol=1e-13)


def test_radiative_flux_analytic(compiled):
    """One radiative row by hand: r_se = n(KrF*(B)) x 1.43e8 s^-1 (the
    Kannari B-state rate; C and Bv radiate into escape, not r_se)."""
    cn, lay = compiled
    dens, te, tg, ne = _random_state(lay, seed=5)
    out = cn.rhs(dens, te, tg, ne, n_ph=0.0, ntot=7.4e25)
    n_upper = dens[lay.index["KrF*(B)"]]
    assert np.allclose(out["r_se"], n_upper * 1.43e8, rtol=1e-12)


def test_stimulated_emission_scales_with_photons(compiled):
    cn, lay = compiled
    dens, te, tg, ne = _random_state(lay, seed=9)
    o1 = cn.rhs(dens, te, tg, ne, n_ph=1e17, ntot=7.4e25)
    o2 = cn.rhs(dens, te, tg, ne, n_ph=2e17, ntot=7.4e25)
    assert np.allclose(o2["r_ste"], 2.0 * o1["r_ste"], rtol=1e-12)
    assert np.allclose(o2["r_x"], 2.0 * o1["r_x"], rtol=1e-12)


def test_three_body_tg_scaling(compiled, loki_ar_fill):
    """The Boichenko Tg^-0.75 set responds to gas temperature; a row
    outside the set does not."""
    cn, lay = compiled
    dens, te, tg, ne = _random_state(lay, seed=13)
    hot = np.full_like(tg, 600.0)
    cold = np.full_like(tg, 300.0)
    o_hot = cn.rhs(dens, te, hot, ne, 0.0, 7.4e25)
    o_cold = cn.rhs(dens, te, cold, ne, 0.0, 7.4e25)
    # AK7 (tg_pow -0.75) feeds Kr2+; its production must fall when hot.
    i = lay.index["Kr2+"]
    assert o_hot["ddens"][i].sum() < o_cold["ddens"][i].sum()


def test_swarm_residual_sink_closes_te_at_the_table_point(loki_arf):
    """0.1.27: on the E/N-direct (circuit) path the electron-energy
    balance must be stationary at the table's own swarm point.

    The bucketed rows book only the level-lumped share of the sigma-set's
    inelastic spend (measured: 12% of e v_d E unbooked at 14.1 Td in the
    Chen fill, rising to 45% at 150 Td) and are E/N-indexed -- Te-FLAT --
    so before 0.1.27 the fluid Te floated x1.67 above the table's own
    <eps> (the 0.1.25 anchor's x1.55 Te row; attributed 2026-08-26). The
    swarm-residual sink is DERIVED from the same table (v_d, <eps>) plus
    the network's own row inventory -- nothing fitted. This gate pins its
    defining property at the pure fill state: dTe/dt = 0 at
    Te = (2/3)<eps>(E/N) exactly (measured 1e-6 eV/s-class against
    +3e9 eV/s-class pre-fix -- teeth margin ~1e9), restoring from both
    sides. Beyond-table channels (stepwise/Penning births) legitimately
    shift the IN-RUN equilibrium below the swarm point; that is dynamics,
    not this identity."""
    from REBUILD.MODULES.network import build_network, mode_optics
    from REBUILD.MODULES.species import (ARF_MODE_SPECIES, StateLayout,
                                         charge_of)
    from REBUILD.MODULES.units import gas_density_m3
    rows = build_network("arf_chen", loki_arf.available_buckets())
    lay = StateLayout(ARF_MODE_SPECIES, n_cells=1)
    cn = CompiledNetwork(rows, lay, loki_arf, mode_optics("arf_chen"))
    cn.set_photon_energy(photon_energy_J(193.0) / EV)
    n_tot = gas_density_m3(3.4, 300.0)

    def dte_dt(te_val, en):
        dens = np.zeros((lay.n_species, 1))
        dens[lay.index["Ne"]] = 0.964 * n_tot
        dens[lay.index["Ar"]] = 0.035 * n_tot
        dens[lay.index["F2"]] = 0.001 * n_tot
        dens[lay.index["Ar+"]] = 1.0e19          # quasineutral seed
        ne = lay.electron_density(dens)
        ch = cn.rhs(dens, np.array([te_val]), np.array([300.0]), ne,
                    0.0, n_tot, en_td=en)
        p_in = (float(ne[0]) * EV * float(loki_arf.drift_velocity(en))
                * en * 1e-21 * n_tot)
        du = float(ch["du_e"][0]) + p_in          # engine adds p_joule
        dne = sum(charge_of(s) * float(ch["ddens"][i, 0])
                  for i, s in enumerate(lay.species))
        return ((2.0 / 3.0) * du / (float(ne[0]) * EV)
                - te_val * dne / float(ne[0]))

    for en in (8.0, 14.1, 40.0):
        te_sw = (2.0 / 3.0) * float(loki_arf.mean_energy_at_en(en))
        assert abs(dte_dt(te_sw, en)) < 1.0       # eV/s; measured 1e-6
        assert dte_dt(0.7 * te_sw, en) > 1e7      # restoring up
        assert dte_dt(1.3 * te_sw, en) < -1e7     # restoring down
