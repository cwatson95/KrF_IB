"""Gate for CompiledNetwork.row_fluxes: the per-row fluxes exposed for the
channel-budget diagnostics must reproduce the RHS exactly — summing each
row's flux times its stoichiometry over all rows equals the rhs ddens when
the non-row terms (stimulated emission, ASE, absorbers) are switched off
(n_ph = 0, ase_rate = 0). A diagnostic that disagrees with the engine it
describes would be worse than none."""

import numpy as np


def _state_mid_run(krf_device, krf_result):
    lay = krf_device.lay
    j = len(krf_result["t_ns"]) // 2            # a beam-on state
    y = krf_result["y_final"].copy()
    # rebuild the mid-run state from the volume-averaged histories (0-D:
    # the averages ARE the state)
    cells = lay.cells_view(y)
    for s, i in lay.index.items():
        cells[i] = krf_result["dens_avg"][s][j]
    ne = lay.electron_density(np.maximum(cells[: lay.n_species], 0.0))
    te = krf_result["Te_eV"][j]
    tg = krf_result["Tg_K"][j]
    return cells, ne, te, tg


def test_row_fluxes_reproduce_rhs(krf_device, krf_result):
    dev = krf_device
    lay = dev.lay
    cells, ne, te, tg = _state_mid_run(dev, krf_result)
    dens = np.maximum(cells[: lay.n_species], 0.0).reshape(lay.n_species, 1)
    ne_a = np.atleast_1d(ne)
    te_a = np.atleast_1d(te)
    tg_a = np.atleast_1d(tg)
    ch = dev.compiled.rhs(dens, te_a, tg_a, ne_a, 0.0, dev.n_tot,
                          en_td=None, ase_rate=0.0)
    fl = dev.compiled.row_fluxes(dens, te_a, tg_a, ne_a, dev.n_tot)
    st = dev.compiled.row_stoichiometry()
    assert set(fl) == set(st) and len(fl) == len(dev.net)
    recon = np.zeros_like(ch["ddens"])
    for rid, flux in fl.items():
        for i_s, c in st[rid]:
            recon[i_s] += c * flux
    scale = np.max(np.abs(ch["ddens"])) + 1e-300
    assert np.max(np.abs(recon - ch["ddens"])) / scale < 1e-10
    # every flux non-negative, and the beam-on state is chemically live
    assert all((f >= 0).all() for f in fl.values())
    assert sum(float(f.sum()) for f in fl.values()) > 0.0
