"""1Dv0_3's transport on a mesh (transport.py) against Fluid1D's.

Needs Python >= 3.11 with numba and Fluid1D, as test_1d.py does."""

import copy
import importlib
import sys

import numpy as np
import pytest

if sys.version_info < (3, 11):
    pytest.skip("Fluid1D needs Python >= 3.11", allow_module_level=True)
pytest.importorskip("numba")


def _m(name):
    return importlib.import_module(f"1Dv0_3.{name}")


bridge = _m("fluid1d")
try:
    bridge.src_dir()
except ImportError as exc:
    pytest.skip(str(exc), allow_module_level=True)

sim1d = _m("sim1d")
meshes = _m("mesh")
tr = _m("transport")
from Fluid1D import Numerics                                   # noqa: E402
from Fluid1D.field import Poisson1DSolver, SemiImplicitField   # noqa: E402
from Fluid1D.plasma_model import DriftDiffusionPlasmaModel     # noqa: E402
from particles.constants import E_CHARGE, EPS0, M_E           # noqa: E402

# Fluid1D's field work: these tests hold the transport to Fluid1D's, and
# the conserving coupling's tests make their own model from the live state
SHORT = dict(pulse_s=20e-9, t_end_s=40e-9, dt_s=25e-12, n_save=11,
             n_profiles=3, verbose=False, field_work="fluid1d")


def _built(tmp_path, steps, **kw):
    """A LightSimulation built and stepped ``steps`` times."""
    p = sim1d.Params1D(**SHORT, output_dir=str(tmp_path), **kw)
    sim = sim1d.simulation_class()(p)
    num = Numerics(steps=1600, dt=p.dt_s, fluid_every=40, fluid_substeps=1,
                   field_heating_cap_factor=None,
                   electron_energy_closure="source_coupled")
    sim.numerics = num
    sim._build(num)
    for i in range(1, steps + 1):
        sim._step(i)
    return sim


def _handle(sim):
    return SemiImplicitField(solver=sim.poisson, phi_left=0.0, phi_right=0.0,
                             cell_centers=sim.cell_centers, prescribed=None,
                             t=0.0, gap_capacitance_f=None, gap_drive_a=0.0)


def _model_args(m):
    return dict(reactions=None, ion_mobility=m.ion_mobility,
                mixture=m.mixture, gas_temp_k=m.gas_temp_k,
                area_m2=m.area_m2, gap_m=m.gap_m,
                secondary_emission_gamma=m.gamma,
                energy_floor_ev=m.eps_floor, energy_ceiling_ev=m.eps_ceiling,
                field_heating_cap_factor=None,
                electron_energy_closure="source_coupled")


def test_uniform_step_is_fluid1d_s_to_the_bit(tmp_path):
    """From a live state (the sheath formed, 300 steps in), one transport
    advance of transport.MeshPlasmaModel on the uniform mesh and one of
    Fluid1D's own model give the same densities, energy, field and
    ledgers to the bit."""
    sim = _built(tmp_path, 300, cells=8, mesh="uniform", name="bit")
    mine = sim.plasma
    theirs = DriftDiffusionPlasmaModel(**_model_args(mine))
    cf_a, cf_b = copy.deepcopy(sim.charged), copy.deepcopy(sim.charged)
    theirs._resolve(sim.fluid, cf_a, sim.coeffs)
    efield = sim._last_efield.copy()
    ha, hb = _handle(sim), _handle(sim)
    ra = theirs.advance(sim.p1.dt_s, efield=efield, neutral_fluid=sim.fluid,
                        charged_fluid=cf_a, coeffs=sim.coeffs, field=ha)
    rb = mine.advance(sim.p1.dt_s, efield=efield, neutral_fluid=sim.fluid,
                      charged_fluid=cf_b, coeffs=sim.coeffs, field=hb)
    assert theirs._last_substeps == mine._last_substeps > 1
    assert np.array_equal(cf_a.density, cf_b.density)
    assert np.array_equal(cf_a.electron_mean_energy_ev,
                          cf_b.electron_mean_energy_ev)
    assert np.array_equal(theirs._last_n_eps, mine._last_n_eps)
    assert np.array_equal(ha.last_efield, hb.last_efield)
    for k in ("field_work_j_m3", "elastic_transfer_j_m3",
              "coefficient_bound_adjustment_j_m3"):
        assert np.array_equal(getattr(ra.energy_closure, k),
                              getattr(rb.energy_closure, k)), k
    assert ra.deposited_power_w == rb.deposited_power_w
    assert ra.conduction_current_a == rb.conduction_current_a


def _poisson_case(mesh, seed=0):
    rng = np.random.default_rng(seed)
    n = mesh.n_cells
    return (rng.normal(0.0, 1e-3, n), rng.uniform(0.0, 1e-3, n),
            np.concatenate([[0.0], rng.normal(0.0, 1e-6, n - 1), [0.0]]),
            25e-12)


def test_field_solve_is_fluid1d_s_on_a_uniform_mesh():
    """MeshPoisson's rows are Fluid1D's own on a uniform mesh (the same
    potential and face fields to the bit)."""
    m = meshes.uniform(40, 1e-2)
    rho, sigma, jf, h = _poisson_case(m)
    ours = tr.MeshPoisson(m)
    theirs = Poisson1DSolver(40, 1e-2 / 40)
    a = theirs.solve_semi_implicit(rho, 0.0, -3.0, sigma_cells=sigma,
                                   j_expl_faces=jf, h=h)
    b = ours.solve_semi_implicit(rho, 0.0, -3.0, sigma_cells=sigma,
                                 j_expl_faces=jf, h=h)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    assert np.array_equal(theirs.last_e_faces, ours.last_e_faces)


def test_field_solve_converges_second_order_on_graded_meshes():
    """A uniform charge between grounded walls, phi = rho x (L - x) / 2 eps0
    exactly: on nested graded meshes the error at the cell centres falls
    fourfold per refinement."""
    rho0 = 1e-4                                       # C/m^3
    length = 1e-2
    errs = []
    m = meshes.graded(48, length, 20e-6)
    for k in range(3):
        mk = meshes.refine(m, k)
        n = mk.n_cells
        phi, _ = tr.MeshPoisson(mk).solve_semi_implicit(
            np.full(n, rho0), 0.0, 0.0, sigma_cells=np.zeros(n),
            j_expl_faces=np.zeros(n + 1), h=0.0)
        x = mk.centres
        exact = rho0 * x * (length - x) / (2 * EPS0)
        errs.append(np.max(np.abs(phi - exact)) / np.max(exact))
    assert errs[0] < 1e-2
    for a, b in zip(errs, errs[1:]):
        assert 3.5 < a / b < 4.5, errs


def test_kernel_field_solve_is_mesh_poisson():
    """The kernel's _solve_field and MeshPoisson assemble the same rows: on
    a graded mesh they agree to rounding."""
    m = meshes.graded(48, 1e-2, 5e-6)
    geo = tr.MeshArrays(m)
    n = m.n_cells
    rng = np.random.default_rng(3)
    dens = np.zeros((n, 2))
    dens[:, 0] = rng.uniform(1e18, 1e20, n)
    dens[:, 1] = dens[:, 0] * rng.uniform(0.98, 1.02, n)
    q = np.array([-1.602176634e-19, 1.602176634e-19])
    mu_e = rng.uniform(0.01, 0.05, n)
    mu_s = np.vstack([mu_e, rng.uniform(1e-4, 2e-4, n)])
    d_face = rng.uniform(0.05, 0.2, (2, n - 1))
    out = [np.empty(k) for k in (n, n + 1, n, n, n + 1, n + 1, n, n, n, n,
                                 n, n, n, n, n)]
    (sigma, jf, rho, a_coef, a_face, e_face, diag, lower, upper, rhs, cp,
     dp, phi, e_axial, eps_eq) = out
    h = 25e-12
    tr._solve_field(dens, np.full(n, 2.45e25), h, geo.widths, geo.spacing,
                    geo.wl, geo.wr, geo.dx_ref, 0.0, -2.0, 0, q, mu_e, mu_s,
                    d_face, False, np.zeros(1), np.zeros(1), 0.03, 60.0,
                    sigma, jf, rho, a_coef, a_face, e_face, diag, lower,
                    upper, rhs, cp, dp, phi, e_axial, eps_eq)
    ref_phi, ref_e = tr.MeshPoisson(m).solve_semi_implicit(
        rho, 0.0, -2.0, sigma_cells=sigma, j_expl_faces=jf, h=h)
    np.testing.assert_allclose(phi, ref_phi, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(e_axial, ref_e, rtol=1e-10, atol=1e-9)


def test_sg_fluxes_hold_the_boltzmann_equilibrium_on_a_graded_mesh():
    """Scharfetter-Gummel is exact for drift against diffusion: an ion in a
    uniform field with n ~ exp(z E x / kT) carries no flux through any face
    of a graded mesh, with its mobility varying cell to cell, as long as
    the Peclet number takes the centre spacing and the drift and the
    diffusion take the same face weights (D = mu kT/e)."""
    m = meshes.graded(48, 1e-2, 10e-6)
    geo = tr.MeshArrays(m)
    n = m.n_cells
    kt = 0.025852                                        # V, kT/e at 300 K
    e0 = -40.0                                           # V/m
    rng = np.random.default_rng(7)
    mu = np.vstack([np.full(n, 0.03), rng.uniform(1e-4, 3e-4, n)])
    d_face = np.empty((2, n - 1))
    tr._face_diffusion(np.full(n, 0.1), mu, 0, kt, geo.wl, geo.wr, d_face)
    dens = np.zeros((n, 2))
    dens[:, 1] = 1e16 * np.exp(e0 * m.centres / kt)     # z = +1
    out = dict(j_left=np.empty(2), j_right=np.empty(2),
               interior=np.empty((2, n - 1)), gam=np.empty((2, n)),
               faces=np.empty(n + 1), rho_before=np.empty(n))
    tr._transport_species(dens, np.full(n, 0.04), np.full(n + 1, e0), 1e-15,
                          geo.widths, geo.spacing, geo.wl, geo.wr, 0,
                          np.array([-1.0, 1.0]),
                          np.array([-1.602176634e-19, 1.602176634e-19]),
                          mu, d_face, np.array([0.0, 400.0]),
                          np.array([1], dtype=np.int64), 0.0,
                          out["j_left"], out["j_right"], out["interior"],
                          out["gam"], out["faces"], out["rho_before"],
                          False, np.ones(n))
    drift = mu[1, :-1] * abs(e0) * 1e16 * np.exp(e0 * m.edges[1:-1] / kt)
    assert np.max(np.abs(out["interior"][1]) / drift) < 1e-11


def test_substep_is_each_cell_s_own_limit():
    """The substep is cfl_safety times the least of every cell's drift,
    thermal-wall and diffusion limits, on its narrowest length: its width
    and the spacings of its two faces."""
    m = meshes.graded(48, 1e-2, 10e-6)
    geo = tr.MeshArrays(m)
    h, d = m.widths, m.spacing
    hmin = h.copy()
    hmin[:-1] = np.minimum(hmin[:-1], d)
    hmin[1:] = np.minimum(hmin[1:], d)
    assert np.array_equal(geo.hmin, hmin)
    n = m.n_cells
    rng = np.random.default_rng(11)
    mu_e = rng.uniform(0.02, 0.05, n)
    mu_s = np.vstack([mu_e, rng.uniform(1e-4, 3e-4, n)])
    # cold electrons and a strong field (the drift binds), and hot ones in
    # a weak field with fast diffusion (the diffusion binds)
    for e_scale, eps_ev, d_scale, binds in ((3e6, 0.03, 0.1, "drift"),
                                            (1e2, 0.03, 3.0, "diffusion")):
        e_axial = rng.normal(0.0, e_scale, n)
        d_e = d_scale * rng.uniform(0.5, 1.0, n)
        d_eps = d_scale * rng.uniform(0.5, 1.5, n)
        eps = np.full(n, eps_ev)
        got = tr._cfl_timestep(e_axial, eps, mu_e, d_e, d_eps, mu_s, 0,
                               geo.hmin, 0.3, 1.0)
        v = np.maximum(5.0 / 3.0 * mu_e * np.abs(e_axial),
                       mu_s[1] * np.abs(e_axial))
        vth = 0.25 * np.sqrt(8 * (2 / 3 * eps_ev * E_CHARGE) / (np.pi * M_E))
        v[[0, -1]] = np.maximum(v[[0, -1]], vth)
        drift = np.min(hmin / v)
        diffusion = np.min(hmin ** 2 / (2 * np.maximum(d_e, d_eps)))
        assert (drift < diffusion) == (binds == "drift")
        assert got == pytest.approx(0.3 * min(drift, diffusion), rel=1e-12,
                                    abs=0), binds


def test_graded_step_conserves_particles(tmp_path):
    """On a graded mesh, one advance from a bump in the gap's middle (no
    charge near the walls, so nothing leaves) keeps every species' count,
    sum n_i h_i, and the electrons' energy, sum n_eps,i h_i (the step moves
    it by its fluxes alone; the field's work goes to the stage), to
    rounding, and refuses the frozen-field path."""
    sim = _built(tmp_path, 1, cells=48, mesh="graded", wall_cell_um=5.0,
                 name="cons")
    m = sim.mesh
    cf = copy.deepcopy(sim.charged)
    x = m.centres / m.length
    # dense enough that Fluid1D's repair of a near-empty cell's mean energy
    # (a blend toward the local field's, weight 1e10 m^-3 / n_e) is nil
    bump = np.exp(-((x - 0.5) / 0.08) ** 2)
    bump[bump < 1e-4] = 0.0
    for s in range(cf.n_species):
        cf.density[:, s] = 1e22 * bump
    cf.electron_mean_energy_ev[:] = 2.0
    before = cf.density.T @ m.widths
    e_before = float(np.sum(cf.density[:, 0] * 2.0 * E_CHARGE * m.widths))
    rep = sim.plasma.advance(1e-12, efield=np.zeros((m.n_cells, 3)),
                             neutral_fluid=sim.fluid, charged_fluid=cf,
                             coeffs=sim.coeffs, field=_handle(sim))
    after = cf.density.T @ m.widths
    np.testing.assert_allclose(after, before, rtol=1e-13, atol=0)
    adjust = rep.energy_closure.coefficient_bound_adjustment_j_m3
    e_after = float(np.sum((sim.plasma._last_n_eps - adjust) * m.widths))
    # the repair's blend, weight 1e10/n_e in the 1e18 m^-3 tail, is 1e-11
    assert e_after == pytest.approx(e_before, rel=1e-9, abs=0)
    assert sim.plasma._last_substeps >= 1
    with pytest.raises(NotImplementedError):
        sim.plasma.advance(1e-12, efield=np.zeros((m.n_cells, 3)),
                           neutral_fluid=sim.fluid, charged_fluid=cf,
                           coeffs=sim.coeffs, field=None)


# ---------------------------------------------------------------------------
# The field's work kept to the field's energy (field_work="conserve")
# ---------------------------------------------------------------------------
def _face_lengths(geo):
    """The faces' lengths: the centres' spacing, half a wall cell at the
    walls (the field's integral over them is the electrodes' difference)."""
    return np.concatenate([[0.5 * geo.widths[0]], geo.spacing,
                           [0.5 * geo.widths[-1]]])


def test_exact_field_is_gauss_s_law_between_grounded_walls():
    """_exact_field holds Gauss's law in every cell, eps0 (E_{i+1/2} -
    E_{i-1/2}) = rho_i w_i, to round-off on a graded mesh, the field's
    integral over the faces' lengths (the grounded electrodes' difference)
    is zero, and _field_energy is sum_f eps0 E_f^2 l_f / 2."""
    m = meshes.graded(48, 1e-2, 10e-6)
    geo = tr.MeshArrays(m)
    n = m.n_cells
    rng = np.random.default_rng(3)
    dens = np.abs(rng.normal(1e16, 3e15, (n, 2)))
    q = np.array([-E_CHARGE, E_CHARGE])
    e = np.empty(n + 1)
    rows = [np.empty(n + 1)] + [np.empty(n) for _ in range(7)]
    tr._exact_field(dens, q, geo.widths, geo.spacing, 0.0, 0.0, *rows, e)
    charge = dens @ q * geo.widths
    np.testing.assert_allclose(EPS0 * np.diff(e), charge, rtol=1e-9,
                               atol=1e-11 * np.max(np.abs(charge)))
    lengths = _face_lengths(geo)
    assert abs(np.sum(e * lengths)) < 1e-12 * np.sum(np.abs(e) * lengths)
    w = tr._field_energy(e, geo.widths, geo.spacing)
    assert w == pytest.approx(0.5 * EPS0 * np.sum(e * e * lengths),
                              rel=1e-13, abs=0)


def test_conserving_work_keeps_the_field_s_energy(tmp_path):
    """From a live state with its sheaths formed (300 steps into a 48-cell
    graded run, wall cells of 10 um), one advance with field_work=
    'conserve': its field's energy at the start and the end are the exact
    fields' of the charge before and after, the faces' work matches their
    difference to round-off (the discrete energy theorem between grounded
    electrodes), the cells' shares, electrons' and ions', sum to the faces'
    total, the ions gain and the electrons lose, and no density goes
    negative, the positive transport having scaled the outflows the wall
    cells could not supply (Fluid1D's floor there made charge and missed
    the theorem by 44 % of the work)."""
    sim = _built(tmp_path, 300, cells=48, mesh="graded", wall_cell_um=10.0,
                 neutrals="still", name="thm")
    cons = tr.MeshPlasmaModel(mesh=sim.mesh, field_work="conserve",
                              **_model_args(sim.plasma))
    cf = copy.deepcopy(sim.charged)
    cons._resolve(sim.fluid, cf, sim.coeffs)
    w_before = cons.field_energy_j_m2(cf)
    rep = cons.advance(sim.p1.dt_s, efield=sim._last_efield.copy(),
                       neutral_fluid=sim.fluid, charged_fluid=cf,
                       coeffs=sim.coeffs, field=_handle(sim))
    w0, w1 = cons._last_field_energy_j_m2
    work = cons._last_field_work_j_m2
    widths = sim.mesh.widths
    fe = rep.energy_closure.field_work_j_m3
    fi = cons._last_ion_work_j_m3
    scale = float(np.sum((np.abs(fe) + np.abs(fi)) * widths))
    assert w0 > 0.0 and scale > 0.0
    assert w0 == pytest.approx(w_before, rel=1e-12, abs=0)
    assert w1 == pytest.approx(cons.field_energy_j_m2(cf), rel=1e-12, abs=0)
    assert abs((w1 - w0) + work) < 1e-8 * scale
    assert abs(float(np.sum((fe + fi) * widths)) - work) < 1e-12 * scale
    assert np.sum(fi * widths) > 0.0 > np.sum(fe * widths)
    assert cons._last_positivity_limits > 0
    assert np.all(cf.density >= 0.0)
    with pytest.raises(NotImplementedError):
        h = _handle(sim)
        h.phi_right = 1.0
        cons.advance(sim.p1.dt_s, efield=sim._last_efield.copy(),
                     neutral_fluid=sim.fluid, charged_fluid=cf,
                     coeffs=sim.coeffs, field=h)


def test_conserving_work_moves_the_electrons_energy_by_their_share(tmp_path):
    """From the same live state (300 steps into the 48-cell graded run),
    one substep with field_work='conserve' and one with Fluid1D's: the
    charges move alike to the bit (no outflow is scaled over so short a
    step), and the electrons' energy before the bounds differs by their
    share of the faces' work, either sign, which Fluid1D's step leaves to
    the stage instead."""
    sim = _built(tmp_path, 300, cells=48, mesh="graded", wall_cell_um=10.0,
                 neutrals="still", name="twin")
    w = sim.mesh.widths
    out = {}
    for mode in ("fluid1d", "conserve"):
        model = tr.MeshPlasmaModel(mesh=sim.mesh, field_work=mode,
                                   **_model_args(sim.plasma))
        cf = copy.deepcopy(sim.charged)
        model._resolve(sim.fluid, cf, sim.coeffs)
        rep = model.advance(1e-13, efield=sim._last_efield.copy(),
                            neutral_fluid=sim.fluid, charged_fluid=cf,
                            coeffs=sim.coeffs, field=_handle(sim))
        assert model._last_substeps == 1, mode
        ec = rep.energy_closure
        out[mode] = (cf.density.copy(),
                     model._last_n_eps - ec.coefficient_bound_adjustment_j_m3,
                     ec.field_work_j_m3.copy(), model)
    dens_f, raw_f, _, _ = out["fluid1d"]
    dens_c, raw_c, fe, cons = out["conserve"]
    assert cons._last_positivity_limits == 0
    np.testing.assert_array_equal(dens_c, dens_f)
    share = float(np.sum(fe * w))
    assert share < 0.0 and np.any(fe > 0.0)
    assert float(np.sum((raw_c - raw_f) * w)) == pytest.approx(share,
                                                               rel=1e-9)

@pytest.mark.parametrize("positive", [True, False])
def test_positive_transport_keeps_the_count_where_the_floor_makes_it(
        positive):
    """An ion drifted out of 1-um cells faster than they hold it (900 m/s
    over 2 ns): with positive the outflows are scaled to the cells'
    content, so nothing goes negative, the last cell gives the wall all it
    holds and no more, and the count changes by the wall fluxes alone,
    exactly. Fluid1D's floor (positive False) sets the negative values to 0
    and so makes particles."""
    m = meshes.uniform(6, 6e-6)
    geo = tr.MeshArrays(m)
    n = m.n_cells
    kt = 0.025852
    mu = np.vstack([np.full(n, 0.03), np.full(n, 3e-4)])
    d_face = np.empty((2, n - 1))
    tr._face_diffusion(np.full(n, 0.1), mu, 0, kt, geo.wl, geo.wr, d_face)
    dens = np.zeros((n, 2))
    dens[:, 1] = [1e16, 1e10, 1e16, 1e16, 1e10, 1e16]
    out = dict(j_left=np.empty(2), j_right=np.empty(2),
               interior=np.empty((2, n - 1)), gam=np.empty((2, n)),
               faces=np.empty(n + 1), rho_before=np.empty(n))
    before = float(dens[:, 1] @ geo.widths)
    held = dens[n - 1, 1] * geo.widths[n - 1]
    h = 2e-9
    _, n_scaled = tr._transport_species(
        dens, np.full(n, 0.04), np.full(n + 1, 3e6), h, geo.widths,
        geo.spacing, geo.wl, geo.wr, 0, np.array([-1.0, 1.0]),
        np.array([-E_CHARGE, E_CHARGE]), mu, d_face, np.array([0.0, 400.0]),
        np.array([1], dtype=np.int64), 0.0, out["j_left"], out["j_right"],
        out["interior"], out["gam"], out["faces"], out["rho_before"],
        positive, np.ones(n))
    after = float(dens[:, 1] @ geo.widths)
    walls = h * (out["j_right"][1] - out["j_left"][1])
    assert np.all(dens >= 0.0)
    if positive:
        assert n_scaled > 0
        assert h * out["j_right"][1] == pytest.approx(held, rel=1e-12, abs=0)
        assert after == pytest.approx(before - walls, rel=1e-14, abs=0)
    else:
        assert n_scaled == 0
        assert after > (before - walls) * (1 + 1e-3)
