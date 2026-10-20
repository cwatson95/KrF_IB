"""1Dv0_3's charged transport on a graded mesh: Fluid1D's step, cell by cell.

Fluid1D (``nathaniel/1DModelKrF``, Nathaniel's, AGPL-3.0) moves the
charged species and the electron energy with one dx everywhere: in every
Scharfetter-Gummel flux, divergence and CFL limit, and in its
semi-implicit Poisson solve (``Fluid1D/transport_kernel.py``). Fluid1D is
imported, never edited. This module is that kernel's coupled substep loop
with the mesh in place of dx, and as a derivative of it it carries the
AGPL. It calls Fluid1D's own geometry-free kernel functions (the mean
energy's recovery, the coefficient lookup, the Bernoulli function, the
gap node and the diagnostic inelastic sink) and replaces those that read
dx:

* a face flux takes the centre-to-centre spacing (its Peclet number and
  D/spacing), and a cell's update its own width;
* a face average (mobility, diffusion, the field's coefficient) weighs
  the two cells by their half-widths, the two half-cells in series to
  first order;
* the field solve's wall faces sit half a wall cell from the first centre;
* the substep limit is each cell's own, from its width and its narrower
  face;
* the gap's mean current and conductivity are length-weighted.

On a uniform mesh each of these is Fluid1D's arithmetic exactly (the
half-width weights are 0.5 and 0.5, and 0.5 a + 0.5 b is 0.5 (a + b) to
the bit), so the step is Fluid1D's to the bit there
(``tests/test_transport.py``).

With ``field_work="conserve"`` (1Dv0_3's default since 2026-10-06) the
field's work is kept to the field's energy. Each substep pairs the face
currents that moved the charge with the field time-centred between the
exact fields of the charge before and after, so between grounded
electrodes the work is minus the field's energy change to round-off (the
discrete energy theorem, above _exact_field). The electrons' share goes
into their energy in the step, either sign, and the ions' is returned for
the gas. The species update then scales any outflow a cell cannot supply
to its content, where Fluid1D's floor sets a negative density to zero and
so makes charge, which would break the theorem. ``"fluid1d"`` is Fluid1D's
step unchanged.
"""

from __future__ import annotations

import warnings

import numpy as np

from .fluid1d import load

load()
from Fluid1D import transport_kernel as _fk                    # noqa: E402
from Fluid1D.field import Poisson1DSolver                      # noqa: E402
from Fluid1D.plasma_model import (DriftDiffusionPlasmaModel,    # noqa: E402
                                  EnergyClosureReport, StepReport)
from particles.constants import E_CHARGE, EPS0, M_E            # noqa: E402

try:
    from numba import njit
except ImportError:                                            # pragma: no cover
    def njit(*args, **kwargs):
        def deco(fn):
            return fn
        return deco(args[0]) if args and callable(args[0]) else deco

_bernoulli = _fk._bernoulli
_recover_energy = _fk._recover_energy
_lookup_transport = _fk._lookup_transport
_advance_gap = _fk._advance_gap
_inelastic_power_cell = _fk._inelastic_power_cell


class MeshArrays:
    """The mesh as the kernel reads it (contiguous float64 arrays)."""

    def __init__(self, mesh):
        self.mesh = mesh
        self.widths = np.ascontiguousarray(mesh.widths)
        self.spacing = np.ascontiguousarray(mesh.spacing)
        wl, wr = mesh.weights
        self.wl, self.wr = np.ascontiguousarray(wl), np.ascontiguousarray(wr)
        n = self.widths.size
        # a cell's narrowest length: its width and its faces' spacings (the
        # wall's flux is the thermal one, limited by the cell's width)
        hmin = self.widths.copy()
        hmin[:-1] = np.minimum(hmin[:-1], self.spacing)
        hmin[1:] = np.minimum(hmin[1:], self.spacing)
        self.hmin = hmin
        self.length = float(mesh.length)
        # the field solve's scale: Fluid1D multiplies its rows by dx; the mean
        # width is that dx on a uniform mesh, and any length on another
        self.dx_ref = float(self.widths[0]) if mesh.uniform else \
            self.length / n
        self.weighted = not mesh.uniform


# ---------------------------------------------------------------------------
# The kernel
# ---------------------------------------------------------------------------
@njit(cache=True)
def _face_diffusion(d_e, mu_s, e_col, kt_gas_ev, wl, wr, d_face_s):
    """Face diffusion by half-width weights (Fluid1D's 0.5/0.5 on a
    uniform mesh, where it equals _lookup_transport's to the bit)."""
    n = d_e.shape[0]
    for s in range(mu_s.shape[0]):
        if s == e_col:
            for i in range(n - 1):
                d_face_s[s, i] = wl[i] * d_e[i] + wr[i] * d_e[i + 1]
        else:
            for i in range(n - 1):
                d_face_s[s, i] = kt_gas_ev * (wl[i] * mu_s[s, i]
                                              + wr[i] * mu_s[s, i + 1])


@njit(cache=True)
def _cfl_timestep(e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col, hmin,
                  cfl_safety, remaining):
    """Each cell's drift, thermal-wall and diffusion limit, and the least
    of them (Fluid1D's dx / v_max and dx^2 / 2 D_max on a uniform mesh)."""
    n = eps.shape[0]
    dt_cfl = 1.0e300
    for i in range(n):
        v = (5.0 / 3.0) * mu_e[i] * abs(e_axial[i])
        for s in range(mu_s.shape[0]):
            if s == e_col:
                continue
            value = mu_s[s, i] * abs(e_axial[i])
            if value > v:
                v = value
        if i == 0 or i == n - 1:
            vth = 0.25 * np.sqrt(
                8.0 * ((2.0 / 3.0) * max(eps[i], 0.0) * E_CHARGE)
                / (np.pi * M_E))
            if vth > v:
                v = vth
        d = d_e[i] if d_e[i] > d_eps[i] else d_eps[i]
        if v > 0.0:
            t = hmin[i] / v
            if t < dt_cfl:
                dt_cfl = t
        if d > 0.0:
            t = hmin[i] * hmin[i] / (2.0 * d)
            if t < dt_cfl:
                dt_cfl = t
    h = cfl_safety * dt_cfl
    if not np.isfinite(h) or h <= 0.0:
        h = remaining
    return min(h, remaining)


@njit(cache=True)
def _solve_field(density, n_gas_safe, h, widths, spacing, wl, wr, dx_ref,
                 phi_l, phi_r, e_col, q, mu_e, mu_s, d_face_s, has_eq, eq_en,
                 eq_eps, eps_lo, eps_hi, sigma_cell, j_expl_f, rho, a_coef,
                 a_face, e_face, diag, lower, upper, rhs, cp, dp, phi_cells,
                 e_axial, eps_eq):
    """The cell-centred semi-implicit Poisson solve on the mesh: each row
    is -[A (phi_{i+1} - phi_i)/d_{i+1/2} - A (phi_i - phi_{i-1})/d_{i-1/2}]
    = h_i rho_i - h (J_{i+1/2} - J_{i-1/2}), times dx_ref, the wall faces
    half a cell from the first centres (Fluid1D's rows on a uniform mesh)."""
    n = density.shape[0]
    s_tot = density.shape[1]
    for i in range(n):
        conductivity = E_CHARGE * density[i, e_col] * mu_e[i]
        charge = 0.0
        for s in range(s_tot):
            if s != e_col:
                conductivity += E_CHARGE * density[i, s] * mu_s[s, i]
            charge += q[s] * density[i, s]
        sigma_cell[i] = conductivity
        rho[i] = charge
    for face in range(n + 1):
        j_expl_f[face] = 0.0
    for s in range(s_tot):
        for i in range(n - 1):
            j_expl_f[i + 1] += q[s] * (
                -d_face_s[s, i] * (density[i + 1, s] - density[i, s])
                / spacing[i])
    for i in range(n):
        a_coef[i] = EPS0 + h * max(sigma_cell[i], 0.0)
    # face conductances, scaled by dx_ref (a_face itself on a uniform mesh,
    # and 2 a_face at the walls)
    a_face[0] = a_coef[0] * (dx_ref / (0.5 * widths[0]))
    a_face[n] = a_coef[n - 1] * (dx_ref / (0.5 * widths[n - 1]))
    for face in range(1, n):
        a_face[face] = (wl[face - 1] * a_coef[face - 1]
                        + wr[face - 1] * a_coef[face]) \
            * (dx_ref / spacing[face - 1])
    for i in range(n):
        diag[i] = a_face[i] + a_face[i + 1]
        lower[i] = -a_face[i]
        upper[i] = -a_face[i + 1]
        rhs[i] = dx_ref * widths[i] * rho[i] \
            - h * dx_ref * (j_expl_f[i + 1] - j_expl_f[i])
    rhs[0] += a_face[0] * phi_l
    rhs[n - 1] += a_face[n] * phi_r
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denominator = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denominator
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denominator
    phi_cells[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        phi_cells[i] = dp[i] - cp[i] * phi_cells[i + 1]
    e_face[0] = -(phi_cells[0] - phi_l) / (0.5 * widths[0])
    e_face[n] = -(phi_r - phi_cells[n - 1]) / (0.5 * widths[n - 1])
    for face in range(1, n):
        e_face[face] = -(phi_cells[face] - phi_cells[face - 1]) \
            / spacing[face - 1]
    for i in range(n):
        e_axial[i] = 0.5 * (e_face[i] + e_face[i + 1])
        if has_eq:
            en_td = abs(e_axial[i]) / n_gas_safe[i] / 1.0e-21
            value = np.interp(en_td, eq_en, eq_eps)
            eps_eq[i] = min(max(value, eps_lo), eps_hi)


# ---------------------------------------------------------------------------
# The field's work, kept to the field's energy (field_work="conserve")
# ---------------------------------------------------------------------------
# Gauss's law per cell, eps0 (E_{i+1/2} - E_{i-1/2}) = rho_i w_i, and the
# charge's continuity, rho_i' = rho_i - h (J_{i+1/2} - J_{i-1/2}) / w_i, give
# eps0 dE_f + h J_f the same on every face. Over the faces' lengths l_f (the
# centres' spacing, half a wall cell at the walls) sum_f E_f l_f is the
# electrodes' potential difference, zero when both are grounded, so the
# field's energy W = sum_f eps0 E_f^2 l_f / 2 changes by exactly
#     W' - W = -h sum_f J_f (E_f + E_f') / 2 l_f,
# the work of the face currents that moved the charge, in the field
# time-centred between the exact fields of the charge before and after. The
# step's own field, solved semi-implicitly from a predicted charge, drives
# the fluxes. The work is booked against these two exact fields.

@njit(cache=True)
def _exact_field(density, q, widths, spacing, phi_l, phi_r, a_face, diag,
                 lower, upper, rhs, cp, dp, phi_cells, e_out):
    """The field [V/m] at the n + 1 faces of the charge as it stands: the
    rows of _solve_field with h = 0 (Gauss's law per cell, the wall faces
    half a cell from the first centres)."""
    n = density.shape[0]
    s_tot = density.shape[1]
    a_face[0] = EPS0 / (0.5 * widths[0])
    a_face[n] = EPS0 / (0.5 * widths[n - 1])
    for face in range(1, n):
        a_face[face] = EPS0 / spacing[face - 1]
    for i in range(n):
        charge = 0.0
        for s in range(s_tot):
            charge += q[s] * density[i, s]
        diag[i] = a_face[i] + a_face[i + 1]
        lower[i] = -a_face[i]
        upper[i] = -a_face[i + 1]
        rhs[i] = widths[i] * charge
    rhs[0] += a_face[0] * phi_l
    rhs[n - 1] += a_face[n] * phi_r
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denominator = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denominator
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denominator
    phi_cells[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        phi_cells[i] = dp[i] - cp[i] * phi_cells[i + 1]
    e_out[0] = -(phi_cells[0] - phi_l) / (0.5 * widths[0])
    e_out[n] = -(phi_r - phi_cells[n - 1]) / (0.5 * widths[n - 1])
    for face in range(1, n):
        e_out[face] = -(phi_cells[face] - phi_cells[face - 1]) \
            / spacing[face - 1]


@njit(cache=True)
def _field_energy(e_faces, widths, spacing):
    """The field's energy per unit area [J/m^2], sum_f eps0 E_f^2 l_f / 2."""
    n = widths.shape[0]
    w = 0.25 * EPS0 * (e_faces[0] * e_faces[0] * widths[0]
                       + e_faces[n] * e_faces[n] * widths[n - 1])
    for face in range(1, n):
        w += 0.5 * EPS0 * e_faces[face] * e_faces[face] * spacing[face - 1]
    return w


@njit(cache=True)
def _face_work(e_prev, e_new, j_left, j_right, interior, q, e_col, widths,
               spacing, pe_cell, pi_cell):
    """The field's work over a substep [W/m^3]: each face's current times
    its time-centred field, J_f (E_f + E_f') / 2, the electrons' and the
    ions' apart, each face's split between its two cells by the half-cell
    each holds of the face's length (all of a wall face to its cell).
    Returns the gap's total, sum_f J_f Ebar_f l_f [W/m^2], which the cells'
    shares times their widths sum to."""
    n = widths.shape[0]
    s_tot = q.shape[0]
    for i in range(n):
        pe_cell[i] = 0.0
        pi_cell[i] = 0.0
    total = 0.0
    for face in range(n + 1):
        e_bar = 0.5 * (e_prev[face] + e_new[face])
        j_e = 0.0
        j_i = 0.0
        for s in range(s_tot):
            if face == 0:
                flux = j_left[s]
            elif face == n:
                flux = j_right[s]
            else:
                flux = interior[s, face - 1]
            if s == e_col:
                j_e += q[s] * flux
            else:
                j_i += q[s] * flux
        p_e = j_e * e_bar
        p_i = j_i * e_bar
        if face > 0:
            pe_cell[face - 1] += 0.5 * p_e
            pi_cell[face - 1] += 0.5 * p_i
        if face < n:
            pe_cell[face] += 0.5 * p_e
            pi_cell[face] += 0.5 * p_i
        if face == 0:
            length = 0.5 * widths[0]
        elif face == n:
            length = 0.5 * widths[n - 1]
        else:
            length = spacing[face - 1]
        total += (p_e + p_i) * length
    return total


@njit(cache=True)
def _transport_species(density, eps, e_face, h, widths, spacing, wl, wr,
                       e_col, z, q, mu_s, d_face_s, ion_vth, pos_ion,
                       gamma_se, j_left, j_right, interior, gam, faces,
                       rho_before, positive, scale_cell):
    """SG and electrode fluxes, and the conservative update of every
    charged density, on the mesh. With positive, the fluxes a cell cannot
    supply are scaled to its content (Fluid1D floors the result instead).
    Returns the largest relative Ampere residual and the scaled cells."""
    n = density.shape[0]
    s_tot = density.shape[1]
    for i in range(n):
        charge = 0.0
        for s in range(s_tot):
            charge += q[s] * density[i, s]
        rho_before[i] = charge
    for s in range(s_tot):
        for i in range(n - 1):
            velocity = z[s] * (wl[i] * mu_s[s, i] + wr[i] * mu_s[s, i + 1]) \
                * e_face[i + 1]
            diffusion = d_face_s[s, i]
            if abs(diffusion) < 1.0e-30:
                diffusion = 1.0e-30
            peclet = velocity * spacing[i] / diffusion
            interior[s, i] = (diffusion / spacing[i]) * (
                _bernoulli(-peclet) * density[i, s]
                - _bernoulli(peclet) * density[i + 1, s])
        if s == e_col:
            vth0 = np.sqrt(8.0 * ((2.0 / 3.0) * max(eps[0], 0.0) * E_CHARGE)
                           / (np.pi * M_E))
            vthl = np.sqrt(8.0 * ((2.0 / 3.0) * max(eps[-1], 0.0) * E_CHARGE)
                           / (np.pi * M_E))
        else:
            vth0 = ion_vth[s]
            vthl = ion_vth[s]
        v0 = z[s] * mu_s[s, 0] * e_face[0]
        vl = z[s] * mu_s[s, n - 1] * e_face[n]
        j_left[s] = (-0.25 * density[0, s] * vth0
                     + (v0 if v0 < 0.0 else 0.0) * density[0, s])
        j_right[s] = (0.25 * density[n - 1, s] * vthl
                      + (vl if vl > 0.0 else 0.0) * density[n - 1, s])
    ion_cathode_flux = 0.0
    for index in range(pos_ion.shape[0]):
        value = j_right[pos_ion[index]]
        if value > 0.0:
            ion_cathode_flux += value
    j_right[e_col] -= gamma_se * ion_cathode_flux

    n_scaled = 0
    for s in range(s_tot):
        faces[0] = j_left[s]
        faces[n] = j_right[s]
        for i in range(n - 1):
            faces[i + 1] = interior[s, i]
        if positive:
            # what a cell's faces would carry out over h beyond its content
            # is scaled to the content, each face by its upwind cell's
            # factor, so the count stays non-negative and its continuity
            # exact (Fluid1D's floor instead sets a negative value to 0,
            # which makes particles and charge)
            for i in range(n):
                out = 0.0
                if faces[i + 1] > 0.0:
                    out += faces[i + 1]
                if faces[i] < 0.0:
                    out -= faces[i]
                out *= h / widths[i]
                if out > density[i, s] and out > 0.0:
                    scale_cell[i] = max(density[i, s], 0.0) / out
                    n_scaled += 1
                else:
                    scale_cell[i] = 1.0
            for face in range(n + 1):
                flux = faces[face]
                if flux > 0.0 and face > 0:
                    faces[face] = flux * scale_cell[face - 1]
                elif flux < 0.0 and face < n:
                    faces[face] = flux * scale_cell[face]
            j_left[s] = faces[0]
            j_right[s] = faces[n]
            for i in range(n - 1):
                interior[s, i] = faces[i + 1]
        for i in range(n):
            gam[s, i] = 0.5 * (faces[i] + faces[i + 1])
            value = density[i, s] - h * (faces[i + 1] - faces[i]) / widths[i]
            density[i, s] = value if value > 0.0 else 0.0

    maximum_residual = 0.0
    for i in range(n):
        rho_after = 0.0
        j_face_l = 0.0
        j_face_r = 0.0
        for s in range(s_tot):
            rho_after += q[s] * density[i, s]
            flux_left = j_left[s] if i == 0 else interior[s, i - 1]
            flux_right = j_right[s] if i == n - 1 else interior[s, i]
            j_face_l += q[s] * flux_left
            j_face_r += q[s] * flux_right
        charge_rate = (rho_after - rho_before[i]) / h
        current_div = (j_face_r - j_face_l) / widths[i]
        ampere_abs = abs(charge_rate + current_div)
        ampere_scale = max(abs(charge_rate), abs(current_div), 1.0e-300)
        relative = ampere_abs / ampere_scale
        if relative > maximum_residual:
            maximum_residual = relative
    return maximum_residual, n_scaled


@njit(cache=True)
def _transport_energy(density, neutral, n_eps, eps, e_face, e_axial, mu_eps,
                      d_eps, mu_e, gam, j_left, j_right, e_col, widths,
                      spacing, wl, wr, h, area, substep_current, v_b_supply,
                      heating_cap_factor, m_gas, eps_gas_ev, closure_mode,
                      key_grid, te_arr, loss_ev, base_loss_k, loss_rate_kind,
                      loss_rate_lut, loss_analytic, loss_react_kind,
                      loss_react_col, loss_react_order, loss_react_count,
                      eps_floor, eps_ceiling, n_e_floor, faces, pf_arr,
                      field_work, elastic_transfer, continuous_sink,
                      bound_adjustment):
    """The electron energy PDE on the mesh, and this substep's ledgers.
    The field's work on a cell's electrons is Fluid1D's, -e Gamma E of the
    cell's mean flux and mean field, booked for the split to deliver
    (closure_mode 0) or applied with the collisions (1). In closure_mode 2
    (field_work="conserve") it is the caller's pf_arr, the electrons' share
    of the faces' work (_face_work), applied here in full, either sign."""
    n = density.shape[0]
    faces[0] = eps[0] * E_CHARGE * j_left[e_col]
    faces[n] = eps[n - 1] * E_CHARGE * j_right[e_col]
    for i in range(n - 1):
        velocity = -(wl[i] * mu_eps[i] + wr[i] * mu_eps[i + 1]) * e_face[i + 1]
        diffusion = wl[i] * d_eps[i] + wr[i] * d_eps[i + 1]
        if abs(diffusion) < 1.0e-30:
            diffusion = 1.0e-30
        peclet = velocity * spacing[i] / diffusion
        faces[i + 1] = (diffusion / spacing[i]) * (
            _bernoulli(-peclet) * n_eps[i] - _bernoulli(peclet) * n_eps[i + 1])

    limiter_activity = 0
    if closure_mode != 2:
        positive_power = 0.0
        negative_power = 0.0
        for i in range(n):
            pf_arr[i] = -E_CHARGE * gam[e_col, i] * e_axial[i]
            if pf_arr[i] > 0.0:
                positive_power += pf_arr[i] * widths[i]
            else:
                negative_power -= pf_arr[i] * widths[i]
        positive_power *= area
        negative_power *= area
        cap_w = heating_cap_factor * abs(area * substep_current * v_b_supply) \
            + 1.0e3
        magnitude = positive_power if positive_power > negative_power \
            else negative_power
        if magnitude > cap_w:
            scale = cap_w / magnitude
            for i in range(n):
                pf_arr[i] *= scale
            limiter_activity += 1

    deposited_dx_j = 0.0
    maximum_jump = 0.0
    maximum_residual = 0.0
    for i in range(n):
        transport_rate = -(faces[i + 1] - faces[i]) / widths[i]
        p_field = pf_arr[i]
        collision_rate = E_CHARGE / (M_E * mu_e[i]) if mu_e[i] > 0.0 else 0.0
        p_elastic = ((2.0 * M_E / m_gas[i]) * collision_rate
                     * density[i, e_col] * (eps[i] - eps_gas_ev) * E_CHARGE)
        p_inelastic = 0.0
        if closure_mode == 1:
            p_inelastic = _inelastic_power_cell(
                i, eps, density, neutral, key_grid, te_arr, loss_ev,
                base_loss_k, loss_rate_kind, loss_rate_lut, loss_analytic,
                loss_react_kind, loss_react_col, loss_react_order,
                loss_react_count)
        old_energy = n_eps[i]
        if closure_mode == 1:
            local_power = p_field - p_elastic - p_inelastic
        elif closure_mode == 2:
            local_power = p_field
        else:
            local_power = 0.0
        raw_energy = old_energy + h * (transport_rate + local_power)
        lower_bound = eps_floor * density[i, e_col] * E_CHARGE
        upper_bound = (eps_ceiling * max(density[i, e_col], 0.0) * E_CHARGE
                       + 1.0e-300)
        next_energy = raw_energy
        if raw_energy < lower_bound:
            next_energy = lower_bound
            limiter_activity += 1
        elif raw_energy > upper_bound:
            next_energy = upper_bound
            limiter_activity += 1
        n_eps[i] = next_energy
        adjustment = next_energy - raw_energy
        field_work[i] += p_field * h
        elastic_transfer[i] += p_elastic * h
        continuous_sink[i] += p_inelastic * h
        bound_adjustment[i] += adjustment
        density_safe = max(density[i, e_col], n_e_floor)
        jump = abs(next_energy - old_energy) / (density_safe * E_CHARGE)
        if jump > maximum_jump:
            maximum_jump = jump
        residual = (next_energy - old_energy
                    - h * (transport_rate + local_power) - adjustment)
        if abs(residual) > maximum_residual:
            maximum_residual = abs(residual)
        deposited_dx_j += p_field * widths[i] * h
    return deposited_dx_j, maximum_jump, limiter_activity, maximum_residual


@njit(cache=True)
def _coupled_loop(density, n_eps, e_axial, neutral, n_gas_safe, m_gas, dt,
                  widths, spacing, wl, wr, hmin, dx_ref, length, weighted,
                  phi_l, phi_r, c_g, i_drive, area, gap, cfl_safety,
                  max_substeps, eps_gas_ev, kt_gas_ev, gamma_se,
                  heating_cap_factor, eps_lo, eps_hi, lmea_ref, n_e_floor,
                  eps_floor, eps_ceiling, e_col, z, q, ion_mu_n, ion_vth,
                  pos_ion, key_grid, mob_arr, dif_arr, emob_arr, edif_arr,
                  has_eq, eq_en, eq_eps, closure_mode, loss_ev, base_loss_k,
                  loss_rate_kind, loss_rate_lut, loss_analytic,
                  loss_react_kind, loss_react_col, loss_react_order,
                  loss_react_count, te_arr, conserve, ion_work):
    """Fluid1D's _coupled_loop on the mesh, in the same sequence. With
    conserve (field_work="conserve", grounded electrodes) the field's work
    is the faces' between the exact fields of the charge before and after
    each substep: the electrons' share goes into their energy here, the
    ions' is summed into ion_work [J/m^3], and the field's energy at the
    start and the end [J/m^2] come back with the total work, which their
    difference matches to round-off."""
    n = density.shape[0]
    s_tot = density.shape[1]

    eps = np.empty(n)
    eps_eq = np.empty(n)
    mu_e = np.empty(n)
    d_e = np.empty(n)
    mu_eps = np.empty(n)
    d_eps = np.empty(n)
    sigma_cell = np.empty(n)
    j_expl_f = np.empty(n + 1)
    rho = np.empty(n)
    a_coef = np.empty(n)
    a_face = np.empty(n + 1)
    e_face = np.empty(n + 1)
    diag = np.empty(n)
    lower = np.empty(n)
    upper = np.empty(n)
    rhs = np.empty(n)
    cp = np.empty(n)
    dp = np.empty(n)
    phi_cells = np.empty(n)
    faces = np.empty(n + 1)
    mu_s = np.empty((s_tot, n))
    d_face_s = np.empty((s_tot, n - 1))
    j_left = np.empty(s_tot)
    j_right = np.empty(s_tot)
    interior = np.empty((s_tot, n - 1))
    gam = np.empty((s_tot, n))
    pf_arr = np.empty(n)
    rho_before = np.empty(n)
    field_work = np.zeros(n)
    elastic_transfer = np.zeros(n)
    continuous_sink = np.zeros(n)
    bound_adjustment = np.zeros(n)
    maximum_energy_jump_ev = 0.0
    limiter_activity = 0
    max_ampere_residual = 0.0
    max_energy_residual = 0.0
    # the conserving work's exact fields and their solve's rows
    e_prev = np.zeros(n + 1)
    e_new = np.zeros(n + 1)
    pi_cell = np.zeros(n)
    xa = np.empty(n + 1)
    xdiag = np.empty(n)
    xlower = np.empty(n)
    xupper = np.empty(n)
    xrhs = np.empty(n)
    xcp = np.empty(n)
    xdp = np.empty(n)
    xphi = np.empty(n)
    w_start = 0.0
    w_end = 0.0
    work_total = 0.0
    scale_cell = np.ones(n)
    positivity_limits = 0
    if conserve:
        _exact_field(density, q, widths, spacing, phi_l, phi_r, xa, xdiag,
                     xlower, xupper, xrhs, xcp, xdp, xphi, e_prev)
        w_start = _field_energy(e_prev, widths, spacing)

    for i in range(n):
        if has_eq:
            en_td = abs(e_axial[i]) / n_gas_safe[i] / 1.0e-21
            v = np.interp(en_td, eq_en, eq_eps)
            eps_eq[i] = min(max(v, eps_lo), eps_hi)
        else:
            eps_eq[i] = 0.0

    deposited_dx_j = 0.0
    jcond_mean = 0.0
    gap_e_abs = 0.0
    gap_v2 = 0.0
    remaining = dt
    steps_done = 0
    capped = False
    h_next = -1.0
    v_b = phi_r

    while True:
        _recover_energy(density, n_eps, eps, eps_eq, e_col, n_e_floor,
                        has_eq, lmea_ref, eps_lo, eps_hi)
        _lookup_transport(eps, n_gas_safe, key_grid, mob_arr, dif_arr,
                          emob_arr, edif_arr, e_col, ion_mu_n, kt_gas_ev,
                          mu_e, d_e, mu_eps, d_eps, mu_s, d_face_s)
        _face_diffusion(d_e, mu_s, e_col, kt_gas_ev, wl, wr, d_face_s)

        if h_next < 0.0:
            h = _cfl_timestep(e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col,
                              hmin, cfl_safety, remaining)
        else:
            h = h_next
        if steps_done >= max_substeps:
            capped = True
            h = remaining

        _solve_field(density, n_gas_safe, h, widths, spacing, wl, wr, dx_ref,
                     phi_l, v_b, e_col, q, mu_e, mu_s, d_face_s, has_eq,
                     eq_en, eq_eps, eps_lo, eps_hi, sigma_cell, j_expl_f, rho,
                     a_coef, a_face, e_face, diag, lower, upper, rhs, cp, dp,
                     phi_cells, e_axial, eps_eq)

        ampere_residual, n_scaled = _transport_species(
            density, eps, e_face, h, widths, spacing, wl, wr, e_col, z, q,
            mu_s, d_face_s, ion_vth, pos_ion, gamma_se, j_left, j_right,
            interior, gam, faces, rho_before, conserve, scale_cell)
        positivity_limits += n_scaled
        if ampere_residual > max_ampere_residual:
            max_ampere_residual = ampere_residual

        jc = 0.0
        sg = 0.0
        for i in range(n):
            jci = 0.0
            for s in range(s_tot):
                jci += q[s] * gam[s, i]
            if weighted:
                jc += jci * widths[i]
                sg += sigma_cell[i] * widths[i]
            else:
                jc += jci
                sg += sigma_cell[i]
        if weighted:
            substep_current = jc / length
            substep_sigma = sg / length
        else:
            substep_current = jc / n
            substep_sigma = sg / n
        jcond_mean += substep_current * (h / dt)
        v_b, v_b_supply, gap_energy = _advance_gap(
            v_b, c_g, i_drive, area, gap, h, substep_current, substep_sigma)
        gap_v2 += v_b_supply * v_b_supply * h
        gap_e_abs += gap_energy

        mode = closure_mode
        if conserve:
            # the charge as it now stands, its exact field, and the faces'
            # work between that and the field the substep began with
            _exact_field(density, q, widths, spacing, phi_l, phi_r, xa,
                         xdiag, xlower, xupper, xrhs, xcp, xdp, xphi, e_new)
            work_total += h * _face_work(e_prev, e_new, j_left, j_right,
                                         interior, q, e_col, widths,
                                         spacing, pf_arr, pi_cell)
            for i in range(n):
                ion_work[i] += pi_cell[i] * h
            mode = 2

        energy_deposited, energy_jump, energy_limiters, energy_residual = \
            _transport_energy(
                density, neutral, n_eps, eps, e_face, e_axial, mu_eps, d_eps,
                mu_e, gam, j_left, j_right, e_col, widths, spacing, wl, wr,
                h, area, substep_current, v_b_supply, heating_cap_factor,
                m_gas, eps_gas_ev, mode, key_grid, te_arr, loss_ev,
                base_loss_k, loss_rate_kind, loss_rate_lut, loss_analytic,
                loss_react_kind, loss_react_col, loss_react_order,
                loss_react_count, eps_floor, eps_ceiling, n_e_floor, faces,
                pf_arr, field_work, elastic_transfer, continuous_sink,
                bound_adjustment)
        deposited_dx_j += energy_deposited
        if energy_jump > maximum_energy_jump_ev:
            maximum_energy_jump_ev = energy_jump
        limiter_activity += energy_limiters
        if energy_residual > max_energy_residual:
            max_energy_residual = energy_residual
        if conserve:
            for face in range(n + 1):
                e_prev[face] = e_new[face]

        steps_done += 1
        remaining -= h
        if remaining <= 1.0e-12 * dt or capped:
            break

        h_next = _cfl_timestep(e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col,
                               hmin, cfl_safety, remaining)

    if conserve:
        w_end = _field_energy(e_prev, widths, spacing)
    return (jcond_mean, deposited_dx_j, steps_done, capped, eps_eq, v_b,
            gap_e_abs, gap_v2, field_work, elastic_transfer, continuous_sink,
            bound_adjustment, maximum_energy_jump_ev, limiter_activity,
            max_ampere_residual, max_energy_residual, w_start, w_end,
            work_total, positivity_limits)


def coupled_advance(model, geo: MeshArrays, dt, e_axial0, neutral_fluid,
                    charged_fluid, packed, *, phi_l, phi_r, c_g=None,
                    i_drive=0.0):
    """One coupled advance on the mesh (Fluid1D's coupled_advance, with the
    mesh's arrays); mutates the charged fluid in place."""
    cf = charged_fluid
    n_gas = neutral_fluid.density.sum(axis=1)
    n_gas_safe = np.where(n_gas > 0.0, n_gas, 1.0)
    masses = model._neutral_masses
    m_gas = np.where(n_gas > 0.0, (neutral_fluid.density @ masses)
                     / n_gas_safe, masses.mean())
    e_col = model._e_col
    n_e0 = cf.density[:, e_col]
    n_eps = np.ascontiguousarray(cf.electron_mean_energy_ev * n_e0 * E_CHARGE)
    e_axial = np.ascontiguousarray(e_axial0, dtype=np.float64).copy()
    density = np.ascontiguousarray(cf.density)
    p = packed
    conserve = bool(getattr(model, "conserve_field_work", False))
    ion_work = np.zeros(density.shape[0])
    (jcond_mean, dep_dx, steps, capped, eps_eq, v_b_final, gap_e_abs, gap_v2,
     field_work, elastic_transfer, continuous_sink, bound_adjustment,
     max_jump, limiter_activity, ampere_residual, energy_residual,
     w_start, w_end, work_total, positivity_limits) = \
        _coupled_loop(
            density, n_eps, e_axial, np.ascontiguousarray(neutral_fluid.density),
            np.ascontiguousarray(n_gas_safe), np.ascontiguousarray(m_gas), dt,
            geo.widths, geo.spacing, geo.wl, geo.wr, geo.hmin, geo.dx_ref,
            geo.length, geo.weighted, phi_l, phi_r,
            np.inf if c_g is None else float(c_g), float(i_drive),
            model.area_m2, model.gap_m, model.cfl_safety, model.max_substeps,
            p["eps_gas_ev"], p["kt_gas_ev"], p["gamma_se"],
            model.field_heating_cap_factor, p["eps_lo"], p["eps_hi"],
            p["lmea_ref"], p["n_e_floor"], p["eps_floor"], p["eps_ceiling"],
            p["e_col"], p["z"], p["q"], p["ion_mu_n"], p["ion_vth"],
            p["pos_ion"], p["key_grid"], p["mob_arr"], p["dif_arr"],
            p["emob_arr"], p["edif_arr"], p["has_eq"], p["eq_en"],
            p["eq_eps"],
            1 if model.electron_energy_closure == "continuous_reconciled"
            else 0,
            p["loss_ev"], p["base_loss_k"], p["loss_rate_kind"],
            p["loss_rate_lut"], p["loss_analytic"], p["loss_react_kind"],
            p["loss_react_col"], p["loss_react_order"], p["loss_react_count"],
            p["te_arr"], conserve, ion_work)
    cf.density[:] = density
    # the conserving coupling's books: the ions' share of the field's work
    # [J/m^3] (for the gas), the field's energy at the advance's start and
    # end and the faces' total work [J/m^2], whose sum is the step's error
    model._last_ion_work_j_m3 = ion_work
    model._last_field_energy_j_m2 = (float(w_start), float(w_end))
    model._last_field_work_j_m2 = float(work_total)
    model._last_positivity_limits = int(positivity_limits)
    return _fk._CoupledAdvanceResult(
        jcond_mean_a_m2=float(jcond_mean),
        deposited_j=float(dep_dx * cf.area_m2),
        n_eps_j_m3=n_eps, e_axial_v_m=e_axial, steps=int(steps),
        capped=bool(capped), equilibrium_energy_ev=eps_eq,
        v_b_final_v=float(v_b_final),
        gap_energy_absorbed_j=float(gap_e_abs),
        gap_voltage_squared_dt_v2_s=float(gap_v2),
        field_work_j_m3=field_work, elastic_transfer_j_m3=elastic_transfer,
        continuous_sink_j_m3=continuous_sink,
        bound_adjustment_j_m3=bound_adjustment,
        maximum_energy_jump_ev=float(max_jump),
        limiter_activity=int(limiter_activity),
        max_ampere_residual=float(ampere_residual),
        max_energy_residual_j_m3=float(energy_residual))


# ---------------------------------------------------------------------------
# The model and the field solver
# ---------------------------------------------------------------------------
class MeshPlasmaModel(DriftDiffusionPlasmaModel):
    """Fluid1D's drift-diffusion model with the coupled step on a mesh
    (``coupled_advance`` above). Only the coupled path exists: a frozen
    field or a prescribed drive is Fluid1D's uniform NumPy path, which
    reads one dx, and is refused here.

    field_work: 'fluid1d' books Fluid1D's -e Gamma E of each cell's mean
    flux and mean field for the split to deliver; 'conserve' applies the
    faces' work between the exact fields of the charge before and after
    each substep (the electrons' share to their energy, the ions' left in
    _last_ion_work_j_m3 for the gas), which keeps the field's energy to
    round-off. It needs both electrodes grounded and no drive."""

    FIELD_WORK = ("fluid1d", "conserve")

    def __init__(self, *, mesh=None, field_work="fluid1d", **kw):
        if field_work not in self.FIELD_WORK:
            raise ValueError(f"field_work is one of {self.FIELD_WORK}")
        super().__init__(**kw)
        self.geo = None if mesh is None else MeshArrays(mesh)
        self.conserve_field_work = field_work == "conserve"
        self._last_ion_work_j_m3 = None
        self._last_field_energy_j_m2 = (0.0, 0.0)
        self._last_field_work_j_m2 = 0.0
        self._last_positivity_limits = 0

    def set_mesh(self, mesh):
        """The mesh, once the orchestrator has made it (before any step)."""
        self.geo = MeshArrays(mesh)

    def field_energy_j_m2(self, charged_fluid) -> float:
        """The field's energy per unit area [J/m^2] of the charge as it
        stands, between grounded electrodes (the kernel's exact field)."""
        g = self.geo
        cf = charged_fluid
        n = g.widths.size
        q = np.array([s.charge_number for s in cf.species],
                     dtype=np.float64) * E_CHARGE
        e = np.empty(n + 1)
        rows = [np.empty(n + 1)] + [np.empty(n) for _ in range(7)]
        _exact_field(np.ascontiguousarray(cf.density, dtype=np.float64), q,
                     g.widths, g.spacing, 0.0, 0.0, *rows, e)
        return float(_field_energy(e, g.widths, g.spacing))

    def advance(self, dt, *, efield, neutral_fluid, charged_fluid, coeffs,
                field=None) -> StepReport:
        if self.geo is None:
            raise RuntimeError("MeshPlasmaModel has no mesh: set_mesh first")
        if self.geo.widths.size != charged_fluid.n_cells:
            raise RuntimeError("the transport's mesh is not the fluid's")
        if not self._resolved:
            self._resolve(neutral_fluid, charged_fluid, coeffs)
        if field is None or getattr(field, "prescribed", None) is not None:
            raise NotImplementedError(
                "1Dv0_3's transport runs the coupled semi-implicit step only")
        return self._advance_on_mesh(dt, efield, neutral_fluid,
                                     charged_fluid, field)

    def _advance_on_mesh(self, dt, efield, neutral_fluid, charged_fluid,
                         field) -> StepReport:
        """Fluid1D's _advance_via_kernel with coupled_advance on the mesh."""
        cf = charged_fluid
        if self.conserve_field_work and not (
                field.phi_left == 0.0 and field.phi_right == 0.0
                and getattr(field, "gap_drive_a", 0.0) == 0.0
                and getattr(field, "gap_capacitance_f", None) is None):
            raise NotImplementedError(
                "field_work='conserve' keeps the field's energy between "
                "grounded electrodes; a potential difference or a drive "
                "needs the circuit's work in its books")
        if self._kernel_packed is None or self._kernel_packed[0] \
                is not self._me_coeffs:
            self._kernel_packed = (
                self._me_coeffs,
                _fk.pack_kernel_args(self, self._me_coeffs, neutral_fluid, cf))
        result = coupled_advance(
            self, self.geo, dt, np.asarray(efield, dtype=np.float64)[:, 0],
            neutral_fluid, cf, self._kernel_packed[1],
            phi_l=field.phi_left, phi_r=field.phi_right,
            c_g=getattr(field, "gap_capacitance_f", None),
            i_drive=getattr(field, "gap_drive_a", 0.0))
        field.v_b_final = result.v_b_final_v
        field.phi_right = result.v_b_final_v
        if result.capped:
            warnings.warn("transport kernel hit max_substeps; its final "
                          "remainder step exceeded the requested CFL target",
                          stacklevel=2)
        eps_final = self._recover_eps(result.n_eps_j_m3,
                                      cf.density[:, self._e_col],
                                      result.equilibrium_energy_ev)
        cf.electron_mean_energy_ev[:] = eps_final
        self._last_n_eps = np.maximum(result.n_eps_j_m3, 0.0)
        self._last_substeps = result.steps
        lef = np.zeros((result.e_axial_v_m.shape[0], 3))
        lef[:, 0] = result.e_axial_v_m
        field.last_efield = lef
        n = cf.n_cells
        return StepReport(
            deposited_power_w=result.deposited_j / dt if dt > 0 else 0.0,
            mean_electron_energy_ev=float(eps_final.mean()),
            dt=dt,
            conduction_current_a=-self.area_m2 * result.jcond_mean_a_m2,
            gap_conductance_s=(max(result.gap_energy_absorbed_j, 0.0)
                               / max(result.gap_voltage_squared_dt_v2_s,
                                     1e-30)),
            gap_energy_absorbed_j=max(result.gap_energy_absorbed_j, 0.0),
            energy_closure=EnergyClosureReport(
                field_work_j_m3=result.field_work_j_m3,
                elastic_transfer_j_m3=result.elastic_transfer_j_m3,
                continuous_inelastic_sink_j_m3=result.continuous_sink_j_m3,
                accepted_reaction_cost_j_m3=np.zeros(n),
                signed_reconciliation_j_m3=np.zeros(n),
                coefficient_bound_adjustment_j_m3=result.bound_adjustment_j_m3,
                maximum_energy_jump_ev=result.maximum_energy_jump_ev,
                limiter_activity=result.limiter_activity,
                max_ampere_residual=result.max_ampere_residual,
                max_energy_residual_j_m3=result.max_energy_residual_j_m3))


class MeshPoisson(Poisson1DSolver):
    """Fluid1D's cell-centred semi-implicit Poisson solve on the mesh: the
    rows of _solve_field above, in NumPy (Fluid1D's own rows on a uniform
    mesh). The orchestrator holds one for its field handle."""

    def __init__(self, mesh, *, eps0: float = EPS0):
        super().__init__(mesh.n_cells, None, eps0=eps0)
        self.geo = MeshArrays(mesh)

    def solve_semi_implicit(self, charge_density_1d, phi_left, phi_right, *,
                            sigma_cells, j_expl_faces, h):
        g = self.geo
        n = self.n_cells
        rho = np.asarray(charge_density_1d, dtype=np.float64)
        a = self.eps0 + h * np.maximum(np.asarray(sigma_cells, float), 0.0)
        jf = np.asarray(j_expl_faces, dtype=np.float64)
        af = np.empty(n + 1)
        af[0] = a[0] * (g.dx_ref / (0.5 * g.widths[0]))
        af[n] = a[-1] * (g.dx_ref / (0.5 * g.widths[-1]))
        af[1:n] = (g.wl * a[:-1] + g.wr * a[1:]) * (g.dx_ref / g.spacing)
        diag = af[:-1] + af[1:]
        lower = np.empty(n)
        upper = np.empty(n)
        lower[0], upper[-1] = 0.0, 0.0
        lower[1:] = -af[1:n]
        upper[:-1] = -af[1:n]
        rhs = g.dx_ref * g.widths * rho - h * g.dx_ref * (jf[1:] - jf[:-1])
        rhs[0] += af[0] * phi_left
        rhs[-1] += af[n] * phi_right
        from Fluid1D.field import _thomas_general
        phi = _thomas_general(lower, diag, upper, rhs)
        e_faces = np.empty(n + 1)
        e_faces[1:n] = -(phi[1:] - phi[:-1]) / g.spacing
        e_faces[0] = -(phi[0] - phi_left) / (0.5 * g.widths[0])
        e_faces[n] = -(phi_right - phi[-1]) / (0.5 * g.widths[-1])
        self.last_e_faces = e_faces
        return phi, 0.5 * (e_faces[:-1] + e_faces[1:])
