"""Numba kernel for the coupled (semi-implicit-field) transport substep loop.

The post-breakdown cost of :meth:`DriftDiffusionPlasmaModel.advance` is the pure-Python
substep loop: at nₑ ~ 10²² the single-cell sheath field forces ~10²–10³ CFL substeps per
outer step, each one ~10² small numpy calls (interpreter overhead) plus the pure-Python
Thomas solve — measured ~2–3 s per 50 ps step. This module runs the **whole substep loop in
one ``@njit`` kernel** — ⟨ε⟩ recovery/blend, LUT lookups (``np.interp`` on packed arrays),
per-species Scharfetter–Gummel fluxes, boundary fluxes + secondary emission, the
semi-implicit σ-corrected tridiagonal field solve, and the electron energy PDE — replicating
the numpy path's arithmetic and substep sequencing
(parity is pinned by ``tests/test_transport_kernel.py``).

Only the coupled path lives here; the ``field=None`` frozen-field path stays numpy (it is
the bit-for-bit reference baseline for the stub tests).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particles.constants import E_CHARGE, EPS0, K_B, M_E

try:
    from numba import njit
except ImportError:  # pragma: no cover - numba is a hard runtime dep, this is a guard
    def njit(*args, **kwargs):  # type: ignore[misc]
        def deco(fn):
            return fn
        return deco(args[0]) if args and callable(args[0]) else deco


@njit(cache=True)
def _bernoulli(x: float) -> float:
    if abs(x) < 1.0e-10:
        return 1.0 - 0.5 * x
    return x / np.expm1(x)


@njit(cache=True)
def _recover_energy(
    density, n_eps, eps, eps_eq, e_col, n_e_floor, has_eq, lmea_ref, eps_lo, eps_hi,
):
    """Recover and bound mean energy, repairing the near-vacuum LMEA ratio in place."""
    for i in range(density.shape[0]):
        n_e_i = density[i, e_col]
        n_e_safe = n_e_i if n_e_i > n_e_floor else n_e_floor
        eps_lmea = n_eps[i] / (n_e_safe * E_CHARGE)
        if has_eq:
            weight = n_e_i / (n_e_i + lmea_ref)
            value = weight * eps_lmea + (1.0 - weight) * eps_eq[i]
        else:
            value = eps_lmea
        if not np.isfinite(value):
            value = eps_lo
        if value < eps_lo:
            value = eps_lo
        elif value > eps_hi:
            value = eps_hi
        eps[i] = value
        n_eps[i] = n_e_safe * value * E_CHARGE


@njit(cache=True)
def _lookup_transport(
    eps, n_gas_safe, key_grid, mob_arr, dif_arr, emob_arr, edif_arr,
    e_col, ion_mu_n, kt_gas_ev, mu_e, d_e, mu_eps, d_eps, mu_s, d_face_s,
):
    """Interpolate electron rows and construct ion mobility/diffusion face data."""
    nk = key_grid.shape[0]
    for i in range(eps.shape[0]):
        value = eps[i]
        j = np.searchsorted(key_grid, value)
        if j <= 0:
            jl = 0
            weight = 0.0
        elif j >= nk:
            jl = nk - 2
            weight = 1.0
        else:
            jl = j - 1
            weight = (value - key_grid[jl]) / (key_grid[jl + 1] - key_grid[jl])
        lower_weight = 1.0 - weight
        mu_e[i] = (lower_weight * mob_arr[jl] + weight * mob_arr[jl + 1]) / n_gas_safe[i]
        d_e[i] = (lower_weight * dif_arr[jl] + weight * dif_arr[jl + 1]) / n_gas_safe[i]
        mu_eps[i] = (lower_weight * emob_arr[jl] + weight * emob_arr[jl + 1]) / n_gas_safe[i]
        d_eps[i] = (lower_weight * edif_arr[jl] + weight * edif_arr[jl + 1]) / n_gas_safe[i]
    for s in range(mu_s.shape[0]):
        if s == e_col:
            for i in range(eps.shape[0]):
                mu_s[s, i] = mu_e[i]
            for i in range(eps.shape[0] - 1):
                d_face_s[s, i] = 0.5 * (d_e[i] + d_e[i + 1])
        else:
            for i in range(eps.shape[0]):
                mu_s[s, i] = ion_mu_n[s] / n_gas_safe[i]
            for i in range(eps.shape[0] - 1):
                d_face_s[s, i] = 0.5 * kt_gas_ev * (mu_s[s, i] + mu_s[s, i + 1])


@njit(cache=True)
def _cfl_timestep(
    e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col, dx, cfl_safety, remaining,
):
    """Return the drift/diffusion/thermal-wall CFL substep for the current state."""
    v_max = 0.0
    d_max = 0.0
    for i in range(eps.shape[0]):
        value = (5.0 / 3.0) * mu_e[i] * abs(e_axial[i])
        if value > v_max:
            v_max = value
        if d_e[i] > d_max:
            d_max = d_e[i]
        if d_eps[i] > d_max:
            d_max = d_eps[i]
    for s in range(mu_s.shape[0]):
        if s == e_col:
            continue
        for i in range(eps.shape[0]):
            value = mu_s[s, i] * abs(e_axial[i])
            if value > v_max:
                v_max = value
    vth0 = 0.25 * np.sqrt(
        8.0 * ((2.0 / 3.0) * max(eps[0], 0.0) * E_CHARGE) / (np.pi * M_E)
    )
    vthl = 0.25 * np.sqrt(
        8.0 * ((2.0 / 3.0) * max(eps[-1], 0.0) * E_CHARGE) / (np.pi * M_E)
    )
    if vth0 > v_max:
        v_max = vth0
    if vthl > v_max:
        v_max = vthl
    dt_cfl = 1.0e300
    if v_max > 0.0:
        dt_cfl = dx / v_max
    if d_max > 0.0 and dx * dx / (2.0 * d_max) < dt_cfl:
        dt_cfl = dx * dx / (2.0 * d_max)
    h = cfl_safety * dt_cfl
    if not np.isfinite(h) or h <= 0.0:
        h = remaining
    return min(h, remaining)


@njit(cache=True)
def _solve_field(
    density, n_gas_safe, h, dx, phi_l, phi_r, e_col, q, mu_e, mu_s, d_face_s,
    has_eq, eq_en, eq_eps, eps_lo, eps_hi, sigma_cell, j_expl_f, rho, a_coef,
    a_face, e_face, diag, lower, upper, rhs, cp, dp, phi_cells, e_axial, eps_eq,
):
    """Assemble and solve the cell-centred semi-implicit Poisson system in place."""
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
                -d_face_s[s, i] * (density[i + 1, s] - density[i, s]) / dx
            )
    for i in range(n):
        a_coef[i] = EPS0 + h * max(sigma_cell[i], 0.0)
    a_face[0] = a_coef[0]
    a_face[n] = a_coef[n - 1]
    for face in range(1, n):
        a_face[face] = 0.5 * (a_coef[face - 1] + a_coef[face])
    for i in range(n):
        diag[i] = a_face[i] + a_face[i + 1]
        lower[i] = -a_face[i]
        upper[i] = -a_face[i + 1]
        rhs[i] = dx * dx * rho[i] - h * dx * (j_expl_f[i + 1] - j_expl_f[i])
    diag[0] = 2.0 * a_face[0] + a_face[1]
    diag[n - 1] = a_face[n - 1] + 2.0 * a_face[n]
    rhs[0] += 2.0 * a_face[0] * phi_l
    rhs[n - 1] += 2.0 * a_face[n] * phi_r
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denominator = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denominator
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denominator
    phi_cells[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        phi_cells[i] = dp[i] - cp[i] * phi_cells[i + 1]
    e_face[0] = -(phi_cells[0] - phi_l) / (0.5 * dx)
    e_face[n] = -(phi_r - phi_cells[n - 1]) / (0.5 * dx)
    for face in range(1, n):
        e_face[face] = -(phi_cells[face] - phi_cells[face - 1]) / dx
    for i in range(n):
        e_axial[i] = 0.5 * (e_face[i] + e_face[i + 1])
        if has_eq:
            en_td = abs(e_axial[i]) / n_gas_safe[i] / 1.0e-21
            value = np.interp(en_td, eq_en, eq_eps)
            eps_eq[i] = min(max(value, eps_lo), eps_hi)


@njit(cache=True)
def _transport_species(
    density, eps, e_face, h, dx, e_col, z, q, mu_s, d_face_s, ion_vth,
    pos_ion, gamma_se, j_left, j_right, interior, gam, faces, rho_before,
):
    """Compute SG/electrode fluxes and conservatively update all charged densities."""
    n = density.shape[0]
    s_tot = density.shape[1]
    for i in range(n):
        charge = 0.0
        for s in range(s_tot):
            charge += q[s] * density[i, s]
        rho_before[i] = charge
    for s in range(s_tot):
        for i in range(n - 1):
            velocity = z[s] * 0.5 * (mu_s[s, i] + mu_s[s, i + 1]) * e_face[i + 1]
            diffusion = d_face_s[s, i]
            if abs(diffusion) < 1.0e-30:
                diffusion = 1.0e-30
            peclet = velocity * dx / diffusion
            interior[s, i] = (diffusion / dx) * (
                _bernoulli(-peclet) * density[i, s]
                - _bernoulli(peclet) * density[i + 1, s]
            )
        if s == e_col:
            vth0 = np.sqrt(
                8.0 * ((2.0 / 3.0) * max(eps[0], 0.0) * E_CHARGE) / (np.pi * M_E)
            )
            vthl = np.sqrt(
                8.0 * ((2.0 / 3.0) * max(eps[-1], 0.0) * E_CHARGE) / (np.pi * M_E)
            )
        else:
            vth0 = ion_vth[s]
            vthl = ion_vth[s]
        v0 = z[s] * mu_s[s, 0] * e_face[0]
        vl = z[s] * mu_s[s, n - 1] * e_face[n]
        j_left[s] = (
            -0.25 * density[0, s] * vth0
            + (v0 if v0 < 0.0 else 0.0) * density[0, s]
        )
        j_right[s] = (
            0.25 * density[n - 1, s] * vthl
            + (vl if vl > 0.0 else 0.0) * density[n - 1, s]
        )
    ion_cathode_flux = 0.0
    for index in range(pos_ion.shape[0]):
        value = j_right[pos_ion[index]]
        if value > 0.0:
            ion_cathode_flux += value
    j_right[e_col] -= gamma_se * ion_cathode_flux

    for s in range(s_tot):
        faces[0] = j_left[s]
        faces[n] = j_right[s]
        for i in range(n - 1):
            faces[i + 1] = interior[s, i]
        for i in range(n):
            gam[s, i] = 0.5 * (faces[i] + faces[i + 1])
            value = density[i, s] - h * (faces[i + 1] - faces[i]) / dx
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
        current_div = (j_face_r - j_face_l) / dx
        ampere_abs = abs(charge_rate + current_div)
        ampere_scale = max(abs(charge_rate), abs(current_div), 1.0e-300)
        relative = ampere_abs / ampere_scale
        if relative > maximum_residual:
            maximum_residual = relative
    return maximum_residual


@njit(cache=True)
def _advance_gap(v_b, c_g, i_drive, area, gap, h, jcond_mean, sigma_mean):
    """Advance the gap-only Lie split and return final voltage, supply voltage, and energy."""
    v_supply = v_b
    i_gap = -area * jcond_mean
    if np.isfinite(c_g):
        conductance = area * sigma_mean / gap
        offset = i_gap - conductance * v_b
        if conductance > 0.0:
            v_eq = (i_drive - offset) / conductance
            v_b = v_eq + (v_b - v_eq) * np.exp(-conductance * h / c_g)
        else:
            v_b = v_b + (i_drive - offset) * h / c_g
        energy_absorbed = 0.5 * c_g * (v_supply * v_supply - v_b * v_b)
    else:
        energy_absorbed = i_gap * v_supply * h
    return v_b, v_supply, energy_absorbed


@njit(cache=True)
def _inelastic_power_cell(
    cell, eps, density, neutral, key_grid, te_arr, loss_ev, base_loss_k,
    loss_rate_kind, loss_rate_lut, loss_analytic, loss_react_kind,
    loss_react_col, loss_react_order, loss_react_count,
):
    """Evaluate the continuous diagnostic reaction-energy sink for one cell."""
    nk = key_grid.shape[0]
    lower = np.searchsorted(key_grid, eps[cell])
    if lower <= 0:
        lower = 0
        weight = 0.0
    elif lower >= nk:
        lower = nk - 2
        weight = 1.0
    else:
        lower -= 1
        weight = (eps[cell] - key_grid[lower]) / (key_grid[lower + 1] - key_grid[lower])
    te = (
        (1.0 - weight) * te_arr[lower] + weight * te_arr[lower + 1]
        if nk > 1 else te_arr[0]
    )
    power = 0.0
    for reaction in range(loss_ev.shape[0]):
        kind = loss_rate_kind[reaction]
        if kind == 1:
            rate = (
                (1.0 - weight) * loss_rate_lut[reaction, lower]
                + weight * loss_rate_lut[reaction, lower + 1]
                if nk > 1 else loss_rate_lut[reaction, 0]
            )
        elif kind == 2:
            k0 = loss_analytic[reaction, 0]
            exponent = loss_analytic[reaction, 1]
            threshold = loss_analytic[reaction, 2]
            rate = k0 * te**exponent * np.exp(-threshold / te) if te > 0.0 else 0.0
        else:
            rate = base_loss_k[reaction]
        progress = rate
        for reactant in range(loss_react_count[reaction]):
            column = loss_react_col[reaction, reactant]
            reactant_density = (
                density[cell, column]
                if loss_react_kind[reaction, reactant] == 0
                else neutral[cell, column]
            )
            progress *= reactant_density ** loss_react_order[reaction, reactant]
        power += loss_ev[reaction] * progress * E_CHARGE
    return power


@njit(cache=True)
def _transport_energy(
    density, neutral, n_eps, eps, e_face, e_axial, mu_eps, d_eps, mu_e, gam,
    j_left, j_right, e_col, dx, h, area, substep_current, v_b_supply,
    heating_cap_factor, m_gas, eps_gas_ev, closure_mode, key_grid, te_arr,
    loss_ev, base_loss_k, loss_rate_kind, loss_rate_lut, loss_analytic,
    loss_react_kind, loss_react_col, loss_react_order, loss_react_count,
    eps_floor, eps_ceiling, n_e_floor, faces, pf_arr, field_work,
    elastic_transfer, continuous_sink, bound_adjustment,
):
    """Advance the energy PDE and return this substep's scalar ledger extrema."""
    n = density.shape[0]
    faces[0] = eps[0] * E_CHARGE * j_left[e_col]
    faces[n] = eps[n - 1] * E_CHARGE * j_right[e_col]
    for i in range(n - 1):
        velocity = -0.5 * (mu_eps[i] + mu_eps[i + 1]) * e_face[i + 1]
        diffusion = 0.5 * (d_eps[i] + d_eps[i + 1])
        if abs(diffusion) < 1.0e-30:
            diffusion = 1.0e-30
        peclet = velocity * dx / diffusion
        faces[i + 1] = (diffusion / dx) * (
            _bernoulli(-peclet) * n_eps[i] - _bernoulli(peclet) * n_eps[i + 1]
        )

    positive_power = 0.0
    negative_power = 0.0
    for i in range(n):
        pf_arr[i] = -E_CHARGE * gam[e_col, i] * e_axial[i]
        if pf_arr[i] > 0.0:
            positive_power += pf_arr[i]
        else:
            negative_power -= pf_arr[i]
    positive_power *= dx * area
    negative_power *= dx * area
    cap_w = heating_cap_factor * abs(area * substep_current * v_b_supply) + 1.0e3
    magnitude = positive_power if positive_power > negative_power else negative_power
    limiter_activity = 0
    if magnitude > cap_w:
        scale = cap_w / magnitude
        for i in range(n):
            pf_arr[i] *= scale
        limiter_activity += 1

    deposited_dx_j = 0.0
    maximum_jump = 0.0
    maximum_residual = 0.0
    for i in range(n):
        transport_rate = -(faces[i + 1] - faces[i]) / dx
        p_field = pf_arr[i]
        collision_rate = E_CHARGE / (M_E * mu_e[i]) if mu_e[i] > 0.0 else 0.0
        p_elastic = (
            (2.0 * M_E / m_gas[i])
            * collision_rate
            * density[i, e_col]
            * (eps[i] - eps_gas_ev)
            * E_CHARGE
        )
        p_inelastic = 0.0
        if closure_mode == 1:
            p_inelastic = _inelastic_power_cell(
                i, eps, density, neutral, key_grid, te_arr, loss_ev, base_loss_k,
                loss_rate_kind, loss_rate_lut, loss_analytic, loss_react_kind,
                loss_react_col, loss_react_order, loss_react_count,
            )
        old_energy = n_eps[i]
        local_power = p_field - p_elastic - p_inelastic if closure_mode == 1 else 0.0
        raw_energy = old_energy + h * (transport_rate + local_power)
        lower_bound = eps_floor * density[i, e_col] * E_CHARGE
        upper_bound = (
            eps_ceiling * max(density[i, e_col], 0.0) * E_CHARGE + 1.0e-300
        )
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
        residual = (
            next_energy - old_energy - h * (transport_rate + local_power) - adjustment
        )
        if abs(residual) > maximum_residual:
            maximum_residual = abs(residual)
        deposited_dx_j += p_field * dx * h
    return deposited_dx_j, maximum_jump, limiter_activity, maximum_residual


@njit(cache=True)
def _coupled_loop(  # noqa: PLR0913 - a kernel takes what a kernel takes
    density,            # (N, S) charged densities, mutated in place
    n_eps,              # (N,) electron energy density J/m^3, mutated in place
    e_axial,            # (N,) axial field V/m, mutated in place (bootstrap in, final out)
    neutral,            # (N, Sn) neutral densities (read-only)
    n_gas_safe,         # (N,)
    m_gas,              # (N,)
    dt, dx, phi_l, phi_r, c_g, i_drive, area, gap,
    cfl_safety, max_substeps, eps_gas_ev, kt_gas_ev, gamma_se, heating_cap_factor,
    eps_lo, eps_hi, lmea_ref, n_e_floor, eps_floor, eps_ceiling,
    e_col, z, q, ion_mu_n, ion_vth, pos_ion,   # species meta: (S,)×4 + (P,)
    key_grid, mob_arr, dif_arr, emob_arr, edif_arr,   # electron LUT rows on key_grid
    has_eq, eq_en, eq_eps,                     # equilibrium <eps>(E/N Td) table
    closure_mode, loss_ev, base_loss_k, loss_rate_kind, loss_rate_lut, loss_analytic,
    loss_react_kind, loss_react_col, loss_react_order, loss_react_count, te_arr,
):
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
    phi_nodes = np.empty(n + 1)
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

    for i in range(n):  # eq-<eps> for the bootstrap field (first recovery's blend)
        if has_eq:
            en_td = abs(e_axial[i]) / n_gas_safe[i] / 1.0e-21
            v = np.interp(en_td, eq_en, eq_eps)
            eps_eq[i] = min(max(v, eps_lo), eps_hi)
        else:
            eps_eq[i] = 0.0

    deposited_dx_j = 0.0   # deposited energy per unit AREA·? — Σ p_field·dx·h; ×area by caller
    jcond_mean = 0.0
    gap_e_abs = 0.0
    gap_v2 = 0.0
    remaining = dt
    steps_done = 0
    capped = False
    h_next = -1.0  # <0 → compute from the CFL at the loop top (first iteration only)
    # Gap voltage is co-advanced per substep for finite c_g so plasma draw reduces the
    # cathode Dirichlet value between circuit exchanges. c_g = inf gives fixed voltage.
    v_b = phi_r

    while True:
        _recover_energy(
            density, n_eps, eps, eps_eq, e_col, n_e_floor,
            has_eq, lmea_ref, eps_lo, eps_hi,
        )
        _lookup_transport(
            eps, n_gas_safe, key_grid, mob_arr, dif_arr, emob_arr, edif_arr,
            e_col, ion_mu_n, kt_gas_ev, mu_e, d_e, mu_eps, d_eps, mu_s, d_face_s,
        )

        # --- substep size (numpy-path sequencing: first from the bootstrap state, then
        #     the value computed at the END of the previous iteration) -------------------
        if h_next < 0.0:
            h = _cfl_timestep(
                e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col,
                dx, cfl_safety, remaining,
            )
        else:
            h = h_next
        if steps_done >= max_substeps:
            capped = True
            h = remaining

        _solve_field(
            density, n_gas_safe, h, dx, phi_l, v_b, e_col, q, mu_e, mu_s,
            d_face_s, has_eq, eq_en, eq_eps, eps_lo, eps_hi, sigma_cell,
            j_expl_f, rho, a_coef, a_face, e_face, diag, lower, upper, rhs,
            cp, dp, phi_nodes, e_axial, eps_eq,
        )

        ampere_residual = _transport_species(
            density, eps, e_face, h, dx, e_col, z, q, mu_s, d_face_s,
            ion_vth, pos_ion, gamma_se, j_left, j_right, interior, gam,
            faces, rho_before,
        )
        if ampere_residual > max_ampere_residual:
            max_ampere_residual = ampere_residual

        jc = 0.0
        sg = 0.0
        for i in range(n):
            jci = 0.0
            for s in range(s_tot):
                jci += q[s] * gam[s, i]
            jc += jci
            sg += sigma_cell[i]
        substep_current = jc / n
        substep_sigma = sg / n
        jcond_mean += substep_current * (h / dt)
        v_b, v_b_supply, gap_energy = _advance_gap(
            v_b, c_g, i_drive, area, gap, h, substep_current, substep_sigma,
        )
        gap_v2 += v_b_supply * v_b_supply * h
        gap_e_abs += gap_energy

        energy_deposited, energy_jump, energy_limiters, energy_residual = _transport_energy(
            density, neutral, n_eps, eps, e_face, e_axial, mu_eps, d_eps, mu_e,
            gam, j_left, j_right, e_col, dx, h, area, substep_current, v_b_supply,
            heating_cap_factor, m_gas, eps_gas_ev, closure_mode, key_grid, te_arr,
            loss_ev, base_loss_k, loss_rate_kind, loss_rate_lut, loss_analytic,
            loss_react_kind, loss_react_col, loss_react_order, loss_react_count,
            eps_floor, eps_ceiling, n_e_floor, faces, pf_arr, field_work,
            elastic_transfer, continuous_sink, bound_adjustment,
        )
        deposited_dx_j += energy_deposited
        if energy_jump > maximum_energy_jump_ev:
            maximum_energy_jump_ev = energy_jump
        limiter_activity += energy_limiters
        if energy_residual > max_energy_residual:
            max_energy_residual = energy_residual

        steps_done += 1
        remaining -= h
        if remaining <= 1.0e-12 * dt or capped:
            break

        h_next = _cfl_timestep(
            e_axial, eps, mu_e, d_e, d_eps, mu_s, e_col,
            dx, cfl_safety, remaining,
        )

    return (jcond_mean, deposited_dx_j, steps_done, capped,
            eps_eq, v_b, gap_e_abs, gap_v2, field_work, elastic_transfer,
            continuous_sink, bound_adjustment, maximum_energy_jump_ev,
            limiter_activity, max_ampere_residual, max_energy_residual)


# ---------------------------------------------------------------------------- caller side
@dataclass(frozen=True)
class _CoupledAdvanceResult:
    """Named Python boundary for the positional Numba kernel ABI."""

    jcond_mean_a_m2: float
    deposited_j: float
    n_eps_j_m3: np.ndarray
    e_axial_v_m: np.ndarray
    steps: int
    capped: bool
    equilibrium_energy_ev: np.ndarray
    v_b_final_v: float
    gap_energy_absorbed_j: float
    gap_voltage_squared_dt_v2_s: float
    field_work_j_m3: np.ndarray
    elastic_transfer_j_m3: np.ndarray
    continuous_sink_j_m3: np.ndarray
    bound_adjustment_j_m3: np.ndarray
    maximum_energy_jump_ev: float
    limiter_activity: int
    max_ampere_residual: float
    max_energy_residual_j_m3: float


def pack_kernel_args(model, me, neutral_fluid, charged_fluid):
    """Pack a resolved :class:`DriftDiffusionPlasmaModel` + mean-energy table for the kernel.

    Everything static across one ``advance`` call: species meta, LUT rows evaluated on the
    table's own key grid (the effective energy pair — including the 5/3 closure when the
    table carries none), and the equilibrium ⟨ε⟩(E/N) inversion.
    """
    cf = charged_fluid
    key = np.ascontiguousarray(me.key_grid, dtype=np.float64)
    mob = np.ascontiguousarray(me.mobility(key), dtype=np.float64)
    dif = np.ascontiguousarray(me.diffusion(key), dtype=np.float64)
    emob = np.ascontiguousarray(me.energy_mobility(key), dtype=np.float64)
    edif = np.ascontiguousarray(me.energy_diffusion(key), dtype=np.float64)
    te_arr = np.ascontiguousarray(me.te(key), dtype=np.float64)
    en = getattr(me, "reduced_field_td", None)
    if en is not None:
        order = np.argsort(np.asarray(en), kind="stable")
        eq_en = np.ascontiguousarray(np.asarray(en, dtype=np.float64)[order])
        eq_eps = np.ascontiguousarray(np.asarray(me.key_grid, dtype=np.float64)[order])
        has_eq = True
    else:
        eq_en = np.zeros(1)
        eq_eps = np.zeros(1)
        has_eq = False

    kt_ev = K_B * model.gas_temp_k / E_CHARGE
    ion_vth = np.zeros(cf.n_species)
    for s in range(cf.n_species):
        if s != model._e_col:
            kt_j = (2.0 / 3.0) * (kt_ev * 1.5) * E_CHARGE
            ion_vth[s] = float(np.sqrt(8.0 * kt_j / (np.pi * model._ion_masses[s])))

    packed = model._inelastic_packed
    if packed is None:
        raise RuntimeError("transport reaction-energy table was not resolved")
    return dict(
        eps_gas_ev=1.5 * kt_ev, kt_gas_ev=kt_ev, gamma_se=model.gamma,
        eps_lo=model._eps_bounds[0], eps_hi=model._eps_bounds[1],
        lmea_ref=model.lmea_ref_density_m3, n_e_floor=1.0,
        eps_floor=model.eps_floor, eps_ceiling=model.eps_ceiling,
        e_col=model._e_col,
        z=np.ascontiguousarray(model._z), q=np.ascontiguousarray(model._q),
        ion_mu_n=np.ascontiguousarray(model._ion_mu_N), ion_vth=ion_vth,
        pos_ion=np.asarray(model._pos_ion_cols, dtype=np.int64),
        key_grid=key, mob_arr=mob, dif_arr=dif, emob_arr=emob, edif_arr=edif,
        te_arr=te_arr,
        has_eq=has_eq, eq_en=eq_en, eq_eps=eq_eps,
        loss_ev=packed["loss_ev"], base_loss_k=packed["base_k"],
        loss_rate_kind=packed["rate_kind"], loss_rate_lut=packed["rate_lut"],
        loss_analytic=packed["analytic"], loss_react_kind=packed["react_kind"],
        loss_react_col=packed["react_col"], loss_react_order=packed["react_order"],
        loss_react_count=packed["react_count"],
    )


def coupled_advance(model, dt, e_axial0, neutral_fluid, charged_fluid, packed, *,
                    phi_l, phi_r, c_g=None, i_drive=0.0):
    """Run one coupled ``advance`` via the kernel; mutates the charged fluid in place.

    Returns current/conductivity, deposited energy, electron energy/field, substep status,
    and final gap state in the same conventions as the NumPy path.
    """
    cf = charged_fluid
    n_gas = neutral_fluid.density.sum(axis=1)
    n_gas_safe = np.where(n_gas > 0.0, n_gas, 1.0)
    masses = model._neutral_masses
    m_gas = np.where(n_gas > 0.0, (neutral_fluid.density @ masses) / n_gas_safe, masses.mean())

    e_col = model._e_col
    n_e0 = cf.density[:, e_col]
    n_eps = np.ascontiguousarray(cf.electron_mean_energy_ev * n_e0 * E_CHARGE)
    e_axial = np.ascontiguousarray(e_axial0, dtype=np.float64).copy()
    density = np.ascontiguousarray(cf.density)

    p = packed
    (jcond_mean, dep_dx, steps, capped, eps_eq,
     v_b_final, gap_e_abs, gap_v2, field_work, elastic_transfer,
     continuous_sink, bound_adjustment, max_jump, limiter_activity,
     ampere_residual, energy_residual) = _coupled_loop(
        density, n_eps, e_axial, np.ascontiguousarray(neutral_fluid.density),
        np.ascontiguousarray(n_gas_safe), np.ascontiguousarray(m_gas),
        dt, cf.dx, phi_l, phi_r,
        np.inf if c_g is None else float(c_g), float(i_drive), model.area_m2, model.gap_m,
        model.cfl_safety, model.max_substeps, p["eps_gas_ev"], p["kt_gas_ev"], p["gamma_se"],
        model.field_heating_cap_factor,
        p["eps_lo"], p["eps_hi"], p["lmea_ref"], p["n_e_floor"], p["eps_floor"],
        p["eps_ceiling"],
        p["e_col"], p["z"], p["q"], p["ion_mu_n"], p["ion_vth"], p["pos_ion"],
        p["key_grid"], p["mob_arr"], p["dif_arr"], p["emob_arr"], p["edif_arr"],
        p["has_eq"], p["eq_en"], p["eq_eps"],
        1 if model.electron_energy_closure == "continuous_reconciled" else 0,
        p["loss_ev"], p["base_loss_k"], p["loss_rate_kind"], p["loss_rate_lut"],
        p["loss_analytic"], p["loss_react_kind"], p["loss_react_col"],
        p["loss_react_order"], p["loss_react_count"], p["te_arr"],
    )
    cf.density[:] = density
    deposited_j = dep_dx * cf.area_m2  # kernel accumulates Σ p_field·dx·h; volume = dx·area
    return _CoupledAdvanceResult(
        jcond_mean_a_m2=float(jcond_mean),
        deposited_j=float(deposited_j),
        n_eps_j_m3=n_eps,
        e_axial_v_m=e_axial,
        steps=int(steps),
        capped=bool(capped),
        equilibrium_energy_ev=eps_eq,
        v_b_final_v=float(v_b_final),
        gap_energy_absorbed_j=float(gap_e_abs),
        gap_voltage_squared_dt_v2_s=float(gap_v2),
        field_work_j_m3=field_work,
        elastic_transfer_j_m3=elastic_transfer,
        continuous_sink_j_m3=continuous_sink,
        bound_adjustment_j_m3=bound_adjustment,
        maximum_energy_jump_ev=float(max_jump),
        limiter_activity=int(limiter_activity),
        max_ampere_residual=float(ampere_residual),
        max_energy_residual_j_m3=float(energy_residual),
    )
