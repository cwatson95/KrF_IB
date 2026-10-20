"""Numba per-cell implicit stiff integrator for combined plasma chemistry.

Each grid cell holds an independent species-density vector; the mass-action ODE
``dn/dt = ν @ R(n)`` is stiff (radiative 1/τ ~ 4e7 s⁻¹ against 3-body rates), so we
take **linearly-implicit (Rosenbrock-)Euler** substeps — one linear solve
``(I − dt·J)·δ = dt·f`` per substep — which is L-stable for the linear part. The step
runs in a ``prange`` over cells (embarrassingly parallel; cells never interact here).
The analytic ``f`` and ``J`` are assembled from the padded reactant tables built by
:class:`~kinetics.chemistry.FluidChemistry`. SciPy Radau (in ``chemistry.py``) is the
single-cell reference/oracle this step is validated against.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from particles.constants import E_CHARGE

from .chemistry import FluidChemistry

try:
    from numba import njit, prange

    HAS_NUMBA = True
except ImportError:  # pragma: no cover - numba is a hard runtime dep, this is a guard
    HAS_NUMBA = False

    def njit(*args, **kwargs):  # type: ignore[misc]
        def deco(fn):
            return fn
        return deco(args[0]) if args and callable(args[0]) else deco

    def prange(*args):  # type: ignore[misc]
        return range(*args)


@njit(cache=True)
def _progress_only(n, k, react_species, react_order, react_count):
    """Rates of progress R[r] = k[r]·Π n^O at the state ``n`` (no Jacobian)."""
    n_react = k.shape[0]
    progress = np.empty(n_react)
    for r in range(n_react):
        cnt = react_count[r]
        prod = 1.0
        for a in range(cnt):
            s = react_species[r, a]
            prod *= n[s] ** react_order[r, a]
        progress[r] = k[r] * prod
    return progress


@njit(cache=True)
def _rhs_and_jac(n, nu, k, react_species, react_order, react_count):
    """Return (f, J, R): f = ν@R(n) (S,), J = ∂f/∂n (S,S), R the progress (R,)."""
    n_species = nu.shape[0]
    n_react = nu.shape[1]
    progress = np.empty(n_react)
    deriv = np.zeros((n_species, n_react))  # ∂R̃[r]/∂n[j]
    for r in range(n_react):
        cnt = react_count[r]
        prod = 1.0
        for a in range(cnt):
            s = react_species[r, a]
            prod *= n[s] ** react_order[r, a]
        progress[r] = k[r] * prod
        for a in range(cnt):
            j = react_species[r, a]
            oj = react_order[r, a]
            others = 1.0
            for b in range(cnt):
                if b != a:
                    s = react_species[r, b]
                    others *= n[s] ** react_order[r, b]
            deriv[j, r] = oj * (n[j] ** (oj - 1)) * others

    f = np.zeros(n_species)
    for i in range(n_species):
        acc = 0.0
        for r in range(n_react):
            acc += nu[i, r] * progress[r]
        f[i] = acc

    jac = np.zeros((n_species, n_species))
    for i in range(n_species):
        for r in range(n_react):
            w = nu[i, r] * k[r]
            if w != 0.0:
                for j in range(n_species):
                    d = deriv[j, r]
                    if d != 0.0:
                        jac[i, j] += w * d
    return f, jac, progress


@njit(cache=True)
def _linearly_implicit_trial(n, nu, k, rs, ro, rc, dt, eth):
    """One positivity-damped Rosenbrock-Euler trial and its inelastic energy cost."""
    f, jac, prog0 = _rhs_and_jac(n, nu, k, rs, ro, rc)
    amat = np.eye(n.shape[0]) - dt * jac
    delta = np.linalg.solve(np.ascontiguousarray(amat), dt * f)
    # One scalar damping factor retains the stoichiometric direction of the update. Per-species
    # clipping creates or destroys atoms and charge whenever an implicit step undershoots.
    alpha = 1.0
    for i in range(n.shape[0]):
        if delta[i] < 0.0 and n[i] + delta[i] < 0.0:
            alpha = min(alpha, n[i] / (-delta[i]))
    candidate = n + alpha * delta
    prog1 = _progress_only(candidate, k, rs, ro, rc)
    rate_cost = 0.0
    for r in range(k.shape[0]):
        if eth[r] > 0.0:
            rate_cost += eth[r] * 0.5 * (prog0[r] + prog1[r])
    return candidate, dt * rate_cost


@njit(cache=True)
def _augmented_linearly_implicit_trial(
    n, energy_ev_m3, nu, k, rs, ro, rc, dt, eth, electron_col, eps_floor,
):
    """One accepted Rosenbrock stage over ``[species, electron energy]``.

    The reaction extent is the same linearized stage that produces ``delta``:
    ``delta = nu @ extent``.  Its threshold debit is therefore not an endpoint quadrature
    estimate.  One damping factor scales the complete species update and energy debit.
    """
    f, jac, prog0 = _rhs_and_jac(n, nu, k, rs, ro, rc)
    amat = np.eye(n.shape[0]) - dt * jac
    delta = np.linalg.solve(np.ascontiguousarray(amat), dt * f)

    # Rosenbrock stage progress R(n) + dR/dn·delta.  This is the reaction-space stage whose
    # stoichiometric projection is exactly the solved species delta.
    stage = prog0.copy()
    for r in range(k.shape[0]):
        for a in range(rc[r]):
            j = rs[r, a]
            order_j = ro[r, a]
            others = 1.0
            for b in range(rc[r]):
                if b != a:
                    s = rs[r, b]
                    others *= n[s] ** ro[r, b]
            deriv = k[r] * order_j * n[j] ** (order_j - 1) * others
            stage[r] += deriv * delta[j]

    raw_cost = 0.0
    for r in range(k.shape[0]):
        raw_cost += eth[r] * dt * stage[r]
    if raw_cost < 0.0:
        # A negative linearized extent is outside mass-action admissibility.  It must not act
        # as an electron-energy source; the species positivity damping below remains visible.
        raw_cost = 0.0

    alpha = 1.0
    for i in range(n.shape[0]):
        if delta[i] < 0.0 and n[i] + delta[i] < 0.0:
            alpha = min(alpha, n[i] / (-delta[i]))

    # A costly stage cannot cross the post-reaction coefficient-table floor.  Threshold-free
    # electron creation (Penning) remains admissible at the floor; the resulting floor energy
    # is a separately ledgered numerical bound adjustment, not free reaction energy.
    if raw_cost > 0.0:
        slack = energy_ev_m3 - eps_floor * max(n[electron_col], 0.0)
        slope = -raw_cost - eps_floor * delta[electron_col]
        if slope < 0.0 and slack + alpha * slope < 0.0:
            alpha = min(alpha, max(slack, 0.0) / (-slope))
    if raw_cost > 0.0 and energy_ev_m3 - alpha * raw_cost < 0.0:
        alpha = min(alpha, energy_ev_m3 / raw_cost)

    alpha = min(max(alpha, 0.0), 1.0)
    candidate = n + alpha * delta
    return candidate, alpha * raw_cost, alpha


@njit(cache=True)
def _implicit_euler_substeps(n, nu, k, rs, ro, rc, dt, n_sub, eth, budget_ev):
    """Advance one cell by ``n_sub`` linearly-implicit Euler substeps of size ``dt`` (in place).

    Returns the cell's **inelastic energy spend** ``Σ_r ε_th,r·∫R_r dt`` (eV/m³) over the whole
    call, accumulated per substep with the trapezoidal rule on the rates of progress —
    ``eth`` is the per-reaction electron energy cost (eV; zero rows are skipped). This is the
    bookkeeping returned by the fixed-coefficient integration path. Source-coupled plasma runs
    use :func:`_energy_coupled_substeps`, which applies this debit inside each substep.
    """
    cost_ev = 0.0
    k_run = k.copy()
    for _ in range(n_sub):
        # Once the accumulated inelastic spend reaches the cell budget, energy-consuming
        # channels stop so the accepted chemistry cannot overspend its supplied energy.
        remaining = max(budget_ev - cost_ev, 0.0)
        if remaining <= 0.0:
            for r in range(k_run.shape[0]):
                if eth[r] > 0.0:
                    k_run[r] = 0.0
        candidate, step_cost = _linearly_implicit_trial(n, nu, k_run, rs, ro, rc, dt, eth)
        if step_cost > remaining and np.isfinite(remaining):
            # Scale only energy-consuming channels and resolve the same implicit step from
            # its original state. Bisection makes the per-cell licence effective *within*
            # the substep instead of allowing the last step to overspend it arbitrarily.
            lo = 0.0
            hi = 1.0
            best = n.copy()
            best_cost = 0.0
            k_trial = k_run.copy()
            for _bisect in range(40):
                mid = 0.5 * (lo + hi)
                for r in range(k_run.shape[0]):
                    k_trial[r] = k_run[r] * mid if eth[r] > 0.0 else k_run[r]
                trial, trial_cost = _linearly_implicit_trial(
                    n, nu, k_trial, rs, ro, rc, dt, eth
                )
                if trial_cost <= remaining:
                    lo = mid
                    best = trial
                    best_cost = trial_cost
                else:
                    hi = mid
            candidate = best
            step_cost = best_cost
        n[:] = candidate
        cost_ev += step_cost
    return cost_ev


@njit(parallel=True, cache=True)
def _advance_cells(densities, nu, k, rs, ro, rc, dt, n_sub, eth, budget_ev, cost_out):
    """Advance every cell's density vector in parallel (rows are independent)."""
    n_cells = densities.shape[0]
    for c in prange(n_cells):
        cost_out[c] = _implicit_euler_substeps(
            densities[c], nu, k, rs, ro, rc, dt, n_sub, eth, budget_ev[c]
        )


@njit(parallel=True, cache=True)
def _advance_cells_percell_k(densities, nu, k_percell, rs, ro, rc, dt, n_sub, eth, budget_ev,
                             cost_out):
    """Advance every cell with its **own** rate-coefficient row ``k_percell[c]`` (LMEA)."""
    n_cells = densities.shape[0]
    for c in prange(n_cells):
        cost_out[c] = _implicit_euler_substeps(
            densities[c], nu, k_percell[c], rs, ro, rc, dt, n_sub, eth, budget_ev[c]
        )


@njit(cache=True)
def _rates_at_energy(
    eps_ev, base_k, dynamic_kind, rate_lut, key_grid, te_lut, analytic_params,
):
    """Return one reaction-coefficient row evaluated at the live mean energy."""
    k = base_k.copy()
    nk = key_grid.shape[0]
    if nk == 1:
        jl = 0
        w = 0.0
    else:
        j = np.searchsorted(key_grid, eps_ev)
        if j <= 0:
            jl = 0
            w = 0.0
        elif j >= nk:
            jl = nk - 2
            w = 1.0
        else:
            jl = j - 1
            w = (eps_ev - key_grid[jl]) / (key_grid[jl + 1] - key_grid[jl])
    if nk == 1:
        te = te_lut[0]
    else:
        te = (1.0 - w) * te_lut[jl] + w * te_lut[jl + 1]
    for r in range(k.shape[0]):
        kind = dynamic_kind[r]
        if kind == 1:
            if nk == 1:
                k[r] = rate_lut[r, 0]
            else:
                k[r] = (1.0 - w) * rate_lut[r, jl] + w * rate_lut[r, jl + 1]
        elif kind == 2:
            k0 = analytic_params[r, 0]
            power = analytic_params[r, 1]
            e0 = analytic_params[r, 2]
            k[r] = k0 * te**power * np.exp(-e0 / te) if te > 0.0 else 0.0
    return k


@njit(cache=True)
def _energy_coupled_substeps(
    n, energy_ev_m3, nu, base_k, dynamic_kind, rate_lut, key_grid, te_lut,
    analytic_params, rs, ro, rc, dt, n_sub, eth, electron_col, eps_floor, eps_ceiling,
    source_ev_s,
):
    """Co-advance one cell's chemistry and electron energy over all implicit substeps.

    The reaction state and its threshold-energy integral are accepted together.  Costly
    channels are reduced only when that same trial would take the post-reaction electron
    population below the resolved energy floor.  Any energy introduced or discarded solely
    by the numerical energy bounds is returned separately instead of being hidden in the
    reaction ledger.

    ``source_ev_s`` is the cell's local power density (eV m^-3 s^-1, field work minus
    elastic exchange, delivered by the orchestrator) added at the top of every substep so
    heating and stiff inelastic cooling coexist inside the stage — the stage then relaxes
    to the true balance point instead of draining its stored buffer to the floor.  A
    negative source (delivery reconciliation) is admissible; the floor clamp absorbs and
    ledgers any overdraw.
    """
    cost_ev = 0.0
    bound_adjustment_ev = 0.0
    maximum_jump_ev = 0.0
    limiter_activity = 0
    energy = energy_ev_m3
    for _ in range(n_sub):
        energy += dt * source_ev_s
        ne_pre = max(n[electron_col], 0.0)
        lo_pre = eps_floor * ne_pre
        hi_pre = eps_ceiling * ne_pre
        bounded_pre = min(max(energy, lo_pre), hi_pre)
        bound_adjustment_ev += bounded_pre - energy
        if bounded_pre != energy:
            limiter_activity += 1
        energy = bounded_pre
        eps = energy / ne_pre if ne_pre > 0.0 else eps_floor
        eps = min(max(eps, eps_floor), eps_ceiling)
        k_run = _rates_at_energy(
            eps, base_k, dynamic_kind, rate_lut, key_grid, te_lut, analytic_params,
        )
        candidate, step_cost, alpha = _augmented_linearly_implicit_trial(
            n, energy, nu, k_run, rs, ro, rc, dt, eth, electron_col, eps_floor,
        )
        if alpha <= 1.0e-14 and energy <= eps_floor * ne_pre * (1.0 + 1.0e-14):
            # At the coefficient floor the costly part has zero admissible extent, but a
            # threshold-free Penning/photoelectron source is still physical.  Re-evaluate the
            # same augmented stage with costly rows disabled; no bisection or independent debit
            # is involved, and any extra floor energy remains a visible bound adjustment.
            k_zero = k_run.copy()
            for r in range(k_zero.shape[0]):
                if eth[r] > 0.0:
                    k_zero[r] = 0.0
            candidate, step_cost, alpha = _augmented_linearly_implicit_trial(
                n, energy, nu, k_zero, rs, ro, rc, dt, eth, electron_col, eps_floor,
            )
        if alpha < 1.0:
            limiter_activity += 1

        raw_energy = max(energy - step_cost, 0.0)
        ne_post = max(candidate[electron_col], 0.0)
        lo_post = eps_floor * ne_post
        hi_post = eps_ceiling * ne_post
        bounded_post = min(max(raw_energy, lo_post), hi_post)
        bound_adjustment_ev += bounded_post - raw_energy
        if bounded_post != raw_energy:
            limiter_activity += 1
        jump = abs(bounded_post - energy) / max(ne_post, 1.0)
        if jump > maximum_jump_ev:
            maximum_jump_ev = jump
        n[:] = candidate
        energy = bounded_post
        cost_ev += step_cost
    return energy, cost_ev, bound_adjustment_ev, maximum_jump_ev, limiter_activity


@njit(parallel=True, cache=True)
def _advance_cells_energy_coupled(
    densities, energy_ev_m3, nu, base_k, dynamic_kind, rate_lut, key_grid, te_lut,
    analytic_params, rs, ro, rc, dt, n_sub, eth, electron_col, eps_floor, eps_ceiling,
    source_ev_s, cost_out, adjustment_out, jump_out, limiter_out,
):
    """Parallel per-cell chemistry/electron-energy advance."""
    for c in prange(densities.shape[0]):
        energy, cost, adjustment, jump, limiter = _energy_coupled_substeps(
            densities[c], energy_ev_m3[c], nu, base_k, dynamic_kind, rate_lut,
            key_grid, te_lut, analytic_params, rs, ro, rc, dt, n_sub, eth,
            electron_col, eps_floor, eps_ceiling, source_ev_s[c],
        )
        energy_ev_m3[c] = energy
        cost_out[c] = cost
        adjustment_out[c] = adjustment
        jump_out[c] = jump
        limiter_out[c] = limiter


class FluidIntegrator:
    """Advance per-cell neutral densities with the numba implicit step (numpy fallback)."""

    def __init__(self, chemistry: FluidChemistry) -> None:
        self.chem = chemistry
        # Per-cell inelastic energy spend of the LAST advance (eV/m³): Σ_r ε_th,r·∫R_r dt.
        # Zero for networks without per-reaction thresholds (e.g. neutral-only chemistry).
        self.last_energy_cost_ev_m3: NDArray[np.float64] | None = None
        # Final co-integrated electron energy density (J/m3) and explicit energy introduced
        # or discarded by coefficient-table bounds.  Both refer to the LAST coupled advance.
        self.last_electron_energy_j_m3: NDArray[np.float64] | None = None
        self.last_energy_bound_adjustment_ev_m3: NDArray[np.float64] | None = None
        self.last_maximum_energy_jump_ev: float = 0.0
        self.last_limiter_activity: int = 0
        self._energy_pack_cache: tuple[object, tuple] | None = None

    def _thresholds(self) -> NDArray[np.float64]:
        eth = getattr(self.chem, "threshold_ev", None)
        if eth is None:
            return np.zeros(self.chem.n_reactions, dtype=np.float64)
        return np.ascontiguousarray(eth, dtype=np.float64)

    def advance(
        self,
        densities: NDArray[np.float64],
        interval: float,
        *,
        substeps: int = 64,
        k_percell: NDArray[np.float64] | None = None,
        energy_budget_ev_m3: NDArray[np.float64] | None = None,
    ) -> NDArray[np.float64]:
        """Advance ``densities`` (n_cells, n_species) by ``interval`` in ``substeps`` steps.

        In place; also returns the array. ``substeps`` (the internal implicit-Euler
        count over the subsample ``interval``) is caller-configurable, not hardcoded.
        ``k_percell`` (n_cells, n_reactions) supplies **per-cell** rate coefficients (the LMEA
        path — electron-impact rates evaluated at each cell's mean energy); when ``None`` the
        chemistry's shared ``k`` array is used for every cell.

        Side effect: ``last_energy_cost_ev_m3`` (n_cells,) records each cell's inelastic
        electron-energy spend over the call (per-reaction ``threshold_ev`` weighting,
        trapezoidal in time) — the operator split debits the electron energy density with it.
        """
        if interval < 0.0:
            raise ValueError("chemistry interval must be non-negative")
        if substeps <= 0:
            raise ValueError("chemistry substeps must be positive")
        d = np.ascontiguousarray(densities, dtype=np.float64)
        if d.ndim != 2 or d.shape[1] != self.chem.n_species:
            raise ValueError("densities must have shape (n_cells, n_species)")
        dt = interval / substeps
        chem = self.chem
        eth = self._thresholds()
        cost = np.zeros(d.shape[0], dtype=np.float64)
        budget = (np.full(d.shape[0], np.inf) if energy_budget_ev_m3 is None
                  else np.ascontiguousarray(energy_budget_ev_m3, dtype=np.float64))
        if budget.shape != (d.shape[0],) or np.any(np.isnan(budget)) or np.any(budget < 0.0):
            raise ValueError("energy_budget_ev_m3 must be one non-negative value per cell")
        percell = k_percell is not None
        if percell:
            k_arr = np.ascontiguousarray(k_percell, dtype=np.float64)
            if k_arr.shape != (d.shape[0], chem.n_reactions):
                raise ValueError("k_percell must have shape (n_cells, n_reactions)")
        if HAS_NUMBA:
            if percell:
                _advance_cells_percell_k(
                    d, chem.nu, k_arr, chem.react_species, chem.react_order,
                    chem.react_count, dt, substeps, eth, budget, cost,
                )
            else:
                _advance_cells(
                    d, chem.nu, chem.k, chem.react_species, chem.react_order,
                    chem.react_count, dt, substeps, eth, budget, cost,
                )
        else:  # pragma: no cover - the same reference kernel works without numba decoration
            for c in range(d.shape[0]):
                k_cell = k_arr[c] if percell else chem.k
                cost[c] = _implicit_euler_substeps(
                    d[c], chem.nu, k_cell, chem.react_species, chem.react_order,
                    chem.react_count, dt, substeps, eth, budget[c],
                )
        self.last_energy_cost_ev_m3 = cost
        if densities is not d:
            densities[...] = d
        return densities

    def advance_energy_coupled(
        self,
        densities: NDArray[np.float64],
        interval: float,
        *,
        electron_energy_j_m3: NDArray[np.float64],
        electron_col: int,
        coeffs,
        energy_bounds_ev: tuple[float, float],
        substeps: int = 64,
        power_source_j_m3_s: NDArray[np.float64] | None = None,
    ) -> NDArray[np.float64]:
        """Co-integrate reactions and their exact threshold-energy debit.

        Electron-impact and analytic-T_e coefficients are refreshed from the live mean
        energy before every internal implicit substep.  The accepted reaction trial and its
        threshold-energy cost are applied together, so chemistry can neither progress without
        paying nor pay for progress that was throttled away.

        ``power_source_j_m3_s`` (per cell, W m^-3, may be negative) is applied inside every
        substep before the trial, so field heating and stiff inelastic cooling share the
        same balance.  The per-cell identity extends to
        ``E_final + E_cost - E_bound_adjustment - source*interval = E_initial``.
        """
        if interval < 0.0:
            raise ValueError("chemistry interval must be non-negative")
        if substeps <= 0:
            raise ValueError("chemistry substeps must be positive")
        if coeffs.key_kind != "mean_energy_eV":
            raise ValueError("energy-coupled chemistry requires a mean-energy-keyed table")
        lo_ev, hi_ev = map(float, energy_bounds_ev)
        if not (0.0 < lo_ev < hi_ev):
            raise ValueError("energy_bounds_ev must be finite, positive, and increasing")

        chem = self.chem
        d = np.ascontiguousarray(densities, dtype=np.float64)
        if d.ndim != 2 or d.shape[1] != chem.n_species:
            raise ValueError("densities must have shape (n_cells, n_species)")
        if not 0 <= electron_col < chem.n_species:
            raise ValueError("electron_col is outside the chemistry species vector")
        energy_j = np.ascontiguousarray(electron_energy_j_m3, dtype=np.float64)
        if energy_j.shape != (d.shape[0],) or np.any(~np.isfinite(energy_j)) or np.any(energy_j < 0.0):
            raise ValueError("electron_energy_j_m3 must be one finite non-negative value per cell")
        if power_source_j_m3_s is None:
            source_ev_s = np.zeros(d.shape[0], dtype=np.float64)
        else:
            source_w = np.ascontiguousarray(power_source_j_m3_s, dtype=np.float64)
            if source_w.shape != (d.shape[0],) or np.any(~np.isfinite(source_w)):
                raise ValueError("power_source_j_m3_s must be one finite value per cell")
            source_ev_s = source_w / E_CHARGE

        if self._energy_pack_cache is not None and self._energy_pack_cache[0] is coeffs:
            key, te_lut, rate_lut, dynamic_kind, analytic = self._energy_pack_cache[1]
        else:
            key = np.ascontiguousarray(coeffs.key_grid, dtype=np.float64)
            te_lut = np.ascontiguousarray(coeffs.electron_temp_ev, dtype=np.float64)
            rate_lut = np.zeros((chem.n_reactions, key.size), dtype=np.float64)
            dynamic_kind = np.zeros(chem.n_reactions, dtype=np.int64)
            analytic = np.zeros((chem.n_reactions, 3), dtype=np.float64)
            if chem.dynamic_mask is not None:
                for r in np.flatnonzero(chem.dynamic_mask):
                    if r in chem.te_analytic:
                        dynamic_kind[r] = 2
                        analytic[r, :] = chem.te_analytic[r]
                    else:
                        dynamic_kind[r] = 1
                        rate_lut[r, :] = np.asarray(
                            coeffs.rate_coeff.get(chem.reaction_ids[r], np.zeros_like(key)),
                            dtype=np.float64,
                        )
            self._energy_pack_cache = (
                coeffs, (key, te_lut, rate_lut, dynamic_kind, analytic),
            )

        energy_ev = energy_j / E_CHARGE
        cost = np.zeros(d.shape[0], dtype=np.float64)
        adjustment = np.zeros(d.shape[0], dtype=np.float64)
        jump = np.zeros(d.shape[0], dtype=np.float64)
        limiter = np.zeros(d.shape[0], dtype=np.int64)
        dt = interval / substeps
        _advance_cells_energy_coupled(
            d, energy_ev, chem.nu, chem.k, dynamic_kind, rate_lut, key, te_lut,
            analytic, chem.react_species, chem.react_order, chem.react_count, dt,
            substeps, self._thresholds(), electron_col, lo_ev, hi_ev, source_ev_s,
            cost, adjustment, jump, limiter,
        )
        self.last_energy_cost_ev_m3 = cost
        self.last_electron_energy_j_m3 = energy_ev * E_CHARGE
        self.last_energy_bound_adjustment_ev_m3 = adjustment
        self.last_maximum_energy_jump_ev = float(np.max(jump)) if jump.size else 0.0
        self.last_limiter_activity = int(np.sum(limiter))
        if densities is not d:
            densities[...] = d
        return densities

    def warmup(self) -> None:
        """JIT-compile the kernels on a 1-cell dummy so the first real step isn't penalized."""
        if not HAS_NUMBA:
            return
        dummy = np.zeros((1, self.chem.n_species), dtype=np.float64)
        self.advance(dummy, 1e-18, substeps=1)
