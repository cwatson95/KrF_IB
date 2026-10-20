"""Two-term Boltzmann solver for the electron coefficients (Hagelaar & Pitchford 2005).

Solves the steady, spatially uniform **two-term** electron Boltzmann equation in energy
space (H&P Eq. 39) for the isotropic EEDF ``F₀(ε)`` at each reduced field E/N, and builds
the same mean-energy-keyed :class:`~Fluid1D.boltzmann.CoefficientTable` the LMEA transport
model consumes — replacing the Maxwellian closure of
:class:`~Fluid1D.eedf_solver.MaxwellianEedfSolver`, whose over-populated tail makes the
discharge break down at roughly half the field Chen et al. (2025) report.

Formulation (H&P Eqs. 39-44, steady state, no spatial gradients):

    d/dε [ W̃(ε)·F₀ − D̃(ε)·dF₀/dε ] = S̃(ε)

    W̃ = −γ ε² Σ_g (2mₑ/M_g)·x_g·σ_{m,g}(ε)                        (elastic cooling)
    D̃ = (γ/3)·(E/N)²·ε/σ̃_m(ε) + γ·(k_B T_g/e)·ε² Σ_g (2mₑ/M_g)·x_g·σ_{m,g}   (field + thermal)
    S̃ = Σ_j γ x_j [ (ε+u_j)σ_j(ε+u_j)F₀(ε+u_j) − ε σ_j(ε)F₀(ε) ]   (non-local inelastic)

with γ = √(2e/mₑ), σ̃_m the effective total momentum-transfer cross-section (mixture elastic
plus the inelastic process cross-sections), normalised ``∫F₀√ε dε = 1``. The energy-space
flux is discretised with the **same Scharfetter–Gummel kernel the plasma transport uses**
(:func:`Fluid1D.transport.sg_interior_flux` — verified identical to H&P's own exponential
scheme, their Eq. 45), and each inelastic process moves electrons **down** by its threshold
``u_j`` with an exactly number-conserving two-cell deposit (the discrete form of H&P's
P/Q operators, Eqs. 46-51).

Growth model (H&P §II.B, temporal): the net production ν̄ = N·(Σ k_ion·x − Σ k_att·x) enters
the balance as the distributed term ``(ν̄/N)·√ε·F₀`` and is found by fixed-point iteration
(F₀ → ν̄ → F₀, a handful of passes). This is **not optional** for an attaching mixture: at low
E/N the no-growth equation would describe only the survivors of attachment — an EEDF pinned
near ~2.4 eV in this Ne/F₂ mix — while the true (decaying) distribution relaxes to the gas
temperature. Ionization deposits the scattered primary at ε−u and the secondary at ε≈0
("one takes all" sharing).

Model approximation:
- σ̃_m adds every inelastic σ on top of the per-gas elastic/effective set; for gases whose
  base set is already "effective" (Ar, F₂ here) the inelastic part is counted twice — a
  few-percent overestimate of σ̃_m near thresholds.

Coefficients (H&P Eqs. 55-56 and the LMEA energy pair):

    μN  = −(γ/3)·∫ (ε/σ̃_m)·dF₀/dε dε          DN  = (γ/3)·∫ (ε/σ̃_m)·F₀ dε
    μ_εN = −(γ/3)·∫ (ε²/σ̃_m)·dF₀/dε dε / ⟨ε⟩   D_εN = (γ/3)·∫ (ε²/σ̃_m)·F₀ dε / ⟨ε⟩
    k_r = γ·∫ ε σ_r(ε)·F₀ dε                    ⟨ε⟩ = ∫ ε^{3/2}·F₀ dε

For a Maxwellian F₀ these reduce exactly to the :class:`MaxwellianEedfSolver` moments (with
μ_ε/μ = D_ε/D = 5/3) — one of the validation anchors in ``tests/test_two_term.py``.
"""

from __future__ import annotations

import numpy as np

from particles.constants import E_CHARGE, K_B, M_E, molecular_mass
from reactions.types import RateForm

from .boltzmann import CoefficientTable, _dynamic_reactions
from .eedf_solver import GAMMA_V_PER_SQRT_EV, MaxwellianEedfSolver, _atoms_of
from .transport import bernoulli


class TwoTermBoltzmannSolver(MaxwellianEedfSolver):
    """Real two-term EEDF solve per E/N (inherits the cross-section plumbing of the
    Maxwellian solver; only :meth:`solve` and the moments differ).

    ``solve`` sweeps the caller's ``en_grid_td`` (the values matter here, unlike the
    Maxwellian solver which is natively Tₑ-keyed), solves F₀ at each point, and returns a
    table keyed on the resulting mean energy — E/N is stored as ``reduced_field_td`` so the
    local-field ⟨ε⟩_eq(E/N) closure works natively.
    """

    def __init__(
        self,
        cross_sections_dir,
        *,
        d_eps_target_ev: float = 0.05,
        n_energy_cells_min: int = 400,
        n_energy_cells_max: int = 4000,
        energy_max_floor_ev: float = 60.0,
        energy_max_cap_ev: float = 400.0,
        extra_loss_process_ids: tuple[str, ...] = (),
        sigma_lut_points: int = 1024,
    ) -> None:
        super().__init__(
            cross_sections_dir,
            extra_loss_process_ids=extra_loss_process_ids,
            sigma_lut_points=sigma_lut_points,
        )
        # Δε must resolve the inelastic "wall" (F₀ plunges over ≲1 eV past the first
        # threshold; everything above it — the ionization integral — is exponentially
        # sensitive to that plunge). The sparse solve makes a few thousand cells cheap.
        self.d_eps_target_ev = float(d_eps_target_ev)
        self.n_energy_cells_min = int(n_energy_cells_min)
        self.n_energy_cells_max = int(n_energy_cells_max)
        self.energy_max_floor_ev = float(energy_max_floor_ev)
        self.energy_max_cap_ev = float(energy_max_cap_ev)

    @property
    def solver_id(self) -> str:
        from .eedf_solver import _cross_section_fingerprint

        return (f"two-term:v2:de{self.d_eps_target_ev:g}:cap{self.energy_max_cap_ev:g}"
                f":floor{self.energy_max_floor_ev:g}:n{self.n_energy_cells_min}-"
                f"{self.n_energy_cells_max}:nlut{self.sigma_lut_points}"
                f":loss{','.join(self.extra_loss_process_ids)}"
                f":xs{_cross_section_fingerprint(self.cross_sections_dir)}")

    def _n_cells(self, eps_max: float) -> int:
        return int(np.clip(round(eps_max / self.d_eps_target_ev),
                           self.n_energy_cells_min, self.n_energy_cells_max))

    # ------------------------------------------------------------------ solve
    def solve(self, mixture, en_grid_td, gas_temp_k=300.0, *, reactions=None) -> CoefficientTable:
        en_grid = np.asarray(en_grid_td, dtype=np.float64)
        # Extend the sweep to low E/N so the ⟨ε⟩-keyed table covers cooled (post-collapse)
        # cells: in a Ne buffer ⟨ε⟩(0.1 Td) is already ~2.4 eV, so without the extension the
        # LMEA clamp would pin every cold cell to that value.
        if en_grid.min() > 1.0e-2:
            low = np.geomspace(1.0e-2, en_grid.min(), 12, endpoint=False)
            en_grid = np.unique(np.concatenate([low, en_grid]))
        processes = self._collect_processes(mixture, reactions)

        # First point serial (pre-warms the σ caches), remainder threaded — the per-point
        # solves are independent and spend most of their time under GIL-releasing
        # scipy/LAPACK calls, so a thread pool gives a real (2-4x) build speedup.
        rows = [self._solve_one(mixture, float(en_grid[0]), gas_temp_k, processes)]
        if en_grid.size > 1:
            import os
            from concurrent.futures import ThreadPoolExecutor

            workers = min(8, max(1, (os.cpu_count() or 2) - 1))
            with ThreadPoolExecutor(max_workers=workers) as pool:
                rows += list(pool.map(
                    lambda en: self._solve_one(mixture, float(en), gas_temp_k, processes),
                    en_grid[1:],
                ))

        mean_energy = np.array([r["mean_energy"] for r in rows])
        order = np.argsort(mean_energy, kind="stable")
        # Strictly increasing key grid (⟨ε⟩(E/N) is physically monotone; enforce numerically).
        key = mean_energy[order]
        key = np.maximum.accumulate(key)
        key += np.arange(key.size) * 1e-12  # break exact ties without moving values

        def col(name):
            return np.array([rows[i][name] for i in order])

        rate_coeff = {}
        for rid in rows[0]["rates"]:
            rate_coeff[rid] = np.array([rows[i]["rates"][rid] for i in order])

        return CoefficientTable(
            key_kind="mean_energy_eV",
            key_grid=key,
            mobility_N=col("mob"),
            diffusion_N=col("dif"),
            drift_velocity=col("drift"),
            mean_energy_ev=key,
            electron_temp_ev=(2.0 / 3.0) * key,
            rate_coeff=rate_coeff,
            mixture=dict(mixture),
            gas_temp_k=float(gas_temp_k),
            energy_mobility_N=col("mob_e"),
            energy_diffusion_N=col("dif_e"),
            reduced_field_td=en_grid[order],
        )

    # ------------------------------------------------------------------ pieces
    def _collect_processes(self, mixture, reactions):
        """Inelastic processes as (id_or_None, x_target, threshold_eV, sigma, is_att, is_ion).

        ``id`` is the reaction id for channels that need a rate-coefficient output; loss-only
        channels (e.g. F₂ vibrational) carry ``None``. Attachment removes the electron
        (negative growth); ionization adds the secondary at ε ≈ 0 (positive growth); every
        channel shifts the scattered primary down by its threshold.
        """
        procs = []
        if reactions is not None:
            for rxn in _dynamic_reactions(reactions):
                if rxn.rate_form is not RateForm.ELECTRON_IMPACT:
                    continue
                key = getattr(rxn, "cross_section", None)
                sig = self._sigma(str(key)) if key else None
                if sig is None:
                    continue  # solve() caller already warned via the Maxwellian path
                # Elastic/effective datasets already form the mixture momentum-transfer
                # cross section.  Adding them again as an inelastic process double-counts
                # momentum transfer and invents a zero-threshold energy-shift channel.
                if sig.record.reaction_kind.value == "elastic":
                    continue
                target_frac = sum(
                    float(mixture.get(name, 0.0)) for _c, name in rxn.reactants
                )
                eth = float(sig.threshold_ev)
                n_e_out = sum(c for c, name in rxn.products if name == "e-")
                n_e_in = sum(c for c, name in rxn.reactants if name == "e-")
                is_ion = n_e_out > n_e_in
                is_att = n_e_out < n_e_in
                procs.append((rxn.reaction_id, target_frac, eth, sig, is_att, is_ion))
        for pid in self.extra_loss_process_ids:
            sig = self._sigma(pid)
            if sig is None:
                continue
            frac = float(mixture.get(sig.record.target, 0.0))
            if frac <= 0.0:
                continue
            procs.append((None, frac, float(sig.threshold_ev), sig, False, False))
        return procs

    def _solve_one(self, mixture, en_td, gas_temp_k, processes):
        """One E/N point: adapt the energy window to ⟨ε⟩, solve F₀, take the moments.

        The window must both resolve the bulk (≳20 cells per ⟨ε⟩) and, once the bulk is warm
        enough to feed the inelastic thresholds, reach past the highest threshold — but a cold
        (sub-threshold) point must NOT be forced to a threshold-wide window, or the thermal
        bulk would be unresolved.
        """
        en_vm2 = en_td * 1.0e-21  # E/N in V·m²
        u_max = max((p[2] for p in processes), default=0.0)
        eps_max = self.energy_max_floor_ev
        result = None
        for _pass in range(6):
            result = self._solve_f0(mixture, en_vm2, gas_temp_k, processes, eps_max)
            f0, eps = result["f0"], result["eps"]
            d_eps = eps[1] - eps[0]
            # Under-resolved thermal bulk: the population collapses into the first cells.
            if float(np.sum(f0[:3] * np.sqrt(eps[:3])) * d_eps) > 0.5:
                eps_max = max(eps_max / 20.0, 30.0 * K_B * gas_temp_k / E_CHARGE)
                continue
            eps1 = result["mean_energy"]
            new_max = 25.0 * eps1
            if u_max > 0.0 and eps1 > 0.5:
                new_max = max(new_max, u_max + 10.0 * eps1)
            new_max = float(np.clip(new_max, 30.0 * K_B * gas_temp_k / E_CHARGE,
                                    self.energy_max_cap_ev))
            if abs(new_max - eps_max) < 0.2 * eps_max:
                break
            eps_max = new_max
        return result

    def _solve_f0(self, mixture, en_vm2, gas_temp_k, processes, eps_max):
        from scipy.sparse import csc_matrix
        from scipy.sparse.linalg import splu

        gamma = GAMMA_V_PER_SQRT_EV
        n = self._n_cells(eps_max)
        d_eps = eps_max / n
        eps = (np.arange(n) + 0.5) * d_eps          # cell centres
        eps_f = np.arange(1, n) * d_eps             # interior faces

        # Elastic recoil weight Σ_g (2mₑ/M_g)·x_g·σ_{m,g}(ε) and the mixture σ_m, on both grids.
        recoil_c = np.zeros(n)
        recoil_f = np.zeros(n - 1)
        sigma_m_c = np.zeros(n)
        for gas, frac in mixture.items():
            if frac <= 0.0:
                continue
            sig = self._elastic_sigma_for(gas)
            if sig is None:
                continue
            try:
                mass = molecular_mass(*_atoms_of(gas))
            except (KeyError, ValueError):
                continue
            w = 2.0 * M_E / mass * float(frac)
            recoil_c += w * np.asarray(sig(eps), dtype=np.float64)
            recoil_f += w * np.asarray(sig(eps_f), dtype=np.float64)
            sigma_m_c += float(frac) * np.asarray(sig(eps), dtype=np.float64)
        # Effective total momentum cross-section: elastic + inelastic processes.
        sigma_tot_c = sigma_m_c.copy()
        for _rid, frac, _eth, sig, _att, _ion in processes:
            sigma_tot_c += frac * np.asarray(sig(eps), dtype=np.float64)
        sigma_tot_c = np.maximum(sigma_tot_c, 1.0e-24)
        sigma_tot_f = 0.5 * (sigma_tot_c[:-1] + sigma_tot_c[1:])

        kt_g_ev = K_B * gas_temp_k / E_CHARGE
        # Energy-space drift (cooling) and diffusion (field heating + thermal), at faces.
        w_face = -gamma * eps_f**2 * recoil_f
        d_face = (gamma / 3.0) * en_vm2**2 * eps_f / sigma_tot_f \
            + gamma * kt_g_ev * eps_f**2 * recoil_f
        d_face = np.maximum(d_face, 1.0e-300)

        # Assemble A·F₀ = 0 (flux divergence + inelastic operator) as a sparse matrix —
        # tridiagonal SG flux (the same exponential-fitting coefficients sg_interior_flux
        # applies; H&P's own Eq. 45 scheme) plus one off-band per inelastic process.
        idx = np.arange(n)
        pe = w_face * d_eps / d_face
        c_lo = (d_face / d_eps) * bernoulli(-pe)    # multiplies F₀_i     in J_{i+1/2}
        c_hi = (d_face / d_eps) * bernoulli(pe)     # multiplies F₀_{i+1} in J_{i+1/2}
        inv = 1.0 / d_eps
        rows = [idx[:-1], idx[:-1], idx[1:], idx[1:]]           # ±J_{i±1/2}/dε
        cols = [idx[:-1], idx[1:], idx[1:], idx[:-1]]
        vals = [inv * c_lo, -inv * c_hi, inv * c_hi, -inv * c_lo]

        growth_sign = []  # +1 ionization, −1 attachment, 0 conservative — for ν̄
        for _rid, frac, eth, sig, is_att, is_ion in processes:
            growth_sign.append(1.0 if is_ion else (-1.0 if is_att else 0.0))
            if frac <= 0.0:
                continue
            nu = gamma * frac * eps * np.asarray(sig(eps), dtype=np.float64)
            act = np.flatnonzero(nu > 0.0)
            if act.size == 0:
                continue
            rows.append(act); cols.append(act); vals.append(nu[act])          # loss
            if is_att:
                continue  # electron removed
            if is_ion:  # secondary at ε ≈ 0 ("one takes all" sharing)
                rows.append(np.zeros(act.size, dtype=np.int64)); cols.append(act)
                vals.append(-nu[act])
            # Scattered electron re-appears at ε_i − u_j (two-cell conservative deposit).
            target = eps[act] - eth
            low = target <= eps[0]
            if low.any():
                al = act[low]
                rows.append(np.zeros(al.size, dtype=np.int64)); cols.append(al)
                vals.append(-nu[al])
            ah = act[~low]
            if ah.size:
                f = target[~low] / d_eps - 0.5
                lo = np.clip(np.floor(f).astype(np.int64), 0, n - 2)
                w_hi = f - lo
                rows.append(lo); cols.append(ah); vals.append(-nu[ah] * (1.0 - w_hi))
                rows.append(lo + 1); cols.append(ah); vals.append(-nu[ah] * w_hi)

        rows_a = np.concatenate(rows)
        cols_a = np.concatenate(cols)
        vals_a = np.concatenate(vals)
        sqrt_eps = np.sqrt(eps)
        # Physical operator (all balance rows) — used by the growth eigen-iteration.
        a_phys = csc_matrix((vals_a, (rows_a, cols_a)), shape=(n, n))
        # Normalised system: the top-energy balance row (F₀ ~ 0 up there) is replaced by
        # the normalisation ∫F₀√ε dε = 1 — used for the final eigenvector refinement.
        keep = rows_a != n - 1
        rows_n = np.concatenate([rows_a[keep], np.full(n, n - 1, dtype=np.int64)])
        cols_n = np.concatenate([cols_a[keep], idx])
        vals_n = np.concatenate([vals_a[keep], sqrt_eps * d_eps])
        a_base = csc_matrix((vals_n, (rows_n, cols_n)), shape=(n, n))
        growth_diag = sqrt_eps.copy()
        growth_diag[n - 1] = 0.0  # never touch the normalisation row
        b = np.zeros(n)
        b[n - 1] = 1.0

        # Temporal-growth eigenvalue (H&P §II.B): the net production ν̄/N enters the balance
        # as the distributed term (ν̄/N)·√ε·F₀ — see the Perron-eigenpair iteration below.
        # Required in an attaching mixture — see the module docstring.
        state = {"col": np.ones(n), "first": True}

        def eval_g(nu_bar: float):
            """Production G(ν̄) and its F₀ — with iterative column scaling: scaling the
            unknowns by the previous iterate keeps them O(1) across the exponentially
            decaying tail, preserving componentwise accuracy through the factorisation."""
            f0 = None
            col = state["col"]
            for _ in range(3 if state["first"] else 1):
                a_it = a_base + csc_matrix(
                    (nu_bar * growth_diag, (idx, idx)), shape=(n, n)
                )
                a_it = a_it.multiply(col[None, :]).tocsr()
                scale = np.abs(a_it).max(axis=1).toarray().ravel()
                scale = np.where(scale > 0.0, scale, 1.0)
                y = splu(a_it.multiply(1.0 / scale[:, None]).tocsc()).solve(b / scale)
                f0 = np.maximum(y, 0.0) * col
                norm = float(np.sum(f0 * sqrt_eps) * d_eps)
                if norm > 0.0:
                    f0 /= norm
                col = np.clip(f0, f0.max() * 1e-14 + 1e-300, None)
            state["col"] = col
            state["first"] = False
            g = 0.0
            for (_rid, frac, _eth, sig, _att, _ion), s in zip(processes, growth_sign):
                if s == 0.0 or frac <= 0.0:
                    continue
                g += s * gamma * frac * float(
                    np.sum(eps * np.asarray(sig(eps), dtype=np.float64) * f0) * d_eps
                )
            return g, f0

        g0, f0 = eval_g(0.0)
        if abs(g0) > 1e-24:
            # The self-consistent (ν̄, F₀) is the **principal (Perron) eigenpair** of the
            # pencil A_phys·F = −ν̄·M·F with M = diag(√ε): F ≥ 0 and ν̄ the algebraically
            # largest growth rate. (Fixed-point/root-finding on the production map is NOT
            # reliable here — attachment heating makes it non-monotone near the
            # attachment↔ionization balance.) Shifted inverse-power iteration: with c > ν̄,
            # B = (A_phys + c·M)⁻¹·M is elementwise non-negative (M-matrix inverse), so
            # power iteration on B converges to the Perron mode; μ = 1/(c − ν̄).
            nu_bar = 0.0
            fvec = np.maximum(f0, 1e-300)
            c = 4.0 * abs(g0)
            for _round in range(4):
                lu = splu((a_phys + csc_matrix((c * sqrt_eps, (idx, idx)),
                                               shape=(n, n))).tocsc())
                mu_prev = None
                stalled = False
                for _it in range(150):
                    rhs = sqrt_eps * fvec
                    f_new = lu.solve(rhs)
                    if not np.all(np.isfinite(f_new)) or f_new.max() <= 0.0:
                        stalled = True
                        break
                    f_new = np.maximum(f_new, 0.0)
                    mu = float(np.sum(f_new * sqrt_eps) / max(np.sum(fvec * sqrt_eps), 1e-300))
                    fvec = f_new / f_new.max()
                    if mu_prev is not None and abs(mu - mu_prev) <= 1e-6 * abs(mu):
                        mu_prev = mu
                        break
                    mu_prev = mu
                if stalled or mu_prev is None or mu_prev <= 0.0:
                    c *= 8.0  # shift sat at/below ν̄ — push it up and retry
                    continue
                nu_est = c - 1.0 / mu_prev
                # Good separation: accept; poor separation (ν̄ ≈ c): retighten the shift.
                if nu_est < 0.75 * c or abs(nu_est) < 1e-24:
                    nu_bar = nu_est
                    break
                c = 4.0 * abs(nu_est) + 1e-24
            # Final eigenvector with the normalised, column-scaled solve at the converged ν̄
            # (the normalisation row resolves the near-singular system exactly at the root).
            _g, f0 = eval_g(nu_bar)

        out = self._moments(eps, d_eps, f0, sigma_tot_c, processes)
        out["drift"] = out["mob"] * en_vm2
        out["f0"] = f0
        out["eps"] = eps
        return out

    @staticmethod
    def _moments(eps, d_eps, f0, sigma_tot_c, processes) -> dict:
        """H&P Eqs. 55-56 moments of an EEDF on a uniform cell grid (validation anchor).

        For a Maxwellian ``f0`` these reduce to closed forms: D/μ = Tₑ exactly (any σ), and
        D_ε/D = μ_ε/μ = 4/3 for constant σ (5/3 for constant collision frequency).
        """
        gamma = GAMMA_V_PER_SQRT_EV
        mean_energy = float(np.sum(eps**1.5 * f0) * d_eps)
        df0 = np.gradient(f0, d_eps)
        mob = -(gamma / 3.0) * float(np.sum(eps / sigma_tot_c * df0) * d_eps)
        dif = (gamma / 3.0) * float(np.sum(eps / sigma_tot_c * f0) * d_eps)
        me_safe = max(mean_energy, 1e-30)
        mob_e = -(gamma / 3.0) * float(np.sum(eps**2 / sigma_tot_c * df0) * d_eps) / me_safe
        dif_e = (gamma / 3.0) * float(np.sum(eps**2 / sigma_tot_c * f0) * d_eps) / me_safe
        rates = {}
        for rid, _frac, _eth, sig, _att, _ion in processes:
            if rid is None:
                continue
            rates[rid] = gamma * float(np.sum(eps * np.asarray(sig(eps)) * f0) * d_eps)
        return dict(
            mean_energy=mean_energy, mob=max(mob, 0.0), dif=max(dif, 0.0),
            mob_e=max(mob_e, 0.0), dif_e=max(dif_e, 0.0), rates=rates,
        )
