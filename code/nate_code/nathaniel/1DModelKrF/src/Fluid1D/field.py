"""1D electrostatic field: discharge-circuit BCs + Poisson solve + field sampling.

The axial electric field over the 1D cavity is obtained from a Poisson solve
``-d²φ/dx² = ρ/ε₀`` with **Dirichlet** electrode potentials as boundary conditions. Those
electrode potentials are supplied each step by an in-loop **discharge-circuit** model
(:class:`DischargeCircuit`). An optional drive exposing ``at_batch`` and ``at_faces`` is
superposed on the Poisson field, and the fluid samples the result as a per-cell vec-3.

Two circuit backends are available:

- :class:`StubCircuit` — fixed/prescribed electrode potentials (a stand-in).
- :class:`MpcCircuit` — the real magnetic-pulse-compression (MPC) driver from Chen et al.,
  Phys. Plasmas 32, 113502 (2025): a 7-state ODE (energy-storage C₀ → magnetic switch L_m →
  peaking C_p → stray L_s → gap, with an R_f–L_f–C_f snubber across the electrodes). It is
  gas-invariant; the plasma enters only through the spatially-averaged conduction current it
  draws. The gap capacitance ``C_g = ε₀A/d`` sits at the cathode node so its voltage V_b is a
  proper ODE state (no DAE), handed back to the Poisson solve as the cathode Dirichlet value.

:class:`Poisson1DSolver` provides the plain Dirichlet solve plus a **semi-implicit
σ-corrected** variant (:meth:`~Poisson1DSolver.solve_semi_implicit`): the ohmic conduction
current responds to the *new* field, which removes the ε₀/σ dielectric-relaxation step limit
that makes the staggered-explicit coupling blow up past breakdown.
:class:`SemiImplicitField` bundles the solver with the electrode potentials so the plasma
model can re-solve the field inside its own substep loop.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from particles.constants import E_CHARGE, EPS0



class DischargeCircuit(ABC):
    """External discharge circuit: advances its ODE state and returns electrode potentials."""

    @abstractmethod
    def advance(
        self,
        dt: float,
        *,
        discharge_current_a: float = 0.0,
        gap_conductance_s: float = 0.0,
    ) -> tuple[float, float]:
        """Advance the circuit one ``dt``; return ``(phi_left, phi_right)`` electrode V.

        ``phi_left`` is the anode potential (grounded, x=0) and ``phi_right`` the cathode
        potential (powered, x=d) — the Dirichlet BCs for the 1D Poisson solve.
        ``discharge_current_a`` is the plasma **conduction** current A·⟨J_cond⟩ (A) drawn at
        the cathode node; the gap **displacement** current is handled internally by the
        circuit's own gap capacitance (see :class:`MpcCircuit`), so passing conduction-only
        avoids double-counting (Chen et al. §VI gap-coupling identity).
        ``gap_conductance_s`` is the optional differential conductance G = ∂I_cond/∂V_b that
        lets an implementation linearise the load about the current state; backends free to
        ignore it (a fixed-potential circuit has no node to load).
        """


class StubCircuit(DischargeCircuit):
    """Fixed / prescribed electrode potentials (no ODE; ignores the discharge current).

    ``waveform`` optionally supplies the right electrode as ``phi_right(t)``.
    """

    def __init__(
        self,
        phi_left: float = 0.0,
        phi_right: float = -20.0e3,
        *,
        waveform: Callable[[float], float] | None = None,
    ) -> None:
        self.phi_left = phi_left
        self.phi_right = phi_right
        self.waveform = waveform
        self._t = 0.0

    def advance(
        self,
        dt: float,
        *,
        discharge_current_a: float = 0.0,
        gap_conductance_s: float = 0.0,
    ) -> tuple[float, float]:
        self._t += dt  # prescribed potentials: the plasma load cannot move them
        if self.waveform is not None:
            return self.phi_left, float(self.waveform(self._t))
        return self.phi_left, self.phi_right


# State-vector indices for the 7-state MPC circuit ODE.
_V_C0, _I_M, _V_P, _I_S, _V_B, _I_F, _V_CF = range(7)


class MpcCircuit(DischargeCircuit):
    """Magnetic-pulse-compression discharge driver — 7-state ODE (Chen et al. 2025, §V).

    Topology (all node potentials vs. the grounded anode):

        C₀ ──L_m──● (V_p, C_p) ──L_s──● (V_b, cathode: C_g ∥ plasma ∥ R_f-L_f-C_f snubber)

    C₀ is pre-charged to ``v0_v`` (negative); energy resonantly transfers C₀→C_p (phase 1)
    then C_p→gap (phase 2, on breakdown). The gap capacitance ``C_g = ε₀·A/d`` (+ parasitic
    ``c_head_f``) at the cathode node makes V_b a proper ODE state — node-B KCL is an honest
    ODE, and the plasma feeds back only the **conduction** current A·⟨J_cond⟩ (displacement is
    C_g·V̇_b), per the §VI gap-coupling identity. Gas-invariant: geometry (``gap_m``,
    ``area_m2``) + the conduction current are the only plasma-side inputs.

    State ``y = [V_C0, I_m, V_p, I_s, V_b, I_f, V_Cf]`` (SI: V, A). The frozen-load linear ODE is
    advanced exactly with an augmented matrix exponential.
    """

    def __init__(
        self,
        *,
        gap_m: float,
        area_m2: float,
        c0_f: float = 10.0e-9,
        lm_h: float = 250.0e-9,
        cp_f: float = 10.0e-9,
        ls_h: float = 5.0e-9,
        rf_ohm: float = 70.71,
        lf_h: float = 100.0e-9,
        cf_f: float = 10.0e-12,
        v0_v: float = -18.0e3,
        r1_ohm: float = 0.02,
        rs_ohm: float = 0.02,
        c_head_f: float = 20.0e-12,
    ) -> None:
        self.c0 = c0_f
        self.lm = lm_h
        self.cp = cp_f
        self.ls = ls_h
        self.rf = rf_ohm
        self.lf = lf_h
        self.cf = cf_f
        self.r1 = r1_ohm
        self.rs = rs_ohm
        # Gap capacitance at the cathode node: physical ε₀A/d plus lumped head/feedthrough.
        self.cg_eff = EPS0 * area_m2 / gap_m + c_head_f
        self._open_gap_transition_cache: dict[float, NDArray[np.float64]] = {}
        # C₀ pre-charged to v0; everything else at rest.
        self.y = np.zeros(7, dtype=np.float64)
        self.y[_V_C0] = v0_v

    def advance(
        self,
        dt: float,
        *,
        discharge_current_a: float = 0.0,
        gap_conductance_s: float = 0.0,
    ) -> tuple[float, float]:
        """Advance the frozen-load 7-state ODE exactly; return (anode=0, cathode=V_b).

        ``discharge_current_a`` is A·⟨J_cond⟩, the plasma conduction current reported for the
        *current* state. ``gap_conductance_s`` is the gap's differential conductance G = ∂I/∂V_b;
        supplying it switches on the **semi-implicit gap linearisation** of Chen et al.'s circuit
        analysis §7,

            I_cond(V_b) ≈ G·V_b + (I_cond^n − G·V_b^n),

        so the load *responds* to V_b inside the ODE instead of acting as a frozen forcing term.
        Node B then relaxes on τ = C_g/G (a stable eigenvalue) rather than being driven open-loop:
        With G = 0 the reported current is a frozen forcing term. Positive G lets the load track
        the node-B voltage and gives the conducting gap a stable ``-G/C_g`` mode.
        """
        g = max(0.0, float(gap_conductance_s))  # a negative conductance would be anti-damping
        # Offset chosen so the linearisation reproduces the reported current at the current V_b.
        i_off = float(discharge_current_a) - g * float(self.y[_V_B])
        # The MPC network with the linearised gap load is a LINEAR constant-coefficient ODE
        # over the step, y' = A·y + b — so advance it **exactly** with the augmented matrix
        # exponential (top block of expm([[A, b], [0, 0]]·dt)). This is unconditionally
        # stable and exact for the stiff gap mode τ = C_g/G, which RK4 sub-stepping had to
        # resolve at ~G/4 substeps per step: at post-breakdown conductances (G spiking to
        # 1e4–1e6 S at the avalanche front) that was either the dominant wall-time cost of a
        # run (resolved honestly) or a V_b → ∞ instability (substeps capped).
        a = np.zeros((8, 8))
        a[_V_C0, _I_M] = -1.0 / self.c0
        a[_I_M, _V_C0] = 1.0 / self.lm
        a[_I_M, _V_P] = -1.0 / self.lm
        a[_I_M, _I_M] = -self.r1 / self.lm
        a[_V_P, _I_M] = 1.0 / self.cp
        a[_V_P, _I_S] = -1.0 / self.cp
        a[_I_S, _V_P] = 1.0 / self.ls
        a[_I_S, _V_B] = -1.0 / self.ls
        a[_I_S, _I_S] = -self.rs / self.ls
        a[_V_B, _I_S] = 1.0 / self.cg_eff
        a[_V_B, _I_F] = -1.0 / self.cg_eff
        a[_V_B, _V_B] = -g / self.cg_eff
        a[_I_F, _V_B] = 1.0 / self.lf
        a[_I_F, _I_F] = -self.rf / self.lf
        a[_I_F, _V_CF] = -1.0 / self.lf
        a[_V_CF, _I_F] = 1.0 / self.cf
        a[_V_B, 7] = -i_off / self.cg_eff  # inhomogeneous term (augmented column)
        if g == 0.0 and i_off == 0.0 and dt in self._open_gap_transition_cache:
            phi = self._open_gap_transition_cache[dt]
        else:
            from scipy.linalg import expm

            phi = expm(a * dt)
            if g == 0.0 and i_off == 0.0:
                # The production Lie split advances this identical open-gap MPC network at
                # one or two fixed timesteps thousands of times.  Cache the exact transition;
                # the state multiplication remains bit-identical to an uncached evaluation.
                self._open_gap_transition_cache[dt] = phi
        y_aug = np.append(self.y, 1.0)
        self.y = (phi @ y_aug)[:7]
        return 0.0, float(self.y[_V_B])

    # -- named accessors (diagnostics) ----------------------------------------
    @property
    def v_c0(self) -> float: return float(self.y[_V_C0])
    @property
    def v_p(self) -> float: return float(self.y[_V_P])
    @property
    def i_s(self) -> float: return float(self.y[_I_S])  # discharge-loop current
    @property
    def v_b(self) -> float: return float(self.y[_V_B])  # cathode / gap voltage
    @property
    def i_f(self) -> float: return float(self.y[_I_F])

    def pin_gap_voltage(self, v_b: float) -> None:
        """Adopt the kernel's co-advanced node-B voltage (gap-draw half of the split)."""
        self.y[_V_B] = float(v_b)

    def stored_energy_j(self) -> float:
        """Total capacitive + inductive energy currently in the circuit (J)."""
        v_c0, i_m, v_p, i_s, v_b, i_f, v_cf = self.y
        cap = 0.5 * (self.c0 * v_c0**2 + self.cp * v_p**2
                     + self.cg_eff * v_b**2 + self.cf * v_cf**2)
        ind = 0.5 * (self.lm * i_m**2 + self.ls * i_s**2 + self.lf * i_f**2)
        return float(cap + ind)


def conduction_current_0d(
    n_e_mean: float, mobility_e: float, area_m2: float, gap_m: float, v_gap: float
) -> float:
    """0D-closure gap conduction current A·⟨J_cond⟩ = G·V_gap, G = e·n_e·μ_e·A/d.

    This is the Chen et al. §VI bulk-conductance closure used by the stub plasma model.
    Drift-diffusion transport reports its load from the charged-species flux integral.
    """
    conductance = E_CHARGE * n_e_mean * mobility_e * area_m2 / gap_m
    return conductance * v_gap


def _thomas_general(
    lower: NDArray[np.float64],
    diag: NDArray[np.float64],
    upper: NDArray[np.float64],
    rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Solve a general tridiagonal system (Thomas algorithm, no pivoting).

    ``lower[r]`` multiplies ``x[r-1]`` (``lower[0]`` unused), ``diag[r]`` multiplies ``x[r]``,
    ``upper[r]`` multiplies ``x[r+1]`` (``upper[-1]`` unused). No pivoting is needed for the
    symmetric-positive-definite, diagonally-dominant systems assembled here.
    """
    m = rhs.shape[0]
    cp = np.empty(m)
    dp = np.empty(m)
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, m):
        denom = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denom
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denom
    x = np.empty(m)
    x[-1] = dp[-1]
    for i in range(m - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


class Poisson1DSolver:
    """Tridiagonal Dirichlet Poisson solve on the 1D mesh: ``-d²φ/dx² = ρ/ε₀``.

    Nodes sit at cell edges ``x_j = j·dx`` (``j = 0..N``); the electrode potentials pin the
    two end nodes (Dirichlet). The interior ``N-1`` nodes are solved with the constant
    ``[-1, 2, -1]`` operator. Cell-centred ``φ`` and ``E = -dφ/dx`` (both length ``N``) are
    returned, with the per-cell charge density mapped to interior nodes by neighbour average.
    """

    def __init__(self, n_cells: int, dx: float, *, eps0: float = EPS0) -> None:
        self.n_cells = int(n_cells)
        self.dx = dx
        self.eps0 = eps0

    def solve(
        self,
        charge_density_1d: NDArray[np.float64],
        phi_left: float,
        phi_right: float,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return ``(phi_cells, e_cells)`` (both length ``N``) for per-cell ρ and Dirichlet BCs."""
        return self.solve_semi_implicit(
            charge_density_1d, phi_left, phi_right,
            sigma_cells=np.zeros(self.n_cells), j_expl_faces=np.zeros(self.n_cells + 1),
            h=0.0,
        )

    def solve_semi_implicit(
        self,
        charge_density_1d: NDArray[np.float64],
        phi_left: float,
        phi_right: float,
        *,
        sigma_cells: NDArray[np.float64],
        j_expl_faces: NDArray[np.float64],
        h: float,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Semi-implicit σ-corrected Dirichlet solve (field prediction over a substep ``h``).

        Solves ``−d/dx[(ε₀ + h·σ(x))·dφ/dx] = ρ − h·dJ_expl/dx`` on the same node layout as
        :meth:`solve`. The ohmic conduction current ``σ·E`` responds to the *new* field —
        exactly the implicit treatment that removes the ``Δt ≲ ε₀/σ`` dielectric-relaxation
        stability limit of the staggered-explicit coupling — while everything else
        (``j_expl_faces``: the face explicit currents, walls zero — in practice the diffusive
        part) stays explicit and keeps the correct ambipolar limit as ``h·σ/ε₀ → ∞``.

        Cell-centred layout (2026-07-31): unknowns at cell centres, Dirichlet walls via dx/2
        half-cells, E native at the N+1 faces (``self.last_e_faces``); the returned
        ``e_cells`` is the face average. With ``σ = 0`` and ``J_expl = 0`` this reduces
        exactly to :meth:`solve`.
        """
        n = self.n_cells
        dx = self.dx
        rho = np.asarray(charge_density_1d, dtype=np.float64)
        a = self.eps0 + h * np.maximum(np.asarray(sigma_cells, dtype=np.float64), 0.0)  # (N,)
        jf = np.asarray(j_expl_faces, dtype=np.float64)  # (N+1,) face currents, walls 0
        # Cell-centred FV: unknowns φ_i at cell centres; interior face coefficient is the
        # arithmetic mean of the adjacent cells; wall faces sit dx/2 from the first/last
        # centre (factor 2). This layout sees the cell-alternating charge mode (the node-CV
        # gather 0.5(ρ_i+ρ_{i+1}) annihilated it EXACTLY — the measured odd–even
        # checkerboard: charge separation between adjacent cells drew no restoring field).
        af = np.empty(n + 1, dtype=np.float64)
        af[1:n] = 0.5 * (a[:-1] + a[1:])
        af[0] = a[0]
        af[n] = a[-1]
        diag = af[:-1] + af[1:]
        diag[0] = 2.0 * af[0] + af[1]
        diag[-1] = af[n - 1] + 2.0 * af[n]
        # _thomas_general convention: lower[r]·x[r−1] (lower[0] unused), upper[r]·x[r+1]
        # (upper[-1] unused) — pad to row length n.
        lower = np.empty(n, dtype=np.float64)
        upper = np.empty(n, dtype=np.float64)
        lower[0] = 0.0
        lower[1:] = -af[1:n]
        upper[-1] = 0.0
        upper[:-1] = -af[1:n]
        rhs = dx * dx * rho - h * dx * (jf[1:] - jf[:-1])
        rhs[0] += 2.0 * af[0] * phi_left
        rhs[-1] += 2.0 * af[n] * phi_right
        phi_cells = _thomas_general(lower, diag, upper, rhs)
        e_faces = np.empty(n + 1, dtype=np.float64)
        e_faces[1:n] = -(phi_cells[1:] - phi_cells[:-1]) / dx
        e_faces[0] = -(phi_cells[0] - phi_left) / (0.5 * dx)
        e_faces[n] = -(phi_right - phi_cells[-1]) / (0.5 * dx)
        self.last_e_faces = e_faces
        e_cells = 0.5 * (e_faces[:-1] + e_faces[1:])
        return phi_cells, e_cells


def sample_efield(
    poisson_e_1d: NDArray[np.float64],
    cell_centers: NDArray[np.float64],
    *,
    axis: int = 0,
    prescribed=None,
    t: float = 0.0,
) -> NDArray[np.float64]:
    """Combine the axial Poisson E with any prescribed field into a per-cell vec-3 ``(N, 3)``.

    ``poisson_e_1d`` is the axial component along ``axis``; ``prescribed`` is an optional
    object with ``at_batch(positions, t)`` sampled at cell centres and superposed.
    """
    n = poisson_e_1d.shape[0]
    efield = np.zeros((n, 3), dtype=np.float64)
    efield[:, axis] = poisson_e_1d
    if prescribed is not None:
        pos = np.zeros((n, 3), dtype=np.float64)
        pos[:, axis] = np.asarray(cell_centers, dtype=np.float64)
        efield += np.asarray(prescribed.at_batch(pos, t), dtype=np.float64)
    return efield


@dataclass
class SemiImplicitField:
    """Per-outer-step field handle: lets the plasma model re-solve E inside its substep loop.

    Bundles the Poisson solver with the electrode Dirichlet potentials — ``phi_right`` is the
    cathode ``V_b``, frozen over the outer step because the circuit owns node B — plus any
    prescribed drive (superposed *after* the solve, exactly as :func:`sample_efield` does; an
    external field is not part of the implicit gap balance). :meth:`solve` runs the
    semi-implicit σ-corrected solve for a substep ``h`` and returns the fresh axial field;
    the full vec-3 field is kept in ``last_efield`` for the orchestrator's diagnostics.
    """

    solver: Poisson1DSolver
    phi_left: float
    phi_right: float
    cell_centers: NDArray[np.float64] | None = None  # only needed when ``prescribed`` is set
    prescribed: object | None = None  # optional at_batch/at_faces field provider
    t: float = 0.0
    axis: int = 0
    last_efield: NDArray[np.float64] | None = None  # (N, 3), refreshed by every solve
    # Optional gap-node co-advance (finite C_g): the kernel sags the cathode Dirichlet value
    # as current is drawn, and reports the final V_b for the circuit to adopt (Lie split of
    # node B: circuit owns the network drive, the kernel owns the gap draw).
    gap_capacitance_f: float | None = None
    gap_drive_a: float = 0.0   # network drive i_s − i_f feeding node B (frozen per exchange)
    v_b_final: float | None = None

    def solve(
        self,
        charge_density_1d: NDArray[np.float64],
        sigma_cells: NDArray[np.float64],
        j_expl_faces: NDArray[np.float64],
        h: float,
    ) -> NDArray[np.float64]:
        """Fresh axial E (V/m, length N) from the semi-implicit solve for substep ``h``."""
        _phi, e1d = self.solver.solve_semi_implicit(
            charge_density_1d, self.phi_left, self.phi_right,
            sigma_cells=sigma_cells, j_expl_faces=j_expl_faces, h=h,
        )
        # Transport uses face-native E while heating/diagnostics use cell-centred E.  Apply
        # an external drive to both representations; otherwise prescribed fields heat the
        # electrons but produce no Scharfetter-Gummel drift or electrode flux.
        self.last_e_faces = np.asarray(self.solver.last_e_faces, dtype=np.float64).copy()
        if self.prescribed is not None:
            face_pos = np.zeros((self.solver.n_cells + 1, 3), dtype=np.float64)
            face_pos[:, self.axis] = np.arange(self.solver.n_cells + 1) * self.solver.dx
            prescribed_faces = np.asarray(
                self.prescribed.at_batch(face_pos, self.t), dtype=np.float64
            )
            self.last_e_faces += prescribed_faces[:, self.axis]
        centers = self.cell_centers
        if centers is None:
            centers = (np.arange(self.solver.n_cells) + 0.5) * self.solver.dx
        self.last_efield = sample_efield(
            e1d, centers, axis=self.axis, prescribed=self.prescribed, t=self.t
        )
        return self.last_efield[:, self.axis]
