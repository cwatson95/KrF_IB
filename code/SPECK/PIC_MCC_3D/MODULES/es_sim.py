"""The minimal electrostatic PIC loop that composes M1 and M2.

Per step, the classic cycle:

    deposit rho  ->  solve phi  ->  E = -grad phi  ->  gather  ->  push v
                 ->  move x     ->  boundaries     ->  COLLIDE (optional)

This exists so the physics gates are REAL: a Poisson solver that matches an
LU oracle proves arithmetic, but only a closed loop can show that deposit,
solve, gather and push are mutually consistent -- which is what the cold
plasma-oscillation gate (measured f/f_pe) actually tests.

Pass ``mcc=`` (a :class:`~PIC_MCC_3D.MODULES.mcc3d.MCC3D`) to make it a
PIC-*MCC*: electrons then collide off the neutral background once per field
step, ionization can create real ion macro-particles, and the per-process
tallies and energy ledger accumulate on the run.

WHAT THIS IS NOT.  No external circuit, no dielectric charging, no Coulomb
collisions, no ion-neutral collisions and no implicit mover: of missing
modules 4-9 (``../PIC_MCC_3D.txt`` section 2.1) only slices are in --
electrode planes take a value or a phi(t) waveform (M4 partial, via
``_bias``), and surface emission lives in ``emission3d`` (M5 partial: the
gamma process + a flux emitter, registered through ``wall_hooks`` /
``injectors``; the G5.1 Turner claim stays open).  The explicit
leapfrog below carries the hard stability requirements

    dx < lambda_D        (or the grid self-heats)
    omega_pe dt < 2      (or the mover is unstable)

and :meth:`ESSim.stability` reports both so a run cannot quietly violate them.

LEAPFROG BOOKKEEPING.  Velocities live at half-steps and positions at whole
steps.  :meth:`prepare` retards v by half a step so a run starts centred; the
kinetic and field energies are therefore offset by dt/2, which is why the
energy history here is a diagnostic and not a conservation claim (an
energy-conserving mover is M7).
"""

from __future__ import annotations

import numpy as np

from . import es_field as EF
from . import species3d as SP
from .poisson3d import PERIODIC, PERIODIC_XY_DIRICHLET_Z, Poisson3D
from .units import E_CHARGE as EF_E_CHARGE
from .units import EPS0, debye_length, plasma_frequency


class ESSim:
    """A 3-D electrostatic PIC run on one grid."""

    def __init__(self, nx: int, ny: int, nz: int, h: float,
                 species: SP.SpeciesSet, dt: float,
                 bc: str = PERIODIC_XY_DIRICHLET_Z,
                 phi_lo=0.0, phi_hi=0.0,
                 background_charge_density: float = 0.0,
                 mcc=None, mcc_species: SP.ParticleSpecies = None,
                 rng: np.random.Generator = None, seed: int = 0) -> None:
        self.shape = (int(nx), int(ny), int(nz))
        self.h = float(h)
        self.dt = float(dt)
        self.bc = bc
        self.species = species
        self.phi_lo = phi_lo
        self.phi_hi = phi_hi
        # a uniform immobile background [C/m^3] -- the neutralising ion jelly
        # of a periodic electron-only test.  Real ion species are particles.
        self.rho_background = float(background_charge_density)
        self.solver = Poisson3D(nx, ny, nz, h, bc=bc)
        self.L = np.array([nx * h, ny * h,
                           (nz * h) if bc == PERIODIC else ((nz - 1) * h)])
        self.step_count = 0
        self.time = 0.0
        self.wall_charge = 0.0
        self.wall_energy_J = 0.0
        # per-electrode attribution: (species name, 'lo'|'hi') -> absorbed
        # real-particle count (sum of macro weights).  The M5 gamma process
        # and any per-electrode current measurement read this.
        self.wall_flux_w = {}
        # M5 hooks: called as hook(pop, idx_lo, idx_hi, self) with the
        # indices of the particles ABOUT to be absorbed, before removal.
        # A hook must not spawn into the population it is called about.
        self.wall_hooks = []
        # injectors: obj.inject(self) runs at the top of every step (M5
        # flux emitter; M9 beam sources can register the same way).
        self.injectors = []
        self.phi = np.zeros(self.shape)
        self.E = (np.zeros(self.shape), np.zeros(self.shape),
                  np.zeros(self.shape))
        self.history = {"t": [], "ke_J": [], "ue_J": [], "px": [],
                        "mode_re": [], "mode_im": [], "n_macro": [],
                        "n_coll": []}
        self._mode_axis = 2
        self._mode_number = 1

        # --- collisions (missing module 3, electron half) ---
        self.mcc = mcc
        self.rng = np.random.default_rng(seed) if rng is None else rng
        self.mcc_species = mcc_species
        if mcc is not None and mcc_species is None:
            self.mcc_species = self._default_mcc_target()
        # immobile charge from ionization/attachment events whose product has
        # NO registered particle species -- deposited on the grid so the
        # solve still sees it.  With ion species registered this stays zero.
        self.rho_heavy = np.zeros(self.shape)
        self.collision_tally = None
        self.collision_energy_eV = None
        self.gas_reservoir_eV = 0.0

    def _default_mcc_target(self) -> SP.ParticleSpecies:
        """The electron population, if there is exactly one."""
        from .units import M_E
        cands = [p for p in self.species if abs(p.species.mass / M_E - 1.0) < 1e-9]
        if len(cands) != 1:
            raise ValueError(
                f"cannot infer the MCC target: found {len(cands)} electron "
                "populations -- pass mcc_species= explicitly")
        return cands[0]

    # ------------------------------------------------------------ the cycle
    def deposit(self) -> np.ndarray:
        rho = np.zeros(self.shape)
        for p in self.species:
            if p.n:
                EF.deposit_rho(self.shape, self.h, p.x[:p.n], p.y[:p.n],
                               p.z[:p.n], p.charge_weight(), self.bc, out=rho)
        if self.rho_background:
            rho += self.rho_background
        if self.mcc is not None:
            rho += self.rho_heavy
        return rho

    def _bias(self, v):
        """Electrode potential right now: a number, or phi(t) if callable
        (the C3/RF path -- a time-varying applied waveform)."""
        return v(self.time) if callable(v) else v

    def solve_fields(self, rho: np.ndarray = None) -> None:
        rho = self.deposit() if rho is None else rho
        self._rho = rho                      # cached for density_mode()
        self.phi = self.solver.solve(rho, self._bias(self.phi_lo),
                                     self._bias(self.phi_hi))
        self.E = EF.efield_from_phi(self.phi, self.h, self.bc)

    def prepare(self) -> None:
        """Retard velocities by half a step so the leapfrog starts centred."""
        self.solve_fields()
        for p in self.species:
            if not p.n:
                continue
            ex, ey, ez = EF.gather_E(*self.E, self.h, p.x[:p.n], p.y[:p.n],
                                     p.z[:p.n], self.bc)
            SP.push_es(p, ex, ey, ez, -0.5 * self.dt)

    def step(self) -> None:
        for inj in self.injectors:
            inj.inject(self)
        self.solve_fields()
        for p in self.species:
            k = max(int(getattr(p, "subcycle", 1)), 1)
            if self.step_count % k:
                continue                       # this species waits its turn
            if not p.n:
                continue
            dt_p = self.dt * k
            ex, ey, ez = EF.gather_E(*self.E, self.h, p.x[:p.n], p.y[:p.n],
                                     p.z[:p.n], self.bc)
            SP.push_es(p, ex, ey, ez, dt_p)
            SP.move(p, dt_p)
            self._apply_boundaries(p)
        self._collide()
        self.step_count += 1
        self.time += self.dt

    def _collide(self) -> None:
        """Run the MCC over the electron population for one field step.

        Order matters: collisions come AFTER the push and move, so the
        operator sees the velocities the field just produced and the
        secondaries it creates start life at a real position inside the
        domain.  Ionization/attachment products go into their registered
        species; anything unregistered is deposited as immobile charge so
        the space-charge accounting stays closed either way.
        """
        if self.mcc is None or self.mcc_species is None:
            return
        p = self.mcc_species
        if p.n == 0:
            return
        ev = self.mcc.collide_adaptive(p, self.dt, self.rng)
        if self.collision_tally is None:
            self.collision_tally = np.zeros_like(ev["tally"])
            self.collision_energy_eV = np.zeros_like(ev["energy_eV"])
        self.collision_tally += ev["tally"]
        self.collision_energy_eV += ev["energy_eV"]
        self.gas_reservoir_eV += ev["gas_reservoir_eV"]
        if ev["ion_w"].size:
            ix, iy, iz = ev["ion_pos"]
            EF.deposit_rho(self.shape, self.h, ix, iy, iz,
                           +EF_E_CHARGE * ev["ion_w"], self.bc,
                           out=self.rho_heavy)
        if ev["att_w"].size:
            ax, ay, az = ev["att_pos"]
            EF.deposit_rho(self.shape, self.h, ax, ay, az,
                           -EF_E_CHARGE * ev["att_w"], self.bc,
                           out=self.rho_heavy)
        for pop in self.species:
            self._apply_boundaries(pop)

    def run(self, n_steps: int, record_every: int = 1) -> dict:
        self.prepare()
        self.record()
        for _ in range(int(n_steps)):
            self.step()
            if self.step_count % max(int(record_every), 1) == 0:
                self.record()
        return {k: np.asarray(v) for k, v in self.history.items()}

    # -------------------------------------------------------- boundaries
    def _apply_boundaries(self, p: SP.ParticleSpecies) -> None:
        n = p.n
        if not n:
            return
        p.x[:n] = np.mod(p.x[:n], self.L[0])
        p.y[:n] = np.mod(p.y[:n], self.L[1])
        if self.bc == PERIODIC:
            p.z[:n] = np.mod(p.z[:n], self.L[2])
            return
        # absorbing electrodes: tally the charge and energy they collect, then
        # remove the particle.  (Surface charge feeding back into the circuit
        # is missing module 4; this only BOOKKEEPS it so nothing is lost
        # silently.)
        lo = p.z[:n] < 0.0
        hi = p.z[:n] > self.L[2]
        gone = lo | hi
        if not gone.any():
            return
        idx = np.nonzero(gone)[0]
        idx_lo = np.nonzero(lo)[0]
        idx_hi = np.nonzero(hi)[0]
        name = p.species.name
        for side, ii in (("lo", idx_lo), ("hi", idx_hi)):
            if ii.size:
                key = (name, side)
                self.wall_flux_w[key] = (self.wall_flux_w.get(key, 0.0)
                                         + float(np.sum(p.w[ii])))
        for hook in self.wall_hooks:
            hook(p, idx_lo, idx_hi, self)
        self.wall_charge += float(p.species.charge * np.sum(p.w[idx]))
        v2 = p.vx[idx] ** 2 + p.vy[idx] ** 2 + p.vz[idx] ** 2
        self.wall_energy_J += float(0.5 * p.species.mass
                                    * np.sum(p.w[idx] * v2))
        p.alive[:n][idx] = False
        p.compact()

    # ------------------------------------------------------- diagnostics
    def set_mode(self, axis: int = 2, mode: int = 1) -> None:
        """Which Fourier mode :meth:`record` tracks (the seeded one)."""
        self._mode_axis = int(axis)
        self._mode_number = int(mode)

    def density_mode(self) -> complex:
        """Complex amplitude of the tracked density mode, from the grid.

        Averages the deposited charge density over the two transverse
        directions and takes the DFT coefficient of the seeded mode along the
        chosen axis.  Signed and linear in the perturbation, so it oscillates
        at omega_pe (the field ENERGY would oscillate at 2 omega_pe).

        Uses the density the last :meth:`solve_fields` already deposited
        rather than re-depositing -- a re-deposit per recorded step doubled
        the cost of every run for an identical answer.
        """
        rho = getattr(self, "_rho", None)
        if rho is None:
            rho = self.deposit()
        ax = self._mode_axis
        other = tuple(a for a in (0, 1, 2) if a != ax)
        line = rho.mean(axis=other)
        n = line.size
        k = 2.0 * np.pi * self._mode_number / n
        phase = np.exp(-1j * k * np.arange(n))
        return complex(np.sum(line * phase) / n)

    def field_energy_J(self) -> float:
        return EF.field_energy_J(*self.E, self.h)

    def record(self) -> None:
        amp = self.density_mode()
        self.history["t"].append(self.time)
        self.history["ke_J"].append(self.species.kinetic_energy_J())
        self.history["ue_J"].append(self.field_energy_J())
        self.history["px"].append(float(self.species.momentum()[2]))
        self.history["mode_re"].append(amp.real)
        self.history["mode_im"].append(amp.imag)
        self.history["n_macro"].append(int(sum(p.n for p in self.species)))
        self.history["n_coll"].append(
            int(self.collision_tally.sum()) if self.collision_tally is not None else 0)

    # --------------------------------------------------------- stability
    def stability(self, n_e_m3: float, Te_eV: float) -> dict:
        """The two explicit-PIC walls, reported (never silently violated).

        Returns dx/lambda_D and omega_pe*dt with the thresholds this package
        is only valid below.  M7 (the implicit mover) is what relaxes them;
        until then a run outside these bounds is not trustworthy and the
        caller is expected to say so.
        """
        ld = debye_length(n_e_m3, Te_eV)
        wpe = plasma_frequency(n_e_m3)
        return {"dx_over_debye": self.h / ld if ld > 0 else np.inf,
                "wpe_dt": wpe * self.dt,
                "resolved": bool(self.h < ld and wpe * self.dt < 2.0),
                "lambda_D_m": ld, "omega_pe": wpe}


def uniform_background_for(pop: SP.ParticleSpecies, volume_m3: float) -> float:
    """The immobile background charge density [C/m^3] that neutralises ``pop``.

    Used by the periodic verification runs, where the ions are a fixed jelly
    rather than particles: the fully periodic Poisson solve has no solution
    unless the total charge is zero (see the Known-bad list).
    """
    return -pop.species.charge * pop.real_count() / volume_m3
