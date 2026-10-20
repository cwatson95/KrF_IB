"""MISSING MODULE 5 -- cathode surface emission (the gamma process + a flux
emitter).

``PIC_MCC_3D.txt`` section 2.1 item 5 / milestone M5: without secondary
emission a discharge is not self-sustaining -- the electron population decays
once the seed is consumed, because volume ionization alone cannot replace the
electrons lost to the walls.  This module closes that loop, exactly as the
plan words it:

    per-step ion areal flux to each electrode x gamma -> WHOLE super-particles
    emitted just inside the wall, directed into the gap, with the emission
    energy distribution stated.

Two emitters, one sampling rule:

* :class:`SecondaryEmitter` -- the gamma process.  Registers as an ``ESSim``
  wall hook; whenever the watched ion population is absorbed by the chosen
  electrode, gamma x (incident real weight) real electrons are queued and
  emitted as whole super-particles of weight ``w_macro`` through a fractional
  accumulator (nothing is rounded away; the remainder carries to the next
  hit).  gamma is a STATED PARAMETER of the surface, not a prediction -- real
  gamma depends on cathode material, state, and ion energy, and every result
  built on this module must say which gamma it assumed.

* :class:`FluxEmitter` -- a constant-current-density cathode (thermionic-style
  source), the piece the Child-Langmuir gate (plan gate G4.1) needs.
  Registers as an ``ESSim`` injector; each step it emits
  j x A x dt / (e x w_macro) super-particles through the same accumulator.

THE EMISSION ENERGY DISTRIBUTION (stated, as the plan requires): emitted
electrons leave the surface with a flux-weighted half-Maxwellian at
temperature ``T_emit_eV`` -- the distribution a thermal source actually emits
through a plane:

    f(v_n)  proportional to  v_n exp(-m v_n^2 / 2 e T)   (into the gap)
    v_n = v_th sqrt(-ln U),   v_th = sqrt(2 e T / m),   U ~ uniform(0,1]
    v_t ~ Normal(0, sqrt(e T / m))                       (each tangential)

Newborns are placed a random fraction of their first free flight inside the
domain (z = U' v_n dt), the standard emitted-during-the-step placement; their
first velocity is not half-step retarded, a stated injection artifact that
vanishes with dt against any gate tolerance used here.

WHAT THIS IS STILL NOT: no electron-impact SEE, no field emission, no
energy-dependent gamma curve, no dielectric charging (M4's surface half), and
ion-neutral collisions (section 2.1 item 3) remain open -- ions arrive at the
cathode by free fall, so gamma here cannot be made ion-energy-resolved yet.
"""

from __future__ import annotations

import numpy as np

from .species3d import ParticleSpecies
from .units import E_CHARGE, EV


def flux_maxwellian(mass_kg: float, T_eV: float, m: int,
                    rng: np.random.Generator, direction: float
                    ) -> tuple:
    """Sample ``m`` emitted velocities: flux-weighted normal, thermal tangent.

    ``direction`` is +1.0 for emission toward +z (a z=0 cathode) and -1.0
    for emission toward -z (a z=Lz cathode).
    """
    if m == 0:
        z = np.zeros(0)
        return z, z, z.copy()
    kT = T_eV * EV
    vth = np.sqrt(2.0 * kT / mass_kg)          # most-probable speed scale
    sig = np.sqrt(kT / mass_kg)                # 1-D thermal sigma
    u = rng.random(m)
    vn = vth * np.sqrt(-np.log(np.clip(1.0 - u, 1e-300, None)))
    vx = rng.normal(0.0, sig, m)
    vy = rng.normal(0.0, sig, m)
    return vx, vy, float(direction) * vn


class _Accumulator:
    """Whole super-particles out of a fractional real-electron budget."""

    def __init__(self, w_macro: float) -> None:
        self.w_macro = float(w_macro)
        if self.w_macro <= 0.0:
            raise ValueError("w_macro must be positive")
        self.frac = 0.0
        self.emitted_macro = 0
        self.emitted_real = 0.0

    def take(self, real_electrons: float) -> int:
        self.frac += float(real_electrons) / self.w_macro
        n = int(np.floor(self.frac))
        self.frac -= n
        self.emitted_macro += n
        self.emitted_real += n * self.w_macro
        return n


class _EmitterBase:
    def __init__(self, pop: ParticleSpecies, wall: str, T_emit_eV: float,
                 w_macro: float) -> None:
        if wall not in ("lo", "hi"):
            raise ValueError(f"wall must be 'lo' or 'hi', got {wall!r}")
        self.pop = pop
        self.wall = wall
        self.T_emit_eV = float(T_emit_eV)
        self.acc = _Accumulator(w_macro)

    def _emit(self, sim, n: int) -> None:
        if n <= 0:
            return
        rng = sim.rng
        direction = +1.0 if self.wall == "lo" else -1.0
        vx, vy, vz = flux_maxwellian(self.pop.species.mass, self.T_emit_eV,
                                     n, rng, direction)
        x = rng.random(n) * sim.L[0]
        y = rng.random(n) * sim.L[1]
        dz = rng.random(n) * np.abs(vz) * sim.dt
        z = dz if self.wall == "lo" else sim.L[2] - dz
        w = np.full(n, self.acc.w_macro)
        self.pop.spawn(x, y, z, vx, vy, vz, w)


class FluxEmitter(_EmitterBase):
    """Constant-current-density electrode source (an ``ESSim`` injector).

    ``j_A_m2`` is the emitted current density magnitude; the emitting area is
    the sim's full periodic cross-section Lx x Ly.  Space-charge limitation is
    then the SIMULATION's job: emit more than Child-Langmuir allows and the
    virtual cathode must turn the surplus back -- that is gate G4.1.
    """

    def __init__(self, pop: ParticleSpecies, j_A_m2: float,
                 wall: str = "lo", T_emit_eV: float = 0.5,
                 w_macro: float = 1.0) -> None:
        super().__init__(pop, wall, T_emit_eV, w_macro)
        self.j_A_m2 = float(j_A_m2)

    def inject(self, sim) -> None:
        area = sim.L[0] * sim.L[1]
        real = self.j_A_m2 * area * sim.dt / E_CHARGE
        self._emit(sim, self.acc.take(real))


class SecondaryEmitter(_EmitterBase):
    """The gamma process (an ``ESSim`` wall hook).

    Watches ``ion_pop`` impacts on the chosen electrode; per incident real
    ion, ``gamma`` real electrons enter the accumulator and come out as whole
    super-particles of weight ``w_macro``.  The hook must never spawn into
    the population it is called about, and does not: it watches ions and
    spawns electrons.
    """

    def __init__(self, pop: ParticleSpecies, ion_pop: ParticleSpecies,
                 gamma: float, wall: str = "lo", T_emit_eV: float = 2.0,
                 w_macro: float = 1.0) -> None:
        super().__init__(pop, wall, T_emit_eV, w_macro)
        if pop is ion_pop:
            raise ValueError("SecondaryEmitter must emit into a different "
                             "population than the one it watches")
        self.ion_pop = ion_pop
        self.gamma = float(gamma)
        self.incident_real = 0.0

    def wall_hit(self, pop: ParticleSpecies, idx_lo: np.ndarray,
                 idx_hi: np.ndarray, sim) -> None:
        if pop is not self.ion_pop:
            return
        idx = idx_lo if self.wall == "lo" else idx_hi
        if idx.size == 0:
            return
        incident = float(np.sum(pop.w[idx]))
        self.incident_real += incident
        self._emit(sim, self.acc.take(self.gamma * incident))
