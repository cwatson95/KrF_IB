"""MISSING MODULE 2 -- mobile, multi-species kinetic particles.

In ``ARPAE/pic3d`` only electrons are kinetic; ions and anions are immobile
charge smeared on the grid ("by scope: ns windows").  That is defensible for a
wave instrument and fatal for a discharge: with frozen ions there is no
ambipolar field, no sheath formation, no ion transit and no ion flux for
secondary emission to act on.

Here every species is an ordinary kinetic macro-particle population carrying
its own charge and mass.  The push takes q/m from the species record, so the
same Boris routine that moves electrons moves Ar+ -- only the mass ratio (and
therefore the sensible timestep) differs.

SUBCYCLING.  Heavy species barely move on an electron timestep, so each
species carries a ``subcycle`` count k: it is advanced once every k field
steps with a step of k*dt, seeing the field held fixed in between.
:func:`ion_subcycle_ratio` suggests k from the mass ratio; it is a suggestion,
and the caller is expected to CFL-check it (a gate compares subcycled and
un-subcycled trajectories).

Weights: ``w`` is the number of real particles a macro-particle represents.
Mixed weights within a species are allowed (ionization will need them); every
diagnostic here is weight-aware.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .units import AMU, E_CHARGE, EV, M_E

_FIELDS = ("x", "y", "z", "vx", "vy", "vz", "w")


@dataclass(frozen=True)
class Species:
    """A particle species: name, signed charge [C], mass [kg]."""

    name: str
    charge: float
    mass: float

    @property
    def q_over_m(self) -> float:
        return self.charge / self.mass

    @classmethod
    def electron(cls) -> "Species":
        return cls("e-", -E_CHARGE, M_E)

    @classmethod
    def ion(cls, name: str, mass_amu: float, charge_number: int = 1
            ) -> "Species":
        """A singly (or multiply) charged atomic/molecular ion."""
        return cls(name, charge_number * E_CHARGE, mass_amu * AMU)


# convenient, provenance-free atomic masses [amu] for the gas set this
# project uses; the CHEMISTRY lives elsewhere -- these are kinematics only
ION_MASSES_AMU = {"Ar+": 39.948, "Kr+": 83.798, "He+": 4.002602,
                  "F2+": 37.9968, "F-": 18.9984}


def ion_subcycle_ratio(species: Species, dt: float, h: float,
                       v_max: float = None, safety: float = 0.5) -> int:
    """Suggest a subcycle count k for a heavy species.

    Two considerations: (a) the species moves sqrt(m_i/m_e) times more slowly
    than an electron at the same energy, so k ~ sqrt(m_i/m_e) is the natural
    scale; (b) a CFL cap, k v_max dt <= safety * h, so a subcycled particle
    still cannot cross a cell in one of its own steps.
    """
    k_mass = int(max(1, round(np.sqrt(species.mass / M_E))))
    if v_max is not None and v_max > 0.0:
        k_cfl = int(max(1, np.floor(safety * h / (v_max * dt))))
        return int(max(1, min(k_mass, k_cfl)))
    return k_mass


class ParticleSpecies:
    """Structure-of-arrays store for one species."""

    def __init__(self, species: Species, capacity: int = 1024) -> None:
        self.species = species
        capacity = max(int(capacity), 1)
        self.capacity = capacity
        self.n = 0
        self.subcycle = 1
        for name in _FIELDS:
            setattr(self, name, np.zeros(capacity, dtype=float))
        self.alive = np.zeros(capacity, dtype=bool)

    # ------------------------------------------------------------- loaders
    @classmethod
    def maxwellian_box(cls, species: Species, n_macro: int, lo, hi,
                       n_real_m3: float, T_eV: float,
                       rng: np.random.Generator) -> "ParticleSpecies":
        """Uniform-in-space, Maxwellian-in-velocity load of a box."""
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        vol = float(np.prod(hi - lo))
        n_macro = int(n_macro)
        p = cls(species, capacity=max(int(1.5 * n_macro) + 16, 16))
        pos = lo[None, :] + rng.random((n_macro, 3)) * (hi - lo)[None, :]
        vth = np.sqrt(max(T_eV, 0.0) * EV / species.mass)
        vel = rng.normal(0.0, vth, size=(n_macro, 3))
        p._install(pos, vel, n_real_m3 * vol / max(n_macro, 1))
        return p

    @classmethod
    def cold_wave(cls, species: Species, n_macro: int, lo, hi,
                  n_real_m3: float, amplitude: float = 0.0,
                  mode: int = 1, axis: int = 2,
                  rng: np.random.Generator = None) -> "ParticleSpecies":
        """Cold, uniform load with a sinusoidal DISPLACEMENT along ``axis``.

        The standard quiet start for a plasma-oscillation test: a density
        perturbation dn/n = -A k cos(k s) is produced by displacing an
        otherwise uniform load, ``s -> s + A sin(k s)``, which seeds one clean
        mode instead of relying on shot noise.  Velocities are exactly zero,
        so the measured oscillation frequency is the cold-plasma omega_pe with
        no thermal (Bohm-Gross) correction to disentangle.
        """
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        vol = float(np.prod(hi - lo))
        n_macro = int(n_macro)
        p = cls(species, capacity=max(int(1.5 * n_macro) + 16, 16))
        # a low-discrepancy (regular) load in each direction keeps the noise
        # floor far below the seeded mode
        per_axis = int(round(n_macro ** (1.0 / 3.0)))
        if per_axis ** 3 != n_macro:
            raise ValueError(
                f"cold_wave wants a perfect cube of macro-particles, got "
                f"{n_macro} (nearest: {per_axis ** 3})")
        g = (np.arange(per_axis) + 0.5) / per_axis
        gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
        pos = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))
        pos = lo[None, :] + pos * (hi - lo)[None, :]
        if amplitude:
            L = float(hi[axis] - lo[axis])
            k = 2.0 * np.pi * mode / L
            pos[:, axis] += amplitude * np.sin(k * (pos[:, axis] - lo[axis]))
        vel = np.zeros((n_macro, 3))
        p._install(pos, vel, n_real_m3 * vol / max(n_macro, 1))
        return p

    def _install(self, pos, vel, weight) -> None:
        n = pos.shape[0]
        self.n = n
        self.x[:n], self.y[:n], self.z[:n] = pos.T
        self.vx[:n], self.vy[:n], self.vz[:n] = vel.T
        self.w[:n] = weight
        self.alive[:n] = True

    # --------------------------------------------------------- diagnostics
    def charge_weight(self) -> np.ndarray:
        """q*w per live macro-particle [C] -- what the deposit consumes."""
        return self.species.charge * self.w[: self.n]

    def real_count(self) -> float:
        return float(np.sum(self.w[: self.n]))

    def energy_eV(self) -> np.ndarray:
        n = self.n
        v2 = self.vx[:n] ** 2 + self.vy[:n] ** 2 + self.vz[:n] ** 2
        return 0.5 * self.species.mass * v2 / EV

    def kinetic_energy_J(self) -> float:
        n = self.n
        v2 = self.vx[:n] ** 2 + self.vy[:n] ** 2 + self.vz[:n] ** 2
        return float(0.5 * self.species.mass * np.sum(self.w[:n] * v2))

    def momentum(self) -> np.ndarray:
        n = self.n
        m = self.species.mass
        w = self.w[:n]
        return m * np.array([float(np.sum(w * self.vx[:n])),
                             float(np.sum(w * self.vy[:n])),
                             float(np.sum(w * self.vz[:n]))])

    def temperature_eV(self) -> float:
        """(2/3)<eps> of the weighted population [eV] (drifting frame kept)."""
        n = self.n
        if n == 0:
            return 0.0
        w = self.w[:n]
        return float((2.0 / 3.0) * np.sum(w * self.energy_eV()) / np.sum(w))

    # ------------------------------------------------------ grow / compact
    def _ensure_capacity(self, extra: int) -> None:
        need = self.n + extra
        if need <= self.capacity:
            return
        new_cap = max(need, 2 * self.capacity)
        for name in _FIELDS:
            arr = getattr(self, name)
            grown = np.zeros(new_cap, dtype=float)
            grown[: self.capacity] = arr
            setattr(self, name, grown)
        grown_alive = np.zeros(new_cap, dtype=bool)
        grown_alive[: self.capacity] = self.alive
        self.alive = grown_alive
        self.capacity = new_cap

    def spawn(self, x, y, z, vx, vy, vz, w) -> None:
        x = np.atleast_1d(np.asarray(x, dtype=float))
        m = x.size
        if m == 0:
            return
        self._ensure_capacity(m)
        s = slice(self.n, self.n + m)
        for name, val in zip(_FIELDS, (x, y, z, vx, vy, vz, w)):
            getattr(self, name)[s] = np.broadcast_to(
                np.asarray(val, dtype=float), (m,))
        self.alive[s] = True
        self.n += m

    def kill(self, mask) -> None:
        self.alive[: self.n][np.asarray(mask, dtype=bool)] = False

    def compact(self) -> None:
        live = self.alive[: self.n]
        if live.all():
            return
        idx = np.nonzero(live)[0]
        m = idx.size
        for name in _FIELDS:
            arr = getattr(self, name)
            arr[:m] = arr[idx]
        self.alive[:m] = True
        self.alive[m: self.n] = False
        self.n = m


class SpeciesSet:
    """The mixture: an ordered collection of :class:`ParticleSpecies`."""

    def __init__(self, populations=()) -> None:
        self.populations = list(populations)

    def __iter__(self):
        return iter(self.populations)

    def __len__(self) -> int:
        return len(self.populations)

    def add(self, pop: ParticleSpecies) -> ParticleSpecies:
        self.populations.append(pop)
        return pop

    def by_name(self, name: str) -> ParticleSpecies:
        for p in self.populations:
            if p.species.name == name:
                return p
        raise KeyError(f"no species named {name!r}")

    def kinetic_energy_J(self) -> float:
        return float(sum(p.kinetic_energy_J() for p in self.populations))

    def momentum(self) -> np.ndarray:
        tot = np.zeros(3)
        for p in self.populations:
            tot += p.momentum()
        return tot

    def net_charge(self) -> float:
        return float(sum(p.species.charge * p.real_count()
                         for p in self.populations))


# --------------------------------------------------------------- the movers
def push_es(p: ParticleSpecies, ex, ey, ez, dt: float) -> None:
    """Electrostatic leapfrog kick: v^{n-1/2} -> v^{n+1/2} under E.

    Positions are NOT advanced here -- the caller moves them, so a future
    charge-conserving deposit can see both endpoints of the flight.
    """
    n = p.n
    if n == 0:
        return
    a = p.species.q_over_m * dt
    p.vx[:n] += a * ex
    p.vy[:n] += a * ey
    p.vz[:n] += a * ez


def push_boris(p: ParticleSpecies, ex, ey, ez, bx, by, bz, dt: float) -> None:
    """Boris rotation under E and B: half kick, magnetic rotation, half kick.

    Exactly speed-conserving when E = 0 (pure gyration), and the rotation
    angle matches 2 atan(omega_c dt / 2) -- both gated.  Identical algebra to
    ``ARPAE/pic3d/particles3d.boris_push`` but with q/m taken from the
    species, so it moves ions as readily as electrons.
    """
    n = p.n
    if n == 0:
        return
    hqmdt = 0.5 * p.species.q_over_m * dt
    vmx = p.vx[:n] + hqmdt * ex
    vmy = p.vy[:n] + hqmdt * ey
    vmz = p.vz[:n] + hqmdt * ez
    tx, ty, tz = hqmdt * bx, hqmdt * by, hqmdt * bz
    t2 = tx * tx + ty * ty + tz * tz
    sx, sy, sz = (2.0 * tx / (1.0 + t2), 2.0 * ty / (1.0 + t2),
                  2.0 * tz / (1.0 + t2))
    vpx = vmx + (vmy * tz - vmz * ty)
    vpy = vmy + (vmz * tx - vmx * tz)
    vpz = vmz + (vmx * ty - vmy * tx)
    vfx = vmx + (vpy * sz - vpz * sy)
    vfy = vmy + (vpz * sx - vpx * sz)
    vfz = vmz + (vpx * sy - vpy * sx)
    p.vx[:n] = vfx + hqmdt * ex
    p.vy[:n] = vfy + hqmdt * ey
    p.vz[:n] = vfz + hqmdt * ez


def move(p: ParticleSpecies, dt: float) -> None:
    """Advance positions by v dt (the drift half of the leapfrog)."""
    n = p.n
    if n == 0:
        return
    p.x[:n] += p.vx[:n] * dt
    p.y[:n] += p.vy[:n] * dt
    p.z[:n] += p.vz[:n] * dt
