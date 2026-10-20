"""MISSING MODULE 3 (electron half) -- null-collision Monte-Carlo collisions.

The MC of PIC-MCC.  Vahedi-Surendra null-collision selection over the LXCat
set, ADAPTED from ``ARPAE/pic3d/mcc3d.py`` (itself adapted from the tested
``pic1d/mcc_ref.py`` @ 4e2a9fc, which matches LoKI-MC swarm rates to x1.6).
The selection, kinematics, angular models, ionization sharing and adaptive
substepping are deliberately identical to the source module, so its gates
still mean something here.

THE ONE REAL UPGRADE: ionization can produce a REAL ION MACRO-PARTICLE.
ARPAE's engine has immobile ions, so its MCC deposits +e of frozen charge at
the birth site.  This package has kinetic ions (missing module 2), so an
ionization event may instead spawn into an actual :class:`ParticleSpecies` --
which is what a sheath, an ion flux to a wall, and therefore secondary
emission (M5) all require.  The immobile-charge path is kept for callers that
do not want to pay for ion particles; both are gated.

ENERGY BOOKKEEPING, and where it stops being exact
--------------------------------------------------
The ELECTRON ledger is exact and gated at 1e-10:

    KE_e(before) - KE_e(after, including secondaries) = sum(ledger)

with elastic -> gas heating, thresholds -> excitation/ionization potential,
attachment -> the removed electron's kinetic energy.

A newborn ion is a separate matter.  It inherits the velocity of the neutral
it was made from, so it arrives carrying THERMAL energy drawn from the gas
reservoir -- energy that was already in the system but not in this model's
particle set.  That is booked separately in ``gas_reservoir_eV`` rather than
folded into the electron closure, because hiding a reservoir term inside an
"exact" invariant is how an energy audit stops meaning anything.  Total energy
is conserved only once that term is counted; the static-gas assumption
(``gas.py``) is what licenses drawing on it.

SCATTERING CAVEAT (inherited, not solved): isotropic scattering off
momentum-transfer cross sections.  ``pic1d/README_LIMITS.md`` §2.2 is the
statement of what that does and does not establish.
"""

from __future__ import annotations

import numpy as np

from .cross_sections import CrossSectionSet
from .gas import GasMixture
from .species3d import ParticleSpecies
from .units import E_CHARGE, EV, K_B, M_E


def _speed_to_energy_eV(speed2: np.ndarray) -> np.ndarray:
    return 0.5 * M_E * speed2 / E_CHARGE


def _energy_to_speed(eps_eV: np.ndarray) -> np.ndarray:
    return np.sqrt(2.0 * np.maximum(eps_eV, 0.0) * E_CHARGE / M_E)


def _isotropic_dirs(rng: np.random.Generator, m: int) -> np.ndarray:
    cos_t = 1.0 - 2.0 * rng.random(m)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t ** 2))
    phi = 2.0 * np.pi * rng.random(m)
    return np.column_stack((sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t))


def _scatter_dirs(vhat: np.ndarray, cos_chi: np.ndarray,
                  eta: np.ndarray) -> np.ndarray:
    sin_chi = np.sqrt(np.maximum(0.0, 1.0 - cos_chi ** 2))
    ref = np.zeros_like(vhat)
    near_x = np.abs(vhat[:, 0]) > 0.9
    ref[:, 0] = 1.0
    ref[near_x, 0] = 0.0
    ref[near_x, 1] = 1.0
    e1 = ref - np.sum(ref * vhat, axis=1)[:, None] * vhat
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(vhat, e1)
    return (cos_chi[:, None] * vhat
            + (sin_chi * np.cos(eta))[:, None] * e1
            + (sin_chi * np.sin(eta))[:, None] * e2)


def sample_cos_chi(model: str, eps: np.ndarray, rng: np.random.Generator,
                   e_screen: float = 27.21) -> np.ndarray:
    """'isotropic' | 'screened' | 'surendra' polar-angle sampling (pic1d)."""
    eps = np.asarray(eps, dtype=float)
    if model == "screened":
        eta = e_screen / np.maximum(eps, 1e-6)
        p = rng.random(eps.size)
        return (eta + 1.0 - (2.0 * eta + 1.0) * p) / (eta + 1.0 - p)
    if model == "surendra":
        p = rng.random(eps.size)
        return (2.0 + eps - 2.0 * (1.0 + eps) ** p) / np.maximum(eps, 1e-12)
    return 1.0 - 2.0 * rng.random(eps.size)


class MCC3D:
    """Null-collision MCC over an electron :class:`ParticleSpecies`.

    ``ion_species`` optionally maps a target gas name (as it appears in the
    LXCat file, e.g. ``"Ar"``) to the :class:`ParticleSpecies` its ions should
    be born into.  Targets absent from the map fall back to immobile-charge
    bookkeeping (event positions/weights are still reported, so the caller can
    deposit them on the grid).  ``anion_species`` does the same for
    attachment.
    """

    def __init__(self, cross_sections: CrossSectionSet,
                 mixture: GasMixture = None,
                 ion_species: dict = None, anion_species: dict = None,
                 ionization_sharing: str = "equal",
                 angular_model: dict = None, e_screen: float = 27.21,
                 eps_split: float = 3.0, p_sub_target: float = 0.25,
                 max_substeps: int = 64, born_ion_thermal: bool = True):
        self.cs = cross_sections
        self.mixture = mixture
        if ionization_sharing not in ("equal", "uniform"):
            raise ValueError(
                f"unknown ionization_sharing {ionization_sharing!r}")
        self.ionization_sharing = ionization_sharing
        self.angular_model = dict(angular_model or {})
        self.e_screen = float(e_screen)
        self.eps_split = float(eps_split)
        self.p_sub_target = float(p_sub_target)
        self.max_substeps = int(max_substeps)
        self.ion_species = dict(ion_species or {})
        self.anion_species = dict(anion_species or {})
        self.born_ion_thermal = bool(born_ion_thermal)
        if self.born_ion_thermal and (self.ion_species or self.anion_species) \
                and mixture is None:
            raise ValueError(
                "born_ion_thermal needs a GasMixture to draw the neutral "
                "temperature from; pass mixture= or set born_ion_thermal=False")

    # ------------------------------------------------------------ selection
    def _sample_cos_chi(self, kind, eps, rng):
        return sample_cos_chi(self.angular_model.get(kind, "isotropic"),
                              eps, rng, self.e_screen)

    def select(self, eps: np.ndarray, dt: float, rng: np.random.Generator,
               nu_prime: float):
        """Null-collision selection on the energy array ``eps`` (pic1d form)."""
        n = eps.size
        p_coll = 1.0 - np.exp(-nu_prime * dt)
        cand_idx = np.nonzero(rng.random(n) < p_coll)[0]
        if cand_idx.size == 0:
            return cand_idx, np.empty(0, dtype=np.intp)
        nu = self.cs.nu_matrix(eps[cand_idx])
        cumulative = np.cumsum(nu, axis=0) / nu_prime
        r = rng.random(cand_idx.size)
        is_real = r < cumulative[-1]
        proc = np.full(cand_idx.size, -1, dtype=np.intp)
        first = (cumulative > r[None, :]).argmax(axis=0)
        proc[is_real] = first[is_real]
        return cand_idx, proc

    # ------------------------------------------------------------- collide
    def collide_indices(self, p: ParticleSpecies, idx: np.ndarray, dt: float,
                        rng: np.random.Generator, nu_prime: float,
                        out: dict) -> None:
        """One MCC substep over live electrons ``idx`` (global indices)."""
        if idx.size == 0:
            return
        vx, vy, vz = p.vx[idx], p.vy[idx], p.vz[idx]
        speed2 = vx * vx + vy * vy + vz * vz
        eps = _speed_to_energy_eV(speed2)
        ci, proc = self.select(eps, dt, rng, nu_prime)
        if ci.size == 0:
            return
        cand = idx[ci]
        real = proc >= 0
        if not real.any():
            return

        v = np.column_stack((p.vx[cand], p.vy[cand], p.vz[cand]))
        speed = np.sqrt(np.maximum(speed2[ci], 1e-300))
        vhat = v / speed[:, None]
        eps_c = eps[ci]
        w_c = p.w[cand]

        kill_mask = np.zeros(cand.size, dtype=bool)
        for k in range(self.cs.n_proc):
            sel = np.nonzero(proc == k)[0]
            if sel.size == 0:
                continue
            out["tally"][k] += int(sel.size)
            kind = self.cs.kinds[k]
            thr = self.cs.thresholds[k]
            mr = self.cs.mass_ratios[k]
            target = self.cs.targets[k]
            cos_chi = self._sample_cos_chi(kind, eps_c[sel], rng)
            eta = 2.0 * np.pi * rng.random(sel.size)

            if kind == "attachment":
                kill_mask[sel] = True
                gi = cand[sel]
                out["energy_eV"][k] += float(np.sum(w_c[sel] * eps_c[sel]))
                out["att"].append((p.x[gi].copy(), p.y[gi].copy(),
                                   p.z[gi].copy(), p.w[gi].copy(), target))
                continue

            if kind == "elastic":
                d_eps = eps_c[sel] * (2.0 * mr * (1.0 - cos_chi))
                eps_new = eps_c[sel] - d_eps
                out["energy_eV"][k] += float(np.sum(w_c[sel] * d_eps))
            elif kind == "excitation":
                eps_new = eps_c[sel] - thr
                out["energy_eV"][k] += float(np.sum(w_c[sel] * thr))
            elif kind == "ionization":
                eps_avail = np.maximum(eps_c[sel] - thr, 0.0)
                if self.ionization_sharing == "uniform":
                    frac = rng.random(sel.size)
                else:
                    frac = np.full(sel.size, 0.5)
                eps_ej = frac * eps_avail
                eps_new = (1.0 - frac) * eps_avail
                out["energy_eV"][k] += float(np.sum(w_c[sel] * thr))
                gi = cand[sel]
                sec_dir = _isotropic_dirs(rng, sel.size)
                sec_v = _energy_to_speed(eps_ej)[:, None] * sec_dir
                out["ion"].append((p.x[gi].copy(), p.y[gi].copy(),
                                   p.z[gi].copy(), p.w[gi].copy(), target))
                out["_spawn"].append((p.x[gi].copy(), p.y[gi].copy(),
                                      p.z[gi].copy(), sec_v, p.w[gi].copy()))
            else:
                continue

            new_dir = _scatter_dirs(vhat[sel], cos_chi, eta)
            new_v = _energy_to_speed(eps_new)[:, None] * new_dir
            gi = cand[sel]
            p.vx[gi] = new_v[:, 0]
            p.vy[gi] = new_v[:, 1]
            p.vz[gi] = new_v[:, 2]

        if kill_mask.any():
            p.alive[cand[kill_mask]] = False

    def collide_adaptive(self, p: ParticleSpecies, dt: float,
                         rng: np.random.Generator) -> dict:
        """Energy-grouped, substepped MCC over the whole population for ``dt``.

        At 1 atm the majorant nu' dt can exceed 1, so several collisions per
        field step are possible.  Electrons are split into cold/hot groups at
        ``eps_split``; each group gets its own majorant and substep count, so
        one hot secondary no longer drags the cold bulk through dozens of
        substeps.
        """
        out = {
            "tally": np.zeros(self.cs.n_proc, dtype=np.int64),
            "energy_eV": np.zeros(self.cs.n_proc),
            "ion": [], "att": [], "_spawn": [],
            "gas_reservoir_eV": 0.0,
            "n_substeps": (0, 0),
        }
        if p.n == 0:
            return self._finish(p, out, rng)
        eps = p.energy_eV()
        cold = np.nonzero(eps < self.eps_split)[0]
        hot = np.nonzero(eps >= self.eps_split)[0]
        subs = []
        for idx, eps_g in ((cold, eps[cold]), (hot, eps[hot])):
            if idx.size == 0:
                subs.append(0)
                continue
            nu_p = self.cs.nu_majorant(float(eps_g.max()) * 1.1 + 0.5)
            k_sub = int(np.ceil(nu_p * dt / self.p_sub_target))
            k_sub = max(1, min(k_sub, self.max_substeps))
            subs.append(k_sub)
            dt_sub = dt / k_sub
            for _ in range(k_sub):
                live = idx[p.alive[idx]]
                if live.size == 0:
                    break
                self.collide_indices(p, live, dt_sub, rng, nu_p, out)
        out["n_substeps"] = tuple(subs)
        return self._finish(p, out, rng)

    # ------------------------------------------------------------- finish
    def _thermal_velocities(self, target: str, m: int,
                            rng: np.random.Generator) -> np.ndarray:
        """Velocity a newborn heavy particle inherits from its neutral."""
        if not self.born_ion_thermal or self.mixture is None:
            return np.zeros((m, 3))
        mass = self.mixture.mass_of(target)
        vth = np.sqrt(K_B * self.mixture.T_gas_K / mass)
        return rng.normal(0.0, vth, size=(m, 3))

    def _birth_heavies(self, records, table, rng, out, sign_note):
        """Spawn ions/anions into their species where one is registered."""
        leftover = []
        for x, y, z, w, target in records:
            pop = table.get(target)
            if pop is None:
                leftover.append((x, y, z, w, target))
                continue
            v = self._thermal_velocities(target, x.size, rng)
            pop.spawn(x, y, z, v[:, 0], v[:, 1], v[:, 2], w)
            # energy drawn from the gas reservoir, booked separately -- it is
            # NOT part of the electron closure (see the module docstring)
            m_t = (self.mixture.mass_of(target) if self.mixture is not None
                   else 0.0)
            if m_t:
                v2 = (v ** 2).sum(axis=1)
                out["gas_reservoir_eV"] += float(0.5 * m_t * np.sum(w * v2) / EV)
        return leftover

    def _finish(self, p: ParticleSpecies, out: dict,
                rng: np.random.Generator) -> dict:
        for sx, sy, sz, sv, sw in out.pop("_spawn"):
            p.spawn(sx, sy, sz, sv[:, 0], sv[:, 1], sv[:, 2], sw)
        if p.n and not p.alive[: p.n].all():
            p.compact()

        ion_left = self._birth_heavies(out["ion"], self.ion_species, rng, out,
                                       "+")
        att_left = self._birth_heavies(out["att"], self.anion_species, rng,
                                       out, "-")

        def _pack(records):
            if not records:
                z = np.zeros(0)
                return (z, z, z), z
            xs = np.concatenate([r[0] for r in records])
            ys = np.concatenate([r[1] for r in records])
            zs = np.concatenate([r[2] for r in records])
            ws = np.concatenate([r[3] for r in records])
            return (xs, ys, zs), ws

        # events whose species is NOT registered stay as immobile-charge
        # records for the caller to deposit (the ARPAE path)
        out["ion_pos"], out["ion_w"] = _pack(ion_left)
        out["att_pos"], out["att_w"] = _pack(att_left)
        out.pop("ion")
        out.pop("att")
        return out

    # --------------------------------------------------------- diagnostics
    def rate_coefficients(self, Te_eV: float, n_sample: int = 200000,
                          rng: np.random.Generator = None) -> dict:
        """Maxwellian-averaged rate coefficients <sigma v> [m^3/s] per kind.

        A scoring diagnostic, not part of the loop: it prices the SAME tables
        the MCC samples, so a disagreement between this and a measured swarm
        rate is a bug in the sampling rather than in the data.
        """
        rng = np.random.default_rng(0) if rng is None else rng
        vth = np.sqrt(max(Te_eV, 1e-12) * EV / M_E)
        v = rng.normal(0.0, vth, size=(int(n_sample), 3))
        eps = 0.5 * M_E * (v ** 2).sum(axis=1) / EV
        nu = self.cs.nu_matrix(eps)                     # (n_proc, n_sample)
        out = {}
        for k in range(self.cs.n_proc):
            kind = self.cs.kinds[k]
            dens = self.cs.density[k]
            if dens <= 0:
                continue
            out[kind] = out.get(kind, 0.0) + float(np.mean(nu[k]) / dens)
        return out
