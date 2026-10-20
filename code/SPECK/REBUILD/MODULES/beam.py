"""E-beam pump: stopping-power deposition, emergent-W yields, conserving split.

PHYSICS ADAPTED from ``PIC_MCC_3D/MODULES/beam3d.py`` (the gated M9 module)
and ``spectre/ebeam.py`` (the conserving per-pair ledger), with provenance:

* **Stopping power is first-principles** — the relativistic Bethe/Møller
  collisional formula with ICRU 37 mean excitation energies (beam3d's
  ``stopping_power``, verified there against NIST ESTAR to 2–6 % at
  600 keV). The 600 kV primary is transported as a CSDA ray losing S(E)·dx
  per cell — never as a particle, and never below the formula's ~1 keV
  validity floor (the remainder deposits locally; sub-percent of the range).
  The march uses the LOCAL cell density, so a density gradient reshapes the
  profile (the MESH.txt density-adaptivity ask, first rung).

* **The W-value is an OUTPUT.** ``spencer_fano_per_process`` runs the
  Fowler/Spencer–Fano yield recursion per process over the same LXCat sets
  LoKI-MC samples; W, the per-species yields, the Platzman subexcitation
  energy and the untracked-channel subfloor all EMERGE, and they satisfy the
  per-pair identity

      W = sum(produced species x E_form) + subfloor + eps_sec

  exactly (gated at 1e-9): the beam injects exactly P_dep — the discipline
  whose absence in phantom11 is PLAN (70)'s open ~18x defect
  (ORESTES_CLASS_REBUILD.txt §0.3).

* **The W-table fallback** (``w_table_yields``) carries spectre's cited Lee
  1989 ledger (W, f_exc, E_ion, E_exc per gas, subexcitation remainder) for
  runs without LXCat data. Which route a run used is recorded on the yields
  object and printed in the run report — a fallback may not masquerade as
  the first-principles result.

Literature W anchors (SCORE, never steer): Ar 26.4, Kr 24.2, Ne 36.6,
He 41.3 eV/ion pair. PIC_MCC_3D's measured misses (Ar +5 %, Kr +8 %,
He +33 %) are the expected class for equal-sharing Spencer–Fano on these
sets; the same misses are expected here and are pinned by tests, not hidden.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cross_sections import CrossSectionSet
from .species import e_form_of
from .units import EV, MEC2_EV, N_AVOGADRO, R_E

# ICRU 37 mean excitation energies I [eV], atomic number Z, molar mass
# [g/mol]. F2 as 2 x atomic F (I_F = 115 eV; the molecular correction is
# small next to the log). Values as in PIC_MCC_3D/beam3d (+ Ne there too,
# + Xe added here for the Chen trace: ICRU 37 I(Xe) = 482 eV).
_ICRU = {
    "He": (41.8, 2, 4.002602),
    "Ne": (137.0, 10, 20.1797),
    "Ar": (188.0, 18, 39.948),
    "Kr": (352.0, 36, 83.798),
    "Xe": (482.0, 54, 131.293),
    "F2": (115.0, 18, 37.9968),
}

# --- spectre's cited Lee-ledger constants (the W-table fallback) -----------
# Kannari 1985 Table I p.4312 beam ledger, verified first-hand 2026-08-04:
# Ar 26.2, Kr 24.3, F2 32 eV/ion pair; beam excitation = 0.28 x ionization
# for BOTH Ar and Kr. Kr stays 24.2 (the ICRU/lit anchor this package
# scores against; Kannari's 24.3 is within the source spread and noted).
# F2 33->32 and F_EXC(Kr) 0.5->0.28 adopted from Kannari (were EST).
W_EV = {"Ne": 36.6, "Ar": 26.2, "Kr": 24.2, "F2": 32.0, "Xe": 21.9}
F_EXC = {"Ne": 0.305, "Ar": 0.284, "Kr": 0.28, "F2": 0.0, "Xe": 0.5}
E_ION_EV = {"Ne": 21.56, "Ar": 15.76, "Kr": 14.0, "F2": 15.69, "Xe": 12.13}
E_EXC_EV = {"Ne": 16.7, "Ar": 11.6, "Kr": 9.92, "Xe": 8.32}
_W_TABLE_CITE = ("Lee 1989 Table I [25][26] W/E_exc ledger as audited into "
                 "spectre ebeam.py; W(Ar) 26.2, W(F2) 32, f_exc 0.28 "
                 "verified first-hand in Kannari 1985 Table I p.4312 "
                 "(2026-08-04)")

BEAM_ION = {"Ne": "Ne+", "Ar": "Ar+", "Kr": "Kr+", "F2": "F2+", "Xe": "Xe+"}
BEAM_EXC = {"Ne": "Ne*", "Ar": "Ar*", "Kr": "Kr*", "Xe": "Xe*"}
BEAM_EXC_HI = {"Ne": "Ne**", "Ar": "Ar**", "Kr": "Kr**"}
# lo/hi manifold split thresholds [eV] (spectre's Ne 3s/3p split at 18.0;
# Ar 4s/4p at 12.9; Kr 5s/5p at 11.0; F2 electronic min 2.5, vib below)
_LO_HI_SPLIT = {"Ne": 18.0, "Ar": 12.9, "Kr": 11.0, "Xe": 9.0}
_F2_ELEC_MIN = 2.5
_F2_D0_EV = 1.6            # F2 dissociation energy (2 x E_form(F) = 0.80)

LITERATURE_W_EV = {"He": 41.3, "Ne": 36.6, "Ar": 26.4, "Kr": 24.2,
                   "N2": 34.8}


# --------------------------------------------------------------- beam pulse
@dataclass(frozen=True)
class BeamSpec:
    """Hibachi-style beam through a foil: trapezoidal pulse.

    DE-AR class defaults: 600 kV, 53 A/cm2 at the gas, ~60 ns flat top
    (docs/DE-AR_targets.txt; spectre ebeam.Beam)."""

    voltage_kV: float = 600.0
    j_A_cm2: float = 53.0
    t_start_ns: float = 0.0
    t_flat_ns: float = 60.0
    t_rise_ns: float = 5.0
    foil_transmission: float = 1.0
    axis: int = 0                   # injection axis (laser axis is z = 2)
    from_low: bool = True

    @property
    def energy_eV(self) -> float:
        return 1.0e3 * float(self.voltage_kV)

    def current_density_A_m2(self, t_s: float) -> float:
        t = t_s * 1e9 - self.t_start_ns
        j0 = self.j_A_cm2 * 1e4 * self.foil_transmission
        r, flat = max(self.t_rise_ns, 0.0), self.t_flat_ns
        if t <= 0.0 or t >= flat + 2 * r:
            return 0.0
        if r > 0.0 and t < r:
            return j0 * t / r
        if t > r + flat:
            return j0 * (flat + 2 * r - t) / r
        return j0

    def power_areal_W_m2(self, t_s: float) -> float:
        return self.current_density_A_m2(t_s) * self.voltage_kV * 1e3


# ---------------------------------------------------------- stopping power
def stopping_power(densities: dict, E_eV: float) -> float:
    """Collisional stopping power -dE/dx [eV/m] of a {gas: n_m3} mixture.

    Relativistic Bethe/Møller for electrons (see PIC_MCC_3D/beam3d for the
    formula's derivation notes); density-effect delta = 0 (negligible at
    ~atm for sub-MeV electrons); bremsstrahlung excluded (<1 % below ~1 MeV
    in these gases, and it does not deposit locally)."""
    E_eV = float(E_eV)
    if E_eV <= 0.0:
        return 0.0
    tau = E_eV / MEC2_EV
    beta2 = 1.0 - 1.0 / (1.0 + tau) ** 2
    if beta2 <= 0.0:
        return 0.0
    total = 0.0
    for gas, n_gas in densities.items():
        if n_gas <= 0.0:
            continue
        if gas not in _ICRU:
            raise KeyError(f"no ICRU stopping data for {gas!r}")
        I_eV, Z, _A = _ICRU[gas]
        n_e = n_gas * Z
        F = (1.0 - beta2
             + (tau ** 2 / 8.0 - (2.0 * tau + 1.0) * np.log(2.0))
             / (1.0 + tau) ** 2)
        L = np.log(tau ** 2 * (tau + 2.0)
                   / (2.0 * (I_eV / MEC2_EV) ** 2)) + F
        if L <= 0.0:
            continue
        total += (2.0 * np.pi * R_E ** 2 * MEC2_EV * n_e / beta2) * L
    return float(total)


def csda_range(densities: dict, E0_eV: float, e_min_eV: float = 1.0e3,
               n_steps: int = 400) -> float:
    """Continuous-slowing-down range [m] (integral dE/S from e_min to E0)."""
    E0_eV = float(E0_eV)
    if E0_eV <= e_min_eV:
        return 0.0
    E = np.linspace(e_min_eV, E0_eV, int(n_steps))
    S = np.array([stopping_power(densities, e) for e in E])
    good = S > 0.0
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(1.0 / S[good], E[good]))


def deposition_profile(spec: BeamSpec, cell_densities, dx: float) -> np.ndarray:
    """Fraction of beam energy deposited in each cell along the beam axis.

    ``cell_densities`` is a sequence of {gas: n_m3} dicts, one per cell IN
    MARCH ORDER (entry face first) — so a density gradient reshapes the
    profile (the march is local). Returns fractions summing to <= 1; the
    remainder exits the far face. Straight-ray, single-energy: no detour
    factor, no straggling, no backscatter — O(10 %) shape corrections,
    stated here rather than discovered later (beam3d's scope note)."""
    E = spec.energy_eV
    n_cells = len(cell_densities)
    dx_k = (np.full(n_cells, float(dx)) if np.isscalar(dx)
            else np.asarray(dx, dtype=float))
    out = np.zeros(n_cells)
    for k in range(n_cells):
        if E <= 1.0e3:
            out[k] += E / spec.energy_eV
            E = 0.0
            break
        dE = min(stopping_power(cell_densities[k], E) * dx_k[k], E)
        out[k] = dE / spec.energy_eV
        E -= dE
    if not spec.from_low:
        out = out[::-1].copy()
    return out


def equal_deposition_edges(spec: BeamSpec, densities: dict, L_m: float,
                           n_cells: int, n_march: int = 800) -> np.ndarray:
    """The MESH.txt refinement rung: axis cell edges that EQUIDISTRIBUTE
    the deposited beam energy — cells crowd where the material actually
    stops the beam (the Bragg peak), coarsen where little deposits.

    A fine CSDA march builds the cumulative deposition F(x); edges sit at
    the n-quantiles of F. Static (computed once per configuration), which
    is the honest v1 of "adaptive": the mesh adapts to the deposition of
    THIS beam in THIS material, not dynamically during the run. Returns
    n_cells + 1 monotone edges spanning [0, L]."""
    xs = np.linspace(0.0, L_m, int(n_march) + 1)
    fine = deposition_profile(spec, [densities] * int(n_march),
                              L_m / int(n_march))
    F = np.concatenate([[0.0], np.cumsum(fine)])
    if F[-1] <= 0.0:
        return np.linspace(0.0, L_m, int(n_cells) + 1)
    targets = np.linspace(0.0, F[-1], int(n_cells) + 1)
    # strictly-monotone F for interpolation (flat tails after full stop)
    F_mono = F + np.arange(F.size) * (1e-12 * max(F[-1], 1e-300))
    edges = np.interp(targets + np.arange(targets.size)
                      * (1e-12 * max(F[-1], 1e-300)), F_mono, xs)
    edges[0], edges[-1] = 0.0, L_m
    return np.maximum.accumulate(edges)


# ------------------------------------------- Spencer–Fano per-process yields
def spencer_fano_per_process(cs: CrossSectionSet, E0_eV: float,
                             n_grid: int = 1200,
                             prompt_floor_eV: float = 0.0):
    """Expected count of every inelastic process per source electron at E0.

    The Fowler/Spencer–Fano recursion of PIC_MCC_3D/beam3d, extended to a
    PER-PROCESS vector (the per-species split needs it): one low→high sweep
    on a uniform grid; an excitation of threshold I continues the electron
    at E−I; an ionization counts one pair and splits the remainder equally
    between primary and secondary at (E−I)/2. Elastic excluded (energy loss
    2m/M per event — the standard approximation).

    ``prompt_floor_eV`` (0.1.12): processes with thresholds BELOW this are
    excluded from the recursion. The recursion is TIME-FREE — it asks
    'which inelastic happens next' with no clock — so below the buffer
    gas's inelastic floor, where only trace-molecular channels stay open,
    it would funnel 100% of the flux through them however slow their
    absolute collision rate is (~ns per event at 0.1-0.5% F2: tens of
    events spanning the whole pulse, NOT prompt). Those channels belong
    to the rate-limited bulk chemistry the network already carries (one
    physical pathway, one owner — the 0.1.11 lesson at the energy level);
    with a floor set, the cascade terminates there and the sub-floor
    energy stays with the electron gas (the Platzman split, which the
    cited w_table convention always used). Pure noble-gas cascades have
    no sub-floor channels, so every W anchor is unchanged bit-for-bit.

    Returns (U (n_proc,), terminal_eV): expected process counts and the
    sub-threshold terminal energy per source electron. Energy identity
    E0 = sum(U_p * thr_p) + terminal holds to round-off by construction.
    """
    idx = [k for k, kind in enumerate(cs.kinds)
           if kind in ("excitation", "ionization")
           and float(cs.thresholds[k]) >= float(prompt_floor_eV)]
    if not idx:
        raise ValueError("no inelastic processes in this cross-section set")
    E0_eV = float(E0_eV)
    ceiling = cs.energy_ceiling_eV
    if E0_eV > ceiling:
        raise ValueError(
            f"source energy {E0_eV:.1f} eV is above the cross-section "
            f"ceiling {ceiling:.1f} eV — every sigma reads ZERO there and "
            "the cascade would silently produce nothing (PIC_MCC_3D "
            "Known-bad). Transport fast electrons by stopping power.")
    thr = cs.thresholds[idx]
    is_ion = np.array([cs.kinds[k] == "ionization" for k in idx])
    n_grid = int(n_grid)
    grid = np.linspace(0.0, E0_eV, n_grid)
    dE = grid[1] - grid[0]
    n_p = len(idx)
    U = np.zeros((n_grid, n_p))

    def gather(e):
        """U interpolated at energy e (uniform grid) -> (n_p,)."""
        if e <= 0.0:
            return np.zeros(n_p)
        pos = min(e / dE, n_grid - 1.000001)
        i0 = int(pos)
        f = pos - i0
        return U[i0] * (1.0 - f) + U[i0 + 1] * f

    for i in range(1, n_grid):
        E = grid[i]
        nu = cs.nu_matrix(np.array([E]))[idx, 0]
        open_ = thr <= E
        tot = float(nu[open_].sum())
        if tot <= 0.0 or not open_.any():
            U[i] = U[i - 1] * 0.0
            continue
        acc = np.zeros(n_p)
        for j in np.nonzero(open_)[0]:
            w = nu[j] / tot
            if w <= 0.0:
                continue
            unit = np.zeros(n_p)
            unit[j] = 1.0
            if is_ion[j]:
                acc += w * (unit + 2.0 * gather(0.5 * (E - thr[j])))
            else:
                acc += w * (unit + gather(E - thr[j]))
        U[i] = acc

    # scatter back to the FULL process list (zeros for elastic/attachment)
    # so callers can index with cs.kinds positions
    U_full = np.zeros(cs.n_proc)
    U_full[idx] = U[-1]
    terminal = E0_eV - float(np.dot(U[-1], thr))
    return U_full, terminal


def default_cascade_energy(cs: CrossSectionSet, margin: float = 0.55) -> float:
    """A source energy safely inside the tables (beam3d's 55 % rule: high
    enough for a real cascade, below where the tables running out bends the
    yield — see differential_W)."""
    return float(margin * cs.energy_ceiling_eV)


def differential_W(cs: CrossSectionSet, e_lo: float = None,
                   e_hi: float = None) -> float:
    """W as dE/dU_ion over a clean 25–55 %-of-ceiling window [eV] — the
    quantity the literature reports (beam3d's rule: the integral ratio
    E0/U_ion includes the sub-excitation bottom of the cascade and always
    sits above the asymptote; the table edge bends it back up)."""
    c = cs.energy_ceiling_eV
    e_lo = 0.25 * c if e_lo is None else float(e_lo)
    e_hi = 0.55 * c if e_hi is None else float(e_hi)
    idx = [k for k, kind in enumerate(cs.kinds) if kind == "ionization"]
    Ua, _ = spencer_fano_per_process(cs, e_lo)
    Ub, _ = spencer_fano_per_process(cs, e_hi)
    d_u = float(Ub[idx].sum() - Ua[idx].sum())
    if d_u <= 0.0:
        return float("inf")
    return (e_hi - e_lo) / d_u


# ------------------------------------------------------------- the yields
@dataclass
class BeamYields:
    """Per-ion-pair beam production, with the conserving identity attached.

    ``per_pair`` maps species -> number produced per ion pair; ``W_eV`` is
    the energy per pair; ``eps_sec_eV`` heats the electron gas;
    ``subfloor_eV`` (untracked sub-threshold channels — F2 vibrations —
    plus the fragment heat of prompt F2 dissociation, 0.1.25) goes to the
    gas. Identity (gated):

        W_eV = sum(per_pair[s] * E_form(s)) + subfloor_eV + eps_sec_eV
    """

    per_pair: dict
    W_eV: float
    eps_sec_eV: float
    subfloor_eV: float
    mode: str                   # 'spencer_fano' | 'w_table'
    cite: str

    def identity_residual_eV(self) -> float:
        spend = sum(n * e_form_of(s) for s, n in self.per_pair.items())
        return self.W_eV - spend - self.subfloor_eV - self.eps_sec_eV


def spencer_fano_yields(densities: dict, E0_eV: float = None) -> BeamYields:
    """First-principles per-species yields from the LXCat sets (see module
    docstring). ``densities``: {gas: n_m3} — composition sets the split.

    The prompt cascade ends at the NOBLE-gas inelastic floor (0.1.12):
    sub-floor trace-molecular channels (F2 electronic/vibrational) are
    rate-limited, not prompt, and the network's own e+F2 rows carry them —
    counting them here spent the Platzman sub-excitation energy at t=0
    (measured: eps_sec 1.81 eV/pair at 0.5% F2 vs ~5 for pure Ar — the
    bulk ran x4 cold and F- collapsed x400 vs the NRL benchmark,
    scorecard.electra_population_rows). Below the floor the secondary IS
    the bulk; its energy stays in the electron gas."""
    cs = CrossSectionSet(densities)
    E0 = default_cascade_energy(cs) if E0_eV is None else float(E0_eV)
    _NOBLE = ("He", "Ne", "Ar", "Kr", "Xe")
    noble_thr = [float(cs.thresholds[k]) for k, kind in enumerate(cs.kinds)
                 if kind in ("excitation", "ionization")
                 and cs.gases[k] in _NOBLE]
    floor = min(noble_thr) if noble_thr else 0.0
    U, terminal = spencer_fano_per_process(cs, E0, prompt_floor_eV=floor)
    ion_rows = [k for k, kind in enumerate(cs.kinds) if kind == "ionization"]
    u_ion = float(U[ion_rows].sum())
    if u_ion <= 0.0:
        raise ValueError("Spencer-Fano cascade produced no ionization")
    W = E0 / u_ion
    per_pair = {}
    subfloor = 0.0
    for k, kind in enumerate(cs.kinds):
        n_per_pair = float(U[k]) / u_ion
        if n_per_pair <= 0.0 or kind not in ("excitation", "ionization"):
            continue
        g = cs.gases[k]
        thr = float(cs.thresholds[k])
        target = None
        if kind == "ionization":
            target = BEAM_ION.get(g)
        elif g == "F2":
            if thr < _F2_ELEC_MIN:
                target = None                      # vibrational: subfloor
            else:
                # F2 electronic excitation IS prompt dissociation (the
                # states are repulsive — the network's 0.1.14 reroute;
                # only C1/H1 at 11.57/13.08 eV clear the noble prompt
                # floor). Booking them as bound F2(a3)/F2(A1) parked the
                # energy in species the Kannari-chemistry modes never
                # decay (found 2026-08-25; ~0.005% of drive at 0.1% F2,
                # growing with F2 fraction). Book 2 F + the excess over
                # D0 as fragment heat — mirroring reroute_repulsive_f2.
                per_pair["F"] = per_pair.get("F", 0.0) + 2.0 * n_per_pair
                subfloor += n_per_pair * (thr - _F2_D0_EV)
                continue
        else:
            split = _LO_HI_SPLIT.get(g)
            if split is not None:
                target = (BEAM_EXC.get(g) if thr < split
                          else BEAM_EXC_HI.get(g, BEAM_EXC.get(g)))
        if target is None:
            subfloor += n_per_pair * thr
            continue
        per_pair[target] = per_pair.get(target, 0.0) + n_per_pair
    # the exported identity uses REGISTRY formation energies (what the
    # species store will hold); the (small) LXCat-threshold vs E_form
    # mismatch lands inside eps_sec, keeping closure exact by construction
    spend = sum(n * e_form_of(s) for s, n in per_pair.items())
    eps_sec = W - spend - subfloor
    if eps_sec < 0.0:
        raise ValueError(
            f"Spencer-Fano ledger over-spent: eps_sec = {eps_sec:.3f} eV — "
            "formation energies and thresholds disagree beyond the identity's "
            "slack; do not proceed with a non-conserving beam split")
    return BeamYields(per_pair, W, eps_sec, subfloor, "spencer_fano",
                      "Fowler/Spencer-Fano recursion over the LXCat sets "
                      f"(E0 = {E0:.0f} eV, {cs.describe()})")


def w_table_yields(fractions: dict) -> BeamYields:
    """The cited W-table fallback (spectre's conserving Lee ledger).

    Per pair in gas g: 1 ion + f_exc excited atoms; the subexcitation
    remainder wbar = W − E_ion − f_exc·E_exc heats the electron gas.
    Mixture-averaged with the STOPPING share as the channel weight (the
    spectre convention)."""
    share = stopping_share(fractions)
    per_pair = {}
    W_mix = eps_mix = 0.0
    for g, fr in share.items():
        if fr <= 0.0 or g not in W_EV:
            continue
        pairs_g = fr / W_EV[g]                 # pairs per eV into gas g
        per_pair[BEAM_ION[g]] = per_pair.get(BEAM_ION[g], 0.0) + pairs_g
        exc = BEAM_EXC.get(g)
        if exc is not None and F_EXC.get(g, 0.0) > 0.0:
            per_pair[exc] = per_pair.get(exc, 0.0) + F_EXC[g] * pairs_g
    tot_pairs = sum(v for s, v in per_pair.items() if s.endswith("+"))
    if tot_pairs <= 0.0:
        raise ValueError("w_table_yields: no ionizing component in mixture")
    per_pair = {s: v / tot_pairs for s, v in per_pair.items()}
    W_mix = 1.0 / tot_pairs                    # eV per pair, mixture
    spend = sum(n * e_form_of(s) for s, n in per_pair.items())
    eps = W_mix - spend
    if eps < 0.0:
        raise ValueError("w_table ledger over-spent — check W/E_form tables")
    return BeamYields(per_pair, W_mix, eps, 0.0, "w_table", _W_TABLE_CITE)


def stopping_share(fractions: dict) -> dict:
    """Fraction of the stopping (≈ deposition) per gas, from the per-atom
    stopping cross sections at the beam energy class (Z-weighted Bethe —
    composition-only, so it needs no absolute density)."""
    contrib = {}
    for g, fr in fractions.items():
        if fr <= 0.0 or g not in _ICRU:
            continue
        contrib[g] = fr * stopping_power({g: 1.0e25}, 6.0e5)
    tot = sum(contrib.values())
    if tot <= 0.0:
        return {}
    return {g: c / tot for g, c in contrib.items()}


# ------------------------------------------------------------- the source
def beam_sources(yields: BeamYields, p_dep_W_m3):
    """Chemistry sources from local deposited power.

    Returns (production {species: rate m^-3 s^-1 (array-shaped like
    p_dep)}, electron-gas heating [W/m^3], gas-heating subfloor [W/m^3]).
    The three outputs exhaust p_dep exactly (the conserving identity)."""
    p = np.asarray(p_dep_W_m3, dtype=float)
    pairs = p / (yields.W_eV * EV)             # pairs m^-3 s^-1
    prod = {s: n * pairs for s, n in yields.per_pair.items()}
    return prod, pairs * yields.eps_sec_eV * EV, pairs * yields.subfloor_eV * EV
