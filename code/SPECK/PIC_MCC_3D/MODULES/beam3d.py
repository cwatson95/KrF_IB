"""MISSING MODULE 9 -- e-beam injection, and the preionization seed.

THE SPLIT THAT MAKES THIS TRACTABLE
-----------------------------------
A 600 kV primary is NOT a PIC particle in this package, and trying to make it
one is the first way an e-beam model goes wrong:

  * it is RELATIVISTIC -- gamma = 2.17 at 600 keV, v = 0.89 c -- so the
    non-relativistic leapfrog in ``species3d`` would be wrong by a factor of
    two in inertia, and
  * it is far above the LXCat tables' ceiling.  ``cross_sections`` tabulates
    to 2 keV; ``Process.sigma_at`` returns ZERO past the last tabulated point.
    A 600 keV macro-particle would therefore sample **zero cross section**,
    collide with nothing, and fly through the gas without stopping -- silently,
    with no error raised.

So the fast primary is transported as a RAY with continuous-slowing-down
(CSDA) energy loss, and only its low-energy PROGENY become PIC particles --
which is exactly the regime the LXCat data and the MCC are valid in.  This is
also what the repo's 0-D model does (``unified/ebeam.py``: stopping-power
deposition + secondary spectrum), so the two are comparable by construction.

STOPPING POWER IS FIRST-PRINCIPLES, NOT A TABLE LOOKUP.  :func:`stopping_power`
evaluates the relativistic Bethe/Moller collisional formula with ICRU mean
excitation energies.  Measured against the NIST ESTAR values that
``unified/ebeam.py`` carries, it agrees to 2-6 % at 600 keV -- and unlike a
single tabulated constant it has the right energy dependence, which is what
sets the deposition profile as the beam slows.

THE W-VALUE IS AN OUTPUT HERE, NOT AN INPUT
-------------------------------------------
The energy cost per ion pair (Ar ~26 eV, Kr ~24 eV, He ~41-46 eV) is the
strongest first-principles check available to this code, and it is only a
check if nothing imposes it.  Two independent routes are provided, both priced
from the SAME LXCat cross sections the MCC samples:

  * :func:`spencer_fano_yields` -- the deterministic Fowler/Spencer-Fano yield
    recursion (the method of ``unified/spencer_fano.py``);
  * :func:`cascade_yields` -- a Monte-Carlo cascade run through the actual
    :class:`~PIC_MCC_3D.MODULES.mcc3d.MCC3D` operator.

They must agree with each other and with the literature.  ``pair_mode="w"``
on :class:`BeamInjector` exists for cheap runs, and is labelled: with it, W is
an input and NO W-value claim may be made from that run.

ON THE SURFACE DISCHARGE (why it is not simulated here)
-------------------------------------------------------
:func:`seed_preionization` enters a preionized sheet as an INITIAL CONDITION
-- a seeded electron/ion population at a stated density and temperature -- and
does not model the surface flashover that produced it.  That is a deliberate
scope call, and it is the repo's existing precedent: PLAN (67)'s Argon_UI puts
"the no-arc surface-discharge sheet ... as an INITIAL CONDITION seeded on Ar+"
for the same reason.  Simulating a surface discharge needs dielectric surface
charging, field emission, photoemission and streamer formation -- i.e. missing
modules 4 and 5 plus physics no module here has -- and every one of those is a
place to introduce a bug that would quietly corrupt the beam answer.  Seeding
the result, stating the density, and letting the gates check what happens NEXT
is the honest and far less bug-prone path.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cross_sections import CrossSectionSet
from .gas import GasMixture
from .species3d import ParticleSpecies, Species
from .units import EV, K_B, M_E

# --- constants for the Bethe/Moller collisional stopping power --------------
R_E = 2.8179403262e-15          # classical electron radius [m]
MEC2_EV = 510998.95             # electron rest energy [eV]
N_A = 6.02214076e23             # 1/mol

# ICRU 37 mean excitation energies I [eV], atomic number Z, molar mass [g/mol].
# F2 is taken as 2 x atomic fluorine (I_F = 115 eV); the molecular correction is
# small next to the log dependence.
_ICRU = {
    "He": (41.8, 2, 4.002602),
    "Ar": (188.0, 18, 39.948),
    "Kr": (352.0, 36, 83.798),
    "F2": (115.0, 18, 37.9968),
    "Ne": (137.0, 10, 20.1797),
    "N2": (82.0, 14, 28.0134),
    "O2": (95.0, 16, 31.9988),
}

# Opal-Peterson-Beaty (1971) secondary-electron spectrum parameter E_bar [eV]:
# dsigma/deps ~ 1 / (1 + (eps/E_bar)^2).  Values are the canonical fits.
OPB_EBAR_EV = {"He": 15.8, "Ne": 19.9, "Ar": 10.0, "Kr": 9.6, "F2": 13.0,
               "N2": 13.0, "O2": 17.4}


# --------------------------------------------------------------- stopping power
def stopping_power(mixture: GasMixture, E_eV: float) -> float:
    """Collisional (electronic) stopping power -dE/dx [eV/m] of the mixture.

    Relativistic Bethe/Moller formula for electrons:

        -dE/dx = 2 pi r_e^2 m_e c^2 n_e / beta^2
                 * [ ln( tau^2 (tau+2) / (2 (I/m_e c^2)^2 ) ) + F(tau) - delta ]

    with tau = E / m_e c^2, F(tau) = 1 - beta^2
    + [tau^2/8 - (2 tau + 1) ln 2] / (tau+1)^2, and n_e the gas ELECTRON
    density (sum of n_i Z_i).

    The density-effect correction ``delta`` is set to zero: at ~1 atm it is
    negligible for the sub-MeV electrons this project uses.  Radiative
    (bremsstrahlung) loss is also excluded -- it is <1 % of the collisional
    term below ~1 MeV in these gases, and it does not deposit locally anyway.
    """
    E_eV = float(E_eV)
    if E_eV <= 0.0:
        return 0.0
    tau = E_eV / MEC2_EV
    beta2 = 1.0 - 1.0 / (1.0 + tau) ** 2
    if beta2 <= 0.0:
        return 0.0
    total = 0.0
    for gas, n_gas in mixture.densities.items():
        if n_gas <= 0.0:
            continue
        if gas not in _ICRU:
            raise KeyError(f"no ICRU stopping data for {gas!r}")
        I_eV, Z, _A = _ICRU[gas]
        n_e = n_gas * Z                                   # electrons / m^3
        F = (1.0 - beta2
             + (tau ** 2 / 8.0 - (2.0 * tau + 1.0) * np.log(2.0))
             / (1.0 + tau) ** 2)
        L = np.log(tau ** 2 * (tau + 2.0)
                   / (2.0 * (I_eV / MEC2_EV) ** 2)) + F
        if L <= 0.0:                    # below the formula's validity, no loss
            continue
        total += (2.0 * np.pi * R_E ** 2 * MEC2_EV * n_e / beta2) * L
    return float(total)


def csda_range(mixture: GasMixture, E0_eV: float, e_min_eV: float = 1.0e3,
               n_steps: int = 400) -> float:
    """Continuous-slowing-down range [m]: integral of dE / S(E) from e_min to E0.

    Integration stops at ``e_min_eV`` (1 keV by default) because the Bethe
    formula loses validity as E approaches the mean excitation energy -- and
    below ~1 keV the electron is inside the LXCat tables' range, where the MCC
    takes over anyway.  The neglected tail is a sub-percent contribution to
    the range of a 600 keV electron.
    """
    E0_eV = float(E0_eV)
    if E0_eV <= e_min_eV:
        return 0.0
    E = np.linspace(e_min_eV, E0_eV, int(n_steps))
    S = np.array([stopping_power(mixture, e) for e in E])
    good = S > 0.0
    return float(np.trapezoid(1.0 / S[good], E[good])) if hasattr(np, "trapezoid") \
        else float(np.trapz(1.0 / S[good], E[good]))


# ------------------------------------------------------------ beam definition
@dataclass
class BeamSpec:
    """A hibachi-style e-beam entering the cell through a foil.

    ``voltage_kV``      beam electron energy at the gas (post-foil)
    ``j_A_cm2``         current density at the gas
    ``pulse_ns``        flat-top duration
    ``rise_ns``         linear rise/fall applied at both ends
    ``foil_transmission`` fraction of the diode current reaching the gas
    ``axis``            injection direction, 0/1/2 -> x/y/z
    ``from_low``        True: enters at the low-coordinate face
    """

    voltage_kV: float = 600.0
    j_A_cm2: float = 20.0
    pulse_ns: float = 50.0
    rise_ns: float = 5.0
    foil_transmission: float = 0.75
    axis: int = 2
    from_low: bool = True

    @property
    def energy_eV(self) -> float:
        return 1.0e3 * float(self.voltage_kV)

    def current_density_A_m2(self) -> float:
        return 1.0e4 * float(self.j_A_cm2) * float(self.foil_transmission)

    def envelope(self, t_s: float) -> float:
        """Temporal shape in [0, 1]: linear rise, flat top, linear fall."""
        t_ns = t_s * 1.0e9
        r = max(float(self.rise_ns), 0.0)
        flat = float(self.pulse_ns)
        if t_ns <= 0.0:
            return 0.0
        if r > 0.0 and t_ns < r:
            return t_ns / r
        if t_ns <= r + flat:
            return 1.0
        if r > 0.0 and t_ns < r + flat + r:
            return max(0.0, 1.0 - (t_ns - r - flat) / r)
        return 0.0

    def power_areal_W_m2(self, t_s: float) -> float:
        """Beam power per unit area entering the gas [W/m^2] at time ``t_s``."""
        return (self.current_density_A_m2() * float(self.voltage_kV) * 1.0e3
                * self.envelope(t_s))


def deposition_profile(mixture: GasMixture, spec: BeamSpec, h: float,
                       n_cells: int) -> np.ndarray:
    """Fraction of the beam energy deposited in each cell along the axis.

    Marches one CSDA ray through ``n_cells`` cells of size ``h``, losing
    S(E)*h each cell until the electron is spent.  Returns a length-``n_cells``
    array summing to <= 1 (the remainder leaves through the far face).

    This is a straight-ray, single-energy treatment: no scattering-induced
    detour factor, no straggling, no backscatter.  For a beam whose range is
    long compared with the cell it is the deposition profile's shape that
    matters, and those corrections change the shape by O(10 %), not its
    character -- stated here rather than discovered later.
    """
    E = spec.energy_eV
    out = np.zeros(int(n_cells))
    if E <= 0.0:
        return out
    for k in range(int(n_cells)):
        if E <= 1.0e3:
            out[k] += E / spec.energy_eV       # dump the remainder locally
            E = 0.0
            break
        dE = min(stopping_power(mixture, E) * h, E)
        out[k] = dE / spec.energy_eV
        E -= dE
    if not spec.from_low:
        out = out[::-1]
    return out


# ------------------------------------------------- secondary-electron spectrum
def sample_opal_peterson(rng: np.random.Generator, n: int, gas: str,
                         eps_max_eV: float) -> np.ndarray:
    """Sample ``n`` secondary energies from the Opal-Peterson-Beaty spectrum.

    dN/deps ~ 1 / (1 + (eps/E_bar)^2) on [0, eps_max], sampled exactly by
    inverting the arctangent CDF (no rejection loop, no tail truncation bias
    beyond the stated eps_max).
    """
    e_bar = OPB_EBAR_EV.get(gas, 13.0)
    u = rng.random(int(n))
    return e_bar * np.tan(u * np.arctan(max(eps_max_eV, 0.0) / e_bar))


# ------------------------------------------------------------ the W-value: two ways
def _inelastic_index(cs: CrossSectionSet):
    return [k for k, kind in enumerate(cs.kinds)
            if kind in ("excitation", "ionization")]


def _check_within_tables(cs: CrossSectionSet, E0_eV: float) -> None:
    """Refuse a source energy the cross-section data cannot support.

    Above ``cs.energy_ceiling_eV`` every sigma reads zero, so a cascade there
    returns "no ionizations" rather than an error -- an answer of exactly
    zero that looks like physics. Fail loudly instead.
    """
    ceiling = cs.energy_ceiling_eV
    if E0_eV > ceiling:
        raise ValueError(
            f"source energy {E0_eV:.1f} eV is above this cross-section set's "
            f"validity ceiling {ceiling:.1f} eV -- every cross section is ZERO "
            "there, so the cascade would silently produce no collisions. Use "
            "beam3d.default_cascade_energy(cs), or transport the fast "
            "electron by stopping power (that is what BeamInjector does).")


LITERATURE_W_EV = {"He": 41.3, "Ne": 36.6, "Ar": 26.4, "Kr": 24.2, "N2": 34.8}


def differential_W(cs: CrossSectionSet, e_lo: float = None,
                   e_hi: float = None) -> float:
    """The W-value as dE / dU_ion over a clean energy window [eV].

    THIS, not ``E0 / U_ion(E0)``, is the quantity the literature reports.
    The integral ratio includes the bottom of the cascade, where the last few
    tens of eV are spent on sub-excitation losses that produce no ion pair, so
    it always sits ABOVE the asymptotic value -- 35 eV at E0 = 100 eV for
    argon, falling toward the high-energy limit as E0 grows.

    The window also has to avoid the TOP of the tables. Approaching
    ``energy_ceiling_eV`` the cross sections run off the end of their tabulated
    range and the yield saturates, which bends W back UP (argon: 28.4 eV at
    600 eV, 31.7 eV at 950 eV -- an artefact of the data ending, not physics).
    Defaults therefore take 25-55 % of the ceiling.
    """
    c = cs.energy_ceiling_eV
    e_lo = 0.25 * c if e_lo is None else float(e_lo)
    e_hi = 0.55 * c if e_hi is None else float(e_hi)
    a = spencer_fano_yields(cs, e_lo)
    b = spencer_fano_yields(cs, e_hi)
    d_u = b["U_ion"] - a["U_ion"]
    if d_u <= 0.0:
        return float("inf")
    return float((e_hi - e_lo) / d_u)


def default_cascade_energy(cs: CrossSectionSet, margin: float = 0.55) -> float:
    """A source energy safely inside the cross-section set's validity range.

    Defaults to 55 % of :attr:`CrossSectionSet.energy_ceiling_eV` -- high
    enough for a real cascade, below the region where the tables running out
    distorts the yield (see :func:`differential_W`). Starting a
    degradation cascade AT the ceiling is the classic silent failure: every
    cross section reads zero there, so the electron never collides and the
    yield comes back as exactly zero rather than as an error.
    """
    return float(margin * cs.energy_ceiling_eV)


def spencer_fano_yields(cs: CrossSectionSet, E0_eV: float = None,
                        n_grid: int = 1200) -> dict:
    """Deterministic Fowler / Spencer-Fano yield recursion.

    Returns ``{'U_ion', 'U_exc', 'W_eV'}`` -- the expected number of ionizations
    and excitations produced by one electron of energy ``E0``, and
    W = E0 / U_ion.

    Method (the one ``unified/spencer_fano.py`` documents): at energy E the
    relative probability of inelastic process p is w_p = n_p sigma_p(E) /
    sum_inel(E).  An excitation of threshold I continues the electron at E - I;
    an ionization of threshold I counts one pair and splits the remainder
    equally between a primary and a secondary, each at (E - I)/2.  Every
    daughter energy is below E, so one low->high sweep on a grid is exact.
    Elastic collisions are excluded: in these heavy gases they degrade energy
    by only 2m/M per event (the standard approximation).
    """
    idx = _inelastic_index(cs)
    if not idx:
        raise ValueError("no inelastic processes in this cross-section set")
    E0_eV = default_cascade_energy(cs) if E0_eV is None else float(E0_eV)
    _check_within_tables(cs, E0_eV)
    thr = cs.thresholds[idx]
    is_ion = np.array([cs.kinds[k] == "ionization" for k in idx])
    grid = np.linspace(0.0, float(E0_eV), int(n_grid))
    U_ion = np.zeros_like(grid)
    U_exc = np.zeros_like(grid)

    def interp(arr, e):
        return np.interp(e, grid, arr, left=0.0, right=arr[-1])

    for i in range(1, grid.size):
        E = grid[i]
        sig = cs.nu_matrix(np.array([E]))[idx, 0]     # nu ~ n sigma v, v cancels
        tot = sig.sum()
        if tot <= 0.0:
            continue
        w = sig / tot
        open_ = thr <= E
        if not open_.any():
            continue
        ui = ue = 0.0
        for j in np.nonzero(open_)[0]:
            if is_ion[j]:
                half = 0.5 * (E - thr[j])
                ui += w[j] * (1.0 + 2.0 * interp(U_ion, half))
                ue += w[j] * (2.0 * interp(U_exc, half))
            else:
                ui += w[j] * interp(U_ion, E - thr[j])
                ue += w[j] * (1.0 + interp(U_exc, E - thr[j]))
        norm = w[open_].sum()
        U_ion[i] = ui / norm if norm > 0 else 0.0
        U_exc[i] = ue / norm if norm > 0 else 0.0

    u_ion = float(U_ion[-1])
    return {"U_ion": u_ion, "U_exc": float(U_exc[-1]),
            "W_eV": float(E0_eV / u_ion) if u_ion > 0 else np.inf}


def cascade_yields(mcc, E0_eV: float = None, n_electrons: int = 200,
                   rng: np.random.Generator = None,
                   max_rounds: int = 4000) -> dict:
    """Monte-Carlo degradation cascade through the ACTUAL MCC operator.

    Launches ``n_electrons`` at ``E0`` into a field-free box and collides them
    (and every secondary they make) until nothing is left above the lowest
    inelastic threshold.  Returns the measured yields and
    W = total energy in / ionizations produced.

    This is the emergent-W measurement: nothing about the answer is imposed,
    and it exercises the same sampling the production loop uses.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    cs = mcc.cs
    E0_eV = default_cascade_energy(cs) if E0_eV is None else float(E0_eV)
    _check_within_tables(cs, E0_eV)
    idx = _inelastic_index(cs)
    ion_rows = [k for k, kind in enumerate(cs.kinds) if kind == "ionization"]
    exc_rows = [k for k, kind in enumerate(cs.kinds) if kind == "excitation"]
    # STOP WHEN NO FURTHER IONIZATION IS POSSIBLE, not when every electron has
    # been ground below the lowest EXCITATION threshold.  Once the fastest
    # electron is under the ionization threshold the ion yield -- and hence W,
    # the quantity being measured -- is final.  Waiting for full exhaustion
    # instead means chasing a straggler sitting just above an excitation
    # threshold with a tiny cross section, cooled elastically at 2m/M ~ 3e-5
    # per collision: >6000 rounds in argon for an answer that stopped changing
    # after a few hundred.  The cost is that U_exc is then a LOWER BOUND,
    # which is reported as ``exc_complete``.
    stop_thr = float(cs.thresholds[ion_rows].min()) if ion_rows else \
        float(cs.thresholds[idx].min())

    n_electrons = int(n_electrons)
    p = ParticleSpecies(Species.electron(), capacity=64 * n_electrons)
    v = np.sqrt(2.0 * float(E0_eV) * EV / M_E)
    dirs = rng.normal(size=(n_electrons, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    zeros = np.zeros(n_electrons)
    p.spawn(zeros, zeros, zeros, v * dirs[:, 0], v * dirs[:, 1],
            v * dirs[:, 2], np.ones(n_electrons))

    n_ion = n_exc = 0
    rounds = 0
    while rounds < max_rounds:
        if p.n == 0:
            break
        eps_max = float(p.energy_eV().max())
        if eps_max < stop_thr:
            break
        # ADAPT the step to the population as it cools.  A step sized once for
        # the source energy is far too short for the tail: the majorant is set
        # by the fastest electron, so once the cascade has ground down to a few
        # eV every round costs a full pass for a vanishing collision
        # probability.  Re-pricing dt against the CURRENT maximum energy cut
        # the argon cascade from >6000 rounds to a few hundred, with the same
        # answer (W agreed to 0.4 % either way).
        dt = 0.5 / max(cs.nu_majorant(eps_max * 1.1), 1.0)
        ev = mcc.collide_adaptive(p, dt, rng)
        n_ion += int(ev["tally"][ion_rows].sum())
        n_exc += int(ev["tally"][exc_rows].sum())
        rounds += 1
    e_in = n_electrons * float(E0_eV)
    return {"U_ion": n_ion / n_electrons, "U_exc": n_exc / n_electrons,
            "W_eV": e_in / n_ion if n_ion else np.inf,
            "rounds": rounds, "converged": rounds < max_rounds,
            "exc_complete": False,      # see stop_thr: excitation may continue
            "n_final": p.n}


# ------------------------------------------------------------- the injector
class BeamInjector:
    """Turns beam energy deposition into PIC particles, step by step.

    ``pair_mode``:
      * ``"kinetic"`` (default) -- the deposited energy is converted into
        secondary electrons sampled from the Opal-Peterson spectrum, and the
        MCC cascades them.  The ion-pair count, and therefore W, EMERGES.
      * ``"w"`` -- pairs are created directly at a supplied W.  Cheap, and W is
        then an INPUT: a run in this mode may not be used to support any
        W-value or efficiency claim.  ``w_value_eV`` must be given explicitly;
        there is deliberately no default, so the assumption cannot be made by
        accident.
    """

    def __init__(self, mixture: GasMixture, spec: BeamSpec,
                 electrons: ParticleSpecies, ions: ParticleSpecies = None,
                 pair_mode: str = "kinetic", w_value_eV: float = None,
                 macro_energy_eV: float = 5.0e4,
                 rng: np.random.Generator = None) -> None:
        if pair_mode not in ("kinetic", "w"):
            raise ValueError(f"unknown pair_mode {pair_mode!r}")
        if pair_mode == "w" and w_value_eV is None:
            raise ValueError(
                "pair_mode='w' makes the W-value an INPUT -- pass w_value_eV "
                "explicitly (and do not quote a W-value from that run)")
        self.mixture = mixture
        self.spec = spec
        self.electrons = electrons
        self.ions = ions
        self.pair_mode = pair_mode
        self.w_value_eV = w_value_eV
        self.macro_energy_eV = float(macro_energy_eV)
        self.rng = np.random.default_rng(0) if rng is None else rng
        self.energy_in_J = 0.0
        self.energy_deposited_J = 0.0
        self.pairs_created = 0.0
        self._profile = None
        self._gas_for_secondaries = max(mixture.fractions,
                                        key=mixture.fractions.get)

    def profile(self, h: float, n_cells: int) -> np.ndarray:
        if self._profile is None or self._profile.size != n_cells:
            self._profile = deposition_profile(self.mixture, self.spec, h,
                                               n_cells)
        return self._profile

    def inject(self, sim, dt: float = None) -> dict:
        """Deposit one step's beam energy into ``sim`` and create the particles.

        Returns a small record of what was injected.
        """
        dt = sim.dt if dt is None else dt
        ax = self.spec.axis
        n_cells = sim.shape[ax]
        h = sim.h
        area = float(np.prod([sim.shape[a] * h for a in range(3) if a != ax]))
        p_areal = self.spec.power_areal_W_m2(sim.time)
        if p_areal <= 0.0:
            return {"energy_in_J": 0.0, "energy_deposited_J": 0.0,
                    "n_macro": 0, "pairs": 0.0}
        e_in = p_areal * area * dt                        # J entering this step
        prof = self.profile(h, n_cells)
        e_cell = e_in * prof                              # J per slab
        e_dep = float(e_cell.sum())                       # J actually stopped
        self.energy_in_J += e_in
        self.energy_deposited_J += e_dep

        if self.pair_mode == "w":
            pairs = e_cell / (self.w_value_eV * EV)
        else:
            # kinetic: the deposit becomes secondary electrons; how many pairs
            # follow is for the MCC to decide, not for this module
            pairs = np.zeros_like(e_cell)

        n_made = 0
        for k in np.nonzero(e_cell > 0.0)[0]:
            n_made += self._seed_slab(sim, ax, k, e_cell[k],
                                      None if self.pair_mode == "kinetic"
                                      else pairs[k])
        self.pairs_created += float(pairs.sum())
        # NOTE the two energies are different and both are reported.  Only the
        # DEPOSITED fraction becomes particles; for a 600 keV beam whose range
        # is metres, a centimetre-scale cell stops well under a percent of what
        # enters.  Returning one number called "energy" invites exactly the
        # wrong comparison.
        return {"energy_in_J": float(e_in),
                "energy_deposited_J": e_dep,
                "n_macro": n_made, "pairs": float(pairs.sum())}

    def _seed_slab(self, sim, ax: int, k: int, e_J: float,
                   pairs) -> int:
        """Place secondaries (and ions, in 'w' mode) in one slab."""
        h = sim.h
        rng = self.rng
        if self.pair_mode == "kinetic":
            n_macro = max(int(np.ceil(e_J / (self.macro_energy_eV * EV))), 1)
            eps = sample_opal_peterson(rng, n_macro,
                                       self._gas_for_secondaries, 200.0)
            eps = np.maximum(eps, 1e-3)
            weight = (e_J / EV) / float(np.sum(eps))      # conserve the energy
            w = np.full(n_macro, weight)
        else:
            n_macro = max(int(np.ceil(pairs / 1.0e10)), 1) if pairs > 0 else 0
            if n_macro == 0:
                return 0
            kT = K_B * self.mixture.T_gas_K / EV
            eps = rng.exponential(kT, n_macro)
            w = np.full(n_macro, pairs / n_macro)

        pos = [rng.random(n_macro) * sim.shape[a] * h for a in range(3)]
        pos[ax] = (k + rng.random(n_macro)) * h
        speed = np.sqrt(2.0 * eps * EV / M_E)
        d = rng.normal(size=(n_macro, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        self.electrons.spawn(pos[0], pos[1], pos[2], speed * d[:, 0],
                             speed * d[:, 1], speed * d[:, 2], w)
        if self.ions is not None and self.pair_mode == "w":
            zeros = np.zeros(n_macro)
            self.ions.spawn(pos[0], pos[1], pos[2], zeros, zeros, zeros, w)
        return n_macro


# ------------------------------------------------------ the preionization seed
def seed_preionization(electrons: ParticleSpecies, ions: ParticleSpecies,
                       n_macro: int, lo, hi, n_e_m3: float, Te_eV: float,
                       rng: np.random.Generator, Ti_eV: float = None) -> float:
    """Seed a quasi-neutral preionized region -- the SURFACE-DISCHARGE STAND-IN.

    Places matched electron and ion populations over the box ``lo..hi`` at
    density ``n_e_m3``, Maxwellian at ``Te_eV`` (ions at ``Ti_eV``, defaulting
    to the gas temperature if a mixture temperature is not supplied).  Returns
    the net charge of what was seeded, which a caller should check is ~0.

    WHY THIS INSTEAD OF SIMULATING THE SURFACE DISCHARGE.  A surface (sliding)
    discharge is a dielectric-surface phenomenon: charge accumulation on the
    insulator, field enhancement at the triple point, field/photo-emission, and
    streamer propagation across the surface.  None of those have a module in
    this package -- they need missing modules 4 and 5 first, and then physics
    beyond them -- so simulating it here would mean inventing several unvalidated
    models at once and letting their errors compound into the answer everyone
    actually wants (the beam-pumped kinetics).

    The repo has already taken this decision once, for the same reason: PLAN
    (67)'s Argon_UI enters "the no-arc surface-discharge sheet ... as an INITIAL
    CONDITION seeded on Ar+".  What the seed density IS remains an input the
    caller must justify from the hardware -- this function makes the assumption
    explicit and checkable rather than buried.
    """
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    vol = float(np.prod(hi - lo))
    n_macro = int(n_macro)
    weight = n_e_m3 * vol / max(n_macro, 1)
    Ti_eV = 0.025 if Ti_eV is None else float(Ti_eV)

    pos = lo[None, :] + rng.random((n_macro, 3)) * (hi - lo)[None, :]
    vth_e = np.sqrt(max(Te_eV, 0.0) * EV / electrons.species.mass)
    ve = rng.normal(0.0, vth_e, size=(n_macro, 3))
    electrons.spawn(pos[:, 0], pos[:, 1], pos[:, 2], ve[:, 0], ve[:, 1],
                    ve[:, 2], np.full(n_macro, weight))

    # ions are co-located with the electrons: seeding them independently would
    # create a spurious space-charge field from sampling noise alone
    vth_i = np.sqrt(max(Ti_eV, 0.0) * EV / ions.species.mass)
    vi = rng.normal(0.0, vth_i, size=(n_macro, 3))
    ions.spawn(pos[:, 0], pos[:, 1], pos[:, 2], vi[:, 0], vi[:, 1], vi[:, 2],
               np.full(n_macro, weight))

    return float(electrons.species.charge * weight * n_macro
                 + ions.species.charge * weight * n_macro)
