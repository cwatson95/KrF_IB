"""The 0-D model: rate equations, electron energy, gas temperature, ledger.

State (densities scaled by the initial gas density N0 so the solver sees
O(1) numbers):

    n_s / N0          every species the active rows touch
    U_e / N0          electron energy [eV per gas particle]; Te = U_e/(1.5 n_e)
    Tg                gas temperature [K] (constant volume: a microsecond
                      pulse is ~45x shorter than the acoustic time of a 1-cm
                      cell, so N stays fixed and the pressure rises)
    ledger / N0       E_dep, elastic heat, captured-electron heat [eV]
    events / N0       time-integrated events of every row

Every energy channel of the report is rebuilt from the event counts times
the per-event energies of network.energy_routing, while the densities come
out of the ODE itself. The closure in report.accounting() therefore tests
the bookkeeping instead of restating it.

Energy input (the notes' 'Input' block): power_W_cm3 into the electron
gas for pulse_s, deposited uniformly and directly ("do not care how the
power is deposited"). It is the only energy input from v0_1 on. No drive is
modeled (v0's EM-wave option is not carried over).
Electrons. No electrons are injected. ne_mode='free', the basic
network's default, starts from ne_cm3 (a seed, the plasma starting
quasineutral with Kr+ = ne_cm3) and lets ionization, attachment and
recombination set n_e, quasineutral: n_e = cations - anions. The seed only
sets the first nanoseconds. ne_mode='clamped' holds n_e = ne_cm3 for the
whole run, as v0 does. It is the 'auto' choice for the notes' list as
written (network 'v0'), which has no electron loss, and it must never run
with attachment or recombination (charge would not be conserved).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import solve_ivp

from . import boltzmann, xsec
from .network import (BANDS, FREE_BANDS, ION_K0, LADDER_KANNARI, bands_read,
                      build, energy_routing)
from .species import ELECTRON, SPECIES, e_form

KB = 1.380649e-23             # J/K
QE = 1.602176634e-19          # J/eV
KB_EV = KB / QE               # eV/K
ATM_PA = 101325.0
AMAGAT_CM3 = 2.6867811e19     # Loschmidt, 273.15 K and 1 atm
EPS0 = 8.8541878128e-12

BASE_SPECIES = ("Kr", "F2", "F", "Kr*", "Kr**", "Kr+", "KrF(B)", "Kr2F*")
# v0_4: E_em is the energy the EM wave deposits (the second input), E_wall
# the energy the charges carry to the electrodes (gap_cm), E_ion_heat the
# ambipolar field's work on the ions, which heats the gas
LEDGER = ("E_dep", "heat_elastic", "heat_capture", "E_em", "E_wall",
          "E_ion_heat", "E_wall_th")
# E_wall_th: the thermal energy (cv k Tg per particle) the charges lost to
# the walls carry out of the gas, which leaves Tg unchanged and is booked
# to the walls apart from their formation energy
WALL_TE_MAX = 30.0            # eV: the electron energy per cation that
#                               drives the ambipolar flow, U_e/(1.5 sum n+),
#                               is capped here. Te = U_e/(1.5 n_e) itself
#                               must not enter the wall term: in the
#                               attaching afterglow n_e is a small
#                               difference of two large ion densities, and
#                               the term's Jacobian explodes and stalls the
#                               stiff solver (2026-10-08)
TE_FLOOR_FIT = 0.026          # te_fit rows are held above this Te [eV]
# the energy a lost electron carries to the wall, in units of Te: the mean
# energy of the Maxwellian flux through a surface (2 Te), and the ambipolar
# potential's work on the ion it leaves with, which the ion gives the gas
# in its collisions (of order Te; 1Dv0_3 found the field's work on the ions
# at 0.2-0.5 % of the input at the bands' edges). Both are the closure of
# a 0-D box, stated here and reported, never fitted
WALL_E_ELECTRON_TE = 2.0
WALL_E_AMBIPOLAR_TE = 1.0
TE_FLOOR_STEPWISE = 0.026     # REBUILD chemistry._TE_MIN_EV
TE_FLOOR_POW = 0.02           # te_pow rows are held below this Te [eV]
NE_FLOOR = 1.0                # cm^-3


@dataclass(frozen=True)
class Params:
    pressure_atm: float = 1.0          # initial fill pressure (0.1 - 2 atm)
    x_f2: float = 0.005                # F2 mole fraction, rest Kr
    tg0_K: float = 300.0
    power_W_cm3: float = 3.0e5         # 100 - 500 kW/cm^3 in the notes
    pulse_s: float = 1.0e-6            # 'a microsecond pulse'
    t_end_s: float = 2.0e-6            # pulse + afterglow (Kr2F* tau 220 ns)
    ne_cm3: float = 1.0e15             # 1e14 - 1e16 in the notes
    ne_mode: str = "auto"              # 'free' | 'clamped' (auto: basic
                                       # free, v0 clamped)
    tau_sp_s: float = 1.0 / 1.43e8     # Kannari's KrF(B) radiative lifetime
    # the trimer's gas-temperature exponent (G5, G5C). None keeps the
    # network's, Shui's calculation (-3.32, computed 300-800 K). It is
    # measured only to 425 K, and pure-Kr laser efficiencies suggest -0.5
    # (Salesky & Kimura 1985-86), so a number here spans that band for a
    # sensitivity run. It never sets the model's own value.
    trimer_tg_pow: float = None
    # v0_2 (2026-09-29): the rates the basic network carries as a band
    # (network.BANDS), as ((name, value), ...). A band left out keeps its
    # default, and the unmeasured ones default to 0. For sensitivity runs
    # (run_v0_3.py --bands), never to match a device.
    rates: tuple = ()
    # v0_4 (2026-10-09): a sensitivity scaling of named rows, ((row id,
    # factor), ...), multiplying the row's rate as compiled (every kind
    # but the elastic row, which it refuses; under the kinetic EEDF the
    # rows the kinetic table carries take it into their cross section or
    # recombination template, 2026-10-14; K76 and K76r take one factor).
    # For propagating a source's quoted error through the model (the
    # accuracy overview), never to match a device. An unknown id is an
    # error
    row_scale: tuple = ()
    network: str = "basic"             # 'basic' (KrF(B,C) set) or 'v0'
    groups: tuple = ()                 # optional groups ('all' = every one)
    # the electron closure. 'auto' gives the basic network the kinetic EEDF
    # (v0_3) and the exact captured energy, and keeps v0 on the notes'
    # Maxwellian with (3/2)Te per captured electron. 'two_t' is v0_2's
    # closure, the Maxwellian with the tabulated tail (2026-09-28, step 3)
    eedf: str = "auto"                 # 'kinetic' | 'two_t' | 'maxwellian'
    capture: str = "auto"              # 'exact' | 'three_halves'
    # v0_3 (2026-09-30): the kinetic EEDF in the loop. Every electron rate,
    # the elastic loss and the captured energies come from the two-term
    # Boltzmann solution with e-e collisions at LIGHT's mean energy
    # (boltzmann.py), renewed when Te has moved 0.12 in ln Te, when n_e or
    # the F2 fraction has moved kin_step (relative), and at least every
    # kin_max_dt_s. At the reference point kin_step 0.05, 5 ns or 1200
    # cells move no band by more than 2e-4 of itself (2026-09-30)
    kin_cells: int = 600               # energy cells of the kinetic grid
    kin_emax_eV: float = 30.0          # its first top; it adapts
    kin_step: float = 0.10
    kin_max_dt_s: float = 10e-9
    # after the pulse the interval grows to this share of the time since
    # the pulse ended (the drift triggers watch every kinetic input)
    kin_afterglow_frac: float = 0.1
    # while the drive is on, and for this long from the start, the kinetic
    # EEDF is advanced in time, a backward-Euler step from the last
    # update's (boltzmann.Solver.step), instead of solved steady. The
    # default, inf, is the whole pulse, and 0 is off (2026-10-01: the
    # heating transient of the first nanosecond, which no steady solution
    # holds). After the pulse the EEDF is steady, as before
    kin_memory_s: float = math.inf
    kin_tol: float = 1e-9
    init_extra_cm3: tuple = ()         # ((species, density), ...) added at t=0
    # v0_4 (2026-10-08). The ladder preset: 'off' leaves every ladder
    # band at its default (the unmeasured ones 0), 'kannari' sets each to
    # Kannari's printed value (network.LADDER_KANNARI) unless `rates`
    # names it. A sensitivity configuration, never a fit.
    ladder: str = "off"
    # the electrode gap [cm] of a 0-D box between grounded plates. None:
    # no walls (the box of v0 to v0_3). With a gap the cations drift and
    # diffuse to the walls at the ambipolar rate D_a / Lambda^2, Lambda =
    # gap/pi, D_a = mu_i k (Te + Tg)/e with the measured reduced
    # mobilities (network.ION_K0), F- at its own free diffusion, and the
    # charges' energy goes to the walls (E_wall) and the gas (E_ion_heat)
    gap_cm: float = None
    # the EM wave (the notes' 'EM wave injection'): the rms field inside
    # the plasma [V/cm] and the frequency [Hz], on during the pulse. Its
    # absorbed power is an output, E_em in the ledger: the Lorentz-gas
    # kernel on the electron-neutral collisions plus Turnbull et al.'s
    # electron-ion inverse bremsstrahlung (ib.py), the latter times the
    # band ib_scale (0.9-1.1, their model's stated error)
    em_E_V_cm: float = 0.0
    em_f_Hz: float = 95e9
    rtol: float = 1e-8
    n_save: int = 801
    method: str = "BDF"

    def replace(self, **kw) -> "Params":
        return replace(self, **kw)


@dataclass
class Compiled:
    row: object
    j: int                 # row index (event-counter slot)
    heavy: tuple           # species indices, with multiplicity
    ne_pow: int
    stoich: tuple          # ((species index, coefficient), ...)
    route: dict
    kind: str
    rate: object           # constant, RateTable, or tuple of parameters
    cap: object = None     # captured energy E(Te) [eV]; None = (3/2)Te
    factor: float = 1.0    # the row's frac bands (v0_4), multiplies k


class Model:
    def __init__(self, p: Params):
        basic = p.network == "basic"
        self.ne_mode = (("free" if basic else "clamped")
                        if p.ne_mode == "auto" else p.ne_mode)
        if self.ne_mode not in ("clamped", "free"):
            raise ValueError(f"ne_mode {p.ne_mode!r}")
        self.eedf = (("kinetic" if basic else "maxwellian")
                     if p.eedf == "auto" else p.eedf)
        self.capture = (("exact" if basic else "three_halves")
                        if p.capture == "auto" else p.capture)
        if self.eedf not in ("kinetic", "two_t", "maxwellian"):
            raise ValueError(f"eedf {p.eedf!r}")
        if self.eedf == "kinetic" and self.capture != "exact":
            raise ValueError("eedf='kinetic' books the captured energy from "
                             "its own EEDF: it needs capture='exact'")
        if self.capture not in ("exact", "three_halves"):
            raise ValueError(f"capture {p.capture!r}")
        self.tail = (xsec.tail_ratio_table() if self.eedf == "two_t"
                     else None)
        self.p = p
        self.rows = build(p.groups, p.network)
        self.rates = dict(p.rates)
        self._row_scale = dict(p.row_scale)
        if any(not v >= 0.0 for v in self._row_scale.values()):
            raise ValueError(f"row_scale factors must be >= 0: {p.row_scale}")
        if p.ladder not in ("off", "kannari"):
            raise ValueError(f"ladder {p.ladder!r}: 'off' or 'kannari'")
        if p.ladder == "kannari":
            for name, v in LADDER_KANNARI.items():
                self.rates.setdefault(name, v)
        read = bands_read(self.rows)
        if set(self.rates) - read - set(FREE_BANDS):
            raise ValueError(f"rates {sorted(set(self.rates) - read)}: no "
                             f"row of this network reads them (bands here: "
                             f"{sorted(read)}; the trimer law is "
                             f"Params.trimer_tg_pow)")
        if any(not v >= 0.0 for v in self.rates.values()):
            raise ValueError(f"rates must be >= 0: {p.rates}")
        for name in ("ii_yield", "ii_f_c"):
            if not 0.0 <= self.rates.get(name, BANDS[name].default) <= 1.0:
                raise ValueError(f"{name} is a fraction: {self.rates[name]}")
        if p.gap_cm is not None and not p.gap_cm > 0.0:
            raise ValueError(f"gap_cm {p.gap_cm}: positive, or None")
        if p.em_E_V_cm < 0.0 or (p.em_E_V_cm > 0.0 and not p.em_f_Hz > 0.0):
            raise ValueError("the EM wave needs em_E_V_cm >= 0 and "
                             "em_f_Hz > 0")
        self.ib_scale = self.rates.get("ib_scale", BANDS["ib_scale"].default)
        species = list(BASE_SPECIES)
        for r in self.rows:
            for s in r.reactants + r.products:
                if s != ELECTRON and s not in species:
                    species.append(s)
        for s, _ in p.init_extra_cm3:
            if s not in species:
                species.append(s)
        self.species = species
        self.idx = {s: i for i, s in enumerate(species)}
        self.ns = len(species)
        self.i_ue = self.ns
        self.i_tg = self.ns + 1
        self.i_led = self.ns + 2
        self.i_ev = self.i_led + len(LEDGER)
        self.size = self.i_ev + len(self.rows)
        self.cations = [self.idx[s] for s in species if SPECIES[s].charge > 0]
        self.anions = [self.idx[s] for s in species if SPECIES[s].charge < 0]
        self.cv = np.array([SPECIES[s].cv for s in species])
        self.eform = np.array([e_form(s) for s in species])
        self.N0 = p.pressure_atm * ATM_PA / (KB * p.tg0_K) * 1e-6    # cm^-3
        self.kel = xsec.elastic_loss_table("Kr")
        # v0_4: the ions that drift to the walls, with their reduced
        # mobilities (every charged species must have one)
        self.wall = []
        if p.gap_cm is not None:
            for s in species:
                q = SPECIES[s].charge
                if q and s not in ION_K0:
                    raise ValueError(f"{s}: no reduced mobility in ION_K0")
                if q:
                    self.wall.append((self.idx[s], q, ION_K0[s]))
        self._em = None
        if p.em_E_V_cm > 0.0:
            from . import ib
            self._em = ib
        unknown = set(self._row_scale) - {r.rid for r in self.rows}
        if unknown:
            raise ValueError(f"row_scale: no row {sorted(unknown)} in this "
                             "network")
        # the elastic exchange is the momentum transfer of the whole EEDF
        # (both closures read the gas's elastic power), not a rate to scale
        elastic = {r.rid for r in self.rows if r.rate[0] == "elastic"}
        if set(self._row_scale) & elastic:
            raise ValueError(f"row_scale: the elastic row(s) "
                             f"{sorted(set(self._row_scale) & elastic)} "
                             "cannot be scaled")
        # K76 and K76r are one cross section and its detailed-balance
        # reverse: the kinetic EEDF builds the reverse from the forward's
        # scaled cross section, so a factor on one alone would mean one
        # thing on each closure (2026-10-14). Give both the same factor
        rids = {r.rid for r in self.rows}
        if {"K76", "K76r"} <= rids and (self._row_scale.get("K76", 1.0)
                                        != self._row_scale.get("K76r", 1.0)):
            raise ValueError("row_scale: K76 and K76r are one cross section "
                             "and its reverse; give them the same factor")
        self.compiled = [self._compile(j, r) for j, r in enumerate(self.rows)]
        # the kinetic EEDF's table of rates, captured energies and elastic
        # power over ln Te (None for the other closures)
        self._tab = None
        self.kin = None
        self.kin_log = []
        if self.eedf == "kinetic":
            self._kinetic_setup()

    # ------------------------------------------------------------ kinetic
    def _kinetic_setup(self):
        """The processes the kinetic equation carries: exactly the electron
        rows the network books, so its energy balance and the ledger close
        the same budget."""
        procs, self._kin_xsec = [], []
        sup = None
        # A row's factor (row_scale, and any band a row reads through its
        # frac) goes into its process here, as K76's band goes into its
        # cross section below, so the kinetic rate the run reads and the
        # collision operator carry it, as the closure's own rate (_k) does.
        # Until 2026-10-14 the table's rate was read without it, and
        # row_scale moved none of these rows under the kinetic EEDF (found
        # in LIGHT_Kr v0_1, 2026-10-13)
        for c in self.compiled:
            r = c.row
            if c.kind == "xsec":
                p = xsec.find(*r.rate[1:4])
                if c.factor != 1.0:
                    p = replace(p, sigma=p.sigma * c.factor)
                procs.append(boltzmann.Process(r.rate[1], p.kind,
                                               p.threshold, p, (r.rid,)))
                self._kin_xsec.append((c.j, p))
            elif c.kind == "superelastic":
                if sup is not None:
                    raise ValueError("more than one superelastic row")
                p = xsec.find(r.rate[1], "excitation", r.rate[2])
                # the Klein-Rosseland rate is linear in g_ratio
                sup = (p, float(r.rate[3]) * c.factor, c.j)
        self._kin_sup = sup
        # Penning ionization: rows that make an electron from heavy
        # reactants and give it their energy
        self._kin_pen = [c for c in self.compiled
                         if c.row.n_electrons_in == 0 and c.row.ne_delta > 0
                         and c.route["to_e"] > 0.0]
        # recombination: te_pow rows, sigma v = C eps^n with the Maxwellian
        # average of that template exactly the row's law k0 (Te/te0)^n
        self._kin_rec = []
        for c in self.compiled:
            if c.kind == "te_pow":
                k0, te0, n_pow = c.rate
                cons = (k0 * te0 ** (-n_pow) * math.gamma(1.5)
                        / math.gamma(1.5 + n_pow) * 1e-6       # m^3/s eV^-n
                        * c.factor)
                ion = [s for s in c.row.heavy]
                self._kin_rec.append((c.j, cons, float(n_pow),
                                      self.idx[ion[0]]))
        if len({n for _, _, n, _ in self._kin_rec}) > 1:
            raise ValueError("recombination rows with different Te laws")
        self._kin_att = [(c.j, p) for c in self.compiled for j, p in
                         self._kin_xsec if j == c.j and p.kind == "attachment"
                         and c.route["capture"]]
        # the elastic energy exchange of the gases the network books it for
        # (G6, Kr): F2's elastic collisions stay in its momentum transfer
        elastic = tuple(c.row.rate[1] for c in self.compiled
                        if c.kind == "elastic")
        if elastic != ("Kr",):
            raise ValueError(f"elastic rows {elastic}: the table's elastic "
                             "power is Kr's alone")
        # v0_4: the 5s <-> 5p electron collisions (K76, K76r) as a
        # population process on Kr*. Kannari's A Te^0.5 exp(-C/Te) is the
        # Maxwellian rate of sigma(eps) = sigma0 (1 - C/eps) above C, with
        # sigma0 = A / (gamma 2/sqrt pi), exactly, so the pseudo cross
        # section reproduces the printed law under a Maxwellian and gives
        # the kinetic EEDF the collision's shape. The energy the electron
        # loses is the ledger's gap (the Process's threshold), and the
        # cross section is 0 below C, so the reverse by detailed balance
        # is LIGHT's K76r. Without it the kinetic energy balance missed
        # LIGHT's by the ladder's share, 3 to 5 % (2026-10-08)
        self._kin_pop = None
        pops = []
        byrid = {c.row.rid: c for c in self.compiled}
        if ("K76" in byrid and "K76r" in byrid
                and byrid["K76"].rate[0] * byrid["K76"].factor > 0.0):
            a_cm3, _, c_ev = byrid["K76"].rate
            # the row's factor (row_scale) on the pseudo cross section, as
            # the closure's own rate (_k) carries it; the reverse follows
            # by detailed balance (2026-10-14)
            a_cm3 *= byrid["K76"].factor
            gap = e_form("Kr**") - e_form("Kr*")
            gamma_cm = math.sqrt(2.0 * QE / boltzmann.ME) * 100.0
            sigma0 = a_cm3 / (gamma_cm * 2.0 / math.sqrt(math.pi))
            e_grid = np.concatenate([[0.0, c_ev], np.geomspace(
                c_ev * (1.0 + 1e-9), 1000.0, 400)])
            sig = np.where(e_grid > c_ev, sigma0 * (1.0 - c_ev / np.maximum(
                e_grid, c_ev)), 0.0)
            pseudo = xsec.Process("excitation", "Kr*", "Kr**(5p)", gap,
                                  e_grid, sig)
            g_ratio = float(byrid["K76r"].row.rate[5])
            pops.append(boltzmann.PopProcess("Kr*", pseudo, g_ratio,
                                             ("K76", "K76r"), upper="Kr**"))
            self._kin_pop = (byrid["K76"].j, byrid["K76r"].j, "Kr*")
        self.kin = boltzmann.Solver(
            procs, n_cells=self.p.kin_cells, emax=self.p.kin_emax_eV,
            superelastic=None if sup is None else sup[:2],
            tol=self.p.kin_tol, elastic=elastic, pops=pops)
        self._kin_t = None
        self._kin_sol = None
        self._kin_pts = None
        self._kin_ref = None

    # the inputs of the kinetic state that renew the EEDF when they drift
    # (_kin_observed), each with the floor below which it cannot matter:
    # Te, n_e, the F2 fraction, Tg, the Penning and recombination
    # frequencies per electron (1e4 /s, nothing next to attachment and
    # energy relaxation) and the Kr* fraction of a superelastic row
    # the Kr* and Kr** fractions' floors are 1e-7 (v0_4): below it their
    # collisions move nothing the table need follow, and a floor of 1e-9
    # had the trigger fire through the first nanoseconds, where Kr** grows
    # from nothing by orders of magnitude a step, until scipy's event
    # root-finder met a bracket without a sign change (2026-10-09)
    KIN_FLOORS = (0.0, 0.0, 1e-7, 0.0, 1e4, 1e4, 1e-7, 1e-7)

    def _kin_observed(self, y):
        """Every input of the kinetic state a state vector sets: Te, n_e,
        the F2 fraction, Tg, the Penning frequency per electron, the
        recombination frequency per electron (its template at Te), the
        Kr* fraction (0 without a superelastic row or the 5s <-> 5p pair)
        and the Kr** fraction (v0_4, 0 without the pair)."""
        n = np.maximum(y[:self.ns], 0.0) * self.N0
        ne = self.electron_density(n)
        te = max(y[self.i_ue] * self.N0 / (1.5 * ne), 1e-3)
        tg = y[self.i_tg]
        n_kr, n_f2 = n[self.idx["Kr"]], n[self.idx["F2"]]
        xf2 = n_f2 / max(n_kr + n_f2, 1.0)
        pen = 0.0
        for c in self._kin_pen:
            flux = self._rate(c, te, tg, n, ne)
            for i in c.heavy:
                flux *= n[i]
            pen += flux / ne
        rec = sum(cons * n[i_ion] * 1e6 * max(te, TE_FLOOR_POW) ** n_pow
                  for _, cons, n_pow, i_ion in self._kin_rec)
        kr_star = (n[self.idx["Kr*"]] / max(n_kr, 1.0)
                   if (self._kin_sup is not None or self._kin_pop is not None)
                   else 0.0)
        kr_pp = (n[self.idx["Kr**"]] / max(n_kr, 1.0)
                 if self._kin_pop is not None else 0.0)
        return te, ne, xf2, tg, pen, rec, kr_star, kr_pp

    def kinetic_injection(self, t, y, on):
        """Electrons a drive injects with a spectrum, as ((rate [cm^-3
        s^-1], e_top [eV], mean [eV]), ...), linear on [0, e_top]. The
        given power heats the electrons as a field and injects none. A
        beam drive returns its secondary electrons (BeamModel in
        LIGHT/docs/survey_2026-09-27/scripts/ebeam_screen.py)."""
        return ()

    def _kinetic_state(self, y, t=0.0, on=False):
        N0 = self.N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        ne = self.electron_density(n)
        te = max(y[self.i_ue] * N0 / (1.5 * ne), 1e-3)
        tg = y[self.i_tg]
        n_kr, n_f2 = n[self.idx["Kr"]], n[self.idx["F2"]]
        n_col = n_kr + n_f2
        pen = []
        for c in self._kin_pen:
            flux = self._rate(c, te, tg, n, ne)
            for i in c.heavy:
                flux *= n[i]
            pen.append((flux / ne, float(c.route["to_e"])))
        rec, rec_pow = 0.0, -0.5
        for j, cons, n_pow, i_ion in self._kin_rec:
            rec += cons * n[i_ion] * 1e6
            rec_pow = n_pow
        kr_star = 0.0
        if self._kin_sup is not None and n_kr > 0.0:
            kr_star = n[self.idx["Kr*"]] / n_kr
        inj = tuple((rate / ne, e_top, mean) for rate, e_top, mean
                    in self.kinetic_injection(t, y, on) if rate > 0.0)
        # v0_4: the cations' Coulomb collisions, and the EM wave while the
        # drive is on (its rms reduced field in Td)
        n_ion = sum(n[i] for i in self.cations) * 1e6
        omega = en_em = 0.0
        if self._em is not None and on:
            omega = 2.0 * math.pi * self.p.em_f_Hz
            en_em = self.p.em_E_V_cm * 100.0 / (n_col * 1e6) / boltzmann.TD
        pops = {}
        if self._kin_pop is not None and n_col > 0.0:
            pops["Kr*"] = n[self.idx["Kr*"]] / n_col
            pops["Kr**"] = n[self.idx["Kr**"]] / n_col
        s = boltzmann.State(x={"Kr": n_kr / n_col, "F2": n_f2 / n_col},
                            tg=tg, n_col=n_col * 1e6, ne=ne * 1e6,
                            pen=tuple(pen), rec=rec, rec_pow=rec_pow,
                            kr_star=kr_star, inj=inj, n_ion=n_ion,
                            omega=omega, en_em_td=en_em, pops=pops)
        return s, te, ne

    # the pulse's table: ln(Te / Te now) at its points, and the drift of Te
    # that renews it (inside the table, so the RHS never reads past its ends)
    KIN_TABLE = (-0.16, -0.08, 0.0, 0.08, 0.16)
    KIN_TABLE_REACH = 0.12
    KIN_NE_MIN = 1e3          # cm^-3: below it the EEDF is not solved

    def _kin_quantities(self, s, sol):
        """One EEDF's kinetic rates [cm^3/s] per row, captured energies [eV]
        per row, and elastic power per electron per Kr atom [eV cm^3/s], as
        the solver exchanges it (Solver.elastic_power)."""
        g, F = sol.grid, sol.F
        k, cap = {}, {}
        for j, p in self._kin_xsec:
            k[j] = g.rate(p, F) * 1e6
        if self._kin_sup is not None:
            p, g_ratio, j = self._kin_sup
            k[j] = g.superelastic_rate(p, g_ratio, F) * 1e6
        if self._kin_pop is not None:                     # v0_4
            j_f, j_r, name = self._kin_pop
            kf, kr = g.pop_rates(name, F)
            k[j_f], k[j_r] = kf * 1e6, kr * 1e6
        for j, cons, n_pow, i_ion in self._kin_rec:
            w = g.ec ** (n_pow + 0.5) * F
            k[j] = cons * float(w.sum()) * g.du * 1e6
            cap[j] = float((g.ec * w).sum() / w.sum())
        for j, p in self._kin_att:
            w = g.ec * boltzmann.sigma(p, g.ec) * F
            if w.sum() > 0.0:
                cap[j] = float((g.ec * w).sum() / w.sum())
        # the EM wave's power per electron, its electron-ion share scaled
        # by the ib_scale band (Turnbull's model error)
        pem = self.kin.em_power(s, sol)
        if pem > 0.0 and self.ib_scale != 1.0:
            atoms = self.kin.em_power(s, sol, ions=False)
            pem = atoms + self.ib_scale * (pem - atoms)
        return (k, cap, self.kin.elastic_power(s, sol, "Kr") * 1e6, pem)

    def _kinetic_update(self, t, y, on):
        """Solve the kinetic EEDF at state y and give the RHS its rates.

        The EEDF LIGHT reads has LIGHT's mean energy, (3/2) Te: the steady
        solution in which a field supplies whatever holds the electrons at
        that mean (the input power in the pulse, the energy LIGHT's
        electrons still carry after it). That is the local mean energy
        approximation of fluid models, with LIGHT's own energy equation
        choosing the mean. Where even the no-field solution is hotter (the
        Penning electrons alone heat it past Te), the field is 0 and the
        solution is the no-field one. The table holds the kinetic rates,
        captured energies and elastic power at five mean energies around
        LIGHT's (Te e^r for r in KIN_TABLE), read in ln Te, so between
        solutions LIGHT's energy equation sees the kinetic Te dependence.
        The kinetic mean always equals LIGHT's where the table is read, so
        no rate is carried from one mean energy to another (factors over a
        Maxwellian, ML's method, did that after the pulse until
        2026-09-30, and where the two means parted they drove LIGHT's Te
        down by a factor 2.5 in 0.1 ns)."""
        s, te, ne = self._kinetic_state(y, t, on)
        warm = self._kin_sol
        self._tab = None
        if ne < self.KIN_NE_MIN:
            # no electrons to speak of (a decay run's floor): their EEDF
            # cannot matter, and the Maxwellian closure stands in
            self._kin_ref = self._kin_observed(y)
            self.kin_log.append({"t": t, "on": on, "te": te, "skipped": True})
            return
        sols = {}
        memory = on and t < self.p.kin_memory_s
        if memory and (warm is None or self._kin_t is None):
            # the start: the electrons are Maxwellian at the gas
            # temperature, which LIGHT's own closure holds exactly. That
            # Maxwellian, with LIGHT's mean energy on the grid, is the
            # history the first time step starts from
            g = self.kin._grid(boltzmann.Solver.EMAX_MIN)
            self.kin.grid = g
            F = g.maxwellian_at_mean(1.5 * te)
            self._kin_sol = boltzmann.Solution(F, 0.0, 0.0, g.mean(F),
                                               False, 0, g)
            self._kin_t = t
            self._kin_ref = self._kin_observed(y)
            self.kin_log.append({"t": t, "on": on, "te": te, "skipped": True})
            return
        if memory:
            # every point is one step from the last centre over the time
            # since it, to its own mean energy (Solver.step)
            dt = max(t - self._kin_t, 1e-18)
            centre = self.kin.step(s, warm, dt, target_mean=1.5 * te)
            sols[0.0] = centre
            for r in (0.08, 0.16, -0.08, -0.16):
                self.kin.grid = centre.grid
                sols[r] = self.kin.step(s, warm, dt,
                                        target_mean=1.5 * te * math.exp(r),
                                        adapt=False)
        else:
            centre = self.kin.eedf(s, 1.5 * te, warm=warm)
            sols[0.0] = centre
            # each point warm-starts from the last table's point at the
            # same offset (else its inner neighbour), on the centre's grid,
            # and without the growth term where the centre has none
            prev = self._kin_pts or {}
            for r, inner in ((0.08, 0.0), (0.16, 0.08), (-0.08, 0.0),
                             (-0.16, -0.08)):
                self.kin.grid = centre.grid       # a stretch is the point's
                sols[r] = self.kin.eedf(s, 1.5 * te * math.exp(r),
                                        warm=prev.get(r, sols[inner]),
                                        growth=not centre.growth_off,
                                        adapt=False)
        self.kin.grid = centre.grid
        self._kin_pts = sols
        # points below the no-field mean are the no-field solution: one
        # entry for them (np.interp holds the ends)
        pts = sorted({round(math.log(v.te), 9): v
                      for v in sols.values()}.items())
        q = [self._kin_quantities(s, v) for _, v in pts]
        lnte = np.array([a for a, _ in pts])
        self._tab = {
            "lnte": lnte,
            "lnk": {j: np.log(np.maximum([qq[0][j] for qq in q], 1e-300))
                    for j in q[0][0]},
            "cap": {j: np.array([qq[1][j] for qq in q]) for j in q[0][1]},
            "pel": np.array([qq[2] for qq in q]),
            # v0_4: the EM wave's absorbed power per electron [eV/s] from
            # the solver's own field term (0 without a wave)
            "pem": np.array([qq[3] for qq in q])}
        # LIGHT's own Maxwellian values at the table's two ends, which carry
        # it beyond them (_tab_end)
        n = np.maximum(y[:self.ns], 0.0) * self.N0
        ends = (math.exp(lnte[0]), math.exp(lnte[-1]))
        byj = {c.j: c for c in self.compiled}    # rows may have been removed
        self._tab["lnkM"] = {
            j: np.log(np.maximum([self._k(byj[j], e, s.tg, n, ne)
                                  for e in ends], 1e-300))
            for j in self._tab["lnk"]}
        self._tab["capM"] = {
            j: np.array([byj[j].cap(e) for e in ends])
            for j in self._tab["cap"] if byj[j].cap is not None}
        self._tab["pelM"] = np.array([self._pel_maxwell(e, s.tg * KB_EV)
                                      for e in ends])
        self._tab["pemM"] = np.array([self._em_maxwell(e, n, ne)
                                      for e in ends])
        self._kin_sol = sol = centre
        t_prev, self._kin_t = self._kin_t, t
        self._kin_ref = self._kin_observed(y)
        ratio = balance = float("nan")
        if on:
            # the drive's electron power, from LIGHT's own equations (what
            # the pulse adds to dUe/dt). A field carries the part not
            # injected as electrons with a spectrum: all of the given
            # power, none of a beam's. The steady-shape kinetic balance has
            # no d<eps>/dt, so the two agree where field power - field part
            # + n_e d<eps>/dt = 0, with n_e d<eps>/dt = dUe/dt - <eps>
            # dn_e/dt from LIGHT's own equations here
            d = self.rhs(t, y, True)
            p_drive = (d[self.i_ue] - self.rhs(t, y, False)[self.i_ue]) \
                * self.N0
            if p_drive > 0.0:
                # the DC field carries the given power alone: not a beam's
                # electrons, not the EM wave's power (v0_4, its own term)
                p_field = p_drive - sum(r * m for r, _, m in
                                        self.kinetic_injection(t, y, on)) \
                    - self.em_power(te, n, ne)
                field_w = self.kin.field_power(s, sol) * ne
                dne = sum(SPECIES[self.species[i]].charge * d[i]
                          for i in range(self.ns)) * self.N0
                shape = d[self.i_ue] * self.N0 - 1.5 * te * dne
                if memory:
                    # a time step's balance holds its own d<eps>/dt, the
                    # change of the mean over the step: what is left is
                    # LIGHT's derivative at t against the step's average
                    shape -= ne * (sol.mean - warm.mean) / (t - t_prev)
                if p_field > 1e-6 * p_drive:
                    ratio = field_w / p_field
                balance = (field_w - p_field + shape) / p_drive
        self.kin_log.append({
            "t": t, "on": on, "te": te, "te_kin": sol.te, "en_td": sol.en_td,
            "reached": sol.reached, "power_ratio": ratio, "balance": balance,
            "emax": sol.grid.emax, "growth_off": sol.growth_off,
            "iters": sol.iters, "nu": sol.nu, "points": len(pts),
            "pem": float(self._tab["pem"][np.searchsorted(lnte, math.log(
                sol.te)) if len(lnte) > 1 else 0]) if s.omega > 0.0 else 0.0,
            "k": {self.rows[j].rid: v for j, v in
                  self._kin_quantities(s, sol)[0].items()}})

    def _drift_events(self):
        ref = self._kin_ref
        step = self.p.kin_step
        te_step = self.KIN_TABLE_REACH if self._tab is not None else step

        def make(i, lim):
            floor = self.KIN_FLOORS[i]
            if ref[i] + floor <= 0.0:
                return None              # an input that is 0 and stays so

            def ev(t, y):
                v = self._kin_observed(y)[i]
                return abs(math.log((v + floor) / (ref[i] + floor))) - lim
            ev.terminal = True
            ev.direction = 1
            return ev
        evs = [make(0, te_step)] + [make(i, step)
                                    for i in range(1, len(ref))]
        return [e for e in evs if e is not None]

    # ------------------------------------------------------------ compile
    def _compile(self, j, r):
        kind = r.rate[0]
        if kind == "const":
            rate = float(r.rate[1])
        elif kind == "inv_tau":
            frac = float(r.rate[2]) if len(r.rate) > 2 else 1.0
            rate = frac / float(getattr(self.p, r.rate[1]))
        elif kind == "param":
            frac = float(r.rate[2]) if len(r.rate) > 2 else 1.0
            rate = frac * float(self.rates.get(r.rate[1],
                                               BANDS[r.rate[1]].default))
        elif kind == "xsec" and self.eedf == "two_t" and xsec.in_tail(
                *r.rate[1:4]):
            kind = "xsec_tail"
            rate = xsec.two_t_table(r.rate[1], r.rate[2], r.rate[3])
        elif kind == "xsec":
            rate = xsec.rate_table(r.rate[1], r.rate[2], r.rate[3])
        elif kind == "superelastic":
            rate = xsec.superelastic_table(r.rate[1], r.rate[2], r.rate[3])
        elif kind in ("stepwise", "dr_sqrt", "flannery", "elastic",
                      "te_pow", "ionion"):
            rate = r.rate[1:]
        elif kind == "balance":
            k_fwd, k_ref, upper, lower = r.rate[1:]
            rate = (float(k_fwd) / float(k_ref),
                    (e_form(upper) - e_form(lower)) / KB_EV)
        elif kind == "te_fit":
            # A Te^B exp(-C/Te); A a number or a band; a reverse row
            # carries A g_ratio and C - dE (detailed balance over the
            # ledger's gap)
            a_spec, b_pow, c_ev = r.rate[1:4]
            a = (self.rates.get(a_spec, BANDS[a_spec].default)
                 if isinstance(a_spec, str) else float(a_spec))
            if len(r.rate) > 4:
                tag, g_ratio, upper, lower = r.rate[4:8]
                if tag != "reverse":
                    raise ValueError(f"{r.rid}: te_fit tag {tag!r}")
                a *= float(g_ratio)
                c_ev = float(c_ev) - (e_form(upper) - e_form(lower))
            rate = (float(a), float(b_pow), float(c_ev))
        else:
            raise ValueError(f"{r.rid}: rate kind {kind!r}")
        factor = 1.0
        for name, how in r.frac:
            v = self.rates.get(name, BANDS[name].default)
            factor *= v if how == "x" else 1.0 - v
        factor *= self._row_scale.get(r.rid, 1.0)
        delta = {}
        for s in r.reactants:
            if s != ELECTRON:
                delta[s] = delta.get(s, 0) - 1
        for s in r.products:
            if s != ELECTRON:
                delta[s] = delta.get(s, 0) + 1
        stoich = tuple((self.idx[s], c) for s, c in delta.items() if c != 0)
        route = energy_routing(r)
        return Compiled(r, j, tuple(self.idx[s] for s in r.heavy),
                        r.n_electrons_in, stoich, route, kind, rate,
                        self._capture_energy(r, kind) if route["capture"]
                        else None, factor)

    def _capture_energy(self, r, kind):
        """The energy one captured electron takes from the electron gas.
        For a Maxwellian it is <eps sigma v>/<sigma v> = (3/2 + d ln k /
        d ln Te) Te, so a Te-independent rate books (3/2)Te, R11p's
        Te^-1/2 template books Te, and an LXCat process reads its own
        table. None keeps (3/2)Te (capture='three_halves', v0's)."""
        if self.capture == "three_halves":
            return None
        if kind == "xsec":
            return xsec.capture_energy_table(r.rate[1], r.rate[2], r.rate[3])
        if kind == "dr_sqrt":
            return lambda te: te if te > 0.02 else 1.5 * te
        if kind == "te_pow":
            k0, te0, n_pow = r.rate[1:]
            return lambda te: ((1.5 + n_pow) * te if te > TE_FLOOR_POW
                               else 1.5 * te)
        if kind in ("const", "inv_tau", "param", "balance", "ionion"):
            return None
        raise ValueError(f"{r.rid}: no captured-energy rule for rate kind "
                         f"{kind!r}")

    # ------------------------------------------------------------ helpers
    def electron_density(self, n):
        if self.ne_mode == "clamped":
            return self.p.ne_cm3
        ne = n[self.cations].sum()
        if self.anions:
            ne -= n[self.anions].sum()
        return max(ne, NE_FLOOR)

    def _k(self, c, te, tg, n, ne):
        kind = c.kind
        if kind in ("const", "inv_tau", "param"):
            k = c.rate
        elif kind in ("xsec", "superelastic"):
            k = float(c.rate(te))
        elif kind == "xsec_tail":
            n_gas = max(float(n.sum()), 1.0)
            r = self.tail(ne / n_gas, n[self.idx["F2"]] / n_gas)
            k = c.rate(te, r * te)
        elif kind == "stepwise":
            t = max(te, TE_FLOOR_STEPWISE)
            k = 1.56e-7 * t ** 0.71 * math.exp(-c.rate[0] / t)
        elif kind == "dr_sqrt":
            k = c.rate[0] / math.sqrt(max(te, 0.02))
        elif kind == "te_pow":
            k = c.rate[0] * (max(te, TE_FLOOR_POW) / c.rate[1]) ** c.rate[2]
        elif kind == "flannery":
            # harmonic mean of a three-body part a N/1.7 and a mobility part
            # a/N; an optional (Tg/300)^n3 acts on the three-body part
            n_am = n.sum() / AMAGAT_CM3
            f3 = ((max(tg, 1.0) / 300.0) ** c.rate[2] if len(c.rate) > 2
                  else 1.0)
            k = c.rate[0] * n_am * f3 / (n_am * n_am * f3 + 1.7) * c.rate[1]
            k *= self._screening(ne, n, tg)
        elif kind == "balance":
            # k_fwd / K(Tg), K(Tg) = K_ref exp(theta (1/Tg - 1/300))
            k = c.rate[0] * math.exp(-c.rate[1] * (1.0 / max(tg, 1.0)
                                                   - 1.0 / 300.0))
        elif kind == "te_fit":
            t = max(te, TE_FLOOR_FIT)
            k = c.rate[0] * t ** c.rate[1] * math.exp(-c.rate[2] / t)
        elif kind == "ionion":
            # Flannery's turnover as the harmonic mean of the three-body
            # limb k3 N (Tg/300)^n3 and the Langevin limb kL / N (N in
            # amagat, both cm3/s at 1 amagat), times the Debye screening
            k3, kl, n3 = c.rate
            n_am = n.sum() / AMAGAT_CM3
            a3 = k3 * n_am * (max(tg, 1.0) / 300.0) ** n3
            al = kl / n_am
            k = a3 * al / (a3 + al) * self._screening(ne, n, tg)
        else:
            raise ValueError(kind)
        tg_pow = c.row.tg_pow
        if self.p.trimer_tg_pow is not None and c.row.rid in ("G5", "G5C"):
            tg_pow = self.p.trimer_tg_pow
        if tg_pow:
            k *= (max(tg, 1.0) / 300.0) ** tg_pow
        return k * c.factor

    def _rate(self, c, te, tg, n, ne):
        """A row's rate constant: the kinetic table's (read in ln Te) where
        the table carries the row, else the closure's own (_k)."""
        tab = self._tab
        if tab is None or c.j not in tab["lnk"]:
            return self._k(c, te, tg, n, ne)
        lt = math.log(te)
        end = self._tab_end(lt)
        if end is None:
            return math.exp(float(np.interp(lt, tab["lnte"],
                                            tab["lnk"][c.j])))
        km = max(self._k(c, te, tg, n, ne), 1e-300)
        return math.exp(tab["lnk"][c.j][end] + math.log(km)
                        - tab["lnkM"][c.j][end])

    # ------------------------------------------------------------ v0_4
    def wall_te(self, ue, n):
        """The electron energy per cation that drives the ambipolar flow
        [eV]: U_e / (1.5 sum n+), which is Te where the plasma is
        electropositive and falls with the electron gas where the anions
        hold the charge (an ion-ion plasma has no ambipolar field to speak
        of). Smooth in the densities, capped at WALL_TE_MAX."""
        n_cat = max(sum(n[i] for i in self.cations), NE_FLOOR)
        return min(ue / (1.5 * n_cat), WALL_TE_MAX)

    def wall_loss(self, ue, tg, n):
        """The charges' drift and diffusion to the electrodes of a gap
        (Params.gap_cm): per charged species the loss frequency [1/s] and
        the particle flux [cm^-3 s^-1]. Cations leave at the ambipolar
        rate D_a / Lambda^2, D_a = mu_i k (Te_w + Tg) / e (the electrons'
        pressure drives them, the standard electropositive closure, with
        Te_w = wall_te), with mu_i = K0 N0 / N from the measured reduced
        mobilities; anions at their free diffusion mu_i k Tg / e /
        Lambda^2, since the ambipolar field holds them in. Lambda = gap /
        pi (the lowest slab mode). The electrons leave with the cations,
        which quasineutrality books by itself. Returns ((species index,
        charge, flux), ...)."""
        if not self.wall:
            return ()
        lam2 = (self.p.gap_cm / math.pi) ** 2
        n_am = n.sum() / AMAGAT_CM3
        tg_ev = tg * KB_EV
        te_w = self.wall_te(ue, n)
        out = []
        for i, q, k0 in self.wall:
            mu = k0 / max(n_am, 1e-12)                  # cm2/V/s
            kt = (te_w + tg_ev) if q > 0 else tg_ev     # eV
            nu = mu * kt / lam2                         # 1/s
            out.append((i, q, nu * n[i]))
        return out

    def _em_maxwell(self, te, n, ne):
        """ib.py's absorbed EM power per electron [eV/s] over a Maxwellian
        at te: the Lorentz-gas kernel on the electron-neutral collisions
        plus Turnbull et al.'s electron-ion absorption (times ib_scale),
        with the Langdon factor on the electron-ion part for a Maxwellian
        closure only (a solved EEDF must not carry it: double counting)."""
        n_kr, n_f2 = n[self.idx["Kr"]], n[self.idx["F2"]]
        n_gas = n_kr + n_f2
        if self._em is None or n_gas <= 0.0 or ne <= 0.0:
            return 0.0
        ni = sum(n[i] for i in self.cations)
        w = 2.0 * math.pi * self.p.em_f_Hz
        args = (te, w, self.p.em_E_V_cm * 100.0, n_gas * 1e6,
                {"Kr": n_kr / n_gas, "F2": n_f2 / n_gas}, ne * 1e6)
        langdon = self.eedf != "kinetic"
        total = self._em.em_power_per_electron(*args, ni * 1e6, langdon=langdon)
        if self.ib_scale == 1.0 or ni <= 0.0:
            return total
        # the band acts on the electron-ion share of the absorbed power,
        # the MARGINAL one (the power with the ions less the power
        # without them), as the kinetic closure's does. The kernel's
        # printed "ei" part is the nu_ei share of the numerator; where
        # nu_en is not << omega (1 atm Kr at 95 GHz: nu_en = 4.4 omega)
        # the two differ, 1.5x there, and only the marginal one is what
        # the ions add
        atoms = self._em.em_power_per_electron(*args, 0.0, langdon=langdon)
        return atoms + self.ib_scale * (total - atoms)

    def em_power(self, te, n, ne):
        """The EM wave's absorbed power [eV cm^-3 s^-1] at this state: the
        kinetic table's value where the EEDF is solved with the wave's
        field term (boltzmann.Solver.em_power, read in ln Te and carried
        past the table's ends by the Maxwellian kernel's own change, as
        the elastic power is), else the Maxwellian kernel (_em_maxwell)."""
        if self._em is None or ne <= 0.0:
            return 0.0
        tab = self._tab
        if tab is not None and "pem" in tab and tab["pem"].any():
            lt = math.log(te)
            end = self._tab_end(lt)
            if end is None:
                per_e = float(np.interp(lt, tab["lnte"], tab["pem"]))
            else:
                per_e = tab["pem"][end] + self._em_maxwell(te, n, ne) \
                    - tab["pemM"][end]
            return max(per_e, 0.0) * ne
        return self._em_maxwell(te, n, ne) * ne

    def _tab_end(self, lt):
        """None inside the kinetic table, else the end LIGHT's ln Te lies
        beyond (0 or -1). Beyond an end the table's values there are
        carried with LIGHT's own Maxwellian Te dependence. Inside it LIGHT's
        mean energy is the kinetic one. Below the bottom, the field is at
        its floor and LIGHT's electrons are colder than any steady solution
        (the no-field one is hotter): a flat table there would leave
        LIGHT's energy equation no loss that falls with its Te, and its Te
        would creep toward the kinetic mean at the electrons' turnover rate
        alone (1.3 against 1.6 eV at Tahoma, 2026-09-30)."""
        lnte = self._tab["lnte"]
        if lt < lnte[0]:
            return 0
        if lt > lnte[-1]:
            return -1
        return None

    def _pel_maxwell(self, te, tg_ev):
        """Kr's elastic power per electron per Kr atom over a Maxwellian
        [eV cm^3/s], LIGHT's own form."""
        return float(self.kel(te)) * (1.0 - tg_ev / te)

    def _screening(self, ne, n, tg):
        """Debye screened-capture factor for ion-ion recombination (REBUILD
        chemistry.ion_ion_screening, cgs in / SI inside)."""
        n_an = n[self.anions].sum() if self.anions else 0.0
        n_scr = max(2.0 * (ne + n_an), 1.0) * 1e6
        tk = max(tg, 1.0)
        lam_d = math.sqrt(EPS0 * KB * tk / (n_scr * QE ** 2))
        r_cap = 2.0 * QE ** 2 / (3.0 * 4.0 * math.pi * EPS0 * KB * tk)
        return 1.0 / (1.0 + r_cap / lam_d)

    def heating_power(self, te, n, ne):
        """Power into the electron gas [eV cm^-3 s^-1] during the pulse:
        the given input power, uniform and direct (the EM wave's absorbed
        power is added in rhs and booked apart, E_em)."""
        return self.p.power_W_cm3 / QE

    # ------------------------------------------------------------ RHS
    def rhs(self, t, y, on):
        N0 = self.N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        ue = y[self.i_ue] * N0
        tg = y[self.i_tg]
        tg_ev = tg * KB_EV
        ne = self.electron_density(n)
        te = max(ue / (1.5 * ne), 1e-3)

        dy = np.zeros(self.size)
        dn = dy[:self.ns]
        heat = 0.0
        p_in = self.heating_power(te, n, ne) if on else 0.0
        p_em = self.em_power(te, n, ne) if on else 0.0
        tab = self._tab
        if tab is not None:            # the kinetic EEDF's table
            lt = math.log(te)
            end = self._tab_end(lt)
            if end is None:
                pel = float(np.interp(lt, tab["lnte"], tab["pel"]))
            else:
                # carried by the Maxwellian's own change, which passes
                # through 0 at Tg (a ratio would not)
                pel = tab["pel"][end] + self._pel_maxwell(te, tg_ev) \
                    - tab["pelM"][end]
            p_el = ne * n[self.idx["Kr"]] * pel
        else:
            p_el = (ne * n[self.idx["Kr"]] * float(self.kel(te))
                    * (1.0 - tg_ev / te))
        due = p_in + p_em - p_el
        heat += p_el
        cap_heat = 0.0
        # v0_4: the charges' loss to the walls (Params.gap_cm). Each lost
        # cation takes its formation energy to the wall, where it
        # recombines, and the electron that leaves with it takes 2 Te from
        # the electron gas to the wall and about Te more through the
        # ambipolar field to the ion, which heats the gas. A lost anion
        # takes its own formation energy to the wall.
        e_wall = e_ion_heat = e_wall_th = 0.0
        dn_wall = np.zeros(self.ns)
        wl = self.wall_loss(ue, tg, n)
        te_w = self.wall_te(ue, n) if wl else 0.0
        for i, q, flux in wl:
            dn[i] -= flux
            dn_wall[i] -= flux
            e_wall += flux * self.eform[i]
            e_wall_th += flux * self.cv[i] * tg_ev      # leaves at Tg
            if q > 0:
                e_wall += flux * WALL_E_ELECTRON_TE * te_w
                e_ion_heat += flux * WALL_E_AMBIPOLAR_TE * te_w
                due -= flux * (WALL_E_ELECTRON_TE + WALL_E_AMBIPOLAR_TE) * te_w
        heat += e_ion_heat
        for c in self.compiled:
            if c.kind == "elastic":
                continue
            flux = self._rate(c, te, tg, n, ne)
            for i in c.heavy:
                flux *= n[i]
            if c.ne_pow:
                flux *= ne ** c.ne_pow
            for i, coeff in c.stoich:
                dn[i] += coeff * flux
            rt = c.route
            due += (rt["to_e"] - rt["from_e"]) * flux
            heat += rt["heat"] * flux
            if rt["capture"]:
                if tab is not None and c.j in tab["cap"]:
                    if end is None or c.j not in tab["capM"]:
                        e_cap = float(np.interp(lt, tab["lnte"],
                                                tab["cap"][c.j]))
                    else:
                        e_cap = tab["cap"][c.j][end] * c.cap(te) \
                            / tab["capM"][c.j][end]
                else:
                    e_cap = 1.5 * te if c.cap is None else c.cap(te)
                cap = e_cap * rt["capture"] * flux
                due -= cap
                cap_heat += cap
            dy[self.i_ev + c.j] = flux
        heat += cap_heat
        # constant-volume gas energy U_th = sum n_s cv_s k Tg, dU_th/dt = heat:
        # particles a reaction creates (F2 -> 2F) take their cv k Tg from
        # it; particles lost to the walls leave at Tg with theirs (v0_4,
        # E_wall_th) and change Tg not at all
        c_v = float((n * self.cv).sum()) * KB_EV
        dy[self.i_tg] = (heat - tg_ev * float(((dn - dn_wall)
                                                * self.cv).sum())) / c_v
        dn /= N0
        dy[self.i_ue] = due / N0
        dy[self.i_led + 0] = p_in / N0
        dy[self.i_led + 1] = p_el / N0
        dy[self.i_led + 2] = cap_heat / N0
        dy[self.i_led + 3] = p_em / N0
        dy[self.i_led + 4] = e_wall / N0
        dy[self.i_led + 5] = e_ion_heat / N0
        dy[self.i_led + 6] = e_wall_th / N0
        dy[self.i_ev:] /= N0
        return dy

    # ------------------------------------------------------------ run
    def initial_state(self):
        p = self.p
        y = np.zeros(self.size)
        N0 = self.N0
        y[self.idx["Kr"]] = (1.0 - p.x_f2) - p.ne_cm3 / N0
        y[self.idx["F2"]] = p.x_f2
        y[self.idx["Kr+"]] = p.ne_cm3 / N0
        for s, val in p.init_extra_cm3:
            y[self.idx[s]] += val / N0
        y[self.i_ue] = 1.5 * p.ne_cm3 * p.tg0_K * KB_EV / N0
        y[self.i_tg] = p.tg0_K
        return y

    def run(self) -> "Result":
        if self.eedf == "kinetic":
            return self._run_kinetic()
        p = self.p
        y0 = self.initial_state()
        t_grid = np.linspace(0.0, p.t_end_s, p.n_save)
        # integrate on-pulse and afterglow separately: the drive switches
        # off discontinuously and a stiff solver must not step across it
        segs = [(0.0, min(p.pulse_s, p.t_end_s), True),
                (p.pulse_s, p.t_end_s, False)]
        segs = [s for s in segs if s[1] > s[0]]
        atol = np.full(self.size, 1e-18)
        atol[self.i_tg] = 1e-8
        ts, ys, nfev = [], [], 0
        for a, b, on in segs:
            te_ = t_grid[(t_grid > a) & (t_grid < b)]
            t_eval = np.concatenate([[a], te_, [b]])
            sol = solve_ivp(lambda t, y: self.rhs(t, y, on), (a, b), y0,
                            method=p.method, t_eval=t_eval, rtol=p.rtol,
                            atol=atol, first_step=1e-15)
            if not sol.success:
                raise RuntimeError(f"solver failed on [{a:.3e}, {b:.3e}] s: "
                                   f"{sol.message}")
            nfev += sol.nfev
            keep = slice(None) if not ts else slice(1, None)   # no duplicate
            ts.append(sol.t[keep])
            ys.append(sol.y[:, keep])
            y0 = sol.y[:, -1]
        return Result(self, np.concatenate(ts), np.concatenate(ys, axis=1),
                      nfev)

    def _run_kinetic(self) -> "Result":
        """The run with the kinetic EEDF. The stiff solver integrates with
        the table held, stops when Te has moved KIN_TABLE_REACH in ln Te or
        n_e or the F2 fraction kin_step in its ln since the last EEDF (or
        kin_max_dt_s has passed), and the EEDF is solved again there. The
        output grid is the same as the other closures' (the stops are not
        saved)."""
        p = self.p
        # rows edited after construction (ML removes or rescales some) are
        # what the kinetic equation carries, and a second run starts cold
        self._kinetic_setup()
        self._kin_sol = self._kin_pts = self._kin_t = None
        self.kin_log = []
        y0 = self.initial_state()
        t_grid = np.linspace(0.0, p.t_end_s, p.n_save)
        segs = [(0.0, min(p.pulse_s, p.t_end_s), True),
                (p.pulse_s, p.t_end_s, False)]
        segs = [s for s in segs if s[1] > s[0]]
        atol = np.full(self.size, 1e-18)
        atol[self.i_tg] = 1e-8
        ts, ys, nfev = [], [], 0
        for a, b, on in segs:
            t, y = a, y0
            self._kinetic_update(t, y, on)
            seg_t, seg_y = [], []
            first = True
            stall = 0
            while t < b:
                cap = p.kin_max_dt_s if on else max(
                    p.kin_max_dt_s, p.kin_afterglow_frac * (t - a))
                t_stop = min(b, t + cap)
                inner = t_grid[(t_grid > t) & (t_grid < t_stop)]
                t_eval = np.concatenate([[t], inner, [t_stop]])
                try:
                    sol = solve_ivp(lambda tt, yy: self.rhs(tt, yy, on),
                                    (t, t_stop), y, method=p.method,
                                    t_eval=t_eval, rtol=p.rtol, atol=atol,
                                    first_step=1e-15 if first else None,
                                    events=self._drift_events())
                except ValueError as exc:
                    if "different signs" not in str(exc):
                        raise
                    # scipy's event root-finder: the step's ends crossed a
                    # trigger that its dense output cannot bracket (an
                    # input jumping, n_e at its floor in an attaching
                    # afterglow; v0_4, 2026-10-09). A short step without
                    # the triggers, and the table is renewed after it
                    t_go = min(b, t + max(1e-3 * cap, 1e-14))
                    sol = solve_ivp(lambda tt, yy: self.rhs(tt, yy, on),
                                    (t, t_go), y, method=p.method,
                                    rtol=p.rtol, atol=atol,
                                    first_step=min(1e-15, t_go - t))
                    if not sol.success:
                        raise RuntimeError(f"solver failed on [{t:.3e}, "
                                           f"{t_go:.3e}] s: {sol.message}")
                    nfev += sol.nfev
                    self.kin_log.append({"t": t_go, "on": on, "te": None,
                                         "skipped": True,
                                         "event_bracket": True})
                    first = False
                    t, y = t_go, sol.y[:, -1]
                    if t < b:
                        self._kinetic_update(t, y, on)
                    continue
                if not sol.success:
                    raise RuntimeError(f"solver failed on [{t:.3e}, "
                                       f"{t_stop:.3e}] s: {sol.message}")
                nfev += sol.nfev
                for k, tt in enumerate(sol.t):
                    if (first and k == 0) or (k > 0 and (tt < t_stop or
                                                        tt == b)):
                        seg_t.append(tt)
                        seg_y.append(sol.y[:, k])
                if sol.status == 1:
                    hit = [(ev[0], ye[0]) for ev, ye in
                           zip(sol.t_events, sol.y_events) if len(ev)]
                    t_end, y_end = min(hit, key=lambda h: h[0])
                    # an event that lands on b ends the segment
                    t_end = min(t_end, b)
                else:
                    t_end, y_end = t_stop, sol.y[:, -1]
                if t_end <= t:
                    # an input passed its trigger inside the solver's first
                    # step, so the event lands where the segment began (a
                    # beam's secondaries heat a seed of electrons within a
                    # femtosecond): a short step without the triggers,
                    # doubled while it repeats, keeps the run moving
                    stall += 1
                    t_go = min(b, t + 1e-15 * 2.0 ** stall)
                    sol = solve_ivp(lambda tt, yy: self.rhs(tt, yy, on),
                                    (t, t_go), y, method=p.method,
                                    rtol=p.rtol, atol=atol,
                                    first_step=min(1e-15, t_go - t))
                    if not sol.success:
                        raise RuntimeError(f"solver failed on [{t:.3e}, "
                                           f"{t_go:.3e}] s: {sol.message}")
                    nfev += sol.nfev
                    t_end, y_end = t_go, sol.y[:, -1]
                else:
                    stall = 0
                first = False
                t, y = t_end, y_end
                if t < b:
                    self._kinetic_update(t, y, on)
            if not seg_t or seg_t[-1] != b:
                seg_t.append(b)
                seg_y.append(y)
            keep = slice(None) if not ts else slice(1, None)
            ts.append(np.array(seg_t)[keep])
            ys.append(np.stack(seg_y, axis=1)[:, keep])
            y0 = y
        res = Result(self, np.concatenate(ts), np.concatenate(ys, axis=1),
                     nfev)
        res.kinetic = list(self.kin_log)
        return res




class Result:
    """Trajectory + ledger in physical units (cm^-3, eV cm^-3, K, eV)."""

    def __init__(self, model: Model, t, y, nfev):
        self.model = model
        self.p = model.p
        self.t = t
        N0 = model.N0
        self.n = {s: y[i] * N0 for s, i in model.idx.items()}
        dens = y[:model.ns] * N0
        self.ne = np.array([model.electron_density(np.maximum(dens[:, k], 0))
                            for k in range(len(t))])
        self.ue = y[model.i_ue] * N0
        self.te = self.ue / (1.5 * self.ne)
        self.tg = y[model.i_tg]
        self.ledger = {name: y[model.i_led + i] * N0
                       for i, name in enumerate(LEDGER)}
        self.events = {r.rid: y[model.i_ev + j] * N0
                       for j, r in enumerate(model.rows)}
        # chemical + electron energy [eV/cm^3]
        self.stored = (model.eform[:, None] * dens).sum(axis=0) + self.ue
        # gas thermal energy [eV/cm^3]; its growth must equal the booked heat
        self.u_thermal = (model.cv[:, None] * dens).sum(axis=0) * KB_EV \
            * self.tg
        n_heavy = dens.sum(axis=0)
        self.pressure_atm = (n_heavy + self.ne) * 1e6 * KB * self.tg / ATM_PA
        self.nfev = nfev

    def at_end(self, key):
        return self.n[key][-1]


def run(params: Params = None, **kw) -> Result:
    """Convenience: run(pressure_atm=2, ne_cm3=1e16) or run(Params(...))."""
    p = params if params is not None else Params()
    if kw:
        p = p.replace(**kw)
    return Model(p).run()
