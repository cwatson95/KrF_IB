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

from . import xsec
from .network import BANDS, build, energy_routing
from .species import ELECTRON, SPECIES, e_form

KB = 1.380649e-23             # J/K
QE = 1.602176634e-19          # J/eV
KB_EV = KB / QE               # eV/K
ATM_PA = 101325.0
AMAGAT_CM3 = 2.6867811e19     # Loschmidt, 273.15 K and 1 atm
EPS0 = 8.8541878128e-12

BASE_SPECIES = ("Kr", "F2", "F", "Kr*", "Kr+", "KrF(B)", "Kr2F*")
LEDGER = ("E_dep", "heat_elastic", "heat_capture")
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
    # (run_v0_2.py --bands), never to match a device.
    rates: tuple = ()
    network: str = "basic"             # 'basic' (KrF(B,C) set) or 'v0'
    groups: tuple = ()                 # optional groups ('all' = every one)
    # the electron closure (2026-09-28, step 3). 'auto' gives the basic
    # network the two-temperature tail and the exact captured energy, and
    # keeps v0 on the notes' Maxwellian with (3/2)Te per captured electron
    eedf: str = "auto"                 # 'two_t' | 'maxwellian'
    capture: str = "auto"              # 'exact' | 'three_halves'
    init_extra_cm3: tuple = ()         # ((species, density), ...) added at t=0
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


class Model:
    def __init__(self, p: Params):
        basic = p.network == "basic"
        self.ne_mode = (("free" if basic else "clamped")
                        if p.ne_mode == "auto" else p.ne_mode)
        if self.ne_mode not in ("clamped", "free"):
            raise ValueError(f"ne_mode {p.ne_mode!r}")
        self.eedf = (("two_t" if basic else "maxwellian") if p.eedf == "auto"
                     else p.eedf)
        self.capture = (("exact" if basic else "three_halves")
                        if p.capture == "auto" else p.capture)
        if self.eedf not in ("two_t", "maxwellian"):
            raise ValueError(f"eedf {p.eedf!r}")
        if self.capture not in ("exact", "three_halves"):
            raise ValueError(f"capture {p.capture!r}")
        self.tail = (xsec.tail_ratio_table() if self.eedf == "two_t"
                     else None)
        self.p = p
        self.rows = build(p.groups, p.network)
        self.rates = dict(p.rates)
        read = {r.rate[1] for r in self.rows if r.rate[0] == "param"}
        if set(self.rates) - read:
            raise ValueError(f"rates {sorted(set(self.rates) - read)}: no "
                             f"row of this network reads them (bands here: "
                             f"{sorted(read)}; the trimer law is "
                             f"Params.trimer_tg_pow)")
        if any(not v >= 0.0 for v in self.rates.values()):
            raise ValueError(f"rates must be >= 0: {p.rates}")
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
        self.compiled = [self._compile(j, r) for j, r in enumerate(self.rows)]

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
                      "te_pow"):
            rate = r.rate[1:]
        elif kind == "balance":
            k_fwd, k_ref, upper, lower = r.rate[1:]
            rate = (float(k_fwd) / float(k_ref),
                    (e_form(upper) - e_form(lower)) / KB_EV)
        else:
            raise ValueError(f"{r.rid}: rate kind {kind!r}")
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
                        else None)

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
        if kind in ("const", "inv_tau", "param", "balance"):
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
        else:
            raise ValueError(kind)
        tg_pow = c.row.tg_pow
        if self.p.trimer_tg_pow is not None and c.row.rid in ("G5", "G5C"):
            tg_pow = self.p.trimer_tg_pow
        if tg_pow:
            k *= (max(tg, 1.0) / 300.0) ** tg_pow
        return k

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
        the given input power, uniform and direct."""
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
        p_el = (ne * n[self.idx["Kr"]] * float(self.kel(te))
                * (1.0 - tg_ev / te))
        due = p_in - p_el
        heat += p_el
        cap_heat = 0.0
        for c in self.compiled:
            if c.kind == "elastic":
                continue
            flux = self._k(c, te, tg, n, ne)
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
                e_cap = 1.5 * te if c.cap is None else c.cap(te)
                cap = e_cap * rt["capture"] * flux
                due -= cap
                cap_heat += cap
            dy[self.i_ev + c.j] = flux
        heat += cap_heat
        # constant-volume gas energy U_th = sum n_s cv_s k Tg, dU_th/dt = heat:
        # particles a reaction creates (F2 -> 2F) take their cv k Tg from it
        c_v = float((n * self.cv).sum()) * KB_EV
        dy[self.i_tg] = (heat - tg_ev * float((dn * self.cv).sum())) / c_v
        dn /= N0
        dy[self.i_ue] = due / N0
        dy[self.i_led + 0] = p_in / N0
        dy[self.i_led + 1] = p_el / N0
        dy[self.i_led + 2] = cap_heat / N0
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
