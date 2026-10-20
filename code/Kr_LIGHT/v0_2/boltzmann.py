"""The kinetic EEDF of v0_3: the two-term Boltzmann equation with
electron-electron collisions, solved during the run.

The scheme is ML's (ML/ml/boltzmann: the website lxcat build's two-term
solver, Hagelaar & Pitchford, PSST 14, 722 (2005), the form LoKI-B solves,
and ML's e-e term, ee.py), written again here on sparse matrices, with the
gas terms assembled once per grid. For the isotropic part F0 of the EEDF,
normalised as int sqrt(eps) F0 d eps = 1, on cells of width du with centres
e_i = (i + 1/2) du, and per gas particle N:

    d/d eps [ W F0 - D dF0/d eps ] = (inelastic, attachment, recombination,
                                      growth and source terms)
    W = -gamma eps^2 sigma_eps + W_ee
    D = gamma/3 (E/N)^2 eps / sigma~ + gamma (k Tg/e) eps^2 sigma_eps + D_ee
    sigma_eps = sum_g x_g (2m/M_g) sigma_el,g   (the gases whose elastic
                                                 loss the network books)
    sigma~    = sum_g x_g [sigma_el,g + sum_c sigma_c] + (nu/N)/(gamma sqrt eps)

The flux through each inner boundary is Scharfetter-Gummel, a collision
removes gamma x e_j sigma(e_j) F_j du from cell j and puts the scattered
electron (both electrons of an ionization, sharing equally) at e_j - u,
split between the two cells whose centres bracket it, and attachment only
removes (boltz_twin.py, whose conventions this follows: cross sections
linear in energy, zero outside their table and at or below a threshold).
LIGHT_Kr v0_2 adds the associative ionization of an excited level
(kind 'assoc', the Hornbeck-Molnar rows): the electron is scattered at
e_j - u as in the excitation, and the second electron is born at rest, in
the first cell.
The first cell's balance is replaced by the normalisation (boltz_twin's
rule). The e-e term is ML's, W_ee = -3 alpha A1 and D_ee = 2 alpha (A2 +
eps^1.5 A3), with the A integrals taken on ln F linear between cell
centres (Grid.ee_integrals), so that a Maxwellian of any temperature is
left alone on the grid, as it is by the exact operator. ML's quadrature
(F constant in each cell) misses that at the lowest boundaries by more
than the elastic coupling to the gas, and an undriven EEDF then condenses
into the first cell.

What v0_3 adds, so that the EEDF belongs to its own plasma and not a
swarm's (a swarm has no sources, and where attachment outruns energy
relaxation its two-term solution fails):
  - the network's electron sources and sinks. Penning ionization (G2,
    PEN_KR2) injects its electrons at the energy the network gives them
    (5.834 and 5.554 eV). Recombination (R11p) removes electrons with
    sigma v ~ eps^n, the template of its Te^n law. Attachment and
    ionization come from their cross sections. The growth term carries the
    net, nu = ionization + Penning - attachment - recombination per
    electron, so every electron is accounted for.
  - LIGHT's mean energy. Every solution has LIGHT's (3/2) Te as its mean
    energy, held by a field E/N >= 0: the input power in the pulse, the
    energy LIGHT's electrons still carry after it (the local mean energy
    approximation of fluid models, with LIGHT's energy equation choosing
    the mean). The field is found inside every step of the iteration
    (_mean_match), since at a fixed field the problem can have two
    solutions. Where even the no-field solution is hotter than the target,
    the field is 0.
  - the growth term where it holds. It keeps the EEDF's shape steady while
    n_e changes at nu. Where the growth solution cannot be as cold as the
    target (attachment of the cold electrons heats the rest), where the
    energy it moves outweighs the heating (growth_valid), or where it has
    no positive solution, the EEDF is the one without it.
  - only the processes the network books enter the collision and energy
    exchange terms, so the kinetic energy balance and LIGHT's ledger close
    the same budget. Every process of a gas's file stays in its momentum
    transfer sigma~ (in LIGHT's Kr/F2 that kept the F2 file's unbooked
    processes there and nowhere else). LIGHT_Kr's one gas is krypton
    (GASES).

LIGHT reads a table of the kinetic rates, captured energies and elastic
power at five mean energies around its own (model.py), so the kinetic mean
and LIGHT's are the same wherever the table is read.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from . import xsec

QE = xsec.QE
ME = xsec.ME
EPS0 = xsec.EPS0
AMU = xsec.AMU
KB = 1.380649e-23
GAMMA = math.sqrt(2.0 * QE / ME)            # v = GAMMA sqrt(eps) [m/s]
TD = 1e-21                                  # V m^2 per Td
C_EE = QE ** 2 / (24.0 * math.pi * EPS0 ** 2)
# the gases the grid assembles, every vendored file (LIGHT_Kr: Kr;
# LIGHT's Kr/F2 model had ("Kr", "F2"))
GASES = tuple(xsec.FILES)
M_KG = {g: xsec.MASS_AMU[g] * AMU for g in GASES}


# ---------------------------------------------------------------------------
# cross sections, as LoKI-MC (and boltz_twin) evaluate them [m^2]
# ---------------------------------------------------------------------------
def sigma(p: xsec.Process, e) -> np.ndarray:
    e = np.asarray(e, dtype=float)
    s = np.interp(e, p.energy, p.sigma * 1e-4)
    s = np.where((e >= p.energy[0]) & (e <= p.energy[-1]), s, 0.0)
    if p.kind not in ("elastic", "effective"):
        s = np.where(e > p.threshold, s, 0.0)
    return s


def elastic_mt(gas: str, e) -> np.ndarray:
    """The elastic momentum-transfer cross section [m^2]. Kr's file gives
    it. A file with an effective (total) block instead (LIGHT's F2) gives
    the elastic part as the effective less every other process on the
    effective table's own energies, clipped at zero
    (EedfGas::elasticFromEffectiveCrossSection, as gases_twin.py does)."""
    procs = xsec.processes(gas)
    el = [p for p in procs if p.kind == "elastic"]
    if el:
        return sigma(el[0], e)
    eff = next(p for p in procs if p.kind == "effective")
    x = eff.energy
    s = eff.sigma * 1e-4
    for p in procs:
        if p.kind not in ("effective", "elastic"):
            s = s - sigma(p, x)
    s = np.where(s < 0.0, 0.0, s)
    e = np.asarray(e, dtype=float)
    out = np.interp(e, x, s)
    return np.where((e >= x[0]) & (e <= x[-1]), out, 0.0)


def bern(x):
    """B(x) = x / (e^x - 1), the Scharfetter-Gummel weight."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-3
    xs = x[small]
    out[small] = 1 - xs / 2 + xs * xs / 12 - xs ** 4 / 720
    xl = x[~small]
    with np.errstate(over="ignore"):
        out[~small] = xl / np.expm1(xl)
    return out


def _deposit(pos):
    """Cells and weights for a deposit at pos (in cell-index units, the
    centre of cell i at i): below the first centre, all in cell 0."""
    pos = np.asarray(pos, dtype=float)
    i0 = np.floor(pos).astype(int)
    w = pos - i0
    low = pos <= 0.0
    i0 = np.where(low, 0, i0)
    w = np.where(low, 0.0, w)
    return i0, 1.0 - w, i0 + 1, w


# ---------------------------------------------------------------------------
# the grid and the gas terms, assembled once
# ---------------------------------------------------------------------------
@dataclass
class Process:
    """One collision the kinetic equation carries."""
    gas: str
    kind: str                 # excitation | ionization | attachment |
    #                           assoc (v0_2: excitation + an electron at rest)
    threshold: float
    p: xsec.Process
    rows: tuple = ()          # the network rows it serves


@dataclass
class PopProcess:
    """v0_4: an electron collision with an excited population, Kr(5s) ->
    Kr(5p) and its reverse. `p` is a Process-like cross section (cm^2 on
    its energy grid) whose threshold is the energy the electron loses (the
    ledger's gap); the reverse is its Klein-Rosseland superelastic with
    g_ratio = g_lower/g_upper. The forward collision acts on the lower
    population (`name`), the reverse on the upper one (`upper`); each
    population's fraction of the colliding gas comes with the State
    (State.pops[name], State.pops[upper]). rows: the network rows
    (forward, reverse) the pair serves."""
    name: str
    p: object
    g_ratio: float
    rows: tuple = ()
    upper: str = ""


class Grid:
    def __init__(self, n: int, emax: float, procs, superelastic=None,
                 pops=()):
        self.n = n
        self.emax = emax
        self.du = du = emax / n
        self.ec = (np.arange(n) + 0.5) * du
        self.eb = np.arange(1, n) * du
        self.edges = np.arange(n + 1) * du
        self.procs = procs
        ec, eb = self.ec, self.eb
        # elastic energy transfer and momentum transfer, per gas
        self.sig_eps = {g: 2.0 * ME / M_KG[g] * elastic_mt(g, eb)
                        for g in GASES}
        sig_m = {g: elastic_mt(g, eb) for g in GASES}
        for g in GASES:                 # every process of the file, as MT
            for p in xsec.processes(g):
                if p.kind not in ("elastic", "effective"):
                    sig_m[g] = sig_m[g] + sigma(p, eb)
        self.sig_m = sig_m
        # the collision operator per gas at unit fraction
        self.op = {}
        for g in GASES:
            r, c, v = [], [], []
            for q in procs:
                if q.gas != g:
                    continue
                s = sigma(q.p, ec)
                j = np.nonzero(s)[0]
                L = GAMMA * ec[j] * s[j] * du
                r.append(j)
                c.append(j)
                v.append(L)
                if q.kind == "attachment":
                    continue
                if q.kind == "ionization":
                    pos = (ec[j] - q.threshold) / 2.0 / du - 0.5
                    mult = 2.0
                else:
                    pos = (ec[j] - q.threshold) / du - 0.5
                    mult = 1.0
                k0, w0, k1, w1 = _deposit(pos)
                r += [k0, k1]
                c += [j, j]
                v += [-mult * w0 * L, -mult * w1 * L]
                if q.kind == "assoc":
                    # v0_2: the electron the excited atom's associative
                    # ionization frees, born at rest (the first cell)
                    r.append(np.zeros_like(j))
                    c.append(j)
                    v.append(-L)
            if r:
                rr, cc, vv = (np.concatenate(r), np.concatenate(c),
                              np.concatenate(v))
                keep = rr < n
                self.op[g] = sp.csc_matrix((vv[keep], (rr[keep], cc[keep])),
                                           shape=(n, n))
            else:
                self.op[g] = sp.csc_matrix((n, n))
        # superelastic Kr* -> Kr (only when the network books it): per unit
        # population of Kr* in Kr, gamma g (e+u) sigma(e+u) du from cell j,
        # deposited at e_j + u (into the top cell above e_max)
        self.op_sup = None
        self.sig_sup_b = None
        if superelastic is not None:
            p, g_ratio = superelastic
            u = p.threshold
            s = g_ratio * (ec + u) * sigma(p, ec + u)
            j = np.nonzero(s)[0]
            L = GAMMA * s[j] * du
            k0, w0, k1, w1 = _deposit((ec[j] + u) / du - 0.5)
            rr = np.concatenate([j, np.minimum(k0, n - 1),
                                 np.minimum(k1, n - 1)])
            cc = np.concatenate([j, j, j])
            vv = np.concatenate([L, -w0 * L, -w1 * L])
            self.op_sup = sp.csc_matrix((vv, (rr, cc)), shape=(n, n))
            with np.errstate(divide="ignore", invalid="ignore"):
                self.sig_sup_b = np.where(
                    eb > 0, g_ratio * (1.0 + u / eb) * sigma(p, eb + u), 0.0)
        # v0_4: the population processes (PopProcess): the forward
        # excitation out of the lower population and its superelastic
        # reverse out of the upper one, each per unit population fraction
        # of the colliding gas, as two operators (the reverse acts on the
        # upper population, which is not the lower one's: scaling both by
        # the lower one returned 7 times too much energy to the electrons
        # at the reference point, 2026-10-08), with their momentum
        # transfer at the boundaries
        self.op_pop, self.sig_pop_b, self.pops = {}, {}, {}
        for q in pops:
            u = q.p.threshold
            s_f = sigma(q.p, ec)                              # forward
            jf = np.nonzero(s_f)[0]
            Lf = GAMMA * ec[jf] * s_f[jf] * du
            k0, w0, k1, w1 = _deposit((ec[jf] - u) / du - 0.5)
            rr = np.concatenate([jf, k0, k1])
            cc = np.concatenate([jf, jf, jf])
            vv = np.concatenate([Lf, -w0 * Lf, -w1 * Lf])
            keep = rr < n
            op_f = sp.csc_matrix((vv[keep], (rr[keep], cc[keep])),
                                 shape=(n, n))
            s_r = q.g_ratio * (ec + u) * sigma(q.p, ec + u)   # reverse
            jr = np.nonzero(s_r)[0]
            Lr = GAMMA * s_r[jr] * du
            k0, w0, k1, w1 = _deposit((ec[jr] + u) / du - 0.5)
            rr = np.concatenate([jr, np.minimum(k0, n - 1),
                                 np.minimum(k1, n - 1)])
            cc = np.concatenate([jr, jr, jr])
            vv = np.concatenate([Lr, -w0 * Lr, -w1 * Lr])
            op_r = sp.csc_matrix((vv, (rr, cc)), shape=(n, n))
            self.op_pop[q.name] = (op_f, op_r)
            with np.errstate(divide="ignore", invalid="ignore"):
                self.sig_pop_b[q.name] = (
                    sigma(q.p, eb),
                    np.where(eb > 0, q.g_ratio * (1.0 + u / eb)
                             * sigma(q.p, eb + u), 0.0))
            self.pops[q.name] = q
        # exact cell integrals of u^1/2 and u^3/2 (the e-e term with F
        # constant in each cell, ML's ee.py)
        self.w1 = 2.0 / 3.0 * np.diff(self.edges ** 1.5)
        self.w2 = 2.0 / 5.0 * np.diff(self.edges ** 2.5)
        self.sqrt_ec = np.sqrt(ec)
        # the e-e integrals with F log-linear between cell centres
        # (ee_integrals): 8-point Gauss-Legendre on each half cell, the
        # first half cell in t = sqrt(u), where u^1/2 is not smooth
        xk, wk = np.polynomial.legendre.leggauss(8)
        xk, wk = 0.5 * (xk + 1.0), 0.5 * wk
        h = 0.5 * du
        u_l = np.arange(n)[:, None] * du + xk[None, :] * h
        u_r = ec[:, None] + xk[None, :] * h
        w_l = np.tile(wk * h, (n, 1))
        t = xk * math.sqrt(h)
        u_l[0] = t * t
        w_l[0] = wk * math.sqrt(h) * 2.0 * t
        u = np.stack([u_l, u_r], axis=1)                     # (n, 2, 8)
        wq = np.stack([w_l, np.tile(wk * h, (n, 1))], axis=1)
        self._ee_off = u - ec[:, None, None]
        self._ee_w = (wq, wq * np.sqrt(u), wq * u ** 1.5)
        # the first cell's balance is replaced by F0 = 1 there (boltz_twin's
        # rule), and the normalisation int sqrt(eps) F0 d eps = 1 is made
        # after (Solver._solve_parts)
        keep = np.ones(n)
        keep[0] = 0.0
        self.keep = sp.diags(keep, format="csr")
        self.pin = sp.csr_matrix(([1.0], ([0], [0])), shape=(n, n))
        self._inj = {}

    def normalise(self, F):
        return F / float(np.sum(self.sqrt_ec * F) * self.du)

    def injection(self, e_top: float, mean: float) -> np.ndarray:
        """Deposit weights, summing to 1, of electrons born with ML's
        sub-excitation spectrum (ml/boltzmann/ee.py injection): linear on
        [0, e_top] with the given mean. The spectrum is sampled at 4000
        points and each is split between the two cell centres around it,
        as a Penning electron is, so the weights sum to 1 and their first
        moment is the mean to 2e-5 of it on grids of 12 to 40 eV (the
        electrons born below the first cell centre are put at it)."""
        key = (round(e_top, 9), round(mean, 9))
        if key not in self._inj:
            b = (e_top / 2.0 - mean) / (e_top / 6.0)
            if not -1.0 <= b <= 1.0:
                raise ValueError(f"injection mean {mean} outside [e_top/3, "
                                 f"2 e_top/3] (e_top {e_top})")
            m = 4000
            x = (np.arange(m) + 0.5) * e_top / m
            w = 1.0 + b * (1.0 - 2.0 * x / e_top)
            w = w / w.sum()
            k0, w0, k1, w1 = _deposit(x / self.du - 0.5)
            q = np.zeros(self.n)
            np.add.at(q, np.minimum(k0, self.n - 1), w * w0)
            np.add.at(q, np.minimum(k1, self.n - 1), w * w1)
            self._inj[key] = q
        return self._inj[key]

    def ee_integrals(self, F, quadrature="loglinear"):
        """A1 = int_0^eps u^1/2 F, A2 = int_0^eps u^3/2 F and A3 =
        int_eps^emax F at the inner boundaries.

        'loglinear' takes ln F linear between cell centres (the end half
        cells extrapolated, by at most a factor e), which is exact for a
        Maxwellian. For a Maxwellian of temperature T the continuous
        identity 3 T A1 = 2 (A2 + eps^3/2 A3) then holds on the grid, so
        W_ee du / D_ee = -du/T, the Scharfetter-Gummel flux of the
        discrete Maxwellian is zero, and the e-e term leaves it alone.
        'cell' takes F constant in each cell (ML's ee.py): where kT is
        not much larger than the cell, the identity fails at the lowest
        boundaries by far more than the elastic coupling to the gas
        (2m/M smaller than e-e at n_e/N ~ 1e-5), and an undriven EEDF
        condenses into the first cell."""
        du = self.du
        if quadrature == "cell":
            return (np.cumsum(self.w1 * F)[:-1], np.cumsum(self.w2 * F)[:-1],
                    (np.sum(F) - np.cumsum(F))[:-1] * du)
        lf = np.log(np.maximum(F, F.max() * 1e-300))
        s = np.diff(lf) / du
        cap = 2.0 / du
        s_l = np.concatenate([[min(max(s[0], -cap), cap)], s])
        s_r = np.concatenate([s, [min(max(s[-1], -cap), cap)]])
        slope = np.stack([s_l, s_r], axis=1)[:, :, None]
        Fq = np.exp(lf[:, None, None] + slope * self._ee_off)
        i0, i1, i2 = ((w * Fq).sum(axis=(1, 2)) for w in self._ee_w)
        return (np.cumsum(i1)[:-1], np.cumsum(i2)[:-1],
                (i0.sum() - np.cumsum(i0))[:-1])

    def mean(self, F) -> float:
        return float(np.sum(self.ec ** 1.5 * F) * self.du)

    def maxwellian(self, te: float) -> np.ndarray:
        return self.normalise(np.exp(-self.ec / max(te, 1e-4)))

    def maxwellian_at_mean(self, mean: float) -> np.ndarray:
        """The grid's Maxwellian with mean energy `mean` [eV] on this grid
        (a cell as wide as kT shifts the grid's own mean from (3/2) kT by
        several percent)."""
        lo, hi = 0.2 * mean, 5.0 * mean
        te = brentq(lambda t: self.mean(self.maxwellian(t)) - mean, lo, hi,
                    xtol=1e-15, rtol=4.0 * np.finfo(float).eps)
        return self.maxwellian(te)

    def remap(self, ec_old, F_old) -> np.ndarray:
        """F from another grid onto this one, keeping its normalisation and
        its mean energy: interpolated, then its energy axis stretched until
        the mean is the old grid's (a time step from F reads that mean as
        its starting energy, and a remap that moved it would be heating
        or cooling of the step's own)."""
        if ec_old.shape == self.ec.shape and np.array_equal(ec_old, self.ec):
            return self.normalise(np.maximum(F_old, 0.0))
        mean = float(np.sum(ec_old ** 1.5 * F_old)
                     / np.sum(np.sqrt(ec_old) * F_old))

        def at(r):
            return self.normalise(np.maximum(
                np.interp(self.ec * r, ec_old, F_old, right=0.0), 0.0))
        r = 1.0
        for _ in range(40):
            m = self.mean(at(r))
            if abs(m / mean - 1.0) < 1e-13:
                break
            r *= m / mean
        return at(r)

    def rate(self, p: xsec.Process, F) -> float:
        """<sigma v> [m^3/s] over F."""
        return GAMMA * float(np.sum(self.ec * sigma(p, self.ec) * F)) \
            * self.du

    def superelastic_rate(self, p: xsec.Process, g_ratio: float, F) -> float:
        """<sigma_sup v> [m^3/s] over F, per de-exciting partner: the
        Klein-Rosseland cross section g (e+u)/e sigma(e+u), as the
        collision operator carries it."""
        u = p.threshold
        return GAMMA * float(np.sum(g_ratio * (self.ec + u)
                                    * sigma(p, self.ec + u) * F)) * self.du

    def pop_rates(self, name: str, F) -> tuple:
        """(forward, reverse) rates [m^3/s] of a population process over
        F, per unit population (v0_4)."""
        q = self.pops[name]
        return (self.rate(q.p, F), self.superelastic_rate(q.p, q.g_ratio, F))


# ---------------------------------------------------------------------------
# one state of the plasma, and the solution there
# ---------------------------------------------------------------------------
@dataclass
class State:
    """What the kinetic equation needs from LIGHT at one moment."""
    x: dict                    # collision fractions {'Kr': 1.0}
    tg: float                  # K
    n_col: float               # colliding gas density [m^-3]
    ne: float                  # electron density [m^-3]
    pen: tuple = ()            # ((frequency per electron [1/s], eV), ...)
    rec: float = 0.0           # recombination template C [Kr2+] [m^3/s eV^-n]
    rec_pow: float = -0.5
    kr_star: float = 0.0       # [Kr*]/[Kr] (superelastic, when booked)
    # electrons born with a spectrum, a beam's secondaries: ((frequency per
    # electron [1/s], e_top [eV], mean [eV]), ...), linear on [0, e_top]
    # (Grid.injection)
    inj: tuple = ()
    # v0_4 (2026-10-08): the cations [m^-3], whose Coulomb collisions
    # enter the momentum transfer (ei), Debye-screened for the DC field
    # that stands in for the given power and with Turnbull et al.'s
    # inverse-bremsstrahlung logarithm (ib.py) for an EM wave of angular
    # frequency omega [rad/s] and rms reduced field en_em_td [Td]. The
    # model's error band (ib_scale) acts on the electron-ion share of the
    # absorbed power LIGHT reads (model.py), not on the cross section
    n_ion: float = 0.0
    ei: bool = True
    omega: float = 0.0
    en_em_td: float = 0.0
    # v0_4: the excited populations the population processes act on, as
    # fractions of the colliding gas {species: n / n_col}, the lower and
    # the upper state of each pair
    pops: dict = field(default_factory=dict)

    def source_rate(self) -> float:
        """The electrons the sources add, per electron [1/s]."""
        return (sum(v for v, _ in self.pen)
                + sum(v for v, _, _ in self.inj))

    def source_power(self) -> float:
        """The energy they bring, per electron [eV/s]."""
        return (sum(v * e for v, e in self.pen)
                + sum(v * m for v, _, m in self.inj))


@dataclass
class Solution:
    F: np.ndarray
    en_td: float
    nu: float                  # growth frequency per electron [1/s]
    mean: float
    growth_off: bool
    iters: int
    grid: Grid = field(repr=False)
    # False where the field fell to zero with the mean still above the
    # target (the target lies below the no-field mean)
    reached: bool = True

    @property
    def te(self) -> float:
        return 2.0 / 3.0 * self.mean


class Solver:
    """The kinetic EEDF for the processes a network books.

    procs: the Process list; superelastic: (xsec.Process, g_ratio) of the
    Kr* -> Kr collision when the network books it (row SUPER), else None."""

    def __init__(self, procs, n_cells=600, emax=40.0, superelastic=None,
                 tol=1e-9, maxit=400, ee_quadrature="loglinear",
                 elastic=("Kr",), pops=()):
        self.procs = procs
        self.superelastic = superelastic
        self.pops = tuple(pops)                # v0_4: PopProcess list
        self.ee_quadrature = ee_quadrature     # Grid.ee_integrals
        # the gases whose elastic energy exchange the equation carries (a
        # network books its own: LIGHT's G6 is Kr's alone). Every gas's
        # momentum transfer stays in sigma~
        self.elastic = tuple(elastic)
        self.n_cells = n_cells
        self.tol = tol
        self.maxit = maxit
        self._grids = {}
        self.emax = emax
        self.grid = self._grid(emax)
        self.trace = None      # a list: solve_joint appends every step to it

    def _grid(self, emax):
        key = round(emax, 6)
        if key not in self._grids:
            self._grids[key] = Grid(self.n_cells, emax, self.procs,
                                    self.superelastic, self.pops)
        return self._grids[key]

    # ---------------------------------------------------------- one solve
    def _frequencies(self, g: Grid, s: State, F):
        """(ionization, attachment, recombination) per electron [1/s].
        Ionization counts the associative ionizations (v0_2)."""
        ion = att = 0.0
        for q in self.procs:
            if q.kind in ("ionization", "attachment", "assoc"):
                k = g.rate(q.p, F) * s.x[q.gas] * s.n_col
                if q.kind == "attachment":
                    att += k
                else:
                    ion += k
        rec = 0.0
        if s.rec > 0.0:
            rec = s.rec * float(np.sum(g.ec ** (s.rec_pow + 0.5) * F)) * g.du
        return ion, att, rec

    def _parts(self, g: Grid, s: State, nu, F, ee=True, mem=0.0):
        """The matrix without the field: the collision operator, the
        recombination and growth sinks on the diagonal, and the flux
        coefficients W and D0 (elastic, thermal and e-e, the e-e ones at
        F). The field adds (E/N)^2 dfield to D (_assemble). mem = 1/dt
        adds a time step's diagonal, the growth term's with nu -> 1/dt
        (its source, F_prev/dt, is solve_joint's)."""
        n, du, eb = g.n, g.du, g.eb
        A = sp.csc_matrix((n, n))
        sig_eps = np.zeros(n - 1)
        sig_t = np.zeros(n - 1)
        for gas, x in s.x.items():
            if x <= 0.0:
                continue
            A = A + x * g.op[gas]
            if gas in self.elastic:
                sig_eps += x * g.sig_eps[gas]
            sig_t += x * g.sig_m[gas]
        if g.op_sup is not None and s.kr_star > 0.0:
            ws = s.x["Kr"] * s.kr_star
            A = A + ws * g.op_sup
            sig_t += ws * g.sig_sup_b
        for name, q in g.pops.items():             # v0_4
            f_lo = s.pops.get(name, 0.0)
            f_up = s.pops.get(q.upper, 0.0)
            op_f, op_r = g.op_pop[name]
            sig_f, sig_r = g.sig_pop_b[name]
            if f_lo > 0.0:
                A = A + f_lo * op_f
                sig_t += f_lo * sig_f
            if f_up > 0.0:
                A = A + f_up * op_r
                sig_t += f_up * sig_r
        nu_n = nu / s.n_col
        sig_t = sig_t + nu_n / (GAMMA * np.sqrt(eb))
        # v0_4: the ions' Coulomb collisions in the momentum transfer. The
        # DC field (the stand-in for the given power) sees them
        # Debye-screened; the EM wave sees Turnbull et al.'s
        # inverse-bremsstrahlung cross section (ib.py), the atoms and the
        # growth term once in each
        sig_hf_atoms = sig_t
        sig_hf = sig_t
        if s.ei and s.n_ion > 0.0:
            x_ion = s.n_ion / s.n_col
            sig_t = sig_t + x_ion * _sigma_ei_dc(eb, s.ne, g.mean(F))
            if s.omega > 0.0 and s.en_em_td > 0.0:
                from . import ib
                sig_hf = sig_hf + x_ion * ib.sigma_ei(eb, s.omega, s.ne,
                                                      s.n_ion)
        hf = (None, 0.0, None)
        if s.omega > 0.0 and s.en_em_td > 0.0:
            # Hagelaar & Pitchford's field term at omega: D_hf = (gamma/3)
            # (E_rms/N)^2 eps sigma~ / (sigma~^2 + (omega/N)^2 / (gamma^2
            # eps)), the DC one at omega = 0
            a2 = (s.omega / s.n_col) ** 2 / (GAMMA ** 2 * eb)
            hf = (GAMMA / 3.0 * eb * sig_hf / (sig_hf ** 2 + a2),
                  s.en_em_td,
                  GAMMA / 3.0 * eb * sig_hf_atoms / (sig_hf_atoms ** 2 + a2))
        W = -GAMMA * eb ** 2 * sig_eps
        D0 = GAMMA * (KB * s.tg / QE) * eb ** 2 * sig_eps
        if ee and s.ne > 0.0:
            alpha = (s.ne / s.n_col) * C_EE * _coulomb_log(s.ne, g.mean(F)) \
                * GAMMA
            a1, a2, a3 = g.ee_integrals(F, self.ee_quadrature)
            W = W - 3.0 * alpha * a1
            D0 = D0 + 2.0 * alpha * (a2 + eb ** 1.5 * a3)
        diag = np.zeros(n)
        if s.rec > 0.0:
            diag += s.rec / s.n_col * g.ec ** (s.rec_pow + 0.5) * du
        diag += (nu_n + mem / s.n_col) * g.sqrt_ec * du
        return A + sp.diags(diag, format="csc"), W, D0, \
            GAMMA / 3.0 * eb / sig_t, hf

    @staticmethod
    def _d_total(parts, en_td):
        """D at the DC field E/N [Td]: thermal, e-e, the DC field's part and
        the EM wave's (parts[4] = (dfield_hf, en_em_td, dfield_hf without
        the ions))."""
        _, _, D0, dfield, hf = parts
        D = D0 + (en_td * TD) ** 2 * dfield
        if hf[0] is not None:
            D = D + (hf[1] * TD) ** 2 * hf[0]
        return D

    def _assemble(self, g: Grid, parts, en_td):
        """The matrix at E/N [Td] from _parts: Scharfetter-Gummel fluxes
        through the inner boundaries."""
        A0, W, D0, dfield, hf = parts
        n, du = g.n, g.du
        D = self._d_total(parts, en_td)
        z = W * du / D
        a = D / du * bern(-z)
        b = -D / du * bern(z)
        i = np.arange(n - 1)
        diag = np.zeros(n)
        diag[:-1] += a
        diag[1:] -= b
        flux = sp.csc_matrix(
            (np.concatenate([b, -a]), (np.concatenate([i, i + 1]),
                                       np.concatenate([i + 1, i]))),
            shape=(n, n))
        return A0 + flux + sp.diags(diag, format="csc")

    def _matrix(self, g: Grid, s: State, en_td, nu, F, ee=True):
        return self._assemble(g, self._parts(g, s, nu, F, ee), en_td)

    def _source(self, g: Grid, s: State):
        q = np.zeros(g.n)
        for nu_p, e in s.pen:
            if nu_p <= 0.0:
                continue
            k0, w0, k1, w1 = _deposit(np.array([e / g.du - 0.5]))
            q[int(k0[0])] += nu_p / s.n_col * float(w0[0])
            if int(k1[0]) < g.n:
                q[int(k1[0])] += nu_p / s.n_col * float(w1[0])
        for nu_i, e_top, mean in s.inj:
            if nu_i > 0.0:
                q = q + nu_i / s.n_col * g.injection(e_top, mean)
        return q

    def solve_at(self, s: State, en_td: float, F0=None, nu0=None,
                 growth=True):
        """The EEDF at a given E/N [Td], e-e iterated to self-consistency
        with the growth frequency. Returns a Solution."""
        g = self.grid
        F = g.maxwellian(max(1.5 * KB * s.tg / QE, 0.03)) if F0 is None \
            else g.normalise(F0)
        q = self._source(g, s)
        ion, att, rec = self._frequencies(g, s, F)
        pen = s.source_rate()
        nu = (ion + pen - att - rec) if growth else 0.0
        if nu0 is not None and growth:
            nu = nu0
        change = np.inf
        for it in range(1, self.maxit + 1):
            Fn = self._one(g, s, F, nu, en_td, q)
            Fn = g.normalise(0.5 * (Fn + F))
            change = float(np.max(np.abs(Fn - F)) / np.max(np.abs(F)))
            F = Fn
            if growth:
                ion, att, rec = self._frequencies(g, s, F)
                nu_new = ion + pen - att - rec
                dnu = abs(nu_new - nu) / max(abs(nu_new), ion + pen + att
                                            + rec, 1.0)
                nu = 0.5 * (nu + nu_new)
            else:
                dnu = 0.0
            if change < self.tol and dnu < self.tol:
                break
        else:
            # the last resort: a slow contraction (a strong net loss with
            # no heating) still gives an EEDF good to its last change
            if max(change, dnu) > 1e-4:
                raise RuntimeError(f"kinetic EEDF: no convergence at "
                                   f"{en_td:.3g} Td ({change:.1e})")
        return Solution(F, en_td, nu, g.mean(F), not growth, it, g)

    def solve(self, s: State, en_td: float, F0=None, nu0=None, growth=True):
        """The EEDF at a given E/N: the joint solve with the field held,
        the damped Picard iteration (solve_at) where it does not converge,
        and without the growth term where neither converges with it or it
        makes F negative (a strong net loss, attachment outrunning energy
        relaxation, the two-term scheme's known failure)."""
        sol = None
        if growth:
            sol = self.solve_joint(s, None, F0, nu0, en_fixed=en_td)
            if sol is None:
                try:
                    sol = self.solve_at(s, en_td, F0, nu0, growth=True)
                except RuntimeError:
                    sol = None
            if sol is not None and sol.F.min() < -1e-8 * sol.F.max():
                sol = None
        if sol is None:
            sol = self.solve_joint(s, None, F0, None, en_fixed=en_td,
                                   growth=False)
            if sol is None:
                sol = self.solve_at(s, en_td, F0, growth=False)
        if sol.F.min() < -1e-8 * sol.F.max():
            raise RuntimeError(f"kinetic EEDF: negative at {en_td:.3g} Td "
                               "even without the growth term")
        # round-off far down the tail (1e-30 of the peak) would give a
        # process that lives there a negative rate
        sol.F = np.maximum(sol.F, 0.0)
        return sol

    # ------------------------------------ the joint solve (the fast route)
    S_MEAN = 0.6          # d ln <eps> / d ln E/N before the first estimate
    EN_MIN, EN_MAX = 1e-4, 3000.0
    REACH_TOL = 1e-3      # a no-field mean this far above the target counts
    MEAN_TOL = 1e-10      # |ln <eps> - ln target| of the field's inner solve

    def _solve_parts(self, g, parts, en_td, q):
        """One linear solve with the first cell's balance replaced by the
        normalisation int sqrt(eps) F0 d eps = 1. The sources q are per
        electron, so with the growth term at its self-consistent nu the
        replaced balance holds by itself, and without the growth term it
        takes up the net gain or loss of electrons (boltz_twin's rule: F
        stays the population LIGHT has, whatever the sources and sinks do
        to its size). Done as boltz_twin pins the first cell, F0 = 1
        there, which keeps the matrix sparse and the far tail clean: one
        factorisation, the pinned solution Fh and the sources' Fp (first
        cell 0), and F = c Fh + Fp with c from the normalisation."""
        A = (g.keep @ self._assemble(g, parts, en_td) + g.pin).tocsc()
        lu = splu(A)
        e0 = np.zeros(g.n)
        e0[0] = 1.0
        Fh = lu.solve(e0)
        if not np.any(q[1:]):
            return g.normalise(Fh)
        qp = q.copy()
        qp[0] = 0.0
        Fp = lu.solve(qp)
        w = g.sqrt_ec * g.du
        c = (1.0 - float(np.sum(w * Fp))) / float(np.sum(w * Fh))
        return g.normalise(c * Fh + Fp)

    def _one(self, g, s, F, nu, en_td, q, mem=0.0):
        """One linear solve with the matrix assembled at F (its e-e
        coefficients) and nu."""
        return self._solve_parts(g, self._parts(g, s, nu, F, mem=mem), en_td,
                                 q)

    def _power_match(self, g, s, F, nu, q, lg, lp, tol, mem=0.0):
        """The field at which the linear problem assembled at (F, nu) takes
        the Joule power e^lp per electron [eV/s], as the solver takes it
        (_split): secant steps on ln P against ln E/N from lg (the power
        goes nearly as the field squared), then a bracket and Brent's
        method. Returns (Fn, ln E/N, at_floor), with Fn None where a solve
        is not finite."""
        parts = self._parts(g, s, nu, F, mem=mem)
        _, W, D0, dfield, _hf = parts
        lg_lo, lg_hi = math.log(self.EN_MIN), math.log(self.EN_MAX)
        seen = {}

        def h(x):
            if x not in seen:
                Fn = self._solve_parts(g, parts, math.exp(x), q)
                if not np.all(np.isfinite(Fn)):
                    raise _NoSolution
                en_v = math.exp(x) * TD
                c = _sg_factor(W, self._d_total(parts, math.exp(x)), g.du)
                p = -en_v ** 2 * s.n_col * float(np.sum(
                    dfield * c * (Fn[1:] - Fn[:-1])))
                seen[x] = (Fn, math.log(p) - lp if p > 0.0 else -50.0)
            return seen[x][1]
        try:
            x0 = min(max(lg, lg_lo), lg_hi)
            h0 = h(x0)
            for _ in range(2):
                if abs(h0) < tol:
                    return seen[x0][0], x0, False
                x1 = min(max(x0 - 0.5 * h0, x0 - 2.0, lg_lo), x0 + 2.0,
                         lg_hi)
                if x1 == x0:
                    break
                x0, h0 = x1, h(x1)
            if abs(h0) < tol:
                return seen[x0][0], x0, False
            above = [x for x, v in seen.items() if v[1] > 0.0]
            if not above:
                if h(lg_hi) <= 0.0:
                    return seen[lg_hi][0], lg_hi, False
                above = [lg_hi]
            b = min(above)
            below = [x for x, v in seen.items() if v[1] < 0.0 and x < b]
            if not below:
                if h(lg_lo) >= 0.0:
                    return seen[lg_lo][0], lg_lo, True
                below = [lg_lo]
            root = brentq(h, max(below), b, xtol=1e-14,
                          rtol=4.0 * np.finfo(float).eps, maxiter=200)
            h(root)
            return seen[root][0], root, False
        except _NoSolution:
            return None, lg, False

    def _mean_match(self, g, s, F, nu, q, lg, lt, slope, tol, mem=0.0):
        """The field at which the linear problem assembled at (F, nu) gives
        the mean energy e^lt. With its coefficients held the mean rises
        with the field, so: two secant steps on ln <eps> against ln E/N
        from lg (enough from a warm start), then a bracket (where the mean
        is still too high, the floor itself: if even no field is too hot
        the answer is the floor, and where it is too low, steps of e^2 up)
        and Brent's method in it. The secant's slope can be useless (near
        no field, where the Penning electrons and not the field set the
        mean, it is 1e-6), and steps with it never arrive. Returns (Fn,
        ln E/N, slope, at_floor), with Fn None where a solve is not finite
        or has no positive mean."""
        parts = self._parts(g, s, nu, F, mem=mem)
        lg_lo, lg_hi = math.log(self.EN_MIN), math.log(self.EN_MAX)
        seen = {}

        def h(x):
            if x not in seen:
                Fn = self._solve_parts(g, parts, math.exp(x), q)
                mn = g.mean(Fn)
                if not (np.all(np.isfinite(Fn)) and mn > 0.0):
                    raise _NoSolution
                seen[x] = (Fn, math.log(mn) - lt)
            return seen[x][1]
        try:
            x0 = min(max(lg, lg_lo), lg_hi)
            h0 = h(x0)
            for _ in range(2):
                if abs(h0) < tol:
                    return seen[x0][0], x0, slope, False
                x1 = min(max(x0 - h0 / slope, x0 - 2.0, lg_lo), x0 + 2.0,
                         lg_hi)
                if x1 == x0:
                    break
                h1 = h(x1)
                est = (h1 - h0) / (x1 - x0)
                if est > 1e-3:
                    slope = est
                x0, h0 = x1, h1
            if abs(h0) < tol:
                return seen[x0][0], x0, slope, False
            below = [x for x, v in seen.items() if v[1] < 0.0]
            above = [x for x, v in seen.items() if v[1] > 0.0]
            if not below:
                if h(lg_lo) >= 0.0:
                    return seen[lg_lo][0], lg_lo, slope, True
                below = [lg_lo]
            x = max(above) if above else max(below)
            while not above:
                if x >= lg_hi:
                    return seen[x][0], x, slope, False
                x = min(x + 2.0, lg_hi)
                if h(x) > 0.0:
                    above = [x]
            a, b = max(below), min(above)
            root = brentq(h, a, b, xtol=1e-14, rtol=4.0 * np.finfo(float).eps,
                          maxiter=200)
            h(root)
            return seen[root][0], root, slope, False
        except _NoSolution:
            return None, lg, slope, False

    def solve_joint(self, s: State, target_mean=None, F0=None, nu0=None,
                    en0=None, growth=True, m=6, maxit=200, en_fixed=0.0,
                    max_resets=12, memory=None, target_power=None):
        """F and the growth frequency as one fixed point, accelerated by
        Anderson mixing. With target_mean, every step solves the linear
        problem at the field that gives exactly that mean energy
        (_mean_match), so the iteration only has the shape to settle: at a
        fixed field the problem can have two solutions (Kr's Ramsauer
        minimum with strong e-e, the thermal runaway), at a fixed mean it
        has one. Without target_mean the field is en_fixed [Td]. m = 0 is
        the plain iteration, damped by half. Returns a Solution, or None
        where it does not converge in maxit or loses positivity (the
        caller then drops the growth term, or without it takes the plain
        iteration). With the growth term, a field at its floor for five
        steps with the mean still above the target ends the solve early
        and unconverged, with reached=False: the growth term cannot give
        that mean energy.

        memory = (energies, F_prev, dt) makes the solve one backward-Euler
        step of the time-dependent equation from F_prev over dt [s]: the
        steady problem plus sqrt(eps) (F - F_prev)/(N dt), a diagonal like
        the growth term's and a source. With target_power [eV/s per
        electron] in place of target_mean the field is the one whose Joule
        power that is."""
        g = self.grid
        mem, Fprev = 0.0, None
        if memory is not None:
            ec_p, F_p, dt = memory
            mem = 1.0 / dt
            Fprev = g.remap(ec_p, F_p)
        F = None
        if F0 is not None:
            F = g.normalise(np.maximum(F0, 0.0))
            if target_mean is not None:
                ratio = g.mean(F) / target_mean
                if abs(math.log(ratio)) > 0.3:
                    # far from the target mean a warm start is a worse
                    # start than a Maxwellian at it
                    F = None
                elif abs(math.log(ratio)) > 1e-3:
                    # its energy axis stretched onto the target mean
                    F = g.normalise(np.interp(g.ec * ratio, g.ec, F,
                                              right=0.0))
        if F is None and Fprev is not None:
            F = Fprev.copy()
        if F is None:
            F = g.maxwellian(2.0 / 3.0 * target_mean) \
                if target_mean is not None \
                else g.maxwellian(max(KB * s.tg / QE, 0.02))
        q = self._source(g, s)
        if Fprev is not None:
            q = q + mem / s.n_col * g.sqrt_ec * g.du * Fprev
        pen = s.source_rate()
        ion, att, rec = self._frequencies(g, s, F)
        nu_scale = max(ion + pen + att + rec, 1.0)
        # a cold start takes the net of the sources and the guessed F's
        # frequencies where it is a gain (a beam's secondaries in a seed of
        # electrons, where nu = 0 would leave the sources nothing to
        # balance), else 0 (as ML's solver does: nu from a guessed F can be
        # far too negative, F2 attaching thermal electrons fast)
        if not growth:
            nu = 0.0
        elif nu0 is not None:
            nu = nu0
        else:
            nu = max(0.0, ion + pen - att - rec)
        field = target_mean is not None or target_power is not None
        lg_lo, lg_hi = math.log(self.EN_MIN), math.log(self.EN_MAX)
        lg = math.log(min(max(en0 or 10.0, self.EN_MIN), self.EN_MAX)) \
            if field else 0.0
        lt = math.log(target_mean) if target_mean is not None else 0.0
        lp = math.log(target_power) if target_power is not None else 0.0
        w = math.sqrt(g.n)
        slope = self.S_MEAN
        tol_in = 1e-6
        floor = 0
        xs, gs = [], []
        resets = 0
        damp = None            # the plain step's size after a reset
        for it in range(1, maxit + 1):
            if target_power is not None:
                Fn, lgn, at_floor = self._power_match(
                    g, s, F, nu, q, lg, lp, tol_in, mem)
            elif field:
                Fn, lgn, slope, at_floor = self._mean_match(
                    g, s, F, nu, q, lg, lt, slope, tol_in, mem)
            else:
                Fn, lgn, at_floor = self._one(g, s, F, nu, en_fixed, q,
                                              mem), 0.0, False
            if Fn is None or not np.all(np.isfinite(Fn)) \
                    or Fn.min() < -1e-8 * Fn.max() or g.mean(Fn) <= 0.0:
                # positivity lost: forget the mixing history and go on
                # with small plain steps (damp). A net loss too large (a nu
                # from an F far from the solution, attachment outrunning
                # the energy relaxation) is pulled toward zero. A gain too
                # small (sources outrunning every sink, so the first cell
                # would take up a net gain) is raised, to the sources' rate
                # at least
                resets += 1
                if resets > max_resets or not growth:
                    return None
                xs.clear()
                gs.clear()
                nu = 0.5 * nu if nu < 0.0 else max(2.0 * nu, pen)
                if nu == 0.0:
                    return None
                damp = 0.3
                continue
            if growth:
                ion, att, rec = self._frequencies(g, s, Fn)
                nun = ion + pen - att - rec
                nu_scale = max(nu_scale, ion + pen + att + rec)
            else:
                nun = 0.0
            if growth and at_floor and target_mean is not None \
                    and math.log(g.mean(Fn)) > lt + 0.01:
                floor += 1
                if floor >= 5:
                    return Solution(Fn, 0.0, nun, g.mean(Fn), False, it, g,
                                    reached=False)
            else:
                floor = 0
            fm = max(float(np.max(F)), 1e-300)
            change = max(float(np.max(np.abs(Fn - F))) / fm,
                         abs(nun - nu) / nu_scale, abs(lgn - lg))
            if self.trace is not None:
                self.trace.append((it, change, g.mean(Fn), lgn, nun))
            # the inner match needs no more precision than the outer step
            tol_in = min(1e-6, max(self.MEAN_TOL, 1e-3 * change))
            if change < self.tol:
                F, nu, lg = Fn, nun, lgn
                break
            x = np.concatenate([F / fm, [w * nu / nu_scale]])
            gx = np.concatenate([Fn / fm, [w * nun / nu_scale]])
            if damp is not None:
                # after a reset nu creeps toward its own consistent value
                # instead of jumping back to the one that broke positivity;
                # the mixing resumes once the steps are small
                new = x + damp * (gx - x)
                if change < 1e-3:
                    damp = None
            elif m == 0:
                new = 0.5 * (x + gx)          # the plain iteration, damped
            else:
                xs.append(x)
                gs.append(gx)
                if len(xs) > m + 1:
                    xs.pop(0)
                    gs.pop(0)
                new = gx
                if len(xs) > 1:
                    R = np.stack([gi - xi for gi, xi in zip(gs, xs)],
                                 axis=1)
                    dR = R[:, 1:] - R[:, :-1]
                    G = np.stack(gs, axis=1)
                    dG = G[:, 1:] - G[:, :-1]
                    coef, *_ = np.linalg.lstsq(dR, gx - x, rcond=None)
                    # multiply and sum, not matmul: numpy 2 on macOS
                    # Accelerate warns spuriously on underflowed zeros
                    # (LIGHT Known-bad)
                    new = gx - (dG * coef[None, :]).sum(axis=1)
                    if not np.all(np.isfinite(new)):
                        new = gx
            F = np.maximum(new[:g.n], 0.0) * fm
            if not np.any(F > 0.0):
                return None
            F = g.normalise(F)
            nu = (min(max(new[g.n] * nu_scale / w, -2.0 * nu_scale),
                      2.0 * nu_scale) if growth else 0.0)
            lg = lgn
        else:
            return None
        en = math.exp(lg) if field else en_fixed
        reached = True
        if field and lg <= lg_lo + 1e-9:
            en = 0.0                      # the target is below the no-field
            reached = target_mean is None or \
                g.mean(F) <= target_mean * (1.0 + self.REACH_TOL)
        return Solution(F, en, nu, g.mean(F), not growth, it, g,
                        reached=reached)

    GROWTH_LIMIT = 10.0

    def growth_valid(self, s: State, sol: Solution) -> bool:
        """The growth term holds the EEDF's shape steady while n_e changes
        at nu. With no heating at all (no field and no Penning electrons
        yet, as at the start of a run) and electrons lost to attachment
        at thermal energies, the two-term swarm solution is an artifact, a
        hot EEDF with nothing to heat it. The term is dropped where the
        energy the population turnover moves, |nu| <eps>, exceeds
        GROWTH_LIMIT times the heating per electron (field and sources),
        unless the turnover is GROWTH_LIMIT times slower than the energy
        relaxation by the gas (relax_freq), which then sets the shape
        whatever the term does (a recombining afterglow, where dropping it
        would put the recombined electrons back cold, below the gas
        temperature). Where turnover and heating are comparable (the
        afterglow, where Penning electrons heat and recombination removes)
        it stays, since without it the kinetic energy balance does not
        close."""
        if sol.nu >= 0.0:
            return True
        if -sol.nu * self.GROWTH_LIMIT <= self.relax_freq(s, sol):
            return True
        heating = (self.field_power(s, sol) + self.em_power(s, sol)
                   + s.source_power())
        return -sol.nu * sol.mean <= self.GROWTH_LIMIT * heating

    def relax_freq(self, s: State, sol: Solution) -> float:
        """The elastic energy-relaxation frequency per electron [1/s],
        N sum_g x_g (2m/M_g) <sigma_m v> over the solution, for the gases
        whose elastic exchange the equation carries."""
        g = sol.grid
        fb = 0.5 * (sol.F[1:] + sol.F[:-1])
        return s.n_col * GAMMA * g.du * sum(
            s.x[gas] * float(np.sum(g.eb * g.sig_eps[gas] * fb))
            for gas in self.elastic if s.x.get(gas, 0.0) > 0.0)

    def eedf(self, s: State, target_mean=None, warm=None, growth=True,
             adapt=True):
        """The kinetic EEDF LIGHT reads at one state, warm-started from
        warm (a previous Solution).

        With the growth term where it converges, gives LIGHT's mean energy
        (reached) and is valid (growth_valid), else without it. After a
        solution without it, the growth term is tried again, on a budget,
        only where the solution without it would pass growth_valid
        (_eedf_on_grid). With adapt, the grid is stretched or shortened so
        that F
        falls 11 to 25 decades across it (_emax_for). Without (the table
        points around a centre, which share its grid), it only stretches,
        where F falls fewer than 8 decades. The growth decision made on
        the first grid holds on the next. A grid already tried is not tried
        again: the best solution so far is returned, and the solver's grid
        is always the returned solution's."""
        lo, hi = (11.0, 25.0) if adapt else (8.0, np.inf)
        tried, best = set(), None
        for _ in range(8):
            tried.add(round(self.grid.emax, 6))
            sol = self._eedf_on_grid(s, target_mean, warm, growth)
            growth = growth and not sol.growth_off
            sol.F = np.maximum(sol.F, 0.0)
            decades = _decades(sol.F)
            # the best: in the window, else long enough, else the longest
            key = (lo <= decades <= hi, decades >= lo, -abs(decades - 15.0))
            if best is None or key > best[0]:
                best = (key, sol)
            if lo <= decades <= hi:
                break
            emax = self._emax_for(sol.grid, sol.F, decades, lo)
            if round(emax, 6) in tried:
                break
            self.grid = self._grid(emax)
            warm = sol
        sol = best[1]
        self.grid = sol.grid
        return sol

    def step(self, s: State, prev: Solution, dt: float, target_mean=None,
             target_power=None, adapt=True):
        """One backward-Euler step of the time-dependent EEDF from the
        Solution prev over dt [s]. The growth term is the population's own
        change, nu from F, with no validity rule: in a step it moves the
        shape by nu dt at most, where a steady solution has to absorb all
        of it. The field gives target_mean or target_power (else none),
        and where even no field leaves the mean above target_mean the
        no-field step stands, not reached. On the solver's grid, stretched
        or shortened as eedf does (_emax_for)."""
        memory = (prev.grid.ec, prev.F, dt)
        lo, hi = (11.0, 25.0) if adapt else (8.0, np.inf)
        tried, best, warm = set(), None, prev
        for _ in range(8):
            tried.add(round(self.grid.emax, 6))
            F0 = np.interp(self.grid.ec, warm.grid.ec, warm.F, right=0.0)
            nu0 = None if warm.growth_off else warm.nu
            kw = dict(memory=memory, target_power=target_power)
            sol = self.solve_joint(s, target_mean, F0, nu0,
                                   warm.en_td or None, True, **kw)
            if sol is None:
                sol = self.solve_joint(s, target_mean, F0, nu0,
                                       warm.en_td or None, True, m=0,
                                       maxit=3000, **kw)
            if sol is None:
                raise RuntimeError(
                    f"kinetic EEDF: no time step of {dt:.3g} s to mean "
                    f"{target_mean} eV (n_e {s.ne:.3g} m^-3, Tg {s.tg:.4g} "
                    f"K, x {s.x})")
            if not sol.reached:
                z = self.solve_joint(s, None, sol.F, sol.nu, en_fixed=0.0,
                                     memory=memory)
                if z is not None:
                    z.reached = False
                    sol = z
            sol.F = np.maximum(sol.F, 0.0)
            decades = _decades(sol.F)
            key = (lo <= decades <= hi, decades >= lo, -abs(decades - 15.0))
            if best is None or key > best[0]:
                best = (key, sol)
            if lo <= decades <= hi:
                break
            emax = self._emax_for(sol.grid, sol.F, decades, lo)
            if round(emax, 6) in tried:
                break
            self.grid = self._grid(emax)
            warm = sol
        sol = best[1]
        self.grid = sol.grid
        return sol

    EMAX_MIN, EMAX_MAX, EMAX_STEP = 12.0, 400.0, 1.1

    def _emax_for(self, g, F, decades, lo, want=15.0):
        """The grid top at which F falls about `want` decades below its
        peak. On a grid long enough, where F reaches that fall; on one too
        short, the fall across the grid's top tenth extrapolated (at most
        a factor 3). The fall is slow in the bulk and steep past the
        inelastic thresholds, so e_max is read off F, not scaled.
        Rounded up onto the lattice 12 x 1.1^k eV, so that each grid is
        built once."""
        lf = np.log10(np.maximum(F, F.max() * 1e-300) / F.max())
        ipk = int(np.argmax(F))
        past = np.nonzero(lf[ipk:] < -want)[0]
        if decades >= lo and past.size:
            e = float(g.ec[ipk + past[0]])
        else:
            k = max(g.n // 10, 2)
            fall = (lf[-k] - lf[-1]) / (g.ec[-1] - g.ec[-k])   # decades/eV
            e = g.emax + (want - decades) / max(fall, 1e-6)
            e = min(max(e, g.emax * 1.1), 3.0 * g.emax)
        k = math.ceil(math.log(max(e, self.EMAX_MIN) / self.EMAX_MIN)
                      / math.log(self.EMAX_STEP) - 1e-9)
        return min(self.EMAX_MIN * self.EMAX_STEP ** k, self.EMAX_MAX)

    def _without_growth(self, s, target_mean, F0, en0):
        """The solution without the growth term: the joint solve, and the
        plain damped iteration (slower, surer) where the mixing fails.
        None where neither converges: where the sources outrun every sink
        (a beam's secondaries in a seed of electrons), the first cell would
        have to take up a net gain, and the population goes negative."""
        sol = self.solve_joint(s, target_mean, F0, None, en0, False)
        if sol is None:
            sol = self.solve_joint(s, target_mean, F0, None, en0, False,
                                   m=0, maxit=3000)
        return sol

    # a growth attempt from a solution without the growth term, or from a
    # cold start: at most this many steps and positivity resets (a
    # converging one takes 40 from there, a failing one would otherwise run
    # the full reset sequence)
    TRY_MAXIT, TRY_RESETS = 60, 2

    def _eedf_on_grid(self, s, target_mean, warm, growth):
        """The growth decision on one grid (eedf's docstring). Where the
        target lies below the growth term's no-field solution (it did not
        reach it), that solution, converged, is kept: it is the EEDF when
        the solution without the growth term does not reach the target
        either, or does not exist."""
        F0 = nu0 = en0 = off = floor = None
        budget = {}
        if warm is not None:
            F0 = np.interp(self.grid.ec, warm.grid.ec, warm.F, right=0.0)
            nu0 = None if warm.growth_off else warm.nu
            en0 = warm.en_td or None
        if growth and warm is not None and warm.growth_off:
            # the growth term was off: the solution without it first, and
            # the growth term tried from it (on a budget) only where that
            # solution passes growth_valid with its own net frequency (the
            # growth term, where nu < 0, needs less field, so it would fail
            # the test the other did)
            off = self._without_growth(s, target_mean, F0, en0)
            if off is not None:
                ion, att, rec = self._frequencies(off.grid, s, off.F)
                nu_est = ion + s.source_rate() - att - rec
                heating = self.field_power(s, off) + s.source_power()
                if nu_est < 0.0 and \
                        -nu_est * off.mean > self.GROWTH_LIMIT * heating \
                        and -nu_est * self.GROWTH_LIMIT > \
                        self.relax_freq(s, off):
                    return off
                F0, en0 = off.F, off.en_td or None
            budget = dict(maxit=self.TRY_MAXIT, max_resets=self.TRY_RESETS)
        elif growth and warm is None:
            # a cold start: the growth term first, on a budget
            budget = dict(maxit=self.TRY_MAXIT, max_resets=self.TRY_RESETS)
        if growth:
            sol = self.solve_joint(s, target_mean, F0, nu0, en0, True,
                                   **budget)
            if sol is not None and sol.reached and \
                    self.growth_valid(s, sol):
                return sol
            if sol is not None and not sol.reached:
                z = self.solve_joint(s, None, sol.F, sol.nu, en_fixed=0.0)
                if z is not None and self.growth_valid(s, z):
                    z.reached = False
                    floor = z
                    # where the population grows, or the gas relaxes the
                    # shape faster than the turnover, this solution stands:
                    # one without the term that reaches a colder mean gets
                    # there only by putting the lost electrons back cold
                    if z.nu >= 0.0 or (-z.nu * self.GROWTH_LIMIT
                                       <= self.relax_freq(s, z)):
                        return z
        if off is None:
            off = self._without_growth(s, target_mean, F0, en0)
        if off is not None and (off.reached or floor is None):
            return off
        if floor is not None:
            return floor
        raise RuntimeError(
            "kinetic EEDF: no solution with or without the growth term at "
            f"mean {target_mean} eV (n_e {s.ne:.3g} m^-3, Tg {s.tg:.4g} K, "
            f"x {s.x})")

    # ------------------------------------------------------ what LIGHT reads
    def _split(self, s: State, sol: Solution):
        """The solution's flux coefficients and c(z) = (z/2) coth(z/2) at
        each inner boundary, z = W du / D. The Scharfetter-Gummel flux the
        solver uses is the central one, W F_b - D dF/deps, with its
        diffusive part times c(z), so split that way (elastic, field and
        e-e alike) its parts sum to the solver's own flux. Central
        differences alone give the solver's own equilibrium an elastic loss
        of 4.7 % of its cooling term on 600 cells to 12 eV at 300 K, where
        the cell is 0.77 kTg and c(z) = 1.049 (2026-09-30)."""
        g = sol.grid
        parts = self._parts(g, s, sol.nu, sol.F)
        _, W, D0, dfield, hf = parts
        return dfield, _sg_factor(W, self._d_total(parts, sol.en_td), g.du), hf

    def elastic_power(self, s: State, sol: Solution, gas: str) -> float:
        """The elastic energy loss per electron per unit fraction of a gas
        [eV m^3/s], as the solver exchanges it: gamma sum_b eps^2 sigma_eps
        (F + c(z) (kTg/e) dF/deps) du at the inner boundaries (_split). For
        a Maxwellian on a fine grid it is (2m/M) <sigma_m v eps> (1 -
        Tg/Te), LIGHT's own form."""
        g = sol.grid
        _, c, _ = self._split(s, sol)
        fb = 0.5 * (sol.F[1:] + sol.F[:-1])
        dfb = (sol.F[1:] - sol.F[:-1]) / g.du
        kt = KB * s.tg / QE
        return GAMMA * float(np.sum(g.eb ** 2 * g.sig_eps[gas]
                                    * (fb + kt * c * dfb))) * g.du

    def field_power(self, s: State, sol: Solution) -> float:
        """Joule power per electron [eV/s], e E w, as the solver takes it
        from the DC field (_split)."""
        if sol.en_td <= 0.0:
            return 0.0
        g = sol.grid
        dfield, c, _ = self._split(s, sol)
        en_v = sol.en_td * TD
        return -en_v ** 2 * s.n_col * float(np.sum(
            dfield * c * (sol.F[1:] - sol.F[:-1])))

    def em_power(self, s: State, sol: Solution, ions=True) -> float:
        """v0_4: the EM wave's absorbed power per electron [eV/s], the part
        of the solver's own flux its field term carries (_split). With
        ions=False the same with the ions' Coulomb collisions left out of
        the wave's momentum transfer, so the difference is the
        electron-ion (inverse bremsstrahlung) share."""
        if s.omega <= 0.0 or s.en_em_td <= 0.0:
            return 0.0
        _, c, hf = self._split(s, sol)
        d = hf[0] if ions else hf[2]
        en_v = hf[1] * TD
        return -en_v ** 2 * s.n_col * float(np.sum(
            d * c * (sol.F[1:] - sol.F[:-1])))


def _sg_factor(W, D, du):
    """c(z) = (z/2) coth(z/2), z = W du / D, at each inner boundary: the
    Scharfetter-Gummel flux is the central W F_b - D dF/deps with its
    diffusive part times c(z)."""
    h = 0.5 * W * du / D
    small = np.abs(h) < 1e-4
    return np.where(small, 1.0 + h * h / 3.0,
                    h / np.tanh(np.where(small, 1.0, h)))


class _NoSolution(Exception):
    """A linear solve that is not finite or has no positive mean."""


def _decades(F) -> float:
    """How many decades F falls from its peak to the grid's top cell."""
    peak = float(F.max())
    return math.log10(peak / max(float(F[-1]), peak * 1e-300))


def _coulomb_log(ne_m3, mean_ev):
    kte = 2.0 / 3.0 * max(mean_ev, 1e-3) * QE
    lam = 12.0 * math.pi * (EPS0 * kte / QE ** 2) ** 1.5 / math.sqrt(ne_m3)
    return math.log(max(lam, 2.0))


HBAR = 1.054571817e-34
E4PE = QE / (4.0 * math.pi * EPS0)          # e/(4 pi eps0) [V m]


def _sigma_ei_dc(eps_ev, ne_m3, mean_ev):
    """v0_4: the electron-ion momentum-transfer cross section of one singly
    charged ion for a static (DC) field [m^2], per electron energy eps
    [eV]: pi b0^2 L, b0 = e/(4 pi eps0 eps) the distance of closest
    approach, L = (1/2) ln(1 + (b_max/b_min)^2) with b_max the electron
    Debye length at the mean energy and b_min the larger of b0/2 and the
    de Broglie half-wavelength hbar/(2 m v) (the classical and quantum
    cutoffs of the Spitzer logarithm, per energy)."""
    eps = np.maximum(np.asarray(eps_ev, dtype=float), 1e-6)
    kte = 2.0 / 3.0 * max(mean_ev, 1e-3)
    lam_d = math.sqrt(EPS0 * kte / (max(ne_m3, 1.0) * QE))
    v = GAMMA * np.sqrt(eps)
    b0 = E4PE / eps
    bmin = np.maximum(0.5 * b0, HBAR / (2.0 * ME * v))
    return math.pi * b0 ** 2 * 0.5 * np.log1p((lam_d / bmin) ** 2)
