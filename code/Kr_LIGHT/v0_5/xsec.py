"""Electron-impact rate coefficients k(Te) = <sigma v> over a Maxwellian.

Gary's notes: "Q = k n_e [Kr], where k is <sigma v>". This module is that
sentence: LXCat cross sections integrated over a Maxwellian EEDF at Te,

    k(Te) = int sigma(eps) v(eps) g(eps; Te) d eps,
    g = (2/sqrt(pi)) Te^-3/2 sqrt(eps) exp(-eps/Te),     v = sqrt(2 e eps / m),

tabulated once on a Te grid and read back log-log. Also here, because they
are the same integral with a different weight:

* the elastic energy-loss coefficient (the Davydov/BOLSIG+ elastic term
  for a Maxwellian, where the gas-temperature correction is exactly
  (1 - Tg/Te)):  P_el = n_e N (2m/M) (1 - Tg/Te) <sigma_m v eps>.

The EEDF. Below the 9.9 eV Kr threshold it is the Maxwellian at Te, which
electron-electron collisions enforce at the laser's ionization degree.
Above it, excitation empties the tail faster than those collisions refill
it, so a Maxwellian k(Te) there is an upper bound. The basic network
integrates its Kr rows over a two-temperature EEDF whose tail temperature
comes from the kinetic EEDF (two_t_table and tail_ratio_table below, since
2026-09-28). LIGHT_Kr reads the kinetic table's F2 = 0 column, pure
krypton; the default closure is the kinetic EEDF itself (boltzmann.py).

Data: vendored copies of the repo's LXCat files, pinned by md5 so an edit
anywhere shows up as a failed test instead of a silent change. LIGHT_Kr
vendors krypton's alone (LIGHT's F2 file is not carried). v0_2
(2026-10-14) reads Biagi's Magboltz 8.97 set, LXCat's 'Biagi' database:
51 excitation levels, each its own block at its own threshold, where
v0_1's Biagi-v7.1 grouped them into four (5s, 5p, d+p, 5d) at the
groups' lowest thresholds. LXCat's note on v7.1 points to this set for
the levels. Every level is a row of v0_2's network (network.py).

v0_3 (2026-10-14) trains a few factors on that file against krypton's
measured swarm data (ML/ml/kr_swarm.py, the owner's "use ML to train
it"): the elastic momentum transfer times a factor read log-linearly in
energy between knots, and one factor on each group of excitation levels
by its Paschen label (1s, 2p, the rest). processes() is the file with the
trained factors (Factors, TRAINED_FILE); raw_processes() is the file as
Biagi published it, and using() runs with other factors (ML's fits; UNIT
gives v0_2's cross sections bit for bit).
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from contextlib import contextmanager
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent / "data"

# gas -> (vendored file, md5, where it was copied from)
FILES = {
    "Kr": ("Kr_Biagi897_LXCat.txt", "faeab9cbbca5a4e1789374b8048414c4",
           "old_files/Kr_LXCat_new.txt (LXCat's 'Biagi' database for Kr, "
           "transcribed from Magboltz 8.97, Sept 2011; retrieved "
           "2026-06-18 by the owner, Kinetic_Modeling commit 36f70d3, "
           "moved unchanged to old_files/ on 2026-08-25; CRLF line ends as "
           "LXCat serves them)"),
}

QE = 1.602176634e-19          # C  (also J per eV)
ME = 9.1093837015e-31         # kg
EPS0 = 8.8541878128e-12       # F/m
AMU = 1.66053906660e-27       # kg
MASS_AMU = {"Kr": 83.798}
V_CM_S = math.sqrt(2.0 * QE / ME) * 100.0     # v [cm/s] = V_CM_S sqrt(eps[eV])

EPS = np.geomspace(1e-5, 1.0e3, 8000)         # quadrature grid [eV]
TE = np.geomspace(0.01, 30.0, 600)            # table grid [eV]
TE_MIN, TE_MAX = float(TE[0]), float(TE[-1])
_K_FLOOR = 1e-300


def _trap_weights(x):
    """Trapezoid weights: int f dx ~= sum(f * w) on the (non-uniform) grid."""
    w = np.zeros_like(x)
    dx = np.diff(x)
    w[:-1] += 0.5 * dx
    w[1:] += 0.5 * dx
    return w


_W = _trap_weights(EPS)


def maxwellian(eps, te):
    """g(eps; Te) [1/eV], normalized to 1 over 0..inf."""
    return (2.0 / math.sqrt(math.pi)) * te ** -1.5 * np.sqrt(eps) \
        * np.exp(-eps / te)


# ---------------------------------------------------------------------------
# LXCat parsing
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Process:
    kind: str            # elastic|effective|excitation|ionization|attachment
    target: str
    product: str         # the LXCat product label, e.g. 'Kr(1S5)'
    threshold: float     # eV (0 for elastic/attachment)
    energy: np.ndarray   # eV
    sigma: np.ndarray    # cm^2
    mass_ratio: float = 0.0   # m/M as printed on elastic/effective blocks

    def sigma_at(self, eps):
        # zero outside the table (tables end near 1 keV; the Maxwellian at
        # Te <= 30 eV carries < e^-30 of its weight up there)
        return np.interp(eps, self.energy, self.sigma, left=0.0, right=0.0)


_KINDS = ("ELASTIC", "EFFECTIVE", "EXCITATION", "IONIZATION", "ATTACHMENT")


def _floats(line):
    try:
        return [float(p) for p in line.split()]
    except ValueError:
        return None


def parse_lxcat(path) -> list:
    lines = Path(path).read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        kw = lines[i].strip()
        if kw not in _KINDS:
            i += 1
            continue
        target, _, product = lines[i + 1].replace("<->", "->").partition("->")
        j, thr, mr = i + 2, 0.0, 0.0
        vals = _floats(lines[j])
        if kw != "ATTACHMENT" and vals:
            if kw in ("EXCITATION", "IONIZATION"):
                thr = vals[0]
            else:
                mr = vals[0]
            j += 1
        while not lines[j].lstrip().startswith("---"):
            j += 1
        j += 1
        e, s = [], []
        while not lines[j].lstrip().startswith("---"):
            a, b = lines[j].split()[:2]
            e.append(float(a))
            s.append(float(b) * 1e4)          # m^2 -> cm^2
            j += 1
        out.append(Process(kw.lower(), target.strip(), product.strip(), thr,
                           np.asarray(e), np.asarray(s), mr))
        i = j + 1
    return out


def md5_of(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def raw_processes(gas: str) -> tuple:
    """Parsed processes of a vendored file, as published; refuses a file
    whose md5 moved."""
    name, md5, _ = FILES[gas]
    path = DATA_DIR / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}: the vendored "
                           "cross sections changed; re-pin deliberately")
    return tuple(parse_lxcat(path))


# ---------------------------------------------------------------------------
# v0_3: the swarm-trained factors on Biagi's cross sections
# ---------------------------------------------------------------------------
# ML's swarm analysis (ML/ml/kr_swarm.py; `python3 ML/run.py kr-train`
# writes the file) fits these factors to the krypton swarm datasets its
# consistency check admits (drift velocity, characteristic energy and the
# ionization coefficient). The ionization cross section is not trained
# (its factor stays 1). The vendored file stays Biagi's, md5-pinned, and
# the factors are pinned the same way.
TRAINED_FILE = ("kr_swarm_trained.json", "d7c48cd865d79c81b9b449ff51c00c69")
EXC_GROUPS = ("1s", "2p", "hi")


@dataclass(frozen=True)
class Factors:
    """Multipliers on Biagi's cross sections: m at m_knots [eV] on the
    elastic momentum transfer (log-linear in energy between the knots,
    the end knots' values beyond them), exc on the 1s (5s), 2p (5p) and
    remaining excitation levels, and ion on the ionization."""
    m_knots: tuple = (1.0,)
    m: tuple = (1.0,)
    exc: tuple = (1.0, 1.0, 1.0)
    ion: float = 1.0

    def m_at(self, e) -> np.ndarray:
        e = np.asarray(e, dtype=float)
        if len(self.m_knots) == 1:
            return np.full(e.shape, float(self.m[0]))
        le = np.log(np.maximum(e, 1e-300))
        return np.exp(np.interp(le, np.log(self.m_knots), np.log(self.m)))


UNIT = Factors()


def exc_group(product: str) -> str:
    """The training group of an excitation level from its LXCat label:
    'Kr(1S5)' -> '1s', 'Kr(2P10)' -> '2p', any other level -> 'hi'."""
    tag = product.split("(", 1)[1].upper() if "(" in product else ""
    if tag.startswith("1S"):
        return "1s"
    if tag.startswith("2P"):
        return "2p"
    return "hi"


def apply_factors(procs, f: Factors) -> tuple:
    """The processes with their factors; a unit factor leaves a process
    untouched (so UNIT returns the file's own arrays)."""
    out = []
    for p in procs:
        if p.kind == "elastic":
            fac = f.m_at(p.energy)
            out.append(p if np.all(fac == 1.0)
                       else replace(p, sigma=p.sigma * fac))
            continue
        if p.kind == "excitation":
            x = f.exc[EXC_GROUPS.index(exc_group(p.product))]
        elif p.kind == "ionization":
            x = f.ion
        else:
            x = 1.0
        out.append(p if x == 1.0 else replace(p, sigma=p.sigma * x))
    return tuple(out)


def read_factors(path) -> Factors:
    """A factors file (JSON): m_knots_eV, m, exc {1s, 2p, hi} and ion,
    with whatever record of the fit the trainer adds."""
    d = json.loads(Path(path).read_text())
    return Factors(tuple(float(v) for v in d["m_knots_eV"]),
                   tuple(float(v) for v in d["m"]),
                   tuple(float(d["exc"][g]) for g in EXC_GROUPS),
                   float(d["ion"]))


@lru_cache(maxsize=1)
def trained() -> Factors:
    """The trained factors, from the pinned file (UNIT until it exists)."""
    name, md5 = TRAINED_FILE
    if md5 is None:
        return UNIT
    path = DATA_DIR / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}: the trained "
                           "factors changed; retrain with ML/run.py "
                           "kr-train and re-pin deliberately")
    return read_factors(path)


_ACTIVE: list = []


def factors() -> Factors:
    """The factors in force: using()'s, or the trained ones."""
    return _ACTIVE[-1] if _ACTIVE else trained()


@contextmanager
def using(f: Factors):
    """Run with other factors (ML's fits; UNIT is v0_2's cross sections).
    Every table built on the cross sections is rebuilt inside and after;
    a model compiled inside keeps the factors it was compiled with."""
    _ACTIVE.append(f)
    _clear_tables()
    try:
        yield f
    finally:
        _ACTIVE.pop()
        _clear_tables()


def _clear_tables():
    for fn in (processes, lumped_5s, rate_table, superelastic_table,
               elastic_loss_table, capture_energy_table, e_break,
               _tail_sums, two_t_table):
        fn.cache_clear()
    _COMMON[0] = _COMMON[1] = None
    ib = sys.modules.get(f"{__package__}.ib")
    if ib is not None:
        ib._MT_CACHE.clear()


@lru_cache(maxsize=None)
def processes(gas: str) -> tuple:
    """The processes the model and its solver use: the vendored file with
    the factors in force (v0_3: the trained ones)."""
    raw = raw_processes(gas)
    return apply_factors(raw, factors()) if gas == "Kr" else raw


# v0_2: the four 5s levels as one cross section at the manifold's energy
# (the 1s5 threshold), for the superelastic return of the lumped Kr*
# (network.py SUPER), as v0_1 read Biagi-v7.1's 5s group. find() returns
# it under this label; it is not one of processes(), so the kinetic
# equation's momentum transfer does not count the 5s levels twice
LUMPED_5S = "Kr(5s manifold)"


@lru_cache(maxsize=1)
def lumped_5s() -> Process:
    lv = [p for p in processes("Kr") if p.kind == "excitation"]
    eb = min(p.threshold for p in lv)
    five_s = [p for p in lv if p.threshold < eb + 1.0]
    if len(five_s) != 4:
        raise RuntimeError(f"{len(five_s)} levels within 1 eV of the "
                           "lowest: the 5s manifold has four")
    e = np.unique(np.concatenate([p.energy for p in five_s]))
    s = sum(p.sigma_at(e) for p in five_s)
    return Process("excitation", "Kr", LUMPED_5S, eb, e, s)


def find(gas: str, kind: str, product: str = None) -> Process:
    if gas == "Kr" and kind == "excitation" and product == LUMPED_5S:
        return lumped_5s()
    hits = [p for p in processes(gas) if p.kind == kind
            and (product is None or p.product == product)]
    if len(hits) != 1:
        raise KeyError(f"{gas} {kind} {product!r}: {len(hits)} matches")
    return hits[0]


# ---------------------------------------------------------------------------
# Maxwellian tables
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _weights():
    """(n_te, n_eps) matrix of g(eps; Te) * trapezoid weight, row-normalized
    so the discrete Maxwellian integrates to exactly 1 (cancels quadrature
    error in the normalization)."""
    g = maxwellian(EPS[None, :], TE[:, None]) * _W[None, :]
    return g / g.sum(axis=1, keepdims=True)


def _mean(f):
    """<f>(Te) on the TE grid. Multiply-and-sum, not matmul: numpy 2 on
    macOS Accelerate raises spurious divide/overflow warnings in a matmul
    whose operands hold underflowed zeros (seen 2026-09-25)."""
    return (_weights() * f[None, :]).sum(axis=1)


class RateTable:
    """k(Te) on the TE grid, read back linearly in (log Te, log k)."""

    def __init__(self, values, label):
        self.values = np.maximum(np.asarray(values, dtype=float), _K_FLOOR)
        self._logk = np.log(self.values)
        self.label = label

    def __call__(self, te):
        lt = np.log(np.clip(te, TE_MIN, TE_MAX))
        return np.exp(np.interp(lt, np.log(TE), self._logk))


@lru_cache(maxsize=None)
def rate_table(gas: str, kind: str, product: str = None) -> RateTable:
    """k(Te) = <sigma v> [cm^3/s] for one LXCat process."""
    p = find(gas, kind, product)
    k = _mean(p.sigma_at(EPS) * V_CM_S * np.sqrt(EPS))
    return RateTable(k, f"{gas} {kind} {p.product}")


@lru_cache(maxsize=None)
def superelastic_table(gas: str, product: str, g_ratio: float) -> RateTable:
    """De-excitation k(Te) [cm^3/s] from the Klein-Rosseland cross section
    sigma_sup(eps) = g_ratio (eps + E)/eps sigma_exc(eps + E), integrated
    directly. (Detailed balance k_exc (g_ratio) exp(E/Te) is the same number
    for a Maxwellian, but evaluating it that way multiplies an underflowed
    table by exp(+381) at room-temperature Te.)"""
    p = find(gas, "excitation", product)
    e_thr = p.threshold
    sig = g_ratio * (EPS + e_thr) / EPS * p.sigma_at(EPS + e_thr)
    k = _mean(sig * V_CM_S * np.sqrt(EPS))
    return RateTable(k, f"{gas} superelastic from {p.product}")


@lru_cache(maxsize=None)
def elastic_loss_table(gas: str) -> RateTable:
    """(2m/M) <sigma_m v eps> [eV cm^3/s]; times n_e N (1 - Tg/Te) gives the
    elastic power [eV cm^-3 s^-1] (exact for a Maxwellian EEDF)."""
    p = find(gas, "elastic")
    k = _mean(p.sigma_at(EPS) * V_CM_S * np.sqrt(EPS) * EPS)
    return RateTable(two_m_over_M(gas) * k, f"{gas} elastic energy loss")


def two_m_over_M(gas: str) -> float:
    """Fractional energy loss per elastic collision, 2 m_e / M, from the
    CODATA electron mass and the atomic mass (the Biagi file prints the
    rounded m/M = 6.5e-6 for Kr; the exact value is 6.5465e-6)."""
    return 2.0 * ME / (MASS_AMU[gas] * AMU)


# ---------------------------------------------------------------------------
# The captured electron's energy (2026-09-28, step 3 of the owner's plan)
# ---------------------------------------------------------------------------
class EnergyTable:
    """E(Te) [eV] on the TE grid, read back linearly in log Te."""

    def __init__(self, values, label):
        self.values = np.asarray(values, dtype=float)
        self.label = label

    def __call__(self, te):
        lt = np.log(np.clip(te, TE_MIN, TE_MAX))
        return float(np.interp(lt, np.log(TE), self.values))


@lru_cache(maxsize=None)
def capture_energy_table(gas: str, kind: str, product: str = None):
    """<eps sigma v> / <sigma v> [eV]: the mean energy of the electrons a
    process removes from a Maxwellian at Te. Differentiating k(Te) under
    the integral gives the same number as (3/2 + d ln k / d ln Te) Te,
    which is how model.py books the analytic rate laws. (In LIGHT's Kr/F2,
    F2 attachment removed slow electrons: 0.54 eV at Te = 1.13 eV, not
    (3/2) Te. No pure-Kr row captures an electron by a cross section.)"""
    p = find(gas, kind, product)
    sv = p.sigma_at(EPS) * V_CM_S * np.sqrt(EPS)
    k = _mean(sv)
    ke = _mean(sv * EPS)
    ok = k > _K_FLOOR
    e = np.where(ok, ke / np.where(ok, k, 1.0), 1.5 * TE)
    return EnergyTable(e, f"{gas} {kind} {p.product}: captured energy")


# ---------------------------------------------------------------------------
# The two-temperature EEDF (2026-09-28, step 3 of the owner's plan)
# ---------------------------------------------------------------------------
# Above the lowest Kr threshold the laser plasma's EEDF is not Maxwellian:
# excitation empties the tail faster than electron-electron collisions
# refill it (ML/ml/eedf.py, the kinetic beam EEDF). LIGHT keeps the
# Maxwellian at Te below the break and lets it fall with a tail temperature
# T_t above it, continuous at the break:
#     g(eps) ~ sqrt(eps) exp(-eps/Te)                          eps < E_b
#     g(eps) ~ sqrt(eps) exp(-E_b/Te) exp(-(eps - E_b)/T_t)    eps >= E_b
# For a process whose threshold is at or above E_b the rate separates into
# 1-D sums on the quadrature grid,
#     k = exp(-E_b/Te) I(T_t) / (Z_b(Te) + exp(-E_b/Te) J(T_t)),
#     I = sum_{eps>=E_b} sigma v sqrt(eps) exp(-(eps-E_b)/T_t) w,
#     J = the same sum without sigma v,  Z_b = sum_{eps<E_b} sqrt(eps)
#     exp(-eps/Te) w,
# so T_t = Te is the Maxwellian table exactly. T_t = r Te, with r read from
# the kinetic table below (tail_ratio_table) at the run's n_e/N, on its
# F2 = 0 column for pure krypton (model.py).
TAIL_GAS = "Kr"
# the kinetic tail ratio: generated by ML (ML/ml/tail.py), pinned by md5 as
# the cross sections are (edit the generator and regenerate, never the file)
TAIL_FILE = ("kr_tail_ratio.txt", "46616e447404c155b3f80ee23e67db8a")
# v0_4: the field-driven pure-Kr table, r(n_e/N, Te), from this model's own
# kinetic EEDF (tail_field.py), the closure's default (Params.tail =
# 'field'); v0_1-v0_3's beam table above stays as Params.tail = 'beam'.
# v0_5: r(n_e/N, Te, x_s), x_s the 5s fraction, its x_s = 0 block v0_4's
# table (Params.tail = 'field0' reads that block alone)
TAIL_FIELD_FILE = ("kr_tail_field.txt", "c04a5d06ba496de23d132f06e4a97525")


@lru_cache(maxsize=1)
def e_break() -> float:
    """The lowest Kr excitation threshold (Kr 1s5, 9.9152 eV): the break."""
    return min(p.threshold for p in processes(TAIL_GAS)
               if p.kind == "excitation")


def in_tail(gas: str, kind: str, product: str = None) -> bool:
    """True for the processes the two-temperature EEDF changes: Kr
    excitation and ionization (thresholds at or above the break)."""
    if gas != TAIL_GAS or kind not in ("excitation", "ionization"):
        return False
    return find(gas, kind, product).threshold >= e_break() - 1e-9


@lru_cache(maxsize=1)
def _tail_sums():
    """(tail mask, log Z_b(Te), log J(T_t)) on the TE grid (T_t uses the
    same grid)."""
    eb = e_break()
    tail = EPS >= eb
    w = _W * np.sqrt(EPS)
    zb = (w[None, ~tail] * np.exp(-EPS[None, ~tail] / TE[:, None])).sum(axis=1)
    j = (w[None, tail] * np.exp(-(EPS[None, tail] - eb) / TE[:, None])).sum(
        axis=1)
    return tail, np.log(zb), np.log(j)


class TwoTTable:
    """k(Te, T_t) [cm^3/s] of one tail process over the two-temperature
    EEDF, from its log I(T_t) on the TE grid."""

    def __init__(self, log_i, label):
        self.log_i = np.asarray(log_i, dtype=float)
        self.label = label

    def __call__(self, te, tt):
        x, ltt, lz = _two_t_common(te, tt)
        li = float(np.interp(ltt, _LOG_TE, self.log_i))
        return max(math.exp(x + li - lz), _K_FLOOR)


_LOG_TE = np.log(TE)
_COMMON = [None, None]


def _two_t_common(te, tt):
    """The parts of a two-temperature rate every tail process shares at
    (Te, T_t): -E_b/Te, ln T_t and ln(Z_b + exp(-E_b/Te) J). v0_2 reads 72
    tail rows at one state, so the last state's are kept (the same
    numbers, computed once)."""
    if _COMMON[0] == (te, tt):
        return _COMMON[1]
    _, log_zb, log_j = _tail_sums()
    lte = math.log(min(max(te, TE_MIN), TE_MAX))
    ltt = math.log(min(max(tt, TE_MIN), TE_MAX))
    x = -e_break() / math.exp(lte)
    lzb = float(np.interp(lte, _LOG_TE, log_zb))
    lj = float(np.interp(ltt, _LOG_TE, log_j))
    # log(Z_b + exp(x) J), without overflow
    a, b = lzb, x + lj
    lz = max(a, b) + math.log1p(math.exp(-abs(a - b)))
    _COMMON[0], _COMMON[1] = (te, tt), (x, ltt, lz)
    return x, ltt, lz


def two_t_from_sigma(sigma, label) -> TwoTTable:
    """The two-temperature table of a cross section sigma(EPS) [cm^2] that
    vanishes below the break."""
    tail, _, _ = _tail_sums()
    if np.any(sigma[~tail] != 0.0):
        raise ValueError(f"{label}: cross section below the break")
    eb = e_break()
    sv = (sigma * V_CM_S * np.sqrt(EPS) * _W * np.sqrt(EPS))[tail]
    i = (sv[None, :] * np.exp(-(EPS[None, tail] - eb) / TE[:, None])).sum(
        axis=1)
    return TwoTTable(np.log(np.maximum(i, _K_FLOOR)), label)


@lru_cache(maxsize=None)
def two_t_table(gas: str, kind: str, product: str = None) -> TwoTTable:
    if not in_tail(gas, kind, product):
        raise ValueError(f"{gas} {kind} {product}: threshold below the break")
    p = find(gas, kind, product)
    return two_t_from_sigma(p.sigma_at(EPS),
                            f"{gas} {kind} {p.product} (two-temperature)")


class TailRatio:
    """r = T_t / Te from the kinetic table, bilinear in (log10 n_e/N, F2
    fraction) and held at the table's edges outside it. The axes are the
    kinetic states' own: n_e/N sets how fast electron-electron collisions
    refill the tail, and F2 how fast the injected electrons are slowed
    below the break (ML/ml/tail.py)."""

    def __init__(self, x, f2, r, header):
        self.x = np.asarray(x, dtype=float)       # log10(n_e/N)
        self.f2 = np.asarray(f2, dtype=float)     # [F2]/N
        self.r = np.asarray(r, dtype=float)       # shape (len x, len f2)
        self.header = header

    def __call__(self, ne_over_n, f2_over_n):
        x = math.log10(max(ne_over_n, 1e-300))
        x = min(max(x, self.x[0]), self.x[-1])
        f = min(max(f2_over_n, self.f2[0]), self.f2[-1])
        i = min(max(int(np.searchsorted(self.x, x)) - 1, 0), len(self.x) - 2)
        j = min(max(int(np.searchsorted(self.f2, f)) - 1, 0),
                len(self.f2) - 2)
        a = (x - self.x[i]) / (self.x[i + 1] - self.x[i])
        b = (f - self.f2[j]) / (self.f2[j + 1] - self.f2[j])
        r = self.r
        return float((1 - a) * (1 - b) * r[i, j] + a * (1 - b) * r[i + 1, j]
                     + (1 - a) * b * r[i, j + 1] + a * b * r[i + 1, j + 1])


def read_tail_ratio(path) -> TailRatio:
    """The table file: '#' header lines, then a line 'F2' followed by the F2
    fractions of the columns, then one line per log10(n_e/N): the value
    and its row of r."""
    header, rows, f2 = [], [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            header.append(line[1:].strip())
        elif line.split() and line.split()[0] == "F2":
            f2 = [float(v) for v in line.split()[1:]]
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    rows = np.asarray(rows)
    return TailRatio(rows[:, 0], f2, rows[:, 1:], "\n".join(header))


class FieldTailRatio:
    """v0_4: r = T_t / Te of field-driven pure Kr (tail_field.py),
    bilinear in (log10 n_e/N, ln Te) on ln r (r spans 1e-3 to 1) and held
    at the table's edges. v0_5: a third axis, x_s = [Kr(5s)]/[Kr], read
    linearly in log10(x_s + XS0); a table without it (v0_4's) reads as
    one block."""

    def __init__(self, x, te, r, header, xs=None, xs0=1e-7):
        self.x = np.asarray(x, dtype=float)       # log10(n_e/N)
        self.te = np.asarray(te, dtype=float)     # eV
        self.lte = np.log(self.te)
        r = np.asarray(r, dtype=float)
        if r.ndim == 2:
            r = r[:, :, None]
            xs = (0.0,) if xs is None else xs
        self.r = r                                # (len x, len te, len xs)
        self.lr = np.log(self.r)
        self.xs = np.asarray(xs, dtype=float)
        self.xs0 = xs0
        self.zs = np.log10(self.xs + xs0)
        self.header = header

    def __call__(self, ne_over_n, te, x_s=0.0):
        x = math.log10(max(ne_over_n, 1e-300))
        x = min(max(x, self.x[0]), self.x[-1])
        lt = math.log(max(te, 1e-300))
        lt = min(max(lt, self.lte[0]), self.lte[-1])
        i = min(max(int(np.searchsorted(self.x, x)) - 1, 0), len(self.x) - 2)
        j = min(max(int(np.searchsorted(self.lte, lt)) - 1, 0),
                len(self.lte) - 2)
        a = (x - self.x[i]) / (self.x[i + 1] - self.x[i])
        b = (lt - self.lte[j]) / (self.lte[j + 1] - self.lte[j])
        if len(self.zs) == 1:
            ks, cs = (0, 0), (1.0, 0.0)
        else:
            z = math.log10(max(x_s, 0.0) + self.xs0)
            z = min(max(z, self.zs[0]), self.zs[-1])
            k = min(max(int(np.searchsorted(self.zs, z)) - 1, 0),
                    len(self.zs) - 2)
            c = (z - self.zs[k]) / (self.zs[k + 1] - self.zs[k])
            ks, cs = (k, k + 1), (1.0 - c, c)
        lr = self.lr
        out = 0.0
        for kk, ck in zip(ks, cs):
            if ck == 0.0:
                continue
            out += ck * ((1 - a) * (1 - b) * lr[i, j, kk]
                         + a * (1 - b) * lr[i + 1, j, kk]
                         + (1 - a) * b * lr[i, j + 1, kk]
                         + a * b * lr[i + 1, j + 1, kk])
        return math.exp(out)


def read_tail_field(path) -> FieldTailRatio:
    """The field table's file: '#' header lines; v0_5: blocks, each a line
    'XS' and its x_s, a line 'TE' and the columns' Te [eV], then one line
    per log10(n_e/N): the value and its row of r (v0_4's file: one block
    without the 'XS' line)."""
    header, blocks, te, xs = [], [], None, []
    cur = None
    for line in Path(path).read_text().splitlines():
        w = line.split()
        if line.startswith("#"):
            header.append(line[1:].strip())
        elif w and w[0] == "XS":
            xs.append(float(w[1]))
            cur = []
            blocks.append(cur)
        elif w and w[0] == "TE":
            te = [float(v) for v in w[1:]]
            if cur is None:
                cur = []
                blocks.append(cur)
        elif w:
            cur.append([float(v) for v in w])
    arrs = [np.asarray(b) for b in blocks]
    x = arrs[0][:, 0]
    r = np.stack([b[:, 1:] for b in arrs], axis=2)
    return FieldTailRatio(x, te, r, "\n".join(header),
                          xs=tuple(xs) if xs else (0.0,))


@lru_cache(maxsize=1)
def tail_field_table() -> FieldTailRatio:
    name, md5 = TAIL_FIELD_FILE
    path = DATA_DIR / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}: the field "
                           "tail table changed; regenerate it with "
                           "'python3 -m v0_5.tail_field' and re-pin "
                           "deliberately")
    return read_tail_field(path)


@lru_cache(maxsize=1)
def tail_ratio_table() -> TailRatio:
    name, md5 = TAIL_FILE
    path = DATA_DIR / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}: the tail "
                           "table changed; regenerate it with ML/ml/tail.py "
                           "and re-pin deliberately")
    return read_tail_ratio(path)
