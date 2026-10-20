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
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
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
def processes(gas: str) -> tuple:
    """Parsed processes of a vendored file; refuses a file whose md5 moved."""
    name, md5, _ = FILES[gas]
    path = DATA_DIR / name
    got = md5_of(path)
    if got != md5:
        raise RuntimeError(f"{path} md5 {got} != pinned {md5}: the vendored "
                           "cross sections changed; re-pin deliberately")
    return tuple(parse_lxcat(path))


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
