"""LIGHT's KrF(B,C) constants against the decay measurements (evidence
script, not LIGHT code).

Two kinds of data, both on pure Kr with F2, no electrons:
  - Eden, Waynant, Searles & Burnham, APL 32, 733 (1978). KrF2 photolysis
    at 193 nm makes vibrationally cold KrF(B) alone. The B-X emission decay
    constant Sigma is measured against Kr (Fig. 3) and F2 (Fig. 5). The
    pump (about 3 ns FWHM with a 2.5-ns exponential tail) was deconvolved
    (JAP 49, 5368).
  - Gadomski, Xu, Setser & Nelson, CPL 189, 153 (1992), Figs. 1 and 2
    (Kansas State). A two-photon Kr + F2 reaction makes hot B and C. The
    slow decay constant is taken from the latest part of the waveform
    against Kr (1 and 2 Torr F2) and F2 (50 and 100 Torr Kr).
The points are figure readings in digitized/ (headers give the
uncertainties).

Eden's Sigma is compared with the intrinsic B(t) of the two-level model
after a B-only start, fitted with one exponential from t = 0 down to 5 %
(what the deconvolution recovers). The Kansas State points are compared
with the slow eigenvalue lam-.

Also Hakuta, Fujino, Nakayama & Takuma, JJAP 29, 1194 (1990), Figs. 2 and 3:
the slow C-A decay after 248-nm photoassociation, against F2 (9 Torr Kr)
and Kr (1 Torr F2), compared with lam-.

Part 1 scores LIGHT's constants now (the joint fit adopted 2026-09-28),
LIGHT's until 2026-09-27, and Gadomski's published set.
Part 2 frees every two-level constant (B and C quench separately, tau_B,
tau_C, K_eq) and profiles the Eden + Kansas State fit over the B -> C
mixing rate k_BC, with the share of KrF in the lasing B state at 1 atm.
Part 3 redoes Gadomski's reduction (K_eq and k_Kr fitted to the Kansas
State decays) with his assumed inputs and with Hakuta's measured ones.
It was the basis of LIGHT's constants until 2026-09-27. Nothing here
changes LIGHT.

    /usr/bin/python3 LIGHT/docs/krf_data/bc_check.py      # about 2 s
"""
import csv
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_LIGHT))
import numpy as np
from scipy.optimize import least_squares

from v0 import two_level as L

T = 3.22e16               # cm^-3 per Torr, the figures' own density axes
DATA = _HERE / "digitized"


def _rows(name):
    with open(DATA / name) as fh:
        rows = list(csv.reader(l for l in fh if not l.startswith("#")))
    return rows[1:]


EDEN_KR = [(float(p), float(s) * 1e7) for p, s in _rows("eden1978_APL_fig3_Kr.csv")]
EDEN_F2 = [(float(p), float(s) * 1e7) for p, s in _rows("eden1978_APL_fig5_F2.csv")]
KSU_2 = [(float(f), float(p), float(s) * 1e7)
         for f, p, s in _rows("gadomski1992_fig2.csv")]
KSU_1 = [(100.0 if s_ == "Kr100_F2" else 50.0, float(p), float(v) * 1e7)
         for s_, p, v in _rows("gadomski1992_fig1.csv") if s_ != "Xe"]
HAK_3 = [(float(p), float(v) * 1e7) for p, v in _rows("hakuta1990_fig3_Kr.csv")]
HAK_2 = [(float(p), float(v) * 1e7) for p, v in _rows("hakuta1990_fig2_F2.csv")]
KRF2_LOSS = 1.4e-10 * 0.5 * T      # Eden's KrF2 quench at his ~0.5 Torr
t = np.arange(0.0, 200e-9, 0.05e-9)


def eden_sigma(c, kr, f2):
    B, _ = L.impulse(c, t, B0=1.0, C0=0.0, kr=kr, f2=f2, extra=KRF2_LOSS)
    sel = B > 0.05
    return -np.polyfit(t[sel], np.log(B[sel]), 1)[0]


def residuals(c):
    r = [math.log(eden_sigma(c, p * T, 0.0) / s) for p, s in EDEN_KR]
    r += [math.log(eden_sigma(c, 0.0, p * T) / s) for p, s in EDEN_F2]
    r += [math.log(L.eigen(c, kr=p * T, f2=f * T)[1] / s) for f, p, s in KSU_2]
    r += [math.log(L.eigen(c, kr=kr * T, f2=p * T)[1] / s)
          for kr, p, s in KSU_1]
    return np.array(r)


def hakuta_residuals(c):
    """Hakuta et al. 1990 Figs. 3 (Kr, 1 Torr F2) and 2 (F2, 9 Torr Kr): the
    slow C-A decay rate, compared with the slow eigenvalue lam-."""
    r = [math.log(L.eigen(c, kr=p * T, f2=1.0 * T)[1] / v) for p, v in HAK_3]
    r += [math.log(L.eigen(c, kr=9.0 * T, f2=p * T)[1] / v) for p, v in HAK_2]
    return np.array(r)


def ksu_residuals(c):
    r = [math.log(L.eigen(c, kr=p * T, f2=f * T)[1] / s) for f, p, s in KSU_2]
    r += [math.log(L.eigen(c, kr=kr * T, f2=p * T)[1] / s)
          for kr, p, s in KSU_1]
    return np.array(r)


def rms(r):
    return 100.0 * math.sqrt(np.mean(np.square(np.expm1(np.asarray(r)))))


def redo_gadomski(fixed):
    """Gadomski's reduction of the Kansas State decays: K_eq and k_Kr free,
    everything else as given. Returns the constants and 1-sigma factors."""
    def f(x):
        return ksu_residuals(fixed.with_(K_eq=math.exp(x[0]),
                                         k_Kr=math.exp(x[1])))
    sol = least_squares(f, [math.log(4.0), math.log(3e-12)], diff_step=1e-3)
    s2 = np.sum(sol.fun ** 2) / (len(sol.fun) - 2)
    err = np.sqrt(np.diag(np.linalg.inv(sol.jac.T @ sol.jac) * s2))
    c = fixed.with_(K_eq=math.exp(sol.x[0]), k_Kr=math.exp(sol.x[1]))
    return c, math.exp(err[0]), math.exp(err[1])


def rms_parts(r):
    n1 = len(EDEN_KR)
    n2 = n1 + len(EDEN_F2)
    n3 = n2 + len(KSU_2)

    def f(a):
        return 100.0 * math.sqrt(np.mean(np.square(np.expm1(a))))
    return f(r[:n1]), f(r[n1:n2]), f(r[n2:n3]), f(r[n3:])


NAMES = ("k_BC", "K_eq", "tau_C", "k_Kr", "k_Kr_C", "k3_Kr", "k3_Kr_C",
         "k_F2", "k_F2_C", "tau_B")
LO = np.log([1e-13, 0.1, 20e-9, 1e-14, 1e-14, 1e-33, 1e-33, 5e-11, 5e-11,
             6.0e-9])
HI = np.log([1e-9, 30.0, 150e-9, 5e-11, 5e-11, 2e-30, 2e-30, 2e-9, 2e-9,
             7.6e-9])


def make(x):
    v = dict(zip(NAMES, np.exp(x)))
    return L.Consts(k_e=0.0, **v)


def b_share_1atm(c):
    """[B]/([B]+[C]) at 1 atm pure Kr with 0.5 % F2 and 300 K, for KrF
    formed in B (the steady state of the B/C block)."""
    n = 760.0 * T
    return L.steady_fraction_B(c, kr=0.995 * n, f2=0.005 * n)


def main():
    c_light = L.Consts.from_light()
    # Gadomski 1992's published set (one quench set for B and C), and
    # LIGHT's until 2026-09-27: the same with Hakuta's tau_C and k_BC
    c_gad = c_light.with_(tau_C=75e-9, k_BC=1.0e-10, K_eq=5.0, k_Kr=4.0e-12,
                          k3_Kr=6.0e-31, k_F2=3.0e-10, k_e=4.5e-8,
                          k_Kr_C=None, k3_Kr_C=None, k_F2_C=None, k_e_C=None)
    c_0927 = c_gad.with_(tau_C=60e-9, k_BC=1.2e-10)
    print("Part 1. rms deviation, %:     Eden Kr  Eden F2  KSU Fig2  KSU Fig1"
          "  Hakuta   B share, 1 atm")
    for label, c in (("LIGHT now (joint fit)", c_light),
                     ("LIGHT until 2026-09-27", c_0927),
                     ("Gadomski 1992", c_gad)):
        print(f"{label:30s}" + "%6.1f   %6.1f   %6.1f    %6.1f"
              % rms_parts(residuals(c)) + f"   {rms(hakuta_residuals(c)):5.1f}"
              + f"      {b_share_1atm(c):.2f}")
    print("\nPart 3 (the basis of LIGHT's constants until 2026-09-27). "
          "Gadomski's reduction of the Kansas\nState decays redone, K_eq and "
          "k_Kr free, k_F2 3.0e-10 and k3 6.0e-31 kept:")
    for label, fx in (("his inputs: tau_C 75 ns, k_BC 1.0e-10", c_gad),
                      ("Hakuta's: tau_C 60 ns, k_BC 1.2e-10", c_0927)):
        c, e_k, e_kr = redo_gadomski(fx)
        print(f"  {label:40s} K_eq {c.K_eq:.2f} (x/ {e_k:.2f})  k_Kr "
              f"{c.k_Kr:.2e} (x/ {e_kr:.2f})  KSU {rms(ksu_residuals(c)):.1f} %"
              f"  Hakuta {rms(hakuta_residuals(c)):.1f} %")
    print()
    x0 = np.log([1e-11, 2.0, 55e-9, 2e-12, 2e-12, 1e-30, 4e-31, 5e-10, 3e-10,
                 7e-9])
    best = least_squares(lambda x: residuals(make(x)), x0, bounds=(LO, HI),
                         diff_step=1e-3, max_nfev=600)
    print("Part 2. All two-level constants free (B and C quench separate),"
          " Eden and KSU only, profiled over k_BC:")
    print("  k_BC      K_eq  tau_C   Eden Kr  Eden F2  KSU Fig2  KSU Fig1"
          "   B share, 1 atm")
    for kbc in (3e-12, 1e-11, 2e-11, 3e-11, 5e-11, 1e-10, 2e-10):
        lo, hi = LO.copy(), HI.copy()
        lo[0], hi[0] = math.log(kbc) - 1e-9, math.log(kbc) + 1e-9
        x = best.x.copy()
        x[0] = math.log(kbc)
        x[1] = math.log(min(29.0, max(0.11, kbc / 8e-12)))
        sol = least_squares(lambda y: residuals(make(y)), x, bounds=(lo, hi),
                            diff_step=1e-3, max_nfev=600)
        c = make(sol.x)
        print(f"  {kbc:.0e}  {c.K_eq:5.2f}  {c.tau_C * 1e9:3.0f} ns  "
              "%6.1f   %6.1f   %6.1f    %6.1f" % rms_parts(sol.fun)
              + f"      {b_share_1atm(c):.2f}")


if __name__ == "__main__":
    main()
