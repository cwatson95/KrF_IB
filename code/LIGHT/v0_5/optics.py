"""The 248-nm optics of the gas (v0_4, 2026-10-08): small-signal gain and
absorption from the densities, species by species.

g0(t) = sigma_B n_B(t) is what a probe measures; alpha(t) = sum_i sigma_i
n_i(t) the absorption the same probe sees (Tahoma's Table 3-3 prints both,
g0 6.39 and alpha0 1.11 %/cm). LIGHT has no cavity and no photon field, so
the absorbers change no density here: they are a diagnostic, the
coefficient the experiments report. In a laser the absorbed photons would
make F + e (F-), Kr+ + Kr (Kr2+), Kr+ + e (Kr**) and heat, which a cavity
model would book.

Cross sections [cm2] at 248 nm, Kannari, Obara & Fujioka, JAP 57, 4309
(1985) Table I p.4312 ("248-nm absorption cross section"), with his
reference numbers; the stimulated-emission cross section from the same
table. A species with no printed cross section is not an absorber here.
"""

from __future__ import annotations

import numpy as np

SIGMA_SE_CM2 = 2.6e-16          # KrF(B-X) stimulated emission, Kannari
#                                 Table I [46] (Tellinghuisen 1976), the
#                                 value ML's forward model and REBUILD use

# species: (sigma [cm2], source). Kannari's Table I block re-read from the
# page (2026-10-08): each entry's reference and its nature, and the one
# direct 248-nm measurement since, Hakuta et al. 1987 for Kr2F*, which
# replaces the compilation's estimate (a measurement before an estimate).
# Kr2*, Ar2* and F are not in his block and are not absorbers here.
ABSORBERS = {
    "F-": (5.6e-18, "F- + hv -> F + e, photodetachment, Kannari 1985 Table "
           "I [82] = Mandl, PRA 3, 251 (1971), a measurement"),
    "Kr2+": (1.6e-18, "Kr2+ + hv -> Kr+ + Kr, Kannari 1985 Table I [62] = "
             "Rokni, Jacob & Mangano, PRA 16, 2216 (1977); Rokni et al., "
             "IEEE JQE 14, 464 (1978) Table IV label the same value "
             "'estimated'"),
    "Kr**": (6.0e-18, "Kr(5p) + hv -> Kr+ + e, photoionization, Kannari "
             "1985 Table I [84] = Hyman, APL 31, 14 (1977), a calculation"),
    "F2": (1.2e-20, "F2 + hv -> 2F, Kannari 1985 Table I [81] = Calvert & "
           "Pitts, Photochemistry (1966), measured"),
    "Kr2F*": (1.0e-18, "Kr2F* + hv, (1.0 +- 0.2)e-18 cm2 at 248 nm, "
              "measured by fluorescence suppression under KrF laser light "
              "(Hakuta, Komori, Mukai & Takuma, JAP 61, 2113 (1987), "
              "abstract: 'about 1/5 of the value assumed in the current "
              "kinetic model'). Kannari 1985 Table I [16] = Johnson & "
              "Hunter 1980 carried 6.0e-18, an estimate, which REBUILD "
              "still does; Geohegan & Eden, CPL 139, 519 (1987) measured "
              "the 4^2 Gamma absorption from 248 to 570 nm; Edwards et al., "
              "APL 38, 843 (1981) give 1.8 +- 1.1e-17 at 257 nm"),
    "Kr*": (1.0e-19, "Kr(5s) + hv -> Kr+ + e, Kannari 1985 Table I [16] = "
            "Johnson & Hunter 1980, a compilation"),
}
# sigma_se: Kannari Table I [46] = Dunning & Hay, APL 28, 649 (1976), a
# calculation; measured 2.5 +- 0.3e-16 from gain and absorption (Rokni et
# al. 1978 p.476) and 2.6 +- 0.3e-16 from the spontaneous spectrum
# (Tellinghuisen et al. 1976)


def gain_cm(res) -> np.ndarray:
    """Small-signal gain g0(t) [1/cm] = sigma_B n_B."""
    return SIGMA_SE_CM2 * res.n["KrF(B)"]


def absorption_cm(res) -> dict:
    """alpha(t) [1/cm] per absorber present in the run, and 'total'."""
    out = {}
    tot = np.zeros_like(res.t, dtype=float)
    for s, (sig, _) in ABSORBERS.items():
        if s in res.n:
            a = sig * res.n[s]
            out[s] = a
            tot = tot + a
    out["total"] = tot
    return out


def window(res, t0=None, t1=None) -> dict:
    """g0, alpha (total and by species) and g0 - alpha averaged over
    [t0, t1] (default: the second half of the pulse), plus their values at
    the time of peak KrF(B), all in %/cm."""
    t = res.t
    if t0 is None:
        t0 = 0.5 * res.p.pulse_s
    if t1 is None:
        t1 = res.p.pulse_s
    win = (t >= t0) & (t <= t1)
    if not win.any():
        win = np.zeros_like(t, dtype=bool)
        win[-1] = True
    g = gain_cm(res)
    al = absorption_cm(res)
    k = int(np.argmax(res.n["KrF(B)"]))
    mean = lambda x: 100.0 * float(np.mean(x[win]))           # noqa: E731
    out = {"t0": float(t[win][0]), "t1": float(t[win][-1]),
           "g0": mean(g), "alpha": mean(al["total"]),
           "net": mean(g - al["total"]),
           "by_species": {s: mean(a) for s, a in al.items()
                          if s != "total"},
           "t_peak": float(t[k]), "g0_peak": 100.0 * float(g[k]),
           "alpha_at_peak": 100.0 * float(al["total"][k])}
    return out


def lines(res) -> list:
    """The report's optics block."""
    w = window(res)
    out = [f"OPTICS at 248 nm, mean over {w['t0'] * 1e9:.0f}-"
           f"{w['t1'] * 1e9:.0f} ns (sigma_B {SIGMA_SE_CM2:.1e} cm2; "
           "absorbers from Kannari 1985 Table I; a diagnostic, no cavity)",
           f"  g0 {w['g0']:.3f} %/cm   alpha {w['alpha']:.3f} %/cm   "
           f"g0 - alpha {w['net']:.3f} %/cm   "
           f"(at peak KrF(B), {w['t_peak'] * 1e9:.0f} ns: g0 "
           f"{w['g0_peak']:.3f}, alpha {w['alpha_at_peak']:.3f})"]
    parts = sorted(w["by_species"].items(), key=lambda kv: -kv[1])
    out.append("  alpha by species: " + ", ".join(
        f"{s} {v:.3f}" for s, v in parts if v > 0.0) + " %/cm")
    return out
