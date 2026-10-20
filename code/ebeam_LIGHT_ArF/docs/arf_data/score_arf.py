"""LIGHT_ArF_ebeam v0_2 (0-D) and 1Dv0_2 (1-D) against ArF measurements.

    python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py          # 0-D: score_0d.json
    FLUID1D_SRC=... ~/.venvs/f2sweep/bin/python \\
        LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --one-d      # 1-D: score_1d.json
    python3 LIGHT_ArF_ebeam/docs/arf_data/score_arf.py --report # SCORE.md

The measurements are SCORING_DATA.md's, their numbers re-read first-hand
where SCORE.md says so; the figure readings are fig_readings/*.csv. Nothing
here sets a rate: every model value comes from the network's defaults (the
bands at their defaults), with LIGHT's own input (a given power into the
electron gas) and, beside it, the model's electron-beam drive (v0_2
drive='beam', beam.py: REBUILD's beam physics with the Opal-Peterson-Beaty
secondaries). In LIGHT_ArF the beam column was the prototype drive
(beam_drive.py, which here only ar_tail.py still uses: the tail table was
made under it and is indexed by n_e/N, not by the drive). The seeded starts
(SRI, Chen & Payne) take the drive's excited atoms per ion and secondaries'
energy at their fill. Every input is otherwise LIGHT_ArF's, so the two
records differ by the drive alone.

What each test computes is in its docstring. Three kinds:
  device    formation efficiency and gain at the gain-measurement pumps
            (held out: no rate came from these experiments)
  shape     relative fluorescence against pressure, pump and F2 (held out)
  yield     SRI's photons per initial excitation after a 3-ns pulse
            (held out; seeded start)
  pure      Morozov's pure-Ar VUV efficiency (held out)
  kinetic   decay rates of Ar*, ArF* and Ar2F* against [F2] and [Ar]: the
            network's constants are means of these same groups, so these
            are consistency checks, not held-out tests
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIGHT_ARF = next(p for p in HERE.parents if (p / "run_v0_2.py").is_file())
if str(LIGHT_ARF) not in sys.path:
    sys.path.insert(0, str(LIGHT_ARF))

from v0_2 import beam as BEAM  # noqa: E402
from v0_2 import network as NET  # noqa: E402
from v0_2.model import QE, Model, Params  # noqa: E402
from v0_2.species import MANIFOLD  # noqa: E402

FIG = HERE / "fig_readings"
TORR = 3.2188e16                 # cm^-3 per Torr at 300 K
ATM = 760.0 * TORR
HV_BX = NET.HV_BX                # 6.414 eV
SIGMA_CALC = 2.75e-16            # cm^2, Lee 1989 Table I (calculated)
OUT0, OUT1 = HERE / "score_0d.json", HERE / "score_1d.json"


def fig(name):
    rows = []
    with open(FIG / name) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(lines):
        rows.append({k: (float(v) if k not in ("species", "series", "figure",
                                               "x_unit", "probe", "band_nm")
                         else v) for k, v in r.items()})
    return rows


# ---------------------------------------------------------------------------
# 0-D runs
# ---------------------------------------------------------------------------
# the electron closure every test runs on: 'auto' is the network's (the
# two-temperature tail); the *_maxwellian sensitivity entries set it, and
# run on LSODA: on the Maxwellian EEDF BDF's finite-difference Jacobian
# stalls at a Jancaitis point for tens of minutes, and the two solvers give
# the same peaks to five digits where both finish
EEDF = "auto"
METHOD = "BDF"


def params(p_atm, x_f2, power, pulse, t_end=None, **kw):
    return Params(pressure_atm=p_atm, x_f2=x_f2, power_W_cm3=power,
                  pulse_s=pulse, t_end_s=t_end or pulse,
                  ne_cm3=kw.pop("ne_cm3", 1e14), n_save=kw.pop("n_save", 801),
                  eedf=kw.pop("eedf", EEDF), method=kw.pop("method", METHOD),
                  **kw)


def run0(p, beam=False):
    m = Model(p.replace(drive="beam") if beam else p)
    return m, m.run()


def formed(m, r):
    """Cumulative ArF(B) formed [cm^-3]."""
    out = np.zeros_like(r.t)
    for c in m.compiled:
        net = (sum(c.row.products.count(s) for s in MANIFOLD)
               - sum(c.row.reactants.count(s) for s in MANIFOLD))
        if net > 0:
            out += net * r.events[c.row.rid]
    return out


def slope(t, y, w):
    sel = (t >= w[0]) & (t <= w[1])
    return float(np.polyfit(t[sel], y[sel], 1)[0])


def device0(p, window, beam):
    """Formation efficiency (ArF(B) formed x 6.414 eV over the deposited
    power) and the mean ArF(B) density in the window."""
    m, r = run0(p, beam)
    eff = slope(r.t, formed(m, r), window) * HV_BX * QE / p.power_W_cm3
    sel = (r.t >= window[0]) & (r.t <= window[1])
    return {"eff": eff, "nB": float(np.mean(r.n["ArF(B)"][sel])),
            "te": float(np.mean(r.te[sel])), "ne": float(np.mean(r.ne[sel])),
            "f2_left": float(r.n["F2"][sel][-1] / r.n["F2"][0])}


def emission(m, r, rid):
    """Emission rate [photons cm^-3 s^-1] of a radiative row: its
    emitter's density times its rate constant. Exact, where the derivative
    of the cumulative event counter picks up spikes from the solver's
    output grid (2026-09-30: two Ar2F* peaks came out equal that way)."""
    c = next(x for x in m.compiled if x.row.rid == rid)
    assert c.kind in ("const", "inv_tau", "param") and c.row.is_rad, rid
    return r.n[c.row.reactants[0]] * c.rate


# ---------------------------------------------------------------------------
# The tests: each returns a list of points {x, measured, model..., note}
# ---------------------------------------------------------------------------
ELECTRA = dict(p_atm=0.82, x_f2=0.003, pulse=120e-9)   # 100-ns flat + ramps
ELECTRA_WIN = (40e-9, 80e-9)       # the flat top's middle: "at temporal
                                   # peak power deposition"


def t_electra_formation(beam):
    """Electra (Wolford 2020 p.2, Obenschain 2020 p.7): 0.82 atm, 0.3 %
    F2, 1.09 MW/cc: g0 3.21 %/cm, Isat 10 MW/cm2, so g0 Isat / P = 0.29, the
    ArF(B) formation efficiency (hv times the B loss rate, which is its
    formation rate at the quasi-steady peak). Petrova's Fig. 7 plots the
    same point at 0.90 MW/cc: 0.36."""
    out = []
    for power, meas in ((1.09e6, 0.0321 * 10.0 / 1.09),
                        (0.90e6, 0.0321 * 10.0 / 0.90)):
        d = device0(params(ELECTRA["p_atm"], ELECTRA["x_f2"], power,
                           ELECTRA["pulse"]), ELECTRA_WIN, beam)
        out.append({"x": power / 1e6, "measured": meas, "model": d["eff"],
                    "state": d})
    return out


def t_electra_yield_band(beam):
    """Sensitivity, not a score of the defaults: Electra's formation
    efficiency with the ArF* yield of Ar* + F2 at the top of its band, 1.0
    (every published ArF model's value), against the measured 0.29. The
    default, 0.5, is Setser's low-pressure measurement."""
    d = device0(params(ELECTRA["p_atm"], ELECTRA["x_f2"], 1.09e6,
                       ELECTRA["pulse"], rates=(("arf_yield", 1.0),)),
                ELECTRA_WIN, beam)
    return [{"x": 1.09, "measured": 0.0321 * 10.0 / 1.09, "model": d["eff"],
             "state": d}]


def t_electra_gain_trend(beam):
    """Electra single pass, 30-cm cathode (NRL ArF report p.8): g0 1.84,
    1.77, 1.29 %/cm at (0.82 atm, 0.395 MW/cc), (1.09, 0.497), (1.36,
    0.583). The ratios to the first point need no cross section; the model's
    are its ArF(B) density ratios. g0 with Lee's calculated 2.75e-16 cm2 is
    given as information."""
    pts = ((0.82, 0.395e6, 1.84), (1.09, 0.497e6, 1.77), (1.36, 0.583e6, 1.29))
    ds = [device0(params(p, 0.003, pw, ELECTRA["pulse"]), ELECTRA_WIN, beam)
          for p, pw, _ in pts]
    out = []
    for (p, pw, g), d in zip(pts, ds):
        out.append({"x": p, "measured": g / pts[0][2],
                    "model": d["nB"] / ds[0]["nB"],
                    "g0_calc_pct_cm": 100 * SIGMA_CALC * d["nB"],
                    "g0_meas_pct_cm": g, "state": d})
    return out


def t_electra_gain_pump(beam):
    """Electra at 0.82 atm: 3.21 %/cm at 1.09 MW/cc (10-cm cathode, double
    pass) against 1.84 %/cm at 0.395 (30-cm, single pass): ratio 1.745. The
    two configurations differ, so this ratio carries their systematics."""
    lo = device0(params(0.82, 0.003, 0.395e6, ELECTRA["pulse"]),
                 ELECTRA_WIN, beam)
    hi = device0(params(0.82, 0.003, 1.09e6, ELECTRA["pulse"]),
                 ELECTRA_WIN, beam)
    return [{"x": 1.09 / 0.395, "measured": 3.21 / 1.84,
             "model": hi["nB"] / lo["nB"]}]


def t_keio_formation(beam):
    """Keio (Lee et al. 1989 text p.2053 and Fig. 5(a)): 1 atm, 0.4 % F2,
    65-ns pulse; g0 4.18 %/cm, Is 3.2 MW/cm2 (figure reading) at the
    printed 0.51 MW/cm3: g0 Is / P = 0.26 (NRL replots the point at 1.21
    MW/cm3, which would make it 0.11)."""
    out = []
    for power in (0.51e6, 1.21e6):
        d = device0(params(1.0, 0.004, power, 65e-9), (30e-9, 60e-9), beam)
        out.append({"x": power / 1e6, "measured": 0.0418 * 3.2e6 / power,
                    "model": d["eff"], "state": d})
    return out


def t_sandia_formation(beam):
    """Sandia (Adams & Dillon 1991, abstract): about 5 MW/cm3, 963 Torr
    (953 Ar + 8 F2, second-hand), g0 12.1 %/cm, Is 6.7 MW/cm2: 0.16. The
    pulse length is not in the abstract: 30 ns, window 10-25 ns."""
    d = device0(params(963 / 760, 8 / 963, 5e6, 30e-9), (10e-9, 25e-9), beam)
    return [{"x": 5.0, "measured": 0.121 * 6.7 / 5.0, "model": d["eff"],
             "state": d}]


def steady_emission(p, beam, window, rid="G4", te=False):
    """Mean emission rate in the window (and, with te=True, the mean Te
    there: LIGHT's own input at a low pump lets n_e decay and Te climb)."""
    m, r = run0(p, beam)
    sel = (r.t >= window[0]) & (r.t <= window[1])
    em = float(np.mean(emission(m, r, rid)[sel]))
    return (em, float(np.mean(r.te[sel]))) if te else em


# W/cm3 per atm of Ar at Rokni's 5 A/cm2 and 150 keV (printed; the power
# density is not): the survey's estimate, kept as LIGHT_ArF has it so that
# only the drive differs between the two records. The drive's Bethe/Moller
# stopping power gives 18.8 on a straight ray, the lower bound the O(10 %)
# detour raises (beam.deposited_power_W_cm3, tests/test_beam.py)
ROKNI_P_PER_ATM = 20e3


def t_rokni_vs_ar(beam):
    """Rokni, Jacob & Mangano 1977 Fig. 4: 193-nm fluorescence at 2 Torr
    F2 and 5 A/cm2 against Ar, 300-ns pulse (steady state). The beam's
    deposition scales with the Ar density, 20 kW/cm3 per atm (an estimate,
    18.8 by the drive's stopping power on a straight ray). Both curves
    are normalised at 1 atm. The Te column is the model's mean Te over
    250-300 ns."""
    data = fig("rokni1977_PRA_fig4_Ar.csv")
    ps = sorted({round(d["ar_atm"], 2) for d in data})
    ref_meas = np.mean([d["I_rel"] for d in data if abs(d["ar_atm"] - 1) < 0.01])
    model, te = {}, {}
    for pa in ps + [1.0]:
        tot = pa + 2 / 760
        model[pa], te[pa] = steady_emission(
            params(tot, (2 / 760) / tot, ROKNI_P_PER_ATM * pa, 300e-9),
            beam, (250e-9, 300e-9), te=True)
    out = []
    for pa in ps:
        meas = np.mean([d["I_rel"] for d in data
                        if round(d["ar_atm"], 2) == pa]) / ref_meas
        out.append({"x": pa, "measured": meas,
                    "model": model[pa] / model[1.0], "te": te[pa]})
    return out


def t_rokni_current_ratio(beam):
    """Rokni 1977 Fig. 5: the 193-nm fluorescence ratio for a x25 beam
    current change, about 25 at every pressure (the ion channel forms ArF*
    with branching near one). The Te columns are the model's mean Te over
    250-300 ns at the full current and at 1/25 of it."""
    data = fig("rokni1977_PRA_fig5_current_ratio.csv")
    out = []
    for d in data[::3]:
        pa = d["ar_atm"]
        tot = pa + 2 / 760
        hi, te_hi = steady_emission(params(tot, (2 / 760) / tot,
                                           ROKNI_P_PER_ATM * pa, 300e-9),
                                    beam, (250e-9, 300e-9), te=True)
        lo, te_lo = steady_emission(params(tot, (2 / 760) / tot,
                                           ROKNI_P_PER_ATM * pa / 25, 300e-9),
                                    beam, (250e-9, 300e-9), te=True)
        out.append({"x": pa, "measured": d["R"], "model": hi / lo,
                    "te": te_hi, "te_low": te_lo})
    return out


JANC_P_PER_A = 0.5e6 / 100.0 * (1500 / 760)   # W/cm3 per A/cm2 at 1500 Torr


def t_jancaitis_scan(beam):
    """Jancaitis & Powell 1983 Fig. 3: peak ArF* and Ar2F* fluorescence
    against current density, 5 Torr F2 + 1500 Torr Ar, 50-ns pulse; the
    deposition is 25 mJ/cm3 per atm at 100 A/cm2 (0.5 MW/cm3 per atm). Each
    curve normalised at its point nearest 74 A/cm2. The Te column is the
    model's mean Te over 20-50 ns, the pulse after its start-up."""
    data = fig("jancaitis1983_fig3_peak_vs_current.csv")
    tot = 1505 / 760
    out = []
    for sp, rid in (("ArF", "G4"), ("Ar2F", "C3")):
        pts = [d for d in data if d["species"] == sp]
        ref = min(pts, key=lambda d: abs(d["J_A_cm2"] - 74))
        peak, te = {}, {}
        for d in pts:
            m, r = run0(params(tot, 5 / 1505, JANC_P_PER_A * d["J_A_cm2"],
                               50e-9, 200e-9), beam)
            peak[d["J_A_cm2"]] = float(np.max(emission(m, r, rid)))
            late = (r.t >= 20e-9) & (r.t <= 50e-9)
            te[d["J_A_cm2"]] = float(np.mean(r.te[late]))
        for d in pts:
            out.append({"x": d["J_A_cm2"], "species": sp,
                        "measured": d["peak_arb"] / ref["peak_arb"],
                        "model": peak[d["J_A_cm2"]] / peak[ref["J_A_cm2"]],
                        "te": te[d["J_A_cm2"]]})
    return out


# Seeded runs outlive their electrons: once n_e sits at its 1 cm^-3 floor
# the electron energy no longer moves anything, BDF's finite-difference
# Jacobian keeps growing that column's step and overflows (Ar 5 atm, 0.5 Torr
# F2 crashed). LSODA gives the same yields to 4 digits wherever both run
# (2026-09-30) and completes those.
SEED_METHOD = "LSODA"


class SeedModel(Model):
    """A 3-ns pulse taken as instantaneous: ions (and electrons carrying
    the secondaries' share each) and excited atoms at t = 0, no power. The
    share and the excited atoms per ion are the beam drive's at the fill
    (beam.yields); its prompt F2 dissociation (0.15 % of the pairs) and
    heat are left out of a seed."""

    def initial_state(self):
        y = super().initial_state()
        y[self.i_ue] = (self.p.ne_cm3 * BEAM.yields(self.p.x_f2).eps_sec
                        / self.N0)
        return y


def seed_params(tot_atm, x_f2, n_i, t_end, n_save=801):
    exc = BEAM.yields(x_f2).ar_star               # excited atoms per ion
    return Params(pressure_atm=tot_atm, x_f2=x_f2, power_W_cm3=0.0,
                  pulse_s=0.0, t_end_s=t_end, ne_cm3=n_i, method=SEED_METHOD,
                  init_extra_cm3=(("Ar*", exc * n_i),), n_save=n_save), exc


def sri_yield(p_ar, f2_torr):
    n_i = 1e15 * p_ar                        # "10^15 [Ar(atm)] cm^-3"
    tot = p_ar + f2_torr / 760
    p, exc = seed_params(tot, (f2_torr / 760) / tot, n_i, 3e-6)
    r = SeedModel(p).run()
    n_x = n_i * (1 + exc)
    return (float(r.events["G4"][-1] / n_x), float(r.events["C3"][-1] / n_x))


def t_sri_vs_ar(beam=False):
    """SRI MP-76-99 Fig. 4: 193- and 290-nm photons per excited state made
    initially (ions + excited atoms, p.25) after a 3-ns Febetron pulse, 1
    Torr F2, against Ar. Seeded start: 1e15 ions per atm (p.37), the beam
    drive's excited atoms per ion at the fill (0.39) and its secondaries'
    share (4.4 eV each) in the electron gas."""
    data = fig("sri1976_fig4_yields_vs_Ar.csv")
    out = []
    for d in data:
        y193, y284 = sri_yield(d["ar_atm"], 1.0)
        out.append({"x": d["ar_atm"], "band": d["band_nm"],
                    "measured": d["yield"],
                    "model": y193 if str(d["band_nm"]) == "193" else y284})
    return out


def t_sri_vs_f2(beam=False):
    """SRI MP-76-99 Fig. 11: the 193-nm yield against F2 at 1, 2, 5 atm."""
    data = fig("sri1976_fig11_ArF_yield_vs_F2.csv")
    out, cache = [], {}
    for d in data:
        key = (d["ar_atm"], d["f2_torr"])
        if key not in cache:
            cache[key] = sri_yield(*key)[0]
        out.append({"x": d["f2_torr"], "ar_atm": d["ar_atm"],
                    "measured": d["yield"], "model": cache[key]})
    return out


def t_morozov_vuv(beam):
    """Morozov et al. 2008: pure-Ar VUV power into 4 pi over deposited
    power, 0.33 +- 0.04 at 500-1400 hPa (printed), independent of power;
    Fig. 10 points. LIGHT_ArF on pure Ar at a low pump (330 W/cm3 for 450
    ns, as ML's pure-Kr check), photons counted to 20 us after the pulse."""
    data = fig("morozov2008_fig10_Ar_efficiency.csv")
    ps = sorted({d["p_hPa"] for d in data})
    out = []
    for ph in ps:
        pa = ph / 1013.25
        p = params(pa, 0.0, 330.0, 450e-9, 450e-9 + 20e-6, ne_cm3=1e8,
                   n_save=2001)
        m, r = run0(p, beam)
        e_dep = float(r.ledger["E_dep"][-1])
        vuv = float(r.events["AR2_RAD"][-1]) * NET.HV_AR2 / e_dep
        meas = np.mean([d["efficiency"] for d in data if d["p_hPa"] == ph])
        out.append({"x": ph, "measured": meas, "model": vuv})
    return out


# --- kinetic consistency (the constants' own sources) -----------------------
K_AR_F2 = NET.BANDS["ar_f2"].default
K8 = NET.BANDS["ar2_form"].default
AR2F_LOSS0 = 1.0 / NET.AR2F_TAU


def t_kinetic(beam=False):
    """Decay rates against the network's constants (analytic: the rates
    the model uses, summed). Ar*: k(Ar*+F2)[F2] + k(Ar*+2Ar)[Ar]^2. Ar2F*:
    1/tau + k(Ar2F*+F2)[F2] (+ Ar2F* + Ar, 0 by default). These sources
    gave the constants (as means), so this is consistency, not accuracy."""
    out = []
    for d in fig("chenpayne1979_fig13_ArF_decay.csv"):
        ar, f2 = d["p_ar_torr"] * TORR, d["f2_torr"] * TORR
        out.append({"set": "Chen & Payne Fig. 13, ArF* (follows Ar*)",
                    "x": d["f2_torr"], "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": K_AR_F2 * f2 + K8 * ar * ar})
    for d in fig("chenpayne1979_fig14_Ar2F_decay.csv"):
        out.append({"set": "Chen & Payne Fig. 14, Ar2F*", "x": d["f2_torr"],
                    "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": AR2F_LOSS0 + NET.K_AR2F_F2 * d["f2_torr"] * TORR})
    for d in fig("jancaitis1983_fig7_Ar2F_decay.csv"):
        f2 = d["p_total_torr"] / (1 + d["ar_to_f2"]) * TORR
        out.append({"set": f"Jancaitis Fig. 7, Ar2F*, Ar/F2 {d['ar_to_f2']:g}",
                    "x": d["p_total_torr"],
                    "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": AR2F_LOSS0 + NET.K_AR2F_F2 * f2})
    for d in fig("ueda1983_fig8_Ar2F_decay.csv"):
        out.append({"set": "Ueda 1983 Fig. 8, Ar2F*",
                    "x": d["f2_1e17_cm3"], "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": AR2F_LOSS0 + NET.K_AR2F_F2 * d["f2_1e17_cm3"] * 1e17})
    for d in fig("ueda1983_fig12_ArMeta_vs_F2.csv"):
        n = 2.5 * ATM
        f2 = d["f2_percent"] / 100 * n
        out.append({"set": "Ueda 1983 Fig. 12, Ar*(3P2) at 2.5 atm",
                    "x": d["f2_percent"], "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": K_AR_F2 * f2 + K8 * (n - f2) ** 2})
    for d in fig("ueda1983_fig10_ArMeta_pureAr.csv"):
        n = d["p_atm"] * ATM
        out.append({"set": "Ueda 1983 Fig. 10, Ar*(3P2) in pure Ar",
                    "x": d["p_atm"], "measured": d["decay_1e7_per_s"] * 1e7,
                    "model": K8 * n * n})
    return out


def t_stern_volmer(beam=False):
    """Stern-Volmer plots, 1/I(193) normalised at the first point: Chen &
    Payne Figs. 16 (F2 at 600 Torr Ar) and 17 (Ar^2 at 0.84 Torr F2), from a
    dilute seeded start (proton excitation, n_e < 1e12); Rokni Fig. 3 (F2
    at 150 Torr Ar) as a steady e-beam run. These gave the k tau products:
    consistency."""
    def seeded(p_ar_torr, f2_torr):
        n_i = 1e11
        tot = (p_ar_torr + f2_torr) / 760
        p, _ = seed_params(tot, f2_torr / 760 / tot, n_i, 4e-6)
        r = SeedModel(p).run()
        return float(r.events["G4"][-1])
    out = []
    sv = fig("chenpayne1979_fig16_17_SternVolmer.csv")
    for figno, key in (("16", "F2 at 600 Torr Ar"),
                       ("17", "Ar^2 at 0.84 Torr F2")):
        pts = [d for d in sv if d["figure"] == figno]
        vals = []
        for d in pts:
            if figno == "16":
                i = seeded(600.0, d["x"])
            else:
                i = seeded(math.sqrt(d["x"] * 1e5), 0.84)
            vals.append(1.0 / i)
        for d, v in zip(pts, vals):
            out.append({"set": f"Chen & Payne Fig. {figno}, 1/I vs {key}",
                        "x": d["x"], "measured": d["inv_I"] / pts[0]["inv_I"],
                        "model": v / vals[0]})
    rk = fig("rokni1977_PRA_fig3_F2.csv")
    xs = sorted({d["f2_torr"] for d in rk})
    model = {}
    for f2 in xs:
        tot = (150 + f2) / 760
        model[f2] = 1.0 / steady_emission(
            params(tot, f2 / 760 / tot, ROKNI_P_PER_ATM * 150 / 760, 300e-9),
            False, (250e-9, 300e-9))
    first = np.mean([d["inv_I"] for d in rk if d["f2_torr"] == xs[0]])
    for f2 in xs:
        meas = np.mean([d["inv_I"] for d in rk if d["f2_torr"] == f2]) / first
        out.append({"set": "Rokni Fig. 3, 1/I vs F2 at 150 Torr Ar",
                    "x": f2, "measured": meas, "model": model[f2] / model[xs[0]]})
    return out


def _maxwellian(fn):
    """fn on the Maxwellian EEDF, LIGHT's own input only: how much of an
    own-input result is the two-temperature tail's (sensitivity)."""
    def run(beam):
        global EEDF, METHOD
        old = EEDF, METHOD
        EEDF, METHOD = "maxwellian", "LSODA"
        try:
            return fn(beam)
        finally:
            EEDF, METHOD = old
    run.__doc__ = ("Sensitivity, not a score of the defaults: the test below "
                   "on the Maxwellian EEDF in place of the two-temperature "
                   "tail, LIGHT's own input only. " + fn.__doc__)
    return run


TESTS = {
    "electra_formation": ("device", t_electra_formation, True),
    "electra_yield_band": ("sensitivity", t_electra_yield_band, True),
    "electra_gain_trend": ("device", t_electra_gain_trend, True),
    "electra_gain_pump": ("device", t_electra_gain_pump, True),
    "keio_formation": ("device", t_keio_formation, True),
    "sandia_formation": ("device", t_sandia_formation, True),
    "rokni_vs_ar": ("shape", t_rokni_vs_ar, True),
    "rokni_current_ratio": ("shape", t_rokni_current_ratio, True),
    "jancaitis_scan": ("shape", t_jancaitis_scan, True),
    "sri_vs_ar": ("yield", t_sri_vs_ar, False),
    "sri_vs_f2": ("yield", t_sri_vs_f2, False),
    "morozov_vuv": ("pure", t_morozov_vuv, True),
    "kinetic_decays": ("kinetic", t_kinetic, False),
    "stern_volmer": ("kinetic", t_stern_volmer, False),
    "electra_formation_maxwellian": ("sensitivity",
                                     _maxwellian(t_electra_formation), False),
    "keio_formation_maxwellian": ("sensitivity",
                                  _maxwellian(t_keio_formation), False),
    "sandia_formation_maxwellian": ("sensitivity",
                                    _maxwellian(t_sandia_formation), False),
    "rokni_vs_ar_maxwellian": ("sensitivity", _maxwellian(t_rokni_vs_ar),
                               False),
    "rokni_current_ratio_maxwellian": (
        "sensitivity", _maxwellian(t_rokni_current_ratio), False),
    "jancaitis_scan_maxwellian": ("sensitivity",
                                  _maxwellian(t_jancaitis_scan), False),
}


def run_0d(only=None):
    out = {"model": "LIGHT_ArF_ebeam v0_2", "chemistry_md5": chem_md5(),
           "tests": {}}
    if only and OUT0.is_file():            # add to the tests already run
        old = json.loads(OUT0.read_text())
        if old.get("chemistry_md5") == out["chemistry_md5"]:
            out["tests"] = old["tests"]
    for name, (kind, fn, has_beam) in TESTS.items():
        if only and name not in only:
            continue
        t0 = time.time()
        res = {"kind": kind, "doc": " ".join(fn.__doc__.split()),
               "own": fn(False)}
        if has_beam:
            res["beam"] = fn(True)
        out["tests"][name] = res
        print(f"  {name:22s} {time.time() - t0:6.1f} s", flush=True)
    OUT0.write_text(json.dumps(out, indent=1, default=float) + "\n")
    return out


CHEM_FILES = ("model.py", "network.py", "species.py", "xsec.py",
              "data/ar_tail_ratio.txt", "beam.py")


def chem_md5(pkg="v0_2"):
    """The chemistry a score was made on (beam.py, the drive, is v0_2's
    only; 1Dv0_2 has none)."""
    import hashlib
    d = LIGHT_ARF / pkg
    return {f: hashlib.md5((d / f).read_bytes()).hexdigest()
            for f in CHEM_FILES if (d / f).is_file()}


# ---------------------------------------------------------------------------
# 1-D runs (1Dv0_2 on Fluid1D): the device and shape tests, on LIGHT's own
# input ('centre', 'gap') and under the electron beam ('centre_beam',
# 'gap_beam')
# ---------------------------------------------------------------------------
DRIVES_1D = (("power", ""), ("beam", "_beam"))


def run_1d(cells=64, only=None):
    sim1d = importlib.import_module("1Dv0_2.sim1d")
    runs = HERE / "runs"
    QE1 = QE

    def one(p_atm, x_f2, power, pulse, name, drive, t_end=None):
        p = sim1d.Params1D(pressure_atm=p_atm, x_f2=x_f2, power_W_cm3=power,
                           drive=drive, pulse_s=pulse,
                           t_end_s=t_end or pulse, ne_cm3=1e14, cells=cells,
                           dt_s=50e-12, n_save=401, n_profiles=3,
                           verbose=False, name=name, output_dir=str(runs))
        return sim1d.run(p)

    def device(res, power, window):
        t = res.t
        sel = (t >= window[0]) & (t <= window[1])
        out = {}
        for where, led, dens in (("centre", "ledger_centre", "centre"),
                                 ("gap", "ledger", "mean")):
            f = res.ledger("arf_formed", led)
            out[where] = {
                "eff": float(np.polyfit(t[sel], f[sel], 1)[0]) * HV_BX * QE1
                / power,
                "nB": float(np.mean(res.series("ArF(B)", dens)[sel]))}
        return out

    def both(make):
        """make(drive, suffix) -> {key: value}; the two drives merged."""
        out = {}
        for drive, suf in DRIVES_1D:
            out.update(make(drive, suf))
        return out

    tests = {}
    t0 = time.time()
    if not only or "electra_formation" in only:
        pts = []
        for power, meas in ((1.09e6, 0.0321 * 10 / 1.09),
                            (0.90e6, 0.0321 * 10 / 0.90)):
            def make(drive, suf):
                d = device(one(0.82, 0.003, power, 120e-9,
                               f"electra_{power:.3g}{suf}", drive),
                           power, ELECTRA_WIN)
                return {"centre" + suf: d["centre"]["eff"],
                        "gap" + suf: d["gap"]["eff"]}
            pts.append({"x": power / 1e6, "measured": meas, **both(make)})
        tests["electra_formation"] = pts
        print(f"  electra_formation 1-D {time.time() - t0:.0f} s", flush=True)
    if not only or "electra_gain_trend" in only:
        nb = {}
        for p, pw, g in ((0.82, 0.395e6, 1.84), (1.09, 0.497e6, 1.77),
                         (1.36, 0.583e6, 1.29)):
            for drive, suf in DRIVES_1D:
                d = device(one(p, 0.003, pw, 120e-9,
                               f"electra_trend_{p:.3g}{suf}", drive),
                           pw, ELECTRA_WIN)
                nb[(p, suf)] = (d["centre"]["nB"], d["gap"]["nB"], g)
        pts = []
        for p in (0.82, 1.09, 1.36):
            q = {"x": p, "measured": nb[(p, "")][2] / nb[(0.82, "")][2]}
            for _, suf in DRIVES_1D:
                q["centre" + suf] = nb[(p, suf)][0] / nb[(0.82, suf)][0]
                q["gap" + suf] = nb[(p, suf)][1] / nb[(0.82, suf)][1]
            pts.append(q)
        tests["electra_gain_trend"] = pts
        print(f"  electra_gain_trend 1-D {time.time() - t0:.0f} s", flush=True)
    if not only or "keio_formation" in only:
        pts = []
        for power in (0.51e6, 1.21e6):
            def make(drive, suf):
                d = device(one(1.0, 0.004, power, 65e-9,
                               f"keio_{power:.3g}{suf}", drive),
                           power, (30e-9, 60e-9))
                return {"centre" + suf: d["centre"]["eff"],
                        "gap" + suf: d["gap"]["eff"]}
            pts.append({"x": power / 1e6, "measured": 0.0418 * 3.2e6 / power,
                        **both(make)})
        tests["keio_formation"] = pts
        print(f"  keio_formation 1-D {time.time() - t0:.0f} s", flush=True)
    if not only or "sandia_formation" in only:
        def make(drive, suf):
            d = device(one(963 / 760, 8 / 963, 5e6, 30e-9, f"sandia{suf}",
                           drive), 5e6, (10e-9, 25e-9))
            return {"centre" + suf: d["centre"]["eff"],
                    "gap" + suf: d["gap"]["eff"]}
        tests["sandia_formation"] = [{"x": 5.0,
                                      "measured": 0.121 * 6.7 / 5.0,
                                      **both(make)}]
        print(f"  sandia_formation 1-D {time.time() - t0:.0f} s", flush=True)
    if not only or "rokni_vs_ar" in only:
        data = fig("rokni1977_PRA_fig4_Ar.csv")
        ps = [0.243, 0.494, 1.0, 2.0, 4.003]
        ref_meas = np.mean([d["I_rel"] for d in data
                            if abs(d["ar_atm"] - 1) < 0.01])
        em = {}
        for pa in ps:
            tot = pa + 2 / 760
            for drive, suf in DRIVES_1D:
                res = one(tot, (2 / 760) / tot, ROKNI_P_PER_ATM * pa, 300e-9,
                          f"rokni_{pa:.3g}{suf}", drive)
                t = res.t
                sel = (t >= 250e-9) & (t <= 300e-9)
                em[(pa, suf)] = {w: float(np.polyfit(
                    t[sel], res.ledger(f"rad:{NET.BX_BAND}", led)[sel], 1)[0])
                    for w, led in (("centre", "ledger_centre"),
                                   ("gap", "ledger"))}
        pts = []
        for pa in ps:
            meas = np.mean([d["I_rel"] for d in data
                            if abs(d["ar_atm"] - pa) < 0.01]) / ref_meas
            q = {"x": pa, "measured": meas}
            for _, suf in DRIVES_1D:
                for w in ("centre", "gap"):
                    q[w + suf] = em[(pa, suf)][w] / em[(1.0, suf)][w]
            pts.append(q)
        tests["rokni_vs_ar"] = pts
        print(f"  rokni_vs_ar 1-D {time.time() - t0:.0f} s", flush=True)
    out = {"model": "LIGHT_ArF_ebeam 1Dv0_2", "cells": cells, "gap_cm": 1.0,
           "chemistry_md5": chem_md5("1Dv0_2"), "tests": tests,
           "wall_s": time.time() - t0}
    OUT1.write_text(json.dumps(out, indent=1, default=float) + "\n")
    return out


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def gm(xs):
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def rms_log(pairs):
    v = [math.log(m / d) for d, m in pairs if d > 0 and m > 0]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--one-d", action="store_true", help="the 1-D runs")
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args(argv)
    if a.report:
        import make_score_report
        make_score_report.write()
        return 0
    if a.one_d:
        sys.path.insert(0, str(LIGHT_ARF))
        run_1d(a.cells, a.only)
    else:
        run_0d(a.only)
    return 0


if __name__ == "__main__":
    sys.exit(main())
