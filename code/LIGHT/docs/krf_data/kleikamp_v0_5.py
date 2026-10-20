"""Kleikamp, Peters and Witteman's coaxial e-beam KrF laser with krypton
as the only buffer (IEEE J. Quantum Electron. 22, 219 (1986)), as a
held-out structure check on LIGHT: where the Kr pressure and the F2
pressure optima of a Kr/F2 laser lie, and the effective three-body
quenching of KrF* by Kr the laser implies.

The measurements (their Figs. 3 and 4, Table I): with F2 at 10 mbar the
output energy of the two-component Kr-F2 system peaks at about 2 atm of
Kr; with Kr at 2 atm the output peaks at about 10 mbar of F2; and from
the optimum, with their steady-state cavity balance, they deduce a
three-body quenching constant of KrF* by Kr, k_R = 4.6e-31 cm6/s (their
Table I; "the quenching constants of Kr and Ne are not known in the
literature"). The diode: 280 kV, an average current density of about
250 A/cm2, 30-ns pulses, electrons of at most 170 keV in the gas, a
1-cm-diameter tube 20 cm long.

A score only, never a fit target (the predictive-only policy). LIGHT has
no cavity. Two proxies for the output are compared with the measured
optimum's location: (i) the net small-signal gain g0 - alpha0 over the
pulse's second half, and (ii) the paper's own steady-state balance of a
saturated laser (their Eqs. (1)-(3)): with the inversion clamped at
threshold, [KrF]_thr = (gamma + alpha)/sigma_s, the output per length is
gamma I = gamma/(gamma + alpha) (P - Q [KrF]_thr), P the KrF(B) formation
rate, Q the loss rate of a KrF(B) molecule (quenching and radiation),
gamma = 5.75e-2 1/cm their outcoupling and sigma_s LIGHT's 2.6e-16 cm2.
P, Q and alpha are LIGHT's, averaged over the window. The optimum of
(ii) balances the pumping, which rises with the density, against the
quenching of the threshold inversion, which rises faster (the trimer
goes as the density squared), so it tests the trimer's rate at the
laser's gas temperature, the largest term of v0_4's propagated
uncertainty.

The deposited power. For a 170-keV electron the collision stopping power
of krypton is about 2.2 MeV cm2/g (Berger and Seltzer's tables, the
source the paper names; the value here is read to +-15 %), which at 2 atm
(6.98e-3 g/cm3) and 250 A/cm2 gives 3.8 MW/cm3. The paper's own pumping
rate B_Kr = 111 cm2 A^-1 s^-1 per atom, read as a W-value conversion of
the same stopping power with W(Kr) = 24.1 eV, gives 5.6 MW/cm3. The scan
runs at their mean, 4.7 MW/cm3 at 2 atm, proportional to the gas
density, and again at 0.7 and 1.4 times it, since the optimum's location
is the observable and its sensitivity to the power is part of the answer.

    LIGHT_VERSION=v0_5 python3 LIGHT/docs/krf_data/kleikamp_v0_5.py [--kinetic] [--json PATH]
"""

import argparse
import importlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIGHT = HERE.parents[1]
KM = LIGHT.parent
sys.path.insert(0, str(KM))
sys.path.insert(0, str(LIGHT))

VERSION = os.environ.get("LIGHT_VERSION", "v0_5")
MEASURED = {"kr_opt_atm": 2.0, "f2_opt_mbar": 10.0, "k_R_cm6": 4.6e-31,
            "current_A_cm2": 250.0, "pulse_ns": 30.0, "e_max_keV": 170.0}
P_2ATM_W_CM3 = 4.7e6                # at 2 atm, see the docstring
P_SPREAD = (0.7, 1.0, 1.4)
KR_SCAN_ATM = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
F2_SCAN_MBAR = (2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0)
WINDOW = (15e-9, 30e-9)             # the pulse's second half
MBAR_PER_ATM = 1013.25
GAMMA_CM = 5.75e-2                  # their outcoupling per unit length


def _run(forward, drive, optics, p_atm, f2_mbar, power_scale, kinetic,
         trimer_scale=1.0):
    from_v = importlib.import_module(f"{VERSION}.model")
    Params = from_v.Params
    x_f2 = f2_mbar / (p_atm * MBAR_PER_ATM)
    closure = {} if kinetic else {"eedf": "two_t"}
    extra = ({} if trimer_scale == 1.0
             else {"row_scale": (("G5", trimer_scale), ("G5C", trimer_scale))})
    p = Params(network="basic", pressure_atm=p_atm, x_f2=x_f2,
               power_W_cm3=P_2ATM_W_CM3 * power_scale * p_atm / 2.0,
               pulse_s=MEASURED["pulse_ns"] * 1e-9,
               t_end_s=MEASURED["pulse_ns"] * 1e-9 + 20e-9, ne_mode="free",
               ne_cm3=1e10, n_save=601, **closure, **extra)
    r = drive.BeamModel(p).run()
    t = r.t
    win = (t >= WINDOW[0]) & (t <= WINDOW[1])
    g0 = 100.0 * optics.SIGMA_SE_CM2 * float(np.mean(r.n["KrF(B)"][win]))
    al = sum(100.0 * sig * float(np.mean(r.n[s][win]))
             for s, (sig, _) in optics.ABSORBERS.items() if s in r.n)
    tg = float(np.mean(r.tg[win]))
    m = r.model
    # the formation rate P of KrF(B) [cm-3 s-1] and the loss rate Q of a
    # KrF(B) molecule [1/s] (every row that makes or consumes it, the
    # radiative ones included), averaged over the window's saved states
    iB = m.idx["KrF(B)"]
    P, Q, nsteps = 0.0, 0.0, 0
    for k in np.where(win)[0]:
        n = np.array([r.n[s][k] for s in m.species])
        te, ne = float(r.te[k]), float(r.ne[k])
        nb = max(float(n[iB]), 1.0)
        for c in m.compiled:
            if c.kind == "elastic":
                continue
            made = c.row.products.count("KrF(B)") - c.row.reactants.count("KrF(B)")
            if made == 0:
                continue
            # the rate constant times the reactants' densities: events
            # per cm3 per s
            rate = m._rate(c, te, float(r.tg[k]), n, ne) * ne ** c.ne_pow
            for i in c.heavy:
                rate *= n[i]
            if made > 0:
                P += made * rate
            else:
                Q += -made * rate / nb
        nsteps += 1
    P /= nsteps
    Q /= nsteps
    alpha_cm = al / 100.0
    thr = (GAMMA_CM + alpha_cm) / optics.SIGMA_SE_CM2
    out_sat = GAMMA_CM / (GAMMA_CM + alpha_cm) * (P - Q * thr)
    # the effective three-body quenching of KrF(B) by Kr at the window's
    # gas temperature: the trimer row G5 with its temperature law
    c5 = next(c for c in m.compiled if c.row.rid == "G5")
    n = np.array([float(np.mean(r.n[s][win])) for s in m.species])
    k_eff = m._k(c5, float(np.mean(r.te[win])), tg, n, float(np.mean(r.ne[win])))
    return {"p_atm": p_atm, "f2_mbar": f2_mbar, "power_scale": power_scale,
            "g0": g0, "alpha": al, "net": g0 - al, "tg_K": tg,
            "te_eV": float(np.mean(r.te[win])), "ne": float(np.mean(r.ne[win])),
            "P_form": P, "Q_loss": Q, "n_thr": thr, "out_sat": out_sat,
            "k_trimer_eff": k_eff, "e_dep_J_cm3": float(r.ledger["E_dep"][-1]) * 1.602176634e-19}


def _peak(xs, ys):
    """The abscissa of the maximum, refined by the parabola through the
    three points around it (the scan's points are not equally spaced)."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    i = int(np.argmax(ys))
    if 0 < i < len(xs) - 1:
        a, b, _ = np.polyfit(xs[i - 1:i + 2], ys[i - 1:i + 2], 2)
        if a < 0:
            return float(-b / (2 * a))
    return float(xs[i])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kinetic", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--trimer-scale", type=float, default=1.0,
                    help="a sensitivity factor on the trimer rows G5, G5C "
                    "(Params.row_scale), e.g. 0.45 brings LIGHT's three-body "
                    "quenching at the laser's Tg to Kleikamp's deduced k_R; "
                    "a comparison, never adopted")
    ap.add_argument("--only-nominal", action="store_true",
                    help="the nominal power only, no F2 scan")
    a = ap.parse_args(argv)
    from ML.ml import forward
    assert forward.LIGHT_VERSION == VERSION, forward.LIGHT_VERSION
    drive = forward._drive()
    optics = importlib.import_module(f"{VERSION}.optics")
    out = {"version": VERSION, "kinetic": a.kinetic, "measured": MEASURED,
           "kr_scan": {}, "f2_scan": []}
    t0 = time.time()
    print(f"LIGHT {VERSION} at Kleikamp 1986's Kr/F2 laser "
          f"({'kinetic' if a.kinetic else 'two_t'} closure):")
    if a.trimer_scale != 1.0:
        print(f"  (trimer rows G5, G5C scaled by {a.trimer_scale:g}: a sensitivity run)")
    for ps in ((1.0,) if a.only_nominal else P_SPREAD):
        rows = [_run(forward, drive, optics, p, 10.0, ps, a.kinetic, a.trimer_scale)
                for p in KR_SCAN_ATM]
        out["kr_scan"][str(ps)] = rows
        nets = [r["net"] for r in rows]
        pk = _peak(np.array(KR_SCAN_ATM), np.array(nets))
        pk_sat = _peak(np.array(KR_SCAN_ATM), np.array([r["out_sat"] for r in rows]))
        print(f"  power x{ps:.1f}: g0-alpha0 over {WINDOW[0]*1e9:.0f}-{WINDOW[1]*1e9:.0f} ns "
              f"peaks at {pk:.2f} atm Kr; the saturated output at {pk_sat:.2f} atm "
              f"(measured output peak ~{MEASURED['kr_opt_atm']} atm)")
        for r in rows:
            print(f"     {r['p_atm']:4.2f} atm: g0 {r['g0']:6.2f} alpha {r['alpha']:5.2f} net {r['net']:6.2f} %/cm"
                  f"  P {r['P_form']:.2e} Q {r['Q_loss']:.2e} out_sat {r['out_sat']:.2e}"
                  f"  Tg {r['tg_K']:5.0f} K  Te {r['te_eV']:4.2f} eV  n_e {r['ne']:.2e}"
                  f"  k_trimer(Tg) {r['k_trimer_eff']:.2e}  E_dep {r['e_dep_J_cm3']:.3f} J/cm3")
    if a.only_nominal:
        at2 = next(r for r in out["kr_scan"]["1.0"] if r["p_atm"] == 2.0)
        print(f"  effective KrF* + 2Kr quenching at 2 atm (Tg {at2['tg_K']:.0f} K): "
              f"{at2['k_trimer_eff']:.2e} cm6/s (Kleikamp {MEASURED['k_R_cm6']:.1e})")
        if a.json:
            Path(a.json).write_text(json.dumps(out, indent=1))
        return 0
    rows = [_run(forward, drive, optics, 2.0, f, 1.0, a.kinetic, a.trimer_scale)
            for f in F2_SCAN_MBAR]
    out["f2_scan"] = rows
    pk = _peak(np.array(F2_SCAN_MBAR), np.array([r["net"] for r in rows]))
    pk_sat = _peak(np.array(F2_SCAN_MBAR), np.array([r["out_sat"] for r in rows]))
    print(f"  F2 at 2 atm Kr: g0-alpha0 peaks at {pk:.1f} mbar F2; the saturated output at "
          f"{pk_sat:.1f} mbar (measured ~{MEASURED['f2_opt_mbar']})")
    for r in rows:
        print(f"     {r['f2_mbar']:5.1f} mbar: g0 {r['g0']:6.2f} alpha {r['alpha']:5.2f} net {r['net']:6.2f} %/cm"
              f"  P {r['P_form']:.2e} Q {r['Q_loss']:.2e} out_sat {r['out_sat']:.2e}  Tg {r['tg_K']:5.0f} K")
    at2 = next(r for r in out["kr_scan"]["1.0"] if r["p_atm"] == 2.0)
    print(f"  effective KrF* + 2Kr quenching at 2 atm, 10 mbar (Tg {at2['tg_K']:.0f} K): "
          f"{at2['k_trimer_eff']:.2e} cm6/s; Kleikamp's deduced k_R {MEASURED['k_R_cm6']:.1e}"
          f" (x{at2['k_trimer_eff'] / MEASURED['k_R_cm6']:.2f})")
    print(f"  {time.time() - t0:.0f} s")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
