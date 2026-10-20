#!/usr/bin/env python3
"""Back-of-the-envelope: inverse-bremsstrahlung coupling vs frequency.

Owner ask (2026-08-21): "does microwave not work at a given frequency?
... what about calculating the inverse bremsstrahlung for a range of
microwave frequency?" — is the reheat failure a FREQUENCY problem
(coupling) or a KINETICS/WALL-PLUG problem?

Physics (electron-neutral inverse bremsstrahlung = the collisional/Drude
absorption the engine's window already uses):

  per-electron heating   P_e = e^2 E_rms^2 nu_m / (m (omega^2 + nu_m^2))
  at fixed intensity I   E_rms^2 = Z0 I  (traveling plane wave)
  absorption coeff       alpha = (omega_p^2 / c) nu_m / (omega^2 + nu_m^2)
                         [valid omega_p^2 << omega^2 + nu_m^2 — flagged]
  critical density       n_c = eps0 m omega^2 / e^2
  quiver/ponderomotive   U_p = e^2 E_rms^2 / (2 m (omega^2 + nu_m^2)) ~
                         the energy scale the field can hold per electron

nu_m comes from the engine's OWN swarm tables (nu = e E / (m v_d) at the
Ar/Kr/F2 fill, 2 atm) — the same quantity the window drive uses. Two
regimes:
  omega << nu_m  (all of 1-300 GHz at >= 1 atm): the AC factor is ~1 AND
                 alpha ~ omega_p^2/(c nu_m) is FREQUENCY-INDEPENDENT —
                 the plasma cannot tell 2.45 GHz from 95 GHz; coupling is
                 total either way. Frequency buys only penetration
                 bookkeeping (n_c) and SOURCE hardware.
  omega >> nu_m  (THz and IR, incl. the CO2 laser at 28.3 THz): heating
                 per unit intensity falls as omega^-2 — IR couples
                 WEAKLY per W, sources just have more W.

The wall-plug inequality that actually decides feasibility (measured
marginals from the 0.1.21-0.1.25 studies: absorbed-EM-joule -> laser-joule
conversion m in [-10 %, +12 %], sign tier-fragile):

  d(eta_wallplug) > 0   <=>   eta_source x m  >  eta_beam_wallplug (~7 %)

so even a PERFECTLY coupled window needs eta_source x m > 7 % — the table
prints eta_source x m_best per frequency. Source efficiencies are stated
literature-class numbers, not measurements. BOTE, scored against nothing.

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/ib_feasibility.py
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.MODULES.loki_tier import rates_for_mode  # noqa: E402
from REBUILD.MODULES.units import (C_LIGHT, E_CHARGE, EPS0, EV, M_E,  # noqa: E402
                                    gas_density_m3)

Z0 = 376.730                       # Ohm, vacuum impedance
M_BEST = 0.12                      # best measured laser marginal (fragile)
ETA_BEAM_WALL = 0.07               # DE-AR: current KrF wall-plug

# (label, f_Hz, source class, stated wall-plug efficiency of the source)
FREQS = [
    ("1 GHz (L-band)", 1e9, "magnetron/solid-state", 0.70),
    ("2.45 GHz (ISM magnetron)", 2.45e9, "magnetron", 0.72),
    ("9.5 GHz (X-band, WR90)", 9.5e9, "magnetron/klystron", 0.60),
    ("35 GHz (Ka gyrotron)", 35e9, "gyrotron", 0.40),
    ("95 GHz (DE-AR gyrotron)", 95e9, "gyrotron", 0.35),
    ("170 GHz (ITER-class gyrotron)", 170e9, "gyrotron", 0.40),
    ("300 GHz", 300e9, "gyrotron (edge)", 0.20),
    ("1 THz", 1e12, "THz source", 0.01),
    ("28.3 THz (CO2 laser, 10.6 um)", 28.3e12, "CO2 laser", 0.12),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pressure-atm", type=float, default=2.0)
    ap.add_argument("--en-td", type=float, default=1.0,
                    help="reference E/N for nu_m from the swarm tables")
    ap.add_argument("--intensity-W-cm2", type=float, default=1e3,
                    help="reference intensity for the per-electron heating "
                         "column (1 kW/cm^2 ~ 250 kW over a 10 cm cube face)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    n_tot = gas_density_m3(args.pressure_atm, 300.0)
    lok = rates_for_mode("krf", {"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
    # nu_m from the engine's own transport at several E/N (the window's
    # operating band) — nu = e E / (m v_d)
    nus = {}
    for en in (0.5, 1.0, 5.0, 30.0):
        e_f = en * 1e-21 * n_tot
        v_d = lok.drift_velocity(en)
        nus[en] = E_CHARGE * e_f / (M_E * v_d)
    nu = nus[args.en_td]
    print(f"[ib_feasibility] {args.pressure_atm:g} atm Ar/Kr/F2 fill; "
          f"nu_m from the swarm tables: "
          + ", ".join(f"{k:g} Td -> {v:.2e} /s" for k, v in nus.items()))
    print(f"reference nu_m = {nu:.2e} /s (at {args.en_td:g} Td) -> the "
          f"collisional/IB crossover f = nu/2pi = {nu / 2 / np.pi / 1e12:.2f} THz\n")
    I = args.intensity_W_cm2 * 1e4                 # W/m^2
    E2 = Z0 * I                                    # E_rms^2
    rows = []
    hdr = (f"{'frequency':>30} {'AC factor':>9} {'heat eV/ns/e-':>13} "
           f"{'U_p eV':>8} {'n_c m^-3':>9} {'L_abs@1e18':>10} "
           f"{'L_abs@1e20':>10} {'eta_src':>7} {'eta*m_best':>10}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for label, f, src, eta_s in FREQS:
        w = 2 * np.pi * f
        fac = nu * nu / (nu * nu + w * w)
        p_e = E_CHARGE ** 2 * E2 * nu / (M_E * (w * w + nu * nu))   # W per e-
        heat = p_e / EV * 1e-9                                      # eV/ns
        u_p = E_CHARGE ** 2 * E2 / (2 * M_E * (w * w + nu * nu)) / EV
        n_c = EPS0 * M_E * w * w / E_CHARGE ** 2
        L18 = L20 = None
        alpha = {}
        for ne, keyname in ((1e18, "1e18"), (1e20, "1e20")):
            wp2 = ne * E_CHARGE ** 2 / (EPS0 * M_E)
            a = (wp2 / C_LIGHT) * nu / (w * w + nu * nu)
            alpha[keyname] = a
            valid = wp2 < 0.3 * (w * w + nu * nu)
            alpha[keyname + "_valid"] = bool(valid)
        L18 = 1.0 / alpha["1e18"]
        L20 = 1.0 / alpha["1e20"]
        payoff = eta_s * M_BEST
        verdict = ("couples fully; source-limited" if fac > 0.5 else
                   "weak coupling per W (omega >> nu)")
        if payoff > ETA_BEAM_WALL:
            verdict += "; eta*m EXCEEDS beam 7% (if the fragile sign holds)"
        rows.append({"label": label, "f_Hz": f, "source": src,
                     "eta_source": eta_s, "ac_factor": fac,
                     "heat_eV_per_ns_per_e": heat, "U_p_eV": u_p,
                     "n_c_m3": n_c, "L_abs_1e18_m": L18, "L_abs_1e20_m": L20,
                     "alpha_valid_1e20": alpha["1e20_valid"],
                     "eta_times_m_best": payoff})
        print(f"{label:>30} {fac:9.4f} {heat:13.3g} {u_p:8.2e} "
              f"{n_c:9.2e} {L18:10.3g} {L20:10.3g} {100 * eta_s:6.0f}% "
              f"{100 * payoff:9.1f}%  {verdict}")
    print(f"""
READING (the BOTE conclusions):
 - Below ~{nu / 2 / np.pi / 1e12:.1f} THz the AC factor is ~1 AND the absorption
   length is frequency-independent (alpha ~ omega_p^2 / (c nu_m)): the
   plasma CANNOT tell the frequencies apart — every GHz band couples
   totally into a beam-made plasma. The reheat failure is NOT a frequency
   or coupling failure; it is kinetics (absorbed joules thermalize; the
   conductive plasma pins the field below excitation thresholds) plus the
   wall-plug inequality.
 - What frequency changes: n_c (penetration bookkeeping), the achievable
   field per watt of source, and ABOVE the crossover (THz, CO2 at
   28.3 THz) the coupling per watt falls as omega^-2 — the CO2 laser
   heats electrons only via its enormous intensity, and its ~12% source
   efficiency times the <=12% (fragile) laser marginal is ~1.4%, far
   below the beam's ~7% wall-plug: CO2 reheat CANNOT raise wall-plug
   efficiency (the owner's hunch, quantified).
 - The wall-plug inequality eta_source x marginal > 7% is satisfied at
   BEST marginally, by high-efficiency magnetrons (1-2.45 GHz, ~70%)
   times the +12% plateau co-pump marginal -> ~8.4% — IF that fragile,
   tier-flipping sign is real (0.1.23), and even then the absolute gain
   is capped at P_MW/P_beam ~ 1e-3 of the pulse. Every other line loses
   outright. The efficiency of this laser lives in the beam-only levers
   (0.1.24: F2 fraction, pressure, optics, Ne buffer).
 - The one EM-positive channel the benchmarks show — afterglow reheat of
   a HIGH-MARGIN cell that is still lasing during the window (Ne-opt
   +9.3% LOKI / Ar-opt +9.5% KINETIC; dE_laser >> dFluor there) — has
   the same mechanism as co-pump (feed a cavity that is extracting) and
   the same fragile-sign, wall-plug-capped arithmetic.""")
    out_dir = (Path(args.out) if args.out
               else ROOT / "rebuild_out" / "mw_benchmarks")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "ib_feasibility.json"
    path.write_text(json.dumps({"pressure_atm": args.pressure_atm,
                                "nu_m_per_s": nus, "rows": rows,
                                "m_best": M_BEST,
                                "eta_beam_wallplug": ETA_BEAM_WALL},
                               indent=2, default=float))
    print(f"\nreport -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
