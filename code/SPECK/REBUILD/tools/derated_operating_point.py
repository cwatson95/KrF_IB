"""Derated operating point for an afterglow-reheat efficiency gain
([0.1.33]): what survives when the coupling is NOT ideal?

Everything reads from recorded artifacts:
  m, eta0, E_dep, E_MW   <- mw_benchmarks_loki.json (ne_optimized,
                             afterglow — the one robust positive)
  nu_m(E/N)              <- ib_feasibility.json (2 atm swarm transport)
  eta_beam_wallplug      <- ib_feasibility.json (0.07: the PROGRAM-
                             stated ~7% current e-beam wall-plug,
                             docs/DE-AR_targets.txt "GOAL" block)

Chain and bookkeeping:
  wall -> eta_src -> RF at the cell -> eta_c (fraction ABSORBED in the
  gas: 1 - reflection - window/mode losses) -> absorbed joules, which
  earn the measured marginal m per joule.

  Blended intrinsic gain (microwave counted at the RF level):
      d_eta ~= (eta_c * m - eta0) * E_RF / E_tot     [exact to 1st order]
      -> intrinsic gain requires  eta_c > eta0 / m
  Wall-plug gain requires        eta_src * eta_c * m > eta_beam_wp

Reflection estimate: normal-incidence Fresnel on a collisional Drude
half-space, eps = 1 - wp^2/(w(w - i nu)) — a BOTE for the launch
budget (real geometry is the full-wave solver's job, stated).
"""
import json
import sys
from pathlib import Path

import numpy as np

for base in Path(__file__).resolve().parents:
    if (base / "REBUILD" / "MODULES").is_dir():
        SPECK = base
        break
else:
    sys.exit("run from inside SPECK")

OUT = SPECK / "rebuild_out" / "mw_benchmarks" / "derated_operating_point.json"

E_CH = 1.602176634e-19
EPS0 = 8.8541878128e-12
M_E = 9.1093837015e-31
C = 299792458.0

bench = json.loads((SPECK / "rebuild_out/mw_benchmarks/mw_benchmarks_loki.json").read_text())
ib = json.loads((SPECK / "rebuild_out/mw_benchmarks/ib_feasibility.json").read_text())

case = next(c for c in bench["cases"] if c["case"] == "ne_optimized")
ag = case["schemes"]["afterglow"]
m = ag["marginal"]                       # +0.1002 per absorbed joule
eta0 = case["eta0"]                      # 0.06194
e_dep = case["E_dep_J"]                  # 125.49 J
e_mw = ag["E_MW_J"]                      # 0.080 J absorbed in the run
eta_wp = ib["eta_beam_wallplug"]         # 0.07 (program-stated)
nu_lo = ib["nu_m_per_s"]["1.0"]          # ~1.95e12 s^-1 (2 atm, low E/N)


def reflectivity(f_hz, ne_m3, nu):
    w = 2 * np.pi * f_hz
    wp2 = ne_m3 * E_CH ** 2 / (EPS0 * M_E)
    eps = 1.0 - wp2 / (w * (w - 1j * nu))
    n = np.sqrt(eps)
    return float(abs((n - 1) / (n + 1)) ** 2)


bands = [("2.45 GHz", 2.45e9), ("9.5 GHz", 9.5e9), ("95 GHz", 95e9)]
nes = [1e19, 1e20, 1e21]
refl = {lab: {f"{ne:.0e}": reflectivity(f, ne, nu_lo) for ne in nes}
        for lab, f in bands}

eta_c_star = eta0 / m                    # intrinsic-gain threshold
# wall-plug: eta_src * eta_c > eta_wp / m
wp_ratio = eta_wp / m
sources = [("magnetron 2.45 GHz", 0.72), ("X-band 9.5 GHz", 0.60),
           ("gyrotron 95 GHz", 0.35)]
wp_table = {lab: {"eta_c_required": wp_ratio / es,
                  "feasible": wp_ratio / es <= 1.0}
            for lab, es in sources}

# blended-gain table: d_eta(pp) over eta_c and E_RF
eta_cs = [0.5, 0.62, 0.75, 0.85, 1.0]
e_rfs = [0.1, 0.5, 2.0]                   # J of RF at the cell
gain = {}
for ec in eta_cs:
    gain[f"{ec:.2f}"] = {
        f"{er:g} J": 100 * (ec * m - eta0) * er / (e_dep + er)
        for er in e_rfs}

print("[derated operating point] built on the ne_optimized afterglow "
      "channel (the one robust positive)")
print(f"  measured per-ABSORBED-joule marginal m = {100 * m:+.2f}%  |  "
      f"machine eta0 = {100 * eta0:.3f}%  (E_dep {e_dep:.1f} J)")
print(f"  first-order identity check: (m - eta0)*E_MW/E_tot = "
      f"{100 * (m - eta0) * e_mw / (e_dep + e_mw):+.5f} pp vs recorded "
      f"{ag['d_eta_pp']:+.5f} pp")
print(f"\n  INTRINSIC-gain threshold: eta_c > eta0/m = "
      f"{100 * eta_c_star:.0f}% of RF energy absorbed in the gas")
print(f"  WALL-PLUG-gain requirement (bar = program-stated "
      f"{100 * eta_wp:.0f}%): eta_src*eta_c > {wp_ratio:.2f}")
for lab, es in sources:
    t = wp_table[lab]
    print(f"    {lab} (eta_src {100 * es:.0f}%): needs eta_c "
          f">= {100 * t['eta_c_required']:.0f}%"
          + ("" if t["feasible"] else "  -> IMPOSSIBLE (>100%)"))
print(f"  bar sensitivity: at a 5% bar eta_src*eta_c > {0.05 / m:.2f}; "
      f"at 8%: > {0.08 / m:.2f}")

print("\n  Fresnel reflection off the collisional plasma face "
      f"(nu = {nu_lo:.2e} s^-1; normal incidence BOTE):")
print(f"    {'band':>10} " + " ".join(f"n_e {ne:.0e}" for ne in nes))
for lab, _f in bands:
    print(f"    {lab:>10} " + " ".join(
        f"{100 * refl[lab][f'{ne:.0e}']:8.0f}%" for ne in nes))

print("\n  blended intrinsic gain d_eta (pp) = (eta_c*m - eta0)*E_RF/E_tot:")
print(f"    {'eta_c':>6} " + " ".join(f"{er:>9g} J" for er in e_rfs))
for ec in eta_cs:
    row = gain[f"{ec:.2f}"]
    print(f"    {ec:6.2f} " + " ".join(
        f"{row[f'{er:g} J']:+11.5f}" for er in e_rfs))

OUT.write_text(json.dumps({
    "inputs": {"m": m, "eta0": eta0, "E_dep_J": e_dep, "E_MW_J": e_mw,
               "eta_beam_wallplug": eta_wp, "nu_m_per_s": nu_lo},
    "intrinsic_threshold_eta_c": eta_c_star,
    "wallplug_requirement_eta_src_x_eta_c": wp_ratio,
    "wallplug_by_source": wp_table,
    "fresnel_reflectivity": refl,
    "blended_gain_pp": gain}, indent=1))
print(f"\nreport -> {OUT}")
