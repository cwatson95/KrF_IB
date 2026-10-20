#!/usr/bin/env python3
"""REBUILD CLI — presets, run report, scorecard, ASE, tier comparison.

Usage (from the Kinetic_Modeling root):

    python3 SPECK/REBUILD/run.py --preset krf_beam_0d
    python3 SPECK/REBUILD/run.py --preset krf_beam_3d --ase
    python3 SPECK/REBUILD/run.py --preset arf_chen_0d --scorecard
    python3 SPECK/REBUILD/run.py --tier-compare
    python3 SPECK/REBUILD/run.py --list

Every preset prints its regime line (grid, reactions, drive route, ledger
residual, wall time) — the no-silent-scope rule.
"""

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))

from REBUILD.MODULES.engine import Device                      # noqa: E402
from REBUILD.MODULES.ledger import format_budget               # noqa: E402
from REBUILD.MODULES.params import engine_parameters           # noqa: E402

PRESETS = {
    "krf_beam_0d": {
        "doc": "DE-AR class KrF cell, e-beam pumped, Ar-dominant, 0-D "
               "column (the v2-class comparison point)",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                   "t_end_ns": 200.0, "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                   "beam_axis": 0},
    },
    "krf_ne_0d": {
        "doc": "Ne-buffered KrF fill (the spectre brief), 0-D",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                   "t_end_ns": 200.0, "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                   "beam_axis": 0},
        "mixture": {"Ne": 0.949, "Kr": 0.05, "F2": 0.001},
    },
    "krf_beam_3d": {
        "doc": "3-D device grid, beam along x with a real deposition "
               "gradient (10 cm path at 3 atm), 240 cells — NB the "
               "gradient is x1.02: an Electra-class beam barely "
               "attenuates, so the plasma is honestly near-uniform",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                   "t_end_ns": 200.0, "nx": 8, "ny": 3, "nz": 10,
                   "Lx_m": 0.10, "Ly_m": 0.0141, "Lz_m": 0.75,
                   "beam_axis": 0},
    },
    "krf_beam_3d_bragg": {
        "doc": "3-D grid with the beam STOPPED inside the cell: 100 kV "
               "into 6 cm at 3 atm (the recorded mesh-refinement Bragg "
               "case) — deposition rises x2.9 to the Bragg peak at "
               "x~4 cm then collapses to ZERO; the spatial-structure "
               "demo, 360 cells",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                   "beam_kV": 100.0, "t_end_ns": 200.0, "nx": 12,
                   "ny": 3, "nz": 10, "Lx_m": 0.06, "Ly_m": 0.0141,
                   "Lz_m": 0.75, "beam_axis": 0},
    },
    "krf_arpa_mw": {
        "doc": "ARPA-E chamber design study (geometry measured into "
               "freecad/arpa_chamber_ports.json): e-beam KrF pulse, then a "
               "WR90 microwave reheat window in the afterglow. Pumped "
               "column = the O190.5 front-aperture beam (area-equivalent "
               "168.8 mm square) through the 264.16 mm chamber depth "
               "(y axis); optical axis z = the chamber top/bottom window "
               "axis, internal mirrors R_OC 0.9 / L_op 0.02 assumed "
               "(research-cell optics; the physical O76/O32 windows would "
               "vignette an external cavity). Vessel-safe fill 2 atm "
               "(1 atm gauge). MW field = the three WR90 feeds' overmoded "
               "cavity taken spatially uniform; the implied absorbed power "
               "is an output — see tools/mw_reheat_experiment.py for the "
               "baseline-vs-reheat marginal-efficiency comparison. Feed "
               "budget: 3x WR90 pulsed magnetron class ~1 MW injected, "
               "field ceiling 30 Td (empty-cavity buildup at Q ~ 1e4).",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 2.0,
                   "t_end_ns": 400.0, "beam_kV": 500.0,
                   "beam_j_A_cm2": 140.0, "beam_flat_ns": 100.0,
                   "beam_rise_ns": 5.0, "beam_axis": 1,
                   "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.1688, "Ly_m": 0.26416, "Lz_m": 0.1688,
                   "r_oc": 0.9, "l_op": 0.02,
                   "mw_power_W": 1.0e6, "mw_en_max_td": 30.0,
                   "mw_t0_ns": 130.0,
                   "mw_flat_ns": 150.0, "mw_rise_ns": 10.0},
    },
    "krf_arpa_mw_3d": {
        "doc": "the ARPA-E chamber study on the 3-D device grid (4x6x4 = "
               "96 cells, beam axis y resolved) — the glow-movie / spatial "
               "showcase of krf_arpa_mw; same physics point",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 2.0,
                   "t_end_ns": 400.0, "beam_kV": 500.0,
                   "beam_j_A_cm2": 140.0, "beam_flat_ns": 100.0,
                   "beam_rise_ns": 5.0, "beam_axis": 1,
                   "nx": 4, "ny": 6, "nz": 4,
                   "Lx_m": 0.1688, "Ly_m": 0.26416, "Lz_m": 0.1688,
                   "r_oc": 0.9, "l_op": 0.02,
                   "mw_power_W": 1.0e6, "mw_en_max_td": 30.0,
                   "mw_t0_ns": 130.0,
                   "mw_flat_ns": 150.0, "mw_rise_ns": 10.0},
    },
    "krf_hyper_mw": {
        "doc": "HYPER-LASER hybrid-pumping frame (docs/DE-AR_targets.txt, "
               "DE-AR0002030): 600 kV beam into the 10 cm pumped cube at "
               "2 atm Ar/Kr/F2, 95 GHz gyrotron booster up to 250 kW "
               "coupled. Default window = afterglow reheat (t0 80 ns; the "
               "beam ends at 70). en_max 30 Td is the ~x5-field-buildup "
               "case (traveling-wave 250 kW over the aperture is ~6 Td at "
               "2 atm) and sits AT the fill's self-sustainment crossing "
               "(~28 Td, ionization = attachment in the swarm tables) — "
               "bracket it with --en-max. Internal optics R_OC 0.95 / "
               "L_op 0.02 stated. The experiment's M2.6 observable is "
               "248.5 nm FLUORESCENCE (prop. to small-signal gain), not "
               "extracted energy — tools/hybrid_timing_study.py reports "
               "BOTH metrics for every timing scheme.",
        "params": {"mode": "krf", "drive": "beam", "pressure_atm": 2.0,
                   "t_end_ns": 400.0, "beam_kV": 600.0,
                   "beam_j_A_cm2": 100.0, "beam_flat_ns": 60.0,
                   "beam_rise_ns": 5.0, "beam_axis": 0,
                   "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.1, "Ly_m": 0.1, "Lz_m": 0.1,
                   "r_oc": 0.95, "l_op": 0.02,
                   "mw_power_W": 2.5e5, "mw_en_max_td": 30.0,
                   "mw_freq_GHz": 95.0,
                   "mw_t0_ns": 80.0, "mw_flat_ns": 200.0,
                   "mw_rise_ns": 10.0},
    },
    "arf_chen_circuit": {
        "doc": "Chen 2025 base case through the REAL circuit (Sec II A "
               "constants + CHEN_ERRATA since 0.1.25; E/N-direct rates, "
               "swarm-drift conduction) — self-consistent breakdown; see "
               "CHANGELOG [0.1.25] for the corrected anchor verdict",
        "params": {"mode": "arf_chen", "drive": "discharge",
                   "discharge_mode": "circuit", "pressure_atm": 3.4,
                   "t_end_ns": 200.0, "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.012, "Ly_m": 0.010, "Lz_m": 0.90},
    },
    "arf_chen_0d": {
        "doc": "Chen 2025 ArF scoring mode (Ne/Ar/F2/Xe, discharge-pumped; "
               "prescribed power waveform stands in for their circuit — an "
               "anchored x3-7 formation-efficiency gap, CHANGELOG 0.1.3)",
        "params": {"mode": "arf_chen", "drive": "discharge",
                   "pressure_atm": 3.4, "t_end_ns": 200.0,
                   "dis_peak_W_m3": 2e12, "dis_t0_ns": 100.0,
                   "dis_fwhm_ns": 25.0,
                   # Chen Table II geometry done RIGHT (found 2026-08-04):
                   # V = 90 cm2 discharge area x 12 mm gap = 108 cm3;
                   # laser axis l_g = 90 cm -> aperture 1.2 cm2. The old
                   # preset read 90 cm2 as the BEAM cross-section (x75
                   # volume error).
                   "nx": 1, "ny": 1, "nz": 1,
                   "Lx_m": 0.012, "Ly_m": 0.010, "Lz_m": 0.90},
    },
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", choices=sorted(PRESETS), default=None)
    ap.add_argument("--list", action="store_true",
                    help="list presets and exit")
    ap.add_argument("--scorecard", action="store_true",
                    help="print the anchors scorecard (S7)")
    ap.add_argument("--ase", action="store_true",
                    help="ASE/parasitics report on the run's gain field")
    ap.add_argument("--tier-compare", action="store_true",
                    help="KINETIC vs LOKI comparison (S4) and exit")
    ap.add_argument("--est-sweep", action="store_true",
                    help="run the S5 EST x1/3..x3 bracket sweep on the "
                         "Ar and Ne fills (~3 min) and exit")
    ap.add_argument("--tier", choices=("loki", "kinetic", "maxwell"), default=None,
                    help="electron-impact rate tier override")
    ap.add_argument("--t-end", type=float, default=None, metavar="NS")
    ap.add_argument("--n-save", type=int, default=121)
    ap.add_argument("--checkpoint", metavar="PATH", default=None,
                    help="write a checkpoint of the final state")
    args = ap.parse_args(argv)

    if args.list:
        for name, spec in sorted(PRESETS.items()):
            print(f"{name:<14s} {spec['doc']}")
        return 0
    if args.tier_compare:
        from REBUILD.MODULES.kinetic_tier import (KineticRates,
                                                  compare_tiers,
                                                  format_comparison)
        from REBUILD.MODULES.loki_tier import rates_for_mode
        kin = KineticRates.load_default()
        lok = rates_for_mode("krf", {"Ar": 0.949, "Kr": 0.05, "F2": 0.001})
        print(format_comparison(compare_tiers(kin, lok)))
        return 0
    if args.est_sweep:
        from REBUILD.MODULES.scorecard import (est_bracket_rows,
                                               format_est_brackets)
        for fill in ({"Ar": 0.949, "Kr": 0.05, "F2": 0.001},
                     {"Ne": 0.949, "Kr": 0.05, "F2": 0.001}):
            print(format_est_brackets(est_bracket_rows(fill)))
            print()
        return 0
    if args.preset is None and args.scorecard:
        from REBUILD.MODULES.scorecard import run_scorecard
        print(run_scorecard(None))
        return 0
    if args.preset is None:
        ap.print_help()
        return 2

    spec = PRESETS[args.preset]
    P = engine_parameters()
    P.set_many(spec["params"], source="problem")
    if args.tier:
        P.set("electron_tier", args.tier, source="cli")
    if args.t_end is not None:
        P.set("t_end_ns", args.t_end, source="cli")
    dev = Device(P, mixture=spec.get("mixture"))
    print(f"[{args.preset}] {spec['doc']}")
    res = dev.run(n_save=args.n_save)
    print(res["regime"])
    print()
    print(format_budget(res["closure"][-1]))
    if args.checkpoint:
        from REBUILD.MODULES.checkpoint import save
        save(args.checkpoint, dev, res["t_ns"][-1] * 1e-9, res["y_final"])
        print(f"checkpoint -> {args.checkpoint}")
    if args.ase:
        from REBUILD.MODULES.ase import ase_report
        rep = ase_report(res)
        print()
        print(f"ASE ({rep['state']}): dlam_eff {rep['dlam_eff_nm']:.2f} nm, "
              f"g0 max {rep['g0_max_m1']:.3f} /m, scalar flux "
              f"{rep['ase_scalar_W_m2']:.3e} W/m2, window "
              f"{rep['ase_window_W_m2']:.3e} W/m2")
        m = rep["parasitic_margins"][rep["worst_parasitic"]]
        print(f"worst parasitic: {rep['worst_parasitic']} — round-trip gain "
              f"{m['round_trip_gain']:.3e} "
              f"({'OSCILLATES' if m['oscillates'] else 'below threshold'})")
    if args.scorecard:
        from REBUILD.MODULES.scorecard import run_scorecard
        print()
        print(run_scorecard(res))
    return 0 if res["success"] else 1


if __name__ == "__main__":
    sys.exit(main())
