"""Child-Langmuir diode -- plan gate G4.1, run as an experiment.

The oldest quantitative test a bounded PIC code can face (Child 1911,
Langmuir 1913): emit MORE electron current than the gap can carry and the
space charge must throttle transmission to

    J_CL = (4/9) eps0 sqrt(2 e / m_e) V^{3/2} / d^2

independent of how much was offered.  This exercises, end to end: the M5
flux emitter, the Dirichlet-z Poisson solve, gather/push, and the absorbing
electrodes with per-wall attribution -- if any of them is wrong, the
transmitted current does not land on the law.

STATED LIMITS of the comparison (all reported in the result dict):

* Emission is a flux-weighted half-Maxwellian at ``T_emit_eV`` (module M5's
  stated distribution); the law assumes cold emission, so finite T biases
  J HIGH by O(sqrt(T/V)) (the Langmuir thermal correction) -- a few percent
  at the defaults (T/V = 4e-4).
* The virtual-cathode layer has the local Debye scale of the emitted cloud
  (~1e-5 m at the defaults).  A uniform grid under-resolves it at small nz;
  ``scan_nz`` exists to SHOW the convergence rather than assert it away.
* Transmitted current is measured from the anode's absorbed-weight series
  (``wall_flux_w``), straight-line fit over the trailing window -- a real
  collector-current measurement, not a formula shortcut.

Run from the Kinetic_Modeling root:

    python3 SPECK/PIC_MCC_3D/tools/diode.py                 # one point + PNG
    python3 SPECK/PIC_MCC_3D/tools/diode.py --sweep         # J(V) vs V^1.5
    python3 SPECK/PIC_MCC_3D/tools/diode.py --scan-nz       # grid convergence
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "PIC_MCC_3D" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("PIC_MCC_3D/MODULES not found above this script")

from PIC_MCC_3D.MODULES import emission3d as EM               # noqa: E402
from PIC_MCC_3D.MODULES import species3d as SP                # noqa: E402
from PIC_MCC_3D.MODULES.es_sim import ESSim                   # noqa: E402
from PIC_MCC_3D.MODULES.units import E_CHARGE, EPS0, M_E      # noqa: E402

OUT_DEFAULT = ROOT / "picmcc_out" / "diode"


def child_langmuir_j(V: float, d_m: float) -> float:
    """The law itself [A/m^2]."""
    return (4.0 / 9.0) * EPS0 * np.sqrt(2.0 * E_CHARGE / M_E) \
        * V ** 1.5 / d_m ** 2


def run_diode(V: float = 250.0, d_m: float = 5.0e-3, nz: int = 129,
              j_factor: float = 3.0, T_emit_eV: float = 0.1,
              transits: float = 10.0, fit_window: float = 0.4,
              n_target: int = 1500, seed: int = 20260814,
              record_every: int = 20) -> dict:
    """One diode at (V, d): emit j_factor x J_CL, measure what crosses."""
    j_cl = child_langmuir_j(V, d_m)
    v_f = np.sqrt(2.0 * E_CHARGE * V / M_E)
    t_cl = 3.0 * d_m / v_f                      # space-charge-limited transit
    h = d_m / (nz - 1)
    vth = np.sqrt(2.0 * E_CHARGE * T_emit_eV / M_E)
    dt = 0.35 * h / (v_f + 4.0 * vth)
    area = (4 * h) * (4 * h)
    # macro weight sized for ~n_target in flight at steady state
    w_macro = max(j_cl * area * t_cl / E_CHARGE / n_target, 1e-12)

    e_pop = SP.ParticleSpecies(SP.Species.electron(),
                               capacity=8 * n_target)
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop]), dt,
                phi_lo=0.0, phi_hi=float(V), seed=seed)
    emitter = EM.FluxEmitter(e_pop, j_factor * j_cl, wall="lo",
                             T_emit_eV=T_emit_eV, w_macro=w_macro)
    sim.injectors.append(emitter)

    n_steps = int(np.ceil(transits * t_cl / dt))
    t_series, w_hi_series, phi_axis = [], [], None
    sim.prepare()
    for k in range(n_steps):
        sim.step()
        if k % record_every == 0:
            t_series.append(sim.time)
            w_hi_series.append(sim.wall_flux_w.get(("e-", "hi"), 0.0))
    phi_axis = sim.phi.mean(axis=(0, 1)).copy()

    t = np.asarray(t_series)
    w_hi = np.asarray(w_hi_series)
    tail = t >= (1.0 - fit_window) * t[-1]
    slope = np.polyfit(t[tail], w_hi[tail], 1)[0]      # real electrons / s
    j_meas = slope * E_CHARGE / area
    return {
        "V": V, "d_m": d_m, "nz": nz, "dt": dt, "j_factor": j_factor,
        "T_emit_eV": T_emit_eV, "transits": transits,
        "j_cl_A_m2": j_cl, "j_meas_A_m2": float(j_meas),
        "ratio": float(j_meas / j_cl),
        "j_emitted_A_m2": j_factor * j_cl,
        "phi_axis_V": phi_axis, "z_m": np.arange(nz) * h,
        "phi_min_V": float(phi_axis.min()),
        "phi_min_z_m": float(np.argmin(phi_axis) * h),
        "t_s": t, "w_hi": w_hi, "n_macro_final": e_pop.n,
        "thermal_bias_note": "finite T_emit biases j HIGH by O(sqrt(T/V)) "
                             f"~ {2.66 * np.sqrt(T_emit_eV / V):.1%}",
    }


def jv_sweep(voltages=(100.0, 175.0, 250.0, 400.0), **kw) -> list:
    return [run_diode(V=v, **kw) for v in voltages]


def scan_nz(nzs=(33, 65, 129, 257), **kw) -> list:
    return [run_diode(nz=n, **kw) for n in nzs]


# ----------------------------------------------------------------- artifacts
def render(res_one: dict, sweep: list, scan: list, out_dir: Path) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    made = []

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    ax = axes[0]
    ax.plot(res_one["z_m"] * 1e3, res_one["phi_axis_V"], lw=1.6)
    ax.axhline(0.0, color="0.7", lw=0.7)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("phi [V]")
    ax.set_title(f"potential; min {res_one['phi_min_V']:.2f} V @ "
                 f"{res_one['phi_min_z_m']*1e6:.0f} um (virtual cathode)")
    ax = axes[1]
    area_note = res_one["ratio"]
    ax.plot(res_one["t_s"] * 1e9, res_one["w_hi"] * 1e-6, lw=1.4)
    ax.set_xlabel("t [ns]")
    ax.set_ylabel("anode-collected real e- [x1e6]")
    ax.set_title(f"collector series; fit ratio J/J_CL = {area_note:.3f}")
    ax = axes[2]
    if sweep:
        vs = np.array([r["V"] for r in sweep])
        js = np.array([r["j_meas_A_m2"] for r in sweep])
        vv = np.linspace(vs.min() * 0.9, vs.max() * 1.05, 100)
        ax.plot(vv, child_langmuir_j(vv, sweep[0]["d_m"]), "-", lw=1.2,
                label="(4/9)eps0 sqrt(2e/m) V^1.5/d^2")
        ax.plot(vs, js, "o", ms=5, label="PIC transmitted")
        ax.set_xlabel("V [V]")
        ax.set_ylabel("J [A/m^2]")
        ax.legend(fontsize=8)
        ax.set_title("gate G4.1: the V^{3/2} law")
    fig.suptitle(f"Child-Langmuir diode, d={res_one['d_m']*1e3:.1f} mm, "
                 f"emit {res_one['j_factor']:.0f}x J_CL at "
                 f"T={res_one['T_emit_eV']} eV")
    fig.tight_layout()
    p = out_dir / "child_langmuir.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    made.append(p)

    if scan:
        fig, ax = plt.subplots(figsize=(5.2, 3.8))
        ns = np.array([r["nz"] for r in scan])
        rs = np.array([r["ratio"] for r in scan])
        ax.axhline(1.0, color="0.6", lw=0.8)
        ax.plot(ns, rs, "o-", ms=5)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("nz (uniform grid)")
        ax.set_ylabel("J / J_CL")
        ax.set_title("grid convergence (virtual-cathode layer ~1e-5 m)")
        fig.tight_layout()
        p = out_dir / "convergence_nz.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        made.append(p)

    np.savez_compressed(
        out_dir / "diode_data.npz",
        z_m=res_one["z_m"], phi_axis_V=res_one["phi_axis_V"],
        t_s=res_one["t_s"], w_hi=res_one["w_hi"],
        sweep_V=np.array([r["V"] for r in sweep] if sweep else []),
        sweep_j=np.array([r["j_meas_A_m2"] for r in sweep] if sweep else []),
        scan_nz=np.array([r["nz"] for r in scan] if scan else []),
        scan_ratio=np.array([r["ratio"] for r in scan] if scan else []),
        meta=np.array([f"V={res_one['V']} d={res_one['d_m']} "
                       f"ratio={res_one['ratio']:.4f}"]))
    made.append(out_dir / "diode_data.npz")
    return made


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--V", type=float, default=250.0)
    ap.add_argument("--d-mm", type=float, default=5.0)
    ap.add_argument("--nz", type=int, default=129)
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--scan-nz", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    res = run_diode(V=args.V, d_m=args.d_mm * 1e-3, nz=args.nz)
    print(f"J_CL      {res['j_cl_A_m2']:9.2f} A/m^2")
    print(f"J_meas    {res['j_meas_A_m2']:9.2f} A/m^2   "
          f"ratio {res['ratio']:.3f}  ({res['thermal_bias_note']})")
    print(f"phi_min   {res['phi_min_V']:.3f} V at {res['phi_min_z_m']*1e6:.0f} um")
    sweep = jv_sweep(d_m=args.d_mm * 1e-3, nz=args.nz) if args.sweep else []
    if sweep:
        for r in sweep:
            print(f"  V={r['V']:6.1f}  J/J_CL = {r['ratio']:.3f}")
    scan = scan_nz(V=args.V, d_m=args.d_mm * 1e-3) if args.scan_nz else []
    if scan:
        for r in scan:
            print(f"  nz={r['nz']:4d}  J/J_CL = {r['ratio']:.3f}")
    made = render(res, sweep, scan, args.out)
    for p in made:
        print("wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
