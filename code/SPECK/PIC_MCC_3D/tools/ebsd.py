"""E-beam-sustained conduction -- the non-self-sustaining i-V characteristic
(C5, kinetic skeleton).

Raizer Fig. 7.1 (read first-hand 2026-08-14): a gap fed by an EXTERNAL
ionization source and biased below breakdown shows (i) a saturation current
when every beam-made pair is collected, then (ii) exponential amplification
i = i_sat exp(alpha d) as the bias climbs -- but it stays non-self-sustaining:
remove the source and the current dies.  That stability is the whole point of
e-beam-sustained discharges (EBSD), and it is reproduced here structurally:
the gamma process is NOT registered, so self-sustainment is impossible by
construction, exactly Raizer's "without secondary emission a discharge is not
self-sustaining".

MACHINERY: beam3d.BeamInjector (kinetic pair mode -- secondaries sampled from
the Opal-Peterson spectrum, cascaded by the MCC; W emerges, it is not an
input) drives a biased Ar/Kr/F2 gap along a periodic axis (uniform volumetric
source); the anode's absorbed-weight series gives the collected current.

STATED LIMITS: ions stay as immobile grid charge (no ion transport, no ion
MCC -- plan section 2.1 items 2-3 for this geometry); windows are kept short
enough that the accumulated ion space charge stays small against the applied
field, and the report prints that ratio.  The full coupled EBSD (chemistry
feedback via SPECK/tools/couple0d) is the named next step.

Run from the Kinetic_Modeling root:

    python3 SPECK/PIC_MCC_3D/tools/ebsd.py --scan     # i-V curve + PNG
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

from PIC_MCC_3D.MODULES import beam3d as BM                    # noqa: E402
from PIC_MCC_3D.MODULES import species3d as SP                 # noqa: E402
from PIC_MCC_3D.MODULES.cross_sections import mixture_set      # noqa: E402
from PIC_MCC_3D.MODULES.es_sim import ESSim                    # noqa: E402
from PIC_MCC_3D.MODULES.gas import GasMixture                  # noqa: E402
from PIC_MCC_3D.MODULES.mcc3d import MCC3D                     # noqa: E402
from PIC_MCC_3D.MODULES.units import E_CHARGE, M_E             # noqa: E402

OUT_DEFAULT = ROOT / "picmcc_out" / "ebsd"


def run_point(V: float, p_torr: float = 38.0, d_m: float = 0.8e-3,
              nz: int = 17, beam_kV: float = 10.0, beam_j_A_cm2: float = 1e-5,
              window_ns: float = 16.0, fit_window: float = 0.5,
              fractions: dict = None, seed: int = 20260814) -> dict:
    """One bias point of the beam-sustained gap."""
    fractions = fractions or {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
    mix = GasMixture(dict(fractions), pressure_atm=p_torr / 760.0,
                     T_gas_K=300.0)
    cs = mixture_set(mix)
    e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=65536)
    mcc = MCC3D(cs, mixture=mix)
    h = d_m / (nz - 1)
    v_cap = np.sqrt(2.0 * E_CHARGE * max(V, 100.0) * 1.5 / M_E)
    dt = 0.35 * h / v_cap
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop]), dt,
                phi_lo=0.0, phi_hi=float(V), mcc=mcc, seed=seed)
    spec = BM.BeamSpec(voltage_kV=beam_kV, j_A_cm2=beam_j_A_cm2,
                       pulse_ns=1.0e5, rise_ns=0.0, foil_transmission=1.0,
                       axis=0, from_low=True)
    inj = BM.BeamInjector(mix, spec, e_pop, macro_energy_eV=2.0e3,
                          rng=np.random.default_rng(seed + 3))
    sim.injectors.append(inj)
    n_steps = int(window_ns * 1e-9 / dt)
    area = sim.L[0] * sim.L[1]
    t_series, w_hi = [], []
    sim.prepare()
    for k in range(n_steps):
        sim.step()
        if k % 50 == 0:
            t_series.append(sim.time)
            w_hi.append(sim.wall_flux_w.get(("e-", "hi"), 0.0))
    t = np.asarray(t_series)
    w = np.asarray(w_hi)
    tail = t >= (1.0 - fit_window) * t[-1]
    slope = np.polyfit(t[tail], w[tail], 1)[0] if tail.sum() > 2 else 0.0
    j_anode = slope * E_CHARGE / area
    # saturation reference: every beam-made pair collected, no multiplication
    pair_rate = inj.pairs_created                       # kinetic mode: 0
    # in kinetic mode pairs EMERGE via MCC ionizations; count them
    ion_k = [k2 for k2 in range(cs.n_proc) if cs.kinds[k2] == "ionization"]
    n_ion = float(sim.collision_tally[ion_k].sum()) if \
        sim.collision_tally is not None else 0.0
    # ion space charge vs applied field (smallness check): the immobile ion
    # charge lives on sim.rho_heavy [C/m^3]; sheet-field bound sigma/2eps0
    cell_vol = sim.h ** 3
    q_ion_m2 = float(np.sum(np.maximum(sim.rho_heavy, 0.0))) * cell_vol / area
    e_applied = V / d_m
    e_spacecharge = q_ion_m2 / 8.854e-12 / 2.0
    return {"V": V, "j_anode_A_m2": float(j_anode),
            "beam_in_J": inj.energy_in_J,
            "beam_dep_J": inj.energy_deposited_J,
            "n_ionizations": n_ion, "n_macro": e_pop.n,
            "spacecharge_over_applied": float(e_spacecharge
                                              / max(e_applied, 1.0)),
            "dt": dt, "steps": n_steps, "d_m": d_m, "p_torr": p_torr}


def scan_bias(voltages=(20.0, 60.0, 120.0, 180.0, 240.0), **kw) -> list:
    rows = []
    for v in voltages:
        r = run_point(v, **kw)
        rows.append(r)
        print(f"  V={v:6.1f}  J_anode={r['j_anode_A_m2']:.4e} A/m^2  "
              f"ionizations={r['n_ionizations']:.0f}  "
              f"dSC/E={r['spacecharge_over_applied']:.2e}")
    return rows


def render(rows: list, out_dir: Path) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    v = np.array([r["V"] for r in rows])
    j = np.array([r["j_anode_A_m2"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(v, j, "o-", ms=6)
    ax.set_xlabel("bias V [V]")
    ax.set_ylabel("anode current density [A/m^2]")
    ax.set_title("beam-sustained gap i-V (Raizer Fig. 7.1 regime):\n"
                 "saturation then exp(alpha d) amplification, "
                 "never self-sustaining (no gamma)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    p = out_dir / "ebsd_iv.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    np.savez_compressed(out_dir / "ebsd_data.npz", V=v, j_anode=j)
    return [p, out_dir / "ebsd_data.npz"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)
    if not args.scan:
        ap.print_help()
        return 0
    rows = scan_bias()
    for p in render(rows, args.out):
        print("wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
