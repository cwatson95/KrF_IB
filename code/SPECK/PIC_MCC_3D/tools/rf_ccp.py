"""RF capacitively-coupled gap -- the machinery for the Turner benchmark
(gate G5.1), BEGUN, not claimed.

WHAT IS IN: the pieces a CCP needs that this package now has -- a phi(t)
electrode waveform (es_sim ``_bias``, leapfrog-gated in tests/test_rf_bias),
the He LXCat set (ships in MODULES/lxcat/), mobile ions + the M5 gamma
process, absorbing electrodes with per-wall attribution.  ``build_rf_gap``
wires them into a runnable RF discharge; ``--smoke`` runs a short scaled
window to prove the loop closes.

WHAT IS NOT CLAIMED: the benchmark itself.  Turner et al., Phys. Plasmas 20,
013507 (2013), DOI 10.1063/1.4775084, prescribe their OWN reaction set,
surface model and numerics per case, and G5.1 means matching the published
n_i(x) profiles within the inter-code spread over thousands of RF cycles
(hours-days of compute).  Blockers, stated: the paper's exact case spec must
be transcribed first-hand (the numbers below are from a reproduction repo,
see CASES_PROVENANCE), and Coulomb collisions (M6) are open -- the paper
discusses where they matter.

CASES_PROVENANCE: transcribed 2026-08-14 from the lase-unb/ccp-benchmark
reproduction (github.com/lase-unb/ccp-benchmark, src/parameters.cpp), NOT
yet verified against the paper's Table -- do that before claiming G5.1.

Run from the Kinetic_Modeling root:

    python3 SPECK/PIC_MCC_3D/tools/rf_ccp.py --smoke     # ~a minute
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

from PIC_MCC_3D.MODULES import species3d as SP                 # noqa: E402
from PIC_MCC_3D.MODULES.cross_sections import mixture_set      # noqa: E402
from PIC_MCC_3D.MODULES.es_sim import ESSim                    # noqa: E402
from PIC_MCC_3D.MODULES.gas import GasMixture                  # noqa: E402
from PIC_MCC_3D.MODULES.mcc3d import MCC3D                     # noqa: E402
from PIC_MCC_3D.MODULES.units import K_B                       # noqa: E402


def turner_cases() -> dict:
    """The four benchmark cases (see CASES_PROVENANCE in the module doc)."""
    base = {"gas": "He", "gap_m": 6.7e-2, "f_Hz": 13.56e6, "T_gas_K": 300.0}
    return {
        1: dict(base, V0=450.0, n_neutral_m3=9.64e20, n_plasma_m3=2.56e14,
                cells=129, steps_per_period=400, ppc=512),
        2: dict(base, V0=200.0, n_neutral_m3=32.1e20, n_plasma_m3=5.12e14,
                cells=257, steps_per_period=800, ppc=256),
        3: dict(base, V0=150.0, n_neutral_m3=96.4e20, n_plasma_m3=5.12e14,
                cells=513, steps_per_period=1600, ppc=128),
        4: dict(base, V0=120.0, n_neutral_m3=321.0e20, n_plasma_m3=3.84e14,
                cells=513, steps_per_period=3200, ppc=64),
    }


def build_rf_gap(case: int = 1, nz: int = None, ppc: int = None,
                 seed: int = 20260814):
    """A runnable RF He gap at a Turner case's physical point.

    nz/ppc may be scaled down for smoke runs; the returned dict says what
    was actually built so nothing scaled is mistaken for the benchmark.
    """
    c = turner_cases()[case]
    nz = int(nz or c["cells"])
    ppc = int(ppc or c["ppc"])
    L = c["gap_m"]
    h = L / (nz - 1)
    dt = 1.0 / (c["f_Hz"] * c["steps_per_period"])
    # neutral background at the case's density (pressure follows from nkT)
    p_atm = c["n_neutral_m3"] * K_B * c["T_gas_K"] / 101325.0
    mix = GasMixture({"He": 1.0}, pressure_atm=p_atm, T_gas_K=c["T_gas_K"])
    cs = mixture_set(mix)
    e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=4 * ppc * nz)
    he_ion = SP.ParticleSpecies(SP.Species.ion("He", 4.0026),
                                capacity=4 * ppc * nz)
    mcc = MCC3D(cs, mixture=mix, ion_species={"He": he_ion})
    n_macro = ppc * (nz - 1)
    rng = np.random.default_rng(seed)
    vol = (4 * h) * (4 * h) * L
    w = c["n_plasma_m3"] * vol / n_macro
    for pop, T_eV in ((e_pop, 2.0), (he_ion, 0.026)):
        vth = np.sqrt(T_eV * 1.602176634e-19 / pop.species.mass)
        pop.spawn(rng.random(n_macro) * 4 * h, rng.random(n_macro) * 4 * h,
                  rng.random(n_macro) * L,
                  rng.normal(0, vth, n_macro), rng.normal(0, vth, n_macro),
                  rng.normal(0, vth, n_macro), np.full(n_macro, w))
    omega = 2.0 * np.pi * c["f_Hz"]
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop, he_ion]), dt,
                phi_lo=0.0, phi_hi=lambda t: c["V0"] * np.sin(omega * t),
                mcc=mcc, seed=seed)
    return {"sim": sim, "e": e_pop, "ion": he_ion, "case": c,
            "built": {"nz": nz, "ppc": ppc, "scaled": (nz != c["cells"]
                                                      or ppc != c["ppc"])},
            "pressure_atm": p_atm}


def smoke(case: int = 1, cycles: float = 2.0, nz: int = 65, ppc: int = 8,
          seed: int = 20260814) -> dict:
    """Prove the loop closes at the case's physical point (scaled numerics)."""
    b = build_rf_gap(case, nz=nz, ppc=ppc, seed=seed)
    sim, e_pop, ion = b["sim"], b["e"], b["ion"]
    c = b["case"]
    n_steps = int(cycles * c["steps_per_period"])
    ne0 = e_pop.real_count()
    stab = sim.stability(c["n_plasma_m3"], 2.0)
    sim.run(n_steps, record_every=max(n_steps // 20, 1))
    return {"built": b["built"], "steps": n_steps,
            "ne0": ne0, "ne1": e_pop.real_count(),
            "ni1": ion.real_count(),
            "collisions": int(sim.collision_tally.sum())
            if sim.collision_tally is not None else 0,
            "wall_keys": sorted(str(k) for k in sim.wall_flux_w),
            "stability": stab}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--case", type=int, default=1, choices=(1, 2, 3, 4))
    args = ap.parse_args(argv)
    if args.smoke:
        r = smoke(args.case)
        print(f"case {args.case} smoke ({'SCALED' if r['built']['scaled'] else 'full'} "
              f"nz={r['built']['nz']} ppc={r['built']['ppc']}), "
              f"{r['steps']} steps:")
        print(f"  n_e {r['ne0']:.3e} -> {r['ne1']:.3e} real; "
              f"ions {r['ni1']:.3e}; collisions {r['collisions']}")
        print(f"  walls hit: {r['wall_keys']}")
        s = r["stability"]
        print(f"  stability: dx/lambda_D={s['dx_over_debye']:.2f} "
              f"wpe*dt={s['wpe_dt']:.3f} resolved={s['resolved']}")
        print("  (benchmark NOT claimed -- see module docstring)")
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
