"""couple0d -- the PIC <-> chemistry coupler, BEGUN (C4): the data contract
plus a working 0-D sub-cycling core, gated in the reproduce-limit.

THE CONTRACT (the reason SPECK/tools exists: tools may import BOTH halves;
the engines never import each other):

    PIC   -> chem : clamps per macro-step -- n_e [m^-3], T_e [eV], T_g [K],
                    E/N [Td] (drives the LOKI swarm-rate tier), n_ph [m^-3]
    chem  -> PIC  : the heavy-species inventory (S,) after the macro-step.
                    Feeding it back into the MCC's neutral background is the
                    NEXT step -- gas.GasMixture is a static background by
                    design today, and that limit is stated, not hidden.

SUB-CYCLING: the chemistry inside one PIC macro-step is stiff (ion-ion
capture at ~1e21 m^-3 turns over in ps), so the clamped compiled network is
advanced with scipy BDF per macro-step, clamps linearly interpolated across
the step -- "sub-cycling" here means an implicit solve, not tiny explicit
steps.

THE GATE THAT MAKES THIS REAL (tests/test_couple0d.py): the reproduce-limit.
Run the REBUILD engine on its own 0-D KrF problem, then re-advance the heavy
species with THIS coupler, clamping (n_e, T_e, T_g, n_ph) to the engine's own
recorded trajectory across an afterglow window (drive strictly zero, checked).
In that limit the coupler and the engine integrate the same physics, so every
major species must land on the engine's answer.  No tuning knob exists here.
MEASURED (2026-08-14): parity error is set by the CLAMP RESOLUTION, not by
drift -- 2 ns clamp spacing left the worst species (F-) at 16%, 1 ns at 4%.
That number is the macro-step guidance for a real PIC loop: clamps must
resolve the fastest clamp-dependent chemistry (here the attachment-fed F-
inventory), or the handoff error is the interpolation, not the physics.

THE DEMO (--demo): a real handoff at a point inside BOTH envelopes: 38 Torr
is REBUILD's validated pressure floor AND PIC-tractable (3 atm is behind
the lambda_D wall, PIC_MCC_3D.txt section 3).  The Ar/Kr/F2 drift gap
reaches quasi-steady state under a seed current; its measured (n_e, T_e,
E/N) clamp a chemistry window and the inventory response is reported.
One-way today; the return leg is named above.

Run from the Kinetic_Modeling root:

    python3 SPECK/tools/couple0d.py --demo
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir() \
            and (base / "PIC_MCC_3D" / "MODULES").is_dir():
        SPECK = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("SPECK root (REBUILD + PIC_MCC_3D) not found above this script")


# ---------------------------------------------------------- chemistry half
class ClampedChemistry:
    """REBUILD's compiled network, advanced under external electron clamps."""

    def __init__(self, mode: str = "krf", pressure_atm: float = 3.0,
                 mixture: dict = None, device=None):
        from REBUILD.MODULES.engine import Device
        from REBUILD.MODULES.params import engine_parameters
        if device is not None:
            self.dev = device
        else:
            p = engine_parameters()
            p.set_many({"mode": mode, "drive": "beam",
                        "pressure_atm": float(pressure_atm),
                        "nx": 1, "ny": 1, "nz": 1}, source="api")
            self.dev = Device(p, mixture=dict(mixture) if mixture else None)
        self.lay = self.dev.lay
        self.compiled = self.dev.compiled
        self.n_tot = self.dev.n_tot

    def initial_dens(self) -> np.ndarray:
        """The neutral fill, shape (S, 1)."""
        y0 = self.dev.initial_state()
        cells = self.lay.cells_view(y0)
        return cells[: self.lay.n_species].copy()

    def advance(self, dens, clamp0: dict, clamp1: dict, dt_s: float,
                rtol: float = 1e-6) -> dict:
        """One macro-step: dens (S,1) -> (S,1) under interpolated clamps.

        clamp keys: ne [m^-3], te [eV], tg [K], n_ph [m^-3] (optional,
        default 0), en_td [Td] (optional, default None -> the tier's own
        operating point).
        """
        from scipy.integrate import solve_ivp
        for c in (clamp0, clamp1):
            for k in ("ne", "te"):
                if k not in c:
                    raise TypeError(f"clamp missing required key {k!r} -- "
                                    "the contract needs ne and te")
        S = self.lay.n_species
        x0 = np.asarray(dens, dtype=float).reshape(S)

        def mixl(key, w, default=None):
            a = clamp0.get(key, default)
            b = clamp1.get(key, default)
            if a is None or b is None:
                return None
            return (1.0 - w) * float(a) + w * float(b)

        dev = self.dev

        def f(t, x):
            w = t / dt_s if dt_s > 0 else 0.0
            d = np.maximum(x.reshape(S, 1), 0.0)
            ase = dev._ase_rate(d) if getattr(dev, "_ase_on", False) else 0.0
            ch = self.compiled.rhs(
                d,
                np.array([mixl("te", w)]),
                np.array([mixl("tg", w, 300.0)]),
                np.array([mixl("ne", w)]),
                float(mixl("n_ph", w, 0.0)),
                self.n_tot,
                en_td=mixl("en_td", w),
                ase_rate=ase)
            return ch["ddens"].ravel()

        sol = solve_ivp(f, (0.0, float(dt_s)), x0, method="BDF",
                        rtol=rtol, atol=1.0e4)
        if not sol.success:
            raise RuntimeError(f"clamped BDF failed: {sol.message}")
        return {"dens": sol.y[:, -1].reshape(S, 1),
                "nfev": int(sol.nfev), "n_steps": len(sol.t) - 1}


# ---------------------------------------------------------------- PIC half
def pic_steady_gap(p_torr: float = 38.0,
                   fractions: dict = None,
                   V: float = 215.0, d_m: float = 0.8e-3, nz: int = 17,
                   j_seed_A_m2: float = 0.05, transits: float = 2.5,
                   seed: int = 20260814) -> dict:
    """An Ar/Kr/F2 drift gap at REBUILD's pressure floor, run to steady state.

    Returns the clamps the chemistry half consumes: bulk-averaged n_e and
    T_e, the applied E/N, and the machinery report (multiplication, walls).
    """
    from PIC_MCC_3D.MODULES import emission3d as EM
    from PIC_MCC_3D.MODULES import species3d as SP
    from PIC_MCC_3D.MODULES.cross_sections import mixture_set
    from PIC_MCC_3D.MODULES.es_sim import ESSim
    from PIC_MCC_3D.MODULES.gas import GasMixture
    from PIC_MCC_3D.MODULES.mcc3d import MCC3D
    from PIC_MCC_3D.MODULES.units import E_CHARGE, M_E

    fractions = fractions or {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
    mix = GasMixture(dict(fractions), pressure_atm=p_torr / 760.0,
                     T_gas_K=300.0)
    cs = mixture_set(mix)
    e_pop = SP.ParticleSpecies(SP.Species.electron(), capacity=65536)
    mcc = MCC3D(cs, mixture=mix)
    h = d_m / (nz - 1)
    v_cap = np.sqrt(2.0 * E_CHARGE * V * 1.5 / M_E)
    dt = 0.35 * h / v_cap
    sim = ESSim(4, 4, nz, h, SP.SpeciesSet([e_pop]), dt,
                phi_lo=0.0, phi_hi=float(V), mcc=mcc, seed=seed)
    emitter = EM.FluxEmitter(e_pop, j_seed_A_m2, wall="lo", T_emit_eV=2.0,
                             w_macro=max(j_seed_A_m2 * (4 * h) ** 2 * dt
                                         / E_CHARGE / 0.15, 1e-12))
    sim.injectors.append(emitter)
    # transit at the drift speed is not known a priori; run a fixed window
    # of "field transits" (ballistic bound) times a generous factor
    t_transit_ballistic = 2.0 * d_m / np.sqrt(2.0 * E_CHARGE * V / M_E)
    n_steps = int(transits * 40.0 * t_transit_ballistic / dt)
    sim.prepare()
    n_hist = []
    for k in range(n_steps):
        sim.step()
        if k % 100 == 0:
            n_hist.append(e_pop.n)
    vol = sim.L[0] * sim.L[1] * sim.L[2]
    eps = e_pop.energy_eV()
    w = e_pop.w[: e_pop.n]
    te = (2.0 / 3.0) * float(np.sum(w * eps) / max(np.sum(w), 1e-300))
    ne = e_pop.real_count() / vol
    en_td = (V / d_m) / mix.total_density_m3 / 1.0e-21
    ion_k = [k2 for k2 in range(cs.n_proc) if cs.kinds[k2] == "ionization"]
    return {"ne": ne, "te": te, "tg": 300.0, "en_td": en_td,
            "p_torr": p_torr, "V": V, "d_m": d_m,
            "n_macro": e_pop.n, "n_hist": np.array(n_hist),
            "ionizations": float(sim.collision_tally[ion_k].sum())
            if sim.collision_tally is not None else 0.0,
            "quasi_steady": bool(len(n_hist) > 10 and
                                 abs(np.mean(n_hist[-5:])
                                     - np.mean(n_hist[-10:-5]))
                                 < 0.35 * max(np.mean(n_hist[-5:]), 1.0))}


# -------------------------------------------------------------------- demo
def run_demo(window_ns: float = 200.0, macro_steps: int = 8) -> dict:
    """PIC gap -> clamps -> chemistry window.  One-way handoff, stated."""
    pic = pic_steady_gap()
    chem = ClampedChemistry(pressure_atm=pic["p_torr"] / 760.0)
    dens = chem.initial_dens()
    clamp = {"ne": pic["ne"], "te": pic["te"], "tg": pic["tg"],
             "en_td": pic["en_td"], "n_ph": 0.0}
    dt = window_ns * 1e-9 / macro_steps
    for _ in range(macro_steps):
        dens = chem.advance(dens, clamp, clamp, dt)["dens"]
    d0 = chem.initial_dens()
    movers = []
    for s, i in chem.lay.index.items():
        a, b = float(d0[i, 0]), float(dens[i, 0])
        if max(a, b) > 1e10 and abs(b - a) > 0.02 * max(a, 1e10):
            movers.append((s, a, b))
    movers.sort(key=lambda r: -abs(r[2] - r[1]))
    return {"pic": pic, "movers": movers, "dens": dens, "chem": chem}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)
    if not args.demo:
        ap.print_help()
        return 0
    r = run_demo()
    pic = r["pic"]
    print(f"PIC gap ({pic['p_torr']} Torr, {pic['V']} V, "
          f"{pic['d_m']*1e3:.1f} mm): n_e={pic['ne']:.3e} m^-3, "
          f"T_e={pic['te']:.2f} eV, E/N={pic['en_td']:.0f} Td, "
          f"quasi_steady={pic['quasi_steady']}, "
          f"ionizations={pic['ionizations']:.0f}")
    print("chemistry window movers (species, before, after):")
    for s, a, b in r["movers"][:12]:
        print(f"  {s:10s} {a:10.3e} -> {b:10.3e}")
    print("(one-way handoff; the chem->PIC gas-update leg is the named "
          "next step)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
