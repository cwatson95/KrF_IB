#!/usr/bin/env python3
"""Export PIC_MCC_3D swarm rate coefficients for REBUILD's KINETIC tier.

THIS SCRIPT IMPORTS PIC_MCC_3D AND THE ENGINE MUST NOT: the coupling between
the packages is this file's OUTPUT (a .npz with provenance), never an import
(ORESTES_CLASS_REBUILD.txt §5; REBUILD's architecture gate enforces it).

Method — a uniform-field swarm driver over PIC_MCC_3D's gated operators:
electron macro-particles accelerated by a constant E field between
null-collision MCC calls (``MCC3D.collide_adaptive``, the ARPAE/pic1d
oracle-matched operator), immobile-ion path, population held in a band by
weight-preserving down/duplication. After a burn-in the per-process event
tallies over the measuring window give

    k_p = events_p / (N_macro_avg * T_window * n_target_p)   [m^3/s]

and the mean energy <eps> is the time-averaged population mean. Buckets
follow the engine's convention (loki_tier._bucket / beam lo-hi splits).

Statistical honesty: per-bucket relative error ~ 1/sqrt(counts) is exported
alongside the rates; buckets with no events export NaN, and the reader
skips them (the LOKI tier remains the fallback for those channels).

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/export_kinetic.py [--quick]
"""

import json
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "PIC_MCC_3D" / "MODULES").is_dir():
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("PIC_MCC_3D not found above this script")

from PIC_MCC_3D.MODULES.beam3d import (cascade_yields,             # noqa: E402
                                       default_cascade_energy)
from PIC_MCC_3D.MODULES.cross_sections import mixture_set          # noqa: E402
from PIC_MCC_3D.MODULES.gas import GasMixture                      # noqa: E402
from PIC_MCC_3D.MODULES.mcc3d import MCC3D                         # noqa: E402
from PIC_MCC_3D.MODULES.species3d import ParticleSpecies, Species  # noqa: E402
from PIC_MCC_3D.MODULES.units import E_CHARGE, EV, M_E             # noqa: E402

# the engine's bucket conventions (kept in sync with loki_tier/beam)
_LO_HI = {"Ne": 18.0, "Ar": 12.9, "Kr": 11.0, "Xe": 9.0}
_F2_ELEC_MIN, _F2_SPLIT = 2.5, 4.0


def bucket_of(kind, gas, thr):
    if kind == "ionization":
        return gas + "+"
    if kind == "attachment":
        return "F-" if gas == "F2" else None
    if kind == "elastic":
        return "el:" + gas
    if kind == "excitation":
        if gas == "F2":
            if thr < _F2_ELEC_MIN:
                return "vib"
            return "F2(a3)" if thr < _F2_SPLIT else "F2(A1)"
        split = _LO_HI.get(gas)
        if split is None:
            return None
        return (gas + "*") if thr < split else (gas + "**")
    return None


def run_point(mix, cs, mcc, e_over_n_td, n_macro, burn_steps, meas_steps,
              seed, te_seed_eV):
    rng = np.random.default_rng(seed)
    n_gas = mix.total_density_m3
    E_field = e_over_n_td * 1e-21 * n_gas          # V/m (1 Td = 1e-21 V m^2)
    a_z = -E_CHARGE * E_field / M_E                # electron acceleration
    p = ParticleSpecies(Species.electron(), capacity=8 * n_macro)
    vth = np.sqrt(te_seed_eV * EV / M_E)
    v = rng.normal(0.0, vth, size=(n_macro, 3))
    zeros = np.zeros(n_macro)
    p.spawn(zeros, zeros, zeros, v[:, 0], v[:, 1], v[:, 2],
            np.ones(n_macro))
    tally = np.zeros(cs.n_proc, dtype=np.int64)
    t_meas = 0.0
    n_avg = 0.0
    eps_sum = 0.0
    n_samples = 0
    for step in range(burn_steps + meas_steps):
        eps_max = float(p.energy_eV().max())
        nu = cs.nu_majorant(eps_max * 1.1 + 0.5)
        dt = 0.25 / nu * 16.0        # collide_adaptive substeps internally
        p.vz[: p.n][p.alive[: p.n]] += a_z * dt    # kick (whole step)
        ev = mcc.collide_adaptive(p, dt, rng)
        p.compact()
        # population control (weight-preserving)
        if p.n > 2 * n_macro:
            keep = rng.random(p.n) < 0.5
            p.kill(~keep)
            p.compact()
            p.w[: p.n] *= 2.0
        elif 0 < p.n < n_macro // 2:
            n = p.n
            p.spawn(p.x[:n].copy(), p.y[:n].copy(), p.z[:n].copy(),
                    p.vx[:n].copy(), p.vy[:n].copy(), p.vz[:n].copy(),
                    p.w[:n].copy() * 0.5)
            p.w[:n] *= 0.5
        if p.n == 0:
            break
        if step >= burn_steps:
            tally += ev["tally"]
            t_meas += dt
            n_avg += p.n * dt
            eps_sum += float(np.mean(p.energy_eV())) * dt
            n_samples += 1
    if t_meas <= 0.0 or n_avg <= 0.0:
        return None
    n_macro_avg = n_avg / t_meas
    eps_mean = eps_sum / t_meas
    k_proc = tally / (n_macro_avg * t_meas * cs.density)
    err_proc = np.where(tally > 0, 1.0 / np.sqrt(np.maximum(tally, 1)),
                        np.nan)
    return eps_mean, k_proc, err_proc, tally


def cascade_w_per_gas(gases, n_electrons, seed):
    """The SECOND W route: a Monte-Carlo degradation cascade through
    PIC_MCC_3D's ACTUAL MCC operator (beam3d.cascade_yields), per pure gas.

    REBUILD's own Spencer–Fano route is deterministic; this one exercises
    independent code (null-collision selection, kinematics, sharing) on the
    same LXCat data. PIC_MCC_3D's in-package pairing of the same two
    methods agrees to 0.4 % (Ar) — the scorecard compares REBUILD-SF
    against THIS number at the recorded E0, so the check crosses package
    boundaries. Integral W = E0/U_ion at E0 = 0.55 x ceiling (the
    differential-window rule does not apply to a single-E0 cascade; both
    routes are evaluated at the SAME E0, so the comparison is like-for-like).
    """
    from PIC_MCC_3D.MODULES.cross_sections import (CrossSectionSet as PicCS,
                                                   GAS_FILES as PIC_GAS_FILES)
    out = {}
    for i, gas in enumerate(gases):
        t0 = time.time()
        mix = GasMixture({gas: 1.0}, pressure_atm=1.0, T_gas_K=300.0)
        if gas in PIC_GAS_FILES:
            cs = mixture_set(mix)
        else:
            # PIC_MCC_3D's loader maps Ar/Kr/F2/He only; Ne's Biagi file is
            # reachable on both packages' search paths but unmapped there.
            # Resolve it with REBUILD's loader and hand PIC_MCC_3D the
            # explicit path — no sibling-package edit (the ARPAE-Loschmidt
            # rule: changing another package as a side effect is not this
            # tool's call).
            from REBUILD.MODULES.cross_sections import gas_file as _rb_file
            cs = PicCS({gas: (_rb_file(gas), mix.density_of(gas))})
        mcc = MCC3D(cs, mixture=mix)
        E0 = default_cascade_energy(cs)
        rng = np.random.default_rng(seed + 1000 + i)
        y = cascade_yields(mcc, E0, n_electrons=n_electrons, rng=rng)
        out[gas] = {
            "E0_eV": float(E0), "W_eV": float(y["W_eV"]),
            "U_ion": float(y["U_ion"]), "n_electrons": int(n_electrons),
            "rounds": int(y["rounds"]), "converged": bool(y["converged"]),
            "stat_err_rel": float(1.0 / np.sqrt(max(
                y["U_ion"] * n_electrons, 1.0))),
        }
        print(f"  cascade W({gas}) @ E0 {E0:.0f} eV -> "
              f"{y['W_eV']:.2f} eV ({y['rounds']} rounds, "
              f"{time.time()-t0:.0f} s)")
    return out


def relaxation_sampling(mix, cs, mcc, n_macro, seed, E0_eV=60.0,
                        max_rounds=20000):
    """The L3 measurement: rates of a field-free RELAXING population.

    Beam-sustained plasmas are cascade/relaxation-shaped, not field-swarm
    shaped, and their mean energies (~0.5-2 eV) sit BELOW the swarm
    tables, where the engine's flat+taper extrapolation currently stands
    in. Here electrons launched at E0 relax collisionally with no field;
    per-process tallies are binned by the population's instantaneous mean
    energy, giving k_p(<eps>) down to ~0.5 eV from the ACTUAL operators.
    Returns (bin centers eV, {bucket: k}, {bucket: rel err})."""
    rng = np.random.default_rng(seed + 5000)
    p = ParticleSpecies(Species.electron(), capacity=8 * n_macro)
    v = np.sqrt(2.0 * E0_eV * EV / M_E)
    d = rng.normal(size=(n_macro, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    z = np.zeros(n_macro)
    p.spawn(z, z, z, v * d[:, 0], v * d[:, 1], v * d[:, 2],
            np.ones(n_macro))
    edges = np.geomspace(0.4, 40.0, 15)
    nb = edges.size - 1
    tally = np.zeros((nb, cs.n_proc))
    exposure = np.zeros(nb)              # sum of N_macro * dt per bin
    for _ in range(max_rounds):
        if p.n == 0:
            break
        eps = p.energy_eV()
        mean = float(np.mean(eps))
        if mean < edges[0]:
            break
        b = int(np.searchsorted(edges, mean) - 1)
        b = min(max(b, 0), nb - 1)
        nu = cs.nu_majorant(float(eps.max()) * 1.1 + 0.5)
        dt = 0.25 / nu * 16.0
        ev = mcc.collide_adaptive(p, dt, rng)
        p.compact()
        tally[b] += ev["tally"]
        exposure[b] += p.n * dt
    centers = np.sqrt(edges[:-1] * edges[1:])
    ks, errs = {}, {}
    for j in range(cs.n_proc):
        gas = cs.targets[j].split("(")[0].strip()
        b = bucket_of(cs.kinds[j], gas, float(cs.thresholds[j]))
        if b is None:
            continue
        ks.setdefault(b, np.zeros(nb))
        errs.setdefault(b, np.zeros(nb))
        with np.errstate(divide="ignore", invalid="ignore"):
            ks[b] += np.where(exposure > 0,
                              tally[:, j] / (exposure * cs.density[j]), 0.0)
        errs[b] += tally[:, j]
    for b in errs:
        with np.errstate(divide="ignore", invalid="ignore"):
            errs[b] = np.where(errs[b] > 0, 1.0 / np.sqrt(errs[b]), np.nan)
    return centers, ks, errs


def deepen(max_rounds: int) -> None:
    """Deepen ONLY the relaxation curves past the Ramsauer stall — the
    recorded L3 follow-up (REBUILD CHANGELOG 0.1.9), now with its
    strongest reason (the F- x0.02 residual vs the NRL population
    benchmark, 0.1.12: the measured F- bins stop at 1.76 eV exactly
    where physical dissociative attachment rises).

    PURE EXTENSION, mechanically enforced: relaxation_sampling seeds its
    own rng (seed + 5000, independent of the swarm/cascade streams), so
    re-running it with the SAME seed and MORE rounds replays rounds
    1..N_old bit-identically and only then continues. The population
    mean decreases monotonically (field-free relaxation), so every bin
    STRICTLY ABOVE the old stall bin receives zero new tallies and must
    reproduce BIT-EXACTLY — asserted below; the stall bin itself gains
    exposure (refines within its statistics) and new bins appear under
    it. The swarm points and W-cascade arrays are carried over UNTOUCHED
    from the shipped npz (their streams never run)."""
    out_dir = _HERE.parent / "MODULES" / "kinetic_tables"
    path = out_dir / "pic_mcc_3d_swarm.npz"
    old = dict(np.load(path, allow_pickle=False))
    meta = json.loads(str(old["meta"]))
    fractions = meta["mixture"]
    seed = int(meta["seed"])
    mix = GasMixture(fractions, pressure_atm=1.0, T_gas_K=300.0)
    cs = mixture_set(mix)
    mcc = MCC3D(cs, mixture=mix)
    print(f"deepen: relax-only, seed {seed}, max_rounds {max_rounds} "
          f"(shipped run used 25000)")
    t0 = time.time()
    r_eps, r_k, r_err = relaxation_sampling(mix, cs, mcc, 500, seed,
                                            max_rounds=max_rounds)
    print(f"  sampled in {time.time()-t0:.0f} s")
    assert np.array_equal(r_eps, old["relax_eps_eV"]), "bin grid moved"
    # --- the pure-extension gate --------------------------------------
    checked = frozen = 0
    for key in list(old):
        if not key.startswith("relaxk_"):
            continue
        b = key[7:]
        k_old = old[key]
        k_new = r_k.get(b.replace("_", ":") if b.startswith("el_")
                        else b, None)
        if k_new is None:
            k_new = r_k.get(b)
        assert k_new is not None, f"bucket {b} vanished in the deepen run"
        covered = np.nonzero(k_old > 0)[0]
        if covered.size == 0:
            continue
        stall = int(covered.min())
        if not np.array_equal(k_new[stall + 1:], k_old[stall + 1:]):
            raise SystemExit(
                f"PURE-EXTENSION VIOLATED for {b}: a bin above the old "
                f"stall bin moved — the stream is not a replay; refusing "
                f"to write")
        checked += 1
        frozen += k_old.size - (stall + 1)
    print(f"  pure-extension gate: {checked} buckets, {frozen} "
          f"above-stall bin values bit-identical")
    # --- merge and write ----------------------------------------------
    meta["relaxation_deepened"] = {
        "max_rounds": int(max_rounds), "date": str(date.today()),
        "note": "same seed, relax stream replayed + extended past the "
                "Ar Ramsauer stall (REBUILD 0.1.13); swarm/cascade "
                "arrays carried over untouched",
    }
    old["meta"] = json.dumps(meta)
    for b, k in r_k.items():
        key = b.replace(":", "_")
        old["relaxk_" + key] = k
        old["relaxerr_" + key] = r_err[b]
    np.savez(path, **old)
    got = sorted(b for b in r_k if np.any(r_k[b] > 0))
    lo = {b: float(r_eps[np.nonzero(r_k[b] > 0)[0].min()])
          for b in got}
    print(f"wrote {path}")
    print("  new coverage floors [eV]:",
          {b: round(v, 2) for b, v in sorted(lo.items())})


def main():
    if "--deepen" in sys.argv:
        i = sys.argv.index("--deepen")
        mr = int(sys.argv[i + 1]) if len(sys.argv) > i + 1 else 250000
        deepen(mr)
        return
    quick = "--quick" in sys.argv
    fractions = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}
    mix = GasMixture(fractions, pressure_atm=1.0, T_gas_K=300.0)
    cs = mixture_set(mix)
    mcc = MCC3D(cs, mixture=mix)
    en_points = ([3.0, 7.0, 12.0, 20.0, 50.0, 100.0, 200.0]
                 if not quick else [50.0, 150.0])
    n_macro = 800 if not quick else 300
    burn, meas = (2500, 4000) if not quick else (400, 600)
    seed = 20260804
    print(f"swarm export: {fractions} | E/N {en_points} Td | "
          f"{n_macro} macros, {burn}+{meas} steps")
    eps_list = []
    k_rows = []
    err_rows = []
    for i, en in enumerate(en_points):
        t0 = time.time()
        te_seed = max(0.8, 0.03 * en)          # start near the answer
        out = run_point(mix, cs, mcc, en, n_macro, burn, meas,
                        seed + i, te_seed)
        if out is None:
            print(f"  {en} Td: population died, skipped")
            continue
        eps_mean, k_proc, err_proc, tally = out
        eps_list.append(eps_mean)
        k_rows.append(k_proc)
        err_rows.append(err_proc)
        print(f"  {en:6.1f} Td -> <eps> {eps_mean:5.2f} eV, "
              f"{int(tally.sum())} events, {time.time()-t0:.0f} s")
    print("second W route (MC cascade through PIC_MCC_3D's operator):")
    w_cascade = cascade_w_per_gas(["Ar"] if quick else ["Ar", "Kr", "Ne"],
                                  150 if quick else 400, seed)
    print("relaxation sampling (L3: rates at beam-mode mean energies):")
    t0 = time.time()
    r_eps, r_k, r_err = relaxation_sampling(
        mix, cs, mcc, 200 if quick else 500, seed,
        max_rounds=4000 if quick else 25000)
    got = [b for b in r_k if np.any(r_k[b] > 0)]
    print(f"  bins {r_eps[0]:.2f}-{r_eps[-1]:.1f} eV, buckets with data: "
          f"{sorted(got)} ({time.time()-t0:.0f} s)")
    order = np.argsort(eps_list)
    eps = np.array(eps_list)[order]
    K = np.array(k_rows)[order]
    ERR = np.array(err_rows)[order]
    # aggregate processes into engine buckets
    buckets = {}
    errs = {}
    for j in range(cs.n_proc):
        gas = cs.targets[j].split("(")[0].strip()   # 'Kr(1S0)' -> 'Kr'
        b = bucket_of(cs.kinds[j], gas, float(cs.thresholds[j]))
        if b is None:
            continue
        buckets.setdefault(b, np.zeros_like(eps))
        errs.setdefault(b, np.full_like(eps, np.nan))
        buckets[b] = buckets[b] + np.nan_to_num(K[:, j])
        e_old = errs[b]
        errs[b] = np.where(np.isnan(e_old), ERR[:, j],
                           np.fmin(e_old, ERR[:, j]))
    out_dir = _HERE.parent / "MODULES" / "kinetic_tables"
    out_dir.mkdir(exist_ok=True)
    path = out_dir / "pic_mcc_3d_swarm.npz"
    meta = {
        "source": "PIC_MCC_3D.MODULES.mcc3d.MCC3D swarm driver "
                  "(this script)",
        "method": "uniform-field null-collision swarm, immobile ions, "
                  "population band-controlled; k = events/(N T n_target)",
        "mixture": fractions, "pressure_atm_run": 1.0,
        "e_over_n_td": [float(e) for e in en_points],
        "n_macro": n_macro, "burn_steps": burn, "meas_steps": meas,
        "seed": seed, "date": str(date.today()),
        "mapping": "<eps> = (3/2) Te lookup, same convention as the LOKI "
                   "tier (see kinetic_tier.py docstring)",
        "w_cascade": w_cascade,
        "w_cascade_method": "beam3d.cascade_yields through MCC3D per pure "
                            "gas; integral W = E0/U_ion at the recorded E0 "
                            "(0.55 x table ceiling); compare REBUILD "
                            "Spencer-Fano AT THE SAME E0",
    }
    meta["relaxation"] = ("field-free relaxing-population rates binned by "
                          "instantaneous <eps> (L3 measurement; E0 = 60 eV, "
                          "bins 0.4-40 eV) — the beam-mode EEDF class the "
                          "swarm tables cannot represent")
    arrays = {"eps_eV": eps, "meta": json.dumps(meta),
              "relax_eps_eV": r_eps}
    for b, k in r_k.items():
        key = b.replace(":", "_")
        arrays["relaxk_" + key] = k
        arrays["relaxerr_" + key] = r_err[b]
    for b, k in buckets.items():
        key = b.replace(":", "_")
        arrays["k_" + key] = k
        arrays["err_" + key] = errs[b]
    np.savez(path, **arrays)
    print(f"wrote {path} ({sorted(buckets)})")


if __name__ == "__main__":
    main()
