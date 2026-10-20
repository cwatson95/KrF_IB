"""The EEDF of the laser plasma, and what it does to LIGHT.

LIGHT integrates every electron rate over a Maxwellian at its Te. At Tahoma
(n_e ~ 5e14 cm^-3 in 2e19 cm^-3 of gas, ionization degree 2e-5) the
kinetic answer is the steady Boltzmann equation of a beam-sustained plasma
(ml/boltzmann/ee.py, solve_beam):
  - electrons enter with the sub-excitation spectrum of the degradation
    cascade, 4.506 eV each (the beam drive's Spencer-Fano value) on
    [0, 9.915 eV];
  - they collide with Kr and F2 (the same md5-pinned LoKI-MC cross sections
    LIGHT reads) and with one another (Coulomb, the term LoKI-MC lacks);
  - they leave by attachment and by recombination (LIGHT's R11p rate per
    electron).
Its n_e and mean energy come out of the solution.

The comparison is made at the same mean energy. The solution's bulk is
Maxwellian to 2-4 % in every rate below 9 eV (attachment, vibrational,
dissociative). Above the Kr thresholds, excitation empties the tail faster
than e-e collisions refill it. The factor
    f_row = k_row(Boltzmann) / k_row(Maxwellian at 2/3 <eps>)
therefore carries the shape and not a temperature. LIGHT keeps its own
energy equation and Te, with each electron row scaled by f_row.

LIGHT's Maxwellian closure (forward.MAXWELL, LIGHT's own until
2026-09-28) removes (3/2)Te with every attached electron. F2 attaches below
1 eV, so at LIGHT's Tahoma state an attaching electron carries off 0.54 eV
(the Boltzmann solution, and a Maxwellian gives about the same), not the
1.69 eV that closure books. The overbooked cooling is 21 % of the
secondary-electron heating. With exact_capture (the default) the loop
books e_att per attachment and Te per R11p recombination (exact for its
Te^-1/2 template, and (3/2 + n) Te for v0_1's Te^n law). The Boltzmann
solve uses the matching sigma v ~ eps^n sink. Since 2026-09-28 LIGHT's basic network books both itself and takes
its tail from ml/tail.py's table. This loop is the reference that closure
is tested against (tests/test_electron.py).

self_consistent iterates LIGHT(f) -> window state -> Boltzmann -> f until f
stops changing. At convergence it checks LIGHT's Te and n_e against the
Boltzmann solution's 2/3 <eps> and n_e. They agree only if LIGHT's energy
equation and the kinetic one close the same budget: with exact capture to
0.4 % and 5 %, with (3/2)Te capture only to 7 % and 2 % (2026-09-27).

The scored outputs (Tahoma g0, formation) never steer anything here. The
EEDF follows from the cross sections and the collision physics alone.
"""
from __future__ import annotations

import numpy as np

from . import ROOT, electron, forward, light
from . import boltzmann  # noqa: F401  (puts the vendored solver on sys.path)
import boltz_twin as BT  # noqa: E402
import ee  # noqa: E402
import gases_twin as GT  # noqa: E402
import lx_lib as L  # noqa: E402

INP = str(ROOT / "LoKI-MC" / "Code" / "Input")
GT.load_tables(INP)

E_TOP = 9.915      # eV: the lowest Kr threshold, top of the sub-excitation band
EPS_SEC = 4.506    # eV per ion pair (LIGHT's beam drive, REBUILD Spencer-Fano)
# LIGHT's electron rows -> (gas, threshold eV) of the solver's collision
ROWS = {"G1": ("Kr", 9.915), "C1a": ("Kr", 11.30), "C1b": ("Kr", 11.998),
        "C1c": ("Kr", 12.75), "C2": ("Kr", 13.996), "ATT": ("F2", 0.0),
        "VIB1": ("F2", 0.111), "VIB2": ("F2", 0.219), "VIB3": ("F2", 0.324),
        "VIB4": ("F2", 0.421), "DIS_a3": ("F2", 3.16), "DIS_A1": ("F2", 4.34)}
_GASES: dict = {}


def r11p_pow() -> float:
    """The Te exponent of LIGHT's Kr2+ + e row: -1/2 for v0's Te^-1/2
    template, the row's own for a te_pow law (v0_1's -0.55)."""
    row = {r.rid: r for r in light("network").build(("all",), "basic")}["R11p"]
    if row.rate[0] == "dr_sqrt":
        return -0.5
    if row.rate[0] == "te_pow":
        return float(row.rate[3])
    raise ValueError(f"R11p: rate kind {row.rate[0]!r}")


def gases(tg):
    key = round(tg, 1)
    if key not in _GASES:
        _GASES[key] = (GT.build("Kr", tg, INP), GT.build("F2", tg, INP))
    return _GASES[key]


def collision(kr, f2, gid, thr):
    gas = kr if gid == "Kr" else f2
    return next(c for c in gas["colls"] if abs(c["thr"] - thr) < 1e-6)


def row_rates(sol, kr, f2):
    """<sigma v> [cm^3/s] of LIGHT's electron rows over sol's F0."""
    return {r: ee.rate(collision(kr, f2, *v), sol) * 1e6
            for r, v in ROWS.items()}


def beam_eedf(dev, ee_on=True, mean_inj=EPS_SEC, N=600, emax=30.0,
              att=(1.0, 1.0), dr=False, tol=1e-9):
    """The Boltzmann EEDF at a device run's window state (forward.run_device
    output). Returns the solution, its row rates, the rates of the
    Maxwellian with the same mean energy, the factors f_row, and the mean
    energy an attaching electron carries off (e_att).
    att: the (s_low, s_high) scales of the F2 attachment resonances.
    dr: recombination as LIGHT's R11p template (sigma ~ 1/eps) instead of
    the same for every electron."""
    n = dev["n"]
    tg = dev["tg"]
    kr, f2 = gases(tg)
    if tuple(att) != (1.0, 1.0):
        from .swarm import f2_modified
        f2 = f2_modified(tg, *att)
    n_gas = (n["Kr"] + n["F2"] + n["F"]) * 1e6
    mix = [(kr, n["Kr"] * 1e6 / n_gas), (f2, n["F2"] * 1e6 / n_gas)]
    pairs = dev["pairs"] * 1e6
    nu_rec = dev["events"]["R11p"] / dev["ne"]
    sol = ee.solve_beam(mix, tg, n_gas, N, emax, pairs,
                        lambda e: ee.injection(e, E_TOP, mean_inj), nu_rec,
                        dev["ne"] * 1e6, ee_on=ee_on,
                        rec_te=dev["te"] if dr else None, tol=tol,
                        rec_pow=r11p_pow())
    k = row_rates(sol, kr, f2)
    ca = collision(kr, f2, "F2", 0.0)
    w = sol["ec"] * L.interp(ca, sol["ec"]) * sol["F"]
    e_att = float(np.sum(sol["ec"] * w) / np.sum(w))
    same = dict(sol, F=ee.maxwellian(sol["ec"], sol["du"],
                                     2.0 / 3.0 * sol["mean"]))
    km = row_rates(same, kr, f2)
    return {"ne": sol["ne"] / 1e6, "mean": sol["mean"], "k": k, "k_max": km,
            "f": {r: k[r] / km[r] for r in ROWS if km[r] > 0}, "sol": sol,
            "e_att": e_att}


def self_consistent(c, cond=None, ee_on=True, att=(1.0, 1.0),
                    exact_capture=True, tol=0.005, maxit=20):
    """LIGHT with its electron rows scaled by the Boltzmann factors, iterated
    to a fixed point. LIGHT runs on its Maxwellian closure here
    (forward.MAXWELL), the reference its own two-temperature closure is
    tested against. att: trained F2 attachment scales, applied to LIGHT's
    ATT table and to the Boltzmann solve alike. exact_capture: book the
    energy an attached electron really carries off (e_att, from the
    Boltzmann solution) instead of (3/2)Te, and (3/2 + n) Te per R11p
    recombination for its Te^n law (Te for v0's Te^-1/2 template), with the
    matching sigma v ~ eps^n sink in the Boltzmann solve. Returns (device result, last beam_eedf,
    history)."""
    tables = ({} if tuple(att) == (1.0, 1.0)
              else {"ATT": electron.att_table(*att)})
    cap = {}

    def run(f):
        return forward.run_device(c, cond, modify=electron.modifier(
            scales=f, tables=tables, capture=cap), closure=forward.MAXWELL)
    f = {}
    hist = []
    for it in range(maxit):
        dev = run(f)
        b = beam_eedf(dev, ee_on=ee_on, att=att, dr=exact_capture)
        if exact_capture:
            new_cap = {"ATT": b["e_att"] / (1.5 * dev["te"]),
                       "R11p": (1.5 + r11p_pow()) / 1.5}
            cap_change = max(abs(np.log(new_cap[r] / cap[r])) if r in cap
                             else np.inf for r in new_cap)
            cap.update({r: (np.sqrt(cap[r] * v) if r in cap else v)
                        for r, v in new_cap.items()})
        else:
            cap_change = 0.0
        hist.append({"g0": dev["g0"], "formation": dev["formation"],
                     "ne": dev["ne"], "te": dev["te"], "ne_B": b["ne"],
                     "te_B": 2.0 / 3.0 * b["mean"]})
        # damped in log space: the factors span 7 decades
        new = {r: (np.sqrt(f[r] * v) if r in f else v)
               for r, v in b["f"].items()}
        change = max(max(abs(np.log(new[r] / f[r])) if r in f else np.inf
                         for r in new), cap_change)
        f = new
        if change < tol:
            break
    else:
        raise RuntimeError(f"LIGHT-Boltzmann loop: change {change:.3f}")
    dev = run(f)
    b = beam_eedf(dev, ee_on=ee_on, att=att, dr=exact_capture)
    hist.append({"g0": dev["g0"], "formation": dev["formation"],
                 "ne": dev["ne"], "te": dev["te"], "ne_B": b["ne"],
                 "te_B": 2.0 / 3.0 * b["mean"]})
    dev["factors"] = f
    dev["capture"] = dict(cap)
    return dev, b, hist
