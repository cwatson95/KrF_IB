"""Command line for LIGHT 1Dv0_2, the 1-D model (launch through LIGHT/run_1Dv0_2.py).

    python LIGHT/run_1Dv0_2.py                   # v0_2's reference point, 1-cm gap
    python LIGHT/run_1Dv0_2.py --compare         # and v0_2's 0-D run of the same inputs
    python LIGHT/run_1Dv0_2.py --cells 256 --gap 2 --dt-ps 25
    python LIGHT/run_1Dv0_2.py --pressure 2 --power 500 --f2-pct 0.3
    python LIGHT/run_1Dv0_2.py --rate kr2_f=3e-10  # one band value (repeatable)
    python LIGHT/run_1Dv0_2.py --profile power.txt # relative power per cell

Needs Python >= 3.11 with numba, as Fluid1D does (~/.venvs/f2sweep/bin/python
on this Mac). Fluid1D is read from Kinetic_Modeling/nathaniel/1DModelKrF/src
(or FLUID1D_SRC). Writes <out>/<name>/result.json and the files Fluid1D read
(fluid1d_inputs/); <out> defaults to LIGHT/1Dv0_2/runs/. v0_2's --bands (each
band's two ends, one run per end) stays with run_v0_2.py: with a uniform
input the gap's interior is the 0-D run.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from .network import BANDS, BASIC_GROUPS
from .sim1d import Params1D, run


def _rate(text):
    """--rate NAME=VALUE, as run_v0_2.py reads it."""
    name, _, value = text.partition("=")
    if name not in BANDS or not value:
        raise SystemExit(f"--rate {text!r}: give NAME=VALUE with NAME one "
                         f"of {', '.join(BANDS)}")
    return name, float(value)


def params_from_args(a) -> Params1D:
    groups = ()
    if a.add:
        groups = (tuple(BASIC_GROUPS) if a.add == "all"
                  else tuple(g.strip() for g in a.add.split(",") if g.strip()))
    profile = None
    if a.profile:
        profile = tuple(float(v) for v in
                        np.loadtxt(a.profile, dtype=float).ravel())
    return Params1D(
        pressure_atm=a.pressure, x_f2=a.f2_pct / 100.0,
        tg0_K=a.tg0, power_W_cm3=a.power * 1e3, pulse_s=a.pulse * 1e-6,
        t_end_s=a.t_end * 1e-6, ne_cm3=a.ne, tau_sp_s=a.tau_sp * 1e-9,
        trimer_tg_pow=a.trimer_tg_pow,
        rates=tuple(_rate(x) for x in a.rate), groups=groups, eedf=a.eedf,
        capture=a.capture, gap_cm=a.gap, area_cm2=a.area, cells=a.cells,
        dt_s=a.dt_ps * 1e-12, power_profile=profile, name=a.name,
        output_dir=a.out, threads=a.threads, verbose=not a.quiet)


def compare_0d(res, p: Params1D) -> str:
    """v0_2's 0-D run of the same inputs (1Dv0_2's own model.py) beside the
    1-D run's centre cell and gap."""
    from .model import run as run_0d
    from .report import row_dict
    d0 = row_dict(run_0d(p.cell_params()))
    last = res.data["records"][-1]
    vuv = "rad:147 nm (Kr2* VUV)"
    led, led_c = last["ledger"], last["ledger_centre"]
    rows = [("248 nm, % of deposited", d0["emission_248"],
             res.share_248("ledger_centre"), res.share_248())]
    if vuv in led:
        rows.append(("147 nm, % of deposited", d0["emission_vuv"],
                     100 * led_c[vuv] / led_c["input"],
                     100 * led[vuv] / led["input"]))
    rows.append(("Tg at the end, K", d0["tg_end_K"], last["tg_centre"],
                 last["tg_mean"]))
    out = ["", "0-D (v0_2's model) vs 1-D         0-D     centre      gap"]
    for label, a0, c, g in rows:
        out.append(f"  {label:30s} {a0:9.4f} {c:9.4f} {g:9.4f}")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="run_1Dv0_2.py", description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pressure", type=float, default=1.0, help="atm")
    ap.add_argument("--f2-pct", type=float, default=0.5,
                    help="F2 mole percent, rest Kr")
    ap.add_argument("--tg0", type=float, default=300.0, help="K")
    ap.add_argument("--power", type=float, default=300.0, help="kW/cm^3")
    ap.add_argument("--pulse", type=float, default=1.0, help="us")
    ap.add_argument("--t-end", type=float, default=2.0, help="us")
    ap.add_argument("--ne", type=float, default=1e15,
                    help="cm^-3, the seed electron (= Kr+) density")
    ap.add_argument("--tau-sp", type=float, default=1e9 / 1.43e8,
                    help="KrF(B) radiative lifetime, ns")
    ap.add_argument("--trimer-tg-pow", type=float, default=None,
                    help="sensitivity only: the trimer's (Tg/300)^n exponent "
                    "(default Shui's -3.32; -0.5 is the pure-Kr efficiency "
                    "bound)")
    ap.add_argument("--rate", action="append", default=[],
                    metavar="NAME=VALUE", help="set one band (repeatable), "
                    "sensitivity only: " + ", ".join(
                        b for b in BANDS if b != "trimer_tg_pow"))
    ap.add_argument("--add", default="",
                    help="optional groups, comma list: "
                    + ",".join(BASIC_GROUPS) + " or all")
    ap.add_argument("--eedf", choices=("auto", "two_t", "maxwellian"),
                    default="auto", help="auto: v0_2's two-temperature tail")
    ap.add_argument("--capture", choices=("auto", "exact", "three_halves"),
                    default="auto")
    ap.add_argument("--gap", type=float, default=1.0, help="cm")
    ap.add_argument("--area", type=float, default=1.0,
                    help="cm^2 (scales totals only)")
    ap.add_argument("--cells", type=int, default=128)
    ap.add_argument("--dt-ps", type=float, default=50.0,
                    help="transport step, ps")
    ap.add_argument("--profile", default=None,
                    help="file of relative power per cell (one number per "
                    "cell, any scale); default uniform, v0_2's rule")
    ap.add_argument("--name", default="1Dv0_2")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--threads", type=int, default=None,
                    help="numba threads")
    ap.add_argument("--compare", action="store_true",
                    help="also run v0_2's 0-D model on the same inputs")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    p = params_from_args(a)
    t0 = time.time()
    res = run(p)
    print(res.summary())
    if a.compare:
        print(compare_0d(res, p))
    root = Path(p.output_dir) if p.output_dir else None
    where = (root / p.name if root else Path(__file__).resolve().parent
             / "runs" / p.name) / "result.json"
    print(f"\nwrote {where} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
