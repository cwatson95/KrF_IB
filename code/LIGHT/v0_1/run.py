"""Command line for LIGHT v0_1 (launch through LIGHT/run_v0_1.py).

    python3 LIGHT/run_v0_1.py                        # the reference point
    python3 LIGHT/run_v0_1.py --pressure 2 --power 500 --ne 1e16
    python3 LIGHT/run_v0_1.py --ne-mode free         # n_e from ionization
    python3 LIGHT/run_v0_1.py --network v0           # the notes' list as written
    python3 LIGHT/run_v0_1.py --add kr2              # promote optional groups
    python3 LIGHT/run_v0_1.py --omissions            # the network vs each group
    python3 LIGHT/run_v0_1.py --scan --csv scan.csv  # the notes' ranges

The only energy input is the given power, deposited uniformly and directly
into the electron gas.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time

from .model import Model, Params, run
from .network import BASIC_GROUPS, GROUPS
from .report import row_dict, summary

SCAN_PRESSURE = (0.1, 0.25, 0.5, 1.0, 2.0)          # atm  (notes: 0.1-2)
SCAN_POWER = (100.0, 300.0, 500.0)                  # kW/cm^3 (notes)
SCAN_NE = (1e14, 1e15, 1e16)                        # cm^-3 (notes)

TABLE_COLS = (
    ("pressure_atm", "p/atm", "{:5.2f}"), ("power_kW_cm3", "kW/cc", "{:5.0f}"),
    ("ne_cm3", "n_e", "{:7.0e}"), ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Krstar", "Kr*", "{:5.1f}"), ("to_Krplus", "Kr+", "{:5.1f}"),
    ("to_KrFB", "KrF(B)", "{:6.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_248", "248nm", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("f2_left_pct", "F2left", "{:6.1f}"),
    ("tg_end_K", "Tg_end", "{:6.0f}"), ("closure", "closure", "{:8.0e}"),
)


def params_from_args(a) -> Params:
    kw = dict(network=a.network, pressure_atm=a.pressure,
              x_f2=a.f2_pct / 100.0,
              power_W_cm3=a.power * 1e3, pulse_s=a.pulse * 1e-6,
              t_end_s=a.t_end * 1e-6, ne_cm3=a.ne, ne_mode=a.ne_mode,
              tau_sp_s=a.tau_sp * 1e-9, rtol=a.rtol,
              trimer_tg_pow=a.trimer_tg_pow,
              groups=tuple(g for g in a.add.split(",") if g))
    return Params(**kw)


def table(rows) -> str:
    head = "  ".join(f"{h:>{len(fmt.format(0))}s}" for _, h, fmt in TABLE_COLS)
    out = [head]
    for r in rows:
        out.append("  ".join(fmt.format(r[k]) for k, _, fmt in TABLE_COLS))
    return "\n".join(out)


def do_scan(base: Params):
    """The notes' pressure x power grid. With n_e held (the notes' list,
    network v0) n_e is a third axis. With n_e free it is only the seed and
    is not scanned."""
    ne_axis = SCAN_NE if Model(base).ne_mode == "clamped" else (base.ne_cm3,)
    rows = []
    for p in SCAN_PRESSURE:
        for pw in SCAN_POWER:
            for ne in ne_axis:
                res = run(base.replace(pressure_atm=p, power_W_cm3=pw * 1e3,
                                       ne_cm3=ne))
                rows.append(row_dict(res))
    return rows


def do_omissions(base: Params):
    rows = []
    allowed = BASIC_GROUPS if base.network == "basic" else GROUPS
    for g in ((),) + tuple((x,) for x in allowed) + (("all",),):
        d = row_dict(run(base.replace(groups=g)))
        rows.append(d)
    return rows


OMIT_COLS = (
    ("groups", "network", "{:>12s}"), ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Krstar", "Kr*", "{:5.1f}"), ("to_Krplus", "Kr+", "{:5.1f}"),
    ("to_KrFB", "KrF(B)", "{:6.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_248", "248nm", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("krf_to_248_pct", "B->248",
                                          "{:6.1f}"),
    ("krstar_to_krf_pct", "Kr*->B", "{:6.1f}"),
    ("f2_left_pct", "F2left", "{:6.1f}"), ("closure", "closure", "{:8.0e}"),
)


def omit_table(rows) -> str:
    def w(fmt):
        return len(fmt.format("x" * 12 if "s}" in fmt else 0))
    head = "  ".join(f"{h:>{w(fmt)}s}" for _, h, fmt in OMIT_COLS)
    out = [head]
    for r in rows:
        out.append("  ".join(fmt.format(r[k]) for k, _, fmt in OMIT_COLS))
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="run_v0_1.py", description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pressure", type=float, default=1.0, help="atm")
    ap.add_argument("--f2-pct", type=float, default=0.5,
                    help="F2 mole percent, rest Kr (not set by the notes)")
    ap.add_argument("--power", type=float, default=300.0, help="kW/cm^3")
    ap.add_argument("--pulse", type=float, default=1.0, help="us")
    ap.add_argument("--t-end", type=float, default=2.0, help="us")
    ap.add_argument("--ne", type=float, default=1e15,
                    help="cm^-3, the initial electron density (a seed)")
    ap.add_argument("--ne-mode", choices=("auto", "free", "clamped"),
                    default="auto", help="auto: free for the basic "
                    "network, held for the notes' list (v0)")
    ap.add_argument("--trimer-tg-pow", type=float, default=None,
                    help="sensitivity only: the trimer's (Tg/300)^n "
                    "exponent (default Shui's -3.32; -0.5 is the "
                    "pure-Kr efficiency bound)")
    ap.add_argument("--tau-sp", type=float, default=1e9 / 1.43e8,
                    help="KrF(B) radiative lifetime, ns (notes: ~4)")
    ap.add_argument("--network", choices=("basic", "v0"), default="basic",
                    help="basic: the KrF(B,C) network (default); v0: the "
                    "notes' six reactions as written")
    ap.add_argument("--add", default="",
                    help="optional groups, comma list: " + ",".join(GROUPS)
                    + " (basic: " + ",".join(BASIC_GROUPS) + ")"
                    + " or all")
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--scan", action="store_true",
                    help="p x power over the notes' ranges (x n_e when "
                    "n_e is held)")
    ap.add_argument("--omissions", action="store_true",
                    help="V0 vs V0 + each optional group vs V0 + all")
    ap.add_argument("--csv", default=None, help="write scan/omissions rows")
    ap.add_argument("--json", default=None, help="write the single-run rows")
    a = ap.parse_args(argv)
    base = params_from_args(a)
    t0 = time.time()
    if a.scan or a.omissions:
        rows = do_scan(base) if a.scan else do_omissions(base)
        print(table(rows) if a.scan else omit_table(rows))
        print(f"\n(percent of deposited energy; {len(rows)} runs in "
              f"{time.time() - t0:.1f} s)")
        if a.csv:
            with open(a.csv, "w", newline="") as fh:
                wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
                wr.writeheader()
                wr.writerows(rows)
        return 0
    res = run(base)
    print(summary(res))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({k: (float(v) if not isinstance(v, str) else v)
                       for k, v in row_dict(res).items()}, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
