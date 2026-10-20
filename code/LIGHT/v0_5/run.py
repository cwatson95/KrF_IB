"""Command line for LIGHT v0_4 (launch through LIGHT/run_v0_5.py).

    python3 LIGHT/run_v0_5.py                        # the reference point
    python3 LIGHT/run_v0_5.py --bands                # each band's two ends
    python3 LIGHT/run_v0_5.py --rate kr2_f=3e-10     # one band value
    python3 LIGHT/run_v0_5.py --pressure 2 --power 500 --ne 1e16
    python3 LIGHT/run_v0_5.py --network v0           # the notes' list as written
    python3 LIGHT/run_v0_5.py --network v0 --add kr2 # v0 and an optional group
    python3 LIGHT/run_v0_5.py --omissions            # the network vs each group
    python3 LIGHT/run_v0_5.py --scan --csv scan.csv  # the notes' ranges

    python3 LIGHT/run_v0_5.py --ladder kannari       # every 5p-rung band at Kannari's value
    python3 LIGHT/run_v0_5.py --k76 kannari          # v0_4's 5s -> 5p fit instead of the measurement
    python3 LIGHT/run_v0_5.py --rate krf_f=7.8e-10   # KrF* + F at Kannari's analog
    python3 LIGHT/run_v0_5.py --gap 1                # a 1-cm gap: the charges drift to the walls
    python3 LIGHT/run_v0_5.py --em-E 1500 --em-f 95  # an EM wave, 1500 V/cm rms at 95 GHz (Turnbull IB)

The given power is deposited uniformly and directly into the electron
gas; an EM wave (--em-E, --em-f) is a second input whose absorbed power
is an output. The electrons default to the kinetic EEDF, and a run then
takes minutes, so --bands, --scan and --omissions take an hour or more.
--eedf two_t runs v0_2's closure, in seconds. The optics block prints g0
and alpha at 248 nm (optics.py).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time

from .model import Model, Params, run
from .network import BANDS, BASIC_GROUPS, GROUPS
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
    ("emission_all", "emit", "{:5.1f}"), ("g0_pct_cm", "g0", "{:6.3f}"),
    ("alpha_pct_cm", "alpha", "{:6.3f}"),
    ("f2_left_pct", "F2left", "{:6.1f}"),
    ("tg_end_K", "Tg_end", "{:6.0f}"), ("closure", "closure", "{:8.0e}"),
)


def params_from_args(a) -> Params:
    kw = dict(network=a.network, pressure_atm=a.pressure,
              x_f2=a.f2_pct / 100.0,
              power_W_cm3=a.power * 1e3, pulse_s=a.pulse * 1e-6,
              t_end_s=a.t_end * 1e-6, ne_cm3=a.ne, ne_mode=a.ne_mode,
              tau_sp_s=a.tau_sp * 1e-9, rtol=a.rtol,
              trimer_tg_pow=a.trimer_tg_pow,
              rates=tuple(_rate(x) for x in a.rate),
              groups=tuple(g for g in a.add.split(",") if g),
              eedf=a.eedf, ladder=a.ladder, k76=a.k76, gap_cm=a.gap,
              em_E_V_cm=a.em_E, em_f_Hz=a.em_f * 1e9)
    return Params(**kw)


def _rate(text):
    name, _, value = text.partition("=")
    if name not in BANDS or not value:
        raise SystemExit(f"--rate {text!r}: give NAME=VALUE with NAME one "
                         f"of {', '.join(BANDS)}")
    return name, float(value)


def with_band(p: Params, name, value) -> Params:
    """p with one band set: the trimer law through its own field, every
    other band through Params.rates."""
    if name == "trimer_tg_pow":
        return p.replace(trimer_tg_pow=value)
    return p.replace(rates=tuple(kv for kv in p.rates if kv[0] != name)
                     + ((name, value),))


def do_bands(base: Params):
    """The network as given, then each band's two ends, one band at a time
    (network.BANDS). An end equal to the band's value in `base` is that
    run again and is not repeated."""
    if base.network != "basic":
        raise SystemExit("--bands needs the basic network, the one that "
                         "carries the bands")
    ref = row_dict(run(base))
    rows = [dict(ref, band="as given", end="", value=float("nan"))]
    for b in BANDS.values():
        now = (base.trimer_tg_pow if b.name == "trimer_tg_pow"
               else dict(base.rates).get(b.name))
        now = b.default if now is None else now
        for end, v in (("low", b.lo), ("high", b.hi)):
            d = ref if v == now else row_dict(run(with_band(base, b.name,
                                                             v)))
            rows.append(dict(d, band=b.name, end=end, value=v))
    return rows


BAND_COLS = (
    ("band", "band", "{:>13s}"), ("end", "end", "{:>4s}"),
    ("value", "value", "{:9.3g}"), ("emission_248", "248nm", "{:6.2f}"),
    ("d248", "change", "{:+6.2f}"), ("emission_vuv", "VUV", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("g0_pct_cm", "g0", "{:6.3f}"),
    ("alpha_pct_cm", "alpha", "{:6.3f}"), ("tg_end_K", "Tg_end", "{:6.0f}"),
)


def band_table(rows) -> str:
    def w(fmt):
        return len(fmt.format("x" if "s}" in fmt else 0.0))
    head = "  ".join(f"{h:>{w(fmt)}s}" for _, h, fmt in BAND_COLS)
    out = [head]
    ref = rows[0]["emission_248"]
    for r in rows:
        r = dict(r, d248=r["emission_248"] - ref)
        out.append("  ".join(fmt.format(r[k]) for k, _, fmt in BAND_COLS))
    out.append("")
    for b in BANDS.values():
        out.append(f"  {b.name:13s} {b.what}, {b.lo:.3g} to {b.hi:.3g} "
                   f"{b.unit}, default {b.default:.3g}")
    return "\n".join(out)


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
    ap = argparse.ArgumentParser(prog="run_v0_5.py", description=__doc__,
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
    ap.add_argument("--rate", action="append", default=[],
                    metavar="NAME=VALUE", help="set one band (repeatable), "
                    "sensitivity only: " + ", ".join(
                        b for b in BANDS if b != "trimer_tg_pow"))
    ap.add_argument("--bands", action="store_true",
                    help="run each band's two ends, one band at a time")
    ap.add_argument("--tau-sp", type=float, default=1e9 / 1.43e8,
                    help="KrF(B) radiative lifetime, ns (notes: ~4)")
    ap.add_argument("--network", choices=("basic", "v0"), default="basic",
                    help="basic: the KrF(B,C) network (default); v0: the "
                    "notes' six reactions as written")
    ap.add_argument("--add", default="",
                    help="optional groups, comma list: " + ",".join(GROUPS)
                    + " (basic: " + ",".join(BASIC_GROUPS) + ")"
                    + " or all")
    ap.add_argument("--eedf", choices=("auto", "kinetic", "two_t",
                                       "maxwellian"), default="auto",
                    help="the electrons: kinetic, the two-term Boltzmann EEDF "
                    "solved during the run (v0_3, the basic network's "
                    "auto); two_t, v0_2's Maxwellian with the tabulated "
                    "tail; maxwellian, the notes'")
    ap.add_argument("--k76", choices=("jung", "kannari"), default="jung",
                    help="v0_5: the 5s -> 5p excitation on Jung et al.'s "
                    "measured cross sections (default) or v0_4's Kannari fit")
    ap.add_argument("--ladder", choices=("off", "kannari"), default="off",
                    help="the 5p rung's bands: off (the unmeasured rows "
                    "at 0, their default) or kannari (each at Kannari "
                    "1985's printed value), a sensitivity preset")
    ap.add_argument("--gap", type=float, default=None,
                    help="electrode gap, cm: the charges drift and diffuse "
                    "to the walls (default: no walls)")
    ap.add_argument("--em-E", type=float, default=0.0,
                    help="EM wave, rms field inside the plasma, V/cm "
                    "(default 0: no wave)")
    ap.add_argument("--em-f", type=float, default=95.0,
                    help="EM wave frequency, GHz")
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
    if a.bands:
        rows = do_bands(base)
        print(band_table(rows))
        print(f"\n(percent of deposited energy; change = 248 nm against the "
              f"first row; {len(rows)} rows in {time.time() - t0:.1f} s)")
        if a.csv:
            with open(a.csv, "w", newline="") as fh:
                wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
                wr.writeheader()
                wr.writerows(rows)
        return 0
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
