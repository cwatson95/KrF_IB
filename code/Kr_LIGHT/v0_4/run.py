"""Command line for LIGHT_Kr v0_4, the pure-krypton model (launch through
LIGHT_Kr/run_v0_4.py).

    python3 LIGHT_Kr/run_v0_4.py                        # the reference point
    python3 LIGHT_Kr/run_v0_4.py --eedf two_t           # LIGHT v0_2's closure
    python3 LIGHT_Kr/run_v0_4.py --bands                # each band's two ends
    python3 LIGHT_Kr/run_v0_4.py --rate kr2_pen=1e-10   # one band value
    python3 LIGHT_Kr/run_v0_4.py --pressure 2 --power 500 --ne 1e16
    python3 LIGHT_Kr/run_v0_4.py --scan --csv scan.csv  # pressure x power
    python3 LIGHT_Kr/run_v0_4.py --eedf two_t --tail beam   # v0_3's tail table
    python3 LIGHT_Kr/run_v0_4.py --afterglow-memory     # the EEDF stepped after the pulse

    python3 LIGHT_Kr/run_v0_4.py --ladder kannari       # every 5p-rung band at Kannari's value
    python3 LIGHT_Kr/run_v0_4.py --k76 kannari          # LIGHT v0_4's 5s -> 5p fit
    python3 LIGHT_Kr/run_v0_4.py --gap 1                # a 1-cm gap: the charges drift to the walls
    python3 LIGHT_Kr/run_v0_4.py --em-E 1500 --em-f 95  # an EM wave, 1500 V/cm rms at 95 GHz (Turnbull IB)

v0_4 carries no optional group (v0_1 to v0_3's 'electron' booked the
lumped Kr*), so --add and --omissions have nothing to add.

The fill is pure krypton. The given power is deposited uniformly and
directly into the electron gas; an EM wave (--em-E, --em-f) is a second
input whose absorbed power is an output. The electrons default to the
kinetic EEDF, and a run then takes minutes, so --bands, --scan and
--omissions take an hour or more. --eedf two_t runs LIGHT v0_2's closure,
in seconds.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time

from .model import Params, run
from .network import BANDS, GROUPS
from .report import row_dict, summary

SCAN_PRESSURE = (0.1, 0.25, 0.5, 1.0, 2.0)          # atm  (notes: 0.1-2)
SCAN_POWER = (100.0, 300.0, 500.0)                  # kW/cm^3 (notes)

TABLE_COLS = (
    ("pressure_atm", "p/atm", "{:5.2f}"), ("power_kW_cm3", "kW/cc", "{:5.0f}"),
    ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Krstar", "Kr*", "{:5.1f}"), ("to_Krplus", "Kr+", "{:5.1f}"),
    ("to_Kr2star", "Kr2*", "{:5.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_vuv", "147nm", "{:5.1f}"),
    ("emission_nir", "nearIR", "{:6.2f}"), ("emission_all", "emit", "{:5.1f}"),
    ("ne_end_cm3", "n_e_end", "{:7.1e}"),
    ("tg_end_K", "Tg_end", "{:6.0f}"), ("closure", "closure", "{:8.0e}"),
)


def params_from_args(a) -> Params:
    kw = dict(pressure_atm=a.pressure,
              power_W_cm3=a.power * 1e3, pulse_s=a.pulse * 1e-6,
              t_end_s=a.t_end * 1e-6, ne_cm3=a.ne, rtol=a.rtol,
              rates=tuple(_rate(x) for x in a.rate),
              groups=tuple(g for g in a.add.split(",") if g),
              eedf=a.eedf, ladder=a.ladder, k76=a.k76, gap_cm=a.gap,
              em_E_V_cm=a.em_E, em_f_Hz=a.em_f * 1e9, tail=a.tail,
              kin_afterglow_memory=a.afterglow_memory)
    return Params(**kw)


def _rate(text):
    name, _, value = text.partition("=")
    if name not in BANDS or not value:
        raise SystemExit(f"--rate {text!r}: give NAME=VALUE with NAME one "
                         f"of {', '.join(BANDS)}")
    return name, float(value)


def with_band(p: Params, name, value) -> Params:
    """p with one band set through Params.rates."""
    return p.replace(rates=tuple(kv for kv in p.rates if kv[0] != name)
                     + ((name, value),))


def do_bands(base: Params):
    """The network as given, then each band's two ends, one band at a time
    (network.BANDS). An end equal to the band's value in `base` is that
    run again and is not repeated."""
    ref = row_dict(run(base))
    rows = [dict(ref, band="as given", end="", value=float("nan"))]
    for b in BANDS.values():
        now = dict(base.rates).get(b.name)
        now = b.default if now is None else now
        for end, v in (("low", b.lo), ("high", b.hi)):
            d = ref if v == now else row_dict(run(with_band(base, b.name,
                                                             v)))
            rows.append(dict(d, band=b.name, end=end, value=v))
    return rows


BAND_COLS = (
    ("band", "band", "{:>9s}"), ("end", "end", "{:>4s}"),
    ("value", "value", "{:9.3g}"), ("emission_vuv", "147nm", "{:6.2f}"),
    ("dvuv", "change", "{:+6.2f}"), ("emission_nir", "nearIR", "{:6.2f}"),
    ("emission_all", "emit", "{:5.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("tg_end_K", "Tg_end", "{:6.0f}"),
)


def band_table(rows) -> str:
    def w(fmt):
        return len(fmt.format("x" if "s}" in fmt else 0.0))
    head = "  ".join(f"{h:>{w(fmt)}s}" for _, h, fmt in BAND_COLS)
    out = [head]
    ref = rows[0]["emission_vuv"]
    for r in rows:
        r = dict(r, dvuv=r["emission_vuv"] - ref)
        out.append("  ".join(fmt.format(r[k]) for k, _, fmt in BAND_COLS))
    out.append("")
    for b in BANDS.values():
        out.append(f"  {b.name:9s} {b.what}, {b.lo:.3g} to {b.hi:.3g} "
                   f"{b.unit}, default {b.default:.3g}")
    return "\n".join(out)


def table(rows) -> str:
    head = "  ".join(f"{h:>{len(fmt.format(0))}s}" for _, h, fmt in TABLE_COLS)
    out = [head]
    for r in rows:
        out.append("  ".join(fmt.format(r[k]) for k, _, fmt in TABLE_COLS))
    return "\n".join(out)


def do_scan(base: Params):
    """The notes' pressure x power grid. n_e is solved, so the seed is not
    an axis."""
    rows = []
    for p in SCAN_PRESSURE:
        for pw in SCAN_POWER:
            res = run(base.replace(pressure_atm=p, power_W_cm3=pw * 1e3))
            rows.append(row_dict(res))
    return rows


def do_omissions(base: Params):
    """The network alone, with each optional group, and with all of them."""
    rows = []
    for g in ((),) + tuple((x,) for x in GROUPS) + (("all",),):
        rows.append(row_dict(run(base.replace(groups=g))))
    return rows


OMIT_COLS = (
    ("groups", "network", "{:>12s}"), ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Krstar", "Kr*", "{:5.1f}"), ("to_Krplus", "Kr+", "{:5.1f}"),
    ("to_Kr2star", "Kr2*", "{:5.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_vuv", "147nm", "{:5.1f}"),
    ("emission_nir", "nearIR", "{:6.2f}"), ("emission_all", "emit", "{:5.1f}"),
    ("krstar_to_kr2_pct", "Kr*->Kr2*", "{:9.1f}"),
    ("kr2_to_vuv_pct", "Kr2*->hv", "{:8.1f}"), ("closure", "closure", "{:8.0e}"),
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
    ap = argparse.ArgumentParser(prog="run_v0_4.py", description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pressure", type=float, default=1.0,
                    help="atm of pure Kr")
    ap.add_argument("--power", type=float, default=300.0, help="kW/cm^3")
    ap.add_argument("--pulse", type=float, default=1.0, help="us")
    ap.add_argument("--t-end", type=float, default=2.0, help="us")
    ap.add_argument("--ne", type=float, default=1e15,
                    help="cm^-3, the initial electron density (a seed)")
    ap.add_argument("--rate", action="append", default=[],
                    metavar="NAME=VALUE", help="set one band (repeatable), "
                    "sensitivity only: " + ", ".join(BANDS))
    ap.add_argument("--bands", action="store_true",
                    help="run each band's two ends, one band at a time")
    ap.add_argument("--add", default="",
                    help="optional groups, comma list: " + ",".join(GROUPS)
                    + " or all")
    ap.add_argument("--eedf", choices=("auto", "kinetic", "two_t",
                                       "maxwellian"), default="auto",
                    help="the electrons: kinetic, the two-term Boltzmann EEDF "
                    "solved during the run (auto); two_t, LIGHT v0_2's "
                    "Maxwellian with the tabulated tail; maxwellian, the "
                    "notes'")
    ap.add_argument("--tail", choices=("field", "beam"), default="field",
                    help="the two-temperature closure's tail table: field "
                    "(v0_4's, field-driven pure Kr) or beam (v0_1-v0_3's)")
    ap.add_argument("--afterglow-memory", action="store_true",
                    help="time-step the kinetic EEDF after the pulse too")
    ap.add_argument("--k76", choices=("jung", "kannari"), default="jung",
                    help="the 5s -> 5p excitation on Jung et al.'s measured "
                    "cross sections (default) or LIGHT v0_4's Kannari fit")
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
                    help="pressure x power over the notes' ranges")
    ap.add_argument("--omissions", action="store_true",
                    help="the network vs the network + each optional group "
                    "vs + all")
    ap.add_argument("--csv", default=None, help="write scan/omissions rows")
    ap.add_argument("--json", default=None, help="write the single-run rows")
    a = ap.parse_args(argv)
    base = params_from_args(a)
    t0 = time.time()
    if a.bands:
        rows = do_bands(base)
        print(band_table(rows))
        print(f"\n(percent of deposited energy; change = 147 nm against the "
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
