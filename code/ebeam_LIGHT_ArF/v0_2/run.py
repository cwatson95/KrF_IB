"""Command line for LIGHT_ArF_ebeam v0_2 (launch through run_v0_2.py).

    python3 LIGHT_ArF_ebeam/run_v0_2.py                  # the reference point
    python3 LIGHT_ArF_ebeam/run_v0_2.py --drive beam     # the same, e-beam
    python3 LIGHT_ArF_ebeam/run_v0_2.py --bands          # each band's two ends
    python3 LIGHT_ArF_ebeam/run_v0_2.py --rate ar2_f=3e-10   # one band value
    python3 LIGHT_ArF_ebeam/run_v0_2.py --pressure 2 --power 500 --f2-pct 0.3
    python3 LIGHT_ArF_ebeam/run_v0_2.py --eedf maxwellian    # no kinetic tail
    python3 LIGHT_ArF_ebeam/run_v0_2.py --omissions      # with/without 'electron'
    python3 LIGHT_ArF_ebeam/run_v0_2.py --scan --csv scan.csv

The Ar/F2 ArF model. The given power is deposited uniformly, either
directly into the electron gas (--drive power, LIGHT's own input) or by an
electron beam (--drive beam, v0_2/beam.py).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time

from .model import Params, run
from .network import ARF_TAU_B, BANDS, BASIC_GROUPS
from .report import row_dict, summary

SCAN_PRESSURE = (0.5, 1.0, 2.0, 3.0)                # atm
SCAN_POWER = (100.0, 300.0, 500.0)                  # kW/cm^3

TABLE_COLS = (
    ("pressure_atm", "p/atm", "{:5.2f}"), ("power_kW_cm3", "kW/cc", "{:5.0f}"),
    ("x_f2", "F2", "{:6.4f}"), ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Arstar", "Ar*", "{:5.1f}"), ("to_Arplus", "Ar+", "{:5.1f}"),
    ("to_ArFB", "ArF(B)", "{:6.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_193", "193nm", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("f2_left_pct", "F2left", "{:6.1f}"),
    ("tg_end_K", "Tg_end", "{:6.0f}"), ("closure", "closure", "{:8.0e}"),
)


def params_from_args(a) -> Params:
    kw = dict(pressure_atm=a.pressure, x_f2=a.f2_pct / 100.0,
              power_W_cm3=a.power * 1e3, drive=a.drive,
              beam_route=a.beam_route, beam_sharing=a.beam_sharing,
              pulse_s=a.pulse * 1e-6,
              t_end_s=a.t_end * 1e-6, ne_cm3=a.ne, eedf=a.eedf,
              tau_sp_s=a.tau_sp * 1e-9, rtol=a.rtol,
              trimer_tg_pow=a.trimer_tg_pow,
              rates=tuple(_rate(x) for x in a.rate),
              groups=tuple(g for g in a.add.split(",") if g))
    return Params(**kw)


def _rate(text):
    name, _, value = text.partition("=")
    if name not in BANDS or name == "trimer_tg_pow" or not value:
        raise SystemExit(f"--rate {text!r}: give NAME=VALUE with NAME one "
                         f"of {', '.join(b for b in BANDS if b != 'trimer_tg_pow')}"
                         " (the trimer law is --trimer-tg-pow)")
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
    ("value", "value", "{:9.3g}"), ("emission_193", "193nm", "{:6.2f}"),
    ("d193", "change", "{:+6.2f}"), ("emission_vuv", "VUV", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("tg_end_K", "Tg_end", "{:6.0f}"),
)


def band_table(rows) -> str:
    def w(fmt):
        return len(fmt.format("x" if "s}" in fmt else 0.0))
    head = "  ".join(f"{h:>{w(fmt)}s}" for _, h, fmt in BAND_COLS)
    out = [head]
    ref = rows[0]["emission_193"]
    for r in rows:
        r = dict(r, d193=r["emission_193"] - ref)
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
    """A pressure x power grid at the given F2 fraction (n_e is free, so
    the seed is not scanned)."""
    rows = []
    for p in SCAN_PRESSURE:
        for pw in SCAN_POWER:
            res = run(base.replace(pressure_atm=p, power_W_cm3=pw * 1e3))
            rows.append(row_dict(res))
    return rows


def do_omissions(base: Params):
    rows = []
    for g in ((),) + tuple((x,) for x in BASIC_GROUPS):
        rows.append(row_dict(run(base.replace(groups=g))))
    return rows


OMIT_COLS = (
    ("groups", "network", "{:>12s}"), ("te_mean_eV", "Te", "{:4.2f}"),
    ("to_Arstar", "Ar*", "{:5.1f}"), ("to_Arplus", "Ar+", "{:5.1f}"),
    ("to_ArFB", "ArF(B)", "{:6.1f}"), ("heat_reactions", "heat", "{:5.1f}"),
    ("elastic", "elast", "{:5.1f}"), ("emission_193", "193nm", "{:5.1f}"),
    ("emission_all", "emit", "{:5.1f}"), ("arf_to_193_pct", "B->193",
                                          "{:6.1f}"),
    ("arstar_to_arf_pct", "Ar*->B", "{:6.1f}"),
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
    ap = argparse.ArgumentParser(prog="run_v0_2.py", description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pressure", type=float, default=1.0, help="atm")
    ap.add_argument("--f2-pct", type=float, default=0.5,
                    help="F2 mole percent, rest Ar")
    ap.add_argument("--power", type=float, default=300.0, help="kW/cm^3")
    ap.add_argument("--drive", choices=("power", "beam"), default="power",
                    help="power: into the electron gas (LIGHT's own input); "
                    "beam: an electron beam (v0_2/beam.py)")
    ap.add_argument("--beam-route", choices=("asymptotic", "integral"),
                    default="asymptotic", help="the beam's yields per pair "
                    "(integral + equal = REBUILD's own, a sensitivity)")
    ap.add_argument("--beam-sharing", choices=("opb", "equal"),
                    default="opb", help="the secondary's energy: "
                    "Opal-Peterson-Beaty or equal sharing (REBUILD)")
    ap.add_argument("--pulse", type=float, default=1.0, help="us")
    ap.add_argument("--t-end", type=float, default=2.0, help="us")
    ap.add_argument("--ne", type=float, default=1e15,
                    help="cm^-3, the initial electron density (a seed; "
                    ">= 1e13 at these powers)")
    ap.add_argument("--eedf", choices=("auto", "two_t", "maxwellian"),
                    default="auto", help="auto: the two-temperature tail")
    ap.add_argument("--trimer-tg-pow", type=float, default=None,
                    help="sensitivity only: the trimer's (Tg/300)^n "
                    "exponent (band trimer_tg_pow)")
    ap.add_argument("--rate", action="append", default=[],
                    metavar="NAME=VALUE", help="set one band (repeatable), "
                    "sensitivity only: " + ", ".join(
                        b for b in BANDS if b != "trimer_tg_pow"))
    ap.add_argument("--bands", action="store_true",
                    help="run each band's two ends, one band at a time")
    ap.add_argument("--tau-sp", type=float, default=ARF_TAU_B * 1e9,
                    help="ArF(B) radiative lifetime, ns")
    ap.add_argument("--add", default="",
                    help="optional groups, comma list: "
                    + ",".join(BASIC_GROUPS) + " or all")
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--scan", action="store_true",
                    help="pressure x power at the given F2")
    ap.add_argument("--omissions", action="store_true",
                    help="the network vs the network + each optional group")
    ap.add_argument("--csv", default=None, help="write scan/omissions rows")
    ap.add_argument("--json", default=None, help="write the single-run rows")
    a = ap.parse_args(argv)
    base = params_from_args(a)
    t0 = time.time()
    if a.bands:
        rows = do_bands(base)
        print(band_table(rows))
        print(f"\n(percent of deposited energy; change = 193 nm against the "
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
