"""Tahoma under the prototype beam drive, for any LIGHT version, with the
248-nm optics (v0_4): g0, alpha0 by species and g0 - alpha0 over the
200-400 ns window against the measured 6.39, -1.11 and 5.28 %/cm (Spectra
Technology 1989 Table 3-3, beam stop in, 330 kW/cm^3), and the pure-Kr
loss (no F2) against the measured -0.67 %/cm at 257.3 nm (Table 3-2; the
cross sections here are 248-nm values, so that one is a rough check).

A score only, never a fit target (the repo's predictive-only policy).

    LIGHT_VERSION=v0_4 python3 LIGHT/docs/krf_data/tahoma_v0_4.py [--kinetic]
    LIGHT_VERSION=v0_3 python3 LIGHT/docs/krf_data/tahoma_v0_4.py

Runs ML's Tahoma configuration (ML/ml/forward.py tahoma_params) under the
survey's BeamModel with the version's own basic network (no constants
imposed), on v0_2's closure by default (seconds) or the kinetic one
(--kinetic, minutes). Versions before v0_4 have no optics module, so for
them alpha is computed here from v0_4's cross sections on their densities.
"""

import argparse
import importlib
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIGHT = HERE.parents[1]
KM = LIGHT.parent
sys.path.insert(0, str(KM))
sys.path.insert(0, str(LIGHT))

VERSION = os.environ.get("LIGHT_VERSION", "v0_4")
MEASURED = {"g0": 6.39, "alpha0": 1.11, "net": 5.28, "pure_kr_257": 0.67}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kinetic", action="store_true",
                    help="the kinetic EEDF (minutes) instead of v0_2's "
                    "closure")
    ap.add_argument("--ladder", default="off", help="v0_4: off | kannari")
    ap.add_argument("--rate", action="append", default=[],
                    metavar="NAME=VALUE", help="v0_4 bands")
    ap.add_argument("--json", default=None,
                    help="also write the numbers to this file")
    a = ap.parse_args(argv)
    from ML.ml import forward                           # binds LIGHT_VERSION
    assert forward.LIGHT_VERSION == VERSION, forward.LIGHT_VERSION
    # the version's own optics module (v0_4 and later carry one; older
    # versions take v0_4's cross sections on their densities)
    try:
        optics = importlib.import_module(f"{VERSION}.optics")
    except ModuleNotFoundError:
        optics = importlib.import_module("v0_4.optics")
    closure = None if a.kinetic else {"eedf": "two_t"}
    p = forward.tahoma_params(closure=closure)
    extra = {}
    if VERSION == "v0_4":
        extra["ladder"] = a.ladder
        if a.rate:
            extra["rates"] = tuple((k, float(v)) for k, v in
                                   (x.split("=") for x in a.rate))
    p = p.replace(**extra)
    drive = forward._drive()
    out = {}
    for tag, pp in (("Kr/F2", p), ("pure Kr", p.replace(x_f2=1e-9))):
        t0 = time.time()
        r = drive.BeamModel(pp).run()
        t = r.t
        win = (t >= forward.WINDOW[0]) & (t <= forward.WINDOW[1])
        g0 = 100.0 * optics.SIGMA_SE_CM2 * float(np.mean(r.n["KrF(B)"][win]))
        al = {s: 100.0 * sig * float(np.mean(r.n[s][win]))
              for s, (sig, _) in optics.ABSORBERS.items() if s in r.n}
        alpha = sum(al.values())
        out[tag] = dict(g0=g0, alpha=alpha, by=al, ne=float(np.mean(r.ne[win])),
                        te=float(np.mean(r.te[win])), tg=float(r.tg[-1]),
                        secs=time.time() - t0)
    o = out["Kr/F2"]
    print(f"LIGHT {VERSION} at Tahoma ({'kinetic' if a.kinetic else 'two_t'}"
          f" closure, ladder {a.ladder}, rates {a.rate or 'defaults'}), "
          f"200-400 ns window:")
    print(f"  g0        {o['g0']:.3f} %/cm   measured {MEASURED['g0']}   "
          f"x{o['g0'] / MEASURED['g0']:.3f}")
    print(f"  alpha0    {o['alpha']:.3f} %/cm   measured {MEASURED['alpha0']}  "
          f"x{o['alpha'] / MEASURED['alpha0']:.3f}")
    print(f"  g0-alpha0 {o['g0'] - o['alpha']:.3f} %/cm   measured "
          f"{MEASURED['net']}   x{(o['g0'] - o['alpha']) / MEASURED['net']:.3f}")
    parts = ", ".join(f"{s} {v:.3f}" for s, v in
                      sorted(o["by"].items(), key=lambda kv: -kv[1]))
    print(f"  alpha by species: {parts} %/cm")
    print(f"  n_e {o['ne']:.3g} cm^-3, Te {o['te']:.2f} eV, Tg end "
          f"{o['tg']:.0f} K  ({o['secs']:.0f} s)")
    q = out["pure Kr"]
    parts = ", ".join(f"{s} {v:.3f}" for s, v in
                      sorted(q["by"].items(), key=lambda kv: -kv[1]))
    print(f"  pure Kr, no F2: alpha {q['alpha']:.3f} %/cm (248-nm cross "
          f"sections) against the measured {MEASURED['pure_kr_257']} at "
          f"257.3 nm: {parts}")
    if a.json:
        import json
        Path(a.json).write_text(json.dumps(
            {"version": VERSION, "kinetic": a.kinetic, "ladder": a.ladder,
             "rates": a.rate, "measured": MEASURED, **out}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
