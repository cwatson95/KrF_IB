"""SCORE.md from score_0d.json and score_1d.json (score_arf.py --report),
beside LIGHT_ArF's record (score_0d_LIGHT_ArF.json, kept as it was)."""
from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TITLES = {
    "electra_formation": "Electra formation efficiency, g0 Isat / P (0.82 atm, 0.3 % F2)",
    "electra_yield_band": "Sensitivity: Electra formation efficiency with the ArF* yield at 1.0",
    "electra_gain_trend": "Electra gain against pressure (single pass), ratio to 0.82 atm",
    "electra_gain_pump": "Electra gain, 1.09 against 0.395 MW/cc at 0.82 atm",
    "keio_formation": "Keio formation efficiency, g0 Is / P (1 atm, 0.4 % F2)",
    "sandia_formation": "Sandia formation efficiency, g0 Is / P (963 Torr)",
    "rokni_vs_ar": "Rokni 1977 Fig. 4: 193-nm fluorescence against Ar (norm. at 1 atm)",
    "rokni_current_ratio": "Rokni 1977 Fig. 5: fluorescence ratio for a x25 current change",
    "jancaitis_scan": "Jancaitis Fig. 3: peak fluorescence against current density (norm. near 74 A/cm2)",
    "sri_vs_ar": "SRI MP-76-99 Fig. 4: photons per initial excitation against Ar (1 Torr F2)",
    "sri_vs_f2": "SRI MP-76-99 Fig. 11: 193-nm photons per initial excitation against F2",
    "morozov_vuv": "Morozov 2008: pure-Ar VUV efficiency against pressure",
    "kinetic_decays": "Decay rates (consistency: the constants' own sources)",
    "stern_volmer": "Stern-Volmer plots, 1/I normalised (consistency)",
    "electra_formation_maxwellian": "Sensitivity: Electra formation efficiency, own input on the Maxwellian EEDF",
    "keio_formation_maxwellian": "Sensitivity: Keio formation efficiency, own input on the Maxwellian EEDF",
    "sandia_formation_maxwellian": "Sensitivity: Sandia formation efficiency, own input on the Maxwellian EEDF",
    "rokni_vs_ar_maxwellian": "Sensitivity: Rokni Fig. 4 against Ar, own input on the Maxwellian EEDF",
    "rokni_current_ratio_maxwellian": "Sensitivity: Rokni Fig. 5 current ratio, own input on the Maxwellian EEDF",
    "jancaitis_scan_maxwellian": "Sensitivity: Jancaitis Fig. 3 against current, own input on the Maxwellian EEDF",
}
XLABEL = {
    "electra_formation": "MW/cc", "electra_yield_band": "MW/cc",
    "electra_gain_trend": "atm",
    "electra_gain_pump": "pump ratio", "keio_formation": "MW/cm3",
    "sandia_formation": "MW/cm3", "rokni_vs_ar": "Ar atm",
    "rokni_current_ratio": "Ar atm", "jancaitis_scan": "A/cm2",
    "sri_vs_ar": "Ar atm", "sri_vs_f2": "F2 Torr", "morozov_vuv": "hPa",
}


READING = """## Reading (written 2026-09-30 from the numbers below)

- **The e-beam drive improves the model overall.** Over the 11 held-out
  tests (device, shape, yield, pure) the mean |ln gm| is 0.314 for
  LIGHT_ArF on its own input, 0.303 here on the own input (the seeded SRI
  starts take the drive's 0.39 excited atoms per ion in place of 0.85) and
  0.256 under the drive, and the rms of the logs over their 87 points
  0.476, 0.460 and 0.392. LIGHT_ArF's prototype drive gave 0.275 and 0.429.
  Under the beam Te sits at 1.5 to 2.2 eV at every device and shape
  point, the state the tail table was made in, where the own input runs to
  21-94 eV at Rokni's pumps.
- **What it fixes.** Rokni's x25 current ratio goes from x1.67 to x0.97,
  every pressure within x0.92 to x1.01, and his pressure curve's rms from
  0.154 to 0.101. Electra's formation efficiency goes from x0.64 to x0.84
  at 1.09 MW/cc (x0.55 to x0.72 at 0.90), its pump ratio from x1.36 to
  x1.13, and the yield-band sensitivity from x1.19 to x1.08. Jancaitis's
  Ar2F* at the top currents goes from x1.67 to x1.34.
- **What it makes worse.** Sandia (abstract only) goes from x1.39 to x1.82
  and Keio at its printed 0.51 MW/cm3 from x0.86 to x1.14. Jancaitis's
  ArF* peak against current goes from rms 0.28 to 0.35: from 74 to 562
  A/cm2 the model's peak rises x5.5 against x9.1 measured. The drive makes
  ArF(B) at 0.25 to 0.30 of the deposited energy at every device pump,
  where the measured values fall from 0.36 at Electra's 0.90 MW/cc to
  0.16 at Sandia's 5 MW/cm3. Neither input follows that fall. Electra's
  gain drop with pressure stays missed (x1.32 at 1.36 atm on both).
- **Morozov is not reached.** The drive takes pure argon from x2.15 to
  x1.75 (x2.07 to x1.70 at 500-1400 hPa, where his 0.33 +- 0.04 is
  printed). From 300 to 1400 hPa it gives 0.54 to 0.56 (below, the
  pulsed test's 20-us count cuts it), the network's one-photon ceiling: each
  ion and each excitation ends as one 127-nm photon, and with the drive's
  25.2 eV per pair and 0.39 excitations per pair that is 0.538, beside the
  0.53 Lorents computed (Morozov's ref. 17). Run as Morozov ran it, a
  continuous beam, the model is flat at 0.54 from 151 to 1000 hPa where
  he measured 0.23 to 0.33. His argon sits at 62 % of that maximum while
  his Kr and Xe sit at 91 and 98 % of theirs (Eckstrom), so no drive
  closes it: either the network lacks an argon VUV loss that has been
  measured, or his argon value is low.
- **REBUILD's yields, found here.** REBUILD books the yields per pair as
  the equal-sharing Spencer-Fano integral at one source energy (516 eV
  with F2). On the IST-Lisbon set that rule makes them oscillate with the
  source energy (0.31 to 1.28 excitations per ion and W 28 to 35 eV
  between 100 and 900 eV), and the prototype sat on one point of it (W
  30.85 eV, 0.85 excitations per pair). With the measured
  Opal-Peterson-Beaty secondary spectrum the cascade converges to W 25.2
  eV, 4.5 % under the measured 26.4, and 0.39 excitations per ion. The
  device scores move by at most 4 % between the two (Electra x0.84 against
  x0.81), Morozov at 1000 hPa from x1.78 to x1.67.
- **SRI's yields per initial excitation** (seeded, the drive's partition)
  are x0.66 at 193 nm and x0.90 at 290 nm against Ar (x0.63 and x0.86 in
  LIGHT_ArF), and x0.84 against F2 (x0.78). At 290 nm the model still
  gives too little at 1 atm (x0.35) and too much at high pressure (x1.4
  at 9 atm), the missing Ar2F* loss LIGHT_ArF found.
- **1Dv0_2 runs both inputs** (64 cells). Under the beam its centre cell
  gives v0_2's beam numbers at every device point to 1.1e-4 (3.9e-4 on the
  own input, both at Sandia), and its gap average lies 2.1 to 2.4 % lower
  (3.1 % on the own input). The own input's exception at Rokni's 2 and 4
  atm, where v0_2's Te passes Fluid1D's 40-eV ceiling and the 1-D centre
  departs (x0.90 and x0.94), is gone under the beam. There v0_2's Te stays
  at 1.5 to 2.2 eV, and the 1-D centre follows v0_2's beam runs at the
  1-D's pressures (0.243 to 4.003 atm) to 1e-4. The own input's 1-D runs
  are unchanged by the drive, bit for bit.
- **Kinetic consistency** is unchanged: the decays their constants came
  from at gm 0.975 (rms 0.12 over 69 points) and the Stern-Volmer plots at
  gm 1.06 (rms 0.17), Chen & Payne's seeded ones now on the drive's
  partition (Fig. 16 x0.88, Fig. 17 x0.98)."""


def gm(xs):
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def rms_log(xs):
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    return math.sqrt(sum(math.log(x) ** 2 for x in xs) / len(xs)) if xs else float("nan")


def fmt(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return ""
    a = abs(v)
    if a == 0:
        return "0"
    if a >= 1e4 or a < 1e-3:
        return f"{v:.3g}"
    return f"{v:.3f}" if a < 10 else f"{v:.2f}"


def load(name):
    p = HERE / name
    return json.loads(p.read_text()) if p.is_file() else None


def ratios(points, key):
    """model/measured, leaving out the normalisation point (ratio 1 by
    construction on both sides)."""
    out = []
    for p in points:
        if key not in p or p.get(key) is None:
            continue
        if abs(p["measured"] - 1.0) < 1e-12 and abs(p[key] - 1.0) < 1e-12:
            continue
        out.append(p[key] / p["measured"])
    return out


def write():
    d0, d1 = load("score_0d.json"), load("score_1d.json")
    old = (load("score_0d_LIGHT_ArF.json") or {}).get("tests", {})
    t1 = (d1 or {}).get("tests", {})
    L = ["# LIGHT_ArF_ebeam against ArF measurements", "",
         "Generated by `score_arf.py --report` from `score_0d.json` "
         "(LIGHT_ArF_ebeam v0_2, 0-D) and `score_1d.json` (1Dv0_2 on "
         f"Fluid1D, {(d1 or {}).get('cells', '-')} cells, 1-cm gap). Every "
         "model value uses the network's defaults (bands at their defaults). "
         "'own' is LIGHT's own input, a given power into the electron gas, "
         "and 'beam' is the model's electron-beam drive (v0_2 drive='beam', "
         "`v0_2/beam.py`: REBUILD's Spencer-Fano cascade on the model's own "
         "cross sections with the Opal-Peterson-Beaty secondaries, W 25.2 eV "
         "and 0.39 excited atoms per ion pair). The 'LIGHT_ArF' columns are "
         "the record this copy started from (`score_0d_LIGHT_ArF.json`): its "
         "own input and its prototype beam drive (REBUILD's equal-sharing "
         "integral yields, W 30.85 eV and 0.85 excited atoms per pair). The "
         "own input differs from LIGHT_ArF's only in the seeded starts "
         "(SRI, Chen & Payne: the drive's 0.39 excited atoms per ion in place "
         "of 0.85) and elsewhere agrees with it to 5e-8. The 1-D "
         "model runs both inputs: 'centre' is its centre cell and 'gap' the "
         "gap average (wall cells included) on the own input, 'centre beam' "
         "and 'gap beam' the same under the electron beam.",
         "", "Ratios are model / measured. For normalised curves the "
         "normalisation point is left out of the statistics. gm = geometric "
         "mean of the ratios, rms = rms of their natural logs.", ""]
    # --- summary --------------------------------------------------------
    L += ["## Summary", "",
          "| test | kind | n | own gm | own rms | beam gm | beam rms | LIGHT_ArF own gm | LIGHT_ArF prototype beam gm | 1-D centre gm | 1-D centre rms | 1-D gap gm | 1-D centre beam gm | 1-D gap beam gm |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    order = [k for k in TITLES if k in d0["tests"]] + [
        k for k in d0["tests"] if k not in TITLES]
    tests = {k: d0["tests"][k] for k in order}
    for name, t in tests.items():
        own = ratios(t["own"], "model")
        beam = ratios(t.get("beam", []), "model") if "beam" in t else []
        o_own = ratios(old.get(name, {}).get("own", []), "model")
        o_beam = ratios(old.get(name, {}).get("beam", []), "model")
        c1 = ratios(t1.get(name, []), "centre")
        g1 = ratios(t1.get(name, []), "gap")
        cb1 = ratios(t1.get(name, []), "centre_beam")
        gb1 = ratios(t1.get(name, []), "gap_beam")
        L.append(f"| {name} | {t['kind']} | {len(own)} | {fmt(gm(own))} | "
                 f"{fmt(rms_log(own))} | {fmt(gm(beam)) if beam else ''} | "
                 f"{fmt(rms_log(beam)) if beam else ''} | "
                 f"{fmt(gm(o_own)) if o_own else ''} | "
                 f"{fmt(gm(o_beam)) if o_beam else ''} | "
                 f"{fmt(gm(c1)) if c1 else ''} | "
                 f"{fmt(rms_log(c1)) if c1 else ''} | "
                 f"{fmt(gm(g1)) if g1 else ''} | "
                 f"{fmt(gm(cb1)) if cb1 else ''} | "
                 f"{fmt(gm(gb1)) if gb1 else ''} |")
    L += ["", READING, ""]
    # --- each test ------------------------------------------------------
    for name, t in tests.items():
        L += [f"## {TITLES.get(name, name)}", "", t["doc"], ""]
        if t["kind"] == "kinetic":
            sets = {}
            for p in t["own"]:
                sets.setdefault(p["set"], []).append(p)
            L += ["| set | n | gm model/measured | rms log |", "|---|---|---|---|"]
            for s, pts in sets.items():
                r = [p["model"] / p["measured"] for p in pts]
                if "1/I" in s:
                    r = ratios(pts, "model")
                L.append(f"| {s} | {len(pts)} | {fmt(gm(r))} | {fmt(rms_log(r))} |")
            L.append("")
            continue
        extra = ""
        base = name.removesuffix("_maxwellian")       # a sensitivity entry
        if base == "jancaitis_scan":
            extra = "species"
        elif base == "sri_vs_ar":
            extra = "band"
        elif base == "sri_vs_f2":
            extra = "ar_atm"
        has_beam = "beam" in t
        has_te = "te" in t["own"][0]
        has_low = "te_low" in t["own"][0]
        pts1 = {round(p["x"], 2): p for p in t1.get(name, [])}
        beam1 = any("centre_beam" in q for q in pts1.values())
        head = f"| {XLABEL.get(base, 'x')} |" + (f" {extra} |" if extra else "") \
            + " measured | own | own/meas |" + (" own Te eV |" if has_te else "") \
            + (" own Te x1/25 eV |" if has_low else "") \
            + (" beam | beam/meas |" if has_beam else "") \
            + (" beam Te eV |" if has_beam and has_te else "") \
            + (" beam Te x1/25 eV |" if has_beam and has_low else "") \
            + (" 1-D centre | 1-D gap |" if pts1 else "") \
            + (" 1-D centre beam | 1-D gap beam |" if beam1 else "")
        L += [head, "|" + "---|" * (head.count("|") - 1)]
        beam = t.get("beam", [None] * len(t["own"]))
        for p, b in zip(t["own"], beam):
            row = f"| {fmt(p['x'])} |" + (f" {p.get(extra)} |" if extra else "")
            row += f" {fmt(p['measured'])} | {fmt(p['model'])} | " \
                   f"{fmt(p['model'] / p['measured'])} |"
            if has_te:
                row += f" {fmt(p['te'])} |"
            if has_low:
                row += f" {fmt(p['te_low'])} |"
            if has_beam:
                row += f" {fmt(b['model'])} | {fmt(b['model'] / p['measured'])} |"
                if has_te:
                    row += f" {fmt(b['te'])} |"
                if has_low:
                    row += f" {fmt(b['te_low'])} |"
            if pts1:
                q = pts1.get(round(p["x"], 2))
                row += (f" {fmt(q['centre'])} | {fmt(q['gap'])} |" if q
                        else " | |")
            if beam1:
                q = pts1.get(round(p["x"], 2))
                row += (f" {fmt(q['centre_beam'])} | {fmt(q['gap_beam'])} |"
                        if q else " | |")
            L.append(row)
        if name == "electra_gain_trend":
            L += ["", "g0 with Lee's calculated 2.75e-16 cm2 (own input): "
                  + ", ".join(f"{fmt(p['x'])} atm {fmt(p['g0_calc_pct_cm'])} "
                              f"%/cm (measured {fmt(p['g0_meas_pct_cm'])})"
                              for p in t["own"])]
        if name.startswith(("electra_formation", "electra_yield_band",
                            "keio_formation", "sandia_formation")):
            for col, label in (("own", "own input"), ("beam", "beam drive")):
                if col not in t:
                    continue
                L += ["", f"Model state in the window ({label}): "
                      + "; ".join(f"{fmt(p['x'])}: Te {fmt(p['state']['te'])} "
                                  f"eV, n_e {fmt(p['state']['ne'])} cm-3, "
                                  f"ArF(B) {fmt(p['state']['nB'])} cm-3, F2 "
                                  f"left {fmt(p['state']['f2_left'])}"
                                  for p in t[col])]
        L.append("")
    (HERE / "SCORE.md").write_text("\n".join(L) + "\n")
    print("wrote", HERE / "SCORE.md")


if __name__ == "__main__":
    write()
