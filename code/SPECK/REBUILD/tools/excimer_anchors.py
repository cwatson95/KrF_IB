#!/usr/bin/env python3
"""excimer_anchors.py -- the MEASURED excimer small-signal-gain anchors.

Ported into the engine's own tools (0.1.40) from the 2026-09-04 validation
audit (~/coding/PIC/Validation/scripts/gain_anchors.py + debug_klimek_spec.py),
so the excimer layer can be scored on every engine change instead of once.
Anchors SCORE; they never steer (Known-bad: no calibration to any anchor).

Small-signal gain is the right observable for the excimer layer: it is linear
in the lasing-state population, measured with a probe below threshold, and
free of the cavity knife-edge that makes device energies exponential in it.
The cavity is opened (R_OC 1e-6, single-pass loss 0.9999) as in the
measurements.

A. Electra ArF small-signal gain vs deposition
   Wolford, Myers & Petrova 2023, NRL/6730/MR--2023/1 Sect. 3.3 (193-nm seed,
   Rigrod analysis, deposition by pressure transducer), 99.7 % Ar / 0.3 % F2,
   500 kV, 100-ns flat top.  Single pass, 30-cm cathode: 1.84 / 1.77 / 1.29 %/cm
   at 0.82 / 1.09 / 1.36 atm (0.395 / 0.497 / 0.583 MW/cm3); double pass,
   10-cm cathode: 3.21 %/cm at 0.82 atm, ~1.09 MW/cm3.  The report does not
   say when in the pulse gain was probed, so both the PEAK and the end-of-
   flat-top (115 ns) values are scored.  Deposition is PRINTED, so it is
   matched exactly (the beam's deposited energy is linear in current density
   and independent of the plasma state -- see `deposited_J`).

B. Klimek et al. 1981 KrF gain / absorption grid (IEEE JQE QE-17, 1847), as
   tabulated first-hand in Kannari, Obara & Fujioka 1985 (JAP 57, 4309)
   Tables II (absorption, F2 EXCLUDED by Kannari's stated convention) and III
   (gain): J_eb = 110 / 55 / 20 A/cm2 x 1 / 2 / 3 amagat x 0.2 / 0.6 % F2,
   steady state.  NOT printed there and therefore ASSUMED: Kr 5 % in Ar, and
   the pump -- matched to Kannari's own printed device calibration (660 kV,
   30 kA over 40 x 4 cm2 = 187.5 A/cm2, "excitation rate is 0.75 MW/cm3 atm",
   p.4314), which is a DIFFERENT device than Klimek's.  The absolute level is
   therefore uncertain (~x2, 2026-09-04 audit Sect. 7); the pump-independent
   discriminators are the firm part and are reported separately:
     - g0/alpha per pressure (both scale with the pump), and
     - the pressure exponent d ln g0 / d ln p at fixed J (measured +0.455
       at 110 A/cm2, 0.2 % F2, 1 -> 3 amagat).

Usage (from the Kinetic_Modeling root):
    python3 SPECK/REBUILD/tools/excimer_anchors.py               # both anchors
    python3 SPECK/REBUILD/tools/excimer_anchors.py --quick       # 4 Electra + 4 Klimek corners
    python3 SPECK/REBUILD/tools/excimer_anchors.py --set petrov17_f2=1 --tag p17
Writes rebuild_out/excimer_anchors/<tag>.json.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for base in [_HERE, *_HERE.parents]:
    if (base / "REBUILD" / "MODULES").is_dir():
        ROOT = base
        sys.path.insert(0, str(base))
        break
else:
    sys.exit("REBUILD/MODULES not found above this script")

from REBUILD.MODULES.engine import Device               # noqa: E402
from REBUILD.MODULES.laser import Cavity                # noqa: E402
from REBUILD.MODULES.network import mode_optics         # noqa: E402
from REBUILD.MODULES.params import engine_parameters    # noqa: E402

AMAGAT_TO_ATM_300K = 300.0 / 273.15   # equal number density at the engine's 300 K

# (atm, MW/cm3, measured g0 %/cm, modelled transverse width m, label)
ELECTRA_GAIN = [
    (0.82, 0.395, 1.84, 0.30, "single pass, 30-cm cathode"),
    (1.09, 0.497, 1.77, 0.30, "single pass, 30-cm cathode"),
    (1.36, 0.583, 1.29, 0.30, "single pass, 30-cm cathode"),
    (0.82, 1.09, 3.21, 0.14, "double pass, 10-cm cathode"),
]
ELECTRA_TREND_MEASURED = 1.29 / 1.84          # 0.82 -> 1.36 atm, 30-cm cathode

KLIMEK_J = (110.0, 55.0, 20.0)
KLIMEK_P = (1.0, 2.0, 3.0)
KLIMEK_F2 = (0.002, 0.006)
# Kannari 1985 Tables II / III "Exp." rows (verified first-hand 2026-09-18;
# also the 2026-09-04 verify_targets.py re-extraction), columns (p, F2)
KLIMEK_ABS = {110.0: [0.86, 0.73, 1.8, 1.3, 2.2, 1.6],
              55.0: [0.41, 0.32, 0.91, 0.60, 1.1, 0.89],
              20.0: [0.21, 0.18, 0.42, 0.26, 0.53, 0.38]}
KLIMEK_GAIN = {110.0: [9.1, 8.3, 14.0, 11.0, 15.0, 11.0],
               55.0: [5.1, 5.3, 6.6, 5.1, 6.6, 5.3],
               20.0: [1.9, 1.4, 2.3, 1.6, 2.2, 1.4]}
KANNARI_MW_CM3_PER_ATM_A = 0.75 / 187.5      # p.4314: 0.75 MW/cm3 atm at 187.5 A/cm2
KLIMEK_KR = 0.05                             # ASSUMED (see B above)
KLIMEK_BOX = (0.04, 0.04, 0.40)              # 4 x 4 cm gain aperture, 40-cm path


def _p(**kw):
    p = engine_parameters()
    p.set_many(kw, source="problem")
    return p


def open_cavity(lz, area, wl, opt):
    return Cavity(l_c_m=lz, l_g_m=lz, K_axial=1e-5, R_OC=1e-6, k_c_loss=0.0,
                  L_op_fixed=0.9999, area_m2=area, wavelength_nm=wl,
                  sigma_se_m2=opt["sigma_se_m2"], absorbers=dict(opt["absorbers"]))


def deposited_J(dev, t_end_ns):
    """Beam energy deposited in the grid over [0, t_end]: the ledger's
    drive_in, evaluated WITHOUT integrating the plasma. The drive term is
    p_areal(t) x beam area x sum(deposition fractions) (engine.rhs), which
    depends on the beam spec and the fill only -- so it is exactly linear in
    the current density and matching a printed deposition needs no
    iteration. Gated against a full run in tests/test_excimer_anchors.py."""
    t = np.linspace(0.0, t_end_ns * 1e-9, 64001)
    pa = np.array([dev.spec.power_areal_W_m2(tt) for tt in t])
    return float(np.trapezoid(pa, t)) * dev._beam_area * float(dev._dep_frac.sum())


def build(mode, mixture, p_atm, kv, j, flat_ns, rise_ns, t_end_ns, box, wl,
          overrides=None, beam_axis=0):
    lx, ly, lz = box
    p = _p(mode=mode, drive="beam", pressure_atm=p_atm, beam_kV=kv,
           beam_j_A_cm2=j, beam_flat_ns=flat_ns, beam_rise_ns=rise_ns,
           t_end_ns=t_end_ns, nx=1, ny=1, nz=1, Lx_m=lx, Ly_m=ly, Lz_m=lz,
           beam_axis=beam_axis)
    if overrides:
        p.set_many(dict(overrides), source="problem")
    opt = mode_optics(mode)
    cav = open_cavity(lz, lx * ly, wl, opt)
    dev = Device(p, mixture=mixture, cavity=cav)
    if dev.beam_yields is None or dev.beam_yields.mode != "spencer_fano":
        raise RuntimeError(
            "beam yields are not on the Spencer-Fano route "
            f"({getattr(dev.beam_yields, 'cite', None)!r}): the LXCat data "
            "are missing and every anchor would be scored on the W-table "
            "fallback (2026-09-04 silent-fallback defect). Refusing to score.")
    return dev, cav, opt


def matched(mode, mixture, p_atm, kv, flat_ns, rise_ns, t_end_ns, box, wl,
            target_J, overrides=None, net_hook=None):
    """Build at a probe current, scale to the printed deposition (exact:
    deposition is linear in j), rebuild, and verify to 1e-9."""
    j0 = 100.0
    dev, _, _ = build(mode, mixture, p_atm, kv, j0, flat_ns, rise_ns,
                      t_end_ns, box, wl, overrides)
    j = j0 * target_J / deposited_J(dev, t_end_ns)
    dev, cav, opt = build(mode, mixture, p_atm, kv, j, flat_ns, rise_ns,
                          t_end_ns, box, wl, overrides)
    dep = deposited_J(dev, t_end_ns)
    if abs(dep / target_J - 1.0) > 1e-9:
        raise RuntimeError(f"deposition match failed: {dep} vs {target_J}")
    if net_hook is not None:
        net_hook(dev)
    return dev, cav, opt, j


def gain_alpha(res, cav, opt, drop=()):
    """Small-signal gain and absorption, 1/m (numerically %/cm), gain-volume
    averaged.  `drop` omits named absorbers to match a measurement's stated
    convention (Kannari's Table II excludes F2, p.4316-17)."""
    up = np.asarray(res["dens_avg"][opt["upper"]], dtype=float)
    g0 = cav.sigma_se_m2 * up
    alpha = np.zeros_like(g0)
    for sp, sig in cav.absorbers.items():
        if sp in res["dens_avg"] and sp not in drop:
            alpha += sig * np.asarray(res["dens_avg"][sp], dtype=float)
    return g0, alpha


def electra(overrides=None, net_hook=None, points=None, log=print):
    mix = {"Ar": 0.997, "F2": 0.003}
    rows = []
    for k, (p_atm, mw_cc, g_meas, lx, label) in enumerate(ELECTRA_GAIN):
        if points is not None and k not in points:
            continue
        box = (lx, 0.30, 1.0)
        target = mw_cc * 1e12 * (lx * 0.30 * 1.0) * 120e-9  # 20-ns ramps + 100-ns flat
        t0 = time.time()
        dev, cav, opt, j = matched("arf_beam", mix, p_atm, 500.0, 100.0, 20.0,
                                   160.0, box, 193.0, target, overrides, net_hook)
        res = dev.run(n_save=41)
        g0, al = gain_alpha(res, cav, opt)
        t = np.asarray(res["t_ns"], dtype=float)
        i_pk = int(np.argmax(g0))
        i_ft = int(np.argmin(np.abs(t - 115.0)))
        row = dict(anchor="electra_arf_gain", p_atm=p_atm, mw_cm3=mw_cc, label=label,
                   measured=g_meas, model_peak=float(g0[i_pk]), model_flat=float(g0[i_ft]),
                   ratio_peak=float(g0[i_pk] / g_meas), ratio_flat=float(g0[i_ft] / g_meas),
                   alpha_at_peak=float(al[i_pk]), t_peak_ns=float(t[i_pk]),
                   Tg_max_K=float(np.max(res["Tg_K"])), j_A_cm2=float(j),
                   deposited_J=float(res["final_channels_J"]["drive_in"]),
                   wall_s=time.time() - t0)
        rows.append(row)
        log(f"  Electra ArF {p_atm:.2f} atm {mw_cc:.3g} MW/cm3: meas {g_meas:.2f}  "
            f"model {row['model_peak']:.2f} (peak) / {row['model_flat']:.2f} (115 ns) "
            f"x{row['ratio_peak']:.2f} / x{row['ratio_flat']:.2f}  "
            f"Tg {row['Tg_max_K']:.0f} K  [{row['wall_s']:.0f} s]")
    return rows


def klimek(overrides=None, net_hook=None, quick=False, kr=KLIMEK_KR, log=print):
    rows = []
    lx, ly, lz = KLIMEK_BOX
    for jj in KLIMEK_J:
        col = 0
        for pa in KLIMEK_P:
            for f2 in KLIMEK_F2:
                c = col
                col += 1
                if quick and not (jj in (110.0, 20.0) and pa in (1.0, 3.0)
                                  and f2 == 0.002):
                    continue
                p_atm = pa * AMAGAT_TO_ATM_300K
                mix = {"Ar": 1.0 - kr - f2, "Kr": kr, "F2": f2}
                target = (KANNARI_MW_CM3_PER_ATM_A * p_atm * jj) * 1e12 \
                    * (lx * ly * lz) * 155e-9            # 5-ns ramps + 150-ns flat
                t0 = time.time()
                dev, cav, opt, j = matched("krf", mix, p_atm, 600.0, 150.0, 5.0,
                                           160.0, KLIMEK_BOX, 248.0, target,
                                           overrides, net_hook)
                res = dev.run(n_save=41)
                g0, al = gain_alpha(res, cav, opt, drop=("F2",))
                t = np.asarray(res["t_ns"], dtype=float)
                i = int(np.argmin(np.abs(t - 150.0)))
                gm, am = KLIMEK_GAIN[jj][c], KLIMEK_ABS[jj][c]
                row = dict(anchor="klimek_krf", J_A_cm2=jj, amagat=pa, f2=f2,
                           measured_gain=gm, model_gain=float(g0[i]),
                           ratio_gain=float(g0[i] / gm),
                           measured_abs=am, model_abs=float(al[i]),
                           ratio_abs=float(al[i] / am),
                           Tg_max_K=float(np.max(res["Tg_K"])),
                           j_engine_A_cm2=float(j), wall_s=time.time() - t0)
                rows.append(row)
                log(f"  Klimek KrF J {jj:5.0f} {pa:.0f} amg F2 {f2 * 100:.1f} %: gain meas "
                    f"{gm:5.1f} model {row['model_gain']:5.2f} x{row['ratio_gain']:.2f} | "
                    f"abs meas {am:.2f} model {row['model_abs']:.2f} x{row['ratio_abs']:.2f}  "
                    f"Tg {row['Tg_max_K']:.0f} K  [{row['wall_s']:.0f} s]")
    return rows


def _gm(xs):
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def summary(e_rows, k_rows):
    out = {}
    if e_rows:
        out["electra_gain_peak"] = _gm([r["ratio_peak"] for r in e_rows])
        out["electra_gain_flat"] = _gm([r["ratio_flat"] for r in e_rows])
        by_p = {r["p_atm"]: r for r in e_rows if r["label"].startswith("single")}
        if 0.82 in by_p and 1.36 in by_p:
            out["electra_trend_model"] = by_p[1.36]["model_peak"] / by_p[0.82]["model_peak"]
            out["electra_trend_measured"] = ELECTRA_TREND_MEASURED
    if k_rows:
        out["klimek_gain"] = _gm([r["ratio_gain"] for r in k_rows])
        out["klimek_abs"] = _gm([r["ratio_abs"] for r in k_rows])
        pexp = {}
        for jj in KLIMEK_J:
            r1 = [r for r in k_rows if r["J_A_cm2"] == jj and r["amagat"] == 1.0 and r["f2"] == 0.002]
            r3 = [r for r in k_rows if r["J_A_cm2"] == jj and r["amagat"] == 3.0 and r["f2"] == 0.002]
            if r1 and r3:
                pexp[f"{jj:.0f}"] = dict(
                    model=math.log(r3[0]["model_gain"] / r1[0]["model_gain"]) / math.log(3.0),
                    measured=math.log(r3[0]["measured_gain"] / r1[0]["measured_gain"]) / math.log(3.0))
        out["klimek_p_exponent"] = pexp
        # pump-independent: (g0/alpha)_model / (g0/alpha)_measured, per pressure
        ga = {}
        for pa in KLIMEK_P:
            rr = [r for r in k_rows if r["amagat"] == pa]
            if rr:
                ga[f"{pa:.0f}amg"] = _gm([(r["model_gain"] / r["model_abs"])
                                          / (r["measured_gain"] / r["measured_abs"])
                                          for r in rr])
        out["klimek_g0_over_alpha"] = ga
    return out


def run_all(overrides=None, net_hook=None, quick=False, log=print):
    t0 = time.time()
    log("Electra ArF small-signal gain (deposition matched) ...")
    e_rows = electra(overrides, net_hook, points=(0, 2) if quick else None, log=log)
    log("Klimek KrF grid (pump matched to Kannari's device calibration) ...")
    k_rows = klimek(overrides, net_hook, quick=quick, log=log)
    s = summary(e_rows, k_rows)
    return dict(summary=s, electra=e_rows, klimek=k_rows,
                overrides=dict(overrides or {}), quick=quick,
                wall_s=time.time() - t0)


def _fmt_summary(s):
    lines = []
    if "electra_gain_peak" in s:
        lines.append(f"  Electra ArF gain   x{s['electra_gain_peak']:.2f} (peak) / "
                     f"x{s['electra_gain_flat']:.2f} (flat top)")
    if "electra_trend_model" in s:
        lines.append(f"  Electra 0.82->1.36 atm  model x{s['electra_trend_model']:.2f}  "
                     f"measured x{s['electra_trend_measured']:.2f}")
    if "klimek_gain" in s:
        lines.append(f"  Klimek KrF gain    x{s['klimek_gain']:.2f}   absorption "
                     f"x{s['klimek_abs']:.2f}")
        for jj, v in s["klimek_p_exponent"].items():
            lines.append(f"  p-exponent @ {jj} A/cm2: model {v['model']:+.2f}  "
                         f"measured {v['measured']:+.2f}")
        for pa, v in s["klimek_g0_over_alpha"].items():
            lines.append(f"  g0/alpha  {pa}: x{v:.2f} of measured")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quick", action="store_true",
                    help="Electra 0.82/1.36 atm + Klimek 110/20 A/cm2 x 1/3 amagat corners")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="engine parameter override (repeatable)")
    ap.add_argument("--tag", default="baseline")
    ap.add_argument("--out", default=str(ROOT / "rebuild_out" / "excimer_anchors"))
    a = ap.parse_args()
    ov = {}
    for kv in a.set:
        k, v = kv.split("=", 1)
        try:
            v = json.loads(v)
        except json.JSONDecodeError:
            pass
        ov[k] = v
    rep = run_all(ov, quick=a.quick)
    print(_fmt_summary(rep["summary"]))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rep["generated"] = time.strftime("%Y-%m-%d %H:%M")
    (out / f"{a.tag}.json").write_text(json.dumps(rep, indent=1, default=float))
    print(f"{rep['wall_s']:.0f} s -> {out / (a.tag + '.json')}")


if __name__ == "__main__":
    main()
