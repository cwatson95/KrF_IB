"""Blind test of LIGHT's basic network against Tahoma (evidence script, not
LIGHT code).

Tahoma is the one full-text, buffer-free Kr/F2 gain dataset found
(Spectra Technology final report to LANL, 19 Jul 1989, OSTI 10136715,
sources/Tahoma_SpectraTech_1989_OSTI10136715.pdf):
  99.6 % Kr/F2, F2 2.8 Torr, 665 Torr total, e-beam, 450 ns (p. 2-2).
  Beam stop in: 330 kW/cm3 deposited at the probed centre (p. 2-3),
  g0 = 6.39 %/cm (Table 3-3), formation efficiency g0*Isat/Pin = 90/330 = 0.27
  (p. 3-16, Isat = 1.4 MW/cm2).

The drive is the PROTOTYPE e-beam source in
../survey_2026-09-27/scripts/ebeam_screen.py (REBUILD's Spencer-Fano yields
for Kr with 0.5 % F2). g0 = sigma_B N_B with sigma_B = 2.6e-16 cm2
(Tellinghuisen 1976 via Kannari 1985 Table I, as in REBUILD mode_optics).
No value is adjusted to Tahoma.

The model's formation efficiency is (KrF formation rate x 5.0 eV) / Pin over
200-400 ns. Tahoma's is g0*Isat/Pin. The two agree only for a simple
two-level medium, so read that comparison as approximate. g0 compares
directly.

Since 2026-09-27 (later) the basic network applies K_eq(Tg) by detailed
balance over the ledger's B-C gap (80 cm-1). The 200 and 300 cm-1 rows
replace that gap with the one Gadomski et al. 1992 p.156 printed (about
300 cm-1, "the collective data favor ... >200").

    /usr/bin/python3 LIGHT/docs/krf_data/tahoma_check.py      # about 1 min
"""
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_LIGHT / "docs" / "survey_2026-09-27" / "scripts"))
sys.path.insert(0, str(_LIGHT))
import numpy as np
from ebeam_screen import BeamModel
from v0.model import Params, QE
from v0.network import KRF_K_BC, KRF_K_EQ

SIGMA_B = 2.6e-16          # cm2
HV_EV = 5.00               # 248 nm
CM1_K = 1.4387769          # K per cm-1
MEAS_G0, MEAS_ETA_F = 6.39, 0.27
LEDGER_GAP_CM1 = 79.849    # E(KrF B) - E(KrF C) in species.py (9.9 meV)


class GapBeamModel(BeamModel):
    """K_eq(Tg) = 5.0 exp(gap (1/Tg - 1/300) hc/k) on the C -> B row."""
    gap_cm1 = 0.0

    def _k(self, c, te, tg, n, ne):
        if c.row.rid == "CB" and self.gap_cm1:
            return KRF_K_BC / keq(self.gap_cm1, tg)
        return super()._k(c, te, tg, n, ne)


def keq(gap_cm1, tg):
    return KRF_K_EQ * math.exp(gap_cm1 * CM1_K * (1.0 / max(tg, 1.0)
                                                  - 1.0 / 300.0))


def check(groups, gap_cm1=0.0, f2_torr=2.8):
    p = Params(network="basic", groups=groups, pressure_atm=665 / 760,
               x_f2=f2_torr / 665, power_W_cm3=3.3e5, pulse_s=450e-9,
               t_end_s=600e-9, ne_mode="free", ne_cm3=1e10, n_save=1201)
    m = GapBeamModel(p)
    m.gap_cm1 = gap_cm1
    r = m.run()
    t = r.t
    man = ("KrF(B)", "KrF(C)")
    formed = np.zeros_like(t)      # KrF made from outside the B/C manifold
    for c in m.compiled:
        net = (sum(c.row.products.count(s) for s in man)
               - sum(c.row.reactants.count(s) for s in man))
        if net > 0:
            formed += net * r.events[c.row.rid]
    win = (t >= 200e-9) & (t <= 400e-9)
    rate = np.polyfit(t[win], formed[win], 1)[0]           # cm-3 s-1
    eta_f = rate * HV_EV * QE / p.power_W_cm3
    g0 = SIGMA_B * float(np.mean(r.n["KrF(B)"][win])) * 100.0   # %/cm
    tg = float(np.mean(r.tg[win]))
    return dict(eta_f=eta_f, g0=g0, ne=float(np.mean(r.ne[win])),
                te=float(np.mean(r.te[win])), tg=tg,
                keq=keq(gap_cm1 or LEDGER_GAP_CM1, tg))


CASES = (
    ("basic + ions + f2e (ledger gap 80)", ("ions", "f2e"), 0.0),
    ("  same, K_eq(Tg) with a 200 cm-1 gap", ("ions", "f2e"), 200.0),
    ("  same, K_eq(Tg) with a 300 cm-1 gap", ("ions", "f2e"), 300.0),
    ("basic + all optional groups", ("all",), 0.0),
    ("basic alone (no recombination)", (), 0.0),
)

if __name__ == "__main__":
    print(f"Tahoma, beam stop in (99.6 % Kr, 665 Torr, 330 kW/cm3): "
          f"g0 {MEAS_G0} %/cm, formation efficiency {MEAS_ETA_F}")
    print(f"{'case':38s} {'g0 %/cm':>8s} {'x meas':>7s} {'form.':>6s} "
          f"{'x meas':>7s} {'n_e':>8s} {'Te':>5s} {'Tg':>5s} {'K_eq':>5s}")
    for label, groups, gap in CASES:
        try:
            d = check(groups, gap)
        except Exception as exc:        # report, do not hide, a solver failure
            print(f"{label:38s} failed: {str(exc)[:60]}")
            continue
        print(f"{label:38s} {d['g0']:8.2f} {d['g0'] / MEAS_G0:7.2f} "
              f"{d['eta_f']:6.3f} {d['eta_f'] / MEAS_ETA_F:7.2f} "
              f"{d['ne']:8.1e} {d['te']:5.2f} {d['tg']:5.0f} {d['keq']:5.2f}")
