"""The prototype ArF e-beam drive: an evidence tool, never part of LIGHT_ArF.

LIGHT's only energy input is a given power into the electron gas. An e-beam
does not heat the electron gas: it makes ions and excited atoms directly and
leaves its secondaries with a few eV each. This module is the survey's KrF
beam drive (LIGHT/docs/survey_2026-09-27/scripts/ebeam_screen.py) with Ar's
per-pair yields, so the ArF tail table (ar_tail.py) and the scoring
(score_arf.py) can run LIGHT_ArF under it beside LIGHT's own input, as the
KrF review ran both.

Per ion pair, REBUILD 0.1.40 MODULES/beam.spencer_fano_yields for Ar with
0.5 % F2 (IST-Lisbon + SIGLO, E0 = 516 eV): W = 30.85 eV, Ar+ 0.9981, Ar**
(4p and above) 0.5474, Ar* (4s) 0.2967, F 0.0044, F2+ 0.0019, and 4.398 eV
of sub-excitation energy (run 2026-09-30; the same call for Kr gives the
KrF drive's 27.07 eV, 0.9989, 0.5119, 0.2634, 4.506 eV exactly). As in the
KrF drive, Ar** cascades promptly to Ar* (REBUILD's Ar** 13.20 eV, the
excess to heat) and the F2 channels (0.6 %) are folded into the electron
share. The energy per pair is rebuilt from LIGHT_ArF's own E_form, so the
ledger closes exactly. The IST-Lisbon set gives W 17 % above the measured
26.4 eV; only the pair rate at a given pump depends on it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

LIGHT_ARF = next(p for p in Path(__file__).resolve().parents
                 if (p / "run_v0_2.py").is_file())
if str(LIGHT_ARF) not in sys.path:
    sys.path.insert(0, str(LIGHT_ARF))

from v0_2 import xsec  # noqa: E402
from v0_2.model import KB_EV, QE, Model  # noqa: E402
from v0_2.species import RG, RG_ION, RG_STAR, e_form  # noqa: E402

Y_ION, Y_PP, Y_S, EPS_SEC = 0.9981, 0.5474, 0.2967, 4.398
E_PP = 13.20
E_B = xsec.e_break()                     # 11.548 eV
CASCADE = Y_PP * (E_PP - e_form(RG_STAR))
W_EFF = Y_ION * e_form(RG_ION) + (Y_PP + Y_S) * e_form(RG_STAR) + CASCADE \
    + EPS_SEC                            # eV per pair
NEX_PER_ION = (Y_PP + Y_S) / Y_ION       # excited atoms per ion, 0.846


class BeamModel(Model):
    """LIGHT_ArF under the prototype beam drive: Params.power_W_cm3 is the
    beam's deposited power, pulse_s its duration."""

    def rhs(self, t, y, on):
        if not on:
            return super().rhs(t, y, False)
        dy = super().rhs(t, y, False)          # chemistry, no power input
        N0 = self.N0
        pairs = self.p.power_W_cm3 / QE / W_EFF
        src = np.zeros(self.ns)
        src[self.idx[RG]] = -pairs * (Y_ION + Y_PP + Y_S)
        src[self.idx[RG_ION]] = pairs * Y_ION
        src[self.idx[RG_STAR]] = pairs * (Y_PP + Y_S)
        dy[:self.ns] += src / N0
        dy[self.i_ue] += pairs * EPS_SEC / N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        tg_ev = y[self.i_tg] * KB_EV
        c_v = float((n * self.cv).sum()) * KB_EV
        dy[self.i_tg] += (pairs * CASCADE
                          - tg_ev * float((src * self.cv).sum())) / c_v
        dy[self.i_led + 0] += pairs * W_EFF / N0          # E_dep
        dy[self.i_led + 2] += pairs * CASCADE / N0        # cascade heat
        return dy
