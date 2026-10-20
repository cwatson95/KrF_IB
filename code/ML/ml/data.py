"""The experimental datasets: values, provenance, error model and role.

Roles encode the repo's predictive-only policy (KrF_and_LoKI/agents/CLAUDE.md,
PLAN_physical_validity.txt). Every rate comes from "a published measurement
of that process", and anchors "score the model after the fact; they never
steer a number":
  train     kinetic measurements of the channels being fitted
  held_out  kinetic measurements left out of the fit, to test it
  device    laser performance (gain, efficiency). These score and never
            train. tests/test_data.py fails if a device set is given the
            train role.

The train set is the largest group of kinetic datasets that agree with one
another within their errors (`run.py consistency`). Hakuta 1990 disagrees
with it and is held out. Changing a role is a modelling decision, and it
belongs in ML/CHANGELOG.md with the reason.

Densities use the figures' own axis, 3.22e16 cm^-3 per Torr (300 K).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field

from . import ROOT

TORR = 3.22e16
DIGITIZED = ROOT / "LIGHT" / "docs" / "krf_data" / "digitized"
SOURCES = ROOT / "LIGHT" / "docs" / "krf_data" / "sources"

RATE_KINDS = ("b_only_decay", "slow_decay")
DEVICE_KINDS = ("device_g0", "device_formation", "device_isat")
ROLES = ("train", "held_out", "device")


@dataclass(frozen=True)
class Dataset:
    name: str
    kind: str              # what the forward model must compute
    role: str
    source: str            # citation with page, table or figure
    points: tuple          # ((conditions dict, measured value), ...)
    sigma: float           # error per point
    sigma_is_fraction: bool = True
    file: str | None = None     # digitized CSV (relative to ROOT), if read
    note: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def values(self):
        return [v for _, v in self.points]


def _rows(name):
    with open(DIGITIZED / name) as fh:
        rows = list(csv.reader(l for l in fh if not l.startswith("#")))
    return rows[1:]


def _rel(name):
    return str((DIGITIZED / name).relative_to(ROOT))


def _eden():
    kr = tuple(({"kr_torr": float(p), "f2_torr": 0.0, "krf2_torr": 0.5},
                float(s) * 1e7) for p, s in _rows("eden1978_APL_fig3_Kr.csv"))
    f2 = tuple(({"kr_torr": 0.0, "f2_torr": float(p), "krf2_torr": 0.5},
                float(s) * 1e7) for p, s in _rows("eden1978_APL_fig5_F2.csv"))
    src = ("Eden, Waynant, Searles & Burnham, APL 32, 733 (1978): KrF2 "
           "photolysis at 193 nm makes cold KrF(B) alone. Single-exponential "
           "B-X decay, pump deconvolved, stated error <= 10 % per point")
    return [Dataset("eden_kr", "b_only_decay", "train", src + ", Fig. 3 "
                    "(against Kr)", kr, 0.10, file=_rel("eden1978_APL_fig3_Kr"
                                                         ".csv")),
            Dataset("eden_f2", "b_only_decay", "train", src + ", Fig. 5 "
                    "(against F2)", f2, 0.10, file=_rel("eden1978_APL_fig5_F2"
                                                         ".csv"))]


def _ksu():
    fig2 = tuple(({"kr_torr": float(p), "f2_torr": float(f)}, float(s) * 1e7)
                 for f, p, s in _rows("gadomski1992_fig2.csv"))
    fig1 = tuple(({"kr_torr": 100.0 if s_ == "Kr100_F2" else 50.0,
                   "f2_torr": float(p)}, float(v) * 1e7)
                 for s_, p, v in _rows("gadomski1992_fig1.csv") if s_ != "Xe")
    src = ("Gadomski, Xu, Setser & Nelson, CPL 189, 153 (1992): two-photon "
           "Kr + F2 makes hot B and C, slow decay from the latest part of the "
           "waveform")
    return [Dataset("ksu_fig2", "slow_decay", "train", src + ", Fig. 2 "
                    "(against Kr, 1 and 2 Torr F2)", fig2, 0.12,
                    file=_rel("gadomski1992_fig2.csv")),
            Dataset("ksu_fig1", "slow_decay", "train", src + ", Fig. 1 "
                    "(against F2, 50 and 100 Torr Kr)", fig1, 0.12,
                    file=_rel("gadomski1992_fig1.csv"))]


def _xu():
    # zero-added-gas intercepts of the He, Ne and Ar series (Figs. 5-7),
    # read by a search agent and checked by eye; the mean of the three
    xu = {50.0: (6.90, 6.76, 6.84), 100.0: (7.69, 7.75, 7.66),
          150.0: (8.52, 8.58, 8.37), 200.0: (9.23, 9.61, 9.13)}
    pts = tuple(({"kr_torr": p, "f2_torr": 2.0}, sum(v) / 3.0 * 1e7)
                for p, v in xu.items())
    return [Dataset("xu_figs5_7", "slow_decay", "held_out",
                    "Xu, Gadomski & Setser, JCP 99, 2591 (1993) Figs. 5-7 "
                    "p.2595-2596: pure Kr/F2 decays at 2 Torr F2 (the "
                    "intercepts of the He, Ne and Ar series)", pts, 0.12,
                    note="Same lab as ksu_fig1/2, a different series. It "
                    "differs from Gadomski Fig. 2b by up to 25 % at 50 Torr.")]


def _hakuta():
    fig3 = tuple(({"kr_torr": float(p), "f2_torr": 1.0}, float(v) * 1e7)
                 for p, v in _rows("hakuta1990_fig3_Kr.csv"))
    fig2 = tuple(({"kr_torr": 9.0, "f2_torr": float(p)}, float(v) * 1e7)
                 for p, v in _rows("hakuta1990_fig2_F2.csv"))
    src = ("Hakuta, Fujino, Nakayama & Takuma, JJAP 29, 1194 (1990): C-A "
           "slow decay after 248-nm photoassociation, 10-ns resolution")
    note = ("Inconsistent with the train set (leave-one-out chi2 about 15 per "
            "point). Its rates run about 1.8x the Kansas State rates at the "
            "same conditions.")
    return [Dataset("hakuta_fig3", "slow_decay", "held_out", src +
                    ", Fig. 3 (against Kr, 1 Torr F2)", fig3, 0.10,
                    file=_rel("hakuta1990_fig3_Kr.csv"), note=note),
            Dataset("hakuta_fig2", "slow_decay", "held_out", src +
                    ", Fig. 2 (against F2, 9 Torr Kr)", fig2, 0.10,
                    file=_rel("hakuta1990_fig2_F2.csv"), note=note)]


def _yu_setser():
    src = ("Yu & Setser, J. Phys. Chem. 94, 2934 (1990) p.2941: equilibrium "
           "C/B emission ratio at high pressure, B-A corrected")
    return [Dataset("ys_ir300", "eq_ratio", "train", src +
                    ", I_r(300 K) = 0.16 +- 0.02 (their eq. 18: I_r = "
                    "1.035 (tau_B/tau_C) K_eq)", (({"T": 300.0}, 0.16),), 0.02,
                    sigma_is_fraction=False),
            Dataset("ys_ratio_230_300", "eq_ratio_T", "held_out", src +
                    ", dE_BC = -180 +- 30 cm^-1 from the 230 and 300 K "
                    "ratios, so I_r(230)/I_r(300) = 1.30 +- 0.06",
                    (({"T1": 230.0, "T2": 300.0}, 1.30),), 0.06,
                    sigma_is_fraction=False,
                    note="Tests the B-C gap, which the fit does not vary "
                    "(LIGHT's ledger gap, 80 cm^-1). Their absolute-ratio "
                    "method gives -110 +- 20 cm^-1.")]


def _burnham_searles():
    return [Dataset("bs_tau_b", "tau_B", "held_out",
                    "Burnham & Searles, JCP 67, 5967 (1977), NRL MR 3753 "
                    "preprint p.56: tau_B = 9.0 +- 0.6 ns from ArF-laser "
                    "photolysis of KrF2 (late decay only)",
                    (({}, 9.0e-9),), 0.6e-9, sigma_is_fraction=False,
                    note="Eden measured 6.8 +- 0.2 ns, and theory gives 7.0.")]


def _tahoma():
    cond = {"pressure_torr": 665.0, "f2_torr": 2.8, "power_W_cm3": 3.3e5,
            "pulse_s": 450e-9}
    src = ("Spectra Technology final report to LANL, 1989 (Tahoma), OSTI "
           "10136715: 99.6 % Kr/F2, e-beam, 330 kW/cm3 at the probed centre")
    return [Dataset("tahoma_g0", "device_g0", "device", src +
                    ", Table 3-3: g0 = 6.39 %/cm (beam stop in)",
                    ((cond, 6.39),), 0.03, note="Scored through the "
                    "prototype beam drive in LIGHT/docs/survey_2026-09-27."),
            Dataset("tahoma_formation", "device_formation", "device", src +
                    ", p.3-16: formation efficiency g0 Isat/Pin = 0.27",
                    ((cond, 0.27),), 0.10, note="The model's value is its own "
                    "g0 x Isat / Pin (ml/formation.py, a 248 nm probe flux "
                    "with Isat fitted at fixed g0), Tahoma's definition."),
            Dataset("tahoma_isat", "device_isat", "device", src +
                    ", p.3-16: saturation flux 1.4 +- 0.1 MW/cm2 (modified "
                    "Rigrod analysis of the extraction data)",
                    ((cond, 1.4),), 0.07, note="The model's Isat from its "
                    "saturated gain at 0.5-4 MW/cm2. It measures the KrF(B) "
                    "lifetime, not the formation.")]


def all_datasets() -> list:
    return (_eden() + _ksu() + _xu() + _hakuta() + _yu_setser()
            + _burnham_searles() + _tahoma())


def by_role(role: str) -> list:
    return [d for d in all_datasets() if d.role == role]


def by_name(*names) -> list:
    table = {d.name: d for d in all_datasets()}
    return [table[n] for n in names]
