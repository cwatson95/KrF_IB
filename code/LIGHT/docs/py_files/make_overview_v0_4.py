"""Numbers and tables for accuracy_overview_v0_4.tex, v0_4's accuracy.

v0_4 is v0_3 with the 5p rung of krypton, the 248-nm optics, KrF(B,C) + F
as a band, ion-ion recombination's Langevin limb and bands, the ions'
drift to the walls and the EM wave with Turnbull et al.'s inverse
bremsstrahlung (v0_4/README.md). The tests without electrons (the
photolysis decays, the C/B ratios) score as v0_2 does, and this script
checks that they do. What v0_4 changes shows at the reference point and
at Tahoma, where the absorption is scored for the first time. No number
in the overview is typed by hand. This script writes two files, each one
\\wfdef{key}{value} per number or a finished table (a macro namespace of
its own, so the overview can hold v0_2's, v0_3's and the map's numbers
beside v0_4's):

  overview_v0_4_fast.tex   minutes: ML's scorecard under v0_4, checked
                           equal to v0_2's, the photolysis decay constants,
                           Tahoma under v0_2's electron closure for v0_4,
                           its variants and v0_3, the saturation intensity
                           and formation efficiency, the wave and the walls
                           at the reference point, and the new terms' own
                           checks (the 5s <-> 5p pseudo cross section, the
                           wave's field term, the electron-ion share)
  overview_v0_4_slow.tex   the reference point and Tahoma under the
                           kinetic EEDF, each a process of its own
  overview_v0_4_unc.tex    the propagated uncertainty: every band at its
                           two ends and every constant with a quoted error
                           at its two ends, the 248-nm share at the
                           reference point and Tahoma's g0, alpha0 and
                           g0 - alpha0 under v0_2's closure, combined in
                           quadrature and as an envelope (~15 min, 2
                           processes for Tahoma)

The reference point under v0_2's closure, its bands and its optics come
from the system map's numbers (system_map_v0_4_numbers.tex), which the
overview inputs. The accuracy numbers come from ML/, which runs v0_4 when
LIGHT_VERSION names it. ML fixes its version when first imported, so
main() sets LIGHT_VERSION=v0_4 before ML is imported, and every worker
is a fresh process that sets its own. v0_4/tests/test_docs.py compares
the files through --check.

    python3 LIGHT/docs/make_overview_v0_4.py               # all, ~40 min
    python3 LIGHT/docs/make_overview_v0_4.py --fast        # fast only, ~15 min
    python3 LIGHT/docs/make_overview_v0_4.py --slow        # slow only, ~12 min
    python3 LIGHT/docs/make_overview_v0_4.py --unc         # the uncertainty only, ~15 min
    python3 LIGHT/docs/make_overview_v0_4.py --check fast  # compare only
    python3 LIGHT/docs/make_energy_ladder.py v0_4
    cd LIGHT/docs && latexmk -pdf accuracy_overview_v0_4.tex
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/, the directory holding run_v0_4.py (no parent count)
LIGHT = next(p for p in HERE.parents if (p / "run_v0_4.py").is_file())
ML_DIR = LIGHT.parent / "ML"
for _p in (str(LIGHT), str(HERE), str(ML_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                            # noqa: E402

import make_overview as MO                                    # noqa: E402
from make_system_map_basic import _k                          # noqa: E402
from v0_4 import __version__                                  # noqa: E402
from v0_4.network import BANDS, KR_5P_TAU_S, LADDER_KANNARI   # noqa: E402

VERSION = "v0_4"
OUT_FAST = HERE / "overview_v0_4_fast.tex"
OUT_SLOW = HERE / "overview_v0_4_slow.tex"
OUT_UNC = HERE / "overview_v0_4_unc.tex"
V02_FAST = HERE / "overview_v0_2_fast.tex"
TAHOMA_SCRIPT = HERE / "krf_data" / "tahoma_v0_4.py"
WORKERS = 2
TRIMER_WEAK = -0.5
TWO_T = {"eedf": "two_t"}
BANDS_SHOWN = ("248", "272", "275", "420", "147", "IR")
BAND_ROWS = {"248": "G4", "272": "G4A", "275": "G4C", "420": "C3",
             "147": "KR2_RAD", "IR": "PRAD"}
# the Tahoma variants under v0_2's closure: label -> Params fields
TAHOMA_VARIANTS = {
    "default": {},
    "kannari": {"ladder": "kannari"},
    "nocycle": {"rates": (("k_5s5p", 0.0),)},
    "kpslo": {"rates": (("kps_kr", BANDS["kps_kr"].lo),)},
    "krff": {"rates": (("krf_f", BANDS["krf_f"].hi),)},
    "kr2f": {"rates": (("kr2_f", BANDS["kr2_f"].hi),)},
    "tri": {"trimer_tg_pow": TRIMER_WEAK},
    "iiyield": {"rates": (("ii_yield", BANDS["ii_yield"].lo),)},
    "iifc": {"rates": (("ii_f_c", BANDS["ii_f_c"].hi),)},
}
SIGMA_K2F_COMPILATION = 6.0e-18     # Kannari Table I [16], the estimate
WAVE = {"em_E_V_cm": 1500.0, "em_f_Hz": 95e9}
WAVE_THZ = {"em_E_V_cm": 3000.0, "em_f_Hz": 1e12}
GAP_CM = 1.0
THIN = {"gap_cm": 0.1, "pressure_atm": 0.1}
# the propagated uncertainty's inputs beyond the bands: the constants the
# model uses at a measured value whose source quotes an error and that no
# band covers (row id -> relative error, the source), and the two cross
# sections of the optics with a quoted error (analytic: g0 is sigma_se
# times n_B, the Kr2F* share of alpha is its sigma times n). The fitted
# B/C set (tau_B, tau_C, k_BC, the B and C quenches, the trimer) is not
# propagated: it is pinned by the training decays and scored above
UNC_ROWS = {
    "R11p": (0.10, r"Shiu and Biondi 1977, $\pm10$\,\%"),
    "K15": (0.09, r"Gr\"ossl et al.\ 1981, $(2.2\pm0.2)\times10^{-31}$"),
    "C3": (0.09, r"Xu et al.\ 1993, $220\pm20$\unit{ns}"),
}
UNC_SIGMA = {
    "se": (0.12, r"Rokni et al.\ 1978, $(2.5\pm0.3)\times10^{-16}$"),
    "k2f": (0.20, r"Hakuta et al.\ 1987, $(1.0\pm0.2)\times10^{-18}$"),
}
UNC_LABELS = {
    "R11p": r"\Krtp{} + \el{} \ra \Krpp{} + Kr",
    "K15": r"\Krp{} + 2Kr \ra \Krtp{} + Kr",
    "C3": r"\KtF{} radiative lifetime",
    "se": r"$\sigma_\mathrm{se}$, 248-nm stimulated emission",
    "k2f": r"$\sigma$(\KtF{}) at 248\unit{nm}",
}


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def _render(n: dict, what: str, tables=()) -> str:
    lines = [
        f"% LIGHT v0_4 accuracy numbers ({what}). GENERATED by "
        "make_overview_v0_4.py",
        f"% from LIGHT v0_4 ({__version__}) and ML/. Do not edit. Rerun the "
        "script instead.",
        r"\providecommand\wfdef[2]{\expandafter\def\csname wf@#1\endcsname{#2}}",
    ]
    lines += [f"\\wfdef{{{k}}}{{{v}}}" for k, v in n.items()]
    lines += list(tables)
    return "\n".join(lines) + "\n"


def _e(v: float, digits: int = 1) -> str:
    """A small number as 1.2e-5 in LaTeX: $1.2{\\times}10^{-5}$."""
    if v == 0.0:
        return "0"
    mant, ex = f"{v:.{digits}e}".split("e")
    return f"${mant}{{\\times}}10^{{{int(ex)}}}$"


def _ml():
    import ml
    if ml.LIGHT_VERSION != VERSION:
        raise RuntimeError(f"ML runs LIGHT {ml.LIGHT_VERSION} in this "
                           "process. Run make_overview_v0_4.py as a program, "
                           "which sets LIGHT_VERSION=v0_4 before ML loads")
    return ml


def _v02_numbers() -> dict:
    return dict(re.findall(r"\\widef\{([^}]*)\}\{(.*)\}$",
                           V02_FAST.read_text(), flags=re.M))


def _sp(s: str) -> str:
    return {"F-": r"\Fm{}", "Kr2+": r"\Krtp{}", "Kr**": r"\Krpp{}",
            "F2": r"\Ft{}", "Kr2F*": r"\KtF{}", "Kr*": r"\Krs{}"}[s]


# ---------------------------------------------------------------------------
# the fast numbers
# ---------------------------------------------------------------------------
def _krf_decay(kr_torr, f2_torr):
    """KrF(B) seeded in Kr/F2 at 300 K with no power, under v0_4's own
    closure (a seed of 1 electron per cm3 leaves the kinetic EEDF to the
    Maxwellian), and the slow decay of KrF(B) + KrF(C) fitted."""
    from v0_4.model import Model, Params
    tot = kr_torr + f2_torr
    p = Params(pressure_atm=tot / 760.0, x_f2=f2_torr / tot,
               power_W_cm3=0.0, pulse_s=0.0, t_end_s=400e-9,
               ne_mode="free", ne_cm3=1.0, n_save=2001,
               init_extra_cm3=(("KrF(B)", 1e12),))
    r = Model(p).run()
    n = r.n["KrF(B)"] + r.n.get("KrF(C)", 0.0)
    rel = n / n[0]
    sel = (r.t > 20e-9) & (rel < 1e-3) & (rel > 1e-8)
    return -np.polyfit(r.t[sel], np.log(n[sel]), 1)[0]


def _tahoma_two_t(**kw) -> dict:
    """Tahoma under the beam drive on v0_2's closure, for v0_4 in this
    process: g0, alpha by species, the window's Te and n_e."""
    _ml()
    from ml import forward
    from v0_4 import optics
    p = forward.tahoma_params(closure=TWO_T).replace(**kw)
    r = forward._drive().BeamModel(p).run()
    win = (r.t >= forward.WINDOW[0]) & (r.t <= forward.WINDOW[1])
    g0 = 100.0 * optics.SIGMA_SE_CM2 * float(np.mean(r.n["KrF(B)"][win]))
    by = {s: 100.0 * sig * float(np.mean(r.n[s][win]))
          for s, (sig, _) in optics.ABSORBERS.items() if s in r.n}
    return {"g0": g0, "alpha": sum(by.values()), "by": by,
            "te": float(np.mean(r.te[win])), "ne": float(np.mean(r.ne[win])),
            "tg": float(r.tg[-1])}


def _tahoma_other_version(version: str) -> dict:
    """Tahoma under v0_2's closure for another LIGHT version, through the
    scoring script in a process of its own (ML binds its version once)."""
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "t.json"
        env = dict(os.environ, LIGHT_VERSION=version,
                   PYTHONDONTWRITEBYTECODE="1")
        r = subprocess.run([sys.executable, "-W", "ignore",
                            str(TAHOMA_SCRIPT), "--json", str(out)],
                           env=env, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stdout + r.stderr)
        return json.loads(out.read_text())["Kr/F2"]


def _isat_two_t() -> dict:
    """ML's Tahoma observables (g0, Isat, g0 Isat / Pin) on v0_4 under
    v0_2's closure: the probe runs at each flux of ML's fit, with
    tahoma_params told the closure."""
    _ml()
    from ml import formation, forward, params
    c = forward.with_rederived_k_e(params.light_now())
    orig = forward.tahoma_params

    def two_t(cond=None, kr2=True, closure=None, afterglow_s=150e-9):
        return orig(cond, kr2, closure or TWO_T, afterglow_s)
    forward.tahoma_params = two_t
    try:
        o = formation.observables(c)
    finally:
        forward.tahoma_params = orig
    return {"g0": o["g0"], "isat": o["isat"], "eta": o["eta"],
            "formation": o["formation"]}


def _reference_two_t(**kw):
    from v0_4.model import Params, run
    from v0_4.report import accounting, row_dict
    r = run(Params(**TWO_T, **kw))
    return r, accounting(r), row_dict(r)


def _wave_numbers() -> dict:
    """The wave at the reference point under v0_2's closure: its absorbed
    energy beside the given one, the band's two ends, the electron-ion
    share of the absorption at mid-pulse, and the THz case."""
    from v0_4 import ib
    from v0_4.model import Model, Params
    n = {}
    r, a, d = _reference_two_t(**WAVE)
    n["wv.E"] = f"{WAVE['em_E_V_cm']:.0f}"
    n["wv.f"] = f"{WAVE['em_f_Hz'] / 1e9:.0f}"
    n["wv.em"] = f"{a['e_em'] * 1.602176634e-19:.3f}"
    n["wv.given"] = f"{a['e_given'] * 1.602176634e-19:.1f}"
    n["wv.share"] = f"{100.0 * a['e_em'] / a['e_dep']:.1f}"
    n["wv.248"] = f"{d['emission_248']:.2f}"
    n["wv.closure"] = _e(abs(a["closure"]))
    ends = {}
    for tag, v in (("lo", BANDS["ib_scale"].lo), ("hi", BANDS["ib_scale"].hi)):
        _, ab, _ = _reference_two_t(rates=(("ib_scale", v),), **WAVE)
        ends[tag] = ab["e_em"]
    n["wv.band"] = f"{100.0 * (ends['hi'] - ends['lo']) / a['e_em']:.2f}"
    # the electron-ion share of the absorption at mid-pulse, Maxwellian
    # kernel at the run's own state, and at the notes' densities
    m = Model(Params(**TWO_T, **WAVE))
    k = int(r.t.searchsorted(0.5 * r.p.pulse_s))
    nvec = np.array([r.n[s][k] for s in m.species])
    ne, te = float(r.ne[k]), float(r.te[k])
    n_gas = nvec[m.idx["Kr"]] + nvec[m.idx["F2"]]
    ni = sum(nvec[i] for i in m.cations)
    w = 2.0 * math.pi * WAVE["em_f_Hz"]
    x_gas = {"Kr": nvec[m.idx["Kr"]] / n_gas, "F2": nvec[m.idx["F2"]] / n_gas}
    args = (te, w, WAVE["em_E_V_cm"] * 100.0, n_gas * 1e6, x_gas, ne * 1e6)
    parts = ib.em_power_per_electron(*args, ni * 1e6, parts=True)
    atoms = ib.em_power_per_electron(*args, 0.0)
    n["wv.mid.te"] = f"{te:.2f}"
    n["wv.mid.ne"] = f"${_k(ne)}$"
    # the electron-ion share is the MARGINAL one (the power with the ions
    # less the power without them), what the band acts on; the kernel's
    # printed "ei" is the nu_ei share of the numerator, beside it
    n["wv.mid.ei"] = f"{100.0 * (parts['total'] - atoms) / parts['total']:.2f}"
    n["wv.mid.ei.share"] = f"{100.0 * parts['ei'] / parts['total']:.2f}"
    n["wv.mid.fL"] = f"{parts['f_L']:.4f}"
    n["wv.mid.nuw"] = f"{parts['nu_en_mean'] / w:.1f}"
    for ne_, tag in ((1e15, "15"), (1e16, "16")):
        for f_, ftag in ((95e9, "ghz"), (1e12, "thz")):
            a_ = (1.5, 2.0 * math.pi * f_, 1.5e5, 2.446e25,
                  {"Kr": 0.995, "F2": 0.005}, ne_ * 1e6)
            pt = ib.em_power_per_electron(*a_, ne_ * 1e6, parts=True)
            at = ib.em_power_per_electron(*a_, 0.0)
            n[f"wv.ei.{tag}.{ftag}"] = f"{100.0 * (pt['total'] - at) / pt['total']:.2f}"
    rt, at, dt = _reference_two_t(power_W_cm3=0.0, **WAVE_THZ)
    n["wv.thz.E"] = f"{WAVE_THZ['em_E_V_cm'] / 1e3:.0f}"
    n["wv.thz.em"] = _e(at["e_em"] * 1.602176634e-19)
    n["wv.err"] = f"{100 * ib.TURNBULL_MODEL_ERROR:.0f}"
    return n


def _wall_numbers() -> dict:
    from v0_4.model import Model, Params, WALL_E_AMBIPOLAR_TE, \
        WALL_E_ELECTRON_TE
    from v0_4.network import ION_K0
    n = {}
    r, a, d = _reference_two_t(gap_cm=GAP_CM)
    n["wl.gap"] = f"{GAP_CM:g}"
    n["wl.share"] = _e(100.0 * a["walls"] / a["e_dep"])          # percent
    n["wl.ionheat"] = _e(100.0 * float(r.ledger["E_ion_heat"][-1]) / a["e_dep"])
    n["wl.248"] = f"{d['emission_248']:.2f}"
    n["wl.closure"] = _e(abs(a["closure"]))
    rt, at, dt = _reference_two_t(**THIN)
    n["wl.thin.gap"] = f"{THIN['gap_cm']:g}"
    n["wl.thin.p"] = f"{THIN['pressure_atm']:g}"
    n["wl.thin.share"] = f"{100.0 * at['walls'] / at['e_dep']:.3f}"
    n["wl.K0.Krp"] = f"{ION_K0['Kr+']:.2f}"
    n["wl.K0.Kr2p"] = f"{ION_K0['Kr2+']:.2f}"
    n["wl.K0.Fm"] = f"{ION_K0['F-']:.2f}"
    n["wl.eTe"] = f"{WALL_E_ELECTRON_TE:g}"
    n["wl.iTe"] = f"{WALL_E_AMBIPOLAR_TE:g}"
    # the ambipolar rate at the reference point's mid-pulse state
    m = r.model
    k = int(r.t.searchsorted(0.5 * r.p.pulse_s))
    nvec = np.array([r.n[s][k] for s in m.species])
    loss = m.wall_loss(float(r.ue[k]), float(r.tg[k]), nvec)
    nu = max(f / nvec[i] for i, q, f in loss if q > 0 and nvec[i] > 0)
    n["wl.nu"] = f"{nu:.0f}"
    # the 1-D model's walls line at the same reference point, from its
    # own generated numbers (1Dv0_2 on its finest uniform mesh; that line
    # also holds its transport step's mean-energy cut in the wall cells,
    # so it bounds this closure from above only). 1Dv0_3 prints no walls
    # line of its own
    d1 = dict(re.findall(r"\\didef\{([^}]*)\}\{(.*)\}$",
                         (HERE / "overview_1Dv0_2_slow.tex").read_text(),
                         flags=re.M))
    n["wl.1d"] = d1["ref.256.walls"]
    n["wl.1d.cells"] = d1["mesh.hi"]
    return n


def _solver_checks() -> dict:
    """The new terms' own checks, as numbers: the 5s <-> 5p pseudo cross
    section against Kannari's law and detailed balance, the wave's field
    term against the Maxwellian kernel, the pops' energy balance in a
    short kinetic pulse, the ion-ion limbs."""
    from v0_4 import boltzmann as KB
    from v0_4 import ib
    from v0_4.model import Model, Params
    n = {}
    m = Model(Params(eedf="kinetic"))
    g = m.kin._grid(40.0)
    g = KB.Grid(2400, 40.0, m.kin.procs, m.kin.superelastic, m.kin.pops)
    q = g.pops["Kr*"]
    worst_k = worst_db = 0.0
    for te in (1.0, 1.5, 2.5):
        F = g.maxwellian(te)
        kf, kr = g.pop_rates("Kr*", F)
        kann = 9.0e-7 * te ** 0.5 * math.exp(-1.51 / te) * 1e-6
        worst_k = max(worst_k, abs(kf / kann - 1.0))
        worst_db = max(worst_db, abs(kr / kf / (q.g_ratio * math.exp(
            (1.51 - 0.125) / te)) - 1.0))
    n["pop.k"] = _e(worst_k)
    n["pop.db"] = _e(worst_db)
    n["pop.cells"] = "2400"
    n["pop.emax"] = "40"
    n["pop.sigma0"] = _k(float(q.p.sigma.max()))
    # the wave's field term: a Maxwellian handed to the solver against the
    # kernel's atom and ion parts
    mw = Model(Params(eedf="kinetic", **WAVE))
    nvec = np.zeros(mw.ns)
    nvec[mw.idx["Kr"]] = 0.995 * mw.N0
    nvec[mw.idx["F2"]] = 0.005 * mw.N0
    nvec[mw.idx["Kr+"]] = 1e15
    y = mw.initial_state()
    y[:mw.ns] = nvec / mw.N0
    te = 1.5
    y[mw.i_ue] = 1.5 * 1e15 * te / mw.N0
    s, _, _ = mw._kinetic_state(y, 0.0, True)
    gw = mw.kin._grid(12.0)
    mw.kin.grid = gw
    F = gw.maxwellian(te)
    sol = KB.Solution(F, 0.0, 0.0, gw.mean(F), False, 0, gw)
    p_all = mw.kin.em_power(s, sol)
    p_at = mw.kin.em_power(s, sol, ions=False)
    w = 2.0 * math.pi * WAVE["em_f_Hz"]
    args = (te, w, WAVE["em_E_V_cm"] * 100.0, s.n_col,
            {"Kr": 0.995, "F2": 0.005}, 1e21)
    parts = ib.em_power_per_electron(*args, 1e21, parts=True, langdon=False)
    atoms = ib.em_power_per_electron(*args, 0.0, langdon=False)
    marginal = parts["total"] - atoms
    n["hf.atoms"] = _e(abs(p_at / parts["en"] - 1.0))
    # the solver's ion part is the marginal one; the kernel's printed
    # share (nu_ei's share of the numerator) is another quantity where
    # nu_en is not << omega
    n["hf.ions"] = _e(abs((p_all - p_at) / marginal - 1.0))
    n["hf.share"] = f"{parts['ei'] / marginal:.2f}"
    n["hf.nuw"] = f"{parts['nu_en_mean'] / w:.1f}"
    n["hf.te"] = f"{te:g}"
    n["hf.ne"] = "$10^{15}$"
    # a short kinetic pulse: the field's share of the input and the gap
    # between the two energy balances with the pair inside the equation
    r = Model(Params(eedf="kinetic", pulse_s=60e-9, t_end_s=70e-9,
                     n_save=141)).run()
    on = [e for e in r.kinetic if e["on"] and not e.get("skipped")
          and e["t"] >= 30e-9]
    ratio = [e["power_ratio"] for e in on if np.isfinite(e["power_ratio"])]
    n["pop.short.t"] = "60"
    n["pop.short.ratio.lo"] = f"{min(ratio):.3f}"
    n["pop.short.ratio.hi"] = f"{max(ratio):.3f}"
    n["pop.short.gap"] = _e(max(abs(e["balance"]) for e in on))
    # the ion-ion limbs at 1 amagat
    from v0_4.network import ION_ION_K3, langevin_k
    k3, kl = ION_ION_K3["Kr+"], langevin_k("Kr+", "F-")
    n["ii.k3"] = _k(k3)
    n["ii.kL"] = _k(kl)
    n["ii.k1"] = _k(k3 * kl / (k3 + kl))
    n["ii.k1.v03"] = _k(2.4 / 2.7 * 3.417e-6)
    n["ii.kL.v03"] = _k(2.4 * 3.417e-6)
    k3b, klb = ION_ION_K3["Kr2+"], langevin_k("Kr2+", "F-")
    n["ii.k1b"] = _k(k3b * klb / (k3b + klb))
    n["ii.kLb"] = _k(klb)
    return n


def _ib_checks() -> str:
    """How many checks test_ib.py collects."""
    r = subprocess.run([sys.executable, "-m", "pytest", "-q",
                        "--collect-only", "-p", "no:cacheprovider",
                        str(LIGHT / "v0_4" / "tests" / "test_ib.py")],
                       capture_output=True, text=True, cwd=str(LIGHT.parent))
    m = re.search(r"(\d+) tests? collected", r.stdout)
    if m is None:
        raise RuntimeError(r.stdout + r.stderr)
    return m.group(1)


def _constants() -> dict:
    from v0_4 import optics
    from v0_4.network import ION_K0, KANNARI_K76, KANNARI_KSTEP_P, \
        KANNARI_KSTEP_S
    n = {"tau5p": f"{KR_5P_TAU_S * 1e9:.1f}",
         "k.KPS": _k(BANDS["kps_kr"].default),
         "k.KPS.lo": _k(BANDS["kps_kr"].lo),
         "k.K27": _k(BANDS["k5p_f2"].default),
         "k.76.A": _k(KANNARI_K76[0]), "k.76.B": f"{KANNARI_K76[1]:.2f}",
         "k.76.C": f"{KANNARI_K76[2]:.2f}",
         "k.stepS.hi": _k(BANDS["kstep_s"].hi),
         "k.stepS.kannari": _k(LADDER_KANNARI["kstep_s"]),
         "k.stepP.hi": _k(BANDS["kstep_p"].hi),
         "k.krff.hi": _k(BANDS["krf_f"].hi),
         "sig.B": _k(optics.SIGMA_SE_CM2),
         "sig.K2F": _k(optics.ABSORBERS["Kr2F*"][0]),
         "sig.K2F.comp": _k(SIGMA_K2F_COMPILATION),
         "sig.Fm": _k(optics.ABSORBERS["F-"][0]),
         "sig.Kr2p": _k(optics.ABSORBERS["Kr2+"][0]),
         "sig.Kpp": _k(optics.ABSORBERS["Kr**"][0]),
         "sig.F2": _k(optics.ABSORBERS["F2"][0]),
         "sig.Krs": _k(optics.ABSORBERS["Kr*"][0]),
         "n.abs": str(len(optics.ABSORBERS)),
         "n.bands": str(len(BANDS)),
         "n.zero": str(sum(1 for b in BANDS.values() if b.default == 0.0)),
         "ii.yield.lo": f"{BANDS['ii_yield'].lo:g}",
         "ii.fc.hi": f"{BANDS['ii_f_c'].hi:g}",
         "ib.lo": f"{BANDS['ib_scale'].lo:g}",
         "ib.hi": f"{BANDS['ib_scale'].hi:g}"}
    for s, k0 in ION_K0.items():
        n[f"K0.{s.replace('+', 'p').replace('-', 'm')}"] = f"{k0:.2f}"
    return n


def fast_numbers() -> tuple:
    """(numbers, tables)."""
    _ml()
    from ml import data, forward, params, score
    n = {}
    v02 = _v02_numbers()
    # ML's scorecard under v0_4, which must be v0_2's: these tests have
    # no electrons and no Kr**, and the new rows are 0 or absent there
    c = forward.with_rederived_k_e(params.light_now())
    same = 0
    for ds in data.all_datasets():
        if ds.kind in data.DEVICE_KINDS:
            continue
        v = score.kinetic_score(ds, c)
        s = MO._pct_signed(v) if len(ds.values) == 1 else f"{v:.1f}"
        n[f"sc.{ds.name}"] = s
        if v02.get(f"sc.{ds.name}") != s:
            raise RuntimeError(f"{ds.name}: v0_4 scores {s}, v0_2 "
                               f"{v02.get(f'sc.{ds.name}')}: the tests "
                               "without electrons should not move")
        same += 1
    n["sc.same"] = str(same)
    rms = [score.kinetic_score(d, c) for d in data.by_role("train")
           if len(d.values) > 1]
    n["sc.train.lo"], n["sc.train.hi"] = f"{min(rms):.1f}", f"{max(rms):.1f}"
    for kr, meas in MO.XU_1993:
        n[f"xu.{kr:.0f}"] = MO._x(_krf_decay(kr, 1.0), meas)
        if n[f"xu.{kr:.0f}"] != v02[f"xu.{kr:.0f}"]:
            raise RuntimeError(f"Xu {kr} Torr: v0_4 {n[f'xu.{kr:.0f}']}, "
                               f"v0_2 {v02[f'xu.{kr:.0f}']}")
    # Tahoma under v0_2's closure: the measured numbers, v0_4 and its
    # variants, v0_3 in its own process
    meas = {d.name: d for d in data.by_role("device")}
    g0m = meas["tahoma_g0"].values[0]
    etam = meas["tahoma_formation"].values[0]
    isatm = meas["tahoma_isat"].values[0]
    alpha_m, net_m = 1.11, 5.28           # Table 3-3, beam stop in
    pure_m = 0.67                         # Table 3-2, 257.3 nm, pure Kr
    n.update({"tah.g0.meas": f"{g0m:.2f}", "tah.alpha.meas": f"{alpha_m:.2f}",
              "tah.net.meas": f"{net_m:.2f}", "tah.pure.meas": f"{pure_m:.2f}",
              "tah.eta.meas": f"{etam:.2f}", "tah.isat.meas": f"{isatm:.1f}"})
    rows = []
    for tag, kw in TAHOMA_VARIANTS.items():
        t = _tahoma_two_t(**kw)
        n[f"tah.{tag}.g0"] = f"{t['g0']:.2f}"
        n[f"tah.{tag}.g0.x"] = MO._x(t["g0"], g0m)
        n[f"tah.{tag}.alpha"] = f"{t['alpha']:.2f}"
        n[f"tah.{tag}.alpha.x"] = MO._x(t["alpha"], alpha_m)
        n[f"tah.{tag}.net"] = f"{t['g0'] - t['alpha']:.2f}"
        n[f"tah.{tag}.net.x"] = MO._x(t["g0"] - t["alpha"], net_m)
        n[f"tah.{tag}.te"] = f"{t['te']:.2f}"
        n[f"tah.{tag}.ne"] = f"${_k(t['ne'])}$"
        if tag == "default":
            by = sorted(t["by"].items(), key=lambda kv: -kv[1])
            n["tah.by"] = ", ".join(f"{_sp(s)} {v:.2f}" for s, v in by)
            n["tah.k2f.share"] = f"{100.0 * t['by']['Kr2F*'] / t['alpha']:.0f}"
            alpha_comp = (t["alpha"] - t["by"]["Kr2F*"]
                          + t["by"]["Kr2F*"] * SIGMA_K2F_COMPILATION
                          / 1.0e-18)
            n["tah.comp.alpha"] = f"{alpha_comp:.2f}"
            n["tah.comp.alpha.x"] = MO._x(alpha_comp, alpha_m)
    # pure krypton's loss, and v0_3 in its own process
    pure = _tahoma_two_t(x_f2=1e-9)
    n["tah.pure.alpha"] = f"{pure['alpha']:.2f}"
    n["tah.pure.x"] = MO._x(pure["alpha"], pure_m)
    n["tah.pure.by"] = ", ".join(f"{_sp(s)} {v:.2f}" for s, v in sorted(
        pure["by"].items(), key=lambda kv: -kv[1]) if v >= 0.005)
    v3 = _tahoma_other_version("v0_3")
    n["tah.v03.g0"] = f"{v3['g0']:.2f}"
    n["tah.v03.g0.x"] = MO._x(v3["g0"], g0m)
    n["tah.v03.alpha"] = f"{v3['alpha']:.2f}"
    n["tah.v03.alpha.x"] = MO._x(v3["alpha"], alpha_m)
    n["tah.v03.net"] = f"{v3['g0'] - v3['alpha']:.2f}"
    n["tah.v03.net.x"] = MO._x(v3["g0"] - v3["alpha"], net_m)
    n["tah.v03.te"] = f"{v3['te']:.2f}"
    n["tah.v03.ne"] = f"${_k(v3['ne'])}$"
    # the saturation intensity and the formation efficiency (ML's
    # observables) under v0_2's closure
    o = _isat_two_t()
    n.update({"tah.isat": f"{o['isat']:.2f}", "tah.isat.x": MO._x(o["isat"], isatm),
              "tah.eta": f"{o['eta']:.3f}", "tah.eta.x": MO._x(o["eta"], etam),
              "tah.form": f"{o['formation']:.3f}",
              "tah.ml.g0": f"{o['g0']:.2f}"})
    # the variants' table
    labels = {"default": "every band at its default",
              "kannari": r"the rung's bands at Kannari's values (\texttt{--ladder kannari})",
              "nocycle": r"the 5s $\leftrightarrow$ 5p cycle off (\texttt{k\_5s5p} = 0)",
              "kpslo": r"\Krpp{} + Kr at the band's low end",
              "krff": r"KrF(B,C) + F at Kannari's analogue",
              "kr2f": r"\Krts{} + F at the top of its band",
              "tri": r"the trimer on $T_g^{-1/2}$",
              "iiyield": r"ion-ion KrF\ensuremath{^*} yield at \wf{ii.yield.lo}",
              "iifc": r"ion-ion C share at \wf{ii.fc.hi}"}
    for tag in TAHOMA_VARIANTS:
        rows.append(f"{labels[tag]} & \\wf{{tah.{tag}.g0}} (\\x{{\\wf{{tah.{tag}.g0.x}}}}) & "
                    f"\\wf{{tah.{tag}.alpha}} (\\x{{\\wf{{tah.{tag}.alpha.x}}}}) & "
                    f"\\wf{{tah.{tag}.net}} (\\x{{\\wf{{tah.{tag}.net.x}}}}) & "
                    f"\\wf{{tah.{tag}.te}} & \\wf{{tah.{tag}.ne}}\\\\")
    tables = ["\\def\\wfTahVariants{" + "\n".join(rows) + "}"]
    n.update(_wave_numbers())
    n.update(_wall_numbers())
    n.update(_solver_checks())
    n["ib.checks"] = _ib_checks()
    n.update(_constants())
    return n, tables


# ---------------------------------------------------------------------------
# the propagated uncertainty (two_t): every input at its two ends
# ---------------------------------------------------------------------------
def _unc_inputs() -> list:
    """(kind, name, lo, hi): every band at its two ends, every row of
    UNC_ROWS at 1 -+ its quoted error."""
    inputs = [("band", name, b.lo, b.hi) for name, b in BANDS.items()]
    inputs += [("row", rid, 1.0 - d, 1.0 + d) for rid, (d, _) in UNC_ROWS.items()]
    return inputs


def _unc_kw(kind, name, v) -> dict:
    """The Params fields that set one input: a band through `rates`, the
    trimer's temperature law through its own field, a row through
    `row_scale`."""
    if kind == "row":
        return {"row_scale": ((name, v),)}
    if name == "trimer_tg_pow":
        return {"trimer_tg_pow": v}
    return {"rates": ((name, v),)}


def task_tahoma_two_t(**kw) -> dict:
    return _tahoma_two_t(**kw)


def _unc_stats(y0: float, pairs: dict) -> dict:
    """Each input's signed half-difference h = (y_hi - y_lo)/2, their
    quadrature sum, their linear sum, and the envelope (y0 plus the sum of
    every input's negative excursions, plus the sum of the positive
    ones)."""
    h = {k: 0.5 * (b - a) for k, (a, b) in pairs.items()}
    quad = math.sqrt(sum(v * v for v in h.values()))
    lin = sum(abs(v) for v in h.values())
    lo = y0 + sum(min(a - y0, b - y0, 0.0) for a, b in pairs.values())
    hi = y0 + sum(max(a - y0, b - y0, 0.0) for a, b in pairs.values())
    return {"h": h, "quad": quad, "lin": lin, "lo": lo, "hi": hi}


def unc_numbers() -> tuple:
    import make_system_map_v0_4 as MAP
    labels = dict(MAP.BAND_LABELS, **UNC_LABELS)
    inputs = _unc_inputs()
    t0 = time.time()
    # the reference point's 248-nm share, in this process
    _, _, d0 = _reference_two_t()
    ref0 = d0["emission_248"]
    ref = {}
    for kind, name, lo, hi in inputs:
        ref[name] = tuple(_reference_two_t(**_unc_kw(kind, name, v))[2]["emission_248"]
                          for v in (lo, hi))
    print(f"reference point: {2 * len(inputs)} runs in {time.time() - t0:.0f} s")
    # Tahoma under v0_2's closure, two processes
    specs = [(VERSION, "tahoma_two_t", {})]
    specs += [(VERSION, "tahoma_two_t", _unc_kw(kind, name, v))
              for kind, name, lo, hi in inputs for v in (lo, hi)]
    ctx = mp.get_context("spawn")
    t1 = time.time()
    with ctx.Pool(WORKERS, maxtasksperchild=1) as pool:
        results = pool.map(_run_task, specs, chunksize=1)
    print(f"Tahoma: {len(specs)} runs in {time.time() - t1:.0f} s on "
          f"{WORKERS} workers")
    tah0, rest = results[0], results[1:]
    tah = {}
    for i, (kind, name, lo, hi) in enumerate(inputs):
        tah[name] = (rest[2 * i], rest[2 * i + 1])
    g0_0, al_0 = tah0["g0"], tah0["alpha"]
    k2f_0 = tah0["by"]["Kr2F*"]
    pairs = {"ref": {k: v for k, v in ref.items()},
             "g0": {k: (a["g0"], b["g0"]) for k, (a, b) in tah.items()},
             "alpha": {k: (a["alpha"], b["alpha"]) for k, (a, b) in tah.items()},
             "net": {k: (a["g0"] - a["alpha"], b["g0"] - b["alpha"])
                     for k, (a, b) in tah.items()}}
    # the cross sections, analytic
    d_se, d_k2f = UNC_SIGMA["se"][0], UNC_SIGMA["k2f"][0]
    pairs["g0"]["se"] = (g0_0 * (1 - d_se), g0_0 * (1 + d_se))
    pairs["alpha"]["k2f"] = (al_0 - d_k2f * k2f_0, al_0 + d_k2f * k2f_0)
    pairs["net"]["se"] = (g0_0 * (1 - d_se) - al_0, g0_0 * (1 + d_se) - al_0)
    pairs["net"]["k2f"] = (g0_0 - al_0 - d_k2f * k2f_0, g0_0 - al_0 + d_k2f * k2f_0)
    y0 = {"ref": ref0, "g0": g0_0, "alpha": al_0, "net": g0_0 - al_0}
    meas = {"g0": 6.39, "alpha": 1.11, "net": 5.28}
    n = {"unc.n": str(len(inputs) + len(UNC_SIGMA)),
         "unc.nband": str(len(BANDS)), "unc.nrow": str(len(UNC_ROWS)),
         "unc.nsig": str(len(UNC_SIGMA)),
         "unc.sig.se": f"{100 * d_se:.0f}", "unc.sig.k2f": f"{100 * d_k2f:.0f}",
         "unc.sig.k2f.alpha": f"{100 * d_k2f * k2f_0 / al_0:.0f}",
         "unc.min": f"{(time.time() - t0) / 60:.0f}"}
    stats = {}
    for o in ("ref", "g0", "alpha", "net"):
        st = _unc_stats(y0[o], pairs[o])
        stats[o] = st
        fmt = "{:.2f}"
        n[f"unc.{o}.y0"] = fmt.format(y0[o])
        n[f"unc.{o}.quad"] = fmt.format(st["quad"])
        n[f"unc.{o}.quad.pct"] = f"{100 * st['quad'] / y0[o]:.0f}"
        n[f"unc.{o}.lin"] = fmt.format(st["lin"])
        n[f"unc.{o}.lin.pct"] = f"{100 * st['lin'] / y0[o]:.0f}"
        n[f"unc.{o}.lo"] = fmt.format(st["lo"])
        n[f"unc.{o}.hi"] = fmt.format(st["hi"])
        # the linear sums of the downward and upward excursions: an
        # overestimate of the spread (the inputs draw on the same KrF(B)
        # and a share cannot fall below 0), printed as such
        n[f"unc.{o}.down"] = f"{st['lo'] - y0[o]:+.2f}"
        n[f"unc.{o}.up"] = f"{st['hi'] - y0[o]:+.2f}"
        top = sorted(st["h"].items(), key=lambda kv: -abs(kv[1]))[:4]
        n[f"unc.{o}.top"] = ", ".join(
            f"{labels[k]} ($\\pm${abs(v):.2f})" for k, v in top)
        n[f"unc.{o}.top1"] = labels[top[0][0]]
        n[f"unc.{o}.top1.pct"] = f"{100 * abs(top[0][1]) / y0[o]:.0f}"
        if o in meas:
            dev = y0[o] - meas[o]
            n[f"unc.{o}.dev"] = f"{dev:+.2f}"
            n[f"unc.{o}.dev.pct"] = f"{100 * dev / meas[o]:+.0f}"
            n[f"unc.{o}.dev.quad"] = f"{abs(dev) / st['quad']:.1f}"
            # within one quadrature unit of the model, or beyond it
            n[f"unc.{o}.covers"] = ("within" if abs(dev) <= st["quad"]
                                   else "beyond")
    # the table: every input, its two ends, and the four outputs at them
    rows = []
    for kind, name, lo, hi in inputs:
        if kind == "band":
            b = BANDS[name]
            unit = MAP.UNITS.get(getattr(b, "unit", ""), "")
            rng = (f"${_k(lo) if lo else '0'}$ to ${_k(hi)}${unit}"
                   if abs(hi) < 1e-2 else f"{lo:.3g} to {hi:.3g}")
        else:
            rng = f"$\\pm{100 * UNC_ROWS[name][0]:.0f}$\\,\\% ({UNC_ROWS[name][1]})"
        cells = [labels[name], rng]
        for o in ("ref", "g0", "alpha", "net"):
            a, b = pairs[o][name]
            cells.append(f"{a:.2f} to {b:.2f}")
        rows.append(" & ".join(cells) + r"\\")
    for key in ("se", "k2f"):
        cells = [labels[key], f"$\\pm{100 * UNC_SIGMA[key][0]:.0f}$\\,\\% ({UNC_SIGMA[key][1]})", "--"]
        for o in ("g0", "alpha", "net"):
            if key in pairs[o]:
                a, b = pairs[o][key]
                cells.append(f"{a:.2f} to {b:.2f}")
            else:
                cells.append("--")
        rows.append(" & ".join(cells) + r"\\")
    tables = ["\\def\\wfUncTable{" + "\n".join(rows) + "}"]
    return n, tables


# ---------------------------------------------------------------------------
# the slow numbers: the kinetic EEDF, one process per run
# ---------------------------------------------------------------------------
def _bands_of(res) -> dict:
    m, e = res.model, res.ledger["E_dep"][-1]
    rows = {x.rid: x for x in m.rows}
    return {k: 100.0 * res.events[rid][-1] * rows[rid].hv_eV / e
            for k, rid in BAND_ROWS.items()}


def task_reference(**kw) -> dict:
    from v0_4.model import Model, Params
    from v0_4.report import accounting, reservoir, row_dict
    from v0_4 import optics
    t0 = time.time()
    m = Model(Params(eedf="kinetic", **kw))
    res = m.run()
    log = [e for e in res.kinetic if not e.get("skipped")]
    on = [e for e in log if e["on"]]
    late = [e for e in on if e["t"] >= 0.5 * m.p.pulse_s]
    krs = reservoir(res, "Kr*")
    kpp = reservoir(res, "Kr**")
    w = optics.window(res)
    return {"time": time.time() - t0, "bands": _bands_of(res),
            "row": {k: v for k, v in row_dict(res).items()
                    if isinstance(v, float)},
            "closure": accounting(res)["closure"],
            "solutions": len(log), "on": len(on),
            "en": (min(e["en_td"] for e in on if e["en_td"] > 0.0),
                   max(e["en_td"] for e in on)),
            "ratio": [e["power_ratio"] for e in late
                      if np.isfinite(e["power_ratio"])],
            "gap": max(abs(e["balance"]) for e in late),
            "brs76": 100.0 * krs["losses"]["K76"] / krs["formed"],
            "brpKPS": 100.0 * kpp["losses"]["KPS"] / kpp["formed"],
            "optics": {k: w[k] for k in ("g0", "alpha", "net")}}


def task_tahoma(**kw) -> dict:
    _ml()
    from ml import forward
    from v0_4 import optics
    t0 = time.time()
    p = forward.tahoma_params().replace(**kw)
    r = forward._drive().BeamModel(p).run()
    win = (r.t >= forward.WINDOW[0]) & (r.t <= forward.WINDOW[1])
    g0 = 100.0 * optics.SIGMA_SE_CM2 * float(np.mean(r.n["KrF(B)"][win]))
    by = {s: 100.0 * sig * float(np.mean(r.n[s][win]))
          for s, (sig, _) in optics.ABSORBERS.items() if s in r.n}
    log = [e for e in r.kinetic if not e.get("skipped")]
    inwin = [e for e in log if forward.WINDOW[0] <= e["t"] <= forward.WINDOW[1]]
    return {"g0": g0, "alpha": sum(by.values()), "by": by,
            "te": float(np.mean(r.te[win])), "ne": float(np.mean(r.ne[win])),
            "gap": max(abs(e["balance"]) for e in inwin),
            "reached": all(e["reached"] for e in inwin),
            "en": max(e["en_td"] for e in inwin),
            "time": time.time() - t0}


TASKS = {"reference": task_reference, "tahoma": task_tahoma,
         "tahoma_two_t": task_tahoma_two_t}


def _run_task(spec):
    version, name, kw = spec
    os.environ["LIGHT_VERSION"] = version
    return TASKS[name](**kw)


def _specs() -> list:
    return [(VERSION, "reference", {}), (VERSION, "tahoma", {})]


def slow_numbers() -> tuple:
    specs = _specs()
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(WORKERS, maxtasksperchild=1) as pool:
        results = pool.map(_run_task, specs, chunksize=1)
    print(f"{len(specs)} runs in {time.time() - t0:.0f} s on {WORKERS} "
          "workers")
    by = {(name, tuple(sorted(kw.items()))): r
          for (_, name, kw), r in zip(specs, results)}
    ref, tah = by[("reference", ())], by[("tahoma", ())]
    n = {}
    for k in BANDS_SHOWN:
        fmt = "{:.2f}" if k in ("248", "420", "IR") else (
            "{:.3f}" if k in ("272", "275") else "{:.1f}")
        n[f"ref.{k}"] = fmt.format(ref["bands"][k])
    for key, row in (("all", "emission_all"), ("tg", "tg_end_K"),
                     ("krstar", "to_Krstar"), ("krplus", "to_Krplus"),
                     ("krfb", "to_KrFB"), ("elastic", "elastic"),
                     ("heat", "heat_reactions"), ("te", "te_mean_eV")):
        fmt = "{:.0f}" if key == "tg" else ("{:.2f}" if key == "te"
                                              else "{:.1f}")
        n[f"ref.{key}"] = fmt.format(ref["row"][row])
    n["ref.closure"] = f"$10^{{{math.ceil(math.log10(abs(ref['closure']) + 1e-300))}}}$"
    n["ref.solutions"] = str(ref["solutions"])
    n["ref.on"] = str(ref["on"])
    n["ref.en.lo"], n["ref.en.hi"] = (f"{ref['en'][0]:.2f}",
                                      f"{ref['en'][1]:.0f}")
    n["ref.ratio.lo"] = f"{min(ref['ratio']):.3f}"
    n["ref.ratio.hi"] = f"{max(ref['ratio']):.3f}"
    n["ref.gap"] = _e(ref["gap"])
    n["ref.brs76"] = f"{ref['brs76']:.1f}"
    n["ref.brpKPS"] = f"{ref['brpKPS']:.1f}"
    n["ref.g0"] = f"{ref['optics']['g0']:.2f}"
    n["ref.alpha"] = f"{ref['optics']['alpha']:.2f}"
    n["ref.net"] = f"{ref['optics']['net']:.2f}"
    n["ref.min"] = f"{ref['time'] / 60:.0f}"
    n["tah.kin.g0"] = f"{tah['g0']:.2f}"
    n["tah.kin.g0.x"] = f"{tah['g0'] / 6.39:.2f}"
    n["tah.kin.alpha"] = f"{tah['alpha']:.2f}"
    n["tah.kin.alpha.x"] = f"{tah['alpha'] / 1.11:.2f}"
    n["tah.kin.net"] = f"{tah['g0'] - tah['alpha']:.2f}"
    n["tah.kin.net.x"] = f"{(tah['g0'] - tah['alpha']) / 5.28:.2f}"
    n["tah.kin.te"] = f"{tah['te']:.2f}"
    n["tah.kin.ne"] = f"${_k(tah['ne'])}$"
    n["tah.kin.gap"] = _e(tah["gap"])
    n["tah.kin.reached"] = "every" if tah["reached"] else "not every"
    n["tah.kin.en"] = f"{tah['en']:.2f}"
    n["tah.kin.min"] = f"{tah['time'] / 60:.0f}"
    by = sorted(tah["by"].items(), key=lambda kv: -kv[1])
    n["tah.kin.by"] = ", ".join(f"{_sp(s)} {v:.2f}" for s, v in by)
    return n, []


def _without_timing(text: str) -> str:
    """The file without its wall-clock keys (*.min: how long the runs
    took), which a comparison of two runs must not read."""
    return "\n".join(l for l in text.split("\n")
                     if not re.match(r"\\wfdef\{[^}]*\.min\}", l))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--fast", action="store_true",
                    help="write the fast numbers only")
    ap.add_argument("--slow", action="store_true",
                    help="write the slow numbers only (the fast file stays)")
    ap.add_argument("--unc", action="store_true",
                    help="write the propagated uncertainty only")
    ap.add_argument("--check", choices=("fast", "slow", "unc"),
                    help="compare one file with a fresh run, write nothing")
    a = ap.parse_args(argv)
    os.environ["LIGHT_VERSION"] = VERSION
    if a.check:
        if a.check == "fast":
            n, tables = fast_numbers()
            fresh = _render(n, "fast", tables)
        elif a.check == "unc":
            n, tables = unc_numbers()
            fresh = _render(n, "uncertainty", tables)
        else:
            fresh = _render(*_slow())
        path = {"fast": OUT_FAST, "slow": OUT_SLOW, "unc": OUT_UNC}[a.check]
        if _without_timing(path.read_text()) != _without_timing(fresh):
            print(f"{path.name} is stale: rerun make_overview_v0_4.py")
            return 1
        print(f"{path.name} is current")
        return 0
    only = a.fast or a.slow or a.unc
    if a.fast or not only:
        n, tables = fast_numbers()
        OUT_FAST.write_text(_render(n, "fast", tables))
        print(f"wrote {OUT_FAST}")
    if a.slow or not only:
        OUT_SLOW.write_text(_render(*_slow()))
        print(f"wrote {OUT_SLOW}")
    if a.unc or not only:
        n, tables = unc_numbers()
        OUT_UNC.write_text(_render(n, "uncertainty", tables))
        print(f"wrote {OUT_UNC}")
    return 0


def _slow():
    n, tables = slow_numbers()
    return n, "slow", tables


if __name__ == "__main__":
    sys.exit(main())
