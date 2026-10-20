"""KrF formation as the target: the chain from the beam's Kr+ and Kr* to
KrF(B,C), each link pinned to measurements, and the quantity Tahoma measured.

The observable. Tahoma's formation efficiency is g0 Isat / Pin = 0.27
(Spectra Technology 1989, OSTI 10136715, p. 3-16), with Isat = 1.4 +- 0.1
MW/cm2 from a Rigrod fit to its extraction data. LIGHT's own formation rate
x 5 eV / Pin is a different number: KrF that forms in, or mixes into, the C
state and decays there never meets the laser field. `observables` puts a
uniform 248 nm intracavity flux I on the model (stimulated emission
sigma_B I / h nu from KrF(B)) and fits g(I) = g0 / (1 + I/Isat) over
0.5-4 MW/cm2 with g0 held at the small-signal gain, as Tahoma held its
measured g0. g0, Isat and g0 Isat / Pin then compare like for like, and
Isat is a second measured number the model must meet.

The chain (Chain). Each link has LIGHT's value as the default (since
2026-09-28 LIGHT carries the measured links, step 2 of the owner's plan):
  k_H      Kr* + F2, total rate                 8.47e-10 (three first-hand)
  gamma_H  KrF* share of Kr* + F2               0.90 (Velazco et al. 1976)
  f_C_H    share of that KrF formed in C        0.38 (Yu & Setser 1990)
  k3       Kr* + 2Kr -> Kr2* at 300 K           1.7e-32 (Ueda 1984)
  k3_tpow  its (Tg/300)^n                       -0.64 (Wieme & Lenaerts)
  kr2_rad  Kr2* radiative rate                  3.7e6 1/s (270 ns)
  k_kr2_f2 Kr2* + F2, total                     2.1e-10 (Quigley & Hughes)
  kr2_krf  KrF* share of Kr2* + F2              0.5 (Kannari, unmeasured)
  k_kr2_f  Kr2* + F -> KrF(B) + Kr              3.0e-10 (Kannari), and 0
           from LIGHT v0_2, which carries it as a band (unmeasured)
  ion_scale ion-ion recombination x             LIGHT 1 (Flannery form)
  y_ion    KrF* share of ion-ion recombination  LIGHT 1
  f_C_ion  share of that KrF formed in C        LIGHT 0
  thermal  thermal-electron Kr excitation and ionization (G1, C1, C2) on
  thermal_upper  only the 5p/4d/5d excitation and ionization (C1a-c, C2):
           the kinetic EEDF empties the tail above 11 eV (ml/eedf.py), and
           in pure Kr at low pump LIGHT's thermal excitation is not physical
           at all (ml/pure.py)
The measured values and their sources are in ml/formation_data.py.

g0 Isat / Pin is a formation number: the fitted and LIGHT's B/C constants
give 0.328 and 0.329 although their g0 differ x2. The model's saturation is
not quite homogeneous. Lasing turns part of the pump into light instead of
heat, the gas runs cooler, and the KrF(B) + 2Kr trimer rate (Tg^-3.3) rises.
A real laser does the same, so Isat carries it.

The B/C constants are ml/params.Consts (the fitted set or LIGHT's), put on
LIGHT's rows by forward.impose. The chain's rows are added to LIGHT's own
network at build time. LIGHT itself is unchanged.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, replace

import numpy as np

from . import LIGHT_DIR, light
from . import forward
from .forward import HV_EV, SIGMA_B

HV_PROBE = 1239.84 / 248.4       # eV, Tahoma's probe line
FLUXES = (0.5, 1.0, 2.0, 4.0)    # MW/cm2, the Isat fit (Tahoma's range)
TAHOMA = {"g0": 6.39, "isat": 1.4, "eta": 0.27}   # %/cm, MW/cm2, 1
MAN = ("KrF(B)", "KrF(C)")


def _light_default(band, value):
    """The LIGHT version's own default for a link it carries as a band
    (network.BANDS, v0_2 on), else `value`, v0's and v0_1's."""
    bands = getattr(light("network"), "BANDS", {})
    return bands[band].default if band in bands else value


@dataclass(frozen=True)
class Chain:
    k_H: float = 8.47e-10
    gamma_H: float = 0.90
    f_C_H: float = 0.38
    k3: float = 1.7e-32
    k3_tpow: float = -0.64
    kr2_rad: float = 3.7e6
    k_kr2_f2: float = 2.1e-10
    kr2_krf: float = 0.5
    k_kr2_f: float = _light_default("kr2_f", 3.0e-10)
    ion_scale: float = 1.0
    y_ion: float = 1.0
    f_C_ion: float = 0.0
    thermal: bool = True
    thermal_upper: bool = True

    def with_(self, **kw):
        return replace(self, **kw)


LIGHT_CHAIN = Chain()


def _rows(chain, I_MW):
    """The chain's extra rows (LIGHT's Row, checked by LIGHT's own
    element, charge and energy rules). The ion-ion shares copy the rate
    form of LIGHT's own K26 and K108 (v0_1 adds a Tg law to it)."""
    network = light("network")
    Row, check_row = network.Row, network.check_row
    own = {r.rid: r.rate for r in network.build(
        forward.tahoma_params().groups, "basic")}
    rows = [
        Row("G3X", ("Kr*", "F2"), ("Kr", "F", "F"), ("const", 0.0),
            "formation", "Kr* + F2 without KrF*: the (1 - gamma_H) share"),
        Row("G3C", ("Kr*", "F2"), ("KrF(C)", "F"), ("const", 0.0),
            "formation", "Kr* + F2 -> KrF(C): the f_C_H share"),
        Row("K26X", ("Kr+", "F-"), ("Kr", "F"), own["K26"],
            "formation", "Kr+ + F- without KrF*: the (1 - y_ion) share"),
        Row("K26C", ("Kr+", "F-"), ("KrF(C)",), own["K26"],
            "formation", "Kr+ + F- -> KrF(C): the f_C_ion share"),
        Row("K108X", ("Kr2+", "F-"), ("Kr", "Kr", "F"), own["K108"],
            "formation", "Kr2+ + F- without KrF*: the (1 - y_ion) share"),
        Row("K108C", ("Kr2+", "F-"), ("KrF(C)", "Kr"), own["K108"],
            "formation", "Kr2+ + F- -> KrF(C) + Kr: the f_C_ion share"),
    ]
    if I_MW > 0:
        k = SIGMA_B * I_MW * 1e6 / (HV_PROBE * forward.QE_SI)
        rows.append(Row("STIM", ("KrF(B)",), ("Kr", "F"), ("const", k),
                        "probe", "248 nm probe flux: sigma_B I / h nu",
                        hv_eV=HV_PROBE, band="248 nm stimulated"))
    for r in rows:
        check_row(r)
    return rows


def build(c, chain=LIGHT_CHAIN, I_MW=0.0, cond=None):
    """LIGHT's Tahoma model (the forward module's drive and conditions)
    with the B/C set c and the formation chain."""
    BeamModel = forward._drive().BeamModel
    M = light("model")
    p = forward.tahoma_params(cond, kr2=chain.k3 > 0)
    extra = _rows(chain, I_MW)
    orig = M.build

    def build_with(g, n):
        base = orig(g, n)
        have = {r.rid for r in base}
        return base + [r for r in extra if r.rid not in have]
    M.build = build_with
    try:
        m = BeamModel(p)
    finally:
        M.build = orig
    forward.impose(m, c)
    _impose_chain(m, chain)
    return m


def _impose_chain(m, ch):
    by = {cc.row.rid: cc for cc in m.compiled}

    def flan(rid, share):
        cc = by[rid]
        cc.rate = ((cc.row.rate[1] * share * ch.ion_scale,)
                   + tuple(cc.row.rate[2:]))
    by["G3"].rate = ch.k_H * ch.gamma_H * (1.0 - ch.f_C_H)
    by["G3C"].rate = ch.k_H * ch.gamma_H * ch.f_C_H
    by["G3X"].rate = ch.k_H * (1.0 - ch.gamma_H)
    for main in ("K26", "K108"):
        flan(main, ch.y_ion * (1.0 - ch.f_C_ion))
        flan(main + "C", ch.y_ion * ch.f_C_ion)
        flan(main + "X", 1.0 - ch.y_ion)
    if "AK8_Kr" in by:          # v0, v0_1: the kr2 group; v0_2: basic
        by["AK8_Kr"].rate = ch.k3
        by["AK8_Kr"].row = replace(by["AK8_Kr"].row, tg_pow=ch.k3_tpow)
        by["KR2_RAD"].rate = ch.kr2_rad
        by["NK7"].rate = ch.k_kr2_f2 * ch.kr2_krf
        by["K2F_BR"].rate = ch.k_kr2_f2 * (1.0 - ch.kr2_krf)
        by["K2F_F"].rate = ch.k_kr2_f
    drop = set()
    if not ch.thermal:
        drop |= {"G1", "C1a", "C1b", "C1c", "C2"}
    if not ch.thermal_upper:
        drop |= {"C1a", "C1b", "C1c", "C2"}
    if drop:
        m.compiled = [cc for cc in m.compiled if cc.row.rid not in drop]


def run(c, chain=LIGHT_CHAIN, I_MW=0.0, cond=None):
    """One run: window means of the gain and of the formation budget."""
    m = build(c, chain, I_MW, cond)
    r = m.run()
    t = r.t
    win = (t >= forward.WINDOW[0]) & (t <= forward.WINDOW[1])
    pairs = forward.pair_rate(m.p)

    def slope(rid):
        return float(np.polyfit(t[win], r.events[rid][win], 1)[0])
    by = {}
    for cc in m.compiled:
        net = (sum(cc.row.products.count(s) for s in MAN)
               - sum(cc.row.reactants.count(s) for s in MAN))
        if net > 0:
            by[cc.row.rid] = net * slope(cc.row.rid) / pairs
    nb = float(np.mean(r.n["KrF(B)"][win]))
    nc = float(np.mean(r.n["KrF(C)"][win]))
    return {"g": SIGMA_B * nb * 100.0, "nb": nb, "nc": nc,
            "tg": float(np.mean(r.tg[win])),
            "per_pair": sum(by.values()), "by": by,
            "formation": sum(by.values()) * pairs * HV_EV * forward.QE_SI
            / m.p.power_W_cm3,
            "events": {cc.row.rid: slope(cc.row.rid) / pairs
                       for cc in m.compiled},
            "pin": m.p.power_W_cm3, "pairs": pairs}


def observables(c, chain=LIGHT_CHAIN, cond=None, fluxes=FLUXES):
    """g0 [%/cm], Isat [MW/cm2], g0 Isat / Pin, the total formation and its
    budget per ion pair. As in Tahoma's analysis g0 is the small-signal gain
    and Isat alone is fitted to the saturated gains: g0/g(I) - 1 = I/Isat,
    least squares through the origin. isat_points gives each flux's own
    value (a homogeneous medium would give one number)."""
    base = run(c, chain, 0.0, cond)
    g0 = base["g"]
    y = np.array([g0 / run(c, chain, I, cond)["g"] - 1.0 for I in fluxes])
    x = np.array(fluxes)
    isat = float(np.sum(x * x) / np.sum(x * y))
    eta = g0 / 100.0 * isat * 1e6 / base["pin"]
    return dict(base, g0=g0, isat=isat, eta=eta, isat_points=x / y)
