"""The electron side of LIGHT: cross sections and the EEDF.

LIGHT integrates its LoKI-MC cross sections (Kr: Biagi-v7.1, F2: SIGLO,
md5-pinned copies in LIGHT/v0/data) over a Maxwellian EEDF at Te. This
module changes that side at run time, for three purposes:
  1. Scale a cross section (sensitivity).
  2. Swap F2 to the TRINITI set that LoKI-MC also ships (a second, independent
     F2 set, and so a measure of how far the F2 cross sections are known).
  3. Replace a row's table with one from a trained cross section
     (att_table: the F2 attachment as ml/swarm.py fits it).
The EEDF itself is ml/eedf.py's business: it scales LIGHT's rows by the
ratio of the beam-sustained Boltzmann EEDF's rates to the Maxwellian's.

What may train these is set by the same policy as the rest of ML. Measured
swarm data (drift velocity, ionization and attachment coefficients) are
measurements of these very cross sections and EEDFs. Laser data only score.
"""
from __future__ import annotations

import numpy as np

from . import ROOT, light

xsec = light("xsec")

LOKI_INPUT = ROOT / "LoKI-MC" / "Code" / "Input"
TRINITI = LOKI_INPUT / "Fluorine2" / "F2_LXCat_TRINITI.txt"

# LIGHT's electron rows and the cross sections behind them
F2_ROWS = {"ATT": ("attachment", None), "DIS_a3": ("excitation", "F2(a3)(3.16eV)"),
           "DIS_A1": ("excitation", "F2(A1)(4.34eV)"),
           "VIB1": ("excitation", "F2(V1)(0.111eV)"),
           "VIB2": ("excitation", "F2(V2)(0.219eV)"),
           "VIB3": ("excitation", "F2(V3)(0.324eV)"),
           "VIB4": ("excitation", "F2(V4)(0.421eV)")}
KR_ROWS = {"G1": "Kr(s)(9.915eV)", "C1a": "Kr(p)(11.30eV)",
           "C1b": "Kr(d+p)(11.998eV)", "C1c": "Kr*(12.75eV)", "C2": "Kr^+"}
GROUPS = {"F2 attachment": ("ATT",), "F2 vibrational": ("VIB1", "VIB2",
                                                        "VIB3", "VIB4"),
          "F2 dissociative excitation": ("DIS_a3", "DIS_A1"),
          "Kr elastic (energy loss)": ("kel",),
          "Kr excitation": ("G1", "C1a", "C1b", "C1c"), "Kr ionization": ("C2",)}


def _table(proc):
    k = xsec._mean(proc.sigma_at(xsec.EPS) * xsec.V_CM_S * np.sqrt(xsec.EPS))
    return xsec.RateTable(k, proc.product)


def triniti_tables() -> dict:
    """LIGHT's F2 rows filled from TRINITI. TRINITI lumps the vibrational
    levels into one 0.11-eV process, which goes into VIB1, and VIB2-4 get
    zero."""
    ps = xsec.parse_lxcat(TRINITI)
    by = {(p.kind, p.product): p for p in ps}
    zero = xsec.RateTable(np.zeros_like(xsec.TE), "none")
    return {"ATT": _table(by[("attachment", "F + F(-)")]),
            "VIB1": _table(by[("excitation", "F2(VIB)")]),
            "VIB2": zero, "VIB3": zero, "VIB4": zero,
            "DIS_a3": _table(by[("excitation", "F2(*1)")]),
            "DIS_A1": _table(by[("excitation", "F2(*2)")])}


def scaled(table, factor):
    return xsec.RateTable(table.values * factor, f"{table.label} x{factor}")


def att_table(s_low=1.0, s_high=1.0):
    """LIGHT's ATT row (Maxwellian) with the attachment cross section's two
    resonances rescaled as ml/swarm.py trains them."""
    from .swarm import att_factor
    p = xsec.find("F2", "attachment", None)
    k = xsec._mean(p.sigma_at(xsec.EPS) * att_factor(xsec.EPS, s_low, s_high)
                   * xsec.V_CM_S * np.sqrt(xsec.EPS))
    return xsec.RateTable(k, f"F2 attachment x({s_low:.3g}, {s_high:.3g})")


def modifier(scales=None, f2_set="SIGLO", tables=None, capture=None):
    """A function that changes a compiled LIGHT model's electron side.
    scales: {row id or 'kel': factor}, applied last. f2_set: 'SIGLO' or
    'TRINITI'. tables: {row id: RateTable} replacing a row's Maxwellian
    table (a trained cross section, att_table). capture: {row id: factor}
    on the energy LIGHT books for each electron the row removes, (3/2)Te on
    its Maxwellian closure (forward.MAXWELL)."""
    scales = scales or {}
    tables = tables or {}
    capture = capture or {}
    tri = triniti_tables() if f2_set == "TRINITI" else None

    def apply(m):
        for cc in m.compiled:
            rid = cc.row.rid
            if cc.kind != "xsec":
                continue
            table = cc.rate
            if tri is not None and rid in tri:
                table = tri[rid]
            if rid in tables:
                table = tables[rid]
            if rid in scales:
                table = scaled(table, scales[rid])
            cc.rate = table
        for cc in m.compiled:
            if cc.row.rid in capture:
                cc.route = dict(cc.route, capture=cc.route["capture"]
                                * capture[cc.row.rid])
        if "kel" in scales:
            m.kel = scaled(m.kel, scales["kel"])
    return apply
