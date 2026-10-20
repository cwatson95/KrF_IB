"""Prototype: explicit Kr** (5p) with electron up-pumping from Kr* (evidence script, not LIGHT code).

Rows (group 'kstar'), values as printed:
  KXP   e + Kr*  -> e + Kr**    Kannari 1985 Table I: 9.00e-7 Te^0.50 exp(-1.51/Te) [76]
  KXPr  e + Kr** -> e + Kr*     J&H 1980 #50: 8.0e-7 (typical, reverse)
  KPPI  e + Kr** -> Kr+ + 2e    Kannari 1985 Table I: 2.00e-7 Te^0.62 exp(-2.53/Te) [79]
  KPS   Kr** + Kr -> Kr* + Kr   REBUILD KPS_Kr class value 1e-10 (EST; Janacaitis list 3.0e-11)
  KRAD  Kr** -> Kr* + hv        Chow 1978 #18: tau 27.6 ns
  K27   Kr** + F2 -> KrF(B) + F REBUILD K27 8.1e-10
Kr** placed 1.51 eV above Kr* (Kannari's gap) so thresholds and stores agree.
"""
import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
import math
from v0 import species, network
from v0 import model as vmodel
from v0.species import Species
from v0.network import Row
species.SPECIES["Kr**"] = Species("Kr**", 9.915 + 1.51, (("Kr", 1),), 0, 1.5, "prototype: Kannari 5s-5p gap")
_orig_opt = network.optional_rows
def optional_rows():
    return _orig_opt() + [
        Row("KXP", ("e", "Kr*"), ("e", "Kr**"), ("te_fit", 9.0e-7, 0.5, 1.51), "kstar", "Kannari [76]"),
        Row("KXPr", ("e", "Kr**"), ("e", "Kr*"), ("const", 8.0e-7), "kstar", "J&H #50", energy_to="electrons"),
        Row("KPPI", ("e", "Kr**"), ("Kr+", "e", "e"), ("te_fit", 2.0e-7, 0.62, 2.53), "kstar", "Kannari [79]"),
        Row("KPS", ("Kr**", "Kr"), ("Kr*", "Kr"), ("const", 1.0e-10), "kstar", "REBUILD KPS_Kr EST"),
        Row("KRAD", ("Kr**",), ("Kr*",), ("const", 1 / 27.6e-9), "kstar", "Chow #18", hv_eV=1.51, band="IR (Kr 5p-5s)"),
        Row("K27", ("Kr**", "F2"), ("KrF(B)", "F"), ("const", 8.1e-10), "kstar", "REBUILD K27"),
    ]
network.optional_rows = optional_rows
network.GROUPS = network.GROUPS + ("kstar",)
_orig_k = vmodel.Model._k
def _k(self, c, te, tg, n, ne):
    if c.kind == "te_fit":
        t = max(te, 0.026)
        return c.rate[0] * t ** c.rate[1] * math.exp(-c.rate[2] / t)
    return _orig_k(self, c, te, tg, n, ne)
vmodel.Model._k = _k
_orig_compile = vmodel.Model._compile
def _compile(self, j, r):
    if r.rate[0] == "te_fit":
        from v0.network import energy_routing
        from v0.model import Compiled
        delta = {}
        for s in r.reactants:
            if s != "e": delta[s] = delta.get(s, 0) - 1
        for s in r.products:
            if s != "e": delta[s] = delta.get(s, 0) + 1
        stoich = tuple((self.idx[s], cc) for s, cc in delta.items() if cc != 0)
        return Compiled(r, j, tuple(self.idx[s] for s in r.heavy), r.n_electrons_in, stoich, energy_routing(r), "te_fit", r.rate[1:])
    return _orig_compile(self, j, r)
vmodel.Model._compile = _compile

from v0.model import Params, run
from v0.report import row_dict, accounting
from ebeam_screen import run_beam
BASE = ("quench", "kr2", "ions", "electron", "bc", "f2e")
def show(label, d, a=None):
    print(f"  {label:28s} Te {d['te_mean_eV']:.2f}  248 {d['emission_248']:5.2f}  KrF(B) {d['to_KrFB']:5.1f}  Kr+ {d['to_Krplus']:5.1f}  heat {d['heat_reactions']:5.1f}  elastic {d['elastic']:4.1f}" + (f"  closure {a['closure']:.0e}" if a else ""))
for ne in (1e15, 1e16):
    print(f"\npower drive, 1 atm, n_e {ne:.0e} held")
    for label, g in (("all groups", BASE), ("all groups + Kr** (5p)", BASE + ("kstar",))):
        r = run(Params(groups=g, ne_cm3=ne)); show(label, row_dict(r), accounting(r))
print("\ne-beam prototype, 1 atm, free n_e")
for label, g in (("all groups", BASE), ("all groups + Kr** (5p)", BASE + ("kstar",))):
    r = run_beam(Params(groups=g, ne_mode="free", ne_cm3=1e10)); show(label, row_dict(r), accounting(r))
