import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
from ebeam_screen import BeamModel
from v0.model import Model, Params
from v0.network import GROUPS
from v0.report import row_dict
def run_with(p, beam, keq5):
    m = (BeamModel if beam else Model)(p)
    if keq5:
        for c in m.compiled:
            if c.row.rid == "CB_Kr": c.rate = 5.0e-10 / 5.0
    return row_dict(m.run())
for beam in (False, True):
    for pres in (1.0, 2.0):
        base = dict(groups=GROUPS, pressure_atm=pres) if not beam else dict(groups=GROUPS, pressure_atm=pres, ne_mode="free", ne_cm3=1e10)
        a = run_with(Params(**base), beam, False)
        b = run_with(Params(**base), beam, True)
        drive = "e-beam" if beam else "power "
        print(f"{drive} {pres} atm, all groups:  K_eq 1.47 -> 248 {a['emission_248']:5.2f} %, KrF->248 {a['krf_to_248_pct']:5.1f} %   |   K_eq 5 -> 248 {b['emission_248']:5.2f} %, KrF->248 {b['krf_to_248_pct']:5.1f} %")
