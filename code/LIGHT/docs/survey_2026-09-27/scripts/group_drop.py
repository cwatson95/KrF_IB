import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
from ebeam_screen import run_beam, BeamModel
from v0.model import Params, run
from v0.network import GROUPS
from v0.report import row_dict
def fmt(d): return f"248 {d['emission_248']:5.1f}  KrF(B) {d['to_KrFB']:5.1f}  heat {d['heat_reactions']:5.1f}  ions {d['stored_ions']:5.1f}  F2left {d['f2_left_pct']:5.1f}"
for drive in ("e-beam", "power"):
    for pres in (1.0, 2.0):
        print(f"\n== {drive} drive, {pres} atm, 300 kW/cm3, 1 us, Kr + 0.5 % F2")
        for label, groups in [("v0 only", ())] + [(f"all minus {g}", tuple(x for x in GROUPS if x != g)) for g in GROUPS] + [("all", GROUPS)]:
            try:
                if drive == "e-beam":
                    p = Params(groups=groups, ne_mode="free", ne_cm3=1e10, pressure_atm=pres)
                    d = row_dict(run_beam(p))
                else:
                    p = Params(groups=groups, pressure_atm=pres)
                    d = row_dict(run(p))
                print(f"  {label:20s} {fmt(d)}")
            except (ValueError, RuntimeError) as exc:
                # expected once: free n_e + attachment drives n_e to its
                # floor and U_e/(1.5 n_e) turns the Jacobian non-finite
                print(f"  {label:20s} solver failed: {str(exc)[:60]}")
