import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
from ebeam_screen import run_beam
from v0.model import Params, run
from v0.report import row_dict
NO_RECOMB = ("quench", "kr2", "electron", "bc")          # everything except ions + f2e
def go(drive, groups, pres):
    for method in ("BDF", "Radau", "LSODA"):
        try:
            if drive == "e-beam":
                return row_dict(run_beam(Params(groups=groups, ne_mode="free", ne_cm3=1e10, pressure_atm=pres, method=method))), method
            return row_dict(run(Params(groups=groups, pressure_atm=pres, method=method))), method
        except Exception as exc:
            last = repr(exc)[:80]
    return None, last
for drive in ("e-beam", "power"):
    for pres in (0.5, 1.0, 2.0):
        a, ma = go(drive, NO_RECOMB, pres)
        b, mb = go(drive, NO_RECOMB + ("ions", "f2e"), pres)
        if a is None or b is None:
            print(f"{drive:6s} {pres} atm: failed ({ma} / {mb})"); continue
        print(f"{drive:6s} {pres:3.1f} atm  without recombination: KrF(B) {a['to_KrFB']:5.1f}  248 {a['emission_248']:5.2f}  ions {a['stored_ions']:5.1f} | "
              f"with: KrF(B) {b['to_KrFB']:5.1f}  248 {b['emission_248']:5.2f}  ions {b['stored_ions']:5.1f}   [{ma}/{mb}]")
