"""Run LIGHT itself against the Xu 1993 / Gadomski 1992 photolysis decays (seed mode, no power, no electrons)."""
import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
import numpy as np
from v0.model import Model, Params
TORR = 133.322368 / (1.380649e-23 * 300.0) * 1e-6        # cm^-3 per Torr at 300 K
CASES = [(50, 1, 5.0e7, "Xu p.2594"), (100, 1, 6.2e7, "Xu p.2594"), (500, 1, 2.63e8, "Gadomski fit curve (agent read)")]
def decay_const(kr_torr, f2_torr, groups, k_bc=None, k_cb=None, t_end=400e-9):
    p = Params(pressure_atm=(kr_torr + f2_torr) / 760.0, x_f2=f2_torr / (kr_torr + f2_torr),
               power_W_cm3=0.0, pulse_s=0.0, t_end_s=t_end, ne_mode="free", ne_cm3=1.0,
               groups=groups, init_extra_cm3=(("KrF(B)", 1e12),), n_save=2001)
    m = Model(p)
    if k_bc is not None:
        for c in m.compiled:
            if c.row.rid == "BC_Kr": c.rate = k_bc
            if c.row.rid == "CB_Kr": c.rate = k_cb
    r = m.run()
    tot = r.n["KrF(B)"] + (r.n["KrF(C)"] if "KrF(C)" in r.n else 0.0)
    rel = tot / tot[0]
    sel = (r.t > 20e-9) & (rel < 1e-2) & (rel > 1e-6)
    slope = np.polyfit(r.t[sel], np.log(tot[sel]), 1)[0]
    return -slope
print("Kr/F2 Torr   measured   v0          v0+quench    v0+quench+bc   +bc with K_eq 5 (k_CB = k_BC/5)   Hakuta k_BC 1.2e-10, K_eq 5")
for kr, f2, meas, src in CASES:
    row = [decay_const(kr, f2, g) for g in ((), ("quench",), ("quench", "bc"))]
    keq5 = decay_const(kr, f2, ("quench", "bc"), k_bc=5.0e-10, k_cb=1.0e-10)
    hak = decay_const(kr, f2, ("quench", "bc"), k_bc=1.2e-10, k_cb=2.4e-11)
    print(f"{kr:3d}/{f2}   {meas:.2e}  " + "  ".join(f"{v:.2e} (x{v/meas:.2f})" for v in row + [keq5, hak]))
