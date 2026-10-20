"""Prototype e-beam drive for LIGHT + leave-one-out screen (evidence script, not LIGHT code).

Per ion pair (REBUILD spencer_fano_yields for Kr/0.5%F2, W = 27.07 eV):
Kr+ 0.9989, Kr** 0.5119 (prompt cascade to Kr* with 11.6-9.915 eV heat),
Kr* 0.2634, secondary-electron energy 4.506 eV into U_e. F2+/F (0.3 %)
folded into the electron share. Energy per pair is rebuilt from LIGHT's own
E_form so the ledger closes exactly.
"""
import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
import time
import numpy as np
from v0.model import Model, Params, QE
from v0.report import row_dict, accounting
KRP, KRS = 13.996, 9.915
Y_ION, Y_PP, Y_S, EPS_SEC = 0.9989, 0.5119, 0.2634, 4.506
CASCADE = Y_PP * (11.6 - KRS)
W_EFF = Y_ION * KRP + (Y_PP + Y_S) * KRS + CASCADE + EPS_SEC   # eV per pair

class BeamModel(Model):
    def rhs(self, t, y, on):
        if not on:
            return super().rhs(t, y, False)
        dy = super().rhs(t, y, False)          # chemistry with no electron heating
        N0 = self.N0
        pairs = self.p.power_W_cm3 / QE / W_EFF  # pairs cm^-3 s^-1
        dy[self.idx["Kr"]] -= pairs * (Y_ION + Y_PP + Y_S) / N0
        dy[self.idx["Kr+"]] += pairs * Y_ION / N0
        dy[self.idx["Kr*"]] += pairs * (Y_PP + Y_S) / N0
        dy[self.i_ue] += pairs * EPS_SEC / N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        tg_ev = y[self.i_tg] * 8.617333262e-5
        c_v = float((n * self.cv).sum()) * 8.617333262e-5
        dn_src = np.zeros(self.ns)
        dn_src[self.idx["Kr"]] = -pairs * (Y_ION + Y_PP + Y_S)
        dn_src[self.idx["Kr+"]] = pairs * Y_ION
        dn_src[self.idx["Kr*"]] = pairs * (Y_PP + Y_S)
        dy[self.i_tg] += (pairs * CASCADE - tg_ev * float((dn_src * self.cv).sum())) / c_v
        dy[self.i_led + 0] += pairs * W_EFF / N0          # E_dep
        dy[self.i_led + 2] += pairs * CASCADE / N0        # cascade heat (booked in the capture slot)
        return dy

    def kinetic_injection(self, t, y, on):
        """The secondary electrons, one per ion pair, born with the
        sub-excitation spectrum (linear on [0, KRS], mean EPS_SEC, as ML's
        solve_beam injects them). Read by a LIGHT version whose kinetic
        EEDF takes electron sources (v0_3, 2026-09-30). Earlier versions
        never call it."""
        if not on:
            return ()
        return ((self.p.power_W_cm3 / QE / W_EFF, KRS, EPS_SEC),)

def run_beam(p, drop=None):
    m = BeamModel(p)
    if drop:
        m.compiled = [c for c in m.compiled if c.row.rid.split("_C")[0] != drop and c.row.rid != drop]
    return m.run()

KEYS = ("emission_248", "to_KrFB", "heat_reactions", "stored_ions", "emission_all")
if __name__ == "__main__":
    base_p = Params(groups=("all",), ne_mode="free", ne_cm3=1e10, power_W_cm3=3e5)
    for pt in (dict(pressure_atm=1.0), dict(pressure_atm=2.0)):
        p = base_p.replace(**pt)
        r = run_beam(p)
        a = accounting(r); d = row_dict(r)
        print(f"\nE-BEAM {pt} closure {a['closure']:.1e}  n_e peak {r.ne.max():.2e}  Te pulse-mean {np.mean(r.te[r.t<=p.pulse_s]):.2f} eV")
        print("  " + "  ".join(f"{k}={d[k]:.1f}" for k in KEYS) + f"  KrF->248={d['krf_to_248_pct']:.1f}  Kr*->KrF={d['krstar_to_krf_pct']:.1f}  F2left={d['f2_left_pct']:.1f}")
        fb = {c.row.rid: float(r.events[c.row.rid][-1]) for c in r.model.compiled}
        krf_src = {rid: v for rid, v in fb.items() if rid in ("G3", "NK7", "K2F_F", "K26", "K108")}
        tot = sum(krf_src.values())
        print("  KrF(B) formation by channel: " + ", ".join(f"{k} {100*v/tot:.1f}%" for k, v in sorted(krf_src.items(), key=lambda kv: -kv[1])))
    # leave-one-out at 1 atm under the beam drive
    p = base_p.replace(pressure_atm=1.0)
    base = row_dict(run_beam(p))
    rids = sorted({c.row.rid.split("_C")[0] if (c.row.rid.endswith("_C") and c.row.rid != "URAD_C") else c.row.rid
                   for c in BeamModel(p).compiled})
    out = []
    t0 = time.time()
    for rid in rids:
        try:
            d = row_dict(run_beam(p, rid))
            out.append((abs(d["emission_248"] - base["emission_248"]), rid, d["emission_248"] - base["emission_248"], d["to_KrFB"] - base["to_KrFB"], d["heat_reactions"] - base["heat_reactions"]))
        except Exception as exc:
            out.append((999, rid, float("nan"), float("nan"), float("nan")))
    print(f"\nLOO under the e-beam drive at 1 atm (base 248 = {base['emission_248']:.2f} %, {time.time()-t0:.0f} s)")
    for m, rid, d248, dkrf, dheat in sorted(out, reverse=True)[:22]:
        print(f"  {rid:8s} d248 {d248:+6.2f}  dKrF {dkrf:+6.2f}  dheat {dheat:+6.2f}")
