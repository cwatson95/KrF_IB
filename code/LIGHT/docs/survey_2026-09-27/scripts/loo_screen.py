"""Leave-one-out importance screen on LIGHT v0 + all optional groups (prototype, not repo code)."""
import sys
from pathlib import Path
_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_LIGHT))
import time
from v0.model import Model, Params
from v0.report import row_dict
POINTS = [dict(pressure_atm=0.1, ne_cm3=1e15), dict(pressure_atm=1.0, ne_cm3=1e14),
          dict(pressure_atm=1.0, ne_cm3=1e15), dict(pressure_atm=1.0, ne_cm3=1e16),
          dict(pressure_atm=2.0, ne_cm3=1e15)]
KEYS = ("emission_248", "to_KrFB", "heat_reactions", "stored_ions", "emission_all")
def run_without(p, drop):
    m = Model(p)
    if drop:
        m.compiled = [c for c in m.compiled if c.row.rid.split("_C")[0] != drop and c.row.rid != drop]
    return row_dict(m.run())
base_rows = [r.rid for r in Model(Params(groups=("all",))).rows]
rids = sorted({r.split("_C")[0] if r.endswith("_C") and r not in ("URAD_C",) else r for r in base_rows})
t0 = time.time()
results = {}
for i, pt in enumerate(POINTS):
    p = Params(groups=("all",), **pt)
    base = run_without(p, None)
    results[f"p{i}"] = {"point": pt, "base": {k: base[k] for k in KEYS}, "loo": {}}
    for rid in rids:
        try:
            d = run_without(p, rid)
            results[f"p{i}"]["loo"][rid] = {k: d[k] - base[k] for k in KEYS}
        except Exception as exc:
            results[f"p{i}"]["loo"][rid] = {"error": repr(exc)[:120]}
    print(f"point {i} {pt} base 248={base['emission_248']:.2f} done at {time.time()-t0:.0f}s", flush=True)
# rank rows by max |d 248| over points
rank = []
for rid in rids:
    vals = [results[f"p{i}"]["loo"][rid] for i in range(len(POINTS))]
    if any("error" in v for v in vals):
        rank.append((999, rid, "error")); continue
    m248 = max(abs(v["emission_248"]) for v in vals)
    mkrf = max(abs(v["to_KrFB"]) for v in vals)
    mheat = max(abs(v["heat_reactions"]) for v in vals)
    rank.append((m248, rid, f"max|d248| {m248:5.2f}  max|dKrF| {mkrf:5.2f}  max|dheat| {mheat:5.2f}  " + " ".join(f"{v['emission_248']:+.2f}" for v in vals)))
for m, rid, s in sorted(rank, reverse=True):
    print(f"{rid:8s} {s}")
