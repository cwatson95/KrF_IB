"""ML runs the LIGHT version the environment names (ml/__init__.py).

The version is fixed when ml is first imported, so each check probes a
fresh process."""

import json
import os
import subprocess
import sys
from pathlib import Path

ML_DIR = Path(__file__).resolve().parents[1]

PROBE = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
import ml
from ml import eedf, formation, forward, params, pure
p = forward.tahoma_params()
drive = forward._drive()
c = forward.with_rederived_k_e(params.light_now())
m = forward.impose(drive.BeamModel(p), c)
rad = {cc.row.rid: cc.rate for cc in m.compiled
       if cc.row.rate[:2] == ("inv_tau", "tau_sp_s")}
rows = {r.rid: r.rate for r in formation._rows(formation.LIGHT_CHAIN, 0.0)}
chain = {cc.row.rid: cc.rate for cc in formation.build(c).compiled}
kr2 = {k: chain.get(k) for k in ("AK8_Kr", "K2F_F", "PEN_KR2")}
refused = False
if ml.LIGHT_VERSION != "v0":
    from ml import tail
    try:
        tail.write_table([])
    except RuntimeError:
        refused = True
print(json.dumps({"version": ml.LIGHT_VERSION, "groups": list(p.groups),
                  "base": drive.BeamModel.__mro__[1].__module__,
                  "rad_rows": sorted(rad),
                  "rad_tau": sum(rad.values()) * c.tau_B,
                  "k26x": list(rows["K26X"]), "k26": list(chain["K26"]),
                  "k108c": list(chain["K108C"]), "r11p_pow": eedf.r11p_pow(),
                  "k15": pure.kr_rates()[1], "tail_refused": refused,
                  "kr2": kr2, "chain_k_kr2_f": formation.LIGHT_CHAIN.k_kr2_f}))
'''


def _probe(version=None):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    env.pop("LIGHT_VERSION", None)
    if version is not None:
        env["LIGHT_VERSION"] = version
    return subprocess.run([sys.executable, "-c", PROBE, str(ML_DIR)],
                          env=env, capture_output=True, text=True)


def _result(version=None):
    r = _probe(version)
    assert r.returncode == 0, r.stderr
    return json.loads(r.stdout.strip().splitlines()[-1])


def test_ml_runs_v0_by_default():
    """With no LIGHT_VERSION, ML builds exactly the Tahoma model it always
    built: v0's basic network with its ion, F2-electron and Kr2* groups,
    B's radiative rate all in G4, v0's ion-ion form and Te^-1/2 template."""
    d = _result()
    assert d["version"] == "v0"
    assert d["groups"] == ["ions", "f2e", "kr2"]
    assert d["base"] == "v0.model"
    assert d["rad_rows"] == ["G4"] and d["rad_tau"] == 1.0
    assert d["k26x"] == ["flannery", 2.4, 3.417e-06]
    assert d["k26"] == [2.4, 3.417e-06] and d["k108c"] == [0.0, 3.096e-06]
    assert d["r11p_pow"] == -0.5
    assert d["k15"] == 2.5e-31
    assert d["kr2"] == {"AK8_Kr": 1.7e-32, "K2F_F": 3.0e-10,
                        "PEN_KR2": 5.0e-10}
    assert d["chain_k_kr2_f"] == 3.0e-10


def test_ml_runs_v0_1_when_named():
    """LIGHT_VERSION=v0_1: its ion and F2-electron rows are basic rows, B's
    radiative rate is shared between G4 and G4A and still sums to 1/tau_B,
    the ion-ion shares keep v0_1's Tg law, the loop uses its Te^-0.55 law,
    and the v0-made tail table cannot be overwritten."""
    d = _result("v0_1")
    assert d["version"] == "v0_1"
    assert d["groups"] == ["kr2"]
    assert d["base"] == "v0_1.model"
    assert d["rad_rows"] == ["G4", "G4A"]
    assert abs(d["rad_tau"] - 1.0) < 1e-12
    assert d["k26x"] == ["flannery", 2.4, 3.417e-06, -3.0]
    assert d["k26"] == [2.4, 3.417e-06, -3.0]
    assert d["k108c"] == [0.0, 3.096e-06, -3.0]
    assert d["r11p_pow"] == -0.55
    assert d["k15"] == 2.4e-31
    assert d["tail_refused"]
    assert d["kr2"] == {"AK8_Kr": 1.7e-32, "K2F_F": 3.0e-10,
                        "PEN_KR2": 5.0e-10}
    assert d["chain_k_kr2_f"] == 3.0e-10


def test_ml_runs_v0_2_when_named():
    """LIGHT_VERSION=v0_2: the Kr2* route is basic, so ML asks for no group
    at all. Its chain keeps the measured links and takes v0_2's 0 for the
    unmeasured Kr2* + F, and Kr2* + Kr2* stays at v0_2's 0. The rest is
    v0_1's."""
    d = _result("v0_2")
    assert d["version"] == "v0_2"
    assert d["groups"] == []
    assert d["base"] == "v0_2.model"
    assert d["rad_rows"] == ["G4", "G4A"]
    assert abs(d["rad_tau"] - 1.0) < 1e-12
    assert d["k26"] == [2.4, 3.417e-06, -3.0]
    assert d["r11p_pow"] == -0.55
    assert d["k15"] == 2.4e-31
    assert d["tail_refused"]
    assert d["kr2"] == {"AK8_Kr": 1.7e-32, "K2F_F": 0.0, "PEN_KR2": 0.0}
    assert d["chain_k_kr2_f"] == 0.0


def test_ml_refuses_a_version_that_does_not_exist():
    r = _probe("v9")
    assert r.returncode != 0
    assert "no LIGHT/run_v9.py" in r.stderr
