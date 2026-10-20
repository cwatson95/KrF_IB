"""ML's krypton swarm analysis (ml/kr_*.py): the registry, the factors on
LIGHT_Kr's cross sections, the error model, the consistency check and the
trained file's free factors. The emulator and the trained file's fit are
checked in test_kr_slow.py."""
import math

import numpy as np
import pytest

from ml import kr_data, kr_light, kr_swarm


def test_registry_files_are_pinned_and_every_set_has_a_doi():
    for key in kr_data.LAPLACE:
        kr_data._pinned(key)                    # raises on a moved md5
    ds = kr_data.load()
    names = [d.name for d in ds]
    assert len(names) == len(set(names))
    for d in ds:
        assert "doi:10." in d.source, d.name
        assert d.role in ("train", "held_out", "check"), d.name
        assert len(d.x) == len(d.y) > 0 and np.all(d.y > 0), d.name
        assert d.kind in kr_swarm.KINDS
    # every point above the forward model's reach is a check point
    for d in ds:
        if d.x.max() > kr_data.ALPHA_TOP_TD:
            assert d.role == "check", d.name
    assert {"a_specht80_total", "a_kruithof40_high"} <= {
        d.name for d in ds if d.role == "check"}


def test_lxcat_units_are_converted():
    """mobility x N (1/(m V s)) times E/N (V m^2) is W (m/s), shown in cm/s;
    LXCat's alpha/N in m^2 is shown in cm^2."""
    hunter = kr_data.by_name("w_hunter88")[0]
    b = [b for b in kr_data.lxcat_blocks(kr_data._pinned("muN"))
         if b["comment"].startswith("Hunter")][0]
    assert hunter.y[0] == pytest.approx(b["y"][0] * b["x"][0] * 1e-21 * 100,
                                        rel=1e-12, abs=0)
    j = kr_data.by_name("a_jacques86")[0]
    raw = [b for b in kr_data.lxcat_blocks(kr_data._pinned("alphaN"))
           if b["comment"].startswith("Jacques")][0]
    k = list(raw["x"]).index(j.x[0])
    assert j.y[0] == pytest.approx(raw["y"][k] * 1e4, rel=1e-12, abs=0)


def test_factors_map_z_onto_lighkr_cross_sections():
    xs = kr_light.module("xsec")
    f0 = kr_swarm.factors(np.zeros(kr_swarm.NP))
    assert f0.m == (1.0,) * 6 and f0.exc == (1.0, 1.0, 1.0) and f0.ion == 1.0
    z = np.zeros(kr_swarm.NP)
    z[kr_swarm.PNAMES.index("x_2p")] = math.log(0.8)
    z[kr_swarm.PNAMES.index("m_0.6")] = math.log(1.2)
    f = kr_swarm.factors(z)
    assert f.exc == pytest.approx((1.0, 0.8, 1.0), rel=1e-12, abs=0)
    assert f.m_at(np.array([0.6]))[0] == pytest.approx(1.2, rel=1e-12, abs=0)
    with xs.using(f):
        p = [q for q in xs.processes("Kr") if q.product == "Kr(2P10)"][0]
    p0 = [q for q in xs.raw_processes("Kr") if q.product == "Kr(2P10)"][0]
    assert np.allclose(p.sigma, 0.8 * p0.sigma, rtol=1e-14, atol=0)


def test_unit_factors_reproduce_the_untrained_swarm():
    xs = kr_light.module("xsec")
    sw = kr_light.module("swarm")
    with xs.using(xs.UNIT):
        a = sw.point(30.0)
    w, al, dt = kr_swarm._solve((tuple(np.zeros(kr_swarm.NP)), 30.0, 300.0))
    assert (w, al) == (a["W"], a["alphaN"])


def test_effective_sigma_adds_the_en_term_through_the_slope():
    d = kr_data.KrSwarm("t", "alphaN", 300.0, np.array([30.0]),
                        np.array([1e-19]), 0.10, "Jacques 1986", "")
    s = kr_swarm.effective_sigma(d, np.array([8.0]))
    assert s[0] == pytest.approx(math.hypot(0.10, kr_swarm.EN_SIGMA * 8.0),
                                 rel=1e-12, abs=0)


def _toy(n_sets=5, bad=None, seed=0):
    """Datasets measured through a linear-in-z toy model; the set `bad` is
    off by x1.6 from every other."""
    rng = np.random.default_rng(seed)
    z_true = 0.1 * rng.standard_normal(kr_swarm.NP) * np.asarray(kr_swarm.FREE)
    sens, ds = {}, []
    for k in range(n_sets):
        x = np.geomspace(1.0, 10.0, 12)
        a = rng.standard_normal((len(x), kr_swarm.NP))
        y0 = np.exp(rng.standard_normal(len(x)))
        y = y0 * np.exp(a @ z_true) * np.exp(0.02 * rng.standard_normal(len(x)))
        if k == bad:
            y = y * 1.6
        d = kr_data.KrSwarm(f"set{k}", "W", 300.0, x, y, 0.03, "Hunter 1988",
                            "", role="train")
        sens[d.name] = (a, y0)
        ds.append(d)

    def model(z, d):
        a, y0 = sens[d.name]
        return y0 * np.exp(a @ np.asarray(z))
    sig = {d.name: np.full(len(d.x), d.sigma) for d in ds}
    return ds, model, sig, z_true


def test_fit_recovers_the_toy_truth_and_refuses_check_data():
    ds, model, sig, z_true = _toy()
    z, cov = kr_swarm.fit(ds, model, sig)
    f = np.asarray(kr_swarm.FREE)
    assert np.max(np.abs(z - z_true)[f] / np.sqrt(np.diag(cov))[f]) < 4.0
    # the factors not trained stay Biagi's: exactly 0 in z, no variance
    assert np.all(z[~f] == 0.0) and np.all(np.diag(cov)[~f] == 0.0)
    chain, acc = kr_swarm.mcmc(ds, model, sig, z, n_walkers=20, n_steps=60,
                               burn=20)
    assert np.all(chain[:, ~f] == 0.0) and 0.05 < acc < 0.95
    ds[0].role = "check"
    with pytest.raises(ValueError):
        kr_swarm.fit(ds, model, sig)


def test_consistency_drops_the_set_that_disagrees():
    ds, model, sig, _ = _toy(bad=2)
    c = kr_swarm.consistency(ds, model, sig, threshold=3.0, log=None
                             or (lambda *a: None))
    assert [n for n, _ in c["dropped"]] == ["set2"]
    assert c["dropped"][0][1] > 50.0
    ds, model, sig, _ = _toy(bad=None)
    c = kr_swarm.consistency(ds, model, sig, log=lambda *a: None)
    assert c["dropped"] == []


def test_the_pinned_record_trains_the_free_factors_only():
    """LIGHT_Kr's factors file was made with kr_swarm.FREE: the factors it
    holds fixed are exactly 1 with no spread, and its variant comparison
    (the evidence for FREE) holds FREE's own variant."""
    import json
    xs = kr_light.module("xsec")
    rec = json.loads((xs.DATA_DIR / xs.TRAINED_FILE[0]).read_text())
    assert rec["free"] == list(kr_swarm.FREE)
    assert rec["names"] == list(kr_swarm.PNAMES)
    f = xs.trained()
    full = list(f.m) + list(f.exc)
    for i, (name, free) in enumerate(zip(kr_swarm.PNAMES, kr_swarm.FREE)):
        if free:
            continue
        assert rec["z"][i] == 0.0 and rec["z_sd"][i] == 0.0, name
        assert full[i] == 1.0, name
        assert rec["posterior"][name] == [1.0, 1.0, 1.0], name
        assert all(d[i] == 0.0 for d in rec["mcmc"]["draws"]), name
    masks = [v["free"] for v in rec["free_variants"].values()]
    assert list(kr_swarm.FREE) in masks
