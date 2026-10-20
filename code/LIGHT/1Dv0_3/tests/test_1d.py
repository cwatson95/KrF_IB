"""The 1-D model on v0_2's chemistry: the bridge to Fluid1D, the stage, the
bands, the 0-D limit, and 1Dv0_3's graded mesh against its uniform one and
against 1Dv0_2.

Fluid1D needs Python >= 3.11 with numba, so run these with such a Python:

    ~/.venvs/f2sweep/bin/python -m pytest LIGHT/1Dv0_3/tests -q

Under another Python (LIGHT's usual /usr/bin/python3) they skip, and the
0-D tests beside them still run. LIGHT_1D_SLOW=1 adds v0_2's whole
reference point in 1-D, on 32 uniform cells (about 50 s) and on the default
graded mesh (about 5 min).
"""

import hashlib
import importlib
import inspect
import math
import os
import sys

import numpy as np
import pytest

if sys.version_info < (3, 11):
    pytest.skip("Fluid1D needs Python >= 3.11", allow_module_level=True)
pytest.importorskip("numba")


def _m(name):
    """1Dv0_3's modules ('1Dv0_3' is not an identifier: no 'from 1Dv0_3 import')."""
    return importlib.import_module(f"1Dv0_3.{name}")


bridge = _m("fluid1d")
try:
    bridge.src_dir()
except ImportError as exc:                       # no Fluid1D beside LIGHT
    pytest.skip(str(exc), allow_module_level=True)

model = _m("model")
xsec = _m("xsec")
sim1d = _m("sim1d")
species = _m("species")


def approx(x, rel):
    # abs=0: pytest.approx's default abs of 1e-12 passes any two rates
    return pytest.approx(x, rel=rel, abs=0)


# every band row at a value: the three unmeasured ones default to 0
BANDS_SET = (("kr2_f", 3e-10), ("kr2_pen", 5e-10), ("kr2f_f", 1e-10))


@pytest.fixture(scope="module", params=[(), BANDS_SET],
                ids=["band defaults", "bands set"])
def chem(request):
    return bridge.CellChemistry(model.Model(model.Params(rates=request.param)))


def _random_state(chem, cells=40, seed=1):
    rng = np.random.default_rng(seed)
    n_gas = 2.45e25
    dens = np.zeros((cells, chem.n_real))
    dens[:, chem.col["Kr"]] = n_gas * 0.995
    dens[:, chem.col["F2"]] = n_gas * 10 ** rng.uniform(-5, -2.3, cells)
    for s in chem.heavy:
        if s not in ("Kr", "F2"):
            dens[:, chem.col[s]] = 10 ** rng.uniform(15, 22, cells)
    dens[:, chem.e_col] = 10 ** rng.uniform(15, 22, cells)
    tg = rng.uniform(250.0, 1500.0, cells)
    eps = 10 ** rng.uniform(math.log10(bridge.EPS_LO), math.log10(50.0),
                            cells)
    return dens, tg, eps


# ---------------------------------------------------------------------------
# The bridge
# ---------------------------------------------------------------------------
def test_registry_matches_species_py(chem):
    bridge.load()
    from particles.constants import molecular_mass
    reg = bridge.registry(chem.heavy)
    assert reg.names == chem.names
    assert reg.get("e-").charge_number == -1
    for s in chem.heavy:
        sp, rec = species.SPECIES[s], reg.get(s)
        assert rec.charge_number == sp.charge
        atoms = sorted(a for a, k in sp.atoms for _ in range(k))
        assert sorted(rec.constituents) == atoms
        assert rec.mass == approx(molecular_mass(*atoms), 1e-15)


# the file holds the electron-impact rows alone, so Fluid1D sees Kr+ and F-
# made and never consumed; the stage consumes them
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cross_sections_round_trip_through_fluid1d(chem, tmp_path):
    """Every process of the two pinned LXCat files, as Fluid1D reads it
    back, is xsec.py's table; and every electron-impact row of the
    reaction file imports cleanly."""
    bridge.load()
    from reactions import import_cross_sections, load_reaction_network
    from reactions.cross_sections import load_cross_section
    xs_dir, toml = bridge.write_inputs(tmp_path, chem.model)
    n = 0
    for gas in xsec.FILES:
        for p in xsec.processes(gas):
            rec = load_cross_section(xs_dir / f"{bridge.dataset_id(gas, p)}"
                                     ".json")
            np.testing.assert_array_equal(rec.energy_ev, p.energy)
            np.testing.assert_allclose(rec.sigma_m2, p.sigma * 1e-4,
                                       rtol=1e-15, atol=0)
            assert rec.target == gas
            kind = "elastic" if p.kind == "effective" else p.kind
            assert rec.reaction_kind.value == kind
            n += 1
    assert n == 17
    compiled = load_reaction_network(
        toml, bridge.registry(chem.heavy)).compile()
    report = import_cross_sections(compiled, xs_dir)
    ids = sorted(i.reaction_id for i in report.imports)
    assert report.n_failed == 0 and len(ids) == 12, report.lines()
    assert ids == sorted(c.row.rid for c in chem.rows
                         if c.row.rate[0] == "xsec")


def test_rates_are_model_k(chem):
    """Per cell, every row's rate is Model._k's at the cell's Te, Tg and
    densities: the Tg, density and ion rows exactly, the Te tables to
    rounding (v0_2's own nodes), the two-temperature tail to 2e-5 (its
    ln(Z_b + exp(-E_b/Te) J) is bilinear between the nodes; 7.8e-6 is the
    largest of 2e5 random points). Below v0_2's floor of 1e-300 cm^3/s
    (Kr excitation at Te = 0.02 eV is ~e^-500) both are zero."""
    dens, tg, eps = _random_state(chem)
    k_static = chem.static_rates(tg, dens)
    log_r = chem.tail_log_ratio(dens)
    k_tab = chem.table_rates(eps, log_r)
    m = chem.model
    for i in range(dens.shape[0]):
        n0 = dens[i, 1:] * 1e-6
        for r, c in enumerate(chem.rows):
            ref = m._k(c, eps[i] / 1.5, tg[i], n0, dens[i, 0] * 1e-6) \
                * chem.si[r]
            got = k_static[i, r] if chem.kind[r] == 0 else k_tab[i, r]
            if ref <= xsec._K_FLOOR * chem.si[r] * (1 + 1e-12):
                assert got <= ref, (c.row.rid, i)
                continue
            rel = {0: 1e-13, 1: 1e-12, 2: 2e-5}[int(chem.kind[r])]
            assert got == approx(ref, rel), (c.row.rid, i)
    kern = bridge.kernel_rates(chem, eps, k_static, log_r)
    np.testing.assert_allclose(kern, k_static + k_tab, rtol=1e-13, atol=0)


def test_param_rows_carry_their_band(chem):
    """v0_2's eight band rows ('param') run at the value Model._compile
    resolved, frac x (Params.rates or the band's default), times their Tg
    law, the same in every cell and constant over the run."""
    bands = _m("network").BANDS
    dens, tg, _ = _random_state(chem, seed=6)
    k = chem.static_rates(tg, dens)
    rates = dict(chem.model.p.rates)
    n = 0
    for r, c in enumerate(chem.rows):
        if c.kind != "param":
            continue
        n += 1
        frac = float(c.row.rate[2]) if len(c.row.rate) > 2 else 1.0
        want = frac * rates.get(c.row.rate[1], bands[c.row.rate[1]].default)
        law = (tg / 300.0) ** c.row.tg_pow if c.row.tg_pow else 1.0
        np.testing.assert_allclose(k[:, r], want * law * chem.si[r],
                                   rtol=1e-13, atol=0)
    assert n == 8


def test_vectorized_helpers_equal_xsec(chem):
    rng = np.random.default_rng(2)
    te = 10 ** rng.uniform(-2.2, 1.6, 300)
    tt = te * 10 ** rng.uniform(-3, 0, 300)
    for r in chem.tail_rows:
        tab = chem.rows[r].rate
        vec = bridge.two_t_rates(tab, te, tt)
        ref = [tab(a, b) for a, b in zip(te, tt)]
        np.testing.assert_allclose(vec, ref, rtol=1e-12, atol=0)
    tail = chem.model.tail
    x = 10 ** rng.uniform(-9, -2, 300)
    f = rng.uniform(0, 0.04, 300)
    np.testing.assert_allclose(bridge.tail_ratio(tail, x, f),
                               [tail(a, b) for a, b in zip(x, f)],
                               rtol=1e-14, atol=0)
    for c in chem.rows:
        if c.route["capture"]:
            ref = [1.5 * t if c.cap is None else c.cap(t) for t in te]
            np.testing.assert_allclose(bridge._capture_energy(c, te), ref,
                                       rtol=1e-13, atol=0)


def test_elastic_power_is_model_rhs(chem):
    """The stage's elastic term is Model.rhs's p_el."""
    dens, tg, eps = _random_state(chem, seed=3)
    got = chem.elastic_power(dens, eps / 1.5, tg)
    m = chem.model
    for i in range(dens.shape[0]):
        te = eps[i] / 1.5
        ref = (dens[i, 0] * 1e-6 * dens[i, chem.col["Kr"]] * 1e-6
               * float(m.kel(te)) * (1 - tg[i] * model.KB_EV / te)) * 1e6
        assert got[i] == approx(ref, 1e-13)


def test_stage_conserves_energy_atoms_and_charge(chem):
    """One stage from random states: per cell, formation energy + electron
    energy + heat + radiation change by exactly the source plus what the
    bounds added; Kr, F and charge are conserved."""
    dens, tg, eps = _random_state(chem, cells=24, seed=4)
    energy = eps * dens[:, 0]
    source = np.random.default_rng(5).uniform(-1e28, 5e30, 24)
    h = 25e-12
    st = chem.stage(dens, energy, eps, tg, source, h, bridge.EPS_LO,
                    bridge.EPS_HI)
    rad = sum(st.counters[f"rad:{b}"] for b in chem.bands)
    lhs = (chem.chemical_energy(st.dens) - chem.chemical_energy(dens)
           + st.energy - energy + st.heat + rad)
    rhs = source * h + st.adjust
    scale = np.abs(source * h) + np.abs(st.heat) + 1e-3 * energy
    assert np.max(np.abs(lhs - rhs) / scale) < 1e-9
    for atom in ("Kr", "F"):
        w = np.array([0.0] + [dict(species.SPECIES[s].atoms).get(atom, 0)
                              for s in chem.heavy])
        np.testing.assert_allclose(st.dens @ w, dens @ w, rtol=1e-12)
    q = np.array([-1.0] + [species.SPECIES[s].charge for s in chem.heavy])
    assert np.max(np.abs(st.dens @ q - dens @ q)) < 1e-9 * np.max(dens[:, 0])
    assert np.all(st.dens >= 0.0)


FLUID1D_HASHES = {
    "Simulation._build": "94d181669d207da9b4ae2b11e548cfeb",
    "Simulation._build_fluids": "cec632b2567b493c16dfc762483c4f4e",
    "Simulation._build_chemistry": "303434fece859c091190b97ffaf5d86e",
    "Simulation._build_plasma_transport": "033f83b538ac518d33e770567aee3928",
    "Simulation._source_exchange": "0a0b8eca1a21d3ce0e94b07405e39773",
    "Simulation._deliverable_heating_j_m3": "54215f93e77df58c3d2df8a7689ede1d",
    "Simulation._advance_electron_source": "664b6873376e971f879954a8e325ac29",
    "Simulation._advance_source_coupled_interval":
        "5025ac606634829b40ed7712afa82d4b",
    "Simulation._step": "a2c7d499bec2436da590409fbe39f3c1",
    "Simulation._exchange": "1d6f3d6462c2dcade7dfab938b5c9769",
    "Simulation._advance_fluid_and_photons":
        "2611a73f13f708db6e455c114e07b94d",
    "Simulation._adopt_interval_state": "3ec003509655f91ee2819f7d85282616",
    "Simulation._advance_optics": "0680a1f45b1393e5c7b9c0c22bfdfde9",
    "Simulation._book_transport_closure": "baa1632b290470a00bf42f04e6fbee6b",
    "Simulation._combined_densities": "4b7febf3b99991a5f2230926d7234ac9",
    "Simulation._adopt_combined_densities":
        "5c1d677890f3aed7294d62e9027c1e29",
    "Simulation.run": "294bbc22bce290419789f5bfe55c27ad",
    "DriftDiffusionPlasmaModel.adopt_electron_energy_density":
        "57832e70ba02ac966223ef44a3ca2aa6",
    "DriftDiffusionPlasmaModel._resolve": "330432c94a424713f3c8d77da9f77fa2",
    "DriftDiffusionPlasmaModel._advance_electron_energy":
        "87a68018e9544bc38494fe57e4d69b30",
    "integrator._augmented_linearly_implicit_trial":
        "cb98c4f6ce32b364e5ac789ecd8b0b8f",
    "integrator._energy_coupled_substeps": "d3f744514812403bcffae9530242277f",
    # 1Dv0_3's mesh: the orchestrator's books it redoes per cell, and the
    # kernel transport.py follows (and the functions of it that it calls)
    "Simulation._initialize_runtime_state": "0fc4238a2cab128fd8eb1804a33da6d6",
    "Simulation._initialize_energetics": "8bba6270417848dc10e0c63765f59b01",
    "Simulation._resync_transport_energy_after_optics":
        "8e38c5a2a4ad419a22e06e01f3650777",
    "DriftDiffusionPlasmaModel.advance": "18e092a2f2a78c7dbefb56aac1d7a717",
    "DriftDiffusionPlasmaModel._advance_via_kernel":
        "386fd76b9f934d8bfa1e6e2bff51c7af",
    "transport_kernel._coupled_loop": "c497232c0a60d5deb0faa61b214c5a58",
    "transport_kernel._cfl_timestep": "1532e3dd1f8bc2d70a8c06ce88448e75",
    "transport_kernel._solve_field": "31ec298350643f2b12c24127e5899a69",
    "transport_kernel._transport_species": "17ac2c7326545370bcaa6a5f85ba2e65",
    "transport_kernel._transport_energy": "b9704204c2bf756734b135e83722b695",
    "transport_kernel.coupled_advance": "cb5a2cefd49910caa9a60c986aaa6937",
    "transport_kernel.pack_kernel_args": "33ff95ff4e7029755c547fafdaa6f810",
    "transport_kernel._recover_energy": "dc9c535d2ba56704257b7d73df20f2cc",
    "transport_kernel._lookup_transport": "a169feac7910b962e6f2ad5f55250902",
    "transport_kernel._advance_gap": "139bf663692e6958cea5f4ad6397c0f7",
    "transport_kernel._inelastic_power_cell":
        "051934d4987834963f0b7d1558c40b3e",
    "transport_kernel._bernoulli": "b1e05b8ab9ab3aece4d423dad36e77de",
    "Poisson1DSolver.solve_semi_implicit": "6e4daa6d5a5b2127349f2a2ff905041f",
    "Poisson1DSolver.solve": "15d5900d89275609cfad04a74941323f",
    "SemiImplicitField.solve": "fae58aca51cf4153cea4dcb9305c423d",
    "field._thomas_general": "74a023d946ac785484c0418c6c209312",
}


def test_fluid1d_methods_1dv0_relies_on_are_unchanged():
    """sim1d.LightSimulation overrides or calls these Fluid1D methods, the
    stage follows _energy_coupled_substeps, and transport.py follows
    Fluid1D's kernel and calls its geometry-free functions. A failure means
    Nathaniel's code moved: re-read the method against sim1d.py,
    fluid1d.py or transport.py (the split's order, the ledger names, the
    trial's signature, the kernel's arithmetic), then record its new hash
    here."""
    bridge.load()
    import Fluid1D.field as ff
    import Fluid1D.plasma_model as pm
    import Fluid1D.simulation as fs
    import Fluid1D.transport_kernel as tk
    import kinetics.integrator as ki
    owners = {"Simulation": fs.Simulation,
              "DriftDiffusionPlasmaModel": pm.DriftDiffusionPlasmaModel,
              "integrator": ki, "transport_kernel": tk,
              "Poisson1DSolver": ff.Poisson1DSolver,
              "SemiImplicitField": ff.SemiImplicitField, "field": ff}
    moved = []
    for key, want in FLUID1D_HASHES.items():
        owner, name = key.split(".")
        f = getattr(owners[owner], name)
        f = getattr(f, "py_func", f)
        got = hashlib.md5(inspect.getsource(f).encode()).hexdigest()
        if got != want:
            moved.append(key)
    assert not moved, f"Fluid1D changed: {moved}"


# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------
# the gas held still and Fluid1D's field work (1Dv0_2's model) unless a
# run says otherwise
SHORT = dict(cells=8, mesh="uniform", neutrals="still", field_work="fluid1d",
             pulse_s=20e-9, t_end_s=40e-9, dt_s=25e-12, n_save=201,
             n_profiles=5, verbose=False)


@pytest.fixture(scope="module")
def short_runs(tmp_path_factory):
    out = str(tmp_path_factory.mktemp("runs1d"))
    p = sim1d.Params1D(**SHORT, output_dir=out, name="short8")
    band = p.replace(rates=(("kr2_f", 1e-9),), name="short8_kr2f")
    res = {"8": sim1d.run(p),
           "16": sim1d.run(p.replace(cells=16, name="short16")),
           "dt50": sim1d.run(p.replace(dt_s=50e-12, name="short8_dt50")),
           "band": sim1d.run(band),
           "u48": sim1d.run(p.replace(cells=48, name="short48")),
           "g48": sim1d.run(p.replace(cells=48, mesh="graded",
                                      wall_cell_um=10.0, name="short48g")),
           "g48m": sim1d.run(p.replace(cells=48, mesh="graded",
                                       wall_cell_um=10.0, neutrals="move",
                                       name="short48gm")),
           "u8m": sim1d.run(p.replace(neutrals="move", name="short8m")),
           "g48c": sim1d.run(p.replace(cells=48, mesh="graded",
                                       wall_cell_um=10.0,
                                       field_work="conserve",
                                       name="short48gc")),
           "g48mc": sim1d.run(p.replace(cells=48, mesh="graded",
                                        wall_cell_um=10.0, neutrals="move",
                                        field_work="conserve",
                                        name="short48gmc"))}
    res["0d"] = model.run(p.cell_params().replace(n_save=4001))
    res["band0d"] = model.run(band.cell_params().replace(n_save=4001))
    return res


def test_centre_cell_is_the_0d_model(short_runs):
    """With uniform power the gap's interior is v0_2's 0-D run: the centre
    cell starts in v0_2's state, and its emission by band, KrF formed,
    species, Tg and (after the first 2 ns) Te follow v0_2's."""
    res, r0 = short_runs["8"], short_runs["0d"]
    report = _m("report")
    first = res.data["records"][0]["centre"]
    for i, s in enumerate(res.species):
        v0 = r0.ne[0] if s == "e-" else r0.n[s][0]
        assert first[i] == approx(v0, 1e-12), s
    d0 = report.row_dict(r0)
    assert res.share_248("ledger_centre") == approx(d0["emission_248"], 1e-3)
    acc = report.accounting(r0)
    led = res.data["records"][-1]["ledger_centre"]
    for b in res.data["bands"]:
        assert led[f"rad:{b}"] / led["input"] == approx(
            acc["bands"][b] / acc["e_dep"], 1e-3), b
    assert led["krf_formed"] == approx(
        report.reservoir(r0, "KrF*")["formed"], 1e-3)
    last = res.data["records"][-1]
    assert last["tg_centre"] == pytest.approx(d0["tg_end_K"], abs=0.5)
    for s in res.species:
        v0 = r0.ne[-1] if s == "e-" else r0.n[s][-1]
        assert res.series(s, "centre")[-1] == approx(v0, 1e-2), s
    te1 = np.array([r["te_centre"] for r in res.data["records"]])
    te0 = np.interp(res.t, r0.t, r0.te)
    late = res.t > 2e-9
    assert np.max(np.abs(te1[late] / te0[late] - 1)) < 0.01


def test_a_band_reaches_the_1d_chemistry(short_runs):
    """Kr2* + F at the top of its band (1e-9): the centre cell follows v0_2's
    0-D run with the same band, and that run differs from the default by
    far more than the 1-D and 0-D runs differ from each other."""
    report = _m("report")
    res, r0, d0 = short_runs["band"], short_runs["band0d"], short_runs["0d"]
    band, dflt = report.row_dict(r0), report.row_dict(d0)
    got = res.share_248("ledger_centre")
    assert got == approx(band["emission_248"], 1e-3)
    led = res.data["records"][-1]["ledger_centre"]
    vuv = 100 * led["rad:147 nm (Kr2* VUV)"] / led["input"]
    assert vuv == approx(band["emission_vuv"], 1e-3)
    assert res.series("Kr2*", "centre")[-1] == approx(r0.n["Kr2*"][-1], 1e-2)
    assert abs(band["emission_248"] / dflt["emission_248"] - 1) > 5e-3
    assert abs(band["emission_vuv"] / dflt["emission_vuv"] - 1) > 3e-2
    assert "bands: kr2_f 1e-09" in res.summary()


def test_energy_ledger_closes(short_runs):
    """The gap's chemical + electron + thermal energy, recomputed from its
    state, equals what the ledger booked (input - radiation + bounds +
    transport), and the fate adds up to 100 %."""
    for key in ("8", "16", "dt50", "band", "u48", "g48", "g48m", "u8m",
                "g48c", "g48mc"):
        res = short_runs[key]
        fate = res.fate()
        assert abs(fate["closure residual"]) < 1e-8, key
        assert sum(fate.values()) == pytest.approx(100.0, abs=1e-6), key


def test_interior_is_uniform_and_the_cut_is_in_the_wall_cells(short_runs):
    """Uniform power: the interior is the centre cell (the two cells beside
    the wall cells to 1e-3, the rest to 1e-6), and the energy the
    mean-energy ceiling cuts is all in the two wall cells, whose electrons
    the walls drain."""
    for key in ("8", "16"):
        res = short_runs[key]
        led = res.data["ledger_cells"]
        cut = -np.array(led["adjust"])
        assert np.all(np.abs(cut[1:-1]) < 1e-9 * np.max(cut)), key
        rad = np.array(led["rad:248 nm (KrF B-X)"])
        mid = rad[len(rad) // 2]
        np.testing.assert_allclose(rad[2:-2], mid, rtol=1e-6)
        np.testing.assert_allclose(rad[[1, -2]], mid, rtol=1e-3)
    f8 = short_runs["8"].fate()["cut by the mean-energy bounds"]
    f16 = short_runs["16"].fate()["cut by the mean-energy bounds"]
    assert 0 < f16 < f8


def test_te_converges_first_order_in_dt(short_runs):
    """Halving dt halves the Te error against the 0-D run (the split and
    the linearly implicit stage are first order)."""
    r0 = short_runs["0d"]

    def err(res):
        te1 = np.array([r["te_centre"] for r in res.data["records"]])
        te0 = np.interp(res.t, r0.t, r0.te)
        return np.max(np.abs(te1[1:] / te0[1:] - 1))
    ratio = err(short_runs["dt50"]) / err(short_runs["8"])
    assert 1.5 < ratio < 2.6


def test_result_file_round_trip(short_runs, tmp_path):
    res = short_runs["8"]
    back = sim1d.Result1D.load(res.save(tmp_path / "r.json"))
    assert back.share_248() == res.share_248()
    assert back.data["provenance"]["chemistry_md5"] == sim1d.chemistry_md5()


def test_cli(tmp_path, capsys):
    main = _m("run").main
    assert main(["--cells", "4", "--mesh", "uniform", "--pulse", "0.002",
                 "--t-end", "0.004",
                 "--out", str(tmp_path), "--name", "cli", "--quiet",
                 "--rate", "kr2_f=1e-10", "--compare"]) == 0
    out = capsys.readouterr().out
    assert "radiated 248 nm (KrF B-X)" in out and "0-D (v0_2's model)" in out
    assert "bands: kr2_f 1e-10;  no measurement, 0: kr2_pen, kr2f_f" in out
    assert "147 nm, % of deposited" in out
    assert (tmp_path / "cli" / "result.json").is_file()
    assert (tmp_path / "cli" / "fluid1d_inputs" / "reactions"
            / "light_1Dv0_3.toml").is_file()
    with pytest.raises(SystemExit):
        main(["--rate", "k2f=3e-10"])


def test_input_is_the_power_over_the_gap_s_volume(short_runs):
    """Every total is a sum over the cells' own volumes: on either mesh the
    gap's mean input is the given power over the pulse."""
    for key in ("8", "u48", "g48"):
        led = short_runs[key].data["records"][-1]["ledger"]
        assert led["input"] * model.QE == approx(3.0e5 * SHORT["pulse_s"],
                                                 1e-12), key


def test_means_and_profiles_weigh_the_cells_by_volume(short_runs):
    """On the graded mesh a gap mean is a volume-weighted mean (the final
    F2 against the last profile's), and a power profile has a
    volume-weighted mean of 1, so it moves power without changing the
    total."""
    g = short_runs["g48"]
    w = g.widths_cm / np.sum(g.widths_cm)
    f2 = np.array(g.data["profiles"][-1]["dens"][g.species.index("F2")])
    assert g.data["profiles"][-1]["t"] == g.t[-1]
    assert g.series("F2")[-1] == approx(float(np.sum(f2 * w)), 1e-12)
    assert abs(float(np.mean(f2)) / float(np.sum(f2 * w)) - 1) > 1e-6
    # twice the power in the ten small cells at each wall (a symmetric
    # profile on a symmetric mesh needs weights that differ from 1/N)
    raised = np.ones(48)
    raised[:10] = raised[-10:] = 2.0
    p = sim1d.Params1D(cells=48, wall_cell_um=10.0,
                       power_profile=tuple(raised))
    m = p.make_mesh(1e-2, 1e-4)
    f = p.profile(m.volumes)
    assert float(np.sum(f * m.volumes) / np.sum(m.volumes)) == approx(
        1.0, 1e-14)
    assert abs(float(np.mean(f)) - 1.0) > 1e-3


def test_graded_mesh_centre_is_the_0d_model(short_runs):
    """The graded mesh keeps the interior: its centre cell gives v0_2's 0-D
    emission by band, KrF formed and Tg, as the uniform mesh's does."""
    res, r0 = short_runs["g48"], short_runs["0d"]
    report = _m("report")
    d0, acc = report.row_dict(r0), report.accounting(r0)
    assert res.share_248("ledger_centre") == approx(d0["emission_248"], 1e-3)
    led = res.data["records"][-1]["ledger_centre"]
    for b in res.data["bands"]:
        assert led[f"rad:{b}"] / led["input"] == approx(
            acc["bands"][b] / acc["e_dep"], 1e-3), b
    assert led["krf_formed"] == approx(
        report.reservoir(r0, "KrF*")["formed"], 1e-3)
    assert res.data["records"][-1]["tg_centre"] == pytest.approx(
        d0["tg_end_K"], abs=0.5)


def test_graded_mesh_resolves_the_walls(short_runs):
    """On the same 48 cells, cells of 10 um at the electrodes (against 208
    um) take the sheath: the energy the walls take, the energy the
    mean-energy ceiling cuts in drained cells and the negative field work
    Fluid1D's split keeps out are small parts of the uniform mesh's, and
    the gap's 248-nm share comes near the centre cell's."""
    u, g = short_runs["u48"], short_runs["g48"]
    fu, fg = u.fate(), g.fate()
    assert g.data["mesh"].startswith("48 graded cells, 10 um at the walls")
    assert np.sum(g.widths_cm) == approx(1.0, 1e-13)
    assert fg["lost to the walls"] < 0.2 * fu["lost to the walls"]
    assert fg["cut by the mean-energy bounds"] \
        < 0.1 * fu["cut by the mean-energy bounds"]

    def kept(r):
        led = r.data["records"][-1]["ledger"]
        return -r.data["field_dropped"] / led["input"]
    assert kept(g) < 2e-3 * kept(u)
    gap_u = u.share_248() - u.share_248("ledger_centre")
    gap_g = g.share_248() - g.share_248("ledger_centre")
    assert abs(gap_g) < 0.1 * abs(gap_u)


def test_uniform_mesh_is_1dv0_2_to_the_bit(short_runs, tmp_path):
    """1Dv0_3 on its uniform mesh runs 1Dv0_2's model: the same short run
    gives the same cells to the bit (densities, Te, Tg, every ledger), and
    the gap's totals to their summation order (1Dv0_3 sums over the cells'
    volumes, 1Dv0_2 times one volume)."""
    try:
        s2 = importlib.import_module("1Dv0_2.sim1d")
    except ImportError:
        pytest.skip("1Dv0_2 is not beside 1Dv0_3")
    res = short_runs["8"]
    same = {k: v for k, v in SHORT.items()                 # 1Dv0_2's own
            if k not in ("mesh", "neutrals", "field_work")}
    ref = s2.run(s2.Params1D(**same, output_dir=str(tmp_path),
                             name="short8_1dv0_2"))
    for k, v in ref.data["ledger_cells"].items():
        assert np.array_equal(np.array(v),
                              np.array(res.data["ledger_cells"][k])), k
    for a, b in zip(ref.data["profiles"], res.data["profiles"]):
        for k in ("dens", "te", "tg"):
            assert np.array_equal(np.array(a[k]), np.array(b[k])), (a["t"], k)
    # e_other ("lost to the walls") sums differences of the gap's whole
    # energy, which dwarfs the input, over every exchange: 1e-12 of it.
    # The closure residual is that rounding itself (1e-11 %)
    fa, fb = ref.fate(), res.fate()
    for k in fa:
        assert fb[k] == pytest.approx(fa[k], rel=1e-9, abs=1e-9), k


def test_the_conserving_field_work_keeps_its_books(short_runs):
    """field_work='conserve' on 48 graded cells: the field's work is not an
    input (none given, none kept out), the electrons give the field what
    the ions take from it, and the electrons' work, the ions' and the
    field's energy at the end, less what the gas changed in it carrying the
    charge (none with the gas still), sum to the theorem's residual, which
    is round-off. The ions' share heats the gas (its own line in the
    fate), the bounds book what the floor holds in the drained wall cells,
    the interior keeps the Fluid1D run's centre, and the summary reports
    the field's books instead of a field work given."""
    for key in ("g48c", "g48mc"):
        res = short_runs[key]
        d = res.data
        led = d["records"][-1]["ledger"]
        assert led["field_work"] == 0.0 and d["field_dropped"] == 0.0, key
        assert led["field_e"] < 0.0 < led["field_i"], key
        assert abs(d["field_defect"]) < 1e-9 * abs(led["field_e"]), key
        books = (led["field_e"] + led["field_i"] + d["field_energy"]
                 - d["field_moved"])
        assert abs(books) < 1e-9 * abs(led["field_e"]), key
        f = res.fate()
        assert f["heat (ions in the field)"] == approx(
            100 * led["field_i"] / led["input"], 1e-12), key
        assert d["positivity_limits"] > 0, key
        assert res.share_248("ledger_centre") == approx(
            short_runs["g48"].share_248("ledger_centre"), 1e-4), key
        text = res.summary()
        assert "the field's work, kept to its energy" in text, key
        assert "field work given" not in text, key
    assert short_runs["g48c"].data["field_moved"] == 0.0
    assert short_runs["g48mc"].data["field_moved"] != 0.0
    assert "the gas carrying the charge" in short_runs["g48mc"].summary()
    assert "the gas carrying the charge" not in short_runs["g48c"].summary()
    # the field's work moves energy inside the gap, never out of it, so the
    # conserving run loses no more to the walls than Fluid1D's, whose walls'
    # line also holds the transport step's own bounds (0.066 % against
    # 0.106 % of the input). The ions' share booked but left out of the gas
    # would show there as 0.94 % more
    assert (short_runs["g48c"].fate()["lost to the walls"]
            < short_runs["g48"].fate()["lost to the walls"])
    cells = short_runs["g48c"].data["ledger_cells"]
    held = np.array(cells["adjust"])
    assert np.all(np.abs(held[1:-1]) < 1e-9 * np.max(np.abs(held)))


def test_cli_graded(tmp_path, capsys):
    main = _m("run").main
    assert main(["--cells", "48", "--mesh", "graded", "--wall-cell", "10",
                 "--pulse", "0.002", "--t-end", "0.004", "--out",
                 str(tmp_path), "--name", "clig", "--quiet"]) == 0
    out = capsys.readouterr().out
    assert "48 graded cells, 10 um at the walls" in out
    assert "transport substeps" in out
    assert "gas: |u| up to" in out                   # the default moves it
    assert main(["--cells", "8", "--mesh", "uniform", "--neutrals", "still",
                 "--pulse", "0.002", "--t-end", "0.004", "--out",
                 str(tmp_path), "--name", "clis", "--quiet"]) == 0
    text = capsys.readouterr().out
    assert "gas: held still (neutrals='still')" in text
    assert "the field's work, kept to its energy" in text   # the default
    assert main(["--cells", "8", "--mesh", "uniform", "--neutrals", "still",
                 "--field-work", "fluid1d", "--pulse", "0.002", "--t-end",
                 "0.004", "--out", str(tmp_path), "--name", "clif",
                 "--quiet"]) == 0
    assert "field work given" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["--mesh", "stretched"])
    with pytest.raises(SystemExit):
        main(["--neutrals", "fly"])
    with pytest.raises(SystemExit):
        main(["--field-work", "rectified"])


def test_the_moving_gas_keeps_the_books_and_the_interior(short_runs):
    """neutrals='move' on 48 graded cells: the books close with the heat the
    electrodes take booked, the centre cell is the still gas's (the waves
    from the walls are microns from it), the gas moves near the walls, and
    the profiles carry its velocity and pressure."""
    m, s = short_runs["g48m"], short_runs["g48"]
    fm = m.fate()
    assert abs(fm["closure residual"]) < 1e-8
    assert sum(fm.values()) == pytest.approx(100.0, abs=1e-6)
    assert m.share_248("ledger_centre") == approx(
        s.share_248("ledger_centre"), 1e-9)
    last = m.data["profiles"][-1]
    u, p = np.abs(np.array(last["u"])), np.array(last["p"])
    mid = len(u) // 2
    assert u.max() > 0.1 and u[mid] < 1e-9 * u.max()
    assert np.all(p > 0.9e5) and np.all(p < 1.2e5)
    assert m.data["records"][-1]["ledger"]["wall_heat"] > 0.0
    assert m.data["gas_substeps"]["mean"] >= 1.0
    assert "u" not in s.data["profiles"][-1] and s.data["gas_substeps"] is None
    assert "gas: |u| up to" in m.summary()
    assert "gas: held still" in s.summary()
    # on 8 uniform cells the wall cells are 1250 um: the moving gas stays
    # inside them, and the interior is the still run's
    u8, s8 = short_runs["u8m"], short_runs["8"]
    assert u8.share_248("ledger_centre") == approx(
        s8.share_248("ledger_centre"), 1e-9)
    assert abs(u8.fate()["closure residual"]) < 1e-8


@pytest.mark.skipif(not os.environ.get("LIGHT_1D_SLOW"),
                    reason="set LIGHT_1D_SLOW=1 (about 50 s)")
def test_reference_point_centre_is_v0_2(tmp_path):
    """v0_2's reference point (1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, 2 us)
    on 32 cells: the centre cell's 248-nm share and Tg are the 0-D run's."""
    p = sim1d.Params1D(cells=32, mesh="uniform", verbose=False,
                       output_dir=str(tmp_path))
    res = sim1d.run(p)
    d0 = _m("report").row_dict(model.run(p.cell_params()))
    assert res.share_248("ledger_centre") == approx(d0["emission_248"], 1e-4)
    assert res.data["records"][-1]["tg_centre"] == pytest.approx(
        d0["tg_end_K"], abs=0.1)
    assert abs(res.fate()["closure residual"]) < 1e-8


@pytest.mark.skipif(not os.environ.get("LIGHT_1D_SLOW"),
                    reason="set LIGHT_1D_SLOW=1 (about 13 min)")
def test_reference_point_on_the_graded_mesh(tmp_path):
    """v0_2's reference point on 1Dv0_3's default, the gas moving on 128
    cells graded from 0.95 um at the walls (the 64-cell mesh from 2 um
    refined once): the centre cell is the 0-D run, the books close with the
    heat the electrodes take, the bounds cut what the drained cells within
    about 10 um of the walls take (0.14 % of the input, where 64 uniform
    cells cut 3.1 %), the field's work is kept to the field's energy (none
    of it input, the electrons giving the ions 0.49 % of the given power,
    the books' residual round-off), and the gap's 248-nm share of the given
    power lies below the interior's by the walls' deficit (4.63 % against
    4.67)."""
    p = sim1d.Params1D(verbose=False, output_dir=str(tmp_path))
    assert (p.mesh, p.cells, p.neutrals, p.field_work) == (
        "graded", 128, "move", "conserve")
    m = sim1d.meshes.refine(sim1d.meshes.graded(64, p.gap_cm * 1e-2, 2e-6))
    assert np.array_equal(p.make_mesh(p.gap_cm * 1e-2, p.area_cm2 * 1e-4)
                          .edges, m.edges)
    res = sim1d.run(p)
    d0 = _m("report").row_dict(model.run(p.cell_params()))
    assert res.share_248("ledger_centre") == approx(d0["emission_248"], 1e-4)
    f = res.fate()
    last = res.data["records"][-1]["ledger"]
    assert abs(f["closure residual"]) < 1e-8
    assert 0.0 < f["cut by the mean-energy bounds"] < 0.5
    assert f["lost to the walls"] < 0.05
    assert last["field_work"] == 0.0
    assert last["field_e"] < 0.0 < last["field_i"] < 0.01 * last["input"]
    assert abs(res.data["field_defect"]) < 1e-9 * abs(last["field_e"])
    assert last["wall_heat"] > 0.0
    assert d0["emission_248"] - 0.15 < res.share_248() < d0["emission_248"]
