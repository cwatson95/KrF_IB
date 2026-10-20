"""The 1-D ArF model on LIGHT_ArF_ebeam v0_2's chemistry: the bridge to
Fluid1D, the stage, the bands, the 0-D limit, and both drives (LIGHT's own
input and the electron beam).

Fluid1D needs Python >= 3.11 with numba, so run these with such a Python:

    ~/.venvs/f2sweep/bin/python -m pytest LIGHT_ArF_ebeam/1Dv0_2/tests -q

Under another Python (LIGHT's usual /usr/bin/python3) they skip, and the
0-D tests beside them still run. LIGHT_1D_SLOW=1 adds v0_2's whole
reference point in 1-D (about a minute).
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
    """1Dv0_2's modules ('1Dv0_2' is not an identifier: no 'from 1Dv0_2 import')."""
    return importlib.import_module(f"1Dv0_2.{name}")


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


# every band row at a value: the four unmeasured ones default to 0
BANDS_SET = (("ar2_f", 3e-10), ("ar2_pen", 5e-10), ("ar_pen", 5e-10),
             ("ar2f_ar", 1e-14), ("arf_yield", 0.8))


@pytest.fixture(scope="module", params=[(), BANDS_SET],
                ids=["band defaults", "bands set"])
def chem(request):
    return bridge.CellChemistry(model.Model(model.Params(rates=request.param)))


def _random_state(chem, cells=40, seed=1):
    rng = np.random.default_rng(seed)
    n_gas = 2.45e25
    dens = np.zeros((cells, chem.n_real))
    dens[:, chem.col["Ar"]] = n_gas * 0.995
    dens[:, chem.col["F2"]] = n_gas * 10 ** rng.uniform(-5, -2.3, cells)
    for s in chem.heavy:
        if s not in ("Ar", "F2"):
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


# the file holds the electron-impact rows alone, so Fluid1D sees Ar+ and F-
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
    assert n == 50                  # 39 Ar (IST-Lisbon) + 11 F2 (SIGLO)
    compiled = load_reaction_network(
        toml, bridge.registry(chem.heavy)).compile()
    report = import_cross_sections(compiled, xs_dir)
    ids = sorted(i.reaction_id for i in report.imports)
    assert report.n_failed == 0 and len(ids) == 45, report.lines()
    assert ids == sorted(c.row.rid for c in chem.rows
                         if c.row.rate[0] == "xsec")


def test_rates_are_model_k(chem):
    """Per cell, every row's rate is Model._k's at the cell's Te, Tg and
    densities: the Tg, density and ion rows exactly, the Te tables to
    rounding (v0_2's own nodes), the two-temperature tail to 2e-5 (its
    ln(Z_b + exp(-E_b/Te) J) is bilinear between the nodes; 7.8e-6 is the
    largest of 2e5 random points). Below v0_2's floor of 1e-300 cm^3/s
    (Ar excitation at Te = 0.02 eV is ~e^-580) both are zero."""
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
    """v0_2's thirteen band rows ('param') run at the value Model._compile
    resolved (network.param_rate: the band's value from Params.rates or its
    default, its modifiers), times their Tg law, the same in every cell and
    constant over the run."""
    net = _m("network")
    dens, tg, _ = _random_state(chem, seed=6)
    k = chem.static_rates(tg, dens)
    m = chem.model
    n = 0
    for r, c in enumerate(chem.rows):
        if c.kind != "param":
            continue
        n += 1
        want = net.param_rate(c.row.rate, m.band_value, m.p.tau_sp_s)
        law = (tg / 300.0) ** c.row.tg_pow if c.row.tg_pow else 1.0
        np.testing.assert_allclose(k[:, r], want * law * chem.si[r],
                                   rtol=1e-13, atol=0)
    assert n == 13


def test_trimer_law_reaches_the_cells():
    """A trimer_tg_pow set in the run (a sensitivity run of the
    unmeasured law) acts on G5 in every cell as in Model._k, and on no
    other row."""
    m = model.Model(model.Params(trimer_tg_pow=-3.9))
    ch = bridge.CellChemistry(m)
    base = bridge.CellChemistry(model.Model(model.Params()))
    dens, tg, _ = _random_state(ch, seed=7)
    k, k0 = ch.static_rates(tg, dens), base.static_rates(tg, dens)
    for r, c in enumerate(ch.rows):
        if ch.kind[r]:
            continue
        if c.row.rid == "G5":
            np.testing.assert_allclose(k[:, r], k0[:, r] * (tg / 300.0)
                                       ** -3.9, rtol=1e-13, atol=0)
        else:
            np.testing.assert_allclose(k[:, r], k0[:, r], rtol=1e-13,
                                       atol=0)


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
        ref = (dens[i, 0] * 1e-6 * dens[i, chem.col["Ar"]] * 1e-6
               * float(m.kel(te)) * (1 - tg[i] * model.KB_EV / te)) * 1e6
        assert got[i] == approx(ref, 1e-13)


def test_stage_conserves_energy_atoms_and_charge(chem):
    """One stage from random states: per cell, formation energy + electron
    energy + heat + radiation change by exactly the source plus what the
    bounds added; Ar, F and charge are conserved."""
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
    for atom in ("Ar", "F"):
        w = np.array([0.0] + [dict(species.SPECIES[s].atoms).get(atom, 0)
                              for s in chem.heavy])
        np.testing.assert_allclose(st.dens @ w, dens @ w, rtol=1e-12)
    q = np.array([-1.0] + [species.SPECIES[s].charge for s in chem.heavy])
    assert np.max(np.abs(st.dens @ q - dens @ q)) < 1e-9 * np.max(dens[:, 0])
    assert np.all(st.dens >= 0.0)


def test_stage_conserves_energy_atoms_and_charge_with_the_beam():
    """The beam's sources in one stage: formation energy + electron energy
    + heat + radiation change by the energy source plus the formation
    energy of what the beam made (plus what the bounds added), and the
    beam's products keep Ar, F and charge (an Ar or F2 ionization is one
    Ar+ and one electron, a prompt F2 dissociation two F)."""
    m = model.Model(model.Params(drive="beam"))
    ch = bridge.CellChemistry(m)
    y = m.beam
    assert ch.beam_vec[ch.e_col] == 1.0
    assert ch.beam_vec[ch.col["Ar+"]] == 1.0
    assert ch.beam_vec[ch.col["Ar*"]] == pytest.approx(y.ar_star, rel=1e-12)
    dens, tg, eps = _random_state(ch, cells=24, seed=6)
    energy = eps * dens[:, 0]
    pairs = np.random.default_rng(7).uniform(1e28, 1e31, 24)   # m^-3 s^-1
    source = pairs * y.eps_sec
    made = pairs[:, None] * ch.beam_vec[None, :]
    h = 25e-12
    st = ch.stage(dens, energy, eps, tg, source, h, bridge.EPS_LO,
                  bridge.EPS_HI, species_source=made)
    rad = sum(st.counters[f"rad:{b}"] for b in ch.bands)
    lhs = (ch.chemical_energy(st.dens) - ch.chemical_energy(dens)
           + st.energy - energy + st.heat + rad)
    rhs = (source + made @ ch.e_form) * h + st.adjust
    scale = (np.abs(source * h) + np.abs(made @ ch.e_form * h)
             + np.abs(st.heat) + 1e-3 * energy)
    assert np.max(np.abs(lhs - rhs) / scale) < 1e-9
    for atom in ("Ar", "F"):
        w = np.array([0.0] + [dict(species.SPECIES[s].atoms).get(atom, 0)
                              for s in ch.heavy])
        np.testing.assert_allclose(st.dens @ w, dens @ w, rtol=1e-12)
    q = np.array([-1.0] + [species.SPECIES[s].charge for s in ch.heavy])
    assert np.max(np.abs(st.dens @ q - dens @ q)) < 1e-9 * np.max(dens[:, 0])
    assert np.all(st.dens >= 0.0)
    # with no source the stage is the own input's, to the last bit
    a = ch.stage(dens, energy, eps, tg, source, h, bridge.EPS_LO,
                 bridge.EPS_HI)
    b = ch.stage(dens, energy, eps, tg, source, h, bridge.EPS_LO,
                 bridge.EPS_HI, species_source=np.zeros_like(made))
    assert np.array_equal(a.dens, b.dens) and np.array_equal(a.energy, b.energy)


def test_beam_electrons_enter_inside_the_energy_bounds():
    """A cell at the mean-energy ceiling (as a wall cell the walls drain
    sits) takes the beam's new electrons with their eps_sec each, far under
    the ceiling: the bounds are read on the electrons after the source, so
    nothing is cut. Read on the electrons before it, the ceiling of the old
    population would cut the new electrons' energy."""
    m = model.Model(model.Params(drive="beam"))
    ch = bridge.CellChemistry(m)
    dens = np.zeros((1, ch.n_real))
    n_gas = 2.45e25
    dens[0, ch.col["Ar"]] = n_gas * 0.995
    dens[0, ch.col["F2"]] = n_gas * 0.005
    dens[0, ch.col["Ar+"]] = dens[0, ch.e_col] = 1e18
    eps = np.array([bridge.EPS_HI])
    energy = eps * dens[:, 0]
    pairs = np.array([1e30])
    st = ch.stage(dens, energy, eps, np.array([300.0]),
                  pairs * m.beam.eps_sec, 25e-12, bridge.EPS_LO,
                  bridge.EPS_HI,
                  species_source=pairs[:, None] * ch.beam_vec[None, :])
    assert st.dens[0, ch.e_col] > 10 * dens[0, ch.e_col]
    assert abs(st.adjust[0]) < 1e-9 * energy[0]


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
}


def test_fluid1d_methods_1dv0_relies_on_are_unchanged():
    """sim1d.LightSimulation overrides or calls these Fluid1D methods, and
    the stage follows _energy_coupled_substeps. A failure means Nathaniel's
    code moved: re-read the method against sim1d.py/fluid1d.py (the
    split's order, the ledger names, the trial's signature), then record
    its new hash here."""
    bridge.load()
    import Fluid1D.plasma_model as pm
    import Fluid1D.simulation as fs
    import kinetics.integrator as ki
    owners = {"Simulation": fs.Simulation,
              "DriftDiffusionPlasmaModel": pm.DriftDiffusionPlasmaModel,
              "integrator": ki}
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
# how far the two bands' tops move the short run's emission (measured on
# the 0-D runs, far above the 1e-3 the 1-D and 0-D runs agree to)
BAND_193, BAND_VUV = 0.1, 1.0


SHORT = dict(cells=8, pulse_s=20e-9, t_end_s=40e-9, dt_s=25e-12, n_save=201,
             n_profiles=5, verbose=False)


@pytest.fixture(scope="module")
def short_runs(tmp_path_factory):
    out = str(tmp_path_factory.mktemp("runs1d"))
    p = sim1d.Params1D(**SHORT, output_dir=out, name="short8")
    band = p.replace(rates=(("arf_yield", 1.0), ("ar2_rad", 6.0e7)),
                     name="short8_bands")
    res = {"8": sim1d.run(p),
           "16": sim1d.run(p.replace(cells=16, name="short16")),
           "dt50": sim1d.run(p.replace(dt_s=50e-12, name="short8_dt50")),
           "band": sim1d.run(band)}
    res["0d"] = model.run(p.cell_params().replace(n_save=4001))
    res["band0d"] = model.run(band.cell_params().replace(n_save=4001))
    beam = p.replace(drive="beam", name="short8_beam")
    res["beam"] = sim1d.run(beam)
    res["beam0d"] = model.run(beam.cell_params().replace(n_save=4001))
    return res


def test_centre_cell_is_the_0d_model(short_runs):
    """With uniform power the gap's interior is v0_2's 0-D run: the centre
    cell starts in v0_2's state, and its emission by band, ArF formed,
    species, Tg and (after the first 2 ns) Te follow v0_2's."""
    res, r0 = short_runs["8"], short_runs["0d"]
    report = _m("report")
    first = res.data["records"][0]["centre"]
    for i, s in enumerate(res.species):
        v0 = r0.ne[0] if s == "e-" else r0.n[s][0]
        assert first[i] == approx(v0, 1e-12), s
    d0 = report.row_dict(r0)
    assert res.share_193("ledger_centre") == approx(d0["emission_193"], 1e-3)
    acc = report.accounting(r0)
    led = res.data["records"][-1]["ledger_centre"]
    for b in res.data["bands"]:
        assert led[f"rad:{b}"] / led["input"] == approx(
            acc["bands"][b] / acc["e_dep"], 1e-3), b
    assert led["arf_formed"] == approx(
        report.reservoir(r0, "ArF*")["formed"], 1e-3)
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
    """Two bands at their tops, the ArF* yield of Ar* + F2 (1.0) and the
    Ar2* radiative rate (full singlet-triplet mixing, 6e7 1/s): the centre
    cell follows v0_2's 0-D run with the same bands, and that run differs
    from the default by far more than the 1-D and 0-D runs differ from each
    other."""
    report = _m("report")
    res, r0, d0 = short_runs["band"], short_runs["band0d"], short_runs["0d"]
    band, dflt = report.row_dict(r0), report.row_dict(d0)
    got = res.share_193("ledger_centre")
    assert got == approx(band["emission_193"], 1e-3)
    led = res.data["records"][-1]["ledger_centre"]
    vuv = 100 * led["rad:127 nm (Ar2* VUV)"] / led["input"]
    assert vuv == approx(band["emission_vuv"], 1e-3)
    assert res.series("Ar2*", "centre")[-1] == approx(r0.n["Ar2*"][-1], 1e-2)
    assert abs(band["emission_193"] / dflt["emission_193"] - 1) > BAND_193
    assert abs(band["emission_vuv"] / dflt["emission_vuv"] - 1) > BAND_VUV
    assert "bands: arf_yield 1, ar2_rad 6e+07" in res.summary()


def test_centre_cell_is_the_0d_model_under_the_beam(short_runs):
    """Under the electron beam the gap's interior is v0_2's 0-D beam run:
    the centre cell's emission by band, ArF formed, the beam's prompt heat,
    species, Tg and (after the first 2 ns) Te follow it as the own input's
    follow v0_2's own-input run."""
    res, r0 = short_runs["beam"], short_runs["beam0d"]
    report = _m("report")
    assert "drive beam" in res.data["model"]
    d0 = report.row_dict(r0)
    assert res.share_193("ledger_centre") == approx(d0["emission_193"], 1e-3)
    acc = report.accounting(r0)
    led = res.data["records"][-1]["ledger_centre"]
    for b in res.data["bands"]:
        assert led[f"rad:{b}"] / led["input"] == approx(
            acc["bands"][b] / acc["e_dep"], 1e-3), b
    assert led["arf_formed"] == approx(
        report.reservoir(r0, "ArF*")["formed"], 1e-3)
    assert led["heat_beam"] / led["input"] == approx(
        r0.ledger["heat_beam"][-1] / r0.ledger["E_dep"][-1], 1e-9)
    last = res.data["records"][-1]
    assert last["tg_centre"] == pytest.approx(d0["tg_end_K"], abs=0.5)
    for s in res.species:
        v0 = r0.ne[-1] if s == "e-" else r0.n[s][-1]
        assert res.series(s, "centre")[-1] == approx(v0, 1e-2), s
    te1 = np.array([r["te_centre"] for r in res.data["records"]])
    te0 = np.interp(res.t, r0.t, r0.te)
    late = res.t > 2e-9
    assert np.max(np.abs(te1[late] / te0[late] - 1)) < 0.01
    # and the beam is not the own input: during the pulse its electrons
    # sit far cooler (after it both cool to the same afterglow), and the
    # 193-nm share differs
    own = short_runs["8"]
    on = (res.t > 5e-9) & (res.t < 20e-9)
    te_own = np.array([r["te_centre"] for r in own.data["records"]])
    assert np.mean(te1[on]) < 0.5 * np.mean(te_own[on])
    assert abs(res.share_193("ledger_centre")
               / own.share_193("ledger_centre") - 1) > 0.05


def test_energy_ledger_closes(short_runs):
    """The gap's chemical + electron + thermal energy, recomputed from its
    state, equals what the ledger booked (input - radiation + bounds +
    transport), and the fate adds up to 100 %, on either drive."""
    for key in ("8", "16", "dt50", "band", "beam"):
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
        rad = np.array(led["rad:193 nm (ArF B-X)"])
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
    assert back.share_193() == res.share_193()
    assert back.data["provenance"]["chemistry_md5"] == sim1d.chemistry_md5()


def test_cli(tmp_path, capsys):
    main = _m("run").main
    assert main(["--cells", "4", "--pulse", "0.002", "--t-end", "0.004",
                 "--out", str(tmp_path), "--name", "cli", "--quiet",
                 "--rate", "ar2_f=1e-10", "--compare"]) == 0
    out = capsys.readouterr().out
    assert "radiated 193 nm (ArF B-X)" in out and "0-D (v0_2's model)" in out
    assert "bands: ar2_f 1e-10;  no measurement, 0: ar2_pen, ar2f_ar, " \
        "ar_pen" in out
    assert "127 nm, % of deposited" in out
    assert (tmp_path / "cli" / "result.json").is_file()
    assert (tmp_path / "cli" / "fluid1d_inputs" / "reactions"
            / "light_1Dv0_2.toml").is_file()
    with pytest.raises(SystemExit):
        main(["--rate", "a2f=3e-10"])
    assert main(["--cells", "4", "--pulse", "0.002", "--t-end", "0.004",
                 "--out", str(tmp_path), "--name", "clibeam", "--quiet",
                 "--drive", "beam", "--compare"]) == 0
    out = capsys.readouterr().out
    assert "by the electron beam" in out and "drive beam" in out
    assert "heat (beam: above 4s, F2 fragments)" in out


@pytest.mark.skipif(not os.environ.get("LIGHT_1D_SLOW"),
                    reason="set LIGHT_1D_SLOW=1 (about a minute)")
def test_reference_point_centre_is_v0_2(tmp_path):
    """v0_2's reference point (1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, 2 us)
    on 32 cells: the centre cell's 193-nm share and Tg are the 0-D run's."""
    p = sim1d.Params1D(cells=32, verbose=False, output_dir=str(tmp_path))
    res = sim1d.run(p)
    d0 = _m("report").row_dict(model.run(p.cell_params()))
    assert res.share_193("ledger_centre") == approx(d0["emission_193"], 1e-4)
    assert res.data["records"][-1]["tg_centre"] == pytest.approx(
        d0["tg_end_K"], abs=0.1)
    assert abs(res.fate()["closure residual"]) < 1e-8
