"""LIGHT v0_3-v0_5's kinetic EEDF on pure krypton: the solver
(boltzmann.py) and the run with it in the loop (eedf='kinetic', the
default).

The solver is checked against ML's dense solver at the same inputs,
against its own particle and energy balances, in the undriven limit (the
Maxwellian at the gas temperature), and on the S curve of krypton's
Ramsauer minimum under strong e-e collisions, where the problem at a fixed
field has several solutions and a fixed-field iteration finds at most the
upper one, but the problem at a fixed mean energy has one. The run is
checked for the ledger's closure and for the agreement of its kinetic
energy balance with LIGHT's.

LIGHT's versions of these tests ran on Kr/F2 laser states. The states here
are those states with the F2 taken out (2026-10-13), and every assertion
was measured again on them. Pure krypton differs in three places: no
attachment, so no state here needs the growth term dropped for it; the S
curve sits near 0.04 Td, not at LIGHT's 0.4 Td, where F2's low-energy
losses had held a cold branch; and the grid search can land on two
different grids in its window from two starts.

The reference point's convergence in the update settings and the grid is
slow (minutes): set LIGHT_KINETIC_SLOW=1."""

import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from v0_4 import boltzmann as KB
from v0_4 import xsec
from v0_4.model import QE, Model, Params
from v0_4.network import build
from v0_4.species import e_form
from v0_4.report import accounting, electrons_line, row_dict

ROOT = Path(__file__).resolve().parents[3]          # Kinetic_Modeling
ML_BOLTZ = ROOT / "ML" / "ml" / "boltzmann"
LOKI_INPUT = ROOT / "LoKI-MC" / "Code" / "Input"
SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_KINETIC_SLOW"),
                          reason="minutes (reference-point runs): set "
                          "LIGHT_KINETIC_SLOW=1")
FAST = dict(pulse_s=0.3e-6, t_end_s=0.6e-6)


def _procs():
    out = []
    for r in build():
        if r.rate[0] == "xsec":
            p = xsec.find(*r.rate[1:4])
            out.append(KB.Process(r.rate[1], p.kind, p.threshold, p,
                                  (r.rid,)))
    return out


def _rec(n_pow=-0.55, k0=1.6e-6, te0=0.025851999786):
    """The R11p template's constant, as model.py makes it [m^3/s eV^-n]."""
    return (k0 * te0 ** (-n_pow) * math.gamma(1.5)
            / math.gamma(1.5 + n_pow) * 1e-6)


def _plasma(pen):
    """A driven krypton plasma (450 K, n_e/N = 3e-6, LIGHT's laser state
    without its 0.42 % F2) with Penning electrons at pen per electron per
    second."""
    return KB.State(x={"Kr": 1.0}, tg=450.0, n_col=2.45e25, ne=7.4e19,
                    pen=((pen, 5.834),), rec=_rec() * 7.2e19,
                    rec_pow=-0.55)


# LIGHT's reference point 0.12 ns into the pulse without its F2: Te 0.2 eV
# and rising, n_e/N 3.4e-5
TRANSIENT = KB.State(x={"Kr": 1.0}, tg=300.0, n_col=2.446e25, ne=8.356e20,
                     pen=((0.221, 5.834),), rec=3.054e6, rec_pow=-0.55)


def _energy(s, st, sol, procs):
    """(gain, loss, net frequency) per electron [eV/s] of a solution: field
    and Penning injection against elastic, inelastic (thresholds), the
    recombined electrons' energy, and nu <eps>. (LIGHT's also took the
    attached electrons'.)"""
    g, F, N = sol.grid, sol.F, st.n_col
    ion, att, rec = s._frequencies(g, st, F)
    assert att == 0.0                               # no attachment in Kr
    inel = sum(st.x[q.gas] * N * g.rate(q.p, F) * q.threshold for q in procs
               if q.kind in ("excitation", "ionization"))
    el = N * sum(st.x[gas] * s.elastic_power(st, sol, gas)
                 for gas in s.elastic)
    wr = g.ec ** (st.rec_pow + 0.5) * F
    gain = s.field_power(st, sol) + st.source_power()
    loss = (inel + el + rec * float((g.ec * wr).sum() / wr.sum())
            + sol.nu * sol.mean)
    return gain, loss, ion + sum(v for v, _ in st.pen) - att - rec


def _loki_input(tmp_path):
    """A LoKI-MC input directory for ML's solver with v0_2's cross sections:
    LXCat's Biagi 8.97 file rewritten in the block form LoKI-MC and ML's
    lx_lib parse (the collision named in a bracketed COMMENT after PARAM,
    LF line ends, plain state names Kr(L01) ... Kr(L51)), beside LoKI-MC's
    own Databases and the two tables gases_twin.load_tables reads."""
    inp = tmp_path / "Input"
    (inp / "Krypton").mkdir(parents=True)
    for d in ("Databases", "Nitrogen", "H2O"):
        (inp / d).symlink_to(LOKI_INPUT / d, target_is_directory=True)
    src = (xsec.DATA_DIR / xsec.FILES["Kr"][0]).read_text()
    lines = src.replace("\r\n", "\n").split("\n")
    out, i, nlev = [], 0, 0
    while i < len(lines):
        kw = lines[i].strip()
        if kw not in ("ELASTIC", "EXCITATION", "IONIZATION"):
            i += 1
            continue
        thr = lines[i + 2]
        j = i + 3
        while not lines[j].lstrip().startswith("---"):
            j += 1
        k = j + 1
        while not lines[k].lstrip().startswith("---"):
            k += 1
        table = lines[j + 1:k]
        if kw == "ELASTIC":
            head = ["ELASTIC", "Kr", thr, "SPECIES: e / Kr",
                    "PROCESS: E + Kr -> E + Kr, Elastic",
                    "PARAM.:  m/M = 0.0000065, complete set",
                    "COMMENT: [Kr(1S0) + e -> Kr(1S0) + e, Elastic]"]
        elif kw == "EXCITATION":
            nlev += 1
            st = f"Kr(L{nlev:02d})"
            head = ["EXCITATION", f"Kr -> {st}", thr, "SPECIES: e / Kr",
                    f"PROCESS: E + Kr -> E + {st}, Excitation",
                    f"PARAM.:  E = {float(thr):.6g} eV, complete set",
                    f"COMMENT: [Kr(1S0) + e -> {st} + e, Excitation]"]
        else:
            head = ["IONIZATION", "Kr -> Kr^+", thr, "SPECIES: e / Kr",
                    "PROCESS: E + Kr -> E + E + Kr+, Ionization",
                    f"PARAM.:  E = {float(thr):.6g} eV, complete set",
                    "COMMENT: [Kr(1S0) + e -> Kr(+,gnd) + e + e, "
                    "Ionization]"]
        out += head + ["UPDATED: 2011-10-18 00:00:00",
                       "COLUMNS: Energy (eV) | Cross section (m2)",
                       "-----------------------------"] + table + [
                       "-----------------------------", ""]
        i = k + 1
    assert nlev == 51
    (inp / "Krypton" / "Kr_LXCat.txt").write_text("\n".join(out) + "\n")
    return inp


# ------------------------------------------------------------- the solver
@pytest.mark.skipif(not (ML_BOLTZ.is_dir() and LOKI_INPUT.is_dir()),
                    reason="ML's solver or LoKI-MC's inputs not present")
@pytest.mark.parametrize("en,ne_n,tg", [(10.0, 1e-5, 300.0),
                                        (20.0, 3e-6, 300.0),
                                        (30.0, 1e-4, 600.0)])
def test_solver_matches_ml_dense_solver(en, ne_n, tg, tmp_path):
    """Every process of krypton's file, a DC field, e-e and the growth
    term: with ML's e-e quadrature the sparse solver is ML's
    (ee.solve_field, numpy.linalg.solve on the whole matrix) to round-off,
    and LIGHT's own quadrature moves F by 1e-4 at most at these driven
    states. v0_2: on v0_2's own Biagi 8.97 file, rewritten for ML's
    parser (_loki_input); v0_1 compared on Biagi-v7.1 as LoKI-MC holds it.
    v0_3: the same file, so v0_3's side runs on it as published
    (xsec.UNIT): the test is of the solver, not of the factors."""
    sys.path.insert(0, str(ML_BOLTZ))
    import ee
    import gases_twin as GT
    inp = str(_loki_input(tmp_path))
    GT.load_tables(inp)
    n_col, N, emax = 2.4e25, 600, 40.0
    kr = GT.build("Kr", tg, inp)
    assert kr["n_blocks"] == 53 and len(kr["colls"]) == 52   # + elastic
    ml = ee.solve_field([(kr, 1.0)], en, tg, N, emax, ne_N=ne_n,
                        ne_m3=ne_n * n_col, tol=1e-11)
    with xsec.using(xsec.UNIT):
        _against_ml(ml, kr, ee, en, tg, ne_n, n_col, N, emax)


def _against_ml(ml, kr, ee, en, tg, ne_n, n_col, N, emax):
    procs = [KB.Process("Kr", p.kind, p.threshold, p)
             for p in xsec.processes("Kr")
             if p.kind in ("excitation", "ionization", "attachment")]
    st = KB.State(x={"Kr": 1.0}, tg=tg, n_col=n_col, ne=ne_n * n_col)
    mine = KB.Solver(procs, n_cells=N, emax=emax, tol=1e-11, maxit=4000,
                     ee_quadrature="cell").solve(st, en)
    assert np.max(np.abs(mine.F - ml["F"])) / ml["F"].max() < 1e-7
    m_ml = float(np.sum(ml["ec"] ** 1.5 * ml["F"]) * ml["du"])
    assert mine.mean == pytest.approx(m_ml, rel=1e-9, abs=0.0)
    assert mine.nu / n_col == pytest.approx(ml["nubar"], rel=1e-7, abs=0.0)
    for thr in (9.9152, 13.9996):
        c = next(c for c in kr["colls"] if abs(c["thr"] - thr) < 1e-6)
        p = next(q.p for q in procs if abs(q.threshold - thr) < 1e-6)
        assert mine.grid.rate(p, mine.F) == pytest.approx(ee.rate(c, ml),
                                                          rel=1e-6, abs=0.0)
    own = KB.Solver(procs, n_cells=N, emax=emax, tol=1e-11,
                    maxit=4000).solve(st, en)
    assert np.max(np.abs(own.F - ml["F"])) / ml["F"].max() < 1e-4


@pytest.mark.parametrize("cells", [600, 2400])
def test_undriven_eedf_is_the_maxwellian_at_the_gas_temperature(cells):
    """No field, no sources, n_e/N = 1e-5: elastic collisions pull toward
    Tg and e-e collisions leave any Maxwellian alone, so the answer is the
    grid's Maxwellian at Tg to round-off. With F constant in each cell in
    the e-e integrals (ML's quadrature) the lowest boundaries' error
    outweighs the elastic coupling and the EEDF condenses into the first
    cell (its mean falls with the cell width)."""
    st = KB.State(x={"Kr": 1.0}, tg=300.0, n_col=2.45e25, ne=2.45e20)
    s = KB.Solver(_procs(), n_cells=cells, emax=12.0)
    sol = s.solve_joint(st, None, None, None, None, True)
    FM = sol.grid.maxwellian(KB.KB * 300.0 / KB.QE)
    bulk = sol.grid.ec < 0.5
    assert np.max(np.abs(sol.F[bulk] / FM[bulk] - 1.0)) < 1e-10
    assert sol.nu == 0.0
    # and it exchanges no energy with the gas, as the solver exchanges it
    g = sol.grid
    cooling = KB.GAMMA * float(np.sum(g.eb ** 2 * g.sig_eps["Kr"] * 0.5
                                      * (sol.F[1:] + sol.F[:-1]))) * g.du
    assert abs(s.elastic_power(st, sol, "Kr")) < 1e-12 * cooling
    cell = KB.Solver(_procs(), n_cells=cells, emax=12.0,
                     ee_quadrature="cell")
    bad = cell.solve_joint(st, None, None, None, None, True)
    assert bad.mean < 0.8 * sol.mean


def test_particle_and_energy_balance():
    """nu is the solution's own net frequency, and the energy the field and
    the Penning electrons bring equals what the collisions, the
    recombination and the growth term take (2e-6 measured). Without the
    growth term (no field at a state the solver drops it at) the net gain
    of electrons comes back in the first cell, and with that the balance
    closes as well (3.5e-4 on 0.05-eV cells). The solver is LIGHT's: Kr's
    elastic loss as the network books it (G6)."""
    procs = _procs()
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    assert Model(Params()).kin.elastic == ("Kr",)
    for st, target in ((_plasma(1.9e7), 4.5), (_plasma(6.4e7), 4.5)):
        sol = s.eedf(st, target)
        assert not sol.growth_off and sol.reached
        gain, loss, nu = _energy(s, st, sol, procs)
        assert sol.nu == pytest.approx(nu, rel=1e-9, abs=0.0)
        assert abs(loss - gain) / gain < 2e-5
    st = _plasma(2.1e7)
    sol = s.eedf(st, None)
    assert sol.growth_off and sol.en_td == 0.0
    gain, loss, net = _energy(s, st, sol, procs)
    injected = -net * sol.grid.ec[0]
    assert abs(loss - gain - injected) / gain < 1e-3


def test_powers_are_the_solvers_own_flux_split():
    """The Scharfetter-Gummel flux through each inner boundary is the
    central W F_b - D dF/deps with its diffusive part times c(z) = (z/2)
    coth(z/2), z = W du/D from all of W and D (elastic, e-e and field). The
    elastic and field powers are its parts, the energy the solver itself
    exchanges. At a thermal state with a weak field c(z) reaches 1.05 on
    600 cells, by which central differences miss."""
    s = KB.Solver(_procs(), n_cells=600, emax=12.0, elastic=("Kr",))
    st = KB.State(x={"Kr": 1.0}, tg=300.0, n_col=2.446e25, ne=2.4e19)
    cold = s.solve_joint(st, None, None, None, 0.0, True)
    sol = s.eedf(st, 1.2 * cold.mean)
    assert sol.reached and sol.en_td > 0.0
    g, F = sol.grid, sol.F
    _, W, D0, dfield, hf = s._parts(g, st, sol.nu, F)
    assert hf[0] is None                                   # no EM wave
    D = D0 + (sol.en_td * KB.TD) ** 2 * dfield
    z = W * g.du / D
    flux = D / g.du * (KB.bern(-z) * F[:-1] - KB.bern(z) * F[1:])
    fb, dfb = 0.5 * (F[1:] + F[:-1]), (F[1:] - F[:-1]) / g.du
    _, c, _ = s._split(st, sol)
    assert np.max(np.abs(flux - (W * fb - D * c * dfb))) \
        < 1e-10 * np.max(np.abs(W * fb))
    assert c.max() > 1.03
    en_v = sol.en_td * KB.TD
    assert s.field_power(st, sol) == pytest.approx(
        -en_v ** 2 * st.n_col * float(np.sum(dfield * c * dfb)) * g.du,
        rel=1e-12, abs=0.0)
    kt = KB.KB * 300.0 / KB.QE
    assert s.elastic_power(st, sol, "Kr") == pytest.approx(
        KB.GAMMA * float(np.sum(g.eb ** 2 * g.sig_eps["Kr"]
                                * (fb + kt * c * dfb))) * g.du,
        rel=1e-12, abs=0.0)


def test_time_step_limits_and_its_energy():
    """Solver.step is one backward-Euler step of the time-dependent EEDF.
    A step too long to matter is the steady growth solution, and one too
    short to matter, with no field, leaves F where it was. From the gas's
    Maxwellian under a Joule power of 1e9 eV/s per electron the field takes
    that power, F stays normalised, nu is its own net frequency, and with
    e-e left out the mean rises by what the field and the sources bring
    less what the collisions, captures and population change take."""
    procs = _procs()
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    st = _plasma(1.9e7)
    steady = s.solve_joint(st, 4.5, None, None, None, True)
    big = s.step(st, steady, 1e3, target_mean=4.5, adapt=False)
    assert np.max(np.abs(big.F - steady.F)) < 1e-10 * steady.F.max()
    assert big.nu == pytest.approx(steady.nu, rel=1e-8, abs=0.0)
    tiny = s.step(st, steady, 1e-20, adapt=False)
    assert np.max(np.abs(tiny.F - steady.F)) < 1e-8 * steady.F.max()
    # a target colder than 0.1 ns can reach is not reached, and the step
    # with no field stands, converged
    low = s.step(st, steady, 1e-10, target_mean=4.0, adapt=False)
    assert not low.reached and low.en_td == 0.0 and low.mean > 4.0
    again = s.solve_joint(st, None, low.F, low.nu, en_fixed=0.0,
                          memory=(steady.grid.ec, steady.F, 1e-10))
    assert np.max(np.abs(again.F - low.F)) < 1e-8 * low.F.max()
    cold_state = KB.State(x=st.x, tg=st.tg, n_col=st.n_col, ne=0.0,
                          pen=st.pen, rec=st.rec, rec_pow=st.rec_pow)
    g = s.grid
    F0 = g.maxwellian_at_mean(1.5 * KB.KB * st.tg / KB.QE)
    cold = KB.Solution(F0, 0.0, 0.0, g.mean(F0), False, 0, g)
    for dt in (1e-12, 1e-10):
        s.grid = g
        sol = s.step(cold_state, cold, dt, target_power=1e9, adapt=False)
        assert s.field_power(cold_state, sol) == pytest.approx(
            1e9, rel=1e-9, abs=0.0)
        assert float(np.sum(g.sqrt_ec * sol.F) * g.du) == pytest.approx(
            1.0, rel=1e-12, abs=0.0)
        gain, loss, nu = _energy(s, cold_state, sol, procs)
        assert sol.nu == pytest.approx(nu, rel=1e-9, abs=0.0)
        rise = (sol.mean - cold.mean) / dt
        assert abs(gain - loss - rise) < 1e-4 * gain


def test_remap_keeps_the_mean():
    """A history moved onto another grid keeps its normalisation and its
    mean energy, and the grid's Maxwellian at a mean has that mean on the
    grid."""
    s = KB.Solver(_procs(), n_cells=600, emax=30.0)
    g30, g12, g40 = s.grid, s._grid(12.0), s._grid(40.0)
    F = g30.maxwellian(0.5)
    for g in (g12, g40):
        Fm = g.remap(g30.ec, F)
        assert g.mean(Fm) == pytest.approx(g30.mean(F), rel=1e-12, abs=0.0)
        assert float(np.sum(g.sqrt_ec * Fm) * g.du) == pytest.approx(
            1.0, rel=1e-12, abs=0.0)
    kt = KB.KB * 300.0 / KB.QE
    assert g12.mean(g12.maxwellian_at_mean(1.5 * kt)) == pytest.approx(
        1.5 * kt, rel=1e-12, abs=0.0)
    assert g12.mean(g12.maxwellian(kt)) < 0.98 * 1.5 * kt


def test_on_the_s_curve_a_field_has_several_solutions_a_mean_has_one():
    """The heating transient's state in pure Kr. At a fixed mean energy
    the solver (Solver.eedf, without the growth term) reaches 0.06, 0.2
    and 0.5 eV from cold, warm and hot starts, one solution and one field
    for each, and the field is not monotonic in the mean: 0.2 eV needs
    less than 0.06 eV or 0.5 eV, the middle branch of the S curve. So at
    the field the 0.06-eV solution needs (0.042 Td), a fixed-field
    iteration from a hot start finds another, upper solution, and from the
    cold start none. (In LIGHT's Kr/F2 state the curve sat at 0.4 Td;
    there F2's low-energy losses held its cold branch.)"""
    s = KB.Solver(_procs(), n_cells=600, emax=12.0, elastic=("Kr",))
    g = s.grid
    en = {}
    for target in (0.06, 0.2, 0.5):
        sols = []
        for te in (0.026, 0.1, 2.0):
            s.grid = g
            F0 = g.maxwellian(te)
            warm = KB.Solution(F0, 0.0, 0.0, g.mean(F0), True, 0, g)
            sols.append(s.eedf(TRANSIENT, target, warm=warm, growth=False))
        for sol in sols:
            assert sol.reached and sol.growth_off
            assert sol.mean == pytest.approx(target, rel=1e-7)   # 8e-9 seen
            assert sol.grid is sols[0].grid
            assert np.max(np.abs(sol.F - sols[0].F)) / sols[0].F.max() < 1e-8
            # the field to 1e-4: on the S curve's middle branch the mean
            # hardly moves with the field, so an EEDF equal to 1e-8 leaves
            # it free to 1e-5 (2026-10-13)
            assert sol.en_td == pytest.approx(sols[0].en_td, rel=1e-4)
        en[target] = sols[0].en_td
    assert en[0.2] < en[0.06] and en[0.2] < en[0.5]
    s.grid = g
    hot = s.solve_joint(TRANSIENT, None, g.maxwellian(2.0), None, None,
                        False, en_fixed=en[0.06])
    assert hot is not None and hot.mean > 5.0 * 0.06
    cold = s.solve_joint(TRANSIENT, None, g.maxwellian(0.026), None, None,
                         False, en_fixed=en[0.06])
    assert cold is None


def test_growth_term_kept_only_where_it_gives_the_mean_energy():
    """In the pulse the growth term is kept and the field gives LIGHT's
    mean energy. 0.12 ns into the pulse the solver solves without it, at
    LIGHT's mean energy all the same. Where the Penning electrons alone
    make the growth solution hotter than the target, the field falls to its
    floor and the solve stops within a few steps, not reached."""
    s = KB.Solver(_procs(), n_cells=600, emax=30.0, elastic=("Kr",))
    sol = s.eedf(_plasma(1.9e7), 4.5)
    assert not sol.growth_off and sol.reached and sol.en_td > 0.0
    assert sol.mean == pytest.approx(4.5, rel=1e-9)
    sol = s.eedf(TRANSIENT, 0.29)
    # the inner match is held to 1e-3 of the outer step's change, so the
    # mean lands within a few 1e-9 (v0_2: 2.2e-9; v0_1 within 1e-9)
    assert sol.growth_off and sol.mean == pytest.approx(0.29, rel=1e-8)
    sol = s.solve_joint(_plasma(6.4e7), 2.7, None, None, None, True)
    assert sol is not None and not sol.reached and sol.en_td == 0.0
    assert sol.mean > 2.7 and sol.iters <= 10


def test_growth_validity_threshold():
    """growth_valid keeps a net loss nu only while the energy it moves,
    |nu| <eps>, is at most ten times the heating, and keeps any net gain.
    With no field the heating is the Penning electrons' alone, which nu
    does not change (a field's power does, through sigma~)."""
    s = KB.Solver(_procs(), n_cells=600, emax=30.0, elastic=("Kr",))
    # the same power brought by Penning electrons and by a beam's
    # secondaries (whose power counts as heating as well)
    beam = KB.State(x={"Kr": 1.0}, tg=450.0, n_col=2.45e25, ne=7.4e19,
                    rec=_rec() * 7.2e19, rec_pow=-0.55,
                    inj=((2.1e7 * 5.834 / 4.506, 9.915, 4.506),))
    for st in (_plasma(2.1e7), beam):
        sol = s.eedf(st, None)
        assert sol.en_td == 0.0 and s.field_power(st, sol) == 0.0
        edge = 10.0 * 2.1e7 * 5.834 / sol.mean
        assert 10.0 * edge > s.relax_freq(st, sol)
        for nu, ok in ((-0.9 * edge, True), (-1.1 * edge, False),
                       (+1.1 * edge, True)):
            sol.nu = nu
            assert s.growth_valid(st, sol) is ok


def test_recombining_afterglow_keeps_the_growth_term():
    """Pure Kr after the pulse, with no heating, Kr2+ recombining the
    electrons thirty times slower than the gas relaxes their energy. The
    growth term stays, since the relaxation and not the turnover sets the
    shape. Where LIGHT's mean energy lies below its no-field solution, that
    solution stands, unreached, though the solution without the term
    reaches the target by putting the recombined electrons back cold
    (2026-09-30: dropping it took a pure-Kr afterglow below the gas
    temperature). So from a cold start and from one without the term."""
    s = KB.Solver(_procs(), n_cells=600, emax=12.0, elastic=("Kr",))
    st = KB.State(x={"Kr": 1.0}, tg=300.0, n_col=2.446e25,
                  ne=2.4e19, rec=_rec() * 1e17, rec_pow=-0.55)
    grow = s.solve_joint(st, None, None, None, 0.0, True)
    cold = s.solve_joint(st, None, None, None, 0.0, False)
    assert grow.nu < 0.0 and -30.0 * grow.nu < s.relax_freq(st, grow)
    assert cold.mean < grow.mean
    # the second target needs almost no field without the term, so the
    # turnover outweighs the heating and only the relaxation keeps it
    for f in (0.5, 0.02):
        target = cold.mean + f * (grow.mean - cold.mean)
        off = s.eedf(st, target, growth=False)
        assert off.growth_off and off.reached
        if f < 0.1:
            ion, att, rec = s._frequencies(off.grid, st, off.F)
            nu = ion + st.source_rate() - att - rec
            assert -nu * off.mean > 10.0 * s.field_power(st, off)
        for warm in (None, off):
            sol = s.eedf(st, target, warm=warm)
            assert not sol.growth_off and not sol.reached
            assert sol.en_td == 0.0
            assert sol.mean == pytest.approx(grow.mean, rel=1e-4)


def test_beam_secondaries_reach_the_grid():
    """A beam's secondary electrons, injected with the sub-excitation
    spectrum (linear on [0, 9.915] eV, mean 4.506 eV, as ML's solve_beam
    injects them), arrive as one electron per injected electron, with the
    injected mean energy, on grids from 12 to 40 eV."""
    for emax in (12.0, 21.26, 40.0):
        g = KB.Grid(600, emax, _procs())
        q = g.injection(9.915, 4.506)
        assert q.sum() == pytest.approx(1.0, rel=1e-12, abs=0.0)
        assert float(np.sum(q * g.ec)) == pytest.approx(4.506, rel=2e-5,
                                                         abs=0.0)
        assert q[g.ec > 9.915 + g.du].sum() == 0.0


# a beam drive as a test FIXTURE: LIGHT's prototype
# (LIGHT/docs/survey_2026-09-27/scripts/ebeam_screen.py), whose per-pair
# yields are REBUILD's Spencer-Fano numbers for Kr with 0.5 % F2 and whose
# Kr** cascade is booked at LIGHT v0's 11.6 eV. It exercises the kinetic
# equation's electron source and the run's stall guard (LIGHT's Known-bad
# list); it is not a pure-Kr beam model
# v0_2: the ledger's own energies (Magboltz 8.97's thresholds), so the
# fixture's books close against v0_2's species; v0_4: its excitations go
# to the lowest 5s level, Kr*(1s5), at the energy v0_3's lumped Kr* had
KRP, KRS = e_form("Kr+"), e_form("Kr*(1s5)")
Y_ION, Y_PP, Y_S, EPS_SEC = 0.9989, 0.5119, 0.2634, 4.506
CASCADE = Y_PP * (11.6 - KRS)
W_EFF = Y_ION * KRP + (Y_PP + Y_S) * KRS + CASCADE + EPS_SEC   # eV per pair


class BeamFixture(Model):
    def rhs(self, t, y, on):
        if not on:
            return super().rhs(t, y, False)
        dy = super().rhs(t, y, False)          # chemistry, no field heating
        N0 = self.N0
        pairs = self.p.power_W_cm3 / QE / W_EFF  # pairs cm^-3 s^-1
        dy[self.idx["Kr"]] -= pairs * (Y_ION + Y_PP + Y_S) / N0
        dy[self.idx["Kr+"]] += pairs * Y_ION / N0
        dy[self.idx["Kr*(1s5)"]] += pairs * (Y_PP + Y_S) / N0
        dy[self.i_ue] += pairs * EPS_SEC / N0
        n = np.maximum(y[:self.ns], 0.0) * N0
        tg_ev = y[self.i_tg] * 8.617333262e-5
        c_v = float((n * self.cv).sum()) * 8.617333262e-5
        dn_src = np.zeros(self.ns)
        dn_src[self.idx["Kr"]] = -pairs * (Y_ION + Y_PP + Y_S)
        dn_src[self.idx["Kr+"]] = pairs * Y_ION
        dn_src[self.idx["Kr*(1s5)"]] = pairs * (Y_PP + Y_S)
        dy[self.i_tg] += (pairs * CASCADE
                          - tg_ev * float((dn_src * self.cv).sum())) / c_v
        dy[self.i_led + 0] += pairs * W_EFF / N0          # E_dep
        dy[self.i_led + 2] += pairs * CASCADE / N0        # the cascade heat
        return dy

    def kinetic_injection(self, t, y, on):
        if not on:
            return ()
        return ((self.p.power_W_cm3 / QE / W_EFF, KRS, EPS_SEC),)


def test_beam_seed_grows_on_its_secondaries():
    """A 1e10 cm^-3 seed at 300 K in 665 Torr of Kr with the fixture's
    secondaries injected at 7.6e12 per electron per second. Without the
    growth term the first cell would have to take up a net gain, so the
    growth term stays, its frequency is the injection's, and the EEDF is
    the injected one, hotter than LIGHT's electrons at the gas temperature
    (not reached: no field can cool it)."""
    pairs = 3.3e5 / QE / W_EFF
    ne = 1e10
    n_col = 665 / 760 * 101325 / (1.380649e-23 * 300.0)
    st = KB.State(x={"Kr": 1.0}, tg=300.0, n_col=n_col, ne=ne * 1e6,
                  inj=((pairs / ne, 9.915, 4.506),))
    s = KB.Solver(_procs(), n_cells=600, emax=30.0, elastic=("Kr",))
    sol = s.eedf(st, 1.5 * 0.025852)
    assert not sol.growth_off and not sol.reached and sol.en_td == 0.0
    assert sol.nu == pytest.approx(pairs / ne, rel=1e-3)
    assert sol.mean == pytest.approx(4.506, rel=1e-2)


def test_grid_is_sized_from_the_tail():
    """From a grid long enough and one far too short the solver returns a
    grid in its window (F falling 11 to 25 decades across it), and the
    solver's grid is always the returned solution's. In pure Kr on v0_1's
    Biagi-v7.1 groups the two starts landed on two different grids of the
    window (19.3 and 25.7 eV), whose rates differed by the grids'
    resolution (0.5 % on the 5s excitation, 2.2 % on the ionization); on
    v0_2's 51 levels both land on one grid, 21.3 eV (2026-10-14)."""
    out = []
    for emax in (30.0, 13.2):
        s = KB.Solver(_procs(), n_cells=600, emax=emax, elastic=("Kr",))
        sol = s.eedf(_plasma(1.9e7), 4.5)
        assert sol.grid is s.grid
        assert 11.0 <= KB._decades(sol.F) <= 25.0
        out.append([sol.grid.rate(q.p, sol.F) for q in _procs()])
    assert out[0] == pytest.approx(out[1], rel=0.03, abs=0.0)


def test_grid_search_from_a_grid_five_times_too_long():
    """v0_1 (strict xfail there): from 60 eV (and 40 eV) the first
    solution's tail reached the noise floor (the top cells clipped to 0),
    the grid of the window it was shortened to repeated that from the long
    grid's warm start, and the search returned the long grid, 300 decades
    on 60 eV, its ionization rate 31 % off the windowed grid's. On v0_2's
    51 levels the search lands on the windowed grid from 60, 40 and 30 eV
    at this state, with one ionization rate (2026-10-14). The solver is
    unchanged, so the rule stays (LIGHT_Kr's Known-bad list): runs start
    at 30 eV and warm-start from the last grid. v0_3's trained cross
    sections and Biagi's file as published both land all three on the
    21.3-eV grid (2026-10-15; a nine-factor trial set had put the 30-eV
    start on a 19.3-eV grid, 1.2 % off in this far-tail rate, inside the
    window: the closure's own spread)."""
    for f in (xsec.factors(), xsec.UNIT):
        with xsec.using(f):
            ion = next(q.p for q in _procs() if q.kind == "ionization")
            got = []
            for emax in (60.0, 40.0, 30.0):
                s = KB.Solver(_procs(), n_cells=600, emax=emax,
                              elastic=("Kr",))
                sol = s.eedf(_plasma(1.9e7), 4.5)
                assert 11.0 <= KB._decades(sol.F) <= 25.0, emax
                got.append(sol.grid.rate(ion, sol.F))
        assert got[0] == pytest.approx(got[2], rel=1e-6, abs=0.0)
        assert got[1] == pytest.approx(got[2], rel=1e-6, abs=0.0)


# ------------------------------------------------------------ the run
def test_beam_drive_feeds_its_secondaries_to_the_kinetic_eedf():
    """Under the beam fixture the kinetic state carries the secondaries
    (one per ion pair, the drive's own numbers), the ledger closes, once
    the seed has grown LIGHT's mean energy is the kinetic one, and the
    kinetic energy balance is LIGHT's with the beam's power injected, not
    carried by the field. As LIGHT v0_4 found in Kr/F2 (1.9e-4 between 10
    and 16 ns), while the Kr* and Kr** populations the 5s <-> 5p pair acts
    on still grow fast between renewals the no-field mean sits above
    LIGHT's Te: in pure Kr up to 1.6e-4, from 15.5 to 21.5 ns, and to 1e-9
    from 22 ns on (2026-10-13)."""
    p = Params(pressure_atm=665 / 760, power_W_cm3=3.3e5, pulse_s=30e-9,
               t_end_s=30e-9, ne_cm3=1e10, n_save=301)
    m = BeamFixture(p)
    inj = m.kinetic_injection(0.0, m.initial_state(), True)
    assert len(inj) == 1
    rate, e_top, mean = inj[0]
    assert rate == pytest.approx(p.power_W_cm3 / QE / W_EFF, rel=1e-12)
    assert (e_top, mean) == (KRS, EPS_SEC)
    assert m.kinetic_injection(0.0, m.initial_state(), False) == ()
    s, _, _ = m._kinetic_state(m.initial_state(), 0.0, True)
    assert s.inj and s.inj[0][0] == pytest.approx(rate / 1e10, rel=1e-12)
    res = m.run()
    assert abs(accounting(res)["closure"]) < 1e-12
    late = [e for e in res.kinetic if e["on"] and not e.get("skipped")
            and e["t"] > 10e-9]
    assert late and all(e["te_kin"] == pytest.approx(e["te"], rel=3e-4)
                        for e in late)
    settled = [e for e in late if e["t"] > 22e-9]
    assert settled and all(e["te_kin"] == pytest.approx(e["te"], rel=1e-6)
                           for e in settled)
    assert max(abs(e["balance"]) for e in late) < 2e-3
    # no field carries a beam's power, so the line gives no ratio
    line = electrons_line(res)
    assert "energy balance gap" in line and "nan" not in line


def _counting(m, limit):
    """Record the times of a model's kinetic updates, and stop the run
    after `limit` of them (a run that stalls would never end)."""
    times, update = [], m._kinetic_update

    def counted(t, y, on):
        times.append((t, y))
        if len(times) > limit:
            raise RuntimeError(f"{limit} kinetic updates, at t = {t:.3e} s")
        return update(t, y, on)
    m._kinetic_update = counted
    return times


def test_beam_run_moves_past_a_trigger_at_its_start():
    """Under a beam, the secondaries heat a seed of 1e10 cm^-3 so fast that
    Te passes its trigger 1.6e-16 s into the run, inside the 8.9e-16 s to
    which scipy locates an event. In pure Kr the event came back at t = 0
    itself, and LIGHT's run renewed its table there 12,405 times without
    moving (2026-09-30). It steps past (the stall guard)."""
    p = Params(pressure_atm=665 / 760, power_W_cm3=3.0e5, pulse_s=3e-15,
               t_end_s=3e-15, ne_cm3=1e10, n_save=3)
    m = BeamFixture(p)
    times = _counting(m, 60)
    res = m.run()
    assert abs(accounting(res)["closure"]) < 1e-12
    assert times[-1][0] > 0.0


def test_rows_dropped_after_construction_leave_the_kinetic_equation():
    """ML drops rows from a model it has built and rescales others. A run's
    kinetic equation carries the rows the model holds when it runs, and the
    table's Maxwellian values at its ends belong to their own rows, though
    a dropped row (C1a_2p10) moves every later row up the list."""
    m = Model(Params(pulse_s=1e-14, t_end_s=1e-14, n_save=3))
    m.compiled = [c for c in m.compiled if c.row.rid != "C1a_2p10"]
    times = _counting(m, 20)
    m.run()
    assert "C1a_2p10" not in {r for pr in m.kin.procs for r in pr.rows}
    tab, y = m._tab, times[-1][1]
    n = np.maximum(y[:m.ns], 0.0) * m.N0
    ne = m.electron_density(n)
    byj = {c.j: c for c in m.compiled}
    ends = np.exp(tab["lnte"][[0, -1]])
    # the other 50 levels, C2, the 5s <-> 5p pairs (v0_4: one per 5s
    # level) and the recombination (the Hornbeck-Molnar rows carry
    # nothing at hm_kr = 0, so are not in the equation)
    levels = [c.row.rid for c in m.compiled
              if c.row.rid.startswith(("G1_", "C1"))]
    assert len(levels) == 50
    pairs = [f"K76{x}_{lv}" for lv in ("1s5", "1s4", "1s3", "1s2")
             for x in ("", "r")]
    assert sorted(m.rows[j].rid for j in tab["lnkM"]) == sorted(
        levels + pairs + ["C2", "R11p"])
    for j, lnk in tab["lnkM"].items():
        want = [math.log(max(m._k(byj[j], e, y[m.i_tg], n, ne), 1e-300))
                for e in ends]
        assert list(lnk) == pytest.approx(want, rel=1e-12, abs=0.0)


def test_table_is_carried_past_its_ends_by_lights_maxwellian():
    """Past the kinetic table's ends LIGHT carries it with its own
    Maxwellian dependence on Te: a rate by its logarithmic slope, a
    captured energy by its ratio and the elastic power by its difference,
    which passes through 0 at Tg. Each meets the table at its end."""
    from v0_4.model import KB_EV
    for memory in (0.0, math.inf):
        _check_table_ends(memory, KB_EV)


def _check_table_ends(memory, KB_EV):
    """The table past its ends, for a steady table (memory 0) and for a
    time step's (one picosecond after the start)."""
    m = Model(Params(kin_memory_s=memory))
    y = m.initial_state()
    y[m.idx["Kr2+"]] = y[m.idx["Kr+"]]            # recombination at work
    n = np.maximum(y[:m.ns], 0.0) * m.N0
    ne = m.electron_density(n)
    tg = y[m.i_tg]
    y[m.i_ue] = 1.5 * ne * 2.0 / m.N0             # the table at Te = 2 eV
    m._kinetic_update(0.0, y, True)
    if memory:
        assert m._tab is None                     # LIGHT's own closure
        m._kinetic_update(1e-12, y, True)
    tab = m._tab
    rows = [c for c in m.compiled if c.j in tab["lnk"]]
    assert len(rows) == 61 and tab["capM"]   # v0_4: four 5s <-> 5p pairs
    lnte = tab["lnte"]
    for end, shift in ((0, -0.3), (-1, 0.3), (0, 0.0), (-1, 0.0)):
        te = math.exp(lnte[end] + shift)
        assert m._tab_end(math.log(te)) == (end if shift else None)
        for c in rows:
            want = tab["lnk"][c.j][end] + math.log(max(
                m._k(c, te, tg, n, ne), 1e-300)) - tab["lnkM"][c.j][end]
            assert math.log(m._rate(c, te, tg, n, ne)) == pytest.approx(
                want, rel=1e-12, abs=1e-9)
        yt = y.copy()
        yt[m.i_ue] = 1.5 * ne * te / m.N0
        d = m.rhs(0.0, yt, True)
        pel = tab["pel"][end] + m._pel_maxwell(te, tg * KB_EV) \
            - tab["pelM"][end]
        assert d[m.i_led + 1] * m.N0 == pytest.approx(
            ne * n[m.idx["Kr"]] * pel, rel=1e-12, abs=0.0)
        cap = 0.0
        for c in m.compiled:
            if c.kind == "elastic" or not c.route["capture"]:
                continue
            flux = m._rate(c, te, tg, n, ne) * np.prod(
                [n[i] for i in c.heavy]) * ne ** (c.ne_pow or 0)
            if c.j in tab["capM"]:
                e_cap = tab["cap"][c.j][end] * c.cap(te) \
                    / tab["capM"][c.j][end]
            elif c.j in tab["cap"]:
                e_cap = tab["cap"][c.j][end]
            else:
                e_cap = 1.5 * te if c.cap is None else c.cap(te)
            cap += e_cap * c.route["capture"] * flux
        assert d[m.i_led + 2] * m.N0 == pytest.approx(cap, rel=1e-12,
                                                       abs=0.0)


def test_first_nanosecond_steps_in_time():
    """The reference pulse's first nanosecond, with the time-dependent step
    (the default in the pulse). The electrons start Maxwellian at the gas
    temperature, under LIGHT's own closure, and every later update is one
    step from the last that reaches LIGHT's mean energy with the growth
    term on. The field then carries the input power to 7 %, and the two
    energy balances agree to 2.1e-3 after 0.1 ns and 7e-4 after 0.2 ns in
    pure Kr (LIGHT's Kr/F2: 2e-3 after 0.1 ns): while the electrons are
    thermal a cell is as wide as kT, and the flux scheme's energy exchange
    holds only to some percent (2026-10-13)."""
    m = Model(Params(t_end_s=1e-9, n_save=11))
    res = m.run()
    assert abs(accounting(res)["closure"]) < 1e-12
    log = res.kinetic
    assert log[0].get("skipped") and log[0]["t"] == 0.0
    steps = log[1:]
    assert len(steps) > 50
    assert all(not e["growth_off"] and e["reached"] for e in steps)
    assert all(e["te_kin"] == pytest.approx(e["te"], rel=1e-6)
               for e in steps)
    assert all(0.9 < e["power_ratio"] < 1.1 for e in steps)
    late = [e for e in steps if e["t"] > 1e-10]
    assert late and all(abs(e["balance"]) < 2.5e-3 for e in late)
    later = [e for e in late if e["t"] > 2e-10]
    assert later and all(abs(e["balance"]) < 1e-3 for e in later)


def test_afterglow_stays_steady(monkeypatch):
    """The time step is the pulse's. After it the EEDF is solved steady,
    since a step there left most of pure krypton's afterglow updates short
    of LIGHT's mean energy (2026-10-01)."""
    stepped, state = [], {}
    step = KB.Solver.step

    def counted(self, *a, **kw):
        stepped.append(state["on"])
        return step(self, *a, **kw)
    monkeypatch.setattr(KB.Solver, "step", counted)
    m = Model(Params(pulse_s=0.2e-9, t_end_s=0.4e-9, n_save=5))
    update = m._kinetic_update

    def tagged(t, y, on):
        state["on"] = on
        return update(t, y, on)
    m._kinetic_update = tagged
    m.run()
    assert stepped and all(stepped)
    assert any(not e["on"] and not e.get("skipped") for e in m.kin_log)


@pytest.fixture(scope="module")
def fast_run():
    return Model(Params(**FAST)).run()


def test_kinetic_run_closes_and_matches_lights_mean_energy(fast_run):
    """The chemical ledger closes, every pulse EEDF has LIGHT's mean energy
    (the table's centre), and in the pulse's second half the kinetic energy
    balance and LIGHT's agree: the field's power over the input power is 1
    less LIGHT's n_e d<eps>/dt over the input power, which the steady-shape
    kinetic balance leaves out."""
    res = fast_run
    assert res.model.eedf == "kinetic" and res.model.capture == "exact"
    assert abs(accounting(res)["closure"]) < 1e-12
    on = [e for e in res.kinetic if e["on"] and not e.get("skipped")]
    assert len(on) > 20
    for e in on:
        assert e["te_kin"] == pytest.approx(e["te"], rel=1e-6)
    late = [e for e in on if e["t"] >= 0.5 * FAST["pulse_s"]]
    assert late and max(abs(e["balance"]) for e in late) < 1e-3
    assert all(not e["growth_off"] for e in late)
    assert "kinetic EEDF" in electrons_line(res)


@SLOW
def test_reference_point_converged_in_the_kinetic_settings():
    """Renewing the EEDF at half the drift (kin_step 0.05) or every 5 ns,
    or doubling the energy cells, moves no emission band by more than 5e-4
    of itself at the reference point (LIGHT's bound on its KrF bands).
    v0_4: the 147-nm band is the two Kr2* states' rows, the near-IR one
    the four 5p -> 5s branches."""
    groups = (("KR2_RAD", "KR2F_RAD"),
              tuple(f"PRAD_{lv}" for lv in ("1s5", "1s4", "1s3", "1s2")))

    def bands(**kw):
        res = Model(Params(**kw)).run()
        m, e = res.model, res.ledger["E_dep"][-1]
        hv = {x.rid: x.hv_eV for x in m.rows}
        return np.array([sum(res.events[r][-1] * hv[r] for r in g) / e
                         for g in groups])
    ref = bands()
    for kw in (dict(kin_step=0.05), dict(kin_max_dt_s=5e-9),
               dict(kin_cells=1200)):
        assert np.max(np.abs(bands(**kw) / ref - 1.0)) < 5e-4, kw
