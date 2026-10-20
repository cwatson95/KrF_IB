"""v0_4's kinetic EEDF: the solver (boltzmann.py) and the run with it in
the loop (eedf='kinetic', the basic network's default).

The solver is checked against ML's dense solver at the same inputs, against
its own particle and energy balances, in the undriven limit (the Maxwellian
at the gas temperature), and at a state of the heating transient where the
problem at a fixed field has two solutions (Kr's Ramsauer minimum with
strong e-e) but the problem at a fixed mean energy has one. The run is
checked for the ledger's closure, for the agreement of its kinetic energy
balance with LIGHT's, and, under v0_2's electron closure, for v0_2's own
numbers.

The reference point's convergence in the update settings and the grid is
slow (minutes): set LIGHT_KINETIC_SLOW=1."""

import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from v0_5 import boltzmann as KB
from v0_5 import xsec
from v0_5.model import Model, Params
from v0_5.network import build
from v0_5.report import accounting, electrons_line, row_dict

ROOT = Path(__file__).resolve().parents[3]          # Kinetic_Modeling
ML_BOLTZ = ROOT / "ML" / "ml" / "boltzmann"
LOKI_INPUT = ROOT / "LoKI-MC" / "Code" / "Input"
SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_KINETIC_SLOW"),
                          reason="minutes (reference-point runs): set "
                          "LIGHT_KINETIC_SLOW=1")
FAST = dict(pulse_s=0.3e-6, t_end_s=0.6e-6)


def _procs():
    out = []
    for r in build((), "basic"):
        if r.rate[0] == "xsec":
            p = xsec.find(*r.rate[1:4])
            out.append(KB.Process(r.rate[1], p.kind, p.threshold, p,
                                  (r.rid,)))
    return out


def _rec(n_pow=-0.55, k0=1.6e-6, te0=0.025851999786):
    """The R11p template's constant, as model.py makes it [m^3/s eV^-n]."""
    return (k0 * te0 ** (-n_pow) * math.gamma(1.5)
            / math.gamma(1.5 + n_pow) * 1e-6)


def _laser(pen):
    """A laser state (0.42 % F2 in Kr, 450 K, n_e/N = 3e-6) with Penning
    electrons at pen per electron per second."""
    return KB.State(x={"Kr": 0.9958, "F2": 0.0042}, tg=450.0, n_col=2.45e25,
                    ne=7.4e19, pen=((pen, 5.834),), rec=_rec() * 7.2e19,
                    rec_pow=-0.55)


# the reference point 0.12 ns into the pulse: Te 0.2 eV and rising, n_e/N
# 3.4e-5, thermal attachment to F2 at 1.3e9 per second
TRANSIENT = KB.State(x={"Kr": 0.995, "F2": 0.005}, tg=300.0, n_col=2.446e25,
                     ne=8.356e20, pen=((0.221, 5.834),), rec=3.054e6,
                     rec_pow=-0.55)


def _energy(s, st, sol, procs):
    """(gain, loss) per electron [eV/s] of a solution: field and Penning
    injection against elastic, inelastic (thresholds), the captured
    energies of attachment and recombination, and nu <eps>."""
    g, F, N = sol.grid, sol.F, st.n_col
    ion, att, rec = s._frequencies(g, st, F)
    inel = sum(st.x[q.gas] * N * g.rate(q.p, F) * q.threshold for q in procs
               if q.kind in ("excitation", "ionization"))
    el = N * sum(st.x[gas] * s.elastic_power(st, sol, gas)
                 for gas in s.elastic)
    p_att = next(q.p for q in procs if q.kind == "attachment")
    w = g.ec * KB.sigma(p_att, g.ec) * F
    wr = g.ec ** (st.rec_pow + 0.5) * F
    gain = s.field_power(st, sol) + st.source_power()
    loss = (inel + el + att * float((g.ec * w).sum() / w.sum())
            + rec * float((g.ec * wr).sum() / wr.sum()) + sol.nu * sol.mean)
    return gain, loss, ion + sum(v for v, _ in st.pen) - att - rec


# ------------------------------------------------------------- the solver
@pytest.mark.skipif(not (ML_BOLTZ.is_dir() and LOKI_INPUT.is_dir()),
                    reason="ML's solver or LoKI-MC's inputs not present")
@pytest.mark.parametrize("xf2,en,ne_n,tg", [(0.005, 10.0, 1e-5, 300.0),
                                            (0.0, 20.0, 3e-6, 300.0),
                                            (0.02, 30.0, 1e-4, 600.0)])
def test_solver_matches_ml_dense_solver(xf2, en, ne_n, tg):
    """Every process of both files, a DC field, e-e and the growth term:
    with ML's e-e quadrature the sparse solver is ML's (ee.solve_field,
    numpy.linalg.solve on the whole matrix) to round-off, and v0_4's own
    quadrature moves F by 1e-5 at most at these driven states."""
    sys.path.insert(0, str(ML_BOLTZ))
    import ee
    import gases_twin as GT
    GT.load_tables(str(LOKI_INPUT))
    n_col, N, emax = 2.4e25, 600, 40.0
    kr, f2 = GT.build("Kr", tg, str(LOKI_INPUT)), GT.build("F2", tg,
                                                           str(LOKI_INPUT))
    mix = [(kr, 1.0 - xf2)] + ([(f2, xf2)] if xf2 > 0 else [])
    ml = ee.solve_field(mix, en, tg, N, emax, ne_N=ne_n,
                        ne_m3=ne_n * n_col, tol=1e-11)
    procs = [KB.Process(g, p.kind, p.threshold, p)
             for g in ("Kr", "F2") for p in xsec.processes(g)
             if p.kind in ("excitation", "ionization", "attachment")]
    st = KB.State(x={"Kr": 1.0 - xf2, "F2": xf2}, tg=tg, n_col=n_col,
                  ne=ne_n * n_col)
    mine = KB.Solver(procs, n_cells=N, emax=emax, tol=1e-11, maxit=4000,
                     ee_quadrature="cell").solve(st, en)
    assert np.max(np.abs(mine.F - ml["F"])) / ml["F"].max() < 1e-7
    m_ml = float(np.sum(ml["ec"] ** 1.5 * ml["F"]) * ml["du"])
    assert mine.mean == pytest.approx(m_ml, rel=1e-9, abs=0.0)
    assert mine.nu / n_col == pytest.approx(ml["nubar"], rel=1e-7, abs=0.0)
    for gid, gas, thr in (("Kr", kr, 9.915), ("Kr", kr, 13.996),
                          ("F2", f2, 0.0), ("F2", f2, 3.16)):
        if gid == "F2" and xf2 == 0:
            continue
        c = next(c for c in gas["colls"] if abs(c["thr"] - thr) < 1e-6)
        p = next(q.p for q in procs if q.gas == gid
                 and abs(q.threshold - thr) < 1e-6)
        assert mine.grid.rate(p, mine.F) == pytest.approx(ee.rate(c, ml),
                                                          rel=1e-6, abs=0.0)
    own = KB.Solver(procs, n_cells=N, emax=emax, tol=1e-11,
                    maxit=4000).solve(st, en)
    assert np.max(np.abs(own.F - ml["F"])) / ml["F"].max() < 1e-4


@pytest.mark.parametrize("cells", [600, 2400])
def test_undriven_eedf_is_the_maxwellian_at_the_gas_temperature(cells):
    """No field, no sources, no attachment, n_e/N = 1e-5: elastic
    collisions pull toward Tg and e-e collisions leave any Maxwellian
    alone, so the answer is the grid's Maxwellian at Tg to round-off. With
    F constant in each cell in the e-e integrals (ML's quadrature) the
    lowest boundaries' error outweighs the elastic coupling and the EEDF
    condenses into the first cell (its mean falls with the cell width)."""
    st = KB.State(x={"Kr": 1.0, "F2": 0.0}, tg=300.0, n_col=2.45e25,
                  ne=2.45e20)
    s = KB.Solver(_procs(), n_cells=cells, emax=12.0)
    sol = s.solve_joint(st, None, None, None, None, True)
    FM = sol.grid.maxwellian(KB.KB * 300.0 / KB.QE)
    bulk = sol.grid.ec < 0.5
    assert np.max(np.abs(sol.F[bulk] / FM[bulk] - 1.0)) < 1e-10
    assert sol.nu == 0.0
    # and it exchanges no energy with the gas, as the solver exchanges it
    # (central differences gave it a loss of 4.7 % of the cooling term on
    # 600 cells, and LIGHT's afterglow electrons fell to 0.93 Tg)
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
    the Penning electrons bring equals what the collisions, the captures
    and the growth term take. Without the growth term (Penning heating
    alone makes the growth solution hotter than the target) the net loss
    of electrons comes back in the first cell, and with that injection
    the balance closes as well. The solver is LIGHT's: Kr's elastic loss
    alone, as the network books it (G6), and F2's elastic collisions in
    its momentum transfer only."""
    procs = _procs()
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    assert Model(Params()).kin.elastic == ("Kr",)
    for st, target, tol in ((_laser(1.9e7), 4.5, 2e-5),
                            (_laser(2.1e7), None, 5e-4)):
        sol = s.eedf(st, target)
        assert not sol.growth_off
        gain, loss, nu = _energy(s, st, sol, procs)
        assert sol.nu == pytest.approx(nu, rel=1e-9, abs=0.0)
        assert abs(loss - gain) / gain < tol
    st = _laser(6.4e7)
    sol = s.eedf(st, 2.7)
    assert sol.growth_off and sol.mean == pytest.approx(2.7, rel=1e-9)
    gain, loss, net = _energy(s, st, sol, procs)
    injected = -net * sol.grid.ec[0]
    # first order in the cell width: 5.9e-4, 2.9e-4 and 6.6e-5 of the gain
    # on cells of 0.048, 0.024 and 0.006 eV (2026-09-30)
    assert abs(loss - gain - injected) / gain < 1e-3


def test_powers_are_the_solvers_own_flux_split():
    """The Scharfetter-Gummel flux through each inner boundary is the
    central W F_b - D dF/deps with its diffusive part times c(z) = (z/2)
    coth(z/2), z = W du/D from all of W and D (elastic, e-e and field). The
    elastic and field powers are its parts, the energy the solver itself
    exchanges. At a thermal state with a weak field c(z) reaches 1.05 on
    600 cells, by which central differences miss."""
    s = KB.Solver(_procs(), n_cells=600, emax=12.0, elastic=("Kr",))
    st = KB.State(x={"Kr": 1.0, "F2": 0.0}, tg=300.0, n_col=2.446e25,
                  ne=2.4e19)
    cold = s.solve_joint(st, None, None, None, 0.0, True)
    sol = s.eedf(st, 1.2 * cold.mean)
    assert sol.reached and sol.en_td > 0.0
    g, F = sol.grid, sol.F
    _, W, D0, dfield, hf = s._parts(g, st, sol.nu, F)      # v0_4: 5 parts
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
    less what the collisions, captures and population change take. (With
    e-e at thermal energies, on cells as wide as kT, the scheme's energy
    exchange holds only to 10 %, and to 0.8 % on four times as many.)"""
    procs = _procs()
    s = KB.Solver(procs, n_cells=600, emax=30.0, elastic=("Kr",))
    st = _laser(1.9e7)
    steady = s.solve_joint(st, 4.5, None, None, None, True)
    big = s.step(st, steady, 1e3, target_mean=4.5, adapt=False)
    assert np.max(np.abs(big.F - steady.F)) < 1e-10 * steady.F.max()
    assert big.nu == pytest.approx(steady.nu, rel=1e-8, abs=0.0)
    tiny = s.step(st, steady, 1e-20, adapt=False)
    assert np.max(np.abs(tiny.F - steady.F)) < 1e-8 * steady.F.max()
    # a target colder than 0.1 ns can reach is not reached, and the step
    # with no field stands, converged (the field's floor exit leaves F 5e-7
    # off it, and at a picosecond's step only 3e-12)
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
    mean energy (a step reads the mean it starts from, and a remap that
    moved it would heat or cool the step), and the grid's Maxwellian at a
    mean has that mean on the grid."""
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


def test_fixed_field_is_bistable_the_mean_match_is_not():
    """0.12 ns into the pulse. At 0.4 Td the problem without the growth
    term has a hot and a cold solution (Kr's Ramsauer minimum under strong
    e-e, the thermal runaway), so no bracket on E/N can find a mean energy
    there. At a fixed mean energy the solution is one, from a hot, a cold
    or a Druyvesteyn start, and the field it needs is not monotonic in the
    mean (0.2 eV needs less than 0.06 eV or 0.5 eV): the middle branch of
    the S curve."""
    s = KB.Solver(_procs(), n_cells=600, emax=12.0, elastic=("Kr",))
    g = s.grid

    def fixed(F0):
        sol = s.solve_joint(TRANSIENT, None, F0, None, None, False,
                            en_fixed=0.4)
        return sol or s.solve_joint(TRANSIENT, None, F0, None, None, False,
                                    en_fixed=0.4, m=0, maxit=3000)
    hot, cold = fixed(g.maxwellian(2.0)), fixed(g.maxwellian(0.026))
    assert hot.mean > 4.0 * cold.mean
    en = {}
    for target in (0.06, 0.2, 0.5):
        a = target * math.gamma(0.75) / math.gamma(1.25)
        starts = (g.maxwellian(2 / 3 * target * math.exp(0.25)),
                  g.maxwellian(2 / 3 * target * math.exp(-0.25)),
                  g.normalise(np.exp(-(g.ec / a) ** 2)))
        sols = [s.solve_joint(TRANSIENT, target, F0, None, None, False)
                for F0 in starts]
        for sol in sols:
            assert sol.mean == pytest.approx(target, rel=1e-9)
            assert np.max(np.abs(sol.F - sols[0].F)) / sols[0].F.max() < 1e-8
        en[target] = sols[0].en_td
    assert en[0.2] < en[0.06] and en[0.2] < en[0.5]


def test_growth_term_kept_only_where_it_gives_the_mean_energy():
    """In the pulse the growth term is kept and the field gives LIGHT's
    mean energy. 0.12 ns into the pulse the growth term moves more energy,
    replacing the electrons F2 attaches, than the field brings, so the
    EEDF is the one without it, at LIGHT's mean energy all the same. Where
    the Penning electrons alone make the growth solution hotter than the
    target, the field falls to its floor and the solve stops within a few
    steps, not reached."""
    s = KB.Solver(_procs(), n_cells=600, emax=30.0, elastic=("Kr",))
    sol = s.eedf(_laser(1.9e7), 4.5)
    assert not sol.growth_off and sol.reached and sol.en_td > 0.0
    assert sol.mean == pytest.approx(4.5, rel=1e-9)
    sol = s.eedf(TRANSIENT, 0.29)
    assert sol.growth_off and sol.mean == pytest.approx(0.29, rel=1e-9)
    sol = s.solve_joint(_laser(6.4e7), 2.7, None, None, None, True)
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
    beam = KB.State(x=_laser(0.0).x, tg=450.0, n_col=2.45e25, ne=7.4e19,
                    rec=_rec() * 7.2e19, rec_pow=-0.55,
                    inj=((2.1e7 * 5.834 / 4.506, 9.915, 4.506),))
    for st in (_laser(2.1e7), beam):
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
    st = KB.State(x={"Kr": 1.0, "F2": 0.0}, tg=300.0, n_col=2.446e25,
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


def test_beam_seed_grows_on_its_secondaries():
    """Tahoma at t = 0 under the beam: a 1e10 cm^-3 seed at 300 K with the
    secondaries injected at 7.6e12 per electron per second. Without the
    growth term the first cell would have to take up a net gain, so the
    growth term stays, its frequency is the injection's, and the EEDF is
    the injected one, hotter than LIGHT's electrons at the gas temperature
    (not reached: no field can cool it)."""
    w_eff = (0.9989 * 13.996 + (0.5119 + 0.2634) * 9.915
             + 0.5119 * (11.6 - 9.915) + 4.506)
    pairs = 3.3e5 / 1.602176634e-19 / w_eff
    ne = 1e10
    n_col = 665 / 760 * 101325 / (1.380649e-23 * 300.0)
    st = KB.State(x={"Kr": 0.99579, "F2": 0.00421}, tg=300.0, n_col=n_col,
                  ne=ne * 1e6, inj=((pairs / ne, 9.915, 4.506),))
    s = KB.Solver(_procs(), n_cells=600, emax=30.0, elastic=("Kr",))
    sol = s.eedf(st, 1.5 * 0.025852)
    assert not sol.growth_off and not sol.reached and sol.en_td == 0.0
    assert sol.nu == pytest.approx(pairs / ne, rel=1e-3)
    assert sol.mean == pytest.approx(4.506, rel=1e-2)


def test_grid_is_sized_from_the_tail():
    """From a grid far too long and one far too short, the same grid, with
    F falling 11 to 25 decades across it, and the same rates (to the
    solver's tolerance: the two reach it from different starts). The fall is
    slow in the bulk and steep past the Kr thresholds, so a grid scaled by
    the decades oscillated (41 and 13 eV at the reference point, and the
    13 eV one, 3 decades, was returned). The solver's grid is always the
    returned solution's."""
    out = []
    for emax in (60.0, 13.2):
        s = KB.Solver(_procs(), n_cells=600, emax=emax,
                      elastic=("Kr",))
        sol = s.eedf(_laser(1.9e7), 4.5)
        assert sol.grid is s.grid
        assert 11.0 <= KB._decades(sol.F) <= 25.0
        out.append((sol.grid.emax, [sol.grid.rate(q.p, sol.F)
                                    for q in _procs()]))
    assert out[0][0] == out[1][0]
    assert out[0][1] == pytest.approx(out[1][1], rel=1e-7, abs=0.0)


# ------------------------------------------------------------ the run
def _beam_model():
    """The prototype beam drive (LIGHT/docs/survey_2026-09-27/scripts/
    ebeam_screen.py) on v0_4, loaded as ML loads it, its two v0 imports
    pointed at v0_4."""
    path = (ROOT / "LIGHT" / "docs" / "survey_2026-09-27" / "scripts"
            / "ebeam_screen.py")
    if not path.is_file():
        pytest.skip("the beam-drive script is not beside LIGHT")
    src = path.read_text()
    assert src.count("from v0.") == 2
    ns = {"__file__": str(path), "__name__": "ebeam_screen_v0_4_test"}
    exec(compile(src.replace("from v0.", "from v0_5."), str(path), "exec"),
         ns)
    return ns


def test_beam_drive_feeds_its_secondaries_to_the_kinetic_eedf():
    """Under the beam drive the kinetic state carries the secondaries (one
    per ion pair, the drive's own numbers), the ledger closes, once the
    seed has grown LIGHT's mean energy is the kinetic one, and the kinetic
    energy balance is LIGHT's with the beam's power injected, not carried
    by the field (6.5e-4 of it after 10 ns, 2026-09-30). v0_4: with the
    5s <-> 5p pair in the kinetic equation the no-field mean sits up to
    1.9e-4 above LIGHT's Te between 10 and 16 ns, while the Kr* and Kr**
    populations the pair acts on still grow fast between renewals, and
    agrees to 1e-9 from 17 ns on (balance 1.2e-3 falling to 5.4e-4)."""
    E = _beam_model()
    p = Params(pressure_atm=665 / 760, x_f2=2.8 / 665, power_W_cm3=3.3e5,
               pulse_s=30e-9, t_end_s=30e-9, ne_mode="free", ne_cm3=1e10,
               n_save=301)
    m = E["BeamModel"](p)
    inj = m.kinetic_injection(0.0, m.initial_state(), True)
    assert len(inj) == 1
    rate, e_top, mean = inj[0]
    assert rate == pytest.approx(p.power_W_cm3 / E["QE"] / E["W_EFF"],
                                 rel=1e-12)
    assert (e_top, mean) == (E["KRS"], E["EPS_SEC"])
    assert m.kinetic_injection(0.0, m.initial_state(), False) == ()
    s, _, _ = m._kinetic_state(m.initial_state(), 0.0, True)
    assert s.inj and s.inj[0][0] == pytest.approx(rate / 1e10, rel=1e-12)
    res = m.run()
    assert abs(accounting(res)["closure"]) < 1e-12
    late = [e for e in res.kinetic if e["on"] and not e.get("skipped")
            and e["t"] > 10e-9]
    assert late and all(e["te_kin"] == pytest.approx(e["te"], rel=1e-3)
                        for e in late)
    settled = [e for e in late if e["t"] > 20e-9]
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
    which scipy locates an event. In pure Kr the event comes back at t = 0
    itself, and the run renewed its table there 12,405 times without moving
    (2026-09-30). It now steps past."""
    E = _beam_model()
    p = Params(pressure_atm=665 / 760, x_f2=0.0, power_W_cm3=3.0e5,
               pulse_s=3e-15, t_end_s=3e-15, ne_mode="free", ne_cm3=1e10,
               groups=("electron",), n_save=3)
    m = E["BeamModel"](p)
    times = _counting(m, 60)
    res = m.run()
    assert abs(accounting(res)["closure"]) < 1e-12
    assert times[-1][0] > 0.0


def test_rows_dropped_after_construction_leave_the_kinetic_equation():
    """ML drops rows from a model it has built (run_beam's drop, its fits)
    and rescales others. A run's kinetic equation carries the rows the model
    holds when it runs, and the table's Maxwellian values at its ends belong
    to their own rows, though a dropped row (C1a, the list's seventh) moves
    every later row up the list."""
    m = Model(Params(pulse_s=1e-14, t_end_s=1e-14, n_save=3))
    m.compiled = [c for c in m.compiled if c.row.rid != "C1a"]
    times = _counting(m, 20)
    m.run()
    assert "C1a" not in {r for pr in m.kin.procs for r in pr.rows}
    tab, y = m._tab, times[-1][1]
    n = np.maximum(y[:m.ns], 0.0) * m.N0
    ne = m.electron_density(n)
    byj = {c.j: c for c in m.compiled}
    ends = np.exp(tab["lnte"][[0, -1]])
    assert len(tab["lnkM"]) >= 10
    for j, lnk in tab["lnkM"].items():
        want = [math.log(max(m._k(byj[j], e, y[m.i_tg], n, ne), 1e-300))
                for e in ends]
        assert list(lnk) == pytest.approx(want, rel=1e-12, abs=0.0)


def test_table_is_carried_past_its_ends_by_lights_maxwellian():
    """Past the kinetic table's ends LIGHT carries it with its own
    Maxwellian dependence on Te: a rate by its logarithmic slope, a
    captured energy by its ratio and the elastic power by its difference,
    which passes through 0 at Tg. Each meets the table at its end. Held flat
    there, the table left LIGHT's energy equation no loss that falls with
    its Te (1.3 against 1.6 eV at Tahoma, 2026-09-30)."""
    from v0_5.model import KB_EV
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
    assert len(rows) >= 10 and tab["capM"]
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
    term on. The field then carries the input power to 10 %, where the
    steady solutions carried 2 to 18 % of it and dropped the growth term at
    every update (2026-10-01). What is left over is the step's average
    d<eps>/dt against LIGHT's instantaneous one, and with it the two energy
    balances agree to 2e-3 after 0.1 ns (while the electrons are thermal a
    cell is as wide as kT, and the flux scheme's energy exchange holds
    only to some percent)."""
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
    assert late and all(abs(e["balance"]) < 2e-3 for e in late)


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


def test_two_t_closure_reproduces_v0_2():
    """Under v0_2's electron closure, with v0_4's 5p rung neutralised (the
    5s <-> 5p cycle off, the 5p relaxed to 5s at once) v0_4 is v0_2 to
    the ion-ion limb change (v0_4): same network otherwise, the emission
    columns within 0.3 % and the ion flows within 1 % (the Langevin limb
    lowers the ion-ion rate 10 to 14 %). v0_4's row_dict carries the
    optics, walls and EM columns on top of v0_2's."""
    from v0_2.model import Params as P2
    from v0_2.model import run as run2
    from v0_2.report import row_dict as row2
    off = (("k_5s5p", 0.0), ("kps_kr", 1e-10))
    a = row_dict(Model(Params(eedf="two_t", rates=off, **FAST)).run())
    b = row2(run2(P2(**FAST)))
    assert set(b) <= set(a)
    for k in b:
        if isinstance(a[k], float) and k not in ("closure",):
            tol = 3e-3 if k in ("emission_248", "to_KrFB",
                                "krf_to_248_pct") else 1e-2
            assert a[k] == pytest.approx(b[k], rel=tol, abs=1e-9), k
        elif not isinstance(a[k], float):
            assert a[k] == b[k], k


@SLOW
def test_reference_point_converged_in_the_kinetic_settings():
    """Renewing the EEDF at half the drift (kin_step 0.05) or every 5 ns,
    or doubling the energy cells, moves no band by more than 5e-4 of
    itself at the reference point."""
    def bands(**kw):
        res = Model(Params(**kw)).run()
        m, e = res.model, res.ledger["E_dep"][-1]
        return np.array([res.events[r][-1] * next(
            x for x in m.rows if x.rid == r).hv_eV / e
            for r in ("G4", "G4A", "G4C", "C3", "KR2_RAD")])
    ref = bands()
    for kw in (dict(kin_step=0.05), dict(kin_max_dt_s=5e-9),
               dict(kin_cells=1200)):
        assert np.max(np.abs(bands(**kw) / ref - 1.0)) < 5e-4, kw
