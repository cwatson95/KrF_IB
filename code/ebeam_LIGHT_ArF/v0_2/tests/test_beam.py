"""The electron-beam drive (beam.py): its yields per ion pair, their
identity and REBUILD parity, W against the measurement, and the model's
closures under the beam."""
import sys
from pathlib import Path

import numpy as np
import pytest

from v0_2 import beam
from v0_2.model import QE, Model, Params, run
from v0_2.network import HV_AR2
from v0_2.report import accounting, reservoir
from v0_2.species import SPECIES

from .test_model import FAST, _check_closures, approx

SPECK = next((p / "SPECK" for p in Path(__file__).resolve().parents
              if (p / "SPECK" / "REBUILD" / "MODULES" / "beam.py").is_file()),
             None)


def _rebuild():
    if SPECK is None:
        pytest.skip("SPECK/REBUILD is not in this checkout")
    if str(SPECK) not in sys.path:
        sys.path.insert(0, str(SPECK))
    from REBUILD.MODULES import beam as rb
    return rb


@pytest.mark.parametrize("x", [0.0, 0.003, 0.005, 0.01])
@pytest.mark.parametrize("route, sharing", [("asymptotic", "opb"),
                                            ("integral", "equal"),
                                            ("asymptotic", "equal")])
def test_yields_close_the_identity_and_conserve_atoms(x, route, sharing):
    y = beam.yields(x, route, sharing)
    assert abs(y.identity_residual()) < 1e-12 * y.W
    assert y.eps_sec > 0.0 and y.heat >= 0.0 and y.ar_star > 0.0
    for atom in ("Ar", "F"):
        net = sum(n * dict(SPECIES[s].atoms).get(atom, 0)
                  for s, n in y.per_pair().items())
        assert net == approx(0.0, abs=1e-12), atom
    if x == 0.0:
        assert y.f2_diss == 0.0 and y.f2_ion == 0.0


def test_rebuild_parity():
    """The equal-sharing integral route is REBUILD's spencer_fano_yields on
    the same (md5-identical) files, and the stopping power is REBUILD's."""
    rb = _rebuild()
    n = 2.446e25
    for x in (0.0, 0.005):
        dens = {"Ar": (1 - x) * n} | ({"F2": x * n} if x else {})
        ref = rb.spencer_fano_yields(dens)
        mine = beam.yields(x, "integral", "equal")
        assert mine.W == approx(ref.W_eV, rel=1e-12)
        exc = sum(v for s, v in ref.per_pair.items() if s in ("Ar*", "Ar**"))
        assert mine.ar_star == approx(exc, rel=1e-12)
    for kev in (12.0, 150.0, 600.0):
        dens = beam.gas_densities_m3(1.0, 0.005)
        assert beam.stopping_power(dens, kev * 1e3) == approx(
            rb.stopping_power(dens, kev * 1e3), rel=1e-12)


def test_w_is_scored_not_set():
    """W is an output: 25.2 eV with the Opal-Peterson-Beaty secondaries,
    4.5 % under the measured 26.4 eV (ICRU 31). Pinned, so a change of the
    data or the recursion shows here; it is never moved to the anchor."""
    y = beam.yields(0.0)
    assert y.W == approx(25.20, rel=0.01)
    assert (y.W - beam.W_LIT_AR) / beam.W_LIT_AR == approx(-0.045, abs=0.01)
    assert y.ar_star == approx(0.389, rel=0.02)


def _sweep_ratios(sharing, energies):
    cas = beam.Cascade(0.0)
    grid, u = beam.sweep(cas, 900.0, sharing, 2000)
    ion = u[:, cas.is_ion].sum(axis=1)
    exc = u[:, ~cas.is_ion].sum(axis=1)
    at = lambda e, a: float(np.interp(e, grid, a))   # noqa: E731
    return ([e / at(e, ion) for e in energies],
            [at(e, exc) / at(e, ion) for e in energies])


def test_equal_sharing_oscillates_and_opal_beaty_does_not():
    """Why the default is not REBUILD's rule: with equal sharing the yields
    per pair swing with the source energy (E/U_ion from 28 to 32 eV and the
    excitations per ion from 0.32 to 0.96 between 300 and 900 eV); with the
    Opal-Peterson-Beaty spectrum they settle smoothly."""
    es = (300.0, 400.0, 500.0, 700.0, 900.0)
    w_eq, x_eq = _sweep_ratios("equal", es)
    w_op, x_op = _sweep_ratios("opb", es)
    assert max(x_eq) / min(x_eq) > 2.0 and max(w_eq) / min(w_eq) > 1.1
    assert max(x_op) / min(x_op) < 1.15 and max(w_op) / min(w_op) < 1.04
    assert all(a > b for a, b in zip(w_op, w_op[1:]))   # falls to its asymptote


def test_asymptotic_route_does_not_depend_on_its_window():
    cas = beam.Cascade(0.005)
    grid, u = beam.sweep(cas, 0.9 * cas.ceiling, "opb")
    out = []
    for lo in (0.15, 0.20, 0.35):
        sel = grid >= lo * cas.ceiling
        d = np.polyfit(grid[sel], u[sel], 1)[0]
        ion = d[cas.is_ion].sum()
        exc = d[[k for k in range(cas.thr.size) if cas.gases[k] == "Ar"
                 and cas.kinds[k] == "excitation"]].sum()
        out.append((1.0 / ion, exc / ion))
    for w, x in out[1:]:
        assert w == approx(out[0][0], rel=2e-3)
        assert x == approx(out[0][1], rel=0.02)


def test_deposited_power_is_the_stopping_power_times_the_current():
    """Rokni's beam, 150 keV and 5 A/cm2 at the gas: 18.8 kW/cm3 per atm of
    Ar on a straight ray (the survey had estimated 20)."""
    p = beam.deposited_power_W_cm3(150.0, 5.0, 1.0, 0.0)
    assert p == approx(18.76e3, rel=2e-3)
    assert beam.deposited_power_W_cm3(150.0, 5.0, 2.0, 0.0) == approx(2 * p)
    assert beam.csda_range_cm(1.0, 0.0, 12e3) == approx(0.33, rel=0.05)


@pytest.mark.parametrize("kw", [
    dict(),
    dict(beam_route="integral", beam_sharing="equal"),
    dict(pressure_atm=0.5, power_W_cm3=5e5, ne_cm3=1e14),
    dict(x_f2=0.0),
    dict(groups=("all",)),
])
def test_beam_runs_close(kw):
    res = run(**FAST, drive="beam", **kw)
    _check_closures(res)
    p = res.p
    assert res.ledger["E_dep"][-1] * QE == approx(p.power_W_cm3 * p.pulse_s,
                                                  rel=1e-9)
    y = res.model.beam
    pairs = res.ledger["E_dep"][-1] / y.W
    assert res.ledger["heat_beam"][-1] == approx(pairs * y.heat, rel=1e-9)


def test_pure_argon_turns_every_pair_and_excitation_into_one_photon():
    """The network's accounting in pure Ar under the beam, the ceiling
    Morozov's 0.33 sits under: every ion and every excitation ends as one
    127-nm photon (the bulk electrons' own excitations add about 1 %)."""
    p = Params(pressure_atm=1000 / 1013.25, x_f2=0.0, power_W_cm3=330.0,
               pulse_s=450e-9, t_end_s=450e-9 + 20e-6, ne_cm3=1e8,
               drive="beam", n_save=401)
    res = Model(p).run()
    y = res.model.beam
    pairs = res.ledger["E_dep"][-1] / y.W
    photons = res.events["AR2_RAD"][-1]
    assert photons / (pairs * (1.0 + y.ar_star)) == approx(1.0, rel=0.03)
    eff = photons * HV_AR2 / res.ledger["E_dep"][-1]
    assert eff == approx(y.one_photon_ceiling(HV_AR2), rel=0.03)
    assert reservoir(res, "Ar*")["balance"] < 1e-6 * reservoir(res, "Ar*")[
        "formed"]


def test_power_drive_is_untouched_and_unknown_drives_are_refused():
    a = accounting(run(**FAST))
    assert "beam: to Ar+ (every pair x 15.76 eV)" not in a["flows"]
    with pytest.raises(ValueError):
        Model(Params(drive="discharge"))
    with pytest.raises(ValueError):
        beam.yields(0.005, "differential")
    with pytest.raises(ValueError):
        beam.yields(0.005, "asymptotic", "half")
