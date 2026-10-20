"""The bridge from 1Dv0_2's chemistry to Fluid1D, Nathaniel's 1-D fluid model.

1Dv0_2 runs LIGHT v0_2's chemistry in every cell of a 1-D mesh: the basic
network with its Kr2* route and its rates (the rates known only as a band
at their defaults unless Params1D.rates sets them), the two-temperature
tail, the captured electrons' exact energies and the gas heating, all read
from this package's own modules (network.py, species.py, xsec.py,
model.py, byte-identical copies of v0_2's when 1Dv0_2 was made). It is
1Dv0's 1-D layer on v0_2's chemistry. Fluid1D supplies what a 0-D model has not
got: the mesh, drift-diffusion transport of the charged species and of the
electron energy (Scharfetter-Gummel fluxes), the space-charge field
(semi-implicit Poisson), the electrode walls, the electron transport
coefficients from the same cross sections, and its linearly implicit
chemistry step. It is imported from ``nathaniel/1DModelKrF/src`` beside
LIGHT (or from ``FLUID1D_SRC``) and never edited: that tree is another
author's (AGPL).

This module holds what does not need a running simulation:

* ``load()`` puts Fluid1D on sys.path and imports it.
* ``registry()`` builds Fluid1D's species records for the network's species.
* ``write_inputs()`` writes the files Fluid1D reads, generated from 1Dv0_2's
  own md5-pinned LXCat files and network rows: the cross sections (JSON +
  CSV) and a reaction file (TOML) with the electron-impact rows, from
  which Fluid1D builds its electron transport table.
* ``CellChemistry`` holds the network as arrays for the per-cell stage:
  stoichiometry with counter columns (heat, emission by band, energy
  returned to the electrons, captured electrons, KrF formed), the per-cell
  rates of every row, and the energy each event costs the electrons.
* ``stage_kernel()`` compiles the per-cell chemistry and electron-energy
  stage. It calls Fluid1D's linearly implicit trial step
  (``kinetics.integrator._augmented_linearly_implicit_trial``) inside an
  outer substep loop that follows Fluid1D's own
  (``_energy_coupled_substeps``), with two changes: the rates are LIGHT's,
  refreshed every substep, including the two-temperature tail rows, whose
  rate also depends on the cell's tail ratio (Fluid1D's loop reads rates of
  the mean energy alone); and they are read at the substep's starting
  energy, before the source is added, which makes a steady Te exact.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from . import xsec
from .model import AMAGAT_CM3, EPS0, KB, KB_EV, QE, TE_FLOOR_POW, Model
from .species import ELECTRON, SPECIES, e_form

PKG_DIR = Path(__file__).resolve().parent
LIGHT_DIR = PKG_DIR.parent
DEFAULT_SRC = LIGHT_DIR.parent / "nathaniel" / "1DModelKrF" / "src"

# The electron mean-energy range [eV] of Fluid1D's transport table and of
# the stage (Te = 2/3 of it). The floor sits below the thermal 3/2 kTg of a
# 300 K gas (0.039 eV) so the electrons can cool to the gas, which v0_2's
# elastic term lets them do. The rate tables themselves sit on xsec.TE, the
# nodes of v0_2's own tables, so the stage reads v0_2's interpolants.
EPS_LO, EPS_HI = 0.03, 60.0
K_FLOOR_SI = 1e-300
SUBSTEP_CHANGE = 0.05       # largest relative change per chemistry substep
MAX_SUBSTEPS = 8000         # Fluid1D's own cap

_KIND_TYPE = {"elastic": "Elastic", "effective": "Effective",
              "excitation": "Excitation", "ionization": "Ionization",
              "attachment": "Attachment"}
# rate kinds whose rate is a function of Te alone: tabulated on the
# mean-energy grid and refreshed every substep
TE_KINDS = ("xsec", "superelastic", "stepwise", "dr_sqrt", "te_pow")
# rate kinds read per cell (Tg, gas density, ion densities). A 'param' row
# (v0_2's bands) is a constant of the run, which Model._compile has already
# resolved to frac x (Params.rates or the band's default)
CONSTANT_KINDS = ("const", "inv_tau", "param")
STATIC_KINDS = CONSTANT_KINDS + ("balance", "flannery")


# ---------------------------------------------------------------------------
# Fluid1D itself
# ---------------------------------------------------------------------------
def src_dir() -> Path:
    """Fluid1D's src/ directory: FLUID1D_SRC, else the nathaniel/ copy
    beside LIGHT."""
    env = os.environ.get("FLUID1D_SRC")
    src = Path(env).expanduser() if env else DEFAULT_SRC
    if not (src / "Fluid1D" / "simulation.py").is_file():
        raise ImportError(
            f"Fluid1D not found at {src}: set FLUID1D_SRC to "
            "Kinetic_Modeling/nathaniel/1DModelKrF/src")
    return src.resolve()


def load():
    """Import Fluid1D from src_dir() and return the package. Refuses a
    Fluid1D that Python found somewhere else."""
    src = src_dir()
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    import Fluid1D
    got = Path(Fluid1D.__file__).resolve()
    if src not in got.parents:
        raise ImportError(f"Fluid1D was imported from {got}, not from {src}")
    return Fluid1D


def ion_mobility_path() -> Path:
    """Fluid1D's cited ion-mobility table (Kr+, Kr2+ and F- in Kr)."""
    return src_dir().parent / "data" / "ion_transport" / "ion_mobility.toml"


def md5_of(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Species
# ---------------------------------------------------------------------------
def fluid_name(name: str) -> str:
    """LIGHT's electron is 'e', Fluid1D's 'e-'; every other name is shared."""
    return "e-" if name == ELECTRON else name


def registry(names):
    """Fluid1D species records: the electron, then ``names`` (LIGHT's heavy
    species) in that order. Masses come from Fluid1D's atomic masses, the
    charge and the atoms from species.py. The category and the excitation
    energy are labels only (Fluid1D's diagnostics); the energy ledger is
    species.py's formation energies."""
    load()
    from particles.constants import molecular_mass
    from particles.species import Species, SpeciesCategory, SpeciesRegistry

    cat = SpeciesCategory
    out = [Species("e-", molecular_mass("e"), -1, cat.ELECTRON,
                   constituents=("e",))]
    e_f = SPECIES["F"].e_form_eV
    for s in names:
        sp = SPECIES[s]
        atoms = tuple(a for a, k in sp.atoms for _ in range(k))
        n_f = sum(k for a, k in sp.atoms if a == "F")
        exc = 0.0 if sp.charge else max(sp.e_form_eV - e_f * n_f, 0.0)
        if sp.charge:
            kind = cat.ION if len(atoms) == 1 else cat.DIMER
        elif exc > 0.0:
            kind = cat.EXCITED if len(atoms) == 1 else cat.EXCIMER
        else:
            kind = cat.NEUTRAL if len(atoms) == 1 else cat.DIMER
        out.append(Species(s, molecular_mass(*atoms), int(sp.charge), kind,
                           excitation_energy_ev=exc, constituents=atoms))
    return SpeciesRegistry.from_species(tuple(out))


# ---------------------------------------------------------------------------
# The files Fluid1D reads
# ---------------------------------------------------------------------------
def dataset_id(gas: str, p) -> str:
    """Fluid1D process id of one LXCat process of 1Dv0_2's vendored files."""
    if p.kind in ("elastic", "effective"):
        return f"light_{gas}_elastic"
    if p.kind in ("attachment", "ionization"):
        return f"light_{gas}_{p.kind}"
    return f"light_{gas}_excitation_{p.threshold:g}eV"


def _xsec_process(row):
    """The LXCat process behind an 'xsec' row, or None."""
    if row.rate[0] != "xsec":
        return None
    return xsec.find(*row.rate[1:4])


def write_cross_sections(out_dir) -> dict:
    """Every process of 1Dv0_2's two vendored LXCat files as a Fluid1D
    JSON + CSV pair (schema cross-section-metadata/1.1, sigma in cm^2 as
    xsec.py parsed it). Returns {process id: json path}."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for gas, (name, md5, origin) in xsec.FILES.items():
        for p in xsec.processes(gas):          # refuses a moved md5
            pid = dataset_id(gas, p)
            if pid in written:
                raise ValueError(f"two {gas} processes map to {pid}")
            if p.kind in ("excitation", "ionization"):
                thr = float(p.threshold)
            elif p.kind == "attachment":
                thr = 0.0
            else:
                thr = None
            reaction = (f"E + {gas} -> {p.product}" if p.product
                        else f"E + {gas} -> E + {gas}")
            meta = {
                "schema": "cross-section-metadata/1.1",
                "data_file": f"{pid}.csv",
                "data_kind": "total",
                "process_id": pid,
                "database": "LIGHT 1Dv0_2 vendored LXCat file",
                "reaction": reaction,
                "type": _KIND_TYPE[p.kind],
                "category": _KIND_TYPE[p.kind],
                "family": "Electron",
                "projectile": "e",
                "target": gas,
                "mass_ratio_m_over_M": float(p.mass_ratio) or None,
                "threshold": thr,
                "n_points": int(p.energy.size),
                "axes": {"x_quantity": "Energy",
                         "y_quantity": "Cross section"},
                "units": {"x": "eV", "y": "cm^2"},
                "comment": (f"generated by LIGHT/1Dv0_2/fluid1d.py from "
                            f"1Dv0_2/data/{name} (md5 {md5}; {origin}); "
                            f"LXCat product label {p.product!r}"),
            }
            (out_dir / f"{pid}.json").write_text(
                json.dumps(meta, indent=1) + "\n")
            lines = [f"# process_id: {pid}", f"# reaction: {reaction}",
                     f"# source: 1Dv0_2/data/{name} (md5 {md5})",
                     "energy_eV,sigma_cm2"]
            lines += [f"{e!r},{s!r}" for e, s in zip(p.energy.tolist(),
                                                     p.sigma.tolist())]
            (out_dir / f"{pid}.csv").write_text("\n".join(lines) + "\n")
            written[pid] = out_dir / f"{pid}.json"
    return written


def _toml_str(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_reactions(path, model: Model) -> Path:
    """The network's electron-impact rows as a Fluid1D reaction file.

    Fluid1D builds its electron transport table from these rows (the
    reduced field that pairs with each mean energy comes from their
    losses). The chemistry itself never runs from this file: every row of
    the network, these included, runs from CellChemistry."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = ["# Generated by LIGHT/1Dv0_2/fluid1d.py from 1Dv0_2/network.py. "
           "Do not edit.",
           "# Only the electron-impact rows. Fluid1D builds its electron "
           "transport table",
           "# from them. The whole network runs from 1Dv0_2's own rates "
           "(fluid1d.CellChemistry).", ""]
    for c in model.compiled:
        p = _xsec_process(c.row)
        if p is None:
            continue
        gas = c.row.rate[1]
        out += ["[[reaction]]",
                f"id = {_toml_str(c.row.rid)}",
                'rate_form = "electron_impact"',
                "reactants = [" + ", ".join(_toml_str(fluid_name(s))
                                             for s in c.row.reactants) + "]",
                "products = [" + ", ".join(_toml_str(fluid_name(s))
                                            for s in c.row.products) + "]",
                f"cross_section = {_toml_str(dataset_id(gas, p))}"]
        if c.route["from_e"] > 0.0:
            out.append(f"threshold_eV = {c.route['from_e']!r}")
        out += [f"reference = {_toml_str('1Dv0_2/network.py row ' + c.row.rid)}",
                ""]
    path.write_text("\n".join(out))
    return path


def write_inputs(out_dir, model: Model):
    """Cross sections and reaction file under ``out_dir``; returns
    (cross-section dir, reaction file)."""
    out_dir = Path(out_dir)
    xs_dir = out_dir / "cross_sections"
    write_cross_sections(xs_dir)
    toml = write_reactions(out_dir / "reactions" / "light_1Dv0_2.toml", model)
    return xs_dir, toml


# ---------------------------------------------------------------------------
# The two-temperature tail, vectorized (xsec.TwoTTable.__call__ on arrays)
# ---------------------------------------------------------------------------
def two_t_rates(table, te, tt):
    """k(Te, T_t) [cm^3/s] of a TwoTTable on arrays; the same arithmetic
    as TwoTTable.__call__ (tests hold them equal)."""
    _, log_zb, log_j = xsec._tail_sums()
    lg = np.log(xsec.TE)
    lte = np.log(np.clip(te, xsec.TE_MIN, xsec.TE_MAX))
    ltt = np.log(np.clip(tt, xsec.TE_MIN, xsec.TE_MAX))
    x = -xsec.e_break() / np.exp(lte)
    li = np.interp(ltt, lg, table.log_i)
    lzb = np.interp(lte, lg, log_zb)
    lj = np.interp(ltt, lg, log_j)
    a, b = lzb, x + lj
    lz = np.maximum(a, b) + np.log1p(np.exp(-np.abs(a - b)))
    return np.maximum(np.exp(x + li - lz), xsec._K_FLOOR)


def tail_ratio(table, ne_over_n, f2_over_n):
    """xsec.TailRatio.__call__ on arrays (same clamping, bilinear)."""
    x = np.log10(np.maximum(ne_over_n, 1e-300))
    x = np.clip(x, table.x[0], table.x[-1])
    f = np.clip(f2_over_n, table.f2[0], table.f2[-1])
    i = np.clip(np.searchsorted(table.x, x) - 1, 0, len(table.x) - 2)
    j = np.clip(np.searchsorted(table.f2, f) - 1, 0, len(table.f2) - 2)
    a = (x - table.x[i]) / (table.x[i + 1] - table.x[i])
    b = (f - table.f2[j]) / (table.f2[j + 1] - table.f2[j])
    r = table.r
    return ((1 - a) * (1 - b) * r[i, j] + a * (1 - b) * r[i + 1, j]
            + (1 - a) * b * r[i, j + 1] + a * b * r[i + 1, j + 1])


def _capture_energy(c, te):
    """The captured electron's energy [eV] of row c at an array of Te
    (Model._capture_energy's rules, vectorized)."""
    te = np.asarray(te, dtype=float)
    if c.cap is None:
        return 1.5 * te
    if isinstance(c.cap, xsec.EnergyTable):
        lt = np.log(np.clip(te, xsec.TE_MIN, xsec.TE_MAX))
        return np.interp(lt, np.log(xsec.TE), c.cap.values)
    if c.kind == "dr_sqrt":
        return np.where(te > 0.02, te, 1.5 * te)
    if c.kind == "te_pow":
        n_pow = c.rate[2]
        return np.where(te > TE_FLOOR_POW, (1.5 + n_pow) * te, 1.5 * te)
    raise ValueError(f"{c.row.rid}: no vectorized captured-energy rule for "
                     f"rate kind {c.kind!r}")


# ---------------------------------------------------------------------------
# The network as per-cell arrays
# ---------------------------------------------------------------------------
class CellChemistry:
    """LIGHT's compiled rows (a 0-D Model's) as arrays for the 1-D stage.

    Columns: Fluid1D's species order (the electron first, then the
    Model's species), then the counters. Densities are SI (m^-3) and rate
    constants SI (m^3/s, m^6/s, s^-1), converted from LIGHT's cgs by
    (1e-6)^(reaction order - 1). The 'elastic' row (G6) is not a reaction
    here: its energy loss is the stage's elastic source term, as in v0_2.
    """

    def __init__(self, model: Model):
        self.model = model
        p = model.p
        if p.network != "basic":
            raise ValueError("1Dv0_2 runs the basic network (network='basic')")
        if model.ne_mode != "free":
            raise ValueError("1Dv0_2 transports the electrons: ne_mode must "
                             "be 'free'")
        self.heavy = list(model.species)
        self.names = ("e-",) + tuple(self.heavy)
        self.col = {n: i for i, n in enumerate(self.names)}
        self.n_real = len(self.names)
        self.e_col = 0
        self.rows = [c for c in model.compiled if c.kind != "elastic"]
        rows = self.rows
        n_r = len(rows)

        # --- counters ---------------------------------------------------
        bands = []
        for c in rows:
            if c.route["rad"] and c.row.band not in bands:
                bands.append(c.row.band)
        self.bands = tuple(bands)
        self.capture_rows = tuple(c.row.rid for c in rows if c.route["capture"])
        counters = (["heat+", "heat-", "to_e", "krf_formed"]
                    + [f"rad:{b}" for b in bands]
                    + [f"cap:{rid}" for rid in self.capture_rows])
        self.counters = tuple(counters)
        self.ctr = {n: self.n_real + i for i, n in enumerate(counters)}
        self.n_cols = self.n_real + len(counters)

        # --- stoichiometry, reactant tables, SI factors --------------------
        nu = np.zeros((self.n_cols, n_r))
        reactants = []
        order = np.zeros(n_r, dtype=int)
        manifold = ("KrF(B)", "KrF(C)")
        for r, c in enumerate(rows):
            for i, coeff in c.stoich:
                nu[1 + i, r] += coeff
            nu[self.e_col, r] += c.row.ne_delta
            counts = {}
            for s in c.row.heavy:
                counts[self.col[s]] = counts.get(self.col[s], 0) + 1
            if c.ne_pow:
                counts[self.e_col] = counts.get(self.e_col, 0) + c.ne_pow
            reactants.append(counts)
            order[r] = sum(counts.values())
            rt = c.route
            nu[self.ctr["heat+"], r] = max(rt["heat"], 0.0)
            nu[self.ctr["heat-"], r] = max(-rt["heat"], 0.0)
            nu[self.ctr["to_e"], r] = rt["to_e"]
            net = (sum(c.row.products.count(s) for s in manifold)
                   - sum(c.row.reactants.count(s) for s in manifold))
            nu[self.ctr["krf_formed"], r] = max(net, 0)
            if rt["rad"]:
                nu[self.ctr[f"rad:{c.row.band}"], r] = rt["rad"]
            if rt["capture"]:
                nu[self.ctr[f"cap:{c.row.rid}"], r] = rt["capture"]
        self.nu = np.ascontiguousarray(nu)
        width = max(len(d) for d in reactants)
        self.rs = np.full((n_r, width), -1, dtype=np.int64)
        self.ro = np.zeros((n_r, width), dtype=np.int64)
        self.rc = np.zeros(n_r, dtype=np.int64)
        for r, d in enumerate(reactants):
            self.rc[r] = len(d)
            for a, (s, o) in enumerate(sorted(d.items())):
                self.rs[r, a] = s
                self.ro[r, a] = o
        self.si = 1e-6 ** (order - 1).astype(float)
        self.from_e = np.array([c.route["from_e"] for c in rows])
        self.capture = np.array([c.route["capture"] for c in rows],
                                dtype=float)
        self.e_form = np.array([e_form(n) if n != "e-" else 0.0
                                for n in self.names])
        self.cv = np.array([0.0] + [SPECIES[s].cv for s in self.heavy])
        self.cations = [self.col[s] for s in self.heavy
                        if SPECIES[s].charge > 0]
        self.anions = [self.col[s] for s in self.heavy
                       if SPECIES[s].charge < 0]
        self.i_kr = self.col["Kr"]
        self.i_f2 = self.col["F2"]
        self.electron_changing = nu[self.e_col] != 0.0

        # --- rate kinds -------------------------------------------------
        self.kind = np.zeros(n_r, dtype=np.int64)
        self.tail_row = np.full(n_r, -1, dtype=np.int64)
        tails = []
        for r, c in enumerate(rows):
            if c.kind in STATIC_KINDS:
                continue
            if c.kind in TE_KINDS or c.kind == "xsec_tail":
                if c.row.tg_pow:
                    raise NotImplementedError(
                        f"{c.row.rid}: a Te-dependent rate with a Tg "
                        "factor has no per-cell table here")
                if c.kind == "xsec_tail":
                    self.kind[r] = 2
                    self.tail_row[r] = len(tails)
                    tails.append(r)
                else:
                    self.kind[r] = 1
                continue
            raise NotImplementedError(f"{c.row.rid}: rate kind {c.kind!r}")
        self.tail_rows = tuple(tails)
        self._build_tables()

    # ------------------------------------------------------------ tables
    def _build_tables(self):
        """The rates of the table rows on v0_2's own Te nodes (xsec.TE), so
        the stage reads v0_2's interpolants.

        Kind 1 (a function of Te): ln k, linear in ln Te between the nodes,
        which is RateTable's rule, so an LXCat row is v0_2's table itself.
        Kind 2 (the two-temperature tail): k = exp(-E_b/Te + h), h =
        ln I(T_t) - ln(Z_b(Te) + exp(-E_b/Te) J(T_t)) on the (Te, T_t)
        nodes. ln I, ln Z_b and ln J are TwoTTable's node values, so h is
        bilinear in (ln Te, ln T_t) where they are linear, and the steep
        -E_b/Te is added back exactly."""
        m = self.model
        te = xsec.TE
        self.lte0 = math.log(xsec.TE_MIN)
        self.lte1 = math.log(xsec.TE_MAX)
        self.inv_dlte = (te.size - 1) / (self.lte1 - self.lte0)
        self.e_break = xsec.e_break() if self.tail_rows else 0.0
        n_r = len(self.rows)
        self.lk = np.zeros((n_r, te.size))
        dummy = np.zeros(m.ns)
        for r, c in enumerate(self.rows):
            if self.kind[r] != 1:
                continue
            if isinstance(c.rate, xsec.RateTable):
                lk = c.rate._logk
            else:
                lk = np.log(np.maximum([m._k(c, float(t), 300.0, dummy, 1.0)
                                        for t in te], xsec._K_FLOOR))
            self.lk[r] = lk + math.log(self.si[r])
        self.h_tail = np.zeros((max(len(self.tail_rows), 1), te.size,
                                te.size))
        if self.tail_rows:
            _, log_zb, log_j = xsec._tail_sums()
            x = -self.e_break / te
            b = x[:, None] + log_j[None, :]               # (Te, T_t)
            a = log_zb[:, None]
            lz = np.maximum(a, b) + np.log1p(np.exp(-np.abs(a - b)))
            for t, r in enumerate(self.tail_rows):
                self.h_tail[t] = (self.rows[r].rate.log_i[None, :] - lz
                                  + math.log(self.si[r]))

    # ------------------------------------------------------------ per cell
    def tail_log_ratio(self, dens):
        """ln r per cell, r = T_t/Te from v0_2's kinetic table at the cell's
        n_e/N and F2/N (N = the heavy species), as Model._k reads it."""
        if not self.tail_rows:
            return np.zeros(dens.shape[0])
        n_gas = np.maximum(dens[:, 1:self.n_real].sum(axis=1), 1.0)
        r = tail_ratio(self.model.tail, dens[:, self.e_col] / n_gas,
                       dens[:, self.i_f2] / n_gas)
        return np.log(r)

    def static_rates(self, tg, dens):
        """(cells, rows) SI rate constants of the rows whose rate needs the
        cell's Tg, gas density or ion densities (zero in the other rows,
        whose rate the stage reads from the tables). Model._k's formulas
        on arrays."""
        p = self.model.p
        tg = np.maximum(np.asarray(tg, dtype=float), 1.0)
        n_cm3 = dens[:, :self.n_real] * 1e-6
        k = np.zeros((dens.shape[0], len(self.rows)))
        for r, c in enumerate(self.rows):
            if self.kind[r]:
                continue
            if c.kind in CONSTANT_KINDS:
                kk = np.full(dens.shape[0], float(c.rate))
            elif c.kind == "balance":
                kk = c.rate[0] * np.exp(-c.rate[1] * (1.0 / tg - 1.0 / 300.0))
            elif c.kind == "flannery":
                n_am = n_cm3[:, 1:].sum(axis=1) / AMAGAT_CM3
                f3 = ((tg / 300.0) ** c.rate[2] if len(c.rate) > 2
                      else np.ones_like(tg))
                kk = (c.rate[0] * n_am * f3 / (n_am * n_am * f3 + 1.7)
                      * c.rate[1] * self.screening(n_cm3, tg))
            else:
                raise NotImplementedError(f"{c.row.rid}: rate kind {c.kind!r}")
            tg_pow = c.row.tg_pow
            if p.trimer_tg_pow is not None and c.row.rid in ("G5", "G5C"):
                tg_pow = p.trimer_tg_pow
            if tg_pow:
                kk = kk * (tg / 300.0) ** tg_pow
            k[:, r] = kk * self.si[r]
        return k

    def screening(self, n_cm3, tg):
        """Model._screening on arrays (n_e is the cell's own)."""
        ne = np.maximum(n_cm3[:, self.e_col], 1.0)
        n_an = n_cm3[:, self.anions].sum(axis=1) if self.anions else 0.0
        n_scr = np.maximum(2.0 * (ne + n_an), 1.0) * 1e6
        lam_d = np.sqrt(EPS0 * KB * tg / (n_scr * QE ** 2))
        r_cap = 2.0 * QE ** 2 / (3.0 * 4.0 * math.pi * EPS0 * KB * tg)
        return 1.0 / (1.0 + r_cap / lam_d)

    def _node(self, lt):
        x = (lt - self.lte0) * self.inv_dlte
        j = np.clip(np.floor(x).astype(np.int64), 0, xsec.TE.size - 2)
        return j, np.clip(x - j, 0.0, 1.0)

    def table_rates(self, eps, log_r):
        """(cells, rows) SI rates of the table rows at mean energies
        ``eps`` and tail ratios exp(log_r) (the stage kernel's arithmetic,
        on arrays)."""
        te = np.asarray(eps, dtype=float) / 1.5
        lte = np.clip(np.log(te), self.lte0, self.lte1)
        j, w = self._node(lte)
        k = np.zeros((te.shape[0], len(self.rows)))
        for r in range(len(self.rows)):
            if self.kind[r] == 1:
                lk = self.lk[r]
                k[:, r] = np.exp(lk[j] + w * (lk[j + 1] - lk[j]))
            elif self.kind[r] == 2:
                h = self.h_tail[self.tail_row[r]]
                ltt = np.clip(np.log(te) + log_r, self.lte0, self.lte1)
                i, v = self._node(ltt)
                lo = (1 - v) * h[j, i] + v * h[j, i + 1]
                hi = (1 - v) * h[j + 1, i] + v * h[j + 1, i + 1]
                k[:, r] = np.exp(-self.e_break / np.exp(lte)
                                 + lo + w * (hi - lo))
        return k

    def energy_cost(self, te):
        """(cells, rows) electron energy [eV] each event takes: the row's
        from_e plus, for a capture row, the captured electron's energy at
        the cell's Te (held over one stage)."""
        eth = np.tile(self.from_e, (te.shape[0], 1))
        for r, c in enumerate(self.rows):
            if c.route["capture"]:
                eth[:, r] += c.route["capture"] * _capture_energy(c, te)
        return eth

    def capture_energies(self, te):
        """{capture row id: captured energy per event [eV] per cell}."""
        return {c.row.rid: c.route["capture"] * _capture_energy(c, te)
                for c in self.rows if c.route["capture"]}

    def progress(self, dens, k):
        """(cells, rows) rates of progress [m^-3 s^-1]."""
        prog = k.copy()
        for r in range(len(self.rows)):
            for a in range(self.rc[r]):
                prog[:, r] *= dens[:, self.rs[r, a]] ** self.ro[r, a]
        return prog

    def frequencies(self, dens, k, mask=None):
        """Per-row fastest (over cells) rate of progress divided by the
        scarcest reactant [s^-1], Fluid1D's progress_frequency rule."""
        worst = 0.0
        for r in range(len(self.rows)):
            if mask is not None and not mask[r]:
                continue
            d = [dens[:, self.rs[r, a]] for a in range(self.rc[r])]
            prod = k[:, r].copy()
            for a in range(self.rc[r]):
                prod *= d[a] ** self.ro[r, a]
            scarce = np.maximum(np.min(np.stack(d), axis=0), 1.0)
            worst = max(worst, float(np.max(prod / scarce)))
        return worst

    def elastic_power(self, dens, te, tg):
        """v0_2's elastic loss [eV m^-3 s^-1]: n_e N_Kr kel(Te) (1 - Tg/Te),
        negative when the gas heats the electrons."""
        kel = self.model.kel
        lt = np.log(np.clip(te, xsec.TE_MIN, xsec.TE_MAX))
        k = np.exp(np.interp(lt, np.log(xsec.TE), kel._logk))   # eV cm^3/s
        te_safe = np.maximum(te, 1e-3)
        return (dens[:, self.e_col] * dens[:, self.i_kr] * k * 1e-6
                * (1.0 - tg * KB_EV / te_safe))

    def chemical_energy(self, dens):
        """Formation energy [eV m^-3] per cell."""
        return dens[:, :self.n_real] @ self.e_form

    def heat_capacity(self, dens):
        """sum_s n_s cv_s k_B [eV m^-3 K^-1] per cell (the electrons have
        their own energy, as in v0_2)."""
        return (dens[:, :self.n_real] @ self.cv) * KB_EV

    # ------------------------------------------------------------ the stage
    def stage(self, dens, energy, eps, tg, source, interval, lo, hi):
        """Every row of the network in every cell over ``interval``, with
        the electron energy.

        dens: (cells, species) m^-3; energy: electron energy [eV m^-3];
        eps: mean energy [eV] (within [lo, hi]); tg: gas temperature [K];
        source: energy into the electrons [eV m^-3 s^-1], constant over
        the interval. Returns a namespace: dens (new densities), energy
        (new electron energy, within the bounds), heat [eV m^-3] (reaction
        heat plus the captured electrons' energy: what raises Tg), counters
        (per counter, per cell: extents times the per-event energies),
        adjust [eV m^-3] (energy the mean-energy bounds added; negative
        where they cut), n_sub, nu_e (fastest electron-changing frequency,
        for Fluid1D's cadence) and freqs (species, electron energy spent,
        source; the substep count's three terms).

        Per cell, to rounding: change of formation energy + change of
        electron energy + heat + radiated = source * interval + adjust.
        """
        cells = dens.shape[0]
        y = np.zeros((cells, self.n_cols))
        y[:, :self.n_real] = dens[:, :self.n_real]
        te = eps / 1.5
        base_k = self.static_rates(tg, dens)
        log_r = self.tail_log_ratio(dens)
        eth = self.energy_cost(te)

        # substeps: at most SUBSTEP_CHANGE of any reactant, of the stored
        # electron energy, or of the energy the source brings. The
        # electron-energy terms are sized on the cells that hold electrons
        # (n_e above 1e-3 of the largest) and, for the source, sit below
        # the mean-energy ceiling: a wall cell whose electrons the walls
        # drain runs to the ceiling, which the clamp enforces at any
        # substep count.
        k_now = base_k + self.table_rates(eps, log_r)
        nu_all = self.frequencies(dens, k_now)
        nu_e = self.frequencies(dens, k_now, self.electron_changing)
        spend = (eth * self.progress(dens, k_now)).sum(axis=1)
        ne = dens[:, self.e_col]
        live = ne >= 1e-3 * max(float(ne.max()), 1.0)
        heating = live & (eps < 0.5 * hi)
        stored = np.maximum(energy, 1.0)
        f_spend = (float(np.max(np.abs(spend[live]) / stored[live]))
                   if live.any() else 0.0)
        f_source = (float(np.max(np.abs(source[heating]) / stored[heating]))
                    if heating.any() else 0.0)
        freq = max(nu_all, f_spend, f_source)
        n_sub = int(min(max(1, math.ceil(interval * freq / SUBSTEP_CHANGE)),
                        MAX_SUBSTEPS))

        e = np.array(energy, dtype=float)                 # a copy
        cost = np.zeros(cells)
        adj = np.zeros(cells)
        jump = np.zeros(cells)
        lim = np.zeros(cells, dtype=np.int64)
        stage_kernel().cells(
            y, e, self.nu, np.ascontiguousarray(base_k), self.kind, self.lk,
            self.h_tail, self.tail_row, kernel_grid(self),
            np.ascontiguousarray(log_r), self.rs, self.ro, self.rc,
            interval / n_sub, n_sub, np.ascontiguousarray(eth), self.e_col,
            lo, hi, np.ascontiguousarray(source, dtype=float), cost, adj,
            jump, lim)

        d = {n: y[:, self.ctr[n]].copy() for n in self.counters}
        e = e + d["to_e"]                          # Penning, superelastic
        ne_post = np.maximum(y[:, self.e_col], 0.0)
        bounded = np.clip(e, lo * ne_post, hi * ne_post)
        adj = adj + (bounded - e)
        heat = d["heat+"] - d["heat-"]
        for rid, e_cap in self.capture_energies(te).items():
            heat = heat + e_cap * d[f"cap:{rid}"]
        return SimpleNamespace(
            dens=np.ascontiguousarray(y[:, :self.n_real]), energy=bounded,
            heat=heat, counters=d, adjust=adj, n_sub=n_sub, nu_e=float(nu_e),
            freqs=(nu_all, f_spend, f_source), limiter=int(lim.sum()))


# ---------------------------------------------------------------------------
# The per-cell stage
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def stage_kernel():
    """Compile the per-cell stage; returns a namespace with ``cells`` and
    ``rates``.

    ``cells(y, energy, nu, base_k, kind, lk, h_tail, tail_row, grid,
    log_r, rs, ro, rc, dt, n_sub, eth, e_col, lo, hi, source, cost, adj,
    jump, lim)`` advances every cell's densities ``y`` (with the counter
    columns) and electron energy [eV m^-3] by n_sub substeps of dt, in
    place, and fills the per-cell energy spent, energy added by the table
    bounds, largest mean-energy jump and limiter count. ``grid`` is
    (ln Te_min, ln Te_max, 1/d ln Te, E_break). ``rates`` is its rate
    evaluation for one cell."""
    load()
    from numba import njit, prange
    from kinetics.integrator import _augmented_linearly_implicit_trial

    trial = _augmented_linearly_implicit_trial

    @njit(cache=False)
    def node(lt, lte0, inv, n):
        x = (lt - lte0) * inv
        if x <= 0.0:
            return 0, 0.0
        if x >= n - 1:
            return n - 2, 1.0
        j = int(x)
        return j, x - j

    @njit(cache=False)
    def rates(eps, log_r, base_k, kind, lk, h_tail, tail_row, grid):
        lte0, lte1, inv, e_break = grid[0], grid[1], grid[2], grid[3]
        n = lk.shape[1]
        k = base_k.copy()
        lt_raw = math.log(eps / 1.5)
        lte = min(max(lt_raw, lte0), lte1)
        j, w = node(lte, lte0, inv, n)
        for r in range(k.shape[0]):
            kd = kind[r]
            if kd == 1:
                a0 = lk[r, j]
                k[r] = math.exp(a0 + w * (lk[r, j + 1] - a0))
            elif kd == 2:
                t = tail_row[r]
                ltt = min(max(lt_raw + log_r, lte0), lte1)
                i, v = node(ltt, lte0, inv, n)
                lo_ = (1.0 - v) * h_tail[t, j, i] + v * h_tail[t, j, i + 1]
                hi_ = (1.0 - v) * h_tail[t, j + 1, i] \
                    + v * h_tail[t, j + 1, i + 1]
                k[r] = math.exp(-e_break / math.exp(lte)
                                + lo_ + w * (hi_ - lo_))
        return k

    @njit(cache=False)
    def cell(n, energy, nu, base_k, kind, lk, h_tail, tail_row, grid, log_r,
             rs, ro, rc, dt, n_sub, eth, e_col, lo, hi, source):
        cost = 0.0
        adjust = 0.0
        jump_max = 0.0
        limiter = 0
        for _ in range(n_sub):
            # the rates come from the substep's starting energy: the
            # balance source = cost then holds at the energy the rates
            # were read at, so a steady Te is exact (read after the source
            # is added, as Fluid1D's loop does, it settles low by
            # dt * source / n_e)
            ne_pre = max(n[e_col], 0.0)
            eps = energy / ne_pre if ne_pre > 0.0 else lo
            eps = min(max(eps, lo), hi)
            energy += dt * source
            bounded = min(max(energy, lo * ne_pre), hi * ne_pre)
            adjust += bounded - energy
            if bounded != energy:
                limiter += 1
            energy = bounded
            k = rates(eps, log_r, base_k, kind, lk, h_tail, tail_row, grid)
            cand, step_cost, alpha = trial(n, energy, nu, k, rs, ro, rc, dt,
                                           eth, e_col, lo)
            if alpha <= 1.0e-14 and energy <= lo * ne_pre * (1.0 + 1.0e-14):
                # at the floor a costly row has no admissible extent;
                # Penning-type electron sources still run (Fluid1D's rule)
                k0 = k.copy()
                for r in range(k0.shape[0]):
                    if eth[r] > 0.0:
                        k0[r] = 0.0
                cand, step_cost, alpha = trial(n, energy, nu, k0, rs, ro, rc,
                                               dt, eth, e_col, lo)
            if alpha < 1.0:
                limiter += 1
            raw = max(energy - step_cost, 0.0)
            ne_post = max(cand[e_col], 0.0)
            bounded = min(max(raw, lo * ne_post), hi * ne_post)
            adjust += bounded - raw
            if bounded != raw:
                limiter += 1
            jump = abs(bounded - energy) / max(ne_post, 1.0)
            if jump > jump_max:
                jump_max = jump
            n[:] = cand
            energy = bounded
            cost += step_cost
        return energy, cost, adjust, jump_max, limiter

    @njit(parallel=True, cache=False)
    def cells(y, energy, nu, base_k, kind, lk, h_tail, tail_row, grid, log_r,
              rs, ro, rc, dt, n_sub, eth, e_col, lo, hi, source, cost, adj,
              jump, lim):
        for c in prange(y.shape[0]):
            e, co, ad, ju, li = cell(
                y[c], energy[c], nu, base_k[c], kind, lk, h_tail, tail_row,
                grid, log_r[c], rs, ro, rc, dt, n_sub, eth[c], e_col, lo, hi,
                source[c])
            energy[c] = e
            cost[c] = co
            adj[c] = ad
            jump[c] = ju
            lim[c] = li

    return SimpleNamespace(cells=cells, rates=rates)


def kernel_grid(chem: CellChemistry) -> np.ndarray:
    return np.array([chem.lte0, chem.lte1, chem.inv_dlte, chem.e_break])


def kernel_rates(chem: CellChemistry, eps, base_k, log_r):
    """The compiled kernel's own rate evaluation for each cell (tests)."""
    rates = stage_kernel().rates
    grid = kernel_grid(chem)
    return np.array([rates(float(e), float(log_r[i]),
                           np.ascontiguousarray(base_k[i]), chem.kind,
                           chem.lk, chem.h_tail, chem.tail_row, grid)
                     for i, e in enumerate(eps)])
