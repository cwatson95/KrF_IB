"""Species registry AS DATA: mass, charge, elemental atoms, formation energy.

Three jobs, all mechanical:

1. **Elemental composition + charge per species** are what make per-reaction
   mass/charge balance a computable check (defect class D2's prevention:
   "mass/charge balance is checkable mechanically per row, and no term can
   hide"). The r115-class sign error — 40 % of dArF, invisible to a parity
   harness — fails a balance check the moment the row is written.

2. **Formation energy per species** is what makes the whole-engine energy
   ledger (S3) close BY CONSTRUCTION: every reaction's energy defect is
   computed from these values, not hand-curated per row, so no reaction can
   create or destroy energy silently. E_form accuracy moves the *split*
   between the chemical-potential store and gas heat; it cannot move the
   closure, which is exact by bookkeeping.

3. **State-vector indexing** over an arbitrary number of device-grid cells.

Electron density is NOT integrated: n_e = sum(cations) - sum(anions), floored
(quasineutrality — the phantom11/spectre convention; the root Known-bad
"`ne` is not `y[0]`" records what reading it any other way cost).

Formation-energy provenance: excitation/ionization thresholds are the LXCat
values the repo already extracted (spectre engine.U_EV / ebeam.E_ION_EV,
tier REPO); dimer/excimer well depths are EST with the stated bracket — they
enter only the heat/potential split (see job 2).
"""

from __future__ import annotations

from dataclasses import dataclass

from .provenance import check_tier
from .units import AMU

# ---------------------------------------------------------------------------
# the registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpeciesDef:
    name: str
    atoms: tuple            # (("Ar", 2), ("F", 1)) — elemental composition
    charge: int             # in units of e (electron = -1 handled separately)
    mass_amu: float
    e_form_eV: float        # energy above the ground-state atoms at rest
    tier: str               # provenance of e_form_eV
    cite: str

    def __post_init__(self):
        check_tier(self.tier)


def _amu(sym_counts) -> float:
    _A = {"Ne": 20.1797, "Ar": 39.948, "Kr": 83.798, "F": 18.998403,
          "Xe": 131.293}
    return sum(_A[s] * c for s, c in sym_counts)


def _sp(name, atoms, charge, e_form, tier, cite):
    return SpeciesDef(name, tuple(atoms), charge, _amu(atoms), e_form, tier,
                      cite)


_THRESH = "LXCat threshold as extracted in spectre engine.U_EV / ebeam.py (REPO)"
_IP = "LXCat ionization potential, spectre ebeam.E_ION_EV (REPO)"
_WELL = ("EST well-depth/emission-energy estimate; enters only the ledger's "
         "heat/potential split (bracket x1/3..x3 harmless to closure)")

REGISTRY = {s.name: s for s in [
    # --- ground-state neutrals (E_form = 0 by definition; F pays half the
    #     F2 bond: D0(F2) = 1.6 eV, the value the repo's own photolysis
    #     adjudication used, spectre CHANGELOG 2026-07-15)
    _sp("Ne", [("Ne", 1)], 0, 0.0, "REPO", "ground state"),
    _sp("Ar", [("Ar", 1)], 0, 0.0, "REPO", "ground state"),
    _sp("Kr", [("Kr", 1)], 0, 0.0, "REPO", "ground state"),
    _sp("Xe", [("Xe", 1)], 0, 0.0, "REPO", "ground state"),
    _sp("F2", [("F", 2)], 0, 0.0, "REPO", "ground state"),
    _sp("F",  [("F", 1)], 0, 0.80, "REPO",
        "D0(F2)/2 = 0.80 eV (D0 = 1.6 eV; spectre CHANGELOG 2026-07-15 "
        "photolysis adjudication)"),
    # --- excited atoms (thresholds = the electron-impact energy costs)
    _sp("Ne*",  [("Ne", 1)], 0, 16.70, "REPO", _THRESH + " (Ne 3s)"),
    _sp("Ne**", [("Ne", 1)], 0, 18.70, "REPO", _THRESH + " (Ne 3p+)"),
    _sp("Ar*",  [("Ar", 1)], 0, 11.60, "REPO", _THRESH + " (Ar 4s)"),
    _sp("Ar**", [("Ar", 1)], 0, 13.20, "REPO", _THRESH + " (Ar 4p+)"),
    _sp("Kr*",  [("Kr", 1)], 0, 9.92,  "REPO", _THRESH + " (Kr 5s)"),
    _sp("Kr**", [("Kr", 1)], 0, 11.60, "REPO", _THRESH + " (Kr 5p+)"),
    _sp("Xe*",  [("Xe", 1)], 0, 8.32,  "REPO", _THRESH + " (Xe 6s)"),
    # --- excited dimers
    _sp("Ne2*", [("Ne", 2)], 0, 16.1, "EST", _WELL + "; Ne2* VUV ~80 nm"),
    _sp("Ar2*", [("Ar", 2)], 0, 10.4, "EST", _WELL + "; Ar2* VUV 126 nm"),
    _sp("Kr2*", [("Kr", 2)], 0, 9.2,  "EST", _WELL + "; Kr2* VUV 147 nm"),
    _sp("ArKr*", [("Ar", 1), ("Kr", 1)], 0, 9.5, "EST",
        _WELL + "; heteronuclear excimer anchored on Kr*(5s) 9.92 eV with "
        "a weak ~0.4 eV well (the Ar-manifold completion, Kannari 1985 / "
        "J&H 1980 Table I, first-hand 2026-08-04)"),
    # --- F2 electronic states (Chen R11/R12 thresholds)
    _sp("F2(a3)", [("F", 2)], 0, 3.16, "REPO", _THRESH + " (Chen R11)"),
    _sp("F2(A1)", [("F", 2)], 0, 4.34, "REPO", _THRESH + " (Chen R12)"),
    # --- atomic/molecular ions
    _sp("Ne+", [("Ne", 1)], +1, 21.56, "REPO", _IP),
    _sp("Ar+", [("Ar", 1)], +1, 15.76, "REPO", _IP),
    _sp("Kr+", [("Kr", 1)], +1, 14.00, "REPO", _IP),
    _sp("Xe+", [("Xe", 1)], +1, 12.13, "REPO", _IP),
    _sp("F2+", [("F", 2)], +1, 15.69, "REPO", _IP),
    _sp("Ne2+", [("Ne", 2)], +1, 20.2, "EST",
        "IP(Ne) - D0(Ne2+) ~ 1.35 eV; " + _WELL),
    _sp("Ar2+", [("Ar", 2)], +1, 14.5, "EST",
        "IP(Ar) - D0(Ar2+) ~ 1.27 eV; " + _WELL),
    _sp("Kr2+", [("Kr", 2)], +1, 12.9, "EST",
        "IP(Kr) - D0(Kr2+) ~ 1.15 eV; " + _WELL),
    _sp("NeKr+", [("Ne", 1), ("Kr", 1)], +1, 13.7, "EST",
        "IP(Kr) - weak heteronuclear bond ~0.3 eV; " + _WELL),
    _sp("ArKr+", [("Ar", 1), ("Kr", 1)], +1, 13.5, "EST",
        "IP(Kr) - heteronuclear well ~0.5 eV (between the Kr2+ 1.15 and "
        "NeKr+ ~0.3 wells; charge sits on the lower-IP Kr). The ion the "
        "Kannari/J&H ArKr+ manifold tracks (0.1.10); " + _WELL),
    _sp("F-", [("F", 1)], -1, -2.60, "REPO",
        "E_form(F) - EA(F) = 0.80 - 3.40 eV (EA: Mandl-class photodetachment "
        "data the repo's absorber table already uses)"),
    # --- excimers (upper laser levels). Reference consistency rule: the
    #     separated fragments carry E_form(F) = 0.80 eV, so
    #     E_form(excimer) = 0.80 + h*nu(emission) + lower-state repulsion —
    #     otherwise the stimulated/spontaneous emission bookkeeping would
    #     leave a spurious negative gas-heat term.
    _sp("KrF*",  [("Kr", 1), ("F", 1)], 0, 6.15, "EST",
        _WELL + "; 0.80 + 248 nm (5.00 eV) + X-wall repulsion ~0.35 eV "
        "(the LUMPED reference state — the KrF mode integrates the "
        "B/C/Bv split below; this entry remains for tooling)"),
    # --- the Kannari 1985 upper-level split (verified first-hand
    #     2026-08-04): E_B - E_C = 80 cm^-1 = 9.9e-3 eV (Julienne &
    #     Krauss via Kannari p.4314); Bv is the vibrationally hot B
    #     manifold (v > 4, ~ a few x 327 cm^-1 above; EST) ---
    _sp("KrF*(B)",  [("Kr", 1), ("F", 1)], 0, 6.15, "LIT",
        "the lasing B state; energy reference as the lumped KrF* "
        "(Kannari 1985 p.4314: E_B - E_C = 80 cm^-1)"),
    _sp("KrF*(C)",  [("Kr", 1), ("F", 1)], 0, 6.1401, "LIT",
        "the C state, 80 cm^-1 BELOW B (Kannari 1985 p.4314; radiates "
        "C->A at 1.33e7 s^-1, does not lase at 248 nm)"),
    _sp("KrF*(Bv)", [("Kr", 1), ("F", 1)], 0, 6.35, "EST",
        "vibrationally hot B manifold (v > 4), ~0.2 eV above B(v<=4) "
        "(a few x Kannari's DeltaE_v = 327 cm^-1); nascent KrF* is born "
        "here and relaxes at his 4e-11 cm3/s"),
    _sp("ArF*",  [("Ar", 1), ("F", 1)], 0, 7.52, "EST",
        _WELL + "; 0.80 + 193 nm (6.42 eV) + X-wall repulsion ~0.3 eV"),
    _sp("Kr2F*", [("Kr", 2), ("F", 1)], 0, 4.10, "EST",
        _WELL + "; 0.80 + broad ~415 nm emission (~3.0 eV) + ~0.3"),
    _sp("Ar2F*", [("Ar", 2), ("F", 1)], 0, 5.45, "EST",
        _WELL + "; 0.80 + broad ~285 nm emission (~4.35 eV) + ~0.3"),
    _sp("ArKrF*", [("Ar", 1), ("Kr", 1), ("F", 1)], 0, 4.8, "EST",
        _WELL + "; placed between the Kr2F* (4.10) and Ar2F* (5.45) "
        "triatomic wells — the Kannari/J&H displacement pair "
        "ArKrF*+Ar <-> Ar2F*+Kr straddles it, so one direction is mildly "
        "endothermic by the EST placement (heat/potential split only)"),
    _sp("NeF*",  [("Ne", 1), ("F", 1)], 0, 12.4, "EST",
        _WELL + "; predissociation-dominated (Lee 1989 [62])"),
]}

ELECTRON_NAME = "e"          # not a registry species: quasineutral, m_e known
PHOTON_PREFIX = "hv"         # photon markers in Chen rows: massless, q = 0


# Species evolved per mode (order = state layout). Same sets as spectre
# (species.py), which the network rows were written against.
KRF_MODE_SPECIES = [
    "Ne", "Ar", "Kr", "F2", "F",
    "Ne*", "Ne**", "Ne2*",
    "Ar*", "Ar**", "Ar2*",
    "Kr*", "Kr**", "Kr2*", "ArKr*",
    "F2(a3)", "F2(A1)",
    "Ne+", "Ne2+", "NeKr+", "Ar+", "Ar2+", "Kr+", "Kr2+", "ArKr+", "F2+",
    "F-",
    "ArF*", "Ar2F*", "KrF*(B)", "KrF*(C)", "KrF*(Bv)", "Kr2F*", "ArKrF*",
    "NeF*",
]

ARF_MODE_SPECIES = [
    "Ne", "Ar", "F2", "F", "Xe",
    "Ne*", "Ne**", "Ne2*",
    "Ar*", "Ar**", "Ar2*",
    "F2(a3)", "F2(A1)", "Xe*",
    "Ne+", "Ne2+", "Ar+", "Ar2+", "F2+", "F-", "Xe+",
    "ArF*", "Ar2F*",
]

DENSITY_FLOOR_M3 = 1.0


def atoms_of(name: str) -> dict:
    """Elemental composition of a species/electron/photon name."""
    if name == ELECTRON_NAME or name.startswith(PHOTON_PREFIX):
        return {}
    return dict(REGISTRY[name].atoms)


def charge_of(name: str) -> int:
    """Charge in units of e (electron −1, photon 0)."""
    if name == ELECTRON_NAME:
        return -1
    if name.startswith(PHOTON_PREFIX):
        return 0
    return REGISTRY[name].charge


def e_form_of(name: str) -> float:
    """Formation energy [eV] (electron 0 — its kinetic energy is tracked
    separately in U_e; photon 0 — its energy is the emitted quantum,
    accounted by the emitting reaction's defect)."""
    if name == ELECTRON_NAME or name.startswith(PHOTON_PREFIX):
        return 0.0
    return REGISTRY[name].e_form_eV


# ---------------------------------------------------------------------------
# state layout over C device-grid cells
# ---------------------------------------------------------------------------
class StateLayout:
    """Indexing for Y = [cells x (species..., U_e, T_g)] ++ [n_ph, ledger...].

    Per cell: species number densities [m^-3] (linear, floored on unpack),
    U_e = electron energy density [J/m^3], T_g = gas temperature [K].
    Global: n_ph = cavity photon density [m^-3] (one axial mode shared by all
    cells — the standard multi-zone laser treatment), then the ledger's named
    channels [J] appended by the engine (see ledger.py).

    n_e is derived per cell from quasineutrality, floored at
    DENSITY_FLOOR_M3.
    """

    def __init__(self, species_names, n_cells: int, n_ledger: int = 0,
                 n_aux: int = 0):
        for s in species_names:
            if s not in REGISTRY:
                raise KeyError(f"unknown species {s!r}")
        self.species = list(species_names)
        self.index = {s: i for i, s in enumerate(self.species)}
        self.n_species = len(self.species)
        self.n_cells = int(n_cells)
        self.per_cell = self.n_species + 2          # + U_e + T_g
        self.i_ue = self.n_species                  # per-cell offset
        self.i_tg = self.n_species + 1
        self.i_ph = self.per_cell * self.n_cells    # global photon index
        self.i_ledger0 = self.i_ph + 1
        self.n_ledger = int(n_ledger)
        # auxiliary GLOBAL states appended after the ledger (the circuit
        # drive's q0/qp/i1/i2/E ride here; 0 for every other mode)
        self.n_aux = int(n_aux)
        self.i_aux0 = self.i_ledger0 + self.n_ledger
        self.size = self.i_aux0 + self.n_aux
        self.positive = [s for s in self.species if s.endswith("+")]
        self.negative = [s for s in self.species if s.endswith("-")]
        self._pos_idx = [self.index[s] for s in self.positive]
        self._neg_idx = [self.index[s] for s in self.negative]

    # --- views (no copies): the engine's RHS works on these -----------------
    def cells_view(self, y):
        """(per_cell, n_cells) view of the cell block of state vector ``y``."""
        return y[: self.per_cell * self.n_cells].reshape(self.n_cells,
                                                         self.per_cell).T

    def density_block(self, y):
        """(n_species, n_cells) densities view."""
        return self.cells_view(y)[: self.n_species]

    def electron_density(self, dens_block):
        """(n_cells,) n_e from quasineutrality, floored."""
        import numpy as np
        ne = dens_block[self._pos_idx].sum(axis=0)
        if self._neg_idx:
            ne = ne - dens_block[self._neg_idx].sum(axis=0)
        return np.maximum(ne, DENSITY_FLOOR_M3)
