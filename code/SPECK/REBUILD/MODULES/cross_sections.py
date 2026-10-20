"""LXCat cross-section loader — the data behind the beam's Spencer–Fano route.

ADAPTED, NOT RE-DERIVED: a port of ``PIC_MCC_3D/MODULES/cross_sections.py``
(itself adapted from ARPAE/pic3d, itself from the tested pic1d @ 4e2a9fc,
oracle-matched against LoKI-MC). Parser, ``Process``/``CrossSectionSet``
semantics, the effective→elastic synthesis and the fixed log-energy lookup
are deliberately identical — the point of adapting a gated module is that
its gates still mean something. Differences, stated:

* the species argument is a plain ``{gas: number_density_m3}`` dict (this
  package's gas state lives on the Device, not in a GasMixture class);
* Ne is in the file map (the Biagi set the root CHANGELOG mandates for Ne —
  never the BSR set beyond 150 Td, never Phelps: both Known-bad);
* the data-root search is marker-based (never parent-counting: the
  PIC_MCC_3D Known-bad records how ``parents[2]`` silently broke inside a
  git worktree), env override ``REBUILD_LXCAT`` first.

CEILING CAVEAT (inherited, load-bearing): ``Process.sigma_at`` returns ZERO
above the last tabulated point — it does not extrapolate and does not raise.
Measured ceilings on this repo's files are ~1 keV; every consumer must stay
below ``energy_ceiling_eV`` (the beam module refuses to cascade above it).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .units import E_CHARGE, M_E

_HERE = Path(__file__).resolve().parent

_REPO_MARKERS = (
    Path("ARPAE") / "pic3d" / "lxcat",        # ARPAE's standalone copy
    Path("LoKI-MC") / "Code" / "Input",       # the repo's LoKI input tree
)

LXCAT_ENV = "REBUILD_LXCAT"


def _repo_data_dirs():
    found = []
    for base in [_HERE, *_HERE.parents]:
        for rel in _REPO_MARKERS:
            cand = base / rel
            if cand.is_dir() and cand not in found:
                found.append(cand)
    return found


def lxcat_roots():
    """Directories searched for LXCat files, in order (first hit wins)."""
    import os
    roots = []
    env = os.environ.get(LXCAT_ENV)
    if env:
        roots.append(Path(env).expanduser())
    roots.append(_HERE / "lxcat")
    roots.extend(_repo_data_dirs())
    return roots


# gas name -> (flat filename, nested repo-tree relative path)
GAS_FILES = {
    "Ar": ("Ar_LXCat.txt", "Argon/Ar_LXCat.txt"),
    "Kr": ("Kr_LXCat.txt", "Krypton/Kr_LXCat.txt"),
    "F2": ("F2_LXCat.txt", "Fluorine2/F2_LXCat.txt"),
    "He": ("He_LXCat.txt", "Helium/He_LXCat.txt"),
    "Ne": ("Ne_LXCat.txt", "Neon/Ne_LXCat.txt"),
}

_KINDS = {"ELASTIC", "EFFECTIVE", "EXCITATION", "IONIZATION", "ATTACHMENT"}


def gas_file(name: str) -> Path:
    if name not in GAS_FILES:
        raise KeyError(f"no LXCat file mapped for gas {name!r}; "
                       f"have {sorted(GAS_FILES)}")
    flat, nested = GAS_FILES[name]
    tried = []
    for root in lxcat_roots():
        for rel in (flat, nested):
            p = root / rel
            tried.append(p)
            if p.is_file():
                return p
    raise FileNotFoundError(
        f"LXCat data for {name!r} not found. Tried:\n  "
        + "\n  ".join(str(t) for t in tried))


@dataclass
class Process:
    kind: str
    target: str
    product: str
    threshold: float          # energy loss [eV]
    mass_ratio: float         # m_e/M (elastic/effective only)
    energy: np.ndarray
    sigma: np.ndarray

    def sigma_at(self, eps):
        return np.interp(eps, self.energy, self.sigma, left=0.0, right=0.0)


def _is_floats(line: str) -> bool:
    parts = line.split()
    if not parts:
        return False
    try:
        for p in parts:
            float(p)
        return True
    except ValueError:
        return False


def parse_lxcat(path):
    lines = Path(path).read_text().splitlines()
    procs = []
    i, n = 0, len(lines)
    while i < n:
        kw = lines[i].strip()
        if kw not in _KINDS:
            i += 1
            continue
        kind = kw.lower()
        target_line = lines[i + 1].strip()
        target, _, product = target_line.replace("<->", "->").partition("->")
        target, product = target.strip(), product.strip()

        j = i + 2
        threshold = 0.0
        mass_ratio = 0.0
        if kind != "attachment" and _is_floats(lines[j]):
            vals = [float(v) for v in lines[j].split()]
            if kind in ("elastic", "effective"):
                mass_ratio = vals[0]
            else:
                threshold = vals[0]
            j += 1

        while j < n and not lines[j].lstrip().startswith("---"):
            j += 1
        j += 1
        e_list, s_list = [], []
        while j < n and not lines[j].lstrip().startswith("---"):
            if _is_floats(lines[j]):
                a, b = lines[j].split()[:2]
                e_list.append(float(a))
                s_list.append(float(b))
            j += 1

        if e_list:
            procs.append(Process(kind=kind, target=target, product=product,
                                 threshold=threshold, mass_ratio=mass_ratio,
                                 energy=np.asarray(e_list, dtype=float),
                                 sigma=np.asarray(s_list, dtype=float)))
        i = j + 1
    return procs


class CrossSectionSet:
    """Collision processes for a mixture, with per-process target densities.

    ``densities`` maps gas name -> number density [m^-3]. Lumped 'effective'
    sets are dropped; a file with only an effective set gets a synthesised
    elastic = effective − sum(inelastic), clipped at zero (the pic1d remedy,
    inherited with its provenance)."""

    def __init__(self, densities: dict):
        self.processes = []
        self.density = []
        self.species_density = {}
        self.target_mass = {}
        for name, dens in densities.items():
            if dens <= 0.0:
                continue
            self.species_density[name] = float(dens)
            procs = parse_lxcat(gas_file(name))
            if not any(p.kind == "elastic" for p in procs):
                eff = next((p for p in procs if p.kind == "effective"), None)
                if eff is not None:
                    sig_inel = np.zeros_like(eff.sigma)
                    for p in procs:
                        if p.kind in ("excitation", "ionization", "attachment"):
                            sig_inel += p.sigma_at(eff.energy)
                    procs.append(Process(
                        kind="elastic", target=eff.target, product="",
                        threshold=0.0, mass_ratio=eff.mass_ratio,
                        energy=eff.energy.copy(),
                        sigma=np.maximum(eff.sigma - sig_inel, 0.0)))
            for proc in procs:
                if proc.kind == "effective":
                    continue
                # remember which GAS the process belongs to (the per-species
                # beam-yield split needs it; LXCat 'target' strings vary)
                proc.gas = name
                self.processes.append(proc)
                self.density.append(float(dens))
                if proc.kind == "elastic" and proc.mass_ratio > 0:
                    self.target_mass[name] = M_E / proc.mass_ratio
        if not self.processes:
            raise ValueError("CrossSectionSet: no processes loaded")
        self.density = np.asarray(self.density, dtype=float)
        self.n_proc = len(self.processes)
        self.kinds = [p.kind for p in self.processes]
        self.thresholds = np.array([p.threshold for p in self.processes])
        self.mass_ratios = np.array([p.mass_ratio for p in self.processes])
        self.gases = [p.gas for p in self.processes]
        self._tab_eps = np.geomspace(1e-3, 2e3, 8192)
        self._tab_log0 = np.log(self._tab_eps[0])
        self._tab_dlog = np.log(self._tab_eps[1] / self._tab_eps[0])
        self._tab_sigma = np.array([p.sigma_at(self._tab_eps)
                                    for p in self.processes])
        self._tab_nu_unit = (self._tab_sigma * self.density[:, None]
                             * np.sqrt(2.0 * self._tab_eps * E_CHARGE
                                       / M_E)[None, :])

    def _tab_idx(self, eps):
        e = np.maximum(np.asarray(eps, dtype=float), self._tab_eps[0])
        idx = ((np.log(e) - self._tab_log0) / self._tab_dlog + 0.5).astype(np.intp)
        return np.clip(idx, 0, self._tab_eps.size - 1)

    def sigma_matrix(self, eps):
        return self._tab_sigma[:, self._tab_idx(eps)]

    def nu_matrix(self, eps):
        """(n_proc, n_eps) collision frequencies n_p sigma_p(eps) v(eps)."""
        return self._tab_nu_unit[:, self._tab_idx(eps)]

    @property
    def energy_ceiling_eV(self) -> float:
        """Highest energy at which the set is still COMPLETE (see module
        docstring — a safety number, not a curiosity). Defined over the
        channels a fast electron needs (elastic/ionization); excitation and
        attachment tables legitimately end where their channel closes."""
        ends = [p.energy[-1] for p in self.processes
                if p.kind in ("elastic", "ionization")]
        return float(min(ends)) if ends else 0.0

    def describe(self) -> str:
        counts = {}
        for k in self.kinds:
            counts[k] = counts.get(k, 0) + 1
        return (f"{self.n_proc} processes over {len(self.species_density)} "
                "species: " + ", ".join(f"{v} {k}"
                                        for k, v in sorted(counts.items())))
