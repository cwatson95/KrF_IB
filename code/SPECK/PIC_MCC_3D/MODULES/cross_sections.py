"""LXCat cross-section loader -- the collision data behind M3.

ADAPTED, NOT RE-DERIVED.  This is a port of ``ARPAE/pic3d/cross_sections.py``
(itself an adapted copy of the tested ``pic1d/cross_sections.py`` @ 4e2a9fc,
oracle-matched against LoKI-MC).  The parser, the Process/CrossSectionSet
semantics, the effective->elastic synthesis and the fixed log-energy lookup
tables are deliberately identical, because the whole point of adapting a
gated module is that its gates still mean something.  ``PIC_MCC_3D.txt`` §7
says so explicitly: adapt ARPAE's tested modules at M3, do not re-derive.

WHAT IS DIFFERENT HERE
  * densities come from a :class:`~PIC_MCC_3D.MODULES.gas.GasMixture` rather
    than from loose (pressure, fraction) floats, so one object owns the
    neutral background;
  * the data root is resolved through a search PATH (see ``lxcat_roots``)
    instead of a single local directory -- this package ships no copy of the
    LXCat files, and duplicating ~6000 lines of tabulated data a third time
    to achieve "standalone" would be a worse trade than a documented lookup.
    Drop a ``MODULES/lxcat/`` directory in and it wins; otherwise ARPAE's
    copy is used, then the repo's LoKI-MC input tree.

LXCat block layout (see the header of any LXCat file):
    line 1: process kind  (ELASTIC | EFFECTIVE | EXCITATION | IONIZATION | ATTACHMENT)
    line 2: target name   (optionally "A -> A*" / "A <-> A*")
    line 3: m/M (elastic/effective) | threshold eV (excitation/ionization) | *absent* (attachment)
    comments ... then a dashed fence, the (energy, sigma) table, a closing fence.

PROVENANCE CAVEAT, carried from the source modules: the elastic entries in
these files are MOMENTUM-TRANSFER cross sections, and the MCC scatters
isotropically off them.  That is the self-consistent reference choice, and it
is also the open question ``pic1d/README_LIMITS.md`` §2.2 names -- the x1.6
LoKI-MC agreement shares the same basis, so it cannot rule out a common
isotropic-on-MT bias.  Do not present isotropic scattering here as validated
against anisotropic data; it is not.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .gas import GasMixture
from .units import E_CHARGE, M_E

_HERE = Path(__file__).resolve().parent

# Relative locations of known data trees, searched under whichever ancestor
# actually contains one.
_REPO_MARKERS = (
    Path("ARPAE") / "pic3d" / "lxcat",        # ARPAE's standalone copy
    Path("LoKI-MC") / "Code" / "Input",       # the repo's LoKI input tree
)

# environment override, checked first -- an explicit path always wins
LXCAT_ENV = "PIC_MCC_3D_LXCAT"


def _repo_data_dirs() -> list:
    """Find the data trees by SEARCHING UPWARD for a marker.

    Deliberately not ``_HERE.parents[2]``: counting parent directories assumes
    this package sits at a fixed depth below the repo root, which is false the
    moment the checkout is relocated -- a git worktree under
    ``.claude/worktrees/<name>/`` adds two levels and the count silently
    resolves to a directory that does not exist.  Walking up until a marker is
    found works in both, and in a standalone copy simply finds nothing (which
    is the correct answer there).
    """
    found = []
    for base in [_HERE, *_HERE.parents]:
        for rel in _REPO_MARKERS:
            cand = base / rel
            if cand.is_dir() and cand not in found:
                found.append(cand)
    return found


def lxcat_roots() -> list:
    """Directories searched for LXCat files, in order (first hit wins)."""
    import os
    roots = []
    env = os.environ.get(LXCAT_ENV)
    if env:
        roots.append(Path(env).expanduser())
    roots.append(_HERE / "lxcat")            # a local copy, if one is dropped in
    roots.extend(_repo_data_dirs())
    return roots


# gas name -> (flat filename, nested repo-tree relative path)
GAS_FILES = {
    "Ar": ("Ar_LXCat.txt", "Argon/Ar_LXCat.txt"),
    "Kr": ("Kr_LXCat.txt", "Krypton/Kr_LXCat.txt"),
    "F2": ("F2_LXCat.txt", "Fluorine2/F2_LXCat.txt"),
    "He": ("He_LXCat.txt", "Helium/He_LXCat.txt"),
}

_KINDS = {"ELASTIC", "EFFECTIVE", "EXCITATION", "IONIZATION", "ATTACHMENT"}


def gas_file(name: str) -> Path:
    """Resolve the LXCat file for ``name`` over :func:`lxcat_roots`."""
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
    kind: str                 # elastic|effective|excitation|ionization|attachment
    target: str               # species name (e.g. 'Kr')
    product: str              # excited/ion product label, or ''
    threshold: float          # energy loss [eV] (0 for elastic/attachment)
    mass_ratio: float         # m_e / M (elastic/effective only; else 0)
    energy: np.ndarray        # tabulated energies [eV], increasing
    sigma: np.ndarray         # cross section [m^2]

    def sigma_at(self, eps: np.ndarray) -> np.ndarray:
        """Linear interpolation, zero below the first point and beyond the last."""
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


def parse_lxcat(path) -> list:
    """Parse every collision process in an LXCat file."""
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
        j += 1                                   # past the opening fence
        e_list, s_list = [], []
        while j < n and not lines[j].lstrip().startswith("---"):
            if _is_floats(lines[j]):
                a, b = lines[j].split()[:2]
                e_list.append(float(a))
                s_list.append(float(b))
            j += 1

        if e_list:
            procs.append(Process(
                kind=kind, target=target, product=product,
                threshold=threshold, mass_ratio=mass_ratio,
                energy=np.asarray(e_list, dtype=float),
                sigma=np.asarray(s_list, dtype=float)))
        i = j + 1
    return procs


class CrossSectionSet:
    """Collision processes for a mixture, with per-process target densities.

    ``species`` maps a name to ``(lxcat_path, number_density_m3)``.  The
    lumped 'effective' cross sections are DROPPED (they would double-count
    with elastic + inelastic); a file tabulating only an effective set gets a
    synthesised elastic = effective - sum(inelastic), clipped at zero.  That
    is the standard remedy and it is inherited, with its provenance, from
    pic1d.
    """

    def __init__(self, species: dict):
        self.processes = []
        self.density = []
        self.species_density = {}
        self.target_mass = {}
        for name, (path, dens) in species.items():
            self.species_density[name] = float(dens)
            procs = parse_lxcat(path)
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
                self.processes.append(proc)
                self.density.append(float(dens))
                if proc.kind == "elastic" and proc.mass_ratio > 0:
                    self.target_mass[name] = M_E / proc.mass_ratio
        self.density = np.asarray(self.density, dtype=float)
        self.n_proc = len(self.processes)
        self.kinds = [p.kind for p in self.processes]
        self.thresholds = np.array([p.threshold for p in self.processes])
        self.mass_ratios = np.array([p.mass_ratio for p in self.processes])
        self.targets = [p.target for p in self.processes]
        # fixed log-energy lookup: 8192 points over 1 meV .. 2 keV (~0.18 %
        # resolution, below the data's own accuracy).  One fancy-index gather
        # replaces n_proc binary searches per MCC call.
        self._tab_eps = np.geomspace(1e-3, 2e3, 8192)
        self._tab_log0 = np.log(self._tab_eps[0])
        self._tab_dlog = np.log(self._tab_eps[1] / self._tab_eps[0])
        self._tab_sigma = np.array([p.sigma_at(self._tab_eps)
                                    for p in self.processes])
        self._tab_nu_unit = (self._tab_sigma * self.density[:, None]
                             * np.sqrt(2.0 * self._tab_eps * E_CHARGE
                                       / M_E)[None, :])
        self._nu_runmax = None

    def _tab_idx(self, eps: np.ndarray) -> np.ndarray:
        e = np.maximum(np.asarray(eps, dtype=float), self._tab_eps[0])
        idx = ((np.log(e) - self._tab_log0) / self._tab_dlog + 0.5).astype(np.intp)
        return np.clip(idx, 0, self._tab_eps.size - 1)

    def sigma_matrix(self, eps: np.ndarray) -> np.ndarray:
        """(n_proc, n_eps) cross sections [m^2]."""
        return self._tab_sigma[:, self._tab_idx(eps)]

    def nu_matrix(self, eps: np.ndarray) -> np.ndarray:
        """(n_proc, n_eps) collision frequencies n_p sigma_p(eps) v(eps) [1/s]."""
        return self._tab_nu_unit[:, self._tab_idx(eps)]

    def nu_total(self, eps: np.ndarray) -> np.ndarray:
        return self.nu_matrix(eps).sum(axis=0)

    def nu_max(self, eps_max: float, n_grid: int = 4000) -> float:
        """Max total collision frequency over [0, eps_max] (the majorant)."""
        eps = np.linspace(1e-4, max(eps_max, 1e-3), n_grid)
        return float(self.nu_matrix(eps).sum(axis=0).max())

    def nu_majorant(self, eps_max: float) -> float:
        """Cached running-max majorant -- the MCC calls this every substep."""
        if self._nu_runmax is None:
            self._nu_runmax = np.maximum.accumulate(
                self._tab_nu_unit.sum(axis=0))
        idx = self._tab_idx(np.array([max(float(eps_max), self._tab_eps[0])]))
        return float(self._nu_runmax[int(idx[0])])

    def nu_momentum(self, eps: np.ndarray) -> np.ndarray:
        """Elastic momentum-transfer collision frequency [1/s].

        The elastic entries ARE momentum-transfer cross sections here, so
        summing the elastic rows is nu_m (see the module docstring caveat).
        """
        el = [k for k, kind in enumerate(self.kinds) if kind == "elastic"]
        return self.nu_matrix(eps)[el].sum(axis=0)

    @property
    def energy_ceiling_eV(self) -> float:
        """Highest energy at which this set is still COMPLETE [eV].

        THIS IS A SAFETY NUMBER, NOT A CURIOSITY.  ``Process.sigma_at``
        returns ZERO above a process's last tabulated point (``np.interp``
        ``right=0.0``) -- it does not extrapolate, and it does not raise.  So
        an electron above the ceiling samples no cross section at all: it
        collides with nothing, never slows, and the run looks fine.  Measured
        on this repo's files the ceiling is ~1000 eV (Ar, He), ~965 eV (Kr),
        ~937 eV (F2) -- far below the ~2 keV the lookup table spans and far
        below a 600 keV beam primary, which is why the beam module transports
        fast electrons by stopping power instead of as particles.

        Defined over the channels a fast electron needs (elastic/effective and
        ionization).  Excitation and attachment tables legitimately END where
        the channel closes, so including them would report a physical cutoff
        as a data limit -- F2's attachment stops at 20 eV.
        """
        ends = [p.energy[-1] for p in self.processes
                if p.kind in ("elastic", "effective", "ionization")]
        return float(min(ends)) if ends else 0.0

    def weight_above_ceiling(self, eps: np.ndarray, w: np.ndarray = None
                             ) -> float:
        """Fraction of ``eps`` (weighted) sitting above the validity ceiling."""
        eps = np.asarray(eps, dtype=float)
        if eps.size == 0:
            return 0.0
        w = np.ones_like(eps) if w is None else np.asarray(w, dtype=float)
        tot = float(np.sum(w))
        return float(np.sum(w[eps > self.energy_ceiling_eV]) / tot) if tot else 0.0

    def describe(self) -> str:
        counts = {}
        for k in self.kinds:
            counts[k] = counts.get(k, 0) + 1
        return (f"{self.n_proc} processes over {len(self.species_density)} "
                f"species: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))


def mixture_set(mixture: GasMixture) -> CrossSectionSet:
    """Build a :class:`CrossSectionSet` for a :class:`GasMixture`."""
    spec = {}
    for name, dens in mixture.densities.items():
        if dens > 0.0:
            spec[name] = (gas_file(name), dens)
    if not spec:
        raise ValueError("mixture_set: the mixture has no components")
    return CrossSectionSet(spec)
