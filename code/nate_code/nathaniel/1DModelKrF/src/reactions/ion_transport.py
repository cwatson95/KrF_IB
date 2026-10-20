"""Ion transport coefficients (mobility + diffusion) for the 1D charged-species fluid.

Electron transport comes from the in-tree two-term Boltzmann solve, which produces no ion
transport because it solves only the *electron* Boltzmann equation (Hagelaar
& Pitchford, Plasma Sources Sci. Technol. 14, 722 (2005), whose title is literally "…electron
transport coefficients"). Chen et al. (2025) supplied ion mobilities externally (COMSOL
requires user-input ion mobility); we do the same, from ion-mobility swarm literature.

Model (the drift-diffusion approximation, Chen Eq. 2 with ``z_j μ_j n_j E``):

- **Constant reduced mobility.** Each (ion, background-gas) pair carries a low-field *standard
  reduced mobility* ``K₀`` (cm²/V·s), referenced to the Loschmidt density
  ``N₀ = 2.6867805e25 m⁻³`` (0 °C, 1 atm). The mobility at an arbitrary gas density ``N`` is
  ``μ = K₀·(N₀/N)`` — i.e. the invariant ``μ·N = K₀·N₀`` is stored and divided by the local
  ``N`` per cell, exactly mirroring how the electron reduced mobility ``μ_e·N`` from the LUT is
  consumed. Mixtures combine available pairs via **Blanc's law**, renormalizing the represented
  neutral fractions when the table has only a subset of mixture pairs.
- **Einstein diffusion.** ``D = μ·k_B·T_gas/e`` — ions are near-thermal at the gas temperature.

Every ``K₀`` value carries a ``source`` citation key and note in the machine-readable
``data/ion_transport/ion_mobility.toml`` table.
"""

from __future__ import annotations

import tomllib
import warnings
from dataclasses import dataclass
from pathlib import Path

from particles.constants import E_CHARGE, K_B

# Loschmidt number density: standard reduced mobility K₀ reference (273.15 K, 101325 Pa).
LOSCHMIDT_M3 = 2.6867805e25
# cm²/(V·s) → m²/(V·s).
_K0_CM2_VS_TO_SI = 1.0e-4


@dataclass(frozen=True)
class IonMobilityRecord:
    """One (ion, gas) low-field reduced-mobility datum, SI-normalised at load.

    ``mu_N`` is the density-invariant reduced mobility ``μ·N = K₀·N₀`` in SI
    ``1/(V·m·s)`` — divide by the local gas number density to get ``μ`` (m²/V·s).
    """

    ion: str
    gas: str
    k0_cm2_Vs: float          # published standard reduced mobility (cm²/V·s)
    mu_N: float               # μ·N invariant (1/(V·m·s)) = K₀·1e-4·N₀
    source: str               # published-source citation key
    note: str = ""


class IonMobilityTable:
    """All (ion, gas) reduced-mobility records + mixture-aware μ/D accessors."""

    def __init__(self, records: tuple[IonMobilityRecord, ...], *, loschmidt_m3: float = LOSCHMIDT_M3) -> None:
        self.records = records
        self.loschmidt_m3 = loschmidt_m3
        # (ion, gas) -> record; and ion -> {gas: record} for Blanc's-law lookups.
        self._by_pair: dict[tuple[str, str], IonMobilityRecord] = {
            (r.ion, r.gas): r for r in records
        }
        self._by_ion: dict[str, dict[str, IonMobilityRecord]] = {}
        for r in records:
            self._by_ion.setdefault(r.ion, {})[r.gas] = r

    def has_ion(self, ion: str) -> bool:
        return ion in self._by_ion

    def reduced_mobility(self, ion: str, mixture: dict[str, float]) -> float:
        """Mixture μ·N invariant (1/(V·m·s)) for ``ion`` via Blanc's law over the neutral fill.

        ``mixture`` maps neutral species -> mole fraction. Only gases for which we hold an
        (ion, gas) record contribute; their fractions are renormalised among themselves.
        Blanc's law: ``1/μ = Σ_g (x_g / μ_g)`` → equivalently for the invariant
        ``1/(μN) = Σ_g (x_g / (μN)_g)`` since ``N`` cancels.
        """
        by_gas = self._by_ion.get(ion)
        if not by_gas:
            raise KeyError(f"no ion-mobility data for ion {ion!r}")
        contribs = {g: x for g, x in mixture.items() if g in by_gas and x > 0.0}
        if not contribs:
            # Fall back to the single most-abundant recorded gas for this ion.
            gas = max(by_gas, key=lambda g: mixture.get(g, 0.0)) if mixture else next(iter(by_gas))
            return by_gas[gas].mu_N
        total = sum(contribs.values())
        inv = sum((x / total) / by_gas[g].mu_N for g, x in contribs.items())
        return 1.0 / inv

    def mobility(self, ion: str, n_gas, mixture: dict[str, float]):
        """Ion mobility μ (m²/V·s) at gas number density ``n_gas`` (scalar or array, m⁻³)."""
        mu_n = self.reduced_mobility(ion, mixture)
        import numpy as np

        n = np.asarray(n_gas, dtype=np.float64)
        return np.where(n > 0.0, mu_n / n, 0.0)

    def diffusion(self, ion: str, n_gas, t_gas: float, mixture: dict[str, float]):
        """Einstein diffusion D = μ·k_B·T_gas/e (m²/s) at gas density ``n_gas``."""
        return self.mobility(ion, n_gas, mixture) * (K_B * t_gas / E_CHARGE)


def load_ion_mobility(path: str | Path) -> IonMobilityTable:
    """Load ``ion_mobility.toml`` into an :class:`IonMobilityTable` (K₀ → μ·N at load)."""
    path = Path(path)
    with path.open("rb") as fh:
        doc = tomllib.load(fh)
    loschmidt = float(doc.get("loschmidt_m3", LOSCHMIDT_M3))
    rows = doc.get("ion", [])
    records: list[IonMobilityRecord] = []
    for row in rows:
        k0 = float(row["K0_cm2_Vs"])
        records.append(
            IonMobilityRecord(
                ion=str(row["ion"]),
                gas=str(row["gas"]),
                k0_cm2_Vs=k0,
                mu_N=k0 * _K0_CM2_VS_TO_SI * loschmidt,
                source=str(row.get("source", "")),
                note=str(row.get("note", "")),
            )
        )
    if not records:
        warnings.warn(f"no ion-mobility records in {path}", stacklevel=2)
    return IonMobilityTable(tuple(records), loschmidt_m3=loschmidt)
