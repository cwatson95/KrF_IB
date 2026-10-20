"""Excimer-laser efficiency from interval-averaged optical and discharge powers.

The live fluid solver records one average power for every chemistry/optics interval. The
initial electron thermal energy is included in the deposited-energy denominator. Pulse-integrated
metrics:

- **intrinsic efficiency** = ∫P_out / (E_initial + ∫P_in) (optical out per unit total
  energy deposited in the discharge, whether stored before t=0 or coupled in by the
  drive during the pulse);
- **wall-plug efficiency** = ∫P_out / E_electrical_total (optical out per unit total
  electrical energy delivered to the device).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray


@dataclass
class LaserEnergetics:
    """Time series of output/input power with pulse-integrated efficiency metrics.

    ``power_out_w`` (total coupler output) is stored explicitly alongside its
    spontaneous/stimulated components; ``record`` maintains the invariant
    ``out = out_spont + out_stim``.
    """

    initial_energy_j: float = 0.0  # energy stored in the plasma before the first step
    times_s: list[float] = field(default_factory=list)
    power_out_w: list[float] = field(default_factory=list)
    power_out_spont_w: list[float] = field(default_factory=list)
    power_out_stim_w: list[float] = field(default_factory=list)
    power_in_w: list[float] = field(default_factory=list)
    interval_s: list[float | None] = field(default_factory=list)

    def record(
        self,
        time_s: float,
        power_out_spont_w: float,
        power_out_stim_w: float,
        power_in_w: float,
        interval_s: float | None = None,
    ) -> None:
        """Append one record; powers are averages over ``interval_s`` when supplied."""
        if interval_s is not None and interval_s < 0.0:
            raise ValueError("interval_s must be non-negative")
        self.times_s.append(time_s)
        self.power_out_spont_w.append(power_out_spont_w)
        self.power_out_stim_w.append(power_out_stim_w)
        self.power_out_w.append(power_out_spont_w + power_out_stim_w)
        self.power_in_w.append(power_in_w)
        self.interval_s.append(interval_s)

    @property
    def t(self) -> NDArray[np.float64]:
        return np.asarray(self.times_s, dtype=np.float64)

    def _integrate(self, series: list[float]) -> float:
        if self.times_s and len(self.interval_s) == len(self.times_s) \
                and all(dt is not None for dt in self.interval_s):
            return float(sum(value * float(dt) for value, dt in zip(series, self.interval_s)))
        if len(self.times_s) < 2:
            return 0.0
        return float(np.trapezoid(series, self.times_s))

    def output_energy_j(self) -> float:
        """Pulse-integrated optical output energy."""
        return self._integrate(self.power_out_w)

    def output_energy_spont_j(self) -> float:
        """Spontaneous-labelled part of the output energy."""
        return self._integrate(self.power_out_spont_w)

    def output_energy_stim_j(self) -> float:
        """Stimulated-labelled part of the output energy."""
        return self._integrate(self.power_out_stim_w)

    def input_energy_j(self) -> float:
        """Pulse-integrated discharge energy ∫P_in dt."""
        return self._integrate(self.power_in_w)

    def deposited_energy_j(self) -> float:
        """Total deposited energy: pre-pulse stored energy + net drive work."""
        return self.initial_energy_j + self.input_energy_j()

    def intrinsic_efficiency(self) -> float:
        """∫P_out / (E_initial + ∫P_in) (0 if no deposited energy)."""
        e_dep = self.deposited_energy_j()
        return self.output_energy_j() / e_dep if e_dep > 0.0 else 0.0

    def wall_plug_efficiency(self, total_electrical_energy_j: float) -> float:
        """∫P_out / total electrical energy delivered (0 if none)."""
        return self.output_energy_j() / total_electrical_energy_j if total_electrical_energy_j > 0.0 else 0.0
