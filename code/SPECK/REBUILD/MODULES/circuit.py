"""The discharge circuit — Chen 2025's C0/Cp/L chain as engine state.

Closes the anchored x3-7 formation-efficiency gap's named cause: the
prescribed-power stand-in had no E/N-driven breakdown kinetics. Here the
gap field is a STATE, conduction follows the swarm drift velocity from the
same LoKI tables the rates come from (natively E/N-indexed — the circuit
path uses k(E/N) directly, with NO <eps> = (3/2)Te mapping at all), and
the avalanche grows or does not grow on its own.

Topology (Chen 2025 Sec. II A — the MPC loop + discharge load loop of the
BASE CASE all results figures use: C0 = 10 nF, Lm = 250 nH, Cp = 10 nF,
Ls = 5 nH, V0 = -18 kV):

    C0 --[ Lm, R ]--> Cp --[ Ld, R ]--> gap (d, A) ++ [ Rf, Lf, Cf ]

States: q0, qp, qf [C], i1, i2, if [A], E [V/m] (the gap field).

    dq0/dt = -i1
    dqp/dt =  i1 - i2
    di1/dt = (q0/C0 - qp/Cp - R i1) / Lm
    di2/dt = (qp/Cp - E d   - R i2) / Ld
    dif/dt = (E d - qf/Cf - Rf if) / Lf      [the paper's RLC filter
    dqf/dt =  if                              branch across the gap,
                                              Rf 70.71 ohm / Lf 100 nH /
                                              Cf 10 pF — "to suppress
                                              circuit oscillations after
                                              gap breakdown"; 0.1.27]
    dE/dt  = (i2 - j_cond A - if) / (eps0 A) [current continuity at the
                                              gap node: circuit =
                                              conduction + displacement
                                              + filter draw]

with j_cond = e n_e v_drift(E/N) sign(E) supplied by the engine per step.

Stated simplification vs the paper (0.1.25): Lm is a fixed inductor at
the paper's stated 250 nH ("a saturable inductor with a low inductance"
— the 5 -> 1 uH switching curve in Sec. II D belongs to the Fig. 5
VALIDATION circuit, whose 10.4 nF/8.1 nF/30 nH components shipped here
as "the Chen circuit" through 0.1.24; see CHANGELOG [0.1.25]). The RLC
filter branch is IN as of 0.1.27 (it was simplification (b)). R is the
series loop resistance (copper/ESR class, a declared physical parameter
with a range, not a fit); Rf/Lf/Cf are the paper's printed filter values.

Energy bookkeeping: the CELL budget boundary is j.E (Joule heating into
the electron gas — the ledger's drive_in). Everything upstream of the gap
(capacitor charge, inductor currents, the gap's own field energy) is
circuit-side and reported via :func:`circuit_energy_J`; with R = 0 and no
plasma the total is conserved (the LC gate).
"""

from __future__ import annotations

from dataclasses import dataclass

from .units import EPS0

N_CIRCUIT_STATES = 7
IQ0, IQP, II1, II2, IE, IIF, IQF = range(N_CIRCUIT_STATES)


@dataclass(frozen=True)
class CircuitSpec:
    c0_F: float
    cp_F: float
    lm_H: float
    ld_H: float
    v0_V: float          # initial charge voltage on C0 (magnitude)
    r_ohm: float
    gap_m: float         # discharge gap d (the field axis)
    area_m2: float       # discharge area A
    # the paper's RLC filter branch across the gap (Sec. II A; 0.1.27)
    rf_ohm: float = 70.71
    lf_H: float = 100e-9
    cf_F: float = 10e-12

    def initial(self):
        import numpy as np
        y = np.zeros(N_CIRCUIT_STATES)
        y[IQ0] = self.c0_F * self.v0_V
        return y

    def derivs(self, y, i_cond_A):
        """d(states)/dt given the conduction current [A] through the gap."""
        import numpy as np
        q0, qp, i1, i2, E, i_f, q_f = (y[IQ0], y[IQP], y[II1], y[II2],
                                       y[IE], y[IIF], y[IQF])
        d = np.empty(N_CIRCUIT_STATES)
        d[IQ0] = -i1
        d[IQP] = i1 - i2
        d[II1] = (q0 / self.c0_F - qp / self.cp_F
                  - self.r_ohm * i1) / self.lm_H
        d[II2] = (qp / self.cp_F - E * self.gap_m
                  - self.r_ohm * i2) / self.ld_H
        d[IIF] = (E * self.gap_m - q_f / self.cf_F
                  - self.rf_ohm * i_f) / self.lf_H
        d[IQF] = i_f
        d[IE] = (i2 - i_cond_A - i_f) / (EPS0 * self.area_m2)
        return d

    def energy_J(self, y) -> float:
        """Total circuit-side energy: capacitors + inductors + gap field
        + the filter branch's Lf/Cf stores (Rf dissipation is a
        circuit-side loss, like R)."""
        q0, qp, i1, i2, E, i_f, q_f = (y[IQ0], y[IQP], y[II1], y[II2],
                                       y[IE], y[IIF], y[IQF])
        return float(q0 * q0 / (2 * self.c0_F) + qp * qp / (2 * self.cp_F)
                     + 0.5 * self.lm_H * i1 * i1
                     + 0.5 * self.ld_H * i2 * i2
                     + 0.5 * self.lf_H * i_f * i_f
                     + q_f * q_f / (2 * self.cf_F)
                     + 0.5 * EPS0 * E * E * self.gap_m * self.area_m2)
