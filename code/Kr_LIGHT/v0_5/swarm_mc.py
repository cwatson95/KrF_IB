"""What a swarm experiment measures that the two-term swarm does not give
(LIGHT_Kr v0_3): the Monte Carlo's ratios.

swarm.py solves the two-term Boltzmann equation with the temporal growth
term: its drift velocity is the flux one and its ionization coefficient
k_i/W. A drift tube times the swarm's centre of mass (the bulk drift
velocity), a steady-state Townsend experiment measures the spatial growth
of its current (the SST coefficient), and the two-term transverse
diffusion leaves out the higher harmonics. LoKI-MC (Dias et al., Comput.
Phys. Commun. 282, 108554 (2023); the vendored LoKI-MC/Code/lokimc), run on
Biagi's 8.97 file as published (xsec.UNIT) at 300 K and 1 Torr with equal
sharing (1e5 electrons; 2e4 at 0.2 Td, 1e4 at 0.5 and 1 Td), gives each
as its value over the two-term's at E/N [Td] (2026-10-14 and 15):

  W       the bulk drift velocity over the two-term's flux one. They part
          as ionization grows (+23 % at 300 Td); the Monte Carlo's flux
          drift velocity is the two-term's to 0.3 % from 10 to 300 Td.
          Below 20 Td the 20-Td value (0.9992): without ionization the bulk
          and the flux drift velocities are one, and the two-term solution
          is compared as it is. Its own error there (the Monte Carlo's
          flux W is 1.8 % above it at 0.2 Td, 0.3 % below at 10 Td; FLUX)
          is absorbed by the trained factors, which are effective for
          LIGHT_Kr's two-term solver.
  alphaN  the SST coefficient over the two-term's k_i/W: the two-term's is
          the Monte Carlo's flux one less 1.5 to 4 %, and the SST one lies
          below the flux one by up to 11 % at 300 Td. Below 40 Td the
          40-Td value: the Monte Carlo counts 29 and 543 ionizations at 20
          and 30 Td (19 and 4 % noise).
  DTmu    the flux characteristic energy D_T/mu over the two-term's: the
          two-term transverse diffusion is 3.6 % high at 0.2 Td, 10 % at
          0.5 Td, 16 % at 1 Td, 46 % at 10 Td and 60 to 80 % from 20 to
          300 Td (DT_LOW below 20 Td; below 0.2 Td the 0.2-Td value,
          where the ratio must reach 1 in the thermal limit: a
          thermal-limit alternative, 1 at 0.005 Td, moves Koizumi's chi2
          by +0.1 to +1.1 per point and keeps the variants' ranking,
          ML/CHANGELOG.md 2026-10-15).

ratio(kind, en) reads them log-linearly in E/N and holds them at the
tables' ends. They are taken on Biagi's file as published. On v0_3's
trained set (xsec.TRAINED_FILE) the Monte Carlo gives W and alphaN ratios
within 1.65 % of these from 30 to 300 Td (TRAINED_CHECK; the most is the
coefficient at 50 Td, whose counting error is about 1 %), where the
data's errors are 3 and 10 %.
"""
from __future__ import annotations

import numpy as np

# E/N [Td], W, alphaN, DTmu
RATIOS = (
    (20.0, 0.99916, 1.03579, 0.62980),
    (30.0, 1.00318, 1.03579, 0.59785),
    (40.0, 1.00595, 1.03579, 0.57769),
    (50.0, 1.01113, 1.02521, 0.57178),
    (75.0, 1.03819, 0.99747, 0.55844),
    (100.0, 1.06468, 0.97928, 0.56065),
    (150.0, 1.11315, 0.95598, 0.58680),
    (200.0, 1.15249, 0.93935, 0.61525),
    (300.0, 1.23178, 0.90618, 0.68065),
)
# E/N [Td], DTmu below 20 Td
DT_LOW = (
    (0.2, 0.96527),
    (0.5, 0.90522),
    (1.0, 0.86223),
    (10.0, 0.68576),
)
# the Monte Carlo's flux values over the two-term's on Biagi's file, E/N
# [Td], W_flux, alpha_flux (the same runs): how close the two-term solution
# itself is; at 20 and 30 Td the coefficient counts 29 and 543 ionizations
FLUX = (
    (0.2, 1.01761, None),
    (1.0, 1.00352, None),
    (10.0, 0.99688, None),
    (20.0, 0.99897, 0.96407),
    (30.0, 1.00115, 1.04330),
    (40.0, 0.99976, 1.03940),
    (50.0, 0.99683, 1.03392),
    (75.0, 0.99881, 1.02000),
    (100.0, 0.99895, 1.01554),
    (150.0, 0.99978, 1.01425),
    (200.0, 0.99874, 1.01678),
    (300.0, 1.00039, 1.01667),
)
# the same on v0_3's trained set (the seven factors), E/N [Td], W, alphaN
# (2026-10-15)
TRAINED_CHECK = (
    (30.0, 1.00220, 1.03445),
    (50.0, 1.01772, 1.00832),
    (100.0, 1.07031, 0.97336),
    (200.0, 1.16533, 0.92985),
    (300.0, 1.24076, 0.90790),
)
KINDS = ("W", "alphaN", "DTmu")


def ratio(kind: str, en) -> np.ndarray:
    """The Monte Carlo's value over the two-term's for one kind at E/N."""
    en = np.atleast_1d(np.asarray(en, dtype=float))
    col = KINDS.index(kind) + 1
    xs = [r[0] for r in RATIOS]
    ys = [r[col] for r in RATIOS]
    if kind == "DTmu":
        xs = [r[0] for r in DT_LOW] + xs
        ys = [r[1] for r in DT_LOW] + ys
    return np.exp(np.interp(np.log(en), np.log(xs), np.log(ys)))


def as_measured(d: dict, en) -> dict:
    """swarm.point's W, alphaN and DTmu as the experiments measure them."""
    return {k: float(d[k] * ratio(k, en)[0]) for k in KINDS}
