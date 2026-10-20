"""0-D excimer laser gain/output optics (Chen et al. 2025 §II.C).

The optical (gain) axis runs transverse to the 1D discharge axis, so the photon field is a
single cavity photon density ``n_ph`` (:class:`PhotonDensity`, Chen Eq. 7) driven by the plasma,
and the laser output follows from the lumped-gain relations (:class:`GainCavity`, Eqs. 10/14/18).
Valid for both KrF* (248 nm) and ArF* (193 nm).
"""

from __future__ import annotations

from .gain import GainCavity
from .photon_density import PhotonDensity

__all__ = (
    "GainCavity",
    "PhotonDensity",
)
