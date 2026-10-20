"""PIC_MCC_3D -- a self-consistent 3-D electrostatic particle-in-cell /
Monte-Carlo-collision code.

The roadmap is ``../PIC_MCC_3D.txt``; this package implements the nine
"missing modules" of its section 2.1, in dependency order.

    MODULES/        every dependency of the simulation (import from here)
    tests/          the falsification gates -- python3 -m pytest PIC_MCC_3D/tests -q
    status.html     the live system map: groups, sub-groups, status, gates
    CHANGELOG.md    the dated record + the permanent Known-bad list -- READ FIRST

Import path (single, explicit -- see MODULES/__init__.py for why there is no
root-level alias):

    from PIC_MCC_3D.MODULES.es_sim import ESSim

LANDED: missing modules 1 (electrostatic Poisson core) and 2 (mobile
multi-species particles).  STILL OPEN: ion-neutral collisions, electrodes
with an external circuit, secondary emission, Coulomb collisions, the
implicit mover, merge/split, beam injection.  There is no MCC here yet.

RELATION TO ``ARPAE/pic3d``: that package is the ELECTROMAGNETIC instrument
(Yee FDTD, horn antenna, wave coupling into a PRESCRIBED plasma).  This one is
the ELECTROSTATIC, self-consistent code.  Neither imports the other yet; when
collisions land (M3) the tested ARPAE MCC and cross-section modules are what
should be adapted, not re-derived.

POLICY (PLAN (26)): predictive-only.  Analytic solutions and published
benchmarks SCORE this code; they never steer a number.
"""

from __future__ import annotations

__all__ = ["MODULES"]
__version__ = "0.2.0"
