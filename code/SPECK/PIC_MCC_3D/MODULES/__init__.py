"""PIC_MCC_3D.MODULES -- the code that makes up the PIC-MCC.

Every dependency of the simulation lives here, so the package root holds only
the things ABOUT the code (CHANGELOG.md, status.html) and the harness that
checks it (tests/).  The canonical import path is explicit:

    from PIC_MCC_3D.MODULES.poisson3d import Poisson3D

There is deliberately NO re-export of these names at the package root.  A
convenience alias would make ``from PIC_MCC_3D import poisson3d`` work while
``from PIC_MCC_3D.poisson3d import Poisson3D`` still failed -- one path that
works and one that does not, for the same module, is worse than one path.

WHAT IS HERE (the landed half of ../PIC_MCC_3D.txt section 2.1):

    units.py        SI constants, plasma scales, the ES/EM timestep ratio
    poisson3d.py    MISSING MODULE 1 -- the separable 3-D Poisson solve
    es_field.py     deposit / gradient / gather (one shape function, shared)
    species3d.py    MISSING MODULE 2 -- multi-species kinetic particles
    es_sim.py       the closed deposit -> solve -> gather -> push loop
    theory.py       analytic overlays that SCORE the code, never steer it

STILL OPEN: ion_mcc (M3), boundaries3d (M4), emission3d (M5), coulomb3d (M6),
implicit3d (M7), resample (M8), beam3d (M9).  ``../status.html`` is the live
board; ``../CHANGELOG.md`` is the record.
"""

from __future__ import annotations

__all__ = ["units", "poisson3d", "es_field", "species3d", "es_sim", "theory"]
