"""LIGHT_Kr v0_3: the lightweight 0-D model of pure krypton (2026-10-15).

v0_3 is v0_2 with its electron cross sections trained on krypton's
measured swarm data by ML (ML/ml/kr_swarm.py, the owner's "use ML to
train it to improve the accuracy"): Biagi's Magboltz 8.97 file stays the
source, and seven factors on it (the elastic momentum transfer at four
knots from 0.6 to 20 eV, Biagi's own at and below 0.2 eV, and one factor
on each group of excitation levels) are fitted to the drift velocities
and ionization coefficients that ML's consistency check admits, with a
neural emulator of the swarm solver, an MCMC and a neural posterior
estimator; the factors are pinned in data/ (xsec.TRAINED_FILE). The swarm
solver also gives the transverse diffusion (swarm.py). Everything else is
v0_2's: every level of Biagi's set its own row, the Hornbeck-Molnar band
hm_kr, LIGHT v0_5's krypton network and machinery. See
LIGHT_Kr/v0_3/README.md.

Canonical imports: ``from v0_3.model import Params, run`` and
``from v0_3.report import summary``. No import-time side effects.
"""

__version__ = "0.3.0-kr"
