"""LIGHT v0_6: v0_5 with a window of its own for the EM wave, so that the
wave can reheat the afterglow with the given power off (2026-10-19).

v0_6 starts as a copy of v0_5. The wave of v0_4 rode with the given
power, on from 0 to pulse_s. Params.em_t_on_s and em_t_off_s set its
window, and Model.segments splits the run where the given power stops and
at both window edges. The default window is v0_5's, and a default run is
v0_5's to the last digit. reheat.py runs the twin pairs, the same pump
with and without the wave in the afterglow, at LIGHT's reference point and
at the envelope notes' worked example (IB/envelope_kr/envelope_kr_v6.tex),
and reports the 248-nm energy each absorbed joule of the wave adds and the
efficiency with the wave counted as input, with the bands. v0 to v0_5
stay as they were. See LIGHT/v0_6/README.md.

Canonical imports: ``from v0_6.model import Params, run`` and
``from v0_6.report import summary``. No import-time side effects.
"""

__version__ = "0.1.6"
