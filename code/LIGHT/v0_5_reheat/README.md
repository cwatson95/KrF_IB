# LIGHT v0_6: reheating the afterglow, tested

v0_6 is v0_5 with a window of its own for the EM wave, and a test of one
claim: that reheating the afterglow of the KrF gas with the wave improves
its efficiency (the owner, 2026-10-19, with the wave of
`IB/envelope_kr/envelope_kr_v6.tex`). FLASH 4.8 cannot test it: its laser
heats an LTE plasma by electron-ion inverse bremsstrahlung, with no
electron-atom absorption, no Kr excitation, no attachment and no KrF
kinetics. LIGHT carries every one of those and needed only the window.
Everything else is v0_5's (`v0_5/README.md`).

```bash
python3 LIGHT/run_v0_6.py --eedf two_t --em-E 33660 --em-f 6200 --em-on 1 --em-off inf --twin   # one pair, seconds
python3 LIGHT/run_v0_6.py reheat --report      # the study's tables again (results/reheat.md)
python3 LIGHT/run_v0_6.py reheat --fast        # the study under two_t, bands included (minutes)
python3 LIGHT/run_v0_6.py reheat --kinetic-bands --further --workers 8   # the whole study (hours); resumes, runs only what is missing
/usr/bin/python3 -m pytest LIGHT/v0_6/tests -q --ignore=LIGHT/v0_6/tests/test_kinetic.py   # test_kinetic.py on its own
```

## The three tests

1. **A window for the wave.** `Params.em_t_on_s` and `em_t_off_s`
   (`--em-on`, `--em-off`, in us; `inf` runs it to the end) say when the
   wave is on. The default is v0_5's wave, on with the pulse, and a
   default run is v0_5's to the bit: densities, times and ledger, checked
   against v0_5 itself under two_t. Under the kinetic EEDF the reference
   point gives v0_5's 5.73 % at 248 nm, 4.82 % VUV and 738 K (to the
   printed digits). `Model.segments` cuts the run where the given power
   stops and at both window edges.
2. **Twin pairs.** The same pump, 300 kW/cm3 for 1 us, run to 2 us,
   without the wave and with it on in the afterglow, at LIGHT's reference
   point (1 atm, 0.5 % F2) and at the envelope notes' worked example with
   fluorine (0.5 atm, 2 % F2). The wave is the notes': E0 = 47.6 kV/cm,
   the amplitude of a 3 MW/cm2 wave (33.7 kV/cm rms inside the plasma).
   It runs at 6.2 THz at the reference point (the notes' pure-Kr
   frequency) and 6.6 THz at the envelope point (the notes' value with 2 %
   F2). Two weaker fields, 1/sqrt(10) and 1/10 of that, give 1/10 and
   1/100 of the power. The windows: the whole afterglow (1-2 us) and its
   first 100 ns. Both closures: the kinetic EEDF, and v0_2's
   two-temperature closure (two_t). A 1 V/cm control twin measures the
   numerical floor.
3. **The marginal and the efficiency, with the bands.** For each pair:
   - the marginal: the 248-nm energy each absorbed joule of the wave
     adds, (E_248,on - E_248,off) / E_wave;
   - the 248-nm efficiency with the wave counted as input, E_248,on /
     (E_given + E_wave), against the pump's own. It rises exactly where
     the marginal beats the pump's efficiency.

   Every band of `network.BANDS` and the three constants the accuracy
   overview propagates (R11p, K15, C3) run at their two ends, one at a
   time, on both twins:
   - under two_t at both points, for the notes' wave and the gentle one;
   - under the kinetic EEDF at the reference point (all 20 inputs) and at
     the envelope point (the three that move its two_t marginal most,
     k_5s5p, krf_f and the trimer law, and kr2_f), for the notes' wave.

   Kr2* + F -> KrF(B), which decides the reference point, also runs at
   seven values across its band. 250 runs in all.

## What it says

**At LIGHT's defaults the reheat lowers the efficiency.** It does so at
both points, under both closures, at every field and in both windows,
with the efficiency counted with the wave as input. The notes' wave adds
248-nm light (the fluorescence rises, as REBUILD found), but each
absorbed joule makes less 248 nm than a joule of the pump does. A gentle
wave takes 248 nm away.

| point | closure | E_wave, J/cm3 (of E_given) | marginal | pump's own | change, pp | source vs pump to break even |
|---|---|---|---|---|---|---|
| reference | kinetic | 0.0093 (3.1 %) | 4.5 % | 5.73 % | -0.037 | 1.27 x |
| reference | two_t | 0.226 (75 %) | 0.44 % | 5.68 % | -2.25 | 12.9 x |
| envelope | kinetic | 2.2e-4 (0.074 %) | 13.9 % | 15.83 % | -0.0014 | 1.14 x |
| envelope | two_t | 0.0018 (0.59 %) | 9.1 % | 14.75 % | -0.033 | 1.62 x |

(The notes' wave over the whole afterglow. The last column is how
efficient the wave's source must be, from electricity to absorbed power,
relative to the pump, for the pair to break even at the wall plug:
eff_off / marginal, REBUILD's inequality eta_source x marginal >
eta_pump.)

**But one rate with no measurement decides the reference point under the
kinetic EEDF.** Kr2* + F -> KrF(B) + Kr (the band kr2_f) is 0 by default,
LIGHT's rule for a rate nobody measured. Its band runs to 1e-9 cm3/s, and
Kannari's 3e-10 traces to a conditional estimate for Ar2* + F. At the
reference point the 0.5 % F2 is burned out (below 10 % of the fill at 815
ns), so the Kr* the wave makes finds no F2 to harpoon. It becomes Kr2*,
and at the default it leaves as VUV (30 % of the wave's energy). Any
Kr2* + F channel hands that Kr2* to the F atoms the burn-out left, and
the reheat starts to pay:

| kr2_f, cm3/s | marginal | pump's own | change, pp | source vs pump | wave to VUV |
|---|---|---|---|---|---|
| 0 (LIGHT's default) | 4.5 % | 5.73 % | -0.037 | 1.27 x | 30 % |
| 1e-12 | 5.2 % | 5.80 % | -0.018 | 1.11 x | 28 % |
| 3e-12 | 6.4 % | 5.93 % | +0.014 | 0.92 x | 25 % |
| 1e-11 | 9.2 % | 6.25 % | +0.077 | 0.68 x | 19 % |
| 3e-11 | 12.3 % | 6.67 % | +0.127 | 0.54 x | 11 % |
| 1e-10 | 14.5 % | 7.14 % | +0.138 | 0.49 x | 4 % |
| 3e-10 (Kannari) | 15.3 % | 7.45 % | +0.131 | 0.49 x | 1.4 % |
| 1e-9 (the band's top) | 15.6 % | 7.67 % | +0.123 | 0.49 x | 0.4 % |

The sign turns at about 2e-12 cm3/s, under 1 % of Kannari's value. Above
it the reheat raises the efficiency by up to 0.14 point (of about 7 %, a
2 % relative gain), with a wave source at least half as efficient as the
pump. A gentle wave still takes 248 nm away there (-0.7 % at 1e-9).

**Where it does not pay at all:**
- **Under two_t.** The wave sustains a hot discharge there, and every
  band end lowers the efficiency at both points, kr2_f's top included
  (-1.37 point at the reference point).
- **At the envelope point.** The 2 % F2 attaches the electrons within
  about 10 ns, so the wave finds almost none to heat. The F2 is not
  burned out, so Kr2* + F has nothing to add: -0.0012 point at kr2_f's
  top, -0.0014 at the default.
- **With the gentle wave**, everywhere.

**So the claim holds in LIGHT only conditionally.** All of these are
needed together:
- an F2-depleted afterglow, as at the reference point;
- a Kr2* + F -> KrF(B) faster than about 2e-12 cm3/s, which no one has
  measured;
- the kinetic EEDF's regime, not two_t's sustained discharge;
- a strong wave: the notes' 3 MW/cm2 wave pays and the gentle one (1/100
  of the power) does not; the 1/10 wave was not run with kr2_f set;
- a THz source at least half as efficient as the pump.

Even then the gain is at most +0.14 point. A measurement of Kr2* + F ->
KrF(B) would settle the chemistry: above about 1e-11 cm3/s the reheat
pays at this point, and below about 1e-12 it does not. A measurement, or
a full kinetic treatment, of the EEDF under the THz field would settle
the closure.

What happens:

- **At the reference point the closures put the wave in different
  regimes.**
  - **two_t: the wave sustains a discharge through the afterglow.** The
    electron-weighted Te in the window rises from 0.61 to 2.2 eV, and
    int n_e dt grows 23-fold. The wave absorbs 75 % as much energy as
    the pump: 85 % of it ends as heat and 7 % as VUV. Tg goes from 742
    to 1120 K, past the trimer law's range and into Kr2+ and Kr2F*
    redissociation, which LIGHT lacks.
  - **Kinetic: it does not.** Te rises to 1.06 eV and int n_e dt
    2.6-fold. The wave absorbs 3.1 % of the pump's energy (at
    kr2_f = 0): 53 % to heat, 30 % to the Kr2* VUV, 11 % stored (F
    atoms, Kr*, Kr2*) and 4.5 % to 248 nm.
  - **Why.** Per electron the wave delivers 5.2e9 eV/s under two_t (at
    2.2 eV) and 2.0e9 eV/s under the kinetic EEDF (at 1.06 eV).
    Krypton's momentum transfer climbs steeply above its Ramsauer
    minimum (the notes, 'The absorption coefficient'), so a hotter EEDF
    absorbs more and is heated further. Under two_t that feedback holds
    the plasma, and the hot electrons spend the wave on the 5s <-> 5p
    cycle, whose quench is heat: 82 % of the wave goes to reaction heat
    and 21 % through Kr* excitation. Under the kinetic EEDF the feedback
    does not hold, and 69 % of the wave goes through Kr* excitation.
  - **Which to trust.** The closures differ in the EEDF alone. two_t's
    is a Maxwellian with v0_2's tail, a ratio calibrated on kinetic
    EEDFs under the prototype e-beam drive (`ML/ml/tail.py`). The
    kinetic one is solved with the wave's own field term in it, and is
    the headline.
- **At the envelope point the electrons are gone.** With 2 % F2 they
  attach within about 10 ns of the pulse's end (n_e 1.35e14 cm-3 at 1 us
  and 2.4e12 at 1.05 us, under two_t). Under two_t the wave absorbs
  0.45 % of the pump's energy in 100 ns and 0.59 % in the whole
  microsecond. It heats the few electrons it finds to 6 eV, and turns
  9 % of its energy into 248 nm against the pump's 14.7 %. Under the
  kinetic EEDF it absorbs 0.074 % in either window (Te 1.3 eV) and turns
  13.9 % of that into 248 nm against the pump's 15.8 %.
- **A gentle reheat costs 248 nm.** At 1/100 of the power (3.4 kV/cm)
  the marginal is negative at both points and under both closures (-0.3
  to -3.1 %), and at kr2_f's top too. Warming the afterglow's electrons
  a little slows the attachment and recombination whose products (F-,
  then KrF* by ion-ion recombination) make the afterglow's KrF. That is
  the finding of REBUILD's no-go audit: "every afterglow formation
  channel wants electrons colder".
- **Gain.** LIGHT has no cavity. The small-signal net gain g0 - alpha in
  the window rises with the notes' wave (kinetic reference point: 0.12
  to 0.20 %/cm), but it is 4 to 6 % of the pulse's 3.1 %/cm, so there is
  little to extract in the afterglow either way.

## The bands

`results/reheat.md` section 3 gives each input's two ends. Below, the
half-differences in quadrature, and the envelope, which adds every
input's excursions each way at once:

| point, wave, closure | marginal | change in efficiency, pp | ends that move the change most |
|---|---|---|---|
| reference, 33.7 kV/cm, kinetic (20 inputs) | 4.5 % +- 7.5 (envelope -17 to 24) | -0.037 +- 0.84 (envelope -2.57 to +0.35) | k_5s5p 0.5 (-1.55), kr2_pen 5e-10 (-0.45), kps_kr 3e-11 (-0.42); raising it: kr2_f 1e-9 (+0.12), kstep_s 1.5e-7 (+0.071), kstep_p 2e-7 (+0.015), k_5s5p 1.4 (+0.002) |
| reference, 33.7 kV/cm, two_t | 0.44 % +- 2.0 (envelope -1.6 to 5.1) | -2.25 +- 2.0 (envelope -8.0 to +2.7) | kstep_s 1.5e-7 (-5.09), k_5s5p 0.5 (-3.49), trimer T^-0.5 (-0.92), kr2_f 1e-9 (-1.37) |
| reference, 3.37 kV/cm, two_t | -1.42 % +- 1.05 | -0.0003 +- 0.0006 | |
| envelope, 33.7 kV/cm, kinetic (3 inputs, and k_5s5p's upper end) | 13.9 % +- 5.4 (envelope -1.0 to 14.3) | -0.0014 +- 0.0004 (envelope -0.0019 to -0.0006) | trimer T^-0.5 (-0.0008), krf_f 7.8e-10 (-0.0019), kr2_f 1e-9 (-0.0012), k_5s5p 1.4 (-0.0014); its 0.5 end did not finish |
| envelope, 33.7 kV/cm, two_t | 9.1 % +- 5.1 (envelope -5.6 to 18.6) | -0.033 +- 0.017 (envelope -0.073 to +0.025) | k_5s5p (-0.036 at 0.5, -0.019 at 1.4), trimer T^-0.5 (-0.020) |
| envelope, 3.37 kV/cm, two_t | -3.1 % +- 3.4 | -0.0001 +- 0.0001 | |

Under two_t every input's end, one at a time, lowers the efficiency at
both points with either wave; only the envelope, which stacks every
favourable end at once, reaches above zero. Under the kinetic EEDF the
envelope point's every end lowers it. At the reference point four ends
raise it, three of them rates with no measurement (kr2_f, kstep_s,
kstep_p) and one the measured 5s -> 5p cross section at its upper error.

## What was found on the way

- **The wave alone must not time-step the EEDF in an attaching
  afterglow.** v0_6's first version advanced the EEDF in time whenever
  a drive was on, the wave alone included. At the envelope point three
  of six runs failed:
  - two found no time step inside the window, the 1 V/cm control one of
    them, so the wave's power was not the cause;
  - one found no steady solution as the window closed at 1.1 us.

  The EEDF is now time-stepped under the given power alone, v0_5's rule.
  With the wave alone it is the steady solution at LIGHT's mean with the
  wave's field term, as after the pulse. Where no steady solution
  exists, a time step from the last one stands in, logged as `stepped`:
  2 of the study's 37,097 kinetic solutions, both in the envelope
  point's early-window run.
- **v0_5's call signatures are part of the interface.** The beam drive
  of `docs/survey_2026-09-27/scripts/ebeam_screen.py` overrides
  `rhs(t, y, on)`, and `test_kinetic.py`'s counters wrap
  `_kinetic_update(t, y, on)`. A fourth argument everywhere broke four
  tests. The wave's flag is now passed only where a window parts it from
  the given power, so every run without a window calls them as v0_5 did.
- **The holding field.** With the wave alone, the kinetic EEDF at
  LIGHT's mean needs a DC field wherever the wave's own steady solution
  would be colder. That field is no input LIGHT books. Its power, as a
  share of the wave's, is logged:
  - median 0.002 with the notes' wave at the reference point, and 0.001
    at the envelope point;
  - about 1 in the first solutions after the pulse, while LIGHT's Te
    still carries the pump's;
  - of order 1 or more with the gentle wave, where the wave's own power
    is small against what holds the afterglow's electrons in v0_5 too.
- **A control twin.** A 1 V/cm wave cuts and renews the run as the wave
  does and absorbs nothing. 248 nm moves by 1e-12 J/cm3, so the pairs'
  differences are physics, not numerics.
- **Running it.** The study's 250 runs ran on 2026-10-19 across two
  sessions and a reboot, on 8 to 14 processes. The envelope point's
  kinetic runs take 20-75 min of compute each: its afterglow renews the
  EEDF every time n_e moves 10 %. The memory reaper stopped one batch,
  with 28 workers beside another session's decompilers. `reheat.py`
  saves each run as it ends, so nothing was lost.
- **One run did not finish.** The envelope point's kinetic pair at
  k_5s5p = 0.5 (its off twin finished) had computed for 3.4 hours, twice
  its siblings' longest, when it was stopped. Whether it stalls or only
  crawls is not known. `python3 LIGHT/run_v0_6.py reheat --kinetic-bands
  --further` runs it again. Its band end is the one the envelope point's
  kinetic band row lacks.

## Limits

- **0-D.** The field inside the plasma is given. The notes' design is
  not modelled (the gas 10 cm deep, absorbing one e-folding of the
  wave), nor reflection, nor how the absorption length changes as n_e
  falls in the afterglow.
- **No cavity.** The 248-nm energy is fluorescence, and the optics
  block is the small-signal diagnostic. REBUILD's one positive afterglow
  case was a cell still lasing through the window; LIGHT cannot see that
  channel.
- **The closures disagree 24-fold** on the absorbed energy at the
  reference point. The kinetic EEDF solves the wave's own field term,
  but with the wave alone it is held at LIGHT's mean by a field LIGHT
  does not book (above), and the EEDF's transient as the wave turns on
  is not followed.
- **Hot gas.** The envelope point runs hot without the wave (Tg 1050 K
  at the end, in both twins), past the trimer law's range and into Kr2+
  and Kr2F* redissociation, which LIGHT lacks. Under two_t the notes'
  wave takes the reference point there too (1120 K).
- **Bands one at a time.** Their envelope is not a joint distribution.
  kr2_f, which moves the result most, was also scanned across its band;
  the others were not.

## Checks

`tests/test_v0_6.py` (15 tests):
- the default window equal to v0_5 to the bit (three cases, against
  v0_5 itself);
- the pieces `segments` cuts, the windows refused, and a window past
  t_end equal to no wave;
- the afterglow wave booking its energy in its window alone, with the
  run before it unchanged, and a window inside the pulse;
- a subclass with v0_5's three arguments;
- the kinetic EEDF seeing the wave in its window alone (about three
  minutes);
- the twin-pair arithmetic, `unc_stats`, and the study's specs;
- the two_t reference pair pinned;
- the CLI.

On 2026-10-19, on the final code, under `/usr/bin/python3` (3.9, numpy
2.0):
- `test_v0_6.py`: 15 passed.
- The other fast files: 327 passed and 1 strict expected failure
  (Hakuta's decays, as in v0_5).
- `test_kinetic.py`: 24 passed and 1 skipped, as in v0_5.

Every saved run of the study closes its ledger to 1e-12 of its input
(the largest residual 8e-13).

## Files

```
v0_6/
  README.md, SOURCES.md    this file; the sources (section 9 is v0_6's)
  model.py       + Params.em_t_on_s, em_t_off_s; Model.segments; the
                 pump's and the wave's flags apart in rhs and the kinetic
                 EEDF (v0_5's call where they agree); the steady-to-
                 time-step fallback
  report.py      the window in the summary, the wave-alone kinetic line,
                 row_dict's given_J_cm3, E248_J_cm3, em_on_us, em_off_us
  reheat.py      the twin pairs: specs, runs (resumable, one process
                 each), the pair's numbers, the bands, FURTHER, the tables
  run.py         + --em-on, --em-off, --twin; `run_v0_6.py reheat`
  results/       reheat_runs.json (every run's record), reheat.md (the
                 tables)
  tests/         test_v0_6.py and v0_5's, re-pointed (test_docs.py left
                 out: v0_6 has no documents of its own)
```
