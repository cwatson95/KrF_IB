# LIGHT v0_4

v0_4 is v0_3 with what the audit of 2026-10-08 found missing against the
measured KrF data and the modern kinetic codes. The owner asked that day
for "the absorbers, quenching, ladder, Turnbull's IB equations (note that
it has a 10 % error), ion drift, and ion-ion recombination". Each is in,
each at a measured value where one exists and as a band where none does,
and the kinetic EEDF of v0_3 carries the new electron collisions. v0_3
stays as it was.

```bash
python3 LIGHT/run_v0_4.py --eedf two_t                 # the reference point in seconds (v0_2's closure)
python3 LIGHT/run_v0_4.py                              # the same with the kinetic EEDF, minutes
python3 LIGHT/run_v0_4.py --bands --eedf two_t         # every band's two ends, 4 min
python3 LIGHT/run_v0_4.py --ladder kannari             # the 5p rung's bands at Kannari's values
python3 LIGHT/run_v0_4.py --rate krf_f=7.8e-10         # KrF* + F at Kannari's analogue
python3 LIGHT/run_v0_4.py --gap 1                      # a 1-cm gap: the charges drift to the walls
python3 LIGHT/run_v0_4.py --em-E 1500 --em-f 95        # an EM wave, 1500 V/cm rms at 95 GHz
LIGHT_VERSION=v0_4 python3 LIGHT/docs/krf_data/tahoma_v0_4.py [--kinetic]   # Tahoma g0, alpha0
python3 -m pytest LIGHT/v0_4/tests -q                  # the tests (test_kinetic.py takes minutes)
```

## What changed

1. **The 5p rung.** A lumped Kr** (the 4p^5 5p manifold and the levels
   that cascade through it, 11.30 eV, the Biagi group's own threshold) is
   made by the higher excitations C1a-c and by Kr2+ + e (Shiu & Biondi saw
   the 5p lines), relaxes to Kr* in collisions with Kr (KPS, the measured
   (1.7 +- 0.2)e-10 cm3/s of Gazeli et al. 2021 for one 5p' level),
   radiates to the 5s levels at the measured manifold lifetime (PRAD, 29.0
   ns, the statistical mean of the ten measured level lifetimes), harpoons
   F2 at the measured (13 +- 2)e-10 (K27, Nelson & Setser 1990 as a search
   snippet renders their abstract: verify), and is pumped from Kr* by the
   electrons on Kannari's fit, 9.0e-7 Te^0.5 exp(-1.51/Te) (K76, Jacob &
   Mangano 1976), and returned by detailed balance (K76r). The two stepwise
   ionizations (KSTEP_S from 5s, KSTEP_P from 5p) have no measurement and
   are bands at 0; `--ladder kannari` sets every ladder band to Kannari's
   printed value at once. In the kinetic EEDF the 5s <-> 5p pair is a
   population process on Kr*: Kannari's A Te^0.5 exp(-C/Te) is exactly the
   Maxwellian rate of sigma(eps) = sigma0 (1 - C/eps), so that pseudo
   cross section (threshold the ledger's 1.385 eV gap, zero below C)
   reproduces the printed law under a Maxwellian, gives the EEDF the
   collision's shape, and makes K76r its Klein-Rosseland reverse
   (`boltzmann.PopProcess`).
2. **KrF(B,C) + F** (QF, QFC): a band, 0 to Kannari's 7.8e-10 'analogous
   to F2', 0 by default. No measurement of it, or of any RgX* + X
   analogue, exists.
3. **The 248-nm optics** (`optics.py`): g0 = sigma_B n_B and alpha =
   sum sigma_i n_i, species by species, from Kannari's absorber block
   (F- 5.6e-18 measured, Kr2+ 1.6e-18 an estimate, Kr** 6.0e-18 a
   calculation, F2 1.2e-20, Kr* 1.0e-19 compiled) with Hakuta et al.'s
   1987 direct measurement of Kr2F* at 248 nm, (1.0 +- 0.2)e-18, in place
   of the compilation's 6.0e-18 that every KrF model carried. A
   diagnostic: LIGHT has no photon field, so the absorbers change no
   density.
4. **Ion-ion recombination** (K26, K108 and their new C and dissociation
   rows): Flannery's turnover as the harmonic mean of REBUILD's three-body
   (Thomson) limb and the Langevin limb of the measured reduced mobilities
   (Ellis et al.; the same table Fluid1D's transport reads), which lowers
   the 1-amagat rate from 3.04e-6 to 2.60e-6 (Kr+) and 2.74e-6 (Kr2+)
   cm3/s. The KrF* yield (ii_yield, 1.0 by default, Rokni's relative 'near
   unity', down to 0.5) and its C share (ii_f_c, 0 as in v0_3, up to
   REBUILD's even split) are bands that multiply the three product rows.
5. **The ions' drift to the walls** (`Params.gap_cm`, `--gap`): in a gap
   the cations leave at the ambipolar rate D_a/Lambda^2, D_a = mu_i k(Te +
   Tg)/e with the measured mobilities, F- at its free diffusion, Lambda =
   gap/pi; the lost charges' formation energy and 2 Te per electron go to
   the walls (a new fate line) and Te per pair to the gas through the
   ambipolar field. None without a gap (the box of v0 to v0_3).
6. **The EM wave** (`--em-E`, `--em-f`): a second energy input whose
   absorbed power is an output (E_em in the ledger). The electron-atom
   part is v0's Lorentz-gas kernel on the Biagi/SIGLO momentum transfer
   (within 1.6 % of the quantum rate at 0.3-4 THz, IB/electron_atom); the
   electron-ion part is Turnbull et al.'s inverse bremsstrahlung (`ib.py`:
   the Sommerfeld logarithm with Pradler & Semmelrock's Gaunt-factor fit,
   ion screening at the ion-sphere radius, the Langdon factor), with
   their model's stated 10 % error as the band ib_scale (0.9 to 1.1) on
   the electron-ion share of the power. Under the kinetic closure the wave
   enters the Boltzmann equation as Hagelaar & Pitchford's high-frequency
   field term with Turnbull's per-energy cross section in its momentum
   transfer, and the absorbed power is the solver's own. Two
   double-counting rules: the Langdon factor is applied only under the
   Maxwellian closures (a solved EEDF already holds the distortion it
   corrects), and the ions' collisions enter each field term once (the DC
   stand-in field sees them Debye-screened, the wave with the IB
   logarithm).
7. **Electron-ion collisions in the kinetic EEDF** at every state, with or
   without a wave: the cations' Coulomb momentum transfer, Debye-screened,
   in the DC field's term (`boltzmann._sigma_ei_dc`).

The ledger closes to 1e-13 in every new mode (test_v0_4.py). With the 5p
rung neutralised (k_5s5p = 0 and the 5p relaxed at once) v0_4 reproduces
v0_3's reference point to 0.03 %, the ion-ion limb change included.

## At the reference point

1 atm, 0.5 % F2, 300 kW/cm^3 for 1 us, run to 2 us, in percent of the
input, under v0_2's closure (two_t) and the kinetic EEDF:

| | v0_3, two_t | v0_4, two_t | v0_4, kinetic |
|---|---|---|---|
| 248 nm | 4.67 | 5.79 | 5.68 |
| 147 nm (Kr2*) | 20.0 | 4.53 | 5.29 |
| 420 nm (Kr2F*) | 4.87 | 3.84 | 4.28 |
| 5p-5s near IR | - | 0.28 | 0.27 |
| heat, reactions | 53.2 | 71.6 | 71.9 |
| g0, 500-1000 ns, %/cm | 2.56 | 4.13 | 3.73 |
| alpha, %/cm | 0.70 | 0.63 | 0.68 |
| Tg end, K | 649 | 742 | 735 |
| Te pulse mean, eV | 4.45 | 5.81 | 2.97 |

(v0_3's kinetic reference gave 4.79 % at 248 nm.) Under the kinetic EEDF
the field carries 1.001 to 1.004 of the input through the pulse and the
kinetic energy balance agrees with LIGHT's to 7.9e-6 of it, with the 5s
<-> 5p pair inside the Boltzmann equation; 463 solutions, 311 in the
pulse, E/N 1.6 to 45 Td.

What moves the result is the 5s <-> 5p cycle. On Kannari's fit the
electrons pump Kr* to the 5p rung faster than F2 harpoons it (81 % of the
Kr* formed goes up, 8 % to KrF), the 5p relaxes to 5s on Kr in 0.2 ns
(KPS, 96 % of its fate) and the 1.385 eV of every cycle heats the gas.
The cycle takes about a fifth of the input, cuts the Kr* population to a
third, and with it the Kr2* route (the VUV from 20 to 4.5 %), while the
248-nm share rises by a quarter as the harpoon and the trimer share what
is left. The near-IR 5p -> 5s band is 0.3 %: the measured relaxation
outruns the 29-ns radiation 120 to 1. With the relaxation band at its low
end (3e-11) nothing changes at 248 nm and the IR band is 1.4 %; with it
at 0, which no measurement supports, the rung could only radiate and 17 %
of the input left as near-IR light (CHANGELOG).

Every band's two ends under the two_t closure (`run_v0_4.py --bands
--eedf two_t`, the 248-nm share and its change against the first row):

```
         band   end      value   248nm  change    VUV   emit      g0   alpha  Tg_end
     as given              nan    5.79   +0.00    4.5   15.5   4.125   0.629     742
        kr2_f   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
        kr2_f  high      1e-09    7.59   +1.80    0.2   13.2   5.881   0.635     757
      kr2_pen   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kr2_pen  high      5e-10    5.80   +0.01    3.4   14.3   4.102   0.626     749
       kr2f_f   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
       kr2f_f  high    1.8e-10    5.93   +0.15    4.5   13.0   4.236   0.260     757
     kr2_form   low    1.7e-32    5.79   +0.00    4.5   15.5   4.125   0.629     742
     kr2_form  high    5.5e-32    5.31   -0.48   10.6   21.5   3.801   0.677     706
       kr2_f2   low    2.1e-10    5.79   +0.00    4.5   15.5   4.125   0.629     742
       kr2_f2  high      5e-10    5.81   +0.02    3.3   14.5   4.107   0.656     747
       kr_pen   low    1.5e-10    6.01   +0.22    5.6   17.0   4.340   0.653     732
       kr_pen  high      5e-10    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kr2f_f2   low    1.2e-10    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kr2f_f2  high    4.3e-10    5.96   +0.17    4.5   13.4   4.223   0.387     754
trimer_tg_pow   low      -3.32    5.79   +0.00    4.5   15.5   4.125   0.629     742
trimer_tg_pow  high       -0.5    2.29   -3.50    4.3   13.3   1.258   0.956     754
        krf_f   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
        krf_f  high    7.8e-10    3.99   -1.80    4.4   12.2   2.423   0.479     761
       kps_kr   low      3e-11    5.79   +0.00    5.5   17.8   4.094   1.053     727
       kps_kr  high    1.7e-10    5.79   +0.00    4.5   15.5   4.125   0.629     742
       k_5s5p   low          0    4.67   -1.12   20.0   30.5   2.562   0.802     649
       k_5s5p  high      9e-07    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kstep_s   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kstep_s  high    1.5e-07    5.13   -0.66    2.3   11.5   3.711   0.546     767
      kstep_p   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
      kstep_p  high      2e-07    5.56   -0.22    3.9   14.3   3.959   0.605     749
       k5p_f2   low          0    5.76   -0.02    4.4   15.3   4.109   0.622     743
       k5p_f2  high    1.3e-09    5.79   +0.00    4.5   15.5   4.125   0.629     742
     ii_yield   low        0.5    5.68   -0.11    4.5   15.3   4.043   0.623     743
     ii_yield  high          1    5.79   +0.00    4.5   15.5   4.125   0.629     742
       ii_f_c   low          0    5.79   +0.00    4.5   15.5   4.125   0.629     742
       ii_f_c  high        0.5    5.76   -0.03    4.5   15.4   4.105   0.628     742
     ib_scale   low        0.9    5.79   +0.00    4.5   15.5   4.125   0.629     742
     ib_scale  high        1.1    5.79   +0.00    4.5   15.5   4.125   0.629     742
```

The largest levers are the trimer's temperature law (-3.5 points at the
pure-Kr efficiency bound), KrF* + F at Kannari's analogue (-1.8), Kr2* + F
at its top (+1.8), the 5s <-> 5p fit off (-1.1) and the 5s stepwise
ionization at REBUILD's template (-0.7). The ion-ion yield and C share
move the line by 0.1 and 0.03 points under this drive; the wall loss of a
1-cm gap and the IB band move nothing here (no wave in the reference run).

## At Tahoma

ML's Tahoma configuration under the prototype beam drive (665 Torr with
2.8 Torr F2, 330 kW/cm^3 for 450 ns, n_e free from 1e10 cm^-3), the
200-400 ns window, against Spectra Technology 1989 Table 3-3 (beam stop
in): g0 6.39, alpha0 1.11, g0 - alpha0 5.28 %/cm. A score, never a fit.

| | measured | v0_3, two_t | v0_4, two_t | v0_4, ladder kannari | v0_4, kinetic |
|---|---|---|---|---|---|
| g0, %/cm | 6.39 | 5.68 (x0.89) | 5.86 (x0.92) | 5.87 (x0.92) | 5.91 (x0.93) |
| alpha0, %/cm | 1.11 | 0.84 (x0.75) | 0.82 (x0.74) | 0.84 (x0.76) | 0.82 (x0.74) |
| g0 - alpha0, %/cm | 5.28 | 4.85 (x0.92) | 5.04 (x0.95) | 5.03 (x0.95) | 5.09 (x0.96) |
| alpha0 by species | | Kr2F* 0.66, F- 0.06 | Kr2F* 0.60, F- 0.09, F2 0.06 | | Kr2F* 0.60, F- 0.09, F2 0.06, Kr** 0.03 |
| n_e, Te in the window | | 6.3e14, 1.60 eV | 4.5e14, 1.18 eV | | 4.4e14, 1.16 eV |

(v0_3's kinetic score was x0.88 on g0.) The kinetic closure moves v0_4's
Tahoma numbers by under 1 % from v0_2's closure, as it did v0_3's.

alpha0 is new in v0_4 and the first absorption LIGHT has scored. With
Kannari's compilation value for Kr2F* (6.0e-18) it came out 3.8 %/cm, 3.4
times the measurement, which sent the audit to the literature: Hakuta,
Komori, Mukai & Takuma measured (1.0 +- 0.2)e-18 at 248 nm by fluorescence
suppression and wrote that it is 'about 1/5 of the value assumed in the
current kinetic model'. With the measurement alpha0 is x0.74 of Tahoma's,
Kr2F* 73 % of it, and the net gain x0.95. The pure-Kr loss Tahoma measured
at 257.3 nm, 0.67 %/cm, LIGHT gives as 0.24 with 248-nm cross sections (a
rough check; Kr2+'s 1.6e-18 is itself an estimate). Under the beam the
5s <-> 5p cycle cools the electrons (1.18 against 1.60 eV) and lowers n_e,
and the gain rises 3 %.

## Found on the way

1. **A 5p rung that cannot relax radiates a sixth of the input.** The
   first build had KPS, the 5p + Kr relaxation, as a band at 0 (no Kr
   measurement was in hand), and the rung's only non-electron exit was the
   30-ns radiation: 17 % of the input left as near-IR light and the VUV
   fell to 14 %. Gazeli et al.'s measured 1.7e-10 (one 5p' level) and
   Whitehead et al.'s level-by-level constants (not retrieved) put the
   relaxation 120 times ahead of the radiation; the IR band is 0.3 %.
2. **The 5s <-> 5p cycle has to be in the kinetic equation.** With K76
   and K76r as LIGHT-side rates alone, the kinetic energy balance missed
   LIGHT's by 0.7 of the input at the reference point (the field carried
   half the power while the cycle took the rest unseen). Kannari's law
   has an exact cross-section equivalent (above), so the pair now sits in
   the Boltzmann equation, each operator scaled by its own population
   (the first build scaled the reverse by the Kr* fraction too, seven
   times too much here, and the gap only fell to 0.4); the field now
   carries 1.001 to 1.004 of the input through the pulse and the two
   balances agree to 7.9e-6 (the full reference run).
8. **The wall term stalled the solver.** Written with Te = U_e/(1.5 n_e)
   it has a Jacobian that explodes in the attaching afterglow, where n_e
   is a small difference of two large ion densities; a 0.3-us run took
   over 10 minutes. The flow is driven by U_e/(1.5 sum n+) instead,
   smooth and capped (`Model.wall_te`), and the lost particles' thermal
   energy leaves at Tg (E_wall_th), outside the Tg equation and counted
   once.
9. **scipy's event root-finder can fail in an attaching afterglow.** Where
   n_e falls to its floor a trigger's sign change between a step's ends
   is not reproduced by the dense output, and solve_ivp raises. The run
   now takes a short step without the triggers there and renews the
   table, as the stall guard does at a segment's start; the Kr* and Kr**
   fraction floors are 1e-7.
3. **Kr2F*'s cross section was the audit's largest single error.** The
   compilation value gave alpha0 3.4 times Tahoma's; the measured one
   gives x0.74 (above).
4. **The IB model's error band must act on the power.** Scaling nu_ei
   inside the Lorentz kernel is linear in the e-i share only for omega >>
   nu; at 1 THz in 1 atm Kr (nu/omega = 0.4) it misses by 10 % of itself
   (ib.py's own tests). ib_scale therefore scales the electron-ion share
   of the absorbed power under both closures.
5. **The solver's wave power is the Maxwellian kernel's, and the ions'
   part has two definitions.** With a Maxwellian handed to it, the
   solver's high-frequency field term gives the Lorentz-gas kernel's atom
   part to 0.7 % and, with the ions taken out of its momentum transfer,
   the kernel's MARGINAL electron-ion part (the power with the ions less
   the power without them) to 0.1 % (test_v0_4.py). That is not the
   kernel's printed "ei" share, nu_ei's share of the numerator: at 1 atm
   and 95 GHz nu_en is 4.4 omega, the wave is collisional, and the
   printed share is 1.5 times the marginal part (the two agree only for
   nu_en << omega). Until 2026-10-09 the test compared the solver's part
   with the printed share at a 50 % tolerance and this README said
   "within its grid"; the overview's own check found the 33 % and the
   definition behind it. The band ib_scale now acts on the marginal
   share under both closures (before, on the printed share under the
   Maxwellian ones), which changes nothing at ib_scale = 1. A full
   kinetic run with the wave closes its ledger with the wave's power
   read from the table.
6. **At 1 THz a kilovolt per centimetre heats nothing.** 3 kV/cm rms at 1
   THz in 1 atm Kr absorbs 1.7e-5 J/cm^3 in a microsecond (P ~ nu/omega^2),
   as the IB project's scans found; 1.5 kV/cm at 95 GHz absorbs 0.08
   J/cm^3 beside the 0.3 J/cm^3 given.
7. **The Langevin limb from the measured mobilities is 30 % below
   REBUILD's.** REBUILD's a scale / N is 8.2e-6 cm3/s amagat for both ion
   pairs; e (K0+ + K0-)/eps0 is 5.66e-6 (Kr+ + F-) and 6.35e-6 (Kr2+ + F-).
   The turnover rate at 1 amagat falls 10 to 14 %, and the 248-nm share at
   the reference point by 0.03 %.

## Limits

- The 5s <-> 5p rate is a fit, not a measurement, and it decides the
  partition. Jung et al. measured the cross sections out of the J = 2
  metastable into the ten 5p levels (PRL 94, 163202 (2005); PRA 73, 022722
  (2006)); their absolute values are not among the sources held. Band
  k_5s5p runs to 0 so the row's whole effect is bracketed (-1.1 points).
- KPS is one level's total loss (5p'[3/2]_2) carried to the manifold, its
  product unresolved; K27 is a search snippet of an abstract. Both are
  flagged in the rows and here.
- The stepwise ionizations (KSTEP_S, KSTEP_P) and the 5p harpoon are
  LIGHT-side rates: the Boltzmann equation does not see them (they are 0
  by default; with the preset on, the kinetic energy balance misses by
  their share). Only the 5s <-> 5p pair is a population process.
- Kr2+ + e makes the 5p rung entirely; Shiu & Biondi say 'mostly'.
- The absorbers are a diagnostic. A probe or a cavity would make F + e,
  Kr+ + Kr and Kr+ + e from the absorbed photons; neither is modelled.
  Kr2+'s and Kr**'s cross sections are an estimate and a calculation.
- The wall loss is the electropositive global-model closure (D_a = mu_i
  kTe/e, 2 Te per lost electron, Te per pair to the gas), stated and
  reported, with F- confined. In a 1-cm gap at the reference point it
  takes under 1e-3 of the input, as 1Dv0_3's resolved sheaths found.
- The wave's field is the one inside the plasma: no propagation,
  reflection or skin depth (ib.skin_hints prints n_e/n_crit and nu/omega).
  Under the two_t closure the Langdon factor is Turnbull's new one at Z =
  1; under the kinetic closure none, and the wave's momentum transfer
  carries the ions once with the IB logarithm.
- A kinetic run takes minutes, as in v0_3; `--eedf two_t` is seconds.

## Checks

`tests/test_v0_4.py` (the new physics), `tests/test_ib.py` (the IB module
against the IB project's own implementations, 219 checks),
`tests/test_network.py` (the rows, bands and limbs) and v0_3's tests
re-pointed, with the three whose expectations the new rows move rewritten
to the neutralised configuration. On 2026-10-09 the fast files give 315
passed and 1 strict expected failure (Hakuta's decays), `test_kinetic.py`
24 passed and 1 skipped (the slow convergence run); run the kinetic file
on its own. They hold:

- the ledger's closure in every new mode (ladder preset, KrF* + F, the
  ion-ion bands, a gap, the wave, all together), to 1e-9;
- the neutralised rung reproducing v0_3's reference point;
- the 5s <-> 5p pair's detailed balance (the electrons alone take Kr* and
  Kr** to g_p/g_s exp(-dE/Te) and pay the gap per net excitation), the
  pseudo cross section's Maxwellian rate against Kannari's law, and the
  rung's routing (C1 -> Kr**, KPS >> PRAD, the stepwise rows at 0);
- KrF* + F at the band's top outrunning the radiation after the burn-out;
- the ion-ion limbs (Thomson at low N, Langevin at high N, the yield and
  C share on the three product rows, the Langevin constant recomputed);
- the wall loss against its closed form, its ledger and its size;
- the optics against the hand sum; the wave as a second input with its
  band moving the electron-ion share alone, the model's kernel against
  v0's, the solver's field term against the Maxwellian kernel and a
  kinetic run with the wave closing;
- the te_fit and ionion rate kinds against their printed forms.

`tests/test_docs.py` holds the accuracy overview and the system map to
the model: every number they print is defined in a generated file, the
energy ladder is v0_4's, the overview's words are held to its numbers,
every citation has its reference and every reference is cited, ML's
Tahoma configuration is the one it names, and (`LIGHT_DOCS_FAST=1`,
`LIGHT_DOCS_SLOW=1`) the generated numbers match a fresh run.

**The propagated uncertainty** (the overview's Section 1, generated by
`make_overview_v0_4.py --unc`): every band at its two ends, the three
constants the model uses at a measured value with a quoted error that no
band covers (Kr2+ + e +-10 %, Kr+ + 2Kr +-9 %, the Kr2F* lifetime +-9 %,
through `Params.row_scale`) and the two cross sections with a quoted
error (sigma_se +-12 %, sigma(Kr2F*) +-20 %, analytic), each alone, on
the 248-nm share at the reference point and on Tahoma's g0, alpha0 and
g0 - alpha0 under v0_2's closure; half the difference between the two
ends is the term, combined in quadrature (indicative: a band's range is
not a standard deviation) and as an envelope (the bound). Not propagated:
the electron-impact cross sections (no quoted error) and the fitted
KrF(B,C) set (pinned by the training decays, scored instead).

## Files

```
v0_4/
  README.md, SOURCES.md    this file; the sources (section 7 is v0_4's)
  species.py     + Kr** (5p, 11.30 eV)
  network.py     + ladder_rows, f_quench_rows, ION_K0, langevin_k, the
                 te_fit and ionion rate kinds, Row.frac, the new bands,
                 LADDER_KANNARI, FREE_BANDS, bands_read
  model.py       + Params.ladder, gap_cm, em_E_V_cm, em_f_Hz, row_scale
                 (a sensitivity factor on one row's rate); wall_loss,
                 em_power, the E_em / E_wall / E_ion_heat ledger slots,
                 the population process and the wave in the kinetic state
  boltzmann.py   + PopProcess, the ions' Coulomb collisions, the
                 high-frequency field term, em_power
  ib.py          Turnbull et al.'s inverse bremsstrahlung and the
                 Lorentz-gas kernel (SI inside)
  optics.py      g0 and alpha at 248 nm
  report.py      + the optics block, the walls line, the EM input
  run.py         + --ladder, --gap, --em-E, --em-f
  tests/         test_v0_4.py, test_ib.py, test_docs.py and v0_3's,
                 re-pointed
docs/
  accuracy_overview_v0_4.tex/.pdf   how accurate v0_4 is (this README's
                 numbers, generated: overview_v0_4_fast/slow/unc.tex by
                 make_overview_v0_4.py, the map's numbers, the ladder),
                 with the errors, the propagated uncertainty and the
                 references
  system_map_v0_4.tex/.pdf          the energy flow with the 5p rung, the
                 bands' table and the optics (make_system_map_v0_4.py)
  energy_ladder_v0_4.tex            make_energy_ladder.py v0_4
  krf_data/tahoma_v0_4.py           Tahoma's g0 and alpha0 for any version
```
