# Reheating the KrF afterglow by inverse bremsstrahlung: the compendium in brief

The quantum derivation in short, the bottom line, and the efficiency that a kinetic model of the laser
gas gives. Summary of *Inverse-bremsstrahlung heating of KrF laser gases* (the compendium, 163 pages),
20 October 2026. The PDF is `compendium_summary.pdf`.

## Bottom line

- **Reheating is possible, and quantum mechanics says so.** The electrons of the afterglow absorb a THz
  wave, and at the rate of classical physics. The quantum rate of free-free absorption between exact
  scattering states reduces, for a photon that is soft against the electron's energy, to the
  momentum-transfer cross section, and the classical absorption built on it is right to **1.6 %** from
  0.3 to 4 THz for electrons from 300 K to 1.08 eV. Stimulated emission takes 0.2 % at the afterglow's
  0.63 eV and 0.67 THz, the exchange of several photons changes no rate, and the nonlinearities of
  Langdon and Silin need 33 times the heating intensity and more. In the afterglow, 10^20 m^-3 at
  0.63 eV, a 1.3-THz wave is absorbed through 10 cm of gas, and 1.4×10^3 to 1.5×10^4 W/cm^2 holds the
  electrons at their temperature.
- **It does not make the laser meaningfully more efficient.** A time-dependent kinetic model of the
  Kr/F2 laser gas, pumped with 300 kW/cm^3 for 1 µs at 1 atm with 0.5 % F2, which turns 5.73 % of the
  pump's energy into 248-nm light, was run with and without the compendium's 3-MW/cm^2 wave in the
  afterglow. With every rate at its default the wave *lowers* the efficiency, by 0.037 percentage
  points: each joule it gives the gas returns 4.5 % as 248-nm light, less than the pump's 5.73 %. It
  raises the efficiency only if Kr2* + F → KrF(B) + Kr, which no one has measured, runs faster than
  about 2×10^-12 cm^3/s, and then by at most 0.14 points, about 2 % of the efficiency, with a wave
  source at least 0.49 times as efficient as the pump's. Across that rate's band the change runs from
  **−0.04 to +0.14 points**, and the model's other uncertain inputs widen it to −2.57 to +0.35. Weaker
  waves lose 248-nm light, and with 2 % F2 the afterglow's electrons last 8.5 ns, which leaves
  −0.0014 points.
- **What would decide it:** a measurement of Kr2* + F → KrF(B) + Kr; the electrons' energy distribution
  under the THz field, on which the model's two closures disagree 24-fold in the absorbed energy; the
  stepwise ionization of Kr* and the electron quench of KrF*; and a model with a cavity, since the
  afterglow's 248-nm light is counted here as fluorescence.

## 1. The question

An electron beam or a field pumps a KrF laser gas for a microsecond. When the pump stops, the electrons
it leaves behind cool and are lost to attachment and recombination: the afterglow. Reheating sends an
electromagnetic wave into the afterglow, the wave heats the free electrons by inverse bremsstrahlung
(IB), and the hotter electrons excite krypton again, which forms KrF* while F2 lasts. The compendium
asks whether that can make the laser more efficient. Its Part I splits the question into five links:
the wave is absorbed through the gain medium (coupling), its absorption is known (absorption), it heats
the electrons (heating), the heat forms KrF* (partition), and the KrF* becomes laser light that pays
for the wave's source (laser and wall plug). The first three are physics that quantum mechanics
settles, Section 2. The last two need the kinetics of the laser gas, Section 3.

## 2. Why the electrons absorb: the derivation in short

Gaussian units, as in the compendium. An electron of momentum ħ**k** and energy E = ħ²k²/2m scatters
on an atom of potential V(r) in a wave E0 cos ωt polarized along **ê**, with u = ħω/k_BT and N atoms
and n_e electrons per unit volume.

**(1) The golden rule.** The rate of absorbing one photon between exact scattering states, E → E + ħω,
divided by the photon flux cE0²/8πħω and the electron flux, is the free-free cross section
(Chapter 9, Eq. 9.14)

$$\sigma_+=\frac{32\pi^2}{3}\,\frac{e^2m^2\omega}{\hbar^3c}\,\frac{S}{k^3k'},\qquad
S=\sum_{l}\big[(l+1)R_{l\to l+1}^2+l\,R_{l\to l-1}^2\big],\qquad
R_{l\to l'}=\frac{1}{m\omega^2}\int_0^\infty u_{k'l'}\,V'(r)\,u_{kl}\,dr,$$

with u_kl the radial functions and ħk' = (ħ²k² + 2mħω)^(1/2). The last form, the acceleration form,
shows that the photon is absorbed where the atom accelerates the electron.

**(2) Emission.** Stimulated emission is the reverse transition, and detailed balance,
k²σ+(k→k') = k'²σ−(k'→k), leaves Maxwell electrons, with g(E) their distribution of energies, the
absorption coefficient (Eqs. 9.17 and 9.20)

$$\kappa=n_eN\,(1-e^{-u})\int_0^\infty g(E)\,v\,\sigma_+(E)\,dE .$$

**(3) Soft photons.** As ħω/E → 0 the force identity ∫ u_{k,l+1} V' u_kl dr =
(ħ²k²/2m) sin(δ_l − δ_{l+1}) (Eq. 9.31) turns S into the phase shifts, and σ+ into the
momentum-transfer cross section σ_m, the low-energy theorem of Low [1] (Eq. 9.34):

$$\sigma_+\to\frac{4\pi}{3}\,\frac{e^2v^2}{c\hbar\omega^3}\,\sigma_m(E),\qquad
\sigma_m=\frac{4\pi}{k^2}\sum_l(l+1)\sin^2(\delta_l-\delta_{l+1}).$$

**(4) The classical limit.** With Eq. (3) and 1 − e^(−u) → u, Eq. (2) becomes the Drude absorption of a
Lorentz gas (Eq. 9.37),

$$\kappa_{\rm cl}=\frac{\omega_p^2}{c\,\omega^2}\,\nu_{\rm eff},\qquad
\nu_{\rm eff}=\frac{mN}{3k_BT}\,\big\langle v^3\sigma_m\big\rangle,\qquad
\omega_p^2=\frac{4\pi n_ee^2}{m},$$

the energy-weighted collision frequency of the compendium's Eq. 4.12. The quantum and the classical
absorption differ by the ratio (Eq. 9.39) R = κ/κ_cl = [(1 − e^(−u))/u] ⟨v³σ_m G⟩/⟨v³σ_m⟩, with G the
free-free Gaunt factor of the atom. Computed from the exact phase shifts of model potentials for
neon, argon and krypton, R lies within 1.6 % of one from 0.3 to 4 THz at every electron temperature
from 300 K to 1.08 eV in both fills (Sec. 9.9). Only far above that band does it matter: at 28.3 THz
and 300 K the quantum absorption is 10 % below the classical one in the neon fill and 28 % below it in
the argon fill.

**(5) Ions.** The same golden rule on the Coulomb potential gives the collision frequency ν_T Λ with
Λ = (π/√3)⟨g⟩(1 − e^(−u))/u (Eq. 4.15), which in the Born approximation is ln(4k_BT/e^γ ħω) for
u ≪ 1: the logarithm needs no cutoff put in by hand (Sec. 8.13). The ions matter only in a thermalized
afterglow. At 300 K and 10^20 m^-3 they absorb 1.4 times as much as the atoms in the neon fill at
1 THz, but while the electrons stay at 0.63 eV the ions catch up with the atoms only above
10^22 m^-3 (Sec. 9.10).

**(6) Strong fields.** For electrons whose energy far exceeds the photon's, the cross section of Kroll
and Watson [2] for exchanging ℓ photons heats at Σ_ℓ ℓ² J_ℓ²(x) = x²/2, the one-photon rate, at any
intensity (Sec. 9.7). The distribution's own nonlinearities need far more than heating does:
Langdon's flattening [3] 33 to 51 times the heating intensity of the afterglow, and Silin's regime a
further 11.5 times (Sec. 6.5).

**(7) Heating.** Each electron gains P = e²E_rms² ν_eff/[m(ω² + ν_eff²)] (Eq. 4.7), and the field sets
the electron temperature once P exceeds the electrons' collisional loss P_loss, at the intensity
I_heat = P_loss/σ_abs with σ_abs = 4πe²ν_eff/mcω² for ω ≫ ν_eff (Eq. 6.3).

**In the afterglow.** In the neon-buffered fill at 0.63 eV and 10^20 m^-3 (Sec. 6.5, Table 6.2), a
0.67-THz photon has u = 0.0044, so stimulated emission takes 1 − (1 − e^(−u))/u = 0.2 %, and the
electrons' collision frequency is 0.12 of the wave's angular frequency. The wave is absorbed in
2.5 mm at 170 GHz and in 2.6 cm at 0.67 THz, so a depth of 10 cm takes 1.3 THz (Sec. 2.1). The
electrons lose 3.6×10^8 eV/s each, and I_heat is 1.4×10^3 W/cm^2 at 170 GHz and 1.5×10^4 W/cm^2 at
0.67 THz, within reach of gyrotrons.

**So:** quantum mechanics does not stand in the way. The afterglow's electrons absorb a THz wave at the
classical rate, to 1.6 %, and the wave heats them at 10^3 to 10^4 W/cm^2.

## 3. How much it raises the efficiency

**The formation chain overstates it.** Counted with the branchings of a KrF kinetic network, a field far
above the collision frequency forms 0.274 to 0.352 KrF* per absorbed joule at the afterglow's density,
1.4 to 1.8 times the electron beam's 0.198 (Secs. 2.5 and 6.8 of the compendium). That count assumes
F2 for every Kr* and leaves out the electron collisions of the excited atoms. A kinetic model that
carries both answers differently.

**The model.** It is a time-dependent model of one uniform box of Kr/F2 at constant volume, with 59
reactions on 12 species and an energy ledger that closes. The electrons' energy distribution is solved
from the two-term Boltzmann equation with electron-electron collisions [4] and the wave's own
high-frequency field term, so that the absorbed power is the solver's own. The atoms absorb through the
classical kernel of Eq. (4) and the ions through the IB model of Turnbull et al. [5, 6]. Nothing in it is
tuned to a device. Against the Tahoma amplifier it gives the net small-signal gain to −2 %, the
saturation intensity to +4 % and the formation efficiency to −3 % (0.261 against 0.27), and its inputs'
propagated uncertainty is ±40 %.

**The test.** The gas, 1 atm of Kr with 0.5 % F2 at 300 K, is pumped with 300 kW/cm^3 put into its
electrons for 1 µs and followed to 2 µs, once without and once with the compendium's wave in the
afterglow: 47.6 kV/cm in amplitude, 3 MW/cm^2, at 6.2 THz (Sec. 4.9). Without the wave the gas turns
5.73 % of the pump's energy into 248-nm light. The *marginal* is the 248-nm energy that each absorbed
joule of the wave adds, and the efficiency counted with the wave as an input rises only where the
marginal beats the pump's 5.73 %. When the wave opens, the gas holds 2.7×10^14 cm^-3 electrons, at
0.60 eV on average over the afterglow.

Table 1. The twin pairs with every rate at its default. The change is that of the 248-nm efficiency
with the wave counted as an input. The break-even source is how efficient the wave's source, from
electricity to absorbed power, must be relative to the pump's for the pair to break even at the wall
plug: the pump's own efficiency over the marginal.

| gas; wave; window | absorbed, % of pump | marginal | pump's own | change, points | break-even source |
|---|---|---|---|---|---|
| 1 atm, 0.5 % F2; 3 MW/cm^2; 1–2 µs | 3.1 | 4.5 % | 5.73 % | −0.037 | 1.27 |
| 1 atm, 0.5 % F2; 3 MW/cm^2; 1–1.1 µs | 1.2 | 3.5 % | 5.73 % | −0.026 | 1.62 |
| 1 atm, 0.5 % F2; 0.3 MW/cm^2; 1–2 µs | 0.055 | 0.81 % | 5.73 % | −0.0027 | 7.1 |
| 1 atm, 0.5 % F2; 0.03 MW/cm^2; 1–2 µs | 0.004 | −0.79 % | 5.73 % | −0.0003 | never |
| 0.5 atm, 2 % F2; 3 MW/cm^2 at 6.6 THz; 1–2 µs | 0.074 | 13.9 % | 15.83 % | −0.0014 | 1.14 |

With 2 % F2 at 0.5 atm the electrons last 8.5 ns after the pump, against 42 ns at 1 atm with 0.5 % F2
(their density integrated over the afterglow, over its value when the wave opens), so the wave finds
few electrons to heat.

**Where the wave's energy goes.** At 1 atm with the strong wave, 4.5 % of the absorbed energy leaves as
248-nm light, 29.5 % as the 147-nm light of Kr2*, 2.3 % in other bands, 52.8 % as heat, and 10.9 % stays
in the gas as F atoms, Kr* and Kr2*. The electrons warm from 0.60 to 1.06 eV, and their
time-integrated density grows 2.6-fold. But the pump has burned the F2 below 10 % of its fill by
815 ns, before the wave opens, so the Kr* the wave makes finds no F2 to harpoon. It forms Kr2*, which
radiates at 147 nm, and more than half of the wave's energy ends as heat. The formation chain counts
neither the burned F2 nor these losses.

**One unmeasured rate decides the sign.** Kr2* + F → KrF(B) + Kr would hand that Kr2* to the F atoms
the burn left. It is 0 in the model by default, the rule for a rate that no one has measured. Kannari
et al. [7] use 3×10^-10 cm^3/s, an estimate carried over from Ar2* + F. Across its band:

Table 2. The strong wave at 1 atm and 0.5 % F2 over the whole afterglow, with Kr2* + F → KrF(B) + Kr
across its band. The rate also raises the pump's own efficiency, since the pump's afterglow holds the
same Kr2*.

| Kr2* + F → KrF(B) + Kr, cm^3/s | marginal | pump's own | change, points | break-even source | wave to 147 nm |
|---|---|---|---|---|---|
| 0 (the default) | 4.5 % | 5.73 % | −0.037 | 1.27 | 29.5 % |
| 1×10^-12 | 5.2 % | 5.80 % | −0.018 | 1.11 | 28.0 % |
| 3×10^-12 | 6.4 % | 5.93 % | +0.014 | 0.92 | 25.4 % |
| 1×10^-11 | 9.2 % | 6.25 % | +0.077 | 0.68 | 18.9 % |
| 3×10^-11 | 12.3 % | 6.67 % | +0.127 | 0.54 | 10.5 % |
| 1×10^-10 | 14.5 % | 7.14 % | +0.138 | 0.49 | 4.0 % |
| 3×10^-10 (Kannari et al.) | 15.3 % | 7.45 % | +0.131 | 0.49 | 1.4 % |
| 1×10^-9 (the band's top) | 15.6 % | 7.67 % | +0.123 | 0.49 | 0.4 % |

The sign turns at about 2×10^-12 cm^3/s, under 1 % of Kannari's value. Above it the reheat gains at
most 0.138 points of 7.14 %, about 2 % of the efficiency, and only with a wave source at least 0.49
times as efficient as the pump's.

**The other inputs.** Each of the model's 20 uncertain inputs at its two ends, one at a time, gives the
change at the defaults as −0.037 ± 0.84 points in quadrature, and the envelope that stacks every
input's excursions each way runs from −2.57 to +0.35 points. The measured 5s→5p cross section at the
lower end of its band alone gives −1.55 points.

**The other closure.** With a Maxwellian closure of two temperatures in place of the solved
distribution, the same wave sustains a discharge through the afterglow: the electrons reach 2.22 eV,
the wave absorbs 75 % of the pump's energy, 85 % of it ends as heat, and the efficiency falls by 2.25
points.

**Limits.** The model is 0-D: the field inside the gas is given, and the gas 10 cm deep, the reflection
and the absorption length as n_e falls are not modelled. It has no cavity. The 248-nm energy is
fluorescence, and the afterglow's small-signal net gain rises with the wave from 0.12 to 0.20 %/cm,
against 3.1 %/cm during the pulse, so there is little to extract after the pump either way. Waves
stronger than 3 MW/cm^2 were not run with the solved distribution. The pure-krypton model of the same
family has not been run with the wave.

**So:** the reheat changes the 248-nm efficiency by −0.037 points at the model's defaults, by at most
+0.14 points if the unmeasured Kr2* + F channel is fast, and by −0.0014 points with 2 % F2. A gain of
at most about 2 % of the efficiency, conditional on an unmeasured rate and on a THz source half as
efficient as the pump, is not a reason to build the reheat.

## 4. What would settle it

1. **Kr2* + F → KrF(B) + Kr.** From about 1×10^-11 cm^3/s up the reheat pays at the model's reference
   point with a source 0.68 times as efficient as the pump's or better, and at 1×10^-12 and below it
   does not pay at all (Table 2).
2. **The electrons' distribution under the THz field.** The model's two closures put the wave in
   different regimes and differ 24-fold in the energy it absorbs. A measurement, or a kinetic treatment
   that follows the distribution as the wave turns on, would settle it.
3. **The electron collisions of the excited atoms:** the stepwise ionization e + Kr* → Kr+ + 2e and the
   electron quench of KrF*, which the compendium names as deciding the sign of a field's gain
   (Sec. 6.8). A THz pulse that sets the electron temperature and reads it back by its transmission can
   measure both.
4. **A cavity**, to count the extraction of the afterglow's KrF*, and the wall plug with a real THz
   source at the frequency of Section 2.

## References

1. F. E. Low, Phys. Rev. **110**, 974 (1958), doi:10.1103/PhysRev.110.974.
2. N. M. Kroll and K. M. Watson, Phys. Rev. A **8**, 804 (1973), doi:10.1103/PhysRevA.8.804.
3. A. B. Langdon, Phys. Rev. Lett. **44**, 575 (1980), doi:10.1103/PhysRevLett.44.575.
4. G. J. M. Hagelaar and L. C. Pitchford, Plasma Sources Sci. Technol. **14**, 722 (2005),
   doi:10.1088/0963-0252/14/4/011.
5. D. Turnbull et al., Phys. Rev. Lett. **130**, 145103 (2023), doi:10.1103/PhysRevLett.130.145103.
6. D. Turnbull et al., Phys. Plasmas **31**, 063304 (2024), doi:10.1063/5.0203446.
7. F. Kannari, M. Obara and T. Fujioka, J. Appl. Phys. **57**, 4309 (1985), doi:10.1063/1.334590.
