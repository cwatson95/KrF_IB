# Building the slab in COMSOL 6.4

Two ways to build the same model. Both read the files in `KrF/` (Kr + 0.5 % F2) or `Kr/` (pure Kr),
and both follow the gas for 1 µs of afterglow from the state LIGHT's kinetic model reaches at the end
of the 1-µs pump.

- **A. Plasma Module + RF Module** (*Microwave Plasma*). The electrons, the ions and the excited
  states are solved by the Plasma interface with the chemistry imported from a file. The wave is
  solved by Electromagnetic Waves, Frequency Domain. The two are coupled by COMSOL's own plasma
  conductivity and electron heat source.
- **B. RF Module + Domain ODEs**. This is the same physics with the chemistry written as ordinary
  differential equations at every mesh point. It is simpler, and it is the exact system the Python
  reference solves, so use it if A will not converge or to check A.

Below, `<gas>` stands for `KrF` or `Kr`, and paths are relative to this folder. Every COMSOL name used
here comes from COMSOL's 6.4 documentation. COMSOL has not been run on these steps, so check the
**Checks** after each part. In what follows, *Gas* means the gas domain, and *Front* and *Back* mean
the outer boundaries at x = −L_vac and x = L_gas + L_vac.

## The physics in one paragraph

A plane wave (out-of-plane E, incident intensity `I0`) is launched from vacuum into gas `L_gas` deep.
The gas is a strip `h_strip` high, whose long sides are perfect magnetic conductors for the wave and
insulating for the plasma, so the problem is one-dimensional. The plasma's permittivity is the Drude
form, with LIGHT's kinetic absorption as its collision frequency `nu_e`. The electrons absorb the
wave, and the absorbed power heats them. Their collisions, rates and energy losses are LIGHT's own,
read from `kinetics_<gas>.txt`. The reaction heat and the elastic power heat the gas (`Tg`), and the
excimers radiate (`P248`, `P147`). Two runs make a twin pair: `I0_MW = 0` (no wave) and
`I0_MW > 0`. Their difference over the energy absorbed is the efficiency the wave buys.

## A. Plasma Module + RF Module

### A1. Model Wizard
1. **File > New > Model Wizard > 2D**.
2. Add **Plasma > Microwave Plasma**. This adds Plasma (plas), Electromagnetic Waves, Frequency
   Domain (emw), and the multiphysics nodes Plasma Conductivity Coupling and Electron Heat Source.
3. Add **Mathematics > ODE and DAE Interfaces > Domain ODEs and DAEs** twice:
   - The first holds the gas temperature. Field name `Tg`, one dependent variable `Tg`.
   - The second holds the accumulated energies. Field name `E`, with dependent variables `E248`,
     `E147` and `Ewave` (KrF), or `E147` and `Ewave` (Kr).
4. Study: **Preset Studies for Selected Multiphysics > Frequency-Transient**.

### A2. Global Definitions
1. **Parameters 1 > Load from File**: `<gas>/params_<gas>.txt`.
2. **Functions > Interpolation**, the closure:
   - Data source: **File**, `<gas>/kinetics_<gas>.txt`.
   - Data format: **Grid**, **Number of arguments: 3**.
   - Function name: `kin_<gas>`, position in file: 1.
   - Interpolation: **Linear**. Extrapolation: **Constant**.
   - Leave the arguments and the function without units.
3. (Kr only) A second **Interpolation** from `Kr/heavy_Kr.txt`:
   - Data format: **Grid**, **2 arguments**.
   - Function name: `heavy_Kr`, position 1.
   - Linear interpolation, Constant extrapolation.

### A3. Geometry
- Three rectangles, all of height `h_strip`, with their lower-left corners at:
  - `(-L_vac, 0)`, width `L_vac`: the vacuum in front;
  - `(0, 0)`, width `L_gas`: the gas;
  - `(L_gas, 0)`, width `L_vac`: the vacuum behind.
- **Form Union**.
- In the gas rectangle's settings, tick **Resulting objects selection** (the gas selection, *Gas*).
- Make explicit boundary selections for *Front*, *Back*, and *Sides* (every horizontal boundary).

### A4. Definitions (Component 1)
1. **Variables 1**:
   - Geometric entity level: **Domain**, selection *Gas*.
   - **Load from File**: `<gas>/variables_<gas>_plasma.txt`.

   Some names in it (`plas.*`, `emw.Qrh`, `Tg`) stay unresolved until the physics below exists. That
   is expected.
2. **Nonlocal Couplings > Integration**: name `intop_gas`, selection *Gas*.

### A5. Materials
- **Blank Material** on all three domains, with:
  - Relative permittivity 1;
  - Relative permeability 1;
  - Electrical conductivity 0.
- In the gas, the Plasma Conductivity Coupling supplies the plasma's conductivity.

### A6. Plasma (plas)
1. Selection: *Gas* only.
2. In the interface settings, tick **Use reduced electron transport properties**.
3. **Plasma Model 1**:
   - Mean electron energy: **Local energy approximation** (the default).
   - Electron transport properties: **Specify mobility only**, with Reduced electron mobility
     μeNn = `muN_e`.
   - Electron energy distribution function: leave it at **Maxwellian**. It is not used, because
     every rate comes from `kin_<gas>`.
   - Model Inputs:
     - Temperature T = `Tg_`;
     - Absolute pressure pA = `N0*k_B_const*Tg_`. The heavy-particle density then stays `N0` as
       the gas heats (constant volume, as in LIGHT).
4. Right-click Plasma, then **Plasma Chemistry Import**:
   - File: `<gas>/chemistry_<gas>.txt`, then **Import**.
   - On 6.1 to 6.3, use the Plasma Chemistry add-in with the same file.

   This creates the Species, the Electron Impact Reactions and the Reactions.
5. **Check the import.**
   - Open the first Electron Impact Reaction. It should show:
     - Collision type **Excitation**;
     - Energy loss **0**;
     - **Specify reaction using: Rate constant**;
     - kf = `kf_G1` (KrF) or `kf_G1_1s5` (Kr), exactly. If the field reads `kf_G1*N_A_const` or
       anything more than the variable's name, delete the extra factor in every reaction. The
       variable already carries N_A and its units, m³/(s·mol).
   - Open one Reaction. Its rate constant should be `kf_<name>` exactly.
6. **Species**:
   - `Kr`: **From mass constraint**.
   - `Kr2+`: **Initial value from electroneutrality constraint**.
   - Every other species: initial mass fraction `w0_<tag>`. The tag is the species name with `+` →
     `_1p` and `-` → `_1m`, so `w0_Kr_1p` and `w0_F_1m`; the other tags are their names (`w0_Krs`,
     `w0_KrFB`, …).
7. **Initial Values 1**: electron density `ne0`, mean electron energy `ebar0`, electric potential 0.
8. **General Power Deposition**: selection *Gas*, General heat source Q = `Qe`. This one term carries
   every electron energy gain and loss of LIGHT's ledger, wave excepted:
   - the inelastic losses;
   - the superelastic and Penning gains;
   - the captured electrons' energy;
   - the elastic power.

   That is why every Electron Impact Reaction has an energy loss of 0.
9. Boundaries:
   - Keep the defaults (Insulation, Zero Charge) on every boundary of the gas.
   - Add **Ground** on the gas's front boundary (x = 0) only, to fix the potential's reference.
   - Add no Wall: the windows are 0.1 mm of plasma away from mattering in 1 µs at 1 atm.

### A7. Electromagnetic Waves, Frequency Domain (emw)
1. Settings: Electric field components solved for: **Out-of-plane vector**.
2. **Perfect Magnetic Conductor** on *Sides*.
3. **Scattering Boundary Condition 1** on *Front*:
   - Incident field: **Wave given by E field**;
   - Incident electric field E0 = (0, 0, `E0`);
   - leave the wave direction at its default, opposite to the boundary normal, so into the domain.
4. **Scattering Boundary Condition 2** on *Back*: no incident field. The transmitted wave leaves
   and is lost. For a mirror behind the gas, use **Perfect Electric Conductor** on *Back* instead.

### A8. Multiphysics
- **Plasma Conductivity Coupling 1**: Electron mobility (the default). COMSOL takes the collision
  frequency from the mobility, ν = e/(m_e μ_e), which is `nu_e`.
- **Electron Heat Source 1**: the default.

### A9. Domain ODEs
1. `Tg`:
   - selection *Gas*;
   - Units: dependent variable **K**, source term **K/s**;
   - Distributed ODE: source term f = `dTg`, damping coefficient 1;
   - Initial value `Tg0`.
2. `E`:
   - selection *Gas*;
   - Units: dependent variables **J/m^3**, source terms **W/m^3**;
   - sources:
     - KrF: `P248`, `P147`, `Qw`;
     - Kr: `P147`, `Qw`;
   - initial values 0.

### A10. Mesh
**Mapped**, with:
- Distribution: on every vertical boundary, **1 element**;
- Size: maximum element size `hmax`, which is λ0/6 (λ0 = 0.23 mm at 1.3 THz, so about 2,200
  elements along 10 cm).

### A11. Study
1. **Frequency-Transient**: frequency `f0`, output times `range(0,5[ns],t_end)`.
2. **Solver Configurations > Time-Dependent Solver**:
   - Relative tolerance **1e-5**;
   - Initial step **1e-15 s**: the electrons cool from 2.4 eV to 0.9 eV in the first nanosecond.
3. **Parametric Sweep** over `I0_MW` with values `0 3` (the twins). Add `0.3 0.03` for the scan the
   reference reports.

### A12. Results (at t = t_end)
Use a **Global Evaluation** for each sweep value.
- Energies per area [J/cm²]:
  - absorbed wave: `intop_gas(Ewave)/h_strip`;
  - 248 nm in the afterglow (KrF): `intop_gas(E248)/h_strip`;
  - 147 nm in the afterglow: `intop_gas(E147)/h_strip`.
- Efficiency of the main band, pump and afterglow together:
  - KrF: `(Epump_248 + intop_gas(E248)/h_strip)/(E_given + intop_gas(Ewave)/h_strip)`;
  - Kr: `(Epump_147 + intop_gas(E147)/h_strip)/(E_given + intop_gas(Ewave)/h_strip)`.
- The pair's numbers:
  - the change d_eff = eff(on) − eff(off);
  - the marginal efficiency = (E_band,on − E_band,off) / E_wave;
  - the source efficiency needed to break even = eff(off) / marginal.
- Compare with `reference/results.md`.

### Checks
1. At t = 0, `plas.ne` = `ne0` and `Te_eV` = `Te0`, and the number densities equal `n0_<tag>`.
2. In the **off twin** (`I0_MW = 0`), against `<gas>/reference_<gas>.csv`:
   - n_e(t) should fall by attachment and recombination;
   - T_e(t) should settle at 0.45 to 0.8 eV;
   - Tg should rise by about 20 K.

   A difference by a factor near 6×10²³ in the rates means N_A was applied twice (A6 step 5).
3. With a thin gas (`L_gas = 0.1[mm]`), `f0 = 6.2[THz]` and `I0_MW = 3`, the twin should reproduce
   the reference's 0-D pair, the kit's closure row in `reference/results.md`.

## B. RF Module + Domain ODEs (no Plasma interface)

The same geometry (A3), parameters, functions (A2), materials and mesh. Then:

1. Variables: `<gas>/variables_<gas>_ode.txt` on *Gas*.
2. **Domain ODEs and DAEs**, three of them on *Gas*:
   - Species, field `n`:
     - one dependent variable per species, named `n<tag>`: `nKr`, `nF2`, `nKrs`, `nKr_1p`,
       `nF_1m`, …;
     - units 1/m^3, source units 1/(m^3*s);
     - sources `dn_<tag>`;
     - initial values `n0_<tag>`.
   - `Tg`, as in A9: source `dTg`, initial value `Tg0`.
   - Energies, field `U`:
     - KrF: `Ue`, `E248`, `E147` and `Ewave`; Kr: `Ue`, `E147` and `Ewave`;
     - units J/m^3, sources W/m^3;
     - sources `Qe+Qw`, then `P248`, `P147` and `Qw` (KrF), or `P147` and `Qw` (Kr);
     - initial values `Ue0`, 0, 0, 0.
3. emw:
   - out-of-plane field;
   - in the gas, the material's relative permittivity is `epsr_D` (complex);
   - permeability 1 and conductivity 0 everywhere;
   - boundaries as in A7.
4. Study: Frequency-Transient as in A11.
5. Results as in A12.
