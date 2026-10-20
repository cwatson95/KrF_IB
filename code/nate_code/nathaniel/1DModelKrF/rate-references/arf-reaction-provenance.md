# ArF network provenance — `data/reactions/arf.toml` vs Chen (2025)

Every reaction in the working ArF network, its adopted rate/dataset, numbered citations
(list at bottom), and whether the value matches **Chen et al., Phys. Plasmas 32, 113502
(2025), Table III** as printed. Where we deviate, the note gives the one-line reason; full
forensic detail lives in `bug-fixes.md` and `changelog.md`. Rates in SI (m³/s, m⁶/s);
electron-impact channels are cross-section datasets resolved through the two-term Boltzmann
solve (Chen used BOLSIG+ over the same SIGLO source).

| ID | Reaction | Rate / dataset | Cites | = Chen? | Difference note |
|---|---|---|---|---|---|
| R1 | e + Ne → Ne\* + e | `e_Ne_excitation_1` (SIGLO 4s-manifold sum) | [1][3] | ~ | Chen cites one lumped level; we sum the full SIGLO 3P2+3P1+1P1 manifold (his source's own content) |
| R2 | e + Ne → Ne\*\* + e | `e_Ne_excitation_2` (SIGLO higher-manifold sum) | [1][3] | ~ | Same manifold-sum treatment (2P+2S3D+3P) |
| R3 | e + Ne → Ne⁺ + 2e | `e_Ne_ionization` | [1][3] | ✓ | — |
| R4 | e + Ne\* → Ne⁺ + 2e | `e_Nes_ionization` (SIGLO Ne(3P2) set) | [1][3] | ~ | Chen prints no threshold; ours estimated 4.95 eV; SIGLO record cross-validated vs TRINITI ≤2% |
| R5 | e + Ar → Ar\* + e | `e_Ar_excitation_1` = SIGLO total × (1−f_BSR(E)) | [1][2][3] | ✗ | Chen's R5/R6 split is unreproducible from SIGLO (one lumped level). We partition the SIGLO total by the BSR ground→4s/4p branching f(E)≈0.25–0.36; R5+R6 ≡ SIGLO exactly |
| R6 | e + Ar → Ar\*\* + e | `e_Ar_excitation_2` = SIGLO total × f_BSR(E) | [1][2][3] | ✗ | Counterpart of R5; threshold 13.074 eV (lowest BSR 4p). Chen's source carries no such channel |
| R7 | e + Ar → Ar⁺ + 2e | `e_Ar_ionization` | [1][3] | ✓ | — |
| R8 | e + Ar\* → Ar⁺ + 2e | `e_Ars_ionization` (TRINITI) | [16][3] | ~ | SIGLO has no excited-Ar ionization; TRINITI is the only measured set in any local database. Threshold estimated 4.21 eV |
| R9 | e + Ar\*\* → Ar⁺ + 2e | Tₑ-fit 1.56×10⁻¹³·Tₑ^0.71·e^(−2.63/Tₑ) | [3][5] | ✓ | Value identical; provenance note: the fit is R&K's Xe\*\* row transplanted by Chen (documented analogy) |
| R10 | e + F₂ → F⁻ + F | `e_F2_attachment` | [1][3] | ✓ | — |
| R11 | e + F₂ → F₂(a3) + e | `e_F2_excitation_a3` | [1][3] | ✓ | — |
| R12 | e + F₂ → F₂(A1) + e | `e_F2_excitation_A1` | [1][3] | ✓ | — |
| R13 | e + F₂(a3) → 2F + e | Tₑ-fit, K₀ = 2.08×10⁻¹⁴ m³/s | [12][3] | ✗ | Chen prints 2.08×10⁻⁸ "m³/s" — ~10⁴× the gas-kinetic ceiling; read as cm³/s |
| R14 | e + F₂(A1) → 2F + e | Tₑ-fit, K₀ = 5.61×10⁻¹⁵ m³/s | [12][3] | ✗ | Same cm³/s reading of Chen's print |
| R15 | e + F₂ → F₂⁺ + 2e | `e_F2_ionization` | [1][3] | ✓ | — |
| R16 | e + Xe → Xe\* + e | `e_Xe_excitation_1` | [1][3] | ✓ | — |
| R17 | e + Xe → Xe⁺ + 2e | `e_Xe_ionization` | [1][3] | ✓ | Xe⁺ has no loss channel in Chen **or** Ma (inherited orphan; documented, not patched) |
| R18 | e + Ne₂\* → Ne₂⁺ + 2e | Tₑ-fit (R&K form) | [5][3] | ✓ | — |
| R19 | e + Ne₂⁺ → Ne + Ne\* | 3.7×10⁻¹⁴·Tₑ^**−0.5** | [5][3] | ✗ | Chen prints Tₑ^+0.5 — a sign flip vs the measured (Philbrick-class) falloff; his form makes recombination rise with Tₑ (unphysical) |
| R20 | 2Ne\* → Ne⁺ + Ne + e | 5×10⁻¹⁶ | [4][3] | ✓ | Chen addition (not in Ma); value sane, kept |
| R21 | Ne\* + Ar → Ar⁺ + Ne + e | 6.5×10⁻¹⁷ | [15][3][4] | ✓ | Penning; matches across chain |
| R21b | Ne\* + F₂ → Ne + 2F | 4.2×10⁻¹⁶ | [11] | **absent** | Measured total quench missing from Chen *and* Ma — a ~15% Ne\* sink at 0.1% F₂; NeF\* is unbound so booked dissociative |
| R22 | F₂⁺ + F⁻ → 3F | 4×10⁻¹⁴ | [4][3] | ✓ | — |
| R23 | Ar\* + Ar + Ne → Ar₂\* + Ne | 1.1×10⁻⁴⁴ m⁶/s | [4][6] | ✗ | Chen prints 2.5×10⁻⁴³ — provably his own R27 (Ar⁺ dimerization) pasted onto the excimer row; Ma's value = the SRI/Kannari literature rate |
| R24 | Ne\* + 2Ne → Ne₂\* + Ne | 4.0×10⁻⁴⁶ m⁶/s | [5][3][4] | ✗ | Chen prints 4×10⁻⁴⁰ (cm⁶→m⁶ slip, ×10⁶: gives Ne\* a 0.4 ps lifetime); Ma prints 8×10⁻⁴⁶ (R&K third-body column transposed, ×2). Primary triple-confirmed: R&K 4.0, Phelps-1959 5, Lee-1989 4.1 (×10⁻³⁴ cm⁶/s) |
| R25 | Ne⁺ + Ar → Ar⁺ + Ne | 1×10⁻¹⁷ | [8][3] | ✓ | Primary-confirmed (Vitols & Oskam in holdings) |
| R26 | Ar⁺ + Ar + Ne → **Ar₂⁺** + Ne | 1.6×10⁻⁴⁴ m⁶/s | [4][3] | ✗ (products) | Both papers print "→ Ne⁺ + 2Ar" — endothermic +5.8 eV and orphans their own Ar₂⁺; species-balance + Chen Fig. 2 prove the dimerization reading. Rate as printed |
| R27 | Ar⁺ + 2Ar → Ar₂⁺ + Ar | 2.5×10⁻⁴³ m⁶/s | [14][3] | ✓ | Kannari/Liu-Conway-confirmed |
| R28 | Ne⁺ + 2Ne → Ne₂⁺ + Ne | 4.4×10⁻⁴⁴ m⁶/s | [5][3] | ✓ | Clean end-to-end (R&K→Ma→Chen) |
| R29 | Ar\* + F₂ → ArF\* + F | 7.5×10⁻¹⁶ | [6][7] | ✗ | Chen prints 4.0×10⁻¹⁶ (Akashi). Lee/Kannari's own compilation books Velazco's full measured quench as ArF\*-forming; adopted (A/B-tested — the split reading collapsed output ×80) |
| R29b | Ar\* + F₂ → Ar + 2F | 3.1×10⁻¹⁶ **DISABLED** | [4][6] | — | Ma's dissociative branch; superseded by [6] (would re-create a >total double count with R29 = 7.5) |
| R29c | Ar\*\* + F₂ → ArF\* + F | 7.5×10⁻¹⁶ | [6] | **absent** | Missing from Chen and Ma: the 4p harpoon (Lee Table I ref 53) — routes the Ar\*\* pool into formation |
| R30 | Ar⁺ + F⁻ → ArF\* | 2.16×10⁻¹² | [10][9][3] | ✗ | Chen prints 9.1×10⁻¹³; measured KrF analogue (Mezyk-Cooper) + Flannery pressure-forms (≈2.5×10⁻¹² at 3.4 atm via Lee Table I) support ~2×10⁻¹². Labeled ANALOGY/THEORY-VALIDATED |
| R31 | Ar₂⁺ + F⁻ → ArF\* + Ar | 1×10⁻¹² | [9][3] | ✗ | Chen prints 9.26×10⁻¹³ 2-body (Ma's ×N form dimensionally ambiguous); Flannery-form value adopted, labeled THEORY |
| R32 | ArF\* + e → Ar + F + e | 2×10⁻¹³ | [14][3][4] | ✓ | — |
| R33 | ArF\* + F₂ → Ar + 3F | 1.9×10⁻¹⁵ | [14][3][4] | ✓ | — |
| R34 | ArF\* + Ar → 2Ar + F | 9×10⁻¹⁸ | [14][3][4] | ✓ | — |
| R35 | ArF\* + Ne → Ar + Ne + F | 1.6×10⁻¹⁸ | [4][3] | ✓ | — |
| R36 | ArF\* + Ar + Ne → Ar₂F\* + Ne | 3.5×10⁻⁴³ m⁶/s | [4][3] | ✗ | Chen prints 3.5×10⁻⁴⁵ — exponent typo vs his cited source; Ma's value = the literature trimer rate (a real ~25%-of-radiative ArF\* sink). Ar₂F\* loss channel absent in both papers (inherited) |
| R37 | ArF\* → Ar + F + hv193 | τ from 2.381×10⁸ s⁻¹ | [13][3] | ✓ | — |
| R38 | ArF\* + hv193 → stimulated | σ_STE = 4×10⁻²⁰ m² | [13][3] | ✓ | Note: Lee [6] uses 2.75×10⁻²⁰ m² — flagged, not adopted |
| R39 | Ar\*\* → **Ar\*** + hv | 10⁹ s⁻¹ | [13][3] | ✗ (products) | Chen prints "→ Ar + hv95" — E1-parity-forbidden (even→even); read as the 4p→4s cascade feeding Ar\*. Rate as printed. (Ma's collisional-only alternative tested and rejected — pinned negative) |
| R40 | Ar\*\* + hv193 → Ar⁺ + e | photo, σ = 10⁻²² m²; e birth 3.76 eV | [13][3][17] | ✓ | Exact frozen-fluence hazard both attenuates photons and consumes Ar\*\*, creates Ar⁺ + e, and books `hv` minus the NIST-anchored excited-state binding energy |
| R41 | F⁻ + hv193 → F + e | photo, σ = 10⁻²¹ m²; e birth 3.02 eV | [13][3][18] | ✓ | Same exact species/photon feedback as R40; birth energy is `hv` minus EA(F). Lee [6] prints 10⁻²² m² for F⁻ — flagged, not adopted |

**Legend** — ✓ matches Chen's print; ~ same source, resolution/threshold detail differs; ✗ deviates (note explains); **absent** = channel Chen's table lacks.

## Cited sources

1. SIGLO database, www.lxcat.net (retrieved 2026-07-30) — Chen's electron-impact cross-section source (his ref 22).
2. BSR database, www.lxcat.net (Zatsarinny & Bartschat B-spline R-matrix; state-resolved Ar ground→4s/4p).
3. J. Chen, L. Zhang, Z. Zhao, Y. Fu, *Phys. Plasmas* **32**, 113502 (2025) — the replication target (Table III).
4. X. Ma, L. Bai, Y. Zhu, X. Jiang, Y. Wu, *Plasma Sources Sci. Technol.* **33**, 075012 (2024) — Chen's ref 13 (heavy-particle set).
5. S. Rauf, M. J. Kushner, *J. Appl. Phys.* **85**, 3460 (1999) — Ne/Xe primary behind Ma's [55]-cited rows.
6. C. H. Lee, T. Matsui, F. Kannari, M. Obara, *IEEE Trans. Electron Devices* **36**, 2054 (1989), Table I — ArF-specific rate compilation.
7. J. E. Velazco, J. H. Kolts, D. W. Setser, *J. Chem. Phys.* **65**, 3468 (1976) — measured Ar(3P2)+F₂ quench.
8. A. P. Vitols, H. J. Oskam — measured Ne⁺ three-body conversion (in holdings).
9. M. R. Flannery, ion–ion recombination theory: *Applied Atomic Collision Physics* (1982) (in holdings) and M. R. Flannery, T. P. Yang, *Appl. Phys. Lett.* **32** (1978).
10. S. M. Mezyk, R. Cooper, D. Sherwell — measured rare-gas–halide ion–ion recombination (KrF analogue).
11. J. H. Brom, J. H. Kolts, D. W. Setser, *Chem. Phys. Lett.* **55**, 44 (1978) — measured Ne(3P2)+F₂ total quench.
12. J. Tennyson et al., Quantemol QDB 2021 release — Chen's ref 38 (F₂ excited-state dissociation fits).
13. H. Akashi, Y. Sakai, H. Tagashira — self-sustained ArF discharge model (Chen's ref 30; radiative/optical constants).
14. A. Demyanov, L. Feenstra, P. Peters, A. Napartovich, W. Witteman et al. — ArF quench set (Chen's cited source on R27/R32–R34).
15. F. Kannari, M. Obara, T. Fujioka — advanced kinetic model (Penning rates).
16. TRINITI database, www.lxcat.net — the only measured e+Ar\*(4s) ionization set in local holdings.
17. NIST Atomic Spectra Database — Ar ionization energy and 3p⁵4p manifold energies used for the
    R40 point photoelectron birth energy.
18. C. Blondel, P. Cacciani, C. Delsart, R. Trainham, *Phys. Rev. A* **40**, 3698 (1989) —
    measured F electron affinity used for the R41 point photoelectron birth energy.

## Open items carried in the table

- Xe⁺ and Ar₂F\* orphans remain inherited from [3][4]. R40/R41 product feedback and point
  birth-energy booking are live; resolved photoelectron birth-energy distributions are not.
- σ_STE and σ(F⁻) discrepancies vs Lee [6] — flagged for a sensitivity pass, not adopted.
- The e+Ar\*(4s)→Ar⁺ record [16] is single-source; an independent measurement/BSR computation remains the top data want.
