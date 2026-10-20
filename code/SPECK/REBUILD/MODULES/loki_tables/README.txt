LOKI-tier source data — provenance
==================================

Copied 2026-08-04 from projects/spectre/LoKI/Curated/Output/<run>/
(lookUpTableRateCoeff.txt + lookUpTableSwarm.txt only; the per-E/N job
directories stay with spectre). These are LoKI-MC (LisbOn KInetics Monte
Carlo, N-PRiME/IST Lisbon) boltzmannMC outputs over the actual device
mixtures — the curated PLAN (45)/(48) sweeps, 2e4 electrons per point:

  ChenArF_NeArF2        Ne/Ar/F2 (Chen 2025 ArF point, 96.4/3.5/0.1)
  V2point_NeKrF2        Kr/F2/Ne 5/0.1/94.9 (the v2 KrF point, Ne fill)
  V2point_NeKrF2_lowEN  same mixture, low-E/N extension
  V2point_ArKrF2        Kr/F2/Ar 5/0.1/94.9 (Ar fill)
  V2point_ArKrF2_lowEN  same mixture, low-E/N extension

Files are byte-verbatim copies (verify with cmp against the spectre tree).
Units: rate coefficients m3/s; MeanE(eV) is the swarm mean electron energy
at each reduced field. The engine looks rates up by <eps> = (3/2) Te — the
documented Maxwellian-mapping approximation the KINETIC tier exists to
score (ORESTES_CLASS_REBUILD.txt L3/S4).

Do not edit these files. To regenerate: LoKI-MC lives in
Kinetic_Modeling/LoKI-MC/ (build/run details in its CLAUDE.md); the input
decks are in the spectre project alongside the outputs.

Added 2026-08-26 (0.1.27):

  ChenArF_NeArF2_lowEN  same Chen mixture, low-E/N extension
                        (0.5/1/2/3 Td — the V2point lowEN convention)

Generated this date with spectre's own toolchain (Code/lokimc, arm64;
2e4 electrons/point, boltzmannMC, equalSharing): the ChenArF deck
byte-identical except reducedElecField and the output folder — the deck
is Input/Mixture/ChenArF_NeArF2_lowEN.in in spectre, job dirs under
Code/Output/krf0d_CHEN_ArF_Ne964_Ar035_F2001_lowEN, curated copy at
Curated/Output/ChenArF_NeArF2_lowEN. Before it, the arf_chen circuit
drive rode the flat+taper stand-in for EVERYTHING under 5 Td — 70% of
the Chen anchor's 200-ns window, including the whole breakdown ramp
(the 0.1.25 en_below_table_frac finding). The loader's overlap guard
verifies grid disjointness at load.
