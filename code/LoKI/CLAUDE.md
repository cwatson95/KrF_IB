# LoKI-MC (C++, vendored)

Upstream C++ Monte Carlo electron-kinetics solver (LisbOn KInetics, N-PRiME group, GPL). Treat as a vendored external tool: the KrF code consumes its outputs as embedded tables in `phantom11.loki_tables`, not at runtime. Only relevant if regenerating swarm tables.

Build from `LoKI-MC/Code/` via the install scripts (needs g++-9, GSL, Boost, CMake, gnuplot):

```bash
cd LoKI-MC/InstallationCommands
./install_linux.sh all        # installs deps + builds; use 'onlyCode' to skip apt installs
```

Run with `./lokimc <setupFile.in> <numThreads>` from `LoKI-MC/Code/` (setup files in `Code/Input/`).

A native binary exists at `LoKI-MC/Code/lokimc` (built 2026-07-05, conda env `loki` libs) for regenerating tables. Raw LoKI outputs land in `KrF_and_LoKI/LoKI Code/Output/`.
