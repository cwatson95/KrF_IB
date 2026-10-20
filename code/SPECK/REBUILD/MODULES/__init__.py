"""The engine modules. Canonical import path: ``REBUILD.MODULES.<name>``.

No aliases on the package root, no import-time side effects, no module-level
mutable state — all three are gated by ``tests/test_architecture.py``.

Layering (low → high; an import that inverts it fails the architecture gate):

    units
    provenance, params
    species, chen2025, cross_sections
    network, loki_tier, kinetic_tier, beam
    ledger, laser, ase, deposition3d, chemistry
    engine, checkpoint, scorecard
"""
