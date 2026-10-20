# data/tests — fixture data for the test suite

Files here are **test fixtures**: the unit/integration tests under
`tests/` load these files by path, and deleting or editing them breaks the
suite. Nothing in this directory is read by the simulation engine or the
examples — experiment inputs live in `data/cross_sections/` (σ datasets
referenced from `data/reactions/*.toml`) and are free to be regenerated,
re-exported, or replaced without touching the tests.

The `example_` prefix marks a file as test-only:

- `example_e_Kr_ionization.{csv,json}`, `example_e_Krs_1s.{csv,json}`,
  `example_e_F2_attachment.{csv,json}` — frozen snapshots of the
  Cross-Section Browser exports (total σ(E); the F₂ set is a cm²-column
  export, deliberately kept to exercise the units-metadata conversion).
  Their JSON `process_id`/`data_file` fields carry the `example_` name so
  key-based resolution is tested against this directory alone.

If a test needs new fixture data, copy it here with the `example_` prefix
(updating `data_file`/`process_id` inside the JSON to match) rather than
pointing the test at `data/cross_sections/` — tests must never depend on
data that experiments are allowed to change. Small synthetic fixtures used
by a single test module may also live in `tests/fixtures/`.
