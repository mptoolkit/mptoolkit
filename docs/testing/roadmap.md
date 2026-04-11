# Roadmap

## Phase 1: Python Runner

Build a small Python engine around `subprocess` with:

- per-test scratch directories
- stable environment defaults such as `LC_ALL=C`, `OMP_NUM_THREADS=1`,
  `MP_BINPATH=<tmp>`
- action recipes
- probe recipes
- typed comparisons

At this stage the test catalog can still be written in Python, but the runtime
objects should already match the declarative model.

## Phase 2: Fixture Cache

Add:

- content-addressed fixture cache
- fixture certification
- read-only cached outputs
- copy-on-use semantics for tests
- locking so parallel workers do not build the same fixture twice

## Phase 3: Declarative Loader

Add YAML loading for:

- recipes
- fixtures
- tests

The YAML surface syntax should lower into the same canonical objects used by
the Python runner.

## Phase 4: Standard Library

Add a small built-in library of common MPToolkit probes and helpers:

- `norm`
- `overlap`
- `attr`
- `history_field`
- `json_field`
- `stdout_last_float`

The aim is that most tests use named probes rather than custom parsing.

## Phase 5: Coverage Growth

Start with fast deterministic tests:

- TEBD imaginary-time norm scaling under `sum_unit(I(0))`
- TEBD with `--normalize`
- real-time norm preservation
- machine-readable JSON output checks
- cross-tool consistency such as `overlap(psi, psi)` versus `norm(psi)`

Then add:

- small exact-reference tests for tiny systems
- regression fixtures for previously fixed bugs
- documentation-derived scenarios from the wiki HOWTO pages

## Phase 6: CI

For pull requests:

- build once
- run the fast integration suite

For nightly or manual jobs:

- run slower exact-reference tests
- run a broader scenario catalog

## Open Questions

- Which tool outputs should gain JSON or scalar modes first?
- Which fixture recipes should become part of the first standard library?
- Whether YAML is the final format, or only the first declarative surface.
