# Prototype Runner

This directory contains a minimal Python prototype for the proposed MPToolkit
integration-test schema.

The current prototype is intentionally narrow. It is designed to validate the
core model rather than provide a complete test framework.

Implemented features:

- YAML suite loading
- raw action commands written in a suite file using the actual MPToolkit CLI
- reusable probe recipes imported from YAML
- shorthand lowering for `run`, `steps`, and probe assertions
- action recipes
- probe recipes
- float and complex scalar extraction
- JSON field extraction
- fixture dependencies via `from`
- fixture certification
- test execution with reusable fixtures
- implicit fixture-scoped `cwd`
- implicit current-state tracking for `state`, `psi1`, and `psi2`
- bare tool names resolved against `--bin-dir`
- fixture dependency validation, including cycle detection
- action output enforcement: listed outputs must appear, and undeclared created
  files are treated as errors
- regex, contains, and exists comparisons in addition to `approx` and `equals`

Not yet implemented:

- persistent content-addressed fixture cache
- parallel execution
- text-pattern extraction beyond scalar helpers
- recipe inheritance
- CLI result formats beyond the current prototype reporter

Usage:

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-expectation.yaml \
  --bin-dir /home/ian/build/main-optimized
```

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-dmrg.yaml \
  --bin-dir /home/ian/build/main-optimized
```

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized
```

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-transforms.yaml \
  --bin-dir /home/ian/build/main-optimized
```

Debugging modes:

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized \
  --explain
```

Shows the execution plan for normal developers:

- fixture dependencies
- cache and scratch directories
- resolved commands
- expected output files or glob patterns

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized \
  --trace
```

Shows execution details for debugging test failures:

- stdout and stderr for each command
- extracted probe values
- assertion pass details

```bash
python3 test-proposal/prototype/run_suite.py \
  test-proposal/examples/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized \
  --dump-ir
```

Shows the fully lowered internal representation. This is intended for
implementers of the test framework, not ordinary test authors.

The first prototype target exercises:

- `spinchain-u1`
- `mp-construct`
- `mp-norm`
- `mp-expectation`

and validates real, complex, `--real`, and `--imag` expectation values on a
small finite Neel product state.

The next prototype target exercises:

- `spinchain-u1`
- `mp-random`
- `mp-dmrg`
- `mp-attr`

and validates a derived fixture chain from a fixed-seed random state to a
ground state, plus the legacy `mp-dmrg -s` sweep count path.

The next tool slice exercises:

- `mp-tebd`
- `mp-itebd`
- `mp-imoments-cross --json`

and validates finite and infinite derived-state workflows, plus a JSON-emitting
probe that reuses an iTEBD output fixture.

The next transform slice exercises:

- `mp-overlap`
- `mp-scale`
- `mp-normalize`
- `mp-conj`

and validates in-place scaling and normalization, plus an out-of-place
conjugation workflow checked via explicit output bindings.

The author-facing syntax is now shorter than the original canonical form:

- fixture and test actions can be written as the real commands users run, for
  example `mp-construct -l lat -o psi --finite ...`
- suites can import a small shared probe library such as
  [../recipes/mptoolkit-probes.yaml](/home/ian/sync/git/main/test-proposal/recipes/mptoolkit-probes.yaml)
- probe recipes can use `probe: ...` and action recipes can use `action: ...`
  instead of spelling out `kind` plus `command`
- recipe-local parameter defaults can use `defaults:` instead of `params:`
  when that reads more naturally
- simple extractors can use `extract: last_float` or
  `extract: {json_path: "{json_path}"}`
- assertions can use probe-name shorthand such as `norm:`, `attr_float:`,
  `expectation:`, or `imoments_cross_json:`
- the working directory is inferred from fixture scope in normal cases, so test
  authors usually do not write `cwd:` at all
- the default state follows the current fixture output, so test authors usually
  do not write `state: psi` either
- single-fixture `use:` can be written as a scalar instead of a one-element
  list
- `use:` can also bind named fixture outputs, for example
  `psi1: left.state` and `psi2: right.state`

The intended split is now:

- Python runner: generic mechanics common to nearly all commands
- suite files: the actual action commands for the workflow being tested
- imported probe recipes: reusable checks such as `norm`, `expectation`, and
  `attr`

Fixture scope rules:

- fixture build and certify actions run in the fixture directory
- tests that `use` a single fixture run steps and probes relative to that
  fixture by default
- mapped `use:` bindings are resolved against the test-local materialized copy
  of each fixture output, not the shared fixture cache
- multi-fixture tests can override scope with `fixture: some_fixture`
- raw `cwd:` still exists as an escape hatch, but is not part of ordinary
  authoring

Default state rules:

- probes and actions that use `state`, `psi1`, or `psi2` default to the
  current fixture's primary state output
- derived fixtures automatically switch that default to the most recent state
  produced by their build step
- tests that `use` a single fixture inherit that default state automatically

Command rules:

- action commands may be written as a shell-like string or as an argv list
- for string commands, the runner tokenizes the template first and then
  substitutes placeholders per token, so placeholder values are not re-split by
  whitespace
- if the command name has no slash and exists in `--bin-dir`, the runner uses
  that binary automatically
- tests should rely on fresh fixture/test directories and unique output names,
  not on `--force`; inconsistent `--force` handling across tools should be
  treated as a separate bug
