# Declarative Test Schema

## Scope

The proposed system is aimed at MPToolkit integration and end-to-end tests.
Tests should describe workflows such as state construction, time evolution,
norm checks, overlap checks, metadata checks, and machine-readable CLI output.

This proposal intentionally treats the declarative syntax as a surface layer.
The underlying model should be explicit enough to implement first in Python and
later load from YAML or another declarative format.

## Design Principles

- Tests depend on fixtures, never on other tests.
- Shared fixtures are immutable and may be cached.
- Each test runs in its own scratch directory.
- Text scraping logic lives in reusable recipes, not in individual tests.
- Common probes such as `norm`, `overlap`, `attr`, and `json_field` should be
  part of the standard library.
- Regex is an escape hatch, not the default authoring interface.

## Core Entities

### Suite

A suite is a collection of recipes, fixtures, and tests.

Suggested top-level shape:

```yaml
version: 0
defaults: ...
recipes: ...
fixtures: ...
tests: ...
```

### Recipe

A recipe is a reusable command specification.

There are two primary kinds:

- `action`: runs a command and produces one or more artifacts.
- `probe`: runs a command and extracts a typed value.

A recipe may inherit from another recipe and override fields.

Suggested canonical fields:

```yaml
recipes:
  norm:
    probe: mp-norm {state}
    env: {}
    cwd: "{cwd}"
    capture:
      stream: stdout
      strip_ansi: true
      whitespace: flexible
    extract: last_float

  tebd_step:
    action: mp-tebd -w {input} -H {hamiltonian} -o {output_prefix} -t {time} -n {steps} -s {save_every}
    defaults:
      input: null
      hamiltonian: null
      output_prefix: "{scratch}/{name}"
      time: null
      steps: 1
      save_every: 1
    outputs:
      state:
        discover:
          kind: glob
          pattern: "{output_prefix}*"
          select: newest
```

Notes:

- `action:` and `probe:` are author-facing shorthands for `kind` plus
  `command`.
- `command` is a template. Parameters are resolved before execution.
- Parameter precedence should be:
  `suite defaults < recipe defaults < per-use overrides`.
- `capture` is relevant mainly to probes, but may also be useful for actions
  whose outputs are reported in text.
- `outputs` tells the runner how to name or discover produced artifacts.

### Fixture

A fixture is a named, reusable artifact or set of artifacts.

Fixtures may be:

- built from one or more action steps
- derived from another fixture
- imported from a checked-in file as a temporary escape hatch

Suggested fields:

```yaml
fixtures:
  psi_up_20:
    steps:
      - run:
          recipe: prepare_product_up
          with:
            sites: 20
        save_as:
          state: state
    certify:
      - probe:
          recipe: norm
          with:
            state: psi
            cwd: "{fixture.dir}"
        compare:
          kind: approx
          value: 1.0
          abs: 1e-12

  psi_up_20_beta_0p1:
    from: psi_up_20
    steps:
      - run:
          recipe: tebd_step
          with:
            input: psi
            hamiltonian: "sum_unit(I(0))"
            time: "0-0.1i"
            output_prefix: "{scratch}/psi"
        save_as:
          state: state
    certify:
      - probe:
          recipe: norm
          with:
            state: psi
            cwd: "{fixture.dir}"
        compare:
          kind: approx
          value: 0.1353352832366127
          abs: 1e-12
```

Fixture semantics:

- A fixture is built at most once per cache key.
- `steps` is the general form. `build` may still exist as shorthand for a
  single action step.
- Certification happens immediately after build.
- A failed certification means the fixture is unusable.
- Cached fixture outputs are treated as read-only.

### Test

A test names the fixtures it uses, runs zero or more steps, and asserts on the
resulting artifacts or probe values.

Suggested fields:

```yaml
tests:
  tebd_identity_normalized:
    use:
      - psi_up_20
    steps:
      - run:
          recipe: tebd_step
          with:
            input: "{fixtures.psi_up_20.state}"
            hamiltonian: "sum_unit(I(0))"
            time: "0-0.1i"
            steps: 1
            save_every: 1
            normalize: true
            output_prefix: "{scratch}/psi"
        save_as:
          state: psi1
    expect:
      - probe:
          recipe: norm
          with:
            state: "{state.psi1}"
        compare:
          kind: approx
          value: 1.0
          abs: 1e-12
```

Test semantics:

- Each test runs in a fresh scratch directory.
- Fixture inputs are copied or linked into the scratch area before use.
- Tests may run in parallel once required fixtures are available.
- No test may rely on another test having already run.

## Canonical Expectation Form

The explicit internal form should separate probing from comparison:

```yaml
expect:
  - probe:
      recipe: norm
      with:
        state: "{state.psi1}"
    compare:
      kind: approx
      value: 1.0
      abs: 1e-12
```

Compact sugar can lower to the same representation:

```yaml
expect:
  - norm:
      target: psi1
      approx: 1.0
      abs: 1e-12
```

Supported comparison kinds should include:

- `equals`
- `approx`
- `matches`
- `contains`
- `exists`

## Output Discovery

Some tools produce outputs at a known path. Others produce paths derived from
prefixes, time steps, or internal naming rules. The schema therefore needs a
standard output discovery mechanism.

Proposed discovery modes:

- `path`: the output path is known exactly.
- `glob`: select from a glob pattern.
- `stdout`: parse a path from text output.

Examples:

```yaml
outputs:
  state:
    discover:
      kind: path
      path: "{output}"
```

```yaml
outputs:
  state:
    discover:
      kind: glob
      pattern: "{output_prefix}*"
      select: newest
```

```yaml
outputs:
  state:
    discover:
      kind: stdout
      pattern: "Saved wavefunction to {path}"
```

## Command Representation

Author-facing suites can use simple command strings:

```yaml
command: mp-expectation {state} lat:Sz(0)
```

The runner should split these without invoking a shell. Argv-list form should
also remain valid when authors need exact token control:

```yaml
command: ["mp-expectation", "{state}", "lat:Sz(0)"]
```

## Text Matching

Text output matching should be line-oriented by default and use a small pattern
language with typed placeholders. See
[pattern-language.md](/home/ian/sync/git/main/test-proposal/pattern-language.md).

The intention is that novice test authors write:

```yaml
stdout: "Energy = {float}"
```

not floating-point regexes.

## Mutability and Isolation

The safe default is:

- fixture cache is immutable
- test scratch directories are writable
- action steps operate on local copies of inputs

Some MPToolkit artifacts have companion-file conventions, so a fixture is often
best treated as a directory bundle rather than a single file.

This removes most hazards around tools that overwrite files in place and makes
parallel execution straightforward.

Implementation detail:

- physical copying can later be optimized with hardlinks or reflinks
- logical semantics should still be copy-on-use

## Parallelism

The execution graph is a fixture DAG plus independent tests.

The scheduler should:

1. Build missing fixtures, with locking per cache key.
2. Certify them once.
3. Run tests in parallel after fixture dependencies are satisfied.

This keeps reuse efficient while preserving order independence.

## Escape Hatches

The declarative layer should allow controlled fallbacks for awkward cases:

- custom extractor implemented in Python
- custom comparator implemented in Python
- imported binary fixture for rare or expensive setups

These should exist, but the standard library should cover the common cases well
enough that most tests never need them.
