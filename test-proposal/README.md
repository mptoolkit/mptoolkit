# MPToolkit Test Proposal

This directory is a design workspace for a modern integration and end-to-end
test system for MPToolkit.

The near-term implementation target is a Python runner built on `subprocess`.
The long-term target is a declarative test format where workflow actions are
written as natural CLI commands and reusable probe recipes provide the common
checks.

The main design goals are:

- Make it easy to add a new test for an existing tool or option.
- Reuse certified wavefunction fixtures without introducing test ordering.
- Support noisy CLI output by parsing typed values rather than exact text.
- Keep the stable contract at the level of MPToolkit workflows, not low-level
  backend implementation details.
- Make the eventual declarative format clean enough for non-experts to write.

Documents in this directory:

- [schema.md](/home/ian/sync/git/main/test-proposal/schema.md): Core object
  model and execution semantics.
- [pattern-language.md](/home/ian/sync/git/main/test-proposal/pattern-language.md):
  Line-oriented matching and typed placeholder syntax for text output.
- [roadmap.md](/home/ian/sync/git/main/test-proposal/roadmap.md): Suggested
  implementation path from Python-first to a stable declarative layer.
- [checklist.md](/home/ian/sync/git/main/test-proposal/checklist.md): Concrete
  rollout checklist for tool coverage, model coverage, and CI growth.
- [prototype/README.md](/home/ian/sync/git/main/test-proposal/prototype/README.md):
  Minimal executable prototype for the first vertical slice.
- [examples/spinchain-expectation.yaml](/home/ian/sync/git/main/test-proposal/examples/spinchain-expectation.yaml):
  First executable suite targeting `spinchain-u1`, `mp-construct`,
  `mp-norm`, and `mp-expectation`.
- [examples/spinchain-dmrg.yaml](/home/ian/sync/git/main/test-proposal/examples/spinchain-dmrg.yaml):
  Executable suite covering `mp-random`, `mp-dmrg`, and `mp-attr`, including a
  derived ground-state fixture and a regression test for the legacy `-s` path.
- [examples/spinchain-tebd.yaml](/home/ian/sync/git/main/test-proposal/examples/spinchain-tebd.yaml):
  Executable suite covering `mp-tebd`, `mp-itebd`, and `mp-imoments-cross
  --json`, including reuse of an iTEBD-derived fixture in a JSON probe.
- [examples/spinchain-transforms.yaml](/home/ian/sync/git/main/test-proposal/examples/spinchain-transforms.yaml):
  Executable suite covering `mp-overlap`, `mp-scale`, `mp-normalize`, and
  `mp-conj`, including explicit output bindings for cross-state checks.
- [examples/spinlessfermion-dmrg.yaml](/home/ian/sync/git/main/test-proposal/examples/spinlessfermion-dmrg.yaml):
  Executable suite covering finite `spinlessfermion-u1`, `mp-construct`,
  `mp-random`, `mp-dmrg`, and sign-sensitive off-diagonal fermion correlators.
- [examples/hubbard-u1u1-dmrg.yaml](/home/ian/sync/git/main/test-proposal/examples/hubbard-u1u1-dmrg.yaml):
  Executable suite covering finite `hubbard-u1u1`, `mp-construct`, `mp-dmrg`,
  and sign-sensitive up-spin hopping correlators in a tiny Hubbard sector.
- [examples/bosehubbard-u1-dmrg.yaml](/home/ian/sync/git/main/test-proposal/examples/bosehubbard-u1-dmrg.yaml):
  Executable suite covering finite `bosehubbard-u1`, `mp-construct`,
  `mp-dmrg`, and a tiny one-boson hopping sector.
- [examples/identity-tebd.yaml](/home/ian/sync/git/main/test-proposal/examples/identity-tebd.yaml):
  A worked example showing fixtures, derived fixtures, probes, and test reuse.

The prototype now has two layers:

- a verbose canonical internal form used by the runner
- a lighter author-facing syntax with raw action commands, imported probe
  recipes, and probe assertion shorthand

The author-facing form now also infers working directory from fixture scope, so
normal tests do not need to mention `cwd`, and follows the current fixture
state so normal tests do not need to mention `state: psi`.

Its current direction is:

- keep only cross-cutting mechanics in Python, such as fixture isolation,
  `--bin-dir` command resolution, implicit fixture scope, and default current
  state
- write fixture and test actions as the real commands a user would run
- keep reusable checks like `norm`, `expectation`, and `attr` as imported probe
  recipes

It also now has two user-facing debug modes:

- `--explain` for normal developers, showing dependencies, files, and resolved
  commands
- `--trace` for execution debugging, showing command output and extracted values

The separate `--dump-ir` mode is intended for implementers of the test runner.

The current direction is to preserve the explicit internal model while making
the author-facing form concise enough for normal contributors without hiding
workflow commands behind a second layer of action-specific sugar.
