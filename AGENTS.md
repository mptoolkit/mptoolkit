# MPToolkit Agent Guidance

This file contains portable repository-level guidance for agents working on
MPToolkit. Machine-local paths, compiler selections, package-management rules,
build-tree names, parallelism defaults, and personal workflow preferences belong
in ignored local override files, not here.

## Status And Authority

- These notes are guidance for development and review; they do not replace the
  source code, tests, or PR review comments.
- If this file disagrees with the code or an explicit user instruction, inspect
  the code and follow the user instruction.
- Prefer concrete local or CI verification over broad claims about the toolkit.
- Keep this file portable. Do not add host paths, user-specific symlink names,
  local package-manager rules, CPU-count assumptions, or scratch-workflow notes
  unless they are true repository requirements.

## Project Context

- MPToolkit is currently transitioning through C++17/C++20 cleanup, Boost
  reduction, and eventual CMake migration.
- Autoconf is still the active build system on `main` unless the current task is
  explicitly on the CMake branch.
- The legacy `LinearAlgebra` tree is expected to be replaced by uni20 later.
  Avoid large cosmetic refactors there unless they are needed for the task.
- The parser stack still uses Boost.Spirit and is expected to be replaced later.
  Do not spend effort modernizing parser internals unless directly requested.
- Boost.Program_options remains an accepted dependency for now.

## General Build Policy

- Use out-of-tree builds. Do not configure, compile, or run integration tests
  from the source root.
- Choose the concrete build directory from explicit user instructions, ignored
  local override files, or a task-specific out-of-tree directory. Do not encode
  machine-specific build paths or symlink names as repository policy.
- Configure from the build directory using the source tree's `configure` script.
- Build from the build directory, for example:
  `make -C <build-dir> <target>`
- Do not run plain `make` from the source root.

## Integration Tests

- The integration suites are currently the main practical regression tests.
- Run integration tests against the build directory, for example:
  `scripts/mptk-test --bin-dir <build-dir> [suite ...]`
- For focused checks, pass suite names explicitly.
- Use temporary work roots unless the task requires preserving outputs.
- Old low-level tests, especially in `linearalgebra/test`, may be useful smoke
  tests but should not be treated as comprehensive coverage.
- When touching dense matrix exponential or time evolution paths, run focused
  TEBD/TDVP-facing suites such as:
  - `spinchain-tebd`
  - `spinchain-matrix-tools`
  - `spinchain-ibc-dynamics`

## Autoconf Outputs

- Let Make regenerate Autoconf outputs when `configure.ac`, `aclocal.m4`, or
  `stamp-h.in` is actually out of date.
- If Autoconf runs from the source tree, clean source-local detritus afterward:
  `autom4te.cache/` and `configure~`.
- Do not hand-edit generated `configure` or `config.h.in`; edit the source
  Autoconf inputs and regenerate.

## Source Tree Hygiene

- If `Makefile`, `config.h`, `*.o`, `*.d`, or toolkit executables appear in the
  source root, treat that as an accidental in-source build and clean it up
  before continuing.
- Do not remove unrelated untracked files unless explicitly asked.
- Before finishing build or test work, check:
  `git status --short --branch`
  and
  `find . -maxdepth 1 \( -name 'Makefile' -o -name 'config.h' -o -name '*.o' -o -name '*.d' -o -name 'autom4te.cache' -o -name 'configure~' \) -print`

## Answering And Review Rules

- State what was actually inspected or tested.
- Separate implemented behavior from roadmap material.
- Do not claim broad numerical or physics correctness from a compile-only check.
- For branch or PR status, verify the actual GitHub base/head before describing
  stacked branch relationships.
- Do not include local machine paths in PR descriptions or comments unless the
  path is essential to explain a local-only diagnostic.
- Resolve review threads only after the code or explanation directly addresses
  the comment.

## Migration Roadmap Cautions

- Do not mix unrelated migration layers in one change unless explicitly asked.
- C++17 cleanup, C++20 cleanup, EXPOKIT removal, CMake migration, and uni20
  replacement work should normally remain separate PRs.
- Do not remove Boost.Spirit or Boost.Program_options opportunistically.
- Do not replace `LinearAlgebra` piecemeal with uni20 unless the task is
  explicitly about that replacement.
