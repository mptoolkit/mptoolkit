# `.codex/` Layout

This directory stores Codex-specific guidance for this checkout.

`AGENTS.md` and `.codex/instructions.md` are intended to be shared repository
guidance. `.codex/instructions.local.md` is local-only.

## Files

- `.codex/instructions.md`: Codex entry point for this repository.
- `.codex/instructions.local.md`: machine-specific build paths, local package
  policy, and preferred commands. This file is ignored and optional.

## What Belongs Where

- Put stable, portable repository rules in `AGENTS.md`.
- Put Codex-specific instruction precedence and routing in `.codex/instructions.md`.
- Put machine-local paths and build tree names in `.codex/instructions.local.md`.
- Do not put local filesystem paths in PR descriptions or committed docs unless
  they are essential to explain a local diagnostic.
- If a local-only note becomes useful to everyone, move the portable part into
  `AGENTS.md` or `.codex/instructions.md` and leave host paths, symlink names,
  package-manager policy, and personal defaults behind.

## Portable Versus Local

Potentially portable guidance:

- out-of-tree builds
- source-tree hygiene
- Autoconf generated-file policy
- integration-test expectations
- migration-roadmap cautions
- PR summary hygiene

Local-only guidance:

- absolute build paths
- local build-tree names and symlink targets
- local compiler/package policy
- local parallelism defaults
- personal scratch-file expectations
