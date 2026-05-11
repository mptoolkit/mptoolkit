# Codex Instructions

These rules are for Codex agents working in this repository. Keep this file
portable; machine-specific paths, compiler selections, package-management rules,
and personal build preferences belong in `.codex/instructions.local.md`, which
is intentionally ignored.

## Instruction Precedence

Before substantive editing, build, or test work:

1. Read `AGENTS.md`.
2. If `.codex/instructions.local.md` exists, read it for this checkout's local
   paths and tool policy.

If there is a conflict, follow:

- user, developer, and system instructions
- `AGENTS.md` for repository build, test, and migration policy
- `.codex/instructions.local.md` for machine-specific details compatible with
  the repository policy above
- this file

## Build And Test Discipline

- Use out-of-tree builds only.
- Use explicit user instructions or the local override file to choose the build
  directory. If neither names one, create a task-specific build directory
  outside the source tree.
- Do not run `configure`, `make`, or integration tests from the source root.
- Run focused integration tests after behavior changes unless explicitly told not
  to.

## Change Scope

- Keep changes focused and minimal.
- Do not touch unrelated code or local scratch files.
- Keep C++17, C++20, CMake, EXPOKIT, Boost, and uni20 migration work separated
  unless the user explicitly asks to combine them.

## GitHub Notes

- Verify current branch, PR head, and PR base before describing PR status.
- Keep PR summaries free of local machine paths.
- After rebases, update stale branch-stacking language in PR descriptions.
