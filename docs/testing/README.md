# Testing Design Notes

This directory contains design and rollout notes for the MPToolkit integration
test system.

Executable test assets live under [tests](/home/ian/sync/git/main/tests). These
documents describe the schema, matching rules, and coverage plan that guide the
runner and suite catalog.

Key documents:

- [schema.md](/home/ian/sync/git/main/docs/testing/schema.md): core suite,
  fixture, recipe, and assertion model
- [pattern-language.md](/home/ian/sync/git/main/docs/testing/pattern-language.md):
  line-oriented typed text matching
- [roadmap.md](/home/ian/sync/git/main/docs/testing/roadmap.md): implementation
  direction and future cleanup items
- [checklist.md](/home/ian/sync/git/main/docs/testing/checklist.md): concrete
  rollout checklist for tools, models, and workflow integration
- [examples/identity-tebd.yaml](/home/ian/sync/git/main/docs/testing/examples/identity-tebd.yaml):
  illustrative design-only example

The executable catalog is currently in:

- [tests/README.md](/home/ian/sync/git/main/tests/README.md)
- [tests/recipes/mptoolkit-probes.yaml](/home/ian/sync/git/main/tests/recipes/mptoolkit-probes.yaml)
- [tests/suites](/home/ian/sync/git/main/tests/suites)

The current direction is:

- keep workflow actions written as the real CLI commands users run
- keep reusable checks as shared probe recipes
- keep fixture isolation, command resolution, and output enforcement in Python
- integrate execution into the normal build workflow through `scripts/mptk-test`
  and `make` targets run from the build directory
