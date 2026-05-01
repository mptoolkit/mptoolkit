# Testing Design Notes

This directory contains design and rollout notes for the MPToolkit integration
test system.

Executable test assets live under [tests](../../tests). These
documents describe the schema, matching rules, and coverage plan that guide the
runner and suite catalog.

For day-to-day test execution, start with
[tests/README.md](../../tests/README.md). The short version is:

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test
```

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test spinchain-tebd
```

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test \
  spinchain-tebd \
  --test itebd_and_itdvp_real_time_ising_match_phase_and_fidelity_per_unit_cell
```

Key documents:

- [schema.md](schema.md): core suite,
  fixture, recipe, and assertion model
- [pattern-language.md](pattern-language.md):
  line-oriented typed text matching
- [roadmap.md](roadmap.md): implementation
  direction and future cleanup items
- [checklist.md](checklist.md): concrete
  rollout checklist for tools, models, and workflow integration
- [examples/identity-tebd.yaml](examples/identity-tebd.yaml):
  illustrative design-only example

The executable catalog is currently in:

- [tests/README.md](../../tests/README.md)
- [tests/recipes/mptoolkit-probes.yaml](../../tests/recipes/mptoolkit-probes.yaml)
- [tests/suites](../../tests/suites)

The current direction is:

- keep workflow actions written as the real CLI commands users run
- keep reusable checks as shared probe recipes
- keep fixture isolation, command resolution, and output enforcement in Python
- integrate execution into the normal build workflow through `scripts/mptk-test`
  and `make` targets run from the build directory
