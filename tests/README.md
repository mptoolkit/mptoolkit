# Integration Tests

This directory contains the declarative MPToolkit integration test system.

Layout:

- `run_suite.py`: the test runner
- `recipes/`: reusable probe recipes
- `suites/`: executable YAML suites
- `data/`: checked-in test data when needed

The runner executes real MPToolkit commands, builds reusable fixtures in
isolated directories, and checks results with typed probes such as norms,
overlaps, expectation values, attributes, JSON fields, and line-pattern
matches.

Typical usage from the source tree:

```bash
python3 tests/run_suite.py \
  tests/suites/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized
```

```bash
python3 tests/run_suite.py \
  tests/suites/spinchain-tebd.yaml \
  --bin-dir /home/ian/build/main-optimized \
  --trace
```

Typical usage from a build directory:

```bash
/home/ian/sync/git/main/scripts/mptk-test
```

```bash
/home/ian/sync/git/main/scripts/mptk-test spinchain-tebd
```

The wrapper defaults `--bin-dir` to the current working directory, so it is
intended to be run from the build tree.

Supported debug modes:

- `--explain`: show fixture dependencies, scratch directories, resolved
  commands, and expected outputs
- `--trace`: show command output, extracted values, and assertion details
- `--dump-ir`: dump the normalized internal suite representation

Representative suites:

- `suites/spinchain-expectation.yaml`: finite construction and scalar probes
- `suites/spinchain-dmrg.yaml`: finite ground-state workflow
- `suites/spinchain-aklt-spt.yaml`: AKLT SPT commutator checks via
  `mp-aux-algebra`
- `suites/spinchain-tebd.yaml`: finite and infinite evolution workflows
- `suites/spinchain-ibc.yaml`: first infinite-boundary-condition workflow checks
- `suites/spinchain-ibc-dynamics.yaml`: IBC TDVP and correlation prefix checks
- `suites/spinchain-ea.yaml`: excitation-ansatz helper and wavepacket checks
- `suites/spinchain-transforms.yaml`: scale/normalize/conjugation behavior
- `suites/spinchain-infinite-analysis.yaml`: moments and spectrum-style output

The runner still has a few planned extensions, notably persistent fixture
caching and parallel execution, but it is now organized as the project’s real
integration-test tree rather than a separate proposal area.
