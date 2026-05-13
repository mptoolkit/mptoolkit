# Integration Tests

This directory contains the declarative MPToolkit integration test system.

## Quick Start

Build MPToolkit first, then run the integration wrapper from the build
directory. The wrapper uses the current directory as `--bin-dir`, so it can find
programs such as `mp-construct`, `mp-itebd`, and `spinchain-u1`.

Run all supported integration suites:

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test
```

By default, the wrapper runs top-level `tests/suites/*.yaml` suites. Optional
contrib suites live under `tests/suites/contrib/`; they can be run explicitly
when working on contrib models, but they are not part of the default CI blocker
set.

Run one suite:

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test spinchain-tebd
```

Run one named test inside one suite:

```bash
cd /path/to/mptoolkit-build
/path/to/mptoolkit-source/scripts/mptk-test \
  spinchain-tebd \
  --test itebd_and_itdvp_real_time_ising_match_phase_and_fidelity_per_unit_cell
```

Run from the source tree by passing the build directory explicitly:

```bash
cd /path/to/mptoolkit-source
scripts/mptk-test --bin-dir /path/to/mptoolkit-build spinchain-tebd
```

`--bin-dir` may point at the legacy flat build directory, or at a build output
root containing executables under `bin/`, `tools/`, `models/`,
`models/contrib/`, `bin/tools/`, `bin/models/`, or `bin/models/contrib/`.
For CMake build trees with per-configuration output directories, pass
`--config Release`, `--config Debug`, or another configuration name to prefer
that configuration's outputs.

Run the lower-level Python runner directly when debugging a single suite:

```bash
cd /path/to/mptoolkit-source
python3 tests/run_suite.py \
  tests/suites/spinchain-tebd.yaml \
  --bin-dir /path/to/mptoolkit-build \
  --trace
```

Useful runner flags:

- `--trace`: print command output, extracted probe values, and assertion details
- `--explain`: show fixture dependencies, scratch directories, and resolved
  commands
- `--config CONFIG`: prefer binaries from that CMake configuration when
  resolving tools in per-configuration output directories
- `--work-root DIR`: keep fixture and test outputs under `DIR` instead of a
  temporary directory
- `--dump-ir`: print the normalized suite representation for one suite

Layout:

- `run_suite.py`: the test runner
- `recipes/`: reusable probe recipes
- `suites/`: supported executable YAML suites
- `suites/contrib/`: optional suites for unsupported contrib models
- `data/`: checked-in test data when needed

The runner executes real MPToolkit commands, builds reusable fixtures in
isolated directories, and checks results with typed probes such as norms,
overlaps, expectation values, attributes, JSON fields, and line-pattern
matches.

For manual overnight stress checks of a single suite, use:

```bash
/path/to/mptoolkit-source/scripts/mptk-stress-test \
  --bin-dir /path/to/mptoolkit-build \
  --forever
```

The default stress target is `spinchain-aklt-spt-uc1-rerun.yaml`.

Common stress-test patterns:

```bash
/path/to/mptoolkit-source/scripts/mptk-stress-test \
  --bin-dir /path/to/mptoolkit-build \
  --forever
```

```bash
/path/to/mptoolkit-source/scripts/mptk-stress-test \
  --bin-dir /path/to/mptoolkit-build \
  --suite spinchain-tebd \
  --test itdvp_projected_amplitude_tracks_time_dependent_imaginary_time_identity \
  --iterations 200 \
  --fail-fast
```

```bash
/path/to/mptoolkit-source/scripts/mptk-stress-test \
  --bin-dir /path/to/mptoolkit-build \
  --keep-all-logs \
  --trace \
  --iterations 20
```

Stress-test behavior:

- default log directory is `/tmp/mptk-stress-<suite>-<timestamp>-<pid>/`
- `summary.txt` records pass/fail by iteration
- failed iteration logs are kept automatically
- passing iteration logs are deleted unless `--keep-all-logs` is given
- `--fail-fast` stops on the first failure

See also:

- [docs/testing/stress-testing.md](../docs/testing/stress-testing.md)

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
