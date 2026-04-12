# Manual Stress Testing

This note documents the manual stress-test wrapper:

- [scripts/mptk-stress-test](/home/ian/sync/git/main/scripts/mptk-stress-test)

It is intended for hunting intermittent failures that are too rare or too slow
to turn into ordinary CI regressions.

Examples:

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --forever
```

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --iterations 1000 \
  --fail-fast
```

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --suite spinchain-tebd \
  --test itdvp_projected_amplitude_tracks_time_dependent_imaginary_time_identity \
  --iterations 200
```

## Purpose

Use the ordinary runner and wrapper for normal integration coverage:

- [tests/run_suite.py](/home/ian/sync/git/main/tests/run_suite.py)
- [scripts/mptk-test](/home/ian/sync/git/main/scripts/mptk-test)

Use `mptk-stress-test` when:

- a failure is suspected to be intermittent
- a debug-only failure needs repeated attempts
- an eigensolver or orthogonalization path may depend on restart history
- a test is too expensive or too nondeterministic to add directly to the main
  green suite

This is intentionally a manual tool. It is useful for overnight runs and
focused debugging, not for ordinary `make test-integration`.

## Default Target

If no suite is given, the script stresses:

- [tests/suites/spinchain-aklt-spt-uc1-rerun.yaml](/home/ian/sync/git/main/tests/suites/spinchain-aklt-spt-uc1-rerun.yaml)

That suite exercises the warm-start rerun of single-site AKLT iDMRG, which is a
good probe of the infinite-wavefunction orthogonalization path.

## How It Works

The script repeatedly invokes the declarative runner on one suite or one named
test.

Per iteration it:

- runs [tests/run_suite.py](/home/ian/sync/git/main/tests/run_suite.py)
- writes the raw output to a per-iteration log file
- records pass/fail in a summary file
- deletes passing logs by default, so failed iterations stand out

At the end it prints:

- total iterations
- passed iterations
- failed iterations
- elapsed time
- the summary file location

If there are failures, it also prints the failed iteration log paths.

## Log Handling

By default the script writes logs under:

- `/tmp/mptk-stress-<suite>-<timestamp>/`

The directory contains:

- `summary.txt`
- `failures.txt`
- `iter-XXXX.log` for failed runs

Passing logs are removed unless `--keep-all-logs` is given.

If you want a stable location for a longer investigation, pass:

```bash
--log-dir /some/path
```

## Common Workflows

Run the default overnight debug stress test:

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --forever
```

Stop on the first failure and keep every log:

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --forever \
  --fail-fast \
  --keep-all-logs
```

Stress one specific test inside a suite:

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --suite spinchain-aklt-spt-uc1-rerun \
  --test rerunning_single_site_aklt_dmrg_preserves_exact_energy_and_spt_projective_phase \
  --iterations 500
```

Add runner verbosity for a shorter debugging session:

```bash
/home/ian/sync/git/main/scripts/mptk-stress-test \
  --bin-dir /home/ian/build/main-debug \
  --iterations 20 \
  --trace \
  --keep-all-logs
```

## Reading Results

Green outcome:

- all iterations pass
- no failed iteration logs remain unless `--keep-all-logs` was used

Interesting failure outcome:

- the script exits nonzero
- `summary.txt` records which iterations failed
- the failing `iter-XXXX.log` files remain for inspection

That gives a stable path from “I think this is random” to a concrete saved
reproducer.

## Relationship To The Main Suite

The stress wrapper does not replace the normal integration suite.

Use:

- [scripts/mptk-test](/home/ian/sync/git/main/scripts/mptk-test) for regular
  development and pre-push checks
- [scripts/mptk-stress-test](/home/ian/sync/git/main/scripts/mptk-stress-test)
  for repeated attempts at one suspicious case

That separation is intentional: the main suite should stay deterministic and
fast enough to run routinely, while the stress tool exists to explore edge
cases that need many repetitions.
