# Coverage Checklist

This document is the concrete rollout checklist for extending the declarative
integration test system into broad MPToolkit coverage.

Scope:

- cover the supported toolkit and model surface in this repository
- eventually cover every model generator in `models/`
- explicitly do **not** treat `models/contrib/` as required coverage
- explicitly do **not** require declarative integration coverage for binaries
  listed under `experimental-tools` in `Makefile.in` unless they are later
  promoted into the supported tool set

The point is not to cross every tool with every model. The point is to build a
test catalog that gives strong confidence in the high-level workflows while the
backend evolves.

## Coverage Levels

Use three coverage levels for model and tool planning.

- `Level A`: smoke / contract
  - generator runs
  - lattice file is usable
  - one simple state can be built
  - one basic probe works
- `Level B`: algorithm workflow
  - one representative ground-state or evolution workflow passes
  - output attributes and simple invariants are checked
- `Level C`: strong reference / regression
  - tiny exact reference where practical
  - or a focused regression for a known bug

The goal is:

- every top-level model in `models/` reaches `Level A`
- representative models in each family reach `Level B`
- selected core workflows reach `Level C`

Current baseline:

- `supported-model-contract` gives every supported top-level model at least
  `Level A` coverage for generator execution, lattice readability, symmetry /
  unit-cell metadata, and key operator names
- focused `Level C` suites cover selected exact tiny-system regressions such as
  Hubbard next-nearest hopping and spin-chain next-nearest exchange aliases

## Coverage Philosophy

The target is not exhaustive flag-permutation testing.

The target is:

- strong coverage of core workflow contracts
- broad smoke coverage across supported models
- explicit regression coverage for real bugs
- focused option coverage for distinct, high-risk behaviors

Guiding rules:

- test each important behavior class, not every syntactic spelling
- prefer model breadth over option explosion on one model
- prefer workflow-level invariants over tool-internal implementation detail
- add a regression test whenever a bug is fixed in a core workflow
- give extra attention to options that:
  - change semantics rather than formatting
  - have broken before
  - affect file mutation or output discovery
  - affect symmetry sector, normalization, or amplitude
  - change finite vs infinite behavior
  - provide machine-readable output

Usually not worth exhaustive direct coverage:

- cosmetic formatting flags
- redundant spellings for the same behavior
- every combination of otherwise independent options

When in doubt, ask:

- does this option create a distinct semantic contract?
- is it risky or historically fragile?
- would a user reasonably depend on it in production workflows?

If the answer is yes, it probably deserves direct coverage.

Tool-scope rule:

- treat the `tools` list in `Makefile.in` as the required integration-test
  surface
- treat the `experimental-tools` list as deferred by default
- if an experimental tool becomes supported, move it into the required list and
  then add coverage here

## High-Priority Risk Areas

These are not optional polish items. They are priority targets because errors
here can survive smoke tests while breaking real physics.

- complex conjugation and related transform behavior
- fermion-sign-sensitive workflows
- symmetry-sector preservation
- normalization and amplitude handling
- finite versus infinite algorithm differences
- output-file mutation and discovery behavior

Immediate implications:

- state-transform coverage is a priority, not a later cleanup task
- fermionic model coverage should start early, not after all spin-only suites
- exact or sharply constrained checks are especially valuable where sign errors
  can otherwise hide behind plausible scalar outputs

## Framework Checklist

- [x] Raw action commands in suites use real CLI syntax.
- [x] Shared probe recipes live outside the runner.
- [x] Fixture outputs can be bound explicitly in `use:`.
- [x] Fixture materialization copies outputs into test-local directories.
- [x] Fixture dependency cycles are detected.
- [x] Action outputs are enforced and undeclared created files are errors.
- [x] String command tokenization is placeholder-safe.
- [x] `approx`, `equals`, `matches`, `contains`, and `exists` comparisons work.
- [x] Add first-class line-pattern extraction with typed placeholders.
- [ ] Add persistent content-addressed fixture cache.
- [ ] Add fixture-build locking for parallel workers.
- [ ] Add parallel test execution.
- [ ] Add better failure summaries for missing/extra output files.
- [ ] Add first-class binding syntax for non-state outputs in examples and docs.
- [ ] Add suite tags or tiers for fast PR vs slow/nightly coverage.
- [x] Add a manual stress runner for intermittent-failure hunting.

Framework policy:

- use fresh fixture/test directories and unique output names rather than
  depending on `--force`
- treat inconsistent `--force` handling as a tool bug, not as required test
  harness behavior

Manual stress-test note:

- use [scripts/mptk-stress-test](../../scripts/mptk-stress-test)
  for overnight or repeated-attempt debugging of a single suspicious suite
- do not treat the stress runner as part of the normal green CI path
- keep ordinary integration coverage deterministic; use stress tests to turn
  suspected flakes into concrete saved failure logs

## Shared Fixture Library Checklist

- [ ] Add a common lattice-fixture catalog for representative models.
- [ ] Add reusable finite product-state fixtures.
- [ ] Add reusable infinite product-state fixtures.
- [ ] Add reusable random-state fixtures with fixed seeds.
- [ ] Add reusable small DMRG ground-state fixtures.
- [ ] Add reusable small TEBD / iTEBD evolved-state fixtures.
- [ ] Add reusable TDVP / iTDVP evolved-state fixtures.
- [ ] Add normalization and canonicalization derived fixtures.
- [ ] Add genuinely complex-state fixtures for conjugation and phase tests.
- [ ] Add fermion-sign-sensitive fixtures on representative fermionic models.
- [ ] Add a small catalog of exact tiny-system references where feasible.

## Shared Probe / Helper Checklist

- [x] `norm`
- [x] `expectation`
- [x] `expectation_real`
- [x] `expectation_imag`
- [x] `attr_float`
- [x] `attr_text`
- [x] `info_text`
- [x] `history_text`
- [x] `lattice_info_text`
- [x] `show_finite_operator_text`
- [x] `show_product_operator_text`
- [x] `imoments_cross_json`
- [x] `overlap`
- [x] `ioverlap`
- [ ] `history_field`
- [ ] `info_field`
- [ ] `moments`
- [x] `imoments`
- [ ] `json_field` as a generic helper
- [x] typed line-pattern extractors for verbose text tools

## Tool-Family Checklist

### Construction And Basic Inspection

- [x] `mp-construct`
- [x] `mp-random`
- [x] `mp-norm`
- [x] `mp-expectation`
- [x] `mp-attr`
- [x] `mp-overlap`
- [x] `mp-info`
- [x] `mp-history`
- [x] `mp-lattice-info`
- [x] `mp-show-operator`

### State Transforms And Canonicalization

- [x] `mp-normalize`
- [x] `mp-scale`
- [x] `mp-conj`
- [x] `mp-reflect`
- [x] `mp-left-canonicalize`
- [x] `mp-right-canonicalize`
- [x] `mp-reorder-symmetry`
- [x] `mp-change-lattice`
- [x] `mp-irepeat`
- [x] `mp-irotate`

### Finite Ground-State Algorithms

- [x] `mp-dmrg`
- [x] `mp-dmrg-2site`
- [x] `mp-dmrg-3s`

### Finite Time Evolution

- [x] `mp-tebd`
- [x] `mp-tdvp`

Finite-algorithm note:

- `mp-dmrg-2site` now has declarative coverage in both optimized and debug
  builds.

### Infinite Ground-State And Evolution

- [ ] `mp-icdmrg`
- [ ] `mp-idmrg-ee`
- [x] `mp-idmrg-s3e`
- [x] `mp-itebd`
- [x] `mp-itdvp`

### Moments, Cross, And Overlap Families

- [x] `mp-overlap`
- [x] `mp-ioverlap`
- [x] `mp-imoments`
- [x] `mp-imoments-cross`
- [ ] `mp-iexpectation` (obsolete; use `mp-expectation`)
- [x] `mp-iexpectation-cross`
- [x] `mp-ies`
- [x] `mp-ies-cross`

### Correlation / Spectrum / Fluctuation Families

- [x] `mp-allcorrelation`
- [x] `mp-icorrelation`
- [x] `mp-ibc-correlation`
- [x] `mp-ispectrum`

### Apply / Matrix / Operator Families

- [x] `mp-apply`
- [x] `mp-iapply`
- [x] `mp-ibc-apply`
- [x] `mp-idivide`
- [ ] `mp-matrix`
- [x] `mp-aux-matrix`
- [x] `mp-aux-algebra`
- [x] `mp-wigner-eckart`

### IBC Families

- [x] `mp-ibc-create`
- [ ] `mp-ibc-dmrg`
- [x] `mp-ibc-overlap`
- [ ] `mp-ibc-splice`
- [x] `mp-ibc-tdvp`
- [x] `mp-ibc-wavepacket`

### Excitations / EA Families

- [ ] `mp-excitation-ansatz`
- [x] `mp-ea-create`
- [ ] `mp-ea-dmrg`
- [x] `mp-ea-extend`
- [x] `mp-ea-change-k`
- [ ] `mp-ea-moments`

### Specialized / Lower Priority

- [x] `mp-coarsegrain`
- [x] `mp-finegrain`
- [x] `mp-finite-create`

### Deferred Experimental Tools

These binaries currently live in `experimental-tools` in `Makefile.in`, so
they are not part of the required declarative integration-test surface unless
they are promoted later.

- `mp-fluctuation`

## Known Gaps Inside Covered Tools

- [ ] `mp-change-lattice` should handle states that project completely out of
  the target local basis without aborting. A manual `bosehubbard-u1` check from
  `N=5` to `N=2` with the state `0:1:3:0` currently triggers a tensor-basis
  precondition failure instead of returning the zero state or a clear error.
- [ ] `mp-matrix` currently aborts in debug builds on simple `spinchain-u1`
  product iMPS inputs, so it is not yet covered declaratively even though
  `mp-aux-matrix` is.
- [ ] `mp-ibc-dmrg` currently segfaults on simple `mp-ibc-create` outputs,
  including identical-boundary `spinchain-u1` IBC states with both empty and
  one-site windows. The crash happens before the initial `Timestep=0` output,
  so the current bug is in the initialization/first-energy path rather than the
  declarative harness.
- [ ] `mp-ibc-splice` currently aborts on a trivial self-splice of identical
  `spinchain-u1` IBC states (`mp-ibc-splice ibc ibc -o out`) with a
  `pheap::pvalue_ptr` null-pointer precondition failure.
- [ ] `mp-ea-moments` currently aborts on the trivial self-overlap of a
  one-site `spinchain-u1` EA state created with `mp-ea-create "lat:Sm(0)"`,
  hitting the ARPACK precondition `n > 0 is false` before producing output.

## Representative Model-Family Checklist

These are the first `Level B` targets, chosen to give symmetry and statistics
coverage without exploding the matrix.

### Spin

- [x] `spinchain-u1`
- [x] `spinchain-su2`
- [x] `spinchain-z2`
- [x] `spinchain-spin2-su2`
  Includes AKLT ground-state and SPT projective-commutator checks.
- [x] `spinladder.cpp`
- [x] `spinladder-su2.cpp`
- [x] `spinladder-z2.cpp`
- [x] `spincylinder.cpp`
- [x] `spincylinder-u1.cpp`
- [x] `spincylinder-su2.cpp`
- [x] `spincylinder-z2.cpp`
- [x] `spinorbitchain.cpp` moved to `models/contrib/`
- [x] `spinorbitchain-u1u1.cpp` moved to `models/contrib/`

### Fermion

- [x] `spinlessfermion-u1.cpp` priority
- [x] `hubbard-u1u1.cpp` priority
- [x] `hubbard-u1su2.cpp` priority
- [x] `hubbard-so4.cpp` priority
- [x] `spinlessfermion-u1.cpp`
- [x] `spinlessfermionladder-u1u1.cpp`
- [x] `hubbard.cpp`
- [x] `hubbard-u1u1.cpp`
- [x] `hubbard-su2.cpp`
- [x] `hubbard-u1su2.cpp`
- [x] `hubbard-so4.cpp`
- [x] `hubbardcylinder-u1su2.cpp`
- [x] `hubbardcylinder-u1su2-k.cpp` moved to `models/contrib/`
- [x] `hubbardcylinder-u1su2-k2.cpp` moved to `models/contrib/`
- [x] `klm-u1u1.cpp`
- [x] `klm-u1su2.cpp`
- [x] `klm-cylinder-u1su2.cpp`

### Boson

- [x] `bosehubbard.cpp`
- [x] `bosehubbard-u1.cpp`
- [x] `bosehubbard-u1-trap.cpp`
- [x] `bosehubbard-ladder-u1.cpp`
- [x] `bosehubbard-flux-2leg-u1.cpp` moved to `models/contrib/`
- [x] `bosehubbard-flux-3leg-u1.cpp` moved to `models/contrib/`
- [x] `bosehubbard-2component-u1u1.cpp` moved to `models/contrib/`
- [x] `bosehubbard-2component-u1z2.cpp` moved to `models/contrib/`

### Discrete

- [x] `clock.cpp`
- [x] `clock-zn.cpp`
- [x] `potts.cpp`
- [x] `potts-zn.cpp`

## Full Top-Level Model Inventory Checklist

Every file directly under `models/` should eventually get at least `Level A`
coverage. This list is the completion target for supported models.

- [x] `bosehubbard-ladder-u1.cpp`
- [x] `bosehubbard-u1-trap.cpp`
- [x] `bosehubbard-u1.cpp`
- [x] `bosehubbard.cpp`
- [x] `clock-zn.cpp`
- [x] `clock.cpp`
- [x] `hubbard-so4.cpp`
- [x] `hubbard-su2.cpp`
- [x] `hubbard-u1su2.cpp`
- [x] `hubbard-u1u1.cpp`
- [x] `hubbard.cpp`
- [x] `hubbardcylinder-u1su2.cpp`
- [x] `klm-cylinder-u1su2.cpp`
- [x] `klm-u1su2.cpp`
- [x] `klm-u1u1.cpp`
- [x] `potts-zn.cpp`
- [x] `potts.cpp`
- [x] `spinchain-spin2-su2.cpp`
- [x] `spinchain-su2.cpp`
- [x] `spinchain-u1.cpp`
- [x] `spinchain-z2.cpp`
- [x] `spinchain.cpp`
- [x] `spincylinder-su2.cpp`
- [x] `spincylinder-u1.cpp`
- [x] `spincylinder-z2.cpp`
- [x] `spincylinder.cpp`
- [x] `spinladder-su2.cpp`
- [x] `spinladder-z2.cpp`
- [x] `spinladder.cpp`
- [x] `spinlessfermion-u1.cpp`
- [x] `spinlessfermionladder-u1u1.cpp`

## Suggested Rollout Order

### Phase 1: Finish The Core Contract Layer

- [x] Add all-supported-model contract coverage.
- [ ] Add `info_field` and `history_field` probes.
- [ ] Extend transform coverage with stronger complex-conjugation checks.
- [x] Add finite `mp-tdvp` suite.
- [x] Add infinite `mp-itdvp` suite.
- [x] Add the first fermion-sign-sensitive suite on `spinlessfermion-u1`.
- [x] Add the first Hubbard-family fermion-sign suite.

### Phase 2: Expand Model Diversity

- [x] Add one spin `su2` suite.
- [x] Add one fermion `u1` suite.
- [x] Add one fermion `su2` suite.
- [x] Add one boson `u1` suite.
- [x] Add one ladder or cylinder suite.

### Phase 3: Add Exact / Stronger References

- [x] Add tiny exact-reference checks for selected finite chains.
- [ ] Add regression tests for every new bug fix in core workflows.
- [ ] Add cross-tool consistency checks between norm/overlap/expectation/info.

### Phase 4: Broaden Coverage

- [x] Bring all top-level `models/*.cpp` to `Level A`.
- [ ] Bring one representative in each family to `Level B`.
- [ ] Promote the most stable and informative suites into fast PR CI.
- [ ] Move broader and slower coverage into nightly CI.

## Explicit Non-Goals

- Do not require full coverage of `models/contrib/`.
- Do not require exhaustive tool × model Cartesian-product coverage.
- Do not depend on hidden action-specific adapters in the runner.
