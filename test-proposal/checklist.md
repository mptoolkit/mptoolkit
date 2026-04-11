# Coverage Checklist

This document is the concrete rollout checklist for extending the declarative
test prototype into broad MPToolkit integration coverage.

Scope:

- cover the supported toolkit and model surface in this repository
- eventually cover every model generator in `models/`
- explicitly do **not** treat `models/contrib/` as required coverage

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

## Framework Checklist

- [x] Raw action commands in suites use real CLI syntax.
- [x] Shared probe recipes live outside the runner.
- [x] Fixture outputs can be bound explicitly in `use:`.
- [x] Fixture materialization copies outputs into test-local directories.
- [x] Fixture dependency cycles are detected.
- [x] Action outputs are enforced and undeclared created files are errors.
- [x] String command tokenization is placeholder-safe.
- [x] `approx`, `equals`, `matches`, `contains`, and `exists` comparisons work.
- [ ] Add first-class line-pattern extraction with typed placeholders.
- [ ] Add persistent content-addressed fixture cache.
- [ ] Add fixture-build locking for parallel workers.
- [ ] Add parallel test execution.
- [ ] Add better failure summaries for missing/extra output files.
- [ ] Add first-class binding syntax for non-state outputs in examples and docs.
- [ ] Add suite tags or tiers for fast PR vs slow/nightly coverage.

## Shared Fixture Library Checklist

- [ ] Add a common lattice-fixture catalog for representative models.
- [ ] Add reusable finite product-state fixtures.
- [ ] Add reusable infinite product-state fixtures.
- [ ] Add reusable random-state fixtures with fixed seeds.
- [ ] Add reusable small DMRG ground-state fixtures.
- [ ] Add reusable small TEBD / iTEBD evolved-state fixtures.
- [ ] Add reusable TDVP / iTDVP evolved-state fixtures.
- [ ] Add normalization and canonicalization derived fixtures.
- [ ] Add a small catalog of exact tiny-system references where feasible.

## Shared Probe / Helper Checklist

- [x] `norm`
- [x] `expectation`
- [x] `expectation_real`
- [x] `expectation_imag`
- [x] `attr_float`
- [x] `attr_text`
- [x] `imoments_cross_json`
- [ ] `overlap`
- [ ] `ioverlap`
- [ ] `history_field`
- [ ] `info_field`
- [ ] `moments`
- [ ] `imoments`
- [ ] `json_field` as a generic helper
- [ ] typed line-pattern extractors for verbose text tools

## Tool-Family Checklist

### Construction And Basic Inspection

- [x] `mp-construct`
- [x] `mp-random`
- [x] `mp-norm`
- [x] `mp-expectation`
- [x] `mp-attr`
- [ ] `mp-overlap`
- [ ] `mp-info`
- [ ] `mp-history`
- [ ] `mp-lattice-info`
- [ ] `mp-show-operator`

### State Transforms And Canonicalization

- [ ] `mp-normalize`
- [ ] `mp-scale`
- [ ] `mp-conj`
- [ ] `mp-reflect`
- [ ] `mp-left-canonicalize`
- [ ] `mp-right-canonicalize`
- [ ] `mp-reorder-symmetry`
- [ ] `mp-change-lattice`
- [ ] `mp-irepeat`
- [ ] `mp-irotate`

### Finite Ground-State Algorithms

- [x] `mp-dmrg`
- [ ] `mp-dmrg-2site`
- [ ] `mp-dmrg-3s`

### Finite Time Evolution

- [x] `mp-tebd`
- [ ] `mp-tdvp`

### Infinite Ground-State And Evolution

- [ ] `mp-icdmrg`
- [ ] `mp-idmrg-ee`
- [ ] `mp-idmrg-s3e`
- [x] `mp-itebd`
- [ ] `mp-itdvp`

### Moments, Cross, And Overlap Families

- [ ] `mp-overlap`
- [ ] `mp-ioverlap`
- [ ] `mp-imoments`
- [x] `mp-imoments-cross`
- [ ] `mp-iexpectation`
- [ ] `mp-iexpectation-cross`
- [ ] `mp-ies`
- [ ] `mp-ies-cross`

### Correlation / Spectrum / Fluctuation Families

- [ ] `mp-allcorrelation`
- [ ] `mp-fluctuation`
- [ ] `mp-icorrelation`
- [ ] `mp-ibc-correlation`
- [ ] `mp-ispectrum`

### Apply / Matrix / Operator Families

- [ ] `mp-apply`
- [ ] `mp-iapply`
- [ ] `mp-ibc-apply`
- [ ] `mp-idivide`
- [ ] `mp-matrix`
- [ ] `mp-aux-matrix`
- [ ] `mp-aux-algebra`
- [ ] `mp-wigner-eckart`

### IBC Families

- [ ] `mp-ibc-create`
- [ ] `mp-ibc-dmrg`
- [ ] `mp-ibc-overlap`
- [ ] `mp-ibc-splice`
- [ ] `mp-ibc-tdvp`
- [ ] `mp-ibc-wavepacket`

### Excitations / EA Families

- [ ] `mp-excitation-ansatz`
- [ ] `mp-ea-create`
- [ ] `mp-ea-dmrg`
- [ ] `mp-ea-extend`
- [ ] `mp-ea-change-k`
- [ ] `mp-ea-moments`

### Specialized / Lower Priority

- [ ] `mp-coarsegrain`
- [ ] `mp-finegrain`
- [ ] `mp-finite-create`

## Representative Model-Family Checklist

These are the first `Level B` targets, chosen to give symmetry and statistics
coverage without exploding the matrix.

### Spin

- [ ] `spinchain-u1`
- [ ] `spinchain-su2`
- [ ] `spinchain-z2`
- [ ] `spinchain-spin2-su2`
- [ ] `spinladder.cpp`
- [ ] `spinladder-su2.cpp`
- [ ] `spinladder-z2.cpp`
- [ ] `spincylinder.cpp`
- [ ] `spincylinder-u1.cpp`
- [ ] `spincylinder-su2.cpp`
- [ ] `spincylinder-z2.cpp`
- [ ] `spinorbitchain.cpp`
- [ ] `spinorbitchain-u1u1.cpp`

### Fermion

- [ ] `spinlessfermion-u1.cpp`
- [ ] `spinlessfermionladder-u1u1.cpp`
- [ ] `hubbard.cpp`
- [ ] `hubbard-u1u1.cpp`
- [ ] `hubbard-su2.cpp`
- [ ] `hubbard-so4.cpp`
- [ ] `hubbard-u1su2.cpp`
- [ ] `hubbardcylinder-u1su2.cpp`
- [ ] `hubbardcylinder-u1su2-k.cpp`
- [ ] `hubbardcylinder-u1su2-k2.cpp`
- [ ] `klm-u1u1.cpp`
- [ ] `klm-u1su2.cpp`
- [ ] `klm-cylinder-u1su2.cpp`

### Boson

- [ ] `bosehubbard.cpp`
- [ ] `bosehubbard-u1.cpp`
- [ ] `bosehubbard-u1-trap.cpp`
- [ ] `bosehubbard-ladder-u1.cpp`
- [ ] `bosehubbard-flux-2leg-u1.cpp`
- [ ] `bosehubbard-flux-3leg-u1.cpp`
- [ ] `bosehubbard-2component-u1u1.cpp`
- [ ] `bosehubbard-2component-u1z2.cpp`

## Full Top-Level Model Inventory Checklist

Every file directly under `models/` should eventually get at least `Level A`
coverage. This list is the completion target for supported models.

- [ ] `bosehubbard-2component-u1u1.cpp`
- [ ] `bosehubbard-2component-u1z2.cpp`
- [ ] `bosehubbard-flux-2leg-u1.cpp`
- [ ] `bosehubbard-flux-3leg-u1.cpp`
- [ ] `bosehubbard-ladder-u1.cpp`
- [ ] `bosehubbard-u1-trap.cpp`
- [ ] `bosehubbard-u1.cpp`
- [ ] `bosehubbard.cpp`
- [ ] `hubbard-so4.cpp`
- [ ] `hubbard-su2.cpp`
- [ ] `hubbard-u1su2.cpp`
- [ ] `hubbard-u1u1.cpp`
- [ ] `hubbard.cpp`
- [ ] `hubbardcylinder-u1su2-k.cpp`
- [ ] `hubbardcylinder-u1su2-k2.cpp`
- [ ] `hubbardcylinder-u1su2.cpp`
- [ ] `klm-cylinder-u1su2.cpp`
- [ ] `klm-u1su2.cpp`
- [ ] `klm-u1u1.cpp`
- [ ] `spinchain-spin2-su2.cpp`
- [ ] `spinchain-su2.cpp`
- [ ] `spinchain-u1.cpp`
- [ ] `spinchain-z2.cpp`
- [ ] `spinchain.cpp`
- [ ] `spincylinder-su2.cpp`
- [ ] `spincylinder-u1.cpp`
- [ ] `spincylinder-z2.cpp`
- [ ] `spincylinder.cpp`
- [ ] `spinladder-su2.cpp`
- [ ] `spinladder-z2.cpp`
- [ ] `spinladder.cpp`
- [ ] `spinlessfermion-u1.cpp`
- [ ] `spinlessfermionladder-u1u1.cpp`
- [ ] `spinorbitchain-u1u1.cpp`
- [ ] `spinorbitchain.cpp`

## Suggested Rollout Order

### Phase 1: Finish The Core Contract Layer

- [ ] Add `overlap`, `ioverlap`, `info_field`, and `history_field` probes.
- [ ] Add transform and canonicalization suites.
- [ ] Add finite `mp-tdvp` suite.
- [ ] Add infinite `mp-itdvp` suite.

### Phase 2: Expand Model Diversity

- [ ] Add one spin `su2` suite.
- [ ] Add one fermion `u1` suite.
- [ ] Add one fermion `su2` suite.
- [ ] Add one boson `u1` suite.
- [ ] Add one ladder or cylinder suite.

### Phase 3: Add Exact / Stronger References

- [ ] Add tiny exact-reference checks for selected finite chains.
- [ ] Add regression tests for every new bug fix in core workflows.
- [ ] Add cross-tool consistency checks between norm/overlap/expectation/info.

### Phase 4: Broaden Coverage

- [ ] Bring all top-level `models/*.cpp` to `Level A`.
- [ ] Bring one representative in each family to `Level B`.
- [ ] Promote the most stable and informative suites into fast PR CI.
- [ ] Move broader and slower coverage into nightly CI.

## Explicit Non-Goals

- Do not require full coverage of `models/contrib/`.
- Do not require exhaustive tool × model Cartesian-product coverage.
- Do not depend on hidden action-specific adapters in the runner.
