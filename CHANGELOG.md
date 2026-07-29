# Changelog

All notable changes to HEBCPF are documented here.

## v5.2.0 - 2026.07.30

- Released the MATLAB suite as `HEBCPF_MEX_v5.2` and
  `HEBCPF_matlab_v5.2`; Python implementations are maintained separately in
  the `pyHEBCPF` repository.
- Replaced dense three-dimensional power-flow quadratic construction with a
  sparse formulation in both the MEX and pure-MATLAB solvers.
- Preserved the existing preprocessing function names and inputs so all serial,
  `parfor`, `parfeval`, and benchmark drivers share one implementation.
- Simplified `quadr_matrix` to return only the three outputs used by the suite:
  `Ma`, `ra`, and `I`; bus-type counts remain internal loop variables.
- Added MATPOWER internal-index normalization for non-consecutive external bus
  numbers, including the bundled `case300` data.
- Corrected PV injection and voltage targets to use bus-based online-generator
  aggregation rather than generator-row ordering.
- Added validation for missing generators, inconsistent voltage setpoints,
  unsupported bus types, and invalid reference-bus counts.
- Corrected quadratic construction for non-symmetric admittance matrices, such
  as networks containing phase-shifting transformers.
- Hardened the MEX `parfor` resume path by initializing missing accumulated
  timing state from older compatible checkpoints.
- Added a 20-case, three-repeat, 23-worker comparison of the pure-MATLAB and
  MEX packages through 57 buses. Both packages completed 3,254 solutions; the
  rounded summary gives a 1.95x aggregate wall-time ratio and a 1.59x
  geometric-mean per-case wall-time factor in favor of MEX.
- Published the V5.2 comparison table and machine-readable metadata while
  explicitly recording that raw repeat records and residual vectors are not
  available for this comparison.
- Updated the suite overview, both user guides, citation metadata, release
  notes, checklist, and third-party notices for a MATLAB-only release.

## v5.1 - Unreleased development line

- V5.1 was not published as a HEBCPF MATLAB release. Its sparse preprocessing
  work is incorporated into v5.2.0. The Python ports developed alongside it
  are released separately as `pyHEBCPF` v5.2.0.

## v5.0.0 - 2026.07.26

V5.5 was an unreleased development package. V5 is documented directly against
the public 2026.07.15 v4 baseline.

- Packaged the live solvers as `HEBCPF_MEX_v5` and `HEBCPF_matlab_v5`.
- Rebuilt `klurf.mexw64` from pinned official SuiteSparse v7.12.2 and added a
  matching corresponding-source, license, build-metadata, and checksum package.
- Added three canonical queue schedulers: v4-compatible `scan`, learned
  productivity `bandit` (default), and bounded projected-distance `novelty`.
  The internal-v5 name `diverse` remains a deprecated alias for `bandit`.
- Replaced auxiliary per-equation pending lists with coverage cursors, reducing
  scheduler bookkeeping from `O(ns*neq)` duplicated lists to `O(neq)` state;
  the required `VBook` coverage matrix remains `O(ns*neq)`.
- Added per-trace anytime instrumentation, scheduler-overhead counters, and a
  resumable three-repeat benchmark of all bundled cases through case57.
- Synchronized scheduler behavior, configuration, wrappers, and batch defaults
  across the MEX and pure-MATLAB packages.
- Centralized the initial/hard solution-encounter caps at 2,000/5,000 through
  `solver_params.m` in every parallel driver.
- Extended checkpoints with policy, learned gains, and trace histories, and
  added a harness that exercises cyclic cross-policy transitions in both packages.
- Added explicit provenance metadata, guarded report publication, and a
  portable raw benchmark archive with a SHA-256 checksum for official dataset
  `20260724_3x23`.
- Reworked the suite overview and both user guides to compare V5 directly
  with v4 and to describe scheduler principles, complexity, and measured
  numerical behavior.
- Adopted official dataset `20260724_3x23` as the V5 scheduler evidence because
  the benchmarked V5.5 development tree and V5 solver source are identical;
  metadata records both the benchmark-time and rebuilt release KLU hashes.

## 2026.07.15

This benchmarked maintenance release packages the v4 solver line:

- Added `HEBCPF_MEX_v4_20260715`, the Windows x64 MEX-accelerated v4 solver.
- Added `HEBCPF_matlab_v4_20260715`, the portable pure MATLAB v4 solver.
- Added tested checkpoint/resume support for long parallel searches through
  `temp_result.mat`.
- Added a completed 20-case queue-`parfeval` benchmark through 57 buses,
  including solution-set cross-checks between the MEX and pure MATLAB v4
  packages.
- Reworked the suite overview and both user guides with the benchmark results
  and an explicit distinction from the 2026.07.14 v4 release.

## 2026.07.14

- Added `HEBCPF_MEX_v4_20260714` and `HEBCPF_matlab_v4_20260714`, performance-optimized
  successors to the MEX and pure-MATLAB v3 releases.
- Replaced the MEX v3 `parfeval` row barrier with a global work queue that preserves
  one in-flight trace per equation; retained the old scheduler as
  `runVBook_hybrid_parfeval_barrier.m`.
- Added cached sparse LU ordering, sparse holomorphic-assembly reuse, and allocation-free
  Pade polynomial evaluation in MEX v4.
- Added the same queue scheduler and holomorphic-kernel optimizations to MATLAB v4; its
  correctors deliberately use MATLAB built-in linear algebra only.
- Standardized `slope_max = 4e5` in the v3 solver parameter files and updated the guides.
- Updated the Suite Overview, README, citation metadata, notices, and User Guides for both
  v4 releases.

## 2026.07.12

Compared with the publicly released `HEBCPF` baseline:

- Added sparse-core Schur-complement bordered Newton solves with optional
  SuiteSparse/KLU symbolic-once refactorization through `klurf`.
- Added `klurf` source, build helper, documentation, and Windows x64 MEX
  binaries for the MEX releases.
- Added keyed bucket solution collection to the serial, `parfor`, and
  `parfeval` workflows for efficient long and resumed runs.
- Reduced solution-collection memory copying through geometric `Zsave` growth.
- Added `MakeJacobianD` sparse-core and cache initialization paths used by the
  KLU corrector.
- Updated the Suite Overview and all User Guides for package 2026.07.12.

## Public Baseline

The preceding public package used the `HEBCPF_matlab_v2.202607`,
`HEBCPF_matlab_v3.202606`, `HEBCPF_MEX_v2.202607`, and
`HEBCPF_MEX_v3.202606` release folders.
