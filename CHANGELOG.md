# Changelog

All notable changes to HEBCPF are documented here.

## 2026.07.14

- Added `HEBCPF_MEX_v4_20260714` and `HEBCPF_matlab_v4_20260714`, performance-optimized
  successors to the MEX and pure-MATLAB v3 releases.
- Replaced the MEX v3 `parfeval` row barrier with a global work queue that preserves
  one in-flight trace per equation; retained the old scheduler as
  `runVBook_hybrid_parfeval_barrier.m`.
- Added cached sparse LU ordering, sparse holomorphic-assembly reuse, and allocation-free
  Pad\'e polynomial evaluation in MEX v4.
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
