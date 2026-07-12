# Changelog

All notable changes to HEBCPF are documented here.

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
