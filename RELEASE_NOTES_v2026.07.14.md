# HEBCPF Solver Suite v2026.07.14

This release adds `HEBCPF_MEX_v4_20260714` and `HEBCPF_matlab_v4_20260714` while
retaining the existing v2 and v3 solver folders.

## Highlights

- Added the MEX v4 global `parfeval` work queue with one in-flight trace per equation.
- Added cached sparse LU ordering, sparse holomorphic-assembly reuse, and allocation-free
  Pade polynomial evaluation in MEX v4.
- Added the same v4 queue scheduler and holomorphic-kernel optimizations to pure MATLAB v4,
  with MATLAB built-in linear algebra retained for the correctors.
- Retained the previous MEX v3 barrier scheduler as
  `runVBook_hybrid_parfeval_barrier.m` for controlled comparisons.
- Standardized `slope_max = 4e5` in the v3 parameter files.
- Updated the Suite Overview, User Guides, citation metadata, notice, and changelog.

## Requirements

- MATLAB R2022a or later
- MATPOWER 7.1
- Parallel Computing Toolbox for `parfor` and `parfeval`

## Notes

- MEX v4 includes Windows x64 MEX binaries, including `klurf.mexw64`.
- The reported MEX v3/v4 `parfeval` reference table is included in the v4 User Guide
  and Suite Overview. It was not re-executed during the package-preparation pass.
- See `CHANGELOG.md` and `HEBCPF_Suite_Overview.pdf` for details.
