# HEBCPF Solver Suite

HEBCPF is a MATLAB solver suite for finding all real-valued solutions of AC
power-flow equations. It combines an ellipsoidal formulation, holomorphic
embedding, Pade approximation, and numerical continuation.

This package release is **2026.07.14**. It updates the publicly released
`HEBCPF` baseline, preserves the v2/v3 solver lines, and adds MEX v4 and pure MATLAB v4.

## Releases

| Folder | Solver line | Orientation | Platform |
| --- | --- | --- | --- |
| `HEBCPF_MEX_v4_20260714` | v4.20260714 | Performance-optimized compiled build | Windows x64 MATLAB, or rebuild MEX |
| `HEBCPF_matlab_v4_20260714` | v4.20260714 | Performance-optimized pure MATLAB build | Cross-platform |
| `HEBCPF_MEX_v3.20260712` | v3.202606 | Robustness-oriented compiled build | Windows x64 MATLAB, or rebuild MEX |
| `HEBCPF_matlab_v3.20260712` | v3.202606 | Robustness-oriented pure MATLAB build | Cross-platform |
| `HEBCPF_MEX_v2.20260712` | v2.202607 | Efficiency-oriented compiled build | Windows x64 MATLAB, or rebuild MEX |
| `HEBCPF_matlab_v2.20260712` | v2.202607 | Efficiency-oriented pure MATLAB build | Cross-platform |

Use v3 for hard, large, meshed, or fold-dense systems. Use v2 when speed is
the priority on well-behaved systems. MEX builds provide precompiled Windows
x64 hot paths; pure MATLAB builds are portable and easy to inspect or modify.

See [HEBCPF_Suite_Overview.pdf](HEBCPF_Suite_Overview.pdf) for the complete
comparison and the `HEBCPF_User_Guide.pdf` inside each release folder for
installation, parameters, and troubleshooting.

## Requirements

- MATLAB R2022a or later. Package 2026.07.14 was validated with MATLAB R2022a.
- MATPOWER 7.1 installed and on the MATLAB path.
- Parallel Computing Toolbox for `run_batch_par.m` and `run_batch_parfeval.m`.
- MATLAB Coder and a supported C compiler only when rebuilding MEX functions.

HEBCPF calls MATPOWER's installed `runpf`, `makeYbus`, `idx_bus`, and
`idx_brch`. It does not redistribute a local `makeY.m` or MATPOWER's
`makeYbus.m`.

## Quick Start

1. Install MATLAB and MATPOWER, then confirm:

   ```matlab
   which runpf
   which makeYbus
   ```

2. Choose exactly one release folder. For example:

   ```matlab
   cd('<path_to_HEBCPF>/HEBCPF_matlab_v4_20260714')
   addpath(pwd)
   [result, solutions] = run_merged_case('case9', 'serial');
   result
   ```

3. For a single programmatic run:

   ```matlab
   r = run_merged_case('case14mod', 'parfeval');
   r.solutions
   r.max_residual
   r.wall_time
   ```

Do not add multiple HEBCPF release folders to the MATLAB path at the same time:
they contain functions and cases with overlapping names.

## KLU Acceleration

The bordered Newton corrector supports a sparse-core Schur-complement solve
with SuiteSparse/KLU symbolic-once refactorization through `klurf`.

- MEX releases include `klurf.mexw64` and enable it automatically on Windows x64.
- MATLAB v2 and MATLAB v3 can use a compatible `klurf` binary when one is placed in
  the release folder. MATLAB v4 deliberately uses MATLAB built-ins only and ignores
  `klurf`, even when it is present on the path.
- The `global USEKLU` override applies to MEX releases and optional-KLU MATLAB v2/v3:

  ```matlab
  USEKLU = 0;   % dense bordered solve
  USEKLU = 1;   % KLU full factorization; requires klu.mexw64
  USEKLU = 2;   % KLU symbolic-once refactor; requires klurf
  USEKLU = [];  % auto: klurf when available, otherwise dense
  ```

## Outputs and Validation

`run_merged_case` returns the number of solutions, the maximum algebraic
residual, and wall-clock time. The batch drivers write ignored CSV summaries:
`results_summary.csv`, `results_par_summary.csv`, and
`results_parfeval_summary.csv`.

The package was checked with MATLAB R2022a and MATPOWER 7.1. The documented MATLAB v4
validation record is in `HEBCPF_matlab_v4_20260714/README_v4_pure.md` and the Suite
Overview. For reproducible long runs, use the checkpoint instructions in the User Guide;
large MATLAB checkpoints should use `-v7.3`.

## Release Notes

See [CHANGELOG.md](CHANGELOG.md). Package 2026.07.14 adds MEX v4 and MATLAB v4,
including the v4 global `parfeval` queue and hot-path allocation/cache improvements.
The KLU bordered-Newton acceleration and keyed solution collection were introduced in
2026.07.12.

## License, Attribution, and Citation

HEBCPF is distributed under the BSD 3-Clause License. See [LICENSE](LICENSE).
MATPOWER and SuiteSparse/KLU remain third-party software under their own terms;
see [NOTICE](NOTICE) and the individual source notices.

If you use HEBCPF, cite the HEBC paper and, when appropriate, this software
package through [CITATION.cff](CITATION.cff).
