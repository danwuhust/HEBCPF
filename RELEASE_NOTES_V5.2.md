# HEBCPF MATLAB Solver Suite v5.2.0

- Version: `5.2.0`
- Git tag: `v5.2.0`
- Release date: 2026-07-30

HEBCPF v5.2.0 is the MATLAB-only successor to the public v5.0.0 release. The
Python implementations now live in the separate
[`pyHEBCPF`](https://github.com/danwuhust/pyHEBCPF) repository.

## Included solvers

- `HEBCPF_MEX_v5.2`: Windows x64 MEX-accelerated MATLAB solver.
- `HEBCPF_matlab_v5.2`: portable pure-MATLAB solver.

The two packages use the same HEBC formulation, MATPOWER cases, numerical
defaults, schedulers, checkpoints, and stopping rule. The MEX package compiles
the KLU corrector, Pade, and holomorphic hot paths.

## V5.2 changes

- Replaced the V5.0 dense cubic-memory quadratic setup with direct real sparse
  assembly in both MATLAB implementations.
- Normalized cases with MATPOWER `ext2int`, including non-consecutive external
  bus labels such as `case300`.
- Corrected PV generation handling by aggregating online generators by internal
  bus and validating shared voltage setpoints.
- Preserved equations for non-symmetric `Ybus` matrices, including
  phase-shifting transformers.
- Simplified `quadr_matrix` to return only the used outputs: `Ma`, `ra`, and
  `I`.
- Hardened MEX `parfor` checkpoint resume when older checkpoints lack the
  accumulated-time bookkeeping variables.
- Separated the MATLAB release from the Python suite and synchronized all
  README, notice, citation, checklist, and PDF documentation.
- Added a matched 20-case pure-MATLAB-versus-MEX comparison through 57 buses.

## Numerical evidence

The V5.2 comparison used queue `parfeval`, the `bandit` scheduler, a shared
23-worker pool, and three repeats per case on the recorded Intel Core
i9-13900K/128 GB/Windows 11 Pro x64 host with MATLAB R2022a. Both packages
returned the same solution count for every case, totaling 3,254 solutions.

From the rounded published table, pure MATLAB required 1,145.3 s in summed
search wall time and MEX required 586.3 s, a 1.95x aggregate factor. The
geometric-mean per-trace kernel factor is 1.97x. See
`MATLAB_BENCHMARK_V5.2.md` for all cases and limitations.

The raw repeated-run records and per-solution residual vectors for this new
comparison are not present, so v5.2.0 makes no new residual-distance claim.
The official V5 three-policy scheduler dataset remains separately documented
in `SCHEDULER_BENCHMARK_v5.md` and `BENCHMARK_METADATA_v5.json`.

## Dependencies and binary redistribution

MATLAB R2022a or later and MATPOWER 7.1 or a compatible release are expected.
MATPOWER is installed separately. Parallel modes require Parallel Computing
Toolbox. The MEX package includes Windows x64 binaries.

The shipped `klurf.mexw64` is the reproducible SuiteSparse v7.12.2 build used
in v5.0.0. The release assets include its matching binary, exact corresponding
source and license archive, benchmark evidence archive, and SHA-256 files.

## Documentation

- `HEBCPF_Suite_Overview.pdf`
- `HEBCPF_MEX_v5.2/HEBCPF_User_Guide.pdf`
- `HEBCPF_matlab_v5.2/HEBCPF_User_Guide.pdf`
- `MATLAB_BENCHMARK_V5.2.md`
- `POLICY_README.md`
- `SCHEDULER_BENCHMARK_v5.md`
