# HEBCPF Solver Suite

HEBCPF is a MATLAB solver suite for finding all real-valued solutions of the
AC power-flow equations using a holomorphic-embedding-based continuation
method. The implementation combines an ellipsoidal formulation, holomorphic
embedding, Pade approximation, and numerical continuation.

The suite contains four packaged releases. They implement the same method and
are intended to return the same solution set on the bundled test cases; the
choice between them is mainly about robustness, speed, and platform support.

## Releases

| Folder | Orientation | Platform |
| --- | --- | --- |
| `HEBCPF_MEX_v3.202606` | Robustness-oriented, compiled build | Windows x64 MATLAB, or rebuild MEX |
| `HEBCPF_matlab_v3.202606` | Robustness-oriented, pure MATLAB build | Cross-platform |
| `HEBCPF_MEX_v2.202607` | Efficiency-oriented, compiled build | Windows x64 MATLAB, or rebuild MEX |
| `HEBCPF_matlab_v2.202607` | Efficiency-oriented, pure MATLAB build | Cross-platform |

Use a v3 release when robustness is the priority, especially for hard, large,
meshed, or fold-dense systems. Use a v2 release when speed is the priority on
well-behaved systems. Use a MEX release on Windows x64 for faster hot-path
routines, or use a pure MATLAB release when portability and editability matter.

For more detail, see `HEBCPF_Suite_Overview.pdf` and the
`HEBCPF_User_Guide.pdf` file inside each release folder.

## Requirements

- MATLAB R2022a or later.
- MATPOWER installed and on the MATLAB path.
- Parallel Computing Toolbox, optional, only for `run_batch_par.m` and
  `run_batch_parfeval.m`.
- MATLAB Coder and a supported C compiler, optional, only if rebuilding the MEX
  files.

The setup stage uses MATPOWER functions such as `runpf`, `makeYbus`, `idx_bus`,
and `idx_brch`. Confirm MATPOWER is available before running HEBCPF:

```matlab
which runpf
which makeYbus
```

## Installation

1. Install MATLAB.
2. Install MATPOWER and add it to the MATLAB path.
3. Choose one HEBCPF release folder.
4. In MATLAB, change into that folder and add it to the path:

```matlab
cd('<path_to_HEBCPF>/HEBCPF_matlab_v3.202606')
addpath(pwd)
```

Do not add multiple HEBCPF release folders at the same time, because they share
many function and case names.

## Quick Start

From an interactive MATLAB session inside one release folder:

```matlab
addpath(pwd)
run_batch           % serial batch
run_batch_par       % parallel parfor, requires Parallel Computing Toolbox
run_batch_parfeval  % parallel parfeval, recommended for large cases
```

Run one case programmatically:

```matlab
r = run_merged_case('case14mod', 'parfeval');
r.solutions
r.max_residual
r.wall_time
```

Run from a system shell:

```powershell
matlab -sd "<path_to_release_folder>" -batch "addpath(pwd); run('run_batch.m')"
```

The parallel drivers start a local `parpool` automatically if none is running.
The one-time pool startup cost is excluded from the reported timing.

## Test Systems

The bundled MATPOWER-format cases include:

| Case | Real solutions |
| --- | ---: |
| `case3` | 6 |
| `case4gs` | 6 |
| `case4BB0` | 14 |
| `case4BBc` | 12 |
| `case5loop` | 10 |
| `case5Salam` | 10 |
| `case6ww` | 6 |
| `case7Salam` | 4 |
| `case9`, `case9Q` | 8 |
| `case14mod` | 30 |
| `case33bw` | 16 |
| `case39` | 176 |
| `case30`, `case_ieee30` | 472 |
| `case57mod` | 606 |
| `case57` | 1322 |
| `case118` | more than 1.5e6, partial enumeration |

Antipodal pairs are collapsed in the reported counts.

## Outputs

Each driver prints per-case progress, solution count, and timing information.
Batch runs also write CSV summaries in the release folder:

| Driver | CSV file |
| --- | --- |
| `run_batch.m` | `results_summary.csv` |
| `run_batch_par.m` | `results_par_summary.csv` |
| `run_batch_parfeval.m` | `results_parfeval_summary.csv` |

The full solution set is available in the MATLAB workspace as `Zsave`, with one
solution per column in the internal variable ordering. Solutions are
sign-normalized and deterministically ordered by `canonicalize_solutions.m` so
serial, `parfor`, and `parfeval` runs can be compared directly.

## MEX Builds

The MEX releases include Windows 64-bit `.mexw64` binaries for the hot-path
Pade and triangular-solve routines. On Windows x64, they are picked up
automatically.

On other platforms, or after editing a hot-path source file, rebuild inside a
MEX release folder:

```matlab
addpath(pwd)
build_mex
```

If rebuilding is not available, use the corresponding pure MATLAB release. The
pure MATLAB and MEX builds are intended to produce numerically identical
results.

## Validation

The packaged default `main.m` case in each of the four release folders has been
smoke-tested on MATLAB R2022a with MATPOWER 7.1. This is a startup and default
workflow check, not a substitute for full numerical validation on every case and
platform.

## File Map

- `main.m`: reference single-case preprocessing script.
- `run_batch.m`: serial batch driver.
- `run_batch_par.m`: `parfor` batch driver.
- `run_batch_parfeval.m`: `parfeval` batch driver.
- `run_merged_case.m`: programmatic single-case harness.
- `get_quadr_mtrx.m`, `quadr_matrix.m`, `MakeEllipse.m`, `Tune_Fac.m`: setup
  and ellipsoidal formulation utilities.
- `hybrid_traceloops_4_*`, `holomorphic_hybrid_5_*`,
  `branch_trace_hybrid_4_*`: branch-tracing and hybrid HEBC workflow.
- `Pade_Apprxmt.m`, `holomorphic_cont_tri.m`, `update_V_Pade.m`: holomorphic
  and Pade approximation routines.
- `canonicalize_solutions.m`: deterministic post-processing of solution sets.
- `cleanup_parallel_artifacts.m`: parallel-worker cache cleanup.
- `case*.m`: bundled MATPOWER-format test systems.
- `solver_params.m`: v3 release configuration source of truth.

## Troubleshooting

- `Undefined function 'runpf'`: MATPOWER is not on the MATLAB path.
- `Undefined function 'makeYbus'`: MATPOWER is not on the MATLAB path or the
  MATPOWER `lib` folder is missing from the path.
- `Undefined function 'caseXXX'`: run `addpath(pwd)` from the selected HEBCPF
  release folder.
- OSQP warning during `runpf`: this can occur during MATPOWER solver feature
  probing; `runpf` may still converge normally.
- Parallel driver appears to hang on first run: the local parallel pool is
  starting.
- MEX function missing or built for a different platform: rebuild with
  `build_mex`, or use a pure MATLAB release.
- Solution count differs from the expected value: check that the load factor is
  1 and that the initial MATPOWER power-flow solve reported success.

## Citation

If you use results produced by this solver, please cite:

```bibtex
@article{wu2019hebc,
  author = {Wu, Dan and Wang, Bin},
  title = {Holomorphic Embedding Based Continuation Method for Identifying Multiple Power Flow Solutions},
  journal = {IEEE Access},
  volume = {7},
  pages = {86843--86853},
  year = {2019},
  doi = {10.1109/ACCESS.2019.2925384}
}
```

Foundational and related works are listed in the user guides.

For citing the software package itself, see `CITATION.cff`.

## License and Third-Party Notices

The original HEBCPF code in this repository is distributed under the BSD
3-Clause License. See `LICENSE`.

This repository also includes MATPOWER-format case files copied from, derived
from, or modified from MATPOWER and other public test-system sources. MATPOWER
is required at runtime and is available from https://matpower.org and
https://github.com/MATPOWER/matpower. See `NOTICE` for attribution details.
