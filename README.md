# HEBCPF Solver Suite 2026.07.15

HEBCPF is a MATLAB solver suite for finding all real-valued solutions of AC
power-flow equations. It combines an ellipsoidal formulation, holomorphic
embedding, Pade approximation, and numerical continuation.

This release contains the v4 solver line in two forms:

| Folder | Solver | Best Use |
| --- | --- | --- |
| `HEBCPF_MEX_v4_20260715` | Windows x64 MEX-accelerated v4 | Fastest release on supported Windows MATLAB installations |
| `HEBCPF_matlab_v4_20260715` | Pure MATLAB v4 | Portable, inspectable, no compiler or MEX binaries required |

Both solvers use the same v4 numerical settings, including consensus Pade-pole
selection, deterministic keyed deduplication, `slope_max = 4e5`, and the v4
asynchronous `parfeval` work queue. The MEX release also includes a retained
barrier-style `parfeval` driver for A/B comparison.

Release 2026.07.15 is the benchmarked maintenance release of the v4 line. It
retains the 2026.07.14 formulation and numerical settings, while adding tested
checkpoint/resume behavior for long parallel searches. Older v2/v3 and v4
packages remain available through their GitHub tags for historical
reproduction; add only one solver folder to the MATLAB path at a time.

## Requirements

- MATLAB R2022a or later. This release was prepared and checked with MATLAB
  R2022a.
- MATPOWER 7.1 or a compatible MATPOWER release on the MATLAB path.
- Parallel Computing Toolbox for `run_batch_par.m`, `run_batch_parfeval.m`,
  `runVBook_hybrid_parallel.m`, and `runVBook_hybrid_parfeval.m`.
- MATLAB Coder and a supported C compiler only if rebuilding MEX binaries.

HEBCPF calls MATPOWER's installed `runpf`, `makeYbus`, `idx_bus`, and
`idx_brch`. It does not redistribute MATPOWER itself.

## Benchmark Results

The released queue-`parfeval` benchmark covers all 20 bundled cases through
57 buses at nominal load, using MATLAB R2022a on Windows x64 with 23 local
workers. Parallel-pool startup was excluded from per-case timing.

| Solver | Cases | Total wall time | Trace wall time | Total solutions |
| --- | ---: | ---: | ---: | ---: |
| MEX v4 | 20 | 718.656 s | 716.062 s | 3256 |
| Pure MATLAB v4 | 20 | 1272.966 s | 1271.222 s | 3256 |

MEX v4 used 1.77x less aggregate wall time and was faster on 14 cases. The
pure MATLAB package was faster on six small cases, where overhead dominates.
Every case had matching solution counts; an order-independent solution-set
check had maximum nearest-solution distance `3.154e-7`, below the shared
deduplication tolerance of `4e-7`. These figures compare the two 2026.07.15
v4 packages only, not earlier releases or other hardware/platforms.

Run `run_benchmark_up_to_57bus_20260715` to reproduce the experiment. The
committed report and CSV are in
`benchmark_results_20260715/20260715_194046/`; full per-case results appear
in `HEBCPF_Suite_Overview.pdf`.

## Quick Start

Use exactly one solver folder on the MATLAB path at a time because the MEX and
pure MATLAB folders contain many functions and cases with the same names.

```matlab
cd('<path_to_HEBCPF>/HEBCPF_MEX_v4_20260715')       % or HEBCPF_matlab_v4_20260715
addpath(pwd)

% one case, recommended large-case mode
[result, solutions] = run_merged_case('case14mod', 'parfeval');
result
```

Available modes are:

- `serial`: no Parallel Computing Toolbox; easiest to debug.
- `parfor`: synchronous parallel batch per starting solution.
- `parfeval`: asynchronous work queue; recommended for larger systems.

Batch entry points are `run_batch.m`, `run_batch_par.m`, and
`run_batch_parfeval.m`.

## Resuming A Ceased Search

The v4 runbook scripts write periodic checkpoints to `temp_result.mat` during
long `parfor` and `parfeval` runs. If MATLAB, the machine, or a long job stops
after this file has appeared, resume from the same solver folder:

```matlab
cd('<path_to_HEBCPF>/HEBCPF_MEX_v4_20260715')       % same folder used originally
addpath(pwd)

% First run main.m or equivalent preprocessing for the same case/load setting.
% Then load the saved solver state and run the same collection driver.
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval.m')                   % or runVBook_hybrid_parallel.m
```

For the serial driver, save the workspace manually before stopping:

```matlab
save('my_checkpoint.mat', '-v7.3')
```

then resume with:

```matlab
load my_checkpoint.mat
wait = 0;
run('runVBook_hybrid_2023.m')
```

Resume checkpoints are framework-specific. Resume a `parfeval` checkpoint with
`runVBook_hybrid_parfeval.m`, a `parfor` checkpoint with
`runVBook_hybrid_parallel.m`, and a serial checkpoint with
`runVBook_hybrid_2023.m`. The case, load factor, and solver folder must match
the original run.

## Documentation

- `HEBCPF_Suite_Overview.pdf`: release-level orientation and solver comparison.
- `HEBCPF_MEX_v4_20260715/HEBCPF_User_Guide.pdf`: MEX solver guide.
- `HEBCPF_matlab_v4_20260715/HEBCPF_User_Guide.pdf`: pure MATLAB solver guide.
- `HEBCPF_MEX_v4_20260715/README_v4.md`: concise MEX v4 notes.
- `HEBCPF_matlab_v4_20260715/README_v4_pure.md`: concise pure MATLAB v4 notes.

## Outputs

`run_merged_case` returns a result struct and the solution matrix. The result
includes the case name, selected mode, number of solutions, maximum algebraic
residual, and trace wall-clock time. Batch drivers write CSV summaries in the
run folder:

- `results_summary.csv`
- `results_par_summary.csv`
- `results_parfeval_summary.csv`

Generated checkpoint files, temporary MATLAB outputs, and CSV timing summaries
are ignored by the release `.gitignore`.

## Release Notes

See `CHANGELOG.md`, `RELEASE_NOTES_20260715.md`, and
`RELEASE_CHECKLIST_20260715.md`. This release records the completed v4
benchmark and adds tested checkpoint/resume support for long searches. The
queue scheduler, cached holomorphic operations, and lower-allocation Pade
evaluation are inherited v4 features.

## License, Attribution, and Citation

HEBCPF is distributed under the BSD 3-Clause License. See `LICENSE`.
MATPOWER and SuiteSparse/KLU remain third-party software under their own terms;
see `NOTICE` and the individual source notices.

If you use HEBCPF, cite the HEBC paper and, when appropriate, this software
package through `CITATION.cff`.
