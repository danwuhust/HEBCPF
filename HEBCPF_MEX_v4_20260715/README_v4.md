# HEBCPF MEX v4 20260715

This folder is the Windows x64 MEX-accelerated v4 solver. It keeps the usual
HEBCPF workflows (`main.m`, `run_batch*.m`, `run_merged_case.m`) and uses the
same case files and solver parameters as the pure MATLAB v4 folder.

## What Is New In v4

| Area | Change | Main Files |
| --- | --- | --- |
| Scheduler | `runVBook_hybrid_parfeval.m` uses a global work queue with one in-flight trace per equation. The older row-barrier scheduler is retained as `runVBook_hybrid_parfeval_barrier.m`. | `runVBook_hybrid_parfeval.m`, `trace_equation_worker.m` |
| Holomorphic LU | Sparse `lu` column ordering is cached per case and reused. | `lu_cached.m`, `holomorphic_hybrid_5_with_mex.m` |
| Holomorphic assembly | Constant network blocks and stacked quadratic operators are reused. | `holomorphic_para_sp.m` |
| Pade evaluation | Voltage numerator/denominator polynomials are evaluated with lower allocation. | `update_V_Pade.m` |
| Resume support | Periodic checkpoints are written to `temp_result.mat`; `VBook/Zsave` checkpoints resume without re-seeding. | `runVBook_hybrid_parallel.m`, `runVBook_hybrid_parfeval.m`, `runVBook_hybrid_parfeval_barrier.m` |

Default numerical settings include deterministic keyed deduplication,
consensus Pade-pole selection, and `slope_max = 4e5`.

The queue scheduler and hot-path changes are inherited from the 2026.07.14
v4 release. This 2026.07.15 package adds the benchmarked checkpoint/resume
workflow; it does not change the underlying v4 formulation or defaults.

## Benchmark Result

In the released queue-`parfeval` benchmark, this package completed all 20
cases through 57 buses in 718.656 s, versus 1272.966 s for the companion pure
MATLAB package (1.77x less aggregate wall time). Both packages returned the
same solution count for every case. The maximum order-independent
nearest-solution distance was `3.154e-7`, below the shared `4e-7` tolerance.
The full table and conditions are in `../HEBCPF_Suite_Overview.pdf`.

## Run One Case

```matlab
cd('<path>/HEBCPF_MEX_v4_20260715')
addpath(pwd)
[result, solutions] = run_merged_case('case14mod', 'parfeval');
```

Use `serial`, `parfor`, or `parfeval` as the second argument. `parfeval` is the
recommended mode for large cases.

## Resume A Ceased Search

For `parfor`, `parfeval`, and the retained `parfeval` barrier driver, long runs
periodically write `temp_result.mat`. Resume from the same folder, same case,
and same load factor:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval.m')          % or runVBook_hybrid_parallel.m
```

For the retained barrier scheduler:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval_barrier.m')
```

The resume path rebuilds runtime-only futures, parallel pools, pool constants,
and keyed deduplication state. Do not rename the checkpoint file inside the
drivers; the periodic checkpoint name is intentionally `temp_result.mat`.

## MEX And KLU Notes

This release includes Windows x64 binaries:

- `Pade_Apprxmt_mex.mexw64`
- `holomorphic_cont_tri_mex.mexw64`
- `klurf.mexw64`

Run `build_mex.m` only if you change the hot-path source or need to rebuild for
your local MATLAB/compiler setup. The solver falls back to non-KLU dense solves
if `klurf` is unavailable.

## Validation Notes

Checkpoint/resume behavior was regression-tested by comparing full runs
against forced checkpoint/resume runs on `case14mod` and `case39` for the MEX
parallel, queue `parfeval`, and barrier `parfeval` drivers. The completed
release benchmark also covers all 20 bundled cases through 57 buses. Solution
sets matched within the shared `4e-7` keyed-deduplication tolerance.
