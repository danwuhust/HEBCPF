# HEBCPF Pure MATLAB v4 20260715

This folder is the portable, no-compiled-code counterpart of
`HEBCPF_MEX_v4_20260715`. It requires MATLAB R2022a or later and MATPOWER on
the MATLAB path. Parallel modes additionally require the Parallel Computing
Toolbox.

## What Is Carried Over From MEX v4

- The v4 solver parameters, including `slope_max = 4e5`.
- Deterministic keyed deduplication.
- Consensus Pade-pole selection.
- Vectorized/cached construction in `holomorphic_para_sp.m`.
- Cached MATLAB `lu` column ordering through `lu_cached.m`.
- Lower-allocation Pade polynomial evaluation in `update_V_Pade.m`.
- The v4 queue-style `parfeval` scheduler with one in-flight trace per
  equation.
- Periodic `temp_result.mat` checkpoints and resume from existing
  `VBook/Zsave` state.

## Pure MATLAB Guarantee

All execution paths use `.m` code and MATLAB built-ins. The correctors in
`branch_trace_hybrid_4_no_mex.m` and `Resolve_no_mex.m` do not call the MEX
hot path. The hot path calls `holomorphic_cont_tri.m` and `Pade_Apprxmt.m`,
never their MEX variants.

## Benchmark Result

In the released queue-`parfeval` benchmark, this package completed all 20
cases through 57 buses in 1272.966 s, versus 718.656 s for the companion MEX
package. The pure MATLAB package was faster on six small cases, but MEX used
1.77x less aggregate wall time. Both packages returned the same solution count
for every case; the maximum order-independent nearest-solution distance was
`3.154e-7`, below the shared `4e-7` tolerance. The full table and conditions
are in `../HEBCPF_Suite_Overview.pdf`.

## Run One Case

```matlab
cd('<path>/HEBCPF_matlab_v4_20260715')
addpath(pwd)
[result, solutions] = run_merged_case('case14mod', 'parfeval');
```

Use `serial`, `parfor`, or `parfeval` as the second argument. `parfeval` is the
recommended mode for large cases. The driver creates `parpool('local')` when no
pool exists. For noninteractive batch runs, pass `true` as the third argument
to close the pool after the case:

```matlab
[result, solutions] = run_merged_case('case3', 'parfeval', true);
```

## Resume A Ceased Search

For `parfor` and `parfeval`, long runs periodically write `temp_result.mat`.
Resume from the same folder, same case, and same load factor:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval.m')          % or runVBook_hybrid_parallel.m
```

For the serial driver, save manually before stopping:

```matlab
save('my_checkpoint.mat', '-v7.3')
```

then resume:

```matlab
load my_checkpoint.mat
wait = 0;
run('runVBook_hybrid_2023.m')
```

The resume path rebuilds runtime-only futures, parallel pools, and keyed
deduplication state. Do not rename the periodic checkpoint inside the drivers;
it is intentionally `temp_result.mat`.

## Validation Notes

The pure MATLAB v4 solver was checked against MEX v4 across the completed
20-case benchmark through 57 buses. Its checkpoint/resume behavior was also
regression-tested on `case14mod` for the parallel and queue `parfeval` drivers.
Solution sets matched within the shared `4e-7` keyed-deduplication tolerance.
