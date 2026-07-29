# HEBCPF Pure MATLAB v5.2

This folder is the portable, no-compiled-code counterpart of
`HEBCPF_MEX_v5.2`. It requires MATLAB R2022a or later and MATPOWER on
the MATLAB path. Parallel modes additionally require the Parallel Computing
Toolbox.

## v4 Foundation And v5 Additions

- The v4 solver parameters, including `slope_max = 4e5`.
- Deterministic keyed deduplication.
- Consensus Pade-pole selection.
- Vectorized/cached construction in `holomorphic_para_sp.m`.
- Cached MATLAB `lu` column ordering through `lu_cached.m`.
- Lower-allocation Pade polynomial evaluation in `update_V_Pade.m`.
- The v4 queue-style `parfeval` execution model, now extended with `scan`,
  `bandit` (default), and bounded `novelty` schedulers.
- Periodic `temp_result.mat` checkpoints and resume from existing
  `VBook/Zsave` state.

## V5.2 Sparse Preprocessing

- `get_quadr_mtrx.m` builds sparse active/reactive quadratic forms directly
  from `real(Ybus)` and `imag(Ybus)`, without dense `2n`-by-`2n`-by-`n` arrays.
- `quadr_matrix.m` keeps the established output ordering while aggregating all
  online generators at each PV bus and validating common voltage setpoints.
- Every entry point uses MATPOWER `ext2int` numbering, so external bus labels do
  not need to be consecutive.
- The formulation supports non-symmetric admittance matrices, including
  phase-shifting transformer models.

## Pure MATLAB Guarantee

All execution paths use `.m` code and MATLAB built-ins. The correctors in
`branch_trace_hybrid_4_no_mex.m` and `Resolve_no_mex.m` do not call the MEX
hot path. The hot path calls `holomorphic_cont_tri.m` and `Pade_Apprxmt.m`,
never their MEX variants.

## Historical v4 Benchmark Result

In the retained 2026.07.15 queue-`parfeval` benchmark, the v4 predecessor completed all 20
cases through 57 buses in 1272.966 s, versus 718.656 s for the companion MEX
package. The pure MATLAB package was faster on six small cases, but MEX used
1.77x less aggregate wall time. Both packages returned the same solution count
for every case; the maximum order-independent nearest-solution distance was
`3.154e-7`, below the shared `4e-7` tolerance. The full table and conditions
are in `../HEBCPF_Suite_Overview.pdf`.

## V5.2 Pure MATLAB Versus MEX

The V5.2 packages were compared on all 20 bundled cases through `case57` with
queue `parfeval`, the `bandit` scheduler, three repeats per case, and a shared
persistent 23-worker pool. Both packages completed 3,254 solutions. From the
rounded table, pure MATLAB used 1,145.3 s in summed search wall time versus
586.3 s for MEX, while the geometric-mean case-wise wall-time factor was 1.59x
in favor of MEX.

See `../MATLAB_BENCHMARK_V5.2.md` for the complete table and environment. Raw
repeat records and residual vectors are unavailable for this comparison, so
it supports a solution-count agreement claim but no new residual-distance
claim.

## Run One Case

```matlab
cd('<path>/HEBCPF_matlab_v5.2')
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

The v5 scheduler implementation was checked against its MEX counterpart on
the smoke suite. The full 20-case scheduler benchmark and solution-set audit
are reported in `../SCHEDULER_BENCHMARK_v5.md`; all audited sets matched the
scan reference within the shared `4e-7` keyed-deduplication tolerance.
