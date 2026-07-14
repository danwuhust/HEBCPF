# HEBCPF pure MATLAB v4

This folder is the portable, no-compiled-code counterpart of
`HEBCPF_MEX_v4_20260714`. It requires MATLAB R2022a+ and MATPOWER 7.1 (or a
compatible MATPOWER release). The parallel modes additionally require the
Parallel Computing Toolbox.

## What is carried over from MEX v4

- The v4 solver parameters, including `slope_max = 4e5` and deterministic
  seven-variable keyed deduplication.
- Vectorized/cached construction in `holomorphic_para_sp.m`.
- Cached MATLAB `lu` column ordering through `lu_cached.m`.
- Allocation-free Padé polynomial evaluation in `update_V_Pade.m`.
- The v4 work-queue `parfeval` scheduler, with one in-flight trace per
  equation to preserve VBook pruning without row barriers.

## Pure-MATLAB guarantee

All execution paths use `.m` code and MATLAB built-ins. The correctors in
`branch_trace_hybrid_4_no_mex.m` and `Resolve_no_mex.m` use MATLAB backslash;
they do not probe for or call KLU/klurf. The hot path calls
`holomorphic_cont_tri.m` and `Pade_Apprxmt.m`, never their MEX variants.

## Run one case

```matlab
cd('<path>/HEBCPF_matlab_v4_20260714')
addpath(pwd)
[result, solutions] = run_merged_case('case9', 'serial');
```

Use `'parfor'` or `'parfeval'` as the second argument to run the corresponding
parallel driver. The driver creates `parpool('local')` when no pool exists.
For noninteractive batch runs, pass `true` as the third argument to close the
pool after the case:

```matlab
[result, solutions] = run_merged_case('case3', 'parfeval', true);
```

`SOLVER_USE_POOL_CONST=true` optionally wraps read-only data in
`parallel.pool.Constant`; the default is direct data transfer, chosen for
desktop MATLAB stability.

## Validation performed on 2026-07-14

| Run | Result |
|---|---|
| `case3`, serial | 6 solutions, max residual `1.69e-14` |
| `case3`, pure v4 vs MEX v4 | 6 vs 6; bidirectional set error `2.42e-14` L1 |
| `case9`, serial | 8 solutions, max residual `7.39e-13` |
| `case9`, pure v4 vs MEX v4 | 8 vs 8; bidirectional set error `2.89e-12` L1 |
| `case14mod`, serial | 30 solutions, max residual `9.25e-13` |
| `case14mod`, pure v4 vs MEX v4 | 30 vs 30; bidirectional set error `2.49e-12` L1 |
| `case3`, `parfeval` | pool started; 6 solutions, max residual `2.35e-14` |

The legacy PDF guide in this folder documents the HEBC method and general
workflow. This README is authoritative for pure v4-specific behavior and
dependencies.
