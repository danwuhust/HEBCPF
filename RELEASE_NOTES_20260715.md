# HEBCPF Solver Suite 2026.07.15

This benchmarked maintenance release packages the v4 solver line in two folders:

- `HEBCPF_MEX_v4_20260715`: Windows x64 MEX-accelerated solver.
- `HEBCPF_matlab_v4_20260715`: portable pure MATLAB solver.

## Highlights

- Queue-style `parfeval` driver with per-equation serialization, inherited from
  the 2026.07.14 v4 release.
- Tested checkpoint/resume support using `temp_result.mat`.
- Standardized `slope_max = 4e5`, inherited from the 2026.07.14 v4 release.
- Consensus Pade-pole strategy and deterministic keyed deduplication, inherited
  v4 numerical features.
- Cached sparse holomorphic operations and lower-allocation Pade evaluation,
  inherited v4 hot-path features.
- Updated user guides for both solver folders and a new v4 suite overview.

## Benchmark Evidence

The released queue-`parfeval` benchmark ran all 20 bundled cases through 57
buses at nominal load using MATLAB R2022a on Windows x64 with 23 local workers.
Parallel-pool startup was excluded from per-case times.

| Solver | Total wall time | Trace wall time | Solutions |
| --- | ---: | ---: | ---: |
| MEX v4 | 718.656 s | 716.062 s | 3256 |
| Pure MATLAB v4 | 1272.966 s | 1271.222 s | 3256 |

MEX v4 used 1.77x less aggregate wall time and was faster on 14 of 20 cases.
Solution counts matched case-by-case. The saved solution sets also matched
one-to-one with maximum nearest-solution distance `3.154e-7`, below the shared
`4e-7` deduplication tolerance. These are package-level results for this
environment, not a speed claim for every platform or execution mode.

The full report and machine-readable table are in
`benchmark_results_20260715/20260715_194046/`. The 2026.07.15 release retains
the 2026.07.14 v4 formulation and numerical defaults; it adds the benchmarked
checkpoint/resume workflow rather than a new mathematical algorithm.

## Resume Reminder

To resume a ceased parallel search, return to the same solver folder, load
`temp_result.mat`, and run the same driver:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval.m')     % or runVBook_hybrid_parallel.m
```

The MEX folder also keeps `runVBook_hybrid_parfeval_barrier.m` for comparison.
Use that same file when resuming a barrier-driver checkpoint.

## Validation

Checkpoint/resume behavior was regression-tested by comparing uninterrupted
runs against forced checkpoint/resume runs:

- MEX v4: `parallel`, queue `parfeval`, and barrier `parfeval` on `case14mod`
  and `case39`.
- Pure MATLAB v4: `parallel` and queue `parfeval` on `case14mod`.

All tested resumed runs matched the uninterrupted solution sets well below the
`4e-7` keyed-deduplication tolerance.

## Notes

- Add only one solver folder to the MATLAB path at a time.
- `parfeval` is the recommended mode for large cases.
- MATPOWER must be installed separately and available on the MATLAB path.
