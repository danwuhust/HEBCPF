# HEBCPF Solver Suite v5.0.0 (Historical Release Notes)

- Version: `5.0.0`
- Git tag: `v5.0.0`
- Release date: 2026-07-26

The historical v5.0.0 package contained the V5 solver line in two folders:

- `HEBCPF_MEX_v5`: Windows x64 MEX-accelerated solver.
- `HEBCPF_matlab_v5`: portable pure MATLAB solver.

The retained 2026.07.15 v4 figures provide the public predecessor baseline.
The V5 scheduler harness writes reproducible runs under
`scheduler_benchmark_v5/<timestamp>/`.

The published V5 numerical and timing tables use only official dataset
`20260724_3x23`, produced on the designated main job machine with MATLAB R2022a
on Windows x64 and 23 workers. Results from other machines are excluded. See
`BENCHMARK_METADATA_v5.json` for provenance and the separately attached
`HEBCPF-v5-benchmark-raw.zip` plus SHA-256 checksum for portable raw evidence.
The run used the source-equivalent V5.5 development tree; metadata records both
the benchmark-time and rebuilt-release KLU hashes. No result from the current
auxiliary machine is included.

The V5 Windows x64 `klurf.mexw64` was rebuilt from official SuiteSparse
`v7.12.2`, commit `42151688813c45846a597edcb601435a0e38f3dd`. The matching
release asset contains the exact component sources, licenses, build metadata,
smoke check, binary, and checksums. Its binary hash matches both copies shipped
in V5.

## Highlights

- Three queue schedulers: rotating `scan`, gain-directed `bandit` (default),
  and projection-based `novelty`.
- Per-trace anytime metrics, including traces and wall time to 50%, 90%, and
  99% of the final solution count.
- Tested checkpoint/resume support using `temp_result.mat`.
- Standardized `slope_max = 4e5`, inherited from the 2026.07.14 v4 release.
- Consensus Pade-pole strategy and deterministic keyed deduplication, inherited
  v4 numerical features.
- Cached sparse holomorphic operations and lower-allocation Pade evaluation,
  inherited v4 hot-path features.
- Updated user guides and suite overview with performance and connectivity
  plots for case14mod and case39.

## Historical v4 Benchmark Evidence

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

The full V5 report is `SCHEDULER_BENCHMARK_v5.md`, with
machine-readable summary tables beside it. The formulation remains the v4
formulation; V5 changes work selection and observability rather than the solved
equations.

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
