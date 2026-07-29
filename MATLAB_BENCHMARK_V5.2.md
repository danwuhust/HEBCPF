# HEBCPF V5.2 Pure MATLAB versus MEX Benchmark

This report preserves the matched V5.2 implementation comparison used in the
suite overview and both user guides. It compares only the pure-MATLAB and MEX
execution paths; it does not replace the official V5 scheduler benchmark in
`SCHEDULER_BENCHMARK_v5.md`.

## Environment and protocol

| Component | Recorded value |
| --- | --- |
| CPU | Intel Core i9-13900K, 24 cores / 32 threads |
| RAM | 128 GB |
| Operating system | Windows 11 Pro x64 |
| MATLAB | R2022a (9.12), `PCWIN64` |
| Parallel execution | Shared persistent local pool, 23 workers |
| Solver mode | Queue `parfeval` |
| Scheduler | `bandit` |
| Repeats | Three per case; table entries are averages |
| Pool startup | Excluded from every per-case time |
| Benchmark date | 2026-07-27 |

The two packages use the same equations, cases, tolerances, scheduler, and
stopping rule. The MEX package compiles the KLU corrector, Pade, and
holomorphic hot paths. Times are machine-specific and should not be compared
directly with the historical v4 or official V5 scheduler timings.

## Complete results

| Case | Buses | Solutions | Pure wall (s) | MEX wall (s) | Pure ms/trace | MEX ms/trace | Trace speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `case3` | 3 | 6 | 0.5 | 0.3 | 119.0 | 68.1 | 1.75x |
| `case3TS` | 3 | 6 | 0.2 | 0.1 | 53.4 | 23.0 | 2.33x |
| `case4gs` | 4 | 4 | 0.3 | 0.2 | 69.8 | 34.1 | 2.05x |
| `case7Salam` | 7 | 4 | 0.3 | 0.2 | 90.3 | 39.5 | 2.29x |
| `case6ww` | 6 | 6 | 0.4 | 0.2 | 74.5 | 33.6 | 2.22x |
| `case9` | 9 | 8 | 0.7 | 0.3 | 118.0 | 54.0 | 2.18x |
| `case9Q` | 9 | 8 | 0.5 | 0.3 | 100.0 | 46.6 | 2.15x |
| `case4BBc` | 4 | 12 | 0.3 | 0.6 | 63.5 | 44.4 | 1.43x |
| `case4BB0` | 4 | 14 | 0.5 | 0.3 | 84.3 | 43.3 | 1.95x |
| `case5loop` | 5 | 10 | 0.6 | 0.3 | 69.3 | 29.3 | 2.36x |
| `case5Salam` | 5 | 10 | 0.6 | 0.5 | 97.5 | 77.3 | 1.26x |
| `case5Salam_mod3` | 5 | 4 | 0.2 | 0.2 | 61.1 | 48.5 | 1.26x |
| `case14mod` | 14 | 30 | 3.1 | 2.2 | 168.0 | 116.0 | 1.46x |
| `case14mod2` | 14 | 68 | 5.3 | 3.8 | 143.0 | 82.2 | 1.74x |
| `case33bw` | 33 | 16 | 7.0 | 3.4 | 247.0 | 113.0 | 2.19x |
| `case39` | 39 | 176 | 50.7 | 28.0 | 200.0 | 93.2 | 2.15x |
| `case30` | 30 | 472 | 85.0 | 48.8 | 159.0 | 72.6 | 2.19x |
| `case_ieee30` | 30 | 472 | 88.1 | 48.6 | 167.0 | 72.0 | 2.31x |
| `case57mod` | 57 | 606 | 293.0 | 139.0 | 212.0 | 80.2 | 2.64x |
| `case57` | 57 | 1322 | 608.0 | 309.0 | 222.0 | 91.8 | 2.42x |

## Summary and interpretation

- Both implementations returned the same solution count for every case,
  totaling 3,254 solutions under the corrected V5.2 preprocessing.
- From the rounded table, summed wall time is 1,145.3 s for pure MATLAB and
  586.3 s for MEX, an aggregate 1.95x factor.
- The geometric-mean case-wise wall factor is 1.59x. The geometric-mean
  per-trace kernel factor is 1.97x.
- Small-case wall time includes scheduler and dispatch overhead; `case4BBc`
  therefore has a faster MEX trace kernel but a slower rounded end-to-end wall
  time. Larger cases show the compiled advantage more clearly.

## Evidence boundary

The preserved release artifact consists of the rounded table and chart in
`matlab_benchmark.tex` and `matlab_benchmark.png`. The raw three-repeat records,
per-solution residual vectors, and an execution commit identifier are not
present in this folder. Consequently, this release claims case-by-case solution
count agreement but does not infer a new residual-distance statistic. Future
replacement benchmarks should retain raw per-repeat CSV/MAT files, residual
summaries, the exact Git commit, dependency versions, and machine identity.
