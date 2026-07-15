# HEBCPF 2026-07-15 benchmark summary

This folder contains the complete benchmark run for all solver folders present in the `HEBCPF_20260715` release folder:

- `HEBCPF_MEX_v4_20260715`
- `HEBCPF_matlab_v4_20260715`

The benchmark covers every bundled case up to and including 57 buses. The parallel pool was started once before the timed case loop, so the per-case timings below do not include parallel-pool startup or shutdown.

## Benchmark artifacts

- Combined CSV: `benchmark_summary_all.csv`
- Per-solver CSVs, per-case console logs, MATLAB snapshots, and the MATLAB record are
  retained locally for validation but are generated artifacts and are not committed to
  the GitHub release.

The committed CSV includes solver/case status, pool workers, case dimensions, staged
timings, trace timing, total case timing, solution counts, coverage status, solver
settings, and solution checksums/norms.

## Run status

All 40 solver-case runs completed successfully.

| Solver | Cases | Total case wall time (s) | Total trace wall time (s) | Slowest case | Slowest case time (s) |
| --- | ---: | ---: | ---: | --- | ---: |
| `MEX_v4_20260715` | 20 | 718.656 | 716.062 | `case57` | 366.985 |
| `matlab_v4_20260715` | 20 | 1272.966 | 1271.222 | `case57` | 693.529 |

## Per-case comparison

The speed ratio is `matlab_v4 total wall time / MEX_v4 total wall time`. Values above 1 mean the MEX solver was faster for that case. For very small cases, overhead dominates, so timing ratios should be interpreted cautiously.

| Case | Buses | Solutions | MEX v4 total (s) | MATLAB v4 total (s) | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| `case3` | 3 | 6 | 1.128 | 0.621 | 0.55x |
| `case3TS` | 3 | 6 | 0.313 | 0.332 | 1.06x |
| `case4BB0` | 4 | 14 | 1.023 | 0.767 | 0.75x |
| `case4BBc` | 4 | 12 | 2.025 | 0.409 | 0.20x |
| `case4gs` | 4 | 6 | 0.352 | 0.356 | 1.01x |
| `case5loop` | 5 | 10 | 0.595 | 0.834 | 1.40x |
| `case5Salam` | 5 | 10 | 0.656 | 0.703 | 1.07x |
| `case5Salam_mod3` | 5 | 4 | 0.424 | 0.307 | 0.72x |
| `case6ww` | 6 | 6 | 0.623 | 0.571 | 0.92x |
| `case7Salam` | 7 | 4 | 0.479 | 0.456 | 0.95x |
| `case9` | 9 | 8 | 0.904 | 0.922 | 1.02x |
| `case9Q` | 9 | 8 | 0.546 | 0.700 | 1.28x |
| `case14mod` | 14 | 30 | 3.110 | 4.495 | 1.45x |
| `case14mod2` | 14 | 68 | 6.134 | 8.640 | 1.41x |
| `case_ieee30` | 30 | 472 | 62.862 | 91.943 | 1.46x |
| `case30` | 30 | 472 | 61.776 | 93.156 | 1.51x |
| `case33bw` | 33 | 16 | 4.855 | 9.713 | 2.00x |
| `case39` | 39 | 176 | 34.870 | 58.694 | 1.68x |
| `case57` | 57 | 1322 | 366.985 | 693.529 | 1.89x |
| `case57mod` | 57 | 606 | 168.996 | 305.818 | 1.81x |

## Numerical agreement

Both packages returned the same solution count for every case. An
order-independent, one-to-one nearest-solution comparison of the saved snapshots gave a
maximum distance of `3.154e-7` over all 20 cases, below the shared `4e-7` keyed-deduplication
tolerance. This confirms agreement of the released MEX and pure MATLAB solution sets within
the package's accepted matching tolerance.

## Notes for reuse

- Use `benchmark_summary_all.csv` as the primary machine-readable benchmark table.
- Retain local per-case logs and snapshots when checking exact trace progression or
  reproducing numerical state; they are intentionally excluded from the GitHub payload.
- The benchmark driver is `../../run_benchmark_up_to_57bus_20260715.m` relative to this result folder.
