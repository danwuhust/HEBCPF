# HEBCPF MATLAB Solver Suite V5.2

HEBCPF is a MATLAB solver suite for enumerating real-valued solutions of AC
power-flow equations by holomorphic-embedding-based continuation (HEBC).

This release contains the two MATLAB implementations of the solver:

| Folder | Implementation | Best use |
| --- | --- | --- |
| `HEBCPF_MEX_v5.2` | Windows x64 MEX-accelerated V5.2 | Fastest supported Windows configuration |
| `HEBCPF_matlab_v5.2` | Pure MATLAB V5.2 | Portable, inspectable, and compiler-free |

The two builds use the same formulation, cases, tolerances, and scheduler
semantics. They returned the same completed solution count on every case in
the published V5.2 comparison; the MEX build compiles selected hot paths. Each
folder contains a full PDF user guide. **Python** implementations of the same
method (pure Python and native Cython) are released separately as the
`pyHEBCPF` v5.2 suite.

## V5.2 preprocessing update

V5.2 retains the V5 tracing, scheduler, checkpoint, and numerical settings but
replaces the dense cubic-memory power-flow matrix setup with sparse assembly.
All implementations construct real sparse quadratic forms directly from the
real and imaginary parts of `Ybus`. Cases are normalized to internal bus order;
active generation is summed across online generators at each PV bus, and shared
voltage setpoints are validated. This also supports non-symmetric `Ybus`
matrices arising from phase-shifting transformers.

V5.2 also adds a matched pure-MATLAB-versus-MEX comparison for all 20 completed
cases through 57 buses. That comparison is reported separately from the
official V5 scheduler dataset because it answers a different question: the
implementation speed difference with the scheduler, cases, and worker pool
held fixed.

## Direct comparison with v4

V5.5 was an unreleased development package and is not a public release
baseline. V5 is therefore compared directly with the 2026.07.15 v4 package:

| Area | v4 (2026.07.15) | V5 |
| --- | --- | --- |
| Queue scheduling | `scan` rotating fairness | `scan`, `bandit` (default), and `novelty` |
| Pending-pair state | Per-equation pending lists | Per-equation cursors; auxiliary scheduler state is `O(neq)` |
| Anytime instrumentation | Final counts and timing | Per-trace discovery, completion time, equation, and worker time |
| Solution encounter cap | Older fixed/default settings | Central 2,000 initial cap and 5,000 hard cap via `solver_params` |
| Checkpoint policy state | Solver state | Solver state, trace history, active policy, and learned gains |
| Benchmarking | MEX-versus-pure v4 benchmark | Three-scheduler, three-repeat V5 MEX benchmark through case57 |

The mathematical formulation, MATPOWER case format, convergence tolerances,
deduplication tolerance, and final completeness criterion remain compatible
with v4. Scheduler choice changes dispatch order and the anytime discovery
profile; it does not intentionally change the completed solution set.

## Requirements

MATLAB implementations:

- MATLAB R2022a or later. V5.2 was prepared with R2022a.
- MATPOWER 7.1 or a compatible release on the MATLAB path.
- Parallel Computing Toolbox for the `parfor` and `parfeval` drivers.
- MATLAB Coder and a supported compiler only when rebuilding MEX binaries.

HEBCPF calls MATPOWER's installed `runpf`, `makeYbus`, `ext2int`, `idx_bus`,
`idx_gen`, and `idx_brch`; MATPOWER is not redistributed in this package.

## MATLAB quick start

Add only one implementation folder to the MATLAB path because both folders
contain functions and cases with the same names.

```matlab
cd('<path_to_HEBCPF>/HEBCPF_MEX_v5.2')  % or HEBCPF_matlab_v5.2
addpath(pwd)

global HEBCPOLICY
HEBCPOLICY = 'bandit';                  % default; also scan or novelty
[result, solutions] = run_merged_case('case14mod','parfeval');
result
```

The legacy v5 name `diverse` is accepted as a deprecated alias for `bandit`.
See `POLICY_README.md` for the algorithms, complexity, checkpoint semantics,
and selection guidance.

Available execution modes are:

- `serial`: no Parallel Computing Toolbox; simplest to debug.
- `parfor`: synchronous parallel batch per starting solution.
- `parfeval`: asynchronous global work queue; the scheduler comparison applies
  to this mode.

Batch entry points are `run_batch.m`, `run_batch_par.m`, and
`run_batch_parfeval.m`.

For the Python implementations of the same method, see the separate `pyHEBCPF`
v5.2 suite.

## Numerical evidence

### V5.2 pure MATLAB versus MEX

Both V5.2 implementations were compared with queue `parfeval`, the `bandit`
scheduler, a shared persistent 23-worker pool, and three repeats per case.
Pool startup was excluded. The run used MATLAB R2022a on the recorded Windows
x64 comparison host (Intel Core i9-13900K, 24 cores/32 threads, 128 GB RAM).

Across the 20 cases, both packages returned the same case-by-case solution
counts, totaling 3,254 solutions under the corrected V5.2 preprocessing. Based
on the rounded published table, pure MATLAB required 1,145.3 s in summed search
wall time and MEX required 586.3 s, an aggregate factor of 1.95x. The
geometric-mean case-wise wall factor is 1.59x, while the geometric-mean
per-trace kernel factor is 1.97x. Small cases can be dominated by dispatch and
bookkeeping overhead; the compiled advantage is clearest on the larger cases.

The full case table, scope, and limitations are in
`MATLAB_BENCHMARK_V5.2.md`. The surviving artifact does not contain
per-solution residual vectors or the raw repeated-run records, so V5.2 makes a
solution-count agreement claim but no new residual-distance claim from this
comparison.

### Historical v4 and V5 scheduler evidence

The retained v4 benchmark is the direct historical baseline. It covered all
20 bundled cases through case57 on MATLAB R2022a/Windows x64 with 23 workers:

| v4 implementation | Cases | Aggregate wall time | Trace wall time | Solutions summed across cases |
| --- | ---: | ---: | ---: | ---: |
| MEX v4 | 20 | 718.656 s | 716.062 s | 3256 |
| Pure MATLAB v4 | 20 | 1272.966 s | 1271.222 s | 3256 |

The v4 MEX and pure-MATLAB packages agreed on every solution count; their
largest order-independent nearest-solution distance was `3.154e-7`, below the
shared `4e-7` deduplication tolerance. These v4 timings are historical and are
not mixed with the V5 scheduler timings.

The V5 scheduler benchmark compares `scan`, `bandit`, and `novelty` on the MEX
implementation for every bundled case through case57, with three repeats per
policy. Its report and CSV are
`SCHEDULER_BENCHMARK_v5.md` and `scheduler_benchmark_v5_summary.csv`.

All published V5 numerical and timing values come solely from official
dataset `20260724_3x23`, produced on the designated main job machine with
MATLAB R2022a on Windows x64 and 23 workers. Results produced on other machines
are not merged into the release evidence. `BENCHMARK_METADATA_v5.json` records
the platform, policy settings, deduplication settings, and provenance; the
exact Git commit and machine identifier were not recorded at execution time and
are therefore reported as unknown rather than inferred afterward. The run used
the source-equivalent V5.5 development tree; metadata records both the
benchmark-time and reproducibly rebuilt release KLU hashes.

All 180 measured runs passed. Summed exhaustive wall time was effectively tied
at 606.947 s (`scan`), 607.184 s (`bandit`), and 606.698 s (`novelty`). In
contrast, the summed trace-to-90% metric fell from 35,975 for `scan` to 16,611
for `bandit` (53.8% lower) and 11,975 for `novelty` (66.7% lower). A separate
40-run order-independent audit matched every adaptive-policy solution set to
its scan reference; the worst distance was `2.66e-9`, below `4e-7`.

The report and PDFs include exhaustive-time, trace-to-90%, anytime-discovery,
and case14mod/case39 connectivity plots. Connectivity panels contain identical
validated solution sets; their colors and edges visualize scheduler traversal
history.

## Checkpoint and resume

Long parallel runs write `temp_result.mat`. Resume in the same implementation
folder and with the same case/load preprocessing:

```matlab
cd('<path_to_HEBCPF>/HEBCPF_MEX_v5.2')
addpath(pwd)
load temp_result.mat

global HEBCPOLICY
HEBCPOLICY = 'novelty';                 % may differ from the saved policy
wait = 0;
run('runVBook_hybrid_parfeval.m')
```

The three canonical schedulers share checkpoint state. Learned equation gains
are reused by `bandit` and `novelty`; policy-local cursors and novelty
projections are rebuilt. The bundled
`test_checkpoint_policy_compatibility_v5.m` harness exercises cross-policy
resume in both implementations; no result produced on a non-designated machine
is included in the published release evidence.

Two environment controls exist for deterministic checkpoint testing:
`HEBCPF_CHECKPOINT_TRACE_INTERVAL` and `HEBCPF_STOP_AFTER_CHECKPOINT=true`.
When unset, the production time-based cadence is unchanged.

## Documentation and outputs

- `HEBCPF_Suite_Overview.pdf`: release-level V5-versus-v4 comparison.
- `RELEASE_NOTES_V5.2.md`: GitHub Release description for tag `v5.2.0`.
- `RELEASE_CHECKLIST_V5.2.md`: final repository, documentation, test, asset,
  and Git verification checklist.
- `HEBCPF_MEX_v5.2/HEBCPF_User_Guide.pdf`: MEX user guide.
- `HEBCPF_matlab_v5.2/HEBCPF_User_Guide.pdf`: pure-MATLAB user guide.
- `MATLAB_BENCHMARK_V5.2.md`: GitHub-readable pure-vs-MEX protocol, complete
  case table, aggregate statistics, and limitations.
- `MATLAB_BENCHMARK_METADATA_V5.2.json`: machine-readable identity and
  provenance boundary for that comparison.
- `matlab_benchmark.tex/png`: pure-vs-MEX benchmark fragment.
- `POLICY_README.md`: scheduler algorithms and complexity.
- `SCHEDULER_BENCHMARK_v5.md`: official V5 three-policy report.
- `scheduler_benchmark_v5_summary.csv`: machine-readable averaged results.
- `scheduler_benchmark_v5_overall.csv`: suite-level scheduler statistics.
- `BENCHMARK_METADATA_v5.json`: official dataset identity and provenance.
- `release_assets/HEBCPF-v5-benchmark-raw.zip`: portable raw CSV, trace
  histories, validation data, and references, attached separately to the
  GitHub Release with its SHA-256 checksum.
- `release_assets/HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip`: exact pinned
  source, licenses, build metadata, smoke check, and V5 `klurf.mexw64`, attached
  beside the standalone V5 MEX and both SHA-256 files.

To produce a future official benchmark on the designated main job machine,
set `HEBCPF_SCHED_BENCH_WORKERS=23` and `HEBCPF_SCHED_BENCH_ROOT` before running
`run_scheduler_benchmark_v5`. After reviewing the completed run and marking
its generated metadata as `official release evidence`, set
`HEBCPF_PUBLISH_BENCHMARK=true` before running
`generate_scheduler_benchmark_report_v5`. The report generator rejects
unapproved metadata, incomplete case-policy matrices, failed validation rows,
and worker-count mismatches.

When rebuilding the LaTeX documentation, run `pdflatex` from the document's
own folder at least twice (a third pass is harmless). The first pass creates
the table-of-contents and cross-reference data; the next pass writes those
entries into the PDF. In TeXstudio, use a build command configured to rerun
LaTeX when references have changed.

`run_merged_case` returns a result struct and solution matrix. In V5 the
result also includes the normalized scheduler name and `trace_metrics`, with
traces to 50%, 90%, and 99% of the final solution count, wall time to 90%,
worker-trace statistics, scheduler selection time, and trace histories.

## License, attribution, and citation

HEBCPF is distributed under the BSD 3-Clause License; see `LICENSE`.
MATPOWER and SuiteSparse/KLU remain third-party software under their own terms;
see `NOTICE` and the individual source notices. Citation metadata is provided
in `CITATION.cff`.
