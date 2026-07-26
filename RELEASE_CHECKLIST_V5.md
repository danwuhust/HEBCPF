# HEBCPF V5 Package Checklist

Release identity:

- Version: `5.0.0`
- Git tag: `v5.0.0`
- Release date: 2026-07-26

Use this checklist before uploading the package as a GitHub release.

## Include

- `README.md` - top-level suite overview, quick start, and benchmark summary.
- `HEBCPF_Suite_Overview.tex` and `HEBCPF_Suite_Overview.pdf` - release-level solver-suite overview.
- `RELEASE_NOTES_V5.md` - GitHub release-note text.
- `CHANGELOG.md`, `CITATION.cff`, `LICENSE`, and `NOTICE`.
- `HEBCPF_MEX_v5/` - Windows x64 MEX solver.
- `HEBCPF_matlab_v5/` - pure MATLAB solver.
- Solver-level README and `HEBCPF_User_Guide.tex/.pdf` files in each solver folder.
- Bundled case files and required support files already present in each solver folder.
- `run_scheduler_benchmark_v5.m` - reproducible three-policy benchmark driver.
- `generate_scheduler_benchmark_report_v5.m` and
  `validate_scheduler_solution_sets_v5.m` - guarded report and validation tools.
- `scheduler_benchmark_v5_summary.csv` and `scheduler_benchmark_v5_overall.csv`
  - published V5 scheduler results.
- `BENCHMARK_METADATA_v5.json` - official dataset identity, platform,
  scheduler parameters, and provenance.
- `release_assets/README.md` - explains the separately uploaded raw-data assets.
- `release_assets/HEBCPF-v5-benchmark-raw.zip` and matching `.sha256` file
  - upload both manually as GitHub Release assets; `.gitignore` intentionally
  keeps these generated files out of the repository commit.
- `release_assets/HEBCPF-v5-klurf-SuiteSparse-v7.12.2.mexw64` and
  matching `.sha256` file - upload together as GitHub Release assets.
- `release_assets/HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip` and matching
  `.sha256` file - upload beside the V5 MEX as its pinned corresponding-source
  and license package.
- `connectivity_case30_v4_reference.{fig,png}` and
  `connectivity_case57_v4_reference.{fig,png}` - retained byte-for-byte as the
  released-v4 connectivity references. Their historical raw v4 snapshots are
  not distributed in V5; `generate_reference_connectivity_v5.m` is an optional
  provenance helper for authors who already hold those snapshots.

## Exclude Generated Run Artifacts

The `.gitignore` file excludes normal local artifacts:

- `temp_result.mat` checkpoint files.
- Manual checkpoint files such as `my_checkpoint.mat` and `*_checkpoint.mat`.
- Generated solver summaries such as `results_parfeval_summary.csv`.
- Benchmark logs, per-case snapshots, `benchmark_records.mat`, and incomplete benchmark attempts.
- LaTeX intermediates such as `.aux`, `.log`, `.out`, and `.toc`.

Do not delete `examples.mat`; it is a required MEX rebuild input and is intentionally not ignored.

## Resume Workflow To Verify Before Release

For either V5 solver, resume an interrupted parallel run from the solver folder with:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval.m')          % queue parfeval
% run('runVBook_hybrid_parallel.m')        % parfor
```

For the MEX solver, the retained barrier scheduler can also be resumed with:

```matlab
load temp_result.mat
wait = 0;
run('runVBook_hybrid_parfeval_barrier.m')
```

## Official Release Evidence

Dataset `20260724_3x23` from the designated 23-worker main job machine supports
the published V5 numerical and timing claims. Its metadata records the
source-equivalent V5.5 execution tree and both KLU hashes. Do not merge results
from another computer into the tracked evidence.

- The retained v4 baseline completed all 20 benchmark cases through 57 buses using queue `parfeval`.
- The V5 timing benchmark completed 180/180 runs (20 cases, three policies,
  three repeats) and the order-independent solution audit passed 40/40 checks.
- Checkpoint and smoke-test harnesses are included for platform-specific use;
  do not publish their outcomes unless they are rerun and recorded on the
  designated main job machine.
- Solution counts matched case-by-case; order-independent solution-set checks remained below `4e-7`.
- All three PDFs were rebuilt from the updated LaTeX sources.
- Documentation was scanned for stale release claims, obsolete checkpoint naming, and old slope-limit defaults.
- MATLAB `checkcode` was run on key drivers; remaining messages are non-blocking style/preallocation notices.

## Report Publication Guard

- Confirm that completed benchmark metadata says
  `"publication_status": "official release evidence"` only after review.
- Confirm `worker_count` is 23 and the metadata lists all 20 cases, all three
  policies, and three repeats.
- Set `HEBCPF_PUBLISH_BENCHMARK=true` only for the reviewed official dataset.
- Run `generate_scheduler_benchmark_report_v5` only for a newly completed and
  approved V5 dataset; it must reject failed rows,
  incomplete matrices, unapproved metadata, and worker-count mismatches.
- Confirm that no absolute local paths occur in tracked documentation or in the
  raw-data release archive.
