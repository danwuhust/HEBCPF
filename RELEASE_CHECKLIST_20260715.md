# HEBCPF 2026.07.15 Release Checklist

Use this checklist before uploading the package as a GitHub release.

## Include

- `README.md` - top-level suite overview, quick start, and benchmark summary.
- `HEBCPF_Suite_Overview.tex` and `HEBCPF_Suite_Overview.pdf` - release-level solver-suite overview.
- `RELEASE_NOTES_20260715.md` - GitHub release-note text.
- `CHANGELOG.md`, `CITATION.cff`, `LICENSE`, and `NOTICE`.
- `HEBCPF_MEX_v4_20260715/` - Windows x64 MEX solver.
- `HEBCPF_matlab_v4_20260715/` - pure MATLAB solver.
- Solver-level README and `HEBCPF_User_Guide.tex/.pdf` files in each solver folder.
- Bundled case files and required support files already present in each solver folder.
- `run_benchmark_up_to_57bus_20260715.m` - reproducible benchmark driver.
- `benchmark_results_20260715/20260715_194046/BENCHMARK_SUMMARY.md` and
  `benchmark_summary_all.csv` - completed release benchmark report and data.

## Exclude Generated Run Artifacts

The `.gitignore` file excludes normal local artifacts:

- `temp_result.mat` checkpoint files.
- Manual checkpoint files such as `my_checkpoint.mat` and `*_checkpoint.mat`.
- Generated solver summaries such as `results_parfeval_summary.csv`.
- Benchmark logs, per-case snapshots, `benchmark_records.mat`, and incomplete benchmark attempts.
- LaTeX intermediates such as `.aux`, `.log`, `.out`, and `.toc`.

Do not delete `examples.mat`; it is a required MEX rebuild input and is intentionally not ignored.

## Resume Workflow To Verify Before Release

For either v4 solver, resume an interrupted parallel run from the solver folder with:

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

## Final Local Checks Used For This Package

- Both solver packages completed all 20 benchmark cases through 57 buses using queue `parfeval`.
- Solution counts matched case-by-case; order-independent solution-set checks remained below `4e-7`.
- All three PDFs were rebuilt from the updated LaTeX sources.
- Documentation was scanned for stale release claims, obsolete checkpoint naming, and old slope-limit defaults.
- MATLAB `checkcode` was run on key v4 drivers; remaining messages are non-blocking style/preallocation notices.
