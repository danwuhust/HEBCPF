# HEBCPF v5.2.0 Release Checklist

Release identity:

- Version: `5.2.0`
- Git tag: `v5.2.0`
- Release date: `2026-07-30`
- Repository: `https://github.com/danwuhust/HEBCPF`

## Repository contents

- [ ] `README.md`, `CHANGELOG.md`, `CITATION.cff`, `LICENSE`, and `NOTICE` are
      staged.
- [ ] `RELEASE_NOTES_V5.2.md` matches the GitHub release text.
- [ ] `HEBCPF_MEX_v5.2/` contains the Windows x64 solver, required `.mexw64`
      files, cases, README, and `HEBCPF_User_Guide.tex/.pdf`.
- [ ] `HEBCPF_matlab_v5.2/` contains the pure-MATLAB solver, cases, README, and
      `HEBCPF_User_Guide.tex/.pdf`.
- [ ] `HEBCPF_Suite_Overview.tex/.pdf` describes only the MATLAB/MEX suite and
      points Python users to the separate `pyHEBCPF` repository.
- [ ] `MATLAB_BENCHMARK_V5.2.md`,
      `MATLAB_BENCHMARK_METADATA_V5.2.json`, and `matlab_benchmark.tex/.png`
      are staged together.
- [ ] V5 scheduler reports, CSV summaries, plots, metadata, and validation
      scripts remain present as historical/official V5 scheduler evidence.
- [ ] `POLICY_README.md`, scheduler guide fragments, checkpoint harnesses, and
      smoke tests are staged.

## Exclusions

- [ ] No `*.asv`, LaTeX intermediates, checkpoints, generated solver summary
      CSV files, incomplete benchmark directories, or absolute local paths are
      staged.
- [x] `examples.mat` remains in the MEX folder; it is a required rebuild input.
- [ ] The ignored `scheduler_benchmark_v5/` raw working tree is not committed;
      its portable reviewed archive is supplied as a release asset.
- [ ] Python solver folders and Python dependency text are not part of this
      MATLAB-only release.

## Numerical and documentation checks

- [x] Confirm both packages use `slope_max = 4e5` and the same scheduler and
      checkpoint defaults.
- [x] Confirm `quadr_matrix` returns only `Ma`, `ra`, and `I` and all callers
      use that interface.
- [x] Run MATLAB `checkcode` on the changed entry points and preprocessing
      functions in both packages.
- [x] Run `smoke_test_v5` with MATPOWER available. Record pass/fail only; do
      not promote timing from an unapproved machine.
- [ ] Run `test_checkpoint_policy_compatibility_v5` when a full checkpoint
      validation is desired.
- [x] Rebuild all three PDFs from their `.tex` sources and inspect every page.
- [x] Validate `CITATION.cff` and both benchmark JSON files.
- [x] Verify the pure/MEX table contains 20 cases, 3,254 total solutions, and
      the documented rounded aggregate statistics.

## Release assets

Upload the six files described by `release_assets/README.md`:

- [ ] `HEBCPF-v5-klurf-SuiteSparse-v7.12.2.mexw64`
- [ ] its `.sha256` file
- [ ] `HEBCPF-v5-klurf-source-SuiteSparse-v7.12.2.zip`
- [ ] its `.sha256` file
- [ ] `HEBCPF-v5-benchmark-raw.zip`
- [ ] its `.sha256` file

The KLU corresponding-source asset must accompany any GitHub release that
redistributes the KLU-linked binary. Verify all three SHA-256 values before
publishing.

## Git verification

- [ ] Review `git status --short`, `git diff --check`, and
      `git diff --cached --stat`.
- [x] Confirm no file exceeds GitHub's normal single-file limit.
- [ ] Commit the reviewed tree, pull/rebase `origin/main`, and push `main`.
- [ ] Create annotated tag `v5.2.0` on the reviewed commit and push it.
- [ ] Create the GitHub Release from `RELEASE_NOTES_V5.2.md`, attach the six
      assets, and verify the three PDFs in the GitHub web interface.
