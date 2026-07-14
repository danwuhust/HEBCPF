# klurf — KLU "symbolic-once" refactor MEX for the HEBCPF corrector

`klurf` accelerates the Phase-I bordered Newton solve inside `branch_trace`. It exposes
KLU's `klu_analyze` (fill-reducing ordering + symbolic factorization, done **once**) and
`klu_refactor` (reuse that ordering, recompute only the numeric values) for the repeated,
same-sparsity-pattern linear systems that arise across Newton iterations and continuation
curves. This keeps the **exact** Jacobian every iteration (quadratic convergence, exact
solution counts preserved) while skipping the AMD ordering + symbolic analysis + pivot
search that a full factorization repeats on every solve.

Measured on the numerical test bed (Intel i9-13900K, Win11, MATLAB R2022a), serial tracing,
solution counts identical to the dense solve on every case:

| system | Nm | refactor speedup vs dense |
|--------|----|---------------------------|
| ≤ case14mod | ≤27 | ~1.0× (overhead-bound, no gain) |
| case30 / case33bw | 59 / 65 | 1.12× / 1.25× |
| case39 | 77 | 1.43× |
| case118 | 235 | ~7.5× |

The advantage scales with the system size Nm; it is essentially free (auto-disabled) on the
small cases. **No gating is applied** — the corrector auto-uses `klurf` whenever the MEX is
present on the MATLAB path and falls back to the dense solve otherwise.

## Usage (from the solver)

Nothing to do: each solver folder that ships `klurf.mexw64` uses it automatically. To force a
specific corrector mode for benchmarking, set the global before running:

```matlab
global USEKLU
USEKLU = 0;   % dense bordered solve  (J1\F1)
USEKLU = 1;   % KLU full factorization each iteration  (needs klu.mexw64, not shipped)
USEKLU = 2;   % KLU symbolic-once refactor  (klurf)   <-- default when klurf is present
USEKLU = [];  % auto: refactor if klurf is present, else dense (default)
```

Direct API: `X = klurf(A, B)` solves `A*X = B` (A real sparse n×n, B real dense n×nrhs).
`klurf('free')` releases the cached factorization.

## Prebuilt binary

`klurf.mexw64` (Windows x64) is included and already copied into the `HEBCPF_MEX_*`
folders. MATLAB v2/v3 can also use a compatible binary placed in their release folder.
MATLAB v4 deliberately keeps its correctors pure MATLAB and does not call `klurf`.

## Rebuilding (Linux / macOS, or from source)

`klurf` links against SuiteSparse/KLU (open source, LGPL-2.1+; Tim Davis). To rebuild:

1. Get SuiteSparse (https://github.com/DrTimothyAldenDavis/SuiteSparse).
2. Copy `klurf_mex.c` and `klurf_make.m` into `SuiteSparse/KLU/MATLAB/`.
3. From that directory in MATLAB: `klurf_make`  (real-valued, no CHOLMOD; uses AMD/COLAMD/BTF).
4. Copy the resulting `klurf.<mexext>` next to the solver (or into the folder you run from).

`klurf_make.m` compiles the needed SuiteSparse sources directly, so no prior SuiteSparse
build/library is required — only a working `mex` compiler.

## Note on correctness

`klu_refactor` reuses the first factorization's pivot ordering. A fixed-order LU is still an
**exact** factorization of the current matrix, so this does not perturb the solution branches:
validated to reproduce the dense-solve solution counts exactly on case9, case14mod, case33bw,
**case39 (ill-conditioned, slopes ~3e9)**, and case30. If a refactor ever hits a zero pivot,
`klurf` transparently falls back to a full factorization for that solve.
