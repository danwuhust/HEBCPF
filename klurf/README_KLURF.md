# klurf - KLU symbolic-once refactor MEX for the HEBCPF corrector

`klurf` accelerates the Phase-I bordered Newton solve inside `branch_trace`. It
exposes KLU's `klu_analyze` operation once for a fixed sparsity pattern, then
uses `klu_refactor` to reuse that ordering while recomputing numerical values
for later Newton iterations and continuation curves.

This retains the current Jacobian at every iteration while avoiding repeated
ordering and symbolic-analysis work. The solver automatically falls back to a
full factorization if refactorization encounters a zero pivot.

## Usage from the solver

No setup is needed when `klurf.mexw64` is beside the MEX solver. To select a
corrector mode explicitly for diagnostics or controlled comparisons, set:

```matlab
global USEKLU
USEKLU = 0;   % MATLAB dense bordered solve (J1\F1)
USEKLU = 1;   % KLU full factorization each iteration (requires klu.mexw64)
USEKLU = 2;   % KLU symbolic-once refactor through klurf
USEKLU = [];  % auto: use klurf when available, otherwise MATLAB solve
```

Direct API: `X = klurf(A, B)` solves `A*X = B`, where `A` is a real sparse
square matrix and `B` is a real dense right-hand side. `klurf('free')`
releases the cached factorization.

## Prebuilt binary

`klurf.mexw64` is provided for Windows x64 in this folder and in
`HEBCPF_MEX_v5.2/`. The pure-MATLAB package deliberately does not call it.

The binary was rebuilt from official SuiteSparse `v7.12.2`, commit
`42151688813c45846a597edcb601435a0e38f3dd`. The corresponding GitHub Release
asset contains the exact required SuiteSparse component source, license texts,
build metadata, wrapper source, build helper, deterministic smoke check, and
binary. Its SHA-256 matches the two repository copies. See
`release_assets/README.md` and the archived `BUILD_INFO.md`.

## Rebuilding

1. Obtain the pinned SuiteSparse `v7.12.2` source or the corresponding-source
   archive attached to the HEBCPF GitHub Release.
2. Copy `klurf_mex.c` and `klurf_make.m` into `SuiteSparse/KLU/MATLAB/`.
3. Run `klurf_make` in that directory from MATLAB with a configured MEX C
   compiler.
4. Copy the resulting `klurf.<mexext>` beside the solver.

`klurf_make.m` compiles the required SuiteSparse sources directly, so no
prebuilt SuiteSparse library is required.

## Correctness scope

Reusing a symbolic ordering does not freeze the numerical Jacobian: numerical
values are refactorized for each solve. Historical V5 validation compared the
KLU and MATLAB corrector paths on the bundled cases. For V5.2, the published
20-case pure-MATLAB-versus-MEX comparison reports matching completed solution
counts. Its raw residual vectors are unavailable, so V5.2 documentation does
not make a new residual-distance claim from that table.
