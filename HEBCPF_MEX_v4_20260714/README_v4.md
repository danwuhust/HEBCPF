# HEBCPF v4_mex

Based on `HEBCPF_MEX_v3.20260712` with four validated performance changes
(compound "ABCE"). Same inputs, same workflow (`main.m`, `run_batch*.m`,
`run_merged_case.m`), and — verified by order-insensitive pairwise matching
on all bundled cases ≤ 57 buses — the identical solution sets.

## Changes vs v3

| Id | Change | Files | Measured effect (23 workers) |
|----|--------|-------|------------------------------|
| A | Parfeval driver: row barrier → global work queue with per-equation serialization (zero redundant traces). Old driver kept as `runVBook_hybrid_parfeval_barrier.m`. | `runVBook_hybrid_parfeval.m` | 1.10–1.33× parfeval wall; shrinks toward parity at ≤8-worker pools |
| B | Sparse LU of the holomorphic matrix H: UMFPACK column analysis cached per case, then 3-output GPLU in the cached column order. | `lu_cached.m`, `holomorphic_hybrid_5_with_mex.m` | LU step 4.1–4.9×; serial 1.08–1.21× |
| C | `holomorphic_para_sp`: P0/Q0 cell loop → persistent stacked-operator mat-vec; spdiags → direct `sparse(i,i,v)`; G/B/dl cached. | `holomorphic_para_sp.m` | serial +2–4% (case39) |
| E | `update_V_Pade`: per-call bus_n×(deg+1) power matrices → two BLAS mat-vecs. Bit-exact. | `update_V_Pade.m` | serial +13–14% (case39) |

A fifth candidate (D: klurf values-only refactor) won its micro-benchmark but
tied end-to-end and was **not** included.

## Reported development measurements (parfeval, 23 workers, vs v3 MEX, 2 reps)

These are development-time reference measurements. They were not re-executed during
the 2026.07.14 package-preparation pass.

| Case | v3 wall (s) | v4 wall (s) | speedup |
|------|---:|---:|---:|
| case9 | 1.61 | 1.33 | 1.21× |
| case14mod | 5.57 | 3.91 | 1.42× |
| case33bw | 8.22 | 5.66 | 1.45× |
| case39 | 49.93 | 38.64 | 1.29× |
| case30 | 101.81 | 63.96 | 1.59× |
| case57mod | 236.43 | 175.70 | 1.35× |

Serial mode benefits from B+C+E only (~1.3× on case39-class systems).

## Compatibility notes

- `run_merged_case(case,'parfeval')` and `run_batch_parfeval.m` use the new
  queue driver transparently; `'serial'` / `'parfor'` paths are unchanged
  apart from the B/C/E kernels.
- To A/B against the old scheduler on identical kernels, run
  `runVBook_hybrid_parfeval_barrier.m` directly.
- All caches (lu_cached, para_sp, MakeJacobianD, quad_form_vals) self-key on
  problem dimensions/values and refresh automatically when the case changes.
