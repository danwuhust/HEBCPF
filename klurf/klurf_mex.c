/* klurf_mex.c  --  KLU "symbolic-once" solve for repeated same-pattern systems.
 *
 *   X = klurf(A, B)     solve A*X = B.  A: real sparse n x n, B: real dense n x nrhs.
 *   klurf('free')       free the cached factorization.
 *
 * Persistent (per-MATLAB-process) Symbolic + Numeric + pattern cache. On each call
 * the incoming pattern (n, nnz, Ap, Ai) is compared to the cache:
 *   - match  -> klu_l_refactor : reuse fill-reducing ordering + pivot sequence,
 *               recompute numeric values only (exact Jacobian, quadratic Newton kept).
 *   - differ -> klu_l_analyze + klu_l_factor : full factorization, refresh cache.
 * If refactor reports a zero-pivot / singular breakdown, it falls back to a full
 * klu_l_factor for that call (safety). Uses the int64 (_l) KLU API to match MATLAB
 * mwIndex on 64-bit; real-valued only.
 *
 * SPDX-License-Identifier: LGPL-2.1+   (uses SuiteSparse/KLU, T. Davis)
 */
#include "mex.h"
#include "klu.h"
#include "SuiteSparse_config.h"
#include <string.h>
#include <stdlib.h>

static klu_l_symbolic *Symbolic = NULL;
static klu_l_numeric  *Numeric  = NULL;
static klu_l_common    Common;
static int      common_init = 0;
static int      atexit_set  = 0;
static int64_t  cached_n = -1, cached_nz = -1;
static int64_t *cached_Ap = NULL;   /* size n+1 */
static int64_t *cached_Ai = NULL;   /* size nz  */

static void free_all(void)
{
    if (Numeric)   { klu_l_free_numeric(&Numeric, &Common);   Numeric = NULL; }
    if (Symbolic)  { klu_l_free_symbolic(&Symbolic, &Common); Symbolic = NULL; }
    if (cached_Ap) { free(cached_Ap); cached_Ap = NULL; }
    if (cached_Ai) { free(cached_Ai); cached_Ai = NULL; }
    cached_n = -1; cached_nz = -1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1) mexErrMsgTxt("usage: X = klurf(A,B)  or  klurf('free')");
    if (mxIsChar(prhs[0])) { free_all(); return; }        /* klurf('free') */
    if (nrhs != 2) mexErrMsgTxt("usage: X = klurf(A,B)");

    const mxArray *A = prhs[0], *B = prhs[1];
    if (!mxIsSparse(A) || mxIsComplex(A)) mexErrMsgTxt("A must be real sparse");
    if ( mxIsSparse(B) || mxIsComplex(B)) mexErrMsgTxt("B must be real dense");

    mwSize n = mxGetM(A);
    if (n != mxGetN(A)) mexErrMsgTxt("A must be square");
    mwIndex *Ap_mw = mxGetJc(A);
    mwIndex *Ai_mw = mxGetIr(A);
    double  *Ax    = mxGetPr(A);
    mwSize   nz    = Ap_mw[n];
    mwSize   bn    = mxGetN(B);
    if (mxGetM(B) != n) mexErrMsgTxt("size(B,1) must equal size(A,1)");

    if (!common_init) {
        /* CRITICAL: in a MEX build SuiteSparse defaults to mxMalloc/mxFree, which are
         * auto-freed when mexFunction returns -> our persistent Symbolic/Numeric would
         * dangle on the next call. Force plain C malloc/free so they survive. */
        SuiteSparse_config_malloc_func_set(malloc);
        SuiteSparse_config_calloc_func_set(calloc);
        SuiteSparse_config_realloc_func_set(realloc);
        SuiteSparse_config_free_func_set(free);
        klu_l_defaults(&Common);
        common_init = 1;
    }
    if (!atexit_set)  { mexAtExit(free_all); atexit_set = 1; }

    /* int64 copies of the CSC pattern (mwIndex is unsigned; klu_l wants int64_t*) */
    int64_t *Ap = (int64_t*) mxMalloc((n+1)*sizeof(int64_t));
    int64_t *Ai = (int64_t*) mxMalloc((nz>0?nz:1)*sizeof(int64_t));
    { mwSize j; for (j=0;j<=n;j++) Ap[j]=(int64_t)Ap_mw[j]; }
    { mwSize k; for (k=0;k<nz;k++) Ai[k]=(int64_t)Ai_mw[k]; }

    int same = (Symbolic && Numeric && cached_Ap && cached_Ai
                && cached_n==(int64_t)n && cached_nz==(int64_t)nz);
    if (same) {
        if (memcmp(cached_Ap, Ap, (n+1)*sizeof(int64_t)) != 0) same = 0;
        else if (nz>0 && memcmp(cached_Ai, Ai, nz*sizeof(int64_t)) != 0) same = 0;
    }

    if (same) {
        int ok = klu_l_refactor(Ap, Ai, Ax, Symbolic, Numeric, &Common);
        if (!ok || Common.status != KLU_OK) {
            /* refactor breakdown: redo a full numeric factor (reuse Symbolic) */
            klu_l_free_numeric(&Numeric, &Common);
            Numeric = klu_l_factor(Ap, Ai, Ax, Symbolic, &Common);
        }
    }
    if (!same || !Numeric) {
        if (Numeric)  klu_l_free_numeric(&Numeric, &Common);
        if (Symbolic) klu_l_free_symbolic(&Symbolic, &Common);
        Symbolic = klu_l_analyze((int64_t)n, Ap, Ai, &Common);
        if (!Symbolic) mexErrMsgTxt("klu_l_analyze failed");
        Numeric = klu_l_factor(Ap, Ai, Ax, Symbolic, &Common);
        if (!Numeric)  mexErrMsgTxt("klu_l_factor failed");
        if (cached_Ap) free(cached_Ap);
        if (cached_Ai) free(cached_Ai);
        cached_Ap = (int64_t*) malloc((n+1)*sizeof(int64_t));
        cached_Ai = (int64_t*) malloc((nz>0?nz:1)*sizeof(int64_t));
        memcpy(cached_Ap, Ap, (n+1)*sizeof(int64_t));
        if (nz>0) memcpy(cached_Ai, Ai, nz*sizeof(int64_t));
        cached_n = (int64_t)n; cached_nz = (int64_t)nz;
    }

    /* solve in place on a copy of B */
    plhs[0] = mxCreateDoubleMatrix(n, bn, mxREAL);
    double *X = mxGetPr(plhs[0]);
    memcpy(X, mxGetPr(B), (size_t)n*bn*sizeof(double));
    if (!klu_l_solve(Symbolic, Numeric, (int64_t)n, (int64_t)bn, X, &Common))
        mexErrMsgTxt("klu_l_solve failed");

    mxFree(Ap); mxFree(Ai);
}
