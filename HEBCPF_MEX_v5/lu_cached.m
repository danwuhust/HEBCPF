function [Lf,Uf,Ppm,Qqm] = lu_cached(H)
%LU_CACHED Sparse LU of H with a cached column ordering (suggestion B).
%
% H's sparsity pattern is fixed for a given power system -- every holo step
% of every trace factors a matrix with the same structure, only the values
% (P0/Q0/V0/W0 diagonals) change. The stock hot path
%     [Lh,Uh,Ph,Qh] = lu(H)
% re-runs UMFPACK's column analysis on every call (measured 17-20% of
% worker self-time). Here the column permutation cv is computed ONCE (first
% call, full 4-output lu), then subsequent calls use the 3-output
%     [L,U,p] = lu(H(:,cv),'vector')
% which factors in the given column order with partial pivoting only -- no
% analysis stage (measured 4-4.9x faster on case14mod/39/57mod, factorization
% exact to machine precision).
%
% Returned Ppm/Qqm satisfy the same identity the callers rely on:
%     H(Ppm, cv) = Lf*Uf   with   Qqm = inverse permutation of cv,
% exactly matching the stock  [Ppm,~]=find(Ph'); [Qqm,~]=find(Qh')  vectors.
%
% Cache key is (n, nnz): O(1) per call, same policy as MakeJacobianD. A
% same-nnz pattern change (theoretically possible if a V0/W0 component is
% exactly zero) only degrades fill/pivoting -- GPLU with partial pivoting
% stays CORRECT for any column order, so results are never wrong.

persistent cv Qqm_c n_c nnz_c
n = size(H,1);
if isempty(cv) || n ~= n_c || nnz(H) ~= nnz_c
    [Lh,Uh,Ph,Qh] = lu(H);
    [Ppm,~,~] = find(Ph');
    [Qqm,~,~] = find(Qh');
    [cvv,~,~] = find(Qh);            % H*Qh = H(:,cvv)
    cv = cvv; Qqm_c = Qqm; n_c = n; nnz_c = nnz(H);
    Lf = full(Lh); Uf = full(Uh);
    return;
end
[Lh,Uh,p] = lu(H(:,cv), 'vector');   % H(p,cv) = Lh*Uh, no column analysis
Ppm = p(:);
Qqm = Qqm_c;
Lf = full(Lh);
Uf = full(Uh);
end
