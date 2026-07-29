function [Lf,Uf,Ppm,Qqm] = lu_cached(H)
%LU_CACHED Sparse LU with a cached column ordering.
% The factorization remains a MATLAB built-in LU; no compiled solver is used.

persistent cv Qqm_c n_c nnz_c
n = size(H,1);
if isempty(cv) || n ~= n_c || nnz(H) ~= nnz_c
    [Lh,Uh,Ph,Qh] = lu(H);
    [Ppm,~,~] = find(Ph');
    [Qqm,~,~] = find(Qh');
    [cvv,~,~] = find(Qh);
    cv = cvv; Qqm_c = Qqm; n_c = n; nnz_c = nnz(H);
    Lf = full(Lh); Uf = full(Uh);
    return;
end
[Lh,Uh,p] = lu(H(:,cv), 'vector');
Ppm = p(:);
Qqm = Qqm_c;
Lf = full(Lh);
Uf = full(Uh);
end
