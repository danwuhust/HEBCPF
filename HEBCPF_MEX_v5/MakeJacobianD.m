function J = MakeJacobianD(x,S,n,action)
% D variant: skip the triplet decode AND the per-call sparse rebuild entirely.
% Vt = sparse(S(:,2),S(:,1),S(:,3), n*n, n) encodes all stacked constraint
% matrices in transposed/blocked form and is built once (persistent).
% Then  J = 2*reshape(Vt*x, n, n)'  reproduces exactly
%   J(c,j) = 2*(M_c' x)_j  =  2*sparse(SJ(:,2),SJ(:,1),SJ(:,3),n,n)'
% from the scalar version. Per call only a sparse mat-vec + reshape remain.
%
% action = 'initialize' : build/refresh the cache only (deep isequal check), return [].
% action = 'sparse'     : return J as a SPARSE Nm x Nm matrix directly (reduced operator).
% action = 'residual'   : return the Jacobian-free residual g = J(x)*x/2 = [x'M_i x]_i.
%
% Hot path (non-initialize) checks DIMENSIONS ONLY -- never isequal(S,cachedS) -- so the
% per-call overhead is O(1). (An earlier sparse-mode rewrite regressed this by adding a
% per-call isequal, which is O(nnz(S)) and cost ~4-5 us/call on case118.)
persistent Vt cachedS nc Nc Vtr rowsJ colsJ
N = size(S,1);
act = ''; if nargin>=4; act = action; end

if strcmp(act,'initialize')
    if isempty(cachedS) || n~=nc || N~=Nc || ~isequal(S,cachedS)
        [Vt,Vtr,rowsJ,colsJ] = build_ops(S,n);
        cachedS = S; nc = n; Nc = N;
    end
    J = [];
    return;
end
if isempty(nc) || n ~= nc || N ~= Nc            % hot path: dimensions only
    [Vt,Vtr,rowsJ,colsJ] = build_ops(S,n);
    cachedS = S; nc = n; Nc = N;
end

if strcmp(act,'sparse')
    vr = Vtr * x;                               % |nz| x 1  (small mat-vec)
    J  = sparse(rowsJ, colsJ, 2*vr, n, n);
elseif strcmp(act,'residual')
    % Jacobian-FREE constraint residual  g = J(x)*x/2 = [x'M_i x]_i, via the reduced
    % operator (no dense reshape or matrix assembly), retained for cache-aware callers.
    vr = Vtr * x;
    J  = accumarray(rowsJ, vr .* x(colsJ), [n 1]);
else
    J = 2 * reshape(Vt*x, n, n)';               % dense (default, unchanged)
end
end

function [Vt,Vtr,rowsJ,colsJ] = build_ops(S,n)
Vt    = sparse(S(:,2), S(:,1), S(:,3), n*n, n);
nz    = find(any(Vt,2));                         % nonzero rows of Vt = structural nnz of J
rowsJ = floor((nz-1)/n) + 1;
colsJ = nz - (rowsJ-1)*n;
Vtr   = Vt(nz,:);                                % reduced operator: |nz| x n
end
