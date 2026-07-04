function J = MakeJacobianD(x,S,n,action)
% D variant: skip the triplet decode AND the per-call sparse rebuild entirely.
% Vt = sparse(S(:,2),S(:,1),S(:,3), n*n, n) encodes all stacked constraint
% matrices in transposed/blocked form and is built once (persistent).
% Then  J = 2*reshape(Vt*x, n, n)'  reproduces exactly
%   J(c,j) = 2*(M_c' x)_j  =  2*sparse(SJ(:,2),SJ(:,1),SJ(:,3),n,n)'
% from the scalar version. Per call only a sparse mat-vec + reshape remain.
persistent Vt cachedS nc Nc
N = size(S,1);
initialize = nargin>=4 && strcmp(action,'initialize');
if initialize
    if isempty(cachedS) || n~=nc || N~=Nc || ~isequal(S,cachedS)
        Vt = sparse(S(:,2), S(:,1), S(:,3), n*n, n);
        cachedS = S;
        nc = n;
        Nc = N;
    end
    J = [];
    return;
end
if isempty(nc) || n ~= nc || N ~= Nc
    Vt = sparse(S(:,2), S(:,1), S(:,3), n*n, n);
    cachedS = S;
    nc = n;
    Nc = N;
end
J = 2 * reshape(Vt*x, n, n)';
end
