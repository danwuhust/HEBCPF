function q = quad_form_vals(x, Ma, n, action)
%QUAD_FORM_VALS  Vectorized q(io) = x' * Ma{io} * x  for io = 1:n.
% Collapses the per-constraint cell-array Hermitian-form loop into a single
% sparse mat-vec, mirroring the persistent-cache trick used in MakeJacobianD.
% Uses its OWN persistent cache because Ma uses a different scaling than the
% MS triplet (so MakeJacobianD's cache must NOT be reused for it).
%
% Note: x' is the conjugate transpose, so this computes the Hermitian form
%   q(io) = sum_a conj(x(a)) * (Ma{io}*x)(a)
% identical to the original loop  q(io) = x'*Ma{io}*x.
persistent Va cachedMa nc
initialize = nargin>=4 && strcmp(action,'initialize');
if initialize
    if isempty(cachedMa) || n~=nc || ~isequal(Ma,cachedMa)
        [Va,cachedMa,nc]=build_operator(Ma,n);
    end
    q=[];
    return;
end
if isempty(nc) || n ~= nc
    [Va,cachedMa,nc]=build_operator(Ma,n);
end
R = reshape(Va*x, n, n);    % column io = Ma{io} * x
q = R.' * conj(x);          % q(io) = conj(x).' * (Ma{io}*x) = x' * Ma{io} * x
end

function [Va,cachedMa,nc]=build_operator(Ma,n)
    rows = cell(n,1); cols = cell(n,1); vals = cell(n,1);
    for io = 1:n
        [jj, bb, vv] = find(Ma{io});   % jj = row, bb = col, vv = value
        rows{io} = jj + (io-1)*n;       % stack constraint io into rows (io-1)*n+1 : io*n
        cols{io} = bb;
        vals{io} = vv;
    end
    Va = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), n*n, n);
    cachedMa = Ma;
    nc = n;
end
