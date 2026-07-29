function [Mp,Mq] = get_quadr_mtrx(Y,n)
%GET_QUADR_MTRX Build sparse active/reactive power quadratic matrices.

if ~isequal(size(Y),[n n])
    error('HEBCPF:get_quadr_mtrx:DimensionMismatch', ...
        'Ybus must be an n-by-n matrix.');
end

Y = sparse(Y);
G = real(Y);
B = imag(Y);
Mp = cell(n,1);
Mq = cell(n,1);

for j = 1:n
    Ap = sparse(2*n,2*n);
    Ap(j,1:n) = G(j,:);
    Ap(j,n+1:2*n) = -B(j,:);
    Ap(j+n,1:n) = B(j,:);
    Ap(j+n,n+1:2*n) = G(j,:);
    Mp{j} = 0.5*(Ap + Ap.');

    Aq = sparse(2*n,2*n);
    Aq(j,1:n) = -B(j,:);
    Aq(j,n+1:2*n) = -G(j,:);
    Aq(j+n,1:n) = G(j,:);
    Aq(j+n,n+1:2*n) = -B(j,:);
    Mq{j} = 0.5*(Aq + Aq.');
end
end
