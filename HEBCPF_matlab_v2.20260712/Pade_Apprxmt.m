function [num,den]=Pade_Apprxmt(Vx,Vy,degree,bus_n)
% Pade approximation based on the Taylor series of complex voltage variables.
%
% STRUCTURAL OPTIMIZATION (equivalent to the original full-matrix solve):
% The original builds, per bus, the (act+1)x(act+1) matrix G=[G1 G2] and solves
% G*nd=-c. With G1=-[eye(p+1);0] the bottom q=denominator rows contain ONLY the
% denominator unknowns, so G is block lower-triangular. We therefore:
%   (1) solve the q x q Toeplitz denominator block  A*d = -c(p+1:p+q),
%   (2) recover the numerator by the convolution  N = (c * [1;d])  truncated to p.
% This is the SAME linear system (no approximation): it replaces a per-bus
% (act+1) dense solve + an inner G2-build loop with a per-bus q-solve and a
% fully vectorised (across buses) numerator step. Results match the original to
% round-off.

p = degree.numerator;        % numerator degree  (coeffs 0..p)
q = degree.denominator;      % denominator degree (coeffs 1..q, D0=1)
M = degree.act;              % highest Taylor coeff index (coeffs 0..M, M=p+q)

ap = Vx(:,1:M+1) + 1i*Vy(:,1:M+1);   % bus_n x (M+1); ap(:,k) = c_{k-1}

% --- Denominator: per-bus q x q Toeplitz solve --------------------------------
% Equation rows i = p+1..p+q (only denominator unknowns appear there):
%   sum_{j=1}^{q} c_{i-j} d_j = -c_i
% A(r,j) = c_{p+r-j} = ap(:, p+r-j+1);  rhs(r) = -c_{p+r} = -ap(:, p+r+1)
[rr,jj] = ndgrid(1:q,1:q);
colIdxA = p + rr - jj + 1;           % q x q index template into the coeff array
rhsIdx  = (p+2):(p+q+1);             % 1 x q  -> c_{p+1..p+q}

dcoef = zeros(bus_n, q) + 0i;
for b = 1:bus_n
    c = ap(b,:);
    A = c(colIdxA);                  % q x q (same shape as colIdxA)
    dcoef(b,:) = (A \ (-c(rhsIdx).')).';
end

den = [ones(bus_n,1), dcoef];        % bus_n x (q+1), leading 1

% --- Numerator: N_k = sum_{j=0}^{min(k,q)} c_{k-j} D_j , k=0..p ----------------
% Vectorised across buses; D_0=1, D_1..D_q = dcoef.
num = ap(:,1:p+1);                   % j=0 term: c_k
for j = 1:q
    k = j:p;                         % k>=j contribute c_{k-j}*d_j
    if isempty(k); break; end
    num(:,k+1) = num(:,k+1) + ap(:,k-j+1) .* dcoef(:,j);
end

end
