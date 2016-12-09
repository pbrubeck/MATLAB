function [D,x] = cheb(N)
% CHEB  Trefethen's Chebyshev spectral collocation scheme over [-1,1].
%
% [D,x] = cheb(N)
%
%  x - Chebyshev collocation points in descending order
%      [column-vector of length N+1]
%
% D - differentiation matrix
%      [square matrix, size N+1]
%
% See also MYCHEB, LEG, PCHEB.

if N==0, D=0; x=1; return, end

x = cos(pi*(0:N)/N)';
c = [2 ones(1,N-1) 2] .* ((-1).^(0:N));  % c = [2 -1 1 -1 ... 1 -1 2]
X = repmat(x,1,N+1);
D_numer = c' * (1./c);
D_denom = (X - X') + eye(N+1);
D = D_numer./D_denom;
D = D - diag(sum(D'));
