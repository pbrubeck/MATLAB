function [D,x] = mycheb(N,xmax,xmin)
%MYCHEB  Chebyshev spectral collocation scheme in one dimension.
%
% [D,x] = mycheb(N,xmax,xmin)
%
%  x - Chebyshev collocation points in ascending order over [xmin,xmax]
%      [column-vector of length N+1]
%
% D - differentiation matrix
%      [square matrix, size N+1]
%
% See also CHEB, LEG, PCHEB.

  if nargin < 3,  xmin = 0;  end
  if nargin < 2,  xmax = 1;  end

  [Ds,s] = cheb(N);
  k = (xmax - xmin)/2;
  x = xmin + k*(1-s);
  D = -Ds/k;
  