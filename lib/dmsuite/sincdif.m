function [x, DM] = sincdif(N, M, h)

%  The function [x, DM] = sincdif(N, M, h) computes sinc the
%  differentiation matrices D1, D2, ..., DM on equidistant points.
%
%  Input:
%  N:    Number of points, i.e., order of differentiation matrix.
%  M:    Number of derivatives required (integer).
%  h:    Step-size (real, positive).
%
%  Note:  0 < M < N-1.
%
%  Output:
%  x:    Vector of nodes.
%  DM:   DM(1:N,1:N,l) contains l-th derivative matrix, l=1..M.

%  J.A.C. Weideman, S.C. Reddy 1998.  Help lines corrected 
%  by JACW, March/April 2003.

k = [1:N-1]';
t = k*pi;
x = h*[-(N-1)/2:(N-1)/2]';

sigma = zeros(size(k));

for ell = 1:M;
    sigma = (-ell*sigma + imag(exp(i*t)*i^ell))./t;
      col = (pi/h)^ell*[imag(i^(ell+1))/(ell+1); sigma];
      row = (-1)^ell*col; row(1) = col(1);
DM(:,:,ell) = toeplitz(col,row);
end
