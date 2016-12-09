function C = barycentric(x,xx)
% BARYCENTRIC  Spectrally-accurate interpolation of collocation points.
%
% C = barycentric(x,xx)
%
%  x - collocation points [vector, length N+1]
% xx - interpolation points [vector, length M]
%  C - interpolation matrix [size M-by-N+1]
%
% See also DIFFBARY.

  x = x(:);  xx = xx(:);
  N = length(x);   M = length(xx);   one = ones(N-1,N);   C = zeros(M,N);

  X1 = repmat(x',N-1,1);   %  N-1 x N matrix with constant columns

  X2 = zeros(N-1,N);
  
  for j = 1:N,
     X2(:,j)  =  x([ 1:(j-1)  (j+1):N ]);     %  x with entry x(j) deleted
  end
  
% dX(k,j)  =  x(j) - x(k),    if  k <  j
%          =  x(j) - x(k+1),  if  j <= k < N.
  
  dX = X1 - X2;
  
% denom(j) =
%   (x(j)-x(1)) (x(j)-x(2)) ... (x(j)-x(j-1)) (x(j)-x(j+1)) ... (x(j)-x(N))
%
% which is the denominator of the j-th Lagrange polynomial.

  denom = prod(dX);

% C(i,j) = 
%   (xx(i)-x(1)) ... (xx(i)-x(j-1)) (xx(i)-x(j+1)) ... (xx(i)-x(N))
%
% which is the numerator of the j-th Lagrange polynomial evaluated at xx(i).
 
  for i = 1:M,
      C(i,:) = prod(one*xx(i) - X2) ./ denom;
  end

