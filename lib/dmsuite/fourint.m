function t = fourint(fk, x)

%  The function t = fourint(fk, x) computes the trigonometric interpolant
%  of the data (xk, fk), where xk are equidistant nodes.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at equidistant points 
%       x(k) = (k-1)*2*pi/N,  k = 1...N.
%  x:   Vector of x-values where interpolant is to be evaluated.
%
%  Output:
%  t:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 46 in
%  P. Henrici, Applied and Computational Complex Analysis III, Wiley, 1986.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)

%  J.A.C. Weideman, S.C. Reddy 1998


 fk = fk(:); x = x(:);           % Make sure data are column vectors.

  N = length(fk); 
  M = length(x);

 xk = (2*pi/N)*[0:N-1]';         % Compute equidistant points
   
  w = (-1).^[0:N-1]';            % Weights for trig interpolation

 x2 = x/2; xk2 = xk/2;
   
  D = x2(:,ones(1,N)) - xk2(:,ones(1,M))'; % Compute quantities x-x(k)

if rem(N,2) == 0                  
  D = 1./tan(D+eps*(D==0));       %  Formula for N even
else
  D = 1./sin(D+eps*(D==0));       %  Formula for N odd
end

  t = D*(w.*fk)./(D*w);            % Evaluate interpolant as 
                                   % matrix-vector products.
