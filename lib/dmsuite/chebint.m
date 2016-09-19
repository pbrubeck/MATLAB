function p = chebint(fk, x)

%  The function p = chebint(fk, x) computes the polynomial interpolant
%  of the data (xk, fk), where xk are the Chebyshev nodes.  
%  Two or more data points are assumed.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at Chebyshev points 
%       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  x:   Vector of x-values where polynomial interpolant is to be evaluated.
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)

%  J.A.C. Weideman, S.C. Reddy 1998

  fk = fk(:); x = x(:);                    % Make sure data are column vectors.

   N = length(fk); 
   M = length(x);
     
  xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    % Compute Chebyshev points.

   w = ones(N,1).*(-1).^[0:N-1]';          % w = weights for Chebyshev formula
w(1) = w(1)/2; w(N) = w(N)/2;
 
   D = x(:,ones(1,N)) - xk(:,ones(1,M))';  % Compute quantities x-x(k)
   D = 1./(D+eps*(D==0));                  % and their reciprocals.
  
   p = D*(w.*fk)./(D*w);                   % Evaluate interpolant as
                                           % matrix-vector products.
