function p = polint(xk, fk, x, alpxk, alpx)

%  The function p = polint(xk, fk, x) computes the polynomial interpolant
%  of the data (xk, fk).  Two or more data points are assumed.
%
%  Input  (constant weight)
%  xk:    Vector of x-coordinates of data (assumed distinct).
%  fk:    Vector of y-coordinates of data.
%  x:     Vector of x-values where polynomial interpolant is to be evaluated.
%
%  In the case of non-constant weights two additional input parameters are
%  needed, as in  p = polint(xk, fk, x, alpxk, alpx)
%
%  Input  (non-constant weight)
%  xk:    Vector of x-coordinates of data (assumed distinct).
%  fk:    Vector of y-coordinates of data.
%  x:     Vector of x-values where polynomial interpolant is to be evaluated.
%  alpxk: Vector of weight values sampled at the points xk.
%  alpx:  Vector of weight values sampled at the points x.
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)
%  Note added May 2003:  Except for certain nice node distributions
%  polynomial interpolation of high-degree is an ill-conditioned
%  problem.  This code does not test for conditioning so use with
%  care. 

%  J.A.C. Weideman, S.C. Reddy 1998.


if nargin == 3
alpx = 1; alpxk = 1;
elseif nargin > 3
fk = fk./alpxk;
end

     x = x(:);                          % Make sure the data are column vectors
    xk = xk(:);  fk = fk(:);
 alpxk = alpxk(:); alpx = alpx(:);

     N = length(xk); 
     M = length(x);
     L = logical(eye(N));

     D = xk(:,ones(1,N))-xk(:,ones(1,N))';  % Compute the weights w(k)
  D(L) = ones(N,1);
     w = 1./prod(D)';                       
 
     D = x(:,ones(1,N)) - xk(:,ones(1,M))'; % Compute quantities x-x(k) 
     D = 1./(D+eps*(D==0));                 % and their reciprocals. 
  
     p = alpx.*(D*(w.*fk)./(D*w));          % Evaluate interpolant as
                                            % matrix-vector products.
