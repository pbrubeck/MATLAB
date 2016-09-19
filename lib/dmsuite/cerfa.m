function y = cerfa(t,N)

% The function y = cerfa(t,N) computes y(t) = exp(t^2) erfc(t)
% for t > 0 using an NxN Chebyshev differentiation matrix.
% The boundary condition is y = 0 at t = infty.   
% The input parameter may be a scalar or vector.

% J.A.C. Weideman, S.C. Reddy 1998

c  = 3.75;                               % Initialize parameter

[x, D] = chebdif(N+1,1);                 % Compute Chebyshev points, 
D = D(2:N+1,2:N+1);                      % assemble differentiation matrix,
x = x(2:N+1);                            % and incorporate boundary condition

A = diag((1-x).^3)*D-diag(4*c^2*(1+x));  % Coefficient matrix
b = 4*c/sqrt(pi)*(x-1);                  % Right-hand side
y = A\b;                                 % Solve system 

y = chebint([0; y], (t-c)./(t+c));       % Interpolate
