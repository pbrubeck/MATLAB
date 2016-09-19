function y = cerfb(t,N)

% The function y = cerfb(t,N) computes y(t) = exp(t^2) erfc(t)
% for t > 0 using an NxN Chebyshev differentiation matrix.
% The boundary condition is y = 1 at t = 0.   
% The input parameter may be a scalar or vector.

% J.A.C. Weideman, S.C. Reddy 1998

c  = 3.75;                               % Initialize parameter

[x, D] = chebdif(N+1,1);                 % Compute Chebyshev points, 

A = diag((1-x).^3)*D-diag(4*c^2*(1+x));  % Coefficient matrix
b = 4*c/sqrt(pi)*(x-1);                  % Right-hand side

a1 = A(1:N,N+1); b = b(1:N);
A = A(1:N,1:N); 
y = A\(b-a1);                            % Solve system 

y = chebint([y; 1], (t-c)./(t+c));       % Interpolate
