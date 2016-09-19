function c = ce0(x, q, N);

%  The function c = ce0(x, q, N) computes the Mathieu cosine-elliptic 
%  function ce0(x, q) using an NxN Fourier differentiation matrix.
%  The input parameter x could be a scalar or a vector.

%  J.A.C. Weideman, S.C. Reddy 1998
 
[t, D] = fourdif(N,2);                  % Assemble Differentiation Matrix.
 [V,E] = eig((q/2)*diag(cos(t))-D);     % Solve Eigenproblem.     
 
 [m,l] = min(diag(E));                  % Determine charcteristic number.
     v = abs(V(:,l))*sqrt(N/2);         % and corresponding eigenfunction.

     c = fourint(v, 2*x);               % Compute function values at t
                                        % with barycentric trigonometric
                                        % interpolation.
