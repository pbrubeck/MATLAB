function r = lagroots(N);

%  The function r = lagroots(N) computes the roots of the 
%  Laguerre polynomial of degree N.

%  J.A.C. Weideman, S.C. Reddy 1998.

J = diag([1:2:2*N-1])-diag([1:N-1],1)-diag([1:N-1],-1);  % Jacobi matrix
r = sort(eig(sparse(J)));                                % Compute eigenvalues
