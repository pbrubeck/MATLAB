function r = herroots(N);

%  The function r = herroots(N) computes the roots of the 
%  Hermite polynomial of degree N.

%  J.A.C. Weideman, S.C. Reddy 1998.
 
J = diag(sqrt([1:N-1]),1)+diag(sqrt([1:N-1]),-1);    % Jacobi matrix
r = sort(eig(sparse(J)))/sqrt(2);                    % Compute eigenvalues
