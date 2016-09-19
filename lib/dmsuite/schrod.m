
%  The script file schrod.m computes the first eigenvalue of the Schrodinger 
%  equation on the half-line using an NxN Laguerre differentiation matrix.

% J.A.C. Weideman, S.C. Reddy 1998. Code modified to display output
% by JACW, May 2003.

N = input(' Order of the differentiation matrix: N = ? ');
b = input(' Scaling parameter of the Laguerre method: b = ? ');

r = 5.08685476; epsi = 0.929852862;   % Parameters for Woods-Saxon potential

[x,D] = lagdif(N+1,2,b);              % Compute Laguerre  derivative matrix,
D2 = D(2:N+1,2:N+1,2);                % and enforce boundary conditions.
x  = x(2:N+1);

Q = diag(1./(1+exp((x-r)/epsi)));     % Woods-Saxon potential
I = eye(size(D2));                    % Identity matrix
 
e = min(eig(-D2+I,Q));                % Compute smallest eigenvalue

disp('Smallest eigenvalue = ')        % real part
disp(e)
