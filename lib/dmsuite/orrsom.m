
%  The script file orrsom.m computes the eigenvalues of the Orr-Sommerfeld
%  equation using NxN Chebyshev differentiation matrices.

% S.C. Reddy, J.A.C. Weideman 1998.  Code modified to display output
% by JACW, May 2003.

N = input(' Order of the differentiation matrix: N = ? ');
R = input(' Reynolds number: R = ? ');
i = sqrt(-1);

[x,DM] = chebdif(N+2,2);                       % Compute second derivative
    D2 = DM(2:N+1,2:N+1,2);                    % Enforce Dirichlet BCs
                                     
[x,D4] = cheb4c(N+2);                          % Compute fourth derivative
     I = eye(size(D4));                        % Identity matrix

A = (D4-2*D2+I)/R-2*i*I-i*diag(1-x.^2)*(D2-I); % Set up A and B matrices
B = D2-I;

e = eig(A,B);                                  % Compute eigenvalues

[m,l] = max(real(e));                          % Find eigenvalue of largest
disp('Eigenvalue with largest real part = ')   % real part
disp(e(l))
