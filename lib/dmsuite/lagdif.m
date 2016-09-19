function [x, DM] = lagdif(N, M, b)

%  The function [x, DM] = lagdif(N, M, b) computes the
%  differentiation matrices D1, D2, ..., DM on Laguerre points.
%
%  Input:
%  N:    Number of points, i.e., order of differentiation matrices (integer).
%  M:    Number of derivatives required (integer).
%  b:    Scaling parameter (real, positive).
%
%  Note:  0 < M < N-1.
%
%  Output:
%  x:    Vector of nodes (zeros of Laguerre polynomial of degree N-1,
%        plus x = 0), all scaled by the parameter b.
%  DM:   DM(1:N,1:N,l) contains l-th derivative matrix, l=1..M.

%  J.A.C. Weideman, S.C. Reddy 1998.
 
x = lagroots(N-1);                    % Compute Laguerre roots
x = [0; x];                           % Add a node at x=0 to facilitate
                                      % the implementation of BCs.

alpha = exp(-x/2);                    % Compute weights.

for ell = 1:M                             % Set up beta matrix s.t.
 beta(ell,:) = (-0.5)^ell*ones(size(x')); % beta(ell,j) is (l^th derivative
end                                       % alpha(x))/alpha(x)
                                          % evaluated at x = x(j).

DM = poldif(x, alpha, beta);          % Compute differentiation matrix (b=1).

x = x/b;                              % Scale nodes by the factor b.

for ell = 1:M                           % Adjust for b not equal to 1.
 DM(:,:,ell) = (b^ell)*DM(:,:,ell);
end
