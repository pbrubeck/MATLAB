%DEMO_QUAD_TWOPAREIG     Demo for quad_twopareig
%
% We solve the quadratic two-parameter eigenvalue problem
%
% (A1 + lambda B1 + mu C1 + lambda^2 D1+ lambda mu E1 + mu^2 F1)x = 0
% (A2 + lambda B2 + mu C2 + lambda^2 D2+ lambda mu E2 + mu^2 F2)y = 0
%
% where
%
% A1=[3 4;6 1]; B1=[1 2;2 1]; C1=[4 1;2 4]; D1=[6 7;5 2]; E1=[1 3;7 1]; F1=[4 1;6 3];
% A2=[1 3;2 1]; B2=[1 4;8 2]; C2=[2 3;4 1]; D2=[2 6;3 1]; E2=[7 2;3 7]; F2=[3 5;5 2];
%
% See also: QUAD_TWOPAREIG

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

A1=[3 4;6 1]; B1=[1 2;2 1]; C1=[4 1;2 4]; D1=[6 7;5 2]; E1=[1 3;7 1]; F1=[4 1;6 3];
A2=[1 3;2 1]; B2=[1 4;8 2]; C2=[2 3;4 1]; D2=[2 6;3 1]; E2=[7 2;3 7]; F2=[3 5;5 2];

% we compute eigenvalues and eigenvectors
[lambda,mu,XR,YR,XL,YL] = quad_twopareig(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2);
eigenvalues = [lambda mu]

% for a check that this are really eigenvalues, we compute minimal singular
% values of the matrices
minsvd = [];
for k = 1:length(lambda)
    minsvd(k,1) = min(svd(A1 + lambda(k)*B1 + mu(k)*C1 + lambda(k)^2*D1 + lambda(k)*mu(k)*E1 + mu(k)^2*F1));
    minsvd(k,2) = min(svd(A2 + lambda(k)*B2 + mu(k)*C2 + lambda(k)^2*D2 + lambda(k)*mu(k)*E2 + mu(k)^2*F2));
end
minsvd

% for a check that we really have eigenpairs, we compute the norms of the
% residuals for right eigenvector
res = [];
for k = 1:length(lambda)
    res(k,1) = norm((A1 + lambda(k)*B1 + mu(k)*C1 + lambda(k)^2*D1 + lambda(k)*mu(k)*E1 + mu(k)^2*F1)*XR(:,k));
    res(k,2) = norm((A2 + lambda(k)*B2 + mu(k)*C2 + lambda(k)^2*D2 + lambda(k)*mu(k)*E2 + mu(k)^2*F2)*YR(:,k));
end
res

