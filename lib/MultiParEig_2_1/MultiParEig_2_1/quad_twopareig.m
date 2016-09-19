function [lambda,mu,XR,YR,XL,YL] = quad_twopareig(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,opts)

%QUAD_TWOPAREIG  Solve a quadratic two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL] = QUAD_TWOPAREIG(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,opts)
% returns eigenvalues and eigenvectors of the quadratic two-parameter eigenvalue problem
%
% (A1 + lambda B1 + mu C1 + lambda^2 D1+ lambda mu E1 + mu^2 F1)x = 0
% (A2 + lambda B2 + mu C2 + lambda^2 D2+ lambda mu E2 + mu^2 F2)y = 0
%
% Output:
%    - lambda, mu : eigenvalues
%    - XR, YR: right eigenvectors
%    - XL, YL: left eigenvectors
%
% Options in opts:
%   - inviter (1) : use inverse iteration for eigenvectors or slow svd (0)
%   - all options of twopareig and its auxiliary functions

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter eigenvalue 
% problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

opts.singular = 1;
if isfield(opts,'inviter'),  inviter = opts.inviter;    else inviter = 1;    end

% if matrices are in vpa format, we use eig and svd instead of qz and lu
if strcmp(class(A1), 'sym') || strcmp(class(A2), 'sym')
    inviter = 0;
end

% Linearization
[LP1, LQ1, LR1] = linearize_quadtwopar(A1,B1,C1,D1,E1,F1);
[LP2, LQ2, LR2] = linearize_quadtwopar(A2,B2,C2,D2,E2,F2);

% Regular eigenvalues of a singular two-parameter eigenvalue problem
[lambda, mu] = twopareig(LP1,-LQ1,-LR1,LP2,-LQ2,-LR2,opts); 

% We compute eigenvectors one by one using inverse iteration. 
% This works only when all eigenvalues are simple. 
neig = max(size(lambda));
m1 = size(A1,1);
m2 = size(A2,1);
XR = zeros(m1,neig); XL = zeros(m1,neig);
YR = zeros(m2,neig); YL = zeros(m2,neig);

for i = 1:neig   
    [XR(:,i),XL(:,i)] = min_sing_vec(A1 + lambda(i)* B1 + mu(i)*C1 + lambda(i)^2 * D1+ lambda(i)*mu(i)*E1 + mu(i)^2 * F1,inviter);
    [YR(:,i),YL(:,i)] = min_sing_vec(A2 + lambda(i)* B2 + mu(i)*C2 + lambda(i)^2 * D2+ lambda(i)*mu(i)*E2 + mu(i)^2 * F2,inviter);
end
