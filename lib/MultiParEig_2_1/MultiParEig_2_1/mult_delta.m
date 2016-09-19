function [G0,G1,G2] = mult_delta(A1,B1,C1,A2,B2,C2,Z,W)

%MULT_DELTA  Projections of Delta matrices for two-parameter eigenvalue problem
%
% [G0,G1,G2] = MULT_DELTA(A1,B1,C1,A2,B2,C2,Z,W) returns matrices G0, G1, G2, 
% such that Gi = Z'*Deltai*W for i=0,1,2, where
% 
% Delta0 = kron(B1,C2) - kron(C1,B2),
% Delta1 = kron(A1,C2) - kron(C1,A2),
% Delta2 = kron(B1,A2) - kron(A1,B2),
%
% without explicitly computing Delta matrices

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

m = size(W,2);
n = size(Z,2);
k1 = size(A1,1);
k2 = size(A2,1);
G0 = zeros(n,m);
G1 = zeros(n,m);
G2 = zeros(n,m);
for j=1:m
    MW = reshape(W(:,j),k2,k1);
    AMW2 = A2*MW;
    BMW2 = B2*MW;
    CMW2 = C2*MW;
       
    D0 = CMW2*B1.' - BMW2*C1.';
    D1 = CMW2*A1.' - AMW2*C1.';
    D2 = AMW2*B1.' - BMW2*A1.';
        
    d0 = reshape(D0, k1*k2, 1);
    d1 = reshape(D1, k1*k2, 1);
    d2 = reshape(D2, k1*k2, 1);
    
    G0(:,j) = Z'*d0;
    G1(:,j) = Z'*d1;
    G2(:,j) = Z'*d2;
end