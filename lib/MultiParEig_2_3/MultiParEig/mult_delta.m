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

% BP 02.12.2016 : modified to be precision-independent
% Last revision: 8.9.2015

narginchk(8,8);
class_t = superiorfloat(A1,B1,C1,A2,B2,C2);

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;
if ~isa(Z,class_t), Z = numeric_t(Z,class_t); end;
if ~isa(W,class_t), W = numeric_t(W,class_t); end;

m = size(W,2);
n = size(Z,2);
k1 = size(A1,1);
k2 = size(A2,1);
G0 = zeros(n,m,class_t);
G1 = zeros(n,m,class_t);
G2 = zeros(n,m,class_t);
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