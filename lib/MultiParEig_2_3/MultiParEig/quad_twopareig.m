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
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% PH 22.11.2016 : added precision-independency.
% PH 26.11.2016 : code simplifications and clean-ups.

% Last revision: 26.11.2016

% Validate number of input parameters.
narginchk(12, 13);

if nargin < 13, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2);
end

opts.singular = 1;
if isfield(opts,'inviter'),  inviter = opts.inviter;    else, inviter = 1;    end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(D1,class_t), D1 = numeric_t(D1,class_t); end;
if ~isa(E1,class_t), E1 = numeric_t(E1,class_t); end;
if ~isa(F1,class_t), F1 = numeric_t(F1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;
if ~isa(D2,class_t), D2 = numeric_t(D2,class_t); end;
if ~isa(E2,class_t), E2 = numeric_t(E2,class_t); end;
if ~isa(F2,class_t), F2 = numeric_t(F2,class_t); end;

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
XR = zeros(m1,neig,class_t); XL = zeros(m1,neig,class_t);
YR = zeros(m2,neig,class_t); YL = zeros(m2,neig,class_t);

% MP: Generate initial vectors (in case of inverse iteration) only
% once. This gives us a bit of speed-up, especially in extended precision
% case, where multiple calls to randn might be noticeable.
if inviter
    x10 = randn(size(A1,1),1,class_t);
    x20 = randn(size(A2,1),1,class_t);    
else
    x10 = numeric_t([],class_t);
    x20 = numeric_t([],class_t);
end;

for i = 1:neig   
    [XR(:,i),XL(:,i)] = min_sing_vec(A1 + lambda(i)* B1 + mu(i)*C1 + lambda(i)^2 * D1+ lambda(i)*mu(i)*E1 + mu(i)^2 * F1,inviter,x10,x10);
    [YR(:,i),YL(:,i)] = min_sing_vec(A2 + lambda(i)* B2 + mu(i)*C2 + lambda(i)^2 * D2+ lambda(i)*mu(i)*E2 + mu(i)^2 * F2,inviter,x20,x20);
end

end