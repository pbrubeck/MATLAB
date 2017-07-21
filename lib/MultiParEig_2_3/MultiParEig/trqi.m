function [lambda,mu,x,y,res,flag,iter] = trqi(A1,B1,C1,A2,B2,C2,x0,y0,maxiter,tol)

%TRQI  Tensor Rayleigh quotient iteration for a two-parameter eigenvalue problem 
%
% [lambda,mu,x,y,res,flag,iter] = TRQI(A1,B1,C1,A2,B2,C2,x0,y0,maxiter,tol)
% applies Tensor Rayleigh quotient iteration to a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x
% A2 y = lambda B2 y + mu C2 y
%
% Input
%   - A1,B1,C1,A2,B2,C2 : matrices of the two-parameter eigenvalue problem
%   - x0, y0 : approximation for the eigenvector
%   - maxiter : maximum number of iterations (default is 5)
%   - tol : tolerance for the residuals (default is 1e6*eps)
%
% Output:
%   - lambda, mu : eigenvalue
%   - x,y : eigenvector
%   - res : norms of the residual
%   - flag : convergence flag (success 1 or fail 0)
%   - iter : number of iterations
%
% See also: TRQI_3P, TRQI_NP.

% Reference: Algorithm 2 (TRQI) from
% B. Plestenjak, A continuation method for a right definite two-parameter
% eigenvalue problem, SIAM J. Matrix Anal. Appl. 21 (2000), 1163-1184.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 04.12.2016: new stopping criteria when residual stops decreasing
% BP 01.12.2016: modified to be precision-independent                 
% BP 07.11.2016: We consider Frobenius norm of the residuals 
% BP 03.09.2015: Parameter tol is used for both equations
% Last revision: 4.12.2016

% Validate number of input parameters
narginchk(8,10);

class_t = superiorfloat(A1,B1,C1,A2,B2,C2,x0,y0);

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end
if ~isa(x0,class_t), x0 = numeric_t(x0,class_t); end
if ~isa(y0,class_t), y0 = numeric_t(y0,class_t); end

if nargin<9 || isempty(maxiter), maxiter = 5; end
if nargin<10 || isempty(tol), tol = numeric_t('1e6*eps',class_t); end

iter = 0;

x = x0/norm(x0);
y = y0/norm(y0);
normres = 1e20;
oldnormres = Inf;

A1x = A1*x; B1x = B1*x; C1x = C1*x;
A2y = A2*y; B2y = B2*y; C2y = C2*y;

while (iter<= maxiter) && (normres>tol) && (normres<0.95*oldnormres)
    
    M = [x'*B1x x'*C1x; y'*B2y y'*C2y];
    b = [x'*A1x; y'*A2y];
    warning off
    tmp = M\b;
    lambda = tmp(1); 
    mu = tmp(2);

    iter = iter + 1;
    V1 = (A1-lambda*B1-mu*C1)\[B1x C1x];  
    V2 = (A2-lambda*B2-mu*C2)\[B2y C2y];
    warning on
    
    if any(any(isnan([V1; V2]))) || any(any(isinf(double([V1; V2]))))
        % in this case we compute vectors from SVD 
        [tilda1,tilda2,tmpV] = svd(A1-lambda*B1-mu*C1); x = tmpV(:,end);
        [tilda1,tilda2,tmpV] = svd(A2-lambda*B2-mu*C2); y = tmpV(:,end);
    else
        M = [x'*V1; y'*V2];
        warning off
        delta = M\[1;1];
        warning on
        xnew = V1*delta;
        ynew = V2*delta;
        x = xnew/norm(xnew);
        y = ynew/norm(ynew);
    end
    
    A1x = A1*x; B1x = B1*x; C1x = C1*x;
    A2y = A2*y; B2y = B2*y; C2y = C2*y;
    
    res1 = norm(A1x-lambda*B1x-mu*C1x);
    res2 = norm(A2y-lambda*B2y-mu*C2y);
    oldnormres = normres;
    normres = norm([res1 res2]);
end    
res = [res1 res2];

if normres<tol, flag=1; else flag=0; end
    
