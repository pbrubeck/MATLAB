function [lambda,mu,eta,x,y,z,res,flag,iter] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,x0,y0,z0,maxiter,tol)

%TRQI_3P   Tensor Rayleigh quotient iteration for a three-parameter eigenvalue problem 
%
% [lambda,mu,x,y,res,flag,iter] = trqi(A1,B1,C1,A2,B2,C2,x0,y0,maxiter,tol)
% applies Tensor Rayleigh quotient iteration to a three-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x + eta D1 x
% A2 y = lambda B2 y + mu C2 y + eta D2 y
% A3 z = lambda B3 z + mu C3 z + eta D3 z
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices of the three-parameter eigenvalue problem
%   - x0, y0, z0 : approximation for the eigenvector
%   - maxiter : maximum number of iterations (default is 5)
%   - tol : tolerance for the residual (default is 1e-10)
%
% Output:
%   - lambda, mu, eta : eigenvalue
%   - x,y,z : eigenvector
%   - res : norms of the residual
%   - flag : convergence flag (success 1 or fail 0)
%   - iter : number of iterations
%
% See also: TRQI, TRQI_NP

% Reference: Algorithm 2 (TRQI) from
% B. Plestenjak, A continuation method for a right definite two-parameter
% eigenvalue problem, SIAM J. Matrix Anal. Appl. 21 (2000), 1163-1184.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<16  || isempty(maxiter), maxiter = 5; end
if nargin<17  || isempty(tol), tol = 1e-10; end

ACell = cell(3,4);
ACell{1,1} = A1; ACell{1,2} = B1; ACell{1,3} = C1; ACell{1,4} = D1;
ACell{2,1} = A2; ACell{2,2} = B2; ACell{2,3} = C2; ACell{2,4} = D2;
ACell{3,1} = A3; ACell{3,2} = B3; ACell{3,3} = C3; ACell{3,4} = D3;

XCell = cell(3);
XCell{1} = x0; XCell{2} = y0; XCell{3} = z0;

[VecLambda,XCell,res,flag,iter] = trqi_np(ACell,XCell,maxiter,tol);

x = XCell{1}; y = XCell{2}; z = XCell{3};
lambda = VecLambda(1); mu = VecLambda(2); eta = VecLambda(3);


