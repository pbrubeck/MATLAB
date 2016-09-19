function [lambda,mu,X,Y,flag] = twopareigs(A1,B1,C1,A2,B2,C2,k,opts)

%TWOPAREIGS   Find a few eigenvalues for a two-parameter eigenvalue problem
%
% [lambda,mu,X,Y] = TWOPAREIGS(A1,B1,C1,A2,B2,C2,k,opts)
% returns k eigenvalues with the smallest |mu| of the two-parameter
% eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using either implicitly restarted inverse Arnoldi (twopareigs_ira) or 
% Krylov-Schur method (twopareigs_ks) on the generalized eigenvalue
% problem Delta2*z = lambda*Delta0*z, where z = kron(x,y) and
% Delta0 = kron(B1,C2) - kron(C1,B2)
% Delta2 = kron(B1,A2) - kron(A1,B2)
% 
% When building Krylov subspace, Bartels-Stewart method is used to solve 
% the Sylvester equation related to linear system Delta2*w = Delta0*z
%
% By default the Krylov-Schur method is used when all matrices are real
% and implicitly restarted Arnoldi otherwise
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : matrices of the two-parameter eigenvalue problem
%   - k: number of eigenvalues (6)
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - X, Y: components of decomposable right eigenvectors
%     (eigenvector is kron(X(:,j),Y(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)*X(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*Y(:,j)=0
%   - flag: convergence (0), no convergence (1)
% 
% Options in opts:
%   - use_ira: (default 0) set to 1 to use IRA for real matrices  
%   - all other options of twopareigs_ira or twopareigs_ks, like 
%       - divA : (default 1) : divide by A1 and A2, or by B1 and B2 (0)
%       - refine : number of TRQ refinement steps in the end (1)
%       - epscluster : distance between eigenvalues in the same cluster (1e-4)
%       - v0 : initial vector 
%       - tol : absolute tolerance 
%       - maxit : maximum number of iterations 
%       - p : minimal size of search space 
%       - disp : display progress 
%
% For Matlab below 2014a package lapack from MatlabCentral is required 
% for faster evaluation 
%
% See also: TWOPAREIG, TWOPAREIGS_IRA, TWOPAREIGS_KS, TWOPAREIGS_JD,
% TWOPAREIGS_SI, THREEPAREIGS

% Reference: K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method 
% for the generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<7, k=6; end      % if not specified, 6 eigenvalues are returned 
if nargin<8, opts=[]; end
if isfield(opts,'use_ira'), use_ira = opts.use_ira; else use_ira = 0; end

n1 = size(A1,1);
n2 = size(A2,1);

% if k>=n1*n2 we compute all eigenvalues 
if k>=n1*n2-1
    [lambda,mu,X,Y] = twopareig(A1,B1,C1,A2,B2,C2);
    flag = 0;
    return
end

if isreal(A1) && isreal(B1) && isreal(C1) && isreal(A2) && isreal(B2) && isreal(C2) && ~use_ira
    [lambda,mu,X,Y,flag] = twopareigs_ks(A1,B1,C1,A2,B2,C2,k,opts);
else
    [lambda,mu,X,Y,flag] = twopareigs_ira(A1,B1,C1,A2,B2,C2,k,opts);
end

