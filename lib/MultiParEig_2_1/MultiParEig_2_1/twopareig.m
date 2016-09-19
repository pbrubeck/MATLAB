function [lambda,mu,XR,YR,XL,YL] = twopareig(A1,B1,C1,A2,B2,C2,opts)

%TWOPAREIG   Solve a two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL] = TWOPAREIG(A1,B1,C1,A2,B2,C2,opts) returns
% eigenvalues and eigenvectors of the two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : real matrices
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR: components of decomposable right eigenvectors
%     (eigenvector is kron(XR(:,j),YR(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL: components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),YL(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%
% Operator determinants Delta0, Delta1, and Delta2 are used, where
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - epscluster (1e-4): relative distance between eigenvalues in a cluster
%   - fast (1): use fast algorithm (can fail for multiple eigenvalues) 
%     or slow algorithm (0) with clustering
%   - inviter (1) : use inverse iteration for eigenvectors or slow svd (0)
%   - all options of auxiliary functions
%
% See also: TWOPAREIGS, TWOPAREIGS_JD, TWOPAREIGS_SI, THREEPAREIG, MULTIPAREIG

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : use extract_regular_part_np
% Last revision: 8.9.2015

if nargin<7, opts=[]; end
if isfield(opts,'epscluster'),  epscluster = opts.epscluster;   else epscluster = 1e-4;   end
if isfield(opts,'fast'),        fast = opts.fast;               else fast = 1;            end
if isfield(opts,'inviter'),     inviter = opts.inviter;         else inviter = 1;         end
if isfield(opts,'singular'),    singular = opts.singular;       else singular = 0;        end

% if matrices are in vpa format, we use eig and svd instead of qz and lu
symb = strcmp(class(A1), 'sym');
if symb
    fast = 1; 
    inviter = 0;
end

% Delta matrices 
[Delta0, Delta1, Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2);

if singular
    DeltaCell = extract_regular_part_np({Delta0,Delta1,Delta2}, opts);
    Delta0 = DeltaCell{1}; Delta1 = DeltaCell{2}; Delta2 = DeltaCell{3}; 
end

if (size(Delta0,1) == 0) || (size(Delta0,1)~=size(Delta0,2))   % no regular part was found
     lambda = []; mu = []; XL = []; XR = []; YR = []; YL = []; % default output
     return
end

if fast
    tmp = Delta0\[Delta1 Delta2];
    n = size(Delta0,1);
    Gamma1 = tmp(:,1:n);
    Gamma2 = tmp(:,n+1:end); 
    if symb 
        [Q1,D1] = eig(Gamma1);
        lambda = diag(D1);
        mu = diag(Q1\Gamma2*Q1);
    else
        [Q1,D1] = schur(Gamma1,'complex');
        lambda = diag(D1);
        mu = diag(Q1'*Gamma2*Q1);
    end
else
    mu = [];
    [S0,S1,Q,Z,order,start,csize,lambda] = clustered_qz(Delta0,Delta1,epscluster);
    S2 = Q*Delta2*Z;
    for k = 1:length(start)
        partS0 = S0(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
        partS2 = S2(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
        partmu = eig(partS2,partS0);
        mu = [mu; partmu];
    end
end

if nargout > 2
    % extraction of eigenvectors (individually using inverse iteration or SVD)
    for k = 1:length(lambda)
        [xr,xl] = min_sing_vec(A1-lambda(k)*B1-mu(k)*C1,inviter);
        XR(:,k) = xr; XL(:,k)=xl;
        [yr,yl] = min_sing_vec(A2-lambda(k)*B2-mu(k)*C2,inviter);
        YR(:,k) = yr; YL(:,k)=yl;
    end   
end

