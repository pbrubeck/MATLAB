function [lambda,mu,eta,X,Y,Z] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,k,opts)

%THREEPAREIGS   Find a few eigenvalues for a three-parameter eigenvalue problem
%
% [lambda,mu,eta,X,Y,Z] = THREEPAREIGS(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues with the smallest |eta| of the three-parameter
% eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using implicit restarted inverse Arnoldi in Matlab routine eigs
% on the related generalized eigenvalue problems Delta3*w = eta*Delta0*w, 
% where Delta0 and Delta3 are corresponding operator determinants
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
%   - k : number of eigenvalues
%   - opts : options (see below)
%
% Output:
%   - lambda , mu, eta: eigenvalues (eigenvalue is (lambda(j),mu(j),eta(j))
%   - X1, X2, X3 : components of decomposable right eigenvectors 
%     (eigenvector is kron(X1(:,j),kron(X2(:,j),X3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)*X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)*X2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)*X3(:,j)=0
%   - Y1, Y2, Y3 : components of decomposable left eigenvectors 
%     (eigenvector is kron(Y1(:,j),kron(Y2(:,j),Y3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)'*Y1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)'*Y2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)'*Y3(:,j)=0
%
% Options in opts:
%   - usesparse : set to 0 if all matrices are full. The default (1) uses 
%     sparse representation of Delta matrices, which is better if some 
%     matrices are diagonal or zero
%   - refine : number of TRQ refinement steps in the end (2)
%   - all options for the eigs
%
% See also: THREEPAREIG, THREEPAREIGS_JD, THREEPAREIGS_SI, TWOPAREIGS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : new option usesparse
% BP 03.09.2015 : refine with trqi_3p
% Last revision: 8.9.2015

% options 
if nargin<14, opts = []; end
if isfield(opts,'usesparse'),  usesparse = opts.usesparse;  else usesparse = 1;  end
if isfield(opts,'refine'),     refine = opts.refine;  else refine = 2;           end

n1 = size(A1,1);
n2 = size(A2,1);
n3 = size(A3,1);

% if k>=n1*n2*n3 we compute all eigenvalues 
if k>=n1*n2*n3-1
    [lambda,mu,eta,X,Y,Z] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
    return
end

if usesparse 
    A1s = sparse(A1); B1s = sparse(B1); C1s = sparse(C1); D1s = sparse(D1);
    A2s = sparse(A2); B2s = sparse(B2); C2s = sparse(C2); D2s = sparse(D2);
    A3s = sparse(A3); B3s = sparse(B3); C3s = sparse(C3); D3s = sparse(D3);
    [Delta0,Delta1,Delta2,Delta3] = threepar_delta(A1s,B1s,C1s,D1s,A2s,B2s,C2s,D2s,A3s,B3s,C3s,D3s);
else
    [Delta0,Delta1,Delta2,Delta3] = threepar_delta(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
end

opts.disp = 0; % turn off diagnostic information level in eigs in earlier versions of Matlab
[W,D] = eigs(Delta3,Delta0,k,'SM',opts);
eta = diag(D);

% lambda and mu are computed from Rayleigh quotients
for i=1:k
    imen = (W(:,i)'*Delta0*W(:,i));
    mu(i,1) = W(:,i)'*Delta2*W(:,i)/imen;
    lambda(i,1) = W(:,i)'*Delta1*W(:,i)/imen;
end

if nargout > 3
    % extraction of eigenvectors (individually using inverse iteration or SVD)
    for k=1:k
        [xr,xl] = min_sing_vec(A1-lambda(k)*B1-mu(k)*C1-eta(k)*D1);
        X(:,k) = xr;
        [xr,xl] = min_sing_vec(A2-lambda(k)*B2-mu(k)*C2-eta(k)*D2);
        Y(:,k) = xr;
        [xr,xl] = min_sing_vec(A3-lambda(k)*B3-mu(k)*C3-eta(k)*D3);
        Z(:,k) = xr;
        if refine>0
            [refl,refu,refe,refx,refy,refz] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,X(:,k),Y(:,k),Z(:,k),refine,eps);
            X(:,k) = refx;  Y(:,k) = refy;  Z(:,k) = refz;
            lambda(k) = refl;  mu(k) = refu;  eta(k) = refe;
        end
    end   
end

