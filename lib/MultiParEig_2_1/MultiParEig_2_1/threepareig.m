function [lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,opts)

%THREEPAREIG   Solve a three-parameter eigenvalue problem
%
% [lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3] = THREEPAREIG(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,opts) 
% finds eigenvalues and eigenvectors of a three-parameter eigenvalue problem
%
% A1 x1 = lambda B1 x1 + mu C1 x1 + eta D1 x1 
% A2 x2 = lambda B2 x2 + mu C2 x2 + eta D2 x2 
% A3 x3 = lambda B3 x3 + mu C3 x3 + eta D3 x3 
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
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
% Operator determinants Delta0, Delta1, Delta2, Delta3 are used, where
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta1 = - | A1 C1 D1; A2 C2 D2; A3 C3 D3 |
% Delta2 =   | B1 A1 D1; B2 A2 D2; B3 A3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - epscluster (1e-4): relative distance between eigenvalues in a cluster
%   - fast (1): use fast algorithm (can fail for multiple eigenvalues) 
%     or slow algorithm (0) with clustering
%   - inviter (1) : use inverse iteration for eigenvectors or slow svd (0)
%
% See also: TWOPAREIG, MULTIPAREIG, THREEPAREIGS, THREEPAREIGS_JD,
% THREEPAREIGS_SI

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 6.9.2105 : support for singular 3EP
% Last revision: 8.9.2015

if nargin<13, opts=[]; end
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
[Delta0,Delta1,Delta2,Delta3] = threepar_delta(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);

if singular
    DeltaCell = extract_regular_part_np({Delta0,Delta1,Delta2,Delta3}, opts);
    Delta0 = DeltaCell{1}; Delta1 = DeltaCell{2}; Delta2 = DeltaCell{3}; Delta3 = DeltaCell{4};
end

if fast
    tmp = Delta0\[Delta1 Delta2 Delta3];
    n = size(Delta0,1);
    Gamma1 = tmp(:,1:n);
    Gamma2 = tmp(:,n+1:2*n); 
    Gamma3 = tmp(:,2*n+1:end); 
    if symb 
        [Q1,L1] = eig(Gamma1);
        lambda = diag(L1);
        mu = diag(Q1\Gamma2*Q1);
        eta = diag(Q1\Gamma3*Q1);
    else
        [Q1,L1] = schur(Gamma1,'complex');
        lambda = diag(L1);
        mu =  diag(Q1'*Gamma2*Q1);
        eta = diag(Q1'*Gamma3*Q1);
    end
else
    % clustering lambda
    [S0,S1,Q,Z,order,start,csize,lambda] = clustered_qz(Delta0,Delta1,epscluster);
    S2 = Q*Delta2*Z;
    S3 = Q*Delta3*Z;
    for k=1:length(start)
         % clustering each block (lambda,mu)
         partS0 = S0(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
         partS2 = S2(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
         partS3 = S3(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
         [partS0,partS2,Q2,Z2,order2,start2,csize2,partmu] = clustered_qz(partS0,partS2,epscluster);
         partS3 = Q2*partS3*Z2;
         mu=[mu; partmu];
         for j=1:length(start2)
             % computing eta
             blockS0 = partS0(start2(j):(start2(j)+csize2(j)-1),start2(j):(start2(j)+csize2(j)-1));
             blockS3 = partS3(start2(j):(start2(j)+csize2(j)-1),start2(j):(start2(j)+csize2(j)-1));
             parteta = eig(blockS3,blockS0);
             eta = [eta; parteta];
         end
     end
end

if nargout > 2
    % extraction of eigenvectors (individually using inverse iteration or SVD)
    for k=1:length(lambda)
        [tx1,ty1] = min_sing_vec(A1-lambda(k)*B1-mu(k)*C1-eta(k)*D1,inviter);
        X1(:,k) = tx1; Y1(:,k)=ty1;
        [tx2,ty2] = min_sing_vec(A2-lambda(k)*B2-mu(k)*C2-eta(k)*D2,inviter);
        X2(:,k) = tx2; Y2(:,k)=ty2;
        [tx3,ty3] = min_sing_vec(A3-lambda(k)*B3-mu(k)*C3-eta(k)*D3,inviter);
        X3(:,k) = tx3; Y3(:,k)=ty3;
    end   
end
