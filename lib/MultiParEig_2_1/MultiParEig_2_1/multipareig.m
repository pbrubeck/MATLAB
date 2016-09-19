function [lambda,X,Y] = multipareig(A,opts)

%MULTIPAREIG   Solve a multiparameter eigenvalue problem
%
% [lambda,X,Y] = MULTIPAREIG(A,opts) returns eigenvalues and eigenvectors 
% of the multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% Input:
%   - A : cell array of size k x (k+1) of matrices A{i,j}, all matrices in
%         the same row have to be square matrices of the same size
%   - opts : options
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - inviter (1) : use inverse iteration for eigenvectors or slow svd (0)
%   - all options of auxiliary functions
%
% Output:
%   - lambda : matrix of size m x k, each row is an eigenvalue
%   - X : cell array of size m x k with right eigenvectors
%   - Y : cell array of size m x k with left eigenvectors
% 
% Method can fail in case of multiple eigenvalue lambda(1). In such case, 
% for k=2 and k=3 use twopareig and threepareig with option fast = 0
%
% See also: TWOPAREIG, THREEPAREIG.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 05.09.2015 : support for singular MEP
% Last revision: 08.09.2015

if nargin<2, opts=[]; end
if isfield(opts,'inviter'),   inviter = opts.inviter;    else inviter = 1;   end
if isfield(opts,'singular'),  singular = opts.singular;  else singular = 0;  end

k = length(A); % number of parameters

% computation of Delta matrices
Delta = multipar_delta(A);

if singular
    Delta = extract_regular_part_np(Delta, opts);
end

if (size(Delta{1},1) == 0) || (size(Delta{1},1)~=size(Delta{1},2))   % no regular part was found
     lambda = []; X = []; Y = []; % default output
     return
end

% computation of eigenvalues, we assume that the unitary transformation
% that puts inv(Delta0)*Delta1 into triangular form, does this to other
% matrices inv(Delta0)*Deltaj as well
Gamma1 = - Delta{1}\Delta{2};
[Q,L] = schur(Gamma1,'complex');
lambda(:,1) = diag(L);
for r = 2:k-1
    Gammar = - Delta{1}\Delta{r+1};
    lambda(:,r) = diag(Q'*Gammar*Q);
end

% extraction of eigenvectors (individually using inverse iteration)
if nargout > 1
    m = size(lambda,1);
    X = cell(m,k-1);
    Y = cell(m,k-1);
    for j = 1:m
        for r = 1:k-1
            M = A{r,1};
            for p = 1:k-1
                M = M - lambda(j,p)*A{r,p+1};
            end
           [zr,zl] = min_sing_vec(M,inviter);
           X{j,r} = zr;
           Y{j,r} = zl;
        end
    end   
end



