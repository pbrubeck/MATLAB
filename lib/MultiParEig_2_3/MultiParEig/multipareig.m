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
%   - rrqr (0): for singular problems only, set to 1 to use rank revealing qr
%   - inviter (1) : use inverse iteration for eigenvectors or slow svd (0)
%   - all options of auxiliary functions
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
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
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 05.09.2015 : support for singular MEP
% PH 01.11.2016 : e-values computation speed-up
% PH 22.11.2016 : precision-independent version 
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications and clean-ups.

% Last revision: 26.11.2016

% Validate number of input parameters
narginchk(1, 2);

% Analyse user supplied options, if any.
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A{:});
end

if isfield(opts,'inviter'),   inviter = opts.inviter;    else,  inviter = 1;  end
if isfield(opts,'singular'),  singular = opts.singular;  else, singular = 0;  end

% Make sure all inputs are of the same numeric type.
for k=1:numel(A)
    if ~isa(A{k}, class_t)
         A{k} = numeric_t(A{k},class_t);
    end
end

k = length(A); % number of parameters + 1

% computation of Delta matrices
Delta = multipar_delta(A);

if singular
    Delta = extract_regular_part_np(Delta, opts);
end

if (size(Delta{1},1) == 0) || (size(Delta{1},1)~=size(Delta{1},2))   % no regular part was found
     lambda = numeric_t([],class_t); X = numeric_t([],class_t); Y = numeric_t([],class_t);   % default output
     return
end

% computation of eigenvalues, we assume that the unitary transformation
% that puts inv(Delta0)*Delta1 into triangular form, does this to other
% matrices inv(Delta0)*Deltaj as well

% Factorize Delta{1} only once and use factorization to solve against
% different righ-hand-sides later on. This boosts the speed slightly.
[L,U,p] = lu(Delta{1},'vector');
Gamma1 = -U\(L\Delta{2}(p,:));

[Q,T] = schur(Gamma1,'complex');
lambda(:,1) = diag(T);
for r = 2:k-1
    % Avoid O(n^3) operations in the loop
    Gammar = (-U\(L\Delta{r+1}(p,:)))*Q;
    for n = 1:size(Q,1)
        lambda(n,r) = Q(:,n)'*Gammar(:,n);
    end
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

end % multipareig