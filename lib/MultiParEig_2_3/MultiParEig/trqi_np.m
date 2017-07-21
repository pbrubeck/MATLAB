function [lambda,x,res,flag,iter] = trqi_np(A,x0,maxiter,tol)

%TRQI_NP  Tensor Rayleigh quotient iteration for a multi-parameter eigenvalue problem 
%
% [lambda,x,res,flag,iter] = TRQI_NP(A,x0,maxiter,tol) applies Tensor 
% Rayleigh quotient iteration to a multiparameter eigenvalue problem

% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% Input:
%   - A : cell array of size k x (k+1) of matrices A{i,j}, all matrices in
%         the same row have to be square matrices of the same size
%   - x0 : cell of size k of vectors x0{i}, their tensor product is an
%         approximation for the eigenvector
%   - maxiter : maximum number of iterations (default is 5)
%   - tol : tolerance for the residuals (default is 1e6*eps)
%
% Output:
%   - lambda : eigenvalue (vector of size k)
%   - x : eigenvector (cell of size k of vectors x{k})
%   - res : norms of the residual
%   - flag : convergence flag (success 1 or fail 0)
%   - iter : number of iterations
%
% See also: TRQI, TRQI_3P

% Reference: Algorithm 2 (TRQI) from
% B. Plestenjak, A continuation method for a right definite two-parameter
% eigenvalue problem, SIAM J. Matrix Anal. Appl. 21 (2000), 1163-1184.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 04.12.2016: new stopping criteria when residual stops decreasing
% BP 01.12.2016: modified to be precision-independent                 
% BP 07.11.2016: We consider Frobenius norm of the residuals 
% Last revision: 04.12.2016

% Validate number of input parameters
narginchk(2,4);

class_t = superiorfloat(A{:},x0{:});

% Make sure all inputs are of the same numeric type.
for k=1:numel(A)
    if ~isa(A{k},class_t)
         A{k} = numeric_t(A{k},class_t);
    end
end
for k=1:numel(x0)
    if ~isa(x0{k}, class_t)
         x0{k} = numeric_t(x0{k},class_t);
    end
end

if nargin<3 || isempty(maxiter), maxiter = 5; end
if nargin<4 || isempty(tol), tol = numeric_t('1e6*eps',class_t); end

iter = 0;

k = length(A)-1; % number of parameters
M = zeros(k,class_t);

for i = 1:k
    x{i} = x0{i}/norm(x0{i});
    res(i) = numeric_t('Inf',class_t);
end

for i = 1:k
    for j = 1:k+1
        Ax{i}(:,j) = A{i,j}*x{i};
    end
end

normres = 1e20;
oldnormres = Inf;

while (iter<= maxiter) && (norm(res)>tol) && (normres<0.95*oldnormres)

    for i = 1:k
        M(i,:) = x{i}'*Ax{i}(:,2:k+1);
        b(i,1) = x{i}'*Ax{i}(:,1);
    end
    lambda = M\b;
    iter = iter + 1;
    warning off
    for i = 1:k
        W{i} = A{i,1};
        for j = 1:k
            W{i} = W{i} - lambda(j)*A{i,j+1};
        end
        V{i} = W{i}\Ax{i}(:,2:k+1);
    end
    warning on
    
    usesvd = 0;
    for i = 1:k
        if any(any(isnan(V{i}))) || any(any(isinf(V{i})))
            usesvd = 1;
        end
    end
    
    if usesvd
        for i = 1:k
            x{i} = min_sing_vec(W{i});
        end
    else
        for i = 1:k
            M(i,:) = x{i}'*V{i};
        end
        warning off
        delta = M\ones(k,1);
        warning on
        for i = 1:k
            xnew{i} = V{i}*delta;
            x{i} = xnew{i}/norm(xnew{i});
        end
    end
    
    for i = 1:k
        for j = 1:k+1
            Ax{i}(:,j) = A{i,j}*x{i};
        end
    end
    
    for i = 1:k
        res(i) = norm(W{i}*x{i});
    end
    oldnormres = normres;
    normres = norm(res);
end    

if norm(res)<tol, flag=1; else flag=0; end
    
