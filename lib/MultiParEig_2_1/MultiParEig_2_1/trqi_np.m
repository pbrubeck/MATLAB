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
%   - tol : tolerance for the residuals (default is 1e-10)
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

% Last revision: 8.9.2015

if nargin<3 || isempty(maxiter), maxiter = 5; end
if nargin<4 || isempty(tol), tol = 1e-10; end

iter = 0;

k = length(A)-1; % number of parameters
M = zeros(k);

for i = 1:k
    x{i} = x0{i}/norm(x0{i});
    res(i) = 1;
end

for i = 1:k
    for j = 1:k+1
        Ax{i}(:,j) = A{i,j}*x{i};
    end
end

while (iter<= maxiter) && (max(abs(res))>tol)
    
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
        
end    

if max(abs(res))<tol, flag=1; else flag=0; end
    
