function [xr,xl] = min_sing_vec(A,fast,xr0,xl0)

%MIN_SING_VEC   Left and right zero singular vector of a singular matrix A 
%
% [xr,xl] = MIN_SING_VEC(A,fast,xr0,xl0) returns the right and 
% the left singular vector for the smallest singular value of matrix A,
% for which we assume that it is almost singular and do one step of inverse 
% iteration
%
% fast = 0 : use SVD (slow)
% fast = 1 : (default) use one step of inverse iteration (fast)
% xr0, xl0 : initial vectors for inverse iteration (otherwise random)

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 03.11.2016: speedup: when xl0 is not given, we use xr0=xl0 (if xl0 is random, we use the same random vector) 
% PH 22.11.2016: fixed bug when initial vectors supplied, code clean-up.

% Last revision: 22.11.2016

if nargin<2 || isempty(fast)
    fast = 1;
end

if ~fast % slow computation via SVD
    [U,S,V] = svd(A); %#ok<*ASGLU>
    xr = V(:,size(A,2));
    xl = U(:,size(A,2));
    return
end

n = size(A,1);

if nargin<3 || isempty(xr0)
    xr0 = randn(n,1,class(A));
end

if nargin<4 || isempty(xl0)
    xl0 = randn(n,1,class(A));
end

steps = 1;

% in inverse iteration we don't want the warnings about rank deficiency or
% singular matrices
warning off 

if ~issparse(A)
    [L,U,p] = lu(A,'vector');
    [tilda,invp] = sort(p);
    opts1.LT = true;
    opts2.UT = true;
    
    xr = xr0;
    for j=1:steps
        s = linsolve(L,xr(p),opts1);
        y0 = linsolve(U,s, opts2);
        xr = y0/norm(y0);
    end
    
    opts1.TRANSA=true;
    opts2.TRANSA=true;
    
    xl = xl0;
    for j=1:steps
        s = linsolve(U,xl,opts2);
        y0 = linsolve(L,s,opts1);
        y0 = y0(invp);
        xl = y0/norm(y0);
    end
else
    y0 = A\xr0;
    xr = y0/norm(y0);
    
    y0 = A'\xl0;
    xl = y0/norm(y0);
end

if ~all(isfinite([xr; xl]))
    % solution for a rare situation where we get NaN or Inf
    [tmpU,tmpS,tmpV] = svd(A);
    xr = tmpV(:,end);
    xl = tmpU(:,end);
end

warning on

   