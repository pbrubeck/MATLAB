function [dof,x,U,xx,C] = solve_nlin_sq(N,xel,xnh,func,tol)
%SOLVE_NLIN_SQ  Solve a nonlinear PDE on a square domain:
%                    Uxx + Uyy = func(U), Dirichlet BCs.

if nargin < 5,  tol = 1e-6;  end

a0 = xel(1);  a1 = xel(end);

[x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
dom  = ~wall;

zero = zeros(size(xv));

u = zero;  u0 = zero; % column vector of length n^2

% construct Poisson operator and forcing function
A = kron(D2,I) + kron(I,D2); 

% mortar BCs - enforce continuity of first derivative
A(m,:) = 0;   A = A + M;

% solve
dof = length(find(dom & ~m));

du = 1.0;  k=0;


while du > tol
    k  = k+1;
    b  = func(u0);  b(m) = 0;
    u(dom) = A(dom,dom)\b(dom);
    du = max(abs(u - u0));  u0 = u;
end

% convert to matrix
U = reshape(u,n,n);
  
disp(['Converged after ',int2str(k),' steps.'])
