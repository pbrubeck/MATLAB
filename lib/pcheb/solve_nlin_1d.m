function [xx,uu,x,u,n] = solve_nlin_1d(N,xel,nh,func,tol)
%SOLVE_NLIN_1D  Solve the nonlinear 1D equation Uxx = func(U).


[x,xx,C,~,D2,m,M] = pcheb(N,xel,nh);  

% 1D Laplace equation
A = D2;

% mortar BCs - continuity of first derivative between elements
A(m,:) = 0;  A = A + M;

% Dirichlet conditions - implement by deleting boundary points
A = A(2:end-1,2:end-1);

du = 1.0;  u0 = zeros(size(x));  n=0;


while du > tol
    n  = n+1;
    b  = func(u0);  b(m) = 0;  b = b(2:end-1);
    u  = [ 0; A\b; 0 ];
    du = max(abs(u - u0));  u0 = u;
end

% interpolate to fine grid
uu = C*u;


end
