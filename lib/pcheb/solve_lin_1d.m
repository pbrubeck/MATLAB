function [x,u,xx,uu] = solve_lin_1d(N,xel,xnh,func)
%SOLVE_LIN_1D  Solve 1D Poisson/Dirichlet equation:
%                Uxx = 0   within  [a0,a1],
%                  U = 0   at x=a0, x=a1
%
%   where [a0,a1] are the end-points of XEL.

[x,xx,C,~,D2,m,M] = pcheb(N,xel,xnh); 

% 1D Laplace equation
A = D2;  b = func(x);  

% mortar BCs - continuity of first derivative between elements
A(m,:) = 0;  A = A + M;  b(m) = 0;

% Dirichlet conditions - implement by deleting boundary points
A = A(2:end-1,2:end-1);  b = b(2:end-1);

% solve and restore boundary points
u = [ 0; A\b; 0 ];   

% interpolate to fine grid
uu = C*u;

end
