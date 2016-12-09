function [lambda,x,Modes,xx,C] = solve_eig_sq(N,xel,xnh,Neig)
%SOLVE_EIG_SQ  Laplace/Dirichlet eigenmodes on a square domain:
%
%     Uxx + Uyy + lambda*U = 0   within  Dom [a0,a1]^2,
%                        U = 0     on    boundary[Dom],
%
%     where lambda is an unknown constant (i.e. eigenvalue).  The Neig
%     leading modes are returned.
%
% [lambda,x,Modes,xx,C] = solve_eig_sq(N,xel,xnh,Neig)
%
% See also SOLVE_EIG_DUCT.

a0 = xel(1);  a1 = xel(end);

[x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
i    = ~wall;
Eqns = spdiag(~m & ~wall); 
% This would also work in the case of Dirichlet BCs:
%    Eqns = spdiag(~m);

A = M + Eqns*(kron(D2,I) + kron(I,D2));
B = double(Eqns);

[modes,Eig] = eigs(A(i,i), B(i,i), Neig,'sm'); 

% eigenvalues are all real, but numerical/roundoff error may generate
%   spurious imaginary part ...
lambda = -real(diag(Eig));

Modes = zeros(n,n,Neig);

for j = 1:Neig,
    Modes(2:n-1,2:n-1,j) = reshape(modes(:,j), n-2,n-2);
end
