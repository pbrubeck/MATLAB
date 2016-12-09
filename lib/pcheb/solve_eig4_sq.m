function [lambda,x,Modes,xx,C] = solve_eig4_sq(N,xel,xnh,Neig)
%SOLVE_EIG_SQ  Biharmonic eigenmodes on [-1,1]^2 with Dirichlet and Neumann
%              boundary conditions.  The eigenvalues
%
%  [lambda,x,Modes,xx,C] = solve_eig4_sq(N,xel,xnh,Neig)
%
%      LAMBDA is a vector of length Neig denoting the (square roots of) the
%      Neig leading eigenvalues, rescaled so that LAMBDA(1) = 1.

% See also SOLVE_LIN4_1D.

a0 = xel(1);  a1 = xel(end);

% choose between Chebyshev (0) and Legendre (1) points
xflag = 0;

% choose symmetry of D-matrices - any one of {0,-1,+1}
dflag = 0;

% must use 4th-order mortar BCs
mflag = 1;

[x,xx,C,~,D2x,mx,Mx,D3x,D4x] = pcheb(N,xel,xnh,xflag,dflag,mflag);

n = length(x);  i = 2:n-1;  I = speye(n);

% Implement Neumann BCs (see SOLVE_LIN4_1D for detailed explanation).
S  = diag([ 0; 1./(1 - x(i).^2); 0; ]); 
D4 = (diag(1-x.^2)*D4x - 8*diag(x)*D3x - 12*D2x)*S;

[m,M,xv,yv] = mort_sq(x,mx,Mx);

wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
i    = ~wall;
Eqns = spdiag(~m & ~wall); 
% This would also work in the case of Dirichlet BCs:
%    Eqns = spdiag(~m);

A = M - Eqns*(kron(D4,I) + kron(I,D4) + 2*kron(D2x,D2x));
B = double(Eqns);

[modes,Eig] = eigs(A(i,i), B(i,i), Neig,'sm'); 

kappa  = -diag(Eig);
lambda = sqrt(kappa/kappa(1));
Modes  = zeros(n,n,Neig);

for j = 1:Neig,
    Modes(2:n-1,2:n-1,j) = reshape(modes(:,j), n-2,n-2);
end
