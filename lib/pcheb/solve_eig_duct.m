function [lambda,x,Modes,xx,C] = solve_eig_duct(N,R1,R2,xnh,Neig)
%SOLVE_EIG_DUCT  Laplace/Dirichlet eigenmodes on a duct, defined as the 
%    square-annular domain [-R2,R2]^2 \ [-R1,R1]^2 for 0 < R1 < R2.
%
% [lambda,x,Modes,xx,C] = solve_eig_duct(N,R1,R2,xnh,Neig)
%
% See also SOLVE_EIG_SQ.

xel = [-R2 -R1 R1 R2];

[x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

xabs = abs(xv);
yabs = abs(yv);

dom = (xabs<R2 & yabs<R2) & (xabs>R1 | yabs>R1);

% wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
% i    = ~wall;
Eqns = spdiag(~m & dom); 
% This would also work in the case of Dirichlet BCs:
%    Eqns = spdiag(~m);

A = M + Eqns*(kron(D2,I) + kron(I,D2));
B = double(Eqns);

[modes,Eig] = eigs(A(dom,dom), B(dom,dom), Neig,'sm'); 

lambda = -diag(Eig);

Mode  = zeros(n,n);
Modes = zeros(n,n,Neig);

for j = 1:Neig,
    
    Mode(dom) = modes(:,j);   Modes(:,:,j) = Mode;
    
end
