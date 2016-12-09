function [dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f,c)
%SOLVE_LIN_SQ_SYM  Solve Poisson/Dirichlet equation on square domain:
%
%       Uxx + Uyy = f(x,y) + c   within  Dom = [a0,a1]^2,
%               U = 0              on    boundary[Dom]
%
%   where [a0,a1] are the end-points of XEL.
%
%   S invokes symmetry to halve the number of degrees of freedom:
%          s = 0  -  do not assume symmetry
%          s > 0  -  assume     symmetry across bisector y = x
%          s < 0  -  assume antisymmetry across bisector y = x.
%
%   If s>0 or s<0, F must be symmetric or antisymmetric respectively.
%
%  SIGNATURE:
%
%  [dof,x,U,xx,C] = solve_lin_sq_sym(N,xel,xnh,s,f)
%
%  See also SOLVE_LIN_SQ_GEN and SOLVE_LIN_RECT.

if isempty(f),  f = @fzero;  end

s = sign(s);  a0 = xel(1);  a1 = xel(end);

[x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

[img,sub,bis] = img_sq(n,0);  sup = img(sub);

wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
dom  = ~wall;

u = zeros(size(xv));  % column vector of length n^2

% construct Poisson operator and forcing function
A = kron(D2,I) ...  %   Dxx
  + kron(I,D2);     % + Dyy

b = f(xv,yv)+c; 

% mortar BCs - enforce continuity of first derivative
A(m,:) = 0;   A = A + M;   b(m) = 0;

% reduce DOF by enforcing function symmetry across bisector
if s,
    dom = dom & (sub | ((s>0) & bis));

    A(:,sub) = A(:,sub) + s*A(:,sup);
end

% solve
dof = length(find(dom & ~m));

tic
u(dom) = A(dom,dom)\b(dom);
toc

% below-bisector values
if s,  u(sup) = s*u(sub);  end

% convert to matrix
U = reshape(u,n,n);
