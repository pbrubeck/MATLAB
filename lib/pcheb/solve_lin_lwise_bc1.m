function [dof,x,U,xx,C] = solve_lin_lwise_bc1(N,L,rho,xnh,s,f,c)
%SOLVE_LIN_LWISE_BC1  Solve Poisson eqn on L-shaped domain with Dirichlet
%   boundary conditions.  The domain comprises the square [-L,L]^2 minus
%   its lower-right quadrant:
%
%            Uxx + Uyy = f(x,y)  within  [-L,L]^2 \ { [0,L] x [-L,0] },
%                    U = 0         on    all boundaries.
%
%   S invokes symmetry to halve the number of degrees of freedom:
%          s = 0  -  do not assume symmetry
%          s > 0  -  assume positive symmetry across bisector y = -x
%          s < 0  -  assume      antisymmetry across bisector y = -x.
%
%   If s>0 or s<0, F must be symmetric or antisymmetric respectively.
%
%   RHO is a (possibly empty) vector in the range (0,1) specifying 
%       intermediate element boundaries.
%
%   See also SOLVE_LIN_LWISE_BC2.

if isempty(f),  f = @fzero;  end

s = sign(s);

xel = symxel(L,rho);

[x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

[img,~,bis,sup] = img_sq(n,1);  sub = img(sup);

u = zeros(size(xv));  % column vector of length n^2

% construct Poisson operator and forcing function
A = kron(D2,I) + kron(I,D2);   b = f(xv,yv)+c; 

% mortar BCs - enforce continuity of first derivative
A(m,:) = 0;   A = A + M;   b(m) = 0;

% knock out boundaries and lower-right quadrant
dom = (abs(xv)<L & abs(yv)<L) & (xv < 0 | yv > 0);

% reduce DOF by enforcing function symmetry across bisector
if s,
    dom = dom & (sup | ((s>0) & bis));

    A(:,sup) = A(:,sup) + s*A(:,sub);
end

dof = length(find(dom & ~m));

u(dom) = A(dom,dom)\b(dom);

% below-bisector values
if s,  u(sub) = s*u(sup);  end

% convert to matrix (zeros in lower-right quadrant)
U = reshape(u,n,n);

end

