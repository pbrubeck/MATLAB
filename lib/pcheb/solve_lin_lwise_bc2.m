function [x,U,xx,C] = solve_lin_lwise_bc2(N,L,rho,xnh,sflag)
%SOLVE_LIN_LWISE_BC2  Solve Laplace eqn on L-shaped domain with mixed
%   (Dirichlet and inhomogeneous Neumman) boundary conditions. The domain
%   comprises the square [-L,L]^2 minus its lower-right quadrant:
%
%         Uxx + Uyy = 0   within   [-L,L]^2 \ { [0,L] x [-L,0] },
%                 U = 0     at     y=0 (x>0)  and  x=0  (y<0) ['ELBOW'],
%                Un = 0     at     x=L (y>0)  and  y=-L (x=0),
%                Un = 1     at       x = -L   and   y = +L,
%
%   where Un denotes the wall-normal derivative.
%
%   SFLAG invokes symmetry to halve the number of degrees of freedom. 
%
%   RHO is a (possibly empty) vector in the range (0,1) specifying 
%       intermediate element boundaries.
%
%   This test problem is drawn from Pathria & Karniadakis (1995):
%
%        "Spectral Element Methods for Elliptic Problems in Nonsmooth
%         Domains",  J. Comp. Physics 112, pp. 83-95 (1995).
%
%   TECHNICAL NOTE:
%
%      To apply Neumman conditions at the elbow, it would be necessary to
%      set DFLAG in PCHEB:
%
%        - With SFLAG set, it would suffice to set DFLAG=+1 to ensure that
%          right-derivatives are evaluated at y=0.
%
%        - Without SFLAG set, the y=0 and x=0 elbow conditions would need
%          to be computed using DFLAG=+1 and DFLAG=-1 respectively.
%
%   See also SOLVE_LIN_LWISE_BC1.

xel = symxel(L,rho);

[x,xx,C,D1,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,xv,yv] = mort_sq(x,mx,Mx);

[img,~,bis,sup] = img_sq(n,1);  sub = img(sup);

Dx = kron(D1,I);  Dy = kron(I,D1);

zero = zeros(size(xv));

u = zero;  b = zero;  % column vectors of length n^2

% construct Poisson operator and forcing function
A = kron(D2,I) + kron(I,D2);

% mortar BCs - enforce continuity of first derivative
A(m,:) = 0;   A = A + M;   b(m) = 0;

% BCs at x = -L [left] and x = +L [right]
%   If SFLAG is set, the left BC is superfluous but harmless.
i = find(abs(xv)==L);  A(i,:) = Dx(i,:);  b(xv==-L) = -1;  b(xv==L) = 0;

% BCs at y = -L [bottom] and y = +L [top]
%   If SFLAG is set, the bottom BC is superfluous but harmless.
i = find(abs(yv)==L);  A(i,:) = Dy(i,:);  b(yv==-L) =  0;  b(yv==L) = 1;

% kill two birds with one stone:
%   - knock out lower-right quadrant
%   - inner-wall BCs - prescribe U=0 at y==0 (x >= 0) and x==0 (y <= 0)
dom = (xv < 0 | yv > 0);

if sflag,    
    dom = dom & (sup | bis);

    A(:,sup) = A(:,sup) + A(:,sub);
end

u(dom) = A(dom,dom)\b(dom);

% below-bisector values
if sflag,  u(sub) = u(sup);  end

% convert to matrix (zeros in lower-right quadrant)
U = reshape(u,n,n);

end

