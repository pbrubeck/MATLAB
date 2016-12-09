function [dof,x,U,xx,C] = ...
    solve_lin_sq_gen(N,xel,xnh,fbv,f,c,f0,c0,fx,cx,fy,cy,fxy,cxy)
%SOLVE_LIN_SQ_GEN  Solve a general 2nd-order linear PDE on a square domain
%                  with prescribed boundary values.
%   Domain:
%      Dom = [a0,a1]^2, where [a0,a1] are the end-points of XEL,
%   PDE:
%       Uxx + Uyy + (fxy+cxy)*Uxy 
%       + (fx+cx)*Ux + (fy+cy)*Uy + (f0+c0)*U = f+c   [within Dom],
%   Boundary:
%              U = fbv    on    x=x0/x=x1/y=x0/y=x1,
%
%   where [fbv,f,f0,fx,fy,fxy] are functions of (x,y), but not of U.
%
%  SIGNATURE:
%
%  [dof,x,U,xx,C] = ...
%          solve_lin_sq_gen(N,xel,xnh,fbv,f,c,f0,c0,fx,cx,fy,cy,fxy,cxy)
%
%  Null functions can be passed as "[]". Excess terms can be omitted.
%  For example, the following code solves a Helmholtz/Dirichlet eqn
%    with wavenumber K:
%
%     [...] = solve_lin_sq_gen(N,xel,xnh,[],f,c,[],k^2);
%
%  See also SOLVE_LIN_SQ_SYM and SOLVE_LIN_RECT.

i=2;  % number of mandatory arguments

i=i+1;  if nargin < i,                  xnh = [];      end
i=i+1;  if nargin < i || isempty(fbv),  fbv = @fzero;  end
i=i+1;  if nargin < i || isempty(f  ),  f   = @fzero;  end
i=i+1;  if nargin < i || isempty(c  ),  c   = 0;       end
i=i+1;  if nargin < i || isempty(f0 ),  f0  = @fzero;  end
i=i+1;  if nargin < i || isempty(c0 ),  c0  = 0;       end
i=i+1;  if nargin < i || isempty(fx ),  fx  = @fzero;  end
i=i+1;  if nargin < i || isempty(cx ),  cx  = 0;       end
i=i+1;  if nargin < i || isempty(fy ),  fy  = @fzero;  end
i=i+1;  if nargin < i || isempty(cy ),  cy  = 0;       end
i=i+1;  if nargin < i || isempty(fxy),  fxy = @fzero;  end
i=i+1;  if nargin < i || isempty(cxy),  cxy = 0;       end

a0 = xel(1);  a1 = xel(end);

[x,xx,C,D,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

[m,M,p,q] = mort_sq(x,mx,Mx);

w   = (p==a0 | p==a1 | q==a0 | q==a1);  % boundary ('wall') points
pde = (~w & ~m);
dof = length(find(pde));

PDE = kron(D2,I) + kron(I,D2)    ...
    + spdiag(fxy(p,q)+cxy)*kron(D,D) ...
    + spdiag( fx(p,q)+cx )*kron(D,I) ...
    + spdiag( fy(p,q)+cy )*kron(I,D) ...
    + spdiag( f0(p,q)+c0 )*speye(n^2); % i.e. spdiag(...)*kron(I,I)

A = M + spdiag(pde)*PDE + spdiag(w);
b =     spdiag(pde)*(f(p,q)+c) + w.*fbv(p,q);

U = reshape(A\b,n,n);
