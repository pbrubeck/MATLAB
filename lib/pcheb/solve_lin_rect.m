function [dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(Nx,Ny,xel,yel,xnh,ynh,f,c,fbv)
%SOLVE_LIN_RECT  Solve Poisson equation on rectangular domain with
%                prescribed boundary values.
%
%       Uxx + Uyy = f(x,y)   within  Dom =  { [a0,a1] x [b0,b1] },
%               U = fbv(x,y)   on    boundary[Dom]
%
%       where [a0,a1] are the end-points of XEL,
%         and [b0,b1] are the end-points of YEL.
%
%  SIGNATURE:
%
%  [dof,x,y,U,xx,yy,Cx,Cy] = solve_lin_rect(Nx,Ny,xel,yel,xnh,ynh,f,fbv)
%
%  See also SOLVE_LIN_SQ_SYM.

i=4;  % number of mandatory arguments

i=i+1;  if nargin < i,  xnh = [];  end
i=i+1;  if nargin < i,  ynh = [];  end

i=i+1;  if nargin < i || isempty(f),  f = @fzero;  end

i=i+1;  if nargin < i,  c = 0;  end

i=i+1;  if nargin < i,  fbv = @fzero;  end

a0 = xel(1);  a1 = xel(end);
b0 = yel(1);  b1 = yel(end);

[x,xx,Cx,~,D2x,mx,Mx] = pcheb(Nx,xel,xnh);  nx = length(x);
[y,yy,Cy,~,D2y,my,My] = pcheb(Ny,yel,ynh);  ny = length(y);

[m,M,xv,yv] = mort_rect(x,y,mx,my,Mx,My);

w   = (xv==a0 | xv==a1 | yv==b0 | yv==b1);  % boundary ('wall') points
dof = length(find(~w & ~m));
Ix  = speye(nx);
Iy  = speye(ny);

% Step 1:  Construct Poisson operator and forcing function
A = kron(D2x,Iy) + kron(Ix,D2y);   

b = f(xv,yv)+c; 

% Step 2a:  Impose mortar BCs
A(m,:) = 0;   A = A + M;   b(m) = 0;

% Step 2b:  Impose boundary values
A(w,:) = 0;   A(w,w) = speye(length(find(w)));  b(w) = fbv(xv(w),yv(w));
%
% NOTES:
%
%  (1) Steps 2a and 2b can be safely interchanged, since MORT_RECT avoids
%      boundary points.
%
%  (2) Steps 1,2a,2b can be combined as follows:
%
%              pde = (~w & ~m);
% 
%              A = M + spdiag(w) ...
%                    + spdiag(pde) * (kron(D2x,Iy) + kron(Ix,D2y));
% 
%              b = spdiag(w)*fbv(xv,yv) + spdiag(pde)*f(xv,yv); 

U = reshape(A\b,ny,nx);


end

