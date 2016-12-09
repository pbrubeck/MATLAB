function [x,U,xx,C] = solve_lin_sq_lrefine(N,xel,rho,xnh,f,c)
%SOLVE_LIN_SQ_LREFINE  Solve Poisson/Dirichlet eqn on square domain with 
%   local refinement of the mesh as follows:
%
%       corner domain (inner region):  { x0 <= x & y <= x1 }
%    principal domain (outer region):  { x1 <= x | y <= xM }
%
%   where x0 = xel(1),
%         x1 = xel(2),
%         xM = xel(end).
%
%   This scheme is classified as non-conformal, meaning that the
%   collocation points at the inner/outer inteface do not line up:
%
%     - The outer domain is discretised using PCHEB(N,xel),
%         i.e. using (M^2 - 1) elements, M = length(xel) - 1.
%       Since the lower-left element (corner region) is to be excised, it
%         is necessary to set DFLAG in PCHEB.
%
%     - The inner (corner) domain is discretised using PCHEB(N,xelc):
%            xelc = x0 + (x1-x0)*[0,rho,1],  0 < rho < 1,
%        i.e. using Mc^2 elements, Mc = length(rho) + 1.

if isempty(f),  f = @fzero;  end

% choose between Chebyshev (0) and Legendre (1) points
xflag = 0;

% matrix D1d (outer domain) requires right-derivatives at interface x1
dflag = +1;

rho = rho(rho > 0 & rho < 1);  rho = [ 0; rho(:); 1 ];

xeld = xel(:);
x0   = xel(1);
x1   = xel(2);
xend = xel(end);
xelc = x0 + (x1 - x0)*rho;   xelc(end) = x1;

y1 = x1;

cxnh = pcheb(N,[x0,x1]);  % the N+1 join-points in the D-domain

% Corner operators (inner region, refined grid)
[xc,~,Cc,D1c,D2c,mxc,Mxc] = pcheb(N,xelc,cxnh,xflag);   dxnh = xc;

% Domain operators (outer region, principal grid)
[xd,~,Cd,D1d,D2d,mxd,Mxd] = pcheb(N,xeld,dxnh,xflag,dflag);

[mc,Mcc,xvc,yvc] = mort_sq(xc,mxc,Mxc);
[md,Mdd,xvd,yvd] = mort_sq(xd,mxd,Mxd);

nc = length(xc);  Ic = speye(nc);  ncsq = nc^2;
nd = length(xd);  Id = speye(nd);  ndsq = nd^2;

xv = [ xvc; xvd ];
yv = [ yvc; yvd ];

wall = (xv==x0 | xv==xend | yv==x0 | yv==xend);

% Will eventually discard corner points in D-domain
ckeep = true(size(xvc));
dkeep = (xvd >= x1 | yvd >= y1);
keep  = [ ckeep; dkeep ];

i1 =  1:ncsq;
i2 = (1:ndsq) + ncsq;

% Cc [size (N+1)-by-nc] - interpolate xc to xd for 0 <= xd <= x1
% Cd [size    nc-by-nd] - interpolate xd to xc

Idj =  Id(xd==x1,:);
Dcj = D1c(xc==x1,:);

ixc1 = find(xvc==x1);  ixd1 = find(xvd==x1 & yvd<=y1);
iyc1 = find(yvc==y1);  iyd1 = find(yvd==y1 & xvd<=x1);

Ddx = kron(D1d, Id);
Ddy = kron( Id,D1d);

Ac = kron(D2c,Ic) + kron(Ic,D2c);
Ad = kron(D2d,Id) + kron(Id,D2d);

Mcd = sparse(ncsq,ndsq);
Mdc = sparse(ndsq,ncsq);

% enforce continuity at x==x1, y==yc  (interpolated to y==yd in D-domain)
i = ixc1;  Mcc(i,i) = -speye(nc);  Mcd(i,:) = kron(Idj,Cd);  mc(i) = true;

% enforce continuity at y==y1, x==xc  (interpolated to x==xd in D-domain)
i = iyc1;  Mcc(i,i) = -speye(nc);  Mcd(i,:) = kron(Cd,Idj);  mc(i) = true;

%    d/dx continuity at x==x1, y==yd  (interpolated to y==yc in C-domain)
i = ixd1;  Mdd(i,:) = -Ddx(i,:);   Mdc(i,:) = kron(Dcj,Cc);

%    d/dy continuity at y==y1, x==xd  (interpolated to x==xc in C-domain)
i = iyd1;  Mdd(i,:) = -Ddy(i,:);   Mdc(i,:) = kron(Cc,Dcj);

m   = [ mc; md ];
int = (~wall & ~m);
b   = int.*(f(xv,yv)+c); 
u   = zeros(size(b));
A   = ...
      spdiag( int ) *   blkdiag(Ac,Ad)     ...
    + spdiag(~wall) * [ Mcc Mcd; Mdc Mdd ] ...
    + spdiag( wall);

u(keep) = A(keep,keep)\b(keep);


Uc = reshape(u(i1), nc,nc);
Ud = reshape(u(i2), nd,nd);

% interpolate Ud to C-domain grid (xc for x<x1, yc for y<y1)
xel = [ xelc; xeld(xeld>x1) ];   

i = find(xd > x1);   

x = [ xc; xd(i) ];

C = [ Cd; Id(i,:) ];   U = expand2d(Ud,C);   U(1:nc,1:nc) = Uc;

% interpolation to fine grid
[C,xx] = pbary(N,xel,xnh);
