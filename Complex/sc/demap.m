function wp = demap(zp,w,beta,z,c,qdat)
%DEMAP  Schwarz-Christoffel exterior map.
%   DEMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%   Christoffel exterior map at the points in vector ZP. The arguments
%   W, BETA, Z, C, and QDAT are as in DEPARAM.  DEMAP returns a vector
%   the same size as ZP.
%
%   DEMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%   answer accurate to within TOL.
%
%   See also DEPARAM, DEPLOT, DEINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: demap.m 298 2009-09-15 14:36:37Z driscoll $

if isempty(zp)
  wp = [];
  return
end

n = length(w);
beta = beta(:);
z = z(:);

% Quadrature data and error tolerance
if nargin < 6
  tol = 1e-8;
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  tol = qdat;
  qdat = scqdata(beta,max(ceil(-log10(tol)),8));
else
  tol = 10^(-size(qdat,1));
end

shape = size(zp);
zp = zp(:);
zprow = zp.';
p = length(zp);
wp = zeros(p,1);

% For each point in zp, find nearest prevertex.
[dist,sing] = min(abs(zprow(ones(n,1),:) - z(:,ones(1,p))));
sing = sing(:);				% indices of prevertices

% Screen out images of prevertices
vertex = (dist(:) < tol);
wp(vertex) = w(sing(vertex));
zero = abs(zp) < tol;
wp(zero) = Inf;
vertex = vertex | zero;

% zs = the starting singularities
zs = z(sing);

wp(~vertex) = w(sing(~vertex));

% Must be careful about the singularity at the origin, since the
% quadrature routine doesn't pay attention to the right endpoint.

abszp = abs(zp); 			% dist to sing at 0

% Integrate for the rest, using vectorization with masking.
unf = ~vertex;        % unfinished cases
zold = zs;            % integration starting points
znew = zold; dist = ones(p,1);     % initial conditions
while any(unf)
  % How far can the integration go?
  dist(unf) = min(1,2*abszp(unf)./abs(zp(unf)-zold(unf)));
  % New integration endpoints
  znew(unf) = zold(unf) + dist(unf).*(zp(unf)-zold(unf));
  wp(unf) = wp(unf) + c*dequad(zold(unf),znew(unf),sing(unf),z,beta,qdat);

  unf = (dist<1);
  zold(unf) = znew(unf);
  sing(:) = 0;        % only first step can have endpoint singularity
end

wp = reshape(wp,shape);
