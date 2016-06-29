function corner = corners(M)
%CORNERS Indices of rectangle/generalized quadrilateral corners.

%   Copyright 1998 by Toby Driscoll.
%   $Id: corners.m 298 2009-09-15 14:36:37Z driscoll $

z = M.prevertex;
tol = 4*eps;

% Find extent of rectangle
K = max(real(z));
Kp = max(imag(z));

% First corner is K + 0i
dif = repmat(z,1,4) - repmat([K K+i*Kp -K+i*Kp -K],length(z),1);
[tmp,corner] = min(abs(dif));

corner = corner(:);
