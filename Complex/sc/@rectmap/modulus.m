function mu = modulus(M)
%MODULUS Conformal modulus of the generalized quadrilateral.
%   Returns the conformal modulus of the polygon in the rectmap (the
%   aspect ratio of the source rectangle).

%   Copyright 1998 by Toby Driscoll.
%   $Id: modulus.m 298 2009-09-15 14:36:37Z driscoll $

z = M.prevertex;
mu = max(imag(z)) / (2*max(real(z)));
