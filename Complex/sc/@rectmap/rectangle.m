function zr = rectangle(M)
%RECTANGLE Return the corners of the rectangle in the fundamental domain.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rectangle.m 298 2009-09-15 14:36:37Z driscoll $

zr = prevertex(M);
zr = zr(corners(M));
