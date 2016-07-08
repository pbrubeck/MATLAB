function r = minus(p,q)
%   Translate a polygon, or subtract the vertices of two polygons.

%   Copyright 1999-2003 by Toby Driscoll.
%   $Id: minus.m 298 2009-09-15 14:36:37Z driscoll $

r = plus(p,-q);
