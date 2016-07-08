function q = uminus(p)
%   Negate the vertices of a polygon.
%   This may have surprising consequences if p is unbounded.

%   Copyright 2003 by Toby Driscoll (driscoll@math.udel.edu).
%   $Id: uminus.m 298 2009-09-15 14:36:37Z driscoll $

q = polygon( -vertex(p), angle(p) );

