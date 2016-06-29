function r = mtimes(p,q)
%   Multiplication of a polygon by a scalar.

%   Copyright 1998 by Toby Driscoll.
%   $Id: mtimes.m 298 2009-09-15 14:36:37Z driscoll $

if isa(q,'polygon')
  if isa(p,'polygon')
    error('Function ''*'' not defined for two polygon objects.')
  end
  tmp = p;
  p = q;
  q = tmp;
end

r = p;
r.vertex = r.vertex*q;
