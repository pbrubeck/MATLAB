function r = mrdivide(p,q)
%   Divide a polygon by a scalar.

%   Copyright 2003 by Toby Driscoll.
%   $Id: mrdivide.m 298 2009-09-15 14:36:37Z driscoll $

if ~isa(q,'double') || length(q)>1
  error('Function ''/'' defined only for a scalar double.')
end

r = p;
r.vertex = r.vertex/q;
