function M = mrdivide(M,a)
%   Divide the image of an SC map by a constant.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mrdivide.m 298 2009-09-15 14:36:37Z driscoll $

if isa(a,'double')
  M = M * (1/a);
else
  error('Cannot divide by an SC map.')
end

